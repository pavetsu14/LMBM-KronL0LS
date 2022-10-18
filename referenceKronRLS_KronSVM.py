from rlscore.learner import CGKronRLS
from rlscore.learner.kron_svm import KronSVM
from rlscore.measure import cindex
from rlscore.kernel import GaussianKernel
from rlscore.utilities.sampled_kronecker_products import sampled_vec_trick

import warnings
warnings.filterwarnings("ignore")

import numpy as np
import time
import pandas as pd
import multiprocessing as mp
import itertools as it

import load_split_data
"""
Functions to compute the kernel matrices either explicitly (naive_kronecker_kernel) 
or implicitly (sampled_vec_trick).
"""
def naive_kronecker_kernel(K1, K2, rows1, cols1, rows2, cols2):
    assert len(rows1) == len(cols1)
    assert len(rows2) == len(cols2)
    o = len(rows1)
    p = len(rows2)
    K = np.zeros((o, p))
    for i in range(o):
        for j in range(p):
            k_ij = K1[rows1[i], rows2[j]]
            g_ij = K2[cols1[i], cols2[j]]
            val = k_ij * g_ij
            K[i,j] = val
    return K

def mv_kronecker(K1, K2, rows1, cols1, rows2, cols2):
    def mv(v):
        return sampled_vec_trick(v, K2, K1, cols1, rows1, cols2, rows2)
    return mv

"""
Callback object so that KronRLS and KronSVM can be run with early stopping.
"""
class CallBack(object):

    def __init__(self, K1, K2, Y, row_inds, col_inds, ESlag = 10):
        self.K1 = K1
        self.K2 = K2
        self.Y = Y
        self.row_inds = row_inds
        self.col_inds = col_inds
        self.iter = 1
        self.maxperf = 0
        self.earlyStopLag = ESlag
        self.best_learner = None

    def callback(self, learner):
        P = learner.predict(self.K1, self.K2, self.row_inds, self.col_inds)
        perf = cindex(self.Y, P)
        # Save the model with the highest c-index.
        if perf > self.maxperf:
            self.maxperf = perf
            self.best_iter = self.iter
            self.best_learner = learner.predictor
        self.iter += 1
        # Raise an error if there has not been any improvement in C-index in ESlag iterations.
        if self.iter > self.best_iter + self.earlyStopLag:
            raise ValueError("Early stopping criteria met", self.best_iter, self.maxperf, self.best_learner)
            
    def finished(self, learner):
        pass

"""
The process of obtaining the reference models with the state-of-the art methods. 
Input is a list of 11 elements that are saved in variables with more descriptive name.
"""
def referenceModels(reference_params):
    Y = reference_params[0]
    # It is not explicitly given whether labels are binary or continuous.
    binary = int(len(set(Y)) == 2)
    drug_inds = reference_params[1]
    target_inds = reference_params[2]
    training_inds = reference_params[3]
    test_inds = reference_params[4][0]
    validation_inds = reference_params[4][1]
    lambdas = reference_params[5]
    MI = reference_params[6]
    XD = reference_params[7]
    XT = reference_params[8]
    data_name = reference_params[9]
    random_seed = reference_params[10]
    zero_limit = 1e-4
    
    # Split the indices of the drugs, targets and Y according to the training-test-validation split.
    train_drug_inds = drug_inds[training_inds]
    train_target_inds = target_inds[training_inds]
    Y_train = Y[training_inds]
    
    test_drug_inds = drug_inds[test_inds]
    test_target_inds = target_inds[test_inds]
    Y_test = Y[test_inds]
    
    validation_drug_inds = drug_inds[validation_inds]
    validation_target_inds = target_inds[validation_inds]
    Y_validation = Y[validation_inds]
    
    # The setting number is not explicitly given so it needs to be found out.
    if set(test_drug_inds).issubset(set(train_drug_inds)):
        if set(test_target_inds).issubset(set(train_target_inds)):
            setting = "S1"
        else:
            setting = "S2"
    else:
        if set(test_target_inds).issubset(set(train_target_inds)):
            setting = "S3"
        else:
            setting = "S4"

    # Create the Kronecker product kernel matrices for drugs and targets.
    early = 1
    kro = 2
    kernel12 = GaussianKernel(XD, gamma=1e-5)
    K12 = kernel12.getKM(XD)
    kernel22 = GaussianKernel(XT, gamma=1e-5)
    K22 = kernel22.getKM(XT)
    if early==1:
        K12_vali = kernel12.getKM(XD)
        K22_vali = kernel22.getKM(XT)
    K12_test = kernel12.getKM(XD)
    K22_test = kernel22.getKM(XT)
    if kro==1:
        krotu = naive_kronecker_kernel(K12, K22, train_drug_inds, train_target_inds, train_drug_inds, train_target_inds)
    else:
        krotu = mv_kronecker(K12, K22, train_drug_inds, train_target_inds, train_drug_inds, train_target_inds)
    n = len(train_drug_inds)

    # Collect the information of the current case in a way that they can be collected as a data frame in the end.
    gathered_infos = pd.Series(data = {'data_set':data_name, 'setting':setting, 'binary':binary, 'max_iterations':MI, 'train_size': len(Y_train), 'test_size':len(Y_test), 'validation_size':len(Y_validation), 'random_seed':random_seed})
    
    """
    KronRLS hyperparameter optimization and reference model.
    """
    start_rls = time.clock()
    kronRLS_performances = []
    kronRLS_MIs = []
    kronRLS_predictors = []
    
    # Hyperparameter optimization.
    for rls_lambda in lambdas:
        cb_validation = CallBack(K12_vali, K22_vali, Y_validation, validation_drug_inds, validation_target_inds, ESlag = 100)
        try:
            CGKronRLS(K1 = K12, K2 = K22, Y=Y_train, label_row_inds = train_drug_inds, label_col_inds = train_target_inds, regparam = rls_lambda, callback = cb_validation, maxiter = MI)
        except ValueError as err:
            kronRLS_MIs.append(err.args[1])
            kronRLS_performances.append(err.args[2])
            kronRLS_predictors.append(err.args[3])
        else:
            kronRLS_MIs.append(cb_validation.best_iter)
            kronRLS_performances.append(cb_validation.maxperf)
            kronRLS_predictors.append(cb_validation.best_learner)
    # Details of the best KronRLS model.
    best_lambda_index = np.argmax(kronRLS_performances)
    rlslam = lambdas[best_lambda_index]
    gathered_infos['lambda_KronRLS'] = rlslam

    # Sparsity and generalization performance of the best KronRLS model.
    learner = kronRLS_predictors[best_lambda_index]
    rlsx = learner.A 
    rlsalpha = np.abs(rlsx)
    rlslas = len(rlsalpha[rlsalpha>zero_limit])
    gathered_infos['nz_KronRLS'] = rlslas
    rlsP = learner.predict(K12_test, K22_test, test_drug_inds, test_target_inds)
    perf_rls_ci = cindex(Y_test, rlsP)
    gathered_infos['CI_KronRLS'] = perf_rls_ci
    criteria_value = np.abs(n-rlslas)/(n-2) + np.abs(0.5-perf_rls_ci)/(0.5)
    gathered_infos['criteria_KronRLS'] = criteria_value

    # Print the runtime of finding the best KronRLS model.
    end_rls = time.clock()
    kesto_rls = end_rls - start_rls
    print("Data %s, setting %s, binary %d, MI %d, Time for finding the optimal dual variable vector for KronRLS %f seconds" % (data_name, setting, int(binary), MI, kesto_rls))
    
    # Save the obtained solution in a text file so that it can be utilized as a starting point in a file LMBMKronL0LS.py.
    refeta = rlsx.tolist()
    if(binary):
        np.savetxt('referenceKronRLS_'+data_name+'_'+setting+'_binary_'+str(random_seed)+'.txt', refeta, delimiter=',')
    else:
        np.savetxt('referenceKronRLS_'+data_name+'_'+setting+'_'+str(random_seed)+'.txt', refeta, delimiter=',')

    """
    KronSVM only for binary labels.
    """
    if(binary):
        print("KronSVM started")
        start_svm = time.clock()
        kronSVM_performances = []
        kronSVM_MIs = []
        kronSVM_predictors = []
        
        # Hyperparameter optimization.
        for svm_lambda in lambdas:
            cb_validation = CallBack(K12_vali, K22_vali, Y_validation, validation_drug_inds, validation_target_inds, ESlag = 100)
            try:
                KronSVM(K1 = K12, K2 = K22, Y=Y_train, label_row_inds = train_drug_inds, label_col_inds = train_target_inds, regparam = svm_lambda, callback = cb_validation, maxiter = MI)
            except ValueError as err:
                kronSVM_MIs.append(err.args[1])
                kronSVM_performances.append(err.args[2])
                kronSVM_predictors.append(err.args[3])
            else:
                kronSVM_MIs.append(cb_validation.best_iter)
                kronSVM_performances.append(cb_validation.maxperf)
                kronSVM_predictors.append(cb_validation.best_learner)
        # Details of the best KronSVM model.
        best_lambda_index = np.argmax(kronSVM_performances)
        svmlam = lambdas[best_lambda_index]
        gathered_infos['lambda_KronSVM'] = svmlam
        
        # Sparsity and generalization performance of the best KronSVM model.
        learner = kronSVM_predictors[best_lambda_index]
        svmx = learner.A 
        svmalpha = np.abs(svmx)
        svmlas = len(svmalpha[svmalpha>zero_limit])
        gathered_infos['nz_KronSVM'] = svmlas
        svmP = learner.predict(K12_test, K22_test, test_drug_inds, test_target_inds)
        perf_svm_ci = cindex(Y_test, svmP)
        gathered_infos['CI_KronSVM'] = perf_svm_ci
        criteria_value = np.abs(n-svmlas)/(n-2) + np.abs(0.5-perf_svm_ci)/(0.5)
        gathered_infos['criteria_KronSVM'] = criteria_value
        
        # Print the runtime of finding the best KronSVM model.
        end_svm = time.clock()
        kesto_svm = end_svm - start_svm
        print("Data %s, setting %s, binary %d, MI %d, Time for finding the optimal dual variable vector for KronSVM %f seconds" % (data_name, setting, int(binary), MI, kesto_svm))
    
    return(gathered_infos)

if __name__ == "__main__":
    random_seeds = [2688385916]
    split_percentage = 0.6
    datasets = ["davis", "metz", "merget"]
    df_list = []
    for ds in datasets:
        for random_seed in random_seeds:
            XD, XT, Y, drug_inds, target_inds, Y_binary = eval('load_split_data.load_'+ds+'()')
            df_indices, training_indices, S1_test_validation, S2_test_validation, S3_test_validation, S4_test_validation = load_split_data.splits(drug_inds, target_inds, split_percentage, random_seed)
            Ys = [Y_binary, Y]
            all_test_validation = [S1_test_validation, S2_test_validation, S3_test_validation, S4_test_validation]
            lambdas = [2.0**(-10), 2.0**(-5), 2.0**(-4), 2.0**(-3), 2.0**(-2), 2.0**(-1), 2.0**(0), 2.0**(1), 2.0**(2), 2.0**(3), 2.0**(4), 2.0**(5), 2.0**(10)]
            
            maxiter = 1000
            parameters = it.product(Ys, [drug_inds], [target_inds], [training_indices], all_test_validation, [lambdas], [maxiter], [XD], [XT], [ds], [random_seed])
            # Compute in parallel.
            pool = mp.Pool(processes = 4)
            output = pool.map(referenceModels, list(parameters))
            pool.close()
            pool.join()
            # Collect the informations of the best KronRLS and KronSVM models in a data frame.
            df = pd.concat(output, ignore_index = False, axis = 1)
            df_list.append(df.transpose())
    # Concatenate the data frames of each data set and save all the information in a csv-file.
    pd.concat(df_list, ignore_index = True).to_csv('./KronRLS_KronSVM_nz_ci_criteria.csv', index = False)