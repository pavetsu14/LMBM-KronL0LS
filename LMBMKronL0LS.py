from rlscore.kernel import GaussianKernel
import warnings
warnings.filterwarnings("ignore")
import lmbmprogram
import load_split_data

import numpy as np
from numpy import linalg
from rlscore.utilities.sampled_kronecker_products import sampled_vec_trick
import time
from heapq import nlargest

from rlscore.measure import cindex
import multiprocessing as mp
import itertools as it
np.set_printoptions(threshold=np.inf)
import pandas as pd

"""
CHOOSE THE PARAMETERS TO MATCH HOW YOU WANT TO RUN THE METHOD!

loss: 
    (0) Do not use any constant coefficient in front of the loss function that is part of the objective function.
    (1) Use 1/2 as a coefficient in front of the loss function that is part of the objective function.

early:
    (0) Do not use early stopping when running LMBM-kronL0LS.
    (1) Use early stopping procedure when running LMBM-kronL0LS.

kro:
    (1) Use explicit Kronecker product kernel matrix.
    (2) Use implicit Kronecker product kernel matrix via generalized vector trick.

main_objective:
    (1) Maintain satisfactory prediction accuracy with respect to the non-sparse reference solution.
    (2) Treat the attempt for sparsity as an equally important goal as the attempt for obtaining a good predictor. 

start_procedure:
    (0) Every optimization problem is started from a vector of ones. 
    (1) For every k, the optimization problem related to the first value of rho is started from a vector of ones. 
        The previous solution is used as a starting point of the optimization problem related to the next value of rho.
    (2) Always inherit the previous solution as a starting point of the next optimization problem.

Other parameters.

random_seed:
    Choose a random seed so that the randomness in splitting the data is controlled and the calculations are reproducable. 

zero_limit:
    A limit value that is used when it is defined which values in the dual vector are close enough to zero.
"""

def run_lmbmkronl0ls(params):
    Y = params[0]
    binary = int(len(set(Y)) == 2)
    drug_inds = params[1]
    target_inds = params[2]
    training_inds = params[3]
    test_inds = params[4][0]
    validation_inds = params[4][1]
    XD = params[5]
    XT = params[6]
    data_name = params[7]
    random_seed = params[8]
    main_objective = params[9]
    start_procedure = params[10]
    decrease = params[11]
    rolista = params[12]
    ESP_lag = params[13]
    loss = 1
    early = 1
    kro = 2
    zero_limit = 1e-4
    
    
    train_drug_inds = drug_inds[training_inds]
    train_target_inds = target_inds[training_inds]
    Y_train = Y[training_inds]
    
    test_drug_inds = drug_inds[test_inds]
    test_target_inds = target_inds[test_inds]
    Y_test = Y[test_inds]
    
    validation_drug_inds = drug_inds[validation_inds]
    validation_target_inds = target_inds[validation_inds]
    Y_validation = Y[validation_inds]
    
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
    """
    Functions for obtaining the kronecker kernel matrix.

    naive_kronecker_kernel:
        Input: the kernel matrices and row and column indices.
        Output: explicit kronecker kernel matrix?

    mv_kronecker:
        Input: the kernel matrices and row and column indices. 
        Output: implicit kronecker kernel matrix?
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
    Functions that are used by LMBM algorithm.

    func:
        Returns the value of the loss function for given k and rho values?

    funcL0:
        Returns the L0 norm of the solution.

    subgra:
        Updates the values in vector g.
    """
    def func(x,n):
        # Explicit Kronecker product kernel matrix
        if kro==1:
            if loss==0:
                lossosa = linalg.norm(np.dot(Y_train - np.matmul(krotu, x), Y_train - np.matmul(krotu, x)))
            else:
                lossosa = (1.0/2.0) * (linalg.norm(np.dot(Y_train - np.matmul(krotu, x), Y_train - np.matmul(krotu, x))))
        # GVT
        else:
            if loss==0:
                lossosa = linalg.norm(np.dot(Y_train - krotu(x), Y_train - krotu(x)))
            else:
                # print("Y type:", type(Y_train[0]), " krotu type:", type(krotu(x)[0]))
                lossosa = (1.0/2.0) * (linalg.norm(np.dot(Y_train - krotu(x), Y_train - krotu(x))))
        absx = np.abs(x)
        # L0-norm
        kla = nlargest(kno, absx)
        regosa = ro * np.sum(absx) - ro * np.sum(kla)
        fuarvo = lossosa + regosa
        return fuarvo


    def funcL0(x):
        # L0-norm
        absaL0 = np.abs(x)
        klaL0 = nlargest(kno, absaL0)
        L0arvo = ro * np.sum(absaL0) - ro * np.sum(klaL0)
        return L0arvo


    def subgra(x,g,n):
        # Explicit Kronecker product kernel matrix
        if kro==1:
            if loss==0:
                lossvek = (-2.0) * (np.matmul(krotu, Y_train - np.matmul(krotu, x)))
            else:
                lossvek = (-1.0) * (np.matmul(krotu, Y_train - np.matmul(krotu, x)))
        # GVT
        else:
            if loss==0:
                lossvek = (-2.0) * (krotu(Y_train - krotu(x)))
            else:
                lossvek = (-1.0) * (krotu(Y_train - krotu(x)))
        # L0-norm
        gl1 = np.zeros(n)
        for i in range(n):
            if x[i] < 0.0:
                gl1[i] = -1.0
            else:
                gl1[i] = 1.0
        # k-pseudonorms
        absx = np.abs(x)
        klai = np.argpartition(absx, -kno)[-kno:]
        glk = np.zeros(n)
        for i in range(len(klai)):
            temp = klai[i]
            if x[temp] < 0.0:
                glk[temp] = -1.0
            else:
                glk[temp] = 1.0
        regvek = ro * (gl1 - glk)
        g_new = lossvek + regvek
        for i in range(n):
            g[i] = g_new[i]
        return


    def vali(x, n):
        if early==0:
            lopa = 6
        else:
            p_validation = sampled_vec_trick(x, K22_vali, K12_vali, validation_target_inds, validation_drug_inds, train_target_inds, train_drug_inds)
            aucoma = cindex(Y_validation, p_validation)
            itsx = np.abs(x)
            harx = len(itsx[itsx>zero_limit])
            # Bi-objective vector that contains the values that are used for calculating the criteria value.
            bov = np.array([harx, aucoma])
            abs_diff_from_worst = np.abs(bov_worst - bov)
            abs_diff_from_worst_nz = float(abs_diff_from_worst[0]) / float(n-2)
            abs_diff_from_worst_ci = float(abs_diff_from_worst[1]) / 0.5

            cri_last = weight_nz * abs_diff_from_worst_nz + weight_ci * abs_diff_from_worst_ci
            x_latest = x.copy()
            # Save the dual vector x_latest in a list listax.
            listax.append(x_latest)
            # Save the criteria value in a list listask.
            listask.append(cri_last)
            # Lenght of the list of criteria values, i.e. number of iterations.
            n_iterations = np.shape(listask)[0]
            # IF you want, you can replace with any number of iteration rounds you wish to force the code to run before starting to observe the progress of the results for validation data
            # Iterate at least 50 times before even thinking about early stopping.
            n_skip = 50
            if n_iterations<n_skip:
                lopa = 0
            # First iteration when the criteria value and the dual vector are saved for later comparisons.
            elif n_iterations==n_skip:
                lopa = 0
                # Save the latest criteria value in variable cri_sol_validation.
                cri_sol_validation = listask[-1]
                # Save the latest dual vector in variable x_sol.
                x_sol = listax[-1]
            # Rest of the iterations go through this part.
            else:
                cri_sol_validation = listask[n_iterations-2-len(no_improvement)]
                x_sol = listax[n_iterations-2-len(no_improvement)]
                # If there is no improvement in the criteria value.
                if cri_sol_validation > cri_last:
                    no_improvement.append("1")
                    # No improvement in the last ESP_lag iterations.
                    if len(no_improvement) == ESP_lag:
                        lopa = 1
                        # Add the best dual vector to a list earlyli.
                        earlyli.append(x_sol)
                    else:
                        lopa = 0
                # Current solution is better than the previous best solution.
                else:
                    # Update the variables that contain the information of the best solution.
                    cri_sol_validation = listask[-1]
                    x_sol = listax[-1]
                    lopa = 0
                    # Initialize the count back to zero
                    del no_improvement[:]
        return lopa
    
    """
    Create Kronecker product kernel matrix.
    """
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
    # The worst possible bi-objective vector (is always the same, does not depend on the reviewed data set).
    bov_worst = np.array([n, 0.5])
    # There is an option to adjust the importance of meeting the objectives with respect to each other.
    # The larger weight the more it affects to the criteria value and thus the more important it is to get the associated
    # objective value (number of nonzero elements / prediction accuracy) as far away as possible from its worst bi-objective vector value.
    # Note that the effects of tuning these weights has not been throughoutly examined and can result in unexpected results.
    weight_nz = 1.0
    weight_ci = 1.0
    """
    The non-sparse reference solution
        The solution is previously found by running file referenceModels.py.
        Calculate the number of non-zero elements in the reference solution.
        Set the small enough values to zero. 
        Use function sampled_vec_trick to obtain test set predictions with the reference solution and calculate a c-index value.
        Calculate a criteria value for the reference solution, 
            i.e. a weighted difference between the reference solution and the worst possible solution proportional to the difference between the best and the worst possible solutions.
            Criteria value = weight of sparsity * relative difference of sparsity + weight of c-index * relative difference of c-indices.
    """
    if(binary):
        file_ref = open('referenceKronRLS_'+data_name+'_'+setting+'_binary_'+str(random_seed)+'.txt', "r")
    else:
        file_ref = open('referenceKronRLS_'+data_name+'_'+setting+'_'+str(random_seed)+'.txt', "r")

    refte = file_ref.readlines()
    x_ref = np.array(refte, dtype='f8')

    """
    Calculate the criteria value for the non-sparse reference solution.
    """
    x_ref_abs = np.abs(x_ref)
    nz_ref = len(x_ref_abs[x_ref_abs>zero_limit])
    x_ref_sparse = x_ref.copy()

    for k in range(n):
        if np.abs(x_ref[k])<=zero_limit:
            x_ref_sparse[k] = 0.0

    P_test_ref_sparse = sampled_vec_trick(x_ref_sparse, K22_test, K12_test, test_target_inds, test_drug_inds, train_target_inds, train_drug_inds)
    ci_ref_sparse = cindex(Y_test, P_test_ref_sparse)
    bov_ref = np.array([nz_ref, ci_ref_sparse])

    diff_from_worst_ref = bov_worst - bov_ref
    abs_diff_from_worst_ref = np.abs(diff_from_worst_ref)
    # The best possible amount of nonzero dual variables = 2, the worst = n
    abs_diff_from_worst_ref_nz = float(abs_diff_from_worst_ref[0]) / float(n-2)
    # The choice of the data set has an effect to the maximal prediction accuracy (the main reason is that with checkerboard
    # data we have an option to add noise to the data
    # In any case, the worst possible prediction accuracy = 0.5 (totally random, no signal to any direction)
    abs_diff_from_worst_ref_ci = float(abs_diff_from_worst_ref[1]) / 0.5

    cri_start = weight_nz * abs_diff_from_worst_ref_nz + weight_ci * abs_diff_from_worst_ref_ci
    # Define the maximum value of k to be the number of nonzero elements in the reference solution.
    k_lim = nz_ref
    print("%s, %s, binary = %r: Reference solution [# nz, CI, criteria] = [%i, %f, %f]" %(data_name, setting, binary, bov_ref[0], bov_ref[1], cri_start))

    """
    The usage of LMBM-KronL0LS algorithm.

    rolista:
        List of penalization weights.

    aucli, harli, muli:
        Initialize empty list where the results for different k and rho values will be saved.

    k_step:
        Step size that is used when increasing the value of k.
        The step size is defined to be 5 % of the training sample size.

    ci_limit:
        The maximum accepted decrease in the value of c-index in the case of MO1.
    """

    print("Number of training samples %f" %n)

    ci_values = [bov_ref[1]]
    nz_values = [bov_ref[0]]
    regularization_coefficients = [0]
    k_list = [n]
    update_info_list = ["ref"]
    criteria_values = [cri_start]
    # Define the step size.
    k_step = int(0.05*n)
    no_updates_done = True
    
    # Initialize the first starting point according to the main objective.
    if main_objective == 1:
        a_start = x_ref.copy()
    else:
        a_start = np.ones(n)

    start_rok = time.clock()
    # Change the rolista and iterate until at least one solution has had a better criteria_value than the reference solution. 
    while(no_updates_done):
        rolkm = len(rolista)
        ci_limit = (1-decrease)*ci_ref_sparse

        for k in range(2, n+2, k_step):
            if k>n:
                kno = n-1
            else:
                kno = k
            if kno > k_lim:
                print("%s, %s, MO%i, SP%i, binary = %r: Run onto the upper limit of k value!" % (data_name, setting, main_objective, start_procedure, binary))
                break
            
            for j in range(len(rolista)):
                if kno > k_lim:
                    print("%s, %s, MO%i, SP%i, binary = %r: Run onto the upper limit of k value!" % (data_name, setting, main_objective, start_procedure, binary))
                    break
                else:
                    ro = rolista[j]
                    earlyli = list()
                    listax = list()
                    listask = list()
                    no_improvement = []
                    # If we wish to use the newest solution point as a starting point of the next optimization problem always when such is available
                    if start_procedure==1:
                        # If we are starting to solve the optimization problem related to the first k&rho pair
                        if k==2 and j==0:
                            x = a_start.copy()
                        else:
                            x = xsol.copy()
                    else:
                        if j==0:
                            if k == 2:
                                x = a_start.copy()
                            else:
                                x = a_small.copy()
                        else:
                            x = xsol.copy()
                    # Initialization of the variables LMBM implementation needs as input.
                    fu = 1.0
                    lopa = 3
                    g = np.ones(n)
                    mc = 7
                    print("LMBM algorithm returns:", lmbmprogram.lmbm_mod.lmbm(x,fu,func,lopa,vali,g,subgra,mc,n))
                    # The best solution returned by LMBM was obtained either by early stopping or is the latest solution.
                    if early==1 and len(earlyli)>0:
                        xsol = earlyli[0]
                    else:
                        xsol = x.copy()
                    x_last = x.copy()
                    absx = np.abs(xsol)
                    n_nz = len(absx[absx>zero_limit])
                    x_sparse = xsol.copy()
                    for i in range(n):
                        if np.abs(xsol[i])<=zero_limit:
                            x_sparse[i] = 0.0
                    P_test_sparse = sampled_vec_trick(x_sparse, K22_test, K12_test, test_target_inds, test_drug_inds, train_target_inds, train_drug_inds)
                    ci_sparse = cindex(Y_test, P_test_sparse)
                    print("%s, %s, MO%i, SP%i, binary = %r: The solution for the optimization problem attached to the current k = %i and rho = %f value is CI = %f, # nz = %i" % (data_name, setting, main_objective, start_procedure, binary, kno, ro, ci_sparse, n_nz))

                    if main_objective == 1:
                        if start_procedure == 2 and k == 2:
                            a_small = x_last.copy()
                        # Test whether the obtained solution is good enough to be saved into the solution lists or whether we should
                        # break and move on to the next optimization problem associated with a new k&rho pair.
                        if ci_sparse < ci_limit:
                            print("The obtained prediction accuracy is too low.")
                            break
                    else:
                        if start_procedure == 2 and k == 2 and j == rolkm-1:
                            a_small = x_last.copy()
                            
                    ci_values.append(ci_sparse)
                    nz_values.append(n_nz)
                    regularization_coefficients.append(ro)
                    k_list.append(k)
                    
                    # Test whether the obtained solution gives us a reason to tighten the upper limit of the traversed k values
                    bov_sol = np.array([n_nz, ci_sparse])
                    diff_from_worst_sol = bov_worst - bov_sol
                    abs_diff_from_worst_sol = np.abs(diff_from_worst_sol)
                    abs_diff_from_worst_sol_nz = float(abs_diff_from_worst_sol[0]) / float(n-2)
                    abs_diff_from_worst_sol_ci = float(abs_diff_from_worst_sol[1]) / 0.5
                    
                    cri_sol = weight_nz * abs_diff_from_worst_sol_nz + weight_ci * abs_diff_from_worst_sol_ci
                    criteria_values.append(cri_sol)
                    # Current solution produces a better criteria value than the reference solution
                    if (cri_sol>=cri_start):
                        # Update the value of k_lim to be the number of non-zero elements in the current solution. 
                        k_lim = n_nz
                        # Update the criteria value to match the best solution found so far
                        cri_start = cri_sol
                        print("%s, %s, MO%i, SP%i, binary = %r: New upper limit for k is %i" %(data_name, setting, main_objective, start_procedure, binary, k_lim))
                        k_ref_updated = k
                        rho_ref_updated = ro
                        no_updates_done = False
                        update_info_list.append("yes")
                    else:
                        update_info_list.append("no")
                    # Let's test whether it is still relevant to continue aka if L0-norm still has a nonzero value
                    normiLMBM = funcL0(x_sparse)
                    if normiLMBM<1e-1:
                        print("No point to continue increasing rho value since the constraint related to the maximum amount of nonzero elements (k kpl) already satisfied!")
                        break

        # Allow a larger decrease in C-index if the method has not yet been able to find any good enough solutions. 
        if(no_updates_done):
            decrease += 0.01
            print("Accepted decrease increased to %f" % decrease)

    end_rok = time.clock()
    print("%s, %s, binary = %r: The amount of time it took to scan through the relevant k&rho value pair space: %f seconds" %(data_name, setting, binary, end_rok - start_rok))
    print("%s, %s, binary = %r: The highest criteria value was obtained with k = %d and rho = %f" %(data_name, setting, binary, k_ref_updated, rho_ref_updated))
    results = pd.DataFrame({'ci':ci_values, 'nz':nz_values, 'rho':regularization_coefficients, 'k':k_list, 'updated_ref':update_info_list, 'criteria_value':criteria_values, 'setting':setting, 'MO':main_objective, 'SP':start_procedure, 'accepted_decrease':decrease, 'binary':binary})
    return(results)

if __name__ == "__main__":
    random_seeds = [2688385916]
    # Choose the percentage of drugs and targets that are defined to be in area S1 and the percentage of area S1 pairs that will be used as training data.
    split_percentage = 0.6
    # Choose the accepted decrease in C-index in the case of MO1.
    decrease = 0.05
    datasets = ["davis", "metz", "merget"]
    main_objectives = [1,2]
    start_procedures = [1,2]
    # The values that the regularization parameter rho can get.
    rho_values = [0.1, 0.5, 1.0, 2.0, 4.0, 12.0, 50.0, 500.0]
    ESP_lag = 200
    for random_seed in random_seeds:
        for ds in datasets:
            XD, XT, Y, drug_inds, target_inds, Y_binary = eval('load_split_data.load_'+ds+'()')
            df_indices, training_indices, S1_test_validation, S2_test_validation, S3_test_validation, S4_test_validation = load_split_data.splits(drug_inds, target_inds, split_percentage, random_seed)
            Ys = [Y_binary, Y]
            all_test_validation = [S1_test_validation, S2_test_validation, S3_test_validation, S4_test_validation]
            
            parameters = it.product(Ys, [drug_inds], [target_inds], [training_indices], all_test_validation, [XD], [XT], [ds], [random_seed], main_objectives, start_procedures, [decrease], [rho_values], [ESP_lag])
            # Compute in parallel.
            pool = mp.Pool(processes = 8)
            output = pool.map(run_lmbmkronl0ls, list(parameters))
            pool.close()
            pool.join()
            pd.concat(output, ignore_index = True).to_csv('LMBMresultsRefKronRLS_'+ds+'_'+str(random_seed)+'_decrease'+str(decrease)+'.csv', index = False)