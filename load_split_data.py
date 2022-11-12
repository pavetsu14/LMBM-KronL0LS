"""
A file containing the functions for importing the data sets and creating the data splits for the 4 different settings.
"""
import numpy as np
import pandas as pd

"""
Function to load an incomplete data set introduced by Metz et al. (2011).
Returns the data matrices and lists of drug and target indices for the known pairs.
"""
def load_metz(threshold = 7.6):
    Y = np.loadtxt("known_drug-target_interaction_affinities_pKi__Metz_et_al.2011.txt")
    XD = np.loadtxt("drug-drug_similarities_2D__Metz_et_al.2011.txt")
    XT = np.loadtxt("target-target_similarities_WS_normalized__Metz_et_al.2011.txt")
    drug_inds, target_inds = np.where(np.isnan(Y)==False)
    Y = Y[drug_inds, target_inds]
    Y_binary = (-1.0)**(Y < threshold)
    return XD, XT, Y, drug_inds.astype(int), target_inds.astype(int), Y_binary

"""    
Function to load a complete data set introduced by Davis et al. (2011).
Returns the data matrices and lists of drug and target indices for the known pairs.
The matrix of drug similarities is multiplied by 100 in order to obtain the same range as in the corresponding matrix of Metz data.
The returned continuous labels are natural logarithm of the Kd values so that the range is again similar to the range of continuous labels in Metz data.
"""
def load_davis(threshold = 30):
    Y = np.loadtxt("drug-target_interaction_affinities_Kd__Davis_et_al.2011.txt")
    XD = np.loadtxt("drug-drug_similarities_2D__Davis_et_al.2011.txt")
    XD = XD*100
    XT = np.loadtxt("target-target_similarities_WS_normalized__Davis_et_al.2011.txt")
    drug_inds, target_inds = np.where(np.isnan(Y)==False)
    Y = Y[drug_inds, target_inds]
    Y_binary = (-1.0)**(Y > threshold)
    Y_log = np.log(Y)
    return XD, XT, Y_log, drug_inds.astype(int), target_inds.astype(int), Y_binary

"""
Function to load an incomplete data set introduced by Merget et al. (2017) and updated by Cichonska et al (2018).
Returns the data matrices and lists of drug and target indices for the known pairs.
The matrices of drug similarities and target similaritied are multiplied by 100 in order to 
obtain the same range as in the corresponding matrices of Metz data.
"""
def load_merget(threshold = 7):
    Y = np.loadtxt("Merget_DTIs_2967com_226kin.txt")
    XD = np.loadtxt("Kd_Tanimoto-shortestpath.txt")
    XD = XD*100
    XT = np.loadtxt("Kp_GS-ATP_L5_Sp4.0_Sc4.0.txt")
    XT = XT*100
    drug_inds, target_inds = np.where(np.isnan(Y)==False)
    Y = Y[drug_inds, target_inds]
    Y_binary = (-1.0)**(Y <= threshold)
    return XD, XT, Y, drug_inds.astype(int), target_inds.astype(int), Y_binary
"""
Function to create the splits for the four settings.
S1: Both the drugs and targets in the test/validation pairs are also included in the training data.
S2: The drugs in the test/validation pairs are included in the training data, but the targets are new.
S3: The targets in the test/validation pairs are included in the training data, but the drugs are new.
S4: Neither the drugs nor the targets in the test/validation pairs are included in the training data.
"""
def splits(drug_inds, target_inds, split_percentage, random_seed):
    n_sample = len(drug_inds)
    drugs = list(set(drug_inds))
    n_drugs = len(drugs)
    targets = list(set(target_inds))
    n_targets = len(targets)

    np.random.seed(random_seed)
    # Shuffle the unique drugs and split them so that a wanted percentage of them are considered as area S1 drugs, 
    # rest of them are split in halves into test and validation drugs.
    np.random.shuffle(drugs)
    index_drugs_S1 = int(np.ceil(n_drugs*split_percentage))
    drugs_S1 = drugs[:index_drugs_S1]
    index_drugs_test = index_drugs_S1 + int(np.ceil((n_drugs-index_drugs_S1)*0.5))
    drugs_test = drugs[index_drugs_S1:index_drugs_test]
    drugs_validation = drugs[index_drugs_test:]

    # Shuffle the unique targets and split them so that a wanted percentage of them are considered as area S1 targets, 
    # rest of them are split in halves into test and validation targets.
    np.random.shuffle(targets)
    index_targets_S1 = int(np.ceil(n_targets*split_percentage))
    targets_S1 = targets[:index_targets_S1]
    index_targets_test = index_targets_S1 + int(np.ceil((n_targets-index_targets_S1)*0.5))
    targets_test = targets[index_targets_S1:index_targets_test]
    targets_validation = targets[index_targets_test:]
    
    S1_indices = []
    # Find the pairs that are in the area S1.
    for i in range(n_sample):
        if drug_inds[i] in drugs_S1 and target_inds[i] in targets_S1:
            S1_indices.append(i)
    
    # Split the set of pair indices in the area S1 into training, S1 test and S2 validation sets.
    np.random.shuffle(S1_indices)
    index_training = int(np.ceil(len(S1_indices)*split_percentage))
    training_indices = S1_indices[:index_training]
    index_S1_test = index_training + int(np.ceil((len(S1_indices)-index_training)*0.5))
    S1_test = S1_indices[index_training:index_S1_test]
    S1_validation = S1_indices[index_S1_test:]
    drugs_training = list(set(drug_inds[training_indices]))
    targets_training = list(set(target_inds[training_indices]))
    
    # Create the test and validation sets of indices for every setting.
    S1_test_indices = []
    S1_validation_indices = []
    S2_test_indices = []
    S2_validation_indices = []
    S3_test_indices = []
    S3_validation_indices = []
    S4_test_indices = []
    S4_validation_indices = []
    for i in range(n_sample):
        if drug_inds[i] in drugs_training:
            if target_inds[i] in targets_training:
                if i in S1_test:
                    S1_test_indices.append(i)
                elif i in S1_validation:
                    S1_validation_indices.append(i)
            elif target_inds[i] in targets_test:
                S2_test_indices.append(i)
            elif target_inds[i] in targets_validation:
                S2_validation_indices.append(i)
        elif drug_inds[i] in drugs_test: 
            if target_inds[i] in targets_training:
                S3_test_indices.append(i)
            elif target_inds[i] in targets_test:
                S4_test_indices.append(i)
        elif drug_inds[i] in drugs_validation:
            if target_inds[i] in targets_training:
                S3_validation_indices.append(i)
            elif target_inds[i] in targets_validation:
                S4_validation_indices.append(i)
    
    # Create a data frame that contains all the information of the created splits for saving purposes.
    df_indices = pd.concat([pd.DataFrame({'roundID':random_seed, 'subset':'training', 'setting':'all', 'index':training_indices, 'drugs':drug_inds[training_indices], 'targets':target_inds[training_indices]}), 
    pd.DataFrame({'roundID':random_seed, 'subset':'test', 'setting':'S1', 'index':S1_test_indices, 'drugs':drug_inds[S1_test_indices], 'targets':target_inds[S1_test_indices]}), 
    pd.DataFrame({'roundID':random_seed, 'subset':'validation', 'setting':'S1', 'index':S1_validation_indices, 'drugs':drug_inds[S1_validation_indices], 'targets':target_inds[S1_validation_indices]}), 
    pd.DataFrame({'roundID':random_seed, 'subset':'test', 'setting':'S2', 'index':S2_test_indices, 'drugs':drug_inds[S2_test_indices], 'targets':target_inds[S2_test_indices]}), 
    pd.DataFrame({'roundID':random_seed, 'subset':'validation', 'setting':'S2', 'index':S2_validation_indices, 'drugs':drug_inds[S2_validation_indices], 'targets':target_inds[S2_validation_indices]}), 
    pd.DataFrame({'roundID':random_seed, 'subset':'test', 'setting':'S3', 'index':S3_test_indices, 'drugs':drug_inds[S3_test_indices], 'targets':target_inds[S3_test_indices]}), 
    pd.DataFrame({'roundID':random_seed, 'subset':'validation', 'setting':'S3', 'index':S3_validation_indices, 'drugs':drug_inds[S3_validation_indices], 'targets':target_inds[S3_validation_indices]}), 
    pd.DataFrame({'roundID':random_seed, 'subset':'test', 'setting':'S4', 'index':S4_test_indices, 'drugs':drug_inds[S4_test_indices], 'targets':target_inds[S4_test_indices]}), 
    pd.DataFrame({'roundID':random_seed, 'subset':'validation', 'setting':'S4', 'index':S4_validation_indices, 'drugs':drug_inds[S4_validation_indices], 'targets':target_inds[S4_validation_indices]})], 
    axis = 0, sort = False)

    # Combine the lists of test and validation indices for each setting.
    S1_test_validation = [S1_test_indices, S1_validation_indices]
    S2_test_validation = [S2_test_indices, S2_validation_indices]
    S3_test_validation = [S3_test_indices, S3_validation_indices]
    S4_test_validation = [S4_test_indices, S4_validation_indices]
    return(df_indices, training_indices, S1_test_validation, S2_test_validation, S3_test_validation, S4_test_validation)

if __name__ == "__main__":
    # Select a seed or multiple seeds for controlling the randomness in creating the splits. 
    random_seeds = [2688385916]
    # Choose the percentage of drugs and targets that are defined to be in area S1 and the percentage of area S1 pairs that will be used as training data.
    split_percentage = 0.6
    datasets = ["davis", "metz", "merget"] 
    for random_seed in random_seeds:
        for ds in datasets:
            # Load the data set in the wanted form.
            XD, XT, Y, drug_inds, target_inds, Y_binary = eval('load_'+ds+'()')
            # Split the data set into training data set and test and validation data sets for each setting.
            df_indices, training_indices, S1_test_validation, S2_test_validation, S3_test_validation, S4_test_validation = splits(drug_inds, target_inds, split_percentage, random_seed)
            # Save the information of the splits as a csv-file.
            df_indices.to_csv('splits_'+ds+'_RS_'+str(random_seed)+'.csv', index = False)