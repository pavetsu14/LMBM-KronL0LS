import numpy as np
from numpy import linalg
from sklearn.metrics import roc_auc_score
from rlscore.kernel import GaussianKernel
from rlscore.utilities.sampled_kronecker_products import sampled_vec_trick
import time
from heapq import nlargest

import lmbmprogram

from data_eplin import generate_chess

from rlscore.learner import CGKronRLS
#from rlscore.learner.kron_svm import KronSVM
#from rlscore.measure import cindex

np.set_printoptions(threshold=np.inf)




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


def func(x,n):
	# Explicit Kronecker product kernel matrix
	if kro==1:
		if loss==0:
			lossosa = linalg.norm(np.dot(Y2th - np.matmul(krotu, x), Y2th - np.matmul(krotu, x)))
		else:
			lossosa = (1.0/2.0) * (linalg.norm(np.dot(Y2th - np.matmul(krotu, x), Y2th - np.matmul(krotu, x))))
	# GVT
	else:
		if loss==0:
			lossosa = linalg.norm(np.dot(Y2th - krotu(x), Y2th - krotu(x)))
		else:
			lossosa = (1.0/2.0) * (linalg.norm(np.dot(Y2th - krotu(x), Y2th - krotu(x))))
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
			lossvek = (-2.0) * (np.matmul(krotu, Y2th - np.matmul(krotu, x)))
		else:
			lossvek = (-1.0) * (np.matmul(krotu, Y2th - np.matmul(krotu, x)))
	# GVT
	else:
		if loss==0:
			lossvek = (-2.0) * (krotu(Y2th - krotu(x)))
		else:
			lossvek = (-1.0) * (krotu(Y2th - krotu(x)))
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
		#ennuvali = sampled_vec_trick(x, K22_vali, K12_vali, tem10, tem9, tem2, tem1)
		# IF you are solving a binarized interaction affinity problem with Ki data, the labels will need to be converted according to a threshold value
		#Y2_valith = Y2_vali.copy()
		#Y2_valith[Y2_valith<7.6] = -1
		#Y2_valith[Y2_valith>=7.6] = 1
		#aucoma = roc_auc_score(Y2_valith, ennuvali)
		##aucoma = roc_auc_score(Y2_vali, ennuvali)
		#print("current prediction accuracy")
		#print(aucoma)
		lopa = 6
	else:
		ennuvali = sampled_vec_trick(x, K22_vali, K12_vali, tem10, tem9, tem2, tem1)
		# IF you are solving a binarized interaction affinity problem with Ki data, the labels will need to be converted according to a threshold value
		Y2_valith = Y2_vali.copy()
		Y2_valith[Y2_valith<7.6] = -1
		Y2_valith[Y2_valith>=7.6] = 1
		aucoma = roc_auc_score(Y2_valith, ennuvali)
		#aucoma = roc_auc_score(Y2_vali, ennuvali)
		itsx = np.abs(x)
		harx = len(itsx[itsx>1e-4])
		vatemp = np.array([harx, aucoma])
		vaero = huve - vatemp
		vapoik = np.abs(vaero)
		vapoikhar = float(vapoik[0]) / float(n-2)
		vapoikauc = float(vapoik[1]) / 0.5
		# IF you are adding noise to training data, extract its amount from 0.5!
		#vapoikauc = float(vapoik[1]) / 0.3
		skasu = harpaino * vapoikhar + aucpaino * vapoikauc
		xuusi = x.copy()
		listax.append(xuusi)
		listask.append(skasu)
		#print("validation results")
		#print(aucoma)
		#print(harx)
		listal = np.shape(listask)[0]
		# IF you want, you can replace with any number of iteration rounds you wish to force the code to run before starting to observe the progress of the results for validation data
		if listal<50:
			lopa = 0
		elif listal==50:
			lopa = 0
			par = listask[-1]
			parx = listax[-1]
		else:
			par = listask[listal-2-len(toisto)]
			parx = listax[listal-2-len(toisto)]
			if par > skasu:
				toisto.append('1')
				if len(toisto) == kipa:
					lopa = 1
					earlyli.append(parx)
				else:
					lopa = 0
			else:
				par = listask[-1]
				parx = listax[-1]
				lopa = 0
				del toisto[:]
	return lopa


# CHOOSE THE PARAMETERS TO MATCH HOW YOU WANT TO RUN THE METHOD!

# Do not use any constant as a coefficient in front of the loss function part of the objective function (0) or use 1/2 as a coefficient (1)
loss = 1
# Apply early stopping procedure (1) or skip it (0)
early = 1
# Use explicit Kronecker product kernel matrix (1) or the implicit version created via generalized vector trick (2)
kro = 2
# Solve a problem related to the Chessboard data (1) or Ki data (2)
datava = 2


if datava==1:
	# The amount of noise you wish to incorporate into the data
	# 0.2 = 20% of the labels will be flipped
	koh = 0.2
	# Parameter that affects to the complexity level of the Chessboard data via determining the "size" of the chessboard we derive the values from
	shakki = 5
	start_da = time.clock()
	# Generate the training data set
	X12, X22, rows2, cols2, Y2 = generate_chess(100, 100, shakki, koh, seed=1)
	kernel12 = GaussianKernel(X12, gamma=1)
	K12 = kernel12.getKM(X12)
	kernel22 = GaussianKernel(X22, gamma=1)
	K22 = kernel22.getKM(X22)
	# Generate the test data set
	X12_test, X22_test, rows2_test, cols2_test, Y2_test = generate_chess(500, 500, shakki, koh, seed=3)
	K12_test = kernel12.getKM(X12_test)
	K22_test = kernel22.getKM(X22_test)
	# Generate the validation data set
	# MAY NOT BE NECESSARY, depending whether you decided to apply early stopping procedure or not
	if early==1:
		X12_vali, X22_vali, rows2_vali, cols2_vali, Y2_vali = generate_chess(500, 500, shakki, koh, seed=2)
		K12_vali = kernel12.getKM(X12_vali)
		K22_vali = kernel22.getKM(X22_vali)
	if kro==1:
		krotu = naive_kronecker_kernel(K12, K22, rows2, cols2, rows2, cols2)
	else:
		krotu = mv_kronecker(K12, K22, rows2, cols2, rows2, cols2)
	n = len(rows2)
	end_da = time.clock()
	kesto_da = end_da - start_da
	print("The amount of time for generating the data sets and forming the Kronecker product kernel matrix ", kesto_da, " seconds")
elif datava==2:
	# THE EXPERIMENTAL SETTING TO USE: A (1), B (2), C (3) or D (4)
	setva = 4
	# SOLVE THE PROBLEM ALSO WITH KronRLS: yes (1) or no (0)
	krorls = 1
	# CALCULATE THE REFERENCE SOLUTION (1) OR READ IT DIRECTLY FROM THE TEXT FILE CREATED SEPARATELY (0)
	refva = 0
	def load_metz():
		Y = np.loadtxt("known_drug-target_interaction_affinities_pKi__Metz_et_al.2011.txt")
		XD = np.loadtxt("drug-drug_similarities_2D__Metz_et_al.2011.txt")
		XT = np.loadtxt("target-target_similarities_WS_normalized__Metz_et_al.2011.txt")
		drug_inds, target_inds = np.where(np.isnan(Y)==False)
		Y = Y[drug_inds, target_inds]
		print("XD dimensions %d %d" %XD.shape)
		print("XT dimensions %d %d" %XT.shape)
		print("Labeled pairs %d, all possible pairs %d" %(len(Y), XD.shape[0]*XT.shape[0]))
		return XD, XT, Y, drug_inds, target_inds
	def settingA_split():
		XD, XT, Y, drug_inds, target_inds = load_metz()
		np.random.seed(77)
		#random split to train/test, corresponds to setting A
		ind = list(range(len(Y)))
		np.random.shuffle(ind)
		train_ind = ind[:50000]
		test_ind = ind[50000:]
		train_drug_inds = drug_inds[train_ind]
		train_target_inds = target_inds[train_ind]
		Y_train = Y[train_ind]
		test_drug_inds = drug_inds[test_ind]
		test_target_inds = target_inds[test_ind]
		Y_test = Y[test_ind]
		return XD, XT, train_drug_inds, train_target_inds, Y_train, test_drug_inds, test_target_inds, Y_test
	def settingAvali_split():
		XD, XT, Y, drug_inds, target_inds = load_metz()
		np.random.seed(77)
		#random split to train/test, corresponds to setting A
		ind = list(range(len(Y)))
		np.random.shuffle(ind)
		train_ind = ind[:50000]
		kytkD=1
		vali_ind = []
		for i in range(len(ind)):
			if ind[i] not in train_ind:
				# If not in training set, split every other to test set and every other to validation set
				if kytkD==1:
					test_ind.append(i)
					kytkD=2
				else:
					vali_ind.append(i)
					kytkD=1
		train_drug_inds = drug_inds[train_ind]
		train_target_inds = target_inds[train_ind]
		Y_train = Y[train_ind]
		test_drug_inds = drug_inds[test_ind]
		test_target_inds = target_inds[test_ind]
		Y_test = Y[test_ind]
		vali_drug_inds = drug_inds[vali_ind]
		vali_target_inds = target_inds[vali_ind]
		Y_vali = Y[vali_ind]
		return XD, XT, train_drug_inds, train_target_inds, Y_train, vali_drug_inds, vali_target_inds, Y_vali, test_drug_inds, test_target_inds, Y_test
	def settingB_split():
		XD, XT, Y, drug_inds, target_inds = load_metz()
		np.random.seed(77)
		#random split to train/test, corresponds to setting B
		drows = list(range(XD.shape[0]))
		np.random.shuffle(drows)
		train_drows = set(drows[:800])
		train_ind = []
		test_ind = []
		for i in range(len(drug_inds)):
			if drug_inds[i] in train_drows:
				train_ind.append(i)
			else:
				test_ind.append(i)
		train_drug_inds = drug_inds[train_ind]
		train_target_inds = target_inds[train_ind]
		Y_train = Y[train_ind]
		test_drug_inds = drug_inds[test_ind]
		test_target_inds = target_inds[test_ind]
		Y_test = Y[test_ind]
		return XD, XT, train_drug_inds, train_target_inds, Y_train, test_drug_inds, test_target_inds, Y_test
	def settingBvali_split():
		XD, XT, Y, drug_inds, target_inds = load_metz()
		np.random.seed(77)
		#random split to train/test, corresponds to setting B
		drows = list(range(XD.shape[0]))
		np.random.shuffle(drows)
		train_drows = set(drows[:700])
		train_ind = []
		test_ind = []
		kytkD=1
		vali_ind = []
		for i in range(len(drug_inds)):
			if drug_inds[i] in train_drows:
				train_ind.append(i)
			else:
				# If not in training set, split every other to test set and every other to validation set
				if kytkD==1:
					test_ind.append(i)
					kytkD=2
				else:
					vali_ind.append(i)
					kytkD=1
		train_drug_inds = drug_inds[train_ind]
		train_target_inds = target_inds[train_ind]
		Y_train = Y[train_ind]
		test_drug_inds = drug_inds[test_ind]
		test_target_inds = target_inds[test_ind]
		Y_test = Y[test_ind]
		vali_drug_inds = drug_inds[vali_ind]
		vali_target_inds = target_inds[vali_ind]
		Y_vali = Y[vali_ind]
		return XD, XT, train_drug_inds, train_target_inds, Y_train, vali_drug_inds, vali_target_inds, Y_vali, test_drug_inds, test_target_inds, Y_test
	def settingC_split():
		XD, XT, Y, drug_inds, target_inds = load_metz()
		np.random.seed(77)
		#random split to train/test, corresponds to setting C
		trows = list(range(XT.shape[0]))
		np.random.shuffle(trows)
		train_trows = set(trows[:80])
		train_ind = []
		test_ind = []
		for i in range(len(target_inds)):
			if target_inds[i] in train_trows:
				train_ind.append(i)
			else:
				test_ind.append(i)
		train_drug_inds = drug_inds[train_ind]
		train_target_inds = target_inds[train_ind]
		Y_train = Y[train_ind]
		test_drug_inds = drug_inds[test_ind]
		test_target_inds = target_inds[test_ind]
		Y_test = Y[test_ind]
		return XD, XT, train_drug_inds, train_target_inds, Y_train, test_drug_inds, test_target_inds, Y_test
	def settingCvali_split():
		XD, XT, Y, drug_inds, target_inds = load_metz()
		np.random.seed(77)
		#random split to train/test, corresponds to setting C
		trows = list(range(XT.shape[0]))
		np.random.shuffle(trows)
		train_trows = set(trows[:70])
		train_ind = []
		test_ind = []
		kytkD=1
		vali_ind = []
		for i in range(len(target_inds)):
			if target_inds[i] in train_trows:
				train_ind.append(i)
			else:
				# If not in training set, split every other to test set and every other to validation set
				if kytkD==1:
					test_ind.append(i)
					kytkD=2
				else:
					vali_ind.append(i)
					kytkD=1
		train_drug_inds = drug_inds[train_ind]
		train_target_inds = target_inds[train_ind]
		Y_train = Y[train_ind]
		vali_drug_inds = drug_inds[vali_ind]
		vali_target_inds = target_inds[vali_ind]
		Y_vali = Y[vali_ind]
		test_drug_inds = drug_inds[test_ind]
		test_target_inds = target_inds[test_ind]
		Y_test = Y[test_ind]
		return XD, XT, train_drug_inds, train_target_inds, Y_train, vali_drug_inds, vali_target_inds, Y_vali, test_drug_inds, test_target_inds, Y_test
	def settingD_split():
		XD, XT, Y, drug_inds, target_inds = load_metz()
		np.random.seed(77)
		#random split to train/test, corresponds to setting D
		drows = list(range(XD.shape[0]))
		np.random.shuffle(drows)
		train_drows = set(drows[:800])
		trows = list(range(XT.shape[0]))
		np.random.shuffle(trows)
		train_trows = set(trows[:80])
		train_ind = []
		test_ind = []
		for i in range(len(target_inds)):
			if drug_inds[i] in train_drows and target_inds[i] in train_trows:
				train_ind.append(i)
			elif drug_inds[i] not in train_drows and target_inds[i] not in train_trows:
				test_ind.append(i)
		train_drug_inds = drug_inds[train_ind]
		train_target_inds = target_inds[train_ind]
		Y_train = Y[train_ind]
		test_drug_inds = drug_inds[test_ind]
		test_target_inds = target_inds[test_ind]
		Y_test = Y[test_ind]
		return XD, XT, train_drug_inds, train_target_inds, Y_train, test_drug_inds, test_target_inds, Y_test
	def settingDvali_split():
		XD, XT, Y, drug_inds, target_inds = load_metz()
		np.random.seed(77)
		drows = list(range(XD.shape[0]))
		np.random.shuffle(drows)
		# Let's split about in three parts
		train_drows = set(drows[:700])
		trows = list(range(XT.shape[0]))
		np.random.shuffle(trows)
		# Let's split about in three parts
		train_trows = set(trows[:70])
		train_ind = []
		test_ind = []
		kytkD=1
		vali_ind = []
		for i in range(len(target_inds)):
			if drug_inds[i] in train_drows and target_inds[i] in train_trows:
				train_ind.append(i)
			elif drug_inds[i] not in train_drows and target_inds[i] not in train_trows:
				# If not in training set, split every other to test set and every other to validation set
				if kytkD==1:
					test_ind.append(i)
					kytkD=2
				else:
					vali_ind.append(i)
					kytkD=1
		train_drug_inds = drug_inds[train_ind]
		train_target_inds = target_inds[train_ind]
		Y_train = Y[train_ind]
		test_drug_inds = drug_inds[test_ind]
		test_target_inds = target_inds[test_ind]
		Y_test = Y[test_ind]
		vali_drug_inds = drug_inds[vali_ind]
		vali_target_inds = target_inds[vali_ind]
		Y_vali = Y[vali_ind]
		return XD, XT, train_drug_inds, train_target_inds, Y_train, vali_drug_inds, vali_target_inds, Y_vali, test_drug_inds, test_target_inds, Y_test
	start_da = time.clock()
	if setva==1:
		if early==1:
			XD, XT, train_drug_inds, train_target_inds, Y2, vali_drug_inds, vali_target_inds, Y2_vali, test_drug_inds, test_target_inds, Y2_test = settingAvali_split()
		else:
			XD, XT, train_drug_inds, train_target_inds, Y2, test_drug_inds, test_target_inds, Y2_test = settingA_split()
	elif setva==2:
		if early==1:
			XD, XT, train_drug_inds, train_target_inds, Y2, vali_drug_inds, vali_target_inds, Y2_vali, test_drug_inds, test_target_inds, Y2_test = settingBvali_split()
		else:
			XD, XT, train_drug_inds, train_target_inds, Y2, test_drug_inds, test_target_inds, Y2_test = settingB_split()
	elif setva==3:
		if early==1:
			XD, XT, train_drug_inds, train_target_inds, Y2, vali_drug_inds, vali_target_inds, Y2_vali, test_drug_inds, test_target_inds, Y2_test = settingCvali_split()
		else:
			XD, XT, train_drug_inds, train_target_inds, Y2, test_drug_inds, test_target_inds, Y2_test = settingC_split()
	elif setva==4:
		if early==1:
			XD, XT, train_drug_inds, train_target_inds, Y2, vali_drug_inds, vali_target_inds, Y2_vali, test_drug_inds, test_target_inds, Y2_test = settingDvali_split()
		else:
			XD, XT, train_drug_inds, train_target_inds, Y2, test_drug_inds, test_target_inds, Y2_test = settingD_split()
	# IF you are solving a binarized interaction affinity problem with Ki data, the labels will need to be converted according to a threshold value
	Y2th = Y2.copy()
	Y2th[Y2th<7.6] = -1
	Y2th[Y2th>=7.6] = 1
	Y2_testth = Y2_test.copy()
	Y2_testth[Y2_testth<7.6] = -1
	Y2_testth[Y2_testth>=7.6] = 1
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
		tem1 = train_drug_inds.astype(np.int32)
		tem2 = train_target_inds.astype(np.int32)
		krotu = mv_kronecker(K12, K22, tem1, tem2, tem1, tem2)
	n = len(train_drug_inds)
	end_da = time.clock()
	kesto_da = end_da - start_da
	print("The amount of time for generating the data sets and forming the Kronecker product kernel matrix ", kesto_da, " seconds")
	if krorls==1:
		start_rls = time.clock()
		# HYPERPARAMETER FOR KronRLS method! NEED TO BE TUNED SEPARATELY OR SELECTED WISELY! Huge affect on the obtained results.
		rlslam = 1.0 / 32.0
		# DEPENDING WHETHER YOU ARE SOLVING A PROBLEM WITH BINARIZED INTERACTION AFFINITY VALUES OR NOT
		#learner = CGKronRLS(K1 = K12, K2 = K22, Y=Y2, label_row_inds = train_drug_inds, label_col_inds = train_target_inds, regparam = rlslam, maxiter=100)
		learner = CGKronRLS(K1 = K12, K2 = K22, Y=Y2th, label_row_inds = train_drug_inds, label_col_inds = train_target_inds, regparam = rlslam, maxiter=100)
		rlsx = learner.A
		rlsalpha = np.abs(rlsx)
		rlslas = len(rlsalpha[rlsalpha>1e-4])
		print("Number of the non-zero components in the KronRLSn solution ", rlslas)
		rlsP = learner.predict(K12_test, K22_test, test_drug_inds, test_target_inds)
		perf = roc_auc_score(Y2_testth, rlsP)
		print("Prediction accuracy for test data with KronRLS solution ", perf)
		end_rls = time.clock()
		kesto_rls = end_rls - start_rls
		print("The amount of time for finding the optimal dual variable vector for KronRLS ", kesto_rls, " seconds")


# Parameters of the LMBM-KronL0LS
print("Number of training samples %f" %n)
rolista = np.array([0.5, 1.0, 2.0, 4.0, 12.0, 50.0, 500.0])
# Alternative options
#rolista = np.array([0.01, 0.1, 0.5, 1.0, 2.0, 4.0, 7.0, 12.0, 50.0])
#rolista = np.array([0.1, 0.5, 1.0, 2.0, 4.0, 12.0, 50.0])
#rolista = np.array([0.5, 1.0, 2.0, 4.0, 12.0, 50.0])
rolkm = len(rolista)

aucli = list()
harli = list()
muli = list()

# The worst possible solution (is always the same, does not depend on the reviewed data set)
huve = np.array([n, 0.5])

# The reference solution is acquired by setting the L0 norm's coefficient to zero
# Corresponds to the situation that we would like to fit the model to the data as well as possible (nevertheless the sparsity level it leads to)
# In practice the large amount of overfitting happens too (since there is no regulation in this phase)
ro = 0.0
kno = n
# Let's run to the bottom
early = 0
# Initialization of the variables LMBM implementation needs as input
x = np.ones(n)
fu = 1.0
# Acts as a switch (its value changes when LMBM has reached the desired accuracy)
lopa = 3
g = np.ones(n)
mc = 7
# Type conversions
tem5 = test_drug_inds.astype(np.int32)
tem6 = test_target_inds.astype(np.int32)
tem9 = vali_drug_inds.astype(np.int32)
tem10 = vali_target_inds.astype(np.int32)
# Determines whether we wish to solve the reference solution (1) or whether we have determined it earlier
# separately, saved it into a text file and want to read the solution directly from it (0)
if refva==1:
	start_poh = time.clock()
	print(lmbmprogram.lmbm_mod.lmbm(x,fu,func,lopa,vali,g,subgra,mc,n))
	refx = x.copy()
	end_poh = time.clock()
	kesto_poh = end_poh - start_poh
	print("The amount of time it took for LMBM to solve the reference solution (no regulation, no early stopping) ", kesto_poh, " seconds")
	refeta = refx.tolist()
	# Let's save the obtained solution to a text file in case we wish to utilize it later
	np.savetxt('referenssiridge.txt', refeta, delimiter=',')
else:
	refti = open("referenssiridge.txt", "r")
	refte = refti.readlines()
	refx = np.array(refte, dtype='f8')
# The amount of "zeros" in the reference solution
refxabs = np.abs(refx)
reflas = len(refxabs[refxabs>1e-4])
refxh = refx.copy()
for k in range(n):
	if np.abs(refx[k])<1e-4:
		refxh[k] = 0.0
# The prediction accuracy of the reference solution
refhennu = sampled_vec_trick(refxh, K22_test, K12_test, tem6, tem5, tem2, tem1)
refhauc = roc_auc_score(Y2_testth, refhennu)
refve = np.array([reflas, refhauc])
# The criteria value for the reference solution
refero = huve - refve
refpoik = np.abs(refero)
# The best possible amount of nonzero dual variables = 2, the worst = n
refpoikhar = float(refpoik[0]) / float(n-2)
# The choice of the data set has an effect to the maximal prediction accuracy (the main reason is that with checkerboard
# data we have an option to add noise to the data
# In any case, the worst possible prediction accuracy = 0.5 (totally random, no signal to any direction)
if datava==1:
	refpoikauc = float(refpoik[1]) / ((1.0-koh)-0.5)
else:
	# In practice, the best possible prediction accuracy = 1.0
	refpoikauc = float(refpoik[1]) / 0.5
# There is an option to adjust the importance of meeting the objectives with respect to each other
# The larger weight the more it affects to the criteria value and thus the more important it is to get the associated
# objective value (number of nonzero elements / prediction accuracy) as far away as possible from its huve vector value
# Note that the affects of tuning these weights has not been throughoutly examined and can result in unexpected results
harpaino = 1.0
aucpaino = 1.0
refpoikyht = harpaino * refpoikhar + aucpaino * refpoikauc
tiukkak = reflas

print("Reference solution:")
print(refve)
print("The worst possible solution:")
print(huve)


# The switch determining the starting point inheriting procedure
# 0: begin to solve each optimization problem (associated to different rho & k pair) from the vector of ones
# 1: inherit the previous solution point to the starting point of the next optimization problem always inside of the fixed k value
# 2: same as 1, but inherit the solution of the previous optimization problem to the starting point of the next optimization
# problem also when moving to the next k value
laky = 2


start_rok = time.clock()

# Depending of the choice of data set, the suitable step size varies
# Principle 1: if the dimension of the problem is large, do not start with too small step size
# -> that might expand the running times
# Principle 2: if you obtain the optimal solution when running through the first k values, consider decreasing the
# step size in order to get even a better solution
for i in range(2, n+2, 100):
#for i in range(2, n+2, 10):
	if i>n:
		kno = n-1
	else:
		kno = i
	if kno > tiukkak:
		print("Run onto the upper limit of k value!")
		break
	print("The current k value: ", kno)
	aucsis = list()
	harsis = list()
	for j in range(len(rolista)):
		ro = rolista[j]
		print(ro)
		if early==1:
			earlyli = list()
			listax = list()
			listask = list()
			toisto = []
			# The amount of iteration rounds we allow the solution process of the optimization problem attached to
			# the current k&rho value pair to continue without any increase of the prediction accuracy for the validation data
			kipa = 100
		if laky==0:
			x = np.ones(n)
		elif laky==1:
			if j==0:
				x = np.ones(n)
			else:
				x[:-1] = alpharva
				x[-1] = betaloppu
		# If we wish to use the newest solution point as a starting point of the next optimization problem always when such is available
		else:
			# If we are starting to solve the optimization problem related to the first k&rho pair
			if i==2 and j==0:
				x = np.ones(n)
				# One option is also to start to solve the first optimization problem from the reference solution point
				# Did not work that well..
				#x = refx
			else:
				x[:-1] = alpharva
				x[-1] = betaloppu
		# Initialization of the variables LMBM implementation needs as input
		fu = 1.0
		lopa = 3
		g = np.ones(n)
		mc = 7
		print(lmbmprogram.lmbm_mod.lmbm(x,fu,func,lopa,vali,g,subgra,mc,n))
		# Note that the execution of the method can interrupt for other reasons too than cause of the early stopping procedure
		if early==1:
			if len(earlyli)>0:
				xloppu = earlyli[0]
			else:
				xloppu = x.copy()
		else:
			xloppu = x.copy()
		absx = np.abs(xloppu)
		laskuri = len(absx[absx>1e-4])
		xharva = xloppu.copy()
		for k in range(n):
			if np.abs(xloppu[k])<1e-4:
				xharva[k] = 0.0
		harvaennu = sampled_vec_trick(xharva, K22_test, K12_test, tem6, tem5, tem2, tem1)
		aucharva = roc_auc_score(Y2_testth, harvaennu)
		print("The solution for the optimization problem attached to the current k and rho value:")
		print(aucharva)
		print(laskuri)
		# Let's test whether the obtained solution is good enough to be saved into the solution lists or whether we should
		# break and move on to the next optimization problem associated with a new k&rho pair
		if refhauc-aucharva>0.05:
			print("The obtained prediction accuracy is too low.")
			break
		aucsis.append(aucharva)
		harsis.append(laskuri)
		muli.append(ro)
		# Let's test whether the obtained solution gives us a reason to tighten the upper limit of the traversed k values
		temp = np.array([laskuri, aucharva])
		erotus = huve - temp
		poik = np.abs(erotus)
		poikhar = float(poik[0]) / float(n-2)
		# The choice of the data set has an effect to the maximal prediction accuracy (the main reason is that with checkerboard
		# data we have an option to add noise to the data
		# In any case, the worst possible prediction accuracy = 0.5 (totally random, no signal to any direction)
		if datava == 1:
			poikauc = float(poik[1]) / ((1.0-koh)-0.5)
		else:
			poikauc = float(poik[1]) / 0.5
		poikyht = harpaino * poikhar + aucpaino * poikauc
		# The larger poikyht value is the farther away we are from the worst possible solution
		if (poikyht>refpoikyht):
			tiukkak = laskuri
			# Let's update the criteria value to match the best solution found so far
			refpoikyht = poikyht
			print("The new upper limit for k:")
			print(tiukkak)
		# Let's test whether it is still relevant to continue aka if L0-norm still has a nonzero value
		normiLMBM = funcL0(xharva)
		if normiLMBM<1e-1:
			print("No point to continue increasing rho value since the constraint related to the maximum amount of nonzero elements (k kpl) already satisfied!")
			break
	print("Save the results obtained with different rho values to the list.")
	sislkm = len(aucsis)
	for r in range(sislkm):
		aucli.append(aucsis[r])
		harli.append(harsis[r])

end_rok = time.clock()
print("The amount of time it took to scan through the relevant k&rho value pair space: %f seconds" %(end_rok - start_rok))

np.savetxt('aucmuro.txt', aucli, delimiter=',')
np.savetxt('harmuro.txt', harli, delimiter=',')
np.savetxt('romuli.txt', muli, delimiter=',')