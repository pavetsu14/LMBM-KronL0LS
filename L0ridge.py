import numpy as np
from numpy import linalg
from sklearn.metrics import roc_auc_score
from rlscore.kernel import GaussianKernel
from rlscore.utilities.sampled_kronecker_products import sampled_vec_trick
import time
from heapq import nlargest

import lmbmprogram

from data_eplin import generate_chess

# Metz data testailuja varten
from rlscore.learner import CGKronRLS
from rlscore.learner.kron_svm import KronSVM
from rlscore.measure import cindex

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
	# Eksplisiittinen kronecker
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
	# L0-normi
	kla = nlargest(kno, absx)
	regosa = ro * np.sum(absx) - ro * np.sum(kla)
	fuarvo = lossosa + regosa
	return fuarvo


# SIIS MUUTTUJA loss VIITTAA NYT SIIHEN, HALUTAANKO LOSS TERMIIN KERROINTA


def funcL0(x):
	# L0-normi
	absaL0 = np.abs(x)
	klaL0 = nlargest(kno, absaL0)
	L0arvo = ro * np.sum(absaL0) - ro * np.sum(klaL0)
	return L0arvo


def subgra(x,g,n):
	# Eksplisiittinen kronecker
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
	# L0-normi
	gl1 = np.zeros(n)
	for i in range(n):
		if x[i] < 0.0:
			gl1[i] = -1.0
		else:
			gl1[i] = 1.0
	# k-pseudonormeihin liittyva osa
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
		# Y2_vali sen sijaan tulee muuttaa labeleiksi threshold arvon mukaan tai muuten meilla on regressio-ongelma
		#Y2_valith = Y2_vali.copy()
		#Y2_valith[Y2_valith<7.6] = -1
		#Y2_valith[Y2_valith>=7.6] = 1
		#aucoma = roc_auc_score(Y2_valith, ennuvali)
		#print("valitulokset")
		#print(aucoma)
		lopa = 6
	else:
		# Ennustus voi olla jatkuva arvo, silla auc perustuu suuruusjarjestykseen
		#ennuvali = sampled_vec_trick(x, K22_vali, K12_vali, tem10, tem9, tem8, tem7)
		#ennuvali = sampled_vec_trick(x, K22_vali, K12_vali, tem10, tem9, tem4, tem3)
		ennuvali = sampled_vec_trick(x, K22_vali, K12_vali, tem10, tem9, tem2, tem1)
		# Y2_vali sen sijaan tulee muuttaa labeleiksi threshold arvon mukaan tai muuten meilla on regressio-ongelma
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
		#vapoikhar = float(vapoik[0]) / 2498.0
		vapoikauc = float(vapoik[1]) / 0.5
		#vapoikauc = float(vapoik[1]) / 0.3
		skasu = harpaino * vapoikhar + aucpaino * vapoikauc
		xuusi = np.zeros(n)
		xuusi = x.copy()
		listax.append(xuusi)
		listask.append(skasu)
		#print("vali tulokset")
		#print(aucoma)
		#print(harx)
		listal = np.shape(listask)[0]
		# Testiversio 1
		if listal<50:
		#if listal<20:
			lopa = 0
		# Testiversio 1
		elif listal==50:
		#elif listal==20:
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


# VALITSE SOPIVAT PARAMETRIT; LUKITTU L0-NORMIIN

# Kaytetaanko loss osan kertoimena ei mitaan (0) vai 1/2sta (1)
loss = 1
# Skipataanko early stopping proseduuri (0) vai otetaanko se mukaan (1)
early = 1
# Kaytetaanko eksplisiittista kroneckerin kerneltulomatriisia (1) vai implisiittista versiota (2)
kro = 2
# Kaytetaanko shakkilautadataa (1) vai metz dataa (2)
datava = 2



if datava==1:
	# Kuinka paljon kohinaa halutaan sisallyttaa dataan eli kuinka suuri prosenttiosuus labeleista flipataan
	# 0.2 hyva vaihtoehto, pitaisi ainakin sen kanssa viela pystya oppimaan datan rakenne
	koh = 0.2
	# Kuinka vaikea tehtava halutaan luoda eli kuinka isosta shakkilaudasta data generoidaan
	# Erittain helpoissa tehtavissa 2, tapion ja antin artikkelissa ollut 100
	shakki = 5
	start_da = time.clock()
	# Generoidaan training datasetti
	# 100 tilalle voi vaihtaa 1000, jos halutaan testata skaalautuvuutta chess-datalla
	X12, X22, rows2, cols2, Y2 = generate_chess(100, 100, shakki, koh, seed=1)
	kernel12 = GaussianKernel(X12, gamma=1)
	K12 = kernel12.getKM(X12)
	kernel22 = GaussianKernel(X22, gamma=1)
	K22 = kernel22.getKM(X22)
	# Generoidaan test datasetti
	# 500 tilalle voi vaihtaa 5000, jos halutaan testata skaalautuvuutta chess-datalla
	X12_test, X22_test, rows2_test, cols2_test, Y2_test = generate_chess(500, 500, shakki, koh, seed=3)
	K12_test = kernel12.getKM(X12_test)
	K22_test = kernel22.getKM(X22_test)
	# Generoidaan validointidatasetti (tarvittaessa)
	# 500 tilalle voi vaihtaa 5000, jos halutaan testata skaalautuvuutta chess-datalla
	if early==1:
		X12_vali, X22_vali, rows2_vali, cols2_vali, Y2_vali = generate_chess(500, 500, shakki, koh, seed=2)
		K12_vali = kernel12.getKM(X12_vali)
		K22_vali = kernel22.getKM(X22_vali)
	# Kroneckerin kerneltulomatriisi
	if kro==1:
		krotu = naive_kronecker_kernel(K12, K22, rows2, cols2, rows2, cols2)
	else:
		krotu = mv_kronecker(K12, K22, rows2, cols2, rows2, cols2)
	n = len(rows2)
	end_da = time.clock()
	kesto_da = end_da - start_da
	print("Datasettien generointiin ja Kroneckerin kerneltulomatriisin muodostamiseen kului aikaa ", kesto_da, " sekuntia")
elif datava==2:
	# Halutaanko metz data jakaa settingin A (1), B (2), C (3) vai D (4) mukaan
	setva = 4
	# Ratkaistaanko ongelma myos KronRLS avulla (1) vai ei (0)
	krorls = 1
	# Lasketaanko referenssipiste (1) vai luetaanko se suoraan tekstitiedostosta (0)
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
		# Setting Assa ollut
		#test_ind = ind[50000:]
		# Oma lisays, tarkoittaa sita etta ensimmainen mahdollinen sample menee testijoukkoon
		kytkD=1
		vali_ind = []
		for i in range(len(ind)):
			if ind[i] not in train_ind:
				# Jos ei oo training joukossa jaetaan erikseen viela validaation ja testin valilla tyylilla joka toinen
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
		# Lisaysta
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
		#test_drug_ind = set(drug_ind[800:])
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
		# Setting Bssa ollut train_drows = set(drows[:800])
		train_drows = set(drows[:700])
		train_ind = []
		test_ind = []
		# Oma lisays, tarkoittaa sita etta ensimmainen mahdollinen sample menee testijoukkoon
		kytkD=1
		vali_ind = []
		for i in range(len(drug_inds)):
			if drug_inds[i] in train_drows:
				train_ind.append(i)
			else:
				#test_ind.append(i)
				# Jos ei oo training joukossa jaetaan erikseen viela validaation ja testin valilla tyylilla joka toinen
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
		# Lisaysta
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
		# Setting Cssa ollut train_trows = set(trows[:80])
		train_trows = set(trows[:70])
		train_ind = []
		test_ind = []
		# Oma lisays, tarkoittaa sita etta ensimmainen mahdollinen sample menee testijoukkoon
		kytkD=1
		vali_ind = []
		for i in range(len(target_inds)):
			if target_inds[i] in train_trows:
				train_ind.append(i)
			else:
				#test_ind.append(i)
				# Jos ei oo training joukossa jaetaan erikseen viela validaation ja testin valilla tyylilla joka toinen
				if kytkD==1:
					test_ind.append(i)
					kytkD=2
				else:
					vali_ind.append(i)
					kytkD=1
		train_drug_inds = drug_inds[train_ind]
		train_target_inds = target_inds[train_ind]
		Y_train = Y[train_ind]
		# Lisaysta
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
		#random split to train/test, corresponds to setting D
		drows = list(range(XD.shape[0]))
		np.random.shuffle(drows)
		# Jaetaankin about kolmeen osaan
		#train_drows = set(drows[:800])
		#train_drows = set(drows[:540])
		train_drows = set(drows[:700])
		trows = list(range(XT.shape[0]))
		np.random.shuffle(trows)
		# Jaetaankin about kolmeen osaan
		#train_trows = set(trows[:80])
		#train_trows = set(trows[:53])
		train_trows = set(trows[:70])
		train_ind = []
		test_ind = []
		# Oma lisays, tarkoittaa sita etta ensimmainen mahdollinen sample menee testijoukkoon
		kytkD=1
		vali_ind = []
		for i in range(len(target_inds)):
			if drug_inds[i] in train_drows and target_inds[i] in train_trows:
				train_ind.append(i)
			elif drug_inds[i] not in train_drows and target_inds[i] not in train_trows:
				#test_ind.append(i)
				# Jos ei oo training joukossa jaetaan erikseen viela validaation ja testin valilla tyylilla joka toinen
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
		# Lisaysta
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
	# Y2 ja Y2_test tulee muuttaa labeleiksi threshold arvon mukaan tai muuten meilla on regressio-ongelma
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
	# Kroneckerin kerneltulomatriisi
	if kro==1:
		krotu = naive_kronecker_kernel(K12, K22, train_drug_inds, train_target_inds, train_drug_inds, train_target_inds)
	else:
		tem1 = train_drug_inds.astype(np.int32)
		tem2 = train_target_inds.astype(np.int32)
		krotu = mv_kronecker(K12, K22, tem1, tem2, tem1, tem2)
	n = len(train_drug_inds)
	end_da = time.clock()
	kesto_da = end_da - start_da
	print("Datasettien generointiin ja Kroneckerin kerneltulomatriisin muodostamiseen kului aikaa ", kesto_da, " sekuntia")
	if krorls==1:
		start_rls = time.clock()
		# rlslam voi olla mika tahansa MIKA ON PARAS?
		rlslam = 1.0 / 32.0
		# TAMA ON regularized least-squares regression with paired-input (dyadic) data and Kronecker kernels
		# TESTAA MITEN TOIMII, JOS KAYTETAAN EI THRESHOLDATTUA ARVOA LABELINA?
		#learner = CGKronRLS(K1 = K12, K2 = K22, Y=Y2, label_row_inds = train_drug_inds, label_col_inds = train_target_inds, regparam = rlslam, maxiter=100)
		learner = CGKronRLS(K1 = K12, K2 = K22, Y=Y2th, label_row_inds = train_drug_inds, label_col_inds = train_target_inds, regparam = rlslam, maxiter=100)
		rlsx = learner.A
		rlsalpha = np.abs(rlsx)
		rlslas = len(rlsalpha[rlsalpha>1e-4])
		print("RLScoren KronRLSn nollasta eroavien komponenttien lkm ", rlslas)
		rlsP = learner.predict(K12_test, K22_test, test_drug_inds, test_target_inds)
		#perfc = cindex(Y2_testth, rlsP)
		perf = roc_auc_score(Y2_testth, rlsP)
		print("RLScoren KronRLSn ennustustarkkuus ", perf)
		end_rls = time.clock()
		kesto_rls = end_rls - start_rls
		print("Optimaalisen duaalimuuttujapisteen loytamiseen RLScoren KronRLSlla kului aikaa ", kesto_rls, " sekuntia")


# LMBM-SVM:n muuttujien ja parametrien alustukset
print("Training samplien lkm on %f" %n)
#rolista = np.array([0.01, 0.1, 0.5, 1.0, 2.0, 4.0, 7.0, 12.0, 50.0])
#rolista = np.array([0.1, 0.5, 1.0, 2.0, 4.0, 12.0, 50.0])
rolista = np.array([0.5, 1.0, 2.0, 4.0, 12.0, 50.0, 500.0])
#rolista = np.array([0.5, 1.0, 2.0, 4.0, 12.0, 50.0])
rolkm = len(rolista)

aucli = list()
harli = list()
muli = list()

# Datasetista riippumatta
huve = np.array([n, 0.5])
# Referenssipisteen ratkaisu on saatu asettamalla L0-termin kerroin 0.0ksi
# Kuvastaa tilannetta, etta haluttaisiin vain sovittaa malli mahdollisimman hyvin dataan valittamatta harvuusasteesta
# (kaytannossa tapahtuu myos overfittingia, jota ei nyt esteta mitenkaan)
# (luultavasti yo. syysta ei saada talla tavalla parasta mahdollista ennustustarkkuutta)
# HUOM. auc arvo on annettu silla tarkkuudella, jolla tietokone sen komentoriville tulostaa
#refve = np.array([2500, 0.77822021597312352])
ro = 0.0
kno = n
# Ajetaan pohjaan taman maarittamiseksi
early = 0
x = np.ones(n)
fu = 1.0
lopa = 3
g = np.ones(n)
mc = 7
tem5 = test_drug_inds.astype(np.int32)
tem6 = test_target_inds.astype(np.int32)
tem9 = vali_drug_inds.astype(np.int32)
tem10 = vali_target_inds.astype(np.int32)
if refva==1:
	start_poh = time.clock()
	# LMBM-SVM
	print(lmbmprogram.lmbm_mod.lmbm(x,fu,func,lopa,vali,g,subgra,mc,n))
	refx = x.copy()
	end_poh = time.clock()
	kesto_poh = end_poh - start_poh
	print("LMBM-SVMlla kului referenssipisteen (ei regularisaatiota, ei early stoppingia) ratkaisemiseen ", kesto_poh, " sekuntia")
	refeta = refx.tolist()
	np.savetxt('referenssiridge.txt', refeta, delimiter=',')
else:
	refti = open("referenssiridge.txt", "r")
	refte = refti.readlines()
	refx = np.array(refte, dtype='f8')
refxabs = np.abs(refx)
reflas = len(refxabs[refxabs>1e-4])
refxh = refx.copy()
for k in range(n):
	if np.abs(refx[k])<1e-4:
		refxh[k] = 0.0
refhennu = sampled_vec_trick(refxh, K22_test, K12_test, tem6, tem5, tem2, tem1)
refhauc = roc_auc_score(Y2_testth, refhennu)
refve = np.array([reflas, refhauc])
# Referenssipistetta vastaava painotettu arvo
refero = huve - refve
refpoik = np.abs(refero)
# Paras mahdollinen harvuus on 2, mutta toisaalta huonoin mahdollinen 2500
refpoikhar = float(refpoik[0]) / float(n-2)
# Maksimaalinen ennustustarkkuus on 0.8, mutta toisaalta 0.5 vastaa taysin randomia (huonoin mahdollinen)
if datava==1:
	refpoikauc = float(refpoik[1]) / ((1.0-koh)-0.5)
else:
	# Kuitenkin muilla datoilla pitaisi periaatteessa pyrkia 1.0 tarkkuuteen
	refpoikauc = float(refpoik[1]) / 0.5
# Kuvastaa tavoitteiden toteutumisen tarkeytta suhteessa toisiinsa
# Mita suurempi paino, sita enemman ko. tavoitteen arvo vaikuttaa painotettuun arvoon ja sita tarkeampaa
# sen saaminen mahd kauaksi huve vektorin ko. tavoitteeseen liittyvasta arvosta on
harpaino = 1.0
aucpaino = 1.0
refpoikyht = harpaino * refpoikhar + aucpaino * refpoikauc
tiukkak = reflas

print("Referenssipiste:")
print(refve)
print("Huonoin mahdollinen tavoitepiste:")
print(huve)


# Kytkin, jonka avulla valitaan
# 0: kaytetaanko jokaisen tehtavan lahtopisteena ykkosvektoria
# 1: hyodynnetaanko fixatun k sisalla edellisen tehtavan ratkaisupistetta yhta isompaan ro arvoon liittyvan
# tehtavan lahtopisteena
# 2: hyodynnetaanko nykyisella k:lla ja suurimmalla lapikaydylla ro arvolla saatua
# ratkaisupistetta seuraavaan (yhden askeleen verran kasvatettuun) k arvoon ja ensimmaiseen (pienimpaan ro arvoon)
# liittyvan optimointiongelman lahtopisteena
laky = 2


start_rok = time.clock()


for i in range(2, n+2, 100):
#for i in range(2, n+2, 10):
	if i>n:
		kno = n-1
	else:
		kno = i
	if kno > tiukkak:
		print("Tormattiin ylarajaan!")
		break
	print("Taman hetkinen kno arvo: ", kno)
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
			kipa = 100
		# Jos halutaan kayttaa samaa lahtopistetta jokaiseen khon liittyvan ensimmaisen tehtavan ratkaisuun
		if laky==0:
			x = np.ones(n)
		elif laky==1:
			if j==0:
				x = np.ones(n)
			else:
				x[:-1] = alpharva
				x[-1] = betaloppu
		# Jos halutaan hyodyntaa edellisella k:lla saatua uusinta optimipistetta seuraavaan k:hon liittyvan
		# ensimmaisen tehtavan lahtopisteena
		else:
			# Mikali ollaan ratkaisemassa ensimmaiseen k:hon liittyvaa ongelmaa ensimmaisella ro:lla
			if i==2 and j==0:
			# VAIHDETAAN TAHAN TILALLE VERSIO, JOSSA JOKAISEEN k JA PIENIMPAAN ro LIITTYVAN TEHTAVAN ALOITUSPISTEENA
			# KAYTETAAN REFERENSSIPISTETTA (sama kuin pienin kaytossa oleva ro olisi nolla)
			#if j==0:
				x = np.ones(n)
				# MUOKATAAN TATA SITEN, ETTA ALOITTAA REFERENSSIPISTEEN RATKAISUPISTEESTA
				#x = refx
			else:
				x[:-1] = alpharva
				x[-1] = betaloppu
		fu = 1.0
		lopa = 3
		g = np.ones(n)
		mc = 7
		# LMBM-SVM
		print(lmbmprogram.lmbm_mod.lmbm(x,fu,func,lopa,vali,g,subgra,mc,n))
		# Suoritus saattaa keskeytya myos muihin kuin early stopping kriteeriin
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
		print("Lopputulos auc LMBM-SVM jalkeen:")
		print(aucharva)
		print(laskuri)
		# Testataan onko saavutettu ratkaisu riittavan hyva tullakseen tallennetuksi listaan
		# Shakkilauta-datan referenssi auc arvo joka on saatu ilman regularisaatiota
		#if aucharva<=0.77822021597312352 and (0.77822021597312352-aucharva>0.05):
		# Laskettu erikseen puolen tunnin ajolla ilman early stoppingia (kiersi 10 000 iteraatiokierrosta)
		#if aucharva<=refhauc and (refhauc-aucharva>0.01):
		# kait sama kuin
		if refhauc-aucharva>0.05:
			print("Ennustustarkkuus tippunut liian alas")
			break
		# Jos ei breikata yo. ehdosta, tallennetaan listaan saadut tulokset
		aucsis.append(aucharva)
		harsis.append(laskuri)
		muli.append(ro)
		# Testataan voidaanko kiristaa k:n ylarajaa nykyisen tuloksen varjolla
		temp = np.array([laskuri, aucharva])
		erotus = huve - temp
		poik = np.abs(erotus)
		poikhar = float(poik[0]) / float(n-2)
		# Paras mahdollinen harvuus on 2, mutta toisaalta huonoin mahdollinen 2500
		#poikhar = float(poik[0]) / 2498.0
		poikauc = float(poik[1]) / 0.5
		# Maksimaalinen ennustustarkkuus on 0.8, mutta toisaalta 0.5 vastaa taysin randomia (huonoin mahdollinen)
		#poikauc = float(poik[1]) / 0.3
		poikyht = harpaino * poikhar + aucpaino * poikauc
		# Pitaa olla suurempi kuin painotettu referenssipisteen arvo silla koitetaan paasta mahd kauaksi huve vektorista
		if (poikyht>refpoikyht):
		# Tai sitten vaihtoehtoisesti, jos ollaan hyvin lahella, kannattaa tutkia hyodyttaisiinko nykyisesta ratkaisusta
		#if (poikyht>refpoikyht) or (np.abs(refpoikyht-poikyht)<1e-4):
			# Testataan kannattaako loydettya ylarajaa muuttaa
			#if laskuri<tiukkak:
				#tiukkak = laskuri
				# Paivitetaan se arvo, johon verrataan, koska muuten paadytaan etsimaan vain optimaalista
				# harvuutta valittamatta siita miten se vaikuttaa AUCiin
				#refpoikyht = poikyht
				#print("k:n ylaraja paivityksen jalkeen")
				#print(tiukkak)
			#else:
				#print("Ylarajan arvo olisi suurentunut")
			# Taytyisiko sallia ylarajan muutos molempiin suuntiin, varsinkin nyt kun kaytossa nain tiukka auc raja?
			tiukkak = laskuri
			# Paivitetaan se arvo, johon verrataan, koska muuten paadytaan etsimaan vain optimaalista
			# harvuutta valittamatta siita miten se vaikuttaa AUCiin
			refpoikyht = poikyht
			print("k:n ylaraja paivityksen jalkeen")
			print(tiukkak)
		# Testataan onko yha relevanttia jatkaa eli onko L0-normilla yha nollasta eroava arvo
		normiLMBM = funcL0(xharva)
		if normiLMBM<1e-1:
			print("Ro breikkiin paasty")
			break
	print("Kay tallentamassa eri ro:illa saadut tulokset listaan")
	sislkm = len(aucsis)
	for r in range(sislkm):
		aucli.append(aucsis[r])
		harli.append(harsis[r])

end_rok = time.clock()
print("Kaikkien yhdistelmien testaamiseen kului aikaa %f sekuntia" %(end_rok - start_rok))

np.savetxt('aucmuro.txt', aucli, delimiter=',')
np.savetxt('harmuro.txt', harli, delimiter=',')
np.savetxt('romuli.txt', muli, delimiter=',')


# HUOM. rows2 = train_drug_inds, cols2 = train_target_inds jne.