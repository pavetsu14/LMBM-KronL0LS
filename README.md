# LMBM-KronL0LS

This repository contains the files needed to run LMBM-KronL0LS method with python. The method can be executed by running either LMBMKronL0LS.py or LMBMKronL0LS_checker.py file. The first one mentioned contains everything needed to apply the method for three different real-world data sets and the second one contains everything needed to apply the method for generate-it-yourself checkerboard data.

The method has been developed for solving pairwise interaction prediction problems with L0-penalized least squares
under four different experimental settings. The method can handle both binary and, often considered as a more realistic approach,
continuous interaction affinity labels. The method utilizes nonsmooth optimization method called Limited Memory Bundle Method,
Haarala et al (2007), for solving the developed bi-objective optimization problem and generalized vec trick, Airola et al (2018),
for forming implicit Kronecker product kernel matrices. As we are using the LMBM to minimize the objective function,
we are able to take an exceptional approach for predicting pairwise interactions. In other words, we will not aim only for
the best possible prediction performance but simultaneously pursue the sparsest solution possible capable of maintaining
satisfactory prediction accuracy. Special thanks to Antti Airola, Napsu Karmitsa and Tapio Pahikkala for co-operation and
the development of the methods and techniques utilized in this work and to Riikka Numminen for your massive efforts in the revision phase.

Pairwise interaction prediction problems are closely related to drug-target prediction problems and thus the developed method
might be useful in the field of drug discovery.

The method has been introduced in the article: REVISION SUBMITTED TO IEEE ACCESS 11.11.2022

**Abstract** Large time, money and resource costs have motivated researchers to constantly develop innovative technologies
for the exploitation of new drugs. However, current knowledge about the drug–target space is limited and the existing
data sets are sparse, i.e. interaction labels have been validated for only
some drug-target pairs. In this paper, a new nonsmooth optimization-based method LMBM-KronL0LS is
introduced for solving large-scale pairwise interaction affinity prediction problems in the framework of
training Kronecker product kernel methods. The proposed method allows balancing between the prediction
accuracy of the learned model and the sparsity of the obtained solution. Here, the level of sparsity refers
to the fraction of zeroes in the solution vector. We apply a least squares approach for loss function and a
continuous formulation of L0-pseudonorm for regularization purposes. In the numerical experiments, we use
three benchmark data sets and two simulated data sets to demonstrate the performance of LMBM-Kronℓ0LS
under four distinct experimental settings, including zero-shot learning, and with both binary and continuous
interaction affinity labels. The results show that the method is capable of finding sparse solutions without
sacrificing too much in the prediction performance.

## F2PY interface

Use f2py to generate lmbmprogram.so package from the lmbmprogram.f95 and lmbmprogram.pyf files. You can find
instructions for doing this e.g. from

https://numpy.org/devdocs/f2py/usage.html

https://www.numfys.net/howto/F2PY/

https://jiffyclub.github.io/numpy/user/c-info.python-as-glue.html

## Data

This repository does not contain all the data used for testing the method. However, we included one of the utilized real-world data sets in here so you can test the method without the need to worry about the data availability. Ki data is related to the experiments done with the affinities found by Metz et al (2011). Ki data has been applied before for the real-valued drug-target binding affinity prediction experiments done in the article by Pahikkala et al (2015).

Ki data is available both in http://staff.cs.utu.fi/~aatapa/data/DrugTarget/ and in this repository through three separate
text files containing the drug-target interaction affinities, the drug-drug structural fingerprint similarities computed
with the 2D Tanimoto coefficients and the target-target sequence similarities computed with the normalized Smith-Waterman (SW) score.

In practice, if you utilize functions in load_split_data.py, please use load_metz function.

In addition to the Ki data, the codes for generating checkerboard data (in generate_nonlinear_data.py file, developed by Antti Airola) have been added to the repository in order to offer different type of data, and most importantly, some data with which running the method is **fast**.

## References:

A. Airola, T. Pahikkala. Fast Kronecker product kernel methods via generalized vec trick. IEEE Transactions on Neural
Networks and Learning Systems, 29(8):3374- 3387, 2018.

N. Haarala, K. Miettinen, M. M. Mäkelä. Globally Convergent Limited Memory Bundle Method for Large-Scale Nonsmooth
Optimization. Mathematical Programming, Vol. 109(1):181-205, 2007.

J. T. Metz, E. F. Johnson, N. B. Soni, P. J. Merta, L. Kifle, P. J. Hajduk. Navigating the kinome. Nature Chemical
Biology, 7:200–202, 2011.

T. Pahikkala, A. Airola, S. Pietilä, S. Shakyawar, A. Szwajda, J. Tang, T. Aittokallio. Toward more realistic drug-target
interaction predictions. Briefings in Bioinformatics, 16(2):325-337, 2015.
