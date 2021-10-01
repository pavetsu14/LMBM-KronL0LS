# LMBM-KronL0LS

This repository contains the files needed to run LMBM-KronL0LS method with python. The method has been developed for solving pairwise interaction prediction problems with L0-penalized least squares under four different experimental settings. The method can handle both binary and, often considered as a more realistic approach, continuous interaction affinity labels. The method utilizes nonsmooth optimization method called Limited Memory Bundle Method, Haarala et al (2007), for solving the developed bi-objective optimization problem and generalized vec trick, Airola et al (2018), for forming implicit Kronecker product kernel matrices. As we are using the LMBM to minimize the objective function, we are able to take an exceptional approach for predicting pairwise interactions. In other words, we will not aim only for the best possible prediction performance but simultaneously pursue the sparsest solution possible capable of maintaining satisfactory prediction accuracy. Special thanks to Antti Airola, Napsu Karmitsa and Tapio Pahikkala for co-operation and the development of the methods and techniques utilized in this work.

Pairwise interaction prediction problems are closely related to drug-target prediction problems and thus the developed method might be useful in the field of drug discovery.

The method has been introduced in the article: TO BE PUBLISHED

**Abstract** Large time, money and resource costs have motivated researchers to constantly develop innovative technologies for the exploitation of new drugs. However, current knowledge about the drug–target space is limited. In this paper, a new nonsmooth optimization based method, LMBM-KronL0LS, for solving large-scale pairwise interaction affinity prediction  problems  in  the  framework  of  training  Kronecker  product  kernel  methods,  is  introduced. The proposed method allows balancing the prediction accuracy with sparsity of the learned model. We apply a least squares approach for loss function and a continuous formulation of L0-pseudonorm for regularization purposes. In addition to that, a novel procedure for tuning the hyperparameters of the developed method is introduced. In our experiments, we will explore the performance of the method under four distinct experimental settings, including zero-shot learning, and with both binary and continuous interaction affinity labels. The obtained results demonstrate the advantages of the proposed method; especially with continuous interaction affinity values. Another major remark is that LMBM-KronL0LS performs best in the most complex situations (with noisy data and highly nonlinear patterns).

## F2PY interface

Use f2py to generate lmbmprogram.so package from the lmbmprogram.f95 and lmbmprogram.pyf files. You can find instructions for doing this e.g. from

https://numpy.org/devdocs/f2py/usage.html

https://www.numfys.net/howto/F2PY/

https://jiffyclub.github.io/numpy/user/c-info.python-as-glue.html

## Data

1. The codes to generate Checkerboard data (in the data_eplin.py file) have been developed by Antti Airola.
2. Ki data is related to the experiments done with the affinities found by Metz et al (2011).

Ki data has been applied before for the real-valued drug-target binding affinity prediction experiments done in the article by Pahikkala et al (2015).

Ki data is available both in http://staff.cs.utu.fi/~aatapa/data/DrugTarget/ and in this repository through three separate text files containing the drug-target interaction affinities, the drug-drug structural fingerprint similarities computed with the 2D Tanimoto coefficients and the target-target sequence similarities computed with the normalized Smith-Waterman (SW) score.

## References:

A. Airola, T. Pahikkala. Fast Kronecker product kernel methods via generalized vec trick. IEEE Transactions on Neural Networks and Learning Systems, 29(8):3374- 3387, 2018.

N. Haarala, K. Miettinen, M. M. Mäkelä. Globally Convergent Limited Memory Bundle Method for Large-Scale Nonsmooth Optimization. Mathematical Programming, Vol. 109(1):181-205, 2007.

J. T. Metz, E. F. Johnson, N. B. Soni, P. J. Merta, L. Kifle, P. J. Hajduk. Navigating the kinome. Nature Chemical Biology, 7:200–202, 2011.

T. Pahikkala, A. Airola, S. Pietilä, S. Shakyawar, A. Szwajda, J. Tang, T. Aittokallio. Toward more realistic drug-target interaction predictions. Briefings in Bioinformatics, 16(2):325-337, 2015.
