# LMBM-KronL0LS

This repository contains the files needed to run LMBM-KronL0LS method with python. The method has been developed for solving pairwise interaction prediction problems with L0-penalized least squares under four different experimental settings. The method can handle both binary and continuous interaction affinity labels. In addition to that, method has been shown to perform well also with noisy and highly nonlinear pairwise data. The method utilizes nonsmooth optimization method called Limited Memory Bundle Method, Haarala et al (2007), for solving the developed bi-objective optimization problem and generalized vec trick for forming implicit Kronecker product kernel matrices.

Pairwise interaction prediction problems are closely related to drug-target problems and the results might be useful e.g. in the field of drug discovery.

The method has been introduced in the article: TO BE PUBLISHED

## Data

1. The codes to generate Checkerboard data (in the data_eplin.py file) have been developed by Antti Airola.
2. Ki data is related to the experiments done with the affinities found by Metz et al (2011).

Ki data has been applied before for the real-valued drug-target binding affinity prediction experiments done in the article by Pahikkala et al (2015).

Ki data is available both in http://staff.cs.utu.fi/~aatapa/data/DrugTarget/ and in this repository through three separate text files containing the drug-target interaction affinities, the drug-drug structural fingerprint similarities computed with the 2D Tanimoto coefficients and the target-target sequence similarities computed with the normalized Smith-Waterman (SW) score.

Special thanks to Antti Airola, Napsu Karmitsa and Tapio Pahikkala for co-operation and the development of the methods and techniques utilized in this work. Great inspiration was derived from the article by Airola et al (2018).

Note. Use f2py to generate lmbmprogram.so package from the lmbmprogram.f95 and lmbmprogram.pyf files. You can find instructions for doing this e.g. from
https://jiffyclub.github.io/numpy/user/c-info.python-as-glue.html
or
https://f2py-users.cens.ioc.narkive.com/c3a0uNn9/wrapping-large-fortran-code-with-f2py


References:

A. Airola, T. Pahikkala. Fast Kronecker product kernel methods via generalized vec trick. IEEE Transactions on Neural Networks and Learning Systems, 29(8):3374- 3387, 2018.

N. Haarala, K. Miettinen, M. M. Mäkelä. Globally Convergent Limited Memory Bundle Method for Large-Scale Nonsmooth Optimization. Mathematical Programming, Vol. 109(1):181-205, 2007.

J. T. Metz, E. F. Johnson, N. B. Soni, P. J. Merta, L. Kifle, P. J. Hajduk. Navigating the kinome. Nature Chemical Biology, 7:200–202, 2011.

T. Pahikkala, A. Airola, S. Pietilä, S. Shakyawar, A. Szwajda, J. Tang, T. Aittokallio. Toward more realistic drug-target interaction predictions. Briefings in Bioinformatics, 16(2):325-337, 2015.
