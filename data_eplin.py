import numpy as np
import random


def generate_chess(size1, size2, chesize, fli, frac=0.25, seed=10):
    np.random.seed(seed)
    random.seed(seed)
    X1 = chesize*np.random.random(size1)
    X2 = chesize*np.random.random(size2)
    paircount = int(frac*size1*size2)
    rowind = np.random.randint(0, size1, paircount)
    colind = np.random.randint(0, size2, paircount)
    #foo = rowind + colind
    #I = np.argsort(foo)
    I = np.argsort(colind)
    rowind = rowind[I]
    colind = colind[I]
    Y = []
    for rind, cind in zip(rowind, colind):
        y1 = int(X1[rind])%2
        y2 = int(X2[cind])%2
        y = np.logical_xor(y1, y2)
        if fli==0.0:
            Y.append(y)
        else:
            if random.random() > fli:   
                Y.append(y)
            else:
                Y.append(y == False)
    Y = 2.*np.array(Y)-1
    X1 = X1.reshape(size1, 1)
    X2 = X2.reshape(size2, 1)
    rowind = np.array(rowind, dtype=np.int32)
    colind = np.array(colind, dtype=np.int32)
    return X1, X2, rowind, colind, Y