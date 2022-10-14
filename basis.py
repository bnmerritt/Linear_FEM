# evalute Lagrange Basis

import unittest
import numpy as np
import math

def evalLagrangeBasis1D(variate,degree,basis_idx):
    step = 2/degree
    xj = np.arange(-1,2,step)
    val = 1 
    for j in range(0,degree+1):
        if j == basis_idx:
            val = val
        else:
            val = val*(variate - xj[j])/(xj[basis_idx] - xj[j])
    
    return val


def evaluateBernsteinBasis1D(variate, degree, basis_idx):
    if variate <= 0:
        variate = abs(-1 - variate) / (2)
        val = math.comb(degree,basis_idx)*(variate**basis_idx)*(1-variate)**(degree-basis_idx)
    else:
        variate = abs(-1 - variate) / (2)
        val = math.comb(degree,basis_idx)*(variate**basis_idx)*(1-variate)**(degree-basis_idx)
    return val
