import unittest
import numpy as np
import math

def evalLagrangeBasis1D(variate,degree,basis_idx,domain):
    step = 2/degree
    xj = np.arange(-1,2,step)
    val = 1 
    for j in range(0,degree+1):
        if j == basis_idx:
            val = val
        else:
            val = val*(variate - xj[j])/(xj[basis_idx] - xj[j])
    return val


# def evalBernsteinBasis1D(variate, degree, basis_idx):
#     if variate <= 0:
#         variate = abs(-1 - variate) / (2)
#         val = math.comb(degree,basis_idx)*(variate**basis_idx)*(1-variate)**(degree-basis_idx)
#     else:
#         variate = abs(-1 - variate) / (2)
#         val = math.comb(degree,basis_idx)*(variate**basis_idx)*(1-variate)**(degree-basis_idx)
#     return val

# def evalBernsteinBasis1D( variate, degree, basis_idx, domain):
#     domain = [-1,1] #if using Gauss-Legendre Quadrature points
#     p = degree
#     i = basis_idx
#     # v = (variate + 1)/2
#     v = (1 / (domain[-1] - domain[0])) * (variate - domain[0])
#     t1 = math.comb(p,i)
#     t2 = v**i
#     t3 = (1-v)**(p-i)
#     val = t1*t2*t3
#     return val

def evalBernsteinBasis1D(variate,degree,basis_idx,domain):
    v = (1/(domain[-1]-domain[0]))*(variate) + 0.5 - ((domain[-1] - domain[0])/(2) + domain[0])
    # print(v)
    # v = (variate + 1)/2 # This work when we are using qp from a basis with a [-1,1] domain
    # v = (1/(domain[-1] - domain[0]))*(variate - domain[0])
    term1 = math.comb(degree,basis_idx)
    term2 = v**basis_idx
    term3 = (1 - v)**(degree-basis_idx)
    val = term1 * term2 * term3 
    return val


def evalLegendreBasis1D(variate,degree,basis_idx,domain):
    if basis_idx == 0:
        val = 1
    elif basis_idx == 1:
        val = variate
    else:
        i = basis_idx - 1
        val = ((i + 1)**(-1)) * ((2*i+1)*variate*evalLegendreBasis1D(variate,degree,i,domain) - i*evalLegendreBasis1D(variate,degree,i-1,domain))
    return val