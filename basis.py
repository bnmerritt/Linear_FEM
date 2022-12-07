import unittest
import numpy as np
import math
import sympy

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
    # v = (1/(domain[-1]-domain[0]))*(variate) + 0.5 - ((domain[-1] - domain[0])/(2) + domain[0])
    if basis_idx < 0 or basis_idx > degree:
        val = 0
    else:
        v = affine_mapping_1D(domain, [0, 1], variate)
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

def evalSplineBasis1D(variate,elem_extraction_operator,basis_idx,domain):
    degree = elem_extraction_operator.shape[0] - 1
    vector = np.zeros(degree+1)
    for bidx in range(0,degree+1):
        vector[bidx] = evalBernsteinBasis1D(variate,degree,bidx,domain)
    val = np.matmul(elem_extraction_operator[basis_idx,:],vector)
    return val

def affine_mapping_1D( domain, target_domain, x ):
    A = np.array( [ [ 1.0, domain[0] ], [ 1.0, domain[1] ] ] )
    b = np.array( [target_domain[0], target_domain[1] ] )
    c = np.linalg.solve( A, b )
    fx = c[0] + c[1] * x
    return fx

def eigenvaluesLegendreBasis( degree ):
    x = sympy.symbols( 'x', real = True )
    poly_fun = sympy.poly( symLegendreBasis( degree, degree, [-1, 1], x ) )
    comp_matrix = computeCompanionMatrix( poly_fun )
    eig_vals = np.sort( np.linalg.eigvals( comp_matrix ) )
    eig_vals = [ float( np.real( val ) ) for val in eig_vals ]
    return eig_vals

def computeCompanionMatrix( poly_fun ):
    coeffs = poly_fun.all_coeffs()
    coeffs.reverse()
    coeffs = [ float( val / coeffs[-1] ) for val in coeffs ]
    coeffs = np.array( coeffs[0:-1] )
    comp_matrix = np.zeros( shape = ( len( coeffs ) , len( coeffs ) ) )
    comp_matrix[:,-1] = -1 * coeffs
    comp_matrix[1:, 0:-1] = np.eye( ( len( coeffs ) - 1 ) )
    return comp_matrix

def symLegendreBasis( degree, basis_idx, domain, variate ):
    x = sympy.symbols( 'x', real = True )
    if degree == 0:
        p = sympy.Poly( 1, x )
    else:
        term_1 = 1.0 / ( ( 2.0 ** degree ) * sympy.factorial( degree ) )
        term_2 = ( ( x**2) - 1.0 ) ** degree 
        term_3 = sympy.diff( term_2, x, degree )
        p = term_1 * term_3
        p = sympy.poly( sympy.simplify( p ) )
    return p

def evalBernsteinBasisDeriv(degree, basis_idx, deriv, domain, variate):
    domain_change = 1 / (domain[-1] - domain[0])
    if deriv == 0:
        val = evalBernsteinBasis1D(variate,degree,basis_idx,domain)
    elif deriv == 1:
        term1 = evalBernsteinBasis1D(variate,degree-1,basis_idx-1,domain)
        term2 = evalBernsteinBasis1D(variate,degree-1,basis_idx,domain)
        val = degree * (term1-term2) * domain_change
    elif deriv == 2:
        term1 = evalBernsteinBasis1D(variate,degree-2,basis_idx-2,domain)
        term2 = 2*evalBernsteinBasis1D(variate,degree-2,basis_idx-1,domain)
        term3 = evalBernsteinBasis1D(variate,degree-2,basis_idx,domain)
        val = degree* (degree-1) * (term1 - term2 + term3) * domain_change**2
    return val

def evalSplineBasisDeriv1D(extraction_operator, basis_idx, deriv, domain, variate):
    # domain_change = 1 / (domain[-1] - domain[0])
    degree = extraction_operator.shape[0] - 1 
    vector = np.zeros(degree+1)
    for bidx in range(0,degree+1):
        vector[bidx] = evalBernsteinBasisDeriv(degree,bidx,deriv,domain,variate)
    val = np.matmul(extraction_operator[basis_idx,:],vector)
    return val

# class test_evalSplineBasis1D( unittest.TestCase ):
#     def test_linear_basis_unit_domain( self ):
#         elem_degree = 1
#         elem_domain = [0, 1]
#         C = np.eye( elem_degree+ 1 )
#         self.assertAlmostEqual( evalSplineBasis1D( variate = 0, elem_extraction_operator = C, basis_idx = 0, domain = elem_domain ), 1.0 )
#         self.assertAlmostEqual( evalSplineBasis1D( variate = 0, elem_extraction_operator = C, basis_idx = 1, domain = elem_domain ), 0.0 )
#         self.assertAlmostEqual( evalSplineBasis1D( variate = 1, elem_extraction_operator = C, basis_idx = 0, domain = elem_domain ), 0.0 )
#         self.assertAlmostEqual( evalSplineBasis1D( variate = 1, elem_extraction_operator = C, basis_idx = 1, domain = elem_domain ), 1.0 )
       
#     def test_linear_basis_biunit_domain( self ):
#         elem_degree = 1
#         elem_domain = [-1, 1]
#         C = np.eye( elem_degree+ 1 )
#         self.assertAlmostEqual( evalSplineBasis1D( variate = -1, elem_extraction_operator = C, basis_idx = 0, domain = elem_domain ), 1.0 )
#         self.assertAlmostEqual( evalSplineBasis1D( variate = -1, elem_extraction_operator = C, basis_idx = 1, domain = elem_domain ), 0.0 )
#         self.assertAlmostEqual( evalSplineBasis1D( variate = 1, elem_extraction_operator = C, basis_idx = 0, domain = elem_domain ), 0.0 )
#         self.assertAlmostEqual( evalSplineBasis1D( variate = 1, elem_extraction_operator = C, basis_idx = 1, domain = elem_domain ), 1.0 )
#     def test_linear_basis_half_unit_domain( self ):
#         elem_degree = 1
#         elem_domain = [0.,0.5]
#         C = np.eye( elem_degree+ 1 )
#         self.assertAlmostEqual( evalSplineBasis1D( variate = 0, elem_extraction_operator = C, basis_idx = 0, domain = elem_domain ), 1.0 )
#         self.assertAlmostEqual( evalSplineBasis1D( variate = 0, elem_extraction_operator = C, basis_idx = 1, domain = elem_domain ), 0.0 )
#         self.assertAlmostEqual( evalSplineBasis1D( variate = 0.5, elem_extraction_operator = C, basis_idx = 0, domain = elem_domain ), 0.0 )
#         self.assertAlmostEqual( evalSplineBasis1D( variate = 0.5, elem_extraction_operator = C, basis_idx = 1, domain = elem_domain ), 1.0 )
       
# unittest.main()