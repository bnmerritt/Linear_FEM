import numpy
import unittest
import basis 
import quadrature

class Test_assembleGramMatrix( unittest.TestCase ):
    def test_quadratic_legendre( self ):
        test_gram_matrix = assembleGramMatrix( domain = [0, 1], degree = 2, solution_basis = basis.evalLegendreBasis1D )
        gold_gram_matrix = numpy.array( [ [1.0, 0.0, 0.0], [0.0, 1.0/3.0, 0.0], [0.0, 0.0, 0.2] ] )
        self.assertTrue( numpy.allclose( test_gram_matrix, gold_gram_matrix ) )
    
    def test_cubic_legendre( self ):
        test_gram_matrix = assembleGramMatrix( domain = [0, 1], degree = 3, solution_basis = basis.evalLegendreBasis1D )
        gold_gram_matrix = numpy.array( [ [1.0, 0.0, 0.0, 0.0], [0.0, 1.0/3.0, 0.0, 0.0], [0.0, 0.0, 0.2, 0.0], [ 0.0, 0.0, 0.0, 1.0/7.0] ] )
        self.assertTrue( numpy.allclose( test_gram_matrix, gold_gram_matrix ) )

    def test_linear_bernstein( self ):
        test_gram_matrix = assembleGramMatrix( domain = [0, 1], degree = 1, solution_basis = basis.evalBernsteinBasis1D )
        gold_gram_matrix = numpy.array( [ [1.0/3.0, 1.0/6.0], [1.0/6.0, 1.0/3.0] ] )
        self.assertTrue( numpy.allclose( test_gram_matrix, gold_gram_matrix ) )

    def test_quadratic_bernstein( self ):
        test_gram_matrix = assembleGramMatrix( domain = [0, 1], degree = 2, solution_basis = basis.evalBernsteinBasis1D )
        gold_gram_matrix = numpy.array( [ [0.2, 0.1, 1.0/30.0], [0.1, 2.0/15.0, 0.1], [1.0/30.0, 0.1, 0.2] ] )
        self.assertTrue( numpy.allclose( test_gram_matrix, gold_gram_matrix ) )
    
    def test_cubic_bernstein( self ):
        test_gram_matrix = assembleGramMatrix( domain = [0, 1], degree = 3, solution_basis = basis.evalBernsteinBasis1D )
        gold_gram_matrix = numpy.array( [ [1.0/7.0, 1.0/14.0, 1.0/35.0, 1.0/140.0], [1.0/14.0, 3.0/35.0, 9.0/140.0, 1.0/35.0], [1.0/35.0, 9.0/140.0, 3.0/35.0, 1.0/14.0], [ 1.0/140.0, 1.0/35.0, 1.0/14.0, 1.0/7.0] ] )
        self.assertTrue( numpy.allclose( test_gram_matrix, gold_gram_matrix ) )

    def test_linear_lagrange( self ):
        test_gram_matrix = assembleGramMatrix( domain = [0, 1], degree = 1, solution_basis = basis.evalLagrangeBasis1D )
        gold_gram_matrix = numpy.array( [ [1.0/3.0, 1.0/6.0], [1.0/6.0, 1.0/3.0] ] )
        self.assertTrue( numpy.allclose( test_gram_matrix, gold_gram_matrix ) )

    def test_quadratic_lagrange( self ):
        test_gram_matrix = assembleGramMatrix( domain = [0, 1], degree = 2, solution_basis = basis.evalLagrangeBasis1D )
        gold_gram_matrix = numpy.array( [ [2.0/15.0, 1.0/15.0, -1.0/30.0], [1.0/15.0, 8.0/15.0, 1.0/15.0], [-1.0/30.0, 1.0/15.0, 2.0/15.0] ] )
        self.assertTrue( numpy.allclose( test_gram_matrix, gold_gram_matrix ) )
    
    def test_cubic_lagrange( self ):
        test_gram_matrix = assembleGramMatrix( domain = [0, 1], degree = 3, solution_basis = basis.evalLagrangeBasis1D )
        gold_gram_matrix = numpy.array( [ [8.0/105.0, 33.0/560.0, -3.0/140.0, 19.0/1680.0], [33.0/560.0, 27.0/70.0, -27.0/560.0, -3.0/140.0], [-3.0/140.0, -27.0/560.0, 27.0/70.0, 33/560.0], [ 19.0/1680.0, -3.0/140.0, 33.0/560.0, 8.0/105.0] ] )

def assembleGramMatrix(domain,degree,solution_basis):
    nodes = degree + 1
    qp, w = quadrature.computeGaussLegendreQuadrature(nodes)
    num_basis_vec = degree + 1
    derivative = (domain[-1] - domain[0]) / 2
    M = numpy.zeros((num_basis_vec,num_basis_vec))
    for A in range(0,num_basis_vec):
        for B in range(0,num_basis_vec):
            for k in range(0,len(qp)):
                M[A,B] += solution_basis(qp[k],degree,A,domain) * solution_basis(qp[k],degree,B,domain) * w[k] * derivative   
    return M

# unittest.main()

# degree = 2
# nodes = degree + 1
# qp, w = quadrature.computeGaussLegendreQuadrature(nodes)
# domain= [0,1]
# solution_basis = basis.evalBernsteinBasis1D
# num_basis_vec = degree + 1
# derivative = (domain[-1] - domain[0]) / 2
# M = numpy.zeros((num_basis_vec,num_basis_vec))
# for A in range(0,num_basis_vec):
#     for B in range(0,num_basis_vec):
#         for k in range(0,len(qp)):
#             M[A,B] += solution_basis(qp[k],degree,A,domain) * solution_basis(qp[k],degree,B,domain) * w[k] * derivative 