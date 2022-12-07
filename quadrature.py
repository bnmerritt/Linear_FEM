import scipy
from scipy import optimize
import numpy
import unittest
# from src import basis
import math
import momentFitting
import basis

def computeGaussLegendreQuadrature( n ):
    M = numpy.zeros( 2*n, dtype = "double" )
    M[0] = 2.0
    x0 = numpy.linspace( -1, 1, n )
    sol = scipy.optimize.least_squares( lambda x : objFun( M, x ), x0, bounds = (-1, 1), ftol = 1e-14, xtol = 1e-14, gtol = 1e-14 )
    qp = sol.x
    w = solveLinearMomentFit( M, qp )
    return qp, w

def quad( fun, domain, num_points ):
    jacobian = ( domain[1] - domain[0] ) / ( 1 - (-1) )
    x_qp, w_qp = getGaussLegendreQuadrature( num_points )
    integral = 0.0
    for qp in range( 0, len( x_qp ) ):
        integral += ( fun( x_qp[qp] ) * w_qp[qp] ) * jacobian
    return integral


def getGaussLegendreQuadrature( num_points ):
    if num_points == 1:
        x = [ 0.0 ]
        w = [ 2.0 ]
    elif num_points == 2:
        x = [ -1.0 / math.sqrt(3), 
              +1.0 / math.sqrt(3) ]

        w = [ 1.0, 
              1.0  ]
    elif num_points == 3:
        x = [ -1.0 * math.sqrt( 3.0 / 5.0 ), 
               0.0, 
              +1.0 * math.sqrt( 3.0 / 5.0 ) ]

        w = [ 5.0 / 9.0, 
              8.0 / 9.0, 
              5.0 / 9.0 ]
    elif num_points == 4:
        x = [ -1.0 * math.sqrt( 3.0 / 7.0 + 2.0 / 7.0 * math.sqrt( 6.0 / 5.0 ) ),
              -1.0 * math.sqrt( 3.0 / 7.0 - 2.0 / 7.0 * math.sqrt( 6.0 / 5.0 ) ),
              +1.0 * math.sqrt( 3.0 / 7.0 - 2.0 / 7.0 * math.sqrt( 6.0 / 5.0 ) ),
              +1.0 * math.sqrt( 3.0 / 7.0 + 2.0 / 7.0 * math.sqrt( 6.0 / 5.0 ) ) ]
        
        w = [ ( 18.0 - math.sqrt( 30.0 ) ) / 36.0,
              ( 18.0 + math.sqrt( 30.0 ) ) / 36.0,
              ( 18.0 + math.sqrt( 30.0 ) ) / 36.0,
              ( 18.0 - math.sqrt( 30.0 ) ) / 36.0 ]
    elif num_points == 5:
        x = [ -1.0 / 3.0 * math.sqrt( 5.0 + 2.0 * math.sqrt( 10.0 / 7.0 ) ),
              -1.0 / 3.0 * math.sqrt( 5.0 - 2.0 * math.sqrt( 10.0 / 7.0 ) ),
               0.0,
              +1.0 / 3.0 * math.sqrt( 5.0 - 2.0 * math.sqrt( 10.0 / 7.0 ) ),
              +1.0 / 3.0 * math.sqrt( 5.0 + 2.0 * math.sqrt( 10.0 / 7.0 ) ) ]
        
        w = [ ( 322.0 - 13.0 * math.sqrt( 70.0 ) ) / 900.0,
              ( 322.0 + 13.0 * math.sqrt( 70.0 ) ) / 900.0,
                128.0 / 225.0,
              ( 322.0 + 13.0 * math.sqrt( 70.0 ) ) / 900.0,
              ( 322.0 - 13.0 * math.sqrt( 70.0 ) ) / 900.0, ]
    elif num_points > 5:
        # x, w = momentFitting.computeQuadrature( num_points, [-1, 1], basis.evalLegendreBasis1D )
        # x, w = computeGaussLegendreQuadratureRule( num_points )
        # A = momentFitting.assembleLinearMomentFitSystem( num_points, basis.symLegendreBasis, [-1, 1], x )
        x = basis.eigenvaluesLegendreBasis( num_points )
        M = momentFitting.computeMomentVector( num_points, basis.symLegendreBasis, [-1, 1] )
        w = momentFitting.solveLinearMomentFit( M, basis.symLegendreBasis, [-1, 1], x )
    else:
        raise( Exception( "num_points_MUST_BE_POSITIVE_INTEGER" ) )
    return x, w

def assembleLinearMomentFitSystem( degree, pts ):
    A = numpy.zeros( shape = ( degree + 1, len( pts ) ), dtype = "double" )
    for m in range(0,degree + 1):
        for n in range(0,len(pts)):
            A[m,n] = evalLegendreBasis1D(m,pts[n])
    return A

def solveLinearMomentFit( M, pts ):
    degree = len( M ) - 1
    A = assembleLinearMomentFitSystem( degree, pts )
    sol = scipy.optimize.lsq_linear( A, M )
    w = sol.x
    return w

def objFun( M, pts ):
    degree = len( M ) - 1
    A = assembleLinearMomentFitSystem( degree, pts )
    w = solveLinearMomentFit( M, pts )
    obj_val = M - numpy.matmul(A,w)
    return obj_val

def evalLegendreBasis1D(degree,variate):
    if degree == 0:
        val = 1
    elif degree == 1:
        val = variate
    else:
        i = degree - 1
        val = ((i + 1)**(-1)) * ((2*i+1)*variate*evalLegendreBasis1D(i,variate) - i*evalLegendreBasis1D(i-1,variate))
    return val

# class Test_computeGaussLegendreQuadrature( unittest.TestCase ):
#     def test_1_pt( self ):
#         qp_gold = numpy.array( [ 0.0 ] )
#         w_gold = numpy.array( [ 2.0 ] )
#         [ qp, w ] = computeGaussLegendreQuadrature( 1 )
#         self.assertAlmostEqual( first = qp, second = qp_gold, delta = 1e-12 )
#         self.assertAlmostEqual( first = w, second = w_gold, delta = 1e-12 )

#     def test_2_pt( self ):
#         qp_gold = numpy.array( [ -1.0/numpy.sqrt(3), 1.0/numpy.sqrt(3) ] )
#         w_gold = numpy.array( [ 1.0, 1.0 ] )
#         [ qp, w ] = computeGaussLegendreQuadrature( 2 )
#         self.assertTrue( numpy.allclose( qp, qp_gold ) )
#         self.assertTrue( numpy.allclose( w, w_gold ) )

#     def test_3_pt( self ):
#         qp_gold = numpy.array( [ -1.0 * numpy.sqrt( 3.0 / 5.0 ),
#                                 0.0,
#                                 +1.0 * numpy.sqrt( 3.0 / 5.0 ) ] )
#         w_gold = numpy.array( [ 5.0 / 9.0,
#                                 8.0 / 9.0,
#                                 5.0 / 9.0 ] )
#         [ qp, w ] = computeGaussLegendreQuadrature( 3 )
#         self.assertTrue( numpy.allclose( qp, qp_gold ) )
#         self.assertTrue( numpy.allclose( w, w_gold ) )

#     def test_4_pt( self ):
#         qp_gold = numpy.array( [ -1.0 * numpy.sqrt( 3.0 / 7.0 + 2.0 / 7.0 * numpy.sqrt( 6.0 / 5.0 ) ),
#                                 -1.0 * numpy.sqrt( 3.0 / 7.0 - 2.0 / 7.0 * numpy.sqrt( 6.0 / 5.0 ) ),
#                                 +1.0 * numpy.sqrt( 3.0 / 7.0 - 2.0 / 7.0 * numpy.sqrt( 6.0 / 5.0 ) ),
#                                 +1.0 * numpy.sqrt( 3.0 / 7.0 + 2.0 / 7.0 * numpy.sqrt( 6.0 / 5.0 ) ) ] )
#         w_gold = numpy.array( [ ( 18.0 - numpy.sqrt( 30.0 ) ) / 36.0,
#                                 ( 18.0 + numpy.sqrt( 30.0 ) ) / 36.0,
#                                 ( 18.0 + numpy.sqrt( 30.0 ) ) / 36.0,
#                                 ( 18.0 - numpy.sqrt( 30.0 ) ) / 36.0 ] )
#         [ qp, w ] = computeGaussLegendreQuadrature( 4 )
#         self.assertTrue( numpy.allclose( qp, qp_gold ) )
#         self.assertTrue( numpy.allclose( w, w_gold ) )

#     def test_5_pt( self ):
#         qp_gold = numpy.array( [ -1.0 / 3.0 * numpy.sqrt( 5.0 + 2.0 * numpy.sqrt( 10.0 / 7.0 ) ),
#                                 -1.0 / 3.0 * numpy.sqrt( 5.0 - 2.0 * numpy.sqrt( 10.0 / 7.0 ) ),
#                                 0.0,
#                                 +1.0 / 3.0 * numpy.sqrt( 5.0 - 2.0 * numpy.sqrt( 10.0 / 7.0 ) ),
#                                 +1.0 / 3.0 * numpy.sqrt( 5.0 + 2.0 * numpy.sqrt( 10.0 / 7.0 ) ) ] )
#         w_gold = numpy.array( [ ( 322.0 - 13.0 * numpy.sqrt( 70.0 ) ) / 900.0,
#                                 ( 322.0 + 13.0 * numpy.sqrt( 70.0 ) ) / 900.0,
#                                 128.0 / 225.0,
#                                 ( 322.0 + 13.0 * numpy.sqrt( 70.0 ) ) / 900.0,
#                                 ( 322.0 - 13.0 * numpy.sqrt( 70.0 ) ) / 900.0, ] )
#         [ qp, w ] = computeGaussLegendreQuadrature( 5 )
#         self.assertTrue( numpy.allclose( qp, qp_gold ) )
#         self.assertTrue( numpy.allclose( w, w_gold ) )
        
# unittest.main()  