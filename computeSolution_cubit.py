import numpy
import unittest
import basis 
import mesh_adapted
# import computeFitError
import scipy
import assembleGramMatrix_cubit
import assembleForceVector_cubit
import readBEXT_JSON
import basis
import quadrature
from scipy import integrate
import matplotlib.pyplot as plt

class Test_computeSolution( unittest.TestCase ):
    def test_cubic_polynomial_target_linear_bspline( self ):
        # print( "POLY TEST" )
        target_fun = lambda x: x**3 - (8/5)*x**2 + (3/5)*x
        # spline_space = { "domain": [0, 1], "degree": [ 1, 1 ], "continuity": [ -1, 0, -1 ] }
        # uspline.make_uspline_mesh( spline_space, "temp_uspline" )
        # uspline_bext = bext.readBEXT( "temp_uspline.json" )
        uspline_bext = readBEXT_JSON.readBEXT( "test_usplineCS1.json" )
        test_sol_coeff = computeSolution( target_fun = target_fun, uspline_bext = uspline_bext )
        plotCompareFunToTestSolution( target_fun, test_sol_coeff, uspline_bext )
        gold_sol_coeff = numpy.array( [ 9.0/160.0, 7.0/240.0, -23.0/480.0 ] )
        abs_err, rel_err = computeFitError( target_fun, test_sol_coeff, uspline_bext )
        self.assertTrue( numpy.allclose( gold_sol_coeff, test_sol_coeff ) )
        self.assertAlmostEqual( first = rel_err, second = 0, delta = 1e0 )

    def test_cubic_polynomial_target_quadratic_bspline( self ):
        # print( "POLY TEST" )
        target_fun = lambda x: x**3 - (8/5)*x**2 + (3/5)*x
        # spline_space = { "domain": [0, 1], "degree": [ 2, 2 ], "continuity": [ -1, 1, -1 ] }
        # uspline.make_uspline_mesh( spline_space, "temp_uspline" )
        # uspline_bext = bext.readBEXT( "temp_uspline.json" )
        uspline_bext = readBEXT_JSON.readBEXT( "test_usplineCS2.json" )
        test_sol_coeff = computeSolution( target_fun = target_fun, uspline_bext = uspline_bext )
        plotCompareFunToTestSolution( target_fun, test_sol_coeff, uspline_bext )
        gold_sol_coeff = numpy.array( [ 1.0/120.0, 9.0/80.0, -1.0/16.0, -1.0/120.0 ] )
        abs_err, rel_err = computeFitError( target_fun, test_sol_coeff, uspline_bext )
        self.assertTrue( numpy.allclose( gold_sol_coeff, test_sol_coeff ) )
        self.assertAlmostEqual( first = rel_err, second = 0, delta = 1e-1 )

    def test_cubic_polynomial_target_cubic_bspline( self ):
        # print( "POLY TEST" )
        target_fun = lambda x: x**3 - (8/5)*x**2 + (3/5)*x
        # spline_space = { "domain": [0, 1], "degree": [ 3, 3 ], "continuity": [ -1, 2, -1 ] }
        # uspline.make_uspline_mesh( spline_space, "temp_uspline" )
        # uspline_bext = bext.readBEXT( "temp_uspline.json" )
        uspline_bext = readBEXT_JSON.readBEXT( "test_usplineCS34.json" )
        test_sol_coeff = computeSolution( target_fun = target_fun, uspline_bext = uspline_bext )
        plotCompareFunToTestSolution( target_fun, test_sol_coeff, uspline_bext )
        gold_sol_coeff = numpy.array( [ 0.0, 1.0/10.0, 1.0/30.0, -1.0/15.0, 0.0 ] )
        abs_err, rel_err = computeFitError( target_fun, test_sol_coeff, uspline_bext )
        self.assertTrue( numpy.allclose( gold_sol_coeff, test_sol_coeff ) )
        self.assertAlmostEqual( first = rel_err, second = 0, delta = 1e-12 )

    def test_sin_target( self ):
        # print( "SIN TEST" )
        target_fun = lambda x: numpy.sin( numpy.pi * x )
        # spline_space = { "domain": [0, 1], "degree": [ 3, 3 ], "continuity": [ -1, 2, -1 ] }
        # uspline.make_uspline_mesh( spline_space, "temp_uspline" )
        # uspline_bext = bext.readBEXT( "temp_uspline.json" )
        uspline_bext = readBEXT_JSON.readBEXT( "test_usplineCS34.json" )
        test_sol_coeff = computeSolution( target_fun = target_fun, uspline_bext = uspline_bext )
        abs_err, rel_err = computeFitError( target_fun, test_sol_coeff, uspline_bext )
        plotCompareFunToTestSolution( target_fun, test_sol_coeff, uspline_bext )
        self.assertAlmostEqual( first = rel_err, second = 0, delta = 1e-2 )

    def test_erfc_target( self ):
        # print( "ERFC TEST" )
        target_fun = lambda x: numpy.real( scipy.special.erfc( x ) )
        # spline_space = { "domain": [-1, 1], "degree": [ 3, 1, 3 ], "continuity": [ -1, 1, 1, -1 ] }
        # uspline.make_uspline_mesh( spline_space, "temp_uspline" )
        # uspline_bext = bext.readBEXT( "temp_uspline.json" )
        uspline_bext = readBEXT_JSON.readBEXT( "test_usplineCS5.json" )
        test_sol_coeff = computeSolution( target_fun = target_fun, uspline_bext = uspline_bext )
        abs_err, rel_err = computeFitError( target_fun, test_sol_coeff, uspline_bext )
        plotCompareFunToTestSolution( target_fun, test_sol_coeff, uspline_bext )
        self.assertAlmostEqual( first = rel_err, second = 0, delta = 1e-2 )

    def test_exptx_target( self ):
        # print( "EXPT TEST" )
        target_fun = lambda x: float( numpy.real( float( x )**float( x ) ) )
        # spline_space = { "domain": [-1, 1], "degree": [ 5, 5, 5, 5 ], "continuity": [ -1, 4, 0, 4, -1 ] }
        # uspline.make_uspline_mesh( spline_space, "temp_uspline" )
        # uspline_bext = bext.readBEXT( "temp_uspline.json" )
        uspline_bext = readBEXT_JSON.readBEXT( "test_usplineCS6.json" )
        test_sol_coeff = computeSolution( target_fun = target_fun, uspline_bext = uspline_bext )
        abs_err, rel_err = computeFitError( target_fun, test_sol_coeff, uspline_bext )
        plotCompareFunToTestSolution( target_fun, test_sol_coeff, uspline_bext )
        self.assertAlmostEqual( first = rel_err, second = 0, delta = 1e-2 )
        
def computeFitError( target_fun, coeff, uspline_bext ):
    num_elems = readBEXT_JSON.getNumElems( uspline_bext )
    abs_error = 0.0
    for elem_idx in range( 0, num_elems ):
        elem_id = readBEXT_JSON.elemIdFromElemIdx( uspline_bext, elem_idx )
        abs_error += computeElementFitError( target_fun, coeff, uspline_bext, elem_id )
    domain = readBEXT_JSON.getDomain( uspline_bext )
    target_fun_norm, _ = scipy.integrate.quad( lambda x: abs( target_fun(x) ), domain[0], domain[1], epsrel = 1e-12, limit = num_elems * 100 )
    rel_error = abs_error / target_fun_norm
    return abs_error, rel_error

def computeElementFitError( target_fun, coeff, uspline_bext, elem_id ):
    domain = readBEXT_JSON.getDomain( uspline_bext )
    elem_domain = readBEXT_JSON.getElementDomain( uspline_bext, elem_id )
    elem_degree = readBEXT_JSON.getElementDegree( uspline_bext, elem_id )
    num_qp = int( numpy.ceil( ( 2*elem_degree + 1 ) / 2.0 ) + 1 )
    abs_err_fun = lambda x : abs( target_fun( basis.affine_mapping_1D( [-1, 1], elem_domain, x ) ) - evaluateSolutionAt( basis.affine_mapping_1D( [-1, 1], elem_domain, x ), coeff, uspline_bext ) )
    abs_error = quadrature.quad( abs_err_fun, elem_domain, num_qp )
    return abs_error

def evaluateSolutionAt( x, coeff, uspline_bext ):
    elem_id = readBEXT_JSON.getElementIdContainingPoint( uspline_bext, x )
    elem_nodes = readBEXT_JSON.getElementNodeIds( uspline_bext, elem_id )
    elem_domain = readBEXT_JSON.getElementDomain( uspline_bext, elem_id )
    elem_degree = readBEXT_JSON.getElementDegree( uspline_bext, elem_id )
    elem_extraction_operator = readBEXT_JSON.getElementExtractionOperator( uspline_bext, elem_id )
    sol = 0.0
    for n in range( 0, len( elem_nodes ) ):
        curr_node = elem_nodes[n]
        sol += coeff[curr_node] * basis.evalSplineBasis1D(  x , elem_extraction_operator, n, elem_domain)
    return sol

def plotCompareFunToTestSolution( fun, test_coeff, uspline_bext ):
    nodes = readBEXT_JSON.getSplineNodes(uspline_bext)[:,0]
    domain = readBEXT_JSON.getDomain( uspline_bext )
    x = numpy.linspace( domain[0], domain[1], 1000 )
    y = numpy.zeros( 1000 )
    yt = numpy.zeros( 1000 )
    for i in range(0, len(x) ):
        y[i] = fun( x[i] )
        yt[i] = evaluateSolutionAt( x[i], test_coeff, uspline_bext )
    plt.plot( x, y )
    plt.plot( x, yt )
    plt.plot(nodes,test_coeff)
    plt.show()

        
def computeSolution(target_fun,uspline_bext):
    M = assembleGramMatrix_cubit.assembleGramMatrix(uspline_bext)
    F = assembleForceVector_cubit.assembleForceVector(target_fun,uspline_bext)
    M_inv = numpy.linalg.inv(M)
    d = numpy.dot(M_inv,F)
    return d

unittest.main()

# target_fun = lambda x: x**3 - (8/5)*x**2 + (3/5)*x
# uspline_bext = readBEXT_JSON.readBEXT( "test_usplineCS1.json" )
# M = assembleGramMatrix_cubit.assembleGramMatrix(uspline_bext)
# F = assembleForceVector_cubit.assembleForceVector(target_fun,uspline_bext)
# M_inv = numpy.linalg.inv(M)
# d = numpy.dot(M_inv,F)
# abs_error, rel_error = computeFitError( target_fun, d, uspline_bext )