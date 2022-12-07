import numpy
import unittest
import basis 
import mesh_adapted
# import computeFitError
import scipy
import readBEXT_JSON
import basis
import quadrature
from scipy import integrate
import matplotlib.pyplot as plt
import assembleForceVector_Stiffness
import assembleGramMatrix_Stiffness

class test_ComputeSolution( unittest.TestCase ):
    def test_simple( self ):
           problem = { "elastic_modulus": 100,
                       "area": 0.01,
                       "length": 1.0,
                       "traction": { "value": 1e-3, "position": 1.0 },
                       "displacement": { "value": 0.0, "position": 0.0 },
                       "body_force": 1e-3 }
           # spline_space = { "domain": [0, problem[ "length" ]], "degree": [ 1, 1, 1 ], "continuity": [ -1, 0, 0, -1 ] }
           # uspline.make_uspline_mesh( spline_space, "temp_uspline" )
           # uspline_bext = bext.readBEXT( "temp_uspline.json" )
           uspline_bext = readBEXT_JSON.readBEXT( "uspline_Final1.json" )
           test_sol_coeff = computeSolution( problem = problem, uspline_bext = uspline_bext )
           gold_sol_coeff = numpy.array( [ 0.0, 11.0 / 18000.0, 1.0 / 900.0, 3.0 / 2000.0 ] )
           self.assertTrue( numpy.allclose( test_sol_coeff, gold_sol_coeff ) )
           # splineBarGalerkin.plotSolution( test_sol_coeff, uspline_bext )
           # splineBarGalerkin.plotCompareFunToExactSolution( problem, test_sol_coeff, uspline_bext )

    def test_textbook_problem( self ):
           problem = { "elastic_modulus": 200e9,
                       "area": 1.0,
                       "length": 5.0,
                       "traction": { "value": 9810.0, "position": 5.0 },
                       "displacement": { "value": 0.0, "position": 0.0 },
                       "body_force": 784800.0 }
           # spline_space = { "domain": [0, problem[ "length" ]], "degree": [ 2, 2 ], "continuity": [ -1, 1, -1 ] }
           # uspline.make_uspline_mesh( spline_space, "temp_uspline" )
           # uspline_bext = bext.readBEXT( "temp_uspline.json" )
           uspline_bext = readBEXT_JSON.readBEXT( "uspline_Final2.json" )
           test_sol_coeff = computeSolution( problem = problem, uspline_bext = uspline_bext )
           gold_sol_coeff = numpy.array( [0.0, 2.45863125e-05, 4.92339375e-05, 4.92952500e-05] )
           self.assertTrue( numpy.allclose( test_sol_coeff, gold_sol_coeff ) )
           # splineBarGalerkin.plotSolution( test_sol_coeff, uspline_bext )
           # splineBarGalerkin.plotCompareFunToExactSolution( problem, test_sol_coeff, uspline_bext )
        
## MAIN CODE
def computeSolution( problem, uspline_bext ):
    M = assembleStiffnessMatrix( problem, uspline_bext )
    F = assembleForceVector( problem, uspline_bext )
    M_final, F_final = applyDisplacement( problem, M, F, uspline_bext )
    M_final_inv = numpy.linalg.inv(M_final)
    coeff = numpy.dot(M_final_inv,F_final)
    coeff_final = numpy.append(problem["displacement"]["value"],coeff)
    return coeff_final

# def assembleSolution( coeff, problem, uspline_bext ):
#     # Your code goes here
#     return coeff

def applyDisplacement( problem, stiffness_matrix, force_vector, uspline_bext ):
    value = problem["displacement"]["value"]
    position = problem["displacement"]["position"]
    stiffness_matrix_displaced = numpy.delete(numpy.delete(stiffness_matrix,0,0),0,1)
    force_vector_displaced = numpy.delete(force_vector,0,0)
    return stiffness_matrix_displaced, force_vector_displaced

# def applyTraction( problem, uspline_bext ):
#     force_tractions = assembleForceVector_Stiffness.applyTractions(problem,uspline_bext)
#     return force_tractions

# def evaluateConstitutiveModel( problem ):
#     # Your code goes here
#     return

def assembleStiffnessMatrix( problem, uspline_bext ):
    stiffness_matrix = assembleGramMatrix_Stiffness.assembleStiffnessMatrix(problem, uspline_bext)
    return stiffness_matrix

def assembleForceVector( problem, uspline_bext ):
    force_vector = assembleForceVector_Stiffness.assembleForceVector(problem,uspline_bext)
    return force_vector

## UTILITY CODE
# def evaluateSolutionAt( x, coeff, uspline_bext ):
#     elem_id = bext.getElementIdContainingPoint( uspline_bext, x )
#     elem_nodes = bext.getElementNodeIds( uspline_bext, elem_id )
#     elem_domain = bext.getElementDomain( uspline_bext, elem_id )
#     elem_degree = bext.getElementDegree( uspline_bext, elem_id )
#     elem_extraction_operator = bext.getElementExtractionOperator( uspline_bext, elem_id )
#     sol = 0.0
#     for n in range( 0, len( elem_nodes ) ):
#         curr_node = elem_nodes[n]
#         sol += coeff[curr_node] * basis.evalSplineBasis1D( extraction_operator = elem_extraction_operator, basis_idx = n, domain = elem_domain, variate = x )
#     return sol

# def computeElementFitError( problem, coeff, uspline_bext, elem_id ):
#     domain = bext.getDomain( uspline_bext )
#     elem_domain = bext.getElementDomain( uspline_bext, elem_id )
#     elem_degree = bext.getElementDegree( uspline_bext, elem_id )
#     num_qp = int( numpy.ceil( ( 2*(elem_degree - 1) + 1 ) / 2.0 ) + 1 )
#     abs_err_fun = lambda x : abs( evaluateExactSolutionAt( problem, basis.affine_mapping_1D( [-1, 1], elem_domain, x ) ) - evaluateSolutionAt( basis.affine_mapping_1D( [-1, 1], elem_domain, x ), coeff, uspline_bext ) )
#     abs_error = quadrature.quad( abs_err_fun, elem_domain, num_qp )
#     return abs_error

# def computeFitError( problem, coeff, uspline_bext ):
#     num_elems = bext.getNumElems( uspline_bext )
#     abs_error = 0.0
#     for elem_idx in range( 0, num_elems ):
#         elem_id = bext.elemIdFromElemIdx( uspline_bext, elem_idx )
#         abs_error += computeElementFitError( problem, coeff, uspline_bext, elem_id )
#     domain = bext.getDomain( uspline_bext )
#     target_fun_norm, _ = scipy.integrate.quad( lambda x: abs( evaluateExactSolutionAt( problem, x ) ), domain[0], domain[1], epsrel = 1e-12, limit = num_elems * 100 )
#     rel_error = abs_error / target_fun_norm
#     return abs_error, rel_error

# def plotCompareGoldTestSolution( gold_coeff, test_coeff, uspline_bext ):
#     domain = bext.getDomain( uspline_bext )
#     x = numpy.linspace( domain[0], domain[1], 1000 )
#     yg = numpy.zeros( 1000 )
#     yt = numpy.zeros( 1000 )
#     for i in range(0, len(x) ):
#         yg[i] = evaluateSolutionAt( x[i], test_coeff, uspline_bext )
#         yt[i] = evaluateSolutionAt( x[i], gold_coeff, uspline_bext )
#     plt.plot( x, yg )
#     plt.plot( x, yt )
#     plt.show()

# def plotCompareFunToExactSolution( problem, test_coeff, uspline_bext ):
#     domain = bext.getDomain( uspline_bext )
#     x = numpy.linspace( domain[0], domain[1], 1000 )
#     ya = numpy.zeros( 1000 )
#     ye = numpy.zeros( 1000 )
#     for i in range(0, len(x) ):
#         ya[i] = evaluateSolutionAt( x[i], test_coeff, uspline_bext )
#         ye[i] = evaluateExactSolutionAt( problem, x[i] )
#     plt.plot( x, ya )
#     plt.plot( x, ye )
#     plt.show()

# def computeConvergenceRate( num_entities, qoi ):
#     def func( x, a, b, c ):
#         return a * numpy.power( x, b ) + c
#     fit = scipy.optimize.curve_fit(func, num_entities, qoi, method=’trf’, bounds = ([-numpy.inf, -numpy.inf, -numpy.inf ], [numpy.inf, 0.0, numpy.inf]) )
#     a,b,c = fit[0]
#     return b

# def plotSolution( sol_coeff, uspline_bext ):
#     domain = bext.getDomain( uspline_bext )
#     x = numpy.linspace( domain[0], domain[1], 1000 )
#     y = numpy.zeros( 1000 )
#     for i in range(0, len(x) ):
#         y[i] = evaluateSolutionAt( x[i], sol_coeff, uspline_bext )
#     plt.plot( x, y )
#     plt.plot( bext.getSplineNodes( uspline_bext )[:,0], sol_coeff, color = "k", marker = "o", markerfacecolor = "k" )
#     plt.show()

# def evaluateExactSolutionAt( problem, x ):
#     term_1 = problem[ "traction" ][ "value" ] / evaluateConstitutiveModel( problem ) * x
#     term_2 = problem[ "displacement" ][ "value" ]
#     term_3 =  ( ( problem[ "length" ]**2.0 * problem[ "body_force" ] / 2 ) / evaluateConstitutiveModel( problem ) ) - ( ( ( problem[ "length" ] - x )**2.0 * problem[ "body_force" ] / 2 ) / evaluateConstitutiveModel( problem ) )
#     sol = term_1 + term_2 + term_3
#     return sol

# def plotExactSolution( problem ):
#     domain = [0, problem[ "length" ] ]
#     x = numpy.linspace( domain[0], domain[1], 1000 )
#     y = numpy.zeros( 1000 )
#     for i in range(0, len(x) ):
#         y[i] = evaluateExactSolutionAt( problem, x[i] )
#     plt.plot( x, y )
#     plt.show()

unittest.main()

# problem = { "elastic_modulus": 200e9,
#                        "area": 1.0,
#                        "length": 5.0,
#                        "traction": { "value": 9810.0, "position": 5.0 },
#                        "displacement": { "value": 0.0, "position": 0.0 },
#                        "body_force": 784800.0 }
# uspline_bext = readBEXT_JSON.readBEXT( "uspline_Final2.json" )
# M = assembleStiffnessMatrix( problem, uspline_bext )
# F = assembleForceVector( problem, uspline_bext )
# M_final, F_final = applyDisplacement( problem, M, F, uspline_bext )
# M_final_inv = numpy.linalg.inv(M_final)
# coeff = numpy.dot(M_final_inv,F_final)
