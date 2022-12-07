import argparse
import sympy
import uspline
import bext
import quadrature

## MAIN FUNCTION
def main( target_fun, spline_space ):
    filename = "temp_uspline"
    uspline.make_uspline_mesh( spline_space, filename )
    uspline_bext = bext.readBEXT( filename + ".json" )
    sol = computeSolution( target_fun, uspline_bext )
    return sol

## SECONDARY FUNCTIONS
def computeSolution():
    # YOUR CODE GOES HERE
    return sol

def assembleGramMatrix():
    # YOUR CODE GOES HERE
    return gram_matrix

def assembleForceVector():
    # YOUR CODE GOES HERE
    return force_vector

## UTILITY FUNCTIONS
def evaluateSolutionAt( x, coeff, uspline_bext ):
    elem_id = bext.getElementIdContainingPoint( uspline_bext, x )
    elem_nodes = bext.getElementNodeIds( uspline_bext, elem_id )
    elem_domain = bext.getElementDomain( uspline_bext, elem_id )
    elem_degree = bext.getElementDegree( uspline_bext, elem_id )
    elem_extraction_operator = bext.getElementExtractionOperator( uspline_bext, elem_id )
    sol = 0.0
    for n in range( 0, len( elem_nodes ) ):
        curr_node = elem_nodes[n]
        sol += coeff[curr_node] * basis.evalSplineBasis1D( extraction_operator = elem_extraction_operator, basis_idx = n, domain = elem_domain, variate = x )
    return sol

def computeElementFitError( target_fun, coeff, uspline_bext, elem_id ):
    domain = bext.getDomain( uspline_bext )
    elem_domain = bext.getElementDomain( uspline_bext, elem_id )
    elem_degree = bext.getElementDegree( uspline_bext, elem_id )
    num_qp = int( numpy.ceil( ( 2*elem_degree + 1 ) / 2.0 ) + 1 )
    abs_err_fun = lambda x : abs( target_fun( basis.affine_mapping_1D( [-1, 1], elem_domain, x ) ) - evaluateSolutionAt( basis.affine_mapping_1D( [-1, 1], elem_domain, x ), coeff, uspline_bext ) )
    abs_error = quadrature.quad( abs_err_fun, elem_domain, num_qp )
    return abs_error

def computeFitError( target_fun, coeff, uspline_bext ):
    num_elems = bext.getNumElems( uspline_bext )
    abs_error = 0.0
    for elem_idx in range( 0, num_elems ):
        elem_id = bext.elemIdFromElemIdx( uspline_bext, elem_idx )
        abs_error += computeElementFitError( target_fun, coeff, uspline_bext, elem_id )
    domain = bext.getDomain( uspline_bext )
    target_fun_norm, _ = scipy.integrate.quad( lambda x: abs( target_fun(x) ), domain[0], domain[1], epsrel = 1e-12, limit = num_elems * 100 )
    rel_error = abs_error / target_fun_norm
    return abs_error, rel_error

## CLI ARGUMENT PARSING
def prepareCommandInputs( target_fun_str, domain, degree, continuity ):
    spline_space = { "domain": domain, "degree": degree, "continuity": continuity }
    target_fun = sympy.parsing.sympy_parser.parse_expr( target_fun_str )
    target_fun = sympy.lambdify( sympy.symbols( "x", real = True ), target_fun )
    return target_fun, spline_space

def parseCommandLineArguments( ):
    parser = argparse.ArgumentParser()
    parser.add_argument( "–function", "-f",   nargs = 1,   type = str,   required = True )
    parser.add_argument( "–domain", "-d",     nargs = 2,   type = float, required = True )
    parser.add_argument( "–degree", "-p",     nargs = ’+’, type = int,   required = True )
    parser.add_argument( "–continuity", "-c", nargs = ’+’, type = int,   required = True )
    args = parser.parse_args( )
    return args.function[0], args.domain, args.degree, args.continuity

## TEST CALLING FROM PYTHON
class Test_python( unittest.TestCase ):
    def test_run( self ):
        target_fun_str = "sin(pi*x)"
        domain = [ 0, 1 ]
        degree = [ 2, 2 ]
        continuity = [ -1, 1, -1 ]
        target_fun, spline_space = prepareCommandInputs( target_fun_str, domain, degree, continuity )
        sol = main( target_fun, spline_space )

## EXAMPLE USAGE FROM CLI
if __name__ == "__main__":
    target_fun_str, domain, degree, continuity = parseCommandLineArguments( )
    target_fun, spline_space = prepareCommandInputs( target_fun_str, domain, degree, continuity )
    main( target_fun, spline_space )