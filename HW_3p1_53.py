####################################################################################### Question 53
import numpy
import unittest

class Test_computeSolution( unittest.TestCase ):
    def test_single_linear_element_poly( self ):
        test_solution, node_coords, ien_array = computeSolution( target_fun = lambda x : x, domain = [-1.0, 1.0 ], num_elems = 1, degree = 1 )
        gold_solution = numpy.array( [ -1.0, 1.0 ] )
        self.assertTrue( numpy.allclose( test_solution, gold_solution ) )
    
    def test_single_quad_element_poly( self ):
        test_solution, node_coords, ien_array = computeSolution( target_fun = lambda x : x**2, domain = [-1.0, 1.0 ], num_elems = 1, degree = 2 )
        gold_solution = numpy.array( [ 1.0, 0.0, 1.0 ] )
        self.assertTrue( numpy.allclose( test_solution, gold_solution ) )
    
    def test_two_linear_element_poly( self ):
        test_solution, node_coords, ien_array = computeSolution( target_fun = lambda x : x**2, domain = [-1.0, 1.0 ], num_elems = 2, degree = 1 )
        gold_solution = numpy.array( [ 1.0, 0.0, 1.0 ] )
        self.assertTrue( numpy.allclose( test_solution, gold_solution ) )
    
    def test_four_quad_element_poly( self ):
        test_solution, node_coords, ien_array = computeSolution( target_fun = lambda x : x**2, domain = [-1.0, 1.0 ], num_elems = 4, degree = 1 )
        gold_solution = numpy.array( [ 1.0, 0.25, 0.0, 0.25, 1.0 ] )
        self.assertTrue( numpy.allclose( test_solution, gold_solution ) )
        
def generateMesh1D(xmin,xmax,num_elems,degree):
    ien_array = []
    if degree == 1:
        node_coords = numpy.linspace(xmin,xmax,num_elems+1)
        for i in range(0,len(node_coords)-1):
            ien = [i,i+1]
            ien_array.append(ien)
        ien_array = numpy.asarray(ien_array)
    elif degree == 2:
        node_coords = numpy.linspace(xmin,xmax,(2*num_elems)+1)
        for i in range(0,len(node_coords)-1,2):
            ien = [i,i+1,i+2]
            ien_array.append(ien)
        ien_array = numpy.asarray(ien_array)
    return node_coords, ien_array
    
def computeSolution(target_fun,domain,num_elems,degree):
    xmin = domain[0]
    xmax = domain[1]
    node_coords = generateMesh1D(xmin,xmax,num_elems,degree)[0]
    ien_array = generateMesh1D(xmin,xmax,num_elems,degree)[1]
    test_solution = []
    for node in range(0,len(node_coords)):
        solution = target_fun(node_coords[node])
        test_solution.append(solution)
    test_solution = numpy.asarray(test_solution)
    return test_solution, node_coords, ien_array
            
unittest.main()