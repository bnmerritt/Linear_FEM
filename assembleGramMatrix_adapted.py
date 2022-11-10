import numpy
import unittest
import basis
import mesh_adapted
import quadrature

class Test_assembleGramMatrix( unittest.TestCase ):
    def test_linear_lagrange( self ):
        domain = [ 0, 1 ]
        degree = [ 1, 1 ]
        node_coords, ien_array = mesh_adapted.generateMesh( domain[0], domain[1], degree )
        test_gram_matrix = assembleGramMatrix( node_coords = node_coords, ien_array = ien_array, solution_basis = basis.evalLagrangeBasis1D )
        gold_gram_matrix = numpy.array( [ [1/6, 1/12, 0], [1/12, 1/3, 1/12], [0, 1/12, 1/6] ] )
        self.assertTrue( numpy.allclose( test_gram_matrix, gold_gram_matrix ) )

    def test_quadratic_lagrange( self ):
        domain = [ 0, 1 ]
        degree = [ 2, 2 ]
        node_coords, ien_array = mesh_adapted.generateMesh( domain[0], domain[1], degree )
        test_gram_matrix = assembleGramMatrix( node_coords = node_coords, ien_array = ien_array, solution_basis = basis.evalLagrangeBasis1D )
        gold_gram_matrix = numpy.array( [ [1/15, 1/30, -1/60, 0, 0 ], [1/30, 4/15, 1/30, 0, 0], [-1/60, 1/30, 2/15, 1/30, -1/60], [ 0, 0, 1/30, 4/15, 1/30], [0, 0, -1/60, 1/30, 1/15] ] )
        self.assertTrue( numpy.allclose( test_gram_matrix, gold_gram_matrix ) )
    
    def test_cubic_lagrange( self ):
        domain = [ 0, 1 ]
        degree = [ 3, 3 ]
        node_coords, ien_array = mesh_adapted.generateMesh( domain[0], domain[1], degree )
        test_gram_matrix = assembleGramMatrix( node_coords = node_coords, ien_array = ien_array, solution_basis = basis.evalLagrangeBasis1D )
        gold_gram_matrix = numpy.array( [ [ 0.03809524,  0.02946429, -0.01071429,  0.00565476,  0.00000000,  0.00000000,  0.00000000 ], 
                                          [ 0.02946429,  0.19285714, -0.02410714, -0.01071429,  0.00000000,  0.00000000,  0.00000000 ], 
                                          [-0.01071429, -0.02410714,  0.19285714,  0.02946429,  0.00000000,  0.00000000,  0.00000000 ], 
                                          [ 0.00565476, -0.01071429,  0.02946429,  0.07619048,  0.02946429, -0.01071429,  0.00565476 ], 
                                          [ 0.00000000,  0.00000000,  0.00000000,  0.02946429,  0.19285714, -0.02410714, -0.01071429 ], 
                                          [ 0.00000000,  0.00000000,  0.00000000, -0.01071429, -0.02410714,  0.19285714,  0.02946429 ], 
                                          [ 0.00000000,  0.00000000,  0.00000000,  0.00565476, -0.01071429,  0.02946429,  0.03809524 ] ] )
        self.assertTrue( numpy.allclose( test_gram_matrix, gold_gram_matrix ) )
    
    def test_linear_bernstein( self ):
        domain = [ 0, 1 ]
        degree = [ 1, 1 ]
        node_coords, ien_array = mesh_adapted.generateMesh( domain[0], domain[1], degree )
        test_gram_matrix = assembleGramMatrix( node_coords = node_coords, ien_array = ien_array, solution_basis = basis.evalBernsteinBasis1D )
        gold_gram_matrix = numpy.array( [ [1/6, 1/12, 0], [1/12, 1/3, 1/12], [0, 1/12, 1/6] ] )
        self.assertTrue( numpy.allclose( test_gram_matrix, gold_gram_matrix ) )

    def test_quadratic_bernstein( self ):
        domain = [ 0, 1 ]
        degree = [ 2, 2 ]
        node_coords, ien_array = mesh_adapted.generateMesh( domain[0], domain[1], degree )
        test_gram_matrix = assembleGramMatrix( node_coords = node_coords, ien_array = ien_array, solution_basis = basis.evalBernsteinBasis1D )
        gold_gram_matrix = numpy.array( [ [1/10, 1/20, 1/60, 0, 0 ], [1/20, 1/15, 1/20, 0, 0 ], [1/60, 1/20, 1/5, 1/20, 1/60], [0, 0, 1/20, 1/15, 1/20], [0, 0, 1/60, 1/20, 1/10] ] )
        self.assertTrue( numpy.allclose( test_gram_matrix, gold_gram_matrix ) )
    
    def test_cubic_bernstein( self ):
        domain = [ 0, 1 ]
        degree = [ 3, 3 ]
        node_coords, ien_array = mesh_adapted.generateMesh( domain[0], domain[1], degree )
        test_gram_matrix = assembleGramMatrix( node_coords = node_coords, ien_array = ien_array, solution_basis = basis.evalBernsteinBasis1D )
        gold_gram_matrix = numpy.array( [ [1/14, 1/28, 1/70, 1/280, 0, 0, 0 ], [1/28, 3/70, 9/280, 1/70, 0, 0, 0 ], [1/70, 9/280, 3/70, 1/28, 0, 0, 0 ], [1/280, 1/70, 1/28, 1/7, 1/28, 1/70, 1/280], [0, 0, 0, 1/28, 3/70, 9/280, 1/70], [0, 0, 0, 1/70, 9/280, 3/70, 1/28], [0, 0, 0, 1/280, 1/70, 1/28, 1/14 ] ] )
        self.assertTrue( numpy.allclose( test_gram_matrix, gold_gram_matrix ) )

def assembleGramMatrix(node_coords, ien_array,solution_basis):
    num_elems = len(ien_array)
    M = numpy.zeros((len(node_coords),len(node_coords)))
    for elem in range(0,num_elems):
        elem_nodes = len(ien_array[elem])
        elem_degree = elem_nodes - 1
        node_idx = ien_array[elem][0]
        elem_domain = [node_coords[ien_array[elem][0]],node_coords[ien_array[elem][-1]]]
        qp, w = quadrature.computeGaussLegendreQuadrature(elem_nodes)
        qp_domain = [-1,1]
        num_basis_vec = elem_degree + 1
        derivative = (elem_domain[-1] - elem_domain[0]) / 2
        for A in range(0,num_basis_vec): #basis_index
            for B in range(0,num_basis_vec): #basis_index
                for k in range(0,len(qp)):
                    M[A+node_idx,B+node_idx] += solution_basis(qp[k],elem_degree,A,qp_domain) * solution_basis(qp[k],elem_degree,B,qp_domain) * w[k] * derivative   
    return M

unittest.main()

# domain = [ 0, 1 ]
# degree = [ 1, 1 ]
# node_coords, ien_array = mesh_adapted.generateMesh( domain[0], domain[1], degree )
# solution_basis = basis.evalLagrangeBasis1D
# num_elems = len(ien_array)
# M = numpy.zeros((len(node_coords),len(node_coords)))
# for elem in range(0,num_elems):
#     elem_nodes = len(ien_array[elem])
#     elem_degree = elem_nodes - 1
#     node_idx = ien_array[elem][0]
#     elem_domain = [node_coords[ien_array[elem][0]],node_coords[ien_array[elem][-1]]]
#     qp, w = quadrature.computeGaussLegendreQuadrature(elem_nodes)
#     num_basis_vec = elem_degree + 1
#     derivative = (elem_domain[-1] - elem_domain[0]) / 2
#     for A in range(0,num_basis_vec): #basis_index
#         for B in range(0,num_basis_vec): #basis_index
#             for k in range(0,len(qp)):
#                 M[A+node_idx,B+node_idx] += solution_basis(qp[k],elem_degree,A,elem_domain) * solution_basis(qp[k],elem_degree,B,elem_domain) * w[k] * derivative 