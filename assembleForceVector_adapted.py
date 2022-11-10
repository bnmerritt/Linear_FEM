import numpy
import unittest
import basis
import mesh_adapted
import quadrature

class Test_assembleForceVector( unittest.TestCase ):
    def test_lagrange_const_force_fun( self ):
        domain = [ 0, 1 ]
        degree = [ 3, 3 ]
        target_fun = lambda x: numpy.pi
        node_coords, ien_array = mesh_adapted.generateMesh( domain[0], domain[1], degree )
        test_force_vector = assembleForceVector( target_fun = target_fun, node_coords = node_coords, ien_array = ien_array, solution_basis = basis.evalLagrangeBasis1D )
        gold_force_vector = numpy.array( [ numpy.pi / 16.0, 3.0 * numpy.pi / 16.0, 3.0 * numpy.pi / 16.0, numpy.pi / 8.0, 3.0 * numpy.pi / 16.0, 3.0 * numpy.pi / 16.0, numpy.pi / 16.0 ] )
        self.assertTrue( numpy.allclose( test_force_vector, gold_force_vector ) )
    
    def test_lagrange_linear_force_fun( self ):
        domain = [ 0, 1 ]
        degree = [ 3, 3 ]
        target_fun = lambda x: 2*x + numpy.pi
        node_coords, ien_array = mesh_adapted.generateMesh( domain[0], domain[1], degree )
        test_force_vector = assembleForceVector( target_fun = target_fun, node_coords = node_coords, ien_array = ien_array, solution_basis = basis.evalLagrangeBasis1D )
        gold_force_vector = numpy.array( [ 0.20468287, 0.62654862, 0.73904862, 0.51769908, 0.81404862, 0.92654862, 0.31301621 ] )
        self.assertTrue( numpy.allclose( test_force_vector, gold_force_vector ) )
    
    def test_lagrange_quadratic_force_fun( self ):
        domain = [ 0, 1 ]
        degree = [ 3, 3 ]
        target_fun = lambda x: x**2.0
        node_coords, ien_array = mesh_adapted.generateMesh( domain[0], domain[1], degree )
        test_force_vector = assembleForceVector( target_fun = target_fun, node_coords = node_coords, ien_array = ien_array, solution_basis = basis.evalLagrangeBasis1D )
        gold_force_vector = numpy.array( [ 1.04166667e-03, 0, 2.81250000e-02, 3.33333333e-02, 6.56250000e-02, 1.50000000e-01, 5.52083333e-02 ] )
        self.assertTrue( numpy.allclose( test_force_vector, gold_force_vector ) )

    def test_bernstein_const_force_fun( self ):
        domain = [ 0, 1 ]
        degree = [ 3, 3 ]
        target_fun = lambda x: numpy.pi
        node_coords, ien_array = mesh_adapted.generateMesh( domain[0], domain[1], degree )
        test_force_vector = assembleForceVector( target_fun = target_fun, node_coords = node_coords, ien_array = ien_array, solution_basis = basis.evalBernsteinBasis1D )
        gold_force_vector = numpy.array( [ numpy.pi / 8.0, numpy.pi / 8.0, numpy.pi / 8.0, numpy.pi / 4.0, numpy.pi / 8.0, numpy.pi / 8.0, numpy.pi / 8.0 ] )
        self.assertTrue( numpy.allclose( test_force_vector, gold_force_vector ) )
    
    def test_bernstein_linear_force_fun( self ):
        domain = [ 0, 1 ]
        degree = [ 3, 3 ]
        target_fun = lambda x: 2*x + numpy.pi
        node_coords, ien_array = mesh_adapted.generateMesh( domain[0], domain[1], degree )
        test_force_vector = assembleForceVector( target_fun = target_fun, node_coords = node_coords, ien_array = ien_array, solution_basis = basis.evalBernsteinBasis1D )
        gold_force_vector = numpy.array( [ 0.41769908, 0.44269908, 0.46769908, 1.03539816, 0.56769908, 0.59269908, 0.61769908 ] )
        self.assertTrue( numpy.allclose( test_force_vector, gold_force_vector ) )
    
    def test_bernstein_quadratic_force_fun( self ):
        domain = [ 0, 1 ]
        degree = [ 3, 3 ]
        target_fun = lambda x: x**2.0
        node_coords, ien_array = mesh_adapted.generateMesh( domain[0], domain[1], degree )
        test_force_vector = assembleForceVector( target_fun = target_fun, node_coords = node_coords, ien_array = ien_array, solution_basis = basis.evalBernsteinBasis1D )
        gold_force_vector = numpy.array( [ 1/480, 1/160, 1/80, 1/15, 1/16, 13/160, 49/480 ] )
        self.assertTrue( numpy.allclose( test_force_vector, gold_force_vector ) )
        
        
def assembleForceVector(target_fun,node_coords,ien_array,solution_basis):
    num_elems = len(ien_array)
    F = numpy.zeros(len(node_coords))
    for elem in range(0,num_elems):
        elem_nodes = len(ien_array[elem])
        elem_degree = elem_nodes - 1
        node_idx = ien_array[elem][0]
        num_basis_vec = elem_degree + 1
        elem_domain = [node_coords[ien_array[elem][0]],node_coords[ien_array[elem][-1]]]
        qp, w = quadrature.computeGaussLegendreQuadrature(elem_nodes)
        qp_domain = [-1,1]
        qp_fun = ((elem_domain[-1] - elem_domain[0])/2)*qp + (elem_domain[0] + elem_domain[-1])/2
        derivative = (elem_domain[-1] - elem_domain[0]) / 2
        for A in range(0,num_basis_vec):
            for k in range(0,len(qp)):
                F[node_idx] += solution_basis(qp[k],elem_degree,A,qp_domain) * target_fun(qp_fun[k]) * w[k] * derivative
            node_idx = node_idx +1
    return F

unittest.main()


# domain = [ 0, 1 ]
# degree = [ 3, 3 ]
# target_fun = lambda x: x**2.0
# node_coords, ien_array = mesh_adapted.generateMesh( domain[0], domain[1], degree )
# solution_basis = basis.evalBernsteinBasis1D
# num_elems = len(ien_array)
# F = numpy.zeros(len(node_coords))
# for elem in range(0,num_elems):
#     elem_nodes = len(ien_array[elem])
#     elem_degree = elem_nodes - 1
#     node_idx = ien_array[elem][0]
#     num_basis_vec = elem_degree + 1
#     elem_domain = [node_coords[ien_array[elem][0]],node_coords[ien_array[elem][-1]]]
#     qp, w = quadrature.computeGaussLegendreQuadrature(elem_nodes)
#     qp_fun = ((elem_domain[-1] - elem_domain[0])/2)*qp + (elem_domain[0] + elem_domain[-1])/2
#     derivative = (elem_domain[-1] - elem_domain[0]) / 2
#     for A in range(0,num_basis_vec):
#         for k in range(0,len(qp)):
#             F[node_idx] += solution_basis(qp[k],elem_degree,A,elem_domain) * target_fun(qp_fun[k]) * w[k] * derivative
#         node_idx = node_idx + 1