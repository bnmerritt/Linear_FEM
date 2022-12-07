import numpy
import unittest
import basis
import mesh_adapted
import quadrature
import readBEXT_JSON

class Test_assembleForceVector( unittest.TestCase ):
    def test_const_force_fun_two_element_linear_bspline( self ):
        target_fun = lambda x: numpy.pi
        # spline_space = { "domain": [-1, 1], "degree": [ 1, 1 ], "continuity": [ -1, 0, -1 ] }
        # uspline.make_uspline_mesh( spline_space, "temp_uspline" )
        # uspline_bext = bext.readBEXT( "temp_uspline.json" )
        uspline_bext = readBEXT_JSON.readBEXT("test_usplineF123.json")
        test_force_vector = assembleForceVector( target_fun = target_fun, uspline_bext = uspline_bext )
        gold_force_vector = numpy.array( [ numpy.pi / 2.0, numpy.pi, numpy.pi / 2.0 ] )
        self.assertTrue( numpy.allclose( test_force_vector, gold_force_vector ) )

    def test_linear_force_fun_two_element_linear_bspline( self ):
        target_fun = lambda x: x
        # spline_space = { "domain": [-1, 1], "degree": [ 1, 1 ], "continuity": [ -1, 0, -1 ] }
        # uspline.make_uspline_mesh( spline_space, "temp_uspline" )
        # uspline_bext = bext.readBEXT( "temp_uspline.json" )
        uspline_bext = readBEXT_JSON.readBEXT("test_usplineF123.json")
        test_force_vector = assembleForceVector( target_fun = target_fun, uspline_bext = uspline_bext )
        gold_force_vector = numpy.array( [ -1.0/3.0, 0.0, 1.0/3.0 ] )
        self.assertTrue( numpy.allclose( test_force_vector, gold_force_vector ) )

    def test_quadratic_force_fun_two_element_linear_bspline( self ):
        target_fun = lambda x: x**2
        # spline_space = { "domain": [-1, 1], "degree": [ 1, 1 ], "continuity": [ -1, 0, -1 ] }
        # uspline.make_uspline_mesh( spline_space, "temp_uspline" )
        # uspline_bext = bext.readBEXT( "temp_uspline.json" )
        uspline_bext = readBEXT_JSON.readBEXT("test_usplineF123.json")
        test_force_vector = assembleForceVector( target_fun = target_fun, uspline_bext = uspline_bext )
        gold_force_vector = numpy.array( [ 1.0/4.0, 1.0/6.0, 1.0/4.0 ] )
        self.assertTrue( numpy.allclose( test_force_vector, gold_force_vector ) )

    def test_const_force_fun_two_element_quadratic_bspline( self ):
        target_fun = lambda x: numpy.pi
        # spline_space = { "domain": [-1, 1], "degree": [ 2, 2 ], "continuity": [ -1, 1, -1 ] }
        # uspline.make_uspline_mesh( spline_space, "temp_uspline" )
        # uspline_bext = bext.readBEXT( "temp_uspline.json" )
        uspline_bext = readBEXT_JSON.readBEXT("test_usplineF456.json")
        test_force_vector = assembleForceVector( target_fun = target_fun, uspline_bext = uspline_bext )
        gold_force_vector = numpy.array( [ numpy.pi/3.0, 2.0*numpy.pi/3.0, 2.0*numpy.pi/3.0, numpy.pi/3.0 ] )
        self.assertTrue( numpy.allclose( test_force_vector, gold_force_vector ) )

    def test_linear_force_fun_two_element_quadratic_bspline( self ):
        target_fun = lambda x: x
        # spline_space = { "domain": [-1, 1], "degree": [ 2, 2 ], "continuity": [ -1, 1, -1 ] }
        # uspline.make_uspline_mesh( spline_space, "temp_uspline" )
        # uspline_bext = bext.readBEXT( "temp_uspline.json" )
        uspline_bext = readBEXT_JSON.readBEXT("test_usplineF456.json")
        test_force_vector = assembleForceVector( target_fun = target_fun, uspline_bext = uspline_bext )
        gold_force_vector = numpy.array( [ -1.0/4.0, -1.0/6.0, 1.0/6.0, 1.0/4.0 ] )
        self.assertTrue( numpy.allclose( test_force_vector, gold_force_vector ) )

    def test_quadratic_force_fun_two_element_quadratic_bspline( self ):
        target_fun = lambda x: x**2
        # spline_space = { "domain": [-1, 1], "degree": [ 2, 2 ], "continuity": [ -1, 1, -1 ] }
        # uspline.make_uspline_mesh( spline_space, "temp_uspline" )
        # uspline_bext = bext.readBEXT( "temp_uspline.json" )
        uspline_bext = readBEXT_JSON.readBEXT("test_usplineF456.json")
        test_force_vector = assembleForceVector( target_fun = target_fun, uspline_bext = uspline_bext )
        gold_force_vector = numpy.array( [ 2.0/10.0, 2.0/15.0, 2.0/15.0, 2.0/10.0 ] )
        self.assertTrue( numpy.allclose( test_force_vector, gold_force_vector ) )
        
        
def assembleForceVector(target_fun,uspline_bext):
    solution_basis = basis.evalSplineBasis1D
    num_elems = readBEXT_JSON.getNumElems(uspline_bext) 
    num_nodes = readBEXT_JSON.getNumNodes(uspline_bext)
    F = numpy.zeros(num_nodes)
    for elem in range(0,num_elems):
        elem_id = readBEXT_JSON.elemIdFromElemIdx(uspline_bext,elem)
        elem_nodes = readBEXT_JSON.getElementNodeIds(uspline_bext, elem_id)
        elem_degree = readBEXT_JSON.getElementDegree(uspline_bext, elem_id) 
        node_idx = readBEXT_JSON.getElementNodeIds(uspline_bext, elem_id)
        num_basis_vec = elem_degree + 1
        elem_domain = readBEXT_JSON.getElementDomain(uspline_bext, elem_id)
        elem_extraction_operator = readBEXT_JSON.getElementExtractionOperator(uspline_bext, elem_id)
        qp, w = quadrature.computeGaussLegendreQuadrature(len(elem_nodes)) 
        qp_domain = [-1,1]
        qp_fun = ((elem_domain[-1] - elem_domain[0])/2)*qp + (elem_domain[0] + elem_domain[-1])/2
        derivative = (elem_domain[-1] - elem_domain[0]) / 2
        for a in range(0,num_basis_vec):
            A = node_idx[a]
            for k in range(0,len(qp)):
                F[A] += solution_basis(qp[k],elem_extraction_operator,a,qp_domain) * target_fun(qp_fun[k]) * w[k] * derivative
    return F

unittest.main()

# target_fun = lambda x: numpy.pi
# uspline_bext = readBEXT_JSON.readBEXT("test_usplineF123.json")  
# solution_basis = basis.evalSplineBasis1D
# num_elems = readBEXT_JSON.getNumElems(uspline_bext) 
# num_nodes = readBEXT_JSON.getNumNodes(uspline_bext)
# F = numpy.zeros(num_nodes)
# for elem in range(0,num_elems):
#     elem_id = readBEXT_JSON.elemIdFromElemIdx(uspline_bext,elem)
#     elem_nodes = readBEXT_JSON.getElementNodeIds(uspline_bext, elem_id)
#     elem_degree = readBEXT_JSON.getElementDegree(uspline_bext, elem_id) 
#     node_idx = readBEXT_JSON.getElementNodeIds(uspline_bext, elem_id)
#     num_basis_vec = elem_degree + 1
#     elem_domain = readBEXT_JSON.getElementDomain(uspline_bext, elem_id)
#     elem_extraction_operator = readBEXT_JSON.getElementExtractionOperator(uspline_bext, elem_id)
#     qp, w = quadrature.computeGaussLegendreQuadrature(len(elem_nodes))
#     qp_domain = [-1,1]
#     qp_fun = ((elem_domain[-1] - elem_domain[0])/2)*qp + (elem_domain[0] + elem_domain[-1])/2
#     derivative = (elem_domain[-1] - elem_domain[0]) / 2
#     for a in range(0,num_basis_vec):
#         A = node_idx[a]
#         for k in range(0,len(qp)):
#             # F[node_idx] += solution_basis(qp[k],elem_degree,A,qp_domain) * target_fun(qp_fun[k]) * w[k] * derivative
#             F[A] += solution_basis(qp[k],elem_extraction_operator,a,qp_domain) * target_fun(qp_fun[k]) * w[k] * derivative