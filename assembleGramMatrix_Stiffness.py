import numpy
import mesh_adapted
import quadrature
import readBEXT_JSON
# import uspline
import unittest
import basis

class test_assembleStressMatrix( unittest.TestCase ):
       def test_one_linear_C0_element( self ):
              problem = { "elastic_modulus": 100,
                     "area": 0.01,
                     "length": 1.0,
                     "traction": { "value": 1e-3, "position": 1.0 },
                     "displacement": { "value": 0.0, "position": 0.0 },
                     "body_force": 0.0 }
              # spline_space = { "domain": [0, problem[ "length" ]], "degree": [ 1 ], "continuity": [ -1, -1 ] }
              # uspline.make_uspline_mesh( spline_space, "temp_uspline" )
              # uspline_bext = readBEXT_JSON.readBEXT( "temp_uspline.json" )
              uspline_bext = readBEXT_JSON.readBEXT( "uspline_SM1.json" )
              test_stiffness_matrix = assembleStiffnessMatrix( problem = problem, uspline_bext = uspline_bext )
              gold_stiffness_matrix = numpy.array( [ [ 1.0, -1.0 ], [ -1.0, 1.0 ] ] )
              self.assertTrue( numpy.allclose( test_stiffness_matrix, gold_stiffness_matrix ) )

       def test_two_linear_C0_element( self ):
              problem = { "elastic_modulus": 100,
                     "area": 0.01,
                     "length": 1.0,
                     "traction": { "value": 1e-3, "position": 1.0 },
                     "displacement": { "value": 0.0, "position": 0.0 },
                     "body_force": 0.0 }
              # spline_space = { "domain": [0, problem[ "length" ]], "degree": [ 1, 1 ], "continuity": [ -1, 0, -1 ] }
              # uspline.make_uspline_mesh( spline_space, "temp_uspline" )
              # uspline_bext = readBEXT_JSON.readBEXT( "temp_uspline.json" )
              uspline_bext = readBEXT_JSON.readBEXT( "uspline_SM2.json" )
              test_stiffness_matrix = assembleStiffnessMatrix( problem = problem, uspline_bext = uspline_bext )
              gold_stiffness_matrix = numpy.array( [ [ 2.0, -2.0, 0.0 ], [ -2.0, 4.0, -2.0 ], [ 0.0, -2.0, 2.0 ] ] )
              self.assertTrue( numpy.allclose( test_stiffness_matrix, gold_stiffness_matrix ) )

       def test_one_quadratic_C0_element( self ):
              problem = { "elastic_modulus": 100,
                     "area": 0.01,
                     "length": 1.0,
                     "traction": { "value": 1e-3, "position": 1.0 },
                     "displacement": { "value": 0.0, "position": 0.0 },
                     "body_force": 0.0 }
              # spline_space = { "domain": [0, problem[ "length" ]], "degree": [ 2 ], "continuity": [ -1, -1 ] }
              # uspline.make_uspline_mesh( spline_space, "temp_uspline" )
              # uspline_bext = readBEXT_JSON.readBEXT( "temp_uspline.json" )
              uspline_bext = readBEXT_JSON.readBEXT( "uspline_SM3.json" )
              test_stiffness_matrix = assembleStiffnessMatrix( problem = problem, uspline_bext = uspline_bext )
              gold_stiffness_matrix = numpy.array( [ [  4.0 / 3.0, -2.0 / 3.0, -2.0 / 3.0 ],
                                                 [ -2.0 / 3.0,  4.0 / 3.0, -2.0 / 3.0 ],
                                                 [ -2.0 / 3.0, -2.0 / 3.0,  4.0 / 3.0 ] ] )
              self.assertTrue( numpy.allclose( test_stiffness_matrix, gold_stiffness_matrix ) )

       def test_two_quadratic_C0_element( self ):
              problem = { "elastic_modulus": 100,
                     "area": 0.01,
                     "length": 1.0,
                     "traction": { "value": 1e-3, "position": 1.0 },
                     "displacement": { "value": 0.0, "position": 0.0 },
                     "body_force": 0.0 }
              # spline_space = { "domain": [0, problem[ "length" ]], "degree": [ 2, 2 ], "continuity": [ -1, 0, -1 ] }
              # uspline.make_uspline_mesh( spline_space, "temp_uspline" )
              # uspline_bext = readBEXT_JSON.readBEXT( "temp_uspline.json" )
              uspline_bext = readBEXT_JSON.readBEXT( "uspline_SM4.json" )
              test_stiffness_matrix = assembleStiffnessMatrix( problem = problem, uspline_bext = uspline_bext )
              gold_stiffness_matrix = numpy.array( [ [  8.0 / 3.0, -4.0 / 3.0, -4.0 / 3.0,  0.0,        0.0 ],
                                                 [ -4.0 / 3.0,  8.0 / 3.0, -4.0 / 3.0,  0.0,        0.0 ],
                                                 [ -4.0 / 3.0, -4.0 / 3.0, 16.0 / 3.0, -4.0 / 3.0, -4.0 / 3.0 ],
                                                 [  0.0,        0.0,       -4.0 / 3.0,  8.0 / 3.0, -4.0 / 3.0 ],
                                                 [  0.0,        0.0,       -4.0 / 3.0, -4.0 / 3.0,  8.0 / 3.0 ] ] )
              self.assertTrue( numpy.allclose( test_stiffness_matrix, gold_stiffness_matrix ) )

       def test_two_quadratic_C1_element( self ):
              problem = { "elastic_modulus": 100,
                     "area": 0.01,
                     "length": 1.0,
                     "traction": { "value": 1e-3, "position": 1.0 },
                     "displacement": { "value": 0.0, "position": 0.0 },
                     "body_force": 0.0 }
              # spline_space = { "domain": [0, problem[ "length" ]], "degree": [ 2, 2 ], "continuity": [ -1, 1, -1 ] }
              # uspline.make_uspline_mesh( spline_space, "temp_uspline" )
              # uspline_bext = readBEXT_JSON.readBEXT( "temp_uspline.json" )
              uspline_bext = readBEXT_JSON.readBEXT( "uspline_SM5.json" )
              test_stiffness_matrix = assembleStiffnessMatrix( problem = problem, uspline_bext = uspline_bext )
              gold_stiffness_matrix = numpy.array( [ [  8.0 / 3.0, -2.0,       -2.0/ 3.0,   0.0 ],
                                                 [ -2.0,        8.0 / 3.0,  0.0,       -2.0 / 3.0 ],
                                                 [ -2.0 / 3.0,  0.0,        8.0 / 3.0, -2.0 ],
                                                 [  0.0,       -2.0 / 3.0, -2.0,        8.0 / 3.0 ] ] )
              self.assertTrue( numpy.allclose( test_stiffness_matrix, gold_stiffness_matrix ) )
             
def assembleStiffnessMatrix(problem, uspline_bext): #need to redefine everything from uspline_bext input
    solution_basis = basis.evalSplineBasisDeriv1D
    num_elems = readBEXT_JSON.getNumElems(uspline_bext) 
    num_nodes = readBEXT_JSON.getNumNodes(uspline_bext)
    M = numpy.zeros((num_nodes,num_nodes))
    for elem in range(0,num_elems): #loop through elements
        elem_id = readBEXT_JSON.elemIdFromElemIdx(uspline_bext,elem)
        elem_nodes = readBEXT_JSON.getElementNodeIds(uspline_bext, elem_id) 
        elem_degree = readBEXT_JSON.getElementDegree(uspline_bext, elem_id) 
        node_idx = readBEXT_JSON.getElementNodeIds(uspline_bext, elem_id) 
        elem_domain = readBEXT_JSON.getElementDomain(uspline_bext, elem_id)
        elem_extraction_operator = readBEXT_JSON.getElementExtractionOperator(uspline_bext, elem_id)
        qp, w = quadrature.computeGaussLegendreQuadrature(len(elem_nodes))
        qp_domain = [-1,1]
        num_basis_vec = elem_degree + 1
        derivative = 2/(elem_domain[-1] - elem_domain[0])
        deriv = 1
        for a in range(0,num_basis_vec): #basis_index
            A = node_idx[a]
            for b in range(0,num_basis_vec): #basis_index
                B = node_idx[b]
                for k in range(0,len(qp)):
                    M[A,B] +=  solution_basis(elem_extraction_operator,a,deriv,qp_domain,qp[k]) * solution_basis(elem_extraction_operator,b,deriv,qp_domain,qp[k]) * w[k] * derivative * problem["area"] * problem["elastic_modulus"]
    return M

unittest.main()


# problem = { "elastic_modulus": 100,
#                      "area": 0.01,
#                      "length": 1.0,
#                      "traction": { "value": 1e-3, "position": 1.0 },
#                      "displacement": { "value": 0.0, "position": 0.0 },
#                      "body_force": 0.0 }
# uspline_bext = readBEXT_JSON.readBEXT( "uspline_SM2.json" )
# solution_basis = basis.evalSplineBasisDeriv1D
# num_elems = readBEXT_JSON.getNumElems(uspline_bext) 
# num_nodes = readBEXT_JSON.getNumNodes(uspline_bext)
# M = numpy.zeros((num_nodes,num_nodes))
# for elem in range(0,num_elems): #loop through elements
#     elem_id = readBEXT_JSON.elemIdFromElemIdx(uspline_bext,elem)
#     elem_nodes = readBEXT_JSON.getElementNodeIds(uspline_bext, elem_id) 
#     elem_degree = readBEXT_JSON.getElementDegree(uspline_bext, elem_id) 
#     node_idx = readBEXT_JSON.getElementNodeIds(uspline_bext, elem_id) 
#     elem_domain = readBEXT_JSON.getElementDomain(uspline_bext, elem_id)
#     elem_extraction_operator = readBEXT_JSON.getElementExtractionOperator(uspline_bext, elem_id)
#     qp, w = quadrature.computeGaussLegendreQuadrature(len(elem_nodes))
#     qp_domain = [-1,1]
#     num_basis_vec = elem_degree + 1
#     derivative = 2/(elem_domain[-1] - elem_domain[0])
#     deriv = 1
#     for a in range(0,num_basis_vec): #basis_index
#         A = node_idx[a]
#         for b in range(0,num_basis_vec): #basis_index
#             B = node_idx[b]
#             for k in range(0,len(qp)):
#                 M[A,B] +=  solution_basis(elem_extraction_operator,a,deriv,qp_domain,qp[k]) * solution_basis(elem_extraction_operator,b,deriv,qp_domain,qp[k]) * w[k] * derivative * problem["area"] * problem["elastic_modulus"]