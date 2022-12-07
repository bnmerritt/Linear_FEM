import numpy
import mesh_adapted
import quadrature
import readBEXT_JSON
# import uspline
import unittest
import basis

class Test_assembleGramMatrix( unittest.TestCase ):
    def test_two_element_linear_bspline( self ):
        target_fun = lambda x: x**0
        # spline_space = { "domain": [0, 2], "degree": [ 1, 1 ], "continuity": [ -1, 0, -1 ] }
        # uspline.make_uspline_mesh( spline_space, "temp_uspline" )
        # uspline_bext = readBEXT_JSON.readBEXT( "temp_uspline.json" )
        uspline_bext = readBEXT_JSON.readBEXT( "temp_uspline1.json" )
        test_gram_matrix = assembleGramMatrix( uspline_bext = uspline_bext )
        gold_gram_matrix = numpy.array( [ [ 1.0/3.0, 1.0/6.0, 0.0 ],
                                          [ 1.0/6.0, 2.0/3.0, 1.0/6.0 ],
                                          [ 0.0, 1.0/6.0, 1.0/3.0 ] ] )
        self.assertTrue( numpy.allclose( test_gram_matrix, gold_gram_matrix ) )

    def test_two_element_quadratic_bspline( self ):
        target_fun = lambda x: x**0
        # spline_space = { "domain": [0, 2], "degree": [ 2, 2 ], "continuity": [ -1, 1, -1 ] }
        # uspline.make_uspline_mesh( spline_space, "temp_uspline" )
        # uspline_bext = readBEXT_JSON.readBEXT( "temp_uspline.json" )
        uspline_bext = readBEXT_JSON.readBEXT( "temp_uspline2.json" )
        test_gram_matrix = assembleGramMatrix( uspline_bext = uspline_bext )
        gold_gram_matrix = numpy.array( [ [ 1.0/5.0, 7.0/60.0, 1.0/60.0, 0.0 ],
                                          [ 7.0/60.0, 1.0/3.0, 1.0/5.0, 1.0/60.0],
                                          [ 1.0/60.0, 1.0/5.0, 1.0/3.0, 7.0/60.0 ],
                                          [ 0.0, 1.0/60.0, 7.0/60.0, 1.0/5.0] ] )
        self.assertTrue( numpy.allclose( test_gram_matrix, gold_gram_matrix ) )

    def test_two_element_cubic_bspline( self ):
        # spline_space = { "domain": [0, 2], "degree": [ 3, 3 ], "continuity": [ -1, 2, -1 ] }
        # uspline.make_uspline_mesh( spline_space, "temp_uspline" )
        # uspline_bext = readBEXT_JSON.readBEXT( "temp_uspline.json" )
        uspline_bext = readBEXT_JSON.readBEXT( "temp_uspline3.json" )
        test_gram_matrix = assembleGramMatrix( uspline_bext = uspline_bext )
        gold_gram_matrix = numpy.array( [ [ 1.0/7.0, 7.0/80.0, 1.0/56.0, 1.0/560.0, 0.0 ],
                                          [ 7.0/80.0, 31.0/140.0, 39.0/280.0, 1.0/20.0, 1.0/560.0 ],
                                          [ 1.0/56.0, 39.0/280.0, 13.0/70.0, 39.0/280.0, 1.0/56.0 ],
                                          [ 1.0/560.0, 1.0/20.0, 39.0/280.0, 31.0/140.0, 7.0/80.0 ],
                                          [ 0.0, 1.0/560.0, 1.0/56.0, 7.0/80.0, 1.0/7.0 ] ] )
        self.assertTrue( numpy.allclose( test_gram_matrix, gold_gram_matrix ) )

# def assembleGramMatrix(node_coords, ien_array, solution_basis):
def assembleGramMatrix(uspline_bext): #need to redefine everything from uspline_bext input
    solution_basis = basis.evalSplineBasis1D
    num_elems = readBEXT_JSON.getNumElems(uspline_bext) 
    # num_elems = len(ien_array)
    # M = numpy.zeros((len(node_coords),len(node_coords)))
    num_nodes = readBEXT_JSON.getNumNodes(uspline_bext)
    M = numpy.zeros((num_nodes,num_nodes))
    for elem in range(0,num_elems): #loop through elements
        elem_id = readBEXT_JSON.elemIdFromElemIdx(uspline_bext,elem)
        elem_nodes = readBEXT_JSON.getElementNodeIds(uspline_bext, elem_id) 
        # elem_nodes = readBEXT_JSON.getVertexConnectivity(uspline_bext)[elem] 
        # elem_nodes = len(ien_array[elem])
        elem_degree = readBEXT_JSON.getElementDegree(uspline_bext, elem_id) 
        # elem_degree = elem_nodes - 1
        # M = numpy.zeros((elem_degree+1,elem_degree+1))
        node_idx = readBEXT_JSON.getElementNodeIds(uspline_bext, elem_id) 
        # node_idx = elem_nodes[0]
        # node_idx = readBEXT_JSON.getElementNodeIds(uspline_bext, elem)[0]
        # node_idx = ien_array[elem][0]
        elem_domain = readBEXT_JSON.getElementDomain(uspline_bext, elem_id)
        # elem_domain = [node_coords[ien_array[elem][0]],node_coords[ien_array[elem][-1]]]
        # uspline = readBEXT_JSON.readBEXT( "quadratic_bspline_match.json" )
        elem_extraction_operator = readBEXT_JSON.getElementExtractionOperator(uspline_bext, elem_id)
        qp, w = quadrature.computeGaussLegendreQuadrature(len(elem_nodes))
        # qp, w = quadrature.computeGaussLegendreQuadrature(elem_nodes)
        qp_domain = [-1,1]
        num_basis_vec = elem_degree + 1
        derivative = (elem_domain[-1] - elem_domain[0]) / 2
        for a in range(0,num_basis_vec): #basis_index
            A = node_idx[a]
            for b in range(0,num_basis_vec): #basis_index
                B = node_idx[b]
                for k in range(0,len(qp)):
                    # M[A,B] +=  solution_basis(qp[k],elem_degree,A,qp_domain) * solution_basis(qp[k],elem_degree,B,qp_domain) * w[k] * derivative
                    M[A,B] +=  solution_basis(qp[k],elem_extraction_operator,a,qp_domain) * solution_basis(qp[k],elem_extraction_operator,b,qp_domain) * w[k] * derivative
    return M

unittest.main()

uspline_bext = readBEXT_JSON.readBEXT( "temp_uspline1.json" )
solution_basis = basis.evalSplineBasis1D
num_elems = readBEXT_JSON.getNumElems(uspline_bext) 
# num_elems = len(ien_array)
# M = numpy.zeros((len(node_coords),len(node_coords)))
num_nodes = readBEXT_JSON.getNumNodes(uspline_bext)
M = numpy.zeros((num_nodes,num_nodes))
for elem in range(0,num_elems): #loop through elements
    elem_id = readBEXT_JSON.elemIdFromElemIdx(uspline_bext,elem)
    elem_nodes = readBEXT_JSON.getElementNodeIds(uspline_bext, elem_id) 
    # elem_nodes = readBEXT_JSON.getVertexConnectivity(uspline_bext)[elem] 
    # elem_nodes = len(ien_array[elem])
    elem_degree = readBEXT_JSON.getElementDegree(uspline_bext, elem_id) 
    # elem_degree = elem_nodes - 1
    # M = numpy.zeros((elem_degree+1,elem_degree+1))
    node_idx = readBEXT_JSON.getElementNodeIds(uspline_bext, elem_id) 
    # node_idx = elem_nodes[0]
    # node_idx = readBEXT_JSON.getElementNodeIds(uspline_bext, elem)[0]
    # node_idx = ien_array[elem][0]
    elem_domain = readBEXT_JSON.getElementDomain(uspline_bext, elem_id)
    # elem_domain = [node_coords[ien_array[elem][0]],node_coords[ien_array[elem][-1]]]
    # uspline = readBEXT_JSON.readBEXT( "quadratic_bspline_match.json" )
    elem_extraction_operator = readBEXT_JSON.getElementExtractionOperator(uspline_bext, elem_id)
    qp, w = quadrature.computeGaussLegendreQuadrature(len(elem_nodes))
    # qp, w = quadrature.computeGaussLegendreQuadrature(elem_nodes)
    qp_domain = [-1,1]
    num_basis_vec = elem_degree + 1
    derivative = (elem_domain[-1] - elem_domain[0]) / 2
    for a in range(0,num_basis_vec): #basis_index
        A = node_idx[a]
        for b in range(0,num_basis_vec): #basis_index
            B = node_idx[b]
            for k in range(0,len(qp)):
                # M[A,B] +=  solution_basis(qp[k],elem_degree,A,qp_domain) * solution_basis(qp[k],elem_degree,B,qp_domain) * w[k] * derivative
                M[A,B] +=  solution_basis(qp[k],elem_extraction_operator,a,qp_domain) * solution_basis(qp[k],elem_extraction_operator,b,qp_domain) * w[k] * derivative
    # M_corrected = numpy.dot(elem_extraction_operator,M)
    # for A in range(0,num_basis_vec): #basis_index
    #     for B in range(0,num_basis_vec): #basis_index 
    #         M_final[A+node_idx,B+node_idx] += M_corrected[A,B]