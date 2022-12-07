import numpy
import unittest
import basis
import mesh_adapted
import quadrature
import readBEXT_JSON
    
# problem = { "elastic_modulus": 100,
#                      "area": 0.01,
#                      "length": 1.0,
#                      "traction": { "value": 1e-3, "position": 1.0 },
#                      "displacement": { "value": 0.0, "position": 0.0 },
#                      "body_force": 0.0 }

# uspline_bext = readBEXT_JSON.readBEXT( "test_usplineF123.json" )

def assembleForceVector(problem,uspline_bext):
    target_fun = lambda x: problem["body_force"]
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
    F_total = applyTractions(problem,uspline_bext) + F
    return F_total

def applyTractions(problem,uspline_bext):
    solution_basis = basis.evalSplineBasis1D
    num_nodes = readBEXT_JSON.getNumNodes(uspline_bext)
    e = readBEXT_JSON.getElementIdContainingPoint( uspline_bext, problem["traction"]["position"] )
    p = readBEXT_JSON.getElementDegree( uspline_bext, e)
    omega = readBEXT_JSON.getElementDomain( uspline_bext, e)
    C = readBEXT_JSON.getElementExtractionOperator( uspline_bext, e)
    Ft = numpy.zeros(num_nodes)
    for a in range(0,p+1):
        A = readBEXT_JSON.getElementNodeIds(uspline_bext, e)[a]
        Ft[A] = solution_basis(problem["traction"]["position"],C,a,omega) * problem["traction"]["value"]
    return Ft
        
# val_t = applyTractions(problem,uspline_bext)
# val = assembleForceVector(problem,uspline_bext)