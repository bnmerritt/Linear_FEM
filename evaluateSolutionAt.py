import numpy
import mesh
import basis

# def evaluateSolutionAt(x,coeff,node_coords,ien_array,eval_basis):
#     for element in range(0,len(ien_array)):
#         if x >= node_coords[ien_array[element,0]] and x <= node_coords[ien_array[element,-1]]:
#             elem_idx = element
#             break
            
#     elem_nodes = ien_array[elem_idx]
#     elem_domain = [node_coords[elem_nodes[0]] , node_coords[elem_nodes[-1]]]
#     param_coord = 2*((x - elem_domain[0]) / (elem_domain[-1] - elem_domain[0])) - 1 #now exists between [-1,1] for Lagrange input
#     sol_at_point = 0
#     for n in range(0,len(elem_nodes)):
#         curr_node = elem_nodes[n]
#         sol_at_point += coeff[curr_node]*basis.evalLagrangeBasis1D(param_coord,len(elem_nodes)-1,n)
#     return sol_at_point

def evaluateSolutionAt(x,coeff,node_coords,ien_array,eval_basis):
    for element in range(0,len(ien_array)):
        if x >= node_coords[ien_array[element][0]] and x <= node_coords[ien_array[element][-1]]:
            elem_idx = element
            break
            
    elem_nodes = ien_array[elem_idx]
    elem_domain = [node_coords[elem_nodes[0]] , node_coords[elem_nodes[-1]]]
    param_coord = 2*((x - elem_domain[0]) / (elem_domain[-1] - elem_domain[0])) - 1 #now exists between [-1,1] for Lagrange input
    sol_at_point = 0
    for n in range(0,len(elem_nodes)):
        curr_node = elem_nodes[n]
        sol_at_point += coeff[curr_node]*basis.evalLagrangeBasis1D(param_coord,len(elem_nodes)-1,n)
    return sol_at_point