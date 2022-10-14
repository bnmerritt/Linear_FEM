####################################################################################### Question 57
import numpy
        
def generateMesh(xmin,xmax,degree):
    num_elems = len(degree)
    elem_length = (xmax - xmin) / num_elems
    element = 0
    node_coords = []
    ien_array = []
    key = []
    num_nodes_elem = degree[element] + 1
    key.append(element)
    xmin_elem = xmin + element*elem_length
    xmax_elem = xmin_elem + elem_length
    node_coords_elem = numpy.linspace(xmin_elem,xmax_elem,num_nodes_elem)
    node_coords_elem = list(node_coords_elem)
    node_coords.extend(node_coords_elem)
    ien = numpy.linspace(0,degree[0],degree[0]+1)
    ien = ien.astype(int)
    ien = list(ien)
    ien_array.append(ien)
    for element in range(1,num_elems):
        num_nodes_elem = degree[element] + 1
        key.append(element)
        xmin_elem = xmin + element*elem_length
        xmax_elem = xmin_elem + elem_length
        node_coords_elem = numpy.linspace(xmin_elem,xmax_elem,num_nodes_elem)
        node_coords_elem = list(node_coords_elem)
        node_coords.extend(node_coords_elem)
        ien = numpy.linspace(ien_array[-1][-1],ien_array[-1][-1] + (num_nodes_elem-1),num_nodes_elem)
        ien = ien.astype(int)
        ien = list(ien)
        ien_array.append(ien)     
    ien_array = dict(zip(key, ien_array))
    node_coords = numpy.unique(numpy.asarray(node_coords))
    return node_coords, ien_array
