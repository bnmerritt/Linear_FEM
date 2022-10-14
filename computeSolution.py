import numpy 
import mesh
import mesh_adapted

# def computeSolution(target_fun,domain,num_elems,degree):
#     xmin = domain[0]
#     xmax = domain[1]
#     node_coords = mesh.generateMesh(xmin,xmax,num_elems,degree)[0]
#     ien_array = mesh.generateMesh(xmin,xmax,num_elems,degree)[1]
#     test_solution = []
#     for node in range(0,len(node_coords)):
#         solution = target_fun(node_coords[node])
#         test_solution.append(solution)
#     test_solution = numpy.asarray(test_solution)
#     return test_solution, node_coords, ien_array

def computeSolution(target_fun,domain,num_elems,degree):
    xmin = domain[0]
    xmax = domain[1]
    node_coords = mesh_adapted.generateMesh(xmin,xmax,degree)[0]
    ien_array = mesh_adapted.generateMesh(xmin,xmax,degree)[1]
    test_solution = []
    for node in range(0,len(node_coords)):
        solution = target_fun(node_coords[node])
        test_solution.append(solution)
    test_solution = numpy.asarray(test_solution)
    return test_solution, node_coords, ien_array