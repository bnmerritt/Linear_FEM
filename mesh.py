import numpy

# def generateMesh(xmin,xmax,num_elems,degree):
#     ien_array = []
#     if degree == 1:
#         node_coords = numpy.linspace(xmin,xmax,num_elems+1)
#         for i in range(0,len(node_coords)-1):
#             ien = [i,i+1]
#             ien_array.append(ien)
#         ien_array = numpy.asarray(ien_array)
#     elif degree == 2:
#         node_coords = numpy.linspace(xmin,xmax,(2*num_elems)+1)
#         for i in range(0,len(node_coords)-1,2):
#             ien = [i,i+1,i+2]
#             ien_array.append(ien)
#         ien_array = numpy.asarray(ien_array)
#     return node_coords, ien_array

def generateMesh(xmin,xmax,num_elems,degree):
    ien_array = []
    node_coords = numpy.linspace(xmin,xmax,(degree*num_elems)+1)
    for i in range(0,len(node_coords)-1,degree):
        ien = numpy.linspace(i,i+degree,degree+1)
        ien = ien.astype(int)
        ien = list(ien)
        ien_array.append(ien)
    ien_array = numpy.asarray(ien_array)
    return node_coords, ien_array

