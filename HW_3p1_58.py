####################################################################################### Question 58
import numpy
import mesh_adapted
import basis
import matplotlib.pyplot as plt


domain = [-1,1]
xmin = domain[0]
xmax = domain[1]
domain2 = [0,4]
xmin2 = domain2[0]
xmax2 = domain2[1]
x_all = numpy.linspace(xmin,xmax,100)
degree = [ 1,2,3,4] 
node_coords, ien_array = mesh_adapted.generateMesh(xmin2,xmax2,degree)
num_elems = len(degree)
eval_basis = basis.evaluateBernsteinBasis1D
xelem = numpy.zeros((len(x_all),num_elems))
Bernstein_basis = numpy.zeros((len(x_all),sum(degree)+len(degree)))
column = 0
element = 0
for element in range(0,num_elems):
    xelem[:,element] = numpy.linspace(node_coords[ien_array[element][0]],node_coords[ien_array[element][-1]],100)
    deg = degree[element]
    for basis_idx in range(0,len(ien_array[element])):
        for xi in range(0,len(x_all)):
            Bernstein_basis[xi,column] = eval_basis(x_all[xi], deg, basis_idx)
        plt.plot(xelem[:,element],Bernstein_basis[:,column],'-')
        column = column + 1
plt.xlim(0,4)
plt.ylim(0,1)
plt.savefig('Q58.png', dpi=600)
