####################################################################################### Question 58
import numpy
import math
import mesh_adapted
import computeSolution
import evaluateSolutionAt
import basis
import matplotlib.pyplot as plt


domain = [-1,1]
xmin = domain[0]
xmax = domain[1]
x_all = numpy.linspace(xmin,xmax,100)
degree = [ 1,2,3,4] 
num_elems = len(degree)
eval_basis = basis.evaluateBernsteinBasis1D
# node_coords = []
# for element in range(0,len(degree)):
#     node_coords, ien_array = mesh_adapted.generateMesh(xmin,xmax,degree[element])

xelem1 = numpy.linspace(0,1,100)
xelem2 = numpy.linspace(1,2,100)
xelem3 = numpy.linspace(2,3,100)
xelem4 = numpy.linspace(3,4,100)

Bernstein1_basis0 = numpy.zeros(len(x_all))
Bernstein1_basis1L = numpy.zeros(len(x_all))
Bernstein2_basis1R = numpy.zeros(len(x_all))
Bernstein2_basis2 = numpy.zeros(len(x_all))
Bernstein2_basis3L = numpy.zeros(len(x_all))
Bernstein3_basis3R = numpy.zeros(len(x_all))
Bernstein3_basis4 = numpy.zeros(len(x_all))
Bernstein3_basis5 = numpy.zeros(len(x_all))
Bernstein3_basis6L = numpy.zeros(len(x_all))
Bernstein4_basis6R = numpy.zeros(len(x_all))
Bernstein4_basis7 = numpy.zeros(len(x_all))
Bernstein4_basis8 = numpy.zeros(len(x_all))
Bernstein4_basis9 = numpy.zeros(len(x_all))
Bernstein4_basis10 = numpy.zeros(len(x_all))

for xi in range(0,len(x_all)):
    Bernstein1_basis0[xi] = eval_basis(x_all[xi], 1, 0)
    Bernstein1_basis1L[xi] = eval_basis(x_all[xi], 1, 1)
    Bernstein2_basis1R[xi] = eval_basis(x_all[xi], 2, 0)
    Bernstein2_basis2[xi] = eval_basis(x_all[xi], 2, 1)
    Bernstein2_basis3L[xi] = eval_basis(x_all[xi], 2, 2)
    Bernstein3_basis3R[xi] = eval_basis(x_all[xi], 3, 0)
    Bernstein3_basis4[xi] = eval_basis(x_all[xi], 3, 1)
    Bernstein3_basis5[xi] = eval_basis(x_all[xi], 3, 2)
    Bernstein3_basis6L[xi] = eval_basis(x_all[xi], 3, 3)
    Bernstein4_basis6R[xi] = eval_basis(x_all[xi], 4, 0)
    Bernstein4_basis7[xi] = eval_basis(x_all[xi], 4, 1)
    Bernstein4_basis8[xi] = eval_basis(x_all[xi], 4, 2)
    Bernstein4_basis9[xi] = eval_basis(x_all[xi], 4, 3)
    Bernstein4_basis10[xi] = eval_basis(x_all[xi], 4, 4)
    
plt.plot(xelem1,Bernstein1_basis0,'-')
plt.plot(xelem1,Bernstein1_basis1L,'-')
plt.plot(xelem2,Bernstein2_basis1R,'-')
plt.plot(xelem2,Bernstein2_basis2,'-')
plt.plot(xelem2,Bernstein2_basis3L,'-')
plt.plot(xelem3,Bernstein3_basis3R,'-')
plt.plot(xelem3,Bernstein3_basis4,'-')
plt.plot(xelem3,Bernstein3_basis5,'-')
plt.plot(xelem3,Bernstein3_basis6L,'-')
plt.plot(xelem4,Bernstein4_basis6R,'-')
plt.plot(xelem4,Bernstein4_basis7,'-')
plt.plot(xelem4,Bernstein4_basis8,'-')
plt.plot(xelem4,Bernstein4_basis9,'-')
plt.plot(xelem4,Bernstein4_basis10,'-')
plt.xlim(0,4)
plt.ylim(0,1)
plt.xlabel('x')
plt.savefig('Q58.png', dpi=600)



