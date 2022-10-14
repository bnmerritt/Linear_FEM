####################################################################################### Question 56
import numpy
import math
import mesh_adapted
import computeSolution
import evaluateSolutionAt
import basis
import matplotlib.pyplot as plt

target_fun1 = lambda x : math.sin(math.pi*x)
target_fun2 = lambda x : math.exp(x)
target_fun3 = lambda x : math.erfc(x)

domain = [-1,1]
xmin = domain[0]
xmax = domain[1]
degree = [ 2, 2, 2]  #quadratic Lagrange
num_elems = len(degree)
eval_basis = basis.evalLagrangeBasis1D

coeff1, node_coords, ien_array = computeSolution.computeSolution(target_fun1,domain,num_elems,degree)
coeff2, node_coords, ien_array = computeSolution.computeSolution(target_fun2,domain,num_elems,degree)
coeff3, node_coords, ien_array = computeSolution.computeSolution(target_fun3,domain,num_elems,degree)

basis_fun1 = numpy.zeros(len(node_coords))
basis_fun2 = numpy.zeros(len(node_coords))
basis_fun3 = numpy.zeros(len(node_coords))
sol_at_point1 = numpy.zeros(len(node_coords))
sol_at_point2 = numpy.zeros(len(node_coords))
sol_at_point3 = numpy.zeros(len(node_coords))
for element in range(0,num_elems):
    ien_subarray = ien_array[element] #basis indeces for Lagrange input
    basis_idx = 0
    for node in range(ien_subarray[0],ien_subarray[-1]):
        x = node_coords[node]
        sol_at_point1[node] = evaluateSolutionAt.evaluateSolutionAt(x,coeff1,node_coords,ien_array,eval_basis)
        sol_at_point2[node] = evaluateSolutionAt.evaluateSolutionAt(x,coeff2,node_coords,ien_array,eval_basis)
        sol_at_point3[node] = evaluateSolutionAt.evaluateSolutionAt(x,coeff3,node_coords,ien_array,eval_basis)
        basis_idx = basis_idx + 1
        
x_all = numpy.linspace(-1,1,100)
ysin = numpy.sin(numpy.pi*x_all)
yexpx = numpy.zeros(len(x_all))
for xi in range(0,len(x_all)):
    yexpx[xi] = math.exp(x_all[xi])
yerfcx = numpy.zeros(len(x_all))
for xi in range(0,len(x_all)):
    yerfcx[xi] = math.erfc(x_all[xi])

xelem1 = numpy.linspace(node_coords[0],node_coords[2],100)
xelem2 = numpy.linspace(node_coords[2],node_coords[4],100)
xelem3 = numpy.linspace(node_coords[4],node_coords[6],100)
Lagrange1_basis0 = numpy.zeros(len(x_all))
Lagrange1_basis1 = numpy.zeros(len(x_all))
Lagrange1_basis2L = numpy.zeros(len(x_all))
Lagrange1_basis2R = numpy.zeros(len(x_all))
Lagrange1_basis3 = numpy.zeros(len(x_all))
Lagrange1_basis4L = numpy.zeros(len(x_all))
Lagrange1_basis4R = numpy.zeros(len(x_all))
Lagrange1_basis5 = numpy.zeros(len(x_all))
Lagrange1_basis6 = numpy.zeros(len(x_all))
Lagrange2_basis0 = numpy.zeros(len(x_all))
Lagrange2_basis1 = numpy.zeros(len(x_all))
Lagrange2_basis2L = numpy.zeros(len(x_all))
Lagrange2_basis2R = numpy.zeros(len(x_all))
Lagrange2_basis3 = numpy.zeros(len(x_all))
Lagrange2_basis4L = numpy.zeros(len(x_all))
Lagrange2_basis4R = numpy.zeros(len(x_all))
Lagrange2_basis5 = numpy.zeros(len(x_all))
Lagrange2_basis6 = numpy.zeros(len(x_all))
Lagrange3_basis0 = numpy.zeros(len(x_all))
Lagrange3_basis1 = numpy.zeros(len(x_all))
Lagrange3_basis2L = numpy.zeros(len(x_all))
Lagrange3_basis2R = numpy.zeros(len(x_all))
Lagrange3_basis3 = numpy.zeros(len(x_all))
Lagrange3_basis4L = numpy.zeros(len(x_all))
Lagrange3_basis4R = numpy.zeros(len(x_all))
Lagrange3_basis5 = numpy.zeros(len(x_all))
Lagrange3_basis6 = numpy.zeros(len(x_all))
for xi in range(0,len(x_all)):
    Lagrange1_basis0[xi] = coeff1[0] * basis.evalLagrangeBasis1D(x_all[xi], 2, 0)
    Lagrange1_basis1[xi] = coeff1[1] * basis.evalLagrangeBasis1D(x_all[xi], 2, 1)
    Lagrange1_basis2L[xi] = coeff1[2] * basis.evalLagrangeBasis1D(x_all[xi], 2, 2)
    Lagrange1_basis2R[xi] = coeff1[2] * basis.evalLagrangeBasis1D(x_all[xi], 2, 0)
    Lagrange1_basis3[xi] = coeff1[3] * basis.evalLagrangeBasis1D(x_all[xi], 2, 1)
    Lagrange1_basis4L[xi] = coeff1[4] * basis.evalLagrangeBasis1D(x_all[xi], 2, 2)
    Lagrange1_basis4R[xi] = coeff1[4] * basis.evalLagrangeBasis1D(x_all[xi], 2, 0)
    Lagrange1_basis5[xi] = coeff1[5] * basis.evalLagrangeBasis1D(x_all[xi], 2, 1)
    Lagrange1_basis6[xi] = coeff1[6] * basis.evalLagrangeBasis1D(x_all[xi], 2, 2)
    Lagrange2_basis0[xi] = coeff2[0] * basis.evalLagrangeBasis1D(x_all[xi], 2, 0)
    Lagrange2_basis1[xi] = coeff2[1] * basis.evalLagrangeBasis1D(x_all[xi], 2, 1)
    Lagrange2_basis2L[xi] = coeff2[2] * basis.evalLagrangeBasis1D(x_all[xi], 2, 2)
    Lagrange2_basis2R[xi] = coeff2[2] * basis.evalLagrangeBasis1D(x_all[xi], 2, 0)
    Lagrange2_basis3[xi] = coeff2[3] * basis.evalLagrangeBasis1D(x_all[xi], 2, 1)
    Lagrange2_basis4L[xi] = coeff2[4] * basis.evalLagrangeBasis1D(x_all[xi], 2, 2)
    Lagrange2_basis4R[xi] = coeff2[4] * basis.evalLagrangeBasis1D(x_all[xi], 2, 0)
    Lagrange2_basis5[xi] = coeff2[5] * basis.evalLagrangeBasis1D(x_all[xi], 2, 1)
    Lagrange2_basis6[xi] = coeff2[6] * basis.evalLagrangeBasis1D(x_all[xi], 2, 2)
    Lagrange3_basis0[xi] = coeff3[0] * basis.evalLagrangeBasis1D(x_all[xi], 2, 0)
    Lagrange3_basis1[xi] = coeff3[1] * basis.evalLagrangeBasis1D(x_all[xi], 2, 1)
    Lagrange3_basis2L[xi] = coeff3[2] * basis.evalLagrangeBasis1D(x_all[xi], 2, 2)
    Lagrange3_basis2R[xi] = coeff3[2] * basis.evalLagrangeBasis1D(x_all[xi], 2, 0)
    Lagrange3_basis3[xi] = coeff3[3] * basis.evalLagrangeBasis1D(x_all[xi], 2, 1)
    Lagrange3_basis4L[xi] = coeff3[4] * basis.evalLagrangeBasis1D(x_all[xi], 2, 2)
    Lagrange3_basis4R[xi] = coeff3[4] * basis.evalLagrangeBasis1D(x_all[xi], 2, 0)
    Lagrange3_basis5[xi] = coeff3[5] * basis.evalLagrangeBasis1D(x_all[xi], 2, 1)
    Lagrange3_basis6[xi] = coeff3[6] * basis.evalLagrangeBasis1D(x_all[xi], 2, 2)
    
plt.figure(0)
plt.plot(x_all, ysin, '-k')
plt.plot(node_coords,sol_at_point1,'ko')
plt.plot(node_coords,coeff1,'bx')
plt.plot(xelem1,Lagrange1_basis0,'-')
plt.plot(xelem1,Lagrange1_basis1,'-')
plt.plot(xelem1,Lagrange1_basis2L,'-')
plt.plot(xelem2,Lagrange1_basis2R,'-')
plt.plot(xelem2,Lagrange1_basis3,'-')
plt.plot(xelem2,Lagrange1_basis4L,'-')
plt.plot(xelem3,Lagrange1_basis4R,'-')
plt.plot(xelem3,Lagrange1_basis5,'-')
plt.plot(xelem3,Lagrange1_basis6,'-')
for i in range(0, len(node_coords)+1, 2):
    plt.axvline(node_coords[i], color = 'k', linestyle = '--')
plt.xlim(domain)
plt.ylim(-1,1)
plt.xlabel('x')
plt.ylabel('target_fun')
plt.title('sin(pi*x)')
plt.savefig('Q56_1.png', dpi=600)

plt.figure(1)
plt.plot(x_all, yexpx, '-k')
plt.plot(node_coords,sol_at_point2,'ko')
plt.plot(node_coords,coeff2,'bx')
plt.plot(xelem1,Lagrange2_basis0,'-')
plt.plot(xelem1,Lagrange2_basis1,'-')
plt.plot(xelem1,Lagrange2_basis2L,'-')
plt.plot(xelem2,Lagrange2_basis2R,'-')
plt.plot(xelem2,Lagrange2_basis3,'-')
plt.plot(xelem2,Lagrange2_basis4L,'-')
plt.plot(xelem3,Lagrange2_basis4R,'-')
plt.plot(xelem3,Lagrange2_basis5,'-')
plt.plot(xelem3,Lagrange2_basis6,'-')
for i in range(0, len(node_coords)+1, 2):
    plt.axvline(node_coords[i], color = 'k', linestyle = '--')
plt.xlim(domain)
plt.ylim(-0.5,3)
plt.xlabel('x')
plt.ylabel('target_fun')
plt.title('exp(x)')
plt.savefig('Q56_2.png', dpi=600)

plt.figure(2)
plt.plot(x_all, yerfcx, '-k')
plt.plot(node_coords,sol_at_point3,'ko')
plt.plot(node_coords,coeff3,'bx')
plt.plot(xelem1,Lagrange3_basis0,'-')
plt.plot(xelem1,Lagrange3_basis1,'-')
plt.plot(xelem1,Lagrange3_basis2L,'-')
plt.plot(xelem2,Lagrange3_basis2R,'-')
plt.plot(xelem2,Lagrange3_basis3,'-')
plt.plot(xelem2,Lagrange3_basis4L,'-')
plt.plot(xelem3,Lagrange3_basis4R,'-')
plt.plot(xelem3,Lagrange3_basis5,'-')
plt.plot(xelem3,Lagrange3_basis6,'-')
for i in range(0, len(node_coords)+1, 2):
    plt.axvline(node_coords[i], color = 'k', linestyle = '--')
plt.xlim(domain)
plt.ylim(-0.5,2)
plt.xlabel('x')
plt.ylabel('target_fun')
plt.title('erfc(x)')
plt.savefig('Q56_3.png', dpi=600)

