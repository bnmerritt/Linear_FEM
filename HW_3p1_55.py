####################################################################################### Question 55
import evaluateSolutionAt
import scipy
from scipy import integrate
import computeSolution
import matplotlib.pyplot as plt
import basis
import numpy
import math


def computeFitError( target_fun, coeff, node_coords, ien_array, eval_basis ):
    num_elems = ien_array.shape[0]
    domain = [ min( node_coords ), max( node_coords ) ]
    abs_err_fun = lambda x : abs( target_fun( x ) - evaluateSolutionAt.evaluateSolutionAt( x, coeff, node_coords, ien_array, eval_basis ) )
    fit_error, residual = scipy.integrate.quad( abs_err_fun, domain[0], domain[1], epsrel = 1e-12, limit = num_elems * 100 )
    return fit_error, residual

####################################################################################### h-refinement 1
target_fun = lambda x : x**2
domain = [0,1]
degree = 1
num_elems_vec = numpy.round(numpy.logspace(0,3,10))
eval_basis = basis.evalLagrangeBasis1D

fit_error = numpy.zeros(len(num_elems_vec))
residual = numpy.zeros(len(num_elems_vec))
num_nodes = []
for element in range(0,len(num_elems_vec)):
    num_elems = int(num_elems_vec[element])
    if degree == 1:
        nodes = num_elems + 1
        num_nodes.append(nodes)
    elif degree == 2:
        nodes = 2*num_elems + 1 #for quadratic
        num_nodes.append(nodes)
    coeff, node_coords, ien_array = computeSolution.computeSolution(target_fun,domain,num_elems,degree)
    fit_error[element],residual[element] = computeFitError( target_fun, coeff, node_coords, ien_array, eval_basis )
num_nodes = numpy.asarray(num_nodes)

plt.figure(0)
plt.plot(num_elems_vec,fit_error,'b*-')
plt.xlim([1, 1000])
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Number of Elements')
plt.ylabel('Fit Error')
plt.title('Linear Elements: x^2')
plt.savefig('Q55_1.png', dpi=600)

plt.figure(1)
plt.plot(num_nodes,fit_error,'b*-')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Number of Nodes')
plt.ylabel('Fit Error')
plt.title('Linear Elements: x^2')
plt.savefig('Q55_2.png', dpi=600)
####################################################################################### h-refinement 2
target_fun = lambda x : x**3
domain = [0,1]
degree = 2
num_elems_vec = numpy.round(numpy.logspace(0,3,10))
eval_basis = basis.evalLagrangeBasis1D

fit_error = numpy.zeros(len(num_elems_vec))
residual = numpy.zeros(len(num_elems_vec))
num_nodes = []
for element in range(0,len(num_elems_vec)):
    num_elems = int(num_elems_vec[element])
    if degree == 1:
        nodes = num_elems + 1
        num_nodes.append(nodes)
    elif degree == 2:
        nodes = 2*num_elems + 1 #for quadratic
        num_nodes.append(nodes)
    coeff, node_coords, ien_array = computeSolution.computeSolution(target_fun,domain,num_elems,degree)
    fit_error[element],residual[element] = computeFitError( target_fun, coeff, node_coords, ien_array, eval_basis )
num_nodes = numpy.asarray(num_nodes)

plt.figure(2)
plt.plot(num_elems_vec,fit_error,'b*-')
plt.xlim([1, 1000])
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Number of Elements')
plt.ylabel('Fit Error')
plt.title('Quadratic Elements: x^3')
plt.savefig('Q55_3.png', dpi=600)

plt.figure(3)
plt.plot(num_nodes,fit_error,'b*-')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Number of Nodes')
plt.ylabel('Fit Error')
plt.title('Quadratic Elements: x^3')
plt.savefig('Q55_4.png', dpi=600)
####################################################################################### h-refinement 3
target_fun = lambda x : math.sin(math.pi*x)
domain = [-1,1]
degree = 1
num_elems_vec = numpy.round(numpy.logspace(0,3,10))
eval_basis = basis.evalLagrangeBasis1D

fit_error = numpy.zeros(len(num_elems_vec))
residual = numpy.zeros(len(num_elems_vec))
num_nodes = []
for element in range(0,len(num_elems_vec)):
    num_elems = int(num_elems_vec[element])
    if degree == 1:
        nodes = num_elems + 1
        num_nodes.append(nodes)
    elif degree == 2:
        nodes = 2*num_elems + 1 #for quadratic
        num_nodes.append(nodes)
    coeff, node_coords, ien_array = computeSolution.computeSolution(target_fun,domain,num_elems,degree)
    fit_error[element],residual[element] = computeFitError( target_fun, coeff, node_coords, ien_array, eval_basis )
num_nodes = numpy.asarray(num_nodes)

plt.figure(4)
plt.plot(num_elems_vec,fit_error,'b*-')
plt.xlim([1, 1000])
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Number of Elements')
plt.ylabel('Fit Error')
plt.title('Linear Elements: sin(pi*x)')
plt.savefig('Q55_5.png', dpi=600)

plt.figure(5)
plt.plot(num_nodes,fit_error,'b*-')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Number of Nodes')
plt.ylabel('Fit Error')
plt.title('Linear Elements: sin(pi*x)')
plt.savefig('Q55_6.png', dpi=600)
####################################################################################### h-refinement 4
target_fun = lambda x : math.sin(math.pi*x)
domain = [-1,1]
degree = 2
num_elems_vec = numpy.round(numpy.logspace(0,3,10))
eval_basis = basis.evalLagrangeBasis1D

fit_error = numpy.zeros(len(num_elems_vec))
residual = numpy.zeros(len(num_elems_vec))
num_nodes = []
for element in range(0,len(num_elems_vec)):
    num_elems = int(num_elems_vec[element])
    if degree == 1:
        nodes = num_elems + 1
        num_nodes.append(nodes)
    elif degree == 2:
        nodes = 2*num_elems + 1 #for quadratic
        num_nodes.append(nodes)
    coeff, node_coords, ien_array = computeSolution.computeSolution(target_fun,domain,num_elems,degree)
    fit_error[element],residual[element] = computeFitError( target_fun, coeff, node_coords, ien_array, eval_basis )
num_nodes = numpy.asarray(num_nodes)

plt.figure(6)
plt.plot(num_elems_vec,fit_error,'b*-')
plt.xlim([1, 1000])
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Number of Elements')
plt.ylabel('Fit Error')
plt.title('Quadratic Elements: sin(pi*x)')
plt.savefig('Q55_7.png', dpi=600)

plt.figure(7)
plt.plot(num_nodes,fit_error,'b*-')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Number of Nodes')
plt.ylabel('Fit Error')
plt.title('Quadratic Elements: sin(pi*x)')
plt.savefig('Q55_8.png', dpi=600)
####################################################################################### p-refinement
target_fun = lambda x : math.sin(math.pi*x)
domain = [-1,1]
degree_vec = [1,2,3,4,5,6,7,8]
num_elems = 2
num_elems_vec = numpy.zeros(len(degree_vec))
for d in range(0,len(num_elems_vec)):
    num_elems_vec[d] = num_elems
eval_basis = basis.evalLagrangeBasis1D
fit_error = numpy.zeros(len(degree_vec))
residual = numpy.zeros(len(degree_vec))
num_nodes = numpy.zeros(len(degree_vec))
for deg in range(0,len(degree_vec)):
    degree = degree_vec[deg]
    nodes_per_element = degree + 1
    nodes_total = (num_elems*degree) + 1
    num_nodes[deg] = nodes_total
    coeff, node_coords, ien_array = computeSolution.computeSolution(target_fun,domain,num_elems,degree)
    fit_error[deg],residual[deg] = computeFitError( target_fun, coeff, node_coords, ien_array, eval_basis )

plt.figure(8)
plt.plot(num_elems_vec,fit_error,'bo-')
plt.xlim([1, 1000])
plt.ylim([1e-9,10])
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Number of Elements')
plt.ylabel('Fit Error')
plt.title('2 Elements, Degrees 1-8: sin(pi*x)')
plt.savefig('Q55_9.png', dpi=600)

plt.figure(9)
plt.plot(num_nodes,fit_error,'bo-')
plt.xlim([1, 1000])
plt.ylim([1e-9,10])
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Number of Nodes')
plt.ylabel('Fit Error')
plt.title('2 Elements, Degrees 1-8: sin(pi*x)')
plt.savefig('Q55_10.png', dpi=600)
####################################################################################### hp-refinement
target_fun = lambda x : math.sin(math.pi*x)
domain = [-1,1]
degree_vec = [1,2,3,4,5,6,7,8]
num_elems_vec = 2 ** numpy.array( range( 0, 7 ) )
eval_basis = basis.evalLagrangeBasis1D
fit_error = numpy.zeros((len(degree_vec),len(num_elems_vec)))
residual = numpy.zeros((len(degree_vec),len(num_elems_vec)))
num_nodes = numpy.zeros((len(degree_vec),len(num_elems_vec)))

for elem in range(0,len(num_elems_vec)):
    num_elems = num_elems_vec[elem]
    for deg in range(0,len(degree_vec)):
        degree = degree_vec[deg]
        nodes_per_element = degree + 1
        nodes_total = (num_elems*degree) + 1
        num_nodes[deg,elem] = nodes_total
        coeff, node_coords, ien_array = computeSolution.computeSolution(target_fun,domain,num_elems,degree)
        fit_error[deg,elem],residual[deg,elem] = computeFitError( target_fun, coeff, node_coords, ien_array, eval_basis )

plt.figure(10)
for elem in range(0,len(num_elems_vec)):
        for deg in range(0,len(degree_vec)):
            if deg == elem - 1:
                plt.plot(num_elems_vec[elem],fit_error[deg,elem],'bo-')
plt.xlim([1, 1000])
plt.ylim([1e-16,10])
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Number of Elements')
plt.ylabel('Fit Error')
plt.title('hp-refinement: sin(pi*x)')
plt.savefig('Q55_11.png', dpi=600)

plt.figure(11)
for elem in range(0,len(num_elems_vec)):
        for deg in range(0,len(degree_vec)):
            if deg == elem - 1:
                plt.plot(num_nodes[deg,elem],fit_error[deg,elem],'bo-')
plt.xlim([1, 1000])
plt.ylim([1e-16,10])
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Number of Nodes')
plt.ylabel('Fit Error')
plt.title('hp-refinement: sin(pi*x)')
plt.savefig('Q55_12.png', dpi=600)