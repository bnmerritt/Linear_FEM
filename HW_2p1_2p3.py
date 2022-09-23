# ME EN 507 HW 2.1-2.3
#9/19/22

import numpy as np
import sympy 
import math
import matplotlib.pyplot as plt

def taylorExpansion( fun, a, order ):
    x = list( fun.atoms( sympy.Symbol ) )[0]
    t = 0
    for i in range( 0, order + 1 ):
       df = sympy.diff( fun, x, i )
       term = ( df.subs( x, a ) / sympy.factorial( i ) ) * ( x - a )**i
       t += term
    return t
  
############################################################# 32.1 sin(pi*x)
xx = np.arange(-1,1.01,0.01)  
xxx = sympy.Symbol('xxx')  
plt.figure(0)
plt.plot(xx,np.sin(np.pi*xx), "-k", label = "sin(x)")
sin_Order0 = taylorExpansion(sympy.sin(np.pi*xxx),0,0)
sin_Orders0 = np.zeros(len(xx))
for j in range(0,len(xx)):
    sin_Orders0[j] = float(sin_Order0.subs(xxx,xx[j]))
plt.figure(0)
plt.plot(xx, sin_Orders0, "m", label = "Order 0")
sin_Orders = np.zeros((len(xx),4))
colors = ["b","r","g","c"]
for i in range(1,5):
    k = 2*i - 1
    sin_Order = taylorExpansion(sympy.sin(np.pi*xxx),0,k)
    for j in range(0,len(xx)):
        sin_Orders[j,i-1] = float(sin_Order.subs(xxx,xx[j]))
    plt.figure(0)
    plt.plot(xx, sin_Orders[:,i-1], colors[i-1], label = f"Order {str(k)}")
plt.legend(loc="upper left")
plt.xlim(-1, 1)
plt.ylim(-3, 3)
plt.xlabel('x')
plt.ylabel('f(x)')
plt.savefig('Q32_1.png', dpi=600)
############################################################# 32.2 e^x
xx = np.arange(-1,1.01,0.01)  
xxx = sympy.Symbol('xxx')  
plt.figure(1)
plt.plot(xx,np.exp(xx), "-k", label = "exp(x)")
ex_Orders = np.zeros((len(xx),5))
colors = ["m","b","r","g","c"]
for i in range(0,5):
    k = i
    ex_Order = taylorExpansion(sympy.exp(xxx),0,k)
    for j in range(0,len(xx)):
        ex_Orders[j,i] = float(ex_Order.subs(xxx,xx[j]))
    plt.figure(1)
    plt.plot(xx, ex_Orders[:,i], colors[i], label = f"Order {str(k)}")
plt.legend(loc="upper left")
plt.xlim(-1, 1)
plt.ylim(0, 3)
plt.xlabel('x')
plt.ylabel('f(x)')
plt.savefig('Q32_2.png', dpi=600)
############################################################# 32.3 erfc(x)
xx = np.arange(-2,2.01,0.01)  
xxx = sympy.Symbol('xxx')  
erfc = np.zeros(len(xx))
for xi in range(len(xx)):
    erfc[xi] = math.erfc(xx[xi])
plt.figure(2)
plt.plot(xx,erfc, "-k", label = "erfc(x)")
erfc_Order0 = taylorExpansion(sympy.erfc(xxx),0,0)
erfc_Orders0 = np.zeros(len(xx))
for j in range(0,len(xx)):
    erfc_Orders0[j] = float(erfc_Order0.subs(xxx,xx[j]))
plt.figure(2)
plt.plot(xx, erfc_Orders0, "m", label = "Order 0")
erfc_Orders = np.zeros((len(xx),6))
colors = ["b","r","g","c","y","k"]
for i in range(1,7):
    k = 2*i - 1
    erfc_Order = taylorExpansion(sympy.erfc(xxx),0,k)
    for j in range(0,len(xx)):
        erfc_Orders[j,i-1] = float(erfc_Order.subs(xxx,xx[j]))
    plt.figure(2)
    plt.plot(xx, erfc_Orders[:,i-1], colors[i-1], label = f"Order {str(k)}")
plt.legend(loc="lower left")
plt.xlim(-2, 2)
plt.ylim(-1, 3)
plt.xlabel('x')
plt.ylabel('f(x)')
plt.savefig('Q32_3.png', dpi=600)
############################################################# 32.1 sin error convergence
xx = np.arange(-1,1.01,0.01)  
sin = np.sin(np.pi*xx)
sin_symb = sympy.sin(np.pi*xxx)
sin_orders = np.arange(0,11,1)
sin_int = np.zeros(len(sin_orders))
for sin_i in range(0,len(sin_orders)):
    sin_Order = taylorExpansion(sympy.sin(np.pi*xxx),0,sin_orders[sin_i])
    sin_int[sin_i] = float(sympy.integrate(abs(sin_symb - sin_Order),(xxx,-1,1)))
plt.figure(3)
plt.plot(sin_orders,sin_int, "*k")
plt.xlim(0, 10)
plt.ylim(0.001, 10)
plt.yscale("log")
plt.xlabel('Order')
plt.ylabel('Error')
plt.savefig('Q32_4.png', dpi=600)
############################################################# 32.2 e^x error convergence
xx = np.arange(-1,1.01,0.01)  
ex = np.exp(xx)
ex_symb = sympy.exp(xxx)
ex_orders = np.arange(0,11,1)
ex_int = np.zeros(len(ex_orders))
for ex_i in range(0,len(ex_orders)):
    ex_Order = taylorExpansion(sympy.exp(xxx),0,ex_orders[ex_i])
    ex_int[ex_i] = float(sympy.integrate(abs(ex_symb - ex_Order),(xxx,-1,1)))
plt.figure(4)
plt.plot(ex_orders,ex_int, "*k")
plt.xlim(0, 10)
plt.ylim(10**-9, 10)
plt.yscale("log")
plt.xlabel('Order')
plt.ylabel('Error')
plt.savefig('Q32_5.png', dpi=600)
############################################################# 32.3 erfc error convergence
xx = np.arange(-2,2.01,0.01)  
erfc_symb = sympy.erfc(xxx)
erfc_orders = np.arange(0,11,1)
erfc_int = np.zeros(len(erfc_orders))
for erfc_i in range(0,len(erfc_orders)):
    erfc_Order = taylorExpansion(sympy.erfc(xxx),0,erfc_orders[erfc_i])
    erfc_int[erfc_i] = float(sympy.integrate(abs(erfc_symb - erfc_Order),(xxx,-2,2)))
plt.figure(5)
plt.plot(erfc_orders,erfc_int, "*k")
plt.xlim(0, 10)
plt.ylim(10**(-9), 10)
plt.yscale("log")
plt.xlabel('Order')
plt.ylabel('Error')
plt.savefig('Q32_6.png', dpi=600)
############################################################# 33 power (monomial) basis 1D w/ test
import unittest

class Test_evaluateMonomialBasis1D( unittest.TestCase ):
    def test_basisAtBounds( self ):
        self.assertAlmostEqual( first = evaluateMonomialBasis1D( degree = 0, variate = 0 ), second = 1.0, delta = 1e-12 )
        for p in range( 1, 11 ):
            self.assertAlmostEqual( first = evaluateMonomialBasis1D( degree = p, variate = 0 ), second = 0.0, delta = 1e-12 )
            self.assertAlmostEqual( first = evaluateMonomialBasis1D( degree = p, variate = 1 ), second = 1.0, delta = 1e-12 )

    def test_basisAtMidpoint( self ):
        for p in range( 0, 11 ):
            self.assertAlmostEqual( first = evaluateMonomialBasis1D( degree = p, variate = 0.5 ), second = 1 / ( 2**p ), delta = 1e-12 )

def evaluateMonomialBasis1D(variate,degree):
    return variate**degree

unittest.main()
############################################################# 34 plot power basis for 10 degrees
interval = np.arange(0,1.01,0.01)
degrees = [0,1,2,3,4,5,6,7,8,9,10]
answer34 = np.zeros((len(interval),len(degrees)))
for p in range(0,len(degrees)):
    for spot in range(0,len(interval)):
        answer34[spot,p] = evaluateMonomialBasis1D(interval[spot],degrees[p])
    plt.figure(6)    
    plt.plot(interval,answer34[:,p],'-k')
plt.xlabel('x')
plt.ylabel('Monomial Basis')
plt.savefig('Q34.png', dpi=600)
############################################################# 36 Legendre recurrence relation
import unittest
import math
import numpy

class Test_evalLegendreBasis1D( unittest.TestCase ):
    def test_basisAtBounds( self ):
        for p in range( 0, 2 ):
            if ( p % 2 == 0 ):
                self.assertAlmostEqual( first = evalLegendreBasis1D( degree = p, variate = -1 ), second = +1.0, delta = 1e-12 )
            else:
                self.assertAlmostEqual( first = evalLegendreBasis1D( degree = p, variate = -1 ), second = -1.0, delta = 1e-12 )
            self.assertAlmostEqual( first = evalLegendreBasis1D( degree = p, variate = +1 ), second = 1.0, delta = 1e-12 )

    def test_constant( self ):
        for x in numpy.linspace( -1, 1, 100 ):
            self.assertAlmostEqual( first = evalLegendreBasis1D( degree = 0, variate = x ), second = 1.0, delta = 1e-12 )

    def test_linear( self ):
        for x in numpy.linspace( -1, 1, 100 ):
            self.assertAlmostEqual( first = evalLegendreBasis1D( degree = 1, variate = x ), second = x, delta = 1e-12 )

    def test_quadratic_at_roots( self ):
        self.assertAlmostEqual( first = evalLegendreBasis1D( degree = 2, variate = -1.0 / math.sqrt(3.0) ), second = 0.0, delta = 1e-12 )
        self.assertAlmostEqual( first = evalLegendreBasis1D( degree = 2, variate = +1.0 / math.sqrt(3.0) ), second = 0.0, delta = 1e-12 )

    def test_cubic_at_roots( self ):
        self.assertAlmostEqual( first = evalLegendreBasis1D( degree = 3, variate = -math.sqrt( 3 / 5 ) ), second = 0.0, delta = 1e-12 )
        self.assertAlmostEqual( first = evalLegendreBasis1D( degree = 3, variate = 0 ), second = 0.0, delta = 1e-12 )
        self.assertAlmostEqual( first = evalLegendreBasis1D( degree = 3, variate = +math.sqrt( 3 / 5 ) ), second = 0.0, delta = 1e-12 )
        
def evalLegendreBasis1D(degree,variate):
    if degree == 0:
        val = 1
    elif degree == 1:
        val = variate
    else:
        i = degree - 1
        val = ((i + 1)**(-1)) * ((2*i+1)*variate*evalLegendreBasis1D(i,variate) - i*evalLegendreBasis1D(i-1,variate))
    return val
    
unittest.main()
############################################################# 38 Lagrange 
import unittest

class Test_evaluateLagrangeBasis1D( unittest.TestCase ):
    def test_linearLagrange( self ):
        self.assertAlmostEqual( first = evaluateLagrangeBasis1D( variate = -1, degree = 1, basis_idx = 0 ), second = 1.0, delta = 1e-12 )
        self.assertAlmostEqual( first = evaluateLagrangeBasis1D( variate = -1, degree = 1, basis_idx = 1 ), second = 0.0, delta = 1e-12 )
        self.assertAlmostEqual( first = evaluateLagrangeBasis1D( variate = +1, degree = 1, basis_idx = 0 ), second = 0.0, delta = 1e-12 )
        self.assertAlmostEqual( first = evaluateLagrangeBasis1D( variate = +1, degree = 1, basis_idx = 1 ), second = 1.0, delta = 1e-12 )

    def test_quadraticLagrange( self ):
        self.assertAlmostEqual( first = evaluateLagrangeBasis1D( variate = -1, degree = 2, basis_idx = 0 ), second = 1.0, delta = 1e-12 )
        self.assertAlmostEqual( first = evaluateLagrangeBasis1D( variate = -1, degree = 2, basis_idx = 1 ), second = 0.0, delta = 1e-12 )
        self.assertAlmostEqual( first = evaluateLagrangeBasis1D( variate = -1, degree = 2, basis_idx = 2 ), second = 0.0, delta = 1e-12 )
        self.assertAlmostEqual( first = evaluateLagrangeBasis1D( variate =  0, degree = 2, basis_idx = 0 ), second = 0.0, delta = 1e-12 )
        self.assertAlmostEqual( first = evaluateLagrangeBasis1D( variate =  0, degree = 2, basis_idx = 1 ), second = 1.0, delta = 1e-12 )
        self.assertAlmostEqual( first = evaluateLagrangeBasis1D( variate =  0, degree = 2, basis_idx = 2 ), second = 0.0, delta = 1e-12 )
        self.assertAlmostEqual( first = evaluateLagrangeBasis1D( variate = +1, degree = 2, basis_idx = 0 ), second = 0.0, delta = 1e-12 )
        self.assertAlmostEqual( first = evaluateLagrangeBasis1D( variate = +1, degree = 2, basis_idx = 1 ), second = 0.0, delta = 1e-12 )
        self.assertAlmostEqual( first = evaluateLagrangeBasis1D( variate = +1, degree = 2, basis_idx = 2 ), second = 1.0, delta = 1e-12 )

def evaluateLagrangeBasis1D(variate,degree,basis_idx):
    step = 2/degree
    xj = np.arange(-1,2,step)
    val = 1 
    for j in range(0,degree+1):
        if j == basis_idx:
            val = val
        else:
            val = val*(variate - xj[j])/(xj[basis_idx] - xj[j])
    
    return val

unittest.main()

############################################################# 41 Bernstein (requires domain [-1,1] instead of [0,1])
import unittest

class Test_evaluateBernsteinBasis1D( unittest.TestCase ):
    def test_linearBernstein( self ):
        self.assertAlmostEqual( first = evaluateBernsteinBasis1D( variate = -1, degree = 1, basis_idx = 0 ), second = 1.0, delta = 1e-12 )
        self.assertAlmostEqual( first = evaluateBernsteinBasis1D( variate = -1, degree = 1, basis_idx = 1 ), second = 0.0, delta = 1e-12 )
        self.assertAlmostEqual( first = evaluateBernsteinBasis1D( variate = +1, degree = 1, basis_idx = 0 ), second = 0.0, delta = 1e-12 )
        self.assertAlmostEqual( first = evaluateBernsteinBasis1D( variate = +1, degree = 1, basis_idx = 1 ), second = 1.0, delta = 1e-12 )

    def test_quadraticBernstein( self ):
        self.assertAlmostEqual( first = evaluateBernsteinBasis1D( variate = -1, degree = 2, basis_idx = 0 ), second = 1.00, delta = 1e-12 )
        self.assertAlmostEqual( first = evaluateBernsteinBasis1D( variate = -1, degree = 2, basis_idx = 1 ), second = 0.00, delta = 1e-12 )
        self.assertAlmostEqual( first = evaluateBernsteinBasis1D( variate = -1, degree = 2, basis_idx = 2 ), second = 0.00, delta = 1e-12 )
        self.assertAlmostEqual( first = evaluateBernsteinBasis1D( variate =  0, degree = 2, basis_idx = 0 ), second = 0.25, delta = 1e-12 )
        self.assertAlmostEqual( first = evaluateBernsteinBasis1D( variate =  0, degree = 2, basis_idx = 1 ), second = 0.50, delta = 1e-12 )
        self.assertAlmostEqual( first = evaluateBernsteinBasis1D( variate =  0, degree = 2, basis_idx = 2 ), second = 0.25, delta = 1e-12 )
        self.assertAlmostEqual( first = evaluateBernsteinBasis1D( variate = +1, degree = 2, basis_idx = 0 ), second = 0.00, delta = 1e-12 )
        self.assertAlmostEqual( first = evaluateBernsteinBasis1D( variate = +1, degree = 2, basis_idx = 1 ), second = 0.00, delta = 1e-12 )
        self.assertAlmostEqual( first = evaluateBernsteinBasis1D( variate = +1, degree = 2, basis_idx = 2 ), second = 1.00, delta = 1e-12 )

def evaluateBernsteinBasis1D(variate, degree, basis_idx):
    if variate <= 0:
        variate = abs(-1 - variate) / (2)
        val = math.comb(degree,basis_idx)*(variate**basis_idx)*(1-variate)**(degree-basis_idx)
    else:
        val = math.comb(degree,basis_idx)*(variate**basis_idx)*(1-variate)**(degree-basis_idx)
    return val

unittest.main()



