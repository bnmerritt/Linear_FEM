################################################################################# Question 33
# import unittest

# class Test_evaluateMonomialBasis1D( unittest.TestCase ):
#     def test_basisAtBounds( self ):
#         self.assertAlmostEqual( first = evaluateMonomialBasis1D( degree = 0, variate = 0 ), second = 1.0, delta = 1e-12 )
#         for p in range( 1, 11 ):
#             self.assertAlmostEqual( first = evaluateMonomialBasis1D( degree = p, variate = 0 ), second = 0.0, delta = 1e-12 )
#             self.assertAlmostEqual( first = evaluateMonomialBasis1D( degree = p, variate = 1 ), second = 1.0, delta = 1e-12 )

#     def test_basisAtMidpoint( self ):
#         for p in range( 0, 11 ):
#             self.assertAlmostEqual( first = evaluateMonomialBasis1D( degree = p, variate = 0.5 ), second = 1 / ( 2**p ), delta = 1e-12 )

# def evaluateMonomialBasis1D(degree,variate):
#     return variate**degree

# unittest.main()


################################################################################# Question 34
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

def evaluateMonomialBasis1D(degree,variate):
    return variate**degree

unittest.main()
################################################################################# Question 36
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