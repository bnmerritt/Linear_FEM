# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 21:38:05 2022

@author: brian
"""

import numpy
import unittest
import math

class Test_getNewtonCotesQuadrature( unittest.TestCase ):
    def test_incorrect_num_points( self ):
        with self.assertRaises( Exception ) as context:
            getNewtonCotesQuadrature( num_points = 0 )
        self.assertEqual( "num_points_MUST_BE_INTEGER_IN_[1,6]", str( context.exception ) )
        with self.assertRaises( Exception ) as context:
            getNewtonCotesQuadrature( num_points = 7 )
        self.assertEqual( "num_points_MUST_BE_INTEGER_IN_[1,6]", str( context.exception ) )

    def test_return_types( self ):
        for num_points in range( 1, 7 ):
            x, w = getNewtonCotesQuadrature( num_points = num_points )
            self.assertIsInstance( obj = x, cls = numpy.ndarray )
            self.assertIsInstance( obj = w, cls = numpy.ndarray )
            self.assertTrue( len( x ) == num_points )
            self.assertTrue( len( w ) == num_points )

class Test_computeNewtonCotesQuadrature( unittest.TestCase ):
    def test_integrate_constant_one( self ):
        constant_one = lambda x : 1 * x**0
        for degree in range( 1, 6 ):
            num_points = degree + 1
            self.assertAlmostEqual( first = computeNewtonCotesQuadrature( fun = constant_one, num_points = num_points ), second = 2.0, delta = 1e-12 )

    def test_exact_poly_int( self ):
        for degree in range( 1, 6 ):
            num_points = degree + 1
            poly_fun = lambda x : ( x + 1.0 ) ** degree
            indef_int = lambda x : ( ( x + 1 ) ** ( degree + 1) ) / ( degree + 1 )
            def_int = indef_int(1.0) - indef_int(-1.0)
            self.assertAlmostEqual( first = computeNewtonCotesQuadrature( fun = poly_fun, num_points = num_points ), second = def_int, delta = 1e-12 )

    def test_integrate_sin( self ):
        sin = lambda x : math.sin(x)
        for num_points in range( 1, 7 ):
            self.assertAlmostEqual( first = computeNewtonCotesQuadrature( fun = sin, num_points = num_points ), second = 0.0, delta = 1e-12 )

    def test_integrate_cos( self ):
        cos = lambda x : math.cos(x)
        self.assertAlmostEqual( first = computeNewtonCotesQuadrature( fun = cos, num_points = 6 ), second = 2*math.sin(1), delta = 1e-4 )
        
# def evaluateLagrangeBasis1D(variate,degree,basis_idx):
#     step = 2/degree
#     xj = numpy.arange(-1,1,step)
#     val = 1 
#     for j in range(0,degree+1):
#         if j == basis_idx:
#             val = val
#         else:
#             val = val*(variate - xj[j])/(xj[basis_idx] - xj[j])
#     return val
            
def getNewtonCotesQuadrature(num_points):
    if num_points <= 0:
        raise ValueError("num_points_MUST_BE_INTEGER_IN_[1,6]")
    elif num_points >= 7:
        raise ValueError("num_points_MUST_BE_INTEGER_IN_[1,6]")
    else:
        # width = 2/num_points
        # w = numpy.full(num_points, width)
        # x = numpy.arange(-1 + width/2, 1, width)
        width = 2/(num_points+1)
        w = numpy.full(num_points+1, width)
        x = numpy.arange(-1, 1 + width, width)
    return x, w 


def computeNewtonCotesQuadrature(fun,num_points):
    x,w = getNewtonCotesQuadrature(num_points)
    y = numpy.zeros(num_points)
    for i in range(0,len(x)):
        y[i] = fun(x[i])
        # y[i] = evaluateLagrangeBasis1D(x[i],len(x)-1,basis_idx)
    areas = numpy.multiply(w,y)
    integral = sum(areas)
    return integral
unittest.main()          