# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 13:07:50 2022

@author: brian
"""
import numpy
import math
import matplotlib.pyplot as plt

alpha_l = [1.57080,4.71239,7.85398,10.9956,14.1372]
points = 100
alpha_x0 = numpy.linspace(0+(alpha_l[4]/100),alpha_l[4],points)
phi = numpy.zeros(len(alpha_x0))
for i in range(0,len(alpha_x0+1)):
    phi[i] = math.cosh(alpha_x0[i]) - math.cos(alpha_x0[i]) - (math.cosh(alpha_x0[i]) + math.cos(alpha_x0[i])) / (math.sinh(alpha_x0[i]) + math.sin(alpha_x0[i]))*(math.sinh(alpha_x0[i])-math.sin(alpha_x0[i]))
plt.plot(alpha_x0,phi)