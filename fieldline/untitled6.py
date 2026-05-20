#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  7 15:08:58 2022

@author: weigel
"""

import numpy as np

from scipy.interpolate import NearestNDInterpolator
import matplotlib.pyplot as plt
rng = np.random.default_rng()
x = np.linspace(0, 1., 10)
y = np.linspace(0, 1., 10)
z = np.linspace(0, 1., 10)

def Field(x, y, z):
    return [x, y, z]

[fx,fy,fz] = Field(x, y, z)

interpx = NearestNDInterpolator(list(zip(x, y, z)), fx)
interpy = NearestNDInterpolator(list(zip(x, y, z)), fy)
interpz = NearestNDInterpolator(list(zip(x, y, z)), fy)

def dXds(s, X):
    #return Field(X[0], X[1], X[2])
    fx = interpx(X[0], X[1], X[2])
    fy = interpy(X[0], X[1], X[2])
    fz = interpz(X[0], X[1], X[2])
    
    return [fx, fy, fz]

from scipy.integrate import solve_ivp

soln = solve_ivp(dXds, [0, 10], [0.1, 0.1, 0.1], t_eval=np.linspace(0, 10, 100)).y.transpose()

print(soln)