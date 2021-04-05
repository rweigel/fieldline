import os
import sys
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

from fieldline import trace

def demo(fname):
    '''
    traces field line with seed points along a line connecting (-1,1,1) and (-1,-1,-1) in GSM
    in resuling plot:
        red   : shows the seed points line, and the location of the origin (0,0,0)
        blue  : shows the result of tracing via scipy's nearest neighbor interpolator on native grid
        green : shows the result of vtk's interpolation.
    '''
    ## set seed points
    IC = np.column_stack([-np.ones(17), np.linspace(-1.,1.,17), np.linspace(-1.,1.,17)])

    ## trace 
    scipy_NearestNeighbor_method = trace.trace_file(IC, fname, method='scipy', debug=False)
    vtk_method                   = trace.trace_file(IC, fname[:-4]+'.vtk', method='vtk', debug=False)

    ## plot the result
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    for i in range(IC.shape[0]):
        ax.plot(scipy_NearestNeighbor_method[i][:,0],
                scipy_NearestNeighbor_method[i][:,1],
                scipy_NearestNeighbor_method[i][:,2],
                color='b')

        ax.plot(vtk_method[i][:,0],
                vtk_method[i][:,1],
                vtk_method[i][:,2],
                color='g')

    ax.plot(IC[:,0],
            IC[:,1],
            IC[:,2],
            color='r')

    ax.plot(np.array([0.]),
            np.array([0.]),
            np.array([0.]),
            color='r',
            marker='o')

    ax.set_xlim(-4.,4.)
    ax.set_ylim(-4.,4.)
    ax.set_zlim(-4.,4.)

    plt.show()


if __name__=='__main__':
    fname = '/home/gary/temp/'+'3d__var_3_e20031120-070000-000.out'
    demo(fname)


