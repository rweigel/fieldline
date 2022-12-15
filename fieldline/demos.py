import os
import sys
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

from fieldline import trace


def generate_northern_hemisphere(R,Nt,Np):
    theta = np.linspace(0,np.pi/2.,Nt+2)[1:-1]
    phi = np.linspace(0,2.*np.pi,Np+2)[1:-1]

    B1, B2 = np.meshgrid(phi, theta)
    B1 = B1.flatten(order='C')
    B2 = B2.flatten(order='C')

    PI = np.pi*np.ones((B1.size, ))
    x = R*np.cos(B1+PI)*np.sin(B2)
    y = R*np.sin(B1+PI)*np.sin(B2)
    z = R*np.cos(B2)
    XYZ = np.column_stack((x, y, z))

    return XYZ

def urlretrieve(url, fname):
    print('\ndownloading %s from %s\n'%(fname,url))
    import urllib3
    import shutil
    http = urllib3.PoolManager()

    with http.request('GET', url, preload_content=False) as req:
        if req.status == 200:
            with open(fname,'wb') as fl:
                shutil.copyfileobj(req, fl)
        else:
            raise ConnectionError ('error in dowloading from '+url)

def demo1(fname=None):
    '''
    inputs: fname
        set to the filename (including path if specified) of the swmf .out file
        e.g. '/home/gary/temp/3d__var_3_e20031120-070000-000.out'
        default: None, which downloads an example file from the internet
    '''
    ## download example data file if none provided
    if fname is None:
        fname = '/tmp/3d__var_2_e20190902-041000-000.out' #!!! doesn't work for windows
        if not os.path.exists(fname):
            urlname = 'http://mag.gmu.edu/git-data/GaryQ-Physics/demodata/3d__var_2_e20190902-041000-000.out'
            urlretrieve(urlname, fname)

    ## set seed points
    IC = np.column_stack([-np.ones(17), np.linspace(-1.,1.,17), np.linspace(-1.,1.,17)])

    ## trace 
    scipy_NearestNeighbor_method = trace.trace_file(IC, fname, method='scipy', debug=False)

    ## plot the result
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    for i in range(IC.shape[0]):
        ax.plot(scipy_NearestNeighbor_method[i][:,0],
                scipy_NearestNeighbor_method[i][:,1],
                scipy_NearestNeighbor_method[i][:,2],
                color='b',marker='.',alpha=0.5)

    ax.plot(IC[:,0],
            IC[:,1],
            IC[:,2],
            color='r',marker='o')

    ax.plot(np.array([0.]),
            np.array([0.]),
            np.array([0.]),
            color='r',
            marker='o')

    ax.set_xlim(-4.,4.)
    ax.set_ylim(-4.,4.)
    ax.set_zlim(-4.,4.)

    plt.show()
    plt.clf()


def demo2(ftag=None, writevtk=False):
    '''
    traces field line with seed points along spaced out on 
    the positive Z (in GSM) half of a sphere of radius 2 R_E.

    inputs:
      ftag:
        set to the filename (including path if specified) of the swmf files, without the extension.
        e.g. '/home/gary/temp/3d__var_3_e20031120-070000-000'
        default: None, which downloads example files from the internet
      writevtk:
        boolean, default False. If True, it writes out a Polyline vtk file for each field line
        that can be loaded into paraview. Otherwise it doesn't, and just plots with matplotlib.

    in resuling plot:
        red   : shows the seed points
        blue  : shows the result of tracing via scipy's nearest neighbor interpolator on native grid
        green : shows the result of vtk's interpolation.
    '''
    from swmf_file_reader import read_swmf_files as rswmf
    from swmf_file_reader.vtk_export_copy import vtk_export

    ## download example data files if none provided
    if ftag is None:
        ftag = '/tmp/3d__var_2_e20190902-041000-000' #!!! doesn't work for windows
        if not os.path.exists(ftag+'.info'):
            urltag = 'http://mag.gmu.edu/git-data/GaryQ-Physics/demodata/3d__var_2_e20190902-041000-000'
            urlretrieve(urltag+'.info', ftag+'.info')
            urlretrieve(urltag+'.tree', ftag+'.tree')
            urlretrieve(urltag+'.out' , ftag+'.out' )

    ## set seed points
    #IC = np.column_stack([-np.ones(17), np.linspace(-2.,2.,17), np.linspace(-2.,2.,17)])
    IC = generate_northern_hemisphere(2.,7,7)

    ## if not already existing, make vtk file from swmf .out file, to be used in tracing
    if not os.path.exists(ftag+'.vtk'):
        print('\n\ngenerating vtk file\n\n')
        rswmf.swmf2vtk(ftag)
    else:
        print('\n\nusing existing vtk file\n\n')

    ## trace 
    scipy_NearestNeighbor_method = trace.trace_file(IC, ftag+'.out', 
                                    method='scipy', integration_direction='northern')
    vtk_method                   = trace.trace_file(IC, ftag+'.vtk',
                                    method='vtk',   integration_direction='northern')

    ## plot the result
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    for i in range(IC.shape[0]):
        ax.plot(scipy_NearestNeighbor_method[i][:,0],
                scipy_NearestNeighbor_method[i][:,1],
                scipy_NearestNeighbor_method[i][:,2],
                color='b',marker='.',alpha=0.5)

        ax.plot(vtk_method[i][:,0],
                vtk_method[i][:,1],
                vtk_method[i][:,2],
                color='g',marker='.',alpha=0.5)

        if writevtk:
            print('i=%d'%(i))
            print('IC[i,:]='+str(IC[i,:]))
            vtk_export('NN_%.2d.vtk'%(i), scipy_NearestNeighbor_method[i],
                                dataset = 'POLYDATA',
                                connectivity = 'LINES',ftype='ASCII')
            vtk_export('VTK_%.2d.vtk'%(i), vtk_method[i],
                                dataset = 'POLYDATA',
                                connectivity = 'LINES',ftype='ASCII')

    ax.scatter(IC[:,0],
               IC[:,1],
               IC[:,2],
               color='r',
               marker='o')

    ax.set_xlim(-5.,5.)
    ax.set_ylim(-5.,5.)
    ax.set_zlim(-5.,5.)

    plt.show()
    plt.clf()


def demo3(fname=None):
    # Not working
    kameleon_method = trace.trace_file(IC, fname+'.cdf', method='kameleon', debug=True)
