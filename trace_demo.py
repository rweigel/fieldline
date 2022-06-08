import numpy as np
from os.path import exists
import matplotlib.pyplot as plt
from urllib.request import urlretrieve

from fieldline import trace

urlbase = 'http://mag.gmu.edu/git-data/swmfio/'
tmpdir = '/tmp/'
filebase = '3d__var_2_e20190902-041000-000'

for ext in ['.tree', '.info', '.out']:
    filename = filebase + ext
    if not exists(tmpdir + filename):
        print("Downloading " + urlbase + filename)
        print("to")
        print(tmpdir + filename)
        urlretrieve(urlbase + filename, tmpdir + filename)

filename = tmpdir + filebase + ".out"

import swmfio as swmfio

import logging
swmfio.logger.setLevel(logging.INFO)

if not exists(tmpdir + filename + ".vtk"):
    swmfio.write_vtk(tmpdir + filebase, logger=swmfio.logger)

## set seed points
IC = np.column_stack([-np.ones(17), np.linspace(-1.,1.,17), np.linspace(-1.,1.,17)])

## trace 
scipy_NearestNeighbor_method = trace.trace_file(IC, filename, method='scipy', debug=False)
vtk_method                   = trace.trace_file(IC, filename[:-4] + '.vtk', method='vtk', debug=False)

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
