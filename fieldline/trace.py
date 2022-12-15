import os
import sys
import numpy as np
import vtk
from vtk.numpy_interface import dataset_adapter as dsa
from vtk.vtkCommonDataModel import vtkDataSet
from vtk.vtkCommonExecutionModel import vtkAlgorithmOutput


#def analytic(IC, Field, integration_direction='backward', debug=False):
def trace(IC, Field, integration_direction='backward', debug=False):

    from scipy.integrate import solve_ivp
    #from scipy.integrate import odeint
    import types
    assert isinstance(Field, types.FunctionType)

    if integration_direction in ['northern', 'negative', 'backward']:
        sign = -1
    elif integration_direction in ['southern', 'positive', 'forward']:
        sign = +1
    else:
        raise ValueError('"' + str(integration_direction)
                         + '" is not a valid integration_direction value')

    def dXds(s, X):
        F = Field(X)
        Fmag = np.linalg.norm(F)
        if 1e-9 < Fmag < 1e+7:
            return (sign/Fmag)*F
        return [0., 0., 0.]

    s_grid = np.arange(0., 10., 0.1)
    max_iterations = 100

    if IC.shape == (3,):
        IC = [IC]
    ret = []
    linenum = 0
    for X0 in list(IC):
        if debug:
            print('linenum = ' + str(linenum))
        done = False
        solns = np.empty((0, 3)) # Combined solutions
        i = 0
        while not done:
            if debug:
                print('i = ' + str(i))

            #soln = odeint(dXds, X0, s_grid, tfirst=True)
            soln = solve_ivp(dXds, [0, 10], X0, t_eval=s_grid).y.transpose()
            if debug: print('hellothere')

            R = soln[:, 0]**2 + soln[:, 1]**2 + soln[:, 2]**2
            # define condition on the field line points
            # Find first location where soln steps out-of-bounds
            #tr = np.where( False == (R >= 1) & (soln[:,0] > -30.) & (np.abs(soln[:, 2]) < 20.) )        
            # Boolean array.

            tr = (R >= 1) & (soln[:,0] > -30.) & (np.abs(soln[:, 2]) < 20.)
            # RuntimeWarning: invalid value encountered in greater_equal

            # Indices where stop conditions satisfied
            tr_out = np.where(tr == False)
            if debug:
                print(tr)
            if tr_out[0].size > 0:
                # Stop condition found at least once. Use solution up to that point.s
                solns = np.vstack((solns, soln[0:tr_out[0][0] + 1, :]))
                done = True
            elif max_iterations == i + 1:
                solns = np.vstack((solns, soln))
                done = True
            else:
                # New initial condition is stop point.
                X0 = soln[-1, :]
                # Append solution but exclude last value, which is the
                # new initial condition.
                solns = np.vstack((solns, soln[0:-1, :]))
            i = i + 1
        ret.append(solns)
        linenum += 1

    if len(ret) == 1:
        return ret[0]
    return ret


#def points(IC, Field, Domain, integration_direction='backward', debug=False):
def interpolate_and_trace(IC, Field, Domain, integration_direction='backward', debug=False):
    from scipy.interpolate import RegularGridInterpolator, NearestNDInterpolator

    Field = np.array(Field)
    Domain = np.array(Domain)

    if Domain.shape == (3,2) or Domain.shape == (2,3):
        #TODO
        #xlims = ... TODO

        #from make_grid import make_axes
        #ax_list = make_axes(xlims,ylims,zlims,d)

        # https://stackoverflow.com/questions/21836067/interpolate-3d-volume-with-numpy-and-or-scipy
        Fx_interp = RegularGridInterpolator(tuple(ax_list), Field[0, :,:,:])
        Fy_interp = RegularGridInterpolator(tuple(ax_list), Field[1, :,:,:])
        Fz_interp = RegularGridInterpolator(tuple(ax_list), Field[2, :,:,:])

        def Fcallable(v):
            return np.array([Fx_interp(v)[0], Fy_interp(v)[0], Fz_interp(v)[0]])

    else:
        assert(Domain.shape == Field.shape)
        if Domain.shape[0]==3 and Domain.shape[1]!=3:
            Domain = Domain.transpose()
            Field = Field.transpose()

        Fx_interp = NearestNDInterpolator(Domain,Field[:,0])
        Fy_interp = NearestNDInterpolator(Domain,Field[:,1])
        Fz_interp = NearestNDInterpolator(Domain,Field[:,2])

        def Fcallable(v):
            return np.array([Fx_interp(v)[0], Fy_interp(v)[0], Fz_interp(v)[0]])

    return trace(IC, Fcallable, integration_direction=integration_direction, debug=debug)

# see magnetosphere/misc/vtk/streamline_from_datafile_demo.py
# def vtk(IC, vtk_object, integration_direction='backward', debug=False, var='b', celldata=True):
def trace_vtk(IC, vtk_object, integration_direction='backward', debug=False, var='b', celldata=True):
    if integration_direction in ['northern', 'negative', 'backward']:
        vtk_int_dir = vtk.VTK_INTEGRATE_BACKWARD # = 1
    elif integration_direction in ['southern', 'positive', 'forward']:
        vtk_int_dir = vtk.VTK_INTEGRATE_FORWARD # = 0
    elif integration_direction in ['both']:
        vtk_int_dir = vtk.VTK_INTEGRATE_BOTH_DIRECTIONS # = 2
    else:
        raise ValueError(str(integration_direction)+' not a valid integration_direction')

    if celldata:
        vtk_object.GetCellData().SetActiveVectors(var)
        c2p = vtk.vtkCellDataToPointData()
        c2p.SetInputDataObject(vtk_object)
        tostreamer = c2p.GetOutputPort()
    else:
        vtk_object.GetPointData().SetActiveVectors(var)
        tostreamer = vtk_object

    if IC.shape == (3,):
        IC = [IC]
    ret = []
    linenum = 0
    for X0 in list(IC):
        # Create integrator
        rk = vtk.vtkRungeKutta45()
        # Create source for streamtubes
        streamer = vtk.vtkStreamTracer()

        if isinstance(tostreamer, vtkDataSet):
            streamer.SetInputDataObject(tostreamer)
        elif isinstance(tostreamer, vtkAlgorithmOutput):
            streamer.SetInputConnection(tostreamer)
            #raise RuntimeWarning ('cannot at present select the apropriate vector field. '\
            #                      +  'will use whatever defaults.')
        else:
            raise ValueError ('not a supported vtk object for streamlines')

        streamer.SetStartPosition(X0) #cannot pass multiple IC's in an array
        streamer.SetMaximumPropagation(400) ###
        #streamer.SetIntegrationStepUnit(2) # apperars overiden by next lines, see https://vtk.org/doc/nightly/html/classvtkStreamTracer.html#afe365e81e110f354065f5adc8401d589
        streamer.SetMinimumIntegrationStep(0.00001)
        streamer.SetMaximumIntegrationStep(0.2)
        streamer.SetInitialIntegrationStep(0.01)
        streamer.SetIntegrationDirection(vtk_int_dir)
        streamer.SetIntegrator(rk)
        streamer.SetRotationScale(0.5)
        streamer.SetMaximumError(1.0e-5)
        #https://stackoverflow.com/questions/38504907/reading-a-vtk-polydata-file-and-converting-it-into-numpy-array
        ## either order works ##
        polydata = streamer.GetOutput()
        streamer.Update() # forces the computation of the stream lines

        ret.append( dsa.WrapDataObject(polydata).Points ) # convert the result to array and put in list
        del polydata 
        del streamer

    return ret

# def file(IC, filename, method='vtk', integration_direction='backward', debug=False):
def trace_file(IC, filename, method='vtk', integration_direction='backward', debug=False):
    if not os.path.exists(filename): raise FileNotFoundError ('no file ' + fname)
    ext = filename[-4:]

    if method == 'scipy':
        if ext == '.out':
            import spacepy.pybats.bats as bats
            data3d = bats.Bats2d(filename)

            X = np.nan*np.empty((data3d['x'].size,3), dtype=np.float32)
            X[:,0] = data3d['x']
            X[:,1] = data3d['y']
            X[:,2] = data3d['z']

            B = np.nan*np.empty(X.shape, dtype=np.float32)
            B[:,0] = data3d['bx']
            B[:,1] = data3d['by']
            B[:,2] = data3d['bz']

            return interpolate_and_trace(IC, B, X, integration_direction=integration_direction, debug=debug)
        elif ext == '.cdf':
            pass # TODO: intepolate using kameleon onto regular grid then pass to interpolate_and_trace
        else:
            raise ValueError

    elif method == 'vtk':
        if ext in ['.vtk','.vtu']:
            ## open file
            if ext == '.vtk':
                reader = vtk.vtkGenericDataObjectReader()
            elif ext == '.vtu':
                reader = vtk.vtkXMLUnstructuredGridReader() # https://stackoverflow.com/questions/54044958/reading-data-from-a-raw-vtk-vtu-file
            reader.SetFileName(filename)

            reader.Update() # forces loading the file into memory
            # get nessesary vtk info from Reader object
            reader_output = reader.GetOutput() #not the "New Pipelin"

            return trace_vtk(IC, reader_output, integration_direction=integration_direction, debug=debug)
        elif ext == '.out':
            VTKfilename = filename[:-4]+'.vtk'
            return traceFile(VTKfilename, IC, method=method, integration_direction=integration_direction, debug=debug)
        elif ext == '.cdf':
            raise RuntimeWarning ('should it be allowed?') #!!!!
        else:
            raise ValueError

    elif method == 'kameleon':
        if ext!='.cdf' : raise ValueError('need a CDF file, currently no way to make it')
        #TODO: use kameleons field line tracer.

