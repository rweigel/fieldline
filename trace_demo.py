def demo1():
    fname = '/home/gary/temp/'+'3d__var_3_e20031120-070000-000.vtu'
    #IC = np.array([1.,1.,1.])
    IC = np.column_stack([np.ones(20), np.ones(20), np.ones(20)])
    print(trace_file(IC, fname, method='vtk', debug=True))
