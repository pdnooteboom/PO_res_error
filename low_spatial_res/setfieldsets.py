from parcels import FieldSet


def set_pop_fieldset(ufiles, dimfiles, dfile, bfile, afile, indices = None):#
    filenames = { 'U': {'lon': bfile,
                        'lat': bfile,
                        'depth': dfile,
                        'data':ufiles},
                'V' : {'lon': bfile,
                        'lat': bfile,
                        'depth': dfile,
                        'data':ufiles},
                'W' : {'lon': bfile,
                        'lat': bfile,
                        'depth': dfile,
                        'data':ufiles},
                'S' : {'lon': bfile,
                        'lat': bfile,
                        'depth': dfile,
                        'data':ufiles},
                'T' : {'lon': bfile,
                        'lat': bfile,
                        'depth': dfile,
                        'data':ufiles}   ,
                'B' : {'lon': bfile,
                        'lat': bfile,
                        'depth': dfile,
                        'data':bfile},
                'cell_areas' : {'lon': bfile,
                        'lat': bfile,
                        'depth': dfile,
                        'data':afile},
                }

    variables = {'U': 'UVEL',
                 'V': 'VVEL',
                 'W': 'WVEL',
                 'cell_areas': 'UAREA',
                 'T': 'TEMP',
                 'S': 'SALT',
                 'B':'Bathymetry'}

    dimensions = {'U':{'lon': 'ULONG', 'lat': 'ULAT', 'depth': 'w_dep', 'time': 'time'},
                  'V': {'lon': 'ULONG', 'lat': 'ULAT', 'depth': 'w_dep', 'time': 'time'},
                    'W': {'lon': 'ULONG', 'lat': 'ULAT', 'depth': 'w_dep', 'time': 'time'},
                  'cell_areas':{'lon': 'ULONG', 'lat': 'ULAT', 'depth': 'w_dep', 'time': 'time'},
                    'T': {'lon': 'ULONG', 'lat': 'ULAT', 'depth': 'w_dep', 'time': 'time'},
                    'S': {'lon': 'ULONG', 'lat': 'ULAT', 'depth': 'w_dep', 'time': 'time'},
                    'B': {'lon': 'ULONG', 'lat': 'ULAT'} }


    if(indices!=None):
        fieldset = FieldSet.from_pop(filenames, variables, dimensions, indices=indices, allow_time_extrapolation=False)
    else:
        fieldset = FieldSet.from_pop(filenames, variables, dimensions, allow_time_extrapolation=False)

    fieldset.U.vmax = 10    # set max of flow to 10 m/s
    fieldset.V.vmax = 10
    fieldset.W.vmax = 10
    fieldset.T.vmin = -5
    fieldset.cell_areas.allow_time_extrapolation = True
    fieldset.B.allow_time_extrapolation = True
    fieldset.cell_areas.set_scaling_factor(0.0001)  # cm^2 to m^2

    return fieldset


def set_pop_fieldset_bolus(ufiles, dimfiles, dfile, bfile, afile, indices=None):#
    filenames = { 'U': {'lon': bfile,
                        'lat': bfile,
                        'depth': dfile,
                        'data':ufiles},
                'V' : {'lon': bfile,
                        'lat': bfile,
                        'depth': dfile,
                        'data':ufiles},
                'W' : {'lon': bfile,
                        'lat': bfile,
                        'depth': dfile,
                        'data':ufiles},
                'Ubolus': {'lon': bfile,
                        'lat': bfile,
                        'depth': dfile,
                        'data':ufiles},
                'Vbolus' : {'lon': bfile,
                        'lat': bfile,
                        'depth': dfile,
                        'data':ufiles},
                'S' : {'lon': bfile,
                        'lat': bfile,
                        'depth': dfile,
                        'data':ufiles},
                'T' : {'lon': bfile,
                        'lat': bfile,
                        'depth': dfile,
                        'data':ufiles}   ,
                'B' : {'lon': bfile,
                        'lat': bfile,
                        'depth': dfile,
                        'data':bfile},
                'cell_areas' : {'lon': bfile,
                        'lat': bfile,
                        'depth': dfile,
                        'data':afile},
                }
    variables = {'U': 'UVEL',
                 'V': 'VVEL',
                 'W': 'WVEL',
                 'cell_areas': 'UAREA',
                 'Ubolus': 'UISOP',
                 'Vbolus': 'VISOP',
                 'T': 'TEMP',
                 'S': 'SALT',
                 'B':'Bathymetry'}
    dimensions = {'U':{'lon': 'ULONG', 'lat': 'ULAT', 'depth': 'w_dep', 'time': 'time'},
                  'V': {'lon': 'ULONG', 'lat': 'ULAT', 'depth': 'w_dep', 'time': 'time'},
                    'W': {'lon': 'ULONG', 'lat': 'ULAT', 'depth': 'w_dep', 'time': 'time'},
                  'cell_areas':{'lon': 'ULONG', 'lat': 'ULAT', 'depth': 'w_dep', 'time': 'time'},
                  'Ubolus':{'lon': 'ULONG', 'lat': 'ULAT', 'depth': 'w_dep', 'time': 'time'},
                  'Vbolus': {'lon': 'ULONG', 'lat': 'ULAT', 'depth': 'w_dep', 'time': 'time'},
                    'T': {'lon': 'ULONG', 'lat': 'ULAT', 'depth': 'w_dep', 'time': 'time'},
                    'S': {'lon': 'ULONG', 'lat': 'ULAT', 'depth': 'w_dep', 'time': 'time'},
                    'B': {'lon': 'ULONG', 'lat': 'ULAT'} }

    if(indices!=None):
        fieldset = FieldSet.from_pop(filenames, variables, dimensions, indices=indices, allow_time_extrapolation=False)
    else:
        fieldset = FieldSet.from_pop(filenames, variables, dimensions, allow_time_extrapolation=False)

    fieldset.U.vmax = 10    # set max of flow to 10 m/s
    fieldset.V.vmax = 10
    fieldset.W.vmax = 10
    fieldset.Ubolus.set_scaling_factor(0.01)  # cm/s to m/s
    fieldset.Vbolus.set_scaling_factor(0.01)  # cm/s to m/s
    fieldset.cell_areas.set_scaling_factor(0.0001)  # cm^2/s to m^2/s
    fieldset.T.vmin = -5
    fieldset.cell_areas.allow_time_extrapolation = True
    fieldset.B.allow_time_extrapolation = True

    return fieldset

