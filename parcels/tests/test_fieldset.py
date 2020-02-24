from parcels import FieldSet, ParticleSet, ScipyParticle, JITParticle, Variable, AdvectionRK4, AdvectionRK4_3D, RectilinearZGrid, ErrorCode
from parcels.field import Field, VectorField
from parcels.tools.converters import TimeConverter, _get_cftime_calendars, _get_cftime_datetimes
from datetime import timedelta as delta
import datetime
import numpy as np
import xarray as xr
import math
import pytest
from os import path
import cftime


ptype = {'scipy': ScipyParticle, 'jit': JITParticle}


def generate_fieldset(xdim, ydim, zdim=1, tdim=1):
    lon = np.linspace(0., 10., xdim, dtype=np.float32)
    lat = np.linspace(0., 10., ydim, dtype=np.float32)
    depth = np.zeros(zdim, dtype=np.float32)
    time = np.zeros(tdim, dtype=np.float64)
    if zdim == 1 and tdim == 1:
        U, V = np.meshgrid(lon, lat)
        dimensions = {'lat': lat, 'lon': lon}
    else:
        U = np.ones((tdim, zdim, ydim, xdim))
        V = np.ones((tdim, zdim, ydim, xdim))
        dimensions = {'lat': lat, 'lon': lon, 'depth': depth, 'time': time}
    data = {'U': np.array(U, dtype=np.float32), 'V': np.array(V, dtype=np.float32)}

    return (data, dimensions)


@pytest.mark.parametrize('xdim', [100, 200])
@pytest.mark.parametrize('ydim', [100, 200])
def test_fieldset_from_data(xdim, ydim):
    """ Simple test for fieldset initialisation from data. """
    data, dimensions = generate_fieldset(xdim, ydim)
    fieldset = FieldSet.from_data(data, dimensions)
    assert len(fieldset.U.data.shape) == 3
    assert len(fieldset.V.data.shape) == 3
    assert np.allclose(fieldset.U.data[0, :], data['U'], rtol=1e-12)
    assert np.allclose(fieldset.V.data[0, :], data['V'], rtol=1e-12)


@pytest.mark.parametrize('ttype', ['float', 'datetime64'])
@pytest.mark.parametrize('tdim', [1, 20])
def test_fieldset_from_data_timedims(ttype, tdim):
    data, dimensions = generate_fieldset(10, 10, tdim=tdim)
    if ttype == 'float':
        dimensions['time'] = np.linspace(0, 5, tdim)
    else:
        dimensions['time'] = [np.datetime64('2018-01-01') + np.timedelta64(t, 'D') for t in range(tdim)]
    fieldset = FieldSet.from_data(data, dimensions)
    for i, dtime in enumerate(dimensions['time']):
        assert fieldset.U.grid.time_origin.fulltime(fieldset.U.grid.time[i]) == dtime


@pytest.mark.parametrize('xdim', [100, 200])
@pytest.mark.parametrize('ydim', [100, 50])
def test_fieldset_from_data_different_dimensions(xdim, ydim, zdim=4, tdim=2):
    """ Test for fieldset initialisation from data using
    dict-of-dict for dimensions. """

    lon = np.linspace(0., 1., xdim, dtype=np.float32)
    lat = np.linspace(0., 1., ydim, dtype=np.float32)
    depth = np.zeros(zdim, dtype=np.float32)
    time = np.zeros(tdim, dtype=np.float64)
    U = np.zeros((xdim, ydim), dtype=np.float32)
    V = np.ones((xdim, ydim), dtype=np.float32)
    P = 2 * np.ones((int(xdim/2), int(ydim/2), zdim, tdim), dtype=np.float32)
    data = {'U': U, 'V': V, 'P': P}
    dimensions = {'U': {'lat': lat, 'lon': lon},
                  'V': {'lat': lat, 'lon': lon},
                  'P': {'lat': lat[0::2], 'lon': lon[0::2], 'depth': depth, 'time': time}}

    fieldset = FieldSet.from_data(data, dimensions, transpose=True)
    assert len(fieldset.U.data.shape) == 3
    assert len(fieldset.V.data.shape) == 3
    assert len(fieldset.P.data.shape) == 4
    assert fieldset.P.data.shape == (tdim, zdim, ydim/2, xdim/2)
    assert np.allclose(fieldset.U.data, 0., rtol=1e-12)
    assert np.allclose(fieldset.V.data, 1., rtol=1e-12)
    assert np.allclose(fieldset.P.data, 2., rtol=1e-12)


@pytest.mark.parametrize('xdim', [100, 200])
@pytest.mark.parametrize('ydim', [100, 200])
def test_fieldset_from_parcels(xdim, ydim, tmpdir, filename='test_parcels'):
    """ Simple test for fieldset initialisation from Parcels FieldSet file format. """
    filepath = tmpdir.join(filename)
    data, dimensions = generate_fieldset(xdim, ydim)
    fieldset_out = FieldSet.from_data(data, dimensions)
    fieldset_out.write(filepath)
    fieldset = FieldSet.from_parcels(filepath)
    assert len(fieldset.U.data.shape) == 3  # Will be 4 once we use depth
    assert len(fieldset.V.data.shape) == 3
    assert np.allclose(fieldset.U.data[0, :], data['U'], rtol=1e-12)
    assert np.allclose(fieldset.V.data[0, :], data['V'], rtol=1e-12)


@pytest.mark.parametrize('calendar, cftime_datetime',
                         zip(_get_cftime_calendars(),
                             _get_cftime_datetimes()))
def test_fieldset_nonstandardtime(calendar, cftime_datetime, tmpdir, filename='test_nonstandardtime.nc', xdim=4, ydim=6):
    filepath = tmpdir.join(filename)
    dates = [getattr(cftime, cftime_datetime)(1, m, 1) for m in range(1, 13)]
    da = xr.DataArray(np.random.rand(12, xdim, ydim),
                      coords=[dates, range(xdim), range(ydim)],
                      dims=['time', 'lon', 'lat'], name='U')
    da.to_netcdf(str(filepath))

    dims = {'lon': 'lon', 'lat': 'lat', 'time': 'time'}
    field = Field.from_netcdf(filepath, 'U', dims)
    assert field.grid.time_origin.calendar == calendar


def test_field_from_netcdf():
    data_path = path.join(path.dirname(__file__), 'test_data/')

    filenames = {'lon': data_path + 'mask_nemo_cross_180lon.nc',
                 'lat': data_path + 'mask_nemo_cross_180lon.nc',
                 'data': data_path + 'Uu_eastward_nemo_cross_180lon.nc'}
    variable = 'U'
    dimensions = {'lon': 'glamf', 'lat': 'gphif'}
    Field.from_netcdf(filenames, variable, dimensions, interp_method='cgrid_velocity')


@pytest.mark.parametrize('indslon', [range(10, 20), [1]])
@pytest.mark.parametrize('indslat', [range(30, 60), [22]])
def test_fieldset_from_file_subsets(indslon, indslat, tmpdir, filename='test_subsets'):
    """ Test for subsetting fieldset from file using indices dict. """
    data, dimensions = generate_fieldset(100, 100)
    filepath = tmpdir.join(filename)
    fieldsetfull = FieldSet.from_data(data, dimensions)
    fieldsetfull.write(filepath)
    indices = {'lon': indslon, 'lat': indslat}
    indices_back = indices.copy()
    fieldsetsub = FieldSet.from_parcels(filepath, indices=indices)
    assert indices == indices_back
    assert np.allclose(fieldsetsub.U.lon, fieldsetfull.U.grid.lon[indices['lon']])
    assert np.allclose(fieldsetsub.U.lat, fieldsetfull.U.grid.lat[indices['lat']])
    assert np.allclose(fieldsetsub.V.lon, fieldsetfull.V.grid.lon[indices['lon']])
    assert np.allclose(fieldsetsub.V.lat, fieldsetfull.V.grid.lat[indices['lat']])

    ixgrid = np.ix_([0], indices['lat'], indices['lon'])
    assert np.allclose(fieldsetsub.U.data, fieldsetfull.U.data[ixgrid])
    assert np.allclose(fieldsetsub.V.data, fieldsetfull.V.data[ixgrid])


@pytest.mark.parametrize('calltype', ['from_data', 'from_nemo'])
def test_illegal_dimensionsdict(calltype):
    error_thrown = False
    try:
        if calltype == 'from_data':
            data, dimensions = generate_fieldset(10, 10)
            dimensions['test'] = None
            FieldSet.from_data(data, dimensions)
        elif calltype == 'from_nemo':
            fname = path.join(path.dirname(__file__), 'test_data', 'mask_nemo_cross_180lon.nc')
            filenames = {'dx': fname, 'mesh_mask': fname}
            variables = {'dx': 'e1u'}
            dimensions = {'lon': 'glamu', 'lat': 'gphiu', 'test': 'test'}
            error_thrown = False
            FieldSet.from_nemo(filenames, variables, dimensions)
    except NameError:
        error_thrown = True
    assert error_thrown


@pytest.mark.parametrize('xdim', [100, 200])
@pytest.mark.parametrize('ydim', [100, 200])
def test_add_field(xdim, ydim, tmpdir, filename='test_add'):
    filepath = tmpdir.join(filename)
    data, dimensions = generate_fieldset(xdim, ydim)
    fieldset = FieldSet.from_data(data, dimensions)
    field = Field('newfld', fieldset.U.data, lon=fieldset.U.lon, lat=fieldset.U.lat)
    fieldset.add_field(field)
    assert fieldset.newfld.data.shape == fieldset.U.data.shape
    fieldset.write(filepath)


@pytest.mark.parametrize('dupobject', ['same', 'new'])
def test_add_duplicate_field(dupobject):
    data, dimensions = generate_fieldset(100, 100)
    fieldset = FieldSet.from_data(data, dimensions)
    field = Field('newfld', fieldset.U.data, lon=fieldset.U.lon, lat=fieldset.U.lat)
    fieldset.add_field(field)
    error_thrown = False
    try:
        if dupobject == 'same':
            fieldset.add_field(field)
        elif dupobject == 'new':
            field2 = Field('newfld', np.ones((2, 2)), lon=np.array([0, 1]), lat=np.array([0, 2]))
            fieldset.add_field(field2)
    except RuntimeError:
        error_thrown = True

    assert error_thrown


@pytest.mark.parametrize('mesh', ['flat', 'spherical'])
def test_fieldset_celledgesizes(mesh):
    data, dimensions = generate_fieldset(10, 7)
    fieldset = FieldSet.from_data(data, dimensions, mesh=mesh)

    fieldset.U.calc_cell_edge_sizes()
    D_meridional = fieldset.U.cell_edge_sizes['y']
    D_zonal = fieldset.U.cell_edge_sizes['x']

    assert np.allclose(D_meridional.flatten(), D_meridional[0, 0])  # all meridional distances should be the same in either mesh
    if mesh == 'flat':
        assert np.allclose(D_zonal.flatten(), D_zonal[0, 0])  # all zonal distances should be the same in flat mesh
    else:
        assert all((np.gradient(D_zonal, axis=0) < 0).flatten())  # zonal distances should decrease in spherical mesh


@pytest.mark.parametrize('dx, dy', [('e1u', 'e2u'), ('e1v', 'e2v')])
def test_fieldset_celledgesizes_curvilinear(dx, dy):
    fname = path.join(path.dirname(__file__), 'test_data', 'mask_nemo_cross_180lon.nc')
    filenames = {'dx': fname, 'dy': fname, 'mesh_mask': fname}
    variables = {'dx': dx, 'dy': dy}
    dimensions = {'dx': {'lon': 'glamu', 'lat': 'gphiu'},
                  'dy': {'lon': 'glamu', 'lat': 'gphiu'}}
    fieldset = FieldSet.from_nemo(filenames, variables, dimensions)

    # explicitly setting cell_edge_sizes from e1u and e2u etc
    fieldset.dx.grid.cell_edge_sizes['x'] = fieldset.dx.data
    fieldset.dx.grid.cell_edge_sizes['y'] = fieldset.dy.data

    A = fieldset.dx.cell_areas()
    assert np.allclose(A, fieldset.dx.data * fieldset.dy.data)


def test_fieldset_write_curvilinear(tmpdir):
    fname = path.join(path.dirname(__file__), 'test_data', 'mask_nemo_cross_180lon.nc')
    filenames = {'dx': fname, 'mesh_mask': fname}
    variables = {'dx': 'e1u'}
    dimensions = {'lon': 'glamu', 'lat': 'gphiu'}
    fieldset = FieldSet.from_nemo(filenames, variables, dimensions)
    assert fieldset.dx.creation_log == 'from_nemo'

    newfile = tmpdir.join('curv_field')
    fieldset.write(newfile)

    fieldset2 = FieldSet.from_netcdf(filenames=newfile+'dx.nc', variables={'dx': 'dx'}, dimensions={'lon': 'nav_lon', 'lat': 'nav_lat'})
    assert fieldset2.dx.creation_log == 'from_netcdf'

    for var in ['lon', 'lat', 'data']:
        assert np.allclose(getattr(fieldset2.dx, var), getattr(fieldset.dx, var))


@pytest.mark.parametrize('mesh', ['flat', 'spherical'])
def test_fieldset_cellareas(mesh):
    data, dimensions = generate_fieldset(10, 7)
    fieldset = FieldSet.from_data(data, dimensions, mesh=mesh)
    cell_areas = fieldset.V.cell_areas()
    if mesh == 'flat':
        assert np.allclose(cell_areas.flatten(), cell_areas[0, 0], rtol=1e-3)
    else:
        assert all((np.gradient(cell_areas, axis=0) < 0).flatten())  # areas should decrease with latitude in spherical mesh
        for y in range(cell_areas.shape[0]):
            assert np.allclose(cell_areas[y, :], cell_areas[y, 0], rtol=1e-3)


@pytest.mark.parametrize('mesh', ['flat', 'spherical'])
def test_fieldset_gradient(mesh):
    data, dimensions = generate_fieldset(5, 3)
    fieldset = FieldSet.from_data(data, dimensions, mesh=mesh)

    # Calculate field gradients for testing against numpy gradients.
    dFdx, dFdy = fieldset.V.gradient()

    # Create numpy fields.
    conv_factor = 6.371e6 * np.pi / 180. if mesh == 'spherical' else 1.
    np_dFdx = np.gradient(fieldset.V.data[0, :, :], (np.diff(fieldset.V.lon) * conv_factor)[0], axis=1)
    np_dFdy = np.gradient(fieldset.V.data[0, :, :], (np.diff(fieldset.V.lat) * conv_factor)[0], axis=0)
    if mesh == 'spherical':
        for y in range(np_dFdx.shape[0]):
            np_dFdx[:, y] /= math.cos(fieldset.V.grid.lat[y] * math.pi / 180.)

    assert np.allclose(dFdx.data, np_dFdx, rtol=5e-2)  # Field gradient dx.
    assert np.allclose(dFdy.data, np_dFdy, rtol=5e-2)  # Field gradient dy.


def addConst(particle, fieldset, time):
    particle.lon = particle.lon + fieldset.movewest + fieldset.moveeast


@pytest.mark.parametrize('mode', ['scipy', 'jit'])
def test_fieldset_constant(mode):
    data, dimensions = generate_fieldset(100, 100)
    fieldset = FieldSet.from_data(data, dimensions)
    westval = -0.2
    eastval = 0.3
    fieldset.add_constant('movewest', westval)
    fieldset.add_constant('moveeast', eastval)
    assert fieldset.movewest == westval

    pset = ParticleSet.from_line(fieldset, size=1, pclass=ptype[mode],
                                 start=(0.5, 0.5), finish=(0.5, 0.5))
    pset.execute(pset.Kernel(addConst), dt=1, runtime=1)
    assert abs(pset[0].lon - (0.5 + westval + eastval)) < 1e-4


@pytest.mark.parametrize('mode', ['scipy', 'jit'])
@pytest.mark.parametrize('swapUV', [False, True])
def test_vector_fields(mode, swapUV):
    lon = np.linspace(0., 10., 12, dtype=np.float32)
    lat = np.linspace(0., 10., 10, dtype=np.float32)
    U = np.ones((10, 12), dtype=np.float32)
    V = np.zeros((10, 12), dtype=np.float32)
    data = {'U': U, 'V': V}
    dimensions = {'U': {'lat': lat, 'lon': lon},
                  'V': {'lat': lat, 'lon': lon}}
    fieldset = FieldSet.from_data(data, dimensions, mesh='flat')
    if swapUV:  # we test that we can freely edit whatever UV field
        UV = VectorField('UV', fieldset.V, fieldset.U)
        fieldset.add_vector_field(UV)

    pset = ParticleSet.from_line(fieldset, size=1, pclass=ptype[mode],
                                 start=(0.5, 0.5), finish=(0.5, 0.5))
    pset.execute(AdvectionRK4, dt=1, runtime=1)
    if swapUV:
        assert abs(pset[0].lon - .5) < 1e-9
        assert abs(pset[0].lat - 1.5) < 1e-9
    else:
        assert abs(pset[0].lon - 1.5) < 1e-9
        assert abs(pset[0].lat - .5) < 1e-9


@pytest.mark.parametrize('mode', ['scipy', 'jit'])
@pytest.mark.parametrize('time_periodic', [True, False])
@pytest.mark.parametrize('dt_sign', [-1, 1])
def test_periodic(mode, time_periodic, dt_sign):
    lon = np.array([0, 1], dtype=np.float32)
    lat = np.array([0, 1], dtype=np.float32)
    depth = np.array([0, 1], dtype=np.float32)
    tsize = 24*60+1
    period = 86400
    time = np.linspace(0, period, tsize, dtype=np.float64)

    def temp_func(time):
        return 20 + 2 * np.sin(time*2*np.pi/period)
    temp_vec = temp_func(time)

    U = np.zeros((2, 2, 2, tsize), dtype=np.float32)
    V = np.zeros((2, 2, 2, tsize), dtype=np.float32)
    V[0, 0, 0, :] = 1e-5
    W = np.zeros((2, 2, 2, tsize), dtype=np.float32)
    temp = np.zeros((2, 2, 2, tsize), dtype=np.float32)
    temp[:, :, :, :] = temp_vec

    data = {'U': U, 'V': V, 'W': W, 'temp': temp}
    dimensions = {'lon': lon, 'lat': lat, 'depth': depth, 'time': time}
    fieldset = FieldSet.from_data(data, dimensions, mesh='flat', time_periodic=time_periodic, transpose=True, allow_time_extrapolation=True)

    def sampleTemp(particle, fieldset, time):
        # Note that fieldset.temp is interpolated at time=time+dt.
        # Indeed, sampleTemp is called at time=time, but the result is written
        # at time=time+dt, after the Kernel update
        particle.temp = fieldset.temp[time+particle.dt, particle.depth, particle.lat, particle.lon]
        # test if we can interpolate UV and UVW together
        (particle.u1, particle.v1) = fieldset.UV[time+particle.dt, particle.depth, particle.lat, particle.lon]
        (particle.u2, particle.v2, w_) = fieldset.UVW[time+particle.dt, particle.depth, particle.lat, particle.lon]

    class MyParticle(ptype[mode]):
        temp = Variable('temp', dtype=np.float32, initial=20.)
        u1 = Variable('u1', dtype=np.float32, initial=0.)
        u2 = Variable('u2', dtype=np.float32, initial=0.)
        v1 = Variable('v1', dtype=np.float32, initial=0.)
        v2 = Variable('v2', dtype=np.float32, initial=0.)

    dt_sign = -1
    pset = ParticleSet.from_list(fieldset, pclass=MyParticle,
                                 lon=[0.5], lat=[0.5], depth=[0.5])
    pset.execute(AdvectionRK4_3D + pset.Kernel(sampleTemp),
                 runtime=delta(hours=51), dt=delta(hours=dt_sign*1))

    if time_periodic:
        t = pset.particles[0].time
        temp_theo = temp_func(t)
    elif dt_sign == 1:
        temp_theo = temp_vec[-1]
    elif dt_sign == -1:
        temp_theo = temp_vec[0]
    assert np.allclose(temp_theo, pset.particles[0].temp, atol=1e-5)
    assert np.allclose(pset.particles[0].u1, pset.particles[0].u2)
    assert np.allclose(pset.particles[0].v1, pset.particles[0].v2)


@pytest.mark.parametrize('fail', [False, pytest.param(True, marks=pytest.mark.xfail(strict=True))])
def test_fieldset_defer_loading_with_diff_time_origin(tmpdir, fail, filename='test_parcels_defer_loading'):
    filepath = tmpdir.join(filename)
    data0, dims0 = generate_fieldset(10, 10, 1, 10)
    dims0['time'] = np.arange(0, 10, 1) * 3600
    fieldset_out = FieldSet.from_data(data0, dims0)
    fieldset_out.U.grid.time_origin = TimeConverter(np.datetime64('2018-04-20'))
    fieldset_out.V.grid.time_origin = TimeConverter(np.datetime64('2018-04-20'))
    data1, dims1 = generate_fieldset(10, 10, 1, 10)
    if fail:
        dims1['time'] = np.arange(0, 10, 1) * 3600
    else:
        dims1['time'] = np.arange(0, 10, 1) * 1800 + (24+25)*3600
    if fail:
        Wtime_origin = TimeConverter(np.datetime64('2018-04-22'))
    else:
        Wtime_origin = TimeConverter(np.datetime64('2018-04-18'))
    gridW = RectilinearZGrid(dims1['lon'], dims1['lat'], dims1['depth'], dims1['time'], time_origin=Wtime_origin)
    fieldW = Field('W', np.zeros(data1['U'].shape), grid=gridW)
    fieldset_out.add_field(fieldW)
    fieldset_out.write(filepath)
    fieldset = FieldSet.from_parcels(filepath, extra_fields={'W': 'W'})
    assert fieldset.U.creation_log == 'from_parcels'
    pset = ParticleSet.from_list(fieldset, pclass=JITParticle, lon=[0.5], lat=[0.5], depth=[0.5],
                                 time=[datetime.datetime(2018, 4, 20, 1)])
    pset.execute(AdvectionRK4_3D, runtime=delta(hours=4), dt=delta(hours=1))


@pytest.mark.parametrize('zdim', [2, 8])
@pytest.mark.parametrize('scale_fac', [0.2, 4, 1])
def test_fieldset_defer_loading_function(zdim, scale_fac, tmpdir, filename='test_parcels_defer_loading'):
    filepath = tmpdir.join(filename)
    data0, dims0 = generate_fieldset(3, 3, zdim, 10)
    data0['U'][:, 0, :, :] = np.nan  # setting first layer to nan, which will be changed to zero (and all other layers to 1)
    dims0['time'] = np.arange(0, 10, 1) * 3600
    dims0['depth'] = np.arange(0, zdim, 1)
    fieldset_out = FieldSet.from_data(data0, dims0)
    fieldset_out.write(filepath)
    fieldset = FieldSet.from_parcels(filepath)

    # testing for combination of deferred-loaded and numpy Fields
    fieldset.add_field(Field('numpyfield', np.zeros((10, zdim, 3, 3)), grid=fieldset.U.grid))

    # testing for scaling factors
    fieldset.U.set_scaling_factor(scale_fac)

    dFdx, dFdy = fieldset.V.gradient()

    dz = np.gradient(fieldset.U.depth)
    DZ = np.moveaxis(np.tile(dz, (fieldset.U.grid.ydim, fieldset.U.grid.xdim, 1)), [0, 1, 2], [1, 2, 0])

    def compute(fieldset):
        # Calculating vertical weighted average
        for f in [fieldset.U, fieldset.V]:
            for tind in f.loaded_time_indices:
                f.data[tind, :] = np.sum(f.data[tind, :] * DZ, axis=0) / sum(dz)

    fieldset.compute_on_defer = compute
    fieldset.computeTimeChunk(1, 1)
    assert np.allclose(fieldset.U.data, scale_fac*(zdim-1.)/zdim)
    assert np.allclose(dFdx.data, 0)

    pset = ParticleSet(fieldset, JITParticle, 0, 0)

    def DoNothing(particle, fieldset, time):
        return ErrorCode.Success

    pset.execute(DoNothing, dt=3600)
    assert np.allclose(fieldset.U.data, scale_fac*(zdim-1.)/zdim)
    assert np.allclose(dFdx.data, 0)


@pytest.mark.parametrize('maxlatind', [3, pytest.param(2, marks=pytest.mark.xfail(strict=True))])
def test_fieldset_from_xarray(maxlatind):
    def generate_dataset(xdim, ydim, zdim=1, tdim=1):
        lon = np.linspace(0., 12, xdim, dtype=np.float32)
        lat = np.linspace(0., 12, ydim, dtype=np.float32)
        depth = np.linspace(0., 20., zdim, dtype=np.float32)
        time = np.linspace(0., 10, tdim, dtype=np.float64)
        Uxr = np.ones((tdim, zdim, ydim, xdim), dtype=np.float32)
        Vxr = np.ones((tdim, zdim, ydim, xdim), dtype=np.float32)
        for t in range(Uxr.shape[0]):
            Uxr[t, :, :, :] = t/10.
        coords = {'lat': lat, 'lon': lon, 'depth': depth, 'time': time}
        dims = ('time', 'depth', 'lat', 'lon')
        return xr.Dataset({'Uxr': xr.DataArray(Uxr, coords=coords, dims=dims),
                           'Vxr': xr.DataArray(Vxr, coords=coords, dims=dims)})

    ds = generate_dataset(3, 3, 2, 10)
    variables = {'U': 'Uxr', 'V': 'Vxr'}
    dimensions = {'lat': 'lat', 'lon': 'lon', 'depth': 'depth', 'time': 'time'}
    indices = {'lat': range(0, maxlatind)}
    fieldset = FieldSet.from_xarray_dataset(ds, variables, dimensions, indices, mesh='flat')
    assert fieldset.U.creation_log == 'from_xarray_dataset'

    pset = ParticleSet(fieldset, JITParticle, 0, 0)

    pset.execute(AdvectionRK4, dt=1)
    assert np.allclose(pset[0].lon, 4.5) and np.allclose(pset[0].lat, 10)


def test_fieldset_from_data_gridtypes(xdim=20, ydim=10, zdim=4):
    """ Simple test for fieldset initialisation from data. """
    lon = np.linspace(0., 10., xdim, dtype=np.float32)
    lat = np.linspace(0., 10., ydim, dtype=np.float32)
    depth = np.linspace(0., 1., zdim, dtype=np.float32)
    depth_s = np.ones((zdim, ydim, xdim))
    U = np.ones((zdim, ydim, xdim))
    V = np.ones((zdim, ydim, xdim))
    dimensions = {'lat': lat, 'lon': lon, 'depth': depth}
    data = {'U': np.array(U, dtype=np.float32), 'V': np.array(V, dtype=np.float32)}
    lonm, latm = np.meshgrid(lon, lat)
    for k in range(zdim):
        data['U'][k, :, :] = lonm * (depth[k]+1) + .1
        depth_s[k, :, :] = depth[k]

    # Rectilinear Z grid
    fieldset = FieldSet.from_data(data, dimensions, mesh='flat')
    pset = ParticleSet(fieldset, ScipyParticle, [0, 0], [0, 0], [0, .4])
    pset.execute(AdvectionRK4, runtime=1, dt=.5)
    plon = [p.lon for p in pset]
    plat = [p.lat for p in pset]
    # sol of  dx/dt = (init_depth+1)*x+0.1; x(0)=0
    assert np.allclose(plon, [0.17173462592827032, 0.2177736932123214])
    assert np.allclose(plat, [1, 1])

    # Rectilinear S grid
    dimensions['depth'] = depth_s
    fieldset = FieldSet.from_data(data, dimensions, mesh='flat')
    pset = ParticleSet(fieldset, ScipyParticle, [0, 0], [0, 0], [0, .4])
    pset.execute(AdvectionRK4, runtime=1, dt=.5)
    assert np.allclose(plon, [p.lon for p in pset])
    assert np.allclose(plat, [p.lat for p in pset])

    # Curvilinear Z grid
    dimensions['lon'] = lonm
    dimensions['lat'] = latm
    dimensions['depth'] = depth
    fieldset = FieldSet.from_data(data, dimensions, mesh='flat')
    pset = ParticleSet(fieldset, ScipyParticle, [0, 0], [0, 0], [0, .4])
    pset.execute(AdvectionRK4, runtime=1, dt=.5)
    assert np.allclose(plon, [p.lon for p in pset])
    assert np.allclose(plat, [p.lat for p in pset])

    # Curvilinear S grid
    dimensions['depth'] = depth_s
    fieldset = FieldSet.from_data(data, dimensions, mesh='flat')
    pset = ParticleSet(fieldset, ScipyParticle, [0, 0], [0, 0], [0, .4])
    pset.execute(AdvectionRK4, runtime=1, dt=.5)
    assert np.allclose(plon, [p.lon for p in pset])
    assert np.allclose(plat, [p.lat for p in pset])
