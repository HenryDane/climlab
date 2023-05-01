def make_volume_domains(num_lev=30, num_lat=30, num_lon=30, water_depth=1.0, **kwargs):
    '''
    Create 3D volumes
    '''
    
    # make axes
    latax = climlab.domain.axis.Axis(axis_type='lat', num_points=num_lat)
    lonax = climlab.domain.axis.Axis(axis_type='lon', num_points=num_lon)
    levax = climlab.domain.axis.Axis(axis_type='lev', num_points=num_lev)
    depax = climlab.domain.axis.Axis(axis_type='depth', bounds=[water_depth, 0.0])
    
    # make ocean
    slab = climlab.domain.domain.SlabOcean(axes={'lat':latax, 'lon':lonax, 'depth':depax}, **kwargs)
    
    # make atmosphere
    atm = climlab.domain.domain.Atmosphere(axes={'lat':latax, 'lon':lonax, 'lev':levax}, **kwargs)
    
    # return the domains
    return slab, atm
    
def volume_state(num_lev=30, num_lat=30, num_lon=30, water_depth=1.0):
    '''
    Create a 3D climlab state
    '''
    
    # make domains
    sfc, atm = make_volume_domains(num_lev=num_lev, num_lat=num_lat, num_lon=num_lon, water_depth=water_depth)
    
    # update elevation value
    num_lev = atm.lev.num_points
    
    # make fields
    Ts = climlab.domain.field.Field(288.*np.ones(sfc.shape), domain=sfc)
    Tinitial = np.tile(np.linspace(200., 288.-10., num_lev), sfc.shape)
    Tatm = climlab.domain.field.Field(Tinitial, domain=atm)
    
    # create state
    state = climlab.utils.attrdict.AttrDict()
    state['Ts'] = Ts
    state['Tatm'] = Tatm
    
    # return the state
    return state    
    
# make state
state = volume_state(num_lat=model_dims['lat'], num_lon=model_dims['lon'], num_lev=model_dims['lev'])

sun = climlab.radiation.DailyInsolation(name='Insolation', domains=state['Ts'].domain) 


