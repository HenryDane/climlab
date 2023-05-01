import climlab
from climlab import constants as const
import numpy as np
import matplotlib.pyplot as plt
import datetime

model_dims = {'lat':20, 'lon':10, 'lev':30}
n_steps    = 1000

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

# create "blank" model
model = climlab.TimeDependentProcess(state=state)

#  Radiation 
rad = climlab.radiation.RRTMG(state=state,
                              albedo=0.3,
                              timestep=const.seconds_per_day)

# Convection scheme 
conv = climlab.convection.ConvectiveAdjustment(name='Convection',
                                               state=state,
                                               adj_lapse_rate=6.5, # this is the key parameter! We'll discuss below
                                               timestep=rad.timestep)

sun = climlab.radiation.DailyInsolation(name='Insolation', domains=state['Ts'].domain) 

# attach processes
model.add_subprocess('Radiation', rad)
#model.add_subprocess('Convection', conv)
model.add_subprocess('Insolation', sun)

# debug model info
print(model)

model.step_forward()    

start = datetime.datetime.now()
for i in range(n_steps):
    model.step_forward()
    print(i, end='\r')
end = datetime.datetime.now()    
print('Integrated {:d} steps in {:.3f} seconds'.format(n_steps, (end-start).total_seconds()))
print('Days elapsed: ', model.time['days_elapsed'])

def lon_temp():
    plt.plot(model.lon, np.squeeze(model.Ts)[5,:], label='Surface')
    plt.plot(model.lon, model.Tatm[5,:,-1], label='Atmosphere')
    plt.gca().ticklabel_format(useOffset=False)
    plt.xlabel('Longitude')
    plt.ylabel('Temperature (K)')
    plt.legend()
    plt.show()

def lat_temp():
    plt.plot(model.lat, np.squeeze(model.Ts)[:,5], label='Surface')
    plt.plot(model.lat, model.Tatm[:,5,-1], label='Atmosphere')
    plt.gca().ticklabel_format(useOffset=False)
    plt.xlabel('Latitude')
    plt.ylabel('Temperature (K)')
    plt.legend()
    plt.show()

def lev_temp():
    plt.plot(model.Tatm[10, 5, :], model.lev)
    plt.gca().invert_yaxis()
    plt.ylabel('Pressure')
    plt.xlabel('Temperature (K)')
    plt.show()

def map_temp():
    plt.imshow(model.Ts, aspect='auto', 
               # might be wrong bounds?
               extent=(np.min(model.lon), np.max(model.lon), np.min(model.lat), np.max(model.lat)))
    plt.colorbar(label='Temperature (K)')
    plt.title('Surface Temperature')
    plt.show()
    
lon_temp()
lat_temp()
lev_temp()
map_temp()


