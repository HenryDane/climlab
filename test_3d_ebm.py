import climlab
import numpy as np
import matplotlib.pyplot as plt

model_dims = {'lat':20, 'lon':10, 'lev':30}

state = climlab.volume_state(num_lat=model_dims['lat'], num_lon=model_dims['lon'], num_lev=model_dims['lev'])

ebm = climlab.EBM(state=state)

ebm.integrate_years(1)

