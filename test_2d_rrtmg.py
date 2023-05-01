import climlab
alb = 0.25
#  State variables (Air and surface temperature)
state = climlab.column_state(num_lev=30)
#  Fixed relative humidity
h2o = climlab.radiation.ManabeWaterVapor(name='WaterVapor', state=state)
#  Couple water vapor to radiation
rad = climlab.radiation.RRTMG(name='Radiation', state=state, specific_humidity=h2o.q, albedo=alb)
#  Convective adjustment
conv = climlab.convection.ConvectiveAdjustment(name='Convection', state=state, adj_lapse_rate=6.5)
#  Couple everything together
rcm = climlab.couple([rad,h2o,conv], name='Radiative-Convective Model')
#  Run the model
rcm.integrate_years(1)
#  Check for energy balance
print(rcm.ASR - rcm.OLR)
