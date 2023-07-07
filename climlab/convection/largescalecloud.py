'''
A process for modelling cloud fraction inspired by simcloud 1.0
'''

from __future__ import absolute_import

import numpy as np
import warnings
from climlab.process import TimeDependentProcess
from climlab.utils.thermo import qsat
from climlab import constants as const

class LargeScaleCloud(TimeDependentProcess):
    def __init__(self, **kwargs):
        super(LargeScaleCloud, self).__init__(**kwargs)
        
        # add diagnostics
        self.add_diagnostic('cldfrac', self.Tatm * 0)
        self.add_diagnostic('r_liq', self.Tatm * 0)
        self.add_diagnostic('r_ice', self.Tatm * 0)
        self.add_diagnostic('r_eff', self.Tatm * 0)
        self.add_diagnostic('clwmr', self.Tatm * 0)
        self.add_diagnostic('clwp', self.Tatm * 0)
        self.add_diagnostic('f_liq', self.Tatm * 0)
        
        # 0.95 @ 1000hPa, 0.85 @ 700hPa, 0.99 @ 200hPa
        self.RH_crit = np.interp(self.lev, [200, 700, 1000], [0.99, 0.85, 0.95])
        
    def _compute(self):
        print(type(self.lev))
        # NOTE: RH is computed as 0-100 by EC but 0-1 by other stuff
        # TODO: need to figure out how to access diagnostics computed by other processes
        if 'relative_humidity' in self.diagnostics:
            RH = np.minimum(1, np.maximum(0, self.diagnostics['relative_humidity'] / 100))
        else:
            RH = np.array([0.0])
        # compute 'a' parameter
        a    = calc_a(self.lev, self.lev.max())
        # compute threshold
        qv   = calc_qv(self.lev)
        # compute freeze-dry adjustment 
        fdc  = calc_fdc(self.state['q'], qv)
        # calc large-scale cloud fraction
        clf  = calc_cs(RH, a) * fdc
        # calculate in-cloud liquid water mixing ratio
        wl   = calc_wl(self.Tatm)
        # calculate cloud liquid water path
        clwp = calc_clwp(self.lev, clf, wl)
        # calculate liquid fraction
        fl   = calc_fl(self.Tatm)
        # calculate effective radius
        reff = calc_re(fl * clf)
        
        # apply diagnostics
        self.cldfrac = clf
        self.r_liq   = 14.0 * fl
        self.r_ice   = 25.0 * (1 - fl)
        self.r_eff   = reff
        self.clwmr   = wl
        self.clwp    = clwp
        self.f_liq   = fl
        
        print(RH.min(), RH.max(), RH.mean())
        
        return {}
        
def calc_a(p, ps, at=13.0, as_=36.0, n=12.0):
    # equation (2)
    return at + (as_ - at) * np.exp(1 - (ps / p)**n)    
    
def calc_cs(RH, a):
    # equation (1)
    return np.minimum(1, np.maximum(0, a * (RH - 1) + 1))
    
def calc_qv(p, q0=6e-3, n=2.5):
    # equation (6)
    return q0 * (p / np.max(p))**n     
    
def calc_fdc(q, qv):
    # last term of RHS of equation (5)
    # q is specific humidity in kg/kg
    return np.maximum(0.15, np.minimum(1.0, q / qv)) 
    
def calc_clwp(p, Cs, wl, g=9.81):
    # equation (13)
    # [0-1] * (g/kg) * [kg/m/s*2] * [s^2/m]
    return np.cumsum(Cs * wl * np.concatenate([[0], np.diff(p)]) / g) * 1e-5         
    
def calc_wl(T, wl0=0.18, w0=3e-4, Tl0=280, T0=220):
    # equation (12)
    return np.maximum(w0, wl0 * np.minimum(1, (T - T0) / (Tl0 - T0)))     

def calc_fl(T, Cs=1.0, Tmin=233.15, Tmax=268.15):
    # equation (10)
    # adjusted to require cloud fraction by me
    return np.maximum(0, np.minimum(1, (T - Tmin) / (Tmax - Tmin))) * Cs 

def calc_re(fl):
    # equation (11)
    return 14.0 * fl + 25.0 * (1.0 - fl)                   
