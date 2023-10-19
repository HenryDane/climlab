'''
A process for modelling cloud fraction inspired by simcloud 1.0
'''

from __future__ import absolute_import

import numpy as np
import warnings
from climlab.process import TimeDependentProcess
from climlab.utils.thermo import qsat, lifting_condensation_level
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
        self.add_diagnostic('ciwp', self.Tatm * 0)
        self.add_diagnostic('f_liq', self.Tatm * 0)
        
        # 0.95 @ 1000hPa, 0.85 @ 700hPa, 0.99 @ 200hPa
        self.RH_crit = np.interp(self.lev, [200, 700, 1000], [0.99, 0.85, 0.95])
        
        # calculate altitude table
        Mair = 28.97; g=9.81; R=8.314; H = (R*250) / (Mair * g); LR = 9.81
        self.alt = np.log(self.lev / self.lev.max()) * H * -1
        
    def _compute(self):
        # NOTE: RH is computed as 0-100 by EC but 0-1 by other stuff
        if 'relative_humidity' in self.diagnostics:
            RH = np.minimum(1, np.maximum(0, self.diagnostics['relative_humidity'] / 100))
        else:
            RH = self.Tatm * 0
        
        # figure out inversion level altitude
        theta = self.Tatm * (self.lev / self.lev.max()) ** (287.052874 / 1004.0)
        dthetadP = np.gradient(theta[:,self.lev > 750], self.lev[self.lev > 750], axis=1) # K/hPa
        mingrad = np.min(dthetadP, axis=1)
        zinvidx = np.argmin(dthetadP, axis=1)
        zinv = self.alt[zinvidx.data] * 1e3 # in meters
        
        # find zlcl
        zlcl = lifting_condensation_level(self.Ts, RH[-1]) # this throws RuntimeWarning (div by 0 in log(RH))
        zlcl = zlcl[:,0] # we only care about LCL computed at the surface
        
        # calculate cloud fraction `C`
        a   = calc_a(self.lev, self.lev.max())
        qv  = calc_qv(self.lev)
        fd  = calc_f(self.q, qv)
        elf = calc_elf(fd, zinv, zlcl)
        #print(fd, zinv, zlcl, elf)
        C   = np.maximum(calc_cs(RH, a) * fd, calc_csc(elf) * (mingrad < -0.08)[:,None])
        
        # calculate cloud properties
        fl  = calc_fl(self.Tatm)
        re  = calc_re(fl * C)
        wl  = calc_wl(self.Tatm)
        cwp = calc_clwp(self.lev, C, wl)
        
        # assign diagnostics
        #'''
        self.cldfrac = C
        self.f_liq   = fl
        self.r_eff   = re
        self.r_liq   = 14.0 * fl
        self.r_ice   = 25.0 * (1 - fl)
        self.clwp    = cwp * fl
        self.ciwp    = cwp * (1 - fl)
        self.clwmr   = wl
        #'''
        
        return {}

def calc_cs(RH, a):
    # equation (1)
    return np.minimum(1, np.maximum(0, a * (RH - 1) + 1))
    
def calc_a(p, ps, at=13.0, as_=36.0, n=12.0):
    # equation (2)
    return at + (as_ - at) * np.exp(1 - (ps / p)**n)
    
def calc_cs_alt(RH, RH_c):
    # equation (3)
    return np.maximum(0, 1 - np.sqrt((1 - RH)/(1 - RH_c)))    
    
def calc_f(q, qv):
    # last term of RHS of equation (5)
    # q is specific humidity in kg/kg
    return np.maximum(0.15, np.minimum(1.0, q / qv))    
    
def calc_qv(p, q0=6e-3, n=2.5):
    # equation (6)
    return q0 * (p / np.max(p))**n    
    
def calc_elf(f, z_inv, z_LCL, delta_zs=2750.0):
    # equation (7)
    return f * (1 - (np.sqrt(z_inv * z_LCL) / delta_zs))    
    
def calc_csc(elf, b=1.3, c=0.1):
    return np.minimum(1, np.maximum(0, b * elf + c))   
    
def calc_fl(T, Cs=1.0, Tmin=233.15, Tmax=268.15):
    # equation (10)
    # adjusted to require cloud fraction 
    return np.maximum(0, np.minimum(1, (T - Tmin) / (Tmax - Tmin))) * Cs     
    
def calc_re(fl):
    # equation (11)
    return 14.0 * fl + 25.0 * (1.0 - fl)
    
def calc_wl(T, wl0=0.18, w0=3e-4, Tl0=280, T0=220):
    # equation (12)
    return np.maximum(w0, wl0 * np.minimum(1, (T - T0) / (Tl0 - T0)))    
    
def calc_clwp(p, Cs, wl, g=9.81):
    # equation (13)
    # [0-1] * (g/kg) * [kg/m/s*2] * [s^2/m]
    return np.cumsum(Cs * wl * np.concatenate([[0], np.diff(p)]) / g, axis=0) * 1e-5    

