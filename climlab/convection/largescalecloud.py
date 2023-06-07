'''
A climlab process for the Simcloud 1.0 routine
'''

from __future__ import absolute_import

import numpy as np
import warnings
from climlab.process import TimeDependentProcess
from climlab.utils.thermo import qsat
from climlab import constants as const

class Simcloud(TimeDependendProcess):
    def __init__(self, **kwargs):
        super(Simcloud, self).__init__(**kwargs)
        
        # is this correct?
        self.add_input('relative_humidity', 0 * self.Tatm)
        
    def _compute(self):
        # todo -- update this???
        tendencies = {}
        for name, value in self.state.items():
            tendencies[name] = value * 0.0
           return tendencies
        
