import pint
import os.path

import membrane
import bioreactor
from parameters import ureg
import random


class MBRModel:
    def __init__(self, state: dict) -> None:
        """ Constructor function for integrated model

        Parameters
        ----------
        state: dict
            See state.py file for all required state variables
        """
        self.state = state
        self.membrane = membrane.Membrane()
        self.bioreactor = bioreactor.Bioreactor()

    def record_state(self, file_path: str = 'data.csv') -> None:
        """Record state variables"""
        if not os.path.isfile(file_path):
            # If file doesn't exist create a new csv file with headers
            self.header_order = list(self.state.keys())
            with open(file_path, 'a+') as f:
                for k in self.header_order:
                    if isinstance(self.state[k], pint.Quantity):
                        units = self.state[k].to_base_units().units
                    else:
                        units = '-'
                    f.write(f'{k} [{units}], ')
                f.write('\n')

        with open(file_path, 'a+') as f:
            for k in self.header_order:
                if isinstance(self.state[k], pint.Quantity):
                    value = self.state[k].to_base_units().magnitude
                else:
                    value = self.state[k]
                f.write(f'{value}, ')
            f.write('\n')

    def vary_flowrate(self, pcnt_change: float) -> pint.Quantity:
        Q = self.state['Q_in']
        Q_max = self.state['Q_max']
        Q_min = self.state['Q_min']

        Q_rng = Q_max - Q_min
        pcnt_change = (2 * pcnt_change) * random.random() - pcnt_change
        Q_new = Q + Q_rng * pcnt_change / 100
        if Q_new > Q_max:
            Q_new = Q_max
        elif Q_new < Q_min:
            Q_new = Q_min
        return(Q_new)

    def step_model(self, t_step: pint.Quantity) -> None:
        """Step the integrated model forward

        Parameters
        ----------
        t_step: pint.Quantity
            The time step
        """
        # self.state['Q_in'] = self.vary_flowrate(1)
        self.state = self.bioreactor.step(t_step, self.state)
        self.state = self.membrane.step(t_step, self.state)
        self.state['time'] = self.state['time'] + t_step
        # self.record_state()
        return(self.state)
