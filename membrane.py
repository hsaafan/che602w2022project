import pint

from typing import Tuple
from math import exp

from parameters import *


class Membrane:
    def air_scouring(self, v_sg: pint.Quantity, X_TSS: pint.Quantity,
                     T_l: pint.Quantity) -> pint.Quantity:
        # Convert variables for empirical model
        v_sg = v_sg.to(ureg.cm / ureg.s).magnitude
        X_TSS = X_TSS.to(ureg.kg / ureg.m ** 3).magnitude
        T_l = T_l.to(ureg.degC).magnitude

        # Calculate model parameters
        def p(a1, a2, a3, a4, a5):
            p = a1 + a2 * X_TSS + a3 * T_l + a4 * X_TSS ** 2 + a5 * X_TSS * T_l
            return(p)
        p1 = p(-9.884e-3, -1.106e-4, 1.256e-5, 1.669e-6, -3.722e-7)
        p2 = p(4.231e-2, 3.862e-4, -9.708e-5, 3.378e-6, 4.288e-6)
        p3 = p(0.2627, 6.695e-3, -5.703e-4, -3.598e-5, -5.445e-5)
        p4 = p(-0.151, -2.212e-3, -4.014e-4, 1.985e-4, 8.685e-7)

        # Calculate shear stress
        tau_w = (p1 * v_sg ** 3 + p2 * v_sg ** 2 + p3 * v_sg + p4) * ureg.Pa
        return(tau_w)

    def resistance_change(self, J: pint.Quantity, m_rback: pint.Quantity,
                          S_UAP: pint.Quantity, S_BAP: pint.Quantity,
                          X_MLSS: pint.Quantity, X_TSS: pint.Quantity,
                          alpha_c: pint.Quantity
                          ) -> Tuple[pint.Quantity, pint.Quantity]:
        R_dot_i = a * k_i * exp(b * J) * J * (S_UAP + S_BAP)
        R_dot_r = alpha_c * (J * X_MLSS - m_rback)
        return(R_dot_i, R_dot_r)

    def membrane_resistance(self, t_step: pint.Quantity, state: dict,
                            J: pint.Quantity) -> dict:
        new_state = state.copy()
        # Interface calculations
        X_EPS = state['X_EPS']
        X_MLSS = state['X_MLSS']
        X_TSS = state['X_TSS']
        T_l = state['temperature']
        v_sg = state['v_sg']
        # alpha_c0 = self.specific_cake_resistance(X_EPS, X_MLSS)
        new_state['tau_w'] = self.air_scouring(v_sg, X_TSS, T_l)

        # Update resistances
        S_UAP = state['S_UAP']
        S_BAP = state['S_BAP']
        alpha_c = state['alpha_c']
        m_rback = state['m_rback']
        R_dot_i, R_dot_r = self.resistance_change(J, m_rback, S_UAP, S_BAP,
                                                  X_MLSS, X_TSS, alpha_c)

        new_state['R_i'] = state['R_i'] + R_dot_i * t_step
        new_state['R_r'] = state['R_r'] + R_dot_r * t_step
        new_state['R_t'] = R_m + new_state['R_i'] + new_state['R_r']

        Delta_P = J * (mu * new_state['R_t'])
        new_state['TMP'] = Delta_P
        # Model taken from (Janus, 2013) with parameters from p.193
        new_state['m_rback'] = back_transport_coefficient * state['X_MLSS']
        # new_state['alpha_c'] = alpha_c0 * (Delta_P / Delta_P_crit) ** 2
        new_state['alpha_c'] = 1.12 * (ureg.m / ureg.kg)  # (Janus, 2013 p.280)
        return(new_state)

    def specific_cake_resistance(self, X_EPS: pint.Quantity,
                                 X_MLSS: pint.Quantity) -> pint.Quantity:
        X_EPS = X_EPS.to_base_units().magnitude
        X_MLSS = X_MLSS.to_base_units().magnitude
        alpha_c0 = (1.966e15 * (X_EPS/X_MLSS) - 2.564e13) * (ureg.m / ureg.kg)
        return(alpha_c0)

    def step(self, t_step: pint.Quantity, state: dict) -> dict:
        """Step the bioreactor model

        Parameters
        ----------
        t_step: pint.Quantity
            The time step
        state: dict
            Contains all state variables
        Output
        ------
        new_state: dict
            Contains all updated state variables
        """
        J = state['Q_out'] / (membrane_density * state['volume'])
        new_state = self.membrane_resistance(t_step, state, J)
        return(new_state)
