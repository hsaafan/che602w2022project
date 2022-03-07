import numpy as np
import pint

from typing import Tuple
from math import exp

from parameters import *


def air_scouring(v_sg: pint.Quantity, X_TSS: pint.Quantity,
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


def resistance_change(J: pint.Quantity, m_rback: pint.Quantity,
                      S_UAP: pint.Quantity, S_BAP: pint.Quantity,
                      X_MLSS: pint.Quantity, X_TSS: pint.Quantity,
                      alpha_c: pint.Quantity
                      ) -> Tuple[pint.Quantity, pint.Quantity]:
    a = mu * alpha_c * X_TSS  # (Janus 2013, p.189)
    R_dot_i = a * k_i * exp(b * J) * J * (S_UAP + S_BAP)
    R_dot_r = alpha_c * (J * X_MLSS - m_rback)
    return(R_dot_i, R_dot_r)


def membrane_resistance(J: pint.Quantity, v_sg: pint.Quantity,
                        S_UAP: pint.Quantity, S_BAP: pint.Quantity,
                        X_EPS: pint.Quantity, X_MLSS: pint.Quantity,
                        X_TSS: pint.Quantity, T_l: pint.Quantity,
                        R_i: pint.Quantity, R_r: pint.Quantity,
                        alpha_c: pint.Quantity, m_rback: pint.Quantity
                        ) -> Tuple[pint.Quantity, pint.Quantity,
                                   pint.Quantity, pint.Quantity]:

    R_dot_i, R_dot_r = resistance_change(J, m_rback, S_UAP, S_BAP, X_MLSS,
                                         X_TSS, alpha_c)
    alpha_c0 = specific_cake_resistance(X_EPS, X_MLSS)
    tau_w = air_scouring(v_sg, X_TSS, T_l)

    R_i += R_dot_i
    R_r += R_dot_r
    R_t = R_m + R_r + R_i

    Delta_P = J / (mu * R_t)
    alpha_c = alpha_c0 * (Delta_P / Delta_P_crit) ** 2
    m_rback = k_r * (tau_w - lambda_m * Delta_P)
    return(R_i, R_r, alpha_c, m_rback)


def specific_cake_resistance(X_EPS: pint.Quantity, X_MLSS: pint.Quantity
                             ) -> pint.Quantity:
    alpha_c0 = 1.966e9 * (X_EPS/X_MLSS) - 2.564e13
    return(alpha_c0)


def membrane_model() -> dict:
    pass
