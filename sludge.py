"""The activated sludge model

All process methods (p1, p2a, p2b, p2c, p3a, p3b, p3c, p4, p5, p6, p7, p8, p9)
will have inputs in terms of concentrations or temperature in pint.Quantity
objects and will output rates in units of ([mass] / [volume] / [time]).

Methods
-------
build_petersen_matrix -> np.ndarray
    Using parameters from the parameters.py file this method will create the
    Petersen matrix.
p1 -> pint.Quantity
    Ammonification
p2a -> pint.Quantity
    Aerobic growth on S_S
p2b -> pint.Quantity
    Aerobic growth on S_BAP
p2c -> pint.Quantity
    Aerobic growth on S_UAP
p3a -> pint.Quantity
    Anoxic growth on S_S
p3b -> pint.Quantity
    Anoxic growth on S_BAP
p3c -> pint.Quantity
    Anoxic growth on S_UAP
p4 -> pint.Quantity
    Decay of heterotrophs
p5 -> pint.Quantity
    Hydrolysis of organic compounds
p6 -> pint.Quantity
    Hydrolysis of organic Nitrogen
p7 -> pint.Quantity
    Hydrolysis of X_EPS
p8 -> pint.Quantity
    Aerobic growth of autotrophs
p9 -> pint.Quantity
    Decay of autotrophs
"""
from math import exp
import pint
from typing import Union
import numpy as np

from parameters import *

MATRIX_COMPONENTS = ["S_I", "S_S", "X_I", "X_S", "X_H", "X_EPS", "S_UAP",
                     "S_BAP", "X_A", "X_P", "S_O", "S_NO", "S_N2", "S_NH",
                     "S_ND", "X_ND", "S_ALK"]
MATRIX_PROCESSES = ["p1", "p2a", "p2b", "p2c", "p3a", "p3b", "p3c", "p4",
                    "p5", "p6", "p7", "p8", "p9"]
MATRIX_COMPOSITION = ["ThOD", "Nitrogen", "Ionic Charge"]


def build_petersen_matrix() -> np.ndarray:
    """Petersen Matrix

    Output
    ------
    m: np.ndarray
        Petersen matrix
    """

    m = np.zeros((13, 17))

    # Ammonification
    m[0, 13] = 1
    m[0, 14] = -1
    m[0, 16] = 1/14

    # Aerobic Growth on S_S
    m[1, 1] = -1 / Y_H
    m[1, 4] = 1 - f_EPSh
    m[1, 5] = f_EPSh
    m[1, 6] = gamma_H / Y_H
    m[1, 10] = x2a
    m[1, 13] = y2a
    m[1, 16] = -i_XB / 14

    # Aerobic Growth on S_BAP
    m[2, 4] = 1 - f_EPSh
    m[2, 5] = f_EPSh
    m[2, 7] = -1 / Y_SMP
    m[2, 10] = x2b
    m[2, 13] = y2b
    m[2, 16] = -i_XB / 14

    # Aerobic Growth on S_UAP
    m[3, 4] = 1 - f_EPSh
    m[3, 5] = f_EPSh
    m[3, 6] = -1 / Y_SMP
    m[3, 10] = x2c
    m[3, 13] = y2c
    m[3, 16] = -i_XB / 14

    # Anoxic Growth on S_S
    m[4, 1] = -1 / Y_H
    m[4, 4] = 1 - f_EPSh
    m[4, 5] = f_EPSh
    m[4, 6] = gamma_H / Y_H
    m[4, 11] = x3a
    m[4, 12] = -x3a
    m[4, 13] = y3a
    m[4, 16] = (1 - Y_H) / (40 * Y_H) - i_XB / 14

    # Anoxic Growth on S_BAP
    m[5, 4] = 1 - f_EPSh
    m[5, 5] = f_EPSh
    m[5, 7] = -1 / Y_SMP
    m[5, 11] = x3b
    m[5, 12] = -x3b
    m[5, 13] = y3b
    m[5, 16] = (1 - Y_H) / (40 * Y_H) - i_XB / 14

    # Anoxic Growth on S_UAP
    m[6, 4] = 1 - f_EPSh
    m[6, 5] = f_EPSh
    m[6, 6] = -1 / Y_SMP
    m[6, 11] = x3c
    m[6, 12] = -x3c
    m[6, 13] = y3c
    m[6, 16] = (1 - Y_H) / (40 * Y_H) - i_XB / 14

    # Decay of heterotrophs
    m[7, 3] = 1 - f_P - f_EPSdh - f_BAP
    m[7, 4] = -1
    m[7, 5] = f_EPSdh
    m[7, 7] = f_BAP
    m[7, 9] = f_P
    m[7, 15] = i_XP - f_P * i_XP

    # Hydrolysis of organic compounds
    m[8, 1] = 1
    m[8, 3] = -1

    # Hydrolysis of organic Nitrogen
    m[9, 14] = 1
    m[9, 15] = -1

    # Hydrolysis of X_EPS
    m[10, 1] = f_S
    m[10, 5] = -1
    m[10, 6] = 1 - f_S
    m[10, 14] = i_XEPS - i_XBAP * (1 - f_S)

    # Aerobic growth of autotrophs
    m[11, 5] = f_EPSa
    m[11, 6] = gamma_A / Y_A
    m[11, 8] = 1 - f_EPSa
    m[11, 10] = -(64 / 14 - Y_A) / Y_A
    m[11, 11] = 1 / Y_A
    m[11, 13] = -i_XB - 1 / Y_A
    m[11, 16] = -i_XB / 14 - 1 / (7 * Y_A)

    # Decay of autotrophs
    m[12, 3] = 1 - f_P - f_EPSda - f_BAP
    m[12, 5] = f_EPSda
    m[12, 7] = f_BAP
    m[12, 8] = -1
    m[12, 9] = f_P
    m[12, 15] = i_XP - f_P * i_XP

    return(m)


def build_composition_matrix():
    """Composition Matrix

    Output
    ------
    m: np.ndarray
        Composition matrix
    """
    m = np.zeros((3, 17))
    m[0, :] = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1, -64/14, -24/14, 0, 0, 0, 0]
    m[1, :] = [0, 0, 0, 0, i_XB, i_XEPS, 0, i_XBAP,
               i_XB, i_XP, 0, 1, 1, 1, 1, 1, 0]
    m[2, :] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1/14, 0, 1/14, 0, 0, -1]
    return(m)


# Process Rate Equations
def p1(S_ND: pint.Quantity, X_H: pint.Quantity) -> pint.Quantity:
    """Ammonification"""
    p = k_a * S_ND * X_H
    return(p)


def p2a(S_O: pint.Quantity, S_S: pint.Quantity,
        X_H: pint.Quantity) -> pint.Quantity:
    """Aerobic growth on S_S"""
    n = S_S * S_O * X_H
    d = (K_S + S_S) * (K_OH + S_O)
    p = mu_H * n / d
    return(p)


def p2b(T: pint.Quantity, S_BAP: pint.Quantity, S_O: pint.Quantity,
        S_ALK: pint.Quantity, X_H: pint.Quantity) -> pint.Quantity:
    """Aerobic growth on S_BAP"""
    T = T.to(ureg.degC).magnitude
    n = S_BAP * S_O * S_ALK * X_H
    d = (K_BAP + S_BAP) * (K_OH + S_O) * (K_ALKH + S_ALK)
    p = exp(-0.069 * (20 - T)) * mu_BAP * n / d
    return(p)


def p2c(T: pint.Quantity, S_UAP: pint.Quantity, S_O: pint.Quantity,
        S_ALK: pint.Quantity, X_H: pint.Quantity) -> pint.Quantity:
    """Aerobic growth on S_UAP"""
    T = T.to(ureg.degC).magnitude
    n = S_UAP * S_O * S_ALK * X_H
    d = (K_UAP + S_UAP) * (K_OH + S_O) * (K_ALKH + S_ALK)
    p = exp(-0.069 * (20 - T)) * mu_UAP * n / d
    return(p)


def p3a(S_O: pint.Quantity, S_NO: pint.Quantity, S_S: pint.Quantity,
        X_H: pint.Quantity) -> pint.Quantity:
    """Anoxic growth on S_S"""
    n = S_S * K_OH * S_NO * X_H
    d = (K_S + S_S) * (K_OH + S_O) * (K_NO + S_NO)
    p = mu_H * eta_g * n / d
    return(p)


def p3b(T: pint.Quantity, S_BAP: pint.Quantity, S_O: pint.Quantity,
        S_NO: pint.Quantity, S_ALK: pint.Quantity,
        X_H: pint.Quantity) -> pint.Quantity:
    """Anoxic growth on S_BAP"""
    T = T.to(ureg.degC).magnitude
    n = S_BAP * K_OH * S_NO * S_ALK * X_H
    d = (K_BAP + S_BAP) * (K_OH + S_O) * (K_NO + S_NO) * (K_ALKH + S_ALK)
    p = exp(-0.069 * (20 - T)) * mu_BAP * eta_g * n / d
    return(p)


def p3c(T: pint.Quantity, S_UAP: pint.Quantity, S_O: pint.Quantity,
        S_NO: pint.Quantity, S_ALK: pint.Quantity,
        X_H: pint.Quantity) -> pint.Quantity:
    """Anoxic growth on S_UAP"""
    T = T.to(ureg.degC).magnitude
    n = S_UAP * K_OH * S_NO * S_ALK * X_H
    d = (K_UAP + S_UAP) * (K_OH + S_O) * (K_NO + S_NO) * (K_ALKH + S_ALK)
    p = exp(-0.069 * (20 - T)) * mu_UAP * eta_g * n / d
    return(p)


def p4(X_H: pint.Quantity) -> pint.Quantity:
    """Decay of heterotrophs"""
    p = b_H * X_H
    return(p)


def p5(S_O: pint.Quantity, S_NO: pint.Quantity, X_S: pint.Quantity,
       X_H: pint.Quantity) -> pint.Quantity:
    """Hydrolysis of organic compounds"""
    n1 = X_S
    n2 = S_O
    n3 = K_OH * S_NO
    d1 = K_X + (X_S / X_H)
    d2 = K_OH + S_O
    d3 = (K_OH + S_O) * (K_NO + S_NO)
    p = k_h * (n1 / d1) * (n2 / d2 + eta_h * n3 / d3)
    return(p)


def p6(S_O: pint.Quantity, S_NO: pint.Quantity, X_S: pint.Quantity,
       X_H: pint.Quantity, X_ND: pint.Quantity) -> pint.Quantity:
    """Hydrolysis of organic Nitrogen"""
    p = p5(S_O, S_NO, X_S, X_H) * (X_ND / X_S)
    return(p)


def p7(T: pint.Quantity, X_EPS: pint.Quantity) -> pint.Quantity:
    """Hydrolysis of X_EPS"""
    T = T.to(ureg.degC).magnitude
    p = exp(-0.11 * (20 - T)) * k_hEPS * X_EPS
    return(p)


def p8(S_NH: pint.Quantity, S_O: pint.Quantity,
       X_A: pint.Quantity) -> pint.Quantity:
    """Aerobic growth of autotrophs"""
    n = S_NH * S_O * X_A
    d = (K_NH + S_NH) * (K_OA + S_O)
    p = mu_A * n / d
    return(p)


def p9(X_A: pint.Quantity) -> pint.Quantity:
    """Decay of autotrophs"""
    p = b_A * X_A
    return(p)
