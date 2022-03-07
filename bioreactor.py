import numpy as np
import pint

from typing import Union

import sludge


def bioreactor_model(T: pint.Quantity,
                     S_I: pint.Quantity, S_S: pint.Quantity,
                     X_I: pint.Quantity, X_S: pint.Quantity,
                     X_H: pint.Quantity, X_EPS: pint.Quantity,
                     S_UAP: pint.Quantity, S_BAP: pint.Quantity,
                     X_A: pint.Quantity, X_P: pint.Quantity,
                     S_O: pint.Quantity, S_NO: pint.Quantity,
                     S_N2: pint.Quantity, S_NH: pint.Quantity,
                     S_ND: pint.Quantity, X_ND: pint.Quantity,
                     S_ALK: pint.Quantity) -> dict:
    # Calculate reaction rates
    p1 = sludge.p1(S_ND, X_H)
    p2a = sludge.p2a(S_O, S_S, X_H)
    p2b = sludge.p2b(T, S_BAP, S_O, S_ALK, X_H)
    p2c = sludge.p2c(T, S_UAP, S_O, S_ALK, X_H)
    p3a = sludge.p3a(S_O, S_NO, S_S, X_H)
    p3b = sludge.p3b(T, S_BAP, S_O, S_NO, S_ALK, X_H)
    p3c = sludge.p3c(T, S_UAP, S_O, S_NO, S_ALK, X_H)
    p4 = sludge.p4(X_H)
    p5 = sludge.p5(S_O, S_NO, X_S, X_H)
    p6 = sludge.p6(S_O, S_NO, X_S, X_H, X_ND)
    p7 = sludge.p7(T, X_EPS)
    p8 = sludge.p8(S_NH, S_O, X_A)
    p9 = sludge.p9(X_A)

    # Store data in matrices for calculation
    p_vector = np.asarray([p1, p2a, p2b, p2c, p3a, p3b,
                           p3c, p4, p5, p6, p7, p8, p9]).reshape((-1, 1))
    petersen_matrix = sludge.build_petersen_matrix()

    # Calculate rates of change
    rates = petersen_matrix.T @ p_vector
    (dS_I, dS_S, dX_I, dX_S, dX_H, dX_EPS, dS_UAP, dS_BAP, dX_A, dX_P,
     dS_O, dS_NO, dS_N2, dS_NH, dS_ND, dX_ND, dS_ALK) = rates

    return(dS_I, dS_S, dX_I, dX_S, dX_H, dX_EPS, dS_UAP, dS_BAP, dX_A, dX_P,
           dS_O, dS_NO, dS_N2, dS_NH, dS_ND, dX_ND, dS_ALK)
