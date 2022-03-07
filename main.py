import numpy as np
import pint

from typing import Union

from parameters import ureg
import parameters as param

influent_dry = {
    'S_I': 30.00 * (ureg.gCOD / ureg.m ** 3),
    'S_S': 69.50 * (ureg.gCOD / ureg.m ** 3),
    'X_I': 51.20 * (ureg.gCOD / ureg.m ** 3),
    'X_S': 202.32 * (ureg.gCOD / ureg.m ** 3),
    'X_H': 28.17 * (ureg.gCOD / ureg.m ** 3),
    'S_NH': 31.56 * (ureg.gN / ureg.m ** 3),
    'S_ND': 6.95 * (ureg.gN / ureg.m ** 3),
    'X_ND': 10.59 * (ureg.gN / ureg.m ** 3),
    'S_ALK': 7.00 * (ureg.molHCO3 / ureg.m ** 3),
    'Q_i,av': 18446.33 * (ureg.m ** 3 / ureg.day),
    'Q_i,max': 32180 * (ureg.m ** 3 / ureg.day)
}

influent_rain = {
    'S_I': 25.96 * (ureg.gCOD / ureg.m ** 3),
    'S_S': 60.13 * (ureg.gCOD / ureg.m ** 3),
    'X_I': 44.30 * (ureg.gCOD / ureg.m ** 3),
    'X_S': 175.05 * (ureg.gCOD / ureg.m ** 3),
    'X_H': 24.37 * (ureg.gCOD / ureg.m ** 3),
    'S_NH': 27.30 * (ureg.gN / ureg.m ** 3),
    'S_ND': 6.01 * (ureg.gN / ureg.m ** 3),
    'X_ND': 9.16 * (ureg.gN / ureg.m ** 3),
    'S_ALK': 7.00 * (ureg.molHCO3 / ureg.m ** 3),
    'Q_i,av': 21319.75 * (ureg.m ** 3 / ureg.day),
    'Q_i,max': 52126 * (ureg.m ** 3 / ureg.day)
}

influent_storm = {
    'S_I': 28.03 * (ureg.gCOD / ureg.m ** 3),
    'S_S': 64.93 * (ureg.gCOD / ureg.m ** 3),
    'X_I': 51.92 * (ureg.gCOD / ureg.m ** 3),
    'X_S': 193.32 * (ureg.gCOD / ureg.m ** 3),
    'X_H': 27.25 * (ureg.gCOD / ureg.m ** 3),
    'S_NH': 29.48 * (ureg.gN / ureg.m ** 3),
    'S_ND': 6.49 * (ureg.gN / ureg.m ** 3),
    'X_ND': 10.24 * (ureg.gN / ureg.m ** 3),
    'S_ALK': 7.00 * (ureg.molHCO3 / ureg.m ** 3),
    'Q_i,av': 19744.72 * (ureg.m ** 3 / ureg.day),
    'Q_i,max': 60000 * (ureg.m ** 3 / ureg.day)
}

def main(time: Union[float, pint.Quantity], influent: dict):
    pass


if __name__ == "__main__":
    main()
