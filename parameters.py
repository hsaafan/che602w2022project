"""
Contains values for model parameters
All units are included in the parenthesis following the parameter value
Units of (1) are dimensionless parameters
"""
import pint

ureg = pint.UnitRegistry()
ureg.load_definitions('custom_units.txt')

# Yield coefficient for heterotrophic growth on S_UAP and S_BAP
Y_SMP = 0.45 * (ureg.gCOD / ureg.gCOD)
# Heterotrophic biomass yield
Y_H = 0.67 / (1 + 0.0924) * (ureg.gCOD / ureg.gCOD)

# Fraction of S_UAP produced during:
gamma_H = 0.0924 * (ureg.gCOD / ureg.gCOD)  # Autotrophic growth
gamma_A = 0 * (ureg.gCOD / ureg.gCOD)       # Heterotrophic growth

# Nitrogen content of [subscript]
i_XBAP = 0.07 * (ureg.gN / ureg.gCOD)
i_XEPS = 0.07 * (ureg.gN / ureg.gCOD)

# Half saturation constant for concentration of [subscript]
K_UAP = 100 * (ureg.gCOD / ureg.m ** 3)
K_BAP = 85 * (ureg.gCOD / ureg.m ** 3)

# Max specific heterotrophic growth rate on concentration of [subscript] at 20C
mu_UAP = 0.45 * (1 / ureg.day)
mu_BAP = 0.15 * (1 / ureg.day)

# Fractions
# # S_S produced during X_EPS hydrolysis
f_S = 0.4 * (ureg.gCOD / ureg.gCOD)
# # Concentration of EPS produced during heterotrophic biomass growth
f_EPSh = 0.10 * (ureg.gCOD / ureg.gCOD)
# # Concentration of EPS produced during heterotrophic biomass decay
f_EPSdh = 0.025 * (ureg.gCOD / ureg.gCOD)
# # Concentration of EPS produced during autotrophic biomass growth
f_EPSa = 0.0 * (ureg.gCOD / ureg.gCOD)
# # Concentration of EPS produced during autotrophic biomass decay
f_EPSda = 0.0 * (ureg.gCOD / ureg.gCOD)
# # Concentration of BAP produced during biomass decay
f_BAP = 0.0215 * (ureg.gCOD / ureg.gCOD)

# Max concentration of EPS hydrolysis rate at 20C
k_hEPS = 0.17 * (1 / ureg.day)

# Clean membrane resistance
R_m = 3.0e12 * (1 / ureg.m)
# Concentration of UAP and concentration of BAP retention on the membrane
f_M = 0.5 * (1)
# Fux dependency coeﬃcient in the irreversible fouling equation
b = 6.8e-2 * (ureg.h * ureg.m ** 2 / ureg.L)
# Threshold pressure below which no cake compression occurs
Delta_P_crit = 30000 * (ureg.Pa)
# Cake compressibility factor
n = 0.25
# Empirical proportionality coeﬃcient in the cake detachment equation
gamma_m = 1500 * (1 / ureg.Pa / ureg.s)
# Static friction coeﬃcient in the cake detachment equation
lambda_m = 2e-6 * (1)

# Parameters from (Henze et al., 2000)
K_S = 20.0 * (ureg.g / ureg.m ** 3)
K_OH = 0.2 * (ureg.g / ureg.m ** 3)
K_NO = 0.5 * (ureg.g / ureg.m ** 3)
K_X = 0.03 * (ureg.g / ureg.g)
K_NH = 1.0 * (ureg.g / ureg.m ** 3)
K_OA = 0.4 * (ureg.g / ureg.m ** 3)
K_ALKH = 6.1 * (ureg.g / ureg.m ** 3)
eta_g = 0.8 * (1)
eta_h = 0.4 * (1)
Y_A = 0.24 * (ureg.g / ureg.g)
Y_H = 0.67 * (ureg.g / ureg.g)
i_XB = 0.086 * (ureg.g / ureg.g)
i_XE = 0.06 * (ureg.g / ureg.g)
i_XP = 0.06 * (ureg.g / ureg.g)
f_P = 0.08 * (1)
b_H = 0.62 * (1 / ureg.day)
b_A = 0.1 * (1 / ureg.day)
# Ammonification rate constant at 20C from (Henze et al., 2000)
k_a = 0.08 * (ureg.m ** 3 / ureg.g / ureg.day)
k_h = 3.0 * (ureg.g / ureg.g / ureg.day)
mu_A = 0.8 * (1 / ureg.day)
mu_H = 6.0 * (1 / ureg.day)
