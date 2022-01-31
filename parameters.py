"""
Contains values for model parameters
"""
import pint

ureg = pint.UnitRegistry()
ureg.load_definitions('custom_units.txt')

# Yield coefficient for heterotrophic growth on S_UAP and S_BAP
Y_smp = 0.45 * (ureg.gCOD / ureg.gCOD)

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
