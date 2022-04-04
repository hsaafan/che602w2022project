import parameters
from parameters import ureg

starting_state = {
    'time': 0 * (ureg.s),
    'temperature': ureg.Quantity(20, ureg.degC),
    'volume': 1500 * (ureg.m ** 3),
    # Starting membrane values
    'R_i': 0 * (1 / ureg.m),
    'R_r': 0 * (1 / ureg.m),
    'R_t': parameters.R_m,
    'v_sg': 1 * (ureg.cm / ureg.s),  # (Janus 2013, p.243)
    'alpha_c': 0 * (ureg.m / ureg.kg),
    'm_rback': 0 * (ureg.kg / ureg.m ** 2 / ureg.s),
    'tau_w': 0 * (ureg.Pa),
    'TMP': 0 * (ureg.Pa),
    # Influent concentrations (From (Janus 2013, p.264))
    'in_S_I': 9.00 * (ureg.gCOD / ureg.m ** 3),
    'in_S_S': 69.50 * (ureg.gCOD / ureg.m ** 3),
    'in_X_I': 51.20 * (ureg.gCOD / ureg.m ** 3),
    'in_X_S': 202.32 * (ureg.gCOD / ureg.m ** 3),
    'in_X_H': 26.76 * (ureg.gCOD / ureg.m ** 3),
    'in_X_EPS': 1.41 * (ureg.gCOD / ureg.m ** 3),
    'in_S_UAP': 0 * (ureg.gCOD / ureg.m ** 3),
    'in_S_BAP': 21.00 * (ureg.gCOD / ureg.m ** 3),
    'in_X_A': 0 * (ureg.gCOD / ureg.m ** 3),
    'in_X_P': 0 * (ureg.gCOD / ureg.m ** 3),
    'in_S_O': 0 * (ureg.gCOD / ureg.m ** 3),
    'in_S_NO': 0 * (ureg.gN / ureg.m ** 3),
    'in_S_N2': 0 * (ureg.gN / ureg.m ** 3),
    'in_S_NH': 31.56 * (ureg.gN / ureg.m ** 3),
    'in_S_ND': 6.95 * (ureg.gN / ureg.m ** 3),
    'in_X_ND': 9.37 * (ureg.gN / ureg.m ** 3),
    'in_S_ALK': 7.00 * (ureg.molHCO3 / ureg.m ** 3),
    'in_X_MLSS': 13400 * (ureg.gCOD / ureg.m ** 3),  # (Lindamulla, 2021)
    # Reactor concentrations (Assume equal to influent to begin with)
    'S_I': 9.00 * (ureg.gCOD / ureg.m ** 3),
    'S_S': 69.50 * (ureg.gCOD / ureg.m ** 3),
    'X_I': 51.20 * (ureg.gCOD / ureg.m ** 3),
    'X_S': 202.32 * (ureg.gCOD / ureg.m ** 3),
    'X_H': 26.76 * (ureg.gCOD / ureg.m ** 3),
    'X_EPS': 1.41 * (ureg.gCOD / ureg.m ** 3),
    'S_UAP': 0 * (ureg.gCOD / ureg.m ** 3),
    'S_BAP': 21.00 * (ureg.gCOD / ureg.m ** 3),
    'X_A': 0 * (ureg.gCOD / ureg.m ** 3),
    'X_P': 0 * (ureg.gCOD / ureg.m ** 3),
    'S_O': 0 * (ureg.gCOD / ureg.m ** 3),
    'S_NO': 0 * (ureg.gN / ureg.m ** 3),
    'S_N2': 0 * (ureg.gN / ureg.m ** 3),
    'S_NH': 31.56 * (ureg.gN / ureg.m ** 3),
    'S_ND': 6.95 * (ureg.gN / ureg.m ** 3),
    'X_ND': 9.37 * (ureg.gN / ureg.m ** 3),
    'S_ALK': 7.00 * (ureg.molHCO3 / ureg.m ** 3),
    'X_MLSS': 13400 * (ureg.gCOD / ureg.m ** 3),  # (Lindamulla, 2021)
    # Flow rates
    'Q_in': 18446.33 * (ureg.m ** 3 / ureg.day),
    'Q_min': 10000.00 * (ureg.m ** 3 / ureg.day),
    'Q_max': 32180.00 * (ureg.m ** 3 / ureg.day),
    'Q_out': 18446.33 * (ureg.m ** 3 / ureg.day)
}

X_TSS = 0.75 * (starting_state['X_S'] + starting_state['X_H']
                + starting_state['X_A'] + starting_state['X_P']
                + starting_state['X_I'] + starting_state['X_EPS'])
starting_state['X_TSS'] = X_TSS
