"""
Y-axis of plots: 
- COD = S_S + S_I +  X_S + X_H + X_A + X_P + X_I + X_EPS + S_UAP + S_BAP mg/L
- S_SMP = S_UAP + S_BAP mg/L
- X_EPS mg/L
- Nitrogen = i_XB * X_H + i_XEPS * X_EPS + i_XBAP * S_BAP + i_XB * X_A 
             + i_XP * X_P + S_NO + S_N2 + S_NH + S_ND + X_ND mg/L
- Membrane Resistance
- Trans-Membrane Pressure kPa 

Varied Values:
- MLSS concentration = [3, 15, 30] g/L
- Oxygen concentration = [0, 1, 4.5] mg/L
- Temperature = [5, 20, 30] C
"""
import matplotlib.pyplot as plt
from parameters import ureg, i_XB, i_XEPS, i_XBAP, i_XP
from state import starting_state
import integrated_model
import pint

# Defining units for shorthand
gL = (ureg.g / ureg.liter)
mgL = (ureg.mg / ureg.liter)
C = ureg.degC
day = ureg.day


def make_plots(data: list, fsz: int = 10) -> None:
    ylabels = ['COD (mg/L)', '$S_{SMP}$ (mg/L)', '$X_{EPS}$ (mg/L)',
               'Total Nitrogen (mg/L)', 'Membrane Resistance (1/m)',
               'TMP (kPa)']
    dict_keys = ['COD (mg/L)', 'S_SMP (mg/L)', 'X_EPS (mg/L)', 'N (mg/L)',
                 'R_t (1/m)', 'TMP (kPa)']
    file_names = ['COD.png', 'SMP.png', 'EPS.png', 'Nitrogen.png',
                  'Resistance.png', 'TMP.png']

    for i in range(len(dict_keys)):
        fig, ax = plt.subplots(nrows=1, ncols=1)
        for data_set in data:
            x_mlss = data_set['X_MLSS (g/L)']
            s_o = data_set['S_O,in (mg/L)']
            temperature = data_set['T (C)']
            data_label = ("$X_{MLSS} = %.1f g/L | " %x_mlss
                          + "S_{O,in} = %.1f mg/L | " %s_o
                          + "T = %.1f C$" %temperature)
            x_data = data_set['t (day)']
            y_data = data_set[dict_keys[i]]
            ax.plot(x_data, y_data, label=data_label)
        ax.set_xlabel('Time (day)', fontsize=fsz)
        ax.set_ylabel(ylabels[i], fontsize=fsz)
        ax.legend(fontsize=fsz)
        fig.set_size_inches(8, 6)
        fig.set_tight_layout(True)
        plt.savefig(f'plots/{file_names[i]}')


def generate_data(time: pint.Quantity = 1 * day) -> None:
    time = time.to(day)
    # The parameters to vary
    X_MLSS = [3 * gL, 15 * gL, 30 * gL]
    oxygen = [0 * mgL, 1 * mgL, 4.5 * mgL]
    temperature = [ureg.Quantity(5, ureg.degC),
                   ureg.Quantity(20, ureg.degC),
                   ureg.Quantity(30, ureg.degC)]
    parameters = ((X_MLSS[0], oxygen[1], temperature[1]),
                  (X_MLSS[1], oxygen[1], temperature[1]),
                  (X_MLSS[2], oxygen[1], temperature[1]),
                  (X_MLSS[1], oxygen[0], temperature[1]),
                  (X_MLSS[1], oxygen[2], temperature[1]),
                  (X_MLSS[1], oxygen[1], temperature[0]),
                  (X_MLSS[1], oxygen[1], temperature[2]))
    n_sims = len(parameters)
    all_data = []
    for i, (x_m, s_o, T) in enumerate(parameters):
        # Generate a new starting state
        state = starting_state.copy()
        state['X_MLSS'] = x_m
        state['in_X_MLSS'] = x_m
        state['in_S_O'] = s_o
        state['temperature'] = T
        # Create a new model and reset timer
        model = integrated_model.MBRModel(state)
        t = 0 * ureg.s
        t_step = 60 * 5 * ureg.s
        data = {
            't (day)': [],
            'COD (mg/L)': [],
            'S_SMP (mg/L)': [],
            'X_EPS (mg/L)': [],
            'N (mg/L)': [],
            'R_t (1/m)': [],
            'TMP (kPa)': [],
            'X_MLSS (g/L)': x_m.to(gL).magnitude,
            'S_O,in (mg/L)': s_o.to(mgL).magnitude,
            'T (C)': T.to(C).magnitude
        }
        while t < time:
            # Simulate model
            days = t.to('day').magnitude
            print(f'\rSimulation {i+1}/{n_sims}| {days:.3f}/{time:.3f} days  ',
                  end='')
            state = model.step_model(t_step)
            data = update_data(data, state)
            t += t_step
        all_data.append(data)
    make_plots(all_data)


def update_data(data: dict, state: dict) -> dict:
    data['t (day)'].append(state['time'].to(day).magnitude)
    data['R_t (1/m)'].append(state['R_t'].to(ureg.m ** -1).magnitude)
    data['TMP (kPa)'].append(state['TMP'].to(ureg.kPa).magnitude)
    data['X_EPS (mg/L)'].append(state['X_EPS'].to(mgL).magnitude)

    COD = (state['S_S'] + state['S_I'] +  state['X_S'] + state['X_H']
           + state['X_A'] + state['X_P'] + state['X_I'] + state['X_EPS']
           + state['S_UAP'] + state['S_BAP'])
    S_SMP = state['S_BAP'] + state['S_UAP']
    N = (i_XB * state['X_H'] + i_XEPS * state['X_EPS']
         + i_XBAP * state['S_BAP'] + i_XB * state['X_A'] + i_XP * state['X_P']
         + state['S_NO'] + state['S_N2'] + state['S_NH'] + state['S_ND']
         + state['X_ND'])

    data['COD (mg/L)'].append(COD.to(mgL).magnitude)
    data['S_SMP (mg/L)'].append(S_SMP.to(mgL).magnitude)
    data['N (mg/L)'].append(N.to(mgL).magnitude)
    return(data)


if __name__ == '__main__':
    generate_data()