import numpy as np
import pint

import sludge


class Bioreactor:
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
        p = self.calculate_process_rates(state)
        m = sludge.build_petersen_matrix()
        rates = self.component_rates(m, p)
        new_state = self.material_balance(t_step, state, rates)
        return(new_state)

    def component_rates(self, m: np.ndarray, p: list) -> np.ndarray:
        rates = [0] * m.shape[1]
        for i in range(m.shape[0]):
            for j in range(m.shape[1]):
                rates[j] += m[i, j] * p[i]
        return(rates)

    def calculate_process_rates(self, state: dict) -> list:
        """Find the rates of all processes in the reactor

        Parameters
        ----------
        state: dict
            Contains all state variables
        Output
        ------
        p: list
            A list containing all process rates in the following order:
            [p1, p2a, p2b, p2c, p3a, p3b, p3c, p4, p5, p6, p7, p8, p9]
        """
        # Unpack state variables
        T = state['temperature']
        S_S = state["S_S"]
        X_S = state["X_S"]
        X_H = state["X_H"]
        X_EPS = state["X_EPS"]
        S_UAP = state["S_UAP"]
        S_BAP = state["S_BAP"]
        X_A = state["X_A"]
        S_O = state["S_O"]
        S_NO = state["S_NO"]
        S_NH = state["S_NH"]
        S_ND = state["S_ND"]
        X_ND = state["X_ND"]
        S_ALK = state["S_ALK"]

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
        p = [p1, p2a, p2b, p2c, p3a, p3b, p3c, p4, p5, p6, p7, p8, p9]

        return(p)

    def material_balance(self, t_step: pint.Quantity, state: dict,
                         rates: np.ndarray) -> dict:
        """Material balance

        Parameters
        ----------
        t_step: pint.Quantity
            The time step
        state: dict
            Contains all state variables
        rates: list
            Contains all rates of change of concentrations due to reactions
        Output
        ------
        new_state: dict
            Contains all updated state variables
        """
        new_state = state.copy()
        V = state['volume']
        Q_in = state['Q_in']
        Q_out = state['Q_out']

        for i, key in enumerate(["S_I", "S_S", "X_I", "X_S", "X_H", "X_EPS",
                                 "S_UAP", "S_BAP", "X_A", "X_P", "S_O", "S_NO",
                                 "S_N2", "S_NH", "S_ND", "X_ND", "S_ALK"]):
            # Update concentrations
            C_0 = state[f'in_{key}']
            C = state[key]
            m_total = C * V + (Q_in * C_0 - Q_out * C + rates[i] * V) * t_step
            if m_total < 0:
                m_total *= 0
            new_state[key] = m_total / V
        # Update MLSS (no reactions which is why its not in loop)
        C_MLSS_0 = state['in_X_MLSS']
        C_MLSS = state['X_MLSS']
        m_MLSS = C_MLSS * V + (Q_in * C_MLSS_0 - Q_out * C_MLSS) * t_step
        new_state['X_MLSS'] = m_MLSS / V

        Delta_V = Q_in - Q_out
        new_state['volume'] = V + Delta_V * t_step

        # Calculate TSS concentration
        new_state['X_TSS'] = 0.75 * (new_state['X_S'] + new_state['X_H']
                                     + new_state['X_A'] + new_state['X_P']
                                     + new_state['X_I'] + new_state['X_EPS'])
        return(new_state)
