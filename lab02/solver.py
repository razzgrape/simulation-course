import numpy as np
from numba import njit

@njit
def solve_heat(L, T0, T_left, T_right,
               rho, c, lam,
               tau, h, t_end):

    N = max(3, round(L / h) + 1)
    h_e = L / (N - 1)

    T = np.full(N, float(T0))
    T[0], T[-1] = float(T_left), float(T_right)

    A_coeff = lam / h_e**2
    C_coeff = lam / h_e**2
    B_coeff = 2.0 * lam / h_e**2 + (rho * c / tau)
    gamma   = rho * c / tau

    n_steps = max(1, round(t_end / tau))

    alpha = np.zeros(N)
    beta  = np.zeros(N)

    for _ in range(n_steps):

        F = -gamma * T

        alpha[1] = A_coeff / B_coeff
        beta[1]  = (C_coeff * T[0] - F[1]) / B_coeff

        for i in range(2, N - 1):
            denom    = B_coeff - C_coeff * alpha[i - 1]
            alpha[i] = A_coeff / denom
            beta[i]  = (C_coeff * beta[i - 1] - F[i]) / denom

        T_new = np.zeros(N)
        T_new[0]  = T_left
        T_new[-1] = T_right

        T_new[-2] = alpha[N - 2] * T_new[-1] + beta[N - 2]

        for i in range(N - 3, 0, -1):
            T_new[i] = alpha[i] * T_new[i + 1] + beta[i]

        T = T_new

    return T, h_e, N