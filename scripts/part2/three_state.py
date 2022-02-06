# I don't remember exactly what this is, but think I used this in the end.

import numpy as np
import matplotlib.pyplot as plt
import time
from scipy.integrate import solve_ivp

from model_equations import model3_grow, model3_kill, famine
from model_equations import n0, S0, mutation_seed

dt_max = 0.01
T_max = 50
f = 10**(-6)
eps = 0


###########################
## Antibiotic parameters ##
###########################
T0 = 4
Tab = 12.0
T = T0 + Tab
p = 1
ab_args = [p, T0, Tab]


###########################
## Defining species grid ##
###########################
lag = 1
omega = 2
delta = 0.1
bac_args = [lag, delta, eps, omega]
gamma = 1


################
## Simulation ##
n_t = np.array([1000, 0, 0, 10**9])

def solve_3state(n_t, r, bac_args, ab_args):
    p, T0, Tab = ab_args
    dt, nt, pt = [], [], []
    t = []
    N = 1

    ode_grow = model3_grow
    ode_kill = model3_kill

    T = T0 + Tab
    t_low = 0

    if r < p:
        # Before AB
        sol_cycle = solve_ivp(ode_grow, [0, T0], n_t, args=bac_args, max_step=dt_max)
        n_t = np.array([sol_cycle.y[i][-1] for i in range(3 * N + 1)])
        dt.append(sol_cycle.y[0])
        nt.append(sol_cycle.y[1])
        pt.append(sol_cycle.y[2])
        t.append(sol_cycle.t)

        # During AB
        sol_cycle = solve_ivp(ode_kill, [T0, T], n_t, args=bac_args, max_step=dt_max)
        n_t = np.array([sol_cycle.y[i][-1] for i in range(3 * N + 1)])
        dt.append(sol_cycle.y[0])
        nt.append(sol_cycle.y[1])
        pt.append(sol_cycle.y[2])
        t.append(sol_cycle.t)

        # updating lower t
        t_low = T

    # before S = 0
    sol_cycle = solve_ivp(ode_grow, [t_low, T_max], n_t, args=bac_args, events=famine, max_step=dt_max)
    dt.append(sol_cycle.y[0])
    nt.append(sol_cycle.y[1])
    pt.append(sol_cycle.y[2])
    t.append(sol_cycle.t)
    d_t = f * np.array([sol_cycle.y[i][-1] + sol_cycle.y[N + i][-1] + sol_cycle.y[2*N + i][-1] for i in range(N)])

    dt = np.concatenate(dt)
    nt = np.concatenate(nt)
    pt = np.concatenate(pt)
    t = np.concatenate(t)

    return t, dt, nt, pt


def a_b(lag, delta):
    D = np.sqrt((lag * (1 - delta) - 1) ** 2 + 4 * lag)

    a = (D - (lag * (1 - delta) - 1)) / (2 * lag)
    b = (D + (lag * (1 - delta) - 1)) / (2 * lag)

    return a, b

def ap_bp(lag, delta):
    D = np.sqrt((lag * (gamma + delta) + 1) ** 2 - 4 * gamma * lag)

    ap = (D + (lag * (gamma + delta) + 1)) / (2 * lag)
    bp = - (D - (lag * (gamma + delta) + 1)) / (2 * lag)

    return ap, bp

def before_3state(T0, d0, bac_args):
    lag, delta, eps, omega = bac_args

    sqrt = np.sqrt((omega * (1 - delta) - 1)**2 + 4 * omega)
    a, b = a_b(omega, delta)
    denom = (-2 * (omega - lag ** 2 + lag * (omega * (1 - delta) - 1)) * sqrt)

    Ag = d0 * ((lag - omega) * sqrt - lag * (omega * (1 + delta) + 1) + omega * (omega * (1 - delta) + 1)) / denom
    Bg = d0 * ((lag - omega) * sqrt + lag * (omega * (1 + delta) + 1) - omega * (omega * (1 - delta) + 1)) / denom
    Dg = d0 * (lag - omega) * (-2 * sqrt) / denom

    Ap = d0 * omega * delta * (lag * sqrt + 2 * omega + lag * (omega * (1 - delta) - 1)) / denom
    Bp = d0 * omega * delta * (lag * sqrt - 2 * omega - lag * (omega * (1 - delta) - 1)) / denom
    Dp = d0 * omega * delta * lag * (-2 * sqrt) / denom

    d_t = d0 * np.exp(-T0 / lag)
    g_t = Ag * np.exp(-a * T0) + Bg * np.exp(b * T0) + Dg * np.exp(-T0 / lag)
    p_t = Ap * np.exp(-a * T0) + Bp * np.exp(b * T0) + Dp * np.exp(-T0 / lag)

    return d_t, g_t, p_t

def test0(t, d0, bac_args):
    lag, delta, eps, omega = bac_args

    a, b = a_b(omega, delta)
    C = (lag - omega) / (omega - lag ** 2 + lag * (omega * (1 - delta) - 1))
    B0 = (d0 / lag) * ((C*(1 - a*lag) + 1) / (a + b))
    A0 = -(d0 / lag) * ((C*(1 + b*lag) + 1) / (a + b))

    p = omega * (b + delta - 1) * B0 * np.exp(b*t) - omega * (a - delta + 1) * A0 * np.exp(-a*t) - omega * d0 * (C + 1 + C*lag*(1-delta)) * np.exp(-t / lag) / lag


    return B0 * np.exp(b*t) + A0 * np.exp(-a*t) + d0 * C * np.exp(-t / lag), p

def test1(T, T0, d0, bac_args):
    lag, delta, eps, omega = bac_args
    a, b = a_b(omega, delta)
    ap, bp = ap_bp(omega, delta)

    C0 = (lag - omega) / (omega - lag ** 2 + lag * (omega * (1 - delta) - 1))
    B0 = (d0 / lag) * ((C0 * (1 - a * lag) + 1) / (a + b))
    A0 = -(d0 / lag) * ((C0 * (1 + b * lag) + 1) / (a + b))

    g_T0 = B0 * np.exp(b * T0) + A0 * np.exp(-a * T0) + d0 * C0 * np.exp(-T0 / lag)
    p_T0 = omega * (b + delta - 1) * B0 * np.exp(b * T0) - omega * (a - delta + 1) * A0 * np.exp(
        -a * T0) - omega * d0 * (C0 + 1 + C0 * lag * (1 - delta)) * np.exp(-T0 / lag) / lag

    Cp = (lag - omega) / (omega + gamma * lag ** 2 - lag * (omega * (gamma + delta) + 1))
    Ap = (b * B0 * np.exp(b * T0) - a * A0 * np.exp(-a * T0) - d0 * (C0 - Cp) * np.exp(
        -T0 / lag) / lag + bp * g_T0 - bp * d0 * Cp * np.exp(-T0 / lag)) * np.exp(ap * T0) / (bp - ap)
    Bp = (g_T0 - d0 * Cp * np.exp(-T0 / lag) - Ap * np.exp(-ap * T0)) * np.exp(bp * T0)

    g_T = Bp * np.exp(-bp * (T-T0)) + Ap * np.exp(-ap * (T-T0)) + d0 * Cp * np.exp(-(T-T0) / lag)
    # p_T = omega * (b + delta - 1) * B0 * np.exp(b*T0) - omega * (a - delta + 1) * A0 * np.exp(-a*T0) - omega * d0 * (C + 1 + C*lag*(1-delta)) * np.exp(-T0 / lag) / lag


    return g_T

t, d_t, g_t, p_t = solve_3state(n_t, 0, bac_args, ab_args)
d, g, p = before_3state(t, 1000, bac_args)
tt = np.linspace(T0, T, 100)
g_test = test1(tt, T0, 1000, bac_args)

fig, ax = plt.subplots(3, 1, sharey=True, sharex=True)
ax[0].plot(t, g_t, label="ode")
ax[1].plot(t, p_t, label="ode")

# ax[0].plot(t, g, '--', label="analytical")
# ax[1].plot(t, p, '--', label="analytical")

ax[0].plot(tt, g_test, '--')
# ax[1].plot(t, p_test)

ax[0].set(yscale="log", ylim=(1, 10**9))
ax[0].legend()
ax[1].legend()
fig.show()
