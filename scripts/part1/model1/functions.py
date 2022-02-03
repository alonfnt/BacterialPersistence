import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from tqdm import tqdm

from model_equations import growth, decay
from model_equations import n0, S0, f, gamma, T_cycle


##########################
## Analytical functions ##
##########################
# computing final population of cycle one
# not used
def analytical_cycle(d0, lag, ab_args, r, extinction=False):
    p, Tab = ab_args
    g = np.zeros_like(d0, dtype=float)
    t, T = 0, 0

    # with antibiotics
    if r < p:
        T = Tab

    # after antibiotics / without antibiotics
    d_T = d0 * np.exp(-T / lag)
    g_T = (d_T - d0 * np.exp(-gamma * T)) / (gamma * lag - 1)

    # removing extinct species
    if extinction:
        g[g_T + d_T < n_min] = 0
    #alive = (g_T + d_T > 0)

    # time that substrate is eaten up
    exp_t = (S0 + sum(g_T)) * np.exp(T) / sum(g_T + d_T / (1 + lag))
    t = np.log(exp_t)

    # population and time when substrate is eaten up
    d = d0 * np.exp(-t / lag)
    g = (g_T + d_T / (1 + lag)) * np.exp(t-T) - d / (1 + lag)

    n_T = g_T + d_T
    n = g + d

    return n, t
    
    
# computing analytical fitness
def analytical_fitness(lag, ab_args):
    p, T0, Tab = ab_args
    T = T0 + Tab
	
    # fitness without antibiotics
    fit_0 = np.log(1 / (1 + lag))
    
    # fitness with antibiotics
    fit_p = np.log(((np.exp(-Tab / lag) - np.exp(-gamma * Tab)) / (gamma * lag - 1) + (np.exp(-Tab / lag) / (1 + lag))) * np.exp(-T0 / lag)
                   + (np.exp(T0) - np.exp(-T0 / lag)) * np.exp(-gamma * Tab) / (1 + lag))
                   
    return (1 - p) * fit_0 + p * (-T + fit_p)


# numerically computing optimal lag time
def optimal_lag(ab_args, resolution=0.1):
    p, T0, Tab = ab_args
    T = T0 + Tab
    lag = np.arange(0, T, resolution) + 10**(-4)

    F = analytical_fitness(lag, ab_args)		# computing analytical fitness
    lag_opt = lag[F == F.max()][0]		        # identifying lag time corresponding to max fitness

    return lag_opt, F.max()



##########################
## Simulation functions ##
##########################
# function for detecting when system runs out of nutrients, used in solve_ivp
def famine(t, g, lag): return g[-1]
famine.terminal = True
famine.direction = -1


# performing one feast-famine cycle
def solve_cycle(n_t, r, lag, ab_args, res):
    p, T0, Tab = ab_args
    dt_max = 1

    T = T0 + Tab
    t_low = 0
    lags = np.array([optimal_lag(ab_args, res)[0], lag])		# lags[0] = optimal lag time, lags[1] = competitor lag time
    S_frac = np.zeros(2, dtype=float)				# empty array for consumption fraction

    # cycles with antibiotics
    if r < p:
        # solving system before antibiotics
        sol_cycle = solve_ivp(growth, [0, T0], n_t, args=(lags, ), max_step=dt_max)
        n_t = np.array([sol_cycle.y[i][-1] for i in range(5)])
        # updating consumption fraction
        S_frac += n_t[2:4] / S0

        # solving system during antibiotics
        sol_cycle = solve_ivp(decay, [T0, T], n_t, args=(lags, ), max_step=dt_max)
        n_t = np.array([sol_cycle.y[i][-1] for i in range(5)])
        # updating consumption fraction
        S_frac -= n_t[2:4] / S0

        # updating lower t
        t_low = T

    # cycles without antibiotics / solving system after antibiotics
    sol_cycle = solve_ivp(growth, [t_low, T_cycle], n_t, args=(lags, ), events=famine, max_step=dt_max)
    n_t = np.array([sol_cycle.y[i][-1] for i in range(5)])
    d_t = f * np.array([n_t[i] + n_t[2 + i] for i in range(2)])	# starvation induced dormancy
    
    S_frac += n_t[2:4] / S0					# updating consumption fraction
    n_t = np.concatenate([d_t, np.zeros(2), [10**9]])		# preparing array for next cycle

    return n_t, S_frac


# looping through lag times and performing several feast-famine cycles for each
def compete(lag_arr, ab_args, sim_args):
    p, T0, Tab = ab_args
    tot_cycles, res = sim_args

    # defining output arrays
    S_frac_optimal = np.zeros([len(lag_arr), tot_cycles])
    S_frac_competitor = np.zeros([len(lag_arr), tot_cycles])

    # looping through lag times
    il = 0
    for lag in tqdm(lag_arr):
        n_t = np.array([n0, n0, 0, 0, 10 ** 9])			# initial populations

	# performing feast-famine cycles
        for cycle in range(tot_cycles):
            r = np.random.rand()
            n_t, sfrac = solve_cycle(n_t, r, lag, ab_args, res)
            S_frac_optimal[il, cycle] = sfrac[0]
            S_frac_competitor[il, cycle] = sfrac[1]

        il += 1

    return S_frac_optimal.mean(axis=1), S_frac_competitor.mean(axis=1)
    
    
####################
## Plot behaviour ##
####################
# setting plotting parameters
font = {'family': 'Times New Roman',
        'weight': 'normal',
        'size': 20}

mpl.rc('font', **font)
fs_text = 14

# plotting behavior of system, both with ode and analytical
# Not working
def plot_cycles(lag, p, r_arr, T_ab, extinction, savefig=False):
    tot_cycles = len(r_arr)
    # integration
    d0 = np.array([n0, n0])
    g0 = [0, 0, S0]
    t_span = [0, T_cycle]

    g1, g2, S = [], [], []
    d1, d2 = [], []
    time = []

    n_ana = np.array([n0/f, n0/f])
    n1_a, n2_a, time_ana = [], [], []

    t_S = np.zeros([2, tot_cycles])


    for ic in range(tot_cycles):
        # diluting
        n_ana *= f
        time_ana.append(ic * T_cycle + 1)
        n1_a.append(n_ana[0])
        n2_a.append(n_ana[1])

        args = (d0, lag, p, r_arr[ic], T_ab, extinction)
        sol_cycle = solve_ivp(constant_rate_ode, t_span, g0, args=args, max_step=0.1)
        n_ana, n_T_ana, t_ana = analytical_cycle(n_ana, lag, p, r_arr[ic], T_ab, extinction)

        # saving data
        time.append(sol_cycle.t + ic * T_cycle)
        d1.append(d0[0] * np.exp(-sol_cycle.t / (lag[0] + 0.1)))
        d2.append(d0[1] * np.exp(-sol_cycle.t / (lag[1] + 0.1)))
        g1.append(sol_cycle.y[0])
        g2.append(sol_cycle.y[1])
        S.append(sol_cycle.y[2])

        t_S[:, ic] = min(sol_cycle.t[sol_cycle.y[2] <= 0]), t_ana

        time_ana.append(T_ab * (r_arr[ic] < p) + ic * T_cycle)
        time_ana.append(t_ana + ic * T_cycle)
        time_ana.append((1+ic) * T_cycle)
        n1_a.append(n_T_ana[0])
        n2_a.append(n_T_ana[1])
        for i in range(2):
            n1_a.append(n_ana[0])
            n2_a.append(n_ana[1])

        # preparing new cycle
        d0 = np.array([f * (sol_cycle.y[0][-1] + d1[-1][-1]), f * (sol_cycle.y[1][-1] + d2[-1][-1])])

    time = np.concatenate(time)
    d1 = np.concatenate(d1)
    d2 = np.concatenate(d2)
    g1 = np.concatenate(g1)
    g2 = np.concatenate(g2)
    S = np.concatenate(S)

    fig, ax = plt.subplots(2, 1, figsize=(8, 3.5), sharex="all", sharey="all")
    fig.suptitle("p = " + str(np.round(p, 2)) + ", lag1 = " + str(np.round(lag[0], 1)) + ", lag2 = " + str(
        np.round(lag[1], 1)) + ", T = " + str(T_ab))

    ax[0].plot(time, S, alpha=0.4, label="substrate")
    ax[1].plot(time, S, alpha=0.4, label="substrate")

    #ax[0].set(title="Integration ODEs")
    ax[0].set(title="Species 1")
    ax[0].plot(time, g1+d1, color="tab:orange", label="ODE")
    ax[1].plot(time, g2+d2, color="tab:green", label="ODE")

    #ax[1].set(title="Analytical approximation")
    ax[1].set(title="Species 2")
    ax[0].plot(time_ana, n1_a, "--", color="k", alpha=0.6, label="analytical")
    ax[1].plot(time_ana, n2_a, "--", color="k", alpha=0.6, label="analytical")


    ax[0].set(ylabel="population")
    ax[1].set(xlabel="time [h]", ylabel="population")
    ax[0].set(yscale="log", ylim=(1, 10 ** 10))
    ax[0].legend(bbox_to_anchor=(1, 1), loc="best")
    ax[1].legend(bbox_to_anchor=(1, 1), loc="best")

    fig.tight_layout()
    fig.show()
    if savefig:
        fig.savefig("figs/single_sample.png", dpi=100)

    return np.array([[d1, d2], [g1, g2]]), t_S, S
