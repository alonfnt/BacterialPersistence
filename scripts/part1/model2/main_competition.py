import numpy as np
import matplotlib.pyplot as plt
import time

from model_equations import a_b, ap_bp, lag_min, delta_max
from simulation_functions import run_competition_in_parallel

np.random.seed(18)


###########################
## Simulation parameters ##
###########################
bac_res = 20 #101                           # resolution in bacterial parameters
ab_res = 20 #100                            # resolution in antibiotic parameters
t_res = 16                              # resolution in time array
tot_cycles = 20
repetitions = 1                         # number of repetitions for ensemble average

threshold = False                       # extinction threshold on population
save_fig = False
save_data = False
data = 'new'                                    # 'new' - generate and plot new data. 'old' plot old data


###########################
## Antibiotic parameters ##
###########################
# for chosing which time parameter to keep constant.
ic = 0                                          # ic = 0 corresponds to constant T0, ic = 1 corresponds to constant Tab
T_const = 5                                     # value of the constant parameter

# defining parameter arrays
T_max = [12, 24]                                # upper bounds on meningful values for T0 and Tab
p_arr = np.linspace(0, 1, ab_res)               # probability of antibiotics
T_arr = np.linspace(0, T_max[1-ic], ab_res)     # time array

T_labels = ['T0', 'Tab']
T_values = [T_const, T_arr]
T0 = T_values[ic]
Tab = T_values[1-ic]



##########################
## Bacterial parameters ##
##########################
# empty arrays parameters
lag = np.ones([2, bac_res, bac_res]) * lag_min
delta = np.zeros([2, bac_res, bac_res])

# competitor parameters
I = np.ones(bac_res)
lag[1] += np.outer(I, np.linspace(0, max(T0 + Tab), bac_res))
delta[1] += np.outer(np.linspace(0, delta_max, bac_res), I)

# transforming parameters to a-b space
a, b = a_b(lag, delta)
ap, bp = ap_bp(lag, delta)



########################
## Running Simulation ##
########################
if data == 'new':
    bac_args = [lag, delta, a, b, ap, bp]
    ab_args = [p_arr, T0, Tab]
    sim_args = [ab_res, bac_res, t_res, tot_cycles, repetitions, save_data]

    tic = time.time()
    if __name__ == '__main__':
        S_frac, lag_opt, del_opt = run_competition_in_parallel(bac_args, ab_args, sim_args)

    toc = time.time()
    print(f'\nTime: {(toc - tic) / 3600}h\n')

    # saving
    if save_data:
        np.savetxt("data/competition_opt/ext_comp_Sfrac-"+T_labels[ic]+str(T_const), S_frac)
        np.savetxt("data/competition_opt/ext_comp_lag-"+T_labels[ic]+str(T_const), lag_opt)
        np.savetxt("data/competition_opt/ext_comp_delta-"+T_labels[ic]+str(T_const), del_opt)



###################
## Plotting data ##
###################
if data == 'old':
    S_frac = np.loadtxt(path+'comp_Sfrac-'+severity+'-'+T_labels[ic]+str(T_const))
    lag_opt = np.loadtxt(path+'comp_lag-'+severity+'-'+T_labels[ic]+str(T_const))
    del_opt = np.loadtxt(path+'comp_delta-'+severity+'-'+T_labels[ic]+str(T_const))


T = np.outer(np.ones_like(T0), T0) + np.outer(np.ones_like(Tab), Tab)
T[0, 0] = 1
time_labels = [r'$T_0$', r'$T_{AB}$']
title = time_labels[ic] + ' = ' + str(T_const)

fig, ax = plt.subplots(1, 3, figsize=(15, 4), sharey=True)
# fig.suptitle(title)

# plotting lag time
ax[0].set(xlabel=time_labels[1-ic] + r'$~[h]$', ylabel=r"$p$")#, title=r"$\lambda^*$")
im0 = ax[0].imshow(lag_opt / T, origin="lower", cmap=lag_cmap, aspect="auto", vmin=0, extent=[0, T_arr.max(), 0, 1])
cbar = fig.colorbar(im0, ax=ax[0], aspect=20, anchor=(-.1, 0.5))
cbar.set_label(r"$\lambda^* / ~T$", rotation=0, labelpad=30)

# plotting delta
ax[1].set(xlabel=time_labels[1-ic] + r'$~[h]$', ylabel=r"$p$")#, title=r"$\delta^*$")
im1 = ax[1].imshow(del_opt, origin="lower", cmap=del_cmap, aspect="auto", vmin=0, vmax=0.075, extent=[0, T_arr.max(), 0, 1])
cbar = fig.colorbar(im1, ax=ax[1], aspect=20, anchor=(-.1, 0.5))
cbar.set_label(r"$\delta^* [h^{-1}]$", rotation=0, labelpad=50)

# ax[2].set(xlabel=time_labels[1-ic] + r'$~[h]$', ylabel=r"$p$")#, title="Consumption fraction")
im2 = ax[2].imshow(S_frac, origin="lower", aspect="auto", extent=[0, max(T_arr), 0, 1])
cbar = fig.colorbar(im2, ax=ax[2], aspect=20, anchor=(-.1, 0.5))
cbar.set_label("Consumption fraction", rotation=0, labelpad=50)

fig.tight_layout()
fig.show()

if save_fig:
    fig.savefig("opt_competitor-"+T_labels[ic]+str(T_const)+"totc30_reps50.png", dpi=100)

