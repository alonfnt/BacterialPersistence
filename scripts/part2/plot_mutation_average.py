# plotting long term average from population average from 'main.py' as function of p
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

font = {'family': 'Times New Roman',
        'weight': 'normal',
        'size': 20}
mpl.rc('font', **font)

saving = True                       # saving figure
model = 2                           # models: 1 (triggered), 2 (spontaneous, one lag time), 3 (spontaneous, two lag times)
mutation = True
extinction = False
mutation_seed = 'optimal' 			# seeds: 'min' , 'max', 'optimal' parameters
mutation_rate = 0.001
transient = 1000                    # index of transient state


###########################
## Antibiotic parameters ##
###########################
T0 = 0 									# application time
Tab = np.array([6.0])					# duration
T = T0 + Tab							# removal time
p_arr = [0.1, 0.2, 0.3, 0.5, 0.7, 0.9] 	# probabilities



####################
## Defining paths ##
####################
data_path1 = "../../data/model" + str(model) 
data_path2 = data_path1 + mutation * ("/mutation" + extinction * "_extinction" + "/mutation_rate-"+'{:.0e}'.format(float(mutation_rate)))
data_path2 += (not mutation) * "/competition_all_species"

fig_path = "../../figures/model" + str(model)
time_parameters = "-T0" + str(T0) + "-T" + str(T) + mutation * ("-" + mutation_seed)
mutation_parameters = mutation * ("_mutation" + extinction * "_extinction")

ab_res = len(np.loadtxt(data_path1 + "/high_resolution/optimal_lag-T0" + str((T0))))
t_index = np.array([int(ab_res * t / 24) for t in Tab])
Nt = len(Tab)


####################
## Importing data ##
####################
T = (T0 + Tab[0])

lag_avrg, delta_avrg = [], []       # average values from part 2
lag_opt, delta_opt = [], []         # optimal values from part 1

for p in p_arr:
    parameters = "p" + str(int(p*100)) + time_parameters
    lag_avrg.append(np.loadtxt(data_path2 + "/average_lag-" + parameters))
    delta_avrg.append(np.loadtxt(data_path2 + "/average_delta-" + parameters))

i = 0
for t in Tab:
    lag_opt.append(np.loadtxt(data_path1 + "/high_resolution/optimal_lag-T0" + str((T0)))[:, t_index[i % Nt]])
    delta_opt.append(np.loadtxt(data_path1 + "/high_resolution/optimal_delta-T0" + str((T0)))[:, t_index[i % Nt]])
    i += 1

lag_p = np.zeros_like(p_arr)
lag_std = np.zeros_like(p_arr)
delta_p = np.zeros_like(p_arr)
delta_std = np.zeros([2, len(p_arr)])

for i in range(len(p_arr)):
    lag_p[i] = np.mean(lag_avrg[i][transient:])
    lag_std[i] = np.std(lag_avrg[i][transient:])
    delta_p[i] = np.mean(delta_avrg[i][transient:])
    delta_std[0, i] = np.std(delta_avrg[i][transient:])
    delta_std[1, i] = np.std(delta_avrg[i][transient:])
    if delta_p[i] - delta_std[0,i] < 0:
        delta_std[0,i] = delta_p[i]


#######################
## PREPARING TO PLOT ##
#######################
p_opt = np.linspace(0, 1, ab_res)
color = ['b']

fig, ax = plt.subplots(1, 2, figsize=(12, 5), sharex=True)
ax[0].set(xlabel=r"$p$", ylabel=r"$\langle \lambda \rangle_c$")
ax[1].set(xlabel=r"$p$", ylabel=r"$\langle \delta \rangle_c$")

for i in range(len(Tab)):
    mask0 = (lag_opt[i] < 0.1)
    mask1 = (lag_opt[i] > 0.1)

    ax[0].plot(p_opt[mask0], lag_opt[i][mask0], '-', color=color[i], fillstyle='none', alpha=0.4, lw=3, markeredgewidth=3)
    ax[0].plot(p_opt[mask1], lag_opt[i][mask1], '-', color=color[i], fillstyle='none', alpha=0.4, lw=3, markeredgewidth=3)
    ax[0].errorbar(p_arr, lag_p, yerr=lag_std, fmt='o', color=color[i], fillstyle='full', elinewidth=1, capsize=3, label=r"$T_0$ = " + str(T0) + r",  $T_{AB}$ = " + str(T-T0))

    ax[1].plot(p_opt[mask0], delta_opt[i][mask0], '-', color=color[i], alpha=0.5, lw=3)
    ax[1].plot(p_opt[mask1], delta_opt[i][mask1], '-', color=color[i], alpha=0.5, lw=3)
    ax[1].errorbar(p_arr, delta_p, yerr=delta_std, fmt='o', color=color[i], fillstyle='full', elinewidth=1, capsize=3, label=r"$T_0$ = " + str(T0) + r",  $T_{AB}$ = " + str(T-T0))

fig.tight_layout(rect=[0, 0, 1, 0.85])
ax[1].legend(loc='upper center', bbox_to_anchor=(-0.18, 1.3), ncol=2, fancybox=True, shadow=False)
fig.show()
if saving:
    fig.savefig("../../figures/model2/mutation_average"+time_parameters+".png")
