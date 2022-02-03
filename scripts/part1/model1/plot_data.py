import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from functions import optimal_lag, analytical_fitness

# setting antibiotic parameters
T0 = 0							# application time
Tab = 6							# duration
p_range = np.array([0.2, 0.3, 0.6, 0.9]) 	# probability of ab application, in percentage
save_fig = True


####################################################
## Importing data and setting plotting parameters ##
####################################################
# importing data
tot_cycles = 20_000
data = []
for p in p_range:
    prob = str(int(100*p))
    # data.append(np.loadtxt("../../../data/model1/consumption_fraction-p"+prob+"-T"+str(Tab)+"-T0"+str(T0)))
    data.append(np.loadtxt("../../../data/model1/cycles_"+str(tot_cycles)+"/consumption_fraction-p"+prob+"-T"+str(Tab)+"-T0"+str(T0)))


# setting plotting parameters
font = {'family': 'Times New Roman',
        'weight': 'normal',
        'size': 20}

mpl.rc('font', **font)
fs_text = 14
cmap = mpl.cm.get_cmap('CMRmap')
colors = [cmap(i/len(p_range)) for i in range(len(p_range))]



##############
## Plotting ##
##############
fig, ax = plt.subplots(2, 1, figsize=(12, 6), sharex=True)
ax[0].grid(False)
ax[1].grid(False)
ax[0].set(ylabel=r"$F({\lambda}^*) - F({\lambda})$")
ax[1].set(xlabel=r"$\lambda$ [h]", ylabel="Consumption fraction")


# Ax 1: plotting analytical fitness
lag_arr = np.linspace(0, Tab, 1000) + 10**(-4)

i = 0
for p in p_range:
    ab_args = [p, T0, Tab]
    lag_opt, F_opt = optimal_lag(ab_args, 0.1)
    F_arr = analytical_fitness(lag_arr, ab_args)
    F_arr[lag_arr == lag_opt] += 10**(2)

    ax[0].plot(lag_arr, F_opt - F_arr, lw=3, color=colors[i], label=r"$p$ = "+str(p))
    i += 1

# ax[0].legend(loc="upper right")
ax[0].set(yscale="log")


# Ax 2: plotting numerical consumption fraction
lag_num = np.linspace(0, Tab, len(data[0])) + 10**(-4)

for i in range(len(p_range)):
	ax[1].plot(lag_num, data[i], 'x-', fillstyle='none', lw=3, ms=10, markeredgewidth=3, color=colors[i], label=r"$p$ = "+str(p_range[i]))
ax[1].set(xscale="linear")
ax[1].hlines(0.5, lag_arr.min(), lag_arr.max(), color="gray", ls="--", lw=2)

fig.tight_layout(rect=[0, 0, 0.84, 1])
handles, labels = ax[0].get_legend_handles_labels()
fig.legend(ax, labels=labels, bbox_to_anchor=[0.98, 0.7], borderaxespad=0.2)

# fig.subplots_adjust(right=0.85)
fig.show()
if save_fig:
	fig.savefig("../../../figures/model1/comparing_numerical_analytical_fitness.png", dpi=100)


