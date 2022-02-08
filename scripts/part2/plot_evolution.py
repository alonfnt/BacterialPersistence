# plotting evolution of average parameters
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# setting parameters..
font = {'family': 'Times New Roman',
		'weight': 'normal',
		'size': 20}
# plt.rcParams.update({"text.usetex": True})
mpl.rc('font', **font)
cmap = mpl.cm.get_cmap('CMRmap')


saving = True                       # saving figure
model = 3                           # models: 1 (triggered), 2 (spontaneous, one lag time), 3 (spontaneous, two lag times)
mutation = True
extinction = False
mutation_seed = 'min'               # seeds: 'min' , 'max', 'optimal' parameters
mutation_rate = 0.001


###########################
## Antibiotic parameters ##
###########################
T0 = 0                              # application time
Tab = 6                             # duration
T = T0 + Tab                        # removal time
p_arr = [0.1, 0.3, 0.5, 0.7, 0.9]   # probabilities


####################
## Defining paths ##
####################
data_path = "../../data/model" + str(model) 
data_path += mutation * ("/mutation" + extinction * "_extinction" + "/mutation_rate-"+'{:.0e}'.format(float(mutation_rate)))
data_path += (not mutation) * "/competition_all_species"

fig_path = "../../figures/model" + str(model)
time_parameters = "-T0" + str(T0) + "-T" + str(T) + mutation * ("-" + mutation_seed)
mutation_parameters = mutation * ("_mutation" + extinction * "_extinction")


####################
## Importing data ##
####################
colors = [cmap(i/len(p_arr)) for i in range(len(p_arr))]

lag_avrg, lag_dom = [], []
delta_avrg, delta_dom = [], []
omega_avrg, omega_dom = [], []
extinctions = []

for p in p_arr:
	parameters = "p" + str(int(p*100)) + time_parameters

	lag_avrg.append(np.loadtxt(data_path + "/average_lag-" + parameters))
	# lag_dom.append(np.loadtxt(data_path + "/dominant_lag-" + parameters))

	if model >= 2:
		delta_avrg.append(np.loadtxt(data_path + "/average_delta-" + parameters))
		# delta_dom.append(np.loadtxt(data_path + "/dominant_delta-" + parameters))

	if model == 3:
		omega_avrg.append(np.loadtxt(data_path + "/average_omega-" + parameters))
		# omega_dom.append(np.loadtxt(data_path + "/dominant_omega-" + parameters))


	# if extinction:
	#	extinctions.append(np.loadtxt(data_path + "/extinctions-" + parameters))

tot_cycles = len(lag_avrg[0])
cycles = np.arange(1, tot_cycles+1)



##############
## Plotting ##
##############
title = "Model " + str(model) + mutation * " with mutation" + extinction * " and extinction"
colors = ['purple', 'red', 'orange', 'yellowgreen', 'green', 'cyan', 'blue', 'black', 'magenta']


# plotting average lag time
if model == 1:
	fig, ax = plt.subplots(1, 1, figsize=(7, 6))
	ax.set(xlabel='Cycle number', ylabel='Average lag time [h]', title=title)
	# ax[1].set(xlabel='Cycle number', ylabel='Dominant lag time [h]', title=title)

	i = 0
	for pi in p_arr:
		ax.plot(cycles, lag_avrg[i], color=colors[i], label="p = " + str(pi))
		#ax[1].plot(cycles, lag_dom[i], color=colors[i], label="p = " + str(pi))
		i += 1


	ax.set(yscale='log', ylim=(10**(-2), 20))
	fig.tight_layout(rect=[0, 0, 1, 0.85])
	ax.legend(loc='upper center', bbox_to_anchor=(0.45, 1.45),
 			ncol=3, fancybox=True, shadow=False)

	fig.show()
	if saving:
		fig.savefig(fig_path + "/evolution_average" + mutation_parameters + time_parameters)


elif model == 2:
    # plotting average lag and delta
    fig, ax = plt.subplots(1, 2, figsize=(12, 6))
    fig.suptitle(title)
    ax[0].set(xlabel='Cycle number', ylabel=r'$\langle\lambda\rangle$ [h]') #, title=title)
    ax[1].set(xlabel='Cycle number', ylabel=r'$\langle\delta\rangle$ [1/h]') #, title=title)

    i = 0
    for pi in p_arr:
        ax[0].plot(cycles, lag_avrg[i], color=colors[i], label="p = "+str(pi))
        # ax[0].hlines(pi * T, 0, 10000, ls="dashed", color=colors[i])
        ax[1].plot(cycles, delta_avrg[i], color=colors[i], label="p = "+str(pi))
        i += 1

    ax[0].set(yscale='log', ylim=(0.009, 20))
    ax[1].set(yscale='log', ylim=(0,0.5))
    fig.tight_layout(rect=[0, 0, 1, 0.9])
    ax[1].legend(loc='upper center', bbox_to_anchor=(-0.15, 1.35),
          ncol=3, fancybox=True, shadow=False)
    fig.show()
    if saving:
        fig.savefig(fig_path + "/evolution_average" + mutation_parameters + time_parameters)

    # plotting dominant lag and delta
    # fig, ax = plt.subplots(1, 2, figsize=(12, 4))
    # ax[0].set(xlabel='Cycle number', ylabel='Dominant lag time [h]', title=title)
    # ax[1].set(xlabel='Cycle number', ylabel='Dominant delta time [h]', title=title)

    # i = 0
    # for pi in p_arr:
    #     ax[0].plot(cycles, lag_dom[i], color=colors[i], label="p = "+str(pi))
    #     ax[1].plot(cycles, delta_dom[i], color=colors[i], label="p = "+str(pi))
    #     i += 1

    # ax[0].set(yscale='log', ylim=(0.009, 10))
    # ax[1].set(yscale='linear')
    # ax[0].legend()
    # fig.tight_layout()
    # # fig.show()
    # if saving:
    #     fig.savefig(fig_path + "/evolution_average_" + mutation_parameters + time_parameters)



elif model == 3:
    xTicks = ['0', r'$10^3$']
    # plotting average lag and delta
    fig, ax = plt.subplots(1, 3, figsize=(12, 4.5))
    fig.suptitle(title)
    ax[0].set(xlabel='Cycle number', ylabel=r'$\langle\lambda\rangle$') #, title=title)
    ax[2].set(xlabel='Cycle number', ylabel=r'$\langle\delta\rangle$') #, title=title)
    ax[1].set(xlabel='Cycle number', ylabel=r'$\langle\omega\rangle$') #, title=title)

    i = 0
    for pi in p_arr:
        ax[0].plot(cycles, lag_avrg[i], color=colors[i], label="p = "+str(pi))
        ax[2].plot(cycles, delta_avrg[i], color=colors[i], label="p = "+str(pi))
        ax[1].plot(cycles, omega_avrg[i], color=colors[i], label="p = "+str(pi))
        i += 1

    ax[0].set(yscale='log', ylim=(0.01, 20))
    ax[1].set(yscale='log', ylim=(0.01, 20))
    ax[2].set(yscale='log', ylim=(5*10**(-4), 0.2))
    # ax[0].legend()
    # fig.tight_layout()
    fig.tight_layout(rect=[0, 0, 1, 0.9])
    fig.subplots_adjust(wspace=0.4, left=0.08, right=0.98)
    ax[1].legend(loc='upper center', bbox_to_anchor=(0.5, 1.35),
                 ncol=5, fancybox=True, shadow=False)
    fig.show()
    if saving:
        fig.savefig(fig_path + "/evolution_average" + mutation_parameters + time_parameters)

    # plotting dominant lag and delta
    # fig, ax = plt.subplots(1, 3, figsize=(12, 4))
    # ax[0].set(xlabel='Cycle number', ylabel='Dominant lag time [h]', title=title)
    # ax[1].set(xlabel='Cycle number', ylabel='Dominant delta', title=title)
    # ax[2].set(xlabel='Cycle number', ylabel='Dominant omega [h]', title=title)

    # i = 0
    # for pi in p_arr:
    #     ax[0].plot(cycles, lag_dom[i], color=colors[i], label="p = "+str(pi))
    #     ax[1].plot(cycles, delta_dom[i], color=colors[i], label="p = "+str(pi))
    #     ax[2].plot(cycles, omega_dom[i], color=colors[i], label="p = "+str(pi))
    #     i += 1

    # ax[0].set(yscale='log', ylim=(0.009, 10))
    # ax[1].set(yscale='linear')
    # ax[0].legend()
    # fig.tight_layout()
    # fig.show()
    # if saving:
    #     fig.savefig(fig_path + "/evolution_average" + mutation_parameters + time_parameters)



# if extinction:
#     fig, ax = plt.subplots(1, 1, figsize=(5, 4))
#     ax.set(xlabel='Cycle number', ylabel='# Extinctions', title=title)

#     i = 0
#     for pi in p_arr:
#         ax.plot(cycles, extinctions[i], color=colors[i], label="p = " + str(pi))
#         i += 1

#     ax.set(xscale='log', yscale='linear')
#     ax.legend()
#     fig.tight_layout()
#     #fig.show()
#     if saving:
#         fig.savefig(fig_path + "/extinctions_" + time_parameters)

