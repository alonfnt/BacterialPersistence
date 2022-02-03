import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

font = {'family': 'Times New Roman',
        'weight': 'normal',
        'size': 20}
mpl.rc('font', **font)

folder = "single_plot"

T_labels = ['T0', 'Tab']
path = 'data/' + folder
ab_res = 400


##############
## SETTINGS ##
##############
p = np.array([0.4, 0.5, 0.6])
T0 = np.array(([3.6, 3.7, 3.8]))
Tab = np.array([10])
Nt0 = len(T0)
t0_index = np.array([int(ab_res * t0 / 12) for t0 in T0])

Nt = len(Tab)
t_index = np.array([int(ab_res * t / 24) for t in Tab])
ts_index = np.array([int(100 * t / 24) for t in Tab])

Np = len(p)
p_index = np.array([int(ab_res * pi)-1 for pi in p])
ps_index = np.array([int(100 * pi) for pi in p])


#######################
## PREPARING TO PLOT ##
#######################
lag_opt_p = np.zeros([Nt*Nt0, 400])
del_opt_p = np.zeros([Nt*Nt0, 400])

lag_opt_t = np.zeros([Np*Nt0, 400])
del_opt_t = np.zeros([Np*Nt0, 400])

lag_com_p = np.zeros([Nt*Nt0, ab_res])
del_com_p = np.zeros([Nt*Nt0, ab_res])

lag_com_t = np.zeros([Np*Nt0, ab_res])
del_com_t = np.zeros([Np*Nt0, ab_res])

i = 0
for t0 in T0:
    for t in Tab:
        lag_opt_p[i] = np.loadtxt('data/single_plot/single_lag-Tab' + str(t))[:, t0_index[i % Nt0]]
        del_opt_p[i] = np.loadtxt('data/single_plot/single_delta-Tab' + str(t))[:, t0_index[i % Nt0]]
        # lag_com_p[i] = np.loadtxt(path + 'comp_lag-all-T0' + str(t0))[:, t_index[i % Nt]]
        # del_com_p[i] = np.loadtxt(path + 'comp_delta-all-T0' + str(t0))[:, t_index[i % Nt]]
        i += 1

i = 0
for t in Tab:
    for pi in p:
        lag_opt_t[i] = np.loadtxt('data/single_plot/single_lag-Tab' + str(t))[p_index[i % Np]]
        del_opt_t[i] = np.loadtxt('data/single_plot/single_delta-Tab' + str(t))[p_index[i % Np]]
        # lag_com_t[i] = np.loadtxt(path + 'comp_lag-all-T0' + str(t0))[p_index[i % Np]]
        # del_com_t[i] = np.loadtxt(path + 'comp_delta-all-T0' + str(t0))[p_index[i % Np]]
        i += 1

p_arr = np.linspace(0, 1, ab_res)
T_arr = np.linspace(0, 24, ab_res)
color_T = ['dodgerblue', 'blue', 'black']
color_p = ['orange', 'red', 'maroon', 'blue']

p_single = np.linspace(0, 1, 400)
T_single = np.linspace(0, 24, 400)

#######################
## Single lines of p ##
# #######################
fig, ax = plt.subplots(1, 2, figsize=(12, 5.5), sharex=True)
ax[0].set(xlabel=r"$T_{AB}$", ylabel=r"$\lambda^*$")
ax[1].set(xlabel=r"$T_{AB}$", ylabel=r"$\delta^*$", ylim=(0, 0.075))
fig.suptitle(r"Optimal persistence strategy for $T_0$ = " + str(T0[0]))

# plotting
for i in range(Np*Nt):
    # ax[0].plot(T_arr, lag_com_t[i], 'o', mew=1.5, alpha=1, color=color_T[i], fillstyle='none',
    #            label=r"$p$ = " + str(p[i % Np]))
    # ax[1].plot(T_arr, del_com_t[i], 'o', mew=1.5, alpha=1, color=color_T[i], fillstyle='none',
    #            label=r"$p$ = " + str(p[i % Np]))

    mask0 = (lag_opt_t[i] < 1)
    mask1 = (lag_opt_t[i] > 0.1)
    ax[0].plot(T_single[mask0], lag_opt_t[i][mask0], '.', alpha=0.7, lw=4, color=color_T[i], fillstyle='none') #,  label=r"$p$ = " + str(p[i % Np]))  # +", T0=" + str(T0[i]))
    ax[0].plot(T_single[mask1], lag_opt_t[i][mask1], '.', alpha=0.7, lw=4, color=color_T[i], fillstyle='none')
    # ax[0].plot(T_single, lag_opt_t[i], '-', color=color_T[i], fillstyle='none')

    ax[1].plot(T_single[mask0], del_opt_t[i][mask0], '.', alpha=0.7, lw=4, color=color_T[i], fillstyle='none') #,  label=r"$p$ = " + str(p[i % Np]))  #+", T0=" + str(T0[i]))
    ax[1].plot(T_single[mask1], del_opt_t[i][mask1], '.', alpha=0.7, lw=4, color=color_T[i], fillstyle='none')
    #ax[1].plot(T_single, del_opt_t[i], '-', color=color_T[i], fillstyle='none')


fig.tight_layout(rect=[0, 0, 1, 0.9])
ax[0].legend(loc='upper center', bbox_to_anchor=(1.08, 1.3), ncol=3, fancybox=True, shadow=False)
ax[0].set(xscale="linear")

# saving
fig.show()
# # fig.savefig("figs/single_lags_p_T0" + str(T0[0]) + "_appendix.png", dpi=100)
# # fig.savefig("figs/compete_lags_p_T0" + str(T0[0]) + ".png", dpi=100)
#

######################
# Single lines of T ##
######################
fig, ax = plt.subplots(1, 2, figsize=(14, 4), sharex=True)
ax[0].set(xlabel=r"$p$", ylabel=r"$\lambda^*$")
ax[1].set(xlabel=r"$p$", ylabel=r"$\delta^*$", ylim=(0, 0.075))
# fig.suptitle(r"Optimal persistence strategy for $T_0$ = " + str(T0[0]))

# plotting
for i in range(Nt*Nt0):
    # ax[0].plot(p_arr, lag_com_p[i], 'o', mew=1.5, alpha=1, color=color_p[i], fillstyle='none', label=r"$T_{AB}$ = " + str(int(Tab[i % Nt])))
    # ax[1].plot(p_arr, del_com_p[i], 'o', mew=1.5, alpha=1, color=color_p[i], fillstyle='none', label=r"$T_{AB}$ = " + str(int(Tab[i % Nt])))

    mask0 = (lag_opt_p[i] < 1)
    mask1 = (lag_opt_p[i] > 0.1) * (lag_opt_p[i] < 10)
    mask2 = (lag_opt_p[i] > 10)
    ax[0].plot(p_single, lag_opt_p[i], '--', alpha=0.5, lw=4, color=color_p[i], fillstyle='none')
    ax[0].plot(p_single[mask0], lag_opt_p[i][mask0], '-', alpha=1, lw=4, color=color_p[i], fillstyle='none', label=r"$T_{0}$ = " + str(np.round(T0[i % Nt0],2)))
    ax[0].plot(p_single[mask1], lag_opt_p[i][mask1], '-', alpha=1, lw=4, color=color_p[i], fillstyle='none')
    ax[0].plot(p_single[mask2], lag_opt_p[i][mask2], '-', alpha=1, lw=4, color=color_p[i], fillstyle='none')

    # ax[0].plot(T_single, lag_opt_t[i], '-', color=color_T[i], fillstyle='none')
    ax[1].plot(p_single, del_opt_p[i], '--', alpha=0.5, lw=4, color=color_p[i], fillstyle='none')
    ax[1].plot(p_single[mask0], del_opt_p[i][mask0], '-', alpha=1, lw=4, color=color_p[i], fillstyle='none', label=r"$T_{0}$ = " + str(np.round(T0[i % Nt0],2)))
    ax[1].plot(p_single[mask1], del_opt_p[i][mask1], '-', alpha=1, lw=4, color=color_p[i], fillstyle='none')
    ax[1].plot(p_single[mask2], del_opt_p[i][mask2], '-', alpha=1, lw=4, color=color_p[i], fillstyle='none')
    #ax[1].plot(T_single, del_opt_t[i], '-', color=color_T[i], fillstyle='none')

#fig.tight_layout(rect=[0, 0, 1, 0.9])
#ax[0].legend(loc='upper center', bbox_to_anchor=(1.08, 1.3),
#          ncol=3, fancybox=True, shadow=False)
ax[1].legend(bbox_to_anchor=(1.02, 0.8))
#          ncol=3, fancybox=True, shadow=False)
ax[0].set(xscale="linear")
fig.tight_layout()
# saving
fig.show()
fig.savefig("/home/silja/Documents/thesis/presentations/lineplot.pdf")
# fig.savefig("figs/single_lags_Tab_T0" + str(T0[0]) + "_appendix.png", dpi=100)
# fig.savefig("figs/compete_lags_Tab_T0" + str(T0[0]) + ".png", dpi=100)


# #####################
# ## Synchronisation ##
# ######################
# # setting up figure
# fig, ax = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
# fig.suptitle(r"Optimal persistence strategy for $T_0$ = 0")
# ax[0].set(xlabel=r"$p}$", ylabel=r"$\lambda^*$")
# ax[1].set(xlabel=r"$T_{AB}$")
#
# # plotting
# for i in range(3):
#     ax[0].plot(p_arr, lag_opt_p[i], 'o', color=color_p[i], fillstyle='full', label=r"$T_{AB}$ = " + str(np.round(Tab[i % Nt], 2)))  # +",
#     ax[1].plot(T_arr, lag_opt_t[i], 'o', color=color_T[i], fillstyle='full', label=r"$p$ = " + str(p[i % Np]))  # +",
#
# ax[0].legend()# bbox_to_anchor=(1, 1), loc="best")
# ax[1].legend()# bbox_to_anchor=(1, 1), loc="best")
#
# ax[0].set(xscale="log")
# ax[1].set(xscale="log")
#
# # saving
# fig.tight_layout()
# fig.show()
# fig.savefig("figs/single_lags_p_Tab_T00", dpi=100)
