# importing libraries
import numpy as np
import matplotlib.pyplot as plt
import time
from tqdm import tqdm

# importing from other files
from functions import compete, plot_cycles
from model_equations import lag_min

# set some basic plotting parameters
plt.rcParams['font.size'] = 12


###########################
## Simulation parameters ##
###########################
save_data = False
save_fig = False

tot_cycles = 100    # total number of cycles
dlag = 0.1			# space between lag times
reps = 1            # repetitions for ensemble average


###########################
## Antibiotic parameters ##
###########################
T0 = 0				# application time
Tab = 6				# duration
p = 0.9				# application probability

prob = str(int(100*p))
ab_args = [p, T0, Tab]

# array of lag times
lag_arr = np.arange(0, T0 + Tab, dlag) + lag_min



################
## Simulation ##
################
sim_args = [tot_cycles, dlag]
comp_frac = np.zeros_like(lag_arr)
optimal_comp_frac = np.zeros_like(lag_arr)

fig, ax = plt.subplots(1, 1, figsize=(5, 4))

for rep in range(reps):
    tic = time.time()
    s_w, s_c = compete(lag_arr, ab_args, sim_args)
    toc = time.time()
    print(f"{toc - tic:2.4f} s")

    optimal_comp_frac += s_w / reps
    comp_frac += s_c / reps

if save_data:
    np.savetxt("../../../data/model1/optimal_consumption_fraction-p"+prob+"-T"+str(Tab)+"-T0"+str(T0), optimal_comp_frac)
    np.savetxt("../../../data/model1/consumption_fraction-p"+prob+"-T"+str(Tab)+"-T0"+str(T0), comp_frac)


ax.set(xlabel=r"$\lambda$ [h]", ylabel="Consumption fraction")
ax.plot(lag_arr, comp_frac, '.', label=r"lag time = $\lambda$")
ax.plot(lag_arr, optimal_comp_frac, '.', label=r"lag time = $\lambda^*$")
ax.legend()

fig.tight_layout()
fig.show()

if save_fig:
        fig.savefig("../../../figures/model1/consumption_fraction-p"+prob+"-T"+str(Tab)+"-T0"+str(T0)+".png")


# plotting a few cycles
# lag = np.array([0.1, 1])
# r_arr = np.random.rand(10)
# plot_cycles(lag, p, r_arr, Tab, False, savefig=True)