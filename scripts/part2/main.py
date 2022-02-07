# Computing average parameters when evolving several species
import numpy as np
import matplotlib.pyplot as plt
import time

from model_equations import initialise_system, evolve_system, savedata


###########################
## Simulation parameters ##
###########################
# np.random.seed(2)
saving = True						# saving data
model = 1							# models: 1 (triggered), 2 (spontaneous, one lag time), 3 (spontaneous, two lag times)
mutation = True
extinction = False
mutation_seed = 'min'				# seeds: 'min' , 'max', 'optimal' parameters
tot_cycles = 10 #10_000				# total number of cycles / lenght of evolution



###########################
## Antibiotic parameters ##
###########################
T0 = 0								# application time
Tab = 6 							# duration
T = T0 + Tab 						# removal time
p_arr = np.linspace(0.1, 0.9, 2)	# application probability



################
## Simulation ##
################
# looping though probabilities
for p in p_arr:
    toc = time.time()

    # ramdom array for determining if applying ab
    r_arr = np.random.rand(tot_cycles)
    sim_args = [model, mutation, extinction, mutation_seed, r_arr]
    ab_args = [p, T0, Tab]

    n_t, bac_args = initialise_system(ab_args, sim_args)
    data = evolve_system(n_t, bac_args, ab_args, sim_args)
    
    tic = time.time()
    print(f'{tic-toc:2.4f} s')

    savedata(ab_args, sim_args, data, saving)

