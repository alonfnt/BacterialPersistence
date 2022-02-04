import numpy as np
import matplotlib.pyplot as plt
import time

from model_equations import initialise_system, evolve_system


###########################
## Simulation parameters ##
###########################
# np.random.seed(2)
saving = True				# saving data

model = 1                   # models: 1 (triggered), 2 (spontaneous, one lag time), 3 (spontaneous, two lag times)
mutation = True
extinction = False
mutation_seed = 'min'		# seeds: min , max, optimal parameters
tot_cycles = 100 #10_000 		# total number of cycles / lenght of evolution

data_path = "model_" + str(model) + "/case_" + case + "/data/seed_" + mut_seed



###########################
## Antibiotic parameters ##
###########################
T0 = 6.0				# application time
Tab = 10.0				# duration
T = T0 + Tab				# removal time
p_arr = np.linspace(0.1, 0.9, 5) 	# application probability



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

    savedata(ab_args, sim_args, saving)





