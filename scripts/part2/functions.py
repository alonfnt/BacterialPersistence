# Functions for computing optimal parameters from functions in model_equations.py
import numpy as np
from tqdm import tqdm
from scipy.integrate import solve_ivp

from model_equations import ODEs_growth, ODEs_decay, model_to_states
from model_equations import dlag, ddelta, domega, delta_max, lag_min, dt_max, T0_max, T_max, n_min
from model_equations import f, S0, n0


#######################
## Solving one cycle ##
#######################
def initialise_system(ab_args, sim_args):
	p, T0, Tab = ab_args
	model, mutation, _, mut_seed, _ = sim_args
	T = T0 + Tab
	n_states = model_to_states[model-1]		# number of states a cell can take
	eps = mutation_rate

	# determining size of system
	Nl = int((T) / dlag)					# number of lag times
	Nd = int(delta_max / ddelta)			# number of deltas
	No = int((T) / domega)					# number of omegas
	if model == 1: Nd = 1
	if model <= 2: No = 1
	N = Nl * Nd * No						# total number of species

	# defining parameter arrays
	lag = dlag * np.arange(0, Nl, 1) + lag_min
	lag = np.tile(np.tile(lag, Nd), No)
	delta = ddelta * np.arange(0, Nd, 1)
	delta = np.tile(np.repeat(delta, Nl), No)
	omega = domega * np.arange(0, No, 1) + lag_min
	omega = np.repeat(np.repeat(omega, Nl), Nd)

	# choosing whether to allow mutation and initialising initial conditions in parameter space
	if mutation == False:
		extinction = False					# no extinction in case without mutation
		eps = 0 							# setting mutation rate to 0
		initial_population = np.concatenate([n0 * np.ones(N), np.zeros((n_states - 1) * N), np.array([S0])])

	else:
		initial_population = np.concatenate([np.zeros(n_states * N), np.array([S0])])

		# initialising simulation from combination of smallest parameters
		if mut_seed == 'min':
			initial_population[0] = n0

		# initialising simulation from combination of biggest parameters
		elif mut_seed == 'max':
			initial_population[len(np.unique(lag)) - 1] = n0

		# initialising simulation from combination of optimal parameters, 
		# identifying optimal parameters from data
		elif mut_seed == 'optimal':
			assert(model > 1)
			ab_res = len(np.loadtxt("../../data/model2/low_resolution/optimal_lag-Tab" + str(int(T))))
			ip = int(p * ab_res)
			it = int(T0 * ab_res / T0_max)
			lag_opt = np.loadtxt("../../data/model2/low_resolution/optimal_lag-Tab" + str(int(T)))[ip, it]
			delta_opt = np.loadtxt("../../data/model2/low_resolution/optimal_delta-Tab" + str(int(T)))[ip, it]

			# corresponding position in parameter space
			ilag = np.where(abs(lag - lag_opt) == min(abs(lag - lag_opt)))[0][0]
			idelta = np.where(abs(delta - delta_opt) == min(abs(delta - delta_opt)))[0][0]

			index = idelta + ilag
			if model == 3:
				index += il * Nl * Nd
		
			initial_population[index] = n0

		else:
			print("Error: Don't recognize mutation seed.")

	# bacterial parameters
	bac_args = (lag, delta, omega, eps)

	return initial_population, bac_args


# function for detecting when system runs out of food, used in solve_ivp()
def famine(t, n_t, lag, delta, omega, eps): return n_t[-1]
famine.terminal = True
famine.direction = -1


def solve_cycle(n_t, r, bac_args, ab_args, sim_args):
	p, T0, Tab = ab_args
	model, _, extinction, _, _ = sim_args
	N = len(bac_args[0])
	n_states = model_to_states[model-1]

	# creating mask of none-extinct species before ab is applied
	populations = n_t[:N]
	for i in range(1, n_states):
		populations += n_t[i*N:(i+1)*N]
	init_alive = (populations >= n_min)

	# chosing functions corresponding to model equations
	func_growth = ODEs_growth[model-1]
	func_decay = ODEs_decay[model-1]

	T = T0 + Tab
	t_low = 0					# lower bound for integration
	c_ext = 0					# counting species that goes extinct during this cycle

	# cycles with antibiotics
	if r < p:
		# solving for period before antibiotics
		sol_cycle = solve_ivp(func_growth, [0, T0], n_t, args=bac_args, max_step=dt_max)
		n_t = np.array([sol_cycle.y[i][-1] for i in range(n_states * N + 1)])

		# solving for period during antibiotics
		sol_cycle = solve_ivp(func_decay, [T0, T], n_t, args=bac_args, max_step=dt_max)
		n_t = np.array([sol_cycle.y[i][-1] for i in range(n_states * N + 1)])

		# checking for extinctions (only for simulations with lower threshold implimented)
		if extinction == True:
			# computing total population of each species
			populations = n_t[:N]
		for i in range(1, n_states):
			populations += n_t[i*N:(i+1)*N]

			# checking for extinction
			for i in range(n_states):
				n_t[i*N:(i+1)*N][populations < n_min] = 0

			c_ext = sum((populations < n_min) * init_alive)		# updating counter
			init_alive = init_alive & (populations > n_min)		# updating mask to contain only species that are alive after ab

		t_low = T												# updating lower bound for integration

	# cycles without antibiotics / solving for period after antibiotics
	sol_cycle = solve_ivp(func_growth, [t_low, T_max], n_t, args=bac_args, events=famine, max_step=dt_max)
	d_t = f * np.array([sol_cycle.y[j][-1] for j in range(N)])

	# returning all species to sleep
	for i in range(1, n_states):
		d_t += f * np.array([sol_cycle.y[i*N + j][-1] for j in range(N)])

	# cheking for extinctions  (only for simulations with lower threshold implimented)
	if extinction == True:
		c_ext += sum((d_t < n_min) * init_alive)	# updating counter
		d_t[d_t < n_min] = 0 						# truncating extinct species

	return d_t, c_ext


def evolve_system(n_t, bac_args, ab_args, sim_args):
	# breaking if mutation rate = 0
	lag, delta, omega, _ = bac_args
	model = sim_args[0]

	# sorting parameters
	r_arr = sim_args[4]
	tot_cycles = len(r_arr)
	N = len(lag)
	n_states = model_to_states[model-1]

	# output arrays
	lag_avrg = np.zeros_like(r_arr)
	delta_avrg = np.zeros_like(r_arr)
	omega_avrg = np.zeros_like(r_arr)
	dominant_species = np.zeros_like(r_arr)
	extinctions = np.zeros_like(r_arr)
	pop_dist = np.zeros([len(r_arr), N])

	# looping through cycles
	for ic in tqdm(range(tot_cycles)):
		d_t, ext = solve_cycle(n_t, r_arr[ic], bac_args, ab_args, sim_args)
		# if all species go extinct
		if sum(d_t) == 0:
			break

		# preparing for new cycle
		g_t = np.zeros((n_states - 1) * N)
		n_t = np.concatenate([d_t, g_t, np.array([S0])])

		lag_avrg[ic] = sum(lag * d_t) / sum(d_t)
		delta_avrg[ic] = sum(delta * d_t) / sum(d_t)
		omega_avrg[ic] = sum(omega * d_t) / sum(d_t)
		dominant_species[ic] = np.where(d_t == d_t.max())[0][0]
		extinctions[ic] = ext
		pop_dist[ic] = d_t / sum(d_t)

	return lag_avrg, delta_avrg, omega_avrg, dominant_species, extinctions, pop_dist


def savedata(ab_args, sim_args, data, saving):
	p, T0, Tab = ab_args
	model, mutation, extinction, mutation_seed, _ = sim_args
	T = T0 + Tab

	if saving:
		path = "../../data/model"+str(model)
		parameters = "p" + str(int(p*100)) + "-T0" + str(T0) + "-T" + str(T)

		if mutation:
			path += "/mutation" + extinction * "_extinction" + "/mutation_rate-"+'{:.0e}'.format(float(mutation_rate))
			parameters += "-" + mutation_seed

		else:
			path += "/competition_all_species"

		np.savetxt(path + "/average_lag-" + parameters, data[0])
		# np.savetxt(path + "/dominant_lag-" + parameters, lag_dom)

		if model >= 2:
			np.savetxt(path + "/average_delta-" + parameters, data[1])
			# np.savetxt(path + "/dominant_delta-" + parameters, delta_dom)

		if model == 3:
			np.savetxt(path + "/average_omega-" + parameters, data[2])
			# np.savetxt(path + "/dominant_omega-" + parameters, omega_dom)


		if extinction == True:
			np.savetxt(path + "/extinctions-" + parameters, data[4])
