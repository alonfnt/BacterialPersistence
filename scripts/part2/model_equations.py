import numpy as np
from tqdm import tqdm
from scipy.integrate import solve_ivp


##########################
## Intrinsic parameters ##
##########################
K = 10**9           		# carrying capacity
beta_max = 1        		# max growth rate
gamma = 1           		# kill rate
eps = 10**(-3)      		# mutation rate


#########################
## External parameters ##
#########################
S0 = K              		# initial substrate
n0 = 1000           		# initial population
f = 10**(-6)        		# dilution fraction
n_min = 1           		# lower threshold
T_max = 50          		# max length of cycle


###########################
## Simulation parameters ##
###########################
dt_max = 1          		# max time step
T0_max = 12			# upper limit on application time
lag_min = 10**(-2)          	# lower limit on lag time
delta_max = 0.1			# upper limit on delta

# resolution in parameter space
dlag = 0.25                  	# wake up rate from dormancy
ddelta = 0.01               	# rate of spontaneous persistence
domega = 0.2			# wake up rate from persistent state



############################
## Differential equations ##
############################
# Growth rate as function of nutrients
def beta(s):
    # return beta_max * (s / (s + K))				# Monod
    return beta_max						# Constant



# Two state model of triggered persistence
def model1_grow(t, n_t, lag, delta, omega, eps):
    n_t[:-1][n_t[:-1] < 0] = 10**(-200)     			# avoid negative populations
    N = len(np.unique(lag))                 			# number of species
    nb = 2                                 			# number of nearest neighbours in parameter space

    dn_dt = np.zeros_like(n_t)             			# empty array for computing derivatives
    g_t = n_t[N:2*N]                        			# isolating growing populations

    if n_t[-1] > 0:
        dn_dt[:N] = -n_t[:N] / lag                                # updating dormant populations
        dn_dt[N:2*N] = -dn_dt[:N] + (1 - nb * eps) * n_t[N:2*N]   # updating growing populations

        # mutation between nearest neightbours in parameter space
        g_pad = np.pad(g_t, 1, mode='edge')
        lag_mut = (np.roll(g_pad, 1) + np.roll(g_pad, -1))[1:N + 1]

        dn_dt[N:2*N] += eps * lag_mut.flatten()			# updating growing populations
        dn_dt[-1] = -dn_dt[N:2*N].sum()                           # updating nutrients

    return dn_dt


def model1_kill(t, n_t, lag, delta, omega, eps):
    n_t[:-1][n_t[:-1] < 0] = 10 ** (-200)               		# avoid negative populations
    N = len(lag)						# number of species

    dn_dt = np.zeros_like(n_t)					# empty array for computing derivatives
    dn_dt[:N] = -n_t[:N] / lag                          		# updating dormant populations
    dn_dt[N:2*N] = -dn_dt[:N] - gamma * n_t[N:2*N]      		# updating growing populations

    return dn_dt


# 2 state model of triggered and spontaneous persistence
def model2_grow(t, n_t, lag, delta, omega, eps):
    n_t[:-1][n_t[:-1] < 0] = 10**(-200)     			# avoid negative populations
    
    Nl = len(np.unique(lag))                			# number of lag times
    Nd = len(np.unique(delta))             			# number of deltas
    N = Nl * Nd                             			# number of species
    nb = 4                                  			# number of nearest neighbours in parameter space

    dn_dt = np.zeros_like(n_t)              			# empty array for computing derivatives
    g_t = n_t[N:2*N].reshape(Nd, Nl)        			# isolating growing populations

    if n_t[-1] > 0:
        dn_dt[:N] = -n_t[:N] / lag + delta * n_t[N:2*N]           # updating dormant populations
        dn_dt[N:2*N] = -dn_dt[:N] + (1 - nb * eps) * n_t[N:2*N]   # updating growing populations

        # mutation between nearest neightbours in parameter space
        g_pad = np.pad(g_t, 1, mode='edge')
        lag_mut = (np.roll(g_pad, 1, axis=1) + np.roll(g_pad, -1, axis=1))[1:Nd+1, 1:Nl+1]
        delta_mut = (np.roll(g_pad, 1, axis=0) + np.roll(g_pad, -1, axis=0))[1:Nd+1, 1:Nl+1]

        dn_dt[N:2*N] += eps * (lag_mut + delta_mut).flatten()  # updating growing populations
        dn_dt[-1] = -dn_dt[N:2*N].sum()                           # updating nutrients

    return dn_dt


def model2_kill(t, n_t, lag, delta, omega, eps):
    n_t[:-1][n_t[:-1] < 0] = 10 ** (-200)               		# avoid negative populations
    
    Nl = len(np.unique(lag))                			# number of lag times
    Nd = len(np.unique(delta))             			# number of deltas
    N = Nl * Nd                             			# number of species

    dn_dt = np.zeros_like(n_t)					# empty array for computing derivatives
    dn_dt[:N] = -n_t[:N] / lag + delta * n_t[N:2*N]     		# updating dormant populations
    dn_dt[N:2*N] = -dn_dt[:N] - gamma * n_t[N:2*N]      		# updating growing populations

    return dn_dt


# 3 state model of triggered and spontaneous persistence
def model3_grow(t, n_t, lag, delta, omega, eps):
    n_t[:-1][n_t[:-1] < 0] = 10**(-200)     			# avoid negative populations

    Nl = len(np.unique(lag))                			# number of lag times
    Nd = len(np.unique(delta))              			# number of re-entering dormancy rates
    No = len(np.unique(omega))					# number of wake-up rates
    N = Nl * Nd * No                        			# number of species
    nb = 6                                  			# number of nearest neighbours in parameter space

    dn_dt = np.zeros_like(n_t)              			# empty array for computing derivatives
    g_t = n_t[N:2*N].reshape(Nd, Nl, No)    			# isolating growing populations

    if n_t[-1] > 0:
        dn_dt[:N] = -n_t[:N] / lag                                				# updating dormant populations
        dn_dt[N:2*N] = -dn_dt[:N] + (1 - delta - nb * eps) * n_t[N:2*N] + n_t[2*N:3*N] / omega     # updating growing populations
        dn_dt[2*N:3*N] = delta * n_t[N:2*N] - n_t[2*N:3*N] / omega  				# updating persistent population

        # mutation between nearest neighbours in parameters space
        g_pad = np.pad(g_t, 1, mode='edge')
        lag_mut = (np.roll(g_pad, 1, axis=1) + np.roll(g_pad, -1, axis=1))[1:Nd+1, 1:Nl+1, 1:No+1]
        delta_mut = (np.roll(g_pad, 1, axis=0) + np.roll(g_pad, -1, axis=0))[1:Nd+1, 1:Nl+1, 1:No+1]
        omega_mut = (np.roll(g_pad, 1, axis=2) + np.roll(g_pad, -1, axis=2))[1:Nd+1, 1:Nl+1, 1:No+1]

        dn_dt[N:2*N] += eps * (lag_mut + delta_mut + omega_mut).flatten()	# updating growing populations
        dn_dt[-1] = -dn_dt[N:2*N].sum()                             				# updating nutrients

    return dn_dt


def model3_kill(t, n_t, lag, delta, omega, eps):
    n_t[:-1][n_t[:-1] < 0] = 10 ** (-200)               					# avoid negative populations

    Nl = len(np.unique(lag))                						# number of lag times
    Nd = len(np.unique(delta))              						# number of re-entering dormancy rates
    No = len(np.unique(omega))								# number of wake-up rates
    N = Nl * Nd * No                        						# number of species
    
    dn_dt = np.zeros_like(n_t)
    dn_dt[:N] = -n_t[:N] / lag    							# updating dormant populations
    dn_dt[N:2*N] = -dn_dt[:N] - (gamma + delta) * n_t[N:2*N] + n_t[2*N:3*N] / omega     	# updating growing populations
    dn_dt[2*N:3*N] = delta * n_t[N:2*N] - n_t[2*N:3*N] / omega				# updating persistent populations

    return dn_dt
    

# collecting model functions in arrays
ODEs_growth = [model1_grow, model2_grow, model3_grow]
ODEs_decay = [model1_kill, model2_kill, model3_kill]
model_to_states = [2, 2, 3]



#######################
## Solving one cycle ##
#######################
def initialise_system(ab_args, sim_args):
	p, T0, Tab = ab_args
	_, mutation, _, mut_seed = sim_args
	T = T0 + Tab
	n_states = model_to_states[model-1]	# number of states a cell can take
	
	# determining size of system
	Nl = int((T) / dlag)        		# number of lag times
	Nd = int(delta_max / ddelta)          	# number of deltas
	No = int((T) / domega)			# number of omegas
	if model == 1: Nd = 1
	if model <= 2: No = 1
	N = Nl * Nd * No                	# total number of species

	# defining parameter arrays
	lag = dlag * np.arange(0, Nl, 1) + lag_min
	lag = np.tile(np.tile(lag, Nd), No)
	delta = ddelta * np.arange(0, Nd, 1)
	delta = np.tile(np.repeat(delta, Nl), No)
	omega = domega * np.arange(0, No, 1) + lag_min
	omega = np.repeat(np.repeat(omega, Nl), Nd)
	
	# choosing whether to allow mutation and initialising initial conditions in parameter space
	if mutation == False:
	    extinction = False			# no extinction in case without mutation
	    eps = 0				# setting mutation rate to 0 
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
            assert(mut_seed > 1)
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
    model, _, extinction, _ = sim_args
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
    	    	
            c_ext = sum((populations < n_min) * init_alive)	# updating counter
            init_alive = init_alive & (populations > n_min)	# updating mask to contain only species that are alive after ab

        t_low = T						# updating lower bound for integration

    # cycles without antibiotics / solving for period after antibiotics
    sol_cycle = solve_ivp(func_growth, [t_low, T_max], n_t, args=bac_args, events=famine, max_step=dt_max)
    d_t = f * np.array([sol_cycle.y[j][-1] for j in range(N)])
    
    # returning all species to sleep
    for i in range(1, n_states):
    	d_t += f * np.array([sol_cycle.y[i*N + j][-1] for j in range(N)])
    
    # cheking for extinctions  (only for simulations with lower threshold implimented)
    if extinction == True:
        c_ext += sum((d_t < n_min) * init_alive)			# updating counter
        d_t[d_t < n_min] = 0					# truncating extinct species

    return d_t, c_ext


def evolve_system(n_t, bac_args, ab_args, sim_args):
    # breaking if mutation rate = 0
    lag, delta, omega, _ = bac_args
    model = sim_args[0]

    # sorting parameters
    r_arr = sim_args[4]
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


def savedata(ab_args, sim_args, save_data):
    p, T0, Tab = ab_args
    model, _, extinction, _ = sim_args
    T = T0 + Tab
    
    if saving:
	p = np.round(p, 2)
	np.savetxt("../../data/model"+str(model)+"/-average_lag-p" + str(p) + "-T0" + str(T0) + "-T" + str(T), data[0])
	np.savetxt("../../data/model"+str(model)+"/-dominant_lag-p" + str(p) + "-T0" + str(T0) + "-T" + str(T), lag_dom)

	if model >= 2:
	    np.savetxt("../../data/model"+str(model)+"/-average_delta-p" + str(p) + "-T0" + str(T0) + "-T" + str(T), data[1])
	    np.savetxt("../../data/model"+str(model)+"/-dominant_delta-p" + str(p) + "-T0" + str(T0) + "-T" + str(T), delta_dom)

	if model == 3:
	    np.savetxt("../../data/model3/"+"-average_omega-p" + str(p) + "-T0" + str(T0) + "-T" + str(T), data[2])
	    np.savetxt("../../data/model3/-dominant_omega-p" + str(p) + "-T0" + str(T0) + "-T" + str(T), omega_dom)


	if extinction == True:
	    np.savetxt("../../data/model"+str(model)+"/-extinctions-p" + str(p) + "-T0" + str(T0) + "-T" + str(T), data[4])
		    
		  
