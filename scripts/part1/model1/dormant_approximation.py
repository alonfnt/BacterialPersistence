# Scipt for investigating whether neglecting dormant population at time Ts: S(Ts) = 0 is a good approximation.
# Plotting the ration of dormant species and growing/awake species at the time Ts in the cycle without antibiotics (because this corresponds to the shortest Ts)

import numpy as np
import matplotlib.pyplot as plt

lag_max = 24			# upper limit on lag time
lag = np.linspace(1, lag_max, 1000)	# array of meaningful lag times

# d(t) / g(t) in cycles without antibiotics
def dt_gt(lag, Ts):
    return lag * np.exp(-Ts / lag) /(np.exp(Ts) + lag * np.exp(-Ts / lag))


fig, ax = plt.subplots(1, 1, figsize=(5, 3))
ax.set(xlabel="Lag time [h]", ylabel=r"đ(t) / g(t)")
ax.set(yscale="log", ylim=(10**(-10), 0.1))

ax.plot(lag, dt_gt(lag, 13), label=r"t = 13")
ax.plot(lag, dt_gt(lag, 15), label=r"t = 15")
ax.plot(lag, dt_gt(lag, lag), "--", color="grey", label=r"t = $\lambda$")
# ax.hlines(10**(-6), lag.min(), lag.max(), ls='dashed', color='black')
# ax.hlines(10**(-5), lag.min(), lag.max(), ls='dashed', color='black')
ax.legend()

fig.tight_layout()
fig.show()
fig.savefig("../../../figures/model1/dormant_approx.png", dpi=100)
