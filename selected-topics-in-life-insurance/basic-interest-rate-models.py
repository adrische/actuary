# Simple interest rate processes - some examples
#
# - Brownian motion
# - AR(1)-model
# - Vasiček-model
# - Cox-Ingersoll-Ross-model

import numpy as np
import matplotlib.pyplot as plt
from math import sqrt


# Brownian motion / Random walk
def sim_brownian(steps, start=0, sigma=1, discretization=1):
    steps_mod = int(steps / discretization)
    noise = sqrt(discretization) * sigma * np.random.standard_normal(steps_mod)
    return(start + np.cumsum(noise))


# AR(1)-model
def sim_ar1(length_out, start=0, sigma=1, phi=1):
    noise = sigma * np.random.standard_normal(length_out)
    ar1 = np.zeros(length_out)
    for i in range(length_out - 1):
        ar1[i + 1] = phi * ar1[i] + noise[i]
    return(start + ar1)


# Vasiček-model
def sim_vasicek(steps, start=0, sigma=1,
                reversion_level=None, reversion_strength=1,
                discretization=1):
    if reversion_level is None:
        reversion_level = start
    steps_mod = int(steps / discretization)
    v = np.zeros(steps_mod)
    v[0] = start
    for i in range(steps_mod - 1):
        noise = sqrt(discretization) * sigma * np.random.standard_normal()
        dv = reversion_strength * (reversion_level - v[i]) * discretization + noise
        v[i + 1] = v[i] + dv
    return(v)


# Cox-Ingersoll-Ross-model
def sim_cir(steps, start=0, sigma=1,
            reversion_level=None, reversion_strength=1,
            discretization=1):
    if reversion_level is None:
        reversion_level = start
    steps_mod = int(steps / discretization)
    v = np.zeros(steps_mod)
    v[0] = start
    for i in range(steps_mod - 1):
        # It is not ensured that v[i] is positive. This true in the 
        # continuous setting, but not necessarily in this discrete setting
        # (especially when 'discretization' is large compared to 'sigma').
        noise = sqrt(discretization * v[i]) * sigma * np.random.standard_normal()
        dv = reversion_strength * (reversion_level - v[i]) * discretization + noise
        v[i + 1] = v[i] + dv
    return(v)


# Parameters
start = 0.03
sigma = 0.05
phi = 0.99
steps = 20
reversion_level = 0.01
reversion_strength = 1.5
discretization = 0.01

steps_mod = int(steps / discretization)

# Simulate
I_Brownian = sim_brownian(steps, start, sigma, discretization)
I_ar1 = sim_ar1(steps_mod, start, sqrt(discretization) * sigma, phi)
I_V = sim_vasicek(steps,
                  start,
                  sigma,
                  reversion_level,
                  reversion_strength,
                  discretization)
I_cir = sim_cir(steps,
                start,
                sigma,
                reversion_level,
                reversion_strength,
                discretization)

# Plot
plt.plot(np.column_stack((I_Brownian, I_ar1, I_V, I_cir)))
plt.legend(['Brownian motion', 'AR(1)-model',
            'Vasiček-model', 'Cox-Ingersoll-Ross-model'])
plt.title('Basic interest rate models')
plt.xlabel('Time')
plt.ylabel('Interest rate')
plt.xticks(range(0, steps_mod, int(1 / discretization)),
           range(0, steps))
plt.show()
