# Basic SIR model
# (susceptible, infected, removed)
# https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology#The_SIR_model_without_vital_dynamics
# Described by ODE:
# dS/dt = -beta*I*S
# dI/dt = beta*I*S - gamma*I
# dR/dt = gamma*I
# Solved with Runge-Kutta

import numpy as np
import matplotlib.pyplot as plt

# Initial conditions
I0 = 0.01
S0 = 1 - I0
R0 = 0

# SIR transition intensities
beta = 0.2    # number of contacts * probability of transmission
gamma = 0.05  # recovery rate

# Step-size
h = 0.01
tmax = 100

# One Runge-Kutta step of order 4
def rkstep(f, t, y, h):
    k1 = f(t, y)
    k2 = f(t + h/2, y + h/2*k1)
    k3 = f(t + h/2, y + h/2*k2)
    k4 = f(t + h, y + h*k3)
    return(y + h*(k1/6 + k2/3 + k3/3 + k4/6))

# The f in y'(t) = f(t, y) for the SIR model
# Three components [S, I, R]
def f_SIR(t, y):
    return( np.array([-beta*y[1]*y[0], 
                      beta*y[1]*y[0] - gamma*y[1], 
                      gamma*y[1]]) )

# Apply Runge-Kutta
nsteps = int(tmax/h)
x = range(nsteps)
y = np.array([S0, I0, R0]).reshape(1,3)
for i in x[:-1]:
    y = np.vstack((y, rkstep(f_SIR, h*(i+1), y[i, :], h) ))

# Plot
plt.plot(x, y)
plt.legend(['Susceptible', 'Infected', 'Removed'])
plt.title('Basic SIR model')
plt.ylabel('As fraction of total population')
plt.xlabel('time')
plt.xticks(range(0,nsteps,10*int(1/h)), range(0,tmax,10))
plt.show()
