# Example 2.4.2 (Disability pension)
# A simple three-state model: alive, disabled, dead

# Initial age
age1 = 30
# Final age
age2 = 111
# Step width
h = 0.001

# One step of Runge-Kutta scheme of order 4
def rkstep(f, t, y, h, add_data_list=None):
    k1 = f(t, y, add_data_list)
    k2 = f(t + h/2, y + h/2*k1, add_data_list)
    k3 = f(t + h/2, y + h/2*k2, add_data_list)
    k4 = f(t + h, y + h*k3, add_data_list)
    return(y + h*(k1/6 + k2/3 + k3/3 + k4/6))

# Transition intensities are given
## alive -> disabled
def sigma(x):
    return(0.0004 + 10**(0.06*x - 5.46))

## alive -> dead
def mu(x):
    return(0.0005 + 10**(0.038*x - 4.12))

# ODE's for conditional probabilities:
# dP(at age 2 in state 2, given at age 1 in state 1)/dt = 

## alive -> alive: 
def f_alive_alive(t, y, add_data_list=None):
    return(-y * (mu(t) + sigma(t)))

## alive -> disabled
def f_alive_disabled(t, y, add_data_list=None):
    return(-y*mu(t) + add_data_list*sigma(t))

## alive -> Dead
def f_alive_Dead(t, y, add_data_list=None):
    return( (add_data_list[0] + add_data_list[1]) * mu(t) )

## disabled -> alive (= 0)
def f_disabled_alive(t, y, add_data_list=None):
    return( 0 )

## disabled -> disabled
def f_disabled_disabled(t, y, add_data_list=None):
    return( -y * mu(t))

## disabled -> Dead
def f_disabled_Dead(t, y, add_data_list=None):
    return( y * mu(t))

## Dead -> Dead (= 1)
def f_Dead_Dead(t, y, add_data_list=None):
    return( 1 )

# Boundary conditions (initial probabilities)
p_alive_alive       = [1]
p_alive_disabled    = [0]
p_alive_Dead        = [0]
p_disabled_alive    = [0]
p_disabled_disabled = [1]
p_disabled_Dead     = [0]
p_Dead_Dead         = [1]

# Apply numerical scheme
nsteps = int((age2 - age1) / h)
x = range(nsteps)
for i in x[:-1]:
    t = age1 + h*(i+1)
    p_alive_alive       .append( rkstep(f_alive_alive, t, p_alive_alive[i], h) )
    p_alive_disabled    .append( rkstep(f_alive_disabled, t, p_alive_disabled[i], h, p_alive_alive[i]) )
    p_alive_Dead        .append( rkstep(f_alive_Dead, t, p_alive_Dead[i], h, [p_alive_disabled[i], p_alive_alive[i]]) )
    p_disabled_alive    .append( 0 )
    p_disabled_disabled .append( rkstep(f_disabled_disabled, t, p_disabled_disabled[i], h) )
    p_disabled_Dead     .append( rkstep(f_disabled_Dead, t, p_disabled_Dead[i], h) )
    p_Dead_Dead         .append(1)

# Plot
import matplotlib.pyplot as plt
plt.plot(x, p_alive_alive, 
         x, p_alive_disabled, 
         x, p_alive_Dead, 
         # x, p_disabled_alive, 
         x, p_disabled_disabled,
         x, p_disabled_Dead,
         x, p_Dead_Dead)
plt.legend(['Still alive', 'Disabled', 'Died', 'Stay disabled', 'Died disabled', 'Stay dead'])
plt.xlabel('Age')
plt.ylabel('Probability')
plt.xticks([int(ind/h) for ind in range(0, age2-age1, 5)], range(age1, age2, 5))
plt.show()

