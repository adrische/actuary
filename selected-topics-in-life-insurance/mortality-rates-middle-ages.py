# Exercise 4.3.1
# Mortality rates in the middle ages

from scipy import optimize
from math import log
from math import exp

# Ages at time of death (given in lecture notes)
ages = [47, 67, 38, 45, 48, 48, 51, 76, 32, 43, 42, 26, 25, 40, 53, 64, 49, 42,
        57, 30, 1, 1, 6, 9, 12, 14, 18]


# Model ansatz
# mu(x) = Probability to die between age x and x+1,
#         conditional on being alive at age x.
def mu(x, a, b, c):
    return(exp(a + b * x + c * x ** 2))


# Negative log-likelihood
# Probability of dying at age x is the product of the probabilities of not dying
# during years before x, and the probability of dying between age x and x+1.
def negll(abc, ages):
    a, b, c = abc
    ll = 0
    for age in ages:
        for x in range(age):
            # not dying before age of death
            ll += log(1 - mu(x, a, b, c))
        # dying at age of death
        ll += a + b * age + c * age ** 2  # same as log(mu(x, a, b, c))
    return(-ll)


# Solving for three parameters in model ansatz
# Initial values are the solution given in lecture notes
print(optimize.minimize(negll, (-4.36, 1.01e-2, 4.08e-4), args=(ages,)))
# This gives [a, b, c] = [-4.48664309e+00, 1.52875034e-02, 3.38010681e-04],
# different from the given solution.
