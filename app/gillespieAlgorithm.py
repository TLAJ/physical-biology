# coding: utf-8

"""
Created on 2015/09/14

@author: Kaoru

Gillespie Algorithm:
Solve the simple case of a promoter expressing mRNA
at a constant rate where the mRNA can in turn decay.

Defining some parameters:
new mRNA is produced approximately every 3.3 seconds
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import poisson

r = 20  # probability of mRNA production per unit time: [min^-1]

# life time of 1.5 minutes measured for the mRNA
gamma = 1 / 1.5  # probability of mRNA decay per mRNA and unit time: [min^-1]

# Iterations of the algorithm: for each iteration we will flip a coin.
MaxT = 500

# Simulations start with no mRNA
m = np.empty(MaxT + 1)
tau = np.empty(MaxT)  # Time to the next reaction
m[0] = 0  # initial condition = 0

for i in range(MaxT):
    k1 = r
    k2 = m[i] * gamma
    k0 = k1 + k2  # Sum of all probabilities per unit time

    # Determine the vector tau
    tau[i] = (
        1 / k0 * np.log(1 / np.random.rand())
    )  # Uniform probability between 0 and 1

    # Flip a coin to determine which one of the reactions will take place
    CoinFlip = np.random.rand()

    # mRNA production / mRNA decay
    if CoinFlip <= k1 / k0:
        m[i + 1] = m[i] + 1
    else:
        m[i + 1] = m[i] - 1


# Calculate the time axis
T = np.empty(MaxT + 1)
T[0] = 0

for i in range(MaxT):
    T[i + 1] = np.sum(tau[0:i])

# Plot results with the deterministic solution
plt.figure()
plt.plot(T, m)
plt.plot(T, r / gamma * (1 - np.exp(-T * gamma)), "-k")
plt.xlabel("time (min)")
plt.ylabel("number of mRNA molecules")

plt.show(block=False)

# Distribution in steady state
r = 20
gamma = 1 / 1.5
MaxT = 10000

# Simulations start with no mRNA
m = np.empty(MaxT + 1)
tau = np.empty(MaxT)  # Time to the next reaction
m[0] = np.round(r / gamma)  # Notice the different initial condition

for i in range(MaxT):
    k1 = r
    k2 = m[i] * gamma
    k0 = k1 + k2  # Sum of all probabilities per unit time

    # Determine the vector tau
    tau[i] = 1 / k0  # Average

    # Flip a coin to determine which one of the reactions will take place
    CoinFlip = np.random.rand()

    # mRNA production / mRNA decay
    if CoinFlip <= k1 / k0:
        m[i + 1] = m[i] + 1
    else:
        m[i + 1] = m[i] - 1

# Check simulation to see if it follows a Poisson distribution
# Variance should be equal to the mean
# Weighting each mRNA occurence using the time over which it existed
MeanM = np.sum(m[0:-1] * tau) / np.sum(tau)
SecondMoment = np.sum(np.power(m[0:-1], 2) * tau) / np.sum(tau)

VarianceM = SecondMoment - MeanM ** 2

# Compute the probability distribution
# Finding a certain number of mRNA molecules is related to the period of time
MaxmRNA = int(np.max(m))

p = np.zeros(MaxmRNA)  # Observed number of mRNA molecules with length MaxmRNA

for i in range(MaxT):  # Time that number existed until next reaction
    p[m[i] - 1] = p[m[i] - 1] + tau[i]
# Normalize the distribution
p = p / np.sum(p)

# Plot histogram of p with expected Poisson distribution
plt.figure()
plt.bar(range(1, MaxmRNA + 1), p)

X = range(MaxmRNA)
Y = poisson.pmf(X, r / gamma)
plt.plot(X, Y, "-r")
plt.xlabel("number of mRNA molecules")
plt.ylabel("probability")
plt.legend(["Poisson distribution", "Simulation"])

plt.show()

"""
The agreement might not be perfect. At least, it is probably not as good
as the agreement in the figure associated with this computational
exploration in the chapter. The main reason for the disagreement is the
fact that we did not perform a coin flip to determine the values of tau
(unlike what we actually did to generate the figure). Instead we only took
its average for each iteration. A second cause of
disagreement is the number of iterations employed in our simulation. Both
of these limitations are easily solvable and are, in fact, taken into
account when this algorithm is used in the context of research.
"""
