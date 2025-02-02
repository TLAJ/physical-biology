# coding: utf-8
import numpy as np
import matplotlib.pyplot as plt

# from sympy import *

N_NS = 10000000.0
T = 298.0
k_B = 1.38e-23
beta = 1.0 / (k_B * T)
delta_ep_ad = -13.12 * k_B * T
delta_ep_pd = -5.3 * k_B * T
ep_ap = 0.0


def Freg(A):
    numerator = 1.0 + (A / N_NS) * np.exp(-beta * delta_ep_ad) * np.exp(-beta * ep_ap)
    denominator = 1.0 + (A / N_NS) * np.exp(-beta * delta_ep_ad)
    return numerator / denominator


def pbound(P, A):
    return 1.0 / (1.0 + (N_NS / (P * Freg(A))) * np.exp(beta * delta_ep_pd))


A_range = np.linspace(0, 100, 100)  # Activator個数の変化範囲
P = 500.0

fig = plt.figure()

ep_ap = -3.0 * k_B * T
plt.plot(A_range, pbound(P, A_range), "r", lw=2)
ep_ap = -4.0 * k_B * T
plt.plot(A_range, pbound(P, A_range), "b", lw=2)
ep_ap = -5.0 * k_B * T
plt.plot(A_range, pbound(P, A_range), "g", lw=2)
plt.xlabel("number of activator molecules")
plt.ylabel(r"$p_{bound}$")
plt.legend(
    [
        r"$\varepsilon_{ap} = -3k_{B}T$",
        r"$\varepsilon_{ap} = -4k_{B}T$",
        r"$\varepsilon_{ap} = -5k_{B}T$",
    ]
)

plt.show()  # 全figure描画
