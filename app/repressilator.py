# coding: utf-8
"""
Created on 2015/09/16

@author: Kaoru
"""
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import odeint


def ODE(conc, time, param):
    x_mRNA = conc[0]
    y_mRNA = conc[1]
    z_mRNA = conc[2]
    x = conc[3]
    y = conc[4]
    z = conc[5]

    alpha = param[0]
    alpha0 = param[1]
    beta = param[2]
    n = param[3]

    dx_mRNA_dt = -x_mRNA + alpha / (1 + y ** n) + alpha0
    dy_mRNA_dt = -y_mRNA + alpha / (1 + z ** n) + alpha0
    dz_mRNA_dt = -z_mRNA + alpha / (1 + x ** n) + alpha0
    dx_dt = beta * (x_mRNA - x)
    dy_dt = beta * (y_mRNA - y)
    dz_dt = beta * (z_mRNA - z)

    return [dx_mRNA_dt, dy_mRNA_dt, dz_mRNA_dt, dx_dt, dy_dt, dz_dt]


conc0 = [4.0, 5.0, 8.0, 4.0, 5.0, 6.0]  # 各mRNA, タンパク質の初期濃度
param = [100.0, 1.0, 5.0, 2.0]  # パラメータ：alpha, alpha0, beta, n
sim_time = np.arange(0.01, 50.0, 0.05)  # シミュレーション時間
timecourse = odeint(ODE, conc0, sim_time, (param,))

# 時系列プロット
plt.figure()
plt.plot(sim_time, timecourse[:, 0:3])
plt.xlabel("time")
plt.ylabel("concentration")
plt.legend(["TetR", "LacI", "lambda cI"])

# 3つのタンパク質の濃度変化を3次元プロット
ax = Axes3D(plt.figure())
ax.plot(timecourse[:, 3], timecourse[:, 4], timecourse[:, 5])
plt.gca().set_xlabel("TetR")
plt.gca().set_ylabel("LacI")
plt.gca().set_zlabel("lambda cI")

plt.show()
