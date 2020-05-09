# coding: utf-8
import numpy as np
import matplotlib.pyplot as plt

def XDiff(x, y, a):
  return a / (1.0+ y **2) - x #Toggle switchモデル

def YDiff(x, y, a):
  return a / (1.0+ x **2) - y

def zeroDiff(range_array, a): #微分方程式=0
  return a / (1.0+ range_array **2)

x  = np.arange(0, 5.0, 0.1) #タンパク質濃度の定義域
a = [2.0, 4.0] #抑制速度パラメータ
nullcline = np.empty(len(a),dtype="O") #共通ヌルクラインを作成する

for ii in range(len(a)):
  nullcline[ii] = zeroDiff(x, a[ii])

for i in range(2):
  plt.figure()
  plt.plot(nullcline[i], x) #xヌルクライン
  plt.plot(x, nullcline[i]) #yヌルクライン
  #ベクトル場
  for j in np.arange(0, 4.1, 0.15):
    for k in np.arange(0, 4.1, 0.15):
      x_direction = XDiff(j, k, a[i])
      y_direction = YDiff(j, k, a[i])
      unit_scale = np.sqrt(x_direction**2 + y_direction**2)
      if unit_scale > 0:
        plt.quiver(j, k, x_direction/unit_scale, y_direction/unit_scale,
                   scale=40, width=0.005, headlength=3, headaxislength=3)
  
  plt.xlim(0, 4)
  plt.ylim(0, 4)

"""
def checkStablePoint(a, x, y): #det(J)が正かどうか確かめる
  if 4* x * y * (a **2) <= ((1.0+ x **2) **2) * ((1.0+ y **2) **2):
    return "ko" #安定
  else:
    return "wo" #不安定
#1つ目の固定点解析
figure()
xlim(0, 4)
ylim(0, 4)
plot(Nullcline[0], x)
plot(x, Nullcline[0])
plot(1.0, nullcline(1.0, a[0]), checkStablePoint(a[0], 1.0, nullcline(1.0, a[0])), markersize =10)
#2つ目の固定点解析
figure()
xlim(0, 4)
ylim(0, 4)
plot(Nullcline[1], x)
plot(x, Nullcline[1])

import sympy
xx = Symbol('xx') #交点を求めるため三次方程式を解く (sympy.solve)
vertex = sympy.solve(xx **5- a[1] * xx **4+2* xx **3-2* a[1] * xx **2  + (a[1] **2 + 1.0) * xx - a[1], xx)
print vertex
plot(vertex[0].evalf(), nullcline(vertex[0].evalf(), a[1]), checkStablePoint(a[1], vertex[0].evalf(), nullcline(vertex[0].evalf(), a[1])), markersize =10)
plot(vertex[2].evalf(), nullcline(vertex[2].evalf(), a[1]), checkStablePoint(a[1], vertex[2].evalf(), nullcline(vertex[2].evalf(), a[1])), markersize =10)
plot(vertex[3].evalf(), nullcline(vertex[3].evalf(), a[1]), checkStablePoint(a[1], vertex[3].evalf(), nullcline(vertex[3].evalf(), a[1])), markersize =10)
"""

plt.show() #全figure描画
