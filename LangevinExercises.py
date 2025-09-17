import numpy as np
import matplotlib.pyplot as plt
import math as mt
import random as random
import statistics

n = 10000

dt = 0.01
m = 1
aota = 1
k = 1
T = .01
mu = 0
e = 1
E = .01
sigma = 1

a = [0] * n
v_half = [0] * n
v_half_KE = [0] * n
x = [0] * n


v_half[0] = 1
v_half_KE[0] = .5 * m * v_half[0]**2
x[0] = 0


for i in range (0, n-1): #remove e*E term below for without external force
        a[i] = -aota * v_half[i] + e*E + mt.sqrt((2*aota * k * T )/m) * random.gauss(mu, sigma)
        v_half[i+1] = v_half[i] + a[i] * dt
        x[i+1] = x[i] + v_half[i+1] * dt



print(statistics.mean(v_half))
print(mt.sqrt( (3 * k * T) / m ))

plt.plot(x)
plt.axhline(y = 1.5 * k * T, color = 'green', label =  "√((3 * k * T) / m)")
plt.plot(v_half, color = 'red', label = "avg. velocity")
plt.legend()


plt.show()