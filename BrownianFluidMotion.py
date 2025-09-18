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

for i in range (0, n-1):
    a[i] = -aota * v_half[i] + e*E + mt.sqrt((2*aota * k * T )/m) * random.gauss(mu, sigma)
    v_half[i+1] = v_half[i] + a[i] * dt
    x[i+1] = x[i] + v_half[i+1] * dt

print(statistics.mean(v_half))
print(mt.sqrt( (3 * k * T) / m ))

fig, ax1 = plt.subplots()
t = np.arange(n)

line_pos, = ax1.plot(t, x, label="Particle (position)")
ax1.set_xlabel("Time step")
ax1.set_ylabel("X position")

ax2 = ax1.twinx()
v_rms = mt.sqrt((3*k*T)/m)
line_vrms = ax2.axhline(y=v_rms, color='red', label="Average velocity (equipartition check)")

ax2.set_ylabel("Velocity")

lines = [line_pos, line_vrms]
labels = [l.get_label() for l in lines]
plt.title("Brownian motion of a particle inside a viscous fluid")
ax1.legend(lines, labels, loc="upper left")
plt.show()
