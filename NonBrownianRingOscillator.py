import numpy as np
import matplotlib.pyplot as plt
import math as mt
import random
import statistics

n = 10000
dt = .001
m1 = 1
m2 = 1
m3 = 1
m4 = 1
pi = mt.pi
k = 1
m = 1
K = 1
T = 1



phi_1 =  [0] * n
dphi_1 = [0] * n
ddphi_1 = [0] * n
x_1 = [0] * n
y_1 = [0] * n

phi_2 = [0] * n
dphi_2 = [0] * n
ddphi_2 = [0] * n
x_2 = [0] * n
y_2 = [0] * n

phi_3 = [0] * n
dphi_3 = [0] * n
ddphi_3 = [0] * n
x_3 = [0] * n
y_3 = [0] * n

phi_4 = [0] * n
dphi_4 = [0] * n
ddphi_4 = [0] * n
x_4 = [0] * n
y_4 = [0] * n

# normal modes occur when all masses are displaced equally, not displaced at all, or m1, m3 are displaced opposite and m2 and m4 are displaced opposite

phi_1[0] = -mt.pi/2
x_1[0] = mt.cos(phi_1[0])
y_1[0] = mt.sin(phi_1[0])

phi_2[0] = mt.pi /3
x_2[0] = mt.cos(phi_2[0])
y_2[0] = mt.sin(phi_2[0])

phi_3[0] = -mt.pi
x_3[0] = mt.cos(phi_3[0])
y_3[0] = mt.sin(phi_3[0])

phi_4[0] = mt.pi/2
x_4[0] = mt.cos(phi_4[0])
y_4[0] = mt.sin(phi_4[0])

dphi_1[0] = 0
dphi_2[0] = 0
dphi_3[0] = 0
dphi_4[0] = 0

for i in range(0,n-1):
    ddphi_1[i] = -(k/m)* (2* phi_1[i] - phi_2[i] - phi_4[i]) + + mt.sqrt(((2 * K * T * dt))/m) * random.gauss(0,1)
    dphi_1[i+1] = dphi_1[i] + ddphi_1[i] * dt
    phi_1[i+1] = (phi_1[i] + dphi_1[i+1] * dt) % (2*pi)
    x_1[i+1] = mt.cos(phi_1[i+1])
    y_1[i+1] = mt.sin(phi_1[i+1])

    ddphi_2[i] = -(k/m) * (2* phi_2[i] - phi_1[i] - phi_3[i]) + + mt.sqrt(((2 * K * T * dt))/m) * random.gauss(0,1)
    dphi_2[i+1] = dphi_2[i] + ddphi_2[i] * dt
    phi_2[i+1] = (phi_2[i] + dphi_2[i+1] * dt) % (2*pi)
    x_2[i+1] = mt.cos(phi_2[i+1])
    y_2[i+1] = mt.sin(phi_2[i+1])

    ddphi_3[i] = -(k/m) * (2* phi_3[i] - phi_2[i] - phi_4[i]) + + mt.sqrt(((2 * K * T * dt))/m) * random.gauss(0,1)
    dphi_3[i+1] = dphi_3[i] + ddphi_3[i] * dt
    phi_3[i+1] = (phi_3[i] + dphi_3[i+1] * dt) % (2*pi)
    x_3[i+1] = mt.cos(phi_3[i+1])
    y_3[i+1] = mt.sin(phi_3[i+1])

    ddphi_4[i] = -(k/m) * (2* phi_4[i] - phi_1[i] - phi_3[i]) + + mt.sqrt(((2 * K * T * dt))/m) * random.gauss(0,1)
    dphi_4[i+1] = dphi_4[i] + ddphi_4[i] * dt
    phi_4[i+1] = (phi_4[i] + dphi_4[i+1] * dt) % (2*pi)
    x_4[i+1] = mt.cos(phi_4[i+1])
    y_4[i+1] = mt.sin(phi_4[i+1])

t = np.arange(n)


p1 = m * np.array(dphi_1)
plt.figure()
plt.plot(phi_1, p1, label="Phase space trajectory")
plt.xlabel("Position")
plt.ylabel("Momentum")
plt.title("Phase space of mass 1")
plt.legend()

plt.figure()
plt.plot(x_1, y_1, color='green', label="Trajectory 1")
plt.plot(x_2, y_2, color='red', label="Trajectory 2")
plt.plot(x_3, y_3, color='blue', label="Trajectory 3")
plt.plot(x_4, y_4, color='purple', label="Trajectory 4")
plt.scatter([x_1[0], x_2[0], x_3[0], x_4[0]],
            [y_1[0], y_2[0], y_3[0], y_4[0]],
            c=['green','red','blue','purple'], s=30)
plt.xlabel("X position")
plt.ylabel("Y position")
plt.title("Cartesian graph of mass position")
plt.gca().set_aspect('equal', 'box')
plt.legend()

plt.figure()
plt.plot(t, phi_1, label="Mass 1", color='green')
plt.plot(t, phi_2, label="Mass 2", color='red')
plt.plot(t, phi_3, label="Mass 3", color='blue')
plt.plot(t, phi_4, label="Mass 4", color='purple')
plt.xlabel("Time steps")
plt.ylabel("Angle (radians)")
plt.title("Polar angle of mass vs. time")
plt.legend()





plt.show()
