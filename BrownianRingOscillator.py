import numpy as np
import matplotlib.pyplot as plt
import math as mt
import random
import statistics


n = 100000
dt = .001 # set to .001 or .01
m1 = 1
m2 = 1
m3 = 1
m4 = 1
pi = mt.pi

k = 1
m = 1
K = 1
T = 1 # Set to 0 for code to remove noise effects
mu = 0
sigma = 1

# normal modes occur when all masses are displaced equally, not displaced at all, or m1, m3 are displaced opposite and m2 and m4 are displaced opposite


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

phi_1[0] = -mt.pi/6
x_1[0] = mt.cos(phi_1[0])
y_1[0] = mt.sin(phi_1[0])

phi_2[0] = -mt.pi/4
x_2[0] = mt.cos(phi_2[0])
y_2[0] = mt.sin(phi_2[0])

phi_3[0] = mt.pi/6

x_3[0] = mt.cos(phi_3[0])
y_3[0] = mt.sin(phi_3[0])

phi_4[0] = mt.pi/4

x_4[0] = mt.cos(phi_4[0])
y_4[0] = mt.sin(phi_4[0])



dphi_1[0] = 0
dphi_2[0] = 0
dphi_3[0] = 0
dphi_4[0] = 0

for i in range(0,n-1):
    
    ddphi_1[i] = -(k/m)* (2* phi_1[i] - phi_2[i] - phi_4[i]) + mt.sqrt((2 * k * T )/m) * random.gauss(mu, sigma)
    dphi_1[i+1] = dphi_1[i-1] + ddphi_1[i] * dt
    phi_1[i+1] = (phi_1[i] + dphi_1[i+1] * dt) % (2*pi)
    x_1[i+1] = mt.cos(phi_1[i+1])
    y_1[i+1] = mt.sin(phi_1[i+1])


    ddphi_2[i] = -(k/m) * (2* phi_2[i] - phi_1[i] - phi_3[i]) + mt.sqrt((2 * k * T )/m) * random.gauss(mu, sigma)
    dphi_2[i+1] = dphi_2[i-1] + ddphi_2[i] * dt
    phi_2[i+1] = (phi_2[i] + dphi_2[i+1] * dt) % (2*pi)
    x_2[i+1] = mt.cos(phi_2[i+1])
    y_2[i+1] = mt.sin(phi_2[i+1])

    ddphi_3[i] = -(k/m) * (2* phi_3[i] - phi_2[i] - phi_4[i]) + mt.sqrt((2 * k * T )/m) * random.gauss(mu, sigma)
    dphi_3[i+1] = dphi_3[i-1] + ddphi_3[i] * dt
    phi_3[i+1] = (phi_3[i] + dphi_3[i+1] * dt) % (2*pi)
    x_3[i+1] = mt.cos(phi_3[i+1])
    y_3[i+1] = mt.sin(phi_3[i+1])

    ddphi_4[i] = -(k/m) * (2* phi_4[i] - phi_1[i] - phi_3[i]) + mt.sqrt((2 * k * T )/m) * random.gauss(mu, sigma)
    dphi_4[i+1] = dphi_4[i-1] + ddphi_4[i] * dt
    phi_4[i+1] = (phi_4[i] + dphi_4[i+1] * dt) % (2*pi)
    x_4[i+1] = mt.cos(phi_4[i+1])
    y_4[i+1] = mt.sin(phi_4[i+1])


print((K * T)/m)


 
t = np.arange(n) * dt

plt.figure()
plt.plot(t, phi_1, label="Mass 1")
plt.plot(t, phi_2, label="Mass 2")
plt.plot(t, phi_3, label="Mass 3")
plt.plot(t, phi_4, label="Mass 4")

plt.xlabel("Time Steps")
plt.ylabel("Angle [rad]")  # amplitude here is the angle in radians
plt.title("Coupled Brownian Motion in Thermal Bath: Angle vs. Time Steps")
plt.legend()
plt.tight_layout()
plt.show()



# plt.plot(x_1, y_1, color='green' )
# plt.scatter(x_1[0], y_1[0], color = 'green')
# plt.plot(x_2, y_2, color='red')
# plt.scatter(x_2[0], y_2[0], color = 'red')
# plt.plot(x_3, y_3, color='blue')
# plt.scatter(x_3[0], y_3[0], color = 'blue')
# plt.plot(x_4, y_4, color='purple')
# plt.scatter(x_4[0], y_4[0], color = 'purple')
# plt.show()




 
