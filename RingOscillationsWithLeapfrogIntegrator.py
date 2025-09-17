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
T = 0


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

phi_1[0] = -mt.pi/2
x_1[0] = mt.cos(phi_1[0])
y_1[0] = mt.sin(phi_1[0])

phi_2[0] = mt.pi /3
x_2[0] = mt.cos(phi_2[0])
y_2[0] = mt.sin(phi_2[0])

phi_3[0] = -mt.pi/2

x_3[0] = mt.cos(phi_3[0])
y_3[0] = mt.sin(phi_3[0])

phi_4[0] = mt.pi/3

x_4[0] = mt.cos(phi_4[0])
y_4[0] = mt.sin(phi_4[0])



dphi_1[0] = 0
dphi_2[0] = 0
dphi_3[0] = 0
dphi_4[0] = 0

# + mt.sqrt(((2 * K * T * dt))/m) * random.gauss(0,1)

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



plt.plot(phi_1, label="Angle of mass 1", color = 'green')
plt.plot(phi_2, label = "Angle of mass 2", color='red')
plt.plot(phi_3, label= "Angle of mass 3", color= "Blue")
plt.plot(phi_4, label="Angle of mass 4", color="Purple")
plt.legend()

dphi_1 = np.array(dphi_1) 
dphi_2 = np.array(dphi_2)
dphi_3 = np.array(dphi_3)
dphi_4 = np.array(dphi_4)


avg_dphi_squared = np.mean(0.5 * m * dphi_1 ** 2 + .5 * m * dphi_2**2 + .5 * m * dphi_3**2 + .5 * m * dphi_4**2 )

# avg_dphi_1_squared = np.mean(.5 * m * dphi_1**2)
# print(avg_dphi_1_squared)

print(avg_dphi_squared)

print(2*K*T)





plt.show()


# plt.plot(x_1, y_1, color='green' )
# plt.scatter(x_1[0], y_1[0])
# plt.plot(x_2, y_2, color='red')
# plt.scatter(x_2[0], y_2[0])
# plt.plot(x_3, y_3, color='blue')
# plt.scatter(x_3[0], y_3[0])
# plt.plot(x_4, y_4, color='purple')
# plt.scatter(x_4[0], y_4[0])
# plt.show()

#phase space for mass 1

# plt.plot(phi_1,m*dphi_1)

plt.show()

 
