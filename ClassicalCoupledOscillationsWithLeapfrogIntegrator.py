import numpy as np
import matplotlib.pyplot as plt
import math as mt

n = 50000
m = 1
dt = .001
k1 = 1
k2 = 1
k3 = 1

E = [0] * n
a_1 = [0] * n
v_1_half = [0] * n
x_1 = [0] * n

a_2 = [0] * n
v_2_half = [0] * n
x_2 = [0] * n

x_1[0] = 5
v_1_half[0] = 0
x_2[0] = 0
v_2_half[0] = 0



def A1(x1, x2):
    ddot_x1 = -1/m * (k1+k2)*x1 + k2 * x2
   
    return ddot_x1
def A2(x1,x2):
     ddot_x2 = -1/m * (k2 + k3) * x2 + k2 * x1
     return ddot_x2
def Energy(v1,v2,x1,x2):
    Energy = .5 * m * (v1**2 + v2**2) + .5 * (k1 * x1**2 + k2*(x1-x2)**2 + k3*x2**2)
    return Energy

for i in range(0, n-1):
    a_1[i] = A1(x_1[i], x_2[i])
    v_1_half[i+1] = v_1_half[i] + a_1[i] * dt
    x_1[i+1] = x_1[i] + v_1_half[i+1] * dt

    a_2[i] = A2(x_1[i], x_2[i])
    v_2_half[i+1] = v_2_half[i] + a_2[i] * dt
    x_2[i+1] = x_2[i] + v_2_half[i+1] * dt
    E[i] = Energy(v_1_half[i+1], v_2_half[i+1], x_1[i+1], x_2[i+1])


plt.plot(x_1, label='x_1', color='blue')  
plt.plot(x_2, label='x_2', color='red')
plt.plot(E, label = "E", color = 'green')   


plt.xlabel('Time Step')
plt.ylabel('Position')
plt.title('Position of Masses over Time')
plt.legend()

# Show the plot
plt.show()