import numpy as np
import matplotlib.pyplot as plt
import math as mt
 
n = 20000
N = 2    
m1 = 1
m2 = 1
G = 1
pi = mt.pi
dt = 1/20.
mu = 1

a_x = [0] * n 
a_y = [0] * n 
x = [0] * n 
y = [0] * n       
x[0] = 1
y[0] = 0.

Vx = [0] * n 
Vy = [0] * n      # relative velocity
Vx[0] = 0.        # ellipse: root(2)>x>1
Vy[0] = 1.3

norm = [0] * n
norm[0] = (x[0]**2 + y[0]**2)**0.5

for i in range(n-1):
    a_x[i] = (G*m1*m2/(norm[i])**3)*x[i]
    a_y[i] = (G*m1*m2/(norm[i])**3)*y[i]
    Vx[i+1] = Vx[i] - (1/mu)*a_x[i]*dt
    Vy[i+1] = Vy[i] - (1/mu)*a_y[i]*dt
    x[i+1] = x[i] + Vx[i+1]*dt
    y[i+1] = y[i] + Vy[i+1]*dt
    norm[i+1] = (x[i+1]**2 + y[i+1]**2)**0.5

plt.figure()

#Trajectory
plt.plot(x, y, linestyle='dashed', label='Trajectory')

#Earth
plt.scatter(x[0], y[0], color='tab:blue', s=60, label='Earth')

#Sun
plt.scatter(0, 0, color='tab:orange', s=80, label='Sun')


plt.annotate('Earth', (x[0], y[0]), textcoords='offset points', xytext=(8, 6))
plt.annotate('Sun', (0, 0), textcoords='offset points', xytext=(8, 6))

plt.xlabel('X-position')
plt.ylabel('Y-position')
plt.legend()
plt.gca().set_aspect('equal', 'box')
plt.title('The Kepler Problem solved by Leapfrog Integration method')
plt.show()
