import numpy as np
import matplotlib.pyplot as plt
import math as mt
import statistics
import random

n = 10000
dt = .005
N = 26
mu = 0.0
T = .003
K = 1
m = 1


d = np.zeros((N,2,n))
r = np.zeros((N,2, n))
v = np.zeros((N,2,n))
a = np.zeros((N,2,n))

# 2 particle system
#r[0,:,0] = [ -.39919503, 0.48733313 ]
#r[1,:,0] = [.39919504, -.48733316]

# 3 particle system

r[0,:,0] = [ -.15321632, 0.81846565 ]
r[1,:,0] =  [0.7854202 -.27654359]
r[2,:,0] = [-0.63220389, -0.54192207]

# 5 particle system


r[0, :, 0] = [-1.32271618,  1.0800052]
r[1, :, 0] = [-0.23916632, -1.69786313]
r[2, :, 0] = [-2.2782518,  -1.46649646]
r[3, :, 0] = [-1.45476528, -2.30722003]
r[4, :, 0] = [ 0.82976178,  1.5004903]
r[5, :, 0] = [-2.03378158,  2.01524623]
r[6, :, 0] = [ 0.12792016, -0.62160437]
r[7, :, 0] = [ 2.33519826,  1.54814811]
r[8, :, 0] = [ 0.91952263, -1.44580872]
r[9, :, 0] = [ 2.5789853,  -0.86226929]
r[10, :, 0] = [ 0.93443053, -2.64124794]
r[11, :, 0] = [-2.61579996,  0.87450883]
r[12, :, 0] = [-1.70629177, -0.06761364]
r[13, :, 0] = [ 1.4573863,   2.47048327]
r[14, :, 0] = [ 1.60403819,  0.60243346]
r[15, :, 0] = [-0.93817645,  2.54184149]
r[16, :, 0] = [-0.32159559, -2.85026801]
r[17, :, 0] = [-0.28059057,  1.45555501]
r[18, :, 0] = [ 2.03584454, -1.93428127]
r[19, :, 0] = [ 1.41812142, -0.47393701]
r[20, :, 0] = [-0.60480888,  0.20213499]
r[21, :, 0] = [-1.09937551, -0.99431796]
r[22, :, 0] = [-2.83691016, -0.3868601]
r[23, :, 0] = [ 2.78970561,  0.32083818]
r[24, :, 0] = [ 0.47594899,  0.4197872]
r[25, :, 0] = [ 0.22538721,  2.71830814]


# r[0, :, 0] = [-.945, -1.05]
# r[1, :, 0] = [-1.3833, .29295]
# r[2, :, 0] = [-.4379,1.344]
# r[3, :, 0] = [.4379,-1.344]
# r[4, :, 0] = [.9453,1.0515]
# r[5, :, 0] = [1.3833,-.29295]
# r[6, :, 0] = [0,0]


for i in range (n-1):
   for k in range(N):

      z = mt.sqrt(((2 * K * T * dt))/m) * random.gauss(0,1)
      w = mt.sqrt(((2 * K * T * dt))/m) * random.gauss(0,1)
      for j in range(1,N):
         distance_contribution = (r[k,:,i] - r[(k+j) % N,:,i])/(np.linalg.norm(r[k,:,i] - r[(k+j) % N,:,i]))**3 
         d[k,:,i] = d[k,:,i] + distance_contribution + np.array([z,w])
      a[k,:,i] = -1*r[k,:,i] + d[k,:,i]
      v[k,:,i+1] = v[k,:,i] + a[k,:,i] * dt
      r[k,:,i+1] = r[k,:,i] + v[k,:,i+1] * dt





time_steps = list(range(n)) 
plt.figure()
plt.plot(time_steps, [r[5, 1, i] for i in range(n)], marker='o', linestyle='-', label='Y coordinate')
plt.ylim(-3, 3)
plt.xlabel("Time steps")
plt.ylabel("Y coordinate")
plt.title("Stability of Y coordinate after perturbation")
plt.grid(True)
plt.legend()

plt.figure()
plt.plot(time_steps, [r[5, 0, i] for i in range(n)], marker='o', linestyle='-', label='X coordinate')
plt.ylim(-3, 3)
plt.xlabel("Time steps")
plt.ylabel("X coordinate")
plt.title("Stability of X coordinate after perturbation")
plt.grid(True)
plt.legend()



print(np.var(r[6,:,:]))
print(np.var(r[24,:,:]))
print(np.var(r[20,:,:]))

#1st shell
#r[6,:,:]
#r[24,:,:]
#r[20,:,:]

U_r = 1/3 * (np.var(r[6,:,:] + np.var(r[24,:,:]) + np.var(r[20,:,:])))

print(U_r)


# r[0,:,0] = [ -.15321632, 0.81846565 ]
# r[1,:,0] = [ 0.7854202, -.27654359 ]
# r[2,:,0] = [-0.63220389, -0.54192207]


time_steps = np.arange(n)


x5 = r[5, 0, :]
y5 = r[5, 1, :]
pad = 0.1
xmin, xmax = x5.min(), x5.max()
ymin, ymax = y5.min(), y5.max()
xr = xmax - xmin if xmax > xmin else 1.0
yr = ymax - ymin if ymax > ymin else 1.0



plt.figure()
plt.plot(time_steps, x5, linewidth=1.2, label="X(t)")
plt.xlabel("Time steps")
plt.ylabel("X-coordinate")
plt.title("Perturbation behavior of X coordinate up close")
plt.grid(True)
plt.legend()

plt.figure()
plt.plot(time_steps, y5, linewidth=1.2, label="Y(t)")
plt.xlabel("Time steps")
plt.ylabel("Y-coordinate")
plt.title("Perturbation behavior of Y coordinate up close")
plt.grid(True)
plt.legend()


plt.figure()
plt.plot(x5, y5, linewidth=1.2, label="Particle 5 trajectory")
plt.scatter(x5[0], y5[0], s=40, label="Start", zorder=3)
plt.xlim(xmin - pad*xr, xmax + pad*xr)
plt.ylim(ymin - pad*yr, ymax + pad*yr)
plt.xlabel("X-coordinate")
plt.ylabel("Y-coordinate")
plt.title("Perturbed motion of a electron in lattice")
plt.grid(True)
plt.legend()


plt.figure(figsize=(8, 8))
for i in range(N):
    plt.scatter(r[i, 0, 0], r[i, 1, 0], s=50,)
for i in range(N):
    plt.plot(r[i, 0, :], r[i, 1, :], alpha=0.7)
plt.xlim(-5, 5)
plt.ylim(-5, 5)
plt.xlabel("X-coordinate")
plt.ylabel("Y-coordinate")
plt.title("Wigner Crystal of " + str(N) + " particles in Langevin heat bath")
plt.grid(True)
plt.legend(ncol=2, fontsize=8)

print(np.var(r[6,:,:]))
print(np.var(r[24,:,:]))
print(np.var(r[20,:,:]))
U_r = (np.var(r[6,:,:]) + np.var(r[24,:,:]) + np.var(r[20,:,:])) / 3.0
print(U_r)

plt.show()
