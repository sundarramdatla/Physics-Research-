import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np
import math as mt
from tqdm import tqdm 
from mpi4py import MPI

comm = MPI.COMM_WORLD
num_processes = comm.Get_size()
rank = comm.Get_rank()

M = 6
N = 2**M
Nf = N//2 +1
L = 2 * mt.pi
x = np.zeros(N)
y = np.zeros(N)
z = np.zeros(N)
for i in range(0, N):
    x[i] = (L* i)/N
    y[i] = (L * i)/N
    z[i] = (L*i)/N



kx = ky = np.fft.fftfreq(N, d=1./N).astype(int)
kz = kx[:Nf].copy()
kz[-1] *= -1

K = np.array(np.meshgrid(kx, ky, kz, indexing = 'ij'), dtype=int)
K2 = np.sum(K*K, 0)
K_over_K2 = K.astype(float) / np.where(K2 == 0, 1, K2).astype(float)
X,Y,Z = np.meshgrid(x,y,z)
fig = plt.figure()
ax = plt.axes(projection = '3d')


U = np.empty((3,N,N,N), dtype=float)
U_hat = np.empty((3, N, N, Nf), dtype=complex)
P = np.empty((3,N,N,N), dtype=float)
P_hat = np.empty((N, N, Nf), dtype=complex)
curl = np.empty((3,N,N,N), dtype=float)

def fftn_mpi(u,fu):
    if num_processes == 1:
        if u.ndim == 4:
            for c in range(3):
                fu[c] = np.fft.rfftn(u[c], axes = (0,1,2))
        else:
            fu[:] = np.fft.rfftn(u, axes = (0,1,2))
    return fu
def ifftn_mpi(fu,u):
    if num_processes == 1:
        if fu.ndim == 4: 
            for c in range(3):
                u[c] = np.fft.irfftn(fu[c], s=(N,N,N), axes=(0,1,2))
        else:             
            u[:] = np.fft.irfftn(fu, s=(N,N,N), axes=(0,1,2))
    return u

U_hat = fftn_mpi(U, U_hat)
U = ifftn_mpi(U_hat, U)

def Cross(a, b, c):
    c[0] = fftn_mpi(a[1] * b[2] - a[2] * b[1], c[0])
    c[1] = fftn_mpi(a[2] * b[0] - a[0] * b[1], c[1])
    c[2] = fftn_mpi(a[0] * b[1] - a[1] * b[0], c[2])
    return c

def Curl(a, c):
    c[2] = ifftn_mpi(1j * (K[0] * a[1] - K[1] * a[0]), c[2] )
    c[1] = ifftn_mpi(1j * (K[2] * a[0] - K[0] * a[2]), c[1])
    c[0] = ifftn_mpi(1j * (K[1] * a[2] - K[2] * a[1]), c[0]) 

    return c

dU = np.empty((3,N,N,Nf), dtype=complex)
dt = .01
nu = .001

kmax_delias = N/3
dealias = np.array((abs(K[0]) < kmax_delias) * (abs(K[1]) < kmax_delias) * (abs(K[2]) < kmax_delias), dtype=bool)


def computeRHS(dU):
    curl2 = Curl(U_hat, curl)
    dU = Cross(U, curl2, dU)
    dU *= dealias
    P_hat[:] = np.sum(dU * K_over_K2, 0, out=None)
    dU -= P_hat * K
    dU -= nu*K2*U_hat
    return dU
t = 0
T = 1.0
while t <= T:
    t+= dt
    U_hat += computeRHS(dU) *dt
    for i in range(3):
        U[i] = ifftn_mpi(U_hat[i], U[i])
nu = 1./1600
dt = .001
U[0] = np.sin(X[0]*np.cos(X[1]) * np.cos(X[2]))
U[1] = - np.cos(X[0]) * np.sin(X[1]) * np.cos(X[2])
U[2] = 0

for i in range(3):
    U_hat[i] = fftn_mpi(U[i], U_hat[i])




step = max(1, N // 8)          
xs = x[::step]; ys = y[::step]; zs = z[::step]
XS, YS, ZS = np.meshgrid(xs, ys, zs, indexing='ij')

def downsample_U(U):
    return (U[0][::step, ::step, ::step],
            U[1][::step, ::step, ::step],
            U[2][::step, ::step, ::step])

Ux_s, Uy_s, Uz_s = downsample_U(U)

fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_zlabel('z')
ax.set_xlim(0, L); ax.set_ylim(0, L); ax.set_zlim(0, L)
ax.set_title('3D velocity vector field')


q = [ax.quiver(XS, YS, ZS, Ux_s, Uy_s, Uz_s, length=L/(N//step*1.5), normalize=True)]


n_frames = 150          
steps_per_frame = 2    

t_anim = [0.0]        

def advance_solver(nsteps):
    for _ in range(nsteps):
        U_hat[:] = U_hat + dt * computeRHS(dU)
        for i in range(3):
            U[i] = ifftn_mpi(U_hat[i], U[i])

def update(frame):
    advance_solver(steps_per_frame)
    t_anim[0] += dt * steps_per_frame
    q[0].remove()
    Ux_s, Uy_s, Uz_s = downsample_U(U)
    q[0] = ax.quiver(XS, YS, ZS, Ux_s, Uy_s, Uz_s,
                     length=L/(N//step*1.5), normalize=True)
    ax.set_title(f'3D velocity vector field  •  t = {t_anim[0]:.3f}')
    return q

ani = animation.FuncAnimation(fig, update, frames=n_frames, blit=False)
plt.show()

