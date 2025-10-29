import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.animation import FuncAnimation
import numpy as np
from scipy.stats import pearsonr
from tqdm import tqdm 
from scipy.fft import rfft, rfftfreq
from scipy.signal import detrend
from scipy.signal import find_peaks


N = 161 # 41 points per axis
domain_size = 5
iterations = 10000
reynolds = 100
radius = .2
kinematic_viscosity = 2*radius/reynolds
truncation = 6000
Fs = 20
height_max = .025
height_min = .02
height_threshold = (height_min, height_max)
density = 1.0
pressure_poisson_iterations = 50
U_in = 1
CFL = .5




def main():
    dt_list = []
    element_length = domain_size / (N - 1) #initialize the grid
    x = np.linspace(0, domain_size, N)
    y = np.linspace(0, domain_size, N)
    X,Y = np.meshgrid(x,y)
    x_c = domain_size / 5
    y_c = domain_size / 2
   
    circle = (X - x_c)**2 + (Y-y_c)**2 <= radius**2

    u_prev = np.zeros_like(X)
    u_prev[:, 0] = U_in           # initialize initial conditions of fluid
    v_prev = np.zeros_like(X)
    p_prev = np.zeros_like(X)
    c = np.zeros(iterations)
    def central_difference_x(f):
        diff = np.zeros_like(f)
        diff[1: -1, 1:-1] = ( f[1: -1, 2: ] - f[1:-1, 0: -2] ) / (2 * element_length)
        return diff
    def central_difference_y(f):
        diff = np.zeros_like(f)
        diff[1: -1, 1: -1] = (f[2:, 1:-1] - f[0: -2, 1: -1 ]) / (2 * element_length)
        return diff
    
    def laplace(f):
        diff = np.zeros_like(f)
        diff[1: -1, 1:-1] = (f[1: -1, 0: -2] + f[0: -2, 1: -1] - 4 * f[1:-1, 1:-1] + f[1: -1, 2: ] + f[2: , 1:-1]) / (element_length **2)
        return diff
    
    def compute_dt(u,v):
        max_u = np.max(u)
        max_v = np.max(v)
        dt = CFL * element_length / max(max_u, max_v, 1e-12)
        dt_list.append(dt)
        return dt



    
    for step in tqdm(range(iterations)): 
        
        du_prev_dx = central_difference_x(u_prev)
        du_prev_dy = central_difference_y(u_prev)
        dv_prev_dx = central_difference_x(v_prev)
        dv_prev_dy = central_difference_y(v_prev)
        laplace_u_prev = laplace(u_prev)
        laplace_v_prev = laplace(v_prev)
        time_step_length = compute_dt(u_prev,v_prev)

        #Find the average velocity magnitude of previous field
        magnitudes = np.sqrt(u_prev[:]**2 + v_prev[:]**2)
        V_prev = np.average(magnitudes)
        w_prev = (u_prev + v_prev).ravel()

        #Solve momentum equation w/o pressure consideration using Chorin's projection method
        u = (u_prev + time_step_length * ( - (u_prev * du_prev_dx + v_prev * du_prev_dy) + kinematic_viscosity * laplace_u_prev ))
        v = (v_prev + time_step_length * (- (u_prev * dv_prev_dx + v_prev * dv_prev_dy) + kinematic_viscosity * laplace_v_prev))
        
        #Velocity BC: homogenous dirichlect BC everywhere 
        # except for horizontal velocity at top, applied here:
        u[circle] = 0
        v[circle] = 0
        u[0, :] = u[1,:]
        u[:, 0] = U_in
        u[:, -1] = u[:, -2]
        u[-1, :] = u[-2,:]
        v[0, :] = 0.0
        v[:, 0] = 0
        v[:, -1] = v[:, -2]
        v[-1, :] = 0
    
        du_dx = central_difference_x(u)
        dv_dy = central_difference_y(v)

        #Compute pressure correction by solving pressure-poisson equation

        rhs = ( density / time_step_length * (du_dx + dv_dy))

        for it in range(pressure_poisson_iterations): #Uses Jacobi method to solve this pressure poisson equation
            p_next = np.zeros_like(p_prev)
            p_next[1:-1, 1:-1] = 1/4 * (+ p_prev[1:-1, 0:-2] + p_prev[0: -2, 1:-1] + p_prev[1:-1, 2:] + p_prev[2:, 1:-1] - element_length**2 * rhs[1: -1, 1:-1])

            #pressure boundary conditions: homogenous neumann BC everywhere
            #except for top, homogenous dirichlect BC given here:
            
            
            p_next[-1, :] = p_next[-2, :]
            p_next[0,  :] = p_next[1,  :]
            p_next[:, 0]  = p_next[:, 1]
            p_next[:, -1] = 0
            
            
            mask_left  = circle & (~np.roll(circle, 1, axis=1))
            ii, jj = np.where(mask_left & (np.arange(N)[None, :] > 0))
            p_next[ii, jj] = p_next[ii, jj - 1]

            mask_right = circle & (~np.roll(circle, -1, axis=1))
            ii, jj = np.where(mask_right & (np.arange(N)[None, :] < N - 1))
            p_next[ii, jj] = p_next[ii, jj + 1]

            mask_down  = circle & (~np.roll(circle, 1, axis=0))
            ii, jj = np.where(mask_down & (np.arange(N)[:, None] > 0))
            p_next[ii, jj] = p_next[ii - 1, jj]

            mask_up    = circle & (~np.roll(circle, -1, axis=0))
            ii, jj = np.where(mask_up & (np.arange(N)[:, None] < N - 1))
            p_next[ii, jj] = p_next[ii + 1, jj]

            
            

            p_prev = p_next

        #Now we apply the pressure correction to velocity to keep fluid incompressible
        dp_next_dx = central_difference_x(p_next)
        dp_next_dy = central_difference_y(p_next)

        u_next = (u - time_step_length / density * dp_next_dx)
        v_next = (v - time_step_length / density * dp_next_dy)

        #Reapply boundary conditions
        u_next[circle] = 0
        v_next[circle] = 0
        u_next[0, :] = 0
        u_next[:, 0] = U_in
        u_next[:, -1] = u_next[:, -2]
        u_next[-1, :] = 0
        v_next[0, :] = 0.0
        v_next[:, 0] = 0
        v_next[:, -1] = v_next[:, -2]
        v_next[-1, :] = 0

        w_next = (u_next + v_next).ravel()

        magnitudes_next = np.sqrt(u_next[:]**2 + v_next[:]**2)
        V_next = np.average(magnitudes_next)
        A = np.minimum(V_prev, V_next) / np.maximum(V_prev, V_next)
        #Creating a signal out of correlation coefficients of each vortice
        correlation_coefficient,p_value = pearsonr(w_prev, w_next)
        c[step] = A * correlation_coefficient
        
        

        #advance time
        u_prev = u_next
        v_prev = v_next
        p_prev = p_next

    fig, ax = plt.subplots()
    cf = ax.contourf(X, Y, p_next)
    fig.colorbar(cf, ax=ax)
    ax.streamplot(X, Y, u_next, v_next, color='black')
    ax.quiver(X, Y, u_next, v_next, pivot='mid', scale=None)

    ax.add_patch(Circle((x_c, y_c), radius,
                        facecolor='black', edgecolor='black', alpha=1, zorder=5))
    ax.set_xlim(0, domain_size); ax.set_ylim(0, domain_size)
    ax.set_aspect('equal')

    #Calculating vorticity
    omega = central_difference_x(v_next) - central_difference_y(u_next)
    omega_masked = omega.copy()
    omega_masked[circle] = np.nan
    fig_w, ax_w = plt.subplots(figsize=(7, 6))
    cs = ax_w.contourf(X, Y, omega_masked, levels=60, cmap="viridis")
    cbar = fig_w.colorbar(cs, ax=ax_w)
    cbar.set_label("vorticity")
    ax_w.add_patch(Circle((x_c, y_c), radius, color="k"))
    ax_w.set_aspect("equal", adjustable="box")
    ax_w.set_xlim(0, domain_size); ax_w.set_ylim(0, domain_size)

    ax_w.set_title("Vorticity field")
    plt.tight_layout()
   

    plt.show() 
    t = np.cumsum(dt_list)


    t_full = np.cumsum(dt_list)        
                           
    t_seg = t_full[truncation:] - t_full[truncation]      
    sig   = c[truncation:].astype(float)
    sig  -= np.mean(sig)                
                            
    tu = np.arange(0.0, t_seg[-1], 1.0/Fs)
    sigu = np.interp(tu, t_seg, sig) # resample the signal at the uniform time steps for FFT to work
    sigu = detrend(sigu, type='linear')
    w = np.hanning(sigu.size)
    Y = rfft(w*sigu)
    X = rfftfreq(sigu.size, d=1.0/Fs)
    shedding_frequency_index, properties = find_peaks(np.abs(Y), height=height_threshold)
    for i in range(len(shedding_frequency_index)):
        print(f"Shedding frequency: {X[shedding_frequency_index[i]]:0.4f}")
        print(f"Strouhal number: {X[shedding_frequency_index[i]] * 2*radius} ")
    



    # Plot
    plt.figure()
    plt.plot(tu, sigu); plt.xlabel("Time (s)"); plt.ylabel("Signal"); plt.title("Uniform 2 Hz resample")

    plt.figure()
    plt.plot(X, np.abs(Y))
    plt.xlim(0, Fs/2)   # Nyquist = 1 Hz
    plt.xlabel("Frequency (Hz)"); plt.ylabel("Amplitude"); plt.title("FFT @ 2 Hz")
    plt.tight_layout(); plt.show()



if __name__ == "__main__":
    main()

