import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm 


N = 41 # 41 points per axis
domain_size = 1
iterations = 500
time_step_length = .001
kinematic_viscosity = .1  #Changing this value can affect stability of simulation, i.e 0.1 to 0.2. 
density = 1.0
top_velocity = 1.0 
pressure_poisson_iterations = 50

def main():
    element_length = domain_size / (N - 1) #initialize the grid
    x = np.linspace(0, domain_size, N)
    y = np.linspace(0, domain_size, N)
    X,Y = np.meshgrid(x,y)

    u_prev = np.zeros_like(X)           # initialize initial conditions of fluid
    v_prev = np.zeros_like(X)
    p_prev = np.zeros_like(X)

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
    
    for _ in tqdm(range(iterations)): 
        du_prev_dx = central_difference_x(u_prev)
        du_prev_dy = central_difference_y(u_prev)
        dv_prev_dx = central_difference_x(v_prev)
        dv_prev_dy = central_difference_y(v_prev)
        laplace_u_prev = laplace(u_prev)
        laplace_v_prev = laplace(v_prev)

        #Solve momentum equation w/o pressure consideration using Chorin's projection method
        u = (u_prev + time_step_length * ( - (u_prev * du_prev_dx + v_prev * du_prev_dy) + kinematic_viscosity * laplace_u_prev ) )
        v = (v_prev + time_step_length * (- (u_prev * dv_prev_dx + v_prev * dv_prev_dy) + kinematic_viscosity * laplace_v_prev))
        
        #Velocity BC: homogenous dirichlect BC everywhere 
        # except for horizontal velocity at top, applied here:

        u[0, :] = 0
        u[:, 0] = 0
        u[:, -1] = 0
        u[-1, :] = top_velocity
        v[0, :] = 0.0
        v[:, 0] = 0
        v[:, -1] = 0
        v[-1, :] = 0
    
        du_dx = central_difference_x(u)
        dv_dy = central_difference_y(v)

        #Compute pressure correction by solving pressure-poisson equation

        rhs = ( density / time_step_length * (du_dx + dv_dy))

        for _ in range(pressure_poisson_iterations):
            p_next = np.zeros_like(p_prev)
            p_next[1:-1, 1:-1] = 1/4 * (+ p_prev[1:-1, 0:-2] + p_prev[0: -2, 1:-1] + p_prev[1:-1, 2:] + p_prev[2:, 1:-1] - element_length**2 * rhs[1: -1, 1:-1])

            #pressure boundary conditions: homogenous neumann BC everywhere
            #except for top, homogenous dirichlect BC given here:

            p_next[:, -1] = p_next[:, -2]
            p_next[0, :] = p_next[1, :]
            p_next[:, 0] = p_next[:, 1]
            p_next[-1, :] = 0

            p_prev = p_next

        #Now we apply the pressure correction to velocity to keep fluid incompressible
        dp_next_dx = central_difference_x(p_next)
        dp_next_dy = central_difference_y(p_next)

        u_next = (u - time_step_length / density * dp_next_dx)
        v_next = (v - time_step_length / density * dp_next_dy)

        #Reapply boundary conditions

        u_next[0, :] = 0
        u_next[:, 0] = 0
        u_next[:, -1] = 0
        u_next[-1, :] = top_velocity
        v_next[0, :] = 0.0
        v_next[:, 0] = 0
        v_next[:, -1] = 0
        v_next[-1, :] = 0

        #advance time
        u_prev = u_next
        v_prev = v_next
        p_prev = p_next

    plt.figure()
    plt.contourf(X,Y, p_next)
    plt.colorbar()
   # plt.quiver(X, Y, u_next, v_next, color = "black")
    plt.streamplot(X, Y, u_next, v_next, color='black')
    plt.show()


if __name__ == "__main__":
    main()

