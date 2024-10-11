import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


fsize = 18
fsize2 = 15

def read_file(filename):
    with open("data/" + filename, "r") as file:
        # ----------- Extract metadata ------------------    
        first_line = file.readline()  # Read first line with metadata
        n_particles, n_timepoints = map(int, first_line.split()) # Split first line
        
        print(f"Filename: {filename}")
        print(f'n_particles: {n_particles} \nn_timepoints: {n_timepoints} \n')
        # -----------------------------------------------

        t = np.zeros(n_timepoints)
        r = np.zeros((n_particles, n_timepoints, 3))
        v = np.zeros((n_particles, n_timepoints, 3))
        
        i = 0 # time-index
        j = 0 # particle-index

        for line in file:
            word = line.split()

            if len(word[:]) == 1: #------- Time - line
                t[i] = np.float64(word[0])
                j = 0
                i += 1

            else: #----------------------- Particles - lines
                r[j,i-1,:] = np.float64(np.array( [word[0],word[1],word[2]] ))
                v[j,i-1,:] = np.float64(np.array( [word[3],word[4],word[5]] ))
                j += 1

    return t, r, v, n_particles
    
def single_z_anl(wz,r0,t_max,N):
    z0 = r0[2]
    t= np.linspace(0, t_max, N) 
    z = z0 * np.cos(wz * t)
    return t, z

def single_xy_anl(wz,w0,r0,v0, phi_p,phi_m,t_max,N):
    v0 = np.sqrt(v0[0]**2 + v0[1]**2)
    x0 = r0[0]
    t = np.linspace(0, t_max, N) 

    omega_plus = (w0 + np.sqrt(w0**2 - 2 * wz**2)) / 2
    omega_minus = (w0 - np.sqrt(w0**2 - 2 * wz**2)) / 2

    # Define A_plus and A_minus based on the provided equations
    A_plus = (v0 + omega_minus * x0) / (omega_minus - omega_plus)
    A_minus = -(v0 + omega_plus * x0) / (omega_minus - omega_plus)

    term1 = A_plus * np.exp(-1j * (omega_plus * t + phi_p))
    term2 = A_minus * np.exp(-1j * (omega_minus * t + phi_m))

    f = term1 + term2
    x = np.real(f)
    y = np.imag(f)
    return t, x, y


def single_particle(r0,v0,w0,wz,phi_m,phi_p):

    # --------------------------- Analytical Solution ---------------------
    N_points = 1001
    t_max = tot_time #50e-6  # s
    t_anl, z_anl = single_z_anl(wz, r0, t_max, N_points)


    filename = "one_particle_int_RK4_time_dep.txt"
    t, r, v, n = read_file(filename)

    # ------------------------------- Plotting motion z -------------------------------
    fig = plt.figure(figsize=(6,6))
    fig.suptitle('Single particle motin in time, z-direction in time', fontsize = fsize)
    ax = fig.add_subplot(1,1,1)

    # plot in micro-seconds/meters
    #t, z, t_anl, z_anl = t*1e6, z*1e6, t_anl*1e6, z_anl*1e6
    ax.plot(t_anl, z_anl, label=r'$z_{analytical}$')
    ax.plot(t, r[0,:,2], label=r'$z$', linestyle='--')

    ax.set_xlabel(r't $[\mu s]$', fontsize = fsize2) 
    ax.set_ylabel(r'z $[\mu m]$', fontsize = fsize2)
    ax.legend()

    plt.tight_layout()
    plt.legend()
    plt.grid()
    plt.savefig("plots/one_tz.pdf")


    # ------------------------------- Plotting Rel. Error -------------------------------
    fig = plt.figure(figsize=(6,6))
    fig.suptitle('Relative Error, Single particle', fontsize = fsize)
    ax = fig.add_subplot(1,1,1)

    # ----- Calculate rel. err -------
    n_list = [4000, 8000, 16000, 32000] # timesteps
    for i in range(len(n_list)):
        N = n_list[i]+1 #should be +1, fix
        t_anl, z_anl = single_z_anl(wz,r0,t_max,N)
        t_anl, x_anl, y_anl = single_xy_anl(wz,w0,r0,v0,phi_p,phi_m,t_max,N)
        r_anl = np.array([x_anl,y_anl,z_anl]).T

        filename = "one_particle_no_int_n=" + str(n_list[i]) + "_RK4.txt"
        t, r, v, n = read_file(filename)

        eps   = 1e-15 # buffer - not divide by zero
        r = r[0,:,:]
        r_err = np.linalg.norm(r-r_anl, axis=1) / (np.linalg.norm(r_anl, axis=1) + eps)

        ax.plot(t_anl, r_err, label=f'Rel. Err, n = {N-1}')


    # plot in micro-seconds/meters
    #ax.plot(t_anl, z_anl, label=r'$z_{analytical}$')
    #ax.plot(t    , r[0,:,2]    , label=r'$z$' )

    ax.set_xlabel(r't $[\mu s]$', fontsize = fsize2) 
    ax.set_ylabel(r'Rel. Error $[]$', fontsize = fsize2)
    ax.legend()

    plt.tight_layout()
    plt.legend()
    plt.grid()
    plt.savefig("plots/one_rel_err.pdf")
    plt.show()


def two_particles(interactions, double_xy, double_xv, double_zv, double_xyz):
    # ------------------------------- Read file -----------------------------
    if interactions:  # with interactions
        #filename = "two_particle_int.txt"
        filename = "two_particles_int_RK4.txt"
        add = "_int" 
    else: # without interactions
        filename = "two_particles_no_int_RK4.txt"
        add = "_no_int"
    # add to .pdf-name to separate interactinon with non-interaction  
    

    t, r, v, n_particles = read_file(filename)



    # ------------------------------- Plotting -------------------------------
    if double_xy:
        fig = plt.figure(figsize=(7,7))
        fig.suptitle('2 Particles, xy-plane', fontsize = fsize)
        ax = fig.add_subplot(1,1,1)

        for n in range(n_particles):
            ax.plot(r[n,:,0],r[n,:,1], label=f'P{n+1}')

        ax.set_xlabel(r'x $[\mu m]$') 
        ax.set_ylabel(r'y $[\mu m]$')
        ax.legend()

        ax.legend()
        ax.grid()
        plt.tight_layout()
        plt.savefig("plots/two_xy" + add + ".pdf")

    
    if double_xv:
        fig = plt.figure(figsize=(7,7))
        fig.suptitle(r'2 Particles, $x v_x$-plane', fontsize = fsize)
        ax = fig.add_subplot(1,1,1)

        # plot in micro-seconds/meters
        #t, z, t_anl, z_anl = t*1e6, z*1e6, t_anl*1e6, z_anl*1e6

        for n in range(n_particles):
            ax.plot(r[n,:,0],v[n,:,0], label=f'P{n+1}')

        ax.set_xlabel(r'x $[\mu m]$') 
        ax.set_ylabel(r'$v_x$ $[\frac{m}{s}]$')
        ax.legend()

        ax.legend()
        ax.grid()
        plt.tight_layout()
        plt.savefig("plots/two_xv" + add + ".pdf")


    if double_zv:
        fig = plt.figure(figsize=(7,7))
        fig.suptitle(r'2 Particles, $z v_z$-plane', fontsize = fsize)
        ax = fig.add_subplot(1,1,1)

        # plot in micro-seconds/meters
        #t, z, t_anl, z_anl = t*1e6, z*1e6, t_anl*1e6, z_anl*1e6

        for n in range(n_particles):
            ax.plot(r[n,:,2],v[n,:,2], label=f'P{n+1}')

        ax.set_xlabel(r'z $[\mu m]$') 
        ax.set_ylabel(r'$v_z$ $[\frac{m}{s}]$')
        ax.legend()

        ax.legend()
        ax.grid()
        plt.tight_layout()
        plt.savefig("plots/two_zv" + add + ".pdf")


    if double_xyz:
        from mpl_toolkits.mplot3d import Axes3D

        fig = plt.figure(figsize=(7, 7))
        fig.suptitle(r'2 Particles, $xyz$-space', fontsize=fsize)
        ax = fig.add_subplot(111, projection='3d')

        for n in range(n_particles):            
            ax.plot(r[n,:,0], r[n,:,1], r[n,:,2], label=f'P{n+1}') 

    
        ax.set_xlabel(r'x $[\mu m]$')
        ax.set_ylabel(r'y $[\mu m]$')
        ax.set_zlabel(r'z $[\mu m]$')

        ax.legend()
        ax.grid()
        plt.tight_layout()
        plt.savefig("plots/two_xyz" + add + ".pdf")



if __name__ == "__main__":

    # -------------------------------------- Constants ------------------------------------
    k_e = 1.38935333e5  # (μm)^3 / (μs)^2 e^2
    T_unit = 9.64852558e1  # μm / (μs * e)
    V_unit = 9.64852558e7  # μm^2 / (μs^2 * e)
    d = 500.0  # μm
    e = 1.0
    m_p = 1.0
    tot_time = 500

    B0 = 1.00  # Tesla
    V0 = 25.0e-3  # Volt
    B0_converted = B0 * T_unit  # μm / (μs * e)
    V0_converted = V0 * V_unit # Volt
  # μm^2 / (μs^2 * e)
    #B0 = B0_converted
    #V0 = V0_converted
    V0_d2 = V0_converted / d**2 # Ratio of V0/d^2

    q = e
    m = 40.078
    wz = np.sqrt(2 * q / m * V0_d2)
    w0 = q * B0_converted / m

    # -------------------- Constants for analytical solution----------------
    
    
    # Constants for Analytical Solution
    x0 = 20
    y0 = 0
    z0 = 20
    r0 = [x0,y0,z0] # micro-meter 

    vx0 = 0
    vy0 = 25
    vz0 = 0
    v0  = [vx0,vy0,vz0] # m/s
    
    phi_p = 0  # Phase phi+
    phi_m = 0  # Phase phi-
    

    # --------------------------------------- PLOTS BOOL -----------------------------------
    single     = True
    double_xy  = True 
    double_xv  = True
    double_zv  = True
    double_xyz = True

    interactions = False  # With interaction = True,   Without interactions = False



    # -------------------------------------   1 Particle --------------------------------
    if single: 
        single_particle(r0,v0,w0,wz,phi_m,phi_p)



    # -------------------------------------   2 Particles ----------------------------------

    # two_particles(interactions, double_xy, double_xv, double_zv, double_xyz)
