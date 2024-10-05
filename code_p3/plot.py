import numpy as np
import matplotlib.pyplot as plt


fsize = 15
fsize2 = 10

def read_file(filename):
    with open("data/" + filename, "r") as file:
        # ----------- Extract metadata ------------------    
        first_line = file.readline()  # Read first line with metadata
        n_particles, n_timesteps = map(int, first_line.split()) # Split first line
        
        print(f"Filename: {filename}")
        print(f'n_particles: {n_particles} \nn_timesteps: {n_timesteps} \n')
        # -----------------------------------------------

        t = np.zeros(n_timesteps)
        r = np.zeros((n_particles, n_timesteps, 3))
        v = np.zeros((n_particles, n_timesteps, 3))
        
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
    



def single_particle(z0, q, m, V0_d2): #tz-plot

    # --------------------------- Analytical Solution ---------------------
    N = 100
    t_max = 50e-6  # s
    wz = np.sqrt(2 * q / m * V0_d2)

    t_anl = np.linspace(0, t_max, N)
    z_anl = z0 * np.cos(wz * t_anl)

   

    filename = "one_particle_int.txt"
    t, r, v, n = read_file(filename)

    # ------------------------------- Plotting -------------------------------
    fig = plt.figure(figsize=(7,7))
    fig.suptitle('Single particle motin in time, z-direction in time', fontsize = fsize)
    ax = fig.add_subplot(1,1,1)

    # plot in micro-seconds/meters
    #t, z, t_anl, z_anl = t*1e6, z*1e6, t_anl*1e6, z_anl*1e6
    ax.plot(t_anl, z_anl, label=r'$z_{analytical}$')
    ax.plot(t    , r[0,:,2]    , label=r'$z}$' )

    ax.set_xlabel(r't $[\mu s]$') 
    ax.set_ylabel(r'z $[\mu m]$')
    ax.legend()


    plt.tight_layout()
    plt.legend()
    plt.grid()
    plt.savefig("plots/one_tz.pdf")


def two_particles(interactions, double_xy, double_xv, double_zv, double_xyz):
    # ------------------------------- Read file -----------------------------
    if interactions:  # with interactions
        #filename = "two_particle_int.txt"
        filename = "two_particles_int.txt"
        add = "_int" 
    else: # without interactions
        filename = "two_particles_no_int.txt"
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
            ax.plot(r[n,:,0], r[n,:,1], v[n,:,2], label=f'P{n+1}') 

    
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
    B0 = 1.00  # Tesla
    V0 = 25.0e-3  # Volt
    B0_converted = B0 * T_unit  # μm / (μs * e)
    V0_converted = V0 * V_unit  # μm^2 / (μs^2 * e)
    d = 500.0  # μm
    e = 1.6e-19 # electron charge C
    m_p = 1.6726219e-27 # Proton mass kg
    V0_d2 = V0_converted / d**2 # Ratio of V0/d^2


    # --------------------------------------- PLOTS BOOL -----------------------------------
    single     = True
    double_xy  = True
    double_xv  = True
    double_zv  = True
    double_xyz = True

    interactions = True  # With interaction = True,   Without interactions = False



    # -------------------------------------   1 Particle --------------------------------
    if single: 
        # Parameters for analytical solution
        z0 = 20 # micro-meter 
        q = e
        m = m_p
        # ----------------------------------

        single_particle(z0, q, m, V0_d2)



    # -------------------------------------   2 Particles ----------------------------------

    two_particles(interactions, double_xy, double_xv, double_zv, double_xyz)
