import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns

fsize = 18
fsize2 = 15


def read_file(filename):
    with open("data/" + filename, "r") as file:
        # ----------- Extract metadata ------------------
        first_line = file.readline()  # Read first line with metadata
        n_particles, n_timepoints = map(int, first_line.split())  # Split first line

        print(f"Filename: {filename}")
        print(f"n_particles: {n_particles} \nn_timepoints: {n_timepoints} \n")
        # -----------------------------------------------

        t = np.zeros(n_timepoints)
        r = np.zeros((n_particles, n_timepoints, 3))
        v = np.zeros((n_particles, n_timepoints, 3))

        i = 0  # time-index
        j = 0  # particle-index

        for line in file:
            word = line.split()

            if len(word[:]) == 1:  # ------- Time - line
                t[i] = np.float64(word[0])
                j = 0
                i += 1

            else:  # ----------------------- Particles - lines
                r[j, i - 1, :] = np.float64(np.array([word[0], word[1], word[2]]))
                v[j, i - 1, :] = np.float64(np.array([word[3], word[4], word[5]]))
                j += 1

    return t, r, v, n_particles


def single_z_anl(wz, r0, t_max, N):
    z0 = r0[2]
    t = np.linspace(0, t_max, N)
    z = z0 * np.cos(wz * t)
    return t, z


def single_xy_anl(wz, w0, r0, v0, phi_p, phi_m, t_max, N):
    v0 = np.sqrt(v0[0] ** 2 + v0[1] ** 2)
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


def single_particle(r0, v0, w0, wz, phi_m, phi_p):

    # --------------------------- Analytical Solution ---------------------
    N_points = 1001
    t_max = tot_time  # 50e-6  # micro-seconds
    t_anl, z_anl = single_z_anl(wz, r0, t_max, N_points)

    filename = "one_particle_int_RK4.txt"
    # change to: "one_particle_int_FE.txt" for Forward Euler method
    t, r, v, n = read_file(filename)

    # ------------------------------- Plotting motion z -------------------------------
    fig = plt.figure(figsize=(6, 6))
    fig.suptitle("Single particle, z-direction Motion vs. Time", fontsize=fsize)
    ax = fig.add_subplot(1, 1, 1)

    ax.plot(t_anl, z_anl, label=r"$z_{analytical}$")
    ax.plot(t, r[0, :, 2], label=r"$z$", linestyle="--")

    ax.set_xlabel(r"t $[\mu s]$", fontsize=fsize2)
    ax.set_ylabel(r"z $[\mu m]$", fontsize=fsize2)
    ax.legend()

    plt.tight_layout()
    plt.legend()
    plt.grid()
    plt.savefig("plots/one_tz_RK4.pdf")

    # ------------------------------- Plotting Rel. Error -------------------------------
    solver_list = ["FE", "RK4"]
    # if only want RK4 set: solver_list = ["RK4"]
    n_list = [4000, 8000, 16000, 32000]  # timesteps


   

    # Create lists for calculating error convergence rate
    h_list = np.zeros((len(solver_list), len(n_list)))
    delta_list = np.zeros((len(solver_list), len(n_list)))
    r_err_con = np.zeros(len(solver_list))

    for i, solver in enumerate(solver_list):
        fig = plt.figure(figsize=(6, 6))
        fig.suptitle(f"Relative Error, Single particle, {solver}", fontsize=fsize)
        ax = fig.add_subplot(1, 1, 1)

        for j in range(len(n_list)):
            N = n_list[j] + 1
            t_anl, z_anl = single_z_anl(wz, r0, t_max, N)
            t_anl, x_anl, y_anl = single_xy_anl(wz, w0, r0, v0, phi_p, phi_m, t_max, N)
            r_anl = np.array([x_anl, y_anl, z_anl]).T


            filename = "one_particle_no_int_n=" + str(n_list[j]) + "_" + solver + ".txt"
            t, r, v, n = read_file(filename)

            # ----------- Calculate rel. err ----------
            r = r[0, :, :]
            r_err = np.linalg.norm(r - r_anl, axis=1) / (np.linalg.norm(r_anl, axis=1))
            
            ax.plot(t_anl, r_err, label=f"Rel. Err, n = {N-1}")

            h_list[i, j] = 50 / n_list[j]  # micro-seconds
            delta_list[i, j] = max(np.linalg.norm(r - r_anl, axis=1))

        

        ax.set_xlabel(r"t $[\mu s]$", fontsize=fsize2)
        ax.set_ylabel(r"Rel. Error $[]$", fontsize=fsize2)
        ax.legend()

        plt.yscale('log')
        plt.tight_layout()
        plt.legend()
        plt.grid()
        plt.savefig("plots/one_rel_err" + solver + ".pdf")
        plt.show()

        # ------------------------------ Error convergence rate ---------------------------

        for j in range(1, len(n_list)):
            r_err_con[i] += (1 / 3 * np.log10(delta_list[i, j] / delta_list[i, j - 1]) / np.log10(h_list[i, j] / h_list[i, j - 1]))

        print(
            "------------------------- Error Convergence Rate, r_err --------------------"
        )
        print(f"Solver: {solver_list[i]} \n   r_err = {r_err_con[i]:.3e} \n \n \n")


def two_particles(interactions, d, double_xy, double_xv, double_zv, double_xyz):

    solver = "RK4"  #change to "FE" for forward Euler, "RK4" for Runge-kutta 4

    if double_xy or double_xv or double_zv or double_xyz:
        # ------------------------------- Read file -----------------------------
        # For RK4 read:  "two_particles_(interaction)_RK4.txt"
        # For FE  read:  "two_particles_(interaction)_FE.txt"
        if interactions:  # with interactions
            filename = "two_particles_int_" + solver + ".txt" 
            add = "_int"
        else:  # without interactions
            filename = "two_particles_no_int_" + solver + ".txt" 
            add = "_no_int"
        # add to .pdf-name to separate interactinon with non-interaction

        t, r, v, n_particles = read_file(filename)


    # ------------------------------- Plotting -------------------------------
    if double_xy:
        fig = plt.figure(figsize=(7, 7))
        if interactions:
            fig.suptitle("2 Particles, xy-plane, with interaction", fontsize=fsize)
        else:
            fig.suptitle("2 Particles, xy-plane, without interaction", fontsize=fsize)
        ax = fig.add_subplot(1, 1, 1)

        for n in range(n_particles):
            ax.plot(r[n, :, 0], r[n, :, 1], label=f"P{n+1}")

        ax.set_xlabel(r"x $[\mu m]$")
        ax.set_ylabel(r"y $[\mu m]$")
        # ax.set_xlim(-d,d)
        # ax.set_ylim(-d,d)
        ax.legend()
        ax.grid()
        plt.tight_layout()
        plt.savefig("plots/two_xy" + add + "_" + solver + ".pdf")

    if double_xv:
        fig = plt.figure(figsize=(7, 7))
        if interactions:
            fig.suptitle(
                r"2 Particles, $x v_x$-plane, with interaction", fontsize=fsize
            )
        else:
            fig.suptitle(
                r"2 Particles, $x v_x$-plane, without interaction", fontsize=fsize
            )
        ax = fig.add_subplot(1, 1, 1)

        for n in range(n_particles):
            ax.plot(r[n, :, 0], v[n, :, 0], label=f"P{n+1}")

        ax.set_xlabel(r"x $[\mu m]$")
        ax.set_ylabel(r"$v_x$ $[\frac{m}{s}]$")
        ax.legend()
        ax.grid()
        plt.tight_layout()
        plt.savefig("plots/two_xv" + add + "_" + solver + ".pdf")

    if double_zv:
        fig = plt.figure(figsize=(7, 7))
        if interactions:
            fig.suptitle(
                r"2 Particles, $z v_z$-plane, with interaction", fontsize=fsize
            )
        else:
            fig.suptitle(
                r"2 Particles, $z v_z$-plane, without interaction", fontsize=fsize
            )
        ax = fig.add_subplot(1, 1, 1)

        for n in range(n_particles):
            ax.plot(r[n, :, 2], v[n, :, 2], label=f"P{n+1}")

        ax.set_xlabel(r"z $[\mu m]$")
        ax.set_ylabel(r"$v_z$ $[\frac{m}{s}]$")
        ax.legend()
        ax.grid()
        plt.tight_layout()
        plt.savefig("plots/two_zv" + add + "_" + solver + ".pdf")

    if double_xyz:
        fig = plt.figure(figsize=(7, 7))
        if interactions:
            fig.suptitle(r"2 Particles, $xyz$-space, with interaction", fontsize=fsize)
        else:
            fig.suptitle(
                r"2 Particles, $xyz$-space, without interaction", fontsize=fsize
            )

        ax = fig.add_subplot(111, projection="3d")

        for n in range(n_particles):
            ax.plot(r[n, :, 0], r[n, :, 1], r[n, :, 2], label=f"P{n+1}")
            ax.scatter(r[n, 0, 0], r[n, 0, 1], r[n, 0, 2], color="red")

        ax.set_xlabel(r"x $[\mu m]$")
        ax.set_ylabel(r"y $[\mu m]$")
        ax.set_zlabel(r"z $[\mu m]$")
        # ax.set_xlim(-d,d)
        # ax.set_ylim(-d,d)
        # ax.set_zlim(-d,d)

        ax.legend()
        ax.grid()
        plt.tight_layout()
        plt.savefig("plots/two_xyz" + add + "_" + solver + ".pdf")


def loss_particles(interactions, solver):
    count_particles = False
    if count_particles:
        # ------------------------------- Read file -----------------------------
        f_list = np.linspace(0.1, 0.7, 3)
        w_list = np.linspace(2.18, 2.30, 7)
        # print(f_list, w_list)

        #f_list = np.linspace(0.1, 0.7, 3)
        #w_list = np.linspace(0.2, 2.48, 115)
        #print(f_list, w_list)


        file_list = [i for i in range(len(f_list) * len(w_list))]
        # file_list = np.zeros(len(f_list) * len(w_list))

        if interactions:  # with interactions
            add = "int_"
            for i in range(len(f_list)):
                for j in range(len(w_list)):
                    file_list[i * len(w_list) + j] = (
                        add + "f_" + str(f_list[i]) + "_w_v_" + str(w_list[j]) + "_.txt"
                    )

        else:  # without interactions
            add = "no_int_"
            for i in range(len(f_list)):
                for j in range(len(w_list)):
                    file_list[i * len(w_list) + j] = (
                        add + "f_" + str(f_list[i]) + "_w_v_" + str(w_list[j]) + "_.txt"
                    )


        # ----- Sum up trapped particles --------
        trapped_list = np.zeros(len(file_list))
        # trapped_list = np.zeros(len(f_list), len(w_list))
        for i, filename in enumerate(file_list):
            t, r, v, n_particles = read_file(filename)

            n_particles = len(r[:, 0, 0])
            particles = np.zeros(n_particles)  # 0 = outside

            for n in range(n_particles):
                # check if |r(t_max)| <= d
                if np.linalg.norm(r[n, -1, :]) <= d: 
                    particles[n] = 1  # 1 = trapped inside

            trapped_list[i] = sum(particles)



        # ------ Plot histogram ----------
        fig = plt.figure(figsize=(7, 7))
        if interactions:
            fig.suptitle(
                f"Fraction of {n_particles} particles trapped after {tot_time}µs \nwith interaction",
                fontsize=fsize,
            )
        else:
            fig.suptitle(
                f"Fraction of {n_particles} particles trapped after {tot_time}µs \nwithout interaction",
                fontsize=fsize,
            )

        ax = fig.add_subplot(111)

        # Create heatmap
        yedges = f_list  # [0.1, 0.4,] # modify these
        xedges = w_list  # [0.2]      # modify these
        trapped_list = trapped_list.reshape(len(w_list), len(f_list))  # reshape for heatmap
        # trapped_list = trapped_list.reshape(1,-1)

        print("File list:  ", file_list)
        print("Trapped list:  ", trapped_list)
        print("Trapped list shape: ", trapped_list.shape)
        #sns.heatmap(trapped_list, annot=True, cmap="Blues", xticklabels=xedges, yticklabels=yedges)
        sns.heatmap(trapped_list.T, annot=False, cmap="Blues", xticklabels=xedges, yticklabels=yedges, linewidths=0)


        ax.set_xlabel(r"Frequency $\omega_v$ [MHz]")
        ax.set_ylabel(r"Amplitude f $[]$")
        ax.grid(False)
        plt.tight_layout()
        plt.savefig("plots/trapped_" + add + solver + ".pdf")

        







    # ----------------- Plot time evolution -------------------
    n_particles = 100 # is defined above and WILL DELITE WHEN FINISHED

    fig = plt.figure(figsize=(6, 6))
    if interactions:  # with interactions
        add = "int_"
        fig.suptitle(f"Time evolution of {n_particles} trapped particles \nwith interaction", fontsize=fsize)
    else:  # without interactions
        add = "no_int_"
        fig.suptitle(f"Time evolution of {n_particles} trapped particles \nwithout interaction", fontsize=fsize)
    ax = fig.add_subplot(1, 1, 1)


    # ---- Selecting chousen values for time evolution ---
    # [[i1,j1] , [i2,j2]]
    #ij_list = [[1,0],[1,1],[1,2]]
    #ij_list = [[1,2]]


    # --------------- DO IT MANUALLY TO CHECK GIVEN FILES -------------------
    f_list = [0.39, 0.4]
    w_list = [2.239, 2.24]

    ij_list = [[0,0], [1,1]]

    if interactions:  # with interactions
        add = "int_"
    else:  # without interactions
        add = "no_int_"
    # -----------------------------------------------------------------------


    for i, j in ij_list:
        filename = add + "f_" + str(f_list[i]) + "_w_v_" + str(w_list[j]) + "_.txt"  
        t, r, v, n_particles = read_file(filename)

        trapped = np.zeros(len(t))
        n_particles = len(r[:, 0, 0])

        for i_t in range(len(t)):
            particles = np.zeros(n_particles)  # 0 = outside

            for n in range(n_particles):
                    # check if |r(t_max)| <= d
                    if np.linalg.norm(r[n, i_t, :]) <= d: 
                        particles[n] = 1       # 1 = trapped inside

            trapped[i_t] = sum(particles)


        half = int(len(t)/5)
        start = int(len(t)/7)
        t = t[start:half]
        trapped = trapped[start:half]
        ax.plot(t, trapped, label=f'f = {f_list[i]}, w_v = {w_list[j]}')

        
    ax.set_xlabel(r"t $[\mu s]$", fontsize=fsize2)
    ax.set_ylabel(r"Trapped particles $[]$", fontsize=fsize2)
    plt.tight_layout()
    plt.legend()
    plt.grid()
    pdfname = 'plots/trapped_time_evolution' + add + solver + ".pdf"
    plt.savefig(pdfname)


    """
    # ------------------------ Plot xyz motion plot ----------------------------------
    fig = plt.figure(figsize=(7, 7))
    if interactions:
        fig.suptitle(r"10 Particles, $xyz$-space, with interaction", fontsize=fsize)
    else:
        fig.suptitle(r"10 Particles, $xyz$-space, without interaction", fontsize=fsize)

    ax = fig.add_subplot(111, projection="3d")

    for n in range(n_particles):
        r_n = r[n,start:half,:]
        ax.plot(r_n[:,0], r_n[:, 1], r_n[:, 2], label=f"P{n+1}")
        #ax.scatter(r[n, 0, 0], r[n, 0, 1], r[n, 0, 2], "r")
        #ax.scatter(r[n, -1, 0], r[n, -1, 1], r[n, -1, 2], "g")

    # -------------------- Plot sphere -------------------
    phi, theta = np.mgrid[0:2*np.pi:100j, 0:np.pi:50j]

    x = d * np.sin(theta) * np.cos(phi)
    y = d * np.sin(theta) * np.sin(phi)
    z = d * np.cos(theta)
    ax.plot_surface(x, y, z, rstride=5, cstride=5, color='r', alpha=0.4)
    # ------------------------------------------------------

    ax.set_xlabel(r"x $[\mu m]$")
    ax.set_ylabel(r"y $[\mu m]$")
    ax.set_zlabel(r"z $[\mu m]$")
    ax.set_xlim(-d,d)
    ax.set_ylim(-d,d)
    ax.set_zlim(-d,d)

    ax.legend()
    ax.grid()
    plt.tight_layout()
    plt.show()
    plt.savefig("plots/" + str(n_particles) + "_xyz" + add + ".pdf")

    # --------------------------------------------- Plot xy plot ------------------------------------
    fig = plt.figure(figsize=(7, 7))
    if interactions:
        fig.suptitle("10 Particles, xy-plane, with interaction", fontsize=fsize)
    else:
        fig.suptitle("10 Particles, xy-plane, without interaction", fontsize=fsize)
    ax = fig.add_subplot(1, 1, 1)

    for n in range(n_particles):
        if n == 2 or n == 3:
            r_n = r[n,start:half,:]
            ax.plot(r_n[:,0], r_n[:, 1], label=f"P{n+1}")


    # -------------------- Plot sphere -------------------
    phi = np.linspace(0,2*np.pi,100)

    x = d * np.cos(phi) 
    y = d * np.sin(phi)
    ax.plot(x, y, color='r', alpha=0.4)
    # ------------------------------------------------------

    ax.set_xlabel(r"x $[\mu m]$")
    ax.set_ylabel(r"y $[\mu m]$")

    margin = 60
    ax.set_xlim(-d-margin,d+margin)
    ax.set_ylim(-d-margin,d+margin)
    ax.legend()
    ax.grid()
    # ax.axis('equal')
    plt.tight_layout()
    plt.savefig("plots/" + str(n_particles) + "_xy" + add + ".pdf")

    """


if __name__ == "__main__":

    # -------------------------------------- Constants ------------------------------------
    k_e = 1.38935333e5  # (μm)^3 / (μs)^2 e^2
    T_unit = 9.64852558e1  # μm / (μs * e)
    V_unit = 9.64852558e7  # μm^2 / (μs^2 * e)
    d = 500.0  # μm
    e = 1  # electron charge = 1.6e-19 C
    m_p = 40.078  # Proton mass = 1.6726219e-27 kg
    tot_time = 50 

    B0 = 1.00  # Tesla
    V0 = 25.0e-3  # Volt
    B0_converted = B0 * T_unit  # μm / (μs * e)
    V0_converted = V0 * V_unit  # μm^2 / (μs^2 * e)

    B0 = B0_converted
    V0 = V0_converted

    V0_d2 = V0_converted / d**2  # Ratio of V0/d^2

    q = e
    m = m_p
    wz = np.sqrt(2 * q / m * V0_d2)
    w0 = q * B0 / m

    # -------------------- Constants for analytical solution----------------

    # Constants for Analytical Solution
    x0 = 20
    y0 = 0
    z0 = 20
    r0 = [x0, y0, z0]  # micro-meter

    vx0 = 0
    vy0 = 25
    vz0 = 0
    v0 = [vx0, vy0, vz0]  # m/s

    phi_p = 0  # Phase phi+
    phi_m = 0  # Phase phi-

    # --------------------------------------- PLOTS BOOL -----------------------------------
    want_single = True

    double_xy = True
    double_xv = True
    double_zv = True
    double_xyz = True

    interactions = False  # With interaction = True,   Without interactions = False

    want_loss_particles = False

    # ------------------------------------   Cyclotron freq -----------------------------
    cy_freq = q/m * B0 #* 1/(2*np.pi) 

    omega_plus = (w0 + np.sqrt(w0**2 - 2 * wz**2)) / 2
    omega_minus = (w0 - np.sqrt(w0**2 - 2 * wz**2)) / 2

    print(f'Cyclotron Frequency Ca+: freq = ', cy_freq)
    print(f'w_0 = ', w0, 'base units')
    print(f'w_z = ', wz, 'base units')
    print(f'w_minus = ', omega_minus, ' base units')
    print(f'w_plus = ', omega_plus, ' base units')
    print()

    # -------------------------------------   1 Particle --------------------------------
    if want_single:
        single_particle(r0, v0, w0, wz, phi_m, phi_p)

    # -------------------------------------   2 Particles ----------------------------------

    two_particles(interactions, d, double_xy, double_xv, double_zv, double_xyz)

    # -------------------------------- Loss of trapped particles ---------------------------
    solver = "RK4"  #only for pdf-name
    if want_loss_particles:
        loss_particles(interactions, solver)

    # Version 22.10-2024
