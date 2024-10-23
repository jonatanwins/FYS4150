import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
import pandas as pd

fsize = 20 # Fontsize for title
fsize2 = 15 # Fontsize for labels



########################################################################

def read_file(filename):
    with open("data/" + filename, "r") as file:
        # ----------- Extract metadata ------------------
        first_line = file.readline()  # Read first line with metadata
        n_particles, n_timepoints = map(int, first_line.split())  # Split first line

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

########################################################################



# single particle analytical solution in z-direction
def single_z_anl(wz, r0, t_max, N):
    z0 = r0[2]
    t = np.linspace(0, t_max, N)
    z = z0 * np.cos(wz * t)
    return t, z

# single particle analytical solution in xy-plane
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


def single_particle(r0, v0, w0, wz, phi_m, phi_p, tot_time):

    # --------------------------- Analytical Solution ---------------------
    N_points = 1001
    t_max = tot_time  # micro-seconds
    t_anl, z_anl = single_z_anl(wz, r0, t_max, N_points)

    filename = "one_particle_int_RK4.txt"
    # change to: "one_particle_int_FE.txt" for Forward Euler method
    t, r, v, n = read_file(filename)

    # ------------------------------- Plotting motion z -------------------------------
    fig = plt.figure(figsize=(7,7))
    fig.suptitle("One particle, z-position over time, RK4", fontsize=fsize)
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
        fig = plt.figure(figsize=(7,7))
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


def two_particles():

    solvers = ["RK4", "FE"]  
    interactions_list = [True, False]
    

    for solver in solvers:
        for interactions in interactions_list:

        
            if interactions:  # with interactions
                filename = "two_particles_int_" + solver + ".txt" 
                add = "_int"
            else:  # without interactions
                filename = "two_particles_no_int_" + solver + ".txt" 
                add = "_no_int"
            # add to .pdf-name to separate interactinon with non-interaction

            t, r, v, n_particles = read_file(filename)

        # ------------------------------- Plotting -------------------------------

            # xy plane
            fig = plt.figure(figsize=(7, 7))
            if interactions:
                fig.suptitle(f"2 Particles, xy-plane, with interaction, {solver}", fontsize=fsize)
            else:
                fig.suptitle(f"2 Particles, xy-plane, without interaction, {solver}", fontsize=fsize)
            ax = fig.add_subplot(1, 1, 1)

            for n in range(n_particles):
                ax.plot(r[n, :, 0], r[n, :, 1], label=f"P{n+1}")

            ax.set_xlabel(r"x $[\mu m]$")
            ax.set_ylabel(r"y $[\mu m]$")
            ax.legend()
            ax.grid()
            plt.tight_layout()
            plt.savefig("plots/two_xy" + add + "_" + solver + ".pdf")

            # x and v_x
            fig = plt.figure(figsize=(7, 7))
            if interactions:
                fig.suptitle(
                    fr"2 Particles, $x v_x$-plane, with interaction, {solver}", fontsize=fsize
                )
            else:
                fig.suptitle(
                    fr"2 Particles, $x v_x$-plane, without interaction, {solver}", fontsize=fsize
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

            # z and v_z
            fig = plt.figure(figsize=(7, 7))
            if interactions:
                fig.suptitle(
                    fr"2 Particles, $z v_z$-plane, with interaction, {solver}", fontsize=fsize
                )
            else:
                fig.suptitle(
                    fr"2 Particles, $z v_z$-plane, without interaction, {solver}", fontsize=fsize
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

            # xyz-space
            fig = plt.figure(figsize=(7, 7))
            if interactions:
                fig.suptitle(fr"2 Particles, $xyz$-space, with interaction, {solver}", fontsize=fsize)
            else:
                fig.suptitle(fr"2 Particles, $xyz$-space, without interaction, {solver}", fontsize=fsize)

            ax = fig.add_subplot(111, projection="3d")

            for n in range(n_particles):
                ax.plot(r[n, :, 0], r[n, :, 1], r[n, :, 2], label=f"P{n+1}")
                ax.scatter(r[n, 0, 0], r[n, 0, 1], r[n, 0, 2], color="red")

            ax.set_xlabel(r"x $[\mu m]$")
            ax.set_ylabel(r"y $[\mu m]$")
            ax.set_zlabel(r"z $[\mu m]$")

            ax.legend()
            ax.grid()
            plt.tight_layout()
            plt.savefig("plots/two_xyz" + add + "_" + solver + ".pdf")


def grid_search_loss():
    f_list = np.array([0.1, 0.4, 0.7])
    w_v_list = np.array([0.2, 0.22, 0.24, 0.26, 0.28, 0.3, 0.32, 0.34, 0.36, 0.38, 0.4, 0.42, 0.44, 0.46, 0.48, 0.5, 0.52, 0.54, 0.56, 0.58,
                        0.6, 0.62, 0.64, 0.66, 0.68, 0.7, 0.72, 0.74, 0.76, 0.78, 0.8, 0.82, 0.84, 0.86, 0.88, 0.9, 0.92, 0.94, 0.96, 0.98,
                        "1", 1.02, 1.04, 1.06, 1.08, 1.1, 1.12, 1.14, 1.16, 1.18, 1.2, 1.22, 1.24, 1.26, 1.28, 1.3, 1.32, 1.34, 1.36, 1.38,
                        1.4, 1.42, 1.44, 1.46, 1.48, 1.5, 1.52, 1.54, 1.56, 1.58, 1.6, 1.62, 1.64, 1.66, 1.68, 1.7, 1.72, 1.74, 1.76, 1.78,
                        1.8, 1.82, 1.84, 1.86, 1.88, 1.9, 1.92, 1.94, 1.96, 1.98, "2", 2.02, 2.04, 2.06, 2.08, 2.1, 2.12, 2.14, 2.16, 2.18,
                        2.2, 2.22, 2.24, 2.26, 2.28, 2.3, 2.32, 2.34, 2.36, 2.38, 2.4, 2.42, 2.44, 2.46, 2.48, 2.5])
    w_v_list = w_v_list.astype(str)


    trapped = [] 

    for f in f_list:
        for w_v in w_v_list:
            filename = f"no_int_f_{f}_w_v_{w_v}_.txt"
            t, r, v, n_particles = read_file(filename)  

            particles = np.zeros(n_particles)  # Array to track trapped particles (0 = outside)
            for n in range(n_particles): 
                if np.linalg.norm(r[n, :]) <= 500.0:
                    particles[n] = 1  # Particle is trapped inside

            info = {"f": f, "w_v": w_v, "trapped": np.sum(particles)}
            trapped.append(info)

    data = {"f": [], "w_v": [], "trapped": []}
    for entry in trapped:
        data["f"].append(entry["f"])
        data["w_v"].append(entry["w_v"])
        data["trapped"].append(entry["trapped"])

    # Create a DataFrame from the dictionary
    df = pd.DataFrame(data)

    df_pivot = df.pivot("f", "w_v", "trapped")

    plt.figure(figsize=(6,7))
    ax = sns.heatmap(df_pivot, annot=False, cmap="crest")
    plt.title(f"Number of trapped particles after \n 500 μs  without interactions", fontsize=fsize)
    plt.xlabel(r"$\omega_V$", fontsize=fsize2)
    plt.ylabel("f", fontsize=fsize2)

    plt.xticks(fontsize=fsize)
    plt.yticks(fontsize=fsize)

    # Adjust number of ticks on x and y axes
    plt.locator_params(axis='x', nbins=15)
    plt.locator_params(axis='y', nbins=10)

    plt.savefig("plots/heatmap_trapped_particles.pdf")


def time_evolution(w_v_middle, interactions = True):

    w_file = np.array([w_v_middle-0.02, w_v_middle-0.01, w_v_middle, w_v_middle+0.01, w_v_middle+0.02])
    w_file = np.round(w_file, 2)  # Round to 2 decimal places
    w_file = w_file.astype(str)

    all_particles = []
        
    for w in w_file:
        if interactions:
            t, r, v, n_particles = read_file(f"int_f_0.4_w_v_{w}_.txt")
        else:
            t, r, v, n_particles = read_file(f"no_int_f_0.4_w_v_{w}_fine.txt")
        t = t[0:51]
        trapped = np.zeros(len(t))
        for i_t in range(len(t)):
            particles = np.zeros(n_particles)  # 0 = outside
            for n in range(n_particles): 
                # check if |r(t_max)| <= d 
                if np.linalg.norm(r[n, i_t, :]) <= 500.0:
                    particles[n] = 1       # 1 = trapped inside
            trapped[i_t] = sum(particles)
        all_particles.append(trapped)

    plt.figure(figsize=(7, 7))

    for i, w in enumerate(w_file):
        plt.plot(t, all_particles[i], label=fr"$\omega_v = {w}$")

    plt.xlabel(r"$Time [\mu s]$", fontsize=fsize2)
    plt.ylabel("Trapped particles", fontsize=fsize2)
    title_info = "with" if interactions else "without"
    plt.title(f"Time evolution of the number of trapped \n particles, {title_info} interactions", fontsize=fsize)
    plt.legend(fontsize=fsize2)
    plt.grid()
    plt.savefig(f"plots/time_evolution_trapped_around_{w_v_middle}.pdf")