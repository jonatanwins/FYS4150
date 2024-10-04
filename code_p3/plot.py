import numpy as np
import matplotlib.pyplot as plt


def single_particle(z0, q, m, V0_d2):

    # -------- Analytical Solution ---------
    N = 100
    t_max = 50e-6  # s
    wz = np.sqrt(2 * q / m * V0_d2)

    t_anl = np.linspace(0, t_max, N)
    z_anl = z0 * np.cos(wz * t_anl)

    with open("data/single_particle.txt", "r") as file:
        # length of .txt
        n = len(file.read())
        t = np.zeros(n)
        z = np.zeros(n)
        for i, line in enumerate(file):
            word = line.split()
            t[i] = word[0]
            z[i] = word[1]

    print("hei")
    plt.plot(t_anl, z_anl)
    plt.savefig("single_particle.png")
    plt.show()


if __name__ == "__main__":

    # -------------------------------------- Constants ------------------------------------
    # Coulomb constant
    k_e = 1.38935333e5  # (μm)^3 / (μs)^2 e^2

    # Magnetic field strength (Tesla) and electric potential (Volt)
    T_unit = 9.64852558e1  # μm / (μs * e)
    V_unit = 9.64852558e7  # μm^2 / (μs^2 * e)

    # Penning trap configuration
    B0 = 1.00  # Tesla
    V0 = 25.0e-3  # Volt

    B0_converted = B0 * T_unit  # μm / (μs * e)
    V0_converted = V0 * V_unit  # μm^2 / (μs^2 * e)

    d = 500.0  # μm

    # Ratio of V0/d^2
    V0_d2 = V0_converted / d**2
    # ----------------------------------------------------------------------------------------

    z0 = 20e-6
    q = 1
    m = 1
    single_particle(z0, q, m, V0_d2)
