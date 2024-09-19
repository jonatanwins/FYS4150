import matplotlib.pyplot as plt

N_values = [5, 10, 15, 20, 30, 40, 50, 100, 150]
iterations_tri = [24, 119, 290, 514, 1199, 2163, 3425, 13864, 31777]  # tridiagonal
iterations_sym = [16, 99, 235, 417, 996, 1827, 2818, 11638, 26510]  # random symmetric

# For comparison
quadratic_fit = [N**2 for N in N_values]

# log-log plot with quadratic fit
plt.figure(figsize=(8, 6))
plt.loglog(
    N_values,
    iterations_tri,
    marker="o",
    linestyle="-",
    color="b",
    label="Rotations (Tridiagonal)",
)
plt.loglog(
    N_values,
    iterations_sym,
    marker="o",
    linestyle="-",
    color="g",
    label="Rotations (Random Symmetric)",
)
plt.loglog(N_values, quadratic_fit, "--r", label=r"$N^2$")
plt.title(r"Log-Log Plot: Rotations per $N$ for $N \times N$ matrices")
plt.xlabel("N (log scale)")
plt.ylabel("Rotations (log scale)")
plt.grid(True, which="both", ls="--")
plt.legend()
plt.savefig("task5_rotations.pdf")
