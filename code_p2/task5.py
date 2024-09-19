import matplotlib.pyplot as plt

N_values = [5, 10, 15, 20, 30, 40, 50, 100, 150]
iterations = [24, 119, 290, 514, 1199, 2163, 3425, 13864, 31777]

# For comparison
quadratic_fit = [N**2 for N in N_values]

# log-log plot with quadratic fit
plt.figure(figsize=(8, 6))
plt.loglog(N_values, iterations, marker="o", linestyle="-", color="b", label="Data")
plt.loglog(N_values, quadratic_fit, "--r", label=r"$N^2$")
plt.title("Log-Log Plot: Iterations vs N with Quadratic Fit")
plt.xlabel("N (log scale)")
plt.ylabel("Iterations (log scale)")
plt.grid(True, which="both", ls="--")
plt.legend()
plt.show()
