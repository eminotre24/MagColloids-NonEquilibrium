import numpy as np
import matplotlib.pyplot as plt

# Parameters
n_simulations = 100  # Number of simulations
n_steps = 1000       # Number of steps in each simulation
dt = 0.01            # Time step
k = 1.0              # Harmonic potential constant
gamma = 0.5          # Friction coefficient (damping)
T = 1.0              # Temperature
beta = 1 / T         # Inverse temperature (Boltzmann constant set to 1)
x0 = 0.0             # Initial position
x_mean = 0.0         # Mean of the stationary distribution
sigma = np.sqrt(T / k)  # Standard deviation of the stationary distribution

# Precomputing constants for the Langevin equation
sqrt_2gammaT_dt = np.sqrt(2 * gamma * T * dt)

def simulate_brownian_harmonic(n_steps, dt, k, gamma, T, x0):
    """Simulate a single Brownian motion trajectory under a harmonic potential."""
    x = np.zeros(n_steps)
    x[0] = x0

    for i in range(1, n_steps):
        # Harmonic force term: -k*x
        force = -k * x[i-1]
        
        # Langevin equation: dx = force * dt + noise
        noise = np.random.normal(0, sqrt_2gammaT_dt)
        x[i] = x[i-1] + (force * dt / gamma) + noise

    return x

def approx(x,val,epsilon):
    if x > val - epsilon and x < val + epsilon:
        return True
    else:
        return False

# Run n simulations
trajectories = np.zeros((n_simulations, n_steps))

for i in range(n_simulations):
    trajectories[i] = simulate_brownian_harmonic(n_steps, dt, k, gamma, T, x0)

# Plotting the results
time = np.arange(n_steps) * dt

plt.figure(figsize=(10, 6))

# Plot a few sample trajectories
for i in range(n_simulations):  # Plot 5 trajectories or fewer
    plt.plot(time, trajectories[i], label='', alpha = 0.5, color='r')

plt.title(f'{n_simulations} Simulations of Brownian Motion under Harmonic Potential')
plt.xlabel('Time')
plt.ylabel('Position (x)')
plt.legend()
plt.grid()
plt.show()

# Plot the distribution of final positions
plt.figure(figsize=(8, 5))
plt.hist(trajectories[:, -1], bins=30, density=True, alpha=0.7, color='b', label='Final positions')
plt.axvline(x=x_mean, color='r', linestyle='--', label='Mean position')
plt.title('Distribution of Final Positions')
plt.xlabel('x')
plt.ylabel('Probability Density')
plt.legend()
plt.grid()
plt.show()

val1 = 1 #At frm1
val2 = 0 #At frm2
frm1 = 125
frm2 = 250
ep = 0.25 #epsilon

for i in range(n_simulations):
    s = trajectories[i]
    d = s[0]
    sn = [x - d for x in s]
    if approx(sn[frm1], val1, ep) and approx(sn[frm2], val2, ep):
        plt.plot(sn, color = "k", alpha=0.1)
plt.errorbar(frm1, val1, yerr=ep, fmt="o", capsize=5, label = "point 1")
plt.errorbar(frm2, val2, yerr=ep, fmt="o", capsize=5, label = "point 2")

plt.title("Trayectoria de las N particulas Normalizada* Filtradas")
plt.xlabel("Num de Frame")
plt.ylabel(r"x $\mu{}m$")
plt.legend()
plt.grid()
plt.show()