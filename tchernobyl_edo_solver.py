#hard coded maphys things bb 

import matplotlib.pyplot as plt
import numpy as np

# Constants
gamma_i = 0.064
gamma_x = 0.004
lambda_i = 2.926e-5
lambda_x = 2.0996e-5
sigma_i = 7e-24
sigma_x = 2.65e-18
maj_sigma_f = 0.0984
maj_sigma_u = 0.1355
psi_de_t = 3e13
end_time = 178200  # seconds

def maj_sigma_x(sigma_x, grandX):
    return sigma_x * grandX

# Euler Method
def euler(derivative, t0, y0, t_end, dt):
    t = t0
    y = y0
    times = [t]
    values = [y]
    while t < t_end:
        y = y + dt * derivative(t, y)
        t += dt
        times.append(t)
        values.append(y)
    return np.array(times), np.array(values)

# Differential Equations
def dérivé_I(t, I):
    return gamma_i * maj_sigma_f * psi_de_t - lambda_i * I - sigma_i * I * psi_de_t

def dérivé_grandX(t, grandX):
    return gamma_x * maj_sigma_f * psi_de_t - lambda_x * grandX - sigma_x * grandX * psi_de_t

# Simulations
def run_simuI():
    t0 = 0
    I0 = 0  # Initial condition for I
    dt = 3600  # Time step in seconds
    times, I = euler(dérivé_I, t0, I0, end_time, dt)
    return times, I

def run_simugrandX():
    t0 = 0
    X0 = 0  # Initial condition for X
    dt = 3600  # Time step in seconds
    times, X = euler(dérivé_grandX, t0, X0, end_time, dt)
    return times, X

# Plotting Functions
def plot_simuI(times, Is):
    plt.plot(times / 3600, Is, label="I (particles)")
    plt.xlabel("Time (hours)")
    plt.ylabel("I")
    plt.title("Simulation of I over time")
    plt.legend()
    plt.grid()
    plt.show()

def plot_simugrandX(times, Xs):
    plt.plot(times / 3600, Xs, label="X (particles)")
    plt.xlabel("Time (hours)")
    plt.ylabel("X")
    plt.title("Simulation of X over time")
    plt.legend()
    plt.grid()
    plt.show()

if __name__ == "__main__":
    times_I, Is = run_simuI()
    times_X, Xs = run_simugrandX()
    plot_simuI(times_I, Is)
    plot_simugrandX(times_X, Xs)
