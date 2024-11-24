#hard beepboop of mozerfoking maphys bb ;)


import matplotlib.pyplot as plt
import numpy as np

# Constante
gamma_i = 0.064
gamma_x = 0.004
lambda_i = 2.926e-5
lambda_x = 2.0996e-5
sigma_i = 7e-24
sigma_x = 2.65e-18
maj_sigma_f = 0.0984
maj_sigma_u = 0.1355
psi_de_t = 3e13
end_time = 178200  

# vu que l'équa diff pour X depend aussi de I on fait euler mais pour un systeme d'équation
def euler_coupled(derivatives, t0, y0, t_end, dt):
    t = t0
    y = np.array(y0)
    times = [t]
    values = [y]
    while t < t_end:
        y = y + dt * np.array(derivatives(t, y))
        t += dt
        times.append(t)
        values.append(y)
    return np.array(times), np.array(values)

# systeme d(équation diff
def derivatives(t, y):
    I, grandX = y  
    dI_dt = gamma_i * maj_sigma_f * psi_de_t - lambda_i * I - sigma_i * I * psi_de_t
    dX_dt = gamma_x * maj_sigma_f * psi_de_t + lambda_i * I - lambda_x * grandX - sigma_x * grandX * psi_de_t
    return [dI_dt, dX_dt]

# Simulation
def run_simulation():
    t0 = 0
    y0 = [0, 0]  
    dt = 900
    times, results = euler_coupled(derivatives, t0, y0, end_time, dt)
    return times, results[:, 0], results[:, 1] 

# Plot des fct
def plot_results(times, Is, Xs):
    plt.plot(times / 3600, Is, label="abondence en Iode")
    plt.plot(times / 3600, Xs, label="abondence en Xénon")
    plt.xlabel("Temps en heure ")
    plt.ylabel("nombre de particule")
    plt.title("Simulation de I et X par rapport au temps")
    plt.legend()
    plt.grid()
    plt.show()

if __name__ == "__main__":
    times, Is, Xs = run_simulation()
    plot_results(times, Is, Xs)
