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
end_time = 178200  
sigma_b_min = 0.1  # Barres relevées
sigma_b_max = 0.2  # Barres enfoncées
flux_target = 1e15  # Flux cible
flux_initial = 1e10  # Flux initial


# vu que l'équa diff pour X depend aussi de I on fait euler explicite mais pour un systeme d'équation diff
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

# systeme d'équation diff
def derivatives(t, y ):
    def adjust_bars(psi_de_t, sigma_b):
        if psi_de_t < flux_target * 0.95:  # Si le flux est trop bas
            sigma_b = max(sigma_b_min, sigma_b - 0.001)  # Relever les barres
        elif psi_de_t > flux_target * 1.05:  # Si le flux est trop haut
            sigma_b = min(sigma_b_max, sigma_b + 0.001)  # Enfoncer les barres
        return sigma_b
    I, grandX , psi_de_t = y  
    maj_sigma_x = sigma_x * grandX 
    dI_dt = gamma_i * maj_sigma_f * psi_de_t - lambda_i * I - sigma_i * I * psi_de_t
    dX_dt = gamma_x * maj_sigma_f * psi_de_t + lambda_i * I - lambda_x * grandX - sigma_x * grandX * psi_de_t
    dpsi_dt = (psi_de_t/1000)*3*(maj_sigma_u - adjust_bars(psi_de_t, sigma_b ) - maj_sigma_x)
    return [dI_dt, dX_dt, dpsi_dt]

# Simulation
def run_simulation():
    t0 = 0
    y0 = [0, 0, 10e10]  
    dt = 900
    times, results = euler_coupled(derivatives, t0, y0, end_time, dt)
    return times, results[:, 0], results[:, 1] 

# Plot de I et X en fct du temps 
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