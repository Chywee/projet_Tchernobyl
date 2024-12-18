#hard beepboop of mozerfoking maphys bb ;)
import matplotlib.pyplot as plt
import numpy as np


# Constantes
gamma_i = 0.064
gamma_x = 0.004
lambda_i = 2.926e-5
lambda_x = 2.0996e-5
sigma_i = 7e-24
sigma_x = 2.65e-18
maj_sigma_f = 0.0984
maj_sigma_u = 0.1355
end_time = 200600  # Durée de la simulation (secondes)
sigma_b_min = 0.1  # Barres complètement relevées
sigma_b_max = 0.2  # Barres complètement enfoncées
flux_target = 1e15  # Flux cible (particules/s)
flux_initial = 1e10  # Flux initial (particules/s)
k = 100



# Simulation avec Euler explicite
def euler_coupled(derivatives, t0, y0, t_end, dt):
    t = t0
    y = np.array(y0)
    sigma_b = sigma_b_min  # Initialisation des barres
    times = [t]
    values = [y]
    sigma_b_values = [sigma_b]

    while t < t_end:
        dy_dt, sigma_b = derivatives(t, y, sigma_b)
        y = y + dt * np.array(dy_dt)
        t += dt
        times.append(t)
        values.append(y)
        sigma_b_values.append(sigma_b)

    return np.array(times), np.array(values), np.array(sigma_b_values)


# Système d'équations différentielles
def derivatives(t, y, sigma_b):
    I, grandX, psi_de_t = y


    #calcul des erreurs pour ajuster les bars de carbone pour stabiliser le flux a 1e15

    err = flux_target - psi_de_t 
    
    if (err/3e15+ 0.066) > 0.1:
        sigma_b = 0.1 
    else :
        sigma_b =  0.2 - err/3e15 - 0.066
    
    
    # Calcul de Σx
    maj_sigma_x = sigma_x * grandX

    # Équations différentielles
    dI_dt = gamma_i * maj_sigma_f * psi_de_t - lambda_i * I - sigma_i * I * psi_de_t
    dX_dt = gamma_x * maj_sigma_f * psi_de_t + lambda_i * I - lambda_x * grandX - sigma_x * grandX * psi_de_t
    dpsi_dt = (psi_de_t / 1000) * 3 * (maj_sigma_u - sigma_b - maj_sigma_x)

    return [dI_dt, dX_dt, dpsi_dt], sigma_b


# Simulation
def run_simulation():
    t0 = 0
    y0 = [0, 0, flux_initial]  # Conditions initiales : I = 0, X = 0, flux initial
    dt = 60  # step de calcul

    # Appel de la simulation
    times, results, sigma_b_values = euler_coupled(derivatives, t0, y0, end_time, dt)

    return times, results[:, 0], results[:, 1], results[:, 2], sigma_b_values


# Tracé des résultats
def plot_results(times, Is, Xs, Psis, sigma_bs):
    times_hours = times / 3600  # Conversion du temps en heures

    # Graphique 1 : Abondances d'iode et de xénon
    plt.figure(figsize=(12, 8))
    plt.subplot(3, 1, 1)
    plt.plot(times_hours, Is, label="Abondance d'iode (I)")
    plt.plot(times_hours, Xs, label="Abondance de xénon (X)")
    plt.xlabel("Temps (heures)")
    plt.ylabel("Abondance")
    plt.title("Évolution des abondances d'iode et de xénon")
    plt.legend()
    plt.grid()

    # Graphique 2 : Flux neutronique
    plt.subplot(3, 1, 2)
    plt.plot(times_hours, Psis, label="Flux neutronique (ψ)", color='orange')
    plt.axhline(flux_target, color='red', linestyle='--', label="Flux cible")
    plt.xlabel("Temps (heures)")
    plt.ylabel("Flux neutronique (ψ)")
    plt.title("Évolution du flux neutronique")
    plt.legend()
    plt.grid()

    # Graphique 3 : Position des barres de contrôle
    plt.subplot(3, 1, 3)
    plt.plot(times_hours, sigma_bs, label="Position des barres (Σb)", color='green')
    plt.xlabel("Temps (heures)")
    plt.ylabel("Position des barres (Σb)")
    plt.title("Évolution de la position des barres de contrôle")
    plt.legend()
    plt.grid()

    plt.tight_layout()
    plt.show()


# Exécution du script
if __name__ == "__main__":
    times, Is, Xs, Psis, sigma_bs = run_simulation()
    plot_results(times, Is, Xs, Psis, sigma_bs)