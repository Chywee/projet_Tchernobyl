
#beepboop of mozerfroking maphys hard coded baby B)


import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
#definition des constente
gamma_i = 0.064
gamma_x = 0.004
lambda_i = 2.926*10 ** (-5)
lambda_x = 2.0996*10 ** (-5)
sigma_i = 7*10 ** (-24)
sigma_x = 2.65*10 ** (-18)
maj_sigma_f = 0.0984
maj_sigma_u = 0.1355
psi_de_t = 3*10 ** (13)
end_time = 178200
def maj_sigma_x (sigma_x,grandX):
    return sigma_x*grandX
#euler
def euler(derivatives, x, y, x_stop, h):
    X = [x]
    Y = [y]
    while x < x_stop:
        y = y + h * derivatives(x, y)
        x = x + h
        X.append(x)
        Y.append(y)
    return np.array(X), np.array(Y)
#edo non résolue
def dérivé_I(I):
    return np.array(gamma_i * maj_sigma_f * psi_de_t - lambda_i * I - sigma_i * I * psi_de_t)

def dérivé_grandX(X, I):
    return np.array(gamma_x * maj_sigma_f * psi_de_t - lambda_i * I - sigma_x * grandX * psi_de_t - lambda_x * grandX)
#ou est le t pour les deux edo non résolue ?? 
def run_simuI():
    #Run simulation and output times and I
    start_time = 0
    I0 = 0
    step = 3600

    times, I = euler(dérivé_I, start_time, I0 , end_time, step)
    Is = I#a modifier


    return times, Is

def run_simugrandX():
    #Run simulation and output times and X
    start_time = 0
    X0 = 0
    step = 3600

    times, grandX = euler(dérivé_grandX, start_time, X0 , end_time, step)
    Xs = grandX#a modifier


    return times, Xs

def plot_simuI(times, Is):
    """Generate plot from the simulation"""
    plt.plot(times, Is[:, 0], label="I")
    
    plt.legend()

    plt.show()

def plot_simugrandX(times, Xs):
    """Generate plot from the simulation"""
    plt.plot(times, Xs[:, 0], label="X")
    
    plt.legend()

    plt.show()



if __name__ == "__main__":
    outI = run_simuI()
    outX = run_simugrandX()
    plot_simugrandX(*outX)
    plot_simuI(*outI)
