

from math import sqrt

import matplotlib
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
#definition des constente
gamma_i = 0.064
gamma_x = 0.004
lambda_i = 2.926*10^-5
lambda_x = 2.0996*10^-5
sigma_i = 7*10^-24
sigma_x = 2.65*10^-18
maj_sigma_f = 0.0984
maj_sigma_u = 0.1355
psi_de_t = float(input("quel est la valeur du flux?  "))
end_time = int(input("après combien de temps"))
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
def dérivé_I(I, t):
    return np.array(gamma_i * maj_sigma_f * psi_de_t - lambda_i * I - sigma_i * I * psi_de_t)

def dérivé_X(X, I, t):
    return np.array(gamma_x * maj_sigma_f * psi_de_t - lambda_i * I - sigma_x * X * psi_de_t - lambda_x * X)

def run_simuI():
    #"""Run simulation and output times and I"""
    start_time = 0
    I0 = 
    step = 3600

    times, I = euler(dérivé_I, start_time, I0 , end_time, step)
    Is = I[:, [0, 5, 10]]


    return times, Is

def run_simuX():
    #"""Run simulation and output times and I"""
    start_time = 0
    X0 = 
    step = 3600

    times, X = euler(dérivé_X, start_time, X0 , end_time, step)
    Xs = X[:, [0, 5, 10]]


    return times, Xs

def plot_simuI(times, Is):
    """Generate plot from the simulation"""
    plt.plot(times, Is[:, 0], label="I")
    
    plt.legend()

    plt.show()

def plot_simuX(times, Xs):
    """Generate plot from the simulation"""
    plt.plot(times, Xs[:, 0], label="X")
    
    plt.legend()

    plt.show()



if __name__ == "__main__":
    outI = run_simuI()
    outX = run_simuX()
    plot_simuX(*outX)
    plot_simuI(*outI)
    
