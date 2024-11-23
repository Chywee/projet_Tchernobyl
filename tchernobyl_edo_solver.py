

from math import sqrt

import matplotlib
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
#definition des constente
gamma_i = 0.064
gamma_x = 0.004
lambda_i = 2.926*10^-5
lambda_x = 2.0996*10^-5
sigma_i = 7*10^-24
sigma_x = 2.65*10¨-18
sigma_f = 0.0984
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



