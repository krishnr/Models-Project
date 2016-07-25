from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from math import tan, cos, sin, pi
from scipy.integrate import odeint

##############
## y0 = yk
## y1 = theta
## y2 = px
## y3 = py
##############

def model(y, t):
    yk, theta, vx, vy = y

    # constants
    m = 0.68 # Mass of pendulum
    k = 300 # Stiffness of Spring
    Y = 0
    R = 0.1 # Torsional Resistance
    L = 0.61  # Length of Pendulum
    g = 9.81 # Gravitional acceleration

    # in between terms
    disp = (yk - Y)

    d_yk = (tan(theta) * vx) + vy
    d_theta = vx / L * cos(theta)
    d_vx = ((-tan(theta)) * k / m * disp) - ((R*vx) / (m * (L * cos(theta))**2))
    d_vy = (-k / m * disp) + g

    return [d_yk, d_theta, d_vx, d_vy]

time = np.linspace(0.0, 20.0, 100)
yinit = [0, pi/4, 0, 0]
y = odeint(model, yinit, time)

plt.plot(time, y[:, 0], label="yk")
plt.xlabel('t')
plt.ylabel('Motion of Spring')
plt.legend()
plt.show()

plt.plot(time, y[:, 1], label="theta")
plt.xlabel('t')
plt.ylabel('Angle of rotation')
plt.legend()
plt.show()

plt.plot(time, y[:, 2], label="px")
plt.xlabel('t')
plt.ylabel('Displacement of Pendulum in Y')
plt.legend()
plt.show()

plt.plot(time, y[:, 3], label="py")
plt.xlabel('t')
plt.ylabel('Displacement of Pendulum in Y')
plt.legend()
plt.show()
