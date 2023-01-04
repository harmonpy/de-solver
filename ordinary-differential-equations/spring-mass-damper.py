"""
This is a spring-mass-damper system solver based on forward Euler method.
The system is a 2nd order ordinary differential (ode) equation. It is written as follow,
	\ddot{x} + 2\zeta \omega_0 v + \omega_0^2 x = 0
where,
x 			: the distance of mass displacement
\zeta 		: damping ratio -> \zeta = \frac{c}{2 \sqrt(km)}
\omega_0 	: natural frequency -> \omega_0 = \sqrt(k/m)

We know that the sistem includes spring constant k and damping constant c.
For simplicity we omit these constant in the code. Instead we will use \omega_0 and \zeta.

Alter the 2nd order ode to a system of 1st order ode,
	\dot{x} = v
	\dot{v} = -2 \zeta \omega_0 v + \omega_0^2 x
Also, we set the initial conditions as x(0) = 2 and v(0) = 0.  

We apply forward euler method to the ode of the system. We have
	x^{n+1} = x^n + \Delta t v^n
	v^{n+1} = v^n + \Delta t (-2 \zeta \omega_0 v^n - \omega_0^2 x^n)
with \zeta = 0.25, \omwga_0 = 2 \pi, \Delta t = 0.01, and the simulation runs at 0 >= t >= 10.
"""

import numpy as np
import matplotlib.pyplot as plt

# parameters
w = 2 * np.pi 	# \omega_0	
z = 0.25 		# \zeta
dt = 0.01 		# \Delta t
T = 10			# the end time of simulation

# resulting parameters
n = int(T/dt)
t = np.linspace(0, 10, n)
x = np.zeros(n)
v = np.zeros(n)

# initial condition
x[0] = 2
v[0] = 0

# forward Euler method
for i in range(n-1):
	x[i+1] = x[i] + dt * v[i]
	v[i+1] = v[i] + dt * (-2 * z * w * v[i] - w**2 * x[i])

# plot the displacement x
plt.plot(t, x)
plt.xlabel('time')
plt.ylabel('displacement')
plt.grid(True)
plt.show()