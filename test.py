import matplotlib.pyplot as plt
from scipy.integrate import odeint
import numpy as np

def pend(y, t, H0, Om_M, Om_R, Om_Lambda):
    a, aprime = y
    dydt = [aprime, -(a/2)*H0**2*(Om_M/a**3+2*Om_R/a**4-2*Om_Lambda)]
    return dydt

# We assume the constants are:
H0=1000/3.086e22*(3600*24*365*1e09)*67
Om_M = 0.3
Om_R = 0
Om_Lambda = 0.7
# Initial conditions
y0 = [1, H0]

t = np.linspace(0, 30, 1000)
t2 = np.linspace(0, -14.15, 1000)

sol = odeint(pend, y0, t, args=(H0, Om_M, Om_R, Om_Lambda))
sol2 = odeint(pend, y0, t2, args=(H0, Om_M, Om_R, Om_Lambda))

#plt.plot(t, sol[:, 0], 'b', t2, sol2[:,0], 'r', label='a(t)')
plt.plot(t, sol[:, 1]/sol[:,0], 'b')
plt.plot(t2,sol2[:,1]/sol2[:,0])
# plt.plot(t, sol[:, 1], 'g', label='omega(t)')
plt.legend(loc='best')
plt.xlabel('t')
plt.ylim(0,30)
plt.grid()
plt.show()