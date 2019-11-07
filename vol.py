import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

################ Définition des fonctions

def vol(theta,v0,tsim,m=0.01,g=9.81,x0=0,z0=0,Cx=0.5,S=np.pi*(0.1)**2):
    N = 10000
    t = np.linspace(0,tsim,N)
    x, z = np.zeros(N),np.zeros(N)
    x[0], z[0] = x0, z0
    dxdt,dzdt = np.zeros(N),np.zeros(N)
    dxdt[0],dzdt[0] = v0*np.cos(theta),v0*np.sin(theta)
    dt = tsim/N

    for i in range(N-1):
        if z[i] < 0:
            break
        elif z[i] <= 11000:
            T = 288.15 - 6.5e-3 *(z[i])
            rho = 1.2*(288.15/T)**(1 + 34.16e-3/(-6.5e-3))
        elif 11000 < z[i] <= 20000:
            T = 216.65
            rho = 0.35654*np.exp(-(36.16e-3 *(z[i] - 11000)/T))
        elif 20000 < z[i] <= 32000:
            T = 288.15 + 1e-3*(z[i] - 20000)
            rho = 0.086261*(216.65/T)**(1 + 34.16e-3/(1e-3))
        elif 32000 < z[i] <= 47000:
            T = 288.15  + 2.8e-3 *(z[i] - 32000)
            rho = 0.012961*(228.65/T)**(1 + 34.16e-3/(2.8e-3))
        elif z[i] > 47000:
            T = 270.65
            rho = 0.0013993*np.exp(-((36.16e-3) *(z[i] - 47000)/T))

        dxdt[i+1] = dxdt[i] -0.5*rho*Cx*S*abs(dxdt[i])*dxdt[i]*dt
        dzdt[i+1] = dzdt[i] + ( -0.5*rho*Cx*S*abs(dzdt[i])*dzdt[i] - m*g*t[i])*dt
        x[i+1] = x[i] + dxdt[i+1]*dt
        z[i+1] = z[i] + dzdt[i+1]*dt
        print(i)
        print(dzdt[i])
    return x,z

############## Définition des variables

theta = np.pi/2
v0 = 2000
tsim = 100
t = np.linspace(0,tsim,10000)
############## Calcul des trajectoires

x,z = vol(theta,v0,tsim)

############## Tracé des graphiques

plt.plot(t,z)
plt.show()