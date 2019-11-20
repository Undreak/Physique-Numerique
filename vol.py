import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

################ Définition des fonctions

def vol(theta,v0,tsim,m=30000,g=9.81,x0=0,z0=0,Cx=0.5,S=np.pi*(0.5)**2,q=166,ve=2.11e3,mcap=1400):
    x, z = np.zeros(N),np.zeros(N)
    x[0], z[0] = x0, z0
    dxdt,dzdt = np.zeros(N),np.zeros(N)
    dxdt[0],dzdt[0] = v0*np.cos(theta),v0*np.sin(theta)
    dt = tsim/N
    dxdt2,dzdt2 = np.zeros(N),np.zeros(N)
    dxdt2[0],dzdt2[0] = 0,0
    Ttab,rhotab = np.zeros(N),np.zeros(N)
    Ttab[0],rhotab[0] = 288.15,1.2*(288.15/288.15)**(1 + 34.16e-3/(-6.5e-3))
    Cpara,Spara = 1.75,0
    prop=q*ve
    for i in range(N-1):
        T, rho = atmostphere(z[i])
        if dzdt[i] < 0:
            Cx,Spara = coeff(dzdt[i],T,z[i])
        if z[i] < 0:
            break
        Ttab[i+1] = T
        rhotab[i+1] = rho

        dxdt2[i+1] = ( prop - 0.5*rho*Cx*S*abs(dxdt[i])*dxdt[i])/m
        dzdt2[i+1] = ( prop - 0.5*rho*abs(dzdt[i])*dzdt[i]*(S*Cx + Spara*Cpara) - m*g)/m

        if m >= mcap:
            m -= q*dt
        if m <= mcap:
            prop = 0
        if dt*i >= 143.5:
            prop = 0
        dxdt[i+1] = dxdt[i] + dxdt2[i+1]*dt
        dzdt[i+1] = dzdt[i] + dzdt2[i+1]*dt

        x[i+1] = x[i] + dxdt[i+1]*dt
        z[i+1] = z[i] + dzdt[i+1]*dt
    return x,z,dxdt,dzdt,dxdt2,dzdt2,Ttab,rhotab

def atmostphere(z):
    if z <= 11000:
        T = 288.15 - 6.5e-3*(z)
        rho = 1.2*(288.15/T)**(1 + 34.16e-3/(-6.5e-3))
    elif 11000 < z <= 20000:
        T = 216.65
        rho = 0.35654*np.exp(-(36.16e-3 *(z - 11000)/T))
    elif 20000 < z <= 32000:
        T = 288.15 + 1e-3*(z - 20000)
        rho = 0.086261*(216.65/T)**(1 + 34.16e-3/(1e-3))
    elif 32000 < z <= 47000:
        T = 288.15  + 2.8e-3 *(z - 32000)
        rho = 0.012961*(228.65/T)**(1 + 34.16e-3/(2.8e-3))
    elif z > 47000:
        T = 270.65
        rho = 0.0013993*np.exp(-((36.16e-3) *(z - 47000)/T))
    return T,rho

def coeff(dzdt,T,z):
    M = dzdt/(1.4*287*T)**2
    S = 0
    if M < 1:
        Cx = 0.9
    else :
        Cx = 1.5
    if z <= 6700:
        S = np.pi*(1)**2
        if z <= 3000:
            S = np.pi*(8)**2
    return Cx,S

############## Définition des variables

theta = np.pi/4
v0 = 0
tsim = 2000
N = 10000
t = np.linspace(0,tsim,N)

############## Calcul des trajectoires

x,z,dxdt,dzdt,dxdt2,dzdt2,Ttab,rhotab = vol(theta,v0,tsim)

if abs(dzdt[np.argmin(z)]) > 20:
    print('CRASH DE LA CAPSULE')
if max(abs(dzdt2))/10 > 50:
    print('ACCÉLÉRATION TROP FORTE')
############## Tracé des graphiques

plt.plot(t,z)
plt.show()
