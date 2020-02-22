##################################################
#
# Objectifs n°1: optimiser le programme pour le rendre plus robuste
#
################ Importation des bibliothèques 

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

################ Définition des fonctions

def vol(data,N,theta,tsim,m=30e3,poussee=0,rot=False,thetaF=0,Trot=0):
    ## Initialisation de tous les vecteurs
    t = int(tsim*N)
    dt = 1/N
    x, z = np.zeros(t),np.zeros(t)
    x[0], z[0] = data[0], data[1]
    vx,vz = np.zeros(t),np.zeros(t)
    vx[0],vz[0] = data[2], data[3]
    ax,az = np.zeros(t),np.zeros(t)
    ax[0],az[0] = data[4],data[5]
    thetaI=theta
    
    for i in range(t-1):
        if z[i] < 0:            ## On arrete le calcul quand la fusée touhe le sol et on ne garde que les dernières données calculé
            x,z,m,vx,vz,ax,az = x[:i],z[:i],m,vx[:i],vz[:i],ax[:i],az[:i]
            break
        if rot == True and i <= Trot/dt:
            theta -= (thetaI - thetaF)*dt/Trot    ## Calcul de la variation d'angle, pi/2 => pi/4 en 24 secondes

        if poussee == 1:    ## Calcul de la variation de masse
            m -= q*dt

        ##Calcul de la vitesse par la méthode de Runge-Kutta d'ordre 4
        k1 = acceleration(z[i],vx[i],vz[i],theta,m,poussee)
        vx1,vz1 = vx[i]+k1[0]*dt/2 , vz[i]+k1[1]*dt/2
        k2 = acceleration(z[i],vx1,vz1,theta,m,poussee)
        vx2,vz2 = vx[i]+k2[0]*dt/2 , vz[i]+k2[1]*dt/2
        k3 = acceleration(z[i],vx2,vz2,theta,m,poussee)
        vx3,vz3 = vx[i]+k3[0]*dt , vz[i]+k3[1]*dt
        k4 = acceleration(z[i],vx3,vz3,theta,m,poussee)
        vx4,vz4 = vx[i]+dt*(k1[0]+2*k2[0]+2*k3[0]+k4[0])/6 , vz[i]+dt*(k1[1]+2*k2[1]+2*k3[1]+k4[1])/6

        vx[i+1] = vx4
        vz[i+1] = vz4

        x[i+1] = x[i] + vx[i]*dt
        z[i+1] = z[i] + vz[i]*dt
        
        ax[i+1],az[i+1] = acceleration(z[i+1],vx[i+1],vz[i+1],theta,m,poussee)

    return np.array((x,z,m,vx,vz,ax,az))

def acceleration(z,vx,vz,theta,m,poussee,S=np.pi*(0.89)**2,Cpara=1.75,Spara=0,Cx=0.5):
    T, rho = atmosphere(z)  ## Calcul de la densité de l'air et de la température en fontion de l'altitude.
    if vz < 0:
        Cx = coeff(vz,vx,T)    ## Calcul du coeff Cx lors de la descente
        Spara = para(z)             ## et de Spara à partir d'une certaine altitude.
    prop = q*ve*poussee
    ax = ( prop - 0.5*rho*abs(vx)*vx*(S*Cx + Spara*Cpara))*np.cos(theta)/m
    az = ( prop - 0.5*rho*abs(vz)*vz*(S*Cx + Spara*Cpara))*np.sin(theta)/m - g

    return ax, az
    
def atmosphere(z):
    B = 34.16e-3
    if z <= 11000:
        lbd = -6.5e-3
        T = 288.15 + lbd*(z)
        rho = 1.2*(288.15/T)**(1 + B/lbd)

    elif 11000 < z <= 20000:
        T = 216.65
        rho = 0.35654*np.exp(-(B *(z - 11000)/T))

    elif 20000 < z <= 32000:
        lbd = 1e-3
        T = 216.65 + lbd*(z - 20000)
        rho = 0.086261*(216.65/T)**(1 + B/lbd)

    elif 32000 < z <= 47000:
        lbd = 2.8e-3
        T = 228.65  + lbd*(z - 32000)
        rho = 0.012961*(228.65/T)**(1 + B/lbd)

    elif z > 47000:
        T = 270.65
        rho = 0.0013993*np.exp(-(B *(z - 47000)/T))

    return T,rho

def coeff(vz,vx,T):
    M = (vz**2 + vx**2)**0.5/(1.4*287*T)**0.5
    if M < 1:
        Cx = 0.9
    else :
        Cx = 1.5
    return Cx

def para(z):
    Spara = 0
    if z <= 6700:
        Spara = np.pi*(1)**2
    if z <= 3000:
        Spara = np.pi*(2)**2
    if z <= 2980:
        Spara = np.pi*(3)**2
    if z <= 2960:
        Spara = np.pi*(4)**2
    if z <= 2940:
        Spara = np.pi*(5)**2
    if z <= 2920:
        Spara = np.pi*(6)**2
    if z <= 2900:
        Spara = np.pi*(7)**2
    if z <= 2900:
        Spara = np.pi*(8)**2
    if z <= 2900:
        Spara = np.pi*(9)**2
    return Spara

############## Définition des constantes

N = 100
g = 9.81
q = 166
ve = 2.11e3
mcap = 1400

################### Programme principal

dataZero = np.zeros(6)
data1 = vol(dataZero,N,np.pi/3,50,poussee=1)
data2 = vol(data1[-1],N,np.pi/3,100)

################### Graphiques

plt.plot(data1[0],data1[1])
plt.plot(data2[0],data2[1])
plt.show()