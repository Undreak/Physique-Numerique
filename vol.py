import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

################ Définition des fonctions

def vol(N,theta,tsim,m=30e3,g=9.81,x0=0,z0=0,vx0=0,vz0=0,ax0=0,az0=0,Cx=0.5,S=np.pi*(0.8)**2,q=166,ve=2.11e3,mcap=1400,poussee=0,thetaF=0):
    ##Initialisation de tous les vecteurs.
    t = int(tsim*N)
    x, z = np.zeros(t),np.zeros(t)
    x[0], z[0] = x0, z0
    vx,vz = np.zeros(t),np.zeros(t)
    vx[0],vz[0] = vx0, vz0
    dt = tsim/t
    ax,az = np.zeros(t),np.zeros(t)
    ax[0],az[0] = ax0,az0
    Ttab,rhotab = np.zeros(t),np.zeros(t)
    Ttab[0],rhotab[0] = atmosphere(z0)
    thetaI=theta
    Cpara,Spara = 1.75,0    ##Coefficient de frottement et surface des parachutes, nulle au départ car pas déployé.
    prop = q*ve*poussee              ##Force de poussée de la fusée
    
    for i in range(t-1):
        T, rho = atmosphere(z[i])  ##Calcul de la densité de l'air et de la température en fontion de l'altitude.
        if vz[i] < 0:
            Cx = coeff(vz[i],T)    ##Calcul du coeff Cx lors de la redescente
            Spara = para(z[i])             ##et de Spara à partir d'une certaine altitude.
        if z[i] < 0:
            x,z,m,vx,vz,ax,az,Ttab,rhotab = x[:i],z[:i],m,vx[:i],vz[:i],ax[:i],az[:i],Ttab[:i],rhotab[:i]
            break
        Ttab[i+1] = T
        rhotab[i+1] = rho
        if thetaF != 0 and i <= 24/dt:
            theta -= (thetaI - thetaF)*dt/24    ##Calcul de la variation d'angle, pi/2 => pi/4 en 24 secondes

        ax[i+1] = ( prop*np.cos(theta) - 0.5*rho*(S*Cx + Spara*Cpara)*S*abs(vx[i])*vx[i])/m
        az[i+1] = ( prop*np.sin(theta) - 0.5*rho*abs(vz[i])*vz[i]*(S*Cx + Spara*Cpara) - m*g)/m

        if poussee == 1:
            m -= q*dt
        if m <= mcap:
            prop = 0

        vx[i+1] = vx[i] + (ax[i+1] + ax[i])*dt/2
        vz[i+1] = vz[i] + (az[i+1] + az[i])*dt/2

        x[i+1] = x[i] + (vx[i+1] + vx[i])*dt/2
        z[i+1] = z[i] + (vz[i+1] + vz[i])*dt/2
    return x,z,m,vx,vz,ax,az,Ttab,rhotab

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
        T = 288.15 + lbd*(z - 20000)
        rho = 0.086261*(216.65/T)**(1 + B/lbd)

    elif 32000 < z <= 47000:
        lbd = 2.8e-3
        T = 288.15  + lbd*(z - 32000)
        rho = 0.012961*(228.65/T)**(1 + B/lbd)

    elif z > 47000:
        T = 270.65
        rho = 0.0013993*np.exp(-(B *(z - 47000)/T))

    return T,rho

def coeff(vz,T):
    M = vz/(1.4*287*T)**0.5
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
    if z <= 2990:
        Spara = np.pi*(3)**2
    if z <= 2980:
        Spara = np.pi*(4)**2
    if z <= 2970:
        Spara = np.pi*(5)**2
    if z <= 2960:
        Spara = np.pi*(6)**2
    if z <= 2950:
        Spara = np.pi*(7)**2
    if z <= 2940:
        Spara = np.pi*(8)**2
    return Spara

############## Définition des variables

theta1 = np.pi/2
theta2 = np.pi/4
tsim1 = 16
tsim2 = 127.5
tsim3 = 700
N = 100

############## Calcul des trajectoires

x1,z1,m1,vx1,vz1,ax1,az1,Ttab1,rhotab1 = vol(N,theta1,tsim1,poussee=1) ##La fusée décolle verticalement
x2,z2,m2,vx2,vz2,ax2,az2,Ttab2,rhotab2 = vol(N,theta1,tsim2,x0=x1[-1],z0=z1[-1],m=m1,poussee=1,vx0=vx1[-1],vz0=vz1[-1],ax0=ax1[-1],az0=az1[-1],thetaF=theta2)   ##Puis prend progressivement un angle de 45 ◦ par rapport à l’horizontale ≈ 16 s après le lancement.
x3,z3,m3,vx3,vz3,ax3,az3,Ttab3,rhotab3 = vol(N,0,4,x0=x2[-1],z0=z2[-1],m=m2-580,vx0=vx2[-1],vz0=vz2[-1],ax0=ax2[-1],az0=az2[-1])    ##Largage de la tour d’éjection d’urgence à la fin de la poussée de la fusée
x4,z4,m4,vx4,vz4,az4,az4,Ttab4,rhotab4 = vol(N,0,tsim3,x0=x3[-1],z0=z3[-1],m=1400,vx0=vx3[-1],vz0=vz3[-1],ax0=ax3[-1],az0=az3[-1])    ##Séparation de la capsule Mercury 4 secondes après la fin de la poussée

############## Tracé des graphiques

fig, ax0 = plt.subplots()

color = 'tab:red'
ax0.set_xlabel('x (km)', color='tab:blue')
ax0.set_ylabel('z (km)', color=color)
ax0.set_title('Vol d\'Alan Shepard')
ax0.plot(x1/1000,z1/1000)
ax0.plot(x2/1000,z2/1000)
ax0.plot(x3/1000,z3/1000, label=int(max(x4))/1000)
ax0.plot(x4/1000,z4/1000, label=int(max(z4))/1000)

ax0.tick_params(axis='y', labelcolor=color)
ax0.tick_params(axis='x', labelcolor='tab:blue')

leg = ax0.legend(loc="center",ncol=2, shadow=True, title="x max et z max (km)", fancybox=True)
leg.get_title().set_color("red")

plt.savefig('Vol d\'Alan Shepard.png',dpi=300)

plt.show()
'''
fig, ax1 = plt.subplots()

color = 'tab:red'
ax1.set_xlabel('time (s)')
ax1.set_ylabel('z', color=color)
ax1.plot(np.linspace(147.5,tsim3+147.5,N),z4, color=color)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()

color = 'tab:blue'
ax2.set_ylabel('g', color=color)
ax2.plot(np.linspace(147.5,tsim3+147.5,N),az4/10, color=color)
ax2.tick_params(axis='y', labelcolor=color)

fig.tight_layout()
plt.savefig('Accélération et postion au cours du temps',dpi=300)
plt.show()
'''