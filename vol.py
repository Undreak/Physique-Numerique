import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

################ Définition des fonctions

def vol(N,theta,tsim,m=29981,x0=0,z0=0,vx0=0,vz0=0,ax0=0,az0=0,Cx=0.5,S=np.pi*(0.89)**2,poussee=0,thetaF=0):
    ## Initialisation de tous les vecteurs.
    t = int(tsim*N)
    dt = 1/N
    x, z = np.zeros(t),np.zeros(t)
    x[0], z[0] = x0, z0
    vx,vz = np.zeros(t),np.zeros(t)
    vx[0],vz[0] = vx0, vz0
    ax,az = np.zeros(t),np.zeros(t)
    ax[0],az[0] = ax0,az0
    Ttab,rhotab = np.zeros(t),np.zeros(t)
    Ttab[0],rhotab[0] = atmosphere(z0)
    thetaI=theta
    Cpara,Spara = 1.75,0    ## Coefficient de frottement et surface des parachutes, nulle au départ car pas déployé.
    prop = q*ve*poussee              ## Force de poussée de la fusée
    
    for i in range(t-1):
        T, rho = atmosphere(z[i])  ## Calcul de la densité de l'air et de la température en fontion de l'altitude.
        if vz[i] < 0:
            Cx = coeff(vz[i],vx[i],T)    ## Calcul du coeff Cx lors de la descente
            Spara = para(z[i])             ## et de Spara à partir d'une certaine altitude.

        if z[i] < 0:            ## On arrete le calcul quand la fusée touhe le sol et on ne garde que les dernières données calculé
            x,z,m,vx,vz,ax,az,Ttab,rhotab = x[:i],z[:i],m,vx[:i],vz[:i],ax[:i],az[:i],Ttab[:i],rhotab[:i]
            break
        Ttab[i+1] = T
        rhotab[i+1] = rho
        if thetaF != 0 and i <= 24/dt:
            theta -= (thetaI - thetaF)*dt/24    ## Calcul de la variation d'angle, pi/2 => pi/4 en 24 secondes
        #ax[i+1] = ( prop - 0.5*rho*abs(vx[i])*vx[i]*(S*Cx + Spara*Cpara))*np.cos(theta)/m
        #az[i+1] = ( prop - 0.5*rho*abs(vz[i])*vz[i]*(S*Cx + Spara*Cpara))*np.sin(theta)/m - g

        if poussee == 1:    ## Calcul de la variation de masse
            m -= q*dt
        if m <= mcap:       ## Simple précaution pour éviter que le programme ne fasse n'importe quoi
            prop = 0
        ##Calcul de la vitesse par la méthode de Runge-Kutta d'ordre 4
        k1 = acceleration(z[i],vx[i],vz[i],theta,m,poussee,N)
        vx1,vz1 = vx[i]+k1[0]*dt/2 , vz[i]+k1[1]*dt/2
        k2 = acceleration(z[i],vx1,vz1,theta,m,poussee,N)
        vx2,vz2 = vx[i]+k2[0]*dt/2 , vz[i]+k2[1]*dt/2
        k3 = acceleration(z[i],vx2,vz2,theta,m,poussee,N)
        vx3,vz3 = vx[i]+k3[0]*dt , vz[i]+k3[1]*dt
        k4 = acceleration(z[i],vx3,vz3,theta,m,poussee,N)
        vx4,vz4 = vx[i]+dt*(k1[0]+2*k2[0]+2*k3[0]+k4[0])/6 , vz[i]+dt*(k1[1]+2*k2[1]+2*k3[1]+k4[1])/6

        vx[i+1] = vx4
        vz[i+1] = vz4

        x[i+1] = x[i] + vx[i]*dt
        z[i+1] = z[i] + vz[i]*dt
        
        ax[i+1] = ( prop - 0.5*rho*abs(vx[i+1])*vx[i+1]*(S*Cx + Spara*Cpara))*np.cos(theta)/m
        az[i+1] = ( prop - 0.5*rho*abs(vz[i+1])*vz[i+1]*(S*Cx + Spara*Cpara))*np.sin(theta)/m - g

    return x,z,m,vx,vz,ax,az,Ttab,rhotab

def acceleration(z,vx,vz,theta,m,poussee,N,S=np.pi*(0.89)**2,Cpara=1.75,Spara=0,Cx=0.5):
    T, rho = atmosphere(z)  ## Calcul de la densité de l'air et de la température en fontion de l'altitude.
    if vz < 0:
        Cx = coeff(vz,vx,T)    ## Calcul du coeff Cx lors de la descente
        Spara = para(z)             ## et de Spara à partir d'une certaine altitude.
    prop = q*ve
    if poussee == 1:    ## Calcul de la variation de masse
        m -= q/N
    if m <= mcap:       ## Simple précaution pour éviter que le programme ne fasse n'importe quoi
        prop = 0
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

g = 9.81
q = 166
ve = 2.11e3
mcap = 1400

############## Définition des variables

theta1 = np.pi/2
theta2 = np.pi/4
tsim1 = 16
tsim2 = 127.5
tsim3 = 700
N = 100

############## Calcul des trajectoires

## La fusée décolle verticalement
x1,z1,m1,vx1,vz1,ax1,az1,Ttab1,rhotab1 = vol(N,theta1,tsim1,poussee=1) 
## Puis prend progressivement un angle de 45 ◦ par rapport à l’horizontale ≈ 16 s après le lancement.
x2,z2,m2,vx2,vz2,ax2,az2,Ttab2,rhotab2 = vol(N,theta1,tsim2,x0=x1[-1],z0=z1[-1],m=m1,vx0=vx1[-1],vz0=vz1[-1],ax0=ax1[-1],az0=az1[-1],poussee=1,thetaF=theta2)
## Largage de la tour d’éjection d’urgence à la fin de la poussée de la fusée
x3,z3,m3,vx3,vz3,ax3,az3,Ttab3,rhotab3 = vol(N,theta2,4,x0=x2[-1],z0=z2[-1],m=m2-580,vx0=vx2[-1],vz0=vz2[-1],ax0=ax2[-1],az0=az2[-1])
## Séparation de la capsule Mercury 4 secondes après la fin de la poussée
x4,z4,m4,vx4,vz4,ax4,az4,Ttab4,rhotab4 = vol(N,np.pi*34/180,tsim3,x0=x3[-1],z0=z3[-1],m=1400,vx0=vx3[-1],vz0=vz3[-1],ax0=ax3[-1],az0=az3[-1])

t1 = np.linspace(0,len(z1)/N,len(z1))
t12 = np.linspace(t1[-1],40,2400)
t2 = np.linspace(t12[-1],143.5,10350)
t3 = np.linspace(t2[-1],t2[-1]+4,len(z3))
t4 = np.linspace(t3[-1],len(z1)/N+len(z2)/N+len(z3)/N+len(z4)/N,len(z4))

para1 = int(np.argwhere(z4<=6700)[0])
para2 = int(np.argwhere(z4<=3000)[0])
############## Tracé des graphiques

#plt.plot(t1,az1) and plt.plot(t2,az2) and plt.plot(t3,az3) and plt.plot(t4,az4)
#plt.plot(t1,ax1) and plt.plot(t2,ax2) and plt.plot(t3,ax3) and plt.plot(t4,ax4)
plt.plot(t1,(az1**2 + ax1**2)**0.5) and plt.plot(t12,(az2[0:2400]**2 + ax2[0:2400]**2)**0.5) and plt.plot(t2,(az2[2400:]**2 + ax2[2400:]**2)**0.5) and plt.plot(t3,(az3**2 + ax3**2)**0.5) and plt.plot(t4,(az4**2 + ax4**2)**0.5)
plt.title('Accélération en fonction du temps')
plt.xlabel('t (s)', color='tab:blue')
plt.ylabel('a (m/s^2)', color='tab:red')
plt.tick_params(axis='y', labelcolor='tab:red')
plt.tick_params(axis='x', labelcolor='tab:blue')
plt.savefig('Accélération en fonction du temps.png',dpi=300)

fig, ax0 = plt.subplots()

color = 'tab:red'
ax0.set_xlabel('x (km)', color='tab:blue')
ax0.set_ylabel('z (km)', color=color)
ax0.set_title('Vol d\'Alan Shepard')
ax0.plot(x1/1000,z1/1000, label='Décollage')
ax0.plot(x2[0:2400]/1000,z2[0:2400]/1000, label='Rotation de la fusée')
ax0.plot(x2[2400:]/1000,z2[2400:]/1000, label='Poussée')
ax0.plot(x3/1000,z3/1000, label='Largage de la tour de sauvetage')
ax0.plot(x4/1000,z4/1000, label='Chute libre')
ax0.plot(x4[para1]/1000,z4[para1]/1000,'.', label='Ouverture du parachute pilote')
ax0.plot(x4[para2]/1000,z4[para2]/1000,'.', label='Ouverture du parachute principal')

print('Temps de vol :',round((len(z1) + len(z2)+ len(z3) + len(z4))/(N*60),2))
print('Distance parcourue :', int(max(x4))/1000)
print('Altitude maximale :', int(max(z4))/1000)

ax0.tick_params(axis='y', labelcolor=color)
ax0.tick_params(axis='x', labelcolor='tab:blue')

leg = ax0.legend(loc="best",ncol=1, shadow=True, title="Phases de vol", fancybox=True)
leg.get_title().set_color("red")

plt.savefig('Vol d\'Alan Shepard.png',dpi=300)

plt.show()
'''
fig, ax1 = plt.subplots()

color = 'tab:red'
ax1.set_xlabel('time (s)')
ax1.set_ylabel('z', color=color)
ax1.plot(z4, color=color)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()

color = 'tab:blue'
ax2.set_ylabel('g', color=color)
ax2.plot(az4, color=color)
ax2.tick_params(axis='y', labelcolor=color)

fig.tight_layout()
plt.savefig('Accélération et postion au cours du temps',dpi=300)
plt.show()
'''