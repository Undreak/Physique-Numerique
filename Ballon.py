import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

############ Définition des fonctions

def tir_num(theta,v0,tsim,g = 9.81,x0=16,z0=0):
    N = 1000
    t = np.linspace(0,tsim,N)
    x,z = np.zeros(N), np.zeros(N)
    x[0], z[0] = x0, z0
    dt = tsim/N
    dxdt = v0*np.cos(theta)
    for i in range(N-1):
        x[i+1] = x[i] + dxdt*dt
        z[i+1] = z[i] + (- g*t[i] + v0*np.sin(theta))*dt
        if z[i+1] < 0 :
            break

    return x,z

def tir_num_f(theta,v0,tsim,m=0.450,g=9.81,x0=16,z0=0,Cx=0.45,S=np.pi*(0.11)**2,rho=1.225):
    N = 1000
    x, z = np.zeros(N),np.zeros(N)
    x[0], z[0] = x0, z0
    dxdt,dzdt = np.zeros(N),np.zeros(N)
    dxdt[0],dzdt[0] = v0*np.cos(theta),v0*np.sin(theta)
    dt = tsim/N
    for i in range(N-1):
        dxdt[i+1] = dxdt[i] -0.5*rho*Cx*S*abs(dxdt[i])*dxdt[i]*dt
        dzdt[i+1] = dzdt[i] + ( -0.5*rho*Cx*S*abs(dzdt[i])*dzdt[i] - m*g)*dt/m 
        x[i+1] = x[i] + dxdt[i+1]*dt
        z[i+1] = z[i] + dzdt[i+1]*dt
        if z[i+1] < 0:
            break
    return x,z

############ Définition des constantes

theta1 = np.pi/6
theta3 = np.pi/4
theta5 = np.pi/2.5
v0 = 40
tsim = 10
x0 = 16

############# Calculs des trajectoires
x, z = tir_num(theta3,v0,tsim)

xf, zf = tir_num_f(theta3,v0,tsim)

############# Recherche de l'angle maximal
thetaTAB = np.zeros(10)
xMAX = np.zeros(10)
for i in range(3,13):
    thetaTAB[i-3] = np.pi/i
    xt,zt = tir_num_f(thetaTAB[i-3],v0,tsim)
    xMAX[i-3] = max(xt)

print('maximum de x avec frottement => theta = pi/6')
print('maximum de x sans frottement => theta = pi/4')
print('xmax = ', max(x))
print('xfmax = ', max(xf))
############# Tracé des graphiques

fig, ax0 = plt.subplots()
ax0.set_ylim(0,55)
p1, = ax0.plot(x,z, label='sans frottements')
p2, = ax0.plot(xf,zf,label='avec frottements',color='red', marker='.', linestyle='dashed',linewidth=1, markersize=0.1)

leg = ax0.legend(loc="upper right",ncol=2, shadow=True, title="Légende", fancybox=True)
leg.get_title().set_color("red")

ax0.set(xlabel='distance (metre)', ylabel='hauteur (metre)', title='comparaison avec et sans frottements')

plt.savefig('frottements.png',dpi=300)
plt.show()
