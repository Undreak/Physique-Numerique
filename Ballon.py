import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def tir(theta, v0, tsim, g = 9.81, Cx=0.45,S=np.pi*(11)**2,m=0.450,x0=16,z0=0):
    t = np.linspace(0,tsim,1000)
    x = v0*t*np.cos(theta) + x0
    z = -(1/2)*g*t**2 + v0*t*np.sin(theta)
    n = 0
    while True:
        if i > len(z)-1:
            break
        if z[n] < 0:
            break
        n += 1
    return x[0:n],z[0:n]

def tir_num(theta,v0,tsim,g = 9.81,x0=16,z0=0):
    N = 577
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

theta1 = np.pi/6
theta3 = np.pi/4
theta5 = np.pi/2.5
v0 = 40
tsim = 10
x0 = 16

x_num, z_num = tir_num(theta3,v0,tsim)
plt.plot(x_num,z_num)

x1,z1 = tir(theta1,v0,tsim)
x3,z3 = tir(theta3,v0,tsim)
x5,z5 = tir(theta5,v0,tsim)

fig, ax0 = plt.subplots()
ax0.set_ylim(0,55)
p1, = ax0.plot(x3,z3, label='résultats analytique')
p2, = ax0.plot(x_num,z_num,label='résultats numérique',color='red', marker='.', linestyle='dashed',linewidth=1, markersize=0.1)

leg = ax0.legend(loc="upper right",ncol=2, shadow=True, title="Légende", fancybox=True)
leg.get_title().set_color("red")

ax0.set(xlabel='distance (metre)', ylabel='hauteur (metre)', title='comparaison entre les résultats analytique et numérique')

plt.savefig('analytique_vs_numerique.png',dpi=300)
plt.show()

fig, ax = plt.subplots()
ax.set_ylim(0,50)

l1, = ax.plot(x1,z1,label=int(x1[-1]))
l3, = ax.plot(x3,z3,label=int(x3[-1]))
l5, = ax.plot(x5,z5,label=int(x5[-1]))

tmax1 = round(x1[-1]/(v0*np.cos(theta1)),1)
tmax3 = round(x3[-1]/(v0*np.cos(theta1)),1)
tmax5 = round(x5[-1]/(v0*np.cos(theta1)),1)

leg1 = ax.legend(loc="upper right",ncol=2, shadow=True, title="Distance maximum (m)", fancybox=True)
leg1.get_title().set_color("red")
leg2 = ax.legend((l1,l3,l5),(tmax1,tmax3,tmax5),loc="upper left",ncol=2, shadow=True, title="Durée du vol (s)", fancybox=True)
leg2.get_title().set_color("green")
leg3 = ax.legend((l1,l3,l5),(round(theta1,1),round(theta3,1),round(theta5,1)),loc="lower center",ncol=2, shadow=True, title="Angle (radians)", fancybox=True)
leg3.get_title().set_color("blue")
ax.add_artist(leg1)
ax.add_artist(leg2)

ax.set(xlabel='distance (metre)', ylabel='hauteur (metre)', title='tir de ballon')

plt.savefig('ballon_analytique_angle.png',dpi=500)
plt.show()
