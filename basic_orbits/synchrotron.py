#!/usr/bin/env python
# coding: utf-8

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

B0 = 1.0
E = 0.00

def gamma(ux,uy,uz):
    g = np.sqrt(1 + ux**2 + uy**2 + uz**2)
    return g

def gamma_cyl(ur,uphi,uz):
    g_cyl = np.sqrt(1 + ur**2 + uphi**2 + uz**2)
    return g_cyl


def syn(state, t, B0):
    
    x  = state[0]
    ux = state[1]
    y  = state[2]
    uy = state[3]
    z  = state[4]
    uz = state[5]
    
    duxdt = uy/gamma(ux, uy, uz)
    duydt = -ux/gamma(ux, uy, uz)
    duzdt = 0.0
    
    dxdt = ux/gamma(ux,uy,uz)
    dydt = uy/gamma(ux,uy,uz)
    dzdt = uz/gamma(ux,uy,uz)
    
    derivs = [dxdt, duxdt, dydt, duydt, dzdt, duzdt]
    
    return derivs

def syn_cyl(state_cyl, t, B0):
    
    r  = state_cyl[0]
    ur = state_cyl[1]
    phi  = state_cyl[2]
    uphi = state_cyl[3]
    z  = state_cyl[4]
    uz = state_cyl[5]
    
    durdt = uphi*B0/gamma(ur, uphi, uz)
    duphidt = -ur*B0//gamma(ur, uphi, uz)
    duzdt = 0.0
    
    drdt = ur/gamma(ur,uphi,uz)
    dphidt = uphi/(r*gamma(ur,uphi,uz))
    dzdt = uz/gamma(ur,uphi,uz)
    
    derivs = [drdt, durdt, dphidt, duphidt, dzdt, duzdt]
    
    return derivs


state0 = [1.0, 10.0, 0.0, 0.0, 0.0, 1.0]

t = np.linspace(0.0, 50.0, 1000)

state = odeint(syn, state0, t, args = (B0, ))

state_cyl = odeint(syn_cyl, state0, t, args = (B0, ))

x  = state[:,0]
ux = state[:,1]
y  = state[:,2]
uy = state[:,3]
z  = state[:,4]
uz = state[:,5]

r  = state_cyl[:,0]
ur = state_cyl[:,1]
phi  = state_cyl[:,2]
uphi = state_cyl[:,3]
z  = state_cyl[:,4]
uz = state_cyl[:,5]


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.plot(x, y, z)
ax.plot(r, phi, z)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.show()



