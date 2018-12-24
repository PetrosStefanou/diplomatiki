#!/usr/bin/env python
# coding: utf-8


from scipy.integrate import odeint
import time
import numpy as np
import speiser_fun as sf
import pulsars



#σύστημα διαφορικών με απώλειες ακτινοβολίας
def speiser3D(state, t, Delta, delta, B_0, gamma1, q = 1):
    
    x, ux, y, uy, z, uz  = state
    
    #αδιάστατες εξισώσεις, οι πραγματικές έχουν διαιρεθεί με c*ω0 = mc/eB0c
    duxdt =  (q*sf.Flor_x(y, z, ux, uy, uz, Delta, delta) -Frad[-1]*ux)*gamma1
    duydt =  (q*sf.Flor_y(y, z, ux, uy, uz, Delta, delta) -Frad[-1]*uy)*gamma1
    duzdt =  (q*sf.Flor_z(y, z, ux, uy, uz, Delta, delta) -Frad[-1]*uz)*gamma1
    
    dxdt = ux/sf.gamma(ux,uy,uz)
    dydt = uy/sf.gamma(ux,uy,uz)
    dzdt = uz/sf.gamma(ux,uy,uz)
    
    Frad.append(sf.F_rad(ux, uy, uz, duxdt, duydt, duzdt, B_0))
    
    derivs = np.array([dxdt, duxdt, dydt, duydt, dzdt, duzdt])
    
#     print(Frad[-1]*ux, Frad[-1]*uy, Frad[-1]*uz, Ez(x))
    
    return derivs

#σύστημα διαφορικών χωρίς απώλειες ακτινοβολίας
def speiser3D_noloss(state1, t, Delta, delta, B_0, gamma1, q = 1):
    
    x, ux, y, uy, z, uz  = state1
    
    #αδιάστατες εξισώσεις, οι πραγματικές έχουν διαιρεθεί με c*ω0 = mc/eB0c
    duxdt =  sf.Flor_x(y, z, ux, uy, uz, Delta, delta)*gamma1 
    duydt =  sf.Flor_y(y, z, ux, uy, uz, Delta, delta)*gamma1
    duzdt =  sf.Flor_z(y, z, ux, uy, uz, Delta, delta)*gamma1
    
    dxdt = ux/sf.gamma(ux,uy,uz)
    dydt = uy/sf.gamma(ux,uy,uz)
    dzdt = uz/sf.gamma(ux,uy,uz)
    
    derivs = np.array([dxdt, duxdt, dydt, duydt, dzdt, duzdt])
      
    # print(sf.gamma(ux, uy, uz), np.sqrt(ux**2 + uy**2 + uz**2), sf.Ez(z,Delta))
    
    return derivs

#ολοκλήρωση

def oloklirosi(gamma0, Delta, delta, B_0, gamma1, t_end = 2., Dt = 10000):

    x, ux, y, uy, z, uz = [0.]*len(gamma0), [0.]*len(gamma0), [0.]*len(gamma0), [0.]*len(gamma0), [0.]*len(gamma0), [0.]*len(gamma0)

    x1, ux1, y1, uy1, z1, uz1 = [0.]*len(gamma0), [0.]*len(gamma0), [0.]*len(gamma0), [0.]*len(gamma0), [0.]*len(gamma0), [0.]*len(gamma0)

    #χρόνος της ολοκλήρωσης, σε μονάδες [qB0/mc]
    t = np.linspace(0.0, t_end*Delta, Dt)
    
    for i in range(0, len(gamma0)-1):
        
        global Frad
        
        Frad = [0.]
        
        #αρχικές συνθήκες
        state0 = np.array([0.0, 0.0, delta, -np.sqrt(gamma0[i]**2 - 1), 0.0, 0.0])
    
        #ολοκλήρωση τροχιάς με απώλειες
        state = odeint(speiser3D, state0, t, args = (Delta, delta, B_0, gamma1, ), full_output=0)
    
        #ολοκλήρωση τροχιάς χωρίς απώλειες
        state1 = odeint(speiser3D_noloss, state0, t, args = (Delta, delta, B_0, gamma1, ), full_output=0)
    
        x[i], ux[i], y[i], uy[i], z[i], uz[i] = state[:,0], state[:,1], state[:,2], state[:,3], state[:,4], state[:,5]  
    
        x1[i], ux1[i], y1[i], uy1[i], z1[i], uz1[i] = state1[:,0], state1[:,1], state1[:,2], state1[:,3], state1[:,4], state1[:,5] 
    
    return np.array([x, ux, y, uy, z, uz, x1, ux1, y1, uy1, z1, uz1])

# x, ux, y, uy, z, uz, 
