#!/usr/bin/env python
# coding: utf-8


from scipy.integrate import odeint
import time
import numpy as np
import speiser_fun as sf
import speiser_fun_cyl as sfc
import pulsars



#σύστημα διαφορικών με απώλειες ακτινοβολίας
def speiser3D(state, t, Delta, delta, B_0, gamma1, q = 1):
    
    x, ux, y, uy, z, uz  = state
    
    #αδιάστατες εξισώσεις, οι πραγματικές έχουν διαιρεθεί με c*ω0 = mc/eB0c
    duxdt =  (q*sf.Flor_x(y, z, ux, uy, uz, Delta, delta) -Frad_cart[-1]*ux)*gamma1
    duydt =  (q*sf.Flor_y(y, z, ux, uy, uz, Delta, delta) -Frad_cart[-1]*uy)*gamma1
    duzdt =  (q*sf.Flor_z(y, z, ux, uy, uz, Delta, delta) -Frad_cart[-1]*uz)*gamma1
    
    dxdt = ux/sf.gamma(ux,uy,uz)
    dydt = uy/sf.gamma(ux,uy,uz)
    dzdt = uz/sf.gamma(ux,uy,uz)
    
    Frad_cart.append(sf.F_rad(ux, uy, uz, duxdt, duydt, duzdt, B_0))
    
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

#σύστημα διαφορικών με απώλειες ακτινοβολίας με ακτίνα καμπυλότητας υπολογισμένη από την τροχιά σε κυλινδρικές
def speiser3D_cyl(state_cyl, t, Delta, delta, B_0, gamma1, q = 1):
    
    r, ur, phi, uphi, z_cyl, uz_cyl = state_cyl
    
    #αδιάστατες εξισώσεις, οι πραγματικές έχουν διαιρεθεί με c*ω0 = mc/eB0c
    durdt =  sfc.Flor_r(r, phi, z_cyl, ur, uphi, uz_cyl, Delta, delta) + uphi**2/(r*sfc.gamma(ur, uphi, uz_cyl)) - Frad_cyl[-1]*ur 
    duphidt =  sfc.Flor_phi(r, phi, z_cyl, ur, uphi, uz_cyl, Delta, delta) - ur*uphi/(r*sfc.gamma(ur, uphi, uz_cyl)) - Frad_cyl[-1]*uphi 
    duzdt =  sfc.Flor_z_cyl(r, phi, z_cyl, ur, uphi, uz_cyl, Delta, delta) - Frad_cyl[-1]*uz_cyl
    
    drdt = ur/sfc.gamma(ur,uphi,uz_cyl)
    dphidt = uphi/(r*sfc.gamma(ur,uphi,uz_cyl))
    dzdt = uz_cyl/sfc.gamma(ur,uphi,uz_cyl)
    
    Frad_cyl.append(sfc.F_rad(r, ur, uphi, uz_cyl, durdt, duphidt, duzdt, B_0, sfc.Ploss))
    
    derivs = np.array([drdt, durdt, dphidt, duphidt, dzdt, duzdt])

    # print(sfc.Rc(r, ur, uphi, uz_cyl, durdt, duphidt, duzdt))  
    # print(sf.gamma(ux, uy, uz), np.sqrt(ux**2 + uy**2 + uz**2), sf.Ez(z,Delta))
    # print('{:1.2E}'.format(Frad_cyl[-1]))
    return derivs

    #σύστημα διαφορικών με απώλειες ακτινοβολίας με ακτίνα καμπυλότητας Rlc σε κυλινδρικές
def speiser3D_cyl_Rlc(state_cyl, t, Delta, delta, B_0, gamma1, q = 1):
    
    r, ur, phi, uphi, z_cyl, uz_cyl = state_cyl
    
    #αδιάστατες εξισώσεις, οι πραγματικές έχουν διαιρεθεί με c*ω0 = mc/eB0c
    durdt =  sfc.Flor_r(r, phi, z_cyl, ur, uphi, uz_cyl, Delta, delta) + uphi**2/(r*sfc.gamma(ur, uphi, uz_cyl)) - Frad_cyl[-1]*ur 
    duphidt =  sfc.Flor_phi(r, phi, z_cyl, ur, uphi, uz_cyl, Delta, delta) - ur*uphi/(r*sfc.gamma(ur, uphi, uz_cyl)) - Frad_cyl[-1]*uphi 
    duzdt =  sfc.Flor_z_cyl(r, phi, z_cyl, ur, uphi, uz_cyl, Delta, delta) - Frad_cyl[-1]*uz_cyl
    
    drdt = ur/sfc.gamma(ur,uphi,uz_cyl)
    dphidt = uphi/(r*sfc.gamma(ur,uphi,uz_cyl))
    dzdt = uz_cyl/sfc.gamma(ur,uphi,uz_cyl)
    
    Frad_cyl.append(sfc.F_rad(r, ur, uphi, uz_cyl, durdt, duphidt, duzdt, B_0, sfc.Ploss_Rlc))
    
    derivs = np.array([drdt, durdt, dphidt, duphidt, dzdt, duzdt])

    # print(sfc.Rc(r, ur, uphi, uz_cyl, durdt, duphidt, duzdt))  
    # print(sf.gamma(ux, uy, uz), np.sqrt(ux**2 + uy**2 + uz**2), sf.Ez(z,Delta))
    # print('{:1.2E}'.format(Frad_cyl[-1]))
    return derivs

#σύστημα διαφορικών χωρίς απώλειες ακτινοβολίας σε κυλινδρικές
def speiser3D_cyl_noloss(state_cyl, t, Delta, delta, B_0, gamma1, q = 1):
    
    r, ur, phi, uphi, z_cyl, uz_cyl = state_cyl
    
    #αδιάστατες εξισώσεις, οι πραγματικές έχουν διαιρεθεί με c*ω0 = mc/eB0c
    durdt =  sfc.Flor_r(r, phi, z_cyl, ur, uphi, uz_cyl, Delta, delta) + uphi**2/(r*sfc.gamma(ur, uphi, uz_cyl)) 
    duphidt =  sfc.Flor_phi(r, phi, z_cyl, ur, uphi, uz_cyl, Delta, delta) - ur*uphi/(r*sfc.gamma(ur, uphi, uz_cyl))
    duzdt =  sfc.Flor_z_cyl(r, phi, z_cyl, ur, uphi, uz_cyl, Delta, delta)
    
    drdt = ur/sfc.gamma(ur,uphi,uz_cyl)
    dphidt = uphi/(r*sfc.gamma(ur,uphi,uz_cyl))
    dzdt = uz_cyl/sfc.gamma(ur,uphi,uz_cyl)
    
    
    derivs = np.array([drdt, durdt, dphidt, duphidt, dzdt, duzdt])

    # print(sfc.Rc(r, ur, uphi, uz_cyl, durdt, duphidt, duzdt))  
    # print(sf.gamma(ux, uy, uz), np.sqrt(ux**2 + uy**2 + uz**2), sf.Ez(z,Delta))
    
    return derivs

#ολοκλήρωση

def oloklirosi(gamma0, Delta, delta, B_0, gamma1, t_end = 2., Dt = 10**4):

    x, ux, y, uy, z, uz = [0]*len(gamma0), [0]*len(gamma0), [0]*len(gamma0), [0]*len(gamma0), [0]*len(gamma0), [0]*len(gamma0)

    x1, ux1, y1, uy1, z1, uz1 = [0]*len(gamma0), [0]*len(gamma0), [0]*len(gamma0), [0]*len(gamma0), [0]*len(gamma0), [0]*len(gamma0)

    r, ur, phi, uphi, z_cyl, uz_cyl = [0]*len(gamma0), [0]*len(gamma0), [0]*len(gamma0), [0]*len(gamma0), [0]*len(gamma0), [0]*len(gamma0)
    
    r1, ur1, phi1, uphi1, z_cyl1, uz_cyl1 = [0]*len(gamma0), [0]*len(gamma0), [0]*len(gamma0), [0]*len(gamma0), [0]*len(gamma0), [0]*len(gamma0)

    r2, ur2, phi2, uphi2, z_cyl2, uz_cyl2 = [0]*len(gamma0), [0]*len(gamma0), [0]*len(gamma0), [0]*len(gamma0), [0]*len(gamma0), [0]*len(gamma0)

    #χρόνος της ολοκλήρωσης, σε μονάδες [qB0/mc]
    t = np.linspace(0.0, t_end*Delta, Dt)
    
    for i in range(0, len(gamma0)-1):
        
        global Frad_cart
        global Frad_cyl
        
        Frad_cart = [0.]
        Frad_cyl = [0.]
        #αρχικές συνθήκες
        state0 = np.array([0.0, 0.0, delta, -np.sqrt(gamma0[i]**2 - 1), 0.0, 0.0])
        state0_cyl = np.array([1.0, 0.0, 0.0, 0.0, delta, -np.sqrt(gamma0[i]**2 - 1)])

        #ολοκλήρωση τροχιάς με απώλειες
        state = odeint(speiser3D, state0, t, args = (Delta, delta, B_0, gamma1, ), mxstep = 2000, full_output=1)
    
        #ολοκλήρωση τροχιάς χωρίς απώλειες
        state1 = odeint(speiser3D_noloss, state0, t, args = (Delta, delta, B_0, gamma1, ), mxstep = 2000, full_output=1)
        
        #ολοκλήρωση τρόχιάς σε κυλινδρικές συντεταγμένες
        state_cyl = odeint(speiser3D_cyl, state0_cyl, t, args = (Delta, delta, B_0, gamma1, ), mxstep = 3000, full_output=1)

        state_cyl1 = odeint(speiser3D_cyl_noloss, state0_cyl, t, args = (Delta, delta, B_0, gamma1, ), mxstep = 3000, full_output=1)

        state_cyl2 = odeint(speiser3D_cyl_Rlc, state0_cyl, t, args = (Delta, delta, B_0, gamma1, ), mxstep = 3000, full_output=1)

        x[i], ux[i], y[i], uy[i], z[i], uz[i] = state[0][:,0], state[0][:,1], state[0][:,2], state[0][:,3], state[0][:,4], state[0][:,5]  
    
        x1[i], ux1[i], y1[i], uy1[i], z1[i], uz1[i] = state1[0][:,0], state1[0][:,1], state1[0][:,2], state1[0][:,3], state1[0][:,4], state1[0][:,5] 

        r[i], ur[i], phi[i], uphi[i], z_cyl[i], uz_cyl[i] = state_cyl[0][:,0], state_cyl[0][:,1], state_cyl[0][:,2], state_cyl[0][:,3], state_cyl[0][:,4], state_cyl[0][:,5]

        r1[i], ur1[i], phi1[i], uphi1[i], z_cyl1[i], uz_cyl1[i] = state_cyl1[0][:,0], state_cyl1[0][:,1], state_cyl1[0][:,2], state_cyl1[0][:,3], state_cyl1[0][:,4], state_cyl1[0][:,5]

        r2[i], ur2[i], phi2[i], uphi2[i], z_cyl2[i], uz_cyl2[i] = state_cyl2[0][:,0], state_cyl2[0][:,1], state_cyl2[0][:,2], state_cyl2[0][:,3], state_cyl2[0][:,4], state_cyl2[0][:,5]
    
    return np.array([x, ux, y, uy, z, uz]), np.array([x1, ux1, y1, uy1, z1, uz1]), np.array([r, ur, phi, uphi, z_cyl, uz_cyl]), np.array([r1, ur1, phi1, uphi1, z_cyl1, uz_cyl1]), np.array([r2, ur2, phi2, uphi2, z_cyl2, uz_cyl2]), [state[1], state1[1], state_cyl[1], state_cyl1[1], state_cyl[1]]



# x, ux, y, uy, z, uz, 
