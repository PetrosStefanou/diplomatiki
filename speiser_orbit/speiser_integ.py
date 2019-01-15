#!/usr/bin/env python
# coding: utf-8


from scipy.integrate import odeint
import time
import numpy as np
import speiser_fun as sf
import speiser_fun_cyl as sfc
import pulsars








 



#ολοκλήρωση

# def oloklirosi1(gamma0, Delta, delta, B_0, gamma1, Rlc, t_end = 2., Dt = 10**4):

#     x, ux, y, uy, z, uz = [0]*len(gamma0), [0]*len(gamma0), [0]*len(gamma0), [0]*len(gamma0), [0]*len(gamma0), [0]*len(gamma0)

#     x1, ux1, y1, uy1, z1, uz1 = [0]*len(gamma0), [0]*len(gamma0), [0]*len(gamma0), [0]*len(gamma0), [0]*len(gamma0), [0]*len(gamma0)

#     r, ur, phi, uphi, z_cyl, uz_cyl = [0]*len(gamma0), [0]*len(gamma0), [0]*len(gamma0), [0]*len(gamma0), [0]*len(gamma0), [0]*len(gamma0)
    
#     r1, ur1, phi1, uphi1, z_cyl1, uz_cyl1 = [0]*len(gamma0), [0]*len(gamma0), [0]*len(gamma0), [0]*len(gamma0), [0]*len(gamma0), [0]*len(gamma0)

#     r2, ur2, phi2, uphi2, z_cyl2, uz_cyl2 = [0]*len(gamma0), [0]*len(gamma0), [0]*len(gamma0), [0]*len(gamma0), [0]*len(gamma0), [0]*len(gamma0)

#     #χρόνος της ολοκλήρωσης, σε μονάδες [qB0/mc]
#     t = np.linspace(0.0, t_end*Delta, Dt)
    
#     for i in range(0, len(gamma0)-1):
        
#         global Frad_cart
#         global Frad_cyl
        
#         Frad_cart = [0.]
#         Frad_cyl = [0.]
#         #αρχικές συνθήκες
#         state0 = np.array([0.0, 0.0, delta, -np.sqrt(gamma0[i]**2 - 1), Rlc - Delta, 0.0])
#         state0_cyl = np.array([Rlc - Delta, 0.0, 0.0, 0.0, delta, -np.sqrt(gamma0[i]**2 - 1)])

#         #ολοκλήρωση τροχιάς με απώλειες
#         state = odeint(speiser3D, state0, t, args = (Delta, delta, B_0, gamma1, ), mxstep = 3000, full_output=1)
    
#         #ολοκλήρωση τροχιάς χωρίς απώλειες
#         state1 = odeint(speiser3D_noloss, state0, t, args = (Delta, delta, B_0, gamma1, ), mxstep = 3000, full_output=1)
        
#         #ολοκλήρωση τρόχιάς σε κυλινδρικές συντεταγμένες
#         state_cyl = odeint(speiser3D_cyl, state0_cyl, t, args = (Delta, delta, B_0, gamma1, ), mxstep = 3000, full_output=1)

#         state_cyl1 = odeint(speiser3D_cyl_noloss, state0_cyl, t, args = (Delta, delta, B_0, gamma1, ), mxstep = 3000, full_output=1)

#         state_cyl2 = odeint(speiser3D_cyl_Rlc, state0_cyl, t, args = (Delta, delta, B_0, gamma1, ), mxstep = 3000, full_output=1)

#         x[i], ux[i], y[i], uy[i], z[i], uz[i] = state[0][:,0], state[0][:,1], state[0][:,2], state[0][:,3], state[0][:,4], state[0][:,5]  
    
#         x1[i], ux1[i], y1[i], uy1[i], z1[i], uz1[i] = state1[0][:,0], state1[0][:,1], state1[0][:,2], state1[0][:,3], state1[0][:,4], state1[0][:,5] 

#         r[i], ur[i], phi[i], uphi[i], z_cyl[i], uz_cyl[i] = state_cyl[0][:,0], state_cyl[0][:,1], state_cyl[0][:,2], state_cyl[0][:,3], state_cyl[0][:,4], state_cyl[0][:,5]

#         r1[i], ur1[i], phi1[i], uphi1[i], z_cyl1[i], uz_cyl1[i] = state_cyl1[0][:,0], state_cyl1[0][:,1], state_cyl1[0][:,2], state_cyl1[0][:,3], state_cyl1[0][:,4], state_cyl1[0][:,5]

#         r2[i], ur2[i], phi2[i], uphi2[i], z_cyl2[i], uz_cyl2[i] = state_cyl2[0][:,0], state_cyl2[0][:,1], state_cyl2[0][:,2], state_cyl2[0][:,3], state_cyl2[0][:,4], state_cyl2[0][:,5]
    
#     return np.array([x, ux, y, uy, z, uz]), np.array([x1, ux1, y1, uy1, z1, uz1]), np.array([r, ur, phi, uphi, z_cyl, uz_cyl]), np.array([r1, ur1, phi1, uphi1, z_cyl1, uz_cyl1]), np.array([r2, ur2, phi2, uphi2, z_cyl2, uz_cyl2]), [state[1], state1[1], state_cyl[1], state_cyl1[1], state_cyl[1]]

def oloklirosi(gamma0, Rlc, Delta, delta, B_0, t_end, Dt, system, coord):
    
    x1, u1, x2, u2, x3, u3 = [0]*len(gamma0), [0]*len(gamma0), [0]*len(gamma0), [0]*len(gamma0), [0]*len(gamma0), [0]*len(gamma0)

    t = np.linspace(0.0, t_end*Delta, Dt)

    for i in range(0, len(gamma0)-1):
        
        global Frad
        Frad = [0.]

        #αρχικές συνθήκες
        if coord == 'cart':
            
            init = np.array([0.0, 0.0, delta, -np.sqrt(gamma0[i]**2 - 1), Rlc - Delta, 0.0])
        
        elif coord == 'cyl':
            
            init = np.array([Rlc - Delta, 0.0, 0.0, 0.0, delta, -np.sqrt(gamma0[i]**2 - 1)])
        
        #ολοκλήρωση τροχιάς
        state = odeint(system, init, t, args = (Rlc, Delta, delta, B_0, Frad, ), mxstep = 3000, full_output=1)

        x1[i], u1[i], x2[i], u2[i], x3[i], u3[i] = state[0][:,0], state[0][:,1], state[0][:,2], state[0][:,3], state[0][:,4], state[0][:,5]

        dic = state[1]

    return np.array([x1, u1, x2, u2, x3, u3]), dic