#!/usr/bin/env python
# coding: utf-8


from scipy.integrate import odeint
import time
import numpy as np
import speiser_fun as sf
import speiser_fun_cyl as sfc
import pulsars



def oloklirosi(gamma0, Rlc, Delta, delta_init, B_0, t, w, q, T, system, coord):
    
    x1, u1, x2, u2, x3, u3 = [0]*len(gamma0), [0]*len(gamma0), [0]*len(gamma0), [0]*len(gamma0), [0]*len(gamma0), [0]*len(gamma0)

    total = 0
    for i in range(0, len(gamma0)-1):
        start_time = time.time()
        global Frad
        Frad = [0.]
        
        u_drift = gamma0[i]/2.
        uphi_0 = u_drift - np.sqrt(gamma0[i]**2/2. - 1)
        uz_0 = -u_drift - np.sqrt(gamma0[i]**2/2. - 1)
        #αρχικές συνθήκες
        if coord == 'cart':
            
            init = np.array([0.0, 0.0, delta_init, -np.sqrt(gamma0[i]**2 - 1), Rlc - Delta, 0.0])
        
        elif coord == 'cyl':
            
            init = np.array([Rlc + w*Delta, 0.0, 0.0, uphi_0, delta_init, uz_0])
        
        #ολοκλήρωση τροχιάς
        perc = [0.]
        state = odeint(system, init, t, args = (Rlc, Delta, delta_init, B_0, Frad, q, T, perc), mxstep = 3000, full_output=1)
        
        x1[i], u1[i], x2[i], u2[i], x3[i], u3[i] = state[0][:,0], state[0][:,1], state[0][:,2], state[0][:,3], state[0][:,4], state[0][:,5]

        dic = state[1]
        
        elapsed = time.time() - start_time
        total += elapsed
        hours, rem = divmod(elapsed, 3600)
        minutes, seconds = divmod(rem, 60)
        print('runtime for gamma0 = {} is {}h {}m {}s'.format(int(gamma0[i]),int(hours), int(minutes), int(seconds)))
    
    hours, rem = divmod(total, 3600)
    minutes, seconds = divmod(rem, 60)
    print('total runtime = {}h {}m {}s'.format(int(hours), int(minutes), int(seconds)))
    
    return np.array([x1, u1, x2, u2, x3, u3]), dic