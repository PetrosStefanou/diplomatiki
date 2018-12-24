#!/usr/bin/env python
# coding: utf-8


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.signal import find_peaks
import numpy as np
import time

from pulsars import Pulsars, c, e_charge, e_mass
import speiser_fun as sf 
import speiser_integ 
import speiser_plots as sp


k = 10**3
pulsar = Pulsars(k)

B_0 = pulsar['crab']['Blc']
gamma1 = 10**0
omegaB = (e_charge*B_0/(gamma1*e_mass*c))    #γυροσυχνότητα
Delta = pulsar['crab']['rlc']/k*omegaB/c

gamma0 = np.array([10000., 1.])    #αρχικός παράγοντας Lorentz
delta = 1000.    #πάχος του φύλλου ρεύματος, αδιάστατο, σε μονάδες [c/ωΒ]

 #έναρξη χονομέτρησης
start_time = time.time()

#ολοκλήρωση
x, ux, y, uy, z, uz, x1, ux1, y1, uy1, z1, uz1 = speiser_integ.oloklirosi(gamma0, Delta, delta, B_0, gamma1, t_end = 200., Dt = 10**6)
# x, ux, y, uy, z, uz, 

#λήξη χρονομέτρησης
elapsed = time.time() - start_time

print('total runtime = {:1.2E} s'.format(elapsed))

#γραφικές παραστάσεις
sp.plot_zy(x, y, z, ux, uy, uz, x1, y1, z1, ux1, uy1, uz1, omegaB, delta, Delta, gamma0)

sp.plot_zx(x, y, z, ux, uy, uz, x1, y1, z1, ux1, uy1, uz1, omegaB, delta, Delta, gamma0)

sp.plot_zyx(x, y, z, ux, uy, uz, x1, y1, z1, ux1, uy1, uz1, omegaB, delta, Delta, gamma0)

sp.plot_gamma(x, y, z, ux, uy, uz, x1, y1, z1, ux1, uy1, uz1, omegaB, delta, Delta, gamma0)

plt.show()