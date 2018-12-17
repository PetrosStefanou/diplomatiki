#!/usr/bin/env python
# coding: utf-8

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.signal import find_peaks
from pulsars import Pulsars, c, e_charge, e_mass
import speiser_fun as sf 
import speiser_integ 
import numpy as np
import time


k = 10**5
pulsar = Pulsars(k)

B_0 = pulsar['crab']['Blc']
omegaB = (e_charge*B_0/(e_mass*c))    #γυροσυχνότητα
Delta = pulsar['crab']['rlc']/k*omegaB/c

gamma0 = np.array([1000., 1.])    #αρχικός παράγοντας Lorentz
delta = 1000    #πάχος του φύλλου ρεύματος, αδιάστατο, σε μονάδες [c/ωΒ]

 #έναρξη χονομέτρησης
start_time = time.time()

#ολοκλήρωση
x, ux, y, uy, z, uz, x1, ux1, y1, uy1, z1, uz1 = speiser_integ.oloklirosi(gamma0, Delta, delta, B_0)

#λήξη χρονομέτρησης
elapsed = time.time() - start_time

print('total runtime = {:1.2E} s'.format(elapsed))

#τροχιά στο y-z επίπεδο
fig1 = plt.figure()
for i in range(0, len(x) - 1):
    plt.plot(z[i]*c/omegaB, y[i]*c/omegaB, label = 'losses, $\gamma_0 = %s$' %int(gamma0[i]))
    plt.plot(z1[i]*c/omegaB, y1[i]*c/omegaB, '--', label = 'no losses, $\gamma_0 = %s$' %int(gamma0[i]))
    peaks, _ = find_peaks(y[i], distance = 1000*c/omegaB)
    # plt.plot(z[i][peaks]*c/omegaB, y[i][peaks]*c/omegaB, '*', color = 'k')
plt.axhline(y = delta*c/omegaB, linestyle = ':', color = 'r')
plt.axhline(y = -delta*c/omegaB, linestyle = ':', color = 'r')
plt.axvline(x = Delta*c/omegaB, linestyle = ':', color = 'r')
plt.xlabel('z (cm)')
plt.ylabel('y (cm)')
plt.legend(loc = 'best')
plt.savefig('troxia zy')

# Εξέλιξη του παράγοντα Lorentz
fig2 = plt.figure()
for i in range(0, len(x)-1):
    plt.plot(z[i]/Delta, sf.gamma(ux[i], uy[i], uz[i]), label = 'losses')
    plt.plot(z1[i]/Delta, sf.gamma(ux1[i], uy1[i], uz1[i]), '--', label = 'no losses')
# plt.yscale('log')
plt.xlabel('z/Delta')
plt.ylabel('gamma')
plt.legend(loc = 'best')
plt.savefig('gamma')

