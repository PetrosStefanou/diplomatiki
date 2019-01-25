#!/usr/bin/env python
# coding: utf-8


import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import time
from scipy.special import kv
from scipy.integrate import quad

from pulsars import Pulsars, c, e_charge, e_mass, h
import speiser_fun as sf 
import speiser_fun_cyl as sfc
import speiser_integ as si
import speiser_plots as sp
import speiser_model as sm


k = 1000.
name = 'crab'
pulsar = Pulsars(k)[name]

gamma1 = 1.
B_0 = pulsar['Blc']
omegaB = (e_charge*B_0/(e_mass*c))    #γυροσυχνότητα
Rlc = pulsar['rlc']*omegaB/c
Delta = pulsar['rlc']/k*omegaB/c

gamma0 = np.linspace(1., 1000., 4)    #αρχικός παράγοντας Lorentz

# N = gamma0**(-1)
delta = 1000.    #πάχος του φύλλου ρεύματος, αδιάστατο, σε μονάδες [c/ωΒ]



#έναρξη χονομέτρησης
start_time = time.time()

T = 60.
N = 6*10**3
t = np.linspace(0.0, T*Delta, N)

(r, ur, phi, uphi, z_cyl, uz_cyl), dic_cyl = si.oloklirosi(gamma0, Rlc, Delta, delta, B_0, t, sm.speiser_cyl, coord = 'cyl')

# #λήξη χρονομέτρησης
elapsed = time.time() - start_time

print('total runtime for cylindrical with losses = {:1.2E} s'.format(elapsed))



fig1, ax1 = plt.subplots()
for i in range(0, len(r)-1):
    ax1.plot(r[i]*c/omegaB, sf.gamma(ur[i], uphi[i], uz_cyl[i]), label = '$\gamma_0$ = {}'.format(int(gamma0[i])))
    ax1.set(xlabel = 'r (cm)', ylabel = '$\gamma$', title = 'Evolution of Lorentz factor')
    ax1.axvline(x = Rlc*c/omegaB, linestyle = ':', color = 'k')
    ax1.axvline(x = (Rlc - Delta)*c/omegaB, linestyle = ':', color = 'k')
    ax1.legend(loc='center left', ncol=1, fancybox=True, shadow=True, bbox_to_anchor=(0.8, 0.5))
#     fig1.savefig('γ(r)_k_{:d}.png'.format(int(k)))

fig2, ax2 = plt.subplots()    
for i in range(0, len(r) - 1):
    ax2.plot(r[i]*c/omegaB, z_cyl[i]*c/omegaB, label = '$\gamma_0$ = {}'.format(int(gamma0[i])))
    ax2.axhline(y = delta*c/omegaB, linestyle = ':', color = 'k')
    ax2.axhline(y = -delta*c/omegaB, linestyle = ':', color = 'k')
    ax2.axvline(x = Rlc*c/omegaB, linestyle = ':', color = 'k')
    ax2.axvline(x = (Rlc - Delta)*c/omegaB, linestyle = ':', color = 'k')
    ax2.set(xlabel = 'r (cm)', ylabel = 'z (cm)', title = 'Orbit in r-z plane')
    ax2.legend(loc='center left', ncol=1, fancybox=True, shadow=True, bbox_to_anchor=(0.8, 0.5))
#     fig2.savefig('epipedo_r-z_k_{:d}.png'.format(int(k))) 
    
fig3, ax3 = plt.subplots()
for i in range(0, len(r) - 1):
    plot1 = ax3.plot(r[i]*np.cos(phi[i])*c/omegaB, r[i]*np.sin(phi[i])*c/omegaB, 
                     label = '$\gamma_0$ = {}'.format(int(gamma0[i])))
    
    ax3.set(xlabel = 'r (cm)', ylabel = 'x (cm)', title = 'Orbit in r-x plane', aspect = 'auto', 
            xlim = [1.59E+8, 1.596E+8], ylim = [-0.01E+8, 0.1E+8])
    
#     ax3.axvline(x = (Rlc - Delta)*c/omegaB*, linestyle = ':', color = 'k')
#     ax3.axvline(x = Rlc*c/omegaB, linestyle = ':', color = 'k')
    ax3.legend(loc='center left', ncol=1, fancybox=True, shadow=True, bbox_to_anchor=(0.8, 0.5))
    
    radius1 = mpatches.Arc((0.0, 0.0), 2*Rlc*c/omegaB, 2*Rlc*c/omegaB, theta1 = 0, theta2 = 4, color = 'k', ls = ':',
                           linewidth=1, fill=False)
    radius2 = mpatches.Arc((0.0, 0.0), 2*(Rlc - Delta)*c/omegaB, 2*(Rlc - Delta)*c/omegaB, theta1 = 0, theta2 = 4,
                           color = 'k', ls = ':', linewidth=1, fill=False)
    ax3.add_patch(radius1)
    ax3.add_patch(radius2)
    
#     fig3.savefig('epipedo_r-rφ_k_{:d}.png'.format(int(k)))

plt.show()

urdot = [0]*len(ur)
uphidot = [0]*len(ur)
uz_cyldot = [0]*len(ur)
los = [0]*len(ur)
nu = [0]*len(ur)
r_curv = [0]*len(ur)
wc = [0]*len(ur)
pow_per_freq = [0]*len(ur)

Frad_cyl = [0.]

for i in range(0, len(ur)):
    
    for j in range(1, len(uphi[i])):
        
        urdot[i] = np.zeros(len(ur[i]))
        uphidot[i] = np.zeros(len(ur[i]))
        uz_cyldot[i] = np.zeros(len(ur[i]))
        los[i] = np.zeros(len(ur[i]))
        nu[i] = np.zeros(len(ur[i]))
        r_curv[i] = np.zeros(len(ur[i]))
        
        urdot[i][j] = sfc.Flor_r(r[i][j], phi[i][j], z_cyl[i][j], ur[i][j], uphi[i][j], uz_cyl[i][j], Rlc, Delta, delta) - uphi[i][j]**2/(r[i][j]*sfc.gamma(ur[i][j], uphi[i][j], uz_cyl[i][j])) - Frad_cyl[-1]*ur[i][j]
        uphidot[i][j] = sfc.Flor_phi(r[i][j], phi[i][j], z_cyl[i][j], ur[i][j], uphi[i][j], uz_cyl[i][j], Rlc, Delta, delta) + ur[i][j]*uphi[i][j]/(r[i][j]*sfc.gamma(ur[i][j], uphi[i][j], uz_cyl[i][j])) - Frad_cyl[-1]*uphi[i][j]
        uz_cyl[i][j] = sfc.Flor_z_cyl(r[i][j], phi[i][j], z_cyl[i][j], ur[i][j], uphi[i][j], uz_cyl[i][j], Rlc, Delta, delta) - Frad_cyl[-1]*uz_cyl[i][j]
        
        nu[i][j] = sfc.nu_crit_curv(r[i][j], ur[i][j], uphi[i][j], uz_cyl[i][j], urdot[i][j], uphidot[i][j], uz_cyldot[i][j])
        
        r_curv[i][j] = sfc.Rc(r[i][j], ur[i][j], uphi[i][j], uz_cyl[i][j], urdot[i][j], uphidot[i][j], uz_cyldot[i][j])
        
#         los[i][j] = sfc.Ploss(r[i][j], ur[i][j], uphi[i][j], uz_cyl[i][j], urdot[i][j], uphidot[i][j], uz_cyldot[i][j], B_0)
        
        Frad_cyl.append(sfc.F_rad(r[i][j], ur[i][j], uphi[i][j], uz_cyl[i][j], urdot[i][j], uphidot[i][j], uz_cyldot[i][j], B_0, sfc.Ploss))


def f(n, nc): 
    I = quad(lambda x: kv(5./3., x), n/nc, +np.inf)[0]
    return I

def sing_en_spec(nc, gamma, rc):
    
    nu = np.linspace(0.2*nc, 5*nc, 100)
    
    p = np.zeros(len(nu))
    
    for j in range(0, len(nu)):
        
        p[j] = (np.sqrt(3)*e_charge**2)/(2*np.pi)*gamma*nu[j]/(rc*nc)*f(nu[j], nc)
    
    return p

nu_crit = [0]*len(gamma0)
pow_per_freq = [0]*len(gamma0)

for i in range(0, len(gamma0)):
    
    nu_crit[i] = np.zeros(len(ur[i]))
    pow_per_freq[i] = np.zeros(len(ur[i]))
    
    for j in range(0, len(ur[i])):
        
        nu_crit[i][j] = 3*c*sfc.gamma(ur[i][j], uphi[i][j], uz_cyl[i][j])**3/2/r_curv[i][j]

        pow_per_freq[i][j] = sing_en_spec(nu_crit[i][j], sfc.gamma(ur[i][j], uphi[i][j], uz_cyl[i][j]), r_curv[i][j])
    
print(pow_per_freq)
#     for n in range(0, len(nu)):
#         pow_per_freq[i][n] = sing_en_spec(nu[n], nu_crit[i], sfc.gamma(ur[i], uphi[i], uz_cyl[i]), r_curv[i])