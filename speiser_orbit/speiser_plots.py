# τροχιά στο y-z επίπεδο
from pulsars import Pulsars, c, e_charge, e_mass
import speiser_fun as sf 
import speiser_integ


import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#γραφική παράσταση της τροχιάς στο z-y επίπεδο
def plot_zy(x, y, z, ux, uy, uz, x1, y1, z1, ux1, uy1, uz1, omegaB, delta, Delta, gamma0, save = 0):
    fig = plt.figure()    
    for i in range(0, len(x) - 1):
        plt.plot(z[i]*c/omegaB, y[i]*c/omegaB, label = 'losses, $\gamma_0 = %s$' %int(gamma0[i]))
        plt.plot(z1[i]*c/omegaB, y1[i]*c/omegaB, '--', label = 'no losses, $\gamma_0 = %s$' %int(gamma0[i]))
        # peaks, _ = find_peaks(y[i], distance = 1000*c/omegaB)
        # plt.plot(z[i][peaks]*c/omegaB, y[i][peaks]*c/omegaB, '*', color = 'k')
        plt.axhline(y = delta*c/omegaB, linestyle = ':', color = 'r')
        plt.axhline(y = -delta*c/omegaB, linestyle = ':', color = 'r')
        plt.axvline(x = Delta*c/omegaB, linestyle = ':', color = 'r')
        plt.xlabel('z (cm)')
        plt.ylabel('y (cm)')
        plt.legend(loc = 'best')
        if save == True:
            plt.savefig('troxia{}_zy'.format(i))
            
#γραφική παράσταση της τροχιάς στο x-y επίπεδο    
def plot_zx(x, y, z, ux, uy, uz, x1, y1, z1, ux1, uy1, uz1, omegaB, delta, Delta, gamma0, save = 0):
    fig = plt.figure()
    for i in range(0, len(x) - 1):
        fig3 = plt.figure()
        plt.plot(z[i]*c/omegaB, x[i]*c/omegaB, label = '$\gamma_0$ = {}'.format(int(gamma0[i])))
        plt.plot(z1[i]*c/omegaB, x1[i]*c/omegaB, '--', label = 'no losses $\gamma_0$ = {}'.format(int(gamma0[i])))
        plt.xlabel('z (cm)')
        plt.ylabel('x (cm)')
        plt.title('Τροχιά στο επίπεδο x-z')
        plt.axvline(x = Delta*c/omegaB, linestyle = ':', color = 'r')
        plt.legend(loc = 'best')
        if save == True:
            fig3.savefig('troxia{}_zx'.format(i))

#γραφική παράσταση της τροχιάς στο z-y-x επίπεδο
def plot_zyx(x, y, z, ux, uy, uz, x1, y1, z1, ux1, uy1, uz1, omegaB, delta, Delta, gamma0, save = 0):
    fig = plt.figure()
    for i in range(0, len(x) - 1):
        fig4 = plt.figure()
        ax = fig4.add_subplot(111, projection='3d')
        ax.plot(z[i]*c/omegaB, x[i]*c/omegaB, y[i]*c/omegaB)
        ax.plot(z1[i]*c/omegaB, x1[i]*c/omegaB, y1[i]*c/omegaB)

        ax.set_xlabel('z')
        ax.set_ylabel('x')
        ax.set_zlabel('y')
        if save == True:
            fig4.savefig('troxia{}_xyz'.format(i))

# Γραφική παράσταση της εξέλιξης του παράγοντα Lorentz με το z
def plot_gamma(x, y, z, ux, uy, uz, x1, y1, z1, ux1, uy1, uz1, omegaB, delta, Delta, gamma0, save = 0):
    fig = plt.figure()
    for i in range(0, len(x1)-1):
        plt.plot(z[i]/Delta, sf.gamma(ux[i], uy[i], uz[i]), label = 'losses, $\gamma_0 = %s$' %int(gamma0[i]))
        plt.plot(z1[i]/Delta, sf.gamma(ux1[i], uy1[i], uz1[i]), '--', label = 'no losses, $\gamma_0 = %s$' %int(gamma0[i]))
        # plt.yscale('log')
        plt.xlabel('z/Delta')
        plt.ylabel('gamma')
        # plt.axes('equal')
        plt.axvline(x = 1*c/omegaB, linestyle = ':', color = 'r')
        plt.legend(loc = 'best')
        if save == True:
            plt.savefig('gamma{}'.format(i))

