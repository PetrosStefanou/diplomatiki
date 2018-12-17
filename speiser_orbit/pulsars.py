#!/usr/bin/env python
# coding: utf-8


import math
import numpy as np

#σταθερές
c = 2.99792458*10**10
r_sur = 1.5*10**6
rest_energy = 9.1093897*10**(-28)*(2.99792458*10**10)**2
e_charge = 4.8032068*10**(-10)
e_mass = 9.1093897*10**(-28)

#παράμετροι
# k = 1000      #multiplicity, αντιστοιχεί σε young pulsar
gamma_inj = 10    #αρχικός παράγοντας λόρεντζ των ποζιτρονίων
beta_rec = 0.1    #παράμετρος που καθορίζει το ηλεκτρικό πεδίο αν το dissipation προέρχεται από reconnection

#pulsars
crab = {'P':0.03339, 'B_sur':3.79*10**12}
vela = {'P':0.08933, 'B_sur':3.38*10**12}
geminga = {'P':0.2371, 'B_sur':1.63*10**12}
old_slow = {'P':1.24442, 'B_sur':3.01*10**12}
old_fast = {'P':0.05920, 'B_sur':4.92*10**9}
ms = {'P':0.01, 'B_sur':10**9}


pulsar = {'crab':crab, 'vela':vela, 'geminga':geminga, 'old_slow':old_slow, 'old_fast':old_fast, 'ms': ms}
#πρόσθεση καινούργιου μεγέθους
def Pulsars(k):
    for i in pulsar:

        #κύλινδρος φωτός
        pulsar[i]['rlc'] =  3*10**10*pulsar[i]['P']/(2*math.pi)

        #μαγνητικό πεδίο στον κύλινδρο φωτός
        pulsar[i]['Blc'] =  pulsar[i]['B_sur']*(r_sur/pulsar[i]['rlc'])**3/2

        #Γ_rrl για ακτινοβολία curvature
        #pulsar[i]['gamma_rrl_curv'] = 4*10**7*(pulsar[i]['B_sur']/10**13)**(0.25)*pulsar[i]['P']**(-0.25)

        #Γ_rrl για ακτινοβολία σύγχροτρον
        pulsar[i]['gamma_rrl_syn'] = 1.2*10**12*pulsar[i]['P']**(3/2)*pulsar[i]['B_sur']**(-0.5)

        #χαρακτηριστικός χρόνος επιτάχυνσης μέχρι το limit της curvature
        #pulsar[i]['t_acc_curv'] = pulsar[i]['gamma_rrl_curv']*rest_energy/(e_charge*3*10**10*pulsar[i]['Blc'])

        #χαρακτηριστικός χρόνος επιτάχυνσης μέχρι το limit της συγχροτρον
        pulsar[i]['t_acc_syn'] = pulsar[i]['gamma_rrl_syn']*rest_energy/(e_charge*3*10**10*pulsar[i]['Blc'])

        #χαρακτηριστικός χρόνος παραμονής στην περιοχή επιτάχυνσης
        pulsar[i]['t_acc_remain'] = pulsar[i]['rlc']/(k*c)

        #γυροπερίοδος των φορτίων που εισέρχονται στο φύλλο με γ_injection
        #pulsar[i]['gperiod_inj'] = 2*math.pi*gamma_inj*e_mass*c/(e_charge*pulsar[i]['Blc'])

        #μέγιστος παράγοντας Lorentz που μπορεί να κερδίσει το σωματίδιο
        pulsar[i]['gamma_max'] = pulsar[i]['Blc']*pulsar[i]['rlc']/k*e_charge/rest_energy
    
#     for i in pulsar:
#         print('{}'.format(i))
#         for j in pulsar[i]:
#             print('\t{} = {:1.2E}\n'.format(j,pulsar[i][j]))
            
    return pulsar






