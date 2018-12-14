#!/usr/bin/env python
# coding: utf-8

# In[1]:


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.signal import find_peaks


# In[2]:


#παράμετροι
get_ipython().run_line_magic('matplotlib', 'notebook')

k = 5E+2    #multiplicity

get_ipython().run_line_magic('run', 'pulsars.ipynb')
Pulsars(k)

B_0 = pulsar['crab']['Blc']    #χαρκτηριτικό μαγνητικό πεδίο
omegaB = (e_charge*B_0/(e_mass*c))    #γυροσυχνότητα

gamma0 = np.array([1000., 1.])    #αρχικός παράγοντας Lorentz

delta = 1000    #πάχος του φύλλου ρεύματος, αδιάστατο, σε μονάδες [c/ωΒ]
Delta = pulsar['crab']['rlc']/k*omegaB/c     #Μήκος του φύλλου ρεύματος, αδιάστατο, σε μονάδες [c/ωΒ]  

    #Μαγνητικο πεδίο στον άξονα y, αδιάστατο, σε μονάδες [Β0]

q = 1    # + ή - για ποσιτρόνιο ή ηλεκτρόνιο

get_ipython().run_line_magic('run', 'speiser_integ.ipynb')


# In[3]:



for j in range(1):
    #χρονομέτσηση ολοκλήρωσης
    start_time = time.time()

    #ολοκλήρωση
    x, ux, y, uy, z, uz, x1, ux1, y1, uy1, z1, uz1 = oloklirosi(gamma0)

    elapsed = time.time() - start_time

    

    delta_i = np.zeros(len(gamma0))
   
    #Τροχια y-z
    fig = plt.figure()
    for i in range(0, len(x) - 1):
        plt.plot(z[i]*c/omegaB, y[i]*c/omegaB, label = 'losses, $\gamma_0 = %s$' %int(gamma0[i]))
    #     plt.plot(z1[i]*c/omegaB, y1[i]*c/omegaB, '--', label = 'no losses, $\gamma_0 = %s$' %int(gamma0[i]))
#     peaks, _ = find_peaks(y[2], distance = 1000*c/omegaB)
        peaks, _ = find_peaks(y[i])
        plt.plot(z[i][peaks]*c/omegaB, y[i][peaks]*c/omegaB, '*', color = 'k')
        plt.axhline(y = delta*c/omegaB, linestyle = ':', color = 'r')
        plt.axhline(y = -delta*c/omegaB, linestyle = ':', color = 'r')
        delta_i[i] = np.mean(y[i][peaks])
    
    plt.axvline(x = Delta*c/omegaB, linestyle = ':', color = 'r')
    plt.xlabel('z (cm)')
    plt.ylabel('y (cm)')
    plt.legend(loc = 'best')
    fig.savefig('troxia{}.jpeg'.format(j))
    plt.close()
    
        
    

    delta = np.mean(delta_i)
    
    print('runtime = {:1.2e} s, delta = {:1.2f}'.format(elapsed, delta*c/omegaB))


# In[6]:


#Τροχια y-z
fig = plt.figure()
for i in range(0, len(x) - 1):
    plt.plot(z[i]*c/omegaB, y[i]*c/omegaB, label = 'losses, $\gamma_0 = %s$' %int(gamma0[i]))
#     plt.plot(z1[i]*c/omegaB, y1[i]*c/omegaB, '--', label = 'no losses, $\gamma_0 = %s$' %int(gamma0[i]))
    peaks, _ = find_peaks(y[i], distance = 1000*c/omegaB)
    plt.plot(z[i][peaks]*c/omegaB, y[i][peaks]*c/omegaB, '*', color = 'k')
plt.axhline(y = delta*c/omegaB, linestyle = ':', color = 'r')
plt.axhline(y = -delta*c/omegaB, linestyle = ':', color = 'r')
plt.axvline(x = Delta*c/omegaB, linestyle = ':', color = 'r')
plt.xlabel('z (cm)')
plt.ylabel('y (cm)')
plt.legend(loc = 'best')


# In[9]:


mean0 = np.mean(y[0][peaks]*c/omegaB)
mean1 = np.mean(y[0][peaks][10:]*c/omegaB)
print(mean0, mean1)
fig = plt.figure()
for t in range(0, len(y[0])-1):
    plt.plot(y[0][t], Bx(y[0][t], z[0][t]))
    


# In[ ]:


delta = np.zeros(len(gamma0))
delta


# In[ ]:


#Εξέλιξη του παράγοντα Lorentz
fig = plt.figure()
for i in range(0, len(x)-1):
    plt.plot(z[i]/Delta, gamma(ux[i], uy[i], uz[i]), label = 'losses')
    plt.plot(z1[i]/Delta, gamma(ux1[i], uy1[i], uz1[i]), '--', label = 'no losses')
# plt.yscale('log')
plt.xlabel('z/Delta')
plt.ylabel('gamma')
plt.legend(loc = 'best')


# In[ ]:



fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

for i in range(0, len(x)-1):
    ax.plot(z[i]/Delta, x[i]*c/omegaB, y[i]*c/omegaB)
    ax.plot(z1[i]/Delta, x1[i]*c/omegaB, y1[i]*c/omegaB, '--')

ax.set_xlabel('z')
ax.set_ylabel('x')
ax.set_zlabel('y')

# fig.savefig('troxia_xyz')


# In[ ]:


print(' Μήκος Δ = {:1.2E}cm \n Γυροσυχνότητα ωΒ = {:1.2E}Hz\n Γυροακτίνα rB = {:1.2E} cm'.format(2*Delta*omegaB/c, omegaB, c/omegaB))


# In[ ]:


fig = plt.figure()
plt.plot(z[0]*c/omegaB,x[0]*c/omegaB, ':', label = '$\gamma_0 = 10$')
# plt.plot(z[1],x[1], '--', label = '$\gamma_0 = 100$')
# plt.plot(z[2],x[2], label = '$\gamma_0 = 1000$')
# plt.plot(z[3],x[3], '-.', label = '$\gamma_0 = 10000$')
plt.xlabel('z')
plt.ylabel('x')
# fig.savefig('troxia_xz')


# In[15]:


xx = np.linspace(0,10,100)
fig = plt.figure()
plt.plot(xx, np.tanh(xx)-1./xx)
plt.plot(xx, np.tanh(xx)-1.*xx)
plt.plot(xx, np.tanh(xx))


# In[ ]:




