#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


# In[2]:


B = 1.0
E = 0.00


# In[3]:


def gamma(ux,uy,uz):
    g = np.sqrt(1 + ux**2 + uy**2 + uz**2)
    return g


# In[4]:


def syn(state, t):
    
    x  = state[0]
    ux = state[1]
    y  = state[2]
    uy = state[3]
    z  = state[4]
    uz = state[5]
    
    duxdt = uy
    duydt = -ux
    duzdt = gamma(ux,uy,uz)*E
    
    dxdt = ux/gamma(ux,uy,uz)
    dydt = uy/gamma(ux,uy,uz)
    dzdt = uz/gamma(ux,uy,uz)
    
    derivs = [dxdt, duxdt, dydt, duydt, dzdt, duzdt]
    
    return derivs


state0 = [0.0, 10.0, 0.0, 0.0, 0.0, 1.0]

t = np.linspace(0.0, 50.0, 1000)

state = odeint(syn, state0, t)

x  = state[:,0]
ux = state[:,1]
y  = state[:,2]
uy = state[:,3]
z  = state[:,4]
uz = state[:,5]


# In[5]:


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.plot(x, y, z)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
# ax.plot(x, z, 'r+', zdir='x')


# In[6]:


plt.plot(x, y)


# In[7]:


r = np.sqrt(x**2 + (y+1)**2)
plt.plot(t, r)
plt.axis([0., 50., 0., 2.0])


# In[ ]:





# In[ ]:




