#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


# In[2]:


q = 1.0
m = 1.0
B = 1.0


# In[3]:


def cycl(state, t):
    x  = state[0]
    vx = state[1]
    y  = state[2]
    vy = state[3]
    z  = state[4]
    vz = state[5]
    
    dvxdt = q*B*vy/m
    dvydt = -q*B*vx/m
    dvzdt = 0
    dxdt = vx
    dydt = vy
    dzdt = vz
    
    derivs = [dxdt, dvxdt, dydt, dvydt, dzdt, dvzdt]
    
    return derivs

state0 = [0.0, 0.0, 0.0, 1.0, 0.0, 1.0]

t = np.linspace(0.0, 50.0, 1000)

state = odeint(cycl, state0, t)


# In[4]:


plt.plot(state[:, 0], state[:, 1])


# In[8]:


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.plot(state[:,0], state[:,2], state[:,4])
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')



# In[ ]:




