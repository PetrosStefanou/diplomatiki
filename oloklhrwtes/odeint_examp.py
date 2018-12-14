#!/usr/bin/env python
# coding: utf-8

# In[10]:


import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def MassSpring(state,t):
    # unpack the state vector
    x = state[0]
    xd = state[1]

    # these are our constants
    k = 2.5 # Newtons per metre
    m = 1.5 # Kilograms
    g = 9.8 # metres per second

    # compute acceleration xdd
    xdd = ((-k*x)/m) + g

    # return the two state derivatives
    return [xd, xdd]

state0 = [0.0, 0.0]
t = np.linspace(0.0, 10.0, 100)

state = odeint(MassSpring, state0, t)

plt.plot(t, state)
# xlabel('TIME (sec)')
# ylabel('STATES')
# title('Mass-Spring System')
# legend(('$x$ (m)', '$\dot{x}$ (m/sec)'))


# In[4]:


xd


# In[ ]:




