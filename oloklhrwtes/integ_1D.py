#!/usr/bin/env python
# coding: utf-8

# In[19]:


#euler method

from pylab import *

def PitchforkODE(r,x):
    return r*x-x**3

def OneDEuler(r,x,f,dt):
    return x + dt*f(r,x)

dt = 0.2

r = 1.0

x1 = [0.1]
x2 = [-0.1]
x3 = [2.1]
x4 = [-2.1]

t = [0.0]

N = 50

for n in range(0,N):
    x1.append(OneDEuler(r,x1[n],PitchforkODE,dt))
    x2.append(OneDEuler(r,x2[n],PitchforkODE,dt))
    x3.append(OneDEuler(r,x3[n],PitchforkODE,dt))
    x4.append(OneDEuler(r,x4[n],PitchforkODE,dt))
    t.append(t[n] + dt)

    
xlabel('Time t')
ylabel('x(t)')
axis([0.0,dt*(N + 1),-2.0,2.0])
plot(t,x1,'b')
plot(t,x2,'r')
plot(t,x3,'g')
plot(t,x4,'m')


# In[20]:


#improved euler method

from pylab import *

def PitchforkODE(r,x):
    return r*x-x**3

def ImprovedOneDEuler(r,x,f,dt):
    xtemp = x + dt*f(r,x)
    return x + dt*(f(r,x)+f(r,xtemp))/2.0
dt = 0.2

r = 1.0

x1 = [0.1]
x2 = [-0.1]
x3 = [2.1]
x4 = [-2.1]

t = [0.0]

N = 50

for n in range(0,N):
    x1.append(ImprovedOneDEuler(r,x1[n],PitchforkODE,dt))
    x2.append(ImprovedOneDEuler(r,x2[n],PitchforkODE,dt))
    x3.append(ImprovedOneDEuler(r,x3[n],PitchforkODE,dt))
    x4.append(ImprovedOneDEuler(r,x4[n],PitchforkODE,dt))
    t.append(t[n] + dt)
    
xlabel('Time t')
ylabel('x(t)')
axis([0.0,dt*(N + 1),-2.0,2.0])
plot(t,x1,'b')
plot(t,x2,'r')
plot(t,x3,'g')
plot(t,x4,'m')


# In[21]:


#runge-kutta

from pylab import *

def PitchforkODE(r,x):
    return r*x-x**3

def RKOneD(r,x,f,dt):
    k1 = dt*f(r,x)
    k2 = dt*f(r,x+k1/2.)
    k3 = dt*f(r,x+k2/2.)
    k4 = dt*f(r,x+k3)
    
    return x + (k1 + 2.*k2 + 2.*k3 + k4)/6.

dt = 0.2

r = 1.0

x1 = [0.1]
x2 = [-0.1]
x3 = [2.1]
x4 = [-2.1]

t = [0.0]

N = 50

for n in range(0,N):
    x1.append(ImprovedOneDEuler(r,x1[n],PitchforkODE,dt))
    x2.append(ImprovedOneDEuler(r,x2[n],PitchforkODE,dt))
    x3.append(ImprovedOneDEuler(r,x3[n],PitchforkODE,dt))
    x4.append(ImprovedOneDEuler(r,x4[n],PitchforkODE,dt))
    t.append(t[n] + dt)
    
xlabel('Time t')
ylabel('x(t)')
axis([0.0,dt*(N + 1),-2.0,2.0])
plot(t,x1,'b')
plot(t,x2,'r')
plot(t,x3,'g')
plot(t,x4,'m')


# In[ ]:




