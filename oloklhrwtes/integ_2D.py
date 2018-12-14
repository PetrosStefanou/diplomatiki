#!/usr/bin/env python
# coding: utf-8

# In[3]:


#runge-kutta 2D

from pylab import *

def VDPXDot(u,x,y):
    return y

def VDPYDot(u,x,y):
    return u*(0.1-x**2)*y - x

def RK2D(r,x,y,f,g,dt):
    k1x = dt*f(r,x,y)
    k1y = dt*g(r,x,y)
    k2x = dt*f(r,x+k1x/2.,y+k1y/2.)
    k2y = dt*g(r,x+k1x/2.,y+k1y/2.)
    k3x = dt*f(r,x+k2x/2.,y+k2y/2.)
    k3y = dt*g(r,x+k2x/2.,y+k2y/2.)
    k4x = dt*f(r,x+k3x,y+k3y)
    k4y = dt*g(r,x+k3x,y+k3y)
    
    x = x + (k1x + 2.*k2x + 2.*k3x + k4x)/6.
    y = y + (k1y + 2.*k2y + 2.*k3y + k4y)/6.
    
    return x,y

dt = 0.05

u = 2.0

x1 = [0.1]
y1 = [0.1]
x2 = [2.0]
y2 = [-2.0]
x3 = [-2.0]
y3 = [2.0]

plot(x1,y1,'bo')
plot(x2,y2,'ro')
plot(x3,y3,'go')

t = [0.0]

N = 1000


for n in range(0,N):
    x,y = RK2D(u,x1[n],y1[n],VDPXDot,VDPYDot,dt)
    x1.append(x)
    y1.append(y)
    x,y = RK2D(u,x2[n],y2[n],VDPXDot,VDPYDot,dt)
    x2.append(x)
    y2.append(y)
    x,y = RK2D(u,x3[n],y3[n],VDPXDot,VDPYDot,dt)
    x3.append(x)
    y3.append(y)
    t.append(t[n] + dt)
    
# Setup the parametric plot
xlabel('x(t)') # set x-axis label
ylabel('y(t)') # set y-axis label
title('4th order Runge-Kutta Method: van der Pol ODE at u = ' + str(u)) # set plot title
axis('equal')
axis([-2.0,2.0,-2.0,2.0])
# Plot the trajectory in the phase plane
plot(x1,y1,'b')
plot(x2,y2,'r')
plot(x3,y3,'g')


# In[ ]:





# In[ ]:




