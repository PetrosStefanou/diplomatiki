#!/usr/bin/env python
# coding: utf-8

# In[1]:


from scipy.integrate import odeint
import time


# In[2]:


#παράμετροι

get_ipython().run_line_magic('run', 'speiser_fun.ipynb')
get_ipython().run_line_magic('matplotlib', 'notebook')


# In[3]:


#σύστημα διαφορικών με απώλειες ακτινοβολίας
def speiser3D(state, t):
    
    x, ux, y, uy, z, uz  = state
    
    #αδιάστατες εξισώσεις, οι πραγματικές έχουν διαιρεθεί με c*ω0 = mc/eB0c
    duxdt = -q*Flor_x(y, z, ux, uy, uz) -Frad[-1]*ux
    duydt =  q*Flor_y(y, z, ux, uy, uz) -Frad[-1]*uy
    duzdt =  q*Flor_z(y, z, ux, uy, uz) -Frad[-1]*uz
    
    dxdt = ux/gamma(ux,uy,uz)
    dydt = uy/gamma(ux,uy,uz)
    dzdt = uz/gamma(ux,uy,uz)
    
    Frad.append(F_rad(ux, uy, uz, duxdt, duydt, duzdt))
    
    derivs = np.array([dxdt, duxdt, dydt, duydt, dzdt, duzdt])
    
#     print(Frad[-1]*ux, Frad[-1]*uy, Frad[-1]*uz, Ez(x))
    
    return derivs


# In[4]:


#σύστημα διαφορικών χωρίς απώλειες ακτινοβολίας
def speiser3D_noloss(state1, t):
    
    x, ux, y, uy, z, uz  = state1
    
    #αδιάστατες εξισώσεις, οι πραγματικές έχουν διαιρεθεί με c*ω0 = mc/eB0c
    duxdt = -Flor_x(y, z, ux, uy, uz) 
    duydt =  Flor_y(y, z, ux, uy, uz)
    duzdt =  Flor_z(y, z, ux, uy, uz)
    
    dxdt = ux/gamma(ux,uy,uz)
    dydt = uy/gamma(ux,uy,uz)
    dzdt = uz/gamma(ux,uy,uz)
    
    derivs = np.array([dxdt, duxdt, dydt, duydt, dzdt, duzdt])
      
#     print(gamma1*Ez(y), Frad)
    
    return derivs


# In[5]:


#ολοκλήρωση

def oloklirosi(gamma0, t_end = 2*Delta, Dt = 10000):

    x, ux, y, uy, z, uz = [0.]*len(gamma0), [0.]*len(gamma0), [0.]*len(gamma0), [0.]*len(gamma0), [0.]*len(gamma0), [0.]*len(gamma0)

    x1, ux1, y1, uy1, z1, uz1 = [0.]*len(gamma0), [0.]*len(gamma0), [0.]*len(gamma0), [0.]*len(gamma0), [0.]*len(gamma0), [0.]*len(gamma0)

    #χρόνος της ολοκλήρωσης, σε μονάδες [qB0/mc]
    t = np.linspace(0.0, t_end, Dt)
    
    for i in range(0, len(gamma0)-1):
        
        global Frad
        
        Frad = [0.]
        
        #αρχικές συνθήκες
        state0 = np.array([0.0, 0.0, delta, -np.sqrt(gamma0[i]**2 - 1), 0.0, 0.0])
    
        #ολοκλήρωση τροχιάς με απώλειες
        state = odeint(speiser3D, state0, t, full_output=0)
    
        #ολοκλήρωση τροχιάς χωρίς απώλειες
        state1 = odeint(speiser3D_noloss, state0, t, full_output=0)
    
        x[i], ux[i], y[i], uy[i], z[i], uz[i] = state[:,0], state[:,1], state[:,2], state[:,3], state[:,4], state[:,5]  
    
        x1[i], ux1[i], y1[i], uy1[i], z1[i], uz1[i] = state1[:,0], state1[:,1], state1[:,2], state1[:,3], state1[:,4], state1[:,5] 
    
    return np.array([x, ux, y, uy, z, uz, x1, ux1, y1, uy1, z1, uz1])


# In[ ]:


get_ipython().run_cell_magic('script', 'false', "#διαγράμματα\nfig = plt.figure()\nfor i in range(0, len(x) - 1):\n    plt.plot(z[i]*(c/(e_charge*B_0/(e_mass*c)))/Delta,y[i], label = 'losses')\n    plt.plot(z1[i]*(c/(e_charge*B_0/(e_mass*c)))/Delta,y1[i], '--', label = 'no losses')\n# plt.plot(z[8],y[8], label = '$\\gamma_0 = 10$')\n# plt.plot(z[1],y[1], '--', label = '$\\gamma_0 = 100$')\n# plt.plot(z[2],y[2], label = '$\\gamma_0 = 1000$')\n# plt.plot(z[3],y[3], '-.', label = '$\\gamma_0 = 10000$')\nplt.xlabel('z/Delta')\nplt.ylabel('y')\nplt.legend(loc = 'best')\n# fig.savefig('troxia_yz_paroysiasi')")


# In[ ]:


get_ipython().run_cell_magic('script', 'false', "fig = plt.figure()\nplt.plot(z[0],x[0], ':', label = '$\\gamma_0 = 10$')\nplt.plot(z[1],x[1], '--', label = '$\\gamma_0 = 100$')\nplt.plot(z[2],x[2], label = '$\\gamma_0 = 1000$')\n# plt.plot(z[3],x[3], '-.', label = '$\\gamma_0 = 10000$')\nplt.xlabel('z')\nplt.ylabel('x')\n# fig.savefig('troxia_xz')")


# In[ ]:


get_ipython().run_cell_magic('script', 'false ', "fig = plt.figure()\n# for i in range(0, len(x) - 1):\n#     plt.plot(x[i],y[i], label = gamma0[i])\n# plt.plot(x[0],y[0], ':', label = '$\\gamma_0 = 10$')\n# plt.plot(x[1],y[1], '--', label = '$\\gamma_0 = 100$')\nplt.plot(x[2],y[2], label = '$\\gamma_0 = 1000$')\n# plt.plot(x[3],y[3], '-.', label = '$\\gamma_0 = 10000$')\nplt.xlabel('x')\nplt.ylabel('y')\n# fig.savefig('troxia_xy')")


# In[ ]:


get_ipython().run_cell_magic('script', 'false', "fig = plt.figure()\nfor i in range(0, len(x)-1):\n    plt.plot(z[i]*(c/(e_charge*B_0/(e_mass*c)))/Delta, gamma(ux[i], uy[i], uz[i]), label = 'losses')\n    plt.plot(z1[i]*(c/(e_charge*B_0/(e_mass*c)))/Delta, gamma(ux1[i], uy1[i], uz1[i]), '--', label = 'no losses')\n# plt.plot(z[8], gamma(ux[8], uy[8], uz[8]))\n# plt.plot(z[2], gamma(ux[2], uy[2], uz[2]), label = '$\\gamma_0 = 1000$')\n# plt.plot(z[3], gamma(ux[3], uy[3], uz[3]), '-.', label = '$\\gamma_0 = 10000$')\nplt.yscale('log')\nplt.xlabel('z/Delta')\nplt.ylabel('gamma')\nplt.legend(loc = 'best')\n# fig.savefig('gamma_paroyusiasi')")


# In[ ]:


get_ipython().run_cell_magic('script', 'false ', "fig = plt.figure()\nax = fig.add_subplot(111, projection='3d')\n\n# ax.plot(z[0], x[0], y[0])\nax.plot(z[1], x[1], y[1])\n# ax.plot(z[2], x[2], y[2])\n# ax.plot(z[3], x[3], y[3])\n\nax.set_xlabel('z')\nax.set_ylabel('x')\nax.set_zlabel('y')\n\nfig.savefig('troxia_xyz')")


# In[ ]:


# nu_crit = [nu_crit(ux[i], uy[i], uz[i], 1.) for i in range(0, len(x) - 1)]
# Lg = [Ploss(ux[i], uy[i], uz[i], Bx(y[i]), gamma0[i]) for i in range(0, len(x) - 1)]


# In[ ]:


get_ipython().run_cell_magic('script', 'false ', "fig = plt.figure()\nplt.xscale('log')\nplt.yscale('log')\nfor i in range(0, len(x)-1):\n    plt.plot(nu_crit(ux[i], uy[i], uz[i], 1.), Ploss(ux[i], uy[i], uz[i], 1.0))")


# In[ ]:


get_ipython().run_cell_magic('script', 'false', 'fig = plt.figure()\nfor i in range(0, len(x)-1):\n    plt.plot(t,Ploss(ux[i], uy[i], uz[i], 1.))')

