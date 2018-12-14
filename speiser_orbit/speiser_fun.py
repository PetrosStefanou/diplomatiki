#!/usr/bin/env python
# coding: utf-8

# In[1]:


#παράγοντας Lorentz
def gamma(ux,uy,uz):
    
    g = np.sqrt(1 + ux**2 + uy**2 + uz**2)
    
    return g

#μαγνητικό πεδίο σε μονάδες Β0 = Βlc
def Bx(y, z):
    if z <= Delta:
        if y < -delta:
            bx = 1.
        elif y >= -delta and y <= delta:
            bx = -np.tanh(y/delta)
        else:
            bx = -1.
    else:
        bx = 0.
    return bx

#ηλεκτρικό πεδίο σε μονάδες Β0 = Βlc
def Ez(z):
    
    if z <= Delta:
        ez = 1.0
    else:
        ez = 30.0
    
    return ez

def By(z):
    if z <= Delta:
        by = 1.
    else:
        by = 30.
        
    return by
        
#γωνία εκπομπής
def theta(ux, uy, uz):
    
    th = np.arctan(uy/yz)
    
    return th

#κρίσιμη συχνότητα
def nu_crit(ux, uy, uz, y):
    
    nu = np.abs(2.8*10E+6*Bx(y)*gamma(ux, uy, uz)**2)
    
    return nu


# In[2]:


#συνιστώσες της δύναμης Lorentz
def Flor_x(y, z, ux, uy, uz):
    
    fx = -By(z)*uz/gamma(ux, uy, uz)
    
    return fx

def Flor_y(y, z, ux, uy, uz):
    
    fy = Bx(y, z)*uz/gamma(ux, uy, uz)
    
    return fy

def Flor_z(y, z, ux, uy, uz):
    
    fz = Ez(z) + (ux*By(z) - uy*Bx(y, z))/gamma(ux, uy, uz)
    
    return fz


# In[3]:


#μέτρο του dudt
def dudt(duxdt, duydt, duzdt):
    
    metro_dudt = (duxdt**2 + duydt**2 + duzdt**2)**0.5
    
    return metro_dudt

#μέτρο του u
def u(ux, uy, uz):
    
    metro_u = (ux**2 + uy**2 + uz**2)**0.5
    
    return metro_u

#μέτρο του εξωτερικού γινομένου ταχύτητας-επιτάχυνσης (σε πραγματικες μονάδες [m**2/s**3])
def u_x_a(ux, uy, uz, duxdt, duydt, duzdt):
    
    ax = duxdt/gamma(ux, uy, uz) - ux*u(ux, uy, uz)*dudt(duxdt, duydt, duzdt)/(c**2*gamma(ux, uy, uz)**3)
    ay = duydt/gamma(ux, uy, uz) - uy*u(ux, uy, uz)*dudt(duxdt, duydt, duzdt)/(c**2*gamma(ux, uy, uz)**3)
    az = duzdt/gamma(ux, uy, uz) - uz*u(ux, uy, uz)*dudt(duxdt, duydt, duzdt)/(c**2*gamma(ux, uy, uz)**3)
    
    cross_product = ((uy*az - uz*ay)**2 + (ux*az - uz*ax)**2 + (ux*ay - uy*ax)**2)**0.5
    
    return cross_product

#ακτίνα καμπυλότητας (Rc = v**3/|vxa| = (u/gamma)**3/(|uxa|/gamma), σε πραγματικές μονάδες [m])
def Rc(ux, uy, uz, duxdt, duydt, duzdt):
    
    rcurv = (u(ux, uy, uz)/gamma(ux, uy, uz))**3/(u_x_a(ux, uy, uz, duxdt, duydt, duzdt)/gamma(ux, uy, uz))
    
    return rcurv

#Ενεργειακές απώλειες λόγω ακτινοβολίας (αδιάστατο, έχει διαιρεθεί με ecB0)
def Ploss(ux, uy, uz, duxdt, duydt, duzdt):

    losses = (2*e_charge*gamma(ux, uy, uz)**4)/(3*B_0*Rc(ux, uy, uz, duxdt, duydt, duzdt)**2)

    return losses

#Radiation reaction force (αδιάστατο αλλά πρέπει να πολλαπλασιαστεί με ui σε κάθε εξίσωση)
def F_rad(ux, uy, uz, duxdt, duydt, duzdt):
    
    frad = Ploss(ux, uy, uz, duxdt, duydt, duzdt)*gamma(ux, uy, uz)/(gamma(ux, uy, uz)**2 - 1)
    
    return frad


# In[ ]:




