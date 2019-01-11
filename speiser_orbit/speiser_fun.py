#!/usr/bin/env python
# coding: utf-8

import numpy as np
from pulsars import c, e_charge, Pulsars

k = 3*10**2
pulsar = Pulsars(k)

#παράγοντας Lorentz
def gamma(ux,uy,uz):
    
    g = np.sqrt(1 + ux**2 + uy**2 + uz**2)
    
    return g

#################################
### Καρτεσιανές Συντεταγμένες ###
#################################

#μαγνητικό πεδίο στη x-διεύθυνση σε μονάδες [Β0 = Βlc]
def Bx(y, z, Delta, delta):
    
    if z <= Delta:
        bx = -np.tanh(y/delta)
    else:
        bx = -np.tanh(y/delta)
    
    return bx

#μαφνητικό πεδίο στην y-διευθυνση σε μονάδες [Β0 = Βlc]
def By(z, Delta):
    
    if z <= Delta:
        by = -1.0
    else:
        by = -0.0
    
    return by

#ηλεκτρικό πεδίο στη z-διευθυνση σε μονάδες [Β0 = Βlc]
def Ez(z, Delta):
    
    if z <= Delta:
        ez = 1.0
    else:
        ez = 0.0
    
    return ez

#συνιστώσες της δύναμης Lorentz
def Flor_x(y, z, ux, uy, uz, Delta, delta):
    
    fx = -By(z, Delta)*uz/gamma(ux, uy, uz)
    
    return fx

def Flor_y(y, z, ux, uy, uz, Delta, delta):
    
    fy = Bx(y, z, Delta, delta)*uz/gamma(ux, uy, uz)
    
    return fy
    
def Flor_z(y, z, ux, uy, uz, Delta, delta):
    
    fz = Ez(z, Delta) + (ux*By(z, Delta) - uy*Bx(y, z, Delta, delta))/gamma(ux, uy, uz)
    
    return fz




#γωνία εκπομπής
def theta(ux, uy, uz):
    
    th = np.arctan(uy/uz)
    
    return th

#κρίσιμη συχνότητα
def nu_crit(y, z, ux, uy, uz, Delta, delta):
    
    nu = np.abs(2.8*10E+6*Bx(y, z, Delta, delta)*gamma(ux, uy, uz)**2)
    
    return nu





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
def Ploss(ux, uy, uz, duxdt, duydt, duzdt, B_0):

    losses = (2*e_charge*gamma(ux, uy, uz)**4)/(3*B_0*Rc(ux, uy, uz, duxdt, duydt, duzdt)**2)
    # losses = (2*e_charge*gamma(ux, uy, uz)**4)/(3*B_0)/pulsar['crab']['rlc']**2
    return losses

#Radiation reaction force (αδιάστατο αλλά πρέπει να πολλαπλασιαστεί με ui σε κάθε εξίσωση)
def F_rad(ux, uy, uz, duxdt, duydt, duzdt, B_0):
    
    frad = Ploss(ux, uy, uz, duxdt, duydt, duzdt, B_0)*gamma(ux, uy, uz)/(gamma(ux, uy, uz)**2 - 1)
    
    return frad

