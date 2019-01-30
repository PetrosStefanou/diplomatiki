import numpy as np
from pulsars import c, e_charge, Pulsars


def gamma(u1,u2,u3):
    
    g = np.sqrt(1 + u1**2 + u2**2 + u3**2)
    
    return g

#################################
### Κυλινδρικές Συντεταγμένες ###
#################################

#μαγνητικό πεδίο στη φ-διεύθυνση σε μονάδες [Β0 = Βlc]
def Bphi(r, z_cyl, Rlc, Delta, delta):
    
    if r <= Rlc and r >= Rlc - Delta:
        bphi = -np.tanh(z_cyl/delta)
    elif r < Rlc - Delta:
        bphi = 0.
    else:
        bphi = -np.tanh(z_cyl/delta)
    
    return bphi

#μαφνητικό πεδίο στην z-διευθυνση σε μονάδες [Β0 = Βlc]
def Bz_cyl(r, Rlc, Delta):
    
    if r <= Rlc and r >= Rlc - Delta:
        bz = -1.0
    elif r < Rlc - Delta:
        bz = -10.0
    else:
        bz = 0.0
    
    return bz

#ηλεκτρικό πεδίο στη r-διευθυνση σε μονάδες [Β0 = Βlc] 
def Er(r, Rlc, Delta):
    
    if r <= Rlc and r >= Rlc - Delta:
        er = 1.0
    else:
        er = 0.0
    
    return er

def Flor_r(r, phi, z_cyl, ur, uphi, uz_cyl, Rlc, Delta, delta):
    
    fr = Er(r, Rlc, Delta) + (uphi*Bz_cyl(r, Rlc, Delta) - uz_cyl*Bphi(r, z_cyl, Rlc, Delta, delta))/gamma(ur, uphi, uz_cyl) 
    
    return fr

def Flor_phi(r, phi, z_cyl, ur, uphi, uz_cyl, Rlc, Delta, delta):
    
    fphi = -Bz_cyl(r, Rlc, Delta)*ur/gamma(ur, uphi, uz_cyl) 
    
    return fphi

def Flor_z_cyl(r, phi, z_cyl, ur, uphi, uz_cyl, Rlc, Delta, delta):
    
    fz = Bphi(r, z_cyl, Rlc, Delta, delta)*ur/gamma(ur, uphi, uz_cyl)
    
    return fz

#μέτρο του dudt
def dudt(r, ur, uphi, uz_cyl, durdt, duphidt, duzdt):
    
    metro_dudt = ((durdt - uphi**2/(r*gamma(ur, uphi, uz_cyl)))**2 + (duphidt + ur*uphi/(r*gamma(ur, uphi, uz_cyl)))**2 + duzdt**2)**0.5
    
    return metro_dudt

#μέτρο του u
def u(u1, u2, u3):
    
    metro_u = (u1**2 + u2**2 + u3**2)**0.5
    
    return metro_u

#μέτρο του εξωτερικού γινομένου ταχύτητας-επιτάχυνσης (σε πραγματικες μονάδες [m**2/s**3])
def u_x_a(r, ur, uphi, uz_cyl, durdt, duphidt, duzdt):
    
    ar = (durdt - uphi**2/(r*gamma(ur, uphi, uz_cyl)))/gamma(ur, uphi, uz_cyl) - ur*u(ur, uphi, uz_cyl)*dudt(r, ur, uphi, uz_cyl, durdt, duphidt, duzdt)/(c**2*gamma(ur, uphi, uz_cyl)**3)
    aphi = (duphidt + ur*uphi/(r*gamma(ur, uphi, uz_cyl)))/gamma(ur, uphi, uz_cyl) - uphi*u(ur, uphi, uz_cyl)*dudt(r, ur, uphi, uz_cyl, durdt, duphidt, duzdt)/(c**2*gamma(ur, uphi, uz_cyl)**3)
    az = duzdt/gamma(ur, uphi, uz_cyl) - uz_cyl*u(ur, uphi, uz_cyl)*dudt(r, ur, uphi, uz_cyl, durdt, duphidt, duzdt)/(c**2*gamma(ur, uphi, uz_cyl)**3)
    
    cross_product = ((uphi*az - uz_cyl*aphi)**2 + (ur*az - uz_cyl*ar)**2 + (ur*aphi - uphi*ar)**2)**0.5
    
    return cross_product

#ακτίνα καμπυλότητας (Rc = v**3/|vxa| = (u/gamma)**3/(|uxa|/gamma), σε πραγματικές μονάδες [m])
def Rc(r, ur, uphi, uz_cyl, durdt, duphidt, duzdt):
    
    rcurv = (u(ur, uphi, uz_cyl)/gamma(ur, uphi, uz_cyl))**3/(u_x_a(r, ur, uphi, uz_cyl, durdt, duphidt, duzdt)/gamma(ur, uphi, uz_cyl))
    
    return rcurv

#Ενεργειακές απώλειες λόγω ακτινοβολίας με ακτίνα καμπυλότητας υπολογισμένη από την τροχιά (αδιάστατο, έχει διαιρεθεί με ecB0)
def Ploss(r, ur, uphi, uz_cyl, durdt, duphidt, duzdt, B_0):

    losses = (2*e_charge*gamma(ur, uphi, uz_cyl)**4)/(3*B_0*Rc(r, ur, uphi, uz_cyl, durdt, duphidt, duzdt)**2)
    # losses = (2*e_charge*gamma(ur, uphi, uz_cyl)**4)/(3*B_0)/pulsar['crab']['rlc']**2
    return losses

#Ενεργειακές απώλειες λόγω ακτινοβολίας με ακτίνα καμπυλότητας την Rlc (αδιάστατο, έχει διαιρεθεί με ecB0)
def Ploss_Rlc(r, ur, uphi, uz_cyl, durdt, duphidt, duzdt, B_0):

    k = 3*10**2
    pulsar = Pulsars(k)
    
    losses = (2*e_charge*gamma(ur, uphi, uz_cyl)**4)/(3*B_0*pulsar['crab']['rlc']**2)

    return losses



#Radiation reaction force (αδιάστατο αλλά πρέπει να πολλαπλασιαστεί με ui σε κάθε εξίσωση)
def F_rad(r, ur, uphi, uz_cyl, durdt, duphidt, duzdt, B_0, rad_losses):
    
    frad = rad_losses(r, ur, uphi, uz_cyl, durdt, duphidt, duzdt, B_0)*gamma(ur, uphi, uz_cyl)/(gamma(ur, uphi, uz_cyl)**2 - 1)
    
    return frad

def nu_crit_curv(r, ur, uphi, uz_cyl, durdt, duphidt, duzdt):
    
    nu = (3*c/2)*gamma(ur, uphi, uz_cyl)**3/Rc(r, ur, uphi, uz_cyl, durdt, duphidt, duzdt)

    return nu
