import numpy as np
import speiser_fun_cyl as sfc
from pulsars import h, c, e_charge
from scipy.special import kv
import scipy.integrate

def rc_nc_pr(r, phi, z_cyl, ur, uphi, uz_cyl, Rlc, Delta, delta, gamma0, B_0): 
    durdt = np.zeros(((len(gamma0) - 1), len(r[0])))
    duphidt = np.zeros(((len(gamma0) - 1), len(r[0])))
    duzdt = np.zeros(((len(gamma0) - 1), len(r[0])))

    r_curv = np.zeros(((len(gamma0) - 1), len(r[0])))
    nu_crit = np.zeros(((len(gamma0) - 1), len(r[0])))
    p_rad = np.zeros(((len(gamma0) - 1), len(r[0])))

    for i in range(len(gamma0) - 1):
        # durdt[i] = np.zeros(len(r[i]))
        # duphidt[i] = np.zeros(len(r[i]))
        # duzdt[i] = np.zeros(len(r[i]))
        
        # r_curv[i] = np.zeros(len(r[i]))
        # nu_crit[i] = np.zeros(len(r[i]))
        # p_rad[i] = np.zeros(len(r[i]))
        
        Frad = [0.]
        for j in range(len(r[i])):
            durdt[i][j] = sfc.Flor_r(r[i][j], phi[i][j], z_cyl[i][j], ur[i][j], uphi[i][j], uz_cyl[i][j], Rlc, Delta, delta) 
            + uphi[i][j]**2/(r[i][j]*sfc.gamma(ur[i][j], uphi[i][j], uz_cyl[i][j])) - Frad[-1]*ur[i][j]*0

            
            duphidt[i][j] = sfc.Flor_phi(r[i][j], phi[i][j], z_cyl[i][j], ur[i][j], uphi[i][j], uz_cyl[i][j], Rlc, Delta, delta) 
            - uphi[i][j]*ur[i][j]/(r[i][j]*sfc.gamma(ur[i][j], uphi[i][j], uz_cyl[i][j])) - Frad[-1]*ur[i][j]*0 - Frad[-1]*ur[i][j]*0

            
            duzdt[i][j] = sfc.Flor_z_cyl(r[i][j], phi[i][j], z_cyl[i][j], ur[i][j], uphi[i][j], uz_cyl[i][j], Rlc, Delta, delta)
            - Frad[-1]*ur[i][j]*0

            
            r_curv[i][j] = sfc.Rc(r[i][j], ur[i][j], uphi[i][j], uz_cyl[i][j], durdt[i][j], duphidt[i][j], duzdt[i][j])
            nu_crit[i][j] = sfc.nu_crit_curv(r[i][j], ur[i][j], uphi[i][j], uz_cyl[i][j], 
                                            durdt[i][j], duphidt[i][j], duzdt[i][j])
            p_rad[i][j] = sfc.Ploss_Rlc(r[i][j], ur[i][j], uphi[i][j], uz_cyl[i][j], 
                                        durdt[i][j], duphidt[i][j], duzdt[i][j], B_0)
            
            Frad.append(sfc.F_rad(r[i][j], ur[i][j], uphi[i][j], uz_cyl[i][j], durdt[i][j], duphidt[i][j], duzdt[i][j], 
                                B_0, sfc.Ploss))

    return (r_curv, nu_crit, p_rad)


def spectrum(nu_crit, p_rad, gamma0, t):

    max_en = max(h*nu_crit[0]*6.25E+11)
    min_en = min(h*nu_crit[0]*6.25E+11)
    for i in range(0, len(gamma0)-1):
        if max(h*nu_crit[i]*6.25E+11) > max_en:
            max_en = max(h*nu_crit[i]*6.25E+11)
        if min(h*nu_crit[i]*6.25E+11) < min_en:
            min_en = h*nu_crit[i]*6.25E+11
            
    # print('{:1.2E}, {:1.2E}'.format(max_en, min_en))

    # en = np.linspace(min_en, max_en, 2*10**2)
    en = [10**6]
    while en[-1] < max_en:
        en.append(en[-1] + 0.1*en[-1])
    print(len(en))
    np.asarray(en)
    print(type(en))
        
    ph_num = np.zeros((len(gamma0) - 1, len(en)))
    ph_en = np.zeros((len(gamma0) - 1, len(en)))
    # print('{:1.2E}, {:1.2E}'.format(max_en, min_en))

    
    ph_num_tot = np.zeros(len(en))
    ph_en_tot = np.zeros(len(en))
    for i in range(0, len(gamma0)-1):
        # ph_num[i] = np.zeros(len(en))
        # ph_en[i] = np.zeros(len(en))
        
        
    #     print(len(ph_num[i]))
        for j in range(0, len(nu_crit[i])):
            for w in range(0, len(en)):
                if h*nu_crit[i][j]*6.25E+11 >= en[w] and h*nu_crit[i][j]*6.25E+11 < en[w+1]:
                    ph_num[i][w] += 1
                    ph_num_tot[w] += 1
                    ph_en[i][w] += p_rad[i][j]
                    ph_en_tot[w] += p_rad[i][j]
                    
    

    return (en, ph_num, ph_en, ph_num_tot, ph_en_tot)

def f(n, wc): 
    I = scipy.integrate.quad(lambda x: kv(5./3., x), n/wc, +np.inf)[0]
    return I

def sing_en_spec(nu, wc, g, rc):

    p = np.zeros(len(nu))
    for i in range(0, len(nu) - 1):
        p[i] = np.sqrt(3)*e_charge**2*g*nu[i]/(2*np.pi*rc*wc)*f(nu[i], wc)

    return p