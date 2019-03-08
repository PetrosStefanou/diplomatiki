import numpy as np
import speiser_fun_cyl as sfc
from pulsars import h, c, e_charge
from scipy.special import kv
import scipy.integrate

def rc_nc_pr(r, phi, z_cyl, ur, uphi, uz_cyl, Rlc, Delta, delta_init, gamma0, B_0): 
    durdt = np.zeros(((len(gamma0) - 1), len(r[0])))
    duphidt = np.zeros(((len(gamma0) - 1), len(r[0])))
    duzdt = np.zeros(((len(gamma0) - 1), len(r[0])))

    r_curv = np.zeros(((len(gamma0) - 1), len(r[0])))
    nu_crit = np.zeros(((len(gamma0) - 1), len(r[0])))
    p_rad = np.zeros(((len(gamma0) - 1), len(r[0])))

    for i in range(len(gamma0) - 1):

        Frad = [0.]
        for j in range(len(r[i])):
            durdt[i][j] = sfc.Flor_r(r[i][j], phi[i][j], z_cyl[i][j], ur[i][j], uphi[i][j], uz_cyl[i][j], Rlc, Delta, delta_init) 
            + uphi[i][j]**2/(r[i][j]*sfc.gamma(ur[i][j], uphi[i][j], uz_cyl[i][j])) - Frad[-1]*ur[i][j]*0

            
            duphidt[i][j] = sfc.Flor_phi(r[i][j], phi[i][j], z_cyl[i][j], ur[i][j], uphi[i][j], uz_cyl[i][j], Rlc, Delta, delta_init) 
            - uphi[i][j]*ur[i][j]/(r[i][j]*sfc.gamma(ur[i][j], uphi[i][j], uz_cyl[i][j])) - Frad[-1]*ur[i][j]*0 - Frad[-1]*ur[i][j]*0

            
            duzdt[i][j] = sfc.Flor_z_cyl(r[i][j], phi[i][j], z_cyl[i][j], ur[i][j], uphi[i][j], uz_cyl[i][j], Rlc, Delta, delta_init)
            - Frad[-1]*ur[i][j]*0

            
            r_curv[i][j] = sfc.Rc(r[i][j], ur[i][j], uphi[i][j], uz_cyl[i][j], durdt[i][j], duphidt[i][j], duzdt[i][j])
            nu_crit[i][j] = sfc.nu_crit_curv(r[i][j], ur[i][j], uphi[i][j], uz_cyl[i][j], 
                                            durdt[i][j], duphidt[i][j], duzdt[i][j])
            p_rad[i][j] = sfc.Ploss_Rlc(r[i][j], ur[i][j], uphi[i][j], uz_cyl[i][j], 
                                        durdt[i][j], duphidt[i][j], duzdt[i][j], B_0)
            
            Frad.append(sfc.F_rad(r[i][j], ur[i][j], uphi[i][j], uz_cyl[i][j], durdt[i][j], duphidt[i][j], duzdt[i][j], 
                                B_0, sfc.Ploss))

    return (r_curv, nu_crit, p_rad)


def spectrum(r, z_cyl, nu_crit, p_rad, gamma0, Rlc, Delta, delta_init, t):

    # max_en = max(h*nu_crit[0]*6.25E+11)
    # min_en = min(h*nu_crit[0]*6.25E+11)
    # for i in range(0, len(gamma0)-1):
    #     if max(h*nu_crit[i]*6.25E+11) > max_en:
    #         max_en = max(h*nu_crit[i]*6.25E+11)
        # if min(h*nu_crit[i]*6.25E+11) < min_en:
        #     min_en = h*nu_crit[i]*6.25E+11
            
    # print('{:1.2E}, {:1.2E}'.format(max_en, min_en))

    min_en_pow = 6
    max_en_pow = 13
    bin_num = 201
    # en = [min_en]
    # while en[-1] < max_en:
    #     en.append(en[-1] + 0.1*en[-1])
    # print(len(en))
    # en = np.asarray(en)

    en = np.logspace(min_en_pow, max_en_pow, bin_num)
        
    ph_num = np.zeros((len(gamma0) - 1, len(en)))
    ph_en = np.zeros((len(gamma0) - 1, len(en)))
    ph_num_out = np.zeros((len(gamma0) - 1, len(en)))
    ph_en_out = np.zeros((len(gamma0) - 1, len(en)))
    ph_num_sep = np.zeros((len(gamma0) - 1, len(en)))
    ph_en_sep = np.zeros((len(gamma0) - 1, len(en)))
    ph_num_out2 = np.zeros((len(gamma0) - 1, len(en)))
    ph_en_out2 = np.zeros((len(gamma0) - 1, len(en)))
    # print('{:1.2E}, {:1.2E}'.format(max_en, min_en))

    delta1 = np.zeros(len(r[0]))
    for j in range(len(r[0])):
        delta1[j] = sfc.delta(r[0][j], Rlc, Delta, delta_init)
        if r[0][j] == Rlc + Delta:
            delta1[j] = np.NaN

    
    ph_num_tot = np.zeros(len(en))
    ph_en_tot = np.zeros(len(en))
    for i in range(0, len(gamma0)-1):
        for j in range(0, len(nu_crit[i])):
            for m in range(0, len(en) - 1):
                if h*nu_crit[i][j]*6.25E+11 >= en[m] and h*nu_crit[i][j]*6.25E+11 < en[m+1]:
                    if r[i][j] <= Rlc + Delta and r[i][j] >= Rlc:
                        ph_num[i][m] += 1
                        ph_en[i][m] += p_rad[i][j]/(en[m]-en[m-1])
                    elif r[i][j] >= Rlc + Delta:
                        if z_cyl[i][j] < delta1 and z_cyl>-delta1:
                            ph_num_out[i][m] += 1
                            ph_en_out[i][m] += p_rad[i][j]/(en[m]-en[m-1])
                        else:
                            ph_num_out2 += 1
                            ph_en_out2[i][m] += p_rad[i][j]/(en[m]-en[m-1])
                    elif z_cyl[i][j] > delta_init or z_cyl[i][j] < -delta_init:
                        ph_num_sep[i][m] += 1
                        ph_en_sep[i][m] += p_rad[i][j]/(en[m]-en[m-1])
                    ph_num_tot[m] += 1
                    ph_en_tot[m] += p_rad[i][j]/(en[m]-en[m-1])
    return (en, ph_num, ph_num_out, ph_num_out2, ph_num_sep, ph_en, ph_en_out, ph_en_out2, ph_en_sep, ph_num_tot, ph_en_tot)

def f(n, wc): 
    I = scipy.integrate.quad(lambda x: kv(5./3., x), n/wc, +np.inf)[0]
    return I

def sing_en_spec(nu, wc, g, rc):

    p = np.zeros(len(nu))
    for i in range(0, len(nu) - 1):
        p[i] = np.sqrt(3)*e_charge**2*g*nu[i]/(2*np.pi*rc*wc)*f(nu[i], wc)

    return p