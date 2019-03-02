import numpy as np
import speiser_fun as sf
import speiser_fun_cyl as sfc
import pulsars

#σύστημα διαφορικών με απώλειες ακτινοβολίας σε καρτεσιανές συντεταγμένες
def speiser_cart(state, t, Rlc, Delta, delta_init, B_0, Frad, q = +1):
    
    x, ux, y, uy, z, uz  = state
    
    #αδιάστατες εξισώσεις, οι πραγματικές έχουν διαιρεθεί με c*ω0 = mc/eB0c
    duxdt =  q*sf.Flor_x(y, z, ux, uy, uz, Rlc, Delta, delta_init) -Frad[-1]*ux*0
    duydt =  q*sf.Flor_y(y, z, ux, uy, uz, Rlc, Delta, delta_init) -Frad[-1]*uy*0
    duzdt =  q*sf.Flor_z(y, z, ux, uy, uz, Rlc, Delta, delta_init) -Frad[-1]*uz*0
    
    dxdt = ux/sf.gamma(ux,uy,uz)
    dydt = uy/sf.gamma(ux,uy,uz)
    dzdt = uz/sf.gamma(ux,uy,uz)
    
    Frad.append(sf.F_rad(ux, uy, uz, duxdt, duydt, duzdt, B_0))
    
    derivs = np.array([dxdt, duxdt, dydt, duydt, dzdt, duzdt])
    
#     print(Frad[-1]*ux, Frad[-1]*uy, Frad[-1]*uz, Ez(x))
    
    return derivs

    #σύστημα διαφορικών χωρίς απώλειες ακτινοβολίας σε καρτεσιανές συντεταγμένες
def speiser_cart_noloss(state, t, Rlc, Delta, delta, B_0, Frad, q = 1):
    
    x, ux, y, uy, z, uz  = state
    
    #αδιάστατες εξισώσεις, οι πραγματικές έχουν διαιρεθεί με c*ω0 = mc/eB0c
    duxdt =  sf.Flor_x(y, z, ux, uy, uz, Rlc, Delta, delta)
    duydt =  sf.Flor_y(y, z, ux, uy, uz, Rlc, Delta, delta)
    duzdt =  sf.Flor_z(y, z, ux, uy, uz, Rlc, Delta, delta)
    
    dxdt = ux/sf.gamma(ux,uy,uz)
    dydt = uy/sf.gamma(ux,uy,uz)
    dzdt = uz/sf.gamma(ux,uy,uz)
    
    derivs = np.array([dxdt, duxdt, dydt, duydt, dzdt, duzdt])
      
    # print(sf.gamma(ux, uy, uz), np.sqrt(ux**2 + uy**2 + uz**2), sf.Ez(z,Delta))
    
    return derivs

    #σύστημα διαφορικών με απώλειες ακτινοβολίας με ακτίνα καμπυλότητας υπολογισμένη από την τροχιά σε κυλινδρικές
def speiser_cyl(state, t, Rlc, Delta, delta_init, B_0, Frad, q, T, perc):
    
    r, ur, phi, uphi, z_cyl, uz_cyl = state
    
    #αδιάστατες εξισώσεις, οι πραγματικές έχουν διαιρεθεί με c*ω0 = mc/eB0c
    durdt =  q*sfc.Flor_r(r, phi, z_cyl, ur, uphi, uz_cyl, Rlc, Delta, delta_init) + uphi**2/(r*sfc.gamma(ur, uphi, uz_cyl)) - Frad[-1]*ur
    duphidt =  q*sfc.Flor_phi(r, phi, z_cyl, ur, uphi, uz_cyl, Rlc, Delta, delta_init) - ur*uphi/(r*sfc.gamma(ur, uphi, uz_cyl)) - Frad[-1]*uphi
    duzdt =  q*sfc.Flor_z_cyl(r, phi, z_cyl, ur, uphi, uz_cyl, Rlc, Delta, delta_init) - Frad[-1]*uz_cyl
    
    drdt = ur/sfc.gamma(ur,uphi,uz_cyl)
    dphidt = uphi/(r*sfc.gamma(ur,uphi,uz_cyl))
    dzdt = uz_cyl/sfc.gamma(ur,uphi,uz_cyl)
    
    Frad.append(sfc.F_rad(r, ur, uphi, uz_cyl, durdt, duphidt, duzdt, B_0, sfc.Ploss))
    # print(r, Rlc)
    derivs = np.array([drdt, durdt, dphidt, duphidt, dzdt, duzdt])
    
    
    # if int(t/(T*Delta)*100) != perc[-1]:
    #     perc.append(int(t/(T*Delta)*100)) 
    #     print('{} %'.format(perc[-1]))
    # print(t, T, Delta)
    # print('{:1.2E}, {:1.2E}, {:1.2E}, {:1.2E}, {:1.2E}'.format(sfc.Flor_r(r, phi, z_cyl, ur, uphi, uz_cyl, Rlc, Delta, delta), uphi**2/(r*sfc.gamma(ur, uphi, uz_cyl)), Frad[-1]*ur, ur, durdt))
    return derivs

#σύστημα διαφορικών με απώλειες ακτινοβολίας με ακτίνα καμπυλότητας Rlc σε κυλινδρικές
def speiser_cyl_Rlc(state, t, Rlc, Delta, delta, B_0, Frad, q = -1):
    
    r, ur, phi, uphi, z_cyl, uz_cyl = state
    
    #αδιάστατες εξισώσεις, οι πραγματικές έχουν διαιρεθεί με c*ω0 = mc/eB0c
    durdt =  q*sfc.Flor_r(r, phi, z_cyl, ur, uphi, uz_cyl, Rlc, Delta, delta) + uphi**2/(r*sfc.gamma(ur, uphi, uz_cyl)) - Frad[-1]*ur 
    duphidt =  q*sfc.Flor_phi(r, phi, z_cyl, ur, uphi, uz_cyl, Rlc, Delta, delta) - ur*uphi/(r*sfc.gamma(ur, uphi, uz_cyl)) - Frad[-1]*uphi 
    duzdt =  q*sfc.Flor_z_cyl(r, phi, z_cyl, ur, uphi, uz_cyl, Rlc, Delta, delta) - Frad[-1]*uz_cyl
    
    drdt = ur/sfc.gamma(ur,uphi,uz_cyl)
    dphidt = uphi/(r*sfc.gamma(ur,uphi,uz_cyl))
    dzdt = uz_cyl/sfc.gamma(ur,uphi,uz_cyl)
    
    Frad.append(sfc.F_rad(r, ur, uphi, uz_cyl, durdt, duphidt, duzdt, B_0, sfc.Ploss_Rlc))
    
    derivs = np.array([drdt, durdt, dphidt, duphidt, dzdt, duzdt])

    return derivs

    #σύστημα διαφορικών χωρίς απώλειες ακτινοβολίας σε κυλινδρικές
def speiser_cyl_noloss(state, t, Rlc, Delta, delta, B_0, Frad, q = -1):
    
    r, ur, phi, uphi, z_cyl, uz_cyl = state
    
    #αδιάστατες εξισώσεις, οι πραγματικές έχουν διαιρεθεί με c*ω0 = mc/eB0c
    durdt =  q*sfc.Flor_r(r, phi, z_cyl, ur, uphi, uz_cyl, Rlc, Delta, delta) + uphi**2/(r*sfc.gamma(ur, uphi, uz_cyl)) 
    duphidt =  q*sfc.Flor_phi(r, phi, z_cyl, ur, uphi, uz_cyl, Rlc, Delta, delta) - ur*uphi/(r*sfc.gamma(ur, uphi, uz_cyl))
    duzdt =  q*sfc.Flor_z_cyl(r, phi, z_cyl, ur, uphi, uz_cyl, Rlc, Delta, delta)
    
    drdt = ur/sfc.gamma(ur,uphi,uz_cyl)
    dphidt = uphi/(r*sfc.gamma(ur,uphi,uz_cyl))
    dzdt = uz_cyl/sfc.gamma(ur,uphi,uz_cyl)
    
    derivs = np.array([drdt, durdt, dphidt, duphidt, dzdt, duzdt])

    return derivs