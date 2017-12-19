# Oliver Angelil 2017
# translation of Paul O'Gorman's matlab code:
# http://www.mit.edu/~pog/src/precip_extremes_scaling.m
# following O'Gorman and Schneider, PNAS, 106, 14773-14777, 2009

# import python modules
import numpy as np

# INPUTS                                                                                                                                                    
# omega: vertical profile of vertical velocity in Pa/s (conditioned on extremes occuring)                                                                   
# temp: vertical profile of temperature in K (conditioned on extremes occuring)                                                                             
# plev: pressure in Pa                                                                                                                                      
# ps: surface pressure in Pa (conditioned on extremes occuring)                                                                                             

# OUTPUT                                                                                                                                                    
# precip: precipitation estimate in kg/m^2/s (note doesn't include precipitation efficiency factor)

# ==================================================================
# FIRST DEFINE FIVE FUNCTIONS: pars, saturation_thermodynamics, sat_deriv, 
# moist_adiabatic_lapse_rate, integrate
# ==================================================================
def pars(parstr):
    if parstr=='gravity':
        par = 9.80665  # in m / s^2
    elif parstr=='gas_constant':
        par = 287.04  # J / kg / K
    elif parstr=='kappa':
        par = 2.0/7.0
    elif parstr=='cp':
        par = pars('gas_constant')/pars('kappa')  # specific heat at constant pressure
    elif parstr=='cp_v':
        par = 1870  # specific heat water vapor [J/kg/K]
    elif parstr=='gas_constant_v':
        par = 461.50  # gas constant water vapor [J/kg/K]
    elif parstr=='latent_heat_v':
        par = 2.5e6  # latent heat of evaporation [J/kg]

    return par



def saturation_thermodynamics(temp, plev, calc_type='simple'):
    # the effective saturation vapor pressure, specific humidity, mixing
    # ratio, and latent heat of vaporization/sublimation. Depending on
    # calc_type, is uses different formulations:
    #    'era': Formulation used by ECMWF based on the modified
    #           Tetens formula.
    #    'simple': constant latent heat of condensation


    Rd = pars('gas_constant')
    Rv = pars('gas_constant_v')
    gc_ratio = Rd/Rv  # ratio of gas constants for dry air and water vapor

    if calc_type=='era':

        # ECMWF formulation described by Simmons et al. (1999: QJRMS, 125,
        # 353--386), which uses saturation over ice for temperatures
        # less than 250 K and a quadratic interpolation between
        # saturation over ice and over liquid water for temperatures
        # between 250 K and 273 K

        # coefficients in Tetens saturation vapor pressure es = es0 * exp(a3 * (T-T0)/(T-a4))
        es0       = 611.21  # saturation vapor pressure at T0 (Pa)
        T0        = 273.16  # (K)
        Ti        = T0 - 23 # (K)

        a3l       = 17.502  # liquid water (Buck 1981)
        a4l       = 32.19   # (K)
         
        a3i       = 22.587  # ice (Alduchov and Eskridge 1996)
        a4i       = -0.7    # (K)

        # saturation vapor pressure over liquid and ice
        esl       = es0 * np.exp(a3l * (temp - T0)/(temp - a4l))
        esi       = es0 * np.exp(a3i * (temp - T0)/(temp - a4i))

        # latent heat of sublimation
        Ls0       = 2.834e6  # latent heat of sublimation  (J / kg) [+- 0.01 error for 173 K < T < 273 K]
        Ls        = Ls0 * np.ones(temp.size)

        # latent heat of vaporization
        Lv0       = 2.501e6  # latent heat of vaporization at triple point (J / kg)
        cpl       = 4190     # heat capacity of liquid water (J / kg / K)
        cpv       = pars('cp_v')     # heat capacity of water vapor (J / kg / K)
        Lv        = Lv0 - (cpl - cpv) * (temp - T0)

        # compute saturation vapor pressure and latent heat over liquid/ice mixture
        iice      = temp <= Ti
        iliquid   = temp >= T0
        imixed    = (temp > Ti) * (temp < T0)
        
        es        = np.ones(temp.size)*np.NaN
        L         = np.ones(temp.size)*np.NaN

        if any(iice):
            es[iice]  = esi[iice]
            L[iice]   = Ls[iice]

        if any(iliquid):
            es[iliquid] = esl[iliquid]
            L[iliquid]  = Lv[iliquid]
 
        if any(imixed):
            a         = ((temp[imixed] - Ti)/(T0 - Ti))**2;
            es[imixed]= (1-a) * esi[imixed] + a * esl[imixed]
            L[imixed] = (1-a) * Ls[imixed] + a * Lv[imixed]

    else:
        # simple calculation assuming constant latent heat of vaporization 
        # and only one vapor-liquid phase transition
        T0        = 273.16
        es0       = 610.78
        L         = pars('latent_heat_v')

        # saturation vapor pressure
        es        = es0 * np.exp(L/Rv*(1.0/T0 - 1.0/temp))

        
    # saturation mixing ratio
    rs          = gc_ratio * es / (plev - es)
        
    # saturation specific humidity
    qs          = rs / (1 + rs)

    return es, qs, rs, L



def sat_deriv(plev, temp):
    # Calculates derivatives of the saturation specific humidity wrt
    # temperature and pressure

    # use finite difference approximation
    dp = 0.1
    dT = 0.01

    es_p_plus, qs_p_plus,_,_ = saturation_thermodynamics(temp, plev+dp, 'era')
    es_p_minus, qs_p_minus,_,_ = saturation_thermodynamics(temp, plev-dp, 'era')
    es_T_plus, qs_T_plus,_,_ = saturation_thermodynamics(temp+dT, plev, 'era')
    es_T_minus, qs_T_minus,_,_ = saturation_thermodynamics(temp-dT, plev, 'era')

    dqsat_dp   = (qs_p_plus-qs_p_minus)/(2.0*dp)
    
    dqsat_dT   = (qs_T_plus-qs_T_minus)/(2.0*dT)
    
    dln_esat_dT = (np.log(es_T_plus)-np.log(es_T_minus))/2/dT

    return dqsat_dp, dqsat_dT, dln_esat_dT



def moist_adiabatic_lapse_rate(temp, plev, calc_type):
    # MOIST_ADIABATIC_LAPSE_RATE Returns saturated moist-adiabatic lapse rate.
    #
    # Units are K / m.

    g               = pars('gravity')
    cpd             = pars('cp')
    cpv             = pars('cp_v')
    Rd              = pars('gas_constant')
    Rv              = pars('gas_constant_v')
    gc_ratio        = Rd / Rv

    es, qs, rs, L = saturation_thermodynamics(temp, plev, calc_type);

    # cf. Holton (p. 503), obtained
    # by using rs<<1 in the more accurate expression below
    if calc_type=='simple':
        lapse_rate = g/cpd * (1 + L*rs / Rd / temp) / (1 + L**2 * rs /(cpd * Rv * temp**2))
    else:
        # cf. Emanuel (p. 131)
        lapse_rate      = g/cpd * (1 + rs) / (1 + cpv/cpd*rs) * (1 + L*rs / Rd / temp) / (1 + L**2 * rs * (1 + rs/gc_ratio)/(Rv * temp**2 * (cpd + rs*cpv)))

    return lapse_rate



def integrate(f, x):
    #INTEGRATE  Computes one-dimensional integral.
    #    INTEGRATE(f, x) computes an approximation to the integral
    #    \int f(x) dx over the range of the input vector x. The input x
    #    must be a vector. The input f can be a column vector or a
    #    matrix. The number of rows of f must be the same as the length
    #    of x. If f is a matrix, a vertical integral is computed for
    #    each column.

    dx1 = np.gradient(x)
    dx1[0] = 0.5 * dx1[1]
    dx1[-1] = 0.5 * dx1[-1]
    dx2 = np.tile(dx1, (1, 1)) # why is this line needed? When would f be a matrix?? Ask Paul O'Gorman. 
    F = np.sum(f * dx2)

    return F


# ==================================================================
# PERFORM FINAL CALCULATION USING THE FIVE FUNCTIONS
# ==================================================================

def scaling(omega, temp, plev, ps):

    # criterion for identifying tropopause
    crit_lapse_rate = 0.002 # (k/m) for tropopause
    plev_mask = 0.05e5 # (Pa) exclude levels above this as a fail-safe

    dqsat_dp, dqsat_dT,_ = sat_deriv(plev, temp)
    es, qsat, rsat, latent_heat = saturation_thermodynamics(temp, plev, 'era')
    lapse_rate = moist_adiabatic_lapse_rate(temp, plev, 'era')

    # virtual temperature
    temp_virtual = temp*(1.0+qsat*(pars('gas_constant_v')/pars('gas_constant')-1.0))

    # density
    rho = plev/pars('gas_constant')/temp_virtual

    dT_dp = lapse_rate/pars('gravity')/rho

    # find derivative of saturation specific humidity with respect to pressure along 
    # a moist adiabat at the given temperature and pressure for each level
    dqsat_dp_total = dqsat_dp+dqsat_dT*dT_dp

    # mask above tropopause using simple lapse rate criterion
    dT_dp_env = np.gradient(temp, plev)
    lapse_rate_env = dT_dp_env*rho*pars('gravity')

    itrop = np.where(lapse_rate_env>crit_lapse_rate)[0]
    if itrop.size!=0:
        if np.max(itrop)+1<len(plev):
            dqsat_dp_total[np.max(itrop)+1:]=0

    # mask above certain level as fail safe
    dqsat_dp_total[plev<plev_mask]=0

    dqsat_dp_total_omega = dqsat_dp_total*omega

    # replaces nans with zeros as subsurface values should not contribute
    # to the column integral
    dqsat_dp_total_omega[np.isnan(dqsat_dp_total_omega)]=0

    # also use surface pressure to zero subsurface values
    kbot = plev>ps
    if any(kbot):
        dqsat_dp_total_omega[kbot]=0

    # integrate in the vertical
    precip = -integrate(-dqsat_dp_total_omega,plev)/pars('gravity')

    return precip













