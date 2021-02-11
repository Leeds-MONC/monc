"""
Functions to calculate physical properties
==============================

"""

import numpy as np
# from __future__ import print_function

def calcAirDensity(temperature, pressure):

    """
    Function to calculate air density from temperature and pressure
    ==============================

    """

        #### EXAMPLE OF USE:
        #### data = calcAirDensity(data['temperature'][:], data['pressure'][:])
        ####        temperature = K
        ####        pressure = hPa

    R = 2.8704  #### hPa kg-1 K-1

    print('Calculating air density profile:')
    print('')
    ### if temperature is 1D
    if np.ndim(temperature) == 1:
        rho = np.zeros([np.size(temperature)])
        for k in range(0,np.size(temperature)):
            rho[k] = pressure[k] / (R * temperature[k])
    ### if temperature is 2D
    elif np.ndim(temperature) == 2:
        rho = np.zeros([np.size(temperature,0), np.size(temperature,1)])
        for k in range(0,np.size(temperature, 1)):
            rho[:,k] = pressure[:,k] / (R * temperature[:,k])

    # print rho

    return rho

def calcThetaE(temperature, pressure, q):

    """
    Function to calculate equivalent potential temperature
    ==============================
    inputs:
    pressure = Pa
    temperature = K
    water vapour mixing ratio = kg/kg

    """

    L_vap = 2.555e6    # J/kg
    L_sub = 2.836e6  # J/kg
    cp = 1004.6      # J/kg.K
    cpd = 1005.7     # J/kg.K

    Rd = 287.04   # dry air J kg^-1 K^-1
    Rv = 461.50

    kd = Rd/cpd     # k dry air

    eps = Rd/Rv

    ### saturation vapour pressue
    evs = polysvp(temperature, 0)

    ### saturation mixing ratio and relative humidity
    qvs = (eps * evs) / (pressure - evs)
    rh = q / qvs
                #### gives RH as a fraction

    print('Calculating theta:')
    theta = temperature * np.power(1e5 / pressure, (Rd/cp))
    print('...')

    print('Calculating theta of dry air:')
    thetad = temperature * np.power(1e5 / (pressure - evs), kd)
    print('...')

    print('Calculating theta_e:')
    tempvar = (-1.0 * kd * q) / eps
    thetaE = thetad * np.power( rh, tempvar ) * np.exp(L_vap * q / (temperature * cpd) )         ###Bryan 2008

    print('...')
    print('Done!')

    return theta, thetaE

def calcThetaVL(temperature, pressure, q, ql, qi, tim, height):

    """
    Function to calculate virtual liquid potential temperature
    ==============================
    inputs:
    pressure = Pa
    temperature = K
    water vapour mixing ratio = kg/kg
    liquid mass mixing ratio = kg/kg
    ice mass mixing ratio = kg/kg

    Note that theta_l is based on 'liquid/frozen water static energy' (= cp.T + g.z - L.ql - Ls.qi)
        rather than potential temperature.
    Theta_vl is a conserved variable that is equal to virtual potential temperature (theta_v) in
        cloud-free air and so is used as a simplified measure of buoyancy

    """

    L_vap = 2.555e6    # J/kg
    L_sub = 2.836e6  # J/kg
    cp = 1004.6      # J/kg.K
    cpd = 1005.7     # J/kg.K

    Rd = 287.04   # dry air J kg^-1 K^-1
    Rv = 461.50

    g = 9.806       # m/s2

    kd = Rd/cpd     # k dry air

    eps = Rd/Rv

    ### total water mixing ratio
    qt = q + ql + qi

    ### saturation vapour pressue
    evs = polysvp(temperature, 0)

    ### saturation mixing ratio and relative humidity
    qvs = (eps * evs) / (pressure - evs)
    rh = q / qvs
                #### gives RH as a fraction

    print('Calculating theta:')
    theta = temperature * np.power(1e5 / pressure, (Rd/cp))
    print('...')

    print('Calculating theta_l:')
    theta_l = temperature - ((L_vap * ql)/cp) - ((L_sub * qi)/cp) + ((g * height)/cp)
    print('...')

    print('Calculating theta_vl:')
    cv = (1/eps) - 1
    theta_vl = theta_l + theta_l * cv * qt
    print('...')

    print('...')
    print('Done!')

    return theta, theta_l, theta_vl

def svp(T):

    """
    Function to calculate saturation vapour pressure
    ==============================
    inputs:
    temperature = K

    """

    tempC = T - 273.15

    satvappres = 6.112 * np.exp( 17.67*tempC / (tempC + 243.5) ) * 100

    return satvappres

def polysvp(t,type):

    """
    Function to calculate saturation vapour pressure
    ==============================
    inputs:
    temperature = K

    POLYSVP RETURNED IN UNITS OF PA.
    TYPE REFERS TO SATURATION WITH RESPECT TO LIQUID (0) OR ICE (1)
    COPIED FROM PARCEL MODEL
    """

    ### ! ICE
    if type == 1:
        dt = np.maximum(-80,t-273.16)
        svp = 6.11147274 + dt * (0.503160820 + dt * (0.188439774e-1 + dt * (0.420895665e-3 + dt *
            (0.615021634e-5 + dt * (0.602588177e-7 + dt * (0.385852041e-9 + dt * (0.146898966e-11 +
            0.252751365e-14 * dt)))))))
        svp = svp*100

    ### ! LIQUID
    elif type == 0:
        dt = np.maximum(-80,t-273.16)
        svp = 6.11239921 + dt * (0.443987641 + dt * (0.142986287e-1 + dt * (0.264847430e-3 + dt *
            (0.302950461e-5 + dt * (0.206739458e-7 + dt * (0.640689451e-10 + dt * (-0.952447341e-13 +
            -0.976195544e-15 * dt)))))))
        svp = svp*100


    return svp

def calcRH(temperature, pressure, q):

    """
    Function to calculate RH from given water vapour mixing ratio
    ==============================
    inputs:
    pressure = Pa
    temperature = K
    water vapour mixing ratio = kg/kg

    """

    L_vap = 2.555e6    # J/kg
    L_sub = 2.836e6  # J/kg
    cp = 1004.6      # J/kg.K
    cpd = 1005.7     # J/kg.K

    Rd = 287.04   # dry air J kg^-1 K^-1
    Rv = 461.50

    kd = Rd/cpd     # k dry air

    eps = Rd/Rv

    ### saturation vapour pressue
    evs = polysvp(temperature, 0)

    ### saturation mixing ratio and relative humidity
    qvs = (eps * evs) / (pressure - evs)
    rh = (q / qvs) * 100
                #### gives RH as a percentage

    return rh

def adiabatic_lwc(temperature, pressure):

    '''
    % [dlwcdz, dqldz, dqdp] = adiabatic_lwc(temperature, pressure)
    %
    %    Calculates the theoretical adiabatic rate of increase of LWC with
    %    height, or pressure, given the cloud base temperature and pressure
    %
    %    Prefers: temperature in K
    %             pressure in Pa
    %
    %    Returns: dlwc/dz in kg m-3 m-1
    %             dql/dz  in kg kg-1 m-1
    %             dql/dp  in kg kg-1 Pa-1
    %
    %    From Brenguier (1991)
    '''

    print('Calculating theoretical adiabatic rate of LWC increase with height:')

    # Define constants
    e  = 0.62198        # ratio of the molecular weight of water vapor to dry air
    g  = -9.81         # acceleration due to gravity (m s-1)
    cp = 1005          # heat capacity of air at const pressure (J kg-1 K-1)
    L  = 2.5e6         # latent heat of evaporation (J kg-1)
    R  = 461.5 * e     # specific gas constant for dry air (J kg-1 K-1)
    drylapse = -g / cp # dry lapse rate (K m-1)

    # This function requires temperature in Kelvin, pressure in Pa,
    # returns saturation vapour pressure (es) in Pa and
    # specific humidity mixing ratio (qs) in kg kg-1
    es = svp(temperature)
    qs=0.62198*es/(pressure-es)

    # calculate the density of the air (kg m-3)
    rhoa = pressure / (R * (1 + (0.6 * qs)) * temperature)

    # These equations expect temperature in K, pressure in
    # Pa, svp in Pa
    # Returns dql/dp in kg kg-1 Pa-1
    dqldp = -(1 - ((cp * temperature) / (L * e))) * (1 / ((cp * temperature / (L * e)) + (L * qs * rhoa / ((pressure - es))))) * (e * es) * (np.power(pressure - es, -2))

    # Returns dql/dz in kg kg-1 m-1
    dqldz = - (1 - (cp * temperature / (L * e))) * (1 / ((cp * temperature / (L * e)) + (L * qs * rhoa / (pressure - es)) )) * (rhoa * g * e * es) * (np.power(pressure - es, -2))

    # Returns dlwc/dz in kg m-3 m-1
    dlwcdz = rhoa * dqldz
    print('...')

    print('...')
    print('Done!')

    return dlwcdz, dqldz, dqldp

def calcTemperature(theta, pressure):

    """
    Function to calculate temperature from theta and pressure
    ==============================
    inputs:
    pressure = Pa
    potential temperature = K

    """

    L_vap = 2.555e6    # J/kg
    L_sub = 2.836e6  # J/kg
    cp = 1004.6      # J/kg.K
    cpd = 1005.7     # J/kg.K

    Rd = 287.04   # dry air J kg^-1 K^-1
    Rv = 461.50

    kd = Rd/cpd     # k dry air

    print('Calculating temperature:')
    temperature = theta / np.power(1e5 / pressure, (Rd/cp))
    print('...')

    print('...')
    print('Done!')

    return temperature
