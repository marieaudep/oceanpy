"""Simple ocean physics functions"""
import numpy as np

# ############################################################
def melting_ice(Q,
                rho0 = 917.,
                L = 333.55e3,
                cpi = 2.05e3,
                dT = 20.):
    """Returns the volume of ice melted by heat energy Q

    Assuming a temperature change of dT has to be overcome

    Paramters
    =========
    Q : float
        supplied amount of heat energy (or flux) [J] ([W])
    rho0 : float
        reference density of ice [kg m-3]
    L : float
        latent heat of fusion [J kg-1]
    cpi : float
        specific heat capacity of ice [J kg-1 K-1]
    dT : float
        temperature change in ice [K] (e.g. 20 means ice has to be heated from -20 to 0 degC)

    Source
    ======
    Johnson et al. (2011, JGR 116.C1)
    """
    V = Q / (rho0 * (L + cpi*dT)) # m3
    return V


# ############################################################
def heating_seawater(Q,
                     V=1.,
                     rho0=1027.,
                     cpw=3986.):
    """Returns the change in temperature of water heated with Q

    A volume V [m3] of water is heated with energy (flux) Q [J]([W]). The water has reference density rho0 [kg m-3] and heat coefficient cpw [J kg-1 K-1]

    Parameters
    ==========
    V : float
        volume of water to be heated [m3]
    rho0 : float
        reference density [kg m-3]
    cpw : float
        specific heat capacity of water [J kg-1 K-1]

    Returns
    =======
    Change in temperature in K

    Source
    ======
    Johnson et al. (2011, JGR 116.C1)
    """
    dT = Q / (V*rho0*cpw)
    return dT


# ############################################################
def seawater_heat_content(V=1.,
                          rho=1027.,
                          dT=1.,
                          cpw=3986.):
    """Returns the energy required to heat sea water by dT degrees.

    Parameters
    ==========
    V : float
        volume [m3]
    dT : float
        change in temperature [K]
    rho : float
        density [kg m-3]
    cpw : float
        specific heat capacity of water [J kg-1 K-1]
        after Johnson et al. (2011, JGR 116.C1)

    Returns
    =======
    Q : heat energy needed to warm water [J]
    """
    Q = V*rho*cpw*dT
    return Q


# ############################################################
def ice_heat_content(V=1.,
                     rho0=917.,
                     L = 333.55e3,
                     cpi = 2.05e3,
                     dT = 20.):
    """Returns the amount of energy required to melt ice

    Assuming a temperature change of dT has to be overcome

    Paramters
    =========
    V : float
        volume to melt [m3]
    rho0 : float
        reference density of ice [kg m-3]
    L : float
        latent heat of fusion [J kg-1]
    cpi : float
        specific heat capacity of ice [J kg-1 K-1]
    dT : float
        temperature change in ice [K] (e.g. 20 means ice has to be heated from -20 to 0 degC)

    Source
    ======
    Johnson et al. (2011, JGR 116.C1)
    """
    Q = V * (rho0 * (L + cpi*dT)) # J
    return Q
