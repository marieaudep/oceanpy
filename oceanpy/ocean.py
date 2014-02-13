"""Formulae for length scales and coefficients in Physical Oceanography"""
import numpy as np
from . import stats


def beta_coefficient(lat):
    """Compute linear coefficient of variation of the coriolis parameter for latitude `lat` (deg)"""
    omega = 7.2921e-5 # rad s-1
    R = 6.371e6 # m
    return 2. * omega * np.cos(np.radians(lat)) / R


def coriolis_parameter(lat0):
    """Compute the coriolis parameter for latitude `lat` (deg)"""
    omega = 7.2921e-5 # rad s-1
    return 2. * omega * np.sin(np.radians(lat0))


def buoyancy_frequency(rho,drhodz):
    """Compute the Brunt-V\"ais\"al\"a frequency for density *rho* and vertical density gradient *drhodz* (SI units)"""
    g = -9.80665 # m s-2    
    return np.sqrt(np.abs(-g / rho * drhodz))


def baroclinic_rossby_radius(lat0,rho,drhodz,H=5000):
    """Compute the baroclinic Rossby radius

    Parameters
    ----------
    lat0 : float
        Reference latitude [deg]
    H : float
        Scale height [m] (default: 5000 m)
    rho : float
        Density [kg m-3]
    drhodz : float
        Vertical density gradient [kg m-3 m-1]
    """
    return np.abs(buoyancy_frequency(rho,drhodz) * H  / coriolis_parameter(lat0))


def baroclinic_rossby_radius_chelton(lat0,rho,drhodz,dz):
    """Compute the first baroclinic Rossby radius after Chelton et al. 1998 J.Phys.Oc."""
    # compute phase speed, equation 2.2
    c1 = 1/np.pi * np.sum(buoyancy_frequency(rho,drhodz)*dz)

    # compute radius
    if np.abs(lat0) >= 5:
        # equation 2.3a
        return c1 / np.abs(coriolis_parameter(lat0))
    else:
        # equation 2.3b
        return np.sqrt(c1/(2*beta_coefficient(lat0)))


def eddy_length_scale(L_Rossby,scale='global'):
    """Eddy length scale in the ocean after Stammer 1997, J.Phys.Oc

    Parameters
    ----------
    L_Rossby : float
        Rossby radius of deformation [m]
    scale : str in ['global','Atlantic']
        Which set of coefficients to use
    """
    if scale == 'global':
        return 2.*(0.8 * L_Rossby + 88e3)
    elif scale == 'Atlantic':
        return 2.*(0.7 * L_Rossby + 88e3)
    else:
        raise NotImplementedError("Only implemented for 'global' and 'Atlantic'")


def munk_scale(lat0,A=4e7):
    """Compute boundary current length scale after Munk 1950 (eq. 24)
    from reference latitude and viscosity

    Parameters
    ----------
    lat0 : float
        Reference latitude to compute Rossby radius [deg]
    A : float
        Viscosity [m2 s-1]

    Returns
    -------
    Munk length scale in [m]
    """
    return 4*np.pi / np.sqrt(3) * (A / beta_coefficient(lat0))**(1./3)


def munk_scale_pop(dx):
    """Compute boundary current length scale after Munk 1950
    from model zonal grid spacing `dx` (m)
    after Jochum et al. 2008, JGR 113, equation A16
    """
    return 4*np.pi / np.sqrt(3) * 0.16**(1/3.) * dx


def munk_scale_pop_highres(lat,dx,AB0=2.7e10):
    """Compute boundary current length scale after Munk 1950
    from latitude `lat` (deg), model zonal grid spacing `dx` (m) and
    biharmonic viscosity at the equator `AB0` (m4 s-1), following
    Maltrud et al. 1998 (doi:10.1029/1998JC900013) and
    Smith et al. 2000 (doi:10.1175/1520-0485(2000)030<1532:NSOTNA>2.0.CO;2)
    """
    return munk_scale(lat,A=(AB0 * np.cos(np.radians(lat))**3 / dx**2))


def density_layer_mask(rho,rho0):
    """Return the interpolation matrix for a density layer at `rho0`

    Parameters
    ----------
    rho : array-like
        3D density field (z,y,x)
    rho0 : float
        Reference density
    """
    nz,ny,nx = rho.shape
    i,j = np.meshgrid(np.arange(nx),np.arange(ny))

    k1 = np.argmin(np.abs(rho-rho0),axis=0)
    
    diff = rho[k1,j,i]-rho0
    k2 = k1 - np.ma.array(np.sign(diff),int)
    k2[k2<0] = np.ma.masked

    c1 = np.abs(diff / (rho[k2,j,i] - rho[k1,j,i]))

    mask = np.ma.zeros((nz,ny,nx))
    mask[k1,j,i] = (1.-c1)
    mask[k2,j,i] = c1

    return mask


def wind_stress_curl(taux,tauy,dx,dy):
    """Compute wind stress curl from wind stress vector components taux,tauy and grid spacing dx,dy
    
    Parameters
    ----------
    taux,tauy : ndarrays [y,x]
        components of wind stress
    dx,dy : ndarrays [y,x]
        grid spacing
    """
    curl =  stats.central_differences(tauy,dx,axis=1)
    curl -= stats.central_differences(taux,dy,axis=0)
    curl = np.ma.masked_where(np.abs(curl)>1e20,curl)
    return curl


def sverdrup_transport(taux,tauy,dx,dy,lat0=30,landmask=None,rho0=1e3,return2D=True):
    """Compute the Sverdrup transport stream function from the wind stress
    
    Parameters
    ----------
    taux, tauy : 2D arrays (y,x)
        surface wind stress
    dx,dy : 2D arrays (y,x)
        grid spacing
    lat0 : float
        reference latitude for computing the beta coefficient
    landmask : 2D array, optional
        the field will be multipiled by this mask before integrating
    rho0 : float
        reference density
    return2D : bool
        whether to return the 2D zonal cumulative integral
        if false, the total zonal integral (1D) is returned

    Returns
    -------
    Sverdrup stream function in Sv (1e6 m3 s-1)
    """
    beta = beta_coefficient(lat0)
    
    if landmask is None:
        landmask = 1
    else:
        landmask = np.asarray(landmask,int)
    
    curl = wind_stress_curl(taux,tauy,dx,dy)
    curl *= landmask
    curl *= dx

    if return2D:
        integral = np.cumsum(curl[:,::-1],axis=1)[:,::-1]
        psi = -1./(beta*rho0) * integral * 1e-6
        psi = np.ma.masked_where(landmask==0,psi)
    else:
        integral = np.sum(curl,axis=1)
        psi = -1./(beta*rho0) * integral * 1e-6
    return psi

