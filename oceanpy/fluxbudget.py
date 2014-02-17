"""Computation of a flux budget for a masked region in a scalar/vector field"""
import numpy as np

# TODO: The _budget_2d_* functions should be made more memory efficient, avoiding to process the entire field and indexing with borderi,borderj from the beginning instead.

_facj,_faci = (np.array([-1,1,0,0]),np.array([0,0,-1,1]))
_faces_desc = ['bottom','top','left','right']


def _faces_idx(idxj,idxi):
    """Get indices for the faces"""
    return (np.tile(idxj,(len(_facj),1)) + _facj[:,np.newaxis],
            np.tile(idxi,(len(_faci),1)) + _faci[:,np.newaxis])


def _pad_2d(a):
    """Returns a copy of `a` with the outer values duplicated"""
    if isinstance(a,np.ma.masked_array):
        apad = np.ma.zeros(np.array(a.shape)+2,a.dtype)
    else:
        apad = np.zeros(np.array(a.shape)+2,a.dtype)
    apad[1:-1,1:-1] = a
    apad[1:-1,[0,-1]] = a[:,[-1,0]]
    return apad


def _outline_idx(mask):
    """Find the grid cells in mask that form the outline of the masked area
    and return their indices
    """
    maskedj,maskedi = np.where(mask)
    maskedjf,maskedif = _faces_idx(maskedj,maskedi)
    maskpad = _pad_2d(mask)
    atborder = ~np.all(maskpad[maskedjf+1,maskedif+1],axis=0)
    return (maskedj[atborder],maskedi[atborder])


def _find_open_faces(borderj,borderi,mask):
    """Find the faces between outline-cells and surrounding
    
    Note
    ----
    The vectors `_facj`,`_faci` define which of the 4 dimensions
    of the returned array corresponds to which face
    """
    borderjf,borderif = _faces_idx(borderj,borderi)
    maskpad = _pad_2d(mask)
    return ~maskpad[borderjf+1,borderif+1]


def _budget_2d_arakawa_b(borderj,borderi,faces,U,V,scalar):
    """Compute flux budget over a given region on an Arakawa-B grid"""
    # interpolate u,v to T grid INTERFACES (effectively an Arakawa C-grid)
    Uint = np.empty_like(U)
    Vint = np.empty_like(V)
    Uint[:-1,:] = (U[1:,:]+U[:-1,:])/2. # interpolate U in y-direction
    Vint[:,:-1] = (V[:,1:]+V[:,:-1])/2. # interpolate V in x-direction
    Uint[-1,:] = Uint[-2,:]
    Vint[:,-1] = Vint[:,-2]

    # pad all arrays by one grid point
    Uint = _pad_2d(Uint)
    Vint = _pad_2d(Vint)
    scalar = _pad_2d(scalar)

    # from here on, all indices to these arrays have to be +1
    borderjp1 = borderj+1
    borderip1 = borderi+1

    # generate an array defining whether the flow is _into_ the grid box
    into_box = np.vstack([
        (Vint[borderjp1-1,borderip1] > 0), # through bottom
        (Vint[borderjp1,borderip1] < 0), # through top
        (Uint[borderjp1,borderip1-1] > 0), # through left
        (Uint[borderjp1,borderip1] < 0), # through right
        ])
    out_of_box = ~into_box

    flux = np.vstack([
        Vint[borderjp1-1,borderip1], # through bottom
        Vint[borderjp1,borderip1], # through top
        Uint[borderjp1,borderip1-1], # through left
        Uint[borderjp1,borderip1], # through right
        ])
    absflux = np.abs(flux)

    borderjf,borderif = _faces_idx(borderj,borderi)
    borderjfp1,borderifp1 = borderjf+1,borderif+1

    sflux_into = scalar[borderjfp1,borderifp1] * absflux * faces * into_box
    sflux_outof = scalar[borderjp1,borderip1] * absflux * faces * out_of_box

    return np.sum(sflux_into - sflux_outof,axis=0)


def _budget_2d_arakawa_a(borderj,borderi,faces,U,V,scalar):
    """Compute flux budget over a given region on an Arakawa-A grid"""
    # interpolate u,v to T grid INTERFACES (effectively an Arakawa C-grid)
    Uint = np.empty_like(U)
    Vint = np.empty_like(V)
    Uint[:,:-1] = (U[:,1:]+U[:,:-1])/2.
    Vint[:-1,:] = (V[1:,:]+V[:-1,:])/2.
    Uint[:,-1] = Uint[:,-2]
    Vint[-1,:] = Vint[-2,:]

    # pad all arrays by one grid point
    Uint = _pad_2d(Uint)
    Vint = _pad_2d(Vint)
    scalar = _pad_2d(scalar)

    # from here on, all indices to these arrays have to be +1
    borderjp1 = borderj+1
    borderip1 = borderi+1

    # generate an array defining whether the flow is _into_ the grid box
    into_box = np.vstack([
        (Vint[borderjp1-1,borderip1] > 0), # through bottom
        (Vint[borderjp1,borderip1] < 0), # through top
        (Uint[borderjp1,borderip1-1] > 0), # through left
        (Uint[borderjp1,borderip1] < 0), # through right
        ])
    out_of_box = ~into_box

    flux = np.vstack([
        Vint[borderjp1-1,borderip1], # through bottom
        Vint[borderjp1,borderip1], # through top
        Uint[borderjp1,borderip1-1], # through left
        Uint[borderjp1,borderip1], # through right
        ])
    absflux = np.abs(flux)

    borderjf,borderif = _faces_idx(borderj,borderi)
    borderjfp1,borderifp1 = borderjf+1,borderif+1

    sflux_into = scalar[borderjfp1,borderifp1] * absflux * faces * into_box
    sflux_outof = scalar[borderjp1,borderip1] * absflux * faces * out_of_box

    return np.sum(sflux_into - sflux_outof,axis=0)


def _budget_2d_single_quantity(borderj,borderi,faces,U,V):
    """Compute flux budget over a given region where all quantities are on the same grid."""
    flux = np.vstack([
        V[borderj,borderi], # through bottom
        -V[borderj,borderi], # through top
        U[borderj,borderi], # through left
        -U[borderj,borderi], # through right
        ])
    flux_net = flux * faces
    return np.sum(flux_net,axis=0)


def budget_over_region_2D(U,V,scalar,mask,grid='ArakawaB'):
    """Compute the flux budget of `scalar` in the vector field `U`,`V` over the region `mask`
    
    Parameters
    ----------
    U,V : 2D arrays
        vector field of volume or mass flux (e.g. in m3/s or kg/s)
    scalar : 2D array
        scalar to advect
    mask : 2D boolean array
        region over which to compute the budget (True marks boxes inside the region)
    grid : str, e.g. 'ArakawaA'
        type of grid
    
    """
    if scalar is None:
        grid = 'single'
    else:
        scalar = np.ma.filled(scalar,0.)

    borderj,borderi = _outline_idx(mask)
    faces = _find_open_faces(borderj,borderi,mask)
 
    if grid == 'single':
        return np.sum(_budget_2d_single_quantity(borderj,borderi,faces,U,V)) or 0.
    elif grid == 'ArakawaA':
        return np.sum(_budget_2d_arakawa_a(borderj,borderi,faces,U,V,scalar)) or 0.
    elif grid == 'ArakawaB':
        return np.sum(_budget_2d_arakawa_b(borderj,borderi,faces,U,V,scalar)) or 0.
    else:
        raise NotImplementedError('No function (yet) for grid \'{}\'.'.format(grid))



def test_2d(U=None,V=None,scalar=None,mask=None,divergent=False,
        plotmask=False,plotborder=True,plotfaces=False,plotbudget=True,
        grid='ArakawaB'):
    """Test function for budget function with simplified flow field
    
    Options
    -------
    divergent : whether the field should be divergent
    plotmask : plot the mask
    plotborder : plot the border
    plotfaces : plot the four faces
    plotbudget : plot a 2D image of the flux budget
    """
    import matplotlib.pyplot as plt ; plt.close('all')

    # mask of region
    if mask is None:
        mask[5:10,10:15] = True
        mask[6:15,15:18] = True
    if plotmask:
        plt.figure()
        plt.pcolormesh(mask)
        plt.xlabel('i')
        plt.ylabel('j')
        plt.show()

    # flux through grid faces
    if U is None and V is None:
        u = np.ones(mask.shape)
        v = np.ones(mask.shape)*2.
        if divergent: u[:,:14] *= -1.
        dx = np.ones(mask.shape)
        dy = np.ones(mask.shape)
        dz = 1.
        U = u*dx*dz
        V = v*dy*dz

    # scalar that is to be transported
    if not scalar:
        scalar = np.ones(mask.shape)

    # get border outline indices
    borderj,borderi = _outline_idx(mask)

    if plotborder:
        border = np.zeros(mask.shape)
        border[borderj,borderi] = True
        plt.figure()
        plt.title('Region border and vector field')
        plt.pcolormesh(border)
        plt.quiver(U,V,color='w')
        plt.xlabel('i')
        plt.ylabel('j')
        plt.show()

    # find open faces
    faces = _find_open_faces(borderj,borderi,mask)
    
    if plotfaces:
        faces2d = np.zeros([4]+list(mask.shape))
        faces2d[:,borderj,borderi] = faces
        for k in xrange(4):
            plt.figure()
            plt.title(_faces_desc[k])
            plt.pcolormesh(faces2d[k])
            plt.xlabel('i')
            plt.ylabel('j')
            plt.show()

    # get indices of faces
    borderjf,borderif = _faces_idx(borderj,borderi)
    
    if grid == 'single':
        budget = _budget_2d_single_quantity(borderj,borderi,faces,U,V)
    elif grid == 'ArakawaA':
        budget = _budget_2d_arakawa_a(borderj,borderi,faces,U,V,scalar)
    elif grid == 'ArakawaB':
        budget = _budget_2d_arakawa_b(borderj,borderi,faces,U,V,scalar)
    else:
        raise NotImplementedError('Grid {0} not implemented.'.format(grid))
    print('Grid is \'{}\''.format(grid))

    if plotbudget:
        budget2d = np.zeros(mask.shape)
        budget2d[borderj,borderi] = budget
        budget = np.ma.masked_where(budget==0,budget)
        plt.figure()
        plt.title('Budget with losses and gains')
        img = plt.pcolormesh(budget2d)
        plt.colorbar(img)
        plt.xlabel('i')
        plt.ylabel('j')
        plt.show()

    print('Total budget: {}'.format(np.sum(budget)))


if __name__ is '__main__':
    
    test_2d(divergent=True,grid='ArakawaA')

