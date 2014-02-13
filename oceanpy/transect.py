"""Tools for defining and working with transects on a model grid"""
import numpy as np
import matplotlib.pyplot as plt
import os.path

from pygeo import geo

from . import varmeta

class ZonalTransect:
    """Zonal transect in a near-regular lon/lat grid

    Parameters
    ----------
    longrid,latgrid : 2D arrays
        lon, lat grid vertices on which the transect lies
    lat0 : float
        reference latitude at or close to which to draw the transect
    lonw, lone : float, optional
        western and eastern boundaries for zonal transect
        if not provided, the mask of latgrid (if existent)
        limits the lateral extent of the transect

    Attributes added
    ----------------
    lon,lat : lon,lat axes of u vertices
    ii,j : indices of transect on lon,lat grid

    Example
    -------
    
        tr = ZonalTransect(longrid,latgrid,lat0,lonw,lone)
        print(tr.ii)
        print(tr.j)
        tr.plot_map(longrid,latgrid)

    """

    def __init__(self,longrid,latgrid,lat0,
                 lonw=None,lone=None):
        self.lat0 = lat0
        
        longrid = np.mod(longrid[:],360)

        # mask the latgrid for mean taking
        if lonw is None or lone is None:
            latgrid_masked = latgrid
        else:
            self.lone = np.mod(lone,360) or 360
            self.lonw = np.mod(lonw,360) or 360

            # limit the selectable region zonally
            if self.lonw <= self.lone:
                latgrid_masked = np.ma.masked_where(
                    ((longrid >= self.lonw) & (longrid <= self.lone)),latgrid[:])
            else:
                latgrid_masked = np.ma.masked_where(
                    ((longrid >= self.lonw) | (longrid <= self.lone)),latgrid[:])

        # find position of best fitting zonal section
        latmean = np.mean(latgrid_masked,axis=1)
        self.j = j = np.argmin(np.abs(latmean-lat0))

        # find the zonal boundaries of the section
        lonband = longrid[j,:]
        if lonw is None or lone is None:
            regmask = (~np.ma.getmaskarray(latgrid_masked))[j,:]
            if regmask.all(): # nothing masked
                iw,ie = 0,len(lonband)-1
            elif regmask[[0,-1]].all(): # region wraps around
                # find largest connected region of False
                cr = np.cumsum(regmask)
                regmask_largest = cr == np.argmax(np.bincount(cr))
                ie,iw = np.where(regmask_largest)[0][[0,-1]] ; ie+=1
            else:
                iw,ie = np.where(regmask)[0][[0,-1]]
        else:
            iw = np.argmin(np.abs(lonband-self.lonw))
            ie = np.argmin(np.abs(lonband-self.lone))

        if iw > ie:
            self.ii = np.concatenate([np.arange(iw,len(lonband)),np.arange(0,ie+1)])
        else:
            self.ii = np.arange(iw,ie+1)

        # make sure longitude axis is monotonically increasing
        lontmp = lonband[self.ii]
        lonbreaks = np.where(np.diff(lontmp)<-180)[0]
        for br in lonbreaks:
            lontmp[:(br+1)] -= 360
        if np.all(lontmp > 180):
            lontmp -= 360
        self.lon = lontmp
        self.lat = np.ma.getdata(latgrid[j,self.ii])

        self.meanlat = np.mean(self.lat)

        # add descriptions
        self.desc = 'zonal transect around {:.2f} N'.format(self.meanlat)
        self.lat0str = '{:.0f}S'.format(-lat0) if lat0 < 0 else '{:.0f}N'.format(lat0)
        self.meanlatstr = '{:.2f}S'.format(-self.meanlat) if lat0 < 0 else '{:.2f}N'.format(self.meanlat)

    def add_grid_bounds(self,lon_p=None,lat_p=None):
        """Generate grid cell boundaries
        either by sub-sampling lon and lat or from input arguments

        Parameters
        ----------
        lon_p,lat_p : 2D arrays, optional
            lon,lat axes of grid cell boundaries

        Attributes added
        ----------------
        lon_p,lat_p : lon,lat axes of p vertices (u boundaries)
        """
        if lon_p is None or lat_p is None:
            lo,la = geo.waypoints_segments(self.lon,self.lat,n=2)
            self.lon_p = lo[1:len(lo):2]
            self.lat_p = la[1:len(la):2]
        else:
            self.lon_p = lon_p[self.j,self.ii]
            self.lat_p = lat_p[self.j,self.ii]

    def add_dx_dist(self,dxgrid=None,dxgrid_u=None):
        """Generate grid spacing and distance axis
        either by calculating distances between lon/lat points or from input arguments

        Parameters
        ----------
        dxgrid : 2D array, optional
            distance between grid vertices (centered at boundaries)
        dxgrid_u : 2D array, optional
            distance between cell boundaries (centered at vertices)

        Attributes added
        ----------------
        dist : cumulative distance along transect from west to east
        dist_p : cumulative distance along transect at u grid boundaries (first point <0)
        """
        # compute transect axis distance in km
        if dxgrid is None:
            self.dx = geo.haversine(self.lon[ :-1],self.lat[ :-1],
                                    self.lon[1:  ],self.lat[1:  ])
        else:
            self.dx = dxgrid[self.j,self.ii[:-1]]
        self.dist = np.concatenate([np.zeros(1),np.cumsum(self.dx)])*1e-3
        
        # generate grid cell boundary axis
        if dxgrid_u is not None:
            self.dx_u = dxgrid_u[j,self.ii]
            self.dist_p = (np.concatenate(
                    [np.zeros(1),np.cumsum(self.dx_u[:-1])])-0.5*self.dx_u[0])*1e-3
        

    def plot_map(self,**plotkwargs):
        """Plot a the transect on a map. 
        
        Parameters
        ----------
        **plotkwargs : keyword arguments, optional
            passed on to plot_zonal_transect_map
        """
        
        return plot_zonal_transect_map(trs=self,**plotkwargs)


def plot_zonal_transect_map(trs,m,meridians,parallels,
                            bathy=None,lon=None,lat=None,
                            figdir=None,prefix='',ext='.pdf'):
    """Plot a map of the transect. If bathy is provided, use as background.

    Parameters
    ----------
    trs : zonal transect object
       Zonal transect data
    m : mpl_tools.basemap instance
       Name of map projection defined in jbntools.maps
    bathy : ndarray, optional
       Bathymetry for background
    lon,lat : ndarrays, optional [required for bathy]
       Grid coordinates
    figdir : str, optional
       Directory where to save figures
    prefix : str, optional
       Figure name prefix
    """

    fig = plt.figure()
    ax = fig.gca()

    m.drawmeridians(meridians,labels=[0,0,0,1],zorder=9)
    m.drawparallels(parallels,labels=[1,0,0,0],zorder=9)

    if bathy is not None and lon is not None and lat is not None:        
        xx,yy = m(lon,lat)
        vm = varmeta.VarMeta('depth')
        cnf = m.contourf(xx,yy,bathy,
                         levels=vm.levels,cmap=vm.cmap,norm=vm.norm,
                         extend='max',zorder=3)
        cb = plt.colorbar(cnf,shrink=0.8)
        cb.ax.set_ylim(cb.ax.get_ylim()[::-1])
    else:
        m.drawcoastlines(zorder=3)

    def _draw_transect(tr,label=None):
        mx,my = m(tr.lon,tr.lat)
        m.plot(mx,my,'rx',zorder=8)

    if isinstance(trs,dict):
        for tr,key in trs.iteritems():
            _draw_transect(tr,label=key)
        fname = 'zonal_transects_{}_map'.format(prefix)
    elif isinstance(trs,list):
        for tr in trs:
            _draw_transect(tr)
        fname = 'zonal_transects_{}_map'.format(prefix)
    else:
        tr = trs
        ax.set_title(tr.desc)
        _draw_transect(tr)
        fname = 'zonal_transect_{}_{}_map'.format(tr.lat0str,prefix)

    if figdir:
        figname = os.path.join(figdir,fname+ext)
        if ext == '.pdf':
            fig.savefig(figname,bbox_inches='tight',pad_inches=0.5)
        elif ext in ['.png','.jpg']:
            fig.savefig(figname,dpi=600)
        else:
            fig.savefig(figname)
        plt.close(fig)
    else:
        fig.show()
    return fig

        
def plot_meridional_velocity_transect(vvel,z,tr,vlines=[],xmode='lon',
                                      figdir=None,prefix='',use_contourf=False):
    """Plot a meridional velocity transect using pcolor

    Parameters
    ----------
    vvel : ndarray or netcdf4 variable
        meridional velocity transect (z,x)
    z : ndarray
        Vertical axis (cell boundaries)
    tr : ZonalTransect instance
        Transect spline data
    vlines : list, optional
        vertical lines to draw at specified positions
    xmode : str, 'lon' or 'dist'
        Plot x axis in coordinates or distance from origin
    figdir : str, optional
        Path to figure directory. If provided, figures are saved
    prefix : str, optional
        Prefix for figure file name (e.g. dataset ID)
    use_contourf : bool
        Use contourf instead of pcolormesh to draw map background
    """
    fig = plt.figure(figsize=(10, 4))
    ax = fig.gca()
    ax.set_title(tr.desc)

    if xmode == 'lon':
        x = tr.lon
    elif xmode == 'dist':
        x = tr.dist
        ax.set_xlabel('Distance in km')
    else:
        raise NotImplementedError

    vm = varmeta.VarMeta('vel')
    if use_contourf:
        cnf = ax.contourf(x,z,vvel,
                          levels=vm.levels,cmap=vm.cmap,norm=vm.norm,
                          extend='both',zorder=3)
    else:
        cnf = ax.pcolormesh(x,z,vvel,
                            cmap=vm.cmap,norm=vm.norm,
                            vmin=vm.vmin,vmax=vm.vmax,
                            zorder=3)
    ax.set_ylabel('Depth in m')
    ax.set_ylim(np.max(z),0)
    ax.set_xlim(x[0],x[-1])
    ax.grid(True)
    
    cb = plt.colorbar(cnf)
    cb.set_label('{} in {}'.format(vm.long,vm.units))

    # add vertical lines
    for vlinekwargs in vlines:
        ax.axvline(zorder=9,**vlinekwargs)

    # save figure if figdir is provided
    if figdir:
        figname = os.path.join(figdir,'zonal_transect_{0}_{1}_vel.pdf'.format(tr.lat0str,prefix))
        fig.savefig(figname,bbox_inches='tight')
        plt.close(fig)
    else:
        fig.show()
    return fig

