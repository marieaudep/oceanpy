"""Meta information about oceanographic variables"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

import prettycpt ; prettycpt.register_package_cmaps()


def _symmetric(a):
    return np.unique(np.concatenate((-np.asarray(a),a)))
def _symmetric0(a):
    return np.unique(np.concatenate((-np.asarray(a),[0],a)))


def num2str(x,mindecs=0):
    """Convert number `x` to string with same number of internal decimal places or `mindecs` trailing 0s

    Example
    -------

        >>> num2str([1, 1., 1.0, 1.2, 1.234])
        ['1', '1', '1', '1.2', '1.234']

    """
    def _stringify(x,mindecs):
        ssplt = str(x).split('.')
        if len(ssplt) == 1 or ssplt[-1] == '0':
            ndigits = mindecs
        else:
            ndigits = np.max((mindecs,len(ssplt[-1])))
        fmtstr = '{:.' + str(ndigits) + 'f}'        
        return fmtstr.format(x)
    return np.vectorize(_stringify)(x,mindecs)[()]


def num2fmtdict(nlst,mindecs=1):
    """Make format-dictionary from a list of numbers
    
    Paramters
    ---------
    nlst : lst
        list of numbers with different number of internal decimal places
    mindecs : int, optional
        can be used to set a minimum number of decimal places
    
    Returns
    -------
    dict
        dictionary with python format-strings
        to be used with pyplot.clabel(fmt=...)
        
    Example
    -------

        >>> num2fmtdict([1.,1.2,1.234])
        {1.0: '1.0', 1.2: '1.2', 1.234: '1.234'}
        
    """
    fmtdct = {}
    for x in nlst:
        ssplt = str(x).split('.')
        if len(ssplt) == 1 or ssplt[-1] == '0':
            ndigits = mindecs
        else:
            ndigits = np.max((mindecs,len(ssplt[-1])))
        fmtstr = '{:.' + str(ndigits) + 'f}'
        fmtdct[x] = fmtstr.format(x)
    return fmtdct



class VarMeta:
    """Meta information about oceanographic variables"""
    varmeta =  {
        'global': {
            'profile': {
                'vel' : {
                    'long': 'velocity',
                    'short': 'vel',
                    'units': 'm s-1',
                    'levels' : np.array([-16,-14,-12,-10,-8,-6,-4,-2,-1,1,2,4,6,8,10,12,14,16])/100.,
                    #'levels' : np.arange(-0.17,0.18,0.02),
                    'ticks' : np.array([-16,-12,-8,-4,-1,1,4,6,8,10,12,14,16])/100.,
                    'cmstr' : 'cpt-city/ncl/amwg/blueyellowred',},
                'velbc' : {
                    'long': 'velocity',
                    'short': 'vel',
                    'units': 'm s-1',
                    'levels' : _symmetric([2,4,6,8,10,15,20,25,30,35,40])/100.,
                    'ticks' : _symmetric([5,10,20,30,40])/100.,
                    'cmstr' : 'RdYlBu_r',},
                'spd' : {
                    'long': 'speed',
                    'short' : 'spd',
                    'units': 'm s-1',
                    'levels' : np.array([0,.5,1,1.5,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10,12,14,16,18,20,30])/100.,
                    'cmstr' : 'Spectral_r',},
                'dvdx' : {
                    'long' : 'velocity gradient',
                    'short' : 'dv/dx',
                    'units' : 'm s-1 m-1',
                    'levels' : _symmetric([5e-8,1e-7,5e-7,1e-6,5e-6,1e-5]),
                    'clabfmt' : '%.1e',
                    'cmstr' : 'RdBu_r',},
                'sal': {
                    'long' : 'salinity',
                    'short' : 'S',
                    'units': '',
                    'levels': np.arange(34.5,37,0.25),
                    'cmstr': 'jet',},
                'temp': {
                    'long': 'temperature',
                    'short' : 'T',
                    'units' : 'degC',
                    'levels' : [-1,-0.5,0,0.5,1.5,2.5,5,10,15,20,25],
                    'cmstr' : 'RdYlBu_r',},
                'depth': {
                    'long': 'depth',
                    'short': 'z',
                    'units': 'm',
                    'levels': [0,50,100,200,400,600,800,1000,
                               1200,1400,1600,1800,2000,
                               2200,2400,2600,2800,3000,
                               3200,3400,3600,3800,4000,
                               4200,4400,4600,4800,5000],
                    'cmstr': 'gist_ncar_r',},
                'psi': {
                    'long': 'stream function',
                    'short': 'psi',
                    'units': 'Sv',
                    'levels': np.arange(-52,52+1,4),
                    'cmstr': 'spectral'},
                'heatflux' : {
                    'long': 'heat flux',
                    'short': 'Qh',
                    'units': '10^6 W m-2',
                    'levels': np.arange(-3,3.1,0.2),
                    'cmstr': 'jet'},
                'saltflux' : {
                    'long': 'salt flux',
                    'short': 'SFL',
                    'units': 'm3 -s PPT',
                    'levels': np.arange(-6.5,6.5001,1.),
                    'cmstr': 'cpt-city/gmt/tn/GMT_panoply'}
                },
            'surface' : {
                'vel' : {
                    'long' : 'velocity',
                    'short' : 'vel',
                    'units': 'm s-1',
                    'levels' : np.arange(-0.20,0.20001,0.04),
                    'cmstr' : 'spectral'},
                'sal': {
                    'long': 'salinity',
                    'short': 'SSS',
                    'units': '',
                    'levels': np.arange(22,38.0001,0.5),
                    'cmstr': 'jet'}
                },
            'barotropic' : {
                'psi' : {
                    'long' : 'stream function',
                    'short' : 'psi',
                    'units': 'Sv',
                    'levels': np.arange(-54,54+1,4),
                    'cmstr' : 'spectral'},
                },
            'deep' : {
                'vel' : {
                    'long' : 'velocity',
                    'short' : 'vel',
                    'units': 'm s-1',
                    'levels' : _symmetric0([0.25,0.5,0.75,1,1.5,2,3,4,5,6])/100.,
                    'cmstr' : 'RdYlBu_r'},
                }},
        'northatlantic': {
            'surface': {
                'sal': {
                    'long': 'salinity',
                    'short': 'SSS',
                    'units': '',
                    'levels': np.arange(29,36+1,0.5),
                    'cmstr': 'jet'},
                'temp': {
                    'long': 'temperature',
                    'short': 'SST',
                    'units': 'degC',
                    'levels': [-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,3,4,5,6,7,8,10,15,20],
                    'clabfmt': '%1.1f',
                    'cmstr': 'RdYlBu_r'},
                'dens': {
                    'long': 'density anomaly',
                    'short': 'sigma',
                    'units': 'kg m-3',
                    'levels': [24,24.5,25,25.5,26,26.5,27,27.5,28],
                    'clabfmt': '%1.1f',
                    'cmstr': 'Spectral_r'},
                'psi': {
                    'long': 'stream function',
                    'short': 'psi',
                    'units': 'Sv',
                    'levels': np.arange(-0.50,0.5001,0.04),
                    'cmstr': 'spectral'},
                }},
        'arctic': {
            'surface': {
                'sal': {
                    'long': 'salinity',
                    'short': 'SSS',
                    'levels': [29,30,31,32,32.5,33,33.5,33.75,34,34.25,34.5,35],
                    'clabfmt': {},
                    'cmstr': 'jet'},
                'temp': {
                    'long': 'temperature',
                    'short': 'SST',
                    'units': 'degC',
                    'levels': np.linspace(-2,3,21),
                    'clabfmt': '%1.1f',
                    'cmstr': 'RdYlBu_r'},
                'dens': {
                    'long': 'density anomaly',
                    'short': 'sigma',
                    'units': 'kg m-3',
                    'levels': [24,24.5,25,25.5,26,26.5,27,27.5,28],
                    'clabfmt': '%1.1f',
                    'cmstr': 'Spectral_r'},
                'sic': {
                    'long': 'sea ice concentration',
                    'short': 'SIC',
                    'levels': np.linspace(0,1,21),
                    'levels_overlay': [0.1,0.2,0.4,0.6,0.8,1.0],
                    'clabfmt': '%1.1f',
                    'cmstr': None,
                    'cmstr_overlay': 'Greys'},
                'hice': {
                    'long': 'sea ice thickness',
                    'short': 'Hice',
                    'units': 'm',
                    'levels': np.arange(0,4+.1,0.25),
                    'levels_overlay': [0,1,2,3],
                    'clamfmt': '%1.1f',
                    'cmstr': None},
                }},
        'fjord': {
            'surface': {
                'sal': {
                    'long': 'salinity',
                    'short': 'SSS',
                    'units': '',
                    'levels': [29,30,31,32,32.5,33,33.5,33.75,34,34.25,34.5,35],
                    'clabfmt': {},
                    'cmstr': 'jet'},
                'temp': {
                    'long': 'temperature',
                    'short': 'SST',
                    'units': 'degC',
                    'levels': np.arange(0,6.001,1),
                    'clabfmt': '%1.1f',
                    'cmstr': 'RdYlBu_r'},
                'dens': {
                    'long': 'density anomaly',
                    'short': 'sigma',
                    'units': 'kg m-3',
                    'levels':[24,24.5,25,25.5,26,26.5,27,27.5,28],
                    'clabfmt': '%1.1f',
                    'cmstr': 'Spectral_r'}},
            'profile': {
                'sal': {
                    'long': 'salinity',
                    'short': 'S',
                    'units': '',
                    'levels': [30,31,32,32.5,33,33.5,33.75,34,34.25,34.5,34.75],
                    'cmstr': 'jet'},
                'temp': {
                    'long': 'temperature',
                    'short': 'T',
                    'units': 'degC',
                    'levels': [0, 1, 2, 3, 4, 5, 6],
                    'clabfmt': '%1.1f',
                    'cmstr': 'RdYlBu_r'},
                'dens': {
                    'long': 'density anomaly',
                    'short': 'sigma',
                    'units': 'kg m-3',
                    'levels': [24,25,26,27,28,29,29.5,30],
                    'clabfmt': '%1.0f',
                    'cmstr': 'Spectral_r'},
                }},
        }

    varnsynonyms = {'salt' : 'sal'}


    def __init__(self,varn=None,region='global',depth='profile',**kwargs):
        """Get meta information about the variable specified in the keyword arguments
        
        Parameters
        ----------
        varn : str
            variable name
        region : str
            geographical region
        depth : str
            vertical extension
        kwargs : dict
            additional keyword arguments are added to the VarMeta instance.

        see VarMeta.varmeta for all available options
        """
        if varn is not None:
            self.query(varn,region,depth)
        
        self.__dict__.update(kwargs)

        try:
            self.update_cmap()
        except:
            pass


    def query(self,varn,region='global',depth='profile'):
        """Query the variable meta information
        
        Parameters
        ----------
        varn : str
            variable name
        region : str
            geographical region
        depth : str
            vertical extension

        see VarMeta.varmeta for all available options
        """
        # get entry of self.varmeta
        try:
            varn,region,depth = map(str.lower,[varn,region,depth])
            try:
                varn = self.varnsynonyms[varn]
            except KeyError:
                pass
            vmdict = self.varmeta[region][depth][varn]
        except:
            raise

        # get all attributes from dict
        self.__dict__.update(vmdict)

        # update ticks and string representations
        self.update_ticks()

        # generate cmap, norm, etc.
        self.update_cmap()


    def update_ticks(self):
        # make ticks field equal to levels if not defined
        if not hasattr(self,'ticks'):
            self.ticks = self.levels

        # convert some fields to arrays and make string lists
        for key in ['levels','ticks','levels_overlay']:
            try:
                setattr(self,key,np.asarray(getattr(self,key)))
                setattr(self,(key+'_str'),num2str(getattr(self,key),mindecs=0))
            except AttributeError:
                pass


    def update_cmap(self,levels=None,basecmap=None,cmstr=None,update_ticks=True):
        """Generate cmap and norm
        
        Parameters
        ----------
        levels : array-like, optional
            updates self.levels before generation
        basecmap : color map to sub-sample, optional
            used instead of the color map defined in the vmdict
        cmapstr : str, optional
            color map name to retrieve using plt.cm.get_cmap()
        update_ticks : bool
            whether to update the ticks afterwards
        """
        if basecmap is None:
            if cmstr is not None:
                basecmap = plt.cm.get_cmap(cmstr)
                self.cmstr = cmstr
            if hasattr(self,'cmap'):
                basecmap = self.cmap
            elif hasattr(self,'cmstr'):
                basecmap = plt.cm.get_cmap(self.cmstr)
            else:
                raise ValueError('missing anything to be used as basecmap')
        if levels is not None:
            self.levels = levels
        colors = basecmap(np.linspace(0, 1, len(self.levels)))
        self.cmap = mcolors.ListedColormap(colors)
        self.norm = mcolors.BoundaryNorm(self.levels,self.cmap.N)
        # store vlim
        self.vmin = np.min(self.levels)
        self.vmax = np.max(self.levels)
        
        if update_ticks:
            self.update_ticks()

    def generate_clabfmt(self):
        self.clabfmt = num2fmtdict(self.levels)

