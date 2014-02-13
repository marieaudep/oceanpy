import numpy as np

def central_differences(y,dx,axis=0):
    """Compute gradient of an array `y` with spacing `dx` 
    using central differences along first dimension"""
    y,dx = np.asarray(y,dtype=np.float),np.asarray(dx)
    yT = np.swapaxes(y,0,axis)
    dxT = np.swapaxes(dx,0,axis)

    # make sure that dx is one column shorter than y
    if yT.shape[0] == dxT.shape[0]:
        dxT = dxT[:-1]

    if isinstance(yT,np.ma.masked_array):
        dydx = np.ma.zeros_like(yT)
    else:
        dydx = np.zeros_like(yT)

    # Second order forward difference for first element
    dydx[0] = -(3*yT[0] - 4*yT[1] + yT[2]) / (2*dxT[0])

    # Central difference interior elements
    dydx[1:-1] = (yT[2:] - yT[0:-2]) / (dxT[1:]+dxT[:-1])

    # Backwards difference final element
    dydx[-1] = (3*yT[-1] - 4*yT[-2] + yT[-3]) / (2*dxT[-1])
    
    # forward / backward differences at edges
    #dydx[0] = (yT[1]-yT[0])/dxT[0]
    #dydx[-1] = (yT[-1]-yT[-2])/dxT[-1]

    # central differences in the interior
    #dydx[1:-1] = (yT[2:] - yT[:-2])/(dxT[1:]+dxT[:-1])
    
    return np.swapaxes(dydx,axis,0)


def smooth(data,width=None,window=None,npasses=2,ckwargs=dict(mode='same')):
    """Smooth a data series by convolution with a windowfunction

    Parameters
    ----------
    data : ndarray
        data vector
    width : int
        width of rectangular window
    window : ndarray, optional
        window to be used instead of rectangular window
    npasses: int
        number of passes
    ckwargs : dict, optional
        keyword arguments to np.convolve
    """
    
    data = np.asanyarray(data)

    if width is None or width <= 1:
        # return unsmoothed signal
        print('Warning: Signal was not smoothed.')
        return data

    if window is None:
        window = np.ones(width)

    if npasses == 1:
        return np.convolve(window/window.sum(),data,**ckwargs)
    else:
        outdata = data.copy()
        for i in xrange(npasses):
            outdata = np.convolve(window/window.sum(),outdata,**ckwargs)
        return outdata



