""" 
Signal analysis tools.
NOTE: naming this module "signal.py" can sometimes create conflict with numpy

"""
import numpy as np
from numpy.random import rand
import pandas as pd


# --- List of available filters
FILTERS=[
    {'name':'Moving average'      , 'param':100 , 'paramName':'Window Size'  , 'paramRange':[1      , 100000] , 'increment':1  , 'digits':0} , 
    {'name':'Low pass 1st order'  , 'param':1.0, 'paramName':'Cutoff Freq.' , 'paramRange':[0.0001 , 100000] , 'increment':0.1, 'digits':4} , 
    {'name':'High pass 1st order' , 'param':0.1 , 'paramName':'Cutoff Freq.' , 'paramRange':[0.0001 , 100000] , 'increment':0.1, 'digits':4} , 
]

SAMPLERS=[
    {'name':'Replace', 'param':[], 'paramName':'New x'},
    {'name':'Insert',  'param':[], 'paramName':'Insert list'},
    {'name':'Remove',  'param':[], 'paramName':'Remove list'},
    {'name':'Every n', 'param':2  , 'paramName':'n'},
    {'name':'Linspace', 'param':[0,1,100]  , 'paramName':'xmin, xmax, n'},
    {'name':'Time-based', 'param':0.01  , 'paramName':'Sample time (s)'},
    {'name':'Delta x', 'param':[0.1,np.nan,np.nan], 'paramName':'dx, xmin, xmax'},
]



def reject_outliers(y, x=None, m = 2., replaceNaN=True):
    """ Reject outliers:
        If replaceNaN is true: they are replaced by NaN 
        Otherwise they are removed
    """
    if m==0: 
        # No rejection...
        pass
    else:
        dd = np.abs(y - np.nanmedian(y))
        mdev = np.nanmedian(dd)
        if mdev:
            ss = dd/mdev 
            b=ss<m
            if replaceNaN:
                y=y.copy()
                y[~b]=np.nan
            else:
                y=y[b]
                if x is not None:
                    x= x[b]
    if x is None:
        return y
    else:
        return x, y


# --------------------------------------------------------------------------------}
# --- Resampling 
# --------------------------------------------------------------------------------{
def multiInterp(x, xp, fp, extrap='bounded'):
    """ 
    Interpolate all the columns of a matrix `fp` based on new values `x`
    INPUTS:
      - x  : array ( n ), new values
      - xp : array ( np ), old values
      - fp : array ( nval, np), matrix values to be interpolated
    """
    # Sanity
    x   = np.asarray(x)
    xp  = np.asarray(xp)
    assert fp.shape[1]==len(xp), 'Second dimension of fp should have the same length as xp'

    j = np.searchsorted(xp, x, 'left') - 1
    dd  = np.zeros(len(x)) #*np.nan
    bOK = np.logical_and(j>=0, j< len(xp)-1)
    jOK = j[bOK]
    dd[bOK] = (x[bOK] - xp[jOK]) / (xp[jOK + 1] - xp[jOK])
    jBef=j 
    jAft=j+1
    # 
    bLower =j<0
    bUpper =j>=len(xp)-1
    # Use first and last values for anything beyond xp
    jAft[bUpper] = len(xp)-1
    jBef[bUpper] = len(xp)-1
    jAft[bLower] = 0
    jBef[bLower] = 0
    if extrap=='bounded':
        #OK
        pass
    elif extrap=='nan':
        # Set values to nan if out of bounds
        bBeyond= np.logical_or(x<np.min(xp), x>np.max(xp))
        dd[bBeyond] = np.nan
    else:
        raise NotImplementedError()

    return (1 - dd) * fp[:,jBef] + fp[:,jAft] * dd

def interpArray(x, xp, fp, extrap='bounded'):
    """ 
    Interpolate all the columns of a matrix `fp` based on one new value `x`
    INPUTS:
      - x  : scalar  new values
      - xp : array ( np ), old values
      - fp : array ( nval, np), matrix values to be interpolated
    """
    # Sanity
    xp  = np.asarray(xp)
    assert fp.shape[1]==len(xp), 'Second dimension of fp should have the same length as xp'

    if fp.shape[1]==0:
        raise Exception('Second dimension of fp should be >0')

    j   = np.searchsorted(xp, x) - 1
    if j<0:
        # Before bounds
        if extrap=='bounded':
            return fp[:,0]
        elif extrap=='nan':
            return fp[:,0]*np.nan
        else:
            raise NotImplementedError()

    elif j>=len(xp)-1:
        # After bounds
        if extrap=='bounded':
            return fp[:,-1]
        elif extrap=='nan':
            return fp[:,-1]*np.nan
        else:
            raise NotImplementedError()
    else:
        # Normal case, within bounds
        dd = (x- xp[j]) / (xp[j+1] - xp[j])
        return (1 - dd) * fp[:,j] + fp[:,j+1] * dd


def interpDF(x_new, xLabel, df, extrap='bounded'):
    """ Resample a dataframe using linear interpolation"""
    x_old = df[xLabel].values
    #x_new=np.sort(x_new)
    # --- Method 1 (pandas)
    #df_new = df_old.copy()
    #df_new = df_new.set_index(x_old)
    #df_new = df_new.reindex(df_new.index | x_new)
    #df_new = df_new.interpolate().loc[x_new]
    #df_new = df_new.reset_index()
    # --- Method 2 interp storing dx
    data_new=multiInterp(x_new, x_old, df.values.T, extrap=extrap)
    df_new = pd.DataFrame(data=data_new.T, columns=df.columns.values)
    df_new[xLabel] = x_new # Just in case this value was replaced by nan..
    return df_new


def resample_interp(x_old, x_new, y_old=None, df_old=None):
    #x_new=np.sort(x_new)
    if df_old is not None:
        # --- Method 1 (pandas)
        #df_new = df_old.copy()
        #df_new = df_new.set_index(x_old)
        #df_new = df_new.reindex(df_new.index | x_new)
        #df_new = df_new.interpolate().loc[x_new]
        #df_new = df_new.reset_index()
        # --- Method 2 interp storing dx
        data_new=multiInterp(x_new, x_old, df_old.values.T)
        df_new = pd.DataFrame(data=data_new.T, columns=df_old.columns.values)
        return df_new

    if y_old is not None:
        return np.interp(x_new, x_old, y_old)


def applySamplerDF(df_old, x_col, sampDict):
    x_old=df_old[x_col].values
    x_new, df_new =applySampler(x_old, y_old=None, sampDict=sampDict, df_old=df_old)
    df_new[x_col]=x_new
    return df_new


def applySampler(x_old, y_old, sampDict, df_old=None):

    param = np.asarray(sampDict['param']).ravel()

    if sampDict['name']=='Replace':
        if len(param)==0:
            raise Exception('Error: At least one value is required to resample the x values with')
        x_new = param
        return x_new, resample_interp(x_old, x_new, y_old, df_old)

    elif sampDict['name']=='Insert':
        if len(param)==0:
            raise Exception('Error: provide a list of values to insert')
        x_new = np.sort(np.concatenate((x_old.ravel(),param)))
        return x_new, resample_interp(x_old, x_new, y_old, df_old)

    elif sampDict['name']=='Remove':
        I=[]
        if len(param)==0:
            raise Exception('Error: provide a list of values to remove')
        for d in param:
            Ifound= np.where(np.abs(x_old-d)<1e-3)[0]
            if len(Ifound)>0:
                I+=list(Ifound.ravel())
        x_new=np.delete(x_old,I)
        return x_new, resample_interp(x_old, x_new, y_old, df_old)

    elif sampDict['name']=='Delta x':
        if len(param)==0:
            raise Exception('Error: provide value for dx')
        dx    = param[0]
        if dx==0:
            raise Exception('Error: `dx` cannot be 0')
        if len(param)==1:
            # NOTE: not using min/max if data loops (like airfoil)
            xmin =  np.nanmin(x_old) 
            xmax =  np.nanmax(x_old) + dx/2
        elif len(param)==3:
            xmin  = param[1]
            xmax  = param[2]
            if np.isnan(xmin):
                xmin =  np.nanmin(x_old) 
            if np.isnan(xmax):
                xmax =  np.nanmax(x_old) + dx/2
        else:
            raise Exception('Error: the sampling parameters should be a list of three values `dx, xmin, xmax`')
        x_new = np.arange(xmin, xmax, dx)
        if len(x_new)==0:
            xmax = xmin+dx*1.1 # NOTE: we do it like this to account for negative dx
            x_new = np.arange(xmin, xmax, dx)
        param = [dx, xmin, xmax]
        return x_new, resample_interp(x_old, x_new, y_old, df_old)

    elif sampDict['name']=='Linspace':
        if len(param)!=3:
            raise Exception('Error: Provide three parameters for linspace: xmin, xmax, n')
        xmin  = float(param[0])
        xmax  = float(param[1])
        n     = int(param[2])
        x_new = np.linspace(xmin, xmax, n)
        return x_new, resample_interp(x_old, x_new, y_old, df_old)

    elif sampDict['name']=='Every n':
        if len(param)==0:
            raise Exception('Error: provide value for n')
        n = int(param[0])
        if n==0:
            raise Exception('Error: |n| should be at least 1')

        x_new=x_old[::n]
        if df_old is not None:
            return x_new, (df_old.copy()).iloc[::n,:]
        if y_old is not None:
            return x_new, y_old[::n]

    elif sampDict['name'] == 'Time-based':
        if len(param) == 0:
            raise Exception('Error: provide value for new sampling time')
        sample_time = float(param[0])
        if sample_time <= 0:
            raise Exception('Error: sample time must be positive')

        time_index = pd.TimedeltaIndex(x_old, unit="S")
        x_new = pd.Series(x_old, index=time_index).resample("{:f}S".format(sample_time)).mean().interpolate().values

        if df_old is not None:
            df_new = df_old.set_index(time_index, inplace=False).resample("{:f}S".format(sample_time)).mean()
            df_new = df_new.interpolate().reset_index(drop=True)
            return x_new, df_new
        if y_old is not None:
            y_new = pd.Series(y_old, index=time_index).resample("{:f}S".format(sample_time)).mean()
            y_new = y_new.interpolate().values
            return x_new, y_new

    else:
        raise NotImplementedError('{}'.format(sampDict))
    pass

# --------------------------------------------------------------------------------}
# --- Filters
# --------------------------------------------------------------------------------{
#     def moving_average(x, w):
#         #t_new    = np.arange(0,Tmax,dt)
#         #nt      = len(t_new)
#         #nw=400
#         #u_new = moving_average(np.floor(np.linspace(0,3,nt+nw-1))*3+3.5, nw)
#         return np.convolve(x, np.ones(w), 'valid') / w
#     def moving_average(x,N,mode='same'):
#        y=np.convolve(x, np.ones((N,))/N, mode=mode)
#        return y
def moving_average(a, n=3) :
    """ 
    perform moving average, return a vector of same length as input

    NOTE: also in kalman.filters
    """
    a   = a.ravel()
    a   = np.concatenate(([a[0]]*(n-1),a)) # repeating first values
    ret = np.cumsum(a, dtype = float)
    ret[n:] = ret[n:] - ret[:-n]
    ret=ret[n - 1:] / n
    return ret

def lowpass1(y, dt, fc=3) :
    """ 
    1st order low pass filter
    """
    tau=1/(2*np.pi*fc)
    alpha=dt/(tau+dt)
    y_filt=np.zeros(y.shape)
    y_filt[0]=y[0]
    for i in np.arange(1,len(y)):
        y_filt[i]=alpha*y[i] + (1-alpha)*y_filt[i-1]
    return y_filt

def highpass1(y, dt, fc=3) :
    """ 
    1st order high pass filter
    """
    tau=1/(2*np.pi*fc)
    alpha=tau/(tau+dt)
    y_filt=np.zeros(y.shape)
    y_filt[0]=0
    for i in np.arange(1,len(y)):
        y_filt[i]=alpha*y_filt[i-1] + alpha*(y[i]-y[i-1])
    m0=np.mean(y)
    m1=np.mean(y_filt)
    y_filt+=m0-m1
    return y_filt


def applyFilter(x, y, filtDict):
    if filtDict['name']=='Moving average':
        return moving_average(y, n=np.round(filtDict['param']).astype(int))
    elif filtDict['name']=='Low pass 1st order':
        dt = x[1]-x[0]
        return lowpass1(y, dt=dt, fc=filtDict['param'])
    elif filtDict['name']=='High pass 1st order':
        dt = x[1]-x[0]
        return highpass1(y, dt=dt, fc=filtDict['param'])
    else:
        raise NotImplementedError('{}'.format(filtDict))

def applyFilterDF(df_old, x_col, options):
    """ apply filter on a dataframe """
    # Brute force loop
    df_new = df_old.copy()
    x = df_new[x_col]
    for (colName, colData) in df_new.iteritems():
        if colName != x_col:
            df_new[colName] = applyFilter(x, colData, options)
    return df_new


# --------------------------------------------------------------------------------}
# ---  
# --------------------------------------------------------------------------------{
def zero_crossings(y, x=None, direction=None, bouncingZero=False):
    """
      Find zero-crossing points in a discrete vector, using linear interpolation.

      direction: 'up' or 'down', to select only up-crossings or down-crossings
      bouncingZero: also returns zeros that are exactly zero and do not change sign

      returns: 
          x values xzc such that y(yzc)==0
          indexes izc, such that the zero is between y[izc] (excluded) and y[izc+1] (included)

      if direction is not provided, also returns:
              sign, equal to 1 for up crossing
    """
    if x is None:
        x=np.arange(len(y))
    else:
        x = np.asarray(x)
    y = np.asarray(y)

    if np.any((x[1:] - x[0:-1]) <= 0.0):
        raise Exception('x values need to be in ascending order')

    # Indices before zero-crossing
    iBef = np.where(y[1:]*y[0:-1] < 0.0)[0]
    
    # Find the zero crossing by linear interpolation
    xzc = x[iBef] - y[iBef] * (x[iBef+1] - x[iBef]) / (y[iBef+1] - y[iBef])
    
    # Selecting points that are exactly 0 and where neighbor change sign
    iZero = np.where(y == 0.0)[0]
    iZero = iZero[np.where((iZero > 0) & (iZero < x.size-1))]
    if not bouncingZero:
        iZero = iZero[np.where(y[iZero-1]*y[iZero+1] < 0.0)] # we only accept zeros that change signs

    # Concatenate 
    xzc  = np.concatenate((xzc, x[iZero]))
    iBef = np.concatenate((iBef, iZero))

    # Sort
    iSort = np.argsort(xzc)
    xzc, iBef = xzc[iSort], iBef[iSort]

    # Return up-crossing, down crossing or both
    sign = np.sign(y[iBef+1]-y[iBef])
    if direction == 'up':
        I= np.where(sign==1)[0]
        return xzc[I],iBef[I]
    elif direction == 'down':
        I= np.where(sign==-1)[0]
        return xzc[I],iBef[I]
    elif direction is not None:
        raise Exception('Direction should be either `up` or `down` or `None`')
    return xzc, iBef, sign


# --------------------------------------------------------------------------------}
# --- Correlation  
# --------------------------------------------------------------------------------{
def correlation(x, nMax=80, dt=1, method='numpy'):
    """ 
    Compute auto correlation of a signal
    """

    def acf(x, nMax=20):
        return np.array([1]+[np.corrcoef(x[:-i], x[i:])[0,1]  for i in range(1, nMax)])


    nvec   = np.arange(0,nMax)
    if method=='manual':
        sigma2 = np.var(x)
        R    = np.zeros(nMax)
        R[0] =1
        for i,nDelay in enumerate(nvec[1:]):
            R[i+1] = np.mean(  x[0:-nDelay] * x[nDelay:]  ) / sigma2
            #R[i+1] = np.corrcoef(x[:-nDelay], x[nDelay:])[0,1] 

    elif method=='numpy':
        R= acf(x, nMax=nMax)
    else:
        raise NotImplementedError()

    tau = nvec*dt
    return R, tau
# Auto-correlation comes in two versions: statistical and convolution. They both do the same, except for a little detail: The statistical version is normalized to be on the interval [-1,1]. Here is an example of how you do the statistical one:
# 
# 
# def autocorr(x):
#     result = numpy.correlate(x, x, mode='full')
#     return result[result.size/2:]





def correlated_signal(coeff, n=1000, seed=None):
    """
    Create a correlated random signal of length `n` based on the correlation coefficient `coeff`
          value[t] = coeff * value[t-1]  + (1-coeff) * random
    """
    if coeff<0 or coeff>1: 
        raise Exception('Correlation coefficient should be between 0 and 1')
    if seed is not None:
        np.random.seed(seed)

    x    = np.zeros(n)
    rvec = rand(n)
    x[0] = rvec[0]
    for m in np.arange(1,n):
        x[m] = coeff*x[m-1] + (1-coeff)*rvec[m] 
    x-=np.mean(x)
    return x


def find_time_offset(t, f, g, outputAll=False):
    """ 
    Find time offset between two signals (may be negative)

    t_offset = find_time_offset(t, f, g)
    f(t+t_offset) ~= g(t)

    """
    import scipy
    from scipy.signal import correlate
    # Remove mean and normalize by std
    f  = f.copy()
    g  = g.copy()
    f -= f.mean() 
    g -= g.mean()
    f /= f.std()
    g /= g.std()

    # Find cross-correlation
    xcorr = correlate(f, g)

    # Lags
    n   = len(f)
    dt  = t[1]-t[0]
    lag = np.arange(1-n, n)*dt

    # Time offset is located at maximum correlation
    t_offset = lag[xcorr.argmax()]

    if outputAll:
        return t_offset, lag, xcorr
    else:
        return t_offset

def amplitude(x, t=None, T = None, mask=None, debug=False):
    """
    Compute signal amplitude (max-min)/2.
    If a frequency is provided, the calculation is the average on each period

    x: signal time series
    mask - time at which transient starts 
    """
    if mask is not None:
        x = x[mask]
        if t is not None:
            t = t[mask]
    #
    if T is not None and t is not None:
        t -= t[0]
        if t[-1]<=T:
            return (np.max(x)-np.min(x))/2
        n = int(t[-1]/T)
        A = 0
        for i in range(n):
            b = np.logical_and(t<=(i+1)*T ,  t>=i*T)
            A+=(np.max(x[b])-np.min(x[b]))/2
        A/=n

        if debug:
            import matplotlib.pyplot as plt
            from welib.tools.colors import python_colors
            fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
            fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
            ax.plot(t,  x-np.mean(x) ,'k-', lw=3, label='Original')
            for i in range(n):
                b = np.logical_and(t<=(i+1)*T ,  t>=i*T)
                A=(np.max(x[b])-np.min(x[b]))/2
                ax.plot(t[b]- i*T,  x[b]-np.mean(x[b]), c=python_colors(i),  label='A={}'.format(A))
                ax.plot([0,0,T,T,0],[-A,A,A,-A,-A] ,'--', c=python_colors(i)  )
            ax.set_xlabel('time')
            ax.set_ylabel('x')
            ax.legend()
        return A

        # split signals into  subsets
    else:
        return (np.max(x)-np.min(x))/2

def phase_shift(A, B, t, omega, tStart=0, deg=True, debug=False):
    """
    A: reference signal
    omega: expected frequency of reference signal
    """
    b =t>=tStart
    t = t[b]
    A = A[b]
    B = B[b]
    t_offset, lag, xcorr = find_time_offset(t, A, B, outputAll=True)
    phi0 = t_offset*omega # Phase offset in radians
    phi  = np.mod(phi0, 2*np.pi)  
    if phi > 0.8 * np.pi:
        phi = phi-2*np.pi
    if deg: 
        phi *=180/np.pi
        phi0*=180/np.pi
        if phi<-190:
            phi+=360
#     if debug:
#         raise NotImplementedError()
    return phi 

def input_output_amplitude_phase(t, u, y, omega_u=None, A_u=None, deg=True, mask=None, debug=False):
    """ 
    Return amplitude ratio and phase shift between a reference input `u` and output `y`
       Typically used when the input is a sinusoidal signal and the output is "similar".

    INPUTS:
     - t: time vector, length nt
     - u: input time series, length nt
     - y: output time series, length nt
     - omega_u: cyclic frequency, required to convert time offset to phase
                when provided, amplitude ratio is computed for each possible periods, and averaged
     - A_u : amplitude of input signal (typically known if y is a sinusoid)
     - deg: phase is returned in degrees
     - mask: mask to be applied to t, u, y
    """
    if mask is not None:
        t = t[mask]
        u = u[mask]
        y = y[mask]
    if omega_u is None:
        raise NotImplementedError()
        T=None
    else:
        T = 2*np.pi/omega_u

    # --- Amplitude ratio
    A_y = amplitude(y, t, T=T, debug=debug)
    if A_u is None:
        A_u = amplitude(u, t, T=T)
    G = A_y/A_u

    # --- Phase shift
    phi = phase_shift(u, y, t, omega_u, deg=deg, debug=debug)

    return G, phi
    

def sine_approx(t, x, method='least_square'):
    """ 
    Sinusoidal approximation of input signal x
    """
    if method=='least_square':
        from welib.tools.curve_fitting import fit_sinusoid
        y_fit, pfit, fitter = fit_sinusoid(t, x)
        omega = fitter.model['coeffs']['omega']
        A     = fitter.model['coeffs']['A']
        phi   = fitter.model['coeffs']['phi']
        x2    = y_fit
    else:
        raise NotImplementedError()


    return x2, omega, A, phi


# --------------------------------------------------------------------------------}
# --- Convolution 
# --------------------------------------------------------------------------------{
def convolution_integral(time, f, g, method='auto'):
    r"""
    Compute convolution integral:
       f * g = \int 0^t f(tau) g(t-tau) dtau  = g * f
    For now, only works for uniform time vector, an exception is raised otherwise

    method=['auto','direct','fft'], 
        see scipy.signal.convolve
        see scipy.signal.fftconvolve
    """
    from scipy.signal import convolve
    dt = time[1]-time[0] 
    if len(np.unique(np.around(np.diff(time)/dt,3)))>1:
        raise Exception('Convolution integral implemented for uniform time vector')

    #np.convolve(f.ravel(), g.ravel() )[:len(time)]*dt
    return convolve(f.ravel(), g.ravel() )[:len(time)]*dt


# --------------------------------------------------------------------------------}
# --- Intervals/peaks 
# --------------------------------------------------------------------------------{
def intervals(b, min_length=1, forgivingJump=True, removeSmallRel=True, removeSmallFact=0.1, mergeCloseRel=False, mergeCloseFact=0.2):
    """
    Describe intervals from a boolean vector where intervals are indicated by True

    INPUT:
      - b         : a logical vector, where 1 means, I'm in an interval.
      - min_length: if provided, do not return intervals of length < min_length
      - forgivingJump: if true, merge intervals that are separated by a distance < min_length
      - removeSmallRel: remove intervals that have a small length compared to the max length of intervals
      - removeSmallFact: factor used for removeSmallRel
      - mergeCloseRel: merge intervals that are closer than a fraction of the typical distance between intervals

    OUTPUTS:
        - IStart : ending  indices
        - IEnd  :  ending  indices
        - Length:  interval lenghts (IEnd-IStart+1)

    IStart, IEnd, Lengths = intervals([False, True, True, False, True, True, True, False])
        np.testing.assert_equal(IStart , np.array([1,4]))
        np.testing.assert_equal(IEnd   , np.array([2,6]))
        np.testing.assert_equal(Lengths, np.array([2,3]))
    """
    b = np.asarray(b)
    total = np.sum(b)

    min_length=max(min_length,1)
    if forgivingJump:
        min_jump=min_length
    else:
        min_jump=1

    if total==0:
        IStart = np.array([])
        IEnd   = np.array([])
        Lengths= np.array([])
        return IStart, IEnd, Lengths
    elif total==1:
        i = np.where(b)[0][0]
        IStart = np.array([i])
        IEnd   = np.array([i])
        Lengths= np.array([1])
    else:
        n = len(b)
        Idx = np.arange(n)[b]
        delta_Idx=np.diff(Idx)
        jumps =np.where(delta_Idx>min_jump)[0]
        if len(jumps)==0:
            IStart = np.array([Idx[0]])
            IEnd   = np.array([Idx[-1]])
        else:
            istart=Idx[0]
            jumps=np.concatenate(([-1],jumps,[len(Idx)-1]))
            IStart = Idx[jumps[:-1]+1] # intervals start right after a jump
            IEnd   = Idx[jumps[1:]]    # intervals stops at jump
        Lengths = IEnd-IStart+1

    # Removing intervals smaller than min_length
    bKeep   = Lengths>=min_length
    IStart  = IStart[bKeep]
    IEnd    = IEnd[bKeep]
    Lengths = Lengths[bKeep]
    # Removing intervals smaller than less than a fraction of the max interval
    if removeSmallRel:
        bKeep   = Lengths>=removeSmallFact*np.max(Lengths)
    IStart  = IStart[bKeep]
    IEnd    = IEnd[bKeep]
    Lengths = Lengths[bKeep]

    # Distances between intervals
    if mergeCloseRel:
        if len(IStart)<=2:
            pass
        else:
            D = IStart[1:]-IEnd[0:-1]
            #print('D',D,np.max(D),int(np.max(D) * mergeCloseFact))
            min_length = max(int(np.max(D) * mergeCloseFact), min_length)
            if min_length<=1:
                pass 
            else:
                #print('Readjusting min_length to {} to accomodate for max interval spacing of {:.0f}'.format(min_length, np.mean(D)))
                return intervals(b, min_length=min_length, forgivingJump=True, removeSmallRel=removeSmallRel, removeSmallFact=removeSmallFact, mergeCloseRel=False)
    return IStart, IEnd, Lengths

def peaks(x, threshold=0.3, threshold_abs=True, method='intervals', min_length=3,
         mergeCloseRel=True, returnIntervals=False):
    """
    Find peaks in a signal, above a given threshold
    INPUTS:
     - x         : 1d-array, signal
     - threshold : scalar, absolute or relative threshold beyond which peaks are looked for
                   relative threshold are proportion of the max-min of the signal (between 0-1)
     - threshold_abs : boolean, specify whether the threshold is absolute or relative
     - method    : string, selects which method is used to find the peaks, between:
                   - 'interval'  : one peak per interval above the threshold
                   - 'derivative': uses derivative to find maxima, may return more than one per interval
     - min_length: 
          - if 'interval'   method is used: minimum interval
          - if 'derivative' method is used: minimum distance between two peaks 

    OPTIONS for interval method:
      - mergeCloseRel: logical, if True, attempts to merge intervals that are close to each other compare to the typical interval spacing
                       set to False if all peaks are wanted
      - returnIntervals: logical, if true, return intervals used for interval method
    OUTPUTS:
      - I : index of the peaks 
      -[IStart, IEnd] if return intervals is true, see function `intervals`


    see also:
       scipy.signal.find_peaks

    """
    if not threshold_abs:
        threshold = threshold * (np.max(y) - np.min(y)) + np.min(y)

    if method =='intervals':
        IStart, IEnd, Lengths = intervals(x>threshold, min_length=min_length, mergeCloseRel=mergeCloseRel)
        I = np.array([iS if L==1 else np.argmax(x[iS:iE+1])+iS for iS,iE,L in zip(IStart,IEnd,Lengths)])
        if returnIntervals:
            return I, IStart, IEnd
        else:
            return I
            
    elif method =='derivative':
        I = indexes(x, thres=threshold, thres_abs=True, min_dist=min_length)
        return I
    else:
        raise NotImplementedError('Method {}'.format(method))



# --------------------------------------------------------------------------------}
# --- Simple signals 
# --------------------------------------------------------------------------------{
def impulse(time, t0=0, A=1, epsilon=None, **kwargs):
    """ 
    returns a dirac  function:
      A/dt    if t==t0
      0    otherwise

    Since the impulse response is poorly defined in discrete space, it's recommended 
    to use a smooth_delta. See the welib.tools.functions.delta
    """
    from .functions import delta
    t=np.asarray(time)-t0
    y= delta(t, epsilon=epsilon, **kwargs)*A
    return y

def step(time, t0=0, valueAtStep=0, A=1):
    """ 
    returns a step function:
      0     if t<t0
      A     if t>t0
      valueAtStep if t==t0

    NOTE: see also welib.tools.functions.Pi
    """
    return np.heaviside(time-t0, valueAtStep)*A

def ramp(time, t0=0, valueAtStep=0, A=1):
    """ 
    returns a ramp function:
      0          if t<t0
      A*(t-t0)   if t>=t0

    NOTE: see also welib.tools.functions.Pi
    """
    t=np.asarray(time)-t0
    y=np.zeros(t.shape)
    y[t>=0]=A*t[t>=0]
    return y


def hat(time, T=1, t0=0, A=1, method='abs'):
    """ 
    returns a hat function:
      A*hat   if |t-t0|<T/2
      0       otherwise 
    T : full time length of hat 
    A : Amplitude of hat

    NOTE: for an integral of 1, one needs: T*A=2
    see also: scipy.signal.windows.triang
    """
    t=np.asarray(time)-t0

    if method == 'abs':
        # use definition in terms of absolute value (recommended)
        y=np.zeros(t.shape)
        b= np.abs(t)<=T/2
        y[b]=A*(1-np.abs(2*t[b]/T))
    elif method == 'sum':
        # For fun, use summation of ramps
        y1= ramp(time, t0=t0-T/2, A=A/T*2)
        y2= ramp(time, t0=t0    , A=-A/T*4)
        y3= ramp(time, t0=t0+T/2, A=A/T*2)
        y= y1+y2+y3
    else:
        raise NotImplementedError()

    return y




if __name__=='__main__':
    import numpy as np
    import matplotlib.pyplot as plt

    # Input
    dt    = 1
    n     = 10000
    coeff = 0.95 # 1:full corr, 00-corr
    nMax  = 180
    # Create a correlated time series
    tvec  = np.arange(0,n)*dt
    ts = correlated_signal(coeff, n)
    # --- Compute correlation coefficient
    R, tau = correlation(x, nMax=nMax)
    fig,axes = plt.subplots(2, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    ax=axes[0]
    # Plot time series
    ax.plot(tvec,ts)
    ax.set_xlabel('t [s]')
    ax.set_ylabel('u [m/s]')
    ax.tick_params(direction='in')
    # Plot correlation
    ax=axes[1]
    ax.plot(tau,  R              ,'b-o', label='computed')
    ax.plot(tau, coeff**(tau/dt) , 'r--' ,label='coeff^{tau/dt}') # analytical coeff^n trend
    ax.set_xlabel(r'$\tau$ [s]')
    ax.set_ylabel(r'$R(\tau)$ [-]')
    ax.legend()
    plt.show()






