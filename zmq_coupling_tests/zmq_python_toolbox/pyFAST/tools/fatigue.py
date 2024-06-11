"""
Tools for fatigue analysis


Main functions:
- equivalent_load: calculate damage equivalent load for a given signal
- find_range_count: returns range and number of cycles for a given signal

Subfunctions:
- eq_load: calculate equivalent loads using one of the two rain flow counting methods
- cycle_matrix: calculates a matrix of cycles (binned on amplitude and mean value)
- eq_load_and_cycles: calculate eq_loads of multiple time series (e.g. life time equivalent load)


Main aglorithms for rain flow counting:
- rainflow_windap: taken from [2], based on [3]
- rainflow_astm: taken from [2], based [4]
- fatpack: using [5]


References:
  [1] Hayman (2012) MLife theory manual for Version 1.00
  [2] Wind energy toolbox, wetb.fatigue_tools, DTU wind energy, Denmark
  [3] "Recommended Practices for Wind Turbine Testing - 3. Fatigue Loads", 2. edition 1990, Appendix A
  [4] Adam Nieslony - Rainflow Counting Algorithm, MATLAB Central File Exchange 
      http://www.mathworks.com/matlabcentral/fileexchange/3026)
  [5] Fatpack - Python package
      https://github.com/Gunnstein/fatpack


"""
import warnings
import numpy as np


__all__  = ['equivalent_load', 'find_range_count']
__all__  += ['rainflow_astm', 'rainflow_windap','eq_load','eq_load_and_cycles','cycle_matrix','cycle_matrix2']


class SignalConstantError(Exception):
    pass


def equivalent_load(time, signal, m=3, Teq=1, bins=100, method='rainflow_windap',
        meanBin=True, binStartAt0=False,
        outputMore=False, debug=False):
    """Equivalent load calculation

    Calculate the damage equivalent load for a given signal and a given Wohler exponent

    INPUTS
     - time : array-like, the time values corresponding to the signal (s)
     - signals : array-like, the load signal
     - m :    Wohler exponent (default is 3)
     - Teq : The equivalent period (Default 1, for 1Hz)
     - bins : Number of bins in rainflow count histogram
     - method: rain flow counting algorithm: 'rainflow_windap', 'rainflow_astm' or 'fatpack'
     - meanBin: if True, use the mean of the ranges within a bin (recommended)
              otherwise use the middle of the bin (not recommended).
     - binStartAt0: if True bins start at zero. Otherwise, start a lowest range
     - outputMore: if True, returns range, cycles and bins as well

    OUTPUTS
     - Leq : the equivalent load for given m and Teq

        or (if outputMore is True )

     - Leq, S, N, bins, DELi: 
        - S: ranges
        - N: cycles
        - bins: bin edges
        - DELi: component 'i' of the DEL (for cycle i)
    """
    time   = np.asarray(time)
    signal = np.asarray(signal)

    # Remove nan, might not be the cleanest
    b = ~np.isnan(signal)
    signal = signal[b]
    time   = time[b]

    try:
        if len(time)<=1:
            raise Exception()
        if type(time[0]) is np.datetime64:
            T = T/np.timedelta64(1,'s') # or T.item().total_seconds()
        else:
            T = time[-1]-time[0] # time length of signal (s). Will fail for signal of length 1
        if T==0:
            raise Exception()

        neq = T/Teq # number of equivalent periods, see Eq. (26) of [1]

        # --- Range (S) and counts (N)
        N, S, bins = find_range_count(signal, bins=bins, method=method, meanBin=meanBin, binStartAt0=binStartAt0)

        # --- get DEL 
        DELi = S**m * N / neq
        Leq = DELi.sum() ** (1/m)     # See e.g. eq. (30) of [1]

    except:
        if outputMore:
            return np.nan, np.nan, np.nan, np.nan, np.nan
        else:
            return np.nan

    if debug:
        for i,(b,n,s,DEL) in enumerate(zip(bins, N, S, DELi)):
            if n>0:
                print('Bin {:3d}: [{:6.1f}-{:6.1f}] Mid:{:6.1f} - Mean:{:6.1f} Counts:{:4.1f} DEL:{:8.1f} Fraction:{:3.0f}%'.format(i,b,bins[i+1],(b+bins[i+1])/2,s,n,DEL,DEL/Leq**m*100))
    if outputMore:
        return Leq, S, N, bins, DELi
    else:
        return Leq

 
def find_range_count(signal, bins, method='rainflow_windap', meanBin=True, binStartAt0=True):
    """
    Returns number of cycles `N` for each range range `S` 
    Equidistant bins are setup based on the min/max of the signal.
    INPUTS:
     - signal: array 
     - bins : 1d-array, int
         If bins is a sequence, left edges (and the rightmost edge) of the bins.
         If bins is an int, a sequence is created dividing the range `min`--`max` of signal into `bins` number of equally sized bins.
    OUTPUTS:
      - N: number of cycles for each bin
      - S: Ranges for each bin
           S is either the center of the bin (meanBin=False) 
               or 
           S is the mean of the ranges within this bin (meanBin=True)
      - S_bin_edges: edges of the bins
    """

    if method in rainflow_func_dict.keys():
        rainflow_func = rainflow_func_dict[method]
        try:
            N, S, S_bin_edges, _, _ = cycle_matrix(signal, ampl_bins=bins, mean_bins=1, rainflow_func=rainflow_func, binStartAt0=binStartAt0)
        except SignalConstantError:
            return np.nan, np.nan, np.nan

        S_bin_edges = S_bin_edges.flatten()
        N           = N.flatten()
        S           = S.flatten()
        S_mid = (S_bin_edges[:-1] + S_bin_edges[1:]) / 2
        if not meanBin:
            S=S_mid

    elif method=='fatpack':
        import fatpack
        # find rainflow ranges
        try:
            ranges = fatpack.find_rainflow_ranges(signal)
        except IndexError:
            # Currently fails for constant signal
            return np.nan, np.nan, np.nan
        # --- Legacy fatpack
        # if (not binStartAt0) and (not meanBin):
        #    N, S = fatpack.find_range_count(ranges, bins)
        # --- Setup bins
        # If binStartAt0 is True, the same bins as WINDAP are used
        S_bin_edges = create_bins(ranges, bins, binStartAt0=binStartAt0)
        # --- Using bin_count to get value at center of bins 
        N, S = bin_count(ranges, S_bin_edges, meanBin=meanBin)

    else:
        raise NotImplementedError('Rain flow algorithm {}'.format(method))

    # Remove NaN
    b = np.isnan(S)
    S[b] = 0
    N[b] = 0

    return N, S, S_bin_edges

def create_bins(x, bins, binStartAt0=False):
    """ 
    Equidistant bins are setup based on the min/max of the x, unless the user provided the bins as a sequence.
    INPUTS:
     - x: array 
     - bins : 1d-array, int
         If bins is a sequence, left edges (and the rightmost edge) of the bins.
         If bins is an int, a sequence is created dividing the range `min`--`max` of x into `bins` number of equally sized bins.
    OUTPUTS: 
      - bins:
    """
    if isinstance(bins, int):
        xmax = np.max(x)
        xmin, xmax = np.min(x), np.max(x)
        if binStartAt0:
            xmin = 0
        else:
            xmin = np.min(x)
            if xmin==xmax:
                # I belive that's what's done by histogram. double check
                xmin=xmin-0.5
                xmax=xmax+0.5
        bins = np.linspace(xmin, xmax, num=bins + 1) 
    return bins


def bin_count(x, bins, meanBin=True):
    """ 
    Return counts of x within bins
    """
    if not meanBin:
        # Use the middle of the bin
        N, bns = np.histogram(x, bins=bins)
        S = bns[:-1] + np.diff(bns) / 2.
    else:
        bins = create_bins(x, bins, binStartAt0=False)
        import pandas as pd
        df = pd.DataFrame(data=x, columns=['x'])
        xmid = (bins[:-1]+bins[1:])/2
        df['x_mid']= pd.cut(df['x'], bins= bins, labels = xmid ) # Adding a column that has bin attribute
        df2        = df.groupby('x_mid').mean()   # Average by bin
        df['N']  = 1
        dfCount       = df[['N','x_mid']].groupby('x_mid').sum()
        df2['N'] = dfCount['N']
        # Just in case some bins are missing (will be nan)
        df2       = df2.reindex(xmid)
        df2 = df2.fillna(0)
        S = df2['x'].values
        N = df2['N'].values
    return N, S

 


def check_signal(signal):
    # check input data validity
    if not type(signal).__name__ == 'ndarray':
        raise TypeError('signal must be ndarray, not: ' + type(signal).__name__)

    elif len(signal.shape) not in (1, 2):
        raise TypeError('signal must be 1D or 2D, not: ' + str(len(signal.shape)))

    if len(signal.shape) == 2:
        if signal.shape[1] > 1:
            raise TypeError('signal must have one column only, not: ' + str(signal.shape[1]))
    if np.min(signal) == np.max(signal):
        raise SignalConstantError("Signal is constant, cannot compute DLC and range")


def rainflow_windap(signal, levels=255., thresshold=(255 / 50)):
    """Windap equivalent rainflow counting


    Calculate the amplitude and mean values of half cycles in signal

    This algorithms used by this routine is implemented directly as described in
    "Recommended Practices for Wind Turbine Testing - 3. Fatigue Loads", 2. edition 1990, Appendix A

    Parameters
    ----------
    Signal : array-like
        The raw signal

    levels : int, optional
        The signal is discretize into this number of levels.
        255 is equivalent to the implementation in Windap

    thresshold : int, optional
        Cycles smaller than this thresshold are ignored
        255/50 is equivalent to the implementation in Windap

    Returns
    -------
    ampl : array-like
        Peak to peak amplitudes of the half cycles

    mean : array-like
        Mean values of the half cycles


    Examples
    --------
    >>> signal = np.array([-2.0, 0.0, 1.0, 0.0, -3.0, 0.0, 5.0, 0.0, -1.0, 0.0, 3.0, 0.0, -4.0, 0.0, 4.0, 0.0, -2.0])
    >>> ampl, mean = rainflow_windap(signal)
    """
    check_signal(signal)
    #type <double> is required by <find_extreme> and <rainflow>
    signal = signal.astype(np.double)
    if np.all(np.isnan(signal)):
        return None
    offset = np.nanmin(signal)
    signal -= offset
    if np.nanmax(signal) > 0:
        gain = np.nanmax(signal) / levels
        signal = signal / gain
        signal = np.round(signal).astype(int)


        # If possible the module is compiled using cython otherwise the python implementation is used


        #Convert to list of local minima/maxima where difference > thresshold
        sig_ext = peak_trough(signal, thresshold)


        #rainflow count
        ampl_mean = pair_range_amplitude_mean(sig_ext)

        ampl_mean = np.array(ampl_mean)
        ampl_mean = np.round(ampl_mean / thresshold) * gain * thresshold
        ampl_mean[:, 1] += offset
        return ampl_mean.T



def rainflow_astm(signal):
    """Matlab equivalent rainflow counting

    Calculate the amplitude and mean values of half cycles in signal

    This implemementation is based on the c-implementation by Adam Nieslony found at
    the MATLAB Central File Exchange http://www.mathworks.com/matlabcentral/fileexchange/3026

    Parameters
    ----------
    Signal : array-like
        The raw signal

    Returns
    -------
    ampl : array-like
        peak to peak amplitudes of the half cycles (note that the matlab implementation
        uses peak amplitude instead of peak to peak)

    mean : array-like
        Mean values of the half cycles


    Examples
    --------
    >>> signal = np.array([-2.0, 0.0, 1.0, 0.0, -3.0, 0.0, 5.0, 0.0, -1.0, 0.0, 3.0, 0.0, -4.0, 0.0, 4.0, 0.0, -2.0])
    >>> ampl, mean = rainflow_astm(signal)
    """
    check_signal(signal)

    # type <double> is reuqired by <find_extreme> and <rainflow>
    signal = signal.astype(np.double)

    # Import find extremes and rainflow.
    # If possible the module is compiled using cython otherwise the python implementation is used

    # Remove points which is not local minimum/maximum
    sig_ext = find_extremes(signal)

    # rainflow count
    ampl_mean = np.array(rainflowcount(sig_ext))

    return np.array(ampl_mean).T


def eq_load(signals, no_bins=46, m=[3, 4, 6, 8, 10, 12], neq=1, rainflow_func=rainflow_windap):
    """Equivalent load calculation

    Calculate the equivalent loads for a list of Wohler exponent and number of equivalent loads

    Parameters
    ----------
    signals : list of tuples or array_like
        - if list of tuples: list must have format [(sig1_weight, sig1),(sig2_weight, sig1),...] where\n
            - sigx_weight is the weight of signal x\n
            - sigx is signal x\n
        - if array_like: The signal
    no_bins : int, optional
        Number of bins in rainflow count histogram
    m : int, float or array-like, optional
        Wohler exponent (default is [3, 4, 6, 8, 10, 12])
    neq : int, float or array-like, optional
        The equivalent number of load cycles (default is 1, but normally the time duration in seconds is used)
    rainflow_func : {rainflow_windap, rainflow_astm}, optional
        The rainflow counting function to use (default is rainflow_windap)

    Returns
    -------
    eq_loads : array-like
        List of lists of equivalent loads for the corresponding equivalent number(s) and Wohler exponents

    Examples
    --------
    >>> signal = np.array([-2.0, 0.0, 1.0, 0.0, -3.0, 0.0, 5.0, 0.0, -1.0, 0.0, 3.0, 0.0, -4.0, 0.0, 4.0, 0.0, -2.0])
    >>> eq_load(signal, no_bins=50, neq=[1, 17], m=[3, 4, 6], rainflow_func=rainflow_windap)
    [[10.311095426959747, 9.5942535021382174, 9.0789213365013932], # neq = 1, m=[3,4,6]
    [4.010099657859783, 4.7249689509841746, 5.6618639965313005]], # neq = 17, m=[3,4,6]

    eq_load([(.4, signal), (.6, signal)], no_bins=50, neq=[1, 17], m=[3, 4, 6], rainflow_func=rainflow_windap)
    [[10.311095426959747, 9.5942535021382174, 9.0789213365013932], # neq = 1, m=[3,4,6]
    [4.010099657859783, 4.7249689509841746, 5.6618639965313005]], # neq = 17, m=[3,4,6]
    """
    try:
        return eq_load_and_cycles(signals, no_bins, m, neq, rainflow_func)[0]
    except TypeError:
        return [[np.nan] * len(np.atleast_1d(m))] * len(np.atleast_1d(neq))


def eq_load_and_cycles(signals, no_bins=46, m=[3, 4, 6, 8, 10, 12], neq=[10 ** 6, 10 ** 7, 10 ** 8], rainflow_func=rainflow_windap):
    """Calculate combined fatigue equivalent load

    Parameters
    ----------
    signals : list of tuples or array_like
        - if list of tuples: list must have format [(sig1_weight, sig1),(sig2_weight, sig1),...] where\n
            - sigx_weight is the weight of signal x\n
            - sigx is signal x\n
        - if array_like: The signal
    no_bins : int, optional
        Number of bins for rainflow counting
    m : int, float or array-like, optional
        Wohler exponent (default is [3, 4, 6, 8, 10, 12])
    neq : int or array-like, optional
        Equivalent number, default is [10^6, 10^7, 10^8]
    rainflow_func : {rainflow_windap, rainflow_astm}, optional
        The rainflow counting function to use (default is rainflow_windap)

    Returns
    -------
    eq_loads : array-like
        List of lists of equivalent loads for the corresponding equivalent number(s) and Wohler exponents
    cycles : array_like
        2d array with shape = (no_ampl_bins, 1)
    ampl_bin_mean : array_like
        mean amplitude of the bins
    ampl_bin_edges
        Edges of the amplitude bins
    """
    cycles, ampl_bin_mean, ampl_bin_edges, _, _ = cycle_matrix(signals, no_bins, 1, rainflow_func)
    if 0:  #to be similar to windap
        ampl_bin_mean = (ampl_bin_edges[:-1] + ampl_bin_edges[1:]) / 2
    cycles, ampl_bin_mean = cycles.flatten(), ampl_bin_mean.flatten()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        #DEL =      [[(          (cycles * ampl_bin_mean ** _m) / _neq)  for _m in np.atleast_1d(m)]  for _neq in np.atleast_1d(neq)]
        eq_loads = [[((np.nansum(cycles * ampl_bin_mean ** _m) / _neq) ** (1. / _m)) for _m in np.atleast_1d(m)]  for _neq in np.atleast_1d(neq)]
    return eq_loads, cycles, ampl_bin_mean, ampl_bin_edges


def cycle_matrix(signals, ampl_bins=10, mean_bins=10, rainflow_func=rainflow_windap, binStartAt0=True):
    """Markow load cycle matrix

    Calculate the Markow load cycle matrix

    Parameters
    ----------
    Signals : array-like or list of tuples
        - if array-like, the raw signal\n
        - if list of tuples, list of (weight, signal), e.g. [(0.1,sig1), (0.8,sig2), (.1,sig3)]\n
    ampl_bins : int or array-like, optional
        if int, Number of amplitude value bins (default is 10)
        if array-like, the bin edges for amplitude
    mean_bins : int or array-like, optional
        if int, Number of mean value bins (default is 10)
        if array-like, the bin edges for mea
    rainflow_func : {rainflow_windap, rainflow_astm}, optional
        The rainflow counting function to use (default is rainflow_windap)
    binStartAt0 : boolean 
        Start the bins at 0. Otherwise, start at the min of ranges

    Returns
    -------
    cycles : ndarray, shape(ampl_bins, mean_bins)
        A bi-dimensional histogram of load cycles(full cycles). Amplitudes are\
        histogrammed along the first dimension and mean values are histogrammed along the second dimension.
    ampl_bin_mean : ndarray, shape(ampl_bins,)
        The average cycle amplitude of the bins
    ampl_edges : ndarray, shape(ampl_bins+1,)
        The amplitude bin edges
    mean_bin_mean : ndarray, shape(ampl_bins,)
        The average cycle mean of the bins
    mean_edges : ndarray, shape(mean_bins+1,)
        The mean bin edges

    Examples
    --------
    >>> signal = np.array([-2.0, 0.0, 1.0, 0.0, -3.0, 0.0, 5.0, 0.0, -1.0, 0.0, 3.0, 0.0, -4.0, 0.0, 4.0, 0.0, -2.0])
    >>> cycles, ampl_bin_mean, ampl_edges, mean_bin_mean, mean_edges = cycle_matrix(signal)
    >>> cycles, ampl_bin_mean, ampl_edges, mean_bin_mean, mean_edges = cycle_matrix([(.4, signal), (.6,signal)])
    """

    if isinstance(signals[0], tuple):
        weights, ampls, means = np.array([(np.zeros_like(ampl)+weight,ampl,mean) for weight, signal in signals for ampl,mean in rainflow_func(signal[:]).T], dtype=np.float64).T
    else:
        ampls, means = rainflow_func(signals[:])
        weights = np.ones_like(ampls)
    if isinstance(ampl_bins, int):
        ampl_bins = create_bins(ampls[weights>0], ampl_bins, binStartAt0=binStartAt0)
    cycles, ampl_edges, mean_edges = np.histogram2d(ampls, means, [ampl_bins, mean_bins], weights=weights)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ampl_bin_sum = np.histogram2d(ampls, means, [ampl_bins, mean_bins], weights=weights * ampls)[0]
        ampl_bin_mean = np.nanmean(ampl_bin_sum / np.where(cycles,cycles,np.nan),1)
        mean_bin_sum = np.histogram2d(ampls, means, [ampl_bins, mean_bins], weights=weights * means)[0]
        mean_bin_mean = np.nanmean(mean_bin_sum / np.where(cycles, cycles, np.nan), 1)
    cycles = cycles / 2  # to get full cycles
    return cycles, ampl_bin_mean, ampl_edges, mean_bin_mean, mean_edges


def cycle_matrix2(signal, nrb_amp, nrb_mean, rainflow_func=rainflow_windap):
    """
    Same as wetb.fatigue.cycle_matrix but bin from min_amp to
    max_amp instead of 0 to max_amp.

    Parameters
    ----------

    Signal : ndarray(n)
        1D Raw signal array

    nrb_amp : int
        Number of bins for the amplitudes

    nrb_mean : int
        Number of bins for the means

    rainflow_func : {rainflow_windap, rainflow_astm}, optional
        The rainflow counting function to use (default is rainflow_windap)

    Returns
    -------

    cycles : ndarray, shape(ampl_bins, mean_bins)
        A bi-dimensional histogram of load cycles(full cycles). Amplitudes are\
        histogrammed along the first dimension and mean values are histogrammed
        along the second dimension.

    ampl_edges : ndarray, shape(no_bins+1,n)
        The amplitude bin edges

    mean_edges : ndarray, shape(no_bins+1,n)
        The mean bin edges

    """
    bins = [nrb_amp, nrb_mean]
    ampls, means = rainflow_func(signal)
    weights = np.ones_like(ampls)
    cycles, ampl_edges, mean_edges = np.histogram2d(ampls, means, bins,
                                                    weights=weights)
    cycles = cycles / 2  # to get full cycles

    return cycles, ampl_edges, mean_edges

# --------------------------------------------------------------------------------}
# --- Rainflowcount_astm.py
# --------------------------------------------------------------------------------{
'''
Created on 27/02/2013

@author: mmpe

How to use:

import_cython("cy_rainflowcount",'cy_rainflowcount.py','')
from cy_rainflowcount import find_extremes,rainflow

ext = find_extremes(np.array([-2,0,1,0,-3,0,5,0,-1,0,3,0,-4,0,4,0,-2]).astype(np.double))
print rainflow(ext)
'''
def find_extremes(signal):  #cpdef find_extremes(np.ndarray[double,ndim=1] signal):
    """return indexes of local minima and maxima plus first and last element of signal"""

    #cdef int pi, i
    # sign of gradient
    sign_grad = np.int8(np.sign(np.diff(signal)))

    # remove plateaus(sign_grad==0) by sign_grad[plateau_index]=sign_grad[plateau_index-1]
    plateau_indexes, = np.where(sign_grad == 0)
    if len(plateau_indexes) > 0 and plateau_indexes[0] == 0:
        # first element is a plateau
        if len(plateau_indexes) == len(sign_grad):
                # All values are equal to crossing level!
                return np.array([0])

        # set first element = first element which is not a plateau and delete plateau index
        i = 0
        while sign_grad[i] == 0:
            i += 1
        sign_grad[0] = sign_grad[i]

        plateau_indexes = np.delete(plateau_indexes, 0)

    for pi in plateau_indexes.tolist():
        sign_grad[pi] = sign_grad[pi - 1]

    extremes, = np.where(np.r_[1, (sign_grad[1:] * sign_grad[:-1] < 0), 1])

    return signal[extremes]


def rainflowcount(sig):  #cpdef rainflowcount(np.ndarray[double,ndim=1] sig):
    """Cython compilable rain ampl_mean count without time analysis


    This implemementation is based on the c-implementation by Adam Nieslony found at
    the MATLAB Central File Exchange http://www.mathworks.com/matlabcentral/fileexchange/3026

    References
    ----------
    Adam Nieslony, "Determination of fragments of multiaxial service loading
    strongly influencing the fatigue of machine components,"
    Mechanical Systems and Signal Processing 23, no. 8 (2009): 2712-2721.

    and is based on the following standard:
    ASTM E 1049-85 (Reapproved 1997), Standard practices for cycle counting in
    fatigue analysis, in: Annual Book of ASTM Standards, vol. 03.01, ASTM,
    Philadelphia, 1999, pp. 710-718.

    Copyright (c) 1999-2002 by Adam Nieslony

    Ported to Cython compilable Python by Mads M Pedersen
    In addition peak amplitude is changed to peak to peak amplitude


    """

    #cdef int sig_ptr, index
    #cdef double ampl
    a = []
    sig_ptr = 0
    ampl_mean = []
    for _ in range(len(sig)):
        a.append(sig[sig_ptr])
        sig_ptr += 1
        while len(a) > 2 and abs(a[-3] - a[-2]) <= abs(a[-2] - a[-1]):
            ampl = abs(a[-3] - a[-2])
            mean = (a[-3] + a[-2]) / 2;
            if len(a) == 3:
                del a[0]
                if ampl > 0:
                    ampl_mean.append((ampl, mean))
            elif len(a) > 3:
                del a[-3:-1]
                if ampl > 0:
                    ampl_mean.append((ampl, mean))
                    ampl_mean.append((ampl, mean))
    for index in range(len(a) - 1):
        ampl = abs(a[index] - a[index + 1])
        mean = (a[index] + a[index + 1]) / 2;
        if ampl > 0:
            ampl_mean.append((ampl, mean))
    return ampl_mean

# --------------------------------------------------------------------------------}
# --- Peak_trough.py
# --------------------------------------------------------------------------------{
# @cython.locals(BEGIN=cython.int, MINZO=cython.int, MAXZO=cython.int, ENDZO=cython.int, \
#                R=cython.int, L=cython.int, i=cython.int, p=cython.int, f=cython.int)
def peak_trough(x, R):  #cpdef np.ndarray[long,ndim=1] peak_trough(np.ndarray[long,ndim=1] x, int R):
    """
    Returns list of local maxima/minima.

    x: 1-dimensional numpy array containing signal
    R: Thresshold (minimum difference between succeeding min and max

    This routine is implemented directly as described in
    "Recommended Practices for Wind Turbine Testing - 3. Fatigue Loads", 2. edition 1990, Appendix A
    """

    BEGIN = 0
    MINZO = 1
    MAXZO = 2
    ENDZO = 3
    S = np.zeros(x.shape[0] + 1, dtype=int)

    L = x.shape[0]
    goto = BEGIN

    while 1:
        if goto == BEGIN:
            trough = x[0]
            peak = x[0]

            i = 0
            p = 1
            f = 0
            while goto == BEGIN:
                i += 1
                if i == L:
                    goto = ENDZO
                    continue
                else:
                    if x[i] > peak:
                        peak = x[i]
                        if peak - trough >= R:
                            S[p] = trough
                            goto = MAXZO
                            continue
                    elif x[i] < trough:
                        trough = x[i]
                        if peak - trough >= R:
                            S[p] = peak
                            goto = MINZO
                            continue

        elif goto == MINZO:
            f = -1

            while goto == MINZO:
                i += 1
                if i == L:
                    goto = ENDZO
                    continue
                else:
                    if x[i] < trough:
                        trough = x[i]
                    else:
                        if x[i] - trough >= R:
                            p += 1
                            S[p] = trough
                            peak = x[i]
                            goto = MAXZO
                            continue
        elif goto == MAXZO:
            f = 1
            while goto == MAXZO:
                i += 1
                if i == L:
                    goto = ENDZO
                    continue
                else:
                    if x[i] > peak:
                        peak = x[i]
                    else:
                        if peak - x[i] >= R:
                            p += 1
                            S[p] = peak
                            trough = x[i]
                            goto = MINZO
                            continue
        elif goto == ENDZO:

            n = p + 1
            if abs(f) == 1:
                if f == 1:
                    S[n] = peak
                else:
                    S[n] = trough
            else:
                S[n] = (trough + peak) / 2
            S = S[1:n + 1]
            return S


# --------------------------------------------------------------------------------}
# --- pair_range.py
# --------------------------------------------------------------------------------{
# @cython.locals(p=cython.int, q=cython.int, f=cython.int, flow=list, k=cython.int, n=cython.int, ptr=cython.int)
def pair_range_amplitude(x):  # cpdef pair_range(np.ndarray[long,ndim=1]  x):
    """
    Returns a list of half-cycle-amplitudes
    x: Peak-Trough sequence (integer list of local minima and maxima)

    This routine is implemented according to
    "Recommended Practices for Wind Turbine Testing - 3. Fatigue Loads", 2. edition 1990, Appendix A
    except that a list of half-cycle-amplitudes are returned instead of a from_level-to_level-matrix
    """

    x = x - np.min(x)
    k = np.max(x)
    n = x.shape[0]
    S = np.zeros(n + 1)

    #A = np.zeros(k+1)
    flow = []
    S[1] = x[0]
    ptr = 1
    p = 1
    q = 1
    f = 0
    # phase 1
    while True:
        p += 1
        q += 1

        # read
        S[p] = x[ptr]
        ptr += 1

        if q == n:
            f = 1
        while p >= 4:
            if (S[p - 2] > S[p - 3] and S[p - 1] >= S[p - 3] and S[p] >= S[p - 2]) \
                or\
                    (S[p - 2] < S[p - 3] and S[p - 1] <= S[p - 3] and S[p] <= S[p - 2]):
                ampl = abs(S[p - 2] - S[p - 1])
                # A[ampl]+=2 #Two half cycles
                flow.append(ampl)
                flow.append(ampl)
                S[p - 2] = S[p]

                p -= 2
            else:
                break

        if f == 0:
            pass
        else:
            break
    # phase 2
    q = 0
    while True:
        q += 1
        if p == q:
            break
        else:
            ampl = abs(S[q + 1] - S[q])
            # A[ampl]+=1
            flow.append(ampl)
    return flow





# @cython.locals(p=cython.int, q=cython.int, f=cython.int, flow=list, k=cython.int, n=cython.int, ptr=cython.int)
def pair_range_from_to(x):  # cpdef pair_range(np.ndarray[long,ndim=1]  x):
    """
    Returns a list of half-cycle-amplitudes
    x: Peak-Trough sequence (integer list of local minima and maxima)

    This routine is implemented according to
    "Recommended Practices for Wind Turbine Testing - 3. Fatigue Loads", 2. edition 1990, Appendix A
    except that a list of half-cycle-amplitudes are returned instead of a from_level-to_level-matrix
    """

    x = x - np.min(x)
    k = np.max(x)
    n = x.shape[0]
    S = np.zeros(n + 1)

    A = np.zeros((k + 1, k + 1))
    S[1] = x[0]
    ptr = 1
    p = 1
    q = 1
    f = 0
    # phase 1
    while True:
        p += 1
        q += 1

        # read
        S[p] = x[ptr]
        ptr += 1

        if q == n:
            f = 1
        while p >= 4:
            #print S[p - 3:p + 1]
            #print S[p - 2], ">", S[p - 3], ", ", S[p - 1], ">=", S[p - 3], ", ", S[p], ">=", S[p - 2], (S[p - 2] > S[p - 3] and S[p - 1] >= S[p - 3] and S[p] >= S[p - 2])
            #print S[p - 2], "<", S[p - 3], ", ", S[p - 1], "<=", S[p - 3], ", ", S[p], "<=", S[p - 2], (S[p - 2] < S[p - 3] and S[p - 1] <= S[p - 3] and S[p] <= S[p - 2])
            #print (S[p - 2] > S[p - 3] and S[p - 1] >= S[p - 3] and S[p] >= S[p - 2]) or (S[p - 2] < S[p - 3] and S[p - 1] <= S[p - 3] and S[p] <= S[p - 2])
            if (S[p - 2] > S[p - 3] and S[p - 1] >= S[p - 3] and S[p] >= S[p - 2]) or \
               (S[p - 2] < S[p - 3] and S[p - 1] <= S[p - 3] and S[p] <= S[p - 2]):
                A[S[p - 2], S[p - 1]] += 1
                A[S[p - 1], S[p - 2]] += 1
                S[p - 2] = S[p]
                p -= 2
            else:
                break

        if f == 1:
            break  # q==n
    # phase 2
    q = 0
    while True:
        q += 1
        if p == q:
            break
        else:
            #print S[q], "to", S[q + 1]
            A[S[q], S[q + 1]] += 1
    return A

# @cython.locals(p=cython.int, q=cython.int, f=cython.int, flow=list, k=cython.int, n=cython.int, ptr=cython.int)
def pair_range_amplitude_mean(x):  # cpdef pair_range(np.ndarray[long,ndim=1]  x):
    """
    Returns a list of half-cycle-amplitudes
    x: Peak-Trough sequence (integer list of local minima and maxima)

    This routine is implemented according to
    "Recommended Practices for Wind Turbine Testing - 3. Fatigue Loads", 2. edition 1990, Appendix A
    except that a list of half-cycle-amplitudes are returned instead of a from_level-to_level-matrix
    """

    x = x - np.min(x)
    k = np.max(x)
    n = x.shape[0]
    S = np.zeros(n + 1)
    ampl_mean = []
    A = np.zeros((k + 1, k + 1))
    S[1] = x[0]
    ptr = 1
    p = 1
    q = 1
    f = 0
    # phase 1
    while True:
        p += 1
        q += 1

                # read
        S[p] = x[ptr]
        ptr += 1

        if q == n:
            f = 1
        while p >= 4:
            if (S[p - 2] > S[p - 3] and S[p - 1] >= S[p - 3] and S[p] >= S[p - 2]) \
                or\
                    (S[p - 2] < S[p - 3] and S[p - 1] <= S[p - 3] and S[p] <= S[p - 2]):
                # Extract two intermediate half cycles
                ampl = abs(S[p - 2] - S[p - 1])
                mean = (S[p - 2] + S[p - 1]) / 2
                ampl_mean.append((ampl, mean))
                ampl_mean.append((ampl, mean))

                S[p - 2] = S[p]

                p -= 2
            else:
                break

        if f == 0:
            pass
        else:
            break
    # phase 2
    q = 0
    while True:
        q += 1
        if p == q:
            break
        else:
            ampl = abs(S[q + 1] - S[q])
            mean = (S[q + 1] + S[q]) / 2
            ampl_mean.append((ampl, mean))
    return ampl_mean


rainflow_func_dict = {'rainflow_windap':rainflow_windap, 'rainflow_astm':rainflow_astm}


# --------------------------------------------------------------------------------}
# --- Unittests
# --------------------------------------------------------------------------------{
import unittest

class TestFatigue(unittest.TestCase):

    def test_leq_1hz(self):
        """Simple test of wetb.fatigue.eq_load using a sine
        signal.
        """
        amplitude = 1
        m = 1
        point_per_deg = 100

        for amplitude in [1,2,3]:
            peak2peak = amplitude * 2
            # sine signal with 10 periods (20 peaks)
            nr_periods = 10
            time = np.linspace(0, nr_periods*2*np.pi, point_per_deg*180)
            neq = time[-1]
            # mean value of the signal shouldn't matter
            signal = amplitude * np.sin(time) + 5
            r_eq_1hz = eq_load(signal, no_bins=1, m=m, neq=neq)[0]
            r_eq_1hz_expected = 2*((nr_periods*amplitude**m)/neq)**(1/m)
            np.testing.assert_allclose(r_eq_1hz, r_eq_1hz_expected)

            # sine signal with 20 periods (40 peaks)
            nr_periods = 20
            time = np.linspace(0, nr_periods*2*np.pi, point_per_deg*180)
            neq = time[-1]
            # mean value of the signal shouldn't matter
            signal = amplitude * np.sin(time) + 9
            r_eq_1hz2 = eq_load(signal, no_bins=1, m=m, neq=neq)[0]
            r_eq_1hz_expected2 = 2*((nr_periods*amplitude**m)/neq)**(1/m)
            np.testing.assert_allclose(r_eq_1hz2, r_eq_1hz_expected2)

            # 1hz equivalent should be independent of the length of the signal
            np.testing.assert_allclose(r_eq_1hz, r_eq_1hz2)

    def test_rainflow_combi(self):
        # Signal with two frequencies and amplitudes
        amplitude = 1
        # peak2peak = amplitude * 2
        m = 1
        point_per_deg = 100

        nr_periods = 10
        time = np.linspace(0, nr_periods*2*np.pi, point_per_deg*180)

        signal = (amplitude*np.sin(time)) + 5 + (amplitude*0.2*np.cos(5*time))
        cycles, ampl_bin_mean, ampl_edges, mean_bin_mean, mean_edges = \
            cycle_matrix(signal, ampl_bins=10, mean_bins=5)

        cycles.sum()



    def test_astm1(self):

        signal = np.array([-2.0, 0.0, 1.0, 0.0, -3.0, 0.0, 5.0, 0.0, -1.0, 0.0, 3.0, 0.0, -4.0, 0.0, 4.0, 0.0, -2.0])

        ampl, mean = rainflow_astm(signal)
        np.testing.assert_array_equal(np.histogram2d(ampl, mean, [6, 4])[0], np.array([[ 0., 1., 0., 0.],
                                                                                                           [ 1., 0., 0., 2.],
                                                                                                           [ 0., 0., 0., 0.],
                                                                                                           [ 0., 0., 0., 1.],
                                                                                                           [ 0., 0., 0., 0.],
                                                                                                           [ 0., 0., 1., 2.]]))

    def test_windap1(self):
        signal = np.array([-2.0, 0.0, 1.0, 0.0, -3.0, 0.0, 5.0, 0.0, -1.0, 0.0, 3.0, 0.0, -4.0, 0.0, 4.0, 0.0, -2.0])
        ampl, mean = rainflow_windap(signal, 18, 2)
        np.testing.assert_array_equal(np.histogram2d(ampl, mean, [6, 4])[0], np.array([[ 0., 0., 1., 0.],
                                                                                       [ 1., 0., 0., 2.],
                                                                                       [ 0., 0., 0., 0.],
                                                                                       [ 0., 0., 0., 1.],
                                                                                       [ 0., 0., 0., 0.],
                                                                                       [ 0., 0., 2., 1.]]))

    def test_eq_load_basic(self):
        import numpy.testing
        signal1 = np.array([-2.0, 0.0, 1.0, 0.0, -3.0, 0.0, 5.0, 0.0, -1.0, 0.0, 3.0, 0.0, -4.0, 0.0, 4.0, 0.0, -2.0])
        try:
            M1=eq_load(signal1, no_bins=50, neq=[1, 17], m=[3, 4, 6], rainflow_func=rainflow_windap)
            doTest=True
        except FloatingPointError as e:
            doTest=False
            print('>>> Floating point error')
        M1_ref=np.array([[10.348414123746581, 9.635653414943068, 9.122399471334054], [4.024613313976801, 4.745357541147315, 5.68897815218057]])
        #M1_ref=np.array([[10.311095426959747, 9.5942535021382174, 9.0789213365013932],[4.010099657859783, 4.7249689509841746, 5.6618639965313005]])
        numpy.testing.assert_almost_equal(M1,M1_ref,decimal=5)
        #signal2 = signal1 * 1.1
        #         print (eq_load(signal1, no_bins=50, neq=17, rainflow_func=rainflow_windap))
        #         print (eq_load(signal1, no_bins=50, neq=17, rainflow_func=rainflow_astm))
        #         # equivalent load for default wohler slopes
        #         # Cycle matrix with 4 amplitude bins and 4 mean value bins
        #         print (cycle_matrix(signal1, 4, 4, rainflow_func=rainflow_windap))
        #         print (cycle_matrix(signal1, 4, 4, rainflow_func=rainflow_astm))
        #         # Cycle matrix where signal1 and signal2 contributes with 50% each
        #         print (cycle_matrix([(.5, signal1), (.5, signal2)], 4, 8, rainflow_func=rainflow_astm))


    def test_equivalent_load(self):
        """ Higher level interface """
        try:
            import fatpack
            hasFatpack=True
        except:
            hasFatpack=False
        dt = 0.1
        f0 = 1  ; 
        A  = 5  ; 
        t=np.arange(0,10,dt);
        y=A*np.sin(2*np.pi*f0*t)

        Leq = equivalent_load(t, y, m=10, bins=100, method='rainflow_windap')
        np.testing.assert_almost_equal(Leq, 9.4714702, 3)

        Leq = equivalent_load(t, y, m=1, bins=100, method='rainflow_windap')
        np.testing.assert_almost_equal(Leq, 9.4625320, 3)

        Leq = equivalent_load(t, y, m=4, bins=10, method='rainflow_windap')
        np.testing.assert_almost_equal(Leq, 9.420937, 3)


        if hasFatpack:
            Leq = equivalent_load(t, y, m=4, bins=10, method='fatpack', binStartAt0=False, meanBin=False)
            np.testing.assert_almost_equal(Leq, 9.584617089, 3)

            Leq = equivalent_load(t, y, m=4, bins=1, method='fatpack', binStartAt0=False, meanBin=False)
            np.testing.assert_almost_equal(Leq, 9.534491302, 3)



    def test_equivalent_load_sines(self):
        # Check analytical formulae for sine of various frequencies
        # See welib.tools.examples.Example_Fatigue.py
        try:
            import fatpack
            hasFatpack=True
        except:
            hasFatpack=False

        # --- Dependency on frequency
        m     = 2   # Wohler slope
        A     = 3   # Amplitude
        nT    = 100 # Number of periods
        nPerT = 100 # Number of points per period
        Teq   = 1  # Equivalent period [s]
        nBins = 10  # Number of bins

        vf =np.linspace(0.1,10,21)
        vT  = 1/vf
        T_max=np.max(vT*nT)
        vomega =vf*2*np.pi
        Leq1    = np.zeros_like(vomega)
        Leq2    = np.zeros_like(vomega)
        Leq_ref = np.zeros_like(vomega)
        for it, (T,omega) in enumerate(zip(vT,vomega)):
            # --- Option 1 - Same number of periods
            time = np.linspace(0, nT*T, nPerT*nT+1)
            signal = A * np.sin(omega*time) # Mean does not matter 
            T_all=time[-1]
            Leq1[it] = equivalent_load(time, signal, m=m, Teq=Teq, bins=nBins, method='rainflow_windap')
            if hasFatpack:
                Leq2[it] = equivalent_load(time, signal, m=m, Teq=Teq, bins=nBins, method='fatpack', binStartAt0=False)
        Leq_ref = 2*A*(vf*Teq)**(1/m)
        np.testing.assert_array_almost_equal(    Leq1/A, Leq_ref/A, 2)
        if hasFatpack:
            np.testing.assert_array_almost_equal(Leq2/A, Leq_ref/A, 2)
        #import matplotlib.pyplot as plt
        #fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
        #fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        #ax.plot(vf, Leq_ref/A,  'kd' , label ='Theory')
        #ax.plot(vf, Leq1   /A,   'o' , label ='Windap m={}'.format(m))
        #if hasFatpack:
        #    ax.plot(vf, Leq2/A,  'k.' , label ='Fatpack')
        #ax.legend()
        #plt.show()

    def test_equivalent_load_sines_sum(self):
        # --- Sum of two sinusoids
        try:
            import fatpack
            hasFatpack=True
        except:
            hasFatpack=False
        bs0 = True # Bin Start at 0
        m     = 2   # Wohler slope
        nPerT = 100 # Number of points per period
        Teq   = 1  # Equivalent period [s]
        nBins = 10 # Number of bins
        nT1   = 10 # Number of periods
        nT2   = 20 # Number of periods
        T1 = 10
        T2 = 5
        A1 = 3 # Amplitude
        A2 = 5 # Amplitude
        # --- Signals
        time1 = np.linspace(0, nT1*T1, nPerT*nT1+1)
        time2 = np.linspace(0, nT2*T2, nPerT*nT2+1)
        signal1 = A1 * np.sin(2*np.pi/T1*time1)
        signal2 = A2 * np.sin(2*np.pi/T2*time2)
        # --- Individual Leq
        #print('----------------- SIGNAL 1')
        DEL1 = (2*A1)**m * nT1/time1[-1]
        Leq_th = (DEL1)**(1/m)
        Leq1 = equivalent_load(time1, signal1, m=m, Teq=Teq, bins=nBins, method='rainflow_windap', binStartAt0=bs0)
        if hasFatpack:
            Leq2 = equivalent_load(time1, signal1, m=m, Teq=Teq, bins=nBins, method='fatpack', binStartAt0=bs0)
            np.testing.assert_array_almost_equal(Leq2, Leq_th, 3)
        np.testing.assert_array_almost_equal(Leq1, Leq_th, 1)
        #print('>>> Leq1   ',Leq1)
        #print('>>> Leq2   ',Leq2)
        #print('>>> Leq TH ',Leq_th)
        #print('----------------- SIGNAL 2')
        DEL2 = (2*A2)**m * nT2/time2[-1]
        Leq_th = (DEL2)**(1/m)
        Leq1 = equivalent_load(time2, signal2, m=m, Teq=Teq, bins=nBins, method='rainflow_windap', binStartAt0=bs0)
        if hasFatpack:
            Leq2 = equivalent_load(time2, signal2, m=m, Teq=Teq, bins=nBins, method='fatpack', binStartAt0=bs0)
            np.testing.assert_array_almost_equal(Leq2, Leq_th, 3)
        np.testing.assert_array_almost_equal(Leq1, Leq_th, 1)
        #print('>>> Leq1   ',Leq1)
        #print('>>> Leq2   ',Leq2)
        #print('>>> Leq TH ',Leq_th)
        # --- Concatenation
        #print('----------------- CONCATENATION')
        signal = np.concatenate((signal1, signal2))
        time   = np.concatenate((time1, time2+time1[-1]))
        T_all=time[-1]
        DEL1 = (2*A1)**m * nT1/T_all  
        DEL2 = (2*A2)**m * nT2/T_all  
        Leq_th = (DEL1+DEL2)**(1/m)
        Leq1 = equivalent_load(time, signal, m=m, Teq=Teq, bins=nBins, method='rainflow_windap', binStartAt0=bs0)
        if hasFatpack:
            Leq2 = equivalent_load(time, signal, m=m, Teq=Teq, bins=nBins, method='fatpack', binStartAt0=bs0)
            np.testing.assert_array_almost_equal(Leq2, Leq_th, 1)
        np.testing.assert_array_almost_equal(Leq1, Leq_th, 1)
        #print('>>> Leq1   ',Leq1)
        #print('>>> Leq2   ',Leq2)
        #print('>>> Leq TH ',Leq_th)



    def test_eqload_cornercases(self):
        try:
            import fatpack
            hasFatpack=True
        except:
            hasFatpack=False
        # Signal of length 1
        time=[0]; signal=[0]
        Leq= equivalent_load(time, signal, m=3, Teq=1, bins=100, method='rainflow_windap')
        np.testing.assert_equal(Leq, np.nan)
            
        # Datetime
        time= [np.datetime64('2023-10-01'), np.datetime64('2023-10-02')]
        signal= [0,1]
        Leq= equivalent_load(time, signal, m=3, Teq=1, bins=100, method='rainflow_windap')
        np.testing.assert_equal(Leq, np.nan)

        # Constant signal
        time =[0,1]
        signal =[1,1]
        Leq= equivalent_load(time, signal, m=3, Teq=1, bins=100, method='rainflow_windap')
        np.testing.assert_equal(Leq, np.nan)
        if hasFatpack:
            Leq= equivalent_load(time, signal, m=3, Teq=1, bins=100, method='fatpack')
            np.testing.assert_equal(Leq, np.nan)



if __name__ == '__main__':
    unittest.main()

