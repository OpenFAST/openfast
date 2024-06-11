#  Tools for spectral analysis of a real valued signal.
#
#  The functions in this file were adapted from the python package scipy according to the following license:
# 
# License: 
# Copyright   2001, 2002 Enthought, Inc.
# All rights reserved.
# 
# Copyright   2003-2013 SciPy Developers.
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
# 
#     Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#     Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
#     Neither the name of Enthought nor the names of the SciPy Developers may be used to endorse or promote products derived from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import numpy as np
import pandas as pd
from six import string_types

__all__  = ['fft_wrap','welch', 'psd', 'fft_amplitude']
__all__ += ['pwelch', 'csd', 'coherence']
__all__ += ['fnextpow2']
__all__ += ['hann','hamming','boxcar','general_hamming','get_window']
__all__ += ['TestSpectral']


# --------------------------------------------------------------------------------}
# --- FFT wrap
# --------------------------------------------------------------------------------{
def fft_wrap(t,y,dt=None, output_type='amplitude',averaging='None',averaging_window='hamming',detrend=False,nExp=None, nPerDecade=None):
    """ 
    Wrapper to compute FFT amplitude or power spectra, with averaging.
    INPUTS:
       output_type      : amplitude, PSD, f x PSD
       averaging : None, Welch, Binning
       averaging_window : Hamming, Hann, Rectangular
    OUTPUTS:
       frq: vector of frequencies
       Y  : Amplitude spectrum, PSD, or f * PSD
       Info: a dictionary of info values
    """

    # Formatting inputs
    output_type      = output_type.lower()
    averaging        = averaging.lower()
    averaging_window = averaging_window.lower()
    t = np.asarray(t)
    y = np.asarray(y)
    n0 = len(y) 
    nt = len(t) 
    if len(t)!=len(y):
        raise Exception('t and y should have the same length')
    y = y[~np.isnan(y)]
    n = len(y) 

    if dt is None:
        dtDelta0 = t[1]-t[0]
        # Hack to use a constant dt
        #dt = (np.max(t)-np.min(t))/(n0-1)
        dt = (t[-1]-t[0])/(n0-1)
        relDiff = abs(dtDelta0-dt)/dt*100
        #if dtDelta0 !=dt:
        if relDiff>0.01:
            print('[WARN] dt from tmax-tmin different from dt from t2-t1 {} {}'.format(dt, dtDelta0) )
    Fs = 1/dt
    if averaging =='none':
        frq, PSD, Info = psd(y, fs=Fs, detrend=detrend, return_onesided=True)
    elif averaging =='binning':
        frq, PSD, Info = psd_binned(y, fs=Fs, detrend=detrend, return_onesided=True, nPerDecade=nPerDecade)
    elif averaging=='welch':
        # --- Welch - PSD
        #overlap_frac=0.5
        #return fnextpow2(np.sqrt(len(x)/(1-overlap_frac)))
        nFFTAll=fnextpow2(n)
        if nExp is None:
            nExp=int(np.log(nFFTAll)/np.log(2))-1
        nPerSeg=2**nExp
        if nPerSeg>n:
            print('[WARN] Power of 2 value was too high and was reduced. Disable averaging to use the full spectrum.');
            nExp=int(np.log(nFFTAll)/np.log(2))-1
            nPerSeg=2**nExp
        if averaging_window=='hamming':
           window = hamming(nPerSeg, True)# True=Symmetric, like matlab
        elif averaging_window=='hann':
           window = hann(nPerSeg, True)
        elif averaging_window=='rectangular':
           window = boxcar(nPerSeg)
        else:
            raise Exception('Averaging window unknown {}'.format(averaging_window))
        frq, PSD, Info = pwelch(y, fs=Fs, window=window, detrend=detrend)
        Info.nExp = nExp
    else:
        raise Exception('Averaging method unknown {}'.format(averaging))

    # --- Formatting output
    if output_type=='amplitude':
        deltaf = frq[1]-frq[0]
        Y = np.sqrt(PSD*2*deltaf)
        # NOTE: the above should be the same as:Y=abs(Y[range(nhalf)])/n;Y[1:-1]=Y[1:-1]*2;
    elif output_type=='psd': # one sided
        Y = PSD
    elif output_type=='f x psd':
        Y = PSD*frq
    else:
        raise NotImplementedError('Contact developer')
    if detrend:
        frq= frq[1:]
        Y  = Y[1:]
    return frq, Y, Info



# --------------------------------------------------------------------------------}
# --- Spectral simple (averaging below) 
# --------------------------------------------------------------------------------{
def fft_amplitude(y, fs=1.0, detrend ='constant', return_onesided=True):
    """ Returns FFT amplitude of signal """
    frq, PSD, Info = psd(y, fs=fs, detrend=detrend, return_onesided=return_onesided)
    deltaf = frq[1]-frq[0]
    Y = np.sqrt(PSD*2*deltaf)
    return frq, Y, Info


def psd_binned(y, fs=1.0, nPerDecade=10, detrend ='constant', return_onesided=True):
    """ 
    Return PSD binned with nPoints per decade
    """
    # --- First  return regular PSD
    frq, PSD, Info = psd(y, fs=fs, detrend=detrend, return_onesided=return_onesided)

    add0=False
    if frq[0]==0:
        add0=True
        f0   = 0
        PSD0 = PSD[0]
        frq=frq[1:]
        PSD=PSD[1:]

    # -- Then bin per decase
    log_f = np.log10(frq)
    ndecades = np.ceil(log_f[-1] -log_f[0])
    xbins = np.linspace(log_f[0], log_f[-1], int(ndecades*nPerDecade))

    # Using Pandas to bin..
    df = pd.DataFrame(data=np.column_stack((log_f,PSD)), columns=['x','y'])
    xmid  = (xbins[:-1]+xbins[1:])/2
    df['Bin'] = pd.cut(df['x'], bins=xbins, labels=xmid ) # Adding a column that has bin attribute
    df2  = df.groupby('Bin', observed=False).mean()                     # Average by bin
    df2  = df2.reindex(xmid)
    log_f_bin = df2['x'].values
    PSD_bin = df2['y'].values
    frq2= 10**log_f_bin
    PSD2= PSD_bin
    if add0:
        frq2=np.concatenate(  ([f0  ], frq2)  )
        PSD2=np.concatenate(  ([PSD0], PSD2)  )
    b = ~np.isnan(frq2)
    frq2 = frq2[b]
    PSD2 = PSD2[b]

    #import matplotlib.pyplot as plt
    #fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    #fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    #ax.plot(log_f, PSD, label='')
    #ax.plot(log_f_bin, PSD_bin, 'o', label='')
    #for x in xbins:
    #    ax.axvline(x, ls=':', c=(0.5,0.5,0.5))
    #ax.set_xlabel('')
    #ax.set_ylabel('')
    #ax.legend()
    #plt.show()

    #Info.df    = frq[1]-frq[0]
    #Info.fMax  = frq[-1]
    #Info.LFreq = len(frq)
    #Info.LSeg  = len(Y)
    #Info.LWin  = len(Y)
    #Info.LOvlp = 0
    #Info.nFFT  = len(Y)
    #Info.nseg  = 1
    Info.nPerDecade  = nPerDecade
    Info.xbins  = xbins

    return frq2, PSD2, Info


def psd(y, fs=1.0, detrend ='constant', return_onesided=True):
    """ Perform PSD without averaging """
    if not return_onesided:
        raise NotImplementedError('Double sided todo')

    if detrend is None:
        detrend=False

    if detrend=='constant' or detrend==True:
        m=np.mean(y);
    else:
        m=0;

    n = len(y) 
    if n%2==0:
        nhalf = int(n/2+1)
    else:
        nhalf = int((n+1)/2)

    frq = np.arange(nhalf)*fs/n;
    Y   = np.fft.rfft(y-m) #Y = np.fft.fft(y) 
    PSD = abs(Y[range(nhalf)])**2 /(n*fs) # PSD
    PSD[1:-1] = PSD[1:-1]*2;
    class InfoClass():
        pass
    Info = InfoClass();
    Info.df    = frq[1]-frq[0]
    Info.fMax  = frq[-1]
    Info.LFreq = len(frq)
    Info.LSeg  = len(Y)
    Info.LWin  = len(Y)
    Info.LOvlp = 0
    Info.nFFT  = len(Y)
    Info.nseg  = 1
    return frq, PSD, Info


# --------------------------------------------------------------------------------}
# --- Windows 
# --------------------------------------------------------------------------------{
"""The suite of window functions."""
def fnextpow2(x):
    return 2**np.ceil( np.log(x)*0.99999999999/np.log(2));

def fDefaultWinLen(x,overlap_frac=0.5):
    return fnextpow2(np.sqrt(len(x)/(1-overlap_frac)))

def fDefaultWinLenMatlab(x):
    return  np.fix((len(x)-3)*2./9.)

def _len_guards(M):
    """Handle small or incorrect window lengths"""
    if int(M) != M or M < 0:
        raise ValueError('Window length M must be a non-negative integer')
    return M <= 1

def _extend(M, sym):
    """Extend window by 1 sample if needed for DFT-even symmetry"""
    if not sym:
        return M + 1, True
    else:
        return M, False

def _truncate(w, needed):
    """Truncate window by 1 sample if needed for DFT-even symmetry"""
    if needed:
        return w[:-1]
    else:
        return w

def general_cosine(M, a, sym=True):
    if _len_guards(M):
        return np.ones(M)
    M, needs_trunc = _extend(M, sym)

    fac = np.linspace(-np.pi, np.pi, M)
    w = np.zeros(M)
    for k in range(len(a)):
        w += a[k] * np.cos(k * fac)

    return _truncate(w, needs_trunc)


def boxcar(M, sym=True):
    """Return a boxcar or rectangular window.

    Also known as a rectangular window or Dirichlet window, this is equivalent
    to no window at all.
    """
    if _len_guards(M):
        return np.ones(M)
    M, needs_trunc = _extend(M, sym)

    w = np.ones(M, float)

    return _truncate(w, needs_trunc)

def hann(M, sym=True): # same as hanning(*args, **kwargs):
    return general_hamming(M, 0.5, sym)


def general_hamming(M, alpha, sym=True):
    r"""Return a generalized Hamming window.
    The generalized Hamming window is constructed by multiplying a rectangular
    window by one period of a cosine function [1]_.
    w(n) = \alpha - \left(1 - \alpha\right) \cos\left(\frac{2\pi{n}}{M-1}\right)
              \qquad 0 \leq n \leq M-1
    """
    return general_cosine(M, [alpha, 1. - alpha], sym)


def hamming(M, sym=True):
    r"""Return a Hamming window.
    The Hamming window is a taper formed by using a raised cosine with
    non-zero endpoints, optimized to minimize the nearest side lobe.
     w(n) = 0.54 - 0.46 \cos\left(\frac{2\pi{n}}{M-1}\right)
               \qquad 0 \leq n \leq M-1
    """
    return general_hamming(M, 0.54, sym)

_win_equiv_raw = {
    ('boxcar', 'box', 'ones', 'rect', 'rectangular'): (boxcar, False),
    ('hamming', 'hamm', 'ham'): (hamming, False),
    ('hanning', 'hann', 'han'): (hann, False),
}

# Fill dict with all valid window name strings
_win_equiv = {}
for k, v in _win_equiv_raw.items():
    for key in k:
        _win_equiv[key] = v[0]

# Keep track of which windows need additional parameters
_needs_param = set()
for k, v in _win_equiv_raw.items():
    if v[1]:
        _needs_param.update(k)


def get_window(window, Nx, fftbins=True):
    """
    Return a window.

    Parameters
    ----------
    window : string, float, or tuple
        The type of window to create. See below for more details.
    Nx : int
        The number of samples in the window.
    fftbins : bool, optional
        If True (default), create a "periodic" window, ready to use with
        `ifftshift` and be multiplied by the result of an FFT (see also
        `fftpack.fftfreq`).
        If False, create a "symmetric" window, for use in filter design.
    """
    sym = not fftbins
    try:
        beta = float(window)
    except (TypeError, ValueError):
        args = ()
        if isinstance(window, tuple):
            winstr = window[0]
            if len(window) > 1:
                args = window[1:]
        elif isinstance(window, string_types):
            if window in _needs_param:
                raise ValueError("The '" + window + "' window needs one or "
                                 "more parameters -- pass a tuple.")
            else:
                winstr = window
        else:
            raise ValueError("%s as window type is not supported." %
                             str(type(window)))

        try:
            winfunc = _win_equiv[winstr]
        except KeyError:
            raise ValueError("Unknown window type.")

        params = (Nx,) + args + (sym,)
    else:
        winfunc = kaiser
        params = (Nx, beta, sym)

    return winfunc(*params)






# --------------------------------------------------------------------------------}
# --- Helpers 
# --------------------------------------------------------------------------------{
def odd_ext(x, n, axis=-1):
    """
    Odd extension at the boundaries of an array
    Generate a new ndarray by making an odd extension of `x` along an axis.
    """
    if n < 1:
        return x
    if n > x.shape[axis] - 1:
        raise ValueError(("The extension length n (%d) is too big. " +
                         "It must not exceed x.shape[axis]-1, which is %d.")
                         % (n, x.shape[axis] - 1))
    left_end = axis_slice(x, start=0, stop=1, axis=axis)
    left_ext = axis_slice(x, start=n, stop=0, step=-1, axis=axis)
    right_end = axis_slice(x, start=-1, axis=axis)
    right_ext = axis_slice(x, start=-2, stop=-(n + 2), step=-1, axis=axis)
    ext = np.concatenate((2 * left_end - left_ext,
                          x,
                          2 * right_end - right_ext),
                         axis=axis)
    return ext


def even_ext(x, n, axis=-1):
    """
    Even extension at the boundaries of an array
    Generate a new ndarray by making an even extension of `x` along an axis.
    """
    if n < 1:
        return x
    if n > x.shape[axis] - 1:
        raise ValueError(("The extension length n (%d) is too big. " +
                         "It must not exceed x.shape[axis]-1, which is %d.")
                         % (n, x.shape[axis] - 1))
    left_ext = axis_slice(x, start=n, stop=0, step=-1, axis=axis)
    right_ext = axis_slice(x, start=-2, stop=-(n + 2), step=-1, axis=axis)
    ext = np.concatenate((left_ext,
                          x,
                          right_ext),
                         axis=axis)
    return ext


def const_ext(x, n, axis=-1):
    """
    Constant extension at the boundaries of an array
    Generate a new ndarray that is a constant extension of `x` along an axis.
    The extension repeats the values at the first and last element of
    the axis.
    """
    if n < 1:
        return x
    left_end = axis_slice(x, start=0, stop=1, axis=axis)
    ones_shape = [1] * x.ndim
    ones_shape[axis] = n
    ones = np.ones(ones_shape, dtype=x.dtype)
    left_ext = ones * left_end
    right_end = axis_slice(x, start=-1, axis=axis)
    right_ext = ones * right_end
    ext = np.concatenate((left_ext,
                          x,
                          right_ext),
                         axis=axis)
    return ext


def zero_ext(x, n, axis=-1):
    """
    Zero padding at the boundaries of an array
    Generate a new ndarray that is a zero padded extension of `x` along
    an axis.
    """
    if n < 1:
        return x
    zeros_shape = list(x.shape)
    zeros_shape[axis] = n
    zeros = np.zeros(zeros_shape, dtype=x.dtype)
    ext = np.concatenate((zeros, x, zeros), axis=axis)
    return ext

def signaltools_detrend(data, axis=-1, type='linear', bp=0):
    """
    Remove linear trend along axis from data.

    Parameters
    ----------
    data : array_like
        The input data.
    axis : int, optional
        The axis along which to detrend the data. By default this is the
        last axis (-1).
    type : {'linear', 'constant'}, optional
        The type of detrending. If ``type == 'linear'`` (default),
        the result of a linear least-squares fit to `data` is subtracted
        from `data`.
        If ``type == 'constant'``, only the mean of `data` is subtracted.
    bp : array_like of ints, optional
        A sequence of break points. If given, an individual linear fit is
        performed for each part of `data` between two break points.
        Break points are specified as indices into `data`.

    Returns
    -------
    ret : ndarray
        The detrended input data.
    """
    if type not in ['linear', 'l', 'constant', 'c']:
        raise ValueError("Trend type must be 'linear' or 'constant'.")
    data = np.asarray(data)
    dtype = data.dtype.char
    if dtype not in 'dfDF':
        dtype = 'd'
    if type in ['constant', 'c']:
        #print('Removing mean')
        ret = data - np.expand_dims(np.mean(data, axis), axis)
        return ret
    else:
        #print('Removing linear?')
        dshape = data.shape
        N = dshape[axis]
        bp = sort(unique(r_[0, bp, N]))
        if np.any(bp > N):
            raise ValueError("Breakpoints must be less than length "
                             "of data along given axis.")
        Nreg = len(bp) - 1
        # Restructure data so that axis is along first dimension and
        #  all other dimensions are collapsed into second dimension
        rnk = len(dshape)
        if axis < 0:
            axis = axis + rnk
        newdims = r_[axis, 0:axis, axis + 1:rnk]
        newdata = reshape(np.transpose(data, tuple(newdims)),
                          (N, _prod(dshape) // N))
        newdata = newdata.copy()  # make sure we have a copy
        if newdata.dtype.char not in 'dfDF':
            newdata = newdata.astype(dtype)
        # Find leastsq fit and remove it for each piece
        for m in range(Nreg):
            Npts = bp[m + 1] - bp[m]
            A = ones((Npts, 2), dtype)
            A[:, 0] = cast[dtype](np.arange(1, Npts + 1) * 1.0 / Npts)
            sl = slice(bp[m], bp[m + 1])
            coef, resids, rank, s = np.linalg.lstsq(A, newdata[sl])
            newdata[sl] = newdata[sl] - dot(A, coef)
        # Put data back in original shape.
        tdshape = take(dshape, newdims, 0)
        ret = np.reshape(newdata, tuple(tdshape))
        vals = list(range(1, rnk))
        olddims = vals[:axis] + [0] + vals[axis:]
        ret = np.transpose(ret, tuple(olddims))
        return ret



# --------------------------------------------------------------------------------}
# --- Spectral Averaging
# --------------------------------------------------------------------------------{
"""Tools for spectral analysis.  """

def welch(x, fs=1.0, window='hann', nperseg=None, noverlap=None, nfft=None,
          detrend='constant', return_onesided=True, scaling='density',
          axis=-1):
    """Interface identical to scipy.signal """

    if detrend==True:
        detrend='constant'

    freqs, Pxx = csd(x, x, fs, window, nperseg, noverlap, nfft, detrend, return_onesided, scaling, axis)
    return freqs, Pxx.real

#>>>>
def pwelch(x, window='hamming', noverlap=None, nfft=None, fs=1.0, nperseg=None, 
          detrend=False, return_onesided=True, scaling='density',
          axis=-1):
    r"""
    NOTE: interface and default options modified to match matlab's implementation
       >> detrend: default to False
       >> window : default to 'hamming'
       >> window: if an integer, use 'hamming(window, sym=True)'


    Estimate power spectral density using Welch's method.

    Welch's method [1]_ computes an estimate of the power spectral
    density by dividing the data into overlapping segments, computing a
    modified periodogram for each segment and averaging the
    periodograms.

    Parameters
    ----------
    x : array_like
        Time series of measurement values
    fs : float, optional
        Sampling frequency of the `x` time series. Defaults to 1.0.
    window : str or tuple or array_like, optional
        Desired window to use. If `window` is a string or tuple, it is
        passed to `get_window` to generate the window values, which are
        DFT-even by default. See `get_window` for a list of windows and
        required parameters. If `window` is array_like it will be used
        directly as the window and its length must be nperseg. Defaults
        to a Hann window.
    nperseg : int, optional
        Length of each segment. Defaults to None, but if window is str or
        tuple, is set to 256, and if window is array_like, is set to the
        length of the window.
    noverlap : int, optional
        Number of points to overlap between segments. If `None`,
        ``noverlap = nperseg // 2``. Defaults to `None`.
    nfft : int, optional
        Length of the FFT used, if a zero padded FFT is desired. If
        `None`, the FFT length is `nperseg`. Defaults to `None`.
    detrend : str or function or `False`, optional
        Specifies how to detrend each segment. If `detrend` is a
        string, it is passed as the `type` argument to the `detrend`
        function. If it is a function, it takes a segment and returns a
        detrended segment. If `detrend` is `False`, no detrending is
        done. Defaults to 'constant'.
    return_onesided : bool, optional
        If `True`, return a one-sided spectrum for real data. If
        `False` return a two-sided spectrum. Note that for complex
        data, a two-sided spectrum is always returned.
    scaling : { 'density', 'spectrum' }, optional
        Selects between computing the power spectral density ('density')
        where `Pxx` has units of V**2/Hz and computing the power
        spectrum ('spectrum') where `Pxx` has units of V**2, if `x`
        is measured in V and `fs` is measured in Hz. Defaults to
        'density'
    axis : int, optional
        Axis along which the periodogram is computed; the default is
        over the last axis (i.e. ``axis=-1``).

    Returns
    -------
    f : ndarray
        Array of sample frequencies.
    Pxx : ndarray
        Power spectral density or power spectrum of x.

    See Also
    --------
    periodogram: Simple, optionally modified periodogram
    lombscargle: Lomb-Scargle periodogram for unevenly sampled data

    Notes
    -----
    An appropriate amount of overlap will depend on the choice of window
    and on your requirements. For the default Hann window an overlap of
    50% is a reasonable trade off between accurately estimating the
    signal power, while not over counting any of the data. Narrower
    windows may require a larger overlap.

    If `noverlap` is 0, this method is equivalent to Bartlett's method
    [2]_.

    .. versionadded:: 0.12.0

    References
    ----------
    .. [1] P. Welch, "The use of the fast Fourier transform for the
           estimation of power spectra: A method based on time averaging
           over short, modified periodograms", IEEE Trans. Audio
           Electroacoust. vol. 15, pp. 70-73, 1967.
    .. [2] M.S. Bartlett, "Periodogram Analysis and Continuous Spectra",
           Biometrika, vol. 37, pp. 1-16, 1950.

    """
    import math
    def fnextpow2(x):
        return 2**math.ceil( math.log(x)*0.99999999999/math.log(2));

    # MANU >>> CHANGE OF DEFAULT OPTIONS
    # MANU - If a length is provided use symmetric hamming window
    if type(window)==int:
        window=hamming(window, True) 
    # MANU - do not use 256 as default
    if isinstance(window, string_types) or isinstance(window, tuple):
        if nperseg is None:
            if noverlap is None:
                overlap_frac=0.5
            elif noverlap == 0:
                overlap_frac=0
            else:
                raise NotImplementedError('TODO noverlap set but not nperseg')
            #nperseg = 256  # then change to default
            nperseg=fnextpow2(math.sqrt(x.shape[-1]/(1-overlap_frac)));

    # MANU accepting true as detrend
    if detrend==True:
        detrend='constant'

    freqs, Pxx, Info = csd(x, x, fs, window, nperseg, noverlap, nfft, detrend,
                     return_onesided, scaling, axis, returnInfo=True)

    return freqs, Pxx.real, Info


def csd(x, y, fs=1.0, window='hann', nperseg=None, noverlap=None, nfft=None,
        detrend='constant', return_onesided=True, scaling='density', axis=-1,
        returnInfo=False
        ):
    r"""
    Estimate the cross power spectral density, Pxy, using Welch's
    method.
    """

    freqs, _, Pxy, Info = _spectral_helper(x, y, fs, window, nperseg, noverlap, nfft,
                                     detrend, return_onesided, scaling, axis,
                                     mode='psd')

    # Average over windows.
    if len(Pxy.shape) >= 2 and Pxy.size > 0:
        if Pxy.shape[-1] > 1:
            Pxy = Pxy.mean(axis=-1)
        else:
            Pxy = np.reshape(Pxy, Pxy.shape[:-1])

    if returnInfo:
        return freqs, Pxy, Info
    else:
        return freqs, Pxy



def coherence(x, y, fs=1.0, window='hann', nperseg=None, noverlap=None,
              nfft=None, detrend='constant', axis=-1):
    r"""
    Estimate the magnitude squared coherence estimate, Cxy, of
    discrete-time signals X and Y using Welch's method.

    ``Cxy = abs(Pxy)**2/(Pxx*Pyy)``, where `Pxx` and `Pyy` are power
    spectral density estimates of X and Y, and `Pxy` is the cross
    spectral density estimate of X and Y.
    """

    freqs, Pxx, Infoxx = welch(x, fs, window, nperseg, noverlap, nfft, detrend, axis=axis)
    _, Pyy, Infoyy     = welch(y, fs, window, nperseg, noverlap, nfft, detrend, axis=axis)
    _, Pxy, Infoxy     = csd(x, y, fs, window, nperseg, noverlap, nfft, detrend, axis=axis, returnInfo=True)

    Cxy = np.abs(Pxy)**2 / Pxx / Pyy

    return freqs, Cxy, Infoxx


def _spectral_helper(x, y, fs=1.0, window='hann', nperseg=None, noverlap=None,
                     nfft=None, detrend='constant', return_onesided=True,
                     scaling='spectrum', axis=-1, mode='psd', boundary=None,
                     padded=False):
    """ Calculate various forms of windowed FFTs for PSD, CSD, etc.  """
    if mode not in ['psd', 'stft']:
        raise ValueError("Unknown value for mode %s, must be one of: "
                         "{'psd', 'stft'}" % mode)
    




    boundary_funcs = {'even': even_ext,
                      'odd': odd_ext,
                      'constant': const_ext,
                      'zeros': zero_ext,
                      None: None}

    if boundary not in boundary_funcs:
        raise ValueError("Unknown boundary option '{0}', must be one of: {1}"
                          .format(boundary, list(boundary_funcs.keys())))

    # If x and y are the same object we can save ourselves some computation.
    same_data = y is x

    if not same_data and mode != 'psd':
        raise ValueError("x and y must be equal if mode is 'stft'")

    axis = int(axis)

    # Ensure we have np.arrays, get outdtype
    x = np.asarray(x)
    if not same_data:
        y = np.asarray(y)
        outdtype = np.result_type(x, y, np.complex64)
    else:
        outdtype = np.result_type(x, np.complex64)

    if not same_data:
        # Check if we can broadcast the outer axes together
        xouter = list(x.shape)
        youter = list(y.shape)
        xouter.pop(axis)
        youter.pop(axis)
        try:
            outershape = np.broadcast(np.empty(xouter), np.empty(youter)).shape
        except ValueError:
            raise ValueError('x and y cannot be broadcast together.')

    if same_data:
        if x.size == 0:
            return np.empty(x.shape), np.empty(x.shape), np.empty(x.shape)
    else:
        if x.size == 0 or y.size == 0:
            outshape = outershape + (min([x.shape[axis], y.shape[axis]]),)
            emptyout = np.rollaxis(np.empty(outshape), -1, axis)
            return emptyout, emptyout, emptyout

    if x.ndim > 1:
        if axis != -1:
            x = np.rollaxis(x, axis, len(x.shape))
            if not same_data and y.ndim > 1:
                y = np.rollaxis(y, axis, len(y.shape))

    # Check if x and y are the same length, zero-pad if necessary
    if not same_data:
        if x.shape[-1] != y.shape[-1]:
            if x.shape[-1] < y.shape[-1]:
                pad_shape = list(x.shape)
                pad_shape[-1] = y.shape[-1] - x.shape[-1]
                x = np.concatenate((x, np.zeros(pad_shape)), -1)
            else:
                pad_shape = list(y.shape)
                pad_shape[-1] = x.shape[-1] - y.shape[-1]
                y = np.concatenate((y, np.zeros(pad_shape)), -1)

    if nperseg is not None:  # if specified by user
        nperseg = int(nperseg)
        if nperseg < 1:
            raise ValueError('nperseg must be a positive integer')

    # parse window; if array like, then set nperseg = win.shape
    win, nperseg = _triage_segments(window, nperseg,input_length=x.shape[-1])

    if nfft is None:
        nfft = nperseg
    elif nfft < nperseg:
        raise ValueError('nfft must be greater than or equal to nperseg.')
    else:
        nfft = int(nfft)

    if noverlap is None:
        noverlap = nperseg//2
    else:
        noverlap = int(noverlap)
    if noverlap >= nperseg:
        raise ValueError('noverlap must be less than nperseg.')
    nstep = nperseg - noverlap

    # Padding occurs after boundary extension, so that the extended signal ends
    # in zeros, instead of introducing an impulse at the end.
    # I.e. if x = [..., 3, 2]
    # extend then pad -> [..., 3, 2, 2, 3, 0, 0, 0]
    # pad then extend -> [..., 3, 2, 0, 0, 0, 2, 3]

    if boundary is not None:
        ext_func = boundary_funcs[boundary]
        x = ext_func(x, nperseg//2, axis=-1)
        if not same_data:
            y = ext_func(y, nperseg//2, axis=-1)

    if padded:
        # Pad to integer number of windowed segments
        # I.e make x.shape[-1] = nperseg + (nseg-1)*nstep, with integer nseg
        nadd = (-(x.shape[-1]-nperseg) % nstep) % nperseg
        zeros_shape = list(x.shape[:-1]) + [nadd]
        x = np.concatenate((x, np.zeros(zeros_shape)), axis=-1)
        if not same_data:
            zeros_shape = list(y.shape[:-1]) + [nadd]
            y = np.concatenate((y, np.zeros(zeros_shape)), axis=-1)

    # Handle detrending and window functions
    if not detrend:
        def detrend_func(d):
            return d
    elif not hasattr(detrend, '__call__'):
        def detrend_func(d):
            return signaltools_detrend(d, type=detrend, axis=-1)
    elif axis != -1:
        # Wrap this function so that it receives a shape that it could
        # reasonably expect to receive.
        def detrend_func(d):
            d = np.rollaxis(d, -1, axis)
            d = detrend(d)
            return np.rollaxis(d, axis, len(d.shape))
    else:
        detrend_func = detrend

    if np.result_type(win,np.complex64) != outdtype:
        win = win.astype(outdtype)

    if scaling == 'density':
        scale = 1.0 / (fs * (win*win).sum())
    elif scaling == 'spectrum':
        scale = 1.0 / win.sum()**2
    else:
        raise ValueError('Unknown scaling: %r' % scaling)

    if mode == 'stft':
        scale = np.sqrt(scale)

    if return_onesided:
        if np.iscomplexobj(x):
            sides = 'twosided'
            #warnings.warn('Input data is complex, switching to ' 'return_onesided=False')
        else:
            sides = 'onesided'
            if not same_data:
                if np.iscomplexobj(y):
                    sides = 'twosided'
                    #warnings.warn('Input data is complex, switching to return_onesided=False')
    else:
        sides = 'twosided'

    if sides == 'twosided':
        raise Exception('NOT IMPLEMENTED')
         #freqs = fftpack.fftfreq(nfft, 1/fs)
    elif sides == 'onesided':
        freqs = np.fft.rfftfreq(nfft, 1/fs)

    # Perform the windowed FFTs
    result = _fft_helper(x, win, detrend_func, nperseg, noverlap, nfft, sides)

    if not same_data:
        # All the same operations on the y data
        result_y = _fft_helper(y, win, detrend_func, nperseg, noverlap, nfft,
                               sides)
        result = np.conjugate(result) * result_y
    elif mode == 'psd':
        result = np.conjugate(result) * result

    result *= scale
    if sides == 'onesided' and mode == 'psd':
        if nfft % 2:
            result[..., 1:] *= 2
        else:
            # Last point is unpaired Nyquist freq point, don't double
            result[..., 1:-1] *= 2

    time = np.arange(nperseg/2, x.shape[-1] - nperseg/2 + 1,
                     nperseg - noverlap)/float(fs)
    if boundary is not None:
        time -= (nperseg/2) / fs

    result = result.astype(outdtype)

    # All imaginary parts are zero anyways
    if same_data and mode != 'stft':
        result = result.real

    # Output is going to have new last axis for time/window index, so a
    # negative axis index shifts down one
    if axis < 0:
        axis -= 1

    # Roll frequency axis back to axis where the data came from
    result = np.rollaxis(result, -1, axis)

    # TODO
    class InfoClass():
        pass
    Info = InfoClass();
    Info.df=freqs[1]-freqs[0]
    Info.fMax=freqs[-1]
    Info.LFreq=len(freqs)
    Info.LSeg=nperseg
    Info.LWin=len(win)
    Info.LOvlp=noverlap
    Info.nFFT=nfft
    Info.nseg=-1
    #print('df:{:.3f} - fm:{:.2f} - nseg:{} - Lf:{:5d} - Lseg:{:5d} - Lwin:{:5d} - Lovlp:{:5d} - Nfft:{:5d} - Lsig:{}'.format(freqs[1]-freqs[0],freqs[-1],-1,len(freqs),nperseg,len(win),noverlap,nfft,x.shape[-1]))
    return freqs, time, result, Info


def _fft_helper(x, win, detrend_func, nperseg, noverlap, nfft, sides):
    """ Calculate windowed FFT """
    # Created strided array of data segments
    if nperseg == 1 and noverlap == 0:
        result = x[..., np.newaxis]
    else:
        # http://stackoverflow.com/a/5568169
        step = nperseg - noverlap
        shape = x.shape[:-1]+((x.shape[-1]-noverlap)//step, nperseg)
        strides = x.strides[:-1]+(step*x.strides[-1], x.strides[-1])
        result = np.lib.stride_tricks.as_strided(x, shape=shape,
                                                 strides=strides)

    # Detrend each data segment individually
    result = detrend_func(result)

    # Apply window by multiplication
    result = win * result

    # Perform the fft. Acts on last axis by default. Zero-pads automatically
    if sides == 'twosided':
        raise Exception('NOT IMPLEMENTED')
        #func = fftpack.fft
    else:
        result = result.real
        func = np.fft.rfft
    result = func(result, n=nfft)

    return result

def _triage_segments(window, nperseg,input_length):
    """
    Parses window and nperseg arguments for spectrogram and _spectral_helper.
    This is a helper function, not meant to be called externally.
    """

    #parse window; if array like, then set nperseg = win.shape
    if isinstance(window, string_types) or isinstance(window, tuple):
        # if nperseg not specified
        if nperseg is None:
            nperseg = 256  # then change to default
        if nperseg > input_length:
            print('nperseg = {0:d} is greater than input length '
                              ' = {1:d}, using nperseg = {1:d}'
                              .format(nperseg, input_length))
            nperseg = input_length
        win = get_window(window, nperseg)
    else:
        win = np.asarray(window)
        if len(win.shape) != 1:
            raise ValueError('window must be 1-D')
        if input_length < win.shape[-1]:
            raise ValueError('window is longer than input signal')
        if nperseg is None:
            nperseg = win.shape[0]
        elif nperseg is not None:
            if nperseg != win.shape[0]:
                raise ValueError("value specified for nperseg is different from"
                                 " length of window")

    return win, nperseg




# --------------------------------------------------------------------------------
# --- Simple implementations to figure out the math
# --------------------------------------------------------------------------------
def DFT(x, method='vectorized'):
    """
    Calculate the Discrete Fourier Transform  (DFT) of real signal x

    Definition:
      for k in 0..N-1 (but defined for k in ZZ, see below for index)

          X_k = sum_n=0^{N-1} x_n  e^{-i 2 pi k n /N}

              = sum_n=0^{N-1} x_n [ cos( 2 pi k n /N) - i sin( 2 pi k n /N)

          Xk are complex numbers. The amplitude/phase are:
              A = |X_k|/N 
              phi = atan2(Im(Xk) / Re(Xk)
    
    Indices:
        The DFT creates a periodic signal of period N
            X[k]  = X[k+N]
            X[-k] = X[N-k]

        Therefore, any set of successive indices could be used.
        For instance, with N=4: [0,1,2,3], [-1,0,1,2], [-2,-1,0,1] (canonical one)
                              [0,..N/2-1],             [-N/2, ..N/2-1]

        If N is even
           0         is the 0th frequency (mean)
           1.....N/2-1  terms corresponds to positive frequencies
           N/2.....N-1  terms corresponds to negative frequencies
        If N is odd
           0         is the 0th frequency (mean)
           1.....(N-1)/2 terms corresponds to positive frequencies
           (N+1)/2..N-1  terms corresponds to negative frequencies

    Frequencies convention: (see np.fft.fftfreq and DFT_freq)
        f = [0, 1, ...,   n/2-1,     -n/2, ..., -1] / (dt*n)   if n is even
        f = [0, 1, ..., (n-1)/2, -(n-1)/2, ..., -1] / (dt*n)   if n is odd

    NOTE: when n is even you could chose to go to +n/2 and start at -n/2+1
          The Python convention goes to n/2-1 and start at -n/2. 

    Properties:
      - if x is a real signal
            X[-k] = X*[k]   (=X[N-k])
      - Parseval's theorem: sum |x_n|^2 = sum |X_k|^2   (energy conservation)


    """
    N = len(x)

    if method=='naive':
        X = np.zeros_like(x, dtype=complex)
        for k in np.arange(N):
            for n in np.arange(N):
                X[k] += x[n] * np.exp(-1j * 2*np.pi * k * n / N)

    elif method=='vectorized':
        n = np.arange(N)
        k = n.reshape((N, 1)) # k*n will be of shape (N x N)
        e = np.exp(-2j * np.pi * k * n / N)
        X = np.dot(e, x)
    elif method=='fft':
        X = np.fft.fft(x)
    elif method=='fft_py':
        X = recursive_fft(x)
    else:
        raise NotImplementedError()
    
    return X

def IDFT(X, method='vectorized'):
    """
    Calculate the Inverse Discrete Fourier Transform  (IDFT) of complex coefficients X

    The transformation DFT->IDFT is fully reversible

    Definition:
       for n in 0..N-1:

          x_n = 1/N  sum_k=0^{N-1} X_k  e^{i 2 pi k n/N}
              = 1/N  sum_k=0^{N-1} X_k [ cos( 2 pi k n/N) + i sin( 2 pi k n/N)
              = 1/N  sum_k=0^{N-1} A_k [ cos( 2 pi k n/N + phi_k) + i sin( 2 pi k n/N + phi_k)

        Xk are complex numbers, that can be written X[k] = A[k] e^{j phi[k]} therefore.

    Properties:
      - if the "X" given as input come from a DFT, then the coefficients are periodic with period N
        Therefore
            X[-k] = X[N-k], 
            X[k]  = X[N+k]
        and therefore (see the discussion in the documentation of DFT), the summation from
            k=0 to N-1 can be interpreted as a summation over any other indices set of length N. 
      - if "X" comes for the DFT of x, where x is a real signal, then:
            X[-k] = X*[k] (=X[N-k])

      - a converse is that, if X has conjugate symmetry (X[k]=X*[N-k]), then the IDFT will be real:

          x_n = 1/N  sum_k=0^{N-1}              A_k cos( 2 pi k n/N + phi_k)
                1/N  sum_k={-(N-1)/2}^{(N-1)/2} A_k cos( 2 pi k n/N + phi_k)

          But remember that the A_k and phi_k need to satisfy the conjugate symmetry, so they are not
          fully independent.
          If we want x to be the sum over "N0" independent components, then we need to do the IDFT 
          of a spectrum "X" of length 2N0-1.

    Indices and frequency (python convention):
      (see np.fft.fftfreq and DFT_freq)
            f = [0, 1, ...,   n/2-1,     -n/2, ..., -1] / (dt*n)   if n is even
            f = [0, 1, ..., (n-1)/2, -(n-1)/2, ..., -1] / (dt*n)   if n is odd
    
          When n is even, we lack symmetry of frequency, so it can potentially
          make sense to enforce that the X[-n/2] component is 0 when generating 
          a signal with IDFT



    """
    N = len(X)

    if method in ['naive', 'manual', 'sum']:
        x = np.zeros_like(X, dtype=complex)
        for k in np.arange(N):
            for n in np.arange(N):
                x[k] += X[n] * np.exp(1j * 2*np.pi * k * n / N)
        x = x/N

    elif method=='vectorized':
        n = np.arange(N)
        k = n.reshape((N, 1)) # k*n will be of shape (N x N)
        e = np.exp(2j * np.pi * k * n / N)
        x = np.dot(e, X) / N

    elif method=='ifft':
        x = np.fft.ifft(X)

    #elif method=='ifft_py':
    #    x = IFFT(X)
    else:
        raise NotImplementedError('IDFT: Method {}'.format(method))
    
    x = np.real_if_close(x)
    
    return x

def DFT_freq(time=None, N=None, T=None, doublesided=True):
    """ Returns the frequencies corresponding to a time vector `time`. 
      The signal "x" and "time" are assumed to have the same length
    INPUTS:
      - time: 1d array of time
      OR
      - N: number of time values
      - T: time length of signal
    """
    if time is not None:
        N = len(time)
        T = time[-1]-time[0]
    dt = T/(N-1)
    df = 1/(dt*N)
    nhalf_pos, nhalf_neg = nhalf_fft(N)
    if doublesided:
        freq_pos = np.arange(nhalf_pos+1)*df
        freq_neg = np.arange(nhalf_neg,0)*df
        freq = np.concatenate((freq_pos, freq_neg))
        assert(len(freq) == N)
    else:
        # single sided
        fMax = nhalf_pos * df
        #freq = np.arange(0, fMax+df/2, df)
        freq = np.arange(nhalf_pos+1)*df
    return freq

def IDFT_time(freq=None, doublesided=True):
    """ Returns the time vector corresponding to a frequency vector `freq`. 

    If doublesided is True
          The signal "x" , "time" and freq are assumed to have the same length
    Note: might lead to some inaccuracies, just use for double checking!

    INPUTS:
      - freq: 1d array of time
    """
    if doublesided:
        N = len(freq)
        time = freq*0
        if np.mod(N,2)==0:
            nhalf=int(N/2)-1
        else:
            nhalf=int((N-1)/2)
        fMax = freq[nhalf] 
        df = (fMax-0)/(nhalf)
        dt = 1/(df*N)
        tMax= (N-1)*dt
        #time = np.arange(0,(N-1)*dt+dt/2, dt)
        time = np.linspace(0,tMax, N)
    else:
        raise NotImplementedError()
    return time

def recursive_fft(x):
    """
    A recursive implementation of the 1D Cooley-Tukey FFT

    Returns the same as DFT (see documentation)

    Input should have a length of power of 2. 
    Reference: Kong, Siauw, Bayen - Python Numerical Methods
    """
    N = len(x)
    if not is_power_of_two(N):
        raise Exception('Recursive FFT requires a power of 2')

    
    if N == 1:
        return x
    else:
        X_even = recursive_fft(x[::2])
        X_odd  = recursive_fft(x[1::2])
        factor = np.exp(-2j*np.pi*np.arange(N)/ N)
        X = np.concatenate([X_even+factor[:int(N/2)]*X_odd, X_even+factor[int(N/2):]*X_odd])
        return X

def nhalf_fft(N):
    """ 
    Follows the convention of fftfreq
      fmax = f[nhalf_pos] = nhalf_pos*df  (fftfreq convention)

      fpos = f[:nhalf_pos+1]
      fneg = f[nhalf_pos+1:]

    """
    if N%2 ==0:
        nhalf_pos =   int(N/2)-1
        nhalf_neg =  -int(N/2)
    else:
        nhalf_pos = int((N-1)/2)
        nhalf_neg = -nhalf_pos
    return nhalf_pos, nhalf_neg


def check_DFT_real(X):
    """ Check that signal X is the DFT of a real signal
    and that therefore IDFT(X) will return a real signal.
    For this to be the case, we need conjugate symmetry:
       X[k] = X[N-k]* 
    """
    from welib.tools.spectral import nhalf_fft
    N = len(X)
    nh, _ = nhalf_fft(N)
    Xpos = X[1:nh+1] # we dont take the DC component [0]
    Xneg = np.flipud(X[nh+1:]) # might contain one more frequency than the pos part

    if np.mod(N,2)==0:
        # We have one extra negative frequency, we check that X is zero there and remove the value.
        if Xneg[-1]!=0:
            raise Exception('check_DFT_real: Component {} (first negative frequency) is {} instead of zero, but it should be zero if N is even.'.format(nh+1, X[nh+1]))
        Xneg = Xneg[:-1]

    notConjugate = Xpos-np.conjugate(Xneg)!=0
    if np.any(notConjugate): 
        nNotConjugate=sum(notConjugate)
        I = np.where(notConjugate)[0][:3] + 1 # +1 for DC component that was removed
        raise Exception('check_DFT_real: {}/{} values of the spectrum are not complex conjugate of there symmetric frequency counterpart. See for instance indices: {}'.format(nNotConjugate, nh, I))

def double_sided_DFT_real(X1, N=None):
    """ 
    Take a single sided part of a DFT (X1) and make it double sided signal X, of length N, 
    ensuring that the IDFT of X will be real.
    This is done by ensuring conjugate symmetry:
       X[k] = X[N-k]* 
    For N even, the first negative frequency component is set to 0 because it has no positive counterpart.

    Calling check_DFT_real(X) should return no Exception.

    INPUTS:
      - X1: array of complex values of length N1
      - N: required length of the output array (2N1-1 or 2N1)
    OUTPUTS:
      - X: double sided spectrum:
            [X1 flip(X1*[1:]) ]
          or
            [X1 [0] flip(X1*[1:]) ]
    """
    if N is None:
        N=2*len(X1)-1  # we make it an odd number to ensure symmetry of frequency
    else:
        if N not in [2*len(X1)-1, 2*len(X1), 2*len(X1)-2]:
            raise Exception('N should be twice the length of the single sided spectrum, or one less.')

    if N % 2 ==0:
        # Even number
        if N == 2*len(X1)-2:
            # rfftfreq
            # TODO, there look into irfft to see the convention
            X = np.concatenate((X1[:-1], [0], np.flipud(np.conjugate(X1[1:-1]))))
        else:
            X = np.concatenate((X1, [0], np.flipud(np.conjugate(X1[1:]))))
    else:
        X = np.concatenate((X1, np.flipud(np.conjugate(X1[1:]))))
    return X


# --------------------------------------------------------------------------------}
# --- Helper functions
# --------------------------------------------------------------------------------{
def is_power_of_two(n):
    """ Uses bit manipulation to figure out if an integer is a power of two"""
    return (n != 0) and (n & (n-1) == 0)

def sinesum(time, As, freqs):
    x =np.zeros_like(time)
    for ai,fi in zip(As, freqs):
        x += ai*np.sin(2*np.pi*fi*time)
    return x


# --------------------------------------------------------------------------------}
# --- Unittests
# --------------------------------------------------------------------------------{
import unittest

class TestSpectral(unittest.TestCase):

    def default_signal(self, time, mean=0):
        freqs=[1,4,7  ] # [Hz]
        As   =[3,1,1/2] # [misc]
        x = sinesum(time, As, freqs) + mean
        return x

    def compare_with_npfft(self, time, x):
        # Compare lowlevels functions with npfft
        # Useful to make sure the basic math is correct
        N    = len(time)
        dt   = (time[-1]-time[0])/(N-1)
        tMax = time[-1]

        # --- Test frequency, dt/df/N-relationships
        f_ref                = np.fft.fftfreq(N, dt)
        nhalf_pos, nhalf_neg = nhalf_fft(N)
        fhalf                = DFT_freq(time, doublesided = False)
        freq                 = DFT_freq(time, doublesided = True)
        df                   = freq[1]-freq[0]
        fmax                 = fhalf[-1]

        np.testing.assert_almost_equal(fhalf      , f_ref[:nhalf_pos+1], 10)
        np.testing.assert_almost_equal(fhalf[-1]  ,  np.max(f_ref), 10)
        np.testing.assert_almost_equal( 1/(dt*df), N)
        np.testing.assert_almost_equal(freq     , f_ref, 10)
        if N%2 == 0:
            np.testing.assert_almost_equal(2*fmax/df, N-2 , 10)
        else:
            np.testing.assert_almost_equal(2*fmax/df, N-1 , 10)

        # --- Test DFT methods
        X0 = DFT(x, method='fft')
        X1 = DFT(x, method='naive')
        X2 = DFT(x, method='vectorized')

        np.testing.assert_almost_equal(X1, X0, 10)
        np.testing.assert_almost_equal(X2, X0, 10)
        if is_power_of_two(N):
            X3 = DFT(x, method='fft_py')
            np.testing.assert_almost_equal(X3, X0, 10)

        # --- Test IDFT methods
        x_back0 = IDFT(X0, method='ifft')
        x_back1 = IDFT(X0, method='naive')
        x_back2 = IDFT(X0, method='vectorized')
        np.testing.assert_almost_equal(x_back1, x_back0, 10)
        np.testing.assert_almost_equal(x_back2, x_back0, 10)

        np.testing.assert_almost_equal(x_back0, x, 10)

    def test_lowlevel_fft_even(self):
        # Test lowlevel functions
        time = np.linspace(0,10,16) # NOTE: need a power of two for fft_py
        x = self.default_signal(time, mean=0)
        self.compare_with_npfft(time, x)

    def test_lowlevel_fft_odd(self):
        # Test lowlevel functions
        time = np.linspace(0,10,17) 
        x = self.default_signal(time, mean=0)
        self.compare_with_npfft(time, x)

    def test_fft_amplitude(self):
        dt=0.1
        t=np.arange(0,10,dt);
        f0=1;
        A=5;
        y=A*np.sin(2*np.pi*f0*t)
        f,Y,_=fft_amplitude(y,fs=1/dt,detrend=False)
        i=np.argmax(Y)
        self.assertAlmostEqual(Y[i],A)
        self.assertAlmostEqual(f[i],f0)

    def test_fft_binning(self):
        dt=0.1
        t=np.arange(0,10,dt);
        f0=1;
        A=5;
        y=A*np.sin(2*np.pi*f0*t)

        f,   Y, Info  = psd_binned(y, fs=1/dt, nPerDecade=10, detrend ='constant')
        f2, Y2, Info2 = psd       (y, fs=1/dt,           detrend ='constant')
        #print(f)
        #print(Y)

        #import matplotlib.pyplot as plt
        #fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
        #fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        #ax.plot( f2, Y2   , label='Full')
        #ax.plot( f,  Y    , label='Binned')
        #ax.set_xlabel('')
        #ax.set_ylabel('')
        #ax.legend()
        #plt.show()
    
if __name__ == '__main__':
    #TestSpectral().test_fft_binning()
    #TestSpectral().test_ifft()
    #TestSpectral().test_lowlevel_fft_even()
    #TestSpectral().test_lowlevel_fft_odd()
    unittest.main()

