"""
Context:

Logarithmic decrement:

    delta = 1/N log [ x(t) / x(t+N T_d)]  = 2 pi zeta / sqrt(1-zeta^2)

Damping ratio:

    zeta = delta / sqrt( 4 pi^2 + delta^2 )

Damped period, frequency:

    Td = 2pi / omega_d

    omegad = omega_0 sqrt(1-zeta**2)

Damping exponent:
     

    alpha = zeta omega_0 = delta/ T_d


"""

import numpy as np

__all__  = ['freqDampEstimator']
__all__ += ['freqDampFromPeaks']
__all__ += ['zetaEnvelop']
__all__ += ['TestDamping']

def indexes(y, thres=0.3, min_dist=1, thres_abs=False):
    """Peak detection routine.

    Finds the numeric index of the peaks in *y* by taking its first order difference. By using
    *thres* and *min_dist* parameters, it is possible to reduce the number of
    detected peaks. *y* must be signed.

    Parameters
    ----------
    y : ndarray (signed)
        1D amplitude data to search for peaks.
    thres : float, defining threshold. Only the peaks with amplitude higher than the
        threshold will be detected.
        if thres_abs is False: between [0., 1.], normalized threshold. 
    min_dist : int
        Minimum distance between each detected peak. The peak with the highest
        amplitude is preferred to satisfy this constraint.
    thres_abs: boolean
        If True, the thres value will be interpreted as an absolute value, instead of
        a normalized threshold.

    Returns
    -------
    ndarray
        Array containing the numeric indexes of the peaks that were detected
    """
    if isinstance(y, np.ndarray) and np.issubdtype(y.dtype, np.unsignedinteger):
        raise ValueError("y must be signed")

    if not thres_abs:
        thres = thres * (np.max(y) - np.min(y)) + np.min(y)
        
    min_dist = int(min_dist)

    # compute first order difference
    dy = np.diff(y)

    # propagate left and right values successively to fill all plateau pixels (0-value)
    zeros,=np.where(dy == 0)
    
    # check if the signal is totally flat
    if len(zeros) == len(y) - 1:
        return np.array([])

    if len(zeros):
        # compute first order difference of zero indexes
        zeros_diff = np.diff(zeros)
        # check when zeros are not chained together
        zeros_diff_not_one, = np.add(np.where(zeros_diff != 1), 1)
        # make an array of the chained zero indexes
        zero_plateaus = np.split(zeros, zeros_diff_not_one)

        # fix if leftmost value in dy is zero
        if zero_plateaus[0][0] == 0:
            dy[zero_plateaus[0]] = dy[zero_plateaus[0][-1] + 1]
            zero_plateaus.pop(0)

        # fix if rightmost value of dy is zero
        if len(zero_plateaus) and zero_plateaus[-1][-1] == len(dy) - 1:
            dy[zero_plateaus[-1]] = dy[zero_plateaus[-1][0] - 1]
            zero_plateaus.pop(-1)

        # for each chain of zero indexes
        for plateau in zero_plateaus:
            median = np.median(plateau)
            # set leftmost values to leftmost non zero values
            dy[plateau[plateau < median]] = dy[plateau[0] - 1]
            # set rightmost and middle values to rightmost non zero values
            dy[plateau[plateau >= median]] = dy[plateau[-1] + 1]

    # find the peaks by using the first order difference
    peaks = np.where((np.hstack([dy, 0.]) < 0.)
                     & (np.hstack([0., dy]) > 0.)
                     & (np.greater(y, thres)))[0]

    # handle multiple peaks, respecting the minimum distance
    if peaks.size > 1 and min_dist > 1:
        highest = peaks[np.argsort(y[peaks])][::-1]
        rem = np.ones(y.size, dtype=bool)
        rem[peaks] = False

        for peak in highest:
            if not rem[peak]:
                sl = slice(max(0, peak - min_dist), peak + min_dist + 1)
                rem[sl] = True
                rem[peak] = False

        peaks = np.arange(y.size)[~rem]

    return peaks


#indexes =indexes(x, thres=0.02/max(x), min_dist=1, thres_abs=true)
def logDecFromThreshold(x, threshold=None, bothSides=False, decay=True):
    """ Detect maxima in a signal, computes the log deg based on it 
    """
    if bothSides:
        ldPos,iTPos,stdPos,IPos,vldPos = logDecFromThreshold( x, threshold=threshold, decay=decay)
        ldNeg,iTNeg,stdNeg,INeg,vldNeg = logDecFromThreshold(-x, threshold=threshold, decay=decay)
        return (ldPos+ldNeg)/2, (iTPos+iTNeg)/2, (stdPos+stdNeg)/2, (IPos,INeg), (vldPos, vldNeg)

    if threshold is None:
        threshold = np.mean(abs(x-np.mean(x)))/3;
    I =indexes(x, thres=threshold, min_dist=1, thres_abs=True)
    # Estimating "index" period
    iT = round(np.median(np.diff(I)));
    vn=np.arange(0,len(I)-1)+1
    # Quick And Dirty Way using one as ref and assuming all periods were found
    if decay:
        # For a decay we take the first peak as a reference
        vLogDec  = 1/vn*np.log( x[I[0]]/x[I[1:]] ) # Logarithmic decrement
    else:
        # For negative damping we take the last peak as a reference
        vLogDec  = 1/vn*np.log( x[I[-2::-1]]/x[I[-1]] ) # Logarithmic decrement
    logdec     = np.mean(vLogDec);
    std_logdec = np.std(vLogDec) ;
    return logdec, iT, std_logdec, I, vLogDec

def logDecTwoTimes(x, t, i1, i2, Td):
    t1, t2 = t[i1], t[i2]
    x1, x2 = x[i1], x[i2]
    N = (t2-t1)/Td
    logdec = 1/N * np.log(x1/x2)
    return logdec

def zetaTwoTimes(x, t, i1, i2, Td):
    logdec = logDecTwoTimes(x, t, i1, i2, Td)
    zeta = logdec/np.sqrt(4*np.pi**2 + logdec**2) # damping ratio
    return zeta

def zetaRange(x, t, IPos, INeg, Td, decay):
    """ 
    Compute naive zeta based on different peak values (first, mid, last)
    """
    def naivezeta(i1, i2):
        zetaP = zetaTwoTimes(  x, t, IPos[i1], IPos[i2], Td)
        zetaN = zetaTwoTimes( -x, t, INeg[i1], INeg[i2], Td)
        return [zetaP, zetaN]
    zetas = []
    # --- Computing naive log dec from first and last peaks
    zetas += naivezeta(0, -1)
    # --- Computing naive log dec from one peak and the middle one
    if len(IPos)>3 and len(INeg)>3:
        if decay:
            i1, i2 = 0, int(len(IPos)/2)
        else:
            i1, i2 = -int(len(IPos)/2), -1
        zetas += naivezeta(i1, i2)
    zetaSup  = np.max(zetas)
    zetaInf  = np.min(zetas)
    zetaMean = np.mean(zetas)
    return zetaSup, zetaInf, zetaMean

def zetaEnvelop(x, t, omega0, zeta, iRef, moff=0):
    """ NOTE: x is assumed to be centered on 0"""
    m = np.mean(x)
    tref = t[iRef]
    Aref = x[iRef]-m
    epos =  Aref*np.exp(-zeta*omega0*(t-tref))+m+moff
    eneg = -Aref*np.exp(-zeta*omega0*(t-tref))+m+moff
    return epos, eneg


def freqDampFromPeaks(x, t, threshold=None, plot=False, refPoint='mid'):
    """ 
    Use Upper and lower peaks to compute log decrements between neighboring peaks
    Previously called logDecFromDecay.
    """
    info = {}
    x = np.array(x).copy()
    m = np.mean(x)
    x = x-m # we remove the mean once and for all
    if threshold is None:
        threshold = np.mean(abs(x))/3
    
    dt = t[1]-t[0] # todo signal with dt not uniform

    # Is it a decay or an exloding signal
    xstart, xend = np.array_split(np.abs(x),2)
    decay= np.mean(xstart)> np.mean(xend)

    # --- Computing log decs from positive and negative side and taking the mean
    logdec,iT,std,(IPos,INeg), (vLogDecPos, vLogDecNeg) = logDecFromThreshold( x, threshold=threshold, bothSides=True, decay=decay)

    # --- Finding damped period
    Td = iT*dt # Period of damped oscillations. Badly estimated due to dt resolution
    # % Better estimate of period
    # [T,~,iCross]=fGetPeriodFromZeroCrossing(x(1:IPos(end)),dt);
    fd = 1/Td           

    # --- Naive ranges
    zetaMax, zetaMin, zetaMean = zetaRange(x, t, IPos, INeg, Td, decay)

    zeta  = logdec/np.sqrt(4*np.pi**2 + logdec**2 ) # Damping Ratio
    fn = fd/np.sqrt(1-zeta**2)
    T0 = 1/fn
    omega0=2*np.pi*fn
    # --- Model
    # Estimating signal params 
    alpha   = logdec/Td    
    omega   = 2*np.pi*fd  
    # Find amplitude at a reference peak
    # (we chose the middle peak of the time series to minimize period drift before and after)
    # We will adjust for phase and time offset later
    i1 = IPos[int(len(IPos)/2)]
    A1 = x[i1]
    t1 = dt*i1
    # --- Find a zero up-crossing around our value of reference for phase determination
    XX=x[i1:]
    ineg = i1+np.where(XX<0)[0][0]
    ipos = ineg-1
    xcross = [x[ipos],x[ineg]]
    icross = [ipos,ineg]
    i0 = np.interp(0,xcross,icross) # precise 0-up-crossing
    t0   = dt*i0
    phi0 = np.mod(2*np.pi- omega*t0+np.pi/2,2*np.pi);
    # --- Model
    A      = A1/(np.exp(-alpha*t1)*np.cos(omega*t1+phi0)); # Adjust for phase and time offset
    x_model = A*np.exp(-alpha*t)*np.cos(omega*t+phi0)+m;
    epos   =  A*np.exp(-alpha*t)+m
    eneg   = -A*np.exp(-alpha*t)+m

    if plot:
        if refPoint=='mid':
            iRef = i1
        elif refPoint=='start':
            iRef = IPos[0]
        else:
            iRef = IPos[-1]
        import matplotlib.pyplot as plt
        print('LogDec.: {:.4f} - Damping ratio: {:.4f} - F_n: {:.4f} - F_d: {:.4f} - T_d:{:.3f} - T_n:{:.3f}'.format(logdec, zeta, fn, fd, Td,T0))
        fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
        fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        ax.plot(t, x+m)
        ax.plot(t[IPos],x[IPos]+m,'o')
        ax.plot(t[INeg],x[INeg]+m,'o')
        epos, eneg = zetaEnvelop(x, t, omega0, zeta, iRef=iRef, moff=m) 
        ax.plot(t ,epos, 'k--', label=r'$\zeta={:.4f}$'.format(zeta))
        ax.plot(t ,eneg, 'k--')
        epos, eneg = zetaEnvelop(x, t, omega0, zetaMax, iRef=iRef, moff=m) 
        ax.plot(t ,epos, 'b:', label=r'$\zeta={:.4f}$'.format(zetaMax))
        ax.plot(t ,eneg, 'b:')
        epos, eneg = zetaEnvelop(x, t, omega0, zetaMin, iRef=iRef, moff=m) 
        ax.plot(t ,epos, 'r:', label=r'$\zeta={:.4f}$'.format(zetaMin))
        ax.plot(t ,eneg, 'r:')
        #ax.plot(t ,x_model,'k:')
        #ax.legend()
        dx = np.max(abs(x-m))
        ax.set_ylim([m-dx*1.1 , m+dx*1.1])

    # We return a dictionary
    info['zeta'] = zeta
    info['fd']   = fd
    info['Td']   = Td
    info['fn']   = fn
    info['omega0'] = omega0
    info['IPos'] = IPos
    info['INeg'] = INeg
    # TODO
    info['x_model'] = x_model
    info['epos']    = epos
    info['eneg']    = eneg
    # 
    info['zeta']     = zeta
    info['zetaMin']  = zetaMin
    info['zetaMax']  = zetaMax
    info['zetaMean'] = zetaMean

    return fn, zeta, info


def freqDampEstimator(x, t, opts):
    """ 
    Estimate natural frequency and damping ratio.
    Wrapper function to use different methods.

    """
    if opts['method']=='fromPeaks':
        fn, zeta, info = freqDampFromPeaks(x, t)
    else:
        raise NotImplementedError()
    return fn, zeta, info



# --------------------------------------------------------------------------------}
# --- Unittests
# --------------------------------------------------------------------------------{
import unittest

class TestDamping(unittest.TestCase):

    def test_logdec_from_peaks(self):
        plot = (__name__ == '__main__')

        for zeta in [0.1, -0.01]:
            T0    = 10                              
            Td    = T0 / np.sqrt(1-zeta**2)
            delta = 2*np.pi*zeta/np.sqrt(1-zeta**2) # logdec
            alpha = delta/Td
            t     = np.linspace(0,30*Td,2000)
            x     = np.cos(2*np.pi/Td*t)*np.exp(-alpha*t)+10;
            fn, zeta_out, info = freqDampFromPeaks(x, t, plot=plot)
            self.assertAlmostEqual(zeta  , zeta_out,4)
            self.assertAlmostEqual(1/T0  , fn      ,2)

        if __name__ == '__main__':
            import matplotlib.pyplot as plt
            plt.show()
    
if __name__ == '__main__':
    unittest.main()
#     import matplotlib.pyplot as plt
#     import pydatview.io as weio
#     df= weio.read('DampingExplodingExample2.csv').toDataFrame()
#     M = df.values
#     x= M[:,1]
#     t= M[:,0]
#     #for zeta in [-0.01, 0.1]:
#     #    T0    = 30                              
#     #    Td    = T0 / np.sqrt(1-zeta**2)
#     #    delta = 2*np.pi*zeta/np.sqrt(1-zeta**2) # logdec
#     #    alpha = delta/Td
#     #    x     = np.cos(2*np.pi/Td*t)*np.exp(-alpha*t)+10;
#     #    df.insert(1,'PureDecay{}'.format(zeta), x)
# 
#     #df.to_csv('DECAY_Example.csv',index=False, sep=',')
# 
# #     Td    = 10                                    
# #     zeta  = -0.01                            # damping ratio (<1)
# #     delta = 2*np.pi*zeta/np.sqrt(1-zeta**2)  # logdec
# #     alpha = delta/Td
# #     t     = np.linspace(0,30*Td,1000)
# #     x     = np.cos(2*np.pi/Td*t)*np.exp(-alpha*t)+10;
# # 
# #     fn, zeta, info = freqDampFromPeaks(x, t, plot=True, refPoint='mid')
#     plt.show()




