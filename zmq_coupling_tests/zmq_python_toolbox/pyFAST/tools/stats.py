""" 
Set of tools for statistics 
  - measures (R^2, RMSE)
  - pdf distributions
  - Binning 

"""
import numpy as np
import pandas as pd

# --------------------------------------------------------------------------------}
# --- Stats measures 
# --------------------------------------------------------------------------------{
def comparison_stats(t1, y1, t2, y2, stats='sigRatio,eps,R2', method='mean', absVal=True):
    """
    y1: ref
    y2: other

    """
    from welib.tools.fatigue import equivalent_load

    sp=stats.split(',')
    stats = {}
    sStats=[]

    t1=np.asarray(t1).astype(float)
    y1=np.asarray(y1).astype(float)
    t2=np.asarray(t2).astype(float)
    y2=np.asarray(y2).astype(float)

    # Loop on statistics requested
    for s in sp:
        s= s.strip().lower()
        if s=='sigratio':
            # Ratio of standard deviation:
            sig_ref = float(np.nanstd(y1))
            sig_est = float(np.nanstd(y2))
            try:
                r_sig = sig_est/sig_ref
            except:
                r_sig = np.nan
            stats = {'sigRatio':r_sig}
            sStats+= [r'$\sigma_\mathrm{est}/\sigma_\mathrm{ref} = $'+r'{:.3f}'.format(r_sig)]

        elif s=='eps':
            # Mean relative error
            eps     = float(mean_rel_err(t1, y1, t2, y2, method=method, absVal=absVal))
            stats['eps'] = eps
            sStats+=['$\epsilon=$'+r'{:.1f}%'.format(eps)]

        elif s=='r2':
            # Rsquare
            R2 = float(rsquare(y2, y1)[0])
            stats['R2'] = R2
            sStats+=[r'$R^2=$'+r'{:.3f}'.format(R2)]

        elif s=='epsleq':
            Leq1 = equivalent_load(t1, y1, m=5, bins=100, method='fatpack')
            Leq2 = equivalent_load(t2, y2, m=5, bins=100, method='fatpack')
            epsLeq = (Leq2-Leq1)/Leq1*100
            stats['epsLeq'] = epsLeq
            sStats+=[r'$\epsilon L_{eq}=$'+r'{:.1f}%'.format(epsLeq)]

        else:
            raise NotImplementedError(s)
    sStats=' - '.join(sStats)
    return stats, sStats



def rsquare(y, f, c = True): 
    """ Compute coefficient of determination of data fit model and RMSE
    [r2 rmse] = rsquare(y,f)
    [r2 rmse] = rsquare(y,f,c)
    RSQUARE computes the coefficient of determination (R-square) value from
    actual data Y and model data F. The code uses a general version of
    R-square, based on comparing the variability of the estimation errors
    with the variability of the original values. RSQUARE also outputs the
    root mean squared error (RMSE) for the user's convenience.
    Note: RSQUARE ignores comparisons involving NaN values.
    INPUTS
      Y       : Actual data
      F       : Model fit
    
    # OPTION
      C       : Constant term in model
                R-square may be a questionable measure of fit when no
              constant term is included in the model.
      [DEFAULT] TRUE : Use traditional R-square computation
               FALSE : Uses alternate R-square computation for model
                     without constant term [R2 = 1 - NORM(Y-F)/NORM(Y)]
    # OUTPUT
      R2      : Coefficient of determination
      RMSE    : Root mean squared error """
    # Sanity
    if not np.all(y.shape == f.shape) :
        raise Exception('Y and F must be the same size')
    y = np.asarray(y).astype(float)
    f = np.asarray(f).astype(float)
    # Check for NaN
    tmp = np.logical_not(np.logical_or(np.isnan(y),np.isnan(f))) 
    y = y[tmp]
    f = f[tmp]
    if c:
        r2 = max(0,1-np.sum((y-f)**2)/np.sum((y-np.mean(y))** 2))
    else:
        r2 = 1 - np.sum((y - f) ** 2) / np.sum((y) ** 2)
        if r2 < 0:
            import warnings
            warnings.warn('Consider adding a constant term to your model')
            r2 = 0
    rmse = np.sqrt(np.mean((y - f) ** 2))
    return r2,rmse

def mean_rel_err(t1=None, y1=None, t2=None, y2=None, method='meanabs', verbose=False, varname='', absVal=True):
    """ 
    return mean relative error in % 

    Methods: 
      'mean'   : 100 * |y1-y2|/mean(y1)
      'meanabs': 100 * |y1-y2|/mean(|y1|)
      'minmax': y1 and y2 scaled between 0.5 and 1.5
                |y1s-y2s|/|y1|
      '0-2': signals are scalled between 0 & 2
    """
    def myabs(y):
        if absVal:
            return np.abs(y)
        else:
            return y


    if t1 is None and t2 is None:
        pass
    else:
        if len(y1)!=len(y2):
            y2=np.interp(t1,t2,y2)
    if method=='mean':
        # Method 1 relative to mean
        ref_val = np.nanmean(y1)
        meanrelerr = np.nanmean(myabs(y2-y1)/ref_val)*100 
    elif method=='meanabs':
        ref_val = np.nanmean(abs(y1))
        meanrelerr = np.nanmean(myabs(y2-y1)/ref_val)*100 
    elif method=='loc':
        meanrelerr = np.nanmean(myabs(y2-y1)/abs(y1))*100 
    elif method=='minmax':
        # Method 2 scaling signals
        Min=min(np.nanmin(y1), np.nanmin(y2))
        Max=max(np.nanmax(y1), np.nanmax(y2))
        y1=(y1-Min)/(Max-Min)+0.5
        y2=(y2-Min)/(Max-Min)+0.5
        meanrelerr = np.nanmean(myabs(y2-y1)/np.abs(y1))*100 
    elif method=='1-2':
        # transform values from 1 to 2
        Min=min(np.nanmin(y1), np.nanmin(y2))
        Max=max(np.nanmax(y1), np.nanmax(y2))
        y1 = (y1-Min)/(Max-Min)+1
        y2 = (y2-Min)/(Max-Min)+1
        meanrelerr = np.nanmean(myabs(y2-y1)/np.abs(y1))*100
    else:
        raise Exception('Unknown method',method)

    if verbose:
        if len(varname)>0:
            print('Mean rel error {:15s} {:7.2f} %'.format(varname, meanrelerr))
        else:
            print('Mean rel error {:7.2f} %'.format( meanrelerr))
    return meanrelerr


# --------------------------------------------------------------------------------}
# --- PDF 
# --------------------------------------------------------------------------------{
def pdf(y, method='histogram', n=50, **kwargs):
    """ 
    Compute the probability density function.
    Wrapper over the different methods present in this package
    """
    if method =='sns':
        xh, yh = pdf_sns(y, nBins=n, **kwargs)
    elif method =='gaussian_kde':
        xh, yh = pdf_gaussian_kde(y, nOut=n, **kwargs)
    elif method =='histogram':
        xh, yh = pdf_histogram(y, nBins=n, **kwargs)
    else:
        raise NotImplementedError(f'pdf method: {method}')
    return xh, yh


def pdf_histogram(y,nBins=50, norm=True, count=False):
    yh, xh = np.histogram(y[~np.isnan(y)], bins=nBins)
    dx   = xh[1] - xh[0]
    xh  = xh[:-1] + dx/2
    if count:
        yh  = yh / (len(n)*dx) # TODO DEBUG /VERIFY THIS
    else:
        yh  = yh / (nBins*dx) 
    if norm:
        yh=yh/np.trapz(yh,xh)
    return xh,yh

def pdf_gaussian_kde(data, bw='scott', nOut=100, cut=3, clip=(-np.inf,np.inf)):
    """ 
    Returns a smooth probability density function (univariate kernel density estimate - kde) 
    Inspired from `_univariate_kdeplot` from `seaborn.distributions`

    INPUTS:
        bw:  float defining bandwidth or method (string) to find it (more or less sigma)   
        cut: number of bandwidth kept for x axis (e.g. 3 sigmas)
        clip: (xmin, xmax) values
    OUTPUTS:
        x, y: where y(x) = pdf(data)
    """
    from scipy import stats
    from six import string_types

    data = np.asarray(data)
    data = data[~np.isnan(data)]
    # Gaussian kde
    kde  = stats.gaussian_kde(data, bw_method = bw)
    # Finding a relevant support (i.e. x values)
    if isinstance(bw, string_types):
        bw_ = "scotts" if bw == "scott" else bw
        bw = getattr(kde, "%s_factor" % bw_)() * np.std(data)
    x_min = max(data.min() - bw * cut, clip[0])
    x_max = min(data.max() + bw * cut, clip[1])
    x = np.linspace(x_min, x_max, nOut)
    # Computing kde on support
    y = kde(x)
    return x, y


def pdf_sklearn(y):
    #from sklearn.neighbors import KernelDensity
    #kde = KernelDensity(kernel='gaussian', bandwidth=0.75).fit(y) #you can supply a bandwidth
    #x=np.linspace(0,5,100)[:, np.newaxis]
    #log_density_values=kde.score_samples(x)
    #density=np.exp(log_density)
    pass

def pdf_sns(y,nBins=50):
    import seaborn.apionly as sns
    hh=sns.distplot(y,hist=True,norm_hist=False).get_lines()[0].get_data()
    xh=hh[0]
    yh=hh[1]
    return xh,yh

# --------------------------------------------------------------------------------}
# --- Binning 
# --------------------------------------------------------------------------------{
def bin_DF(df, xbins, colBin, stats='mean'):
    """ 
    Perform bin averaging of a dataframe
    INPUTS:
      - df   : pandas dataframe
      - xBins: end points delimiting the bins, array of ascending x values
      - colBin: column name (string) of the dataframe, used for binning 
    OUTPUTS:
       binned dataframe, with additional columns 'Counts' for the number 

    """
    if colBin not in df.columns.values:
        raise Exception('The column `{}` does not appear to be in the dataframe'.format(colBin))
    xmid      = (xbins[:-1]+xbins[1:])/2
    df['Bin'] = pd.cut(df[colBin], bins=xbins, labels=xmid ) # Adding a column that has bin attribute
    if stats=='mean':
        df2       = df.groupby('Bin', observed=False).mean()                     # Average by bin
    elif stats=='std':
        df2       = df.groupby('Bin', observed=False).std()                     # std by bin
    # also counting
    df['Counts'] = 1
    dfCount=df[['Counts','Bin']].groupby('Bin', observed=False).sum()
    df2['Counts'] = dfCount['Counts']
    # Just in case some bins are missing (will be nan)
    df2       = df2.reindex(xmid)
    return df2

def bin_signal(x, y, xbins=None, stats='mean', nBins=None):
    """ 
    Perform bin averaging of a signal
    INPUTS:
      - x: x-values 
      - y: y-values, signal values
      - xBins: end points delimiting the bins, array of ascending x values
    OUTPUTS:
      - xBinned, yBinned

    """
    if xbins is None:
        xmin, xmax = np.min(x), np.max(x)
        dx = (xmax-xmin)/nBins
        xbins=np.arange(xmin, xmax+dx/2, dx)
    df = pd.DataFrame(data=np.column_stack((x,y)), columns=['x','y'])
    df2 = bin_DF(df, xbins, colBin='x', stats=stats)
    return df2['x'].values, df2['y'].values



def bin2d_signal(x, y, z, xbins=None, ybins=None, nXBins=None, nYBins=None):
    """ 
    Bin signal z based on x and y values using xbins and ybins

    """
    if xbins is None:
        xmin, xmax = np.min(x), np.max(x)
        dx = (xmax-xmin)/nXBins
        xbins=np.arange(xmin, xmax+dx/2, dx)
    if ybins is None:
        ymin, ymax = np.min(y), np.max(y)
        dy = (ymax-ymin)/nYBins
        ybins=np.arange(ymin, ymax+dy/2, dy)

    x = np.asarray(x).flatten()
    y = np.asarray(y).flatten()
    z = np.asarray(z).flatten()

    Counts = np.zeros((len(xbins)-1, len(ybins)-1))
    XMean  = np.zeros((len(xbins)-1, len(ybins)-1))*np.nan
    YMean  = np.zeros((len(xbins)-1, len(ybins)-1))*np.nan
    ZMean  = np.zeros((len(xbins)-1, len(ybins)-1))*np.nan
    ZStd   = np.zeros((len(xbins)-1, len(ybins)-1))*np.nan

    xmid = xbins[:-1] + np.diff(xbins)/2
    ymid = ybins[:-1] + np.diff(ybins)/2
    YMid, XMid = np.meshgrid(ymid, xmid)

    for ixb, xb in enumerate(xbins[:-1]):
        print(ixb)
        bX = np.logical_and(x >= xb, x <= xbins[ixb+1]) # TODO decide on bounds
        for iyb, yb in enumerate(ybins[:-1]):
            bY = np.logical_and(y >= yb, y <= ybins[iyb+1]) # TODO decide on bounds

            bXY = np.logical_and(bX, bY)
            Counts[ixb, iyb] = sum(bXY)
            if Counts[ixb,iyb]>0:
                ZMean [ixb, iyb] = np.mean(z[bXY])
                ZStd  [ixb, iyb] = np.std( z[bXY])
                XMean [ixb, iyb] = np.mean(x[bXY])
                YMean [ixb, iyb] = np.mean(y[bXY])

    return XMean, YMean, ZMean, ZStd, Counts, XMid, YMid






def azimuthal_average_DF(df, psiBin=np.arange(0,360+1,10), colPsi='Azimuth_[deg]', tStart=None, colTime='Time_[s]'):
    """ 
    Average a dataframe based on azimuthal value
    Returns a dataframe with same amount of columns as input, and azimuthal values as index
    """
    if tStart is not None:
        if colTime not in df.columns.values:
            raise Exception('The column `{}` does not appear to be in the dataframe'.format(colTime))
        df=df[ df[colTime]>tStart].copy()

    dfPsi= bin_DF(df, psiBin, colPsi, stats='mean')
    if np.any(dfPsi['Counts']<1):
        print('[WARN] some bins have no data! Increase the bin size.')

    return dfPsi


def azimuthal_std_DF(df, psiBin=np.arange(0,360+1,10), colPsi='Azimuth_[deg]', tStart=None, colTime='Time_[s]'):
    """ 
    Average a dataframe based on azimuthal value
    Returns a dataframe with same amount of columns as input, and azimuthal values as index
    """
    if tStart is not None:
        if colTime not in df.columns.values:
            raise Exception('The column `{}` does not appear to be in the dataframe'.format(colTime))
        df=df[ df[colTime]>tStart].copy()

    dfPsi= bin_DF(df, psiBin, colPsi, stats='std')
    if np.any(dfPsi['Counts']<1):
        print('[WARN] some bins have no data! Increase the bin size.')

    return dfPsi




