import pandas as pd
import numpy as np
import re


def pd_interp1(x_new, xLabel, df):
    """ Interpolate a panda dataframe based on a set of new value
    This function assumes that the dataframe is a simple 2d-table
    """
    from pyFAST.tools.signal_analysis import multiInterp
    x_old = df[xLabel].values
    data_new=multiInterp(x_new, x_old, df.values.T)
    return pd.DataFrame(data=data_new.T, columns=df.columns.values)
    #nRow,nCol = df.shape
    #nRow = len(xnew)
    #data = np.zeros((nRow,nCol))
    #xref =df[xLabel].values.astype(float)
    #for col,i in zip(df.columns.values,range(nCol)):
    #    yref = df[col].values
    #    if yref.dtype!=float:
    #        raise Exception('Wrong type for yref, consider using astype(float)')
    #    data[:,i] = np.interp(xnew, xref, yref)
    #return pd.DataFrame(data=data, columns = df.columns)

def create_dummy_dataframe(size):
    return pd.DataFrame(data={'col1': np.linspace(0,1,size), 'col2': np.random.normal(0,1,size)})



def remap_df(df, ColMap, bColKeepNewOnly=False, inPlace=False, dataDict=None, verbose=False):
    """ 
    NOTE: see welib.fast.postpro

    Add/rename columns of a dataframe, potentially perform operations between columns

    dataDict: dictionary of data to be made available as "variable" in the column mapping
         'key' (new) : value (old)

    Example:

        ColumnMap={
          'WS_[m/s]'         : '{Wind1VelX_[m/s]}'             , # create a new column from existing one
          'RtTSR_[-]'        : '{RtTSR_[-]} * 2  +  {RtAeroCt_[-]}'    , # change value of column
          'RotSpeed_[rad/s]' : '{RotSpeed_[rpm]} * 2*np.pi/60 ', # new column [rpm] -> [rad/s]
          'q_p' :  ['Q_P_[rad]', '{PtfmSurge_[deg]}*np.pi/180']  # List of possible matches
        }
        # Read
        df = weio.read('FASTOutBin.outb').toDataFrame()
        # Change columns based on formulae, potentially adding new columns
        df = fastlib.remap_df(df, ColumnMap, inplace=True)

    """
    # Insert dataDict into namespace
    if dataDict is not None:
        for k,v in dataDict.items():
            exec('{:s} = dataDict["{:s}"]'.format(k,k))


    if not inPlace:
        df=df.copy()
    ColMapMiss=[]
    ColNew=[]
    RenameMap=dict()
    # Loop for expressions
    for k0,v in ColMap.items():
        k=k0.strip()
        if type(v) is not list:
            values = [v]
        else:
            values = v
        Found = False
        for v in values:
            v=v.strip()
            if Found:
                break # We avoid replacing twice
            if v.find('{')>=0:
                # --- This is an advanced substitution using formulae
                search_results = re.finditer(r'\{.*?\}', v)
                expr=v
                if verbose:
                    print('Attempt to insert column {:15s} with expr {}'.format(k,v))
                # For more advanced operations, we use an eval
                bFail=False
                for item in search_results:
                    col=item.group(0)[1:-1]
                    if col not in df.columns:
                        ColMapMiss.append(col)
                        bFail=True
                    expr=expr.replace(item.group(0),'df[\''+col+'\']')
                #print(k0, '=', expr)
                if not bFail:
                    df[k]=eval(expr)
                    ColNew.append(k)
                else:
                    print('[WARN] Column not present in dataframe, cannot evaluate: ',expr)
            else:
                #print(k0,'=',v)
                if v not in df.columns:
                    ColMapMiss.append(v)
                    if verbose:
                        print('[WARN] Column not present in dataframe: ',v)
                else:
                    if k in RenameMap.keys():
                        print('[WARN] Not renaming {} with {} as the key is already present'.format(k,v))
                    else:
                        RenameMap[k]=v
                        Found=True

    # Applying renaming only now so that expressions may be applied in any order
    for k,v in RenameMap.items():
        if verbose:
            print('Renaming column {:15s} > {}'.format(v,k))
        k=k.strip()
        iCol = list(df.columns).index(v)
        df.columns.values[iCol]=k
        ColNew.append(k)
    df.columns = df.columns.values # Hack to ensure columns are updated

    if len(ColMapMiss)>0:
        print('[FAIL] The following columns were not found in the dataframe:',ColMapMiss)
        #print('Available columns are:',df.columns.values)

    if bColKeepNewOnly:
        ColNew = [c for c,_ in ColMap.items() if c in ColNew]# Making sure we respec order from user
        ColKeepSafe = [c for c in ColNew if c in df.columns.values]
        ColKeepMiss = [c for c in ColNew if c not in df.columns.values]
        if len(ColKeepMiss)>0:
            print('[WARN] Signals missing and omitted for ColKeep:\n       '+'\n       '.join(ColKeepMiss))
        df=df[ColKeepSafe]
    return df

def changeUnits(df, flavor='SI', inPlace=True):
    """ Change units of a dataframe

    # TODO harmonize with dfToSIunits in welib.fast.tools.lin.py !
    """
    def splitunit(s):
        iu=s.rfind('[')
        if iu>0:
            return s[:iu], s[iu+1:].replace(']','')
        else:
            return s, ''
    def change_units_to_WE(s, c):
        """ 
        Change units to wind energy units
        s: channel name (string) containing units, typically 'speed_[rad/s]'
        c: channel (array)
        """
        svar, u = splitunit(s)
        u=u.lower()
        scalings = {}
        #        OLD      =     NEW
        scalings['rad/s'] =  (30/np.pi,'rpm') # TODO decide
        scalings['rad' ]  =   (180/np.pi,'deg')
        scalings['n']     =   (1e-3, 'kN')
        scalings['nm']    =   (1e-3, 'kNm')
        scalings['n-m']   =   (1e-3, 'kNm')
        scalings['n*m']   =   (1e-3, 'kNm')
        scalings['w']     =   (1e-3, 'kW')
        if u in scalings.keys():
            scale, new_unit = scalings[u]
            s = svar+'['+new_unit+']'
            c *= scale
        return s, c

    def change_units_to_SI(s, c):
        """ 
        Change units to SI units
        TODO, a lot more units conversion needed...will add them as we go
        s: channel name (string) containing units, typically 'speed_[rad/s]'
        c: channel (array)
        """
        svar, u = splitunit(s)
        u=u.lower()
        scalings = {}
        #        OLD      =     NEW
        scalings['rpm']   =  (np.pi/30,'rad/s') 
        scalings['rad' ]  =   (180/np.pi,'deg')
        scalings['deg/s' ] =   (np.pi/180,'rad/s')
        scalings['kn']     =   (1e3, 'N')
        scalings['knm']    =   (1e3, 'Nm')
        scalings['kn-m']   =   (1e3, 'Nm')
        scalings['kn*m']   =   (1e3, 'Nm')
        scalings['kw']     =   (1e3, 'W')
        if u in scalings.keys():
            scale, new_unit = scalings[u]
            s = svar+'['+new_unit+']'
            c *= scale
        return s, c

    if not inPlace:
        raise NotImplementedError()

    if flavor == 'WE':
        cols = []
        for i, colname in enumerate(df.columns):
            colname_new, df.iloc[:,i] = change_units_to_WE(colname, df.iloc[:,i])
            cols.append(colname_new)
        df.columns = cols
    elif flavor == 'SI':
        cols = []
        for i, colname in enumerate(df.columns):
            colname_new, df.iloc[:,i] = change_units_to_SI(colname, df.iloc[:,i])
            cols.append(colname_new)
        df.columns = cols
    else:
        raise NotImplementedError(flavor)
    return df

