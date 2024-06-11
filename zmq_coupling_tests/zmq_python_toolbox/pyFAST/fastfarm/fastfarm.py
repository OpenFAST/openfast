import os
import glob
import numpy as np
import pandas as pd
from pyFAST.input_output.fast_input_file import FASTInputFile
from pyFAST.input_output.fast_output_file import FASTOutputFile
from pyFAST.input_output.turbsim_file import TurbSimFile
import pyFAST.postpro as fastlib

# --------------------------------------------------------------------------------}
# --- Small helper functions
# --------------------------------------------------------------------------------{
def insertTN(s,i,nWT=1000, noLeadingZero=False):
    """ insert turbine number in name """
    if nWT<10:
        fmt='{:d}'
    elif nWT<100:
        fmt='{:02d}'
    else:
        fmt='{:03d}'

    if noLeadingZero:
        fmt='{:d}'

    if s.find('T1')>=0:
        s=s.replace('T1','T'+fmt.format(i))
    elif s.find('T0')>=0:
        print('this should not be printed')
        s=s.replace('T0','T'+fmt.format(i))
    else:
        sp=os.path.splitext(s)
        s=sp[0]+'_T'+fmt.format(i)+sp[1]
    return s
def forceCopyFile (sfile, dfile):
    # ---- Handling error due to wrong mod
    if os.path.isfile(dfile):
        if not os.access(dfile, os.W_OK):
            os.chmod(dfile, stat.S_IWUSR)
    #print(sfile, ' > ', dfile)
    shutil.copy2(sfile, dfile)

# --------------------------------------------------------------------------------}
# --- Tools to create fast farm simulations
# --------------------------------------------------------------------------------{
def writeFSTandDLL(FstT1Name, nWT):
    """ 
    Write FST files for each turbine, with different ServoDyn files and DLL 
    FST files, ServoFiles, and DLL files will be written next to their turbine 1
    files, with name Ti. 

    FstT1Name: absolute or relative path to the Turbine FST file
    """ 

    FstT1Full = os.path.abspath(FstT1Name).replace('\\','/')
    FstDir  = os.path.dirname(FstT1Full)

    fst=FASTInputFile(FstT1Name)
    SrvT1Name    = fst['ServoFile'].strip('"')
    SrvT1Full    = os.path.join(FstDir, SrvT1Name).replace('\\','/')
    SrvDir       = os.path.dirname(SrvT1Full)
    SrvT1RelFst  = os.path.relpath(SrvT1Full,FstDir)
    if os.path.exists(SrvT1Full):
        srv=FASTInputFile(SrvT1Full)
        DLLT1Name = srv['DLL_FileName'].strip('"')
        DLLT1Full = os.path.join(SrvDir, DLLT1Name)
        if os.path.exists(DLLT1Full):
            servo=True
        else:
            print('[Info] DLL file not found, not copying servo and dll files ({})'.format(DLLT1Full))
            servo=False
    else:
        print('[Info] ServoDyn file not found, not copying servo and dll files ({})'.format(SrvT1Full))
        servo=False

    #print(FstDir)
    #print(FstT1Full)
    #print(SrvT1Name)
    #print(SrvT1Full)
    #print(SrvT1RelFst)

    for i in np.arange(2,nWT+1):
        FstName     = insertTN(FstT1Name,i,nWT)
        if servo:
            # TODO handle the case where T1 not present
            SrvName     = insertTN(SrvT1Name,i,nWT)
            DLLName     = insertTN(DLLT1Name,i,nWT)
            DLLFullName = os.path.join(SrvDir, DLLName)

        print('')
        print('FstName: ',FstName)
        if servo:
            print('SrvName: ',SrvName)
            print('DLLName: ',DLLName)
            print('DLLFull: ',DLLFullName)

        # Changing main file
        if servo:
            fst['ServoFile']='"'+SrvName+'"'
        fst.write(FstName)
        if servo:
            # Changing servo file
            srv['DLL_FileName']='"'+DLLName+'"'
            srv.write(SrvName)
            # Copying dll
            forceCopyFile(DLLT1Full, DLLFullName)



def rectangularLayoutSubDomains(D,Lx,Ly):
    """ Retuns position of turbines in a rectangular layout 
    TODO, unfinished function parameters
    """
    # --- Parameters
    D          = 112  # turbine diameter [m]
    Lx         = 3840 # x dimension of precusor
    Ly         = 3840 # y dimension of precusor
    Height     = 0    # Height above ground, likely 0 [m]
    nDomains_x = 2    # number of domains in x
    nDomains_y = 2    # number of domains in y
    # --- 36 WT
    nx         = 3    # number of turbines to be placed along x in one precursor domain
    ny         = 3    # number of turbines to be placed along y in one precursor domain
    StartX     = 1/2  # How close do we start from the x boundary
    StartY     = 1/2  # How close do we start from the y boundary
    # --- Derived parameters
    Lx_Domain = Lx * nDomains_x   # Full domain size
    Ly_Domain = Ly * nDomains_y
    DeltaX = Lx / (nx)          # Turbine spacing
    DeltaY = Ly / (ny)
    xWT = np.arange(DeltaX*StartX,Lx_Domain,DeltaX) # Turbine positions 
    yWT = np.arange(DeltaY*StartY,Ly_Domain,DeltaY)

    print('Full domain size [D]  :  {:.2f} x {:.2f}  '.format(Lx_Domain/D, Ly_Domain/D))
    print('Turbine spacing  [D]  : {:.2f} x  {:.2f} '.format(DeltaX/D,DeltaX/D))
    print('Number of turbines    : {:d} x {:d} = {:d}'.format(len(xWT),len(yWT),len(xWT)*len(yWT)))

    XWT,YWT=np.meshgrid(xWT,yWT)
    ZWT=XWT*0+Height

    # --- Export coordinates only
    M=np.column_stack((XWT.ravel(),YWT.ravel(),ZWT.ravel()))
    np.savetxt('Farm_Coordinates.csv', M, delimiter=',',header='X_[m], Y_[m], Z_[m]')
    print(M)

    return XWT, YWT, ZWT


def fastFarmTurbSimExtent(TurbSimFilename, hubHeight, D, xWT, yWT, Cmeander=1.9, chord_max=3, extent_X=1.1, extent_YZ=1.1, meanUAtHubHeight=False):
    """ 
    Determines "Ambient Wind" box parametesr for FastFarm, based on a TurbSimFile ('bts')

    Implements the guidelines listed here:
        https://openfast.readthedocs.io/en/dev/source/user/fast.farm/ModelGuidance.html

    INPUTS:
     - TurbSimFilename: name of the BTS file used in the FAST.Farm simulation
     - hubHeight      : Hub height [m]
     - D              : turbine diameter [m]
     - xWT            : vector of x positions of the wind turbines (e.g. [0,300,600])
     - yWT            : vector of y positions of the wind turbines (e.g. [0,0,0])
     - Cmeander       : parameter for meandering used in FAST.Farm [-]
     - chord_max      : maximum chord of the wind turbine blade. Used to determine the high resolution 
     - extent_X       : x-extent of high res box in diamter around turbine location
     - extent_YZ      : y-extent of high res box in diamter around turbine location

    """
    # --- TurbSim data
    ts = TurbSimFile(TurbSimFilename)

    if meanUAtHubHeight:
        # Use Hub Height to determine convection velocity
        iy,iz   = ts.closestPoint(y=0,z=hubHeight)
        meanU   = ts['u'][0,:,iy,iz].mean()
    else:
        # Use middle of the box to determine convection velocity
        zMid, meanU =  ts.midValues()

    return fastFarmBoxExtent(ts.y, ts.z, ts.t, meanU, hubHeight, D, xWT, yWT, Cmeander=Cmeander, chord_max=chord_max, extent_X=extent_X, extent_YZ=extent_YZ) 

def fastFarmBoxExtent(yBox, zBox, tBox, meanU, hubHeight, D, xWT, yWT, 
        Cmeander=1.9, chord_max=3, extent_X=1.1, extent_YZ=1.1, 
        extent_wake=8, LES=False):
    """ 
    Determines "Ambient Wind" box parametesr for FastFarm, based on turbulence box parameters
    INPUTS:
     - yBox  : y vector of grid points of the box
     - zBox  : z vector of grid points of the box
     - tBox  : time vector of the box
     - meanU : mean velocity used to convect the box
     - hubHeight      : Hub height [m]
     - D              : turbine diameter [m]
     - xWT            : vector of x positions of the wind turbines (e.g. [0,300,600])
     - yWT            : vector of y positions of the wind turbines (e.g. [0,0,0])
     - Cmeander       : parameter for meandering used in FAST.Farm [-]
     - chord_max      : maximum chord of the wind turbine blade. Used to determine the high resolution 
     - extent_X       : x-extent of high-res box (in diameter) around turbine location
     - extent_YZ      : y-extent of high-res box (in diameter) around turbine location
     - extent_wake    : extent of low-res box (in diameter) to add beyond the "last" wind turbine
     - LES: False for TurbSim box, true for LES. Perform additional checks for LES.
    """
    if LES:
        raise NotImplementedError()
    # --- Box resolution and extents
    dY_Box      = yBox[1]-yBox[0]
    dZ_Box      = zBox[1]-zBox[0]
    dT_Box      = tBox[1]-tBox[0]
    dX_Box      = dT_Box * meanU
    Z0_Box      = zBox[0]
    LY_Box      = yBox[-1]-yBox[0]
    LZ_Box      = zBox[-1]-zBox[0]
    LT_Box      = tBox[-1]-tBox[0]    
    LX_Box      = LT_Box * meanU

    # --- Desired resolution, rules of thumb
    dX_High_desired = chord_max
    dX_Low_desired  = Cmeander*D*meanU/150.0
    dY_Low_desired  = dX_Low_desired
    dZ_Low_desired  = dX_Low_desired
    dT_Low_desired  = Cmeander*D/(10.0*meanU)

    # --- Suitable resolution for high res
    dX_High = int(dX_High_desired/dX_Box)*dX_Box
    if dX_High==0: raise Exception('The x-resolution of the box ({}) is too large and cannot satisfy the requirements for the high-res domain of dX~{} (based on chord_max). Reduce DX (or DT) of the box.'.format(dX_Box, dX_High_desired))
    dY_High = dY_Box  # TODO? 
    dZ_High = dZ_Box  # TODO? 
    dT_High = dT_Box  # TODO? 

    # --- Suitable resolution for Low res
    dT_Low = int(dT_Low_desired/dT_Box )*dT_Box
    dX_Low = int(dX_Low_desired/dX_High)*dX_High
    dY_Low = int(dY_Low_desired/dY_High)*dY_High
    dZ_Low = int(dZ_Low_desired/dZ_High)*dZ_High
    if dT_Low==0: raise Exception('The time-resolution of the box ({}) is too large and cannot satisfy the requirements for the low-res domain of dT~{} (based on D & U). Reduce the DT of the box.'.format(dT_Box, dT_Low_desired))
    if dX_Low==0: raise Exception('The X-resolution of the box ({}) is too large and cannot satisfy the requirements for the low-res domain of dX~{} (based on D & U). Reduce the DX of the box.'.format(dX_Box, dX_Low_desired))
    if dY_Low==0: raise Exception('The Y-resolution of the box ({}) is too large and cannot satisfy the requirements for the low-res domain of dY~{} (based on D & U). Reduce the DY of the box.'.format(dY_Box, dY_Low_desired))
    if dZ_Low==0: raise Exception('The Z-resolution of the box ({}) is too large and cannot satisfy the requirements for the low-res domain of dZ~{} (based on D & U). Reduce the DZ of the box.'.format(dZ_Box, dZ_Low_desired))

    # --- Low-res domain
    # NOTE: more work is needed to make sure the domain encompass the turbines
    #       Also, we need to know the main flow direction to add a buffere with extent_wake
    # Origin
    nD_Before = extent_X/2 # Diameters before the first turbine to start the domain
    X0_Low = np.floor( (min(xWT)-nD_Before*D-dX_Low)) # Starting on integer value for esthetics. With a dX_Low margin.
    Y0_Low = np.floor( -LY_Box/2                    ) # Starting on integer value for esthetics
    Z0_Low = zBox[0] # we start at lowest to include tower
    if LES:
        if Y0_Low > min(yWT)-3*D:
            Y0_Low = np.floor(min(yWT)-3*D) 
    # Extent NOTE: this assumes main flow about x. Might need to be changed
    
    XMax_Low = max(xWT) + extent_wake*D
    LX_Low = XMax_Low-X0_Low
    LY_Low = LY_Box 
    LZ_Low = LZ_Box 
    # Number of points
    nX_Low = int(np.ceil(LX_Low/dX_Low))
    nY_Low = int(np.ceil(LY_Low/dY_Low))
    nZ_Low = int(np.ceil(LZ_Low/dZ_Low))
    # Make sure we don't exceed box in Y and Z #rt: this essentially gives us 1 less grid point than what is on the inp/bst files
    if (nY_Low*dY_Low>LY_Box): nY_Low=nY_Low-1 
    if (nZ_Low*dZ_Low>LZ_Box): nZ_Low=nZ_Low-1 

    # --- High-res domain extent and number of points
    ZMax_High = hubHeight+extent_YZ*D/2
    Z0_High   = zBox[0] # we start at lowest to include tower
    LX_High =  extent_X*D        
    LY_High =  min(LY_Box, extent_YZ*D      ) # Bounding to not exceed the box dimension
    LZ_High =  min(LZ_Box, ZMax_High-Z0_High) # Bounding to not exceed the box dimension
    nX_High = int(np.ceil(LX_High/dX_High))
    nY_High = int(np.ceil(LY_High/dY_High))
    nZ_High = int(np.ceil(LZ_High/dZ_High))
    # Make sure we don't exceed box in Y and Z
    if (nY_High*dY_High>LY_Box): nY_High=nY_High-1 
    if (nZ_High*dZ_High>LZ_Box): nZ_High=nZ_High-1 

    # --- High-res location per turbine 
    X0_desired = np.asarray(xWT)-LX_High/2 # high-res is centered on turbine location
    Y0_desired = np.asarray(yWT)-LY_High/2 # high-res is centered on turbine location
    X0_High    = X0_Low + np.floor((X0_desired-X0_Low)/dX_High)*dX_High
    Y0_High    = Y0_Low + np.floor((Y0_desired-Y0_Low)/dY_High)*dY_High

    d = dict()
    d['DT_Low']  = np.around(dT_Low ,4)
    d['DT_High'] = np.around(dT_High,4)
    d['NX_Low']  = nX_Low
    d['NY_Low']  = nY_Low
    d['NZ_Low']  = nZ_Low
    d['X0_Low']  = np.around(X0_Low,4)
    d['Y0_Low']  = np.around(Y0_Low,4)
    d['Z0_Low']  = np.around(Z0_Low,4)
    d['dX_Low']  = np.around(dX_Low,4)
    d['dY_Low']  = np.around(dY_Low,4)
    d['dZ_Low']  = np.around(dZ_Low,4)
    d['NX_High'] = nX_High
    d['NY_High'] = nY_High
    d['NZ_High'] = nZ_High
    # --- High extent info for turbine outputs
    d['dX_High'] = np.around(dX_High,4)
    d['dY_High'] = np.around(dY_High,4)
    d['dZ_High'] = np.around(dZ_High,4)
    d['X0_High'] = np.around(X0_High,4)
    d['Y0_High'] = np.around(Y0_High,4)
    d['Z0_High'] = np.around(Z0_High,4)
    # --- Misc
    d['dX_des_High'] = dX_High_desired
    d['dX_des_Low']  = dX_Low_desired
    d['DT_des']      = dT_Low_desired
    d['U_mean']      = meanU

    # --- Sanity check: check that the high res is at "almost" an integer location
    X_rel = (np.array(d['X0_High'])-d['X0_Low'])/d['dX_High']
    Y_rel = (np.array(d['Y0_High'])-d['Y0_Low'])/d['dY_High']
    dX = X_rel - np.round(X_rel) # Should be close to zero
    dY = Y_rel - np.round(Y_rel) # Should be close to zero
    if any(abs(dX)>1e-3):
        print('Deltas:',dX)
        print('Exception has been raise. I put this print statement instead. Check with EB.')
        print('Exception: Some X0_High are not on an integer multiple of the high-res grid')
        #raise Exception('Some X0_High are not on an integer multiple of the high-res grid')
    if any(abs(dY)>1e-3):
        print('Deltas:',dY)
        print('Exception has been raise. I put this print statement instead. Check with EB.')
        print('Exception: Some Y0_High are not on an integer multiple of the high-res grid')
        #raise Exception('Some Y0_High are not on an integer multiple of the high-res grid')

    return d


def writeFastFarm(outputFile, templateFile, xWT, yWT, zWT, FFTS=None, OutListT1=None, noLeadingZero=False):
    """ Write FastFarm input file based on a template, a TurbSimFile and the Layout
    
    outputFile: .fstf file to be written
    templateFile: .fstf file that will be used to generate the output_file
    XWT,YWT,ZWT: positions of turbines
    FFTS: FastFarm TurbSim parameters as returned by fastFarmTurbSimExtent
    """
    # --- Read template fast farm file
    fst=FASTInputFile(templateFile)
    # --- Replace box extent values
    if FFTS is not None:
        fst['Mod_AmbWind'] = 2
        ModVars = ['DT_Low', 'DT_High', 'NX_Low', 'NY_Low', 'NZ_Low', 'X0_Low', 'Y0_Low', 'Z0_Low', 'dX_Low', 'dY_Low', 'dZ_Low', 'NX_High', 'NY_High', 'NZ_High']
        for k in ModVars:
            if isinstance(FFTS[k],int):
                fst[k] = FFTS[k] 
            else:
                fst[k] = np.around(FFTS[k],3)
        fst['WrDisDT'] = FFTS['DT_Low']

    # --- Set turbine names, position, and box extent
    nWT = len(xWT)
    fst['NumTurbines'] = nWT
    if FFTS is not None:
        nCol= 10
    else:
        nCol = 4
    ref_path = fst['WindTurbines'][0,3]
    WT = np.array(['']*nWT*nCol,dtype='object').reshape((nWT,nCol))
    for iWT,(x,y,z) in enumerate(zip(xWT,yWT,zWT)):
        WT[iWT,0]=x
        WT[iWT,1]=y
        WT[iWT,2]=z
        WT[iWT,3]=insertTN(ref_path,iWT+1,nWT,noLeadingZero=noLeadingZero)
        if FFTS is not None:
            WT[iWT,4]=FFTS['X0_High'][iWT]
            WT[iWT,5]=FFTS['Y0_High'][iWT]
            WT[iWT,6]=FFTS['Z0_High']
            WT[iWT,7]=FFTS['dX_High']
            WT[iWT,8]=FFTS['dY_High']
            WT[iWT,9]=FFTS['dZ_High']
    fst['WindTurbines']=WT

    fst.write(outputFile)
    if OutListT1 is not None:
        setFastFarmOutputs(outputFile, OutListT1)

def setFastFarmOutputs(fastFarmFile, OutListT1):
    """ Duplicate the output list, by replacing "T1" with T1->Tn """
    fst = FASTInputFile(fastFarmFile)
    nWTOut = min(fst['NumTurbines'],9) # Limited to 9 turbines
    OutList=['']
    for s in OutListT1:
        s=s.strip('"')  
        if 'T1' in s:
            OutList+=['"'+s.replace('T1','T{:d}'.format(iWT+1))+'"' for iWT in np.arange(nWTOut) ]
        elif 'W1VAmb' in s: # special case for ambient wind
            OutList+=['"'+s.replace('1','{:d}'.format(iWT+1))+'"' for iWT in np.arange(nWTOut) ]
        elif 'W1VDis' in s: # special case for disturbed wind
            OutList+=['"'+s.replace('1','{:d}'.format(iWT+1))+'"' for iWT in np.arange(nWTOut) ]
        else:
            OutList+='"'+s+'"'
    fst['OutList']=OutList
    fst.write(fastFarmFile)


def plotFastFarmSetup(fastFarmFile, grid=True, fig=None, D=None, plane='XY', hubHeight=None, showLegend=True):
    """ """
    import matplotlib.pyplot as plt

    def col(i): 
        Colrs=plt.rcParams['axes.prop_cycle'].by_key()['color']
        return Colrs[ np.mod(i,len(Colrs)) ]
    def boundingBox(x, y):
        """ return x and y coordinates to form a box marked by the min and max of x and y"""
        x_bound = [x[0],x[-1],x[-1],x[0] ,x[0]]
        y_bound = [y[0],y[0] ,y[-1],y[-1],y[0]]
        return x_bound, y_bound


    # --- Read FAST.Farm input file
    fst=FASTInputFile(fastFarmFile)

    if fig is None:
        fig = plt.figure(figsize=(13.5,8))
        ax  = fig.add_subplot(111,aspect="equal")

    WT=fst['WindTurbines']
    xWT = WT[:,0].astype(float)
    yWT = WT[:,1].astype(float)
    zWT = yWT*0 
    if hubHeight is not None:
        zWT += hubHeight

    if plane == 'XY':
        pass
    elif plane == 'XZ':
        yWT = zWT
    elif plane == 'YZ':
        xWT = yWT
        yWT = zWT
    else:
        raise Exception("Plane should be 'XY' 'XZ' or 'YZ'")

    if fst['Mod_AmbWind'] == 2:
        x_low = fst['X0_Low'] + np.arange(fst['NX_Low']+1)*fst['DX_Low']
        y_low = fst['Y0_Low'] + np.arange(fst['NY_Low']+1)*fst['DY_Low']
        z_low = fst['Z0_Low'] + np.arange(fst['NZ_Low']+1)*fst['DZ_Low']
        if plane == 'XZ':
            y_low = z_low
        elif plane == 'YZ':
            x_low = y_low
            y_low = z_low
        # Plot low-res box
        x_bound_low, y_bound_low =  boundingBox(x_low, y_low)
        ax.plot(x_bound_low, y_bound_low ,'--k',lw=2,label='Low-res')
        # Plot Low res grid lines
        if grid:
            ax.vlines(x_low, ymin=y_low[0], ymax=y_low[-1], ls='-', lw=0.3, color=(0.3,0.3,0.3))
            ax.hlines(y_low, xmin=x_low[0], xmax=x_low[-1], ls='-', lw=0.3, color=(0.3,0.3,0.3))

        X0_High = WT[:,4].astype(float)
        Y0_High = WT[:,5].astype(float)
        Z0_High = WT[:,6].astype(float)
        dX_High = WT[:,7].astype(float)[0]
        dY_High = WT[:,8].astype(float)[0]
        dZ_High = WT[:,9].astype(float)[0]
        nX_High = fst['NX_High']
        nY_High = fst['NY_High']
        nZ_High = fst['NZ_High']

        # high-res boxes
        for wt in range(len(xWT)):
            x_high = X0_High[wt] + np.arange(nX_High+1)*dX_High
            y_high = Y0_High[wt] + np.arange(nY_High+1)*dY_High
            z_high = Z0_High[wt] + np.arange(nZ_High+1)*dZ_High
            if plane == 'XZ':
                y_high = z_high
            elif plane == 'YZ':
                x_high = y_high
                y_high = z_high

            x_bound_high, y_bound_high =  boundingBox(x_high, y_high)
            ax.plot(x_bound_high, y_bound_high, '-', lw=2, c=col(wt))
            # Plot High res grid lines
            if grid:
                ax.vlines(x_high, ymin=y_high[0], ymax=y_high[-1], ls='--', lw=0.4, color=col(wt))
                ax.hlines(y_high, xmin=x_high[0], xmax=x_high[-1], ls='--', lw=0.4, color=col(wt))

    # Plot turbines
    for wt in range(len(xWT)):
        ax.plot(xWT[wt], yWT[wt], 'x', ms=8, mew=2, c=col(wt),label="WT{}".format(wt+1))
        if plane=='XY' and D is not None:
            ax.plot([xWT[wt],xWT[wt]], [yWT[wt]-D/2,yWT[wt]+D/2], '-', lw=2, c=col(wt))
        elif plane=='XZ' and D is not None and hubHeight is not None:
            ax.plot([xWT[wt],xWT[wt]], [yWT[wt]-D/2,yWT[wt]+D/2], '-', lw=2, c=col(wt))
        elif plane=='YZ' and D is not None and hubHeight is not None:
            theta = np.linspace(0,2*np.pi, 40)
            x = xWT[wt] + D/2*np.cos(theta)
            y = yWT[wt] + D/2*np.sin(theta)
            ax.plot(x, y, '-', lw=2, c=col(wt))

    #plt.legend(bbox_to_anchor=(1.05,1.015),frameon=False)
    if showLegend:
        ax.legend()
    if plane=='XY':
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")
    elif plane=='XZ':
        ax.set_xlabel("x [m]")
        ax.set_ylabel("z [m]")
    elif plane=='YZ':
        ax.set_xlabel("y [m]")
        ax.set_ylabel("z [m]")
    fig.tight_layout
    # fig.savefig('FFarmLayout.pdf',bbox_to_inches='tight',dpi=500)

    return fig

# --------------------------------------------------------------------------------}
# --- Tools for postpro 
# --------------------------------------------------------------------------------{

def spanwiseColFastFarm(Cols, nWT=9, nD=9):
    """ Return column info, available columns and indices that contain AD spanwise data"""
    FFSpanMap=dict()
    for i in np.arange(nWT):
        FFSpanMap[r'^CtT{:d}N(\d*)_\[-\]'.format(i+1)]='CtT{:d}_[-]'.format(i+1)
    for i in np.arange(nWT):
        for k in np.arange(nD):
            FFSpanMap[r'^WkDfVxT{:d}N(\d*)D{:d}_\[m/s\]'.format(i+1,k+1) ]='WkDfVxT{:d}D{:d}_[m/s]'.format(i+1, k+1)  
    for i in np.arange(nWT):
        for k in np.arange(nD):
            FFSpanMap[r'^WkDfVrT{:d}N(\d*)D{:d}_\[m/s\]'.format(i+1,k+1) ]='WkDfVrT{:d}D{:d}_[m/s]'.format(i+1, k+1)  

    return fastlib.find_matching_columns(Cols, FFSpanMap)

def diameterwiseColFastFarm(Cols, nWT=9):
    """ Return column info, available columns and indices that contain AD spanwise data"""
    FFDiamMap=dict()
    for i in np.arange(nWT):
        for x in ['X','Y','Z']:
            FFDiamMap[r'^WkAxs{}T{:d}D(\d*)_\[-\]'.format(x,i+1)]   ='WkAxs{}T{:d}_[-]'.format(x,i+1) 
    for i in np.arange(nWT):
        for x in ['X','Y','Z']:
            FFDiamMap[r'^WkPos{}T{:d}D(\d*)_\[m\]'.format(x,i+1)]   ='WkPos{}T{:d}_[m]'.format(x,i+1)
    for i in np.arange(nWT):
        for x in ['X','Y','Z']:
            FFDiamMap[r'^WkVel{}T{:d}D(\d*)_\[m/s\]'.format(x,i+1)] ='WkVel{}T{:d}_[m/s]'.format(x,i+1) 
    for i in np.arange(nWT):
        for x in ['X','Y','Z']:
            FFDiamMap[r'^WkDiam{}T{:d}D(\d*)_\[m\]'.format(x,i+1)]  ='WkDiam{}T{:d}_[m]'.format(x,i+1)
    return fastlib.find_matching_columns(Cols, FFDiamMap)

def SensorsFARMRadial(nWT=3,nD=10,nR=30,signals=None):
    """ Returns a list of FASTFarm sensors that are used for the radial distribution
    of quantities (e.g. Ct, Wake Deficits).
    If `signals` is provided, the output is the list of sensors within the list `signals`.
    """
    WT  = np.arange(nWT)
    r   = np.arange(nR)
    D   = np.arange(nD)
    sens=[]
    sens+=['CtT{:d}N{:02d}_[-]'.format(i+1,j+1) for i in WT for j in r]
    sens+=['WkDfVxT{:d}N{:02d}D{:d}_[m/s]'.format(i+1,j+1,k+1) for i in WT for j in r for k in D]
    sens+=['WkDfVrT{:d}N{:02d}D{:d}_[m/s]'.format(i+1,j+1,k+1) for i in WT for j in r for k in D]
    if signals is not None:
        sens = [c for c in sens if c in signals]
    return sens

def SensorsFARMDiam(nWT,nD):
    """ Returns a list of FASTFarm sensors that contain quantities at different downstream diameters
     (e.g. WkAxs, WkPos, WkVel, WkDiam)
    If `signals` is provided, the output is the list of sensors within the list `signals`.
    """
    WT  = np.arange(nWT)
    D   = np.arange(nD)
    XYZ = ['X','Y','Z']
    sens=[]
    sens+=['WkAxs{}T{:d}D{:d}_[-]'.format(x,i+1,j+1) for x in XYZ for i in WT for j in D]
    sens+=['WkPos{}T{:d}D{:d}_[m]'.format(x,i+1,j+1) for x in XYZ for i in WT for j in D]
    sens+=['WkVel{}T{:d}D{:d}_[m/s]'.format(x,i+1,j+1) for x in XYZ for i in WT for j in D]
    sens+=['WkDiam{}T{:d}D{:d}_[m]'.format(x,i+1,j+1) for x in XYZ for i in WT for j in D]
    if signals is not None:
        sens = [c for c in sens if c in signals]
    return sens


def extractFFRadialData(fastfarm_out,fastfarm_input,avgMethod='constantwindow',avgParam=30,D=1,df=None):
    # LEGACY
    return spanwisePostProFF(fastfarm_input,avgMethod=avgMethod,avgParam=avgParam,D=D,df=df,fastfarm_out=fastfarm_out)


def spanwisePostProFF(fastfarm_input,avgMethod='constantwindow',avgParam=30,D=1,df=None,fastfarm_out=None):
    """ 
    Opens a FASTFarm output file, extract the radial data, average them and returns spanwise data

    D: diameter TODO, extract it from the main file

    See faslibt.averageDF for `avgMethod` and `avgParam`.
    """
    # --- Opening ouputfile
    if df is None:
        df=FASTOutputFile(fastfarm_out).toDataFrame()

    # --- Opening input file and extracting inportant variables
    if fastfarm_input is None:
        # We don't have an input file, guess numbers of turbine, diameters, Nodes...
        cols, sIdx = fastlib.find_matching_pattern(df.columns.values, r'T(\d+)')
        nWT = np.array(sIdx).astype(int).max()
        cols, sIdx = fastlib.find_matching_pattern(df.columns.values, r'D(\d+)')
        nD = np.array(sIdx).astype(int).max()
        cols, sIdx = fastlib.find_matching_pattern(df.columns.values, r'N(\d+)')
        nr = np.array(sIdx).astype(int).max()
        vr=None
        vD=None
        D=0
        main = None
    else:
        main=FASTInputFile(fastfarm_input)
        iOut    = main['OutRadii']
        dr      = main['dr']              # Radial increment of radial finite-difference grid (m)
        OutDist = main['OutDist']         # List of downstream distances for wake output for an individual rotor
        WT     = main['WindTurbines']
        nWT    = len(WT)
        vr     = dr*np.array(iOut)
        vD     = np.array(OutDist)
        nr=len(iOut)
        nD=len(vD)


    # --- Extracting time series of radial data only
    colRadial = SensorsFARMRadial(nWT=nWT,nD=nD,nR=nr,signals=df.columns.values)
    colRadial=['Time_[s]']+colRadial
    dfRadialTime = df[colRadial] # TODO try to do some magic with it, display it with a slider

    # --- Averaging data
    dfAvg = fastlib.averageDF(df,avgMethod=avgMethod,avgParam=avgParam)

    # --- Extract radial data
    ColsInfo, nrMax = spanwiseColFastFarm(df.columns.values, nWT=nWT, nD=nD)
    dfRad        = fastlib.extract_spanwise_data(ColsInfo, nrMax, df=None, ts=dfAvg.iloc[0])
    #dfRad       = fastlib.insert_radial_columns(dfRad, vr)
    if dfRad is not None:
        dfRad.insert(0, 'i_[#]', np.arange(nrMax)+1) # For all, to ease comparison
        if vr is not None: 
            dfRad.insert(0, 'r_[m]', vr[:nrMax]) # give priority to r_[m] when available
        dfRad['i/n_[-]']=np.arange(nrMax)/nrMax

    # --- Extract downstream data
    ColsInfo, nDMax = diameterwiseColFastFarm(df.columns.values, nWT=nWT)
    dfDiam       = fastlib.extract_spanwise_data(ColsInfo, nDMax, df=None, ts=dfAvg.iloc[0])
    if dfDiam is not None:
        dfDiam.insert(0, 'i_[#]', np.arange(nDMax)+1) # For all, to ease comparison
        if vD is not None:
            dfDiam.insert(0, 'x_[m]', vD[:nDMax])
        dfDiam['i/n_[-]'] = np.arange(nDMax)/nDMax
    return dfRad, dfRadialTime, dfDiam

