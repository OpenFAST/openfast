""" 
Example usage for the TurbSimFile class.

- Read a TurbSim file and display main properties
- Extract time series at a given y, z location and plot it
- Extract a horizontal plane and plot it
- Compute vertical profile/shear and plot it
- Fit vertical profiel with a power law
- Compute cross corelation in y and z directions
- Modify field (add a constant velocity in the streamwise direction) and write to disk
- Write to Mann Box format
- Write time series at given locations to a CSV file

NOTE: this example uses an extremely small TurbSim box. 
      Results will be more "physical" on a more realstic box.
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from pyFAST.input_output import TurbSimFile

def main():
    MyDir = os.path.dirname(__file__)

    # --- Read a TurbSim file and display main properties
    filename = os.path.join(MyDir, '../tests/example_files/TurbSim_WithTwr.bts')
    ts = TurbSimFile(filename)
    print(ts)

    # --- Extract time series at a given y, z location and plot it
    # Method 1 - use object method
    u,v,w = ts.valuesAt(y=0, z=90, method='nearest')
    # Method 2 - use data directly
    iy, iz = ts.closestPoint(y=0, z=90)
    u2,v2,w2 = ts['u'][0, :, iy, iz], ts['u'][1, :, iy, iz], ts['u'][2, :, iy, iz]

    fig,ax = plt.subplots(1, 1)
    ax.plot(ts.t, u, label='u')
    ax.plot(ts.t, v, label='v')
    ax.plot(ts.t, w, label='w')
    ax.plot(ts.t, u2, 'k--')
    ax.plot(ts.t, v2, 'k--')
    ax.plot(ts.t, w2, 'k--')
    ax.legend()
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Velocity [m/s]')
    ax.set_title('Velocity at y=0 z=90')

    # --- Extract a horizontal plane and plot it
    U, V, W = ts.horizontalPlane(z=90)
    T, Y = np.meshgrid(ts.t, ts.y)
    fig,ax = plt.subplots(1, 1)
    ax.contourf(T, Y, U.T)
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('y [m]')
    ax.set_title('Velocity at z=90')

    # --- Compute vertical profile/shear and plot it
    # NOTE: the example file has only three points in y&z
    z, u_mean, u_std = ts.vertProfile(y_span='full')
    # Fit a power law
    u_fit, pfit, model, z_ref =  ts.fitPowerLaw()
    print('Power law: alpha={:.5f}, u_ref={:.5f}, z_ref={:.5f}'.format(pfit[1],pfit[0],z_ref))
    print('Formula: {} '.format(model['formula']))

    fig,ax = plt.subplots(1, 1)
    ax.plot(u_mean[0,:], z, label='u')
    ax.plot(u_mean[1,:], z, label='v')
    ax.plot(u_mean[2,:], z, label='w')
    ax.plot(u_fit      , z, 'k--', label='u fit (power law)')
    if 'uTwr' in ts:
        ax.plot(np.mean(ts['uTwr'][0,:,:], axis=0), ts['zTwr'], label='u on tower')
    ax.legend()
    ax.set_xlabel('Velocity [m/s]')
    ax.set_ylabel('z [m]')
    ax.set_title('Vertical profiles (averaged over y and time)')





    # --- Compute cross corelation in y and z directions
    # NOTE: the example file has only three points in y&z
    iy0, iz0 = ts.iMid # Index at middle of the box
    y, rho_uu_y, rho_vv_y, rho_ww_y = ts.crosscorr_y(iy0, iz0)
    z, rho_uu_z, rho_vv_z, rho_ww_z = ts.crosscorr_z(iy0, iz0)

    fig,ax = plt.subplots(1, 1)
    ax.plot(y, rho_uu_y, label=r'rho_uu}')
    ax.plot(y, rho_vv_y, label=r'rho_vv}')
    ax.plot(y, rho_ww_y, label=r'rho_ww}')
    ax.set_xlabel('y [m]')
    ax.set_ylabel('Cross correlation')
    ax.set_title('Cross correlation in y direction at middle of the box')
    ax.legend()

    # --- Convert to "DataFrame" 
    # Contains relevant time series like vertical profile, midline, coherence, cross correlation
    #    dfs = ts.toDataFrame()

    # --- Modify field and write to disk
    ts['u'][0,:,:,:] += 1 # Adding 1 m/s in the streamwise
    ts.makePeriodic()     # Make the field periodic by mirroring it in the streamwise direction
    ts.write('_MyNewTurbBox.bts')

    # --- Write to Mann Box format
    ts.toMannBox()

    # --- Write time series at given locations to a CSV file
    ts.writeProbes('_Probes.csv', yProbe=[0], zProbe=[65,115])

if __name__ == "__main__":
    main()
    plt.show()

if __name__=='__test__':
    main()
    try:
        os.remove('_MyNewTurbBox.bts')
        os.remove('_MyNewTurbBox_198x3x4.u')
        os.remove('_MyNewTurbBox_198x3x4.v')
        os.remove('_MyNewTurbBox_198x3x4.w')
        os.remove('_Probes.csv')
    except:
        pass
