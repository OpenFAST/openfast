"""
Tools to work with OLAF the vortex code implemented in openfast
"""
import numpy as np


def OLAFParams(omega_rpm, U0, R, a=0.3, aScale=1.2,
          deltaPsiDeg=6, nPerRot=None,
          targetFreeWakeLengthD=1,
          targetWakeLengthD=4.,
          nNWrot=8, 
          nNWrotFree=1,
          nFWrot=0, nFWrotFree=0,
          verbose=True, dt_glue_code=None):
    """
    Computes recommended time step and wake length for OLAF based on:

    INPUTS:
     - omega_rpm: rotational speed [RPM]
     - U0: mean wind speed [m/s]
     - R: rotor radius [m]

    OPTIONS FOR TIME STEP:
      - either:
         - deltaPsiDeg : target azimuthal discretization [deg]
              or
         - nPerRot     : number of time step per rotations.
                deltaPsiDeg  -  nPerRot
                     5            72
                     6            60
                     7            51.5
                     8            45
     - dt_glue_code: glue code time step. If provided, the time step of OLAF will be approximated
                     such that it is a multiple of the glue-code time step.

    OPTIONS FOR WAKE LENGTH:
     - a: average axial induction factor at the rotor [-]
     - aScale: scaling factor to estimate induction, such that the wake convection velocity is:
               Uc=U0(1-aScale*a)
     - targetWakeLengthD: target wake length in diameter [D]
     - nNWrot     : minimum number of near wake rotations
     - nFWrot     : minimum number of far wake rotations
     - nFWrotFree : minimum number of far wake rotations (free panels)

    """
    def myprint(*args, **kwargs):
        if verbose:
            print(*args, **kwargs)

    # Rotational period
    omega = omega_rpm*2*np.pi/60
    T = 2*np.pi/omega
    # Convection velocity
    Uc = U0 * (1-aScale*a)

    # Desired time step
    if nPerRot is not None:
        dt_wanted    = np.around(T/nPerRot,5)
        deltaPsiDeg  = np.around(omega*dt_wanted*180/np.pi ,2)
    else:
        dt_wanted    = np.around(deltaPsiDeg/(6*omega_rpm),5)
        nPerRot = int(2*np.pi /(deltaPsiDeg*np.pi/180))

    # Adapting desired time step based on glue code time step
    if dt_glue_code is not None:
        dt_rounded = round(dt_wanted/dt_glue_code)*dt_glue_code
        deltaPsiDeg2 = np.around(omega*dt_rounded *180/np.pi ,2)
        myprint('>>> To satisfy glue-code dt:')
        myprint('    Rounding dt   from {} to {}'.format(dt_wanted, dt_rounded    ))
        myprint('    Changing dpsi from {} to {}'.format(deltaPsiDeg, deltaPsiDeg2))
        dt_fvw   = dt_rounded
        deltaPsiDeg = deltaPsiDeg2
        nPerRot = int(2*np.pi /(deltaPsiDeg*np.pi/180))
    else:
        dt_fvw = dt_wanted

    # Useful functions
    n2L = lambda n: (n * dt_fvw * Uc)/(2*R)  # convert number of panels to distance
    n2R = lambda n:  n * dt_fvw / T          # convert number of panels to number of rotations
    n2t = lambda n:  n * dt_fvw              # convert number of panels to time

    # All Wake (AW) panels - Wake length from mean wind speed
    targetWakeLength = targetWakeLengthD * 2 * R
    nAWPanels_FromU0 = int(targetWakeLength / (Uc*dt_fvw))
    # Free near wake panels (based on distance)
    targetFreeWakeLength = targetFreeWakeLengthD * 2 * R
    nNWPanelsFree_FromU0 = int(targetFreeWakeLength / (Uc*dt_fvw))
    nNWPanelsFree_FromRot = int(nNWrotFree*nPerRot)
    # Far wake (FW) panels, always from number of rotations
    nFWPanels     = int(nFWrot*nPerRot)
    nFWPanelsFree = int(nFWrotFree*nPerRot)
    # Wake length from rotational speed and number of rotations
    nAWPanels_FromRot = int(nNWrot*nPerRot) # Total number of panels NW+FW

    # Below we chose between criteria on number of rotation or donwstream distance
    # This can be adapted/improved
    myprint('Number of panels (NW free) from wind speed and distance:{:15d}'.format(nNWPanelsFree_FromU0))
    myprint('Number of panels (NW free) from number of rotations    :{:15d}'.format(nNWPanelsFree_FromRot))
    myprint('Number of panels (NW+FW)   from wind speed and distance:{:15d}'.format(nAWPanels_FromU0))
    myprint('Number of panels (NW+FW)   from number of rotations    :{:15d}'.format(nAWPanels_FromRot))
    myprint('Number of panels (NW+FW)   from average between two    :{:15d}'.format(int((nAWPanels_FromRot+nAWPanels_FromU0)/2)))
    if nAWPanels_FromRot>nAWPanels_FromU0:
        # Criteria based on rotation wins:
        myprint('[INFO] Using number of rotations to setup number of panels')
        nAWPanels = nAWPanels_FromRot # Total number of panels NW+FW
    else:
        myprint('[INFO] Using wind speed and distance to setup number of panels')
        # Wake distance wins, we keep the nFW from rot but increase nNW
        nAWPanels = nAWPanels_FromU0  # Total number of panels NW+FW
    if nNWPanelsFree_FromRot>nNWPanelsFree_FromU0:
        # Criteria based on rotation wins:
        myprint('[INFO] Using number of rotations to setup number of free panels')
        nNWPanelsFree = nNWPanelsFree_FromRot
    else:
        myprint('[INFO] Using wind speed and distance to setup number of free panels')
        # Wake distance wins, we keep the nFW from rot but increase nNW
        nNWPanelsFree = nNWPanelsFree_FromU0

    nNWPanels = nAWPanels - nFWPanels # nNW = All-Far Wake

    # See "free" near wake
    if nNWPanelsFree>nNWPanels:
        nNWPanelsFree=nNWPanels
        myprint('[INFO] Capping number of free NW panels to max.')
    if nNWPanelsFree<nNWPanels and nFWPanelsFree>0:
        nFWPanelsFree=0
        myprint('[INFO] Setting number of Free FW panels to zero because a frozen near wake is used')

    # Transient time (twice the time to develop the full wake extent)
    # This is the minimum recommended time before convergence of the wake is expected
    # (might be quite long)
    tMin = 2 * dt_fvw*nAWPanels
    if verbose:
        myprint('')
        myprint('{:15.2f} Transient time   ({:5.1f} rot)'.format(tMin, tMin/T))
        myprint('{:15d} nAWPanels        ({:5.1f} rot, {:5.1f}D)'.format(nAWPanels, n2R(nAWPanels), n2L(nAWPanels)))
        myprint('')
        myprint('OLAF INPUT FILE:')
        myprint('----------------------- GENERAL OPTIONS ---------------------')
        myprint('{:15.6f} DTFVW        (delta psi = {:5.1f}deg)'.format(dt_fvw, deltaPsiDeg))
        myprint('--------------- WAKE EXTENT AND DISCRETIZATION --------------')
        myprint('{:15d} nNWPanels     ({:5.1f} rot, {:5.1f}D)'.format(nNWPanels    , n2R(nNWPanels    ), n2L(nNWPanels    )))
        myprint('{:15d} nNWPanelsFree ({:5.1f} rot, {:5.1f}D {:5.1f})'.format(nNWPanelsFree, n2R(nNWPanelsFree), n2L(nNWPanelsFree), n2t(nNWPanelsFree)))
        myprint('{:15d} nFWPanels     ({:5.1f} rot, {:5.1f}D)'.format(nFWPanels    , n2R(nFWPanels    ), n2L(nFWPanels    )))
        myprint('{:15d} nFWPanelsFree ({:5.1f} rot, {:5.1f}D)'.format(nFWPanelsFree, n2R(nFWPanelsFree), n2L(nFWPanelsFree)))

    return dt_fvw, tMin, nNWPanels, nNWPanelsFree, nFWPanels, nFWPanelsFree



def OLAFParamsRPM(omega_rpm, deltaPsiDeg=6, nNWrot=2, nFWrot=10, nFWrotFree=3, nPerRot=None,  verbose=True, dt_glue_code=None, totalRot=None):
    """ 
    Computes recommended time step and wake length based on the rotational speed in RPM

    INPUTS:
     - omega_rpm: rotational speed in RPM
     - deltaPsiDeg : azimuthal discretization in deg
     - nNWrot : number of near wake rotations
     - nFWrot : total number of far wake rotations
     - nFWrotFree : number of far wake rotations that are free

        deltaPsiDeg  -  nPerRot
             5            72    
             6            60    
             7            51.5  
             8            45    
    """
    omega_rpm = np.asarray(omega_rpm)
    omega = omega_rpm*2*np.pi/60
    T = 2*np.pi/omega
    if nPerRot is not None:
        dt_wanted    = np.around(T/nPerRot,5)
        deltaPsiDeg  = np.around(omega*dt_wanted*180/np.pi ,2)
    else:
        dt_wanted    = np.around(deltaPsiDeg/(6*omega_rpm),5)
        nPerRot = int(2*np.pi /(deltaPsiDeg*np.pi/180))
    if dt_glue_code is not None:
        dt_rounded = round(dt_wanted/dt_glue_code)*dt_glue_code
        deltaPsiDeg2 = np.around(omega*dt_rounded *180/np.pi ,2)
        print('>>> To satisfy glue-code dt:')
        print('    Rounding dt   from {} to {}'.format(dt_wanted, dt_rounded    ))
        print('    Changing dpsi from {} to {}'.format(deltaPsiDeg, deltaPsiDeg2))
        dt_wanted   = dt_rounded
        deltaPsiDeg = deltaPsiDeg2
        nPerRot = int(2*np.pi /(deltaPsiDeg*np.pi/180))

    nNWPanel     = int(nNWrot*nPerRot)
    nFWPanel     = int(nFWrot*nPerRot)
    nFWPanelFree = int(nFWrotFree*nPerRot)

    if totalRot is None:
        totalRot = (nNWrot + nFWrot)*3 # going three-times through the entire wake

    tMax = dt_wanted*nPerRot*totalRot

    if verbose:
        print(dt_wanted              , '  dt')
        print(int      (nNWPanel    ), '  nNWPanel          ({} rotations)'.format(nNWrot))
        print(int      (nFWPanel    ), '  FarWakeLength     ({} rotations)'.format(nFWrot))
        print(int      (nFWPanelFree), '  FreeFarWakeLength ({} rotations)'.format(nFWrotFree))
        print(tMax              , '  Tmax ({} rotations)'.format(totalRot))

    return dt_wanted, tMax, nNWPanel, nFWPanel, nFWPanelFree






if __name__ == '__main__':
    OLAFParams(omega_rpm = 4.87558, U0=8, R=63, deltaPsiDeg=6, verbose=True)
