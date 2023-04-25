Working with OLAF
=================


.. _Running-OLAF:

Running OLAF
~~~~~~~~~~~~

As OLAF is a module of OpenFAST, the process of downloading, compiling,
and running OLAF is the same as that for OpenFAST. Such instructions are
available in the :ref:`installation` documentation.

.. note::
   To improve the speed of FVW module, the user may wish to compile with
   `OpenMP`.  To do so, add the `-DOPENMP=ON` option with CMake.


.. _Guidelines-OLAF:

Guidelines
~~~~~~~~~~

Most options of OLAF can be left to default. The results will depend on the time discretization, wake length, and regularization parameters. 
We provide guidelines for these parameters in this section, together with a simple python code to compute these parameters.
Please check this section again as we might further refine our guidelines with time.

We provide a python script at the end of the this section to programmatically set the main parameters.


**Time step**
We recommend to set OLAF's time step (**DTfvw**) such that it corresponds to :math:`\Delta \psi = 6` degrees of one rotor revolution:

.. math::
   
    \Delta t_\text{FVW}
    = \frac{\Delta \psi_\text{rad}}{\Omega_\text{rad/s}}
    = \frac{\Delta \psi_\text{deg}}{6 \times \Omega_\text{RPM}}

If the structural solver requires a smaller time step, the time step for the glue code can be set to a different value than **DTfvw** as long as **DTfvw** is a multiple of the glue code time step.



**Wake length and number of panels**

Each vortex element in the wake contributes to the induced velocity on the lifting line of the rotor. 
The longer the wake, the more contributions will be present. 
There is yet a trade-off to reach, as an increased wake length leads to more vortex elements, and therefore a higher computational time. 


If the wake consisted of a vortex cylinder of constant vorticity distribution, it can be shown that a wake length of 4D, 5D and 6D provide 99.2%, 99.5% and 99.7% of the induced velocity respectively. 
We therefore recommend to have a total wake length equal to at least 4D.
As an approximation, the distance  convected by the wake length is a function of the mean wind speed :math:`U_0` and the induced velocity within the wake. 
The induced velocity vary with downstream distance from :math:`-aU_0` at the rotor, where :math:`a` is the average axial induction factor at the rotor, to :math:`-2aU_0` once the wake has reached equilibrium (according to momentum theory). 
As an approximation, it can be assumed that the convection velocity is :math:`U_c = U_0(1-k_a a)`, with :math:`k_a\approx1.2`. 
We note that viscous diffusion is not accounted for, so for a simulation without turbulence, the wake convection velocity is not expected to recover to the freestream, 
and a "large" value of :math:`k_a` should be used (e.g. :math:`k_a\approx1.5`)
For a simulation with turbulence, meandering will diffuse the wake and a small value of :math:`k_a`  should be used (e.g. :math:`k_a\approx1.0`).
The axial induction factor is a function of the operating condition and design of the rotor. For estimates below, we will use :math:`a\approx1/3`.
The approximate time needed for the wake to reach a desired downstream distance :math:`d_\text{target}` is therefore:


.. math:: 
   
    T_\text{target} = \frac{d_\text{target}}{U_0(1-k_a a)}

This time corresponds to a number of time steps (i.e. the total number of wake panels) equal to:

.. math:: :label: ntargetDist

   
    n_\text{target,distance} 
       =  \frac{T_\text{target}}{\Delta t_\text{FVW}} 
       =  \frac{d_\text{target}}{U_0(1-k_a a) \Delta t_\text{FVW}} 
       \approx \frac{d_\text{target}}{0.5 U_0 \Delta t_\text{FVW}} 
       \approx \frac{12 D}{U_0 \Delta t_\text{FVW}} 
       \qquad \text{(integer)}

where the first approximation uses :math:`k_a a\approx 0.5` and the second approximation assumes a target distance of 6D.
It is also possible to define the number of near-wake panels based on a total number of revolution, :math:`n_\text{rot}`, leading to:

.. math:: :label: ntargetRot
   
    n_\text{target, rotations} 
       = \frac{n_\text{rot} T_\text{rot}}{dt_\text{FVW}}
       = \frac{n_\text{rot} 2 \pi }{\Omega dt_\text{FVW}}
       = \frac{n_\text{rot} 60 }{\Omega_\text{RPM} dt_\text{FVW}}
       \qquad \text{(integer)}


The wake of OLAF consists of two regions defined as "near wake" (NW) and "far wake" (FW), where the far consists of rolled-up tip and root vortices. 
The far wake has reduced accuracy, therefore  velocity profiles crossing the far wake will likely be inaccurate. 
The role of the far wake is to reduce the computational time. For increased accuracy, it is recommended to use a longer near wake and a shorter far wake. 
The far wake may be removed altogether.

The near wake and far wake are further decomposed into two regions, a region where the vortex filaments are free to convect, and another one where the filaments convect with a frozen, averaged induced velocity. 
The advantage of having a "frozen" wake region, is that it mitigates the impact of wake truncation which is an erroneous boundary condition (vortex lines cannot end in the fluid). If the wake is truncated while still being "free", then the vorticity will rollup in this region. 
Another advantage is that in the absence of diffusion, the wake tends to become excessively distorted downstream, reaching limit on the validity of the vortex filament representation. 
It is therefore useful to have a frozen far-wake region. 
The total number of wake panels is equal to the number of free near-wake panels, frozen near-wake panels, free far-wake panels and frozen far-wake panels:

.. math::
   
    n_\text{target} = n_\text{NW} + n_\text{FW} = n_\text{NW,Free} + n_\text{NW,Frozen} + n_\text{FW,Free} + n_\text{FW,Frozen}

OLAF input file defines
:math:`n_\text{NW}`      (**nNWPanels**),
:math:`n_\text{NW,Free}` (**nNWPanelsFree**),
:math:`n_\text{FW}`      (**nFWPanels**), and
:math:`n_\text{FW,Free}` (**nFWPanelsFree**).

If a "frozen" near-wake region is used (:math:`n_\text{NW,Free}<n_\text{NW}`) then the "free" far-wake region needs to be of zero length (:math:`n_\text{FW,Free}=0`).


We currently recommend:

- a total wake length of at least 4D (see Eq. :eq:`ntargetDist`), and
- a total wake length corresponding to at least 10 rotor revolutions (see Eq. :eq:`ntargetRot`)
- (depending on the operating conditions, one of the two conditions above will dominate, the largest wake length between the two is used).
- a near-wake extent corresponding to at least 8 revolutions
- a free near-wake extent corresponding to at least 1D
- a zero far-wake extent, or a frozen far-wake extent corresponding to 2 revolutions

The python script provided at the end of this section implements these guidelines.

General considerations:

- If you obtain a power coefficient above the Betz limit, most likely your number of panels is not set correctly according to these guidelines (number of panels too small). 
- If a far wake is used, do not set it as "free" for more than half of the length (i.e. nFWPanelsFree =< nFWPanels/2)
- The near wake is the most accurate. If computational time is not much of a concern, a long near wake is preferred, with a short or no far wake.
- For now, it's recommended to always have a frozen near wake or frozen far wake to avoid the error introduced by the wake truncation. 
- Wake velocity profiles may be erroneous within the far wake. 
- Write the wake (**WrVTK=1, or 2**) at regular interval (e.g. every 100 time steps) for visual inspection



**Regularization parameters**

One critical parameter of vortex methods is the regularization parameter, also referred to as core radius. We currently recommend to set the regularization parameter as a fraction of the spanwise discretization (:math:`\Delta r`), that is: **RegDeterMethod=3** , **WakeRegFactor=0.6**, **WingRegFactor=0.6**.
When the RegFactors are set as function of the spanwise discretization, we expect the factors to be somewhere between 0.25 and 3. 


We also recommend to have the regularization increasing with downstream distance:
**WakeRegMethod=3**.  
The factor with which the regularization parameter will increase with downstream distance can be set as **CoreSpreadEddyVisc=1000** for modern multi-MW turbines. 
When plotting the wake, (**WrVTK**), if the wake appears to be "too distorted" for a steady state simulation, increase the **CoreSpreadEddyVisc** parameter to "smoothen" the wake.





**Python script**

The following python script computes the parameters according to these guidelines.
(Check here for updates: `olaf.py <https://github.com/ebranlard/welib/blob/dev/welib/fast/olaf.py>`_)


.. code::

   def OLAFParams(omega_rpm, U0, R, a=0.3, aScale=1.2,
             deltaPsiDeg=6, nPerRot=None,
             targetFreeWakeLengthD=1, 
             targetWakeLengthD=4., 
             nNWrot=8, nFWrot=0, nFWrotFree=0,
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

       # All Wake (AW) panels - Wake length from mean wind speed
       targetWakeLength = targetWakeLengthD * 2 * R
       nAWPanels_FromU0 = int(targetWakeLength / (Uc*dt_fvw))
       # Free near wake panels (based on distance)
       targetFreeWakeLength = targetFreeWakeLengthD * 2 * R
       nNWPanelsFree_FromU0 = int(targetFreeWakeLength / (Uc*dt_fvw))
       # Far wake (FW) panels, always from number of rotations
       nFWPanels     = int(nFWrot*nPerRot)
       nFWPanelsFree = int(nFWrotFree*nPerRot)
       # Wake length from rotational speed and number of rotations
       nAWPanels_FromRot = int(nNWrot*nPerRot) # Total number of panels NW+FW

       # Below we chose between criteria on number of rotation or donwstream distance
       # This can be adapted/improved
       myprint('Number of panels (NW free) from wind speed and distance:{:15d}'.format(nNWPanelsFree_FromU0))
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
       nNWPanels = nAWPanels - nFWPanels # nNW = All-Far Wake

       # See "free" near wake
       nNWPanelsFree = nNWPanelsFree_FromU0
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
           myprint('{:15d} nNWPanelsFree ({:5.1f} rot, {:5.1f}D)'.format(nNWPanelsFree, n2R(nNWPanelsFree), n2L(nNWPanelsFree)))
           myprint('{:15d} nFWPanels     ({:5.1f} rot, {:5.1f}D)'.format(nFWPanels    , n2R(nFWPanels    ), n2L(nFWPanels    )))
           myprint('{:15d} nFWPanelsFree ({:5.1f} rot, {:5.1f}D)'.format(nFWPanelsFree, n2R(nFWPanelsFree), n2L(nFWPanelsFree)))

       return dt_fvw, tMin, nNWPanels, nNWPanelsFree, nFWPanels, nFWPanelsFree


