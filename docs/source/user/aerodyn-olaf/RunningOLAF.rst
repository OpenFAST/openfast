
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


Guidelines
~~~~~~~~~~

Most options of OLAF can be left to default. The results will depend on the time discretization, wake length, and regularization parameters. We provide guidelines for these parameters in this section, together with a simple python code to compute these parameters.
Please check this section again as we might further refine our guidelines with time.


**Time step and wake length**
We recommend to set OLAF's time step (**DTfvw**) such that it corresponds to :math:`\Delta \psi = 6` degrees of one rotor revolution:

.. math::
   
    \Delta t
    = \frac{\Delta \psi_\text{rad}}{\Omega_\text{rad/s}}
    = \frac{\Delta \psi_\text{deg}}{6 \times \Omega_\text{RPM}}

If the structural solver requires a smaller time step, the time step for the glue code can be set to a different value than **DTfvw** as long as **DTfvw** is a multiple of the glue code time step.


We recommend to set the near wake length to the number of time steps needed to reach two rotor revolutions. For the far wake, we recommend 10 rotor revolutions. 
For the free far-wake, we recommend to set the distance to a value somewhere between 25% and 50% of the far wake length, (e.g. 3 revolutions).

The following python script computes the parameters according to these guidelines.

.. code::

   def OLAFParams(omega_rpm, deltaPsiDeg=6, nNWrot=2, nFWrot=10, nFWrotFree=3, nPerAzimuth=None):
       """ 
       Computes recommended time step and wake length based on the rotational speed in RPM

       INPUTS:
        - omega_rpm: rotational speed in RPM
        - deltaPsiDeg : azimuthal discretization in deg
        - nNWrot : number of near wake rotations
        - nFWrot : total number of far wake rotations
        - nFWrotFree : number of far wake rotations that are free

           deltaPsiDeg  -  nPerAzimuth
                5            72    
                6            60    
                7            51.5  
                8            45    
       """
       omega = omega_rpm*2*np.pi/60
       T = 2*np.pi/omega
       if nPerAzimuth is not None:
           dt_wanted    = np.around(T/nPerAzimuth,3)
       else:
           dt_wanted    = np.around(deltaPsiDeg/(6*omega_rpm),3)
           nPerAzimuth = int(2*np.pi /(deltaPsiDeg*np.pi/180))

       nNWPanel     = nNWrot*nPerAzimuth
       nFWPanel     = nFWrot*nPerAzimuth
       nFWPanelFree = nFWrotFree*nPerAzimuth

       print(dt_wanted              , '  DTfvw')
       print(int      (nNWPanel    ), '  nNWPanel      ')
       print(int      (nFWPanel    ), '  WakeLength    ')
       print(int      (nFWPanelFree), '  FreeWakeLength')

       return dt_wanted, nNWPanel, nFWPanel, nFWPanelFree


**Regularization parameters**

One critical parameter of vortex methods is the regularization parameter, also referred to as core radius. We currently recommend to set the regularization parameter as a fraction of the spanwise discretization, that is: **RegDetMethod=3** , **WakeRegFactor=0.6**, **WingRegFactor=0.6**.
We will likely update these guidelines in the future.


We also recommend to have the regularization increasing with downstream distance:
**WakeRegMethod=3**. 

The factor with which the regularization parameter will increase with downstream distance can be set as
**CoreSpreadEddyVisc=1000** for modern multi-MW turbines. Further guidelines will follow for this parameter in the future. 




