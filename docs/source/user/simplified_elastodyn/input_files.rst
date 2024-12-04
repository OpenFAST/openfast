.. _sed_input-files:

Input and Output Files
======================


Units
-----

SED uses the SI system (kg, m, s, N). 

.. _sed_input-file:

Input file
----------


.. code::
    
    ------- SIMPLIFIED ELASTODYN INPUT FILE ----------------------------------------
   Comment
    ---------------------- SIMULATION CONTROL --------------------------------------
    False         Echo        - Echo input data to "<RootName>.ech" (flag)
              3   Method      - Integration method: {1: RK4, 2: AB4, or 3: ABM4} (-)
    "default"     DT          - Integration time step (s)
    ---------------------- DEGREES OF FREEDOM --------------------------------------
    True          GenDOF      - Generator DOF (flag)
    ---------------------- INITIAL CONDITIONS --------------------------------------
              0   Azimuth     - Initial azimuth angle for blades (degrees)
              0   BlPitch     - Blades initial pitch (degrees)
            0.0   RotSpeed    - Initial or fixed rotor speed (rpm)
              0   NacYaw      - Initial or fixed nacelle-yaw angle (degrees)
              0   PtfmPitch   - Fixed pitch tilt rotational displacement of platform (degrees)
    ---------------------- TURBINE CONFIGURATION -----------------------------------
              3   NumBl       - Number of blades (-)
             63   TipRad      - The distance from the rotor apex to the blade tip (meters)
            1.5   HubRad      - The distance from the rotor apex to the blade root (meters)
           -2.5   PreCone     - Blades cone angle (degrees)
        -5.0191   OverHang    - Distance from yaw axis to rotor apex [3 blades] or teeter pin [2 blades] (meters)
             -5   ShftTilt    - Rotor shaft tilt angle (degrees)
        1.96256   Twr2Shft    - Vertical distance from the tower-top to the rotor shaft (meters)
           87.6   TowerHt     - Height of tower above ground level [onshore] or MSL [offshore] (meters)
    ---------------------- MASS AND INERTIA ----------------------------------------
         115926   RotIner     - Rot inertia about rotor axis [blades + hub] (kg m^2)
        534.116   GenIner     - Generator inertia about HSS (kg m^2)
    ---------------------- DRIVETRAIN ----------------------------------------------
            100   GBoxEff     - Gearbox efficiency (%)
             97   GBRatio     - Gearbox ratio (-)
    ---------------------- OUTPUT --------------------------------------------------
                  OutList     - The next line(s) contains a list of output parameters.  See OutListParameters.xlsx for a listing of available output channels, (-)
    "Azimuth"                 - Blades azimuth angle
    "RotSpeed"                - Low-speed shaft rotational speed
    "RotAcc"                  - Low-speed shaft rotational acceleration
    END of input file (the word "END" must appear in the first 3 columns of this last OutList line)
    -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------





.. _sed_outputs:

Outputs
-------

The write outputs are:
 -  "Azimuth" : Blades azimuth angle (deg) between 0 and 2pi
 -  "RotSpeed": Low-speed shaft rotational speed (rpm)
 -  "RotAcc": Low-speed shaft rotational acceleration (rad/s^2)
 -  "GenSpeed": High-speed shaft rotational speed (rpm)
 -  "GenAcc": High-speed shaft rotational acceleration (rad/s^2)
