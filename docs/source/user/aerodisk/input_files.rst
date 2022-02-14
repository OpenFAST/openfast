.. _adsk_input-files:

Input and Output Files
======================


Units
-----

AeroDisk uses the SI system (kg, m, s, N). 

.. _adsk_input-file:

Input file
----------


.. code::
    
   --- AERO DISK INPUT FILE -------
   Sample actuator disk input file
   --- ENVIRONMENTAL CONDITIONS ---
         1.225   AirDens         - Air density (kg/m^3) (or "default")
   --- ACTUATOR DISK PROPERTIES ---
        63.0     RotRad          - Rotor radius (m) (or "default")
   "RtSpd,VRel,Skew,Pitch" InColNames      - Input column headers (string) {may include a combination of "TSR, RtSpd, VRel, Pitch, Skew"} (up to 4 columns) [choose TSR or RtSpd,VRel; if Skew is absent, Skew is modeled as (COS(Skew))^2]
   9,1,1,1       InColDims       - Number of values in each column (-) (must have same number of columns as InColName) [each >=2]
   RtSpd     VRel     Skew    Pitch   C_Fx   C_Fy   C_Fz   C_Mx   C_My   C_Mz
   (rpm)     (m/s)    (deg)   (deg)   (-)    (-)    (-)    (-)    (-)    (-)
   -20.0      0.0       0.0    0.0
   -15.0      0.0       0.0    0.0
   -10.0      0.0       0.0    0.0
    -5.0      0.0       0.0    0.0
     0.0      0.0       0.0    0.0
     5.0      0.0       0.0    0.0
    10.0      0.0       0.0    0.0
    15.0      0.0       0.0    0.0
    20.0      0.0       0.0    0.0
   --- OUTPUTS --------------------
                 OutList         - The next line(s) contains a list of output parameters.  See OutListParameters.xlsx for a listing of available output channels, (-)
   END of input file (the word "END" must appear in the first 3 columns of this last OutList line)
   --------------------------------






.. _adsk_outputs:

Outputs
-------

The write outputs are:
 -  "ADSpeed":    Actuator disk rotational speed               (rpm)   
 -  "ADTSR":      Actuator disk tip-speed ratio                (-)   
 -  "ADPitch":    Actuator-disk collective blade-pitch angle   (deg)   
 -  "ADVWindx, ADVWindy, ADVWindz":    Actuator-disk-averaged wind velocity in the local coordinate system  (m/s)
 -  "ADSTVx, ADSTVy, ADSTVz":          Actuator-disk structural translational velocity in the local coordinate system  (m/s)
 -  "ADVRel":     Actuator-disk-averaged relative wind speed    (m/s)
 -  "ADSkew":     Actuator-disk inflow-skew angle               (deg)
 -  "ADCp, ADCt, ADCq":   Actuator-disk power, thrust, and torque coefficients  (-)
 -  "ADFx, ADFy, ADFz":   Actuator disk aerodynamic force loads in the local coordinate system (N)
 -  "ADMx, ADMy, ADMz":   Actuator disk aerodynamic moment loads in the local coordinate system (N-m)
 -  "ADPower":   Actuator disk power   (W)

