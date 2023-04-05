module ifw_test_tools

    use NWTC_IO
    use NWTC_Library
    use NWTC_Library_Types

    implicit none

contains

    function getInputFileData()

        INTEGER                             :: ErrStat
        CHARACTER(ErrMsgLen)                :: ErrMsg
        TYPE(FileInfoType)                  :: getInputFileData
        CHARACTER(1024), DIMENSION(69)      :: data = (/ &
            '------- InflowWind v3.01.* INPUT FILE -------------------------------------------------------------------------                                                                    ', &
            'Steady 8 m/s winds with no shear for FAST CertTests #20 and #25                                                                                                                    ', &
            '---------------------------------------------------------------------------------------------------------------                                                                    ', &
            '       false  Echo           - Echo input data to <RootName>.ech (flag)                                                                                                            ', &
            '          1   WindType       - switch for wind file type (1=steady; 2=uniform; 3=binary TurbSim FF; 4=binary Bladed-style FF; 5=HAWC format; 6=User defined; 7=native Bladed FF)   ', &
            '          0   PropagationDir - Direction of wind propagation (meteoroligical rotation from aligned with X (positive rotates towards -Y) -- degrees)                                ', &
            '          0   VFlowAng       - Upflow angle (degrees) (not used for native Bladed format WindType=7)                                                                               ', &
            '      false   VelInterpCubic - Use cubic interpolation for velocity in time (false=linear, true=cubic) [Used with WindType=2,3,4,5,7]                                              ', &
            '          1   NWindVel       - Number of points to output the wind velocity    (0 to 9)                                                                                            ', &
            '          0   WindVxiList    - List of coordinates in the inertial X direction (m)                                                                                                 ', &
            '          0   WindVyiList    - List of coordinates in the inertial Y direction (m)                                                                                                 ', &
            '         90   WindVziList    - List of coordinates in the inertial Z direction (m)                                                                                                 ', &
            '================== Parameters for Steady Wind Conditions [used only for WindType = 1] =========================                                                                    ', &
            '        8.0   HWindSpeed     - Horizontal windspeed                            (m/s)                                                                                               ', &
            '         90   RefHt          - Reference height for horizontal wind speed      (m)                                                                                                 ', &
            '        0.1   PLexp          - Power law exponent                              (-)                                                                                                 ', &
            '================== Parameters for Uniform wind file   [used only for WindType = 2] ============================                                                                    ', &
            '"Wind/08ms.wnd"    FileName_Uni   - Filename of time series data for uniform wind field.      (-)                                                                                  ', &
            '         90   RefHt_Uni      - Reference height for horizontal wind speed                (m)                                                                                       ', &
            '     125.88   RefLength      - Reference length for linear horizontal and vertical sheer (-)                                                                                       ', &
            '================== Parameters for Binary TurbSim Full-Field files   [used only for WindType = 3] ==============                                                                    ', &
            '"Wind/08ms.wnd"    filename_bts       - name of the full field wind file to use (.bts)                                                                                             ', &
            '================== Parameters for Binary Bladed-style Full-Field files   [used only for WindType = 4] =========                                                                    ', &
            '"unused"      FilenameRoot   - Rootname of the full-field wind file to use (.wnd, .sum)                                                                                            ', &
            'False         TowerFile      - Have tower file (.twr) (flag)                                                                                                                       ', &
            '================== Parameters for HAWC-format binary files  [Only used with WindType = 5] =====================                                                                    ', &
            '"wasp\Output\basic_5u.bin"    FileName_u     - name of the file containing the u-component fluctuating wind (.bin)                                                                 ', &
            '"wasp\Output\basic_5v.bin"    FileName_v     - name of the file containing the v-component fluctuating wind (.bin)                                                                 ', &
            '"wasp\Output\basic_5w.bin"    FileName_w     - name of the file containing the w-component fluctuating wind (.bin)                                                                 ', &
            '         64   nx             - number of grids in the x direction (in the 3 files above) (-)                                                                                       ', &
            '         32   ny             - number of grids in the y direction (in the 3 files above) (-)                                                                                       ', &
            '         32   nz             - number of grids in the z direction (in the 3 files above) (-)                                                                                       ', &
            '         16   dx             - distance (in meters) between points in the x direction    (m)                                                                                       ', &
            '          3   dy             - distance (in meters) between points in the y direction    (m)                                                                                       ', &
            '          3   dz             - distance (in meters) between points in the z direction    (m)                                                                                       ', &
            '         90   RefHt_HAWC     - reference height; the height (in meters) of the vertical center of the grid (m)                                                                     ', &
            ' -------------   Scaling parameters for turbulence   ---------------------------------------------------------                                                                     ', &
            '          1   ScaleMethod    - Turbulence scaling method   [0 = none, 1 = direct scaling, 2 = calculate scaling factor based on a desired standard deviation]                      ', &
            '          1   SFx            - Turbulence scaling factor for the x direction (-)   [ScaleMethod=1]                                                                                 ', &
            '          1   SFy            - Turbulence scaling factor for the y direction (-)   [ScaleMethod=1]                                                                                 ', &
            '          1   SFz            - Turbulence scaling factor for the z direction (-)   [ScaleMethod=1]                                                                                 ', &
            '         12   SigmaFx        - Turbulence standard deviation to calculate scaling from in x direction (m/s)    [ScaleMethod=2]                                                     ', &
            '          8   SigmaFy        - Turbulence standard deviation to calculate scaling from in y direction (m/s)    [ScaleMethod=2]                                                     ', &
            '          2   SigmaFz        - Turbulence standard deviation to calculate scaling from in z direction (m/s)    [ScaleMethod=2]                                                     ', &
            '  -------------   Mean wind profile parameters (added to HAWC-format files)   ---------------------------------                                                                    ', &
            '          5   URef           - Mean u-component wind speed at the reference height (m/s)                                                                                           ', &
            '          2   WindProfile    - Wind profile type (0=constant;1=logarithmic,2=power law)                                                                                            ', &
            '          0   PLExp_HAWC     - Power law exponent (-) (used for PL wind profile type only)                                                                                         ', &
            '       0.03   Z0             - Surface roughness length (m) (used for LG wind profile type only)                                                                                   ', &
            '          0   XOffset        - Initial offset in +x direction (shift of wind box)                                                                                                  ', &
            '  ---------------- LIDAR Parameters ---------------------------------------------------------------------------                                                                    ', &
            '          0 SensorType          - Switch for lidar configuration (0 = None, 1 = Single Point Beam(s), 2 = Continuous, 3 = Pulsed)                                                  ', &
            '          0 NumPulseGate        - Number of lidar measurement gates (used when SensorType = 3)                                                                                     ', &
            '         30 PulseSpacing        - Distance between range gates (m) (used when SensorType = 3)                                                                                      ', &
            '          0 NumBeam             - Number of lidar measurement beams (0-5)(used when SensorType = 1)                                                                                ', &
            '       -200 FocalDistanceX      - Focal distance co-ordinates of the lidar beam in the x direction (relative to hub height) (only first coordinate used for SensorType 2 and 3) (m)', &
            '          0 FocalDistanceY      - Focal distance co-ordinates of the lidar beam in the y direction (relative to hub height) (only first coordinate used for SensorType 2 and 3) (m)', &
            '          0 FocalDistanceZ      - Focal distance co-ordinates of the lidar beam in the z direction (relative to hub height) (only first coordinate used for SensorType 2 and 3) (m)', &
            '0.0 0.0 0.0 RotorApexOffsetPos  - Offset of the lidar from hub height (m)                                                                                                          ', &
            '         17 URefLid             - Reference average wind speed for the lidar[m/s]                                                                                                  ', &
            '       0.25 MeasurementInterval - Time between each measurement [s]                                                                                                                ', &
            '      False LidRadialVel        - TRUE => return radial component, FALSE => return "x" direction estimate                                                                          ', &
            '          1 ConsiderHubMotion   - Flag whether to consider the hub motion impact on Lidar measurements                                                                             ', &
            '====================== OUTPUT ==================================================                                                                                                   ', &
            'False         SumPrint     - Print summary data to <RootName>.IfW.sum (flag)                                                                                                       ', &
            '              OutList      - The next line(s) contains a list of output parameters.  See OutListParameters.xlsx for a listing of available output channels, (-)                    ', &
            '"Wind1VelX,Wind1VelY,Wind1VelZ"     - Wind velocity at point WindVxiList(1),WindVyiList(1),WindVziList(1).  X, Y, and Z direction components.                                      ', &
            'END of input file (the word "END" must appear in the first 3 columns of this last OutList line)                                                                                    ', &
            '---------------------------------------------------------------------------------------                                                                                            ' &
        /)

        CALL InitFileInfo(data, getInputFileData, ErrStat, ErrMsg)

    end function

    function getInputFileDataWindType2()

        INTEGER                             :: ErrStat
        CHARACTER(ErrMsgLen)                :: ErrMsg
        TYPE(FileInfoType)                  :: getInputFileDataWindType2
        CHARACTER(1024), DIMENSION(69)      :: data = (/ &
            '------- InflowWind v3.01.* INPUT FILE -------------------------------------------------------------------------                                                                    ', &
            'Steady 8 m/s winds with no shear for FAST CertTests #20 and #25                                                                                                                    ', &
            '---------------------------------------------------------------------------------------------------------------                                                                    ', &
            '       true   Echo           - Echo input data to <RootName>.ech (flag)                                                                                                            ', &
            '          2   WindType       - switch for wind file type (1=steady; 2=uniform; 3=binary TurbSim FF; 4=binary Bladed-style FF; 5=HAWC format; 6=User defined; 7=native Bladed FF)   ', &
            '          0   PropagationDir - Direction of wind propagation (meteoroligical rotation from aligned with X (positive rotates towards -Y) -- degrees)                                ', &
            '          0   VFlowAng       - Upflow angle (degrees) (not used for native Bladed format WindType=7)                                                                               ', &
            '      false   VelInterpCubic - Use cubic interpolation for velocity in time (false=linear, true=cubic) [Used with WindType=2,3,4,5,7]                                              ', &
            '          1   NWindVel       - Number of points to output the wind velocity    (0 to 9)                                                                                            ', &
            '          0   WindVxiList    - List of coordinates in the inertial X direction (m)                                                                                                 ', &
            '          0   WindVyiList    - List of coordinates in the inertial Y direction (m)                                                                                                 ', &
            '         90   WindVziList    - List of coordinates in the inertial Z direction (m)                                                                                                 ', &
            '================== Parameters for Steady Wind Conditions [used only for WindType = 1] =========================                                                                    ', &
            '        8.0   HWindSpeed     - Horizontal windspeed                            (m/s)                                                                                               ', &
            '         90   RefHt          - Reference height for horizontal wind speed      (m)                                                                                                 ', &
            '        0.1   PLexp          - Power law exponent                              (-)                                                                                                 ', &
            '================== Parameters for Uniform wind file   [used only for WindType = 2] ============================                                                                    ', &
            '"Wind/08ms.wnd"    FileName_Uni   - Filename of time series data for uniform wind field.      (-)                                                                                  ', &
            '         90   RefHt_Uni      - Reference height for horizontal wind speed                (m)                                                                                       ', &
            '     125.88   RefLength      - Reference length for linear horizontal and vertical sheer (-)                                                                                       ', &
            '================== Parameters for Binary TurbSim Full-Field files   [used only for WindType = 3] ==============                                                                    ', &
            '"Wind/08ms.wnd"    FileName_BTS          - Filename of time series data for uniform wind field.      (-)                                                                           ', &
            '================== Parameters for Binary Bladed-style Full-Field files   [used only for WindType = 4] =========                                                                    ', &
            '"unused"      FilenameRoot   - Rootname of the full-field wind file to use (.wnd, .sum)                                                                                            ', &
            'False         TowerFile      - Have tower file (.twr) (flag)                                                                                                                       ', &
            '================== Parameters for HAWC-format binary files  [Only used with WindType = 5] =====================                                                                    ', &
            '"wasp\Output\basic_5u.bin"    FileName_u     - name of the file containing the u-component fluctuating wind (.bin)                                                                 ', &
            '"wasp\Output\basic_5v.bin"    FileName_v     - name of the file containing the v-component fluctuating wind (.bin)                                                                 ', &
            '"wasp\Output\basic_5w.bin"    FileName_w     - name of the file containing the w-component fluctuating wind (.bin)                                                                 ', &
            '         64   nx             - number of grids in the x direction (in the 3 files above) (-)                                                                                       ', &
            '         32   ny             - number of grids in the y direction (in the 3 files above) (-)                                                                                       ', &
            '         32   nz             - number of grids in the z direction (in the 3 files above) (-)                                                                                       ', &
            '         16   dx             - distance (in meters) between points in the x direction    (m)                                                                                       ', &
            '          3   dy             - distance (in meters) between points in the y direction    (m)                                                                                       ', &
            '          3   dz             - distance (in meters) between points in the z direction    (m)                                                                                       ', &
            '         90   RefHt_HAWC     - reference height; the height (in meters) of the vertical center of the grid (m)                                                                     ', &
            ' -------------   Scaling parameters for turbulence   ---------------------------------------------------------                                                                     ', &
            '          1   ScaleMethod    - Turbulence scaling method   [0 = none, 1 = direct scaling, 2 = calculate scaling factor based on a desired standard deviation]                      ', &
            '          1   SFx            - Turbulence scaling factor for the x direction (-)   [ScaleMethod=1]                                                                                 ', &
            '          1   SFy            - Turbulence scaling factor for the y direction (-)   [ScaleMethod=1]                                                                                 ', &
            '          1   SFz            - Turbulence scaling factor for the z direction (-)   [ScaleMethod=1]                                                                                 ', &
            '         12   SigmaFx        - Turbulence standard deviation to calculate scaling from in x direction (m/s)    [ScaleMethod=2]                                                     ', &
            '          8   SigmaFy        - Turbulence standard deviation to calculate scaling from in y direction (m/s)    [ScaleMethod=2]                                                     ', &
            '          2   SigmaFz        - Turbulence standard deviation to calculate scaling from in z direction (m/s)    [ScaleMethod=2]                                                     ', &
            '  -------------   Mean wind profile parameters (added to HAWC-format files)   ---------------------------------                                                                    ', &
            '          5   URef           - Mean u-component wind speed at the reference height (m/s)                                                                                           ', &
            '          2   WindProfile    - Wind profile type (0=constant;1=logarithmic,2=power law)                                                                                            ', &
            '          0   PLExp_HAWC     - Power law exponent (-) (used for PL wind profile type only)                                                                                         ', &
            '       0.03   Z0             - Surface roughness length (m) (used for LG wind profile type only)                                                                                   ', &
            '          0   XOffset        - Initial offset in +x direction (shift of wind box)                                                                                                  ', &
            '  ---------------- LIDAR Parameters ---------------------------------------------------------------------------                                                                    ', &
            '          0 SensorType          - Switch for lidar configuration (0 = None, 1 = Single Point Beam(s), 2 = Continuous, 3 = Pulsed)                                                  ', &
            '          0 NumPulseGate        - Number of lidar measurement gates (used when SensorType = 3)                                                                                     ', &
            '         30 PulseSpacing        - Distance between range gates (m) (used when SensorType = 3)                                                                                      ', &
            '          0 NumBeam             - Number of lidar measurement beams (0-5)(used when SensorType = 1)                                                                                ', &
            '       -200 FocalDistanceX      - Focal distance co-ordinates of the lidar beam in the x direction (relative to hub height) (only first coordinate used for SensorType 2 and 3) (m)', &
            '          0 FocalDistanceY      - Focal distance co-ordinates of the lidar beam in the y direction (relative to hub height) (only first coordinate used for SensorType 2 and 3) (m)', &
            '          0 FocalDistanceZ      - Focal distance co-ordinates of the lidar beam in the z direction (relative to hub height) (only first coordinate used for SensorType 2 and 3) (m)', &
            '0.0 0.0 0.0 RotorApexOffsetPos  - Offset of the lidar from hub height (m)                                                                                                          ', &
            '         17 URefLid             - Reference average wind speed for the lidar[m/s]                                                                                                  ', &
            '       0.25 MeasurementInterval - Time between each measurement [s]                                                                                                                ', &
            '      False LidRadialVel        - TRUE => return radial component, FALSE => return "x" direction estimate                                                                          ', &
            '          1 ConsiderHubMotion   - Flag whether to consider the hub motion impact on Lidar measurements                                                                             ', &
            '====================== OUTPUT ==================================================                                                                                                   ', &
            'False         SumPrint     - Print summary data to <RootName>.IfW.sum (flag)                                                                                                       ', &
            '              OutList      - The next line(s) contains a list of output parameters.  See OutListParameters.xlsx for a listing of available output channels, (-)                    ', &
            '"Wind1VelX,Wind1VelY,Wind1VelZ"     - Wind velocity at point WindVxiList(1),WindVyiList(1),WindVziList(1).  X, Y, and Z direction components.                                      ', &
            'END of input file (the word "END" must appear in the first 3 columns of this last OutList line)                                                                                    ', &
            '---------------------------------------------------------------------------------------                                                                                            ' &
        /)

        CALL InitFileInfo(data, getInputFileDataWindType2, ErrStat, ErrMsg)

    end function

end module
