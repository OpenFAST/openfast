def default_StC_vt():
    StC = {}

    StC['Echo'] = True

    #---------------------- StC DEGREES OF FREEDOM ----------------------------------
    StC['StC_DOF_MODE'] = 0
    StC['StC_X_DOF']    = False
    StC['StC_Y_DOF']    = False
    StC['StC_Z_DOF']    = False

    # ---------------------- StC LOCATION ------------------------------------------- [relative to the reference origin of component attached to]
    StC['StC_P_X']  = 0     #      - At rest X position of StC (m)
    StC['StC_P_Y']  = 0     #      - At rest Y position of StC (m)
    StC['StC_P_Z']  = 0     #      - At rest Z position of StC (m)
    # ---------------------- StC INITIAL CONDITIONS --------------------------------- [used only when StC_DOF_MODE=1 or 2]
    StC['StC_X_DSP'] = 0    #    - StC X initial displacement (m) [relative to at rest position]
    StC['StC_Y_DSP'] = 0    #    - StC Y initial displacement (m) [relative to at rest position]
    StC['StC_Z_DSP'] = 0    #    - StC Z initial displacement (m) [relative to at rest position; used only when StC_DOF_MODE=1 and StC_Z_DOF=TRUE]
    StC['StC_Z_PreLd'] = "none"    #    - StC Z initial displacement (m) [relative to at rest position; used only when StC_DOF_MODE=1 and StC_Z_DOF=TRUE]
    # ---------------------- StC CONFIGURATION -------------------------------------- [used only when StC_DOF_MODE=1 or 2]
    StC['StC_X_PSP'] = 0    #    - Positive stop position (maximum X mass displacement) (m)
    StC['StC_X_NSP'] = 0    #    - Negative stop position (minimum X mass displacement) (m)
    StC['StC_Y_PSP'] = 0    #    - Positive stop position (maximum Y mass displacement) (m)
    StC['StC_Y_NSP'] = 0    #    - Negative stop position (minimum Y mass displacement) (m)
    StC['StC_Z_PSP'] = 0    #    - Positive stop position (maximum Z mass displacement) (m) [used only when StC_DOF_MODE=1 and StC_Z_DOF=TRUE]
    StC['StC_Z_NSP'] = 0    #    - Negative stop position (minimum Z mass displacement) (m) [used only when StC_DOF_MODE=1 and StC_Z_DOF=TRUE]
    #---------------------- StC MASS, STIFFNESS, & DAMPING ------------------------- [used only when StC_DOF_MODE=1 or 2]
    StC['StC_X_M']      = 0     #  - StC X mass (kg) [must equal StC_Y_M for StC_DOF_MODE = 2]
    StC['StC_Y_M']      = 0     #  - StC Y mass (kg) [must equal StC_X_M for StC_DOF_MODE = 2]
    StC['StC_Z_M']      = 0     #  - StC Z mass (kg) [used only when StC_DOF_MODE=1 and StC_Z_DOF=TRUE]
    StC['StC_XY_M']     = 0     #  - StC Z mass (kg) [used only when StC_DOF_MODE=2]
    StC['StC_X_K']      = 0     #  - StC X stiffness (N/m)
    StC['StC_Y_K']      = 0     #  - StC Y stiffness (N/m)
    StC['StC_Z_K']      = 0     #  - StC Z stiffness (N/m) [used only when StC_DOF_MODE=1 and StC_Z_DOF=TRUE]
    StC['StC_X_C']      = 0     #  - StC X damping (N/(m/s))
    StC['StC_Y_C']      = 0     #  - StC Y damping (N/(m/s))
    StC['StC_Z_C']      = 0     #  - StC Z damping (N/(m/s)) [used only when StC_DOF_MODE=1 and StC_Z_DOF=TRUE]
    StC['StC_X_KS']     = 0     #  - Stop spring X stiffness (N/m)
    StC['StC_Y_KS']     = 0     #  - Stop spring Y stiffness (N/m)
    StC['StC_Z_KS']     = 0     #  - Stop spring Z stiffness (N/m) [used only when StC_DOF_MODE=1 and StC_Z_DOF=TRUE]
    StC['StC_X_CS']     = 0     #  - Stop spring X damping (N/(m/s))
    StC['StC_Y_CS']     = 0     #  - Stop spring Y damping (N/(m/s))
    StC['StC_Z_CS']     = 0     #  - Stop spring Z damping (N/(m/s)) [used only when StC_DOF_MODE=1 and StC_Z_DOF=TRUE]
    #---------------------- StC USER-DEFINED SPRING FORCES ------------------------- [used only when StC_DOF_MODE=1 or 2]
    StC['Use_F_TBL']    = False   # - Use spring force from user-defined table (flag)
    StC['NKInpSt']      = 2       # - Number of spring force input stations
    #---------------------- StC SPRING FORCES TABLE -------------------------------- [used only when StC_DOF_MODE=1 or 2]
    StC['SpringForceTable'] = {}
    StC['SpringForceTable']['X']        =  [-6,6]
    StC['SpringForceTable']['F_X']      =  [-4.8e6,4.8e6]
    StC['SpringForceTable']['Y']        =  [-6,6]
    StC['SpringForceTable']['F_Y']      =  [-4.8e6,4.8e6]
    StC['SpringForceTable']['Z']        =  [-6,6]
    StC['SpringForceTable']['F_Z']      =  [-4.8e6,4.8e6]
    # ---------------------- StructCtrl CONTROL -------------------------------------------- [used only when StC_DOF_MODE=1 or 2]
    StC['StC_CMODE']        = 0       # - Control mode (switch) {0:none; 1: Semi-Active Control Mode; 2: Active Control Mode}
    StC['StC_SA_MODE']      = 1       # - Semi-Active control mode {1: velocity-based ground hook control; 2: Inverse velocity-based ground hook control; 3: displacement-based ground hook control 4: Phase difference Algorithm with Friction Force 5: Phase difference Algorithm with Damping Force} (-)
    StC['StC_CChan']        = 0       # - Control mode (switch) {0:none; 1: Semi-Active Control Mode; 2: Active Control Mode}
    StC['StC_X_C_HIGH']     = 0       # - StC X high damping for ground hook control
    StC['StC_X_C_LOW']      = 0       # - StC X low damping for ground hook control
    StC['StC_Y_C_HIGH']     = 0       # - StC Y high damping for ground hook control
    StC['StC_Y_C_LOW']      = 0       # - StC Y low damping for ground hook control
    StC['StC_Z_C_HIGH']     = 0       # - StC Z high damping for ground hook control [used only when StC_DOF_MODE=1 and StC_Z_DOF=TRUE]
    StC['StC_Z_C_LOW']      = 0       # - StC Z low damping for ground hook control  [used only when StC_DOF_MODE=1 and StC_Z_DOF=TRUE]
    StC['StC_X_C_BRAKE']    = 0       # - StC X high damping for braking the StC (Don't use it now. should be zero)
    StC['StC_Y_C_BRAKE']    = 0       # - StC Y high damping for braking the StC (Don't use it now. should be zero)
    StC['StC_Z_C_BRAKE']    = 0       # - StC Z high damping for braking the StC (Don't use it now. should be zero) [used only when StC_DOF_MODE=1 and StC_Z_DOF=TRUE]
    # ---------------------- TLCD --------------------------------------------------- [used only when StC_DOF_MODE=3]
    StC['L_X']              = 0     # - X TLCD total length (m)
    StC['B_X']              = 0     # - X TLCD horizontal length (m)
    StC['area_X']           = 0     # - X TLCD cross-sectional area of vertical column (m^2)
    StC['area_ratio_X']     = 0     # - X TLCD cross-sectional area ratio (vertical column area divided by horizontal column area) (-)
    StC['headLossCoeff_X']  = 0     # - X TLCD head loss coeff (-)
    StC['rho_X']            = 0     # - X TLCD liquid density (kg/m^3)
    StC['L_Y']              = 0     # - Y TLCD total length (m)
    StC['B_Y']              = 0     # - Y TLCD horizontal length (m)
    StC['area_Y']           = 0     # - Y TLCD cross-sectional area of vertical column (m^2)
    StC['area_ratio_Y']     = 0     # - Y TLCD cross-sectional area ratio (vertical column area divided by horizontal column area) (-)
    StC['headLossCoeff_Y']  = 0     # - Y TLCD head loss coeff (-)
    StC['rho_Y']            = 0     # - Y TLCD liquid density (kg/m^3)
    # ---------------------- PRESCRIBED TIME SERIES --------------------------------- [used only when StC_DOF_MODE=4]
    StC['PrescribedForcesCoord']    = 0                 # - Prescribed forces are in global or local coordinates (switch) {1: global; 2: local}
    StC['PrescribedForcesFile']     = "unused"          #  - Time series force and moment (7 columns of time, FX, FY, FZ, MX, MY, MZ)

    return StC


if __name__=="__main__":
    
    a = default_StC_vt()
    print('here')

