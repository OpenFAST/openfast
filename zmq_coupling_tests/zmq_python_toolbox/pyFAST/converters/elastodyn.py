import numpy as np
import pandas as pd


def elastodynBlade2Hawc2_raw(dfBld, R, E=3e10):
    """ 
    INPUTS:
     - dfBld: dataframe resulting from reading an ElastoDyn blade file
              with columns: [ 'BlFract_[-]', 'PitchAxis_[-]', 'StrcTwst_[deg]', 'BMassDen_[kg/m]', 'FlpStff_[Nm^2]', 'EdgStff_[Nm^2]'
      - R : rotor radius [m]
    OUTPUTS:
      - dfMeanLine: hawc2 c2def as DataFrame with columns ['x_[m]','y_[m]','z_[m]','twist_[deg]']
      - dfStructure: hawc2 stdata as Dataframe
    """
    r      = dfBld['BlFract_[-]'].values*R
    twist  = dfBld['StrcTwst_[deg]'].values
    m      = dfBld['BMassDen_[kg/m]'].values
    EIflap = dfBld['FlpStff_[Nm^2]'].values
    EIedge = dfBld['EdgStff_[Nm^2]'].values


    # --- ElastoDyn 2 Hawc2 Structural data
    columns=['r_[m]','m_[kg/m]','x_cg_[m]','y_cg_[m]','ri_x_[m]','ri_y_[m]','x_sh_[m]','y_sh_[m]','E_[N/m^2]','G_[N/m^2]','I_x_[m^4]','I_y_[m^4]','I_p_[m^4]','k_x_[-]','k_y_[-]','A_[m^2]','pitch_[deg]','x_e_[m]','y_e_[m]']
    nr = len(r)
    dfStructure = pd.DataFrame(data=np.zeros((nr,len(columns))), columns=columns)
    dfStructure['r_[m]']     = r
    dfStructure['m_[kg/m]'] = m
    # 'x_cg_[m]','y_cg_[m]' = 0
    # 'ri_x_[m]','ri_y_[m]' # TODO
    # 'x_sh_[m]','y_sh_[m]','E_[N/m^2]',
    # dfStructure['x_sh_[m]'] = 0
    # dfStructure['y_sh_[m]']= 0
    dfStructure['E_[N/m^2]'] = E                 # Arbitrary value
    dfStructure['G_[N/m^2]'] = R/2/(1+0.3)*100   # Arbitrary value, but large (ElastoDyn stiff in torsion)
    dfStructure['I_x_[m^4]'] = EIedge / E
    dfStructure['I_y_[m^4]'] = EIflap / E
    dfStructure['I_p_[m^4]'] = np.mean(EIflap+EIedge)/2*E*100 # Ip approximated as mean value of Ix+Iy, large value (ElastoDyn stiff in torsion)
    dfStructure['k_x_[-]']   = 0
    dfStructure['k_y_[-]']   = 0
    dfStructure['A_[m^2]']   = 1       # Arbitrary value
    dfStructure['pitch_[deg]']= -twist  # TODO Taken as structural twist
    # 'x_e_[m]','y_e_[m]' =0
    # Hawc2 = BeamDyn
    #     x_cg    = -y_G
    #     y_cg    = x_G
    #     x_sh    = -y_S
    #     y_sh    = x_S
    #     x_e     = -y_C
    #     y_e     = x_C
    #     I_y     = EIx/E            # [m^4] Hawc2 Iy is wrt to principal bending ye axis
    #     I_x     = EIy/E            # [m^4] Hawc2 Ix is wrt to principal bending xe axis
    #     I_p     = GKt/G            # [m^4]
    #     k_y     = kxsGA/(G*A)
    #     k_x     = kysGA/(G*A)
    #     pitch   = theta_p*180/np.pi # [deg] NOTE: could use theta_p, theta_i or theta_s
    #     if np.all(np.abs(m)<1e-16):
    #         ri_y    = m*0
    #         ri_x    = m*0
    #     else:
    #         ri_y    = np.sqrt(Ixi/m)    # [m]
    #         ri_x    = np.sqrt(Iyi/m)    # [m]
    #     # Curvilinear position of keypoints (only used to get max radius...)
    #     dr= np.sqrt((kp_x[1:]-kp_x[0:-1])**2 +(kp_y[1:]-kp_y[0:-1])**2 +(kp_z[1:]-kp_z[0:-1])**2)
    #     r_p= np.concatenate(([0],np.cumsum(dr)))
    #     r=r_bar * r_p[-1]


    # --- Defining "c2def" meanline (NOTE: ElastoDyn has no prebend or presweep, so x&y are set to zero for now
    MMeanLine = np.column_stack((r*0, r*0, r, -twist))
    dfMeanLine = pd.DataFrame(data = MMeanLine, columns=['x_[m]','y_[m]','z_[m]','twist_[deg]'])
    # BeamDyn:
    #     Y_H2     = kp_x
    #     Z_H2     = kp_z
    #     twist_H2 = - twist*180/np.pi # -[deg]
    #     columns=['x_[m]', 'y_[m]', 'z_[m]', 'twist_[deg]']
    #     data = np.column_stack((X_H2, Y_H2, Z_H2, twist_H2))
    #     dfMeanLine = pd.DataFrame(data=data, columns=columns)

    return dfMeanLine, dfStructure
