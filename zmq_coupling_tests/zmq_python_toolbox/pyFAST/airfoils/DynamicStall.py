import numpy as np
import pandas as pd
from .Polar import Polar as Pol


# --------------------------------------------------------------------------------}
# --- Wagner function 
# --------------------------------------------------------------------------------{
A1_Jones, A2_Jones, b1_Jones, b2_Jones = 0.165 , 0.335 , 0.0455 , 0.3
A1_FAST,  A2_FAST,  b1_FAST,  b2_FAST  = 0.3   , 0.7   , 0.14   , 0.53

def wagner(tau_t, constants=None, A1=None, A2=None, b1=None, b2=None):
    """ 
    Lift coefficient, Cl, from Wagner function
    INPUTS:
    - tau_t: dimensionless time
    - constants: string in ['Jones', 'OpenFAST'] or None 
    - A1, A2, b1, b2 : wagner constants, should be provided if constants is None

    Reference:  Wagner - R.T Jones approximation (Jones 1938)
    """
    if constants in ['Jones','HAWC2']: # R.T Jones approximation to Wagner's function (Jones 1938)
        A1, A2, b1, b2 = A1_Jones, A2_Jones, b1_Jones, b2_Jones
    elif constants=='OpenFAST':
        A1, A2, b1, b2 = A1_FAST, A2_FAST, b1_FAST, b2_FAST
    elif constants is None:
        if all([A1, A2, b1, b2]):
            pass # all good
        else:
            raise Exception('Provide A1, A2, b1, b2 if constants is None')


    else:
        raise NotImplementedError('Constants {}'.format(constants))

    Cl = 1-A1_Jones*np.exp(-b1_Jones*tau_t)-A2_Jones*np.exp(-b2_Jones*tau_t)

    return Cl


# --------------------------------------------------------------------------------}
# --- MHH dynamic stall 
# --------------------------------------------------------------------------------{
def dynstall_mhh_sim(time, u, p, x0=None, prefix='', method='continuous'):
    """ Perform simulation using the MHH/HGM dynamic stall model
    Reference:
       Hansen Gaunaa Madsen - Riso-R-1354 (2004) A Beddoes-Lieshman type dynamic stall model

    INPUTS:
     - time: time vector
     - u: dictionary of input functions of time
     - p: dictionary of parameters
     - x0: initial conditions for the 4 states of the model. If None, the steady steady values are used
     - method: 'continuous' or 'discrete' to chose a formulation
     - prefix: prefix used for channel names of the dataframe. Use 'AB1N001' to match OpenFAST.
    OUTPUTS:
     - df: dataframe with outputs similar to UA module of OpenFAST
    """
    from scipy.integrate import solve_ivp

    # --- Initial conditions for states
    if x0 is None:
        # x0 = [0,0,0,0]
        x0 = dynstall_mhh_steady(0,u,p)

    # --- Time Integration of states
    if method=='continuous':
        sol = solve_ivp(lambda t,x: dynstall_mhh_dxdt(t,x,u,p), t_span=[time[0],time[-1]], y0=x0, t_eval=time)
        y = sol.y
    elif method=='discrete':
        y  = np.zeros((8,len(time)))
        xd = np.zeros(8)
        xd[:4] = x0
        xd[4]  = u['alpha_34'](time[0])
        xd[5]  = 0    # Cl_p
        xd[6]  = 1.0  # fp
        xd[7]  = u['U'](time[0])  # U
        y[:,0] = xd
        for it,t in enumerate(time[1:]):
            dt = t - time[it] # Note: time[it] is in fact t-dt
            xd = dynstall_mhh_update_discr(t, dt, xd, u, p)
            y[:,it+1] = xd
    else:
        raise NotImplementedError('Method: {}'.format(method))

    # --- Compute outputs
    df=pd.DataFrame()
    #print('>>> HACK time 0.002')
    df['Time_[s]'] = time
    df[prefix + 'Vrel_[m/s]']     = np.zeros(len(time))
    df[prefix + 'alpha_34_[deg]'] = np.zeros(len(time))
    df[prefix + 'Cl_[-]']         = np.zeros(len(time))
    df[prefix + 'Cd_[-]']         = np.zeros(len(time))
    df[prefix + 'Cm_[-]']         = np.zeros(len(time))
    df[prefix + 'Tu_[-]']         = np.zeros(len(time))
    df[prefix + 'alphaE_[deg]']   = np.zeros(len(time))
    df[prefix + 'alphaF_[deg]']   = np.zeros(len(time))
    df[prefix + 'Clp_[-]']   = np.zeros(len(time))
    df[prefix + 'fs_aE_[-]']   = np.zeros(len(time))
    df[prefix + 'fs_aF_[-]']   = np.zeros(len(time))
    df['torsrate']  = np.zeros(len(time))
    df['alpha_34']  = np.zeros(len(time))
    df['U']  = np.zeros(len(time))
    df['T_0']  = np.zeros(len(time))
    df['x1'] = y[0,:]
    df['x2'] = y[1,:]
    df['x3'] = y[2,:]
    df['x4'] = y[3,:]
    df['alphaE']   = np.zeros(len(time))
    df['alphaF']   = np.zeros(len(time))
    df['ClP']   = np.zeros(len(time))
    for it,t in enumerate(time):
        Cl, Cd, Cm, alphaE, Tu, fs_aE, Cl_fs, alpha_34, omega, U, alphaF, Clp, fs_aF = dynstall_mhh_outputs(t, y[:,it], u, p, calcOutput=True)
        df.loc[it,prefix + 'Vrel_[m/s]']     = U
        df.loc[it,prefix + 'alpha_34_[deg]'] = alpha_34*180/np.pi
        df.loc[it,prefix + 'Cl_[-]']         = Cl
        df.loc[it,prefix + 'Cd_[-]']         = Cd
        df.loc[it,prefix + 'Cm_[-]']         = Cm
        df.loc[it,prefix + 'Tu_[-]']         = Tu
        df.loc[it,prefix + 'alphaE_[deg]']   = alphaE*180/np.pi
        df.loc[it,prefix + 'alphaF_[deg]']   = alphaF*180/np.pi
        df.loc[it,prefix + 'omega_[deg/s]']  = omega*180/np.pi
        df.loc[it,prefix + 'Clp_[-]']   = Clp
        df.loc[it,prefix + 'fs_aF_[-]'] = fs_aF
        df.loc[it,prefix + 'fs_aE_[-]'] = fs_aE
        df.loc[it,'U']        = U
        df.loc[it,'alpha_34'] = alpha_34 
        df.loc[it,'T_0']      = Tu
        df.loc[it,'torsrate'] = omega
        df.loc[it,'alphaE']   = alphaE
        df.loc[it,'ClP']      = Clp
        df.loc[it,'alphaF']   = alphaF
    df[prefix + 'x1_[rad]'] = y[0,:]
    df[prefix + 'x2_[rad]'] = y[1,:]
    df[prefix + 'x3_[-]']   = y[2,:]
    df[prefix + 'x4_[-]']   = y[3,:]
    return df


def dynstall_mhh_param_from_polar(P, chord, Tf0=6.0, Tp0=1.5, A1=A1_Jones, A2=A2_Jones, b1=b1_Jones, b2=b2_Jones, constants='Jones', p=None):
    if not isinstance(P,Pol):
        raise Exception('Input should be an instance of the `Polar` class')
    if not P._radians :
        raise Exception('MHH dynamic stall implemented for polars in radians only')

    if constants in ['Jones','HAWC2']: 
        A1, A2, b1, b2 = A1_Jones, A2_Jones, b1_Jones, b2_Jones
        Tf0 = 6.0
        Tp0 = 1.5
    elif constants=='OpenFAST': 
        A1, A2, b1, b2 = A1_FAST, A2_FAST, b1_FAST, b2_FAST
        Tf0 = 3.0
        Tp0 = 1.7
    else:
        raise NotImplementedError('Constants {}'.format(constants))

    if p is None:
        p=dict()
    # Airfoil parameters
    p['alpha0']     = P._alpha0 # TODO TODO requires compute params
    p['Cla']        = P._linear_slope
    if p['alpha0'] is None:
        raise Exception('>>>> TODO need to compute params on polar for MHH dyn stall model')
    if p['Cla'] is None:
        raise Exception('>>>> TODO need to compute params on polar for MHH dyn stall model')
    p['chord']      = chord
    # Polar functions
    p['F_st']  = P.fs_interp
    p['Cl_fs'] = P.cl_fs_interp
    p['Cl']    = P.cl_interp
    p['Cd']    = P.cd_interp
    p['Cm']    = P.cm_interp
    # Dynamics constants
    p['Tf0'] = Tf0
    p['Tp0'] = Tp0
    p['A1']  = A1
    p['A2']  = A2
    p['b1']  = b1
    p['b2']  = b2
    p['alpha0_in_x1x2']  = True
    p['U_in_x1x2']       = False
    p['scale_x1_x2']     = False
    p['old_ClCd_dyn']    = True
    return p

def dynstall_mhh_dxdt(t,x,u,p):
    """ Time derivative of states for continous formulation """
    # Inputs
    U         = u['U'](t)
    U_dot     = u['U_dot'](t)
    omega     = u['omega'](t)
    alpha_34  = u['alpha_34'](t)
    return dynstall_mhh_dxdt_simple(t, x, U, U_dot, omega, alpha_34, p)

def dynstall_mhh_dxdt_simple(t, x, U, U_dot, omega, alpha_34, p):
    """ Time derivative of states for continous formulation """
    # States
    x1=x[0] # Downwash memory term 1
    x2=x[1] # Downwash memory term 2
    x3=x[2] # Clp', Lift coefficient with a time lag to the attached lift coeff
    x4=x[3] # f'' , Final separation point function
    # Parameters
    alpha0 = p['alpha0']
    Cla    = p['Cla']
    c      = p['chord']
    A1     = p['A1']
    A2     = p['A2']
    b1     = p['b1']
    b2     = p['b2']
    F_st   = p['F_st']
    # Variables derived from inputs
    U  = max(U, 0.01)
    Tu = max(c/(2*U), 1e-4)                                     # Eq. 23
    Tf     = p['Tf0']*Tu  # OLD was twice: Tf = p['Tf0']*c/U
    Tp     = p['Tp0']*Tu  # OLD was twice: Tp = p['Tp0']*c/U
    # Variables derived from states
    if p['alpha0_in_x1x2']:
        alphaE  = alpha_34*(1-A1-A2)+ x1 + x2  # Eq. 12
    else:
        alphaE  = (alpha_34-alpha0)*(1-A1-A2)+ x1 + x2  + alpha0  # Eq. 12

#     alphaE = u['alphaE'](t) # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< HACK HACK TODO TODO TODO TODO TODO

    Clp     = Cla * (alphaE-alpha0) + np.pi * Tu * omega      # Eq. 13
    alphaF  = x3/Cla+alpha0                                   # p. 13
    fs_aF   = F_st(alphaF)                                    # p. 13
    if(fs_aF<0):
        print('Problematic fs:',fs_aF)
    x4 = np.clip(x4, 1e-16, 1.0) # Constraining x4 between 0 and 1 increases numerical stability
    # State equation
    xdot = [0]*4
    if p['alpha0_in_x1x2']:
        xdot[0] = -1/Tu * (b1 + c * U_dot/(2*U**2)) * x1 + b1 * A1 / Tu * alpha_34
        xdot[1] = -1/Tu * (b2 + c * U_dot/(2*U**2)) * x2 + b2 * A2 / Tu * alpha_34
    else:
        xdot[0] = -1/Tu * (b1 + c * U_dot/(2*U**2)) * x1 + b1 * A1 / Tu * (alpha_34-alpha0)
        xdot[1] = -1/Tu * (b2 + c * U_dot/(2*U**2)) * x2 + b2 * A2 / Tu * (alpha_34-alpha0)
    xdot[2] = -1/Tp                             * x3 + 1/Tp * Clp
    xdot[3] = -1/Tf                             * x4 + 1/Tf * fs_aF
    return xdot



def dynstall_mhh_update_discr(t, dt, xd_old, u, p):
    """ Update discrete states 
    NOTE: discrete states include additional discrete states
    """
    # States
    x1_old       = xd_old[0] # Downwash memory term 1
    x2_old       = xd_old[1] # Downwash memory term 2
    x3_old       = xd_old[2] # Clp', Lift coefficient with a time lag to the attached lift coeff
    x4_old       = xd_old[3] # f'' , Final separation point function
    alpha_34_old = xd_old[4] # 
    Cl_p_old     = xd_old[5] # 
    fs_aF_old    = xd_old[6] # 
    U_old        = xd_old[7] # 
    xd = xd_old.copy()
    # Inputs
    U         = u['U'](t)
    U_dot     = u['U_dot'](t)
    omega     = u['omega'](t)
    alpha_34  = u['alpha_34'](t)
    # Parameters
    alpha0 = p['alpha0']
    Cla    = p['Cla']
    c      = p['chord']
    A1     = p['A1']
    A2     = p['A2']
    b1     = p['b1']
    b2     = p['b2']
    F_st   = p['F_st']
    Cl_fs  = p['Cl_fs']
    Cd     = p['Cd']
    # Variables derived from inputs
    U  = max(U, 0.01)
    Tu = max(c/(2*U), 1e-4)                                     # Eq. 23
    T1 = Tu/b1
    T2 = Tu/b2
    Tp = p['Tp0']*Tu
    Tf = p['Tf0']*Tu
    # Temporarily remove alpha0 (not necessary)
    if p['alpha0_in_x1x2']:
        alphaQS_old = alpha_34_old
        alphaQS     = alpha_34    
    else:
        alphaQS_old = alpha_34_old - alpha0
        alphaQS     = alpha_34     - alpha0

    eps=1e-4
    exp_val1=np.exp( np.clip(-dt/T1, np.log(eps), 0 ))
    exp_val2=np.exp( np.clip(-dt/T2, np.log(eps), 0 ))
    exp_val3=np.exp( np.clip(-dt/Tp, np.log(eps), 0 ))
    exp_val4=np.exp( np.clip(-dt/Tf, np.log(eps), 0 ))
    if ['scale_x1_x2']:
        xd[0] = x1_old*exp_val1 + (alpha_34*U-alpha_34_old*U_old) * A1/b1*Tu/dt*(1-exp_val1) * x4_old
        xd[1] = x2_old*exp_val2 + (alpha_34*U-alpha_34_old*U_old) * A2/b2*Tu/dt*(1-exp_val2) * x4_old
        # x1_ = x1_old*exp      + (alpha*U   -alpha_old*U_old  )   *A1/b1*Tu/dt*(1-exp)*x4

    else:
        if p['U_in_x1x2']:
            xd[0] = x1_old*exp_val1 + 0.5*(alphaQS_old+alphaQS) * A1*U*(1-exp_val1)
            xd[1] = x2_old*exp_val2 + 0.5*(alphaQS_old+alphaQS) * A2*U*(1-exp_val2)
        else:
            xd[0] = x1_old*exp_val1 + 0.5*(alphaQS_old+alphaQS) * A1*(1-exp_val1)
            xd[1] = x2_old*exp_val2 + 0.5*(alphaQS_old+alphaQS) * A2*(1-exp_val2)
    # Effective angle of attack
    if ['scale_x1_x2']:
        alphaE = alpha_34 - (xd[0]+xd[1])/U
    else:
        if p['U_in_x1x2']:
            alphaE = alphaQS*(1-A1-A2) + (xd[0]+xd[1])/U
        else:
            alphaE = alphaQS*(1-A1-A2) + (xd[0]+xd[1])  
    if p['alpha0_in_x1x2']:
        pass
    else:
        alphaE  += alpha0
        alphaQS += alpha0

#     alphaE = u['alphaE'](t) # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< HACK HACK TODO TODO TODO TODO TODO

    Cl_p   = Cla*(alphaE-alpha0)+np.pi*Tu*omega
    xd[2]  = x3_old*exp_val3 + 0.5*(Cl_p_old+Cl_p)*(1-exp_val3)
    alphaF = xd[2]/Cla + alpha0
    fs_aF  = F_st(alphaF)                                    # p. 13
    xd[3]  = x4_old*exp_val4 + 0.5*(fs_aF_old + fs_aF)*(1-exp_val4)
    # Store "old" values
    xd[4] = alpha_34
    xd[5] = Cl_p
    xd[6] = fs_aF
    xd[7] = U
    return xd

def dynstall_mhh_steady(t, u, p):
    """ Return steady state values for the 4 states of the MHH/HGM model"""
    # Inputs
    U         = u['U'](t)
    alpha_34  = u['alpha_34'](t)
    return dynstall_mhh_steady_simple(U, alpha_34, p)

def dynstall_mhh_steady_simple(U, alpha_34, p):
    # Parameters
    c      = p['chord']
    alpha0 = p['alpha0']
    Cla    = p['Cla']
    A1     = p['A1']
    A2     = p['A2']
    b1     = p['b1']
    b2     = p['b2']
    F_st   = p['F_st']
    # Variables derived from inputs
    U  = max(U, 0.01)
    Tu = max(c/(2*U), 1e-4)                                     # Eq. 23
    # Steady states
    if p['alpha0_in_x1x2']:
        x1     = A1*alpha_34
        x2     = A2*alpha_34
        alphaE = alpha_34*(1-A1-A2) + x1 + x2                   # Eq. 12
    else:
        x1     = A1*(alpha_34 - alpha0)
        x2     = A2*(alpha_34 - alpha0)
        alphaE = (alpha_34-alpha0)*(1-A1-A2) + x1 + x2 + alpha0 # Eq. 12
    x3     = Cla * (alphaE-alpha0)
    alphaF = x3/Cla+alpha0               # p. 13
    x4     = F_st(alphaF)
    return [x1,x2,x3,x4]

def dynstall_mhh_outputs(t, x, u, p, calcOutput=False):
    # Inputs
    U         = u['U'](t)
    U_dot     = u['U_dot'](t)
    alpha     = u['alpha'](t)
    omega     = u['omega'](t)
    alpha_34  = u['alpha_34'](t)
    return  dynstall_mhh_outputs_simple(t, x, U, U_dot, omega, alpha_34, p, calcOutput=calcOutput)

def dynstall_mhh_outputs_simple(t, x, U, U_dot, omega, alpha_34, p, calcOutput=False):
    # States
    x1=x[0] # Downwash memory term 1
    x2=x[1] # Downwash memory term 2
    x3=x[2] # Clp', Lift coefficient with a time lag to the attached lift coeff
    x4=x[3] # f'' , Final separation point function
    # Parameters
    alpha0 = p['alpha0']
    Cla    = p['Cla']
    c      = p['chord']
    A1     = p['A1']
    A2     = p['A2']
    b1     = p['b1']
    b2     = p['b2']
    F_st   = p['F_st']
    Cl_fs  = p['Cl_fs']
    Cl     = p['Cl']
    Cd     = p['Cd']
    Cm     = p['Cm']

    #Cd0 = fCd(alpha0)
    #a_st = ??
    # Variables derived from inputs
    U  = max(U, 0.01)
    Tu = max(c/(2*U), 1e-4)                                     # Eq. 23

    # Variables derived from states
    if p['scale_x1_x2']:
        alphaE = alpha_34 - (x[0]+x[1])/U
    else:
        if p['alpha0_in_x1x2']:
            if p['U_in_x1x2']:
                alphaE = alpha_34*(1-A1-A2)+ (x1 + x2)/U                  # Eq. 12
            else:
                alphaE = alpha_34*(1-A1-A2)+ (x1 + x2)                  # Eq. 12
        else:
            if p['U_in_x1x2']:
                alphaE = (alpha_34-alpha0)*(1-A1-A2) + (x1 + x2)/U + alpha0# Eq. 12
            else:
                alphaE = (alpha_34-alpha0)*(1-A1-A2) + x1 + x2 + alpha0# Eq. 12

#     alphaE = u['alphaE'](t) # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< HACK HACK TODO TODO TODO TODO TODO

    fs_aE    = F_st(alphaE)
    Cl_sep_e = Cl_fs(alphaE)
    Cd_e     = Cd(alphaE)
    Cm_e     = Cm(alphaE)
    x4 = np.clip(x4,0,1)
    DeltaCdfpp = (np.sqrt(fs_aE)-np.sqrt(x4))/2 - (fs_aE-x4) /4

    #ast_x4  = (fCm(x4)  - fCm(alpha0))/Cl(x4)
    #ast_faE = (fCm(fs_aE) - fCm(alpha0))/Cl(fs_aE)
    #DeltaCmfpp = (fa_st(x4) - fa_st(fs_aE))
    DeltaCmfpp = 0 # <<<<<<<<<TODO
    Cl_att_e = max(min(Cla * (alphaE-alpha0),5), -5) # Clip between -5 and 5
    if p['old_ClCd_dyn']:
        #Cd_ind  = (alpha_34-alphaE)*Cl_dyn  # <<< TODO alpha_34 or alpha_ac
        #Cd_tors =  0                        # Old
        # Outputs
        Cl_dyn =  Cl_att_e*x4 + Cl_sep_e*(1-x4) + np.pi*Tu*omega      # <<< ADDED ALPHA0?
        Cd_dyn =  Cd_e + (alpha_34-alphaE)*Cl_dyn + (Cd_e-Cd(alpha0))*DeltaCdfpp # <<< TODO alpha_34 or alpha_ac
    else:
        # Cl components
        Cl_circ = Cl_att_e*x4 + Cl_sep_e*(1-x4)
        Cl_tors = np.pi*Tu*omega
        Cl_acc  = 0 # -np.pi * Tu * yddot/U  # <<< TODO TODO for VAWT
        # Cd components
        #Cd_ind  = (alpha_34-alphaE)*Cl_dyn  # <<< TODO alpha_34 or alpha_ac
        #Cd_tors =  0                        # Old
        Cd_ind  = (alpha_34-alphaE)*Cl_circ  # New Riso-E-0171
        Cd_tors =  Cl_circ * Tu * omega      # New Riso-E-0171
        Cd_sep  = (Cd_e-Cd(alpha0))*DeltaCdfpp # <<< TODO alpha_34 or alpha_ac
        # Outputs
        Cl_dyn =  Cl_circ +  Cl_tors + Cl_acc
        Cd_dyn =  Cd_e + Cd_ind + Cd_sep + Cd_tors 
    #Cd_dyn =  Cd(alphaE) + (alpha-alphaE)*Cl(alphaE)
    #Cd_dyn =  Cd(alphaE) + Tu*omega
    Cm_dyn =  Cm_e + Cl_dyn*DeltaCmfpp - np.pi/2*Tu*omega    
    if calcOutput:
        alphaF = x3/Cla + alpha0
        Clp     = Cla * (alphaE-alpha0) + np.pi * Tu * omega      # Eq. 13
        alphaF  = x3/Cla+alpha0                                   # p. 13
        fs_aF   = F_st(alphaF)                                    # p. 13


        return Cl_dyn, Cd_dyn, Cm_dyn, alphaE, Tu, fs_aE, Cl_sep_e, alpha_34, omega, U, alphaF, Clp, fs_aF
    else:
        return Cl_dyn, Cd_dyn, Cm_dyn




# --------------------------------------------------------------------------------}
# --- Oye's dynamic stall 
# --------------------------------------------------------------------------------{
def dynstall_oye_param_from_polar(P,tau=None,tau_chord=None, p=None):
    if tau_chord is None and tau is None:
        raise Exception('Provide `tau` or provide `tau_chord`')
    if p is None:
        p=dict()
    p['tau']   = 3*tau_chord if tau is None else tau
    p['F_st']  = P.fs_interp
    p['Clinv'] = P.cl_inv_interp
    p['Clfs']  = P.cl_fs_interp
    return p

def dynstall_oye_dxdt(t,fs,u,p):
    """ d(fs)/dt = 1/tau (fs_st - fs) """
    alpha   = u['alpha'](t)
    f_st    = p['F_st'](alpha)
    return 1/p['tau'] * (f_st - fs)

def dynstall_oye_dxdt_simple(fs, fs_alpha, tau):
    """ d(fs)/dt = 1/tau (fs_st - fs) """
    return 1/tau * (fs_alpha - fs)

def dynstall_oye_steady(alpha, p):
    """ """
    return p['F_st'](alpha)

def dynstall_oye_output(t,fs,u,p):
    alpha   = u['alpha'](t)
    Clfs    = p['Clfs'](alpha)
    Clinv   = p['Clinv'](alpha)
    Cl      = fs*Clinv+(1-fs)*Clfs               
    return Cl

def dynstall_oye_output_simple(fs, Clfs, Clinv, Cl_qs, Cd_qs, Cm_qs):
    Cl      = fs*Clinv+(1-fs)*Clfs               
    return Cl, Cd_qs, Cm_qs
