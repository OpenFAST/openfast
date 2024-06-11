""" This module contains:
  - Polar: class to represent a polar (computes steady/unsteady parameters, corrections etc.)
  - blend: function to blend two polars
  - thicknessinterp_from_one_set: interpolate polars at different thickeness based on one set of polars 
"""

import os
import numpy as np
from .polar_file import loadPolarFile # For IO

class NoCrossingException(Exception):
    pass
class NoStallDetectedException(Exception):
    pass


class Polar(object):
    """
    Defines section lift, drag, and pitching moment coefficients as a
    function of angle of attack at a particular Reynolds number.
    Different parameters may be computed and different corrections applied.

    Available routines:
        - cl_interp         : cl at given alpha values
        - cd_interp         : cd at given alpha values
        - cm_interp         : cm at given alpha values
        - cn_interp         : cn at given alpha values
        - fs_interp       : separation function (compared to fully separated polar)
        - cl_fs_interp      : cl fully separated at given alpha values
        - cl_inv_interp     : cl inviscid at given alpha values
        - correction3D      : apply 3D rotatational correction
        - extrapolate       : extend polar data set using Viterna's method
        - unsteadyParams    : computes unsteady params e.g. needed by AeroDyn15
        - unsteadyparam     : same but (old)
        - plot              : plots the polar
        - alpha0            : computes and returns alpha0, also stored in _alpha0
        - linear_region     : determines the alpha and cl values in the linear region
        - cl_max            : cl_max
        - cl_linear_slope   : linear slope and the linear region
        - cl_fully_separated: fully separated cl
        - toAeroDyn: write AeroDyn file 
    """

    def __init__(self, filename=None, alpha=None, cl=None, cd=None, cm=None, Re=None, 
            compute_params=False, radians=None, cl_lin_method='max',
            fformat='auto', verbose=False):
        """Constructor

        Parameters
        ----------
        filename: string
            If provided, the polar will be read from filename using fformat
        Re : float
            Reynolds number
        alpha : ndarray (deg)
            angle of attack
        cl : ndarray
            lift coefficient
        cd : ndarray
            drag coefficient
        cm : ndarray
            moment coefficient
        fformat: string
            file format to be used when filename is provided

        """
        # --- Potentially-locked properties
        # Introducing locks so that some properties become readonly if prescribed by user
        self._fs_lock     = False
        self._cl_fs_lock  = False
        self._cl_inv_lock = False
        # TODO lock _alpha0 and cl_alpha
        self.fs = None  # steady separation function
        self.cl_fs = None  # cl_fully separated
        self.cl_inv = None  # cl inviscid/linear/potential flow
        self._alpha0 = None
        self._linear_slope = None

        # Read polar according to fileformat, if filename provided
        if filename is not None:
            df, Re = loadPolarFile(filename, fformat=fformat, to_radians=radians, verbose=verbose)
            alpha = df['Alpha'].values
            cl    = df['Cl'].values
            cd    = df['Cd'].values
            cm    = df['Cm'].values
            if 'fs' in df.keys():
                print('[INFO] Using separating function from input file.')
                self.fs = df['fs'].values
                self._fs_lock = True
            if 'Cl_fs' in df.keys():
                print('[INFO] Using Cl fully separated from input file.')
                self.cl_fs =df['Cl_fs'].values
                self._cl_fs_lock = True
            if 'Cl_inv' in df.keys():
                print('[INFO] Using Cl inviscid from input file.')
                self.cl_inv = df['Cl_inv'].values
                self._cl_inv_lock = True
                # TODO we need a trigger if cl_inv provided, we should get alpha0 and slope from it
            nLocks = sum([self._fs_lock, self._cl_fs_lock, self._cl_inv_lock])
            if nLocks>0 and nLocks<3:
                raise Exception("For now, input files are assumed to have all or none of the columns: (fs, cl_fs, and cl_inv). Otherwise, we\'ll have to ensure consitency, and so far we dont...")

        self.Re = Re
        self.alpha = np.array(alpha)
        if cl is None:
            cl = np.zeros_like(self.alpha)
        if cd is None:
            cd = np.zeros_like(self.alpha)
        if cm is None:
            cm = np.zeros_like(self.alpha)
        self.cl = np.array(cl)
        self.cd = np.array(cd)
        self.cm = np.array(cm)
        if radians is None:
            # If the max alpha is above pi, most likely we are in degrees
            self._radians = np.mean(np.abs(self.alpha)) <= np.pi / 2
        else:
            self._radians = radians

        # NOTE: method needs to be in harmony for linear_slope and the one used in cl_fully_separated
        if compute_params:
            self._linear_slope, self._alpha0 = self.cl_linear_slope(method=cl_lin_method)
            if not self._cl_fs_lock:
                self.cl_fully_separated(method=cl_lin_method)
            if not self._cl_inv_lock:
                self.cl_inv = self._linear_slope * (self.alpha - self._alpha0)

    def __repr__(self):
        s='<{} object>:\n'.format(type(self).__name__)
        sunit = 'deg'
        if self._radians:
            sunit = 'rad'
        s+='Parameters:\n'
        s+=' - alpha, cl, cd, cm  : arrays of size {}\n'.format(len(self.alpha))
        s+=' - Re     :            {} \n'.format(self.Re)
        s+=' - _radians:           {} (True if alpha in radians)\n'.format(self._radians)
        s+=' - _alpha0:            {} [{}]\n'.format(self._alpha0, sunit)
        s+=' - _linear_slope:      {} [1/{}]\n'.format(self._linear_slope, sunit)
        s+='Derived parameters:\n'
        s+=' * cl_inv             : array of size {} \n'.format(len(self.alpha))
        s+=' * cl_fs              : array of size {} \n'.format(len(self.alpha))
        s+=' * fs                 : array of size {} \n'.format(len(self.alpha))
        s+=' * cl_lin (UNSURE)    : array of size {} \n'.format(len(self.alpha))
        s+='Functional parameters:\n'
        s+=' * alpha0 :            {} [{}]\n'.format(self.alpha0(),sunit)
        s+=' * cl_linear_slope :   {} [1/{}]\n'.format(self.cl_linear_slope()[0],sunit)
        s+=' * cl_max :            {} \n'.format(self.cl_max())
        s+=' * unsteadyParams :    {} \n'.format(self.unsteadyParams())
        s+='Useful functions:   cl_interp, cd_interp, cm_interp, fs_interp \n'
        s+='                    cl_fs_interp, cl_inv_interp,                 \n'
        s+='                    interpolant \n'
        s+='                    plot, extrapolate\n'
        return s


    # --- Potential read only properties
    @property
    def cl_inv(self):
        if self._cl_inv is None:
            self.cl_fully_separated() # computes cl_fs, cl_inv and fs
        return self._cl_inv
    @cl_inv.setter
    def cl_inv(self, cl_inv):
        if self._cl_inv_lock:
            raise Exception('Cl_inv was set by user, cannot modify it')
        else:
            self._cl_inv = cl_inv

    @property
    def cl_fs(self):
        if self._cl_fs is None:
            self.cl_fully_separated() # computes cl_fs, cl_inv and fs
        return self._cl_fs
    @cl_fs.setter
    def cl_fs(self, cl_fs):
        if self._cl_fs_lock:
            raise Exception('cl_fs was set by user, cannot modify it')
        else:
            self._cl_fs = cl_fs

    @property
    def fs(self):
        if self._fs is None:
            self.cl_fully_separated() # computes fs, cl_inv and fs
        return self._fs
    @fs.setter
    def fs(self, fs):
        if self._fs_lock:
            raise Exception('fs was set by user, cannot modify it')
        else:
            self._fs = fs

    
    # --- Interpolants
    def cl_interp(self, alpha):
        return np.interp(alpha, self.alpha, self.cl)

    def cd_interp(self, alpha):
        return np.interp(alpha, self.alpha, self.cd)

    def cm_interp(self, alpha):
        return np.interp(alpha, self.alpha, self.cm)

    def cn_interp(self, alpha):
        return np.interp(alpha, self.alpha, self.cn)

    def fs_interp(self, alpha):
        return np.interp(alpha, self.alpha, self.fs)

    def cl_fs_interp(self, alpha):
        return np.interp(alpha, self.alpha, self.cl_fs)

    def cl_inv_interp(self, alpha):
        return np.interp(alpha, self.alpha, self.cl_inv)

    def interpolant(self, variables=['cl', 'cd', 'cm'], radians=None):
        """ 
        Create an interpolant `f` for a set of requested variables with alpha as input variable:
            var_array = f(alpha)

        This is convenient to quickly interpolate multiple polar variables at once.
        The interpolant returns an array corresponding to the interpolated values of the
        requested `variables`, in the same order as they are requested.

        When alpha is a scalar, f(alpha) is of length nVar = len(variables)
        When alpha is an array of length n, f(alpha) is of shape (nVar x n)

        INPUTS:
          - variables: list of variables that will be returned by the interpolant
                    Allowed values: ['alpha', 'cl', 'cd', 'cm', 'fs', 'cl_inv', 'cl_fs']
          - radians: enforce whether `alpha` is in radians or degrees

        OUTPUTS:
         - f: interpolant

        """
        from scipy.interpolate import interp1d

        MAP = {'alpha':self.alpha, 'cl':self.cl, 'cd':self.cd, 'cm':self.cm,
                'cl_inv':self.cl_inv, 'cl_fs':self.cl_fs, 'fs':self.fs}

        if radians is None:
            radians = self._radians

        # Create a Matrix with columns requested by user
        #polCols = polar.columns.values[1:]
        M = self.alpha # we start by alpha for convenience
        for v in variables:
            v = v.lower().strip()
            if v not in MAP.keys():
                raise Exception('Polar: cannot create an interpolant for variable `{}`, allowed variables: {}'.format(v, MAP.keys()))
            M = np.column_stack( (M, MAP[v]) )
        # We remove alpha
        M = M[:,1:]
        # Determine the "x" value for the interpolant (alpha in rad or deg)
        if radians == self._radians:
            alpha = self.alpha # the user requested the same as what we have
        else:
            if radians:
                alpha = np.radians(self.alpha)
            else:
                alpha = np.degrees(self.alpha)
        # Create the interpolant for requested variables with alpha as "x" axis
        f = interp1d(alpha, M.T)
        return f

    @property
    def cn(self):
        """ returns  : Cl cos(alpha) +  Cd      sin(alpha)
                  NOT: Cl cos(alpha) + (Cd-Cd0) sin(alpha)
        """ 
        if self._radians:
            return self.cl * np.cos(self.alpha) + self.cd * np.sin(self.alpha)
        else:
            return self.cl * np.cos(self.alpha * np.pi / 180) + self.cd * np.sin(self.alpha * np.pi / 180)

    @property
    def cl_lin(self): # TODO consider removing
        print('[WARN] Polar: cl_lin is a bit of a weird property. Not sure if it will be kept')
        if self.cl_inv is None:
            self.cl_fully_separated() # computes cl_fs, cl_inv and fs
        return self.cl_inv
        #if (self._linear_slope is None) and (self._alpha0 is None):
        #    self._linear_slope,self._alpha0=self.cl_linear_slope()
        #return self._linear_slope*(self.alpha-self._alpha0)

    @classmethod
    def fromfile(cls, filename, fformat='auto', compute_params=False, to_radians=False):
        """Constructor based on a filename
        # NOTE: this is legacy
        """
        print('[WARN] Polar: "fromfile" is depreciated and will be removed in a future release')
        return cls(filename, fformat=fformat, compute_params=compute_params, radians=to_radians)

    def correction3D(
        self,
        r_over_R,
        chord_over_r,
        tsr,
        lift_method="DuSelig",
        drag_method="None",
        blending_method="linear_25_45",
        max_cl_corr=0.25,
        alpha_max_corr=None,
        alpha_linear_min=None,
        alpha_linear_max=None,
    ):
        """Applies 3-D corrections for rotating sections from the 2-D data.

        Parameters
        ----------
        r_over_R : float
            local radial position / rotor radius
        chord_over_r : float
            local chord length / local radial location
        tsr : float
            tip-speed ratio
        lift_method : string, optional
            flag switching between Du-Selig and Snel corrections
        drag_method : string, optional
            flag switching between Eggers correction and None
        blending_method: string:
             blending method used to blend from 3D to 2D polar. default 'linear_25_45'
        max_cl_corr: float, optional
             maximum correction allowed, default is 0.25.
        alpha_max_corr : float, optional (deg)
            maximum angle of attack to apply full correction
        alpha_linear_min : float, optional (deg)
            angle of attack where linear portion of lift curve slope begins
        alpha_linear_max : float, optional (deg)
            angle of attack where linear portion of lift curve slope ends

        Returns
        -------
        polar : Polar
            A new Polar object corrected for 3-D effects

        Notes
        -----
        The Du-Selig method :cite:`Du1998A-3-D-stall-del` is used to correct lift, and
        the Eggers method :cite:`Eggers-Jr2003An-assessment-o` is used to correct drag.

        """

        if alpha_max_corr == None and alpha_linear_min == None and alpha_linear_max == None:
            alpha_linear_region, _, cl_slope, alpha0 = self.linear_region()
            alpha_linear_min = alpha_linear_region[0]
            alpha_linear_max = alpha_linear_region[-1]
            _, alpha_max_corr = self.cl_max()
            find_linear_region = False
        elif alpha_max_corr * alpha_linear_min * alpha_linear_max == None:
            raise Exception(
                "Define all or none of the keyword arguments alpha_max_corr, alpha_linear_min, and alpha_linear_max"
            )
        else:
            find_linear_region = True

        # rename and convert units for convenience
        if self._radians:
            alpha = alpha
        else:
            alpha = np.radians(self.alpha)
        cl_2d = self.cl
        cd_2d = self.cd
        alpha_max_corr = np.radians(alpha_max_corr)
        alpha_linear_min = np.radians(alpha_linear_min)
        alpha_linear_max = np.radians(alpha_linear_max)

        # parameters in Du-Selig model
        a = 1
        b = 1
        d = 1
        lam = tsr / (1 + tsr ** 2) ** 0.5  # modified tip speed ratio
        if np.abs(r_over_R)>1e-4:
            expon = d / lam / r_over_R
        else:
            expon = d / lam / 1e-4

        # find linear region with numpy polyfit
        if find_linear_region:
            idx = np.logical_and(alpha >= alpha_linear_min, alpha <= alpha_linear_max)
            p = np.polyfit(alpha[idx], cl_2d[idx], 1)
            cl_slope = p[0]
            alpha0 = -p[1] / cl_slope
        else:
            cl_slope = np.degrees(cl_slope)
            alpha0 = np.radians(alpha0)

        if lift_method == "DuSelig":
            # Du-Selig correction factor
            if np.abs(cl_slope)>1e-4:
                fcl = (
                    1.0
                    / cl_slope
                    * (1.6 * chord_over_r / 0.1267 * (a - chord_over_r ** expon) / (b + chord_over_r ** expon) - 1)
                )
                # Force fcl to stay non-negative
                if fcl < 0.:
                    fcl = 0.
            else:
                fcl=0.0
        elif lift_method == "Snel":
            # Snel correction
            fcl = 3.0 * chord_over_r ** 2.0
        else:
            raise Exception("The keyword argument lift_method (3d correction for lift) can only be DuSelig or Snel.")

        # 3D correction for lift
        cl_linear = cl_slope * (alpha - alpha0)
        cl_corr = fcl * (cl_linear - cl_2d)
        # Bound correction +/- max_cl_corr
        cl_corr = np.clip(cl_corr, -max_cl_corr, max_cl_corr)
        # Blending
        if blending_method == "linear_25_45":
            # We adjust fully between +/- 25 deg, linearly to +/- 45
            adj_alpha = np.radians([-180, -45, -25, 25, 45, 180])
            adj_value = np.array([0, 0, 1, 1, 0, 0])
            adj = np.interp(alpha, adj_alpha, adj_value)
        elif blending_method == "heaviside":
            # Apply (arbitrary!) smoothing function to smoothen the 3D corrections and zero them out away from alpha_max_corr
            delta_corr = 10
            y1 = 1.0 - smooth_heaviside(alpha, k=1, rng=(alpha_max_corr, alpha_max_corr + np.deg2rad(delta_corr)))
            y2 = smooth_heaviside(alpha, k=1, rng=(0.0, np.deg2rad(delta_corr)))
            adj = y1 * y2
        else:
            raise NotImplementedError("blending :", blending_method)
        cl_3d = cl_2d + cl_corr * adj

        # Eggers 2003 correction for drag
        if drag_method == "Eggers":
            delta_cd = cl_corr * (np.sin(alpha) - 0.12 * np.cos(alpha)) / (np.cos(alpha) + 0.12 * np.sin(alpha)) * adj
        elif drag_method == "None":
            delta_cd = 0.0
        else:
            raise Exception("The keyword argument darg_method (3d correction for drag) can only be Eggers or None.")

        cd_3d = cd_2d + delta_cd

        return type(self)(Re=self.Re, alpha=np.degrees(alpha), cl=cl_3d, cd=cd_3d, cm=self.cm, radians=False)

    def extrapolate(self, cdmax, AR=None, cdmin=0.001, nalpha=15):
        """Extrapolates force coefficients up to +/- 180 degrees using Viterna's method
        :cite:`Viterna1982Theoretical-and`.

        Parameters
        ----------
        cdmax : float
            maximum drag coefficient
        AR : float, optional
            aspect ratio = (rotor radius / chord_75% radius)
            if provided, cdmax is computed from AR
        cdmin: float, optional
            minimum drag coefficient.  used to prevent negative values that can sometimes occur
            with this extrapolation method
        nalpha: int, optional
            number of points to add in each segment of Viterna method

        Returns
        -------
        polar : Polar
            a new Polar object

        Notes
        -----
        If the current polar already supplies data beyond 90 degrees then
        this method cannot be used in its current form and will just return itself.

        If AR is provided, then the maximum drag coefficient is estimated as

        >>> cdmax = 1.11 + 0.018*AR


        """

        if cdmin < 0:
            raise Exception("cdmin cannot be < 0")

        # lift coefficient adjustment to account for assymetry
        cl_adj = 0.7

        # estimate CD max
        if AR is not None:
            cdmax = 1.11 + 0.018 * AR
        self.cdmax = max(max(self.cd), cdmax)

        # extract matching info from ends
        alpha_high = np.radians(self.alpha[-1])
        cl_high = self.cl[-1]
        cd_high = self.cd[-1]
        cm_high = self.cm[-1]

        alpha_low = np.radians(self.alpha[0])
        cl_low = self.cl[0]
        cd_low = self.cd[0]

        if alpha_high > np.pi / 2:
            raise Exception("alpha[-1] > pi/2")
            return self
        if alpha_low < -np.pi / 2:
            raise Exception("alpha[0] < -pi/2")
            return self

        # parameters used in model
        sa = np.sin(alpha_high)
        ca = np.cos(alpha_high)
        self.A = (cl_high - self.cdmax * sa * ca) * sa / ca ** 2
        self.B = (cd_high - self.cdmax * sa * sa) / ca

        # alpha_high <-> 90
        alpha1 = np.linspace(alpha_high, np.pi / 2, nalpha)
        alpha1 = alpha1[1:]  # remove first element so as not to duplicate when concatenating
        cl1, cd1 = self.__Viterna(alpha1, 1.0)

        # 90 <-> 180-alpha_high
        alpha2 = np.linspace(np.pi / 2, np.pi - alpha_high, nalpha)
        alpha2 = alpha2[1:]
        cl2, cd2 = self.__Viterna(np.pi - alpha2, -cl_adj)

        # 180-alpha_high <-> 180
        alpha3 = np.linspace(np.pi - alpha_high, np.pi, nalpha)
        alpha3 = alpha3[1:]
        cl3, cd3 = self.__Viterna(np.pi - alpha3, 1.0)
        cl3 = (alpha3 - np.pi) / alpha_high * cl_high * cl_adj  # override with linear variation

        if alpha_low <= -alpha_high:
            alpha4 = []
            cl4 = []
            cd4 = []
            alpha5max = alpha_low
        else:
            # -alpha_high <-> alpha_low
            # Note: this is done slightly differently than AirfoilPrep for better continuity
            alpha4 = np.linspace(-alpha_high, alpha_low, nalpha)
            alpha4 = alpha4[1:-2]  # also remove last element for concatenation for this case
            cl4 = -cl_high * cl_adj + (alpha4 + alpha_high) / (alpha_low + alpha_high) * (cl_low + cl_high * cl_adj)
            cd4 = cd_low + (alpha4 - alpha_low) / (-alpha_high - alpha_low) * (cd_high - cd_low)
            alpha5max = -alpha_high

        # -90 <-> -alpha_high
        alpha5 = np.linspace(-np.pi / 2, alpha5max, nalpha)
        alpha5 = alpha5[1:]
        if alpha_low == -alpha_high:
            alpha5 = alpha5[:-1]
        cl5, cd5 = self.__Viterna(-alpha5, -cl_adj)

        # -180+alpha_high <-> -90
        alpha6 = np.linspace(-np.pi + alpha_high, -np.pi / 2, nalpha)
        alpha6 = alpha6[1:]
        cl6, cd6 = self.__Viterna(alpha6 + np.pi, cl_adj)

        # -180 <-> -180 + alpha_high
        alpha7 = np.linspace(-np.pi, -np.pi + alpha_high, nalpha)
        cl7, cd7 = self.__Viterna(alpha7 + np.pi, 1.0)
        cl7 = (alpha7 + np.pi) / alpha_high * cl_high * cl_adj  # linear variation

        alpha = np.concatenate((alpha7, alpha6, alpha5, alpha4, np.radians(self.alpha), alpha1, alpha2, alpha3))
        cl = np.concatenate((cl7, cl6, cl5, cl4, self.cl, cl1, cl2, cl3))
        cd = np.concatenate((cd7, cd6, cd5, cd4, self.cd, cd1, cd2, cd3))

        cd = np.maximum(cd, cdmin)  # don't allow negative drag coefficients

        # Setup alpha and cm to be used in extrapolation
        cm1_alpha = np.floor(self.alpha[0] / 10.0) * 10.0
        cm2_alpha = np.ceil(self.alpha[-1] / 10.0) * 10.0
        if cm2_alpha == self.alpha[-1]:
            self.alpha = self.alpha[:-1]
            self.cm = self.cm[:-1]
        alpha_num = abs(int((-180.0 - cm1_alpha) / 10.0 - 1))
        alpha_cm1 = np.linspace(-180.0, cm1_alpha, alpha_num)
        alpha_cm2 = np.linspace(cm2_alpha, 180.0, int((180.0 - cm2_alpha) / 10.0 + 1))
        alpha_cm = np.concatenate(
            (alpha_cm1, self.alpha, alpha_cm2)
        )  # Specific alpha values are needed for cm function to work
        cm1 = np.zeros(len(alpha_cm1))
        cm2 = np.zeros(len(alpha_cm2))
        cm_ext = np.concatenate((cm1, self.cm, cm2))
        if np.count_nonzero(self.cm) > 0:
            cmCoef = self.__CMCoeff(cl_high, cd_high, cm_high)  # get cm coefficient
            cl_cm = np.interp(alpha_cm, np.degrees(alpha), cl)  # get cl for applicable alphas
            cd_cm = np.interp(alpha_cm, np.degrees(alpha), cd)  # get cd for applicable alphas
            alpha_low_deg = self.alpha[0]
            alpha_high_deg = self.alpha[-1]
            for i in range(len(alpha_cm)):
                cm_new = self.__getCM(i, cmCoef, alpha_cm, cl_cm, cd_cm, alpha_low_deg, alpha_high_deg)
                if cm_new is None:
                    pass  # For when it reaches the range of cm's that the user provides
                else:
                    cm_ext[i] = cm_new
        cm = np.interp(np.degrees(alpha), alpha_cm, cm_ext)
        return type(self)(self.Re, np.degrees(alpha), cl, cd, cm)




    def __Viterna(self, alpha, cl_adj):
        """private method to perform Viterna extrapolation"""

        alpha = np.maximum(alpha, 0.0001)  # prevent divide by zero

        cl = self.cdmax / 2 * np.sin(2 * alpha) + self.A * np.cos(alpha) ** 2 / np.sin(alpha)
        cl = cl * cl_adj

        cd = self.cdmax * np.sin(alpha) ** 2 + self.B * np.cos(alpha)

        return cl, cd

    def __CMCoeff(self, cl_high, cd_high, cm_high):
        """private method to obtain CM0 and CMCoeff"""

        found_zero_lift = False

        for i in range(len(self.cm) - 1):
            if abs(self.alpha[i]) < 20.0 and self.cl[i] <= 0 and self.cl[i + 1] >= 0:
                p = -self.cl[i] / (self.cl[i + 1] - self.cl[i])
                cm0 = self.cm[i] + p * (self.cm[i + 1] - self.cm[i])
                found_zero_lift = True
                break

        if not found_zero_lift:
            p = -self.cl[0] / (self.cl[1] - self.cl[0])
            cm0 = self.cm[0] + p * (self.cm[1] - self.cm[0])
        self.cm0 = cm0
        alpha_high = np.radians(self.alpha[-1])
        XM = (-cm_high + cm0) / (cl_high * np.cos(alpha_high) + cd_high * np.sin(alpha_high))
        cmCoef = (XM - 0.25) / np.tan((alpha_high - np.pi / 2))
        return cmCoef

    def __getCM(self, i, cmCoef, alpha, cl_ext, cd_ext, alpha_low_deg, alpha_high_deg):
        """private method to extrapolate Cm"""

        cm_new = 0
        if alpha[i] >= alpha_low_deg and alpha[i] <= alpha_high_deg:
            return
        if alpha[i] > -165 and alpha[i] < 165:
            if abs(alpha[i]) < 0.01:
                cm_new = self.cm0
            else:
                if alpha[i] > 0:
                    x = cmCoef * np.tan(np.radians(alpha[i]) - np.pi / 2) + 0.25
                    cm_new = self.cm0 - x * (
                        cl_ext[i] * np.cos(np.radians(alpha[i])) + cd_ext[i] * np.sin(np.radians(alpha[i]))
                    )
                else:
                    x = cmCoef * np.tan(-np.radians(alpha[i]) - np.pi / 2) + 0.25
                    cm_new = -(
                        self.cm0
                        - x * (-cl_ext[i] * np.cos(-np.radians(alpha[i])) + cd_ext[i] * np.sin(-np.radians(alpha[i])))
                    )
        else:
            if alpha[i] == 165:
                cm_new = -0.4
            elif alpha[i] == 170:
                cm_new = -0.5
            elif alpha[i] == 175:
                cm_new = -0.25
            elif alpha[i] == 180:
                cm_new = 0
            elif alpha[i] == -165:
                cm_new = 0.35
            elif alpha[i] == -170:
                cm_new = 0.4
            elif alpha[i] == -175:
                cm_new = 0.2
            elif alpha[i] == -180:
                cm_new = 0
            else:
                print("Angle encountered for which there is no CM table value " "(near +/-180 deg). Program will stop.")
        return cm_new

    def unsteadyParams(self, window_offset=None, nMin=720):
        """compute unsteady aero parameters used in AeroDyn input file

                TODO Questions to solve:
                  - Is alpha 0 defined at zero lift or zero Cn?
                  - Are Cn1 and Cn2 the stall points of Cn or the regions where Cn deviates from the linear region?
                  - Is Cd0 Cdmin?
                  - Should Cd0 be used in cn?
                  - Should the TSE points be used?
                  - If so, should we use the linear points or the points on the cn-curve
                  - Should we prescribe alpha0cn when determining the slope?
                NOTE:
                  alpha0Cl and alpha0Cn are usually within 0.005 deg of each other, less thatn 0.3% difference, with alpha0Cn > alpha0Cl. The difference increase thought towards the root of the blade

                  Using the f=0.7 points doesnot change much for the lower point
                            but it has quite an impact on the upper point
        %

                Parameters
                ----------
                window_dalpha0: the linear region will be looked for in the region alpha+window_offset

                Returns
                -------
                alpha0   : lift or 0 cn (TODO TODO) angle of attack (deg)
                alpha1   : angle of attack at f=0.7 (approximately the stall angle) for AOA>alpha0 (deg)
                alpha2   : angle of attack at f=0.7 (approximately the stall angle) for AOA<alpha0 (deg)
                cnSlope  : slope of 2D normal force coefficient curve (1/rad)
                Cn1      : Critical value of C0n at leading edge separation. It should be extracted from airfoil data at a given Mach and Reynolds number. It can be calculated from the static value of Cn at either the break in the pitching moment or the loss of chord force at the onset of stall. It is close to the condition of maximum lift of the airfoil at low Mach numbers.
                Cn2      : As Cn1 for negative AOAs.
                Cd0      : Drag coefficient at zero lift TODO
                Cm0      : Moment coefficient at zero lift TODO


        """
        if window_offset is None:
            dwin = np.array([-5, 10])
            if self._radians:
                dwin = np.radians(dwin)
        cl = self.cl
        cd = self.cd
        cl[np.abs(cl)<1e-10]=0
        alpha = self.alpha

        if len(alpha)<nMin:
            #print('[INFO] Polar: unsteady params, interpolating polar data to have sufficient number of points')
            # we interpolate
            alpha_lin = np.linspace(np.min(alpha), np.max(alpha), nMin)
            alpha = np.unique(np.sort(np.concatenate((alpha, alpha_lin))))
            cl = self.cl_interp(alpha)
            cd = self.cd_interp(alpha)

        if self._radians:
            cn = cl * np.cos(alpha) + cd * np.sin(alpha)
        else:
            cn = cl * np.cos(alpha * np.pi / 180) + cd * np.sin(alpha * np.pi / 180)

        # --- Zero lift
        alpha0 = self.alpha0()
        cd0 = self.cd_interp(alpha0)
        cm0 = self.cm_interp(alpha0)

        # --- Zero cn
        if self._radians:
            window = [np.radians(-20), np.radians(20)]
        else:
            window = [-20, 20]
        try:
            alpha0cn = _find_alpha0(alpha, cn, window, direction='up', value_if_constant = 0.)
        except NoCrossingException:
            print("[WARN] Polar: Cn unsteady, cannot find zero crossing with up direction, trying down direction")
            alpha0cn = _find_alpha0(alpha, cn, window, direction='down')

        # checks for inppropriate data (like cylinders)
        if len(np.unique(cl)) == 1:
            return (alpha0, 0.0, 0.0, 0.0, 0.0, 0.0, cd0, cm0)

        # --- cn "inflection" or "Max" points
        # These point are detected from slope changes of cn, positive of negative inflections
        # The upper stall point is the first point after alpha0 with a "hat" inflection
        # The lower stall point is the first point below alpha0 with a "v" inflection
        try:
            a_MaxUpp, cn_MaxUpp, a_MaxLow, cn_MaxLow = _find_max_points(alpha, cn, alpha0, method="inflections")
        except NoStallDetectedException:
            print('[WARN] Polar: Cn unsteady, cannot find stall based on inflections, using min and max')
            a_MaxUpp, cn_MaxUpp, a_MaxLow, cn_MaxLow = _find_max_points(alpha, cn, alpha0, method="minmax")

        # --- cn slope
        # Different method may be used. The max method ensures the the curve is always below its tangent
        # Leastsquare fit in the region alpha0cn+window_offset
        cnSlope_poly, a0cn_poly = _find_slope(alpha, cn, window=alpha0cn + dwin, method="leastsquare", x0=alpha0cn)
        #cnSlope_poly, a0cn_poly = _find_slope(alpha, cn, window=alpha0cn + dwin, method="leastsquare")
        # Max (KEEP ME)
        # cnSlope_max,a0cn_max = _find_slope(alpha, cn, window=[alpha0cn,a_StallUpp], method='max', xi=alpha0cn)
        # Optim
        # cnSlope_optim,a0cn_optim = _find_slope(alpha, cn, window=[alpha0-5,alpha0+20], method='optim', x0=alpha0cn)
        ## FiniteDiff
        # cnSlope_FD,a0cn_FD = _find_slope(alpha, cn, method='finitediff_1c', xi=alpha0cn)
        # slopesRel=np.array([cnSlope_poly,cnSlope_max,cnSlope_optim,cnSlope_FD])*180/np.pi/(2*np.pi)
        cnSlope = cnSlope_poly

        # --- cn at "stall onset" (Trailling Edge Separation) locations, when cn deviates from the linear region
        a_TSELow, a_TSEUpp = _find_TSE_region(alpha, cn, cnSlope, alpha0cn, deviation=0.05)
        cn_TSEUpp_lin = cnSlope * (a_TSEUpp - alpha0cn)
        cn_TSELow_lin = cnSlope * (a_TSELow - alpha0cn)
        cn_TSEUpp = np.interp(a_TSEUpp, alpha, cn)
        cn_TSELow = np.interp(a_TSELow, alpha, cn)

        # --- cn at points where f=0.7
        cn_f = cnSlope * (alpha - alpha0cn) * ((1 + np.sqrt(0.7)) / 2) ** 2
        xInter, _ = _intersections(alpha, cn_f, alpha, cn)
        #fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
        #fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        #ax.plot(alpha, cn_f    , label='cn_f')
        #ax.plot(alpha, cn    , label='cn')
        #ax.set_xlabel('')
        #ax.set_ylabel('')
        #ax.legend()
        #plt.show()
        if len(xInter) == 3:
            a_f07_Upp = xInter[2]
            a_f07_Low = xInter[0]
        else:
            print('[WARN] Polar: Cn unsteady, cn_f does not intersect cn 3 times. Intersections:{}.'.format(xInter))
            a_f07_Upp =  abs(xInter[0]) 
            a_f07_Low = -abs(xInter[0])

        # --- DEBUG plot
        #         import matplotlib.pyplot as plt
        #         plt.plot(alpha, cn,label='cn')
        #         plt.xlim([-50,50])
        #         plt.ylim([-3,3])
        #         plt.plot([alpha0-5,alpha0-5]  ,[-3,3],'k--')
        #         plt.plot([alpha0+10,alpha0+10],[-3,3],'k--')
        #         plt.plot([alpha0,alpha0],[-3,3],'r-')
        #         plt.plot([alpha0cn,alpha0cn],[-3,3],'b-')
        #
        #         plt.plot(alpha, cn_f,label='cn_f')
        #         plt.plot(a_f07_Upp,self.cn_interp(a_f07_Upp),'d',label='Cn f07 Up')
        #         plt.plot(a_f07_Low,self.cn_interp(a_f07_Low),'d',label='Cn f07 Low')
        #         plt.plot(a_TSEUpp,cn_TSEUpp,'o',label='Cn TSEUp')
        #         plt.plot(a_TSELow,cn_TSELow,'o',label='Cn TSELow')
        #         plt.plot(a_TSEUpp,cn_TSEUpp_lin,'+',label='Cn TSEUp lin')
        #         plt.plot(a_TSELow,cn_TSELow_lin,'+',label='Cn TSELow lin')
        #         plt.plot(alpha,cnSlope *(alpha-alpha0cn),'--',  label ='Linear')
        # #         plt.plot(a_MaxUpp,cnMaxUpp,'o',label='Cn MaxUp')
        # #         plt.plot(a_MaxLow,cnMaxLow,'o',label='Cn MaxLow')
        # #         plt.plot(alpha,cnSlope_poly *(alpha-a0cn_poly),'--',  label ='Polyfit   '+sSlopes[0])
        # #         plt.plot(alpha,cnSlope_max  *(alpha-a0cn_max),'--',   label ='Max       '+sSlopes[1])
        # #         plt.plot(alpha,cnSlope_optim*(alpha-a0cn_optim),'--', label ='Optim     '+sSlopes[2])
        # #         plt.plot(alpha,cnSlope_FD   *(alpha-a0cn_FD),'--',    label ='FiniteDiff'+sSlopes[3])
        # # #         plt.plot(alpha      , np.pi/180*cnSlope*(alpha-alpha0),label='cn lin')
        # # #         plt.plot(alpha1, np.pi/180*cnSlope*(alpha1-alpha0),'o',label='cn Stall')
        # # #         plt.plot(alpha2, np.pi/180*cnSlope*(alpha2-alpha0),'o',label='cn Stall')
        #         plt.legend()
        #         mng=plt.get_current_fig_manager()
        #         mng.full_screen_toggle()
        #         plt.show()
        #         raise Exception()

        # --- Deciding what we return
        # Critical value of C0n at leading edge separation
        #         cn1   = cn_TSEUpp_lin
        #         cn2   = cn_TSELow_lin
        cn1 = cn_MaxUpp
        cn2 = cn_MaxLow
        # Alpha at f=0.7
        #         alpha1= a_TSEUpp
        #         alpha2= a_TSELow
        alpha1 = a_f07_Upp
        alpha2 = a_f07_Low

        #
        if self._radians:
            alpha0 = np.degrees(alpha0)
            alpha1 = np.degrees(alpha1)
            alpha2 = np.degrees(alpha2)
            cnSlope = cnSlope
        else:
            cnSlope = cnSlope * 180 / np.pi
        # --- Sanity checks performed by OpenFAST
        deltaAlpha = 5
        if alpha0<alpha2:
            print('[WARN] Polar: alpha0<alpha2, changing alpha2..')
            alpha2 = alpha0 - deltaAlpha
            #raise Exception('alpha0 must be greater than alpha2')

        return (alpha0, alpha1, alpha2, cnSlope, cn1, cn2, cd0, cm0)

    def unsteadyparam(self, alpha_linear_min=-5, alpha_linear_max=5):
        """compute unsteady aero parameters used in AeroDyn input file

        Parameters
        ----------
        alpha_linear_min : float, optional (deg)
            angle of attack where linear portion of lift curve slope begins
        alpha_linear_max : float, optional (deg)
            angle of attack where linear portion of lift curve slope ends

        Returns
        -------
        aerodynParam : tuple of floats
            (control setting, stall angle, alpha for 0 cn, cn slope,
            cn at stall+, cn at stall-, alpha for min CD, min(CD))

        """

        alpha = np.radians(self.alpha)
        cl = self.cl
        cd = self.cd

        alpha_linear_min = np.radians(alpha_linear_min)
        alpha_linear_max = np.radians(alpha_linear_max)

        cn = cl*np.cos(alpha) + cd*np.sin(alpha)

        # find linear region
        idx = np.logical_and(alpha >= alpha_linear_min,
                             alpha <= alpha_linear_max)

        # checks for inppropriate data (like cylinders)
        if len(idx) < 10 or len(np.unique(cl)) < 10:
            return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0

        # linear fit
        p = np.polyfit(alpha[idx], cn[idx], 1)
        m = p[0]
        alpha0 = -p[1]/m

        # find cn at "stall onset" locations, when cn deviates from the linear region
        alphaUpper    = np.radians(np.arange(40.0))
        alphaLower    = np.radians(np.arange(5.0, -40.0, -1))
        cnUpper       = np.interp(alphaUpper, alpha, cn)
        cnLower       = np.interp(alphaLower, alpha, cn)
        cnLinearUpper = m*(alphaUpper - alpha0)
        cnLinearLower = m*(alphaLower - alpha0)
        deviation     = 0.05                  # threshold for cl in detecting stall

        alphaU = np.interp(deviation, cnLinearUpper-cnUpper, alphaUpper)
        alphaL = np.interp(deviation, cnLower-cnLinearLower, alphaLower)

        # compute cn at stall according to linear fit
        cnStallUpper = m*(alphaU-alpha0)
        cnStallLower = m*(alphaL-alpha0)

        # find min cd
        minIdx = cd.argmin()

        # return: control setting, stall angle, alpha for 0 cn, cn slope,
        #         cn at stall+, cn at stall-, alpha for min CD, min(CD)
        return (0.0, np.degrees(alphaU), np.degrees(alpha0), m,
                cnStallUpper, cnStallLower, alpha[minIdx], cd[minIdx])

    def plot(self):
        """plot cl/cd/cm polar

        Returns
        -------
        figs : list of figure handles

        """
        import matplotlib.pyplot as plt

        p = self

        figs = []

        # plot cl
        fig = plt.figure()
        figs.append(fig)
        ax = fig.add_subplot(111)
        plt.plot(p.alpha, p.cl, label="Re = " + str(p.Re / 1e6) + " million")
        ax.set_xlabel("angle of attack (deg)")
        ax.set_ylabel("lift coefficient")
        ax.legend(loc="best")

        # plot cd
        fig = plt.figure()
        figs.append(fig)
        ax = fig.add_subplot(111)
        ax.plot(p.alpha, p.cd, label="Re = " + str(p.Re / 1e6) + " million")
        ax.set_xlabel("angle of attack (deg)")
        ax.set_ylabel("drag coefficient")
        ax.legend(loc="best")

        # plot cm
        fig = plt.figure()
        figs.append(fig)
        ax = fig.add_subplot(111)
        ax.plot(p.alpha, p.cm, label="Re = " + str(p.Re / 1e6) + " million")
        ax.set_xlabel("angle of attack (deg)")
        ax.set_ylabel("moment coefficient")
        ax.legend(loc="best")

        return figs

    def alpha0(self, window=None):
        """ Finds alpha0, angle of zero lift """
        if window is None:
            if self._radians:
                window = [np.radians(-30), np.radians(30)]
            else:
                window = [-30, 30]
        window = _alpha_window_in_bounds(self.alpha, window)
        # print(window)
        # print(self.alpha)
        # print(self._radians)
        # print(self.cl)
        # print(window)

        return _find_alpha0(self.alpha, self.cl, window)

    def linear_region(self, delta_alpha0=4, method_linear_fit="max"):
        cl_slope, alpha0 = self.cl_linear_slope()
        alpha_linear_region = np.asarray(_find_TSE_region(self.alpha, self.cl, cl_slope, alpha0, deviation=0.05))
        cl_linear_region = (alpha_linear_region - alpha0) * cl_slope

        return alpha_linear_region, cl_linear_region, cl_slope, alpha0

    def cl_max(self, window=None):
        """ Finds cl_max , returns (Cl_max,alpha_max) """
        if window is None:
            if self._radians:
                window = [np.radians(-40), np.radians(40)]
            else:
                window = [-40, 40]

        # Constant case or only one value
        if np.all(self.cl == self.cl[0]) or len(self.cl) == 1:
            return self.cl, self.alpha

        # Ensuring window is within our alpha values
        window = _alpha_window_in_bounds(self.alpha, window)

        # Finding max within window
        iwindow = np.where((self.alpha >= window[0]) & (self.alpha <= window[1]))
        alpha = self.alpha[iwindow]
        cl = self.cl[iwindow]
        i_max = np.argmax(cl)
        if i_max == len(iwindow):
            raise Exception(
                "Max cl is at the window boundary ([{};{}]), increase window (TODO automatically)".format(
                    window[0], window[1]
                )
            )
            pass
        cl_max = cl[i_max]
        alpha_cl_max = alpha[i_max]
        #         alpha_zc,i_zc = _zero_crossings(x=alpha,y=cl,direction='up')
        #         if len(alpha_zc)>1:
        #             raise Exception('Cannot find alpha0, {} zero crossings of Cl in the range of alpha values: [{} {}] '.format(len(alpha_zc),window[0],window[1]))
        #         elif len(alpha_zc)==0:
        #             raise Exception('Cannot find alpha0, no zero crossing of Cl in the range of alpha values: [{} {}] '.format(window[0],window[1]))
        #
        #         alpha0=alpha_zc[0]
        return cl_max, alpha_cl_max


    def cl_linear_slope(self, window=None, method="optim", radians=False):
        """Find slope of linear region
        Outputs: a 2-tuplet of:
           slope (in inverse units of alpha, or in radians-1 if radians=True)
           alpha_0 in the same unit as alpha, or in radians if radians=True
        """
        return cl_linear_slope(self.alpha, self.cl, window=window, method=method, inputInRadians=self._radians, radians=radians)

    def cl_fully_separated(self, method='max'):
        alpha0 = self.alpha0()
        cla,_, = self.cl_linear_slope(method=method)
        if cla == 0:
            cl_fs = self.cl  # when fs ==1
            fs = self.cl * 0
        else:
            cl_ratio = self.cl / (cla * (self.alpha - alpha0))
            cl_ratio[np.where(cl_ratio < 0)] = 0
            if self._fs_lock:
                fs = self.fs
            else:
                fs = (2 * np.sqrt(cl_ratio) - 1) ** 2
                fs[np.where(fs < 1e-15)] = 0
            # Initialize to linear region (in fact only at singularity, where fs=1)
            cl_fs = self.cl / 2.0  # when fs ==1
            # Region where fs<1, merge
            I = np.where(fs < 1)
            cl_fs[I] = (self.cl[I] - cla * (self.alpha[I] - alpha0) * fs[I]) / (1.0 - fs[I])
            # Outside region, use steady data
            iHig = np.ma.argmin(np.ma.MaskedArray(fs, self.alpha < alpha0))
            iLow = np.ma.argmin(np.ma.MaskedArray(fs, self.alpha > alpha0))
            cl_fs[0 : iLow + 1] = self.cl[0 : iLow + 1]
            cl_fs[iHig + 1 : -1] = self.cl[iHig + 1 : -1]

        # Ensuring everything is consistent (but we cant with user provided values..)
        cl_inv = cla * (self.alpha - alpha0)
        if not self._fs_lock:
            fs = (self.cl - cl_fs) / (cl_inv - cl_fs + 1e-10)
            fs[np.where(fs < 1e-15)] = 0
            fs[np.where(fs > 1)] = 1
            # Storing
            self.fs = fs
        self.cl_fs = cl_fs
        if not self._cl_inv_lock:
            self.cl_inv = cl_inv
        return cl_fs, fs

    def dynaStallOye_DiscreteStep(self, alpha_t, tau, fs_prev, dt):
        # compute aerodynamical force from aerodynamic data
        # interpolation from data
        fs = self.fs_interp(alpha_t)
        Clinv = self.cl_inv_interp(alpha_t)
        Clfs = self.cl_fs_interp(alpha_t)
        # dynamic stall model
        fs_dyn = fs + (fs_prev - fs) * np.exp(-dt / tau)
        Cl = fs_dyn * Clinv + (1 - fs_dyn) * Clfs
        return Cl, fs_dyn

    def toAeroDyn(self, filenameOut=None, templateFile=None, Re=1.0, comment=None, unsteadyParams=True):
        from pyFAST.input_output.fast_input_file import ADPolarFile
        cleanComments=comment is not None
        # Read a template file for AeroDyn polars
        if templateFile is None:
            MyDir=os.path.dirname(__file__)
            templateFile = os.path.join(MyDir,'../../data/NREL5MW/5MW_Baseline/Airfoils/Cylinder1.dat')
            cleanComments=True
        if isinstance(templateFile, ADPolarFile):
            ADpol = templateFile
        else:
            ADpol = ADPolarFile(templateFile)



        # --- Updating the AD polar file 
        ADpol['Re'] = Re # TODO UNKNOWN

        # Compute unsteady parameters
        if unsteadyParams:
            (alpha0,alpha1,alpha2,cnSlope,cn1,cn2,cd0,cm0)=self.unsteadyParams()

            # Setting unsteady parameters
            if np.isnan(alpha0):
                ADpol['alpha0'] = 0
            else:
                ADpol['alpha0'] = np.around(alpha0, 4)
            ADpol['alpha1']    = np.around(alpha1, 4) # TODO approximate
            ADpol['alpha2']    = np.around(alpha2, 4) # TODO approximate
            ADpol['C_nalpha']  = np.around(cnSlope ,4)
            ADpol['Cn1']       = np.around(cn1, 4)    # TODO verify
            ADpol['Cn2']       = np.around(cn2, 4)
            ADpol['Cd0']       = np.around(cd0, 4)
            ADpol['Cm0']       = np.around(cm0, 4)

        # Setting polar 
        PolarTable = np.column_stack((self.alpha, self.cl, self.cd, self.cm))
        ADpol['NumAlf'] = self.cl.shape[0]
        ADpol['AFCoeff'] = np.around(PolarTable, 5)

        # --- Comment
        if cleanComments:
            ADpol.comment='' # remove comment from template
        if comment is not None:
            ADpol.comment=comment

        if filenameOut is not None:
            ADpol.write(filenameOut)
        return ADpol




def blend(pol1, pol2, weight):
    """Blend this polar with another one with the specified weighting

    Parameters
    ----------
    pol1:  (class Polar or array) first polar
    pol2:  (class Polar or array) second polar
    weight: (float)  blending parameter between 0 (first polar) and 1 (second polar)

    Returns
    -------
    polar : (class Polar or array) a blended Polar
    """
    bReturnObject = False
    if hasattr(pol1, "cl"):
        bReturnObject = True
        alpha1 = pol1.alpha
        M1 = np.zeros((len(alpha1), 4))
        M1[:, 0] = pol1.alpha
        M1[:, 1] = pol1.cl
        M1[:, 2] = pol1.cd
        M1[:, 3] = pol1.cm
    else:
        alpha1 = pol1[:, 0]
        M1 = pol1
    if hasattr(pol2, "cl"):
        bReturnObject = True
        alpha2 = pol2.alpha
        M2 = np.zeros((len(alpha2), 4))
        M2[:, 0] = pol2.alpha
        M2[:, 1] = pol2.cl
        M2[:, 2] = pol2.cd
        M2[:, 3] = pol2.cm
    else:
        alpha2 = pol2[:, 0]
        M2 = pol2
    # Define range of alpha, merged values and truncate if one set beyond the other range
    alpha = np.union1d(alpha1, alpha2)
    min_alpha = max(alpha1.min(), alpha2.min())
    max_alpha = min(alpha1.max(), alpha2.max())
    alpha = alpha[np.logical_and(alpha >= min_alpha, alpha <= max_alpha)]
    # alpha = np.array([a for a in alpha if a >= min_alpha and a <= max_alpha])

    # Creating new output matrix to store polar
    M = np.zeros((len(alpha), M1.shape[1]))
    M[:, 0] = alpha

    # interpolate to new alpha and linearly blend
    for j in np.arange(1, M.shape[1]):
        v1 = np.interp(alpha, alpha1, M1[:, j])
        v2 = np.interp(alpha, alpha2, M2[:, j])
        M[:, j] = (1 - weight) * v1 + weight * v2
    if hasattr(pol1, "Re"):
        Re = pol1.Re + weight * (pol2.Re - pol1.Re)
    else:
        Re = np.nan

    if bReturnObject:
        return type(pol1)(Re=Re, alpha=M[:, 0], cl=M[:, 1], cd=M[:, 2], cm=M[:, 3])
    else:
        return M


def thicknessinterp_from_one_set(thickness, polarList, polarThickness):
    """Returns a set of interpolated polars from one set of polars at known thicknesses and a list of thickness
    The nearest polar is used when the thickness is beyond the range of values of the input polars.
    """
    thickness = np.asarray(thickness)
    polarThickness = np.asarray(polarThickness)
    polarList = np.asarray(polarList)
    tmax_in = np.max(thickness)
    tmax_pol = np.max(polarThickness)
    if (tmax_in > 1.2 and tmax_pol <= 1.2) or (tmax_in <= 1.2 and tmax_pol > 1.2):
        raise Exception(
            "Thicknesses of polars and input thickness need to be both in percent ([0-120]) or in fraction ([0-1.2])"
        )

    # sorting thickness
    Isort = np.argsort(polarThickness)
    polarThickness = polarThickness[Isort]
    polarList = polarList[Isort]

    polars = []
    for it, t in enumerate(thickness):
        ihigh = len(polarThickness) - 1
        for ip, tp in enumerate(polarThickness):
            if tp > t:
                ihigh = ip
                break
        ilow = 0
        for ip, tp in reversed(list(enumerate(polarThickness))):
            if tp < t:
                ilow = ip
                break

        if ihigh == ilow:
            polars.append(polarList[ihigh])
            print("[WARN] Using nearest polar for section {},   t={} , t_near={}".format(it, t, polarThickness[ihigh]))
        else:
            if (polarThickness[ilow] > t) or (polarThickness[ihigh] < t):
                raise Exception("Implementation Error")
            weight = (t - polarThickness[ilow]) / (polarThickness[ihigh] - polarThickness[ilow])
            # print(polarThickness[ilow],'<',t,'<',polarThickness[ihigh],'Weight',weight)
            pol = blend(polarList[ilow], polarList[ihigh], weight)
            polars.append(pol)
            # import matplotlib.pyplot as plt
            # fig=plt.figure()
            # plt.plot(polarList[ilow][: ,0],polarList[ilow][: ,2],'b',label='thick'+str(polarThickness[ilow]))
            # plt.plot(pol[:,0],pol[:,2],'k--',label='thick'+str(t))
            # plt.plot(polarList[ihigh][:,0],polarList[ihigh][:,2],'r',label='thick'+str(polarThickness[ihigh]))
            # plt.legend()
            # plt.show()
    return polars


def _alpha_window_in_bounds(alpha, window):
    """Ensures that the window of alpha values is within the bounds of alpha
    Example: alpha in [-30,30], window=[-20,20] => window=[-20,20]
    Example: alpha in [-10,10], window=[-20,20] => window=[-10,10]
    Example: alpha in [-30,30], window=[-40,10] => window=[-40,10]
    """
    IBef = np.where(alpha <= window[0])[0]
    if len(IBef) > 0:
        im = IBef[-1]
    else:
        im = 0
    IAft = np.where(alpha >= window[1])[0]
    if len(IAft) > 0:
        ip = IAft[0]
    else:
        ip = len(alpha) - 1
    window = [alpha[im], alpha[ip]]
    return window


def _find_alpha0(alpha, coeff, window, direction='up', value_if_constant = np.nan):
    """Finds the point where coeff(alpha)==0 using interpolation.
    The search is narrowed to a window that can be specified by the user. The default window is yet enough for cases that make physical sense.
    The angle alpha0 is found by looking at a zero up crossing in this window, and interpolation is used to find the exact location.
    """
    # Constant case or only one value
    if np.all(abs((coeff - coeff[0])<1e-8)) or len(coeff) == 1:
        if coeff[0] == 0:
            return 0
        else:
            return value_if_constant
    # Ensuring window is within our alpha values
    window = _alpha_window_in_bounds(alpha, window)

    # Finding zero up-crossing within window
    iwindow = np.where((alpha >= window[0]) & (alpha <= window[1]))
    alpha = alpha[iwindow]
    coeff = coeff[iwindow]
    alpha_zc, i_zc, s_zc = _zero_crossings(x=alpha, y=coeff, direction=direction)

    if len(alpha_zc) > 1:
        print('WARN: Cannot find alpha0, {} zero crossings of Coeff in the range of alpha values: [{} {}] '.format(len(alpha_zc),window[0],window[1]))
        print('>>> Using second zero')
        alpha_zc=alpha_zc[1:]
        #raise Exception('Cannot find alpha0, {} zero crossings of Coeff in the range of alpha values: [{} {}] '.format(len(alpha_zc),window[0],window[1]))
    elif len(alpha_zc) == 0:
        raise NoCrossingException('Cannot find alpha0, no zero crossing of Coeff in the range of alpha values: [{} {}] '.format(window[0],window[1]))

    alpha0 = alpha_zc[0]
    return alpha0


def _find_TSE_region(alpha, coeff, slope, alpha0, deviation):
    """Find the Trailing Edge Separation points, when the coefficient separates from its linear region
    These points are defined as the points where the difference is equal to +/- `deviation`
    Typically deviation is about 0.05 (absolute value)
    The linear region is defined as coeff_lin = slope (alpha-alpha0)

    returns:
       a_TSE: values of alpha at the TSE point (upper and lower)

    """
    # How off are we from the linear region
    DeltaLin = slope * (alpha - alpha0) - coeff

    # Upper and lower regions
    bUpp = alpha >= alpha0
    bLow = alpha <= alpha0

    # Finding the point where the delta is equal to `deviation`
    a_TSEUpp = np.interp(deviation, DeltaLin[bUpp], alpha[bUpp])
    a_TSELow = np.interp(-deviation, DeltaLin[bLow], alpha[bLow])
    return a_TSELow, a_TSEUpp


def _find_max_points(alpha, coeff, alpha0, method="inflections"):
    """Find upper and lower max points in `coeff` vector.
    if `method` is "inflection":
       These point are detected from slope changes of `coeff`, positive of negative inflections
       The upper stall point is the first point after alpha0 with a "hat" inflection
       The lower stall point is the first point below alpha0 with a "v" inflection
    """
    if method == "inflections":
        dC = np.diff(coeff)
        IHatInflections = np.where(np.logical_and.reduce((dC[1:] < 0, dC[0:-1] > 0, alpha[1:-1] > alpha0)))[0]
        IVeeInflections = np.where(np.logical_and.reduce((dC[1:] > 0, dC[0:-1] < 0, alpha[1:-1] < alpha0)))[0]
        if len(IHatInflections) <= 0:
            raise NoStallDetectedException("Not able to detect upper stall point of curve")
        if len(IVeeInflections) <= 0:
            raise NoStallDetectedException("Not able to detect lower stall point of curve")
        a_MaxUpp = alpha[IHatInflections[0] + 1]
        c_MaxUpp = coeff[IHatInflections[0] + 1]
        a_MaxLow = alpha[IVeeInflections[-1] + 1]
        c_MaxLow = coeff[IVeeInflections[-1] + 1]
    elif method == "minmax":
        iMax = np.argmax(coeff)
        iMin = np.argmin(coeff)
        a_MaxUpp = alpha[iMax]
        c_MaxUpp = coeff[iMax]
        a_MaxLow = alpha[iMin]
        c_MaxLow = coeff[iMin]
    else:
        raise NotImplementedError()
    return (a_MaxUpp, c_MaxUpp, a_MaxLow, c_MaxLow)




# --------------------------------------------------------------------------------}
# --- Low-level functions 
# --------------------------------------------------------------------------------{
def fn_fullsep(alpha, dclda, alpha0, alpha_sl_neg, alpha_sl_pos, valpha, vCl):
    """ Function that is zero when f=0 from the Kirchhoff theory """
    cl_linear = cl_lin(alpha, dclda, alpha0, alpha_sl_neg, alpha_sl_pos, valpha, vCl)
    return cl_linear - 0.25* dclda*(alpha-alpha0)

def cl_lin(alpha, dclda, alpha0, alpha_sl_neg, alpha_sl_pos, valpha, vCl):
    """ Linear Cl """
    if alpha > alpha_sl_neg and alpha < alpha_sl_pos :
        cl = dclda*(alpha-alpha0)
    else:
        cl=np.interp(alpha, valpha, vCl)
    return cl

def cl_fullsep(alpha, dclda, alpha0, alpha_sl_neg, alpha_sl_pos, valpha, vCl, alpha_fs_l, alpha_fs_u):
    """ fully separated lift coefficient"""
    cl_linear = cl_lin(alpha, dclda, alpha0, alpha_sl_neg, alpha_sl_pos, valpha, vCl)
    if alpha > alpha_sl_neg and alpha < alpha_sl_pos :
        cl = cl_linear*0.5
    else:
        fp = f_point(alpha, dclda, alpha0, alpha_sl_neg, alpha_sl_pos, valpha, vCl, alpha_fs_l, alpha_fs_u)
        cl=(cl_linear- dclda*(alpha-alpha0)*fp)/(1-fp)
    return cl

def f_point(alpha, dclda, alpha0, alpha_sl_neg, alpha_sl_pos, valpha, vCl, alpha_fs_l, alpha_fs_u):
    """ separation function
    # TODO harmonize with cl_fully_separated maybe?
    """
    if dclda==0:
        return 0
    if alpha < alpha_fs_l or alpha > alpha_fs_u:
        return 0
    else:
        cl_linear = cl_lin(alpha, dclda, alpha0, alpha_sl_neg, alpha_sl_pos, valpha, vCl)
        xx=cl_linear/(dclda*(alpha-alpha0) + np.sign(cl_linear)*1e-15)
        return (2*np.sqrt(xx)-1)**2


def polar_params(alpha, cl, cd, cm):
    """
    - alpha in radians
    """
    #  Treat zero-lift sections separately
    zero_offset=1e-15
    x1 = -7*np.pi/180
    x2 = 7*np.pi/180
    y1 = np.interp(x1, alpha, cl)
    y2 = np.interp(x2, alpha, cl)

    if (y1<zero_offset)and(y2<zero_offset):
       alpha0        = 0.
       Cd0          = np.interp(alpha0, alpha, cd)
       Cm0          = np.interp(alpha0, alpha, cm)
       dclda        = 0.
       alpha_fs_u   = pi
       alpha_fs_l   = -pi
       alpha_sl_neg = 0.
       alpha_sl_pos = 0.
    else:
        alpha0 = _find_alpha0(alpha, cl, [np.radians(-5), np.radians(20)])

        # Find positive angle of attack stall limit alpha_sl_pos
        alpha_sl_pos=20.*np.pi/180
        move_up=True
        maxclp=0.
        nobs=50
        while move_up :
            dalpha=(alpha_sl_pos-alpha0)/nobs
            y2=np.interp(alpha_sl_pos, alpha, cl)
            dclda=y2/(alpha_sl_pos-alpha0)
            if  dclda>maxclp:
                alpha_maxclp=alpha_sl_pos
                maxclp=dclda
            relerr=0.
            for k in np.arange(nobs):
                x1 = alpha0+(k+1)*dalpha
                y1 = np.interp(x1, alpha, cl)
                y2 = dclda*(x1-alpha0)
                relerr=relerr+(y1-y2)/y2
            relerr=relerr/nobs
            move_up= relerr>1e-2
            if move_up:
                alpha_sl_pos=alpha_sl_pos-dalpha
        alpha_sl_pos=max(alpha_maxclp, alpha_sl_pos)

        # Find negative angle of attack stall limit alpha_sl_neg
        alpha_sl_neg=-20*np.pi/180
        move_down=True
        maxclp=0.
        while move_down:
            dalpha=(alpha0-alpha_sl_neg)/nobs
            y1=np.interp(alpha_sl_neg, alpha, cl)
            dclda=y1/(alpha_sl_neg-alpha0)
            if  dclda>maxclp:
                alpha_maxclp=alpha_sl_neg
                maxclp=dclda
            relerr=0.
            for k in np.arange(nobs):
                x1=alpha_sl_neg+(k+1)*dalpha
                y1 = np.interp(x1, alpha, cl)
                y2=dclda*(x1-alpha0)
                relerr=relerr+(y1-y2)/y2
            relerr=relerr/nobs
            move_down=relerr>1e-2
            if move_down:
                alpha_sl_neg=alpha_sl_neg+dalpha
        alpha_sl_neg=min(alpha_maxclp, alpha_sl_neg)

        # Compute the final alpha and dclda (linear lift coefficient slope) values
        y1=np.interp(alpha_sl_neg, alpha, cl)
        y2=np.interp(alpha_sl_pos, alpha, cl)
        alpha0=(y1*alpha_sl_pos-y2*alpha_sl_neg)/(y1-y2)
        dclda =(y1-y2)/(alpha_sl_neg-alpha_sl_pos)

        # Find Cd0 and Cm0
        Cd0          = np.interp(alpha0, alpha, cd)
        Cm0          = np.interp(alpha0, alpha, cm)
        # Find upper surface fully stalled angle of attack alpha_fs_u (Upper limit of the Kirchhoff flat plate solution)
        y1=-1.
        y2=-1.
        delta = np.pi/180
        x2=alpha_sl_pos + delta
        while  y1*y2>0. and x2+delta<np.pi:
            x1=x2
            x2=x1+np.pi/180
            y1=fn_fullsep(x1, dclda, alpha0, alpha_sl_neg, alpha_sl_pos, alpha, cl)
            y2=fn_fullsep(x2, dclda, alpha0, alpha_sl_neg, alpha_sl_pos, alpha, cl)
        if y1*y2<0:
            alpha_fs_u=(0.-y2)/(y1-y2)*x1+(0.-y1)/(y2-y1)*x2
        else:
            alpha_fs_u=np.pi
        # Find lower surface fully stalled angle of attack alpha_fs_l (lower limit of the Kirchhoff flat plate solution)
        y1=-1.
        y2=-1.
        x2=alpha_sl_neg-delta
        while y1*y2>0. and x2-delta>-np.pi:
            x1=x2
            x2=x1-delta
            y1=fn_fullsep(x1, dclda, alpha0, alpha_sl_neg, alpha_sl_pos, alpha, cl)
            y2=fn_fullsep(x2, dclda, alpha0, alpha_sl_neg, alpha_sl_pos, alpha, cl)
        if y1*y2<0:
            alpha_fs_l=(0.-y2)/(y1-y2)*x1+(0.-y1)/(y2-y1)*x2
        else:
            alpha_fs_l=-np.pi

        # --- Compute values at all angle of attack
        Cl_fully_sep = np.zeros(alpha.shape)
        fs           = np.zeros(alpha.shape)
        Cl_linear    = np.zeros(alpha.shape)
        for i,al in enumerate(alpha):
            Cl_fully_sep[i] = cl_fullsep(al, dclda, alpha0, alpha_sl_neg, alpha_sl_pos, alpha, cl, alpha_fs_l, alpha_fs_u)
            fs          [i] = f_point   (al, dclda, alpha0, alpha_sl_neg, alpha_sl_pos, alpha, cl, alpha_fs_l, alpha_fs_u)
            Cl_linear   [i] = cl_lin    (al, dclda, alpha0, alpha_sl_neg, alpha_sl_pos, alpha, cl)

    p=dict()
    p['alpha0'] = alpha0
    p['Cd0'] = Cd0
    p['Cm0'] = Cm0
    p['dclda'] = dclda
    p['alpha_fs_l'] = alpha_fs_l
    p['alpha_fs_u'] = alpha_fs_u
    p['alpha_sl_neg'] = alpha_sl_neg
    p['alpha_sl_pos'] = alpha_sl_pos

    return p, Cl_linear, Cl_fully_sep, fs



def cl_linear_slope(alpha, cl, window=None, method="max", nInterp=721, inputInRadians=False, radians=False):
    """ 
    Find slope of linear region
    Outputs: a 2-tuplet of:
       slope (in inverse units of alpha, or in radians-1 if radians=True)
       alpha_0 in the same unit as alpha, or in radians if radians=True


    INPUTS:
      - alpha: angle of attack in radians
      - Cl   : lift coefficient
      - window: [alpha_min, alpha_max]: region when linear slope is sought
      - method: 'max', 'optim', 'leastsquare', 'leastsquare_constraint'
      
    OUTPUTS:
      - Cl_alpha, alpha0: lift slope (1/rad) and angle of attack (rad) of zero lift
    """
    # --- Return function
    def myret(sl, a0):
        # wrapper function to return degrees or radians # TODO this should be a function of self._radians
        if radians:
            if inputInRadians:
                return sl, a0
            else:
                return np.rad2deg(sl), np.deg2rad(a0) # NOTE: slope needs rad2deg, alpha needs deg2rad
        else:
            return sl, a0

    # finding alpha0 # TODO TODO TODO THIS IS NOT NECESSARY
    if inputInRadians:
        windowAlpha0 = [np.radians(-30), np.radians(30)]
    else:
        windowAlpha0 = [-30, 30]
    windowAlpha0 = _alpha_window_in_bounds(alpha, windowAlpha0)
    alpha0 = _find_alpha0(alpha, cl, windowAlpha0)

    # Constant case or only one value
    if np.all(cl == cl[0]) or len(cl) == 1:
        return myret(0, alpha0)

    if window is None:
        if np.nanmin(cl) > 0 or np.nanmax(cl) < 0:
            window = [alpha[0], alpha[-1]]
        else:
            # define a window around alpha0
            if inputInRadians:
                window = alpha0 + np.radians(np.array([-5, +20]))
            else:
                window = alpha0 +            np.array([-5, +20])

    # Ensuring window is within our alpha values
    window = _alpha_window_in_bounds(alpha, window)

    if method in ["max", "leastsquare"]:
        slope, off = _find_slope(alpha, cl, xi=alpha0, window=window, method=method)

    elif method == "leastsquare_constraint":
        slope, off = _find_slope(alpha, cl, x0=alpha0, window=window, method="leastsquare")

    elif method == "optim":
        # Selecting range of values within window
        idx = np.where((alpha >= window[0]) & (alpha <= window[1]) & ~np.isnan(cl))[0]
        cl, alpha = cl[idx], alpha[idx]
        # Selecting within the min and max of this window to improve accuracy
        imin = np.where(cl == np.min(cl))[0][-1]
        idx = np.arange(imin, np.argmax(cl) + 1)
        window = [alpha[imin], alpha[np.argmax(cl)]]
        cl, alpha = cl[idx], alpha[idx]
        # Performing minimization of slope
        slope, off = _find_slope(alpha, cl, x0=alpha0, window=None, method="optim")

    else:
        raise Exception("Method unknown for lift slope determination: {}".format(method))

    # --- Safety checks
    if len(cl) > 10:
        # Looking at slope around alpha 0 to see if we are too far off
        slope_FD, off_FD = _find_slope(alpha, cl, xi=alpha0, window=window, method="finitediff_1c")
        if abs(slope - slope_FD) / slope_FD * 100 > 50:
            #raise Exception('Warning: More than 20% error between estimated slope ({:.4f}) and the slope around alpha0 ({:.4f}). The window for the slope search ([{} {}]) is likely wrong.'.format(slope,slope_FD,window[0],window[-1]))
            print('[WARN] More than 20% error between estimated slope ({:.4f}) and the slope around alpha0 ({:.4f}). The window for the slope search ([{} {}]) is likely wrong.'.format(slope,slope_FD,window[0],window[-1]))
#         print('slope ',slope,' Alpha range: {:.3f} {:.3f} - nLin {}  nMin {}  nMax {}'.format(alpha[iStart],alpha[iEnd],len(alpha[iStart:iEnd+1]),nMin,len(alpha)))
    return myret(slope, off)

# --------------------------------------------------------------------------------}
# --- Generic curve handling functions
# --------------------------------------------------------------------------------{
def _find_slope(x, y, xi=None, x0=None, window=None, method="max", opts=None, nInterp=721):
    """Find the slope of a curve at x=xi based on a given method.
    INPUTS:
    x: array of x values
    y: array of y values
    xi: point where the slope is to be computed
    x0: point where y(x0)=0
        if provided the constraint y(x0)=0 is added.
    window:
        If a `window` is provided the search is restrained to this region of x values.
        Typical windows for airfoils are: window=[alpha0,Clmax], or window=[-5,5]+alpha0
        If window is None,  the whole extent is used (window=[min(x),max(x)])

    The methods available are:
        'max'    : returns the maximum slope within the window. Needs `xi`
        'leastsquare': use leastsquare (or polyfit), to fit the curve within the window
        'finitediff_1c': first order centered finite difference. Needs `xi`
        'optim': find the slope by looking at all possible slope values, and try to find an optimal where the length of linear region is maximized.

    returns:
        (a,x0): such that the slope is a(x-x0)
                (x0=-b/a where y=ax+b)
    """
    if window is not None:
        x_=x
        y_=y
        if nInterp is not None:
            x_ = np.linspace(x[0], x[-1], max(nInterp,len(x))) # using 0.5deg resolution at least
            y_ = np.interp(x_, x, y)
        I = np.where(np.logical_and(x_>=window[0],x_<=window[1]))
        x = x_[I]
        y = y_[I]


    if len(y) <= 1:
        raise Exception('Cannot find slope, two points needed ({} after window selection)'.format(len(y)))


    if len(y)<4 and method=='optim':
        method='leastsquare'
        #print('[WARN] Not enought data to find slope with optim method, using leastsquare')


    if method == "max":
        if xi is not None:
            I = np.nonzero(x - xi)
            yi = np.interp(xi, x, y)
            a = max((y[I] - yi) / (x[I] - xi))
            x0 = xi - yi / a
        else:
            raise Exception("For now xi needs to be set to find a slope with the max method")

    elif method == "finitediff_1c":
        # First order centered finite difference
        if xi is not None:
            # First point strictly before xi
            im = np.where(x < xi)[0][-1]
            dx = x[im + 1] - x[im - 1]
            if np.abs(dx) > 1e-7:
                a = (y[im + 1] - y[im - 1]) / dx
                yi = np.interp(xi, x, y)
                x0 = xi - yi / a
            else:
                a = np.inf
                x0 = xi
            #print('a',a)
            #print('x0',x0)
            #print('yi',yi)
            dx=(x[im+1]-x[im])
            if np.abs(dx)>1e-7:
                a = ( y[im+1] - y[im] ) / dx
                yi = np.interp(xi,x,y)
                x0 = xi - yi/a
            else:
                a=np.inf
                x0 = xi
            #print('a',a)
            #print('x0',x0)
            #print('yi',yi)
        else:
            raise Exception("For now xi needs to be set to find a slope with the finite diff method")

    elif method == "leastsquare":
        if x0 is not None:
            try:
                a = np.linalg.lstsq((x - x0).reshape((-1, 1)), y.reshape((-1, 1)), rcond=None)[0][0][0]
            except:
                a = np.linalg.lstsq((x - x0).reshape((-1, 1)), y.reshape((-1, 1)))[0][0][0]
        else:
            p = np.polyfit(x, y, 1)
            a = p[0]
            x0 = -p[1] / a
    elif method == "optim":
        if opts is None:
            nMin = max(3, int(len(x) / 2))
        else:
            nMin = opts["nMin"]

        a, x0, iStart, iEnd = _find_linear_region(x, y, nMin, x0)

    else:
        raise NotImplementedError()
    return a, x0

def _find_linear_region(x, y, nMin, x0=None):
    """Find a linear region by computing all possible slopes for all possible extent.
    The objective function tries to minimize the error with the linear slope
    and maximize the length of the linear region.
    nMin is the mimum number of points to be present in the region
    If x0 is provided, the function a*(x-x0) is fitted

    returns:
        slope :
        offset:
        iStart: index of start of linear region
        iEnd  : index of end of linear region
    """
    if x0 is not None:
        x = x.reshape((-1, 1)) - x0
        y = y.reshape((-1, 1))
    n = len(x) - nMin + 1
    err = np.zeros((n, n)) * np.nan
    slp = np.zeros((n, n)) * np.nan
    off = np.zeros((n, n)) * np.nan
    spn = np.zeros((n, n)) * np.nan
    for iStart in range(n):
        for j in range(iStart, n):
            iEnd = j + nMin
            if x0 is not None:
                sl = np.linalg.lstsq(x[iStart:iEnd], y[iStart:iEnd], rcond=None)[0][0]
                slp[iStart, j] = sl
                off[iStart, j] = x0
                y_lin = x[iStart:iEnd] * sl
            else:
                coefs = np.polyfit(x[iStart:iEnd], y[iStart:iEnd], 1)
                slp[iStart, j] = coefs[0]
                off[iStart, j] = -coefs[1] / coefs[0]
                y_lin = x[iStart:iEnd] * coefs[0] + coefs[1]
            err[iStart, j] = np.mean((y[iStart:iEnd] - y_lin) ** 2)
            spn[iStart, j] = iEnd - iStart
    spn = 1 / (spn - nMin + 1)
    err = (err) / (np.nanmax(err))
    obj = np.multiply(spn, err)
    obj = err
    (iStart, j) = np.unravel_index(np.nanargmin(obj), obj.shape)
    iEnd = j + nMin - 1  # note -1 since we return the index here
    return slp[iStart, j], off[iStart, j], iStart, iEnd


def _zero_crossings(y, x=None, direction=None):

    """
    Find zero-crossing points in a discrete vector, using linear interpolation.
    direction: 'up' or 'down', to select only up-crossings or down-crossings
    Returns:
        x values xzc such that y(yzc)==0
        indexes izc, such that the zero is between y[izc] (excluded) and y[izc+1] (included)
    if direction is not provided, also returns:
            sign, equal to 1 for up crossing
    """
    y = np.asarray(y)
    if x is None:
        x = np.arange(len(y))

    deltas = x[1:] - x[0:-1]
    if np.any( deltas == 0.0):
        I=np.where(deltas==0)[0]
        print("[WARN] Some x values are repeated at index {}. Removing them.".format(I))
        x=np.delete(x,I)
        y=np.delete(x,I)
    if np.any(deltas<0):
        raise Exception("x values need to be in ascending order")

    # Indices before zero-crossing
    iBef = np.where(y[1:] * y[0:-1] < 0.0)[0]

    # Find the zero crossing by linear interpolation
    xzc = x[iBef] - y[iBef] * (x[iBef + 1] - x[iBef]) / (y[iBef + 1] - y[iBef])

    # Selecting points that are exactly 0 and where neighbor change sign
    iZero = np.where(y == 0.0)[0]
    iZero = iZero[np.where((iZero > 0) & (iZero < x.size - 1))]
    iZero = iZero[np.where(y[iZero - 1] * y[iZero + 1] < 0.0)]

    # Concatenate
    xzc = np.concatenate((xzc, x[iZero]))
    iBef = np.concatenate((iBef, iZero))

    # Sort
    iSort = np.argsort(xzc)
    xzc, iBef = xzc[iSort], iBef[iSort]

    # Return up-crossing, down crossing or both
    sign = np.sign(y[iBef + 1] - y[iBef])
    if direction == "up":
        I = np.where(sign == 1)[0]
        return xzc[I], iBef[I], sign[I]
    elif direction == "down":
        I = np.where(sign == -1)[0]
        return xzc[I], iBef[I], sign[I]
    elif direction is not None:
        raise Exception("Direction should be either `up` or `down`")
    return xzc, iBef, sign


def _intersections(x1, y1, x2, y2, plot=False, minDist=1e-6, verbose=False):
    """
    INTERSECTIONS Intersections of curves.
    Computes the (x,y) locations where two curves intersect.  The curves
    can be broken with NaNs or have vertical segments.

    Written by: Sukhbinder, https://github.com/sukhbinder/intersection
    adapted by E.Branlard to allow for minimum distance between points
    License: MIT
    usage:
        x,y=intersection(x1,y1,x2,y2)

     Example:
     a, b = 1, 2
     phi = np.linspace(3, 10, 100)
     x1 = a*phi - b*np.sin(phi)
     y1 = a - b*np.cos(phi)

     x2=phi
     y2=np.sin(phi)+2
     x,y=intersections(x1,y1,x2,y2)

     plt.plot(x1,y1,c='r')
     plt.plot(x2,y2,c='g')
     plt.plot(x,y,'*k')
     plt.show()

    """

    def _rect_inter_inner(x1, x2):
        n1 = x1.shape[0] - 1
        n2 = x2.shape[0] - 1
        X1 = np.c_[x1[:-1], x1[1:]]
        X2 = np.c_[x2[:-1], x2[1:]]
        S1 = np.tile(X1.min(axis=1), (n2, 1)).T
        S2 = np.tile(X2.max(axis=1), (n1, 1))
        S3 = np.tile(X1.max(axis=1), (n2, 1)).T
        S4 = np.tile(X2.min(axis=1), (n1, 1))
        return S1, S2, S3, S4

    def _rectangle_intersection_(x1, y1, x2, y2):
        S1, S2, S3, S4 = _rect_inter_inner(x1, x2)
        S5, S6, S7, S8 = _rect_inter_inner(y1, y2)

        C1 = np.less_equal(S1, S2)
        C2 = np.greater_equal(S3, S4)
        C3 = np.less_equal(S5, S6)
        C4 = np.greater_equal(S7, S8)

        ii, jj = np.nonzero(C1 & C2 & C3 & C4)
        return ii, jj

    ii, jj = _rectangle_intersection_(x1, y1, x2, y2)
    n = len(ii)

    dxy1 = np.diff(np.c_[x1, y1], axis=0)
    dxy2 = np.diff(np.c_[x2, y2], axis=0)

    T = np.zeros((4, n))
    AA = np.zeros((4, 4, n))
    AA[0:2, 2, :] = -1
    AA[2:4, 3, :] = -1
    AA[0::2, 0, :] = dxy1[ii, :].T
    AA[1::2, 1, :] = dxy2[jj, :].T

    BB = np.zeros((4, n))
    BB[0, :] = -x1[ii].ravel()
    BB[1, :] = -x2[jj].ravel()
    BB[2, :] = -y1[ii].ravel()
    BB[3, :] = -y2[jj].ravel()

    for i in range(n):
        try:
            T[:, i] = np.linalg.solve(AA[:, :, i], BB[:, i])
        except:
            T[:, i] = np.NaN

    in_range = (T[0, :] >= 0) & (T[1, :] >= 0) & (T[0, :] <= 1) & (T[1, :] <= 1)

    xy0 = T[2:, in_range]
    xy0 = xy0.T

    x = xy0[:, 0]
    y = xy0[:, 1]

    # --- Remove "duplicates"
    if minDist is not None:
        pointKept=[(x[0],y[0])]
        pointSkipped=[]
        for p in zip(x[1:],y[1:]):
            distances = np.array([np.sqrt((p[0]-pk[0])**2 + (p[1]-pk[1])**2) for pk in pointKept])
            if all(distances>minDist):
                pointKept.append((p[0],p[1]))
            else:
                pointSkipped.append((p[0],p[1]))
        if verbose:
            if len(pointSkipped)>0:
                print('Polar:Intersection:Point Kept    :', pointKept)
                print('Polar:Intersection:Point Skipped:', pointSkipped)

        M = np.array(pointKept)
        x = M[:,0]
        y = M[:,1]
    if plot:
        import matplotlib.pyplot as plt
        plt.plot(x1,y1,'.',c='r')
        plt.plot(x2,y2,'',c='g')
        plt.plot(x,y,'*k')

    return x, y


def smooth_heaviside(x, k=1, rng=(-np.inf, np.inf), method="exp"):
    r"""
    Smooth approximation of Heaviside function where the step occurs between rng[0] and rng[1]:
       if rng[0]<rng[1]: then  f(<rng[0])=0, f(>=rng[1])=1
       if rng[0]>rng[1]: then  f(<rng[1])=1, f(>=rng[0])=0
    exp:
       rng=(-inf,inf):  H(x)=[1 + exp(-2kx)            ]^-1
       rng=(-1,1):      H(x)=[1 + exp(4kx/(x^2-1)      ]^-1
       rng=(0,1):       H(x)=[1 + exp(k(2x-1)/(x(x-1)) ]^-1
    INPUTS:
        x  : scalar or vector of real x values \in ]-infty; infty[
        k  : float >=1, the higher k the "steeper" the heaviside function
        rng: tuple of min and max value such that f(<=min)=0  and f(>=max)=1.
             Reversing the range makes the Heaviside function from 1 to 0 instead of 0 to 1
        method: smooth approximation used (e.g. exp or tan)
    NOTE: an epsilon is introduced in the denominator to avoid overflow of the exponentail
    """
    if k < 1:
        raise Exception("k needs to be >=1")
    eps = 1e-2
    mn, mx = rng
    x = np.asarray(x)
    H = np.zeros(x.shape)
    if mn < mx:
        H[x <= mn] = 0
        H[x >= mx] = 1
        b = np.logical_and(x > mn, x < mx)
    else:
        H[x <= mx] = 1
        H[x >= mn] = 0
        b = np.logical_and(x < mn, x > mx)
    x = x[b]
    if method == "exp":
        if np.abs(mn) == np.inf and np.abs(mx) == np.inf:
            # Infinite support
            x[k * x > 100] = 100.0 / k
            x[k * x < -100] = -100.0 / k
            if mn < mx:
                H[b] = 1 / (1 + np.exp(-k * x))
            else:
                H[b] = 1 / (1 + np.exp(k * x))
        elif np.abs(mn) != np.inf and np.abs(mx) != np.inf:
            n = 4.0
            # Compact support
            s = 2.0 / (mx - mn) * (x - (mn + mx) / 2.0)  # transform compact support into ]-1,1[
            x = -n * s / (s ** 2 - 1.0)  # then transform   ]-1,1[  into ]-inf,inf[
            x[k * x > 100] = 100.0 / k
            x[k * x < -100] = -100.0 / k
            H[b] = 1.0 / (1 + np.exp(-k * x))
        else:
            raise NotImplementedError("Heaviside with only one bound infinite")
    else:
        # TODO tan approx
        raise NotImplementedError()
    return H


if __name__ == "__main__":
    pass
