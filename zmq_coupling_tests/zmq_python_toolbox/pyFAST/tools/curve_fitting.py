"""
Set of tools to fit a model to data.

The quality of a fit is usually a strong function of the initial guess. 
Because of this this package contains different kind of "helpers" and "wrapper" tools.

FUNCTIONS
---------

This package can help fitting using:
    1) High level functions, e.g. fit_sinusoid
OR using the `model_fit` function that handles:
    2) User defined "eval" model, e.g. the user sets a string '{a}*x + {b}*x**2'
    3) Predefined models, e.g. Gaussian, logarithmic, weibull_pdf, etc.
    4) Predefined fitters, e.g. SinusoidFitter, DiscretePolynomialFitter, ContinuousPolynomialFitter

1) The high level fitting functions available are:
    - fit_sinusoid
    - fit_polynomial
    - fit_gaussian

2) User defined model, using the `model_fit_function`:
    - model_fit('eval: {a} + {b}*x**3 + {c}*x**5', x, y)
    - model_fit('eval: {u_ref}*(x/{z_ref})**{alpha}', x, y, p0=(8,9,0.1), bounds=(0.001,100))
  User defined models, will require the user to provide an initial guess and potentially bounds

3) Fitting using predefined models using the `model_fit` function :
    - model_fit('predef: gaussian', x, y)
    - model_fit('predef: gaussian-yoff', x, y)
    - model_fit('predef: powerlaw_alpha', x, y, p0=(0.1), **fun_kwargs)
    - model_fit('predef: powerlaw_u_alpha', x, y, **fun_kwargs)
    - model_fit('predef: expdecay', x, y)
    - model_fit('predef: weibull_pdf', x, y)
  Predefined models have default values for bounds and guesses that can be overriden.

4) Predefined fitters, wrapped with the `model_fit` function:
    - model_fit('fitter: sinusoid', x, y) 
    - model_fit('fitter: polynomial_discrete', x, y, exponents=[0,2,4])
    - model_fit('fitter: polynomial_continuous', x, y, order=3)
  Predefined fitters can handle bounds/initial guess better

INPUTS:
--------
All functions have the following inputs:
    - x: array on the x-axis
    - y: values on the y-axis (to be fitted against a model)
Additionally some functions have the following inputs:
    - p0: initial values for parameters, either a string or a dict: 
         - string: the string is converted to a dictionary, assuming key value pairs
                 example: 'a=0, b=1.3'
         - dictionary, then keys should corresponds to the parameters of the model
               example: {'a':0, 'b':1.3}
    - bounds: bounds for each parameters, either a string or a dictionary. 
            NOTE: pi and inf are available to set bounds
       - if a string, the string is converted to a dictionary assuming key value pairs
               example: 'a=(0,3), b=(-inf,pi)'
       - if a dictionary, the keys should corresponds to the parameters of the model
               example: {'a':(0,3), 'b':(-inf,pi)}

OUTPUTS:
--------
All functions returns the same outputs:
    - y_fit : the fit to the y data
    - pfit  : the list of parameters used
    - fitter: a `ModelFitter` object useful to manipulate the fit, in particular: 
         - fitter.model: dictionary with readable versions of the parameters, formula, 
                        function to reevaluate the fit on a different x, etc.
         - fitter.data: data used for the fit
         - fitter.fit_data: perform another fit using different data

MISC
----
High-level fitters, predefined models or fitters can be added to this class.

"""
import numpy as np
import scipy.optimize as so
import scipy.stats as stats
import string
import re
from collections import OrderedDict
from numpy import sqrt, pi, exp, cos, sin, log, inf, arctan # for user convenience
import six

# --------------------------------------------------------------------------------}
# --- High level fitters 
# --------------------------------------------------------------------------------{
def fit_sinusoid(x,y,physical=False):
    """ Fits a sinusoid to y with formula:
     if physical is False: y_fit=A*sin(omega*x+phi)+B 
     if physical is True:  y_fit=A*sin(2*pi(f+phi/360))+B """
    y_fit, pfit, fitter = model_fit('fitter: sinusoid', x, y, physical=physical)
    return y_fit, pfit, fitter 

def fit_polynomial(x, y, order=None, exponents=None):
    """ Fits a polynomial to y, either:
     - full up to a given order:             y_fit= {a_i} x^i   , i=0..order
     - or using a discrete set of exponents: y_fit= {a_i} x^e[i], i=0,..len(exponents)
    OPTIONAL INPUTS:
     - order: integer
           Maximum order of polynomial, e.g. 2: for a x**0 + b x**1 + c x**2
     - exponents: array-like
           Exponents to be used. e.g. [0,2,5] for a x**0 + b x**2 + c x**5
    """
    if order is not None:
        y_fit, pfit, fitter = model_fit('fitter: polynomial_continuous', x, y, order=order)
    else:
        y_fit, pfit, fitter = model_fit('fitter: polynomial_discrete', x, y, exponents=exponents)
    return y_fit, pfit, fitter 

def fit_gaussian(x, y, offset=False):
    """ Fits a gaussin to y, with the following formula:
    offset is True :  '1/({sigma}*sqrt(2*pi)) * exp(-1/2 * ((x-{mu})/{sigma})**2)'
    offset is False:  '1/({sigma}*sqrt(2*pi)) * exp(-1/2 * ((x-{mu})/{sigma})**2) + {y0}'
    """
    if offset:
        return model_fit('predef: gaussian-yoff', x, y)
    else:
        return model_fit('predef: gaussian', x, y)

# --------------------------------------------------------------------------------}
# --- Simple mid level fitter  
# --------------------------------------------------------------------------------{
def fit_polynomial_continuous(x, y, order):
    """Fit a polynomial with a continuous set of exponents up to a given order

    Parameters
    ----------
    x,y: see `model_fit`
    order: integer
        Maximum order of polynomial, e.g. 2: for a x**0 + b x**1 + c x**2

    Returns
    -------
    see `model_fit`
    """
    pfit  = np.polyfit(x,y,order)
    y_fit = np.polyval(pfit,x)

    # coeffs_dict, e.g. {'a':xxx, 'b':xxx}, formula = 'a*x + b'
    variables    = string.ascii_lowercase[:order+1]
    coeffs_dict  = OrderedDict([(var,coeff) for i,(coeff,var) in enumerate(zip(pfit,variables))])
    formula      = ' + '.join(['{}*x**{}'.format(var,order-i) for i,var in enumerate(variables)])
    formula      = _clean_formula(formula)
    
    return y_fit,pfit,{'coeffs':coeffs_dict,'formula':formula,'fitted_function':lambda xx : np.polyval(pfit,xx)}

def fit_polynomial_discrete(x, y, exponents):
    """Fit a polynomial with a discrete set of exponents

    Parameters
    ----------
    x,y: see `model_fit`
    exponents: array-like
        Exponents to be used. e.g. [0,2,5] for a x**0 + b x**2 + c x**5

    Returns
    -------
    see `model_fit`
    """
    #exponents=-np.sort(-np.asarray(exponents))
    X_poly=np.array([])
    for i,e in enumerate(exponents):
        if i==0:
            X_poly = np.array([x**e])
        else:
            X_poly = np.vstack((X_poly,x**e))
    try:
        pfit = np.linalg.lstsq(X_poly.T, y, rcond=None)[0]
    except:
        pfit = np.linalg.lstsq(X_poly.T, y)
    y_fit= np.dot(pfit, X_poly)

    variables    = string.ascii_lowercase[:len(exponents)]
    coeffs_dict  = OrderedDict([(var,coeff) for i,(coeff,var) in enumerate(zip(pfit,variables))])
    formula      = ' + '.join(['{}*x**{}'.format(var,e) for var,e in zip(variables,exponents)])
    formula      = _clean_formula(formula)

    return y_fit,pfit,{'coeffs':coeffs_dict,'formula':formula}


def fit_powerlaw_u_alpha(x, y, z_ref=100, p0=(10,0.1)):
    """ 
    p[0] : u_ref
    p[1] : alpha
    """
    pfit, _ = so.curve_fit(lambda x, *p : p[0] * (x / z_ref) ** p[1], x, y, p0=p0)
    y_fit = pfit[0] * (x / z_ref) ** pfit[1]
    coeffs_dict=OrderedDict([('u_ref',pfit[0]),('alpha',pfit[1])])
    formula = '{u_ref} * (z / {z_ref}) ** {alpha}'
    fitted_fun = lambda xx: pfit[0] * (xx / z_ref) ** pfit[1]
    return y_fit, pfit, {'coeffs':coeffs_dict,'formula':formula,'fitted_function':fitted_fun}


def polyfit2d(x, y, z, kx=3, ky=3, order=None):
    '''
    Two dimensional polynomial fitting by least squares.
    Fits the functional form f(x,y) = z.

    Notes
    -----
    Resultant fit can be plotted with:
    np.polynomial.polynomial.polygrid2d(x, y, soln.reshape((kx+1, ky+1)))

    Parameters
    ----------
    x, y: array-like, 1d
        x and y coordinates.
    z: np.ndarray, 2d
        Surface to fit.
    kx, ky: int, default is 3
        Polynomial order in x and y, respectively.
    order: int or None, default is None
        If None, all coefficients up to maxiumum kx, ky, ie. up to and including x^kx*y^ky, are considered.
        If int, coefficients up to a maximum of kx+ky <= order are considered.

    Returns
    -------
    Return paramters from np.linalg.lstsq.

    soln: np.ndarray
        Array of polynomial coefficients.
    residuals: np.ndarray
    rank: int
    s: np.ndarray

    # The resultant fit can be visualised with:
    # 
    # fitted_surf = np.polynomial.polynomial.polyval2d(x, y, soln.reshape((kx+1,ky+1)))
    # plt.matshow(fitted_surf


    '''

    # grid coords
    x, y = np.meshgrid(x, y)
    # coefficient array, up to x^kx, y^ky
    coeffs = np.ones((kx+1, ky+1))

    # solve array
    a = np.zeros((coeffs.size, x.size))

    # for each coefficient produce array x^i, y^j
    for index, (j, i) in enumerate(np.ndindex(coeffs.shape)): # TODO should it be i,j
        # do not include powers greater than order
        if order is not None and i + j > order:
            arr = np.zeros_like(x)
        else:
            arr = coeffs[i, j] * x**i * y**j
        a[index] = arr.ravel()

    # do leastsq fitting and return leastsq result
    return np.linalg.lstsq(a.T, np.ravel(z), rcond=None)



# --------------------------------------------------------------------------------}
# --- Predifined functions NOTE: they need to be registered in variable `MODELS`
# --------------------------------------------------------------------------------{
def gaussian(x, p):
    """ p = (mu,sigma) """
    return 1/(p[1]*np.sqrt(2*np.pi)) * np.exp(-1/2*((x-p[0])/p[1])**2)

def gaussian_w_offset(x, p):
    """ p = (mu,sigma,y0) """
    return 1/(p[1]*np.sqrt(2*np.pi)) * np.exp(-1/2*((x-p[0])/p[1])**2) + p[2]

def logarithmic(x, p):
    """ p = (a,b) """
    return p[0]*np.log(x)+p[1]

def powerlaw_all(x, p):
    """ p = (alpha,u_ref,z_ref) """
    return p[1] * (x / p[2]) ** p[0]

def powerlaw_alpha(x, p, u_ref=10, z_ref=100):
    """ p = alpha """
    return u_ref * (x / z_ref) ** p[0]

def powerlaw_u_alpha(x, p, z_ref=100):
    """ p = (alpha, u_ref) """
    return p[1] * (x / z_ref) ** p[0]

def expdecay(x, p, z_ref=100):
    """ p = (A, k, B) formula: {A}*exp(-{k}*x)+{B} """,
    return p[0]* np.exp(-p[1]*x) + p[2]

def weibull_pdf(x, p, z_ref=100):
    """ p = (A, k) formula: {k}*x**({k}-1) / {A}**{k} * np.exp(-x/{A})**{k} """,
    # NOTE: if x is 0, a divide by zero error is incountered if p[1]-1<0
    p=list(p)
    return  p[1] * x ** (p[1] - 1) / p[0] ** p[1] * np.exp(-(x / p[0]) ** p[1])

def sinusoid(x, p):
    """ p = (A,omega,phi,B) """
    return p[0]*np.sin(p[1]*x+p[2]) + p[3]
def sinusoid_f(x, p):
    """ p = (A,f,phi_deg,B) """
    return p[0]*np.sin(2*pi*(p[1]*x+p[2]/360)) + p[3]



def secondorder_impulse(t, p):
    """ p = (A, omega0, zeta, B, t0) """
    A, omega0, zeta, B, t0 = p
    omegad = omega0 * sqrt(1-zeta**2)
    phi    = np.arctan2(zeta, sqrt(1-zeta**2))
    x  = np.zeros(t.shape)
    bp = t>=t0
    t  = t[bp]-t0
    x[bp] += A * sin(omegad * t) * exp(-zeta * omega0 * t)
    x+=B
    return x

def secondorder_step(t, p):
    """ p = (A, omega0, zeta, B, t0) """
    A, omega0, zeta, B, t0 = p
    omegad = omega0 * sqrt(1-zeta**2)
    phi    = np.arctan2(zeta, sqrt(1-zeta**2))
    x  = np.zeros(t.shape)
    bp = t>=t0
    t  = t[bp]-t0
    x[bp] += A * ( 1- exp(-zeta*omega0 *t)/sqrt(1-zeta**2) * cos(omegad*t - phi))
    x+=B
    return x


def gentorque(x, p):
    """ 
    INPUTS:
     x: generator or rotor speed
     p= (RtGnSp, RtTq  , Rgn2K , SlPc , SpdGenOn)
     RtGnSp  Rated generator speed for simple variable-speed generator control (HSS side) (rpm) 
     RtTq    Rated generator torque/constant generator torque in Region 3 for simple variable-speed generator control (HSS side) (N-m) 
     Rgn2K   Generator torque constant in Region 2 for simple variable-speed generator control (HSS side) (N-m/rpm^2) 
     SlPc    Rated generator slip percentage in Region 2 1/2 for simple variable-speed generator control (%) 

     OUTPUTS:
       GenTrq: Generator torque [Nm]

     """

    # Init
    RtGnSp, RtTq  , Rgn2K , SlPc, SpdGenOn = p
    GenTrq=np.zeros(x.shape)

    xmin,xmax=np.min(x), np.max(x) 
#     if RtGnSp<(xmin+xmax)*0.4:
#         return GenTrq

    # Setting up different regions
    xR21_Start = RtGnSp*(1-SlPc/100)
    bR0      = x<SpdGenOn
    bR2      = np.logical_and(x>SpdGenOn    , x<xR21_Start)
    bR21     = np.logical_and(x>=xR21_Start , x<=RtGnSp)
    bR3      = x>RtGnSp
    # R21
    y1, y2 = Rgn2K*xR21_Start**2, RtTq
    x1, x2 = xR21_Start            , RtGnSp
    m=(y2-y1)/(x2-x1)
    GenTrq[bR21] =  m*(x[bR21]-x1) + y1  # R21
    GenTrq[bR2] =  Rgn2K * x[bR2]**2  # R2
    GenTrq[bR3] =  RtTq               # R3
    return GenTrq


MODELS =[
#     {'label':'User defined model',
#          'name':'eval:',
#          'formula':'{a}*x**2 + {b}', 
#          'coeffs':None,
#          'consts':None,
#          'bounds':None },
{'label':'Gaussian', 'handle':gaussian,'id':'predef: gaussian',
'formula':'1/({sigma}*sqrt(2*pi)) * exp(-1/2 * ((x-{mu})/{sigma})**2)',
'coeffs' :'mu=0, sigma=1', # Order Important
'consts' :None,
'bounds' :None},
{'label':'Gaussian with y-offset','handle':gaussian_w_offset,'id':'predef: gaussian-yoff',
'formula':'1/({sigma}*sqrt(2*pi)) * exp(-1/2 * ((x-{mu})/{sigma})**2) + {y0}',
'coeffs' :'mu=0, sigma=1, y0=0', #Order Important
'consts' :None,
'bounds' :'sigma=(-inf,inf), mu=(-inf,inf), y0=(-inf,inf)'},
{'label':'Exponential', 'handle': expdecay, 'id':'predef: expdecay',
'formula':'{A}*exp(-{k}*x)+{B}',
'coeffs' :'A=1, k=1, B=0',  # Order Important
'consts' :None,
'bounds' :None},
{'label':'Logarithmic', 'handle': logarithmic, 'id':'predef: logarithmic',
'formula':'{a}*log(x)+{b}',
'coeffs' :'a=1, b=0',  # Order Important
'consts' :None,
'bounds' :None},
{'label':'2nd order impulse/decay (manual)', 'handle': secondorder_impulse, 'id':'predef: secondorder_impulse',
'formula':'{A}*exp(-{zeta}*{omega}*(x-{x0})) * sin({omega}*sqrt(1-{zeta}**2))) +{B}',
'coeffs' :'A=1, omega=1, zeta=0.001, B=0, x0=0',  # Order Important
'consts' :None,
'bounds' :'A=(-inf,inf), omega=(0,100), zeta=(0,1), B=(-inf,inf), x0=(-inf,inf)'},
{'label':'2nd order step (manual)', 'handle': secondorder_step, 'id':'predef: secondorder_step',
'formula':'{A}*(1-exp(-{zeta}*{omega}*(x-{x0}))/sqrt(1-{zeta}**2) * cos({omega}*sqrt(1-{zeta}**2)-arctan({zeta}/sqrt(1-{zeta}**2)))) +{B}',
'coeffs' :'A=1, omega=1, zeta=0.001, B=0, x0=0',  # Order Important
'consts' :None,
'bounds' :'A=(-inf,inf), omega=(0,100), zeta=(0,1), B=(-inf,inf), x0=(-inf,inf)'},

# --- Wind Energy
{'label':'Power law (alpha)', 'handle':powerlaw_alpha, 'id':'predef: powerlaw_alpha',
'formula':'{u_ref} * (z / {z_ref}) ** {alpha}',
'coeffs' : 'alpha=0.1',          # Order important
'consts' : 'u_ref=10, z_ref=100',
'bounds' : 'alpha=(-1,1)'},
{'label':'Power law (alpha,u)', 'handle':powerlaw_u_alpha, 'id':'predef: powerlaw_u_alpha',
'formula':'{u_ref} * (z / {z_ref}) ** {alpha}',
'coeffs': 'alpha=0.1, u_ref=10', # Order important
'consts': 'z_ref=100',
'bounds': 'u_ref=(0,inf), alpha=(-1,1)'},
# 'powerlaw_all':{'label':'Power law (alpha,u,z)', 'handle':powerlaw_all, # NOTE: not that useful
#         'formula':'{u_ref} * (z / {z_ref}) ** {alpha}',
#         'coeffs': 'alpha=0.1, u_ref=10, z_ref=100',
#         'consts': None,
#         'bounds': 'u_ref=(0,inf), alpha=(-1,1), z_ref=(0,inf)'},
{'label':'Weibull PDF', 'handle': weibull_pdf, 'id':'predef: weibull_pdf',
'formula':'{k}*x**({k}-1) / {A}**{k} * np.exp(-x/{A})**{k}',
'coeffs' :'A=1, k=1',  # Order Important
'consts' :None,
'bounds' :'A=(0.1,inf), k=(0,5)'},
{'label':'Generator Torque', 'handle': gentorque, 'id':'predef: gentorque',
'formula': '{RtGnSp} , {RtTq}  , {Rgn2K} , {SlPc} , {SpdGenOn}',
'coeffs' : 'RtGnSp=100 , RtTq=1000  , Rgn2K=0.01 ,SlPc=5 , SpdGenOn=0',  # Order Important
'consts' :None,
'bounds' :'RtGnSp=(0.1,inf) , RtTq=(1,inf), Rgn2K=(0.0,0.1) ,SlPc=(0,20) , SpdGenOn=(0,inf)'}
]

# --------------------------------------------------------------------------------}
# --- Main function wrapper
# --------------------------------------------------------------------------------{
def model_fit(func, x, y, p0=None, bounds=None, **fun_kwargs):
    """
    Parameters
    ----------
    func: string or function handle
        - function handle
        - string starting with "fitter: ": (see  variable FITTERS)
            - "fitter: polynomial_continuous 5'    : polyfit order 5
            - "fitter: polynomial_discrete  0 2 3 ': fit polynomial of exponents 0 2 3
        - string providing an expression to evaluate, e.g.: 
            - "eval: {a}*x + {b}*x**2     " 
        - string starting with "predef": (see  variable MODELS)
            - "predef: powerlaw_alpha"  :
            - "predef: powerlaw_all"    :
            - "predef: gaussian "       :

    x: array of x values
    y: array of y values
    p0: initial values for parameters, either a string or a dict: 
       - if a string: the string is converted to a dictionary, assuming key value pairs
               example: 'a=0, b=1.3'
       - if a dictionary, then keys should corresponds to the parameters of the model
               example: {'a':0, 'b':1.3}
    bounds: bounds for each parameters, either a string or a dictionary. 
            NOTE: pi and inf are available to set bounds
       - if a string, the string is converted to a dictionary assuming key value pairs
               example: 'a=(0,3), b=(-inf,pi)'
       - if a dictionary, the keys should corresponds to the parameters of the model
               example: {'a':(0,3), 'b':(-inf,pi)}
      
    Returns
    -------
    y_fit:  array with same shape as `x`
        fitted data.
    pfit : fitted parameters
    fitter: ModelFitter object
    """

    if isinstance(func,six.string_types) and func.find('fitter:')==0:
        # --- This is a high level fitter, we call the class
        # The info about the class are storred in the global variable FITTERS
        # See e.g. SinusoidFitter, DiscretePolynomialFitter
        predef_fitters=[m['id'] for m in FITTERS]
        if func not in predef_fitters:
            raise Exception('Function `{}` not defined in curve_fitting module\n Available fitters: {}'.format(func,predef_fitters))
        i          = predef_fitters.index(func)
        FitterDict = FITTERS[i]
        consts     = FITTERS[i]['consts']
        args, missing = set_common_keys(consts, fun_kwargs)
        if len(missing)>0:
            raise Exception('Curve fitting with `{}` requires the following arguments {}. Missing: {}'.format(func,consts.keys(),missing))
        # Calling the class
        fitter = FitterDict['handle'](x=x, y=y, p0=p0, bounds=bounds, **fun_kwargs)
    else:
        fitter = ModelFitter(func, x, y, p0=p0, bounds=bounds, **fun_kwargs)

    pfit   = [v for _,v in fitter.model['coeffs'].items()]
    return fitter.data['y_fit'], pfit , fitter


# --------------------------------------------------------------------------------}
# --- Main Class 
# --------------------------------------------------------------------------------{
class ModelFitter():
    def __init__(self,func=None, x=None, y=None, p0=None, bounds=None, **fun_kwargs):

        self.model={
            'name':None, 'model_function':None, 'consts':fun_kwargs, 'formula': 'unavailable', # model signature
            'coeffs':None, 'formula_num':'unavailable', 'fitted_function':None,  'coeffs_init':p0, 'bounds':bounds,  # model fitting
            'R2':None,
        }
        self.data={'x':x,'y':y,'y_fit':None}

        if func is None:
            return
        self.set_model(func, **fun_kwargs)

        # Initialize function if present
        # Perform fit if data and function is present
        if x is not None and y is not None:
            self.fit_data(x,y,p0,bounds)

    def set_model(self,func, **fun_kwargs):
        if callable(func):
            # We don't have much additional info
            self.model['model_function'] = func
            self.model['name']           = func.__name__
            pass

        elif isinstance(func,six.string_types):
            if func.find('predef:')==0:
                # --- Minimization from a predefined function
                predef_models=[m['id'] for m in MODELS]
                if func not in predef_models:
                    raise Exception('Predefined function `{}` not defined in curve_fitting module\n Available functions: {}'.format(func,predef_models))
                i = predef_models.index(func)
                ModelDict = MODELS[i]
                self.model['model_function'] = ModelDict['handle']
                self.model['name']           = ModelDict['label']
                self.model['formula']        = ModelDict['formula']
                self.model['coeffs']         = extract_key_num(ModelDict['coeffs'])
                self.model['coeffs_init']    = self.model['coeffs'].copy()
                self.model['consts']         = extract_key_num(ModelDict['consts'])
                self.model['bounds']         = extract_key_tuples(ModelDict['bounds'])

            elif func.find('eval:')==0:
                # --- Minimization from a eval string 
                formula=func[5:]
                # Extract coeffs {a} {b} {c}, replace by p[0]
                variables, formula_eval = extract_variables(formula)
                nParams=len(variables)
                if nParams==0:
                    raise Exception('Formula should contains parameters in curly brackets, e.g.: {a}, {b}, {u_1}. No parameters found in {}'.format(formula))

                # Check that the formula evaluates
                x=np.array([1,2,5])*np.sqrt(2) # some random evaluation vector..
                p=[np.sqrt(2)/4]*nParams         # some random initial conditions
                try:
                    y=eval(formula_eval)
                    y=np.asarray(y)
                    if y.shape!=x.shape:
                        raise Exception('The formula does not return an array of same size as the input variable x. The formula must include `x`: {}'.format(formula_eval))
                except SyntaxError:
                    raise Exception('The formula does not evaluate, syntax error raised: {}'.format(formula_eval))
                except ZeroDivisionError:
                    pass

                # Creating the actual function
                def func(x, p):
                    return eval(formula_eval)

                self.model['model_function'] = func
                self.model['name']           = 'user function'
                self.model['formula']        = formula
                self.model['coeffs']         = OrderedDict([(k,v) for k,v in zip(variables,p)])
                self.model['coeffs_init']    = self.model['coeffs'].copy()
                self.model['consts']         = {}
                self.model['bounds']         = None

            else:
                raise Exception('func string needs to start with `eval:` of `predef:`, func: {}'.format(func))
        else:
            raise Exception('func should be string or callable')

        if fun_kwargs is None:
            return
        if len(fun_kwargs)==0:
            return
        if self.model['consts'] is None:
            raise Exception('Fun_kwargs provided, but no function constants were defined')

        self.model['consts'], missing = set_common_keys(self.model['consts'],  fun_kwargs )
        if len(missing)>0:
            raise Exception('Curve fitting with function `{}` requires the following arguments {}. Missing: {}'.format(func.__name__,consts.keys(),missing))

    def setup_bounds(self, bounds, nParams):
        if bounds is not None:
            self.model['bounds']=bounds # store in model
        bounds=self.model['bounds'] # usemodel bounds as default
        if bounds is not None:
            if isinstance(bounds ,six.string_types): 
                bounds=extract_key_tuples(bounds)

            if isinstance(bounds ,dict): 
                if len(bounds)==0 or 'all' in bounds.keys():
                    bounds=([-np.inf]*nParams,[np.inf]*nParams)
                elif self.model['coeffs'] is not None:
                    b1=[]
                    b2=[]
                    for k in self.model['coeffs'].keys():
                        if k in bounds.keys():
                            b1.append(bounds[k][0])
                            b2.append(bounds[k][1])
                        else:
                            # TODO merge default bounds
                            raise Exception('Bounds dictionary is missing the key: `{}`'.format(k))
                    bounds=(b1,b2)
                else:
                    raise NotImplementedError('Bounds dictionary with no known model coeffs.')
            else:
                # so.curve_fit needs a 2-tuple 
                b1,b2=bounds[0],bounds[1]
                if not hasattr(b1,'__len__'):
                    b1=[b1]*nParams
                if not hasattr(b2,'__len__'):
                    b2=[b2]*nParams
                bounds=(b1,b2)
        else:
            bounds=([-np.inf]*nParams,[np.inf]*nParams)

        self.model['bounds']=bounds # store in model

    def setup_guess(self, p0, bounds, nParams):
        """ 
        Setup initial parameter values for the fit, based on what the user provided, and potentially the bounds

        INPUTS:
         - p0: initial parameter values for the fit
             - if a string (e.g. " a=1, b=3"), it's converted to a dict 
             - if a dict, the ordered keys of model['coeffs'] are used to sort p0
         - bounds: tuple of lower and upper bounds for each parameters.
                   Parameters are ordered as function of models['coeffs']
                   bounds[0]: lower bounds or all parameters
                   bounds[1]: upper bounds or all parameters

        We can assume that the bounds are set
        """
        def middleOfBounds(i):
            """ return middle of bounds for parameter `i`"""
            bLow  = bounds[0][i]
            bHigh = bounds[0][2]
            if (bLow,bHigh)==(-np.inf,np.inf):
                p_i=0
            elif bLow==-np.inf:
                p_i = -abs(bHigh)*2
            elif bHigh== np.inf:
                p_i =  abs(bLow)*2
            else:
                p_i = (bLow+bHigh)/2
            return p_i

        if isinstance(p0 ,six.string_types): 
            p0=extract_key_num(p0)
            if len(p0)==0:
                p0=None

        if p0 is None:
            # There is some tricky logic here between the priority of bounds and coeffs
            if self.model['coeffs'] is not None:
                # We rely on function to give us decent init coefficients
                p0 = ([v for _,v in self.model['coeffs'].items()])
            elif bounds is None:
                p0 = ([0]*nParams)
            else:
                # use middle of bounds
                p0 = [0]*nParams
                for i,(b1,b2) in enumerate(zip(bounds[0],bounds[1])):
                    p0[i] = middleOfBounds(i)
                p0 = (p0)
        elif isinstance(p0,dict):
            # User supplied a dictionary, we use the ordered keys of coeffs to sort p0
            p0_dict=p0.copy()
            if self.model['coeffs'] is not None:
                p0=[]
                for k in self.model['coeffs'].keys():
                    if k in p0_dict.keys():
                        p0.append(p0_dict[k])
                    else:
                        raise Exception('Guess dictionary is missing the key: `{}`'.format(k))
            else:
                raise NotImplementedError('Guess dictionary with no known model coeffs.')


        if not hasattr(p0,'__len__'):
            p0=(p0,)

        # --- Last check that p0 is within bounds
        if bounds is not None:
            for p,k,lb,ub in zip(p0, self.model['coeffs'].keys(), bounds[0], bounds[1]):
                if p<lb:
                    raise Exception('Parameter `{}` has the guess value {}, which is smaller than the lower bound ({})'.format(k,p,lb))
                if p>ub:
                    raise Exception('Parameter `{}` has the guess value {}, which is larger than the upper bound ({})'.format(k,p,ub))
                # TODO potentially set it as middle of bounds

        # --- Finally, store the initial guesses in the model
        self.model['coeffs_init'] = p0

    def fit(self, func, x, y, p0=None, bounds=None, **fun_kwargs):
        """ Fit model defined by a function to data (x,y) """
        # Setup function
        self.set_model(func, **fun_kwargs)
        # Fit data to model
        self.fit_data(x, y, p0, bounds)

    def clean_data(self,x,y):
        x=np.asarray(x)
        y=np.asarray(y)
        bNaN=~np.isnan(y)
        y=y[bNaN]
        x=x[bNaN]
        bNaN=~np.isnan(x)
        y=y[bNaN]
        x=x[bNaN]
        self.data['x']=x
        self.data['y']=y
        return x,y

    def fit_data(self, x, y, p0=None, bounds=None):
        """ fit data, assuming a model is already setup"""
        if self.model['model_function'] is None:
            raise Exception('Call set_function first')

        # Cleaning data, and store it in object
        x,y=self.clean_data(x,y)

        # nParams
        if isinstance(p0 ,six.string_types): 
            p0=extract_key_num(p0)
            if len(p0)==0:
                p0=None
        if p0 is not None:
            if hasattr(p0,'__len__'):
                nParams=len(p0)
            else:
                nParams=1
        elif self.model['coeffs'] is not None:
            nParams=len(self.model['coeffs'])
        else:
            raise Exception('Initial guess `p0` needs to be provided since we cant infer the size of the model coefficients.')
        if self.model['coeffs'] is not None:
            if len(self.model['coeffs'])!=nParams:
                raise Exception('Inconsistent dimension between model guess (size {}) and the model parameters (size {})'.format(nParams,len(self.model['coeffs'])))

        # Bounds
        self.setup_bounds(bounds,nParams)

        # Initial conditions
        self.setup_guess(p0,self.model['bounds'],nParams)

        # Fitting
        minimize_me = lambda x, *p : self.model['model_function'](x, p, **self.model['consts'])
        pfit, pcov = so.curve_fit(minimize_me, x, y, p0=self.model['coeffs_init'], bounds=self.model['bounds']) 

        # --- Reporting information about the fit (after the fit)
        y_fit = self.model['model_function'](x, pfit, **self.model['consts'])
        self.store_fit_info(y_fit, pfit)

        # --- Return a fitted function
        self.model['fitted_function'] = lambda xx: self.model['model_function'](xx, pfit, **self.model['consts'])

    def store_fit_info(self, y_fit, pfit):
        # --- Reporting information about the fit (after the fit)
        self.data['y_fit']=y_fit
        self.model['R2'] = rsquare(self.data['y'], y_fit)
        if self.model['coeffs'] is not None:
            if not isinstance(self.model['coeffs'], OrderedDict):
                raise Exception('Coeffs need to be of type OrderedDict')
            for k,v in zip(self.model['coeffs'].keys(), pfit):
                self.model['coeffs'][k]=v

        # Replace numerical values in formula
        if self.model['formula'] is not None:
            formula_num=self.model['formula']
            for k,v in self.model['coeffs'].items():
                formula_num = formula_num.replace('{'+k+'}',str(v))
            for k,v in self.model['consts'].items():
                formula_num = formula_num.replace('{'+k+'}',str(v))
            self.model['formula_num'] = formula_num

    def formula_num(self, fmt=None):
        """ return formula with coeffs and consts evaluted numerically"""
        if fmt is None:
            fmt_fun = lambda x: str(x)
        elif isinstance(fmt,six.string_types):
            fmt_fun = lambda x: ('{'+fmt+'}').format(x)
        elif callable(fmt):
            fmt_fun = fmt
        formula_num=self.model['formula']
        for k,v in self.model['coeffs'].items():
            formula_num = formula_num.replace('{'+k+'}',fmt_fun(v))
        for k,v in self.model['consts'].items():
            formula_num = formula_num.replace('{'+k+'}',fmt_fun(v))
        return formula_num



    def plot(self, x=None, fig=None, ax=None):
        if x is None:
            x=self.data['x']

        sFormula = _clean_formula(self.model['formula'],latex=True)

        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches

        if fig is None:
            fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
            fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)

        ax.plot(self.data['x'], self.data['y'], '.', label='Data')
        ax.plot(x, self.model['fitted_function'](x), '-', label='Model ' + sFormula)

        # Add extra info to the legend
        handles, labels = ax.get_legend_handles_labels() # get existing handles and labels
        empty_patch = mpatches.Patch(color='none', label='Extra label') # create a patch with no color
        for k,v in self.model['coeffs'].items():
            handles.append(empty_patch)  # add new patches and labels to list
            labels.append(r'${:s}$ = {}'.format(pretty_param(k),pretty_num_short(v)))
        handles.append(empty_patch)  # add new patches and labels to list
        labels.append('$R^2$ = {}'.format(pretty_num_short(self.model['R2'])))
        ax.legend(handles, labels)


        #ax.set_xlabel('')
        #ax.set_ylabel('')
        return fig,ax

    def print_guessbounds(self):
        s=''
        p0     = self.model['coeffs_init']
        bounds = self.model['bounds']
        for i,(k,v) in enumerate(self.model['coeffs'].items()):
            print( (pretty_num(bounds[0][i]),pretty_num(p0[i]), pretty_num(bounds[1][i])) )
            s+='{:15s}: {:10s} < {:10s} < {:10s}\n'.format(k, pretty_num(bounds[0][i]),pretty_num(p0[i]), pretty_num(bounds[1][i]))
        print(s)
            

    def __repr__(self):
        s='<{} object> with fields:\n'.format(type(self).__name__)
        s+=' - data, dictionary with keys: \n'
        s+='   - x: [{} ... {}], n: {} \n'.format(self.data['x'][0],self.data['x'][-1],len(self.data['x']))
        s+='   - y: [{} ... {}], n: {} \n'.format(self.data['y'][0],self.data['y'][-1],len(self.data['y']))
        s+=' - model, dictionary with keys: \n'
        for k,v in self.model.items():
            s=s+'   - {:15s}: {}\n'.format(k,v)
        return s


# --------------------------------------------------------------------------------}
# --- Wrapper for predefined fitters 
# --------------------------------------------------------------------------------{
class PredefinedModelFitter(ModelFitter):
    def __init__(self, x=None, y=None, p0=None, bounds=None, **kwargs):
        ModelFitter.__init__(self,x=None, y=None, p0=p0, bounds=bounds) # NOTE: not passing data

        self.kwargs=kwargs

        if x is not None and y is not None:
            self.fit_data(x,y,p0,bounds)

    def setup_model(self):
        """
        Setup model:
         - guess/coeffs_init: return params in format needed for curve_fit (p0,p1,p2,p3)
         - bound            :  bounds in format needed for curve_fit ((low0,low1,low2), (high0, high1))
         - coeffs           :  OrderedDict, necessary for user print
         - formula          :  necessary for user print
        """
        #self.model['coeffs']  = OrderedDict([(var,1) for i,var in enumerate(variables)])
        #self.model['formula'] = ''
        #self.model['coeffs_init']=p_guess
        #self.model['bounds']=bounds_guess
        raise NotImplementedError('To be implemented by child class')

    def model_function(self, x, p):
        raise NotImplementedError('To be implemented by child class')

    def fit_data(self, x, y, p0=None, bounds=None):
        # Cleaning data
        x,y=self.clean_data(x,y)

        # --- setup model
        # guess initial parameters, potential bounds, and set necessary data
        self.setup_model()

        # --- Minimization
        minimize_me = lambda x, *p : self.model_function(x, p)
        if self.model['bounds'] is None:
            pfit, pcov = so.curve_fit(minimize_me, x, y, p0=self.model['coeffs_init'])
        else:
            pfit, pcov = so.curve_fit(minimize_me, x, y, p0=self.model['coeffs_init'], bounds=self.model['bounds']) 
        # --- Reporting information about the fit (after the fit)
        # And Return a fitted function
        y_fit = self.model_function(x,  pfit)
        self.model['fitted_function']=lambda xx : self.model_function(xx, pfit)
        self.store_fit_info(y_fit, pfit)

    def plot_guess(self, x=None, fig=None, ax=None):
        """ plotthe guess values"""
        if x is None:
            x=self.data['x']
        import matplotlib.pyplot as plt
        if fig is None:
            fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
            fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)

        p_guess = self.model['coeffs_init']

        ax.plot(self.data['x'], self.data['y']   , '.', label='Data')
        ax.plot(x, self.model_function(x,p_guess), '-', label='Model at guessed parameters')
        ax.legend()


# --------------------------------------------------------------------------------}
# --- Predefined fitters
# --------------------------------------------------------------------------------{
class SecondOrderFitterImpulse(PredefinedModelFitter):

    def model_function(self, x, p):
        return secondorder_impulse(x, p)

    def setup_model(self):
        """ p = (A, omega0, zeta, B, t0) """
        self.model['coeffs'] = OrderedDict([('A',1),('omega',1),('zeta',0.01),('B',0),('t0',0)])
        self.model['formula'] = '{A}*exp(-{zeta}*{omega}*(x-{x0}))*sin({omega}*sqrt(1-{zeta}**2)))+{B}'

        # --- Guess Initial values
        x, y = self.data['x'],self.data['y']
        # TODO use signal
        dt       = x[1]-x[0]
        omega0   = main_frequency(x,y) 
        A        = np.max(y) - np.min(y)
        B        = np.mean(y)
        zeta     = 0.1
        y_start  = y[0]+0.01*A
        bDeviate = np.argwhere(abs(y-y_start)>abs(y_start-y[0]))[0]
        t0       = x[bDeviate[0]]
        p_guess  = np.array([A, omega0, zeta, B, t0])
        self.model['coeffs_init'] = p_guess
        # --- Set Bounds
        T  = x[-1]-x[0]
        dt = x[1]-x[0]
        om_min = 2*np.pi/T/2
        om_max = 2*np.pi/dt/2
        b_A    = (A*0.1,A*3)
        b_om   = (om_min,om_max)
        b_zeta = (0,1)
        b_B    = (np.min(y),np.max(y))
        b_x0   = (np.min(x),np.max(x))
        self.model['bounds'] = ((b_A[0],b_om[0],b_zeta[0],b_B[0],b_x0[0]),(b_A[1],b_om[1],b_zeta[1],b_B[1],b_x0[1]))
        #self.plot_guess(); import matplotlib.pyplot as plt; plt.show()
        #self.print_guessbounds(); 

class SecondOrderFitterStep(PredefinedModelFitter):

    def model_function(self, x, p):
        return secondorder_step(x, p)

    def setup_model(self):
        """ p = (A, omega0, zeta, B, t0) """
        self.model['coeffs'] = OrderedDict([('A',1),('omega',1),('zeta',0.01),('B',0),('t0',0)])
        self.model['formula'] ='{A}*(1-exp(-{zeta}*{omega}*(x-{x0}))/sqrt(1-{zeta}**2) * cos({omega}*sqrt(1-{zeta}**2)-arctan({zeta}/sqrt(1-{zeta}**2)))) +{B}'
        # --- Guess Initial values
        x, y = self.data['x'],self.data['y']
        # TODO use signal
        omega0   = main_frequency(x,y) 
        A        = np.max(y) - np.min(y)
        B        = y[0]
        zeta     = 0.1
        y_start  = y[0]+0.01*A
        bDeviate = np.argwhere(abs(y-y_start)>abs(y_start-y[0]))[0]
        t0       = x[bDeviate[0]]
        p_guess  = np.array([A, omega0, zeta, B, t0])
        self.model['coeffs_init'] = p_guess
        # --- Set Bounds
        T  = x[-1]-x[0]
        dt = x[1]-x[0]
        om_min = 2*np.pi/T/2
        om_max = 2*np.pi/dt/2
        b_A    = (A*0.1,A*3)
        b_om   = (om_min,om_max)
        b_zeta = (0,1)
        b_B    = (np.min(y),np.max(y))
        b_x0   = (np.min(x),np.max(x))
        self.model['bounds'] = ((b_A[0],b_om[0],b_zeta[0],b_B[0],b_x0[0]),(b_A[1],b_om[1],b_zeta[1],b_B[1],b_x0[1]))
        #self.plot_guess(); import matplotlib.pyplot as plt; plt.show()
        #self.print_guessbounds(); 

# --------------------------------------------------------------------------------}
# --- Predefined fitter  
# --------------------------------------------------------------------------------{
class ContinuousPolynomialFitter(ModelFitter):
    def __init__(self,order=None, x=None, y=None, p0=None, bounds=None):
        ModelFitter.__init__(self,x=None, y=None, p0=p0, bounds=bounds)
        self.setOrder(int(order))
        if order is not None and x is not None and y is not None:
            self.fit_data(x,y,p0,bounds)

    def setOrder(self, order):
        self.order=order
        if order is not None:
            variables= string.ascii_lowercase[:order+1]
            self.model['coeffs']  = OrderedDict([(var,1) for i,var in enumerate(variables)])
            formula  = ' + '.join(['{}*x**{}'.format('{'+var+'}',order-i) for i,var in enumerate(variables)])
            self.model['formula']  = _clean_formula(formula)

    def fit_data(self, x, y, p0=None, bounds=None):
        if self.order is None:
            raise Exception('Polynomial Fitter not set, call function `setOrder` to set order')
        # Cleaning data
        x,y=self.clean_data(x,y)

        nParams=self.order+1
        # Bounds
        self.setup_bounds(bounds, nParams) # TODO
        # Initial conditions
        self.setup_guess(p0, bounds, nParams) # TODO

        # Fitting
        pfit  = np.polyfit(x,y,self.order)

        # --- Reporting information about the fit (after the fit)
        y_fit = np.polyval(pfit,x)
        self.store_fit_info(y_fit, pfit)

        # --- Return a fitted function
        self.model['fitted_function']=lambda xx : np.polyval(pfit,xx)


class DiscretePolynomialFitter(ModelFitter):
    def __init__(self,exponents=None, x=None, y=None, p0=None, bounds=None):
        ModelFitter.__init__(self,x=None, y=None, p0=p0, bounds=bounds)
        self.setExponents(exponents)
        if exponents is not None and x is not None and y is not None:
            self.fit_data(x,y,p0,bounds)

    def setExponents(self, exponents):
        self.exponents=exponents
        if exponents is not None:
            #exponents=-np.sort(-np.asarray(exponents))
            self.exponents=exponents
            variables= string.ascii_lowercase[:len(exponents)]
            self.model['coeffs']  = OrderedDict([(var,1) for i,var in enumerate(variables)])
            formula  = ' + '.join(['{}*x**{}'.format('{'+var+'}',e) for var,e in zip(variables,exponents)])
            self.model['formula']  = _clean_formula(formula)

    def fit_data(self, x, y, p0=None, bounds=None):
        if self.exponents is None:
            raise Exception('Polynomial Fitter not set, call function `setExponents` to set exponents')
        # Cleaning data, and store it in object
        x,y=self.clean_data(x,y)

        nParams=len(self.exponents)
        # Bounds
        self.setup_bounds(bounds, nParams) # TODO
        # Initial conditions
        self.setup_guess(p0, bounds, nParams) # TODO

        X_poly=np.array([])
        for i,e in enumerate(self.exponents):
            if i==0:
                X_poly = np.array([x**e])
            else:
                X_poly = np.vstack((X_poly,x**e))
        try:
            pfit = np.linalg.lstsq(X_poly.T, y, rcond=None)[0]
        except:
            pfit = np.linalg.lstsq(X_poly.T, y)

        # --- Reporting information about the fit (after the fit)
        y_fit= np.dot(pfit, X_poly)
        self.store_fit_info(y_fit, pfit)

        # --- Return a fitted function
        def fitted_function(xx):
            y=np.zeros(xx.shape)
            for i,(e,c) in enumerate(zip(self.exponents,pfit)):
                y += c*xx**e
            return y
        self.model['fitted_function']=fitted_function


class SinusoidFitter(ModelFitter):
    def __init__(self, physical=False, x=None, y=None, p0=None, bounds=None):
        ModelFitter.__init__(self, x=None, y=None, p0=p0, bounds=bounds)
        #self.setOrder(int(order))
        self.physical=physical
        if physical:
            self.model['coeffs']  = OrderedDict([('A',1),('f',1),('phi',0),('B',0)])
            self.model['formula']  = '{A} * sin(2*pi*({f}*x + {phi}/360)) + {B}'
        else:
            self.model['coeffs']  = OrderedDict([('A',1),('omega',1),('phi',0),('B',0)])
            self.model['formula']  = '{A} * sin({omega}*x + {phi}) + {B}'

        if x is not None and y is not None:
            self.fit_data(x,y,p0,bounds)

    def fit_data(self, x, y, p0=None, bounds=None):
        # Cleaning data
        x,y=self.clean_data(x,y)

        # TODO use signal
        guess_freq= main_frequency(x,y)/(2*np.pi) # [Hz]
        guess_amp = np.std(y) * 2.**0.5
        guess_offset = np.mean(y)
        if self.physical:
            guess = np.array([guess_amp, guess_freq, 0., guess_offset])
            minimize_me = lambda x, *p : sinusoid_f(x, p)
        else:
            guess = np.array([guess_amp, 2.*np.pi*guess_freq, 0., guess_offset])
            minimize_me = lambda x, *p : sinusoid(x, p)
        self.model['coeffs_init'] = guess

        pfit, pcov = so.curve_fit(minimize_me, x, y, p0=guess)

        # --- Reporting information about the fit (after the fit)
        # And Return a fitted function
        if self.physical:
            y_fit = sinusoid_f(x,  pfit)
            self.model['fitted_function']=lambda xx : sinusoid_f(xx, pfit)
        else:
            y_fit = sinusoid(x,  pfit)
            self.model['fitted_function']=lambda xx : sinusoid(xx, pfit)
        self.store_fit_info(y_fit, pfit)



class GeneratorTorqueFitter(ModelFitter):
    def __init__(self,x=None, y=None, p0=None, bounds=None):
        ModelFitter.__init__(self,x=None, y=None, p0=p0, bounds=bounds)

#         RtGnSp, RtTq  , Rgn2K , SlPc , SpdGenOn = p
#         {'label':'Generator Torque', 'handle': gentorque, 'id':'predef: gentorque',
#         'formula': '{RtGnSp} , {RtTq}  , {Rgn2K} , {SlPc} , {SpdGenOn}',
        self.model['coeffs']= extract_key_num('RtGnSp=100 , RtTq=1000  , Rgn2K=0.01 ,SlPc=5 , SpdGenOn=0')
#         'consts' :None,
#         'bounds' :'RtGnSp=(0.1,inf) , RtTq=(1,inf), Rgn2K=(0.0,0.1) ,SlPc=(0,20) , SpdGenOn=(0,inf)'}
        if x is not None and y is not None:
            self.fit_data(x,y,p0,bounds)

    def fit_data(self, x, y, p0=None, bounds=None):
        #nParams=5
        ## Bounds
        #self.setup_bounds(bounds,nParams) # TODO
        ## Initial conditions
        #self.setup_guess(p0,bounds,nParams) # TODO

        # Cleaning data, and store it in object
        x,y=self.clean_data(x,y)

        I = np.argsort(x)
        x=x[I]
        y=y[I]

        # Estimating deltas
        xMin, xMax=np.min(x),np.max(x)
        yMin, yMax=np.min(y),np.max(y)
        DeltaX = (xMax-xMin)*0.02
        DeltaY = (yMax-yMin)*0.02

        # Binning data
        x_bin=np.linspace(xMin,xMax,min(200,len(x)))
        x_lin=x_bin[0:-1]+np.diff(x_bin)
        #y_lin=np.interp(x_lin,x,y) # TODO replace by bining
        y_lin = np.histogram(y, x_bin, weights=y)[0]/ np.histogram(y, x_bin)[0]
        y_lin, _, _ = stats.binned_statistic(x, y, statistic='mean', bins=x_bin)
        x_lin, _, _ = stats.binned_statistic(x, x, statistic='mean', bins=x_bin)
        bNaN=~np.isnan(y_lin)
        y_lin=y_lin[bNaN]
        x_lin=x_lin[bNaN]

        # --- Find good guess of parameters based on data
        # SpdGenOn
        iOn = np.where(y>0)[0][0]
        SpdGenOn_0    =  x[iOn]
        SpdGenOn_Bnds = (max(x[iOn]-DeltaX,xMin), min(x[iOn]+DeltaX,xMax))
        # Slpc
        Slpc_0    = 5
        Slpc_Bnds = (0,10)
        # RtTq
        RtTq_0    = yMax
        RtTq_Bnds = (yMax-DeltaY, yMax+DeltaY)
        # RtGnSp
        iCloseRt = np.where(y>yMax*0.50)[0][0]
        RtGnSp_0    = x[iCloseRt]
        RtGnSp_Bnds = ( RtGnSp_0 -DeltaX*2, RtGnSp_0+DeltaX*2)
        # Rgn2K
        #print('>>>',SpdGenOn_0, RtGnSp_0)
        bR2=np.logical_and(x>SpdGenOn_0, x<RtGnSp_0)
        exponents=[2]
        _, pfit, _ = fit_polynomial_discrete(x[bR2], y[bR2], exponents)
        #print(pfit)
        Rgn2K_0   =pfit[0]
        Rgn2K_Bnds=(pfit[0]/2, pfit[0]*2)
#         import matplotlib.pyplot as plt
#         fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
#         fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
#         ax.plot(x,y ,'-'   , label='')
#         ax.plot(x[bR2],y[bR2],'ko', label='')
#         ax.plot(x_lin,y_lin,'bd', label='')
#         ax.set_xlabel('')
#         ax.set_ylabel('')
#         ax.tick_params(direction='in')
#         plt.show()
        def minimize_me(p):
            RtGnSp, RtTq  , Rgn2K , SlPc , SpdGenOn = p
            y_model=np.array([gentorque(x_lin, (RtGnSp, RtTq  , Rgn2K , SlPc , SpdGenOn))])
            eps = np.mean((y_lin-y_model)**2)
#             print(eps,p)
            return  eps
        bounds = (RtGnSp_Bnds, RtTq_Bnds, Rgn2K_Bnds, Slpc_Bnds, SpdGenOn_Bnds)
        p0     = [RtGnSp_0, RtTq_0, Rgn2K_0, Slpc_0, SpdGenOn_0]
        #print('Bounds',bounds)
        #print('p0',p0)
        res = so.minimize(minimize_me, x0=p0, bounds= bounds, method='SLSQP')
        pfit=res.x

        # --- Reporting information about the fit (after the fit)
        y_fit= gentorque(x, pfit)
        self.store_fit_info(y_fit, pfit)
        # --- Return a fitted function
        self.model['fitted_function']=lambda x: gentorque(x,pfit)



# --- Registering FITTERS. The formula info is redundant, only used by pyDatView, should be removed
FITTERS= [
{'label':'Polynomial (full)'   ,'id':'fitter: polynomial_continuous', 'handle': ContinuousPolynomialFitter,
'consts':{'order':3}, 'formula': '{a_i} x^i'},
{'label':'Polynomial (partial)','id':'fitter: polynomial_discrete'  , 'handle': DiscretePolynomialFitter  ,
'consts':{'exponents':[0,2,3]},'formula': '{a_i} x^j'},
{'label':'Sinusoid','id':'fitter: sinusoid'  , 'handle': SinusoidFitter  ,
'consts':{'physical':True},'formula': '{A}*sin({omega or 2 pi f}*x+{phi or phi_deg}) + {B} '},
{'label':'2nd order impulse/decay (auto)','id':'fitter: secondorder_impulse', 'handle': SecondOrderFitterImpulse  ,
'consts':{},'formula': '{A}*exp(-{zeta}*{omega}*(x-{x0}))*sin({omega}*sqrt(1-{zeta}**2)))+{B}'},
{'label':'2nd order step (auto)','id':'fitter: secondorder_step', 'handle': SecondOrderFitterStep  ,
'consts':{},'formula':'{A}*(1-exp(-{zeta}*{omega}*(x-{x0}))/sqrt(1-{zeta}**2)*cos({omega}*sqrt(1-{zeta}**2)-arctan({zeta}/sqrt(1-{zeta}**2))))+{B}'},
# {'label':'Generator Torque','id':'fitter: gentorque'  , 'handle': GeneratorTorqueFitter  ,
# 'consts':{},'formula': ''}
]

# --------------------------------------------------------------------------------}
# --- Helper functions  
# --------------------------------------------------------------------------------{
def extract_variables(sFormula):
    """ Extract variables in expression, e.g.  {a}*x + {b} -> ['a','b']
    The variables are replaced with p[0],..,p[n] in order of appearance
    """
    regex = r"\{(.*?)\}"
    matches = re.finditer(regex, sFormula, re.DOTALL)
    formula_eval=sFormula
    variables=[]
    ivar=0
    for i, match in enumerate(matches):
        for groupNum in range(0, len(match.groups())):
            var = match.group(1)
            if var not in variables: 
                variables.append(var)
                formula_eval = formula_eval.replace('{'+match.group(1)+'}','p[{:d}]'.format(ivar))
                ivar+=1
    return variables, formula_eval


def extract_key_tuples(text):
    """
    all=(0.1,-2),b=(inf,0), c=(-inf,0.3e+10)
    """
    if text is None:
        return {}
    regex = re.compile(r'(?P<key>[\w\-]+)=\((?P<value1>[0-9+epinf.-]*?),(?P<value2>[0-9+epinf.-]*?)\)($|,)')
    return  {match.group("key"): (float(match.group("value1")),float(match.group("value2"))) for match in regex.finditer(text.replace(' ',''))}

def extract_key_num(text):
    """
    all=0.1, b=inf, c=-0.3e+10
    """
    if text is None:
        return {}
    regex = re.compile(r'(?P<key>[\w\-]+)=(?P<value>[0-9+epinf.-]*?)($|,)')
    return OrderedDict([(match.group("key"), float(match.group("value"))) for match in regex.finditer(text.replace(' ',''))])

def extract_key_miscnum(text):
    """
    all=0.1, b=(inf,0), c=[-inf,0.3e+10,10,11])
    """
    def isint(s):
        try:
            int(s)
            return True
        except:
            return False

    if text is None:
        return {}
    sp=re.compile('([\w]+)=').split(text.replace(' ',''))
    if len(sp)<3:
        return {}
    sp=sp[1:]
    keys   = sp[0::2]
    values = sp[1::2]
    d={}
    for (k,v) in zip(keys,values):
        if v.find('(')>=0:
            v=v.replace('(','').replace(')','')
            v=v.split(',')
            vect=tuple([float(val) for val in v if len(val.strip())>0])
        elif v.find('[')>=0:
            v=v.replace('[','').replace(']','')
            v=v.split(',')
            vect=[int(val) if isint(val) else float(val) for val in v if len(val.strip())>0] # NOTE returning lists
        elif v.find('True')>=0:
            v=v.replace(',','').strip()
            vect=True
        elif v.find('False')>=0:
            v=v.replace(',','').strip()
            vect=False
        else:
            v=v.replace(',','').strip()
            vect=int(v) if isint(v) else float(v)
        d[k]=vect
    return d

def set_common_keys(dict_target, dict_source):
    """ Set a dictionary using another one, missing keys in source dictionary are reported"""
    keys_missing=[]
    for k in dict_target.keys():
        if k in dict_source.keys():
            dict_target[k]=dict_source[k]
        else:
            keys_missing.append(k)
    return dict_target, keys_missing

def _clean_formula(s, latex=False):
    s = s.replace('+-','-').replace('**1','').replace('*x**0','')
    s = s.replace('np.','')
    if latex:
        #s = s.replace('{','$').replace('}','$')
        s = s.replace('phi',r'\phi')
        s = s.replace('alpha',r'\alpha')
        s = s.replace('beta' ,r'\alpha')
        s = s.replace('zeta' ,r'\zeta')
        s = s.replace('mu'   ,r'\mu'   )
        s = s.replace('pi'   ,r'\pi'   )
        s = s.replace('sigma',r'\sigma')
        s = s.replace('omega',r'\omega')
        s = s.replace('_ref',r'_{ref}') # make this general
        s = s.replace(r'(',r'{(')
        s = s.replace(r')',r')}')
        s = s.replace(r'**',r'^')
        s = s.replace(r'*', '')
        s = s.replace('sin',r'\sin')
        s = s.replace('exp',r'\exp')
        s = s.replace('sqrt',r'\sqrt')
        s = r'$'+s+r'$'
    else:
        s = s.replace('{','').replace('}','')
    return s


def main_frequency(t,y):
    """ 
    Returns main frequency of a signal
    NOTE: this tool below to welib.tools.signal_analysis, but put here for convenience
    """
    dt       = t[1]-t[0]  # assume uniform spacing of time and frequency
    om       = np.fft.fftfreq(len(t), (dt))*2*np.pi 
    Fyy      = abs(np.fft.fft(y))
    omega    = abs(om[np.argmax(Fyy[1:])+1]) # exclude the zero frequency (mean)
    return omega

def rsquare(y, f): 
    """ Compute coefficient of determination of data fit model and RMSE
    [r2] = rsquare(y,f)
    RSQUARE computes the coefficient of determination (R-square) value from
    actual data Y and model data F. 
    INPUTS
      y       : Actual data
      f       : Model fit
    OUTPUT
      R2      : Coefficient of determination
    """
    # Compare inputs
    if not np.all(y.shape == f.shape) :
        raise Exception('Y and F must be the same size')
    # Check for NaN
    tmp = np.logical_not(np.logical_or(np.isnan(y),np.isnan(f))) 
    y = y[tmp]
    f = f[tmp]
    R2 = max(0,1-np.sum((y-f)**2)/np.sum((y-np.mean(y))** 2))
    return R2

def pretty_param(s):
    if s in ['alpha','beta','delta','gamma','epsilon','zeta','lambda','mu','nu','pi','rho','sigma','phi','psi','omega']:
        s = r'\{}'.format(s)
    s = s.replace('_ref',r'_{ref}') # make this general..
    return s

def pretty_num(x):
    if abs(x)<1000 and abs(x)>1e-4:
        return "{:9.4f}".format(x)
    else:
        return '{:.3e}'.format(x)

def pretty_num_short(x,digits=3):
    if digits==4:
        if abs(x)<1000 and abs(x)>1e-1:
            return "{:.4f}".format(x)
        else:
           return "{:.4e}".format(x)
    elif digits==3:
        if abs(x)<1000 and abs(x)>1e-1:
            return "{:.3f}".format(x)
        else:
           return "{:.3e}".format(x)
    elif digits==2:
        if abs(x)<1000 and abs(x)>1e-1:
            return "{:.2f}".format(x)
        else:
           return "{:.2e}".format(x)


if __name__ == '__main__':
    # --- Writing example models to file for pyDatView tests
    a,b,c = 2.0, 3.0, 4.0
    u_ref,z_ref,alpha=10,12,0.12
    mu,sigma=0.5,1.2
    x = np.linspace(0.1,30,20)
    A,k,B=0.5,1.2,10
    y_exp=expdecay(x,(A,k,B))
    A, k = 10, 2.3,
    y_weib=weibull_pdf(x,(A,k))
    y_log=logarithmic(x,(a,b))
    exponents=[0,3,5]
    y_poly = a + b*x**3 + c*x**5
    y_power=powerlaw_all(x,(alpha,u_ref,z_ref))
    y_gauss=gaussian(x,(mu,sigma))
    A= 101; B= -200.5; omega = 0.4; phi = np.pi/3
    y_sin=sinusoid(x,(A,omega,phi,B)) + np.random.normal(0, 0.1, len(x))
    M=np.column_stack((x,y_poly,y_power,y_gauss,y_gauss+10,y_weib,y_exp,y_log,y_sin))
    np.savetxt('../TestFit.csv',M,header='x,poly,power,gauss,gauss_off,weib,expdecay,log,sin',delimiter=',')
