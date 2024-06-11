""" 
Example to generate an AeroDyn polar file from a set of Cl-Cd data
 - The various parameters (e.g. unsteady parameters) are computed and updated
 - The AD file is written
"""
import numpy as np
import os

from pyFAST.airfoils.Polar import Polar
from pyFAST.input_output.fast_input_file import FASTInputFile
from pyFAST.input_output.csv_file import CSVFile

# Get current directory so this script can be called from any location
scriptDir=os.path.dirname(__file__)


def main_ReWriteADFile():
    """ 
    Example 1: 
      - open an existing AeroDyn polar file
      - rewrite it (unsteady parameters are recomputed)
    """
    AD_polarFile_in = os.path.join(scriptDir,'../../../data/NREL5MW/5MW_Baseline/Airfoils/Cylinder1.dat')
    AD_polarFile_out = '_Polar_out.dat.ignore'

    # Open an existing AeroDyn polar file
    polar = Polar(AD_polarFile_in, fformat='ADPolar')
    # Rewrite it (unsteady parameters are recomputed)
    # NOTE: you can provide templateFile=AD_polarFile_in to the function below
    #       to ensure that the file will look the same
    comment = 'Cylinder at Re=6Million\nUpdated unsteady parameters' # Optional comment
    Re      = 6 # Reynolds number in Million (Optional)
    ADpol = polar.toAeroDyn(AD_polarFile_out, Re=6, comment=comment)

    return ADpol, polar


def main_WriteADPolar():
    """ 
    Example 2: 
      - Open a tabulated file with alpha,Cl,Cd,Cm 
      - Write an AeroDyn file from it (unsteady parameters are computed)
    """
    polarFile_in     = os.path.join(scriptDir,'../data/DU21_A17.csv')
    polarFile_AD_out = '_Polar_out.dat.ignore'

    # Open a tabulated file with alpha,Cl,Cd,Cm 
    polar = Polar(polarFile_in, fformat='delimited')
    # Write an AeroDyn file from it (unsteady parameters are computed)
    # NOTE: you can provide templateFile='ADTemplate.dat' to the function below
    #       to ensure that the AD file will look the same as the template.
    ADpol = polar.toAeroDyn(polarFile_AD_out)

    return ADpol, polar

def main_WriteADPolarLowLevel():
    """
    Example 3: Same as Example 2, but with low level interface.
    """
    # --- Reading an existing AD file, just as a template, we'll replace things in it
    templateADFile = os.path.join(scriptDir,'../../../data/NREL5MW/5MW_Baseline/Airfoils/Cylinder1.dat')
    ADpol = FASTInputFile(templateADFile)

    # --- Creating a Polar object from Cl-Cd data
    polarFile = os.path.join(scriptDir,'../data/DU21_A17.csv')
    p=CSVFile(polarFile).toDataFrame().values
    polar= Polar(alpha=p[:,0],cl=p[:,1],cd=p[:,2],cm=p[:,3])
    (alpha0,alpha1,alpha2,cnSlope,cn1,cn2,cd0,cm0)=polar.unsteadyParams()

    # --- Updating the AD polar file 
    # Setting unsteady parameters
    ADpol['Re'] = 1.0000 # TODO UNKNOWN
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
    PolarTable = np.column_stack((polar.alpha,polar.cl,polar.cd,polar.cm))
    ADpol['NumAlf'] = polar.cl.shape[0]
    ADpol['AFCoeff'] = np.around(PolarTable, 5)

    filename='_Polar_out.dat.ignore'
    ADpol.write(filename)
    #print('Writing polar to file:',filename,' thick={}'.format(t))

    return ADpol, polar



ADpol,polar = main_ReWriteADFile()
ADpol,polar = main_WriteADPolar()
ADpol,polar = main_WriteADPolarLowLevel()

if __name__ == '__main__':
    print(polar)
    import matplotlib.pyplot as plt
    plt.plot(polar.alpha   ,polar.cl     , '-' , label= 'cl')
    plt.plot(polar.alpha   ,polar.cl_lin , '--', label= 'cl_lin')
    plt.ylim([-2,2])
    plt.legend()
    plt.show()

if __name__ == '__test__':
    ADpol,polar = main_ReWriteADFile()
    ADpol,polar = main_WriteADPolar()
    ADpol,polar = main_WriteADPolar()
    try:
        os.remove('_Polar_out.dat.ignore')
    except:
        pass
