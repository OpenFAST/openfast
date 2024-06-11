"""
 - Create an AeroDyn polar from a CSV file
 - (optional) Read the AeroDyn polar back 
 - Plot unsteady aerodynamic parameters
"""
import os
import matplotlib.pyplot as plt

# Get current directory so this script can be called from any location
scriptDir=os.path.dirname(__file__)

# --- Create an AeroDyn polar from a CSV file
from pyFAST.airfoils.Polar import Polar
polarFile_in = os.path.join(scriptDir,'../data/DU21_A17.csv')
polar = Polar(polarFile_in, fformat='delimited')
ADpol = polar.toAeroDyn('_AeroDyn_Polar_DU21_A17.dat')


fig = plt.figure()
plt.plot(polar.alpha, polar.cl, label='Cl')
plt.plot(polar.alpha, polar.cl_lin, label='Cl_lin')
plt.xlabel(r'Angle of attack $\alpha$ [deg]')
plt.ylabel(r'Lift coefficient $C_l$ [-]')
plt.ylim([-1.5,1.9])
plt.legend()

# --- Optional: read the AeroDyn polar
from pyFAST.input_output import FASTInputFile
ADpol = FASTInputFile('_AeroDyn_Polar_DU21_A17.dat')

# --- Plot important data for unsteady aerodynamics
df = ADpol.toDataFrame()
print(df.keys())
fig = plt.figure()
plt.plot(df['Alpha_[deg]'], df['Cn_[-]']    , label='Cn')
plt.plot(df['Alpha_[deg]'], df['Cn_pot_[-]'], label='Cn_lin')
plt.plot(df['Alpha_[deg]'], df['Cn_012_[-]'],'--', label='Cn_2,0,1')
plt.xlabel(r'Angle of attack $\alpha$ [deg]')
plt.ylabel(r'Normal coefficient $C_n$ [-]')
plt.ylim([-1.5,1.9])
plt.legend()


# --- Alternative using class Polar (latest python toolbox)
try:
    fig = plt.figure()
    plt.plot(polar.alpha, polar.cn    , label='Cn')
    plt.plot(polar.alpha, polar.cn_lin, label='Cn_lin')
    plt.plot(polar.alpha_201, polar.cn_201 ,'o', label='Cn_2,0,1')
    plt.xlabel(r'Angle of attack $\alpha$ [deg]')
    plt.ylabel(r'Normal coefficient $C_n$ [-]')
    plt.ylim([-1.5,1.9])
    plt.legend()
except:
    pass


if __name__ == '__main__':
    plt.show()

if __name__ == '__test__':
    try:
        os.remove('_AeroDyn_Polar_DU21_A17.dat')
    except:
        pass
