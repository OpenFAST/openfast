import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
from pyFAST.airfoils.Polar import Polar

MyDir=os.path.dirname(__file__)


def main_correction3D(test=False):
    polarFile_in = os.path.join(MyDir,'../data/DU21_A17.csv')

    r_over_R     = 0.2   # spanwise location [-]
    chord_over_r = 3./5. # chord divided by local radius [-]
    tsr          = 10    # tip speed ratio [-]

    polar = Polar(polarFile_in, compute_params=True, verbose=False)
    #ADpol = polar.toAeroDyn(polarFile_AD) # Optional, write to AeroDyn format
    polar3D= polar.correction3D(r_over_R, chord_over_r, tsr)

    # --- Plot
    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    ax.plot(polar.alpha, polar.cl     ,'k-'  , label= r'2D polar')
    ax.plot(polar3D.alpha, polar3D.cl , '-'  , label= r'3D corrected')
    ax.plot(polar.alpha, polar.cl_inv ,'k--' , label= r'inviscid')
    ax.tick_params(direction='in', top=True, right=True)
    ax.set_xlabel(r'Angle of attack, $\alpha$ [deg]')
    ax.set_ylabel(r'Lift coefficient, $C_l$ [-]')
    ax.set_title(r'Airfoils - 3D correction')
    ax.set_xlim([-50,50])
    ax.set_ylim([-1.5,2])
    ax.legend()

    return polar, polar3D

polar,polar3D = main_correction3D()

if __name__ == '__main__':
    plt.show()
if __name__ == '__test__':
    pass
