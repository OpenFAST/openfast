#! /usr/bin/env python
# -*- coding: utf-8 -*-


from mapsys import *
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.ticker as mtick
from matplotlib import rcParams
import numpy as np 
np.set_printoptions(precision=2)
np.set_printoptions(suppress=True)
rcParams.update({'figure.autolayout': True})


# user function to plot the mooring profile and footprint
def plot_mooring_system(mooring_data):
    # plot the mooring profile
    fig = plt.figure(1)
    ax = Axes3D(fig)
    colors = ['b','g','r','c']
    for i in xrange(mooring_data.size_lines()):
        x = mooring_data.plot_x( i, 20 ) # i is the the line number, and 20 is the number of points plotted on the line 
        y = mooring_data.plot_y( i, 20)
        z = mooring_data.plot_z( i, 20)        
        ax.plot(x,y,z,colors[i]+'-')     
    ax.set_xlabel('X [m]')
    ax.set_ylabel('Y [m]')
    ax.set_zlabel('Z [m]')        


def start():
    """
    Step 1) First initialize an instance of a mooring system

    Step 2) Assume that (X, Y, Z, phi, theta, psi) are the translation and rotation displacement of the vessel. 
    These displacements are fed into MAP as an argument to displace the fairlead. With the fairlead(s) 
    rigidly connected to the vessel, the (X, Y, Z, phi, theta, psi) directly manifests into the fairlead
    position in the global frame. 
    
    For the time being, assume a sinusoidal displacement of the vessel

    Step 3) This for-loop emulates the time-stepping procedure. You want to loop through the length of 
    the arrays (X,Y,Z,phi,theta,psi) to retrieve the fairlead force

    Step 4) update the MAP state. The arguments in displace_vessel are the displace displacements and rotations about the reference origin. 
    In this case, the reference origin is (0,0,0). 
    They can be set to a different potision using a run-time argument (this is an advanced feature).

    Step 5) get the fairlead tension. The get_fairlead_force_3d returns the fairlead force in 
    X, Y Z coordinates. This must be called at each time-step, and then stored into an array. We append 
    the empty lists created on lines 84-88.

    .. Note::

       MAP does *NOT* return the mooring restoring moment, The user must calculate this 
       themself using the cross-product between the WEC reference origin and the mooring attached
       point, i.e.,
    
       :math:`\mathbf{Moment} = \mathbf{r} \times \mathbf{F}`
    """

    # Step 1
    mooring = Map() 
    mooring.map_set_sea_depth(120)        # m
    mooring.map_set_gravity(9.81)         # m/s^2
    mooring.map_set_sea_density(1020.0)   # kg/m^2
    mooring.read_file('../test/baseline_1.map')   # input file
    mooring.summary_file('summary_file.sum.txt') # output summary file name at the conclusion of initialization
    
    mooring.init()                # solve the cable equilibrium profile
    plot_mooring_system(mooring) # Optional: call the user function to illustrate the mooring equilibrium profile
    

    # initialize list to zero (this is artificial. This would be prescribed the by vessel program)
    X,Y,Z,phi,theta,psi = ([0.0 for i in xrange(500)] for _ in xrange(6))
    time = []

    # variable to specify the amplitude of surge oscillation and period factor
    dt = 0.1
    amplitude = 10.0

    # prescribe artificial surge and pitch displacement. Again, this should be supplied based on the WEC motion or from time-marching routine
    for i in xrange(len(X)):
        time.append(i*dt)
        X[i] = (amplitude)*(math.sin(i*0.05))
        theta[i] = (amplitude)*(math.sin(i*0.025))

    # Optional: plot the vessel displacement (surge=X and pitch=theta) as a function of time
    plt.figure(2)
    plt.plot(time,X,lw=2,label='Surge displacement')
    plt.plot(time,theta,lw=2,label='Pitch displacement')
    plt.title('Vessel Translation/Rotation')
    plt.ylabel('Amplitude [m,deg]')
    plt.xlabel('Time [sec]')
    plt.legend()

    # create an empty list of the line tension. We will store result from MAP in these lists
    line1_fx, line1_fy, line1_fz = ([] for _ in xrange(3))
    line2_fx, line2_fy, line2_fz = ([] for _ in xrange(3))
    line3_fx, line3_fy, line3_fz = ([] for _ in xrange(3))
    line4_fx, line4_fy, line4_fz = ([] for _ in xrange(3))

    # Step 3) 
    for i in xrange(len(X)):        
        # Step 4) 

        # displace the vessel, X,Y,X are in units of m, and phi, theta, psi are in units of degrees
        mooring.displace_vessel(X[i], Y[i], Z[i], phi[i], theta[i], psi[i]) 

        # first argument is the current time. Second argument is the coupling interval (used in FAST)
        mooring.update_states(time[i], 0)                                   

        # Step 5) 
        # line 1 tensions in X, Y and Z. Note that python is indexed started at zero
        fx, fy, fz = mooring.get_fairlead_force_3d(0) # arugment is the line number
        line1_fx.append(fx)
        line1_fy.append(fy)
        line1_fz.append(fz)        

        # line 2 tensions in X, Y and Z.
        fx, fy, fz = mooring.get_fairlead_force_3d(1)
        line2_fx.append(fx)
        line2_fy.append(fy)
        line2_fz.append(fz)        

        # line 3 tensions in X, Y and Z.
        fx, fy, fz = mooring.get_fairlead_force_3d(2)
        line3_fx.append(fx)
        line3_fy.append(fy)
        line3_fz.append(fz)        

        # line 4 tensions in X, Y and Z.
        fx, fy, fz = mooring.get_fairlead_force_3d(3)
        line4_fx.append(fx)
        line4_fy.append(fy)
        line4_fz.append(fz)        

    # Optional: plot line tension time history
    plt.figure(3)
    ax=plt.subplot(3,1,1)
    plt.plot(time,line1_fx,label='Line 1')
    plt.plot(time,line2_fx,label='Line 2')
    plt.plot(time,line3_fx,label='Line 3')
    plt.plot(time,line4_fx,label='Line 4')
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
    plt.ylabel('X Fairlead Force [N]')
    plt.legend()

    ax = plt.subplot(3,1,2)
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
    plt.plot(time,line1_fy)
    plt.plot(time,line2_fy)
    plt.plot(time,line3_fy)
    plt.plot(time,line4_fy)
    plt.ylabel('Y Fairlead Force [N]')

    ax = plt.subplot(3,1,3)
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
    plt.plot(time,line1_fz)
    plt.plot(time,line2_fz)
    plt.plot(time,line3_fz)
    plt.plot(time,line4_fz)
    plt.ylabel('Z Fairlead Force [N]')
    plt.xlabel('Time [sec]')

    plt.show()

if __name__ == '__main__':      
    start()
