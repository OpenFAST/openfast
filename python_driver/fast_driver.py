#! /usr/bin/env python
# -*- coding: utf-8 -*-

'''
  Copyright (C) 2014 mdm                                     
  http://www.apache.org/licenses/LICENSE-2.0                 

  use like this:
  $ ./fast_driver.py -f ../test/Test22.out -n T[5]
'''  

from fast_driver_support import *


def get_line_tension(mooring, vessel, i):
    fairlead_number = 1
    mooring.displace_vessel(vessel.x[i], vessel.y[i], vessel.z[i], vessel.phi[i], vessel.the[i], vessel.psi[i])
    mooring.update_states(vessel.time[i], 0)
    fx,fy,fz = mooring.get_fairlead_force_3d(fairlead_number)
    return np.sqrt(fx**2 + fy**2 + fz**2)


if __name__ == '__main__':      
    # command line argument processing
    usage = 'usage: %prog [options]'
    parser = OptionParser(usage)
    parser.add_option('-f', '--file', dest='file_name', help='read data from FILENAME')
    parser.add_option('-n', '--output stream', dest='chanel', help='read output stream from file')
    (options, args) = parser.parse_args()

    # index is the surge, ... , yaw column number in the FAST output file
    index = get_vessel_column_index(options.file_name,options.chanel) 

    mooring = Map()
    mooring.map_set_sea_depth(220) # 150 for barge
    mooring.map_set_gravity(9.81) 
    mooring.map_set_sea_density(1020.0)
    mooring.read_file('Mooring.map') # 100 m depth
    mooring.summary_file('barge.sum.txt')
    mooring.init()

    # read column data in FAST output file
    with open(options.file_name, 'rb') as csv_file:        
        reader = csv.reader(csv_file,delimiter=' ',skipinitialspace=True)
        table = list(reader)


    # Set vessel displacement time series based on FAST output file
    vessel = set_vessel_prescribed_motion(table,index)

    # Now displace the vessel for each time step and calculate the fairlead tension.
    # This is looped over the entire vessel-displacement time series. Tension is an
    # array. 
    tension = [get_line_tension(mooring,vessel,i) for i in range(0,len(vessel.time))]
    
    mooring.end()
    
    plt.figure(1)
    plt.title('Tension T[1]')
    plt.xlabel('Time [sec]')
    plt.ylabel('Force [kN]')
    plt.plot(vessel.time,vessel.out,'b',lw=3,alpha=0.5,label='FAST Output')
    plt.plot(vessel.time,tension,'k',label='MAP++ Python Driver')
    plt.legend(loc=2)

    plt.figure(2)
#    plt.title('Tension T[5]')
    plt.xlabel('Time [sec]')
    plt.ylabel('Displacement [m]')
    plt.plot(vessel.time,vessel.x,'r',label='surge')
    plt.plot(vessel.time,vessel.x,'g',label='sway')
    plt.plot(vessel.time,vessel.z,'b',label='heave')
    plt.legend(loc=2)

    plt.show()    
