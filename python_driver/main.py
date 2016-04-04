#! /usr/bin/env python
# -*- coding: utf-8 -*-


#   Copyright (C) 2014 mdm                                     
#   map[dot]plus[dot]plus[dot]help[at]gmail                     
#                                                              
# Licensed to the Apache Software Foundation (ASF) under one   
# or more contributor license agreements.  See the NOTICE file 
# distributed with this work for additional information        
# regarding copyright ownership.  The ASF licenses this file   
# to you under the Apache License, Version 2.0 (the            
# "License"); you may not use this file except in compliance   
# with the License.  You may obtain a copy of the License at   
#                                                              
#   http://www.apache.org/licenses/LICENSE-2.0                 
#                                                              
# Unless required by applicable law or agreed to in writing,   
# software distributed under the License is distributed on an  
# "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY       
# KIND, either express or implied.  See the License for the    
# specific language governing permissions and limitations            
# under the License.                                             


if __name__ == '__main__':      
    from mapsys import *
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    import numpy as np
    np.set_printoptions(formatter={'float': '{: 13.1f}'.format},linewidth=100)

    mooring_1 = Map( )
    
    mooring_1.map_set_sea_depth(320)
    mooring_1.map_set_gravity(9.81)
    mooring_1.map_set_sea_density(1025.0)
    
    # mooring_1.read_file("input.map")                                         # 350 m depth        
    # mooring_1.read_file("inwind_new.map")
    # mooring_1.read_file("../test/test_for_nrel.map")
    # mooring_1.read_file("../test/baseline_1.map")                                        # 120 m depth    
    # mooring_1.read_file("../test/baseline_2.map")                                       # 350 m depth    
    # mooring_1.read_file("../test/baseline_5.map")                                       # 80 m depth
    # mooring_1.read_file("../test/NRELOffshrBsline5MW_Platform_OC3Hywind.map")           # 320 m depth
    # mooring_1.read_file("../test/NRELOffshrBsline5MW_Platform_OC3Hywind_segmented.map")   # 320 m depth
    mooring_1.read_file("../test/NRELOffshrBsLine5MW_OC4.map")                            # 200 m depth
    # mooring_1.read_file("../test/NRELOffshrBsLine5MW_TLP.map")                          # 200 m depth

    # mooring_1.summary_file('name_me.txt')
    mooring_1.init( )

    epsilon = 1e-5
    K = mooring_1.linear(epsilon)    
    print "\nHere is the linearized stiffness matrix with zero vessel displacement:"
    print np.array(K)

    #mooring_1.displace_vessel(5,0,0,0,0,0)
    #mooring_1.update_states(0.0,0)

    #mooring_1.displace_vessel(17,0,0,0,0,0)
    #mooring_1.update_states(0.0,0)

    # H,V = mooring_1.get_fairlead_force_2d(0)    
    # print H, "  ", V
      
    # fx,fy,fz = mooring_1.get_fairlead_force_3d(0)    
    # print fx, "  ", fy, "  ", fz
    # 
    # ''' 
    # function residual at (hopefully) the solution
    # '''
    # 
    # print mooring_1.funch(0) 
    # print mooring_1.funcl(0)
    # 
    # '''
    # derivatives at solution
    # '''
    # print mooring_1.dxdh(0)
    # print mooring_1.dxdv(0)    
    # print mooring_1.dzdh(0)
    # print mooring_1.dzdv(0)
    # 
    # print mooring_1.dxdh(1)
    # print mooring_1.dxdv(1)    
    # print mooring_1.dzdh(1)
    # print mooring_1.dzdv(1)

    fig = plt.figure()
    ax = Axes3D(fig)
    for i in range(0,mooring_1.size_lines()):
        x = mooring_1.plot_x( i, 20 )
        y = mooring_1.plot_y( i, 20 )
        z = mooring_1.plot_z( i, 20 )        
        ax.plot(x,y,z,'b-')
     
    ax.set_xlabel('X [m]')
    ax.set_ylabel('Y [m]')
    ax.set_zlabel('Z [m]')        
    #ax.set_xlim([-3.0,3])        
    #ax.set_ylim([-3.0,3])        
    #ax.set_zlim([-3.0,0])        
     
    plt.show()
    
    mooring_1.end( )
