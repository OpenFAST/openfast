import numpy as np
import os, shutil, platform

# Path to XFoil
path2xfoil      = 'Xfoil/bin/xfoil'

## inputs
if platform.system() == 'Windows':
    FAST_directory  = os.path.dirname( os.path.dirname( os.path.dirname( os.path.dirname( os.path.realpath(__file__) ) ) ) ) + os.sep + 'reg_tests' + os.sep + 'r-tests' + os.sep + 'glue-codes' + os.sep + 'openfast' + os.sep + 'IEA_LB_RWT-AeroAcoustics'
else:
    FAST_directory  = os.path.dirname( os.path.dirname( os.path.dirname( os.path.dirname( os.path.realpath(__file__) ) ) ) ) + os.sep + 'openfast' + os.sep + 'reg_tests' + os.sep + 'r-tests' + os.sep + 'glue-codes' + os.sep + 'openfast' + os.sep + 'IEA_LB_RWT-AeroAcoustics'
folder_inputs   = FAST_directory + os.sep + 'Airfoils'
folder_outputs  = folder_inputs
n_stations      = 30
aoa             = np.linspace(-5., 25., 30)             # List of angles of attack
Re              = np.array([0.5e+06, 1.e+06, 5.e+06, 10.e+06])  # List of Reynolds numbers
TEoffset        = 1 # Number of nodes away from the TE where the BL properties are extracted, TEoffset = 0 is prone to convergence issues 
trip            = 0 # Flag to set free (0) or force (1)
TEAngle         = 10.*np.ones(n_stations) # Distribution of trailing edge angles for BPM trailing edge bluntness noise model
TEThick         = 0.005*np.ones(n_stations) # Distribution of trailing edge angles for BPM trailing edge bluntness noise model

## code

inputs_list  = []
outputs_list = []
for i in range(n_stations):
        if i < 10: 
            inputs_list.append(folder_inputs   + os.sep + 'AF0' + str(i) + '_Coords.txt')    # List of files containing the airfoil coordinates
        else:
            inputs_list.append(folder_inputs   + os.sep + 'AF' + str(i) + '_Coords.txt')    # List of files containing the airfoil coordinates
        outputs_list.append(folder_outputs + os.sep + 'AF' + str(i) + '_BL.txt')      # List of files containing the boundary layer characteristics


for id in range(5,len(inputs_list)):
    print('Compute BL properties station ' + str(id))
    filename    = inputs_list[id]
    coord       = np.loadtxt(filename, skiprows = 9)
    np.savetxt('airfoil.dat', coord)
    Ue_Vinf_SS  = np.zeros((len(aoa),len(Re)))
    Ue_Vinf_PS  = np.zeros((len(aoa),len(Re)))
    Dstar_SS    = np.zeros((len(aoa),len(Re)))
    Dstar_PS    = np.zeros((len(aoa),len(Re)))
    Theta_SS    = np.zeros((len(aoa),len(Re)))
    Theta_PS    = np.zeros((len(aoa),len(Re)))
    Cf_SS       = np.zeros((len(aoa),len(Re)))
    Cf_PS       = np.zeros((len(aoa),len(Re)))
    H_SS        = np.zeros((len(aoa),len(Re)))
    H_PS        = np.zeros((len(aoa),len(Re)))
    
    
    for k in range(len(Re)):
        fid = open('inputxfoil.vbs','w')
        fid.write('\n')
        fid.write('load airfoil.dat\n')
        fid.write('pane\n')
        fid.write('plop\n')
        fid.write('g\n')
        fid.write('%f\n'    % 0.8)
        fid.write('\n')
        fid.write('oper\n')
        fid.write('visc %1.3e\n' % (Re[k]))
        fid.write('iter\n')
        fid.write('%u\n'    % 200)

        if trip == 1:
            fid.write('vpar\n')
            fid.write('xtr 0.05 0.1\n')
            fid.write('\n')
        
        for j in range(len(aoa)):
            fid.write('alfa %3.2f\n'        % aoa[j])
            fid.write('dump airfoil.bl%u\n' % j)
        fid.write('\n')
        fid.write('quit\n')
        fid.close()
        
        os.system(path2xfoil + ' < inputxfoil.vbs')
        
        os.remove('inputxfoil.vbs')

        
        for j in range(len(aoa)):
            s_counter=0
            i_TE_PS = []
            bl_data = np.array([])
            text_file = 'airfoil.bl' + str(j)
            with open(text_file) as xfile:
                xfile.readline()
                for line in xfile:
                    data = np.array([float(x) for x in line.strip().split()])
                    if s_counter==0:
                        bl_data = data
                    else:
                        if len(data) == 12:
                            bl_data = np.vstack((bl_data,data))
                        elif len(data) == 8:
                            data = np.hstack((data, np.zeros(4)))
                            if i_TE_PS == []:
                                i_TE_PS = s_counter - 1
                        else:
                            pass

                    s_counter +=1
            
            
            Ue_Vinf_SS[j,k] = bl_data[0  + TEoffset , 3]
            Ue_Vinf_PS[j,k] = bl_data[i_TE_PS - TEoffset , 3]
            Dstar_SS[j,k]   = bl_data[0  + TEoffset , 4]
            Dstar_PS[j,k]   = bl_data[i_TE_PS - TEoffset , 4]
            Theta_SS[j,k]   = bl_data[0  + TEoffset , 5]
            Theta_PS[j,k]   = bl_data[i_TE_PS - TEoffset , 5]
            Cf_SS[j,k]      = bl_data[0  + TEoffset , 6]
            Cf_PS[j,k]      = bl_data[i_TE_PS - TEoffset , 6]
            H_SS[j,k]       = bl_data[0  + TEoffset , 7]
            H_PS[j,k]       = bl_data[i_TE_PS - TEoffset , 7]
            
            os.remove('airfoil.bl' + str(j))
    
    os.remove('airfoil.dat')
    
    # Compute Delta (nominal boundary layer thickness) from Dstar (boundary layer displacement thickness), Theta (boundary layer momentum thickness), and H (kinematic shape factor)
    Delta_SS = Theta_SS*(3.15+1.72/(H_SS-1.))+Dstar_SS
    Delta_PS = Theta_PS*(3.15+1.72/(H_PS-1.))+Dstar_PS

    fid=open(outputs_list[id],'w')
    fid.write('! Boundary layer characteristics at the trailing edge for the airfoil coordinates of %s\n' % filename)
    fid.write('! Legend: aoa - angle of attack (deg), Re - Reynolds number (-, millions), PS - pressure side, SS - suction side,  Ue_Vinf - edge velocity (-), Dstar - displacement thickness (-), Delta - nominal boundary layer thickness (-) Cf  - friction coefficient (-)\n')
    fid.write('%u \t ReListBL   -  Number of Reynolds numbers (it corresponds to the number of tables)\n' % len(Re))
    fid.write('%u \t aoaListBL  -  Number of angles of attack (it corresponds to the number of rows in each table)\n' % len(aoa))
    for k in range(len(Re)):
        fid.write('%1.2f \t   -  Re\n' % (Re[k]*1.e-6))
        fid.write('aoa \t \t Ue_Vinf_SS \t Ue_Vinf_PS \t  Dstar_SS \t \t  Dstar_PS \t \t  Delta_SS \t \t  Delta_PS \t \t    Cf_SS \t \t    Cf_PS\n')
        fid.write('(deg) \t \t \t (-) \t \t \t (-) \t \t \t (-) \t \t \t (-) \t \t \t (-) \t \t \t (-) \t \t \t (-) \t \t \t (-)\n')
        for j in range(len(aoa)):
            fid.write('%1.5f \t %1.5e \t %1.5e \t %1.5e \t %1.5e \t %1.5e \t %1.5e \t %1.5e \t %1.5e \n' % (aoa[j], Ue_Vinf_SS[j,k], Ue_Vinf_PS[j,k], Dstar_SS[j,k], Dstar_PS[j,k], Delta_SS[j,k], Delta_PS[j,k], Cf_SS[j,k], Cf_PS[j,k]))
    
    fid.write('\n')
    fid.write('! Inputs to trailing edge bluntness noise model\n')
    fid.write('%1.5f \t TEAngle    - Angle of the trailing edge (deg)\n'%TEAngle[id])
    fid.write('%1.5f \t TEThick    - Finite thickness of the trailing edge (deg)\n'%TEThick[id])

    fid.close()
    
    