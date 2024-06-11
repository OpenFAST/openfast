""" 
Reads a VTK file containing a velocity Grid/Plane.
The format is typically used as inputs/outputs of FAST.Farm or outputs of OLAF
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
# Local 
import pyFAST.input_output as io 
import pyFAST.input_output.vtk_file 

# Get current directory so this script can be called from any location
MyDir=os.path.dirname(__file__)
# Read Plane
vtkFileName = os.path.join(MyDir,'../../../data/example_files/Plane.vtk')
vtk = io.vtk_file.VTKFile(vtkFileName)
# Extract field into individual components
u = vtk.point_data_grid['Velocity'][:,:,:,0]
v = vtk.point_data_grid['Velocity'][:,:,:,1]
w = vtk.point_data_grid['Velocity'][:,:,:,2]


if __name__ == '__main__':
    # Print useful information
    print(vtk)
    # Plot a cross section
    fig,ax = plt.subplots()
    im=ax.contourf(vtk.xp_grid, vtk.zp_grid, u[:,0,:].T)
    fig.colorbar(im)
    ax.set_xlabel('x [m]')
    ax.set_ylabel('z [m]')
    ax.set_title('Streamwise velocity in vertical plane (simple shear)')
    plt.show()

if __name__=='__test__':
    np.testing.assert_array_equal(u.shape, (4,1,6))
    np.testing.assert_almost_equal(u[0,-1,3],10.183994 )
    pass
