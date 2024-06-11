""" 
Example usage for the MannBox class.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from pyFAST.input_output.mannbox_file import MannBoxFile


scriptDir = os.path.dirname(__file__)

# --- Inputs
file_u = os.path.join(scriptDir, '../../../data/example_files/MannBox_32x4x8.bin')
dy = 3  # horizontal spacing [m]
dz = 10 # vertical spacing [m]

# --- Read Mann Box file and display main properties
mb = MannBoxFile(file_u, dy=dy, dz=dz) # create a MannBoxFile object, store it in mb
print(mb)                        # print misc info about the object
print(mb['field'].shape)         # print the shape of your field. 
nx, ny, nz =  mb['field'].shape  # Store the field dimension


# --- Plot values of the field at given indices. Note this example has only 2 values in x
# NOTE: can also be obtained directly using: 
# u = mb.valuesAt(y=5, z=50)
t = mb.t(dx=1, U=10)
u = mb['field'][:, 3, 5] # Value at indices iy=3, iz=5
plt.figure(1)
plt.plot(t, u, label='Input') 
plt.xlabel('time [s]') 
plt.ylabel('u [m/s]')
plt.title('Value of field at chosen location as function of x/t')


# --- Plot vertical profile averaged over all x and y
# NOTE: can also be obtained directly using: 
#    z, means, stds = mb.vertProfile
u_mean_vertical = np.mean(np.mean(mb['field'][:,:,:], axis=0), axis=0)
plt.figure(2)
plt.plot(u_mean_vertical, mb.z, label='Input') 
plt.xlabel('u mean [m/s]') 
plt.ylabel('z [m]')
plt.title('Average vertical profile for all x and y')

# --- Compute the mean over all "x" values for each points in the y-z plane
UmeanYZ = np.mean(mb['field'][:,:,:],axis=0) # This has shape ny x nz

# --- Modify the field  (some examples)
mb['field'] -= UmeanYZ         # remove the mean of all datapoints along x 
mb['field'][:,:,:] *= 2        # multiply the field by 10 everywhere
mb['field'][:,:,:] += UmeanYZ  # add the mean again
mb['field'][:,:,:] += 0.5      # add 0.5 everywhere in the field
mb['field'][:,:,0] = 0         # set the velocity to be zero for the first z index


# --- Plot value of the field again
u = mb['field'][:, 3, 5] # Value at indices iy=3, iz=5
plt.figure(1)
plt.plot(t, u, label='Modified') 
plt.xlabel('time [s]') 
plt.ylabel('u [m/s]')

# --- Plot the vertical profile again
u_mean_vertical2 = np.mean(np.mean(mb['field'][:,:,:], axis=0), axis=0)
plt.figure(2)
plt.plot(u_mean_vertical2, mb.z, '--', label='Modified') 
plt.xlabel('u mean [m/s]') 
plt.ylabel('z [m]')




# --- Write the modified field to disk with another filename
outFile = file_u.replace('.bin','_modified.bin')
mb.write(outFile) 
print('File written: ', outFile)



if __name__ == "__main__":
    plt.show()

if __name__=='__test__':
    try:
        os.remove(outFile)
    except:
        pass

