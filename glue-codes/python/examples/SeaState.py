

from OpynFAST.seastate import SeaStateLib
import matplotlib.pyplot as plt
import numpy as np

project_root = '/Users/rmudafor/Development/openfast'
library_path = project_root + '/build/modules/seastate/libseastate_c_binding.dylib'

if __name__=="__main__":

    dt = 1.0
    time_steps = 10

    seastatelib = SeaStateLib(
        library_path,
        "NRELOffshrBsline5MW_OC4DeepCwindSemi_SeaState_WaveMod5.dat"
    )

    seastatelib.init(
        time_interval=dt,
        n_steps=time_steps,
    )

    seastate_outputs = np.zeros((time_steps, seastatelib.num_outs))
    for i in range(time_steps):
        seastatelib.calc_output(i)
        seastate_outputs[i] = seastatelib.output_values
        print(i, [f"{value:3.4f} - " for value in seastate_outputs[i]])
    seastatelib.end()

    # Plot the results
    # plt.figure(figsize=(10, 6))
    # plt.plot(seastate_outputs[:, 0])
    # plt.plot(seastate_outputs[:, 1])
    # plt.xlabel('Time Step')
    # plt.ylabel('Value')
    # plt.title('Sea State Outputs')
    # plt.legend()
    # plt.show()
