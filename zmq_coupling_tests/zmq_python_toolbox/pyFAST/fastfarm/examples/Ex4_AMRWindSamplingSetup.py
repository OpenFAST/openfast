"""
Set up sampling planes for AMR-Wind to use for inflow winds.
"""

from pyFAST.fastfarm.AMRWindSimulation import AMRWindSimulation

def main():
    # -----------------------------------------------------------------------------
    # USER INPUT: Modify these
    # -----------------------------------------------------------------------------

    # ----------- Wind farm
    wts  = {
            0 :{'x':1280.0,  'y':2560,  'z':0.0,  'D':126.9, 'zhub':86.5, 'cmax':5, 'fmax':10/6, 'Cmeander':1.9, 'name':'T0'},
            1 :{'x':1280.0,  'y':3200,  'z':0.0,  'D':126.9, 'zhub':86.5, 'cmax':5, 'fmax':10/6, 'Cmeander':1.9, 'name':'T1'},
            2 :{'x':1280.0,  'y':3840,  'z':0.0,  'D':126.9, 'zhub':86.5, 'cmax':5, 'fmax':10/6, 'Cmeander':1.9, 'name':'T2'},
            }

    # ----------- AMR-Wind parameters
    fixed_dt = 0.25
    prob_lo = (0.0, 0.0, 0.0)
    prob_hi = (2560.0, 6400.0, 1280.0)
    n_cell = (256, 640, 128)
    max_level = 1  # Number of grid refinement levels

    incflo_velocity_hh = (0.0, 10.0, 0.0)  # Hub-height velocity
    postproc_name = 'sampling'

    # ----------- Wake model (1: Polar; 2: Curl; 3: Cartesian)
    mod_wake = 1

    # -----------------------------------------------------------------------------
    # END OF USER INPUT
    # -----------------------------------------------------------------------------

    # Initial setup
    amr = AMRWindSimulation(wts, fixed_dt, prob_lo, prob_hi,
                           n_cell, max_level, incflo_velocity_hh,
                           postproc_name,
                           mod_wake = wake_mod)

    # You can print the complete instance by simply
    print(amr)

    # The output of AMRWindSimulation can be directly used as input to FFCaseCreation:
    # FFCaseCreation(..., dt_high_les=amr.dt_high_les, ds_high_les=amr.ds_high_les,
    #                     dt_low_les=amr.dt_low_les, ds_low_les=amr.ds_low_les,
    #                     extent_high=amr.extent_high, extent_low=amr.extent_low,
    #                       ...)

if __name__ == '__main__':
    main()
