# -*- mode: yaml -*-
#
# C++ glue-code for OpenFAST - Example input file
#

n_turbines_glob: 3       # Total number of turbines in the simulation

debug: False             # Enable debug outputs if set to true

dry_run: False           # The simulation will not run if dryRun is set to true

sim_start: init          # Flag indicating whether the simulation starts from scratch or restart
                         # [init | trueRestart | restartDriverInitFAST]

coupling_mode: strong    #  Coupling mode
                         # [strong | classic]

t_start: 0.0             # Start time of the simulation

t_end: 1.0               # End time of the simulation; tEnd <= tMax.

t_max: 4.0               # Max time of the simulation

dt_FAST: 0.00625         # Time step for FAST. All turbines should have the same time step.

n_substeps: 1            # Number of substeps per timestep of the glue-code

n_checkpoint: 160        # Restart files will be written every so many time steps

set_exp_law_wind: false  # Set velocity at the the turbine using an exponential law profile.

Turbine0:

  turbine_base_pos: [ 0.0, 0.0, 0.0 ]  # The position of the turbine base for actuator-line simulations

  num_force_pts_blade: 0               # The number of actuator points along each blade for actuator-line simulations

  num_force_pts_tower: 0               # The number of actuator points along the tower for actuator-line simulations.

  restart_filename: "banana"           # The checkpoint file for this turbine when restarting a simulation

  FAST_input_filename: "t1_Test05.fst" # The FAST input file for this turbine

  turb_id:  1                          # A unique turbine id for each turbine

Turbine1:
  turbine_base_pos: [ 0.0, 0.0, 0.0 ]
  num_force_pts_blade: 0
  num_force_pts_tower: 0
  restart_filename: "banana"
  FAST_input_filename: "t2_Test05.fst"
  turb_id:  2

Turbine2:
  turbine_base_pos: [ 0.0, 0.0, 0.0 ]
  num_force_pts_blade: 0
  num_force_pts_tower: 0
  restart_filename: "banana"
  FAST_input_filename: "t3_Test05.fst"
  turb_id:  3
