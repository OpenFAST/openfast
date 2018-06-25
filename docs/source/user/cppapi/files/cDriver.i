# -*- mode: yaml -*-
#
# C++ glue-code for OpenFAST - Example input file
#

#Total number of turbines in the simulation
n_turbines_glob: 3
#Enable debug outputs if set to true
debug: False
#The simulation will not run if dryRun is set to true
dry_run:  False
#Flag indicating whether the simulation starts from scratch or restart
sim_start: init # init/trueRestart/restartDriverInitFAST
#Coupling mode
coupling_mode: strong # strong/classic
#Start time of the simulation
t_start:  0.0
#End time of the simulation. tEnd <= tMax
t_end:    1.0
#Max time of the simulation
t_max:    4.0
#Time step for FAST. All turbines should have the same time step.
dt_FAST:  0.00625
#Number of substeps per timestep of the glue-code
n_substeps: 1
#Restart files will be written every so many time steps
n_checkpoint: 160
#Set velocity at the the turbine using an exponential law profile.
set_exp_law_wind: false

Turbine0:
  #The position of the turbine base for actuator-line simulations
  turbine_base_pos: [ 0.0, 0.0, 0.0 ]
  #The number of actuator points along each blade for actuator-line simulations
  num_force_pts_blade: 0
  #The number of actuator points along the tower for actuator-line simulations.
  num_force_pts_tower: 0
  #The checkpoint file for this turbine when restarting a simulation
  restart_filename: "banana"
  #The FAST input file for this turbine
  FAST_input_filename: "t1_Test05.fst"
  #A unique turbine id for each turbine
  turb_id:  1

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
