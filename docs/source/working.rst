.. _working_with_OF:

Working with OpenFAST
=====================

This section provides support for some of the typical use cases of OpenFAST. 
It assumes that the user has a working version of OpenFAST (see :ref:`installation` for installation).
Links to scripts and software that can help the analyses are given.


Running and troubleshooting a simulation
----------------------------------------


Running OpenFAST
~~~~~~~~~~~~~~~~

OpenFAST is usually run from the command line. The simplest method to run OpenFAST is to copy the executable and dlls (e.g. `MAP.dll`) into your working directory, and then open a command line into that directory. You can navigate between folders using the command `cd` on all platforms.
You can visualize the command line by running the following

.. code-block:: bash
    
    cd path/to/working/directory
    ./openfast /h

This will display the version of OpenFAST, the compilation options, and the syntax to run OpenFAST.
It is important to keep track of the version of OpenFAST you are using, since the input files might change from version to version (see :ref:`api_change`).  Note the `./` at the beginning of the command to indicate that the executable is located in the current directory. To avoid having to copy the executable in your working directory, you can place the executable and dlls into a folder and add this folder to your system path. If you chose this method, are restart your terminal, you should be able to run `openfast /h` from any folder, and this time, `./` is not needed. In the following we will omit the `./`, but you might need to use it if you haven't added the executable to your path. 

The typical syntax to run a simulation is: 

.. code-block:: bash
    
    openfast InputFile.fst
    
where `InputFile.fst` is a main OpenFAST input file. The extension `.fst` is recommended to maintain a convention. On certain platform (like Windows), you can drag and drop an input file to the OpenFAST executable in your file explorer, and this will run the simulation. If an error occurs using this method, you will not be able to see the error message. When a simulation runs successfully, OpenFAST generates a `.out` or `.outb` output file. To visualize the outputs, see :ref:`visualizing_input_output_OF`.  If an error is displayed, see :ref:`troubleshooting`.



.. _troubleshooting:

Simple troubleshooting
~~~~~~~~~~~~~~~~~~~~~~

When an error is caught during a simulation, OpenFAST will abort with a message of the following kind:

.. code-block:: bash

    FAST encountered an error during module initialization.
    Simulation error level: FATAL ERROR

    Aborting OpenFAST.

The lines above this message will reveal the nature of the error, and this information can be used to troubleshoot your simulation.
Some typical errors and solutions are listed below:

- *The input file "FILE" was not found*: As indicated, the input file is not found. Linux and Mac platforms are case sensitive, and require forward slashes. OpenFAST accepts relative or absolute path. Relative paths are expressed with respect to the file where they are referenced.

- *Invalid input in file "FILE" while trying to read VAR*: Such errors typically occurs at initialization, when reading one of the input file of OpenFAST. It can be that the variable in the input file has a wrong type (integer instead of logical, float instead of string, etc.). Very often though, such error  indicates that the input file is not at the right version for the version of OpenFAST. Most likely your file is outdated. Lines are often added to the OpenFAST input files. You can have a look at :ref:`api_change` to see what lines have changed between versions of OpenFAST, or look at the `r-test <https://github.com/openfast/r-test>`__ to find working examples of input files for the latest release and dev version of OpenFAST.

- *A fatal error occurred when parsing data from "FILE". The variable "VAR" was not found on line #II*. Such errors are similar to the one described above. Check that your file has the proper format for the version of OpenFAST you are using. 

Similar messages indicate user-input errors (when selected options are not available or compatible).
Such error messages are usually explicit enough. You can have a look at the comments in the input file for some guidance, and refer to the user guide for more details on individual inputs of each module: :ref:`user_guide`.

Some warnings might occasionally occur from different modules (typically the aerodynamic modules) and be reported to the command window. 

 - *SkewedWakeCorrection encountered a large value of chi*: indicates that the turbine is highly yawed/titled. Could happen when the turbine undergoes important motions. 
 - *The BEM solution is being turned off due to low TSR.*: indicate that the instantaneous rotor speed is close to zero, or the relative wind speed is large (check the outputs `RtSpeed` and `RtVavgx`).

The warnings can sometimes be ignored, but they can also indicate an issue in the model, since some of the model reach a point beyond their expected physical range. See the next section of advanced troubleshooting.



Advanced troubleshooting
~~~~~~~~~~~~~~~~~~~~~~~~

In some cases, simulations may abort during the simulation (*FAST encountered an error at simulation time T*), or they may run through but have empty or "NaN" inputs after few time steps (as little as one time steps). Such errors are typically due to the model being unphysical.
In such case, you might see error messages of the following kind in the command window:

- *Small angle assumption violated* or *Angles in GetSmllRotAngs() are larger than 0.4 radians*: such warnings indicate that part of the structure is undergoing large rotations, whereas some module of OpenFAST are only valid under the small angle approximation. 
- *Denominator is zero in GetSmllRotAngs()*

Typically, when a simulation aborts or has unrealistic or NaN values, it is likely that there are errors in the model (the structure is too stiff, too soft, the inflow is incorrect, the initial conditions are incorrect, the controller is behaving unexpectedly, OLAF regularization parameters are set wrong, etc.).

Below are some steps you can take to troubleshoot your model, in particular trying to isolate the problem to a given module and then input:


- Simplify the model by using simple environmental conditions: steady uniform inflow, still water.

- Remove the controller: Turn `GenDOF` to False in ElastoDyn, and set `CompServo` to 0 in the main input file. The rotor will spin at constant RPM.

- Simplify your model by turning off most degrees of freedom in your ElastoDyn input file. You can start by keeping all degrees of freedom off, and progressively adding more degrees of freedom. This might indicat if the issue comes from the blade, nacelle, tower or substructure. Some degrees of freedom that are often problematic are the drive train torsion (`DrTrDOF`), and the yaw degree of freedom (`YawDOF`). The drive train stiffness and damping values in ElastoDyn are often set wrong. A common issues with yaw, is when `NacYaw` (in ElastoDyn) and `YawNeut` (in ServoDyn), are in disagreement, or, when the yaw spring and damping `YawSpr` and `YawDamp` are not physical. For offshore simulations, if `YawDOF` and `PtfmYDOF` are on, the model needs to have a realistic `PtfmYIner` present, otherwise these degrees of freedom will be ill-defined in ElastoDyn. The Y-inertia, should represent the torsional inertia of the platform/TP (if any) and the one from the tower.

- Simplify the physical models: use ElastoDyn (`CompElast=1`) over BeamDyn, use BEM (`WakeMod=1`) over OLAF, use 0 Craig-Bampton modes in SubDyn.

- Visualize the outputs (see :ref:`visualizing_input_output_OF`). Add relevant displacement outputs to your model for instance: PtfmSurge, PtfmSway, PtfmHeave, PtfmRoll, PtfmPitch, PtfmYaw, NacYaw, TTDspFA, TTDspSS, RotSpeed, OoPDefl1, IPDefl1 and RtSkew. It is likely that the turbine has some large displacements due to some errors in the model. 

- Adjust your initial conditions. As mentioned above, `NacYaw` (ElastoDyn) and `YawNeut` (ServoDyn) need to match when the yaw degrees of freedom is on. If the structural is at an initial position that is unrealistic given the environmental condition, it is likely to overshoot (e.g. high wind speed but pitch too low). 

- Visualize the inputs (see :ref:`visualizing_input_output_OF`). Check that the mass and stiffness distributions of the blade and tower are as expected.

- Verify the masses and stiffness of your system. The Blade mass and tower-top mass are shown in the ElastoDyn summary file. The equivalent 6x6 matrix of the substructure is found in the SubDyn summary file.

- If you have isolated the problem to a given module, check the information provided in the summary file of this module. Most module have a flag at the end of their input file called `SumPrint` or similar, so that the summary file is written to disk. 

- Reduce the time step. The simulation time step needs to be adjusted based on the frequencies that are modelled in the system (typically the time step needs to be at least half or even a tenth of the fastest frequency). Modules like BeamDyn and SubDyn usually require fine time steps.
  Instead of reducing the time step, it is often equivalent to introduce 1 correction step (`NumCrctn`). When corrections are used the Jacobian need to be updated regularly, for instance setting `DT_UJac` to 100 time steps.

- Perform a linearization of your structure in vacuum (`CompInflow=0`, `CompAero=0`) and in standstill (`RotSpeed=0`) (see :ref:`linearization_analysis_OF`) and check that the frequencies and damping are within the range you expect. Adjust your structural inputs otherwise.

- Generate VTK outputs for visualization of the turbine and the various meshes used by OpenFAST. VTK outputs are activated using `WrVTK=1` or `WrVTK=2`. The VTK are written in folders `vtk*` in the main directory, and can be visualized using Paraview (see :ref:`visualizing_input_output_OF`).


Additional guidelines may be added in the future. In particular, troubleshooting dedicated to each module.






Scripting
---------

NREL maintains several repositories of scripts to work with OpenFAST. 
The scripts can for instance be used to read the input and outputs of OpenFAST, visualize them, and generate multiple simulation inputs, and postprocess them. Some of these applications will be detailed in the following sections.


The repositories maintained by NREL are the following:

- `python-toolbox <https://github.com/OpenFAST/python-toolbox>`__:  collection of python tools to work with OpenFAST.

- `matlab-toolbox <https://github.com/OpenFAST/matlab-toolbox>`__: collection of matlab tools to work with OpenFAST. , e.g. `FAST2Matlab.m` and `Matlab2FAST.m`
  
- `WEIS <https://github.com/WEIS>`__ : Wind Energy with Integrated Servo-control, performs multifidelity co-design of wind turbines. WEIS is a framework that combines multiple NREL-developed tools to enable design optimization of floating offshore wind turbines.

The users are invited to consult the documentations of the individual repository, and discuss related issues on their individual github pages. Contribution by the community to the NREL repositories are welcome and encouraged.



Additional repositories maintained by NREL are listed below:

- `WISDEM <https://github.com/WISDEM/WISDEM>`__: models for assessing overall wind plant cost of energy (COE), also contains file IO, (DLC) case generation, polar manipulations, visualization, and much more! 
- `ROSCO_toolbox <https://github.com/NREL/ROSCO_toolbox>`__: tools to work with the `ROSCO <https://github.com/NREL/ROSCO>`__ controller that is supported by OpenFAST



Repositories maintained by third-parties are listed below:

- `WindEnergyToolbox <https://gitlab.windenergy.dtu.dk/toolbox/WindEnergyToolbox>`_: library developed by DTU, providing some support for different file formats

- `FASTTool <https://github.com/TUDelft-DataDrivenControl/FASTTool>`_ : NREL FASTv8, MATLAB GUI and Simulink integration developed by TUDelft




 
.. _models_OF:

Open-source OpenFAST models
---------------------------

Open-source OpenFAST wind turbine models can be found here:

- `openfast-turbine-models <https://github.com/NREL/openfast-turbine-models>`_: open source wind turbine models (in development)
- `r-test <https://github.com/OpenFAST/r-test>`__: regression tests for OpenFAST, contains models for OpenFAST and its drivers (AeroDyn, SubDyn, HydroDyn, etc.). This repository is not intended to be used as a "database" of models, but it has the advantage that the input files are always up to date with the latest `format specifications <https://openfast.readthedocs.io/en/master/source/user/api_change.html>`_






.. _visualizing_input_output_OF:

Visualizing inputs and outputs files
------------------------------------



To visualize the input and output files of OpenFAST the following graphical interface tool can be used:

- `pyDatView <https://github.com/ebranlard/pyDatView>`_ : tool to plot the input and output files of OpenFAST, CSV-files, and other files from other wind energy software (Hawc2, Flex, Bladed). Multiple files can be opened at once to compare results from different simulations.

The VTK visualization files that are written by OpenFAST can be opened using:

- `paraview <https://www.paraview.org/>`_ : tool to open the VTK files generated by OpenFAST, i.e. velocity fields and turbine geometry.


For advanced cases, the user may want to script the reading and plotting of the input files.
Python and Matlab tools are respectively being provided in the `python-toolbox <https://github.com/OpenFAST/python-toolbox>`_ and `matlab-toolbox <https://github.com/OpenFAST/matlab-toolbox>`_. 
In the matlab toolbox, the scripts `FAST2Matlab.m` and `Matlab2FAST.m` are used to read and write input files, the script `ReadFASTbinary` is used to open binary (`.outb`) output files. 
The README files of these repositories points to examples and more documentation.
  



.. _running_multiple_OF:

Running parametric studies and design load cases (DLC)
------------------------------------------------------

Parametric studies can be run by using the scripts to read and write OpenFAST input files provided in the `matlab-toolbox <https://github.com/OpenFAST/matlab-toolbox>`__
and 
`python-toolbox <https://github.com/OpenFAST/python-toolbox>`__
.  The python-toolbox provides dedicated scripts and examples to automatize the process (see the README of the repository for more).
The `AeroelasticSE` module of `WEIS <https://github.com/WEIS>`__ can generate input files for the design load cases specified in the standards. In the future, `AeroelasticSE` will be likely ported to the python-toolbox.





.. _linearization_analysis_OF:

Performing linearization analyses
---------------------------------



Background
~~~~~~~~~~

Many applications require a linear model of a system: eigenvalue analyses, frequency domain analysis, linear state space models for observers, etc. Most models of OpenFAST are non-linear, and a linearization of the underlying system is therefore required. 
Linearization is done about a given operating point, which corresponds to the set of values of the states and inputs of the system (typically, a given time of a simulation). 
The output of the linearization is a linear state space model (four matrices relating states, inputs and outputs) valid in the neighborhood of the operating point.

When the rotor is spinning, the linearized model obtained at a given time cannot be directly interpreted since it is effectively in a rotating frame of reference. It is necessary to linearize at different operating points over a period of revolution (i.e. at different azimuthal positions) to be able to obtain a linear state space model of the system that is in a fixed (non-rotating) frame of reference. The linearized outputs at different azimuthal positions are combined to form a linearized system in a fixed frame of reference. This manipulation is referred to as the multi-blade coordinate transformation (MBC). For this approach to be valid, the system must be in a periodic steady state.

One of the output of the linearization is the state matrix (`A`) which relates the system states to their time derivatives. An eigenvalue analysis of the state matrix provide the frequencies and damping of the system. Unlike a linear finite-element software, OpenFAST does not have a notion of stiffness and mass matrix. The underlying system of equation is non-linear, the frequencies of the system will vary with the operating conditions (e.g. wind speed, rotational speed).
The sections below detail the process of obtaining a linear model with OpenFAST, and will focus on its application to obtain the frequencies and damping of the system modes.




Linearized models for one simulation (manually)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This section describes the key steps to generate a linearized model of the system with OpenFAST. 

The steps to perform simple linearization analysis are given below:

- Edit the main `.fst` file, set `Linearize=True`

- There are two main methods to determine at which times the linearization will be made: 1) using `CalcSteady=False`, the user prescribes the times where linearization is to occur using `NLinTimes` and `LinTimes` (it is the responsibility of the user to provide times where the system is in equilibrium or at a periodic steady state, i.e. sufficiently long time); 2) `CalcSteady=True` (recommended approach), OpenFAST will automatically start the linearization when the system is at a periodic steady state and will perform `NLinTimes` linearizations over a rotor revolution.
  When a controller is used the option `CalcSteady` will also adjust the controller inputs (either Pitch, Yaw, or Generator Torque, based on the input `TrimCase`) such as to reach the rotational speed indicated by the initial condition.  The `TrimGain` and `TrimTol` might need to be adjusted. 

- Chose the number of linearizations. For a standstill case, `NLinTimes=1`, for a rotating case, it is recommended to use `NLinTimes=36`. If `CalcSteady=False`, the user needs to input `LinTimes` and give values that corresponds to the rotor being at 36 unique azimuthal position (based on the rotor speed).

- For a typical linearization, the user may set `LinInputs=1`, `LinOutputs=1`, `LinOutJac=False`, `LinOutMod=False`, `Twr_Kdmp=0`, `Bld_Kdmp=0`. The standard set of linearization inputs inherent in the linearized system are available when `LinInputs=1`. This includes e.g. collective blade pitch. With `LinOutputs = 1`, every output various `OutList` sections of each module are included in the linearized system. For instance, `GenSpeed` can be included by including `GenSpeed` in the `OutList` of ElastoDyn.

- Set the output format `OutFmt` to "ES20.11E3". The output files will be written with this high resolution, which is required for accurate eigenvalue analyses.

- Currently linearization of unsteady aerodynamics are not fully supported. It is therefore recommended to use `FrozenWake=True` in the AeroDyn input file. In the future, the user will be able to select a special `DBEMT` and `UAMod` options to account for unsteady aerodynamics while linearizing.

- Run OpenFAST on this `.fst` file. OpenFAST will report errors if some options are activated that are not supported with linearization (see section below). Adapt your input files accordingly. When the simulation runs OpenFAST will display a message when performing each individual linearizations, and individual files with the `.lin` extension will be written.

- It is recommended to check the regular output file `.out` or `.outb`. If `CalcSteady=False`, the user should look to see whether the turbine had indeed reached a steady state (or periodic steady state) at the time where linearization was run. If `CalcSteady=True` and a controller is used, the user can check that the rotational speed has indeed converged to the desired RPM, and potentially chose to adjust `TrimGain` and `TrimTol` for future runs.

The linearization files `*.lin` are then to be postprocessed using the python or matlab tools provided.


Postprocessing 
~~~~~~~~~~~~~~

To obtain the eigen-frequencies of the system the user can open a `.lin` file, extract the state matrix `A` and perform a eigenvalue analysis. For a spinning rotors, all lin-files generated from a simulation at different azimuthal positions need to be opened, and converted using the MBC-transformation. We provide scripts for such cases.

When only one linearization file is to be used (e.g. at standstill), the script `postproLin_OneLinFile_NoRotation` can be used. Is is found in `matlab-toolbox/Campbell/example` or `python-toolbox/pyFAST/linearization/examples/`.

When several linearization files are to be postprocessed (in particular several files corresponding to different azimuthal positions), the script `postproLin_MultiLinFile_Campbel` can be used, located in the same folders mentioned above.
The script can also be used if linearizations were performed at different wind speed and RPM (via different OpenFAST calls). Displaying the frequencies and damping at these different wind turbine operating conditions is referred to as Campbell diagram.



Campbell diagrams
~~~~~~~~~~~~~~~~~

To avoid the manual process of editing input files for different wind turbine operating conditions, we provide the script `runCampbell`, found in `matlab-toolbox/Campbell/example` or `python-toolbox/pyFAST/linearization/examples/`.
The script relies on a template folder which a reference "fst" file. The folder is duplicated, files are created for each wind turbine operating conditions wind speed/rpm), OpenFAST is run, and the linearization files are postprocessed.

The script `runCampbell` generates either a set of CSV files or an Excel file. The script attempts to identify the modes (for instance 1st tower fore-aft mode, 1st flap mode, etc.), but a manual process is usually required to fully identify the mode. This process can be difficult and tedious. It is recommended to proceed first with simulations in vacuum, and with few operating points, to get familiar with the system.

The manual identification process consists in changing the CSV file `Campbell_ModesID.csv` (or the Excel spread sheet `ModesID` if Excel output is used). To avoid having this file rewritten when rerunning `runCampbell`, it is recommended to rename this file as `Campbell_ModesID_Manual.csv`.  The part of the script `runCampbell` that plots the Campbell diagram can be adjusted so as to use the "Manual" file. 
It is recommended to use the CSV format since this is the method compatible with Python and MacOS.

The manual identification process consists in attributing indexes in the table of modes, where the index corresponds to the list of sorted mode frequencies.

For instance, opening the CSV file in excel, the `ModeID` file might look as follows:

.. code::

    Mode Number Table      
    Wind Speed (mps)   2.0   5.0   8.0
    1st Tower FA        0     0     0
    1st Tower SS        1     0     0

In this example, we  assume that linearizations were run at 2, 5 and 8m/s. "0" in the table indicates that a mode was not identified. You can look at the file `Campbell_Summary.txt` to have a look at the frequencies, damping and "modal content" for each mode and operating point. For more details, you can open the individual CSV files for each operating point. (If you used the Excel format, these are in different sheets).
You might find that for 2 and 5m/s, the tower Fore-Aft is the second frequency, and the tower side-side is the first frequency that shows up in the list of modes. At 8m/s you might find that the opposite occurs. In that case, you will edit the file such that it is as follows:

.. code::

    Mode Number Table      
    Wind Speed (mps)   2.0   5.0   8.0
    1st Tower FA        2     2     1
    1st Tower SS        1     1     2


The main question is how to determine which mode is which. There is no true solution to this question, here are some elements to help the identifications: 

 - The system frequencies are usually easy to determine at 0 m/s and 0 rpm. The system frequencies will vary progressively from this reference point as the RPM/WS/pitch changes. Blade modes will typically display a "splitting" equal to +/- the rotational speed frequency as the rotational speed increases.

 - Blade flap modes are typically highly damped 

 - From an operating point to the next, the damping will not change drastically.

 - Tower modes are not strongly affected by the change of operating conditions

 - You will need to look at the "mode content", to see where the energy is for each mode. The file `Campbell_Summary.csv` displays a summary of the mode content. In some cases, there is no clear maximum (the keyword `NoMax` is shown). In that case, identifying the mode might be difficult. A similar content is found in the individual operating point files.

 - Visualization of the modes can help identify them (see the next section). The process can yet be lengthy.

Once the identification table is set. Save the file, and plot the Campbell diagram. The process may be iterative until a satisfying diagram is obtained. There should be no need to close Excel in this process.

We are aware that the process is lengthy,  we thank you for your patience while we attempt to streamline this process.



Mode shape visualization
~~~~~~~~~~~~~~~~~~~~~~~~

Mode shape visualization is currently possible. It requires a generation of viz files for each simulations, and rerunning OpenFAST to generate VTK files. The matlab script `runCampbell` assists in this process, but for now limited support and documentation is provided.

The user is invited to consult the following example:
-  https://github.com/OpenFAST/r-test/tree/main/glue-codes/openfast/5MW_Land_ModeShapes

And it's associated documentation:
- https://github.com/OpenFAST/r-test/blob/main/glue-codes/openfast/5MW_Land_ModeShapes/vtk-visualization.md


Additional references
~~~~~~~~~~~~~~~~~~~~~

Some linearization issues have been discussed in the forum and as github issues:

- https://wind.nrel.gov/forum/wind/
  
- https://github.com/OpenFAST/openfast/issues/480

Thank you for your patience while we attempt to streamline this process.





