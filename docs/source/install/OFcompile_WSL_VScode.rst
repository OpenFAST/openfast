Guide: Setting Up OpenFAST on Windows with WSL and VS Code
============================================================

.. note::
   Adapted from guide written by `skvibimigger <https://github.com/skvibimigger>`_ in an issue comment: `OpenFAST Issue #2904 <https://github.com/OpenFAST/openfast/issues/2904#issuecomment-3056252458>`_

This guide provides a comprehensive, step-by-step process for compiling OpenFAST from source on a Windows machine using the Windows Subsystem for Linux (WSL). This method is highly reliable and avoids common issues with native Windows compilers.

Part 1: Setting Up the Linux Environment (WSL)
-----------------------------------------------

1. Install Windows Subsystem for Linux (WSL)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

WSL allows you to run a Linux distribution directly on Windows.

- Open Windows PowerShell or Command Prompt as an Administrator.
- Run the installation command:

  .. code-block:: bash

     wsl --install

- Restart your computer when prompted. After rebooting, an Ubuntu terminal will open to complete the setup. You will be asked to create a Unix username and password.

2. Install Build Tools and Dependencies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Once your Ubuntu terminal is running, you need to install the necessary compilers and libraries.

- Update Package Lists:

  .. code-block:: bash

     sudo apt-get update

- Install Compilers and Build System:

  .. code-block:: bash

     sudo apt-get install -y build-essential gfortran cmake git

- Install Required Math Libraries:

  .. code-block:: bash

     sudo apt-get install -y libblas-dev liblapack-dev

- (Optional) Install Documentation Dependencies:

  .. warning::
     The ``texlive-full`` package is several gigabytes and may take a long time to install.

  .. code-block:: bash

     sudo apt-get install -y graphviz texlive-full doxygen
     pip install sphinx

Part 2: Compiling OpenFAST
---------------------------

1. Clone the OpenFAST Repository
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Download the main source code:

.. code-block:: bash

   git clone https://github.com/OpenFAST/openfast.git

2. Download the Test Cases (Submodule)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The regression test files (r-test) are in a separate repository (a Git submodule) and must be downloaded explicitly.

- Navigate into the ``openfast`` directory:

  .. code-block:: bash

     cd openfast

- Initialize and download the submodule:

  .. code-block:: bash

     git submodule init
     git submodule update --progress

3. Build with CMake and Make
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Configure the build using CMake and compile the project with make.

- Create a build directory and navigate into it:

  .. code-block:: bash

     mkdir build
     cd build

- Run CMake to configure the build:

  .. code-block:: bash

     cmake .. -DBUILD_SHARED_LIBS=ON -DBUILD_TESTING=ON -DBUILD_DOCUMENTATION=ON

- Compile the project:

  .. code-block:: bash

     make

After this completes without errors, OpenFAST is successfully built.

Part 3: Setting Up the Python Environment
------------------------------------------

1. Install a Compatible Python Version
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The default Python version in new Ubuntu releases can be too new for some packages. Install Python 3.11:

- Add the "deadsnakes" repository:

  .. code-block:: bash

     sudo add-apt-repository ppa:deadsnakes/ppa

- Install Python 3.11 and its virtual environment module:

  .. code-block:: bash

     sudo apt update
     sudo apt install python3.11 python3.11-venv

2. Create and Activate a Virtual Environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Create a self-contained environment for your Python packages.

- Navigate to your home directory:

  .. code-block:: bash

     cd ~

- Create the virtual environment:

  .. code-block:: bash

     python3.11 -m venv openfast_env

- Activate the environment:

  .. code-block:: bash

     source openfast_env/bin/activate

3. Install Python Packages
~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Install ``pyOpenFAST`` from local source:

  .. code-block:: bash

     cd ~/openfast/glue-codes/python/
     pip install .

- Install Jupyter and other libraries:

  .. code-block:: bash

     pip install jupyterlab pandas matplotlib

Part 4: Configuring and Using VS Code
--------------------------------------

1. Install VS Code Extensions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Open VS Code on Windows.
- Go to the Extensions tab.
- Install the following official Microsoft extensions:

  - WSL
  - Jupyter
  - Python

2. Connect VS Code to WSL
~~~~~~~~~~~~~~~~~~~~~~~~~

- Click the green >< button in the bottom-left corner of the VS Code window.
- Select "Connect to WSL" from the menu. The window will reload.

3. Create and Run a Notebook
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- In the connected VS Code window, open your project folder (e.g., ``/home/YOUR_USERNAME``).
- Create a new file named ``simulation.ipynb``.
- When the notebook opens, click "Select Kernel" in the top-right corner.
- Choose "Python Environments..." and select the kernel associated with your virtual environment: ``openfast_env (Python 3.11.x)``.

Part 5: Example Test Script
----------------------------

Paste this code into a cell in your notebook to run a test simulation and plot the results.

.. code-block:: python

   # Import the correct class from the correct file
   from pyOpenFAST.fast import FastLibAPI
   import os
   import matplotlib.pyplot as plt

   print("--- Step 1: Setting up paths ---")

   # The OpenFAST library file
   openfast_library_path = os.path.expanduser('~/openfast/build/modules/openfast-library/libopenfastlib.so')

   # Use the minimal example test case
   sample_input_file = '/home/YOUR_USERNAME/openfast/reg_tests/r-test/glue-codes/openfast/MinimalExample/Main.fst'

   print("\n--- Step 2: Running OpenFAST ---")
   try:
           # Initialize the OpenFAST library interface
           fast_api = FastLibAPI(
                   library_path=openfast_library_path,
                   input_file_name=sample_input_file
           )

           # Run the simulation
           fast_api.run()
           print("✅ OpenFAST simulation ran successfully!")

           # --- Step 3: Post-processing ---
           print("\n--- Step 3: Accessing results ---")

           results = fast_api.output_values
           channels = fast_api.output_channel_names

           print(f"Simulation produced {results.shape[0]} time steps and {results.shape[1]} output channels.")
           print(f"Available channels: {channels}")

           # Define the channel to plot (choose one from the "Available channels" list)
           channel_to_plot = 'OoPDefl1'

           # Find the index of the channel
           try:
                   channel_idx = channels.index(channel_to_plot)
                   time = results[:, 0]
                   channel_data = results[:, channel_idx]

                   # Plot the results
                   plt.figure(figsize=(10, 5))
                   plt.plot(time, channel_data)
                   plt.xlabel('Time (s)')
                   plt.ylabel(f'{channel_to_plot} (m)')
                   plt.grid(True)
                   plt.title(f'{channel_to_plot} vs. Time')
                   plt.show()
                   print("✅ Plot generated successfully!")

           except ValueError:
                   print(f"❌ Could not find '{channel_to_plot}' in the output channels.")

   except Exception as e:
           print(f"❌ An error occurred during the simulation: {e}")

Notes
-----

- Replace ``YOUR_USERNAME`` with your actual Linux username where applicable.
- Ensure all dependencies are installed correctly to avoid runtime errors.
- For additional help, refer to the `OpenFAST GitHub repository <https://github.com/OpenFAST/openfast>`_.
