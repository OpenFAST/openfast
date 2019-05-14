Performance Profiling and Optimization
======================================
The OpenFAST team has been engaged in performance profiling and optimization
work in an effort to improve the time-to-solution performance for the most
computationally expensive use cases. This work is supported by Intel through
its designation of NREL as a Parallel Computing Center (IPCC). Further,
Envision Energy has continuously contributed code and expertise in this area.

The general methods targeted for performance improvements in OpenFAST are:

- Intel tech stack (compiler, math library)
- algorithmic improvements
- memory access optimization
- multhithreading

Test cases
----------
Two OpenFAST test cases have been chosen to provide meaningful and
realistic timing benchmarks. In addition to real-world turbine and
atmospheric models, these cases are computationally expensive and expose
the areas where performance improvements would make a difference.

5MW_Land_BD_DLL_WTurb
~~~~~~~~~~~~~~~~~~~~~
Download files `here <https://github.com/OpenFAST/r-test/tree/dev/glue-codes/openfast/5MW_Land_BD_DLL_WTurb>`__.

The physics modules used in this case are:

- BeamDyn
- InflowWind
- AeroDyn 15
- ServoDyn

This is a land based NREL 5MW turbine simulation using BeamDyn as the
structural module. It simulates 20 seconds with a timestep size of 0.001
seconds and executes in `3m 55s <https://my.cdash.org/testDetails.php?test=40171217&build=1649048>`__
on NREL's peregrine supercomputer.

5MW_OC4Jckt_DLL_WTurb_WavesIrr_MGrowth
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Download files `here <https://github.com/OpenFAST/r-test/tree/dev/glue-codes/openfast/5MW_OC4Jckt_DLL_WTurb_WavesIrr_MGrowth>`__.

This is an offshore, fixed-bottom NREL 5MW turbine simulation with the majority
of the computational expense occuring in the HydroDyn wave dynamics
calculation.

The physics modules used in this case are:

- ElastoDyn
- InflowWind
- AeroDyn 15
- ServoDyn
- HydroDyn
- SubDyn

It simulates 60 seconds with a timestep size of 0.01 seconds and executes in
`20m 27s <https://my.cdash.org/testDetails.php?test=40171219&build=1649048>`__
on NREL's peregrine supercomputer.

Profiling
---------
The OpenFAST test cases were profiled with Intel VTune Amplifier to
identify performance hotspots. This tools provides a clear picture of
bottlenecks in a simulation.

LAPACK
~~~~~~
In the offshore case, the LAPACK usage was found the be the bottleneck.

<SHOW VTUNE SCREENSHOT>

BeamDyn
~~~~~~~
While BeamDyn yields a high fidelity blade calculation, it is a computationally
expensive module. Initial profiling highlighted the `bd_elementmatrixga2`
subroutine, in particular, as a hotspot. However, initial attempts to improve
performance in BeamDyn highlighted needs for alogirthmic improvements
and refinements to the module's data structures.

Results
-------
Though work in ongoing, OpenFAST time-to-solution performance has improved
and the performance potential is better understood.

Specifically, some keys outcomes from the first year are:

- Use of Intel compiler and MKL library provides dramatic speedup over GCC
  and LAPACK

  - Additional significant gains through MKL threading for offshore
    configuration

- Offshore-wind-turbine simulation configuration is poorly load balanced
  across modules

  - Land-based-turbine configuration better balanced
  - OpenMP Tasks are employed to achieve better load-balancing

- OpenMP module-level parallelism provides significant, but limited speed
  up due to imbalance across different module tasks
- Core algorithms need significant modification to enable OpenMP and SIMD
  benefits


Speedup - Intel Tech Stack
~~~~~~~~~~~~~~~~~~~~~~~~~~
By employing the standard Intel developer tools tech stack, a performance
improvement over GNU tools was demonstrated:

======== ================ ===================== ======================================
Compiler Math Library     5MW_Land_BD_DLL_WTurb 5MW_OC4Jckt_DLL_WTurb_WavesIrr_MGrowth
======== ================ ===================== ======================================
GNU      LAPACK           2265 s (1.0x)         673 s (1.0x)
Intel 17 LAPACK           1650 s (1.4x)         251 s (2.7x)
Intel 17 MKL              1235 s (1.8x)         ---
Intel 17 MKL Multitheaded 722 s (3.1x)          ---
======== ================ ===================== ======================================


Speedup - OpenMP at FAST_Solver
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
A performance improvement was domenstrated by adding OpenMP directives to the
`FAST_Solver` module. Although the solution scheme is not well balanced,
parallelizing mesh mapping and calculation routines resulted in the following
speedup:

======== =============== ===================== ======================================
Compiler Math Library    5MW_Land_BD_DLL_WTurb 5MW_OC4Jckt_DLL_WTurb_WavesIrr_MGrowth
======== =============== ===================== ======================================
Intel 17 MKL - 1 thread  1073 s (2.1x)         100 s (6.7x)
Intel 17 MKL - 8 threads 597 s (3.8x)          ---
======== =============== ===================== ======================================


Ongoing Work
------------
The next phase of the OpenFAST performance improvements are focused in two key
areas:

1. Implementing the outcomes from previous work throughout OpenFAST modules and
   glue codes
2. Preparing OpenFAST for efficient execution on Intel's next generation
   platforms

.. Year 2 stuff:

.. Furthermore, NREL is optimizing OpenFAST for the future through profiling on
.. Intel next generation platform (NGP) simulators.

.. bd_5MW_dynamic
.. ~~~~~~~~~~~~~~
.. Download files `here <https://github.com/OpenFAST/r-test/tree/dev/modules/beamdyn/bd_5MW_dynamic>`__.

.. This is a standalone BeamDyn case of the NREL 5MW wind turbine. It simulates 30
.. seconds with a timestep size of 0.002 seconds and executes in 24s on NREL's
.. peregrine supercomputer.

.. BeamDyn dynamic solve

.. Performance Improvements
.. ------------------------
.. BeamDyn chosen as the module to improve from year 1

.. How to improve vectorization

.. BeamDyn Memory Alignment
.. ~~~~~~~~~~~~~~~~~~~~~~~~
.. Work accomplished to align beamdyn types in the dervive types module
.. - Ultimately, this needs to be done in the registry

.. Multithreading
.. ~~~~~~~~~~~~~~
.. OpenMP at the highest level
.. OpenMP added to BeamDyn dynamic solve

.. Speedup
.. -------

.. These are the areas where we have demonstrated performance improvements

.. BeamDyn Dynamic
.. ---------------
.. This improved beamdyn's time to solution by XX%

.. - VTune / Advisor
.. - Vectorization report
.. - SIMD report

.. Optimization Reports
.. The optimization reports provided by the Intel fortran compiler give a static
.. analysis of code optimization. Specifically, the vectorization and openmp
.. reports were analyzed to determine
