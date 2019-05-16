Performance-Profiling and Optimization
======================================
The OpenFAST team has been engaged in performance-profiling and optimization
work in an effort to improve the time-to-solution performance for the most
computationally expensive use cases. This work is supported by Intel® through
its designation of NREL as an
`Intel® Parallel Computing Center (IPCC) <https://software.intel.com/en-us/ipcc>`_.

After initial profiling and hotspot analysis, specific subroutines in the
physics modules of OpenFAST were targeted for optimization. Among other
takeaways, it was learned that the memory alignment of the derived data
types could yield a significant increase in performance. Ultimately, tuning
the Intel® tools to perform best on NREL's hardware and adding high level
multithreading yielded a maximum 3.8x time-to-solution improvement for one
of the benchmark cases.

Approach
--------
The general mechanisms identified for performance improvements in OpenFAST are:

- Intel® compiler suite and Intel® Math Kernel Library (Intel® MKL)
- Algorithmic improvements
- Memory-access optimization enabling more efficient cache usage
- Data type alignment allowing for SIMD vectorization
- Multithreading with OpenMP

To establish a path forward with any of these options, OpenFAST was first
profiled with Intel® VTune™ Amplifier which provides a clear breakdown of
time spent in the simulation. Then, the optimization report generated from the
Intel® Fortran compiler was analyzed to determine area which were not
autovectorized. Finally, Intel® Advisor was used to highlight areas of the code
which the compiler identified as potentially improved with multithreading.

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

This is a land based NREL 5-MW turbine simulation using BeamDyn as the
structural module. It simulates 20 seconds with a time step size of 0.001
seconds and executes in `3m 55s <https://my.cdash.org/testDetails.php?test=40171217&build=1649048>`__
on NREL's `Peregrine <https://www.nrel.gov/hpc/peregrine-system.html>`__
supercomputer.

5MW_OC4Jckt_DLL_WTurb_WavesIrr_MGrowth
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Download files `here <https://github.com/OpenFAST/r-test/tree/dev/glue-codes/openfast/5MW_OC4Jckt_DLL_WTurb_WavesIrr_MGrowth>`__.

This is an offshore, fixed-bottom NREL 5-MW turbine simulation with the
majority of the computational expense occurring in the HydroDyn wave-dynamics
calculation.

The physics modules used in this case are:

- ElastoDyn
- InflowWind
- AeroDyn 15
- ServoDyn
- HydroDyn
- SubDyn

It simulates 60 seconds with a time step size of 0.01 seconds and executes in
`20m 27s <https://my.cdash.org/testDetails.php?test=40171219&build=1649048>`__
on NREL's `Peregrine <https://www.nrel.gov/hpc/peregrine-system.html>`__
supercomputer.

Profiling
---------
The OpenFAST test cases were profiled with Intel® VTune™ Amplifier to
identify performance hotspots. Being that the two test cases exercise
difference portions of the OpenFAST software, different hotspots were
identified. In all cases and environment settings, the majority of the
CPU time was spent in `fast_solution` loop which is a high-level subroutine
that coordinates the solution calculation from each physics module.

LAPACK
~~~~~~
In the offshore case, the LAPACK usage was identified as a performance load.
Within the `fast_solution` loop, the calls to the LAPACK function `dgetrs`
consume 3.3% of the total CPU time.

.. figure:: images/offshore_lapack.png
   :width: 100%
   :align: center

BeamDyn
~~~~~~~
While BeamDyn provides a high-fidelity blade-response calculation, it is a
computationally expensive module. Initial profiling highlighted the
`bd_elementmatrixga2` subroutine, in particular, as a hotspot. However, initial
attempts to improve performance in BeamDyn highlighted needs for algorithmic
improvements and refinements to the module's data structures.

Results
-------
Though work is ongoing, OpenFAST time-to-solution performance has improved
and the performance potential is better understood.

Some keys outcomes from the first year of the IPCC project are as follows:

- Use of Intel® compiler and MKL library provides dramatic speedup over GCC
  and LAPACK

  - Additional significant gains are possible through MKL threading for
    offshore simulations

- Offshore-wind-turbine simulations are poorly load balanced
  across modules

  - Land-based-turbine configuration better balanced
  - OpenMP Tasks are employed to achieve better load-balancing

- OpenMP module-level parallelism provides significant, but limited speed
  up due to imbalance across different module tasks
- Core algorithms need significant modification to enable OpenMP and SIMD
  benefits


Speedup - Intel® Compiler and MKL
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
By employing the standard Intel® developer tools tech stack, a performance
improvement over GNU tools was demonstrated:

========= ================= ===================== ======================================
Compiler  Math Library      5MW_Land_BD_DLL_WTurb 5MW_OC4Jckt_DLL_WTurb_WavesIrr_MGrowth
========= ================= ===================== ======================================
GNU       LAPACK            2265 s (1.0x)         673 s (1.0x)
Intel® 17 LAPACK            1650 s (1.4x)         251 s (2.7x)
Intel® 17 MKL               1235 s (1.8x)         ---
Intel® 17 MKL Multithreaded 722 s (3.1x)          ---
========= ================= ===================== ======================================


Speedup - OpenMP at FAST_Solver
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
A performance improvement was domenstrated by adding OpenMP directives to the
`FAST_Solver` module. Although the solution scheme is not well balanced,
parallelizing mesh mapping and calculation routines resulted in the following
speedup:

========= =============== ===================== ======================================
Compiler  Math Library    5MW_Land_BD_DLL_WTurb 5MW_OC4Jckt_DLL_WTurb_WavesIrr_MGrowth
========= =============== ===================== ======================================
Intel® 17 MKL - 1 thread  1073 s (2.1x)         100 s (6.7x)
Intel® 17 MKL - 8 threads 597 s (3.8x)          ---
========= =============== ===================== ======================================


Ongoing Work
------------
The next phase of the OpenFAST performance improvements are focused in two key
areas:

1. Implementing the outcomes from previous work throughout OpenFAST modules and
   glue codes
2. Preparing OpenFAST for efficient execution on Intel®'s next generation
   platforms

.. Year 2 stuff:

.. Further, `Envision Energy USA, Ltd <http://www.envision-group.com/en/energy.html>`_
.. has continuously contributed code and expertise in this area.


.. Furthermore, NREL is optimizing OpenFAST for the future through profiling on
.. Intel next generation platform (NGP) simulators.

.. bd_5MW_dynamic
.. ~~~~~~~~~~~~~~
.. Download files `here <https://github.com/OpenFAST/r-test/tree/dev/modules/beamdyn/bd_5MW_dynamic>`__.

.. This is a standalone BeamDyn case of the NREL 5MW wind turbine. It simulates 30
.. seconds with a time step size of 0.002 seconds and executes in 24s on NREL's
.. Peregrine supercomputer.

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
