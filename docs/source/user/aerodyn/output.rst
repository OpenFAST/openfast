.. _ad_output:

Output Files
============

AeroDyn produces three types of output files: an echo file, a summary
file, and a time-series results file. The following sections detail the
purpose and contents of these files.

Echo Files
----------

If you set the ``Echo`` flag to ``TRUE`` in the AeroDyn driver file or the
AeroDyn primary input file, the contents of those files will be echoed
to a file with the naming conventions, *OutFileRoot.ech* for the
driver input file and *OutFileRoot.AD.ech* for the AeroDyn primary
input file. ``OutFileRoot`` is either specified in the I/O SETTINGS
section of the driver input file when running AeroDyn standalone, or by
the FAST program when running a coupled simulation. The echo files are
helpful for debugging your input files. The contents of an echo file
will be truncated if AeroDyn encounters an error while parsing an input
file. The error usually corresponds to the line after the last
successfully echoed line.

.. _sec:ad_SumFile:

Summary File
------------

AeroDyn generates a summary file with the naming convention,
*OutFileRoot.AD.sum* if the ``SumPrint`` parameter is set to ``TRUE``.
``OutFileRoot`` is either specified in the I/O SETTINGS section of the
driver input file when running AeroDyn standalone, or by the FAST
program when running a coupled simulation. This file summarizes key
information about your aerodynamics model, including which features have
been enabled and what outputs have been selected.

Results Files
-------------

In standalone mode, the AeroDyn time-series results (a separate file for
each case) are written to text-based files with the naming convention
*OutFileRoot.#.out*, where ``OutFileRoot`` is specified in the I/O
SETTINGS section of the driver input file and the ‘\ *#*\ ’ character is
an integer number corresponding to a test case line found in the
COMBINED-CASE ANALYSIS section. If AeroDyn is coupled to FAST, then FAST
will generate a master results file that includes the AeroDyn results
and AeroDyn will not write out its own results. The results are in table
format, where each column is a data channel (the first column always
being the simulation time), and each row corresponds to a simulation
output time step. The data channels are specified in the OUTPUTS section
of the AeroDyn primary input file. The column format of the
AeroDyn-generated files is specified using the ``OutFmt`` parameter of
the driver input file.
