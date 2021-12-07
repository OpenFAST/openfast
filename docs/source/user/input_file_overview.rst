.. _input_file_overview:


Input file formats
==================

OpenFAST uses two primary input file formats: *value column* where the first
value on the line is read, and *key+value* where a value and keyword pair are
read. Both formats are line number based where a specific input is expected on a
specific line, with some exceptions. 

.. _sec_value_column:

Value column input files
------------------------

Only the first column in a *value column* based input file is read.  This is the
historical format used by OpenFAST and its predecessors (the keyword was often
referenced in the source code and documentation, but OpenFAST did not process
the keyword or description). Everything after the
first value read is simply ignored by the code.  This allowed the user to keep
old values while modifying things.  So for example, and input line like

::

 2       20   TMax            - Total run time (s)

would be read as `2` and the `20` and everything after it ignored.

This format and associated parsing methodology is somewhat limited in informing
the user of errors in parsing, and limited the ability to pass entire inpute
files as text strings from another code (such as a Python driver code).


.. _sec_format_key_value:

Key + Value input files
-----------------------

The first two columns are read in *key + value* input files.  One of these two
columns must contain the **exact** keyword, and the other must contain the value
that corresponds to it.  For example, an input line 

::

         20   TMax            - Total run time (s)

is equivalent to

::

   TMax         20            - Total run time (s)

One additional feature of this input file format is the ability to add an
arbitrary number of comment lines wherever the user wishes.  Any line starting
with `!`, `#`, or `%` will be treated as a comment line and ignored.  For
example,


::

  ! This is a comment line that will be skipped
          %  and this is also a comment line that will be skipped
  # as is this comment line
         20   TMax            - Total run time (s)
  ! the first two columns in the above line will be read as the value + key

The parser for this format of input file also tracks which lines were comments,
and which lines contained the value and key pair.  If a keyname is not found the
parser will return an error with information about which line it was reading
from.


Modules using Key + Value Format
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following modules use the *key + value* format input files (all other
modules use the *value column* format):

============== ==========================================================
 Module         Input file           
============== ==========================================================
AeroDyn         Main AD15 input file 
AeroDyn         Airfoil files
HydroDyn        Main HD input file
InflowWind      Main IfW input file
InflowWind      Uniform wind input file
InflowWind      Bladed wind summary file
ServoDyn        Main ServoDyn input file
ServoDyn        Structural control submodule input file
ServoDyn        Structural control sumbodule prescribed force input file
SubDyn          SubDyn SSI matrix input files
============== ==========================================================

Note that key + value format and value column input files can be identical
if the value is stated before the key.

Reasons for change
~~~~~~~~~~~~~~~~~~

The main reason for the change in the input file parsing was to allow for the
passing of a complete input file in memory from a wrapper code into OpenFAST or
a module.  For example, when including the AeroDyn module into a Python code,
the input file can be passed in directly in memory without writing to disk
first.  This helps reduce the IO overhead in optimization loops where the module
might be called many times sequentially with very small changes to the input
file.  *NOTE: this is still a work in progress, so not all modules can be linked
this way yet*.

To accomplish this, the file parser written by Marshall Buhl for parsing airfoil
tables in AeroDyn 15 in FAST8 was used.  This parser included the more robust
*key + value* input format.



.. _sec_troubleshoot_input_file:

Troubleshooting input files
---------------------------

When troubleshooting an input file error, try the following procedure:

1. An error message containing a line number and variable name, the file format
   being parsed is a *key + value* format.  Check that the key is spelled
   exactly as the input file.  See :numref:`sec_format_key_value` below.
2. An error message containing only the variable name but no line number is a
   *value column* input file format.  See :numref:`sec_value_column` below.
3. Turn on `echo` option in the input file and check the resulting `.ech` for
   which line the file parsing stopped at. This may help isolate where the input
   file parsing failed when no line number is given in the error message.
4. Compare the problematic input file with an input file of the same type from
   the regression test suite distributed with OpenFAST.  See section
   :numref:`testing` for details on the regression tests, or check the
   repository at `r-test <https://github.com/openfast/r-test>`__ .


..
   Input file type by module
   -------------------------
   ============== ====================== =====================
    Module         Input file             Type
   ============== ====================== =====================
   OpenFAST        Main .fst input file   Value column
   OpenFAST        Matlab mode shape      Value column
   OpenFAST        Mode shape             Value column
   OpenFAST        Checkpoint file        Binary
   ============== ====================== =====================
