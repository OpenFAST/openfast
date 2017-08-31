

Add fast v8 to openfast transition info here. the link is already created in the main page side menu



The rst code below is just an example of the various formatting stuff

link: 
`r-test <https://github.com/openfast/r-test>`__

header (something like <h2>): 
Testing OpenFAST
================

smaller header (maybe <h4>?)
Dependencies
------------

code:
``reg_tests``

list:
Dependencies
------------
- Python 3+
- Numpy
- CMake and CTest (Optional in regression test)
- matplotlib (Optional in regression test)

list of links... table of contents style
Test specific documentation
---------------------------
.. toctree::
   :maxdepth: 2

   Unit Test <test/unit_test.rst>
   Regression Test <test/regression_test.rst>
   