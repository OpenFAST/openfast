.. _code_style:

Code Style
~~~~~~~~~~
OpenFAST and its underlying modules are mostly written in Fortran adhering to
the 2003 standard, but modules can be written in C or C++. The
`NWTC Programmer's Handbook <https://nwtc.nrel.gov/system/files/ProgrammingHandbook_Mod20130717.pdf>`__
is the definitive reference for all questions related to working with the
FAST Framework and adding code to OpenFAST.

Generally, code should be written such that it is straightforward to read.
Syntactic sugar or brevity should not detract from readability. The exception
to this is in situations where performance requires poorly readable code.
Here, comment blocks should be used to describe what is not readily apparent
in the code. Indentation is typically three spaces and no tabs.
