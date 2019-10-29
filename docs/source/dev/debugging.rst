.. _debugging:

Debugging OpenFAST
==================

Being a Fortran project, OpenFAST can be challenging to debug and the process
is unique for each system and environment. However, a common requirement for
all systems is to compile OpenFAST in debug mode.

Keep in mind that some OpenFAST cases can be quite large in their memory
footprint and may take a long time to reach the point you're targetting in
the code. Choosing a minimal test case could save significant time.

It may by helpful to write a small fortran program to verify that your
debugging tools are set up properly before diving in to OpenFAST. Be sure to
simulate a bug by doing something like accessing an array element that is not
allocated and verify that you can catch the bug with your tools.

Debugging on Windows
--------------------
Windows developers using Intel tools can use Visual Studio for debugging. This
is a straightforward process with lots of support from Intel.

Otherwise, Windows developers compiling with CygWin or MinGW should proceed to
the section for debugging with linux.

Debugging on Linux and macOS
----------------------------
First, compile OpenFAST in debug mode by setting CMAKE_BUILD_TYPE to Debug.
You can do this on the command line with

.. code-block:: bash

    cmake .. -D CMAKE_BUILD_TYPE=Debug

or by using ccmake to open the command line cmake gui to change it.

The GNU debugger, GDB, works well for debugging compiled code. It has a
comprehensive command line interface which enables developers to add
breakpoints and inspect variables.

Driving the debugger through an IDE can make inspecting the code much more
efficient. One IDE known to work well is Visual Studio Code with the Native
Debug extension. You can set up a launch configuration in VS Code so that
you can debug a particular OpenFAST case through the IDE. To do this, open
the launch configuration and add a block similar to this:

.. code-block:: json

        {
            "name": "AOC_WSt",
            "type": "gdb",
            "request": "launch",
            "printCalls": false,
            "showDevDebugOutput": false,
            "valuesFormatting": "prettyPrinters",
            "gdbpath": "gdb",
            "target": "${workspaceRoot}/build/glue-codes/openfast/openfast",
            "cwd": "${workspaceRoot}/build/reg_tests/glue-codes/openfast/AOC_WSt/",
            "arguments": "${workspaceRoot}/build/reg_tests/glue-codes/openfast/AOC_WSt/AOC_WSt.fst",
        },

macOS configuration
~~~~~~~~~~~~~~~~~~~
GDB on macOS needs some configuration before the system allows it to take
over a process. It is recommended that gdb be installed with homebrew

.. code-block:: bash

    brew info gdb
    brew install gdb

After that completes, be sure to follow the caveats to finish the installation.
For gdb 8.2.1, it looks like this:

.. code-block:: bash

    ==> Caveats
    gdb requires special privileges to access Mach ports.
    You will need to codesign the binary. For instructions, see:

    https://sourceware.org/gdb/wiki/BuildingOnDarwin

    On 10.12 (Sierra) or later with SIP, you need to run this:

    echo "set startup-with-shell off" >> ~/.gdbinit

For Native Debug on macOS, you have to sort of hack the extension to allow
breakpoints in fortran files by adding this line to your settings.json:

.. code-block:: json

    "debug.allowBreakpointsEverywhere": true,
