.. _regression_test_windows:

Windows with Visual Studio regression test
==========================================

1) Clone the openfast repo and initialize the testing database

    a) Open a git command shell window (like git bash)

    b) Change your working directory to the location above where you want your local repo to be located (the repo will be placed into a folder called openfast at this location)

    c. Type:  ``git clone https://github.com/openfast/openfast.git`` (this creates a local version of the openfast repo on your computer).
    You should see something like this:

    :: 

          Cloning into 'openfast'...
          remote: Counting objects: 23801, done.
          remote: Compressing objects: 100% (80/80), done.
          remote: Total 23801 (delta 73), reused 102 (delta 50), pack-reused 23670
          Receiving objects: 100% (23801/23801), 92.10 MiB  18.99 MiB/s, done.
          Resolving deltas: 100% (13328/13328), done.
          Checking connectivity... done.


    d) Type: ``cd openfast``  (change your working directory to the openfast folder)

    e) Type: ``git checkout dev`` (this places your local repo on the correct branch of the openfast repo)

    f) Type: ``git submodule update --init --recursive`` (this downloads the testing database to your computer)
       You should see something like this:

    ::

          Submodule 'reg_tests/r-test' (https://github.com/openfast/r-test.git) registered for path 'reg_tests/r-test'
          Cloning into 'reg_tests/r-test'...
          remote: Counting objects: 3608, done.
          remote: Compressing objects: 100% (121/121), done.
          remote: Total 3608 (delta 22), reused 161 (delta 21), pack-reused 3442
          Receiving objects: 100% (3608/3608), 154.52 MiB  26.29 MiB/s, done.
          Resolving deltas: 100% (2578/2578), done.
          Checking connectivity... done.
          Submodule path 'reg_tests/r-test': checked out 'b808f1f3c1331fe5d03c5aaa4167532c2492d378'


2) Build The Regression Testing DISCON DLLs

    a) Open the Visual Studio Solution (``Discon.sln``) located in ``openfast\vs-build\Discon`` folder

    b) Choose Release and x64 for the Solutions Configuration and Solutions Platform, respectively

    c) From the menu bar select ``Build->Build Solution``

    d) You should now see the files ``Discon.dll``, ``Discon_ITIBarge.dll``, and ``Discon_OC3Hywind.dll`` in your ``openfast\reg_tests\r-test\glue-codes\fast\5MW_Baseline\ServoData`` folder.

3) Build OpenFAST using Visual Studio

    a) Open the Visual Studio Solution (``FAST.sln``) located in ``openfast\vs-build\FAST`` folder

    b) Choose Release_Double and x64 for the Solutions Configuration and Solutions Platform, respectively

    c) From the menu bar select ``Build->Build Solution``

       i)  If this is the first time you have tried to build openfast, you will get build errors!!! [continue to steps (ii) and (iii), otherwise if FAST builds successfully, continue to step (3d) ]

       ii) Cancel build using the menubar ``Build->Cancel``
            [ VS is confused about the build-order/dependency of the project files in FASTlib., but canceling and restarting VS, it somehow as enough info from the partial build to get this right, now]

       iii) Close your Visual Studio and then Repeat Steps (a) through (c)

    d) You should now see the file ``openfast_x64_Double.exe`` in your ``openfast\build\bin`` folder.


4) Prepare regression tests

    a) Create a subdirectory called ``reg_tests`` in your ``openfast\build`` folder.

    b) Copy the contents of ``openfast\reg_tests\r-test`` to ``openfast\build\reg_tests``.


5) Execute the OpenFAST regression Tests

    a) Open a command prompt which is configured for Python [ like Anaconda3 ]
 
    b) Change your working directory to ``openfast\reg_tests``

    c) Type: ``python manualRegressionTest.py ..\build\bin\openfast_x64_Double.exe 2.0 1.9`` 
         You should see this: ``executing AWT_YFix_WSt``

    d) The tests will continue to execute one-by-one until you finally see something like this:

    ::

      executing AWT_YFix_WSt                           PASS
      executing AWT_WSt_StartUp_HighSpShutDown         PASS
      executing AWT_YFree_WSt                          PASS
      executing AWT_YFree_WTurb                        PASS
      executing AWT_WSt_StartUpShutDown                PASS
      executing AOC_WSt                                PASS
      executing AOC_YFree_WTurb                        PASS
      executing AOC_YFix_WSt                           PASS
      executing UAE_Dnwind_YRamp_WSt                   PASS
      executing UAE_Upwind_Rigid_WRamp_PwrCurve        PASS
      executing WP_VSP_WTurb_PitchFail                 PASS
      executing WP_VSP_ECD                             PASS
      executing WP_VSP_WTurb                           PASS
      executing WP_Stationary_Linear                   PASS
      executing SWRT_YFree_VS_EDG01                    PASS
      executing SWRT_YFree_VS_EDC01                    PASS
      executing SWRT_YFree_VS_WTurb                    PASS
      executing 5MW_Land_DLL_WTurb                     PASS
      executing 5MW_OC3Mnpl_DLL_WTurb_WavesIrr         PASS
      executing 5MW_OC3Trpd_DLL_WSt_WavesReg           PASS
      executing 5MW_OC4Jckt_DLL_WTurb_WavesIrr_MGrowth PASS
      executing 5MW_ITIBarge_DLL_WTurb_WavesIrr        PASS
      executing 5MW_TLP_DLL_WTurb_WavesIrr_WavesMulti  PASS
      executing 5MW_OC3Spar_DLL_WTurb_WavesIrr         PASS
      executing 5MW_OC4Semi_WSt_WavesWN                PASS
      executing 5MW_Land_BD_DLL_WTurb                  PASS
      
    e) If an individual test succeeds you will see ``PASS`` otherwise you will see ``FAIL`` after that test's name
