# OpenFAST release process

## Prep
### pull request
1. Create release changelog.md
2. Post PR with contents of the changelog/release notes
3. Get reviews and address all issues
4. Add changelog for the release -- include all PR's that will be included

### PR branch updates
1. checkout branch to merge in and verify builds on VS if any changes for VS or new files added
2. Update the documentation version in docs/conf.py
3. Update the versions in docs/source/user/api_change.rst
4. Verify readthedocs builds correctly
5. Update `openfast_io/pyproject.toml`
6. Update `glue-codes/python/pyproject.toml` (for `pyOpenFAST`)

****

## Posting the release
### r-test
1. Merge and add annotated tag
2. Update pointer on main OF repository (not always necessary)

### Main repository
1. Merge PR
2. Create release with new tag
   * Copy `Changelog` section down from the changlog.md file to release notes
   * add short intro section at top with 2 sentence synopsis (see prior release for this)
   * copy `Precompiled Windows Binaries` section from prior release, and update as needed into the release notes
   * check the `create discussion` box
   * Post
3. delete `rc-` branch if merging from one

### Windows executables build and upload
After posting and tagging release
1. Pull main and tags
   * `git fetch --tags OpenFAST`
   * `git fetch OpenFAST main:main`
   * `git checkout main`
2. Delete `vs-build` and checkout again
   * `rm -rf vs-build`
   * `git checkout vs-build`
3. Set a couple of VS files to not track changes on files that VS wants to update Windows related stuff in
   ```
   git update-index --assume-unchanged vs-build/MAPlib/MAP_dll.vcxproj vs-uild/Registry/FAST_Registry.vcxproj
   ```

4. Compile executables for Windows builds (manual process - use GH actions `deploy` if possible)
   * Run one of the executables and check the version info. Muck about with VS if there is an issue.
   * Also run `dumpbin.exe /dependents <exe>.exe` to check static linking
   * NOTE: build the simulink last -- it messes up some things otherwise
    - [ ] `AeroDisk_Driver_x64.exe`
    - [ ] `AeroDyn_Driver_x64.exe`
    - [ ] `AeroDyn_Driver_x64_OpenMP.exe`
    - [ ] `AeroDyn_Inflow_c_binding_x64.dll`
    - [ ] `AeroDyn_Inflow_c_binding_x64_OpenMP.dll`
    - [ ] `BeamDyn_Driver_x64.exe`
    - [ ] `DISCON.dll` (x64)
    - [ ] `DISCON_ITIBarge.dll` (x64)
    - [ ] `DISCON_OC3Hywind.dll` (x64)
    - [ ] `FAST.Farm_x64.exe`
    - [ ] `FAST.Farm_x64_OMP.exe`
    - [ ] `FAST_SFunc.mexw64` -- build from MATLAB
    - [ ] `HydroDynDriver_x64.exe`
    - [ ] `HydroDyn_c_binding_x64.dll`
    - [ ] `InflowWind_c_binding_x64.dll`
    - [ ] `InflowWind_Driver_x64.exe`
    - [ ] `InflowWind_Driver_x64_OpenMP.exe`
    - [ ] `MoorDyn_Driver_x64.exe`
    - [ ] `MoorDyn_c_binding_x64.dll`
    - [ ] `OpenFAST-Simulink_x64.dll` -- change `additional dependencies` in the `OpenFAST-Simulink` project in `FAST` to point to correct install of MATLAB
    - [ ] `openfast_x64.exe`
    - [ ] `SeaStateDriver_x64.exe`
    - [ ] `SeaState_c_binding_x64.dll`
    - [ ] `SimpleElastoDyn_x64.exe`
    - [ ] `SubDyn_x64.exe`
    - [ ] `Turbsim_x64.exe`
    - [ ] `UnsteadyAero_x64.exe`
5. Upload all filesUnset the no tracking of files
   ```
   git ls-files -v | grep "^[a-z]"
   git update-index --no-assume-unchanged <files-from-cmd-above>
   ```
## Post-release
### Docker Image push to ghcr.io
1. Build latest `OpenFAST/main` image locally (GH actions fails due to memory usage)
2. Push image to ghcr.io/openfast/openfast using tags `latest` and `<version>`

