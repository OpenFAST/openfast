from pathlib import Path

import bs4
import copy

formatter = bs4.formatter.HTMLFormatter(indent=4)

options_debug_release = {"Debug": {}, "Release": {}}

cfg_names = [
    "Debug|x64",
    "Debug_Double|x64",
    "Debug_Matlab|x64",
    "Release|x64",
    "Release_Double|x64",
    "Release_Matlab|x64",
    "Release_OpenMP|x64",
    "Release_Double_OpenMP|x64",
]

for path in Path(".").rglob("*.vfproj"):

    print(path)
    with open(path) as fp:
        soup = bs4.BeautifulSoup(fp, "xml")
    cfgs = soup.find("Configurations")
    cfg_map = {
        "Debug|x64": cfgs.find("Configuration", Name="Debug|x64"),
        "Release|x64": cfgs.find("Configuration", Name="Release|x64"),
    }
    cfgs.clear()
    for cfg_name in cfg_names:
        if "Debug" in cfg_name:
            cfg = copy.copy(cfg_map["Debug|x64"])
        else:
            cfg = copy.copy(cfg_map["Release|x64"])
        cfg["Name"] = cfg_name

        # Get tool elements
        compiler_tool = cfg.find("Tool", Name="VFFortranCompilerTool")
        linker_tool = cfg.find("Tool", Name="VFLinkerTool")
        prebuild_tool = cfg.find("Tool", Name="VFPreBuildEventTool")

        # Compiler tool settings
        compiler_tool["Preprocess"] = "preprocessYes"
        compiler_tool["MultiProcessorCompilation"] = "true"
        compiler_tool["UseMkl"] = "mklSequential"
        compiler_tool["WarnUnusedVariables"] = "false"
        if "Debug" in cfg["Name"]:
            compiler_tool["RuntimeLibrary"] = "rtMultiThreadedDebug"
        else:
            compiler_tool["RuntimeLibrary"] = "rtMultiThreaded"

        # Determine project type (static lib, shared lib, executable)
        if cfg.attrs.get("ConfigurationType", "") == "typeStaticLibrary":
            cfg["OutputDirectory"] = "..\\..\\build\\lib"
        elif cfg.attrs.get("ConfigurationType", "") == "typeDynamicLibrary":
            cfg["OutputDirectory"] = "..\\..\\build\\bin"
            if 'Debug' in cfg_name:
                compiler_tool["FloatingPointExceptionHandling"] = "fpe0"
            linker_tool["StackReserveSize"] = "9999999"
        elif linker_tool != None and linker_tool["SubSystem"] == "subSystemConsole":
            cfg["OutputDirectory"] = "..\\..\\build\\bin"
            if 'Debug' in cfg_name:
                compiler_tool["FloatingPointExceptionHandling"] = "fpe0"
            linker_tool["StackReserveSize"] = "9999999"
            linker_tool['GenerateManifest'] = "false"
        else:
            print("unknown project type")
            continue

        # Set intermediate build directory
        cfg["IntermediateDirectory"] = (
            "..\\..\\build\\$(Configuration)_$(Platform)\\$(ProjectName)\\"
        )

        # Preprocessor defines
        defines = []

        # Project specific settings
        if "NWTC" in str(path):
            # defines.append("HAS_FORTRAN2008_FEATURES")
            pass
        if "VersionInfo" in str(path):
            defines.append("GIT_INCLUDE_FILE='..\\gitVersionInfo.h'")
            prebuild_tool["CommandLine"] = "..\\CreateGitVersion.bat"

        # Configuration spectific settings
        if "Double" in cfg["Name"]:
            compiler_tool["RealKIND"] = "realKIND8"
            compiler_tool["DoublePrecisionKIND"] = "doublePrecisionKIND8"
            if "NWTC" in str(path):
                defines.append("OPENFAST_DOUBLE_PRECISION")
        if "OpenMP" in cfg["Name"]:
            compiler_tool["OpenMP"] = "OpenMPParallelCode"
            compiler_tool["EnableOpenMPSupport"] = "OpenMPParallelCodeIFX"
        if "Matlab" in cfg["Name"]:
            defines.append("COMPILE_SIMULINK")
            defines.append("CONSOLE_FILE")

        # Preprocessor defines
        compiler_tool["PreprocessorDefinitions"] = ";".join(defines)

        # Add config to configs
        cfgs.append(cfg)

    # Update registry file configurations
    for f in soup.find_all("File"):
        fcs = f.find_all("FileConfiguration")
        if len(fcs) == 0:
            continue
        fc_base = copy.copy(fcs[0])
        for fc in f.find_all("FileConfiguration"):
            fc.decompose()
        for cfg_name in cfg_names:
            fc = copy.copy(fc_base)
            fc["Name"] = cfg_name
            f.append(fc)

    # Write file
    with open(path, "w") as fp:
        for line in soup.prettify().splitlines():
            try:
                n = line.index("<")
            except:
                n = 0
            line = ("\t" * n) + line[n:] + "\n"
            fp.write(line)
