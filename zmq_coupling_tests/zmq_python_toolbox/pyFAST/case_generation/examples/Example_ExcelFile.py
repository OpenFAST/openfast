"""
Generate, run and postprocess openfast cases using an Excelfile which defines the parameters to change for each simulation
"""
import os
import pandas as pd

import pyFAST.case_generation.case_gen as case_gen
import pyFAST.case_generation.runner as runner
import pyFAST.input_output.postpro as postpro
import pyFAST.input_output as io

# Get current directory so this script can be called from any location
scriptDir=os.path.dirname(__file__)


def main(run=True):
    # --- Main Parameters
    ref_dir        = os.path.join(scriptDir, '../../../data/NREL5MW/')  # Folder where the fast input files are located (will be copied)
    FAST_EXE       = os.path.join(scriptDir, '../../../data/openfast.exe') # Location of a FAST exe (and dll)
    main_file      = 'Main_Onshore.fst'          # Main file in ref_dir, used as a template
    work_dir       = os.path.join(scriptDir, '_NREL5MW_ParametricExcel/')     # Output folder (will be created)
    parametricFile = os.path.join(scriptDir, 'ParametricExcel.xlsx')         # Excel file containing set of parameters

    # --- Reading Excel file, converting it to a list of dictionaries, and generate input files
    dfs    = io.excel_file.ExcelFile(parametricFile).toDataFrame()
    if isinstance(dfs, dict):
        df = dfs[list(dfs.keys())[0]]
    else:
        df = dfs
    PARAMS = df.to_dict('records')
    fastFiles=case_gen.templateReplace(PARAMS, ref_dir, outputDir=work_dir, removeRefSubFiles=True, removeAllowed=False, main_file=main_file)

    if run:
        # --- Running fast simulations
        print('>>> Running {} simulations in {} ...'.format(len(fastFiles), work_dir))
        runner.writeBatch(os.path.join(work_dir,'_RUN_ALL.bat'),fastFiles,fastExe=FAST_EXE)
        runner.run_fastfiles(fastFiles, showOutputs=True, fastExe=FAST_EXE, nCores=4)

        # --- Postpro - Computing averages at the end of the simluation
        print('>>> Postprocessing...')
        outFiles = [os.path.splitext(f)[0]+'.outb' for f in fastFiles]
        ColKeepStats  = ['RotSpeed_[rpm]','BldPitch1_[deg]','RtAeroCp_[-]','RtAeroCt_[-]','Wind1VelX_[m/s]']
        result = postpro.averagePostPro(outFiles,avgMethod='constantwindow',avgParam=5,ColKeep=ColKeepStats,ColSort='RotSpeed_[rpm]')
        csv_file = '_ParametricExcel_Summary.csv'
        result.to_csv(csv_file ,sep='\t',index=False)
        print('Average values saved to _ParametricExcel_Summary.csv')
    return work_dir

if __name__=='__main__':
    main() 

if __name__=='__test__':
    work_dir = main(run=False) # Need openfast.exe, not running
    import shutil
    shutil.rmtree(work_dir)


