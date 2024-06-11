import glob
import os
import copy
from wetb.hawc2.hawc2_pbs_file import JESS_WINE32_HAWC2MB
from wetb.hawc2.htc_file import HTCFile
from wetb.utils.cluster_tools.pbsfile import PBSMultiRunner


class HTCFileSet():
    def __init__(self, model_path, htc_lst="**/*.htc"):
        self.model_path = model_path

        if not isinstance(htc_lst, list):
            htc_lst = [htc_lst]

        self.htc_files = []
        for htc_path in htc_lst:
            if os.path.isfile(htc_path):
                self.htc_files.append(htc_path)
            else:
                if not os.path.isabs(htc_path):
                    htc_path = os.path.join(model_path, htc_path)
                for filename in glob.iglob(htc_path, recursive=True):
                    self.htc_files.append(filename)

    def pbs_files(self, hawc2_path, hawc2_cmd, queue='workq', walltime=None,
                  input_files=None, output_files=None, copy_turb=(True, True)):

        return (HTCFile(htc).pbs_file(hawc2_path, hawc2_cmd, queue=queue, walltime=walltime,
                                      input_files=copy.copy(input_files),
                                      output_files=copy.copy(output_files),
                                      copy_turb=copy_turb) for htc in self.htc_files)

    def save_pbs_files(self, hawc2_path=None, hawc2_cmd=JESS_WINE32_HAWC2MB, queue='workq', walltime=None,
                       input_files=None, output_files=None, copy_turb=(True, True)):
        for pbs in self.pbs_files(hawc2_path, hawc2_cmd, queue=queue, walltime=walltime,
                                  input_files=input_files, output_files=output_files,
                                  copy_turb=copy_turb):
            pbs.save(self.model_path)


if __name__ == '__main__':
    #model_path = r'R:\HAWC2_tests\v12.6_mmpe3\win32\simple1'
    model_path = "w:/simple1"
    pbs_files = HTCFileSet(model_path).pbs_files(
        hawc2_path=r"R:\HAWC2_tests\v12.6_mmpe3\hawc2\win32", hawc2_cmd=JESS_WINE32_HAWC2MB, input_files=['data/*'])
    import pandas as pd
    time_overview = pd.read_excel(
        r'C:\mmpe\programming\Fortran\HAWC2_git\HAWC2\pytest_hawc2\release_tests\Time_overview.xlsx')
    for pbs in pbs_files:
        f = pbs.filename

        pbs.walltime = time_overview.loc[f[:-3].replace("pbs_in/", 'simple1/')]['mean'] * 24 * 3600
        pbs.save(model_path)
    PBSMultiRunner(model_path, nodes=1, ppn=10).save()
