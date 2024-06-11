""" 
Wrapper around wetb to read/write hawc2 st files.
"""
from .file import File, WrongFormatError

import numpy as np
import pandas as pd
import os

from .wetb.hawc2.st_file import StFile

class HAWC2StFile(File):

    @staticmethod
    def defaultExtensions():
        return ['.st','.dat']

    @staticmethod
    def formatName():
        return 'HAWC2 st file'

    def __init__(self,filename=None, **kwargs):
        self.filename = None
        if filename:
            self.read(filename, **kwargs)

    def _read(self):
        # --- Sanity check read first few lines, check if some start with # and $
        nLinesMax=12
        hasPound=False
        hasDollar=False
        with open(self.filename,'r') as fid:
            for i, line in enumerate(fid):
                hasPound  = hasPound or  line.startswith('#') 
                hasDollar = hasDollar or line.startswith('$')
                if i==nLinesMax:
                    break
        if (not hasPound) or (not hasDollar):
            raise WrongFormatError('The first line of the file dont have `#` or `$`, likely not a st file.')
        #  --- Actual reading
        self.data = StFile(self.filename)

    def _write(self):
        self.data.save(self.filename, precision='%15.07e', encoding='utf-8')

    def __repr__(self):
        s='<{} object> with attribute `data`\n'.format(type(self).__name__)
        return s

    def toDataFrame(self, extraCols=True):
        col_reg=['r_[m]','m_[kg/m]','x_cg_[m]','y_cg_[m]','ri_x_[m]','ri_y_[m]', 'x_sh_[m]','y_sh_[m]','E_[N/m^2]','G_[N/m^2]','I_x_[m^4]','I_y_[m^4]','I_p_[m^4]','k_x_[-]','k_y_[-]','A_[m^2]','pitch_[deg]','x_e_[m]','y_e_[m]']
        col_fpm=['r_[m]','m_[kg/m]','x_cg_[m]','y_cg_[m]','ri_x_[m]','ri_y_[m]','pitch_[deg]','x_e_[m]','y_e_[m]','K11','K12','K13','K14','K15','K16','K22','K23','K24','K25','K26','K33','K34','K35','K36','K44','K45','K46','K55','K56','K66']

        dfs ={}
        nm = len(self.data.main_data_sets)
        for mset in self.data.main_data_sets.keys():
            for iset in self.data.main_data_sets[mset].keys():
                tab = self.data.main_data_sets[mset][iset]
                if tab.shape[1]==19:
                    FPM = False
                    col = col_reg
                else:
                    FPM = True
                    col = col_fpm
                df = pd.DataFrame(data =tab, columns=col )
                if extraCols:
                    df['Ixi_[kg.m]'] = df['ri_x_[m]']**2 * df['m_[kg/m]']
                    df['Iyi_[kg.m]'] = df['ri_y_[m]']**2 * df['m_[kg/m]']
                    df['StrcTwst_[deg]'] = -df['pitch_[deg]']
                    if not FPM:
                        df['EIx_[Nm^2]'] = df['E_[N/m^2]']*df['I_x_[m^4]']
                        df['EIy_[Nm^2]'] = df['E_[N/m^2]']*df['I_y_[m^4]']
                        df['GKt_[Nm^2]'] = df['G_[N/m^2]']*df['I_p_[m^4]']
                        df['EA_[N]']     = df['E_[N/m^2]']*df['A_[m^2]']
                        df['GA_[N]']     = df['G_[N/m^2]']*df['A_[m^2]']
                    df['r_bar_[-]']  = df['r_[m]']/df['r_[m]'].values[-1]

                dfs['{}_{}'.format(mset,iset)] = df


        return dfs
