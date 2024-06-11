import unittest
import os
import numpy as np
from shutil import copyfile
import pyFAST.converters.beamdyn as bd

class Test(unittest.TestCase):

    def test_BAR0(self):

        # Get current directory so this script can be called from any location
        MyDir=os.path.dirname(__file__)
        # --- BeamDynToHawc2
        htc_template   = os.path.join(MyDir,'../../../data/templates/hawc2_template.htc') # readonly
        BD_mainfile    = os.path.join(MyDir,'../../../data/BAR0/BAR0_BeamDyn.dat')      # readonly
        BD_bladefile   = os.path.join(MyDir,'../../../data/BAR0/BAR0_BeamDyn_Blade.dat') # readonly
        H2_htcfile_new = '_BAR0.htc' # will be created
        H2_stfile      = '_BAR0_st.st' # will be created
        copyfile(htc_template,  H2_htcfile_new) # Backup template
        df_c2, df_st = bd.beamDynToHawc2(BD_mainfile, BD_bladefile, H2_htcfile_new, H2_stfile, 'blade1', FPM=True, verbose=True)

        # NOTE: NREL5MW is too simple of a tests since straight
        np.testing.assert_almost_equal(df_c2['x_[m]'].values[10]      ,  0.0   ,3)
        np.testing.assert_almost_equal(df_c2['y_[m]'].values[-1]      , -4.0    ,3)
        np.testing.assert_almost_equal(df_c2['z_[m]'].values[-1]      ,99.9963   ,3)
        np.testing.assert_almost_equal(df_c2['twist_[deg]'].values[-1],2.57534   ,3)
        np.testing.assert_almost_equal(df_c2['twist_[deg]'].values[0] ,-20.0020  ,3)


        np.testing.assert_almost_equal(df_st['m_[kg/m]'].values[10]   ,694.6133  ,3)

        np.testing.assert_almost_equal(df_st['x_cg_[m]'].values[10]   ,-0.74337  ,3)
        np.testing.assert_almost_equal(df_st['y_cg_[m]'].values[10]   , 0.01759  ,3)

        np.testing.assert_almost_equal(df_st['x_e_[m]'].values[10]    ,-0.75435  ,3)
        np.testing.assert_almost_equal(df_st['y_e_[m]'].values[10]    , 0.01374  ,3)

        np.testing.assert_almost_equal(df_st['pitch_[deg]'].values[10],-1.48000  ,3)

        np.testing.assert_almost_equal(df_st['ri_x_[m]'].values[10]   , 0.57088  ,3)
        np.testing.assert_almost_equal(df_st['ri_y_[m]'].values[10]   , 1.4513  ,3)

        np.testing.assert_almost_equal(df_st['K11'].values[10]/1e8    ,5.694503  ,3)
        np.testing.assert_almost_equal(df_st['K22'].values[10]/1e8    ,2.467721  ,3)
        np.testing.assert_almost_equal(df_st['K33'].values[10]/1e8    ,125.1373  ,3)
        np.testing.assert_almost_equal(df_st['K44'].values[10]/1e8    ,44.0502   ,3)
        np.testing.assert_almost_equal(df_st['K55'].values[10]/1e8    ,274.3825  ,3)
        np.testing.assert_almost_equal(df_st['K66'].values[10]/1e8    ,8.232813  ,3)

        np.testing.assert_almost_equal(df_st['K12'].values[10]/1e6    , 1.20013  ,3)
        np.testing.assert_almost_equal(df_st['K16'].values[10]/1e6    ,-33.9432  ,3)
        np.testing.assert_almost_equal(df_st['K26'].values[10]/1e6    ,211.6052  ,3)

        os.remove(H2_htcfile_new)
        os.remove(H2_stfile)


if __name__ == '__main__':

    unittest.main()
