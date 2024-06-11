import unittest
import os
import numpy as np
from pyFAST.input_output.tests.helpers_for_test import MyDir, reading_test 
from pyFAST.input_output.flex_profile_file import FLEXProfileFile
from pyFAST.input_output.flex_blade_file import FLEXBladeFile 
from pyFAST.input_output.flex_wavekin_file import FLEXWaveKinFile 
from pyFAST.input_output.flex_doc_file import FLEXDocFile 

    
import pandas as pd

class Test(unittest.TestCase):
 
    #def test_001_read_all(self, DEBUG=True):
    #    reading_test('FLEX*.*', weio.read)

    def DF(self,FN):
        """ Reads a file with weio and return a dataframe """ 
        pass
        #return Flex(os.path.join(MyDir,FN)).toDataFrame()

    def test_FLEXProfiles(self):
        df = FLEXProfileFile(os.path.join(MyDir,'FLEXProfile.pro')).toDataFrame()
        self.assertAlmostEqual(df['pc_set_2_t_57.0'].values[2,2],0.22711022)

    def test_FLEXBlade(self):
        Bld=FLEXBladeFile(os.path.join(MyDir,'FLEXBlade002.bld')).toDataFrame()
        self.assertAlmostEqual(Bld['r_[m]'].values[-1],61.5)
        self.assertAlmostEqual(Bld['Mass_[kg/m]'].values[-1],10.9)
        self.assertAlmostEqual(Bld['Chord_[m]'].values[3],3.979815059)

    def test_FLEXWaves(self):
        wk = FLEXWaveKinFile(os.path.join(MyDir, 'FLEXWaveKin.wko'))
        self.assertEqual(wk['MaxLongiVel'],2.064)
        self.assertEqual(wk['Tp']         ,12.54)
        self.assertEqual(len(wk['RelDepth']),12)
        self.assertEqual(wk['data']['Time_[s]'].values[-1],3.0)
        self.assertEqual(wk['data']['a_z=20.0_x=0.0_[m/s^2]'].values[-1],0.06)

    def test_FLEXDoc(self):
        doc = FLEXDocFile(os.path.join(MyDir, 'FLEXDocFile.out'))
        self.assertAlmostEqual(doc['RNA']['Mass'], 2.85e-6)
        self.assertAlmostEqual(doc['Tower']['Length'], 1.0)
        self.assertAlmostEqual(doc['Tower']['SectionData'].shape[0], 11)
        self.assertAlmostEqual(doc['Tower']['SectionData'].shape[1], 9)
        self.assertAlmostEqual(doc['Tower']['ShapeFunction_DOF1_Shape'].shape[0], 12)
        self.assertAlmostEqual(doc['Foundation']['Mass'], 900000)
        self.assertAlmostEqual(doc['Foundation']['ShapeFunction_DOF1_Shape'].shape[0], 101)
        self.assertAlmostEqual(doc['Foundation']['ShapeFunction_DOF1_Shape']['H_[m]'].values[-1], 100)
        self.assertAlmostEqual(doc['Foundation']['ShapeFunction_DOF1_Shape']['U_[m]'].values[-1], 1)
        self.assertAlmostEqual(doc['Foundation']['ShapeFunction_DOF2_Shape']['U_[m]'].values[-1], 0)
        self.assertEqual(type(doc['Blade']['ShapeFunction_DOF1_Shape']) is pd.DataFrame, True)


if __name__ == '__main__':
#     Test().test_FLEXWaves()
#     Test().test_FLEXDoc()
    unittest.main()
