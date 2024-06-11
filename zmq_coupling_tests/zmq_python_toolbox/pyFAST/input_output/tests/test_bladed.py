import os
import numpy as np
import re
import pandas as pd
import unittest
from .helpers_for_test import MyDir, reading_test 

from pyFAST.input_output.bladed_out_file import BladedFile



class Test(unittest.TestCase):
 
    def test_001_read_all(self):
        reading_test('Bladed_out_*.*', BladedFile)   
    

    def test_Bladed(self):
        ## check for binary
        F = BladedFile(os.path.join(MyDir,'Bladed_out_binary.$41')) 
        #F = BladedFile(os.path.join(MyDir,'Bladed_out_binary.$41'))
        DF = F.toDataFrame()
        self.assertAlmostEqual(DF['0.0m-Blade 1 Fx (Root axes) [N]'].values[0],146245.984375)
        self.assertAlmostEqual(DF['0.0m-Blade 1 Fx (Root axes) [N]'].values[-1],156967.484375)

        ## check for ASCII
        F = BladedFile(os.path.join(MyDir,'Bladed_out_ascii.$41'))
        DF = F.toDataFrame()
        self.assertAlmostEqual(DF['0.0m-Blade 1 Fx (Root axes) [N]'].values[0],146363.8)
        self.assertAlmostEqual(DF['0.0m-Blade 1 Fx (Root axes) [N]'].values[-1],156967.22)

    def test_Bladed_case2_project(self):
        F = BladedFile(os.path.join(MyDir,'Bladed_out_binary_case2.$PJ')) 
        DF = F.toDataFrame()
        #print(DFS.keys())
        #DF=DFS['Misc']

        #print(DF.shape)
        #print(DF.columns)
        #print(DF.columns[0])
        #print(DF.columns[50])
        self.assertEqual(DF.shape, (10, 89))
        self.assertEqual(DF.columns[0]  , 'Time [s]')
        self.assertEqual(DF.columns[1]  , 'Time from start of simulation [s]')
        self.assertEqual(DF.columns[27] , '26.41m-DPMOM1 [Nm/m]')
        self.assertEqual(DF.columns[69], '38.75m-Blade 1 y-position [m]')
        self.assertEqual(DF.columns[88], 'Foundation Fz [N]')
        self.assertAlmostEqual(DF['Time from start of simulation [s]'][0]    ,  7.0  ) 
        self.assertAlmostEqual(DF['26.41m-DPMOM1 [Nm/m]'][0]       , -226.85083, 5  ) 
        self.assertAlmostEqual(DF['38.75m-Blade 1 y-position [m]'].values[0], -27.949090957, 5 ) 
        self.assertAlmostEqual(DF['38.75m-Blade 1 y-position [m]'].values[-1], -39.96076965, 5 ) 
        self.assertAlmostEqual(DF['Foundation Fz [N]'][0]             , -1092165.5  ) 
        self.assertAlmostEqual(DF['Foundation Fz [N]'].values[-1]     , -1093664.75  ) 
        
        self.assertFalse(DF.isnull().values.any())

    def test_Bladed_case2_indiv(self):
        F = BladedFile(os.path.join(MyDir,'Bladed_out_binary_case2.$12')) 
        DF=F.toDataFrame()
        self.assertEqual(DF.shape, (10, 14))
        self.assertEqual(DF.columns[0]  , 'Time [s]')
        self.assertEqual(DF.columns[1]  , 'POW2 [W]')
        self.assertAlmostEqual(DF['POW2 [W]'].values[-1] , 1940463.0) 

        F = BladedFile(os.path.join(MyDir,'Bladed_out_binary_case2.$25')) 
        DF=F.toDataFrame()
        self.assertEqual(DF.shape, (10, 17))
        self.assertEqual(DF.columns[0]  , 'Time [s]')
        self.assertEqual(DF.columns[1]  , '-15.0m-MXT [Nm]')
        self.assertAlmostEqual(DF['-15.0m-MXT [Nm]'].values[-1], 1587526.625)

        F = BladedFile(os.path.join(MyDir,'Bladed_out_binary_case2.$69')) 
        DF=F.toDataFrame()
        self.assertEqual(DF.shape, (10, 7))
        self.assertEqual(DF.columns[0]  , 'Time [s]')
        self.assertEqual(DF.columns[1]  , 'Foundation Mx [Nm]')
        self.assertAlmostEqual(DF['Foundation Mx [Nm]'].values[-1], 1587236.375)


        F = BladedFile(os.path.join(MyDir,'Bladed_out_binary_case2.$37')) 
        DF=F.toDataFrame()
        self.assertEqual(DF.shape, (522, 6))
        self.assertEqual(DF.columns[0]  , 'Time [s]')
        self.assertEqual(DF.columns[1]  , 'Simulation Time [s]')
        self.assertAlmostEqual(DF['State with largest error [N]'].values[-1], 9.0)

        # NOTE: this binary file is detected as ascii, and the reading fails..
        F = BladedFile(os.path.join(MyDir,'Bladed_out_binary_case2_fail.$55')) 
        DF=F.toDataFrame()
        self.assertEqual(DF.shape, (50, 1))
        self.assertEqual(DF.columns[0]  , 'Step size histogram [N]')
        #self.assertTrue(np.isnan(DF['Step size histogram [N]'].values[-1]))
        self.assertEqual(DF['Step size histogram [N]'].values[-1],0.0)

        # NOTE: this one is properly dected as binary
        F = BladedFile(os.path.join(MyDir,'Bladed_out_binary_case2.$55')) 
        DF=F.toDataFrame()
        self.assertEqual(DF.shape, (50, 1))
        self.assertEqual(DF.columns[0]  , 'Step size histogram [N]')
        self.assertEqual(DF['Step size histogram [N]'].values[-1], 0)

        F = BladedFile(os.path.join(MyDir,'Bladed_out_binary_case2.$46')) 
        DF=F.toDataFrame()
        self.assertEqual(DF.shape, (10, 13))
        self.assertEqual(DF.columns[1]  , 'Node 1-Water particle velocity in X direction [m/s]')
        self.assertEqual(DF['Node 1-Water particle velocity in X direction [m/s]'].values[-1],-0.25)

        F = BladedFile(os.path.join(MyDir,'Bladed_out_binary_case2.$06')) 
        DF=F.toDataFrame()
        self.assertEqual(DF.shape, (10, 5))
        self.assertEqual(DF.columns[1]  , 'Generator torque [Nm]')
        self.assertEqual(DF['Generator torque [Nm]'].values[-1],12852.1953125)

        F = BladedFile(os.path.join(MyDir,'Bladed_out_binary_case2.$23'))
        DF=F.toDataFrame()
        self.assertEqual(DF.shape, (10, 9))
        self.assertEqual(DF.columns[1]  , 'Stationary hub Mx [Nm]')
        self.assertEqual(DF['Stationary hub Mx [Nm]'].values[-1],1112279.375)


if __name__ == '__main__':
     unittest.main()
     #Test().test_001_read_all()
     #Test().test_Bladed()
     #Test().test_Bladed_case2_project()
     #Test().test_Bladed_case2_indiv()
    
