import unittest
import os
import numpy as np
from pyFAST.input_output.tests.helpers_for_test import MyDir, reading_test 
from pyFAST.input_output import FASTSummaryFile

class Test(unittest.TestCase):

    def test_001_read_all(self, DEBUG=True):
        reading_test('FASTSum*.*', FASTSummaryFile)

    def test_FASTSum(self):
        f = FASTSummaryFile(os.path.join(MyDir, 'FASTSum_Pendulum.SD.sum.yaml'))
        np.testing.assert_almost_equal(f['CB_frequencies'].ravel(),[2.571561E-02,5.154897E+00,3.448768E+01,3.639185E+01,9.826435E+01], 5)

        # Test toDataFrame
        df=f.toDataFrame(sortDim=2)
        np.testing.assert_almost_equal(df['z_[m]'].values,[-6,-1,0])
        np.testing.assert_almost_equal(df['GuyanMode1x_[m]'].values[0],0.6)

        # Test toJSON
        dJSON=f.toJSON('_test.json')
        np.testing.assert_almost_equal(dJSON['Connectivity'], [[0,1],[1,2]])
        try:
            os.remove('_test.json')
        except:
            pass


    def test_FASTSumGraph(self):
        f = FASTSummaryFile(os.path.join(MyDir, 'FASTSum_Pendulum.SD.sum.yaml'))
        graph = f.toGraph()
       # print(graph)
        self.assertEqual(len(graph.Nodes), 3)
        self.assertEqual(len(graph.Elements), 2)
        self.assertEqual(len(graph.Modes), 11)
        np.testing.assert_almost_equal(graph.Modes[10]['freq'], 98.26435)


if __name__ == '__main__':
#     Test().test_000_debug()
    unittest.main()
