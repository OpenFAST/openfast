import unittest
import os
import numpy as np
from pyFAST.input_output.tests.helpers_for_test import MyDir, reading_test
from pyFAST.input_output.vtk_file import VTKFile



class Test(unittest.TestCase):

    def test_001_read_all(self, DEBUG=True):
        reading_test('VTK*.*', VTKFile)

    def test_VTKStruct(self):
        f = VTKFile(os.path.join(MyDir, 'VTKStructuredPointsPointData.vtk'))
        np.testing.assert_almost_equal(f.points,f.point_data['DisXY'])

        np.testing.assert_almost_equal(f.xp_grid,[0,20,40])
        np.testing.assert_almost_equal(f.point_data_grid['DisXY'][:,0,0,0],[0,20,40])


if __name__ == '__main__':
#     Test().test_000_debug()
    unittest.main()
