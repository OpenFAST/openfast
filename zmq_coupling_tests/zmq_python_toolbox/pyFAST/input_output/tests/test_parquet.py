import unittest
import os
from pyFAST.input_output.tests.helpers_for_test import MyDir, reading_test 
from pyFAST.input_output.parquet_file import ParquetFile

class Test(unittest.TestCase):

    def test_001_read_all(self, DEBUG=True):
        reading_test('ParquetFile*.*', ParquetFile)

    def DF(self,FN):
        """ Reads a file and return a dataframe """ 
        return ParquetFile(os.path.join(MyDir,FN)).toDataFrame()
 
    def test_ParquetFile(self):
        df=self.DF('ParquetFile_test.parquet')
        self.assertListEqual(list(df.columns),["Column1","Column 2","Column Str"])
        self.assertEqual(df.shape,(3,3))
        self.assertEqual(df.loc[0,"Column Str"],'abc')
        self.assertEqual(df.loc[0, "Column1"], 1)
 


if __name__ == '__main__':
#     Test().test_000_debug()
    unittest.main()
