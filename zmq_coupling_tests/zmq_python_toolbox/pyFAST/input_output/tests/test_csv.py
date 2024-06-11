import unittest
import os
import numpy as np
from pyFAST.input_output.tests.helpers_for_test import MyDir, reading_test
from pyFAST.input_output import CSVFile

class Test(unittest.TestCase):

    def test_001_read_all(self, DEBUG=True):
        reading_test('CSV*.*', CSVFile)

    def DF(self,FN):
        """ Reads a file and return a dataframe """ 
        return CSVFile(os.path.join(MyDir,FN)).toDataFrame()
 
    def test_CSV(self):
        self.assertEqual(self.DF('CSVAutoCommentChar.txt').shape,(11,6))
 
        DF=self.DF('CSVColInHeader.csv')
        self.assertTrue(all(DF.columns.values==['ColA','ColB','ColC']))
        self.assertEqual(DF.shape,(2,3))
 
        DF=self.DF('CSVColInHeader2.csv')
        self.assertTrue(all(DF.columns.values==['ColA','ColB','ColC']))
        self.assertEqual(DF.shape,(2,3))
 
        DF=self.DF('CSVColInHeader3.csv')
        self.assertTrue(all(DF.columns.values==['ColA','ColB','ColC']))
        self.assertEqual(DF.shape,(2,3))

        #DF=self.DF('CSVComma_UTF16.csv') # TODO encoding
        #self.assertEqual(DF.shape,(4,3))
 
        self.assertEqual(self.DF('CSVComma.csv').shape,(4,2))
        self.assertEqual(self.DF('CSVDateNaN.csv').shape,(11,2))
        self.assertEqual(self.DF('CSVNoHeader.csv').shape,(4,2))
        self.assertEqual(self.DF('CSVSemi.csv').shape,(3,2))
        self.assertEqual(self.DF('CSVSpace_ExtraCol.csv').shape,(5,4))
        self.assertEqual(self.DF('CSVTab.csv').shape,(5,2))
 
        DF = self.DF('CSVTwoLinesHeaders.txt')
        self.assertEqual(DF.columns.values[-1],'GenTq_(kN m)')
        self.assertEqual(DF.shape,(9,6))

    def test_CSV_string(self):
        DF=self.DF('CSVxIsString.csv')
        self.assertEqual(DF.shape,(7,2))
        self.assertEqual(DF.columns.values[0],'Label_[-]')


if __name__ == '__main__':
    #Test().test_CSV()
    unittest.main()
