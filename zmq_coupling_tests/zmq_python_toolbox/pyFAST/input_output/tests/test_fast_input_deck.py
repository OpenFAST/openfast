import unittest
import os
import numpy as np
from pyFAST.input_output.tests.helpers_for_test import MyDir, reading_test 
from pyFAST.input_output.fast_input_deck import FASTInputDeck

class Test(unittest.TestCase):

    def test_deck_driver(self):
        F=FASTInputDeck(os.path.join(MyDir,'input_decks/Main_EllipticalWingInf_OLAF.dvr'))
        #F.test_ascii(bCompareWritesOnly=True,bDelete=True)
        self.assertEqual(F.fst['NumTurbines'],1)

        self.assertEqual(F.version,'AD_driver')
        self.assertEqual(F.ADversion,'AD15')
        self.assertTrue(F.fst is not None)
        self.assertTrue(F.IW is not None)
        self.assertTrue(F.AD is not None)
        self.assertTrue(F.AD.Bld1 is not None)



if __name__ == '__main__':
    #Test().test_FASTIn()
    unittest.main()
