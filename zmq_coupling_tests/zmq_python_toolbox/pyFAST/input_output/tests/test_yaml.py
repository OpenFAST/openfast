import unittest
import numpy as np    
import os as os
from pyFAST.input_output.mini_yaml import *


class Test(unittest.TestCase):
    def test_scalar(self):
        # Scalar integer
        D = yaml_read(text='k:  10 # comment')
        self.assertDictEqual(D, {'k':10})
        self.assertIsInstance(D['k'], int)
        # Scalar float
        D = yaml_read(text='k:   10.# comment')
        self.assertDictEqual(D, {'k':10})
        self.assertIsInstance(D['k'], float)
        # Scalar string
        D = yaml_read(text='k:s# comment')
        self.assertDictEqual(D, {'k':'s'})
        self.assertIsInstance(D['k'], str)

    def test_lists(self):
        # list integer
        D = yaml_read(text='k:[ 10 ,] # comment')
        self.assertDictEqual(D, {'k':[10]})
        self.assertTrue(np.issubdtype(D['k'][0], np.integer))
        #self.assertIsInstance(D['k'][0], np.int32)
        # list float
        D = yaml_read(text='k: [ 10., ]# comment')
        self.assertDictEqual(D, {'k':[10]})
        self.assertIsInstance(D['k'][0], float)
        # list string
        D = yaml_read(text='k: [ s ,] #comment')
        self.assertDictEqual(D, {'k': ['s']})
        self.assertIsInstance(D['k'][0], str)
        # empty
        D = yaml_read(text='k:[ ] # comment')
        self.assertIsInstance(D['k'], np.ndarray)
        self.assertTrue(len(D['k'])==0)

        # list integers
        D = yaml_read(text='k: [ 10, 20, 30,] # comment')
        np.testing.assert_array_equal(D['k'], [10,20,30])
        self.assertTrue(np.issubdtype(D['k'].dtype, np.integer))
        # list float
        D = yaml_read(text='k: [ 10., 20.5, 30,] # comment')
        np.testing.assert_array_equal(D['k'], [10,20.5,30])
        self.assertTrue(np.issubdtype(D['k'].dtype, float))
        # list str
        D = yaml_read(text='k: [ a  ,  b   ,   c,] # comment')
        np.testing.assert_array_equal(D['k'], ['a','b','c'])
        self.assertTrue(np.issubdtype(D['k'].dtype, str))

        # list mixed, for now all as strings
        D = yaml_read(text='k: [ 1  ,  b   ,   0,] # comment')
        np.testing.assert_array_equal(D['k'], ['1','b','0'])
        self.assertTrue(np.issubdtype(D['k'].dtype, str))

    def test_arrays(self):
        # Empty array
        D = yaml_read(text="""k: # comment
        - [ ] """)
        self.assertIsInstance(D['k'], np.ndarray)
        self.assertEqual(D['k'].shape, (1,0)) # that's a choice...

        # Float array
        D = yaml_read(text="""k:# comment
          - [ 1, 3.5 ]# comment """)
        self.assertIsInstance(D['k'], np.ndarray)
        self.assertEqual(D['k'].shape, (1,2))
        self.assertTrue(np.issubdtype(D['k'].dtype, float))
        np.testing.assert_array_equal(D['k'], [[1,3.5]])

        # int array
        D = yaml_read(text="""k:# comment
          - [ 1 ,   ]# comment
          - [2]

          """)
        self.assertIsInstance(D['k'], np.ndarray)
        self.assertEqual(D['k'].shape, (2,1))
        self.assertTrue(np.issubdtype(D['k'].dtype, np.integer))
        np.testing.assert_array_equal(D['k'], [[1],[2]])

        # string array
        D = yaml_read(text="""k:# comment
          - [ a , b ,c   ]
          - [d, e  ,f  , ]   
        # comment
          """)
        self.assertIsInstance(D['k'], np.ndarray)
        self.assertEqual(D['k'].shape, (2,3))
        self.assertTrue(np.issubdtype(D['k'].dtype, str))
        np.testing.assert_array_equal(D['k'], [['a','b','c'],['d','e','f']])

        # Mixed array
        D = yaml_read(text="""k:# comment
          - [ a, 3.5 ]# comment """)
        self.assertIsInstance(D['k'], np.ndarray)
        self.assertEqual(D['k'].shape, (1,2))
        self.assertTrue(np.issubdtype(D['k'].dtype, str))
        np.testing.assert_array_equal(D['k'], [['a','3']])


if __name__=='__main__':
    #Test().test_arrays()
    unittest.main()
