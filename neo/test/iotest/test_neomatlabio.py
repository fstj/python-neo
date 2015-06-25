# -*- coding: utf-8 -*-
"""
Tests of neo.io.neomatlabio
"""

# needed for python 3 compatibility
from __future__ import absolute_import, division

try:
    import unittest2 as unittest
except ImportError:
    import unittest

from neo.test.iotest.common_io_test import BaseTestIO
from neo.io.neomatlabio import NeoMatlabIO, HAVE_SCIPY
import numpy as np

from neo.test.iotest.generate_datasets import generate_one_simple_block


@unittest.skipUnless(HAVE_SCIPY, "requires scipy")
class TestNeoMatlabIO(BaseTestIO, unittest.TestCase):
    """
    
    """
    ioclass = NeoMatlabIO
    files_to_test = []
    files_to_download = []

    def setUp(self):
        """
        
        """
        BaseTestIO.setUp(self)
        self.supported_objects = NeoMatlabIO.supported_objects
        self.test_file = self.get_filename_path('test_mat_neo.mat')

    def test_write(self):
        """
        Test that writes then reads and check it gets the same data

        """
        test_block = generate_one_simple_block(
            supported_objects = self.supported_objects)

        # 1 - write test_block
        writer = NeoMatlabIO(filename = self.test_file)
        writer.write_block(test_block)

        # 2 - read test_block
        read_block = writer.read_block(cascade=True, lazy=False)

        # 3 - compare blocks
        # FIXME: assertEqual doesn't seem to handle numpy arrays comparison,
        # modify the line below to really compare the 2 blocks
        self.assertEqual(test_block, read_block)

    #=========================================================================
    # DEFINITION TEST
    #=========================================================================

    # FIXME: The tests commented below don't seem to use assertRaises as it 
    # should be used. The tests should be rewritten to correct them.

#    def test_create_struct_from_obj(self):
#        """
#        DEF TEST : create_struct_from_obj
#        
#        """
#        # raise error if entry is not a Neo object
#        w1_test = {'a':1827, 'b': 'test', 'j': 0.4}
#        w2_test = ['test','test','test']
#        w3_test = 'test'
#        self.assertRaises(NeoMatlabIO.create_struct_from_obj(w1_test))
#        self.assertRaises(NeoMatlabIO.create_struct_from_obj(w2_test))
#        self.assertRaises(NeoMatlabIO.create_struct_from_obj(w3_test))

#    def test_create_ob_from_struct(self):
#        """
#        DEF TEST : create_ob_from_struct
#
#        """
#        # raise error if entry is not a struct format
#        w1_test = {'a':1827, 'b': 'test', 'j': 0.4}
#        w2_test = ['test','test','test']
#        w3_test = 'test'
#        self.assertRaises(NeoMatlabIO.create_ob_from_struct(w1_test))
#        self.assertRaises(NeoMatlabIO.create_ob_from_struct(w2_test))
#        self.assertRaises(NeoMatlabIO.create_ob_from_struct(w3_test))


#    def test_create_dict_from_struct(self):
#        """
#        DEF TEST : create_dict_from_struct

#        """
#        # raise error if entry is not a mat struct
#        w1_test = {'a':1827, 'b': 'test', 'j': 0.4}
#        w2_test = ['test','test','test']
#        w3_test = 'test'
#        self.assertRaises(NeoMatlabIO.create_dict_from_struct(w1_test))
#        self.assertRaises(NeoMatlabIO.create_dict_from_struct(w2_test))
#        self.assertRaises(NeoMatlabIO.create_dict_from_struct(w3_test))

#    def create_struct_array_from_container(self):
#        """
#        DEF TEST : create_struct_from_container

#        """
#        # raise error if entry is not a container (list, dict or np.array)
#        w1_test = 'test'
#        w2_test = 18284
#        w3_test = 0.4
#        self.assertRaises(NeoMatlabIO.create_struct_array_from_container(
#            w1_test))
#        self.assertRaises(NeoMatlabIO.create_struct_array_from_container(
#            w2_test))
#        self.assertRaises(NeoMatlabIO.create_struct_array_from_container(
#            w3_test))

#        # verif out is a numpy array
#        w1_test = {'a':1827, 'b': 'test', 'j': 0.4}
#        w2_test = ['test','test','test']
#        w3_test = np.array(['test','test','test'], dtype='S')
#        self.assertIsInstance(NeoMatlabIO.create_struct_array_from_container(
#            w1_test), np.ndarray)
#        self.assertIsInstance(NeoMatlabIO.create_struct_array_from_container(
#            w2_test), np.ndarray)
#        self.assertIsInstance(NeoMatlabIO.create_struct_array_from_container(
#            w3_test), np.ndarray)

#    def test_key_verif(self):
#        """
#        DEF TEST : key_verif

#        """
#        # raise error if key is not str or int
#        w1_test = {'a':1827, 'b': 'test', 'j': 0.4}
#        w2_test = 0.4
#        w3_test = np.array(['test','test','test'], dtype='S')
#        self.assertRaises(NeoMatlabIO.key_verif(w1_test))
#        self.assertRaises(NeoMatlabIO.key_verif(w2_test))
#        self.assertRaises(NeoMatlabIO.key_verif(w3_test))

#        # raise error if key contain invalid character
#        w1_test = '*$#@02938hfjfkzj_-()'
#        w2_test = 'test2&'
#        w3_test = 'test3\00x0'
#        self.assertRaises(NeoMatlabIO.key_verif(w1_test))
#        self.assertRaises(NeoMatlabIO.key_verif(w2_test))
#        self.assertRaises(NeoMatlabIO.key_verif(w3_test))

#        # raise error if key contain str key is too long
#        w1_test = 'test'*63
#        self.assertRaises(NeoMatlabIO.key_verif(w1_test))


if __name__ == "__main__":
    unittest.main()

