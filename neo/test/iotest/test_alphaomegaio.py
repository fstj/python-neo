# -*- coding: utf-8 -*-
"""
Tests of neo.io.alphaomegaio
"""

# needed for python 3 compatibility
from __future__ import absolute_import, division

try:
    import unittest2 as unittest
except ImportError:
    import unittest

from neo.io import AlphaOmegaIO
from neo.test.iotest.common_io_test import BaseTestIO
import scipy


class TestAlphaOmegaIO(BaseTestIO, unittest.TestCase):
    """
    
    """

    ioclass = AlphaOmegaIO

    files_to_test = ['File_AlphaOmega_1.map',
                     'File_AlphaOmega_2.map']
    files_to_download = files_to_test

    def setUp(self):
        """
        
        """
        BaseTestIO.setUp(self)
        self.AOreader = AlphaOmegaIO(
            filename=self.get_filename_path(self.files_to_test[0]))
        self._change_options()

    def _change_options(self, lazy = False, cascade = True,
                        read_level = 'AnalogSignal', read_waveform = False):
        """
        Change options and read the Neo block in the test file

        """
        self.AOreader.option_read_level = read_level
        self.AOreader.option_read_waveform = read_waveform
        self.test_block =  self.AOreader.read_block(lazy = lazy,
                                                    cascade = cascade)

    #=========================================================================
    # Objects test
    #=========================================================================

    def test_lenObject(self):
        """
        Test the number of neo object expected

        """
        # NOT LAZY

        self._change_options(lazy = False)
        # Test len(segments)
        segments = self.test_block.segments
        self.assertEqual(len(segments), 1)
        # Test len(AnalogSignal)
        analogsignals = segments[0].analogsignals
        self.assertEqual(len(analogsignals), 80)
        # Test len(SpikeTrain)
        spiketrains = segments[0].spiketrains
        self.assertEqual(len(spiketrains), 0)
        # Test len(EventArray)
        eventarrays = segments[0].eventarrays
        self.assertEqual(len(eventarrays), 52)

        # LAZY

        self._change_options(lazy = True)
        # Test len(segments)
        segments = self.test_block.segments
        self.assertEqual(len(segments), 1)
        # Test len(AnalogSignal)
        analogsignals = segments[0].analogsignals
        self.assertEqual(len(analogsignals), 80)
        # Test len(SpikeTrain)
        spiketrains = segments[0].spiketrains
        self.assertEqual(len(spiketrains), 0)
        # Test len(EventArray)
        eventarrays = segments[0].eventarrays
        self.assertEqual(len(eventarrays), 52)

    #=========================================================================
    # Options tests
    #=========================================================================
    
    def test_read_level(self):
        """
        OPTIONS TEST : option_read_level
        Test read_level options

        """
                    
        # 1 - 'AnalogSignal' option
        self._change_options(read_level = 'AnalogSignal')
        segments = self.test_block.segments
        analogsignals = segments[0].analogsignals
        self.assertEqual(len(analogsignals), 80)

        # 2 - 'SpikeTrain' option
        self._change_options(read_level = 'SpikeTrain')
        segments = self.test_block.segments
        spiketrains = segments[0].spiketrains
        self.assertEqual(len(spiketrains), 16)

        # 3 - Wrong option
        self.assertRaises(self._change_options, read_level = 'WRONG')

    def test_read_waveform(self):
        """
        OPTIONS TEST : option_read_waveform
        Test the presence or absence of spiktrains waveforms

        """

        # 1 - option is False, return must be None
        self._change_options(read_level = 'SpikeTrain', read_waveform = False)
        
        for spkt in self.test_block.segments[0].spiketrains:
            self.assertIsNone(spkt.waveforms)

        # 2 - option is True, return must a 2D array
        self._change_options(read_level = 'SpikeTrain', read_waveform = True)

        for spkt in self.test_block.segments[0].spiketrains:
            self.assertIsNotNone(spkt.waveforms)

        # 3 - Wrong option
        self.assertRaises(self._change_options, read_level = 'SpikeTrain',
            read_waveform = 'WRONG')

    #=========================================================================
    # OBLIGATORY ATTRIBUTES PRESENCE Tests
    #=========================================================================

    def test_groupsPresent(self):
        """
        Test the presence/absence of dict "groups" in block annotations

        """
        # not lazy
        self._change_options(lazy = False)
        self.assertIn('groups', self.test_block.annotations)

        # lazy
        self._change_options(lazy = True)
        self.assertIn('groups', self.test_block.annotations)

        # not cascade
        self._change_options(cascade = False)
        self.assertNotIn('groups', self.test_block.annotations)

    #=========================================================================
    # FILE SPECIFICITY Tests
    #=========================================================================

    def test_stopData(self):
        """
        Test specific file with lost datas

        """

    def test_verifAttr(self):
        """
        
        """
    
    def test_comparFile(self):
        """
        Test capability if AlphaOmegaIO to get same informations than the
        convert file in MAT format, from Mapfile.exe software
        
        """
        # FIXME: Only the metadata are compared below, the data in analogsignals
        # and eventarrays should also be compared.
                
        # 1- Reading of MAP file
        self._change_options()
        
        # 2- Reading of MAT file
        # FIXME: The mat-file File_AlphaOmega_1.mat is not yet on the gnode
        # it should be added
        mat_file = self.get_filename_path('File_AlphaOmega_1.mat')
        mat_struct = scipy.io.loadmat(mat_file, 
                                  struct_as_record=False, 
                                  squeeze_me=True)
        
        # 3- Comparison
        segments = self.test_block.segments
        # AnalogSignal
        analogsignals = segments[0].analogsignals
        for anasig in analogsignals:
            if anasig.name in mat_struct:
                self.assertEqual(len(anasig), len(mat_struct(anasig.name)))
                self.assertAlmostEqual(anasig.sampling_rate, 
                                 mat_struct(anasig.name + "_KHz"))
                self.assertAlmostEqual(anasig.t_start, 
                                 mat_struct(anasig.name + "_TimeBegin"))
        # Eventarrays
        eventarrays = segments[0].eventarrays
        for ea in eventarrays:
            if ea.name in mat_struct:
                self.assertEqual(len(ea), len(mat_struct(ea.name)))
                self.assertAlmostEqual(ea.times, len(mat_struct(ea.name)))
                self.assertAlmostEqual(anasig.annotations('sampling_rate'), 
                    mat_struct(anasig.name + "_KHz"))
                self.assertAlmostEqual(anasig.annotations('t_start'), 
                    mat_struct(anasig.name + "_TimeBegin"))
        
if __name__ == "__main__":
    unittest.main()

