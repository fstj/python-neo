# -*- coding: utf-8 -*-
"""

Class for reading data from Alpha Omega .map files.

==============================================================================
This code is written from the incomplete file specifications available in:

[1] AlphaMap Data Acquisition System User's Manual Version 10.1.1
Section 5 APPENDIX B: ALPHAMAP FILE STRUCTURE, pages 120-140
Edited by ALPHA OMEGA Home Office: P.O. Box 810, Nazareth Illit 17105, Israel
http://www.alphaomega-eng.com/

[2] AlphaMap User's Manual Version 1.2.1
Section 6: ALPHALAB SNR FILE FORMAT, pages 82-89
Edited by ALPHA OMEGA

and from the source code of a C software for conversion of .map files to
.eeg elan software files :

[3] alphamap2eeg 1.0, 12/03/03, Anne CHEYLUS - CNRS ISC UMR 5015

==============================================================================
Supported : Read

@author : sgarcia, Lea Andre, Florent Jaillet

"""

# NOTE: For some specific types of comments, the following convention is used:
# - "TODO:" Desirable future evolution
# - "WARNING:" Information about code that is based on broken or missing
#   specifications and that might be wrong
# - "FIXME:" Something that is not working as it should or not implemented as it 
#   should


# Main limitations of this reader:
# - This reader hasn't been thoroughly tested and is therefore probably buggy in
#   some places.
# - Loading multichannel signals as AnalogSignalArrays is not supported (they 
#   are loaded as seperated AnalogSignals).
# - This reader only reads the whole data in the file (returned in a Neo Block),
#   without providing any way to load only part of the data.

# needed for python 3 compatibility
from __future__ import absolute_import, division

# specific imports
import datetime
import logging
import os
import struct
from sets import Set

# file no longer exists in Python3
try:
    file
except NameError:
    import io
    file = io.BufferedReader

# note neo.core needs only numpy and quantities
import numpy as np
import quantities as pq

from neo.io.baseio import BaseIO
from neo.core import Block, Segment, AnalogSignal, EventArray, SpikeTrain
from neo.io.tools import create_many_to_one_relationship
from neo.io.tools import populate_RecordingChannel


class AlphaOmegaIO(BaseIO):
    """
    Class for reading data from Alpha Omega .map files

    Usage:
        >>> from neo import io
        >>> r = io.AlphaOmegaIO(filename='File_AlphaOmega_1.map')
        >>> r.option_read_level('SpikeTrain')
        >>> r.option_read_waveform = True
        >>> blck = r.read_block(lazy=False, cascade=True)
        >>> print blck.segments[0].analogsignals

    :var str option_read_level: name of neo type used to store level data
        2 options :
            - 'AnalogSignal' :
               DEFAULT. level signal are loaded in AnalogSignal neo object
            - 'SpikeTrain' :
               level signal are loaded in SpikeTrain neo object
    :var boolean option_read_waveform: if True, waveforms are read for
        SpikeTrain object
    """

    # FIXME: When using cascade=False, the returned block doesn't contain all 
    # the annotations that are obtained when using cascade=True, because some of
    # them are gathered while reading the signals in the file.
    # For example informations on channels groups are missing (see _print_group
    # for information about the corresponding data).
    # The code should be restructured to get all the annotations even when
    # cascade=False.
    
    is_readable = True  # This is a reading only class
    is_writable = False  # writting is not supported

    # This class is able to directly or inderectly read the following kind of
    # objects
    supported_objects = [Block, Segment, AnalogSignal, EventArray, SpikeTrain]
    # TODO: Add support for other objects that should be extractable from .map
    # files (AnalogSignalArray, Event, Epoch?, Epoch Array?, Spike?)

    readable_objects = [Block]  # This class can only return a Block
    # TODO: create readers for different type of objects (Segment,
    # AnalogSignal,...)

    writeable_objects = []  # This class is not able to write objects

    # This is for GUI stuff : a definition for parameters when reading.
    read_params = {Block: []}

    write_params = None  # Writing is not supported, so no GUI stuff

    name = 'AlphaOmega'
    extensions = ['map', 'mpx']
    mode = 'file'

    # write is not supported so I do not overload write method from BaseIO

    def __init__(self, filename=None, option_read_level="SpikeTrain",
                 option_read_waveform=False):
        """
        :param str filename: The Alpha Omega file name. DEFAULT = None.
        :param str option_read_level:  name of neo type used to store level
            data. Options:
                - 'SpikeTrain' : DEFAULT
                - 'AnalogSignal'
        :param boolean option_read_waveform: if True, waveforms are read for
            SpikeTrain object

        """
        BaseIO.__init__(self)
        self.filename = filename
        self.option_read_level = option_read_level
        self.option_read_waveform = option_read_waveform

        # Logger definition
        # NOTE: Since 10 April, self.logger is already present with
        # BaseIO.__init__ method, but with actual version, it isn't present
        # TODO: remove this logger part if present in BaseIO class

        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.WARNING)
        ch = logging.StreamHandler()
        ch.setLevel(logging.WARNING)
        formatter = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
        ch.setFormatter(formatter)
        self.logger.addHandler(ch)

    #=========================================================================
    #   Public options
    #=========================================================================

    def option_read_level(self, option_read_level):
        """
        option_read_level setter.

        :param str option_read_level: name of neo object to store level data
            options :
                - 'SpikeTrain' : DEFAULT
                - 'AnalogSignal'

        """
        options = ['AnalogSignal', 'SpikeTrain']

        try:
            option_read_level in options
        except:
            raise Exception("Option %s is not available." % (option_read_level))
        else:
            self.option_read_level = option_read_level

    #=========================================================================
    #   Private read data method
    #=========================================================================

    def _read_analog_data(self, fid, list_blocks, file_blocks, sampling_rate,
                          amplitude):
        """
        Read data present in block 5, and return data and information

        Data type read:
            - Continuous Data
            - Level Data

        :param file fid: file object, ready to read
        :param list list_blocks: list of indexes of data blocks
        :param list file_blocks: list of blocks
        :param float sampling_rate: sample rate for the read channel
        :param float amplitude: signal amplitude

        :return np.array signal_array: array fill with the analog signal
        :return int ind_start: start index of analog signal

        """
        # NOTE: we assume that the data block are chronologically ordered, so 
        # the first block contains the ind_start and with the last block, we 
        # can find the ind_stop of the signal.

        # Start index
        count = self._count_samples(file_blocks[list_blocks[0]]['m_length'])
        fid.seek(file_blocks[list_blocks[0]]['pos'] + 6 + (count * 2))
        ind_start = struct.unpack('I', fid.read(4))[0]

        # End index
        count = self._count_samples(
            file_blocks[list_blocks[len(list_blocks) - 1]]['m_length'])
        fid.seek(file_blocks[list_blocks[len(list_blocks) - 1]]['pos'] +
                 6 + (count * 2))
        ind_stop = struct.unpack('I', fid.read(4))[0] + count

        # Initialize signal array
        chan_signal_len = ind_stop - ind_start
        signal_array = np.empty(chan_signal_len)
        signal_array.fill(np.nan)

        for ind_block in list_blocks:
            # Count number of samples
            count = self._count_samples(file_blocks[ind_block]['m_length'])

            # Find start of block
            fid.seek(file_blocks[ind_block]['pos'] + 6 + (count * 2))
            start = struct.unpack('I', fid.read(4))[0]

            # Data reading
            start_block = start - ind_start
            end_block = start - ind_start + count
            fid.seek(file_blocks[ind_block]['pos'] + 6)
            signal_array[start_block:end_block] = np.fromfile(
                fid,
                dtype=np.int16,
                count=count)

        # Treatment for ind_start and signal_array
        t_start = (ind_start / sampling_rate)
        quantum = float(amplitude / (2 ** 15))
        signal_array *= quantum

        return signal_array, t_start

    def _read_digital_data(self, fid, list_blocks, file_blocks, chan_len,
                          ind_chan, ind, sampling_rate):
        """
        Read data present in block 5, and return 2 temp_array

        Data type read:
            - Digital Data

        :param file fid: file object, ready to read
        :param list list_blocks: list of indexes of data blocks
        :param list file_blocks: list of blocks
        :param np.array chan_len: array of len(total data) for a channel
        :param int ind_chan: index of channel in 'chan_len' vector
        :param int ind: index in the data vector
        :param float sampling_rate: sample rate for the read channel

        :return np.array temp_array_up: array fill with times for UP event
        :return np.array temp_array_down: array fill with times for DOWN event

        """
        # Initialize array
        temp_array = np.empty(chan_len[ind_chan], dtype=np.float64)
        labels = np.empty(chan_len[ind_chan], dtype=np.int16)

        for ind_block in list_blocks:
            # Count number of samples
            count = self._count_samples(file_blocks[ind_block]['m_length'])

            # Position in file for times data
            fid.seek(file_blocks[ind_block]['pos'] + 6)
            temp_array[ind:ind + count] = np.fromfile(
                fid,
                dtype=np.int32,
                count=count)

            # Position in file for labels data
            fid.seek(file_blocks[ind_block]['pos'] + 6 + 4)
            labels[ind:ind + count] = np.fromfile(
                fid,
                dtype=np.int16,
                count=count)
            ind += count

        temp_array *= pq.ms

        # Separate in Up and Down temp_array
        temp_array_up = temp_array[labels == 1]
        temp_array_down = temp_array[labels == 0]
        
        # Treatment
        times_up = temp_array_up / sampling_rate
        times_down = temp_array_down / sampling_rate

        return times_up, times_down

    def _read_spike_data(self, fid, list_blocks, file_blocks, sampling_rate,
                         amplitude):
        """
        Read data present in block 5, and return objects for SpikeTrain
        creation

        Data type read:
            - Level Data

        :param file fid: file object, ready to read
        :param list list_blocks: list of indexes of data blocks
        :param list : list of blocks
        :param float smpling_rate : sample rate of signal 
        :param float amplitude: signal amplitude

        :return np.array times: array fill with time (float)
        :return int t_start: start time of data
        :return int t_stop: stop time of data
        :return Quantity 3D waveforms:  object containing spike waveform

        """
        waveforms = np.array([])
        quantum = float(amplitude / (2 ** 15))

        # Start time
        count = self._count_samples(file_blocks[list_blocks[0]]['m_length'])
        fid.seek(file_blocks[list_blocks[0]]['pos'] + 6 + (count * 2))
        ind_start = struct.unpack('I', fid.read(4))[0]

        # End time
        count = self._count_samples(
            file_blocks[list_blocks[len(list_blocks) - 1]]['m_length'])
        fid.seek(file_blocks[list_blocks[len(list_blocks) - 1]]['pos'] +
                 6 + (count * 2))
        ind_stop = struct.unpack('I', fid.read(4))[0] + count

        # Initialize times and waveforms array
        count = self._count_samples(file_blocks[list_blocks[0]]['m_length'])
        times = np.empty(len(list_blocks))
        times.fill(np.nan)
        if self.option_read_waveform:
            waveforms = np.empty((len(list_blocks), 1, count))
            waveforms.fill(np.nan)

        for ind, ind_block in enumerate(list_blocks):
            # Count number of samples
            count = self._count_samples(file_blocks[ind_block]['m_length'])

            # Times data reading
            fid.seek(file_blocks[ind_block]['pos'] + 6 + (count * 2))
            time = struct.unpack('I', fid.read(4))[0]
            times[ind] = time

            # Waveforms data reading
            if self.option_read_waveform:
                # Waveform reading
                fid.seek(file_blocks[ind_block]['pos'] + 6)
                waveforms[ind, 0, :] = np.fromfile(
                    fid, dtype=np.int16, count=count)
                waveforms[ind, 0, :] *= quantum

        # Treatment 
        times /= sampling_rate
        t_start = (ind_start / sampling_rate)
        t_stop = (ind_stop / sampling_rate)
        
        return times, t_start, t_stop, waveforms

    def _read_port_data(self, fid, list_blocks, file_blocks, chan_len,
                       ind_chan, ind, sampling_rate):
        """
        Read data present in block 5, and return 2 temp_array

        Data type read:
            - Port Data

        :param file fid: file object, ready to read
        :param list list_blocks: list of indexes of data blocks
        :param list file_blocks: list of blocks
        :param np.array chan_len: array of len(total data) for a channel
        :param int ind_chan: index of channel in 'chan_len' vector
        :param int ind: index in the data vector

        :return np.array temp_array: array fill with time (float)
        :return np.array labels: array fill with labels (int)

        """
        # Initialize time and label array
        # NOTE: array are normally created with np.empty method, but it cause
        # a problem for neomatlab reading. It actually work with np.zeros
        # method.
        temp_array = np.zeros(chan_len[ind_chan], dtype=np.float64)
        labels = np.zeros(chan_len[ind_chan], dtype=np.int16)

        for ind_block in list_blocks:
            # Count number of samples
            count = self._count_samples(file_blocks[ind_block]['m_length'])

            # Position in file for events or labels data
            fid.seek(file_blocks[ind_block]['pos'] + 6)
            labels[ind:ind + count] = np.fromfile(
                fid,
                dtype=np.int16,
                count=count)

            # Position in file for times data
            fid.seek(file_blocks[ind_block]['pos'] + 6 + 2)
            temp_array[ind:ind + count] = np.fromfile(
                fid,
                dtype=np.int32,
                count=count)
            ind += count

        # Treatment
        times = temp_array / sampling_rate

        return times, labels

#==============================================================================
#   Private annotate block method
#==============================================================================

    def _annotate_block(self, neo_object, block):
        """
        Annotate neo_object with data present in information block

        Block supported:
            - 2
            - b

        :param Neo Object neo_object: EventArray, AnalogSignal, ...
        :param dict block: dictionnary which contains block information

        """
        block_available = ['2', 'b']

        try:
            block['m_TypeBlock'] in block_available
        except:
            raise Exception("Annotate is not available for block type %s"
                % (block['m_TypeBlock']))
        else:
            if block['m_TypeBlock'] == '2':
                self._annotate_block2(neo_object, block)
            elif block['m_TypeBlock'] == 'b':
                self._annotate_blockb(neo_object, block)

    def _annotate_block2(self, neo_object, block):
        """
        Annotate neo_object with data present in block 2

        :param Neo Object neo_object: EventArray, AnalogSignal, ...
        :param dict block: dictionnary which contains block information

        """
        neo_object.annotate(isInput=block['m_isInput'])
        neo_object.annotate(numColor=block['m_numColor'])
        neo_object.annotate(channel_type=block['type_subblock'])

        if block['type_subblock'] != 'digital':
            neo_object.annotate(Amplitude=block['m_Amplitude'] * pq.V)

        if block['type_subblock'] == 'digital':
            neo_object.annotate(SaveTrigger=block['m_SaveTrigger'])
            neo_object.annotate(Duration=block['m_Duration'] * pq.ms)
            neo_object.annotate(PreviousStatus=block['m_PreviousStatus'])

        if block['type_subblock'] == 'level':
            neo_object.annotate(LevelValue=block['m_LevelValue'] * pq.V)
            neo_object.annotate(TrgMode=block['m_TrgMode'])
            neo_object.annotate(YesRMS=block['m_YesRMS'])
            neo_object.annotate(bAutoScale=block['m_bAutoScale'])

        if (block['type_subblock'] == 'level' or
            block['type_subblock'] == 'external_trigger'):
                neo_object.annotate(nSpikeCount=block['m_nSpikeCount'])
                neo_object.annotate(
                    nPreTrigmSec=block['m_nPreTrigmSec'] * pq.ms)
                neo_object.annotate(
                    nPostTrigmSec=block['m_nPostTrigmSec'] * pq.ms)

        if (block['type_subblock'] == 'external_trigger' or
            block['type_subblock'] == 'digital'):
                neo_object.annotate(channel_index=int(block['m_numChannel']))

    def _annotate_blockb(self, neo_object, block):
        """
        Annotate neo_object with data present in block b

        :param Neo Object neo_object: EventArray, AnalogSignal, ...
        :param dict block: dictionnary which contains block information

        """
        neo_object.annotate(BoardNumber=block['m_BoardNumber'])
        neo_object.annotate(channel_index=block['m_numChannel'])
        neo_object.annotate(channel_type='port')

#==============================================================================
#   Private specific method
#==============================================================================

    def _count_samples(self, m_length):
        """
        Count the number of signal samples available in a type 5 data block

        :param int m_length: length of block 5
        :return int count: number of samples in block 5

        """
        # INFO: for information about type 5 data block, see [1].
        # -6 corresponds to the header of block 5, and the -2 take into
        # account the fact that last 2 values are not available as the 4
        # corresponding bytes are coding the time stamp of the beginning
        # of the block
        count = int(((m_length - 6) / 2) - 2)
        return count

    def _print_group(self, groups):
        """
        Pretty print groups data

        :param list groups: list with groups and subgroups informations

        """
        # INFO: groups information structure:
        #       groups = [group, group, group, ...]
        #       group = {
        #           'name': group_name,
        #           'number': group_number,
        #           'Z_Order': appearence_order,
        #           'subgroups': [subgroup, subgroups, subgroup, ...]}
        #       subgroup = {
        #           'name': group_name,
        #           'number': group_number,
        #           'Z_Order': appearence_order,
        #           'TypeOverLap': typeoverlap,
        #           'channels': [channel, channel, channel, ...]}
        #       general and groups information are stored in
        #       annotations of Block.

        print("")
        for group in groups:
            print("GROUP: ", group['name'])
            print("NB SUBGROUP: ",len(group['subgroups']))
            for subgroups in group['subgroups']:
                print("---SUBGROUP: ", subgroups['name'])
                print("---NB CHANNELS: ", len(subgroups['channels']))
                print("---CHANNELS: ", subgroups['channels'])
        print("")

    def _select_available_channel(self, groups, list_available_channel):
        """
        Sort on groups and subgroups with available channel.

        :param list groups: list of group (dict)
        :param list list_available_channel: list of available channel

        :return list groups: list of available group (dict)

        Groups information structure
        -----------------------------
        ""groups = [group, group, group, ...]
        group = {
            'name': group_name,
            'number': group_number,
            'Z_Order': appearence_order,
            'subgroups': [subgroup, subgroups, subgroup, ...]}
        subgroup = {
            'name': group_name,
            'number': group_number,
            'Z_Order': appearence_order,
            'TypeOverLap': typeoverlap,
            'channels': [channel, channel, channel, ...]}""

        """

        # Step 1: delete unavailable channels
        for group in groups:
            for subgroup in group['subgroups']:
                new_channels = []
                for channel in subgroup['channels']:
                    if channel in list_available_channel:
                        new_channels.append(channel)
                subgroup['channels'] = new_channels

        # Step 2: delete empty subgroups (without available channels)
        for group in groups:
            new_subgroups = []
            for subgroup in group['subgroups']:
                if len(subgroup['channels']) != 0:
                    new_subgroups.append(subgroup)
            group['subgroups'] = new_subgroups

        # Step 3: delete empty groups (without available subgroups)
        new_groups = []
        for group in groups:
            if len(group['subgroups']) != 0:
                new_groups.append(group)

        groups = new_groups

        return groups

#==============================================================================
#   Public read method
#==============================================================================

    def read_block(self, lazy=False, cascade=True):
        """
        Return a Block.

        :param boolean lazy: if True, the numerical data is not loaded, only
            properties and metadata
        :param boolean cascade: if False, only a single object is loaded,
            without its references to other objects

        :return Block blck: Neo Block

        """

        # create the neo Block that will be returned at the end
        blck = Block(file_origin=os.path.basename(self.filename))

        fid = open(self.filename, 'rb')

        #=====================================================================
        # Step 1: read the headers of all the data blocks to load the file
        #         structure
        #=====================================================================
        # NOTE: in the following, the word "block" is used in the sense
        # used in the alpha-omega specifications (ie a data chunk in the
        # file), rather than in the sense of the usual Block object in neo

        pos_block = 0  # position of the current block in the file
        file_blocks = []  # list of data blocks available in the file
        list_type_block = []  # list of data blocks type in the file

        information = {}  # dict containing general information
        groups = []  # list containing groups information
        out_info = [  # info not returned in information dict
            "pos",
            "blank",
            "blank2",
            "m_length",
            "m_TypeBlock",
            "m_nextBlock",
            "m_EraseCount",
            "m_Reserved",
            "m_placeMainWindow"]
        dap_info = [  # info for log
            "DAP_buffers_filling",
            "RMS_value",
            "num_channel",
            "sample_count"]
        i_group = -1 # index of group
        i_subgroup = 0  # index of subgroup

        if not cascade:  # we read only the main header

            m_length, m_TypeBlock = struct.unpack('Hcx', fid.read(4))
            # m_TypeBlock should be 'h', as we read the first block
            block = HeaderReader(
                fid,
                dict_header_type.get(m_TypeBlock.decode("utf-8"),
                                     Type_Unknown)).read_f()
            block.update({
                'm_length': m_length,
                'm_TypeBlock': m_TypeBlock,
                'pos': pos_block})
            file_blocks.append(block)
            list_type_block.append(m_TypeBlock)

            # retrieve date and time information
            if m_TypeBlock == 'h':
                blck.rec_datetime = datetime.datetime(
                    block['m_date_year'],
                    block['m_date_month'],
                    block['m_date_day'],
                    block['m_time_hour'],
                    block['m_time_minute'],
                    block['m_time_second'],
                    10000 * block['m_time_hsecond'])
                    # the 10000 is here to convert m_time_hsecond from
                    # centisecond to microsecond

        else:  # cascade == True

            seg = Segment(file_origin=os.path.basename(self.filename))
            blck.segments.append(seg)

            while True:
                first_4_bytes = fid.read(4)

                if len(first_4_bytes) < 4:
                    # we have reached the end of the file
                    break
                else:
                    m_length, m_TypeBlock = struct.unpack('Hcx', first_4_bytes)
                
                # for Python 3 compatibility we need to convert m_TypeBlock to 
                # utf-8
                try:
                    m_TypeBlock = m_TypeBlock.decode("utf-8")
                except:
                    continue

                block = HeaderReader(
                    fid,
                    dict_header_type.get(m_TypeBlock,
                                         Type_Unknown)).read_f()

                block.update({
                    'm_length': m_length,
                    'm_TypeBlock': m_TypeBlock,
                    'pos': pos_block})
                list_type_block.append(m_TypeBlock)

                #========================================
                # a - Read of subblock for block 2 and 7
                #========================================

                if m_TypeBlock == '2':
                    # The beginning of the block of type '2' is identical for
                    # all types of channels, but the following part depends on
                    # the type of channel. So we need a special case here.

                    # WARNING: How to check the type of channel is not
                    # described in the documentation. But the following tests 
                    # seem to work (adapted from C code [3])

                    type_subblock = 'unknown_channel_type'
                    description = Type2_SubBlockUnknownChannels
                    block.update({'m_Name': 'unknown_name'})

                    if block['m_isAnalog'] == 0:
                        # digital channel
                        type_subblock = 'digital'
                        description = Type2_SubBlockDigitalChannels
                    elif block['m_isAnalog'] == 1:
                        # analog channel
                        
                        # NOTE: analog channel block have a different
                        # header than digital channel.
                        
                        # Here is the reading of the common part for analog 
                        # channels
                        description = Type2_SubBlockAnalogChannels
                        subblock = HeaderReader(fid, description).read_f()
                        block.update(subblock)

                        if block['m_Mode'] == 1:
                            # level channel
                            type_subblock = 'level'
                            description = Type2_SubBlockLevelChannels
                        elif block['m_Mode'] == 2:
                            # external trigger channel
                            type_subblock = 'external_trigger'
                            description = Type2_SubBlockExtTriggerChannels
                        else:
                            # continuous channel
                            type_subblock = 'continuous (Mode %i)' % (
                                block['m_Mode'])
                            description = Type2_SubBlockContinuousChannels

                    subblock = HeaderReader(fid, description).read_f()
                    block.update(subblock)
                    block.update({'type_subblock': type_subblock})

                if m_TypeBlock == '7':
                    # The beggining of the block '7' is identical, data
                    # next are different according to the documentation [1].

                    description = Type7_DataSubblockUnknown

                    if block['FINT'] == -111:
                        description = Type7_DataSubblockDAPBuffers
                    elif block['FINT'] == -222:
                        description = Type7_DataSubblockRMS
                    elif block['FINT'] == -444:
                        description = Type7_DataSubblockRestart
                    elif block['FINT'] == -333:
                        description = Type7_DataSubblockDataLoss
                    elif block['FINT'] == -557 or block['FINT'] == -558:
                        description = Type7_DataSubblockStopAnalogOutput

                    subblock = HeaderReader(fid, description).read_f()
                    block.update(subblock)

                #=============================================================
                #  b - Annotate the Neo Block and log information
                #=============================================================
                #  INFO: groups information structure:
                #       groups = [group, group, group, ...]
                #       group = {
                #           'name': group_name,
                #           'number': group_number,
                #           'Z_Order': appearence_order,
                #           'subgroups': [subgroup, subgroup, subgroup, ...]}
                #       subgroup = {
                #           'name': subgroup_name,
                #           'number': subgroup_number,
                #           'Z_Order': appearence_order,
                #           'TypeOverLap': typeoverlap,
                #           'channels': [channel, channel, channel, ...]}
                #       general and groups information is stored in
                #       annotations of Block.

                # retrieve date and time information
                elif m_TypeBlock == 'h':
                    blck.rec_datetime = datetime.datetime(
                        block['m_date_year'],
                        block['m_date_month'],
                        block['m_date_day'],
                        block['m_time_hour'],
                        block['m_time_minute'],
                        block['m_time_second'],
                        10000 * block['m_time_hsecond'])
                        # the 10000 is here to convert m_time_hsecond from
                        # centisecond to microsecond

                    seg.rec_datetime = blck.rec_datetime.replace()
                    # I couldn't find a simple copy function for datetime,
                    # using replace without arguments is a twisted way to make
                    # a copy

                # retrieve general information
                if m_TypeBlock in ['h', '0']:
                    for info in block:
                        if (info in out_info or
                            info.startswith("m_MainWindow_") or
                            info.startswith("m_date_") or
                            info.startswith("m_time_")):
                                pass
                        else:
                            if info in info_quantities:
                                information[info[2:]] = (
                                    block[info] * info_quantities[info])
                            else:
                                information[info[2:]] = block[info]

                if block['m_TypeBlock'] == '1':
                    for info in block:
                        if info not in out_info:
                            if info in info_quantities:
                                information["Board_%i_%s" % (
                                 block["m_Number"],
                                 info[2:])] = (
                                     block[info] * info_quantities[info])
                            else:
                                information["Board_%i_%s" % (
                                 block["m_Number"],
                                 info[2:])] = block[info]

                # retrieve groups information
                if m_TypeBlock == '3':
                    i_group += 1
                    i_subgroup = 0
                    dict_group = {
                        'name': block['m_nameGroup'],
                        'number': block['m_Number'],
                        'Z_Order': block['m_Z_Order'],
                        'subgroups': []}
                    groups.append(dict_group)
                    

                if m_TypeBlock == '4':
                    dict_subgroup = {
                        'name': block['m_Name'],
                        'number': block['m_Number'],
                        'Z_Order': block['m_Z_Order'],
                        'TypeOverlap': block['m_TypeOverlap'],
                        'channels': []}
                    groups[i_group]['subgroups'].append(dict_subgroup)
                    for info in block:
                        if (info.startswith('m_numChannel') and
                            block[info] != 0):
                                groups[i_group]['subgroups'] \
                                      [i_subgroup]['channels'] \
                                      .append(block[info])
                    i_subgroup += 1

                # display warnings about DAP information
                if m_TypeBlock == '7':
                    error = "\nBoard %i - Error %i" % (
                        block['m_numBoard'],
                        block['FINT'])
                    for info in block:
                        if info in dap_info:
                            error += " - %s = %s" % (info, block[info])
                    if block['FINT'] == -111:
                        error = "DAP Buffers - " + error
                    elif block['FINT'] == -222:
                        error = "RMS - " + error
                    elif block['FINT'] == -444:
                        error = "DAP Restart - " + error
                    elif block['FINT'] == -333:
                        error = "Data Loss " + error
                        error += " - Loss between %i - %i" % (
                            block["first_lost_sample"],
                            block["last_lost_sample"])
                    elif block['FINT'] == -666:
                        error = "Start Analog Output - " + error
                    elif block['FINT'] == -557 or block['FINT'] == -558:
                        error = "Stop Analog Output - " + error

                    self.logger.info(error)

                file_blocks.append(block)
                pos_block += m_length
                fid.seek(pos_block)

            #==============================================================
            # Step 2: find the available channels
            #==============================================================
            #         NOTE: Block 2 contains information for continuous,
            #         digital, level and trigger channels. Block b contains
            #         information for Port channels, except the sampling_rate
            #         present in corresponding digital channels (number 11033
            #         to 11040)

            list_chan = []  # list containing indexes of channel blocks
            available_port_channel = list(range(11033, 11041))
            port_sampling_rate = None

            for ind_block, block in enumerate(file_blocks):
                if block['m_TypeBlock'] == '2' or block['m_TypeBlock'] == 'b':
                    list_chan.append(ind_block)
                    if block['m_numChannel'] in available_port_channel:
                        if port_sampling_rate is None:
                            port_sampling_rate = block['m_SampleRate']    
                        if block['m_SampleRate'] != port_sampling_rate:
                            self.logger.warning(
                                'Sample_Rate Error : %s is different of %s' %
                                (block['m_SampleRate'], port_sampling_rate))

            #================================================================
            # Step 3: find blocks containing data for the available channels
            #================================================================

            list_data = []  # list of lists of indexes of data blocks
                            # corresponding to each channel
            list_available_channel = Set()  # list of available numChannel

            for ind_chan, chan in enumerate(list_chan):
                list_data.append([])
                num_chan = file_blocks[chan]['m_numChannel']
                for ind_block, block in enumerate(file_blocks):
                    if block['m_TypeBlock'] == '5':
                        if block['m_numChannel'] == num_chan:
                            list_available_channel.add(block['m_numChannel'])
                            list_data[ind_chan].append(ind_block)

            #================================================================
            # Step 4: compute the length (number of samples) of the channels
            #================================================================

            # NOTE: we compute the length only to know if there is data in
            # block 5 for continous ans level signal. The length to create
            # np.array is calculate later.
            chan_len = np.zeros(len(list_data), dtype=np.int)
            ind_valid_chan = Set()

            for ind_chan, list_blocks in enumerate(list_data):
                for ind_block in list_blocks:
                    chan_len[ind_chan] += self._count_samples(
                        file_blocks[ind_block]['m_length'])
                    if self._count_samples(
                        file_blocks[ind_block]['m_length']) != 0:
                            ind_valid_chan.add(ind_chan)

            #=================================================================
            # Step 6: load the data
            #     TODO give the possibility to load data as AnalogSignalArrays
            #=================================================================

            for ind_chan in ind_valid_chan:
                ind = 0  # index in the data vector
                list_blocks = list_data[ind_chan]

                # find the signal type
                if file_blocks[list_chan[ind_chan]]['m_TypeBlock'] == '2':
                    type_signal = file_blocks[list_chan[ind_chan]]['type_subblock']
                elif file_blocks[list_chan[ind_chan]]['m_TypeBlock'] == 'b':
                    type_signal = 'port'
                sampling_rate = \
                    file_blocks[list_chan[ind_chan]]['m_SampleRate']

                if ('continuous' in type_signal or type_signal == 'level'):
                    # read time stamp for the beginning of the signal
                    form = '<l'  # reading format
                    ind_block = list_blocks[0]

                    count = self._count_samples(
                        file_blocks[ind_block]['m_length'])
                    fid.seek(file_blocks[ind_block]['pos'] + 6 + (count * 2))
                    buf = fid.read(struct.calcsize(form))
                    val = struct.unpack(form, buf)
                    start_index = val[0]
                    
                    if sampling_rate != 0.:
                        t_start = start_index / sampling_rate
                    else:
                        t_start = start_index

                #=============================================================
                # CONTINUOUS signal and LEVEL signal with ANALOGSIGNAL option
                #=============================================================
                if ('continuous' in type_signal or
                    (type_signal == 'level' and
                     self.option_read_level == 'AnalogSignal')):

                    if lazy:
                        # FIXME: Check if we should add the quantity (so for 
                        # example V for signals) when doing lazy reading.
                        # If yes, add the quantity here.
                        signal_array = np.array([])

                    else:
                        amplitude = file_blocks[
                            list_chan[ind_chan]]['m_Amplitude']
                        (signal_array, t_start) = self._read_analog_data(
                            fid,
                            list_blocks,
                            file_blocks,
                            sampling_rate,
                            amplitude)

                    # NOTE: The reason why we have'nt stock immediately
                    # signal values in the AnalogSignal object is because
                    # it is longer to put values one after the other than
                    # put the total signal_array in.
                    ana_sig = AnalogSignal(
                        signal_array * pq.V,
                        sampling_rate=sampling_rate * pq.kHz,
                        t_start=t_start * pq.ms,
                        name=file_blocks[list_chan[ind_chan]]['m_Name'],
                        file_origin=os.path.basename(self.filename),
                        channel_index=int(file_blocks
                            [list_chan[ind_chan]]['m_numChannel']))
                    
                    if lazy:
                        ana_sig.lazy_shape = chan_len[ind_chan]
                    self._annotate_block(
                        ana_sig,
                        file_blocks[list_chan[ind_chan]])

                    seg.analogsignals.append(ana_sig)

                #=====================================
                # LEVEL signal with SPIKETRAIN option
                #=====================================
                if (type_signal == 'level' and
                     self.option_read_level == 'SpikeTrain'):

                    if lazy:
                        times = np.array([])
                        t_stop = 0
                        waveforms = np.array([])

                    else:
                        amplitude = file_blocks[
                            list_chan[ind_chan]]['m_Amplitude']
                        (times, t_start, t_stop, waveforms) = \
                            self._read_spike_data(
                                fid,
                                list_blocks,
                                file_blocks,
                                sampling_rate,
                                amplitude)
                    
                    if self.option_read_waveform:                    
                        spike_train = SpikeTrain(
                            times * pq.ms,
                            sampling_rate=sampling_rate * pq.kHz,
                            t_start=t_start * pq.ms,
                            t_stop=t_stop * pq.ms,
                            waveforms=waveforms * pq.V,
                            left_sweep=file_blocks[
                                list_chan[ind_chan]]['m_nPreTrigmSec'] * pq.ms,
                            name=file_blocks[list_chan[ind_chan]]['m_Name'],
                            file_origin=os.path.basename(self.filename))
                    else:
                        spike_train = SpikeTrain(
                            times * pq.ms,
                            sampling_rate=sampling_rate * pq.kHz,
                            t_start=t_start * pq.ms,
                            t_stop=t_stop * pq.ms,
                            left_sweep=file_blocks[
                                list_chan[ind_chan]]['m_nPreTrigmSec'] * pq.ms,
                            name=file_blocks[list_chan[ind_chan]]['m_Name'],
                            file_origin=os.path.basename(self.filename))    
                        
                    

                    if lazy:
                        spike_train.lazy_shape = chan_len[ind_chan]

                    spike_train.channel_index = int(
                        file_blocks[list_chan[ind_chan]]['m_numChannel'])
                    self._annotate_block(
                        spike_train,
                        file_blocks[list_chan[ind_chan]])
                    seg.spiketrains.append(spike_train)

                #==================
                #  DIGITAL signal
                #==================
                elif type_signal == 'digital':
                    if lazy:
                        times_up = np.array([])
                        times_down = np.array([])
                        labels_up = np.array([], dtype='S')
                        labels_down = np.array([], dtype='S')

                    else:
                        (times_up,
                         times_down) = self._read_digital_data(
                            fid,
                            list_blocks,
                            file_blocks,
                            chan_len,
                            ind_chan,
                            ind,
                            sampling_rate)
                        labels_up = np.zeros(times_up.shape[0],
                                             dtype='S')
                        labels_down = np.zeros(times_down.shape[0],
                                               dtype='S')
                        
                    dig_sig_up = EventArray(
                        times_up * pq.ms,
                        labels=labels_up,
                        name=file_blocks[list_chan[ind_chan]]['m_Name'] \
                            + "_Up",
                        file_origin=os.path.basename(self.filename))
                    dig_sig_down = EventArray(
                        times_down * pq.ms,
                        labels=labels_down,
                        name=file_blocks[list_chan[ind_chan]]['m_Name'] \
                            + "_Down",
                        file_origin=os.path.basename(self.filename))
                       
                    if lazy:
                        dig_sig_up.lazy_shape = chan_len[ind_chan] \
                            - chan_len[ind_chan]//2
                        dig_sig_down.lazy_shape = chan_len[ind_chan]//2

                    self._annotate_block(
                        dig_sig_up,
                        file_blocks[list_chan[ind_chan]])
                    self._annotate_block(
                        dig_sig_down,
                        file_blocks[list_chan[ind_chan]])

                    dig_sig_up.annotate(sampling_rate=sampling_rate * pq.kHz)
                    dig_sig_down.annotate(sampling_rate=sampling_rate * pq.kHz)

                    seg.eventarrays.append(dig_sig_up)
                    seg.eventarrays.append(dig_sig_down)

                #==================
                #  PORT signal
                #==================
                elif type_signal == 'port':

                    if lazy:
                        times = np.array([])
                        labels = np.array([], dtype='S')

                    else:
                        # WARNING: sampling_rate is not present in
                        # the block b which refer to the channel (value of 0.0)
                        # the sampling_rate is present in block 2 for
                        # corresponding digital channels
                        sampling_rate = port_sampling_rate
                        times, labels = self._read_port_data(
                            fid,
                            list_blocks,
                            file_blocks,
                            chan_len,
                            ind_chan,
                            ind,
                            sampling_rate)

                    port_sig = EventArray(
                        times * pq.ms,
                        labels=labels,
                        name=file_blocks[list_chan[ind_chan]]['m_Name'],
                        file_origin=os.path.basename(self.filename))

                    self._annotate_block(
                        port_sig,
                        file_blocks[list_chan[ind_chan]])
                    port_sig.annotate(sampling_rate=sampling_rate * pq.kHz)
                    t_start /= sampling_rate
                    port_sig.annotate(t_start=t_start * pq.ms)
                    seg.eventarrays.append(port_sig)

        fid.close()

        #===============================
        # Step 7: annotate the Neo block
        #===============================

        # select group and subgroup with available data
        if cascade:
            groups = self._select_available_channel(groups,
                                                   list_available_channel)
            # annotations
            for info in information:
                blck.annotations[info] = information[info]
            blck.annotate(groups=groups)

            populate_RecordingChannel(blck, remove_from_annotation=False)

            create_many_to_one_relationship(blck)

        return blck

"""
Information for special types in [1]
------------------------------------

_dostime_t type definition:
struct dos_time_t
{
    unsigned char hour; /* hours (0-23) */
    unsigned char minute; /* minutes (0-59) */
    unsigned char second; /* seconds (0-59) */
    unsigned char hsecond; /* seconds/ 100 (0-99) */
}

_dosdate_t type definition:
struct _dosdate_t
{
    unsigned char day;       /* day of month( 1-31) */
    unsigned char month;     /* month (1-12) */
    unsigned int year;       /* year (1980-2099) */
    unsigned char dayofweek; /* day of week (0 = Sunday) */
}

WINDOWPLACEMENT16 type definition (according to WINE source code):
typedef struct
{
    UINT16   length;
    UINT16   flags;
    UINT16   showCmd;
    POINT16  ptMinPosition;
    POINT16  ptMaxPosition;
    RECT16   rcNormalPosition;
} WINDOWPLACEMENT16,*LPNONCLIENTMETRICS16;

POINT16 struct
{
    INT16 x;
    INT16 y;
}

RECT16 Struct
{
    INT16 bottom;
    INT16 left;
    INT16 right;
    INT16 top;
}

"""

max_string_len = '80s'

# maximal length of variable length strings in the file
# WARNING: I don't know what is the real value here. The maximum value actually
# found is 80 (the longest string is in the subgroup names, as it is built by
# combining the names of the channels in the subgroup, with up to 8 channels per
# subgroup).

# FIXME: A cleaner way to handle strings reading is suitable. Currently I
# read a buffer of max_string_len bytes and look for the C "end of string"
# character ('\x00'). It would be better either to read characters until
# reaching '\x00' or to read the exact number of characters needed, if the
# length of a string can be deduced from the lentgh of the block and the number
# of bytes already read (it seems possible, at least for certain block types).

# The name of the keys in the folowing dicts are chosen to match as closely as
# possible the names in document [1]

TypeH_Header = [
    ('m_nextBlock', 'l'),
    ('m_version', 'h'),
    ('m_time_hour', 'B'),
    ('m_time_minute', 'B'),
    ('m_time_second', 'B'),
    ('m_time_hsecond', 'B'),
    ('m_date_day', 'B'),
    ('m_date_month', 'B'),
    ('m_date_year', 'H'),
    ('m_date_dayofweek', 'B'),
    ('blank', 'x'),  # one byte blank for 2 bytes alignement
    ('m_minTime', 'd'),
    ('m_maxTime', 'd'),
    ('m_EraseCount', 'l'),
    ('m_mapVersion', 'b'),
    ('m_ApplicationName', '10s'),
    ('m_ResourceVersion', '4s'),  # WARNING: present only in version 65,
    ('blank2', 'x'),
    ('m_Reserved', 'l')]  # WARNING: Mpx version. Must be checked


Type0_SetBoards = [
    ('m_nextBlock', 'l'),
    ('m_BoardCount', 'h'),
    ('m_GroupCount', 'h'),
    ('m_MainWindow_length', 'H'),
    ('m_MainWindow_flags', 'H'),
    ('m_MainWindow_showCmd', 'H'),
    ('m_MainWindow_ptMinPosition_x', 'H'),
    ('m_MainWindow_ptMinPosition_y', 'H'),
    ('m_MainWindow_ptMaxPosition_x', 'H'),
    ('m_MainWindow_ptMaxPosition_y', 'H'),
    ('m_MainWindow_rcNormalPosition_bottom', 'H'),
    ('m_MainWindow_rcNormalPosition_left', 'H'),
    ('m_MainWindow_rcNormalPosition_right', 'H'),
    ('m_MainWindow_rcNormalPosition_top', 'H')]  # WARNING: order unknown

Type1_Boards = [
    ('m_nextBlock', 'l'),
    ('m_Number', 'h'),
    ('m_countChannel', 'h'),
    ('m_countAnIn', 'h'),
    ('m_countAnOut', 'h'),
    ('m_countDigIn', 'h'),
    ('m_countDigOut', 'h'),
    ('m_TrigCount', 'h'),  # not defined in 5.3.3 but appears in 5.5.1 and
                           # seems to really exist in files
    ('m_Amplitude', 'f'),
    ('m_cSampleRate', 'f'),  # sample rate seems to be given in kHz
    ('m_Duration', 'f'),
    ('m_nPreTrigmSec', 'f'),
    ('m_nPostTrigmSec', 'f'),
    ('m_TrgMode', 'h'),
    ('m_LevelValue', 'h'),  # after this line, 5.3.3 is wrong,
                            # check example in 5.5.1 for the right fields
    ('m_nSamples', 'h'),
    ('m_LevelFactorRMS', 'f'),
    ('m_ScaleFactor', 'f'),
    ('m_DapTime', 'f'),
    ('m_nameBoard', '35s'),
    ('m_DiscMaxValue', 'b'),  # WARNING: should this exist?
    ('m_DiscMinValue', 'b'),  # WARNING: should this exist?
    ('blank', 'x')]  # WARNING : may be before m_DiscMaxValue or after.

Type2_DefBlocksChannels = [
    # common parameters for all types of channels
    ('m_nextBlock', 'l'),
    ('m_isAnalog', 'h'),
    ('m_isInput', 'h'),
    ('m_numChannel', 'h'),
    ('m_numColor', 'h')]

Type2_SubBlockAnalogChannels = [
    ('blank', '2x'),  # WARNING : this is not in the specs but it seems needed
    ('m_Mode', 'h')]

Type2_SubBlockContinuousChannels = [
    # continuous channels parameters
    ('m_Amplitude', 'f'),
    ('m_SampleRate', 'f'),
    ('m_ContBlkSize', 'h'),
    ('m_ModeSpike', 'h'),  # WARNING: the C code [3] uses usigned short here
    ('m_Duration', 'f'),
    ('m_bAutoScale', 'h'),
    ('m_Name', max_string_len)]

Type2_SubBlockLevelChannels = [
    # level channels parameters
    ('m_Amplitude', 'f'),
    ('m_SampleRate', 'f'),
    ('m_nSpikeCount', 'h'),
    ('m_ModeSpike', 'h'),
    ('m_nPreTrigmSec', 'f'),
    ('m_nPostTrigmSec', 'f'),
    ('m_LevelValue', 'h'),
    ('m_TrgMode', 'h'),
    ('m_YesRMS', 'h'),
    ('m_bAutoScale', 'h'),
    ('m_Name', max_string_len)]

Type2_SubBlockExtTriggerChannels = [  # WARNING: untested
    # external trigger channels parameters
    ('m_Amplitude', 'f'),
    ('m_SampleRate', 'f'),
    ('m_nSpikeCount', 'h'),
    ('m_ModeSpike', 'h'),
    ('m_nPreTrigmSec', 'f'),
    ('m_nPostTrigmSec', 'f'),
    ('m_TriggerNumber', 'h'),
    ('m_Name', max_string_len)]

Type2_SubBlockDigitalChannels = [
    # digital channels parameters
    ('m_Mode', 'h'),
    ('m_SampleRate', 'f'),
    ('m_SaveTrigger', 'h'),
    ('m_Duration', 'f'),
    ('m_PreviousStatus', 'h'),
    ('m_Name', max_string_len)]

Type2_SubBlockUnknownChannels = [
    # WARNING: We have a mode that doesn't appear in our spec, so we don't
    # know what are the fields.
    # It seems that for non-digital channels the beginning is
    # similar to continuous channels. Let's hope we're right...
    ('blank', '2x'),
    ('m_Amplitude','f'),
    ('m_SampleRate','f')]
    # there are probably other fields after...

Type6_DefBlockTrigger = [
    ('m_nextBlock', 'l'),
    ('m_Number', 'h'),
    ('m_countChannel', 'h'),
    ('m_StateChannels', 'h'),
    ('m_numChannel1', 'h'),
    ('m_numChannel2', 'h'),
    ('m_numChannel3', 'h'),
    ('m_numChannel4', 'h'),
    ('m_numChannel5', 'h'),
    ('m_numChannel6', 'h'),
    ('m_numChannel7', 'h'),
    ('m_numChannel8', 'h'),
    ('m_Name', max_string_len)]

Type3_DefBlockGroup = [
    ('m_nextBlock', 'l'),
    ('m_Number', 'h'),
    ('m_Z_Order', 'h'),
    ('m_countSubGroups', 'h'),
    ('m_MainWindow_length', 'H'),
    ('m_MainWindow_flags', 'H'),
    ('m_MainWindow_showCmd', 'H'),
    ('m_MainWindow_ptMinPosition_x', 'H'),
    ('m_MainWindow_ptMinPosition_y', 'H'),
    ('m_MainWindow_ptMaxPosition_x', 'H'),
    ('m_MainWindow_ptMaxPosition_y', 'H'),
    ('m_MainWindow_rcNormalPosition_bottom', 'H'),
    ('m_MainWindow_rcNormalPosition_left', 'H'),
    ('m_MainWindow_rcNormalPosition_right', 'H'),
    ('m_MainWindow_rcNormalPosition_top', 'H'),  # WARNING: order unknown
    ('m_NetLoc', 'h'),
    ('m_locatMax', '4I'),
    ('m_nameGroup', max_string_len)]  # 'c' in documentation

Type4_DefBlockSubgroup = [
    ('m_nextBlock', 'l'),
    ('m_Number', 'h'),
    ('m_TypeOverlap', 'h'),
    ('m_Z_Order', 'h'),
    ('m_countChannel', 'h'),
    ('m_NetLoc', 'h'),
    ('m_location', '4I'),
    ('m_bIsMaximized', 'h'),
    ('m_numChannel1', 'H'),
    ('m_numChannel2', 'H'),
    ('m_numChannel3', 'H'),
    ('m_numChannel4', 'H'),
    ('m_numChannel5', 'H'),
    ('m_numChannel6', 'H'),
    ('m_numChannel7', 'H'),
    ('m_numChannel8', 'H'),
    ('m_Name', max_string_len)]

Type5_DataBlockOneChannel = [
    ('m_numChannel', 'h')]
    # WARNING: 'm_numChannel' (called 'm_Number' in 5.4.1 of [1]) is supposed
    # to be uint according to 5.4.1 but it seems to be a short in the files
    # (or should it be ushort ?)

# WARNING: In 5.1.1 page 121 of [1], they say "Note: 5 is used for demo
# purposes, 7 is used for real data", but looking at some real datafiles,
# it seems that block of type 5 are also used for real data...

Type7_DataBlockMultipleChannels = [
    ('m_numBoard', 'H'),  # WARNING: unknown true type
    ('FINT', 'h')]

Type7_DataSubblockDAPBuffers = [
    # FINT = - 111
    ('DAP_buffers_filling', 'h')]

Type7_DataSubblockRMS = [  # WARNING: untested
    # FINT = - 222
    ('RMS_value', 'b')]

Type7_DataSubblockRestart = [  # WARNING: untested
    # FINT = - 444
    ('num_channel', 'h')]

Type7_DataSubblockDataLoss = [
    # FINT = - 333
    ('num_channel', 'H'),
    ('first_lost_sample', 'i'),
    ('last_lost_sample', 'i')]

Type7_DataSubblockStartAnalogOutput = [  # WARNING: untested
    # FINT = - 666
    ]

Type7_DataSubblockStopAnalogOutput = [  # WARNING: untested
    # FINT = - 557 or - 558
    ('sample_count', 'i')]

Type7_DataSubblockUnknown = [  # WARNING: untested
    ('value', 'h')]

TypeP_DefBlockPeriStimHist = [  # WARNING: present in version 100
    ('m_numChannel', 'h'),
    ('m_Position', '4I'),
    ('m_isStatVisible', 'h'),
    ('m_DurationSec', 'f'),
    ('m_Rows', 'i'),
    ('m_DurationSecPre', 'f'),
    ('m_Bins', 'i'),
    ('m_NoTrigger', 'h')]

TypeF_DefBlockFRTachogram = [  # WARNING: untested
    ('m_numChannel', 'h'),
    ('m_Position', '4I'),
    ('m_isStatVisible', 'h'),
    ('m_DurationSec', 'f'),
    ('m_AutoManualScale', 'i'),
    ('m_Max', 'i')]

TypeR_DefBlockRaster = [
    ('m_numChannel', 'h'),
    ('m_Position', '4I'),
    ('m_isStatVisible', 'h'),
    ('m_DurationSec', 'f'),
    ('blank', '2x'),  # WARNING: position not sure
    ('m_Rows', 'i'),
    ('m_NoTrigger', 'h')]

TypeI_DefBlockISIHist = [  # WARNING: untested
    ('m_numChannel', 'h'),
    ('m_Position', '4I'),
    ('m_isStatVisible', 'h'),
    ('m_DurationSec', 'f'),
    ('m_Bins', 'i'),
    ('m_TypeScale', 'i')]

Type8_MarkerBlock = [  # WARNING: untested
    ('m_numChannel', 'h'),
    ('m_Time', 'l')]  # WARNING: check what's the right type here.
    # It seems that the size of time_t type depends on the system typedef,
    # I put long here but I couldn't check if it is the right type

Type9_ScaleBlock = [  # WARNING: untested
    ('m_numChannel', 'h'),
    ('m_Scale', 'f')]

Typeb_DigPortDef = [
    ('m_BoardNumber', 'i'),
    ('m_numChannel', 'h'),
    ('m_SampleRate', 'f'),
#    ('m_PrevValue','H'),    # WARNING : seems absent in file
    ('m_Name', max_string_len)]

TypeS_DefStreamData = [  # WARNING: untested
    ('m_nextBlock', 'l'),
    ('m_Number', 'h'),
    ('m_SampleRate', 'f'),
    ('m_Name', max_string_len)]

TypeE_StreamDataBlock = [  # WARNING: untested
    ('m_TypeEvent', 'c'),
    ('m_uTimeStamp', 'L'),
    ('m_StreamData', max_string_len)]

Type_Unknown = []

dict_header_type = {
    'h': TypeH_Header,
    '0': Type0_SetBoards,
    '1': Type1_Boards,
    '2': Type2_DefBlocksChannels,
    '6': Type6_DefBlockTrigger,
    '3': Type3_DefBlockGroup,
    '4': Type4_DefBlockSubgroup,
    '5': Type5_DataBlockOneChannel,
    '7': Type7_DataBlockMultipleChannels,
    'P': TypeP_DefBlockPeriStimHist,
    'F': TypeF_DefBlockFRTachogram,
    'R': TypeR_DefBlockRaster,
    'I': TypeI_DefBlockISIHist,
    '8': Type8_MarkerBlock,
    '9': Type9_ScaleBlock,
    'b': Typeb_DigPortDef,
    'S': TypeS_DefStreamData,
    'E': TypeE_StreamDataBlock,
    }

info_quantities = {
    "m_minTime": pq.s,
    "m_maxTime": pq.s,
    "m_Amplitude": pq.V,
    "m_cSampleRate": pq.kHz,
    "m_Duration": pq.ms,
    "m_nPreTrigmSec": pq.ms,
    "m_nPostTrigmSec": pq.ms,
    "m_LevelValue": pq.V,
    }


class HeaderReader():
    '''
    Class for metadata reading in AlphaOmega.map files

    '''

    def __init__(self, fid, description):
        '''
        :param file fid: file ready for reading
        :param dict description: dictionnary with :
            * key (string) = value_name
            * value (string) = C type format reading

        '''
        self.fid = fid
        self.description = description

    def read_f(self, offset=None):
        '''
        Method for reading file with unpack method

        :param int offset: DEFAULT = None.
        :return dict d: dictionnary with metadata

        '''
        if offset is not None:
            self.fid.seek(offset)

        d = {}  # dict with metadata
        for key, fmt in self.description:
            fmt = '<' + fmt  # insures use of standard sizes
            buf = self.fid.read(struct.calcsize(fmt))
            # value reading
            if len(buf) != struct.calcsize(fmt):
                return None
            val = list(struct.unpack(fmt, buf))
            # value treatment
            # NOTE: in Python 3, the string data are read in bytes format.
            for i, ival in enumerate(val):
                if hasattr(ival, 'split'):
                    try:  # python 2
                        val[i] = ival.split('\x00')[0]
                    except:  # python 3
                        val[i] = ival.split(b'\x00')[0].decode("utf-8")
            if len(val) == 1:
                val = val[0]
            if isinstance(val, bytes):
                val = val.decode("utf-8")
            # NOTE: in python 2, some value are read and return in unicode
            # unless str format. Neo object can't take unicode object, so
            # they are convert in str format.
            try:  # python 2
                if isinstance(val, unicode):
                    val = str(val)
            except:  # python 3
                pass
            # value storage
            d[key] = val
        return d
