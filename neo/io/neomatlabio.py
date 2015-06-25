# -*- coding: utf-8 -*-
"""
Module for reading/writing Neo objects in MATLAB format (.mat) versions 5 to
7.2.

This module is a bridge for MATLAB users who want to adopt the Neo object
representation. The nomenclature is the same but using Matlab structs and
cell arrays. With this module MATLAB users can use neo.io to read a format and
convert it to .mat.

=============================================================================
Supported : Read/Write

Author: sgarcia, Lea Andre
"""

from datetime import datetime
from time import time
from distutils import version
import re

import numpy as np
import quantities as pq

import copy

# check scipy
try:
    import scipy.io
    import scipy.version
except ImportError as err:
    HAVE_SCIPY = False
    SCIPY_ERR = err
else:

    if version.LooseVersion(scipy.version.version) < '0.8':
        HAVE_SCIPY = False
        SCIPY_ERR = ImportError("""your scipy version is too old to support
                                MatlabIO, you need at least 0.8.
                                You have %s""" % (scipy.version.version))
    else:
        HAVE_SCIPY = True
        SCIPY_ERR = None


from neo.io.baseio import BaseIO
from neo.core import (Block, Segment, RecordingChannelGroup, RecordingChannel,
                      Unit, AnalogSignal, AnalogSignalArray, Event, EventArray,
                      Spike, SpikeTrain, Epoch, EpochArray,
                      IrregularlySampledSignal)
from neo.io.tools import create_many_to_one_relationship
from neo import description


classname_lower_to_upper = {}
for k in description.class_by_name.keys():
    classname_lower_to_upper[k.lower()] = k


class NeoMatlabIO(BaseIO):
    """
    Class for reading/writing Neo objects in MATLAB format (.mat) versions 5 to
    7.2.

    This module is a bridge for MATLAB users who want to adopt the Neo object
    representation. The nomenclature is the same but using Matlab structs and
    cell arrays. With this module MATLAB users can use neo.io to read a format
    and convert it to .mat.

    Rules of conversion:
      * Neo classes are converted to MATLAB structs.
        e.g., a Block is a struct with attributes "name", "file_datetime", ...
      * Neo one_to_many relationships are cellarrays in MATLAB.
        e.g., ``seg.analogsignals[2]`` in Python Neo will be
        ``seg.analogsignals{3}`` in MATLAB.
      * Quantity attributes are represented by 2 fields in MATLAB.
        e.g., ``anasig.t_start = 1.5 * s`` in Python will be
        ``anasig.t_start = 1.5`` and ``anasig.t_start_unit = 's'`` in MATLAB.
      * classes that inherit from Quantity (AnalogSignal, SpikeTrain, ...) in
        Python will have 2 fields (array and units) in the MATLAB struct.
        e.g.: ``AnalogSignal( [1., 2., 3.], 'V')`` in Python will be
        ``anasig.array = [1. 2. 3]`` and ``anasig.units = 'V'`` in MATLAB.

    1 - **Scenario 1: create data in MATLAB and read them in Python**

        This MATLAB code generates a block::

            block = struct();
            block.segments = { };
            block.name = 'my block with matlab';
            for s = 1:3
                seg = struct();
                seg.name = strcat('segment ',num2str(s));
                seg.analogsignals = { };
                for a = 1:5
                    anasig = struct();
                    anasig.array = rand(100,1);
                    anasig.units = 'mV';
                    anasig.t_start = 0;
                    anasig.t_start_units = 's';
                    anasig.sampling_rate = 100;
                    anasig.sampling_rate_units = 'Hz';
                    seg.analogsignals{a} = anasig;
                end
                seg.spiketrains = { };
                for t = 1:7
                    sptr = struct();
                    sptr.array = rand(30,1)*10;
                    sptr.units = 'ms';
                    sptr.t_start = 0;
                    sptr.t_start_units = 'ms';
                    sptr.t_stop = 10;
                    sptr.t_stop_units = 'ms';
                    seg.spiketrains{t} = sptr;
                end

                block.segments{s} = seg;
            end
            save 'myblock.mat' block -V7


        This code reads it in Python::

            import neo
            r = neo.io.NeoMatlabIO(filename='myblock.mat')
            bl = r.read_block()
            print bl.segments[1].analogsignals[2]
            print bl.segments[1].spiketrains[4]


    2 - **Scenario 2: create data in Python and read them in MATLAB**

        This Python code generates the same block as in the previous scenario:

            import neo
            import quantities as pq
            from scipy import rand

            bl = neo.Block(name='my block with neo')
            for s in range(3):
                seg = neo.Segment(name='segment' + str(s))
                bl.segments.append(seg)
                for a in range(5):
                    anasig = neo.AnalogSignal(
                        rand(100),
                        units='mV',
                        t_start=0*pq.s,
                        sampling_rate=100*pq.Hz)
                    seg.analogsignals.append(anasig)
                for t in range(7):
                    sptr = neo.SpikeTrain(
                        rand(30),
                        units='ms',
                        t_start=0*pq.ms,
                        t_stop=10*pq.ms)
                    seg.spiketrains.append(sptr)

        w = neo.io.NeoMatlabIO(filename='myblock.mat')
        w.write_block(bl)


        This MATLAB code reads it::

            load 'myblock.mat'
            block.name
            block.segments{2}.analogsignals{3}.array
            block.segments{2}.analogsignals{3}.units
            block.segments{2}.analogsignals{3}.t_start
            block.segments{2}.analogsignals{3}.t_start_units


    3 - **Scenario 3: conversion**

        This Python code converts a Spike2 file to MATLAB::

            from neo import Block
            from neo.io import Spike2IO, NeoMatlabIO

            r = Spike2IO(filename='myspike2file.smr')
            w = NeoMatlabIO(filename='convertedfile.mat')
            seg = r.read_segment()
            bl = Block(name='a block')
            bl.segments.append(seg)
            w.write_block(bl)

    """
    is_readable = True
    is_writable = True

    supported_objects = [Block, Segment, RecordingChannelGroup,
                         RecordingChannel, Unit, AnalogSignal,
                         AnalogSignalArray, Event, EventArray, Spike,
                         SpikeTrain, Epoch, EpochArray,
                         IrregularlySampledSignal]
    readable_objects = [Block, ]
    writeable_objects = [Block, ]

    has_header = False
    is_streameable = False
    read_params = {Block: []}
    write_params = {Block: []}

    name = 'neomatlab'
    extensions = ['mat']

    mode = 'file'

    def __init__(self, filename=None):
        """
        This class reads/writes neo objects in matlab 5 to 7.2 format.

        :param str filename : the filename to read

        """
        if not HAVE_SCIPY:
            raise SCIPY_ERR
        BaseIO.__init__(self)
        self.filename = filename
        
        # FIXME: Some types allowed in annotations are not yet supported below.
        # Their support should be added.
        # For details of allowed types, see the variable 
        # ALLOWED_ANNOTATION_TYPES and the function _check_annotations() in
        # neo/core/baseneo.py.       
        
        self.available_types = [int, float, complex, pq.Quantity, str, time,
                                datetime]
        self.available_containers = [list, dict, np.ndarray]

    def read_block(self, cascade=True, lazy=False):
        """
        Method to read Neo Block in a file

        :param boolean lazy: if True, the numerical data is not loaded, only
            properties and metadata
        :param boolean cascade: if False, only a single object is loaded,
            without its references to other objects

        """
        d = scipy.io.loadmat(self.filename,
                             struct_as_record=False,
                             squeeze_me=True)
        assert 'block' in d, 'no block in' + self.filename
        bl_struct = d['block']
        bl = self.create_ob_from_struct(bl_struct,
                                        'Block',
                                        cascade=cascade,
                                        lazy=lazy)
        
        create_many_to_one_relationship(bl)

        
        
                        
        return bl

    def write_block(self, bl):
        """
        Method to write Neo Block in a file

        :param Block bl : the block to save

        """
        bl_struct = self.create_struct_from_obj(bl)

        for seg in bl.segments:
            seg_struct = self.create_struct_from_obj(seg)
            bl_struct['segments'].append(seg_struct)
            for anasig in seg.analogsignals:
                anasig_struct = self.create_struct_from_obj(anasig)
                seg_struct['analogsignals'].append(anasig_struct)

            for anasigarray in seg.analogsignalarrays:
                anasigarray_struct = self.create_struct_from_obj(anasigarray)
                seg_struct['analogsignalarrays'].append(anasigarray_struct)

            for event in seg.events:
                event_struct = self.create_struct_from_obj(event)
                seg_struct['events'].append(event_struct)

            for ea in seg.eventarrays:
                ea_struct = self.create_struct_from_obj(ea)
                seg_struct['eventarrays'].append(ea_struct)

            for sp in seg.spikes:
                sp_struct = self.create_struct_from_obj(sp)
                seg_struct['spikes'].append(sp_struct)

            for sptr in seg.spiketrains:
                sptr_struct = self.create_struct_from_obj(sptr)
                seg_struct['spiketrains'].append(sptr_struct)

            for ep in seg.epochs:
                ep_struct = self.create_struct_from_obj(ep)
                seg_struct['epochs'].append(ep_struct)

            for eparray in seg.epocharrays:
                eparray_struct = self.create_struct_from_obj(eparray)
                seg_struct['epocharrays'].append(eparray_struct)

            for iss in seg.irregularlysampledsignals:
                iss_struct = self.create_struct_from_obj(iss)
                seg_struct['irregularlysampledsignals'].append(iss_struct)

        for recChGr in bl.recordingchannelgroups:
            recChGr_struct = self.create_struct_from_obj(recChGr)
            bl_struct['recordingchannelgroups'].append(recChGr_struct)

            for reCh in recChGr.recordingchannels:
                reCh_struct = self.create_struct_from_obj(reCh)
                recChGr_struct['recordingchannels'].append(reCh_struct)

                # NOTE: recordingchannelgroups are present, but we don't take
                # them. It makes a infinity circle of recordingchannel and
                # recordingchannelgroups.
                # FIXME: a good solution to handle many to many relationship 
                # should be found, but it is not simple as MATLAB doesn't really 
                # support references to objects. 

#                for recChGr in reCh.recordingchannelgroups:
#                    recChGr_struct = self.create_struct_from_obj(recChGr)
#                    reCh_struct['recordingchannelgroups'].append(recChGr_struct)

                # NOTE: Seems analogsignals and irregularlysampledsignals are
                # not present in description.one_to_many_relationship in
                # RecordingChannel object.

#                for anasig in reCh.analogsignals:
#                    anasig_struct = self.create_struct_from_obj(anasig)
#                    reCh_struct['analogsignals'].append(anasig_struct)
#
#                for iss in reCh.irregularlysampledsignals:
#                    iss_struct = self.create_struct_from_obj(iss)
#                    reCh_struct['irregularlysampledsignals'].append(iss_struct)

            for unit in recChGr.units:
                unit_struct = self.create_struct_from_obj(unit)
                recChGr_struct['units'].append(unit_struct)

                for sp in unit.spikes:
                    sp_struct = self.create_struct_from_obj(sp)
                    unit_struct['spikes'].append(sp_struct)

                for sptr in unit.spiketrains:
                    sptr_struct = self.create_struct_from_obj(sptr)
                    unit_struct['spiketrains'].append(sptr_struct)

            for anasigarray in recChGr.analogsignalarrays:
                anasigarray_struct = self.create_struct_from_obj(anasigarray)
                recChGr_struct['analogsignalarrays'].append(anasigarray_struct)

        scipy.io.savemat(self.filename,
                         {'block': bl_struct},
                         long_field_names=True,
                         do_compression=True,
                         oned_as='row',)

    def create_struct_from_obj(self, ob):
        """
        Method to create a dictionnary (writable for matlab format) from a Neo
        object

        :param Neo objet ob : The Neo object to transform in struct format
        :return struct struct : the struct of the neo object

        """
        classname = ob.__class__.__name__
        struct = {}

        # relationship
        rel = copy.deepcopy(description.one_to_many_relationship)
        # NOTE: The update method is switching some one_to_many_relationship
        # information. To addition dictionnary, we use an other method.
        for cl in description.many_to_many_relationship:
            if cl in rel:
                for neoobj in description.many_to_many_relationship[cl]:
                    if neoobj not in rel[cl]:
                        rel[cl].append(neoobj)
            else:
                rel[cl] = description.many_to_many_relationship[cl]

        if classname in rel:
            for childname in rel[classname]:
                if description.class_by_name[
                    childname] in self.supported_objects:
                        struct[childname.lower() + 's'] = []
        
        # TODO: below only the necessary and recommend attributes are handled,
        # it would be better to handle all compatible attributes of the objects
        # even if they are not in the necessary and recommend attributes lists.

        # attributes
        necess = description.classes_necessary_attributes[classname]
        recomm = description.classes_recommended_attributes[classname]
        attributes = necess + recomm

        for i, attr in enumerate(attributes):

            attrname, attrtype = attr[0], attr[1]

            #~ if attrname =='':
                #~ struct['array'] = ob.magnitude
                #~ struct['units'] = ob.dimensionality.string
                #~ continue

            if (classname in description.classes_inheriting_quantities and
                description.classes_inheriting_quantities[
                    classname] == attrname):
                        struct[attrname] = ob.magnitude
                        struct[attrname + '_units'] = ob.dimensionality.string
                        continue

            if not(attrname in ob.annotations or hasattr(ob, attrname)):
                continue

            if getattr(ob, attrname) is None:
                continue

            if attrtype == pq.Quantity:
                #ndim = attr[2]
                struct[attrname] = getattr(ob, attrname).magnitude
                struct[attrname +
                    '_units'] = getattr(ob, attrname).dimensionality.string

            elif attrtype == datetime:
                struct[attrname] = str(getattr(ob, attrname))

            else:
                struct[attrname] = getattr(ob, attrname)

        # annotations
        if len(ob.annotations) != 0:
            dict_annot = ob.annotations
            struct['annotations'] = self.create_struct_annot(dict_annot)

        return struct

    def create_ob_from_struct(self, struct, classname, cascade=True,
                              lazy=False):
        """
        Method to create a Neo object from a mat struct format

        :param struct struct : structure from which we get the Neo object
        :param string classname : classname of the Neo object
        :param boolean lazy: if True, the numerical data is not loaded, only
            properties and metadata
        :param boolean cascade: if False, only a single object is loaded,
            without its references to other objects

        :return Neo object ob: the Neo object to return from struct

        """
        
        cl = description.class_by_name[classname]
        
        if classname in description.classes_inheriting_quantities:
            quantity_attr = description.classes_inheriting_quantities[
                classname]
            arr = getattr(struct, quantity_attr)
            data_complement = dict(
                units=str(getattr(struct, quantity_attr + '_units')))

            if "sampling_rate" in (at[0] for at in \
            description.classes_necessary_attributes[classname]):
                # FIXME: put fake value for now, should be replaced by correct 
                # value
                data_complement["sampling_rate"] = 0 * pq.kHz

            if "t_stop" in (at[0] for at in \
            description.classes_necessary_attributes[classname]):
                if len(arr) > 0:
                    data_complement["t_stop"] = arr.max()
                else:
                    data_complement["t_stop"] = 0.0

            if "t_start" in (at[0] for at in \
            description.classes_necessary_attributes[classname]):
                if len(arr) > 0:
                    data_complement["t_start"] = arr.min()
                else:
                    data_complement["t_start"] = 0.0

            if lazy:
                ob = cl([], **data_complement)
                ob.lazy_shape = arr.shape
            else:
                ob = cl(arr, **data_complement)
        else:
            ob = cl()

        for attrname in struct._fieldnames:

            # check children
            rel = copy.deepcopy(description.one_to_many_relationship)
            if (classname in rel and
                    attrname[:-1] in [r.lower() for r in rel[classname]]):
                try:
                    for c in range(len(getattr(struct, attrname))):
                        if cascade:
                            child = self.create_ob_from_struct(
                                getattr(struct, attrname)[c],
                                classname_lower_to_upper[attrname[:-1]],
                                cascade=cascade, lazy=lazy)
                            getattr(ob, attrname.lower()).append(child)

                except TypeError:
                    # strange behavior in scipy.io:
                    # if len is 1 so there is no len()
                    if cascade:
                        child = self.create_ob_from_struct(
                            getattr(struct, attrname),
                            classname_lower_to_upper[attrname[:-1]],
                            cascade=cascade, lazy=lazy)
                        getattr(ob, attrname.lower()).append(child)

                continue

####

            # check children
            rel = description.many_to_many_relationship
            
            if (classname in rel and
                attrname[:-1] in [r.lower() for r in rel[classname]]):
                    try:
                        for c in range(len(getattr(struct, attrname))):
                            if cascade:
                                child = self.create_ob_from_struct(
                                    getattr(struct, attrname)[c],
                                    classname_lower_to_upper[attrname[:-1]],
                                    cascade=cascade, lazy=lazy)
                                getattr(ob, attrname.lower()).append(child)

                    except TypeError:
                        # strange behavior in scipy.io:
                        # if len is 1 so there is no len()
                        if cascade:
                            child = self.create_ob_from_struct(
                                getattr(struct, attrname),
                                classname_lower_to_upper[attrname[:-1]],
                                cascade=cascade, lazy=lazy)
                            getattr(ob, attrname.lower()).append(child)

                    continue


#####
            # attributes
            if attrname.endswith('_units') or attrname == 'units':
                # linked with another field
                continue

            if classname in description.classes_inheriting_quantities and \
                    description.classes_inheriting_quantities[classname] == \
                    attrname:
                continue
            item = getattr(struct, attrname)

            # put the good type
            necess = description.classes_necessary_attributes[classname]
            recomm = description.classes_recommended_attributes[classname]
            attributes = necess + recomm

            dict_attributes = dict([(a[0], a[1:]) for a in attributes])

            # retrieve data necess and recomm
            if attrname in dict_attributes:
                attrtype = dict_attributes[attrname][0]

                if attrtype == datetime:
                    m = '(\d+)-(\d+)-(\d+) (\d+):(\d+):(\d+).(\d+)'
                    r = re.findall(m, str(item))
                    if len(r) == 1:
                        item = datetime(*[int(e) for e in r[0]])
                    else:
                        item = None

                elif attrtype == np.ndarray:
                    dt = dict_attributes[attrname][2]
                    if lazy:
                        item = np.array([], dtype=dt)
                        ob.lazy_shape = item.shape
                    else:
                        try:
                            item = item.astype(dt)
                        except:
                            item = item

                elif attrtype == pq.Quantity:
                    ndim = dict_attributes[attrname][1]
                    units = str(getattr(struct, attrname + '_units'))
                    if ndim == 0:
                        item = pq.Quantity(item, units)
                    else:
                        if lazy:
                            item = pq.Quantity([], units)
                            item.lazy_shape = item.shape
                        else:
                            item = pq.Quantity(item, units)

                else:
                    item = attrtype(item)

            # retrieve annotations
            if attrname == 'annotations':
                item = self.create_annot_dict_from_struct(item)

            setattr(ob, attrname, item)

        return ob

    def create_annot_dict_from_struct(self, item):
        """
        Method to create a annotation dictionnary from a mat struct

        :param struct item : a mat object read in a matlab file
        :return dict obj_dict : the dict annotations

        """
        obj_dict = {}

        for attrname in item._fieldnames:
            obj = None
            sub_item = getattr(item, attrname)

            if ((type(sub_item) in self.available_types or
                 type(sub_item) in self.available_containers) and
                (not attrname.endswith('_units') or
                attrname != 'units')):
                    obj = sub_item

            # tuple and list
            if (isinstance(sub_item, list) or isinstance(sub_item, tuple)):
                obj = []
                for s_sub_item in sub_item:
                    if type(sub_item) == scipy.io.matlab.mio5_params.mat_struct:
                        s_sub_item = self.create_annot_dict_from_struct(
                            s_sub_item)
                    obj.append(s_sub_item)

            # dict
            if isinstance(sub_item, dict):
                obj = {}
                for key, value in sub_item:
                    if type(sub_item) == scipy.io.matlab.mio5_params.mat_struct:
                           value = self.create_annot_dict_from_struct(
                                value)
                    obj[key] = value

            # array
            if isinstance(sub_item, np.ndarray):
                obj_shape = sub_item.shape
                try:  # Case 1 : array is not a scalar
                    dim_total = 1
                    for dim in sub_item.shape:
                        dim_total *= dim
                    sub_item = sub_item.reshape((dim_total))
                    obj = np.empty(sub_item.shape, dtype=np.object)
                    obj.fill(np.nan)
                except:  # Case 2 : array contains a scalar and shape is ()
                         # here we have a 0d array

                    obj = sub_item.item()
                else:
                    # NOTE: with reshape, we have a 1 dimension array, that
                    # will be reshape with good dimension after verification
                    for ind, value in enumerate(sub_item):
                        if type(sub_item) == \
                                scipy.io.matlab.mio5_params.mat_struct:
                            value = self.create_annot_dict_from_struct(value)
                        obj[ind] = value
                obj = obj.reshape(obj_shape)

            # mat_struct
            if type(sub_item) == scipy.io.matlab.mio5_params.mat_struct:
                obj = self.create_annot_dict_from_struct(sub_item)

            obj_dict[attrname] = obj

        new_obj_dict = {}

        # processing to get quantities
        for key in obj_dict:
            if key.endswith('_units'):
                unit = obj_dict[key]
                obj = obj_dict[key[:-6]]
                new_obj_dict[key[:-6]] = pq.Quantity(obj, unit)
            else:
                if key + '_units' not in obj_dict:
                    new_obj_dict[key] = obj_dict[key]

        return new_obj_dict

    def create_struct_annot(self, container):
        """
        Method to create an available struct (for matlab writing) from a 
        container

        :param container: (list, dict or np.ndarray)
        :return struct struct: the struct return from the container

        """
        struct = None

        try:
            type(container) in self.available_containers
        except:
            raise Exception("Invalid container %s" % (container))

        # LIST verification
        if isinstance(container, list):
            struct = []
            for value in container:
                # VALUE verification
                new_value = self.value_verif(value)
                struct.append(new_value)
            struct = np.array(struct, dtype=object)

        # DICT verification
        if isinstance(container, dict):
            # dtype contruction
            dtype = []
            for key in container:
                # KEY verification
                new_key = self.key_verif(key)
                dt = (new_key, np.object)
                dtype.append(dt)
                # Quantity case
                if isinstance(container[key], pq.Quantity):
                    dt = (new_key + '_units', np.object)
                    dtype.append(dt)
            if not dtype:
                dtype = np.object
            struct = np.empty([1], dtype=dtype)

            for key in container:
                # KEY verification
                new_key = self.key_verif(key)
                # VALUE verification
                new_value = self.value_verif(container[key])

                struct[0][new_key] = new_value
                # Quantity case
                if isinstance(container[key], pq.Quantity):
                    new_value_unit = new_value.dimensionality.string
                    struct[0][new_key + '_units'] = new_value_unit

        # NDARRAY verification
        if isinstance(container, np.ndarray):
            array_shape = container.shape
            total_dim = 1
            for dim in container.shape:
                total_dim *= dim
            # NOTE: with reshape, we have a 1 dimension array, that will be
            # reshape with good dimension after verification
            container = container.reshape((totaldim))
            struct = np.empty(container.shape)
            struct.fill(np.nan)
            for ind, value in enumerate(container):
                # VALUE verification
                new_value = self.value_verif(value)
                struct[ind] = new_value
            struct = struct.reshape(array_shape)

        return struct

    def value_verif(self, value):
        """
        Verification of value for mat struct creation

        :param value : value to verif
        :return value : value after verification

        """

        # 1 - value type must be in self.available_types
        try:
            type(value) in self.available_types
        except:
            raise Exception("Invalid value %s in dict" % (value))
        # 2 - if value is a container, verif and recup his struct
        if type(value) in self.available_containers:
            value = self.create_struct_annot(value)
        return value

    def key_verif(self, key):
        """
        Verification of dict key for mat struct creation

        4 rules for dict key:
            1. must be string type
            2. must start with a letter
            3. must only contain single word character : [a-z][A-Z][0-9] or '_'
            4. must have a correct length, less than ""keylengthmax""

        :param (strin, int) key : key from a dict to verif

        """
        keylengthmax = 63  # value of 'namelengthmax' (MATLAB variable)

        # 1 - key must be str
        if not isinstance(key, str):
            raise Exception("Invalid key %s in dict" % (key))
        # 2 - key must be MATLAB variables names format
        if isinstance(key, str):
            match_format = re.search('\A[a-zA-Z][\w]*', key)
            if not match_format:
                # don't start with letter
                match_first = re.search('\A[^a-zA-Z]', key)
                if match_first:
                    raise Exception("Invalid key %s in dict. " % (key) +
                                    "The key must start with a" + 
                                    "letter ")
                # invalid character inside
                match_invalid = re.search('[\W]', key)
                if match_invalid:
                    raise Exception("Invalid key %s in dict. " % (key) +
                                    "The key must only contain the following " +
                                    "characters: [a-z], [A-Z], [0-9] or '_' ")
                # else
                if not match_first and not match_invalid:
                    raise Exception("Invalid key %s in dict" % (key))
        # 3 - key must be less than keylengthmax
            try:
                len(key) <= keylengthmax
            except:
                raise Exception("Invalid key %s in dict. "% (key) +
                                "The length of the key must be lower " +
                                "than %i" % (keylengthmax))
        return key
