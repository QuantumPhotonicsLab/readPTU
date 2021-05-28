from __future__ import print_function
from __future__ import absolute_import

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  4 18:36:54 2017

@author: raphaelproux and guillemballesteros

Largely inspired from PicoQuant examples:
    https://github.com/PicoQuant/PicoQuant-Time-Tagged-File-Format-Demos
and (in particular, to read the header) from a jupyter notebook
by tritemio on GitHub:
    https://gist.github.com/tritemio/734347586bc999f39f9ffe0ac5ba0e66

Please note:
- this library has been tested only with T2 mode of a Hydraharp v2, but it
  should work with a Hydraharp v1 or a Picoharp file (written, untested),
- does not support T3 mode for now (not written), but should be easy to
  implement

"""

from matplotlib import pyplot as plt
import numpy as np
import os
import mmap
import struct
import time
import collections as coll

from _readTTTRRecords_HHT2_HH2 import ffi, lib
from _readTTTRRecords_HHT2_HH2 import lib as HHT2_HH2_lib
from _readTTTRRecords_HHT2_HH1 import lib as HHT2_HH1_lib
from _readTTTRRecords_PHT2 import lib as PHT2_lib

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx,array[idx]

class PTUfile():
    """
    PTUfile() handles Picoquant PTU files.

    It takes care of  opening and closing, and header reading to extract
    information like the record type, measurement time and number of records.

    More detailed analysis like timetraces and g2 are done using the
    PTUmeasurement class.

    Attributes:
        acq_time (int): measurement acquisition time in picoseconds
        c_filehandle (C FILE *): C file handle used for possible analysis in C
        c_rec_num (C uint64_t *): C pointer to record number (used when
            analysing the file in C)
        end_header_offset (int): Position in bytes of the beginning of the
            records section in the file
        filehandle (TYPE): Python filehandle used only to read the header
        FileTagEnd (str): Constant string, tag to detect end of header
        globres (float): resolution of timetags in seconds
        mm (mmap object): memory map object created to read the file header
        num_records (int): number of records in the file
        rec_type (dict): reference dictionary of record types identifiers
        rec_type_r (dict): reversed dictionary of rec_type (allows inverse
            referencing)
        record_type (int): string identifier of the record type for the current
            file. Note you should send rec_type[record_type] to use the int
            identifier (ie in C).
        tag_type (dict): reference dictionary of header tag types identifiers
        tag_type_r (dict): reversed dictionary of tag_type
        tags (ordered dict): ordered dictionary of the header tags for the
            current file
        TTTRTagAcquisitionTime (str): header tag string identifier for
            acquisition time
        TTTRTagGlobRes (str): header tag string identifier for global
            resolution
        TTTRTagNumRecords (str): header tag string identifier for number of
            records
        TTTRTagRes (str): header tag string identifier for delay time
            resolution (T3 mode)
        TTTRTagTTTRRecType (str): header tag string identifier for record type
    """

    # Constants
    TTTRTagTTTRRecType = "TTResultFormat_TTTRRecType"
    TTTRTagNumRecords  = "TTResult_NumberOfRecords"  # Number of TTTR Records in the File;
    TTTRTagRes         = "MeasDesc_Resolution"       # Resolution for the Dtime (T3 Only)
    TTTRTagGlobRes     = "MeasDesc_GlobalResolution"  # Global Resolution of TimeTag(T2) /NSync (T3)
    TTTRTagAcquisitionTime = "MeasDesc_AcquisitionTime"
    FileTagEnd         = "Header_End"                # Always appended as last tag (BLOCKEND)

    # Tag Types
    tag_type = dict(
        tyEmpty8      = 0xFFFF0008,
        tyBool8       = 0x00000008,
        tyInt8        = 0x10000008,
        tyBitSet64    = 0x11000008,
        tyColor8      = 0x12000008,
        tyFloat8      = 0x20000008,
        tyTDateTime   = 0x21000008,
        tyFloat8Array = 0x2001FFFF,
        tyAnsiString  = 0x4001FFFF,
        tyWideString  = 0x4002FFFF,
        tyBinaryBlob  = 0xFFFFFFFF)

    # Record Types
    rec_type = dict(
        rtPicoHarpT3     = 0x00010303,  # (SubID = $00 ,RecFmt: $01) (V1), T-Mode: $03 (T3), HW: $03 (PicoHarp)
        rtPicoHarpT2     = 0x00010203,  # (SubID = $00 ,RecFmt: $01) (V1), T-Mode: $02 (T2), HW: $03 (PicoHarp)
        rtHydraHarpT3    = 0x00010304,  # (SubID = $00 ,RecFmt: $01) (V1), T-Mode: $03 (T3), HW: $04 (HydraHarp)
        rtHydraHarpT2    = 0x00010204,  # (SubID = $00 ,RecFmt: $01) (V1), T-Mode: $02 (T2), HW: $04 (HydraHarp)
        rtHydraHarp2T3   = 0x01010304,  # (SubID = $01 ,RecFmt: $01) (V2), T-Mode: $03 (T3), HW: $04 (HydraHarp)
        rtHydraHarp2T2   = 0x01010204,  # (SubID = $01 ,RecFmt: $01) (V2), T-Mode: $02 (T2), HW: $04 (HydraHarp)
        rtTimeHarp260NT3 = 0x00010305,  # (SubID = $00 ,RecFmt: $01) (V1), T-Mode: $03 (T3), HW: $05 (TimeHarp260N)
        rtTimeHarp260NT2 = 0x00010205,  # (SubID = $00 ,RecFmt: $01) (V1), T-Mode: $02 (T2), HW: $05 (TimeHarp260N)
        rtTimeHarp260PT3 = 0x00010306,  # (SubID = $00 ,RecFmt: $01) (V1), T-Mode: $03 (T3), HW: $06 (TimeHarp260P)
        rtTimeHarp260PT2 = 0x00010206,  # (SubID = $00 ,RecFmt: $01) (V1), T-Mode: $02 (T2), HW: $06 (TimeHarp260P)
    )

    def __init__(self, filename):
        """
        Open the file and leave it opened for possible analysis with a
        PTUmeasurement object.

        Please note PTUfile supports context manager and it is strongly
        advised to create it using a "with ... as ..." statement.

        Args:
            filename (str): path + filename to the measurement file to open
        """
        # Reverse mappings of the tags and record dictionaries
        self.tag_type_r = {v: k for k, v in self.tag_type.items()}
        self.rec_type_r = {v: k for k, v in self.rec_type.items()}

        # var initialization
        self.tags = coll.OrderedDict()

#        print('Size of file:', os.path.getsize(filename))
        self.filename = filename
        self.filehandle = open(filename, 'rb')

        # with mmap, we will be able to use the file as a string
        self.mm = mmap.mmap(self.filehandle.fileno(),
                            0, access=mmap.ACCESS_READ)

#        magic = self.mm[:8].rstrip(b'\0')
#        version = self.mm[8:16].rstrip(b'\0')
#        print(magic, version)

        self._read_header(self.mm)

        self.mm.close()
        self.filehandle.flush()

        # open the file in C
        self.c_filehandle = lib.fdopen(os.dup(self.filehandle.fileno()),
                                       "rb".encode('ascii'))
        self.reset_rec_num()

    def __enter__(self):
        """Enter with section."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Exit with section."""
        self.close()

    def close(self):
        """Close the file both in C and in Python."""
        try:
            lib.fclose(self.c_filehandle)
        except:
            pass
        finally:
            self.filehandle.close()

    def _ptu_read_tag(self, s, offset):
        """
        Private function which reads a tag from the file header.

        Args:
            s (str): String read from the file (can be a memory map object)
            offset (int): Bytes offset for the reader in the string s

        Returns:
            3-tuple: tag identifier string, tag content and new offset for the
            string reader.
        """
        # Get the header struct as a tuple
        # Struct fields: 32-char string, int32, uint32, int64
        tag_struct = struct.unpack('32s i I q', s[offset:offset + 48])
        offset += 48
        # and save it into a dict
        tagname = tag_struct[0].rstrip(b'\0').decode()
        keys = ('idx', 'type', 'value')
        tag = {k: v for k, v in zip(keys, tag_struct[1:])}
        # Recover the name of the type (a string)
        tag['type'] = self.tag_type_r[tag['type']]

        # Some tag types need conversion
        if tag['type'] == 'tyFloat8':
            tag['value'] = np.int64(tag['value']).view('float64')
        elif tag['type'] == 'tyBool8':
            tag['value'] = bool(tag['value'])
        elif tag['type'] == 'tyTDateTime':
            TDateTime = np.uint64(tag['value']).view('float64')
            t = time.gmtime(self._ptu_TDateTime_to_time_t(TDateTime))
            tag['value'] = time.strftime("%Y-%m-%d %H:%M:%S", t)

        # Some tag types have additional data
        if tag['type'] == 'tyAnsiString':
            tag['data'] = s[offset: offset + tag['value']].rstrip(b'\0').decode()
            offset += tag['value']
        elif tag['type'] == 'tyFloat8Array':
            tag['data'] = np.frombuffer(s, dtype='float', count=tag['value']/8)
            offset += tag['value']
        elif tag['type'] == 'tyWideString':
            # WideString use type WCHAR in the original C++ demo code.
            # WCHAR size is not fixed by C++ standard, but on windows
            # is 2 bytes and the default encoding is UTF-16.
            # I'm assuming this is what the PTU requires.
            tag['data'] = s[offset: offset + tag['value'] * 2].decode('utf16')
            offset += tag['value']
        elif tag['type'] == 'tyBinaryBlob':
            tag['data'] = s[offset: offset + tag['value']]
            offset += tag['value']

        return tagname, tag, offset

    def _ptu_TDateTime_to_time_t(self, TDateTime):
        """
        Private function which converts a time in python format.

        (days from 01/01/1900) to timestamp in seconds since 01/01/1970

        Args:
            TDateTime (int): number of days since 01/01/1970

        Returns:
            int: timestamp corresponding to the date, number of seconds since
                 01/01/1970.
        """
        EpochDiff = 25569  # days between 30/12/1899 and 01/01/1970
        SecsInDay = 86400  # number of seconds in a day
        return (TDateTime - EpochDiff) * SecsInDay

    def _read_header(self, s):
        """
        Private function which reads the file header given a string read from
        the file (can be a memory-map object).
        Will store important tags in object properties.

        Args:
            s (str): string read from the file (can be a memory-map object).
        """
        offset = 16
        tag_end_offset = s.find(self.FileTagEnd.encode())

        tagname, tag, offset = self._ptu_read_tag(s, offset)
        self.tags[tagname] = tag
        while offset < tag_end_offset:
            tagname, tag, offset = self._ptu_read_tag(s, offset)
            self.tags[tagname] = tag

        # some mandatory tags and vars for the data processing
        self.record_type = self.rec_type_r[self.tags[self.TTTRTagTTTRRecType]['value']]
        self.num_records = self.tags[self.TTTRTagNumRecords]['value']
        self.globres = self.tags[self.TTTRTagGlobRes]['value']
        self.acq_time = int(self.tags[self.TTTRTagAcquisitionTime]['value'] * 1e9)  # in picoseconds
        self.end_header_offset = offset + 48  # + 48 for the end-of-header tag

    def print_header(self):
        """Print the header tags identifiers and values."""
        line = '{:30s} %s {:8}  {:12} '
        for key, tag in self.tags.items():
            value_fmt = '{:>20}'
            if tag['type'] == 'tyFloat8':
                value_fmt = '{:20.4g}'
            endline = '\n'
            if tag['type'] == 'tyAnsiString':
                endline = tag['data'] + '\n'
            print((line % value_fmt).format(key,
                                            tag['value'],
                                            tag['idx'],
                                            tag['type'],
                                            end=endline))

    def reset_rec_num(self, at=None):
        """
        Reset the variables and file readers to a given position in the file.

        By default we go to the beggining of the record section.

        Args:
            at (int, optional): offset in record numbers from the beginning of
            the records section.
        """
        self.c_rec_num = ffi.new("uint64_t *")
        if at is not None:
            self.c_rec_num[0] = at
        else:
            self.c_rec_num[0] = 0

        offset = self.c_rec_num[0] * 4 + self.end_header_offset
        lib.c_fseek(self.c_filehandle, offset)


class PTUmeasurement():

    """
    PTUmeasurement() analyses a PTUfile object to extract meaningful data like
    timetraces or g2 measurements.

    Attributes:
        meas (PTUfile object): the PTUfile object, corresponding to a
                               measurement file to analyse.
    """

    def __init__(self, ptu_file):
        """
        Constructs a PTUmeasurement object which allows for analysis of a PTUfile object.

        Args:
            ptu_file (PTUfile object): the PTUfile object to analyse.
                                       Should be open when analysing.
        """
        self.meas = ptu_file

    def _select_record_library(self):
        record_type = self.meas.rec_type[self.meas.record_type]
        # rtPicoHarpT3     = 0x00010303  rtPicoHarpT2     = 0x00010203
        # rtHydraHarpT3    = 0x00010304  rtHydraHarpT2    = 0x00010204
        # rtHydraHarp2T3   = 0x01010304  rtHydraHarp2T2   = 0x01010204
        # rtTimeHarp260NT3 = 0x00010305  rtTimeHarp260NT2 = 0x00010205
        # rtTimeHarp260PT3 = 0x00010306  rtTimeHarp260PT2 = 0x00010206
        if record_type in [0x01010204, 0x00010205, 0x00010206]:
            rec_lib = HHT2_HH2_lib
        elif record_type == 0x00010204:
            rec_lib = HHT2_HH1_lib
        elif record_type == 0x00010203:
            rec_lib = PHT2_lib
        else:
            print('Not implemented record type!')
            assert(0)

        return rec_lib

    def timetrace(self, resolution=1, record_range=[0, None], time_range=None, channel=-1, n_threads=2):
        """ 
        resolution is the time bin of the timetrace in seconds,
        record_range[0] is the record number of the position where to start the timetrace, default is 0
        record_range[1] is the record number of the position where to stop the timetrace, default is None (until the end of the measurement)
        time_range is the time length between stop_record and start_record, default is None (will take the full measurement time)
        """
        tt_f = self._select_record_library().timetrace

        # resolution in "timetag unit" i.e. number of globres
        # globres is in seconds
        resolution = int(float(resolution) / self.meas.globres)

        if record_range is None:
            record_range = [0, None]

        if record_range[1] is None:
            record_range[1] = self.meas.num_records

        record_range = [int(record_range[0]), int(record_range[1])]

        assert(0 <= record_range[0] < record_range[1] <= self.meas.num_records)

        if time_range is None:
            time_range = self.meas.acq_time * 1.e-12  # acq_time is in picoseconds

        self.meas.reset_rec_num(at=record_range[0])

        nb_of_bins = int(float(time_range) / self.meas.globres / resolution)
        c_time_trace = ffi.new("int[{}]".format(nb_of_bins))
        c_rec_num_trace = ffi.new("uint64_t[{}]".format(nb_of_bins))

        filepath = ffi.new("char[]", self.meas.filename.encode('ascii'))
        tt_f(filepath,
             self.meas.end_header_offset,
             record_range[0],
             record_range[1] - record_range[0],
             resolution,
             c_time_trace,
             c_rec_num_trace,
             channel,
             nb_of_bins,
             n_threads)

        time_vector = np.arange(nb_of_bins, dtype='float') * resolution * self.meas.globres
        time_trace = np.array([element for element in c_time_trace])
        rec_num_trace = np.array([element for element in c_rec_num_trace])

        return (time_vector, time_trace, rec_num_trace)

    def calculate_g2(self, correlation_window, resolution=None,
                     post_selec_ranges=None, channel_start=0, channel_stop=1,
                     mode='ring', n_threads=2):
        """
        Return the g2 calculated from the file, given the start and stop
        channels and the record number range to analyse.

        Args:
            correlation_window (int): correlation window length in seconds
            resolution (int): length of one time bin in seconds
            post_selec_ranges (list, optional): 2 levels list (eg [[0,100]]).
                Each element of the first level is a 2-element list with a
                start record number and a stop record number.
                By default, will take all the measurement
            channel_start (int, optional): channel number of the start photons
                (default 0, sync)
            channel_stop (int, optional): channel number of the stop photons
                (default 1)
            mode (fast, ring, classics, symmetric): fast is the naive algorithm
                                          classic uses a linked list
                                          ring uses a ring buffer
                                          symmetric uses ring buffer and
                                          also looks at negative deltas.

        Returns:
            2-tuple: numpy array vector of times (beginning of each time-bin),
                     numpy array vector of histogram (number of start-stop
                        photon couples per delay time bin)
        """
        g2_f = self._select_record_library().calculate_g2

        # Resolution in globres units, typically picoseconds
        if mode == 'symmetric':
            correlation_window /= 2.

        correlation_window = int(float(correlation_window) / self.meas.globres)
        if resolution is None:
            nb_of_bins = 1024
            resolution = int(correlation_window / float(nb_of_bins))
        else:
            resolution = int(float(resolution) / self.meas.globres)

        nb_of_bins = int(np.floor(float(correlation_window) / resolution))
        if mode == 'symmetric':
            nb_of_bins *= 2

        if mode == 'fast':
            mode_idx = 0
        elif mode == 'ring':
            mode_idx = 1
        elif mode == 'classic':
            mode_idx = 2
        elif mode == 'symmetric':
            mode_idx = 3
        else:
            print('Non-Existent Mode!')
            assert(0)

        # Make sure the correlation window is a multiple of the resolution
        correlation_window = nb_of_bins * resolution

        histogram = np.zeros(nb_of_bins)

        # Prepare the post selection ranges
        if post_selec_ranges is None:
            records_per_thread = self.meas.num_records / n_threads
            post_selec_ranges = [[int(i*records_per_thread), int((i+1)*records_per_thread)]
                                 for i in range(n_threads)]

        for post_selec_range in post_selec_ranges:
            if post_selec_range[1] > self.meas.num_records:
                post_selec_range[1] = self.meas.num_records

        # Initialize output arrays
        c_histogram = ffi.new("int[{}]".format(nb_of_bins))
        for i in range(nb_of_bins):
            c_histogram[i] = 0

        # Initialize ranges C arrays
        nb_of_ranges = len(post_selec_ranges)
        c_post_select_starts = ffi.new("uint64_t[{}]".format(nb_of_ranges))
        c_post_select_stops = ffi.new("uint64_t[{}]".format(nb_of_ranges))
        for i in range(nb_of_ranges):
            c_post_select_starts[i] = post_selec_ranges[i][0]
            c_post_select_stops[i] = post_selec_ranges[i][1]

        filepath = ffi.new("char[]", self.meas.filename.encode('ascii'))
        buffer_size = 2**4

        # Calculate
        g2_f(filepath,                     # file to analyze
             self.meas.end_header_offset,  # header offset
             c_post_select_starts,       # array with starting record of ranges
             c_post_select_stops,        # array with stopping record of ranges
             nb_of_ranges,                 # number of postselection ranges
             correlation_window,           # length of g2
             c_histogram,                  # C vector with the histogram output
             nb_of_bins,                   # number of bins in histogram
             channel_start,                # start channel for g2 algorithm
             channel_stop,                 # stop channel for g2 algorithm
             buffer_size,                  # size of ring buffer for g2_ring
             n_threads,                    # number of parallelization threads
             mode_idx)                     # select g2 algorithm

        for i in range(nb_of_bins):
            histogram[i] = c_histogram[i]

        if mode == 'symmetric':
            time_vector = np.linspace(-nb_of_bins/2., nb_of_bins/2., nb_of_bins, dtype='float') * resolution * self.meas.globres
        else:
            time_vector = np.arange(nb_of_bins, dtype='float') * resolution * self.meas.globres
        return (time_vector, np.array(histogram))


def construct_postselect_vector(timetrace_y, timetrace_recnum, threshold, above=True):
    """Constructs a postselection vector based on a threshold condition, selecting items above or below as specified by user.
    
    Args:
        timetrace_y (numpy array): number of photons in the time bin.
        timetrace_recnum (numpy array): array of record numbers corresponding to timebin start (for each time bin of the time trace)
        threshold (number): number of photons threshold used for selection.
        above (bool, optional): if True, will return time ranges where timetrace_y > threshold. If False, will return time ranges where timetrace_y < threshold.
    
    Returns:
        (list, list): first: List of 2-item lists being [start_index, stop_index] array indices in timetrace_y for each post-selected range of points.
                      second: corresponding record numbers allowing for direct use as post-selection array in the calculate_g2() function.
    """

    if above:
        select = timetrace_y > threshold
        if timetrace_y[0] > threshold:
            add_first_time = True
        else:
            add_first_time = False
    else:
        select = timetrace_y < threshold
        if timetrace_y[0] < threshold:
            add_first_time = True
        else:
            add_first_time = False

    nselect = ~select
    post_selec_ranges = ((select[:-1] & nselect[1:]) | (nselect[:-1] & select[1:])).nonzero()[0]
    if add_first_time:
        post_selec_ranges = np.insert(post_selec_ranges, 0, 0).astype(int)
    if len(post_selec_ranges) % 2 != 0:  # odd length means we need to add the end time
        post_selec_ranges = np.append(post_selec_ranges, len(timetrace_y) - 2).astype(int)  # - 2 because 1 will be added below

    # # take into account the end of the interval (otherwise stops 1 time bin before)
    # for i, value in enumerate(post_selec_ranges[1::2]):
    #     post_selec_ranges[2*i - 1] += 1

    # # the previous step will introduce [0,11], [11,15] kind of situation (where ranges are contiguous). We remove these intermediate duplicates.
    # double_index = (post_selec_ranges[:-1] != post_selec_ranges[1:])
    # double_index = (np.append(double_index, True) & np.insert(double_index, 0, True))
    # post_selec_ranges = post_selec_ranges[double_index]

    post_selec_ranges = post_selec_ranges.reshape((-1, 2))

    recnum_post_selec_ranges = [[timetrace_recnum[post_selec_range[0]], timetrace_recnum[post_selec_range[1]]] for post_selec_range in post_selec_ranges]

    return post_selec_ranges, recnum_post_selec_ranges


if __name__ == '__main__':

    timetrace_resolution = 10    # in seconds
    g2_resolution = 600 * 1e-12  # picoseconds * 1e-12 to change to seconds
    g2_window = 50000 * 1e-12   # picoseconds * 1e-12 to change to seconds
    threshold = 0 # in unit of counts per unit timetrace resolution
    constraint_above = True # post-selection counts > threshold for the given timetrace resolution

    filename = r'/Users/garfield/Downloads/test_big.ptu'

    with PTUfile(filename) as ptu_file:
        ptu_file.print_header()
        ptu_meas = PTUmeasurement(ptu_file)

        start_time = time.time()
        timetrace_x, timetrace_y, timetrace_recnum =\
            ptu_meas.timetrace(resolution=timetrace_resolution, n_threads=4)
        stop_time = time.time()
        read_speed = os.path.getsize(filename)/float(stop_time - start_time)/1024./1024./1024.
        print('timetrace calculation took', stop_time - start_time, 's')
        print('processing speed:', read_speed, 'GBps')
        print('Total number of photons: ', np.sum(timetrace_y))

        plt.figure()
        plt.plot(timetrace_x, timetrace_y)
        plt.xlabel('Time (s)')
        plt.ylabel('Counts/{} s'.format(timetrace_resolution))
        plt.title('Timetrace')

        plt.figure()
        plt.plot(timetrace_x, timetrace_recnum)
        plt.xlabel('Time (s)')
        plt.ylabel('Record number')
        plt.title('Record number vs measurement time')

        # to post-selection photon counts  > (or <) threshold
        post_selec_ranges, recnum_post_selec_ranges = construct_postselect_vector(timetrace_y,timetrace_recnum,threshold,constraint_above)

        # to post-select counts within a certain time (say from time 0 to 10 s)
        # recnum_post_selec_ranges = [[0, timetrace_recnum[find_nearest(timetrace_x,10)[0]]]]


        if (len(recnum_post_selec_ranges))==1:
            recnum_post_selec_ranges = None

        start_time = time.time()
        hist_x_ring, hist_y_ring = ptu_meas.calculate_g2(g2_window, g2_resolution,
                                                         post_selec_ranges=recnum_post_selec_ranges,
                                                         n_threads=4,
                                                         mode='symmetric')
        stop_time = time.time()
        read_speed = os.path.getsize(filename)/float(stop_time - start_time)/1024./1024./1024.
        print('g2 calculation took', stop_time - start_time, 's')
        print('processing speed:', read_speed, 'GBps')
        plt.figure()
        plt.plot(hist_x_ring * 1e9, hist_y_ring)
        plt.xlabel('Delay (ns)')
        plt.title('G2 measurements')

        # start_time = time.time()
        # print('\nclassic ALGORITHM')
        # hist_x_classic, hist_y_classic = ptu_meas.calculate_g2(g2_window, g2_resolution,
        #                                        post_selec_ranges=None,
        #                                        n_threads=1,
        #                                        mode='classic')
        # stop_time = time.time()
        # print('g2 calculation took', stop_time - start_time, 's')
        # plt.figure()
        # plt.plot(hist_x_classic * 1e9, hist_y_classic)
        # plt.xlabel('Delay (ns)')
        # plt.title('G2 measurements (Classic algorithm)')

        # plt.figure()
        # plt.plot(hist_x_classic * 1e9, hist_y_classic-hist_y_ring)
        # plt.xlabel('Delay (ns)')
        # plt.title('Difference')

        plt.show()
