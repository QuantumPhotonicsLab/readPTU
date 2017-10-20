#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  4 18:36:54 2017

@author: raphaelproux

Largely inspired from PicoQuant examples:
    https://github.com/PicoQuant/PicoQuant-Time-Tagged-File-Format-Demos
and (in particular, to read the header) from a jupyter notebook
by tritemio on GitHub:
    https://gist.github.com/tritemio/734347586bc999f39f9ffe0ac5ba0e66

Please note:
- this library has been tested only with T2 mode of a Hydraharp v2, but it should work with a Hydraharp v1 or a Picoharp file (written, untested),
- does not support T3 mode for now (not written), but should be easy to implement

"""

import pylab as pl
import os
import sys
import struct
import mmap
import time
import collections as coll
from _readTTTRRecords import ffi, lib


class PTUfile():
    """
    PTUfile() handles Picoquant PTU files opening and closing, and header 
    reading to extract information like the record type, measurement time and
    number of records. More detailed analysis like timetraces and g2 are done
    using the PTUmeasurement class.
    
    Attributes:
        acq_time (int): measurement acquisition time in picoseconds
        c_filehandle (C FILE *): C file handle used for possible analysis in C
        c_rec_num (C uint64_t *): C pointer to record number (used when analysing the file in C)
        end_header_offset (int): Position in bytes of the beginning of the records section in the file
        filehandle (TYPE): Python filehandle used only to read the header
        FileTagEnd (str): Constant string, tag to detect end of header
        globres (float): resolution of timetags in seconds
        mm (mmap object): memory map object created to read the file header
        num_records (int): number of records in the file
        rec_type (dict): reference dictionary of record types identifiers
        rec_type_r (dict): reversed dictionary of rec_type (allows inverse referencing)
        record_type (int): string identifier of the record type for the current file. Note you should send rec_type[record_type] to use the int identifier (ie in C).
        tag_type (dict): reference dictionary of header tag types identifiers
        tag_type_r (dict): reversed dictionary of tag_type
        tags (ordered dict): ordered dictionary of the header tags for the current file
        TTTRTagAcquisitionTime (str): header tag string identifier for acquisition time
        TTTRTagGlobRes (str): header tag string identifier for global resolution
        TTTRTagNumRecords (str): header tag string identifier for number of records
        TTTRTagRes (str): header tag string identifier for delay time resolution (T3 mode)
        TTTRTagTTTRRecType (str): header tag string identifier for record type
    
    """

    # Constants
    TTTRTagTTTRRecType = "TTResultFormat_TTTRRecType"
    TTTRTagNumRecords  = "TTResult_NumberOfRecords"  # Number of TTTR Records in the File;
    TTTRTagRes         = "MeasDesc_Resolution"       # Resolution for the Dtime (T3 Only)
    TTTRTagGlobRes     = "MeasDesc_GlobalResolution" # Global Resolution of TimeTag(T2) /NSync (T3)
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
        tyBinaryBlob  = 0xFFFFFFFF,
        )
    
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
        Constructs a PTUfile object given a filename and number of measurement channels. 
        Will open the file and leave it opened for possible analysis with a PTUmeasurement object.
        Please note PTUfile supports context manager and it is strongly advised to create it using
        a "with ... as ..." statement.
        
        Args:
            filename (str): path + filename to the measurement file to open
        """
        # Reverse mappings of the tags and record dictionaries
        self.tag_type_r = {v: k for k, v in self.tag_type.items()}
        self.rec_type_r = {v: k for k, v in self.rec_type.items()}
        
        # var initialization
        self.tags = coll.OrderedDict()

#        print('Size of file:', os.path.getsize(filename))

        self.filehandle = open(filename, 'rb')
        self.mm = mmap.mmap(self.filehandle.fileno(), 0, access=mmap.ACCESS_READ)  # with mmap, we will be able to use the file as a string

#        magic = self.mm[:8].rstrip(b'\0')
#        version = self.mm[8:16].rstrip(b'\0')
#        print(magic, version)

        self._read_header(self.mm)

        self.mm.close()
        self.filehandle.flush()

        # open the file in C
        self.c_filehandle = lib.fdopen(os.dup(self.filehandle.fileno()), "rb".encode('ascii'))
        self.reset_rec_num()


    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def close(self):
        """
        Closes the file both in C and in Python.
        """
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
            3-tuple: tag identifier string, tag content and new offset for the string reader.
        """
        # Get the header struct as a tuple
        # Struct fields: 32-char string, int32, uint32, int64
        tag_struct = struct.unpack('32s i I q', s[offset:offset+48])
        offset += 48
        # and save it into a dict
        tagname = tag_struct[0].rstrip(b'\0').decode()
        keys = ('idx', 'type', 'value')
        tag = {k: v for k, v in zip(keys, tag_struct[1:])}
        # Recover the name of the type (a string)
        tag['type'] = self.tag_type_r[tag['type']]
        
        # Some tag types need conversion
        if tag['type'] == 'tyFloat8':
            tag['value'] = pl.int64(tag['value']).view('float64')
        elif tag['type'] == 'tyBool8':
            tag['value'] = bool(tag['value'])
        elif tag['type'] == 'tyTDateTime':
            TDateTime = pl.uint64(tag['value']).view('float64')
            t = time.gmtime(self._ptu_TDateTime_to_time_t(TDateTime))
            tag['value'] = time.strftime("%Y-%m-%d %H:%M:%S", t)
            
        # Some tag types have additional data
        if tag['type'] == 'tyAnsiString':
            tag['data'] = s[offset: offset + tag['value']].rstrip(b'\0').decode()
            offset += tag['value']
        elif tag['type'] == 'tyFloat8Array':
            tag['data'] = pl.frombuffer(s, dtype='float', count=tag['value']/8)
            offset += tag['value']
        elif tag['type'] == 'tyWideString':
            # WideString use type WCHAR in the original C++ demo code.
            # WCHAR size is not fixed by C++ standard, but on windows 
            # is 2 bytes and the default encoding is UTF-16.
            # I'm assuming this is what the PTU requires.
            tag['data'] = s[offset: offset + tag['value']*2].decode('utf16')
            offset += tag['value']
        elif tag['type'] == 'tyBinaryBlob':
            tag['data'] = s[offset: offset + tag['value']]
            offset += tag['value']
            
        return tagname, tag, offset

    def _ptu_TDateTime_to_time_t(self, TDateTime):
        """
        Private function which converts a time in python format (days from 01/01/1900) to timestamp in seconds since 01/01/1970
        
        Args:
            TDateTime (int): number of days since 01/01/1970
        
        Returns:
            int: timestamp corresponding to the date, number of seconds since 01/01/1970.
        """
        EpochDiff = 25569  # days between 30/12/1899 and 01/01/1970
        SecsInDay = 86400  # number of seconds in a day
        return (TDateTime - EpochDiff) * SecsInDay
        
    def _read_header(self, s):
        """
        Private function which reads the file header given a string read from the file (can be a memory-map object).
        Will store important tags in object properties.
        
        Args:
            s (str): string read from the file (can be a memory-map object).
        """
        offset = 16
        tag_end_offset = s.find(self.FileTagEnd.encode())
        
        tagname, tag, offset = self._ptu_read_tag(s, offset)
        self.tags[tagname] = tag
        while offset < tag_end_offset:
            tagname, tag, offset =self._ptu_read_tag(s, offset)
            self.tags[tagname] = tag
            
        # some mandatory tags and vars for the data processing
        self.record_type = self.rec_type_r[self.tags[self.TTTRTagTTTRRecType]['value']]
        self.num_records = self.tags[self.TTTRTagNumRecords]['value']
        self.globres = self.tags[self.TTTRTagGlobRes]['value']
        self.acq_time = int(self.tags[self.TTTRTagAcquisitionTime]['value'] * 1e9)  # in picoseconds
        self.end_header_offset = offset + 48  # + 48 for the end-of-header tag
            
    def print_header(self):
        """
        Prints the header tags identifiers and values.
        """
        line = '{:30s} %s {:8}  {:12} '
        for key, tag in self.tags.items():
            value_fmt = '{:>20}'
            if tag['type'] == 'tyFloat8':
                value_fmt = '{:20.4g}'
            endline = '\n'
            if tag['type'] == 'tyAnsiString':
                endline = tag['data'] + '\n'
            print((line % value_fmt).format(key, tag['value'], tag['idx'], tag['type']), 
                  end=endline)
            
    def reset_rec_num(self, at=None):
        """
        Resets the variables and file readers to a given position in the file (by default, the beginning of the records section).
        
        Args:
            at (int, optional): offset in record numbers from the beginning of the records section.
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
    PTUmeasurement() analyses a PTUfile object to extract meaningful data like timetraces or g2 measurements.

    Attributes:
        meas (PTUfile object): the PTUfile object, corresponding to a measurement file to analyse.
    """
    
    def __init__(self, ptu_file):
        """
        Constructs a PTUmeasurement object which allows for analysis of a PTUfile object.

        Args:
            ptu_file (PTUfile object): the PTUfile object to analyse. Should be open when analysing.
        """
        self.meas = ptu_file

    def time_trace(self, resolution=1):
        """
        Returns the timetrace, a tuple where first item is a time vector 
        and second item is the vector of number of photons per resolution bin.

        Parameters:
            resolution: the time resolution in seconds

        Returns:
            3-tuple: numpy array vector of times, vector of intensity per bin, vector of record numbers locating the end of each bin in the file.
        """
        resolution *= 1e12 # now in picoseconds
        self.meas.reset_rec_num()
        
        nb_of_bins = int(self.meas.acq_time / resolution) ;
        c_time_vector = ffi.new("uint64_t[{}]".format(nb_of_bins))
        c_time_trace = ffi.new("int[{}]".format(nb_of_bins))
        c_rec_num_trace = ffi.new("uint64_t[{}]".format(nb_of_bins))

        lib.timetrace(self.meas.c_filehandle,
                      self.meas.rec_type[self.meas.record_type],
                      self.meas.end_header_offset,
                      self.meas.c_rec_num,
                      self.meas.num_records,
                      int(resolution),
                      c_time_vector,
                      c_time_trace,
                      c_rec_num_trace,
                      nb_of_bins)

        time_vector = [element for element in c_time_vector]
        time_trace = [element for element in c_time_trace]
        rec_num_trace = [element for element in c_rec_num_trace]

        return pl.array(time_vector), pl.array(time_trace), pl.array(rec_num_trace)

    def calculate_g2(self, correlation_window, resolution, post_selec_ranges=None, channel_start=0, channel_stop=1, fast=True):
        """
        Returns the g2 calculated from the file, given the start and stop channels and the record number range to analyse (which allows post-selection)

        Args:
            correlation_window (int): correlation window length in number of global resolutions (typically picoseconds)
            resolution (int): length of one time bin in number of global resolutions (typically picoseconds)
            post_selec_ranges (list, optional): 2 levels list (eg [[0,100]]). Each element of the first level is a 2-element list 
                with a start record number and a stop record number. By default, will take all the measurement (post_selec_ranges=None)
            channel_start (int, optional): channel number of the start photons (default 0, sync)
            channel_stop (int, optional): channel number of the stop photons (default 1)
            fast (bool, optional): in fast mode (default, fast=True), the g2 is calculated using subsequent pairs of start-stop photons reading 
                them chronologically along the file (start-stop -> start-stop -> etc.) and discarding other photons. This algorithm will produce an 
                exponential decay artefact on long time scales or with high photon rates. In not fast mode (fast=False), the g2 is calculated 
                considering all start-stop photon combinations, therefore not exhibiting any artefact.

        Returns:
            2-tuple: numpy array vector of times (beginning of each time-bin), 
                     numpy array vector of histogram (number of start-stop photon couples per delay time bin)
        """
        if fast:
            calc_g2 = lib.calculate_g2_fast
        else:
            calc_g2 = lib.calculate_g2

        nb_of_bins = int(pl.floor(float(correlation_window) / resolution))
        histogram = pl.zeros(nb_of_bins)

        if post_selec_ranges is None:
            post_selec_ranges = [[0, self.meas.num_records]]

        for post_selec_range in post_selec_ranges:
            if post_selec_range[1] > self.meas.num_records:
                post_selec_range[1] = self.meas.num_records

            print(post_selec_range)

            c_time_vector = ffi.new("uint64_t[{}]".format(nb_of_bins+1))
            for i in range(nb_of_bins+1):
                c_time_vector[i] = i * resolution

            c_histogram = ffi.new("int[{}]".format(nb_of_bins))
            for i in range(nb_of_bins):
                c_histogram[i] = 0

#            print(nb_of_bins)
            calc_g2(self.meas.c_filehandle,
                    self.meas.rec_type[self.meas.record_type],
                    self.meas.end_header_offset,
                    self.meas.c_rec_num,
                    self.meas.num_records,
                    post_selec_range[0],
                    post_selec_range[1],
                    c_time_vector,
                    c_histogram,
                    nb_of_bins,
                    channel_start,
                    channel_stop)

            for i in range(nb_of_bins):
                histogram[i] += c_histogram[i]

        time_vector = [time for time in c_time_vector]

        return pl.array(time_vector[:-1]), pl.array(histogram)

    def calculate_g2_ring(self, correlation_window, resolution,
                          post_selec_ranges=None, channel_start=0,
                          channel_stop=1, buffer_size=2**10):
        """
        Return g2 using ring buffer algorithm.

        Calculated from the file, given the start and stop channels and the
        record number range to analyse (which allows post-selection)

        Args:
            correlation_window (int): correlation window length in number of
                global resolutions (typically picoseconds)
            resolution (int): length of one time bin in number of global
                resolutions (typically picoseconds)
            post_selec_ranges (list, optional): 2 levels list (eg [[0,100]]).
                Each element of the first level is a 2-element list with a
                start record number and a stop record number. By default,
                will take all the measurement (post_selec_ranges=None)
            channel_start (int, optional): channel number of the start photons (default 0, sync)
            channel_stop (int, optional): channel number of the stop photons (default 1)
            fast (bool, optional): in fast mode (default, fast=True), the g2 is
                calculated using subsequent pairs of start-stop photons reading
                them chronologically along the file (start-stop -> start-stop -> etc.) and discarding other photons. This algorithm will produce an 
                exponential decay artefact on long time scales or with high photon rates. In not fast mode (fast=False), the g2 is calculated 
                considering all start-stop photon combinations, therefore not exhibiting any artefact.

        Returns:
            2-tuple: numpy array vector of times (beginning of each time-bin),
                     numpy array vector of histogram (number of start-stop photon couples per delay time bin)
        """
        calc_g2 = lib.calculate_g2_ring

        nb_of_bins = int(pl.floor(float(correlation_window) / resolution))
        histogram = pl.zeros(nb_of_bins)

        if post_selec_ranges is None:
            post_selec_ranges = [[0, self.meas.num_records]]

        for post_selec_range in post_selec_ranges:
            if post_selec_range[1] > self.meas.num_records:
                post_selec_range[1] = self.meas.num_records

            print(post_selec_range)

            c_time_vector = ffi.new("uint64_t[{}]".format(nb_of_bins + 1))
            for i in range(nb_of_bins + 1):
                c_time_vector[i] = i * resolution

            c_histogram = ffi.new("int[{}]".format(nb_of_bins))
            for i in range(nb_of_bins):
                c_histogram[i] = 0

#            print(nb_of_bins)
            calc_g2(self.meas.c_filehandle,
                    self.meas.rec_type[self.meas.record_type],
                    self.meas.end_header_offset,
                    self.meas.c_rec_num,
                    self.meas.num_records,
                    post_selec_range[0],
                    post_selec_range[1],
                    c_time_vector,
                    c_histogram,
                    nb_of_bins,
                    channel_start,
                    channel_stop,
                    buffer_size)

            for i in range(nb_of_bins):
                histogram[i] += c_histogram[i]

        time_vector = [time for time in c_time_vector]

        return pl.array(time_vector[:-1]), pl.array(histogram)


if __name__ == '__main__':

    timetrace_resolution = 1  # in seconds

    g2_resolution = 1000  # in picoseconds
    g2_coincidence_window = 1e6  # in picoseconds
    g2_post_selec_ranges = [[4e5,1e6]]  # in record numbers

    filename = r'C:\Users\QPL\Desktop\temp_measurements\default.ptu'

    # filename = r'/Users/raphaelproux/Desktop/TTTR/t2htr2a1loc2.ptu'

    with PTUfile(filename) as ptu_file:
        ptu_file.print_header()
        ptu_meas = PTUmeasurement(ptu_file)

        start_time = time.time()
        timetrace_x, timetrace_y, timetrace_recnum = ptu_meas.time_trace(resolution=timetrace_resolution)
        stop_time = time.time()
        print('timetrace calculation took', stop_time - start_time, 's')

        pl.savetxt('timetrace.txt', pl.array([timetrace_x, timetrace_y, timetrace_recnum]).transpose(), delimiter='\t')

        start_time = time.time()
        histo_x, histo_y = ptu_meas.calculate_g2(1000000,10000, post_selec_ranges=[[0,100000000]])  #, fast=False)
        stop_time = time.time()
        print('g2 calculation took', stop_time - start_time, 's')

        pl.savetxt('timetrace.txt', pl.array([timetrace_x, timetrace_y, timetrace_recnum]).transpose(), delimiter='\t')

        pl.figure()
        pl.plot(timetrace_x, timetrace_y)
        pl.xlabel('Time (ps)')
        pl.ylabel('Counts/{} s'.format(timetrace_resolution))
        pl.title('Timetrace')
        pl.figure()
        pl.plot(timetrace_x, timetrace_recnum)
        pl.xlabel('Time (ps)')
        pl.ylabel('Record number')
        pl.title('Record number vs measurement time')
        pl.figure()
        pl.plot(histo_x, histo_y)
        pl.xlabel('Delay (ps)')
        pl.ylabel('Coincidence counts/{} ps'.format(g2_resolution))
        pl.title('g2 measurement')
        pl.show()
