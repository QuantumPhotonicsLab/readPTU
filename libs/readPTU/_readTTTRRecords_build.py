#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 16:16:03 2017

@author: raphaelproux
"""

from cffi import FFI
from os import remove
import setuptools  # necessary magic import for windows -- don't think, just accept it!

# To add a new parsers just add the appropiate parser in parsers.c with the
# name convention Parse____ where the gap corresponds to the record
# type format name. Then add the same name to the list below and you are done!

# You will also have to add the new parser in readPTU.py

# Finally as it is now the records are assumed to be 32 bits long. Using
# the same conidtional compiling strategy as below this could be easily
# modifiable for other types of records.

parsers = ['HHT2_HH1', 'HHT2_HH2', 'PHT2']

with open('readTTTRRecords-for-import.c', 'r') as myfile:
    c_code = myfile.read()
    codes_dict = {k: c_code.replace("##parser##", k) for k in parsers}

prototypes = r"""
    void timetrace(char filepath[], int end_of_header, uint64_t *RecNum,
                   uint64_t NumRecords, uint64_t time_bin_length,
                   int time_trace[], uint64_t RecNum_trace[], int nb_of_bins,
                   int n_threads);
    void calculate_g2(char filepath[], int end_of_header,
                      uint64_t *RecNum_start, uint64_t *RecNum_stop,
                      int nb_of_ranges, uint64_t max_time, int histogram[],
                      int nb_of_bins, int channel_start, int channel_stop,
                      int buffer_size, int n_threads, int mode);
    FILE *fdopen(int, const char *);   // from the C <stdio.h>
    int fclose(FILE *);
    int c_fseek(FILE *, long int);"""

if __name__ == "__main__":
    for code in codes_dict:
        ffibuilder = FFI()
        ffibuilder.cdef(prototypes)
        ffibuilder.set_source("_readTTTRRecords_{}".format(code),
                              codes_dict[code])
        print('\nCompiling version {}'.format(code))
        ffibuilder.compile(verbose=True)
        remove('_readTTTRRecords_{}.c'.format(code))
        remove('_readTTTRRecords_{}.o'.format(code))
