#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 16:16:03 2017

@author: raphaelproux
"""

from cffi import FFI
from os import remove
import setuptools  # necessary magic import for windows -- don't think, just accept it!

with open('readTTTRRecords-for-import.c', 'r') as myfile:
    c_code = myfile.read()
    c_code_HHT2_HH1 = c_code.replace("##parser##", "HHT2_HH1")
    c_code_HHT2_HH2 = c_code.replace("##parser##", "HHT2_HH2")
    c_code_PHT2 = c_code.replace("##parser##", "PHT2")


prototypes = r"""
    void timetrace(char filepath[], int end_of_header, uint64_t *RecNum, uint64_t NumRecords,
                   uint64_t time_bin_length, int time_trace[],
                   uint64_t RecNum_trace[], int nb_of_bins, int n_threads);
    void calculate_g2(char filepath[], int end_of_header, uint64_t *RecNum_start,
                      uint64_t *RecNum_stop, int nb_of_ranges, uint64_t max_time, int histogram[], int nb_of_bins,
                      int channel_start, int channel_stop, int buffer_size, int n_threads, int mode);
    FILE *fdopen(int, const char *);   // from the C <stdio.h>
    int fclose(FILE *);
    int c_fseek(FILE *, long int);"""

if __name__ == "__main__":
    ffibuilder = FFI()
    ffibuilder.cdef(prototypes)
    ffibuilder.set_source("_readTTTRRecords_HHT2_HH2", c_code_HHT2_HH2)
    print('\nCompiling version HHT2_HH2')
    ffibuilder.compile(verbose=True)
    remove('_readTTTRRecords_HHT2_HH2.c')
    remove('_readTTTRRecords_HHT2_HH2.o')

    ffibuilder = FFI()
    ffibuilder.cdef(prototypes)
    ffibuilder.set_source("_readTTTRRecords_HHT2_HH1", c_code_HHT2_HH1)
    print('\nCompiling version HHT2_HH1')
    ffibuilder.compile(verbose=True)
    remove("_readTTTRRecords_HHT2_HH1.c")
    remove("_readTTTRRecords_HHT2_HH1.o")

    ffibuilder = FFI()
    ffibuilder.cdef(prototypes)
    ffibuilder.set_source("_readTTTRRecords_PHT2", c_code_PHT2)
    print('\nCompiling version PHT2')
    ffibuilder.compile(verbose=True)
    remove("_readTTTRRecords_PHT2.c")
    remove("_readTTTRRecords_PHT2.o")
