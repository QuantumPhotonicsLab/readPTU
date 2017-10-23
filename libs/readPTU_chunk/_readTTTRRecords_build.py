#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 16:16:03 2017

@author: raphaelproux
"""

from cffi import FFI
import setuptools  # necessary magic import for windows -- don't think, just accept it!
ffibuilder = FFI()

with open('readTTTRRecords-for-import.c', 'r') as myfile:
    c_code=myfile.read()

# ffibuilder.cdef(r"""
#     void timetrace(FILE* filehandle, long long record_type, int end_of_header, uint64_t *RecNum, uint64_t NumRecords, uint64_t time_bin_length, uint64_t *time_vector, int *time_trace, uint64_t *RecNum_trace, int nb_of_bins);
#     void calculate_g2_ring(FILE* filehandle, long long record_type, int end_of_header, uint64_t *RecNum, uint64_t NumRecords, uint64_t RecNum_start, uint64_t RecNum_stop, uint64_t *time_vector, int *histogram, int nb_of_bins, int channel_start, int channel_stop, int buffer_size);
#     void calculate_g2(FILE* filehandle, long long record_type, int end_of_header, uint64_t *RecNum, uint64_t NumRecords, uint64_t RecNum_start, uint64_t RecNum_stop, uint64_t *time_vector, int *histogram, int nb_of_bins, int channel_start, int channel_stop);
#     void calculate_g2_fast(FILE* filehandle, long long record_type, int end_of_header, uint64_t *RecNum, uint64_t NumRecords, uint64_t RecNum_start, uint64_t RecNum_stop, uint64_t *time_vector, int *histogram, int nb_of_bins, int channel_start, int channel_stop);
#     FILE *fdopen(int, const char *);   // from the C <stdio.h>
#     int fclose(FILE *);
#     int c_fseek(FILE *, long int);""")

ffibuilder.cdef(r"""
    void timetrace(FILE* filehandle, long long record_type, int end_of_header, uint64_t *RecNum, uint64_t NumRecords, uint64_t time_bin_length, uint64_t *time_vector, int *time_trace, uint64_t *RecNum_trace, int nb_of_bins);
    void calculate_g2_ring(FILE* filehandle, long long record_type, int end_of_header, uint64_t *RecNum, uint64_t NumRecords, uint64_t RecNum_start, uint64_t RecNum_stop, uint64_t *time_vector, int *histogram, int nb_of_bins, int channel_start, int channel_stop, int buffer_size);
    FILE *fdopen(int, const char *);   // from the C <stdio.h>
    int fclose(FILE *);
    int c_fseek(FILE *, long int);""")
#void calculate_g2_fast(FILE* filehandle, long long record_type, int end_of_header, uint64_t *RecNum, uint64_t NumRecords, uint64_t RecNum_start, uint64_t RecNum_stop, uint64_t *time_vector, int *histogram, int nb_of_bins, int channel_start, int channel_stop);
ffibuilder.set_source("_readTTTRRecords", c_code)

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
