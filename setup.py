from setuptools.command.install import install
from setuptools.command.develop import develop

from cffi import FFI
from os import remove
import setuptools  # necessary magic import for windows -- don't think, just accept it!

# To add a new parsers just add the appropiate parser in parsers.c with the
# name convention Parse____ where the gap corresponds to the record
# type format name. Then add the same name to the list below and you are done!

# You will also have to add the new parser in readPTU.py in the method
# _select_record_library of PTUmeasurement.

# Finally as it is now the records are assumed to be 32 bits long. Using
# the same conidtional compiling strategy as below this could be easily
# modifiable for other types of records.


def compile_library():
    parsers = ['HHT2_HH1', 'HHT2_HH2', 'PHT2']

    with open('./readPTU/readTTTRRecords-for-import.c', 'r') as myfile:
        c_code = myfile.read()
        codes_dict = {k: c_code.replace("##parser##", k) for k in parsers}

    prototypes = r"""
        void timetrace(char filepath[], int end_of_header, uint64_t RecNum_start,
                       uint64_t NumRecords, uint64_t time_bin_length,
                       int time_trace[], uint64_t RecNum_trace[], int select_channel,
                       int nb_of_bins, int n_threads);
        void calculate_g2(char filepath[], int end_of_header,
                          uint64_t *RecNum_start, uint64_t *RecNum_stop,
                          int nb_of_ranges, uint64_t max_time, int histogram[],
                          int nb_of_bins, int channel_start, int channel_stop,
                          int buffer_size, int n_threads, int mode);
        FILE *fdopen(int, const char *);   // from the C <stdio.h>
        int fclose(FILE *);
        int c_fseek(FILE *, long int);"""

    for code in codes_dict:
        ffibuilder = FFI()
        ffibuilder.cdef(prototypes)
        ffibuilder.set_source("readPTU._readTTTRRecords_{}".format(code),
                              codes_dict[code],
                              extra_compile_args=["-O3"])
        print('\nCompiling version {}'.format(code))
        ffibuilder.compile(verbose=True)
        remove('./readPTU/_readTTTRRecords_{}.c'.format(code))
        # remove('./readPTU/_readTTTRRecords_{}.o'.format(code))


class Build(install):
    """Custom handler for the 'install' command."""

    def run(self):
        compile_library()
        install.run(self)

class Build_dev(develop):
    """Custom handler for the 'install' command."""

    def run(self):
        compile_library()
        develop.run(self)


setuptools.setup(
    name='readPTU',
    version='0.1',
    description='Read and Analyze PicoQuant file formats (ptu, pt2, pt3 ...)',
    url='https://github.com/QuantumPhotonicsLab/readPTU',
    author='QPL Lab',
    author_email='brian.gerardot@hw.ac.uk',
    license='MIT',
    packages=['readPTU'],
    zip_safe=False,
    setup_requires=["cffi>=1.11.2"],
    install_requires=["cffi>=1.11.2"],
    cmdclass={'install': Build, 'develop': Build_dev},
    package_data={'': ['*.so', '*.c', '*.o', '*.pyd']}) # Note that we also have to include .pyd files for binaries to be copied over on windows
