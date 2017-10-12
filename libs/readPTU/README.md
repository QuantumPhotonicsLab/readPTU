# readPTU library

PicoQuant uses a specific homemade file format called PTU for its time-tag time-resolved measurements. This format is subsequent to the former .pt2 or .pt3 and can handle both T2 and T3 acquisition modes. 

For now the library handles Hydraharp V2 T2 mode files (tested) but it should be able to handle Hydraharp V1 and Picoharp in T2 mode (untested but implemented).

## Use example
To start, you should open the PTU file using a PTUfile() object. Constructing a PTUfile() object automatically opens the file and leaves it open for further analysis. As a first step, you can then print the header which will bring useful information on the measurement, like measurement time, number of records, record format type, etc.
```
ptu_file = PTUfile(r'path/to/the/file')
ptu_file.print_header()
```

To analyse the file, e.g. to calculate a timetrace or a g2, you should use the PTUmeasurement() class. To construct a PTUmeasurement() object, you need to provide an open PTUfile object (```ptu_file``` in our example). For example, to calculate a timetrace with a time bin resolution of 1 second:
```
timetrace_x, timetrace_y, timetrace_recnum = ptu_meas.time_trace(resolution=1)
```

For a g2 measurement, you can use:
```
histo_x, histo_y = ptu_meas.calculate_g2(correlation_window=1000000, resolution=10000, fast=True)
```

The g2 measurement in particular has two available algorithms which may be selected using the ```fast``` parameter:

* ```fast=True``` is faster to compute (typical time: a few minutes for a 10 GB file). It consists in finding chronologically the sequence of start-stop — start-stop — start-stop... photons and record the associated delays in a histogram. This is very similar to the histogram mode of most correlation cards, which will do this live during the measurement. This algorithm presents the disadvantage of keeping only  the first stop photon and remove all subsequent stop photons which should count in a proper g2 but are lost by the algorithm. The concrete effect is an exponential decay artefact which will appear obviously when using long correlation windows or performing a measurement with high photon rates.
* ```fast=False``` is typically an order of magnitude slower than fast mode. This algorithm will consider all possible pairs of start-stop photons with a delay smaller than the correlation window. The main advantage is that for a given start photon, you count all the subsequent stop photons, removing the exponential decay artefact. You usually end up with more coincidence counts as well.


## Dynamic library compilation

readPTU uses the cffi package which handles building and interfacing of C dynamic libraries. The C code is in the _readTTTRRecords-for-import.c_ file, which is imported in Python and compiled by executing the _\_readTTTRRecords_build.py_ file.

If there is no dynamic library matching your system, you should simply execute _\_readTTTRRecords_build.py_ in python. You still need to have a compiler installed. If it does not work, you should explore cffi's [documentation](http://cffi.readthedocs.io/en/latest/overview.html?highlight=compile#abi-versus-api). Some hints for Windows: [here](https://stackoverflow.com/questions/16787649/how-to-configure-python-cffi-library-to-use-mingw) and [here](http://preshing.com/20141108/how-to-install-the-latest-gcc-on-windows/). On Linux, simply installing gcc should work.



## Changelog
Raphaël Proux - 12/10/2017 - first release with T2 mode only and support of Hydraharp v2 tested, Hydraharp v1 and Picohard untested (but implemented).
