# readPTU library

PicoQuant uses a specific homemade file format called PTU for its time-tag time-resolved measurements. This format is subsequent to the former .pt2 or .pt3 and can handle both T2 and T3 acquisition modes for a variety of devices (Hydraharp, Picoharp, Timeharp, etc.). 

For now the library handles Hydraharp V2 and Picoharp files in T2 mode only — which are the only devices we have available for testing.

## Installation

To install this package, you need to have a compiler present on your system (`gcc`on Linux, `clang` on macOS or `vc` on Windows). Explicitly, for Windows, this looks like (since it isn't as straightforward as using a package manager):

1. Go to the [Visual Studio downloads page](https://visualstudio.microsoft.com/downloads)

2. Scroll to the bottom and expand the "Tools for Visual Studio 20XX" tab

3. Download "Build Tools for Visual Studio 20XX"

4. Install the "Desktop development with C++" package
	- This should include the 3 core libraries ("C++ Build Tools core features", "C++ 20XX Redistributable Update", "C++ core desktop features") as well as "MSVC build tools", "Windows 10 SDK", and "C++ CMake tools".

5. Restart if required

You can run the `setup.py` file using:
```
python setup.py install
```
The package should compile the libraries when installed (creates `.pyd` files on windows, `.so` files on linux, haven't tested on mac), and copy them over to the proper install location.

## Use example
To start, you should import the readPTU library you should open the PTU file using a PTUfile() object. Constructing a PTUfile() object automatically opens the file and leaves it open for further analysis. As a first step, you can then print the header which will bring useful information on the measurement, like measurement time, number of records, record format type, etc.
```
ptu_file = PTUfile(r'path/to/the/file')
ptu_file.print_header()
```

To analyse the file, e.g. to calculate a timetrace or a g2, you should use the PTUmeasurement() class. To construct a PTUmeasurement() object, you need to provide an open PTUfile object (```ptu_file``` in our example).
```
ptu_meas = PTUmeasurement(ptu_file)
```

The PTUmeasurement() object has methods to analyse the measurement. For example, to calculate a timetrace with a time bin resolution of 1 second:
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

If there is no dynamic library matching your system, you should simply execute _\_readTTTRRecords_build.py_ in python. You still need to have a compiler installed. If it does not work, you should explore cffi's documentation for [platform-specific installation instructions](https://cffi.readthedocs.io/en/latest/installation.html#windows-regular-32-bit). Some hints for Windows: set up visual studio from the website of Microsoft and using ```import setuptools``` has worked... You can have a look [here](https://stackoverflow.com/questions/16787649/how-to-configure-python-cffi-library-to-use-mingw) and [here](http://preshing.com/20141108/how-to-install-the-latest-gcc-on-windows/) as well. On Mac or Linux, simply installing gcc should work, if it is not already there.



## Changelog
Raphaël Proux - 12/10/2017 - first release with T2 mode only and support of Hydraharp v2 tested, Hydraharp v1 and Picohard untested (but implemented).

Raphaël Proux - 18/10/2017 - tested on windows 7 64 bits! It works!! (you need to add ```import setuptools``` in the python build file though).

Guillem Ballesteros-Garcia and Raphaël Proux - 18/12/2017 - very optimised version with new algorithms (ring, symmetric) and new functionnalities (timetrace channel selection, post-selection). Also was put in a package style to import it easily on new machines.
