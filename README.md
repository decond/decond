Decond
=====

[![Join the chat at https://gitter.im/kmtu/decond](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/kmtu/decond?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
Decond is a framework to analyze the transport properties in electrolyte solution with various decomposition methods.

Currently, the decomposition by species and by spatial correlation is implemented and is based on [J. Chem. Phys. 141, 044126 (2014)](http://dx.doi.org/10.1063/1.4890741).

Project structure
-----
This project consists of two parts:

1. Decomposition of the velocity time correlation functions of ions calculated from the [Gromacs](http://www.gromacs.org/) MD trajectories. It is done throught the Fortran program `fortran/decompose_mpi`.
2. Analysis of the decomposed correlation data. It is done through the Python scripts in `python/` folder.


Dependency
-----
Fortran:
   - [XTC library](http://www.gromacs.org/Developer_Zone/Programming_Guide/XTC_Library):
     Please use the customized version in `fortran/lib/xdrfile-1.1.1-d`, which has been modified to double precision.
   - [HDF5](http://www.hdfgroup.org/HDF5/):
     You can use the system built-in HDF5 if available, or download it from the [HDF5 group page](http://www.hdfgroup.org/HDF5/).

Python:
   - Version >= 3.0
   - [SciPy](http://www.scipy.org/)

Using the [Anaconda Python Distribution](http://continuum.io/downloads#34) is the easiest way, which includes many packages useful for scientific calculations including SciPy.


General setup instructions
-----
The general idea of the setup process is:

1. Compile and install the customized XTC library `fortran/lib/xdrfile-1.1.1-d`.
2. Compile and install the HDF5 library, if it is not available on your system already.
3. Manually edit the file `fortran/Makefile` to suit your environment and then compile with `make`.
4. Download and install [Anaconda Python Distribution](http://continuum.io/downloads#34). (Rememeber to download the Python 3+ version)

More details can be found on the [wiki page](https://github.com/kmtu/decond/wiki).

----------
This project follows [semantic versioning](http://semver.org/)
