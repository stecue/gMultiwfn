# gMultiwfn
*The (unofficial) gfortran port of Multiwfn*

## About
gMultiwfn is an unofficial and (maybe) enhanced gfortran port of the popular wavefunction analyzing software [Multiwfn](http://sobereva.com/multiwfn) developed by Tian Lu. This gfortran port is maintained by Xing Yin (stecue@gmail.com). Email Xing or [open an issue on the github](https://github.com/stecue/gMultiwfn/issues) (__*strongly preferred!*__) on the github if you find a bug or need a new additional feature.

## Installation
RPM builds for openSUSE, Fedora and CentOS will be released soon. Please follow the following steps to install from the source for now.

### Download the package
The source tarball can be found [here](http://sobereva.com/multiwfn). 

### Compile from source
`gMultiwfn` uses the standard GNU autotools. Make sure you have lapack/blas and their development files (usually named as `lapack-devel` and `blas-devel` (or `lapack-dev` and `blas-dev`) in your distro's repository) installed before trying to build `gMultiwfn`.
Install a binary package

## Switch to a faster lapack/blas implementation
`LD_PRELOAD=/path/to/libopenblas.so`

## (Must Read) Differences between Multiwfn and gMultiwfn
