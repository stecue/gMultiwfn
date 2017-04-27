# gMultiwfn

## About
gMultiwfn is the (unofficial) gfortran port of the popular wavefunction analyzing software [Multiwfn](http://sobereva.com/multiwfn) developed by Tian Lu. This gfortran port is maintained by Xing Yin (stecue@gmail.com). Email Xing or [open an issue on the github (*strongly preferred!*)](https://github.com/stecue/gMultiwfn/issues) on the github if you find any bug or need a new additional feature. You can also 

## Download
The source tarball can be found [here](http://sobereva.com/multiwfn). RPM builds for openSUSE, Fedora and CentOS will be released soon.

## Install
`gMultiwfn` uses the standard GNU autotools. Make sure you have lapack/blas and their development files (usually named as `lapack-devel` and `blas-devel` (or `lapack-dev` and `blas-dev`) in your distro's repository) installed before trying to build `gMultiwfn`.

## Switch to a faster lapack/blas implimentation.
`LD_PRELOAD=/path/to/libblas.so.3`
