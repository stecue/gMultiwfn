# About *gMultiwfn*
gMultiwfn is an unofficial and (maybe) enhanced gfortran port of the popular wavefunction analyzing software *[Multiwfn](http://sobereva.com/multiwfn)* developed by Tian Lu. This gfortran port is maintained by Xing Yin (stecue@gmail.com). Email Xing or [open an issue on the github](https://github.com/stecue/gMultiwfn/issues) (__*strongly preferred!*__) on the github if you find a bug or need a new additional feature. You can also discuss related topics in the [official Chinese forum](http://bbs.keinsci.com/forum-112-1.html) for `Multiwfn`.

# Binary Packages
RPM builds for openSUSE, Fedora and CentOS will be released soon. Please follow the following steps to install from the source for now.

# Compile from source
`gMultiwfn` uses the standard [GNU Build System](https://en.wikipedia.org/wiki/GNU_Build_System). If you are familiar with `./configure` and `make && make install`, building gMultiwfn is very easy. The following steps sevre as a general guideline and you can make changes accordingly if you are an expert on GNU autotools or have special needs.

## Download the package
The source tarball can be found [here](https://github.com/stecue/gMultiwfn/releases). 

## Requirements
* You can use either Intel Fortran Compiler (`ifort`) or `gfortran`.
* If you choose `ifort`, make sure Intel Math Kernel Library (Intel MKL) is installed. If not sure, you can just try to continue building first because MKL is usually bundled and installed with `ifort` by default.
* If you choose `gfortran`, make sure you have lapack/blas and their development files (usually named as `lapack-devel` and `blas-devel` *or* `lapack-dev` and `blas-dev` in your distro's repository) installed. The optimized LAPACK/BLAS implementations (see below) will not be searched and used by the `configure` script during building.

## Basic build procedure
1. Open a terminal and go to the directory where the tarball is downloaded, or move the tarball to your current directory.
2. Unzip and extract files from the source tarball. If the name of the tarball is `gMultiwfn-3.3.9-1.tar.gz`, the command would be:
```
tar -xvzf gMultiwfn-3.3.9-0.tar.gz
```
3. Make a separate build direcotry under the original "root" source directory:
```
cd gMultiwfn-3.3.9-0
mkdir build
```
4. Run the `configure` script to configure the installation path and generate the make files. The binary (`Multiwfn`) will always be installed to a `bin` directory but you can specify the where the `bin` direcory is located. For `$HOME/bin`, type:
```
../configure --prefix=$HOME
```
5. make and install. You can specify the number of processes to be used by `make`. For example, to fully utilize the power of a 4-core/8-thread CPU, type:
```
make -j8 && make install
```
6. Now `Multiwfn` should be found at `$HOME/bin`. Add `$HOME/bin` to your `PATH` environment variable if not added before. Type `Multiwfn` to see if the building is successful.

# Switch to a faster LAPACK/BLAS implementation
`gMultiwfn` is dynamically linked to `lapack` and `blas` by default. The reference implementations of `lapack` and `blas` are usually the slowest and in a lot of cases they can be safely replaced by optimized implementations such as `OpenBLAS` and `ATLAS` using the steps descripted below. Note that installing OpenBLAS or ATLAS is beyond the scope of this document and please refer to your distro's manual on that information.

## Use `update-alternatives` to make a system-wide change
Use `man` or refer to your disto's manual on the usage. You may need root user's privilege to run this command.

## Use `LD_PRELOAD` to change just for `gMultiwfn`
`LD_PRELOAD` is the enviroment variable to force the dynamic linker in Linux to use a certain version of shared libraries (.so files). It provides a quick solution to try a new libarary without being asked for root user privilege. Assuming `OpenBLAS` is installed to `/path/to/libopenblas.so`, `LD_PRELOAD=/path/to/libopenblas.so Multiwfn` will start `Multiwfn` with the optimized OpenBLAS. To avoid typing the extra letters every time, simply add `alias gMultiwfn='LD_PRELOAD=/path/to/libopenblas.so Multiwfn'` to your bash initialization file (usually `~/.bashrc`). `gMultiwfn` will be equivalent to `LD_PRELOAD=/path/to/libopenblas.so Multiwfn'` for all *new* termnal windows.

# (Must Read) Differences between *Multiwfn* and *gMultiwfn*
1. `gMultiwfn` does not contain the GUI which is base on the closed-source DISLIN library.
2. `gMultiwfn` supports the `OMP_NUM_THREADS` enviroment variable and the traditional `nthreads` parameter from `settings.ini` or the interactive input (the hidden option `1000` at the main menu). `nthreads` has the higher priority. Only when `nthreads` is set to `0` (which is the default setting without a `settings.ini` file), `OMP_NUM_THREADS` determines the number of threads to use. For example, to use 8 threads without `settings.ini`, just run
```
OMP_NUM_THREADS=8 gMultiwfn
```
3. Depending on the version of OpenBLAS, you may need to use the serial version or the OpenMP version of OpenBLAS with `OMP_NUM_THREADS`. See discussions [here](https://groups.google.com/forum/#!topic/openblas-users/W6ehBvPsKTw) and [here](https://github.com/xianyi/OpenBLAS/issues/208).
