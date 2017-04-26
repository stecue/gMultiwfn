AUTOMAKE_OPTIONS = foreign
bin_PROGRAMS = Multiwfn
Multiwfn_SOURCES = Lebedev-Laikov.F DFTxclib.F define.f90 util.f90 function.f90 sub.f90 fileIO.f90 spectrum.f90 DOS.f90 wfn.f90 population.f90 compana.f90 bondana.f90 topology.f90 excittransana.f90 otherfunc.f90 otherfunc2.f90 surfana.f90 procgriddata.f90 AdNDP.f90 fuzzyana.f90 CDA.f90 basinana.f90 atmraddens.f90
Multiwfn_FCFLAGS = -fopenmp -ffree-line-length-0 -std=legacy
Multiwfn_FFLAGS = -fopenmp -ffree-line-length-0 -std=legacy
#Multiwfn_FCFLAGS = -O2 -fopenmp -diag-disable 8290,8291,6371 -threads -openmp-link=static
#Multiwfn_FFLAGS = -O2 -fopenmp -diag-disable 8290,8291,6371 -threads -openmp-link=static
#Multiwfn_LDFLAGS = -mkl -static-intel 
Multiwfn_LDFLAGS = -llapack -lblas
