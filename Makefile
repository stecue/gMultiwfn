OPT = -O2 -openmp -diag-disable 8290,8291,6371 -threads -multiple-processes=4 -openmp-link=static
LIB = -mkl -static-intel
FC = ifort
EXE = Multiwfn

default:
	$(FC) $(OPT) Lebedev-Laikov.F DFTxclib.F define.f90 util.f90 function.f90 sub.f90 fileIO.f90 spectrum.f90 DOS.f90 wfn.f90 population.f90 compana.f90 bondana.f90 topology.f90 excittransana.f90 otherfunc.f90 otherfunc2.f90 surfana.f90 procgriddata.f90 AdNDP.f90 fuzzyana.f90 CDA.f90 basinana.f90 atmraddens.f90 $(LIB) -o $(EXE)

clean:
	rm -f $(EXE) *.o *.mod
