!!----------------- plot TDOS/PDOS/OPDOS, works for Gaussian output file and .fch/.molden file, but not for .wfn
!For .out or plain text file, only one type of spin MOs will be loaded and processed, and then we will not consider spin type
!For .fch/.molden, user can switch spin type anytime
subroutine plotdos
use defvar
use util
use function
implicit real*8 (a-h,o-z)
integer,parameter :: nfragmax=10
integer,parameter :: num2Dpoints=200 !The number of points constituting the X-axis of 2D LDOS
real*8 :: curvexpos(num1Dpoints),TDOScurve(num1Dpoints),OPDOScurve(num1Dpoints),PDOScurve(num1Dpoints,nfragmax),LDOScurve(num1Dpoints)
real*8 :: LDOSxpos(num2Dpoints)
!All ?DOSliney share DOSlinex(:) as X axis
real*8,allocatable :: str(:),FWHM(:),DOSlinex(:),TDOSliney(:),PDOSliney(:,:),OPDOSliney(:),LDOSliney(:)
real*8,allocatable :: compfrag(:,:) !i,k element is the MPA composition of fragment k in MO i
real*8,allocatable :: OPfrag12(:) !Overlap population between fragment 1 and 2
real*8,allocatable :: LDOScomp(:) !Composition at a point of each orbital
real*8,allocatable :: LDOSptscomp(:,:) !Composition of each MO, ipt in a given line
real*8,allocatable :: LDOS2Dmap(:,:) !LDOS curve, ipt in a given line
integer :: nfragDOS(nfragmax)
integer,allocatable :: fragDOS(:,:) !The index of basis functions in fragments. nfragDOS is the number of basis functions in them (0=undefined)
real*8,pointer :: tmpmat(:,:)
real*8 HOMOlevx(2),HOMOlevy(2)
real*8,allocatable :: MOene_dos(:),MOocc_dos(:) !Using the ene/occ in this to plot DOS, the values are scaled when changing between a.u. and eV. The original MOene is remain unchanged
integer,allocatable :: selorbarr(:)
character clegend*960 !(10+2) lines * 80 character per line
character unitstr*5,c200tmp*200,c80tmp*80
character :: TDOSstring*80="TDOS",OPDOSstring*80="OPDOS"
character(len=80), dimension(nfragmax) :: PDOSstring = [character(len=80):: "PDOS frag.1","PDOS frag.2","PDOS frag.3","PDOS frag.4", &
                 "PDOS frag.5","PDOS frag.6","PDOS frag.7","PDOS frag.8","PDOS frag.9","PDOS frag.10"]
integer :: ishowPDOSline(nfragmax),ishowPDOScurve(nfragmax)
integer :: iclrPDOS(nfragmax)=(/ 1,3,10,14,12,9,13,11,6,7 /)
defFWHM=0.05D0 !Default FWHM
ipopmethod=1 !The method calculated OPDOS, =1 Mulliken, =3 SCPA, stout-politzer is too bad so don't consider it
ibroadfunc=2 !Default is Gaussian function
scalecurve=0.1D0 !Multiply curves with this value
enelow=-0.8D0 !Energy range, a.u.
enehigh=0.2D0
stepx=0.1D0
stepy=2
gauweigh=0.5D0 !The weight of Gaussian in Pseudo-Voigt function
nfragDOS=0
ishowTDOScurve=1
ishowTDOSline=1
ishowPDOSline=0
ishowPDOScurve=0
ishowOPDOScurve=0
ishowOPDOSline=0
ishowlegend=1
ishowHOMOlev=1
iunitx=1
unitstr=" a.u."
ispin=0 !restricted wavefunction
if (wfntype==1.or.wfntype==4) ispin=1 !Unrestricted wavefunction, output alpha part by default
iusersetYleft=0 !If user has set lower and upper range of Y axis by himself
Yrightscafac=0.5D0 !Scale factor relative to left Y-axis of OPDOS (right Y-axis)
yxratio=1D0

ireadgautype=1
if (ifiletype==0) then
    !Read energy level information from text file, the first number in first row define how many energy levels
    !in there, the second number in first row if equals to 1, means below data are only energies, if equals to 2,
    !means both strength and FWHM also present.
    open(10,file=filename,status="old")
    call loclabel(10,"Gaussian",igauout)
    rewind(10)
    if (igauout==1) then
        write(*,*) "This is Gaussian output file"
        if (ireadgautype==1) then !Read energy level from Gaussian output
            call loclabel(10,"NBsUse=")
            read(10,*) c200tmp,nbasis
            nmo=nbasis
            allocate(MOene(nmo),MOocc(nmo),str(nmo),FWHM(nmo))
            call loclabel(10,"Orbital energies and kinetic energies",ifound) !First assume this is close-shell
            if (ifound==1) then
                read(10,*)
                read(10,"(a)") c200tmp
                do i=1,nbasis
                    read(10,"(a21)",advance="no") c200tmp
                    read(10,*) MOene(i)
                    MOocc(i)=0
                    if (index(c200tmp,'O')/=0) MOocc(i)=2
                end do
                call loclabel(10,"Orbital energies and kinetic energies (beta)",ifound)
                if (ifound==1) then
                    where(MOocc==2) MOocc=1
                    write(*,*) "Read which type of orbitals? 1:alpha 2:beta"
                    read(*,*) inp
                    if (inp==2) then !Read beta energies, overlay read alpha counterpart
                        read(10,*)
                        read(10,*)
                        do i=1,nbasis
                            read(10,"(a21)",advance="no") c200tmp
                            read(10,*) MOene(i)
                            MOocc(i)=0
                            if (index(c200tmp,'O')/=0) MOocc(i)=1
                        end do
                    end if
                end if            
                write(*,*) "Read orbital energy from the file"
            else
                write(*,"(a)") " Error: Cannot find orbital energies from this file, don't forget using pop=full keyword"
                write(*,*)
                return
            end if
        end if
        str=1D0
        FWHM=defFWHM
    else !Plain text file
        read(10,*) nmo,inp
        allocate(MOene(nmo),MOocc(nmo),str(nmo),FWHM(nmo))
        if (inp==1) then
            do imo=1,nmo
                read(10,*) MOene(imo),MOocc(imo)
            end do
            str=1D0
            FWHM=defFWHM
        else if (inp==2) then
            do imo=1,nmo
                read(10,*) MOene(imo),MOocc(imo),str(imo),FWHM(imo)
            end do
        end if
    end if
    close(10)
    allocate(DOSlinex(3*nmo),TDOSliney(3*nmo))
else if (ifiletype==1.or.ifiletype==9) then
    allocate(str(nbasis),FWHM(nbasis)) !I assume the number of MOs taken into account is equal to nbasis
    !Allocate all arrays that may be used, don't consider if they will actually be used, because memory consuming is very little
    allocate(DOSlinex(3*nbasis),TDOSliney(3*nbasis),PDOSliney(3*nbasis,nfragmax),OPDOSliney(3*nbasis),LDOSliney(3*nbasis))
    allocate(compfrag(nbasis,nfragmax),OPfrag12(nbasis))
    allocate(fragDOS(nbasis,nfragmax+1)) !The last slot is used to exchange fragment
    allocate(LDOScomp(nbasis))
    str=1D0
    FWHM=defFWHM
else
    write(*,*) "ERROR: This input file is unsupported!"
    return
end if

!======Set from where to where are active energy levels
if (ispin==0) imoend=nmo !Text file or restricted .fch
if (ispin/=0) imoend=nbasis !For unrestricted fch or Gaussian output file, from 1 to imoend is the length of one type of spin orbitals

if (allocated(MOene_dos)) deallocate(MOene_dos,MOocc_dos) !MOene_dos is the working horse
allocate(MOene_dos(nmo),MOocc_dos(nmo))
MOene_dos=MOene
MOocc_dos=MOocc

do while(.true.) !!!!! main loop
idoPDOS=0
if (any(nfragDOS>0)) idoPDOS=1
idoOPDOS=0
if (all(nfragDOS(1:2)>0)) idoOPDOS=1

!Unknow text file doesn't contains wavefunction info, couldn't define fragment
write(*,*) "          ================ Plot density-of-states ==============="
write(*,*) "-10 Return to main menu"
write(*,*) "-5 Customize energy levels, occupations, strengths and FWHMs for specific MOs"
write(*,*) "-4 Show all orbital information"
write(*,*) "-3 Export energy levels, occupations, strengths and FWHMs to plain text file"
if (ifiletype==1.or.ifiletype==9) write(*,*) "-1 Define fragments"
if (idoOPDOS==1) then
    write(*,*) "0 Draw TDOS+PDOS+OPDOS graph!"
else if (idoPDOS==1) then
    write(*,*) "0 Draw TDOS+PDOS graph!"
else
    write(*,*) "0 Draw TDOS graph!" !Reading text file can only draw spinless TDOS, because they impossible to define fragment
end if
if (ibroadfunc==1) write(*,*) "1 Select broadening function, current: Lorentzian"
if (ibroadfunc==2) write(*,*) "1 Select broadening function, current: Gaussian"
if (ibroadfunc==3) write(*,*) "1 Select broadening function, current: Pseudo-Voigt"
write(*,"(a,f14.5,a,f14.5,a)") " 2 Set energy range, current:",enelow," to",enehigh,unitstr
if (maxval(FWHM)==minval(FWHM)) then
    write(*,"(a,f10.5,a)") " 3 Set full width at half maximum (FWHM), current:",FWHM(1),unitstr
else
    write(*,"(a,f10.5)") " 3 Set full width at half maximum (FWHM), current: Orbital dependent"
end if
write(*,"(a,f10.5)") " 4 Set scale ratio for DOS curve, current:",scalecurve
if (ibroadfunc==3) write(*,"(a,f10.5)") " 5 Set Gaussian-weighting coefficient, current:",gauweigh
if (ispin==1) write(*,*) "6 Switch spin, current: Alpha"
if (ispin==2) write(*,*) "6 Switch spin, current: Beta"
if (ipopmethod==1.and.(ifiletype==1.or.ifiletype==9)) write(*,*) "7 Switch method for calculating PDOS, current: Mulliken"
if (ipopmethod==3.and.(ifiletype==1.or.ifiletype==9)) write(*,*) "7 Switch method for calculating PDOS, current: SCPA"
write(*,*) "8 Switch unit between a.u. and eV, current:"//unitstr
write(*,*) "10 Draw local DOS for a point"
write(*,*) "11 Draw local DOS along a line"

read(*,*) isel

if (isel==-10) then
    exit
else if (isel==-5) then
    do while(.true.)
        write(*,*)
        write(*,*) "0 Return"
        write(*,*) "1 Set orbital energies for specific orbitals"
        write(*,*) "2 Set occupation numbers for specific orbitals"
        write(*,*) "3 Set strengths for specific orbitals"
        write(*,*) "4 Set FWHMs for specific orbitals"
        read(*,*) isel
        
        if (isel==0) then
            exit
        else 
            write(*,"(a)") " Input orbital indices. e.g. 1,3-6,8,10-11 means the orbital 1,3,4,5,6,8,10,11 will be selected"
            read(*,"(a)") c200tmp
            call str2arr(c200tmp,nselorb)
            if (allocated(selorbarr)) deallocate(selorbarr)
            allocate(selorbarr(nselorb))
            call str2arr(c200tmp,nselorb,selorbarr)
            if (isel==1) then
                write(*,*) "Set their energies to which value? e.g. -0.13"
                write(*,*) "Note: The value should be given in"//unitstr
                read(*,*) enetmp
                if (ispin==2) selorbarr=selorbarr+nbasis
                do imoidx=1,nselorb
                    imo=selorbarr(imoidx)
                    MOene_dos(imo)=enetmp
                end do
            else if (isel==2) then
                write(*,*) "Set their occupation numbers to which value? e.g. 2.0"
                read(*,*) occtmp
                if (ispin==2) selorbarr=selorbarr+nbasis
                do imoidx=1,nselorb
                    imo=selorbarr(imoidx)
                    MOocc_dos(imo)=occtmp
                end do
            else if (isel==3) then
                write(*,*) "Set their strength to which value? e.g. 1.0"
                read(*,*) strtmp
                do imoidx=1,nselorb
                    imo=selorbarr(imoidx)
                    str(imo)=strtmp
                end do
            else if (isel==4) then
                write(*,*) "Set their FWHM to which value? e.g. 0.05"
                read(*,*) FWHMtmp
                do imoidx=1,nselorb
                    imo=selorbarr(imoidx)
                    FWHM(imo)=FWHMtmp
                end do
            end if
        end if
    end do
else if (isel==-4) then
    iFermi=0
    if (ispin==1) write(*,*) "Below orbitals are Alpha type"
    if (ispin==2) write(*,*) "Below orbitals are Beta type"
    do imo=1,imoend
        irealmo=imo
        if (ispin==2) irealmo=imo+nbasis
        write(*,"('#',i5,' Energy(',a,'):',f14.6,' Occ:',f8.5,' Str:',f8.5,' FWHM:',f8.5)") imo,trim(unitstr),MOene_dos(irealmo),MOocc_dos(irealmo),str(imo),FWHM(imo)
        if (imo>1) then
            if (MOocc_dos(irealmo)==0D0.and.MOocc_dos(irealmo-1)/=0D0) iFermi=irealmo-1
        end if
    end do
    if (iFermi/=0) write(*,"(' Fermi energy level:',f12.6,a)") MOene_dos(iFermi),unitstr
    write(*,*)
else if (isel==-3) then
    open(10,file="orginfo.txt",status="replace")
    write(10,"(2i6)") imoend,2
    do imo=1,imoend
        irealmo=imo
        if (ispin==2) irealmo=imo+nbasis
        write(10,"(f14.6,3f12.6)") MOene_dos(irealmo),MOocc_dos(irealmo),str(imo),FWHM(imo)
    end do
    close(10)
    write(*,"(a)") " The energy levels, occupation numbers, strengths, FWHMs have been exported to orginfo.txt in current directory, &
    you can modify it and then load it into Multiwfn again."
    write(*,*) "Note: The energy unit of energy levels and FWHMs in this file is in"//unitstr
    write(*,*)
else if (isel==-1) then
    write(*,*) "           ----------------- Define fragments -----------------"
    write(*,"(a)") " Note: Up to 10 fragments can be defined for plotting PDOS, but OPDOS will only be plotted for fragments 1 and 2"
    do while(.true.)
        do ifrag=1,nfragmax
            if (nfragDOS(ifrag)==0) then
                write(*,"(' Fragment',i5,', has not been defined')") ifrag
            else
                write(*,"(' Fragment',i5,', the number of basis functions:',i6)") ifrag,nfragDOS(ifrag)
            end if
        end do
        write(*,*) "Input fragment index to define it, e.g. 2"
        write(*,*) "Input a negative index can unset the fragment, e.g. -2"
        write(*,*) "Input two indices can exchange the two fragments, e.g. 1,4"
        write(*,*) "Input ""e"" can export current fragment setting to DOSfrag.txt in current folder"
        write(*,*) "Input ""i"" can import fragment setting from DOSfrag.txt in current folder"
        write(*,*) "To return to the last menu, input 0"
        read(*,"(a)") c80tmp
        if (index(c80tmp(1:len_trim(c80tmp)),' ')/=0.or.c80tmp==" ") then
            write(*,*) "Input error!"
            write(*,*)
        else if (c80tmp=='e') then
            open(10,file="DOSfrag.txt",status="replace")
            do ifrag=1,nfragmax
                write(10,*)
                write(10,"(' #Fragment:',i4,'   nbasis:',i8)") ifrag,nfragDOS(ifrag)
                write(10,"(8i8)") fragDOS(1:nfragDOS(ifrag),ifrag)
            end do
            close(10)
            write(*,*) "Export finished!"
            write(*,*)
        else if (c80tmp=='i') then
            open(10,file="DOSfrag.txt",status="old")
            do ifrag=1,nfragmax
                read(10,*)
                read(10,*) c80tmp,inouse,c80tmp,nfragDOS(ifrag)
                read(10,"(8i8)") fragDOS(1:nfragDOS(ifrag),ifrag)
            end do
            close(10)
            write(*,*) "Import finished!"
            write(*,*)
        else if (index(c80tmp,',')==0) then !Set specific fragment
            read(c80tmp,*) ifragsel
            if (ifragsel==0) then
                exit
            else if (ifragsel>nfragmax) then
                write(*,*) "ERROR: The index exceeded upper limit!"
            else if (ifragsel<0) then
                nfragDOS(abs(ifragsel))=0
            else !deffrag routine is only able to deal with global array frag1/2, so we use frag1 as intermediate array
                allocate(frag1(nfragDOS(ifragsel)))
                frag1(:)=fragDOS(1:nfragDOS(ifragsel),ifragsel)
                call deffrag(1)
                if (allocated(frag1)) then
                    nfragDOS(ifragsel)=size(frag1)
                    fragDOS(1:nfragDOS(ifragsel),ifragsel)=frag1(:)
                    deallocate(frag1)
                else
                    nfragDOS(ifragsel)=0
                end if
            end if
        else !Exchange fragments
            read(c80tmp,*) ifragsel,jfragsel
            ntmp=nfragDOS(jfragsel)
            nfragDOS(jfragsel)=nfragDOS(ifragsel)
            nfragDOS(ifragsel)=ntmp
            fragDOS(:,nfragmax+1)=fragDOS(:,jfragsel)
            fragDOS(:,jfragsel)=fragDOS(:,ifragsel)
            fragDOS(:,ifragsel)=fragDOS(:,nfragmax+1)
        end if
    end do
else if (isel==1) then
    write(*,*) "1 Lorentzian"
    write(*,*) "2 Gaussian"
    write(*,*) "3 Pseudo-Voigt"
    read(*,*) ibroadfunc
else if (isel==2) then
    if (iunitx==1) then
        write(*,*) "Input lower, upper limits and stepsize between legends (in a.u.)"
        write(*,*) "e.g. -1.5,0.2,0.3"
    else if (iunitx==2) then
        write(*,*) "Input lower, upper limits and stepsize between legends (in eV)"
        write(*,*) "e.g. -20,5,2"
    end if
    read(*,*) enelow,enehigh,stepx
else if (isel==3) then
    write(*,*) "Input a value"
    read(*,*) FWHMtmp
    if (FWHMtmp<0D0) write(*,*) "Error: The value should larger than zero, input again"
    FWHM=FWHMtmp
else if (isel==4) then
    write(*,*) "Input a value"
    read(*,*) scalecurve
    if (scalecurve<0D0) write(*,*) "Error: The value should larger than zero, input again"
else if (isel==5) then
    write(*,*) "Input a value"
    read(*,*) gauweigh
    if (gauweigh<0D0) write(*,*) "Error: The value should larger than zero, input again"    
else if (isel==6) then
    if (ispin==1) then
        ispin=2
    else
        ispin=1
    end if
else if (isel==7) then
    if (ipopmethod==1) then
        ipopmethod=3
    else
        ipopmethod=1
    end if
else if (isel==8) then
    if (iunitx==1) then !a.u.->eV
        iunitx=2
        MOene_dos=MOene_dos*au2eV
        FWHM=FWHM*au2eV
        enelow=enelow*au2eV
        enehigh=enehigh*au2eV
        unitstr=" eV"
        !After change the unit, in principle, the curve (and hence Y-range) will be automatically reduced by 27.2114.& 
        !str should also be reduced by 27.2114 so that the discrete line can be properly shown in the graph range &
        !To compensate the reduce of str, scalecurve thus be augmented by corresponding factor
        str=str/au2eV
        scalecurve=scalecurve*au2eV
        stepx=stepx*au2eV
        stepy=stepy/au2eV
    else !eV->a.u.
        iunitx=1
        MOene_dos=MOene_dos/au2eV
        FWHM=FWHM/au2eV
        enelow=enelow/au2eV
        enehigh=enehigh/au2eV
        unitstr=" a.u."
        str=str*au2eV
        scalecurve=scalecurve/au2eV
        stepx=stepx/au2eV
        stepy=stepy*au2eV
    end if
    
    
else if (isel==0.or.isel==10) then

    if (isel==10) then
        write(*,*) "Input coordinate (in Bohr), e.g. 1.0,1.5,0.2"
        read(*,*) x,y,z
    end if

    !======Generate fragment composition and overlap population=======
    tmpmat=>cobasa
    if (ispin==2) tmpmat=>cobasb
    
    !Reset display setting
    ishowPDOScurve=0
    ishowPDOSline=0
    do ifrag=1,nfragmax
        if (nfragDOS(ifrag)>0) then
            ishowPDOScurve(ifrag)=1
            ishowPDOSline(ifrag)=1
        end if
    end do
    ishowOPDOScurve=0
    ishowOPDOSline=0
    if (idoOPDOS==1) then
        ishowOPDOScurve=1
        ishowOPDOSline=1
    end if
    ishowLDOScurve=0
    ishowLDOSline=0
    if (isel==10) then
        ishowLDOScurve=1
        ishowLDOSline=1
    end if
    
    if (idoPDOS==1.or.isel==10) then
        write(*,*) "Calculating orbital composition, please wait..."
        do istart=1,nbasis
            enetmp=MOene_dos(istart)
            if (ispin==2) enetmp=MOene_dos(istart+nbasis)
            if (enetmp>enelow-3*FWHM(istart)) exit
        end do
        do iend=nbasis,1,-1
            enetmp=MOene_dos(iend)
            if (ispin==2) enetmp=MOene_dos(iend+nbasis)
            if (enetmp<enehigh+3*FWHM(iend)) exit
        end do
        write(*,"(' Note: MOs from',i7,' to',i7,' are taken into account')") istart,iend
        compfrag=0
        LDOScomp=0
        OPfrag12=0
        if (idoPDOS==1) then
!$OMP PARALLEL DO SHARED(compfrag,OPfrag12) PRIVATE(ifrag,imo,allsqr,i,j,ibas,jbas) schedule(dynamic) NUM_THREADS( nthreads  )
            do imo=istart,iend !Cycle each orbital, don't use nmo, because for unrestricted wavefunction nmo=2*nbasis
                if (ipopmethod==3) allsqr=sum(tmpmat(:,imo)**2)
                do ifrag=1,nfragmax
                    if (nfragDOS(ifrag)==0) cycle
                    do i=1,nfragDOS(ifrag) !Cycle each basis in the fragment
                        ibas=fragDOS(i,ifrag)
                        if (ipopmethod==3) then !SCPA
                            compfrag(imo,ifrag)=compfrag(imo,ifrag)+tmpmat(ibas,imo)**2/allsqr
                        else !Mulliken
                            do jbas=1,nbasis !Cycle all basis, included inner&external cross term and local term (when ibas==jbas)
                                compfrag(imo,ifrag)=compfrag(imo,ifrag)+tmpmat(ibas,imo)*tmpmat(jbas,imo)*Sbas(ibas,jbas)
                            end do
                        end if
                    end do
                end do
                 !Calculate Overlap population between frag 1&2
                if (idoOPDOS==1) then
                    do i=1,nfragDOS(1)
                        ibas=fragDOS(i,1)
                        do j=1,nfragDOS(2)
                            jbas=fragDOS(j,2)
                            OPfrag12(imo)=OPfrag12(imo)+2*tmpmat(ibas,imo)*tmpmat(jbas,imo)*Sbas(ibas,jbas)
                        end do
                    end do
                end if
            end do
!$OMP END PARALLEL DO
        end if
        if (isel==10) then !Calculate LDOS for a point
            do imo=istart,iend !Cycle each orbital
                LDOScomp(imo)=fmo(x,y,z,imo)**2
            end do
        end if
    end if
    

    !======Set X position of curves==========
    enestep=(enehigh-enelow)/(num1Dpoints-1) 
    do i=1,num1Dpoints
        curvexpos(i)=enelow+(i-1)*enestep
    end do
    
    !======Generate energy levels line=======
    do imo=1,imoend
        inow=3*(imo-1)
        irealmo=imo
        if (ispin==2) irealmo=imo+nbasis
        DOSlinex(inow+1:inow+3)=MOene_dos(irealmo)
        if (isel==0) then
            TDOSliney(inow+1)=0D0
            TDOSliney(inow+2)=str(imo)
            TDOSliney(inow+3)=0D0
            do ifrag=1,nfragmax
                if (nfragDOS(ifrag)>0) then
                    PDOSliney(inow+1,ifrag)=0D0
                    PDOSliney(inow+2,ifrag)=str(imo)*compfrag(imo,ifrag)
                    PDOSliney(inow+3,ifrag)=0D0
                end if
            end do
            if (idoOPDOS==1) then
                OPDOSliney(inow+1)=0D0
                OPDOSliney(inow+2)=str(imo)*OPfrag12(imo)
                OPDOSliney(inow+3)=0D0
            end if
        else if (isel==10) then
            LDOSliney(inow+1)=0D0
            LDOSliney(inow+2)=str(imo)*LDOScomp(imo)
            LDOSliney(inow+3)=0D0
        end if
    end do
    
    !======Generate DOS curve=======
    TDOScurve=0D0
    PDOScurve=0D0
    OPDOScurve=0D0
    LDOScurve=0D0
    if (ibroadfunc==1.or.ibroadfunc==3) then !Lorentzian function, see http://mathworld.wolfram.com/LorentzianFunction.html
        do imo=1,imoend !Cycle each orbital
            irealmo=imo
            if (ispin==2) irealmo=imo+nbasis
            preterm=str(imo)*0.5D0/pi*FWHM(imo)
            do ipoint=1,num1Dpoints !Broaden imo as curve
                tmp=preterm/( (curvexpos(ipoint)-MOene_dos(irealmo))**2+0.25D0*FWHM(imo)**2 )
                if (isel==0) then
                    TDOScurve(ipoint)=TDOScurve(ipoint)+tmp
                    do ifrag=1,nfragmax
                        if (nfragDOS(ifrag)>0) PDOScurve(ipoint,ifrag)=PDOScurve(ipoint,ifrag)+tmp*compfrag(imo,ifrag)
                    end do
                    if (idoOPDOS==1) OPDOScurve(ipoint)=OPDOScurve(ipoint)+tmp*OPfrag12(imo)
                else if (isel==10) then
                    LDOScurve(ipoint)=LDOScurve(ipoint)+tmp*LDOScomp(imo)
                end if
            end do
        end do
    end if
    if (ibroadfunc==2.or.ibroadfunc==3) then
        if (ibroadfunc==3) TDOScurve=(1-gauweigh)*TDOScurve
        do imo=1,imoend !Cycle each orbital
            irealmo=imo
            if (ispin==2) irealmo=imo+nbasis
            gauss_c=FWHM(imo)/2D0/sqrt(2*dlog(2D0))
            gauss_a=str(imo)/(gauss_c*sqrt(2D0*pi))
            do ipoint=1,num1Dpoints !Broaden imo as curve
                !Gaussian function, see http://en.wikipedia.org/wiki/Gaussian_function
                tmp=gauss_a*dexp( -(curvexpos(ipoint)-MOene_dos(irealmo))**2/(2*gauss_c**2) )
                if (ibroadfunc==3) tmp=gauweigh*tmp !Combine Lorentizan and Gaussian function
                if (isel==0) then
                    TDOScurve(ipoint)=TDOScurve(ipoint)+tmp
                    do ifrag=1,nfragmax
                        if (nfragDOS(ifrag)>0) PDOScurve(ipoint,ifrag)=PDOScurve(ipoint,ifrag)+tmp*compfrag(imo,ifrag)
                    end do
                    if (idoOPDOS==1) OPDOScurve(ipoint)=OPDOScurve(ipoint)+tmp*OPfrag12(imo)
                else if (isel==10) then
                    LDOScurve(ipoint)=LDOScurve(ipoint)+tmp*LDOScomp(imo)
                end if
            end do
        end do
    end if
    TDOScurve=TDOScurve*scalecurve
    PDOScurve=PDOScurve*scalecurve
    OPDOScurve=OPDOScurve*scalecurve
    LDOScurve=LDOScurve*scalecurve
    
    idraw=1
    isavepic=0
    if (iusersetYleft==0) then !Y axis range was not set by user, we use recommended value
        if (isel==0) then
            yupperleft=1.1D0*maxval(TDOScurve)
            ylowerleft=0
            if (idoPDOS==1) ylowerleft=minval(PDOScurve(:,:)) !PDOS may be negative
            if (idoOPDOS==1) ylowerleft=-yupperleft/2 !OPDOS may be large negative value
            if (ylowerleft>0) ylowerleft=0D0 !Don't allow lower plotting limit >0
            stepy=nint((yupperleft-ylowerleft)*10)/100D0
        else if (isel==10) then
            ylowerleft=0
            yupperleft=1.1D0*maxval(LDOScurve)
            stepy=(yupperleft-ylowerleft)/10
        end if
    end if
    ylowerright=ylowerleft*Yrightscafac !Lower and upper limit for OPDOS
    yupperright=yupperleft*Yrightscafac
    
    do while(.true.)
        if (idraw==1.and.isilent==0) then
            !======Draw DOS now=======
            if (idoOPDOS==1) then
            else
            end if
            numleg=0
            ileg=1

            if (isel==0) then
                !Set legends
                if (ishowTDOScurve==1.or.ishowTDOSline==1) numleg=numleg+1
                do ifrag=1,nfragmax
                    if (nfragDOS(ifrag)==0) cycle
                    if (ishowPDOScurve(ifrag)==1.or.ishowPDOSline(ifrag)==1) numleg=numleg+1
                end do
                if (ishowOPDOScurve==1.or.ishowOPDOSline==1) numleg=numleg+1

                !Draw TDOS
                if (ishowTDOScurve==1.or.ishowTDOSline==1) then
                    ileg=ileg+1
                end if
                
                !Draw a vertical dashed line to highlight HOMO level
                if (ishowHOMOlev==1) then
                    do imo=1,imoend
                        irealmo=imo
                        if (ispin==2) irealmo=imo+nbasis
                        if (imo>1) then
                            if (MOocc_dos(irealmo)==0D0.and.MOocc_dos(irealmo-1)/=0D0) iFermi=irealmo-1
                        end if
                    end do
                    HOMOlevx=MOene_dos(iFermi)
                    HOMOlevy(1)=ylowerleft
                    HOMOlevy(2)=yupperleft
                end if
                
                !Draw PDOS of each defined fragment
                do ifrag=1,nfragmax
                    if (nfragDOS(ifrag)>0) then
                        if (ishowPDOScurve(ifrag)==1.or.ishowPDOSline(ifrag)==1) then
                            ileg=ileg+1
                        end if
                    end if
                end do

                !Draw OPDOS
                if (ishowOPDOScurve==1.or.ishowOPDOSline==1) then
                end if
            
            !Draw LDOS
            else if (isel==10) then
            end if

            if (isavepic==1) write(*,*) "Graph file has been saved to current folder with ""DISLIN"" prefix"
        end if
        idraw=0

        !======Submenu=======
        write(*,*)
        write(*,*) "                    ======== Post-process menu ========"
        if (isel==0) then !T/P/OPDOS
            write(*,*) "0 Return"
            write(*,*) "1 Show graph again"
            write(*,*) "2 Save the picture to current folder"
            write(*,*) "3 Export curve data to plain text file in current folder"
            if (idoPDOS==1) then
                write(*,"(' 4 Set Y-axis range and step for TDOS+PDOS, current:',f8.2,' to',f8.2,',',f6.2)") ylowerleft,yupperleft,stepy
            else
                write(*,"(' 4 Set Y-axis range and step for TDOS, current:',f8.2,' to',f8.2,',',f6.2)") ylowerleft,yupperleft,stepy
            end if
            if (ishowTDOScurve==1) write(*,*) "5 Toggle showing TDOS curve, current: Yes"
            if (ishowTDOScurve==0) write(*,*) "5 Toggle showing TDOS curve, current: No"
            if (ishowTDOSline==1) write(*,*) "6 Toggle showing TDOS line, current: Yes"
            if (ishowTDOSline==0) write(*,*) "6 Toggle showing TDOS line, current: No"
            if (idoPDOS==1) then
                write(*,*) "7 Toggle showing PDOS curves"
                write(*,*) "8 Toggle showing PDOS lines"
            end if
            if (idoOPDOS==1) then
                if (ishowOPDOScurve==1) write(*,*) "9 Toggle showing OPDOS curve, current: Yes"
                if (ishowOPDOScurve==0) write(*,*) "9 Toggle showing OPDOS curve, current: No"
                if (ishowOPDOSline==1) write(*,*) "10 Toggle showing OPDOS line, current: Yes"
                if (ishowOPDOSline==0) write(*,*) "10 Toggle showing OPDOS line, current: No"
            end if
            if (idoPDOS==1) write(*,*) "11 Set color for PDOS curves and lines"
            if (ishowlegend==1) write(*,*) "13 Toggle showing legend, current: Yes"
            if (ishowlegend==0) write(*,*) "13 Toggle showing legend, current: No"
            if (idoOPDOS==1) write(*,"(a,f10.5)") " 14 Set scale factor of Y-axis for OPDOS, current:",Yrightscafac
            write(*,*) "15 Toggle showing vertical dashed line to highlight HOMO level"
            write(*,*) "16 Set the texts in the legends"
            read(*,*) isel2

            if (isel2==0) then
                exit
            else if (isel2==1) then
                idraw=1
                isavepic=0
            else if (isel2==2) then
                idraw=1
                isavepic=1
            else if (isel2==3) then
                open(10,file="DOS_line.txt",status="replace")
                if (idoOPDOS==1) then
                    do i=1,3*imoend
                        write(10,"(f10.5,12f12.6)") DOSlinex(i),TDOSliney(i),OPDOSliney(i),PDOSliney(i,:)
                    end do
                else if (idoPDOS==1) then
                    do i=1,3*imoend
                        write(10,"(f10.5,11f12.6)") DOSlinex(i),TDOSliney(i),PDOSliney(i,:)
                    end do
                else
                    do i=1,3*imoend
                        write(10,"(f10.5,f12.6)") DOSlinex(i),TDOSliney(i)
                    end do
                end if
                close(10)
                open(10,file="DOS_curve.txt",status="replace")
                if (idoOPDOS==1) then
                    do i=1,num1Dpoints
                        write(10,"(f10.5,12(1PE15.6))") curvexpos(i),TDOScurve(i),OPDOScurve(i),PDOScurve(i,:)
                    end do
                else if (idoPDOS==1) then
                    do i=1,num1Dpoints
                        write(10,"(f10.5,11(1PE15.6))") curvexpos(i),TDOScurve(i),PDOScurve(i,:)
                    end do
                else
                    do i=1,num1Dpoints
                        write(10,"(f10.5,1PE15.6)") curvexpos(i),TDOScurve(i)
                    end do
                end if
                close(10)
                write(*,*) "Curve data have been written to DOS_curve.txt in current folder"
                write(*,*) "Discrete line data have been written to DOS_line.txt in current folder"
                write(*,*) "Column 1: Energy ("//trim(unitstr)//")"
                write(*,*) "Column 2: TDOS"
                if (idoOPDOS==1) then
                    write(*,*) "Column 3: OPDOS"
                    write(*,*) "Column 4-13: PDOS of fragment 1-10"
                else if (idoPDOS==1) then
                    write(*,*) "Column 3-12: PDOS of fragment 1-10"
                end if
            else if (isel2==4) then
                write(*,*) "Input lower and upper limit as well as stepsize, e.g. 0.0,2.4,0.3"
                read(*,*) ylowerleft,yupperleft,stepy
                iusersetYleft=1
                ylowerright=ylowerleft*Yrightscafac !Lower and upper limit for OPDOS. Set it here aims for immediately make effect
                yupperright=yupperleft*Yrightscafac
            else if (isel2==5) then
                if (ishowTDOScurve==0) then
                    ishowTDOScurve=1
                else
                    ishowTDOScurve=0
                end if
            else if (isel2==6) then
                if (ishowTDOSline==0) then
                    ishowTDOSline=1
                else
                    ishowTDOSline=0
                end if
            else if (isel2==7) then
                do while(.true.)
                    write(*,*)
                    write(*,*) "0 Return"
                    do ifrag=1,nfragmax
                        if (nfragDOS(ifrag)>0) then
                            if (ishowPDOScurve(ifrag)==1) write(*,"(i2,' Toggle showing PDOS curve of fragment',i3,', current: Yes')") ifrag,ifrag
                            if (ishowPDOScurve(ifrag)==0) write(*,"(i2,' Toggle showing PDOS curve of fragment',i3,', current: No')") ifrag,ifrag
                        end if
                    end do
                    read(*,*) iselfrag
                    if (iselfrag==0) exit
                    if (ishowPDOScurve(iselfrag)==0) then
                        ishowPDOScurve(iselfrag)=1
                    else
                        ishowPDOScurve(iselfrag)=0
                    end if
                end do
            else if (isel2==8) then
                do while(.true.)
                    write(*,*)
                    write(*,*) "0 Return"
                    do ifrag=1,nfragmax
                        if (nfragDOS(ifrag)>0) then
                            if (ishowPDOSline(ifrag)==1) write(*,"(i2,' Toggle showing PDOS line of fragment',i3,', current: Yes')") ifrag,ifrag
                            if (ishowPDOSline(ifrag)==0) write(*,"(i2,' Toggle showing PDOS line of fragment',i3,', current: No')") ifrag,ifrag
                        end if
                    end do
                    read(*,*) iselfrag
                    if (iselfrag==0) exit
                    if (ishowPDOSline(iselfrag)==0) then
                        ishowPDOSline(iselfrag)=1
                    else
                        ishowPDOSline(iselfrag)=0
                    end if
                end do
            else if (isel2==9) then
                if (ishowOPDOScurve==0) then
                    ishowOPDOScurve=1
                else
                    ishowOPDOScurve=0
                end if
            else if (isel2==10) then
                if (ishowOPDOSline==0) then
                    ishowOPDOSline=1
                else
                    ishowOPDOSline=0
                end if
            else if (isel2==11) then
                do while(.true.)
                    write(*,*)
                    write(*,*) "0 Return"
                    do ifrag=1,nfragmax
                        if (nfragDOS(ifrag)>0) write(*,"(' Set color for fragment',i3,', current: ',a)") ifrag,colorname(iclrPDOS(ifrag))
                    end do
                    read(*,*) iselfrag
                    if (iselfrag<0.or.iselfrag>nfragmax) then
                        write(*,*) "ERROR: Index exceeded valid range!"
                    else if (iselfrag==0) then
                        exit
                    else
                        write(*,*) "Use which color?"
                        write(*,*) "1/2/3/4/5 = Red/Green/Blue/White/Black"
                        write(*,*) "6/7/8/9/10 = Gray/Cyan/Yellow/Orange/Magenta"
                        write(*,*) "11/12/13/14 = Crimson/Dark green/Purple/Brown"
                        read(*,*) iclrPDOS(iselfrag)
                    end if
                end do
            else if (isel2==13) then
                if (ishowlegend==0) then
                    ishowlegend=1
                else
                    ishowlegend=0
                end if
            else if (isel2==14) then
                write(*,*) "Input scale factor, e.g. 0.15"
                read(*,*) Yrightscafac
                ylowerright=ylowerleft*Yrightscafac !Lower and upper limit for OPDOS. Set it here aims for immediately make effect
                yupperright=yupperleft*Yrightscafac
            else if (isel2==15) then
                if (ishowHOMOlev==0) then
                    ishowHOMOlev=1
                else
                    ishowHOMOlev=0
                end if
            else if (isel2==16) then
                do while(.true.)
                    write(*,*) "Select a term to change its legend, e.g. 3"
                    write(*,"(' -2 TDOS, current text: ',a)") trim(TDOSstring)
                    if (idoOPDOS==1) write(*,"(' -1 OPDOS, current text: ',a)") trim(OPDOSstring)
                    write(*,*) " 0 Return"
                    do ifrag=1,nfragmax
                        if (nfragDOS(ifrag)>0) write(*,"(i3,' PDOS of frag',i2,', current text: ',a)") ifrag,ifrag,trim(PDOSstring(ifrag))
                    end do
                    read(*,*) iseltmp
                    if (iseltmp==0) exit
                    write(*,*) "Input text for the legend"
                    read(*,"(a)") c80tmp
                    if (iseltmp==-2) then
                        TDOSstring=c80tmp
                    else if (iseltmp==-1) then
                        OPDOSstring=c80tmp
                    else
                        PDOSstring(iseltmp)=c80tmp
                    end if
                end do
            end if
            
        else if (isel==10) then !LDOS in 1D
            write(*,*) "0 Return"
            write(*,*) "1 Show graph again"
            write(*,*) "2 Save the picture to current folder"
            write(*,*) "3 Export curve data to plain text file in current folder"
            if (ishowLDOSline==1) write(*,*) "6 Toggle showing LDOS line, current: Yes"
            if (ishowLDOSline==0) write(*,*) "6 Toggle showing LDOS line, current: No"
            read(*,*) isel2

            if (isel2==0) then
                exit
            else if (isel2==1) then
                idraw=1
                isavepic=0
            else if (isel2==2) then
                idraw=1
                isavepic=1
            else if (isel2==3) then
                open(10,file="LDOS_line.txt",status="replace")
                do i=1,3*imoend
                    write(10,"(f10.5,1PE15.6)") DOSlinex(i),LDOSliney(i)
                end do
                close(10)
                open(10,file="LDOS_curve.txt",status="replace")
                do i=1,num1Dpoints
                    write(10,"(f10.5,1PE15.6)") curvexpos(i),LDOScurve(i)
                end do
                close(10)
                write(*,*) "Curve data have been written to LDOS_curve.txt in current folder"
                write(*,*) "Discrete line data have been written to LDOS_line.txt in current folder"
                write(*,*) "Column 1: Energy ("//trim(unitstr)//")"
                write(*,*) "Column 2: LDOS"
            else if (isel2==6) then
                if (ishowLDOSline==0) then
                    ishowLDOSline=1
                else
                    ishowLDOSline=0
                end if
            end if
        end if
    end do


!Plot local DOS along a line (2D color-filled map)
else if (isel==11) then
    if (allocated(LDOS2Dmap)) deallocate(LDOS2Dmap,LDOSptscomp)
    write(*,*) "Input the starting point (in Bohr), e.g. 1.0,0,0.2"
    read(*,*) orgx,orgy,orgz
    write(*,*) "Input the end point (in Bohr), e.g. 2.0,0,0.4"
    read(*,*) endx,endy,endz
    write(*,*) "Input the number of points along the line"
    read(*,*) numLDOSpt
    allocate(LDOS2Dmap(num2Dpoints,numLDOSpt),LDOSptscomp(nbasis,numLDOSpt))
    write(*,*) "Calculating orbital composition, please wait..."
    do istart=1,nbasis
        enetmp=MOene_dos(istart)
        if (ispin==2) enetmp=MOene_dos(istart+nbasis)
        if (enetmp>enelow-3*FWHM(istart)) exit
    end do
    do iend=nbasis,1,-1
        enetmp=MOene_dos(iend)
        if (ispin==2) enetmp=MOene_dos(iend+nbasis)
        if (enetmp<enehigh+3*FWHM(iend)) exit
    end do
!     istart=1
!     iend=nbasis
    write(*,"(' Note: MOs from',i7,' to',i7,' are taken into account')") istart,iend
    xlen=endx-orgx
    dx=xlen/(numLDOSpt-1)
    ylen=endy-orgy
    dy=ylen/(numLDOSpt-1)
    zlen=endz-orgz
    dz=zlen/(numLDOSpt-1)
    totlen=dsqrt(xlen**2+ylen**2+zlen**2)
    dlen=dsqrt(dx**2+dy**2+dz**2)
    LDOSptscomp=0D0
    do ipt=1,numLDOSpt
        x=orgx+dx*(ipt-1)
        y=orgy+dy*(ipt-1)
        z=orgz+dz*(ipt-1)
        do imo=istart,iend
            LDOSptscomp(imo,ipt)=fmo(x,y,z,imo)**2
        end do
    end do
    
    enestep=(enehigh-enelow)/(num2Dpoints-1) 
    do i=1,num2Dpoints
        LDOSxpos(i)=enelow+(i-1)*enestep
    end do
    
    LDOS2Dmap=0D0
    if (ibroadfunc==1.or.ibroadfunc==3) then !Lorentzian function, see http://mathworld.wolfram.com/LorentzianFunction.html
        do ilinept=1,numLDOSpt !Cycle each point in the line
            do imo=1,imoend !Cycle each orbital
                irealmo=imo
                if (ispin==2) irealmo=imo+nbasis
                preterm=str(imo)*0.5D0/pi*FWHM(imo)
                do imappt=1,num2Dpoints !Broaden imo as curve
                    tmp=preterm/( (LDOSxpos(imappt)-MOene_dos(irealmo))**2+0.25D0*FWHM(imo)**2 )
                    LDOS2Dmap(imappt,ilinept)=LDOS2Dmap(imappt,ilinept)+tmp*LDOSptscomp(imo,ilinept)
                end do
            end do
        end do
    end if
    if (ibroadfunc==2.or.ibroadfunc==3) then
        if (ibroadfunc==3) LDOS2Dmap=(1-gauweigh)*LDOS2Dmap
        do ilinept=1,numLDOSpt !Cycle each point in the line
            do imo=1,imoend !Cycle each orbital
                irealmo=imo
                if (ispin==2) irealmo=imo+nbasis
                gauss_c=FWHM(imo)/2D0/sqrt(2*dlog(2D0))
                gauss_a=str(imo)/(gauss_c*sqrt(2D0*pi))
                do imappt=1,num2Dpoints !Broaden imo as curve
                    !Gaussian function, see http://en.wikipedia.org/wiki/Gaussian_function
                    tmp=gauss_a*dexp( -(LDOSxpos(imappt)-MOene_dos(irealmo))**2/(2*gauss_c**2) )
                    if (ibroadfunc==3) tmp=gauweigh*tmp !Combine Lorentizan and Gaussian function
                    LDOS2Dmap(imappt,ilinept)=LDOS2Dmap(imappt,ilinept)+tmp*LDOSptscomp(imo,ilinept)
                end do
            end do
        end do
    end if
    LDOS2Dmap=LDOS2Dmap*scalecurve
    write(*,*)
    write(*,"(' Energy range: from',f12.5,' to',f12.5,1x,a)") enelow,enehigh,trim(unitstr)
    write(*,"(' Starting point:',3f12.6,' Bohr')") orgx,orgy,orgz
    write(*,"(' End point:     ',3f12.6,' Bohr')") endx,endy,endz
    write(*,"(' Stepsize:',f12.6,' Bohr, total length:',f12.6,' Bohr')") dlen,totlen
    
    idraw=1
    isavepic=0
    do while(.true.)
        if (isilent==0.and.idraw==1) then
            lengthx=2300
            if (isavepic==0) then
            else if (isavepic==1) then
            end if
        end if
        
        idraw=0
        isavepic=0
        write(*,*)
        write(*,*) "                    ======== Post-process menu ========"
        write(*,*) "0 Return"
        write(*,*) "1 Show graph again"
        write(*,*) "2 Save the picture to current folder"
        write(*,*) "3 Export curve data to plain text file in current folder"
        write(*,"(a,f8.3)") " 4 Set ratio of Y-axis relative to X-axis, current:",yxratio
        read(*,*) isel2

        if (isel2==0) then
            exit
        else if (isel2==1) then
            idraw=1
            isavepic=0
        else if (isel2==2) then
            idraw=1
            isavepic=1
        else if (isel2==3) then
            open(10,file="LDOS.txt",status="replace")
            do imappt=1,num2Dpoints
                do ipt=1,numLDOSpt
                    write(10,"(f8.3,f10.4,1PE16.8)") LDOSxpos(imappt),dlen*(ipt-1),LDOS2Dmap(imappt,ipt)
                end do
            end do
            close(10)
            write(*,*) "LDOS data have been written to LDOS.txt in current folder"
            write(*,*) "Column 1: Energy ("//trim(unitstr)//")"
            write(*,*) "Column 2: Coordinate in the line (Bohr)"
            write(*,*) "Column 3: LDOS value"
        else if (isel2==4) then
            write(*,*) "Input the ratio, e.g. 1.4"
            read(*,*) yxratio
        end if
    end do
end if




end do !End of main loop
end subroutine
