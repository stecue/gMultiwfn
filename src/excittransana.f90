!!------- Analyze or visualize hole-electron distribution, transition dipole moment and transition density
!ROHF,RODFT are assumed to be impossible to be ground state
!itype=1 is normal mode
!itype=2 is specifically used to calculate delta_r
!itype=3 is specifically used to generate NTOs
!itype=4 is specifically used to Calculate inter-fragment charger transfer
subroutine hetransdipdens(itype)
use defvar
use util
use function
implicit real*8 (a-h,o-z)
integer :: itype,idomag=0
logical,allocatable :: skippair(:)
integer,allocatable :: excdir(:) !excdir=1 means ->, =2 means <-
integer,allocatable :: orbleft(:),orbright(:) !denote the actual MO at the left/right side in the excitation data (1:nbasis=alpha/spatial, nbasis+1:2*nbasis=beta)
real*8,allocatable :: exccoeff(:),exccoeffbackup(:) !Coefficient of an orbital pair transition. exccoeffbackup is used to backup, because users can modify the coefficients
real*8,allocatable :: holegrid(:,:,:),elegrid(:,:,:),holeeleovlp(:,:,:),transdens(:,:,:),holecross(:,:,:),elecross(:,:,:),Cele(:,:,:),Chole(:,:,:),magtrdens(:,:,:,:)
real*8,allocatable :: dipcontri(:,:) !(1/2/3,iexc) contribution of orbital pairs "iexc" to transition dipole moment in X/Y/Z
character c200tmp*200,c80tmp*80,transmodestr*200,leftstr*80,rightstr*80,strtmp1*10,strtmp2*10,strtmp3*10,selectyn
character,save :: excitfilename*200=" "
integer walltime1,walltime2
real*8 time_end,time_endtmp,orbval(nmo),wfnderv(3,nmo)
real*8,allocatable :: tmparr1(:),tmparr2(:) !Arrays for temporary use
real*8,allocatable :: tdmata(:,:),tdmatb(:,:) !Transition density matrix in basis function representation of alpha/total and beta electrons
real*8,allocatable :: bastrdip(:,:) !Contribution from each basis function to transition dipole moment, the first index 1/2/3=X/Y/Z
real*8,allocatable :: trdipmatbas(:,:,:),trdipmatatm(:,:,:) !Transition dipole moment matrix in basis function / atom representation, the first index 1/2/3=X/Y/Z
real*8,allocatable :: bastrpopa(:),bastrpopb(:) !Transition population of each basis function
!Store dipole moment integral between all GTFs, the matrix is symmetry hence stored as 1D-array to save memory. The first index is 1/2/3=X/Y/Z
real*8,allocatable :: GTFdipint(:,:)
!Store information of another excitation
integer,allocatable :: orbleft2(:),orbright2(:),excdir2(:)
real*8,allocatable :: exccoeff2(:)
!Other
real*8,allocatable :: cubx(:),cuby(:),cubz(:)
! real*8 eigvecmat(nbasis,nbasis),eigvalarr(nbasis)

if (.not.allocated(CObasa)) then
    write(*,*) "Error: The input file does not contain basis function information!"
    write(*,*) "Please read manual to make clear which kinds of input file could be used!"
    write(*,*) "Press ENTER to exit"
    read(*,*)
    return
end if

if (excitfilename==" ") then
    write(*,"(a)") " Input the path of the Gaussian/ORCA output file or plain text file containing excitation data, e.g. c:\a.out"
    do while(.true.)
        read(*,"(a)") excitfilename
        inquire(file=excitfilename,exist=alive)
        if (alive) exit
        write(*,*) "Cannot find this file, input again"
    end do
else
    write(*,"(' Loading ',a)") trim(excitfilename)
end if
open(10,file=excitfilename,status="old")
call loclabel(10,"Gaussian",igauout,1,50)
call loclabel(10,"O   R   C   A",iORCAout,1,50)

if (igauout==1) then !Gaussian output file
    write(*,*) "Analyzing the file..."
    call loclabel(10,"Excitation energies and oscillator strengths:")
    read(10,*)
    nstates=0 !The number of transition modes
    do while(.true.)
        call loclabel(10,"Excited State",ifound,0)
        if (ifound==1) then
            nstates=nstates+1
            read(10,*)
        else
            exit
        end if
    end do
    write(*,"(' There are',i5,' transition modes, analyze which one?  e.g. 2')") nstates
    do while(.true.)
        read(*,*) iexcitmode
        if (iexcitmode>0.and.iexcitmode<=nstates) exit
        write(*,*) "Error: The index exceeded valid range, input again"
    end do
    call loclabel(10,"Excitation energies and oscillator strengths:")
    do itmp=1,iexcitmode !Move to the iexcitmode field
        call loclabel(10,"Excited State",ifound,0)
        read(10,*)
    end do
    backspace(10)
    read(10,"(a)") transmodestr
    iexcmulti=1 !Multiplicity of the excited state
    if (index(transmodestr,"Triplet")/=0) iexcmulti=3
    do i=10,70
        if (transmodestr(i:i+1)=="eV") exit
    end do
    read(transmodestr(i-10:i-1),*) excenergy
else if (iORCAout==1) then !ORCA output file
    write(*,*) "Analyzing the file..."
    call loclabel(10,"Number of roots to be determined",ifound)
    if (ifound==0) then
        write(*,*) "Error: It seems that this is not a output file of CIS/TDA/TD task"
        write(*,*) "Press Enter to exit"
        read(*,*)
        return
    end if
    read(10,"(50x,i7)") nstates
    write(*,"(' There are',i5,' transition modes, analyze which one?  e.g. 2')") nstates
    do while(.true.)
        read(*,*) iexcitmode
        if (iexcitmode>0.and.iexcitmode<=nstates) exit
        write(*,*) "Error: The index exceeded valid range, input again"
    end do
    iexcmulti=1 !Multiplicity of the excited state
    if (wfntype==0.or.wfntype==3) then
        call loclabel(10,"Generation of triplets")
        read(10,"(a)") c80tmp
        if (index(c80tmp," on ")/=0) then
            write(*,*) "Load which kind of excited states?"
            write(*,*) "1: Singlet   3: Triplet"
            read(*,*) iexcmulti
        end if
    end if
    write(*,*) "Loading data, please wait..."
    call loclabel(10,"the weight of the individual excitations are printed")
    if (iexcmulti==3) then !When triplets=on, ORCA calculate both singlet and triplet excited state, now move to the latter
        read(10,*)
        call loclabel(10,"the weight of the individual excitations are printed",ifound,0)
    end if
    do itmp=1,iexcitmode !Move to the iexcitmode field
        call loclabel(10,"STATE",ifound,0)
        read(10,*)
    end do
    backspace(10)
    read(10,"(a)") transmodestr
    do i=10,70
        if (transmodestr(i:i+1)=="eV") exit
    end do
    read(transmodestr(i-10:i-1),*) excenergy
else
    rewind(10)
    read(10,*) iexcmulti,excenergy !The first line should be multiplicity of excited state and excitation energy
end if
write(*,"(' Multiplicity of this excited state is',i4)") iexcmulti
write(*,"(' Excitation energy is',f12.7,' eV')") excenergy

!Count how many orbital pairs are involved in this transition mode
nexcitorb=0
do while(.true.)
    read(10,"(a)",iostat=ierror) c80tmp
    if ((index(c80tmp,'<-')==0.and.index(c80tmp,'->')==0).or.ierror/=0) exit
    nexcitorb=nexcitorb+1
end do
allocate(excdir(nexcitorb),orbleft(nexcitorb),orbright(nexcitorb),exccoeff(nexcitorb),exccoeffbackup(nexcitorb))
if (igauout==1.or.iORCAout==1) then !Rewind to head of the entry
    call loclabel(10,trim(transmodestr))
    read(10,*)
else
    rewind(10)
    read(10,*)
end if
write(*,"(a,i8,a)") " There are",nexcitorb," orbital pairs in this transition mode"

!Load transition detail. Notice that for unrestricted case, A and B are separately recorded in input file, &
!however after loading, they are combined as single index, namely if orbital index is larger than nbasis, then it is B, else A
if (iORCAout==1) then
    !Worthnotingly, in at least ORCA 4.0, de-excitation is not separately outputted as <-, but combined into ->
    !Here we still determine <-, because hopefully Neese may change the convention of ORCA output in the future...
    do itmp=1,nexcitorb
        read(10,"(a)") c80tmp
        excdir(itmp)=-1
        if (index(c80tmp,'->')/=0) excdir(itmp)=1 !means ->
        if (index(c80tmp,'<-')/=0) excdir(itmp)=2 !means <-
        if (excdir(itmp)==-1) then
            write(*,"(a)") " Error: This output file does not correspond to a CIS or TD or TDA task, configuration information cannot be found!"
            write(*,*) "Press Enter to exit"
            read(*,*)
            return
        end if
        do isign=1,80 !Find position of <- or ->
            if (c80tmp(isign:isign)=='-'.or.c80tmp(isign:isign)=='<') exit
        end do
        !Process left side of <- or ->
        read(c80tmp(:isign-1),"(a)") leftstr
        read(leftstr(:len_trim(leftstr)-1),*) orbleft(itmp)
        orbleft(itmp)=orbleft(itmp)+1 !ORCA counts orbital from 0 rather than 1!!!
        if (index(leftstr,'b')/=0) orbleft(itmp)=orbleft(itmp)+nbasis
        !Process right side of <- or ->
        read(c80tmp(isign+2:),*) rightstr
        read(rightstr(:len_trim(rightstr)-1),*) orbright(itmp)
        orbright(itmp)=orbright(itmp)+1
        if (index(rightstr,'b')/=0) orbright(itmp)=orbright(itmp)+nbasis
        iTDA=index(c80tmp,'c=')
        if (iTDA/=0) then !CIS, TDA task, configuration coefficients are presented
            read(c80tmp(iTDA+2:iTDA+13),*) exccoeff(itmp)
        else !TD task, configuration coefficients are not presented. Contribution of i->a and i<-a are summed up and outputted as i->a
            if (itmp==1) then
                write(*,"(a)") " Warning: For TD task, ORCA does not print configuration coefficients but only print corresponding contributions of each orbital pair, &
                in this case Multiwfn determines configuration coefficients simply as square root of contribution values. However, this treatment is &
                evidently inappropriate and the result is nonsense when de-excitation is significant (In this situation you have to use TDA-DFT instead)"
                write(*,*) "If you really want to proceed, press ENTER to continue"
                read(*,*)
            end if
            read(c80tmp(23:32),*) tmpval
            if (tmpval<0) excdir(itmp)=2 !Negative contribution is assumed to be de-excitation (of course this is not strict since -> and <- have been combined together)
            exccoeff(itmp)=dsqrt(abs(tmpval))
        end if
        !Although for closed-shell ground state, ORCA still outputs coefficients as normalization to 100%, &
        !However, in order to follow the Gaussian convention, we change the coefficient as normalization to 50%
        if (wfntype==0.or.wfntype==3) exccoeff(itmp)=exccoeff(itmp)/dsqrt(2D0)
    end do
else !Gaussian output file or plain text file
    do itmp=1,nexcitorb
        read(10,"(a)") c80tmp
        excdir(itmp)=1 !means ->
        if (index(c80tmp,'<-')/=0) excdir(itmp)=2 !means <-
        do isign=1,80 !Find position of <- or ->
            if (c80tmp(isign:isign)=='-'.or.c80tmp(isign:isign)=='<') exit
        end do
        !Process left side of <- or ->
        read(c80tmp(:isign-1),"(a)") leftstr
        ilefttype=0 !close
        if (index(leftstr,'A')/=0) ilefttype=1 !Alpha
        if (index(leftstr,'B')/=0) ilefttype=2 !Beta
        if (ilefttype==0) then
            read(leftstr,*) orbleft(itmp)
        else
            read(leftstr(:len_trim(leftstr)-1),*) orbleft(itmp)
            if (ilefttype==2) orbleft(itmp)=orbleft(itmp)+nbasis
        end if
        !Process right side of <- or ->
        read(c80tmp(isign+2:),"(a)") rightstr
        irighttype=0 !close
        if (index(rightstr,'A')/=0) irighttype=1 !Alpha
        if (index(rightstr,'B')/=0) irighttype=2 !Beta
        if (irighttype==0) then
            read(rightstr,*) orbright(itmp),exccoeff(itmp)
        else
            do isplit=1,80
                if (rightstr(isplit:isplit)=='A'.or.rightstr(isplit:isplit)=='B') exit
            end do
            read(rightstr(:isplit-1),*) orbright(itmp)
            read(rightstr(isplit+1:),*) exccoeff(itmp)
            if (irighttype==2) orbright(itmp)=orbright(itmp)+nbasis
        end if
    end do
end if

!Test sum of square of the coefficients
sumsqrexc=0
sumsqrdeexc=0
do itmp=1,nexcitorb
    if (excdir(itmp)==1) sumsqrexc=sumsqrexc+exccoeff(itmp)**2
    if (excdir(itmp)==2) sumsqrdeexc=sumsqrdeexc-exccoeff(itmp)**2
end do
sumsqrall=sumsqrexc+sumsqrdeexc
write(*,"(' The sum of square of excitation coefficients:',f10.6)") sumsqrexc
write(*,"(' The negative of the sum of square of de-excitation coefficients:',f10.6)") sumsqrdeexc
write(*,"(' The sum of above two values',f10.6)") sumsqrall
close(10)
exccoeffbackup=exccoeff


1 do while(.true.)
    if (itype==1) then
        write(*,*)
        if (idomag==0) write(*,*) "-1 Toggle calculating transit. magnetic dip. density in option 1, current: No"
        if (idomag==1) write(*,*) "-1 Toggle calculating transit. magnetic dip. density in option 1, current: Yes"
        write(*,*) "0 Return"
        write(*,"(a)") " 1 Visualize and analyze hole, electron and transition density and so on"
        write(*,*) "2 Show contribution of MO pairs to transition dipole moment"
        write(*,*) "3 Show contribution of each MO to hole and electron distribution"
        write(*,*) "4 Generate and export transition density matrix"
        write(*,*) "5 Decompose transition dipole moment to basis function and atom contributions"
        write(*,*) "6 Calculate Mulliken atomic transition charges"
        !Note: 7 is not utilized
        write(*,*) "10 Modify or check excitation coefficients"
        read(*,*) isel
    else if (itype==2) then
        call calcdelta_r(nexcitorb,orbleft,orbright,excdir,exccoeff)
        return
    else if (itype==3) then
        call NTO(nexcitorb,orbleft,orbright,excdir,exccoeff)
        return
    else if (itype==4) then
        call excfragCT(nexcitorb,orbleft,orbright,excdir,exccoeff)
        return
    end if
    
    if (isel==-1) then
        if (idomag==0) then
            idomag=1
        else if (idomag==1) then
            idomag=0
        end if
    else if (isel==0) then
        return
    else if (isel==1) then
        exit
    !Show contribution of MO pairs to transition dipole moment
    else if (isel==2) then
        write(*,"(a)") " Input the threshold for printing, e.g. 0.01 means the orbital pairs having contribution &
        to any component of transition dipole moment >= 0.01 will be shown"
        read(*,*) printthres
        write(*,"(a)") " Input the threshold for calculating contribution, e.g. 0.005 means the configurations with &
        absolute value of coefficient smaller than 0.005 will be ignored. The smaller the value, the higher the computational cost"
        read(*,*) calcthres
        write(*,*) "Generating dipole moment integral matrix..."
        nsize=nprims*(nprims+1)/2
        allocate(GTFdipint(3,nsize))
        call genGTFDmat(GTFdipint,nsize)
        allocate(dipcontri(3,nexcitorb))
        dipcontri=0
        fac=1
        if (wfntype==0.or.wfntype==3) fac=2
        write(*,*) "Calculating contribution to transition dipole moment, please wait..."
nthreads=getNThreads()
!$OMP PARALLEL DO SHARED(dipcontri) PRIVATE(iGTF,jGTF,ides,iexcitorb,imo,jmo) schedule(dynamic) NUM_THREADS(nthreads)
        do iexcitorb=1,nexcitorb
            if (abs(exccoeff(iexcitorb))<calcthres) cycle
!             write(*,*) iexcitorb,nexcitorb
            imo=orbleft(iexcitorb)
            jmo=orbright(iexcitorb)
            do iGTF=1,nprims
                do jGTF=1,nprims
                    if (iGTF>=jGTF) then
                        ides=iGTF*(iGTF-1)/2+jGTF
                    else
                        ides=jGTF*(jGTF-1)/2+iGTF
                    end if
                    dipcontri(:,iexcitorb)=dipcontri(:,iexcitorb)+co(imo,iGTF)*co(jmo,jGTF)*GTFdipint(:,ides)
                end do
            end do
            dipcontri(:,iexcitorb)=dipcontri(:,iexcitorb)*exccoeff(iexcitorb)*fac
        end do
!$OMP END PARALLEL DO
        deallocate(GTFdipint)
        xdipsum=sum(dipcontri(1,:))
        ydipsum=sum(dipcontri(2,:))
        zdipsum=sum(dipcontri(3,:))
        ishownpair=0
        write(*,*) "   #Pair                  Coefficient   Transition dipole X/Y/Z (au)"
        do iexcitorb=1,nexcitorb
            if ( any(abs(dipcontri(:,iexcitorb))>printthres) ) then
                ishownpair=ishownpair+1
                imo=orbleft(iexcitorb)
                jmo=orbright(iexcitorb)
                strtmp1=" ->"
                if (excdir(iexcitorb)==2) strtmp1=" <-"
                if (wfntype==0.or.wfntype==3) then
                    write(*,"(i8,i7,a,i7,f12.6,3f11.6)") iexcitorb,imo,trim(strtmp1),jmo,exccoeff(iexcitorb),dipcontri(:,iexcitorb)
                else
                    strtmp2="A"
                    strtmp3="A"
                    if (imo>nbasis) then
                        imo=imo-nbasis
                        strtmp2="B"
                    end if
                    if (jmo>nbasis) then
                        jmo=jmo-nbasis
                        strtmp3="B"
                    end if
                    write(*,"(i8,i6,a,a,i6,a,f12.6,3f11.6)") iexcitorb,imo,trim(strtmp2),trim(strtmp1),jmo,trim(strtmp3),exccoeff(iexcitorb),&
                    dipcontri(:,iexcitorb)
                end if
            end if
        end do
        if (printthres==0) then
            write(*,"(' Sum:                                ',3f11.6)") xdipsum,ydipsum,zdipsum
        else
            write(*,"(' Sum (including the ones not shown): ',3f11.6)") xdipsum,ydipsum,zdipsum
!             write(*,"(i8,' orbital pairs are shown above')") ishownpair
        end if
        oscillstr=2D0/3D0*excenergy/au2eV*(xdipsum**2+ydipsum**2+zdipsum**2)
        write(*,"(' Norm of total transition dipole moment:',f11.6,' a.u.')") dsqrt(xdipsum**2+ydipsum**2+zdipsum**2)
        write(*,"(' Oscillator strength:',f12.7)") oscillstr
        if ((naelec==nbelec).and.iexcmulti==3) write(*,"(a)") " Note: Since the spin multiplicity between the ground state and excited state is different, &
        the transition dipole moment and thus oscillator strength should be exactly zero in principle"
        
        write(*,*)
        write(*,*) "If output all orbital pairs to transdip.txt in current folder? (y/n)"
        read(*,*) c80tmp
        if (c80tmp(1:1)=='y'.or.c80tmp(1:1)=='Y') then
            open(10,file="transdip.txt",status="replace")
            write(10,*) "   #Pair                  Coefficient   Transition dipole X/Y/Z (au)"
            do iexcitorb=1,nexcitorb
                imo=orbleft(iexcitorb)
                jmo=orbright(iexcitorb)
                strtmp1=" ->"
                if (excdir(iexcitorb)==2) strtmp1=" <-"
                if (wfntype==0.or.wfntype==3) then
                    write(10,"(i8,i7,a,i7,f12.6,3f11.6)") iexcitorb,imo,trim(strtmp1),jmo,exccoeff(iexcitorb),dipcontri(:,iexcitorb)
                else
                    strtmp2="A"
                    strtmp3="A"
                    if (imo>nbasis) then
                        imo=imo-nbasis
                        strtmp2="B"
                    end if
                    if (jmo>nbasis) then
                        jmo=jmo-nbasis
                        strtmp3="B"
                    end if
                    write(10,"(i8,i6,a,a,i6,a,f12.6,3f11.6)") iexcitorb,imo,trim(strtmp2),trim(strtmp1),jmo,trim(strtmp3),exccoeff(iexcitorb),dipcontri(:,iexcitorb)
                end if
            end do
            write(10,"(' Sum:                                ',3f11.6)") xdipsum,ydipsum,zdipsum
            write(10,"(' Norm of transition dipole moment:',f12.7,' a.u.')") dsqrt(xdipsum**2+ydipsum**2+zdipsum**2)
            write(10,"(' Oscillator strength:',f12.7)") oscillstr
            close(10)
            write(*,*) "Done! Output finished"
        end if
        deallocate(dipcontri)
    
    !Show contribution of each MO to hole and electron distribution
    else if (isel==3) then
        write(*,*) "Input printing criterion"
        write(*,*) "e.g. 0.005 means only the MOs having contribution >=0.005 will be printed"
        read(*,*) printcrit
        allocate(tmparr1(nmo),tmparr2(nmo))
        tmparr1=0 !Record contribution of each MO to hole
        tmparr2=0 !Record contribution of each MO to electron
        do iexcitorb=1,nexcitorb !Cycle each excitation pair
            imo=orbleft(iexcitorb)
            jmo=orbright(iexcitorb)
            excwei=exccoeff(iexcitorb)**2
            if (excdir(iexcitorb)==1) then ! ->
                tmparr1(imo)=tmparr1(imo)+excwei
                tmparr2(jmo)=tmparr2(jmo)+excwei
            else ! <-
                tmparr1(imo)=tmparr1(imo)-excwei
                tmparr2(jmo)=tmparr2(jmo)-excwei
            end if
        end do
        if (wfntype==0.or.wfntype==3) then !Ground state is close-shell
            tmparr1=tmparr1*2
            tmparr2=tmparr2*2
            do imo=1,nmo
                if (tmparr1(imo)>=printcrit.or.tmparr2(imo)>=printcrit) &
                write(*,"(' MO',i8,', Occ:',f10.5,'    Hole:',f9.5,'     Electron:',f9.5)") imo,MOocc(imo),tmparr1(imo),tmparr2(imo)
            end do
        else !Ground state is open-shell
            do imo=1,nmo
                if (tmparr1(imo)>=printcrit.or.tmparr2(imo)>=printcrit) then
                    if (imo<=nbasis) then
                        write(*,"(' MO',i7,'A, Occ:',f10.5,'    Hole:',f9.5,'     Electron:',f9.5)") imo,MOocc(imo),tmparr1(imo),tmparr2(imo)
                    else
                        write(*,"(' MO',i7,'B, Occ:',f10.5,'    Hole:',f9.5,'     Electron:',f9.5)") imo-nbasis,MOocc(imo),tmparr1(imo),tmparr2(imo)
                    end if
                end if
            end do
        end if
        write(*,"(' Sum of hole:',f9.5,'    Sum of electron:',f9.5)") sum(tmparr1),sum(tmparr2)
        deallocate(tmparr1,tmparr2)
    
    !==4: Generating transition density matrix (TDM)
    !==5: Decompose transition electric/magnetic dipole moment as basis function and atom contributions
    !==6: Calculate Mulliken atomic transition charges
    else if (isel==4.or.isel==5.or.isel==6) then
        if (isel==5) then
            write(*,*) "Decompose which type of transition dipole moment?"
            write(*,*) "1: Electric      2: Magnetic"
            read(*,*) idecomptype
        end if
        !There are two ways to construct TDM for TD method, see eqs. 22~24 in JCP,66,3460 for detail.
        !The first one is correct for transition electric dipole moment, excitation and de-excitation are not distinguished
        if (isel==4.or.(isel==5.and.idecomptype==1).or.isel==6) iTDMtype=1
        !The second one is correct for transition velocity/magnetic dipole moment, excitation and de-excitation are considered individually 
        if (isel==5.and.idecomptype==2) iTDMtype=2
        write(*,*) "Generating transition density matrix..."
        if (.not.allocated(tdmata)) allocate(tdmata(nbasis,nbasis))
        if ((wfntype==1.or.wfntype==4).and.(.not.allocated(tdmatb))) allocate(tdmatb(nbasis,nbasis))
        iprog=0
nthreads=getNThreads()
!$OMP parallel do shared(tdmata,tdmatb,iprog) private(ibas,jbas,iexcitorb,imo,jmo,tmpval,tmpval2) num_threads(nthreads) SCHEDULE(DYNAMIC)
        do ibas=1,nbasis
            do jbas=1,nbasis
                tmpval=0
                tmpval2=0
                if (iTDMtype==1) then
                    do iexcitorb=1,nexcitorb
                        imo=orbleft(iexcitorb)
                        jmo=orbright(iexcitorb)
                        if (MOtype(imo)==0.or.MOtype(imo)==1) then !Close-shell or alpha part of open-shell
                            tmpval=tmpval+exccoeff(iexcitorb)*cobasa(ibas,imo)*cobasa(jbas,jmo)
                        else !Beta part of open-shell
                            tmpval2=tmpval2+exccoeff(iexcitorb)*cobasb(ibas,imo-nbasis)*cobasb(jbas,jmo-nbasis)
                        end if
                    end do
                else if (iTDMtype==2) then
                    do iexcitorb=1,nexcitorb
                        imo=orbleft(iexcitorb)
                        jmo=orbright(iexcitorb)
                        if (MOtype(imo)==0.or.MOtype(imo)==1) then !Close-shell or alpha part of open-shell
                            if (excdir(iexcitorb)==1) tmpval=tmpval+exccoeff(iexcitorb)*cobasa(ibas,imo)*cobasa(jbas,jmo)
                            if (excdir(iexcitorb)==2) tmpval=tmpval-exccoeff(iexcitorb)*cobasa(ibas,imo)*cobasa(jbas,jmo)
                        else !Beta part of open-shell
                            if (excdir(iexcitorb)==1) tmpval2=tmpval2+exccoeff(iexcitorb)*cobasb(ibas,imo-nbasis)*cobasb(jbas,jmo-nbasis)
                            if (excdir(iexcitorb)==2) tmpval2=tmpval2-exccoeff(iexcitorb)*cobasb(ibas,imo-nbasis)*cobasb(jbas,jmo-nbasis)
                        end if
                    end do
                end if
                tdmata(ibas,jbas)=tmpval
                if (wfntype==1.or.wfntype==4) tdmatb(ibas,jbas)=tmpval2
            end do
!$OMP CRITICAL
            iprog=iprog+1
            if (nbasis>800) write(*,"(' Progress:',i6,'  /',i6)") iprog,nbasis
!$OMP END CRITICAL
        end do
!$OMP END PARALLEL do
        
        if (wfntype==0.or.wfntype==3) tdmata=tdmata*2 !Close-shell, should double the TDM
        
        !! Below codes are used to check transition properties based on transition density matrix and corresponding integral matrix
        !BEWARE THAT CARTESIAN BASIS FUNCTIONS MUST BE USED, SO DON'T FORGET 6D 10F KEYWORDS IN DUE CASES! (However, pure type is OK if you have
        !set igenDbas/igenMagbas in settings.ini, because the Cartesian integral matrix generated when loading has already been converted to pure case)
        !Check transition eletric dipole moment. The result is correct when iTDMtype==1, see above
!         if (.not.allocated(Dbas)) then
!             allocate(Dbas(3,nbasis,nbasis))
!             call genDbas
!         end if
!         Teledipx=sum(tdmata*Dbas(1,:,:))
!         Teledipy=sum(tdmata*Dbas(2,:,:))
!         Teledipz=sum(tdmata*Dbas(3,:,:))
!         write(*,"(' Transition electric dipole moment:',3f12.6)") Teledipx,Teledipy,Teledipz
!         !Check transition velocity dipole moment. The result is correct when iTDMtype==2, see above
!         if (.not.allocated(Velbas)) then
!             allocate(Velbas(3,nbasis,nbasis))
!             call genvelbas
!         end if
!         Tvdipx=sum(tdmata*Velbas(1,:,:))
!         Tvdipy=sum(tdmata*Velbas(2,:,:))
!         Tvdipz=sum(tdmata*Velbas(3,:,:))
!         write(*,"(' Transition velocity dipole moment:',3f12.6)") Tvdipx,Tvdipy,Tvdipz
!         !Check transition magnetic dipole moment. The result is correct when iTDMtype==2, see above
!         if (.not.allocated(Magbas)) then
!             allocate(Magbas(3,nbasis,nbasis))
!             call genMagbas
!         end if
!         call showmatgau(tdmata,"tdmata matrix",1)
!         Tmagdipx=sum(tdmata*Magbas(1,:,:))
!         Tmagdipy=sum(tdmata*Magbas(2,:,:))
!         Tmagdipz=sum(tdmata*Magbas(3,:,:))
!         write(*,"(' Transition magnetic dipole moment:',3f12.6)") Tmagdipx,Tmagdipy,Tmagdipz
!         !Below formula is absolutely correct, however the printed result is not correct here because 
!         !we cannot obtain correct transition electric and magnetic dipole moments based on the same TDM.
!         !If we directly use the value outputted by Gaussian, you will find the result is in line with the rotatory strength printed by Gaussian
!         !Four points:
!         !1) The negative sign: Because we ignored the negative sign when evaluating magnetic integrals
!         !2) Diveded by two: Necessary by definition
!         !3) 2.54174619D-018: convert electric dipole moment from a.u. to cgs, see gabedit build-in converter
!         !4) 1.85480184D-020: convert magnetic dipole moment from a.u. to cgs, see gabedit build-in converter
!         Rlen=-(Teledipx*Tmagdipx+Teledipy*Tmagdipy+Teledipz*Tmagdipz)/2D0*2.54174619D-018*1.85480184D-020 *1D40
!         write(*,"(' Rotatory strength in length representation:',f14.8,' 10^-40 cgs')") Rlen
!         cycle
        
        if (isel==4) then !Export transition density matrix
            write(*,*) "If symmetrizing the transition density matrix? (y/n)"
            read(*,*) c80tmp
            if (c80tmp(1:1)=='y'.or.c80tmp(1:1)=='Y') then
                do ibas=1,nbasis
                    do jbas=ibas,nbasis
                        tdmata(ibas,jbas)=(tdmata(ibas,jbas)+tdmata(jbas,ibas))/dsqrt(2D0)
                        tdmata(jbas,ibas)=tdmata(ibas,jbas)
                        if (wfntype==1.or.wfntype==4) then
                            tdmatb(ibas,jbas)=(tdmatb(ibas,jbas)+tdmatb(jbas,ibas))/dsqrt(2D0)
                            tdmatb(jbas,ibas)=tdmatb(ibas,jbas)
                        end if
                    end do
                end do
            end if
            !Diagonalize transition density matrix, only for test purpose
!             call diagsymat(tdmata,eigvecmat,eigvalarr,istat)
!             call diaggemat(tdmata,eigvecmat,eigvalarr,istat)
!             write(*,"(f12.6)") eigvalarr(:)
!             write(*,*) sum(eigvalarr),sum(eigvalarr,eigvalarr>0),sum(eigvalarr,eigvalarr<0)
            if (wfntype==0.or.wfntype==3) then
                open(10,file="tdmat.txt",status="replace")
                call showmatgau(tdmata,"Transition density matrix",0,"f14.8",10)
                close(10)
                write(*,"(a)") " Done! The transition density matrix has been outputted to tdmat.txt in current folder"
                write(*,"(a)") " Hint: You can use function 2 of main function 18 to plot the transition density matrix by using this file as input file"
            else
                open(10,file="tdmata.txt",status="replace")
                call showmatgau(tdmata,"Alpha transition density matrix",0,"f14.8",10)
                close(10)
                open(10,file="tdmatb.txt",status="replace")
                call showmatgau(tdmatb,"Beta transition density matrix",0,"f14.8",10)
                close(10)
                write(*,"(a)") " Done! The alpha and beta transition density matrix has been outputted to tdmata.txt and tdmatb.txt in current folder, respectively"
                write(*,"(a)") " Hint: You can use function 2 of main function 18 to plot the transition density matrix by using these files as input file"
            end if
        
        else if (isel==5) then !Decompose transition electric/magnetic dipole moment as basis function and atom contributions
            if (idecomptype==1.and.(.not.allocated(Dbas))) then
                write(*,"(a)") " ERROR: Electric dipole moment integral matrix is not presented. You should set parameter ""igenDbas"" in settings.ini to 1, &
                so that the matrix will be generated when loading wavefunction file"
            else if (idecomptype==2.and.(.not.allocated(Magbas))) then
                write(*,"(a)") " ERROR: Magnetic dipole moment integral matrix is not presented. You should set parameter ""igenMagbas"" in settings.ini to 1, &
                so that the matrix will be generated when loading wavefunction file"
            else
                allocate(trdipmatbas(3,nbasis,nbasis),bastrdip(3,nbasis),trdipmatatm(3,ncenter,ncenter))
                if (idecomptype==1) then
                    if (wfntype==1.or.wfntype==4) then
                        write(*,*) "Analyze alpha or beta part?   1=Alpha  2=Beta"
                        read(*,*) iseltmp
                        if (iseltmp==1) then
                            trdipmatbas(1,:,:)=tdmata(:,:)*Dbas(1,:,:)
                            trdipmatbas(2,:,:)=tdmata(:,:)*Dbas(2,:,:)
                            trdipmatbas(3,:,:)=tdmata(:,:)*Dbas(3,:,:)
                        else
                            trdipmatbas(1,:,:)=tdmatb(:,:)*Dbas(1,:,:)
                            trdipmatbas(2,:,:)=tdmatb(:,:)*Dbas(2,:,:)
                            trdipmatbas(3,:,:)=tdmatb(:,:)*Dbas(3,:,:)
                        end if
                    else !Close-shell
                        trdipmatbas(1,:,:)=tdmata(:,:)*Dbas(1,:,:)
                        trdipmatbas(2,:,:)=tdmata(:,:)*Dbas(2,:,:)
                        trdipmatbas(3,:,:)=tdmata(:,:)*Dbas(3,:,:)
                    end if
                else if (idecomptype==2) then
                    if (wfntype==1.or.wfntype==4) then
                        write(*,*) "Analyze alpha or beta part?   1=Alpha  2=Beta"
                        read(*,*) iseltmp
                        if (iseltmp==1) then
                            trdipmatbas(1,:,:)=tdmata(:,:)*Magbas(1,:,:)
                            trdipmatbas(2,:,:)=tdmata(:,:)*Magbas(2,:,:)
                            trdipmatbas(3,:,:)=tdmata(:,:)*Magbas(3,:,:)
                        else
                            trdipmatbas(1,:,:)=tdmatb(:,:)*Magbas(1,:,:)
                            trdipmatbas(2,:,:)=tdmatb(:,:)*Magbas(2,:,:)
                            trdipmatbas(3,:,:)=tdmatb(:,:)*Magbas(3,:,:)
                        end if
                    else !Close-shell
                        trdipmatbas(1,:,:)=tdmata(:,:)*Magbas(1,:,:)
                        trdipmatbas(2,:,:)=tdmata(:,:)*Magbas(2,:,:)
                        trdipmatbas(3,:,:)=tdmata(:,:)*Magbas(3,:,:)
                    end if
                end if
                open(10,file="trdipcontri.txt",status="replace")
                ides=10
                write(ides,"(a)") " Contribution from basis functions to transition dipole moment (X,Y,Z, in a.u.) derived by Mulliken partition:"
                bastrdip=0
                do ibas=1,nbasis !Note that trdipmatbas is not strictly a symmetry matrix, so we get average value for off-diagonal elements
                    do jbas=1,nbasis
                        if (ibas==jbas) then
                            bastrdip(:,ibas)=bastrdip(:,ibas)+trdipmatbas(:,ibas,jbas)
                        else
                            bastrdip(:,ibas)=bastrdip(:,ibas)+(trdipmatbas(:,ibas,jbas)+trdipmatbas(:,jbas,ibas))/2
                        end if
                    end do
                    write(ides,"(i5,'#  Shell:',i5,' Cen:',i5,'(',a2,') Type:',a,3f11.5)") &
                    ibas,basshell(ibas),bascen(ibas),a(bascen(ibas))%name,GTFtype2name(bastype(ibas)),bastrdip(:,ibas)
                end do
                write(ides,*)
                write(ides,"(a)") " Contribution from atoms to transition dipole moment (X,Y,Z, in a.u.) derived by Mulliken partition:"
                do iatm=1,ncenter
                    write(ides,"(i5,'(',a2,') :',3f12.5)") iatm,a(iatm)%name,&
                    sum(bastrdip(1,basstart(iatm):basend(iatm))),sum(bastrdip(2,basstart(iatm):basend(iatm))),sum(bastrdip(3,basstart(iatm):basend(iatm)))
                end do
                write(ides,*)
                write(ides,"(a,3f12.6,a)") " Transition dipole moment in X/Y/Z",sum(bastrdip(1,:)),sum(bastrdip(2,:)),sum(bastrdip(3,:))," a.u."
                close(10)
                write(*,*) "Done! The result has been outputted to trdipcontri.txt in current folder"
                
                write(*,*)
                write(*,"(a)") " Would you also like to output atom-atom contribution matrix to AAtrdip.txt in current folder? (y/n)"
                read(*,*) selectyn
                if (selectyn=='y') then
                    do iatm=1,ncenter
                        do jatm=1,ncenter
                            trdipmatatm(1,iatm,jatm)=sum( trdipmatbas(1,basstart(iatm):basend(iatm),basstart(jatm):basend(jatm)) )
                            trdipmatatm(2,iatm,jatm)=sum( trdipmatbas(2,basstart(iatm):basend(iatm),basstart(jatm):basend(jatm)) )
                            trdipmatatm(3,iatm,jatm)=sum( trdipmatbas(3,basstart(iatm):basend(iatm),basstart(jatm):basend(jatm)) )
                            !Use the first slot to record sum of square
                            trdipmatatm(1,iatm,jatm)=trdipmatatm(1,iatm,jatm)**2+trdipmatatm(2,iatm,jatm)**2+trdipmatatm(3,iatm,jatm)**2
                        end do
                    end do
                    open(10,file="AAtrdip.txt",status="replace")
                    call showmatgau(trdipmatatm(1,:,:),"Atom-Atom contribution matrix",0,"f14.8",10)
                    close(10)
                    write(*,"(a)") " Done! The matrix has been outputted to AAtrdip.txt in current folder"
                    write(*,"(a)") " This file can be plotted as colored matrix via subfunction 2 of main function 18"
                end if
                deallocate(trdipmatbas,bastrdip,trdipmatatm)
            end if
            
        else if (isel==6) then !Calculate Mulliken atomic transition charges
            allocate(bastrpopa(nbasis))
            bastrpopa=0
            if (wfntype==1.or.wfntype==4) then
                allocate(bastrpopb(nbasis))
                bastrpopb=0
                open(10,file="atmtrchga.chg",status="replace")
                open(11,file="atmtrchgb.chg",status="replace")
            else
                open(10,file="atmtrchg.chg",status="replace")
            end if
            do ibas=1,nbasis !Note that tdmat is not strictly a symmetry matrix, so we get average value for off-diagonal elements
                do jbas=1,nbasis
                    if (ibas==jbas) then
                        bastrpopa(ibas)=bastrpopa(ibas)+tdmata(ibas,jbas)*Sbas(ibas,jbas)
                        if (wfntype==1.or.wfntype==4) bastrpopb(ibas)=bastrpopb(ibas)+tdmatb(ibas,jbas)*Sbas(ibas,jbas)
                    else
                        bastrpopa(ibas)=bastrpopa(ibas)+(tdmata(ibas,jbas)+tdmata(jbas,ibas))*Sbas(ibas,jbas)/2
                        if (wfntype==1.or.wfntype==4) bastrpopb(ibas)=bastrpopb(ibas)+(tdmatb(ibas,jbas)+tdmatb(jbas,ibas))*Sbas(ibas,jbas)/2
                    end if
                end do
            end do
            do iatm=1,ncenter
                write(10,"(1x,a,4f12.6)") a(iatm)%name,a(iatm)%x*b2a,a(iatm)%y*b2a,a(iatm)%z*b2a,sum(bastrpopa(basstart(iatm):basend(iatm)))
                if (wfntype==1.or.wfntype==4) write(11,"(1x,a,4f12.6)") a(iatm)%name,a(iatm)%x*b2a,a(iatm)%y*b2a,a(iatm)%z*b2a,sum(bastrpopb(basstart(iatm):basend(iatm)))
            end do
            close(10)
            if (wfntype==1.or.wfntype==4) then
                close(11)
                write(*,"(a)") " Done! Mulliken atomic transition charges for alpha and beta electrons have been outputted to atmtrchga.chg and atmtrchgb.chg in current folder, respectively"
            else
                write(*,"(a)") " Done! Mulliken atomic transition charges have been outputted to atmtrchg.chg in current folder"
            end if
            deallocate(bastrpopa)
            if (wfntype==1.or.wfntype==4) deallocate(bastrpopb)
        end if
        deallocate(tdmata)
        if (wfntype==1.or.wfntype==4) deallocate(tdmatb)
        
    !Modify or check excitation coefficients    
    else if (isel==10) then
        write(*,"(' Note: There are',i7,' orbital pairs')") nexcitorb
        do while(.true.)
            sumsqrexc=0
            sumsqrdeexc=0
            do itmp=1,nexcitorb
                if (excdir(itmp)==1) sumsqrexc=sumsqrexc+exccoeff(itmp)**2
                if (excdir(itmp)==2) sumsqrdeexc=sumsqrdeexc-exccoeff(itmp)**2
            end do
            write(*,*)
            write(*,"(' Sum of current coeff.^2 of excitation and de-excitation:',2f11.6)") sumsqrexc,sumsqrdeexc
            nzerocoeff=count(exccoeff==0)
            if (nzerocoeff>0) write(*,"(' The number of pairs having zero coefficient:',i7)") nzerocoeff
!             if (wfntype==0.or.wfntype==3) then
!                 write(*,"(' The MO index ranges from',i7,' to',i7,', the first',i7,' are occupied')") 1,nmo,nint(naelec)
!             else
!                 write(*,"(' The alpha MO index ranges from',i7,' to',i7,', the first',i7,' are occupied')") 1,nbasis,nint(naelec)
!                 write(*,"(' The beta MO index ranges from ',i7,' to',i7,', the first',i7,' are occupied')") 1,nbasis,nint(nbelec)
!             end if
            write(*,*)
            write(*,*) "-2 Print coefficient (and contribution) of some orbital pairs"
            write(*,*) "-1 Reset coefficient of all orbital pairs"
            write(*,*) "0 Return"
            write(*,*) "1 Set coefficient of an orbital pair"
            write(*,*) "2 Set coefficient for specific range of orbital pairs"
            write(*,*) "3 Clean all coefficients but for a specific orbital pair"
            read(*,*) isel
            if (isel==-2) then
                write(*,*) "Input the threshold for printing, e.g. 0.01"
                write(*,"(a)") " Note: The orbital pairs whose absolute value of coefficient >= this value will be printed. &
                If input 0, all pairs will be printed."
                read(*,*) printthres
                ntmp=0
                do iexcitorb=1,nexcitorb
                    if (abs(exccoeff(iexcitorb))<printthres) cycle
                    imo=orbleft(iexcitorb)
                    jmo=orbright(iexcitorb)
                    strtmp1=" ->"
                    contrisign=1
                    if (excdir(iexcitorb)==2) then
                        strtmp1=" <-"
                        contrisign=-1
                    end if
                    if (wfntype==0.or.wfntype==3) then
                        write(*,"(i8,i7,a,i7,'   Coeff.:',f10.4,'   Contri.:',f10.4,'%')") iexcitorb,imo,trim(strtmp1),&
                        jmo,exccoeff(iexcitorb),contrisign*exccoeff(iexcitorb)**2/0.5D0*100
                    else
                        strtmp2="A"
                        strtmp3="A"
                        if (imo>nbasis) then
                            imo=imo-nbasis
                            strtmp2="B"
                        end if
                        if (jmo>nbasis) then
                            jmo=jmo-nbasis
                            strtmp3="B"
                        end if
                        write(*,"(i8,i6,a,a,i6,a,'   Coeff.:',f10.4,'   Contri.:',f10.4,'%')") iexcitorb,imo,trim(strtmp2),trim(strtmp1),&
                        jmo,trim(strtmp3),exccoeff(iexcitorb),contrisign*exccoeff(iexcitorb)**2*100
                    end if
                    ntmp=ntmp+1
                end do
                write(*,"(' There are',i8,' orbital pairs met the criterion you defined')") ntmp
            else if (isel==-1) then
                exccoeff=exccoeffbackup
                write(*,*) "All coefficents have been reset"
            else if (isel==0) then
                exit
            else if (isel==1.or.isel==3) then
                if (wfntype==0.or.wfntype==3) then
                    write(*,*) "Input the index of the MOs orginally occupied and unoccupied, e.g. 12,26"
                    read(*,*) iMO,jMO
                else
                    write(*,*) "Input the index of the MOs orginally occupied and unoccupied"
                    write(*,*) "e.g. 12A,46A    or    23B,39B"
                    read(*,*) strtmp1,strtmp2
                    read(strtmp1(1:len_trim(strtmp1)-1),*) iMO
                    read(strtmp2(1:len_trim(strtmp2)-1),*) jMO
                    if (index(strtmp1,'B')/=0) then
                        iMO=iMO+nbasis
                        jMO=jMO+nbasis
                    end if
                end if
                idir=1
                if (any(excdir==2)) then
                    write(*,*) "Which is the type of the orbital pair?"
                    if (isel==1) write(*,*) "1: Excitation (namely ->)  2: De-excitation (namely <-)"
                    if (isel==3) write(*,*) "1: Excitation (namely ->)  2: De-excitation (namely <-)  3: Both"
                    read(*,*) idir
                end if
                if (isel==1) then !Set coefficient for an orbital pair
                    do iexcitorb=1,nexcitorb
                        if (orbleft(iexcitorb)==iMO.and.orbright(iexcitorb)==jMO.and.excdir(iexcitorb)==idir) then
                            ipair=iexcitorb
                            exit
                        end if
                    end do
                    if (iexcitorb==nexcitorb+1) then
                        write(*,*) "Warning: Cannot find corresponding orbital pair!"
                    else
                        write(*,"(' Original value is',f11.7,', now input the new value, e.g. 0.0143')") exccoeff(ipair)
                        read(*,*) exccoeff(ipair)
                    end if
                else if (isel==3) then !Clean all coefficients but conserve one
                    ipair1=0
                    ipair2=0
                    do iexcitorb=1,nexcitorb
                        if (orbleft(iexcitorb)==iMO.and.orbright(iexcitorb)==jMO) then
                            if ( excdir(iexcitorb)==1.and.(idir==1.or.idir==3) ) then
                                ipair1=iexcitorb
                                tmp1=exccoeff(ipair1)
                            else if ( excdir(iexcitorb)==2.and.(idir==2.or.idir==3) ) then
                                ipair2=iexcitorb
                                tmp2=exccoeff(ipair2)
                            end if
                        end if
                    end do
                    if ((idir==1.and.ipair1==0).or.(idir2==2.and.ipair2==0).or.(idir==3.and.ipair1==0.and.ipair2==0)) then
                        write(*,*) "Error: Cannot find corresponding orbital pair!"
                    else
                        exccoeff=0
                        if (idir==1.or.(idir==3.and.ipair1/=0)) exccoeff(ipair1)=tmp1
                        if (idir==2.or.(idir==3.and.ipair2/=0)) exccoeff(ipair2)=tmp2
                        write(*,*) "Done!"
                    end if
                end if
            else if (isel==2) then
                if (wfntype==0.or.wfntype==3) then
                    write(*,*) "Input the index range of the MOs orginally occupied, e.g. 14,17"
                else
                    write(*,*) "Input the index range of the MOs orginally occupied"
                    write(*,*) "e.g. 12A,14A    or    22B,28B"
                end if
                write(*,*) "Note: Simply input 0 can select all orginally occupied MOs"
                read(*,"(a)") c80tmp
                if (c80tmp(1:1)=='0') then
                    istart=1
                    iend=nmo
                else
                    if (wfntype==0.or.wfntype==3) then
                        read(c80tmp,*) istart,iend
                    else
                        read(c80tmp,*) strtmp1,strtmp2
                        read(strtmp1(1:len_trim(strtmp1)-1),*) istart
                        read(strtmp2(1:len_trim(strtmp2)-1),*) iend
                        if (index(strtmp1,'B')/=0) then
                            istart=istart+nbasis
                            iend=iend+nbasis
                        end if
                    end if
                end if
                
                if (wfntype==0.or.wfntype==3) then
                    write(*,*) "Input the index range of the MOs orginally unoccupied, e.g. 72,93"
                else
                    write(*,*) "Input the index range of the MOs orginally unoccupied"
                    write(*,*) "e.g. 72A,84A    or    62B,98B"
                end if
                write(*,*) "Note: Simply input 0 can select all orginally unoccupied MOs"
                read(*,"(a)") c80tmp
                if (c80tmp(1:1)=='0') then
                    jstart=1
                    jend=nmo
                else
                    if (wfntype==0.or.wfntype==3) then
                        read(c80tmp,*) jstart,jend
                    else
                        read(c80tmp,*) strtmp1,strtmp2
                        read(strtmp1(1:len_trim(strtmp1)-1),*) jstart
                        read(strtmp2(1:len_trim(strtmp2)-1),*) jend
                        if (index(strtmp1,'B')/=0) then
                            jstart=jstart+nbasis
                            jend=jend+nbasis
                        end if
                    end if
                end if
                idir=1
                if (any(excdir==2)) then
                    write(*,*) "Which is the type of these orbital pairs?"
                    write(*,*) "1: Excitation (namely ->)  2: De-excitation (namely <-)  3: Any"
                    read(*,*) idir
                end if
                write(*,*) "Set the coefficients to which value?  e.g. 0.00148"
                read(*,*) tmpval
                iset=0
                do iexcitorb=1,nexcitorb
                    if (orbleft(iexcitorb)>=istart.and.orbleft(iexcitorb)<=iend.and.orbright(iexcitorb)>=jstart.and.orbright(iexcitorb)<=jend) then
                        if (idir==3.or.excdir(iexcitorb)==idir) then
                            exccoeff(iexcitorb)=tmpval
                            iset=iset+1
                        end if
                    end if
                end do
                write(*,"(' Coefficient of',i7,' orbital pairs are set to',f12.7)") iset,tmpval
            end if
        end do
    end if
end do


!Below we will calculate grid data
!Set up grid first
write(*,*)
call setgrid(0,igridsel)
if (allocated(holegrid)) deallocate(holegrid,elegrid,holeeleovlp,transdens,holecross,elecross)
allocate(holegrid(nx,ny,nz),elegrid(nx,ny,nz),holeeleovlp(nx,ny,nz),transdens(nx,ny,nz),holecross(nx,ny,nz),elecross(nx,ny,nz))
holegrid=0D0
elegrid=0D0
holecross=0D0
elecross=0D0
transdens=0D0
if (idomag==1) then !Will also calculate transition magnetic dipole moment density
    if (allocated(magtrdens)) deallocate(magtrdens)
    allocate(magtrdens(nx,ny,nz,3))
    magtrdens=0D0
end if
write(*,"(a)") " Note: During the calculation of cross term of electron, the orbital pairs whose absolute value of coefficient <0.01 will be ignored to save time"
allocate(skippair(nexcitorb))
skippair=.false.
do iexcitorb=1,nexcitorb
    if (abs(exccoeff(iexcitorb))<0.01D0) skippair(iexcitorb)=.true.
end do
write(*,*) "Calculating grid data..."
call walltime(walltime1)
CALL CPU_TIME(time_begin)
ifinish=0
nthreads=getNThreads()
!$OMP PARALLEL DO SHARED(ifinish,holegrid,elegrid,transdens,holecross,elecross,magtrdens) &
!$OMP PRIVATE(i,j,k,tmpx,tmpy,tmpz,orbval,wfnderv,imo,jmo,excwei,iexcitorb,jexcitorb,ileft,jleft,iright,jright,tmpleft,tmpright,idir,jdir,tmpval) &
!$OMP schedule(dynamic) NUM_THREADS(nthreads)
do k=1,nz
    tmpz=orgz+(k-1)*dz
    do j=1,ny
        tmpy=orgy+(j-1)*dy
        do i=1,nx
            tmpx=orgx+(i-1)*dx
            if (idomag==1) then !Study transition magnetic dipole moment density requests orbital 1st-derivative
                call orbderv(2,1,nmo,tmpx,tmpy,tmpz,orbval,wfnderv)
            else
                call orbderv(1,1,nmo,tmpx,tmpy,tmpz,orbval)
            end if
            !Calculate local term of hole and electron
            do iexcitorb=1,nexcitorb
                imo=orbleft(iexcitorb)
                jmo=orbright(iexcitorb)
                excwei=exccoeff(iexcitorb)**2
                if (excdir(iexcitorb)==1) then ! ->
                    holegrid(i,j,k)=holegrid(i,j,k)+excwei*orbval(imo)**2
                    elegrid(i,j,k)=elegrid(i,j,k)+excwei*orbval(jmo)**2
                else ! <-
                    holegrid(i,j,k)=holegrid(i,j,k)-excwei*orbval(imo)**2
                    elegrid(i,j,k)=elegrid(i,j,k)-excwei*orbval(jmo)**2
!                     holegrid(i,j,k)=holegrid(i,j,k)+excwei*orbval(jmo)**2
!                     elegrid(i,j,k)=elegrid(i,j,k)+excwei*orbval(imo)**2
                end if
                transdens(i,j,k)=transdens(i,j,k)+exccoeff(iexcitorb)*orbval(imo)*orbval(jmo)
                if (idomag==1) then
                    if (excdir(iexcitorb)==1) then
                        magtrdens(i,j,k,1)=magtrdens(i,j,k,1)+exccoeff(iexcitorb)*orbval(imo)*(tmpy*wfnderv(3,jmo)-tmpz*wfnderv(2,jmo))
                        magtrdens(i,j,k,2)=magtrdens(i,j,k,2)+exccoeff(iexcitorb)*orbval(imo)*(tmpz*wfnderv(1,jmo)-tmpx*wfnderv(3,jmo))
                        magtrdens(i,j,k,3)=magtrdens(i,j,k,3)+exccoeff(iexcitorb)*orbval(imo)*(tmpx*wfnderv(2,jmo)-tmpy*wfnderv(1,jmo))
                    else !The de-excitation has important influence on the transition magnetic dipole moment density, so must be considered explicitly
                        magtrdens(i,j,k,1)=magtrdens(i,j,k,1)-exccoeff(iexcitorb)*orbval(imo)*(tmpy*wfnderv(3,jmo)-tmpz*wfnderv(2,jmo))
                        magtrdens(i,j,k,2)=magtrdens(i,j,k,2)-exccoeff(iexcitorb)*orbval(imo)*(tmpz*wfnderv(1,jmo)-tmpx*wfnderv(3,jmo))
                        magtrdens(i,j,k,3)=magtrdens(i,j,k,3)-exccoeff(iexcitorb)*orbval(imo)*(tmpx*wfnderv(2,jmo)-tmpy*wfnderv(1,jmo))
                    end if
                end if
            end do
            !Calculate cross term of hole and electron
            do iexcitorb=1,nexcitorb
                !Below cases are skipped:
                !i->l,i->l and i->l,i<-l and i<-l,i<-l, since ileft==jleft.and.iright==jright
                !i->l,j->m, since ileft/=jleft.and.iright/=jright
                !i->l,i<-m and i<-l,j->l, since excdir(iexcitorb)/=excdir(jexcitorb)
                !**If i->l,i<-l should be taken into account is unsolved
                ! Currently only take below cases into account:
                ! Cross term of hole (do <i|j>):     i->l,j->l substract i<-l,j<-l
                ! Cross term of electron (do <l|m>): i->l,i->m substract i<-l,i<-m
                if (skippair(iexcitorb) .eqv. .true.) cycle
                ileft=orbleft(iexcitorb)
                iright=orbright(iexcitorb)
                tmpleft=exccoeff(iexcitorb)*orbval(ileft) !Use temporary variable to save the time for locating element
                tmpright=exccoeff(iexcitorb)*orbval(iright)
                idir=excdir(iexcitorb)
                do jexcitorb=1,nexcitorb
                    if (skippair(jexcitorb) .eqv. .true.) cycle
                    jleft=orbleft(jexcitorb)
                    jright=orbright(jexcitorb)
                    jdir=excdir(jexcitorb)
                    if (idir/=jdir) cycle
                    if (ileft==jleft) then !do <l|m>
                        if (iright==jright) cycle
                        tmpval=tmpright*exccoeff(jexcitorb)*orbval(jright) !Originally virtual orbital
                        if (idir==1) then !->
                            elecross(i,j,k)=elecross(i,j,k)+tmpval
                        else !<-
                            elecross(i,j,k)=elecross(i,j,k)-tmpval
                        end if
                    else if (iright==jright) then !do <i|j>
                        tmpval=tmpleft*exccoeff(jexcitorb)*orbval(jleft) !Originally occupied orbital
                        if (idir==1) then !->
                            holecross(i,j,k)=holecross(i,j,k)+tmpval
                        else !<-
                            holecross(i,j,k)=holecross(i,j,k)-tmpval
                        end if
                    end if
                end do
            end do
        end do
    end do
    ifinish=ifinish+1
    write(*,"(' Finished:',i5,'/',i5)") ifinish,nz
end do
!$OMP END PARALLEL DO
deallocate(skippair)
if (wfntype==0.or.wfntype==3) then !For close-shell wavefunction, the weights are normalized to 0.5 (or say the orbitals are doubly occupied), so correct it
    holegrid=holegrid*2
    elegrid=elegrid*2
    transdens=transdens*2
    holecross=holecross*2
    elecross=elecross*2
    if (idomag==1) magtrdens=magtrdens*2
end if
!Combine local term and cross term of hole to holegrid, that of electron to elegrid. Then local term will be holegrid-holecross and elegrid-elecross
holegrid=holegrid+holecross
elegrid=elegrid+elecross
! holeeleovlp=holegrid*elegrid*min(abs(holegrid),abs(elegrid)) / max(abs(holegrid),abs(elegrid))
! holeeleovlp=holegrid*elegrid
holeeleovlp=min(holegrid,elegrid)

CALL CPU_TIME(time_end)
call walltime(walltime2)
write(*,"(' Totally took up CPU time',f12.2,'s, wall clock time',i10,'s',/)") time_end-time_begin,walltime2-walltime1

!Check normalization
dvol=dx*dy*dz
rnormhole=0
rnormele=0
rnormovlp=0
rinttransdens=0
rtransdipx=0
rtransdipy=0
rtransdipz=0
centholex=0
centholey=0
centholez=0
centelex=0
centeley=0
centelez=0
rtransmagx=0
rtransmagy=0
rtransmagz=0
do k=1,nz
    tmpz=orgz+(k-1)*dz
    do j=1,ny
        tmpy=orgy+(j-1)*dy
        do i=1,nx
            tmpx=orgx+(i-1)*dx
            rnormhole=rnormhole+holegrid(i,j,k)
            rnormele=rnormele+elegrid(i,j,k)
            rnormovlp=rnormovlp+holeeleovlp(i,j,k)
            rinttransdens=rinttransdens+transdens(i,j,k)
            rtransdipx=rtransdipx-tmpx*transdens(i,j,k)
            rtransdipy=rtransdipy-tmpy*transdens(i,j,k)
            rtransdipz=rtransdipz-tmpz*transdens(i,j,k)
            centholex=centholex+holegrid(i,j,k)*tmpx
            centholey=centholey+holegrid(i,j,k)*tmpy
            centholez=centholez+holegrid(i,j,k)*tmpz
            centelex=centelex+elegrid(i,j,k)*tmpx
            centeley=centeley+elegrid(i,j,k)*tmpy
            centelez=centelez+elegrid(i,j,k)*tmpz
            if (idomag==1) then
                rtransmagx=rtransmagx+magtrdens(i,j,k,1)
                rtransmagy=rtransmagy+magtrdens(i,j,k,2)
                rtransmagz=rtransmagz+magtrdens(i,j,k,3)
            end if
        end do
    end do
end do
rnormhole=rnormhole*dvol
rnormele=rnormele*dvol
rnormovlp=rnormovlp*dvol
rinttransdens=rinttransdens*dvol
rtransdipx=rtransdipx*dvol
rtransdipy=rtransdipy*dvol
rtransdipz=rtransdipz*dvol
centholex=centholex*dvol/rnormhole
centholey=centholey*dvol/rnormhole
centholez=centholez*dvol/rnormhole
centelex=centelex*dvol/rnormele
centeley=centeley*dvol/rnormele
centelez=centelez*dvol/rnormele
rtransmagx=rtransmagx*dvol
rtransmagy=rtransmagy*dvol
rtransmagz=rtransmagz*dvol
write(*,"(' Integral of hole:',f12.6)") rnormhole
write(*,"(' Integral of electron:',f12.6)") rnormele
write(*,"(' Integral of transition density:',f12.6)") rinttransdens
write(*,"(' Transition dipole moment in X/Y/Z:',3f11.6,' a.u.')") rtransdipx,rtransdipy,rtransdipz
if (idomag==1) write(*,"(' Transition magnetic dipole moment in X/Y/Z:',3f10.6,' a.u.')") rtransmagx,rtransmagy,rtransmagz
write(*,"(' Integral of overlap of hole-electron distribution:',f12.7)") rnormovlp
write(*,"(' Centroid of hole in X/Y/Z:    ',3f12.6,' Angstrom')") centholex*b2a,centholey*b2a,centholez*b2a
write(*,"(' Centroid of electron in X/Y/Z:',3f12.6,' Angstrom')") centelex*b2a,centeley*b2a,centelez*b2a
write(*,"(' Distance between centroid of hole and electron in X/Y/Z:')")
disx=abs(centelex-centholex)
disy=abs(centeley-centholey)
disz=abs(centelez-centholez)
disnorm=dsqrt(disx**2+disy**2+disz**2)
write(*,"(3f12.6,' Angstrom   Norm:',f12.6' Angstrom')") disx*b2a,disy*b2a,disz*b2a,disnorm*b2a
!Of course, if relaxed density it used, we cannot evaluate the variation of dipole moment here exactly
write(*,"(' Variation of dipole moment with respect to ground state in X/Y/Z:')")
avgtransval=(rnormele+rnormhole)/2D0 !Ideal value is 1.0
dipvarx=-(centelex-centholex)*avgtransval
dipvary=-(centeley-centholey)*avgtransval
dipvarz=-(centelez-centholez)*avgtransval
write(*,"(3f12.6,' a.u.   Norm:',f12.6' a.u.')") dipvarx,dipvary,dipvarz,dsqrt(dipvarx**2+dipvary**2+dipvarz**2)

sigxele=0D0
sigyele=0D0
sigzele=0D0
sigxhole=0D0
sigyhole=0D0
sigzhole=0D0
do k=1,nz
    tmpz=orgz+(k-1)*dz
    do j=1,ny
        tmpy=orgy+(j-1)*dy
        do i=1,nx
            tmpx=orgx+(i-1)*dx
            sigxele=sigxele+elegrid(i,j,k)*(tmpx-centelex)**2
            sigyele=sigyele+elegrid(i,j,k)*(tmpy-centeley)**2
            sigzele=sigzele+elegrid(i,j,k)*(tmpz-centelez)**2
            sigxhole=sigxhole+holegrid(i,j,k)*(tmpx-centholex)**2
            sigyhole=sigyhole+holegrid(i,j,k)*(tmpy-centholey)**2
            sigzhole=sigzhole+holegrid(i,j,k)*(tmpz-centholez)**2
        end do
    end do
end do
sigxele=dsqrt(sigxele*dvol/rnormele)
sigyele=dsqrt(sigyele*dvol/rnormele)
sigzele=dsqrt(sigzele*dvol/rnormele)
signormele=dsqrt(sigxele**2+sigyele**2+sigzele**2)
sigxhole=dsqrt(sigxhole*dvol/rnormhole)
sigyhole=dsqrt(sigyhole*dvol/rnormhole)
sigzhole=dsqrt(sigzhole*dvol/rnormhole)
signormhole=dsqrt(sigxhole**2+sigyhole**2+sigzhole**2)
write(*,"(' RMSD of electron in x,y,z:',3f8.3,'   Total:',f8.3,' Angstrom')") sigxele*b2a,sigyele*b2a,sigzele*b2a,signormele*b2a
write(*,"(' RMSD of hole in x,y,z:    ',3f8.3,'   Total:',f8.3,' Angstrom')") sigxhole*b2a,sigyhole*b2a,sigzhole*b2a,signormhole*b2a
delta_sigCT=dsqrt((sigxele-sigxhole)**2+(sigyele-sigyhole)**2+(sigzele-sigzhole)**2)
write(*,"(' Difference between RMSD of hole and electron:',f8.3,' Angstrom')") delta_sigCT*b2a
Hx=(sigxele+sigxhole)/2D0
Hy=(sigyele+sigyhole)/2D0
Hz=(sigzele+sigzhole)/2D0
Hnorm=dsqrt(Hx**2+Hy**2+Hz**2)
write(*,"(' H index in x,y,z:',3f8.3,'   Norm:',f8.3,' Angstrom')") Hx*b2a,Hy*b2a,Hz*b2a,Hnorm*b2a
tx=disx-Hx
ty=disy-Hy
tz=disz-Hz
tnorm=dsqrt(tx**2+ty**2+tz**2)
write(*,"(' t index in x,y,z:',3f8.3,'   Norm:',f8.3,' Angstrom')") tx*b2a,ty*b2a,tz*b2a,tnorm*b2a

!---- Calculate ghost state diagnostic index proposed by Adamo (DOI: 10.1002/jcc.24862)
ghostp2=1/disnorm !in a.u.
!Definition 1 (the original paper definition). The original paper doesn't clearly mention how to deal with de-excitation, so I treat it as usual
sumC=sum(exccoeff(1:nexcitorb))
ghostp1=0
do itmp=1,nexcitorb
    ghostp1=ghostp1+exccoeff(itmp)/sumC*(-MOene(orbleft(itmp))-MOene(orbright(itmp)))
end do
ghostidx=ghostp1-ghostp2
write(*,"(' Ghost-hunter index (defin. 1):',f8.3,' eV, 1st/2nd terms:',2f9.3,' eV')") ghostidx*au2eV,ghostp1*au2eV,ghostp2*au2eV
!Definition 2 (defined by Tian Lu, more reasonable)
sumCsqr=0
do itmp=1,nexcitorb
    if (excdir(itmp)==2) cycle !Skip de-excitation
    sumCsqr=sumCsqr+exccoeff(itmp)**2
end do
ghostp1=0
do itmp=1,nexcitorb
    if (excdir(itmp)==2) cycle !Skip de-excitation
    ghostp1=ghostp1+exccoeff(itmp)**2/sumCsqr*(-MOene(orbleft(itmp))-MOene(orbright(itmp)))
end do
ghostidx=ghostp1-ghostp2
write(*,"(' Ghost-hunter index (defin. 2):',f8.3,' eV, 1st/2nd terms:',2f9.3,' eV')") ghostidx*au2eV,ghostp1*au2eV,ghostp2*au2eV
write(*,"(' Excitation energy of this state:',f10.3,' eV')") excenergy
if (excenergy<ghostidx*au2eV) write(*,*) "Warning: Probably this is a ghost state"

if (allocated(Cele)) deallocate(Cele,Chole)
allocate(Cele(nx,ny,nz),Chole(nx,ny,nz))
do i=1,nx
    do j=1,ny
        do k=1,nz
            rnowx=orgx+(i-1)*dx
            rnowy=orgy+(j-1)*dy
            rnowz=orgz+(k-1)*dz
            Cele(i,j,k)=exp( -(rnowx-centelex)**2/(2*sigxele**2) -(rnowy-centeley)**2/(2*sigyele**2) -(rnowz-centelez)**2/(2*sigzele**2))
            Chole(i,j,k)=exp( -(rnowx-centholex)**2/(2*sigxhole**2) -(rnowy-centholey)**2/(2*sigyhole**2) -(rnowz-centholez)**2/(2*sigzhole**2))
        end do
    end do
end do
Cele=Cele*rnormele/(sum(Cele)*dvol)
Chole=Chole*rnormhole/(sum(Chole)*dvol)

do while(.true.)
    write(*,*)
    write(*,*) "0 Return"
    write(*,*) "1 Show isosurface of hole distribution"
    write(*,*) "2 Show isosurface of electron distribution"
    write(*,*) "3 Show isosurface of hole and electron distribution simultaneously"
    write(*,*) "4 Show isosurface of overlap of hole-electron"
    write(*,*) "5 Show isosurface of transition density"
    write(*,*) "6 Show isosurface of transition dipole moment density"
    write(*,*) "7 Show isosurface of charge density difference"
    write(*,*) "8 Show isosurface of Cele and Chole functions simultaneously"
    if (idomag==1) write(*,*) "9 Show isosurface of transition magnetic dipole moment density"
    write(*,*) "10 Output cube file of hole distribution to current folder"
    write(*,*) "11 Output cube file of electron distribution to current folder"
    write(*,*) "12 Output cube file of overlap of hole-electron to current folder"
    write(*,*) "13 Output cube file of transition density to current folder"
    write(*,*) "14 Output cube file of transition dipole moment density to current folder"
    write(*,*) "15 Output cube file of charge density difference to current folder"
    write(*,*) "16 Output cube file of Cele and Chole functions to current folder"
    write(*,*) "17 Output cube file of transition magnetic dipole moment density"
    write(*,*) "18 Calculate hole-electron Coulomb attractive energy"
    read(*,*) isel
    if (isel==0) then
        goto 1
    else if (isel==1.or.isel==2.or.isel==4.or.isel==5.or.isel==6.or.isel==7.or.isel==9) then
         if (allocated(cubmat)) deallocate(cubmat)
         allocate(cubmat(nx,ny,nz))
         if (isel==1) then
             write(*,*) "Select the type of hole distribution to be shown. 1 is commonly used"
             write(*,*) "1 Total (local term + cross term)"
             write(*,*) "2 Local term only"
             write(*,*) "3 Cross term only"
             read(*,*) isel2
             if (isel2==1) then
                 cubmat=holegrid
             else if (isel2==2) then
                 cubmat=holegrid-holecross
             else if (isel2==3) then
                 cubmat=holecross
             end if
             sur_value=0.002D0
         else if (isel==2) then
             write(*,*) "Select the type of electron distribution to be shown. 1 is commonly used"
             write(*,*) "1 Total (local term + cross term)"
             write(*,*) "2 Local term only"
             write(*,*) "3 Cross term only"
             read(*,*) isel2
             if (isel2==1) then
                 cubmat=elegrid
             else if (isel2==2) then
                 cubmat=elegrid-elecross
             else if (isel2==3) then
                 cubmat=elecross
             end if
             sur_value=0.002D0
         else if (isel==4) then
             cubmat=holeeleovlp
             sur_value=0.002D0
         else if (isel==5) then
             cubmat=transdens
            isosur1style=4
             sur_value=0.001D0
        else if (isel==6) then
            write(*,*) "Select the component of transition dipole moment density"
            write(*,*) "1: X component  2: Y component  3: Z component"
             read(*,*) ifac
            isosur1style=4
             sur_value=0.001D0
             do k=1,nz
                tmpz=orgz+(k-1)*dz
                do j=1,ny
                    tmpy=orgy+(j-1)*dy
                    do i=1,nx
                        tmpx=orgx+(i-1)*dx
                        if (ifac==1) cubmat(i,j,k)=-tmpx*transdens(i,j,k)
                        if (ifac==2) cubmat(i,j,k)=-tmpy*transdens(i,j,k)
                        if (ifac==3) cubmat(i,j,k)=-tmpz*transdens(i,j,k)
                    end do
                end do
            end do
         else if (isel==7) then
             cubmat=elegrid-holegrid
             sur_value=0.002D0
         else if (isel==9) then
            write(*,*) "Select the component of transition magnetic dipole moment density"
            write(*,*) "1: X component  2: Y component  3: Z component"
             read(*,*) ifac
            isosur1style=4
             sur_value=0.007D0
            if (ifac==1) cubmat=magtrdens(:,:,:,1)
            if (ifac==2) cubmat=magtrdens(:,:,:,2)
            if (ifac==3) cubmat=magtrdens(:,:,:,3)
         end if
        deallocate(cubmat)
    else if (isel==3.or.isel==8) then
        if (isel==3) write(*,"(a)") " Note: Blue and green isosurfaces represent hole and electron distributions, respectively"
        if (isel==8) write(*,"(a)") " Note: Blue and green isosurfaces represent Chole and Cele functions, respectively"
        sur_value=0.002D0
        isosursec=1
        clrRcub2sameold=clrRcub2same !Backup previous color setting
        clrGcub2sameold=clrGcub2same
        clrBcub2sameold=clrBcub2same
        clrRcub2samemeshptold=clrRcub2samemeshpt
        clrGcub2samemeshptold=clrGcub2samemeshpt
        clrBcub2samemeshptold=clrBcub2samemeshpt
        clrRcub2same=0.3D0
        clrGcub2same=0.45D0
        clrBcub2same=0.9D0
        clrRcub2samemeshpt=0.3D0
        clrGcub2samemeshpt=0.45D0
        clrBcub2samemeshpt=0.9D0
        if (allocated(cubmat)) deallocate(cubmat)
        if (allocated(cubmattmp)) deallocate(cubmattmp)
        allocate(cubmat(nx,ny,nz),cubmattmp(nx,ny,nz))
        if (isel==3) then
            cubmat=elegrid
            cubmattmp=holegrid
        else if (isel==8) then
            cubmat=Cele
            cubmattmp=Chole
        end if
        isosur1style=4
        isosur2style=4
        deallocate(cubmat,cubmattmp)
        clrRcub2same=clrRcub2sameold !Recover previous color setting
        clrGcub2same=clrGcub2sameold
        clrBcub2same=clrBcub2sameold
        clrRcub2samemeshpt=clrRcub2samemeshptold
        clrGcub2samemeshpt=clrGcub2samemeshptold
        clrBcub2samemeshpt=clrBcub2samemeshptold
    else if (isel==10) then
         write(*,*) "Select the type of hole distribution to be exported. 1 is commonly used"
         write(*,*) "1 Total (local term + cross term)"
         write(*,*) "2 Local term only"
         write(*,*) "3 Cross term only"
         read(*,*) isel2
        write(*,*) "Outputting hole distribution to hole.cub in current folder"
        open(10,file="hole.cub",status="replace")
        if (isel2==1) call outcube(holegrid,nx,ny,nz,orgx,orgy,orgz,dx,dy,dz,10)
        if (isel2==2) call outcube(holegrid-holecross,nx,ny,nz,orgx,orgy,orgz,dx,dy,dz,10)
        if (isel2==3) call outcube(holecross,nx,ny,nz,orgx,orgy,orgz,dx,dy,dz,10)
        close(10)
        write(*,*) "Done!"
    else if (isel==11) then
         write(*,*) "Select the type of electron distribution to be exported. 1 is commonly used"
         write(*,*) "1 Total (local term + cross term)"
         write(*,*) "2 Local term only"
         write(*,*) "3 Cross term only"
         read(*,*) isel2
        write(*,*) "Outputting electron distribution to electron.cub in current folder"
        open(10,file="electron.cub",status="replace")
        if (isel2==1) call outcube(elegrid,nx,ny,nz,orgx,orgy,orgz,dx,dy,dz,10)
        if (isel2==2) call outcube(elegrid-elecross,nx,ny,nz,orgx,orgy,orgz,dx,dy,dz,10)
        if (isel2==3) call outcube(elecross,nx,ny,nz,orgx,orgy,orgz,dx,dy,dz,10)
        close(10)
        write(*,*) "Done!"
    else if (isel==12) then
        write(*,"(a)") " Outputting overlap of hole-electron to holeeleovlp.cub in current folder"
        open(10,file="holeeleovlp.cub",status="replace")
        call outcube(holeeleovlp,nx,ny,nz,orgx,orgy,orgz,dx,dy,dz,10)
        close(10)
        write(*,*) "Done!"
    else if (isel==13) then
        write(*,"(a)") " Outputting transition density to transdens.cub in current folder"
        open(10,file="transdens.cub",status="replace")
        call outcube(transdens,nx,ny,nz,orgx,orgy,orgz,dx,dy,dz,10)
        close(10)
        write(*,*) "Done!"
     else if (isel==14) then
        write(*,*) "Select the component of transition dipole moment density"
        write(*,*) "1: X component  2: Y component  3: Z component"
         read(*,*) ifac
         if (allocated(cubmat)) deallocate(cubmat)
         allocate(cubmat(nx,ny,nz))
         do k=1,nz
            tmpz=orgz+(k-1)*dz
            do j=1,ny
                tmpy=orgy+(j-1)*dy
                do i=1,nx
                    tmpx=orgx+(i-1)*dx
                    if (ifac==1) cubmat(i,j,k)=-tmpx*transdens(i,j,k)
                    if (ifac==2) cubmat(i,j,k)=-tmpy*transdens(i,j,k)
                    if (ifac==3) cubmat(i,j,k)=-tmpz*transdens(i,j,k)
                end do
            end do
        end do
        write(*,*) "Outputting to transdipdens.cub in current folder"
        open(10,file="transdipdens.cub",status="replace")
        call outcube(cubmat,nx,ny,nz,orgx,orgy,orgz,dx,dy,dz,10)
        close(10)
        write(*,*) "Done!"
        deallocate(cubmat)
    else if (isel==15) then
        write(*,*) "Outputting charge density difference to CDD.cub in current folder"
        open(10,file="CDD.cub",status="replace")
        call outcube(elegrid-holegrid,nx,ny,nz,orgx,orgy,orgz,dx,dy,dz,10)
        close(10)
        write(*,*) "Done!"
    else if (isel==16) then
        open(10,file="Cele.cub",status="replace")
        call outcube(Cele,nx,ny,nz,orgx,orgy,orgz,dx,dy,dz,10)
        close(10)
        write(*,"('Cele function has been outputted to ""Cele.cub"" in current folder')")
        open(10,file="Chole.cub",status="replace")
        call outcube(Chole,nx,ny,nz,orgx,orgy,orgz,dx,dy,dz,10)
        close(10)
        write(*,"('Chole function has been outputted to ""Chole.cub"" in current folder')")
     else if (isel==17) then
        write(*,*) "Select the component of transition magnetic dipole moment density"
        write(*,*) "1: X component  2: Y component  3: Z component"
         read(*,*) ifac
        write(*,*) "Outputting to magtrdipdens.cub in current folder"
        open(10,file="magtrdipdens.cub",status="replace")
        if (ifac==1) call outcube(magtrdens(:,:,:,1),nx,ny,nz,orgx,orgy,orgz,dx,dy,dz,10)
        if (ifac==2) call outcube(magtrdens(:,:,:,2),nx,ny,nz,orgx,orgy,orgz,dx,dy,dz,10)
        if (ifac==3) call outcube(magtrdens(:,:,:,3),nx,ny,nz,orgx,orgy,orgz,dx,dy,dz,10)
        close(10)
        write(*,*) "Done!"
     else if (isel==18) then
        call walltime(iwalltime1)
        CALL CPU_TIME(time_begin)
        write(*,*) "Calculating, please wait..."
        allocate(cubx(nx),cuby(ny),cubz(nz))
        do i=1,nx
            cubx(i)=orgx+(i-1)*dx
        end do
        do i=1,ny
            cuby(i)=orgy+(i-1)*dy
        end do
        do i=1,nz
            cubz(i)=orgz+(i-1)*dz
        end do
         coulene=0
        do i=1,nx
            do j=1,ny
                do k=1,nz
                    if (Cele(i,j,k)<1D-6) cycle !Typically leads to error at 0.001 magnitude
nthreads=getNThreads()
!$OMP parallel shared(coulene) private(ii,jj,kk,distx2,disty2,distz2,dist,coulenetmp) num_threads(nthreads)
                    coulenetmp=0
!$OMP do schedule(DYNAMIC)
                    do ii=1,nx
                        distx2=(cubx(i)-cubx(ii))**2
                        do jj=1,ny
                            disty2=(cuby(j)-cuby(jj))**2
                            do kk=1,nz
                                if (i==ii.and.j==jj.and.k==kk) cycle
                                distz2=(cubz(k)-cubz(kk))**2
                                dist=dsqrt(distx2+disty2+distz2)
                                coulenetmp=coulenetmp+Cele(i,j,k)*Chole(ii,jj,kk)/dist
                            end do
                        end do
                    end do
!$OMP END DO
!$OMP CRITICAL
                    coulene=coulene+coulenetmp
!$OMP END CRITICAL
!$OMP END PARALLEL
                end do
            end do
            write(*,"(' Finished:',i4,' /',i4)") i,nx
        end do
        dvol=dx*dy*dz
        coulene=-coulene*dvol*dvol
        CALL CPU_TIME(time_end)
        call walltime(iwalltime2)
        write(*,"(' Calculation took up CPU time',f12.2,'s, wall clock time',i10,'s')") time_end-time_begin,iwalltime2-iwalltime1
        write(*,*)
        write(*,"(' Coulomb attractive energy:',f12.6,' a.u.  (',f12.6,' eV )')") coulene,coulene*au2eV
        deallocate(cubx,cuby,cubz)
    end if
end do
end subroutine




!!--- Calculate delta_r index, see J. Chem. Theory Comput., 9, 3118 (2013)
!excitation and de-excitation cofficients are summed together, according to Eq.9 of the original paper
subroutine calcdelta_r(nexcitorb,orbleft,orbright,excdir,exccoeff)
use defvar
implicit real*8 (a-h,o-z)
integer nexcitorb,orbleft(nexcitorb),orbright(nexcitorb),excdir(nexcitorb)
real*8 exccoeff(nexcitorb)
real*8,allocatable :: exccoefftot(:) !Store the coefficient combined from excitation and de-excitation of the same pair
real*8,allocatable :: orbcenx(:),orbceny(:),orbcenz(:) !Store centroid of MOs
real*8,allocatable :: GTFdipint(:,:)
character strtmp1,strtmp2,selectyn
write(*,*) "Calculating, please wait..."
allocate(exccoefftot(nexcitorb))
exccoefftot=exccoeff
coeffsumsqr=0D0
do itmp=1,nexcitorb
    if (excdir(itmp)==1) then !->, find corresponding <- pair and sum its coefficient to here
        do jtmp=1,nexcitorb
            if (excdir(jtmp)==2.and.orbleft(itmp)==orbleft(jtmp).and.orbright(itmp)==orbright(jtmp)) exccoefftot(itmp)=exccoefftot(itmp)+exccoeff(jtmp)
        end do
    end if
    coeffsumsqr=coeffsumsqr+exccoefftot(itmp)**2
end do
!Calculate dipole moment integral matrix
nsize=nprims*(nprims+1)/2
allocate(GTFdipint(3,nsize))
call genGTFDmat(GTFdipint,nsize)
!Calculate centroid of MOs
allocate(orbcenx(nmo),orbceny(nmo),orbcenz(nmo))
orbcenx=0
orbceny=0
orbcenz=0
nthreads=getNThreads()
!$OMP PARALLEL DO SHARED(orbcenx,orbceny,orbcenz) PRIVATE(iGTF,jGTF,ides,imo) schedule(dynamic) NUM_THREADS(nthreads)
do imo=1,nmo
    do iGTF=1,nprims
        do jGTF=1,nprims
            if (iGTF>=jGTF) then
                ides=iGTF*(iGTF-1)/2+jGTF
            else
                ides=jGTF*(jGTF-1)/2+iGTF
            end if
            orbcenx(imo)=orbcenx(imo)+co(imo,iGTF)*co(imo,jGTF)*GTFdipint(1,ides)
            orbceny(imo)=orbceny(imo)+co(imo,iGTF)*co(imo,jGTF)*GTFdipint(2,ides)
            orbcenz(imo)=orbcenz(imo)+co(imo,iGTF)*co(imo,jGTF)*GTFdipint(3,ides)
        end do
    end do
end do
!$OMP END PARALLEL DO
deallocate(GTFdipint)
delta_r=0
do itmp=1,nexcitorb
    if (excdir(itmp)==2) cycle
    imo=orbleft(itmp)
    jmo=orbright(itmp)
    delta_r=delta_r+exccoefftot(itmp)**2/coeffsumsqr *dsqrt((orbcenx(imo)-orbcenx(jmo))**2+(orbceny(imo)-orbceny(jmo))**2+(orbcenz(imo)-orbcenz(jmo))**2)
end do
write(*,"(' Delta_r=',f12.6,' Bohr,',f12.6,' Angstrom')") delta_r,delta_r*b2a
write(*,*)
write(*,*) "If print orbital pair contribution to delta_r? (y/n)"
read(*,*) selectyn
if (selectyn=='y'.or.selectyn=='Y') then
    write(*,*) "Input threshold for printing e.g. 0.05"
    write(*,"(a)") " Note: If input -1, then all contributions will be exported to delta_r.txt in curren folder"
    read(*,*) printthres
    iout=6
    if (printthres==-1) then
        iout=10
        open(10,file="delta_r.txt",status="replace")
    end if
    write(iout,"(a)") " Note: The transition coefficients shown below have combined both excitation and de-excitation parts"
    write(iout,"(' Sum of square of transition coefficient:',f12.6)") coeffsumsqr
    write(iout,*) "   #Pair     Orbitals      Coefficient     Contribution (Bohr and Angstrom)"
    do itmp=1,nexcitorb
        if (excdir(itmp)==2) cycle
        imo=orbleft(itmp)
        jmo=orbright(itmp)
        contrival=exccoefftot(itmp)**2/coeffsumsqr *dsqrt((orbcenx(imo)-orbcenx(jmo))**2+(orbceny(imo)-orbceny(jmo))**2+(orbcenz(imo)-orbcenz(jmo))**2)
        if (contrival<printthres) cycle
        if (wfntype==0.or.wfntype==3) then
            write(iout,"(i8,2i7,f16.7,3x,2f16.7)") itmp,imo,jmo,exccoefftot(itmp),contrival,contrival*b2a
        else
            strtmp1="A"
            strtmp2="A"
            if (imo>nbasis) then
                imo=imo-nbasis
                strtmp1="B"
            end if
            if (jmo>nbasis) then
                jmo=jmo-nbasis
                strtmp2="B"
            end if
            write(iout,"(i8,i6,a,i6,a,f16.7,3x,2f16.7)") itmp,imo,strtmp1,jmo,strtmp2,exccoefftot(itmp),contrival,contrival*b2a
        end if
    end do
    if (printthres==-1) then
        close(10)
        write(*,*) "Done, outputting finished"
    end if
end if
deallocate(exccoefftot,orbcenx,orbceny,orbcenz)
write(*,*)
end subroutine




!!---- Generate NTOs, original paper of NTO: J. Chem. Phys., 118, 4775 (2003)
!The NTO eigenvalues for unrestricted wavefunction is 1/2 of the ones outputted by Gaussian, &
!which is incorrect (i.e. the sum is 2.0 rather than 1.0), the result must be Gaussian used incorrect factor
subroutine NTO(nexcitorb,orbleft,orbright,excdir,exccoeff)
use defvar
use util
implicit real*8 (a-h,o-z)
real*8,allocatable :: T_MO(:,:),TT(:,:),NTOvec(:,:),NTOval(:)
integer nexcitorb,orbleft(nexcitorb),orbright(nexcitorb),excdir(nexcitorb)
real*8 exccoeff(nexcitorb),tmparr(nbasis)
character c200tmp*200
write(*,*)
NTOvalcoeff=2
if (allocated(CObasb)) then
    NTOvalcoeff=1
    write(*,*) "Result of Alpha part:"
end if
!*** Alpha part or restricted case
nocc=nint(naelec)
nvir=nbasis-nocc
allocate(T_MO(nocc,nvir))
T_MO=0
do iexcitorb=1,nexcitorb
    iocc=orbleft(iexcitorb)
    ivir=orbright(iexcitorb)-nocc
    if (iocc>nbasis) cycle !Here we only process alpha part, index of Beta orbitals are higher than nbasis
    !Ignoring de-excitation, this treatment makes the result identical to that ouputted by Gaussian
    if (excdir(iexcitorb)==1) T_MO(iocc,ivir)=T_MO(iocc,ivir)+exccoeff(iexcitorb)
end do
!Occupied part
allocate(TT(nocc,nocc),NTOvec(nocc,nocc),NTOval(nocc))
TT=matmul(T_MO,transpose(T_MO))
call diagsymat(TT,NTOvec,NTOval,ierror)
NTOval=NTOval*NTOvalcoeff
MOene(1:nocc)=NTOval !By default, the diagsymat gives result from low to high
CObasa(:,1:nocc)=matmul(CObasa(:,1:nocc),NTOvec)
if (nocc>10) then
    write(*,*) "The highest 10 eigenvalues of occupied NTOs:"
    write(*,"(5f12.6)") MOene(nocc-9:nocc)
else
    write(*,*) "Eigenvalues of occupied NTOs:"
    write(*,"(5f12.6)") MOene(1:nocc)
end if
deallocate(TT,NTOvec,NTOval)
!Virtual part
allocate(TT(nvir,nvir),NTOvec(nvir,nvir),NTOval(nvir))
TT=matmul(transpose(T_MO),T_MO)
call diagsymat(TT,NTOvec,NTOval,ierror)
NTOval=NTOval*NTOvalcoeff
MOene(nocc+1:nbasis)=NTOval
CObasa(:,nocc+1:nbasis)=matmul(CObasa(:,nocc+1:nbasis),NTOvec)
do itmp=1,int(nvir/2D0) !Exchange array, so that the sequence will be high->low rather than the default low->high
    i=nocc+itmp
    j=nbasis+1-itmp
    tmpval=MOene(i)
    MOene(i)=MOene(j)
    MOene(j)=tmpval
    tmparr=CObasa(:,i)
    CObasa(:,i)=CObasa(:,j)
    CObasa(:,j)=tmparr
end do
write(*,*)
if (nvir>10) then
    write(*,*) "The highest 10 eigenvalues of virtual NTOs:"
    write(*,"(5f12.6)") MOene(nocc+1:nocc+10)
else
    write(*,*) "Eigenvalues of virtual NTOs:"
    write(*,"(5f12.6)") MOene(nocc+1:nbasis)
end if
deallocate(TT,NTOvec,NTOval,T_MO)

!*** Beta part
if (allocated(CObasb)) then
    write(*,*)
    write(*,*) "Result of Beta part:"
    nocc=nint(nbelec)
    nvir=nbasis-nocc
    allocate(T_MO(nocc,nvir))
    T_MO=0
    do iexcitorb=1,nexcitorb
        iocc=orbleft(iexcitorb)
        ivir=orbright(iexcitorb)-nocc
        if (iocc>nbasis) then !Only process transition of Beta orbitals
            iocc=iocc-nbasis
            ivir=ivir-nbasis
            if (excdir(iexcitorb)==1) T_MO(iocc,ivir)=T_MO(iocc,ivir)+exccoeff(iexcitorb)
        end if
    end do
    !Occupied part
    allocate(TT(nocc,nocc),NTOvec(nocc,nocc),NTOval(nocc))
    TT=matmul(T_MO,transpose(T_MO))
    call diagsymat(TT,NTOvec,NTOval,ierror)
    NTOval=NTOval*NTOvalcoeff
    MOene(nbasis+1:nbasis+nocc)=NTOval !By default, the diagsymat gives result from low to high
    CObasb(:,1:nocc)=matmul(CObasb(:,1:nocc),NTOvec)
    if (nocc>10) then
        write(*,*) "The highest 10 eigenvalues of occupied NTOs:"
        write(*,"(5f12.6)") MOene(nbasis+nocc-9:nbasis+nocc)
    else
        write(*,*) "Eigenvalues of occupied NTOs:"
        write(*,"(5f12.6)") MOene(nbasis+1:nbasis+nocc)
    end if
    deallocate(TT,NTOvec,NTOval)
    !Virtual part
    allocate(TT(nvir,nvir),NTOvec(nvir,nvir),NTOval(nvir))
    TT=matmul(transpose(T_MO),T_MO)
    call diagsymat(TT,NTOvec,NTOval,ierror)
    NTOval=NTOval*NTOvalcoeff
    MOene(nbasis+nocc+1:nbasis+nbasis)=NTOval
    CObasb(:,nocc+1:nbasis)=matmul(CObasb(:,nocc+1:nbasis),NTOvec)
    do itmp=1,int(nvir/2D0) !Exchange array, so that the sequence will be high->low rather than the default low->high
        i=nocc+itmp
        j=nbasis+1-itmp
        tmpval=MOene(nbasis+i)
        MOene(nbasis+i)=MOene(nbasis+j)
        MOene(nbasis+j)=tmpval
        tmparr=CObasb(:,i)
        CObasb(:,i)=CObasb(:,j)
        CObasb(:,j)=tmparr
    end do
    write(*,*)
    if (nvir>10) then
        write(*,*) "The highest 10 eigenvalues of virtual NTOs:"
        write(*,"(5f12.6)") MOene(nbasis+nocc+1:nbasis+nocc+10)
    else
        write(*,*) "Eigenvalues of virtual NTOs:"
        write(*,"(5f12.6)") MOene(nbasis+nocc+1:nbasis+nbasis)
    end if
    deallocate(TT,NTOvec,NTOval,T_MO)
end if
write(*,*)
write(*,*) "0 Return"
write(*,*) "1 Output NTO orbitals to .molden file"
write(*,*) "2 Output NTO orbitals to .fch file"
read(*,*) iselNTO
if (iselNTO==1) then
    write(*,*) "Input the file path to output, e.g. C:\S1.molden"
    read(*,"(a)") c200tmp
    call outmolden(c200tmp,10)
    write(*,*) "Now you can load the newly generated .molden file to visualize NTOs"
else if (iselNTO==2) then
    write(*,*) "Input the file path to output, e.g. C:\S1.fch"
    read(*,"(a)") c200tmp
    call outfch(c200tmp,10)
    write(*,*) "Now you can load the newly generated .fch file to visualize NTOs"
end if
write(*,*)
write(*,*) "Reloading the initial file to recover status..."
call dealloall
call readinfile(firstfilename,1)
write(*,*) "Loading finished!"
write(*,*)
end subroutine




!!---------- Calculate all transition dipole moments between all excited states
!Note that this function cannot borrow the code in hetransdipdens for loading data, because this function &
!need transition information of all states as well as excitation energies
subroutine exctransdip
use defvar
use util
use function
implicit real*8 (a-h,o-z)
integer itype
!The last index in each arrays is the state index
integer,allocatable :: excdir(:,:) !excdir=1 means ->, =2 means <-
integer,allocatable :: orbleft(:,:),orbright(:,:) !denote the actual MO (have already considered alpha/beta problem) at the left/right side in the excitation data
real*8,allocatable :: exccoeff(:,:) !Coefficient of an orbital pair transition
real*8,allocatable :: excenergy(:) !Excitation energies
integer,allocatable :: excmulti(:) !Multiplicity of the states
integer,allocatable :: nexcitorb(:) !The number of MO pairs in the states
character c80tmp*80,transmodestr*200,leftstr*80,rightstr*80
character,save :: excitfilename*200=" "
real*8,allocatable :: GTFdipint(:,:) !Dipole moment integral between GTFs, use compressed index. The first index is x,y,z
real*8,allocatable :: MOdipint(:,:,:) !Dipole moment integral between all MOs. The first index is x,y,z
real*8,allocatable :: tdvecmat(:,:,:)
real*8 tdvec(3),grounddip(3),tmpvec(3)

! excitfilename="c:\gtest\acetic_acid.out"
if (excitfilename==" ") then
    write(*,*) "Input the path of the Gaussian/ORCA output file or plain text file"
    write(*,*) "e.g. c:\lovelive\sunshine\yosoro.out"
    do while(.true.)
        read(*,"(a)") excitfilename
        inquire(file=excitfilename,exist=alive)
        if (alive) exit
        write(*,*) "Cannot find this file, input again"
    end do
else
    write(*,"(' Loading ',a)") trim(excitfilename)
end if
open(10,file=excitfilename,status="old")
call loclabel(10,"Gaussian",igauout,maxline=50)
call loclabel(10,"O   R   C   A",iORCAout,maxline=50)
rewind(10)

!Determine the number of excited states, so that proper size of arrays can be allocated
nstates=0
if (igauout==1) then !Gaussian output file
    write(*,*) "This is a Gaussian output file"
    call loclabel(10,"Excitation energies and oscillator strengths:")
    do while(.true.)
        call loclabel(10,"Excited State",ifound,0)
        if (ifound==1) then
            nstates=nstates+1
            read(10,*)
        else
            exit
        end if
    end do
else if (iORCAout==1) then !ORCA output file
    write(*,*) "This is an ORCA output file"
    call loclabel(10,"Number of roots to be determined",ifound)
    read(10,"(50x,i7)") nstates
else !Plain text file
    do while(.true.)
        call loclabel(10,"Excited State",ifound,0)
        if (ifound==1) then
            read(10,*)
            nstates=nstates+1
        else
            exit
        end if
    end do
end if
write(*,"(' There are',i5,' transitions, loading...')") nstates

allocate(excenergy(nstates),excmulti(nstates),nexcitorb(nstates),tdvecmat(3,nstates,nstates))
nexcitorb=0

!Load excitation energy, multiplicity, the number of MO pairs of each excited state
if (igauout==1) then !Gaussian output file
    call loclabel(10,"Excitation energies and oscillator strengths:")
    do iexc=1,nstates
        call loclabel(10,"Excited State",ifound,0)
        read(10,"(a)") transmodestr
        excmulti(iexc)=0 !Multiplicity of the excited state
        if (index(transmodestr,"Singlet")/=0) excmulti(iexc)=1
        if (index(transmodestr,"Triplet")/=0) excmulti(iexc)=3
        do i=10,70
            if (transmodestr(i:i+1)=="eV") exit
        end do
        read(transmodestr(i-10:i-1),*) excenergy(iexc)
        !Count how many orbital pairs are involved in this transition mode
        do while(.true.)
            read(10,"(a)",iostat=ierror) c80tmp
            if (index(c80tmp,'<-')==0.and.index(c80tmp,'->')==0) exit
            nexcitorb(iexc)=nexcitorb(iexc)+1
        end do
    end do
else if (iORCAout==1) then !ORCA output file
    excmulti=1 !Multiplicity of all excited states are assumed to be singlet
    if (wfntype==0.or.wfntype==3) then
        call loclabel(10,"Generation of triplets")
        read(10,"(a)") c80tmp
        if (index(c80tmp," on ")/=0) then
            write(*,*) "Load which kind of excited states?"
            write(*,*) "1: Singlet   3: Triplet"
            read(*,*) iexcmulti
            excmulti=iexcmulti
        end if
    end if
    call loclabel(10,"the weight of the individual excitations are printed")
    if (iexcmulti==3) then !When triplets=on, ORCA calculate both singlet and triplet excited state, now move to the latter
        read(10,*)
        call loclabel(10,"the weight of the individual excitations are printed",ifound,0)
    end if
    do iexc=1,nstates
        call loclabel(10,"STATE",ifound,0)
        read(10,"(a)") transmodestr
        do i=10,70
            if (transmodestr(i:i+1)=="eV") exit
        end do
        read(transmodestr(i-10:i-1),*) excenergy(iexc)
        !Count how many orbital pairs are involved in this transition mode
        do while(.true.)
            read(10,"(a)",iostat=ierror) c80tmp
            if (index(c80tmp,'<-')==0.and.index(c80tmp,'->')==0) exit
            nexcitorb(iexc)=nexcitorb(iexc)+1
        end do
    end do
else
    rewind(10)
    do iexc=1,nstates
        call loclabel(10,"Excited State",ifound,0)
        read(10,*) c80tmp,c80tmp,inouse,excmulti(iexc),excenergy(iexc)
        !Count how many orbital pairs are involved in this transition mode
        do while(.true.)
            read(10,"(a)",iostat=ierror) c80tmp
            if ((index(c80tmp,'<-')==0.and.index(c80tmp,'->')==0).or.ierror/=0) exit
            nexcitorb(iexc)=nexcitorb(iexc)+1
        end do
    end do
end if
maxpair=maxval(nexcitorb)

allocate(excdir(maxpair,nstates),orbleft(maxpair,nstates),orbright(maxpair,nstates),exccoeff(maxpair,nstates))

!Load MO transition coefficients and direction
if (igauout==1) then
    call loclabel(10,"Excitation energies and oscillator strengths:")
else if (iORCAout==1) then
    call loclabel(10,"the weight of the individual excitations are printed")
    if (iexcmulti==3) then !When triplets=on, ORCA calculate both singlet and triplet excited state, now move to the latter
        read(10,*)
        call loclabel(10,"the weight of the individual excitations are printed",ifound,0)
    end if
else
    rewind(10)
end if
!Notice that for unrestricted case, A and B are separately recorded in input file, &
!however after loading, they are combined as single index, namely if orbital index is larger than nbasis, then it is B, else A
if (iORCAout==1) then !ORCA output file
    !Worthnotingly, in at least ORCA 4.0, de-excitation is not separately outputted as <-, but combined into ->
    !Here we still determine <-, because hopefully Neese may change the convention of ORCA output in the future...
    do iexc=1,nstates
        call loclabel(10,"STATE",ifound,0)
        read(10,*)
        do itmp=1,nexcitorb(iexc)
            read(10,"(a)") c80tmp
            excdir(itmp,iexc)=-1
            if (index(c80tmp,'->')/=0) excdir(itmp,iexc)=1 !means ->
            if (index(c80tmp,'<-')/=0) excdir(itmp,iexc)=2 !means <-
!              write(*,*) trim(c80tmp)
            do isign=1,80 !Find position of <- or ->
                if (c80tmp(isign:isign)=='-'.or.c80tmp(isign:isign)=='<') exit
            end do
            !Process left side of <- or ->
            read(c80tmp(:isign-1),"(a)") leftstr
            read(leftstr(:len_trim(leftstr)-1),*) orbleft(itmp,iexc)
            orbleft(itmp,iexc)=orbleft(itmp,iexc)+1 !ORCA counts orbital from 0 rather than 1!!!
            if (index(leftstr,'b')/=0) orbleft(itmp,iexc)=orbleft(itmp,iexc)+nbasis
            !Process right side of <- or ->
            read(c80tmp(isign+2:),*) rightstr
            read(rightstr(:len_trim(rightstr)-1),*) orbright(itmp,iexc)
            orbright(itmp,iexc)=orbright(itmp,iexc)+1
            if (index(rightstr,'b')/=0) orbright(itmp,iexc)=orbright(itmp,iexc)+nbasis
            iTDA=index(c80tmp,'c=')
            if (iTDA/=0) then !CIS, TDA task, configuration coefficients are presented
                read(c80tmp(iTDA+2:iTDA+13),*) exccoeff(itmp,iexc)
            else !TD task, configuration coefficients are not presented. Contribution of i->a and i<-a are summed up and outputted as i->a
                if (iexc==1.and.itmp==1) then
                    write(*,"(a)") " Warning: For TD task, ORCA does not print configuration coefficients but only print corresponding contributions of each orbital pair, &
                    in this case Multiwfn determines configuration coefficients simply as square root of contribution values. However, this treatment is &
                    evidently inappropriate and the result is nonsense when de-excitation is significant (In this situation you have to use TDA-DFT instead)"
                    write(*,*) "If you really want to proceed, press ENTER to continue"
                    read(*,*)
                end if
                read(c80tmp(23:32),*) tmpval
                if (tmpval<0) excdir(itmp,iexc)=2 !Negative contribution is assumed to be de-excitation (of course this is not strict since -> and <- have been combined together)
                exccoeff(itmp,iexc)=dsqrt(abs(tmpval))
            end if
            !Although for closed-shell ground state, ORCA still outputs coefficients as normalization to 100%, &
            !However, in order to follow the Gaussian convention, we change the coefficient as normalization to 50%
            if (wfntype==0.or.wfntype==3) exccoeff(itmp,iexc)=exccoeff(itmp,iexc)/dsqrt(2D0)
        end do
    end do
else !Gaussian output or plain text file
    do iexc=1,nstates
        call loclabel(10,"Excited State",ifound,0)
        read(10,*)
        do itmp=1,nexcitorb(iexc)
            read(10,"(a)") c80tmp
            excdir(itmp,iexc)=1 !means ->
            if (index(c80tmp,'<-')/=0) excdir(itmp,iexc)=2 !means <-
            do isign=1,80 !Find position of <- or ->
                if (c80tmp(isign:isign)=='-'.or.c80tmp(isign:isign)=='<') exit
            end do
            !Process left side of <- or ->
            read(c80tmp(:isign-1),"(a)") leftstr
            ilefttype=0 !close
            if (index(leftstr,'A')/=0) ilefttype=1 !Alpha
            if (index(leftstr,'B')/=0) ilefttype=2 !Beta
            if (ilefttype==0) then
                read(leftstr,*) orbleft(itmp,iexc)
            else
                read(leftstr(:len_trim(leftstr)-1),*) orbleft(itmp,iexc)
                if (ilefttype==2) orbleft(itmp,iexc)=orbleft(itmp,iexc)+nbasis
            end if
            !Process right side of <- or ->
            read(c80tmp(isign+2:),"(a)") rightstr
            irighttype=0 !close
            if (index(rightstr,'A')/=0) irighttype=1 !Alpha
            if (index(rightstr,'B')/=0) irighttype=2 !Beta
            if (irighttype==0) then
                read(rightstr,*) orbright(itmp,iexc),exccoeff(itmp,iexc)
            else
                do isplit=1,80
                    if (rightstr(isplit:isplit)=='A'.or.rightstr(isplit:isplit)=='B') exit
                end do
                read(rightstr(:isplit-1),*) orbright(itmp,iexc)
                read(rightstr(isplit+1:),*) exccoeff(itmp,iexc)
                if (irighttype==2) orbright(itmp,iexc)=orbright(itmp,iexc)+nbasis
            end if
        end do
    end do
end if
close(10)

!Test sum of square of the coefficients
write(*,*) "Summary:"
write(*,*) "Exc.state#     Exc.energy(eV)     Multi.   N_pairs    Sum coeff.^2"
do iexc=1,nstates
    sumsqrexc=0
    sumsqrdeexc=0
    do itmp=1,nexcitorb(iexc)
        if (excdir(itmp,iexc)==1) sumsqrexc=sumsqrexc+exccoeff(itmp,iexc)**2
        if (excdir(itmp,iexc)==2) sumsqrdeexc=sumsqrdeexc-exccoeff(itmp,iexc)**2
    end do
    sumsqrall=sumsqrexc+sumsqrdeexc
    write(*,"(i8,f18.5,4x,i8,i11,f16.6)") iexc,excenergy(iexc),excmulti(iexc),nexcitorb(iexc),sumsqrall
end do
write(*,*)
write(*,*) "Please select the destination of the output:"
write(*,*) "1 Output transition dipole moments to screen"
write(*,*) "2 Output transition dipole moments to transdipmom.txt in current folder"
write(*,*) "3 Generate input file of SOS module of Multiwfn as SOS.txt in current folder"
read(*,*) isel
write(*,*)
write(*,*) "Stage 1: Calculating dipole moment integrals between all GTFs..."
nsize=nprims*(nprims+1)/2
allocate(GTFdipint(3,nsize))
call genGTFDmat(GTFdipint,nsize)

call walltime(iwalltime1)
write(*,*) "Stage 2: Calculating dipole moment integrals between all MOs..."
allocate(MOdipint(3,nmo,nmo))
!MOdipint will record dipole moment integrals between all MOs, including all occ+vir alpha and occ+vir beta ones
iprog=0
nthreads=getNThreads()
!$OMP PARALLEL DO SHARED(MOdipint,MOdipintb,iprog) PRIVATE(imo,jmo,iGTF,jGTF,ides,tmpvec) schedule(dynamic) NUM_THREADS(nthreads)
do imo=1,nmo
    do jmo=imo,nmo
        tmpvec=0
        do iGTF=1,nprims
            do jGTF=1,nprims
                if (iGTF>=jGTF) then
                    ides=iGTF*(iGTF-1)/2+jGTF
                else
                    ides=jGTF*(jGTF-1)/2+iGTF
                end if
                tmpvec=tmpvec+co(imo,iGTF)*co(jmo,jGTF)*GTFdipint(:,ides)
            end do
        end do
        MOdipint(:,imo,jmo)=tmpvec
    end do
    if (nprims>300) then
!$OMP CRITICAL
        iprog=iprog+1
        write(*,"(' Finished:',i6,'  /',i6)") iprog,nmo
!$OMP END CRITICAL
    end if
end do
!$OMP END PARALLEL DO
!Fill lower triangle part
do imo=1,nmo
    do jmo=imo+1,nmo
        MOdipint(:,jmo,imo)=MOdipint(:,imo,jmo)
    end do
end do
call walltime(iwalltime2)
write(*,"(' (Stage 2 took up wall clock time',i10,'s)')") iwalltime2-iwalltime1

if (isel==1) then
    iout=6
else if (isel==2) then
    iout=10
    open(iout,file="transdipmom.txt",status="replace")
else if (isel==3) then
    iout=10
    open(iout,file="SOS.txt",status="replace")
end if
fac=1
if (wfntype==0.or.wfntype==3) fac=2
grounddip=0
do imo=1,nmo
    grounddip=grounddip+MOocc(imo)*MOdipint(:,imo,imo)
end do

!Output index, excitation energies and transition dipole moment between ground state and excited states, the format is in line with SOS module
if (isel==3) then
    write(iout,*) nstates
    do i=1,nstates
        write(iout,"(i6,f12.6)") i,excenergy(i)
    end do
    do iexc=0,nstates
        if (iexc==0) then
            write(iout,"(2i6,3(1PE15.6))") 0,iexc,grounddip
        else
            tdvec=0
            do ipair=1,nexcitorb(iexc)
                imo=orbleft(ipair,iexc)
                lmo=orbright(ipair,iexc)
                wei=exccoeff(ipair,iexc)
                tdvec=tdvec+wei*MOdipint(:,imo,lmo)
            end do
            write(iout,"(2i6,3(1PE15.6))") 0,iexc,tdvec*fac
        end if
    end do
end if

write(*,*) "Stage 3: Calculating transition dipole moment between excited states..."
call walltime(iwalltime1)
iprog=0
nthreads=getNThreads()
!$OMP PARALLEL DO SHARED(tdvecmat,iprog) PRIVATE(iexc,jexc,tdvec,imo,lmo,jmo,kmo,wei) schedule(dynamic) NUM_THREADS(nthreads)
do iexc=1,nstates
    do jexc=iexc,nstates
        tdvec=0
        do ipair=1,nexcitorb(iexc)
            imo=orbleft(ipair,iexc)
            lmo=orbright(ipair,iexc)
            do jpair=1,nexcitorb(jexc)
                jmo=orbleft(jpair,jexc)
                kmo=orbright(jpair,jexc)
                wei=exccoeff(ipair,iexc)*exccoeff(jpair,jexc)
                if (excdir(ipair,iexc)==excdir(jpair,jexc)) then !If don't apply this condition, for TD may be evidently deviate to actual result
                    if (excdir(ipair,iexc)==1) then ! -> case
                        if (imo==jmo.and.lmo/=kmo) then
                            tdvec=tdvec+wei*MOdipint(:,lmo,kmo)
                        else if (imo/=jmo.and.lmo==kmo) then
                            tdvec=tdvec-wei*MOdipint(:,imo,jmo)
                        else if (imo==jmo.and.lmo==kmo) then
                            tdvec=tdvec+wei*(grounddip-MOdipint(:,imo,imo)+MOdipint(:,lmo,lmo))
                        end if
!                     else ! <- case. Frankly speaking, I don't know how to exactly deal with this circumstance.&
                    !So I simply ignore them and I found this is the best way currently. The error must be very small, since de-excitation coefficients are often quite small
!                         if (imo==jmo.and.lmo/=kmo) then
!                             tdvec=tdvec-wei*MOdipint(:,lmo,kmo)
!                         else if (imo/=jmo.and.lmo==kmo) then
!                             tdvec=tdvec+wei*MOdipint(:,imo,jmo)
!                         else if (imo==jmo.and.lmo==kmo) then
!                             tdvec=tdvec+wei*(grounddip+MOdipint(:,imo,imo)-MOdipint(:,lmo,lmo))
!                         end if
                    end if
                end if
            end do
        end do
        tdvecmat(:,iexc,jexc)=tdvec*fac
    end do
    
    if (nprims>300) then
!$OMP CRITICAL
        iprog=iprog+1
        write(*,"(' Finished:',i6,'  /',i6)") iprog,nstates
!$OMP END CRITICAL
    end if
end do
!$OMP END PARALLEL DO
call walltime(iwalltime2)
write(*,"(' (Stage 3 took up wall clock time',i10,'s)',/)") iwalltime2-iwalltime1

if (isel==1.or.isel==2) then
    write(iout,"(' Ground state dipole moment in X,Y,Z:',3f12.6,' a,u,',/)") grounddip
    write(iout,"(' Transition dipole moment between excited states (a.u.):')")
    write(iout,*) "    i     j         X             Y             Z        Diff.(eV)   Oscil.str"
end if

do iexc=1,nstates
    do jexc=iexc,nstates
        if (isel==1.or.isel==2) then
            enediff=abs(excenergy(jexc)-excenergy(iexc))
            oscillstr=2D0/3D0*enediff/au2eV*sum(tdvecmat(:,iexc,jexc)**2)
            write(iout,"(2i6,3f14.7,2f12.5)") iexc,jexc,tdvecmat(:,iexc,jexc),enediff,oscillstr
        else if (isel==3) then
            write(iout,"(2i6,3(1PE15.6))") iexc,jexc,tdvecmat(:,iexc,jexc)
        end if
    end do
end do

if (isel==2) then
    close(iout)
    write(*,*) "Done! The result has been outputted to transdipmom.txt in current folder"
else if (isel==3) then
    close(iout)
    write(*,*) "Done! The result has been outputted to SOS.txt in current folder"
end if
write(*,*)
end subroutine




!!!------------- Analyze charge transfer
!See J. Chem. Theory Comput., 7, 2498
subroutine CTanalyze
use defvar
implicit real*8(a-h,o-z)
real*8,allocatable :: Cpos(:,:,:),Cneg(:,:,:),tmpmat(:,:,:)
sumpos=0D0
sumneg=0D0
Rxpos=0D0
Rypos=0D0
Rzpos=0D0
Rxneg=0D0
Ryneg=0D0
Rzneg=0D0
do i=1,nx
    do j=1,ny
        do k=1,nz
            if (cubmat(i,j,k)>0) then
                sumpos=sumpos+cubmat(i,j,k)
                Rxpos=Rxpos+cubmat(i,j,k)*(orgx+(i-1)*dx)
                Rypos=Rypos+cubmat(i,j,k)*(orgy+(j-1)*dy)
                Rzpos=Rzpos+cubmat(i,j,k)*(orgz+(k-1)*dz)
            else
                sumneg=sumneg+cubmat(i,j,k)
                Rxneg=Rxneg+cubmat(i,j,k)*(orgx+(i-1)*dx)
                Ryneg=Ryneg+cubmat(i,j,k)*(orgy+(j-1)*dy)
                Rzneg=Rzneg+cubmat(i,j,k)*(orgz+(k-1)*dz)
            end if
        end do
    end do
end do
dvol=dx*dy*dz
sumpos=sumpos*dvol
sumneg=sumneg*dvol
Rxpos=Rxpos*dvol/sumpos
Rypos=Rypos*dvol/sumpos
Rzpos=Rzpos*dvol/sumpos
Rxneg=Rxneg*dvol/sumneg
Ryneg=Ryneg*dvol/sumneg
Rzneg=Rzneg*dvol/sumneg
disx=abs(Rxpos-Rxneg)
disy=abs(Rypos-Ryneg)
disz=abs(Rzpos-Rzneg)
disnorm=dsqrt(disx**2+disy**2+disz**2)
write(*,"(' Transferred charge (positive and negative parts):',2f8.3)") sumpos,sumneg
write(*,"(' Barycenter of positive part in x,y,z (Angstrom):',3f8.3)") Rxpos*b2a,Rypos*b2a,Rzpos*b2a
write(*,"(' Barycenter of negative part in x,y,z (Angstrom):',3f8.3)") Rxneg*b2a,Ryneg*b2a,Rzneg*b2a
write(*,"(' Distance of CT in x,y,z (Angstrom):',3f8.3,' Norm:',f8.3)") disx*b2a,disy*b2a,disz*b2a,disnorm*b2a
dipx=-(Rxpos-Rxneg)*sumpos
dipy=-(Rypos-Ryneg)*sumpos
dipz=-(Rzpos-Rzneg)*sumpos
dipnorm=disnorm*sumpos
write(*,"(' Dipole moment variation (a.u.) :',3f8.3,' Norm:',f8.3)") dipx,dipy,dipz,dipnorm
write(*,"(' Dipole moment variation (Debye):',3f8.3,' Norm:',f8.3)") dipx*au2debye,dipy*au2debye,dipz*au2debye,dipnorm*au2debye

sigxpos=0D0
sigypos=0D0
sigzpos=0D0
sigxneg=0D0
sigyneg=0D0
sigzneg=0D0
do i=1,nx
    do j=1,ny
        do k=1,nz
            if (cubmat(i,j,k)>0) then
                sigxpos=sigxpos+cubmat(i,j,k)*((orgx+(i-1)*dx)-Rxpos)**2
                sigypos=sigypos+cubmat(i,j,k)*((orgy+(j-1)*dy)-Rypos)**2
                sigzpos=sigzpos+cubmat(i,j,k)*((orgz+(k-1)*dz)-Rzpos)**2
            else
                sigxneg=sigxneg+cubmat(i,j,k)*((orgx+(i-1)*dx)-Rxneg)**2
                sigyneg=sigyneg+cubmat(i,j,k)*((orgy+(j-1)*dy)-Ryneg)**2
                sigzneg=sigzneg+cubmat(i,j,k)*((orgz+(k-1)*dz)-Rzneg)**2
            end if
        end do
    end do
end do
sigxpos=dsqrt(sigxpos/(sumpos/dvol))
sigypos=dsqrt(sigypos/(sumpos/dvol))
sigzpos=dsqrt(sigzpos/(sumpos/dvol))
signormpos=dsqrt(sigxpos**2+sigypos**2+sigzpos**2)
sigxneg=dsqrt(sigxneg/(sumneg/dvol))
sigyneg=dsqrt(sigyneg/(sumneg/dvol))
sigzneg=dsqrt(sigzneg/(sumneg/dvol))
signormneg=dsqrt(sigxneg**2+sigyneg**2+sigzneg**2)
write(*,"(' RMSD of positive part in x,y,z (Angstrom):',3f8.3,' Tot:',f7.3)") sigxpos*b2a,sigypos*b2a,sigzpos*b2a,signormpos*b2a
write(*,"(' RMSD of negative part in x,y,z (Angstrom):',3f8.3,' Tot:',f7.3)") sigxneg*b2a,sigyneg*b2a,sigzneg*b2a,signormneg*b2a
delta_sigCT=dsqrt((sigxpos-sigxneg)**2+(sigypos-sigyneg)**2+(sigzpos-sigzneg)**2)
write(*,"(' Difference between RMSD of positive and negative parts (Angstrom):',f8.3)") delta_sigCT*b2a
Hx=(sigxpos+sigxneg)/2D0
Hy=(sigypos+sigyneg)/2D0
Hz=(sigzpos+sigzneg)/2D0
Hnorm=dsqrt(Hx**2+Hy**2+Hz**2)
write(*,"(' H index in x,y,z (Angstrom):',3f8.3,' Norm:',f8.3)") Hx*b2a,Hy*b2a,Hz*b2a,Hnorm*b2a
tx=disx-Hx
ty=disy-Hy
tz=disz-Hz
tnorm=dsqrt(tx**2+ty**2+tz**2)
write(*,"(' t index in x,y,z (Angstrom):',3f8.3,' Norm:',f8.3)") tx*b2a,ty*b2a,tz*b2a,tnorm*b2a

allocate(Cpos(nx,ny,nz),Cneg(nx,ny,nz),tmpmat(nx,ny,nz))
do i=1,nx
    do j=1,ny
        do k=1,nz
            rnowx=orgx+(i-1)*dx
            rnowy=orgy+(j-1)*dy
            rnowz=orgz+(k-1)*dz
            Cpos(i,j,k)=exp( -(rnowx-Rxpos)**2/(2*sigxpos**2) -(rnowy-Rypos)**2/(2*sigypos**2) -(rnowz-Rzpos)**2/(2*sigzpos**2))
            Cneg(i,j,k)=exp( -(rnowx-Rxneg)**2/(2*sigxneg**2) -(rnowy-Ryneg)**2/(2*sigyneg**2) -(rnowz-Rzneg)**2/(2*sigzneg**2))
        end do
    end do
end do
Cpos=Cpos*sumpos/(sum(Cpos)*dvol)
Cneg=Cneg*sumneg/(sum(Cneg)*dvol)

write(*,"(' Overlap integral between C+ and C-:',f12.6)") sum(dsqrt(Cpos/sumpos)*dsqrt(Cneg/sumneg))*dvol

sur_value=0.001
do while(.true.)
    write(*,*)
    write(*,*) "0 Return"
    write(*,*) "1 Show isosurface of C+ and C- functions simultaneously"
    write(*,*) "2 Export C+ and C- functions to cube file in current folder"
    read(*,*) isel

    if (isel==0) then
        return
    else if (isel==1) then
        isosursec=1
        clrRcub1sameold=clrRcub1same !Backup previous color setting
        clrGcub1sameold=clrGcub1same
        clrBcub1sameold=clrBcub1same
        clrRcub2oppoold=clrRcub2oppo
        clrGcub2oppoold=clrGcub2oppo
        clrBcub2oppoold=clrBcub2oppo
        clrRcub2samemeshptold=clrRcub2samemeshpt
        clrGcub2samemeshptold=clrGcub2samemeshpt
        clrBcub2samemeshptold=clrBcub2samemeshpt
        clrRcub2oppomeshptold=clrRcub2oppomeshpt
        clrGcub2oppomeshptold=clrGcub2oppomeshpt
        clrBcub2oppomeshptold=clrBcub2oppomeshpt
        clrRcub1same=0.3D0 !Set color to that Cpos is green, Cneg is blue. (Cpos/Cneg function is positive/negative everywhere due to sumpos and sumneg)
        clrGcub1same=0.75D0
        clrBcub1same=0.3D0
        clrRcub2oppo=0.3D0
        clrGcub2oppo=0.45D0
        clrBcub2oppo=0.9D0
        clrRcub2samemeshpt=0.3D0
        clrGcub2samemeshpt=0.75D0
        clrBcub2samemeshpt=0.3D0
        clrRcub2oppomeshpt=0.3D0
        clrGcub2oppomeshpt=0.45D0
        clrBcub2oppomeshpt=0.9D0
        if (allocated(cubmattmp)) deallocate(cubmattmp)
        allocate(cubmattmp(nx,ny,nz))
        tmpmat=cubmat
        cubmat=Cpos
        cubmattmp=Cneg
        !Since drawisosurgui is mainly designed for plotting one isosurface, while here we use it to draw both cubmat and cubmattmp, so disable change of style to avoid troubles
        cubmat=tmpmat !Recover original data, namely density difference between two states
        deallocate(cubmattmp)
        clrRcub1same=clrRcub1sameold !Recover previous color setting
        clrGcub1same=clrGcub1sameold
        clrBcub1same=clrBcub1sameold
        clrRcub2oppo=clrRcub2oppoold
        clrGcub2oppo=clrGcub2oppoold
        clrBcub2oppo=clrBcub2oppoold
        clrRcub2samemeshpt=clrRcub2samemeshptold
        clrGcub2samemeshpt=clrGcub2samemeshptold
        clrBcub2samemeshpt=clrBcub2samemeshptold
        clrRcub2oppomeshpt=clrRcub2oppomeshptold
        clrGcub2oppomeshpt=clrGcub2oppomeshptold
        clrBcub2oppomeshpt=clrBcub2oppomeshptold
    else if (isel==2) then
        open(10,file="Cpos.cub",status="replace")
        call outcube(Cpos,nx,ny,nz,orgx,orgy,orgz,dx,dy,dz,10)
        close(10)
        write(*,"('C+ function has been outputted to ""Cpos.cub"" in current folder')")
        open(10,file="Cneg.cub",status="replace")
        call outcube(Cneg,nx,ny,nz,orgx,orgy,orgz,dx,dy,dz,10)
        close(10)
        write(*,"('C- function has been outputted to ""Cneg.cub"" in current folder')")
    end if
end do
end subroutine





!--------------------------------------------------------------------------
!----------- Plot transition density matrix as color-filled map -----------
!--------------------------------------------------------------------------
subroutine plottransdensmat
use defvar
use util
implicit real*8 (a-h,o-z)
character tdmatfilename*200
real*8 tdmatbas(nbasis,nbasis),tdmatatmtmp(ncenter,ncenter),tdmatatm(ncenter,ncenter)
real*8,allocatable :: tdmatatmnoh(:,:)
!Load density transition matrix in basis expansion
tdmatbas=0D0
write(*,"(a)") " Input the path of the Gaussian output file or plain text file containing transition density matrix, e.g. c:\a.out"
do while(.true.)
    read(*,"(a)") tdmatfilename
    inquire(file=tdmatfilename,exist=alive)
    if (alive) exit
    write(*,*) "Error: Cannot find the file, please input again"
end do

open(10,file=tdmatfilename,status="old")
if (index(tdmatfilename,"AAtrdip.txt")/=0) then !Directly load and use atom-atom matrix (commonly generated by option 5 of hole-electron module)
    write(*,*) "Loading atom-atom contribution matrix..."
    call readmatgau(10,tdmatatm,0,"f14.8",6,5,1)
else !Load basis-basis matrix and then condense to atom-atom matrix
    call loclabel(10,"Gaussian",igauout,maxline=100)
    rewind(10)
    if (igauout==1) then !Gaussian output file
        write(*,*) "This is a Gaussian output file"
        call loclabel(10,"Alpha Density Matrix:",ifound)
        if (ifound==1) then !Open-shell
            write(*,*) "Use which type of transition density matrix?"
            write(*,*) "1=Alpha    2=Beta"
            read(*,*) iTDMtype
            write(*,*) "Loading transition density matrix..."
            if (iTDMtype==1) then
                call readmatgau(10,tdmatbas,1,"f10.5",21,5,1)
            else if  (iTDMtype==2) then
                call loclabel(10,"Beta Density Matrix:",ifound)
                call readmatgau(10,tdmatbas,1,"f10.5",21,5,1)
            end if
        else !Close-shell
            call loclabel(10,"Density Matrix:",ifound)
            if (ifound==0) call loclabel(10,"DENSITY MATRIX.",ifound)
            if (ifound==0) then
                write(*,"(a,/)") "Error: Cannot found transition density matrix information from the Gaussian output file, please check if iop(6/8=3) has been specified"
                return
            end if
            write(*,*) "Loading transition density matrix..."
            call readmatgau(10,tdmatbas,1,"f10.5",21,5,1)
        end if
    else !Plain text file with transition density matrix in basis function
        write(*,*) "Loading transition density matrix..."
        call readmatgau(10,tdmatbas,0,"f14.8",6,5,1)
    end if
    !Use tdmatbas to construct tdmatatm
    tdmatatm=0D0
    do iatm=1,ncenter
        do jatm=1,ncenter
            do ibas=basstart(iatm),basend(iatm)
                do jbas=basstart(jatm),basend(jatm)
                    tdmatatm(iatm,jatm)=tdmatatm(iatm,jatm)+tdmatbas(ibas,jbas)**2
    !                 tdmatatm(iatm,jatm)=tdmatatm(iatm,jatm)+tdmatbas(ibas,jbas)
    !                 tdmatatm(iatm,jatm)=tdmatatm(iatm,jatm)+abs(tdmatbas(ibas,jbas))
                end do
            end do
        end do
    end do
end if
close(10)

tdmatatmtmp=tdmatatm
!Contract tdmatatm to tdmatatmnoh, tdmatatmtmp is used as a intermediate
!tdmatatmnoh is the tdmatatm without hydrogens
ncenreal=count(a(:)%index/=1)
write(*,"(' The number of non-hydrogen atoms:',i10)") ncenreal
allocate(tdmatatmnoh(ncenreal,ncenreal))
itmp=0
do iatm=1,ncenter
    if (a(iatm)%index/=1) then
        itmp=itmp+1
        tdmatatmtmp(:,itmp)=tdmatatmtmp(:,iatm)
        tdmatatmtmp(itmp,:)=tdmatatmtmp(iatm,:)
    end if
end do
tdmatatmnoh(:,:)=tdmatatmtmp(1:ncenreal,1:ncenreal)

ifhydrogen=0
clrlimlow=minval(tdmatatmnoh)
clrlimhigh=maxval(tdmatatmnoh)
ninterpo=10
nstepsize=5
ifnormsum=0
! clrlimlow=minval(tdmatbas)
! clrlimhigh=maxval(tdmatbas)
! write(*,*) clrlimlow,clrlimhigh
facnorm=1D0 !Normalization factor, default is 1, namely don't do normalization

do while(.true.)
    write(*,*)
    write(*,*) "0 Return"
    write(*,*) "1 Show transition density matrix map"
    write(*,*) "2 Save transition density matrix map to a graphical file in current folder"
    write(*,*) "3 Export data to plain text file"
    if (ifhydrogen==0) write(*,*) "4 Switch if take hydrogens into account, current: No"
    if (ifhydrogen==1) write(*,*) "4 Switch if take hydrogens into account, current: Yes"
    write(*,"(a,f7.4,a,f7.4)") " 5 Change lower and upper limit of color scale, current:",clrlimlow," to",clrlimhigh
    write(*,"(a,i3)") " 6 Set the number of interpolation steps between grid, current:",ninterpo
    write(*,"(a,i3)") " 7 Set stepsize between labels, current:",nstepsize
    if (ifnormsum==0) write(*,*) "8 Switch if normalized the sum of all elements to unity, current: No"
    if (ifnormsum==1) write(*,*) "8 Switch if normalized the sum of all elements to unity, current: Yes"
    read(*,*) isel

    if (isel==0) then
        return
    else if (isel==1.or.isel==2) then
        if (isel==2) isavepic=1
        if (ifhydrogen==1) then
            nspc=(ncenter-1)/nlabstep
        else if (ifhydrogen==0) then
            nspc=(ncenreal-1)/nlabstep
        end if
        if (isel==2) then
            isavepic=0
            write(*,*) "Done, the picture has been saved to current folder with ""DISLIN"" prefix"
        end if
    else if (isel==3) then
        open(10,file="tdmat.txt",status="replace")
        if (ifhydrogen==0) then 
            do iatm=1,ncenreal
                do jatm=1,ncenreal
                    write(10,"(2i8,f12.6)") iatm,jatm,tdmatatmnoh(iatm,jatm)/facnorm
                end do
            end do
        else if (ifhydrogen==1) then
            do iatm=1,ncenter
                do jatm=1,ncenter
                    write(10,"(2i8,f12.6)") iatm,jatm,tdmatatm(iatm,jatm)/facnorm
                end do
            end do
        end if
        close(10)
        write(*,*) "Done, the data have been exported to tdmat.txt in current folder"
    else if (isel==4) then
        if (ifhydrogen==0) then
            ifhydrogen=1
            clrlimlow=minval(tdmatatm)
            clrlimhigh=maxval(tdmatatm)
        else if (ifhydrogen==1) then
            ifhydrogen=0
            clrlimlow=minval(tdmatatmnoh)
            clrlimhigh=maxval(tdmatatmnoh)
        end if
    else if (isel==5) then
        write(*,*) "Input lower and upper limits, e.g. 0,1.5"
        read(*,*) clrlimlow,clrlimhigh
    else if (isel==6) then
        write(*,*) "Please input a number"
        write(*,"(a)") " Note: Larger value gives rise to more smooth graph, 1 means don't do interpolation"
        read(*,*) ninterpo
    else if (isel==7) then
        write(*,*) "Please input a number, e.g. 5"
        read(*,*) nstepsize
    else if (isel==8) then
        if (ifnormsum==0) then
            ifnormsum=1
            facnorm=sum(tdmatatm(:,:)) !The sum of all elements
            write(*,*) "The normalization factor is", facnorm
        else if (ifnormsum==1) then
            ifnormsum=0
            facnorm=1D0
        end if
    end if
end do
end subroutine





!---------------------------------------------------------------------------------------
!---------- Calculate interfragment charger transfer in electronic excitation ----------
!---------------------------------------------------------------------------------------
subroutine excfragCT(nexcitorb,orbleft,orbright,excdir,exccoeff)
use defvar
use util
implicit real*8 (a-h,o-z)
integer nexcitorb,orbleft(nexcitorb),orbright(nexcitorb),excdir(nexcitorb)
real*8 exccoeff(nexcitorb),CTmat(ncenter,ncenter),CTmattmp(ncenter,ncenter),atmcomp(ncenter,nmo)
character c2000tmp*2000

inquire(file="orbcomp.txt",exist=alive)
if (alive) then !Load atomic contribution from orbcomp.txt, which may be outputted by option -4 of Hirshfeld/Becke composition analysis
    write(*,"(a)") " orbcomp.txt was found in current folder, now load atomic contribution to all orbitals from this file..."
    open(10,file="orbcomp.txt",status="old")
    do imo=1,nmo
        read(10,*)
        do iatm=1,ncenter
            read(10,*) inouse,atmcomp(iatm,imo)
        end do
    end do
    close(10)
    atmcomp=atmcomp/100
else
    write(*,*) "Calculating atomic contribution to all orbitals by SCPA method..."
    do imo=1,nmo
        if (MOtype(imo)==0.or.MOtype(imo)==1) then !Close-shell or alpha part of open-shell
            sumsqr=sum(cobasa(:,imo)**2)
            do iatm=1,ncenter
                atmcomp(iatm,imo)=sum(cobasa(basstart(iatm):basend(iatm),imo)**2)/sumsqr
            end do
        else !Beta part of open-shell
            iimo=imo-nbasis
            sumsqr=sum(cobasb(:,iimo)**2)
            do iatm=1,ncenter
                atmcomp(iatm,imo)=sum(cobasb(basstart(iatm):basend(iatm),iimo)**2)/sumsqr
            end do
        end if
    end do
end if

write(*,*) "Constructing inter-atomic CT matrix..."
CTmat=0
fac=1
if (wfntype==0.or.wfntype==3) fac=2 !Since for closed-shell case, program only print one-half part
do iexc=1,nexcitorb
    imo=orbleft(iexc)
    jmo=orbright(iexc)
    !Calculate transferred electron from iatm(row) to jatm(column)
    do iatm=1,ncenter
        do jatm=1,ncenter
            CTmattmp(iatm,jatm)=atmcomp(iatm,imo)*atmcomp(jatm,jmo)
        end do
    end do
    if (excdir(iexc)/=1) CTmattmp=-CTmattmp
    CTmat=CTmat+CTmattmp*fac*exccoeff(iexc)**2
end do
write(*,*) "Preparation is finished!"

do while(.true.)
    write(*,*)
    write(*,*) "Input atom list for fragment 1, e.g. 3,5-8,15-20"
    write(*,*) "Note: Input 0 can exit"
    read(*,"(a)") c2000tmp
    if (c2000tmp(1:1)=="0") return
    call str2arr(c2000tmp,nterm1)
    if (allocated(frag1)) deallocate(frag1)
    allocate(frag1(nterm1))
    call str2arr(c2000tmp,nterm1,frag1)
    write(*,*) "Input atom list for fragment 2, e.g. 1,2,4,9-14"
    read(*,"(a)") c2000tmp
    call str2arr(c2000tmp,nterm2)
    if (allocated(frag2)) deallocate(frag2)
    allocate(frag2(nterm2))
    call str2arr(c2000tmp,nterm2,frag2)

    CTval1=0
    CTval2=0
    varpop1=0
    varpop2=0
    do idx=1,nterm1
        iatm=frag1(idx)
        do jdx=1,nterm2
            jatm=frag2(jdx)
            CTval1=CTval1+CTmat(iatm,jatm)
            CTval2=CTval2+CTmat(jatm,iatm)
        end do
        varpop1=varpop1+sum(CTmat(:,iatm))-sum(CTmat(iatm,:))
    end do
    do jdx=1,nterm2
        jatm=frag2(jdx)
        varpop2=varpop2+sum(CTmat(:,jatm))-sum(CTmat(jatm,:))
    end do
    write(*,"(a,f10.5)") " Electron transferred from fragment 1 to 2:",CTval1
    write(*,"(a,f10.5)") " Electron transferred from fragment 2 to 1:",CTval2
    write(*,"(a,f10.5)") " Electron net transferred from fragment 1 to 2:",CTval1-CTval2
    write(*,"(a,f10.5)") " Variation of electron population of fragment 1:",varpop1
    write(*,"(a,f10.5)") " Variation of electron population of fragment 2:",varpop2
end do

end subroutine
