subroutine bondana
use defvar
use util
implicit real*8 (a-h,o-z)
character c2000tmp*2000
do while(.true.)
    write(*,*) "           ================ Bond order analysis ==============="
    if (ifiletype==1.or.ifiletype==9) then
        if (allocated(frag1)) then
            write(*,*) "-1 Redefine fragment 1 and 2 for option 1,3,4,7,8"
        else
            write(*,*) "-1 Define fragment 1 and 2 for option 1,3,4,7,8 (to be defined)"
        end if
    end if
    write(*,*) "0 Return"
    write(*,*) "1 Mayer bond order analysis"
    write(*,*) "2 Multicenter bond order analysis (up to 12 centers)"
    write(*,*) "-2 Multicenter bond order analysis in NAO basis (up to 10 centers)"
!     write(*,*) "22 Multicenter bond order analysis in Lowdin orthogonalized basis (up to 10 centers)"
    write(*,*) "3 Wiberg bond order analysis in Lowdin orthogonalized basis"
    write(*,*) "4 Mulliken bond order analysis"
    write(*,*) "5 Decompose Mulliken bond order between two atoms to orbital contributions"
    write(*,*) "6 Orbital occupancy-perturbed Mayer bond order"
    write(*,*) "7 Fuzzy bond order analysis"
    write(*,*) "8 Laplacian bond order"
    read(*,*) ibondana
    if (.not.allocated(shtype).and.(ibondana>=1.and.ibondana<=6)) then
        write(*,"(a)") " ERROR: The input file you used does not contain basis function information! Please check Section 2.5 of the manual for explanation"
        return
    else if (.not.allocated(b).and.(ibondana==7.or.ibondana==8)) then
        write(*,"(a)") " ERROR: The input file you used does not contain GTF information! Please check Section 2.5 of the manual for explanation"
        return
    end if
    
    if (ibondana==-1) then
        !Define frag1, the size just accomodates content
        if (allocated(frag1)) then
            write(*,*) "Atoms in current fragment 1:"
            write(*,"(13i6)") frag1
            write(*,"(a)") " Input 0 to keep unchanged, or redefine fragment 1, e.g. 1,3-6,8,10-11 means the atoms 1,3,4,5,6,8,10,11 will constitute the fragment 1"
        else
            write(*,"(a)") " Input atomic indices to define fragment 1. e.g. 1,3-6,8,10-11 means the atoms 1,3,4,5,6,8,10,11 will constitute the fragment 1"
        end if
        read(*,"(a)") c2000tmp
        if (c2000tmp(1:1)/='0') then
            if (allocated(frag1)) deallocate(frag1)
            call str2arr(c2000tmp,nfragatm)
            allocate(frag1(nfragatm))
            call str2arr(c2000tmp,nfragatm,frag1)
        end if
        !Define frag2, the size just accomodates content
        if (allocated(frag2)) then
            write(*,*) "Atoms in current fragment 2:"
            write(*,"(13i6)") frag2
            write(*,"(a)") " Input 0 to keep unchanged, or redefine fragment 2, e.g. 1,3-6,8,10-11 means the atoms 1,3,4,5,6,8,10,11 will constitute fragment 2"
        else
            write(*,"(a)") " Input atomic indices to define fragment 2. e.g. 1,3-6,8,10-11 means the atoms 1,3,4,5,6,8,10,11 will constitute the fragment 2"
        end if
        read(*,"(a)") c2000tmp
        if (c2000tmp(1:1)/='0') then
            if (allocated(frag2)) deallocate(frag2)
            call str2arr(c2000tmp,nfragatm)
            allocate(frag2(nfragatm))
            call str2arr(c2000tmp,nfragatm,frag2)
        end if
        if (any(frag1>ncenter).or.any(frag2>ncenter)) then
            write(*,*) "Error: Some atomic indices exceeded valid range! Please define again"
            write(*,*)
            deallocate(frag1,frag2)
            cycle
        end if
        write(*,*) "Setting is saved"
        write(*,*) "Now the atoms in fragment 1 are"
        write(*,"(13i6)") frag1
        write(*,*) "Now the atoms in fragment 2 are"
        write(*,"(13i6)") frag2
        if (any(frag1<=0).or.any(frag1>ncenter).or.any(frag2<=0).or.any(frag2>ncenter)) write(*,*) "Warning: Indices of some atoms exceed valid range! Please redefine fragment"
        do i=1,size(frag1)
            if (any(frag2==frag1(i))) then
                write(*,"(a)") "Warning: Indices of some atoms are duplicated in the two fragments! Please redefine them"
                exit
            end if
        end do
        write(*,*)
        
    else if (ibondana==0) then
        if (allocated(frag1)) deallocate(frag1)
        if (allocated(frag2)) deallocate(frag2)
        exit
    else if (ibondana==1) then
        write(*,*) "Please wait..."
        call mayerbndord
    else if (ibondana==2) then
        call multicenter
    else if (ibondana==-2) then
        call multicenterNAO
    else if (ibondana==22) then
        write(*,*) "Performing Lowdin orthogonalization..."
         call symmortho
        call multicenter
    else if (ibondana==3) then
        !In symmortho the density matrix, cobasa/b and Sbas will change, so backup them
        if (allocated(Cobasb)) then !Open-shell
            allocate(Cobasa_org(nbasis,nmo/2),Cobasb_org(nbasis,nmo/2),Sbas_org(nbasis,nbasis))
            Cobasb_org=Cobasb
        else
            allocate(Cobasa_org(nbasis,nmo),Sbas_org(nbasis,nbasis)) 
        end if
        Cobasa_org=Cobasa
        Sbas_org=Sbas
        write(*,*) "Performing Lowdin orthogonalization..."
         call symmortho
        write(*,*) "Calculating Wiberg bond order..."
        call mayerbndord
        write(*,*) "Regenerating density matrix..."
        write(*,*)
        Cobasa=Cobasa_org
        Sbas=Sbas_org
        deallocate(Cobasa_org,Sbas_org)
        if (allocated(Cobasb_org)) then
            Cobasb=Cobasb_org
            deallocate(Cobasb_org)
        end if
        call genP
    else if (ibondana==4) then
        call mullikenbndord
    else if (ibondana==5) then
        call decompMBO
    else if (ibondana==6) then
        call OrbPertMayer
    else if (ibondana==7) then
        call intatomspace(1)
    else if (ibondana==8) then
        call intatomspace(2)
    end if
end do
end subroutine


!! ----------------- Mayer/Generalized Wiberg 2-c bond order analysis
! Mayer bond order analysis and Generalized Wiberg bond order (GWBO) analysis
! If Lowdin orthogonalization has been performed, that is carry out Wiberg bond order analysis in Lowdin orthogonalized basis
! Note: For close-shell, two methods give the same result. For open-shell, Mayer bond order for all electrons is the sum of
! alpha and beta bond order, while GWBO directly use total density matrix to generate
! total bond order, the "Mayer bond order" in gaussian is actually GWBO!
subroutine mayerbndord
use defvar
use util
implicit real*8 (a-h,o-z)
real*8 :: bndmata(ncenter,ncenter),bndmatb(ncenter,ncenter),bndmattot(ncenter,ncenter),&
PSmata(nbasis,nbasis),PSmatb(nbasis,nbasis),PSmattot(nbasis,nbasis)
! integer myfrag1(6),myfrag2(6),myfrag3(6) !For my silicon cluster study
character selectyn

bndmata=0D0
bndmatb=0D0
bndmattot=0D0
!Calculate total bond order for restricted close-shell wavefunction (for open-shell do GWBO, P=Palpha+Pbeta)
PSmattot=matmul(Ptot,Sbas)
do i=1,ncenter
    do j=i+1,ncenter
        accum=0D0
        do ii=basstart(i),basend(i)
            do jj=basstart(j),basend(j)
                accum=accum+PSmattot(ii,jj)*PSmattot(jj,ii)
            end do
        end do
        bndmattot(i,j)=accum
    end do
end do
bndmattot=bndmattot+transpose(bndmattot) !Because we only filled one triangular region, copy it to another
do i=1,ncenter
    bndmattot(i,i)=sum(bndmattot(i,:))
end do

if (wfntype==1.or.wfntype==2.or.wfntype==4) then
    PSmata=matmul(Palpha,Sbas)
    PSmatb=matmul(Pbeta,Sbas)
    do i=1,ncenter
        do j=i+1,ncenter
            accuma=0D0
            accumb=0D0
            do ii=basstart(i),basend(i)
                do jj=basstart(j),basend(j)
                    accuma=accuma+PSmata(ii,jj)*PSmata(jj,ii)
                    accumb=accumb+PSmatb(ii,jj)*PSmatb(jj,ii)
                end do
            end do
            bndmata(i,j)=accuma
            bndmatb(i,j)=accumb
        end do
    end do
    bndmata=2*(bndmata+transpose(bndmata))
    bndmatb=2*(bndmatb+transpose(bndmatb))
    do i=1,ncenter
        bndmata(i,i)=sum(bndmata(i,:))
        bndmatb(i,i)=sum(bndmatb(i,:))
    end do
end if

write(*,"(' The total bond order >=',f10.6)") bndordthres
itmp=0
if (wfntype==1.or.wfntype==2.or.wfntype==4) then
    do i=1,ncenter
        do j=i+1,ncenter
            if (bndmata(i,j)+bndmatb(i,j)>=bndordthres) then
                itmp=itmp+1
                write(*,"(' #',i5,':',i5,a,i5,a,' Alpha: ',f10.6,' Beta:',f10.6,' Total:',f10.6)") &
                itmp,i,'('//a(i)%name//')',j,'('//a(j)%name//')',bndmata(i,j),bndmatb(i,j),bndmata(i,j)+bndmatb(i,j)
            end if
        end do
    end do
    write(*,*)
    write(*,"(' Bond order from mixed alpha&beta density matrix >=',f10.6)") bndordthres
end if
itmp=0
do i=1,ncenter
    do j=i+1,ncenter
        if (bndmattot(i,j)>=bndordthres) then
            itmp=itmp+1
            write(*,"(' #',i5,':',5x,i5,a,i5,a,f14.8)") itmp,i,'('//a(i)%name//')',j,'('//a(j)%name//')',bndmattot(i,j)
        end if
    end do
end do

write(*,*)
write(*,*) "Total valences and free valences defined by Mayer:"
do i=1,ncenter
    accum=0D0
    accum2=0D0
    do ii=basstart(i),basend(i)
        accum=accum+2*PSmattot(ii,ii)
        do jj=basstart(i),basend(i)
            accum2=accum2+PSmattot(ii,jj)*PSmattot(jj,ii)
        end do
    end do
    freeval=accum-accum2-(bndmata(i,i)+bndmatb(i,i))
    if (wfntype==0.or.wfntype==3) freeval=0D0
    write(*,"(' Atom',i6,'(',a,') :',2f14.8)") i,a(i)%name,accum-accum2,freeval
end do

!Between fragment
if (allocated(frag1)) then
    bndordfraga=0
    bndordfragb=0
    bndordfragtot=0
    do i=1,size(frag1)
        do j=1,size(frag2)
            bndordfraga=bndordfraga+bndmata(frag1(i),frag2(j))
            bndordfragb=bndordfragb+bndmatb(frag1(i),frag2(j))
            bndordfragtot=bndordfragtot+bndmattot(frag1(i),frag2(j))
        end do
    end do
    write(*,*)
    if (wfntype==1.or.wfntype==2.or.wfntype==4) then
        write(*,"(' The bond order between fragment 1 and 2:')")
        write(*,"(' Alpha:',f10.6,' Beta:',f10.6,' Total:',f10.6,' Mixed Alpha&Beta:',f10.6)") bndordfraga,bndordfragb,bndordfraga+bndordfragb,bndordfragtot
    else
        write(*,"(' The bond order between fragment 1 and 2:',f12.6)") bndordfragtot
    end if
end if
write(*,*)

write(*,*) "If output bond order matrix to bndmat.txt in current folder? (y/n)"
read(*,*) selectyn
if (selectyn=='y'.or.selectyn=='Y') then
    open(10,file="bndmat.txt",status="replace")
    write(10,*) "Note: The diagonal elements are the sum of corresponding row elements"
    if (wfntype==0.or.wfntype==3) then
        call showmatgau(bndmattot,"Bond order matrix",0,"f14.8",10)
    else if (wfntype==1.or.wfntype==2.or.wfntype==4) then
        call showmatgau(bndmata,"Bond order matrix for alpha electrons",0,"f14.8",10)
        call showmatgau(bndmatb,"Bond order matrix for beta electrons",0,"f14.8",10)
        call showmatgau(bndmata+bndmatb,"Bond order matrix for all electrons",0,"f14.8",10)
        call showmatgau(bndmattot,"Bond order matrix from mixed density",0,"f14.8",10)
    end if
    close(10)
    write(*,*) "Result have been outputted to bndmat.txt in current folder"
    write(*,*)
end if
end subroutine



!! ----------- Multi-center bond order analysis
subroutine multicenter
use defvar
use util
implicit real*8 (a-h,o-z)
real*8 :: PSmat(nbasis,nbasis),PSmata(nbasis,nbasis),PSmatb(nbasis,nbasis),PSmattot(nbasis,nbasis),maxbndord
integer cenind(12),maxcenind(12) !maxcenind is used to store the combination that has the maximum bond order in automatical search
integer itype
character c1000tmp*1000
write(*,*) "Please wait..."

ntime=1 !close-shell
PSmattot=matmul(Ptot,Sbas)
if (wfntype==1.or.wfntype==2.or.wfntype==4) then !open-shell or close-shell with symmetry breaking
    ntime=3
    PSmata=matmul(Palpha,Sbas)
    PSmatb=matmul(Pbeta,Sbas)
end if

do while(.true.)
    write(*,*) "Input atom indices, e.g. 3,4,7,8,10    (up to 12 atoms)"
    write(*,*) "Input -3/-4/-5/-6 can search all possible three/four/five/six-center bonds"
    write(*,*) "Input 0 can return to upper level menu"
    read(*,"(a)") c1000tmp

    if (c1000tmp(1:1)=='0') then
        Return
    else if (c1000tmp(1:1)/='-') then
        call str2arr(c1000tmp,nbndcen,cenind)
        
        do itime=1,ntime
            if (wfntype==0.or.wfntype==3) then
                PSmat=PSmattot
            else if (wfntype==1.or.wfntype==2.or.wfntype==4) then
                if (itime==1) PSmat=PSmattot
                if (itime==2) PSmat=PSmata
                if (itime==3) PSmat=PSmatb
            end if
            if (nbndcen>=8) write(*,*) "Please wait..."
            call calcmultibndord(nbndcen,cenind,PSmat,nbasis,accum)
            
            if (itime==1) bndordmix=accum
            if (itime==2) bndordalpha=accum*2**(nbndcen-1)
            if (itime==3) bndordbeta=accum*2**(nbndcen-1)
        end do !end ntime
        if (wfntype==0.or.wfntype==3) then
            write(*,"(a,f16.10)") " The multicenter bond order:",accum
            !Normalized multicenter bond order, see Electronic Aromaticity Index for Large Rings DOI: 10.1039/C6CP00636A
            !When it is negative, first obtain **(1/n) using its absolute value, then multiply it by -1
            write(*,"(a,f16.10)") " The normalized multicenter bond order:",accum/abs(accum) * (abs(accum)**(1D0/nbndcen))
            
        else if (wfntype==1.or.wfntype==2.or.wfntype==4) then
            write(*,"(a,f13.7)") " The multicenter bond order from alpha density matrix:",bndordalpha
            write(*,"(a,f13.7)") " The multicenter bond order from beta density matrix: ",bndordbeta
            totbndorder=bndordalpha+bndordbeta
            write(*,"(a,f13.7)") " The sum of multicenter bond order from alpha and beta parts:    ",totbndorder
            write(*,"(a,f13.7)") " Above result in normalized form:",totbndorder/abs(totbndorder) * (abs(totbndorder)**(1D0/nbndcen))
            write(*,"(a,f13.7)") " The multicenter bond order from mixed alpha&beta density matrix:",bndordmix
            write(*,"(a,f13.7)") " Above result in normalized form:",bndordmix/abs(bndordmix) * (abs(bndordmix)**(1D0/nbndcen))
        end if
        write(*,*)
        
    else if (c1000tmp(1:1)=='-') then !Automatically search
!         do iatm=1,ncenter
!             do jatm=1,ncenter
!             if (jatm==iatm) cycle
!                 do katm=1,ncenter
!                 if (katm==iatm.or.katm==jatm) cycle
!                     do latm=1,ncenter
!                     if (latm==iatm.or.latm==jatm.or.latm==katm) cycle
!                         do matm=1,ncenter
!                         if (matm==iatm.or.matm==jatm.or.matm==katm.or.matm==latm) cycle
!                             do natm=1,ncenter
!                             if (natm==iatm.or.natm==jatm.or.natm==katm.or.natm==latm.or.natm==matm) cycle
        read(c1000tmp,*) nbndcen
        nbndcen=abs(nbndcen)
        PSmat=PSmattot
        !Search all combinations. Owing to simplicity and efficiency consideration, for open-shell system, compulsory to use mixed alpha&beta density matrix
        if (wfntype==1.or.wfntype==2.or.wfntype==4) write(*,"(a)") "Note: The bond order considered here comes from mixed alpha&beta density matrix"
        write(*,*)
        write(*,*) "The bond order larger than what value should be outputted?"
        write(*,*) "Input a value, e.g. 0.03"
        read(*,*) thres
        
        nfound=0
        maxbndord=0D0
        if (nbndcen/=3) write(*,*) "Note: The search may be not exhaustive. Please wait..."
        if (nbndcen==3) then
            !All combinations
            do iatm=1,ncenter
                do jatm=iatm+1,ncenter
                    do katm=jatm+1,ncenter
                        cenind(1)=iatm
                        cenind(2)=jatm
                        cenind(3)=katm
                        !clockwise and anticlockwise
                        do iseq=1,2
                            if (iseq==2) call invarr(cenind,1,nbndcen)
                            call calcmultibndord(nbndcen,cenind,PSmat,nbasis,accum)
                            if (accum>=thres) then
                                write(*,"(3i8,'    Result:'f16.10)") cenind(1:nbndcen),accum
                                nfound=nfound+1
                            end if
                            if (accum>maxbndord) then
                                maxbndord=accum
                                maxcenind=cenind
                            end if
                        end do
                    end do
                end do
            end do
        else if (nbndcen==4) then
        nthreads=getNThreads()
!$OMP PARALLEL DO private(iatm,jatm,katm,latm,cenind,iseq,accum) shared(nfound) schedule(dynamic) NUM_THREADS(nthreads)
            do iatm=1,ncenter
                do jatm=iatm+1,ncenter
                    do katm=jatm+1,ncenter
                        do latm=katm+1,ncenter
                            cenind(1)=iatm
                            cenind(2)=jatm
                            cenind(3)=katm
                            cenind(4)=latm
                            do iseq=1,2
                                if (iseq==2) call invarr(cenind,1,nbndcen)
                                call calcmultibndord(nbndcen,cenind,PSmat,nbasis,accum)
                                if (accum>=thres) then
                                    write(*,"(4i8,'    Result:'f16.10)") cenind(1:nbndcen),accum
                                    nfound=nfound+1
                                end if
                                if (accum>maxbndord) then
                                    maxbndord=accum
                                    maxcenind=cenind
                                end if
                            end do
                        end do
                    end do
                end do
            end do
!$OMP end parallel do
        else if (nbndcen==5) then
             do iatm=1,ncenter
        nthreads=getNThreads()
!$OMP PARALLEL DO private(jatm,katm,latm,matm,cenind,iseq,accum) shared(nfound) schedule(dynamic) NUM_THREADS(nthreads)
                do jatm=iatm+1,ncenter
                    do katm=jatm+1,ncenter
                        do latm=katm+1,ncenter
                            do matm=latm+1,ncenter
                                cenind(1)=iatm
                                cenind(2)=jatm
                                cenind(3)=katm
                                cenind(4)=latm
                                cenind(5)=matm
                                do iseq=1,2
                                    if (iseq==2) call invarr(cenind,1,nbndcen)
                                    call calcmultibndord(nbndcen,cenind,PSmat,nbasis,accum)
                                    if (accum>=thres) then
                                        write(*,"(5i8,'    Result:'f16.10)") cenind(1:nbndcen),accum
                                        nfound=nfound+1
                                    end if
                                    if (accum>maxbndord) then
                                        maxbndord=accum
                                        maxcenind=cenind
                                    end if
                                end do
                            end do
                        end do
                    end do
                end do
!$OMP end parallel do
            end do
        else if (nbndcen==6) then
            do iatm=1,ncenter
        nthreads=getNThreads()
!$OMP PARALLEL DO private(jatm,katm,latm,matm,cenind,iseq,accum) shared(nfound) schedule(dynamic) NUM_THREADS(nthreads)
                do jatm=iatm+1,ncenter
                    do katm=jatm+1,ncenter
                        do latm=katm+1,ncenter
                            do matm=latm+1,ncenter
                                do natm=matm+1,ncenter
                                    cenind(1)=iatm
                                    cenind(2)=jatm
                                    cenind(3)=katm
                                    cenind(4)=latm
                                    cenind(5)=matm
                                    cenind(6)=natm
                                    do iseq=1,2
                                        if (iseq==2) call invarr(cenind,1,nbndcen)
                                        call calcmultibndord(nbndcen,cenind,PSmat,nbasis,accum)
                                        if (accum>=thres) then
                                            write(*,"(6i8,'    Result:'f16.10)") cenind(1:nbndcen),accum
                                            nfound=nfound+1
                                        end if
                                        if (accum>maxbndord) then
                                            maxbndord=accum
                                            maxcenind=cenind
                                        end if
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
!$OMP end parallel do
                write(*,"(a,i5,a,i5)") " Finished:",iatm,"/",ncenter
            end do
        end if
        if (nfound==0) write(*,*) "No multi-center bonds above criteria were found"
        write(*,*)
        write(*,*) "The maximum bond order is"
        if (nbndcen==3) write(*,"(3i8,'    Result:'f16.10)") maxcenind(1:nbndcen),maxbndord
        if (nbndcen==4) write(*,"(4i8,'    Result:'f16.10)") maxcenind(1:nbndcen),maxbndord
        if (nbndcen==5) write(*,"(5i8,'    Result:'f16.10)") maxcenind(1:nbndcen),maxbndord
        if (nbndcen==6) write(*,"(6i8,'    Result:'f16.10)") maxcenind(1:nbndcen),maxbndord
        write(*,*)
    end if
end do
end subroutine

!!!Directly calculate multi-center bond order without complex things
! Return result
subroutine calcmultibndord(nbndcen,cenind,PSmat,matdim,result)
use defvar
implicit real*8(a-h,o-z)
real*8 PSmat(matdim,matdim)
integer nbndcen,cenind(12)
result=0D0
!Two centers
if (nbndcen==2) then
    do ib=basstart(cenind(2)),basend(cenind(2))
        do ia=basstart(cenind(1)),basend(cenind(1))
    result=result+PSmat(ia,ib)*PSmat(ib,ia)
        end do
    end do
!Three centers
else if (nbndcen==3) then
    do ic=basstart(cenind(3)),basend(cenind(3))
        do ib=basstart(cenind(2)),basend(cenind(2))
            do ia=basstart(cenind(1)),basend(cenind(1))
    result=result+PSmat(ia,ib)*PSmat(ib,ic)*PSmat(ic,ia)
            end do
        end do
    end do
!Four centers
else if (nbndcen==4) then
    do id=basstart(cenind(4)),basend(cenind(4))
        do ic=basstart(cenind(3)),basend(cenind(3))
            do ib=basstart(cenind(2)),basend(cenind(2))
                do ia=basstart(cenind(1)),basend(cenind(1))
    result=result+PSmat(ia,ib)*PSmat(ib,ic)*PSmat(ic,id)*PSmat(id,ia)
                end do
            end do
        end do
    end do
!Five centers
else if (nbndcen==5) then
    do ie=basstart(cenind(5)),basend(cenind(5))
        do id=basstart(cenind(4)),basend(cenind(4))
            do ic=basstart(cenind(3)),basend(cenind(3))
                do ib=basstart(cenind(2)),basend(cenind(2))
                    do ia=basstart(cenind(1)),basend(cenind(1))
    result=result+PSmat(ia,ib)*PSmat(ib,ic)*PSmat(ic,id)*PSmat(id,ie)*PSmat(ie,ia)
                    end do
                end do
            end do
        end do
    end do
!Six centers, you can easily extend this analysis to more than six centers, by simply adding cycles
else if (nbndcen==6) then
    do i_f=basstart(cenind(6)),basend(cenind(6))
        do ie=basstart(cenind(5)),basend(cenind(5))
            do id=basstart(cenind(4)),basend(cenind(4))
                do ic=basstart(cenind(3)),basend(cenind(3))
                    do ib=basstart(cenind(2)),basend(cenind(2))
                        do ia=basstart(cenind(1)),basend(cenind(1))
    result=result+PSmat(ia,ib)*PSmat(ib,ic)*PSmat(ic,id)*PSmat(id,ie)*PSmat(ie,i_f)*PSmat(i_f,ia)
                        end do
                    end do
                end do
            end do
        end do
    end do
else if (nbndcen==7) then
    do i_g=basstart(cenind(7)),basend(cenind(7))
        do i_f=basstart(cenind(6)),basend(cenind(6))
            do ie=basstart(cenind(5)),basend(cenind(5))
                do id=basstart(cenind(4)),basend(cenind(4))
                    do ic=basstart(cenind(3)),basend(cenind(3))
                        do ib=basstart(cenind(2)),basend(cenind(2))
                            do ia=basstart(cenind(1)),basend(cenind(1))
    result=result+PSmat(ia,ib)*PSmat(ib,ic)*PSmat(ic,id)*PSmat(id,ie)*PSmat(ie,i_f)*PSmat(i_f,i_g)*PSmat(i_g,ia)
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
else if (nbndcen==8) then
    do i_h=basstart(cenind(8)),basend(cenind(8))
        do i_g=basstart(cenind(7)),basend(cenind(7))
            do i_f=basstart(cenind(6)),basend(cenind(6))
                do ie=basstart(cenind(5)),basend(cenind(5))
                    do id=basstart(cenind(4)),basend(cenind(4))
                        do ic=basstart(cenind(3)),basend(cenind(3))
                            do ib=basstart(cenind(2)),basend(cenind(2))
                                do ia=basstart(cenind(1)),basend(cenind(1))
    result=result+PSmat(ia,ib)*PSmat(ib,ic)*PSmat(ic,id)*PSmat(id,ie)*PSmat(ie,i_f)*PSmat(i_f,i_g)*PSmat(i_g,i_h)*PSmat(i_h,ia)
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
else if (nbndcen==9) then
    do i_i=basstart(cenind(9)),basend(cenind(9))
        do i_h=basstart(cenind(8)),basend(cenind(8))
            do i_g=basstart(cenind(7)),basend(cenind(7))
                do i_f=basstart(cenind(6)),basend(cenind(6))
                    do ie=basstart(cenind(5)),basend(cenind(5))
                        do id=basstart(cenind(4)),basend(cenind(4))
                            do ic=basstart(cenind(3)),basend(cenind(3))
                                do ib=basstart(cenind(2)),basend(cenind(2))
                                    do ia=basstart(cenind(1)),basend(cenind(1))
    result=result+PSmat(ia,ib)*PSmat(ib,ic)*PSmat(ic,id)*PSmat(id,ie)*PSmat(ie,i_f)*PSmat(i_f,i_g)*PSmat(i_g,i_h)*PSmat(i_h,i_i)*PSmat(i_i,ia)
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
else if (nbndcen==10) then
    itmp=0
    ntot=basend(cenind(10))-basstart(cenind(10))+1
    do i_j=basstart(cenind(10)),basend(cenind(10))
        itmp=itmp+1
        write(*,"(' Finished:',i5,'  /',i5)") itmp,ntot
        do i_i=basstart(cenind(9)),basend(cenind(9))
            do i_h=basstart(cenind(8)),basend(cenind(8))
                do i_g=basstart(cenind(7)),basend(cenind(7))
                    do i_f=basstart(cenind(6)),basend(cenind(6))
                        do ie=basstart(cenind(5)),basend(cenind(5))
                            do id=basstart(cenind(4)),basend(cenind(4))
                                do ic=basstart(cenind(3)),basend(cenind(3))
                                    do ib=basstart(cenind(2)),basend(cenind(2))
                                        do ia=basstart(cenind(1)),basend(cenind(1))
    result=result+PSmat(ia,ib)*PSmat(ib,ic)*PSmat(ic,id)*PSmat(id,ie)*PSmat(ie,i_f)*PSmat(i_f,i_g)*PSmat(i_g,i_h)*PSmat(i_h,i_i)*PSmat(i_i,i_j)*PSmat(i_j,ia)
                                        end do
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
else if (nbndcen==11) then
    itmp=0
    ntot=( basend(cenind(11))-basstart(cenind(11))+1 ) * ( basend(cenind(10))-basstart(cenind(10))+1 )
    do i_k=basstart(cenind(11)),basend(cenind(11))
        do i_j=basstart(cenind(10)),basend(cenind(10))
            itmp=itmp+1
            write(*,"(' Finished:',i5,'  /',i5)") itmp,ntot
            do i_i=basstart(cenind(9)),basend(cenind(9))
                do i_h=basstart(cenind(8)),basend(cenind(8))
                    do i_g=basstart(cenind(7)),basend(cenind(7))
                        do i_f=basstart(cenind(6)),basend(cenind(6))
                            do ie=basstart(cenind(5)),basend(cenind(5))
                                do id=basstart(cenind(4)),basend(cenind(4))
                                    do ic=basstart(cenind(3)),basend(cenind(3))
                                        do ib=basstart(cenind(2)),basend(cenind(2))
                                            do ia=basstart(cenind(1)),basend(cenind(1))
    result=result+PSmat(ia,ib)*PSmat(ib,ic)*PSmat(ic,id)*PSmat(id,ie)*PSmat(ie,i_f)*PSmat(i_f,i_g)*PSmat(i_g,i_h)*PSmat(i_h,i_i)*PSmat(i_i,i_j)*PSmat(i_j,i_k)*PSmat(i_k,ia)
                                            end do
                                        end do
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
else if (nbndcen==12) then
    itmp=0
    ntot=( basend(cenind(12))-basstart(cenind(12))+1 ) * ( basend(cenind(11))-basstart(cenind(11))+1 ) * ( basend(cenind(10))-basstart(cenind(10))+1 )
    do i_l=basstart(cenind(12)),basend(cenind(12))
        do i_k=basstart(cenind(11)),basend(cenind(11))
            do i_j=basstart(cenind(10)),basend(cenind(10))
                itmp=itmp+1
                write(*,"(' Finished:',i5,'  /',i5)") itmp,ntot
                do i_i=basstart(cenind(9)),basend(cenind(9))
                    do i_h=basstart(cenind(8)),basend(cenind(8))
                        do i_g=basstart(cenind(7)),basend(cenind(7))
                            do i_f=basstart(cenind(6)),basend(cenind(6))
                                do ie=basstart(cenind(5)),basend(cenind(5))
                                    do id=basstart(cenind(4)),basend(cenind(4))
                                        do ic=basstart(cenind(3)),basend(cenind(3))
                                            do ib=basstart(cenind(2)),basend(cenind(2))
                                                do ia=basstart(cenind(1)),basend(cenind(1))
    result=result+PSmat(ia,ib)*PSmat(ib,ic)*PSmat(ic,id)*PSmat(id,ie)*PSmat(ie,i_f)*PSmat(i_f,i_g)*PSmat(i_g,i_h)*PSmat(i_h,i_i)*PSmat(i_i,i_j)*PSmat(i_j,i_k)*PSmat(i_k,i_l)*PSmat(i_l,ia)
                                                end do
                                            end do
                                        end do
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
end if
end subroutine





!!------ Multicenter bond order analysis in NAO basis
!Load NBO output file with DMNAO keyword
subroutine multicenterNAO
use defvar
use util
implicit real*8 (a-h,o-z)
integer cenind(12)
integer,allocatable :: NAOinit(:),NAOend(:)
character :: c80tmp*80,c80tmp2*80,c1000tmp*1000 !type2char(0:2)=(/"Cor","Val","Ryd"/)
character,allocatable :: centername(:)*2
real*8,allocatable :: DMNAO(:,:),DMNAOb(:,:),DMNAOtot(:,:),tmpmat(:,:)

open(10,file=filename,status="old")

call loclabel(10,"NAO density matrix:",ifound,1)
if (ifound==0) then
    write(*,"(a)") "Error: Cannot found density matrix in NAO basis in the input file! Please read manual carefully"
    write(*,*)
else !Acquire number of NAOs and centers
    call loclabel(10,"NATURAL POPULATIONS",ifound,1)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    ilastspc=0
    do while(.true.) !Find how many centers and how many NAOs. We need to carefully check where is ending
        read(10,"(a)") c80tmp
        if (c80tmp==' '.or.index(c80tmp,"low occupancy")/=0.or.index(c80tmp,"Population inversion found")/=0.or.index(c80tmp,"effective core potential")/=0) then
            if (ilastspc==1) then
                ncenter=iatm
                numNAO=inao
                exit
            end if
            ilastspc=1 !last line is space
        else
            read(c80tmp,*) inao,c80tmp2,iatm
            ilastspc=0
        end if
    end do
    write(*,"(' The number of atoms:',i10)") ncenter
    write(*,"(' The number of NAOs: ',i10)") numNAO
    write(*,*)
    allocate(NAOinit(ncenter),NAOend(ncenter),centername(ncenter))
    !Get relationship between center and NAO indices
    call loclabel(10,"NATURAL POPULATIONS",ifound,1)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    ilastspc=1
    do while(.true.)
        read(10,"(a)") c80tmp
        if (c80tmp/=' ') then
            read(c80tmp,*) inao,c80tmp2,iatm
            if (ilastspc==1) NAOinit(iatm)=inao
            ilastspc=0
        else
            NAOend(iatm)=inao
            centername(iatm)=c80tmp2
            if (iatm==ncenter) exit
            ilastspc=1
        end if
    end do
end if

!Determine if this file is close or open shell calculation
call loclabel(10,"*******         Alpha spin orbitals         *******",ispinmode,1)
if (ispinmode==0) write(*,*) "This is close-shell calculation"
if (ispinmode==1) write(*,*) "This is open-shell calculation"

!! Read density matrix in NAO basis
!Once some linearly dependent functions are eliminated, the actual matrix dimension may be larger than numNAO,
!in this case we load actual size of DMNAO as ndim*ndim, and copy the first numNAO*numNAO block to DMNAO array
!In particular
!NBO3: The ndim >= numNAO
!NBO6: The ndim always equals to numNAO
call loclabel(10,"NAO density matrix:",ifound,1)
read(10,*)
read(10,*)
read(10,*)
read(10,*)
ndim=0
do while(.true.)
    read(10,"(a)") c80tmp
    if (c80tmp==" ") exit
    ndim=ndim+1
end do
allocate(tmpmat(ndim,ndim)) !Temporary matrix to load the entire DMNAO matrix
allocate(DMNAO(numNAO,numNAO)) !Actually used DMNAO
if (ispinmode==0) write(*,*) "Loading NAO density matrix"
if (ispinmode==1) write(*,*) "Loading NAO density matrix of Alpha electrons"
call loclabel(10,"NAO density matrix:",ifound,1)
call readmatgau(10,tmpmat,0,"f8.4 ",16,8,3)
DMNAO=tmpmat(1:numNAO,1:numNAO)
!Read beta density matrix in NAO basis
if (ispinmode==1) then
    write(*,*) "Loading NAO density matrix of Beta electrons"
    allocate(DMNAOb(numNAO,numNAO),DMNAOtot(numNAO,numNAO))
    call loclabel(10,"NAO density matrix:",ifound,0)
    call readmatgau(10,tmpmat,0,"f8.4 ",16,8,3)
    DMNAOb=tmpmat(1:numNAO,1:numNAO)
    DMNAOtot=DMNAO+DMNAOb
end if
close(10)
write(*,*)

allocate(basstart(ncenter),basend(ncenter))
basstart=NAOinit
basend=NAOend

! call showmatgau(DMNAO,"Density matrix in NAO basis ",0,"f14.8",6)
do while(.true.)
    write(*,*) "Input atom indices, e.g. 3,4,7,8,10    (Up to 12 atoms)"
    write(*,*) "Input 0 can exit"
    read(*,"(a)") c1000tmp
    if (c1000tmp(1:1)=='0') then
        deallocate(basstart,basend)
        return
    else
        call str2arr(c1000tmp,nbndcen,cenind)
        if (nbndcen>=7) write(*,*) "Please wait..."
        if (ispinmode==0) then
            call calcmultibndord(nbndcen,cenind,DMNAO,numNAO,bndord)
            write(*,"(a,f16.10)") " The bond order is",bndord
        else
            call calcmultibndord(nbndcen,cenind,DMNAO,numNAO,bndordalpha)
            write(*,"(a,f16.10)") " The bond order from alpha density matrix:",bndordalpha
            call calcmultibndord(nbndcen,cenind,DMNAOb,numNAO,bndordbeta)
            write(*,"(a,f16.10)") " The bond order from beta density matrix: ",bndordbeta
            call calcmultibndord(nbndcen,cenind,DMNAOtot,numNAO,bndordmix)
            write(*,"(a,f16.10)") " The bond order from mixed alpha&beta density matrix: ",bndordmix
            write(*,"(a,f16.10)") " The sum of bond order from alpha&beta density matrix:",bndordalpha+bndordbeta
        end if
        write(*,*)
    end if
end do
end subroutine





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Calculate Mulliken bond order
subroutine mullikenbndord
use defvar
use util
implicit real*8 (a-h,o-z)
real*8 :: PSmata(nbasis,nbasis),PSmatb(nbasis,nbasis)
real*8 :: bndmattot(ncenter,ncenter),bndmatb(ncenter,ncenter),bndmata(ncenter,ncenter)
character selectyn
bndmattot=0D0
if (wfntype==0.or.wfntype==3) then
    PSmata=Sbas*Ptot !Condensed to basis function matrix
    do i=1,ncenter !Contract PSmata to Condensed to "Condensed to atoms" 
        do j=i+1,ncenter
            bndmattot(i,j)=sum(PSmata(basstart(i):basend(i),basstart(j):basend(j)))
        end do
    end do
    bndmattot=2*(bndmattot+transpose(bndmattot))
    forall (i=1:ncenter) bndmattot(i,i)=sum(bndmattot(i,:))
else if (wfntype==1.or.wfntype==2.or.wfntype==4) then
    bndmata=0D0
    bndmatb=0D0
    PSmata=Palpha*Sbas
    PSmatb=Pbeta*Sbas
    do i=1,ncenter
        do j=i+1,ncenter
            bndmata(i,j)=sum(PSmata(basstart(i):basend(i),basstart(j):basend(j)))
            bndmatb(i,j)=sum(PSmatb(basstart(i):basend(i),basstart(j):basend(j)))
        end do
    end do
    bndmata=2*(bndmata+transpose(bndmata))
    bndmatb=2*(bndmatb+transpose(bndmatb))
    forall (i=1:ncenter) bndmata(i,i)=sum(bndmata(i,:))
    forall (i=1:ncenter) bndmatb(i,i)=sum(bndmatb(i,:))
    bndmattot=bndmata+bndmatb
end if

write(*,"(' The absolute value of bond order >=',f10.6)") bndordthres
itmp=0
do i=1,ncenter
    do j=i+1,ncenter
        if (wfntype==0.or.wfntype==3) then
            if (abs(bndmattot(i,j))>=bndordthres) then
                itmp=itmp+1
                write(*,"(' #',i5,':',5x,i5,a,i5,a,f14.8)") itmp,i,'('//a(i)%name//')',j,'('//a(j)%name//')',bndmattot(i,j)
            end if
        else if (wfntype==1.or.wfntype==2.or.wfntype==4) then
            if (abs(bndmata(i,j)+bndmatb(i,j))>=bndordthres) then
                itmp=itmp+1
                write(*,"(' #',i5,':',i5,a,i5,a,' Alpha: ',f10.6,' Beta:',f10.6,' Total:',f10.6)") &
                itmp,i,'('//a(i)%name//')',j,'('//a(j)%name//')',bndmata(i,j),bndmatb(i,j),bndmattot(i,j)
            end if
        end if
    end do
end do
write(*,*)

!Between fragment
if (allocated(frag1)) then
    bndordfraga=0
    bndordfragb=0
    bndordfragtot=0
    do i=1,size(frag1)
        do j=1,size(frag2)
            bndordfraga=bndordfraga+bndmata(frag1(i),frag2(j))
            bndordfragb=bndordfragb+bndmatb(frag1(i),frag2(j))
            bndordfragtot=bndordfragtot+bndmattot(frag1(i),frag2(j))
        end do
    end do
    if (wfntype==1.or.wfntype==2.or.wfntype==4) then
        write(*,"(' The Mulliken bond order between fragment 1 and 2:')")
        write(*,"(' Alpha:',f12.6,' Beta:',f12.6,' Total:',f12.6)") bndordfraga,bndordfragb,bndordfragtot
    else if (wfntype==0.or.wfntype==3) then
        write(*,"(' The Mulliken bond order between fragment 1 and 2:',f12.6)") bndordfragtot
    end if
    write(*,*)
end if

write(*,*) "If output bond order matrix to bndmat.txt in current folder? (y/n)"
read(*,*) selectyn
if (selectyn=='y'.or.selectyn=='Y') then
    open(10,file="bndmat.txt",status="replace")
    write(10,*) "Note:The diagonal elements are the sum of corresponding row elements"
    if (wfntype==0.or.wfntype==3) then
        call showmatgau(bndmattot,"Mulliken bond order matrix",0,"f14.8",10)
    else if (wfntype==1.or.wfntype==2.or.wfntype==4) then
        call showmatgau(bndmata,"Mulliken bond order matrix for alpha electrons",0,"f14.8",10)
        call showmatgau(bndmatb,"Mulliken bond order matrix for beta electrons",0,"f14.8",10)
        call showmatgau(bndmattot,"Mulliken bond order matrix all electrons",0,"f14.8",10)
    end if
    close(10)
    write(*,*) "Result have been outputted to bndmat.txt in current folder"
    write(*,*)
end if
end subroutine



!!--------- Decompose Mulliken bond order to MO contribution
subroutine decompMBO
use defvar
use util
implicit real*8 (a-h,o-z)
real*8,pointer :: ptmat(:,:)

do while(.true.)
    write(*,*) "Input index of two atom (e.g. 3,5)"
    write(*,*) "Note: Input 0,0 can return to upper level menu"
    read(*,*) ind1,ind2

    if (ind1==0.and.ind2==0) exit
    bndorda=0D0
    bndordb=0D0
    do itime=1,2
        if (itime==1) ptmat=>cobasa
        if (itime==2) ptmat=>cobasb
        if (itime==1.and.(wfntype==1.or.wfntype==4)) write(*,*) "Alpha orbitals:"
        if (itime==2.and.(wfntype==1.or.wfntype==4)) write(*,*) "Beta orbitals:"
        do imo=1,nbasis
            if (itime==1) irealmo=imo
            if (itime==2) irealmo=imo+nbasis
            if (MOocc(irealmo)==0D0) cycle
            accum=0D0
            do i=basstart(ind1),basend(ind1)
                do j=basstart(ind2),basend(ind2)
                    accum=accum+MOocc(irealmo)*ptmat(i,imo)*ptmat(j,imo)*Sbas(i,j)
                end do
            end do
            if (itime==1) bndorda=bndorda+accum*2
            if (itime==2) bndordb=bndordb+accum*2
            write(*,"(' Orbital',i6,' Occ:',f10.6,' Energy:',f12.6,' contributes',f14.8)") imo,MOocc(irealmo),MOene(irealmo),accum*2
        end do
        if (wfntype==0.or.wfntype==2.or.wfntype==3) then
            write(*,"(' Total Mulliken bond order:',f14.8,/)") bndorda
            exit
        else if (wfntype==1.or.wfntype==4) then
            if (itime==1) write(*,"(' Mulliken bond order of all alpha electrons:',f14.8,/)") bndorda
            if (itime==2) write(*,"(' Mulliken bond order of all beta electrons:',f14.8,/)") bndordb
            if (itime==2) write(*,"(' Total Mulliken bond order:',f14.8,/)") bndorda+bndordb
        end if
    end do
end do
end subroutine



!!----- Orbital occupancy-perturbed Mayer bond order (Decompose Mayer bond-order between two atoms to orbital contributions)
!--- J. Chem. Theory Comput. 2012, 8, 908, 914
!For simplicity, this routine only calculate Mayer bond for alpha and beta and then sum them up, don't concern mixed alpha+beta cases
subroutine OrbPertMayer
use defvar
implicit real*8 (a-h,o-z)
character orbtypechar*2
real*8 :: bndmata(ncenter,ncenter),bndmatb(ncenter,ncenter),bndmattot(ncenter,ncenter)
real*8,allocatable :: PSmattot(:,:),Ptottmp(:,:)
real*8,allocatable :: PSmata(:,:),PSmatb(:,:),Palphatmp(:,:),Pbetatmp(:,:)
bndmata=0D0
bndmatb=0D0
do while(.true.)
write(*,*) "Input indices of two atoms, e.g. 3,5"
    read(*,*) iatm,jatm
    if (iatm>=1.and.iatm<=ncenter.and.jatm>=1.and.jatm<=ncenter.and.iatm/=jatm) exit
    write(*,*) "Error: Invalid input, please input again"
end do
if (wfntype==0.or.wfntype==3) then !Close shell
    allocate(PSmattot(nbasis,nbasis),Ptottmp(nbasis,nbasis))
    sumupvar=0D0
    do imo=0,nmo !Cycle all MOs
        Ptottmp=Ptot !Don't use Ptot to make troubles, because Ptot is a global array
        if (imo/=0) then !Calculate perturbed density. At the first time (imo=1), we don't pertube density matrix to yield original Mayer bond order
            if (MOocc(imo)<=1D-10) cycle
            do ibas=1,nbasis
                do jbas=1,nbasis
                    Ptottmp(ibas,jbas)=Ptottmp(ibas,jbas)-MOocc(imo)*CObasa(ibas,imo)*CObasa(jbas,imo)
                end do
            end do
        end if
        PSmattot=matmul(Ptottmp,Sbas) !Calculate Mayer bond order based on Ptottmp
        bndordtot=0D0
        do ii=basstart(iatm),basend(iatm)
            do jj=basstart(jatm),basend(jatm)
                bndordtot=bndordtot+PSmattot(ii,jj)*PSmattot(jj,ii)
            end do
        end do
        if (imo==0) then
            beforepert=bndordtot
            write(*,"(' Mayer bond order before orbital occupancy-perturbation:',f12.6)") beforepert
            write(*,*)
            write(*,"(' Mayer bond order after orbital occupancy-perturbation:')")
            write(*,*) "Orbital     Occ      Energy    Bond order   Variance"
        else
            bndordvar=bndordtot-beforepert
            write(*,"(i6,f12.5,f11.5,2f12.6)") imo,MOocc(imo),MOene(imo),bndordtot,bndordvar
            sumupvar=sumupvar+bndordvar
        end if
    end do
    write(*,"(' Summing up occupancy perturbation from all orbitals:',f10.5)") sumupvar
    
else if (wfntype==1.or.wfntype==2.or.wfntype==4) then !Open shell
    sumupvar=0D0
    allocate(PSmata(nbasis,nbasis),PSmatb(nbasis,nbasis),Palphatmp(nbasis,nbasis),Pbetatmp(nbasis,nbasis))
    do imo=0,nmo
        Palphatmp=Palpha
        Pbetatmp=Pbeta
        if (imo/=0) then !The first time, we calculate actual Mayer bond order
            if (MOocc(imo)<=1D-10) cycle
            if (wfntype==1.or.wfntype==4) then
                if (imo<=nbasis) then !Alpha orbitals
                    do ibas=1,nbasis
                        do jbas=1,nbasis
                            Palphatmp(ibas,jbas)=Palphatmp(ibas,jbas)-MOocc(imo)*CObasa(ibas,imo)*CObasa(jbas,imo)
                        end do
                    end do
                else !Beta orbitals, between nbasis+1 and nmo
                    do ibas=1,nbasis
                        do jbas=1,nbasis
                            Pbetatmp(ibas,jbas)=Pbetatmp(ibas,jbas)-MOocc(imo)*CObasb(ibas,imo-nbasis)*CObasb(jbas,imo-nbasis)
                        end do
                    end do
                end if
            else if (wfntype==2) then !ROHF
                if (MOtype(imo)==0) then !Doubly occupied orbitals
                    do ibas=1,nbasis
                        do jbas=1,nbasis
                            Palphatmp(ibas,jbas)=Palphatmp(ibas,jbas)-1D0*CObasa(ibas,imo)*CObasa(jbas,imo)
                            Pbetatmp(ibas,jbas)=Pbetatmp(ibas,jbas)-1D0*CObasa(ibas,imo)*CObasa(jbas,imo) !For ROHF, Cobasb==Cobasa, and hence Cobasb is not allocated
                        end do
                    end do
                else if (MOtype(imo)==1) then !Alpha orbitals
                    do ibas=1,nbasis
                        do jbas=1,nbasis
                            Palphatmp(ibas,jbas)=Palphatmp(ibas,jbas)-1D0*CObasa(ibas,imo)*CObasa(jbas,imo)
                        end do
                    end do                
                end if
            end if
        end if
        PSmata=matmul(Palphatmp,Sbas)
        PSmatb=matmul(Pbetatmp,Sbas)
        bndorda=0D0
        bndordb=0D0
        do ii=basstart(iatm),basend(iatm)
            do jj=basstart(jatm),basend(jatm)
                bndorda=bndorda+PSmata(ii,jj)*PSmata(jj,ii)
                bndordb=bndordb+PSmatb(ii,jj)*PSmatb(jj,ii)
            end do
        end do
        bndorda=bndorda*2
        bndordb=bndordb*2
        if (imo==0) then
            beforepert=bndorda+bndordb
            write(*,"(' Mayer bond order before orbital occupancy-perturbation:')") 
            write(*,"(' Alpha:',f12.6,'  Beta:',f12.6,'  Total:',f12.6)") bndorda,bndordb,bndorda+bndordb
            write(*,*)
            write(*,"(' Mayer bond order after orbital occupancy-perturbation:')")
            write(*,*) "Orbital     Occ      Energy  Type     Alpha      Beta     Total      Variance"
        else
            bndordvar=bndorda+bndordb-beforepert
            if (MOtype(imo)==0) orbtypechar="AB"
            if (MOtype(imo)==1) orbtypechar="A "
            if (MOtype(imo)==2) orbtypechar="B "        
            write(*,"(i6,f12.6,f11.5,2x,a,2x,3f10.5,3x,f10.5)") imo,MOocc(imo),MOene(imo),orbtypechar,bndorda,bndordb,bndorda+bndordb,bndordvar
            sumupvar=sumupvar+bndordvar
        end if
    end do
    write(*,"(' Summing up occupancy perturbation from all orbitals:',f10.5)") sumupvar
end if
write(*,*)
end subroutine
