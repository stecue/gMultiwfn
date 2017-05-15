subroutine popana
use defvar
use util
implicit real*8 (a-h,o-z)
character c2000tmp*2000

if (ifragcontri==1) then
    write(*,*) "Population analysis function couldn't be used combining with self-defined fragment"
else
    do while(.true.)
        imodwfnold=imodwfn
        write(*,*) "                ============== Population analysis =============="
        if (.not.allocated(frag1)) then
            write(*,*) "-1 Define fragment"
        else
            write(*,"(a,i5)") "-1 Redefine fragment, current atoms:",size(frag1)
        end if
        write(*,*) "0 Return"
        write(*,*) "1 Hirshfeld population"
        write(*,*) "2 Voronoi deformation density (VDD) population"
    !         write(*,*) "3 Integrate electron density in voronoi cell"
    !         write(*,*) "4 Adjusted method 3 by Rousseau et al."
        if (ifiletype==1.or.ifiletype==9) then
            write(*,*) "5 Mulliken atom & basis function population analysis"
            write(*,*) "6 Lowdin population analysis"
            write(*,*) "7 Modified Mulliken population defined by Ros & Schuit (SCPA)"
            write(*,*) "8 Modified Mulliken population defined by Stout & Politzer"
            write(*,*) "9 Modified Mulliken population defined by Bickelhaupt"
        end if
        write(*,*) "10 Becke atomic charge with atomic dipole moment correction"
        write(*,*) "11 Atomic dipole corrected Hirshfeld population (ADCH)"
        write(*,*) "12 CHELPG ESP fitting charge"
        write(*,*) "13 Merz-Kollmann (MK) ESP fitting charge"
        write(*,*) "14 AIM charge"
        read(*,*) ipopsel
        
        if (ipopsel==0) then
            if (allocated(frag1)) deallocate(frag1)
            return
        else if (ipopsel==-1) then
            if (allocated(frag1)) then
                write(*,*) "Atoms in current fragment:"
                write(*,"(13i6)") frag1
                write(*,"(a)") " Input 0 to keep unchanged, or redefine fragment, e.g. 1,3-6,8,10-11 means the atoms 1,3,4,5,6,8,10,11 will constitute the fragment"
            else
                write(*,"(a)") " Input atomic indices to define the fragment. e.g. 1,3-6,8,10-11 means the atoms 1,3,4,5,6,8,10,11 will constitute the fragment"
            end if
            read(*,"(a)") c2000tmp
            if (c2000tmp(1:1)/='0') then
                if (allocated(frag1)) deallocate(frag1)
                call str2arr(c2000tmp,nfragatm)
                allocate(frag1(nfragatm))
                call str2arr(c2000tmp,nfragatm,frag1)
                if (any(frag1>ncenter)) then
                    write(*,*) "Error: Some atomic indices exceeded valid range! Please define again"
                    write(*,*)
                    deallocate(frag1)
                end if
            end if
        else if (ipopsel==1) then
            write(*,*) "Citation: Theor. Chim. Acta. (Berl), 44, 129-138 (1977)"
            call spacecharge(1)
        else if (ipopsel==2) then
            write(*,*) "Citation: Organomet., 15, 2923-2931"
            write(*,*) "Citation: J. Comput. Chem., 25, 189-210 (2004)"
            call spacecharge(2)
        else if (ipopsel==3) then
            call spacecharge(3)
        else if (ipopsel==4) then
            write(*,*) "Citation: J. Mol. Struct.(Theochem), 538, 235-238 (2001)"
            call spacecharge(4)
        else if (ipopsel==5) then
            write(*,*) "0 Return"
            write(*,*) "1 Output Mulliken charges and decompose them to MO's contribution"
            write(*,*) "2 Output gross atomic population matrix and decompose it"
            write(*,*) "3 Output gross basis function population matrix and decompose it"
            read(*,*) ipopsel2
            if (ipopsel2==0) cycle
            call MPA(ipopsel2)
        else if (ipopsel==6) then
            call symmortho
            call MPA(1)
            call dealloall
            call readinfile(firstfilename,1) !Current wavefunction has been altered, recover the initial state
        else if (ipopsel==7) then
            call MMPA(1)
        else if (ipopsel==8) then
            call MMPA(2)
        else if (ipopsel==9) then
            call Bickelhaupt
        else if (ipopsel==10) then
            call spacecharge(5)
        else if (ipopsel==11) then
            write(*,*) "Citation of ADCH: Tian Lu, Feiwu Chen, J. Theor. Comput. Chem., 11, 163 (2011)"
            call spacecharge(6)
        else if (ipopsel==12) then
            call fitESP(1)
        else if (ipopsel==13) then
            call fitESP(2)
        else if (ipopsel==14) then
            write(*,"(a)") " NOTE: AIM charges cannot be calculated in present module but can be calculated in basin analysis module, &
            please check the example given in Section 4.17.1 of the manual on how to do this"
        end if        
        if (imodwfnold==1.and.(ipopsel==1.or.ipopsel==2.or.ipopsel==6.or.ipopsel==11)) then !1,2,6,11 is the method needed to reload the initial wavefunction
            write(*,"(a)") " Note: The wavefunction file has been reloaded, your previous modifications on occupation number will be ignored"
        end if
    end do
end if
end subroutine



!!---------- Modified Mulliken population analysis defined by Bickelhaupt
subroutine Bickelhaupt
use defvar
use util
implicit real*8 (a-h,o-z)
character selectyn,chgfilename*80
real*8,target :: atmeletot(ncenter),atmelea(ncenter)
real*8,pointer :: tmpmat(:,:),tmpele(:)
do itime=1,2
    if (itime==1) then
        tmpele=>atmeletot
        tmpmat=>Ptot
    else if (itime==2) then
        tmpele=>atmelea
        tmpmat=>Palpha
    end if
    tmpele=0D0
    do i=1,ncenter
        do j=basstart(i),basend(i)
            cross=0D0
            do k=1,nbasis !This method equalvalent to use diagonal element of density matrix to partition nondiagonal element of P*S matrix
                if (k/=j) then
                    if (tmpmat(j,j)+tmpmat(k,k)/=0D0) then
                        cross=cross+tmpmat(j,j)/(tmpmat(j,j)+tmpmat(k,k))*tmpmat(j,k)*Sbas(j,k)
                    else
                        cross=cross+tmpmat(j,k)*Sbas(j,k)/2D0  !Use equivalent partition, when denominator is zero
                    end if
                end if
            end do
            tmpele(i)=tmpele(i)+tmpmat(j,j)+2*cross !Plus electrons localized in basis and partitioned cross term
        end do
    end do
    if (wfntype==0.or.wfntype==3) exit
end do

if (wfntype==0.or.wfntype==3) then
    do iatm=1,ncenter
        write(*,"(' Atom',i6,'(',a2,')','  Population:',f10.5,'  Atomic charge:',f10.5)") iatm,a(iatm)%name,atmeletot(iatm),a(iatm)%charge-atmeletot(iatm)
    end do
    write(*,"(' Total net charge:',f10.5)") sum(a(:)%charge)-sum(atmeletot(:))
else if (wfntype==1.or.wfntype==2.or.wfntype==4) then
    write(*,*) "    Atom      Alpha_pop.   Beta_pop.    Spin_pop.     Atomic charge"
    totbetapop=0D0
    do iatm=1,ncenter
        betapop=atmeletot(iatm)-atmelea(iatm)
        totbetapop=totbetapop+betapop
        write(*,"(i6,'(',a2,')',3f13.5,f16.5)") iatm,a(iatm)%name,atmelea(iatm),betapop,atmelea(iatm)-betapop,a(iatm)%charge-atmeletot(iatm)
    end do
    write(*,"(' Total net charge:',f10.5,'      Total spin electrons:',f10.5)") sum(a(:)%charge)-sum(atmeletot(:)),sum(atmelea(:))-totbetapop
end if

!Show fragment charge
if (allocated(frag1)) write(*,"(/,' Fragment charge:',f12.6)") sum(a(frag1)%charge-atmeletot(frag1))

call path2filename(firstfilename,chgfilename)
write(*,*)
write(*,"(a)") " If output atom with charges to "//trim(chgfilename)//".chg in current folder? (y/n)"
read(*,*) selectyn
if (selectyn=="y".or.selectyn=="Y") then
    open(10,file=trim(chgfilename)//".chg",status="replace")
    do i=1,ncenter
        write(10,"(a4,4f12.6)") a(i)%name,a(i)%x*b2a,a(i)%y*b2a,a(i)%z*b2a,a(i)%charge-atmeletot(i)
    end do
    close(10)
    write(*,"(a)") " Result have been saved to "//trim(chgfilename)//".chg in current folder"
    write(*,"(a)") " Columns ranging from 1 to 5 are name,X,Y,Z,charge, respectively. The ith row corresponds to the ith atom. The unit is Angstrom"
end if
end subroutine



!!---------- Modified Mulliken population analysis
! isel=1 :Defined by Ros & Schuit (SCPA)"
! isel=2 :Defined by Stout & Politzer
subroutine MMPA(isel)
use defvar
use util
implicit real*8 (a-h,o-z)
integer isel
character selectyn,chgfilename*80
real*8,target :: atmelea(ncenter),atmeleb(ncenter)
real*8,pointer :: tmpmat(:,:),tmpele(:)
atmelea=0D0
atmeleb=0D0
do itime=1,2
    do imo=1,nbasis
        if (itime==1) then
            irealmo=imo
            tmpele=>atmelea
            tmpmat=>cobasa
        else if (itime==2) then
            irealmo=imo+nbasis
            tmpele=>atmeleb
            tmpmat=>cobasb
        end if
        if (MOocc(irealmo)==0D0) cycle
        if (isel==1) allbassqr=sum(tmpmat(:,imo)**2)
        do i=1,ncenter
            atmbassqr=sum(tmpmat(basstart(i):basend(i),imo)**2)
            if (isel==1) then
                tmpele(i)=tmpele(i)+MOocc(irealmo)*atmbassqr/allbassqr
            else if (isel==2) then
                cross=0D0
                do ii=basstart(i),basend(i)
                    do jj=1,nbasis
                        denomin=tmpmat(ii,imo)**2+tmpmat(jj,imo)**2
                        if (jj/=ii.and.denomin>=1D-120) cross=cross+tmpmat(ii,imo)**2/denomin*tmpmat(ii,imo)*tmpmat(jj,imo)*Sbas(ii,jj)
                    end do
                end do
                tmpele(i)=tmpele(i)+MOocc(irealmo)*(atmbassqr+2*cross)
            end if
        end do
    end do
    if (wfntype==0.or.wfntype==3) exit
end do
if (wfntype==0.or.wfntype==3) then
    do iatm=1,ncenter
        write(*,"(' Atom',i6,'(',a2,')','  Population:',f10.5,'  Atomic charge:',f10.5)") iatm,a(iatm)%name,atmelea(iatm),a(iatm)%charge-atmelea(iatm)
    end do
    write(*,"(' Total net charge:',f10.5)") sum(a(:)%charge)-sum(atmelea(:))
    if (allocated(frag1)) write(*,"(/,' Fragment charge:',f12.6)") sum(a(frag1)%charge-atmelea(frag1))
else if (wfntype==1.or.wfntype==2.or.wfntype==4) then
    write(*,*) "    Atom      Alpha_pop.   Beta_pop.    Spin_pop.     Atomic charge"
    do iatm=1,ncenter
        write(*,"(i6,'(',a2,')',3f13.5,f16.5)") iatm,a(iatm)%name,atmelea(iatm),atmeleb(iatm),atmelea(iatm)-atmeleb(iatm),a(iatm)%charge-atmelea(iatm)-atmeleb(iatm)
    end do
    write(*,"(' Total net charge:',f10.5,'      Total spin electrons:',f10.5)") sum(a(:)%charge)-sum(atmelea(:))-sum(atmeleb(:)),sum(atmelea(:))-sum(atmeleb(:))
    if (allocated(frag1)) write(*,"(/,' Fragment charge:',f12.6)") sum(a(frag1)%charge-atmelea(frag1)-atmeleb(frag1))
end if

call path2filename(firstfilename,chgfilename)
write(*,*)
write(*,"(a)") " If output atom coordinates with charges to "//trim(chgfilename)//".chg file in current folder? (y/n)"
read(*,*) selectyn
if (selectyn=="y".or.selectyn=="Y") then
    open(10,file=trim(chgfilename)//".chg",status="replace")
    do i=1,ncenter
        if (wfntype==0.or.wfntype==3) write(10,"(a4,4f12.6)") a(i)%name,a(i)%x*b2a,a(i)%y*b2a,a(i)%z*b2a,a(i)%charge-atmelea(i)
        if (wfntype==1.or.wfntype==2.or.wfntype==4) write(10,"(a4,4f12.6)") a(i)%name,a(i)%x*b2a,a(i)%y*b2a,a(i)%z*b2a,a(i)%charge-atmelea(i)--atmeleb(i)
    end do
    close(10)
    write(*,"(a)") "Result have been saved to "//trim(chgfilename)//".chg in current folder"
    write(*,"(a)") "Columns ranging from 1 to 5 are name,X,Y,Z,charge respectively, unit is Angstrom"
end if
end subroutine



!!--------- Mulliken/Lowdin population analysis & decompose to MO contribution
! isel=1 Output Mulliken/Lowdin charge and decompose it to MO contribution"
! isel=2 Output gross atomic population matrix and decompose it"
! isel=3 Output gross basis function population matrix and decompose it"
!Note: If doing Lowdin population, density matrix and overlap matrix should be transformed first before invoking this routine
subroutine MPA(isel)
use defvar
use util
implicit real*8 (a-h,o-z)
real*8 MOcenmat(nbasis,ncenter),groatmmat(ncenter+1,ncenter),atmele(ncenter),charge(ncenter)
real*8,pointer :: ptmat(:,:)
real*8,allocatable :: tmpmat(:,:),basmata(:,:),angorbpop(:,:),angorbpopa(:,:),angorbpopb(:,:)
character selectyn,corbnum*6,cOcc*12,chgfilename*80
integer isel
! integer :: atmarr(38)=(/69,70,66,63,64,60,57,58,54,51,52,47,48,44,41,42,37,38,34,31,32,27,28,24,21,22,17,18,14,11,12,8,5,6,1,71,73,74/)
! integer :: atmarrtype(38)=(/0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0/)

if (isel==1.or.isel==2) then
    allocate(tmpmat(nbasis,nbasis),basmata(nbasis,nbasis)) !basmata stores basis gross population of alpha part
    do itime=1,3 !total elec, alpha elec, beta elec
        if (itime==1) tmpmat=Ptot*Sbas
        if (itime==2) then
            tmpmat=Palpha*Sbas
            basmata=tmpmat !Backup
        end if
        if (itime==3) tmpmat=Pbeta*Sbas
        !Calculate gross atomic population matrix
        do i=1,ncenter
            do j=1,ncenter
                accum=0D0
                do ii=basstart(i),basend(i)
                    do jj=basstart(j),basend(j)
                        accum=accum+tmpmat(ii,jj)
                    end do
                end do
                groatmmat(i,j)=accum
            end do
        end do
        do i=1,ncenter !Stored atom populations to the last row
            groatmmat(ncenter+1,i)=sum(groatmmat(1:ncenter,i))
        end do
        totelec=0D0
        
        if (isel==2) then !Output gross atomic population matrix
            if (itime==1) call showmatgau(groatmmat,"Total gross atomic population matrix",0,"f14.8")
            if (itime==2) call showmatgau(groatmmat,"Alpha gross atomic population matrix",0,"f14.8")
            if (itime==3) call showmatgau(groatmmat,"Beta gross atomic population matrix",0,"f14.8")
        else if (isel==1) then !Contract gross atomic population matrix and output population in each basis function/shell/atom
            if (wfntype==0.or.wfntype==3) then !Notice that only perform once (itime=1)
                allocate(angorbpop(ncenter,0:5)) !Record the population number in each angular moment orbitals, up to H
                angorbpop=0D0
                write(*,*) "Population of basis functions:"
                write(*,"('  Basis Type    Atom    Shell   Population')")
                do ibas=1,nbasis
                    write(*,"(i6,3x,a,i5,a,i5,f13.5)") ibas,GTFtype2name(bastype(ibas)),bascen(ibas),'('//a(bascen(ibas))%name//')',basshell(ibas),sum(tmpmat(ibas,:))
                end do
                write(*,*)
                write(*,*) "Population of shells of basis functions:"
                do ish=1,nshell
                    shellpop=0D0
                    do ibas=1,nbasis
                        if (basshell(ibas)==ish) then
                            iatm=bascen(ibas) !Which atom this shell attribute to
                            shellpop=shellpop+sum(tmpmat(ibas,:))
                        end if
                    end do
                    write(*,"(' Shell',i6,' Type: ',a,'    in atom',i5,'(',a,') :',f9.5)") ish,shelltype2name(shtype(ish)),iatm,a(iatm)%name,shellpop
                    iangtmp=abs(shtype(ish))
                    angorbpop(iatm,iangtmp)=angorbpop(iatm,iangtmp)+shellpop
                end do
                write(*,*)
                write(*,*) "Population of each type of angular moment orbitals:"
                do iatm=1,ncenter
                    write(*,"(' Atom',i6,'(',a2,')',' s:',f7.4,' p:',f7.4,' d:',f7.4,' f:',f7.4,' g:',f7.4,' h:',f7.4)") iatm,a(iatm)%name,angorbpop(iatm,:)
                end do
                write(*,"(' Sum  s:',f9.4,' p:',f9.4,' d:',f9.4,' f:',f9.4,' g:',f9.4,' h:',f9.4)") &
                sum(angorbpop(:,0)),sum(angorbpop(:,1)),sum(angorbpop(:,2)),sum(angorbpop(:,3)),sum(angorbpop(:,4)),sum(angorbpop(:,5))
                write(*,*)
                write(*,*) "Population of atoms:"
                do iatm=1,ncenter
                    charge(iatm)=a(iatm)%charge-groatmmat(ncenter+1,iatm)
                    write(*,"(' Atom',i6,'(',a2,')','    Population:',f11.5,'    Net charge:',f11.5)") iatm,a(iatm)%name,groatmmat(ncenter+1,iatm),charge(iatm)
                end do
                write(*,"(' Total net charge:',f10.5)") sum(a(:)%charge)-sum(groatmmat(ncenter+1,:))
            else if ((wfntype==1.or.wfntype==2.or.wfntype==4).and.itime==3) then !For unrestrict wfn, at last "itime" cycle print result
                allocate(angorbpopa(ncenter,0:5),angorbpopb(ncenter,0:5))
                angorbpopa=0D0
                angorbpopb=0D0
                write(*,*) "Population of basis functions:"
                write(*,"('  Basis Type    Atom    Shell   Alpha_pop.  Beta_pop.   Total_pop.  Spin_pop.')")
                do ibas=1,nbasis !Note: Currently, tmpmat stores basis gross population of beta part, basmata stores alpha part
                    baspopa=sum(basmata(ibas,:))
                    baspopb=sum(tmpmat(ibas,:))
                    write(*,"(i6,3x,a,i5,a,i5,1x,4f12.5)") ibas,GTFtype2name(bastype(ibas)),bascen(ibas),&
                    '('//a(bascen(ibas))%name//')',basshell(ibas),baspopa,baspopb,baspopa+baspopb,baspopa-baspopb
                end do
                write(*,*)
                write(*,*) "Population of shells:"
                write(*,*) "Shell  Type     Atom     Alpha_pop.  Beta_pop.   Total_pop.  Spin_pop."
                do ish=1,nshell
                    shellpopa=0D0
                    shellpopb=0D0
                    do ibas=1,nbasis
                        if (basshell(ibas)==ish) then
                            iatm=bascen(ibas) !Which atom this shell attribute to
                            shellpopa=shellpopa+sum(basmata(ibas,:))
                            shellpopb=shellpopb+sum(tmpmat(ibas,:))
                        end if
                    end do
                    write(*,"(i5,5x,a,i7,'(',a,')' ,4f12.5)") ish,shelltype2name(shtype(ish)),iatm,a(iatm)%name,shellpopa,shellpopb,shellpopa+shellpopb,shellpopa-shellpopb
                    iangtmp=abs(shtype(ish))
                    angorbpopa(iatm,iangtmp)=angorbpopa(iatm,iangtmp)+shellpopa
                    angorbpopb(iatm,iangtmp)=angorbpopb(iatm,iangtmp)+shellpopb
                end do
                write(*,*)
                write(*,*) "Population of each type of angular moment orbitals:"
                write(*,*) "    Atom    Type   Alpha_pop.   Beta_pop.    Total_pop.   Spin_pop."
                do iatm=1,ncenter
                    if (angorbpopa(iatm,0)/=0D0.or.angorbpopb(iatm,0)/=0D0) write(*,"(i6,'(',a2,')    s',4f13.5)") &
                    iatm,a(iatm)%name,angorbpopa(iatm,0),angorbpopb(iatm,0),angorbpopa(iatm,0)+angorbpopb(iatm,0),angorbpopa(iatm,0)-angorbpopb(iatm,0)
                    if (angorbpopa(iatm,1)/=0D0.or.angorbpopb(iatm,1)/=0D0) write(*,"('              p',4f13.5)") &
                    angorbpopa(iatm,1),angorbpopb(iatm,1),angorbpopa(iatm,1)+angorbpopb(iatm,1),angorbpopa(iatm,1)-angorbpopb(iatm,1)
                    if (angorbpopa(iatm,2)/=0D0.or.angorbpopb(iatm,2)/=0D0) write(*,"('              d',4f13.5)") &
                    angorbpopa(iatm,2),angorbpopb(iatm,2),angorbpopa(iatm,2)+angorbpopb(iatm,2),angorbpopa(iatm,2)-angorbpopb(iatm,2)
                    if (angorbpopa(iatm,3)/=0D0.or.angorbpopb(iatm,3)/=0D0) write(*,"('              f',4f13.5)") &
                    angorbpopa(iatm,3),angorbpopb(iatm,3),angorbpopa(iatm,3)+angorbpopb(iatm,3),angorbpopa(iatm,3)-angorbpopb(iatm,3)
                    if (angorbpopa(iatm,4)/=0D0.or.angorbpopb(iatm,4)/=0D0) write(*,"('              g',4f13.5)") &
                    angorbpopa(iatm,4),angorbpopb(iatm,4),angorbpopa(iatm,4)+angorbpopb(iatm,4),angorbpopa(iatm,4)-angorbpopb(iatm,4)
                    if (angorbpopa(iatm,5)/=0D0.or.angorbpopb(iatm,5)/=0D0) write(*,"('              h',4f13.5)") &
                    angorbpopa(iatm,5),angorbpopb(iatm,5),angorbpopa(iatm,5)+angorbpopb(iatm,5),angorbpopa(iatm,5)-angorbpopb(iatm,5)
                end do
                write(*,*)
                if (sum(angorbpopa(:,0))/=0D0.or.sum(angorbpopb(:,0))/=0D0) write(*,"('     Total    s',4f13.5)") &
                sum(angorbpopa(:,0)),sum(angorbpopb(:,0)),sum(angorbpopa(:,0))+sum(angorbpopb(:,0)),sum(angorbpopa(:,0))-sum(angorbpopb(:,0))
                if (sum(angorbpopa(:,1))/=0D0.or.sum(angorbpopb(:,1))/=0D0) write(*,"('              p',4f13.5)") &
                sum(angorbpopa(:,1)),sum(angorbpopb(:,1)),sum(angorbpopa(:,1))+sum(angorbpopb(:,1)),sum(angorbpopa(:,1))-sum(angorbpopb(:,1))
                if (sum(angorbpopa(:,2))/=0D0.or.sum(angorbpopb(:,2))/=0D0) write(*,"('              d',4f13.5)") &
                sum(angorbpopa(:,2)),sum(angorbpopb(:,2)),sum(angorbpopa(:,2))+sum(angorbpopb(:,2)),sum(angorbpopa(:,2))-sum(angorbpopb(:,2))
                if (sum(angorbpopa(:,3))/=0D0.or.sum(angorbpopb(:,3))/=0D0) write(*,"('              f',4f13.5)") &
                sum(angorbpopa(:,3)),sum(angorbpopb(:,3)),sum(angorbpopa(:,3))+sum(angorbpopb(:,3)),sum(angorbpopa(:,3))-sum(angorbpopb(:,3))
                if (sum(angorbpopa(:,4))/=0D0.or.sum(angorbpopb(:,4))/=0D0) write(*,"('              g',4f13.5)") &
                sum(angorbpopa(:,4)),sum(angorbpopb(:,4)),sum(angorbpopa(:,4))+sum(angorbpopb(:,4)),sum(angorbpopa(:,4))-sum(angorbpopb(:,4))
                if (sum(angorbpopa(:,5))/=0D0.or.sum(angorbpopb(:,5))/=0D0) write(*,"('              h',4f13.5)") &
                sum(angorbpopa(:,5)),sum(angorbpopb(:,5)),sum(angorbpopa(:,5))+sum(angorbpopb(:,5)),sum(angorbpopa(:,5))-sum(angorbpopb(:,5))
                write(*,*)
                write(*,*) "Population of atoms:"
                write(*,*) "    Atom      Alpha_pop.   Beta_pop.    Spin_pop.     Atomic charge"
                totspinpop=0D0
                do iatm=1,ncenter
                    alphaele=atmele(iatm)
                    betaele=groatmmat(ncenter+1,iatm)
                    charge(iatm)=a(iatm)%charge-(alphaele+betaele)
                    write(*,"(i6,'(',a2,')',3f13.5,f16.5)") iatm,a(iatm)%name,alphaele,betaele,alphaele-betaele,charge(iatm)
                    totspinpop=totspinpop+alphaele-betaele
                end do
                write(*,"(' Total net charge:',f10.5,'      Total spin electrons:',f10.5)") sum(charge),totspinpop
            end if
            if (itime==2) atmele(:)=groatmmat(ncenter+1,:) !Store alpha occupation of each atom to a temporary array
        end if
        if (wfntype==0.or.wfntype==3) exit !RHF or ROHF, don't continue to process alpha & beta respectively
    end do
    
    !Show fragment charge
    if (allocated(frag1)) write(*,"(/,' Fragment charge:',f12.6)") sum(charge(frag1))

    !Asking user if decompose result
    write(*,*)
    if (isel==1) then
        write(*,*) "Decompose atomic charges to MO contribution? (y/n)"
    else if (isel==2) then
        write(*,*) "The last row is the sum of corresponding column elements (Atomic population)"
        write(*,*)
        write(*,"(a)") "Decompose to MO contribution and write to groatmdcp.txt in current folder? (y/n)"
    end if
    read(*,*) selectyn
    if (selectyn=='y'.or.selectyn=='Y') then
        if (isel==1) then
            write(*,*) "The (i,j) elements means the contribution to the jth atoms from the ith MO"
        else if (isel==2) then
            open(10,file="groatmdcp.txt",status="replace")
            write(10,*) "The last row is the sum of corresponding column elements"
            write(10,*)
        end if
        do itime=1,2
            MOcenmat=0D0
            if (itime==1) ptmat=>cobasa
            if (itime==2) ptmat=>cobasb
            do imo=1,nbasis
                if (itime==1) irealmo=imo
                if (itime==2) irealmo=imo+nbasis
                write(corbnum,"(i6)") imo
                write(cOcc,"(f12.8)") MOocc(irealmo)
                if (MOocc(irealmo)==0D0) cycle

                do i=1,ncenter
                    do j=1,ncenter
                        accum=0D0
                        do ii=basstart(i),basend(i)
                            do jj=basstart(j),basend(j)
                                accum=accum+MOocc(irealmo)*ptmat(ii,imo)*ptmat(jj,imo)*Sbas(ii,jj)
                            end do
                        end do
                        groatmmat(i,j)=accum
                    end do
                end do
                do i=1,ncenter
                    groatmmat(ncenter+1,i)=sum(groatmmat(1:ncenter,i))
                end do
                if (isel==1) then
                    MOcenmat(imo,:)=groatmmat(ncenter+1,:)
                else if (isel==2) then
                    if (wfntype==0.or.wfntype==2.or.wfntype==3) then
                        call showmatgau(groatmmat,"Orbital"//corbnum//"  Occ:"//cOcc,0,"f14.8",10)
                    else if (itime==1.and.(wfntype==1.or.wfntype==4)) then
                        call showmatgau(groatmmat,"Alpha Orbital"//corbnum//"  Occ:"//cOcc,0,"f14.8",10)
                    else if (itime==2.and.(wfntype==1.or.wfntype==4)) then
                        call showmatgau(groatmmat,"Beta Orbital"//corbnum//"  Occ:"//cOcc,0,"f14.8",10)
                    end if
                    write(10,*)
                end if
            end do
            if (isel==1) then
                if (wfntype==0.or.wfntype==2) call showmatgau(MOcenmat,"",0,"f14.8",6,nint(naelec))
                if (itime==1.and.wfntype==1) call showmatgau(MOcenmat,"Alpha part",0,"f14.8",6,nint(naelec))
                if (itime==2.and.wfntype==1) call showmatgau(MOcenmat,"Beta part",0,"f14.8",6,nint(nbelec))
                if (wfntype==3) call showmatgau(MOcenmat,"",0,"f14.8")
                if (itime==1.and.wfntype==4) call showmatgau(MOcenmat,"Alpha part",0,"f14.8")
                if (itime==2.and.wfntype==4) call showmatgau(MOcenmat,"Beta part",0,"f14.8")
            end if
            if (wfntype==0.or.wfntype==2.or.wfntype==3) exit !ROHF needn't to separate to alpha and beta
        end do
        if (isel==2) close(10)
        if (isel==2) write(*,*) "Done!"
    end if
    
    !If calculating atomic charges, asking user if output it
    if (isel==1) then
        call path2filename(firstfilename,chgfilename)
        write(*,*)
        write(*,"(a)") " If output atom with charges to "//trim(chgfilename)//".chg in current folder? (y/n)"
        read(*,*) selectyn
        if (selectyn=="y".or.selectyn=="Y") then
            open(10,file=trim(chgfilename)//".chg",status="replace")
            do i=1,ncenter
                if (wfntype==0.or.wfntype==3) write(10,"(a4,4f12.6)") a(i)%name,a(i)%x*b2a,a(i)%y*b2a,a(i)%z*b2a,charge(i)
                if (wfntype==1.or.wfntype==2.or.wfntype==4) write(10,"(a4,4f12.6)") a(i)%name,a(i)%x*b2a,a(i)%y*b2a,a(i)%z*b2a,charge(i)
            end do
            close(10)
            write(*,"(a)") " Result have been saved to "//trim(chgfilename)//".chg in current folder"
            write(*,"(a)") " Columns ranging from 1 to 5 are name,X,Y,Z,charge respectively, unit is Angstrom"
        end if
    end if
        
else if (isel==3) then
    call showmatgau(Ptot*Sbas,"Total gross basis function population matrix",0,"f14.8")
    if (wfntype==1.or.wfntype==2.or.wfntype==4) then
        call showmatgau(Palpha*Sbas,"Alpha gross basis function population matrix",0,"f14.8")
        call showmatgau(Pbeta*Sbas,"Beta gross basis function population matrix",0,"f14.8")
    end if
    write(*,*) "The (i,j) element means ¡Æ[imo] Occ(imo)*C(i,imo)*C(j,imo)*S(i,j)"
    write(*,*)
    write(*,"(a)") "Decompose to MO contribution and write to grobasdcp.txt in current folder? (y/n)"
    read(*,*) selectyn
    if (selectyn=='y'.or.selectyn=='Y') then
        open(10,file="grobasdcp.txt",status="replace")
        write(10,*) "Notes:"
        write(10,*) "The (i,j) element means C(i,imo)*C(j,imo)*S(i,j)"
        write(10,*) "The last row is the sum of corresponding column elements"
        write(10,*)
        allocate(tmpmat(nbasis+1,nbasis))
        do itime=1,2
            if (itime==1) ptmat=>cobasa
            if (itime==2) ptmat=>cobasb
            do imo=1,nbasis
                if (itime==1) irealmo=imo
                if (itime==2) irealmo=imo+nbasis
                write(corbnum,"(i6)") imo
                write(cOcc,"(f12.8)") MOocc(irealmo)
                if (MOocc(irealmo)==0D0) cycle

                do ii=1,nbasis
                    do jj=1,nbasis
                        tmpmat(ii,jj)=MOocc(irealmo)*ptmat(ii,imo)*ptmat(jj,imo)*Sbas(ii,jj)
                    end do
                end do
                do j=1,nbasis
                    tmpmat(nbasis+1,j)=sum(tmpmat(1:nbasis,j))
                end do
                if (wfntype==0.or.wfntype==2.or.wfntype==3) then
                    call showmatgau(tmpmat,"Orbital"//corbnum//"  Occ:"//cOcc,0,"f14.8",10)
                else if (itime==1.and.(wfntype==1.or.wfntype==4)) then
                    call showmatgau(tmpmat,"Alpha Orbital"//corbnum//"  Occ:"//cOcc,0,"f14.8",10)
                else if (itime==2.and.(wfntype==1.or.wfntype==4)) then
                    call showmatgau(tmpmat,"Beta Orbital"//corbnum//"  Occ:"//cOcc,0,"f14.8",10)
                end if
                write(10,*)
            end do
            if (wfntype==0.or.wfntype==2.or.wfntype==3) exit
        end do
        close(10)
        write(*,*) "Done!"
    end if
end if
end subroutine




!!-------------- Calculate charge based on space partition method
subroutine spacecharge(chgtype)
!1=Hirshfeld, 2=VDD, 3=Integrate electron density in voronoi cell
!4=Adjusted method 3 by Rousseau et al., 5= Becke with/without ADC, 6= ADCH
use defvar
use function
use util
implicit real*8(a-h,o-z)
! integer :: atmarr(38)=(/69,70,66,63,64,60,57,58,54,51,52,47,48,44,41,42,37,38,34,31,32,27,28,24,21,22,17,18,14,11,12,8,5,6,1,71,73,74/)
! integer :: atmarrtype(38)=(/0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0/)
integer chgtype
real*8 molrho(radpot*sphpot),promol(radpot*sphpot),tmpdens(radpot*sphpot),beckeweigrid(radpot*sphpot),selfdens(radpot*sphpot)
type(content) gridatm(radpot*sphpot),gridatmorg(radpot*sphpot)
real*8 atmdipx(ncenter),atmdipy(ncenter),atmdipz(ncenter),charge(ncenter)
real*8 :: covr_becke(0:nelesupp) !covalent radii used for Becke population
character selectyn,chgfilename*800
character radfilename*800
integer :: nbeckeiter=3

if (chgtype==5) then !Select atomic radii for Becke population
    covr_becke=covr_TianLu
    iraddefine=2
    do while(.true.)
        write(*,*) "-1 Return"
        write(*,*) "0 Start calculation of Becke charge!"
        if (iraddefine==0) write(*,*) "1 Select the definition of atomic radii, current: Custom"
        if (iraddefine==1) write(*,*) "1 Select the definition of atomic radii, current: CSD"
        if (iraddefine==2) write(*,*) "1 Select the definition of atomic radii, current: Modified CSD"
        if (iraddefine==3) write(*,*) "1 Select the definition of atomic radii, current: Pyykko"
        if (iraddefine==4) write(*,*) "1 Select the definition of atomic radii, current: Suresh"
        if (iraddefine==5) write(*,*) "1 Select the definition of atomic radii, current: Hugo"
        write(*,"(a,i2)") " 2 Set the number of iterations for defining Becke atomic space, current:",nbeckeiter
        write(*,*) "10 Read radii from external file"
        write(*,*) "11 Modify current radii by manual input"
        write(*,*) "12 Print current radii list"
        read(*,*) isel
        if (isel==-1) then
            return
        else if (isel==0) then
            exit
        else if (isel==1) then
            write(*,*) "1 Use CSD radii (Dalton Trans., 2008, 2832-2838)"
            write(*,*) "2 Use the modified version of CSD radii defined by Tian Lu (Recommended)"
            write(*,*) "3 Use Pyykko radii (Chem. Eur.-J., 15, 186-197)"
            write(*,*) "4 Use Suresh radii (J. Phys. Chem. A, 105, 5940-5944)"
            write(*,*) "5 Use Hugo radii (Chem. Phys. Lett., 480, 127-131)"
            read(*,*) iselrad
            if (iselrad==1) then
                covr_becke=covr
                iraddefine=1
            else if (iselrad==2) then
                covr_becke=covr_TianLu
                iraddefine=2
            else if (iselrad==3) then
                covr_becke=covr_pyy
                iraddefine=3
            else if (iselrad==4) then
                covr_becke=covr_Suresh
                iraddefine=4
            else if (iselrad==5) then
                covr_becke=radii_hugo
                iraddefine=5
            end if
        else if (isel==2) then
            write(*,*) "Input a number, e.g. 3"
            read(*,*) nbeckeiter
        else if (isel==10) then
            iraddefine=0
            write(*,"(a)") " About the file format:"
            write(*,"(a)") " The first line should be the number of elements you want to modify, followed by element indices and radii (in Angstrom), for example:"
            write(*,"(a)") " 4"
            write(*,"(a)") " 1 0.35"
            write(*,"(a)") " 4 1.2"
            write(*,"(a)") " 5 1.12"
            write(*,"(a)") " 14 1.63"
            write(*,*)
            write(*,*) "Input filename"
            read(*,"(a)") radfilename
            inquire(file=radfilename,exist=alive)
            if (alive.eqv..true.) then
                open(10,file=radfilename,status="old")
                read(10,*) nmodrad
                do irad=1,nmodrad
                    read(10,*) indtmp,radtmp
                    covr_becke(indtmp)=radtmp/b2a
                end do
                close(10)
                write(*,*) "Done!"
            else
                write(*,*) "Error: File cannot be found"
            end if
        else if (isel==11) then
            iraddefine=0
            write(*,*) "Input element index and radius (in Angstrom), e.g. 5,0.84"
            read(*,*) indtmp,radtmp
            covr_becke(indtmp)=radtmp/b2a
            write(*,*) "Done!"
        else if (isel==12) then
            do irad=0,nelesupp
                write(*,"(' Element:',i5,'(',a,')   Radius:',f8.3,' Angstrom')") irad,ind2name(irad),covr_becke(irad)*b2a
            end do
            write(*,*)
        end if
    end do
end if

!Generate quadrature point and weighs by combination of Gauss-Chebyshev and Lebedev grids
call gen1cintgrid(gridatmorg,iradcut)

!***** 1=Hirshfeld, 2=VDD, 6=ADCH
if (chgtype==1.or.chgtype==2.or.chgtype==6) then
    write(*,*) "This task requests atomic densities, please select how to obtain them"
    write(*,*) "1 Use build-in sphericalized atomic densities in free-states (more convenient)"
    write(*,"(a)") " 2 Provide wavefunction file of involved elements by yourself or invoke Gaussian to automatically calculate them"
    read(*,*) iatmdensmode
    if (iatmdensmode==2) call setpromol !In this routine reload first molecule at the end
    write(*,"(' Radial grids:',i5,'    Angular grids:',i5,'   Total:',i10)") radpot,sphpot,radpot*sphpot
    write(*,*) "Calculating, please wait..."
    write(*,*)
    call walltime(nwalltime1)
    do iatm=1,ncenter !Cycle each atom to calculate their charges and dipole
        call delvirorb(0) !For faster calculation, remove virtual MOs in whole system, will not affect result
        atmx=a(iatm)%x
        atmy=a(iatm)%y
        atmz=a(iatm)%z
        gridatm%value=gridatmorg%value !Weight in this grid point
        gridatm%x=gridatmorg%x+atmx !Move quadrature point to actual position in molecule
        gridatm%y=gridatmorg%y+atmy
        gridatm%z=gridatmorg%z+atmz
        !Calculate molecular density first
!$OMP parallel do shared(molrho) private(i) num_threads( nthreads  )
        do i=1,radpot*sphpot
            molrho(i)=fdens(gridatm(i)%x,gridatm(i)%y,gridatm(i)%z)
        end do
!$OMP end parallel do
        !Calc free atomic density to obtain promolecule density
        promol=0D0
        if (iatmdensmode==1) then
            do jatm=1,ncenter
!$OMP parallel do shared(tmpdens) private(ipt) num_threads( nthreads  )
                do ipt=1,radpot*sphpot
                    tmpdens(ipt)=calcatmdens(jatm,gridatm(ipt)%x,gridatm(ipt)%y,gridatm(ipt)%z,0)
                end do
!$OMP end parallel do
                promol=promol+tmpdens
                if (jatm==iatm) selfdens=tmpdens
            end do
        else if (iatmdensmode==2) then
            do jatm=1,ncenter
                call dealloall
                call readwfn(custommapname(jatm),1)
!$OMP parallel do shared(tmpdens) private(ipt) num_threads( nthreads  )
                do ipt=1,radpot*sphpot
                    tmpdens(ipt)=fdens(gridatm(ipt)%x,gridatm(ipt)%y,gridatm(ipt)%z)
                end do
!$OMP end parallel do
                promol=promol+tmpdens
                if (jatm==iatm) selfdens=tmpdens
            end do
            call dealloall
            call readinfile(firstfilename,1) !Retrieve to first loaded file(whole molecule) to calc real rho again
        end if
        !Now we have needed data in hand, calculate atomic charges and atomic dipole moments
        tmpcharge=0D0
        dipx=0D0
        dipy=0D0
        dipz=0D0
        if (chgtype==1.or.chgtype==6) then !Hirshfeld, ADCH charge
            do i=1,radpot*sphpot
                if (promol(i)/=0D0) then
                    tmpv=selfdens(i)/promol(i)*molrho(i)*gridatm(i)%value
                    tmpcharge=tmpcharge-tmpv
                    dipx=dipx-(gridatm(i)%x-atmx)*tmpv
                    dipy=dipy-(gridatm(i)%y-atmy)*tmpv
                    dipz=dipz-(gridatm(i)%z-atmz)*tmpv
                end if
            end do
            charge(iatm)=a(iatm)%charge+tmpcharge
        else if (chgtype==2) then !VDD charge
            do i=1,radpot*sphpot !Cycle each grid point of iatm, if the distance between the grid point and other atom is shorter than iatm, weight=0
                vddwei=1D0
                discen2=(gridatm(i)%x-atmx)**2+(gridatm(i)%y-atmy)**2+(gridatm(i)%z-atmz)**2 !Distance between this grid and current center atom
                do jatm=1,ncenter_org !Note: Current wfn is atomic wfn, so use _org suffix
                    if (jatm==iatm) cycle
                    disother2=(gridatm(i)%x-a_org(jatm)%x)**2+(gridatm(i)%y-a_org(jatm)%y)**2+(gridatm(i)%z-a_org(jatm)%z)**2
                    if (disother2<discen2) then
                        vddwei=0D0 !Using this weight is equivalent to using Voronoi cell
                        exit
                    end if
                end do
                tmpv=vddwei*(molrho(i)-promol(i))*gridatm(i)%value
                tmpcharge=tmpcharge-tmpv
                dipx=dipx-(gridatm(i)%x-atmx)*tmpv
                dipy=dipy-(gridatm(i)%y-atmy)*tmpv
                dipz=dipz-(gridatm(i)%z-atmz)*tmpv
            end do
            charge(iatm)=tmpcharge
        end if
        atmdipx(iatm)=dipx
        atmdipy(iatm)=dipy
        atmdipz(iatm)=dipz
        if (chgtype==1.or.chgtype==6) write(*,"(' Hirshfeld charge of atom ',i5,'(',a2,')',' is',f12.6)") iatm,a_org(iatm)%name,charge(iatm)
        if (chgtype==2) write(*,"(' VDD charge of atom ',i5,'(',a2,')',' is',f12.6)") iatm,a_org(iatm)%name,charge(iatm)
    end do
    
!***** 3=Integrate electron density in Voronoi cell, 4=Adjusted method 3 by Rousseau et al
else if (chgtype==3.or.chgtype==4) then
    write(*,"(' Radial grids:',i5,'    Angular grids:',i5,'   Total:',i10)") radpot,sphpot,radpot*sphpot
    write(*,*) "Calculating, please wait..."
    write(*,*)
    call walltime(nwalltime1)
    if (chgtype==4) then !vdW radius From J.Mol.Stru.(Theo.) 538,235-238 is not identical to original definition
        vdwr(1)=0.68D0/b2a
        !B,C,N,O,F
        vdwr(5)=1.46D0/b2a
        vdwr(6)=1.46D0/b2a
        vdwr(7)=1.39D0/b2a
        vdwr(8)=1.35D0/b2a
        vdwr(9)=1.29D0/b2a
        !P S Cl
        vdwr(15)=1.78D0/b2a
        vdwr(16)=1.74D0/b2a
        vdwr(17)=1.69D0/b2a
    end if
    do iatm=1,ncenter
        tmpcharge=0D0
        dipx=0D0
        dipy=0D0
        dipz=0D0
        atmx=a(iatm)%x
        atmy=a(iatm)%y
        atmz=a(iatm)%z
        gridatm%value=gridatmorg%value
        gridatm%x=gridatmorg%x+atmx !Move quadrature point with center of current atom
        gridatm%y=gridatmorg%y+atmy
        gridatm%z=gridatmorg%z+atmz
        do i=1,radpot*sphpot
            vorwei=1.0D0
            discen2=(gridatm(i)%x-atmx)**2+(gridatm(i)%y-atmy)**2+(gridatm(i)%z-atmz)**2 !Distance between this grid and current center atom
            do jatm=1,ncenter !Determine the boundary of cell
                if (jatm==iatm) cycle
                disother2=(gridatm(i)%x-a(jatm)%x)**2+(gridatm(i)%y-a(jatm)%y)**2+(gridatm(i)%z-a(jatm)%z)**2
                if (chgtype==3) then
                    if (disother2<discen2) then
                        vorwei=0.0D0 !Use this weights equivalent to use voronoi cell
                        exit
                    end if
                else if (chgtype==4) then !Adjusted voronoi
                    vdwra=vdwr(a(iatm)%index)
                    vdwrb=vdwr(a(jatm)%index)
                    RAB=distmat(iatm,jatm)
                    rhoval=(RAB**2+discen2-disother2)/2.0D0/RAB
                    rhoa=vdwra/(vdwra+vdwrb)*RAB
                    if (rhoval>rhoa) then
                        vorwei=0.0D0
                        exit
                    end if
                end if
            end do
            if (vorwei/=0.0D0) then
                tmpv=vorwei*fdens(gridatm(i)%x,gridatm(i)%y,gridatm(i)%z)*gridatm(i)%value
                tmpcharge=tmpcharge-tmpv
                dipx=dipx-(gridatm(i)%x-atmx)*tmpv
                dipy=dipy-(gridatm(i)%y-atmy)*tmpv
                dipz=dipz-(gridatm(i)%z-atmz)*tmpv
            end if
        end do
        charge(iatm)=tmpcharge+a(iatm)%charge
        atmdipx(iatm)=dipx
        atmdipy(iatm)=dipy
        atmdipz(iatm)=dipz
        write(*,"(' The charge of atom ',i5,'(',a2,')',' is',f12.6)") iatm,a(iatm)%name,charge(iatm)
    end do
    
!***** Becke population
else if (chgtype==5) then
    write(*,"(' Radial grids:',i5,'    Angular grids:',i5,'   Total:',i10)") radpot,sphpot,radpot*sphpot
    write(*,*) "Calculating, please wait..."
    write(*,*)
    call walltime(nwalltime1)
    do iatm=1,ncenter !Cycle each atom to calculate their charges and dipole
        gridatm%value=gridatmorg%value !Weight in this grid point
        gridatm%x=gridatmorg%x+a(iatm)%x !Move quadrature point to actual position in molecule
        gridatm%y=gridatmorg%y+a(iatm)%y
        gridatm%z=gridatmorg%z+a(iatm)%z
!$OMP parallel do shared(tmpdens) private(i) num_threads( nthreads  )
        do i=1,radpot*sphpot !Calc molecular density first
            tmpdens(i)=fdens(gridatm(i)%x,gridatm(i)%y,gridatm(i)%z)
        end do
!$OMP end parallel do
        call gen1cbeckewei(iatm,iradcut,gridatm,beckeweigrid)
        tmpcharge=0D0
        dipx=0D0
        dipy=0D0
        dipz=0D0
        do i=1+iradcut*sphpot,radpot*sphpot
            tmpv=tmpdens(i)*beckeweigrid(i)*gridatm(i)%value
            tmpcharge=tmpcharge-tmpv
            dipx=dipx-(gridatm(i)%x-a(iatm)%x)*tmpv
            dipy=dipy-(gridatm(i)%y-a(iatm)%y)*tmpv
            dipz=dipz-(gridatm(i)%z-a(iatm)%z)*tmpv
        end do
        charge(iatm)=tmpcharge+a(iatm)%charge
        atmdipx(iatm)=dipx
        atmdipy(iatm)=dipy
        atmdipz(iatm)=dipz
        write(*,"(' Becke charge of atom ',i5,'(',a2,')',' is',f12.6)") iatm,a(iatm)%name,charge(iatm)
    end do
end if

write(*,*)
write(*,*) "Atomic dipole moments (a.u.):"
do iatm=1,ncenter
    write(*,"(' Atom ',i5,'(',a2,')',' in X/Y/Z:',3f11.6,' Norm:',f11.6)") iatm,a(iatm)%name,atmdipx(iatm),atmdipy(iatm),atmdipz(iatm),dsqrt(atmdipx(iatm)**2+atmdipy(iatm)**2+atmdipz(iatm)**2)
end do
xmoldip=0.0D0
ymoldip=0.0D0
zmoldip=0.0D0
do i=1,ncenter
    xmoldip=xmoldip+a(i)%x*charge(i)
    ymoldip=ymoldip+a(i)%y*charge(i)
    zmoldip=zmoldip+a(i)%z*charge(i)
end do
write(*,*)
write(*,"(' Summing up all charges:',f15.8)") sum(charge)
if (allocated(b_EDF).and.sum(charge)<=-9.5) then
    write(*,"(a)") " Warning: One or more atomic charges may be evidently incorrect. If you have used pseudopotential, &
    set ""readEDF"" in settings.ini to 0, reboot Multiwfn and redo the calculation, then you will get correct result"
    pause
end if
totdip=dsqrt(xmoldip**2+ymoldip**2+zmoldip**2)
write(*,"(' Total dipole from atomic charges:',f12.6,' a.u.')") totdip
write(*,"(' X/Y/Z of dipole from atomic charge:',3f12.6,' a.u.')") xmoldip,ymoldip,zmoldip
totatmdip=dsqrt(sum(atmdipx)**2+sum(atmdipy)**2+sum(atmdipz)**2)
write(*,"(' Total atomic dipole:',f12.6,' a.u.')") totatmdip
write(*,"(' X/Y/Z of total atomic dipole:',3f12.6,' a.u.')") sum(atmdipx),sum(atmdipy),sum(atmdipz)
corrdipx=xmoldip+sum(atmdipx)
corrdipy=ymoldip+sum(atmdipy)
corrdipz=zmoldip+sum(atmdipz)
realdip=dsqrt(corrdipx**2+corrdipy**2+corrdipz**2)
if (chgtype/=6) then !Avoid confusing users what does "corrected" means
    write(*,"(' Total corrected dipole:',f12.6,' a.u.')") realdip
    write(*,"(' X/Y/Z of corrected dipole:',3f12.6,' a.u.')") corrdipx,corrdipy,corrdipz
end if
write(*,*)
call walltime(nwalltime2)
write(*,"(' Calculation took up',i8,' seconds wall clock time')")  nwalltime2-nwalltime1

if (chgtype==5) call doADC(atmdipx,atmdipy,atmdipz,charge,realdip,5)
if (chgtype==6) call doADC(atmdipx,atmdipy,atmdipz,charge,realdip,6)

!Show fragment charge
if (allocated(frag1)) write(*,"(/,' Fragment charge:',f12.6)") sum(charge(frag1))

call path2filename(firstfilename,chgfilename)
write(*,*)
write(*,"(a)") " If output atoms with charges to "//trim(chgfilename)//".chg in current folder? (y/n)"
read(*,*) selectyn
if (selectyn=="y".or.selectyn=="Y") then
    open(10,file=trim(chgfilename)//".chg",status="replace")
    do i=1,ncenter
        write(10,"(a4,4f12.6)") a(i)%name,a(i)%x*b2a,a(i)%y*b2a,a(i)%z*b2a,charge(i)
    end do
    close(10)
    write(*,"(a)") " Result have been saved to "//trim(chgfilename)//".chg in current folder"
    write(*,"(a)") " Columns ranging from 1 to 5 are name,X,Y,Z,charge respectively, unit is Angstrom"
end if
end subroutine


!!------ Calculate atomic dipole moment corrected charge based on existing atomic charge (charge) and atomic dipole moments (dipx/y/z)
!This routine is previously specific for ADCH, but can be extended to any other types of atomic charges
!chgtype 5= Becke with/without ADC, 6= ADCH, used to determine outputting which kind of hints
subroutine doADC(dipx,dipy,dipz,charge,realdip,chgtype)
use defvar
use util
implicit real*8 (a-h,o-z)
integer chgtype
real*8 gammamat(3,3),mat(3,3),avgr(3,1),avgrr(3,3),r(3,1),dip(3,1),tmp(1,1),eigval(3),eigvecmat(3,3)
real*8 w(ncenter),chargecorr(ncenter),charge(ncenter)
real*8 dipx(ncenter),dipy(ncenter),dipz(ncenter),realdip

write(*,*)
write(*,*) "Now calculating atomic dipole moment corrected charge..."
write(*,*)
chargecorr=charge

do i=1,ncenter
    if (ADCtransfer==1) write(*,"('Atom: 'i4,a)") i,a(i)%name !ADCtransfer==1 means output detail of charge transferation process during atomic dipole moment correction
    !Initialize variables
    totq=0.0D0
    tottmpdipx=0.0D0
    tottmpdipy=0.0D0
    tottmpdipz=0.0D0
    avgr=0.0D0
    avgrr=0.0D0
    dip(1,1)=dipx(i)
    dip(2,1)=dipy(i)
    dip(3,1)=dipz(i)

    !Calculate weight of every atom
    do j=1,ncenter
        r(1,1)=a(j)%x-a(i)%x
        r(2,1)=a(j)%y-a(i)%y
        r(3,1)=a(j)%z-a(i)%z
        r2=r(1,1)**2+r(2,1)**2+r(3,1)**2
        distij=dsqrt(r2)
        
        !Use modified Becke weight function with vdW radii criterion
        rmaxdist=vdwr(a(i)%index)+vdwr(a(j)%index)
            tr=distij/(rmaxdist/2.0D0)-1 !Transform variable so that it can in 0~rmaxdist range
            tr=1.5*(tr)-0.5*(tr)**3
            tr=1.5*(tr)-0.5*(tr)**3
            w(j)=0.5*(1-(1.5*tr-0.5*tr**3))
        if (distij>rmaxdist) w(j)=0.0D0

        avgr=avgr+w(j)*r
        avgrr=avgrr+w(j)*matmul(r,transpose(r))
    end do

    wtot=sum(w)
    avgr=avgr/wtot !Get <r>
    avgrr=avgrr/wtot !Get <rr+>
    gammamat=avgrr-matmul(avgr,transpose(avgr))
    call Diagmat(gammamat,eigvecmat,eigval,500,1D-10)

    rmaxv=maxval(eigval)
    mat=0.0D0
    tmpmin=1D-5
    mat(1,1)=1.0D0/(eigval(1)+tmpmin*(rmaxv+tmpmin))
    mat(2,2)=1.0D0/(eigval(2)+tmpmin*(rmaxv+tmpmin))
    mat(3,3)=1.0D0/(eigval(3)+tmpmin*(rmaxv+tmpmin))

    !Use transform matrix to transform r in old coordinate to r' in new coordinate, and transform P to P'
    !r=matmul(eigvecmat,r'), so r'=matmul(eigvecmat^(-1),r), because eigvecmat is unitary matrix, r'=matmul(transepose(eigvecmat),r)
    avgr=matmul(transpose(eigvecmat),avgr)
    dip=matmul(transpose(eigvecmat),dip)

    !All values need have been calculated, now calculate final result
    do j=1,ncenter
        r(1,1)=a(j)%x-a(i)%x
        r(2,1)=a(j)%y-a(i)%y
        r(3,1)=a(j)%z-a(i)%z
        r=matmul(transpose(eigvecmat),r) ! Get r(i,j) vector in new coordinate
        tmp=w(j)/wtot*matmul(matmul(transpose(r-avgr),mat) ,dip) !delta q, namely the charge which atom i gives atom j
        chargecorr(j)=chargecorr(j)+tmp(1,1) !Charge after corrected
        if (ADCtransfer==1) write(*,"(' Give atom ',i4,a4,f15.12,'  Weight',2f15.12)") j,a(j)%name,tmp(1,1),w(j)
        totq=totq+tmp(1,1)
        tottmpdipx=tottmpdipx+(a(j)%x-a(i)%x)*tmp(1,1)
        tottmpdipy=tottmpdipy+(a(j)%y-a(i)%y)*tmp(1,1)
        tottmpdipz=tottmpdipz+(a(j)%z-a(i)%z)*tmp(1,1)
    end do
    if (ADCtransfer==1) write(*,*)
end do

write(*,*) "   ======= Summary of atomic dipole moment corrected (ADC) charges ======="
do i=1,ncenter
    write(*,"(' Atom: ',i4,a,'  Corrected charge:',f12.6,'  Before:',f12.6)") i,a(i)%name,chargecorr(i),charge(i)
end do
! write(*,"(' Summing up charges before correction',f12.7)") sum(charge)
! write(*,"(' Summing up charges after correction',f12.7)") sum(chargecorr)
if (chgtype==5) write(*,"(a)") " Note: The values shown after ""Corrected charge"" are atomic dipole moment corrected Becke charges, the ones after ""Before"" are normal Becke charges"
if (chgtype==6) write(*,"(a)") " Note: The values shown after ""Corrected charge"" are ADCH charges, the ones after ""Before"" are Hirshfeld charges"
ADCdipx=sum(a%x*chargecorr)
ADCdipy=sum(a%y*chargecorr)
ADCdipz=sum(a%z*chargecorr)
ADCdip=sqrt(ADCdipx**2+ADCdipy**2+ADCdipz**2)
write(*,*)
write(*,"(' Total dipole from ADC charge(a.u.)',f11.7,'  Error:',f11.7)") ADCdip,abs(ADCdip-realdip)
write(*,"(' X/Y/Z of dipole moment from the charge(a.u.)',3f11.7)") ADCdipx,ADCdipy,ADCdipz
charge=chargecorr !Overlay charge array, then return to Hirshfeld module and output result to .chg file
end subroutine


!!------------ Calculate charges by fitting ESP, currently CHELPG grid and MK grid are used
!itype=1:CHELPG   itype=2:MK
subroutine fitESP(itype)
use util
use defvar
use function
implicit real*8 (a-h,o-z)
character*200 addcenfile,extptfile
character selectyn,chgfilename*80
integer itype
integer :: nlayer=4 !Number of layers of points for MK
real*8 :: espfitvdwr(0:nelesupp)=-1D0,sclvdwlayer(100)=(/1.4D0,1.6D0,1.8D0,2.0D0,(0D0,i=5,100)/)
real*8,allocatable :: ESPptval(:),ESPptx(:),ESPpty(:),ESPptz(:),Bmat(:),Amat(:,:),Amatinv(:,:),atmchg(:)
real*8,allocatable :: fitcenx(:),fitceny(:),fitcenz(:),fitcenvdwr(:),disptcen(:),origsphpt(:,:)

! Read ESP and coordinates of fitting points from Gaussian Iop(6/33=2) output. Corresponding wavefunction file must be loaded to provide atom coordinates
! naddcen=0
! open(10,file="C:\gtest\benzene.out",status="old")
! call loclabel(10,"NAtoms")
! read(10,"(8x,i6)") nfitcen
! allocate(fitcenx(nfitcen),fitceny(nfitcen),fitcenz(nfitcen))
! call loclabel(10,"Atomic Center    1 is at")
! do i=1,nfitcen
!     read(10,"(32x,3f10.6)") fitcenx(i),fitceny(i),fitcenz(i)
! end do
! fitcenx=fitcenx/b2a
! fitceny=fitceny/b2a
! fitcenz=fitcenz/b2a
! call loclabel(10,"points will be used for",ifound)
! read(10,*) nESPpt
! write(*,"('Number of fitting points used:',i10)") nESPpt
! allocate(ESPptval(nESPpt),ESPptx(nESPpt),ESPpty(nESPpt),ESPptz(nESPpt))
! matdim=nfitcen+1
! call loclabel(10,"ESP Fit Center",ifound)
! if (ifound==0) write(*,*) "Cannot locate ""ESP Fit Center"" field"
! do ipt=1,nESPpt
!     read(10,"(32x,3f10.6)") ESPptx(ipt),ESPpty(ipt),ESPptz(ipt)
! end do
! ESPptx=ESPptx/b2a !Convert to Bohr unit
! ESPpty=ESPpty/b2a
! ESPptz=ESPptz/b2a
! call loclabel(10," Fit ",ifound,0)
! if (ifound==0) write(*,*) "Cannot locate "" Fit "" field"
! do ipt=1,nESPpt
!     read(10,"(14x,f10.6)") ESPptval(ipt)
! end do
! ! goto 332
! goto 333

fitspc=0.3D0/b2a !Spacing between grid for CHELPG
extdis=2.8D0/b2a !Extend 2.8 Angstrom to each side for CHELPG
densperarea=6D0*b2a**2 !Point density per Angstrom**2 for MK, in order to convert to Bohr**2, multiply by b2a**2

iaddcen=0 !If give Additional center
iuseextpt=0 !If use external points
iskipespcalc=0 !If read ESP from external file directly rather than calculate here

do while(.true.)
    write(*,*)
    write(*,*) "Note: All units in this module are in a.u."
    write(*,*) "-2 Load additional fitting centers from external file"
    if (iuseextpt==0) write(*,*) "-1 Use fitting points recorded in external file instead of generating them"
    write(*,*) "0 Return"
    write(*,*) "1 Start calculation!"
    if (iuseextpt==0) then
        if (itype==1) then
            write(*,"(' 2 Set grid spacing, current:',f7.3,' Bohr (',f7.3,' Angstrom)')") fitspc,fitspc*b2a
            write(*,"(' 3 Set box extension, current:',f7.3,' Bohr (',f7.3,' Angstrom)')") extdis,extdis*b2a
        else if (itype==2) then
            write(*,"(' 2 Set number of points per Angstrom^2, current:',f10.3)") densperarea/b2a**2 !Temporary convert to Angstrom**2 for convention
            write(*,"(' 3 Set number of layers per atom, current:',i4)") nlayer
            write(*,"(' 4 Set the value times van der Waals radius in each layer')")
        end if
    end if
    read(*,*) isel
    
    if (isel==-2) then
        iaddcen=1
        write(*,*) "Input the name of the file recording coordinates of additional fitting centers"
        read(*,"(a)") addcenfile
        write(*,*) "Done!"
    else if (isel==-1) then
        iuseextpt=1
        write(*,*) "Input the name of the file recording coordinates of ESP fitting points"
        read(*,"(a)") extptfile
        write(*,*) "OK, the points recorded in this file will be used as fitting points"
    else if (isel==0) then
        Return
    else if (isel==1) then
        exit
    else if (isel==2) then
        write(*,*) "Input new value"
        if (itype==1) read(*,*) fitspc
        if (itype==2) then
            read(*,*) densperarea
            densperarea=densperarea*b2a**2
        end if
    else if (isel==3) then
        write(*,*) "Input new value"
        if (itype==1) read(*,*) extdis
        if (itype==2) read(*,*) nlayer
    else if (isel==4.and.itype==2) then
        write(*,*) "Current values:"
        do ilayer=1,nlayer
            write(*,"(' Layer',i4,' :',f8.4)") ilayer,sclvdwlayer(ilayer)
        end do
        write(*,*)
        do ilayer=1,nlayer
            write(*,"(a,i4)") " Input value for layer",ilayer
            read(*,*) sclvdwlayer(ilayer)
        end do
    end if
end do

!Set vdW radius
if (itype==1) then !For CHELPG
    espfitvdwr(1:2)=1.45D0 !vdW radius copied from GetvdW routine (utilam), some of them are given in CHELPG original paper
    espfitvdwr(3:6)=1.5D0
    espfitvdwr(7:10)=1.7D0
    espfitvdwr(11:18)=2D0
    espfitvdwr(1:18)=espfitvdwr(1:18)/b2a !to Bohr
else if (itype==2) then !For MK, copied from GetvdW routine (utilam)
    espfitvdwr(1:17)=(/1.20d0,1.20d0,1.37d0,1.45d0,1.45d0,1.50d0,1.50d0,&
    1.40d0,1.35d0,1.30d0,1.57d0,1.36d0,1.24d0,1.17d0,1.80d0,1.75d0,1.70d0/)
    espfitvdwr(1:17)=espfitvdwr(1:17)/b2a
end if
write(*,*) "Atomic radii used:"
do ielem=1,nelesupp
    if (any(a%index==ielem).and.espfitvdwr(ielem)/=-1D0) write(*,"(' Element:',a,'     vdW radius (Angstrom):',f6.3)") ind2name(ielem),espfitvdwr(ielem)*b2a
end do

!Check sanity and complete vdW radius table for all involved elements
do iatm=1,ncenter
    if (espfitvdwr(a(iatm)%index)==-1D0) then
        write(*,"(' vdW radius used in fitting for element ',a,' is missing, input the radius (Bohr)')") ind2name(a(iatm)%index)
        write(*,"(a)") " Hint: If you don't know how to deal with the problem, simply input 3.4. (However, the radius of 3.4 Bohr may be not very appropriate for current element)" 
        read(*,*) espfitvdwr(a(iatm)%index)
    end if
end do

!Check total number of fitting centers
naddcen=0
if (iaddcen==1) then
    open(10,file=addcenfile,status="old")
    read(10,*) naddcen
end if
nfitcen=ncenter+naddcen
allocate(fitcenx(nfitcen),fitceny(nfitcen),fitcenz(nfitcen),fitcenvdwr(nfitcen),disptcen(nfitcen))

!Generate information of fitting centers
do iatm=1,ncenter
    fitcenx(iatm)=a(iatm)%x
    fitceny(iatm)=a(iatm)%y
    fitcenz(iatm)=a(iatm)%z
    fitcenvdwr(iatm)=espfitvdwr(a(iatm)%index) !vdW radius for each fitting center
end do
if (iaddcen==1) then
    do icen=ncenter+1,ncenter+naddcen
        read(10,*) fitcenx(icen),fitceny(icen),fitcenz(icen)
        fitcenvdwr(icen)=0D0
    end do
    close(10)
end if

write(*,*)
if (iuseextpt==0) then !Count number and generate coordinates of fitting points
    if (itype==1) then !CHELPG
        xlow=minval(fitcenx)-extdis
        xhigh=maxval(fitcenx)+extdis
        ylow=minval(fitceny)-extdis
        yhigh=maxval(fitceny)+extdis
        zlow=minval(fitcenz)-extdis
        zhigh=maxval(fitcenz)+extdis
        xlen=xhigh-xlow
        ylen=yhigh-ylow
        zlen=zhigh-zlow
        nxfit=int(xlen/fitspc)+1
        nyfit=int(ylen/fitspc)+1
        nzfit=int(zlen/fitspc)+1
        nESPpt=0
        do ix=0,nxfit
            do iy=0,nyfit
                do iz=0,nzfit
                    tmpx=xlow+ix*fitspc
                    tmpy=ylow+iy*fitspc
                    tmpz=zlow+iz*fitspc
                    do icen=1,nfitcen
                        disptcen(icen)=dsqrt( (fitcenx(icen)-tmpx)**2+(fitceny(icen)-tmpy)**2+(fitcenz(icen)-tmpz)**2 )
                        if (disptcen(icen)<=fitcenvdwr(icen)) exit
                        if (icen==nfitcen.and.any(disptcen<=extdis)) nESPpt=nESPpt+1
                    end do
                end do
            end do
        end do
        write(*,"(' Number of fitting points used:',i10)") nESPpt
        allocate(ESPptval(nESPpt),ESPptx(nESPpt),ESPpty(nESPpt),ESPptz(nESPpt))
        iESPpt=0
        do ix=0,nxfit
            do iy=0,nyfit
                do iz=0,nzfit
                    tmpx=xlow+ix*fitspc
                    tmpy=ylow+iy*fitspc
                    tmpz=zlow+iz*fitspc
                    do icen=1,nfitcen
                        disptcen(icen)=dsqrt( (fitcenx(icen)-tmpx)**2+(fitceny(icen)-tmpy)**2+(fitcenz(icen)-tmpz)**2 )
                        if (disptcen(icen)<=fitcenvdwr(icen)) exit
                        if (icen==nfitcen.and.any(disptcen<=extdis)) then
                            iESPpt=iESPpt+1
                            ESPptx(iESPpt)=tmpx
                            ESPpty(iESPpt)=tmpy
                            ESPptz(iESPpt)=tmpz
                        end if
                    end do
                end do
            end do
        end do
        
    else if (itype==2) then !MK, don't reserve spatial space for additional center
        cutinnerscl=minval(sclvdwlayer(1:nlayer))
        write(*,"(' Note: If distance between a ESP point and any atom is smaller than',f6.3,' multiplied by corresponding vdW radius, then the point will be discarded')") cutinnerscl
        nESPpt=0
        maxsphpt=nint(4D0*pi*(maxval(fitcenvdwr)*maxval(sclvdwlayer))**2 *densperarea) !Find maximal possible number of points in unit sphere to allocate temporary origsphpt
        allocate(origsphpt(3,maxsphpt))
        do icen=1,ncenter !Rather than nfitcen.   Count how many possible ESP points in total
            do ilayer=1,nlayer
                numsphpt=nint(4D0*pi*(fitcenvdwr(icen)*sclvdwlayer(ilayer))**2 *densperarea)
                nESPpt=nESPpt+numsphpt
            end do
        end do
        allocate(ESPptval(nESPpt),ESPptx(nESPpt),ESPpty(nESPpt),ESPptz(nESPpt)) !Currently nESPpt is upper limit
        iESPpt=0
        do icen=1,ncenter
            do ilayer=1,nlayer
                radius=fitcenvdwr(icen)*sclvdwlayer(ilayer)
                numsphpt=nint(4D0*pi*radius**2 *densperarea)
                call unitspherept(origsphpt,numsphpt) !Input expected number of point in unit sphere, return actual number of points
                origsphpt(:,1:numsphpt)=origsphpt(:,1:numsphpt)*radius
                origsphpt(1,1:numsphpt)=origsphpt(1,1:numsphpt)+fitcenx(icen) !Move unit sphere to atomic center
                origsphpt(2,1:numsphpt)=origsphpt(2,1:numsphpt)+fitceny(icen)
                origsphpt(3,1:numsphpt)=origsphpt(3,1:numsphpt)+fitcenz(icen)
                do ipt=1,numsphpt
                    tmpx=origsphpt(1,ipt)
                    tmpy=origsphpt(2,ipt)
                    tmpz=origsphpt(3,ipt)
                    iok=1
                    do icen2=1,ncenter
                        if (icen2==icen) cycle
                        disptcensq=(fitcenx(icen2)-tmpx)**2+(fitceny(icen2)-tmpy)**2+(fitcenz(icen2)-tmpz)**2 !distance between point and center
                        if (disptcensq<(fitcenvdwr(icen2)*cutinnerscl)**2) then !Less than vdW RADIUS*cutinner of atom icen2, it should be ommitted
                            iok=0
                            exit
                        end if
                    end do
                    if (iok==1) then
                        iESPpt=iESPpt+1
                        ESPptx(iESPpt)=tmpx
                        ESPpty(iESPpt)=tmpy
                        ESPptz(iESPpt)=tmpz
                    end if
                end do
            end do
        end do
        nESPpt=iESPpt
        deallocate(origsphpt)
        write(*,"(' Number of fitting points used:',i10)") nESPpt
    end if
    
else if (iuseextpt==1) then !Directly use external fitting points
    open(10,file=extptfile,status="old")
    read(10,*) nESPpt
    if (nESPpt<0) then
        iskipespcalc=1 !If the number of fitting points is negative, that means the fourth column records ESP value and needn't to be recalculated
        write(*,*) "ESP value of all fitting points are read from external file directly"
    end if
    nESPpt=abs(nESPpt)
    write(*,"(' Number of fitting points used:',i10)") nESPpt
    allocate(ESPptval(nESPpt),ESPptx(nESPpt),ESPpty(nESPpt),ESPptz(nESPpt))
    do i=1,nESPpt
        if (iskipespcalc==0) read(10,*) ESPptx(i),ESPpty(i),ESPptz(i)
        if (iskipespcalc==1) read(10,*) ESPptx(i),ESPpty(i),ESPptz(i),ESPptval(i)
    end do
    close(10)
end if

332 continue
!Generate ESP value of fitting points
if (iskipespcalc==0) then
    write(*,*) "Calculating ESP, please wait..."
    itmp=1
    do ipt=1,nESPpt
        if (ipt>=itmp*300) then
            write(*,"(' Finished:',i10,'  /',i10)") ipt,nESPpt
            itmp=itmp+1
        end if
        ESPptval(ipt)=totesp((ESPptx(ipt)),(ESPpty(ipt)),(ESPptz(ipt)))
    end do
    write(*,*) "Done!"
end if

333 continue !We just read ESP points from Gaussian output directly, so jump to here
matdim=nfitcen+1
allocate(Bmat(matdim),Amat(matdim,matdim),Amatinv(matdim,matdim),atmchg(matdim))
!See original paper of MK for detail of algorithem
!Forming Amat
Amat=0D0
do icen=1,nfitcen
    do jcen=icen,nfitcen
        do ipt=1,nESPpt
            dis1=dsqrt( (ESPptx(ipt)-fitcenx(icen))**2 + (ESPpty(ipt)-fitceny(icen))**2 + (ESPptz(ipt)-fitcenz(icen))**2 )
            dis2=dsqrt( (ESPptx(ipt)-fitcenx(jcen))**2 + (ESPpty(ipt)-fitceny(jcen))**2 + (ESPptz(ipt)-fitcenz(jcen))**2 )
            Amat(icen,jcen)=Amat(icen,jcen)+1D0/dis1/dis2
        end do
    end do
end do
Amat=Amat+transpose(Amat)
do i=1,nfitcen
    Amat(i,i)=Amat(i,i)/2D0
end do
Amat(matdim,:)=1D0
Amat(:,matdim)=1D0
Amat(matdim,matdim)=0D0
!Forming Bmat
Bmat=0D0
do icen=1,nfitcen
    do ipt=1,nESPpt
        dis=dsqrt( (ESPptx(ipt)-fitcenx(icen))**2 + (ESPpty(ipt)-fitceny(icen))**2 + (ESPptz(ipt)-fitcenz(icen))**2 )
        Bmat(icen)=Bmat(icen)+ESPptval(ipt)/dis
    end do
end do
Bmat(matdim)=sum(a(:)%charge)-nelec !Net charge
Amatinv=invmat(Amat,matdim)
atmchg=matmul(Amatinv,Bmat)

!Output summary
write(*,*) " Center        X           Y           Z            Charge"
do i=1,ncenter
    write(*,"(i6,a,3f12.6,f16.6)") i,ind2name(a(i)%index),fitcenx(i),fitceny(i),fitcenz(i),atmchg(i)
end do
do i=ncenter+1,ncenter+naddcen
    write(*,"(i6,2x,3f12.6,f16.6)") i,fitcenx(i),fitceny(i),fitcenz(i),atmchg(i)
end do
write(*,"(' Sum of charges:',f12.6)") sum(atmchg(1:nfitcen))
!Calculate RMSE and RRMSE
RMSE=0D0
do ipt=1,nESPpt
    atmchgesp=0D0
    do icen=1,nfitcen
        dis=dsqrt( (ESPptx(ipt)-fitcenx(icen))**2 + (ESPpty(ipt)-fitceny(icen))**2 + (ESPptz(ipt)-fitcenz(icen))**2 )
        atmchgesp=atmchgesp+atmchg(icen)/dis
    end do
    RMSE=RMSE+(ESPptval(ipt)-atmchgesp)**2
end do
RRMSE=dsqrt(RMSE/sum(ESPptval(1:nESPpt)**2))
RMSE=dsqrt(RMSE/nESPpt)
write(*,"(' RMSE:',f12.6,'   RRMSE:',f12.6)") RMSE,RRMSE

!Show fragment charge
if (allocated(frag1)) write(*,"(/,' Fragment charge:',f12.6)") sum(atmchg(frag1))

write(*,*)
write(*,"(a)") " If output coordinates and ESP value of all fitting points to ESPfitpt.txt in current folder? (y/n)"
read(*,*) selectyn
if (selectyn=='y'.or.selectyn=="Y") then
    open(10,file="ESPfitpt.txt",status="replace")
    write(10,*) nESPpt
    do ipt=1,nESPpt
        write(10,"(3f12.6,f14.8)") ESPptx(ipt),ESPpty(ipt),ESPptz(ipt),ESPptval(ipt)
    end do
    write(*,*) "Data have been outputted to ESPfitpt.txt in current folder"
    write(*,"(a)") " All units are in a.u. The first line shows the number of fitting points, the first three columns are X,Y,Z coordinates, the last column corresponds to ESP value"
    close(10)
end if

call path2filename(firstfilename,chgfilename)
write(*,"(a)") " If output atom coordinates with charges to "//trim(chgfilename)//".chg file in current folder? (y/n)"
read(*,*) selectyn
if (selectyn=='y'.or.selectyn=="Y") then
    open(10,file=trim(chgfilename)//".chg",status="replace")
    do i=1,ncenter
        write(10,"(a4,4f12.6)") a(i)%name,a(i)%x*b2a,a(i)%y*b2a,a(i)%z*b2a,atmchg(i)
    end do
    do i=ncenter+1,ncenter+naddcen
        write(10,"(a4,4f12.6)") " Bq",a(i)%x*b2a,a(i)%y*b2a,a(i)%z*b2a,atmchg(i)
    end do
    close(10)
    write(*,"(a)") " Result have been saved to "//trim(chgfilename)//".chg in current folder"
    write(*,"(a)") " Columns from 1 to 5 are name,X,Y,Z,charge respectively, unit is Angstrom"
end if

!Output fitting points to pdb file for visualizing, only for debug purpose
! open(10,file="ESPfitpt.pdb",status="replace")
! do ipt=1,nESPpt
!     write(10,"(a6,i5,1x,a4,1x,a3, 1x,a1,i4,4x,3f8.3,2f6.2,10x,a2)") &
!                     "HETATM",i,' '//"O "//' ',"MOL",'A',1,ESPptx(ipt),ESPpty(ipt),ESPptz(ipt),1.0,ESPptval(ipt)*au2kcal,"O "
! end do
! write(*,*) "Data have been outputted to ESPfitpt.pdb in current folder"
! close(10)

end subroutine


!!!------------- Generate numpt points scattered evenly in an unit sphere
! Input argument numpt is the expected number of points, while the return value is actual number
! ptcrd store coordinates of the points 
subroutine unitspherept(ptcrd,numpt)
implicit real*8 (a-h,o-z)
real*8 ptcrd(3,numpt)
integer numpt
pi=3.141592653589793D0
!The average number of equator points in all XY layes is numequ*2/pi, and there are numvert=numequ/2 layers
!Solve (numequ*2/pi)*numequ/2=numpt one can get numequ=sqrt(numpt*pi)
numequ=int(sqrt(numpt*pi)) !Maximal number of point in each XY layer
numvert=numequ/2
ipt=0
do ivert=0,numvert
    angz=dfloat(ivert)/numvert*pi
    scalexy=sin(angz)
    z=cos(angz)
    numxy=int(numequ*scalexy)
    if (numxy==0) numxy=1
    do ihori=1,numxy
        ipt=ipt+1
        if (ipt>numpt) then
            numpt=ipt-1
            return
        end if
        angxy=2D0*pi*ihori/numxy
        ptcrd(1,ipt)=cos(angxy)*scalexy
        ptcrd(2,ipt)=sin(angxy)*scalexy
        ptcrd(3,ipt)=z
    end do
end do
numpt=ipt
end subroutine
