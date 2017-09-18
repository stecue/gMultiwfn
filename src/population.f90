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
        !Not available because integration error of below two methods by means of Becke integration are too large
    !         write(*,*) "3 Integrate electron density in voronoi cell"
    !         write(*,*) "4 Adjusted method 3 by Rousseau et al."
        if (allocated(cobasa)) then
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
        write(*,*) "15 Hirshfeld-I charge"
        write(*,*) "16 CM5 charge"
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
        else if (ipopsel==15) then
            if (iautointgrid==1) then
                radpot=30
                sphpot=170
                if (any(a%index>18)) radpot=40
                if (any(a%index>36)) radpot=50
                if (any(a%index>54)) radpot=60
            end if
            call Hirshfeld_I_wrapper(1)
        else if (ipopsel==16) then
            call spacecharge(7)
        end if
        if (imodwfnold==1.and.(ipopsel==1.or.ipopsel==2.or.ipopsel==6.or.ipopsel==11)) then !1,2,6,11 are the methods need to reload the initial wavefunction
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
character selectyn,chgfilename*200
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
character selectyn,chgfilename*200
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
character selectyn,corbnum*6,cOcc*12,chgfilename*200
integer isel

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
                    write(*,"(' Shell',i6,' Type: ',a,'    in atom',i5,'(',a,') :',f9.5)") ish,shtype2name(shtype(ish)),iatm,a(iatm)%name,shellpop
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
                    write(*,"(i5,5x,a,i7,'(',a,')' ,4f12.5)") ish,shtype2name(shtype(ish)),iatm,a(iatm)%name,shellpopa,shellpopb,shellpopa+shellpopb,shellpopa-shellpopb
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
!1=Hirshfeld, 2=VDD, 3=Integrate electron density in voronoi cell
!4=Adjusted method 3 by Rousseau et al., 5= Becke with/without ADC, 6= ADCH, 7= CM5
subroutine spacecharge(chgtype)
use defvar
use function
use util
implicit real*8(a-h,o-z)
integer chgtype
real*8 molrho(radpot*sphpot),promol(radpot*sphpot),tmpdens(radpot*sphpot),beckeweigrid(radpot*sphpot),selfdens(radpot*sphpot)
type(content) gridatm(radpot*sphpot),gridatmorg(radpot*sphpot)
real*8 atmdipx(ncenter),atmdipy(ncenter),atmdipz(ncenter),charge(ncenter)
real*8 :: covr_becke(0:nelesupp) !covalent radii used for Becke population
character selectyn,chgfilename*200
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

!***** 1=Hirshfeld, 2=VDD, 6=ADCH, 7=CM5
if (chgtype==1.or.chgtype==2.or.chgtype==6.or.chgtype==7) then
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
nthreads=getNThreads()
!$OMP parallel do shared(molrho) private(i) num_threads(nthreads)
        do i=1,radpot*sphpot
            molrho(i)=fdens(gridatm(i)%x,gridatm(i)%y,gridatm(i)%z)
        end do
!$OMP end parallel do
        !Calc free atomic density to obtain promolecule density
        promol=0D0
        if (iatmdensmode==1) then
            do jatm=1,ncenter
nthreads=getNThreads()
!$OMP parallel do shared(tmpdens) private(ipt) num_threads(nthreads)
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
nthreads=getNThreads()
!$OMP parallel do shared(tmpdens) private(ipt) num_threads(nthreads)
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
        if (chgtype==1.or.chgtype==6.or.chgtype==7) then !Hirshfeld, ADCH charge, CM5 charge
            do i=1,radpot*sphpot
                if (promol(i)/=0D0) then
                    tmpv=selfdens(i)/promol(i)*molrho(i)*gridatm(i)%value
                    tmpcharge=tmpcharge-tmpv
                    dipx=dipx-(gridatm(i)%x-atmx)*tmpv
                    dipy=dipy-(gridatm(i)%y-atmy)*tmpv
                    dipz=dipz-(gridatm(i)%z-atmz)*tmpv
                end if
            end do
            if (nEDFelec==0) charge(iatm)=a(iatm)%charge+tmpcharge
            if (nEDFelec>0) charge(iatm)=a(iatm)%index+tmpcharge !EDF is provided
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
        if (chgtype==1.or.chgtype==6.or.chgtype==7) write(*,"(' Hirshfeld charge of atom ',i5,'(',a2,')',' is',f12.6)") iatm,a_org(iatm)%name,charge(iatm)
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
nthreads=getNThreads()
!$OMP parallel do shared(tmpdens) private(i) num_threads(nthreads)
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
        if (nEDFelec==0) charge(iatm)=a(iatm)%charge+tmpcharge
        if (nEDFelec>0) charge(iatm)=a(iatm)%index+tmpcharge !EDF is provided
        atmdipx(iatm)=dipx
        atmdipy(iatm)=dipy
        atmdipz(iatm)=dipz
        write(*,"(' Becke charge of atom ',i5,'(',a2,')',' is',f12.6)") iatm,a(iatm)%name,charge(iatm)
    end do
end if

write(*,"(' Summing up all charges:',f15.8)") sum(charge)
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
totdip=dsqrt(xmoldip**2+ymoldip**2+zmoldip**2)
write(*,"(' Total dipole moment from atomic charges:',f12.6,' a.u.')") totdip
write(*,"(' X/Y/Z of dipole from atomic charge:',3f12.6,' a.u.')") xmoldip,ymoldip,zmoldip
totatmdip=dsqrt(sum(atmdipx)**2+sum(atmdipy)**2+sum(atmdipz)**2)
write(*,"(' Total atomic dipole moment:',f12.6,' a.u.')") totatmdip
write(*,"(' X/Y/Z of total atomic dipole:',3f12.6,' a.u.')") sum(atmdipx),sum(atmdipy),sum(atmdipz)
corrdipx=xmoldip+sum(atmdipx) !Corresponding to actual molecular dipole moment derived from molecular density
corrdipy=ymoldip+sum(atmdipy)
corrdipz=zmoldip+sum(atmdipz)
realdip=dsqrt(corrdipx**2+corrdipy**2+corrdipz**2)

if (chgtype==5) call doADC(atmdipx,atmdipy,atmdipz,charge,realdip,5)
if (chgtype==6) call doADC(atmdipx,atmdipy,atmdipz,charge,realdip,6)
if (chgtype==7) call doCM5(charge)

!Show fragment charge
if (allocated(frag1)) write(*,"(/,' Fragment charge:',f12.6)") sum(charge(frag1))

write(*,*)
call walltime(nwalltime2)
write(*,"(' Calculation took up',i8,' seconds wall clock time')")  nwalltime2-nwalltime1

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
!The "charge" is inputted Hirshfeld charge, finally it is replaced by ADC charge 
!chgtype 5= Becke with/without ADC, 6= ADCH
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
    if (ishowchgtrans==1) write(*,"('Atom: 'i4,a)") i,a(i)%name !ishowchgtrans==1 means output detail of charge transferation process during atomic dipole moment correction
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
        if (ishowchgtrans==1) write(*,"(' Give atom ',i4,a4,f15.12,'  Weight',2f15.12)") j,a(j)%name,tmp(1,1),w(j)
        totq=totq+tmp(1,1)
        tottmpdipx=tottmpdipx+(a(j)%x-a(i)%x)*tmp(1,1)
        tottmpdipy=tottmpdipy+(a(j)%y-a(i)%y)*tmp(1,1)
        tottmpdipz=tottmpdipz+(a(j)%z-a(i)%z)*tmp(1,1)
    end do
    if (ishowchgtrans==1) write(*,*)
end do

write(*,*) "   ======= Summary of atomic dipole moment corrected (ADC) charges ======="
do i=1,ncenter
    write(*,"(' Atom: ',i4,a,'  Corrected charge:',f12.6,'  Before:',f12.6)") i,a(i)%name,chargecorr(i),charge(i)
end do
write(*,"(' Summing up all corrected charges:',f12.7)") sum(chargecorr)
if (chgtype==5) write(*,"(a)") " Note: The values shown after ""Corrected charge"" are atomic dipole moment corrected Becke charges, the ones after ""Before"" are normal Becke charges"
if (chgtype==6) write(*,"(a)") " Note: The values shown after ""Corrected charge"" are ADCH charges, the ones after ""Before"" are Hirshfeld charges"
ADCdipx=sum(a%x*chargecorr)
ADCdipy=sum(a%y*chargecorr)
ADCdipz=sum(a%z*chargecorr)
ADCdip=sqrt(ADCdipx**2+ADCdipy**2+ADCdipz**2)
write(*,*)
write(*,"(' Total dipole from ADC charges (a.u.)',f11.7,'  Error:',f11.7)") ADCdip,abs(ADCdip-realdip)
write(*,"(' X/Y/Z of dipole moment from the charge (a.u.)',3f11.7)") ADCdipx,ADCdipy,ADCdipz
charge=chargecorr !Overlay charge array, then return to Hirshfeld module and output result to .chg file
end subroutine


!!--------- Calculate CM5 charge based on Hirshfeld charge
subroutine doCM5(charge)
use defvar
implicit real*8 (a-h,o-z)
real*8 charge(ncenter),CMcharge(ncenter),radius(118),Dparm(118)
alpha=2.474D0
Dparm=0D0
Dparm(1)=0.0056D0
Dparm(2)=-0.1543D0
Dparm(4)=0.0333D0
Dparm(5)=-0.1030D0
Dparm(6)=-0.0446D0
Dparm(7)=-0.1072D0
Dparm(8)=-0.0802D0
Dparm(9)=-0.0629D0
Dparm(10)=-0.1088D0
Dparm(11)=0.0184D0
Dparm(13)=-0.0726D0
Dparm(14)=-0.0790D0
Dparm(15)=-0.0756D0
Dparm(16)=-0.0565D0
Dparm(17)=-0.0444D0
Dparm(18)=-0.0767D0
Dparm(19)=0.0130D0
Dparm(31)=-0.0512D0
Dparm(32)=-0.0557D0
Dparm(33)=-0.0533D0
Dparm(34)=-0.0399D0
Dparm(35)=-0.0313D0
Dparm(36)=-0.0541D0
Dparm(37)=0.0092D0
Dparm(49)=-0.0361D0
Dparm(50)=-0.0393D0
Dparm(51)=-0.0376D0
Dparm(52)=-0.0281D0
Dparm(53)=-0.0220D0
Dparm(54)=-0.0381D0
Dparm(55)=0.0065D0
Dparm(81)=-0.0255D0
Dparm(82)=-0.0277D0
Dparm(83)=-0.0265D0
Dparm(84)=-0.0198D0
Dparm(85)=-0.0155D0
Dparm(86)=-0.0269D0
Dparm(87)=0.0046D0
Dparm(113)=-0.0179D0
Dparm(114)=-0.0195D0
Dparm(115)=-0.0187D0
Dparm(116)=-0.0140D0
Dparm(117)=-0.0110D0
Dparm(118)=-0.0189D0
!As shown in CM5 paper, the covalent radii used in CM5 equation are tabulated in CRC book 91th, where they are obtained as follows:
!For Z=1~96, the radii are the average of CSD radii (For Fe, Mn, Co the low-spin is used) and Pyykko radii
!For Z=97~118, the radii are Pyykko radii
radius(1:96)=(covr(1:96)+covr_pyy(1:96))/2D0
radius(97:118)=covr_pyy(97:118)
radius=radius*b2a !Because the radii have already been converted to Bohr, so we convert them back to Angstrom

if (ishowchgtrans==1) write(*,"(/,a)") " Details of CM5 charge correction:"

do iatm=1,ncenter
    CMcorr=0
    iZ=a(iatm)%index
    do jatm=1,ncenter
        if (iatm==jatm) cycle
        jZ=a(jatm)%index
        Bval=exp( -alpha*(distmat(iatm,jatm)*b2a-radius(iZ)-radius(jZ)) )
        if (iZ==1.and.jZ==6) then
            Tval=0.0502D0
        else if (iZ==6.and.jZ==1) then
            Tval=-0.0502D0
        else if (iZ==1.and.jZ==7) then
            Tval=0.1747D0
        else if (iZ==7.and.jZ==1) then
            Tval=-0.1747D0
        else if (iZ==1.and.jZ==8) then
            Tval=0.1671D0
        else if (iZ==8.and.jZ==1) then
            Tval=-0.1671D0
        else if (iZ==6.and.jZ==7) then
            Tval=0.0556D0
        else if (iZ==7.and.jZ==6) then
            Tval=-0.0556D0
        else if (iZ==6.and.jZ==8) then
            Tval=0.0234D0
        else if (iZ==8.and.jZ==6) then
            Tval=-0.0234D0
        else if (iZ==7.and.jZ==8) then
            Tval=-0.0346D0
        else if (iZ==8.and.jZ==7) then
            Tval=0.0346D0
        else
            Tval=Dparm(iZ)-Dparm(jZ)
        end if
        CMcorr=CMcorr+Tval*Bval
        if (ishowchgtrans==1) then
            write(*,"(i4,a,i4,a,'  B_term:',f10.5,'  T_term:',f10.5,'  Corr. charge:',f10.5)") iatm,a(iatm)%name,jatm,a(jatm)%name,Bval,Tval,Tval*Bval
        end if
    end do
    CMcharge(iatm)=charge(iatm)+CMcorr
end do
write(*,*)
write(*,*) "                    ======= Summary of CM5 charges ======="
do i=1,ncenter
    write(*,"(' Atom: ',i4,a,'  CM5 charge:',f12.6,'  Hirshfeld charge:',f12.6)") i,a(i)%name,CMcharge(i),charge(i)
end do
write(*,"(' Summing up all CM5 charges:',f15.8)") sum(CMcharge)
CM5dipx=sum(a%x*CMcharge)
CM5dipy=sum(a%y*CMcharge)
CM5dipz=sum(a%z*CMcharge)
CM5dip=sqrt(CM5dipx**2+CM5dipy**2+CM5dipz**2)
write(*,*)
write(*,"(' Total dipole moment from CM5 charges',f12.7,' a.u.')") CM5dip
write(*,"(' X/Y/Z of dipole moment from CM5 charges',3f10.5, ' a.u.')") CM5dipx,CM5dipy,CM5dipz
charge=chargecorr
end subroutine



!!------------ Calculate charges by fitting ESP, currently CHELPG grid and MK grid are used
!itype=1:CHELPG   itype=2:MK
subroutine fitESP(itype)
use util
use defvar
use function
implicit real*8 (a-h,o-z)
character*200 addcenfile,extptfile
character selectyn,chgfilename*200
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








!!============================ Hirshfeld-I ============================!!
!!============================ Hirshfeld-I ============================!!
!!============================ Hirshfeld-I ============================!!
!!============================ Hirshfeld-I ============================!!
!!============================ Hirshfeld-I ============================!!
!Wrapper of Hirshfeld-I module to automatically set radpot and sphpot to proper values
!itype=1: Normal population analysis =2: Only used to generate proper atomic space (i.e. Don't do unnecessary things)
subroutine Hirshfeld_I_wrapper(itype)
use defvar
implicit real*8 (a-h,o-z)
radpotold=radpot
sphpotold=sphpot
if (iautointgrid==1) then
    radpot=30
    sphpot=170
    if (any(a%index>18)) radpot=40
    if (any(a%index>36)) radpot=50
    if (any(a%index>54)) radpot=60
end if
call Hirshfeld_I(itype)
if (iautointgrid==1) then
    radpot=radpotold
    sphpot=sphpotold
end if
end subroutine

!!--------- Calculate Hirshfeld-I charge and yield final atomic radial density
!I've compared this module with hipart, this module is faster than hipart, and the accuracy under default setting is at least never lower than hipart
subroutine Hirshfeld_I(itype)
use defvar
use function
use util
implicit real*8 (a-h,o-z)
integer itype
type(content) gridatm(radpot*sphpot),gridatmorg(radpot*sphpot)
real*8 molrho(radpot*sphpot),promol(radpot*sphpot),tmpdens(radpot*sphpot),selfdens(radpot*sphpot),molrhoall(ncenter,radpot*sphpot)
real*8 charge(ncenter),lastcharge(ncenter) !Atomic charge of current iter. and last iter.
real*8 radrholow(200),radrhohigh(200)
character sep,c80tmp*80,chgfilename*200,selectyn
character*2 :: statname(-4:4)=(/ "-4","-3","-2","-1","_0","+1","+2","+3","+4" /)
integer :: maxcyc=50,ioutmedchg=0
real*8 :: crit=0.0002D0
real*8,external :: fdens_rad
!Used for mode 2. e.g. atmstatgrid(iatm,igrid,jatm,-1) means density of jatm with -1 charge state at igrid around iatm
real*8,allocatable :: atmstatgrid(:,:,:,:)
ntotpot=radpot*sphpot

!Mode 1 use very low memory but expensive, because most data is computed every iteration
!Mode 2 use large memory but fast, because most data is only computed once at initial stage
!The result of the two modes differ with each other marginally, probably because in mode 1 radial density is related to max(npthigh,nptlow), which is not involved in mode 2
!In principle, result of mode 2 is slightly better
imode=2

!Ignore jatm contribution to iatm centered grids if distance between iatm and jatm is larger than 1.5 times of sum of their vdwr
!This can decrease lots of time for large system, the lose of accuracy can be ignored (error is ~0.0001 per atom)
ignorefar=1

write(*,*)
if (itype==1) write(*,*) "     =============== Iterative Hirshfeld (Hirshfeld-I) ==============="
if (itype==2) write(*,*) "     ============== Generate Hirshfeld-I atomic weights =============="
do while(.true.)
    if (imode==1) write(*,*) "-2 Switch algorithm, current: Slow & low memory requirement"
    if (imode==2) write(*,*) "-2 Switch algorithm, current: Fast & large memory requirement"
    if (itype==1) then
        if (ioutmedchg==0) write(*,*) "-1 Switch if output intermediate results, current: No"
        if (ioutmedchg==1) write(*,*) "-1 Switch if output intermediate results, current: Yes"
        write(*,*) "0 Return"
    end if
    write(*,*) "1 Start calculation!"
    write(*,"(a,i4)") " 2 Set the maximum number of iterations, current:",maxcyc
    write(*,"(a,f10.6)") " 3 Set convergence criterion of atomic charges, current:",crit
    read(*,*) isel
    if (isel==-2) then
        if (imode==1) then
            imode=2
        else
            imode=1
            crit=0.001 !mode 1 is more time-consuming, use loose criterion
        end if
    else if (isel==-1) then
        if (ioutmedchg==1) then
            ioutmedchg=0
        else
            ioutmedchg=1
        end if
    else if (isel==0) then
        return
    else if (isel==1) then
        exit
    else if (isel==2) then
        write(*,*) "Input maximum number of iterations, e.g. 30"
        read(*,*) maxcyc
    else if (isel==3) then
        write(*,*) "Input convergence criterion of atomic charges, e.g. 0.001"
        read(*,*) crit
    end if
end do

!Generate all needed .rad files
call genatmradfile

!====== Start calculation ======!
call walltime(iwalltime1)
CALL CPU_TIME(time_begin)

!Currently all atoms share the same radial points
nradpt=200
itmp=0
do irad=nradpt,1,-1
    radx=cos(irad*pi/(nradpt+1))
    itmp=itmp+1
    atmradpos(itmp)=(1+radx)/(1-radx)
end do

!Generate single center integration grid
call gen1cintgrid(gridatmorg,iradcut)
write(*,*)
write(*,"(' Radial grids:',i4,'  Angular grids:',i5,'  Total:',i7,'  After pruning:',i7)") radpot,sphpot,radpot*sphpot,radpot*sphpot-iradcut*sphpot

!Calculate molecular density
write(*,*) "Calculating density of actual molecule for all grids..."
nthreads=getNThreads()
!$OMP parallel do shared(molrhoall) private(iatm,ipt,gridatm) num_threads(nthreads)
do iatm=1,ncenter
    gridatm%x=gridatmorg%x+a(iatm)%x !Move quadrature point to actual position in molecule
    gridatm%y=gridatmorg%y+a(iatm)%y
    gridatm%z=gridatmorg%z+a(iatm)%z
    do ipt=1+iradcut*sphpot,ntotpot
        molrhoall(iatm,ipt)=fdens(gridatm(ipt)%x,gridatm(ipt)%y,gridatm(ipt)%z)
    end do
end do
!$OMP end parallel do

if (allocated(atmradnpt)) deallocate(atmradnpt)
if (allocated(atmradrho)) deallocate(atmradrho)
allocate(atmradnpt(ncenter),atmradrho(ncenter,200))
sep='/' !Separation symbol of directory
if (isys==1) sep='\'

!Calculate contribution of all atoms in every state to each atomic centered grids
if (imode==2) then
    allocate(atmstatgrid(ncenter,ntotpot,ncenter,-2:2))
    atmstatgrid=0
    write(*,*) "Calculating atomic density contribution to grids..."
    do iatm=1,ncenter !The center of grids
        write(*,"(' Progress:',i5,' /',i5)") iatm,ncenter
        gridatm%value=gridatmorg%value !Weight in this grid point
        gridatm%x=gridatmorg%x+a(iatm)%x !Move quadrature point to actual position in molecule
        gridatm%y=gridatmorg%y+a(iatm)%y
        gridatm%z=gridatmorg%z+a(iatm)%z
        do istat=-2,2 !Charge state
            do jatm=1,ncenter
                if (ignorefar==1) then
                    atmdist=dsqrt( (a(iatm)%x-a(jatm)%x)**2+(a(iatm)%y-a(jatm)%y)**2+(a(iatm)%z-a(jatm)%z)**2 )
                    if (atmdist>(vdwr(iatm)+vdwr(jatm))*1.5D0) cycle
                end if
                if (a(jatm)%index==1.and.istat==1) cycle !H+ doesn't contains electron and cannot compute density
                c80tmp="atmrad"//sep//trim(a(jatm)%name)//statname(istat)//".rad"
                inquire(file=c80tmp,exist=alive)
                if (alive.eqv. .false.) cycle
                open(10,file=c80tmp,status="old")
                read(10,*) atmradnpt(jatm)
                do ipt=1,atmradnpt(jatm)
                    read(10,*) rnouse,atmradrho(jatm,ipt)
                end do
                close(10)
                do ipt=1+iradcut*sphpot,ntotpot
                    atmstatgrid(iatm,ipt,jatm,istat)=fdens_rad(jatm,gridatm(ipt)%x,gridatm(ipt)%y,gridatm(ipt)%z)
                end do
            end do
        end do
    end do
end if

!Set atomic initial radial density as neutral state, which is loaded from corresponding .rad file
atmradrho=0
do iatm=1,ncenter
    open(10,file="atmrad"//sep//trim(a(iatm)%name)//"_0.rad",status="old")
    read(10,*) atmradnpt(iatm)
    do ipt=1,atmradnpt(iatm)
        read(10,*) rnouse,atmradrho(iatm,ipt)
    end do
    close(10)
end do

write(*,*)
write(*,*) "Performing Hirshfeld-I iteration to refine atomic spaces..."
lastcharge=0
!Cycle each atom to calculate their charges
do icyc=1,maxcyc
    if (ioutmedchg==1) write(*,*)
    if (icyc==1) then
        write(*,"(' Cycle',i5)") icyc
    else
        write(*,"(' Cycle',i5,'   Maximum change:',f10.6)") icyc,varmax
    end if
    
    do iatm=1,ncenter
        gridatm%value=gridatmorg%value !Weight in this grid point
        gridatm%x=gridatmorg%x+a(iatm)%x !Move quadrature point to actual position in molecule
        gridatm%y=gridatmorg%y+a(iatm)%y
        gridatm%z=gridatmorg%z+a(iatm)%z
        
        !Molecular density
        molrho=molrhoall(iatm,:)
        
        !Calculate promolecular and proatomic density 
        promol=0D0
        do jatm=1,ncenter
            if (ignorefar==1) then
                atmdist=dsqrt( (a(iatm)%x-a(jatm)%x)**2+(a(iatm)%y-a(jatm)%y)**2+(a(iatm)%z-a(jatm)%z)**2 )
                if (atmdist>(vdwr(iatm)+vdwr(jatm))*1.5D0) cycle
            end if
            if (imode==1) then
nthreads=getNThreads()
!$OMP parallel do shared(tmpdens) private(ipt) num_threads(nthreads)
                do ipt=1+iradcut*sphpot,ntotpot
                    tmpdens(ipt)=fdens_rad(jatm,gridatm(ipt)%x,gridatm(ipt)%y,gridatm(ipt)%z)
                end do
!$OMP end parallel do
            else if (imode==2) then
                if (icyc==1) then
                    tmpdens=atmstatgrid(iatm,:,jatm,0)
                else
                    ichglow=floor(lastcharge(jatm))    
                    ichghigh=ceiling(lastcharge(jatm))
                    tmpdens=(lastcharge(jatm)-ichglow)*atmstatgrid(iatm,:,jatm,ichghigh) + (ichghigh-lastcharge(jatm))*atmstatgrid(iatm,:,jatm,ichglow)
                end if
            end if
            promol=promol+tmpdens
            if (jatm==iatm) selfdens=tmpdens
        end do
        
        !Calculate atomic charge
        electmp=0D0
        do ipt=1+iradcut*sphpot,ntotpot
            if (promol(ipt)/=0D0) electmp=electmp+selfdens(ipt)/promol(ipt)*molrho(ipt)*gridatm(ipt)%value
        end do
        if (nEDFelec==0) charge(iatm)=a(iatm)%charge-electmp
        if (nEDFelec>0) charge(iatm)=a(iatm)%index-electmp !Assume EDF information provides inner-core electrons for all atoms using ECP
        if (ioutmedchg==1) write(*,"(' Charge of atom',i5,'(',a2,')',': ',f12.6,'  Delta:',f12.6)") &
        iatm,a(iatm)%name,charge(iatm),charge(iatm)-lastcharge(iatm)
    end do
    
    !Check convergence
    varmax=maxval(abs(charge-lastcharge))
    if (varmax<crit) then
        if (itype==1) then
            write(*,"(a,f10.6)") " All atomic charges have converged to criterion of",crit
            write(*,"(' Sum of all charges:',f14.8)") sum(charge)
            !Normalize to get rid of integration inaccuracy
            totnumelec=sum(a%charge-charge)
            facnorm=nelec/totnumelec
            do iatm=1,ncenter
                charge(iatm)=a(iatm)%charge-facnorm*(a(iatm)%charge-charge(iatm))
            end do
            write(*,*)
            write(*,*) "Final atomic charges, after normalization to actual number of electrons"
            do iatm=1,ncenter
                write(*,"(' Atom',i5,'(',a2,')',': ',f12.6)") iatm,a(iatm)%name,charge(iatm)
            end do
            exit
        else
            write(*,*) "Hirshfeld-I atomic spaces converged successfully!"
            write(*,*)
            return
        end if
    else
        if (icyc==maxcyc) then
            write(*,"(/,' Convergence failed within',i4,' cycles!')") maxcyc
            exit
        end if
    end if
    
    !Update atomic radial density by means of interpolation of adjacent charge state
    do iatm=1,ncenter
        !Read radial density of lower limit state
        ichglow=floor(charge(iatm))
        radrholow=0
        c80tmp="atmrad"//sep//trim(a(iatm)%name)//statname(ichglow)//".rad"
        inquire(file=c80tmp,exist=alive)
        if (alive.eqv. .false.) then
            write(*,"(' Error: ',a,' was not prepared!')") trim(c80tmp)
            return
        end if
        open(10,file=c80tmp,status="old")
        read(10,*) nptlow
        do ipt=1,nptlow
            read(10,*) rnouse,radrholow(ipt)
        end do
        close(10)
        !Read radial density of upper limit state
        ichghigh=ceiling(charge(iatm))
        radrhohigh=0
        c80tmp="atmrad"//sep//trim(a(iatm)%name)//statname(ichghigh)//".rad"
        inquire(file=c80tmp,exist=alive)
        if (alive.eqv. .false.) then
            write(*,"(' Error: ',a,' was not prepared!')") trim(c80tmp)
            return
        end if
        open(10,file=c80tmp,status="old")
        read(10,*) npthigh
        do ipt=1,npthigh
            read(10,*) rnouse,radrhohigh(ipt)
        end do
        close(10)
        !Update current radial density
        atmradrho(iatm,:)=(charge(iatm)-ichglow)*radrhohigh(:) + (ichghigh-charge(iatm))*radrholow(:)
        atmradnpt(iatm)=max(npthigh,nptlow)
    end do
    
    lastcharge=charge
end do

if (allocated(frag1)) write(*,"(/,' Fragment charge:',f12.6)") sum(charge(frag1))
CALL CPU_TIME(time_end)
call walltime(iwalltime2)
write(*,"(' Calculation took up CPU time',f12.2,'s, wall clock time',i10,'s')") time_end-time_begin,iwalltime2-iwalltime1

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


!!------- Generate atomic radial density files at different states, used for such as Hirshfeld-I
!"atmrad" in current folder is used as working directory
!-2,-1,0,+1,+2 charge states of each element will be calculated to produce atomic .wfn file by Gaussian, predefined ground state multiplicity is used
!After that, radial density file (.rad) is generated for each state of each element
!If atomic wfn file is already existed, calculation will be skipped
!Radial distance values are the same as built-in atomic density, i.e. those in atmraddens.f90
subroutine genatmradfile
use defvar
use util
implicit real*8 (a-h,o-z)
character c80tmp*80,c200tmp*200,calclevel*80,radname*200,sep
character*2 :: statname(-3:3)=(/ "-3","-2","-1","_0","+1","+2","+3" /)
integer :: chgmulti(nelesupp,-3:3)=0 !Ground state multiplicity of each charge state of each element. If value=0, means undefined

!Define chgmulti for elements for possible states
!H,Li,Na,K,Rb,Cs
chgmulti(1,0)=2
chgmulti(1,1)=1
chgmulti(1,-1)=1
chgmulti(3,:)=chgmulti(1,:)
chgmulti(11,:)=chgmulti(1,:)
chgmulti(19,:)=chgmulti(1,:)
chgmulti(37,:)=chgmulti(1,:)
chgmulti(55,:)=chgmulti(1,:)
!He,Ne,Ar,Kr,Xe,Rn
chgmulti(2,0)=1
chgmulti(2,1)=2
chgmulti(2,-1)=2
chgmulti(10,:)=chgmulti(2,:)
chgmulti(18,:)=chgmulti(2,:)
chgmulti(36,:)=chgmulti(2,:)
chgmulti(54,:)=chgmulti(2,:)
chgmulti(86,:)=chgmulti(2,:)
!Be,Mg,Ca,Sr,Ba
chgmulti(4,0)=1
chgmulti(4,1)=2
chgmulti(4,2)=1
chgmulti(4,-1)=2
chgmulti(12,:)=chgmulti(4,:)
chgmulti(20,:)=chgmulti(4,:)
chgmulti(38,:)=chgmulti(4,:)
chgmulti(56,:)=chgmulti(4,:)
!B,Al,Ga,In,Tl
chgmulti(5,0)=2
chgmulti(5,1)=1
chgmulti(5,2)=2
chgmulti(5,-1)=3
chgmulti(5,-2)=4
chgmulti(13,:)=chgmulti(5,:)
chgmulti(31,:)=chgmulti(5,:)
chgmulti(49,:)=chgmulti(5,:)
chgmulti(81,:)=chgmulti(5,:)
!C,Si,Ge,Sn,Pb
chgmulti(6,0)=3
chgmulti(6,1)=2
chgmulti(6,2)=1
chgmulti(6,-1)=4
chgmulti(6,-2)=3
chgmulti(14,:)=chgmulti(6,:)
chgmulti(32,:)=chgmulti(6,:)
chgmulti(50,:)=chgmulti(6,:)
chgmulti(82,:)=chgmulti(6,:)
!N,P,As,Sb,Bi
chgmulti(7,0)=4
chgmulti(7,1)=3
chgmulti(7,2)=2
chgmulti(7,-1)=3
chgmulti(7,-2)=2
chgmulti(15,:)=chgmulti(7,:)
chgmulti(33,:)=chgmulti(7,:)
chgmulti(51,:)=chgmulti(7,:)
chgmulti(83,:)=chgmulti(7,:)
!O,S,Se,Te,Po
chgmulti(8,0)=3
chgmulti(8,1)=4
chgmulti(8,2)=3
chgmulti(8,-1)=2
chgmulti(8,-2)=1
chgmulti(16,:)=chgmulti(8,:)
chgmulti(34,:)=chgmulti(8,:)
chgmulti(52,:)=chgmulti(8,:)
chgmulti(84,:)=chgmulti(8,:)
!F,Cl,Br,I,At
chgmulti(9,0)=2
chgmulti(9,1)=3
chgmulti(9,2)=4
chgmulti(9,-1)=1
chgmulti(17,:)=chgmulti(9,:)
chgmulti(35,:)=chgmulti(9,:)
chgmulti(53,:)=chgmulti(9,:)
chgmulti(85,:)=chgmulti(9,:)
!Spin multiplicity of transition metal for each state is determined by chemical intuition as well as a few single point energy data
!For simplicity, I assume that later elements in each row has identical configuration, of course this is incorrect but not too bad
!Sc (3d1,4s2)
chgmulti(21,0)=2
chgmulti(21,1)=3
chgmulti(21,2)=2
chgmulti(21,-1)=3
chgmulti(39,:)=chgmulti(21,:) !Y
chgmulti(57,:)=chgmulti(21,:) !La
!Ti (3d2,4s2)
chgmulti(22,0)=3
chgmulti(22,1)=4
chgmulti(22,2)=3
chgmulti(22,-1)=4
chgmulti(40,:)=chgmulti(22,:) !Zr
chgmulti(72,:)=chgmulti(22,:) !Hf
!V  (3d3,4s2)
chgmulti(23,0)=4
chgmulti(23,1)=5
chgmulti(23,2)=4
chgmulti(23,-1)=5
chgmulti(41,:)=chgmulti(23,:) !Nb
chgmulti(73,:)=chgmulti(23,:) !Ta
!Cr (3d5,4s1)
chgmulti(24,0)=7
chgmulti(24,1)=6
chgmulti(24,2)=5
chgmulti(24,-1)=6
chgmulti(42,:)=chgmulti(24,:) !Mo
chgmulti(74,:)=chgmulti(24,:) !W
!Mn (3d5,4s2)
chgmulti(25,0)=6
chgmulti(25,1)=7
chgmulti(25,2)=6
chgmulti(25,-1)=5
chgmulti(43,:)=chgmulti(25,:) !Tc
chgmulti(75,:)=chgmulti(25,:) !Re
!Fe (3d6,4s2)
chgmulti(26,0)=5
chgmulti(26,1)=6
chgmulti(26,2)=7
chgmulti(26,-1)=4
chgmulti(44,:)=chgmulti(26,:) !Ru
chgmulti(76,:)=chgmulti(26,:) !Os
!Co (3d7,4s2)
chgmulti(27,0)=4
chgmulti(27,1)=5
chgmulti(27,2)=4
chgmulti(27,-1)=3
chgmulti(45,:)=chgmulti(27,:) !Rh
chgmulti(77,:)=chgmulti(27,:) !Ir
!Ni (3d8,4s2)
chgmulti(28,0)=3
chgmulti(28,1)=4
chgmulti(28,2)=3
chgmulti(28,-1)=2
chgmulti(46,:)=chgmulti(28,:) !Pd
chgmulti(78,:)=chgmulti(28,:) !Pt
!Cu (3d10,4s1)
chgmulti(29,0)=2
chgmulti(29,1)=1
chgmulti(29,2)=2
chgmulti(29,-1)=1
chgmulti(47,:)=chgmulti(29,:) !Ag
chgmulti(79,:)=chgmulti(29,:) !Au
!Zn (3d10,4s2)
chgmulti(30,0)=1
chgmulti(30,1)=2
chgmulti(30,2)=1
chgmulti(30,-1)=2
chgmulti(48,:)=chgmulti(30,:) !Cd
chgmulti(80,:)=chgmulti(30,:) !Hg

sep='/' !Separation symbol of directory
if (isys==1) sep='\'
calclevel=" "
! calclevel="B3LYP/def2SVP"

!Cycle each charge state of each atom. Each element is only calculated once. If the file is existing, don't recalculate again
do iatm=1,ncenter
    iele=a(iatm)%index
    do istat=-3,3
        if (chgmulti(iele,istat)==0) cycle !Undefined state
        radname="atmrad"//sep//trim(a(iatm)%name)//statname(istat)//".wfn"
        inquire(file=radname,exist=alive)
        if (alive) cycle
        
        !Check Gaussian path
        inquire(file=gaupath,exist=alive)
        if (.not.alive) then
            write(*,*) "Couldn't find Gaussian path defined in ""gaupath"" variable in settings.ini"
            if (isys==1) write(*,*) "Input the path of Gaussian executable file, e.g. ""d:\study\g09w\g09.exe"""
            if (isys==2.or.isys==3) write(*,*) "Input the path of Gaussian executable file, e.g. ""/sob/g09/g09"""
            do while(.true.)
                read(*,"(a)") gaupath
                inquire(file=gaupath,exist=alive)
                if (alive) exit
                write(*,*) "Couldn't find Gaussian executable file, input again"
            end do
        end if
        
        !Input calculation level
        if (calclevel==" ") then
            write(*,*) "Some atomic .wfn files are not found in ""atmrad"" folder in current directory"
            write(*,"(a)") " Now input the level for calculating these .wfn files, e.g. B3LYP/def2SVP"
            write(*,"(a)") " You can also add other keywords at the same time, e.g. M062X/6-311G(2df,2p) scf=xqc int=ultrafine"
            read(*,"(a)") calclevel
        end if
        
        !Generate .gjf file 
        inquire(file="./atmrad/.",exist=alive)
        if (alive.eqv. .false.) call system("mkdir atmrad")
        c200tmp="atmrad"//sep//trim(a(iatm)%name)//statname(istat)//".gjf"
        open(10,file=c200tmp,status="replace")
        write(10,"(a)") "# "//trim(calclevel)//" out=wfn"
        write(10,*)
        write(10,"(a)") trim(a(iatm)%name)//statname(istat)
        write(10,*)
        write(10,"(2i3)") istat,chgmulti(iele,istat)
        write(10,"(a)") a(iatm)%name
        write(10,*)
        c200tmp="atmrad"//sep//trim(a(iatm)%name)//statname(istat)//".wfn"
        write(10,"(a)") trim(c200tmp)
        write(10,*)
        write(10,*)
        close(10)
        
        !Start calculation
        c80tmp="atmrad"//sep//trim(a(iatm)%name)//statname(istat)
        write(*,*) "Running: "//trim(Gaupath)//' "'//trim(c80tmp)//'.gjf" "'//trim(c80tmp)//'"'
        call system(trim(Gaupath)//' "'//trim(c80tmp)//'.gjf" "'//trim(c80tmp)//'"')
        
        !Check if Gaussian task was successfully finished
        if (isys==1) then
            inquire(file=trim(c80tmp)//".out",exist=alive)
        else
            inquire(file=trim(c80tmp)//".log",exist=alive)
        end if
        if (alive) then
            if (isys==1) then
                open(10,file=trim(c80tmp)//".out",status="old")
            else
                open(10,file=trim(c80tmp)//".log",status="old")
            end if
            call loclabel(10,"Normal termination",igaunormal)
            close(10)
            if (igaunormal==0) then
                write(*,"(a)") " Gaussian running may be failed! Please manually check Gaussian input and output files in atmrad folder"
                write(*,*) "Press ENTER to continue"
                read(*,*)
            end if
        else
            write(*,"(a)") " Gaussian running may be failed! Please manually check Gaussian input and output files in atmrad folder"
            write(*,*) "Press ENTER to continue"
            read(*,*)
        end if
    end do
end do

!All element wfn files have been generated, now calculate corresponding radial density file (.rad)
!Existing .rad file will not be re-calculated
write(*,*)
write(*,*) "Generating atomic radial density from atomic wfn file..."
do iatm=1,ncenter
    iele=a_org(iatm)%index
    do istat=-3,3
        if (chgmulti(iele,istat)==0) cycle !Undefined state
        c80tmp="atmrad"//sep//trim(a_org(iatm)%name)//statname(istat)
        inquire(file=trim(c80tmp)//".rad",exist=alive)
        if (alive) cycle
        inquire(file=trim(c80tmp)//".wfn",exist=alive)
        if (alive.eqv. .false.) then
            write(*,"(' Error: ',a,' was not found!')") trim(c80tmp)//".wfn"
            write(*,*) "If you want to skip, press ENTER directly"
            read(*,*)
            cycle
        end if
        write(*,"(' Converting ',a,' to ',a)") trim(c80tmp)//".wfn",trim(c80tmp)//".rad"
        call atmwfn2atmrad(trim(c80tmp)//".wfn",trim(c80tmp)//".rad")
    end do
end do

!Recover to the firstly loaded file
call dealloall
call readinfile(firstfilename,1)
end subroutine


!!----- Generate atomic radial density from atomic wfn file
!The code is adapted from sphatmraddens
subroutine atmwfn2atmrad(infile,outfile)
use defvar
use function
implicit real*8 (a-h,o-z)
character(len=*) infile,outfile
real*8,allocatable :: potx(:),poty(:),potz(:),potw(:),radpos(:),sphavgval(:)
call dealloall
call readinfile(infile,1)
truncrho=1D-8
rlow=0D0
rhigh=12
nsphpt=974
nradpt=200 !Totally 200 radial points, but the number of point is truncated at truncrho (because the interpolation routine doesn't work well for very low value)
allocate(potx(nsphpt),poty(nsphpt),potz(nsphpt),potw(nsphpt),radpos(nradpt),sphavgval(nradpt))
sphavgval=0
call Lebedevgen(nsphpt,potx,poty,potz,potw)
nthreads=getNThreads()
!$OMP PARALLEL DO SHARED(sphavgval,radpos) PRIVATE(irad,radx,radr,isph,rnowx,rnowy,rnowz) schedule(dynamic) NUM_THREADS(nthreads)
do irad=1,nradpt
    radx=cos(irad*pi/(nradpt+1))
    radr=(1+radx)/(1-radx) !Becke transform
    radpos(irad)=radr
    do isph=1,nsphpt
        rnowx=potx(isph)*radr
        rnowy=poty(isph)*radr
        rnowz=potz(isph)*radr
        sphavgval(irad)=sphavgval(irad)+fdens(rnowx,rnowy,rnowz)*potw(isph)
    end do
end do
!$OMP END PARALLEL DO
open(10,file=outfile,status="replace")
write(10,*) count(sphavgval>truncrho)
do irad=nradpt,1,-1
    if (sphavgval(irad)>truncrho) write(10,"(f20.12,E18.10)") radpos(irad),sphavgval(irad)
end do
close(10)
end subroutine


!!---- Calculate density at a point for iatm based on loaded atomic radial density
real*8 function fdens_rad(iatm,x,y,z)
use defvar
use util
integer iatm,npt
real*8 x,y,z,r,rnouse
npt=atmradnpt(iatm)
r=dsqrt((a(iatm)%x-x)**2+(a(iatm)%y-y)**2+(a(iatm)%z-z)**2)
call lagintpol(atmradpos(1:npt),atmradrho(iatm,1:npt),npt,r,fdens_rad,rnouse,rnouse,1)
end function
