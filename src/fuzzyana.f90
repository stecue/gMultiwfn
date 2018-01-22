!------ Integrate fuzzy atomic space
!Normally iwork=0. If iwork=1, directly choose isel==4 to calculate delocalization index in fuzzy atomic spase (namely fuzzy bond order, see statement in JPCA,110,5108 below Eq.9) and then return
!If iwork=2, directly choose isel==8 to calculate Laplacian bond order and then return
!
!The integration grid is directly controlled by sphpot and radpot in settings.ini, since integrand may be not proportional to electron density,
!the grid will not be adjusted automatically as proposed by Becke for more efficient integration of XC functional
subroutine intatomspace(iwork)
use function
use util
use topo
implicit real*8 (a-h,o-z)
integer iwork
integer atmcalclist(ncenter),natmcalclist !The atoms to be calculated will be recorded in the array
real*8 potx(sphpot),poty(sphpot),potz(sphpot),potw(sphpot)
type(content) gridatm(radpot*sphpot)
real*8 smat(ncenter,ncenter),Pvec(ncenter),rintval(ncenter,10),funcval(radpot*sphpot),atmspcweight(radpot*sphpot) !rintval store integral, can record 10 integrand at the same time
real*8 rintvalp(ncenter,10) !Private for each OpenMP thread
real*8 promol(radpot*sphpot),atomdens(radpot*sphpot),selfdens(radpot*sphpot),selfdensgrad(3,radpot*sphpot) !For Hirshfeld partition. selfdensgrad2 is only used in Shubin's project
real*8 specrho(radpot*sphpot),specrhograd2(radpot*sphpot) !Density and its gradient^2 of atom in specific state (user-provided atomic wavefunction). Used for taking Hirshfeld as reference to calculate relative Shannon and Fisher entropy
real*8 :: covr_becke(0:nelesupp)=0D0 !covalent radii used for Becke partition
real*8 DI(ncenter,ncenter),DIa(ncenter,ncenter),DIb(ncenter,ncenter) !Delocalization index matrix
real*8 LI(ncenter),LIa(ncenter),LIb(ncenter) !Localization index array
real*8 CLRK(ncenter,ncenter) !Condensed linear response kernel
real*8 ovlpinttot(ncenter,ncenter),ovlpintpos(ncenter,ncenter),ovlpintneg(ncenter,ncenter),ovlpintpostmp(ncenter,ncenter),ovlpintnegtmp(ncenter,ncenter) !Integration between fuzzy atoms, store positive part and negative part respectively
real*8 atmmono(ncenter) !Atomic monopole, filled during multipole integration task
integer :: ifunc=3,ipartition=1,PDIatom(6),FLUatom(ncenter),FLUorb(nmo),PLRatom(6)
real*8,allocatable :: AOM(:,:,:),AOMa(:,:,:),AOMb(:,:,:),AOMsum(:,:),AOMsuma(:,:),AOMsumb(:,:) !AOM(i,j,k) means overlap matrix of MO i,j in atom k space
real*8 :: FLUref(nelesupp,nelesupp)=-1D0
real*8 :: AOMtmp(nmo,nmo),orbval(nmo)
integer :: iraddefine=-1 !-1= Specific for Laplacian bond order. 0=Custom 1=CSD 2=Pyykko 3=Suresh
integer :: nbeckeiter=3,sphpotold,radpotold
integer :: cenind(10) !Record atom index for multicenter DI
real*8 hess(3,3),rhogradw(3)
character :: radfilename*200,selectyn,c80inp*80,specatmfilename*80,c200tmp*200
real*8,external :: fdens_rad
if (ispecial==2) then
    ipartition=2 !Use Hirshfeld for shubin's 2nd project
    expcutoff=1 !Full accuracy
end if

do i=1,ncenter !Initialize the list of the atoms to be integrated
    atmcalclist(i)=i
end do
natmcalclist=ncenter

if (all(covr_becke==0D0)) then
    if (iraddefine==-1) covr_becke=covr_tianlu !The first time
    if (iraddefine==1) covr_becke=covr
    if (iraddefine==2) covr_becke=covr_pyy
    if (iraddefine==3) covr_becke=covr_Suresh
end if
if (all(FLUref==-1D0)) then !If =-1, means missing reference value
    FLUref(6,6)=1.468D0 !Calculated for benzene under HF/6-31G* opted
    FLUref(6,7)=1.566D0 !Pyridine
    FLUref(7,6)=1.566D0
    FLUref(5,7)=1.260D0 !Borazine
    FLUref(7,5)=1.260D0
end if
!Backup original grid setting, because calculating bond order may use lower grid quality
radpotold=radpot
sphpotold=sphpot


!==== Interface loop ====!
!==== Interface loop ====!
do while(.true.) 

!For some functions, e.g. calculate DI, it is safe to use relatively low grid quality for saving time,
!so sphpot and radpot may be adjusted automatically, but each time enter main interface we recover the one set by users
radpot=radpotold
sphpot=sphpotold
if (iwork==0) then
    write(*,*) "                ======== Fuzzy atomic space analysis ========"
    if (numcp>0) write(*,*) "-11 Choose a critical point as reference point"
    write(*,"(a,3f10.6,' Bohr')") " -10 Set X,Y,Z of reference point, current: ",refx,refy,refz
    if (natmcalclist==ncenter) write(*,*) "-5 Define the atoms to be calculated in functions 1 and 2, current: all atoms"
    if (natmcalclist/=ncenter) write(*,"(a,i5,a)") " -5 Define the atoms to be calculated in functions 1 and 2, current:",natmcalclist," atoms"
    write(*,*) "-4 Adjust reference parameter for FLU"
    if (ipartition==1) then !For Becke
        write(*,"(' -3 Set the number of iterations for Becke partition, current:',i3)") nbeckeiter
        if (iraddefine==-1) write(*,*) "-2 Select radius definition for Becke partition, current: Modified CSD"
        if (iraddefine==0) write(*,*) "-2 Select radius definition for Becke partition, current: Custom"
        if (iraddefine==1) write(*,*) "-2 Select radius definition for Becke partition, current: CSD"
        if (iraddefine==2) write(*,*) "-2 Select radius definition for Becke partition, current: Pyykko"
        if (iraddefine==3) write(*,*) "-2 Select radius definition for Becke partition, current: Suresh"
        if (iraddefine==4) write(*,*) "-2 Select radius definition for Becke partition, current: Hugo"
    end if
    if (ipartition==1) write(*,"(a)") " -1 Select method for partitioning atomic space, current: Becke"
    if (ipartition==2) write(*,"(a)") " -1 Select method for partitioning atomic space, current: Hirshfeld"
    if (ipartition==3) write(*,"(a)") " -1 Select method for partitioning atomic space, current: Hirshfeld*"
    if (ipartition==4) write(*,"(a)") " -1 Select method for partitioning atomic space, current: Hirshfeld-I"
    write(*,*) "0 Return"
    write(*,*) "1 Perform integration in fuzzy atomic spaces for a real space function"
    write(*,*) "2 Calculate atomic multipole moments"
    write(*,*) "3 Calculate and output atomic overlap matrix to AOM.txt in current folder"
    write(*,*) "4 Calculate localization and delocalization index (Fuzzy bond order)"
    write(*,*) "5 Calculate PDI (Para-delocalization index)"
    write(*,*) "6 Calculate FLU (Aromatic fluctuation index)"
    write(*,*) "7 Calculate FLU-pi"
    write(*,*) "8 Perform integration in fuzzy overlap region for a real space functions"
    if (allocated(CObasa)) write(*,*) "9 Calculate condensed linear response kernel (CLRK)" !Need virtual orbital informations
    if (allocated(CObasa)) write(*,*) "10 Calculate PLR (Para linear response index)" !Need virtual orbital informations
    write(*,*) "11 Calculate multi-center delocalization index" !Only can be used for HF/DFT close-shell wavefunction
!     write(*,*) "101 Integrating real space function in Hirshfeld atom with molecular grid"
    if (ispecial==2) then
        write(*,*) "99 Calculate relative Shannon and Fisher entropy and 2nd-order term"
        write(*,"(a)") " 100 Calculate relative Shannon and Fisher entropy of specific state w.r.t. Hirshfeld density"
        write(*,*) "102 Obtain quadratic and cubic Renyi entropy"
        write(*,*) "103 Obtain quadratic and cubic Renyi relative entropy"
    end if
    read(*,*) isel
    
else if (iwork==1) then
    isel=4 !Directly calculate delocalization index
else if (iwork==2) then
    isel=8 !Directly calculate Laplacian bond order
end if


!!===================================
!!--------- Adjust settings ---------
!!===================================
if (isel==0) then
    exit
    
else if (isel==101) then
    call intHirsh_molgrid
    
else if (isel==-11) then
    if (numcp>0) then
        write(*,*) "Summary of found CPs:"
        write(*,*) " Index              Coordinate               Type"
        do icp=1,numcp
            write(*,"(i6,3f12.6,3x,a)") icp,CPpos(:,icp),CPtyp2lab(CPtype(icp))
        end do
        write(*,*) "Select a CP by inputting its index, e.g. 5"
        read(*,*) icp
        refx=CPpos(1,icp)
        refy=CPpos(2,icp)
        refz=CPpos(3,icp)
    else
        write(*,*) "Error: No CPs have been found"
    end if
    
else if (isel==-10) then
    write(*,*) "Input X,Y,Z of reference point, e.g. 3.0,-4.12,0.0"
    read(*,*) refx,refy,refz
    write(*,*) "You inputted coordinate is in which unit?  1:Bohr  2:Angstrom"
    read(*,*) iunit
    if (iunit==2) then
        refx=refx/b2a
        refy=refy/b2a
        refz=refz/b2a
    end if
    
else if (isel==-5) then
    do while(.true.)
        write(*,*) "Input atom indices, e.g. 1,3-7,9,12"
        read(*,"(a)") c200tmp
        call str2arr(c200tmp,natmcalclist,atmcalclist)
        if (any(atmcalclist(1:natmcalclist)>ncenter).or.any(atmcalclist(1:natmcalclist)<=0)) then
            write(*,*) "One or more atoms exceeded valid range!"
        else
            exit
        end if
    end do
    write(*,*) "Done! The atoms you chose:"
    write(*,"(10i6)") atmcalclist(1:natmcalclist)
    write(*,*)
    
else if (isel==-4) then
    write(*,*) "Current FLU reference paramters:"
    do iref=1,nelesupp
        do jref=iref,nelesupp
            if (FLUref(iref,jref)/=-1) write(*,"(' ',a,a,a,a,f10.5)") ind2name(iref),'-',ind2name(jref),':',FLUref(iref,jref)
        end do
    end do
    do while(.true.)
        write(*,*) "Input two element indices and a new reference parameter"
        write(*,*) "e.g. 6,7,1.35  means set reference parameter for C-N to 1.35"
        write(*,*) "(Input q can return)"
        read(*,"(a)") c80inp
        if (c80inp(1:1)=='q'.or.c80inp(1:1)=='Q') exit
        read(c80inp,*) itmp,jtmp,refval
        FLUref(itmp,jtmp)=refval
        FLUref(jtmp,itmp)=refval
        write(*,*) "Done!"
    end do
    
else if (isel==-3) then
    nbeckeiterold=nbeckeiter
    write(*,*) "Do how many times of iteration? e.g. 3"
    write(*,*) "Note: Larger value gives rise to sharper atomic boundary"
    read(*,*) nbeckeiter
    if (nbeckeiter/=nbeckeiterold) then
        if (allocated(AOM)) deallocate(AOM,AOMsum)
        if (allocated(AOMa)) deallocate(AOMa,AOMb,AOMsuma,AOMsumb)
    end if
    
else if (isel==-2) then
    if (allocated(AOM)) deallocate(AOM,AOMsum)
    if (allocated(AOMa)) deallocate(AOMa,AOMb,AOMsuma,AOMsumb)
    do while(.true.)
        write(*,*)
        write(*,*) "-1 Use the modified version of CSD radii defined by Tian Lu"
        write(*,*) "0 Return"
        write(*,*) "1 Use CSD radii (Dalton Trans., 2008, 2832-2838)"
        write(*,*) "2 Use Pyykko radii (Chem. Eur.-J., 15, 186-197)"
        write(*,*) "3 Use Suresh radii (J. Phys. Chem. A, 105, 5940-5944)"
        write(*,*) "4 Use Hugo radii (Chem. Phys. Lett., 480, 127-131)"
        write(*,*) "10 Read radii from external file"
        write(*,*) "11 Modify current radii by manual input"
        write(*,*) "12 Print current radii list"
        
        read(*,*) iselrad
        if (iselrad==-1) then
            covr_becke=covr_TianLu
            iraddefine=-1
            write(*,*) "Done!"
        else if (iselrad==0) then
            exit
        else if (iselrad==1) then
            covr_becke=covr
            iraddefine=1
            write(*,*) "Done!"
        else if (iselrad==2) then
            covr_becke=covr_pyy
            iraddefine=2
            write(*,*) "Done!"
        else if (iselrad==3) then
            covr_becke=covr_Suresh
            iraddefine=3
            write(*,*) "Done!"
        else if (iselrad==4) then
            covr_becke=radii_hugo
            iraddefine=4
            write(*,*) "Done!"
        else if (iselrad==10) then
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
        else if (iselrad==11) then
            iraddefine=0
            write(*,*) "Input element index and radius (in Angstrom), e.g. 5,0.84"
            read(*,*) indtmp,radtmp
            covr_becke(indtmp)=radtmp/b2a
            write(*,*) "Done!"
        else if (iselrad==12) then
            do irad=0,nelesupp
                write(*,"(' Element:',i5,'(',a,')   Radius:',f8.3,' Angstrom')") irad,ind2name(irad),covr_becke(irad)*b2a
            end do
        end if
    end do
    
else if (isel==-1) then
    ipartitionold=ipartition
    write(*,*) "Select atomic space partition method"
    write(*,*) "1 Becke"
    write(*,*) "2 Hirshfeld"
    write(*,*) "3 Hirshfeld* (preferred over 2)"
    write(*,*) "4 Hirshfeld-I"
    write(*,"(a)") " Note: (2) uses atomic .wfn files to calculate Hirshfeld weights, they must be provided by yourself or let Multiwfn automatically &
    invoke Gaussian to generate them. (3) evaluates the weights based on built-in radial atomic densities, thus is more convenient than (2)"
    read(*,*) ipartition
    if (imodwfn==1.and.(ipartition==2.or.ipartition==4)) then !These two modes need reloading firstly loaded file, so they cannot be already modified
        write(*,"(a)") " Error: Since the wavefunction has been modified by you or by other functions, present function is unable to use. &
        Please reboot Multiwfn and reload the file"
        ipartition=ipartitionold
        cycle
    end if
    if (ipartition/=ipartitionold) then
        if (allocated(AOM)) deallocate(AOM,AOMsum)
        if (allocated(AOMa)) deallocate(AOMa,AOMb,AOMsuma,AOMsumb)
    end if
    if (ipartition==4) then !Generate radial density of all atoms by Hirshfeld-I
        call Hirshfeld_I(2)
    end if
end if
if (isel==101.or.isel<0) cycle


!!=======================================
!!--------- Prepare calculation ---------
!!=======================================

if (isel==1.or.isel==8) then !Select which function to be integrated in single atomic space or overlap between two atomic spaces
    if (isel==8.and.ipartition/=1) then
        write(*,"(a)") " Error: Only the fuzzy atomic space defined by Becke can be used together with this function"
        cycle
    end if
    if (iwork==2) then  !When ==2, means calculate Laplacian bond order
        ifunc=3 !Laplacian of rho
        if (iautointgrid==1) then !Allow change integration grid adapatively. Do not use sphpot=230/266, the result will be frantic, I don't know why
            radpot=45
            sphpot=302
        end if
    else
        write(*,*) "-2 Deformation density"
        call selfunc_interface(ifunc)
    end if
else if (isel==2) then !Multipole moment integral need electron density
    ifunc=1
    write(*,*) "Note: All units below are in a.u."
    write(*,*)
    
!AOM,LI/DI,PDI,FLU/-pi/CLRK/PLR/Multicenter DI. Note: MO values will be generated when collecting data
else if (isel==3.or.isel==4.or.isel==5.or.isel==6.or.isel==7.or.isel==9.or.isel==10.or.isel==11) then
    !Even (30,110) can be used for fuzzy bond order, so (45,170) is absolutely enough
    if (iautointgrid==1) then !Allow change integration grid adapatively
        radpot=45
        sphpot=170
    end if
    !FLU,FLU-pi,PDI,PLR,CLRK are only applied to close-shell system
    if (isel==5.or.isel==6.or.isel==7.or.isel==9.or.isel==10) then
        if (wfntype/=0.and.wfntype/=3) then
            write(*,*) "Error: This function is only available for close-shell system!"
            cycle
        end if
    else if (isel==11.and.wfntype/=0) then
        write(*,"(a)") " Error: This function is only available for single-determinant close-shell system!"
        cycle
    end if
    !Allocate space for AOM
    if ((isel==9.or.isel==10).and.allocated(AOM).and.size(AOM,1)==nmo) then !The AOM calculated for FLU/PDI is smaller than nmo, since virtual orbitals are not taken into account
        goto 10
    else if ((isel==5.or.isel==6.or.isel==7).and.(allocated(AOM).or.allocated(AOMa))) then !Don't calculate AOM again. If AOM is calculated for PLR, the AOM is also applicable since it is more than necessary
        write(*,*) "Note: AOM has already been generated before, so skipping its calculation"
        goto 10
    else !Haven't calculate AOM, hence it is needed to be calculated this time
        if (wfntype==0.or.wfntype==2.or.wfntype==3) then !RHF,ROHF,R-post-HF
            if (wfntype==3.or.isel==9.or.isel==10) then !R-post-HF or CLRK or PLR, need to consider all orbitals
                nmatsize=nmo
            else !RHF,ROHF
                !High-lying virtual orbitals will be deleted, especially for .fch case
                !Notice that occupation number may be not contiguous, some low-lying orbital may have
                !zero occupation due to modification by users, so we can't simply use nelec to determine matrix size
                do nmatsize=nmo,1,-1
                    if (MOocc(nmatsize)/=0) exit
                end do
                if (nmo-nmatsize>0) write(*,"(' Note: The highest',i6,' virtual orbitals will not be taken into account')") nmo-nmatsize
            end if
            if (allocated(AOM)) deallocate(AOM,AOMsum) !For PLR, the previous AOM and AOMsum allocated by PDI/FLU is too small, so here should be released
            allocate(AOM(nmatsize,nmatsize,ncenter),AOMsum(nmatsize,nmatsize))
            AOM=0
            AOMsum=0
        else if (wfntype==1.or.wfntype==4) then !UHF,U-post-HF
            do iendalpha=nmo,1,-1
                if (MOtype(iendalpha)==1) exit
            end do
            if (wfntype==4) then !U-post-HF
                nmatsizea=iendalpha !Total number of alpha orbitals
                nmatsizeb=nmo-nmatsizea !Total number of beta orbitals
            else !UHF
                do nmatsizea=iendalpha,1,-1
                    if (MOocc(nmatsizea)/=0D0) exit
                end do
                if (nint(nbelec)==0) then
                    nmatsizeb=0
                else
                    do nmatsizeb=nmo,iendalpha+1,-1
                        if (MOocc(nmatsizeb)/=0D0) exit
                    end do
                    nmatsizeb=nmatsizeb-iendalpha
                end if
                if (iendalpha-nmatsizea>0) write(*,"(' Note: The highest',i6,' alpha virtual orbitals will not be taken into account')") iendalpha-nmatsizea
                if (nmo-iendalpha-nmatsizeb>0) write(*,"(' Note: The highest',i6,' beta virtual orbitals will not be taken into account')") nmo-iendalpha-nmatsizeb
            end if
            if (allocated(AOMa)) deallocate(AOMa,AOMb,AOMsuma,AOMsumb)
            allocate( AOMa(nmatsizea,nmatsizea,ncenter),AOMb(nmatsizeb,nmatsizeb,ncenter) )
            allocate( AOMsuma(nmatsizea,nmatsizea),AOMsumb(nmatsizeb,nmatsizeb) )
            AOMa=0
            AOMb=0
            AOMsuma=0
            AOMsumb=0
        end if
    end if
end if


!!=======================================
!!---------------------------------------
!!--------- Start calculation -----------
!!---------------------------------------
!!=======================================
rintval=0D0 !Clean accumulated variables
xintacc=0D0
yintacc=0D0
zintacc=0D0
ovlpintpos=0D0
ovlpintneg=0D0
if (ipartition==2.or.ifunc==-2) call setpromol !In this routine reload first molecule at the end
write(*,"(' Radial points:',i5,'    Angular points:',i5,'   Total:',i10,' per center')") radpot,sphpot,radpot*sphpot
write(*,*) "Please wait..."
write(*,*)
call walltime(nwalltime1)

call Lebedevgen(sphpot,potx,poty,potz,potw)

do iatm=1,ncenter !! Cycle each atom
    !Show progress for integrating function. For electric multipole moment integration the process is not shown, because it print data in the process
    if (isel/=2) write(*,"(' Progress:',i6,'   /',i6)") iatm,ncenter

    if ( (isel==1.or.isel==2).and.all(atmcalclist(1:natmcalclist)/=iatm) ) cycle

    !Prepare grid points on current center
    iradcut=0 !Before where the radial points will be cut
    parm=1D0
    do i=1,radpot !Combine spherical point&weights with second kind Gauss-Chebyshev method for radial part
        radx=cos(i*pi/(radpot+1))
        radr=(1+radx)/(1-radx)*parm !Becke transform
        radw=2*pi/(radpot+1)*parm**3 *(1+radx)**2.5D0/(1-radx)**3.5D0 *4*pi
        gridatm( (i-1)*sphpot+1:i*sphpot )%x=radr*potx
        gridatm( (i-1)*sphpot+1:i*sphpot )%y=radr*poty
        gridatm( (i-1)*sphpot+1:i*sphpot )%z=radr*potz
        gridatm( (i-1)*sphpot+1:i*sphpot )%value=radw*potw
        if (radcut/=0D0.and.iradcut==0.and.radr<radcut) iradcut=i-1
    end do
!     if (radcut/=0) write(*,"(i8,' points per center were discarded')") iradcut*sphpot
    gridatm%x=gridatm%x+a(iatm)%x !Move quadrature point to actual position in molecule
    gridatm%y=gridatm%y+a(iatm)%y
    gridatm%z=gridatm%z+a(iatm)%z
    
    !For integrating real space function (1,8), calculate selected function value at each point here
    !For multipole moment integration (2), calculate electron density here
    !For AOM, LI/DI, FLU/PDI calculation, wavefunction value will be computed in integration stage rather than here
    if (isel==1.or.isel==2.or.isel==8.or.isel==102) then
nthreads=getNThreads()
!$OMP parallel do shared(funcval) private(rnowx,rnowy,rnowz,i) num_threads(nthreads)
        do i=1+iradcut*sphpot,radpot*sphpot
            if (ifunc==-2.or.isel==102) then
                funcval(i)=fdens(gridatm(i)%x,gridatm(i)%y,gridatm(i)%z)
            else
                funcval(i)=calcfuncall(ifunc,gridatm(i)%x,gridatm(i)%y,gridatm(i)%z)
            end if
        end do
!$OMP end parallel do
        
        !Calculate deformation density. We've calculated total density, thus now minusing it by each atomic density in free-state
        if ((isel==1.or.isel==8).and.ifunc==-2) then
            do jatm=1,ncenter_org !Calc free atomic density
                call dealloall
                call readwfn(custommapname(jatm),1)
nthreads=getNThreads()
!$OMP parallel do shared(atomdens) private(ipt) num_threads(nthreads)
                do ipt=1+iradcut*sphpot,radpot*sphpot
                    atomdens(ipt)=fdens(gridatm(ipt)%x,gridatm(ipt)%y,gridatm(ipt)%z)
                end do
!$OMP end parallel do
                funcval=funcval-atomdens
            end do
            call dealloall
            call readinfile(firstfilename,1) !Retrieve to first loaded file(whole molecule) to calc real rho again
        end if
    end if
    
    !Calculate "iatm" atomic space weight at all points around it (recorded in atmspcweight), which will be used later
    !Also integrate fuzzy overlap region here (only available for Becke partition)
    if (ipartition==1) then !Becke
!$OMP parallel shared(atmspcweight,ovlpintpos,ovlpintneg) private(i,rnowx,rnowy,rnowz,smat,&
!$OMP ii,ri,jj,rj,rmiu,chi,uij,aij,tmps,iter,Pvec,tmpval,tmpval2,ovlpintpostmp,ovlpintnegtmp) num_threads(nthreads)
        ovlpintpostmp=0D0
        ovlpintnegtmp=0D0
!$OMP do schedule(dynamic)
        do i=1+iradcut*sphpot,radpot*sphpot
            rnowx=gridatm(i)%x
            rnowy=gridatm(i)%y
            rnowz=gridatm(i)%z
            !Calculate weight, by using Eq. 11,21,13,22 in Becke's paper (JCP 88,15)
            smat=1.0D0
            do ii=1,ncenter
                ri=dsqrt( (rnowx-a(ii)%x)**2+(rnowy-a(ii)%y)**2+(rnowz-a(ii)%z)**2 )
                do jj=1,ncenter
                    if (ii==jj) cycle
                    rj=dsqrt( (rnowx-a(jj)%x)**2+(rnowy-a(jj)%y)**2+(rnowz-a(jj)%z)**2 )
                    rmiu=(ri-rj)/distmat(ii,jj)
                    !Adjust for heteronuclear
                    chi=covr_becke(a(ii)%index)/covr_becke(a(jj)%index)
                    uij=(chi-1)/(chi+1)
                    aij=uij/(uij**2-1)
                    if (aij>0.5D0) aij=0.5D0
                    if (aij<-0.5D0) aij=-0.5D0
                    rmiu=rmiu+aij*(1-rmiu**2)
                    
                    tmps=rmiu
                    do iter=1,nbeckeiter
                        tmps=1.5D0*(tmps)-0.5D0*(tmps)**3
                    end do
                    smat(ii,jj)=0.5D0*(1-tmps)
                end do
            end do
            Pvec=1.0D0
            do ii=1,ncenter
                Pvec=Pvec*smat(:,ii)
            end do
            Pvec=Pvec/sum(Pvec)
            atmspcweight(i)=Pvec(iatm) !Normalized Pvec, Pvec contain partition weight of each atom in current point, namely i
            
            if (isel==8) then !Integration between two fuzzy atoms
                tmpval=Pvec(iatm)*funcval(i)*gridatm(i)%value
                do ii=1,ncenter !Note, ovlpint is lower triangular matrix, will be convert to full matrix during statistic stage
                    do jj=ii,ncenter
                        tmpval2=Pvec(jj)*Pvec(ii)*tmpval
                        if (tmpval>0) then !ovlpinttot will be summed up to single value in statistic stage
                            ovlpintpostmp(jj,ii)=ovlpintpostmp(jj,ii)+tmpval2
                        else
                            ovlpintnegtmp(jj,ii)=ovlpintnegtmp(jj,ii)+tmpval2
                        end if
                    end do
                end do
            end if
        end do
!$OMP end do
!$OMP CRITICAL
            ovlpintpos=ovlpintpos+ovlpintpostmp
            ovlpintneg=ovlpintneg+ovlpintnegtmp
!$OMP end CRITICAL
!$OMP end parallel
    else if (ipartition==2) then !Hirshfeld based on atomic .wfn files
        promol=0D0
        do jatm=1,ncenter_org !Calculate free atomic density of each atom and promolecular density
            call dealloall
            call readwfn(custommapname(jatm),1)
nthreads=getNThreads()
!$OMP parallel do shared(atomdens,selfdensgrad) private(ipt) num_threads(nthreads)
            do ipt=1+iradcut*sphpot,radpot*sphpot
                atomdens(ipt)=fdens(gridatm(ipt)%x,gridatm(ipt)%y,gridatm(ipt)%z)
                if (jatm==iatm.and.(isel==99.or.isel==100)) then !SPECIAL: Calculate rho and its gradient for free atom
                    call calchessmat_dens(1,gridatm(ipt)%x,gridatm(ipt)%y,gridatm(ipt)%z,selfdens(ipt),selfdensgrad(1:3,ipt),hess)
                end if
            end do
!$OMP end parallel do
            promol=promol+atomdens
            if (jatm==iatm) selfdens=atomdens
        end do
        do i=1+iradcut*sphpot,radpot*sphpot !Get Hirshfeld weight of present atom
            if (promol(i)/=0.0D0) then
                atmspcweight(i)=selfdens(i)/promol(i)
            else
                atmspcweight(i)=0D0
            end if
        end do
        call dealloall
        call readinfile(firstfilename,1) !Retrieve the firstly loaded file(whole molecule) in order to calculate real rho later
    else if (ipartition==3.or.ipartition==4) then !Hirshfeld or Hirshfeld-I
        promol=0D0
        if (ipartition==3) then !Hirshfeld based on interpolation of built-in atomic radius density
            do jatm=1,ncenter !Calculate free atomic density of each atom and promolecular density
nthreads=getNThreads()
!$OMP parallel do shared(atomdens) private(ipt) num_threads(nthreads)
                do ipt=1+iradcut*sphpot,radpot*sphpot
                    atomdens(ipt)=calcatmdens(jatm,gridatm(ipt)%x,gridatm(ipt)%y,gridatm(ipt)%z,0)
                end do
!$OMP end parallel do
                promol=promol+atomdens
                if (jatm==iatm) selfdens=atomdens
            end do
        else !Hirshfeld-I based on refined atomic radial density
            do jatm=1,ncenter !Calculate free atomic density of each atom and promolecular density
nthreads=getNThreads()
!$OMP parallel do shared(atomdens) private(ipt) num_threads(nthreads)
                do ipt=1+iradcut*sphpot,radpot*sphpot
                    atomdens(ipt)=fdens_rad(jatm,gridatm(ipt)%x,gridatm(ipt)%y,gridatm(ipt)%z)
                end do
!$OMP end parallel do
                promol=promol+atomdens
                if (jatm==iatm) selfdens=atomdens
            end do
        end if
        do i=1+iradcut*sphpot,radpot*sphpot !Get Hirshfeld weight of present atom
            if (promol(i)/=0D0) then
                atmspcweight(i)=selfdens(i)/promol(i)
            else
                atmspcweight(i)=0D0
            end if
        end do
    end if
    
    !SPECIAL CASE: Calculate density of atom in specific-state, for Shubin's idea
    if (isel==100) then
        write(specatmfilename,"(a,i4.4,a)") "specwfn/",iatm,".wfn"
        write(*,*) "Prodessing "//trim(specatmfilename)
        call dealloall
        call readwfn(specatmfilename,1)
        a=a_org(iatm) !Set atom position to actual atom position
nthreads=getNThreads()
!$OMP parallel do shared(specrho,specrhograd2) private(ipt) num_threads(nthreads)
        do ipt=1+iradcut*sphpot,radpot*sphpot
            specrho(ipt)=fdens(gridatm(ipt)%x,gridatm(ipt)%y,gridatm(ipt)%z)
            specrhograd2(ipt)=fgrad(gridatm(ipt)%x,gridatm(ipt)%y,gridatm(ipt)%z,'t')**2
        end do
!$OMP end parallel do
        call dealloall
        call readinfile(firstfilename,1) !Retrieve the first loaded file(whole molecule) to calc real rho later
    end if
    
    !Perform integration on single center
    if (isel==1) then
        do i=1+iradcut*sphpot,radpot*sphpot
            rintval(iatm,1)=rintval(iatm,1)+atmspcweight(i)*funcval(i)*gridatm(i)%value
        end do
    else if (isel==99.or.isel==100.or.isel==103) then !SPECIAL SPECIAL SPECIAL
        !=99:  Calculate relative Shannon and Fisher entropy and 2nd-order term
        !=100: Calculate relative Shannon/Fisher entropy by taking Hirshfeld density as reference
        !=103: Calculate quadratic and cubic Renyi relative entropy
nthreads=getNThreads()
!$OMP parallel shared(rintval) private(i,rnowx,rnowy,rnowz,rhow,rhogradw,rhograd2w,rintvalp,tmpx,tmpy,tmpz) num_threads(nthreads)
        rintvalp=0D0
!$OMP do schedule(dynamic)
        do i=1+iradcut*sphpot,radpot*sphpot
            rnowx=gridatm(i)%x
            rnowy=gridatm(i)%y
            rnowz=gridatm(i)%z
            if (isel==99.or.isel==100) then
                call calchessmat_dens(1,rnowx,rnowy,rnowz,rhow,rhogradw(:),hess)
                rhow=atmspcweight(i)*rhow !rhoA at current point
                rhogradw=atmspcweight(i)*rhogradw
                rhograd2w=sum(rhogradw(:)**2) !|grad_rhoA|^2 at current point
            else
                rhow=atmspcweight(i)*fdens(rnowx,rnowy,rnowz) !rhoA at current point
            end if
            if (isel==99) then
                !Relative Shannon entropy w.r.t. free-state
                rintvalp(iatm,1)=rintvalp(iatm,1)+rhow*log(rhow/selfdens(i))*gridatm(i)%value
                !Relative Fisher information entropy w.r.t. free-state (old formula, incorrect)
                rintvalp(iatm,2)=rintvalp(iatm,2)+(rhograd2w/rhow-sum(selfdensgrad(1:3,i)**2)/selfdens(i))*gridatm(i)%value
                !Atomic Shannon
                rintvalp(iatm,3)=rintvalp(iatm,3)-rhow*log(rhow)*gridatm(i)%value
                !Atomic Fisher
                rintvalp(iatm,4)=rintvalp(iatm,4)+(rhograd2w/rhow)*gridatm(i)%value
                !1st-order term: rhoA-rho0
                rintvalp(iatm,5)=rintvalp(iatm,5)+(rhow-selfdens(i))*gridatm(i)%value
                !2nd-order term
                rintvalp(iatm,6)=rintvalp(iatm,6)+(rhow-selfdens(i))**2/rhow/2D0 *gridatm(i)%value
                !Relative Fisher information entropy w.r.t. free-state (new formula, correct)
                tmpx=rhogradw(1)/rhow-selfdensgrad(1,i)/selfdens(i)
                tmpy=rhogradw(2)/rhow-selfdensgrad(2,i)/selfdens(i)
                tmpz=rhogradw(3)/rhow-selfdensgrad(3,i)/selfdens(i)
                rintvalp(iatm,7)=rintvalp(iatm,7)+rhow*(tmpx**2+tmpy**2+tmpz**2)*gridatm(i)%value
            else if (isel==100) then
                !Relative Shannon entropy of specific atomic state with Hirshfeld density as reference
                rintvalp(iatm,1)=rintvalp(iatm,1)+specrho(i)*log(specrho(i)/rhow)*gridatm(i)%value
                !Relative Fisher information entropy of specific atomic state with Hirshfeld density as reference
                rintvalp(iatm,2)=rintvalp(iatm,2)+(specrhograd2(i)/specrho(i)-rhograd2w/rhow)*gridatm(i)%value
            else if (isel==103) then
                !Quadratic Renyi relative entropy 
                rintvalp(iatm,1)=rintvalp(iatm,1)+rhow**2/selfdens(i) *gridatm(i)%value
                !Cubic Renyi relative entropy 
                rintvalp(iatm,2)=rintvalp(iatm,2)+rhow**3/selfdens(i)**2 *gridatm(i)%value
            end if
        end do
!$OMP end do
!$OMP CRITICAL
            rintval=rintval+rintvalp
!$OMP end CRITICAL
!$OMP end parallel
    else if (isel==102) then !SPECIAL SPECIAL SPECIAL: Obtain quadratic and cubic Renyi entropy
        do i=1+iradcut*sphpot,radpot*sphpot
            rintval(iatm,1)=rintval(iatm,1)+atmspcweight(i)*funcval(i)**2*gridatm(i)%value
            rintval(iatm,2)=rintval(iatm,2)+atmspcweight(i)*funcval(i)**3*gridatm(i)%value
        end do
    else if (isel==2) then !Integrate multipole momentS
        eleint=0D0
        xint=0D0
        yint=0D0
        zint=0D0
        xxint=0D0
        yyint=0D0
        zzint=0D0
        xyint=0D0
        yzint=0D0
        xzint=0D0
        xxxint=0D0
        yyyint=0D0
        zzzint=0D0
        yzzint=0D0
        xzzint=0D0
        xxzint=0D0
        yyzint=0D0
        xxyint=0D0
        xyyint=0D0
        xyzint=0D0
        rrxint=0D0
        rryint=0D0
        rrzint=0D0
        do i=1+iradcut*sphpot,radpot*sphpot
            rx=gridatm(i)%x-a(iatm)%x
            ry=gridatm(i)%y-a(iatm)%y
            rz=gridatm(i)%z-a(iatm)%z
            tmpmul=atmspcweight(i)*funcval(i)*gridatm(i)%value
            eleint=eleint-tmpmul !monopole
            xint=xint+rx*tmpmul
            yint=yint+ry*tmpmul
            zint=zint+rz*tmpmul
            xxint=xxint+rx*rx*tmpmul
            yyint=yyint+ry*ry*tmpmul
            zzint=zzint+rz*rz*tmpmul
            xyint=xyint+rx*ry*tmpmul
            yzint=yzint+ry*rz*tmpmul
            xzint=xzint+rx*rz*tmpmul
            !Used for octopole moments
            xxxint=xxxint+rx*rx*rx*tmpmul
            yyyint=yyyint+ry*ry*ry*tmpmul
            zzzint=zzzint+rz*rz*rz*tmpmul
            yzzint=yzzint+ry*rz*rz*tmpmul
            xzzint=xzzint+rx*rz*rz*tmpmul
            xxzint=xxzint+rx*rx*rz*tmpmul
            yyzint=yyzint+ry*ry*rz*tmpmul
            xxyint=xxyint+rx*rx*ry*tmpmul
            xyyint=xyyint+rx*ry*ry*tmpmul
            xyzint=xyzint+rx*ry*rz*tmpmul
        end do
        rrint=xxint+yyint+zzint
        rrxint=xxxint+xyyint+xzzint
        rryint=xxyint+yyyint+yzzint
        rrzint=xxzint+yyzint+zzzint
        atmchgtmp=a(iatm)%charge+eleint
        write(*,"('              *****  Atomic multipole moments of',i6,'(',a2,')  *****')") iatm,a(iatm)%name
        write(*,"(' Atomic monopole moment (electron):',f12.6,'   Atomic charge:',f12.6)") eleint,atmchgtmp
        write(*,"(' Atomic dipole moments:')")
        write(*,"(' X=',f12.6,'  Y=',f12.6,'  Z=',f12.6,'  Norm=',f12.6)") -xint,-yint,-zint,dsqrt(xint**2+yint**2+zint**2)
        write(*,"(' Contribution to molecular dipole moment:')")
        contridipx=atmchgtmp*a(iatm)%x-xint
        contridipy=atmchgtmp*a(iatm)%y-yint
        contridipz=atmchgtmp*a(iatm)%z-zint
        write(*,"(' X=',f12.6,'  Y=',f12.6,'  Z=',f12.6,'  Norm=',f12.6)") contridipx,contridipy,contridipz,dsqrt(contridipx**2+contridipy**2+contridipz**2)
        write(*,"(' Atomic quadrupole moments (Cartesian form):')")
        QXX=(-3*xxint+rrint)/2
        QYY=(-3*yyint+rrint)/2
        QZZ=(-3*zzint+rrint)/2
        write(*,"(' QXX=',f12.6,'  QXY=',f12.6,'  QXZ=',f12.6)") QXX,(-3*xyint)/2,(-3*xzint)/2
        write(*,"(' QYX=',f12.6,'  QYY=',f12.6,'  QYZ=',f12.6)") (-3*xyint)/2,QYY,(-3*yzint)/2
        write(*,"(' QZX=',f12.6,'  QZY=',f12.6,'  QZZ=',f12.6)") (-3*xzint)/2,(-3*yzint)/2,QZZ
        write(*,"( ' The magnitude of atomic quadrupole moment (Cartesian form):',f12.6)") sqrt(2D0/3D0*(QXX**2+QYY**2+QZZ**2))
        R20=-(3*zzint-rrint)/2D0 !Notice that the negative sign, because electrons carry negative charge
        R2n1=-dsqrt(3D0)*yzint
        R2p1=-dsqrt(3D0)*xzint
        R2n2=-dsqrt(3D0)*xyint
        R2p2=-dsqrt(3D0)/2D0*(xxint-yyint)
        write(*,"(' Atomic quadrupole moments (Spherical harmonic form):')")
        write(*,"(' Q_2,0 =',f11.6,'   Q_2,-1=',f11.6,'   Q_2,1=',f11.6)") R20,R2n1,R2p1
        write(*,"(' Q_2,-2=',f11.6,'   Q_2,2 =',f11.6)") R2n2,R2p2
        write(*,"( ' Magnitude: |Q_2|=',f12.6)") dsqrt(R20**2+R2n1**2+R2p1**2+R2n2**2+R2p2**2)
        R30=-(5*zzzint-3*rrzint)/2D0
        R3n1=-dsqrt(3D0/8D0)*(5*yzzint-rryint)
        R3p1=-dsqrt(3D0/8D0)*(5*xzzint-rrxint)
        R3n2=-dsqrt(15D0)*xyzint
        R3p2=-dsqrt(15D0)*(xxzint-yyzint)/2D0
        R3n3=-dsqrt(5D0/8D0)*(3*xxyint-yyyint)
        R3p3=-dsqrt(5D0/8D0)*(xxxint-3*xyyint)
        write(*,"(' Atomic octopole moments (Spherical harmonic form):')")
        write(*,"(' Q_3,0 =',f11.6,'  Q_3,-1=',f11.6,'  Q_3,1 =',f11.6)") R30,R3n1,R3p1
        write(*,"(' Q_3,-2=',f11.6,'  Q_3,2 =',f11.6,'  Q_3,-3=',f11.6,'  Q_3,3 =',f11.6)") R3n2,R3p2,R3n3,R3p3
        write(*,"( ' Magnitude: |Q_3|=',f12.6)") dsqrt(R30**2+R3n1**2+R3p1**2+R3n2**2+R3p2**2+R3n3**2+R3p3**2)
        write(*,*)
        atmmono(iatm)=eleint
        xintacc=xintacc+xint
        yintacc=yintacc+yint
        zintacc=zintacc+zint
    
    !Calculate atomic overlap matrix (AOM) for all tasks that require it
    else if (isel==3.or.isel==4.or.isel==5.or.isel==6.or.isel==7.or.isel==9.or.isel==10.or.isel==11) then
        if (wfntype==0.or.wfntype==2.or.wfntype==3) then !RHF,ROHF,R-post-HF
nthreads=getNThreads()
!$OMP parallel shared(AOM) private(i,imo,jmo,AOMtmp,orbval) num_threads(nthreads)
            AOMtmp=0D0
!$OMP do schedule(dynamic)
            do i=1+iradcut*sphpot,radpot*sphpot
                call orbderv(1,1,nmatsize,gridatm(i)%x,gridatm(i)%y,gridatm(i)%z,orbval) !Calculate orbital wavefunction value of all MOs in current position and store to orbval
                do imo=1,nmatsize
                    do jmo=imo,nmatsize    
                        AOMtmp(imo,jmo)=AOMtmp(imo,jmo)+atmspcweight(i)*orbval(imo)*orbval(jmo)*gridatm(i)%value
                    end do
                end do
            end do
!$OMP end do
!$OMP CRITICAL
                AOM(:,:,iatm)=AOM(:,:,iatm)+AOMtmp(1:nmatsize,1:nmatsize)
!$OMP end CRITICAL
!$OMP end parallel
            AOM(:,:,iatm)=AOM(:,:,iatm)+transpose(AOM(:,:,iatm))
            do imo=1,nmatsize
                AOM(imo,imo,iatm)=AOM(imo,imo,iatm)/2D0
            end do
        else if (wfntype==1.or.wfntype==4) then !UHF,U-post-HF
            !Alpha part
            do i=1+iradcut*sphpot,radpot*sphpot
                call orbderv(1,1,nmatsizea,gridatm(i)%x,gridatm(i)%y,gridatm(i)%z,orbval)
                do imo=1,nmatsizea
                    do jmo=imo,nmatsizea
                        AOMa(imo,jmo,iatm)=AOMa(imo,jmo,iatm)+atmspcweight(i)*orbval(imo)*orbval(jmo)*gridatm(i)%value
                    end do
                end do
            end do
            AOMa(:,:,iatm)=AOMa(:,:,iatm)+transpose(AOMa(:,:,iatm))
            do imo=1,nmatsizea
                AOMa(imo,imo,iatm)=AOMa(imo,imo,iatm)/2D0
            end do
            !Beta part
            if (nmatsizeb>0) then !If there is no beta orbital, then nmatsizeb=0, and we will do nothing
                MOinit=iendalpha+1
                MOend=iendalpha+nmatsizeb
                do i=1+iradcut*sphpot,radpot*sphpot
                    call orbderv(1,MOinit,MOend,gridatm(i)%x,gridatm(i)%y,gridatm(i)%z,orbval)
                    do imo=MOinit,MOend
                        imotmp=imo-iendalpha !So that the index start from 1 to nbelec
                        do jmo=imo,MOend
                            jmotmp=jmo-iendalpha
                            AOMb(imotmp,jmotmp,iatm)=AOMb(imotmp,jmotmp,iatm)+atmspcweight(i)*orbval(imo)*orbval(jmo)*gridatm(i)%value
                        end do
                    end do
                end do                    
                AOMb(:,:,iatm)=AOMb(:,:,iatm)+transpose(AOMb(:,:,iatm))
                do imo=1,nmatsizeb
                    AOMb(imo,imo,iatm)=AOMb(imo,imo,iatm)/2D0
                end do
            end if
        end if
    end if
    
end do !End cycling atoms

call walltime(nwalltime2)
write(*,"(' Calculation took up',i8,' seconds wall clock time')") nwalltime2-nwalltime1



!==== Check sanity of AOM ====!
if (isel==3.or.isel==4.or.isel==5.or.isel==6.or.isel==7.or.isel==9.or.isel==10.or.isel==11) then
    if (wfntype==0.or.wfntype==2.or.wfntype==3) then !RHF,ROHF,R-post-HF
        do iatm=1,ncenter
            AOMsum=AOMsum+AOM(:,:,iatm)
        end do
        AOMerror=identmaterr(AOMsum)/ncenter
        write(*,"(' Error of AOM is',f14.8)") AOMerror
        if (AOMerror>0.001D0) write(*,"(a)") " Warning: The integration is not very accurate, in the settings.ini, &
        you need to set iautointgrid to 0 and proper set radpot and sphpot parameters"        
!         call showmatgau(AOMsum,"",1,"f14.8",7)
    else if (wfntype==1.or.wfntype==4) then !UHF,U-post-HF
        AOMerrorb=0D0
        do iatm=1,ncenter
            AOMsuma=AOMsuma+AOMa(:,:,iatm)
            if (nmatsizeb>0) AOMsumb=AOMsumb+AOMb(:,:,iatm)
        end do
        AOMerrora=identmaterr(AOMsuma)/ncenter
        if (nmatsizeb>0) AOMerrorb=identmaterr(AOMsumb)/ncenter
        write(*,"(' Error of alpha AOM is',f14.8)") AOMerrora
        if (nmatsizeb>0) write(*,"(' Error of Beta AOM is ',f14.8)") AOMerrorb
        if (AOMerrora>0.001D0.or.AOMerrorb>0.001D0)  write(*,"(a)") " Warning: The integration is not very accurate, in the settings.ini, &
        you need to set iautointgrid to 0 and proper set radpot and sphpot parameters"        
    end if
end if

!==== Generate DI, LI or condensed linear response kernel (CLRK) ====!
!DI-pi will be calculated for FLU-pi at later stage
!Multicenter DI will be calculated at later stage
10    if (isel==4.or.isel==5.or.isel==6) then !For LI/DI, PDI and FLU
    if (any(MOocc<0)) then
        where(MOocc<0) MOocc=0
        write(*,"(a)") " Note: Some occupation numbers are negative. In order to make the calculation feasible, they have been set to zero"
        write(*,*) "Press ENTER to continue"
        read(*,*)
    end if
    !RHF,R-post-HF, DI_A,B=2¡Æ[i,j]dsqrt(n_i*n_j)*S_i,j_A * S_i,j_B     where i and j are non-spin orbitals
    if (wfntype==0.or.wfntype==3) then
        DI=0D0
        do iatm=1,ncenter
            do jatm=iatm,ncenter
                do iorb=1,nmatsize
                    do jorb=1,nmatsize
                        DI(iatm,jatm)=DI(iatm,jatm)+dsqrt(MOocc(iorb)*MOocc(jorb))*AOM(iorb,jorb,iatm)*AOM(iorb,jorb,jatm)
                    end do
                end do
            end do
            LI(iatm)=DI(iatm,iatm)
        end do
        DI=2*(DI+transpose(DI))
        do iatm=1,ncenter !Diagonal terms are the sum of corresponding row or column
            DI(iatm,iatm)=0D0
            DI(iatm,iatm)=sum(DI(iatm,:))
        end do
    else if (wfntype==2) then !ROHF
        DIa=0D0
        DIb=0D0
        do nmoclose=nmatsize,1,-1
            if (MOtype(nmoclose)==0) exit
        end do
        do iatm=1,ncenter
            do jatm=iatm,ncenter
                !Alpha
                do iorb=1,nmatsize !The number of close or alpha orbitals needed to be concerned
                    occi=MOocc(iorb)
                    if (MOtype(iorb)==0) occi=occi/2D0
                    do jorb=1,nmatsize
                        occj=MOocc(jorb)
                        if (MOtype(jorb)==0) occj=occj/2D0
                        DIa(iatm,jatm)=DIa(iatm,jatm)+dsqrt(occi*occj)*AOM(iorb,jorb,iatm)*AOM(iorb,jorb,jatm)
                    end do
                end do
                !Beta
                do iorb=1,nmoclose !The number of close orbitals needed to be concerned
                    do jorb=1,nmoclose
                        DIb(iatm,jatm)=DIb(iatm,jatm)+dsqrt(MOocc(iorb)/2D0*MOocc(jorb)/2D0)*AOM(iorb,jorb,iatm)*AOM(iorb,jorb,jatm)
                    end do
                end do
            end do
            LIa(iatm)=DIa(iatm,iatm)
            LIb(iatm)=DIb(iatm,iatm)
        end do
        DIa=2*(DIa+transpose(DIa))
        DIb=2*(DIb+transpose(DIb))
        do iatm=1,ncenter !Diagonal terms are the sum of corresponding row or column
            DIa(iatm,iatm)=0D0
            DIb(iatm,iatm)=0D0
            DIa(iatm,iatm)=sum(DIa(iatm,:))
            DIb(iatm,iatm)=sum(DIb(iatm,:))
        end do
        !Combine alpha and Beta to total
        DI=DIa+DIb
        LI=LIa+LIb
    !UHF,U-post-HF   DI(A,B)=2¡Æ[i,j]dsqrt(n_i*n_j)*S_i,j_A * S_i,j_B   where i and j are spin orbitals
    else if (wfntype==1.or.wfntype==4) then
        !Alpha
        DIa=0D0
        do iatm=1,ncenter
            do jatm=iatm,ncenter
                do iorb=1,nmatsizea
                    do jorb=1,nmatsizea
                        DIa(iatm,jatm)=DIa(iatm,jatm)+dsqrt(MOocc(iorb)*MOocc(jorb))*AOMa(iorb,jorb,iatm)*AOMa(iorb,jorb,jatm)
                    end do
                end do
            end do
            LIa(iatm)=DIa(iatm,iatm)
        end do
        DIa=2*(DIa+transpose(DIa))
        !Beta
        if (nmatsizeb>0) then
            DIb=0D0
            MOinit=iendalpha+1 !Index range of beta orbitals
            MOend=iendalpha+nmatsizeb
            do iatm=1,ncenter
                do jatm=iatm,ncenter
                    do iorb=MOinit,MOend
                        iorbtmp=iorb-iendalpha
                        do jorb=MOinit,MOend
                            jorbtmp=jorb-iendalpha
                            DIb(iatm,jatm)=DIb(iatm,jatm)+dsqrt(MOocc(iorb)*MOocc(jorb))*AOMb(iorbtmp,jorbtmp,iatm)*AOMb(iorbtmp,jorbtmp,jatm)
                        end do
                    end do
                end do
                LIb(iatm)=DIb(iatm,iatm)
            end do
            DIb=2*(DIb+transpose(DIb))
        end if
        do iatm=1,ncenter !Diagonal terms are the sum of corresponding row or column
            DIa(iatm,iatm)=0D0
            DIb(iatm,iatm)=0D0
            DIa(iatm,iatm)=sum(DIa(iatm,:))
            DIb(iatm,iatm)=sum(DIb(iatm,:))
        end do
        !Combine alpha and Beta to total
        DI=DIa+DIb
        LI=LIa+LIb
    end if
    
else if (isel==9.or.isel==10) then !Calculate condensed linear response kernel, PLR also uses it
    CLRK=0D0
    do iatm=1,ncenter
        do jatm=iatm,ncenter
            do iorb=1,nmo !Occupied MOs
                if (nint(MOocc(iorb))==2D0) then
                    do jorb=idxHOMO+1,nmo !Virtual MOs
                        if (nint(MOocc(jorb))==0D0) CLRK(iatm,jatm)=CLRK(iatm,jatm)+AOM(iorb,jorb,iatm)*AOM(iorb,jorb,jatm)/(MOene(iorb)-MOene(jorb))
                    end do
                end if
            end do
        end do
    end do
    CLRK=CLRK*4D0
    CLRK=CLRK+transpose(CLRK)
    do iatm=1,ncenter
        CLRK(iatm,iatm)=CLRK(iatm,iatm)/2D0
    end do
end if




!!====================================================
!!------- Statistic results or post-processing -------
!!====================================================
write(*,*)
if (isel==1) then
    sumval=sum(rintval(:,1))
    sumabsval=sum(abs(rintval(:,1)))
    write(*,*) "  Atomic space        Value                % of sum            % of sum abs"
    if (any(rintval(:,1)>1D9).or.all(rintval(:,1)<1D-7))then
        do iatm=1,ncenter
            write(*,"(i6,'(',a2,')  ',D20.10,1x,f20.6,1x,f20.6)") iatm,a(iatm)%name,rintval(iatm,1),rintval(iatm,1)/sumval*100,rintval(iatm,1)/sumabsval*100
        end do
        write(*,"(' Summing up above values:',D20.10)") sumval
        write(*,"(' Summing up absolute value of above values:',D20.10)") sumabsval
        else
        do iatm=1,ncenter
            write(*,"(i6,'(',a2,')  ',f20.8,1x,f20.6,1x,f20.6)") iatm,a(iatm)%name,rintval(iatm,1),rintval(iatm,1)/sumval*100,rintval(iatm,1)/sumabsval*100
        end do
        write(*,"(' Summing up above values:',f20.8)") sumval
        write(*,"(' Summing up absolute value of above values:',f20.8)") sumabsval
    end if
    write(*,*)
    
else if (isel==99) then !SPECIAL: Relative Shannon and Fisher entropy and 2nd-order term
    write(*,*) "Relative Shannon entropy and relative Fisher information w.r.t. its free-state"
    write(*,*) "   Atom           Rel.Shannon       Rel.Fisher(old)   Rel.Fisher(new)"
    do iatm=1,ncenter
        write(*,"(i6,'(',a2,')  ',3f18.8)") iatm,a(iatm)%name,rintval(iatm,1),rintval(iatm,2),rintval(iatm,7)
    end do
    write(*,"(' Summing up above values:',3f18.8)") sum(rintval(:,1)),sum(rintval(:,2)),sum(rintval(:,7))
    write(*,*)
    write(*,*) "Shannon and Fisher information entropy of each atom"
    write(*,*) "   Atom             Shannon            Fisher"
    do iatm=1,ncenter
        write(*,"(i6,'(',a2,')  ',2f18.8)") iatm,a(iatm)%name,rintval(iatm,3),rintval(iatm,4)
    end do
    write(*,"(' Summing up above values:',2f22.8)") sum(rintval(:,3)),sum(rintval(:,4))
    write(*,*)
    write(*,*) "1st and 2nd-order terms of each atom"
    write(*,*) "   Atom           1st           2nd"
    do iatm=1,ncenter
        write(*,"(i6,'(',a2,')  ',2f14.8)") iatm,a(iatm)%name,rintval(iatm,5),rintval(iatm,6)
    end do
    write(*,"(' Summing up above values:',2f16.8)") sum(rintval(:,5)),sum(rintval(:,6))
    write(*,*)
else if (isel==100) then !SPECIAL: Relative Shannon/Fisher by taking Hirshfeld density as reference
    write(*,*) "Relative Shannon and Fisher entropy of specific state w.r.t. Hirshfeld density"
    write(*,*) "   Atom         Relat_Shannon      Relat_Fisher"
    do iatm=1,ncenter
        write(*,"(i6,'(',a2,')  ',2f18.8)") iatm,a(iatm)%name,rintval(iatm,1),rintval(iatm,2)
    end do
    write(*,*)
else if (isel==102) then !SPECIAL: Quadratic and cubic Renyi entropy
    write(*,*) "Atomic contribution to int(rho^2) and int(rho^3) under Hirshfeld partition:"
    write(*,*) "   Atom            Quadratic             Cubic"
    do iatm=1,ncenter
        write(*,"(i6,'(',a2,')  ',2(1PE20.8))") iatm,a(iatm)%name,rintval(iatm,1),rintval(iatm,2)
    end do
    write(*,"('    Total   ',2(1PE20.8))") sum(rintval(:,1)),sum(rintval(:,2))
    write(*,*)
    write(*,"(' Molecular quadratic Renyi entropy:',f18.8)") -log10(sum(rintval(:,1)))
    write(*,"(' Molecular cubic Renyi entropy:    ',f18.8)") -log10(sum(rintval(:,2)))/2
    write(*,*)
else if (isel==103) then !SPECIAL: Quadratic and cubic Renyi relative entropy
    write(*,"(a)") " Note: rhoA=w_A(r)*rho(r) is density of A in molecule, rhoA0 is density of A in its free-state"
    write(*,*) "   Atom        int(rhoA^2/rhoA0)   int(rhoA^3/rhoA0^2)"
    do iatm=1,ncenter
        write(*,"(i6,'(',a2,')  ',2(1PE20.8))") iatm,a(iatm)%name,rintval(iatm,1),rintval(iatm,2)
    end do
    write(*,"('    Total   ',2(1PE20.8))") sum(rintval(:,1)),sum(rintval(:,2))
    write(*,*)
    write(*,"(' Molecular quadratic Renyi relative entropy:',f18.8)") -log10(sum(rintval(:,1)))
    write(*,"(' Molecular cubic Renyi relative entropy:    ',f18.8)") -log10(sum(rintval(:,2)))
    write(*,*)
        
else if (isel==2) then !Multipole moment integration
    write(*,"(' Total number of electrons:',f14.6)") -sum(atmmono)
    xmoldip=-xintacc+sum(a(:)%x*(atmmono(:)+a(:)%charge))
    ymoldip=-yintacc+sum(a(:)%y*(atmmono(:)+a(:)%charge))
    zmoldip=-zintacc+sum(a(:)%z*(atmmono(:)+a(:)%charge))
    write(*,"(' Molecular dipole moment (a.u.): ',3f14.6)") xmoldip,ymoldip,zmoldip
    write(*,"(' Molecular dipole moment (Debye):',3f14.6)") xmoldip*au2debye,ymoldip*au2debye,zmoldip*au2debye
    xmoldipmag=sqrt(xmoldip**2+ymoldip**2+zmoldip**2)
    write(*,"(' Magnitude of molecular dipole moment (a.u.&Debye):',2f14.6)") xmoldipmag,xmoldipmag*au2debye
    write(*,*)
    
else if (isel==3) then !Output AOM
    open(10,file="AOM.txt",status="replace")
    if (wfntype==0.or.wfntype==2.or.wfntype==3) then
        do iatm=1,ncenter
            write(10,"('Atomic overlap matrix of',i6,'(',a2,')')") iatm,a(iatm)%name
            call showmatgau(AOM(:,:,iatm),"",1,"f14.8",10)
            write(10,*)
        end do
    else if (wfntype==1.or.wfntype==4) then
        do iatm=1,ncenter
            write(10,"('Alpha part of atomic overlap matrix of',i6,'(',a2,')')") iatm,a(iatm)%name
            call showmatgau(AOMa(:,:,iatm),"",1,"f14.8",10)
            if (nmatsizeb>0) then
                write(10,"('Beta part of atomic overlap matrix of',i6,'(',a2,')')") iatm,a(iatm)%name
                call showmatgau(AOMb(:,:,iatm),"",1,"f14.8",10)
            end if
            write(10,*)
        end do
    end if
    close(10)
    write(*,*) "Done, AOM have been exported to AOM.txt in current folder"
    write(*,*)
    
else if (isel==4) then !Show LI and DI or fuzzy bond order
    if (iwork==0) then !Output LI and DI
        write(*,"(a)") " Note: Delocalization index in fuzzy atomic space is also known as fuzzy bond order"
        !The strict definition of atomic valence in fuzzy space is Eq.18 in CPL,368,375, however in close-shell case free valence is zero, so sum of bond order is just atomic valence
        write(*,"(a)") " Note: Diagonal terms are the sum of corresponding row or column elements, for close-shell cases, they are also known as atomic valence"
        write(*,*)
        selectyn='n'
        ioutid=6
        do while(.true.)
            if (wfntype==1.or.wfntype==2.or.wfntype==4) then !UHF,ROHF,U-post-HF, output each spin component first
                !Alpha
                call showmatgau(DIa,"Delocalization index matrix for alpha spin",0,"f14.8",ioutid)
                write(ioutid,*)
                if (iwork/=1) then
                    write(ioutid,*) "Localization index for alpha spin:"
                    do iatm=1,ncenter
                        write(ioutid,"(i5,'(',a,'):',f7.3)",advance='no') iatm,ind2name(a(iatm)%index),LIa(iatm)
                        if (mod(iatm,4)==0) write(ioutid,*)
                    end do
                    write(ioutid,*)
                    write(ioutid,*)
                end if
                !Beta
                call showmatgau(DIb,"Delocalization index matrix for beta spin",0,"f14.8",ioutid)
                write(ioutid,*)
                if (iwork/=1) then
                    write(ioutid,*) "Localization index for beta spin:"
                    do iatm=1,ncenter
                        write(ioutid,"(i5,'(',a,'):',f7.3)",advance='no') iatm,ind2name(a(iatm)%index),LIb(iatm)
                        if (mod(iatm,4)==0) write(ioutid,*)
                    end do
                    write(ioutid,*)
                    write(ioutid,*)
                end if
            end if
            !Alpha+Beta
            call showmatgau(DI,"Total delocalization index matrix",0,"f14.8",ioutid)
            write(ioutid,*)
            if (iwork/=1) then
                write(ioutid,*) "Localization index:"
                do iatm=1,ncenter
                    write(ioutid,"(i5,'(',a,'):',f7.3)",advance='no') iatm,ind2name(a(iatm)%index),LI(iatm)
                    if (mod(iatm,4)==0) write(ioutid,*)
                end do
                write(ioutid,*)
                write(ioutid,*)
            end if
            if (selectyn=='n') then !Just output result to screen, choose if output result to plain text file
                write(*,*) "If also output LI and DI to LIDI.txt in current folder? (y/n)"
                read(*,*) selectyn
                if (selectyn=='y') then
                    open(10,file="LIDI.txt",status="replace")
                    ioutid=10
                else if (selectyn=='n') then
                    exit
                end if
            else if (selectyn=='y') then !Have already outputted result to LIDI.txt, exit cycle
                write(*,*) "Done, the LI and DI have been outputted to LIDI.txt in current folder"
                close(10)
                exit
            end if
        end do
        
    else if (iwork==1) then !Output fuzzy bond order
        write(*,"('The total bond order >=',f10.6)") bndordthres
        itmp=0
        if (wfntype==1.or.wfntype==2.or.wfntype==4) then
            do i=1,ncenter
                do j=i+1,ncenter
                    if (DIa(i,j)+DIb(i,j)>=bndordthres) then
                        itmp=itmp+1
                        write(*,"('#',i5,':',i5,a,i5,a,' Alpha: ',f10.6,' Beta:',f10.6,' Total:',f10.6)") &
                        itmp,i,'('//a(i)%name//')',j,'('//a(j)%name//')',DIa(i,j),DIb(i,j),DIa(i,j)+DIb(i,j)
                    end if
                end do
            end do
        else if (wfntype==0.or.wfntype==3) then
            itmp=0
            do i=1,ncenter
                do j=i+1,ncenter
                    if (DI(i,j)>=bndordthres) then
                        itmp=itmp+1
                        write(*,"('#',i5,':',5x,i5,a,i5,a,f14.8)") itmp,i,'('//a(i)%name//')',j,'('//a(j)%name//')',DI(i,j)
                    end if
                end do
            end do
        end if
        !Output bond order sequence defined by oliarr
!         do itmp=1,size(oliarr)-1
!             write(*,"(2i6,f12.6)") oliarr(itmp),oliarr(itmp+1),DI(oliarr(itmp),oliarr(itmp+1))
!         end do
        if (allocated(frag1)) then
            bndordfraga=0
            bndordfragb=0
            bndordfragtot=0
            do i=1,size(frag1)
                do j=1,size(frag2)
                    if (wfntype==1.or.wfntype==2.or.wfntype==4) then 
                        bndordfraga=bndordfraga+DIa(frag1(i),frag2(j))
                        bndordfragb=bndordfragb+DIb(frag1(i),frag2(j))
                    else if (wfntype==0.or.wfntype==3) then
                        bndordfragtot=bndordfragtot+DI(frag1(i),frag2(j))
                    end if
                end do
            end do
            write(*,*)
            if (wfntype==1.or.wfntype==2.or.wfntype==4) then
                write(*,"('The bond order between fragment 1 and 2:')")
                write(*,"('Alpha:',f10.6,' Beta:',f10.6,' Total:',f10.6)") bndordfraga,bndordfragb,bndordfraga+bndordfragb
            else if (wfntype==0.or.wfntype==3) then
                write(*,"('The bond order between fragment 1 and 2:',f12.6)") bndordfragtot
            end if
        end if
        write(*,*)
        write(*,*) "If output bond order matrix? (y/n)"
        read(*,*) selectyn
        if (selectyn=='y'.or.selectyn=='Y') then
            open(10,file="bndmat.txt",status="replace")
            if (wfntype==1.or.wfntype==2.or.wfntype==4) then !UHF,ROHF,U-post-HF, output each spin component first
                call showmatgau(DIa,"Delocalization index matrix for alpha spin",0,"f14.8",10)
                write(10,*)
                call showmatgau(DIb,"Delocalization index matrix for beta spin",0,"f14.8",10)
                write(10,*)
            end if
            call showmatgau(DI,"Total delocalization index matrix",0,"f14.8",10)
            write(10,*)
            close(10)
            write(*,*) "Done, bond order matrix has been outputted to bndmat.txt in current folder"
            write(*,"(a)") " Note: Diagonal terms in the bond order matrix are the sum of corresponding row or column elements, for close-shell cases, they are also known as atomic valence"
        end if
        radpot=radpotold
        sphpot=sphpotold
        return !Fuzzy bond order has been shown, now (normally) return to bond order analysis interface
    end if
    
else if (isel==5) then !PDI
    call showmatgau(DI,"Delocalization index matrix",0,"f14.8")
    write(*,"(a)") " Note: Diagonal terms are the sum of corresponding row or column elements, for close-shell cases, these also known as atomic valence"
    write(*,*)
    do while(.true.)
        write(*,"(a)") " Input indices of the six atoms constituting the ring, in clockwise or anti-clockwise. e.g. 4,5,6,7,8,2"
        write(*,*) "(Input q can return)"
        read(*,"(a)") c80inp
        if (c80inp(1:1)=='q'.or.c80inp(1:1)=='Q') exit
        read(c80inp,*) PDIatom(:)
        write(*,"(' Delocalization index of ',i5,'(',a,')   --',i5,'(',a,'):',f12.6)") PDIatom(1),a(PDIatom(1))%name,PDIatom(4),a(PDIatom(4))%name,DI(PDIatom(1),PDIatom(4))
        write(*,"(' Delocalization index of ',i5,'(',a,')   --',i5,'(',a,'):',f12.6)") PDIatom(2),a(PDIatom(2))%name,PDIatom(5),a(PDIatom(5))%name,DI(PDIatom(2),PDIatom(5))
        write(*,"(' Delocalization index of ',i5,'(',a,')   --',i5,'(',a,'):',f12.6)") PDIatom(3),a(PDIatom(3))%name,PDIatom(6),a(PDIatom(6))%name,DI(PDIatom(3),PDIatom(6))
        write(*,"(' PDI value is',f12.6)") ( DI(PDIatom(1),PDIatom(4))+DI(PDIatom(2),PDIatom(5))+DI(PDIatom(3),PDIatom(6)) )/3D0
        write(*,*)
    end do
    
else if (isel==6) then !FLU
    call showmatgau(DI,"Delocalization index matrix",0,"f14.8")
    write(*,"(a)") " Note: Diagonal terms are the sum of corresponding row or column elements, for close-shell cases, these also known as atomic valence"
    write(*,*)
    write(*,*) "Current FLU reference parameters:"
    do iref=1,nelesupp
        do jref=iref,nelesupp
            if (FLUref(iref,jref)/=-1) write(*,"(' ',a,a,a,a,f10.5)") ind2name(iref),'-',ind2name(jref),':',FLUref(iref,jref)
        end do
    end do
    write(*,*)
    do while(.true.)
        write(*,"(a)") " Input indices of the atoms in the ring, in clockwise or anti-clockwise"
        write(*,*) "e.g. 4,7,8,1,2,3      (Input q can exit)"
        read(*,"(a)") c80inp
        if (c80inp(1:1)=='q'.or.c80inp(1:1)=='Q') exit
        call str2arr(c80inp,nFLUatom,FLUatom)
        FLUval=0D0
        write(*,*) "        Atom pair         Contribution          DI"
        do iidx=1,nFLUatom
            jidx=iidx+1
            if (iidx==nFLUatom) jidx=1 !Return to the first element of the ring after a cycle
            iatm=FLUatom(iidx) !Actual atom index in present system
            jatm=FLUatom(jidx)
            iatmeleidx=a(iatm)%index !Index in periodic table
            jatmeleidx=a(jatm)%index
            refval=FLUref(iatmeleidx,jatmeleidx)
            if (refval==-1D0) then
                write(*,"(' Error: Missing reference parameter for',a,'-',a)") ind2name(iatmeleidx),ind2name(jatmeleidx)
                exit
            end if
            valenratio=DI(iatm,iatm)/DI(jatm,jatm) !DI(iatm,iatm) is the sum of corresponding row or column elements, namely atomic valence, rather than LI*2 of iatm
            if (valenratio<1) valenratio=1D0/valenratio
            FLUpair=(valenratio*( (DI(iatm,jatm)-refval)/refval ))**2/nFLUatom
            write(*,"(i5,'(',a,')  --',i5,'(',a,'):',f15.6,f16.6)") iatm,ind2name(iatmeleidx),jatm,ind2name(jatmeleidx),FLUpair,DI(iatm,jatm)
            FLUval=FLUval+FLUpair
        end do
        write(*,"(' FLU value is',f12.6)") FLUval
        write(*,*)
    end do
    
else if (isel==7) then !FLU-pi
    write(*,*) "Which occupied orbitals are pi orbitals? Input their serials, e.g. 17,20,21"
    read(*,"(a)") c80inp
    call str2arr(c80inp,nFLUorb,FLUorb)
    !Generate DI for pi orbitals. DI_A,B=2¡Æ[i,j]dsqrt(n_i*n_j)*S_i,j_A * S_i,j_B     where i and j are non-spin orbital
    DI=0D0
    do iatm=1,ncenter
        do jatm=iatm+1,ncenter
            tmpval=0D0
            do iidx=1,nFLUorb
                iorb=FLUorb(iidx)
                do jidx=1,nFLUorb
                    jorb=FLUorb(jidx)
                    tmpval=tmpval+dsqrt(MOocc(iorb)*MOocc(jorb))*AOM(iorb,jorb,iatm)*AOM(iorb,jorb,jatm)
                end do
            end do
            DI(iatm,jatm)=tmpval
        end do
    end do
    DI=2*(DI+transpose(DI))
    do iatm=1,ncenter !Calculate atomic valence
        DI(iatm,iatm)=sum(DI(iatm,:))
    end do
    call showmatgau(DI,"Delocalization index matrix for pi electrons",0,"f14.8")
    write(*,"(a)") " Note: Diagonal terms are the sum of corresponding row or column elements"
    write(*,*)
    do while(.true.)
        write(*,"(a)") " Input indices of the atoms in the ring, in clockwise or anti-clockwise"
        write(*,*) "e.g. 4,7,8,1,2,3      (Input q can exit)"
        read(*,"(a)") c80inp
        if (c80inp(1:1)=='q'.or.c80inp(1:1)=='Q') exit
        call str2arr(c80inp,nFLUatom,FLUatom)
        !Calculate average of DI-oi first
        avgDI=0D0
        do iidx=1,nFLUatom
            jidx=iidx+1
            if (iidx==nFLUatom) jidx=1
            avgDI=avgDI+DI(FLUatom(iidx),FLUatom(jidx))
        end do
        avgDI=avgDI/nFLUatom
        write(*,"(' Average of DI-pi is',f12.6)") avgDI
        FLUval=0D0
        write(*,*) "        Atom pair         Contribution          DI"
        do iidx=1,nFLUatom
            jidx=iidx+1
            if (iidx==nFLUatom) jidx=1
            iatm=FLUatom(iidx) !Actual atom index in present system
            jatm=FLUatom(jidx)
            iatmeleidx=a(iatm)%index !Index in periodic table
            jatmeleidx=a(jatm)%index
            valenratio=DI(iatm,iatm)/DI(jatm,jatm)
            if (valenratio<1) valenratio=1D0/valenratio
            FLUpair=(valenratio*(DI(iatm,jatm)-avgDI)/avgDI)**2/nFLUatom
            write(*,"(i5,'(',a,')  --',i5,'(',a,'):',f15.6,f16.6)") iatm,ind2name(iatmeleidx),jatm,ind2name(jatmeleidx),FLUpair,DI(iatm,jatm)
            FLUval=FLUval+FLUpair
        end do
        write(*,"(' FLU-pi value is',f12.6)") FLUval
        write(*,*)
    end do
    
else if (isel==8) then !Integral in overlap region
    ovlpintpos=ovlpintpos+transpose(ovlpintpos)
    ovlpintneg=ovlpintneg+transpose(ovlpintneg)
    sumdiagpos=0D0
    sumdiagneg=0D0
    do i=1,ncenter
        ovlpintpos(i,i)=ovlpintpos(i,i)/2D0
        sumdiagpos=sumdiagpos+ovlpintpos(i,i)
        ovlpintneg(i,i)=ovlpintneg(i,i)/2D0
        sumdiagneg=sumdiagneg+ovlpintneg(i,i)
    end do
    ovlpinttot=ovlpintpos+ovlpintneg
    if (iwork==2) then !Output Laplacian bond order
        write(*,"('The bond order >=',f10.6)") bndordthres
        itmp=0
        do i=1,ncenter
            do j=i+1,ncenter
                if (-10*ovlpintneg(i,j)>=bndordthres) then
                    itmp=itmp+1
                    write(*,"('#',i5,':',i5,a,i5,a,':',f10.6)") &
                    itmp,i,'('//a(i)%name//')',j,'('//a(j)%name//')',-10*ovlpintneg(i,j)
                end if
            end do
        end do
        if (allocated(frag1)) then !Output interfragment bond order
            bndordfragtot=0
            do i=1,size(frag1)
                do j=1,size(frag2)
                    bndordfragtot=bndordfragtot-10*ovlpintneg(frag1(i),frag2(j))
                end do
            end do
            write(*,*)
            write(*,"('The bond order between fragment 1 and 2:',f12.6)") bndordfragtot
        end if
        write(*,*)
        write(*,*) "If output bond order matrix? (y/n)"
        read(*,*) selectyn
        if (selectyn=='y'.or.selectyn=='Y') then
            open(10,file="bndmat.txt",status="replace")
            do i=1,ncenter
                ovlpintneg(i,i)=sum(ovlpintneg(i,:))-ovlpintneg(i,i) !Make diagonal terms are the sum of corresponding row elements, namely valence
            end do
            call showmatgau(-10*ovlpintneg,"Laplacian bond order matrix",0,"f14.8",10)
            close(10)
            write(*,*) "Done, bond order matrix has been outputted to bndmat.txt in current folder"
            write(*,"(a)") " Note: Diagonal terms in the bond order matrix are the sum of corresponding row or column elements"
            write(*,*)
        end if
        radpot=radpotold
        sphpot=sphpotold
        return !Laplacian bond order has been shown, now (normally) return to bond order analysis interface
    else
        call showmatgau(ovlpintpos,"Integration of positive values in overlap region",0,"f14.8")
        sumovlppos=sum(ovlpintpos)
        write(*,"('Summing up diagonal matrix elements:     ',f20.8)") sumdiagpos
        write(*,"('Summing up non-diagonal, matrix elements:',f20.8)") sumovlppos-sumdiagpos
        write(*,"('Summing up all matrix elements:          ',f20.8)") sumovlppos
        write(*,*)
        sumovlpneg=sum(ovlpintneg)
        call showmatgau(ovlpintneg,"Integration of negative values in overlap region",0,"f14.8")
        write(*,"('Summing up diagonal matrix elements:     ',f20.8)") sumdiagneg
        write(*,"('Summing up non-diagonal, matrix elements:',f20.8)") sumovlpneg-sumdiagneg
        write(*,"('Summing up all matrix elements:          ',f20.8)") sumovlpneg
        write(*,*)
        sumovlptot=sum(ovlpinttot)
        call showmatgau(ovlpinttot,"Integration of all values in overlap region",0,"f14.8")
        write(*,"('Summing up diagonal matrix elements:     ',f20.8)") sumdiagpos+sumdiagneg
        write(*,"('Summing up non-diagonal, matrix elements:',f20.8)") sumovlptot-sumdiagpos-sumdiagneg
        write(*,"('Summing up all matrix elements:          ',f20.8)") sumovlptot
        write(*,*)
        write(*,*) "If also output above matrices to intovlp.txt in current folder? (y/n)"
        read(*,*) selectyn
        if (selectyn=='y'.or.selectyn=='Y') then
            open(10,file="intovlp.txt",status="replace")
            call showmatgau(ovlpintpos,"Integration of positive values in overlap region",0,"f14.8",10)
            write(10,*)
            call showmatgau(ovlpintneg,"Integration of negative values in overlap region",0,"f14.8",10)
            write(10,*)
            call showmatgau(ovlpinttot,"Integration of all values in overlap region",0,"f14.8",10)
            write(10,*)
            close(10)
            write(*,*) "Done, the matrices have been outputted to intovlp.txt in current folder"
            write(*,*)
        end if
    end if
    
else if (isel==9) then !CLRK
    call showmatgau(CLRK,"Condensed linear response kernel (CLRK) matrix",0,"f14.8")
    write(*,*)
    write(*,*) "If also output CLRK to CLRK.txt in current folder? (y/n)"
    read(*,*) selectyn
    if (selectyn=='y'.or.selectyn=='Y') then
        open(10,file="CLRK.txt",status="replace")
        call showmatgau(CLRK,"Condensed linear response kernel (CLRK) matrix",0,"f14.8",10)
        close(10)
        write(*,*) "Done, the CLRK matrix has been outputted to CLRK.txt in current folder"
        write(*,*)
    end if
    
else if (isel==10) then !PLR
    call showmatgau(CLRK,"Condensed linear response kernel (CLRK) matrix",0,"f14.8")
    write(*,*)
    do while(.true.)
        write(*,"(a)") " Input indices of the six atoms constituting the ring, in clockwise or anti-clockwise. e.g. 4,5,6,7,8,2"
        write(*,*) "(Input q can return)"
        read(*,"(a)") c80inp
        if (c80inp(1:1)=='q'.or.c80inp(1:1)=='Q') exit
        read(c80inp,*) PLRatom(:)
        write(*,"(' CLRK of ',i5,'(',a,')   --',i5,'(',a,'):',f12.6)") PLRatom(1),a(PLRatom(1))%name,PLRatom(4),a(PLRatom(4))%name,CLRK(PLRatom(1),PLRatom(4))
        write(*,"(' CLRK of ',i5,'(',a,')   --',i5,'(',a,'):',f12.6)") PLRatom(2),a(PLRatom(2))%name,PLRatom(5),a(PLRatom(5))%name,CLRK(PLRatom(2),PLRatom(5))
        write(*,"(' CLRK of ',i5,'(',a,')   --',i5,'(',a,'):',f12.6)") PLRatom(3),a(PLRatom(3))%name,PLRatom(6),a(PLRatom(6))%name,CLRK(PLRatom(3),PLRatom(6))
        write(*,"(' PLR index is',f12.6)") ( CLRK(PLRatom(1),PLRatom(4))+CLRK(PLRatom(2),PLRatom(5))+CLRK(PLRatom(3),PLRatom(6)) )/3D0
        write(*,*)
    end do
    
else if (isel==11) then !Multicenter DI
    do while(.true.)
        write(*,*) "Input atom indices, e.g. 3,4,7,8,10    (Up to 10 atoms)"
        write(*,*) "Input q can return to upper level menu"
        read(*,"(a)") c80inp
        if (c80inp(1:1)=='q') then
            exit
        else
            call str2arr(c80inp,nDIcen,cenind)
        end if
        DImulti=0D0
        write(*,*) "Please wait..."
        if (nDIcen==3) then
            do iorb=1,nmatsize
                if (MOocc(iorb)==0D0) cycle
                do jorb=1,nmatsize
                    if (MOocc(jorb)==0D0) cycle
                    do korb=1,nmatsize
                        if (MOocc(korb)==0D0) cycle
            DImulti=DImulti+&
            AOM(iorb,jorb,cenind(1))*AOM(jorb,korb,cenind(2))*AOM(korb,iorb,cenind(3))
                    end do
                end do
            end do
        else if (nDIcen==4) then
            do iorb=1,nmatsize
                if (MOocc(iorb)==0D0) cycle
                do jorb=1,nmatsize
                    if (MOocc(jorb)==0D0) cycle
                    do korb=1,nmatsize
                        if (MOocc(korb)==0D0) cycle
                        do lorb=1,nmatsize
                            if (MOocc(lorb)==0D0) cycle
            DImulti=DImulti+&
            AOM(iorb,jorb,cenind(1))*AOM(jorb,korb,cenind(2))*AOM(korb,lorb,cenind(3))*AOM(lorb,iorb,cenind(4))
                        end do
                    end do
                end do
            end do
        else if (nDIcen==5) then
            do iorb=1,nmatsize
                if (MOocc(iorb)==0D0) cycle
                do jorb=1,nmatsize
                    if (MOocc(jorb)==0D0) cycle
                    do korb=1,nmatsize
                        if (MOocc(korb)==0D0) cycle
                        do lorb=1,nmatsize
                            if (MOocc(lorb)==0D0) cycle
                            do morb=1,nmatsize
                                if (MOocc(morb)==0D0) cycle
            DImulti=DImulti+&
            AOM(iorb,jorb,cenind(1))*AOM(jorb,korb,cenind(2))*AOM(korb,lorb,cenind(3))*AOM(lorb,morb,cenind(4))*AOM(morb,iorb,cenind(5))
                            end do
                        end do
                    end do
                end do
            end do                
        else if (nDIcen==6) then
            do iorb=1,nmatsize
                if (MOocc(iorb)==0D0) cycle
                do jorb=1,nmatsize
                    if (MOocc(jorb)==0D0) cycle
                    do korb=1,nmatsize
                        if (MOocc(korb)==0D0) cycle
                        do lorb=1,nmatsize
                            if (MOocc(lorb)==0D0) cycle
                            do morb=1,nmatsize
                                if (MOocc(morb)==0D0) cycle
                                do norb=1,nmatsize
                                    if (MOocc(norb)==0D0) cycle
            DImulti=DImulti+& !dsqrt(MOocc(iorb)*MOocc(jorb)*MOocc(korb)*MOocc(lorb)*MOocc(morb)*MOocc(norb))*
            AOM(iorb,jorb,cenind(1))*AOM(jorb,korb,cenind(2))*AOM(korb,lorb,cenind(3))*AOM(lorb,morb,cenind(4))*AOM(morb,norb,cenind(5))*AOM(norb,iorb,cenind(6))
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        else if (nDIcen==7) then
            do iorb=1,nmatsize
                write(*,"(' Finished',i8,'/',i8)") iorb,nmatsize
                if (MOocc(iorb)==0D0) cycle
                do jorb=1,nmatsize
                    if (MOocc(jorb)==0D0) cycle
                    do korb=1,nmatsize
                        if (MOocc(korb)==0D0) cycle
                        do lorb=1,nmatsize
                            if (MOocc(lorb)==0D0) cycle
                            do morb=1,nmatsize
                                if (MOocc(morb)==0D0) cycle
                                do norb=1,nmatsize
                                    if (MOocc(norb)==0D0) cycle
                                    do iiorb=1,nmatsize
                                        if (MOocc(iiorb)==0D0) cycle
            DImulti=DImulti+&
            AOM(iorb,jorb,cenind(1))*AOM(jorb,korb,cenind(2))*AOM(korb,lorb,cenind(3))*AOM(lorb,morb,cenind(4))*AOM(morb,norb,cenind(5))*&
            AOM(norb,iiorb,cenind(6))*AOM(iiorb,iorb,cenind(7))
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        else if (nDIcen==8) then
            do iorb=1,nmatsize
                write(*,"(' Finished',i8,'/',i8)") iorb,nmatsize
                if (MOocc(iorb)==0D0) cycle
                do jorb=1,nmatsize
                    if (MOocc(jorb)==0D0) cycle
                    do korb=1,nmatsize
                        if (MOocc(korb)==0D0) cycle
                        do lorb=1,nmatsize
                            if (MOocc(lorb)==0D0) cycle
                            do morb=1,nmatsize
                                if (MOocc(morb)==0D0) cycle
                                do norb=1,nmatsize
                                    if (MOocc(norb)==0D0) cycle
                                    do iiorb=1,nmatsize
                                        if (MOocc(iiorb)==0D0) cycle
                                        do jjorb=1,nmatsize
                                            if (MOocc(jjorb)==0D0) cycle
            DImulti=DImulti+&
            AOM(iorb,jorb,cenind(1))*AOM(jorb,korb,cenind(2))*AOM(korb,lorb,cenind(3))*AOM(lorb,morb,cenind(4))*AOM(morb,norb,cenind(5))*&
            AOM(norb,iiorb,cenind(6))*AOM(iiorb,jjorb,cenind(7))*AOM(jjorb,iorb,cenind(8))
                                        end do
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        else if (nDIcen==9) then
            do iorb=1,nmatsize
                write(*,"(' Finished',i8,'/',i8)") iorb,nmatsize
                if (MOocc(iorb)==0D0) cycle
                do jorb=1,nmatsize
                    if (MOocc(jorb)==0D0) cycle
                    do korb=1,nmatsize
                        if (MOocc(korb)==0D0) cycle
                        do lorb=1,nmatsize
                            if (MOocc(lorb)==0D0) cycle
                            do morb=1,nmatsize
                                if (MOocc(morb)==0D0) cycle
                                do norb=1,nmatsize
                                    if (MOocc(norb)==0D0) cycle
                                    do iiorb=1,nmatsize
                                        if (MOocc(iiorb)==0D0) cycle
                                        do jjorb=1,nmatsize
                                            if (MOocc(jjorb)==0D0) cycle
                                            do kkorb=1,nmatsize
                                                if (MOocc(kkorb)==0D0) cycle
            DImulti=DImulti+&
            AOM(iorb,jorb,cenind(1))*AOM(jorb,korb,cenind(2))*AOM(korb,lorb,cenind(3))*AOM(lorb,morb,cenind(4))*AOM(morb,norb,cenind(5))*&
            AOM(norb,iiorb,cenind(6))*AOM(iiorb,jjorb,cenind(7))*AOM(jjorb,kkorb,cenind(8))*AOM(kkorb,iorb,cenind(9))
                                            end do
                                        end do
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        else if (nDIcen==10) then
            do iorb=1,nmatsize
                write(*,"(' Finished',i8,'/',i8)") iorb,nmatsize
                if (MOocc(iorb)==0D0) cycle
                do jorb=1,nmatsize
                    if (MOocc(jorb)==0D0) cycle
                    do korb=1,nmatsize
                        if (MOocc(korb)==0D0) cycle
                        do lorb=1,nmatsize
                            if (MOocc(lorb)==0D0) cycle
                            do morb=1,nmatsize
                                if (MOocc(morb)==0D0) cycle
                                do norb=1,nmatsize
                                    if (MOocc(norb)==0D0) cycle
                                    do iiorb=1,nmatsize
                                        if (MOocc(iiorb)==0D0) cycle
                                        do jjorb=1,nmatsize
                                            if (MOocc(jjorb)==0D0) cycle
                                            do kkorb=1,nmatsize
                                                if (MOocc(kkorb)==0D0) cycle
                                                do llorb=1,nmatsize
                                                    if (MOocc(llorb)==0D0) cycle
            DImulti=DImulti+&
            AOM(iorb,jorb,cenind(1))*AOM(jorb,korb,cenind(2))*AOM(korb,lorb,cenind(3))*AOM(lorb,morb,cenind(4))*AOM(morb,norb,cenind(5))*&
            AOM(norb,iiorb,cenind(6))*AOM(iiorb,jjorb,cenind(7))*AOM(jjorb,kkorb,cenind(8))*AOM(kkorb,llorb,cenind(9))*AOM(llorb,iorb,cenind(10))
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
        DImulti=DImulti*2**(nDIcen-1)
        write(*,"(' Multicenter DI:',f13.7,/)") DImulti
        write(*,"(' Multicenter DI in normalized form: ',f13.7,/)") DImulti**(1D0/nDIcen)
    end do
end if
    
end do !End interface loop

end subroutine





