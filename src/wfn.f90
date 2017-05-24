program multiwfn
use defvar
use util
use function
use topo
implicit real*8(a-h,o-z)
character nowdate*20,nowtime*20,inpstring*80,c200tmp*200,c2000tmp*2000,outcubfile*200,selectyn,lovername*80,settingpath*200
real*8 :: inx,iny,inz,tmpvec(3)
integer :: iprintfunc=1 !The default function whose gradient and Hessian will be outputted at a point by main function 1
integer,allocatable :: exclfragatm(:),tmparrint(:)
integer walltime1,walltime2
real*8,allocatable :: d1add(:,:),d1min(:,:),d2add(:,:),d2min(:,:),d1addtmp(:,:),d1mintmp(:,:),d2addtmp(:,:),d2mintmp(:,:) !Store temporary data for drawing gradient map
real*8,allocatable :: planemat_cust(:,:) !For storing temporary data of doing custom map
real*8,allocatable :: planemat_bk(:,:) !Used to backup plane data
! real*8,allocatable :: tmpmat(:,:)

call getarg(1,filename)
call getarg(2,cmdarg2)
11 call loadsetting
if (isys==1) write(*,*) "Multiwfn -- A Multifunctional Wavefunction Analyzer (for Windows 64bit)"
if (isys==2) write(*,*) "Multiwfn -- A Multifunctional Wavefunction Analyzer (for Linux 64bit)"
if (isys==3) write(*,*) "Multiwfn -- A Multifunctional Wavefunction Analyzer (for MacOS)"
write(*,*) "Version 3.4(dev), release date: 2017-May-24"
write(*,"(a)") " Project leader: Tian Lu (Beijing Kein Research Center for Natural Sciences)"
write(*,*) "Citation of Multiwfn: Tian Lu, Feiwu Chen, J. Comput. Chem. 33, 580-592 (2012)"
write(*,*) "Multiwfn official website: http://sobereva.com/multiwfn"
write(*,*) "Multiwfn official forum (in Chinese): http://bbs.keinsci.com"
write(*,*) "Bug reporting, question and suggestion, please contact: Sobereva@sina.com"

!!!!!!!!
!if (isys==1) call KMP_SET_STACKSIZE_S(ompstacksize) !For Linux/MacOSX version, it seems the only way to set stacksize of each thread is to define KMP_STACKSIZE environment variable
!!!!!!!!

nthreads=getNThreads()
call date_and_time(nowdate,nowtime)
write(*,"(' ( The number of threads:',i3,'   Current date: ',a,'-',a,'-',a,'   Time: ',a,':',a,':',a,' )')") &
nthreads,nowdate(1:4),nowdate(5:6),nowdate(7:8),nowtime(1:2),nowtime(3:4),nowtime(5:6)
write(*,*)

! call system("echo %GAUSS_EXEDIR%")
! call selfileGUI !Doesn't work, I don't know why
if (trim(filename)=="") then !Haven't defined filename variable
    call mylover(lovername)
    write(*,"(a,a,a)") " Input file path, for example E:\",trim(lovername),".wfn"
    write(*,*) "(Supported types: .wfn/.wfx/.fch/.31/.chg/.pdb/.xyz/.cub/.grd/.molden, etc.)"
    write(*,"(a)") " Hint: To reload the file last time used, simply input the letter ""o"". Input such as ?miku.fch can open miku.fch in the same folder of the file last time used"
    do while(.true.)
        read(*,"(a)") filename
        ltmp=len_trim(filename)
        if (ltmp==0) cycle
        if (filename=='o') then
            write(*,"(' The file last time used: ',a)") trim(lastfile)
            filename=lastfile
        end if
        !Remove the first and the last " or  'symbol, because directly dragging file into the window will result in " or ' symbol, which is unrecognized by Multwifn
        if (filename(1:1)=='"'.or.filename(1:1)=="'") filename(1:1)=" "
        if (filename(ltmp:ltmp)=='"'.or.filename(ltmp:ltmp)=="'") filename(ltmp:ltmp)=" "
        if (filename(1:1)=='?') then
            do itmp=len_trim(lastfile),1,-1
                if (isys==1.and.lastfile(itmp:itmp)=='\') exit
                if ((isys==2.or.isys==3).and.lastfile(itmp:itmp)=='/') exit
            end do
            filename=lastfile(1:itmp)//trim(filename(2:))
        end if
        inquire(file=filename,exist=alive)
        if (alive.eqv..true.) exit
        write(*,"('""',a,'"" ',a)") trim(filename),"cannot be found, input again"
    end do
    !Write current opened file to "lastfile" in settings.ini
    inquire(file="settings.ini",exist=alive)
    if (alive .eqv. .true.) then
        settingpath="settings.ini"
    else if (alive .eqv. .false.) then
        call getenv("Multiwfnpath",c200tmp)
        if (isys==1) then
            settingpath=trim(c200tmp)//"\settings.ini"
        else if (isys==2.or.isys==3) then
            settingpath=trim(c200tmp)//"/settings.ini"
        end if
    end if
    inquire(file=settingpath,exist=alive)
    if (alive) then
        open(20,file=settingpath,status="old")
        call loclabel(20,"lastfile")
        write(20,"(a)") "lastfile= "//trim(filename)
        close(20)
    end if
else
    inquire(file=filename,exist=alive)
    if (alive.eqv..false.) then
        write(*,*) "File not found, exit program..."
        pause
        stop
    end if
end if
call readinfile(filename,0)

!!-- Backup various information of first loaded (meanwhile unmodified) molecule
firstfilename=filename
allocate(a_org(ncenter))
allocate(b_org(nprims))
allocate(CO_org(nmo,nprims))
allocate(MOocc_org(nmo))
a_org=a
b_org=b
CO_org=CO
MOocc_org=MOocc
nprims_org=nprims
nmo_org=nmo
ncenter_org=ncenter

!!-- Initialize fragment
nfragatmnum=ncenter !Default fragment is the whole molecule
nfragatmnumbackup=ncenter
allocate(fragatm(nfragatmnum),fragatmbackup(nfragatmnum))
forall (i=1:nfragatmnum) fragatm(i)=i
forall (i=1:nfragatmnum) fragatmbackup(i)=i
ifragcontri=0

!!-- Call some routine only once
if (ncenter>5000) then
    write(*,"(a)") "Warning: There are too many atoms, the distance matrix cannot be generated! Some functions may not work properly"
else
    call gendistmat !Generate distance matrix
end if
!Convert prebuild radii from Angstrom to Bohr. But some radii such as radii_hugo will remain unchanged since it is recorded as Bohr
if (ifirstMultiwfn==1) then
    vdwr=vdwr/b2a
    vdwr_tianlu=vdwr_tianlu/b2a
    covr=covr/b2a
    covr_Suresh=covr_Suresh/b2a
    covr_pyy=covr_pyy/b2a
    covr_tianlu=covr_tianlu/b2a
end if
!Only get into dislin level 1 to get width and height of screen in pixels, don't do any other things
!-- Show related molecular information
if (ifiletype/=0.and.ifiletype/=8) then
    call showformula
    totmass=sum(atmwei(a%index))
    write(*,"(' Molecule weight:',f16.5)") totmass
end if

write(*,"(/,3a)") " Loaded ",trim(filename)," successfully!"
! call sys1eprop !Show some system 1e properties, only works when Cartesian basis functions are presented

! call intisosurface


!!!--------------------- Now everything start ---------------------!!!
!!!--------------------- Now everything start ---------------------!!!
!!!--------------------- Now everything start ---------------------!!!
do while(.true.) !Main loop

10 write(*,*)
if (allocated(cubmat)) write(*,*) "Note: A set of grid data presents in memory"
write(*,*) "                   ------------ Main function menu ------------"
! write(*,*) "-11 Load a new file"
! write(*,*) "-10 Exit program"
if (ifragcontri/=1) write(*,*) "-4 Exclude some atoms contribution to property"
if (ifragcontri/=1) write(*,*) "-3 Obtain a fragment contribution to property"
if (ifiletype/=7.and.ifiletype/=8) write(*,*) "0 Show molecular structure and view orbitals"
if (ifiletype==7.or.ifiletype==8) write(*,*) "0 Show molecular structure and view isosurface"
write(*,*) "1 Output all properties at a point"
write(*,*) "2 Topology analysis"
write(*,*) "3 Output and plot specific property in a line"
write(*,*) "4 Output and plot specific property in a plane"
write(*,*) "5 Output and plot specific property within a spatial region (calc. grid data)"
write(*,*) "6 Check & modify wavefunction"
write(*,*) "7 Population analysis"
write(*,*) "8 Orbital composition analysis"
write(*,*) "9 Bond order analysis"
write(*,*) "10 Plot Total/Partial/Overlap population density-of-states (DOS)"
write(*,*) "11 Plot IR/Raman/UV-Vis/ECD/VCD spectrum"
write(*,*) "12 Quantitative analysis of molecular surface"
if (allocated(cubmat)) write(*,*) "13 Process grid data"
if (.not.allocated(cubmat)) write(*,*) "13 Process grid data (No grid data is presented currently)"
write(*,*) "14 Adaptive natural density partitioning (AdNDP) analysis"
write(*,*) "15 Fuzzy atomic space analysis"
write(*,*) "16 Charge decomposition analysis (CDA) and extended CDA (ECDA)"
write(*,*) "17 Basin analysis"
write(*,*) "18 Electron excitation analysis"
write(*,*) "100 Other functions (Part1)"
write(*,*) "200 Other functions (Part2)"
! write(*,*) "1000 Set some parameters"
read(*,*) infuncsel1

!!! Setting various content before implement formal functions

if (infuncsel1==-10) then !Exit program
    stop
else if (infuncsel1==-11) then !Load a new file
    call dealloall
    filename=""
    deallocate(a_org,b_org,CO_org,MOocc_org,fragatm,fragatmbackup)
    ifirstMultiwfn=0
    goto 11
else if (infuncsel1==-3.or.infuncsel1==-4) then
    deallocate(fragatm) !fragatm has been defined previously by default, fragatm contains all atoms
    if (infuncsel1==-3) then
        ! "fragatm" is convertion relationship from fragment to the whole,
        ! e.g. fragatm(4) is the actual atom index corresponding the 4th atom in fragment list
        write(*,"(a)") " Input atomic indices to define the fragment, e.g. 1,3-6,8,10-11 means the atoms 1,3,4,5,6,8,10,11 will constitute the fragment"
        read(*,"(a)") c2000tmp
        call str2arr(c2000tmp,nfragatmnum)
        allocate(fragatm(nfragatmnum))
        call str2arr(c2000tmp,nfragatmnum,fragatm)
        call sorti4(fragatm,"val")
    else if(infuncsel1==-4) then
        write(*,*) "How many atoms will be excluded?"
        write(*,*) "e.g. 1,3-6,8,10-11 means the atoms 1,3,4,5,6,8,10,11 will be excluded"
        read(*,"(a)") c2000tmp
        call str2arr(c2000tmp,nexclatm)
        nfragatmnum=ncenter-nexclatm
        allocate(fragatm(nfragatmnum),exclfragatm(nexclatm))
        call str2arr(c2000tmp,nexclatm,exclfragatm)
        j=0
        do i=1,ncenter
            if (all(exclfragatm/=i)) then
                j=j+1
                fragatm(j)=i
            end if
        end do
    end if
    j=0
    do i=1,nprims
        if (any(fragatm==b(i)%center)) then
            j=j+1      !Move function in the fragment to head of list
            CO(:,j)=CO(:,i)
            b(j)=b(i)
        end if
    end do
    ifragcontri=1 !We have defined fragment
    write(*,"(' Done,',i8,' GTFs have been discarded,',i8,' GTFs reserved')") nprims-j,j
    nprims=j !Cut list at j, all functions after j seem non exist
    if (infuncsel1==-4) deallocate(exclfragatm)

    !Modification of wavefunction has finished, now reduce size of b, CO... to current nprims and nmo to avoid potential problems
    if (allocated(b)) then !Only for input file contains wavefunctions
        call resizebynmo(nmo,nprims) !Reduce size of CO, MOene, MOocc, MOtype
        allocate(b_tmp(nprims))
        b_tmp(:)=b(1:nprims)
        deallocate(b)
        allocate(b(nprims))
        b=b_tmp
        deallocate(b_tmp)
    end if
end if


!!! Every actual thing start from now on...

!!!---------------------------------------
!1!!------------------- Show system structure and view isosurface of MOs or the grid data read from cube file
if (infuncsel1==0) then !
    if (ncenter>0) write(*,*) "Nucleus list:"
    do i=1,ncenter
        write(*,"(i5,'(',a2,')',' --> Charge:',f10.6,'  x,y,z(Bohr):',3f11.6)") i,a(i)%name,a(i)%charge,a(i)%x,a(i)%y,a(i)%z
    end do
    if (allocated(CObasa).and.imodwfn==0) then !fch and occupation number hasn't been modified
        if (wfntype==0) then
            write(*,"(' Note: Orbital',i6,' is HOMO, energy:',f12.6,' a.u.',f12.6,' eV')") nint(nelec/2),MOene(nint(nelec/2)),MOene(nint(nelec/2))*au2eV
            if (nint(nelec/2)+1<=nmo) then
                write(*,"('       Orbital',i6,' is LUMO, energy:',f12.6' a.u.',f12.6,' eV')") nint(nelec/2)+1,MOene(nint(nelec/2)+1),MOene(nint(nelec/2)+1)*au2eV
                gapene=MOene(nint(nelec/2)+1)-MOene(nint(nelec/2))
                write(*,"('       LUMO/HOMO gap:',f12.6,' a.u.',f12.6,' eV',f14.6,' kJ/mol')") gapene,gapene*au2eV,gapene*au2kJ
            end if
        else if (wfntype==1) then
            write(*,"(' Note: Orbital',i6,' is HOMO of alpha spin, orbital',i6,' is HOMO of beta spin')") nint(naelec),nbasis+nint(nbelec)
            if (nbasis>=nint(naelec)+1) then
                gapenea=MOene(nint(naelec)+1)-MOene(nint(naelec))
                write(*,"('       LUMO/HOMO gap of alpha orbitals:',f12.6,' a.u.',f12.6,' eV')") gapenea,gapenea*au2eV
                gapeneb=MOene(nbasis+nint(nbelec)+1)-MOene(nbasis+nint(nbelec))
                write(*,"('       LUMO/HOMO gap of beta orbitals: ',f12.6,' a.u.',f12.6,' eV')") gapeneb,gapeneb*au2eV
            end if
        else if (wfntype==2) then
            write(*,"(' Index of SOMO orbitals:',7i6)") (i,i=nint(nbelec+1),nint(naelec))
        end if
    end if
    if (ifiletype==7.or.ifiletype==8) then !cube file
    else
    end if

!!!---------------------------------------
!1!!------------------- Output properties at a point
else if (infuncsel1==1) then
    do while(.true.)
        write(*,"(a)") " Input x,y,z, divided by space or comma, e.g. 3.3,2.0,-0.3"
        write(*,"(a)") " or input e.g. ""a5"" to use nuclear position of atom 5"
        write(*,"(a)") "    input e.g. ""o8"" to select orbital 8, whose wavefunction value will be shown"
        write(*,"(a)") "    input e.g. ""f3"" to select function 3, whose gradient and Hessian will be shown, input ""allf"" can print all available functions"        
        write(*,"(a)") "    input ""q"" to return"
        read(*,"(a)") inpstring
        inpstring=adjustl(inpstring)
        if (inpstring(1:1)=='q') then
            exit
        else if (inpstring(1:4)=='allf') then
            write(*,*) "1 Electron density (Analytical Hessian)"
            write(*,*) "3 Laplacian of electron density"
            write(*,*) "4 Value of orbital wavefunction"
            if (ELFLOL_type==0) write(*,*) "9 Electron localization function(ELF)"
            if (ELFLOL_type==1) write(*,*) "9 Electron localization function(ELF) defined by Tsirelson" 
            if (ELFLOL_type==2) write(*,*) "9 Electron localization function(ELF) defined by Lu, Tian" 
            if (ELFLOL_type==0) write(*,*) "10 Localized orbital locator(LOL)"
            if (ELFLOL_type==1) write(*,*) "10 Localized orbital locator(LOL) defined by Tsirelson" 
            if (ELFLOL_type==2) write(*,*) "10 Localized orbital locator(LOL) defined by Lu, Tian"
            write(*,*) "12 Total electrostatic potential"
            write(*,*) "100 User defined real space function"
        else if (inpstring(1:1)=='f') then    
            read(inpstring(2:),*) iprintfunc
        else if (inpstring(1:1)=='a') then
            read(inpstring(2:),*) iatm
            if (iatm>0.and.iatm<=ncenter) then
                write(*,*) "Note: Unless otherwise specified, all units are in a.u."
                write(*,"(' Atom',i6,'  X,Y,Z(Bohr):',3f14.8)") iatm,a(iatm)%x,a(iatm)%y,a(iatm)%z
                call showptprop(a(iatm)%x,a(iatm)%y,a(iatm)%z,iprintfunc,6)
            else
                write(*,*) "The index exceeds valid range"
            end if
        else if (inpstring(1:1)=='o') then
            read(inpstring(2:),*) iorbseltmp
            if (iorbseltmp>0.and.iorbseltmp<=nmo) then
                iorbsel=iorbseltmp
            else
                write(*,*) "The index exceeds valid range"
            end if
        else
            read(inpstring,*) inx,iny,inz
            write(*,*) "You inputted coordinate is in which unit?  1:Bohr  2:Angstrom"
            read(*,*) iunit
            if (iunit==2) then
                inx=inx/b2a
                iny=iny/b2a
                inz=inz/b2a
            end if
            write(*,*) "Note: Unless otherwise specified, all units are in a.u."
            call showptprop(inx,iny,inz,iprintfunc,6)
        end if
        write(*,*)
    end do

!!!---------------------------------------
!2!!------------------- Topology analysis
else if (infuncsel1==2) then
    call delvirorb(1) !Don't need virtual orbitals, delete them for faster calculation
    call topomain




!!!---------------------------------------
!!!---------------------------------------
!3!------------------- Draw property in a line
!!!---------------------------------------
!!!---------------------------------------
else if (infuncsel1==3) then
    ncustommap=0 !Clean custom operation setting that possibly defined by other modules
    if (allocated(custommapname)) deallocate(custommapname)
    if (allocated(customop)) deallocate(customop)
    write(*,*) "-10 Return to main menu"
    write(*,*) "-2 Obtain deformation property"
    write(*,*) "-1 Obtain promolecule property"
    write(*,*) "0 Set custom operation"
300 call selfunc_interface(infuncsel2)

    if (infuncsel2==0.or.infuncsel2==-1.or.infuncsel2==-2) then
        if (infuncsel2==0) call customplotsetup
        if (infuncsel2==-1) then
            ipromol=1
        else
            ipromol=0
        end if
        if (infuncsel2==-1.or.infuncsel2==-2) call setPromol
        write(*,*) "-10 Return to main menu"
        goto 300
    else if (infuncsel2==-10) then
        cycle
    else if (infuncsel2==111) then
        write(*,*) "Input indices of two atoms, e.g. 1,4, or input an atom and zero, e.g. 5,0"
        read(*,*) iatmbecke1,iatmbecke2
    end if

301    write(*,"(a,f8.4,a)") " 0 Set extension distance for mode 1, current:",aug1D," Bohr"
    write(*,*) "1 Input index of two nuclei to define a line"
    write(*,*) "2 Input coordinate of two points to define a line"
    read (*,*) infuncsel3
    
    if (infuncsel3==0) then
        write(*,*) "Input augment distance (in Bohr, e.g. 2.5)"
        read(*,*) aug1D
        goto 301
    else if (infuncsel3==1) then
        do while(.true.)
            write(*,*) "Input two number to select two nuclei (e.g. 1,3)"
            read(*,*) iselatm1,iselatm2
            if (iselatm1/=iselatm2.and.min(iselatm1,iselatm2)>=1.and.max(iselatm1,iselatm2)<=ncenter) exit
            write(*,*) "Invalid input"
        end do
        write(*,"(' Nucleus',i5,'  Charge:',f6.2,'  X,Y,Z:',3f12.6)") iselatm1,a(iselatm1)%charge,a(iselatm1)%x,a(iselatm1)%y,a(iselatm1)%z
        write(*,"(' Nucleus',i5,'  Charge:',f6.2,'  X,Y,Z:',3f12.6)") iselatm2,a(iselatm2)%charge,a(iselatm2)%x,a(iselatm2)%y,a(iselatm2)%z
        torgx=a(iselatm1)%x
        torgy=a(iselatm1)%y
        torgz=a(iselatm1)%z
        tendx=a(iselatm2)%x
        tendy=a(iselatm2)%y
        tendz=a(iselatm2)%z
        ratio=dsqrt((tendx-torgx)**2+(tendy-torgy)**2+(tendz-torgz)**2)/aug1D
        orgx1D=torgx-(tendx-torgx)/ratio
        orgy1D=torgy-(tendy-torgy)/ratio
        orgz1D=torgz-(tendz-torgz)/ratio
        endx1D=tendx+(tendx-torgx)/ratio
        endy1D=tendy+(tendy-torgy)/ratio
        endz1D=tendz+(tendz-torgz)/ratio
        totdist=dsqrt((endx1D-orgx1D)**2+(endy1D-orgy1D)**2+(endz1D-orgz1D)**2)
        atomr1=aug1D
        atomr2=totdist-aug1D
    else if (infuncsel3==2) then
        write(*,*) "Input x1,y1,z1,x2,y2,z2 to define two points (in Bohr)"
        write(*,*) "e.g. 0,0,3.2,-1,-0.26,2.8"
        read(*,*) x1,y1,z1,x2,y2,z2
        orgx1D=x1
        orgy1D=y1
        orgz1D=z1
        endx1D=x2
        endy1D=y2
        endz1D=z2
        atomr1=0D0
        atomr2=0D0
    end if
    npointcurve=num1Dpoints  !The number of data to plot
    if (infuncsel2==12) then
        npointcurve=num1Dpoints/6 !Calculate ESP is time consuming, so decrease the number of points
        write(*,*) "Please wait..."
    end if
    if (infuncsel2/=12.and.expcutoff<0) write(*,"(' Note: All exponential functions exp(x) with x<',f8.3,' is be ignored ')") expcutoff

    if (allocated(curvex)) deallocate(curvex)
    if (allocated(curvey)) deallocate(curvey)
    allocate(curvex(npointcurve))
    allocate(curvey(npointcurve))
    if (allocated(curveytmp)) deallocate(curveytmp)
    if (ncustommap/=0) allocate(curveytmp(npointcurve))
    transx=(endx1D-orgx1D)/npointcurve
    transy=(endy1D-orgy1D)/npointcurve
    transz=(endz1D-orgz1D)/npointcurve
    transr=dsqrt(transx**2+transy**2+transz**2)
    
    write(*,*)
    write(*,"(' Original point in X,Y,Z:    ',3f10.5)") orgx1D,orgy1D,orgz1D !Output grid information
    write(*,"(' End point in X,Y,Z:         ',3f10.5)") endx1D,endy1D,endz1D
    write(*,"(' Translation vector in X,Y,Z:',3f10.5,'   Norm:',f10.5)") transx,transy,transz,transr
    write(*,"(' Number of points:',i10)") npointcurve

    icustom=0
    curvey=0D0
    if (ipromol==1) goto 311
310    continue
!$OMP parallel do shared(curvex,curvey) private(i,rnowx,rnowy,rnowz) num_threads( nthreads  )
    do i=1,npointcurve  !Calculate data for line plot
        rnowx=orgx1D+i*transx
        rnowy=orgy1D+i*transy
        rnowz=orgz1D+i*transz
        curvex(i)=i*transr
        if (infuncsel2==111) then
            curvey(i)=beckewei(rnowx,rnowy,rnowz,iatmbecke1,iatmbecke2)
        else
            curvey(i)=calcfuncall(infuncsel2,rnowx,rnowy,rnowz)
        end if
    end do
!$OMP end parallel do

311    if (ncustommap/=0) then !Calculate data for custom map
        if (icustom==0) then
            curveytmp=curvey
        else if (icustom/=0) then
            if (customop(icustom)=='+') curveytmp=curveytmp+curvey
            if (customop(icustom)=='-') curveytmp=curveytmp-curvey
            if (customop(icustom)=='x'.or.customop(icustom)=='*') curveytmp=curveytmp*curvey
            if (customop(icustom)=='/') curveytmp=curveytmp/curvey
        end if
        if (icustom/=ncustommap) then
            icustom=icustom+1
            filename=custommapname(icustom)
            call dealloall
            write(*,"(' Loading:  ',a)") trim(filename)
            call readinfile(filename,1)
            !Generate temporary fragatm
            deallocate(fragatm)
            nfragatmnum=ncenter
            allocate(fragatm(nfragatmnum))
            do iatm=1,ncenter
                fragatm(iatm)=iatm
            end do
            !Input the MO index for current file. Since the MO index may be not the same as the first loaded one
            if (infuncsel2==4) then
                write(*,"(' Input the index of the orbital to be calculated for ',a,'   e.g. 3')") trim(filename)
                read(*,*) iorbsel
            end if
            goto 310
        else
            curvey=curveytmp
            call dealloall
            write(*,"(' Reloading:  ',a)") trim(firstfilename)
            call readinfile(firstfilename,1)
            !Recovery user defined fragatm from the backup
            deallocate(fragatm)
            nfragatmnum=nfragatmnumbackup
            allocate(fragatm(nfragatmnum))
            fragatm=fragatmbackup
        end if
    end if

    write(*,"(' Minimal/Maximum value:',2D16.8)") minval(curvey),maxval(curvey)
    write(*,"(' Summing up all values:',D18.8,'  Integration value:',D18.8)") sum(curvey),sum(curvey)*transr
    exty=(maxval(curvey)-minval(curvey))/10
    if (exty<maxval(abs(curvey))/100) exty=maxval(abs(curvey))/10 !Sometimes the difference between maximal and minimal values are too small, so do special treatment
    curveymin=minval(curvey)-exty
    curveymax=maxval(curvey)+exty
    steplabx=maxval(curvex)/10
    steplaby=exty
    i=-1
    
    do while(.true.)
        if (isilent==0.and.i==-1) then
            if (atomr1==atomr2) then
            else !Draw the two atom positions
            end if
        end if
        write(*,*)
        write(*,*) "-1 Show the graph again"
        write(*,*) "0 Return to main menu"
        write(*,*) "1 Save the graph to a file in current folder"
        write(*,*) "2 Export the data to line.txt in current folder"
        write(*,"(a,1PE14.6,a,1PE14.6)") " 3 Change range of Y axis, current: from",curveymin," to",curveymax
        if (icurve_vertlinex==0) write(*,*) "4 Draw a vertical line with specific X"
        if (icurve_vertlinex==1) write(*,*) "4 Delete the vertical line"
        write(*,"(' 5 Change the ratio of X and Y length, current:',f10.5)") curvexyratio
        write(*,"(' 6 Find the positions of local minimum and maximum')")
        write(*,"(' 7 Find the positions where function value equals to specified value')")
        if (ilog10y==0) write(*,*) "8 Use logarithmic scaling of Y axis"
        if (ilog10y==1) write(*,*) "8 Use linear scaling of Y axis"
        write(*,*) "9 Change the line color"
        if (ilog10y==0) write(*,"(a,f8.3,1PE14.5)") " 10 Set stepsize in X and Y axes, current:",steplabx,steplaby
        if (ilog10y==1) write(*,"(a,f8.3)") " 10 Set the stepsize in X axis, current:",steplabx
        if (ilenunit1D==1) write(*,*) "11 Change length unit of the graph to Angstrom"
        if (ilenunit1D==2) write(*,*) "11 Change length unit of the graph to Bohr"

        read(*,*) i
        if (i==-1) then
            cycle
        else if (i==0) then
            exit
        else if (i==1) then
            write(*,"(a,a,a)") " Graph have been saved to ",trim(graphformat)," file with ""DISLIN"" prefix in current directory"
        else if (i==2) then        
            open(14,file="line.txt",status="replace")
            do i=1,npointcurve  !Output result to file
                rnowx=orgx1D+i*transx
                rnowy=orgy1D+i*transy
                rnowz=orgz1D+i*transz
                curvex(i)=i*transr
                write(14,"(4f12.6,1PE18.10)") rnowx*b2a,rnowy*b2a,rnowz*b2a,curvex(i)*b2a,curvey(i)
            end do
            close (14)
            write(*,*) "Results have been saved to line.txt in current folder"
            write(*,"(a)") " Unit is Angstrom. The first three columns are actual coordinates, the fourth column &
            is X position in the curve graph, the fifth column is function value"
        else if (i==3) then
            if (ilog10y==0) write(*,*) "Input minimum and maximum value of Y axis  e.g. -0.1,2"
            if (ilog10y==1)    write(*,*) "Input minimum and maximum value of Y axis  e.g. -1,5 means from 10^-1 to 10^5"
            read(*,*) curveymin,curveymax
        else if (i==4.and.icurve_vertlinex==0) then
            write(*,*) "Input X (in Bohr) e.g. 1.78"
            read(*,*) curve_vertlinex
            icurve_vertlinex=1
        else if (i==4.and.icurve_vertlinex==1) then
            icurve_vertlinex=0
        else if (i==5) then
            write(*,*) "Input a value"
            read(*,*) curvexyratio
        else if (i==6) then
            numlocmin=0
            numlocmax=0
            do ipoint=2,npointcurve-1
                gradold=curvey(ipoint)-curvey(ipoint-1)
                gradnew=curvey(ipoint+1)-curvey(ipoint)
                if (gradold*gradnew<0D0) then
                    if (gradold>gradnew) then
                        numlocmax=numlocmax+1
                        write(*,"(' Local maximum X (Bohr):',f12.6,'  Value:',D18.8)") curvex(ipoint),curvey(ipoint)
                    else if (gradold<gradnew) then
                        numlocmin=numlocmin+1
                        write(*,"(' Local minimum X (Bohr):',f12.6,'  Value:',D18.8)") curvex(ipoint),curvey(ipoint)
                    end if
                end if
            end do
            write(*,"(' Totally found',i5,' local minimum,',i5,' local maximum')") numlocmin,numlocmax
        else if (i==7) then
            write(*,*) "Input a value"
            read(*,*) specvalue
            numfind=0
            do ipoint=1,npointcurve-1
                if ( (specvalue>curvey(ipoint).and.specvalue<curvey(ipoint+1)) .or.&
                 (specvalue<curvey(ipoint).and.specvalue>curvey(ipoint+1)) ) then
                    !Use linear interpolation to evaluate the X position
                    tmpratio=abs(specvalue-curvey(ipoint))/abs(curvey(ipoint+1)-curvey(ipoint))
                    specvaluex=curvex(ipoint)+tmpratio*transr
                    numfind=numfind+1
                    write(*,"(' #',i5,' X (Bohr):',f12.6)") numfind,specvaluex
                end if
            end do
            if (numfind==0) write(*,*) "Found nothing"
        else if (i==8) then
            if (ilog10y==1) then
                ilog10y=0
                !Recover default limit
                curveymin=minval(curvey)-(maxval(curvey)-minval(curvey))/10
                curveymax=maxval(curvey)+(maxval(curvey)-minval(curvey))/10
            else if (ilog10y==0) then
                ilog10y=1
                write(*,*) "Input minimum and maximum value of Y axis  e.g. -1,5 means from 10^-1 to 10^5"
                read(*,*) curveymin,curveymax
            end if
        else if (i==9) then
            write(*,*) "Use which color?"
            write(*,*) "1/2/3/4/5 = Red/Green/Blue/White/Black"
            write(*,*) "6/7/8/9/10 = Gray/Cyan/Yellow/Orange/Magenta"
            write(*,*) "11/12/13/14 = Crimson/Dark green/Purple/Brown"
            read(*,*) iclrcurve
        else if (i==10) then
            if (ilog10y==0) then 
                write(*,*) "Input the step size between the labels in X and Y axes, respectively"
                write(*,*) "e.g. 1.5,20"
                read(*,*) steplabx,steplaby
            else if (ilog10y==1) then
                write(*,*) "Input the step size between the labels in X axis, e.g. 1.5"
                read(*,*) steplabx
            end if
        else if (i==11) then
            if (ilenunit1D==1) then
                ilenunit1D=2 !Angstrom
            else if (ilenunit1D==2) then
                ilenunit1D=1 !Bohr
            end if
        end if
    end do



!!!----------------------------------------
!!!----------------------------------------
!4!----------------------- Draw plane graph
!!!----------------------------------------
!!!----------------------------------------
else if (infuncsel1==4) then
    ncustommap=0 !Clean custom operation setting that possibly defined by other modules
    if (allocated(custommapname)) deallocate(custommapname)
    if (allocated(customop)) deallocate(customop)
    if (allocated(tmparrint)) deallocate(tmparrint)
    write(*,*) "-10 Return to main menu"
    write(*,*) "-2 Obtain of deformation property"
    write(*,*) "-1 Obtain of promolecule property"
    write(*,*) "0 Set custom operation"
400    call funclist
    read(*,*) infuncsel2
    if (infuncsel2==4) then !4=MO value, 111=becke weighting function
        write(*,"(a,i10)") " Input the orbital index that you want to plot, should between 1 and",nmo
        write(*,"(a)") " If you want to plot contour map for two orbitals at the same time, input two indices and separate them by comma, e.g. 8,10"
        read(*,"(a)") c200tmp
        if (index(c200tmp,',')/=0) then !Inputted two orbitals
            read(c200tmp,*) iorbsel,iorbsel2
        else
            read(c200tmp,*) iorbsel
        end if
    else if (infuncsel2==111) then !Calculate Becke weighting function
        write(*,*) "Input indices of two atoms to calculate Becke overlap weight, e.g. 1,4"
        write(*,*) "or input index of an atom and zero to calculate Becke atomic weight, e.g. 5,0"
        read(*,*) iatmbecke1,iatmbecke2
    else if (infuncsel2==112) then !Calculate Hirshfeld weighting function
        write(*,*) "Input index of the atoms you are interested in, e.g. 2,3,7-10"
        read(*,"(a)") c2000tmp
        call str2arr(c2000tmp,ntmp)
        allocate(tmparrint(ntmp))
        call str2arr(c2000tmp,ntmp,tmparrint)
        write(*,"(a)") " How to generate the atomic densities that used in the calculation of Hirshfeld weight?"
        write(*,*) "1 Based on atomic .wfn files"
        write(*,*) "2 Based on built-in atomic densities (see Appendix 3 of the manual for detail)"
        read(*,*) iHirshdenstype
        if (iHirshdenstype==1) call setpromol
    else if (infuncsel2==500.or.infuncsel2==510.or.infuncsel2==511.or.infuncsel2==512) then !Calculate rho(A)*ln[rho(A)/rho0(A)], or rho(A), or rho0(A)
        call setpromol
        allocate(tmparrint(1))
        write(*,*) "Input index of the atom you are interested in, e.g. 4"
        read(*,*) tmparrint
    else if (infuncsel2==501.or.infuncsel2==502.or.infuncsel2==503) then !sum[rhoA*ln(rhoA/rho0A), sum[(rhoA-rho0A)/rhoA], difference between relative entropy and deformation density
        call setpromol
    end if
    
    inucespplot=0
    if (infuncsel2==8) inucespplot=1 !For deal with plotting nucesp property, special treatment is needed
    imarkrefpos=0
    if (infuncsel2==17.or.infuncsel2==19.or.(infuncsel2==100.and.iuserfunc==24)) imarkrefpos=1 !Only for correlation hole/factor, source function, linear response kernel... reference point will be marked on contour map

    if (infuncsel2==0.or.infuncsel2==-1.or.infuncsel2==-2) then
        if (infuncsel2==0) call customplotsetup
        if (infuncsel2==-1) then
            ipromol=1
        else
            ipromol=0
        end if
        if (infuncsel2==-1.or.infuncsel2==-2) call setPromol
        write(*,*) "-10 Return to main menu"
        goto 400
    else if (infuncsel2==-10) then
        cycle
    end if
    
    if (iorbsel2/=0) then
        idrawtype=2 !Only contour line map is available when another orbital is needed to be plotted together
    else if (iorbsel2==0) then
        write(*,*) " -10 Return to main menu"
        write(*,*) "Draw which type of map?"
        write(*,*) "1 Color-filled map"
        write(*,*) "2 Contour line map"
        write(*,*) "3 Relief map"
        write(*,*) "4 Shaded surface map"
        write(*,*) "5 Shaded surface map with projection"
        write(*,*) "6 Gradient lines map with/without contour lines"
        write(*,*) "7 Vector field map with/without contour lines"
        read(*,*) idrawtype
    end if

    if (idrawtype==1.or.idrawtype==2.or.idrawtype==6.or.idrawtype==7) then  !Initialize contour line setting
        if (idrawtype==2.or.idrawtype==6) idrawcontour=1
        if (idrawtype==1.or.idrawtype==7) idrawcontour=0
        if (idrawtype/=1) iatom_on_contour=1
        ilabel_on_contour=0
        ctrval=0.0D0
        if (infuncsel2==9.or.infuncsel2==10) then  !A special contour setting suitable for ELF and LOL
            lastctrval=21
            do i=1,21
                ctrval(i)=(i-1)*0.05D0
            end do
        else  !General contour setting for other real space functions
            lastctrval=62
            ctrval(1)=1D-3
            do i=0,9 !set the value of contour line
                ctrval(3*i+2)=2*1D-3*10**i
                ctrval(3*i+3)=4*1D-3*10**i
                ctrval(3*i+4)=8*1D-3*10**i
            end do
            ctrval(32:62)=-ctrval(1:31)
        end if
    else if (idrawtype==-10) then
        cycle
    end if
    
    write(*,*) " -10 Return to main menu"
    write(*,*) "How many grids in the two dimensions respectively?"
    if (idrawtype==1.or.idrawtype==2.or.idrawtype==6) then
        if (infuncsel2==12.or.(infuncsel2==100.and.iuserfunc==60).or.iuserfunc==39.or.iuserfunc==101.or.iuserfunc==102) then
            write(*,*) "(100,100 is recommended)" !Because calculating ESP is very time consuming, so use lower grid
        else
            write(*,*) "(200,200 is recommended)"
        end if
    else if (idrawtype==3.or.idrawtype==4.or.idrawtype==5) then
        write(*,*) "(100,100 is recommended)"
    else if (idrawtype==7) then
        write(*,*) "(80,80 is recommended)"
    end if
    write(*,*) "Hint: You can press ENTER button directly to use recommended value"
    read(*,"(a)") c200tmp
    if (c200tmp==' ') then !Press enter directly
        if (idrawtype==1.or.idrawtype==2.or.idrawtype==6) then
            if (infuncsel2==12.or.(infuncsel2==100.and.iuserfunc==60).or.iuserfunc==39.or.iuserfunc==101.or.iuserfunc==102) then
                ngridnum1=100
                ngridnum2=100
            else
                ngridnum1=200
                ngridnum2=200
            end if
        else if (idrawtype==3.or.idrawtype==4.or.idrawtype==5) then
            ngridnum1=100
            ngridnum2=100
        else if (idrawtype==7) then
            ngridnum1=80
            ngridnum2=80
        end if
    else if (c200tmp=='-10') then
        cycle
    else
        read(c200tmp,*) ngridnum1,ngridnum2
    end if

    if (allocated(planemat)) deallocate(planemat)
    if (allocated(planemattmp)) deallocate(planemattmp)
    allocate(planemat(ngridnum1,ngridnum2),planemattmp(ngridnum1,ngridnum2)) !planemattmp is used in many cases below
    if (ncustommap/=0) then
        if (allocated(planemat_cust)) deallocate(planemat_cust)
        if (allocated(d1addtmp)) deallocate(d1addtmp)
        if (allocated(d1mintmp)) deallocate(d1mintmp)
        if (allocated(d2addtmp)) deallocate(d2addtmp)
        if (allocated(d2mintmp)) deallocate(d2mintmp)
        allocate(planemat_cust(ngridnum1,ngridnum2))
        allocate(d1addtmp(ngridnum1,ngridnum2))
        allocate(d1mintmp(ngridnum1,ngridnum2))
        allocate(d2addtmp(ngridnum1,ngridnum2))
        allocate(d2mintmp(ngridnum1,ngridnum2))
    end if
    if (idrawtype==6.or.idrawtype==7) then !Draw gradient lines
        if (allocated(d1add)) deallocate(d1add)
        if (allocated(d1min)) deallocate(d1min)
        if (allocated(d2add)) deallocate(d2add)
        if (allocated(d2min)) deallocate(d2min)
        if (allocated(gradd1)) deallocate(gradd1)
        if (allocated(gradd2)) deallocate(gradd2)
        allocate(d1add(ngridnum1,ngridnum2))
        allocate(d1min(ngridnum1,ngridnum2))
        allocate(d2add(ngridnum1,ngridnum2))
        allocate(d2min(ngridnum1,ngridnum2))
        allocate(gradd1(ngridnum1,ngridnum2))
        allocate(gradd2(ngridnum1,ngridnum2))
    end if

    write(*,*) " -10 Return to main menu" 
401 write(*,*) "Define the plane to be plotted"
    write(*,*) "1:XY 2:XZ 3:YZ 4:Define by three atoms 5:Define by three points"
    write(*,*) "6:Input origin and translation vector (For expert user)"
    write(*,*) "7:Parallel to a bond and meantime normal to a plane defined by three atoms"
    write(*,"(a,f8.4,a)") " 0:Set extension distance for plane type 1~5, current:",aug2D," Bohr"
    read(*,*) plesel
    aug2D2=aug2D !If don't draw gradient line map, needn't make the augment in the two dimension different
    orgx2D=0D0
    orgy2D=0D0
    orgz2D=0D0
    v1x=0D0
    v1y=0D0
    v1z=0D0
    v2x=0D0
    v2y=0D0
    v2z=0D0
    if (plesel==-10) then
        cycle
    else if (plesel==0) then
        write(*,*) "Input extension distance in Bohr, e.g. 4.5"
        read(*,*) aug2D
        goto 401
    else if (plesel==-1) then
        write(*,*) "Input rotation angle of the plotting plane in degree, e.g. 30.5"
        read(*,*) rot2Dple
        goto 401
    else if (plesel==1) then
        write(*,*) "Input Z value in Bohr"
        write(*,*) "Note: If the unit is in Angstrom, add ""a"" suffix, e.g. -1.6a"
        read(*,*) c200tmp
        if (index(c200tmp,'a')/=0) then
            read(c200tmp(1:len_trim(c200tmp)-1),*) orgz2D
            orgz2D=orgz2D/b2a
        else
            read(c200tmp,*) orgz2D
        end if
        if (ncenter==0) then
            orgx2D=orgx
            orgy2D=orgy
            v1x=(endx-orgx)/(ngridnum1-1)
            v2y=(endy-orgy)/(ngridnum2-1)
        else
            !Adjust aug2D/aug2D2 to make the length in two direction same
            !Because of bug in stream() of dislin 9.5D, if the two length unequal, some region won't be shown
            if (idrawtype==6.or.idrawtype==7) then  
                sup=(maxval(a%x)-minval(a%x))-(maxval(a%y)-minval(a%y))
                if (sup>0) aug2D2=aug2D+sup/2
                if (sup<0) aug2D=aug2D-sup/2
            end if
            orgx2D=minval(a%x)-aug2D
            orgy2D=minval(a%y)-aug2D2
            v1x=(maxval(a%x)+aug2D-orgx2D)/ngridnum1
            v2y=(maxval(a%y)+aug2D2-orgy2D)/ngridnum2
        end if
    else if (plesel==2) then
        write(*,*) "Input Y value in Bohr"
        write(*,*) "Note: If the unit is in Angstrom, add ""a"" suffix, e.g. -1.6a"
        read(*,*) c200tmp
        if (index(c200tmp,'a')/=0) then
            read(c200tmp(1:len_trim(c200tmp)-1),*) orgy2D
            orgy2D=orgy2D/b2a
        else
            read(c200tmp,*) orgy2D
        end if
        if (ncenter==0) then
            orgx2D=orgx
            orgz2D=orgz
            v1x=(endx-orgx)/(ngridnum1-1)
            v2z=(endz-orgz)/(ngridnum2-1)
        else
            if (idrawtype==6.or.idrawtype==7) then !adjust aug2D/aug2D2
                sup=(maxval(a%x)-minval(a%x))-(maxval(a%z)-minval(a%z))
                if (sup>0) aug2D2=aug2D+sup/2
                if (sup<0) aug2D=aug2D-sup/2
            end if
            orgx2D=minval(a%x)-aug2D
            orgz2D=minval(a%z)-aug2D2
            v1x=(maxval(a%x)+aug2D-orgx2D)/ngridnum1
            v2z=(maxval(a%z)+aug2D2-orgz2D)/ngridnum2
        end if
    else if (plesel==3) then
        write(*,*) "Input X value in Bohr"
        write(*,*) "Note: If the unit is in Angstrom, add ""a"" suffix, e.g. -1.6a"
        read(*,*) c200tmp
        if (index(c200tmp,'a')/=0) then
            read(c200tmp(1:len_trim(c200tmp)-1),*) orgx2D
            orgx2D=orgx2D/b2a
        else
            read(c200tmp,*) orgx2D
        end if
        if (ncenter==0) then
            orgy2D=orgy
            orgz2D=orgz
            v1y=(endy-orgy)/(ngridnum1-1)
            v2z=(endz-orgz)/(ngridnum2-1)
        else
            if (idrawtype==6.or.idrawtype==7) then !adjust aug2D/aug2D2
                sup=(maxval(a%y)-minval(a%y))-(maxval(a%z)-minval(a%z))
                if (sup>0) aug2D2=aug2D+sup/2
                if (sup<0) aug2D=aug2D-sup/2
            end if
            orgy2D=minval(a%y)-aug2D
            orgz2D=minval(a%z)-aug2D2
            v1y=(maxval(a%y)+aug2D-orgy2D)/ngridnum1
            v2z=(maxval(a%z)+aug2D2-orgz2D)/ngridnum2
        end if
    else if (plesel==4.or.plesel==5) then
410        if (plesel==4) then
            write(*,*) "Input the number of three atoms (e.g. 3,6,7)"
            read(*,*) i1,i2,i3
            if (i1==i2.or.i1==i3.or.i2==i3.or.min(i1,i2,i3)<1.or.max(i1,i2,i3)>ncenter) then
                if (min(i1,i2,i3)<1.or.max(i1,i2,i3)>ncenter) write(*,*) "Atom indices are out of valid range, please input again"
                if (i1==i2.or.i1==i3.or.i2==i3) write(*,*) "Atom indices are duplicated, please input again"
                goto 410
            end if
            a1x=a(i1)%x
            a1y=a(i1)%y
            a1z=a(i1)%z
            a2x=a(i2)%x
            a2y=a(i2)%y
            a2z=a(i2)%z
            a3x=a(i3)%x
            a3y=a(i3)%y
            a3z=a(i3)%z
        else if (plesel==5) then
            write(*,*) "Input x,y,z of point 1 (in Bohr, e.g. 1.2,1.3,0.0)"
            read(*,*) a1x,a1y,a1z
            write(*,*) "Input x,y,z of point 2 (in Bohr)"
            read(*,*) a2x,a2y,a2z
            write(*,*) "Input x,y,z of point 3 (in Bohr)"
            read(*,*) a3x,a3y,a3z
        end if
        v1x=a1x-a2x
        v1y=a1y-a2y
        v1z=a1z-a2z
        v2x=a3x-a2x
        v2y=a3y-a2y
        v2z=a3z-a2z
        rnorm1=dsqrt(v1x**2+v1y**2+v1z**2) !Norm of vector 1
        rnorm2=dsqrt(v2x**2+v2y**2+v2z**2)
        if (abs(v1x*v2x+v1y*v2y+v1z*v2z)/(rnorm1*rnorm2)>0.999) then
            if (plesel==4) write(*,*) "The three atoms should not lie in the same line!"
            if (plesel==5) write(*,*) "The three points should not lie in the same line!"
            goto 410
        end if
        rangle=acos( abs(v1x*v2x+v1y*v2y+v1z*v2z)/(rnorm1*rnorm2) )
        if (idrawtype==6.or.idrawtype==7) then !adjust aug2D/aug2D2
            sup=rnorm1-rnorm2*cos(pi/2-rangle)
            if (sup>0) aug2D2=aug2D+sup/2
            if (sup<0) aug2D=aug2D-sup/2
        end if
        dist1=rnorm1+2*aug2D !Total length in direction 1
        dist2=rnorm2*cos(pi/2-rangle)+2*aug2D2 !Total length in direction 2
!        write(*,*) "angle",acos(cos(pi/2-rangle))/pi*180
        d1=dist1/ngridnum1  !Transitional step length in direction 1
        d2=dist2/ngridnum2
        v1x=v1x*d1/rnorm1  !Make the norm of v1 equal to expected step lengh (d1)
        v1y=v1y*d1/rnorm1
        v1z=v1z*d1/rnorm1
        schmit=(v1x*v2x+v1y*v2y+v1z*v2z)/(v1x**2+v1y**2+v1z**2) !Use schmit method to make v2 ortho to v1
        v2x=v2x-schmit*v1x
        v2y=v2y-schmit*v1y
        v2z=v2z-schmit*v1z
        rnorm2=dsqrt(v2x**2+v2y**2+v2z**2)
        v2x=v2x*d2/rnorm2   !Make the norm of v2 equal to expected step lengh (d2)
        v2y=v2y*d2/rnorm2
        v2z=v2z*d2/rnorm2
!         write(*,*) "test ortho",v1x*v2x+v1y*v2y+v1z*v2z
        orgx2D=a2x-aug2D/d1*v1x-aug2D2/d2*v2x  !aug2D/d1*v1x=aug2D*(v1x/d1), v1x/d1 correspond the x component of unit vector in v1x direction
        orgy2D=a2y-aug2D/d1*v1y-aug2D2/d2*v2y
        orgz2D=a2z-aug2D/d1*v1z-aug2D2/d2*v2z
    else if (plesel==6) then
        write(*,*) "Input origin of x,y,z in Bohr, e.g. 3.5,-1,0.2"
        read(*,*) orgx2D,orgy2D,orgz2D
        write(*,*) "Input x,y,z of transitional vector 1 in Bohr, e.g. 0.08,0.03,0"
        read(*,*) v1x,v1y,v1z
        write(*,*) "Input x,y,z of transitional vector 2 in Bohr, should be orthogonal to vector 1"
        read(*,*) v2x,v2y,v2z
        d1=dsqrt(v1x**2+v1y**2+v1z**2)
        d2=dsqrt(v2x**2+v2y**2+v2z**2)
        dist1=d1*ngridnum1
        dist2=d2*ngridnum2
        a1x=orgx2D !Although a1x...a3z is no use for generate data, but these are critical for plotting atom label (subroutine drawplane)
        a1y=orgy2D
        a1z=orgz2D
        a2x=orgx2D+ngridnum1*v1x
        a2y=orgy2D+ngridnum1*v1y
        a2z=orgz2D+ngridnum1*v1z
        a3x=orgx2D+ngridnum2*v2x
        a3y=orgy2D+ngridnum2*v2y
        a3z=orgz2D+ngridnum2*v2z
    else if (plesel==7) then
        write(*,*) "Input two atoms to define the bond, e.g. 4,5"
        read(*,*) iatm1,iatm2
        write(*,*) "Input three atoms to define a plane, e.g. 1,4,7"
        read(*,*) jatm1,jatm2,jatm3
        write(*,*) "Input length of X-axis in Bohr e.g. 10"
        read(*,*) dist1
        write(*,*) "Input length of Y-axis in Bohr e.g. 8"
        read(*,*) dist2
        v1x=a(iatm2)%x-a(iatm1)%x
        v1y=a(iatm2)%y-a(iatm1)%y
        v1z=a(iatm2)%z-a(iatm1)%z
        call pointABCD(a(jatm1)%x,a(jatm1)%y,a(jatm1)%z,a(jatm2)%x,a(jatm2)%y,a(jatm2)%z,a(jatm3)%x,a(jatm3)%y,a(jatm3)%z,v2x,v2y,v2z,tmpD)
        schmit=(v1x*v2x+v1y*v2y+v1z*v2z)/(v1x**2+v1y**2+v1z**2) !Use schmit method to make v2 ortho to v1
        v2x=v2x-schmit*v1x
        v2y=v2y-schmit*v1y
        v2z=v2z-schmit*v1z
        d1=dist1/(ngridnum1-1)
        d2=dist2/(ngridnum2-1)
        rnorm1=dsqrt(v1x**2+v1y**2+v1z**2)
        v1x=v1x/rnorm1*d1
        v1y=v1y/rnorm1*d1
        v1z=v1z/rnorm1*d1
        rnorm2=dsqrt(v2x**2+v2y**2+v2z**2)
        v2x=v2x/rnorm2*d2
        v2y=v2y/rnorm2*d2
        v2z=v2z/rnorm2*d2
        orgx2D=(a(iatm1)%x+a(iatm2)%x)/2 - v1x*ngridnum1/2 - v2x*ngridnum2/2
        orgy2D=(a(iatm1)%y+a(iatm2)%y)/2 - v1y*ngridnum1/2 - v2y*ngridnum2/2
        orgz2D=(a(iatm1)%z+a(iatm2)%z)/2 - v1z*ngridnum1/2 - v2z*ngridnum2/2
        a1x=orgx2D !Although a1x...a3z is no use for generate data, but these are critical for plotting atom label (subroutine drawplane)
        a1y=orgy2D
        a1z=orgz2D
        a2x=orgx2D+ngridnum1*v1x
        a2y=orgy2D+ngridnum1*v1y
        a2z=orgz2D+ngridnum1*v1z
        a3x=orgx2D+ngridnum2*v2x
        a3y=orgy2D+ngridnum2*v2y
        a3z=orgz2D+ngridnum2*v2z
    end if
    
    write(*,*)
    write(*,"(' X/Y/Z of origin of the plane:',3f10.5,' Bohr')") orgx2D,orgy2D,orgz2D
    endx2D=orgx2D+v1x*(ngridnum1-1)+v2x*(ngridnum2-1)
    endy2D=orgy2D+v1y*(ngridnum1-1)+v2y*(ngridnum2-1)
    endz2D=orgz2D+v1z*(ngridnum1-1)+v2z*(ngridnum2-1)
    write(*,"(' X/Y/Z of end of the plane:   ',3f10.5,' Bohr')") endx2D,endy2D,endz2D
    write(*,"(' X/Y/Z of translation vector 1:',3f10.5,' Bohr')") v1x,v1y,v1z
    write(*,"(' X/Y/Z of translation vector 2:',3f10.5,' Bohr')") v2x,v2y,v2z
    write(*,*)
    rnorm1=dsqrt(v1x**2+v1y**2+v1z**2) !The final length of vector 1
    rnorm2=dsqrt(v2x**2+v2y**2+v2z**2) !The final length of vector 2
    diff=1D-5
    diffv1x=diff*v1x/rnorm1 !The infinitesimal in each direction for gradient plot
    diffv1y=diff*v1y/rnorm1
    diffv1z=diff*v1z/rnorm1
    diffv2x=diff*v2x/rnorm2
    diffv2y=diff*v2y/rnorm2
    diffv2z=diff*v2z/rnorm2
    
    !Don't directly use Multiwfn to calculate the plane data, but load from external file
    if (iplaneextdata==1) then
        open(10,file="planept.txt",status="replace")
        open(11,file="cubegenpt.txt",status="replace")
        write(10,*) ngridnum1*ngridnum2
        do ipt=1,ngridnum1
            do jpt=1,ngridnum2
                rnowx=orgx2D+(ipt-1)*v1x+(jpt-1)*v2x
                rnowy=orgy2D+(ipt-1)*v1y+(jpt-1)*v2y
                rnowz=orgz2D+(ipt-1)*v1z+(jpt-1)*v2z
                write(10,"(3f16.8)") rnowx,rnowy,rnowz
                write(11,"(3f16.8)") rnowx*b2a,rnowy*b2a,rnowz*b2a
            end do
        end do
        close(10)
        close(11)
        write(*,"(a)") " The coordinate of all points needed to be calculated have been outputted to plane.txt in current folder, the unit is in Bohr"
        write(*,"(a)") " cubegenpt.txt is also outputted, which is similar to plane.txt, but the unit is in Angstrom, &
        and there is no first line (the number of points). It can be directly utilize by cubegen"
        write(*,"(a)") " For example ""cubegen 0 potential CNT.fch result.cub -5 h < cubegenpt.txt"""
        write(*,*)
        write(*,"(a)") " Now input the path of the file containing function values, e.g. c:\t.txt, whose format should be identical to plane.txt, but with function value in the fourth column"
        write(*,"(a)") " Note: If the suffix is .cub, then the file will be recognized and loaded as output file of cubegen"
        do while(.true.)
            read(*,"(a)") c200tmp
            inquire(file=c200tmp,exist=alive)
            if (alive) exit
            write(*,*) "File cannot be found, input again"
        end do
        open(10,file=c200tmp,status="old")
        if (index(c200tmp,".cub")/=0) then !cubegen output
            do iskip=1,6+ncenter
                read(10,*)
            end do
        else
            read(10,*)
        end if
        do ipt=1,ngridnum1
            do jpt=1,ngridnum2
                read(10,*) rnouse,rnouse,rnouse,planemat(ipt,jpt)
            end do
        end do
        close(10)
        goto 430
    end if

    !!! Start calculation of plane data now!
    if (infuncsel2/=4) call delvirorb(1) !Delete high-lying virtual orbitals for faster calculation, but don't do this for analyzing MO
    write(*,*) "Please wait..."    
    if (infuncsel2/=12.and.expcutoff<0) write(*,"(' Note: All exponential functions exp(x) with x<',f8.3,' will be ignored ')") expcutoff
    CALL CPU_TIME(time_begin)
    call walltime(walltime1)
    icustom=0
    planemat=0D0 !For promolecular property, first clean up planemat, if not planemat may not clean
    if (ipromol==1) goto 421 ! To obtain promolecule property, pass the first loaded molecule
    !Note: If the task refers to plotting ESP map, the function of drawing gradient line will be disabled
    
420    if (infuncsel2==12.and.idrawtype/=6.and.idrawtype/=7) then
        call planeesp(max(ngridnum1,ngridnum2)) !Special treatment to calculate ESP
    else if (infuncsel2==112) then !Calculate atomic Hirshfeld weighting function
        call genhirshplanewei(tmparrint,size(tmparrint),iHirshdenstype)
        ncustommap=0
    else if (infuncsel2==500.or.infuncsel2==510.or.infuncsel2==511.or.infuncsel2==512) then
    !500: Calculate rho(A)*ln[rho(A)/rho0(A)], 510: Calculate rho(A), 511: Calculate rho0(A), 512: other
        call genhirshplanewei(tmparrint,size(tmparrint),1)
        ncustommap=0
        do i=1,ngridnum1 !Now planemat is Hirshfeld weight of iatmentropy, and planemattmp is its density in free-state
            do j=1,ngridnum2
                rnowx=orgx2D+(i-1)*v1x+(j-1)*v2x
                rnowy=orgy2D+(i-1)*v1y+(j-1)*v2y
                rnowz=orgz2D+(i-1)*v1z+(j-1)*v2z
                rhoA=planemat(i,j)*fdens(rnowx,rnowy,rnowz)
                if (infuncsel2==500) planemat(i,j)=rhoA*log(rhoA/planemattmp(i,j))
                if (infuncsel2==510) planemat(i,j)=rhoA
                if (infuncsel2==511) planemat(i,j)=planemattmp(i,j)
                if (infuncsel2==512) planemat(i,j)=log(rhoA/planemattmp(i,j))
            end do
        end do
    else if (infuncsel2==501) then !Calculate sum{rho(A)*ln[rho(A)/rho0(A)]}
        call genentroplane(1)
        ncustommap=0
    else if (infuncsel2==502) then !Calculate sum(x), where x=[rho(A)-rho0(A)]/rho(A)
        call genentroplane(2)
        ncustommap=0
    else if (infuncsel2==503) then !Calculate difference between total relative Shannon entropy and deformation density 
        call genentroplane(3)
        ncustommap=0
    else
nthreads=getNThreads()
!$OMP PARALLEL DO private(i,j,rnowx,rnowy,rnowz) shared(planemat,d1add,d1min,d2add,d2min) schedule(dynamic) NUM_THREADS(nthreads)
        do i=1,ngridnum1
            do j=1,ngridnum2
                rnowx=orgx2D+(i-1)*v1x+(j-1)*v2x
                rnowy=orgy2D+(i-1)*v1y+(j-1)*v2y
                rnowz=orgz2D+(i-1)*v1z+(j-1)*v2z
                if (infuncsel2==111) then
                    planemat(i,j)=beckewei(rnowx,rnowy,rnowz,iatmbecke1,iatmbecke2)
                else
                    planemat(i,j)=calcfuncall(infuncsel2,rnowx,rnowy,rnowz)
                    if (infuncsel2==4.and.iorbsel2/=0) planemattmp(i,j)=fmo(rnowx,rnowy,rnowz,iorbsel2) !Calculate another orbital together
                end if
                if (idrawtype==6.or.idrawtype==7) then !Generate two vector to plot vector field line graph
                    d1add(i,j)=calcfuncall(infuncsel2,rnowx+diffv1x,rnowy+diffv1y,rnowz+diffv1z)
                    d1min(i,j)=calcfuncall(infuncsel2,rnowx-diffv1x,rnowy-diffv1y,rnowz-diffv1z)
                    d2add(i,j)=calcfuncall(infuncsel2,rnowx+diffv2x,rnowy+diffv2y,rnowz+diffv2z)
                    d2min(i,j)=calcfuncall(infuncsel2,rnowx-diffv2x,rnowy-diffv2y,rnowz-diffv2z)
                end if
            end do
        end do
!$OMP END PARALLEL DO
    end if

421    if (ncustommap/=0) then !Calculate data for custom map
        if (icustom==0) then !The first time
            planemat_cust=planemat
            if (idrawtype==6.or.idrawtype==7) then
                d1addtmp=d1add
                d1mintmp=d1min
                d2addtmp=d2add
                d2mintmp=d2min
            end if
        else if (icustom/=0) then
            if (customop(icustom)=='+') then
                planemat_cust=planemat_cust+planemat
                if (idrawtype==6.or.idrawtype==7) then
                    d1addtmp=d1addtmp+d1add
                    d1mintmp=d1mintmp+d1min
                    d2addtmp=d2addtmp+d2add
                    d2mintmp=d2mintmp+d2min
                end if
            else if (customop(icustom)=='-') then
                planemat_cust=planemat_cust-planemat
                if (idrawtype==6.or.idrawtype==7) then
                    d1addtmp=d1addtmp-d1add
                    d1mintmp=d1mintmp-d1min
                    d2addtmp=d2addtmp-d2add
                    d2mintmp=d2mintmp-d2min
                end if
            else if (customop(icustom)=='x'.or.customop(icustom)=='*') then
                planemat_cust=planemat_cust*planemat
                if (idrawtype==6.or.idrawtype==7) then
                    d1addtmp=d1addtmp*d1add
                    d1mintmp=d1mintmp*d1min
                    d2addtmp=d2addtmp*d2add
                    d2mintmp=d2mintmp*d2min
                end if
            else if (customop(icustom)=='/') then
                planemat_cust=planemat_cust/planemat
                if (idrawtype==6.or.idrawtype==7) then
                    d1addtmp=d1addtmp/d1add
                    d1mintmp=d1mintmp/d1min
                    d2addtmp=d2addtmp/d2add
                    d2mintmp=d2mintmp/d2min
                end if
            end if
        end if
        if (icustom/=ncustommap) then !Not the final time
            icustom=icustom+1
            filename=custommapname(icustom)
            call dealloall
            write(*,"(' Loading:  ',a)") trim(filename)
            call readinfile(filename,1)
            if (infuncsel2/=4) call delvirorb(0)
            !Generate temporary fragatm
            deallocate(fragatm)
            nfragatmnum=ncenter
            allocate(fragatm(nfragatmnum))
            do iatm=1,ncenter
                fragatm(iatm)=iatm
            end do
            !Input the MO index for current file. Since the MO index may be not the same as the first loaded one
            if (infuncsel2==4) then
                write(*,"(' Input the index of the orbital to be calculated for ',a,'   e.g. 3')") trim(filename)
                read(*,*) iorbsel
            end if
            goto 420
        else !The final time, reload the first loaded system
            planemat=planemat_cust
            if (idrawtype==6.or.idrawtype==7) then
                d1add=d1addtmp
                d1min=d1mintmp
                d2add=d2addtmp
                d2min=d2mintmp
            end if
            call dealloall
            write(*,"(' Reloading:  ',a)") trim(firstfilename)
            call readinfile(firstfilename,1)
            !Recovery user defined fragatm from the backup
            deallocate(fragatm)
            nfragatmnum=nfragatmnumbackup
            allocate(fragatm(nfragatmnum))
            fragatm=fragatmbackup
        end if
    end if

    if (idrawtype==6.or.idrawtype==7) then !Finish the finite differential to yield gradient
        gradd1=(d1add-d1min)/2/diff
        gradd2=(d2add-d2min)/2/diff
    end if
    CALL CPU_TIME(time_end)
    call walltime(walltime2)

!     open(14,file="data(h2o-O1).txt",status="old") !read data from external file, normal user don't play with this
!     do i=1,ngridnum1
!         do j=1,ngridnum2
!             read(14,*) rabbish,rabbish,rabbish,planemat(i,j)
!         end do
!     end do
!     close(14)

430    if (isys==1) write(*,"(' Calculation took up CPU time',f12.2,'s, wall clock time',i8,'s')") time_end-time_begin,walltime2-walltime1
    write(*,*) "The minimum of data:",minval(planemat)
    write(*,*) "The maximum of data:",maxval(planemat)
    !! Set default lower and upper of Z for plot
    surcolorzmin=-3
    surcolorzmax=3
    if (infuncsel2==1) then
        drawlowlim=0.0D0
        drawuplim=0.65D0
    else if (infuncsel2==2) then
        drawlowlim=0.0D0
        drawuplim=0.65D0
    else if (infuncsel2==3) then 
        drawlowlim=-8.0D0
        drawuplim=15.0D0
    else if (infuncsel2==4) then
        drawlowlim=-0.8D0
        drawuplim=0.8D0
    else if (infuncsel2==5) then
        drawlowlim=-0.1D0
        drawuplim=0.1D0
    else if (infuncsel2==8.and.ifiletype/=4) then !Nuclear ESP
        drawlowlim=0.0D0
        drawuplim=50D0
        surcolorzmin=0
        surcolorzmax=50
    else if (infuncsel2==8.and.ifiletype==4) then !Atomic charge ESP
        drawlowlim=-0.4D0
        drawuplim=0.4D0
    else if (infuncsel2==9) then
        drawlowlim=0.0D0
        drawuplim=1.0D0
    else if (infuncsel2==10) then
        drawlowlim=0.0D0
        drawuplim=0.8D0
    else if (infuncsel2==11) then
        drawlowlim=0.0D0
        drawuplim=0.1D0
    else if (infuncsel2==12) then
        drawlowlim=-0.1D0
        drawuplim=0.1D0
    else if (infuncsel2==13.or.infuncsel2==14) then
        drawlowlim=0D0
        drawuplim=1D0
    else if (infuncsel2==15.or.infuncsel2==16) then
        drawlowlim=-0.65D0
        drawuplim=0.65D0
    else if (infuncsel2==17) then
        drawlowlim=-0.5D0
        drawuplim=0.1D0
    else if (infuncsel2==18) then
        drawlowlim=0D0
        drawuplim=2D0
    else if (infuncsel2==111.or.infuncsel2==112) then !Becke/Hirshfeld weight
        drawlowlim=0D0
        drawuplim=1D0
    else if (infuncsel2==100.and.iuserfunc==20) then !DORI
        drawlowlim=0D0
        drawuplim=1D0
    else ! Include infuncsel2==100
        drawlowlim=0.0D0
        drawuplim=5.0D0
    end if
    !Set up range of X and Y axes
    if (plesel==1) then
        axlow1=orgx2D
        axhigh1=endx2D
        axlow2=orgy2D
        axhigh2=endy2D
    else if (plesel==2) then
        axlow1=orgx2D
        axhigh1=endx2D
        axlow2=orgz2D
        axhigh2=endz2D
    else if (plesel==3) then
        axlow1=orgy2D
        axhigh1=endy2D
        axlow2=orgz2D
        axhigh2=endz2D
    else if (plesel==4.or.plesel==5.or.plesel==6.or.plesel==7) then
        axlow1=0D0
        axhigh1=dist1-d1
        axlow2=0D0
        axhigh2=dist2-d2
    end if
    !Step size between labels
    planestpx=(axhigh1-axlow1)/7
    planestpy=(axhigh2-axlow2)/7
    planestpz=(drawuplim-drawlowlim)/10

    XVU=150.0D0 !Reinitialize view
    YVU=30.0D0
    ZVU=7.0D0 !More suitable than 6.0D0 for drawing plane
    i=-1

    idrawintbasple=0 !Refresh (3,-1) information
    nple3n1path=0
    cp2ple3n1path=0
    iatmcontri=0 !=0 means haven't define atomic contribution
    if (allocated(ple3n1path)) deallocate(ple3n1path)

    do while(.true.)
    !! We have calculated planemat or x/y grad, now use the data to draw graph
    !! Because the max cycle is ngridnum1/2 - 1, so the upper coordinate likes dist1-d1 rather than dist1
        if ((i==-1.and.isilent==0).or.isavepic==1) then
            if (idrawtype==3.or.idrawtype==4.or.idrawtype==5) then !Draw 3D plane, first use drawplaneGUI to setup GUI, which then invokes drawplane
            else
            end if
            if (isavepic==1) write(*,"(a,a,a)") " Graph have been saved to ",trim(graphformat)," file with ""DISLIN"" prefix in current directory"
            isavepic=0
        end if
        
        !! After show the plot once, ask user what to do next. 0 and negative options are general options
        write(*,*)
!         write(*,*) "-10 Multiply data by the data in a plane text file"
        if (iatmcontri==0) write(*,*) "-9 Only plot the data around certain atoms"
        if (iatmcontri==1) write(*,*) "-9 Recovery original plane data"
        if (ilenunit2D==1) write(*,*) "-8 Change length unit of the graph to Angstrom"
        if (ilenunit2D==2) write(*,*) "-8 Change length unit of the graph to Bohr"        
        write(*,*) "-7 Multiply data by a factor"
        write(*,*) "-6 Export calculated plane data to plane.txt in current folder"
        write(*,*) "-5 Return to main menu"
        write(*,*) "-4 Switch ON/OFF of reversing ticks"
        write(*,"(a,i3)") " -3 Set the number of ticks between the labels, current:",iticks
        if (idrawtype==1.or.idrawtype==4.or.idrawtype==5) then
            write(*,"(a,2f7.3,f10.5)") " -2 Set stepsize in X,Y,Z axes, current:",planestpx,planestpy,planestpz
        else
            write(*,"(a,2f7.3)") " -2 Set stepsize in X and Y axes, current:",planestpx,planestpy
        end if
        write(*,*) "-1 Show the graph again"
        write(*,*) "0 Save the graph to a file"

        if (idrawtype==1) then !Color-filled map
            if (abs(drawlowlim)<1000000.and.abs(drawuplim)<1000000) then
                write(*,"(a,2f15.7)") " 1 Set lower&upper limit of color scale, current:",drawlowlim,drawuplim
            else
                write(*,"(a,2(1PE15.6))") " 1 Set lower&upper limit of color scale, current:",drawlowlim,drawuplim
            end if
            if (idrawcontour==1) write(*,*) "2 Disable showing contour lines"
            if (idrawcontour==0) write(*,*) "2 Enable showing contour lines"
            write(*,*) "3 Change contour line setting"
            if (iatom_on_contour==0) write(*,*) "4 Enable showing atom labels and reference point"
            if (iatom_on_contour==1) write(*,*) "4 Disable showing atom labels and reference point"
            if (idrawplanevdwctr==0.and.iorbsel2==0.and.allocated(b)) write(*,*) "15 Draw a contour line of vdW surface (electron density=0.001)" !meaningless if custom operation is performed
            if (idrawplanevdwctr==1) write(*,*) "15 Don't draw a contour line of vdW surface"
            if (idrawplanevdwctr==1) write(*,*) "16 Set label size, style and color of the contour line of vdW surface" !When iorbsel2/=0, that means plot another orbital, cubmattmp will be pre-occupied
            if (iatom_on_contour==1) write(*,"(a,f7.3)") " 17 Set the distance criterion for showing atom labels, current:",disshowlabel
            if (iatmlabtype==1) write(*,*) "18 Change style of atomic labels: Only plot element symbol"
            if (iatmlabtype==2) write(*,*) "18 Change style of atomic labels: Only plot atomic index"
            if (iatmlabtype==3) write(*,*) "18 Change style of atomic labels: Plot both element symbol and atomic index"
        else if (idrawtype==2.or.idrawtype==6.or.idrawtype==7) then !contour map, gradient line, vector field with/without contour
            if (iatom_on_contour==0) write(*,*) "1 Enable showing atom labels and reference point"
            if (iatom_on_contour==1) write(*,*) "1 Disable showing atom labels and reference point"
            if (ilabel_on_contour==0) write(*,*) "2 Enable showing isovalue on contour lines"
            if (ilabel_on_contour==1) write(*,*) "2 Disable showing isovalue on contour lines"
            write(*,*) "3 Change contour line setting"
            if (numcp>0.or.numpath>0) write(*,*) "4 Set marks of critical points and paths"
            if (idrawcontour==1) write(*,*) "5 Disable showing contour lines"
            if (idrawcontour==0) write(*,*) "5 Enable showing contour lines"
            if (numcp>0.and.idrawintbasple==0) write(*,*) "6 Generate and show interbasin paths"
            if (numcp>0.and.idrawintbasple==1) write(*,*) "6 Delete interbasin paths"
            if (numcp>0.and.idrawintbasple==0) write(*,*) "7 Set stepsize and maximal iteration for interbasin path generation" 

            if (idrawtype==6) then
                if (igrad_arrow==0) write(*,*) "10 Show arrows"
                if (igrad_arrow==1) write(*,*) "10 Don't show arrows"
                write(*,"(a,f8.4)") " 11 Set integration step for gradient line, current:",gradplotstep
                write(*,"(a,f8.4)") " 12 Set interstice between gradient line, current:",gradplotdis
                write(*,"(a,f8.4)") " 13 Set test value for drawing a new gradient line, current:",gradplottest
                write(*,*) "14 Set color and line width for gradient lines"
            else if (idrawtype==7) then
                write(*,"(a,f8.4)") " 10 Set upper limit of absolute value for scaling, current:",cutgradvec
                if (icolorvecfield==0) write(*,*) "11 Map color to arrows"
                if (icolorvecfield==1) write(*,*) "11 Don't Map color to arrows"
                if (icolorvecfield==0) write(*,*) "12 Set color for arrow heads" !If color map was set, the color set by user themselves is nulified
                if (iinvgradvec==0) write(*,*) "13 Invert gradient vectors"
                if (iinvgradvec==1) write(*,*) "13 Don't Invert gradient vectors"
            end if            
            if (idrawplanevdwctr==0.and.iorbsel2==0.and.allocated(b)) write(*,*) "15 Draw a contour line of vdW surface (electron density=0.001)" !meaningless if custom operation is performed
            if (idrawplanevdwctr==1) write(*,*) "15 Don't draw a contour line of vdW surface"
            if (idrawplanevdwctr==1) write(*,*) "16 Set label size, style and color of the contour line of vdW surface"
            if (iatom_on_contour==1) write(*,"(a,f7.3)") " 17 Set the distance criterion for showing atom labels, current:",disshowlabel
            if (iatmlabtype==1) write(*,*) "18 Change style of atomic labels: Only plot element symbol"
            if (iatmlabtype==2) write(*,*) "18 Change style of atomic labels: Only plot atomic index"
            if (iatmlabtype==3) write(*,*) "18 Change style of atomic labels: Plot both element symbol and atomic index"
        else if (idrawtype==4.or.idrawtype==5) then !Colored relief map with/without projected color-filled map
            write(*,*) "1 Set color scale range for filling color"
            write(*,"(a,a)") " 2 Toggle drawing mesh on the surface, current: ",drawsurmesh
        else if (idrawtype==3) then
            continue
        end if
        
        read(*,*) i
        
        !Below are general options
        if (i==-9) then
            if (iatmcontri==0) then
                allocate(planemat_bk(ngridnum1,ngridnum2))
                planemat_bk=planemat
                if (allocated(tmparrint)) deallocate(tmparrint)
                write(*,*) "Input index of the atoms you are interested in, e.g. 2,3,7-10"
                read(*,"(a)") c2000tmp
                call str2arr(c2000tmp,ntmp)
                allocate(tmparrint(ntmp))
                call str2arr(c2000tmp,ntmp,tmparrint)
                write(*,*) "Updating plane data, please wait..."
                do i=1,ngridnum1 !First calculate promolecular density and store it to planemat
                    do j=1,ngridnum2
                        rnowx=orgx2D+(i-1)*v1x+(j-1)*v2x
                        rnowy=orgy2D+(i-1)*v1y+(j-1)*v2y
                        rnowz=orgz2D+(i-1)*v1z+(j-1)*v2z
                        densall=0
                        densfrag=0
                        do iatm=1,ncenter
                            tmpval=calcatmdens(iatm,rnowx,rnowy,rnowz,0)
                            densall=densall+tmpval
                            if (any(tmparrint==iatm)) densfrag=densfrag+tmpval
                        end do
                        planemat(i,j)=planemat(i,j)*densfrag/densall
                    end do
                end do
                write(*,*) "Done! The data have been updated, you can replot it"
                deallocate(tmparrint)
                iatmcontri=1
            else !Recovery the backed up data
                planemat=planemat_bk
                deallocate(planemat_bk)
                iatmcontri=0
            end if
        else if (i==-8) then
            if (ilenunit2D==1) then
                ilenunit2D=2
            else if (ilenunit2D==2) then
                ilenunit2D=1
            end if
        else if (i==-10) then !Load plane data in another plain text file and operate to current plane data, the plane settings must be identical
            write(*,*) "Input file name, e.g. C:\plane.txt"
            read(*,"(a)") c200tmp
            write(*,*) "How many columns? (4 or 6. The data in the last column will be loaded)"
            read(*,*) ncol
            open(10,file=c200tmp,status="old")
                do i=0,ngridnum1-1
                    do j=0,ngridnum2-1
                        if (ncol==4) then
                            read(10,*) tmpv,tmpv,tmpv,planemattmp(i+1,j+1)
                        else
                            read(10,*) tmpv,tmpv,tmpv,tmpv,tmpv,planemattmp(i+1,j+1)
                        end if
                    end do
                end do
            close(10)
            write(*,*) "Which operation? +,-,x,/"
            read(*,*) c200tmp(1:1)
            if (c200tmp(1:1)=="+") then
                planemat=planemat+planemattmp
            else if (c200tmp(1:1)=="-") then
                planemat=planemat-planemattmp
            else if (c200tmp(1:1)=="x") then
                planemat=planemat*planemattmp
            else if (c200tmp(1:1)=="/") then
                planemat=planemat/planemattmp
            end if
            write(*,*) "Done!"
        else if (i==-7) then        
            write(*,*) "Input a value, e.g. 0.3"
            read(*,*) scaleval
            planemat=planemat*scaleval
            write(*,*) "Done!"
        else if (i==-6) then
            open(10,file="plane.txt",status="replace")
            do i=0,ngridnum1-1
                do j=0,ngridnum2-1
                    rnowx=orgx2D+i*v1x+j*v2x
                    rnowy=orgy2D+i*v1y+j*v2y
                    rnowz=orgz2D+i*v1z+j*v2z
!                     if (planemat(i+1,j+1)<0.00098.or.planemat(i+1,j+1)>0.00103) cycle
                    if (plesel==4.or.plesel==5.or.plesel==6.or.plesel==7) then
                        write(10,"(5f10.5,1PE18.10)") rnowx*b2a,rnowy*b2a,rnowz*b2a,i*d1*b2a,j*d2*b2a,planemat(i+1,j+1)
                    else !Plane is vertical, the coordinate in a direction is zero
                        write(10,"(3f10.5,1PE18.10)") rnowx*b2a,rnowy*b2a,rnowz*b2a,planemat(i+1,j+1)
                    end if
                end do
            end do
            close(10)
            if (plesel==4.or.plesel==5.or.plesel==6.or.plesel==7) then
                write(*,"(a)") " The column 1,2,3 correspond to actual coordinates respectively"
                write(*,"(a)") " The column 4,5,6 correspond to X,Y coordinates in the graph and function value respectively"
            else
                write(*,"(a)") " The column 1,2,3,4 correspond to X,Y,Z and function value respectively"
            end if
            write(*,"(a)") "Result has been ouputted to plane.txt in current folder, length unit is in Angstrom"
            if (idrawtype/=3.and.idrawtype/=4.and.idrawtype/=5) then
                if ( numcp>0.or.(numpath>0.and.imarkpath==1).or.(nple3n1path>0.and.idrawintbasple==1) ) then
                    write(*,"(/,a)") " If also output critical points and topology/basin paths to plain text files in current folder? (y/n)"
                    read(*,*) selectyn
                    if (selectyn=='y'.or.selectyn=='Y') then
                        iplaneoutall=1 !Global variable, which tells drawplane routine to export topology data
                        iplaneoutall=0
                        if (numcp>0) write(*,"(a)") " Critical points have been outputted to planeCP.txt in current folder. The third column is type: 1=(3,-3), 2=(3,-1), 3=(3,+1), 4=(3,+3)"
                        if (numpath>0.and.imarkpath==1) write(*,*) "Topology paths have been outputted to planepath.txt in current folder"
                        if (nple3n1path>0.and.idrawintbasple==1) write(*,*) "Interbasin paths have been outputted to planeinterbasin.txt in current folder"
                        write(*,"(a)") " The first two columns correspond to X,Y coordinates in the graph, the unit is in Angstrom"
                    end if
                end if
            end if
        else if (i==-5) then
            deallocate(planemat,planemattmp)
            if (allocated(planemat_bk)) deallocate(planemat_bk)
            idrawplanevdwctr=0
            iorbsel2=0
            exit
        else if (i==-4) then
            if (itickreverse==0) then
                itickreverse=1
            else
                itickreverse=0
            end if
        else if (i==-3) then
            write(*,*) "How many ticks between labels do you want?"
            read(*,*) iticks
            iticks=iticks+1
        else if (i==-2) then
            if (idrawtype==1.or.idrawtype==4.or.idrawtype==5) then
                write(*,"(a)") " Input the stepsize between the labels in X,Y,Z axes, e.g. 1.5,2.0,0.1"
                read(*,*) planestpx,planestpy,planestpz
            else
                write(*,"(a)") " Input the stepsize between the labels in X and Y axes, e.g. 1.5,2.0"
                read(*,*) planestpx,planestpy
            end if
        else if (i==-1) then
            cycle
        else if (i==0) then
            isavepic=1
        end if

        !Shared options for idrawtype 1,2,6,7 are given first
        if (idrawtype==1.or.idrawtype==2.or.idrawtype==6.or.idrawtype==7) then
            if (i==15) then
                if (idrawplanevdwctr==0) then
                    idrawplanevdwctr=1
                    write(*,*) "Please wait..."
nthreads=getNThreads()
!$OMP PARALLEL DO private(ipt,jpt,rnowx,rnowy,rnowz) shared(planemattmp) schedule(dynamic) NUM_THREADS(nthreads)
                    do ipt=0,ngridnum1-1
                        do jpt=0,ngridnum2-1
                            rnowx=orgx2D+ipt*v1x+jpt*v2x
                            rnowy=orgy2D+ipt*v1y+jpt*v2y
                            rnowz=orgz2D+ipt*v1z+jpt*v2z
                            planemattmp(ipt+1,jpt+1)=fdens(rnowx,rnowy,rnowz)
                        end do
                    end do
!$OMP END PARALLEL DO
                    write(*,*) "Done, now you can replot the graph to check effect"
                else if (idrawplanevdwctr==1) then
                    idrawplanevdwctr=0
                end if
            else if (i==16) then
                write(*,*) "Input the size of the label, e.g. 30"
                write(*,*) "Note: If you input 0, then the label will not be shown"
                read(*,*) ivdwctrlabsize
                write(*,*) "Select color of the contour line:"
                write(*,*) "1/2/3/4/5 = Red/Green/Blue/White/Black"
                write(*,*) "6/7/8/9/10 = Gray/Cyan/Yellow/Orange/Magenta"
                write(*,*) "11/12/13/14 = Crimson/Dark green/Purple/Brown"
                read(*,*) ivdwclrindctr
                write(*,*) "Input the width of the contour line, e.g. 10"
                read(*,*) iwidthvdwctr
                write(*,*) "Input length of line segment and interstice"
                write(*,*) "e.g. 1,0 means solid line; 1,10 means DOT; 10,15 means DASH"
                write(*,*) "     10,25 means DASH with larger interstice"
                read(*,*) vdwctrstyle
                write(*,*) "Done, now you can replot the graph to check effect"
            else if (i==17) then
                write(*,"(a)") " Note: If the distance between an atom nucleus/critical point and the plane of interest is less than this criterion, &
                the atom/critical point will be labelled on the graph"
                write(*,*) "Input the criterion value in Bohr, e.g. 0.5"
                read(*,*) disshowlabel
                write(*,"(a)") " If also show the labels of the atoms that beyond this criterion as light face type? (y/n)"
                read(*,*) selectyn
                if (selectyn=='y'.or.selectyn=='Y') then
                    iatom_on_contour_far=1
                else
                    iatom_on_contour_far=0
                end if
            else if (i==18) then
                write(*,*) "1: Only plot element symbol"
                write(*,*) "2: Only plot atomic index"
                write(*,*) "3: Plot both element symbol and atomic index"
                write(*,*) "Note that the default value can be set by ""iatmlabtype"" in settings.ini"
                read(*,*) iatmlabtype
            end if
            
            !Options only for idrawtype 1=====================
            if (idrawtype==1) then 
                if (i==1) then
                    write(*,*) "Input lower & upper limit of Z (e.g. -0.3,0.3)"
                    read(*,*) drawlowlim,drawuplim
                else if (i==2) then
                    if (idrawcontour==1) then
                        idrawcontour=0
                    else if (idrawcontour==0) then
                        idrawcontour=1
                    end if
                else if (i==3) then  !! Change isovalues of contour line
                    call setctr
                else if (i==4) then
                    if (iatom_on_contour==1) then
                        iatom_on_contour=0
                    else if (iatom_on_contour==0) then
                        iatom_on_contour=1
                        write(*,*) "Use which color for labelling atoms?"
                        write(*,*) "1/2/3/4/5 = Red/Green/Blue/White/Black"
                        write(*,*) "6/7/8/9/10 = Gray/Cyan/Yellow/Orange/Magenta"
                        write(*,*) "11/12/13/14 = Crimson/Dark green/Purple/Brown"
                        read(*,*) iclrindatmlab
                    end if
                end if
            !Options only for idrawtype 2,6,7========================
            else if (idrawtype==2.or.idrawtype==6.or.idrawtype==7) then
                !General option for idrawtype 2,6,7
                if (i==1) then
                    if (iatom_on_contour==1) then
                        iatom_on_contour=0
                    else if (iatom_on_contour==0) then
                        iatom_on_contour=1
                        write(*,*) "Use which color?"
                        write(*,*) "1/2/3/4/5 = Red/Green/Blue/White/Black"
                        write(*,*) "6/7/8/9/10 = Gray/Cyan/Yellow/Orange/Magenta"
                        write(*,*) "11/12/13/14 = Crimson/Dark green/Purple/Brown"
                        read(*,*) iclrindatmlab
                    end if
                else if (i==2) then
                    if (ilabel_on_contour==1) then
                        ilabel_on_contour=0
                    else if (ilabel_on_contour==0) then
                        ilabel_on_contour=1
                        write(*,*) "Input label size   e.g. 30"
                        read(*,*) ictrlabsize
                        write(*,"(a)") " Hint: The number of digits after the decimal point of label on contour lines can be set by ""numdigctr"" in settings.ini"
                    end if
                else if (i==3) then  !! Change isovalues of contour line
                    call setctr
                else if (i==4) then
                    call settopomark
                else if (i==5) then
                    if (idrawcontour==1) then
                        idrawcontour=0
                    else if (idrawcontour==0) then
                        idrawcontour=1
                    end if
                else if (i==6) then !Generate paths from (3,-1)
                    if (idrawintbasple==0) then    
                        do icp=1,numcp
                            if (CPtype(icp)==2) then
                                cpx=CPpos(1,icp)
                                cpy=CPpos(2,icp)
                                cpz=CPpos(3,icp)
                                if (plesel==1) then
                                    if (abs(cpz-orgz2D) > disshowlabel) cycle
                                else if (plesel==2) then
                                    if (abs(cpy-orgy2D) > disshowlabel) cycle
                                else if (plesel==3) then
                                    if (abs(cpx-orgx2D) > disshowlabel) cycle
                                else if (plesel==4.or.plesel==5.or.plesel==6.or.plesel==7) then
                                    call pointprjple(a1x,a1y,a1z,a2x,a2y,a2z,a3x,a3y,a3z,cpx,cpy,cpz,prjx,prjy,prjz)
                                    if ( (cpx-prjx)**2+(cpy-prjy)**2+(cpz-prjz)**2 > disshowlabel**2) cycle
                                end if
                                nple3n1path=nple3n1path+1
                                cp2ple3n1path(icp)=nple3n1path !Default is zero, means this CP is not on the given plane and has no corresponding interbasin path
                            end if
                        end do
                        if (nple3n1path>0) then
                            idrawintbasple=1
                            write(*,"(' Found',i8,' (3,-1) CPs in the plane')") nple3n1path
                            allocate(ple3n1path(3,n3n1plept,2,nple3n1path))
                            write(*,*) "Generating interbasin paths from (3,-1) CPs, Please wait..."
nthreads=getNThreads()
!$OMP PARALLEL DO SHARED(numcp) PRIVATE(icp) schedule(dynamic) NUM_THREADS(nthreads)
                            do icp=1,numcp
                                if (cp2ple3n1path(icp)/=0) call gen3n1plepath(ifunctopo,icp,cp2ple3n1path(icp))
        !                         write(*,"('Finished the in-plane path generation from (3,-1)',i8)") icp
                            end do
!$OMP END PARALLEL DO
                        else
                            write(*,*) "No (3,-1) CP is closed to the plane"
                        end if
                    else if (idrawintbasple==1) then
                        idrawintbasple=0
                        nple3n1path=0
                        cp2ple3n1path=0
                        if (allocated(ple3n1path)) deallocate(ple3n1path)
                    end if
                else if (i==7) then
                    write(*,*) "Input stepsize (Bohr) and maximal iterations"
                    write(*,"(a,f8.5,',',i6)") " Note: Current values:",ple3n1pathstpsiz,n3n1plept
                    read(*,*) ple3n1pathstpsiz,n3n1plept
                end if
                
                !Option only for idrawtype 6
                if (idrawtype==6) then
                    if (i==10) then
                        if (igrad_arrow==1) then
                            igrad_arrow=0
                        else if (igrad_arrow==0) then
                            igrad_arrow=1
                        end if
                    else if (i==11) then
                        write(*,*) "Input a value, default value is 0.002D0"
                        read(*,*) gradplotstep
                    else if (i==12) then
                        write(*,"(a,f8.4)") "Input a value, should be larger than",gradplotstep
                        read(*,*) gradplotdis
                    else if (i==13) then
                        write(*,*) "Input a value, default value is 0.2"
                        read(*,*) gradplottest
                    else if (i==14) then
                        write(*,*) "Use which color?"
                        write(*,*) "1/2/3/4/5 = Red/Green/Blue/White/Black"
                        write(*,*) "6/7/8/9/10 = Gray/Cyan/Yellow/Orange/Magenta"
                        write(*,*) "11/12/13/14 = Crimson/Dark green/Purple/Brown"
                        read(*,*) iclrindgradline
                        write(*,*) "Input line width, e.g. 5  (default value is 1)"
                        read(*,*) iwidthgradline
                    end if
                !Individual option for idrawtype 7    
                else if (idrawtype==7) then
                    if (i==10) then
                        write(*,*) "Input a value"
                        read(*,*) cutgradvec
                    else if (i==11) then
                        if (icolorvecfield==1) then
                            icolorvecfield=0
                        else if (icolorvecfield==0) then
                            icolorvecfield=1
                        end if
                    else if (i==12) then
                        write(*,*) "Input color index, e.g. 1 =black, 50 =blue, 150 =green, 250 =red"
                        read(*,*) vecclrind
                    else if (i==13) then
                        if (iinvgradvec==1) then
                            iinvgradvec=0
                        else if (iinvgradvec==0) then
                            iinvgradvec=1
                        end if
                    end if
                end if
            end if
            
        !Options for idrawtype 4,5=================
        else if (idrawtype==4.or.idrawtype==5) then
            if (i==1) then
                write(*,*) "Input lower & upper limits of color scale for shading the surface"
                write(*,*) "e.g. -0.1,0.3"
                read(*,*) surcolorzmin,surcolorzmax
                if (idrawtype==5) then
                    write(*,*) "Input lower & upper limits of color scale for the projected map"
                    write(*,*) "e.g. -0.1,0.3"
                    read(*,*) drawlowlim,drawuplim
                end if
            else if (i==2) then
                if (drawsurmesh=="ON ") then
                    drawsurmesh="OFF"
                else
                    drawsurmesh="ON"
                end if
            end if
        else if (idrawtype==3) then !No options for map type 3 currently
            continue
        end if
    end    do



!!!--------------------------------------------------------
!!!--------------------------------------------------------
!5!------------------- Calculate, show and output grid file
!!!--------------------------------------------------------
!!!--------------------------------------------------------

else if (infuncsel1==5) then
    ncustommap=0 !Clean custom operation setting that possibly defined by other modules
    if (allocated(custommapname)) deallocate(custommapname)
    if (allocated(customop)) deallocate(customop)
    write(*,*) "-10 Return to main menu"
    write(*,*) "-2 Obtain of deformation property"
    write(*,*) "-1 Obtain of promolecule property"
    write(*,*) "0 Set custom operation"
500    call selfunc_interface(infuncsel2)

    if (infuncsel2==0.or.infuncsel2==-1.or.infuncsel2==-2) then
        if (infuncsel2==0) call customplotsetup
        if (infuncsel2==-1) then
            ipromol=1
        else
            ipromol=0
        end if
        if (infuncsel2==-1.or.infuncsel2==-2) call setPromol
        write(*,*) "-10 Return to main menu"
        goto 500
    else if (infuncsel2==-10) then
        cycle
    else if (infuncsel2==111) then !Calculate Becke weighting function
        write(*,*) "Input indices of two atoms to calculate Becke overlap weight, e.g. 1,4"
        write(*,*) "or input index of an atom and zero to calculate Becke atomic weight, e.g. 5,0"
        read(*,*) iatmbecke1,iatmbecke2
    else if (infuncsel2==112) then !Calculate Hirshfeld weighting function
        if (allocated(tmparrint)) deallocate(tmparrint)
        write(*,*) "Input index of the atoms you are interested in, e.g. 2,3,7-10"
        read(*,"(a)") c2000tmp
        call str2arr(c2000tmp,ntmp)
        allocate(tmparrint(ntmp))
        call str2arr(c2000tmp,ntmp,tmparrint)
        write(*,"(a)") " How to generate the atomic densities that used in the calculation of Hirshfeld weight?"
        write(*,*) "1 Based on atomic .wfn files"
        write(*,*) "2 Based on built-in atomic densities (see Appendix 3 of the manual for detail)"
        read(*,*) iHirshdenstype
        if (iHirshdenstype==1) call setpromol
    end if
    
    call setgrid(1,igridsel)
    
    if (igridsel==100) then !Calculate value on a set of points loaded from external file
        if (allocated(extpttmp)) deallocate(extpttmp)
        if (ncustommap/=0) allocate(extpttmp(numextpt)) !!! temp file for difference cube
        if (infuncsel2/=4) call delvirorb(1) !Delete high-lying virtual orbitals for faster calculation
        icustom=0
        extpt(:,4)=0D0
        call walltime(walltime1)
        CALL CPU_TIME(time_begin)
        if (ipromol==1) goto 509 !Calculate promolecular property, so skip the first time calculation (namely for the whole system)
    508 continue
nthreads=getNThreads()
!$OMP PARALLEL DO SHARED(extpt) PRIVATE(iextpt) schedule(dynamic) NUM_THREADS(nthreads)
        do iextpt=1,numextpt !Calculate function value
            extpt(iextpt,4)=calcfuncall(infuncsel2,extpt(iextpt,1),extpt(iextpt,2),extpt(iextpt,3))
        end do
!$OMP END PARALLEL DO
    509    if (ncustommap/=0) then !cycling will stop when all the file have been dealed
            if (icustom==0) then
        !Note: For promolecular property, x,y,z hasn't been saved in extpt at first time, while after calculation of atoms, extpt already has %x,%y,%z
                extpttmp(:)=extpt(:,4) !first time
            else if (icustom/=0) then !not first time
                if (customop(icustom)=='+') extpttmp(:)=extpttmp(:)+extpt(:,4)
                if (customop(icustom)=='-') extpttmp(:)=extpttmp(:)-extpt(:,4)
                if (customop(icustom)=='x'.or.customop(icustom)=='*') extpttmp(:)=extpttmp(:)*extpt(:,4)
                if (customop(icustom)=='/') extpttmp(:)=extpttmp(:)/extpt(:,4)
            end if
            if (icustom/=ncustommap) then
                icustom=icustom+1
                filename=custommapname(icustom)
                call dealloall
                write(*,"(' Loading:  ',a)") trim(filename)
                call readinfile(filename,1)
                if (infuncsel2/=4) call delvirorb(0)
                !Generate temporary fragatm
                deallocate(fragatm)
                nfragatmnum=ncenter
                allocate(fragatm(nfragatmnum))
                do iatm=1,ncenter
                    fragatm(iatm)=iatm
                end do
                !Input the MO index for current file. Since the MO index may be not the same as the first loaded one
                if (infuncsel2==4) then
                    write(*,"(' Input the index of the orbital to be calculated for ',a,'   e.g. 3')") trim(filename)
                    read(*,*) iorbsel
                end if
                goto 508
            else if (icustom==ncustommap) then !last time
                extpt(:,4)=extpttmp(:)
                call dealloall
                write(*,"(' Reloading:  ',a)") trim(firstfilename)
                call readinfile(firstfilename,1)
                !Recovery user defined fragatm from the backup
                deallocate(fragatm)
                nfragatmnum=nfragatmnumbackup
                allocate(fragatm(nfragatmnum))
                fragatm=fragatmbackup
            end if
        end if
        CALL CPU_TIME(time_end)
        call walltime(walltime2)
        if (isys==1) write(*,"(' Calculation is finished, took up CPU time',f10.2,'s, wall clock time',i9,'s')") time_end-time_begin,walltime2-walltime1
        write(*,"(a)") " Output the points with function values to which file? e.g. c:\ltwd.txt"
        read(*,"(a)") c200tmp
        open(10,file=c200tmp,status="replace")
        write(10,"(i10)") numextpt
        do iextpt=1,numextpt
            write(10,"(3f13.7,E20.10)") extpt(iextpt,1:3),extpt(iextpt,4)
        end do
        close(10)
        write(*,"(a)") " Done! In this file the first line is the number of points, Column 1~4 correspond to X,Y,Z coordinates and function values, respectively. All units are in a.u."
    
    else !Calculate grid data
        if (allocated(cubmat)) deallocate(cubmat)
        allocate(cubmat(nx,ny,nz))
        if (allocated(cubmattmp)) deallocate(cubmattmp)
        if (ncustommap/=0) allocate(cubmattmp(nx,ny,nz)) !!! temp file for difference cube
        if (infuncsel2/=4) call delvirorb(1) !Delete high-lying virtual orbitals for faster calculation
        !!!!! Calculate grid data
        if (infuncsel2==111) then !Becke's weight
            do k=1,nz
                tmpz=orgz+(k-1)*dz
                do j=1,ny
                    tmpy=orgy+(j-1)*dy
                    do i=1,nx
                        tmpx=orgx+(i-1)*dx
                        cubmat(i,j,k)=beckewei(tmpx,tmpy,tmpz,iatmbecke1,iatmbecke2)
                    end do
                end do
            end do
        else if (infuncsel2==112) then !Hirshfeld weight
            ncustommap=0
            call genhirshcubewei(tmparrint,size(tmparrint),iHirshdenstype)
        else if (infuncsel2==120) then !Calculate and output three components of Steric force to plain text file
            open(20,file="stericforce.txt",status="replace")
            do k=1,nz
                write(*,"(' Finished:',i5,'  /',i5)") k,nz
                tmpz=orgz+(k-1)*dz
                do j=1,ny
                    tmpy=orgy+(j-1)*dy
                    do i=1,nx
                        tmpx=orgx+(i-1)*dx
                        call stericderv(tmpx,tmpy,tmpz,tmpvec)
                        write(20,"(7f12.6)") (orgx+(i-1)*dx)*b2a,(orgy+(j-1)*dy)*b2a,(orgz+(k-1)*dz)*b2a,-tmpvec, dsqrt(sum(tmpvec**2))
                    end do
                end do
            end do
            close(20)
            write(*,*) "Done, the results have been outputted to stericforce.txt in current folder"
            write(*,"(a)") " Columns 1,2,3 correspond to X,Y,Z coordinates, 4,5,6 correspond to steric force component in X,Y,Z. The last column denotes magnitude of steric force"
            write(*,*)
            pause
        else !Common cases
            cubmat=0D0
            icustom=0
            if (ipromol==1) goto 511 !Calculate promolecular property, so skip the first time calculation (namely for the whole system)
        510    call savecubmat(infuncsel2,0,iorbsel) !Save data to cubmat matrix
        511    if (ncustommap/=0) then !Calculate data for custom cube, cycling stop when all the file have been dealed
                if (icustom==0) then
                !Note: For promolecular property, x,y,z hasn't been saved in cubmat at first time, while after calculation of atoms, cubmat already has %x,%y,%z
                    cubmattmp=cubmat !first time
                else if (icustom/=0) then !not first time
                    if (customop(icustom)=='+') cubmattmp=cubmattmp+cubmat
                    if (customop(icustom)=='-') cubmattmp=cubmattmp-cubmat
                    if (customop(icustom)=='x'.or.customop(icustom)=='*') cubmattmp=cubmattmp*cubmat
                    if (customop(icustom)=='/') cubmattmp=cubmattmp/cubmat
                end if
                if (icustom/=ncustommap) then
                    icustom=icustom+1
                    filename=custommapname(icustom)
                    call dealloall
                    write(*,"(' Loading:  ',a)") trim(filename)
                    call readinfile(filename,1)
                    if (infuncsel2/=4) call delvirorb(0)
                    !Generate temporary fragatm
                    deallocate(fragatm)
                    nfragatmnum=ncenter
                    allocate(fragatm(nfragatmnum))
                    do iatm=1,ncenter
                        fragatm(iatm)=iatm
                    end do
                    !Input the MO index for current file. Since the MO index may be not the same as the first loaded one
                    if (infuncsel2==4) then
                        write(*,"(' Input the index of the orbital to be calculated for ',a,'   e.g. 3')") trim(filename)
                        read(*,*) iorbsel
                    end if
                    goto 510
                else if (icustom==ncustommap) then !last time
                    cubmat=cubmattmp
                    call dealloall
                    write(*,"(' Reloading:  ',a)") trim(firstfilename)
                    call readinfile(firstfilename,1)
                    !Recovery user defined fragatm from the backup
                    deallocate(fragatm)
                    nfragatmnum=nfragatmnumbackup
                    allocate(fragatm(nfragatmnum))
                    fragatm=fragatmbackup
                end if
            end if
        end if

        !! Output result
        outcubfile="griddata.cub" !General name
        if (infuncsel2==1) then
            outcubfile="density.cub"
            dipx=0
            dipy=0
            dipz=0
            do k=1,nz
                do j=1,ny
                    do i=1,nx
                        dipx=dipx-cubmat(i,j,k)*(orgx+(i-1)*dx)
                        dipy=dipy-cubmat(i,j,k)*(orgy+(j-1)*dy)
                        dipz=dipz-cubmat(i,j,k)*(orgz+(k-1)*dz)
                    end do
                end do
            end do
            dipx=sum(a%charge*a%x)+dipx*dx*dy*dz
            dipy=sum(a%charge*a%y)+dipy*dx*dy*dz
            dipz=sum(a%charge*a%z)+dipz*dx*dy*dz
            write(*,*)
            write(*,*) "System dipole moment in a.u. (e/Bohr) and Debye, respectively:"
            write(*,"(' X component is',2f12.6)") dipx,dipx*8.47835281D-30*2.99792458D+29
            write(*,"(' Y component is',2f12.6)") dipy,dipy*8.47835281D-30*2.99792458D+29
            write(*,"(' Z component is',2f12.6)") dipz,dipz*8.47835281D-30*2.99792458D+29
            write(*,"(' Total magnitude is   ',2f12.6)") dsqrt(dipx**2+dipy**2+dipz**2),dsqrt(dipx**2+dipy**2+dipz**2)*8.47835281D-30*2.99792458D+29
            write(*,*)
        else if (infuncsel2==2) then
            outcubfile="gradient.cub"
        else if (infuncsel2==3) then
            outcubfile="laplacian.cub"
        else if (infuncsel2==4) then
            outcubfile="MOvalue.cub"
        else if (infuncsel2==5) then
            outcubfile="spindensity.cub"
        else if (infuncsel2==6) then
            outcubfile="K(r).cub"
        else if (infuncsel2==7) then
            outcubfile="G(r).cub"
        else if (infuncsel2==8) then
            outcubfile="nucleiesp.cub"
        else if (infuncsel2==9) then
            outcubfile="ELF.cub"
            sur_value=0.7
        else if (infuncsel2==10) then
            outcubfile="LOL.cub"
            sur_value=0.5
        else if (infuncsel2==11) then
            outcubfile="infoentro.cub"
        else if (infuncsel2==12) then
            outcubfile="totesp.cub"
        else if (infuncsel2==13) then
            sur_value=0.5
            outcubfile="RDG.cub"
        else if (infuncsel2==14) then
            sur_value=0.4
            outcubfile="RDGprodens.cub"
        else if (infuncsel2==15) then
            outcubfile="signlambda2rho.cub"
        else if (infuncsel2==16) then
            outcubfile="signlambda2rhoprodens.cub"
        else if (infuncsel2==17) then
            outcubfile="fermihole.cub"
        else if (infuncsel2==18) then
            outcubfile="avglocion.cub"
        else if (infuncsel2==19) then
            outcubfile="srcfunc.cub"
        else if (infuncsel2==100) then
            outcubfile="userfunc.cub"
        else if (infuncsel2==111) then
            outcubfile="Becke.cub"
        else if (infuncsel2==112) then
            outcubfile="Hirshfeld.cub"
        end if

        temp=minval(cubmat)
        call findvalincub(cubmat,temp,i,j,k)
        write(*,"(' The minimum is',D16.8,' at',3f10.5,' (Bohr)')") temp,orgx+(i-1)*dx,orgy+(j-1)*dy,orgz+(k-1)*dz
        temp=maxval(cubmat)
        call findvalincub(cubmat,temp,i,j,k)
        write(*,"(' The maximum is',D16.8,' at',3f10.5,' (Bohr)')") temp,orgx+(i-1)*dx,orgy+(j-1)*dy,orgz+(k-1)*dz
        write(*,"(' Summing up all value and multiply differential element:')") 
        write(*,*) sum(cubmat)*dx*dy*dz
        write(*,"(' Summing up positive value and multiply differential element:')")
        write(*,*) sum(cubmat,mask=cubmat>0)*dx*dy*dz
        write(*,"(' Summing up negative value and multiply differential element:')")
        write(*,*) sum(cubmat,mask=cubmat<0)*dx*dy*dz

        !Reinitialize plot parameter
        bondcrit=1.15D0
        textheigh=30.0D0
        ratioatmsphere=1.0D0
        bondradius=0.2D0
        ishowatmlab=1
        ishowaxis=1
        idrawmol=1
        isosurshowboth=1
        ishowdatarange=0
        
        do while(.true.)
            write(*,*)
            write(*,*) "-1 Show isosurface graph"
            write(*,*) "0 Return to main menu"
            write(*,*) "1 Save graph of isosurface to file in current folder"
            write(*,*) "2 Export data to Gaussian cube file in current folder"
            write(*,*) "3 Export data to formatted text file in current folder"
            write(*,"(a,f10.5)") " 4 Set the value of isosurface to be shown, current:",sur_value
            write(*,*) "5 Multiply all grid data by a factor"
            write(*,*) "6 Divide all grid data by a factor"
            write(*,*) "7 Add a value to all grid data"
            write(*,*) "8 Substract a value from all grid data"
            read(*,*) i
            
            if (i==-1) then
            else if (i==0) then
                exit
            else if (i==1) then
                idrawisosur=1
                isavepic=1
                isavepic=0
                write(*,*) "Graph has been saved to current folder with ""DISLIN"" prefix"
            else if (i==2) then
                open(10,file=outcubfile,status="replace")
                call outcube(cubmat,nx,ny,nz,orgx,orgy,orgz,dx,dy,dz,10)
                close(10)
                write(*,"(' New grid file has been outputted to: ',a35)") adjustl(outcubfile)
            else if (i==3) then
                open(10,file="output.txt",status="replace")
                write(*,*) "Outputting output.txt in current directory..."
                do i=1,nx
                    do j=1,ny
                        do k=1,nz
                            write(10,"(3f12.6,2x,1PD15.8)") (orgx+(i-1)*dx)*b2a,(orgy+(j-1)*dy)*b2a,(orgz+(k-1)*dz)*b2a,cubmat(i,j,k)
                        end do
                    end do
                end do
                close(10)
                write(*,*) "Output finished, column 1/2/3/4 correspond to X/Y/Z/Value, unit is Angstrom"
            else if (i==4) then
                 write(*,*) "Input the value of isosurface"
                read(*,*) sur_value
            else if (i==5) then
                write(*,*) "Input a value"
                read(*,*) tmpval
                cubmat=cubmat*tmpval
            else if (i==6) then
                write(*,*) "Input a value"
                read(*,*) tmpval
                cubmat=cubmat/tmpval
            else if (i==7) then
                write(*,*) "Input a value"
                read(*,*) tmpval
                cubmat=cubmat+tmpval
            else if (i==8) then
                write(*,*) "Input a value"
                read(*,*) tmpval
                cubmat=cubmat-tmpval
            end if
        end do
    end if

!!!---------------------------------------
!6!!------------------- Check & Modify wavefunction or show GTF/Orbital information
else if (infuncsel1==6) then
    call modwfn


!!!---------------------------------------
!7!!------------------- Population analysis
else if (infuncsel1==7) then
    call popana


!!!---------------------------------------
!8!!------------------- Orbital composition analysis
else if (infuncsel1==8) then
    call compana


!!!---------------------------------------
!9!!------------------- Bond order analysis
else if (infuncsel1==9) then
    call bondana


!!!---------------------------------------
!10!!------------------- Plot DOS
else if (infuncsel1==10) then
    if (ifiletype/=0.or.ifiletype/=1.or.ifiletype/=9) then
        call plotdos
    else
        write(*,"(a,/)") " Error: This function is only available for input file containing basis function information &
        (e.g. .fch/.molden/.gms) and plain text file with energy levels!"
    end if
    
    
!!!---------------------------------------
!11!!------------------- Plot spectrums
else if (infuncsel1==11) then
    call plotspectrum


!!!---------------------------------------
!12!!------------------- Molecular surface analysis
else if (infuncsel1==12) then
    call surfana

    
!!!---------------------------------------
!13!!------------------- Process grid data
else if (infuncsel1==13) then
    if (.not.allocated(cubmat)) then
        write(*,"(a)") " Error: Grid data was not loaded or generated! If you want to load a grid data now, input its path, e.g. C:\nico.cub, else input 0 to exit"
        do while(.true.)
            read(*,"(a)") c200tmp
            if (c200tmp(1:1)=='0') then
                exit
            else
                inquire(file=c200tmp,exist=alive)
                if (alive) then
                    inamelen=len_trim(c200tmp)
                    !Only load grid data, do not pertube other variables
                    if (c200tmp(inamelen-2:inamelen)=="grd") then
                        call readgrd(c200tmp,1,1)
                    else if (c200tmp(inamelen-2:inamelen)=="cub".or.c200tmp(inamelen-3:inamelen)=="cube") then
                        call readcube(c200tmp,1,1)
                    else
                        write(*,*) "Error: Unknown file type, input again"
                        cycle
                    end if
                    call procgriddata
                    exit
                else
                    write(*,*) "Cannot find the file, input again"
                end if
            end if
        end do
    else
        call procgriddata
    end if


!!!---------------------------------------
!14!!------------------- Adaptive natural density partitioning (AdNDP)
else if (infuncsel1==14) then
    call AdNDP
    
    
!!!---------------------------------------
!15!!------------------- Integrate fuzzy atomic space
else if (infuncsel1==15) then
    call intatomspace(0)
    

!!!---------------------------------------
!16!!------------------- Charge decomposition analysis
else if (infuncsel1==16) then
    write(*,*) "Citation of original GCDA and CDA used in Multiwfn, respectively:"
    write(*,"(a)") " Meng Xiao, Tian Lu, Generalized Charge Decomposition Analysis (GCDA) Method, J. Adv. Phys. Chem., 4, 111-124 (2015), http://dx.doi.org/10.12677/JAPC.2015.44013"
    write(*,"(a)") " Stefan Dapprich, Gernot Frenking, J. Phys. Chem., 99, 9352-9362 (1995)"
    call CDA


!!!---------------------------------------
!17!!------------------- Basin integration
else if (infuncsel1==17) then
    call basinana


!!!---------------------------------------
!18!!------------------- Electron excitation analysis
else if (infuncsel1==18) then
    do while(.true.)
        write(*,*) "            ------------ Electron excitation analyses ------------ "
        write(*,*) "0 Return"
        write(*,"(a)") " 1 Analyze and visualize hole-electron distribution, transition dipole moment and transition density"
        write(*,*) "2 Plot transition density matrix as color-filled map"
        write(*,*) "3 Analyze charge-transfer based on density difference grid data (JCTC,7,2498)"
        write(*,*) "4 Calculate delta_r index to measure charge-transfer length (JCTC,9,3118)"
        write(*,*) "5 Calculate transition dipole moments between all excited states"
        
        read(*,*) isel
        if (isel==0) then
            goto 10
        else if (isel==1) then
            call hetransdipdens(1)
        else if (isel==2) then
            call plottransdensmat
        else if (isel==3) then
            if (allocated(cubmat)) then
                call CTanalyze
            else
                write(*,"(a,/)") " Error: Grid data of electron density difference must be calculated by main function 5 or loaded from external file first!"
            end if
        else if (isel==4) then
            call hetransdipdens(2)
        else if (isel==5) then
            call exctransdip
        end if
    end do
    
    
!!!---------------------------------------
!100!!------------------- Misc and some not important functions, Part 1
else if (infuncsel1==100) then
    do while(.true.)
        write(*,*) "              ------------ Other functions (Part 1) ------------ "
        write(*,*) "0 Return"
        write(*,*) "1 Draw scatter graph between two functions and generate their cube files"
        write(*,*) "2 Export .pdb/.xyz/.wfn/.wfx/.molden/.fch/.47 or Gaussian/GAMESS-US input file"
        write(*,*) "3 Calculate molecular van der Waals Volume"
        write(*,*) "4 Integrate a function in whole space"
        write(*,*) "5 Show overlap integral between alpha and beta orbitals"
        write(*,*) "6 Monitor SCF convergence process of Gaussian"
        write(*,*) "7 Generate Gaussian input file with initial guess from converged wavefunction"
        write(*,*) "8 Generate Gaussian input file with initial guess from fragment wavefunctions"
        write(*,*) "9 Evaluate coordination number for all atoms"
        write(*,*) "11 Calculate overlap and centroid distance between two orbitals"
        write(*,*) "13 Calculate HOMA and Bird aromaticity index"
        write(*,*) "14 Calculate LOLIPOP (LOL Integrated Pi Over Plane)"
        write(*,*) "15 Calculate intermolecular orbital overlap"
!         write(*,*) "16 Calculate intermolecular tranfer integral by direct coupling method"
         write(*,*) "18 Yoshizawa's electron transport route analysis"
         write(*,*) "19 Generate promolecular .wfn file from fragment wavefunctions"
         write(*,*) "20 Calculate Hellmann-Feynman forces"
         write(*,*) "21 Calculate properties based on geometry information for specific atoms"
         write(*,*) "22 Detect pi orbitals and set occupation numbers"
         write(*,*) "23 Fit function distribution to atomic value"
         write(*,*) "24 Obtain NICS_ZZ for non-planar system"
         write(*,*) "25 Calculate area and perimeter for a ring"

        read(*,*) infuncsel2
        if (infuncsel2==0) then
            goto 10
        else if (infuncsel2==1) then
            call funcvsfunc
        else if (infuncsel2==2) then
            write(*,*) "1 Output current structure to .pdb file"
            write(*,*) "2 Output current structure to .xyz file"
            write(*,*) "4 Output current wavefunction as .wfx file"
            write(*,*) "5 Output current wavefunction as .wfn file"
            write(*,*) "6 Output current wavefunction as Molden input file (.molden)"
            write(*,*) "7 Output current wavefunction as .fch file"
            write(*,*) "8 Output current wavefunction as .47 file"
            write(*,*) "10 Output current structure to Gaussian input file"
            write(*,*) "11 Output current structure to GAMESS-US input file"
            read(*,*) itmp
            if (itmp==1) then
                write(*,*) "Input the path for pdb file, e.g. c:\ltwd.pdb"
                read(*,"(a)") c200tmp
                call outpdb(c200tmp,10)
            else if (itmp==2) then
                write(*,*) "Input the path for xyz file, e.g. c:\ltwd.xyz"
                read(*,"(a)") c200tmp
                call outxyz(c200tmp,10)
            else if (itmp==4) then
                if (.not.allocated(b)) then
                    write(*,*) "Error: The input file you used does not contain GTF information!"
                else
                    write(*,*) "Input the path, e.g. c:\ltwd.wfx"
                    read(*,"(a)") c200tmp
                    call outwfx(c200tmp,1,10)
                    write(*,*) "Done!"
                end if
            else if (itmp==5) then
                if (.not.allocated(b)) then
                    write(*,*) "Error: The input file you used does not contain GTF information!"
                else
                    write(*,*) "Input the path, e.g. c:\ltwd.wfn"
                    read(*,"(a)") c200tmp
                    call outwfn(c200tmp,1,1,10)
                    write(*,*) "Done!"
                end if
            else if (itmp==6) then
                if (.not.allocated(CObasa)) then
                    write(*,*) "Error: This function works only when input file contains basis function information"
                else
                    write(*,*) "Input the path, e.g. c:\ltwd.molden"
                    read(*,"(a)") c200tmp
                    call outmolden(c200tmp,10)
                end if
            else if (itmp==7) then
                if (.not.allocated(CObasa)) then
                    write(*,*) "Error: This function works only when input file contains basis function information"
                else
                    write(*,*) "Input the path, e.g. c:\ltwd.fch"
                    read(*,"(a)") c200tmp
                    call outfch(c200tmp,10)
                end if
            else if (itmp==8) then
                if (.not.allocated(CObasa)) then
                    write(*,*) "Error: This function works only when input file contains basis function information"
                else
                    write(*,*) "Input the path, e.g. c:\ltwd.47"
                    read(*,"(a)") c200tmp
                    call out47(c200tmp,10)
                end if
            else if (itmp==10) then
                write(*,*) "Input the path, e.g. c:\ltwd.gjf"
                read(*,"(a)") c200tmp
                call outgjf(c200tmp,10)
            else if (itmp==11) then
                write(*,*) "Input the path, e.g. c:\ltwd.inp"
                read(*,"(a)") c200tmp
                call outGAMESSinp(c200tmp,10)
            end if
        else if (infuncsel2==3) then
            if (MCvolmethod==1) then
                write(*,*) "100*2^i points will be used to evaluate the volume by Monte Carlo method"
                write(*,*) "The volume is defined as superposition of vdW sphere of atoms"
                write(*,*)
                do while(.true.)
                    write(*,*) "Please input i, generally 10 is recommended, big system requires bigger i"
                    write(*,*) "Input 0 can return"
                    read(*,*) pointexp
                    if (pointexp==0) exit
                    call calcvolume(1,pointexp,0D0,1D0)
                end do
            else if (MCvolmethod==2) then
                write(*,*) "100*2^i points will be used to evaluate the volume by Monte Carlo method."
                write(*,"(a)") " The volume is defined as the region encompassed by the isosurface of density equals to x, &
                The box used in Monte Carlo procedure will be enlarged by k multiples vdW radius in each side."
                write(*,"(a)") " Hint: For evaluating the volume encompassed by 0.001 isosurface of density for small molecule, &
                            we suggest you input 9,0.001,1.7"
                write(*,*)
                do while(.true.)
                    write(*,*) "Please input i,x,k (Input 0,0,0 can return)"
                    read(*,*) pointexp,tmpisoval,enlarbox
                    if (pointexp==0.and.tmpisoval==0.and.enlarbox==0) exit
                    call calcvolume(2,pointexp,tmpisoval,enlarbox)
                end do
            end if
        else if (infuncsel2==4) then
            if (ispecial/=1) then
                call selfunc_interface(ifunc)
                call intfunc(ifunc)
            else if (ispecial==1) then
                call intfunc(1)
            end if
        else if (infuncsel2==-4) then
            call intdiffsqr
        else if (infuncsel2==5) then
            call aboverlap
        else if (infuncsel2==6) then
            call monitorscf
        else if (infuncsel2==7) then
            call gengauguess
        else if (infuncsel2==8) then
            call fragguess
        else if (infuncsel2==9) then
            call coordnum
        else if (infuncsel2==11) then
            call ovlpdistorb
        else if (infuncsel2==13) then
            call HOMA_Bird
        else if (infuncsel2==14) then
            call LOLIPOP
        else if (infuncsel2==15) then
            call intmolovlp
        else if (infuncsel2==16) then
            call intmoltransint
        else if (infuncsel2==18) then
            call Yoshieletrans
        else if (infuncsel2==19) then
            call genpromolwfn
        else if (infuncsel2==20) then
            call hellmann_feynman
        else if (infuncsel2==21) then
            call calcgeomprop
        else if (infuncsel2==22) then
            call procpiorb
        else if (infuncsel2==23) then
            call fitfunc
        else if (infuncsel2==24) then
            call utilNICS_ZZ
        else if (infuncsel2==25) then
            call calcringsize
        end if
        write(*,*)
    end do
    
    
!!!---------------------------------------
!200!!------------------- Misc and not important functions, Part 2
else if (infuncsel1==200) then
    do while(.true.)
        write(*,*) "              ------------ Other functions (Part 2) ------------ "
        write(*,*) "0 Return"
        write(*,*) "1 Weak interaction analysis for fluctuation environment by RDG method"
        write(*,*) "2 Calculate atomic and bond dipole moments in Hilbert space"
        write(*,*) "3 Generate cube file for multiple orbital wavefunctions"
        write(*,*) "4 Generate iso-chemical shielding surfaces (ICSS) and related quantities"
        write(*,*) "5 Plot radial distribution function for a real space function"
        write(*,*) "6 Analyze correspondence between orbitals in two wavefunctions"
         write(*,*) "7 Parse output of (hyper)polarizability task of Gaussian"
        write(*,*) "8 Calculate (hyper)polarizability by sum-over-states (SOS) method"
        write(*,*) "9 Calculate average bond length and average coordinate number"
        write(*,*) "10 Output various kinds of integral between orbitals"
        write(*,*) "11 Calculate center, the first and second moments of a real space function"
        write(*,*) "12 Calculate energy index (EI) or bond polarity index (BPI)"
        write(*,*) "13 Pipek-Mezey orbital localization"
        write(*,*) "14 Perform integration within isosurfaces of a real space function"
!         write(*,*) "20 Calculate electronic coupling matrix element by FCD or GMH method"
        read(*,*) infuncsel2
        if (infuncsel2==0) then
            goto 10
        else if (infuncsel2==1) then
            call RDG_MD
        else if (infuncsel2==2) then
            call atmbonddip
        else if (infuncsel2==3) then
            call genmultiorbcube
        else if (infuncsel2==4) then
            call ICSS
        else if (infuncsel2==5) then
            call plotraddis
        else if (infuncsel2==6) then
            call orbcorres
        else if (infuncsel2==7) then
            call parseGauPolar
        else if (infuncsel2==8) then
            call SOS
        else if (infuncsel2==9) then
            call atmavgdist
        else if (infuncsel2==10) then
            call outorbint
        else if (infuncsel2==11) then
            call funcmoment
        else if (infuncsel2==12) then
            call calcEIBPI
        else if (infuncsel2==13) then
            call pipek_mezey
        else if (infuncsel2==14) then
            call intisosurface
        else if (infuncsel2==20) then
            call FCD
        end if
        write(*,*)
    end do
        
else if (infuncsel1==98) then
    call sphatmraddens
! else if (infuncsel1==99) then
!     call AIMbasinint


!!!---------------------------------------
!1000!!------------------- Set special parameters
else if (infuncsel1==1000) then
1000 write(*,*)
    write(*,*) "0 Return to main menu"
    write(*,"(a,3f12.6)") " 1 Set reference point, current(Bohr):",refx,refy,refz
    write(*,"(a,i6)") " 2 Set iuserfunc, current:",iuserfunc
    write(*,"(a,i6)") " 3 Set iskipnuc, current:",iskipnuc
    if (pleA==0D0.and.pleB==0D0.and.pleC==0D0.and.pleD==0D0) then
        write(*,"(a)") " 4 Set the plane for user-defined function 38 (Not defined)"
    else
        write(*,"(a)") " 4 Set the plane for user-defined function 38 (Defined)"
    end if
    write(*,"(a,1PD18.8)") " 5 Set global temporary variable, current:",globaltmp
    write(*,"(a,i3)") " 10 Set the number of threads, current:", getNThreads()
    write(*,*) "100 Check the sanity of present wavefunction"
    read(*,*) i
    if (i==1) then
        write(*,*) "Input x,y,z in Bohr, e.g. 3.0,0.0,1.3"
        read(*,*) refx,refy,refz
        write(*,*) "Done!"
    else if (i==2) then
        write(*,*) "Input an integer, e.g. 24"
        read(*,*) iuserfunc
        write(*,*) "Done!"
    else if (i==3) then
        write(*,*) "Input the index of the nucleus, e.g. 24"
        read(*,*) iskipnuc
        write(*,*) "Done!"
    else if (i==4) then
        write(*,*) "1 Input index of three atoms to define the plane"
        write(*,*) "2 Input XYZ coordinate of three points to define the plane"
        read(*,*) iseldef
        if (iseldef==1) then
            write(*,*) "Input three indices, e.g. 2,4,5"
            read(*,*) i1,i2,i3
            call pointABCD(a(i1)%x,a(i1)%y,a(i1)%z,a(i2)%x,a(i2)%y,a(i2)%z,a(i3)%x,a(i3)%y,a(i3)%z,pleA,pleB,pleC,pleD)
        else if (iseldef==2) then
            write(*,*) "Input coordinate for point 1 (in Bohr), e.g. 1.0,-0.2,0.3"
            read(*,*) xtmp1,ytmp1,ztmp1
            write(*,*) "Input coordinate for point 2 (in Bohr), e.g. 2.0,-0.3,0.1"
            read(*,*) xtmp2,ytmp2,ztmp2
            write(*,*) "Input coordinate for point 3 (in Bohr), e.g. 1.3,-1.2,0.33"
            read(*,*) xtmp3,ytmp3,ztmp3
            call pointABCD(xtmp1,ytmp1,ztmp1,xtmp2,ytmp2,ztmp2,xtmp3,ytmp3,ztmp3,pleA,pleB,pleC,pleD)
        end if
        tmpval=dsqrt(pleA**2+pleB**2+pleC**2)
        write(*,"(' The unit vector normal to the plane is:',3f10.5)") pleA/tmpval,pleB/tmpval,pleC/tmpval
        goto 1000
    else if (i==5) then
        write(*,*) "Input the value"
        read(*,*) globaltmp
    else if (i==10) then
        write(*,*) "Input an integer, e.g. 8"
        read(*,*) iniNThreads
        write(*,*) "Done!"
    else if (i==100) then
        call wfnsanity
    end if

end if

end do !End main cycle




contains
!===================================================================================================
!===================================================================================================
!!!!!!!!!!!!!!!! Below routines are contained in and closely related to main program !!!!!!!!!!!!!!!
!===================================================================================================
!===================================================================================================




!!!--------------------------  Set content of custom plot
subroutine customplotsetup 
integer i
write(*,*) "How many files to deal with? (Excluding the file that has been loaded)"
read(*,*) ncustommap
if (allocated(custommapname)) deallocate(custommapname)
if (allocated(customop)) deallocate(customop)
allocate(custommapname(ncustommap))
allocate(customop(ncustommap))
write(*,*) "The avaliable operator: +,-,*,/"  !* can also be intputted as x
write(*,*) "e.g. -,sob.wfn means minus the property of sob.wfn from the first file"
do i=1,ncustommap
    write(*,"('Input the operator and filename of',i5)") i
    do while(.true.)
        read(*,"(a1,1x,a)") customop(i),custommapname(i)
        inquire(file=custommapname(i),exist=alive)
        if (alive.eqv..true.) then
            exit
        else
            write(*,*) "File not found, input again"
        end if
    end do
end do
end subroutine


!!!----------------------Set contour line
subroutine setctr
integer i,j
character outfilename*80,selectyn*1
do while(.true.)
    if (j/=6) then !If last operation is not saving file
        write(*,*) "Current isovalue line:"
        do j=1,lastctrval  !The lastest contour line serial
            write(*,"(i3,':',f19.8,'   ')",advance='no') j,ctrval(j)
            if (mod(j,3)==0) write(*,*)
        end do
        if (mod(lastctrval,3)/=0) write(*,*)
    end if
    if (lastctrval==0) write(*,*) "None"
    if (allocated(boldlinelist)) write(*,*) "Bolded line index:",boldlinelist
    write(*,*) "1 Save setting and return"
    write(*,*) "2 Replace a value of contour line by inputting"
    write(*,*) "3 Add a new contour line"
    write(*,*) "4 Delete a contour line"
    write(*,*) "5 Delete a range of contour lines"
    write(*,*) "6 Save contour setting to external file"
    write(*,*) "7 Load contour setting from external file"
    write(*,*) "8 Generate contour value by arithmetic progression"
    write(*,*) "9 Generate contour value by geometric series"
    if (.not.allocated(boldlinelist)) write(*,*) "10 Enable some contour lines bold style"
    if (allocated(boldlinelist)) write(*,*) "10 Disable some contour lines bold style"
    write(*,*) "11 Set color for positive contour lines"
    write(*,*) "12 Set line style and width for positive contour lines"
    write(*,*) "13 Set color for negative contour lines"
    write(*,*) "14 Set line style and width for negative contour lines"
    read(*,*) j 
    if (j==1) then
        exit
    else if (j==2) then
        write(*,*) "Replace which value of contour line by what value?"
        write(*,*) "(e.g. Input 4,0.015 to replace the fourth value of contour line by 0.015)"
        read(*,*) j,selfctrval
        if (j>=1.and.j<=lastctrval) then
            ctrval(j)=selfctrval
        else
            write(*,*) "The number exceed valid range"
        end if
    else if (j==3) then
        if (lastctrval<size(ctrval)) then
            write(*,*) "Input the value of contour line you want to add"
            read(*,*) selfctrval
            lastctrval=lastctrval+1
            ctrval(lastctrval)=selfctrval
        else
            write(*,"(a,i8)") " Error: The number of contour line couldn't exceed",size(ctrval)
        end if
    else if (j==4) then
        write(*,*) "Delete which contour line? Input a number"
        read(*,*) j
        if (j>=1.and.j<=lastctrval) then
            lastctrval=lastctrval-1
            ctrval(j:lastctrval)=ctrval(j+1:lastctrval+1)
        else
            write(*,*) "The number exceed valid range"
        end if
    else if (j==5) then
        write(*,*) "Delete the range of contour line number (e.g. 3,10)"
        write(*,*) "Note: You can also input 0,0 to delete all"
        read(*,*) j1,j2
        if (j1==0.and.j2==0) then
            lastctrval=0
        else if (j2>j1.and.j1>=1.and.j2<=lastctrval) then
            lastctrval=lastctrval-(j2-j1+1)
            ctrval(j1:lastctrval)=ctrval(j2+1:lastctrval+(j2-j1+1))
        else
            write(*,*) "Invalid input, input excced current range"
        end if
    else if (j==6) then
        write(*,*) "Input the filename for outputting current setting  e.g. c:\ltwd.txt"
        read(*,*) outfilename
        open(10,file=outfilename)
        write(10,"(f19.8)") (ctrval(i),i=1,lastctrval)
        close(10)
        write(*,"(a,a)") " Contour setting has been saved to ",trim(outfilename)
        write(*,*)
    else if (j==7) then
        write(*,*) "Input filename, e.g. C:\ctr.txt"
        do while(.true.)
            read(*,"(a)") extctrsetting
            inquire(file=extctrsetting,exist=alive)
            if (alive .eqv. .true.) exit
            write(*,*) "File not found, input again"
        end do
        open(10,file=extctrsetting,access="sequential",status="old")
        ierror=0
        lastctrval=0
        do i=1,size(ctrval)
            read(10,*,iostat=ierror) temp
            if (ierror/=0) exit
            ctrval(i)=temp
            lastctrval=i
        end do
        close(10)
    else if (j==8.or.j==9) then
        write(*,*) "Input start value, step and total number"
        if (j==8) write(*,*) "e.g. 0.4,0.1,10 generate 0.4,0.5,0.6,0.7 ... 1.3"
        if (j==9) write(*,*) "e.g. 2,3,10 generate 2,6,18,54 ... 39366"
        read(*,*) ctrgenstrval,ctrgenstep,igennum
        write(*,*) "If remove existed contour lines? (y/n)"
        read(*,*) selectyn
        if (selectyn=='y'.or.selectyn=='Y') then
            lastctrval=0
        else if (selectyn=='n'.or.selectyn=='N') then !Append to existed contour lines
            if (lastctrval+igennum>size(ctrval)) then
                igennum=size(ctrval)-lastctrval
                write(*,*) "Warning: The total number of coutour lines exceeded the upper limit!"
            end if
        end if
        do jj=1,igennum
            if (j==8) ctrval(lastctrval+jj)=ctrgenstrval+ctrgenstep*(jj-1)
            if (j==9) ctrval(lastctrval+jj)=ctrgenstrval*ctrgenstep**(jj-1)
        end do
        lastctrval=lastctrval+igennum
    else if (j==10) then
        if (allocated(boldlinelist)) then !Disable bold line
            deallocate(boldlinelist)
            write(*,*) "No line is bolded now, you can select this function again to set bolded lines"
        else
            write(*,*) "Set how many lines bold?"
            read(*,*) numboldline
            allocate(boldlinelist(numboldline))
            do itmp=1,numboldline
                write(*,"(' Input the index of bolded line',i4)") itmp
                read(*,*) boldlinelist(itmp)
            end do
        end if
    else if (j==11) then
        write(*,*) "Use which color?"
        write(*,*) "1/2/3/4/5 = Red/Green/Blue/White/Black"
        write(*,*) "6/7/8/9/10 = Gray/Cyan/Yellow/Orange/Magenta"
        write(*,*) "11/12/13/14 = Crimson/Dark green/Purple/Brown"
        read(*,*) iclrindctrpos
    else if (j==12) then
        write(*,*) "Input length of line segment and interstice"
        write(*,*) "e.g. 1,0 means solid line; 1,10 means DOT; 10,10 means DASH"
        write(*,*) "     10,15 means DASH with larger interstice"
        read(*,*) ctrposstyle(1),ctrposstyle(2)
        write(*,*) "Input line width, e.g. 2"
        read(*,*) iwidthposctr
    else if (j==13) then
        write(*,*) "Use which color?"
        write(*,*) "1/2/3/4/5 = Red/Green/Blue/White/Black"
        write(*,*) "6/7/8/9/10 = Gray/Cyan/Yellow/Orange/Magenta"
        write(*,*) "11/12/13/14 = Crimson/Dark green/Purple/Brown"
        read(*,*) iclrindctrneg
    else if (j==14) then
        write(*,*) "Input length of line segment and interstice"
        write(*,*) "e.g. 1,0 means solid line; 1,10 means DOT; 10,10 means DASH"
        write(*,*) "     10,15 means DASH with larger interstice"
        read(*,*) ctrnegstyle(1),ctrnegstyle(2)
        write(*,*) "Input line width, e.g. 2"
        read(*,*) iwidthnegctr
    end if
end do
end subroutine


!!------------------ Set marker of CPs and paths on contour/gradient map
subroutine settopomark
implicit real*8 (a-h,o-z)
do while(.true.)
    write(*,*) "0 Return"
    if (imark3n3==0) write(*,*) "1 Show (3,-3) CPs"
    if (imark3n3==1) write(*,*) "1 Don't show (3,-3) CPs"
    if (imark3n1==0) write(*,*) "2 Show (3,-1) CPs"
    if (imark3n1==1) write(*,*) "2 Don't show (3,-1) CPs"
    if (imark3p1==0) write(*,*) "3 Show (3,+1) CPs"
    if (imark3p1==1) write(*,*) "3 Don't show (3,+1) CPs"
    if (imark3p3==0) write(*,*) "4 Show (3,+3) CPs"
    if (imark3p3==1) write(*,*) "4 Don't show (3,+3) CPs"
    if (imarkpath==0) write(*,*) "5 Show paths"
    if (imarkpath==1) write(*,*) "5 Don't show paths"
    write(*,*) "10 Set size of markers of CPs"
    write(*,*) "11 Set thickness of paths"
    write(*,*) "12 Set color of paths"
    write(*,*) "13 Set thickness of the interbasin paths derived from (3,-1)"
    write(*,*) "14 Set color of the interbasin paths derived from (3,-1)"
    read(*,*) isel

    if (isel==0) then
        exit
    else if (isel==1) then
        if (imark3n3==1) then
            imark3n3=0
        else
            imark3n3=1
        end if
    else if (isel==2) then
        if (imark3n1==1) then
            imark3n1=0
        else
            imark3n1=1
        end if
    else if (isel==3) then
        if (imark3p1==1) then
            imark3p1=0
        else
            imark3p1=1
        end if
    else if (isel==4) then
        if (imark3p3==1) then
            imark3p3=0
        else
            imark3p3=1
        end if
    else if (isel==5) then
        if (imarkpath==1) then
            imarkpath=0
        else
            imarkpath=1
        end if
    else if (isel==10) then
        write(*,*) "Input a value, e.g. 30"
        read(*,*) sizemarkcp
    else if (isel==11) then
        write(*,*) "Input a value, e.g. 5"
        read(*,*) sizemarkpath
    else if (isel==12) then
        write(*,*) "Input R,G,B value, between 0.0 to 1.0. e.g. 0.5,0.7.1.0"
        write(*,*) "Note: 0,0,0 corresponds to black, 1,1,1 corresponds to white"
        read(*,*) clrRpath,clrGpath,clrBpath
    else if (isel==13) then
        write(*,*) "Input a value, e.g. 5"
        read(*,*) sizemark3n1path
    else if (isel==14) then
        write(*,*) "Input R,G,B value, between 0.0 to 1.0. e.g. 0.5,0.7.1.0"
        write(*,*) "Note: 0,0,0 corresponds to black, 1,1,1 corresponds to white"
        read(*,*) clrR3n1path,clrG3n1path,clrB3n1path
    end if
end do
end subroutine




end program
