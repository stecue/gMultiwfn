!!!============== Molecular surface analysis, output a variety of properties
subroutine surfana
use defvar
use surfvertex
use util
use function
implicit real*8(a-h,o-z)
integer*2,allocatable :: corpos(:,:,:) !corner position
logical,allocatable :: ifbndcub(:,:,:) !if true, means this is a boundary cub
integer,allocatable :: mergerelat(:),HirBecatm(:)
integer tmpintarr3(3)
character pdbfilename*200,c200tmp*200,c80tmp*80,c2000tmp*2000,c10000tmp*10000,selectyn,grdfilename*200,char1tmp
real*8 fragsurarea(ncenter,3),fragsuravg(ncenter,3),fragsurvar(ncenter,3) !Area, average value and variance of each atom surface. 1,2,3 corresponds to all,positive,negative part
real*8 fragsurmax(ncenter),fragsurmin(ncenter),fragsurchgsep(ncenter)
integer surfrag(ncenter),ifatmfrag(ncenter) !User defined fragment contain which atoms; ifatmfrag(iatm)=1/0 means iatm belong / doesn't belong to user defined fragment
integer,allocatable :: surtrifrag(:) !Each surface triangle belongs to which fragment
real*8 nucchgbackup(ncenter) !Backup nuclear charge, because which may be flushed by .chg file
integer isurftype !How to define the surface. 1=Isosurface of electron density, 2=A certain real space function, 5/6=Hirshfeld/Becke surface, 10=Isosurface of existing grid data
real*8 smat(ncenter,ncenter),Pvec(ncenter),tmprarr(ncenter)
isurftype=1
surfisoval=0.001D0
imapfunc=1
ifelim=1 !If elimnate redundant surface vertices
ireadextmapval=0
grdspc=0.25D0
spcmergeratio=0.5D0
critmerge=grdspc*spcmergeratio !If the distance between two surface vertices smaller than this value, merge them
vdwmulti=1.7D0
nbisec=3

surfanaloop: do while(.true.)
do while(.true.)
    tetravol0=0D0
    tetravol1=0D0
    tetravol2=0D0
    tetravol3=0D0
    write(*,*)
    write(*,"(a)") " If this module is used in your research, one should also cite below paper, which described the basic algorithm of this module"
    write(*,*) "Tian Lu, Feiwu Chen, J. Mol. Graph. Model., 38, 314-323 (2012)"
    write(*,*)
    write(*,*) "     ============= Quantitative Molecular Surface Analysis ============="
    write(*,*) "-1 Return to main menu"
    write(*,*) "0 Start analysis now!"
    if (isurftype==1) write(*,"(a,f9.5)") " 1 Select the way to define surface, current: Electron density, iso:",surfisoval
    if (isurftype==2) write(*,"(a,i4,a,f9.5)") " 1 Select the way to define surface, current: Function",ifuncintp,", iso:",surfisoval
    if (isurftype==5) write(*,"(a,i5,a)") " 1 Select the way to define surface, current: Hirshfeld surface, for",nHirBecatm," atoms"
    if (isurftype==6) write(*,"(a,i5,a)") " 1 Select the way to define surface, current: Becke surface, for",nHirBecatm," atoms"
    if (isurftype==10) write(*,"(a,f9.5)") " 1 Select the way to define surface, current: Existing grid data, iso:",surfisoval
    
    if (imapfunc==-1) write(*,*) "2 Select mapped function, current: User defined real space function"
    if (imapfunc==0) write(*,*) "2 Select mapped function, current: A function loaded from external file"
    if (imapfunc==1) write(*,*) "2 Select mapped function, current: Electrostatic potential"
    if (imapfunc==2) write(*,*) "2 Select mapped function, current: Average local ionization energy"
    if (imapfunc==3) write(*,"(a)") " 2 Select mapped function, current: Electrostatic potential from atomic charge"
    if (imapfunc==4) write(*,*) "2 Select mapped function, current: Local electron affinity"
    if (imapfunc==10) write(*,*) "2 Select mapped function, current: Pair density"
    if (imapfunc==11) write(*,*) "2 Select mapped function, current: Electron density"
    if (imapfunc==12) write(*,*) "2 Select mapped function, current: Sign(lambda2)*rho"
    if (imapfunc==20) write(*,*) "2 Select mapped function, current: d_i"
    if (imapfunc==21) write(*,*) "2 Select mapped function, current: d_e"
    if (imapfunc==22) write(*,*) "2 Select mapped function, current: d_norm"
    
    if (isurftype/=10) write(*,"(a,f10.6)") " 3 Spacing of grid points for generating molecular surface:",grdspc
    write(*,*) "4 Advanced options"
    if (ireadextmapval==0) write(*,*) "5 If load mapped function values from external file, current: No"
    if (ireadextmapval==1) write(*,*) "5 If load mapped function values from external file, current: Scatter points"
    if (ireadextmapval==2) write(*,*) "5 If load mapped function values from external file, current: Cubegen output"
    if (ireadextmapval==3) write(*,*) "5 If load mapped function values from external file, current: Cube file"
    
    read(*,*) isel
    if (isel==-1) then
        exit surfanaloop
    else if (isel==0) then
        exit
    else if (isel==1) then
        write(*,*) "How to define the surface?"
        write(*,*) "1: The isosurface of electron density"
        write(*,*) "2: The isosurface of a specific real space function "
        write(*,*) "5: Hirshfeld surface (isosurface of Hirshfeld weight) of a fragment"
        write(*,*) "6: Becke surface (isosurface of Becke weight) of a fragment"
        write(*,*) "10: The isosurface of a grid data loaded from external file"
        if (allocated(cubmat)) write(*,*) "11: The isosurface of the grid data in memory"
        read(*,*) isurftypetmp
        if (isurftypetmp==1) then
            isurftype=1
        else if (isurftypetmp==2) then
            isurftype=2
            call selfunc_interface(ifuncintp)
        else if (isurftypetmp==5.or.isurftypetmp==6) then
            if (isurftypetmp==5) isurftype=5
            if (isurftypetmp==6) isurftype=6
            write(*,"(a)") " Input atomic indices. e.g. 1,3-6,8,10-11 means the atoms 1,3,4,5,6,8,10,11 will be selected"
            read(*,"(a)") c10000tmp
            call str2arr(c10000tmp,nHirBecatm)
            if (allocated(HirBecatm)) deallocate(HirBecatm)
            allocate(HirBecatm(nHirBecatm))
            call str2arr(c10000tmp,nHirBecatm,HirBecatm)
            imapfunc=22
        else if (isurftypetmp==10) then
            isurftype=10
            write(*,*) "Input the path of .cub or .grd file, e.g. c:\ltwd.cub"
            do while(.true.)
                read(*,"(a)") grdfilename
                inquire(file=grdfilename,exist=alive)
                if (alive) exit
                write(*,*) "Cannot find the file, input again"
            end do
            inamelen=len_trim(grdfilename)
            if (grdfilename(inamelen-2:inamelen)=="cub".or.grdfilename(inamelen-3:inamelen)=="cube") then
                call readcube(grdfilename,1,1)
            else if (grdfilename(inamelen-2:inamelen)=="grd") then
                call readgrd(grdfilename,1,1)
            end if
            if (dx/=dy.or.dy/=dz) write(*,*) "Warning: The grid in this file is not cubic! The result may be problematic"
            critmerge=min(dx,dy,dz)*spcmergeratio
        else if (isurftypetmp==11) then
            isurftype=10
            critmerge=min(dx,dy,dz)*spcmergeratio
        end if
        !Set isovalue
        if (isurftype==5.or.isurftype==6) then
            surfisoval=0.5D0 !Hirshfeld and Becke surface the isoval must be 0.5
        else
            write(*,*) "Input the isovalue for defining the isosurface, e.g. 0.5"
            if (isurftype==1) write(*,*) "Hint: Isovalue of 0.001 is commonly used to define molecular vdW surface"
            read(*,*) surfisoval
        end if
        !Set additional parameters
        if (isurftype==1) then !Only remove redundant vertices for electron density isosurface, for other types, elimination method may work poorly
            ifelim=1
            critmerge=grdspc*0.5D0
            nbisec=3
        else if (isurftype==2) then
            ifelim=0
            nbisec=3
        else if (isurftype==5.or.isurftype==6.or.isurftype==10) then
            ifelim=1
            nbisec=0 !Linear interpolation
        end if
        if (.not.allocated(b).and.isurftype/=5.and.isurftype/=6) imapfunc=-1
    else if (isel==2) then
        if (imapfunc==3) a%charge=nucchgbackup !If .chg file is loaded previously, recovery actual nuclear charges
        write(*,*) "Select to real space function to be mapped on the molecular surface"
        write(*,*) "-1 User defined function (determined by ""iuserfunc"" in settings.ini)"
        write(*,*) "0 Certain function loaded from external file"
        write(*,*) "1 Electrostatic potential" !If we have loaded chg file, the atomic charges will be disappear
        write(*,*) "2 Average local ionization energy"
        write(*,*) "3 Electrostatic potential from atomic charge"
        write(*,*) "4 Local electron affinity"
!         write(*,*) "10 Pair density"
        write(*,*) "11 Electron density"
        write(*,*) "12 Sign(lambda2)*rho"
        if (nHirBecatm>0) then
            write(*,*) "20 d_i: distance from the nearest nucleus inside the surface"
            write(*,*) "21 d_e: distance from the nearest nucleus outside the surface"
            write(*,*) "22 d_norm: Normalized contact distance"
        end if
        read(*,*) imapfunc
        
        if (imapfunc==0) then
            ireadextmapval=1
            grdspc=0.2D0
        else
            ireadextmapval=0
            if (imapfunc==1.or.imapfunc==11.or.imapfunc==12) grdspc=0.25D0
            if (imapfunc==2.or.imapfunc==3.or.imapfunc==4.or.imapfunc==-1.or.imapfunc==10) grdspc=0.2D0
        end if
        critmerge=grdspc*spcmergeratio
        if (imapfunc==3) then
            write(*,"(a)") " Please input the path of the .chg file which contains the atomic charges of present system. e.g. c:\t.chg"
            do while(.true.)
                read(*,*) c200tmp
                inquire(file=c200tmp,exist=alive)
                if (alive) exit
                write(*,*) "Couldn't find the file, input again"
            end do
            nucchgbackup=a%charge
            open(10,file=c200tmp,status="old")
            do icen=1,ncenter
                read(10,*) c80tmp,rnouse,rnouse,rnouse,a(icen)%charge
                if (trim(c80tmp)/=trim(a(icen)%name)) write(*,"(' Warning: The name of atom',i7,' in .chg file is not consistent with present system',i7)") icen
            end do
            close(10)
            write(*,*) "The atomic charges have been loaded"
        end if
        if (imapfunc==20.or.imapfunc==21) write(*,*) "NOTE: ALL VALUES OF THIS FUNCTION SHOWN IN LATER STAGE WILL BE BOHR!"
    else if (isel==3) then
        write(*,*) "Input a value (in Bohr)"
        if (imapfunc==0.or.imapfunc==1.or.imapfunc==20.or.imapfunc==21.or.imapfunc==22) write(*,*) "Note: In general 0.25 is enough. For higher accuracy, 0.15~0.20 is recommended"
        if (imapfunc==2.or.imapfunc==3.or.imapfunc==4.or.imapfunc==-1) write(*,*) "Note: In general 0.20 is enough. For higher accuracy, 0.13~0.17 is recommended"
        read(*,*) grdspc
        critmerge=grdspc*spcmergeratio
    else if (isel==4) then
        do while(.true.)
            write(*,*) "0 Return to upper level menu"
            write(*,"(a,f7.4)") " 1 The ratio of vdW radius used to extend spatial region of cubic grids:",vdwmulti
            if (ifelim==0) write(*,*) "2 If eliminate redundant vertices: No"
            if (ifelim==1) write(*,"(a,f6.3,a)") " 2 If elimnate redundant vertices: Yes, criteria:",critmerge," Bohr"
            write(*,"(' 3 Number of bisections before linear interpolation, current:',i5)") nbisec
            read(*,*) isel2
            
            if (isel2==0) then
                exit
            else if (isel2==1) then
                write(*,*) "Input a value"
                write(*,"(a)") "Note: 1.7 is enough for the case of isovalue=0.001, for lower isovalue, a larger value is needed"
                read(*,*) vdwmulti
            else if (isel2==2) then
                if (ifelim==1) then
                    ifelim=0
                else if (ifelim==0) then
                    write(*,"(a)") "Input a value (in Bohr). If the distance between any two vertices is smaller than this value, one of them will be merged to the other."
                    write(*,"(a)") "Hint: 0.5 times grid spacing is recommended in general"
                    read(*,*) critmerge
                    ifelim=1
                end if
            else if (isel2==3) then
                write(*,*) "Perform how many times biections before linear interpolation?"
                write(*,*) "Note: Input 0 means do interpolation directly"
                read(*,*) nbisec
            end if
        end do
    else if (isel==5) then
        write(*,*) "0 Do not load mapped function but directly calculate by Multiwfn"
        write(*,*) "1 Load mapped function at all surface vertices from plain text file"
        write(*,*) "2 Similar to 1, but specific for the case of using cubegen utility of Gaussian"
        write(*,*) "3 Interpolate mapped function from a external cube file"
        read(*,*) ireadextmapval
    end if
end do

!!!!!!!!!!!!!!!!!!! Delete high-lying virtual orbitals to speed up calculation
if (imapfunc/=0.and.imapfunc/=4.and.imapfunc/=20.and.imapfunc/=21.and.imapfunc/=22) call delvirorb(1)

call walltime(iclktime1)
if (isurftype==1.or.isurftype==2.or.isurftype==5.or.isurftype==6) then !Calculate grid data for determining isosurface
    endx=maxval( a(:)%x+vdwmulti*vdwr(a(:)%index) )
    endy=maxval( a(:)%y+vdwmulti*vdwr(a(:)%index) )
    endz=maxval( a(:)%z+vdwmulti*vdwr(a(:)%index) )
    orgx=minval( a(:)%x-vdwmulti*vdwr(a(:)%index) )
    orgy=minval( a(:)%y-vdwmulti*vdwr(a(:)%index) )
    orgz=minval( a(:)%z-vdwmulti*vdwr(a(:)%index) )
    ! This determination may be not find actual min and max range, when different types of atoms are laying on the same plane
    ! endx=maxval(a%x)+vdwmulti*vdwr(a(maxloc(a%x,1))%index)
    ! endy=maxval(a%y)+vdwmulti*vdwr(a(maxloc(a%y,1))%index)
    ! endz=maxval(a%z)+vdwmulti*vdwr(a(maxloc(a%z,1))%index)
    ! orgx=minval(a%x)-vdwmulti*vdwr(a(minloc(a%x,1))%index)
    ! orgy=minval(a%y)-vdwmulti*vdwr(a(minloc(a%y,1))%index)
    ! orgz=minval(a%z)-vdwmulti*vdwr(a(minloc(a%z,1))%index)
    write(*,"(' Spatial range of grid data:')")
    write(*,"(' X is from',f10.4,'  to',f10.4,' Bohr')") orgx,endx
    write(*,"(' Y is from',f10.4,'  to',f10.4,' Bohr')") orgy,endy
    write(*,"(' Z is from',f10.4,'  to',f10.4,' Bohr')") orgz,endz
    xlength=endx-orgx
    ylength=endy-orgy
    zlength=endz-orgz
    dx=grdspc
    dy=grdspc
    dz=grdspc
    nx=nint(xlength/dx)+1
    ny=nint(ylength/dy)+1
    nz=nint(zlength/dz)+1
    write(*,"(' The number of point in x,y,z:',3i6,'  Total:',i10)") nx,ny,nz,nx*ny*nz
    write(*,*)
    
    if (isurftype==1.or.isurftype==2) then
        if (isurftype==1) write(*,*) "Calculating grid data of electron density..."
        if (isurftype==1) write(*,*) "Calculating grid data of the real space function..."
        if (allocated(cubmat)) deallocate(cubmat)
        allocate(cubmat(nx,ny,nz))
        if (isurftype==1) then
            call savecubmat(1,0,1) !Calculate electron density data
        else if (isurftype==2) then
            call savecubmat(ifuncintp,0,iorbsel)
        end if
        !Calculate MICC(molecular intrinsic characteristic contour), however it seems that the result obtained in this manner is not in line with the one in original paper
        ! iuserfunc=10
        ! call savecubmat(100,1,1)
        
    else if (isurftype==5) then !Hirshfeld analysis, calculate weighting distribution of a specific set of atom
        write(*,*) "Hirshfeld analysis requests atomic densities, please select how to obtain them"
        write(*,*) "1 Use build-in sphericalized atomic densities in free-states (recommended)"
        write(*,"(a)") " 2 Provide wavefunction file of involved elements by yourself or invoke Gaussian to automatically calculate them"
        read(*,*) ihirshmode
        if (allocated(cubmat)) deallocate(cubmat)
        if (allocated(cubmattmp)) deallocate(cubmattmp)
        allocate(cubmat(nx,ny,nz),cubmattmp(nx,ny,nz))
        cubmat=0D0
        cubmattmp=0D0
        if (ihirshmode==1) then !Doesn't work well currently, because interpolation of density at long range is problematic by Lagrange method
            do iatm=1,ncenter
                write(*,"(' Processing atom',i6,'(',a,') /',i6)") iatm,a_org(iatm)%name,ncenter_org
                !$OMP PARALLEL DO SHARED(cubmat,cubmattmp) PRIVATE(i,j,k,tmpx,tmpy,tmpz,denstmp) schedule(dynamic) NUM_THREADS(rtNThreads())
                do k=1,nz
                    tmpz=orgz+(k-1)*dz
                    do j=1,ny
                        tmpy=orgy+(j-1)*dy
                        do i=1,nx
                            tmpx=orgx+(i-1)*dx
                            denstmp=calcatmdens(iatm,tmpx,tmpy,tmpz,18)
                            cubmattmp(i,j,k)=cubmattmp(i,j,k)+denstmp
                            if (any(HirBecatm==iatm)) cubmat(i,j,k)=cubmat(i,j,k)+denstmp !Density of specified fragment
                        end do
                    end do
                end do
                !$OMP END PARALLEL DO
            end do
        else if (ihirshmode==2) then
            call setpromol
            do iatm=1,ncenter
                call dealloall
                write(*,"(' Processing atom',i6,'(',a,') /',i6)") iatm,a_org(iatm)%name,ncenter_org
                call readwfn(custommapname(iatm),1)
                !$OMP PARALLEL DO SHARED(cubmat,cubmattmp) PRIVATE(i,j,k,tmpx,tmpy,tmpz,denstmp) schedule(dynamic) NUM_THREADS(rtNThreads())
                do k=1,nz
                    tmpz=orgz+(k-1)*dz
                    do j=1,ny
                        tmpy=orgy+(j-1)*dy
                        do i=1,nx
                            tmpx=orgx+(i-1)*dx
                            denstmp=fdens(tmpx,tmpy,tmpz)
                            cubmattmp(i,j,k)=cubmattmp(i,j,k)+denstmp
                            if (any(HirBecatm==iatm)) cubmat(i,j,k)=cubmat(i,j,k)+denstmp !Density of specified fragment
                        end do
                    end do
                end do
                !$OMP END PARALLEL DO
            end do
            call dealloall
            write(*,"(' Reloading ',a)") trim(firstfilename)
            call readinfile(firstfilename,1) !Retrieve to first loaded file(whole molecule)
        end if
        cubmat=cubmat/cubmattmp

    else if (isurftype==6) then !Becke analysis, calculate weighting distribution of a specific set of atom
        if (allocated(cubmat)) deallocate(cubmat)
        allocate(cubmat(nx,ny,nz))
        cubmat=0D0
        ifinish=0
        !We calculate Becke weight for all atoms, but only summing up the value of we selected atoms to cubmat
        !$OMP PARALLEL DO SHARED(cubmat) PRIVATE(i,j,k,tmpx,rnowx,rnowy,rnowz,smat,ii,jj,tmprarr,rmiu,chi,uij,aij,tmps,iter,Pvec) schedule(dynamic) NUM_THREADS(rtNThreads())
        do k=1,nz
            rnowz=orgz+(k-1)*dz
            do j=1,ny
                rnowy=orgy+(j-1)*dy
                do i=1,nx
                    rnowx=orgx+(i-1)*dx
                    !Calculate Becke weight with modified CSD radii, by using Eq. 11,21,13,22 in Becke's paper (JCP 88,15)
                    smat=1.0D0
                    do ii=1,ncenter
                        tmprarr(ii)=dsqrt( (rnowx-a(ii)%x)**2+(rnowy-a(ii)%y)**2+(rnowz-a(ii)%z)**2 )
                    end do
                    do ii=1,ncenter
                        do jj=1,ncenter
                            if (ii==jj) cycle
                            rmiu=(tmprarr(ii)-tmprarr(jj))/distmat(ii,jj)
                            chi=vdwr_tianlu(a(ii)%index)/vdwr_tianlu(a(jj)%index)
                            uij=(chi-1)/(chi+1)
                            aij=uij/(uij**2-1)
                            if (aij>0.5D0) aij=0.5D0
                            if (aij<-0.5D0) aij=-0.5D0
                            rmiu=rmiu+aij*(1-rmiu**2)
                            tmps=rmiu
                            do iter=1,3
                                tmps=1.5D0*(tmps)-0.5D0*(tmps)**3
                            end do
                            smat(ii,jj)=0.5D0*(1-tmps)
                        end do
                    end do
                    Pvec=1.0D0
                    do ii=1,ncenter
                        Pvec=Pvec*smat(:,ii)
                    end do
                    Pvec=Pvec/sum(Pvec) !Normalized Pvec, Pvec contain partition weight of each atom in current point
                    cubmat(i,j,k)=cubmat(i,j,k)+sum(Pvec(HirBecatm(:)))
                end do
            end do
            !$OMP CRITICAL
            ifinish=ifinish+1
            write(*,"(' Finished:',i4,'   /',i4)") ifinish,nz
            !$OMP end CRITICAL
        end do
        !$OMP END PARALLEL DO
    end if
end if

allocate(corpos(nx,ny,nz),ifbndcub(nx-1,ny-1,nz-1))
write(*,*)

!Initialize
numcubx=nx-1
numcuby=ny-1
numcubz=nz-1
ifbndcub=.false.
nsurvtx=0
nsurtri=0
tetravol1=0D0 !Volume in generated surface tetrahedron of type 1,2,3, will be accumulated in polygonization process
tetravol2=0D0
tetravol3=0D0
allocate(abs2suridx(nx,ny,nz))
abs2suridx=0
surlocminidx=0
surlocmaxidx=0

!Which corners(nx*ny*nz in total) is in internal (=0) or external (=1) of specified isosurface
do ix=1,nx
    do iy=1,ny
        do iz=1,nz
            if (cubmat(ix,iy,iz)>=surfisoval) then !internal corner
                corpos(ix,iy,iz)=0
            else !external corner
                corpos(ix,iy,iz)=1
            end if
        end do
    end do
end do

nintcub=0
nextcub=0
nbndcub=0
icorind=0
do ix=1,numcubx
    do iy=1,numcuby
        do iz=1,numcubz
            icubtest=corpos(ix,iy,iz)+corpos(ix+1,iy,iz)+corpos(ix,iy+1,iz)+corpos(ix,iy,iz+1)+corpos(ix+1,iy+1,iz)+corpos(ix,iy+1,iz+1)+corpos(ix+1,iy,iz+1)+corpos(ix+1,iy+1,iz+1)
            if (icubtest==0) then !internal cube
                nintcub=nintcub+1
            else if (icubtest==8) then !external cube
                nextcub=nextcub+1
            else
                ifbndcub(ix,iy,iz)=.true.
                nbndcub=nbndcub+1 !boundary cube
                !Give each corner of boundary cube a unique index
                if (abs2suridx(ix,iy,iz)==0) then
                    icorind=icorind+1
                    abs2suridx(ix,iy,iz)=icorind
                end if
                if (abs2suridx(ix+1,iy,iz)==0) then
                    icorind=icorind+1
                    abs2suridx(ix+1,iy,iz)=icorind
                end if
                if (abs2suridx(ix,iy+1,iz)==0) then
                    icorind=icorind+1
                    abs2suridx(ix,iy+1,iz)=icorind
                end if
                if (abs2suridx(ix,iy,iz+1)==0) then
                    icorind=icorind+1
                    abs2suridx(ix,iy,iz+1)=icorind
                end if
                if (abs2suridx(ix+1,iy+1,iz)==0) then
                    icorind=icorind+1
                    abs2suridx(ix+1,iy+1,iz)=icorind
                end if
                if (abs2suridx(ix,iy+1,iz+1)==0) then
                    icorind=icorind+1
                    abs2suridx(ix,iy+1,iz+1)=icorind
                end if
                if (abs2suridx(ix+1,iy,iz+1)==0) then
                    icorind=icorind+1
                    abs2suridx(ix+1,iy,iz+1)=icorind
                end if
                if (abs2suridx(ix+1,iy+1,iz+1)==0) then
                    icorind=icorind+1
                    abs2suridx(ix+1,iy+1,iz+1)=icorind
                end if
            end if
        end do
    end do
end do
numsurfcubcor=count(abs2suridx/=0)
deallocate(corpos)
write(*,"(' The number of boundary cubes:        ',i10)") nbndcub  !debug info
write(*,"(' The number of corners of boundary cubes:',i10)") numsurfcubcor  !debug info

nassumsurfvtx=7*nbndcub !I assume can generate up to 7*nbndcub surface vertices. In principle, each cube can generate up to 19 vertices, since there are 19 edges in main-axis decomposition
nassumfacet=12*nbndcub !each boundary cube can generate up to 2*6=12 facets in principle
allocate(survtx(nassumsurfvtx)) !, so in principle should be 19*nbndcub... 
allocate(surtriang(nassumfacet))
allocate(surfcor2vtx(numsurfcubcor,14),surcor2vtxpos(numsurfcubcor)) !14 is max number that each corner can link other corner for main-axis decomposition of cube
surfcor2vtx(:,:)%athcor=0
surcor2vtxpos=0
allocate(vtxconn(nassumsurfvtx,50),vtxconnpos(nassumsurfvtx)) !I assume that each surface vertex can connect up to 50 other surface vertices
vtxconn=0
vtxconnpos=0

CALL CPU_TIME(time_begin)
write(*,*) "Generating isosurface by Marching Tetrahedra algorithm, please wait..."
!DO NOT use parallel method in this step, otherwise the numbering will be disordered, and we can't read external mapped function data, since indices cannot be consistent
do ix=1,numcubx
    do iy=1,numcuby
        do iz=1,numcubz
            if (ifbndcub(ix,iy,iz) .eqv. .true.) then !Numbering of cube corner is identical to figure 3 of WFA original paper
                call marchtetra(ix,iy,iz)
!                 call marchcube(ix,iy,iz)
            end if
        end do
    end do
!     write(*,"('Finished:',f5.1,'%')") dfloat(ix)/numcubx*100D0
end do
CALL CPU_TIME(time_end)
deallocate(ifbndcub,abs2suridx,surfcor2vtx,surcor2vtxpos) !Discard arrays used in polygonization
write(*,"(' Polygonization took up CPU time',f12.2,' seconds')") time_end-time_begin

nsuredge=sum(vtxconnpos(1:nsurvtx))/2
write(*,"(' The number of surface vertices (V): ',i10)") nsurvtx
write(*,"(' The number of surface edges (E):    ',i10)") nsuredge
write(*,"(' The number of triangular facets (F):',i10)") nsurtri
iVFEtest=nsurvtx+nsurtri-nsuredge-2
if (imapfunc/=20.and.imapfunc/=21.and.imapfunc/=22) then !The surface of Hirshfeld/Becke is usually open and thus must break this rule, so don't check this
    write(*,"(' V+F-E-2=',i10)") iVFEtest  !debug info
    if (iVFEtest/=0) write(*,"(' Warning: V+F-E-2=0 is violated! Probably grid spacing is too large or the isosurface is not closed')")
end if
volcub=dx*dy*dz
totintvol=nintcub*volcub !Volume of internal cubes
totvol=tetravol0+tetravol1+tetravol2+tetravol3+totintvol !Volume of boundary tetrahedra
! write(*,"('Vol. of type-0,1,2,3:',4f11.7,' Ang.^3')") tetravol0*b2a**3,tetravol1*b2a**3,tetravol2*b2a**3,tetravol3*b2a**3
write(*,"(' Volume enclosed by the isosurface:',f12.5,' Bohr^3  (',f10.5,' Angstrom^3)')") totvol,totvol*b2a**3
write(*,*) "Among all surface vertices:"
write(*,"(' Min-X:',f12.4,'  Max-X:',f10.4,' Angstrom')") minval(survtx(1:nsurvtx)%x)*b2a,maxval(survtx(1:nsurvtx)%x)*b2a
write(*,"(' Min-Y:',f12.4,'  Max-Y:',f10.4,' Angstrom')") minval(survtx(1:nsurvtx)%y)*b2a,maxval(survtx(1:nsurvtx)%y)*b2a
write(*,"(' Min-Z:',f12.4,'  Max-Z:',f10.4,' Angstrom')") minval(survtx(1:nsurvtx)%z)*b2a,maxval(survtx(1:nsurvtx)%z)*b2a
write(*,*)

!Merge redundant vertices
allocate(elimvtx(nsurvtx),elimtri(nint(nsurtri*1.1D0)),mergerelat(nsurvtx)) !The reason for *1.1D0, is because during elimination of vertices with three connections, we need add new triangle
elimvtx=0 !If it is 1, means this surface vertex has been eliminated
elimtri=0
forall(i=1:nsurvtx) mergerelat(i)=i !If mergerelat(i)=j, means surface vertex i has been absorbed into vertex j
if (ifelim==1) then
    write(*,*) "Eliminating redundant surface vertices..."
    do ivt1=1,nsurvtx
        if (elimvtx(ivt1)==1) cycle
        vt1x=survtx(ivt1)%x
        vt1y=survtx(ivt1)%y
        vt1z=survtx(ivt1)%z
        do j=1,vtxconnpos(ivt1)
            ivt2=vtxconn(ivt1,j)
            if (elimvtx(ivt2)==1) cycle
            vt2x=survtx(ivt2)%x
            vt2y=survtx(ivt2)%y
            vt2z=survtx(ivt2)%z
            dist2=(vt1x-vt2x)**2+(vt1y-vt2y)**2+(vt1z-vt2z)**2
            if (dist2<critmerge**2) then !merge vt2 to vt1
                !Scan and update connectivity of neighbour of vt2 (except vt1), also add new connections to vt1 from vt2
                do k=1,vtxconnpos(ivt2)
                    inei=vtxconn(ivt2,k)
                    if (elimvtx(inei)==1) cycle
                    if (inei==ivt1) cycle
                    ihasvt1=0
                    do itmp=1,vtxconnpos(inei)
                        ineinei=vtxconn(inei,itmp)
                        if (elimvtx(ineinei)==1) cycle
                        if (ineinei==ivt1) ihasvt1=1
                    end do
                    if (ihasvt1==0) then
                        vtxconnpos(inei)=vtxconnpos(inei)+1 !if this neighbour vertex of vt2 hasn't connected to vt1, add the connection
                        vtxconn(inei,vtxconnpos(inei))=ivt1
                        vtxconnpos(ivt1)=vtxconnpos(ivt1)+1
                        vtxconn(ivt1,vtxconnpos(ivt1))=inei
                    end if
                end do
                mergerelat(ivt2)=ivt1
                elimvtx(ivt2)=1 !eliminate ivt2
                !Change coordinate of ivtx1 to the average of ivtx1 and ivtx2
                survtx(ivt1)%x=(vt1x+vt2x)/2D0
                survtx(ivt1)%y=(vt1y+vt2y)/2D0
                survtx(ivt1)%z=(vt1z+vt2z)/2D0
                vt1x=survtx(ivt1)%x
                vt1y=survtx(ivt1)%y
                vt1z=survtx(ivt1)%z
            end if
        end do
    end do
!     do i=1,nsurvtx
!         write(1,*) i,mergerelat(i),elimvtx(i)
!     end do
    
    !Update vertices of facets
    do itri=1,nsurtri
        do itmp=1,3
            do while(.true.) !Update vertex indices of each facet according to mergence relation map
                if (elimvtx(surtriang(itri)%idx(itmp))==1) then
                    surtriang(itri)%idx(itmp)=mergerelat(surtriang(itri)%idx(itmp))
                else
                    exit
                end if
            end do
        end do
outdo:    do itmp=1,3 !If a facet has two identical vertices, then delete this facet
            do jtmp=itmp+1,3
                if (surtriang(itri)%idx(itmp)==surtriang(itri)%idx(jtmp)) then
                    elimtri(itri)=1
                    exit outdo
                end if
            end do
        end do outdo
    end do
!     write(*,"('Eliminated surface vertices in merge step:',i10)") count(elimvtx==1)
!     write(*,"('Eliminated surface facets in merge step:',i10)") count(elimtri==1)
    
    !Eliminate surface vertices with two neighbours
    nelimtwoconvtx=0
    do ivt1=1,nsurvtx
        if (elimvtx(ivt1)==1) cycle
        nlink=0
        do j=1,vtxconnpos(ivt1)
            ivt2=vtxconn(ivt1,j)
            if (elimvtx(ivt2)==0) nlink=nlink+1
            if (nlink>2) exit
        end do
        if (nlink==2) then
            nelimtwoconvtx=nelimtwoconvtx+1
            elimvtx(ivt1)=1
        end if
    end do
    !Eliminate surface vertices with three neighbours
    nelimthreeconvtx=0
    do ivt1=1,nsurvtx
        if (elimvtx(ivt1)==1) cycle
        nlink=0
        do j=1,vtxconnpos(ivt1)
            ivt2=vtxconn(ivt1,j)
            if (elimvtx(ivt2)==0) then
                nlink=nlink+1
                if (nlink>3) exit
                tmpintarr3(nlink)=ivt2
            end if
        end do
        if (nlink==3) then
            nelimthreeconvtx=nelimthreeconvtx+1
            elimvtx(ivt1)=1
            nsurtri=nsurtri+1
            surtriang(nsurtri)%idx(:)=tmpintarr3(:)
            elimtri(nsurtri)=0
        end if
    end do
    !Remove the triangles which has at least one vertex has been discarded
!     iz=0
    do itri=1,nsurtri
        if (elimtri(itri)==1) cycle
        do itmp=1,3
            if (elimvtx(surtriang(itri)%idx(itmp))==1) then
                elimtri(itri)=1
!                 iz=iz+1
!                 write(*,"(5i10)") iz,surtriang(itri)%idx(itmp),surtriang(itri)%idx(:)
                exit
            end if
        end do
    end do
!     write(*,"('Eliminated surface vertices with two connections:',i10)") nelimtwoconvtx
!     write(*,"('Eliminated surface vertices with three connections:',i10)") nelimthreeconvtx
    
    !Elimination has finished, count how many surface vertices, edges and facets currently
    ncurrtri=count(elimtri(1:nsurtri)==0)
    ncurrvtx=count(elimvtx(1:nsurvtx)==0)
    ncurredge=0
    do ivtx=1,nsurvtx
        if (elimvtx(ivtx)==1) cycle
        do itmp=1,vtxconnpos(ivtx)
            if (elimvtx(vtxconn(ivtx,itmp))==1) cycle
            ncurredge=ncurredge+1
        end do
    end do
    ncurredge=ncurredge/2
    write(*,"(' After elimination, V=',i9,',  E=',i9,',  F=',i9,',  V+F-E-2=',i6)") ncurrvtx,ncurredge,ncurrtri,ncurrvtx+ncurrtri-ncurredge-2
!     write(*,*) maxval(vtxconnpos(1:nsurvtx))
else
    ncurrtri=nsurtri !If not do vertices elimnation
    ncurrvtx=nsurvtx
    ncurredge=nsuredge
end if


!Check deviation of surface vertices from specified isovalue
! valtotdev=0D0
! valRMSD=0D0
! valmaxdev=0D0
! num=0
! do i=1,nsurvtx
!     if (elimvtx(i)==1) cycle
!     num=num+1
!     devtmp=abs(fdens(survtx(i)%x,survtx(i)%y,survtx(i)%z)-surfisoval)
!     valtotdev=valtotdev+devtmp
!     valRMSD=valRMSD+devtmp**2
!     if (devtmp>valmaxdev) valmaxdev=devtmp
! end do
! valtotdev=valtotdev/num
! valRMSD=dsqrt(valRMSD/num)
! write(*,"('Average deviation from specified isovalue:',f20.16)") valtotdev
! write(*,"('RMSD from specified isovalue:',f20.16)") valRMSD
! write(*,"('Maximum deviation from specified isovalue:',f20.16)") valmaxdev
! read (*,*)

!Calculate isosurface area and geometry center of each facet
surfareaall=0D0
triareamin=100D0
triareamax=0D0
do icyc=1,nsurtri
    if (elimtri(icyc)==1) cycle
    idx1=surtriang(icyc)%idx(1)
    idx2=surtriang(icyc)%idx(2)
    idx3=surtriang(icyc)%idx(3)
    xidx1=survtx(idx1)%x
    yidx1=survtx(idx1)%y
    zidx1=survtx(idx1)%z
    xidx2=survtx(idx2)%x
    yidx2=survtx(idx2)%y
    zidx2=survtx(idx2)%z
    xidx3=survtx(idx3)%x
    yidx3=survtx(idx3)%y
    zidx3=survtx(idx3)%z
    surtriang(icyc)%area=gettriangarea( xidx1,yidx1,zidx1,xidx2,yidx2,zidx2,xidx3,yidx3,zidx3 )
    surfareaall=surfareaall+surtriang(icyc)%area
    if (surtriang(icyc)%area>triareamax) triareamax=surtriang(icyc)%area
    if (surtriang(icyc)%area<triareamin) triareamin=surtriang(icyc)%area
end do
! write(*,"('Minimum and maximum facet area:',D18.10,f15.10,' Angstrom^2')") triareamin*b2a*b2a,triareamax*b2a*b2a  !debug info
write(*,"(' Isosurface area:',f12.5,' Bohr^2  (',f10.5,' Angstrom^2)')") surfareaall,surfareaall*b2a*b2a


!Calculate mapped function at each surface vertex, save to survtx%value
if (ireadextmapval==0) then
    write(*,*)
    iprogstep=5
    if (imapfunc==1.or.imapfunc==3) then
        write(*,*) "Calculating electrostatic potential at surface vertices, please wait patiently"
        if (nprims>500) then
            write(*,"(a)") " NOTE: This system is relatively large, electrostatic potential calculation will consume very long time!!!"
            iprogstep=1
        end if
    else if (imapfunc==2) then
        write(*,*) "Calculating average local ionization energy at surface vertices..."
    else if (imapfunc==4) then
        write(*,*) "Calculating local electron affinity at surface vertices..."
    else if (imapfunc==-1) then
        write(*,*) "Calculating user defined real space function at surface vertices..."
    else
        write(*,*) "Calculating mapped function value at surface vertices..."
    end if
    ii=0
    iprog=0
    CALL CPU_TIME(time_begin)
    !$OMP PARALLEL DO SHARED(survtx,indsurfmax,indsurfmin,iprog) PRIVATE(icyc) schedule(dynamic) NUM_THREADS(rtNThreads())
    do icyc=1,nsurvtx
        if (elimvtx(icyc)==1) cycle
        if (imapfunc==1) then
            survtx(icyc)%value=totesp(survtx(icyc)%x,survtx(icyc)%y,survtx(icyc)%z)
        else if (imapfunc==2) then
            survtx(icyc)%value=avglocion(survtx(icyc)%x,survtx(icyc)%y,survtx(icyc)%z)
        else if (imapfunc==3) then
            survtx(icyc)%value=nucesp(survtx(icyc)%x,survtx(icyc)%y,survtx(icyc)%z)
        else if (imapfunc==4) then
            survtx(icyc)%value=loceleaff(survtx(icyc)%x,survtx(icyc)%y,survtx(icyc)%z)
        else if (imapfunc==-1) then
            survtx(icyc)%value=userfunc(survtx(icyc)%x,survtx(icyc)%y,survtx(icyc)%z)
        else if (imapfunc==10) then
            survtx(icyc)%value=pairfunc(survtx(icyc)%x,survtx(icyc)%y,survtx(icyc)%z)  !Calculate pair density, not for normal users
        else if (imapfunc==11) then
            survtx(icyc)%value=fdens(survtx(icyc)%x,survtx(icyc)%y,survtx(icyc)%z)  !Calculate electron density
        else if (imapfunc==12) then
            survtx(icyc)%value=signlambda2rho(survtx(icyc)%x,survtx(icyc)%y,survtx(icyc)%z)  !Calculate Sign(lambda2)*rho
        else if (imapfunc==20) then
            survtx(icyc)%value=surfana_di(survtx(icyc)%x,survtx(icyc)%y,survtx(icyc)%z,nHirBecatm,HirBecatm)
        else if (imapfunc==21) then
            survtx(icyc)%value=surfana_de(survtx(icyc)%x,survtx(icyc)%y,survtx(icyc)%z,nHirBecatm,HirBecatm)
        else if (imapfunc==22) then
            survtx(icyc)%value=surfana_norm(survtx(icyc)%x,survtx(icyc)%y,survtx(icyc)%z,nHirBecatm,HirBecatm)
        end if
        iprog=iprog+1
        progress=dfloat(iprog)/ncurrvtx*100
        if (progress>ii) then
            ii=ii+iprogstep
            if (imapfunc/=20.and.imapfunc/=21.and.imapfunc/=22)  write(*,"(' Finished ',f6.1,'%')") progress !Calculate d is very fast, don't show prog
        end if
    end do
    !$OMP END PARALLEL DO
    if (imapfunc/=20.and.imapfunc/=21.and.imapfunc/=22)  write(*,"(' Finished ',f6.1,'%')") 100D0
    CALL CPU_TIME(time_end)
    write(*,"(' Calculation took up CPU time',f12.2,' seconds')") time_end-time_begin
else if (ireadextmapval==1) then
    open(10,file="surfptpos.txt",status="replace")
    write(10,"(i10)") ncurrvtx
    do icyc=1,nsurvtx
        if (elimvtx(icyc)==1) cycle
        write(10,"(3f13.7)") survtx(icyc)%x,survtx(icyc)%y,survtx(icyc)%z
    end do
    close(10)
    write(*,"(/,a)") "The X/Y/Z coordinates of the surface vertices have been exported to surfptpos.txt in current folder, unit is in Bohr. &
    Now you can use your favourite program to calculate mapped function values at this points"
    write(*,"(/,a)") " Now input the path of the file containing calculated mapped function values. The data will be read in free format. The first row is the number of points, &
    in following lines the 1/2/3/4 column should correspond to X,Y,Z and mapped function value, respectively. All units must be in a.u."
    do while(.true.)
        read(*,"(a)") c200tmp
        inquire(file=c200tmp,exist=alive)
        if (alive) exit
        write(*,*) "File cannot be found, input again"
    end do
    open(10,file=c200tmp,status="old")
    read(10,*) npttmp !How many points in the file
    if (npttmp/=ncurrvtx) then
        write(*,"(a,i10,a,i10)") " Error: The number of points in this file",npttmp," is inconsistent with the number of vertices on the molecular surface",ncurrvtx
        exit surfanaloop
    end if
    do icyc=1,nsurvtx
        if (elimvtx(icyc)==1) cycle
        read(10,*) rnouse,rnouse,rnouse,survtx(icyc)%value
    end do
    close(10)
    write(*,*) "Loading data finished!"
else if (ireadextmapval==2) then
    open(10,file="cubegenpt.txt",status="replace")
    do icyc=1,nsurvtx
        if (elimvtx(icyc)==1) cycle
        write(10,"(3f13.7)") survtx(icyc)%x*b2a,survtx(icyc)%y*b2a,survtx(icyc)%z*b2a
    end do
    close(10)
    write(*,"(/,a)") " The coordinate of the surface vertices has been outputted to cubegenpt.txt"
    write(*,"(a,/)") " Use for example ""cubegen 0 potential CNT.fch result.cub -5 h < cubegenpt.txt"" to get cubegen output file"
    write(*,*) "Now input the file name of cubegen output file, e.g. d:\g09\result.cub"
    do while(.true.)
        read(*,"(a)") c200tmp
        inquire(file=c200tmp,exist=alive)
        if (alive) exit
        write(*,*) "File cannot be found, input again"
    end do
    open(10,file=c200tmp,status="old")
    do iskip=1,6+ncenter
        read(10,*)
    end do
    do icyc=1,nsurvtx
        if (elimvtx(icyc)==1) cycle
        read(10,*) rnouse,rnouse,rnouse,survtx(icyc)%value
    end do
    close(10)
    write(*,*) "Loading data finished!"
else if (ireadextmapval==3) then !Will calculate mapped function by interpolating from external cube file (cubmattmp)
    write(*,*)
    write(*,*) "Outputting template cube file..."
    open(10,file="template.cub",status="replace")
    call outcube(cubmat,nx,ny,nz,orgx,orgy,orgz,dx,dy,dz,10)
    close(10)
    write(*,"(/,a)") " The template cube file has been outputted to template.cub in current folder"
    write(*,"(a)") " Now input the name of the cube file representing mapped function, e.g. c:\t.cub. &
    The grid setting in this cube file must be exactly identical to the template cube file"
    do while(.true.)
        read(*,"(a)") c200tmp
        inquire(file=c200tmp,exist=alive)
        if (alive) exit
        write(*,*) "File cannot be found, input again"
    end do
    call readcubetmp(c200tmp,inconsis)
    write(*,*) "Loading data finished!"
    if (inconsis==1) then
        write(*,"(a)") " Warning: The grid setting in the cube file you inputted is not identical to the template cube file! The analysis result may be meaningless!"
        read (*,*)
    end if
    do icyc=1,nsurvtx
        if (elimvtx(icyc)==1) cycle
        survtx(icyc)%value=linintp3d(survtx(icyc)%x,survtx(icyc)%y,survtx(icyc)%z,2)
    end do
end if
!Find maximum and minimum value of surface vertex
itime=0
do icyc=1,nsurvtx
    if (elimvtx(icyc)==1) cycle
    itime=itime+1
    if (itime==1) then
        indsurfmin=icyc
        indsurfmax=icyc
    end if
    if (survtx(icyc)%value>survtx(indsurfmax)%value) indsurfmax=icyc
    if (survtx(icyc)%value<survtx(indsurfmin)%value) indsurfmin=icyc
end do
!Approximately evaluate mapped function value at each facet, use average of values at three verices as the value of facet 
do icyc=1,nsurtri
    if (elimtri(icyc)==1) cycle
    idx1=surtriang(icyc)%idx(1)
    idx2=surtriang(icyc)%idx(2)
    idx3=surtriang(icyc)%idx(3)
    surtriang(icyc)%value=(survtx(idx1)%value+survtx(idx2)%value+survtx(idx3)%value)/3D0
end do
if (size(survtx)>0) then !The number of vertex is zero for empty grid data
    write(*,"(' Global surface minimum:',f10.6,' a.u. at',3f11.6,' Ang')") survtx(indsurfmin)%value,survtx(indsurfmin)%x*b2a,survtx(indsurfmin)%y*b2a,survtx(indsurfmin)%z*b2a
    write(*,"(' Global surface maximum:',f10.6,' a.u. at',3f11.6,' Ang')") survtx(indsurfmax)%value,survtx(indsurfmax)%x*b2a,survtx(indsurfmax)%y*b2a,survtx(indsurfmax)%z*b2a
end if
write(*,*)


! Find local minimum of mapped function on surface
nsurlocmin=0
cmin: do ivtx=1,nsurvtx
    if (elimvtx(ivtx)==1) cycle
    vali=survtx(ivtx)%value
    do jtmp=1,vtxconnpos(ivtx)
        jvtx=vtxconn(ivtx,jtmp)
        if (elimvtx(jvtx)==1) cycle
        if (vali>=survtx(jvtx)%value) cycle cmin
        do ktmp=1,vtxconnpos(jvtx)
            kvtx=vtxconn(jvtx,ktmp)
            if (elimvtx(kvtx)==1) cycle
            if (kvtx==ivtx) cycle
            if (vali>=survtx(kvtx)%value) cycle cmin
        end do
    end do
    nsurlocmin=nsurlocmin+1
    surlocminidx(nsurlocmin)=ivtx
end do cmin
write(*,"(a,i6)") " The number of surface minima:",nsurlocmin
if (imapfunc==1.or.imapfunc==2.or.imapfunc==3.or.imapfunc==4) then
    write(*,*) "  #       a.u.         eV      kcal/mol           X/Y/Z coordinate(Angstrom)"
else
    write(*,*) "  #             Value           X/Y/Z coordinate(Angstrom)"
end if
do i=1,nsurlocmin
    idx=surlocminidx(i)
    char1tmp=' '
    if (idx==indsurfmin) char1tmp='*' !Mark minimum term by asterisk
    if (imapfunc==1.or.imapfunc==2.or.imapfunc==3.or.imapfunc==4) then
        write(*,"(a,i5,f12.8,f12.6,f12.6,4x,3f11.6)") char1tmp,i,survtx(idx)%value,survtx(idx)%value*au2eV,survtx(idx)%value*au2kcal,survtx(idx)%x*b2a,survtx(idx)%y*b2a,survtx(idx)%z*b2a
    else
        write(*,"(a,i5,f18.6,4x,3f11.6)") char1tmp,i,survtx(idx)%value,survtx(idx)%x*b2a,survtx(idx)%y*b2a,survtx(idx)%z*b2a
    end if
end do
write(*,*)

! Find local maximum of mapped function on surface
nsurlocmax=0
cmax: do ivtx=1,nsurvtx
    if (elimvtx(ivtx)==1) cycle
    vali=survtx(ivtx)%value
    do jtmp=1,vtxconnpos(ivtx)
        jvtx=vtxconn(ivtx,jtmp)
        if (elimvtx(jvtx)==1) cycle
        if (vali<=survtx(jvtx)%value) cycle cmax
        do ktmp=1,vtxconnpos(jvtx)
            kvtx=vtxconn(jvtx,ktmp)
            if (elimvtx(kvtx)==1) cycle
            if (kvtx==ivtx) cycle
            if (vali<=survtx(kvtx)%value) cycle cmax
        end do
    end do
    nsurlocmax=nsurlocmax+1
    surlocmaxidx(nsurlocmax)=ivtx
end do cmax
write(*,"(a,i6)") " The number of surface maxima:",nsurlocmax
if (imapfunc==1.or.imapfunc==2.or.imapfunc==3.or.imapfunc==4) then
    write(*,*) "  #       a.u.         eV      kcal/mol           X/Y/Z coordinate(Angstrom)"
else
    write(*,*) "  #             Value           X/Y/Z coordinate(Angstrom)"
end if
do i=1,nsurlocmax
    idx=surlocmaxidx(i) !Convert to absolute surface vertice index
    char1tmp=' '
    if (idx==indsurfmax) char1tmp='*' !Mark maximum term by asterisk
    if (imapfunc==1.or.imapfunc==2.or.imapfunc==3.or.imapfunc==4) then
        write(*,"(a,i5,f12.8,f12.6,f12.6,4x,3f11.6)") char1tmp,i,survtx(idx)%value,survtx(idx)%value*au2eV,survtx(idx)%value*au2kcal,survtx(idx)%x*b2a,survtx(idx)%y*b2a,survtx(idx)%z*b2a
    else
        write(*,"(a,i5,f18.6,4x,3f11.6)") char1tmp,i,survtx(idx)%value,survtx(idx)%x*b2a,survtx(idx)%y*b2a,survtx(idx)%z*b2a
    end if
end do


!============ Perform statistic
write(*,*)
write(*,*) "      ================= Summary of surface analysis ================="
write(*,*)
totvalang3=totvol*b2a**3
totmass=sum(atmwei(a%index))
densmass=totmass/avogacst*(1D24/totvalang3)
write(*,"(' Volume:',f12.5,' Bohr^3  (',f10.5,' Angstrom^3)')") totvol,totvalang3
write(*,"(' Estimated density according to mass and volume:',f10.4,' g/cm^3')") densmass
if (imapfunc==1.or.imapfunc==3) then !ESP or ESP from atomic charges
    write(*,"(' Minimal value:',f13.5,' kcal/mol   Maximal value:',f13.5,' kcal/mol')") survtx(indsurfmin)%value*au2kcal,survtx(indsurfmax)%value*au2kcal
else if (imapfunc==2.or.imapfunc==4) then !ALIE or LEA
    write(*,"(' Minimal value:',f13.5,' eV,   Maximal value:',f13.5,' eV')") survtx(indsurfmin)%value*au2eV,survtx(indsurfmax)%value*au2eV
else
    write(*,"(' Minimal value:',f16.8,'    Maximal value:',f16.8)") survtx(indsurfmin)%value,survtx(indsurfmax)%value
end if
surfareapos=0D0
surfareaneg=0D0
do i=1,nsurtri
    if (elimtri(i)==1) cycle
    if (surtriang(i)%value>=0) then
        surfareapos=surfareapos+surtriang(i)%area
    else
        surfareaneg=surfareaneg+surtriang(i)%area
    end if
end do
write(*,"(' Overall surface area:      ',f12.5,' Bohr^2  (',f10.5,' Angstrom^2)')") surfareaall,surfareaall*b2a*b2a
write(*,"(' Positive surface area:     ',f12.5,' Bohr^2  (',f10.5,' Angstrom^2)')") surfareapos,surfareapos*b2a*b2a
write(*,"(' Negative surface area:     ',f12.5,' Bohr^2  (',f10.5,' Angstrom^2)')") surfareaneg,surfareaneg*b2a*b2a

sumposval=0D0
sumnegval=0D0
do i=1,nsurtri
    if (elimtri(i)==1) cycle
    valtmp=surtriang(i)%value
    if (valtmp>=0D0) then
        sumposval=sumposval+valtmp*surtriang(i)%area
    else
        sumnegval=sumnegval+valtmp*surtriang(i)%area
    end if
end do
sumallval=sumposval+sumnegval
avgpos=sumposval/surfareapos
avgneg=sumnegval/surfareaneg
avgall=sumallval/surfareaall
if (imapfunc==1.or.imapfunc==3) then !ESP or ESP from atomic charges
    write(*,"(' Overall average value: ',f13.8,' a.u. (',f13.5,' kcal/mol)')") avgall,avgall*au2kcal
    write(*,"(' Positive average value:',f13.8,' a.u. (',f13.5,' kcal/mol)')") avgpos,avgpos*au2kcal
    write(*,"(' Negative average value:',f13.8,' a.u. (',f13.5,' kcal/mol)')") avgneg,avgneg*au2kcal
else if (imapfunc==2) then !ALIE is always positive
    write(*,"(' Average value: ',f13.8,' a.u. (',f13.5,' eV,',f13.5,' kcal/mol)')") avgpos,avgpos*au2eV,avgpos*au2kcal
else if (imapfunc==4) then !LEA
    write(*,"(' Overall average value: ',f11.7,' a.u. (',f10.5,' eV,',f13.5,' kcal/mol)')") avgall,avgall*au2eV,avgall*au2kcal
    write(*,"(' Positive average value:',f11.7,' a.u. (',f10.5,' eV,',f13.5,' kcal/mol)')") avgpos,avgpos*au2eV,avgpos*au2kcal
    write(*,"(' Negative average value:',f11.7,' a.u. (',f10.5,' eV,',f13.5,' kcal/mol)')") avgneg,avgneg*au2eV,avgneg*au2kcal
else
    write(*,"(' Overall average value: ',f20.8)") avgall
    write(*,"(' Positive average value:',f20.8)") avgpos
    write(*,"(' Negative average value:',f20.8)") avgneg
end if

chgsep=0D0
varipos=0D0
varineg=0D0
do i=1,nsurtri
    if (elimtri(i)==1) cycle
    valtmp=surtriang(i)%value
    chgsep=chgsep+abs(valtmp-avgall)*surtriang(i)%area
    if (valtmp>=0D0) then
        varipos=varipos+surtriang(i)%area*(valtmp-avgpos)**2
    else
        varineg=varineg+surtriang(i)%area*(valtmp-avgneg)**2
    end if
end do
if (surfareapos/=0D0) varipos=varipos/surfareapos
if (surfareaneg/=0D0) varineg=varineg/surfareaneg
variall=varipos+varineg
chgsep=chgsep/surfareaall
if (imapfunc==1.or.imapfunc==3) then
    balencechg=varipos*varineg/(varipos+varineg)**2
    write(*,"(' Overall variance (sigma^2_tot):',f12.8,' a.u.^2 (',f12.5,' (kcal/mol)^2)')") variall,variall*au2kcal**2
    write(*,"(' Positive variance:     ',f13.8,' a.u.^2 (',f13.5,' (kcal/mol)^2)')") varipos,varipos*au2kcal**2
    write(*,"(' Negative variance:     ',f13.8,' a.u.^2 (',f13.5,' (kcal/mol)^2)')") varineg,varineg*au2kcal**2
    write(*,"(' Balance of charges (miu):',f13.8)") balencechg
    write(*,"(' Product of sigma^2_tot and miu: ',f12.8,' a.u.^2 (',f11.5,' (kcal/mol)^2)')") balencechg*variall,balencechg*variall*au2kcal**2
    write(*,"(' Internal charge separation (Pi):',f13.8,' a.u. (',f13.5,' kcal/mol)')") chgsep,chgsep*au2kcal
else if (imapfunc==2) then
    write(*,"(' Variance:  ',f13.8,' a.u.^2  (',f13.5,' eV^2,',E13.5,' kcal/mol^2)')") variall,variall*au2eV**2,variall*au2kcal**2
else if (imapfunc==4) then
    write(*,"(' Overall variance: ',f10.6,' a.u.^2  (',f10.5,' eV^2,',E12.5,' kcal/mol^2)')") variall,variall*au2eV**2,variall*au2kcal**2
    write(*,"(' Positive variance:',f10.6,' a.u.^2  (',f10.5,' eV^2,',E12.5,' kcal/mol^2)')") varipos,varipos*au2eV**2,varipos*au2kcal**2
    write(*,"(' Negative variance:',f10.6,' a.u.^2  (',f10.5,' eV^2,',E12.5,' kcal/mol^2)')") varineg,varineg*au2eV**2,varineg*au2kcal**2
else
    write(*,"(' Overall variance:  ',f20.10)") variall
    write(*,"(' Positive variance: ',f20.10)") varipos
    write(*,"(' Negative variance: ',f20.10)") varineg
end if
write(*,*)
write(*,*) "Surface analysis finished!"

call walltime(iclktime2)
write(*,"(' The total wall clock time passed during the task:',i6,'s')") iclktime2-iclktime1


!============= Post-processing step
textheigh=30.0D0
ratioatmsphere=1.0D0
bondradius=0.2D0
ishowatmlab=1
ishowaxis=1
ishowlocminlab=0
ishowlocmaxlab=0
ishowlocminpos=1
ishowlocmaxpos=1
do while(.true.)
    write(*,*)
    write(*,*) "                   ========== Post-process interface =========="
    write(*,*) "-3 Visualize the surface"
    write(*,*) "-2 Export the grid data to surf.cub in current folder"
    write(*,*) "-1 Return to upper level menu"
    write(*,*) "0 View molecular structure, surface minima and maxima"
    write(*,*) "1 Export surface extrema as plain text file"
    write(*,*) "2 Export surface extrema as pdb file"
    write(*,*) "3 Discard surface minima in certain value range"
    write(*,*) "4 Discard surface maxima in certain value range"
    write(*,*) "5 Export molecule as pdb format file"
    write(*,*) "6 Export all surface vertices to vtx.pdb in current folder" !if 66, also output connectivity
    write(*,*) "7 Export all surface vertices to vtx.txt in current folder"
!     write(*,*) "8 Export center of surface facets as pdb file" !Can also output to xyz file
    write(*,*) "9 Output surface area in specific value range of mapped function"
    write(*,*) "10 Output the closest and farthest distance between the surface and a point"
    write(*,*) "11 Output surface properties of each atom"
    write(*,*) "12 Output surface properties of specific fragment" !if -12, can define additional geometry rule
    if (imapfunc==22) write(*,*) "20 Fingerprint plot analysis"
    read(*,*) isel
    
    if (isel==-3) then
        sur_value=surfisoval
        
    else if (isel==-2) then
        open(10,file="surf.cub",status="replace")
        call outcube(cubmat,nx,ny,nz,orgx,orgy,orgz,dx,dy,dz,10)
        close(10)
        write(*,*) "Done, the grid data has been exported to surf.cub in current folder"
        
    else if (isel==-1) then
        exit
        
    else if (isel==0.or.isel==1) then
        if (isel==0) then
            ides=6
        else if (isel==1) then
            ides=10
            open(ides,file="surfanalysis.txt",status="replace")
        end if
        write(ides,"(a,i6)") "Number of surface minima:",count(surlocminidx(1:nsurlocmin)/=0)
        if (imapfunc==1.or.imapfunc==2.or.imapfunc==3.or.imapfunc==4) then
            write(ides,*) "  #       a.u.         eV      kcal/mol           X/Y/Z coordinate(Angstrom)"
        else
            write(ides,*) "  #             Value           X/Y/Z coordinate(Angstrom)"
        end if
        do i=1,nsurlocmin
            idx=surlocminidx(i)
            if (idx==0) cycle
            char1tmp=" "
            if (idx==indsurfmin) char1tmp="*"
            if (imapfunc==1.or.imapfunc==2.or.imapfunc==3.or.imapfunc==4) then
                write(ides,"(a,i5,f12.8,f12.6,f12.6,4x,3f11.6)") char1tmp,i,survtx(idx)%value,&
                survtx(idx)%value*au2eV,survtx(idx)%value*au2kcal,survtx(idx)%x*b2a,survtx(idx)%y*b2a,survtx(idx)%z*b2a
            else
                write(ides,"(a,i5,f18.6,4x,3f11.6)") char1tmp,i,survtx(idx)%value,survtx(idx)%x*b2a,survtx(idx)%y*b2a,survtx(idx)%z*b2a
            end if
        end do
        write(ides,*)
        write(ides,"(a,i6)") "Number of surface maxima:",count(surlocmaxidx(1:nsurlocmax)/=0)
        if (imapfunc==1.or.imapfunc==2.or.imapfunc==3.or.imapfunc==4) then
            write(ides,*) "  #       a.u.         eV      kcal/mol           X/Y/Z coordinate(Angstrom)"
        else
            write(ides,*) "  #             Value           X/Y/Z coordinate(Angstrom)"
        end if
        do i=1,nsurlocmax
            idx=surlocmaxidx(i)
            if (idx==0) cycle
            char1tmp=" "
            if (idx==indsurfmax) char1tmp="*"
            if (imapfunc==1.or.imapfunc==2.or.imapfunc==3.or.imapfunc==4) then
                write(ides,"(a,i5,f12.8,f12.6,f12.6,4x,3f11.6)") char1tmp,i,survtx(idx)%value,&
                survtx(idx)%value*au2eV,survtx(idx)%value*au2kcal,survtx(idx)%x*b2a,survtx(idx)%y*b2a,survtx(idx)%z*b2a
            else
                write(ides,"(a,i5,f18.6,4x,3f11.6)") char1tmp,i,survtx(idx)%value,survtx(idx)%x*b2a,survtx(idx)%y*b2a,survtx(idx)%z*b2a
            end if
        end do
        if (isel==0) then
        else if (isel==1) then
            close(10)
            write(*,"(a)") " Results have been outputted to surfanalysis.txt in current folder"
        end if
        
    else if (isel==2) then
        !Output positions of local maximum (as carbon) and local minimum (as oxygen) to pdb file
        open(10,file="surfanalysis.pdb",status="replace")
        do i=1,nsurlocmax
            idx=surlocmaxidx(i)
            if (idx==0) cycle
            if (imapfunc==1.or.imapfunc==3) then
                tmpval=survtx(idx)%value*au2kcal
            else if (imapfunc==2.or.imapfunc==4) then
                tmpval=survtx(idx)%value*au2eV
            else
                tmpval=survtx(idx)%value
            end if
            write(10,"(a6,i5,1x,a4,1x,a3, 1x,a1,i4,4x,3f8.3,2f6.2,10x,a2)") &
            "HETATM",i,' '//"C "//' ',"MOL",'A',1,survtx(idx)%x*b2a,survtx(idx)%y*b2a,survtx(idx)%z*b2a,1.0,tmpval,"C "
        end do
        do i=1,nsurlocmin
            idx=surlocminidx(i)
            if (idx==0) cycle
            if (imapfunc==1.or.imapfunc==3) then
                tmpval=survtx(idx)%value*au2kcal
            else if (imapfunc==2.or.imapfunc==4) then
                tmpval=survtx(idx)%value*au2eV
            else
                tmpval=survtx(idx)%value
            end if
            write(10,"(a6,i5,1x,a4,1x,a3, 1x,a1,i4,4x,3f8.3,2f6.2,10x,a2)") &
            "HETATM",i,' '//"O "//' ',"MOL",'A',1,survtx(idx)%x*b2a,survtx(idx)%y*b2a,survtx(idx)%z*b2a,1.0,tmpval,"O "
        end do
        write(10,"('END')")
        close(10)
        write(*,"(a)") " Results have been outputted to surfanalysis.pdb in current folder"
        write(*,"(a)",advance='no') " Note: Carbons and oxygens correspond to local maximum and minimum points respectively, "
        if (imapfunc==1.or.imapfunc==3) then
            write(*,"(a)") "function values (kcal/mol) are recorded in B-factor field"
        else if (imapfunc==2.or.imapfunc==4) then
            write(*,"(a)") "function values (eV) are recorded in B-factor field"
        else
            write(*,"(a)") "function values (in original unit) are recorded in B-factor field"
        end if
        
    else if (isel==3) then
        write(*,*) "Input value range (in a.u.),  e.g. 0.2,999"
        read(*,*) vallowlim,valhighlim
        do i=1,nsurlocmin
            idx=surlocminidx(i)
            if (idx==0) cycle !This slot has been discarded before
            if (survtx(idx)%value>vallowlim.and.survtx(idx)%value<valhighlim) surlocminidx(i)=0
        end do
        write(*,"(' Surface minima with value from',f10.5,' to',f10.5,' have been discarded')") vallowlim,valhighlim
        
    else if (isel==4) then
        write(*,*) "Input value range (in a.u.),  e.g. 0.2,999"
        read(*,*) vallowlim,valhighlim
        do i=1,nsurlocmax
            idx=surlocmaxidx(i)
            if (idx==0) cycle !This slot has been discarded before
            if (survtx(idx)%value>vallowlim.and.survtx(idx)%value<valhighlim) surlocmaxidx(i)=0
        end do
        write(*,"(' Surface maxima with value from',f10.5,' to',f10.5,' have been discarded')") vallowlim,valhighlim
        
    else if (isel==5) then
        write(*,*) "Input the filename you want to save to, e.g. c:\K-ON\Mio.pdb"
        read(*,"(a)") pdbfilename
        call outpdb(pdbfilename,10)
        
    else if (isel==6.or.isel==66) then !if 66, also output connectivity
        if (imapfunc==11.or.imapfunc==12) then
            !The density value in intermolecular region is always very small, and beta column in .pdb file is very narrow, we need scale it by a factor
            write(*,*) "Multiply the density by what value? e.g. 100"
            read(*,*) tmpfac
        end if
        open(10,file="vtx.pdb",status="replace")
        write(10,"('REMARK   Generated by Multiwfn, totally',i10,' surface vertices')") nsurvtx
        do i=1,nsurvtx
            if (elimvtx(i)==0) then
                if (imapfunc==1.or.imapfunc==3) then
                    tmpfuncval=survtx(i)%value*au2kcal
                else if (imapfunc==2.or.imapfunc==4) then
                    tmpfuncval=survtx(i)%value*au2eV
                else if (imapfunc==11.or.imapfunc==12) then
                    tmpfuncval=survtx(i)%value*tmpfac
                else
                    tmpfuncval=survtx(i)%value
                end if
                if (tmpfuncval>999.99D0) tmpfuncval=999.99D0 !Avoid excess limit then become, because B-factor field only have three integer position
                if (tmpfuncval<-99.99D0) tmpfuncval=-99.99D0
                write(10,"(a6,i5,1x,a4,1x,a3, 1x,a1,i4,4x,3f8.3,2f6.2,10x,a2)") &
                "HETATM",i,' '//"C "//' ',"MOL",'A',1,survtx(i)%x*b2a,survtx(i)%y*b2a,survtx(i)%z*b2a,1.0,tmpfuncval,"C "
            else !has been eliminated, B-factor is the minimal value at global surface minima
                if (isel==66) then
                    if (imapfunc==-1.or.imapfunc==0) tmpfuncval=survtx(indsurfmin)%value
                    if (imapfunc==1.or.imapfunc==3) tmpfuncval=survtx(indsurfmin)%value*au2kcal
                    if (imapfunc==2.or.imapfunc==4) tmpfuncval=survtx(indsurfmin)%value*au2eV
                    write(10,"(a6,i5,1x,a4,1x,a3, 1x,a1,i4,4x,3f8.3,2f6.2,10x,a2)") &
                    "HETATM",i,' '//"O "//' ',"MOL",'A',1,survtx(i)%x*b2a,survtx(i)%y*b2a,survtx(i)%z*b2a,1.0,tmpfuncval,"O "
                end if
            end if
        end do
        if (isel==66) then
            do i=1,nsurvtx !Output connectivity
                if (elimvtx(i)==1) cycle
                write(10,"('CONECT',i6)",advance='no') i 
                do j=1,vtxconnpos(i)
                    if (elimvtx(vtxconn(i,j))==0) write(10,"(i6)",advance='no') vtxconn(i,j)
                end do
                write(10,*)
            end do
        end if
        write(10,"('END')")
        close(10)
        write(*,"(a)") " Surface vertices have been outputted to vtx.pdb in current folder"
        if (imapfunc==1.or.imapfunc==3) then
            write(*,"(a)") " B-factor field records mapped function value in kcal/mol"
        else if (imapfunc==2.or.imapfunc==4) then
            write(*,"(a)") " B-factor field records mapped function value in eV"
        else
            write(*,"(a)") " B-factor field records mapped function value in original unit"
        end if
        if (isel==6) then
            write(*,"(a)") " Note: The eliminated vertices and connectivity are not outputted. If you would like to output them you should select option 66 (a hidden option)"
        else if (isel==66) then
            write(*,"(a)") " Note: Carbons and oxygens correspond to actually used and eliminated vertices respecitvely. CONECT field records connectivity between the vertices."
        end if
        
        !Output vertices with three connnections
!         open(10,file="vtx3.pdb",status="replace")
!         do i=1,nsurvtx
!             if (elimvtx(i)==1) cycle
!             numconn=0
!             do j=1,vtxconnpos(i)
!                 if (elimvtx(vtxconn(i,j))==0) numconn=numconn+1
!             end do
!             if (numconn==3) write(10,"(a6,i5,1x,a4,1x,a3, 1x,a1,i4,4x,3f8.3,2f6.2,10x,a2)") &
!                 "HETATM",i,' '//"N "//' ',"MOL",'A',1,survtx(i)%x*b2a,survtx(i)%y*b2a,survtx(i)%z*b2a,1.0,0.0,"N "
!         end do
!         close(10)

    else if (isel==7) then
        open(10,file="vtx.txt",status="replace")
        write(10,"(i10)") ncurrvtx
        do i=1,nsurvtx
            if (elimvtx(i)==0) then
                if (imapfunc==1.or.imapfunc==2.or.imapfunc==3.or.imapfunc==4) then
                    write(10,"(3f13.7,4x,3f16.10)") survtx(i)%x,survtx(i)%y,survtx(i)%z,survtx(i)%value,survtx(i)%value*au2eV,survtx(i)%value*au2kcal
                else
                    write(10,"(3f13.7,4x,f18.10)") survtx(i)%x,survtx(i)%y,survtx(i)%z,survtx(i)%value
                end if
            end if
        end do
        write(*,"(a)") " Done, all surface vertices (not including the ones have been eliminated) have been outputted to vtx.txt in current folder"
        if (imapfunc==1.or.imapfunc==2.or.imapfunc==3.or.imapfunc==4) then
            write(*,"(a)") " The first line is the number of points. Column 1,2,3,4,5,6 respectively correspond to coordinate of X/Y/Z in Bohr, mapped function in a.u., eV and kcal/mol"
        else
            write(*,"(a)") " The first line is the number of points. Column 1,2,3,4 respectively correspond to coordinate of X/Y/Z in Bohr and mapped function in original unit"
        end if
        close(10)
        
    else if (isel==8) then
        !Output facets via xyz file for debugging
!         open(10,file="fac.xyz",status="replace")
!         write(10,*) count(elimtri==0)
!         write(10,*)
!         do i=1,nsurtri
!             if (elimtri(i)==1) cycle
!             write(10,"(a,3f14.8)") "O   ",surtriang(i)%cenx*b2a,surtriang(i)%ceny*b2a,surtriang(i)%cenz*b2a
!         end do
!         close(10)
!         write(*,*) "Center of surface facets have been outputted to fac.xyz in current folder"

        open(10,file="tri.pdb",status="replace")
        write(10,"('REMARK   Generated by Multiwfn, totally',i10,' surface triangles')") nsurtri
        do itri=1,nsurtri
            if (elimtri(itri)/=0) cycle
            surtrix=sum(survtx(surtriang(itri)%idx(1:3))%x)/3D0 !Center of triangles
            surtriy=sum(survtx(surtriang(itri)%idx(1:3))%y)/3D0
            surtriz=sum(survtx(surtriang(itri)%idx(1:3))%z)/3D0
            if (imapfunc==-1.or.imapfunc==0) tmpfuncval=surtriang(itri)%value
            if (imapfunc==1.or.imapfunc==3) tmpfuncval=surtriang(itri)%value*au2kcal
            if (imapfunc==2.or.imapfunc==4) tmpfuncval=surtriang(itri)%value*au2eV
            if (tmpfuncval>999.99D0) tmpfuncval=999.99D0 !Avoid excess limit then become, because B-factor field only have three integer position
            if (tmpfuncval<-99.99D0) tmpfuncval=-99.99D0
            write(10,"(a6,i5,1x,a4,1x,a3, 1x,a1,i4,4x,3f8.3,2f6.2,10x,a2)") &
            "HETATM",itri,' '//"C "//' ',"MOL",'A',1,surtrix*b2a,surtriy*b2a,surtriz*b2a,1.0,tmpfuncval,"C "
        end do
        write(*,*) "Center of surface triangles have been outputted to tri.pdb in current folder"
        close(10)
        
    else if (isel==9) then
        write(*,"(a)") " Input atomic indices to define the fragment. e.g. 1,3-6,8,10-11 means the atoms 1,3,4,5,6,8,10,11 will constitute the fragment"
        write(*,*) "If input ""all"", then the whole system will be considered"
        read(*,"(a)") c2000tmp
        if (allocated(surtrifrag)) deallocate(surtrifrag)
        allocate(surtrifrag(nsurtri))
        if (index(c2000tmp,'all')/=0) then
            surtrifrag=1 !All vertices are taken into account
        else
            call str2arr(c2000tmp,nsurfragatm,surfrag) !surfrag contains nsurfragatm elements, they are the member of the user-defined fragment
            do iatm=1,ncenter !Generate table ifatmfrag, determine if an atom is belong to the user-defined fragment
                if ( any(surfrag(1:nsurfragatm)==iatm) ) then
                    ifatmfrag(iatm)=1
                else
                    ifatmfrag(iatm)=0
                end if
            end do
            !Determine the triangles belong to which atom, use Voronoi-like partition
            surtrifrag=0 !The eliminated vertices will be attributed to null (0)
            do itri=1,nsurtri
                if (elimtri(itri)==1) cycle
                effdistmax=-1D50
                surtrix=sum(survtx(surtriang(itri)%idx(1:3))%x)/3D0 !Center of triangles
                surtriy=sum(survtx(surtriang(itri)%idx(1:3))%y)/3D0
                surtriz=sum(survtx(surtriang(itri)%idx(1:3))%z)/3D0
                do iatm=1,ncenter
                    disttmp=dsqrt((surtrix-a(iatm)%x)**2+(surtriy-a(iatm)%y)**2+(surtriz-a(iatm)%z)**2)
                     if (imolsurparmode==1) then
                         vdwrtmp=1D0 !Use original Voronoi partition
                     else
                        vdwrtmp=vdwr(a(iatm)%index)
                    end if
                    effdisttmp=1-disttmp/vdwrtmp
                    if (effdisttmp>effdistmax) then
                        surtrifrag(itri)=iatm
                        effdistmax=effdisttmp
                    end if
                end do
                surtrifrag(itri)=ifatmfrag(surtrifrag(itri)) ! Then surtrifrag(itri)=1/0 means this itri belongs / doesn't belong to user-defined fragment
            end do
        end if
        write(*,*) "Input range of mapped function value, e.g. -45,50.5"
        read(*,*) rangelow,rangehigh
        write(*,*) "Input the number of intervals, e.g. 5"
        read(*,*) nintval
        convunit=1D0
        if (imapfunc==1.or.imapfunc==2.or.imapfunc==3.or.imapfunc==4) then
            write(*,*) "The inputted range is in which unit?"
            write(*,*) "1: a.u.     2: eV     3: kcal/mol"
            read(*,*) iunit
            if (iunit==2) convunit=au2eV
            if (iunit==3) convunit=au2kcal
        end if
        step=(rangehigh-rangelow)/nintval
        write(*,*) "Note: Area unit is in Angstrom^2"
        write(*,*) "     Begin        End       Center       Area         %"
        areatmptot=0D0
        do istep=1,nintval
            tmplow=rangelow+(istep-1)*step
            tmphigh=rangelow+istep*step
            areatmp=0D0
            do itri=1,nsurtri
                if (elimtri(itri)==1) cycle
                if (surtrifrag(itri)==0) cycle
                if (surtriang(itri)%value*convunit>=tmplow.and.surtriang(itri)%value*convunit<tmphigh) areatmp=areatmp+surtriang(itri)%area
            end do
            areatmp=areatmp*b2a**2
            areatmptot=areatmptot+areatmp
            write(*,"(5f12.4)") tmplow,tmphigh,(tmphigh+tmplow)/2D0,areatmp,areatmp/(surfareaall*b2a**2)*100
        end do
        write(*,"(' Sum:',31x,2f12.4)") areatmptot,areatmptot/(surfareaall*b2a**2)*100
        
    else if (isel==10) then
        write(*,*) "Input XYZ coordinate of the point (in Angstrom), e.g. 0.2,-4.3,1.66"
        write(*,*) "or input atomic index to use its nuclear position as the point, e.g. 5" 
        write(*,*) "or input ""g"" to set geometry center of present system as the point"
        write(*,"(a)") " If input ""f"", then the farthest distance between all surface points will be outputted"
        read(*,"(a)") c200tmp
        if (index(c200tmp,'f')/=0) then
            write(*,*) "Calculating, please wait..."
            distmax2=0
            do ivtx=1,nsurvtx
                do jvtx=ivtx+1,nsurvtx
                    dist2=(survtx(ivtx)%x-survtx(jvtx)%x)**2+(survtx(ivtx)%y-survtx(jvtx)%y)**2+(survtx(ivtx)%z-survtx(jvtx)%z)**2
                    if (dist2>distmax2) distmax2=dist2
                end do
            end do
            write(*,"(' The farthest distance is:',f12.6,' Bohr (',f12.6,' Angstrom)')") dsqrt(distmax2),dsqrt(distmax2)*b2a
        else
            if (index(c200tmp,',')/=0) then
                read(c200tmp,*) tmpx,tmpy,tmpz
                tmpx=tmpx/b2a
                tmpy=tmpy/b2a
                tmpz=tmpz/b2a
            else if (index(c200tmp,'g')/=0) then
                tmpx=sum(a(:)%x)/ncenter
                tmpy=sum(a(:)%y)/ncenter
                tmpz=sum(a(:)%z)/ncenter
            else
                read(c200tmp,*) iatm
                tmpx=a(iatm)%x
                tmpy=a(iatm)%y
                tmpz=a(iatm)%z
            end if
            do i=1,nsurvtx
                dist=dsqrt( (survtx(i)%x-tmpx)**2+(survtx(i)%y-tmpy)**2+(survtx(i)%z-tmpz)**2 )
                if (i==1.or.dist<distmin) distmin=dist
                if (i==1.or.dist>distmax) distmax=dist
            end do
            if (index(c200tmp,',')==0) write(*,"(' The XYZ coordinate of the point you chosen:',/,3f12.6,' Angstrom')") tmpx*b2a,tmpy*b2a,tmpz*b2a
            write(*,"(a,f12.6,a,f12.6,a)") " The closest distance to the point:",distmin," Bohr (",distmin*b2a," Angstrom)"
            write(*,"(a,f12.6,a,f12.6,a)") " The farthest distance to the point:",distmax," Bohr (",distmax*b2a," Angstrom)"
        end if
        write(*,*)
        
    else if (isel==11.or.isel==12.or.isel==-12) then !Separate the properties of the whole molecular surface into atomic or fragment local surface
        if (isel==11) then
            nsurfrag=ncenter !Each fragment corresponds to an atom
        else if (abs(isel)==12) then
            nsurfrag=1 !User define only one fragment
            write(*,"(a)") " Input atomic indices to define the fragment. e.g. 1,3-6,8,10-11 means the atoms 1,3,4,5,6,8,10,11 will constitute the fragment"
            read(*,"(a)") c2000tmp
            call str2arr(c2000tmp,nsurfragatm,surfrag) !surfrag contains nsurfragatm elements, they are the member of the user-defined fragment
            do iatm=1,ncenter !Generate table ifatmfrag, determine if an atom is belong to the user-defined fragment
                if ( any(surfrag(1:nsurfragatm)==iatm) ) then
                    ifatmfrag(iatm)=1
                else
                    ifatmfrag(iatm)=0
                end if
            end do
            if (isel==-12) then !Can input additional geometry rule to determine the local surface. Currently not used
!                 write(*,*) "Input additional geometry rule"
!                 read(*,*) geomrule
            end if
        end if
        !Determine the triangles belong to which atom, use Voronoi-like partition
        if (allocated(surtrifrag)) deallocate(surtrifrag)
        allocate(surtrifrag(nsurtri))
        surtrifrag=0 !The eliminated vertices will be attributed to null (0)
        do itri=1,nsurtri
            if (elimtri(itri)==1) cycle
            effdistmax=-1D50
            surtrix=sum(survtx(surtriang(itri)%idx(1:3))%x)/3D0 !Center of triangles
            surtriy=sum(survtx(surtriang(itri)%idx(1:3))%y)/3D0
            surtriz=sum(survtx(surtriang(itri)%idx(1:3))%z)/3D0
            
            ! Additional condition for determining local surface!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!             if (a(surfrag(1))%name=="Si") then
!                 tmp=2.5D0/b2a
!                 if (isel==12.and.(surtriz>tmp.or.surtriz<-tmp)) cycle  !equatorial
!                 if (isel==-12.and.(surtriz<tmp.and.surtriz>-tmp)) cycle  !end
!             end if

            do iatm=1,ncenter
                disttmp=dsqrt((surtrix-a(iatm)%x)**2+(surtriy-a(iatm)%y)**2+(surtriz-a(iatm)%z)**2)
                 if (imolsurparmode==1) then
                     vdwrtmp=1D0 !Use original Voronoi partition
                 else
                    vdwrtmp=dsqrt(vdwr(a(iatm)%index))
                end if
                effdisttmp=1-disttmp/vdwrtmp
                if (effdisttmp>effdistmax) then
                    surtrifrag(itri)=iatm
                    effdistmax=effdisttmp
                end if
            end do
            if (abs(isel)==12) surtrifrag(itri)=ifatmfrag(surtrifrag(itri)) ! Then surtrifrag(itri)=1/0 means this itri belongs / doesn't belong to user-defined fragment
        end do
!         do ifrag=1,nsurfrag    !Show the number of constituent facets for each fragment
!             write(*,*) ifrag,count(surtrifrag==ifrag)
!         end do
        !Obtain atomic or fragmental values
        fragsurarea=0D0
        fragsuravg=0D0
        fragsurvar=0D0
        fragsurmin=1D50
        fragsurmax=-1D50
        do itri=1,nsurtri !Get area and average
            if (elimtri(itri)==1) cycle
            iattfrag=surtrifrag(itri)
            if (iattfrag==0) cycle
            tmpval=surtriang(itri)%value
            tmpvalmin=minval(survtx(surtriang(itri)%idx(1:3))%value) !Find minimal and maximal value from its three vertices
            tmpvalmax=maxval(survtx(surtriang(itri)%idx(1:3))%value)
            fragsurarea(iattfrag,1)=fragsurarea(iattfrag,1)+surtriang(itri)%area
            fragsuravg(iattfrag,1)=fragsuravg(iattfrag,1)+surtriang(itri)%area*tmpval
            if (tmpval>0) then
                fragsurarea(iattfrag,2)=fragsurarea(iattfrag,2)+surtriang(itri)%area
                fragsuravg(iattfrag,2)=fragsuravg(iattfrag,2)+surtriang(itri)%area*tmpval
            else if (tmpval<0) then
                fragsurarea(iattfrag,3)=fragsurarea(iattfrag,3)+surtriang(itri)%area
                fragsuravg(iattfrag,3)=fragsuravg(iattfrag,3)+surtriang(itri)%area*tmpval
            end if
            if (tmpvalmax>fragsurmax(iattfrag)) fragsurmax(iattfrag)=tmpvalmax
            if (tmpvalmin<fragsurmin(iattfrag)) fragsurmin(iattfrag)=tmpvalmin
        end do
        fragsuravg=fragsuravg/fragsurarea
        fragsurchgsep=0D0
        do itri=1,nsurtri !Get variance
            if (elimtri(itri)==1) cycle
            iattfrag=surtrifrag(itri)
            if (iattfrag==0) cycle
            tmpval=surtriang(itri)%value
            if (tmpval>0) then
                fragsurvar(iattfrag,2)=fragsurvar(iattfrag,2)+surtriang(itri)%area*(tmpval-fragsuravg(iattfrag,2))**2
            else if (tmpval<0) then
                fragsurvar(iattfrag,3)=fragsurvar(iattfrag,3)+surtriang(itri)%area*(tmpval-fragsuravg(iattfrag,3))**2
            end if
            fragsurchgsep(iattfrag)=fragsurchgsep(iattfrag)+surtriang(itri)%area*abs(tmpval-fragsuravg(iattfrag,1))
        end do
        fragsurvar=fragsurvar/fragsurarea
        fragsurvar(:,1)=fragsurvar(:,2)+fragsurvar(:,3)
        fragsurchgsep=fragsurchgsep/fragsurarea(:,1)
        !Print atomic values
        write(*,*)
        if (isel==11) write(*,*) "Note: The atoms having zero surface area are not shown below"
        if (abs(isel)==12) write(*,*) "Properties on the surface of this fragment:"
        if (imapfunc==1.or.imapfunc==3) then !ESP
            if (isel==11) then
                write(*,*) "Note: Minimal and maximal value below are in kcal/mol"
                write(*,*) " Atom#    All/Positive/Negative area (Ang^2)  Minimal value   Maximal value"
                do ifrag=1,nsurfrag
                    if (fragsurarea(ifrag,1)==0D0) cycle
                    write(*,"(i7,1x,3f12.5,2f16.8)") ifrag,fragsurarea(ifrag,:)*b2a*b2a,fragsurmin(ifrag)*au2kcal,fragsurmax(ifrag)*au2kcal
                end do
                write(*,*)
                write(*,*) "Note: Average and variance below are in kcal/mol and (kcal/mol)^2 respectively"
                write(*,*) " Atom#    All/Positive/Negative average       All/Positive/Negative variance"
                do ifrag=1,nsurfrag
                    if (fragsurarea(ifrag,1)==0D0) cycle
                    write(*,"(i7,1x,3f11.5,3x,3f11.5)") ifrag,fragsuravg(ifrag,:)*au2kcal,fragsurvar(ifrag,:)*au2kcal**2
                end do
                write(*,*)
                write(*,*) "Note: Internal charge separation (Pi) is in kcal/mol, miu = Balance of charges"
                write(*,*) " Atom#           Pi              miu         miu*sigma^2"
                do ifrag=1,nsurfrag
                    if (fragsurarea(ifrag,1)==0D0) cycle
                    balencechg=fragsurvar(ifrag,2)*fragsurvar(ifrag,3)/fragsurvar(ifrag,1)**2
                    write(*,"(i7,1x,3f16.6)") ifrag,fragsurchgsep(ifrag)*au2kcal,balencechg,balencechg*fragsurvar(ifrag,1)*au2kcal**2
                end do
            else if (abs(isel)==12) then
                balencechg=fragsurvar(1,2)*fragsurvar(1,3)/fragsurvar(1,1)**2
                write(*,"('Minimal value:',f13.6,' kcal/mol   Maximal value:',f13.6,' kcal/mol')") fragsurmin(1)*au2kcal,fragsurmax(1)*au2kcal
                write(*,"('Overall surface area:      ',f12.5,' Bohr^2  (',f10.5,' Angstrom^2)')") fragsurarea(1,1),fragsurarea(1,1)*b2a*b2a
                write(*,"('Positive surface area:     ',f12.5,' Bohr^2  (',f10.5,' Angstrom^2)')") fragsurarea(1,2),fragsurarea(1,2)*b2a*b2a
                write(*,"('Negative surface area:     ',f12.5,' Bohr^2  (',f10.5,' Angstrom^2)')") fragsurarea(1,3),fragsurarea(1,3)*b2a*b2a
                write(*,"('Overall average value: ',f13.8,' a.u. (',f14.8,' kcal/mol)')") fragsuravg(1,1),fragsuravg(1,1)*au2kcal
                write(*,"('Positive average value:',f13.8,' a.u. (',f14.8,' kcal/mol)')") fragsuravg(1,2),fragsuravg(1,2)*au2kcal
                write(*,"('Negative average value:',f13.8,' a.u. (',f14.8,' kcal/mol)')") fragsuravg(1,3),fragsuravg(1,3)*au2kcal
                write(*,"('Overall variance (sigma^2_tot):',f12.8,' a.u.^2 (',f13.7,' (kcal/mol)^2)')") fragsurvar(1,1),fragsurvar(1,1)*au2kcal**2
                write(*,"('Positive variance:     ',f13.8,' a.u.^2 (',f14.8,' (kcal/mol)^2)')") fragsurvar(1,2),fragsurvar(1,2)*au2kcal**2
                write(*,"('Negative variance:     ',f13.8,' a.u.^2 (',f14.8,' (kcal/mol)^2)')") fragsurvar(1,3),fragsurvar(1,3)*au2kcal**2
                write(*,"('Balance of charges (miu):',f13.8)") balencechg
                write(*,"('Product of sigma^2_tot and miu: ',f12.8,' a.u.^2 (',f12.7,' (kcal/mol)^2)')") balencechg*fragsurvar(1,1),balencechg*fragsurvar(1,1)*au2kcal**2
                write(*,"('Internal charge separation (Pi):',f13.8,' a.u. (',f13.8,' kcal/mol)')") fragsurchgsep(1),fragsurchgsep(1)*au2kcal
            end if
        else if (imapfunc==2) then !ALIE
            if (isel==11) then
                write(*,*) "Minimal, maximal and average value are in eV, variance is in eV^2"
                write(*,*) " Atom#      Area(Ang^2)  Min value   Max value       Average        Variance"
                do ifrag=1,nsurfrag
                    if (fragsurarea(ifrag,2)==0D0) cycle
                    write(*,"(i7,f16.5,2f12.6,2f15.6)") ifrag,fragsurarea(ifrag,2)*b2a*b2a,fragsurmin(ifrag)*au2eV,fragsurmax(ifrag)*au2eV,fragsuravg(ifrag,2)*au2eV,fragsurvar(ifrag,2)*au2eV**2
                end do
            else if (abs(isel)==12) then
                write(*,"('Minimal value:',f13.6,' eV   Maximal value:',f13.6,' eV')") fragsurmin(1)*au2eV,fragsurmax(1)*au2eV
                write(*,"('Surface area:      ',f12.5,' Bohr^2  (',f10.5,' Angstrom^2)')") fragsurarea(1,2),fragsurarea(1,2)*b2a*b2a
                write(*,"('Average value: ',f13.8,' a.u. (',f13.8,' eV,',f14.8,' kcal/mol)')") fragsuravg(1,2),fragsuravg(1,2)*au2eV,fragsuravg(1,2)*au2kcal
                write(*,"('Variance:  ',f13.8,' a.u.^2  (',f13.8,' eV^2,',E14.6,' kcal/mol^2)')") fragsurvar(1,2),fragsurvar(1,2)*au2eV**2,fragsurvar(1,2)*au2kcal**2
            end if
        else if (imapfunc==4) then
            if (isel==11) then
                write(*,*) "Note: Below minimal and maximal values are in eV"
                write(*,*) " Atom#    All/Positive/Negative area (Ang^2)  Minimal value   Maximal value"
                do ifrag=1,nsurfrag
                    if (fragsurarea(ifrag,1)==0D0) cycle
                    write(*,"(i7,1x,3f12.5,2f16.8)") ifrag,fragsurarea(ifrag,:)*b2a*b2a,fragsurmin(ifrag)*au2eV,fragsurmax(ifrag)*au2eV
                end do
                write(*,*)
                write(*,*) "Note: Average and variance below are in eV and eV^2 respectively"
                write(*,*) " Atom#    All/Positive/Negative average       All/Positive/Negative variance"
                do ifrag=1,nsurfrag
                    if (fragsurarea(ifrag,1)==0D0) cycle
                    write(*,"(i7,1x,3f11.5,3x,3f11.5)") ifrag,fragsuravg(ifrag,:)*au2eV,fragsurvar(ifrag,:)*au2eV**2
                end do
                write(*,*)
            else if (abs(isel)==12) then
                write(*,"('Minimal value:',f14.8,' eV,   Maximal value:',f14.8,' eV')") fragsurmin(1)*au2eV,fragsurmax(1)*au2eV
                write(*,"('Overall surface area: ',f12.5,' Bohr^2  (',f10.5,' Angstrom^2)')") fragsurarea(1,1),fragsurarea(1,1)*b2a*b2a
                write(*,"('Positive surface area:',f12.5,' Bohr^2  (',f10.5,' Angstrom^2)')") fragsurarea(1,2),fragsurarea(1,2)*b2a*b2a
                write(*,"('Negative surface area:',f12.5,' Bohr^2  (',f10.5,' Angstrom^2)')") fragsurarea(1,3),fragsurarea(1,3)*b2a*b2a
                write(*,"('Overall average value: ',f13.8,' a.u. (',f14.8,' eV)')") fragsuravg(1,1),fragsuravg(1,1)*au2eV
                write(*,"('Positive average value:',f13.8,' a.u. (',f14.8,' eV)')") fragsuravg(1,2),fragsuravg(1,2)*au2eV
                write(*,"('Negative average value:',f13.8,' a.u. (',f14.8,' eV)')") fragsuravg(1,3),fragsuravg(1,3)*au2eV
                write(*,"('Overall variance: ',f13.8,' a.u.^2 (',f14.8,' eV^2)')") fragsurvar(1,1),fragsurvar(1,1)*au2eV**2
                write(*,"('Positive variance:',f13.8,' a.u.^2 (',f14.8,' eV^2)')") fragsurvar(1,2),fragsurvar(1,2)*au2eV**2
                write(*,"('Negative variance:',f13.8,' a.u.^2 (',f14.8,' eV^2)')") fragsurvar(1,3),fragsurvar(1,3)*au2eV**2
            end if
        else !Other or unknown real space function provided by user
            if (isel==11) then
                write(*,*) " Atom#    All/Positive/Negative area (Ang^2)  Minimal value   Maximal value"
                do ifrag=1,nsurfrag
                    if (fragsurarea(ifrag,1)==0D0) cycle
                    write(*,"(i7,1x,3f12.5,2f16.9)") ifrag,fragsurarea(ifrag,:)*b2a*b2a,fragsurmin(ifrag),fragsurmax(ifrag)
                end do
                write(*,*)
                write(*,*) " Atom#   All/Positive/Negative average       All/Positive/Negative variance"
                do ifrag=1,nsurfrag
                    if (fragsurarea(ifrag,1)==0D0) cycle
                    write(*,"(i7,3(1PE13.5),3(1PE11.4))") ifrag,fragsuravg(ifrag,:),fragsurvar(ifrag,:)
                end do
            else if (abs(isel)==12) then
                write(*,"('Minimal value:',f16.8,'    Maximal value:',f16.8)") fragsurmin(1),fragsurmax(1)
                write(*,"('Overall surface area:      ',f12.5,' Bohr^2  (',f10.5,' Angstrom^2)')") fragsurarea(1,1),fragsurarea(1,1)*b2a*b2a
                write(*,"('Positive surface area:     ',f12.5,' Bohr^2  (',f10.5,' Angstrom^2)')") fragsurarea(1,2),fragsurarea(1,2)*b2a*b2a
                write(*,"('Negative surface area:     ',f12.5,' Bohr^2  (',f10.5,' Angstrom^2)')") fragsurarea(1,3),fragsurarea(1,3)*b2a*b2a
                write(*,"('Overall average value: ',f16.8)") fragsuravg(1,1)
                write(*,"('Positive average value:',f16.8)") fragsuravg(1,2)
                write(*,"('Negative average value:',f16.8)") fragsuravg(1,3)
                write(*,"('Overall variance:  ',f16.8)") fragsurvar(1,1)
                write(*,"('Positive variance: ',f16.8)") fragsurvar(1,2)
                write(*,"('Negative variance: ',f16.8)") fragsurvar(1,3)
            end if
        end if
        ! write(*,"(' Sum of area of total/pos./neg.:',3f12.6,' Bohr^2')") sum(fragsurarea(:,1)),sum(fragsurarea(:,2)),sum(fragsurarea(:,3))
        write(*,"(/,a)") " If output the surface facets to locsurf.pdb in current folder? By which you can visualize local surface via third-part visualization program such as VMD (y/n)"
        read(*,*) selectyn
        if (selectyn=='y'.or.selectyn=='Y') then
            open(10,file="locsurf.pdb",status="replace")
            write(10,"('REMARK   Generated by Multiwfn, totally',i10,' surface triangles')") nsurtri
            do itri=1,nsurtri
                if (elimtri(itri)==1) cycle
                surtrix=sum(survtx(surtriang(itri)%idx(1:3))%x)/3D0 !Center of triangles
                surtriy=sum(survtx(surtriang(itri)%idx(1:3))%y)/3D0
                surtriz=sum(survtx(surtriang(itri)%idx(1:3))%z)/3D0
                write(10,"(a6,i5,1x,a4,1x,a3, 1x,a1,i4,4x,3f8.3,2f6.2,10x,a2)") &
                "HETATM",itri,' '//"C "//' ',"MOL",'A',1,surtrix*b2a,surtriy*b2a,surtriz*b2a,1.0,dfloat(surtrifrag(itri)),"C "
            end do
            write(*,"(a)") " Coordinate of the geometry center of the surface facets have been outputted to locsurf.pdb in current folder"
            if (isel==11) write(*,"(a)") " B-factor column records the attribution of the surface facets"
            if (abs(isel)==12) write(*,"(a)") " In this file the atom with B-factor = 1/0 means corresponding surface facet belongs / does not belong to the fragment you defined"
            close(10)
        end if
        
    else if (isel==20) then !Fingerprint plot
        call fingerprt(HirBecatm,nHirBecatm)
    end if
end do

deallocate(elimvtx,elimtri,mergerelat)
deallocate(survtx,surtriang,vtxconn,vtxconnpos)
end do surfanaloop
end subroutine


!!------ Fingerprint plot analysis
subroutine fingerprt(HirBecatm,nHirBecatm)
use surfvertex
use util
use function
implicit real*8 (a-h,o-z)
integer nHirBecatm,HirBecatm(nHirBecatm),tmparr(ncenter),ifcontactvtx(nsurvtx)
integer,parameter :: nval=200
real*8 :: mat(nval,nval),vtxdnorm(nsurvtx),rlow=0.6D0,rhigh=2.6D0 !Angstrom
integer,allocatable :: inarr(:),outarr(:),notHirBecatm(:)
real*8,allocatable :: surval1(:),surval2(:)
character c2000tmp*2000,c2tmp*2
!Set default inside and outside fragment
ninarr=nHirBecatm
nnotHirBecatm=ncenter-nHirBecatm !The number of atoms do not belong to Hirshfeld/Becke fragment
noutarr=nnotHirBecatm
allocate(inarr(ninarr),outarr(noutarr),notHirBecatm(nnotHirBecatm))
inarr=HirBecatm
itmp=0
do i=1,ncenter
    if (all(HirBecatm/=i)) then
        itmp=itmp+1
        notHirBecatm(itmp)=i
    end if
end do
outarr=notHirBecatm

do while(.true.)
    write(*,*) "-1 Return"
    write(*,*) "0 Start analysis"
    write(*,"(a,i8)") " 1 Set the inside atoms, current number is",ninarr
    write(*,"(a,i8)") " 2 Set the outside atoms, current number is",noutarr
    read(*,*) isel

    if (isel==-1) then
        return
    else if (isel==1) then
        tmparr=0 !If tmparr(i)=1, i will be in this fragment
        write(*,*) "The index of current inside atoms"
        write(*,"(10i7)") inarr
        write(*,*) "Now input two conditions, the inside atoms will be their intersection"
        write(*,*)
        write(*,"(a)") " First, input index range, e.g. 1,3-6,8,10-11 means the atoms 1,3,4,5,6,8,10,11 are selected"
        write(*,"(a)") " Note: If press ENTER directly, all atoms in the Hirshfeld/Becke fragment will be taken into consideration"
        read(*,"(a)") c2000tmp
        if (c2000tmp==" ") then
            nelement=nHirBecatm
            tmparr(1:nelement)=HirBecatm
        else
            call str2arr(c2000tmp,nelement)
            call str2arr(c2000tmp,nelement,tmparr(1:nelement))
        end if
        write(*,*) "Next, input element, e.g. Cl"
        write(*,*) "Note: If press ENTER directly, the element filter will be ignored"
        read(*,"(a)") c2tmp
        deallocate(inarr)
        if (c2tmp==" ") then
            allocate(inarr(nelement))
            inarr=tmparr(1:nelement)
            ninarr=nelement
        else
            ninarr=count(a(tmparr(1:nelement))%name==c2tmp)
            allocate(inarr(ninarr))
            itmp=0
            do i=1,nelement
                if (a(tmparr(i))%name==c2tmp) then
                    itmp=itmp+1
                    inarr(itmp)=tmparr(i)
                end if
            end do
        end if
        write(*,*) "The inside atoms you chosen are"
        write(*,"(10i7)") inarr
        write(*,*)
        do i=1,ninarr
            if (all(HirBecatm/=inarr(i))) write(*,"(' Warning: atom',i7,' does not belong to the Hirshfeld/Becke fragment you previously set')") i
        end do
    else if (isel==2) then
        tmparr=0 !If tmparr(i)=1, i will be in this fragment
        write(*,*) "The index of current outside atoms"
        write(*,"(10i7)") outarr
        write(*,*) "Now input two conditions, the outside atoms will be their intersection"
        write(*,*)
        write(*,"(a)") " First, input index range, e.g. 1,3-6,8,10-11 means the atoms 1,3,4,5,6,8,10,11 are selected"
        write(*,"(a)") " Note: If press ENTER directly, all atoms do not belong to the Hirshfeld/Becke fragment will be taken into consideration"
        read(*,"(a)") c2000tmp
        if (c2000tmp==" ") then
            nelement=nnotHirBecatm
            tmparr(1:nelement)=notHirBecatm
        else
            call str2arr(c2000tmp,nelement)
            call str2arr(c2000tmp,nelement,tmparr(1:nelement))
        end if
        write(*,*) "Next, select element, e.g. Cl"
        write(*,*) "Note: If press ENTER directly, the element filter will be ignored"
        read(*,"(a)") c2tmp
        deallocate(outarr)
        if (c2tmp==" ") then
            allocate(outarr(nelement))
            outarr=tmparr(1:nelement)
            noutarr=nelement
        else
            noutarr=count(a(tmparr(1:nelement))%name==c2tmp)
            allocate(outarr(noutarr))
            itmp=0
            do i=1,nelement
                if (a(tmparr(i))%name==c2tmp) then
                    itmp=itmp+1
                    outarr(itmp)=tmparr(i)
                end if
            end do
        end if
        write(*,*) "The outside atoms you chosen are"
        write(*,"(10i7)") outarr
        write(*,*)
    else if (isel==0) then
        write(*,*) "Calculating contact surface points..."
        ifcontactvtx=0
        do icyc=1,2 !The first time count how many contact point are there, the second time record information
            if (icyc==2) allocate(surval1(ncontactvtx),surval2(ncontactvtx))
            ncontactvtx=0
            ncurrvtx=0
            do ivtx=1,nsurvtx
                if (elimvtx(ivtx)==1) cycle
                ncurrvtx=ncurrvtx+1
                dist2minin=1D100
                dist2minout=1D100
                iminin=0
                iminout=0
                do iatm=1,ncenter
                    dist2=(a(iatm)%x-survtx(ivtx)%x)**2+(a(iatm)%y-survtx(ivtx)%y)**2+(a(iatm)%z-survtx(ivtx)%z)**2
                    if (any(HirBecatm==iatm)) then !Find the closest atom in Hirshfeld/Becke fragment to the surface point
                        if (dist2<dist2minin) then
                            dist2minin=dist2
                            iminin=iatm
                        end if
                    else !Find the closest atom that does not belong to Hirshfeld/Becke fragment to the surface point
                        if (dist2<dist2minout) then
                            dist2minout=dist2
                            iminout=iatm
                        end if
                    end if
                end do
                if (any(inarr==iminin).and.any(outarr==iminout)) then
                    ncontactvtx=ncontactvtx+1
                    if (icyc==2) then
                        surval1(ncontactvtx)=dsqrt(dist2minin) !d_i
                        surval2(ncontactvtx)=dsqrt(dist2minout) !d_e
                        ifcontactvtx(ivtx)=1 !This is a vertex in contact surface
                        vtxdnorm(ivtx)=surfana_norm(survtx(ivtx)%x,survtx(ivtx)%y,survtx(ivtx)%z,nHirBecatm,HirBecatm)
                    end if
                end if
            end do
        end do
        write(*,"(' The number of contact surface points is',i10)") ncontactvtx
        write(*,"(' The number of total Hirshfeld/Becke surface points is',i10)") ncurrvtx
        write(*,"(' They occupy',f10.3,'% of total Hirshfeld/Becke surface points')") dfloat(ncontactvtx)/ncurrvtx*100D0
        
!         rhigh=max(maxval(surval1),maxval(surval2)) !Automatically determines the axis range
!         rlow=min(minval(surval1),minval(surval2))
!         range=rhigh-rlow
!         shift=range/10D0
!         rlow=(rlow-shift)*b2a !All length variables in this part are in Angstrom
!         rhigh=(rhigh+shift)*b2a
        write(*,*) "Calculating the distribution of surface vertices..."
        write(*,*)
        call xypt2densmat(surval1*b2a,surval2*b2a,ncontactvtx,mat,nval,nval,rlow,rhigh,rlow,rhigh)
        spc=(rhigh-rlow)/(nval-1)
        clrlow=1D-5
        clrhigh=maxval(mat)/2D0
        do while(.true.)
            write(*,*) "-1 Return"
            write(*,*) "0 Show fingerprint plot on screen"
            write(*,*) "1 Save fingerprint plot to current folder"
            write(*,*) "2 Export original data of fingerprint plot to finger.txt in current folder"
            write(*,"(a,f10.3,a,f10.3)") " 3 Set color scale of fingerprint plot, current: from",clrlow," to",clrhigh
            write(*,"(a,f10.3,a,f10.3)") " 4 Set the range of axes, current: from",rlow," to",rhigh
            write(*,*) "5 Export the points in fingerprint plot to finger.pdb in current folder"
            write(*,*) "10 Draw scatter map of surface points between d_i and d_e"
            write(*,*) "11 Export d_i and d_e of surface points to di_de.txt in current folder"
            read(*,*) isel3
            if (isel3==-1) then
                return
            else if (isel3==0) then
                write(*,*) "Note: X and Y axes (in Angstrom) correspond to d_i and d_e, respectively"
            else if (isel3==1) then
                isavepic=1
                isavepic=0
                write(*,*) "The graph has been saved to current folder with DISLIN prefix"
                write(*,*) "Note: X and Y axes (in Angstrom) correspond to d_i and d_e, respectively"
            else if (isel3==2) then
                open(10,file="finger.txt",status="replace")
                do i=1,nval
                    xval=xmin+(i-1)*spc
                    do j=1,nval
                        yval=ymin+(j-1)*spc !Already in Angstrom
                        write(10,"(2f11.5,f12.3)") xval,yval,mat(i,j)
                    end do
                end do
                close(10)
                write(*,"(a)") " Done! The distribution has been exported to finger.txt in current folder. &
                The first and second columns (in Angstrom) correspond to d_i and d_e, respectively"
            else if (isel3==3) then
                write(*,*) "Input lower and upper limit of color scale, e.g. 0.00001,6.5"
                read(*,*) clrlow,clrhigh
                if (clrlow==0D0) clrlow=1D-5 !Ensure the vacant region is black
            else if (isel3==4) then
                write(*,*) "Input lower and upper limit of the axes (in Angstrom), e.g. 0.5,2.8"
                read(*,*) rlow,rhigh
                spc=(rhigh-rlow)/(nval-1)
            else if (isel3==5) then
                open(10,file="finger.pdb",status="replace")
                write(10,"('REMARK   Generated by Multiwfn, totally',i10,' contact surface vertices')") ncontactvtx
                do ivtx=1,nsurvtx
                    if (ifcontactvtx(ivtx)==0) cycle
                    write(10,"(a6,i5,1x,a4,1x,a3, 1x,a1,i4,4x,3f8.3,2f6.2,10x,a2)") &
                    "HETATM",ivtx,' '//"C "//' ',"MOL",'A',1,survtx(ivtx)%x*b2a,survtx(ivtx)%y*b2a,survtx(ivtx)%z*b2a,1.0,vtxdnorm(ivtx),"C "
                end do
                write(10,"('END')")
                close(10)
                write(*,"(a)") " The vertices in contact surface have been outputted to finger.pdb in current folder"
            else if (isel3==10) then !Scatter map between d_i and d_e
                isavepic=1
                isavepic=0
                write(*,"(a)") " The graph has been saved to current folder with DISLIN prefix. X and Y axes correspond to d_i and d_e (in Angstrom), respectively"
            else if (isel3==11) then
                open(10,file="di_de.txt",status="replace")
                do ivtx=1,ncontactvtx
                    write(10,"(2f12.6)") surval1(ivtx)*b2a,surval2(ivtx)*b2a
                end do
                close(10)
                write(*,"(a)") " Done! The data has been exported to di_de.txt in current folder. &
                The first and second columns (in Angstrom) correspond to d_i and d_e, respectively"
            end if
            write(*,*)
        end do
    end if
end do
end subroutine



!!-------- Use Marching Tetrahedra algorithem, decompose cube to several tetrahedra
subroutine marchtetra(ix,iy,iz)
implicit real*8(a-h,o-z)
integer ix,iy,iz
!    7-------8
!   /|      /|
!  6-+-----4 |
!  | |     | |
!  | 3-----+-1
!  |/      |/
!  5-------2
!
!   Z
!   |
!   0---Y    
!  / 
! X

! Five tetrahedra, may lead too big spacing somewhere
! call genvertex(ix,iy,iz,1,2,5,4) 
! call genvertex(ix,iy,iz,4,5,6,7)
! call genvertex(ix,iy,iz,1,3,5,7)
! call genvertex(ix,iy,iz,1,4,7,8)
! call genvertex(ix,iy,iz,1,4,7,5)

! WFA original paper
! call genvertex(ix,iy,iz,3,2,1,8) 
! call genvertex(ix,iy,iz,3,2,4,8)
! call genvertex(ix,iy,iz,3,7,4,8)
! call genvertex(ix,iy,iz,3,4,5,2)
! call genvertex(ix,iy,iz,3,4,5,7)
! call genvertex(ix,iy,iz,6,4,5,7)

! Main-axis decomposition, all tetrahedra share 4-3 axis, as in http://paulbourke.net/geometry/polygonise/
call genvertex(ix,iy,iz,4,3,5,2) 
call genvertex(ix,iy,iz,4,3,5,6)
call genvertex(ix,iy,iz,4,3,7,6)
call genvertex(ix,iy,iz,4,3,7,8)
call genvertex(ix,iy,iz,4,3,1,8)
call genvertex(ix,iy,iz,4,3,1,2)
end subroutine

!!-------- Use marching cube algorithem, interpolate each edge. Unfinished routine, can only generated surface vertices but not connecitivity
! subroutine marchcube(ix,iy,iz)
! use surfvertex
! use defvar
! implicit real*8(a-h,o-z)
! integer ix,iy,iz
! itestc1=0
! itestc2=0
! itestc3=0
! itestc4=0
! itestc5=0
! itestc6=0
! itestc7=0
! itestc8=0
! !Test function value at each corner of cube, =0/1 means lower/higher than isovalue
! if (cubmat(ix,iy+1,iz)>=surfisoval) itestc1=1
! if (cubmat(ix+1,iy+1,iz)>=surfisoval) itestc2=1
! if (cubmat(ix,iy,iz)>=surfisoval) itestc3=1
! if (cubmat(ix+1,iy+1,iz+1)>=surfisoval) itestc4=1
! if (cubmat(ix+1,iy,iz)>=surfisoval) itestc5=1
! if (cubmat(ix+1,iy,iz+1)>=surfisoval) itestc6=1
! if (cubmat(ix,iy,iz+1)>=surfisoval) itestc7=1
! if (cubmat(ix,iy+1,iz+1)>=surfisoval) itestc8=1
! !interpolate each cube edge
! ! 1-2,2-5,5-3,3-1
! if (itestc1+itestc2==1) call vertexinterpolate(ix,iy+1,iz,ix+1,iy+1,iz)
! if (itestc2+itestc5==1) call vertexinterpolate(ix+1,iy+1,iz,ix+1,iy,iz)
! if (itestc5+itestc3==1) call vertexinterpolate(ix+1,iy,iz,ix,iy,iz)
! if (itestc3+itestc1==1) call vertexinterpolate(ix,iy,iz,ix,iy+1,iz)
! ! 1-8,2-4,5-6,3-7
! if (itestc1+itestc8==1) call vertexinterpolate(ix,iy+1,iz,ix,iy+1,iz+1)
! if (itestc2+itestc4==1) call vertexinterpolate(ix+1,iy+1,iz,ix+1,iy+1,iz+1)
! if (itestc5+itestc6==1) call vertexinterpolate(ix+1,iy,iz,ix+1,iy,iz+1)
! if (itestc3+itestc7==1) call vertexinterpolate(ix,iy,iz,ix,iy,iz+1)
! ! 8-4,4-6,6-7,7-8
! if (itestc8+itestc4==1) call vertexinterpolate(ix,iy+1,iz+1,ix+1,iy+1,iz+1)
! if (itestc4+itestc6==1) call vertexinterpolate(ix+1,iy+1,iz+1,ix+1,iy,iz+1)
! if (itestc6+itestc7==1) call vertexinterpolate(ix+1,iy,iz+1,ix,iy,iz+1)
! if (itestc7+itestc8==1) call vertexinterpolate(ix,iy,iz+1,ix,iy+1,iz+1)
! end subroutine


!!---------- Generate surface vertices from tetrahedra. inum is the vertex index within each tetrahedron 
subroutine genvertex(baseix,baseiy,baseiz,inum1,inum2,inum3,inum4)
use util
use defvar
use surfvertex
implicit real*8(a-h,o-z)
integer baseix,baseiy,baseiz,inum1,inum2,inum3,inum4
integer itestv(4),vix(4),viy(4),viz(4),newvtxind(4) !1~4 is the four vertices index within current tetrahedron. newvtxind is absolute surface vertex index
real*8 vtxval(4)
!Convert index within cube to absolute cubic index
call getvertind(vix(1),viy(1),viz(1),baseix,baseiy,baseiz,inum1)
call getvertind(vix(2),viy(2),viz(2),baseix,baseiy,baseiz,inum2)
call getvertind(vix(3),viy(3),viz(3),baseix,baseiy,baseiz,inum3)
call getvertind(vix(4),viy(4),viz(4),baseix,baseiy,baseiz,inum4)
!Values of vertices in current tetrahedron
vtxval(1)=cubmat(vix(1),viy(1),viz(1))
vtxval(2)=cubmat(vix(2),viy(2),viz(2))
vtxval(3)=cubmat(vix(3),viy(3),viz(3))
vtxval(4)=cubmat(vix(4),viy(4),viz(4))
!Test function value at each vertex of tetrahedron, =0/1 means lower/higher than isovalue
itestv=0
where (vtxval>surfisoval) itestv=1 !If the tetrahedron vertex is inside isosurface, then its itestv=1, else 0
!    2
!    |
!    |
!    1
!   / \
!  /   \
! 3     4
!Interpolate surface vertices, if it is new, then accumulate it to "survtx" array
itesttot=sum(itestv)
newvtxind=0
inewvtxind=0
if (itesttot==0) then ! external tetrahedron, return directly
    return
else if (itesttot==4) then !Internal tetrahedron, calculate its volume and then return
    pax=orgx+(vix(1)-1)*dx
    pay=orgy+(viy(1)-1)*dy
    paz=orgz+(viz(1)-1)*dz
    pbx=orgx+(vix(2)-1)*dx
    pby=orgy+(viy(2)-1)*dy
    pbz=orgz+(viz(2)-1)*dz
    pcx=orgx+(vix(3)-1)*dx
    pcy=orgy+(viy(3)-1)*dy
    pcz=orgz+(viz(3)-1)*dz
    pdx=orgx+(vix(4)-1)*dx
    pdy=orgy+(viy(4)-1)*dy
    pdz=orgz+(viz(4)-1)*dz
    tetravol0=tetravol0+gettetravol(pax,pay,paz,pbx,pby,pbz,pcx,pcy,pcz,pdx,pdy,pdz) !Volume type 0
else if (itesttot==1) then
    do iint=1,4 !Find the only internal tetrahedron vertex
        if (itestv(iint)==1) exit
    end do
    do iplt=1,4
        if (iplt==iint) cycle
        inewvtxind=inewvtxind+1
        call vertexinterpolate(vix(iint),viy(iint),viz(iint),vix(iplt),viy(iplt),viz(iplt),newvtxind(inewvtxind))
    end do
    call addvtxconn(newvtxind(1),newvtxind(2))
    call addvtxconn(newvtxind(1),newvtxind(3))
    call addvtxconn(newvtxind(2),newvtxind(3))
    nsurtri=nsurtri+1
    surtriang(nsurtri)%idx(1:3)=newvtxind(1:3)
    !Calculate volume of newly generated tetrahedron
    
    pax=orgx+(vix(iint)-1)*dx
    pay=orgy+(viy(iint)-1)*dy
    paz=orgz+(viz(iint)-1)*dz
    pbx=survtx(newvtxind(1))%x
    pby=survtx(newvtxind(1))%y
    pbz=survtx(newvtxind(1))%z
    pcx=survtx(newvtxind(2))%x
    pcy=survtx(newvtxind(2))%y
    pcz=survtx(newvtxind(2))%z
    pdx=survtx(newvtxind(3))%x
    pdy=survtx(newvtxind(3))%y
    pdz=survtx(newvtxind(3))%z
    tetravol1=tetravol1+gettetravol(pax,pay,paz,pbx,pby,pbz,pcx,pcy,pcz,pdx,pdy,pdz) !Volume type 1
else if (itesttot==3) then
    do iext=1,4 !Find the only external tetrahedron vertex
        if (itestv(iext)==0) exit
    end do
    do iplt=1,4
        if (iplt==iext) cycle
        inewvtxind=inewvtxind+1
        call vertexinterpolate(vix(iext),viy(iext),viz(iext),vix(iplt),viy(iplt),viz(iplt),newvtxind(inewvtxind))
    end do
    call addvtxconn(newvtxind(1),newvtxind(2))
    call addvtxconn(newvtxind(1),newvtxind(3))
    call addvtxconn(newvtxind(2),newvtxind(3))
    nsurtri=nsurtri+1
    surtriang(nsurtri)%idx(1:3)=newvtxind(1:3)
    !Calculate volume of newly generated tetrahedron
    pax=orgx+(vix(iext)-1)*dx
    pay=orgy+(viy(iext)-1)*dy
    paz=orgz+(viz(iext)-1)*dz
    pbx=survtx(newvtxind(1))%x
    pby=survtx(newvtxind(1))%y
    pbz=survtx(newvtxind(1))%z
    pcx=survtx(newvtxind(2))%x
    pcy=survtx(newvtxind(2))%y
    pcz=survtx(newvtxind(2))%z
    pdx=survtx(newvtxind(3))%x
    pdy=survtx(newvtxind(3))%y
    pdz=survtx(newvtxind(3))%z
    voidvol=gettetravol(pax,pay,paz,pbx,pby,pbz,pcx,pcy,pcz,pdx,pdy,pdz)
    !Calculate the entire volume of original tetrahedron
    
    pax=orgx+(vix(1)-1)*dx
    pay=orgy+(viy(1)-1)*dy
    paz=orgz+(viz(1)-1)*dz
    pbx=orgx+(vix(2)-1)*dx
    pby=orgy+(viy(2)-1)*dy
    pbz=orgz+(viz(2)-1)*dz
    pcx=orgx+(vix(3)-1)*dx
    pcy=orgy+(viy(3)-1)*dy
    pcz=orgz+(viz(3)-1)*dz
    pdx=orgx+(vix(4)-1)*dx
    pdy=orgy+(viy(4)-1)*dy
    pdz=orgz+(viz(4)-1)*dz
    wholetetravol=gettetravol(pax,pay,paz,pbx,pby,pbz,pcx,pcy,pcz,pdx,pdy,pdz)
    tetravol2=tetravol2+(wholetetravol-voidvol) !Volume type 2
else if (itesttot==2) then
    itv1=0
    itv2=0
    do iint=1,4
        if (itestv(iint)==1) then  !Find the two internal tetrahedron vertices
            if (itv1==0) then
                itv1=iint
            else
                itv2=iint
            end if
            do iother=1,4 !Find external vertices (namely itestv=0)
                if (itestv(iother)==1) cycle !If this is also interal one, skip it
                inewvtxind=inewvtxind+1
                call vertexinterpolate(vix(iint),viy(iint),viz(iint),vix(iother),viy(iother),viz(iother),newvtxind(inewvtxind))
            end do
        end if
    end do
    ! 1--2
    ! |  |
    ! 3--4
    ! 1-2-4 is the first triangle, 1-3-4 is the second one
    call addvtxconn(newvtxind(1),newvtxind(2))
    call addvtxconn(newvtxind(1),newvtxind(3))
    call addvtxconn(newvtxind(2),newvtxind(4))
    call addvtxconn(newvtxind(3),newvtxind(4))
    call addvtxconn(newvtxind(1),newvtxind(4))
    nsurtri=nsurtri+1
    surtriang(nsurtri)%idx(1:2)=newvtxind(1:2)
    surtriang(nsurtri)%idx(3)=newvtxind(4)
    nsurtri=nsurtri+1
    surtriang(nsurtri)%idx(1)=newvtxind(1)
    surtriang(nsurtri)%idx(2:3)=newvtxind(3:4)
    !Now split the widget with six vertices to three tetrahedra, so that volume can be calculated
    !decompose to tv1-sv1-sv3-sv4,tv1-tv2-sv1-sv4,tv2-sv1-sv2-sv4  (tv: tetrahedron vertex with value larger than isovalue;  sv: square vertex)
!
!           tv1-tv2
!     
!      sv1-------sv2
!      /         /
!    sv3-------sv4

    tv1x=orgx+(vix(itv1)-1)*dx
    tv1y=orgy+(viy(itv1)-1)*dy
    tv1z=orgz+(viz(itv1)-1)*dz
    tv2x=orgx+(vix(itv2)-1)*dx
    tv2y=orgy+(viy(itv2)-1)*dy
    tv2z=orgz+(viz(itv2)-1)*dz
    sv1x=survtx(newvtxind(1))%x
    sv1y=survtx(newvtxind(1))%y
    sv1z=survtx(newvtxind(1))%z
    sv2x=survtx(newvtxind(2))%x
    sv2y=survtx(newvtxind(2))%y
    sv2z=survtx(newvtxind(2))%z
    sv3x=survtx(newvtxind(3))%x
    sv3y=survtx(newvtxind(3))%y
    sv3z=survtx(newvtxind(3))%z
    sv4x=survtx(newvtxind(4))%x
    sv4y=survtx(newvtxind(4))%y
    sv4z=survtx(newvtxind(4))%z
    tetravol3=tetravol3+gettetravol(tv1x,tv1y,tv1z,sv1x,sv1y,sv1z,sv3x,sv3y,sv3z,sv4x,sv4y,sv4z) !Volume type 3
    tetravol3=tetravol3+gettetravol(tv1x,tv1y,tv1z,tv2x,tv2y,tv2z,sv1x,sv1y,sv1z,sv4x,sv4y,sv4z)
    tetravol3=tetravol3+gettetravol(tv2x,tv2y,tv2z,sv1x,sv1y,sv1z,sv2x,sv2y,sv2z,sv4x,sv4y,sv4z)
end if
    
! if (itestv(1)+itestv(2)==1) call vertexinterpolate(vix(1),viy(1),viz(1),vix(2),viy(2),viz(2),newvtxindt)
! if (itestv(1)+itestv(3)==1) call vertexinterpolate(vix(1),viy(1),viz(1),vix(3),viy(3),viz(3),newvtxindt)
! if (itestv(1)+itestv(4)==1) call vertexinterpolate(vix(1),viy(1),viz(1),vix(4),viy(4),viz(4),newvtxindt)
! if (itestv(2)+itestv(3)==1) call vertexinterpolate(vix(2),viy(2),viz(2),vix(3),viy(3),viz(3),newvtxindt)
! if (itestv(2)+itestv(4)==1) call vertexinterpolate(vix(2),viy(2),viz(2),vix(4),viy(4),viz(4),newvtxindt)
! if (itestv(3)+itestv(4)==1) call vertexinterpolate(vix(3),viy(3),viz(3),vix(4),viy(4),viz(4),newvtxindt)
end subroutine

!!--------- Convert numbering of each vertex of tetrahedron to absolute corner-index
subroutine getvertind(indx,indy,indz,baseix,baseiy,baseiz,inum)
implicit real*8(a-h,o-z)
integer indx,indy,indz,baseix,baseiy,baseiz,inum
indx=baseix
indy=baseiy
indz=baseiz
if (inum==1) then
    indy=baseiy+1
else if (inum==2) then
    indx=baseix+1
    indy=baseiy+1
else if (inum==3) then
    continue
else if (inum==4) then
    indx=baseix+1
    indy=baseiy+1
    indz=baseiz+1
else if (inum==5) then
    indx=baseix+1
else if (inum==6) then
    indx=baseix+1
    indz=baseiz+1
else if (inum==7) then
    indz=baseiz+1
else if (inum==8) then
    indy=baseiy+1
    indz=baseiz+1
end if
end subroutine

!!--------- Interpolate surface vertex from input index of vertex, then save to "survtx" array
!newind is the index of this surface vertex
!ifuncint define use which real space function to do the interpolation
subroutine vertexinterpolate(iax,iay,iaz,ibx,iby,ibz,newind)
use defvar
use surfvertex
use function
implicit real*8(a-h,o-z)
integer iax,iay,iaz,ibx,iby,ibz,newind
icora=abs2suridx(iax,iay,iaz)
icorb=abs2suridx(ibx,iby,ibz)
do i=1,surcor2vtxpos(icora)
    if (surfcor2vtx(icora,i)%athcor==icorb) then !Have already interpolated for these two corner, return interpolated vertex index directly
        newind=surfcor2vtx(icora,i)%itpvtx
        return
    end if
end do

nsurvtx=nsurvtx+1
newind=nsurvtx
surcor2vtxpos(icora)=surcor2vtxpos(icora)+1
surfcor2vtx(icora,surcor2vtxpos(icora))%athcor=icorb
surfcor2vtx(icora,surcor2vtxpos(icora))%itpvtx=nsurvtx
surcor2vtxpos(icorb)=surcor2vtxpos(icorb)+1
surfcor2vtx(icorb,surcor2vtxpos(icorb))%athcor=icora
surfcor2vtx(icorb,surcor2vtxpos(icorb))%itpvtx=nsurvtx

vala=cubmat(iax,iay,iaz)
valb=cubmat(ibx,iby,ibz)

if (nbisec==0) then !Linear interpolate directly
    scl=(vala-surfisoval)/(vala-valb)
    survtx(nsurvtx)%x=(1-scl)*(orgx+(iax-1)*dx)+scl*(orgx+(ibx-1)*dx)
    survtx(nsurvtx)%y=(1-scl)*(orgy+(iay-1)*dy)+scl*(orgy+(iby-1)*dy)
    survtx(nsurvtx)%z=(1-scl)*(orgz+(iaz-1)*dz)+scl*(orgz+(ibz-1)*dz)
else !First perform bisect several times
    ax=orgx+(iax-1)*dx
    ay=orgy+(iay-1)*dy
    az=orgz+(iaz-1)*dz
    bx=orgx+(ibx-1)*dx
    by=orgy+(iby-1)*dy
    bz=orgz+(ibz-1)*dz
    do i=1,nbisec
        bisecx=(ax+bx)/2D0
        bisecy=(ay+by)/2D0
        bisecz=(az+bz)/2D0
        bisecval=calcfuncall(ifuncintp,bisecx,bisecy,bisecz)
        if ((bisecval-surfisoval)*(vala-surfisoval)<0) then
            bx=bisecx
            by=bisecy
            bz=bisecz
            valb=bisecval
        else
            ax=bisecx
            ay=bisecy
            az=bisecz
            vala=bisecval
        end if
    end do
    !Last time, use linear interpolation
    if ((bisecval-surfisoval)*(vala-surfisoval)<0) then !interpolate for a and bisec
        scl=(vala-surfisoval)/(vala-bisecval)
        survtx(nsurvtx)%x=(1-scl)*ax+scl*bisecx
        survtx(nsurvtx)%y=(1-scl)*ay+scl*bisecy
        survtx(nsurvtx)%z=(1-scl)*az+scl*bisecz
    else !interpolate for b and bisec
        scl=(valb-surfisoval)/(valb-bisecval)
        survtx(nsurvtx)%x=(1-scl)*bx+scl*bisecx
        survtx(nsurvtx)%y=(1-scl)*by+scl*bisecy
        survtx(nsurvtx)%z=(1-scl)*bz+scl*bisecz
    end if
end if
end subroutine


!!----------Add a connection between two surface vertices
subroutine addvtxconn(ind1,ind2)
use surfvertex
implicit real*8(a-h,o-z)
integer ind1,ind2
do i=1,vtxconnpos(ind1)
    if (vtxconn(ind1,i)==ind2) return !Have already exist this connection entry
end do
vtxconnpos(ind1)=vtxconnpos(ind1)+1
vtxconn(ind1,vtxconnpos(ind1))=ind2
vtxconnpos(ind2)=vtxconnpos(ind2)+1
vtxconn(ind2,vtxconnpos(ind2))=ind1
end subroutine
