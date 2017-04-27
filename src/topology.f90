!!!============== Find critical point of real space functions
subroutine topomain
use topo
use defvar
use util
use function
implicit real*8(a-h,o-z)
character c200*200,c1000*1000,ctmp1*20,ctmp2*20,ctmp3*20,ctmp4*20,icp1text*12,icp2text*12
real*8,allocatable :: randptx(:),randpty(:),randptz(:) !x,y,z of the points in the sphere
real*8,allocatable :: bassurpathtmp(:,:,:,:) !Used for temporary store bassurpath
real*8,allocatable :: shanCPrho(:) !For calculate shannon aromaticity
integer :: shanCPind(100),searchcenlist(500)
real*8 :: hesstmp(3,3),gradtmp(3) !Temporarily used for calculating curvature

!Initialize searching and plotting parameters
toposphrad=3D0 !Radius of searching sphere is 3 Bohr
numsearchpt=1000 !1000 points randomly scattered in the sphere
sphcenx=0D0 !Position of the sphere center
sphceny=0D0
sphcenz=0D0
!Initialize plot parameter
ZVU=5.0D0 !Closer than other case
idrawisosur=0
ishow3n3=1
ishow3n1=1
ishow3p1=1
ishow3p3=1
idrawpath=1
bondradius=0.07D0 !For showing CPs, we use thiner default bond to avoid overlay the (3,-1)
ratioatmsphere=0.6D0 !Use smaller radius of atom
textheigh=36
ishowCPlab=0
ishowatmlab=0
ishowpathlab=0
ishowsearchlevel=0
ishowattlab=0 !Don't show the result of basin analysis

write(*,*) "         !!! Note: All length units in this module are Bohr !!!"

do while(.true.)
    write(*,*)
    write(*,"(a)") "            ================ Topology analysis ==============="
    write(*,"(a,i3)") " -11 Delete results and reselect real space function, current:",ifunctopo
    write(*,*) "-10 Return"
    write(*,*) "-9 Measure distances, angles and dihedral angles between CPs or atoms"
    write(*,*) "-5 Modify or print detail or export paths, or calculate property along a path"
    write(*,*) "-4 Modify or export CPs (critical points)"
    write(*,*) "-3 Set interbasin surface generating parameters"
    write(*,*) "-2 Set path searching parameters"
    write(*,*) "-1 Set CP searching parameters"
    write(*,*) "0 Print and visualize all generated CPs, paths and surfaces"
    write(*,*) "1 Search CPs from a given starting point"
    write(*,*) "2 Search CPs from nuclear positions"
    write(*,*) "3 Search CPs from midpoint of atom pairs"
    write(*,*) "4 Search CPs from triangle center of three atoms"
    write(*,*) "5 Search CPs from pyramid center of four atoms"
    write(*,*) "6 Search CPs from a batch of points within a sphere" !cube, random in sphere
    write(*,*) "7 Show real space function values at specific CP or all CPs"
    write(*,*) "8 Generate the path connected (3,-3) and (3,-1)"
    write(*,*) "9 Generate the path connected (3,+1) and (3,+3)"
    write(*,*) "10 Add or delete interbasin surfaces"
    if (ifunctopo==1) write(*,*) "20 Calculate Shannon aromaticity index"
    write(*,*) "21 Calculate density curvature perpendicular to a specific plane at a point"
    read(*,*) isel

    if (isel==-11) then
        write(*,*) "0 Return"
        write(*,*) "1 Electron density (Analytical Hessian)"
        write(*,*) "3 Laplacian of electron density"
        write(*,*) "4 Value of orbital wavefunction"
        if (ELFLOL_type==0) write(*,*) "9 Electron localization function(ELF)"
        if (ELFLOL_type==1) write(*,*) "9 Electron localization function(ELF) defined by Tsirelson" 
        if (ELFLOL_type==2) write(*,*) "9 Electron localization function(ELF) defined by Lu, Tian" 
        if (ELFLOL_type==0) write(*,*) "10 Localized orbital locator(LOL)"
        if (ELFLOL_type==1) write(*,*) "10 Localized orbital locator(LOL) defined by Tsirelson" 
        if (ELFLOL_type==2) write(*,*) "10 Localized orbital locator(LOL) defined by Lu, Tian"
        write(*,"(a,i5)") " 100 User defined real space function, iuserfunc=",iuserfunc
!         write(*,*) "12 Total electrostatic potential"
        read(*,*) ifunctopo
        if (ifunctopo==4) then
            write(*,"(a,i10)") " Input orbital index, between 1 and",nmo
            read(*,*) iorbsel
        end if
        if (ifunctopo/=0) then
            numcp=0 !Clean number
            numpath=0
            nple3n1path=0
            numbassurf=0
            CPtype=0 !Clean relationship
            cp2surf=0
            cp2ple3n1path=0
            if (allocated(topopath)) deallocate(topopath)
            if (allocated(bassurpath)) deallocate(bassurpath)
            if (allocated(ple3n1path)) deallocate(ple3n1path)
            write(*,*) "Note: All found CPs, paths, surfaces have been clean"
            !Set special parameters for specific real space functions
            if (ifunctopo==1.or.ifunctopo==4) then !High criteria for full analytical functions
                gradconv=1D-7
                dispconv=1D-8
            else !User lower criteria for those functions with numerical gradient and Hessian
                gradconv=1D-5
                dispconv=1D-6
            end if
            if (ifunctopo==1) then
                toposphrad=3D0
                maxpathpt=451 !Default parameters for rho
                pathstepsize=0.03D0
                numsearchpt=2000
                nsurfpathpercp=60
                nsurfpt=100
                surfpathstpsiz=0.03D0
            else
                toposphrad=3D0
                if (ifunctopo==4) toposphrad=1.5D0 !CPs for orbital wavefunction is very close to nuclei, so use smaller radii
                maxpathpt=901
                pathstepsize=0.015D0 !Curvature is much larger than paths for rho, so smaller differentiate step must be used
                numsearchpt=1000
                if (ifunctopo==4) numsearchpt=10000 !Use large value since evaluation of orbital wavefunction is very fast, but hard to locate CPs
                nsurfpathpercp=200
                nsurfpt=100
                surfpathstpsiz=0.008D0
            end if
        end if
    else if (isel==-10) then
        exit
!-9 -9 -9 -9 -9 -9 -9
    else if (isel==-9) then
        write(*,*) "q = quit"
        write(*,*) "Selection method: a? = Atom?, c? = Critical point ?"
        write(*,*) "e.g. ""a1 c3"" returns the distance between atom1 and CP3"
        write(*,*) "     ""a4 a2"" returns the distance between atom4 and atom2"
        write(*,*) "     ""c6 a2 a5"" returns the angle of CP6-atom2-atom5"
        write(*,*) "     ""c2 c4 a3 c7"" returns the dihedral angle of CP2-CP4-atom3-CP7"
        do while(.true.)
            read(*,"(a)") c200
            c200=adjustl(c200)
            imeasure=0
            do ichar=1,len_trim(c200) !imeasure=1/2/3: measure distance,angle,dihedral
                if (c200(ichar:ichar)==','.or.c200(ichar:ichar)==' ') imeasure=imeasure+1
            end do
            if (c200(1:1)=='q') then
                exit
            else if (imeasure==1.or.imeasure==2.or.imeasure==3) then
                if (imeasure==1) read(c200,*) ctmp1,ctmp2 !Read two terms
                if (imeasure==2) read(c200,*) ctmp1,ctmp2,ctmp3 !Read three terms
                if (imeasure==3) read(c200,*) ctmp1,ctmp2,ctmp3,ctmp4 !Read four terms
                
                if (ctmp1(1:1)=='a') then
                    read(ctmp1(2:),*) iatm
                    tmpx1=a(iatm)%x
                    tmpy1=a(iatm)%y
                    tmpz1=a(iatm)%z
                else if (ctmp1(1:1)=='c') then
                    read(ctmp1(2:),*) icp
                    tmpx1=CPpos(1,icp)
                    tmpy1=CPpos(2,icp)
                    tmpz1=CPpos(3,icp)
                end if
                if (ctmp2(1:1)=='a') then
                    read(ctmp2(2:),*) iatm
                    tmpx2=a(iatm)%x
                    tmpy2=a(iatm)%y
                    tmpz2=a(iatm)%z
                else if (ctmp2(1:1)=='c') then
                    read(ctmp2(2:),*) icp
                    tmpx2=CPpos(1,icp)
                    tmpy2=CPpos(2,icp)
                    tmpz2=CPpos(3,icp)
                end if
                if (imeasure==1) write(*,"(' The distance is',f12.6,' Bohr (',f12.6 ' Angstrom)')") &
                dsqrt((tmpx1-tmpx2)**2+(tmpy1-tmpy2)**2+(tmpz1-tmpz2)**2),dsqrt((tmpx1-tmpx2)**2+(tmpy1-tmpy2)**2+(tmpz1-tmpz2)**2)*b2a

                if (imeasure==2.or.imeasure==3) then !Analyze one more term, then print angle
                    if (ctmp3(1:1)=='a') then
                        read(ctmp3(2:),*) iatm
                        tmpx3=a(iatm)%x
                        tmpy3=a(iatm)%y
                        tmpz3=a(iatm)%z
                    else if (ctmp3(1:1)=='c') then
                        read(ctmp3(2:),*) icp
                        tmpx3=CPpos(1,icp)
                        tmpy3=CPpos(2,icp)
                        tmpz3=CPpos(3,icp)
                    end if
                end if
                if (imeasure==2) write(*,"(' The angle is',f12.6,' degree')") xyz2angle(tmpx1,tmpy1,tmpz1,tmpx2,tmpy2,tmpz2,tmpx3,tmpy3,tmpz3)
                
                if (imeasure==3) then !Analyze one more term, then print dihedral angle
                    if (ctmp4(1:1)=='a') then
                        read(ctmp4(2:),*) iatm
                        tmpx4=a(iatm)%x
                        tmpy4=a(iatm)%y
                        tmpz4=a(iatm)%z
                    else if (ctmp4(1:1)=='c') then
                        read(ctmp4(2:),*) icp
                        tmpx4=CPpos(1,icp)
                        tmpy4=CPpos(2,icp)
                        tmpz4=CPpos(3,icp)
                    end if
                    write(*,"(' The dihedral angle is',f12.6,' degree')") xyz2dih(tmpx1,tmpy1,tmpz1,tmpx2,tmpy2,tmpz2,tmpx3,tmpy3,tmpz3,tmpx4,tmpy4,tmpz4)
                end if
            else
                write(*,*) "Input error"
            end if
        end do
        
!-5 -5 -5 -5 -5 -5 -5
!-5 -5 -5 -5 -5 -5 -5
!-5 -5 -5 -5 -5 -5 -5
    else if (isel==-5) then
        do while(.true.)
            write(*,*)
            write(*,*) "                      ======== Process paths ========"
            write(*,*) "0 Return"
            write(*,*) "1 Print summary of paths"
            write(*,*) "2 Print detail of a path"
            write(*,*) "3 Delete some paths"
            write(*,*) "4 Save points of all paths to paths.txt in current folder"
            write(*,*) "5 Load paths from an external file"
            write(*,*) "6 Export paths as paths.pdb file in current folder"
            write(*,*) "7 Calculate specific real space funtion along the path"
            read(*,*) isel2
            
            if (isel2==0) then
                exit
                
            else if (isel2==1) then
                if (numpath>0) then
                    do i=1,numpath
                        call path_cp(i,icp1,icp2,ipathtype)
                        if (icp1==0) then
                            icp1text="   Unknown  "
                        else
                            write(icp1text,"(i5,1x,a)") icp1,CPtyp2lab(CPtype(icp1))
                        end if
                        if (icp2==0) then
                            icp2text="  Unknown   "
                        else
                            write(icp2text,"(i5,1x,a)") icp2,CPtyp2lab(CPtype(icp2))
                        end if
                        write(*,"('#',i5,5x,'CP:',a,' --->',' CP:',a,'   Length:',f9.5)") i,icp1text,icp2text,(pathnumpt(i)-1)*pathstepsize
                    end do
                else
                    write(*,*) "No paths have been found"
                end if
                
            else if (isel2==2) then
                write(*,*) "Input the index of the path"
                read(*,*) ipath
                if (ipath>numpath.or.ipath<=0) then
                    write(*,*) "Invalid index"
                else
                    write(*,"(a,i6,a,f10.5,a,i5)") "Path:",ipath,"   Length:",(pathnumpt(ipath)-1)*pathstepsize,"   Total points:",pathnumpt(ipath)
                    write(*,"('From',3f18.12,/,'to  ',3f18.12)") topopath(:,1,ipath),topopath(:,pathnumpt(ipath),ipath)
                    write(*,*) "The X/Y/Z coordinate (Bohr) and length of points in the path:"
                    do ipt=1,pathnumpt(ipath)
                        write(*,"(i6,3f16.10,f9.4)") ipt,topopath(:,ipt,ipath),(ipt-1)*pathstepsize
                    end do
                end if
                
            else if (isel2==3) then
                write(*,*) "Input the index range that will be deleted, e.g. 3,10"
                write(*,*) "Note: Input 0,0 can delete all paths"
                read(*,*) idelstart,idelend

                if (idelstart==0.and.idelend==0) then
                    numpath=0
                    if (allocated(topopath)) deallocate(topopath)
                else if ( (idelstart>idelend).or.idelend>numpath.or.idelstart<=0) then
                    write(*,*) "Invalid input"
                else
                    numdel=idelend-idelstart+1
                    nafter=numpath-idelend
                    topopath(:,:,idelstart:idelstart+nafter-1)=topopath(:,:,idelend+1:numpath)
                    pathnumpt(idelstart:idelstart+nafter-1)=pathnumpt(idelend+1:numpath)
                    numpath=numpath-numdel
                    write(*,"(' Now there are',i6,' paths left')") numpath
                end if
                
            else if (isel2==4) then
                open(10,file="paths.txt",status="replace")
                write(10,"(2i10)") numpath,maxpathpt
                do ipath=1,numpath
                    write(10,"(/,'Path index:',i10)") ipath
                    write(10,"(i10)") pathnumpt(ipath)
                    do ipt=1,pathnumpt(ipath)
                        write(10,"(3E20.12)") topopath(:,ipt,ipath)
                    end do
                end do
                close(10)
                write(*,*) "Done, path information have been saved to paths.txt in current folder"
                write(*,"(a)") " Units are in Bohr. The first two numbers respectively denote the number of paths, and the maximal number of points in the paths"
                
            else if (isel2==5) then
                write(*,"(a)") " Note: The format of the input file must be identical to the one outputted by option 4"
                if (numpath>0) write(*,*) "Note: After loading the file, all current paths will be clean"
                write(*,*) "Input filename, e.g. c:\paths.txt"
                read(*,*) c200
                inquire(file=c200,exist=alive)
                if (alive.eqv..false.) then
                    write(*,*) "File not found"
                else
                    open(10,file=c200,status="old")
                    read(10,*) numpath,maxpathpt
                    if (allocated(topopath)) deallocate(topopath)
                    allocate(topopath(3,maxpathpt,numpath))
                    do ipath=1,numpath
                        read(10,*)
                        read(10,*)
                        read(10,*) pathnumpt(ipath)
                        do ipt=1,pathnumpt(ipath) 
                            read(10,*) topopath(:,ipt,ipath)
                        end do
                    end do
                    close(10)
                    write(*,*) "Done, path information have been recovered from the file"
                end if
                
            else if (isel2==6) then
                open(10,file="paths.pdb",status="replace")
                itmp=0
                do ipath=1,numpath
                    do ipt=1,pathnumpt(ipath)
                        itmp=itmp+1
                        write(10,"(a6,i5,1x,a4,1x,a3, 1x,a1,i4,4x,3f8.3,2f6.2,10x,a2)") "HETATM",itmp,' '//"C "//' ',"PTH",'A',ipath,topopath(:,ipt,ipath)*b2a,1.0,0.0,"C "
                    end do
                    write(10,"('TER')")
                end do
                close(10)
                write(*,*) "Done, path information have been saved to paths.pdb in current folder"
                
            else if (isel2==7) then
                do while(.true.)
                    write(*,*)
                    write(*,*) "Input index of a path, e.g. 3"
                    write(*,*) "Input ""q"" can return"
                    write(*,"(a)") " Hint: If input index of two paths (e.g. 6,7) emitted from the same (3,-1) CP, &
                    then the real space function along the combined paths will be outputted"
                    if (allocated(MOsym)) write(*,"(a)") " Hint: You can input ""s"" to choose which irreducible representations will be taken into account in the real space function evaluation"
                    read(*,"(a)") c200
                    
                    if (index(c200,'q')/=0) then
                        exit
                    else if (index(c200,'s')/=0) then
                        call SelMO_IRREP
                        cycle
                    end if
                    
                    itwopath=0
                    if (index(c200,',')==0) then !Only one path
                        read(c200,*) ipath
                    else !Two paths
                        itwopath=1
                        read(c200,*) ipath,jpath
                    end if
                    if (ipath>numpath.or.ipath<=0 .or. (itwopath==1.and.(jpath>numpath.or.jpath<=0)) ) then
                        write(*,*) "Error: Invalid index"
                        cycle
                    end if
                    if (itwopath==1.and.( topopath(1,1,ipath)/=topopath(1,1,jpath) .or. topopath(2,1,ipath)/=topopath(2,1,jpath) .or. topopath(3,1,ipath)/=topopath(3,1,jpath) )) then
                        write(*,*) "Error: The two paths are not emitted from the same (3,-1) critical point!"
                        cycle
                    end if
                    write(*,*) "Select the real space function to be calculated along the path"
                    call selfunc_interface(iselfunc)
                    !Store the values along the path into curvex and curvey, show them as text and curve map
                    npointcurve=pathnumpt(ipath)
                    if (itwopath==1) npointcurve=pathnumpt(ipath)+pathnumpt(jpath)-1
                    if (allocated(curvex)) deallocate(curvex)
                    if (allocated(curvey)) deallocate(curvey)
                    allocate(curvex(npointcurve),curvey(npointcurve))
                    if (itwopath==0) then
                        do ipt=1,pathnumpt(ipath)
                            curvex(ipt)=(ipt-1)*pathstepsize
                            curvey(ipt)=calcfuncall(iselfunc,topopath(1,ipt,ipath),topopath(2,ipt,ipath),topopath(3,ipt,ipath))
                        end do
                    else if (itwopath==1) then
                        ipttmp=0
                        do ipt=pathnumpt(ipath),1,-1
                            ipttmp=ipttmp+1
                            curvex(ipttmp)=(ipttmp-1)*pathstepsize
                            curvey(ipttmp)=calcfuncall(iselfunc,topopath(1,ipt,ipath),topopath(2,ipt,ipath),topopath(3,ipt,ipath))
                        end do
                        do jpt=2,pathnumpt(jpath) !The first point of jpath is (3,-1), which has been included in ipath above
                            ipttmp=ipttmp+1
                            curvex(ipttmp)=(ipttmp-1)*pathstepsize
                            curvey(ipttmp)=calcfuncall(iselfunc,topopath(1,jpt,jpath),topopath(2,jpt,jpath),topopath(3,jpt,jpath))
                        end do
                    end if    
                    write(*,*) " Index           X/Y/Z Coordinate (Bohr)          Dist.       Value"
                    if (itwopath==0) then
                        do ipt=1,pathnumpt(ipath)
                            write(*,"(i6,3f14.8,f9.4,E16.8)") ipt,topopath(:,ipt,ipath),curvex(ipt),curvey(ipt)
                        end do
                        icurve_vertlinex=0
                    else if (itwopath==1) then
                        ipttmp=0
                        do ipt=pathnumpt(ipath),1,-1
                            ipttmp=ipttmp+1
                            write(*,"(i6,3f14.8,f9.4,E16.8)") ipttmp,topopath(:,ipt,ipath),curvex(ipttmp),curvey(ipttmp)
                        end do
                        ibcptmp=ipttmp
                        do jpt=2,pathnumpt(jpath) !The first point of jpath is (3,-1), which has been included in ipath above
                            ipttmp=ipttmp+1
                            write(*,"(i6,3f14.8,f9.4,E16.8)") ipttmp,topopath(:,jpt,jpath),curvex(ipttmp),curvey(ipttmp)
                        end do
                        write(*,"(' Note: The point',i5,' corresponds to (3,-1) critical point')") ibcptmp
                        icurve_vertlinex=1 !Draw a vertical line to highlight BCP point
                        curve_vertlinex=curvex(ibcptmp)
                        write(*,*) "The dash line corresponds to the position of BCP"
                    end if
                    !Show the result as curve map
                    if (minval(curvey)>=0) then
                        curveymin=0
                        curveymax=1.1D0*maxval(curvey)
                        steplaby=curveymax/11
                    else
                        exty=(maxval(curvey)-minval(curvey))/10
                        curveymin=minval(curvey)-exty
                        curveymax=maxval(curvey)+exty
                        steplaby=exty
                    end if
                    steplabx=0.5D0
                    ilog10y=0
                    do while(.true.)
                        write(*,*)
                        write(*,*) "0 Return"
                        write(*,*) "1 Plot the graph again"
                        write(*,*) "2 Save the graph in current folder"
                        if (ilog10y==0) then
                            write(*,"(' 3 Set range of Y-axis, current:',f13.5,' to',f13.5,' Step:',f12.5)") curveymin,curveymax,steplaby
                            write(*,*) "4 Switch Y-axis type, current: linear scaling"
                        else if (ilog10y==1) then
                            write(*,"(' 3 Set range of Y-axis, current:',1PE14.6,' to',1PE14.6)") 10**curveymin,10**curveymax
                            write(*,*) "4 Switch Y-axis type, current: logarithmic scaling"
                        end if
                        write(*,*) "5 Export the data of the path shown above as pathvalue.txt in current folder"
                        read(*,*) iselcurve
                        
                        if (iselcurve==0) then
                            exit
                        else if (iselcurve==1) then
                        else if (iselcurve==2) then
                        else if (iselcurve==3) then
                            if (ilog10y==0) then
                                write(*,*) "Input lower and upper limits and step between labels, e.g. 0,1.5,0.2"
                                read(*,*) curveymin,curveymax,steplaby
                            else if (ilog10y==1) then
                                write(*,*) "Input minimum and maximum value of Y axis  e.g. -1,5 means from 10^-1 to 10^5"
                                read(*,*) curveymin,curveymax
                            end if
                        else if (iselcurve==4) then
                            if (ilog10y==0) then
                                ilog10y=1
                                curveymin=log10(curveymin)
                                curveymax=log10(curveymax)
                            else
                                ilog10y=0
                                curveymin=10**curveymin
                                curveymax=10**curveymax
                            end if
                        else if (iselcurve==5) then                
                            open(10,file="pathvalue.txt",status="replace")
                            if (itwopath==0) then
                                do ipt=1,pathnumpt(ipath)
                                    write(10,"(i6,3f14.8,f9.4,E16.8)") ipt,topopath(:,ipt,ipath),curvex(ipt),curvey(ipt)
                                end do
                            else if (itwopath==1) then
                                ipttmp=0
                                do ipt=pathnumpt(ipath),1,-1
                                    ipttmp=ipttmp+1
                                    write(10,"(i6,3f14.8,f9.4,E16.8)") ipttmp,topopath(:,ipt,ipath),curvex(ipttmp),curvey(ipttmp)
                                end do
                                do jpt=2,pathnumpt(jpath) !The first point of jpath is (3,-1), which has been included in ipath above
                                    ipttmp=ipttmp+1
                                    write(10,"(i6,3f14.8,f9.4,E16.8)") ipttmp,topopath(:,jpt,jpath),curvex(ipttmp),curvey(ipttmp)
                                end do
                            end if
                            close(10)
                            write(*,*) "Done! The data has been outputted to pathvalue.txt in current folder"
                        end if
                    end do
                    icurve_vertlinex=0
                    deallocate(curvex,curvey)
                end do
            end if
        end do
        
!-4 -4 -4 -4 -4 -4 -4
!-4 -4 -4 -4 -4 -4 -4
!-4 -4 -4 -4 -4 -4 -4
    else if (isel==-4) then
        do while(.true.)
            write(*,*) "               ============ Modify or export found CPs ============"
            write(*,*) "-1 Print summary of CPs in Angstrom"
            write(*,*) "0 Return"
            write(*,*) "1 Print summary of CPs (in Bohr)"
            write(*,*) "2 Delete some CPs"
            write(*,*) "3 Add a CP artificially"
            write(*,*) "4 Save CPs to CPs.txt in current folder"
            write(*,*) "5 Load CPs from a file"
            write(*,*) "6 Export CPs as CPs.pdb file in current folder"
            read(*,*) isel2
            
            if (isel2==0) then
                exit
            else if (isel2==1.or.isel2==-1) then
                if (numcp>0) then
                    write(*,*) "Summary of found CPs:"
                    write(*,*) " Index                    Coordinate                     Type"
                    do icp=1,numcp
                        if (isel2==1) write(*,"(i6,3f16.9,3x,a)") icp,CPpos(:,icp),CPtyp2lab(CPtype(icp))
                        if (isel2==-1) write(*,"(i6,3f16.9,3x,a)") icp,CPpos(:,icp)*b2a,CPtyp2lab(CPtype(icp))
                    end do
                else
                    write(*,*) "No CPs have been found"
                end if
            else if (isel2==2) then
                if (numbassurf>0) then
                    write(*,"(a)") "Warning: If one or more CPs are deleted, all interbasin surfaces will be lost,  continue? 0/1=No/Yes"
                    read(*,*) ifok
                    if (ifok==0) then
                        cycle
                    else
                        numbassurf=0
                        cp2surf=0
                        deallocate(bassurpath)
                    end if
                end if
                write(*,*) "Input the index range that will be deleted, e.g. 3,10"
                write(*,*) "Note 1: Input 0,0 can delete all CPs"
                write(*,*) "Note 2: If you want to delete CPs with too low electron density, select -1"
                read(*,*) idelstart,idelend
                if (idelstart==0.and.idelend==0) then
                    numcp=0
                else if ( (idelstart>idelend).or.idelend>numcp.or.idelstart<=0) then
                    write(*,*) "Invalid input"
                else
                    numdel=idelend-idelstart+1
                    nafter=numcp-idelend
                    CPpos(:,idelstart:idelstart+nafter-1)=CPpos(:,idelend+1:numcp)
                    CPtype(idelstart:idelstart+nafter-1)=CPtype(idelend+1:numcp)
                    numcp=numcp-numdel
                end if
                write(*,"(' Now there are',i6,' CPs left')") numcp
            else if (isel2==3) then
                numcp=numcp+1
                write(*,*) "Input coordinate, e.g. 1.2,0.0,3.44"
                read(*,*) CPpos(1,numcp),CPpos(2,numcp),CPpos(3,numcp)
                write(*,*) "Input CP type, 1=(3,-3) 2=(3,-1) 3=(3,+1) 4=(3,+3)"
                read(*,*) CPtype(numcp)
            else if (isel2==4) then
                open(10,file="CPs.txt",status="replace")
                write(10,"(i6)") numcp
                do icp=1,numcp
                    write(10,"(i6,3f12.6,3x,i4)") icp,CPpos(:,icp),CPtype(icp)
                end do
                close(10)
                write(*,*) "Done, CP information have been saved to CPs.txt in current folder"
                write(*,*) "Note: The last column is CP type, 1=(3,-3) 2=(3,-1) 3=(3,+1) 4=(3,+3)"
            else if (isel2==5) then
                write(*,*) "Input filename, e.g. C:\ltwd\CPs.txt"
                write(*,"(a)") " (The format of the file must be identical to the one outputted by option 4)"
                if (numcp>0) write(*,*) "Note: After loading the file, all found CPs will be clean"
                read(*,*) c200
                inquire(file=c200,exist=alive)
                if (alive.eqv..false.) then
                    write(*,*) "File not found"
                else
                    open(10,file=c200,status="old")
                    read(10,*) numcp
                    do icp=1,numcp
                        read(10,*) nouse,CPpos(:,icp),CPtype(icp)
                    end do
                    close(10)
                    write(*,*) "Done, CP information have been recovered from the file"
                end if
            else if (isel2==6) then
                open(10,file="CPs.pdb",status="replace")
                do icp=1,numcp
                    if (CPtype(icp)==1) write(10,"(a6,i5,1x,a4,1x,a3, 1x,a1,i4,4x,3f8.3,2f6.2,10x,a2)") "HETATM",icp,' '//"C "//' ',"CPS",'A',1,CPpos(:,icp)*b2a,1.0,0.0,"C "
                    if (CPtype(icp)==2) write(10,"(a6,i5,1x,a4,1x,a3, 1x,a1,i4,4x,3f8.3,2f6.2,10x,a2)") "HETATM",icp,' '//"N "//' ',"CPS",'A',1,CPpos(:,icp)*b2a,1.0,0.0,"N "
                    if (CPtype(icp)==3) write(10,"(a6,i5,1x,a4,1x,a3, 1x,a1,i4,4x,3f8.3,2f6.2,10x,a2)") "HETATM",icp,' '//"O "//' ',"CPS",'A',1,CPpos(:,icp)*b2a,1.0,0.0,"O "
                    if (CPtype(icp)==4) write(10,"(a6,i5,1x,a4,1x,a3, 1x,a1,i4,4x,3f8.3,2f6.2,10x,a2)") "HETATM",icp,' '//"F "//' ',"CPS",'A',1,CPpos(:,icp)*b2a,1.0,0.0,"F "
                end do
                write(*,*) "Done, CP information have been saved to CPs.pdb in current folder"
                write(*,*) "Note: Element C/N/O/F correspond to (3,-3)/(3,-1)/(3,+1)/(3,+3) respectively"
                close(10)
            end if
            write(*,*)
        end do
!-3 -3 -3 -3 -3 -3 -3
    else if (isel==-3) then
        do while(.true.)
            write(*,*) "============ Set interbasin surface generating parameters ============"
            write(*,"(a)")      " 0 Return"
            write(*,"(a,i6)")   " 1 Number of paths in each interbasin surface",nsurfpathpercp
            write(*,"(a,i6)")   " 2 Number of points in each interbasin surface path, current:",nsurfpt
            write(*,"(a,f8.4)") " 3 Stepsize, current:",surfpathstpsiz
            read(*,*) isel2

            if (isel2==0) then
                exit
            else
                if (allocated(bassurpath)) then
                    write(*,*) "If the parameter is changed, all already generated surfaces will be clean, OK?"
                    write(*,*) "0/1=No/Yes"
                    read(*,*) ifok
                    if (ifok==1) then
                        deallocate(bassurpath)
                        numbassurf=0
                        CP2surf=0
                    end if
                end if
                if (.not.allocated(bassurpath)) then
                    write(*,*) "Input the value"
                    if (isel2==1) read(*,*) nsurfpathpercp
                    if (isel2==2) read(*,*) nsurfpt
                    if (isel2==3) read(*,*) surfpathstpsiz
                end if
            end if
        end do
!-2 -2 -2 -2 -2 -2 -2
    else if (isel==-2) then
        do while(.true.)
            write(*,*) "           ============ Set path searching parameters ============"
            write(*,"(a)")      " 0 Return"
            write(*,"(a,i4)")   " 1 Maximal number of points of a path, current:",maxpathpt
            write(*,"(a,f8.4)") " 2 Stepsize, current:",pathstepsize
            write(*,"(a,f8.4)") " 3 Stop generation if distance to any CP is smaller than:",discritpathfin
!             write(*,"(a,i4)")   " 4 Maximal number of bisection, current:",npathtry
            read(*,*) isel2

            if (isel2==0) then
                exit
            else if (isel2==1.or.isel2==2) then
                iok=1
                if (numpath>0) then
                    write(*,*) "Warning: All found paths will be clean, OK? 1=Yes 2=No"
                    read(*,*) ifok
                end if
                if (iok==1) then
                    write(*,*) "Input a value"
                    if (isel2==1) read(*,*) maxpathpt
                    if (isel2==2) read(*,*) pathstepsize
                    numpath=0
                    if (allocated(topopath)) deallocate(topopath)
                end if
            else if (isel2==3) then
                write(*,*) "Input a value"
                read(*,*) discritpathfin
            else if (isel2==4) then
                write(*,*) "Input an integer, must >=2"
                read(*,*) npathtry
            end if
        end do
!-1 -1 -1 -1 -1 -1 -1
    else if (isel==-1) then
        do while(.true.)
            write(*,*)"              ============ Set CP searching parameters ============"
            write(*,"(a)") " -1 Set to default"
            write(*,"(a)") " 0 Return"
            write(*,"(a,i5)") " 1 Set maximal iterations:",topomaxcyc
            write(*,"(a,f12.6)") " 2 Set scale factor for stepsize:",CPstepscale
            write(*,"(a,1PE12.5)") " 3 Criteria for gradient-norm convergence:",gradconv
            write(*,"(a,1PE12.5)") " 4 Criteria for displacement convergence:",dispconv
            write(*,"(a,f12.6)") " 5 Minimal distance between CPs:",minicpdis
            write(*,"(a,f8.2)") " 6 Skip search if distance between atoms is longer than the sum of their vdW radius multiplied by:",vdwsumcrit
            if (ishowsearchlevel==0) write(*,"(a)") " 7 If print details of CP searching procedure: No"
            if (ishowsearchlevel==1) write(*,"(a)") " 7 If print details of CP searching procedure: Some detail"
            if (ishowsearchlevel==2) write(*,"(a)") " 7 If print details of CP searching procedure: All detail"
            write(*,"(a,1PE15.8)") " 8 Criteria for determining if Hessian matrix is singular:",singularcrit
            if (CPsearchlow==CPsearchhigh) then
                write(*,*) "9 Select rule for reserving CPs, current: reserve All CPs"
            else
                write(*,"(a,1PE12.4,a,1PE12.4)") " 9 Select rule for reserving CPs, current: between",CPsearchlow,' and',CPsearchhigh
            end if
            read(*,*) isel2

            if (isel2==-1) then
                topomaxcyc=120
                CPstepscale=1D0
                gradonv=1D-7
                dispconv=1D-9
                minicpdis=0.03D0
                vdwsumcrit=1.2D0
                ishowsearchlevel=0
            else if (isel2==0) then
                exit
            else if (isel2==1) then
                write(*,*) "Input an integer"
                read(*,*) topomaxcyc
            else if (isel2==2) then
                write(*,*) "Input a value, default is 1.0"
                read(*,*) CPstepscale
            else if (isel2==3) then
                write(*,*) "Input a value"
                read(*,*) gradconv
            else if (isel2==4) then
                write(*,*) "Input a value"
                read(*,*) dispconv
            else if (isel2==5) then
                write(*,*) "Input a value"
                read(*,*) minicpdis
            else if (isel2==6) then
                write(*,*) "Input a value"
                read(*,*) vdwsumcrit
            else if (isel2==7) then
                write(*,*) "0 Don't print details"
                write(*,*) "1 Print some details"
                write(*,*) "2 Print all details"
                read(*,*) ishowsearchlevel
                if ( nthreads  >1) write(*,*) "Warning: The printed details may be messed up since parallel mode is enabled!"
            else if (isel2==8) then
                write(*,"(a)") "Input a value, if absolute value of determinant of Hessiant matrix is lower than this value, then it will be regarded as singular, e.g. 1D-21"
                read(*,*) singularcrit
            else if (isel2==9) then
                write(*,"(a)") " Input lower and upper limits. For example, if you input 0.05,0.22, &
                then during CP searching, only when the real space function of a new CP is between 0.05 and 0.22 then it will be reserved"
                write(*,*) "Note: If the two values are identical, all CPs will be reserved (default case)"
                read(*,*) CPsearchlow,CPsearchhigh
            end if
        end do
!0000000000000000
    else if (isel==0) then
        if (numpath>0) then
            write(*,*) "Summary of found paths:"
            do i=1,numpath
                call path_cp(i,icp1,icp2,ipathtype)
                if (icp1==0) then
                    icp1text="   Unknown  "
                else
                    write(icp1text,"(i5,1x,a)") icp1,CPtyp2lab(CPtype(icp1))
                end if
                if (icp2==0) then
                    icp2text="  Unknown   "
                else
                    write(icp2text,"(i5,1x,a)") icp2,CPtyp2lab(CPtype(icp2))
                end if
                write(*,"('#',i5,5x,'CP:',a,' --->',' CP:',a,'   Length:',f9.5)") i,icp1text,icp2text,(pathnumpt(i)-1)*pathstepsize
            end do
        else
            write(*,*) "No paths have been found"
        end if
        write(*,*)
        if (numcp>0) then
            write(*,*) "Summary of found CPs:"
            write(*,*) " Index               XYZ Coordinate (Bohr)               Type"
            do icp=1,numcp
                icptype=CPtype(icp)
                if (ifunctopo==1.and.icptype==1) then
                    do iatm=1,ncenter
                        disttmp2=(CPpos(1,icp)-a(iatm)%x)**2+(CPpos(2,icp)-a(iatm)%y)**2+(CPpos(3,icp)-a(iatm)%z)**2
                        if (disttmp2<0.01D0) then
                            write(*,"(i6,3f16.9,3x,a,'  Nuc:',i5,'(',a')')") icp,CPpos(:,icp),CPtyp2lab(icptype),iatm,a(iatm)%name
                            exit
                        end if
                    end do
                    if (iatm==ncenter+1) write(*,"(i6,3f16.9,3x,a,'  Nuc:   Unknown')") icp,CPpos(:,icp),CPtyp2lab(icptype)
                else
                    write(*,"(i6,3f16.9,3x,a)") icp,CPpos(:,icp),CPtyp2lab(icptype)
                end if
            end do
        else
            write(*,*) "No CPs have been found"
        end if
        NumCPtype1=count(CPtype(1:numcp)==1)
        NumCPtype2=count(CPtype(1:numcp)==2)
        NumCPtype3=count(CPtype(1:numcp)==3)
        NumCPtype4=count(CPtype(1:numcp)==4)
        write(*,*) "The number of critical points of each type:"
        write(*,"(' (3,-3):',i6,',   (3,-1):',i6,',   (3,+1):',i6,',   (3,+3):',i6)") NumCPtype1,NumCPtype2,NumCPtype3,NumCPtype4
        
        itestPH=NumCPtype1-NumCPtype2+NumCPtype3-NumCPtype4 !Poincare-Hopf relationship
        write(*,"(' Poincare-Hopf relationship verification:',i5,'  -',i5,'  +',i5,'  -',i5,'  =',i4)") NumCPtype1,NumCPtype2,NumCPtype3,NumCPtype4,itestPH
        if (itestPH/=1) write(*,*) "Warning: Poincare-Hopf relationship is not satisfied, some CPs may be missing"
        if (itestPH==1) write(*,*) "Fine, Poincare-Hopf relationship is satisfied, all CPs may have been found"
        
        if (numbassurf>0) write(*,"(' The number of generated interbasin surfaces:',i8)") numbassurf
        if (numpath>0) idrawmol=0 !Avoid atom and bond covered paths
!111111111111111111111
    else if (isel==1) then
        numcpold=numcp
        write(*,*) "Input X,Y,Z of starting point (in bohr, e.g. 2.0,3.1,-0.5)"
        write(*,"(a)") " You can also input two atomic indices (e.g. 4,5), then midpoint of corresponding two atoms will be taken as starting point"
        read(*,"(a)") c200
        read(c200,*,iostat=ierror) x,y,z
        if (ierror/=0) then
            read(c200,*) iatm,jatm
            x=(a(iatm)%x+a(jatm)%x)/2
            y=(a(iatm)%y+a(jatm)%y)/2
            z=(a(iatm)%z+a(jatm)%z)/2
        end if
        call findcp(x,y,z,ifunctopo,0)
        if (numcp==numcpold) then
            write(*,*) "No new critical point was found"
        else
            write(*,"(' Find',i5,' new critical points')") numcp-numcpold
        end if
!222222222222222222222
    else if (isel==2) then
        numcpold=numcp
        do iatm=1,ncenter
            write(*,"('#',i5,' /',i5,a,i5,'(',a,')')") iatm,ncenter,": Trying from nuclear position of ",iatm,a(iatm)%name
            if (a(iatm)%index>10) then
                call findcp(a(iatm)%x,a(iatm)%y,a(iatm)%z,ifunctopo,1) !For heavy atoms, use lower criteria, because the cusp of electron density is sharp so hard to locate
            else
                call findcp(a(iatm)%x,a(iatm)%y,a(iatm)%z,ifunctopo,0)
            end if
        end do
        if ((numcp-numcpold)/=0) then
            call sortCP(numcpold+1)
            write(*,*) "                            ==== Summary ===="
            write(*,*) " Index              Coordinate               Type"
            do icp=numcpold+1,numcp
                write(*,"(i6,3f12.6,3x,a)") icp,CPpos(:,icp),CPtyp2lab(CPtype(icp))
            end do
        end if
        write(*,"(' Totally find',i6,' new critical points')") numcp-numcpold
        if (ifunctopo==1.and.count(CPtype(1:numcp)==1)<ncenter) write(*,*) "Warning: Some (3,-3) may missing, try to search again with different parameters"
!333333333333333333333
    else if (isel==3) then
        numcpold=numcp
        itime=0
        ntime=0
        do iatm=1,ncenter !Test how many iterations will be done
            do jatm=iatm+1,ncenter
                if ( distmat(iatm,jatm) <= vdwsumcrit*(vdwr(a(iatm)%index)+vdwr(a(jatm)%index)) ) ntime=ntime+1
            end do
        end do
nthreads=getNThreads()
!$OMP PARALLEL DO SHARED(itime) PRIVATE(iatm,jatm) schedule(dynamic) NUM_THREADS(nthreads)
        do iatm=1,ncenter
            do jatm=iatm+1,ncenter
                if ( distmat(iatm,jatm) > vdwsumcrit*(vdwr(a(iatm)%index)+vdwr(a(jatm)%index)) ) cycle
!$OMP CRITICAL
                    itime=itime+1
!$OMP end CRITICAL
                write(*,"('#',i5,' /',i5,a,i5,'(',a,')',a,i5,'(',a,')')") &
                itime,ntime,": Trying from midpoint between ",iatm,a(iatm)%name," and",jatm,a(jatm)%name
                call findcp( (a(iatm)%x+a(jatm)%x)/2D0,(a(iatm)%y+a(jatm)%y)/2D0,(a(iatm)%z+a(jatm)%z)/2D0, ifunctopo,0)
            end do
        end do    
!$OMP END PARALLEL DO
        if ((numcp-numcpold)/=0) then
            call sortCP(numcpold+1)
            write(*,*) "                            ==== Summary ===="
            write(*,*) " Index              Coordinate               Type"
            do icp=numcpold+1,numcp
                write(*,"(i6,3f12.6,3x,a)") icp,CPpos(:,icp),CPtyp2lab(CPtype(icp))
            end do
        end if
        write(*,"(' Totally find',i6,' new critical points')") numcp-numcpold
!4444444444444444444
    else if (isel==4) then
        numcpold=numcp
        itime=0
        ntime=0
        do iatm=1,ncenter !Test how many iterations will be done
            do jatm=iatm+1,ncenter
                if ( distmat(iatm,jatm) > vdwsumcrit*(vdwr(a(iatm)%index)+vdwr(a(jatm)%index)) ) cycle
                do katm=jatm+1,ncenter
                    if ( distmat(katm,iatm) > vdwsumcrit*(vdwr(a(katm)%index)+vdwr(a(iatm)%index)).or.&
                    distmat(katm,jatm) > vdwsumcrit*(vdwr(a(katm)%index)+vdwr(a(jatm)%index)) ) cycle
                    ntime=ntime+1
                end do
            end do
        end do
nthreads=getNThreads()
!$OMP PARALLEL DO SHARED(itime) PRIVATE(iatm,jatm,katm) schedule(dynamic) NUM_THREADS(nthreads)
        do iatm=1,ncenter
            do jatm=iatm+1,ncenter
                if ( distmat(iatm,jatm) > vdwsumcrit*(vdwr(a(iatm)%index)+vdwr(a(jatm)%index)) ) cycle
                do katm=jatm+1,ncenter
                    if ( distmat(katm,iatm) > vdwsumcrit*(vdwr(a(katm)%index)+vdwr(a(iatm)%index)).or.&
                    distmat(katm,jatm) > vdwsumcrit*(vdwr(a(katm)%index)+vdwr(a(jatm)%index)) ) cycle
!$OMP CRITICAL
                    itime=itime+1
!$OMP end CRITICAL
                    write(*,"('#',i5,' /',i5,a ,i5,'(',a,')' ,i5,'(',a,')' ,i5,'(',a,')' )") &
                    itime,ntime,": Trying from triangle center of ",iatm,a(iatm)%name,jatm,a(jatm)%name,katm,a(katm)%name
                    call findcp( (a(iatm)%x+a(jatm)%x+a(katm)%x)/3D0,(a(iatm)%y+a(jatm)%y+a(katm)%y)/3D0,(a(iatm)%z+a(jatm)%z+a(katm)%z)/3D0, ifunctopo,0)
                end do
            end do
        end do
!$OMP END PARALLEL DO
        if ((numcp-numcpold)/=0) then
            call sortCP(numcpold+1)
            write(*,*) "                            ==== Summary ===="
            write(*,*) " Index              Coordinate               Type"
            do icp=numcpold+1,numcp
                write(*,"(i6,3f12.6,3x,a)") icp,CPpos(:,icp),CPtyp2lab(CPtype(icp))
            end do
        end if
        write(*,"(' Totally find',i6,' new critical points')") numcp-numcpold
!5555555555555555555
    else if (isel==5) then
        numcpold=numcp
        itime=0 
        ntime=0
        do iatm=1,ncenter !Test how many iterations will be done; ij,jk,kl,li,lj,ik
            do jatm=iatm+1,ncenter
                if ( distmat(iatm,jatm) > vdwsumcrit*(vdwr(a(iatm)%index)+vdwr(a(jatm)%index)) ) cycle
                do katm=jatm+1,ncenter
                    if ( distmat(katm,jatm) > vdwsumcrit*(vdwr(a(katm)%index)+vdwr(a(jatm)%index)) ) cycle
                    do latm=katm+1,ncenter
                        if ( distmat(latm,katm) > vdwsumcrit*(vdwr(a(latm)%index)+vdwr(a(katm)%index)).or.&
                        distmat(latm,iatm) > vdwsumcrit*(vdwr(a(latm)%index)+vdwr(a(iatm)%index)).or.&
                        distmat(latm,jatm) > vdwsumcrit*(vdwr(a(latm)%index)+vdwr(a(jatm)%index)).or.&
                        distmat(iatm,katm) > vdwsumcrit*(vdwr(a(iatm)%index)+vdwr(a(katm)%index))) cycle
                        ntime=ntime+1
                    end do
                end do
            end do
        end do
nthreads=getNThreads()
!$OMP PARALLEL DO SHARED(itime) PRIVATE(iatm,jatm,katm,latm) schedule(dynamic) NUM_THREADS(nthreads)
        do iatm=1,ncenter !Test how many iterations will be done; ij,jk,kl,li,lj,ik
            do jatm=iatm+1,ncenter
                if ( distmat(iatm,jatm) > vdwsumcrit*(vdwr(a(iatm)%index)+vdwr(a(jatm)%index)) ) cycle
                do katm=jatm+1,ncenter
                    if ( distmat(katm,jatm) > vdwsumcrit*(vdwr(a(katm)%index)+vdwr(a(jatm)%index)) ) cycle
                    do latm=katm+1,ncenter
                        if ( distmat(latm,katm) > vdwsumcrit*(vdwr(a(latm)%index)+vdwr(a(katm)%index)).or.&
                        distmat(latm,iatm) > vdwsumcrit*(vdwr(a(latm)%index)+vdwr(a(iatm)%index)).or.&
                        distmat(latm,jatm) > vdwsumcrit*(vdwr(a(latm)%index)+vdwr(a(jatm)%index)).or.&
                        distmat(iatm,katm) > vdwsumcrit*(vdwr(a(iatm)%index)+vdwr(a(katm)%index))) cycle
!$OMP CRITICAL
                        itime=itime+1
!$OMP end CRITICAL
                        write(*,"('#',i5,' /',i5,a ,i5,'(',a,')' ,i5,'(',a,')' ,i5,'(',a,')' ,i5,'(',a,')')") &
                        itime,ntime,": Trying from center of ",&
                        iatm,a(iatm)%name,jatm,a(jatm)%name,katm,a(katm)%name,latm,a(latm)%name
                        call findcp( (a(iatm)%x+a(jatm)%x+a(katm)%x+a(latm)%x)/4D0,&
                        (a(iatm)%y+a(jatm)%y+a(katm)%y+a(latm)%y)/4D0,(a(iatm)%z+a(jatm)%z+a(katm)%z+a(latm)%z)/4D0, ifunctopo,0)
                    end do
                end do
            end do
        end do
!$OMP END PARALLEL DO
        if ((numcp-numcpold)/=0) then
            call sortCP(numcpold+1)
            write(*,*) "                            ==== Summary ===="
            write(*,*) " Index              Coordinate               Type"
            do icp=numcpold+1,numcp
                write(*,"(i6,3f12.6,3x,a)") icp,CPpos(:,icp),CPtyp2lab(CPtype(icp))
            end do
        end if
        write(*,"(' Totally found',i6,' new critical points')") numcp-numcpold
!6666666666666666666
    else if (isel==6) then
        do while(.true.)
            write(*,*) "         ============= Set the starting points ============="
            write(*,"(' Center:',3f10.5,' Radius:',f6.2,' Points:',i8)") sphcenx,sphceny,sphcenz,toposphrad,numsearchpt
            write(*,"(a)") " -9 Return"
            write(*,"(a)") " -2 Start the search using some nuclei as sphere center in turn"
            write(*,"(a)") " -1 Start the search using each nucleus as sphere center in turn"
            write(*,"(a)") " 0 Start the search using the defined sphere center"
            write(*,"(a)") " 1 Input coordinate of the sphere center"
            write(*,"(a)") " 2 Set the sphere center at a nuclear position"
            write(*,"(a)") " 3 Set the sphere center at midpoint between two atoms"
            write(*,"(a)") " 4 Set the sphere center at triangle center of three atoms"
            write(*,"(a)") " 5 Set the sphere center at a CP"
            write(*,"(a)") " 6 Set the sphere center at midpoint between two CPs"            
            write(*,"(a)") " 10 Set the sphere radius"
            write(*,"(a)") " 11 Set the number of points in the sphere"
            read(*,*) isel2

            if (isel2==-9) then
                exit
            else if (isel2==0.or.isel2==-1.or.isel2==-2) then
                if (isel2==-2) then
                    write(*,*) "Input the indices of the atoms, e.g. 3,4,5,12, at most 1000 characters"
                    read(*,"(a)") c1000
                    call str2arr(c1000,nsearchcen,searchcenlist)
                end if
                if (isel2==0) nsearchcen=1
                if (isel2==-1) nsearchcen=ncenter
                !4.189=4/3*pi, this assess how many points need to be generated in the cube
                !so that numsearchpt points could in the sphere
                numcpold=numcp
                numsearchpt_tmp=nint(8D0/4.189D0*numsearchpt)
                allocate(randptx(numsearchpt_tmp),randpty(numsearchpt_tmp),randptz(numsearchpt_tmp))
                
                itime=0
                ioutcount=0
                do icenidx=1,nsearchcen
                    icen=icenidx !isel==0.or.isel==-1
                    if (isel2==-2) icen=searchcenlist(icenidx)
                    if (isel2==-1.or.isel2==-2) then !Cycle each atom center
                        sphcenx=a(icen)%x
                        sphceny=a(icen)%y
                        sphcenz=a(icen)%z
                    end if
                    CALL RANDOM_NUMBER(randptx)
                    CALL RANDOM_NUMBER(randpty)
                    CALL RANDOM_NUMBER(randptz)
                    randptx=randptx*2*toposphrad+(sphcenx-toposphrad) !Move distribution center of random point to sphere center
                    randpty=randpty*2*toposphrad+(sphceny-toposphrad)
                    randptz=randptz*2*toposphrad+(sphcenz-toposphrad)
                    randptx(1)=sphcenx !The first try point is set to sphere center, this is faciliate to locate CP at nuclei
                    randpty(1)=sphceny
                    randptz(1)=sphcenz
nthreads=getNThreads()
!$OMP PARALLEL DO SHARED(itime) PRIVATE(i) schedule(dynamic) NUM_THREADS(nthreads)
                    do i=1,numsearchpt_tmp
                        dispt_cen=dsqrt( (randptx(i)-sphcenx)**2+(randpty(i)-sphceny)**2+(randptz(i)-sphcenz)**2 )
                        if (dispt_cen>toposphrad) cycle
                        call findcp(randptx(i),randpty(i),randptz(i),ifunctopo,0)
!$OMP CRITICAL
                        itime=itime+1
                        if (itime>ioutcount+99.or.ifunctopo==12) then
                            write(*,"('#',i10,' /',i10)") itime,numsearchpt*nsearchcen
                            ioutcount=ioutcount+100
                        end if
!$OMP end CRITICAL
                    end do
!$OMP END PARALLEL DO
                end do
                deallocate(randptx,randpty,randptz)
                
                if ((numcp-numcpold)/=0) then
!                     call sortCP(numcpold+1) !Senseless here, because the guessing points occur randomly
                    write(*,*) "                            ==== Summary ===="
                    write(*,*) " Index              Coordinate               Type"
                    do icp=numcpold+1,numcp
                        write(*,"(i6,3f12.6,3x,a)") icp,CPpos(:,icp),CPtyp2lab(CPtype(icp))
                    end do
                end if
                write(*,"(' Totally find',i6,' new critical points')") numcp-numcpold
                write(*,*)
            else if (isel2==1) then
                write(*,*) "Input x,y,z   e.g.  1.2,0.2,-0.44"
                read(*,*) sphcenx,sphceny,sphcenz
            else if (isel2==2) then
                write(*,*) "Input atom index"
                read(*,*) iatm
                if (iatm>ncenter.or.iatm<=0) then
                    write(*,*) "Invalid input"
                else
                    sphcenx=a(iatm)%x
                    sphceny=a(iatm)%y
                    sphcenz=a(iatm)%z
                end if
            else if (isel2==3) then
                write(*,*) "Input index of the two atoms,  e.g.  3,7"
                read(*,*) iatm,jatm
                if ( iatm>ncenter.or.iatm<=0.or.jatm>ncenter.or.jatm<=0 ) then
                    write(*,*) "Invalid input"
                else
                    sphcenx=(a(iatm)%x+a(jatm)%x)/2D0
                    sphceny=(a(iatm)%y+a(jatm)%y)/2D0
                    sphcenz=(a(iatm)%z+a(jatm)%z)/2D0
                end if
            else if (isel2==4) then
                write(*,*) "Input index of the three atoms,  e.g.  2,3,7"
                read(*,*) iatm,jatm,katm
                if ( iatm>ncenter.or.iatm<=0.or.jatm>ncenter.or.jatm<=0.or.katm>ncenter.or.katm<=0 ) then
                    write(*,*) "Invalid input"
                else
                    sphcenx=(a(iatm)%x+a(jatm)%x+a(katm)%x)/3D0
                    sphceny=(a(iatm)%y+a(jatm)%y+a(katm)%y)/3D0
                    sphcenz=(a(iatm)%z+a(jatm)%z+a(katm)%z)/3D0
                end if
            else if (isel2==5) then
                write(*,*) "Input CP index"
                read(*,*) icp
                if (icp>numcp.or.icp<=0) then
                    write(*,*) "Invalid input"
                else
                    sphcenx=CPpos(1,icp)
                    sphceny=CPpos(2,icp)
                    sphcenz=CPpos(3,icp)
                end if
            else if (isel2==6) then
                write(*,*) "Input index of the two CPs,  e.g.  3,7"
                read(*,*) icp,jcp
                if ( icp>numcp.or.icp<=0.or.jcp>numcp.or.jcp<=0 ) then
                    write(*,*) "Invalid input"
                else
                    sphcenx=(CPpos(1,icp)+CPpos(1,jcp))/2D0
                    sphceny=(CPpos(2,icp)+CPpos(2,jcp))/2D0
                    sphcenz=(CPpos(3,icp)+CPpos(3,jcp))/2D0
                end if            
            else if (isel2==10) then
                write(*,*) "Input a radius"
                read(*,*) toposphrad
            else if (isel2==11) then
                write(*,*) "Input a number"
                read(*,*) numsearchpt
            end if
        end do
!7777777777777777777
    else if (isel==7) then
        write(*,*) "Input the index of the CP that you are interested in"
        write(*,"(a)") " Note: If input 0, then properties of all CPs will be outputted to CPprop.txt in current folder &
        (and if you feel the output speed is slow, you can input -1 to  avoid outputting ESP, which is the most expensive one)"
        read(*,*) indcp
        if (indcp==0.or.indcp==-1) then
            write(*,*) "Please wait..."
            iback=ishowptESP
            if (indcp==-1) ishowptESP=0 !Avoid outputting ESP
            open(10,file="CPprop.txt",status="replace")
            do icp=1,numcp
                write(*,"(' Outputting CP',i6,'  /',i6)") icp,numcp
                write(10,"(' ================   CP',i6,',     Type ',a,'   ================')") icp,CPtyp2lab(CPtype(icp))
                write(10,"(' Position (Bohr):',3f20.14)") CPpos(:,icp)
                call showptprop(CPpos(1,icp),CPpos(2,icp),CPpos(3,icp),ifunctopo,10)
                write(10,*)
            end do
            close(10)
            if (indcp==-1) ishowptESP=iback
            write(*,*) "Done! The results have been outputted to CPprop.txt in current folder"
            write(*,*) "Note: Unless otherwise specified, all units are in a.u."
        else if (indcp>0.and.indcp<=numcp) then
            write(*,*) "Note: Unless otherwise specified, all units are in a.u."
            write(*,"(' CP Position:',3f20.14)") CPpos(:,indcp)
            write(*,"(' CP type: ',a)") CPtyp2lab(CPtype(indcp))
            call showptprop(CPpos(1,indcp),CPpos(2,indcp),CPpos(3,indcp),ifunctopo,6)
        else
            write(*,*) "Error: Invalid input"
        end if
!8888888888888888888
    else if (isel==8) then
        numpathold=numpath
        do i=1,numcp
            if (CPtype(i)==2) call findpath(i,1,ifunctopo)
        end do
        write(*,"(' Totally found',i6,' new paths')") numpath-numpathold
!9999999999999999999
    else if (isel==9) then
        numpathold=numpath
        do i=1,numcp
            if (CPtype(i)==3) call findpath(i,2,ifunctopo)
        end do
        write(*,"(' Totally found',i6,' new paths')") numpath-numpathold
!10 10 10 10 10 10 10
    else if (isel==10.and.count(CPtype(1:numcp)==2)==0) then
        write(*,*) "Error: You have to find at least one (3,-1) critical point"
    else if (isel==10) then
        do while(.true.)
            write(*,*)
            write(*,*) "Generate or delete interbasin surface for which (3,-1)?"
            write(*,*) "e.g. 5 means generate the surface from the (3,-1) with index of 5"
            if (numbassurf> 0) write(*,*) "     -3 means delete the surface from the (3,-1) with index of 3" 
            if (numbassurf> 0) write(*,*) "If input 0, surfaces from all (3,-1) will be deleted"
            if (numbassurf==0) write(*,*) "If input 0, surfaces from all (3,-1) will be generated"
            if (numbassurf> 0) write(*,*) "To list all generated surfaces, input the letter ""l"""
            if (numbassurf> 0) write(*,*) "To export the surfaces from the (3,-1) with index of 4, input ""o 4"""
            write(*,*) "To return, input ""q"""
            read(*,"(a)") c200

            if (c200(1:1)=='q') then
                exit
            else if (c200(1:1)=='l') then
                do isurf=1,numbassurf
                    do icp=1,numcp
                        if (cp2surf(icp)==isurf) exit
                    end do
                    write(*,"('Index of surface:',i8,'     Index of corresponding (3,-1):',i8)") isurf,icp
                end do
            else if (c200(1:1)=='o') then
                read(c200(3:),*) icp
                if (icp>numcp.or.icp<=0) then
                    write(*,*) "The index of the surface is nonexisted"
                else if (cp2surf(icp)==0) then
                    write(*,*) "This CP is not (3,-1), input again"
                else
                    open(10,file="surpath.txt",status="replace")
                    do ipath=1,nsurfpathpercp
                        write(10,"('Path',i8)") ipath
                        do ipt=1,nsurfpt
                            write(10,"(i6,3f14.8)") ipt,bassurpath(:,ipt,ipath,cp2surf(icp))
                        end do
                    end do
                    close(10)
                    write(*,"(a)") "The coordinates of the paths of the surface have been exported to surpath.txt in current folder"
                end if
            else
                read(c200,*) isel2
                if (isel2==0) then
                    if (numbassurf>0) then
                        numbassurf=0
                        cp2surf=0 !If cps2surf(i)==0, means the (3,-1) with total index of i hasn't been given surface
                        deallocate(bassurpath)
                    else if (numbassurf==0) then !Generate all interbasin surfaces
                        numbassurf=count(CPtype(1:numcp)==2)
                        allocate(bassurpath(3,nsurfpt,nsurfpathpercp,numbassurf))
                        isurf=1
                        write(*,"(i8,' surfaces will be generated, please wait...')") numbassurf
                        do icp=1,numcp
                            if (CPtype(icp)==2) then
                                cp2surf(icp)=isurf
                                call genbassurf(icp,isurf,ifunctopo)
                                isurf=isurf+1
                                write(*,"(a,i6)") " Finished the surface generation from (3,-1) with index of",icp
                            end if
                        end do
                    end if
                else if (abs(isel2)>numcp) then
                    write(*,*) "Error: This CP is nonexisted"
                else if (CPtype(abs(isel2))/=2) then
                    write(*,*) "This CP is not (3,-1), input again"
                else if (isel2>0) then !Add a surface
                    if (cp2surf(isel2)/=0) then
                        write(*,*) "The interbasin surface has already been generated" 
                    else
                        write(*,*) "Please wait..."
                        if (numbassurf>0) then
                            allocate(bassurpathtmp(3,nsurfpt,nsurfpathpercp,numbassurf))
                            bassurpathtmp=bassurpath
                            deallocate(bassurpath)
                            numbassurf=numbassurf+1
                            allocate(bassurpath(3,nsurfpt,nsurfpathpercp,numbassurf))
                            bassurpath(:,:,:,1:numbassurf-1)=bassurpathtmp(:,:,:,:)
                            deallocate(bassurpathtmp)
                        else if (numbassurf==0) then
                            numbassurf=numbassurf+1
                            allocate(bassurpath(3,nsurfpt,nsurfpathpercp,numbassurf))
                        end if
                        call genbassurf(isel2,numbassurf,ifunctopo)
                        cp2surf(isel2)=numbassurf
                        write(*,*) "Done!"
                    end if
                else if (isel2<0) then !Delete a surface
                    idelsurf=cp2surf(abs(isel2))
                    if (idelsurf==0) then
                        write(*,*) "Surface has not been generated from this CP"
                    else
                        allocate(bassurpathtmp(3,nsurfpt,nsurfpathpercp,numbassurf))
                        bassurpathtmp=bassurpath
                        deallocate(bassurpath)
                        numbassurf=numbassurf-1
                        allocate(bassurpath(3,nsurfpt,nsurfpathpercp,numbassurf))
                        bassurpath(:,:,:,1:idelsurf-1)=bassurpathtmp(:,:,:,1:idelsurf-1)
                        bassurpath(:,:,:,idelsurf:)=bassurpathtmp(:,:,:,idelsurf+1:)
                        deallocate(bassurpathtmp)
                        cp2surf(abs(isel2))=0
                        where(cp2surf>idelsurf) cp2surf=cp2surf-1
                    end if
                end if
            end if
        end do
    else if (isel==20.and.ifunctopo==1) then
        do while(.true.)
            write(*,*) "Input the indices of the CPs in the ring, e.g. 22,23,25,28,32,11"
            write(*,*) "(Input q can exit)"
            read(*,"(a)") c200
            if (c200(1:1)=='q'.or.c200(1:1)=='Q') exit
            call str2arr(c200,nshanaromat,shanCPind)
            allocate(shanCPrho(nshanaromat))
            totdens=0D0
            do ishan=1,nshanaromat
                shanCPrho(ishan)=fdens(CPpos(1,shanCPind(ishan)),CPpos(2,shanCPind(ishan)),CPpos(3,shanCPind(ishan)))
            end do
            totrho=sum(shanCPrho(1:nshanaromat))
            shant=0D0
            do ishan=1,nshanaromat
                shanentropy=-shanCPrho(ishan)/totrho*log(shanCPrho(ishan)/totrho)
                write(*,"(' Electron density at CP',i3,':',f15.10,'  Local entropy:',f15.10)") shanCPind(ishan),shanCPrho(ishan),shanentropy
                shant=shant+shanentropy
            end do
            shanmax=log(dfloat(nshanaromat))
            write(*,"(' Total electron density:',f15.10)") totrho
            write(*,"(' Total Shannon entropy:',f15.10)") shant
            write(*,"(' Expected maximum Shannon entropy:',f15.10)") shanmax        
            write(*,*)
            write(*,"(' Shannon aromaticity index:',f16.10)") shanmax-shant 
            deallocate(shanCPrho)
        end do
    else if (isel==21.and.ifunctopo==1) then
        write(*,*) "Input the coordinate, e.g. 2.0,2.4,1.1     or input indices of a CP, e.g. 4"
        read(*,"(a)") c200
        if ( index(c200,',')==0 .and. index(trim(c200),' ')==0 ) then
            read(c200,*) ithisCP
            tmpx=CPpos(1,ithisCP)
            tmpy=CPpos(2,ithisCP)
            tmpz=CPpos(3,ithisCP)
        else
            read(c200,*) tmpx,tmpy,tmpz
            write(*,*) "The coordinate you inputted is in which unit?  1=Bohr  2=Angstrom"
            read(*,*) iunit
            if (iunit==2) then
                tmpx=tmpx/b2a
                tmpy=tmpy/b2a
                tmpz=tmpz/b2a
            end if
        end if
        write(*,*) "Input indices of three atoms to define a plane, e.g. 3,4,9"
        read(*,*) iatm1,iatm2,iatm3
        call pointABCD(a(iatm1)%x,a(iatm1)%y,a(iatm1)%z,a(iatm2)%x,a(iatm2)%y,a(iatm2)%z,a(iatm3)%x,a(iatm3)%y,a(iatm3)%z,xnor,ynor,znor,rnouse) !Normal vector is (xnor,ynor,znor)
        facnorm=sqrt(xnor**2+ynor**2+znor**2)
        xnor=xnor/facnorm !Normalize normal vector, then (xnor,ynor,znor) is the unit vector normal to the plane defined by iatm1,iatm2,iatm3
        ynor=ynor/facnorm
        znor=znor/facnorm
        if (allocated(b)) then
            call gencalchessmat(2,1,tmpx,tmpy,tmpz,densvalue,gradtmp,hesstmp)
            densgrad=xnor*gradtmp(1)+ynor*gradtmp(2)+znor*gradtmp(3)
            denscurvature=xnor*xnor*hesstmp(1,1)+xnor*ynor*hesstmp(1,2)+xnor*znor*hesstmp(1,3)+&
                          ynor*xnor*hesstmp(2,1)+ynor*ynor*hesstmp(2,2)+ynor*znor*hesstmp(2,3)+&
                          znor*xnor*hesstmp(3,1)+znor*ynor*hesstmp(3,2)+znor*znor*hesstmp(3,3)
            write(*,"(' The unit normal vector is',3f14.8)") xnor,ynor,znor
            write(*,"(' Electron density is          ',f30.10)") densvalue
            write(*,"(' Electron density gradient is ',f30.10)") densgrad
            write(*,"(' Electron density curvature is',f30.10)") denscurvature
        end if
        write(*,*)
        write(*,"(a)") " BTW: The X,Y,Z coordinate (row) of current point, the points below and above 1 Angstrom of the plane from current point, respectively (in Angstrom)."
        write(*,"(3f16.10)") tmpx*b2a,tmpy*b2a,tmpz*b2a
        write(*,"(3f16.10)") (tmpx-xnor/b2a)*b2a,(tmpy-ynor/b2a)*b2a,(tmpz-znor/b2a)*b2a
        write(*,"(3f16.10)") (tmpx+xnor/b2a)*b2a,(tmpy+ynor/b2a)*b2a,(tmpz+znor/b2a)*b2a
        write(*,*)
    end if
end do
end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!--------- Generate interbasin surface from (3,-1)
subroutine genbassurf(ithisCP,ithissurf,ifunc)
use defvar
use function
use util
use topo
implicit real*8(a-h,o-z)
integer ifunc,ithisCP,ithissurf
real*8 initstpsize,hess(3,3),grad(3),eigvecmat(3,3),eigval(3),basvec1(3),basvec2(3),k1(3),k2(3)
initstpsize=surfpathstpsiz/4D0 !smaller than pathstepsize(0.02)

call gencalchessmat(2,ifunc,CPpos(1,ithisCP),CPpos(2,ithisCP),CPpos(3,ithisCP),value,grad,hess)
call diagmat(hess,eigvecmat,eigval,300,1D-12)
!Generate normalized basis vectors from two negative eignvectors of hessian matrix
itmp=0
do i=1,3
    if (eigval(i)<0) then
        if (itmp==0) then
            rnorm=dsqrt(eigvecmat(1,i)**2+eigvecmat(2,i)**2+eigvecmat(3,i)**2)
            basvec1=eigvecmat(:,i)/rnorm
            itmp=1
        else if (itmp==1) then
            rnorm=dsqrt(eigvecmat(1,i)**2+eigvecmat(2,i)**2+eigvecmat(3,i)**2)
            basvec2=eigvecmat(:,i)/rnorm
        end if
    end if
end do
!Use the two basis vectors to generate a circle of initial points around (3,-1)
angstp=360D0/nsurfpathpercp
do i=1,nsurfpathpercp
    ang=(i-1)*angstp/180D0*pi !Convert to arc unit
    bassurpath(:,1,i,ithissurf)=initstpsize*( sin(ang)*basvec1(:)+cos(ang)*basvec2(:) ) + CPpos(:,ithisCP)
end do

!Generate gradient path from each initial point to comprise interbasin surface
do ipath=1,nsurfpathpercp
    do ipt=2,nsurfpt
        !Move point, RK2 method. Only calculate function value and gradient
        xtmp=bassurpath(1,ipt-1,ipath,ithissurf)
        ytmp=bassurpath(2,ipt-1,ipath,ithissurf)
        ztmp=bassurpath(3,ipt-1,ipath,ithissurf)
        call gencalchessmat(1,ifunc,xtmp,ytmp,ztmp,value,grad,hess)
        k1=grad/dsqrt(sum(grad**2))
        call gencalchessmat(1,ifunc,xtmp+surfpathstpsiz/2*k1(1),ytmp+surfpathstpsiz/2*k1(2),ztmp+surfpathstpsiz/2*k1(3),value,grad,hess)
        k2=grad/dsqrt(sum(grad**2))
        bassurpath(:,ipt,ipath,ithissurf)=bassurpath(:,ipt-1,ipath,ithissurf)-surfpathstpsiz*k2
    end do
end do
end subroutine



!!--------- Generate interbasin path from (3,-1) on a given plane
!ifunc=which real space function, ithisCP: The total index of (3,-1), ithispath: Store to which slot of ple3n1path array
subroutine gen3n1plepath(ifunc,ithisCP,ithispath)
use defvar
use topo
use function
use util
implicit real*8(a-h,o-z)
integer ifunc,ithisCP,ithispath
real*8 initstpsize,hess(3,3),grad(3),eigvecmat(3,3),eigval(3),initvec(3),plenormvec(3),k1(3),k2(3)
initstpsize=ple3n1pathstpsiz/4D0 !smaller than pathstepsize(0.02)

call gencalchessmat(2,ifunc,CPpos(1,ithisCP),CPpos(2,ithisCP),CPpos(3,ithisCP),value,grad,hess)
call diagmat(hess,eigvecmat,eigval,300,1D-12)
do iposvec=1,3
    if (eigval(iposvec)>0) exit
end do
if (plesel==1) then !XY
    plenormvec(1:2)=0D0
    plenormvec(3)=1D0
else if (plesel==2) then !XZ
    plenormvec(1:3)=0D0
    plenormvec(2)=1D0
else if (plesel==3) then !YZ
    plenormvec(2:3)=0D0
    plenormvec(1)=1D0
else
    call pointABCD(a1x,a1y,a1z,a2x,a2y,a2z,a3x,a3y,a3z,plenormvec(1),plenormvec(2),plenormvec(3),tmp)
end if

call vecprod(plenormvec(1),plenormvec(2),plenormvec(3), eigvecmat(1,iposvec),eigvecmat(2,iposvec),eigvecmat(3,iposvec), initvec(1),initvec(2),initvec(3))
rnorm=dsqrt(sum(initvec**2))
initvec=initvec/rnorm
ple3n1path(:,1,1,ithispath)=CPpos(:,ithisCP)+initvec(:)*initstpsize
ple3n1path(:,1,2,ithispath)=CPpos(:,ithisCP)-initvec(:)*initstpsize

!Generate gradient path from the given (3,-1) to comprise interbasin path on the plane
do idir=1,2
    do ipt=2,n3n1plept
        !Move point, RK2 method. Only calculate function value and gradient, don't calculate Hessian for saving time
        xtmp=ple3n1path(1,ipt-1,idir,ithispath)
        ytmp=ple3n1path(2,ipt-1,idir,ithispath)
        ztmp=ple3n1path(3,ipt-1,idir,ithispath)
        call gencalchessmat(1,ifunc,xtmp,ytmp,ztmp,value,grad,hess)
        k1=grad/dsqrt(sum(grad**2))
        call gencalchessmat(1,ifunc,xtmp+surfpathstpsiz/2*k1(1),ytmp+surfpathstpsiz/2*k1(2),ztmp+surfpathstpsiz/2*k1(3),value,grad,hess)
        k2=grad/dsqrt(sum(grad**2))
        ple3n1path(:,ipt,idir,ithispath)=ple3n1path(:,ipt-1,idir,ithispath)-surfpathstpsiz*k2
    end do
end do
end subroutine



!!!--------- Find path from a critical point at x,y,z
!itype=1: from (3,-1) to (3,-3)  =2: from (3,+1) to (3,+3)  =3: between (3,+1) and (3,-1)
subroutine findpath(ithisCP,itype,ifunc)
use topo
use function
use util
implicit real*8(a-h,o-z)
integer ifunc,ithisCP,foundind(npathtry) !foundind records that in pathtmp, which saved newly found path, =1 means saved
real*8 grad(3),hess(3,3),k1(3),k2(3)
real*8 eigvecmat(3,3),eigval(3),pathtmp(3,maxpathpt,npathtry) !Will maximally try to find npathtry paths
real*8,allocatable :: tmparr(:,:,:)
!Note: Each time invoke this routine, for itype=1 and 2, search two times; for itype=3, search npathtry times
!The new path first store to pathtmp, then enlarge topopath and pass the new path information to it

foundind=0
noldpath=numpath
ipath=0
!Determine eigenvalue and eigenvector of Hessian at initial point
call gencalchessmat(2,ifunc,CPpos(1,ithisCP),CPpos(2,ithisCP),CPpos(3,ithisCP),value,grad,hess)
call diagmat(hess,eigvecmat,eigval,300,1D-12)

if (itype==1.or.itype==2) then
    pathtmp(:,1,1)=CPpos(:,ithisCP) !Set first point as input coordinate
    pathtmp(:,1,2)=CPpos(:,ithisCP)
    if (itype==1) then
        do iposi=1,3
            if (eigval(iposi)>0) exit !Find positive eigenvalue of Hessian of (3,-1)
        end do
    else
        do iposi=1,3
            if (eigval(iposi)<0) exit !Find negative eigenvalue of Hessian of (3,+1)
        end do
    end if
iterdir:    do idir=1,2
        if (idir==1) write(*,"(' Go forward from CP: ',i6,1x,a,' Position:',3f12.6)") ithisCP,CPtyp2lab(CPtype(ithisCP)),CPpos(1:3,ithisCP)
        if (idir==2) write(*,"(' Go backward from CP:',i6,1x,a,' Position:',3f12.6)") ithisCP,CPtyp2lab(CPtype(ithisCP)),CPpos(1:3,ithisCP)
        posvecnorm=dsqrt(sum(eigvecmat(:,iposi)**2))
        if (idir==1) pathtmp(:,2,idir)=pathtmp(:,1,idir)+pathstepsize*eigvecmat(:,iposi)/posvecnorm !Move forwards along eigenvector with positive eigenvalue
        if (idir==2) pathtmp(:,2,idir)=pathtmp(:,1,idir)-pathstepsize*eigvecmat(:,iposi)/posvecnorm !Move backwards along eigenvector with positive eigenvalue
        !Check if the path has already presented by comparing corresponding first two points
        if (allocated(topopath)) then
            do ickpath=1,noldpath
!                  write(*,*) ickpath,numpath,size(topopath,3)
                if ( dsqrt(sum( (topopath(:,1,ickpath)-pathtmp(:,1,idir))**2 ))<0.01D0 &
                .and.dsqrt(sum( (topopath(:,2,ickpath)-pathtmp(:,2,idir))**2 ))<0.01D0) then
                    write(*,*) "The path has already presented, skip search"
                    write(*,*)
                    cycle iterdir
                end if
            end do
        end if

iterpt:    do ipt=2,maxpathpt
            !Check if the distance between current point and any other CP is smaller than threshold (discritpathfin)
            do icp=1,numcp
                if (icp==ithisCP) cycle
                distcp=dsqrt(sum( (pathtmp(:,ipt,idir)-CPpos(:,icp))**2 ))
                if (distcp<discritpathfin) then
                    numpath=numpath+1
                    foundind(idir)=1
                    pathnumpt(numpath)=ipt
                    write(*,"(a,i6,a,f8.4)") " Found new path after",ipt," iterations, path length:",(ipt-1)*pathstepsize
                    write(*,"(' Reached CP',i6,1x,a,' Position:',3f12.6)") icp,CPtyp2lab(CPtype(icp)),CPpos(:,icp)
!                     do i=1,ipt
!                         write(*,"(i6,3f12.6)") i,pathtmp(:,i,idir)
!                     end do
                    exit iterpt
                end if
            end do
            if (ipt==maxpathpt) then
                write(*,*) "Search failed, exceeded upper limit of number of iterations"
                exit
            end if
            
            !Move point, RK2 method. Only calculate function value and gradient
            xtmp=pathtmp(1,ipt,idir)
            ytmp=pathtmp(2,ipt,idir)
            ztmp=pathtmp(3,ipt,idir)
            call gencalchessmat(1,ifunc,xtmp,ytmp,ztmp,value,grad,hess)
            k1=grad/dsqrt(sum(grad**2))
            call gencalchessmat(1,ifunc,xtmp+pathstepsize/2*k1(1),ytmp+pathstepsize/2*k1(2),ztmp+pathstepsize/2*k1(3),value,grad,hess)
            k2=grad/dsqrt(sum(grad**2))
            if (itype==1) pathtmp(:,ipt+1,idir)=pathtmp(:,ipt,idir)+pathstepsize*k2
            if (itype==2) pathtmp(:,ipt+1,idir)=pathtmp(:,ipt,idir)-pathstepsize*k2
        end do iterpt

        write(*,*)
    end do iterdir
else if (itype==3) then
end if

!Enlarge old topopath and then add the path found at this time to it
numnewpath=numpath-noldpath !How many new path found in this time
if (numnewpath/=0) then
    itmp=1
    if (allocated(topopath)) then
        allocate(tmparr(3,maxpathpt,noldpath))
        tmparr=topopath
        deallocate(topopath)
        allocate(topopath(3,maxpathpt,numpath))
        topopath(:,:,1:noldpath)=tmparr
        do i=1,npathtry !The number of foundind(i)==1 equals numnewpath
            if (foundind(i)==1) then
                topopath(:,:,noldpath+itmp)=pathtmp(:,:,i)
                itmp=itmp+1
            end if
        end do
        deallocate(tmparr)
    else !This first time invoke this routine
        allocate(topopath(3,maxpathpt,numnewpath))
        do i=1,npathtry
            if (foundind(i)==1) then
                topopath(:,:,itmp)=pathtmp(:,:,i)
                itmp=itmp+1
            end if
        end do
    end if
end if
end subroutine



!!--------- Determine path type and connected which two CPs, 0 means unknown (not connected a found CP)
!ipathtype =0: other   =1: (3,-1)->(3,-3) =2: (3,+1)->(3,+3) =3: (3,-1)<-->(3,+1)
subroutine path_cp(ipath,icp1,icp2,ipathtype)
use topo
implicit real*8(a-h,o-z)
integer ipath,icp1,icp2
iunknown=0
do icp1=1,numcp
    if (sum( (topopath(:,1,ipath)-CPpos(:,icp1))**2 )<discritpathfin**2) exit !Test the first point in the path
    if (icp1==numcp) iunknown=1
end do
if (iunknown==1) icp1=0
iunknown=0
do icp2=1,numcp
    if (sum( (topopath(:,pathnumpt(ipath),ipath)-CPpos(:,icp2))**2 )<discritpathfin**2) exit !Test the last point
    if (icp2==numcp) iunknown=0
end do
if (iunknown==1) icp2=0

if (icp1==0.or.icp2==0) then
    ipathtype=0
else
    icp1type=CPtype(icp1)
    icp2type=CPtype(icp2)
    if (icp1type==2.and.icp2type==1) then
        ipathtype=1
    else if (icp1type==3.and.icp2type==4) then
        ipathtype=2
    else if ( (icp1type==2.and.icp2type==3).or.(icp1type==3.and.icp2type==2) ) then
        ipathtype=3
    else
        ipathtype=0
    end if
end if
end subroutine



!!!----------- Find critical points from initial guess at X,Y,Z
!ifunc is index of real space functions
!If ilowcrit==1, use lower critiera, because for heavy atoms, cusp of electron density at nuclear position is very sharp hence hard to locate by default criteria
!ishowsearchlevel=0/1/2:  Print none/some detail/all detail. Notice that in parallel mode, the outputted details are messed up
subroutine findcp(x,y,z,ifunc,ilowcrit)
use topo
use function
use util
implicit real*8(a-h,o-z)
integer ifunc
integer ilowcrit
real*8 x,y,z
real*8 coord(3,1),grad(3,1),hess(3,3),disp(3,1)
real*8 eigvecmat(3,3),eigval(3) !,tmpmat(3,3)
if (ilowcrit==1) then
    realdispconv=dispconv*10000D0
    realgradconv=gradconv*10000D0
else !For most cases use default criteria
    realdispconv=dispconv
    realgradconv=gradconv
end if
coord(1,1)=x
coord(2,1)=y
coord(3,1)=z
if (ishowsearchlevel>=1) write(*,"(' Starting point:',3f12.6)") coord(1:3,1)

do i=1,topomaxcyc
    call gencalchessmat(2,ifunc,coord(1,1),coord(2,1),coord(3,1),value,grad(1:3,1),hess)
!     call showmatgau(hess)
!     write(*,*) detmat(hess)
    singulartest=abs(detmat(hess))
    if (singulartest<singularcrit) then
        if (ishowsearchlevel>=1) then
            write(*,*) "Hessian matrix is singular at current position, stop iteration"
            write(*,"(' Absolute of determinant of Hessian matrix:',E18.10)") singulartest
            write(*,"(' Criterion for detecting singular:',E18.10)") singularcrit
        end if
        exit
    end if
    disp=-matmul(invmat(hess,3),grad)
    coord=coord+CPstepscale*disp
    disperr=dsqrt(sum(disp**2))
    graderr=dsqrt(sum(grad**2))

    if (ishowsearchlevel==2) then
        write(*,"(/,' Step',i5,'  Function Value:',f18.10)") i,value
        write(*,"(' Displacement vector:',3f18.10)") disp
        write(*,"(' Current coordinate :',3f18.10)") coord
        write(*,"(' Current gradient   :',3E18.10)") grad
        write(*,"(' Norm of displacement:',E18.8,'  Norm of gradient:',E18.8)") disperr,graderr
        write(*,"(' Goal: |disp|<',E18.8,'  |Grad|<',E18.8)") realdispconv,realgradconv
    end if
!     tmpmat=hess
!     call diagmat(tmpmat,eigvecmat,eigval,300,1D-12)
!     write(*,"('Eigenvalue of Hessian   :  ',3D16.8)") eigval

    if (disperr<realdispconv.and.graderr<realgradconv) then
        if (ishowsearchlevel>=1) write(*,"(' After',i6,' iterations')") i
        if (ishowsearchlevel==2) write(*,*) "====================== Iteration ended ======================"
        inewcp=1
        do icp=1,numcp
            r=dsqrt(sum( (coord(:,1)-CPpos(:,icp))**2 ))
            if (r<=minicpdis) then
                if (ishowsearchlevel>=1) write(*,"(a,i6,a)") " This CP is too close to CP",icp,", ignored..."
                inewcp=0
                exit
            end if
        end do
        if (CPsearchlow/=CPsearchhigh.and.(value<CPsearchlow.or.value>CPsearchhigh)) then
            write(*,"(a,1PD12.5,a)") " The value of this CP is ",value,", which exceeded user-defined range and thus ignored"
            inewcp=0
        end if
        if (inewcp==1) then
!$OMP CRITICAL
            numcp=numcp+1
            CPpos(:,numcp)=coord(:,1)
            call diagmat(hess,eigvecmat,eigval,300,1D-15)
!             call diagsymat(hess,eigvecmat,eigval,idiagok)
!             if (idiagok/=0) write(*,*) "Diagonization of Hessian matrix failed"
            igt0=count(eigval>0)
            if (igt0==3) then
                if (ishowsearchlevel>=1) write(*,"(' Found new (3,+3) at',3f15.10)") coord
                CPtype(numcp)=4
            else if (igt0==2) then
                if (ishowsearchlevel>=1) write(*,"(' Found new (3,+1) at',3f15.10)") coord
                CPtype(numcp)=3
            else if (igt0==1) then
                if (ishowsearchlevel>=1) write(*,"(' Found new (3,-1) at',3f15.10)") coord
                CPtype(numcp)=2
                call sort(eigval)
                if (ishowsearchlevel>=1) write(*,"(' Bond ellipticity is',f15.10)") eigval(1)/eigval(2)-1.0D0
            else if (igt0==0) then
                if (ishowsearchlevel>=1) write(*,"(' Found new (3,-3) at',3f15.10)") coord
                CPtype(numcp)=1
            end if
            if (ishowsearchlevel>=1) write(*,"(' Eigenvalue:',3f20.10)") eigval
!$OMP end CRITICAL
        end if
        exit
    end if
    if (i==topomaxcyc.and.(ishowsearchlevel>=1)) write(*,*) "!! Exceeded maximum cycle until find stationary point !!"
end do
if (ishowsearchlevel>=1) write(*,*)
end subroutine


!Sort newly found CPs according to coordinates
!numcpoldp1: The number of CPs before this search + 1
subroutine sortCP(numcpoldp1)
use topo
implicit real*8(a-h,o-z)
integer numcpoldp1,typetmp
real*8 tmparr(3)
do itmp=numcpoldp1,numcp
    do jtmp=itmp+1,numcp
        tmpvali=CPpos(1,itmp)*0.234134D0+CPpos(2,itmp)*1.9837322D0-CPpos(3,itmp)*0.5413578924D0 !Use three random number to generate unique code
        tmpvalj=CPpos(1,jtmp)*0.234134D0+CPpos(2,jtmp)*1.9837322D0-CPpos(3,jtmp)*0.5413578924D0 
        if (tmpvali>tmpvalj) then
            typetmp=CPtype(itmp)
            CPtype(itmp)=CPtype(jtmp)
            CPtype(jtmp)=typetmp
            tmparr=CPpos(:,itmp)
            CPpos(:,itmp)=CPpos(:,jtmp)
            CPpos(:,jtmp)=tmparr
        end if
    end do
end do
end subroutine
