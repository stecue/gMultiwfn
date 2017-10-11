!!------ Forming basin for specific real space function, and integrate real space functions in the basins
!!------ The method is adapted from J. Phys.: Condens. Matter 21 (2009) 084204
!Grid must be ortho
subroutine basinana
use defvar
use basinintmod
use util
implicit real*8(a-h,o-z)
integer walltime1,walltime2
integer :: igridmethod=3
integer,allocatable :: attconvold(:),realattconv(:),usercluslist(:)
character basinfilename*200,selectyn,c80tmp*80,ctmp1*20,ctmp2*20,ctmp3*20,ctmp4*20,c1000tmp*1000,c200tmp*200
real*8 :: threslowvalatt=1D-5
ishowattlab=1
ishowatmlab=0
ishowatt=1
!Don't show topology analysis results to avoid confusion
ishowCPlab=0
ishowpathlab=1
ishow3n3=0
ishow3n1=0
ishow3p1=0
ishow3p3=0
!When enter this module, we reset the whole state, which is signified by numatt. Because the user may calculate grid data at external functions,
!and hence messed up grid setting (orgx/y/z,dx/y/z...), so all of the functions in this module that rely on grid setting will be totally wrong
numatt=0

call delvirorb(1)

do while(.true.)
    write(*,*)
    write(*,*) "                 ============= Basin analysis ============="
    write(*,*) "-10 Return to main menu"
    if (numatt>0) write(*,*) "-6 Set parameter for attractor clustering or manually perform clustering"
    if (numatt==0) write(*,*) "-6 Set parameter for attractor clustering"
    if (numatt>0) write(*,*) "-5 Export basins as cube file"
    if (numatt>0) write(*,*) "-4 Export attractors as pdb file"
    if (numatt>0) write(*,*) "-3 Show information of attractors"
    if (numatt>0) write(*,*) "-2 Measure distances, angles and dihedral angles between attractors or atoms"
    if (igridmethod==1) write(*,*) "-1 Select the method for generating basins, current: On-grid"
    if (igridmethod==2) write(*,*) "-1 Select the method for generating basins, current: Near-grid"
    if (igridmethod==3) write(*,*) "-1 Select the method for generating basins, current: Near-grid with refinement"
    if (numatt>0) write(*,*) " 0 Visualize attractors and basins"
    if (numatt>0) then
        write(*,*) " 1 Regenerate basins and relocate attractors"
    else
        write(*,*) " 1 Generate basins and locate attractors"
    end if
    if (numatt>0) then
        write(*,*) " 2 Integrate real space functions in the basins"
        write(*,*) " 3 Calculate electric multipole moments in the basins"
        write(*,*) " 4 Calculate localization index and delocalization index for the basins"
        write(*,*) " 5 Output orbital overlap matrix in basins to BOM.txt in current folder"
!         write(*,*) " 6 Integrate real space functions in the basins with multi-level refinement"
        if (ifuncbasin==1) then
            write(*,*) " 7 Integrate real space functions in AIM basins with mixed type of grids"
            write(*,*) " 8 Calculate electric multipole moments in AIM basins with mixed type of grids"
            write(*,*) " 9 Obtain atomic contribution to population of external basins"
        end if
    end if
    read(*,*) isel
    if (numatt==0.and.(isel==-5.or.isel==-4.or.isel==-3.or.isel==-2.or.isel==0.or.isel==2.or.isel==3.or.isel==4.or.isel==5)) then
        write(*,*) "Error: You should use option 1 to generate basins first!"
        cycle
    end if

    if (isel==-10) then
        idrawbasinidx=-10 !Don't display interbasin in any other GUI
        ishowattlab=0 !Don't show attractors in other GUIs
        ishowatt=1
        return
        
    else if (isel==-6) then
        do while(.true.)
            write(*,*) "0 Return"
            write(*,"(a,f14.5,a)") " 1 Set relative value difference criterion, current:",valcritclus*100,"%"
            write(*,"(a,i3)") " 2 Set the multiplier for distance criterion, current:",mergeattdist
            if (numatt>=2) write(*,*) "3 Cluster specified attractors"
            read(*,*) isel2
            if (isel2==0) then
                exit
            else if (isel2==1) then
                write(*,"(a)") " Input a value in percentage, e.g. 0.5"
                write(*,"(a)") " Note: Inputting X means if the value difference between two attractors is less than X%, &
                then they will be regarded as degenerate and may be clustered according to distance criterion"
                read(*,*) valcritclus
                valcritclus=valcritclus/100D0
            else if (isel2==2) then
                write(*,"(a)") " Input a value, e.g. 6"
                write(*,"(a)") " Note: Inputting P means for any two attractors that have relative value difference less than the one set by option 1, &
                if their interval is less than P*sqrt(dx^2+dy^2+dz^2), &
                where dx,dy,dz are grid spacing in X,Y,Z,then they will be clustered together. If you want to nullify the clustering step after generating &
                basins, simply set this value to 0"
                read(*,*) mergeattdist
            else if (isel2==3) then
                if (numatt<2) then
                    write(*,*) "Error: At least two attractors must be presented!"
                    cycle
                end if
                write(*,*) "Input the attractors you want to cluster, e.g. 4,5,8,9,22,21"
                read(*,"(a)") c1000tmp
                call str2arr(c1000tmp,ntobeclus) !Find how many terms
                if (allocated(attconvold)) deallocate(attconvold,realattconv,usercluslist)
                allocate(attconvold(-2:numatt),realattconv(-2:numrealatt),usercluslist(ntobeclus))
                realattconv(-2)=-2
                realattconv(-1)=-1
                realattconv(0)=0
                call str2arr(c1000tmp,ntobeclus,usercluslist)
                call sorti4(usercluslist)
                attconvold=attconv
                idesrealatt=usercluslist(1)
                do itmp=2,ntobeclus
                    irealatt=usercluslist(itmp)
                    do jtmp=1,nrealatthas(irealatt)
                        iatt=realatttable(irealatt,jtmp)
                        attconv(iatt)=idesrealatt
                    end do
                end do
                call clusdegenatt(1)
                realattconv(1)=1
                do itmp=1,numatt !Build conversion relationship between previous real attractors and new real attractors
                    if (itmp/=1) then
                        if (attconvold(itmp)/=attconvold(itmp-1)) realattconv(attconvold(itmp))=attconv(itmp)
                    end if
                end do
                do iz=2,nz-1
                    do iy=2,ny-1
                        do ix=2,nx-1
                            gridbas(ix,iy,iz)=realattconv(gridbas(ix,iy,iz))
                        end do
                    end do
                end do
                call detectinterbasgrd(6) !Detect interbasin grids
                write(*,*) "Done!"
                write(*,*)
            end if
        end do
    
    else if (isel==-5) then !Outout cube file for basins for visualizing interbasin surface
        write(*,*) "Input the index range of basin"
        write(*,"(a)") " Inputting 4,6 will output basin 4,5,6 to basin0004.cub, basin0005.cub, basin0006.cub to current folder. Inputting 5,5 will only output basin 5"
        write(*,"(a)") " Inputting 0,0, will output all basins into basins.cub, the grid value corresponds to basin index"
        read(*,*) ilow,ihigh
        if (ilow/=0) then
            write(*,*) "If output internal region of the basin? (y/n)"
            read(*,*) selectyn
        end if
        if (allocated(cubmattmp)) deallocate(cubmattmp)
        allocate(cubmattmp(nx,ny,nz))
        do ibas=ilow,ihigh
            if (ilow/=0) then
                write(basinfilename,"('basin',i4.4,'.cub')") ibas
                write(*,"(' Outputting basin',i8,' as ',a)") ibas,trim(basinfilename)
                cubmattmp=gridbas(:,:,:)
                where(cubmattmp/=ibas)
                    cubmattmp=0
                elsewhere
                    cubmattmp=1
                end where
                if (selectyn=='n') where (interbasgrid .eqv. .false.) cubmattmp=0
            else
                write(*,*) "Outputting basin.cub..."
                write(basinfilename,"('basin.cub')")
                cubmattmp=gridbas(:,:,:)
            end if
            open(10,file=basinfilename,status="replace")
            write(10,"(' Generated by Multiwfn')")
            write(10,"(' Totally ',i12,' grid points')") nx*ny*nz
            write(10,"(i5,3f12.6)") ncenter,orgx,orgy,orgz
            write(10,"(i5,3f12.6)") nx,dx,0.0,0.0
            write(10,"(i5,3f12.6)") ny,0.0,dy,0.0
            write(10,"(i5,3f12.6)") nz,0.0,0.0,dz
            do icenter=1,ncenter
                write(10,"(i5,4f12.6)") a(icenter)%index,a(icenter)%charge,a(icenter)%x,a(icenter)%y,a(icenter)%z
            end do
            do ix=1,nx
                do iy=1,ny
                    write(10,"(6(1PE13.5))",advance="no") cubmattmp(ix,iy,1:nz)
                    write(10,*)
                end do
            end do
            close(10)
        end do
        deallocate(cubmattmp)
        if (ilow/=0) then
            write(*,"(a)") " Done! Cube files for the selected basins have been outputted to current folder"
            write(*,"(a)") " The value 1 and 0 in the files denote the corresponding point belongs and not belongs to the basin, & 
            respectively. You can plot such as the isosurface with isovalue=0.5 to visualize the basins"
        else
            write(*,"(a)") " Done! basin.cub has been outputted to current folder. The grid values correspond to basin index"
        end if
        
    else if (isel==-4) then
        open(10,file="attractors.pdb",status="replace")
        do iatt=1,numatt
            write(10,"(a6,i5,1x,a4,1x,a3, 1x,a1,i4,4x,3f8.3,2f6.2,10x,a2)") "HETATM",iatt,' '//"C "//' ',"ATT",'A',attconv(iatt),attxyz(iatt,:)*b2a,1.0,0.0,"C "
        end do
        close(10)
        write(*,*) "Done, all attractors have been exported to attractors.pdb in current folder"
        write(*,*) "Note: The residue indices in the pdb file correspond to attractor indices"
        
    else if (isel==-3) then
        write(*,*) "  Attractor       X,Y,Z coordinate (Angstrom)            Value"
        do irealatt=1,numrealatt
            write(*,"(i8,3f14.8,1E18.9)") irealatt,realattxyz(irealatt,1:3)*b2a,realattval(irealatt)
        end do
        do irealatt=1,numrealatt
            if (nrealatthas(irealatt)/=1) then
                write(*,*)
                write(*,"(' The members of degenerate attractor',i6,':')") irealatt
                do itmp=1,nrealatthas(irealatt)
                    iatt=realatttable(irealatt,itmp)
                    write(*,"(i8,'  XYZ:',3f13.7,'  Value:',1E18.9)") iatt,attxyz(iatt,:)*b2a,attval(iatt)
                end do
            end if
        end do
        
    else if (isel==-2) then
        write(*,*) "q = Quit. Selection method:"
        write(*,*) "a? = Atom ?"
        write(*,*) "c? = Attractor ? (If the attractor is degenerate, average coordinate is used)"
        write(*,*) "d? = Attractor ? (Only for degenerate attractors, cycle each of its members)"
        write(*,*) "For example:"
        write(*,*) """a1 c3"" returns the distance between atom1 and att.3"
        write(*,*) """a4 a2"" returns the distance between atom4 and atom2"
        write(*,*) """c6 a2 a5"" returns the angle of att.6-atom2-atom5"
        write(*,*) """c2 c4 a3 c7"" returns the dihedral angle of att.2-att.4-atom3-att.7"
        write(*,*) """d3 a1"" returns the distance between members of att.3 and atom1"
        write(*,"(a)") " ""d3 a1 c2"" returns the the angle of (members of att.3)--atom1--att.2, &
        meanwhile outputs the vertical distance from (members of att.3) to the line linking atom1 and att.2"
        write(*,*)
        do while(.true.)
            read(*,"(a)") c80tmp
            c80tmp=adjustl(c80tmp)
            imeasure=0
            do ichar=1,len_trim(c80tmp) !imeasure=1/2/3: measure distance,angle,dihedral
                if (c80tmp(ichar:ichar)==','.or.c80tmp(ichar:ichar)==' ') imeasure=imeasure+1
            end do
            nelement=0 !Validate "d" selection
            do ichar=1,len_trim(c80tmp)
                if (c80tmp(ichar:ichar)=='d') nelement=nelement+1
            end do
            if (nelement/=0) then
                if (c80tmp(1:1)/='d'.or.nelement>1) then
                    write(*,*) "Error: ""d"" type of selection can only occur once and must be the first term"
                    cycle
                else if (imeasure==3) then
                    write(*,*) "Error: Dihedral angle calculation doesn't support ""d"" type of selection"
                    cycle                    
                end if
            end if
                            
            if (c80tmp(1:1)=='q') then
                exit
            else if (imeasure==1.or.imeasure==2.or.imeasure==3) then
                if (imeasure==1) read(c80tmp,*) ctmp1,ctmp2 !Read two terms
                if (imeasure==2) read(c80tmp,*) ctmp1,ctmp2,ctmp3 !Read three terms
                if (imeasure==3) read(c80tmp,*) ctmp1,ctmp2,ctmp3,ctmp4 !Read four terms
                
                if (ctmp1(1:1)=='a') then
                    read(ctmp1(2:),*) iatm
                    tmpx1=a(iatm)%x
                    tmpy1=a(iatm)%y
                    tmpz1=a(iatm)%z
                else if (ctmp1(1:1)=='c') then
                    read(ctmp1(2:),*) irealatt
                    tmpx1=realattxyz(irealatt,1)
                    tmpy1=realattxyz(irealatt,2)
                    tmpz1=realattxyz(irealatt,3)
                else if (ctmp1(1:1)=='d') then
                    read(ctmp1(2:),*) irealattde
                end if
                if (ctmp2(1:1)=='a') then
                    read(ctmp2(2:),*) iatm
                    tmpx2=a(iatm)%x
                    tmpy2=a(iatm)%y
                    tmpz2=a(iatm)%z
                else if (ctmp2(1:1)=='c') then
                    read(ctmp2(2:),*) irealatt
                    tmpx2=realattxyz(irealatt,1)
                    tmpy2=realattxyz(irealatt,2)
                    tmpz2=realattxyz(irealatt,3)
                end if
                
                if (imeasure==1) then
                    if (ctmp1(1:1)/='d') then
                        write(*,"(' The distance is',f12.6,' Bohr (',f12.6 ' Angstrom)')") &
                        dsqrt((tmpx1-tmpx2)**2+(tmpy1-tmpy2)**2+(tmpz1-tmpz2)**2),dsqrt((tmpx1-tmpx2)**2+(tmpy1-tmpy2)**2+(tmpz1-tmpz2)**2)*b2a
                    else if (ctmp1(1:1)=='d') then
                        distmin=9999999
                        distmax=-1
                        distavg=0
                        do itmp=1,nrealatthas(irealattde)
                            iatt=realatttable(irealattde,itmp)
                            tmpx1=attxyz(iatt,1)
                            tmpy1=attxyz(iatt,2)
                            tmpz1=attxyz(iatt,3)
                            distnow=dsqrt((tmpx1-tmpx2)**2+(tmpy1-tmpy2)**2+(tmpz1-tmpz2)**2)
                            distavg=distavg+distnow
                            if (distnow<distmin) distmin=distnow
                            if (distnow>distmax) distmax=distnow
                        end do
                        distavg=distavg/nrealatthas(irealattde)
                        write(*,"(' Note: Attractor',i6,' has',i6,' members')") irealattde,nrealatthas(irealattde)
                        write(*,"(' The minimum distance is',f12.6,' Bohr (',f12.6 ' Angstrom)')") distmin,distmin*b2a
                        write(*,"(' The maximum distance is',f12.6,' Bohr (',f12.6 ' Angstrom)')") distmax,distmax*b2a
                        write(*,"(' The average distance is',f12.6,' Bohr (',f12.6 ' Angstrom)')") distavg,distavg*b2a
                    end if
                end if

                if (imeasure==2.or.imeasure==3) then !Analyze one more term, then print angle
                    if (ctmp3(1:1)=='a') then
                        read(ctmp3(2:),*) iatm
                        tmpx3=a(iatm)%x
                        tmpy3=a(iatm)%y
                        tmpz3=a(iatm)%z
                    else if (ctmp3(1:1)=='c') then
                        read(ctmp3(2:),*) irealatt
                        tmpx3=realattxyz(irealatt,1)
                        tmpy3=realattxyz(irealatt,2)
                        tmpz3=realattxyz(irealatt,3)
                    end if
                end if
                if (imeasure==2) then
                    if (ctmp1(1:1)/='d') then
                        write(*,"(' The angle is',f12.6,' degree')") xyz2angle(tmpx1,tmpy1,tmpz1,tmpx2,tmpy2,tmpz2,tmpx3,tmpy3,tmpz3)
                    else if (ctmp1(1:1)=='d') then
                        angmin=180
                        angmax=-1
                        angavg=0
                        distmin=999999
                        distmax=-1
                        distavg=0
                        prjx=0
                        prjy=0
                        prjz=0
                        do itmp=1,nrealatthas(irealattde)
                            iatt=realatttable(irealattde,itmp)
                            tmpx1=attxyz(iatt,1)
                            tmpy1=attxyz(iatt,2)
                            tmpz1=attxyz(iatt,3)
                            angnow=xyz2angle(tmpx1,tmpy1,tmpz1,tmpx2,tmpy2,tmpz2,tmpx3,tmpy3,tmpz3)
                            angavg=angavg+angnow
                            if (angnow<angmin) angmin=angnow
                            if (angnow>angmax) angmax=angnow
                            distnow=potlinedis(tmpx1,tmpy1,tmpz1,tmpx2,tmpy2,tmpz2,tmpx3,tmpy3,tmpz3)
                            distavg=distavg+distnow
                            if (distnow<distmin) distmin=distnow
                            if (distnow>distmax) distmax=distnow
                            call pointprjline(tmpx1,tmpy1,tmpz1,tmpx2,tmpy2,tmpz2,tmpx3,tmpy3,tmpz3,prjxtmp,prjytmp,prjztmp)
                            prjx=prjx+prjxtmp
                            prjy=prjy+prjytmp
                            prjz=prjz+prjztmp
                        end do
                        angavg=angavg/nrealatthas(irealattde)
                        distavg=distavg/nrealatthas(irealattde)
                        prjx=prjx/nrealatthas(irealattde)
                        prjy=prjy/nrealatthas(irealattde)
                        prjz=prjz/nrealatthas(irealattde)
                        write(*,"(' Note: Attractor',i6,' has',i6,' members')") irealattde,nrealatthas(irealattde)
                        write(*,"(' The minimum angle is',f12.6,' degree')") angmin
                        write(*,"(' The maximum angle is',f12.6,' degree')") angmax
                        write(*,"(' The average angle is',f12.6,' degree')") angavg
                        write(*,"(' The minimum distance is',f12.6,' Bohr (',f12.6 ' Angstrom)')") distmin,distmin*b2a
                        write(*,"(' The maximum distance is',f12.6,' Bohr (',f12.6 ' Angstrom)')") distmax,distmax*b2a
                        write(*,"(' The average distance is',f12.6,' Bohr (',f12.6 ' Angstrom)')") distavg,distavg*b2a
                        write(*,"(' The average X,Y,Z coordinate by projecting the members of ',a,' to the line linking ',a,' and ',a)") trim(ctmp1),trim(ctmp2),trim(ctmp3)
                        write(*,"(3f14.8,'   Angstrom')") prjx*b2a,prjy*b2a,prjz*b2a
                    end if
                end if
                
                if (imeasure==3) then !Analyze one more term, then print dihedral angle
                    if (ctmp4(1:1)=='a') then
                        read(ctmp4(2:),*) iatm
                        tmpx4=a(iatm)%x
                        tmpy4=a(iatm)%y
                        tmpz4=a(iatm)%z
                    else if (ctmp4(1:1)=='c') then
                        read(ctmp4(2:),*) irealatt
                        tmpx4=realattxyz(irealatt,1)
                        tmpy4=realattxyz(irealatt,2)
                        tmpz4=realattxyz(irealatt,3)
                    end if
                    write(*,"(' The dihedral angle is',f12.6,' degree')") xyz2dih(tmpx1,tmpy1,tmpz1,tmpx2,tmpy2,tmpz2,tmpx3,tmpy3,tmpz3,tmpx4,tmpy4,tmpz4)
                end if
            else
                write(*,*) "Error: Invalid input"
            end if
        end do
        
    else if (isel==-1) then
        write(*,*) "1: On-grid method, Comput .Mat. Sci., 36, 354 (2006)"
        write(*,*) "2: Near-grid method, J. Phys.: Condens. Matter, 21, 08420 (2009)"
        write(*,*) "3: Near-grid method with boundary refinement step"
        write(*,"(a)") " Note: Near-grid method (adapted by Tian Lu) is more accurate than On-grid method and thus is more recommended; with the boundary refinement step, the result will be better"
        read(*,*) igridmethod
        
    else if (isel==0) then
        ioldtextheigh=textheigh
        textheigh=40 !Default textheigh is too small
        textheigh=ioldtextheigh
        
    else if (isel==1) then
        isourcedata=1 !Default status is calculating grid data here
        if (allocated(cubmat)) then
            write(*,"(a)") " Note: There has been a grid data in the memory, please select generating the basins by which manner"
            write(*,*) "0 Return"
            write(*,*) "1 Generate the basins by selecting a real space function"
            write(*,*) "2 Generate the basins by using the grid data stored in memory"
            read(*,*) isourcedata
        end if
        if (isourcedata==0) then
            cycle
        else if (isourcedata==1) then
            write(*,*) "Select the real space function to be integrated"
            call selfunc_interface(ifuncbasin)
            call setgridforbasin(ifuncbasin)
            if (allocated(cubmat)) deallocate(cubmat)
            allocate(cubmat(nx,ny,nz))
            call savecubmat(ifuncbasin,0,iorbsel)
        else if (isourcedata==2) then
            ifuncbasin=1 !I assume that the file loaded records electron density
            write(*,"('The range of x is from ',f12.6,' to ',f12.6,' Bohr,' i5,' points')") ,orgx,orgx+(nx-1)*dx,nx
            write(*,"('The range of y is from ',f12.6,' to ',f12.6,' Bohr,',i5,' points')") ,orgy,orgy+(ny-1)*dy,ny
            write(*,"('The range of z is from ',f12.6,' to ',f12.6,' Bohr,',i5,' points')") ,orgz,orgz+(nz-1)*dz,nz
            write(*,"('Total number of grid points is ',i10)") nx*ny*nz
            write(*,"('Grid spacing in X,Y,Z (Bohr):',3f12.6)") dx,dy,dz
        end if
        
        allocate(grdposneg(nx,ny,nz)) !Record which grids have negative value
        grdposneg=.true.
        do iz=1,nz
            do iy=1,ny
                do ix=1,nx
                    if (cubmat(ix,iy,iz)<0D0) then
                        grdposneg(ix,iy,iz)=.false.
                        cubmat(ix,iy,iz)=-cubmat(ix,iy,iz)
                    end if
                end do
            end do
        end do
!         where (cubmat<0D0)  !DO NOT USE THIS, because I found that "where" will consuming vary large amount of memory!
!             grdposneg=.false.
!             cubmat=-cubmat !Invert negative values to positive, after basins are generated the values will be recovered
!         end where
        if (allocated(gridbas)) deallocate(gridbas)
        allocate(gridbas(nx,ny,nz))
        call setupmovevec
        write(*,*)
        write(*,*) "Generating basins, please wait..."
        call walltime(walltime1)
        CALL CPU_TIME(time_begin)
        call generatebasin(igridmethod)
        CALL CPU_TIME(time_end)
        call walltime(walltime2)
        write(*,"(' Generating basins took up CPU time',f12.2,'s, wall clock time',i10,'s')") time_end-time_begin,walltime2-walltime1    
        numunassign=count(gridbas(2:nx-1,2:ny-1,2:nz-1)==0)
        write(*,"(' The number of unassigned grids:',i12)") numunassign
        numgotobound=count(gridbas(2:nx-1,2:ny-1,2:nz-1)==-1)
        write(*,"(' The number of grids travelled to box boundary:',i12)") numgotobound
        where (grdposneg .eqv. .false.) cubmat=-cubmat !Recover original grid data
        deallocate(grdposneg)
        
        do iatt=1,numatt !Eliminate the attractors with very low value
            if ( abs(cubmat(attgrid(iatt,1),attgrid(iatt,2),attgrid(iatt,3)))<threslowvalatt ) then
                write(*,"(' Note: There are attractors having very low absolute value (<',1PE8.2,') and thus insignificant, how to deal with them?')") threslowvalatt
                write(*,*) "1 Do nothing"
                write(*,*) "2 Set corresponding grids as unassigned status"
                write(*,*) "3 Assign corresponding grids to the nearest significant attractors"
                write(*,*) "Hint: For most cases, option 3 is recommended"
                read(*,*) isel2
                call elimlowvalatt(threslowvalatt,isel2)
                exit
            end if
        end do
        
        !Generate actual coordinate of attractors, which will be used to plot in the local GUI, and in any other external GUIs
        if (allocated(attxyz)) deallocate(attxyz,attval)
        allocate(attxyz(numatt,3),attval(numatt))
        do iatt=1,numatt
            ix=attgrid(iatt,1)
            iy=attgrid(iatt,2)
            iz=attgrid(iatt,3)
            attxyz(iatt,1)=orgx+(ix-1)*dx
            attxyz(iatt,2)=orgy+(iy-1)*dy
            attxyz(iatt,3)=orgz+(iz-1)*dz
            attval(iatt)=cubmat(ix,iy,iz)
        end do
        !Currently attractors have been finally determined, one shouldn't perturb them further more
        
        call clusdegenatt(0) !Cluster degenerate attractors as "real attractors" and calculate average coordinate and value for the real attractors

        call detectinterbasgrd(6) !Detect interbasin grids
        numinterbas=count(interbasgrid .eqv. .true.)
        write(*,"(' The number of interbasin grids:',i12)") numinterbas    
        
    else if (isel==2) then
        call integratebasin
        
    else if (isel==3.or.isel==4.or.isel==5.or.isel==7.or.isel==-7.or.isel==8) then
        if (.not.allocated(b)) then
            write(*,"(a)") " Note: No GTF (gauss type function) information is available in your input file, please input the file &
            containing GTF information of your system, such as .wfn/.wfx and .fch file. e.g. c:\abc.wfn"
            read(*,"(a)") c200tmp
            call readinfile(c200tmp,1)
        end if
        if (isel==3) then
            call multipolebasin
        else if (isel==4) then
            call LIDIbasin(0)
        else if (isel==5) then
            call LIDIbasin(1)
        else if (isel==7) then
            write(*,*) "0 Return"
            write(*,*) "1 Integrate a specific function with atomic-center + uniform grids"
            write(*,*) "2 The same as 1, but with exact refinement of basin boundary"
            write(*,*) "3 The same as 2, but with approximate refinement of basin boundary"
            write(*,*) "Hint:"
            write(*,*) "Accuracy: 2>=3>>1     Time spent: 2>3>>1     Memory requirement: 3>2=1"
    !         write(*,*) "Robost: 1=2>3"
            read(*,*) iseltmp
            call integratebasinmix(iseltmp)
        else if (isel==-7) then
            call integratebasinmix_LSB
        else if (isel==8) then
            call integratebasinmix(10)
        end if
        
!     else if (isel==6) then
!         call integratebasinrefine
    else if (isel==9) then
        call atmpopinbasin
    end if
end do

end subroutine



!!---- setup move vector and move length for the 26 neighbour
subroutine setupmovevec
use defvar
use basinintmod
vec26x=0
vec26y=0
vec26z=0
!The nearest neighbours:
len26(1:2)=dx
vec26x(1)=1
vec26x(2)=-1
len26(3:4)=dy
vec26y(3)=1
vec26y(4)=-1
len26(5:6)=dz
vec26z(5)=1
vec26z(6)=-1
!On the edges:
len26(7:10)=dsqrt(dx*dx+dy*dy)
vec26x(7)=1
vec26y(7)=1
vec26x(8)=-1
vec26y(8)=1
vec26x(9)=-1
vec26y(9)=-1
vec26x(10)=1
vec26y(10)=-1
len26(11:14)=dsqrt(dx*dx+dz*dz)
vec26x(11)=1
vec26z(11)=1
vec26x(12)=-1
vec26z(12)=1
vec26x(13)=-1
vec26z(13)=-1
vec26x(14)=1
vec26z(14)=-1
len26(15:18)=dsqrt(dy*dy+dz*dz)
vec26y(15)=1
vec26z(15)=1
vec26y(16)=1
vec26z(16)=-1
vec26y(17)=-1
vec26z(17)=-1
vec26y(18)=-1
vec26z(18)=1
!At the vertices:
len26(19:26)=dsqrt(dx*dx+dy*dy+dz*dz)
vec26z(19:22)=1
vec26x(19)=1
vec26y(19)=1
vec26x(20)=-1
vec26y(20)=1
vec26x(21)=-1
vec26y(21)=-1
vec26x(22)=1
vec26y(22)=-1
vec26z(23:26)=-1
vec26x(23)=1
vec26y(23)=-1
vec26x(24)=1
vec26y(24)=1
vec26x(25)=-1
vec26y(25)=1
vec26x(26)=-1
vec26y(26)=-1
gridbas=0 !Unassigned state
gridbas(1,:,:)=-2 !Set status for box boundary grids
gridbas(nx,:,:)=-2
gridbas(:,1,:)=-2
gridbas(:,ny,:)=-2
gridbas(:,:,1)=-2
gridbas(:,:,nz)=-2
end subroutine



!!------- Generate basins from regular grid
! igridmethod=1: On grid, Comput.Mat.Sci.,36,354    =2: Near-grid method (slower, but more accurate), see J.Phys.:Condens.Matter,21,084204
! =3: Near-grid with refinement
! The near-grid method is improved by Tian Lu, namely at later stage automatically switch to on-grid method to guarantee convergence
subroutine generatebasin(igridmethod)
use defvar
use basinintmod
implicit real*8(a-h,o-z)
integer,parameter :: nmaxtrjgrid=3000
integer igridmethod
integer ntrjgrid !Recording trjgrid contains how many elements now
integer trjgrid(nmaxtrjgrid,3) !The trajectory contains which grids, record their indices sequentially, trjgrid(i,1/2/3) = ix/iy/iz of the ith grid
! real*8 trjval(nmaxtrjgrid),gradmaxval(nmaxtrjgrid) !******For debugging******
if (allocated(attgrid)) deallocate(attgrid)
allocate(attgrid(nint(nx*ny*nz/20D0),3)) !I think the number of attractors in general is impossible to exceeds nx*ny*nz/20
numatt=0
write(*,*) "  Attractor       X,Y,Z coordinate (Angstrom)                Value"
!Cycle all grids, but the box boundary ones are ignored (due to gridbas=-2), since can't evalute its gradients in all directions by finite difference

nsteplimit=min( nmaxtrjgrid,nint(dsqrt(dfloat(nx*nx+ny*ny+nz*nz))*2) )
nstepdiscorr=nint(nsteplimit/2D0)
if (igridmethod==1) then
1    continue
nthreads=getNThreads()
!$OMP PARALLEL DO private(ix,iy,iz,ntrjgrid,inowx,inowy,inowz,trjgrid,valnow,imove,gradtmp,igradmax,gradmax,iatt,itrjgrid,idtmp) &
!$OMP shared(gridbas,numatt,attgrid) schedule(DYNAMIC) NUM_THREADS(nthreads)
    do iz=2,nz-1
        do iy=2,ny-1
            do ix=2,nx-1
                if (gridbas(ix,iy,iz)/=0) cycle
!                 if (interbasgrid(ix,iy,iz) .eqv. .false.) cycle
                ntrjgrid=0
                inowx=ix
                inowy=iy
                inowz=iz
                do while(.true.) !Steepest ascent
                    ntrjgrid=ntrjgrid+1
                    if (ntrjgrid>nsteplimit) exit !Unconverged. The gridbas for these unassigned grids will still be 0
                    trjgrid(ntrjgrid,1)=inowx
                    trjgrid(ntrjgrid,2)=inowy
                    trjgrid(ntrjgrid,3)=inowz
                    !Test all 26 directions to find out the maximal gradient direction
                    valnow=cubmat(inowx,inowy,inowz)
                    do imove=1,26
                        gradtmp=(cubmat(inowx+vec26x(imove),inowy+vec26y(imove),inowz+vec26z(imove))-valnow)/len26(imove)
                        if (imove==1.or.gradtmp>gradmax) then
                            igradmax=imove
                            gradmax=gradtmp
                        end if
                    end do
                    !Test if this is an attractor, if yes, assign all grid in this trajectory
                    if (gradmax<=0) then !Equal sign is important, because when system has symmetry, adjacent grid may be degenerate about mirrow plane, now the ascent should be terminated
                        if (valnow==0D0) exit !The region far beyond system, the value may be exactly zero due to cutoff of exponent, these grids shouldn't be regarded as attractors
!$OMP CRITICAL
cyciatt3:                do iatt=1,numatt
                            if (inowx==attgrid(iatt,1).and.inowy==attgrid(iatt,2).and.inowz==attgrid(iatt,3)) then
                                do itrjgrid=1,ntrjgrid
                                    gridbas(trjgrid(itrjgrid,1),trjgrid(itrjgrid,2),trjgrid(itrjgrid,3))=iatt
                                end do
                                exit cyciatt3
                            end if
                            do imove=1,26
                                if (inowx+vec26x(imove)==attgrid(iatt,1).and.inowy+vec26y(imove)==attgrid(iatt,2).and.inowz+vec26z(imove)==attgrid(iatt,3)) then
                                    do itrjgrid=1,ntrjgrid
                                        gridbas(trjgrid(itrjgrid,1),trjgrid(itrjgrid,2),trjgrid(itrjgrid,3))=iatt
                                    end do
                                    exit cyciatt3
                                end if
                            end do
                        end do cyciatt3
                        if (iatt>numatt) then !A new attractor, iatt=numatt+1 currently
                            numatt=numatt+1
                            do itrjgrid=1,ntrjgrid
                                gridbas(trjgrid(itrjgrid,1),trjgrid(itrjgrid,2),trjgrid(itrjgrid,3))=numatt
                            end do
                            attgrid(numatt,1)=inowx
                            attgrid(numatt,2)=inowy
                            attgrid(numatt,3)=inowz
                            if (grdposneg(inowx,inowy,inowz) .eqv. .true.) then
                                write(*,"(i8,3f14.8,f20.8)") numatt,(orgx+(inowx-1)*dx)*b2a,(orgy+(inowy-1)*dy)*b2a,(orgz+(inowz-1)*dz)*b2a,cubmat(inowx,inowy,inowz)
                            else !This grid should has negative value
                                write(*,"(i8,3f14.8,f20.8)") numatt,(orgx+(inowx-1)*dx)*b2a,(orgy+(inowy-1)*dy)*b2a,(orgz+(inowz-1)*dz)*b2a,-cubmat(inowx,inowy,inowz)
                            end if
                        end if
!$OMP end CRITICAL
                        exit
                    end if
                    !Move to next grid
                    inowx=inowx+vec26x(igradmax)
                    inowy=inowy+vec26y(igradmax)
                    inowz=inowz+vec26z(igradmax)
                    !Test if this grid has already been assigned
                    idtmp=gridbas(inowx,inowy,inowz)
                    if (idtmp>0) then
                        do itrjgrid=1,ntrjgrid
                            gridbas(trjgrid(itrjgrid,1),trjgrid(itrjgrid,2),trjgrid(itrjgrid,3))=idtmp
                        end do
                        exit
                    end if
                    !Test if encountered box boundary
                    if (inowx==1.or.inowx==nx.or.inowy==1.or.inowy==ny.or.inowz==1.or.inowz==nz) then
                        do itrjgrid=1,ntrjgrid
                            gridbas(trjgrid(itrjgrid,1),trjgrid(itrjgrid,2),trjgrid(itrjgrid,3))=-1
                        end do
                        exit
                    end if
                end do !End ascent
            end do !End cycle x
        end do !End cycle y
    end do !End cycle z
!$OMP END PARALLEL DO
    
else if (igridmethod==2.or.igridmethod==3) then
nthreads=getNThreads()
!$OMP PARALLEL DO private(ix,iy,iz,corrx,corry,corrz,ntrjgrid,inowx,inowy,inowz,trjgrid,valnow,imove,gradtmp,igradmax,gradmax,iatt,&
!$OMP itrjgrid,idtmp,icorrx,icorry,icorrz,gradx,grady,gradz,sclgrad,ineiidx) shared(gridbas,numatt,attgrid) schedule(DYNAMIC) NUM_THREADS(nthreads)
    do iz=2,nz-1
        do iy=2,ny-1
            do ix=2,nx-1
                if (gridbas(ix,iy,iz)/=0) cycle
                ntrjgrid=0
                inowx=ix
                inowy=iy
                inowz=iz
                corrx=0D0 !Correction vector
                corry=0D0
                corrz=0D0
                do while(.true.) !Steepest ascent
                    ntrjgrid=ntrjgrid+1
                    trjgrid(ntrjgrid,1)=inowx
                    trjgrid(ntrjgrid,2)=inowy
                    trjgrid(ntrjgrid,3)=inowz
                    if (ntrjgrid>nsteplimit) then !These gridbas for these unconverged grids will still be 0
!                         do itmp=1,ntrjgrid-1 !******For debugging******
!                             write(*,"(i5,3i6,2f16.10)") itmp,trjgrid(itmp,1:3),trjval(itmp),gradmaxval(itmp)
!                         end do
                        exit
                    end if
                    !Test all 26 directions to find out the maximal gradient direction
                    valnow=cubmat(inowx,inowy,inowz)
                    do imove=1,26
                        gradtmp=(cubmat(inowx+vec26x(imove),inowy+vec26y(imove),inowz+vec26z(imove))-valnow)/len26(imove)
                        if (imove==1.or.gradtmp>gradmax) then
                            igradmax=imove
                            gradmax=gradtmp
                        end if
                    end do
!                      write(*,"(3i5,2f14.8)") inowx,inowy,inowz,gradmax,valnow !Trace the trajectory !******For debugging******
                    if (gradmax<=0D0) then !Equal sign is important, because when system has symmetry, adjacent grid may be degenerate about mirrow plane, now the ascent should be terminated
                        if (valnow==0D0) exit !The region far beyond system, the value may be exactly zero due to cutoff of exponent, these grids shouldn't be regarded as attractors
!$OMP CRITICAL
cyciatt:                do iatt=1,numatt
                            !Test if current grid is attractor iatt, is yes, assign all grid in this trajectory
                            if (inowx==attgrid(iatt,1).and.inowy==attgrid(iatt,2).and.inowz==attgrid(iatt,3)) then
                                do itrjgrid=1,ntrjgrid
                                    gridbas(trjgrid(itrjgrid,1),trjgrid(itrjgrid,2),trjgrid(itrjgrid,3))=iatt
                                end do
                                exit cyciatt
                            end if
                            !Test if neighbour grid (+/-x,+/-y,+/-z) is attractor iatt, is yes, assign all grid in this trajectory
                            !The reason I do this is because when the points are symmetric to Cartesian plane, many adjacent and degenerate attractors will occur,&
                            !while this treatment can combine them as a single one
                            do imove=1,26
                                if (inowx+vec26x(imove)==attgrid(iatt,1).and.inowy+vec26y(imove)==attgrid(iatt,2).and.inowz+vec26z(imove)==attgrid(iatt,3)) then
                                    do itrjgrid=1,ntrjgrid
                                        gridbas(trjgrid(itrjgrid,1),trjgrid(itrjgrid,2),trjgrid(itrjgrid,3))=iatt
                                    end do
                                    exit cyciatt
                                end if
                            end do
                        end do cyciatt
                        !A new attractor
                        if (iatt==numatt+1) then
                            numatt=numatt+1
                            do itrjgrid=1,ntrjgrid
                                gridbas(trjgrid(itrjgrid,1),trjgrid(itrjgrid,2),trjgrid(itrjgrid,3))=numatt
                            end do
                            attgrid(numatt,1)=inowx
                            attgrid(numatt,2)=inowy
                            attgrid(numatt,3)=inowz
                            if (grdposneg(inowx,inowy,inowz) .eqv. .true.) then
                                write(*,"(i8,3f14.8,f20.8)") numatt,(orgx+(inowx-1)*dx)*b2a,(orgy+(inowy-1)*dy)*b2a,(orgz+(inowz-1)*dz)*b2a,cubmat(inowx,inowy,inowz)
                            else !This grid should has negative value
                                write(*,"(i8,3f14.8,f20.8)") numatt,(orgx+(inowx-1)*dx)*b2a,(orgy+(inowy-1)*dy)*b2a,(orgz+(inowz-1)*dz)*b2a,-cubmat(inowx,inowy,inowz)
                            end if
                        end if
!$OMP end CRITICAL
                        exit
                    end if
                    !Correction step may lead to oscillator when encountering circular ELF/LOL attractor, so if the current number of step is already large
                    !(larger than half of upper limit of step number), then correction will be disabled (namely switch to on-grid method) to guarantee convergence.
                    if ( ntrjgrid<nstepdiscorr .and. (abs(corrx)>(dx/2D0).or.abs(corry)>(dy/2D0).or.abs(corrz)>(dz/2D0)) ) then !This time we do correction step
                        if (abs(corrx)>(dx/2D0)) then
                            icorrx=nint(corrx/abs(corrx)) !Get sign of corrx
                            inowx=inowx+icorrx
                            corrx=corrx-icorrx*dx
                        end if
                        if (abs(corry)>(dy/2D0)) then
                            icorry=nint(corry/abs(corry))
                            inowy=inowy+icorry
                            corry=corry-icorry*dy
                        end if
                        if (abs(corrz)>(dz/2D0)) then
                            icorrz=nint(corrz/abs(corrz))
                            inowz=inowz+icorrz
                            corrz=corrz-icorrz*dz
                        end if
                    else !Move to next grid according to maximal gradient and then update correction vector
                        !Calculate true gradient
                        gradx=(cubmat(inowx+1,inowy,inowz)-cubmat(inowx-1,inowy,inowz))/(2*dx)
                        grady=(cubmat(inowx,inowy+1,inowz)-cubmat(inowx,inowy-1,inowz))/(2*dy)
                        gradz=(cubmat(inowx,inowy,inowz+1)-cubmat(inowx,inowy,inowz-1))/(2*dz)
                        sclgrad=min(dx/abs(gradx),dy/abs(grady),dz/abs(gradz))
                        inowx=inowx+vec26x(igradmax)
                        inowy=inowy+vec26y(igradmax)
                        inowz=inowz+vec26z(igradmax)
                        corrx=corrx+gradx*sclgrad-vec26x(igradmax)*dx
                        corry=corry+grady*sclgrad-vec26y(igradmax)*dy
                        corrz=corrz+gradz*sclgrad-vec26z(igradmax)*dz
!                         if (abs(corrx)>(dx/2D0).or.abs(corry)>(dy/2D0).or.abs(corrz)>(dz/2D0)) cycle
                    end if
                    !Test if this grid has already been assigned and all of the neighbours were also assigned to the same attractor
                    idtmp=gridbas(inowx,inowy,inowz)
                    if (idtmp>0) then
                        do imove=1,26
                            ineiidx=gridbas(inowx+vec26x(imove),inowy+vec26y(imove),inowz+vec26z(imove))
                            if (ineiidx/=idtmp) exit
                        end do
                        if (imove==27) then
                            do itrjgrid=1,ntrjgrid
                                gridbas(trjgrid(itrjgrid,1),trjgrid(itrjgrid,2),trjgrid(itrjgrid,3))=idtmp
                            end do
                            exit
                        end if
                    end if
                    !Test if encountered box boundary                
                    if (inowx==1.or.inowx==nx.or.inowy==1.or.inowy==ny.or.inowz==1.or.inowz==nz) then
                        do itrjgrid=1,ntrjgrid
                            gridbas(trjgrid(itrjgrid,1),trjgrid(itrjgrid,2),trjgrid(itrjgrid,3))=-1
                        end do
                        exit
                    end if
                end do !End ascent
            end do !End cycle x
        end do !End cycle y
    end do !End cycle z
!$OMP END PARALLEL DO
    
    !!!!   Refining the basin boundary
    if (igridmethod==3) then
!         do itime=1,3 !Refine one time in general is sufficient
        write(*,*) "Detecting boundary grids..."
        call detectinterbasgrd(6) !It seems that using 26 directions to determine boundary grids doesn't bring evident benefit
        write(*,"(' There are',i12,' grids at basin boundary')") count(interbasgrid .eqv. .true.)
        write(*,*) "Refining basin boundary..."
        !Below code is the adapted copy of above near-grid code
nthreads=getNThreads()
!$OMP PARALLEL DO private(ix,iy,iz,corrx,corry,corrz,ntrjgrid,inowx,inowy,inowz,valnow,imove,gradtmp,igradmax,gradmax,iatt,&
!$OMP ineiidx,idtmp,icorrx,icorry,icorrz,gradx,grady,gradz,sclgrad) shared(gridbas,attgrid) schedule(DYNAMIC) NUM_THREADS(nthreads)
        do iz=2,nz-1
            do iy=2,ny-1
                do ix=2,nx-1
                    if (interbasgrid(ix,iy,iz) .eqv. .false.) cycle
                    if (gridbas(ix,iy,iz)<=0) cycle !Ignored the ones unassigned or gone to box boundary
                    ntrjgrid=0
                    inowx=ix
                    inowy=iy
                    inowz=iz
                    corrx=0D0 !Correction vector
                    corry=0D0
                    corrz=0D0
                    do while(.true.) !Steepest ascent
                        ntrjgrid=ntrjgrid+1
                        if (ntrjgrid>nsteplimit) exit
                        valnow=cubmat(inowx,inowy,inowz)
                        do imove=1,26
                            gradtmp=(cubmat(inowx+vec26x(imove),inowy+vec26y(imove),inowz+vec26z(imove))-valnow)/len26(imove)
                            if (imove==1.or.gradtmp>gradmax) then
                                igradmax=imove
                                gradmax=gradtmp
                            end if
                        end do
                        if (gradmax<=0) then !Equal sign is important, because when system has symmetry, adjacent grid may be degenerate about mirrow plane, now the ascent should be terminated
cyciatt2:                    do iatt=1,numatt
                                if (inowx==attgrid(iatt,1).and.inowy==attgrid(iatt,2).and.inowz==attgrid(iatt,3)) then
                                    gridbas(ix,iy,iz)=iatt
                                    exit cyciatt2
                                end if
                                do imove=1,26 !Test if neighbour grid (+/-x,+/-y,+/-z) is attractor iatt
                                    if (inowx+vec26x(imove)==attgrid(iatt,1).and.inowy+vec26y(imove)==attgrid(iatt,2).and.inowz+vec26z(imove)==attgrid(iatt,3)) then
                                        gridbas(ix,iy,iz)=iatt
                                        exit cyciatt2
                                    end if
                                end do
                            end do cyciatt2
                            if (iatt>numatt) then
                                write(*,*) "Warning: Found new attractor at refining process!"
                            end if
                            exit
                        end if
                        if ( ntrjgrid<nstepdiscorr .and. (abs(corrx)>(dx/2D0).or.abs(corry)>(dy/2D0).or.abs(corrz)>(dz/2D0)) ) then !This time we do correction step
                            if (abs(corrx)>(dx/2D0)) then
                                icorrx=nint(corrx/abs(corrx)) !Get sign of corrx
                                inowx=inowx+icorrx
                                corrx=corrx-icorrx*dx
                            end if
                            if (abs(corry)>(dy/2D0)) then
                                icorry=nint(corry/abs(corry))
                                inowy=inowy+icorry
                                corry=corry-icorry*dy
                            end if
                            if (abs(corrz)>(dz/2D0)) then
                                icorrz=nint(corrz/abs(corrz))
                                inowz=inowz+icorrz
                                corrz=corrz-icorrz*dz
                            end if
                        else !Move to next grid according to maximal gradient and then update correction vector
                            gradx=(cubmat(inowx+1,inowy,inowz)-cubmat(inowx-1,inowy,inowz))/(2*dx)
                            grady=(cubmat(inowx,inowy+1,inowz)-cubmat(inowx,inowy-1,inowz))/(2*dy)
                            gradz=(cubmat(inowx,inowy,inowz+1)-cubmat(inowx,inowy,inowz-1))/(2*dz)
                            sclgrad=min(dx/abs(gradx),dy/abs(grady),dz/abs(gradz))
                            inowx=inowx+vec26x(igradmax)
                            inowy=inowy+vec26y(igradmax)
                            inowz=inowz+vec26z(igradmax)
                            corrx=corrx+gradx*sclgrad-vec26x(igradmax)*dx
                            corry=corry+grady*sclgrad-vec26y(igradmax)*dy
                            corrz=corrz+gradz*sclgrad-vec26z(igradmax)*dz
                        end if
                        idtmp=gridbas(inowx,inowy,inowz)
                        if (ntrjgrid>60.and.idtmp>0) then !If enable this doesn't affect result detectably, but enabling it will evidently reduce computational cost
                            do imove=1,26
                                ineiidx=gridbas(inowx+vec26x(imove),inowy+vec26y(imove),inowz+vec26z(imove))
                                if (ineiidx/=idtmp) exit
                            end do
                            if (imove==27) then
                                gridbas(ix,iy,iz)=idtmp
                                exit
                            end if
                        end if
                        if (inowx==1.or.inowx==nx.or.inowy==1.or.inowy==ny.or.inowz==1.or.inowz==nz) exit
                    end do !End ascent
                end do !End cycle x
            end do !End cycle y
        end do !End cycle z
!$OMP END PARALLEL DO
!         end do
    end if
    
end if
end subroutine



!!------- Eliminate low value attractors and the corresponding basin
! imode=1 Do nothing
! imode=2 Set corresponding grids as unassigned status
! imode=3 Assign corresponding grids to the nearest significant attractors
subroutine elimlowvalatt(threslowvalatt,imode)
use defvar
use basinintmod
implicit real*8(a-h,o-z)
integer imode
real*8 threslowvalatt
integer highvalatt(numatt)
integer att2att(-2:numatt) !Convert indices of old attractors to new ones
integer atttypelist(numatt)
if (imode==1) then
    return
else
    icounthigh=0
    atttypelist=0
    do iatt=1,numatt !Generate list
        if (abs(cubmat(attgrid(iatt,1),attgrid(iatt,2),attgrid(iatt,3)))>=threslowvalatt) then
            icounthigh=icounthigh+1
            highvalatt(icounthigh)=iatt
            atttypelist(iatt)=1 !high
        end if
    end do
    nlowvalatt=numatt-icounthigh
    if (imode==2) then
        do iz=2,nz-1
            do iy=2,ny-1
                do ix=2,nx-1
                    if (atttypelist(gridbas(ix,iy,iz))==0) gridbas(ix,iy,iz)=0
                end do
            end do
        end do
    else if (imode==3) then
        do iz=2,nz-1
            rnowz=orgz+(iz-1)*dz
            do iy=2,ny-1
                rnowy=orgy+(iy-1)*dy
                do ix=2,nx-1
                    rnowx=orgx+(ix-1)*dx
                    if (atttypelist(gridbas(ix,iy,iz))==0) then !Find out which grid belongs to insignificant attractors
                        shortmaxsqr=1D20
                        do iatt=1,numatt !Will be attribute to the nearest significant attractors
                            if (atttypelist(iatt)==1) then
                                xatt=orgx+(attgrid(iatt,1)-1)*dx
                                yatt=orgy+(attgrid(iatt,2)-1)*dy
                                zatt=orgz+(attgrid(iatt,3)-1)*dz
                                disttmpsqr=(xatt-rnowx)**2+(yatt-rnowy)**2+(zatt-rnowz)**2
                                if (disttmpsqr<shortmaxsqr) then
                                    ishortmax=iatt
                                    shortmaxsqr=disttmpsqr
                                end if
                            end if
                        end do
                        gridbas(ix,iy,iz)=ishortmax
                    end if
                end do
            end do
        end do
    end if
    !Build conversion table between old and new attractors, so that the indices of attractors could be contiguous
    att2att=-3
    att2att(-2)=-2
    att2att(-1)=-1
    att2att(0)=0
    icount=0
    do iatt=1,numatt
        if (atttypelist(iatt)==0) cycle
        icount=icount+1
        att2att(iatt)=icount
    end do
    numatt=numatt-nlowvalatt
    !Update attgrid
    do iatt=1,numatt
        attgrid(iatt,:)=attgrid(highvalatt(iatt),:)
    end do
    !Final update of grid attribution
    do iz=2,nz-1
        do iy=2,ny-1
            do ix=2,nx-1
                gridbas(ix,iy,iz)=att2att(gridbas(ix,iy,iz))
            end do
        end do
    end do
    write(*,"(i6,' insignificant attractors have been eliminated')") nlowvalatt
end if
end subroutine



!!------- !Detect interbasin grids for real attractors
!ndir can be 6 and 26, meaning if only primary 6 directions or all 26 directions are used to determine boundary grids
subroutine detectinterbasgrd(ndir)
use defvar
use basinintmod
implicit real*8(a-h,o-z)
integer ndir
if (allocated(interbasgrid)) deallocate(interbasgrid)
allocate(interbasgrid(nx,ny,nz))
interbasgrid=.false.
do iz=2,nz-1
    do iy=2,ny-1
        do ix=2,nx-1
            idxnow=gridbas(ix,iy,iz)
            do imove=1,ndir !To reduce the number of interbasin grids and thus speed up displaying, I don't use full 26 direction as criterion
                idxmove=gridbas(ix+vec26x(imove),iy+vec26y(imove),iz+vec26z(imove))
                if ((idxmove/=idxnow).and.idxmove/=-2) then !This is a grid between two or more basins, and is not adjacent to box boundary
                    interbasgrid(ix,iy,iz)=.true.
                    exit
                end if
            end do
        end do
    end do
end do
end subroutine



!!------- Cluster degenerate and adjacent attractors together
!!Nearest neighbour method is used to cluster the attractors, see Leach book p493. Each cluster corresponds to a final attractor
!If imode=0, run the whole subroutine; if imode=1, then external attconv is provided, and only some code in this routine is used to make the index contiguous
subroutine clusdegenatt(imode)
use defvar
use basinintmod
implicit real*8(a-h,o-z)
integer imode
integer neleatt,eleatt(numatt) !eleatt(i)=j means the ith element in clustering stage corresponds to actual jth attractor. nclustlist is its total current elements
integer passlist(numatt) !if the ith element is 1, means i attractor has passed clustering stage
integer ncluster,ncluele(numatt),cluele(numatt,numatt) !cluele(i,4/5/9...) means cluster i has element 4,5,9..., ncluele(i) is current size of cluster i, ncluster is total number of cluster
real*8 attdismat(numatt,numatt),attvalmat(numatt,numatt) !Distance matrix and relative value difference matrix of the original attractors
distcrit=mergeattdist*dsqrt(dx*dx+dy*dy+dz*dz) !Distance criterion
irefined=0 !If 1, that means this combining process takes effects
idebugclusdegenatt=0 !If output debug information

if (imode==1) goto 2

if (allocated(attconv)) deallocate(attconv)
allocate(attconv(-2:numatt))
do iatt=-2,numatt !First assume that each attractor is a real attractor
    attconv(iatt)=iatt
end do

!Generate difference matrix for distance and relative value 
attdismat=0D0
attvalmat=0D0
do iatt=1,numatt
    valiatt=attval(iatt)
    xiatt=attxyz(iatt,1)
    yiatt=attxyz(iatt,2)
    ziatt=attxyz(iatt,3)
    do jatt=iatt+1,numatt
        valjatt=attval(jatt)
        xjatt=attxyz(jatt,1)
        yjatt=attxyz(jatt,2)
        zjatt=attxyz(jatt,3)
        distdiff=dsqrt( (xiatt-xjatt)**2+(yiatt-yjatt)**2+(ziatt-zjatt)**2 )
        attdismat(iatt,jatt)=distdiff
        valdiff=abs((valiatt-valjatt)/valiatt)
        attvalmat(iatt,jatt)=valdiff
        if (valdiff<valcritclus.and.distdiff<distcrit) irefined=1
!         write(10,"(2i5,f16.10)") iatt,jatt,attvalmat(iatt,jatt)
    end do
end do
attdismat=attdismat+transpose(attdismat)
attvalmat=attvalmat+transpose(attvalmat)
do iatt=1,numatt
    attdismat(iatt,iatt)=attdismat(iatt,iatt)/2D0
    attvalmat(iatt,iatt)=attvalmat(iatt,iatt)/2D0
end do
if (irefined==0) then !Nothing to do. Construct real attractor information and then directly return
    numrealatt=numatt
    if (allocated(nrealatthas)) deallocate(nrealatthas,realattval,realattxyz,realatttable)
    allocate(nrealatthas(numrealatt),realattval(numrealatt),realattxyz(numrealatt,3),realatttable(numrealatt,1))
    nrealatthas=1
    do iatt=1,numatt
        realattxyz(iatt,:)=attxyz(iatt,:)
        realattval(iatt)=attval(iatt)
        realatttable(iatt,1)=iatt
    end do
    return 
else
    write(*,*) "Degenerate attractors detected, clustering them..."
    write(*,"(' Criterion for clustering: Interval <',f8.5,' Bohr, value difference <',f8.5,'%')") distcrit,valcritclus*100
end if
! do i=1,numatt !Print distance matrix
!     do j=1,numatt
!         write(10,*) i,j,attdismat(i,j)
!     end do
! end do

!Clustering via nearest neighbour method, and accordingly update conversion list
passlist=0
do iatt=1,numatt
    if (passlist(iatt)==1) cycle
    if (count(attdismat(iatt,:)<distcrit)==1) then !This attractor has no neighbour
        if (passlist(iatt)==1) cycle
        cycle
    end if
    !Construct mapping table of attractor, capacity is neleatt
    neleatt=0 !The number of degenerate attractors in this batch
    do jatt=iatt,numatt
        if (attvalmat(iatt,jatt)<valcritclus) then
            neleatt=neleatt+1
            eleatt(neleatt)=jatt
            passlist(jatt)=1
        end if
    end do
    if (neleatt==1) cycle !Although this attractor has one or more neighbours, but no other attractor is degenerate with itself, so pass
    if (idebugclusdegenatt==1) then
        write(*,"(a,i5)") " Attractors in this time of clustering triggered by attractor:",iatt
        write(*,"(' The number of initial clusters:',i5,', corresponding to these attractor:')") neleatt
        write(*,"(15i5)") eleatt(1:neleatt)
    end if

    !Start clustering. Now you should forget actual index of attractors. Only consider element 1,2,3...
    !Initialize, assume that each element in this batch of degenerate attractors is a cluster
    ncluster=neleatt
    ncluele(:)=1
    do iclu=1,ncluster
        cluele(iclu,1)=iclu
    end do
    ntimeclu=0
    do while(.true.) !Gradually clustering until closet distance between any cluster pair is larger than criteria
        ntimeclu=ntimeclu+1
        imergeclu=0
        do iclu=1,ncluster !Note that ncluster remain unchanged, eliminated cluster is simply emptyed but still reserve its index
            if (ncluele(iclu)==0) cycle !Has been incorporated to others
            shortmax=9999999D0
            ishortmax=0
            do jclu=iclu+1,ncluster
                if (ncluele(jclu)==0) cycle
                !Calculate shortest distance between cluster i and j
                shortmaxtmp=9999999D0
                do ieletmp=1,ncluele(iclu)
                    do jeletmp=1,ncluele(jclu)
                        iele=cluele(iclu,ieletmp)
                        jele=cluele(jclu,jeletmp)
                        distele=attdismat(eleatt(iele),eleatt(jele))
                        if (distele<shortmaxtmp) shortmaxtmp=distele
                    end do
                end do
                if (shortmaxtmp<shortmax) then !Find out the closest cluster to cluster iclu
                    shortmax=shortmaxtmp
                    ishortmax=jclu
                end if
            end do
!             write(*,"(i3,f16.6,4i5,i2)") ntimeclu,shortmax,iclu,ishortmax,abs(shortmax<distcrit)
            if (shortmax<distcrit) then !Combine cluster ishortmax to cluster iclu
                imergeclu=1
                newlen=ncluele(iclu)+ncluele(ishortmax)
                cluele(iclu,ncluele(iclu)+1:newlen)=cluele(ishortmax,1:ncluele(ishortmax))
                ncluele(iclu)=newlen
                ncluele(ishortmax)=0 !Empty
            end if
        end do
        if (imergeclu==0) exit !Between all cluster pair the shortmax is longer than distcrit, therefore exit
    end do
    do iclu=1,ncluster !Update attractor conversion list
        if (ncluele(iclu)==0) cycle
        if (idebugclusdegenatt==1) then
            write(*,"(' Clustering finished after',i5,' times of cycle')") ntimeclu
            write(*,"(' Cluster',i5,', deriving from  attractor',i5,', contains:')") iclu,eleatt(iclu)
             write(*,"(' Element:  ',15i5)") cluele(iclu,1:ncluele(iclu))
            write(*,"(' Attractor:',15i5)") eleatt(cluele(iclu,1:ncluele(iclu)))
        end if
        do itmp=1,ncluele(iclu)
            ielimele=cluele(iclu,itmp)
            ielimatt=eleatt(ielimele)
            itargetatt=eleatt(iclu)
            attconv(ielimatt)=itargetatt
        end do
    end do
end do

! do iatt=1,numatt !Print conversion list just after clustering
!     write(*,*) iatt,attconv(iatt)
! end do
!Update attractor conversion list to make the attractor index contiguous, and hence get "real" attractor
!e.g. old conversion list, after slash is the new one
!   1           1/1
!   2           2/2
!   ...         2/2
!  13           2/2
!  14          14/3 <---After above clustering, 14,16,17 will convert to 14. The first new entry at right column always matched corresponding left column, this is guaranteed by above clustering code
!  15           2/2
!  16          14/3
!  17          14/3
!  18          18/4
2 numrealatt=0
passlist=0 !Record which attractor has already been converted
do iatt=1,numatt
    if (passlist(iatt)==1) cycle
    numrealatt=numrealatt+1 !Of course, numrealatt always <=iatt, so will not unexpectly overlap data during update conversion list
    idestmp=attconv(iatt)
    do jatt=iatt,numatt
        if (attconv(jatt)==idestmp) then
            passlist(jatt)=1
            attconv(jatt)=numrealatt
        end if
    end do
end do
! do iatt=1,numatt !Print conversion list after reordered
!     write(*,*) iatt,attconv(iatt)
! end do

call updaterealattprop

if (imode==1) return
!Update basin attribution
do iz=2,nz-1
    do iy=2,ny-1
        do ix=2,nx-1
            gridbas(ix,iy,iz)=attconv(gridbas(ix,iy,iz))
        end do
    end do
end do
!Output final attractors
write(*,*) "The attractors after clustering:"
write(*,*) "   Index      Average X,Y,Z coordinate (Angstrom)               Value"
do irealatt=1,numrealatt
    write(*,"(i8,3f15.8,f20.8)") irealatt,realattxyz(irealatt,:)*b2a,realattval(irealatt)
end do
end subroutine



!!------- Update real attractors' properties (average value, xyz, the number of members)
subroutine updaterealattprop
use defvar
use basinintmod
implicit real*8(a-h,o-z)
if (allocated(nrealatthas)) deallocate(nrealatthas,realattval,realattxyz,realatttable)
allocate(nrealatthas(numrealatt),realattval(numrealatt),realattxyz(numrealatt,3),realatttable(numrealatt,numatt))
nrealatthas=0
realattval=0D0
realattxyz=0D0
do iatt=1,numatt
    irealatt=attconv(iatt)
    nrealatthas(irealatt)=nrealatthas(irealatt)+1
    realatttable(irealatt,nrealatthas(irealatt))=iatt
    realattxyz(irealatt,:)=realattxyz(irealatt,:)+attxyz(iatt,:)
    realattval(irealatt)=realattval(irealatt)+attval(iatt)
end do
do irealatt=1,numrealatt
    realattval(irealatt)=realattval(irealatt)/nrealatthas(irealatt)
    realattxyz(irealatt,:)=realattxyz(irealatt,:)/nrealatthas(irealatt)
end do
end subroutine



!!------------------ Set grid for generating basin, adapted from the subroutine "setgrid"
subroutine setgridforbasin(ifuncsel)
use defvar
implicit real*8 (a-h,o-z)
integer :: iselexttype=3,ifuncsel
real*8 :: molxlen,molylen,molzlen,tmpx,tmpy,tmpz,rhocrit=1D-6
real*8 :: gridextdist=5D0,enlarbox=2.1D0,spclowqual=0.2D0,spcmedqual=0.1D0,spchighqual=0.06D0,spclunaqual=0.04D0,tmparr6(6)
character c80tmp*80,cubefilename*200
if (.not.allocated(b)) iselexttype=2 !Unable to set box using rho
do while(.true.)
    if (iselexttype==1) then
        orgx=minval(a%x)-gridextdist
        orgy=minval(a%y)-gridextdist
        orgz=minval(a%z)-gridextdist
        endx=maxval(a%x)+gridextdist
        endy=maxval(a%y)+gridextdist
        endz=maxval(a%z)+gridextdist
    else if (iselexttype==2) then
        orgx=minval( a(:)%x-enlarbox*vdwr_tianlu(a(:)%index) )
        orgy=minval( a(:)%y-enlarbox*vdwr_tianlu(a(:)%index) )
        orgz=minval( a(:)%z-enlarbox*vdwr_tianlu(a(:)%index) )
        endx=maxval( a(:)%x+enlarbox*vdwr_tianlu(a(:)%index) )
        endy=maxval( a(:)%y+enlarbox*vdwr_tianlu(a(:)%index) )
        endz=maxval( a(:)%z+enlarbox*vdwr_tianlu(a(:)%index) )
    end if
    molxlen=endx-orgx
    molylen=endy-orgy
    molzlen=endz-orgz
    ntotlow=(nint(molxlen/spclowqual)+1)*(nint(molylen/spclowqual)+1)*(nint(molzlen/spclowqual)+1)
    ntotmed=(nint(molxlen/spcmedqual)+1)*(nint(molylen/spcmedqual)+1)*(nint(molzlen/spcmedqual)+1)
    ntothigh=(nint(molxlen/spchighqual)+1)*(nint(molylen/spchighqual)+1)*(nint(molzlen/spchighqual)+1)
    ntotluna=(nint(molxlen/spclunaqual)+1)*(nint(molylen/spclunaqual)+1)*(nint(molzlen/spclunaqual)+1)
    
    write(*,*) "Please select a method for setting up grid"
    if (iselexttype==1) write(*,"(a,f10.5,a)") " -10 Set grid extension distance for mode 1~6, current: Fixed,",gridextdist," Bohr"
    if (iselexttype==2) write(*,"(a)") " -10 Set grid extension distance for mode 1~6, current: Adaptive"
    if (iselexttype==3) write(*,"(a)") " -10 Set grid extension distance for mode 1~6, current: Detect rho isosurface"
    if (iselexttype==1.or.iselexttype==2) then
        write(*,"(a,f4.2,a,i14)") " 1 Low quality grid, spacing=",spclowqual," Bohr, number of grids:    ",ntotlow
        write(*,"(a,f4.2,a,i14)") " 2 Medium quality grid, spacing=",spcmedqual," Bohr, number of grids: ",ntotmed
        write(*,"(a,f4.2,a,i14)") " 3 High quality grid, spacing=",spchighqual," Bohr, number of grids:   ",ntothigh
        write(*,"(a,f4.2,a,i14)") " 4 Lunatic quality grid, spacing=",spclunaqual," Bohr, number of grids:",ntotluna
    else
        write(*,"(a,f4.2,a,i14)") " 1 Low quality grid, spacing=",spclowqual," Bohr, cost: 1x"
        write(*,"(a,f4.2,a,i14)") " 2 Medium quality grid, spacing=",spcmedqual," Bohr, cost: 8x"
        write(*,"(a,f4.2,a,i14)") " 3 High quality grid, spacing=",spchighqual," Bohr, cost: 36x"
        write(*,"(a,f4.2,a,i14)") " 4 Lunatic quality grid, spacing=",spclunaqual," Bohr, cost: 120x"
    end if
    write(*,*) "5 Only input grid spacing, automatically set other parameters"
    write(*,*) "6 Only input the number of points in X,Y,Z, automatically set other parameters"
    write(*,*) "7 Input original point, translation vector and the number of points"
    write(*,*) "8 Set center position, grid spacing and box length"
    write(*,*) "9 Use grid setting of another cube file"
    
    read(*,*) igridsel
    if (igridsel/=-10) then
        exit
    else
        write(*,*) "Please select the type of extension distance:"
        write(*,*) "1 Use fixed value in each direction"
        write(*,*) "2 Adaptively determined according to van der Waals radius"
        write(*,"(a,1PE7.1)") " 3 Detect rho to make the box just accommodate its isosurface, isoval: ",rhocrit
        write(*,*) "Hint: In general mode 3 is recommended, mode 1 should be avoided to be used"
        read(*,*) iselexttype
        if (iselexttype==1) then
            do while(.true.)
                write(*,*) "Input extension distance (Bohr) e.g. 6.5"
                read(*,*) gridextdist
                if (gridextdist>0D0) exit
                write(*,*) "Error: The value must be larger than 0!" 
            end do
        else if (iselexttype==2) then
            do while(.true.)
                write(*,*) "Input the ratio for scaling vdW radii, e.g. 1.7"
                write(*,"(' Current value is',f12.6)") enlarbox
                read(*,*) enlarbox
                if (enlarbox>0D0) exit
                write(*,*) "Error: The value must be larger than 0!" 
            end do
        else if (iselexttype==3) then
            write(*,*) "Select the isovalue of rho"
            write(*,*) "1: 0.0001       (10^-4)"
            write(*,*) "2: 0.00001      (10^-5)"
            write(*,*) "3: 0.000005   (5*10^-6)"
            write(*,*) "4: 0.000003   (3*10^-6)"
            write(*,*) "5: 0.000001     (10^-6) In general recommended"
            write(*,*) "6: 0.0000005  (5*10^-7)"
            write(*,*) "7: 0.0000001    (10^-7)"
            write(*,*) "8: 0.00000005 (5*10^-8)"
            write(*,*) "9: 0.00000001   (10^-8)"
            write(*,*) "10: Input by yourself"
            read(*,*) iselrho
            if (iselrho==1) rhocrit=1D-4
            if (iselrho==2) rhocrit=1D-5
            if (iselrho==3) rhocrit=5D-6
            if (iselrho==4) rhocrit=3D-6
            if (iselrho==5) rhocrit=1D-6
            if (iselrho==6) rhocrit=5D-7
            if (iselrho==7) rhocrit=1D-7
            if (iselrho==8) rhocrit=5D-8
            if (iselrho==9) rhocrit=1D-8
            if (iselrho==10) then
                write(*,*) "Input the value, e.g. 0.000005 or 5D-6"
                read(*,*) rhocrit
            end if
        end if
    end if
end do

if (igridsel==1.or.igridsel==2.or.igridsel==3.or.igridsel==4.or.igridsel==5.or.igridsel==6) then
    if (iselexttype==1) then
        aug3D=gridextdist !Set aug3D, which only affects GUI display, so that the scope of the axis is wide enough to contain the whole grid region
    else if (iselexttype==2) then
        aug3D=maxval(enlarbox*vdwr_tianlu(a(:)%index))*1.35D0
    else if (iselexttype==3) then
        call detectboxrho(rhocrit) !determine orgx/y/z and endx/y/z
        molxlen=endx-orgx
        molylen=endy-orgy
        molzlen=endz-orgz
        tmparr6(1)=abs(orgx-minval(a%x))
        tmparr6(2)=abs(orgy-minval(a%y))
        tmparr6(3)=abs(orgz-minval(a%z))
        tmparr6(4)=abs(endx-maxval(a%x))
        tmparr6(5)=abs(endy-maxval(a%y))
        tmparr6(6)=abs(endz-maxval(a%z))
        aug3D=maxval(tmparr6)+1
    end if
end if

!Note: orgx,orgy,orgz,endx,endy,endz as well as molx/y/zlen for igridsel==1~6 have already been set above
if (igridsel==1.or.igridsel==2.or.igridsel==3.or.igridsel==4.or.igridsel==5) then
    if (igridsel==1) dx=spclowqual
    if (igridsel==2) dx=spcmedqual
    if (igridsel==3) dx=spchighqual
    if (igridsel==4) dx=spclunaqual
    if (igridsel==5) then
        write(*,*) "Input the grid spacing (bohr)  e.g. 0.08"
        read(*,*) dx
    end if
    dy=dx
    dz=dx
    nx=nint(molxlen/dx)+1
    ny=nint(molylen/dy)+1
    nz=nint(molzlen/dz)+1
else if (igridsel==6) then
    write(*,*) "Input the number of grid points in X,Y,Z direction   e.g. 139,59,80"
    read(*,*) nx,ny,nz
    dx=molxlen/(nx-1)
    dy=molylen/(ny-1)
    dz=molzlen/(nz-1)
else if (igridsel==7) then
    write(*,*) "Input X,Y,Z coordinate of original point (Bohr) e.g. 0.1,4,-1"
    read(*,*) orgx,orgy,orgz
    write(*,*) "Input X,Y,Z component of translation vector (Bohr) e.g. 0.1,0.1,0.15"
    read(*,*) dx,dy,dz
    write(*,*) "Input the number of points in X,Y,Z direction e.g. 139,59,80"
    read(*,*) nx,ny,nz
    endx=orgx+dx*(nx-1)
    endy=orgy+dy*(ny-1)
    endz=orgz+dz*(nz-1)
else if (igridsel==8) then
    write(*,*) "Input X,Y,Z coordinate of box center (in Angstrom)"
    write(*,*) "or input such as a8 to take the coordinate of atom 8 as box center"
    write(*,*) "or input such as a3,a7 to take the midpoint of atom 3 and atom 7 as box center"
    read(*,"(a)") c80tmp
    if (c80tmp(1:1)=='a') then
        do ich=1,len_trim(c80tmp)
            if (c80tmp(ich:ich)==',') exit
        end do
        if (ich==len_trim(c80tmp)+1) then
            read(c80tmp(2:),*) itmp
            cenx=a(itmp)%x
            ceny=a(itmp)%y
            cenz=a(itmp)%z
        else
            read(c80tmp(2:ich-1),*) itmp
            read(c80tmp(ich+2:),*) jtmp            
            cenx=(a(itmp)%x+a(jtmp)%x)/2D0
            ceny=(a(itmp)%y+a(jtmp)%y)/2D0
            cenz=(a(itmp)%z+a(jtmp)%z)/2D0
        end if
    else
        read(c80tmp,*) cenx,ceny,cenz
        cenx=cenx/b2a
        ceny=ceny/b2a
        cenz=cenz/b2a
    end if
    write(*,*) "Input the grid spacing (bohr)  e.g. 0.08"
    read(*,*) dx
    dy=dx
    dz=dx
    write(*,*) "Input the box lengths in X,Y,Z direction (Bohr) e.g. 8.0,8.0,13.5"
    read(*,*) molxlen,molylen,molzlen
    orgx=cenx-molxlen/2D0
    orgy=ceny-molylen/2D0
    orgz=cenz-molzlen/2D0
    endx=orgx+molxlen
    endy=orgy+molylen
    endz=orgz+molzlen
    nx=nint(molxlen/dx)+1
    ny=nint(molylen/dy)+1
    nz=nint(molzlen/dz)+1
else if (igridsel==9) then
    write(*,*) "Input filename of a cube file"
    do while(.true.)
        read(*,"(a)") cubefilename
        inquire(file=cubefilename,exist=alive)
        if (alive) then
            open(10,file=cubefilename,status="old")
            read(10,*)
            read(10,*)
            read(10,*) nouse,orgx,orgy,orgz
            read(10,*) nx,dx
            read(10,*) ny,rnouse,dy
            read(10,*) nz,rnouse,rnouse,dz
            close(10)
            exit
        else
            write(*,*) "Error: File cannot be found, input again"
        end if
    end do
    endx=orgx+dx*(nx-1)
    endy=orgy+dy*(ny-1)
    endz=orgz+dz*(nz-1)
end if

!If system is symmetric to Cartesian plane, then slightly adjust grid setting, so that the distribution of grids are symmetric to Cartesian plane
!This treatment will make the integrals have much better symmetry
diffx=abs(orgx+endx)
diffy=abs(orgy+endy)
diffz=abs(orgz+endz)
if (igridsel>=1.and.igridsel<=6) then
    if (diffx<0.05D0) then !The system is symmetry to YZ plane
        distxmin=1D10
        do ix=1,nx
            rnowx=orgx+(ix-1)*dx
            if (rnowx>=0D0.and.rnowx<distxmin) distxmin=rnowx
        end do
        !1D-12 is a perturbation to avoid the grids are degenerate with respect to Cartesian plane, which results in cumbersome degenerate attractors
        !However if it is introduced, the integrals will not so close to symmetry
        orgx=orgx+(dx/2D0-distxmin) !+1D-15
    end if
    if (diffy<0.05D0) then !The system is symmetry to XZ plane
        distymin=1D10
        do iy=1,ny
            rnowy=orgy+(iy-1)*dy
            if (rnowy>=0D0.and.rnowy<distymin) distymin=rnowy
        end do
        orgy=orgy+(dy/2D0-distymin) !+1D-15
    end if
    if (diffz<0.05D0) then !The system is symmetry to XY plane
        distzmin=1D10
        do iz=1,nz
            rnowz=orgz+(iz-1)*dz
            if (rnowz>=0D0.and.rnowz<distzmin) distzmin=rnowz
        end do
        orgz=orgz+(dz/2D0-distzmin) !+1D-15
    end if
end if

write(*,"(' Coordinate of origin in X,Y,Z is   ',3f12.6)") orgx,orgy,orgz
write(*,"(' Coordinate of end point in X,Y,Z is',3f12.6)") endx,endy,endz
write(*,"(' Spacing in X,Y,Z is',3f11.6)") dx,dy,dz
write(*,"(' Number of points in X,Y,Z is',3i5,'   Total',i10)") nx,ny,nz,nx*ny*nz
! do i=1,nx
!     write(c80tmp,"(D20.13)") orgx+(i-1)*dx
!     read(c80tmp,*) tmpval
!     write(*,*) i," x ",tmpval
! end do
! do i=1,ny
!     write(c80tmp,"(D20.13)") orgy+(i-1)*dy
!     read(c80tmp,*) tmpval
!     write(*,*) i," y ",tmpval
! end do
! do i=1,nz
!     write(c80tmp,"(D20.13)") orgz+dfloat(i-1)*dz
!     read(c80tmp,*) tmpval
!     write(*,*) i," z ",tmpval
! end do
! read(*,*)
end subroutine


!!-- Detect proper box (namely set orgx/y/z,endx/y/z) so that the box can just enclose the isosurface of rho=rhocrit
subroutine detectboxrho(rhocrit)
use defvar
use function
implicit real*8(a-h,o-z)
real*8 rhocrit
tmpfac=1.0D0
orgxtmp=minval( a(:)%x-tmpfac*vdwr_tianlu(a(:)%index) ) !Define a too small box, extend each side (to obtain orgx/y/z,endx/y/z) by detecting rho
orgytmp=minval( a(:)%y-tmpfac*vdwr_tianlu(a(:)%index) )
orgztmp=minval( a(:)%z-tmpfac*vdwr_tianlu(a(:)%index) )
endxtmp=maxval( a(:)%x+tmpfac*vdwr_tianlu(a(:)%index) )
endytmp=maxval( a(:)%y+tmpfac*vdwr_tianlu(a(:)%index) )
endztmp=maxval( a(:)%z+tmpfac*vdwr_tianlu(a(:)%index) )
scanstp=0.3D0 !Spacing of scan
nstpx=nint((endxtmp-orgxtmp)/scanstp)+1 !The number of detection grid in X direction
nstpy=nint((endytmp-orgytmp)/scanstp)+1 !The number of detection grid in Y direction
nstpz=nint((endztmp-orgztmp)/scanstp)+1 !The number of detection grid in Z direction
distmove=0.2D0
write(*,*) "Detecting proper box size..."
!Detect Z low. Scan XY plane, if rho at a point is larger than "rhocrit", lower down the plane by "distmove" and check again, until all points have rho<rhocrit
orgz=orgztmp
do while(.true.)
    iok=1
zl:do ix=1,nstpx
        rnowx=orgxtmp+(ix-1)*scanstp
        do iy=1,nstpy
            rnowy=orgytmp+(iy-1)*scanstp
            if (fdens(rnowx,rnowy,orgz)>rhocrit) then
                iok=0
                exit zl
            end if
        end do
    end do zl
    if (iok==1) exit
    orgz=orgz-0.2D0 !Lower down the layer
end do
!Detect Z high
endz=endztmp
do while(.true.)
    iok=1
zh:do ix=1,nstpx
        rnowx=orgxtmp+(ix-1)*scanstp
        do iy=1,nstpy
            rnowy=orgytmp+(iy-1)*scanstp
            if (fdens(rnowx,rnowy,endz)>rhocrit) then
                iok=0
                exit zh
            end if
        end do
    end do zh
    if (iok==1) exit
    endz=endz+0.2D0
end do
!Detect X low. Scan YZ plane
orgx=orgxtmp
do while(.true.)
    iok=1
xl:do iy=1,nstpy
        rnowy=orgytmp+(iy-1)*scanstp
        do iz=1,nstpz
            rnowz=orgztmp+(iz-1)*scanstp
            if (fdens(orgx,rnowy,rnowz)>rhocrit) then
                iok=0
                exit xl
            end if
        end do
    end do xl
    if (iok==1) exit
    orgx=orgx-0.2D0
end do
!Detect X high
endx=endxtmp
do while(.true.)
    iok=1
xh:do iy=1,nstpy
        rnowy=orgytmp+(iy-1)*scanstp
        do iz=1,nstpz
            rnowz=orgztmp+(iz-1)*scanstp
            if (fdens(endx,rnowy,rnowz)>rhocrit) then
                iok=0
                exit xh
            end if
        end do
    end do xh
    if (iok==1) exit
    endx=endx+0.2D0
end do
!Detect Y low. Scan XZ plane
orgy=orgytmp
do while(.true.)
    iok=1
yl:do ix=1,nstpx
        rnowx=orgxtmp+(ix-1)*scanstp
        do iz=1,nstpz
            rnowz=orgztmp+(iz-1)*scanstp
            if (fdens(rnowx,orgy,rnowz)>rhocrit) then
                iok=0
                exit yl
            end if
        end do
    end do yl
    if (iok==1) exit
    orgy=orgy-0.2D0
end do
!Detect Y high
endy=endytmp
do while(.true.)
    iok=1
yh:do ix=1,nstpx
        rnowx=orgxtmp+(ix-1)*scanstp
        do iz=1,nstpz
            rnowz=orgztmp+(iz-1)*scanstp
            if (fdens(rnowx,endy,rnowz)>rhocrit) then
                iok=0
                exit yh
            end if
        end do
    end do yh
    if (iok==1) exit
    endy=endy+0.2D0
end do

end subroutine




!!------- Integrate a real space function in the basins already partitioned
subroutine integratebasin
use defvar
use util
use function
use basinintmod
implicit real*8(a-h,o-z)
real*8 intval(-1:numatt),basinvol(-1:numatt),intvalpriv(-1:numatt),basinvolpriv(-1:numatt)
integer walltime1,walltime2
character grdfilename*200

write(*,*) "Please select integrand:"
write(*,*) "-2 Return"
write(*,*) "-1 The values of the grid data stored in an external file (.cub/.grd)"
write(*,*) "0 The values of the grid data stored in memory"
call selfunc_interface(ifuncint)
if (ifuncint==-1) then
    do while(.true.)
        write(*,*) "Input another .cub or .grd file name"
        read(*,"(a)") grdfilename
        inquire(file=grdfilename,exist=alive)
        if (alive) exit
        write(*,*) "File not found, input again"
        write(*,*)
    end do
    inamelen=len_trim(grdfilename)
    if (grdfilename(inamelen-2:inamelen)=="cub".or.grdfilename(inamelen-3:inamelen)=="cube") then
        call readcubetmp(grdfilename,inconsis)
    else if (grdfilename(inamelen-2:inamelen)=="grd") then
        call readgrdtmp(grdfilename,inconsis)
    end if
    if (inconsis==1) return
    write(*,*)
else if (ifuncint==-2) then
    return
end if

call walltime(walltime1)
CALL CPU_TIME(time_begin)
write(*,*) "Integrating, please wait..."
intval=0D0
basinvol=0D0
ifinish=0
nthreads=getNThreads()
!$OMP PARALLEL private(ix,iy,iz,irealatt,rnowx,rnowy,rnowz,tmpval,intvalpriv,basinvolpriv) shared(intval,basinvol,ifinish) NUM_THREADS(nthreads)
intvalpriv=0D0
basinvolpriv=0D0
!$OMP do schedule(DYNAMIC)
do iz=2,nz-1
    do iy=2,ny-1
        do ix=2,nx-1
            rnowx=orgx+(ix-1)*dx
            rnowy=orgy+(iy-1)*dy
            rnowz=orgz+(iz-1)*dz
            if (ifuncint==-1) then
                tmpval=cubmattmp(ix,iy,iz)
            else if (ifuncint==0) then
                tmpval=cubmat(ix,iy,iz)
            else
                tmpval=calcfuncall(ifuncint,rnowx,rnowy,rnowz)
            end if
            irealatt=gridbas(ix,iy,iz)
            intvalpriv(irealatt)=intvalpriv(irealatt)+tmpval
            basinvolpriv(irealatt)=basinvolpriv(irealatt)+1
        end do
    end do
    ifinish=ifinish+1
    if (ifuncint>=1) write(*,"(' Finished:',i5,'/',i5)") ifinish,nz
end do
!$OMP end do
!$OMP CRITICAL
    intval=intval+intvalpriv
    basinvol=basinvol+basinvolpriv
!$OMP end CRITICAL
!$OMP END PARALLEL
dvol=dx*dy*dz
intval=intval*dvol
basinvol=basinvol*dvol !Basin volume
write(*,*) "  #Basin        Integral(a.u.)      Volume(a.u.^3)"
do irealatt=1,numrealatt
    write(*,"(i8,f22.10,f20.8)") irealatt,intval(irealatt),basinvol(irealatt)
end do
write(*,"(' Sum of above values:',f20.8)") sum(intval(1:numrealatt))
if (any(gridbas(2:nx-1,2:ny-1,2:nz-1)==0)) write(*,"(' Integral of unassigned grids:',f20.8)") intval(0)
if (any(gridbas(2:nx-1,2:ny-1,2:nz-1)==-1)) write(*,"(' Integral of the grids travelled to box boundary:',f20.8)") intval(-1)

CALL CPU_TIME(time_end)
call walltime(walltime2)
write(*,"(' Integrating basins took up CPU time',f12.2,'s, wall clock time',i10,'s')") time_end-time_begin,walltime2-walltime1
end subroutine



!!------- Obtain atomic population in a basin region (e.g. ELF bond basin) via AIM partition
subroutine atmpopinbasin
use defvar
use basinintmod
implicit real*8(a-h,o-z)
inquire(file="basin.cub",exist=alive)
if (alive .eqv. .false.) then
    write(*,*) "Error: basin.cub is not existed in current folder!"
    return
else
    if (allocated(cubmattmp)) deallocate(cubmattmp)
    call readcubetmp("basin.cub",inconsis)
    if (inconsis==1) then
        write(*,*) "Error: The grid setting of basin.cub is inconsistent with present grid data"
        return
    end if
end if
do while(.true.)
    write(*,*)
    write(*,*) "Study population of which atom? Input the index of corresponding attractor"
    write(*,*) "For example, 3"
    write(*,*) "Input 0 can exit"
    read(*,*) iatt
    if (iatt==0) then
        exit
    else if (iatt<0.or.iatm>numrealatt) then
        write(*,*) "Error: Attractor index exceeded valid range!"
        return
    end if
    write(*,*) "Study its population in which basin of the basin.cub file? e.g. 6"
    read(*,*) ibasin
    atmpop=0
    do iz=2,nz-1
        do iy=2,ny-1
            do ix=2,nx-1
                if (nint(cubmattmp(ix,iy,iz))==ibasin.and.gridbas(ix,iy,iz)==iatt) atmpop=atmpop+cubmat(ix,iy,iz)
            end do
        end do
    end do
    atmpop=atmpop*dx*dy*dz
    write(*,"(' Population of attractor',i4,' in external basin',i4,' is',f12.5)") iatt,ibasin,atmpop
end do
deallocate(cubmattmp)
end subroutine



!!----------------- Calculate multipole moment in the basins
subroutine multipolebasin
use defvar
use basinintmod
use function
implicit real*8 (a-h,o-z)
do while(.true.)
    write(*,*) "Input the index of the basin in question, e.g. 5"
    write(*,"(a)") " Note: -1 means printing all basin results on screen, -2 means printing to multipol.txt in current folder. Input 0 can return"
    read(*,*) itmp
    if (itmp==-1.or.itmp==-2) then
        ibegin=1
        iend=numrealatt
        if (itmp==-2) then
            ioutid=10
            open(10,file="multipol.txt",status="replace")
        end if
    else if (itmp==0) then
        return
    else
        ibegin=itmp
        iend=itmp
        ioutid=6
    end if
    write(*,*) "Calculating, please wait..."
    write(ioutid,*) "Note: All units shown below are in a.u.!"
    write(ioutid,*)
    
    dipelextot=0D0
    dipeleytot=0D0
    dipeleztot=0D0
    eleinttot=0D0
    do ibas=ibegin,iend
    !     xcen=0 !Use centroid of basin as center, but this is a bad idea
    !     ycen=0
    !     zcen=0
    !     nbasgrid=0
    !     do iz=2,nz-1
    !         do iy=2,ny-1
    !             do ix=2,nx-1
    !                 if (gridbas(ix,iy,iz)==ibas) then
    !                     nbasgrid=nbasgrid+1
    !                     rnowx=orgx+(ix-1)*dx
    !                     rnowy=orgy+(iy-1)*dy
    !                     rnowz=orgz+(iz-1)*dz
    !                     tmpdens=ELF_LOL(rnowx,rnowy,rnowz,"ELF")
    !                     xcen=xcen+rnowx*tmpdens
    !                     ycen=ycen+rnowy*tmpdens
    !                     zcen=zcen+rnowz*tmpdens
    !                 end if
    !             end do
    !         end do
    !     end do
    !     xcen=xcen/nbasgrid
    !     ycen=ycen/nbasgrid
    !     zcen=zcen/nbasgrid
        xcen=realattxyz(ibas,1)
        ycen=realattxyz(ibas,2)
        zcen=realattxyz(ibas,3)
        dvol=dx*dy*dz
    !     write(*,"(' The X,Y,Z of the center of the basin:')")
    !     write(*,"(3f14.8,' Bohr',/)") xcen,ycen,zcen
    
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
nthreads=getNThreads()
!$OMP PARALLEL private(ix,iy,iz,rnowx,rnowy,rnowz,tmpmul,rx,ry,rz,eleintp,xintp,yintp,zintp,xxintp,yyintp,zzintp,xyintp,yzintp,xzintp) NUM_THREADS(nthreads)
        eleintp=0D0
        xintp=0D0
        yintp=0D0
        zintp=0D0
        xxintp=0D0
        yyintp=0D0
        zzintp=0D0
        xyintp=0D0
        yzintp=0D0
        xzintp=0D0
!$OMP do schedule(DYNAMIC)
        do iz=2,nz-1
            rnowz=orgz+(iz-1)*dz
            do iy=2,ny-1
                rnowy=orgy+(iy-1)*dy
                do ix=2,nx-1
                    if (gridbas(ix,iy,iz)==ibas) then
                        rnowx=orgx+(ix-1)*dx
                        if (ifuncbasin==1) then
                            tmpmul=dvol*cubmat(ix,iy,iz) !The cubmat currently is just electron density
                        else
                            tmpmul=dvol*fdens(rnowx,rnowy,rnowz)
                        end if
                        rx=rnowx-xcen
                        ry=rnowy-ycen
                        rz=rnowz-zcen
                        eleintp=eleintp+tmpmul !monopole
                        xintp=xintp+rx*tmpmul
                        yintp=yintp+ry*tmpmul
                        zintp=zintp+rz*tmpmul
                        xxintp=xxintp+rx*rx*tmpmul
                        yyintp=yyintp+ry*ry*tmpmul
                        zzintp=zzintp+rz*rz*tmpmul
                        xyintp=xyintp+rx*ry*tmpmul
                        yzintp=yzintp+ry*rz*tmpmul
                        xzintp=xzintp+rx*rz*tmpmul
                    end if
                end do
            end do
        end do
!$OMP end do
!$OMP CRITICAL
        eleint=eleint+eleintp
        xint=xint+xintp
        yint=yint+yintp
        zint=zint+zintp
        xxint=xxint+xxintp
        yyint=yyint+yyintp
        zzint=zzint+zzintp
        xyint=xyint+xyintp
        yzint=yzint+yzintp
        xzint=xzint+xzintP
!$OMP end CRITICAL
!$OMP END PARALLEL
        if (itmp==-1.or.itmp==-2) write(ioutid,"(' Basin',i8)") ibas
        if (itmp==-2) write(*,"(' Outputting basin',i8)") ibas
        write(ioutid,"(' Basin electric monopole moment:',f12.6)") -eleint
        write(ioutid,"(' Basin electric dipole moment:')") 
        write(ioutid,"('  X=',f12.6,'  Y=',f12.6,'  Z=',f12.6,'  Magnitude=',f12.6)") -xint,-yint,-zint,sqrt(xint**2+yint**2+zint**2)
        eleinttot=eleinttot+eleint
        dipelex=-eleint*xcen+(-xint) !Contribution to molecular total dipole moment
        dipeley=-eleint*ycen+(-yint)
        dipelez=-eleint*zcen+(-zint)
        dipelextot=dipelextot+dipelex
        dipeleytot=dipeleytot+dipeley
        dipeleztot=dipeleztot+dipelez
        write(ioutid,"(' Basin electron contribution to molecular dipole moment:')")
        write(ioutid,"('  X=',f12.6,'  Y=',f12.6,'  Z=',f12.6,'  Magnitude=',f12.6)") dipelex,dipeley,dipelez,sqrt(dipelex**2+dipeley**2+dipelez**2)
        write(ioutid,"(' Basin electric quadrupole moment (Cartesian form):')")
        rrint=xxint+yyint+zzint
        QXX=-(3*xxint-rrint)/2
        QYY=-(3*yyint-rrint)/2
        QZZ=-(3*zzint-rrint)/2
        write(ioutid,"(' QXX=',f12.6,'  QXY=',f12.6,'  QXZ=',f12.6)") QXX,-(3*xyint)/2,-(3*xzint)/2
        write(ioutid,"(' QYX=',f12.6,'  QYY=',f12.6,'  QYZ=',f12.6)") -(3*xyint)/2,QYY,-(3*yzint)/2
        write(ioutid,"(' QZX=',f12.6,'  QZY=',f12.6,'  QZZ=',f12.6)") -(3*xzint)/2,-(3*yzint)/2,QZZ
        write(ioutid,"( ' The magnitude of electric quadrupole moment (Cartesian form):',f12.6)") dsqrt(2D0/3D0*(QXX**2+QYY**2+QZZ**2))
        R20=-(3*zzint-rrint)/2D0 !Notice that the negative sign, because electrons carry negative charge
        R2n1=-dsqrt(3D0)*yzint
        R2p1=-dsqrt(3D0)*xzint
        R2n2=-dsqrt(3D0)*xyint
        R2p2=-dsqrt(3D0)/2D0*(xxint-yyint)
        write(ioutid,"(' Electric quadrupole moments (Spherical harmonic form):')")
        write(ioutid,"(' Q_2,0 =',f11.6,'   Q_2,-1=',f11.6,'   Q_2,1=',f11.6)") R20,R2n1,R2p1
        write(ioutid,"(' Q_2,-2=',f11.6,'   Q_2,2 =',f11.6)") R2n2,R2p2
        write(ioutid,"( ' Magnitude: |Q_2|=',f12.6)") dsqrt(R20**2+R2n1**2+R2p1**2+R2n2**2+R2p2**2)
        write(ioutid,*)
    end do
    if (itmp==-1.or.itmp==-2) then !Output overall properties, most users are not interested in them
        dipnucx=sum(a(:)%x*a(:)%charge)
        dipnucy=sum(a(:)%y*a(:)%charge)
        dipnucz=sum(a(:)%z*a(:)%charge)
        write(ioutid,"( ' Molecular net charge:',f12.6)") sum(a%charge)-eleinttot
        write(ioutid,"( ' Nuclear contribution to molecular dipole moment:')") 
        write(ioutid,"(' X=',f12.6,' Y=',f12.6,' Z=',f12.6,' Norm=',f12.6)") dipnucx,dipnucy,dipnucz,sqrt(dipnucx**2+dipnucy**2+dipnucz**2)
        write(ioutid,"( ' Electron contribution to molecular dipole moment:')") 
        write(ioutid,"(' X=',f12.6,' Y=',f12.6,' Z=',f12.6,' Norm=',f12.6)") dipelextot,dipeleytot,dipeleztot,sqrt(dipelextot**2+dipeleytot**2+dipeleztot**2)
        dipmolx=dipnucx+dipelextot
        dipmoly=dipnucy+dipeleytot
        dipmolz=dipnucz+dipeleztot
        write(ioutid,"( ' Molecular dipole moment:')")
        write(ioutid,"(' X=',f12.6,' Y=',f12.6,' Z=',f12.6,' Norm=',f12.6)") dipmolx,dipmoly,dipmolz,sqrt(dipmolx**2+dipmoly**2+dipmolz**2)
        write(ioutid,*)
    end if
    if (itmp==-2) then
        close(10)
        write(*,*) "Outputting finished!"
        write(*,*)
    end if
end do
end subroutine


!!------------------ Calculate localized index and delocalization index within and between basins
!Most of the codes are adapted from fuzzyana.f90
!itask==0 means calculate BOM and LI,DI, itask==1 means calculate and output BOM to BOM.txt in current folder
!BOM means overlap matrix of all occupied orbitals in basins
subroutine LIDIbasin(itask)
use defvar
use basinintmod
use function
use util
implicit real*8 (a-h,o-z)
real*8 DI(numrealatt,numrealatt),DIa(numrealatt,numrealatt),DIb(numrealatt,numrealatt) !Delocalization index matrix
real*8 LI(numrealatt),LIa(numrealatt),LIb(numrealatt) !Localization index array
real*8,allocatable :: BOM(:,:,:),BOMa(:,:,:),BOMb(:,:,:),BOMsum(:,:),BOMsuma(:,:),BOMsumb(:,:) !BOM(i,j,k) means overlap matrix of MO i,j in atom k space
real*8 :: BOMtmp(nmo,nmo),orbval(nmo)
character selectyn
integer itask

!Allocate space for BOM
if (wfntype==0.or.wfntype==2.or.wfntype==3) then !RHF,ROHF,R-post-HF
    if (wfntype==3) then !R-post-HF, need to consider all orbitals
        nmatsize=nmo
    else !RHF,ROHF
        !High-lying virtual orbitals will be deleted, especially for .fch case
        !Notice that occupation number may be not contiguous, some low-lying orbital may have
        !zero occupation due to modification by users, so we can't simply use nelec to determine matrix size
        do nmatsize=nmo,1,-1
            if (MOocc(nmatsize)/=0) exit
        end do
        if (nmo-nmatsize>0) write(*,"('Note: The highest',i6,' virtual orbitals will not be taken into account')") nmo-nmatsize
    end if
    if (allocated(BOM)) deallocate(BOM,BOMsum)
    allocate(BOM(nmatsize,nmatsize,numrealatt),BOMsum(nmatsize,nmatsize))
    BOM=0
    BOMsum=0
else if (wfntype==1.or.wfntype==4) then !UHF,U-post-HF
    do iendalpha=nmo,1,-1
        if (MOtype(iendalpha)==1) exit
    end do
    if (wfntype==4) then
        nmatsizea=iendalpha !Total number of alpha orbitals
        nmatsizeb=nmo-nmatsizea !Total number of beta orbitals
    else
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
        if (iendalpha-nmatsizea>0) write(*,"('Note: The highest',i6,' alpha virtual orbitals will not be taken into account')") iendalpha-nmatsizea
        if (nmo-iendalpha-nmatsizeb>0) write(*,"('Note: The highest',i6,' beta virtual orbitals will not be taken into account')") nmo-iendalpha-nmatsizeb
    end if
    if (allocated(BOMa)) deallocate(BOMa,BOMb,BOMsuma,BOMsumb)
    allocate( BOMa(nmatsizea,nmatsizea,numrealatt),BOMb(nmatsizeb,nmatsizeb,numrealatt) )
    allocate( BOMsuma(nmatsizea,nmatsizea),BOMsumb(nmatsizeb,nmatsizeb) )
    BOMa=0
    BOMb=0
    BOMsuma=0
    BOMsumb=0
end if

!Calculate BOM
dvol=dx*dy*dz
call walltime(nwalltime1)
if (wfntype==0.or.wfntype==2.or.wfntype==3) then !RHF,ROHF,R-post-HF
    do ibas=1,numrealatt
        write(*,"(' Generating orbital overlap matrix for basin',i6,'  of',i6,' ......')") ibas,numrealatt
nthreads=getNThreads()
!$OMP parallel shared(BOM) private(ix,iy,iz,rnowx,rnowy,rnowz,imo,jmo,BOMtmp,orbval) num_threads(nthreads)
        BOMtmp=0D0
!$OMP do schedule(DYNAMIC)
        do iz=2,nz-1
            do iy=2,ny-1
                do ix=2,nx-1
                    if (gridbas(ix,iy,iz)==ibas) then
                        rnowx=orgx+(ix-1)*dx
                        rnowy=orgy+(iy-1)*dy
                        rnowz=orgz+(iz-1)*dz
                        call orbderv(1,1,nmatsize,rnowx,rnowy,rnowz,orbval)
                        do imo=1,nmatsize
                            do jmo=imo,nmatsize
                                BOMtmp(imo,jmo)=BOMtmp(imo,jmo)+orbval(imo)*orbval(jmo)*dvol
                            end do
                        end do
                    end if
                end do
            end do
        end do
!$OMP end do
!$OMP CRITICAL
        BOM(:,:,ibas)=BOM(:,:,ibas)+BOMtmp(1:nmatsize,1:nmatsize)
!$OMP end CRITICAL
!$OMP end parallel
        BOM(:,:,ibas)=BOM(:,:,ibas)+transpose(BOM(:,:,ibas))
        do imo=1,nmatsize
            BOM(imo,imo,ibas)=BOM(imo,imo,ibas)/2D0
        end do
    end do
else if (wfntype==1.or.wfntype==4) then !UHF,U-post-HF
    do ibas=1,numrealatt
        !Alpha part
        write(*,"(' Generating orbital overlap matrix for basin',i6,'  of',i6,' ......')") ibas,numrealatt
nthreads=getNThreads()
!$OMP parallel shared(BOMa) private(ix,iy,iz,rnowx,rnowy,rnowz,imo,jmo,BOMtmp,orbval) num_threads(nthreads)
        BOMtmp=0D0
!$OMP do schedule(DYNAMIC)
        do iz=2,nz-1
            do iy=2,ny-1
                do ix=2,nx-1
                    if (gridbas(ix,iy,iz)==ibas) then
                        rnowx=orgx+(ix-1)*dx
                        rnowy=orgy+(iy-1)*dy
                        rnowz=orgz+(iz-1)*dz
                        call orbderv(1,1,nmatsizea,rnowx,rnowy,rnowz,orbval)
                        do imo=1,nmatsizea
                            do jmo=imo,nmatsizea
                                BOMtmp(imo,jmo)=BOMtmp(imo,jmo)+orbval(imo)*orbval(jmo)*dvol
                            end do
                        end do
                    end if
                end do
            end do
        end do
!$OMP end do
!$OMP CRITICAL
        BOMa(:,:,ibas)=BOMa(:,:,ibas)+BOMtmp(1:nmatsizea,1:nmatsizea)
!$OMP end CRITICAL
!$OMP end parallel
        BOMa(:,:,ibas)=BOMa(:,:,ibas)+transpose(BOMa(:,:,ibas))
        do imo=1,nmatsizea
            BOMa(imo,imo,ibas)=BOMa(imo,imo,ibas)/2D0
        end do
        
        !Beta part
        if (nmatsizeb>0) then !If there is no beta orbital, then nmatsizeb=0, and we will do nothing
            MOinit=iendalpha+1
            MOend=iendalpha+nmatsizeb
!             write(*,*) MOinit,MOend,nmatsizeb
nthreads=getNThreads()
!$OMP parallel shared(BOMb) private(ix,iy,iz,rnowx,rnowy,rnowz,imo,imotmp,jmo,jmotmp,BOMtmp,orbval) num_threads(nthreads)
            BOMtmp=0D0
!$OMP do schedule(DYNAMIC)
            do iz=2,nz-1
                do iy=2,ny-1
                    do ix=2,nx-1
                        if (gridbas(ix,iy,iz)==ibas) then
                            rnowx=orgx+(ix-1)*dx
                            rnowy=orgy+(iy-1)*dy
                            rnowz=orgz+(iz-1)*dz
                            call orbderv(1,MOinit,MOend,rnowx,rnowy,rnowz,orbval)
                            do imo=MOinit,MOend
                                imotmp=imo-iendalpha !So that the index start from 1 to nbelec
                                do jmo=imo,MOend
                                    jmotmp=jmo-iendalpha
                                    BOMtmp(imotmp,jmotmp)=BOMtmp(imotmp,jmotmp)+orbval(imo)*orbval(jmo)*dvol
                                end do
                            end do
                        end if
                    end do
                end do
            end do
!$OMP end do
!$OMP CRITICAL
            BOMb(:,:,ibas)=BOMb(:,:,ibas)+BOMtmp(1:nmatsizeb,1:nmatsizeb)
!$OMP end CRITICAL
!$OMP end parallel
            BOMb(:,:,ibas)=BOMb(:,:,ibas)+transpose(BOMb(:,:,ibas))
            do imo=1,nmatsizeb
                BOMb(imo,imo,ibas)=BOMb(imo,imo,ibas)/2D0
            end do
        end if
    end do
end if

write(*,*)
!Check sanity of BOM
if (wfntype==0.or.wfntype==2.or.wfntype==3) then !RHF,ROHF,R-post-HF
    do ibas=1,numrealatt
        BOMsum=BOMsum+BOM(:,:,ibas)
    end do
    BOMerror=identmaterr(BOMsum)/numrealatt
    write(*,"(' Error of BOM is',f14.8)") BOMerror
    if (BOMerror>0.02D0) write(*,"(a)") " Warning: The integration is not very accurate"
!     call showmatgau(BOMsum,"BOMsum",0,"f14.8",6)
else if (wfntype==1.or.wfntype==4) then !UHF,U-post-HF
    BOMerrorb=0D0
    do ibas=1,numrealatt
        BOMsuma=BOMsuma+BOMa(:,:,ibas)
        if (nmatsizeb>0) BOMsumb=BOMsumb+BOMb(:,:,ibas)
    end do
    BOMerrora=identmaterr(BOMsuma)/numrealatt
    if (nmatsizeb>0) BOMerrorb=identmaterr(BOMsumb)/numrealatt
    write(*,"(' Error of alpha BOM is',f14.8)") BOMerrora
    if (nmatsizeb>0) write(*,"(' Error of Beta BOM is ',f14.8)") BOMerrorb
    if (BOMerrora>0.02D0.or.BOMerrorb>0.02D0) write(*,"(a)") " Warning: The integration is not very accurate"
!     call showmatgau(BOMsumb,"BOMsumb",0,"f14.8",6)
end if

!Output BOM
if (itask==1) then
    open(10,file="BOM.txt",status="replace")
    if (wfntype==0.or.wfntype==2.or.wfntype==3) then
        do ibas=1,numrealatt
            write(10,"('Orbital overlap matrix of basin',i6)") ibas
            call showmatgau(BOM(:,:,ibas),"",1,"f14.8",10)
            write(10,*)
        end do
    else if (wfntype==1.or.wfntype==4) then
        do ibas=1,numrealatt
            write(10,"('Alpha part of orbital overlap matrix of basin',i6)") ibas
            call showmatgau(BOMa(:,:,ibas),"",1,"f14.8",10)
            if (nmatsizeb>0) then
                write(10,"('Beta part of orbital overlap matrix of basin',i6)") ibas
                call showmatgau(BOMb(:,:,ibas),"",1,"f14.8",10)
            end if
            write(10,*)
        end do
    end if
    close(10)
    write(*,*)
    write(*,*) "Done, the matrices have been exported to BOM.txt in current folder"
    return
end if

!Generate DI and LI
!RHF,R-post-HF, DI_A,B=2¡Æ[i,j]dsqrt(n_i*n_j)*S_i,j_A * S_i,j_B     where i and j are non-spin orbitals
write(*,*) "Generating LI and DI..."
if (wfntype==0.or.wfntype==3) then
    DI=0D0
    do ibas=1,numrealatt
        do jbas=ibas,numrealatt
            do iorb=1,nmatsize
                do jorb=1,nmatsize
                    DI(ibas,jbas)=DI(ibas,jbas)+dsqrt(MOocc(iorb)*MOocc(jorb))*BOM(iorb,jorb,ibas)*BOM(iorb,jorb,jbas)
                end do
            end do
        end do
        LI(ibas)=DI(ibas,ibas)
    end do
    DI=2*(DI+transpose(DI))
    do ibas=1,numrealatt !Diagonal terms are the sum of corresponding row or column
        DI(ibas,ibas)=0D0
        DI(ibas,ibas)=sum(DI(ibas,:))
    end do
else if (wfntype==2) then !ROHF
    DIa=0D0
    DIb=0D0
    do nmoclose=nmatsize,1,-1
        if (MOtype(nmoclose)==0) exit
    end do
    do ibas=1,numrealatt
        do jbas=ibas,numrealatt
            !Alpha
            do iorb=1,nmatsize !The number of close or alpha orbitals needed to be concerned
                occi=MOocc(iorb)
                if (MOtype(iorb)==0) occi=occi/2D0
                do jorb=1,nmatsize
                    occj=MOocc(jorb)
                    if (MOtype(jorb)==0) occj=occj/2D0
                    DIa(ibas,jbas)=DIa(ibas,jbas)+dsqrt(occi*occj)*BOM(iorb,jorb,ibas)*BOM(iorb,jorb,jbas)
                end do
            end do
            !Beta
            do iorb=1,nmoclose !The number of close orbitals needed to be concerned
                do jorb=1,nmoclose
                    DIb(ibas,jbas)=DIb(ibas,jbas)+dsqrt(MOocc(iorb)/2D0*MOocc(jorb)/2D0)*BOM(iorb,jorb,ibas)*BOM(iorb,jorb,jbas)
                end do
            end do
        end do
        LIa(ibas)=DIa(ibas,ibas)
        LIb(ibas)=DIb(ibas,ibas)
    end do
    DIa=2*(DIa+transpose(DIa))
    DIb=2*(DIb+transpose(DIb))
    do ibas=1,numrealatt !Diagonal terms are the sum of corresponding row or column
        DIa(ibas,ibas)=0D0
        DIb(ibas,ibas)=0D0
        DIa(ibas,ibas)=sum(DIa(ibas,:))
        DIb(ibas,ibas)=sum(DIb(ibas,:))
    end do
    !Combine alpha and Beta to total
    DI=DIa+DIb
    LI=LIa+LIb
!UHF,U-post-HF   DI(A,B)=2¡Æ[i,j]dsqrt(n_i*n_j)*S_i,j_A * S_i,j_B   where i and j are spin orbitals
else if (wfntype==1.or.wfntype==4) then
    !Alpha
    DIa=0D0
    do ibas=1,numrealatt
        do jbas=ibas,numrealatt
            do iorb=1,nmatsizea
                do jorb=1,nmatsizea
                    DIa(ibas,jbas)=DIa(ibas,jbas)+dsqrt(MOocc(iorb)*MOocc(jorb))*BOMa(iorb,jorb,ibas)*BOMa(iorb,jorb,jbas)
                end do
            end do
        end do
        LIa(ibas)=DIa(ibas,ibas)
    end do
    DIa=2*(DIa+transpose(DIa))
    !Beta
    if (nmatsizeb>0) then
        DIb=0D0
        MOinit=iendalpha+1 !Index range of beta orbitals
        MOend=iendalpha+nmatsizeb
        do ibas=1,numrealatt
            do jbas=ibas,numrealatt
                do iorb=MOinit,MOend
                    iorbtmp=iorb-iendalpha
                    do jorb=MOinit,MOend
                        jorbtmp=jorb-iendalpha
                        DIb(ibas,jbas)=DIb(ibas,jbas)+dsqrt(MOocc(iorb)*MOocc(jorb))*BOMb(iorbtmp,jorbtmp,ibas)*BOMb(iorbtmp,jorbtmp,jbas)
                    end do
                end do
            end do
            LIb(ibas)=DIb(ibas,ibas)
        end do
        DIb=2*(DIb+transpose(DIb))
    end if
    do ibas=1,numrealatt !Diagonal terms are the sum of corresponding row or column
        DIa(ibas,ibas)=0D0
        DIb(ibas,ibas)=0D0
        DIa(ibas,ibas)=sum(DIa(ibas,:))
        DIb(ibas,ibas)=sum(DIb(ibas,:))
    end do
    !Combine alpha and Beta to total
    DI=DIa+DIb
    LI=LIa+LIb
end if

call walltime(nwalltime2)
write(*,"(' Calculation took up',i8,' seconds wall clock time',/)")  nwalltime2-nwalltime1

!Output LI and DI
ioutid=6
100 write(ioutid,"(a,/)") " Note: Diagonal terms are the sum of corresponding row or column elements"
if (wfntype==1.or.wfntype==2.or.wfntype==4) then !UHF,ROHF,U-post-HF, output each spin component first
    !Alpha
    call showmatgau(DIa,"Delocalization index matrix for alpha spin",0,"f14.8",ioutid)
    write(ioutid,*)
    write(ioutid,*) "Localization index for alpha electron:"
    do ibas=1,numrealatt
        write(ioutid,"(i6,':',f7.3)",advance='no') ibas,LIa(ibas)
        if (mod(ibas,5)==0) write(ioutid,*)
    end do
    write(ioutid,*)
    write(ioutid,*)
    !Beta
    call showmatgau(DIb,"Delocalization index matrix for beta spin",0,"f14.8",ioutid)
    write(ioutid,*)
    write(ioutid,*) "Localization index for beta spin:"
    do ibas=1,numrealatt
        write(ioutid,"(i6,':',f7.3)",advance='no') ibas,LIb(ibas)
        if (mod(ibas,5)==0) write(ioutid,*)
    end do
    write(ioutid,*)
    write(ioutid,*)
end if
!Alpha+Beta
call showmatgau(DI,"Total delocalization index matrix",0,"f14.8",ioutid)
write(ioutid,*)
write(ioutid,*) "Total localization index:"
do ibas=1,numrealatt
    write(ioutid,"(i6,':',f7.3)",advance='no') ibas,LI(ibas)
    if (mod(ibas,5)==0) write(ioutid,*)
end do
if (ioutid==10) then
    write(*,*) "Done!"
    close(10)
    return
end if
write(*,*)
write(*,*)
write(*,*) "If also output the result to LIDI.txt in current folder? y/n"
read(*,*) selectyn
if (selectyn=='y'.or.selectyn=='Y') then
    open(10,file="LIDI.txt",status="replace")
    ioutid=10
    goto 100
end if
end subroutine



!!------- Integrate a real space function in the basins with multi-level refinement, indenspensible for Laplacian
!DEPRECATED, since mixed type of grids perform better and faster
subroutine integratebasinrefine
use defvar
use util
use function
use basinintmod
implicit real*8(a-h,o-z)
character c200tmp*200
real*8 intval(-1:numatt),basinvol(-1:numatt),intvalpriv(-1:numatt),basinvolpriv(-1:numatt)
integer walltime1,walltime2
real*8 :: critlevel1=0.1D0,critlevel2=0.5D0,critlevel3=1D0
integer :: nrefine1=1,nrefine2=2,nrefine3=5,nrefine4=7
if (ifuncbasin/=1) then
    write(*,*) "Error: This function is only applicable to AIM basins!"
    return
end if

do while(.true.)
    write(*,*) "-2 Return"
    write(*,*) "-1 Print and set parameters for multi-level refinement"
    call selfunc_interface(ifuncint)
    if (ifuncint==-2) then
        return
    else if (ifuncint==-1) then
        write(*,"(a)") " Note: A number n means each grid in corresponding value range will be transformed &
        as n^3 grids around it during the integration to gain a higher integration accuracy. n=1 means the grids will remain unchanged. &
        The ""value range"" referred here is the value range of the function used to generate basins"
        write(*,*) "The value range and the times of refinement:"
        write(*,"(' Smaller than ',f10.5,' :',i4)") critlevel1,nrefine1
        write(*,"(' Between',f10.5,' and ',f10.5,' :',i4)") critlevel1,critlevel2,nrefine2
        write(*,"(' Between',f10.5,' and ',f10.5,' :',i4)") critlevel2,critlevel3,nrefine3
        write(*,"(' Larger than  ',f10.5,' :',i4)") critlevel3,nrefine4
        write(*,*)
        write(*,*) "Please input three thresholds to define the four ranges, e.g. 0.1,0.5,1"
        write(*,*) "Note: Press ENTER button can retain current values unchanged"
        read(*,"(a)") c200tmp
        if (c200tmp/=' ') read(c200tmp,*) critlevel1,critlevel2,critlevel3
        write(*,*) "Input the times of refinement for the grids in the four ranges, e.g. 1,2,5,7"
        write(*,*) "Note: Press ENTER button can retain current values unchanged"
        read(*,"(a)") c200tmp
        if (c200tmp/=' ') read(c200tmp,*) nrefine1,nrefine2,nrefine3,nrefine4
        write(*,*) "Done!"
        write(*,*)
    end if
end do

call walltime(walltime1)
CALL CPU_TIME(time_begin)
write(*,*) "Integrating, please wait..."
intval=0D0
basinvol=0D0
ifinish=0
nthreads=getNThreads()
!$OMP PARALLEL private(ix,iy,iz,ixref,iyref,izref,ndiv,irealatt,rnowx,rnowy,rnowz,rnowxtmp,rnowytmp,rnowztmp,orgxref,orgyref,orgzref,dxref,dyref,dzref,&
!$OMP tmpval,tmpvalrefine,intvalpriv,basinvolpriv,nrefine) shared(intval,basinvol,ifinish) NUM_THREADS(nthreads)
intvalpriv=0D0
basinvolpriv=0D0
!$OMP do schedule(DYNAMIC)
do iz=2,nz-1
    do iy=2,ny-1
        do ix=2,nx-1
            if (cubmat(ix,iy,iz)<critlevel1) then
                nrefine=nrefine1 !The number of point to represent each edge
            else if (cubmat(ix,iy,iz)<critlevel2) then
                nrefine=nrefine2
            else if (cubmat(ix,iy,iz)<critlevel3) then
                nrefine=nrefine3
            else
                nrefine=nrefine4
            end if
            ndiv=nrefine**3
            rnowx=orgx+(ix-1)*dx
            rnowy=orgy+(iy-1)*dy
            rnowz=orgz+(iz-1)*dz
            orgxref=rnowx-dx/2 !Take corner position as original point of microcycle
            orgyref=rnowy-dy/2
            orgzref=rnowz-dz/2
            dxref=dx/nrefine
            dyref=dy/nrefine
            dzref=dz/nrefine
            tmpval=0D0
            do ixref=1,nrefine
                do iyref=1,nrefine
                    do izref=1,nrefine
                        rnowxtmp=orgxref+(ixref-0.5D0)*dxref
                        rnowytmp=orgyref+(iyref-0.5D0)*dyref
                        rnowztmp=orgzref+(izref-0.5D0)*dzref
                        if (ifuncint==-1) then
                            tmpvalrefine=cubmattmp(ix,iy,iz)
                        else if (ifuncint==0) then
                            tmpvalrefine=cubmat(ix,iy,iz)
                        else
                            tmpvalrefine=calcfuncall(ifuncint,rnowxtmp,rnowytmp,rnowztmp)
                        end if
                        tmpval=tmpval+tmpvalrefine/ndiv
                    end do
                end do
            end do
            irealatt=gridbas(ix,iy,iz)
            intvalpriv(irealatt)=intvalpriv(irealatt)+tmpval
            basinvolpriv(irealatt)=basinvolpriv(irealatt)+1
        end do
    end do
    ifinish=ifinish+1
    write(*,"(' Finished:',i5,'/',i5)") ifinish,nz
end do
!$OMP end do
!$OMP CRITICAL
    intval=intval+intvalpriv
    basinvol=basinvol+basinvolpriv
!$OMP end CRITICAL
!$OMP END PARALLEL
dvol=dx*dy*dz
intval=intval*dvol
basinvol=basinvol*dvol !Basin volume
write(*,*) "  #Basin          Integral        Volume(a.u.^3)"
do irealatt=1,numrealatt
    write(*,"(i8,f22.10,f20.8)") irealatt,intval(irealatt),basinvol(irealatt)
end do
write(*,"(' Sum of above values:',f20.8)") sum(intval(1:numrealatt))
if (any(gridbas(2:nx-1,2:ny-1,2:nz-1)==0)) write(*,"(' Integral of unassigned grids:',f20.8)") intval(0)
if (any(gridbas(2:nx-1,2:ny-1,2:nz-1)==-1)) write(*,"(' Integral of the grids travelled to box boundary:',f20.8)") intval(-1)

CALL CPU_TIME(time_end)
call walltime(walltime2)
write(*,"(' Integrating basins took up CPU time',f12.2,'s, wall clock time',i10,'s')") time_end-time_begin,walltime2-walltime1
end subroutine



!!------- Integrate AIM basins using mixed atomic-center and uniform grids
! itype=1: Integrate specific real space function
! itype=2: Integrate specific real space function with exact refinement of basin boundary
! itype=3: Integrate specific real space function with approximate refinement of basin boundary by exact+linear interpolation
! itype=10: Produce electric multipole moments
!NNA and ECP are supported
!Notice that the grid generated must cover the whole space!
subroutine integratebasinmix(itype)
use defvar
use util
use function
use basinintmod
use topo
implicit real*8(a-h,o-z)
!p suffix means private variable for parallel mode
real*8 intval(-1:numrealatt,20),intvalp(-1:numrealatt,20) !Up to 20 functions can be evaluated and stored simultaneously
real*8 basinvol(-1:numrealatt),basinvolp(-1:numrealatt),basinvdwvol(-1:numrealatt),basinvdwvolp(-1:numrealatt) !vdW is used to obtain the basin volume enclosed by 0.001 isosurface of rho
real*8 trustrad(numrealatt),intbasinthread(numrealatt),intbasin(numrealatt)
real*8 dens,grad(3),hess(3,3),k1(3),k2(3),k3(3),k4(3),xarr(nx),yarr(ny),zarr(nz)
real*8,allocatable :: potx(:),poty(:),potz(:),potw(:)
real*8,allocatable :: rhogrid(:,:,:),rhograd2grid(:,:,:),prorhogrid(:,:,:) !rhogrid and rhograd2grid are used for shubin's 2nd project, prorhogrid for integrating deformation density
type(content),allocatable :: gridatt(:) !Record correspondence between attractor and grid
integer att2atm(numrealatt) !The attractor corresponds to which atom. If =0, means this is a NNA
real*8 eleint(-1:numrealatt),xint(-1:numrealatt),yint(-1:numrealatt),zint(-1:numrealatt),&
xxint(-1:numrealatt),yyint(-1:numrealatt),zzint(-1:numrealatt),xyint(-1:numrealatt),yzint(-1:numrealatt),xzint(-1:numrealatt)
real*8 eleintp(-1:numrealatt),xintp(-1:numrealatt),yintp(-1:numrealatt),zintp(-1:numrealatt),&
xxintp(-1:numrealatt),yyintp(-1:numrealatt),zzintp(-1:numrealatt),xyintp(-1:numrealatt),yzintp(-1:numrealatt),xzintp(-1:numrealatt)
integer walltime1,walltime2,radpotAIM,sphpotAIM
real*8 gridval(100000,3) !Used for shubin's 2nd project
character c80tmp*80,selectyn
nbeckeiter=8
if (ifuncbasin/=1) then
    write(*,"(a)") " Error: This function is only applicable to AIM basins! That means in option 1, you should select electron density to construct the basins."
    return
end if

if (ispecial==0) then
    if (itype==1.or.itype==2.or.itype==3) then
        write(*,*) "Please select integrand:"
        write(*,*) "-2 Return"
        write(*,*) "-1 Deformation density"
        call selfunc_interface(ifuncint)
        if (ifuncint==-2) then
            return
        else if (ifuncint==-1) then
            call setpromol
        end if
    end if
else if (ispecial==1) then
    continue !Don't let user to select integrand
else if (ispecial==2) then
    call setpromol !Don't let user to select integrand
    expcutoff=1 !Use full accuracy for shubin's 2nd project
end if

call walltime(walltime1)
CALL CPU_TIME(time_begin)
intval=0D0
basinvol=0D0
basinvdwvol=0D0
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
numcp=0

att2atm=0
!Determine trust radius and then integrate in the trust sphere
!We use CPpos to record the center of the trust sphere
write(*,*) "Integrating in trust sphere..."
do iatt=1,numrealatt !Cycle each attractors
    !Refine the crude position of attractor by exact newton method
    do iatm=1,ncenter
        disttest=dsqrt( (realattxyz(iatt,1)-a(iatm)%x)**2+(realattxyz(iatt,2)-a(iatm)%y)**2+(realattxyz(iatt,3)-a(iatm)%z)**2 )
        if (disttest<0.3D0) then
            att2atm(iatt)=iatm
            write(*,"(' Attractor',i6,' corresponds to atom',i6,' (',a,')')") iatt,iatm,a(iatm)%name
            numcpold=numcp
            call findcp(a(iatm)%x,a(iatm)%y,a(iatm)%z,1,0)
            if (numcp==numcpold) then
                write(*,*) "Note: Unable to locate exact CP position! Use nuclear position"
                numcp=numcp+1
                CPpos(1,numcp)=a(iatm)%x
                CPpos(2,numcp)=a(iatm)%y
                CPpos(3,numcp)=a(iatm)%z
            end if
!             write(*,"(' Coordinate after refinement:',3f16.8)") CPpos(:,numcp)
            exit
        end if
    end do
    if (att2atm(iatt)==0) then
        write(*,"(a,i6,a)") " Warning: Unable to determine the attractor",iatt," belongs to which atom!"
        write(*,"(a)") " If this is a non-nuclear attractor, simply press ENTER button to continue. If you used pseudopotential &
        and this attractor corresponds to the cluster of all maxima of its valence electron, then input the index of this atom (e.g. 9). &
        Else you should input q to return and regenerate basins with smaller grid spacing"
        read(*,"(a)") c80tmp
        if (c80tmp=='q') then
            return
        else if (c80tmp==" ") then
            numcpold=numcp
            call findcp(realattxyz(iatt,1),realattxyz(iatt,2),realattxyz(iatt,3),1,0)
            if (numcp==numcpold) then
                write(*,*) "Unable to locate exact CP position! Exit..."
                return
            end if
        else !ECP, input the corresponding atom by user
            read(c80tmp,*) iatmtmp
            att2atm(iatt)=iatmtmp
            numcp=numcp+1
            CPpos(1,numcp)=a(iatmtmp)%x
            CPpos(2,numcp)=a(iatmtmp)%y
            CPpos(3,numcp)=a(iatmtmp)%z
        end if
    end if
    
    !Determine trust radius and set integration points and weight
    radpotAIM=200
    parm=1
    isettrustrad=0
    nintgrid=0 !Then number of integration grids within trust radius
    if (allocated(gridatt)) deallocate(gridatt) !Used to record grids in trust sphere of this attractor
    allocate(gridatt(radpotAIM*500))
    do ish=1,radpotAIM !Cycle each radial shell. Radius distance is from near to far
        if (isettrustrad==1) exit !The trust radius has been finally determined in last shell cycle
        !Becke, namely the second-kind Gauss-Chebyshev
        itmp=radpotAIM+1-ish !Invert ish to make radr from near to far
        radx=cos(itmp*pi/(radpotAIM+1D0))
        radr=(1+radx)/(1-radx)*parm
        radw=2*pi/(radpotAIM+1)*parm**3 *(1+radx)**2.5D0/(1-radx)**3.5D0 *4*pi
!         !Handy, also known as Euler-Maclaurin. See: Murray, C. W.; Handy, N. C.; Laming, G. J. Mol Phys 1993, 78, 997
!         radx=dfloat(ish)/(radpotAIM+1D0)
!         radr=radx**2/(1-radx)**2*parm
!         radw=2*radx**5/dfloat(radpotAIM+1)/(1-radx)**7*parm**3 *4*pi
        
        !Set Lebedev grids according to shell radius
        !For more inner shell, the distribution is more akin to spherically symmetric, therefore lower number of grids could be used
        if (att2atm(iatt)==0) then !NNA
            sphpotAIM=302
        else
            radtmp=covr(a(att2atm(iatt))%index)
            if (radr<0.2D0*radtmp) then
                sphpotAIM=26
            else if (radr<0.5D0*radtmp) then
                sphpotAIM=74
            else if (radr<0.8D0*radtmp) then
                sphpotAIM=146
            else
                sphpotAIM=194
            end if
        end if
        if (allocated(potx)) deallocate(potx,poty,potz,potw)
        allocate(potx(sphpotAIM),poty(sphpotAIM),potz(sphpotAIM),potw(sphpotAIM))
        call Lebedevgen(sphpotAIM,potx,poty,potz,potw)
        !Combine radial point and weights with angular part, and make them centered at current attractor
        gridatt( nintgrid+1:nintgrid+sphpotAIM )%x=radr*potx+CPpos(1,numcp)
        gridatt( nintgrid+1:nintgrid+sphpotAIM )%y=radr*poty+CPpos(2,numcp)
        gridatt( nintgrid+1:nintgrid+sphpotAIM )%z=radr*potz+CPpos(3,numcp)
        gridatt( nintgrid+1:nintgrid+sphpotAIM )%value=radw*potw
        !Find out trust radius for present attractor
        !If in a shell, the angle between "linking line between nucleus and a shell point" and "gradient vector of this point" &
        !is larger than 45 degree, then this shell is trust radius
        angmax=0
        if (att2atm(iatt)==0) then
            radinit=0
        else
            radrinit=0.15D0
            if (a(att2atm(iatt))%index>2) radrinit=0.5D0
        end if
        if (isettrustrad==0.and.radr>radrinit) then
            do isphpt=1,sphpotAIM
                xtmp=gridatt(nintgrid+isphpt)%x
                ytmp=gridatt(nintgrid+isphpt)%y
                ztmp=gridatt(nintgrid+isphpt)%z
                call calchessmat_dens(1,xtmp,ytmp,ztmp,dens,grad,hess) !Only value and gradient
                dirx=CPpos(1,numcp)-xtmp
                diry=CPpos(2,numcp)-ytmp
                dirz=CPpos(3,numcp)-ztmp
                angtmp=vecang(dirx,diry,dirz,grad(1),grad(2),grad(3))
                if (angtmp>angmax) angmax=angtmp
                if (angtmp>45) then
                    write(*,"(' The trust radius of attractor',i6,' is',f10.3,' Bohr',/)") iatt,trustrad(iatt)
                    isettrustrad=1 !The radius of last shell should be the final trust radius. Now exit
                    exit
                end if
            end do
            if (isettrustrad==0) trustrad(iatt)=radr !Passed this shell and temporarily set the radius as trust radius. Continue to enlarge the trust radius, until reached angmax>45 degree
        end if
        nintgrid=nintgrid+sphpotAIM
    end do
    if (isettrustrad==0) then !Trust radius was not set after run over all shells
        trustrad(iatt)=1000 !Infinite, for isolated atom
        if (ispecial==2) trustrad(iatt)=20 !For Shubin's 2nd project, should not be as large as 1000, because for a point very far from nucleus the relative entropy cannot be evaluated
        write(*,"(' The trust radius of attractor',i6,' is',f10.3,' Bohr',/)") iatt,trustrad(iatt)
    end if
    
    !Use DFT integration algorithm to integrate the region inside trust radius
nthreads=getNThreads()
!$OMP PARALLEL private(ipt,ptx,pty,ptz,rx,ry,rz,dist,tmps,iter,switchwei,intvalp,&
!$OMP eleintp,xintp,yintp,zintp,xxintp,yyintp,zzintp,xyintp,yzintp,xzintp,tmpval,tmpval2,tmpval3) shared(intval,gridval) NUM_THREADS(nthreads)
    intvalp=0D0
    eleintp=0D0
    xintp=0D0
    yintp=0D0
    zintp=0D0
    xxintp=0D0
    yyintp=0D0
    zzintp=0D0
    xyintp=0D0
    yzintp=0D0
    xzintp=0D0
!$OMP do schedule(DYNAMIC)
    do ipt=1,nintgrid
        ptx=gridatt(ipt)%x
        pty=gridatt(ipt)%y
        ptz=gridatt(ipt)%z
        rx=ptx-CPpos(1,numcp) !The relative distance between current point to corresponding attractor
        ry=pty-CPpos(2,numcp)
        rz=ptz-CPpos(3,numcp)
        !Calculate switching function
        dist=dsqrt(rx*rx+ry*ry+rz*rz)
        tmps=dist-trustrad(iatt)
        if (tmps>1) then
            switchwei=0
        else if (tmps<-1) then
            switchwei=1
        else
            do iter=1,nbeckeiter
                tmps=1.5D0*(tmps)-0.5D0*(tmps)**3
            end do
            switchwei=0.5D0*(1-tmps)
        end if
        gridval(ipt,3)=switchwei
!         if (dist>trustrad(iatt)) cycle !Discrete separation between atomic-center and uniform integration
        if (switchwei<1D-7) cycle !For saving computational time
        
        if (itype==1.or.itype==2.or.itype==3) then !Integrate a function
            if (ifuncint==-1) then !Deformation density, store molecular density of present center temporarily
                gridval(ipt,1)=fdens(ptx,pty,ptz)
            else if (ispecial==0) then !Normal case
                tmpval=calcfuncall(ifuncint,ptx,pty,ptz)
                intvalp(iatt,1)=intvalp(iatt,1)+gridatt(ipt)%value*tmpval*switchwei
            else if (ispecial==1) then
                tmpval=infoentro(2,ptx,pty,ptz) !Shannon entropy density, see JCP,126,191107 for example
                tmpval2=Fisherinfo(1,ptx,pty,ptz) !Fisher information density, see JCP,126,191107 for example
                tmpval3=weizsacker(ptx,pty,ptz) !Steric energy
                tmpval4=fdens(ptx,pty,ptz) !Electron density
                intvalp(iatt,1)=intvalp(iatt,1)+gridatt(ipt)%value*tmpval*switchwei
                intvalp(iatt,2)=intvalp(iatt,2)+gridatt(ipt)%value*tmpval2*switchwei
                intvalp(iatt,3)=intvalp(iatt,3)+gridatt(ipt)%value*tmpval3*switchwei
                intvalp(iatt,4)=intvalp(iatt,4)+gridatt(ipt)%value*tmpval4*switchwei
            else if (ispecial==2) then !Temporarily for Shubin's 2nd project
                gridval(ipt,1)=fdens(ptx,pty,ptz)
                gridval(ipt,2)=fgrad(ptx,pty,ptz,'t')**2
            end if
        else if (itype==10) then !Calculate multipole moment
            tmpval=gridatt(ipt)%value*fdens(ptx,pty,ptz)*switchwei
            eleintp(iatt)=eleintp(iatt)+tmpval
            xintp(iatt)=xintp(iatt)+rx*tmpval
            yintp(iatt)=yintp(iatt)+ry*tmpval
            zintp(iatt)=zintp(iatt)+rz*tmpval
            xxintp(iatt)=xxintp(iatt)+rx*rx*tmpval
            yyintp(iatt)=yyintp(iatt)+ry*ry*tmpval
            zzintp(iatt)=zzintp(iatt)+rz*rz*tmpval
            xyintp(iatt)=xyintp(iatt)+rx*ry*tmpval
            yzintp(iatt)=yzintp(iatt)+ry*rz*tmpval
            xzintp(iatt)=xzintp(iatt)+rx*rz*tmpval
        end if
    end do
!$OMP end do
!$OMP CRITICAL
    if (itype==1.or.itype==2.or.itype==3) then
        intval=intval+intvalp
    else if (itype==10) then
        eleint=eleint+eleintp
        xint=xint+xintp
        yint=yint+yintp
        zint=zint+zintp
        xxint=xxint+xxintp
        yyint=yyint+yyintp
        zzint=zzint+zzintp
        xyint=xyint+xyintp
        yzint=yzint+yzintp
        xzint=xzint+xzintp
    end if
!$OMP end CRITICAL
!$OMP END PARALLEL
    
    !Some special cases:
    if (ispecial==2) then !Shubin's 2nd project, integrate relative Shannon and Fisher entropy. density and gradient^2 have been stored in gridval(1/2)
        call dealloall
        write(*,"(' Loading ',a,/)") trim(custommapname(att2atm(iatt)))
        call readwfn(custommapname(att2atm(iatt)),1)
        do ipt=1,nintgrid
            switchwei=gridval(ipt,3)
            if (switchwei<1D-7) cycle
            ptx=gridatt(ipt)%x
            pty=gridatt(ipt)%y
            ptz=gridatt(ipt)%z
            prodens=fdens(ptx,pty,ptz) !rho0_A
            prodensgrad2=fgrad(ptx,pty,ptz,'t')**2
            tmpval=gridval(ipt,1)*log(gridval(ipt,1)/prodens)
            intval(iatt,1)=intval(iatt,1)+gridatt(ipt)%value*tmpval*switchwei !Relative Shannon, integrate rho*log(rho/rho0_A)]
            tmpval=gridval(ipt,2)/gridval(ipt,1)-prodensgrad2/prodens
            intval(iatt,2)=intval(iatt,2)+gridatt(ipt)%value*tmpval*switchwei !Relative Fisher, integrate grad2rho/rho-grad2rho0_A/rho0_A
        end do
        call dealloall
        call readinfile(firstfilename,1) !Retrieve to first loaded file(whole molecule) to calc real rho again
    else if (ifuncint==-1) then !Integrate deformation density
        gridval(:,2)=0D0
        do iatm=1,ncenter_org !Cycle each atom to calculate deformation density at all integration grid
            call dealloall
            call readwfn(custommapname(iatm),1)
            do ipt=1,nintgrid
                switchwei=gridval(ipt,3)
                if (switchwei<1D-7) cycle
                gridval(ipt,2)=gridval(ipt,2)+fdens(gridatt(ipt)%x,gridatt(ipt)%y,gridatt(ipt)%z)
            end do
        end do
        do ipt=1,nintgrid !Now gridval(ipt,1/2/3) records actual, promolecular density and weight at ipt
            defdens=gridval(ipt,1)-gridval(ipt,2)
            intval(iatt,1)=intval(iatt,1)+gridatt(ipt)%value*defdens*gridval(ipt,3)
        end do
        call dealloall
        call readinfile(firstfilename,1) !Retrieve to first loaded file(whole molecule) to calc real rho again
    end if
end do !End cycle attractors

if (itype==1.or.itype==2.or.itype==3) then
    write(*,*) "Integration result inside trust spheres"
    if (ispecial/=2) then
        write(*,*) "  #Sphere       Integral(a.u.)"
        do iatt=1,numrealatt
            write(*,"(i8,f22.10)") iatt,intval(iatt,1)
        end do
        write(*,"(' Sum of above values:',f20.8)") sum(intval(1:numrealatt,1))
    else if (ispecial==2) then !Shubin's 2nd project
        write(*,*) "  #Sphere     Relat_Shannon        Relat_Fisher"
        do iatt=1,numrealatt
            write(*,"(i8,2f20.10)") iatt,intval(iatt,1:2)
        end do
        write(*,"(' Sum of relat_Shannon:',f20.8)") sum(intval(1:numrealatt,1))
        write(*,"(' Sum of relat_Fisher: ',f20.8)") sum(intval(1:numrealatt,2))
    end if
end if

!Set coordinate of uniform grids
dvol=dx*dy*dz
do ix=1,nx
    xarr(ix)=orgx+(ix-1)*dx
end do
do iy=1,ny
    yarr(iy)=orgy+(iy-1)*dy
end do
do iz=1,nz
    zarr(iz)=orgz+(iz-1)*dz
end do

!--------- Integrating uniform grids, basin boundary grids will be calculated in the later stage
!NOTE: For shubin's 2nd project or deformation density, do not integrate here. At next stage boundary grids will be updated, and then at the next stage,&
!they will be integrated by a special module. Because at each grid if we reload atomic .wfn file will be too time consuming
if (ispecial==2.or.ifuncint==-1) goto 10 
write(*,*)
write(*,*) "Integrating uniform grids..."
ifinish=0
nthreads=getNThreads()
!$OMP PARALLEL private(ix,iy,iz,iatt,icp,rnowx,rnowy,rnowz,rx,ry,rz,tmpval,tmpval2,tmpval3,intvalp,basinvolp,basinvdwvolp,dist,tmps,iter,switchwei,&
!$OMP eleintp,xintp,yintp,zintp,xxintp,yyintp,zzintp,xyintp,xzintp,yzintp) shared(intval,basinvol,basinvdwvol,ifinish) NUM_THREADS(nthreads)
intvalp=0D0
basinvolp=0D0
basinvdwvolp=0D0
eleintp=0D0
xintp=0D0
yintp=0D0
zintp=0D0
xxintp=0D0
yyintp=0D0
zzintp=0D0
xyintp=0D0
yzintp=0D0
xzintp=0D0
!$OMP do schedule(DYNAMIC)
do iz=2,nz-1
    rnowz=zarr(iz)
    do iy=2,ny-1
        rnowy=yarr(iy)
        do ix=2,nx-1
            rnowx=xarr(ix)
            if ((itype==2.or.itype==3).and.interbasgrid(ix,iy,iz) .eqv. .true.) cycle !If refine boundary grid at next stage, we don't calculate them at present stage
!             if (iz==nint(nz/2D0)) write(1,"('C',3f14.8)") rnowx*b2a,rnowy*b2a,rnowz*b2a !Examine grid distribution
            iatt=gridbas(ix,iy,iz)
!             do icp=1,numcp
!                 dist=dsqrt( (rnowx-CPpos(1,icp))**2+(rnowy-CPpos(2,icp))**2+(rnowz-CPpos(3,icp))**2 )
!                 if (disttest<trustrad(icp)) cycle cycix !The function inside trust radius is integrated by DFT integration
!             end do
            rx=rnowx-CPpos(1,iatt) !The relative distance between current point to corresponding attractor
            ry=rnowy-CPpos(2,iatt)
            rz=rnowz-CPpos(3,iatt)
            !Calculate switching function at current grid
            dist=dsqrt(rx*rx+ry*ry+rz*rz)
            tmps=dist-trustrad(iatt)
            if (tmps>1) then
                switchwei=0
            else if (tmps<-1) then
                switchwei=1
            else
                do iter=1,nbeckeiter
                    tmps=1.5D0*(tmps)-0.5D0*(tmps)**3
                end do
                switchwei=0.5D0*(1-tmps)
            end if
            switchwei=1-switchwei
            basinvolp(iatt)=basinvolp(iatt)+1 !Calculate basin volume
            if (cubmat(ix,iy,iz)>0.001D0) basinvdwvolp(iatt)=basinvdwvolp(iatt)+1
            if (switchwei<1D-7) cycle !For saving time
            if (itype==1.or.itype==2.or.itype==3) then
                if (ispecial==0) then
                    if (ifuncint==1) then !Electron density on each grid has already been calculated
                        tmpval=cubmat(ix,iy,iz)
                    else
                        tmpval=calcfuncall(ifuncint,rnowx,rnowy,rnowz)
                    end if
                    intvalp(iatt,1)=intvalp(iatt,1)+tmpval*switchwei
                else if (ispecial==1) then !For RCY
                    tmpval=infoentro(2,rnowx,rnowy,rnowz) !Shannon entropy density, see JCP,126,191107 for example
                    tmpval2=Fisherinfo(1,rnowx,rnowy,rnowz) !Fisher information density, see JCP,126,191107 for example
                    tmpval3=weizsacker(rnowx,rnowy,rnowz) !Steric energy
                    tmpval4=fdens(rnowx,rnowy,rnowz) !Electron density
                    intvalp(iatt,1)=intvalp(iatt,1)+tmpval*switchwei
                    intvalp(iatt,2)=intvalp(iatt,2)+tmpval2*switchwei
                    intvalp(iatt,3)=intvalp(iatt,3)+tmpval3*switchwei
                    intvalp(iatt,4)=intvalp(iatt,4)+tmpval4*switchwei
                end if
            else if (itype==10) then
                tmpval=cubmat(ix,iy,iz)*switchwei
                eleintp(iatt)=eleintp(iatt)+tmpval
                xintp(iatt)=xintp(iatt)+rx*tmpval
                yintp(iatt)=yintp(iatt)+ry*tmpval
                zintp(iatt)=zintp(iatt)+rz*tmpval
                xxintp(iatt)=xxintp(iatt)+rx*rx*tmpval
                yyintp(iatt)=yyintp(iatt)+ry*ry*tmpval
                zzintp(iatt)=zzintp(iatt)+rz*rz*tmpval
                xyintp(iatt)=xyintp(iatt)+rx*ry*tmpval
                yzintp(iatt)=yzintp(iatt)+ry*rz*tmpval
                xzintp(iatt)=xzintp(iatt)+rx*rz*tmpval
            end if
        end do
    end do
    ifinish=ifinish+1
    if (ifuncint/=1) write(*,"(' Integrating uniform grids, finished:',i5,'/',i5)") ifinish,nz
end do
!$OMP end do
!$OMP CRITICAL
if (itype==1.or.itype==2.or.itype==3) then
    intval=intval+intvalp*dvol
    basinvol=basinvol+basinvolp*dvol
    basinvdwvol=basinvdwvol+basinvdwvolp*dvol
else if (itype==10) then
    eleint=eleint+eleintp*dvol
    xint=xint+xintp*dvol
    yint=yint+yintp*dvol
    zint=zint+zintp*dvol
    xxint=xxint+xxintp*dvol
    yyint=yyint+yyintp*dvol
    zzint=zzint+zzintp*dvol
    xyint=xyint+xyintp*dvol
    yzint=yzint+yzintp*dvol
    xzint=xzint+xzintp*dvol
end if
!$OMP end CRITICAL
!$OMP END PARALLEL

10 continue
!---- Exact refinement with/without multi-level splitting of boundary grids
if (itype==2.or.itype==3) then
    if (itype==3) then !Calculate grid data of gradient of electron density used to linear interpolation to obtain the value at any point
        write(*,*)
        call gengradmat
    end if
    write(*,*)
    nrk4lim=100
    nrk4gradswitch=40
    hsizeinit=0.25D0
    ifinish=0
nthreads=getNThreads()
!$OMP PARALLEL private(ix,iy,iz,iatt,rnowx,rnowy,rnowz,rx,ry,rz,tmpval,tmpval2,tmpval3,intvalp,basinvolp,basinvdwvolp,dist,tmps,iter,switchwei,&
!$OMP rnowxtmp,rnowytmp,rnowztmp,orgxref,orgyref,orgzref,dxref,dyref,dzref,ixref,iyref,izref,nrefine,ndiv,&
!$OMP k1,k2,k3,k4,dens,denshold,grad,hess,iattref,xtmp,ytmp,ztmp,irk4,hsize,ixtest,iytest,iztest,tmpdist) shared(intval,basinvol,basinvdwvol,ifinish) NUM_THREADS(nthreads)
    intvalp=0D0
    basinvolp=0D0
    basinvdwvolp=0D0
!$OMP do schedule(DYNAMIC)
    do iz=2,nz-1
        rnowz=zarr(iz)
        do iy=2,ny-1
            rnowy=yarr(iy)
            do ix=2,nx-1
                rnowx=xarr(ix)
                if (interbasgrid(ix,iy,iz) .eqv. .false.) cycle
!                 if (cubmat(ix,iy,iz)>0.001D0) then
!                     nrefine=2 !3 is the best
!                 else if (cubmat(ix,iy,iz)>0.001D0) then !0.0001 is the best
!                     nrefine=1
!                 else
!                     nrefine=1
!                 end if
                 nrefine=1
                ndiv=nrefine**3
                orgxref=rnowx-dx/2 !Take corner position as original point of microcycle
                orgyref=rnowy-dy/2
                orgzref=rnowz-dz/2
                dxref=dx/nrefine
                dyref=dy/nrefine
                dzref=dz/nrefine
                do ixref=1,nrefine
                    do iyref=1,nrefine
                        do izref=1,nrefine
                            rnowxtmp=orgxref+(ixref-0.5D0)*dxref !Coordinate of current refined grid
                            rnowytmp=orgyref+(iyref-0.5D0)*dyref
                            rnowztmp=orgzref+(izref-0.5D0)*dzref
                            if (cubmat(ix,iy,iz)<=0.001D0) then !Only refine the boundary inside vdW surface
                                iattref=gridbas(ix,iy,iz)
                            else
                                xtmp=rnowxtmp !This point will continuously move in the iteration
                                ytmp=rnowytmp
                                ztmp=rnowztmp
                                hsize=hsizeinit
                                densold=0D0
                                !** Tracing steepest ascent trajectory using 4-order Runge-Kutta (RK4)
        cycrk4:                    do irk4=1,nrk4lim
                                    !For full accuracy refinement, or the first step, or when interpolation gradient works worse,&
                                    !namely has not converge until nrk4gradswitch, use exactly evaluated gradient
                                    if (itype==2.or.irk4==1.or.irk4==2.or.irk4>nrk4gradswitch) then 
                                        if (itype==3.and.irk4==nrk4gradswitch+1) then !Interpolated gradient doesn't work well, switch to full accuracy, reset the coordinate
                                            xtmp=rnowxtmp
                                            ytmp=rnowytmp
                                            ztmp=rnowztmp
                                            hsize=hsizeinit
                                        end if
                                        call calchessmat_dens(1,xtmp,ytmp,ztmp,dens,grad,hess) !Only value and gradient
                                        if (dens<densold-1D-10) then
                                            hsize=hsize*0.75D0 !Reduce step size if density decrease
                                        else if (dens>densold+1D-10) then
                                            hsize=hsizeinit !Recover to initial step size
                                        end if
                                        denshold=dens
                                        k1=grad/dsqrt(sum(grad**2))
                                        call calchessmat_dens(1,xtmp+hsize/2*k1(1),ytmp+hsize/2*k1(2),ztmp+hsize/2*k1(3),dens,grad,hess) !Only value and gradient
                                        k2=grad/dsqrt(sum(grad**2))
                                        call calchessmat_dens(1,xtmp+hsize/2*k2(1),ytmp+hsize/2*k2(2),ztmp+hsize/2*k2(3),dens,grad,hess) !Only value and gradient
                                        k3=grad/dsqrt(sum(grad**2))
                                        call calchessmat_dens(1,xtmp+hsize*k3(1),ytmp+hsize*k3(2),ztmp+hsize*k3(3),dens,grad,hess) !Only value and gradient
                                        k4=grad/dsqrt(sum(grad**2))
                                    else !Using the gradients evaluated by trilinear interpolation from pre-calculated grid data to save computational time
                                        call linintp3dvec(xtmp,ytmp,ztmp,grad) !Only value and gradient
!                                         if (dens<densold-1D-10) then
!                                             hsize=hsize*0.75D0
!                                         else if (dens>densold+1D-10) then
!                                             hsize=hsizeinit !Recover to initial step size
!                                         end if
!                                         denshold=dens
                                        k1=grad/dsqrt(sum(grad**2))
                                        call linintp3dvec(xtmp+hsize/2*k1(1),ytmp+hsize/2*k1(2),ztmp+hsize/2*k1(3),grad)
                                        k2=grad/dsqrt(sum(grad**2))
                                        call linintp3dvec(xtmp+hsize/2*k2(1),ytmp+hsize/2*k2(2),ztmp+hsize/2*k2(3),grad)
                                        k3=grad/dsqrt(sum(grad**2))
                                        call linintp3dvec(xtmp+hsize*k3(1),ytmp+hsize*k3(2),ztmp+hsize*k3(3),grad)
                                        k4=grad/dsqrt(sum(grad**2))
                                    end if
                                    xtmp=xtmp+hsize/6*(k1(1)+2*k2(1)+2*k3(1)+k4(1)) !Update current coordinate
                                    ytmp=ytmp+hsize/6*(k1(2)+2*k2(2)+2*k3(2)+k4(2))
                                    ztmp=ztmp+hsize/6*(k1(3)+2*k2(3)+2*k3(3)+k4(3))
                                    !Check if current position has entered trust radius of an attractor
                                    do iatttmp=1,numrealatt
                                        dist=dsqrt( (xtmp-CPpos(1,iatttmp))**2+(ytmp-CPpos(2,iatttmp))**2+(ztmp-CPpos(3,iatttmp))**2 )
                                        if (dist<trustrad(iatttmp)) then
                                            iattref=iatttmp
                                            exit cycrk4
                                        end if
                                    end do
                                    !Check if the closest grid and its 26 neighbours have the same attribution, if yes, employ its attribution then exit
                                    do ixtest=2,nx-1
                                        tmpdist=abs(xtmp-xarr(ixtest))
                                        if (tmpdist<dx/2D0) exit
                                    end do
                                    do iytest=2,ny-1
                                        tmpdist=abs(ytmp-yarr(iytest))
                                        if (tmpdist<dy/2D0) exit
                                    end do
                                    do iztest=2,nz-1
                                        tmpdist=abs(ztmp-zarr(iztest))
                                        if (tmpdist<dz/2D0) exit
                                    end do
                                    iattref=gridbas(ixtest,iytest,iztest)
                                    do imove=1,26
                                        if ( gridbas(ixtest+vec26x(imove),iytest+vec26y(imove),iztest+vec26z(imove))/=iattref ) exit
                                    end do
                                    if (imove==27) exit !Successfully passed neighbour test
                                end do cycrk4
                                if (irk4==nrk4lim+1) then !Didn't enter trust radius or didn't approach a grid who and whose neighbour have the same attribution
                                    write(*,*) "Warning: Exceeded the step limit of steepest ascent process!"
                                    iattref=gridbas(ix,iy,iz) !Use its original attribution
                                end if
!                                 write(*,*) irk4
                            end if
                            gridbas(ix,iy,iz)=iattref !Update attribution of boundary grids
                            !Calculate switching function at current grid
                            rx=rnowxtmp-CPpos(1,iattref) !The relative distance between current point to corresponding attractor
                            ry=rnowytmp-CPpos(2,iattref)
                            rz=rnowztmp-CPpos(3,iattref)
                            dist=dsqrt(rx*rx+ry*ry+rz*rz)
                            tmps=dist-trustrad(iattref)
                            if (tmps>1) then
                                switchwei=0
                            else if (tmps<-1) then
                                switchwei=1
                            else
                                do iter=1,nbeckeiter
                                    tmps=1.5D0*(tmps)-0.5D0*(tmps)**3
                                end do
                                switchwei=0.5D0*(1-tmps)
                            end if
                            switchwei=1-switchwei
                            basinvolp(iattref)=basinvolp(iattref)+1D0/ndiv !Calculate boundary basin volume
                            if (cubmat(ix,iy,iz)>0.001D0) basinvdwvolp(iattref)=basinvdwvolp(iattref)+1D0/ndiv
                            if (ispecial==2.or.ifuncint==-1) then
                                continue !Don't calculate function value, but only update attribution of boundary grids at this stage
                            else if (ispecial==0) then
                                tmpval=calcfuncall(ifuncint,rnowxtmp,rnowytmp,rnowztmp)
                                intvalp(iattref,1)=intvalp(iattref,1)+tmpval*switchwei/ndiv
                            else if (ispecial==1) then
                                tmpval=infoentro(2,rnowxtmp,rnowytmp,rnowztmp) !Shannon entropy density, see JCP,126,191107 for example
                                tmpval2=Fisherinfo(1,rnowxtmp,rnowytmp,rnowztmp) !Fisher information density, see JCP,126,191107 for example
                                tmpval3=weizsacker(rnowxtmp,rnowytmp,rnowztmp) !Steric energy
                                tmpval4=fdens(rnowxtmp,rnowytmp,rnowztmp) !Electron density
                                intvalp(iattref,1)=intvalp(iattref,1)+tmpval*switchwei/ndiv
                                intvalp(iattref,2)=intvalp(iattref,2)+tmpval2*switchwei/ndiv
                                intvalp(iattref,3)=intvalp(iattref,3)+tmpval3*switchwei/ndiv
                                intvalp(iattref,4)=intvalp(iattref,4)+tmpval4*switchwei/ndiv
                            end if
                        end do !End refine grid
                    end do
                end do
                
            end do !End cycle ix grid
        end do
        ifinish=ifinish+1
        write(*,"(' Integrating grids at basin boundary, finished:',i5,'/',i5)") ifinish,nz
    end do
!$OMP end do
!$OMP CRITICAL
    intval=intval+intvalp*dvol
    basinvol=basinvol+basinvolp*dvol
    basinvdwvol=basinvdwvol+basinvdwvolp*dvol
!$OMP end CRITICAL
!$OMP END PARALLEL
    call detectinterbasgrd(6)
    write(*,*) "Basin boundary has been updated"
    numinterbas=count(interbasgrid .eqv. .true.)
    write(*,"(' The number of interbasin grids:',i12)") numinterbas    
end if

!Below are special modules for the cases when density of atom in free-state are involved
if (ispecial==2) then !Shubin's 2nd project, integrate relative Shannon and Fisher entropy
    allocate(rhogrid(nx,ny,nz),rhograd2grid(nx,ny,nz))
    ifinish=0
    write(*,*)
    write(*,*) "Calculating electron density and its gradient for actual system at each grid"
nthreads=getNThreads()
!$OMP PARALLEL DO SHARED(rhogrid,rhograd2grid,ifinish) PRIVATE(ix,iy,iz,ptx,pty,ptz) schedule(dynamic) NUM_THREADS(nthreads)
    do iz=2,nz-1
        ptz=zarr(iz)
        do iy=2,ny-1
            pty=yarr(iy)
            do ix=2,nx-1
                ptx=xarr(ix)
                rhogrid(ix,iy,iz)=fdens(ptx,pty,ptz)
                rhograd2grid(ix,iy,iz)=fgrad(ptx,pty,ptz,'t')**2
            end do
        end do
!$OMP CRITICAL
        ifinish=ifinish+1
        write(*,"(' Finished',i6,' /',i6)") ifinish,nz-1
!$OMP end CRITICAL
    end do
!$OMP end PARALLEL DO
    write(*,*)
    write(*,*) "Calculating electron density and its gradient for free-state atom at each grid"
    do iatt=1,numrealatt !Cycle each attractors
        write(*,"(' Processing ',a)") trim(custommapname(att2atm(iatt)))
        call dealloall
        call readwfn(custommapname(att2atm(iatt)),1)
nthreads=getNThreads()
!$OMP PARALLEL private(intvalp,ix,iy,iz,ptx,pty,ptz,rx,ry,rz,dist,tmps,switchwei,prodens,prodensgrad2,tmpval1,tmpval2) shared(intval) NUM_THREADS(nthreads)
        intvalp=0D0
!$OMP do schedule(DYNAMIC)
        do iz=2,nz-1
            ptz=zarr(iz)
            do iy=2,ny-1
                pty=yarr(iy)
                do ix=2,nx-1
                    ptx=xarr(ix)
                    if (gridbas(ix,iy,iz)==iatt) then
                        !Calculate switching function at current grid
                        rx=ptx-CPpos(1,iatt) !The relative distance between current point to corresponding attractor
                        ry=pty-CPpos(2,iatt)
                        rz=ptz-CPpos(3,iatt)
                        dist=dsqrt(rx*rx+ry*ry+rz*rz)
                        tmps=dist-trustrad(iatt)
                        if (tmps>1) then
                            switchwei=0
                        else if (tmps<-1) then
                            switchwei=1
                        else
                            do iter=1,nbeckeiter
                                tmps=1.5D0*(tmps)-0.5D0*(tmps)**3
                            end do
                            switchwei=0.5D0*(1-tmps)
                        end if
                        switchwei=1-switchwei
                        prodens=fdens(ptx,pty,ptz)
                        prodensgrad2=fgrad(ptx,pty,ptz,'t')**2
                        tmpval1=rhogrid(ix,iy,iz)*log(rhogrid(ix,iy,iz)/prodens)
                        tmpval2=rhograd2grid(ix,iy,iz)/rhogrid(ix,iy,iz)-prodensgrad2/prodens
                        intvalp(iatt,1)=intvalp(iatt,1)+tmpval1*switchwei
                        intvalp(iatt,2)=intvalp(iatt,2)+tmpval2*switchwei
                    end if
                end do
            end do
        end do
!$OMP end do
!$OMP CRITICAL
        intval=intval+intvalp*dvol
!$OMP end CRITICAL
!$OMP END PARALLEL
    end do
    deallocate(rhogrid,rhograd2grid)
    call dealloall
    write(*,"(' Reloading ',a)") trim(firstfilename)
    call readinfile(firstfilename,1) !Retrieve to first loaded file(whole molecule)
else if (ifuncint==-1) then !Deformation density
    allocate(prorhogrid(nx,ny,nz))
    write(*,*)
    write(*,*) "Calculating promolecular density at each grid"
    prorhogrid=0D0
    do iatm=1,ncenter_org !Cycle each atom
        write(*,"(' Processing atom',i6,a,'...')") iatm,a_org(iatm)%name
        call dealloall
        call readwfn(custommapname(iatm),1)
nthreads=getNThreads()
!$OMP PARALLEL DO SHARED(prorhogrid) PRIVATE(ix,iy,iz) schedule(dynamic) NUM_THREADS(nthreads)
        do iz=2,nz-1
            do iy=2,ny-1
                do ix=2,nx-1
                    prorhogrid(ix,iy,iz)=prorhogrid(ix,iy,iz)+fdens(xarr(ix),yarr(iy),zarr(iz))
                end do
            end do
        end do
!$OMP end PARALLEL DO
    end do
    do iatt=1,numrealatt !Cycle each attractors
        do iz=2,nz-1
            do iy=2,ny-1
                do ix=2,nx-1
                    if (gridbas(ix,iy,iz)==iatt) then
                        !Calculate switching function at current grid
                        rx=xarr(ix)-CPpos(1,iatt) !The relative distance between current point to corresponding attractor
                        ry=yarr(iy)-CPpos(2,iatt)
                        rz=zarr(iz)-CPpos(3,iatt)
                        dist=dsqrt(rx*rx+ry*ry+rz*rz)
                        tmps=dist-trustrad(iatt)
                        if (tmps>1) then
                            switchwei=0
                        else if (tmps<-1) then
                            switchwei=1
                        else
                            do iter=1,nbeckeiter
                                tmps=1.5D0*(tmps)-0.5D0*(tmps)**3
                            end do
                            switchwei=0.5D0*(1-tmps)
                        end if
                        switchwei=1-switchwei
                        defdens=cubmat(ix,iy,iz)-prorhogrid(ix,iy,iz)
                        intval(iatt,1)=intval(iatt,1)+defdens*switchwei*dvol
                    end if
                end do
            end do
        end do
    end do
    deallocate(prorhogrid)
    call dealloall
    call readinfile(firstfilename,1) !Retrieve to first loaded file(whole molecule) to calc real rho again
end if

!!----------- Output SUMMARY
if (itype==1.or.itype==2.or.itype==3) then !Integrate specific real space function(s)
    write(*,*)
    if (ifuncint==1) then
        write(*,*) "Total result:"
        write(*,*) "  #Basin        Integral(a.u.)      Vol(Bohr^3)    Vol(rho>0.001)"
        do iatt=1,numrealatt
            write(*,"(i8,f22.10,2f16.3)") iatt,intval(iatt,1),basinvol(iatt),basinvdwvol(iatt)
        end do
        write(*,"(' Sum of above integrals:',f20.8)") sum(intval(1:numrealatt,1))
        write(*,"(' Sum of basin volumes (rho>0.001):',f12.3,' Bohr^3')") sum(basinvdwvol(1:numrealatt))
        if (any(gridbas(2:nx-1,2:ny-1,2:nz-1)==0)) write(*,"(' Integral of unassigned grids:',f20.8)") intval(0,1)
        if (any(gridbas(2:nx-1,2:ny-1,2:nz-1)==-1)) write(*,"(' Integral of the grids travelled to box boundary:',f20.8)") intval(-1,1)
        write(*,*)
        rnormfac=sum(intval(1:numrealatt,1))/(nelec+nEDFelec) !The electrons represented by EDF must be taken into account!
        write(*,"(' Normalization factor of the integral of electron density is',f12.6)") rnormfac
        write(*,*) "The atomic charges after normalization and atomic volumes:"
        do iatm=0,ncenter
            do iatt=1,numrealatt
                if (att2atm(iatt)==iatm.and.iatm==0) then
                    write(*,"(i7,' (NNA)   Charge:',f12.6,'     Volume:',f10.3,' Bohr^3')") iatt,-intval(iatt,1)/rnormfac,basinvdwvol(iatt)
                else if (att2atm(iatt)==iatm.and.iatm/=0) then
                    if (nEDFelec==0) then !Normal case, all electron basis or using pseudopotential but not accompanied by EDF
                        write(*,"(i7,' (',a,')    Charge:',f12.6,'     Volume:',f10.3,' Bohr^3')") iatm,a(iatm)%name,a(iatm)%charge-intval(iatt,1)/rnormfac,basinvdwvol(iatt)
                    else !EDF is used, so using a(iatm)%index instead of a(iatm)%charge
                        write(*,"(i7,' (',a,')    Charge:',f12.6,'     Volume:',f10.3,' Bohr^3')") iatm,a(iatm)%name,a(iatm)%index-intval(iatt,1)/rnormfac,basinvdwvol(iatt)
                    end if
                end if
            end do
        end do
    else
        if (ispecial==0) then
            write(*,*) "Total result:"
            write(*,*) "    Atom       Basin       Integral(a.u.)   Vol(Bohr^3)   Vol(rho>0.001)"
            do iatm=0,ncenter
                do iatt=1,numrealatt
                    if (att2atm(iatt)==iatm.and.iatm==0) then
                        write(*,"('      NNA   ',i8,f20.8,2f14.3)") iatt,intval(iatt,1),basinvol(iatt),basinvdwvol(iatt)
                    else if (att2atm(iatt)==iatm.and.iatm/=0) then
                        write(*,"(i7,' (',a,')',i8,f20.8,2f14.3)") iatm,a(iatm)%name,iatt,intval(iatt,1),basinvol(iatt),basinvdwvol(iatt)
                    end if
                end do
            end do
            write(*,"(' Sum of above integrals:',f23.8)") sum(intval(1:numrealatt,1))
            write(*,"(' Sum of basin volumes (rho>0.001):',f12.3,' Bohr^3')") sum(basinvdwvol(1:numrealatt))
            if (any(gridbas(2:nx-1,2:ny-1,2:nz-1)==0)) write(*,"(' Integral of unassigned grids:',f20.8)") intval(0,1)
            if (any(gridbas(2:nx-1,2:ny-1,2:nz-1)==-1)) write(*,"(' Integral of the grids travelled to box boundary:',f20.8)") intval(-1,1)
        else if (ispecial==1) then
            write(*,*) "Total result:"
            write(*,*) "    Atom       Basin     Shannon       Fisher      Steric ene  Vol(rho>0.001)"
            do iatm=0,ncenter
                do iatt=1,numrealatt
                    if (att2atm(iatt)==iatm.and.iatm==0) then
                        write(*,"('      NNA   ',i8,3f14.7,f13.5)") iatt,intval(iatt,1),intval(iatt,1:3),basinvdwvol(iatt)
                    else if (att2atm(iatt)==iatm.and.iatm/=0) then
                        write(*,"(i7,' (',a,')',i8,3f14.7,f13.5)") iatm,a(iatm)%name,iatt,intval(iatt,1:3),basinvdwvol(iatt)
                    end if
                end do
            end do
            write(*,"(' Total Shannon entropy:   ',f23.8)") sum(intval(1:numrealatt,1))
            write(*,"(' Total Fisher information:',f23.8)") sum(intval(1:numrealatt,2))
            write(*,"(' Total steric energy:     ',f23.8)") sum(intval(1:numrealatt,3))
            write(*,"(' Sum of basin volumes (rho>0.001):',f12.3,' Bohr^3')") sum(basinvdwvol(1:numrealatt))
            rnormfac=sum(intval(1:numrealatt,4))/(nelec+nEDFelec) !The electrons represented by EDF must be taken into account!
            write(*,"(/,' Normalization factor of the integral of electron density is',f12.6)") rnormfac
            write(*,*) "The atomic charges after normalization and atomic volumes:"
            do iatm=0,ncenter
                do iatt=1,numrealatt
                    if (att2atm(iatt)==iatm.and.iatm==0) then
                        write(*,"(i7,' (NNA)   Charge:',f12.6,'     Volume:',f10.3,' Bohr^3')") iatt,-intval(iatt,1)/rnormfac,basinvdwvol(iatt)
                    else if (att2atm(iatt)==iatm.and.iatm/=0) then
                        if (nEDFelec==0) then !Normal case, all electron basis or using pseudopotential but not accompanied by EDF
                            write(*,"(i7,' (',a,')    Charge:',f12.6,'     Volume:',f10.3,' Bohr^3')") iatm,a(iatm)%name,a(iatm)%charge-intval(iatt,4)/rnormfac,basinvdwvol(iatt)
                        else !EDF is used, so using a(iatm)%index instead of a(iatm)%charge
                            write(*,"(i7,' (',a,')    Charge:',f12.6,'     Volume:',f10.3,' Bohr^3')") iatm,a(iatm)%name,a(iatm)%index-intval(iatt,4)/rnormfac,basinvdwvol(iatt)
                        end if
                    end if
                end do
            end do
        else if (ispecial==2) then
            write(*,*) "Total result:"
            write(*,*) "     Atom       Basin         Relat_Shannon           Relat_Fisher"
            do iatm=1,ncenter
                do iatt=1,numrealatt
                    if (att2atm(iatt)==iatm) write(*,"(i8,' (',a,')',i8,2f23.8)") iatm,a(iatm)%name,iatt,intval(iatt,1:2)
                end do
            end do
            write(*,"(' Sum of relat_Shannon:',f23.8)") sum(intval(1:numrealatt,1))
            write(*,"(' Sum of relat_Fisher: ',f23.8)") sum(intval(1:numrealatt,2))
        end if
    end if
    write(*,*)
else if (itype==10) then !Electric multipole moment
    ioutid=6
101    eleinttot=0D0
    dipelextot=0D0
    dipeleytot=0D0
    dipeleztot=0D0
    if (ioutid==6) write(*,*)
    write(ioutid,*) "Note: All units shown below are in a.u.!"
    write(ioutid,*)
    do iatm=0,ncenter
        do iatt=1,numrealatt
            if (att2atm(iatt)/=iatm) cycle
            if (iatm==0) then
                write(ioutid,"(' Result of NNA',i6)") iatt
            else
                write(ioutid,"(' Result of atom',i6,' (',a,')')") iatm,a(iatm)%name
            end if
            write(ioutid,"(' Basin electric monopole moment:',f12.6)") -eleint(iatt)
            write(ioutid,"(' Basin electric dipole moment:')") 
            write(ioutid,"('  X=',f12.6,'  Y=',f12.6,'  Z=',f12.6,'  Magnitude=',f12.6)") -xint(iatt),-yint(iatt),-zint(iatt),sqrt(xint(iatt)**2+yint(iatt)**2+zint(iatt)**2)
            eleinttot=eleinttot+eleint(iatt)
            dipelex=-eleint(iatt)*CPpos(1,iatt)+(-xint(iatt)) !Contribution to molecular total dipole moment
            dipeley=-eleint(iatt)*CPpos(2,iatt)+(-yint(iatt))
            dipelez=-eleint(iatt)*CPpos(3,iatt)+(-zint(iatt))
            dipelextot=dipelextot+dipelex
            dipeleytot=dipeleytot+dipeley
            dipeleztot=dipeleztot+dipelez
            write(ioutid,"(' Basin electron contribution to molecular dipole moment:')")
            write(ioutid,"('  X=',f12.6,'  Y=',f12.6,'  Z=',f12.6,'  Magnitude=',f12.6)") dipelex,dipeley,dipelez,sqrt(dipelex**2+dipeley**2+dipelez**2)
            write(ioutid,"(' Basin electric quadrupole moment (Cartesian form):')")
            rrint=xxint(iatt)+yyint(iatt)+zzint(iatt)
            QXX=-(3*xxint(iatt)-rrint)/2
            QYY=-(3*yyint(iatt)-rrint)/2
            QZZ=-(3*zzint(iatt)-rrint)/2
            write(ioutid,"(' QXX=',f12.6,'  QXY=',f12.6,'  QXZ=',f12.6)") QXX,-(3*xyint(iatt))/2,-(3*xzint(iatt))/2
            write(ioutid,"(' QYX=',f12.6,'  QYY=',f12.6,'  QYZ=',f12.6)") -(3*xyint(iatt))/2,QYY,-(3*yzint(iatt))/2
            write(ioutid,"(' QZX=',f12.6,'  QZY=',f12.6,'  QZZ=',f12.6)") -(3*xzint(iatt))/2,-(3*yzint(iatt))/2,QZZ
            write(ioutid,"( ' The magnitude of electric quadrupole moment (Cartesian form):',f12.6)") dsqrt(2D0/3D0*(QXX**2+QYY**2+QZZ**2))
            R20=-(3*zzint(iatt)-rrint)/2D0 !Notice that the negative sign, because electrons carry negative charge
            R2n1=-dsqrt(3D0)*yzint(iatt)
            R2p1=-dsqrt(3D0)*xzint(iatt)
            R2n2=-dsqrt(3D0)*xyint(iatt)
            R2p2=-dsqrt(3D0)/2D0*(xxint(iatt)-yyint(iatt))
            write(ioutid,"(' Electric quadrupole moments (Spherical harmonic form):')")
            write(ioutid,"(' Q_2,0 =',f11.6,'   Q_2,-1=',f11.6,'   Q_2,1=',f11.6)") R20,R2n1,R2p1
            write(ioutid,"(' Q_2,-2=',f11.6,'   Q_2,2 =',f11.6)") R2n2,R2p2
            write(ioutid,"( ' Magnitude: |Q_2|=',f12.6)") dsqrt(R20**2+R2n1**2+R2p1**2+R2n2**2+R2p2**2)
            write(ioutid,*)
        end do
    end do
    !Output overall electric properties
    dipnucx=sum(a(:)%x*a(:)%charge)
    dipnucy=sum(a(:)%y*a(:)%charge)
    dipnucz=sum(a(:)%z*a(:)%charge)
    write(ioutid,"( ' Molecular net charge:',f12.6)") sum(a%charge)-eleinttot
    write(ioutid,"( ' Nuclear contribution to molecular dipole moment:')") 
    write(ioutid,"(' X=',f12.6,' Y=',f12.6,' Z=',f12.6,' Norm=',f12.6)") dipnucx,dipnucy,dipnucz,sqrt(dipnucx**2+dipnucy**2+dipnucz**2)
    write(ioutid,"( ' Electron contribution to molecular dipole moment:')") 
    write(ioutid,"(' X=',f12.6,' Y=',f12.6,' Z=',f12.6,' Norm=',f12.6)") dipelextot,dipeleytot,dipeleztot,sqrt(dipelextot**2+dipeleytot**2+dipeleztot**2)
    dipmolx=dipnucx+dipelextot
    dipmoly=dipnucy+dipeleytot
    dipmolz=dipnucz+dipeleztot
    write(ioutid,"( ' Molecular dipole moment:')")
    write(ioutid,"(' X=',f12.6,' Y=',f12.6,' Z=',f12.6,' Norm=',f12.6)") dipmolx,dipmoly,dipmolz,sqrt(dipmolx**2+dipmoly**2+dipmolz**2)
    if (ioutid==6) write(*,*)
    if (ioutid==10) then
        close(10)
        write(*,*) "Done!"
        return
    end if
end if

CALL CPU_TIME(time_end)
call walltime(walltime2)
write(*,"(' Integrating basins took up CPU time',f12.2,'s, wall clock time',i10,'s')") time_end-time_begin,walltime2-walltime1

if (itype==10) then
    write(*,*)
    write(*,*) "If also output the result to multipol.txt in current folder? y/n"
    read(*,*) selectyn
    if (selectyn=='y'.or.selectyn=='Y') then
        ioutid=10
        open(10,file="multipol.txt",status="replace")
        goto 101
    end if
end if
end subroutine





!!------ Generate grid data of gradient of electron density, used to refine basin boundary
subroutine gengradmat
use defvar
use function
implicit real*8 (a-h,o-z)
real*8 wfnval(nmo),wfnderv(3,nmo),gradrho(3),EDFgrad(3),sumgrad2
if (allocated(cubmatvec)) then
    if (size(cubmatvec,2)==nx.and.size(cubmatvec,3)==ny.and.size(cubmatvec,4)==nz) return !Already generated
    deallocate(cubmatvec)
end if
allocate(cubmatvec(3,nx,ny,nz))
ifinish=0
nthreads=getNThreads()
!$OMP PARALLEL DO SHARED(cubmatvec,ifinish) PRIVATE(ix,iy,iz,tmpx,tmpy,tmpz,wfnval,wfnderv,gradrho,imo) schedule(dynamic) NUM_THREADS(nthreads)
do iz=1,nz
    tmpz=orgz+(iz-1)*dz
    do iy=1,ny
        tmpy=orgy+(iy-1)*dy
        do ix=1,nx
            tmpx=orgx+(ix-1)*dx
            call orbderv(2,1,nmo,tmpx,tmpy,tmpz,wfnval,wfnderv)
            gradrho=0D0
            do imo=1,nmo
                gradrho(:)=gradrho(:)+MOocc(imo)*wfnval(imo)*wfnderv(:,imo)
            end do
            cubmatvec(:,ix,iy,iz)=2*gradrho(:)
        end do
    end do
    ifinish=ifinish+1
    write(*,"(' Generating grid data of gradient, finished:',i5,'  /',i5)") ifinish,nz
end do
!$OMP END PARALLEL DO
end subroutine
