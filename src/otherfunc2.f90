!!-------- RDG analysis for MD, use promolecular approximation
subroutine RDG_MD
use defvar
use util
use function
implicit real*8 (a-h,o-z)
!The first index of avggrad and the first two indices of avghess correspond to components of gradient and Hessian, respectively
real*8,allocatable :: avgdens(:,:,:),avggrad(:,:,:,:),avghess(:,:,:,:,:)
real*8 denstmp,gradtmp(3),hesstmp(3,3),eigvecmat(3,3),eigval(3)
real*8,allocatable :: avgRDG(:,:,:),thermflu(:,:,:),avgsl2r(:,:,:) !Final result
real*8,allocatable :: scatterx(:),scattery(:)
integer walltime1,walltime2
real*8 time_end,time_endtmp

write(*,*) "Input the range of the frames to be analyzed, e.g. 150,400"
write(*,*) "Note: The index of frame starts from 1"
read(*,*) ifpsstart,ifpsend
call setgrid(0,igridsel)

nfps=ifpsend-ifpsstart+1
write(*,"(' Totally',i8,' frames, from frame',i8,' to',i8,' will be processed')") nfps,ifpsstart,ifpsend

!Generate averaged density, averaged gradient and averaged Hessian
allocate(avgdens(nx,ny,nz),avggrad(3,nx,ny,nz),avghess(3,3,nx,ny,nz))
avgdens=0D0
avggrad=0D0
avghess=0D0
call walltime(walltime1)
CALL CPU_TIME(time_begin)
open(10,file=filename,access="sequential",status="old")
write(*,*)
write(*,*) "Now calculate averaged density, density gradient and density Hessian"
do ifps=1,ifpsend
    call readxyz(filename,1,0)
    if (ifps<ifpsstart) cycle
    write(*,"(' Processing frame',i8,' ...')") ifps
    !fragatm must be defined for each frame, because each frame may have different content, and calchessmat_prodens will use it
    nfragatmnum=ncenter
    if (allocated(fragatm)) deallocate(fragatm)
    allocate(fragatm(nfragatmnum))
    do itmp=1,nfragatmnum
        fragatm(itmp)=itmp
    end do
nthreads=getNThreads()
nthreads=getNThreads()
!$OMP PARALLEL DO SHARED(avgdens,avggrad,avghess) PRIVATE(i,j,k,tmpx,tmpy,tmpz,denstmp,gradtmp,hesstmp) schedule(dynamic) NUM_THREADS(nthreads)
    do k=1,nz
        tmpz=orgz+(k-1)*dz
        do j=1,ny
            tmpy=orgy+(j-1)*dy
            do i=1,nx
                tmpx=orgx+(i-1)*dx
                call calchessmat_prodens(tmpx,tmpy,tmpz,denstmp,gradtmp,hesstmp)
                avgdens(i,j,k)=avgdens(i,j,k)+denstmp
                avggrad(:,i,j,k)=avggrad(:,i,j,k)+gradtmp
                avghess(:,:,i,j,k)=avghess(:,:,i,j,k)+hesstmp
            end do
        end do
    end do
!$OMP END PARALLEL DO
end do
close(10)
avgdens=avgdens/nfps
avggrad=avggrad/nfps
avghess=avghess/nfps

!Use averaged density, averaged gradient and averaged Hessian to compute averaged RDG and averaged Sign(lambda2)*rho
write(*,*)
write(*,*) "Calculating averaged RDG and averaged Sign(lambda2)*rho..."
allocate(avgRDG(nx,ny,nz),avgsl2r(nx,ny,nz))
do k=1,nz
    tmpz=orgz+(k-1)*dz
    do j=1,ny
        tmpy=orgy+(j-1)*dy
        do i=1,nx
            tmpx=orgx+(i-1)*dx

            avggradnormtmp=dsqrt(sum(avggrad(:,i,j,k)**2))
            if (RDGprodens_maxrho/=0D0.and.avgdens(i,j,k)>=RDGprodens_maxrho) then
                avgRDG(i,j,k)=100D0
            else if (avggradnormtmp==0D0.or.avgdens(i,j,k)==0D0) then
                avgRDG(i,j,k)=999D0
            else
                avgRDG(i,j,k)=0.161620459673995D0*avggradnormtmp/avgdens(i,j,k)**(4D0/3D0) !0.161620459673995D0=1/(2*(3*pi**2)**(1/3))
            end if
            
            call diagmat(avghess(:,:,i,j,k),eigvecmat,eigval,100,1D-6)
            call sort(eigval)
            if (eigval(2)/=0D0) then
                avgsl2r(i,j,k)=avgdens(i,j,k)*eigval(2)/abs(eigval(2)) !At nuclei of single atom system, Hessian returned may be zero matrix
            else
                avgsl2r(i,j,k)=-avgdens(i,j,k) !Around nuclei, eigval(2)/abs(eigval(2)) always be negative
            end if
            
        end do
    end do
end do
CALL CPU_TIME(time_end)
call walltime(walltime2)
write(*,"(' Calculation totally took up CPU time',f12.2,'s, wall clock time',i10,'s')") time_end-time_begin,walltime2-walltime1

deallocate(avghess) !avghess will not be used further, so release its memory
allocate(scatterx(nx*ny*nz),scattery(nx*ny*nz))
ii=1
do k=1,nz
    do j=1,ny
        do i=1,nx
            scatterx(ii)=avgsl2r(i,j,k)
            scattery(ii)=avgRDG(i,j,k)
            ii=ii+1
        end do
    end do
end do
!Default axis range of scatter plot
xmin=-RDGprodens_maxrho
xmax=RDGprodens_maxrho
if (RDGprodens_maxrho==0.0D0) xmin=-2.0D0
if (RDGprodens_maxrho==0.0D0) xmax=2.0D0
ymin=0.0D0
ymax=1.0D0

write(*,*)
do while (.true.)
    write(*,*) "0 Return"
    write(*,*) "1 Draw scatter graph between averaged RDG and Sign(lambda2)*rho"
    write(*,*) "2 Save the scatter graph to file"
    write(*,*) "3 Change range of X-axis of the scatter graph"
    write(*,*) "4 Change range of Y-axis of the scatter graph"
    write(*,*) "5 Export scatter points to output.txt"
    write(*,*) "6 Export averaged RDG and Sign(lambda2)*rho to cube files in current folder"
    write(*,*) "7 Compute and export thermal fluctuation index to cube file in current folder"    
    write(*,*) "8 Show isosurface of averaged RDG"
    read(*,*) isel
    
    if (isel==0) then
        exit
    else if (isel==1) then
        write(*,*) "Drawing graph, please wait..."
    else if (isel==2) then
        isavepic=1
        isavepic=0
        write(*,"(a,a,a)") " Graph have been saved to ",trim(graphformat)," file with ""DISLIN"" prefix in current directory"
    else if (isel==3) then
        write(*,*) "Input lower limit and upper limit of X axis e.g. 0,1.5"
        read(*,*) xmin,xmax
    else if (isel==4) then
        write(*,*) "Input lower limit and upper limit of Y axis e.g. 0,1.5"
        read(*,*) ymin,ymax
    else if (isel==5) then
        open(10,file="output.txt",status="replace")
        write(*,*) "Outputting output.txt in current folder..."
        write(10,"(3f11.6,2E16.8)") ((( (orgx+(i-1)*dx),(orgy+(j-1)*dy),(orgz+(k-1)*dz),avgRDG(i,j,k),avgsl2r(i,j,k),i=1,nx),j=1,ny),k=1,nz)
        close(10)
        write(*,"(a)") " Finished, column 1/2/3/4/5 = X/Y/Z/averaged RDG/averaged Sign(lambda2)*rho, unit is Bohr"
    else if (isel==6) then
        write(*,*) "Outputting averaged reduced density gradient to avgRDG.cub in current folder"
        open(10,file="avgRDG.cub",status="replace")
        call outcube(avgRDG,nx,ny,nz,orgx,orgy,orgz,dx,dy,dz,10)
        close(10)
        write(*,*) "Done!"
        write(*,*)
        write(*,*) "Outputting averaged Sign(lambda2)*rho to avgsl2r.cub in current folder"
        open(10,file="avgsl2r.cub",status="replace")
        call outcube(avgsl2r,nx,ny,nz,orgx,orgy,orgz,dx,dy,dz,10)
        close(10)
        write(*,*) "Done!"
    else if (isel==7) then
        call walltime(walltime1)
        CALL CPU_TIME(time_begin)
        write(*,*) "Calculating thermal fluctuation index..."
        if (allocated(thermflu)) deallocate(thermflu)
        allocate(thermflu(nx,ny,nz))
        thermflu=0D0
        open(10,file=filename,access="sequential",status="old")
        !Calculate fluctuation square term of density and temporarily store it to thermflu array
        do ifps=1,ifpsend
            call readxyz(filename,1,0)
            if (ifps<ifpsstart) cycle
            write(*,"(' Processing frame',i8,' ...')") ifps
nthreads=getNThreads()
nthreads=getNThreads()
!$OMP PARALLEL DO SHARED(thermflu) PRIVATE(i,j,k,tmpx,tmpy,tmpz,denstmp) schedule(dynamic) NUM_THREADS(nthreads)
            do k=1,nz
                tmpz=orgz+(k-1)*dz
                do j=1,ny
                    tmpy=orgy+(j-1)*dy
                    do i=1,nx
                        tmpx=orgx+(i-1)*dx
                        call calchessmat_prodens(tmpx,tmpy,tmpz,denstmp)
                        thermflu(i,j,k)=thermflu(i,j,k)+(denstmp-avgdens(i,j,k))**2
                    end do
                end do
            end do
!$OMP END PARALLEL DO
        end do
        close(10)
        thermflu=dsqrt(thermflu/nfps)/avgdens
        CALL CPU_TIME(time_end)
        call walltime(walltime2)
        write(*,"(' Totally took up CPU time',f12.2,'s, wall clock time',i10,'s',/)") time_end-time_begin,walltime2-walltime1
        write(*,*) "Outputting thermal fluctuation index to thermflu.cub in current folder"
        open(10,file="thermflu.cub",status="replace")
        call outcube(thermflu,nx,ny,nz,orgx,orgy,orgz,dx,dy,dz,10)
        close(10)
        write(*,*) "Done!"
    else if (isel==8) then
        write(*,*) "Please wait..."
         sur_value=0.3D0
         if (allocated(cubmat)) deallocate(cubmat)
         allocate(cubmat(nx,ny,nz))
         cubmat=avgRDG
        deallocate(cubmat)
    end if
    write(*,*)
end do

end subroutine


!!----- Calculate atomic and bond dipole moments in Hilbert space
!For derivation, see Ideas of Quantum Chemistry, p634
subroutine atmbonddip
use defvar
use util
implicit real*8 (a-h,o-z)
real*8 xdipmat(nbasis,nbasis),ydipmat(nbasis,nbasis),zdipmat(nbasis,nbasis),Ptottmp(nbasis,nbasis)
character c80tmp*80
!!!!Beware that the dipole moment integral has taken the negative sign of electron charge into account!
if (igenDbas==0) then !Haven't calculated dipole moment integral matrix, so reload the input file and calculate it
    Ptottmp=Ptot !Backup Ptot, which may have already been modified by users (via modifying occ), otherwise it will be flushed when loading
    igenDbas=1
    write(*,*) "Reloading input file and meantime generating dipole moment integral matrix..."
    call dealloall
    call readinfile(firstfilename,1)
    Ptot=Ptottmp
end if
xdipmat=Dbas(1,:,:)
ydipmat=Dbas(2,:,:)
zdipmat=Dbas(3,:,:)

! call showmatgau(xdipmat,"dip x",0,"f14.8",6)
! call showmatgau(ydipmat,"dip y",0,"f14.8",6)
! call showmatgau(zdipmat,"dip z",0,"f14.8",6)

!Calculate total dipole moment
xnucdip=sum(a(:)%charge*a(:)%x)
ynucdip=sum(a(:)%charge*a(:)%y)
znucdip=sum(a(:)%charge*a(:)%z)
write(*,"(' Molecular nuclear dipole moment (a.u.):')")
write(*,"('  X=',f12.6,'  Y=',f12.6,'  Z=',f12.6,'  Norm=',f12.6)") xnucdip,ynucdip,znucdip,dsqrt(xnucdip**2+ynucdip**2+znucdip**2)
xeledip=sum(Ptot*xdipmat)
yeledip=sum(Ptot*ydipmat)
zeledip=sum(Ptot*zdipmat)
write(*,"(' Molecular electron dipole moment (a.u.):')")
write(*,"('  X=',f12.6,'  Y=',f12.6,'  Z=',f12.6,'  Norm=',f12.6)") xeledip,yeledip,zeledip,dsqrt(xeledip**2+yeledip**2+zeledip**2)
xmoldip=xnucdip+xeledip
ymoldip=ynucdip+yeledip
zmoldip=znucdip+zeledip
write(*,"(' Molecular dipole moment (a.u.):')")
write(*,"('  X=',f12.6,'  Y=',f12.6,'  Z=',f12.6,'  Norm=',f12.6)") xmoldip,ymoldip,zmoldip,dsqrt(xmoldip**2+ymoldip**2+zmoldip**2)

do while(.true.)
    write(*,*)
    write(*,*) "         ----- Atomic and bond dipole moments in Hilbert space -----"
    write(*,*) "0 Return"
    write(*,*) "1 Output atomic dipole moment of specific atom"
    write(*,*) "2 Output bond dipole moment of specific atom pair"
    write(*,*) "3 Output atomic overall dipole moment of specific atom (Mulliken partition)"
    write(*,*) "10 Export entire dipole moment matrix"
    read(*,*) isel
    if (isel==0) then
        exit
    else if (isel==1) then
        do while(.true.)
            write(*,*) "Input the atom index, e.g. 5"
            write(*,*) "Note: Input 0 can return, input -1 can output result for all atoms"
            read(*,*) isel2
            if (isel2==0) then
                exit
            else if (isel2==-1) then
                iatmstart=1
                iatmend=ncenter
            else
                if (isel2>ncenter) then
                    write(*,*) "Atom index exceeded valid range! Input again"
                    cycle
                end if
                iatmstart=isel2
                iatmend=isel2
            end if
            do iatm=iatmstart,iatmend
                write(*,"(' Result of atom',i8,' (',a2,')')") iatm,a(iatm)%name
                istart=basstart(iatm)
                iend=basend(iatm)
                atmelepop=sum(Ptot(istart:iend,istart:iend)*Sbas(istart:iend,istart:iend))
                xatmelediptot=sum(Ptot(istart:iend,istart:iend)*xdipmat(istart:iend,istart:iend))
                yatmelediptot=sum(Ptot(istart:iend,istart:iend)*ydipmat(istart:iend,istart:iend))
                zatmelediptot=sum(Ptot(istart:iend,istart:iend)*zdipmat(istart:iend,istart:iend))
                xatmeledip=xatmelediptot+atmelepop*a(iatm)%x
                yatmeledip=yatmelediptot+atmelepop*a(iatm)%y
                zatmeledip=zatmelediptot+atmelepop*a(iatm)%z
                xatmnucdip=a(iatm)%charge*a(iatm)%x
                yatmnucdip=a(iatm)%charge*a(iatm)%y
                zatmnucdip=a(iatm)%charge*a(iatm)%z
                xatmdiptot=xatmnucdip+xatmelediptot
                yatmdiptot=yatmnucdip+yatmelediptot
                zatmdiptot=zatmnucdip+zatmelediptot
                write(*,"(' Atomic local population number:',f12.6)") atmelepop
                write(*,"(' Atomic dipole moment (a.u.):')")
                write(*,"('  X=',f12.6,'  Y=',f12.6,'  Z=',f12.6,'  Norm=',f12.6)") xatmeledip,yatmeledip,zatmeledip,dsqrt(xatmeledip**2+yatmeledip**2+zatmeledip**2)
                write(*,"(' Contribution to system dipole moment due to nuclear charge (a.u.):')")
                write(*,"('  X=',f12.6,'  Y=',f12.6,'  Z=',f12.6,'  Norm=',f12.6)") xatmnucdip,yatmnucdip,zatmnucdip,dsqrt(xatmnucdip**2+yatmnucdip**2+zatmnucdip**2)
                write(*,"(' Contribution to system dipole moment due to electron (a.u.):')")
                write(*,"('  X=',f12.6,'  Y=',f12.6,'  Z=',f12.6,'  Norm=',f12.6)") xatmelediptot,yatmelediptot,zatmelediptot,dsqrt(xatmelediptot**2+yatmelediptot**2+zatmelediptot**2)
                write(*,"(' Contribution to system dipole moment (a.u.):')")
                write(*,"('  X=',f12.6,'  Y=',f12.6,'  Z=',f12.6,'  Norm=',f12.6)") xatmdiptot,yatmdiptot,zatmdiptot,dsqrt(xatmdiptot**2+yatmdiptot**2+zatmdiptot**2)
                write(*,*)
            end do
        end do
    else if (isel==2) then
        do while(.true.)
            write(*,*) "Input the index of two atoms, e.g. 5,8"
            write(*,*) "Note: Input q can return. Input b can output result for all bonds"
            read(*,"(a)") c80tmp
            if (index(c80tmp,'q')/=0) then
                exit
            else if (index(c80tmp,'b')/=0) then
                write(*,*) "Notice that the bonds are determined according to distance criterion"
                write(*,*)
                bondcritval=1.15D0
            else
                read(c80tmp,*) iatomsel1,iatomsel2
                if (iatomsel1>ncenter.or.iatomsel2>ncenter) then
                    write(*,*) "Atom index exceeded valid range! Input again"
                    cycle
                end if
                if (iatomsel1>iatomsel2) then
                    itmp=iatomsel2
                    iatomsel2=iatomsel1
                    iatomsel1=itmp
                end if
            end if
            do iatm=1,ncenter
                do jatm=iatm+1,ncenter
                    bonddist=dsqrt((a(iatm)%x-a(jatm)%x)**2+(a(iatm)%y-a(jatm)%y)**2+(a(iatm)%z-a(jatm)%z)**2)
                    if (index(c80tmp,'b')/=0) then
                        if (bonddist>( covr(a(iatm)%index)+covr(a(jatm)%index) )*bondcritval) cycle
                    else
                        if (iatm/=iatomsel1.or.jatm/=iatomsel2) cycle
                    end if
                    xcen=(a(iatm)%x+a(jatm)%x)/2D0
                    ycen=(a(iatm)%y+a(jatm)%y)/2D0
                    zcen=(a(iatm)%z+a(jatm)%z)/2D0
                    write(*,"(' Result between atom',i7,' (',a2,')  and atom',i7,' (',a2,'), distance:',f10.5,' Ang')") iatm,a(iatm)%name,jatm,a(jatm)%name,bonddist*b2a
                    istart=basstart(iatm)
                    iend=basend(iatm)
                    jstart=basstart(jatm)
                    jend=basend(jatm)
                    bondpop=2*sum(Ptot(istart:iend,jstart:jend)*Sbas(istart:iend,jstart:jend)) !The matrix is symmetical, so multiplied by 2
                    xbonddiptot=2*sum(Ptot(istart:iend,jstart:jend)*xdipmat(istart:iend,jstart:jend))
                    ybonddiptot=2*sum(Ptot(istart:iend,jstart:jend)*ydipmat(istart:iend,jstart:jend))
                    zbonddiptot=2*sum(Ptot(istart:iend,jstart:jend)*zdipmat(istart:iend,jstart:jend))
                    xbonddip=xbonddiptot+bondpop*xcen
                    ybonddip=ybonddiptot+bondpop*ycen
                    zbonddip=zbonddiptot+bondpop*zcen
                    write(*,"(' Bond population number (Overlap population):',f12.6)") bondpop
                    write(*,"(' Bond dipole moment (a.u.):')")
                    write(*,"('  X=',f12.6,'  Y=',f12.6,'  Z=',f12.6,'  Norm=',f12.6)") xbonddip,ybonddip,zbonddip,dsqrt(xbonddip**2+ybonddip**2+zbonddip**2)
                    write(*,"(' Contribution to system dipole moment (a.u.):')")
                    write(*,"('  X=',f12.6,'  Y=',f12.6,'  Z=',f12.6,'  Norm=',f12.6)") xbonddiptot,ybonddiptot,zbonddiptot,dsqrt(xbonddiptot**2+ybonddiptot**2+zbonddiptot**2)
                    write(*,*)
                end do
            end do
        end do
    else if (isel==3) then
        do while(.true.)
            write(*,*) "Input the atom index, e.g. 5"
            write(*,*) "Note: Input 0 can return, input -1 can output result for all atoms"
            read(*,*) isel2
            if (isel2==0) then
                exit
            else if (isel2==-1) then
                iatmstart=1
                iatmend=ncenter
            else
                if (isel2>ncenter) then
                    write(*,*) "Atom index exceeded valid range! Input again"
                    cycle
                end if
                iatmstart=isel2
                iatmend=isel2
            end if
            do iatm=iatmstart,iatmend
                write(*,"(' Result of atom',i8,' (',a2,')')") iatm,a(iatm)%name
                atmelepop=0D0
                xatmelediptot=0D0
                yatmelediptot=0D0
                zatmelediptot=0D0
                istart=basstart(iatm)
                iend=basend(iatm)
                do jatm=1,ncenter
                    jstart=basstart(jatm)
                    jend=basend(jatm)
                    atmelepop=atmelepop+sum(Ptot(istart:iend,jstart:jend)*Sbas(istart:iend,jstart:jend))
                    xatmelediptot=xatmelediptot+sum(Ptot(istart:iend,jstart:jend)*xdipmat(istart:iend,jstart:jend))
                    yatmelediptot=yatmelediptot+sum(Ptot(istart:iend,jstart:jend)*ydipmat(istart:iend,jstart:jend))
                    zatmelediptot=zatmelediptot+sum(Ptot(istart:iend,jstart:jend)*zdipmat(istart:iend,jstart:jend))
                end do
                xatmeledip=xatmelediptot+atmelepop*a(iatm)%x
                yatmeledip=yatmelediptot+atmelepop*a(iatm)%y
                zatmeledip=zatmelediptot+atmelepop*a(iatm)%z
                xatmnucdip=a(iatm)%charge*a(iatm)%x
                yatmnucdip=a(iatm)%charge*a(iatm)%y
                zatmnucdip=a(iatm)%charge*a(iatm)%z
                xatmdiptot=xatmnucdip+xatmelediptot
                yatmdiptot=yatmnucdip+yatmelediptot
                zatmdiptot=zatmnucdip+zatmelediptot
                write(*,"(' Atomic Mulliken population number:',f12.6)") atmelepop
                write(*,"(' Atomic overall dipole moment (a.u.):')")
                write(*,"('  X=',f12.6,'  Y=',f12.6,'  Z=',f12.6,'  Norm=',f12.6)") xatmeledip,yatmeledip,zatmeledip,dsqrt(xatmeledip**2+yatmeledip**2+zatmeledip**2)
                write(*,"(' Contribution to system dipole moment due to nuclear charge (a.u.):')")
                write(*,"('  X=',f12.6,'  Y=',f12.6,'  Z=',f12.6,'  Norm=',f12.6)") xatmnucdip,yatmnucdip,zatmnucdip,dsqrt(xatmnucdip**2+yatmnucdip**2+zatmnucdip**2)
                write(*,"(' Contribution to system dipole moment due to electron (a.u.):')")
                write(*,"('  X=',f12.6,'  Y=',f12.6,'  Z=',f12.6,'  Norm=',f12.6)") xatmelediptot,yatmelediptot,zatmelediptot,dsqrt(xatmelediptot**2+yatmelediptot**2+zatmelediptot**2)
                write(*,"(' Contribution to system dipole moment (a.u.):')")
                write(*,"('  X=',f12.6,'  Y=',f12.6,'  Z=',f12.6,'  Norm=',f12.6)") xatmdiptot,yatmdiptot,zatmdiptot,dsqrt(xatmdiptot**2+yatmdiptot**2+zatmdiptot**2)
                write(*,*)
            end do
        end do
    else if (isel==10) then
        open(10,file="dipmatx.txt",status="replace")
        call showmatgau(Ptot*xdipmat,"",1,"f14.8",10)
        close(10)
        open(10,file="dipmaty.txt",status="replace")
        call showmatgau(Ptot*ydipmat,"",1,"f14.8",10)
        close(10)
        open(10,file="dipmatz.txt",status="replace")
        call showmatgau(Ptot*zdipmat,"",1,"f14.8",10)
        close(10)
        write(*,"(a)") " X, Y and Z components of electron dipole moment matrix have been outputted to dipmatx, dipmaty and dipmatz.txt in current folder, respectively"
    end if
end do
end subroutine


!!------------ Generate cube file for multiple orbitals
subroutine genmultiorbcube
use defvar
use util
implicit real*8 (a-h,o-z)
integer orbsellist(nmo)
integer tmparr(nmo+1)
character c1000tmp*1000,cubname*20
real*8,allocatable :: orbcubmat(:,:,:,:)
write(*,"(a)") " Input orbital index. e.g. 1,3-6,8,10-11 denotes 1,3,4,5,6,8,10,11"
write(*,*) "Input q can return"
read(*,"(a)") c1000tmp
if (index(c1000tmp,'q')/=0) return
call str2arr(c1000tmp,norbsel,orbsellist)
if ( any(orbsellist(1:norbsel)<1) .or. any(orbsellist(1:norbsel)>nmo) ) then
    write(*,*) "Error: The orbitals you selected exceeded valid range!"
    return
end if
call setgrid(0,itmp)
if (allocated(cubmat)) deallocate(cubmat)
allocate(cubmat(nx,ny,nz))
write(*,*)
write(*,*) "1 Output the grid data of these orbitals as separate cube files"
write(*,*) "2 Output the grid data of these orbitals as a single cube file"
read(*,*) ioutmode

if (ioutmode==1) then
    do iorbidx=1,norbsel
        iorb=orbsellist(iorbidx)
        write(cubname,"('orb',i6.6,'.cub')") iorb
        write(*,"(' Calculating and exporting orbital',i6)") iorb
        call savecubmat(4,1,iorb)
        open(10,file=cubname,status="replace")
        call outcube(cubmat,nx,ny,nz,orgx,orgy,orgz,dx,dy,dz,10)
        close(10)
        write(*,"(' Orbital',i7,' has been exported to ',a,' in current folder',/)") iorb,trim(cubname)
    end do
    
else if (ioutmode==2) then
    allocate(orbcubmat(nx,ny,nz,norbsel))
    do iorbidx=1,norbsel
        iorb=orbsellist(iorbidx)
        write(*,"(a,i6,a)") " Calculating grid data for orbital",iorb,"..."
        call savecubmat(4,1,iorb)
        orbcubmat(:,:,:,iorbidx)=cubmat(:,:,:)
    end do
    where (abs(orbcubmat)<=1D-99) orbcubmat=0D0 !Diminish too small value, otherwise the symbol "E" cannot be shown by 1PE13.5 format e.g. 9.39376-116, 
    write(*,*)
    write(*,*) "Exporting cube file, please wait..."
    open(10,file="orbital.cub",status="replace")
    write(10,"(' Generated by Multiwfn')")
    write(10,"(' Totally ',i12,' grid points')") nx*ny*nz
    write(10,"(i5,3f12.6)") -ncenter,orgx,orgy,orgz
    write(10,"(i5,3f12.6)") nx,dx,0.0,0.0
    write(10,"(i5,3f12.6)") ny,0.0,dy,0.0
    write(10,"(i5,3f12.6)") nz,0.0,0.0,dz
    do i=1,ncenter
        write(10,"(i5,4f12.6)") a(i)%index,a(i)%charge,a(i)%x,a(i)%y,a(i)%z
    end do
    tmparr(1)=norbsel
    tmparr(2:norbsel+1)=orbsellist(1:norbsel)
    write(10,"(10i5)") tmparr(1:norbsel+1)
    do ix=1,nx
        do iy=1,ny
            write(10,"(6(1PE13.5))",advance="no") ((orbcubmat(ix,iy,iz,iorbidx),iorbidx=1,norbsel),iz=1,nz)
            write(10,*)
        end do
    end do
    close(10)
    write(*,*) "The grid data of the orbitals have been stored to orbital.cub in current folder"
    deallocate(orbcubmat)
end if

deallocate(cubmat)
end subroutine




!!---------- Generate grid data of iso-chemical shielding surfaces (ICSS) or related quantities
subroutine ICSS
use defvar
use util
implicit real*8 (a-h,o-z)
character c200tmp*200,gauinpfile*200,gauoutfile*200,selectyn,suffix*4
character,allocatable :: gauinpcontent(:)*79
!Set grid for calculating NICS
call setgrid(0,itmp)
numbqper=NICSnptlim-ncenter
write(*,"(' The number of Bq per batch:',i10)") numbqper
npttot=nx*ny*nz
nfile=ceiling(dfloat(npttot)/numbqper)
!Generate Gaussian input file
write(*,*)
write(*,*) "If skip generating Gaussian input file of NMR task? (y/n)"
read(*,*) selectyn
if (selectyn=='n'.or.selectyn=='N') then
    write(*,*) "Input the path of template Gaussian input file, e.g. c:\ltwd.gjf"
    do while(.true.)
        read(*,"(a)") gauinpfile
        inquire(file=gauinpfile,exist=alive)
        if (alive) exit
        write(*,*) "Cannot find corresponding files, input again"
    end do
    open(10,file=gauinpfile,status="old")
    numgauline=totlinenum(10,2)
    allocate(gauinpcontent(numgauline))
    numblank=0
    iendcoord=numgauline !Which line is the last line recording coordinates
    do i=1,numgauline
        read(10,"(a)") gauinpcontent(i)
        if (gauinpcontent(i)==" ") then
            numblank=numblank+1
            if (numblank==3) iendcoord=i-1
        end if
        if (index(gauinpcontent(i),'#')/=0) then
            gauinpcontent(i)=trim(gauinpcontent(i))//" NMR"
            if (index(gauinpcontent(i),'conn')==0) gauinpcontent(i)=trim(gauinpcontent(i))//" geom=connectivity"
        end if
!         write(*,"(a)") gauinpcontent(i)
    end do
    close(10)

    gauinpfile="NICS"
    do ifile=1,nfile
        write(c200tmp,"(a,i4.4,a)") trim(gauinpfile),ifile,".gjf"
        open(10,file=c200tmp,status="replace")
        write(*,"(a,a,a)") " Outputting ",trim(c200tmp)," to current folder..."
        !Write head part
        do i=1,iendcoord
            if (ifile>1.and.index(gauinpcontent(i),'#')/=0) then
                write(10,"(a)") trim(gauinpcontent(i))//" guess=read"
            else
                write(10,"(a)") gauinpcontent(i)
            end if
        end do
        !Write Bq information
        itmp=0
        do i=1,nx
            do j=1,ny
                do k=1,nz
                    itmp=itmp+1
                    rnowx=orgx+(i-1)*dx
                    rnowy=orgy+(j-1)*dy
                    rnowz=orgz+(k-1)*dz
                    if (itmp<=(ifile-1)*numbqper) cycle
                    if (itmp>ifile*numbqper) exit
                    write(10,"('Bq ',3f12.6)") rnowx*b2a,rnowy*b2a,rnowz*b2a
                end do
            end do
        end do
        write(10,*)
        !Write connectivity explicitly
        if (ifile/=nfile) then
            do i=1,NICSnptlim
                write(10,"(i7)") i
            end do
        else if (ifile==nfile) then
            do i=1,npttot-(ifile-1)*numbqper+ncenter
                write(10,"(i7)") i
            end do
        end if
        !Write remaining part
        do i=iendcoord+1,numgauline
            write(10,"(a)") gauinpcontent(i)
        end do
        close(10)
    end do
end if
write(*,*) "Now please run these input files by Gaussian"
write(*,*)
!Load NICS from Gaussian output file
if (allocated(cubmat)) deallocate(cubmat)
allocate(cubmat(nx,ny,nz))
write(*,*) "Input the path of Gaussian output file of NMR task"
write(*,"(a)") " Assume that you input ""c:\ltwd\NICS"", then c:\ltwd\NICS0001.out, c:\ltwd\NICS0002.out, c:\ltwd\NICS0003.out... will be loaded"
suffix=".out"
do while(.true.)
    read(*,"(a)") gauoutfile
    inquire(file=trim(gauoutfile)//"0001"//suffix,exist=alive)
    if (alive) exit
    if (.not.alive) then
        suffix=".log"
        inquire(file=trim(gauoutfile)//"0001"//suffix,exist=alive)
    end if
    if (alive) exit
    write(*,*) "Cannot find corresponding files, input again"
end do
100 write(*,*) "Load which term?"
write(*,*) "1: Isotropic  2: Anisotropy  3: XX component  4: YY component  5: ZZ component"
read(*,*) iload
do ifile=1,nfile
    write(c200tmp,"(a,i4.4,a)") trim(gauoutfile),ifile,suffix
    open(10,file=c200tmp,status="old")
    write(*,"(' Loading ',a,'...')") trim(c200tmp)
    call loclabel(10,"GIAO Magnetic shielding tensor")
    read(10,*)
    !Detect format. The NMR output format changes since G09 D.01 to leave more space for atomic index
    read(10,"(a80)") c200tmp
    backspace(10)
    iformat=1
    if (c200tmp(25:25)=='=') iformat=2 !Since G09 D.01
    
    do i=1,ncenter !Skip atom's result
        read(10,*)
        read(10,*)
        read(10,*)
        read(10,*)
        read(10,*)
    end do
    itmp=0
    do i=1,nx
        do j=1,ny
            do k=1,nz
                itmp=itmp+1
                if (itmp<=(ifile-1)*numbqper) cycle
                if (itmp>ifile*numbqper) exit
                if (iload==1.or.iload==2) then
!                     read(10,"(a80)") c200tmp
!                     write(*,"(a80)") c200tmp
!                     backspace(10)
                    if (iformat==1) then
                        if (iload==1) read(10,"(22x,f10.4)") cubmat(i,j,k)
                        if (iload==2) read(10,"(48x,f10.4)") cubmat(i,j,k)
                    else if (iformat==2) then
                        if (iload==1) read(10,"(26x,f10.4)") cubmat(i,j,k)
                        if (iload==2) read(10,"(52x,f10.4)") cubmat(i,j,k)
                    end if
                    read(10,*)
                    read(10,*)
                    read(10,*)
                    read(10,*)
                else
                    read(10,*)
                    if (iload==3) then
                        read(10,"(8x,f10.4)") cubmat(i,j,k)
                        read(10,*)
                        read(10,*)
                    else if (iload==4) then
                        read(10,*)
                        read(10,"(24x,f10.4)") cubmat(i,j,k)
                        read(10,*)
                    else if (iload==5) then
                        read(10,*)
                        read(10,*)
                        read(10,"(42x,f10.4)") cubmat(i,j,k)
                    end if
                    read(10,*)
                end if
            end do
        end do
    end do
    close(10)
end do
write(*,*) "Loading finished!"
do while(.true.)
    write(*,*)
    write(*,*) "-1 Load another property"
    write(*,*) "0 Return to main menu"
    if (iload==1) write(*,*) "1 Visualize iso-chemical shielding surface"
    if (iload==2) write(*,*) "1 Visualize aniso-chemical shielding surface"
    if (iload==3) write(*,*) "1 Visualize XX-chemical shielding surface"
    if (iload==4) write(*,*) "1 Visualize YY-chemical shielding surface"
    if (iload==5) write(*,*) "1 Visualize ZZ-chemical shielding surface"
    if (iload==1) write(*,*) "2 Export the grid data to ICSS.cub current folder"
    if (iload==2) write(*,*) "2 Export the grid data to AICSS.cub current folder"
    if (iload==3) write(*,*) "2 Export the grid data to ICSSXX.cub current folder"
    if (iload==4) write(*,*) "2 Export the grid data to ICSSYY.cub current folder"
    if (iload==5) write(*,*) "2 Export the grid data to ICSSZZ.cub current folder"
    read(*,*) isel
    if (isel==-1) then
        goto 100
    else if (isel==0) then
        exit
    else if (isel==1) then
    else if (isel==2) then
        if (iload==1) open(10,file="ICSS.cub",status="replace")
        if (iload==2) open(10,file="AICSS.cub",status="replace")
        if (iload==3) open(10,file="ICSSXX.cub",status="replace")
        if (iload==4) open(10,file="ICSSYY.cub",status="replace")
        if (iload==5) open(10,file="ICSSZZ.cub",status="replace")
        call outcube(cubmat,nx,ny,nz,orgx,orgy,orgz,dx,dy,dz,10)
        close(10)
        write(*,"(a)") " The cube file has been exported to current folder"
    end if
end do
end subroutine



!!----- Plot radial distribution function for a real space function
subroutine plotraddis
use defvar
use function
implicit real*8 (a-h,o-z)
real*8,allocatable :: potx(:),poty(:),potz(:),potw(:),radval(:),radpos(:),intradval(:),sphavgval(:)
ifunc=1
cenx=0D0
ceny=0D0
cenz=0D0
rlow=0D0
rhigh=5D0/b2a !5 Angstrom
nsphpt=2030
nradpt=500
do while(.true.)
    write(*,*)
    write(*,*) "  ====== Plot radial distribution function for a real space function ======"
    write(*,*) "-1 Exit"
    write(*,*) "0 Calculate radial distribution function and its integration curve"
    write(*,"(a,i4)") " 1 Select real space function, current:",ifunc
    write(*,"(a,3f10.4,' Ang')") " 2 Set sphere center, current",cenx*b2a,ceny*b2a,cenz*b2a
    write(*,"(a,2f9.4,' Ang')") " 3 Set lower and upper limit of radial distance, current:",rlow*b2a,rhigh*b2a
    write(*,"(a,i6)") " 4 Set the number of integration point in each shell, current:",nsphpt
    write(*,"(a,i6)") " 5 Set the number of radial points, current:",nradpt
    read(*,*) isel
    if (isel==-1) then
        return
    else if (isel==1) then
        call selfunc_interface(ifunc)
    else if (isel==2) then
        write(*,*) "Input sphere center (Angstrom), e.g. 0.0,1.2,-0.4"
        read(*,*) cenx,ceny,cenz
        cenx=cenx/b2a !Convert to Bohr
        ceny=ceny/b2a
        cenz=cenz/b2a
    else if (isel==3) then
        write(*,*) "Input lower and upper limit (Angstrom), e.g. 0.0,8.0"
        read(*,*) rlow,rhigh
        rlow=rlow/b2a !Convert to Bohr
        rhigh=rhigh/b2a
    else if (isel==4) then
        write(*,"(a)") " Input the number of integration point in each shell, the value must be one of &
        110/170/230/266/302/434/590/770/974/1454/1730/2030/2354/2702/3074/3470/3890/4334/4802/5294/5810"
        read(*,*) nsphpt
    else if (isel==5) then
        write(*,*) "Input the number of radial points, e.g. 800"
        read(*,*) nradpt
    else if (isel==0) then
        allocate(potx(nsphpt),poty(nsphpt),potz(nsphpt),potw(nsphpt))
        allocate(radval(nradpt),radpos(nradpt),intradval(nradpt),sphavgval(nradpt)) !radval records RDF, radpot records r position, intradval records integration of RDF
        call Lebedevgen(nsphpt,potx,poty,potz,potw)
        radval=0D0
        radstp=(rhigh-rlow)/(nradpt-1)
        ifinish=0
        iprogstp=20
        iprogcrit=iprogstp
        write(*,*) "Calculating..."
nthreads=getNThreads()
nthreads=getNThreads()
!$OMP PARALLEL DO SHARED(radval,radpos,ifinish,iprogcrit) PRIVATE(irad,radnow,isph,rnowx,rnowy,rnowz,tmpval) schedule(dynamic) NUM_THREADS(nthreads)
        do irad=1,nradpt
            radnow=rlow+(irad-1)*radstp
            radpos(irad)=radnow
            tmpval=0
            do isph=1,nsphpt
                rnowx=potx(isph)*radnow+cenx
                rnowy=poty(isph)*radnow+ceny
                rnowz=potz(isph)*radnow+cenz
                tmpval=tmpval+calcfuncall(ifunc,rnowx,rnowy,rnowz)*potw(isph)
            end do
            radval(irad)=4*pi*tmpval*radnow**2 !Multiplied by 4*pi is because the Lebedev integration routine produces unity rather than 4*pi
            sphavgval(irad)=tmpval !Spherically average function
            ifinish=ifinish+1
            if (ifinish==iprogcrit) then
                write(*,"(' Finished:',i6,'  /',i6)") ifinish,nradpt
                iprogcrit=iprogcrit+iprogstp
            end if
        end do
!$OMP END PARALLEL DO
        
        !Calculate integration of RDF
        intradval(1)=0D0
        do irad=2,nradpt
            intradval(irad)=intradval(irad-1)+radval(irad-1)*radstp
        end do
        write(*,"(a,f22.10)") " Integrating the RDF in the specified range is",intradval(nradpt)
        valrange=maxval(radval)-minval(radval)
        valrangeint=maxval(intradval)-minval(intradval)
        ilenunit1D=2
        do while(.true.)
            write(*,*)
            if (ilenunit1D==1) write(*,*) "-1 Switch the length unit for plotting, current: Bohr"
            if (ilenunit1D==2) write(*,*) "-1 Switch the length unit for plotting, current: Angstrom"
            write(*,*) "0 Return"
            write(*,*) "1 Plot the radial distribution function"
            write(*,*) "2 Plot integration curve of the RDF"
            write(*,*) "3 Save the radial distribution function map in current folder"
            write(*,*) "4 Save integration curve of the RDF map in current folder"
            write(*,*) "5 Export the radial distribution function to RDF.txt in current folder"
            write(*,*) "6 Export integration curve of the RDF to intRDF.txt in current folder"
            write(*,*) "7 Export the spherically averaged function"
            read(*,*) isel
            if (isel==-1) then
                if (ilenunit1D==1) then
                    ilenunit1D=2
                else if (ilenunit1D==2) then
                    ilenunit1D=1
                end if
            else if (isel==0) then
                deallocate(potx,poty,potz,potw,radval,radpos,intradval,sphavgval)
                exit
            else if (isel==1.or.isel==3) then
                ylow=minval(radval)-0.1D0*valrange
                yhigh=maxval(radval)+0.1D0*valrange
                if (isel==1) then
                else
                    write(*,"(a,a,a)") " Graph have been saved to ",trim(graphformat)," file with ""DISLIN"" prefix in current directory"
                end if
            else if (isel==2.or.isel==4) then
                ylow=minval(intradval)-0.1D0*valrangeint
                yhigh=maxval(intradval)+0.1D0*valrangeint
                if (isel==2) then
                else if (isel==4) then
                    write(*,"(a,a,a)") " Graph have been saved to ",trim(graphformat)," file with ""DISLIN"" prefix in current directory"
                end if
            else if (isel==5) then
                open(10,file="RDF.txt",status="replace")
                do irad=1,nradpt
                    write(10,"(i7,f12.4,f22.10)") irad,radpos(irad)*b2a,radval(irad)
                end do
                close(10)
                write(*,*) "The result has been output to RDF.txt in current folder"
                write(*,*) "The second column is radial distance (Angstrom), the third column is value"
            else if (isel==6) then
                open(10,file="intRDF.txt",status="replace")
                do irad=1,nradpt
                    write(10,"(i7,f12.4,f22.10)") irad,radpos(irad)*b2a,intradval(irad)
                end do
                close(10)
                write(*,*) "The result has been output to intRDF.txt in current folder"
                write(*,*) "The second column is radial distance (Angstrom), the third column is value"
            else if (isel==7) then
                open(10,file="sphavgval.txt",status="replace")
                do irad=1,nradpt
                    write(10,"(i7,f12.6,f18.7)") irad,radpos(irad),sphavgval(irad)
                end do
                close(10)
                write(*,*) "The result has been output to sphavgval.txt in current folder"
                write(*,*) "The second column is radial distance (Bohr), the third column is value"
            end if
        end do
    
    end if
end do
end subroutine



!!--------- Analyze correspondence between orbitals in two wavefunctions
subroutine orbcorres
use defvar
use function
use util
implicit real*8 (a-h,o-z)
real*8,allocatable :: convmat(:,:) !(i,j) is the coefficient of j MO of the second wavefunction in i MO of current wavefunction
real*8,allocatable :: MOvalgrd(:,:),MOvalgrd2(:,:) !MOvalgrd(j,n),MOvalgrd2(j,n) means the the value of the nth MO of the first/second wavefunction at the jth grid
real*8,allocatable :: comparr(:),beckeweigrid(:)
integer,allocatable :: comparridx(:)
character filename2*200,c80tmp*80
type(content),allocatable :: gridatm(:),gridatmorg(:)
if (iautointgrid==1) then
    sphpotold=sphpot
    radpotold=radpot
    sphpot=230
    radpot=50
end if
allocate(gridatm(radpot*sphpot),gridatmorg(radpot*sphpot),beckeweigrid(radpot*sphpot))

do isep=nmo,1,-1
    if (MOtype(isep)==1) exit
end do
if (wfntype==1.or.wfntype==4) write(*,"(' Note: The orbitals from',i6,' to',i6,' are alpha; from',i6,' to',i6,' are beta')") 1,isep,isep+1,nmo
write(*,"(a)") " Input the range of the orbitals of present wavefunction to be considered, e.g. 2,9. &
If press ENTER directly, all orbitals will be taken into account"
read(*,"(a)") c80tmp
if (c80tmp==" ") then
    istart1=1
    iend1=nmo
else
    read(c80tmp,*) istart1,iend1
end if

write(*,*)
write(*,*) "Input filename of the second wavefunction, e.g. c:\ltwd.fch"
do while(.true.)
    read(*,"(a)") filename2
    inquire(file=filename2,exist=alive)
    if (alive) exit
    write(*,*) "Cannot find the file, input again"
end do
call dealloall
call readinfile(filename2,1) !Get some knowledge about the second wavefunction
nmo2=nmo !The number of MOs in the wfn2
iwfntype2=wfntype
do isep=nmo2,1,-1
    if (MOtype(isep)==1) exit
end do
write(*,*)
if (iwfntype2==1.or.iwfntype2==4) write(*,"(' Note: The orbitals from',i6,' to',i6,' are alpha; from',i6,' to',i6,' are beta')") 1,isep,isep+1,nmo
write(*,"(a)") " Input the range of the orbitals of the second wavefunction to be considered, e.g. 2,9. &
If press ENTER directly, all orbitals will be taken into account"
read(*,"(a)") c80tmp
if (c80tmp==" ") then
    istart2=1
    iend2=nmo2
else
    read(c80tmp,*) istart2,iend2
end if

call dealloall
call readinfile(firstfilename,1)
allocate(MOvalgrd(radpot*sphpot,nmo),MOvalgrd2(radpot*sphpot,nmo2),convmat(nmo,nmo2))
convmat=0D0

write(*,"(' Radial points:',i5,'    Angular points:',i5,'   Total:',i10,' per center')") radpot,sphpot,radpot*sphpot
write(*,*) "Calculating, please wait..."
call gen1cintgrid(gridatmorg,iradcut)

call walltime(iwalltime1)
CALL CPU_TIME(time_begin)

do iatm=1,ncenter
    write(*,"(' Progress: ',i5,' /',i5)") iatm,ncenter
    gridatm%x=gridatmorg%x+a(iatm)%x !Move quadrature point to actual position in molecule
    gridatm%y=gridatmorg%y+a(iatm)%y
    gridatm%z=gridatmorg%z+a(iatm)%z
    
    !Calculate value of all MOs of the first and second wavefunction at all grids
    call dealloall
    call readinfile(filename2,1) !Load wfn2
nthreads=getNThreads()
!$OMP parallel do shared(MOvalgrd2) private(ipt) num_threads(nthreads) schedule(dynamic)
    do ipt=1+iradcut*sphpot,radpot*sphpot
        call orbderv(1,istart2,iend2,gridatm(ipt)%x,gridatm(ipt)%y,gridatm(ipt)%z,MOvalgrd2(ipt,:))
    end do
!$OMP end parallel do
    call dealloall
    call readinfile(firstfilename,1) !Retrieve to wfn1
nthreads=getNThreads()
!$OMP parallel do shared(MOvalgrd) private(ipt) num_threads(nthreads) schedule(dynamic)
    do ipt=1+iradcut*sphpot,radpot*sphpot
        call orbderv(1,istart1,iend1,gridatm(ipt)%x,gridatm(ipt)%y,gridatm(ipt)%z,MOvalgrd(ipt,:))
    end do
!$OMP end parallel do

    !Calculate Becke weight at all grids
    call gen1cbeckewei(iatm,iradcut,gridatm,beckeweigrid)
    
nthreads=getNThreads()
!$OMP parallel do shared(convmat) private(imo,jmo,tmpval,ipt) num_threads(nthreads) schedule(dynamic)
    do imo=istart1,iend1
        do jmo=istart2,iend2
            tmpval=0D0
            do ipt=1+iradcut*sphpot,radpot*sphpot
                tmpval=tmpval+beckeweigrid(ipt)*gridatmorg(ipt)%value*MOvalgrd(ipt,imo)*MOvalgrd2(ipt,jmo)
            end do
            convmat(imo,jmo)=convmat(imo,jmo)+tmpval
        end do
    end do
!$OMP end parallel do
end do
CALL CPU_TIME(time_end)
call walltime(iwalltime2)
write(*,"(' Calculation took up CPU time',f12.2,'s, wall clock time',i10,'s',/)") time_end-time_begin,iwalltime2-iwalltime1

! call showmatgau(convmat,"convmat",0,"f12.3")

devmax=0D0
idevmax=1
! write(*,*) "The sum of composition of each orbital of current wavefunction"
if ((wfntype==0.or.wfntype==2).and.(iwfntype2==0.or.iwfntype2==2)) then !Both of the two wavefunctions are R or RO types
    do imo=istart1,iend1 !Check normalization
        totcomp=sum(convmat(imo,:)**2)*100D0
    !     write(*,"(i6,':',f12.5,' %')") imo,totcomp
        if (abs(totcomp-100D0)>devmax) then
            devmax=abs(totcomp-100D0)
            idevmax=imo
        end if
    end do
end if
write(*,"(' The maximum deviation to normalization condition is',f8.3,' % (Orbital',i6,')')") devmax,idevmax
write(*,"(a)") " Note: The first column below is the index of the orbitals in present wavefunction, the largest five contributions from &
the orbitals in the second wavefunction are shown at right side. If the dominative index is inconsistent to the first column, the row will be marked by asterisk"
write(*,*)
allocate(comparr(nmo2),comparridx(nmo2))
do imo=istart1,iend1
    !Sort the composition array from small to large
    comparr(:)=convmat(imo,:)
    do itmp=1,nmo2
        comparridx(itmp)=itmp
    end do
    do i=istart2,iend2
        do j=i+1,iend2
            if (abs(comparr(i))>abs(comparr(j))) then
                temp=comparr(i)
                comparr(i)=comparr(j)
                comparr(j)=temp
                itemp=comparridx(i)
                comparridx(i)=comparridx(j)
                comparridx(j)=itemp
            end if
        end do
    end do
    if (comparridx(iend2)/=imo) then
        write(*,"('*',i5,':  ')",advance="no") imo
    else
        write(*,"(' ',i5,':  ')",advance="no") imo
    end if
    do i=iend2,iend2-4,-1
        write(*,"(i5,'(',f6.2,'%)')",advance="no") comparridx(i),comparr(i)**2*100D0
    end do
    write(*,*)
end do

write(*,*)
do while(.true.)
    write(*,*) "Input the orbital index to print detail compositions and coefficients, e.g. 5"
    write(*,*) "input -1 can output all overlap integrals between the chosen orbitals"
    write(*,*) "Input 0 can exit"
    read(*,*) imo
    if (imo==0) then
        exit
    else if (imo==-1) then
        open(10,file="convmat.txt",status="replace")
        do i=istart1,iend1
            do j=istart2,iend2
                write(10,"(2i7,f12.6)") i,j,convmat(i,j)
            end do
        end do
        close(10)
        write(*,*) "The overlap integrals have been outputted to convmat.txt in current folder"
        write(*,"(a,/)") " The first and second columns correspond to the orbital indices in present and in the second wavefunctions"
    else if (imo<istart1.or.imo>iend1) then
        write(*,"(a,i6,a,i6)") "Error: Exceed valid range! The value should within",istart1," and",iend1
    else
        tmpval=0D0
        do jmo=istart2,iend2
            tmpval=tmpval+convmat(imo,jmo)**2*100D0
            write(*,"(i6,'   Contribution:',f10.3,' %    Coefficient:',f12.6)") jmo,convmat(imo,jmo)**2*100D0,convmat(imo,jmo)
        end do
        write(*,"(' Total:',f10.3,' %')") tmpval
        write(*,*)
    end if
end do

if (iautointgrid==1) then
    radpot=radpotold
    sphpot=sphpotold
end if
end subroutine





!!-------- Parse the output of (hyper)polarizability task of Gaussian to make it more understandable
subroutine parseGauPolar
use defvar
use util
implicit real*8 (a-h,o-z)
real*8 dipole(3),poltens(3,3),hypoltens(3,3,3) !hypoltens2(3,3,3,3)
real*8 eigvecmat(3,3),eigval(3),freqval(100000)
character c200tmp*200,sepchar,c210tmp*210
character*20 :: form,formau="(a,f16.6)",formother="(a,1PE16.6)"
integer :: irdfreq=0,ides=6,iunit=1
poltens=0D0
hypoltens=0D0
write(*,*) "Note: This function only works for ""polar"" tasks of Gaussian 09 with #P"
do while(.true.)
    write(*,"(a,/)") " Select the type of your Gaussian task by option 1~6. 2,4,6 only give polarizability (alpha), &
    the others also give hyperpolarizability (beta). -1 can be chosen only if CPHF=RdFreq or polar=DCSHG is used"
    if (iunit==1) write(*,*) "-3 Set the unit in the output, current: a.u."
    if (iunit==2) write(*,*) "-3 Set the unit in the output, current: SI"
    if (iunit==3) write(*,*) "-3 Set the unit in the output, current: esu"
    if (ides==6) write(*,*) "-2 Set the output destination, current: Output to screen"
    if (ides==11) write(*,*) "-2 Set the output destination, current: polar.txt in current folder"
    if (irdfreq==1) write(*,*) "-1 Toggle if load frequency-dependent result for option 1, current: Yes"
    if (irdfreq==0) write(*,*) "-1 Toggle if load frequency-dependent result for option 1, current: No"
    write(*,*) "0 Return"
    write(*,*) "1 ""Polar"" + analytic 3-order deriv. (HF/DFT/Semi-empirical)"
    write(*,*) "2 ""Polar"" + analytic 2-order deriv. (MP2...)"
    write(*,*) "3 ""Polar=Cubic"" + analytic 2-order deriv."
    write(*,*) "4 ""Polar"" + analytic 1-order deriv. (CISD,QCISD,CCSD,MP3,MP4(SDQ)...)"
    write(*,*) "5 ""Polar=DoubleNumer"" or ""Polar=EnOnly"" + analytic 1-order deriv."
    write(*,*) "6 ""Polar"" + energy only (CCSD(T),QCISD(T),MP4(SDTQ),MP5...)"
    read(*,*) isel
    if (isel==-1) then
        if (irdfreq==1) then
            irdfreq=0
        else
            irdfreq=1
        end if
    else if (isel==-2) then
        write(*,*) "6: Output to screen"
        write(*,*) "11: Output to polar.txt in current folder"
        read(*,*) ides
    else if (isel==-3) then
        write(*,*) "1: Atomic unit"
        write(*,*) "2: SI unit (C^2*m^2/J for alpha, C^3*m^3/J for beta)"
        write(*,*) "3: esu"
        read(*,*) iunit
    else if (isel==0) then
        return
    else
        exit
    end if
end do
if (irdfreq==1.and.isel>=2) then
    write(*,*) "ERROR: Frequency-dependent values are only available for HF/DFT/semi-empirical!"
    return
end if

open(10,file=filename,status="old")
if (ides==11) open(ides,file="polar.txt",status="replace")

if (iunit==1) then
    write(ides,*) "Note: All units shown below are in a.u."
    form=formau
else if (iunit==2) then
    write(ides,*) "Note: All units shown below are in SI unit (C^2*m^2/J for alpha, C^3*m^3/J for beta)"
    form=formother
else if (iunit==3) then
    write(ides,*) "Note: All units shown below are in esu unit"
    form=formother
end if
write(ides,*)


! Dipole moment part (miu)
call loclabel(10,"Dipole moment (field-independent basis, Debye)",ifound)
if (ifound==1) then
    read(10,*)
    read(10,*) c200tmp,xtmp,c200tmp,ytmp,c200tmp,ztmp
    dipole(1)=xtmp/au2debye !Convert to a.u.
    dipole(2)=ytmp/au2debye
    dipole(3)=ztmp/au2debye
    if (iunit==2) dipole=dipole*8.47835D-30
    if (iunit==3) dipole=dipole*2.54175D-18    
    dipnorm=dsqrt(sum(dipole**2))
    write(ides,*) "Dipole moment:"
    if (iunit==1) then
        write(ides,"(' X,Y,Z=',3f12.6,'   Norm=',f12.6)") dipole(:),dipnorm
    else
        write(ides,"(' X,Y,Z=',3(1PE15.6),'   Norm=',1PE15.6)") dipole(:),dipnorm
    end if
    write(ides,*)
end if

!Selecting the result at which frequency will be loaded
if (irdfreq==1) then
    nfreqval=0
    rewind(10)
    do while(.true.)
        call loclabel(10,"-- Alpha(-w,w) frequency",ifound,0) !Not rewind
        if (ifound==1) then
            read(10,"(46x,f12.6)") tmpval
            if ( any(freqval(1:nfreqval)==tmpval) ) exit !The output in stage 2 and stage 3 are identical, so determine if has entered stage 3
            nfreqval=nfreqval+1
            freqval(nfreqval)=tmpval
        else
            exit
        end if
    end do
    write(*,*) "Load which one? Input the index"
    do i=1,nfreqval
        if (freqval(i)>0) write(*,"(i8,'   w=',f12.6,' (',f12.2,'nm )')") i,freqval(i),1240.7011D0/(freqval(i)*au2eV)
        if (freqval(i)==0) write(*,"(i8,'   w=',f12.6,' (     Static    )')") i,freqval(i)
    end do
    read(*,*) ifreq
    if (freqval(ifreq)>0) write(*,"(' Note: All Alpha and Beta below correspond to w=',f12.6,' (',f12.2,'nm )',/)") freqval(ifreq),1240.7011D0/(freqval(ifreq)*au2eV)
    if (freqval(ifreq)==0) write(*,"(' Note: All Alpha and Beta below correspond to w=',f12.6,' ( Static )',/)") freqval(ifreq)
    rewind(10) !Move to the beginning
end if


!!!! Polarizability part (Alpha)
if (isel==1.or.isel==4.or.isel==5) then
    if (isel==1) then
        if (irdfreq==0) then 
            call loclabel(10,"SCF Polarizability",ifound)
        else if (irdfreq==1) then
            do i=1,ifreq
                call loclabel(10,"SCF Polarizability",ifound,0) !Not rewind
                read(10,*)
            end do
            backspace(10)
        end if
    else if (isel==4.or.isel==5) then
        call loclabel(10,"Isotropic polarizability",ifound)
    end if
    call skiplines(10,2)
    read(10,*) rnouse,poltens(1,1)
    read(10,*) rnouse,poltens(2,1),poltens(2,2)
    read(10,*) rnouse,poltens(3,1),poltens(3,2),poltens(3,3)
else if (isel==2.or.isel==3) then
    call loclabel(10,"Exact polarizability")
    read(10,"(23x,6f8.3)") poltens(1,1),poltens(2,1),poltens(2,2),poltens(3,1),poltens(3,2),poltens(3,3)
else if (isel==6) then !Find result from archive part
    sepchar="\"
    if (isys==1) sepchar="|"
    do while(.true.)
        read(10,"(1x,a70)") c210tmp(1:70) !Combine two lines
        read(10,"(1x,a70)") c210tmp(71:140)
        read(10,"(1x,a70)") c210tmp(141:210)
        if (index(c210tmp(1:140),sepchar//"Polar=")/=0) exit
        backspace(10)
        backspace(10)
    end do
    do istart=1,210
        if (c210tmp(istart:istart+6)==sepchar//"Polar=") exit
    end do
    do iend=istart+1,210
        if (c210tmp(iend:iend)=="|".or.c200tmp(iend:iend)=="\") exit
    end do
    c210tmp(1:istart+6)=" " !Clean other information so that the data can be read in free format
    c210tmp(iend:)=" "
    read(c210tmp,*) poltens(1,1),poltens(2,1),poltens(2,2),poltens(3,1),poltens(3,2),poltens(3,3)
end if
poltens(1,2)=poltens(2,1)
poltens(1,3)=poltens(3,1)
poltens(2,3)=poltens(3,2)
!Convert to other unit
if (iunit==2) poltens=poltens*1.6488D-41
if (iunit==3) poltens=poltens*1.4819D-25
if (irdfreq==0) write(ides,*) "Static polarizability:"
if (irdfreq==1) write(ides,*) "Frequency-dependent polarizability:"
write(ides,form) " XX=",poltens(1,1)
write(ides,form) " XY=",poltens(1,2)
write(ides,form) " YY=",poltens(2,2)
write(ides,form) " XZ=",poltens(1,3)
write(ides,form) " YZ=",poltens(2,3)
write(ides,form) " ZZ=",poltens(3,3)
alphaiso=(poltens(1,1)+poltens(2,2)+poltens(3,3))/3D0
write(ides,form) ' Isotropic average polarizability:',alphaiso
term1=(poltens(1,1)-poltens(2,2))**2 + (poltens(1,1)-poltens(3,3))**2 + (poltens(2,2)-poltens(3,3))**2
term2=6*(poltens(1,2)**2+poltens(1,3)**2+poltens(2,3)**2)
alphaani1=dsqrt((term1+term2)/2D0)
write(ides,form) ' Polarizability anisotropy (definition 1):',alphaani1
alphaani2=dsqrt(term1/2D0)
write(ides,form) ' Polarizability anisotropy (definition 2):',alphaani2
call diagmat(poltens,eigvecmat,eigval,300,1D-10)
call sortr8(eigval)
if (iunit==1) then
    write(ides,"(a,3f13.5)") ' Eigenvalues of polarizability tensor:',eigval
else
    write(ides,"(a,3(1PE13.5))") ' Eigenvalues of polarizability tensor:',eigval
end if
alphaani3=eigval(3)-(eigval(1)+eigval(2))/2D0
write(ides,form) ' Polarizability anisotropy (definition 3):',alphaani3
write(ides,*)


!!!! First hyperpolarizability part (Beta)
if (isel==1.or.isel==3.or.isel==5) then
    if (irdfreq==0) then
        if (isel==1) then
            call loclabel(10,"SCF Static Hyperpolarizability",ifound)
        else if (isel==3) then
            call loclabel(10,"Final packed hyperpolarizability",ifound)
        else if (isel==5) then
            call loclabel(10,"Static Hyperpolarizability",ifound)
        end if
        call skiplines(10,3)
        read(10,*) rnouse,hypoltens(1,1,1) !XXX
        call skiplines(10,2)
        read(10,*) rnouse,hypoltens(1,1,2) !XXY
        read(10,*) rnouse,hypoltens(1,2,2),hypoltens(2,2,2) !XYY,YYY
        call skiplines(10,2)
        read(10,*) rnouse,hypoltens(1,1,3) !XXZ
        read(10,*) rnouse,hypoltens(1,2,3),hypoltens(2,2,3) !XYZ,YYZ
        read(10,*) rnouse,hypoltens(1,3,3),hypoltens(2,3,3),hypoltens(3,3,3) !XZZ,YZZ,ZZZ
        !XYX=YXX    =XXY
        hypoltens(1,2,1)=hypoltens(1,1,2)
        hypoltens(2,1,1)=hypoltens(1,1,2)
        !YXY=YYX    =XYY
        hypoltens(2,1,2)=hypoltens(1,2,2)
        hypoltens(2,2,1)=hypoltens(1,2,2)
        !ZXX=XZX    =XXZ
        hypoltens(3,1,1)=hypoltens(1,1,3)
        hypoltens(1,3,1)=hypoltens(1,1,3)
        !XZY=YXZ=ZYX=YXZ=YZX   =XYZ
        hypoltens(1,3,2)=hypoltens(1,2,3)
        hypoltens(2,1,3)=hypoltens(1,2,3)
        hypoltens(3,2,1)=hypoltens(1,2,3)
        hypoltens(2,1,3)=hypoltens(1,2,3)
        hypoltens(2,3,1)=hypoltens(1,2,3)
        !ZYY=YZY    =YYZ
        hypoltens(3,2,2)=hypoltens(2,2,3)
        hypoltens(2,3,2)=hypoltens(2,2,3)
        !ZXZ,ZZX    =XZZ
        hypoltens(3,1,3)=hypoltens(1,3,3)
        hypoltens(3,3,1)=hypoltens(1,3,3)
        !ZYZ=ZZY    =YZZ
        hypoltens(3,2,3)=hypoltens(2,3,3)
        hypoltens(3,3,2)=hypoltens(2,3,3)
        write(ides,"(a,/)") " Note: It is well known that the sign of hyperpolarizability of Gaussian 09 should be multiplied by -1, the outputs below have already been corrected"
        hypoltens=-hypoltens
        !Convert to other unit
        if (iunit==2) hypoltens=hypoltens*3.20636D-53
        if (iunit==3) hypoltens=hypoltens*8.63922D-33
        write(ides,*) "Static first hyperpolarizability:"
        write(ides,form) " XXX=",hypoltens(1,1,1)
        write(ides,form) " XXY=",hypoltens(1,1,2)
        write(ides,form) " XYY=",hypoltens(1,2,2)
        write(ides,form) " YYY=",hypoltens(2,2,2)
        write(ides,form) " XXZ=",hypoltens(1,1,3)
        write(ides,form) " XYZ=",hypoltens(1,2,3)
        write(ides,form) " YYZ=",hypoltens(2,2,3)
        write(ides,form) " XZZ=",hypoltens(1,3,3)
        write(ides,form) " YZZ=",hypoltens(2,3,3)
        write(ides,form) " ZZZ=",hypoltens(3,3,3)
        write(ides,*)
        betaX=hypoltens(1,1,1)+hypoltens(1,2,2)+hypoltens(1,3,3)
        betaY=hypoltens(2,2,2)+hypoltens(2,1,1)+hypoltens(2,3,3)
        betaZ=hypoltens(3,3,3)+hypoltens(3,2,2)+hypoltens(3,1,1)
        if (iunit==1) then
            write(ides,"(' Beta_X=',f15.5,'  Beta_Y=',f15.5,'  Beta_Z=',f15.5)") betaX,betaY,betaZ
        else
            write(ides,"(' Beta_X=',1PE15.5,'  Beta_Y=',1PE15.5,'  Beta_Z=',1PE15.5)") betaX,betaY,betaZ
        end if
        write(ides,form) " Magnitude of first hyperpolarizability:",dsqrt(betaX**2+betaY**2+betaZ**2)
        betaprj=(betaX*dipole(1)+betaY*dipole(2)+betaZ*dipole(3))/dipnorm
        write(ides,form) " Projection of beta on dipole moment:",betaprj
        write(ides,form) " Beta ||     :",betaprj/5D0*3D0
        write(ides,form) " Beta ||(z)  :",betaZ/5D0*3D0
        write(ides,form) " Beta _|_(z) :",betaZ/5D0
        
        !---For easier inspection in special use
!         write(*,*)
!         write(*,*) "============="
!         write(ides,form) ' Isotropic average polarizability:',alphaiso
!         write(ides,form) ' Polarizability anisotropy (definition 2):',alphaani2
!         write(ides,form) " Magnitude of first hyperpolarizability:",dsqrt(betaX**2+betaY**2+betaZ**2)
!         write(ides,form) " Beta ||     :",betaprj/5D0*3D0
!         pause
        
    else if (irdfreq==1) then !Frequency-dependent hyperpolarizability, only available for HF/DFT/semi-empirical
        write(*,*) "Loading which type of hyperpolarizability?"
        write(*,*) "1: Beta(-w;w,0)   2: Beta(-2w;w,w)"
        write(*,*) "Note: Option 2 is meaningless if ""DCSHG"" keyword was not used" 
        read(*,*) ibeta
        rewind(10)
        
        if (ibeta==1) then !Beta(-w;w,0) case
            call loclabel(10,"-- Beta(-w,w,0) frequency",ifound)
            do i=1,ifreq
                call loclabel(10,"-- Beta(-w,w,0) frequency",ifound,0)
                read(10,*)
            end do
            read(10,*)
            read(10,*) c200tmp,hypoltens(1,1,1)
            read(10,*) c200tmp,hypoltens(2,1,1)
            read(10,*) c200tmp,hypoltens(2,2,1)
            read(10,*) c200tmp,hypoltens(3,1,1)
            read(10,*) c200tmp,hypoltens(3,2,1)
            read(10,*) c200tmp,hypoltens(3,3,1)
            hypoltens(1,2,1)=hypoltens(2,1,1)
            hypoltens(1,3,1)=hypoltens(3,1,1)
            hypoltens(2,3,1)=hypoltens(3,2,1)
            read(10,*) c200tmp,hypoltens(1,1,2)
            read(10,*) c200tmp,hypoltens(2,1,2)
            read(10,*) c200tmp,hypoltens(2,2,2)
            read(10,*) c200tmp,hypoltens(3,1,2)
            read(10,*) c200tmp,hypoltens(3,2,2)
            read(10,*) c200tmp,hypoltens(3,3,2)
            hypoltens(1,2,2)=hypoltens(2,1,2)
            hypoltens(1,3,2)=hypoltens(3,1,2)
            hypoltens(2,3,2)=hypoltens(3,2,2)
            read(10,*) c200tmp,hypoltens(1,1,3)
            read(10,*) c200tmp,hypoltens(2,1,3)
            read(10,*) c200tmp,hypoltens(2,2,3)
            read(10,*) c200tmp,hypoltens(3,1,3)
            read(10,*) c200tmp,hypoltens(3,2,3)
            read(10,*) c200tmp,hypoltens(3,3,3)
            hypoltens(1,2,3)=hypoltens(2,1,3)
            hypoltens(1,3,3)=hypoltens(3,1,3)
            hypoltens(2,3,3)=hypoltens(3,2,3)
            write(ides,"(a,/)") " Note: It is well known that the sign of hyperpolarizability of Gaussian 09 should be multiplied by -1, the outputs below have already been corrected"
            hypoltens=-hypoltens
            !Convert to other unit
            if (iunit==2) hypoltens=hypoltens*3.20636D-53
            if (iunit==3) hypoltens=hypoltens*8.63922D-33
            write(ides,*) "Frequency-dependent first hyperpolarizability Beta(-w;w,0)"
            write(ides,form) " XXX=     ",hypoltens(1,1,1)
            write(ides,form) " XYX= YXX=",hypoltens(1,2,1)
            write(ides,form) " YYX=     ",hypoltens(2,2,1)
            write(ides,form) " XZX= ZXX=",hypoltens(1,3,1)
            write(ides,form) " YZX= ZYX=",hypoltens(2,3,1)
            write(ides,form) " ZZX=     ",hypoltens(3,3,1)
            write(ides,form) " XXY=     ",hypoltens(1,1,2)
            write(ides,form) " XYY= YXY=",hypoltens(1,2,2)
            write(ides,form) " YYY=     ",hypoltens(2,2,2)
            write(ides,form) " XZY= ZXY=",hypoltens(1,3,2)
            write(ides,form) " YZY= ZYY=",hypoltens(2,3,2)
            write(ides,form) " ZZY=     ",hypoltens(3,3,2)
            write(ides,form) " XXZ=     ",hypoltens(1,1,3)
            write(ides,form) " XYZ= YXZ=",hypoltens(1,2,3)
            write(ides,form) " YYZ=     ",hypoltens(2,2,3)
            write(ides,form) " XZZ= ZXZ=",hypoltens(1,3,3)
            write(ides,form) " YZZ= ZYZ=",hypoltens(2,3,3)
            write(ides,form) " ZZZ=     ",hypoltens(3,3,3)
            write(ides,*)
            
        else if (ibeta==2) then !used DCSHG, parsing Beta(-2w;w,w)
            call loclabel(10,"-- Beta(w,w,-2w) frequency",ifound)
            do i=1,ifreq
                call loclabel(10,"-- Beta(w,w,-2w) frequency",ifound,0)
                read(10,*)
            end do
            read(10,*)
            read(10,*) c200tmp,hypoltens(1,1,1)
            read(10,*) c200tmp,hypoltens(1,1,2)
            hypoltens(1,2,1)=hypoltens(1,1,2)
            read(10,*) c200tmp,hypoltens(1,2,2)
            read(10,*) c200tmp,hypoltens(1,1,3)
            hypoltens(1,3,1)=hypoltens(1,1,3)
            read(10,*) c200tmp,hypoltens(1,2,3)
            hypoltens(1,3,2)=hypoltens(1,2,3)
            read(10,*) c200tmp,hypoltens(1,3,3)
            read(10,*) c200tmp,hypoltens(2,1,1)
            read(10,*) c200tmp,hypoltens(2,2,1)
            hypoltens(2,1,2)=hypoltens(2,2,1)
            read(10,*) c200tmp,hypoltens(2,2,2)
            read(10,*) c200tmp,hypoltens(2,1,3)
            hypoltens(2,3,1)=hypoltens(2,1,3)
            read(10,*) c200tmp,hypoltens(2,2,3)
            hypoltens(2,3,2)=hypoltens(2,2,3)
            read(10,*) c200tmp,hypoltens(2,3,3)
            read(10,*) c200tmp,hypoltens(3,1,1)
            read(10,*) c200tmp,hypoltens(3,1,2)
            hypoltens(3,2,1)=hypoltens(3,1,2)
            read(10,*) c200tmp,hypoltens(3,2,2)
            read(10,*) c200tmp,hypoltens(3,1,3)
            hypoltens(3,3,1)=hypoltens(3,1,3)
            read(10,*) c200tmp,hypoltens(3,2,3)
            hypoltens(3,3,2)=hypoltens(3,2,3)
            read(10,*) c200tmp,hypoltens(3,3,3)
            write(ides,"(a,/)") " Note: It is well known that the sign of hyperpolarizability of Gaussian 09 should be multiplied by -1, the outputs below have already been corrected"
            hypoltens=-hypoltens
            !Convert to other unit
            if (iunit==2) hypoltens=hypoltens*3.20636D-53
            if (iunit==3) hypoltens=hypoltens*8.63922D-33
            if (ibeta==2) write(ides,*) "Frequency-dependent first hyperpolarizability Beta(-2w;w,w)"
            write(ides,form) " XXX=     ",hypoltens(1,1,1)
            write(ides,form) " YXX=     ",hypoltens(2,1,1)
            write(ides,form) " ZXX=     ",hypoltens(3,1,1)
            write(ides,form) " XYX= XXY=",hypoltens(1,2,1)
            write(ides,form) " YYX= YXY=",hypoltens(2,2,1)
            write(ides,form) " ZYX= ZXY=",hypoltens(3,2,1)
            write(ides,form) " XYY=     ",hypoltens(1,2,2)
            write(ides,form) " YYY=     ",hypoltens(2,2,2)
            write(ides,form) " ZYY=     ",hypoltens(3,2,2)
            write(ides,form) " XZX= XXZ=",hypoltens(1,3,1)
            write(ides,form) " YZX= YXZ=",hypoltens(2,3,1)
            write(ides,form) " ZZX= ZXZ=",hypoltens(3,3,1)
            write(ides,form) " XZY= XYZ=",hypoltens(1,3,2)
            write(ides,form) " YZY= YYZ=",hypoltens(2,3,2)
            write(ides,form) " ZZY= ZYZ=",hypoltens(3,3,2)
            write(ides,form) " XZZ=     ",hypoltens(1,3,3)
            write(ides,form) " YZZ=     ",hypoltens(2,3,3)
            write(ides,form) " ZZZ=     ",hypoltens(3,3,3)
            write(ides,*)
        end if
        
        betaX=0
        betaY=0
        betaZ=0
        do j=1,3
            betaX=betaX+(hypoltens(1,j,j)+hypoltens(j,j,1)+hypoltens(j,1,j))/3
            betaY=betaY+(hypoltens(2,j,j)+hypoltens(j,j,2)+hypoltens(j,2,j))/3
            betaZ=betaZ+(hypoltens(3,j,j)+hypoltens(j,j,3)+hypoltens(j,3,j))/3
        end do
        if (iunit==1) then
            write(ides,"(' Beta_X=',f15.5,'  Beta_Y=',f15.5,'  Beta_Z=',f15.5)") betaX,betaY,betaZ
        else
            write(ides,"(' Beta_X=',1PE15.5,'  Beta_Y=',1PE15.5,'  Beta_Z=',1PE15.5)") betaX,betaY,betaZ
        end if
        write(ides,form) " Magnitude of first hyperpolarizability:",dsqrt(betaX**2+betaY**2+betaZ**2)
        betaprj=(betaX*dipole(1)+betaY*dipole(2)+betaZ*dipole(3))/dipnorm
        write(ides,form) " Projection of beta on dipole moment:",betaprj
        write(ides,form) " Beta ||     :",betaprj*3D0/5D0
        write(ides,form) " Beta ||(z)  :",betaZ*3D0/5D0
        beta_per=0
        do j=1,3
            beta_per=beta_per+(2*hypoltens(3,j,j)+2*hypoltens(j,j,3)-3*hypoltens(j,3,j))/5
        end do
        write(ides,form) " Beta _|_(z) :",beta_per
        write(*,*)
    end if
end if

close(10)
if (ides==11) close(ides)
end subroutine



!!--------- Sum-over-states (SOS) calculation for (hyper)polarizability
!Programmed based on the formulae in Sasagane et al. J. Chem. Phys., 99, 3738 (1993)
subroutine SOS
use defvar
use util
implicit real*8 (a-h,o-z)
character transmodestr*80,c80tmp*80,c200tmp*80
character :: dirlab(3)=(/ "X","Y","Z" /)
real*8,allocatable :: trandip(:,:,:) !Transition dipole moment between i and j in X,Y,Z. 0 corresponds to ground state
real*8,allocatable :: excene(:) !Excitation energy
real*8 :: alpha(3,3),beta(3,3,3),gamma(3,3,3,3),delta(3,3,3,3,3)
real*8 eigval(3),eigvecmat(3,3),tmpw(5)
real*8,allocatable :: freqlist(:,:) !Store the frequency to be calculated for beta and gamma
integer tmpdir(5),arrb(6,3),arrg(24,4),arrd(120,5),dir1,dir2,dir3,dir4,dir5

write(*,*) "Loading data..."
open(10,file=filename,status="old")
call loclabel(10,"Gaussian",igauout)
rewind(10)
if (igauout==1) then !Load excitation energies and <0|r|n>
    write(*,*) "This is a Gaussian output file"
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
    allocate(trandip(0:nstates,0:nstates,3),excene(0:nstates))
    trandip=0
    excene=0
    call loclabel(10,"Ground to excited state transition electric dipole moments")
    read(10,*)
    read(10,*)
    do istat=1,nstates
        read(10,*) inouse,trandip(0,istat,:) !Transition dipole moment from ground state (0) to excited "istate"
    end do
    call loclabel(10,"Excitation Energies [eV] at current iteration:",ifound)
    if (ifound==1) then !Read the excitation energies shown in the last step of iteration process will be more accurate, but #P must be used
        ncyc=1
        do while(.true.) !Find the position of the last iteration
            call loclabel(10,"Excitation Energies [eV] at current iteration:",ifound,0)
            if (ifound==0) exit
            read(10,*)
            ncyc=ncyc+1
        end do
        rewind(10)
        do icyc=1,ncyc
            call loclabel(10,"Excitation Energies [eV] at current iteration:",ifound,0)
        end do
        read(10,*)
        do istat=1,nstates
            read(10,"(14x)",advance="no")
            read(10,*) excene(istat)
        end do
    else !Read excitation energies from normal output
        call loclabel(10,"Excitation energies and oscillator strengths:")
        do istat=1,nstates
            call loclabel(10,"Excited State",ifound,0)
            read(10,"(a)") transmodestr
            do i=10,70
                if (transmodestr(i:i+1)=="eV") exit
            end do
            read(transmodestr(i-10:i-1),*) excene(istat) !Read as eV
        end do
    end if
    call loclabel(10,"Dipole moment (field-independent basis, Debye)")
    read(10,*)
    read(10,*) c80tmp,xtmp,c80tmp,ytmp,c80tmp,ztmp
    trandip(0,0,1)=xtmp/au2debye
    trandip(0,0,2)=ytmp/au2debye
    trandip(0,0,3)=ztmp/au2debye
    ionlyalpha=1
else !Load excitation energies and all <m|r|n> from plain text file
    ionlyalpha=0
    read(10,*) nstates
    if (nstates<0) ionlyalpha=1 !Only calculate polarizability, not hyperpolarizability, so will not read transition dipole moments among excited states
    nstates=abs(nstates)
    allocate(trandip(0:nstates,0:nstates,3),excene(0:nstates))
    do i=1,nstates !Read as eV
        read(10,*) inouse,excene(i)
    end do
    do i=0,nstates !i-i corresponds to dipole moment of corresponding state
        do j=i,nstates
            read(10,*) inouse,inouse,trandip(i,j,:)
        end do
        if (ionlyalpha==1) exit
    end do
end if
close(10)
excene(0)=0
excene=excene/au2eV
!Transition dipole moments are loaded as upper trigonal matrix, now we convert it as symmetry matrix
do i=0,nstates
    do j=i+1,nstates
        trandip(j,i,:)=trandip(i,j,:)
    end do
end do

!! Output some information
! write(*,*) "  State#      Exc.ene.(a.u.)     Transition dipole moment in X,Y,Z (a.u.)"
! do istat=1,nstates
!     write(*,"(i8,f20.6,3f15.6)") istat,excene(istat),trandip(0,istat,:)
! end do
do istat=1,nstates
    write(*,"(' State',i6,'   Excitation energy:',f12.5,' a.u.',f14.6,' eV')") istat,excene(istat),excene(istat)*au2eV
end do
write(*,"(' There are',i6,' excited states')") nstates
write(*,"(' Dipole moment of ground state:'3f12.5,' a.u.')") trandip(0,0,:)
write(*,*) "NOTE: All units used in this function is a.u."
!Gaussian output file is impossible to provide <m|r|n>, even if alltransitiondensities is used for CIS, it doesn't output <m|r|m>
do while(.true.)
write(*,*)
write(*,*) "0 Return"
write(*,*) "1 Calculate polarizability (alpha)"
if (ionlyalpha==0) write(*,*) "2 Calculate first hyperpolarizability (beta)"
if (ionlyalpha==0) write(*,*) "3 Calculate second hyperpolarizability (gamma)"
if (ionlyalpha==0) write(*,*) "4 Calculate third hyperpolarizability (delta)"
write(*,*) "5 Show the variation of alpha w.r.t. the number of states in consideration"
if (ionlyalpha==0) write(*,*) "6 Show the variation of beta w.r.t. the number of states in consideration"
if (ionlyalpha==0) write(*,*) "7 Show the variation of gamma w.r.t. the number of states in consideration"
write(*,*) "15 Calculate alpha in a range of frequencies"
if (ionlyalpha==0) write(*,*) "16 Calculate beta in a range of frequencies"
if (ionlyalpha==0) write(*,*) "17 Calculate gamma in a range of frequencies"
read(*,*) isel

if (isel==0) then
    return
else if (isel==1.or.isel==5.or.isel==15) then ! Calculate polarizability
    if (isel==1.or.isel==5) then
        write(*,*) "Input frequency of external field w for alpha(-w;w)"
        write(*,*) "e.g. 0.25"
        write(*,*) "Note: Negative value means using nm as unit, e.g. -693.5"
        read(*,*) freqbeg
        if (freqbeg<0) freqbeg=1240.7011D0/au2eV/abs(freqbeg)
        freqend=freqbeg
        freqstep=1
        if (freqbeg/=0) then
            wavlen=1240.7011D0/(freqbeg*au2eV)
            write(*,"(' Wavelength of w:',f12.3,' nm',/)") wavlen
        end if
    else if (isel==15) then
        write(*,*) "Input lower and upper limits as well as stepsize of w for alpha(-w;w)"
        write(*,*) "e.g. 0.2,0.5,0.01"
        read(*,*) freqbeg,freqend,freqstep
    end if
    
    if (isel==1) then
        istart=nstates
        iend=nstates
    else if (isel==5) then
        istart=1
        iend=nstates
        open(10,file="alpha_n.txt",status="replace")
    else if (isel==15) then
        istart=nstates
        iend=nstates
        open(10,file="alpha_w.txt",status="replace")
    end if
    
    write(*,*) "Please wait..."
    do numstat=istart,iend
    do freq=freqbeg,freqend,freqstep
        do idir=1,3
            do jdir=1,3
                tmpval=0
                do istat=1,numstat
                    tmpval1=trandip(0,istat,idir)*trandip(istat,0,jdir)/(excene(istat)-freq)
                    tmpval2=trandip(0,istat,jdir)*trandip(istat,0,idir)/(excene(istat)+freq)
                    tmpval=tmpval+tmpval1+tmpval2
                end do
                alpha(idir,jdir)=tmpval
            end do
        end do
        alphaiso=(alpha(1,1)+alpha(2,2)+alpha(3,3))/3D0
        term1=(alpha(1,1)-alpha(2,2))**2 + (alpha(1,1)-alpha(3,3))**2 + (alpha(2,2)-alpha(3,3))**2
        term2=6*(alpha(1,2)**2+alpha(1,3)**2+alpha(2,3)**2)
        alphaani1=dsqrt((term1+term2)/2D0)
        call diagmat(alpha,eigvecmat,eigval,300,1D-10)
        call sortr8(eigval)
        alphaani2=eigval(3)-(eigval(1)+eigval(2))/2D0
        
        if (isel==1) then
            write(*,*) "Polarizability tensor:"
            write(*,*) "             1              2              3"
            do idir=1,3
                write(*,"(i3,3f15.6)") idir,alpha(idir,:)
            end do
            write(*,"(' Isotropic average polarizability:',f15.6)") alphaiso
            write(*,"(' Polarizability anisotropy (definition 1):',f15.6)") alphaani1
            write(*,"(' Eigenvalues:',3f15.6)") eigval(:)
            write(*,"(' Polarizability anisotropy (definition 2):',f15.6)") alphaani2
        else if (isel==5) then
            write(10,"(i6,9f15.6)") numstat,alphaiso,alphaani1,alphaani2,alpha(1,1),alpha(2,2),alpha(3,3),alpha(1,2),alpha(1,3),alpha(2,3)
        else if (isel==15) then
            write(10,"(f12.6,9(1PE14.5))") freq,alphaiso,alphaani1,alphaani2,alpha(1,1),alpha(2,2),alpha(3,3),alpha(1,2),alpha(1,3),alpha(2,3)
        end if
    end do
    end do
    if (isel==5.or.isel==15) then
        close(10)
        if (isel==5) write(*,*) "Done! The result has been outputted to alpha_n.txt in current folder"
        if (isel==15) write(*,*) "Done! The result has been outputted to alpha_w.txt in current folder"
        write(*,*) "The correspondence between columns and information in this file is as follows"
        if (isel==5) write(*,*) "Column 1:  The number of states in consideration"
        if (isel==15)  write(*,*) "Column 1:  Frequency of external field"
        write(*,*) "Column 2:  Isotropic average polarizability"
        write(*,*) "Column 3:  Polarizability anisotropy (definition 1)"
        write(*,*) "Column 4:  Polarizability anisotropy (definition 2)"
        write(*,*) "Column 5:  XX"
        write(*,*) "Column 6:  YY"
        write(*,*) "Column 7:  ZZ"
        write(*,*) "Column 8:  XY"
        write(*,*) "Column 9:  XZ"
        write(*,*) "Column 10: YZ"
    end if

! Calculate first hyperpolarizability
else if (isel==2.or.isel==6.or.isel==16) then
    if (allocated(freqlist)) deallocate(freqlist)
    if (isel==2.or.isel==6) then
        nfreq=1
        allocate(freqlist(1,2))
        write(*,*) "Input frequency of external field w1 and w2 for beta(-(w1+w2);w1,w2)"
        write(*,*) "e.g. 0.25,0.13"
        write(*,*) "Note: Negative values mean using nm as unit, e.g. -693,0"
        read(*,*) freqlist(1,:)
        where(freqlist<0) freqlist=1240.7011D0/au2eV/abs(freqlist)
        if (freqlist(1,1)/=0) then
            wavlen1=1240.7011D0/(freqlist(1,1)*au2eV)
            write(*,"(' Wavelength of w1:',f12.3,' nm')") wavlen1
        end if
        if (freqlist(1,2)/=0) then
            wavlen2=1240.7011D0/(freqlist(1,2)*au2eV)
            write(*,"(' Wavelength of w2:',f12.3,' nm')") wavlen2
        end if
        write(*,*)
    else if (isel==16) then
        write(*,*) "Input the file recording frequency list, e.g. C:\freqlist.txt"
        do while(.true.)
            read(*,"(a)") c200tmp
            inquire(file=c200tmp,exist=alive)
            if (alive) exit
            write(*,*) "Cannot find the file, input again"
        end do
        open(10,file=c200tmp,status="old")
        nfreq=totlinenum(10,1)
        allocate(freqlist(nfreq,2))
        do ifreq=1,nfreq
            read(10,*) freqlist(ifreq,:)
        end do
        close(10)
        write(*,*) "The frequencies loaded:"
        do ifreq=1,nfreq
            write(*,"(' #',i5,'  w1=',f10.5,' a.u.',f10.3,' nm   w2=',f10.5' a.u.',f10.3,' nm')") &
            ifreq,freqlist(ifreq,1),1240.7011D0/(freqlist(ifreq,1)*au2eV),freqlist(ifreq,2),1240.7011D0/(freqlist(ifreq,2)*au2eV)
        end do
        write(*,*)
    end if

    if (isel==2) then
        istart=nstates
        iend=nstates
    else if (isel==6) then
        istart=1
        iend=nstates
        open(10,file="beta_n.txt",status="replace")
    else if (isel==16) then
        istart=nstates
        iend=nstates
        open(10,file="beta_w.txt",status="replace")
    end if
    
    write(*,*) "Please wait..."
    call fullarrange(arrb,6,3) !Generate full arrangement matrix (3!=6 :3) for beta, each row corresponds to one permutation, e.g. 231
    do numstat=istart,iend
    do ifreq=1,nfreq
        freq1=freqlist(ifreq,1)
        freq2=freqlist(ifreq,2)
        freqtot=freq1+freq2
        tmpw(1:3)=(/ -freqtot,freq1,freq2 /)
        do idir=1,3
            do jdir=1,3
                do kdir=1,3

                    tmpval=0
                    tmpdir(1:3)=(/ idir,jdir,kdir /)
                    do iper=1,6 !Do permutation, arrb(1,:)=1,2,3
                        dir1=tmpdir(arrb(iper,1))
                        dir2=tmpdir(arrb(iper,2))
                        dir3=tmpdir(arrb(iper,3))
                        w0=tmpw(arrb(iper,1))
                        w2=tmpw(arrb(iper,3))
                        do istat=1,numstat
                            do jstat=1,numstat
                                cen=trandip(istat,jstat,dir2)
                                if (istat==jstat) cen=cen-trandip(0,0,dir2)
                                tmpval=tmpval+trandip(0,istat,dir1)*cen*trandip(jstat,0,dir3)/(excene(istat)+w0)/(excene(jstat)-w2)
                            end do
                        end do
                    end do
                    beta(idir,jdir,kdir)=tmpval
                    
                    !Below is the code manually considering each permutation, for teaching purpose, but foolish
!                     tmpval=0
!                     do istat=1,numstat
!                         do jstat=1,numstat
!                             ! Consider six permutations, i=-freqtot, j=freq1, k=freq2,  the denominator is (+1th),(-3th)
!                             ! 1_3
!                             ! ijk; -freqtot,-freq2
!                             ! ikj; -freqtot,-freq1
!                             ! jik; +freq1  ,-freq2
!                             ! jki; +freq1  ,+freqtot
!                             ! kij; +freq2  ,-freq1
!                             ! kji; +freq2  ,+freqtot
!                             t1c=trandip(istat,jstat,jdir)
!                             t2c=trandip(istat,jstat,kdir)
!                             t3c=trandip(istat,jstat,idir)
!                             t4c=trandip(istat,jstat,kdir)
!                             t5c=trandip(istat,jstat,idir)
!                             t6c=trandip(istat,jstat,jdir)
!                             if (istat==jstat) then
!                                 t1c=t1c-trandip(0,0,jdir)
!                                 t2c=t2c-trandip(0,0,kdir)
!                                 t3c=t3c-trandip(0,0,idir)
!                                 t4c=t4c-trandip(0,0,kdir)
!                                 t5c=t5c-trandip(0,0,idir)
!                                 t6c=t6c-trandip(0,0,jdir)
!                             end if
!                             t1=trandip(0,istat,idir)*t1c*trandip(jstat,0,kdir) /(excene(istat)-freqtot)/(excene(jstat)-freq2)
!                             t2=trandip(0,istat,idir)*t2c*trandip(jstat,0,jdir) /(excene(istat)-freqtot)/(excene(jstat)-freq1)
!                             t3=trandip(0,istat,jdir)*t3c*trandip(jstat,0,kdir) /(excene(istat)+freq1)  /(excene(jstat)-freq2)
!                             t4=trandip(0,istat,jdir)*t4c*trandip(jstat,0,idir) /(excene(istat)+freq1)  /(excene(jstat)+freqtot)
!                             t5=trandip(0,istat,kdir)*t5c*trandip(jstat,0,jdir) /(excene(istat)+freq2)  /(excene(jstat)-freq1)
!                             t6=trandip(0,istat,kdir)*t6c*trandip(jstat,0,idir) /(excene(istat)+freq2)  /(excene(jstat)+freqtot)
!                             tmpval=tmpval+t1+t2+t3+t4+t5+t6
!                         end do
!                     end do
!                     beta(idir,jdir,kdir)=tmpval
                end do
            end do
        end do

        betaX=0
        betaY=0
        betaZ=0
        do j=1,3
            betaX=betaX+(beta(1,j,j)+beta(j,j,1)+beta(j,1,j))/3
            betaY=betaY+(beta(2,j,j)+beta(j,j,2)+beta(j,2,j))/3
            betaZ=betaZ+(beta(3,j,j)+beta(j,j,3)+beta(j,3,j))/3
        end do
        betatot=dsqrt(betaX**2+betaY**2+betaZ**2)
        dipx=trandip(0,0,1)
        dipy=trandip(0,0,2)
        dipz=trandip(0,0,3)
        dipnorm=dsqrt(dipx**2+dipy**2+dipz**2)
        betaprj=(betaX*dipx+betaY*dipy+betaZ*dipz)/dipnorm
        beta_per=0
        do j=1,3
            beta_per=beta_per+(2*beta(3,j,j)+2*beta(j,j,3)-3*beta(j,3,j))/5
        end do
        
        if (isel==2) then
            write(*,*) "First hyperpolarizability tensor:"
            do jdir=1,3
                do kdir=1,3
                    do idir=1,3
                        write(*,"(2x,3a,'=',f17.5,2x)",advance="no") dirlab(idir),dirlab(jdir),dirlab(kdir),beta(idir,jdir,kdir)
                        if (idir==3) write(*,*)
                    end do
                end do
            end do
            write(*,"(/,' Beta_X:',f17.5,'  Beta_Y:',f17.5,'  Beta_Z:',f17.5)") betaX,betaY,betaZ
            write(*,"(a,f17.5)") " Magnitude of beta:",betatot
            write(*,"(a,f17.5)") " Projection of beta on dipole moment:",betaprj
            write(*,"(a,f17.5)") " Beta ||     :",betaprj*3D0/5D0
            write(*,"(a,f17.5)") " Beta ||(z)  :",betaZ*3D0/5D0
            write(*,"(a,f17.5)") " Beta _|_(z) :",beta_per
        else if (isel==6) then
            write(10,"(i6,8(1PE14.5))") numstat,betaX,betaY,betaZ,betatot,betaprj,betaprj*3D0/5D0,betaZ*3D0/5D0,beta_per
        else if (isel==16) then
            write(10,"(2f12.6,8(1PE14.5))") freq1,freq2,betaX,betaY,betaZ,betatot,betaprj,betaprj*3D0/5D0,betaZ*3D0/5D0,beta_per
        end if
    end do
    end do
    
    if (isel==6.or.isel==16) then
        close(10)
        if (isel==6) then
            write(*,*) "Done! The result has been outputted to beta_n.txt in current folder"
            write(*,*) "The correspondence between columns and information in this file is as follows"
            write(*,*) "Column 1:  The number of states in consideration"
            write(*,*) "Column 2:  Beta_X"
            write(*,*) "Column 3:  Beta_Y"
            write(*,*) "Column 4:  Beta_Z"
            write(*,*) "Column 5:  Magnitude of hyperpolarizability"
            write(*,*) "Column 6:  Hyperpolarizability component along dipole moment direction"
            write(*,*) "Column 7:  Beta ||"
            write(*,*) "Column 8:  Beta ||(z)"
            write(*,*) "Column 9:  Beta _|_(z)"
        else if (isel==16) then
            write(*,*) "Done! The result has been outputted to beta_w.txt in current folder"
            write(*,*) "The correspondence between columns and information in this file is as follows"
            write(*,*) "Column 1:  Frequency of the first external field (w1)"
            write(*,*) "Column 2:  Frequency of the second external field (w2)"
            write(*,*) "Column 3:  Beta_X"
            write(*,*) "Column 4:  Beta_Y"
            write(*,*) "Column 5:  Beta_Z"
            write(*,*) "Column 6:  Magnitude of hyperpolarizability"
            write(*,*) "Column 7:  Hyperpolarizability component along dipole moment direction"
            write(*,*) "Column 8:  Beta ||"
            write(*,*) "Column 9:  Beta ||(z)"
            write(*,*) "Column 10: Beta _|_(z)"
        end if 
    end if
    
! Calculate second hyperpolarizability (gamma)
else if (isel==3.or.isel==7.or.isel==17) then
    if (allocated(freqlist)) deallocate(freqlist)
    if (isel==3.or.isel==7) then
        nfreq=1
        allocate(freqlist(1,3))
        write(*,*) "Input frequency of external fields w1, w2, w3 for gamma(-(w1+w2+w3);w1,w2,w3)"
        write(*,*) "e.g. 0.13,0.13,0"
        write(*,*) "Note: Negative values mean using nm as unit, e.g. -693,0,0"
        read(*,*) freqlist(1,:)
        where(freqlist<0) freqlist=1240.7011D0/au2eV/abs(freqlist)
        do i=1,3
            if (freqlist(1,i)/=0) then
                wavlen=1240.7011D0/(freqlist(1,i)*au2eV)
                write(*,"(' Wavelength of w',i1,':',f12.3,' nm')") i,wavlen
            end if
        end do
        write(*,*)
    else if (isel==17) then
        write(*,*) "Input the file recording frequency list, e.g. C:\freqlist.txt"
        do while(.true.)
            read(*,"(a)") c200tmp
            inquire(file=c200tmp,exist=alive)
            if (alive) exit
            write(*,*) "Cannot find the file, input again"
        end do
        open(10,file=c200tmp,status="old")
        nfreq=totlinenum(10,1)
        allocate(freqlist(nfreq,3))
        do ifreq=1,nfreq
            read(10,*) freqlist(ifreq,:)
        end do
        close(10)
        write(*,*) "The frequencies loaded:"
        do ifreq=1,nfreq
            write(*,"(' #',i5,'   w1=',f10.5,' a.u.    w2=',f10.5,' a.u.    w3=',f10.5,' a.u.')") ifreq,freqlist(ifreq,:)
        end do
        write(*,*)
    end if
    
    if (isel==3.or.isel==17) then
        write(*,"(' Consider how many states? Should be <=',i6)") nstates
        read(*,*) istart
        iend=istart
        nstatstep=1
    else if (isel==7) then
        write(*,*) "Input upper limit and stepsize of the number of states"
        write(*,*) "e.g. 150,2"
        read(*,*) iend,nstatstep
        istart=1
        if (iend>nstates) iend=nstates
    end if
    if (isel==7) open(10,file="gamma_n.txt",status="replace")
    if (isel==17) open(10,file="gamma_w.txt",status="replace")
    
    write(*,*) "Please wait..."
    call walltime(iwalltime1)
    call fullarrange(arrg,24,4) !Generate full arrangement matrix (4!=24 :4) for gamma, each row corresponds to one permutation, e.g. 2341
!     do i=1,24
!         write(*,"(i4,4i1)") arrg(i,:) 
!     end do
    iprog=0
    do numstat=istart,iend,nstatstep
    iprog=iprog+nstatstep
    if (isel==7) write(*,"(' Finished:',i5,' /',i5)") iprog,iend-istart+1
    do ifreq=1,nfreq
        freq1=freqlist(ifreq,1)
        freq2=freqlist(ifreq,2)
        freq3=freqlist(ifreq,3)
        freqtot=freq1+freq2+freq3
        tmpw(1:4)=(/ -freqtot,freq1,freq2,freq3 /)
    
        do idir=1,3 !Cycle direction component
            do jdir=1,3
                do kdir=1,3
                    do ldir=1,3
                    
                        gamma1=0
                        gamma2=0
                        tmpdir(1:4)=(/ idir,jdir,kdir,ldir /)
                        do iper=1,24 !Do permutation, arrg(1,:)=1,2,3,4
                            dir1=tmpdir(arrg(iper,1))
                            dir2=tmpdir(arrg(iper,2))
                            dir3=tmpdir(arrg(iper,3))
                            dir4=tmpdir(arrg(iper,4))
                            w0=tmpw(arrg(iper,1))
                            w1=tmpw(arrg(iper,2))
                            w2=tmpw(arrg(iper,3))
                            w3=tmpw(arrg(iper,4))
nthreads=getNThreads()
nthreads=getNThreads()
!$OMP PARALLEL SHARED(gamma1,gamma2) PRIVATE(istat,jstat,kstat,t1c,t2c,p1,p2,g1t,g2t) NUM_THREADS(nthreads)
                            g1t=0
                            g2t=0
!$OMP DO schedule(dynamic)
                            do istat=1,numstat
                                do jstat=1,numstat
                                    !Gamma 1
                                    do kstat=1,numstat
                                        t1c=trandip(istat,jstat,dir2)
                                        if (istat==jstat) t1c=t1c-trandip(0,0,dir2)
                                        t2c=trandip(jstat,kstat,dir3)
                                        if (jstat==kstat) t2c=t2c-trandip(0,0,dir3)
                                        p1=trandip(0,istat,dir1)*t1c*t2c*trandip(kstat,0,dir4)
                                        p2=(excene(istat)+w0)*(excene(jstat)-w2-w3)*(excene(kstat)-w3)
                                        g1t=g1t+p1/p2
                                    end do
                                    !Gamma 2
                                    p1=trandip(0,istat,dir1)*trandip(istat,0,dir2)*trandip(0,jstat,dir3)*trandip(jstat,0,dir4)
                                    p2=(excene(istat)+w0)*(excene(istat)-w1)*(excene(jstat)-w3) !(excene(jstat)-w3) can also be (excene(jstat)+w2), they are equivalent
                                    g2t=g2t+p1/p2
                                end do
                            end do
!$OMP END DO
!$OMP CRITICAL
                            gamma1=gamma1+g1t
                            gamma2=gamma2+g2t
!$OMP END CRITICAL
!$OMP END PARALLEL
                        end do
                        gamma(idir,jdir,kdir,ldir)=gamma1-gamma2
                        
                    end do
                end do
            end do
        end do
        
        gammaX=0
        gammaY=0
        gammaZ=0
        do i=1,3
            gammaX=gammaX+gamma(1,i,i,1)+gamma(1,i,1,i)+gamma(1,1,i,i)
            gammaY=gammaY+gamma(2,i,i,2)+gamma(2,i,2,i)+gamma(2,2,i,i)
            gammaZ=gammaZ+gamma(3,i,i,3)+gamma(3,i,3,i)+gamma(3,3,i,i)
        end do
        gammaX=gammaX/15
        gammaY=gammaY/15
        gammaZ=gammaZ/15
        gammatot=dsqrt(gammaX**2+gammaY**2+gammaZ**2)
        gammaavg1=gammaX+gammaY+gammaZ
        gammaavg2=( gamma(1,1,1,1)+gamma(2,2,2,2)+gamma(3,3,3,3) + gamma(1,1,2,2)+gamma(1,1,3,3)+gamma(2,2,3,3) + gamma(2,2,1,1)+gamma(3,3,1,1)+gamma(3,3,2,2) )/5
        
        if (isel==3) then
            write(*,*) "Second hyperpolarizability tensor:"
            do jdir=1,3
                do kdir=1,3
                    do ldir=1,3
                        do idir=1,3
                            write(*,"(2x,4a,'=',1PE14.5,2x)",advance="no") dirlab(idir),dirlab(jdir),dirlab(kdir),dirlab(ldir),gamma(idir,jdir,kdir,ldir)
                            if (idir==3) write(*,*)
                        end do
                    end do
                end do
            end do
            write(*,"(/,' Gamma_X:',1PE14.5,'  Gamma_Y:',1PE14.5,'  Gamma_Z:',1PE14.5)") gammaX,gammaY,gammaZ
            write(*,"(a,1PE14.5)") " Magnitude of gamma:",gammatot
            write(*,"(a,1PE14.5)") " Average of gamma (definition 1):",gammaavg1
            write(*,"(a,1PE14.5)") " Average of gamma (definition 2):",gammaavg2
            write(*,*)
        else if (isel==7) then
            write(10,"(i6,6(1PE14.5))") numstat,gammaX,gammaY,gammaZ,gammatot,gammaavg1,gammaavg2
        else if (isel==17) then
            write(10,"(3f12.6,6(1PE14.5))") freq1,freq2,freq3,gammaX,gammaY,gammaZ,gammatot,gammaavg1,gammaavg2
        end if
    end do !end cycle freqlist
    end do !end cycle the number of states
    
    call walltime(iwalltime2)
    write(*,"(' Calculation took up wall clock time',i10,'s')") iwalltime2-iwalltime1
    if (isel==7.or.isel==17) then
        close(10)
        if (isel==7) then
            write(*,*) "Done! The result has been outputted to gamma_n.txt in current folder"
            write(*,*) "The correspondence between columns and information in this file is as follows"
            write(*,*) "Column 1:  The number of states in consideration"
            write(*,*) "Column 2:  Gamma_X"
            write(*,*) "Column 3:  Gamma_Y"
            write(*,*) "Column 4:  Gamma_Z"
            write(*,*) "Column 5:  Magnitude of gamma"
            write(*,*) "Column 6:  Average of gamma (definition 1)"
            write(*,*) "Column 7:  Average of gamma (definition 2)"
        else if (isel==17) then
            write(*,*) "Done! The result has been outputted to gamma_w.txt in current folder"
            write(*,*) "The correspondence between columns and information in this file is as follows"
            write(*,*) "Column 1:  Frequency of the first external field (w1)"
            write(*,*) "Column 2:  Frequency of the second external field (w2)"
            write(*,*) "Column 3:  Frequency of the third external field (w3)"
            write(*,*) "Column 4:  Gamma_X"
            write(*,*) "Column 5:  Gamma_Y"
            write(*,*) "Column 6:  Gamma_Z"
            write(*,*) "Column 7:  Magnitude of gamma"
            write(*,*) "Column 8:  Average of gamma (definition 1)"
            write(*,*) "Column 9:  Average of gamma (definition 2)"
        end if 
    end if

! Calculate third hyperpolarizability
else if (isel==4) then
    write(*,*) "Input w1,w2,w3,w4 for delta(-(w1+w2+w3+w4);w1,w2,w3,w4)"
    write(*,*) "e.g. 0.13,0.13,0,-0.13"
    write(*,*) "Note: Negative values mean using nm as unit, e.g. -693,0,0"
    read(*,*) freq1,freq2,freq3,freq4
    if (freq1<0) freq1=1240.7011D0/au2eV/abs(freq1)
    if (freq2<0) freq2=1240.7011D0/au2eV/abs(freq2)
    if (freq3<0) freq3=1240.7011D0/au2eV/abs(freq3)
    if (freq4<0) freq4=1240.7011D0/au2eV/abs(freq4)
    freqtot=freq1+freq2+freq3+freq4
    tmpw(1:5)=(/ -freqtot,freq1,freq2,freq3,freq4 /)
    if (freq1/=0) then
        wavlen1=1240.7011D0/(freq1*au2eV)
        write(*,"(' Wavelength of w1:',f12.3,' nm')") wavlen1
    end if
    if (freq2/=0) then
        wavlen2=1240.7011D0/(freq2*au2eV)
        write(*,"(' Wavelength of w2:',f12.3,' nm')") wavlen2
    end if
    if (freq3/=0) then
        wavlen3=1240.7011D0/(freq3*au2eV)
        write(*,"(' Wavelength of w3:',f12.3,' nm')") wavlen3
    end if
    if (freq4/=0) then
        wavlen4=1240.7011D0/(freq4*au2eV)
        write(*,"(' Wavelength of w4:',f12.3,' nm')") wavlen4
    end if
    write(*,*)
    write(*,"(' Consider how many states? Should <=',i6)") nstates
    read(*,*) numstat
    
    write(*,*) "Please wait patiently..."
    call walltime(iwalltime1)
    call fullarrange(arrd,120,5) !Generate full arrangement matrix (4!=24 :4) for gamma, each row corresponds to one permutation, e.g. 2341
    iprog=0
    do idir=1,3 !Cycle direction component
        do jdir=1,3
            do kdir=1,3
                do ldir=1,3
                    iprog=iprog+1
                    write(*,"(' Finished:',i4,' /',i4)") iprog,81
                    do mdir=1,3
                
                        delta1=0
                        delta2=0
                        delta3=0
                        tmpdir(1:5)=(/ idir,jdir,kdir,ldir,mdir /)
                        do iper=1,120 !Do permutation, arrd(1,:)=1,2,3,4,5
                            dir1=tmpdir(arrd(iper,1))
                            dir2=tmpdir(arrd(iper,2))
                            dir3=tmpdir(arrd(iper,3))
                            dir4=tmpdir(arrd(iper,4))
                            dir5=tmpdir(arrd(iper,5))
                            w0=tmpw(arrd(iper,1))
                            w1=tmpw(arrd(iper,2))
                            w2=tmpw(arrd(iper,3))
                            w3=tmpw(arrd(iper,4))
                            w4=tmpw(arrd(iper,5))
nthreads=getNThreads()
nthreads=getNThreads()
!$OMP PARALLEL SHARED(delta1,delta2,delta3) PRIVATE(istat,jstat,kstat,lstat,t1c,t2c,t3c,p1,p2,p3,p4,d1t,d2t,d3t) NUM_THREADS(nthreads)
                            d1t=0
                            d2t=0
                            d3t=0
!$OMP DO schedule(dynamic)
                            do istat=1,numstat
                                do jstat=1,numstat
                                    do kstat=1,numstat
                                        !Delta 1
                                        do lstat=1,numstat
                                            t1c=trandip(istat,jstat,dir2)
                                            if (istat==jstat) t1c=t1c-trandip(0,0,dir2)
                                            t2c=trandip(jstat,kstat,dir3)
                                            if (jstat==kstat) t2c=t2c-trandip(0,0,dir3)
                                            t3c=trandip(kstat,lstat,dir4)
                                            if (kstat==lstat) t3c=t3c-trandip(0,0,dir4)
                                            p1=trandip(0,istat,dir1)*t1c*t2c*t3c*trandip(lstat,0,dir5)
                                            p2=(excene(istat)+w0)*(excene(jstat)+w0+w1)*(excene(kstat)-w3-w4)*(excene(lstat)-w4)
                                            d1t=d1t+p1/p2
                                        end do
                                        !Delta 2
                                        t3c=trandip(jstat,kstat,dir4)
                                        if (jstat==kstat) t3c=t3c-trandip(0,0,dir4)
                                        p1=trandip(0,istat,dir1)*trandip(istat,0,dir2)*trandip(0,jstat,dir3)*t3c*trandip(kstat,0,dir5)
                                        p2=(excene(jstat)+w2)*(excene(kstat)-w4)
                                        p3=1/(excene(istat)+w0)+1/(excene(istat)-w1)
                                        p4=1/(excene(jstat)-w3-w4)+1/(excene(kstat)+w2+w3)
                                        d2t=d2t+p1/p2*p3*p4
                                        !Delta 3, p1 is identical to delta 2 counterpart
                                        p2=(excene(istat)+w0)*(excene(istat)-w1)
                                        p3=1/(excene(jstat)-w3-w4)/(excene(kstat)-w4)
                                        p4=1/(excene(jstat)+w2)/(excene(kstat)+w2+w3)
                                        d3t=d3t+p1/p2*(p3+p4)
                                    end do
                                end do
                            end do
!$OMP END DO
!$OMP CRITICAL
                            delta1=delta1+d1t
                            delta2=delta2+d2t
                            delta3=delta3+d3t
!$OMP END CRITICAL
!$OMP END PARALLEL
                        end do
                        delta(idir,jdir,kdir,ldir,mdir)=delta1-delta2/2-delta3/2
                        
                    end do
                end do
            end do
        end do
    end do
    
    write(*,*)
    write(*,*) "Third hyperpolarizability tensor:"
    do jdir=1,3
        do kdir=1,3
            do ldir=1,3
                do mdir=1,3
                    do idir=1,3
                        write(*,"(2x,5a,'=',1PE14.5,2x)",advance="no") &
                        dirlab(idir),dirlab(jdir),dirlab(kdir),dirlab(ldir),dirlab(mdir),delta(idir,jdir,kdir,ldir,mdir)
                        if (idir==3) write(*,*)
                    end do
                end do
            end do
        end do
    end do
    call walltime(iwalltime2)
    write(*,"(' Calculation took up wall clock time',i10,'s')") iwalltime2-iwalltime1
    
end if
end do
end subroutine


!!---------- Calculate average bond length between two elements and average coordinate number
subroutine atmavgdist
use defvar
use util
implicit real*8 (a-h,o-z)
character elesel1*2,elesel2*2,selectyn
if (all(a%name==a(1)%name)) then !Only contain one kind of atom
    elesel1=a(1)%name
    elesel2=a(1)%name
else
    write(*,*) "Input two elements for which their average bond length will be calculated"
    write(*,*) "For example, B,Al"
    read(*,*) elesel1,elesel2
    call lc2uc(elesel1(1:1))
    call uc2lc(elesel1(2:2))
    call lc2uc(elesel2(1:1))
    call uc2lc(elesel2(2:2))
end if
write(*,*) "Input distance cutoff in Angstrom, e.g. 3.2"
read(*,*) discrit
discrit=discrit/b2a

call gendistmat
avgdist=0
iwithin=0
distmax=0
distmin=1000000
do iatm=1,ncenter
    do jatm=iatm+1,ncenter
        if ((a(iatm)%name==elesel1.and.a(jatm)%name==elesel2).or.(a(jatm)%name==elesel1.and.a(iatm)%name==elesel2)) then
            dist=distmat(iatm,jatm)
            if (dist<=discrit) then
                iwithin=iwithin+1
                write(*,"(i6,'#   ',i6,'(',a')   --',i6,'(',a,')    Length:',f12.6,' Angstrom')") iwithin,iatm,a(iatm)%name,jatm,a(jatm)%name,dist*b2a
                avgdist=avgdist+dist
                if ((iatm==1.and.jatm==iatm+1).or.dist>distmax) then
                    distmax=dist
                    idistmax1=iatm
                    idistmax2=jatm
                end if
                if ((iatm==1.and.jatm==iatm+1).or.dist<distmin) then
                    distmin=dist
                    idistmin1=iatm
                    idistmin2=jatm
                end if
            end if
        end if
    end do
end do
if (iwithin>0) then
    avgdist=avgdist/iwithin
    write(*,"(' Average bond length between ',a,1x,a,' is',f12.6,' Angstrom')") elesel1,elesel2,avgdist*b2a
    write(*,"(' Minimum length is',f12.6,' Angstrom, between',i6,'(',a,') and',i6,'(',a,')')") distmin*b2a,idistmin1,elesel1,idistmin2,elesel2
    write(*,"(' Maximum length is',f12.6,' Angstrom, between',i6,'(',a,') and',i6,'(',a,')')") distmax*b2a,idistmax1,elesel1,idistmax2,elesel2
else
    write(*,*) "No bond that satisfied your criterion is found"
    return
end if
write(*,*)
write(*,*) "If also calculate average coordinate number? (y/n)"
read(*,*) selectyn
if (selectyn=='y'.or.selectyn=='Y') then
    ncoordtot=0
    ntmp=0
    do iatm=1,ncenter
        if (a(iatm)%name/=elesel1) cycle
        ncoord=0
        do jatm=1,ncenter
            if (iatm==jatm.or.a(jatm)%name/=elesel2) cycle
            if (distmat(iatm,jatm)<=discrit) ncoord=ncoord+1
        end do
        write(*,"(' The coordinate number of',i6,'(',a,') due to ',a,' - ',a,' bond:',i5)") iatm,elesel1,elesel1,elesel2,ncoord
        ncoordtot=ncoordtot+ncoord
        ntmp=ntmp+1
    end do
    write(*,"(/,' The average coordinate number of ',a,' due to ',a,' - ',a,' bond:',f10.5)") elesel1,elesel1,elesel2,dfloat(ncoordtot)/ntmp
end if
end subroutine


        
!!!------- Calculate electric/magnetic/velocity... integral between orbitals
subroutine outorbint
use defvar
use util
implicit real*8 (a-h,o-z)
real*8,allocatable :: GTFint(:),GTFvecint(:,:)
real*8 vecint(3),vecinttmp(3),intval
write(*,*) "Output which kind of integral between orbitals?"
write(*,*) "1: Electric dipole moment  2: Magnetic dipole moment  3: Velocity"
write(*,*) "4: Kinetic energy   5: Overlap"
read(*,*) itype
if (wfntype==0.or.wfntype==1.or.wfntype==2) then
    write(*,*) "Output the integrals between which orbitals?"
    if (wfntype==0.or.wfntype==1.or.wfntype==2) then
        write(*,*) "1 Between all occupied orbitals"
        write(*,*) "2 Between all occupied and all unoccupied orbitals"
    end if
    write(*,*) "3 Between all orbitals"
    write(*,*) "4 Between the same orbitals"
    write(*,*) "5 Between specific range of orbitals"
    write(*,*) "6 Between two orbitals"
    read(*,*) irange
end if
if (irange==5) then
    write(*,*) "Input orbital range of the first index, e.g. 25,30"
    read(*,*) ibeg,iend
    write(*,*) "Input orbital range of the second index, e.g. 25,30"
    read(*,*) jbeg,jend
else if (irange==6) then
100    write(*,*) "Input index for the two orbitals, e.g. 144,340"
    write(*,*) "Input 0,0 can exit"
    read(*,*) iMOsel,jMOsel
    if (iMOsel==0.and.jMOsel==0) return
end if

call walltime(iwalltime1)

nsize=nprims*(nprims+1)/2
if (itype==1.or.itype==2.or.itype==3) then
    allocate(GTFvecint(3,nsize))
else
    allocate(GTFint(nsize))
end if
if (itype==1) then
    call genGTFDmat(GTFvecint,nsize)
else if (itype==2) then
    call genGTFMmat(GTFvecint,nsize) !Notice that this only generate lower triangular matrix, (i,j)=-(j,i)
else if (itype==3) then    
    call genGTFVelmat(GTFvecint,nsize)
else if (itype==4) then    
    call genGTFTmat(GTFint,nsize)
else if (itype==5) then    
    call genGTFSmat(GTFint,nsize)
end if

if (irange/=6) open(10,file="orbint.txt",status="replace")
do imo=1,nmo
    do jmo=1,nmo
        if (irange==1) then
            if (MOocc(imo)==0D0.or.MOocc(jmo)==0D0) cycle
        else if (irange==2) then
            if (MOocc(imo)==0D0.or.MOocc(jmo)/=0D0) cycle
        else if (irange==4) then
            if (imo/=jmo) cycle
        else if (irange==5) then
            if (imo<ibeg.or.imo>iend.or.jmo<jbeg.or.jmo>jend) cycle
        else if (irange==6) then
            if (imo/=iMOsel.or.jmo/=jMOsel) cycle
        end if
        
        if (itype==1.or.itype==2.or.itype==3) then !Vector integral
            vecint=0D0
nthreads=getNThreads()
nthreads=getNThreads()
!$OMP PARALLEL SHARED(vecint) PRIVATE(iGTF,jGTF,ides,vecinttmp) NUM_THREADS(nthreads)
            vecinttmp=0D0
!$OMP DO schedule(dynamic)
            do iGTF=1,nprims
                do jGTF=1,nprims
                    if (iGTF>=jGTF) then
                        ides=iGTF*(iGTF-1)/2+jGTF
                    else
                        ides=jGTF*(jGTF-1)/2+iGTF
                    end if
                    if ((itype==2.or.itype==3).and.iGTF>jGTF) then !Magnetic and velocity operators are Hermitean, so inverse sign
                        vecinttmp=vecinttmp-co(imo,iGTF)*co(jmo,jGTF)*GTFvecint(:,ides)
                    else
                        vecinttmp=vecinttmp+co(imo,iGTF)*co(jmo,jGTF)*GTFvecint(:,ides)
                    end if
                end do
            end do
!$OMP END DO
!$OMP CRITICAL
            vecint=vecint+vecinttmp
!$OMP END CRITICAL
!$OMP END PARALLEL
            if (irange==6) then
                write(*,"(' Result:',3f18.10)") vecint(:)
            else
                write(10,"(2i8,4x,3f18.10)") imo,jmo,vecint(:)
            end if
        else !Scalar integral
            intval=0D0
            do iGTF=1,nprims
                do jGTF=1,nprims
                    if (iGTF>=jGTF) then
                        ides=iGTF*(iGTF-1)/2+jGTF
                    else
                        ides=jGTF*(jGTF-1)/2+iGTF
                    end if
                    intval=intval+co(imo,iGTF)*co(jmo,jGTF)*GTFint(ides)
                end do
            end do
            if (irange==6) then
                write(*,"(' Result:',f18.10)") intval
            else
                write(10,"(2i8,4x,f18.10)") imo,jmo,intval
            end if
        end if
    end do
    if (irange/=6) write(*,"(' Finished',i7,' /',i7)") imo,nmo
end do
if (irange==6) then
    if (allocated(GTFvecint)) deallocate(GTFvecint)
    if (allocated(GTFint)) deallocate(GTFint)
    write(*,*)
    goto 100
else
    close(10)
    call walltime(iwalltime2)
    write(*,"(' Done! Calculation took up wall clock time',i10,'s',/)") iwalltime2-iwalltime1
    write(*,*) "The integrals has been outputted to orbint.txt in current folder"
    if (itype==1.or.itype==2.or.itype==3) then
        write(*,"(a)") " The first and the second columns correspond to orbital indices. The last three columns correspond to the integral in X/Y/Z (a.u.)"
    else
        write(*,"(a)") " The first and the second columns correspond to orbital indices. The last column corresponds to the integral value (a.u.)"
    end if
end if
end subroutine



!!----------- Calculate center, first and second moments of a real space function
subroutine funcmoment
use defvar
use util
use function
implicit real*8 (a-h,o-z)
real*8 intval,moment1(3),moment2(3,3),moment2nuc(3,3),funcval(radpot*sphpot),beckeweigrid(radpot*sphpot),eigvecmat(3,3),eigval(3)
type(content) gridatm(radpot*sphpot),gridatmorg(radpot*sphpot)
character selectyn
ifunc=1
cenx=0
ceny=0
cenz=0
do while(.true.)
    write(*,*) "0 Return"
    write(*,*) "1 Calculate the first and second moments of the function"
    write(*,*) "2 Calculate the center and integral of the function"
    write(*,"(a,i5)") " 3 Select the function to be studied, current:",ifunc
    write(*,"(a,3f11.5,' Ang')") " 4 Set the center for option 1, current:",cenx*b2a,ceny*b2a,cenz*b2a
    read(*,*) isel
    
    if (isel==0) then
        return
    else if (isel==3) then
        call selfunc_interface(ifunc)
        cycle
    else if (isel==4) then
        write(*,*) "Input X,Y,Z of the center in Angstrom, e.g. 2.0,0,1.5"
        read(*,*) cenx,ceny,cenz
        cenx=cenx/b2a !To Bohr
        ceny=ceny/b2a
        cenz=cenz/b2a
        cycle
    end if

    write(*,"(' Radial points:',i5,'    Angular points:',i5,'   Total:',i10,' per center')") radpot,sphpot,radpot*sphpot
    call gen1cintgrid(gridatmorg,iradcut)

    call walltime(iwalltime1)
    CALL CPU_TIME(time_begin)

    intval=0
    moment1=0
    moment2=0
    realcenx=0
    realceny=0
    realcenz=0
    do iatm=1,ncenter
        write(*,"(' Processing center',i6,'(',a2,')   /',i6)") iatm,a(iatm)%name,ncenter
        gridatm%x=gridatmorg%x+a(iatm)%x !Move quadrature point to actual position in molecule
        gridatm%y=gridatmorg%y+a(iatm)%y
        gridatm%z=gridatmorg%z+a(iatm)%z
nthreads=getNThreads()
!$OMP parallel do shared(funcval) private(i) num_threads(nthreads)
        do i=1+iradcut*sphpot,radpot*sphpot
            funcval(i)=calcfuncall(ifunc,gridatm(i)%x,gridatm(i)%y,gridatm(i)%z)
        end do
!$OMP end parallel do
        call gen1cbeckewei(iatm,iradcut,gridatm,beckeweigrid)
        do i=1+iradcut*sphpot,radpot*sphpot
            tmpval=funcval(i)*gridatmorg(i)%value*beckeweigrid(i)
            xtmp=gridatm(i)%x-cenx
            ytmp=gridatm(i)%y-ceny
            ztmp=gridatm(i)%z-cenz
            intval=intval+tmpval
            moment1(1)=moment1(1)+xtmp*tmpval
            moment1(2)=moment1(2)+ytmp*tmpval
            moment1(3)=moment1(3)+ztmp*tmpval
            moment2(1,1)=moment2(1,1)+xtmp*xtmp*tmpval
            moment2(2,2)=moment2(2,2)+ytmp*ytmp*tmpval
            moment2(3,3)=moment2(3,3)+ztmp*ztmp*tmpval
            moment2(1,2)=moment2(1,2)+xtmp*ytmp*tmpval
            moment2(2,3)=moment2(2,3)+ytmp*ztmp*tmpval
            moment2(1,3)=moment2(1,3)+xtmp*ztmp*tmpval
            realcenx=realcenx+gridatm(i)%x*tmpval
            realceny=realceny+gridatm(i)%y*tmpval
            realcenz=realcenz+gridatm(i)%z*tmpval
        end do
    end do
    CALL CPU_TIME(time_end)
    call walltime(iwalltime2)
    write(*,"(' Calculation took up CPU time',f12.2,'s, wall clock time',i10,'s',/)") time_end-time_begin,iwalltime2-iwalltime1
    
    if (isel==1) then
        moment2(3,1)=moment2(1,3)
        moment2(2,1)=moment2(1,2)
        moment2(3,2)=moment2(2,3)
        write(*,*) "Note: All units below are in a.u."
        write(*,"(/,' Integral:',1PE16.8,/)") intval
        write(*,"(' The first moment:')")
        write(*,"(' X= ',1PE16.8,'   Y= ',1PE16.8,'   Z= ',1PE16.8)") moment1
        write(*,"(' Norm= ',1PE16.8,/)") sum(moment1**2)
        write(*,"(' The second moment:')")
        write(*,"(' XX=',1PE16.8,'   XY=',1PE16.8,'   XZ=',1PE16.8)") moment2(1,:)
        write(*,"(' YX=',1PE16.8,'   YY=',1PE16.8,'   YZ=',1PE16.8)") moment2(2,:)
        write(*,"(' ZX=',1PE16.8,'   ZY=',1PE16.8,'   ZZ=',1PE16.8)") moment2(3,:)

        call diagmat(moment2,eigvecmat,eigval,300,1D-10)
        call sortr8(eigval)
        write(*,"(a,3(1PE16.8))") ' Eigenvalues:',eigval
        write(*,"(' Anisotropy:',1PE16.8,/)") eigval(3)-(eigval(1)+eigval(2))/2D0
        write(*,"(' Radius of gyration:',1PE16.8,/)") dsqrt((moment2(1,1)+moment2(2,2)+moment2(3,3))/intval)

        if (ifunc==1) then
            moment2nuc=0
            do iatm=1,ncenter
                xtmp=a(iatm)%x-cenx
                ytmp=a(iatm)%y-ceny
                ztmp=a(iatm)%z-cenz
                tmpval=a(iatm)%charge
                moment2nuc(1,1)=moment2nuc(1,1)+xtmp*xtmp*tmpval
                moment2nuc(2,2)=moment2nuc(2,2)+ytmp*ytmp*tmpval
                moment2nuc(3,3)=moment2nuc(3,3)+ztmp*ztmp*tmpval
                moment2nuc(1,2)=moment2nuc(1,2)+xtmp*ytmp*tmpval
                moment2nuc(2,3)=moment2nuc(2,3)+ytmp*ztmp*tmpval
                moment2nuc(1,3)=moment2nuc(1,3)+xtmp*ztmp*tmpval
            end do
            moment2nuc(3,1)=moment2nuc(1,3)
            moment2nuc(2,1)=moment2nuc(1,2)
            moment2nuc(3,2)=moment2nuc(2,3)
            write(*,*)
            write(*,"(' The quadrupole moment of nuclear charges:')")
            write(*,"(' XX=',f16.8,'   XY=',f16.8,'   XZ=',f16.8)") moment2nuc(1,:)
            write(*,"(' YX=',f16.8,'   YY=',f16.8,'   YZ=',f16.8)") moment2nuc(2,:)
            write(*,"(' ZX=',f16.8,'   ZY=',f16.8,'   ZZ=',f16.8)") moment2nuc(3,:)
            write(*,*)
            write(*,"(' The quadrupole moment of the system:')")
            write(*,"(' XX=',f16.8,'   XY=',f16.8,'   XZ=',f16.8)") moment2nuc(1,:)-moment2(1,:)
            write(*,"(' YX=',f16.8,'   YY=',f16.8,'   YZ=',f16.8)") moment2nuc(2,:)-moment2(2,:)
            write(*,"(' ZX=',f16.8,'   ZY=',f16.8,'   ZZ=',f16.8)") moment2nuc(3,:)-moment2(3,:)
        end if
        
    else if (isel==2) then
        write(*,"(' Integral:',1PE16.8,/)") intval
        realcenx=realcenx/intval
        realceny=realceny/intval
        realcenz=realcenz/intval
        write(*,"(' The center of the function:')")
        write(*,"(' X=',f16.8,' Y=',f16.8,' Z=',f16.8,' Angstrom',/)") realcenx*b2a,realceny*b2a,realcenz*b2a
        write(*,*) "Use this center for succeeding calculations? (y/n)"
        read(*,*) selectyn
        if (selectyn=='y') then
            cenx=realcenx
            ceny=realceny
            cenz=realcenz
        end if
    end if
end do
end subroutine



!!----------- Calculate energy index (EI) or bond polarity index (BPI)
!!J. Phys. Chem., 94, 5602-5607 and J. Phys. Chem.,96, 157-164
subroutine calcEIBPI
use defvar
implicit real*8 (a-h,o-z)
if (wfntype==3.or.wfntype==4) then
    write(*,*) "Error: Post-HF wavefunction is not supported yet!"
    return
end if
ninnerele=0
do i=1,ncenter
    if (int(a(i)%charge)/=a(i)%index) cycle !For the atom used ECP, skip it
    if (a(i)%index>2.and.a(i)%index<=10) ninnerele=ninnerele+2
    if (a(i)%index>10.and.a(i)%index<=18) ninnerele=ninnerele+10
    if (a(i)%index>18.and.a(i)%index<=36) ninnerele=ninnerele+18
    if (a(i)%index>36.and.a(i)%index<=54) ninnerele=ninnerele+36
    if (a(i)%index>54.and.a(i)%index<=86) ninnerele=ninnerele+54
    if (a(i)%index>86) ninnerele=ninnerele+86
end do
write(*,"(' The number of inner electrons is assumed to be',i5,/)") ninnerele
do while(.true.)
    write(*,*) "Calculate EI index for which atom? e.g. 5"
    write(*,*) "Input 0 can exit"
    read(*,*) iatm
    if (iatm==0) exit
    val1=0D0
    val2=0D0
    do imo=ninnerele/2+1,nbasis
        if (MOocc(imo)==0D0) exit
        compos=0
        do ibas=basstart(iatm),basend(iatm)
            compos=compos+CObasa(ibas,imo)**2
            do jbas=1,nbasis
                if (jbas==ibas) cycle
                compos=compos+CObasa(ibas,imo)*CObasa(jbas,imo)*Sbas(ibas,jbas)
            end do
        end do
!         write(*,"(i5,3f12.6)") imo,MOocc(imo),MOene(imo),compos
        val1=val1+MOocc(imo)*MOene(imo)*compos
        val2=val2+MOocc(imo)*compos
    end do
    if (wfntype==1) then !beta part
        do imo=nbasis+ninnerele/2+1,nmo
            if (MOocc(imo)==0D0) exit
            compos=0
            do ibas=basstart(iatm),basend(iatm)
                compos=compos+CObasb(ibas,imo-nbasis)**2
                do jbas=1,nbasis
                    if (jbas==ibas) cycle
                    compos=compos+CObasb(ibas,imo-nbasis)*CObasb(jbas,imo-nbasis)*Sbas(ibas,jbas)
                end do
            end do
!             write(*,"(i5,3f12.6)") imo,MOocc(imo),MOene(imo),compos
            val1=val1+MOocc(imo)*MOene(imo)*compos
            val2=val2+MOocc(imo)*compos
        end do
    end if
    write(*,"(' The numerator:  ',f12.6,' a.u.')") val1
    write(*,"(' The denominator:',f12.6,' a.u.')") val2 !Corresponding to Mulliken occupation number of this atom in valence MOs
    write(*,"(' The EI index:   ',f12.6,' a.u.',/)") val1/val2
end do
end subroutine



!!--------- Calculate electronic coupling matrix element by FCD method
!FCDtest.fch :1-16,17-31,32-47
subroutine FCD
use defvar
use util
implicit real*8 (a-h,o-z)
integer :: itranstype=1   !1=hole transfer, 2=electron transfer
integer :: iFCDorbtype=1   !For UHF, 1=alpha  2=beta
integer :: imethod=1  !1=FCD   2=GMH
integer,allocatable :: donor(:),bridge(:),acceptor(:),donorMO(:),bridgeMO(:),acceptorMO(:)
real*8 :: compthresFCD=0.89D0,miu1(3),miu2(3),miu12(3)
real*8 :: MOcomp(nbasis,3) !Composition of donor,bridge,acceptor in each MO
character c2000tmp*2000,c80tmp*80
if (wfntype/=0) then
    write(*,*) "Error: This module is only applicable to RHF/RKS wavefunction"
    return
end if
ndonoratm=0
nbridgeatm=0
nacceptoratm=0
ndonorMO=1
nbridgeMO=0
nacceptorMO=1

!For debug, GTG base stack example
ndonoratm=16
nbridgeatm=15
nacceptoratm=16
allocate(donor(ndonoratm))
allocate(bridge(nbridgeatm))
allocate(acceptor(nacceptoratm))
do i=1,16
    donor(i)=i
    if (i/=16) bridge(i)=i+16
    acceptor(i)=i+31
end do
ndonorMO=1
nbridgeMO=1
nacceptorMO=1
!For debug, GG base stack example
! ndonoratm=16
! nbridgeatm=0
! nacceptoratm=16
! allocate(donor(ndonoratm))
! allocate(bridge(nbridgeatm))
! allocate(acceptor(nacceptoratm))
! do i=1,16
!     donor(i)=i
!     acceptor(i)=i+16
! end do
! ndonorMO=1
! nbridgeMO=0
! nacceptorMO=1

write(*,*) "Citation of FCD method: J. Chem. Phys., 117, 5607 (2002)"
130 do while(.true.)
    write(*,*)
    write(*,*) "       ========== Fragment charge difference method (FCD) =========="
    if (itranstype==1) write(*,*) "-1 Show information of some highest occupied MOs"
    if (itranstype==2) write(*,*) "-1 Show information of some lowest unoccupied MOs"
    write(*,*) "0 Start calculation"
    write(*,"(a,i4)") " 1 Define donor fragment, number of atoms:",ndonoratm
    write(*,"(a,i4)") " 2 Define bridge fragment (if any), number of atoms:",nbridgeatm
    write(*,"(a,i4)") " 3 Define acceptor fragment, number of atoms:",nacceptoratm
    write(*,"(a,f5.1,a)") " 4 Set composition threshold for determining attribution, current:",compthresFCD*100,"%"
    write(*,"(a,3i2)") " 5 The number of MOs attributed to donor, bridge and acceptor, current:",ndonorMO,nbridgeMO,nacceptorMO
    if (itranstype==1) write(*,*) "6 Switch type of electron transfer, current: Hole transfer"
    if (itranstype==2) write(*,*) "6 Switch type of electron transfer, current: Electron transfer"
!     if (wfntype==1.and.iFCDorbtype==1) write(*,*) "7 Switch molecular orbital type, current: Alpha"
!     if (wfntype==1.and.iFCDorbtype==2) write(*,*) "7 Switch molecular orbital type, current: Beta"
    if (imethod==1) write(*,*) "8 Switch method between FCD and GMH, current: FCD"
    if (imethod==2) write(*,*) "8 Switch method between FCD and GMH, current: GMH"
    read(*,*) isel
    if (isel==-1) then
        if (itranstype==1) write(*,*) "Show how many highest occupied MOs? e.g. 30"
        if (itranstype==2) write(*,*) "Show how many lowest unoccupied MOs? e.g. 30"
        write(*,*) "If press ENTER directly, 50 MOs will be shown"
        read(*,"(a)") c80tmp
        if (c80tmp==" ") then
            nshowMO=50
        else
            read(c80tmp,*) nshowMO
        end if
        exit
    else if (isel==0) then
        if (ndonoratm==0.or.nacceptoratm==0) then
            write(*,*) "Error: You should define donor and acceptor first!"
            cycle
        end if
        exit
    else if (isel==1) then
        write(*,*) "Input the atom indices constituting the donor"
        write(*,*) "e.g. 1,3-5,9"
        read(*,"(a)") c2000tmp
        if (allocated(donor)) deallocate(donor)
        call str2arr(c2000tmp,ndonoratm)
        allocate(donor(ndonoratm))
        call str2arr(c2000tmp,ndonoratm,donor)
        write(*,*) "Done!"
    else if (isel==2) then
        write(*,*) "Input the atom indices constituting the bridge"
        write(*,*) "e.g. 1,3-5,9"
        write(*,*) "(If directly press ENTER then bridge will not be defined)"
        read(*,"(a)") c2000tmp
        if (c2000tmp==" ") then
            nbridgeatm=0
            cycle
        end if
        if (allocated(bridge)) deallocate(bridge)
        call str2arr(c2000tmp,nbridgeatm)
        allocate(bridge(nbridgeatm))
        call str2arr(c2000tmp,nbridgeatm,bridge)
        write(*,*) "Done!"
    else if (isel==3) then
        write(*,*) "Input the atom indices constituting the acceptor"
        write(*,*) "e.g. 1,3-5,9"
        read(*,"(a)") c2000tmp
        if (allocated(acceptor)) deallocate(acceptor)
        call str2arr(c2000tmp,nacceptoratm)
        allocate(acceptor(nacceptoratm))
        call str2arr(c2000tmp,nacceptoratm,acceptor)
        write(*,*) "Done!"
    else if (isel==4) then
        write(*,*) "Input composition threshold (%), e.g. 89"
        read(*,*) compthresFCD
        compthresFCD=compthresFCD/100
    else if (isel==5) then
        write(*,*) "Input the number of MOs attributed to donor, bridge and acceptor, respectively"
        write(*,*) "e.g. 2,2,1"
        read(*,*) ndonorMO,nbridgeMO,nacceptorMO
    else if (isel==6) then
        if (itranstype==1) then
            itranstype=2
        else if (itranstype==2) then
            itranstype=1
        end if
    else if (isel==7) then
        if (iFCDorbtype==1) then
            iFCDorbtype=2
        else if (iFCDorbtype==2) then
            iFCDorbtype=1
        end if
    else if (isel==8) then
        if (imethod==1) then
            if (igenDbas==0) then !Haven't calculated dipole moment integral matrix
                write(*,"(a)") " Warning: To use GMH method, you should set ""igenDbas"" in settings.ini to 1, &
                and then restart Multiwfn, so that dipole moment integral matrix could be generated when loading file"
                write(*,"(a)") " Alternatively, now you can press ENTER to reload input file and meantime generate dipole moment integral matrix"
                pause
                call dealloall
                write(*,*) "Reloading input file..."
                igenDbas=1
                call readinfile(firstfilename,1)
            end if
            imethod=2
        else if (imethod==2) then
            imethod=1
        end if
    end if
end do

iHOMO=nelec/2
iLUMO=iHOMO+1
write(*,"(' HOMO is',i6)") iHOMO
write(*,"(' LUMO is',i6)") iLUMO

!Determine attribution of complex MOs
if (isel==-1) then
    write(*,*) "    MO   Energy(eV)    Donor(%)   Bridge(%)  Acceptor(%)"
    icount=0
else
    allocate(donorMO(ndonorMO),bridgeMO(nbridgeMO),acceptorMO(nacceptorMO))
end if

idonor=1
ibridge=1
iacceptor=1
if (itranstype==1) then !hole transfer, use HOMO and below MOs
    iMOstart=iHOMO
    iMOend=1
    iMOstep=-1
else !electron transfer, use LUMO and above MOs
    iMOstart=iLUMO
    iMOend=nmo
    iMOstep=1
end if
do imo=iMOstart,iMOend,iMOstep
    compdonor=0
    do iatmidx=1,ndonoratm
        iatm=donor(iatmidx)
        do ibas=basstart(iatm),basend(iatm)
            do jbas=1,nbasis
                compdonor=compdonor+CObasa(ibas,imo)*CObasa(jbas,imo)*Sbas(ibas,jbas)
            end do
        end do
    end do
    compbridge=0
    do iatmidx=1,nbridgeatm
        iatm=bridge(iatmidx)
        do ibas=basstart(iatm),basend(iatm)
            do jbas=1,nbasis
                compbridge=compbridge+CObasa(ibas,imo)*CObasa(jbas,imo)*Sbas(ibas,jbas)
            end do
        end do
    end do
    compacceptor=0
    do iatmidx=1,nacceptoratm
        iatm=acceptor(iatmidx)
        do ibas=basstart(iatm),basend(iatm)
            do jbas=1,nbasis
                compacceptor=compacceptor+CObasa(ibas,imo)*CObasa(jbas,imo)*Sbas(ibas,jbas)
            end do
        end do
    end do
    MOcomp(imo,1)=compdonor
    MOcomp(imo,2)=compbridge
    MOcomp(imo,3)=compacceptor
    if (isel==-1) then !Examine composition of some MOs
        write(*,"(i8,f10.3,3f12.2)") imo,MOene(imo)*au2eV,MOcomp(imo,:)*100
        icount=icount+1
        if (icount==nshowMO) exit !Show up to "nshowMO" MOs
        cycle
    end if
    if (idonor<=ndonorMO.and.compdonor>compthresFCD) then
        donorMO(idonor)=imo
        idonor=idonor+1
        write(*,"(' MO',i4,' (',f9.3,'eV) is attributed to donor with composition of',f5.1'%')") imo,MOene(imo)*au2eV,compdonor*100
    else if (iacceptor<=nacceptorMO.and.compacceptor>compthresFCD) then
        acceptorMO(iacceptor)=imo
        iacceptor=iacceptor+1
        write(*,"(' MO',i4,' (',f9.3,'eV) is attributed to acceptor with composition of',f5.1'%')") imo,MOene(imo)*au2eV,compacceptor*100
    else if (ibridge<=nbridgeMO.and.compbridge>compthresFCD) then
        bridgeMO(ibridge)=imo
        ibridge=ibridge+1
        write(*,"(' MO',i4,' (',f9.3,'eV) is attributed to bridge with composition of',f5.1'%')") imo,MOene(imo)*au2eV,compbridge*100
    end if
    if (idonor>ndonorMO.and.ibridge>nbridgeMO.and.iacceptor>nacceptorMO) exit
end do

if (isel==-1) goto 130

write(*,*)
write(*,"(' Donor atoms:')")
write(*,"(15i5)") donor
if (nbridgeatm>0.and.nbridgeMO>0) then
    write(*,"(' Bridge atoms:')")
    write(*,"(15i5)") bridge
end if
write(*,"(' Acceptor atoms:')")
write(*,"(15i5)") acceptor

write(*,*)
write(*,*)
write(*,*) "             ======== Results between donor and acceptor ========"
write(*,*)
rmsH_DA=0
avgE_DA=0
do idx=1,ndonorMO
    imo=donorMO(idx)
    E1=MOene(imo)
    q1D=MOcomp(imo,1)
    q1A=MOcomp(imo,3)
    dq1=q1D-q1A
    do jdx=1,nacceptorMO
        jmo=acceptorMO(jdx)
        E2=MOene(jmo)
        enediff=abs(E2-E1) !Always take positive energy difference
        q2D=MOcomp(jmo,1)
        q2A=MOcomp(jmo,3)
        dq2=q2D-q2A
        q12_1=0
        write(*,"('                 **** MO',i4,' (donor) --- MO',i4,' (acceptor) ****')") imo,jmo
        write(*,"(' MO energy difference:',f12.4,' eV')") enediff*au2eV
        if (imethod==1) then
            do iatmidx=1,ndonoratm
                iatm=donor(iatmidx)
                do ibas=basstart(iatm),basend(iatm)
                    do jbas=1,nbasis
                        q12_1=q12_1+( CObasa(ibas,imo)*CObasa(jbas,jmo)+CObasa(ibas,jmo)*CObasa(jbas,imo) )*Sbas(ibas,jbas)
                    end do
                end do
            end do
            q12_1=q12_1/2
            q12_2=0
            do iatmidx=1,nacceptoratm
                iatm=acceptor(iatmidx)
                do ibas=basstart(iatm),basend(iatm)
                    do jbas=1,nbasis
                        q12_2=q12_2+( CObasa(ibas,imo)*CObasa(jbas,jmo)+CObasa(ibas,jmo)*CObasa(jbas,imo) )*Sbas(ibas,jbas)
                    end do
                end do
            end do
            q12_2=q12_2/2
            dq12=q12_1-q12_2
            H_DA=enediff*abs(dq12)/dsqrt((dq1-dq2)**2+4*dq12**2)
            write(*,"(' q1(donor):',f10.6,'   q1(acceptor):',f10.6,'   delta_q1:',f10.6)") q1D,q1A,dq1
            write(*,"(' q2(donor):',f10.6,'   q2(acceptor):',f10.6,'   delta_q2:',f10.6)") q2D,q2A,dq2
            write(*,"(' delta_q1 - delta_q2:',f10.6)") dq1-dq2
            write(*,"(' q12(donor):',f10.6,'  q12(acceptor):',f10.6,'  delta_q12:',f10.6)") q12_1,q12_2,dq12
        else if (imethod==2) then
            miu1=0
            miu2=0
            miu12=0
            do ibas=1,nbasis !Note that Dbas has already taken the negative sign of electron charge into account
                do jbas=1,nbasis
                    miu1(:)=miu1(:)+CObasa(ibas,imo)*CObasa(jbas,imo)*Dbas(:,ibas,jbas)
                    miu2(:)=miu2(:)+CObasa(ibas,jmo)*CObasa(jbas,jmo)*Dbas(:,ibas,jbas)
                    miu12(:)=miu12(:)+( CObasa(ibas,imo)*CObasa(jbas,jmo)+CObasa(ibas,jmo)*CObasa(jbas,imo) )*Dbas(:,ibas,jbas)
                end do
            end do
            miu12=miu12/2
            valmiu1=dsqrt(sum(miu1**2))
            valmiu2=dsqrt(sum(miu2**2))
            valmiu12=dsqrt(sum(miu12**2))
            write(*,"(' miu1 (X/Y/Z): ',3f11.6,'   Norm:',f11.6,' a.u.')") miu1,valmiu1
            write(*,"(' miu2 (X/Y/Z): ',3f11.6,'   Norm:',f11.6,' a.u.')") miu2,valmiu2
            write(*,"(' miu1-miu2 (X/Y/Z):',3f11.6,'   Norm:',f11.6,' a.u.')") miu1-miu2,dsqrt(sum((miu1-miu2)**2))
            write(*,"(' miu12 (X/Y/Z):    ',3f11.6,'   Norm:',f11.6,' a.u.')") miu12,valmiu12
            H_DA=enediff*valmiu12/dsqrt(sum((miu1-miu2)**2)+4*valmiu12**2)
        end if
        write(*,"(' H_DA:',f14.8,' eV')") H_DA*au2eV
        write(*,*)
        rmsH_DA=rmsH_DA+H_DA**2
        avgE_DA=avgE_DA+enediff
    end do
end do
rmsH_DA=dsqrt(rmsH_DA/ndonorMO/nacceptorMO)
avgE_DA=avgE_DA/(ndonorMO*nacceptorMO)
write(*,"(' Average MO energy difference E_DA:',f14.8,'eV')") avgE_DA*au2eV
write(*,"(' Root mean square H_DA:',f14.8,'eV')") rmsH_DA*au2eV

if (nbridgeatm>0.and.nbridgeMO>0) then
    write(*,*)
    write(*,*)
    write(*,*) "             ======== Results between donor and bridge ========"
    write(*,*)
    rmsH_DB=0
    avgE_DB=0
    do idx=1,ndonorMO
        imo=donorMO(idx)
        E1=MOene(imo)
        q1D=MOcomp(imo,1)
        q1B=MOcomp(imo,2)
        dq1=q1D-q1B
        do jdx=1,nbridgeMO
            jmo=bridgeMO(jdx)
            E2=MOene(jmo)
            enediff=abs(E2-E1)  !Always take positive energy difference
            q2D=MOcomp(jmo,1)
            q2B=MOcomp(jmo,2)
            dq2=q2D-q2B
            q12_1=0
            write(*,"('                 **** MO',i4,' (donor) --- MO',i4,' (bridge) ****')") imo,jmo
            write(*,"(' MO energy difference:',f12.4,' eV')") enediff*au2eV
            if (imethod==1) then
                do iatmidx=1,ndonoratm
                    iatm=donor(iatmidx)
                    do ibas=basstart(iatm),basend(iatm)
                        do jbas=1,nbasis
                            q12_1=q12_1+( CObasa(ibas,imo)*CObasa(jbas,jmo)+CObasa(ibas,jmo)*CObasa(jbas,imo) )*Sbas(ibas,jbas)
                        end do
                    end do
                end do
                q12_1=q12_1/2
                q12_2=0
                do iatmidx=1,nbridgeatm
                    iatm=bridge(iatmidx)
                    do ibas=basstart(iatm),basend(iatm)
                        do jbas=1,nbasis
                            q12_2=q12_2+( CObasa(ibas,imo)*CObasa(jbas,jmo)+CObasa(ibas,jmo)*CObasa(jbas,imo) )*Sbas(ibas,jbas)
                        end do
                    end do
                end do
                q12_2=q12_2/2
                dq12=q12_1-q12_2
                H_DB=enediff*abs(dq12)/dsqrt((dq1-dq2)**2+4*dq12**2)
                write(*,"(' q1(donor):',f10.6,'   q1(bridge):',f10.6,'   delta_q1:',f10.6)") q1D,q1B,dq1
                write(*,"(' q2(donor):',f10.6,'   q2(bridge):',f10.6,'   delta_q2:',f10.6)") q2D,q2B,dq2
                write(*,"(' delta_q1 - delta_q2:',f10.6)") dq1-dq2
                write(*,"(' q12(donor):',f10.6,'  q12(bridge):',f10.6,'  delta_q12:',f10.6)") q12_1,q12_2,dq12
            else if (imethod==2) then
                miu1=0
                miu2=0
                miu12=0
                do ibas=1,nbasis !Note that Dbas has already taken the negative sign of electron charge into account
                    do jbas=1,nbasis
                        miu1(:)=miu1(:)+CObasa(ibas,imo)*CObasa(jbas,imo)*Dbas(:,ibas,jbas)
                        miu2(:)=miu2(:)+CObasa(ibas,jmo)*CObasa(jbas,jmo)*Dbas(:,ibas,jbas)
                        miu12(:)=miu12(:)+( CObasa(ibas,imo)*CObasa(jbas,jmo)+CObasa(ibas,jmo)*CObasa(jbas,imo) )*Dbas(:,ibas,jbas)
                    end do
                end do
                miu12=miu12/2
                valmiu1=dsqrt(sum(miu1**2))
                valmiu2=dsqrt(sum(miu2**2))
                valmiu12=dsqrt(sum(miu12**2))
                write(*,"(' miu1 (X/Y/Z): ',3f11.6,'   Norm:',f11.6,' a.u.')") miu1,valmiu1
                write(*,"(' miu2 (X/Y/Z): ',3f11.6,'   Norm:',f11.6,' a.u.')") miu2,valmiu2
                write(*,"(' miu1-miu2 (X/Y/Z):',3f11.6,'   Norm:',f11.6,' a.u.')") miu1-miu2,dsqrt(sum((miu1-miu2)**2))
                write(*,"(' miu12 (X/Y/Z):    ',3f11.6,'   Norm:',f11.6,' a.u.')") miu12,valmiu12
                H_DB=enediff*valmiu12/dsqrt(sum((miu1-miu2)**2)+4*valmiu12**2)
            end if
            write(*,"(' H_DB:',f14.8,' eV')") H_DB*au2eV
            write(*,*)
            rmsH_DB=rmsH_DB+H_DB**2
            avgE_DB=avgE_DB+enediff
        end do
    end do
    rmsH_DB=dsqrt(rmsH_DB/ndonorMO/nbridgeMO)
    avgE_DB=avgE_DB/(ndonorMO*nbridgeMO)
    write(*,"(' Average MO energy difference E_DB:',f14.8,'eV')") avgE_DB*au2eV
    write(*,"(' Root mean square H_DB:',f14.8,'eV')") rmsH_DB*au2eV

    write(*,*)
    write(*,*)
    write(*,*) "             ======== Results between bridge and acceptor ========"
    write(*,*)
    rmsH_BA=0
    avgE_BA=0
    do idx=1,nbridgeMO
        imo=bridgeMO(idx)
        E1=MOene(imo)
        q1B=MOcomp(imo,2)
        q1A=MOcomp(imo,3)
        dq1=q1B-q1A
        do jdx=1,nacceptorMO
            jmo=acceptorMO(jdx)
            E2=MOene(jmo)
            enediff=abs(E2-E1)   !Always take positive energy difference
            q2B=MOcomp(jmo,2)
            q2A=MOcomp(jmo,3)
            dq2=q2B-q2A
            q12_1=0
            write(*,"('                 **** MO',i4,' (bridge) --- MO',i4,' (acceptor) ****')") imo,jmo
            write(*,"(' MO energy difference:',f12.4,' eV')") enediff*au2eV
            if (imethod==1) then
                do iatmidx=1,nbridgeatm
                    iatm=bridge(iatmidx)
                    do ibas=basstart(iatm),basend(iatm)
                        do jbas=1,nbasis
                            q12_1=q12_1+( CObasa(ibas,imo)*CObasa(jbas,jmo)+CObasa(ibas,jmo)*CObasa(jbas,imo) )*Sbas(ibas,jbas)
                        end do
                    end do
                end do
                q12_1=q12_1/2
                q12_2=0
                do iatmidx=1,nacceptoratm
                    iatm=acceptor(iatmidx)
                    do ibas=basstart(iatm),basend(iatm)
                        do jbas=1,nbasis
                            q12_2=q12_2+( CObasa(ibas,imo)*CObasa(jbas,jmo)+CObasa(ibas,jmo)*CObasa(jbas,imo) )*Sbas(ibas,jbas)
                        end do
                    end do
                end do
                q12_2=q12_2/2
                dq12=q12_1-q12_2
                H_BA=enediff*abs(dq12)/dsqrt((dq1-dq2)**2+4*dq12**2)
                write(*,"(' q1(bridge):',f10.6,'   q1(acceptor):',f10.6,'   delta_q1:',f10.6)") q1D,q1A,dq1
                write(*,"(' q2(bridge):',f10.6,'   q2(acceptor):',f10.6,'   delta_q2:',f10.6)") q2D,q2A,dq2
                write(*,"(' delta_q1 - delta_q2:',f10.6)") dq1-dq2
                write(*,"(' q12(bridge):',f10.6,'  q12(acceptor):',f10.6,'  delta_q12:',f10.6)") q12_1,q12_2,dq12
            else if (imethod==2) then
                miu1=0
                miu2=0
                miu12=0
                do ibas=1,nbasis !Note that Dbas has already taken the negative sign of electron charge into account
                    do jbas=1,nbasis
                        miu1(:)=miu1(:)+CObasa(ibas,imo)*CObasa(jbas,imo)*Dbas(:,ibas,jbas)
                        miu2(:)=miu2(:)+CObasa(ibas,jmo)*CObasa(jbas,jmo)*Dbas(:,ibas,jbas)
                        miu12(:)=miu12(:)+( CObasa(ibas,imo)*CObasa(jbas,jmo)+CObasa(ibas,jmo)*CObasa(jbas,imo) )*Dbas(:,ibas,jbas)
                    end do
                end do
                miu12=miu12/2
                valmiu1=dsqrt(sum(miu1**2))
                valmiu2=dsqrt(sum(miu2**2))
                valmiu12=dsqrt(sum(miu12**2))
                write(*,"(' miu1 (X/Y/Z): ',3f11.6,'   Norm:',f11.6,' a.u.')") miu1,valmiu1
                write(*,"(' miu2 (X/Y/Z): ',3f11.6,'   Norm:',f11.6,' a.u.')") miu2,valmiu2
                write(*,"(' miu1-miu2 (X/Y/Z):',3f11.6,'   Norm:',f11.6,' a.u.')") miu1-miu2,dsqrt(sum((miu1-miu2)**2))
                write(*,"(' miu12 (X/Y/Z):    ',3f11.6,'   Norm:',f11.6,' a.u.')") miu12,valmiu12
                H_BA=enediff*valmiu12/dsqrt(sum((miu1-miu2)**2)+4*valmiu12**2)
            end if
            write(*,"(' H_BA:',f14.8,' eV')") H_BA*au2eV
            write(*,*)
            rmsH_BA=rmsH_BA+H_BA**2
            avgE_BA=avgE_BA+enediff
        end do
    end do
    rmsH_BA=dsqrt(rmsH_BA/nbridgeMO/nacceptorMO)
    avgE_BA=avgE_BA/(nbridgeMO*nacceptorMO)
    write(*,"(' Average MO energy difference E_BA:',f14.8,'eV')") avgE_BA*au2eV
    write(*,"(' Root mean square H_BA:',f14.8,'eV')") rmsH_BA*au2eV
    write(*,*)
!     write(*,*) "                       ======== Final result ========"
!     write(*,*)
!     H_SE=rmsH_DB*rmsH_BA/(avgE_DA-avgE_DB/2)
!     write(*,"(' H_SE (Superexchange):   ',f14.8,'eV')") H_SE
!     H_DBA=rmsH_DA+H_SE
!     write(*,"(' H_DBA (H_DA via bridge):',f14.8,'eV')") H_DBA
end if

deallocate(donorMO,bridgeMO,acceptorMO)

goto 130
end subroutine



!---------- Perform Pipek-Mezey orbital localization
!The algorithm can be found in J. Comput. Chem., 14, 736 (1993)
!Performing E2 analysis based on LMO oribtal is not appropriate, the non-diagonal elements of Fock matrix are almost zero, &
!this phenonmenon can also be seen in NLMO Fock matrix ($NBO FNLMO $END). This is because we don't allow mixure between occupied and unoccupied MOs
!
!The final wavefunction can be exported as .fch, I don't select .molden because it doesn't record atomic charge, this will be problematic when ECP is used
subroutine pipek_mezey
use defvar
use util
implicit real*8 (a-h,o-z)
integer :: maxcyc=80,ireload=1,idoene=0,domark(4)
real*8 :: arrayi(nbasis),arrayj(nbasis),crit=1D-4,tmpbasarr(nbasis),tmpprimarr(nprims)
real*8,pointer :: Cmat(:,:)
real*8,allocatable :: FmatA(:,:),FmatB(:,:),FLMOA(:,:),FLMOB(:,:)
character c200tmp*200,typestr*4
if (wfntype==2.or.wfntype==3.or.wfntype==4) then
write(*,*) "Error: This function only works for restricted or unrestricted SCF wavefunction!"
    write(*,*) "Press ENTER to return"
    pause
    return
end if

do while(.true.)
    if (idoene==1) write(*,*) "-4 If calculate and print orbital energies, current: Yes"
    if (idoene==0) write(*,*) "-4 If calculate and print orbital energies, current: No"
    if (ireload==1) write(*,*) "-3 If finally reload newly generated .fch file, current: Yes"
    if (ireload==0) write(*,*) "-3 If finally reload newly generated .fch file, current: No"
    write(*,"(a,f12.8)") " -2 Set criterion of convergence, current:",crit
    write(*,"(a,i4)") " -1 Set maximum cycles, current:",maxcyc
    write(*,*) "0 Return"
    write(*,*) "1 Localizing occupied orbitals only"
    write(*,*) "2 Localizing both occupied and unoccupied orbitals separately"
    read(*,*) isel
    if (isel==0) then
        return
    else if (isel==-1) then
        write(*,*) "Input maximum cycles, e.g. 30"
        read(*,*) maxcyc
    else if (isel==-2) then
        write(*,*) "Input criterion of convergence, e.g. 0.0001"
        read(*,*) crit
    else if (isel==-3) then
        if (ireload==1) then
            ireload=0
        else if (ireload==0) then
            ireload=1
        end if
    else if (isel==-4) then
        if (idoene==1) then
            idoene=0
        else if (idoene==0) then
            write(*,"(a)") " Input the file recording Fock matrix in original basis functions in lower triangular form, e.g. C:\fock.txt"
            read(*,"(a)") c200tmp
            inquire(file=c200tmp,exist=alive)
            if (alive.eqv. .false.) then
                write(*,*) "Error: Unable to find this file!"
                cycle
            end if
            open(10,file=c200tmp,status="old")
            if (.not.allocated(FmatA)) allocate(FmatA(nbasis,nbasis))
            read(10,*) ((FmatA(i,j),j=1,i),i=1,nbasis)
            do i=1,nbasis !Fill upper triangular part
                 do j=i+1,nbasis
                     FmatA(i,j)=FmatA(j,i)
                 end do
            end do
            if (wfntype==1) then
                if (.not.allocated(FmatB)) allocate(FmatB(nbasis,nbasis))
                read(10,*) ((FmatB(i,j),j=1,i),i=1,nbasis)
                do i=1,nbasis
                     do j=i+1,nbasis
                         FmatB(i,j)=FmatB(j,i)
                     end do
                end do
            end if
            close(10)
            write(*,*) "Fock matrix loaded successfully!"
            write(*,*)
            idoene=1
        end if
    else if (isel==1.or.isel==2) then
        exit
    end if
end do

call walltime(iwalltime1)
CALL CPU_TIME(time_begin)

domark=0
if (wfntype==0) then
    domark(1)=1
    if (isel==2) domark(2)=1
else if (wfntype==1) then
    domark(1)=1
    domark(3)=1
    if (isel==2) domark=1
end if

!Alpha-occ,Alpha-vir,Beta-occ,Beta-vir
do itime=1,4
    if (domark(itime)==0) cycle
    if (itime<=2) then
        Cmat=>CObasa
    else
        Cmat=>CObasb
    end if
    if (itime==1) then
        nmobeg=1
        nmoend=naelec
    else if (itime==2) then
        nmobeg=naelec+1
        nmoend=nbasis
    else if (itime==3) then
        nmobeg=1
        nmoend=nbelec
    else if (itime==4) then
        nmobeg=nbelec+1
        nmoend=nbasis
    end if
    
    if (wfntype==0) then
        if (itime==1) write(*,"(/,a)") " Localizing occupied orbitals..."
        if (itime==2) write(*,"(/,a)") " Localizing unoccupied orbitals..."
    else if (wfntype==1) then
        if (itime==1) write(*,"(/,a)") " Localizing alpha occupied orbitals..."
        if (itime==2) write(*,"(/,a)") " Localizing alpha unoccupied orbitals..."
        if (itime==3) write(*,"(/,a)") " Localizing beta occupied orbitals..."
        if (itime==4) write(*,"(/,a)") " Localizing beta unoccupied orbitals..."
    end if

    Pvalold=0
    do icyc=1,maxcyc
        do imo=nmobeg,nmoend-1
            do jmo=imo+1,nmoend
                Aval=0
                Bval=0
nthreads=getNThreads()
!$OMP parallel shared(Aval,Bval) private(Avalpriv,Bvalpriv,Qij,Qii,Qjj) num_threads(nthreads)
                Avalpriv=0
                Bvalpriv=0
!$OMP do schedule(DYNAMIC)
                do iatm=1,ncenter
                    Qij=Qval(Cmat,imo,jmo,iatm)
                    Qii=Qval(Cmat,imo,imo,iatm)
                    Qjj=Qval(Cmat,jmo,jmo,iatm)
                    Avalpriv=Avalpriv+( Qij**2-(Qii-Qjj)**2/4D0 )
                    Bvalpriv=Bvalpriv+( Qij*(Qii-Qjj) )
                end do
!$OMP END DO
!$OMP CRITICAL
                Aval=Aval+Avalpriv
                Bval=Bval+Bvalpriv
!$OMP END CRITICAL
!$OMP END PARALLEL
                if (Aval==0.and.Bval==0) cycle
                gamma=sign(1D0,Bval)*acos(-Aval/dsqrt(Aval**2+Bval**2))/4D0
                arrayi=cos(gamma)*Cmat(:,imo)+sin(gamma)*Cmat(:,jmo)
                arrayj=-sin(gamma)*Cmat(:,imo)+cos(gamma)*Cmat(:,jmo)
                Cmat(:,imo)=arrayi
                Cmat(:,jmo)=arrayj
    !             write(*,*) imo,jmo,gamma,Aval,Bval
    !             pause
            end do
        end do
        Pval=0
        do imo=nmobeg,nmoend
            do iatm=1,ncenter
                Pval=Pval+Qval(Cmat,imo,imo,iatm)**2
            end do
        end do
        deltaPval=Pval-Pvalold
        write(*,"(' Cycle:',i5,'  P:',f16.8,'  Delta P:',f16.8)") icyc,Pval,deltaPval
        if (abs(deltaPval)<crit) exit
        Pvalold=Pval
    end do
    if (icyc==maxcyc+1) then
        write(*,*) "Warning: Convergence failed!"
    else
        write(*,"(a)") " Successfully converged!"
    end if
end do

CALL CPU_TIME(time_end)
call walltime(iwalltime2)
write(*,"(/,' Calculation took up CPU time',f12.2,'s, wall clock time',i10,'s')") time_end-time_begin,iwalltime2-iwalltime1

!Print orbital energies, sort orbital according to energies and do E2 analysis
if (idoene==1) then
    !Do Alpha part or closed-shell orbitals
    allocate(FLMOA(nbasis,nbasis))
    FLMOA=matmul(matmul(transpose(CObasa),FmatA),CObasa)
    if (isel==1) nmoend=naelec
    if (isel==2) nmoend=nbasis
    do iorb=1,nmoend
        MOene(iorb)=FLMOA(iorb,iorb)
    end do
    do iorb=1,nmoend
        do jorb=iorb+1,nmoend
            if (MOene(iorb)>MOene(jorb)) then    
                tmpbasarr=CObasa(:,iorb)
                CObasa(:,iorb)=CObasa(:,jorb)
                CObasa(:,jorb)=tmpbasarr
                tmpprimarr=CO(iorb,:)
                CO(iorb,:)=CO(jorb,:)
                CO(jorb,:)=tmpprimarr
                tmpene=MOene(iorb)
                MOene(iorb)=MOene(jorb)
                MOene(jorb)=tmpene
            end if
        end do
    end do
    !Do beta part
    if (wfntype==1) then
        allocate(FLMOB(nbasis,nbasis))
        FLMOB=matmul(matmul(transpose(CObasb),FmatB),CObasb)
        if (isel==1) nmoend=nbelec
        if (isel==2) nmoend=nbasis
        do iorb=1,nmoend
            MOene(nbasis+iorb)=FLMOB(iorb,iorb)
        end do
        do iorb=1,nmoend
            do jorb=iorb+1,nmoend
                if (MOene(nbasis+iorb)>MOene(nbasis+jorb)) then    
                    tmpbasarr=CObasb(:,iorb)
                    CObasb(:,iorb)=CObasb(:,jorb)
                    CObasb(:,jorb)=tmpbasarr
                    tmpprimarr=CO(nbasis+iorb,:)
                    CO(nbasis+iorb,:)=CO(nbasis+jorb,:)
                    CO(nbasis+jorb,:)=tmpprimarr
                    tmpene=MOene(nbasis+iorb)
                    MOene(nbasis+iorb)=MOene(nbasis+jorb)
                    MOene(nbasis+jorb)=tmpene
                end if
            end do
        end do
    end if
    !Print energies, those not localized are not printed
    write(*,*) "Energies of localized orbitals:"
    do iorb=1,nmo
        if (isel==1) then
            if (wfntype==0) then
                if (iorb>nint(naelec)) cycle
            else if (wfntype==1) then
                if (MOtype(iorb)==1.and.iorb>nint(naelec)) cycle
                if (MOtype(iorb)==2.and.iorb>nbasis+nint(nbelec)) cycle
            end if
        end if
        if (MOtype(iorb)==0) typestr="A+B"
        if (MOtype(iorb)==1) typestr="A"
        if (MOtype(iorb)==2) typestr="B"
        write(*,"(i6,'   Energy (a.u.):',f16.8,'    Type: ',a,'  Occ:',f4.1)") iorb,MOene(iorb),typestr,MOocc(iorb)
    end do
    if (isel==1) write(*,*) "Energies of unoccupied orbitals are not updated since they were not localized"
    
    !Second-order perturbation analysis between occupied and virtual orbitals, this only works when both of them have been localized
    !This part is commented since it don't print any useful E(2), because Fock element between occupied and virtual orbitals are always nearly zero
!     if (isel==1) then
!         write(*,*) " Note: E(2) analysis is skipped since virtual orbitals were not localized"
!     else if (isel==2) then
!         write(*,*) "Second-order perturbation theory analysis of interaction energy:"
!         !Regenerated Fock matrix in LMO, since they have been sorted
!         FLMOA=matmul(matmul(transpose(CObasa),FmatA),CObasa)
!         call showmatgau(FLMOA)
!         pause
!         coeff=2
!         do iocc=1,naelec
!             do ivir=naelec+1,nbasis
!                 E2val=-coeff* FLMOA(iocc,ivir)**2/(MOene(ivir)-MOene(iocc))
!                 qCT=coeff* ( FLMOA(iocc,ivir)/(MOene(ivir)-MOene(iocc)) )**2
!                 if (abs(E2val*au2kcal)>0.2D0) then
!                     write(*,"(' Donor:',i5,'  -  Acceptor:',i5,'   E(2):',f7.2,' kcal/mol   q_CT:',f10.5)") iocc,ivir,E2val*au2kcal,qCT
!                     write(*,"(3f16.10)") FLMOA(iocc,ivir),MOene(ivir)-MOene(iocc)
!                 end if
!             end do
!         end do
!     end if
end if

write(*,*)
call outfch("new.fch",10)
write(*,*) "The localized orbitals have been exported as new.fch in current folder"
if (ireload==1) then
    call dealloall
    write(*,*) "Loading new.fch..."
    call readfch("new.fch",1)
    write(*,"(a)") " Loading finished, now you can use main function 0 to visualize them as isosurface"
end if
end subroutine

!------ Calculate Q value used in Pipek-Mezey localization
real*8 function Qval(Cmat,imo,jmo,iatm)
use defvar
implicit real*8 (a-h,o-z)
real*8 Cmat(nbasis,nbasis)
Qval=0
do ibas=basstart(iatm),basend(iatm)
    Qval=Qval+sum( (Cmat(:,imo)*Cmat(ibas,jmo)+Cmat(ibas,imo)*Cmat(:,jmo)) *Sbas(ibas,:) )
end do
Qval=Qval/2D0
end function





!!------------ Integrate any real space function within isosurface of a real space function
!I use the same data structure as basin analysis to illustrate definition of isosurfaces
subroutine intisosurface
use defvar
use util
use function
use basinintmod !Use its vec26 array
implicit real*8 (a-h,o-z)
integer :: ifunciso=13,ifuncint=1
integer,allocatable :: mergelist(:),grididx(:,:,:),dogrid(:,:)
character :: defdomain*20="<0.5",c1000tmp*1000

if (allocated(gridxyz)) deallocate(gridxyz,domainsize,domaingrid)
do while(.true.)
    write(*,*)
    if (allocated(cubmat)) write(*,*) "-1 Directly use the grid data in memory to yield domains"
    write(*,*) "0 Return"
    write(*,*) "1 Calculate grid data and yield domains"
    write(*,"(a,i5)") " 2 Select real space function to be calculated, current:",ifunciso
    write(*,"(a,a)") " 3 Set criterion for defining domain, current: ",trim(defdomain)
    read(*,*) isel
    if (isel==0) then
        return
    else if (isel==1.or.isel==-1) then
        exit
    else if (isel==2) then
        call selfunc_interface(ifunciso)
    else if (isel==3) then
        write(*,*) "Input the definition, e.g. <0.05"
        write(*,*) "Note: The first character must be < or >"
        read(*,"(a)") defdomain
    end if
end do

!Set grid and generate grid data
if (isel==1) then
    call setgridfixspc
    dvol=dx*dy*dz
    if (allocated(cubmat)) deallocate(cubmat)
    allocate(cubmat(nx,ny,nz))
    call savecubmat(ifunciso,0,iorbsel)
end if

!Count the number of grids satisfying the criterion
read(defdomain(2:),*) valiso
if (defdomain(1:1)=='<') then
    ngrid=count(cubmat<valiso)
else if (defdomain(1:1)=='>') then
    ngrid=count(cubmat>valiso)
else
    write(*,*) "Error: Parsing of the domain definition failed!"
    return
end if
write(*,"(/,' The number of grids satisfied the criterion:',i10)") ngrid

!Clustering grids that satisfied criterion to domain
!The idea is very clever: I assign each grid with a different index, and frequently perform iteration, in each cycle all grids (that meets isovalue criterion) &&
!are run over, and index of each grid is set to that of one of the nearest 26 grids as long as its index is larger than current grid. Finally, &
!grids in each domian will have identical indices
write(*,*) "Clustering domains..."
call walltime(iwalltime1)
CALL CPU_TIME(time_begin)
!dogrid: XYZ index of useful grids
!gridxyz: XYZ coordinate of useful grids
!grididx: currently records initial index of all gris
allocate(dogrid(ngrid,3),grididx(nx,ny,nz),gridxyz(ngrid,3))
grididx=-1 !Irrelevant grids have very small value
idx=0
do iz=1,nz
    do iy=1,ny
        do ix=1,nx
            if ((defdomain(1:1)=='<'.and.cubmat(ix,iy,iz)<valiso).or.(defdomain(1:1)=='>'.and.cubmat(ix,iy,iz)>valiso)) then
                idx=idx+1
                dogrid(idx,1)=ix
                dogrid(idx,2)=iy
                dogrid(idx,3)=iz
                grididx(ix,iy,iz)=idx
                gridxyz(idx,1)=orgx+(ix-1)*dx
                gridxyz(idx,2)=orgy+(iy-1)*dy
                gridxyz(idx,3)=orgz+(iz-1)*dz
            end if
        end do
    end do
end do

call setupmovevec
icyc=0
do while(.true.)
    iupdate=0
    icyc=icyc+1
!     write(*,*) icyc
    do itmp=1,ngrid
        ix=dogrid(itmp,1)
        iy=dogrid(itmp,2)
        iz=dogrid(itmp,3)
        do imove=1,26
            ixtmp=ix+vec26x(imove)
            iytmp=iy+vec26y(imove)
            iztmp=iz+vec26z(imove)
            if (ixtmp<1.or.ixtmp>nx.or.iytmp<1.or.iytmp>ny.or.iztmp<1.or.iztmp>nz) then
                write(*,"(a)") " Error: Some domains may be truncated! You should enlarge extension size or completely re-defined grid to avoid this circumstance!"
                write(*,*) "Press Enter to exit"
                pause
                return
            end if
            if (grididx(ix,iy,iz)<grididx(ixtmp,iytmp,iztmp)) then
                grididx(ix,iy,iz)=grididx(ixtmp,iytmp,iztmp)
                iupdate=1
                exit
            end if
        end do
    end do
    if (iupdate==0) exit
end do
!     do itmp=1,ngrid
!         ix=dogrid(itmp,1)
!         iy=dogrid(itmp,2)
!         iz=dogrid(itmp,3)
!         write(*,*) itmp,grididx(ix,iy,iz)
!     end do
!     pause
ndone=0
ndomain=0
do while(.true.)
    ndomain=ndomain+1
    mintmp=1000000
    do itmp=1,ngrid
        idxtmp=grididx(dogrid(itmp,1),dogrid(itmp,2),dogrid(itmp,3))
        if (idxtmp<ndomain) cycle
        if (idxtmp<mintmp) mintmp=idxtmp
    end do
    do itmp=1,ngrid
        idxtmp=grididx(dogrid(itmp,1),dogrid(itmp,2),dogrid(itmp,3))
        if (idxtmp==mintmp) then
            grididx(dogrid(itmp,1),dogrid(itmp,2),dogrid(itmp,3))=ndomain
            ndone=ndone+1
        end if
    end do
    if (ndone==ngrid) exit
!     write(*,*) ndomain,ndone
end do
!Generate domainsize and domaingrid (grid index that contained in each domain)
allocate(domainsize(ndomain),domaingrid(ndomain,ngrid))
do idom=1,ndomain
    j=0
    do itmp=1,ngrid
        if (grididx(dogrid(itmp,1),dogrid(itmp,2),dogrid(itmp,3))==idom) then
            j=j+1
            domaingrid(idom,j)=itmp
        end if
    end do
    domainsize(idom)=j
end do

write(*,*) "Clustering domains finished!"
CALL CPU_TIME(time_end)
call walltime(iwalltime2)
write(*,"(' Clustering took up CPU time',f12.2,'s, wall clock time',i10,'s')") time_end-time_begin,iwalltime2-iwalltime1

do idom=1,ndomain
    write(*,"(' Domain:',i6,'     Grids:',i8)") idom,domainsize(idom)
end do

do while(.true.)
    write(*,*)
    write(*,*) "-1 Merge specific domains"
    write(*,*) "0 Exit"
    write(*,*) "1 Perform integration for a domain"
    write(*,*) "2 Perform integration for all domains"
    write(*,*) "3 Visualize domains"
    write(*,"(a,i5)") " 4 Select the real space function to be integrated, current:",ifuncint
    write(*,*) "5 Calculate q_bind index for a domain"
    read(*,*) isel2
    if (isel2==0) then
        return
    else if (isel2==-1) then
        if (ndomain<2) then
            write(*,*) "Error: At least two domains must be presented!"
            cycle
        end if
        write(*,*) "Input indices of the domains you want to merge, e.g. 4,5,8-10"
        read(*,"(a)") c1000tmp
        call str2arr(c1000tmp,nmerge) !Find how many terms
        allocate(mergelist(nmerge))
        call str2arr(c1000tmp,nmerge,mergelist)
        call sorti4(mergelist)
        idom=mergelist(1)
        do jdx=nmerge,2,-1
            jdom=mergelist(jdx)
            nsizei=domainsize(idom)
            nsizej=domainsize(jdom)
            domainsize(idom)=nsizei+nsizej
            domaingrid(idom,nsizei+1:nsizei+nsizej)=domaingrid(jdom,1:nsizej)
            ndomain=ndomain-1
            domainsize(jdom:ndomain)=domainsize(jdom+1:ndomain+1)
            domaingrid(jdom:ndomain,:)=domaingrid(jdom+1:ndomain+1,:)
        end do
        deallocate(mergelist)
        write(*,"(a,i6)") " Done! The domains you selected have been merged as domain",idom
        write(*,*) "Size of current domains:"
        do idom=1,ndomain
            write(*,"(' Domain:',i6,'     Grids:',i8)") idom,domainsize(idom)
        end do
    else if (isel2==1) then
        write(*,*) "Input the index of the domain to be integrated, e.g. 3"
        read(*,*) intdom
        if (intdom<1.or.intdom>ndomain) then
            write(*,"(a)") " Error: The index of the domain to be integrated is incorrect"
            cycle
        end if
        valint=0
        volint=0
        valmin=1D100
        valmax=-1D100
        do igrd=1,domainsize(intdom)
            idx=domaingrid(intdom,igrd)
            tmpval=calcfuncall(ifuncint,gridxyz(idx,1),gridxyz(idx,2),gridxyz(idx,3))
            valint=valint+tmpval
            volint=volint+1
            if (tmpval<valmin) valmin=tmpval
            if (tmpval>valmax) valmax=tmpval
        end do
        avgval=valint/domainsize(intdom)
        valint=valint*dvol
        volint=volint*dvol
        write(*,"(' Integration result:',E20.10,' a.u.')") valint
        write(*,"(' Volume:',f12.6,' Bohr^3  (',f12.6,' Angstrom^3 )')") volint,volint*b2a**3
        write(*,"(' Average:',E20.10)") avgval
        write(*,"(' Maximum:',E20.10,'   Minimum:',E20.10)") valmax,valmin
    else if (isel2==2) then
        write(*,*) "Domain    Integral (a.u.)     Volume (Bohr^3)      Average"
        valinttot=0
        volinttot=0
        do intdom=1,ndomain
            valint=0
            volint=0
            do igrd=1,domainsize(intdom)
                idx=domaingrid(intdom,igrd)
                valint=valint+calcfuncall(ifuncint,gridxyz(idx,1),gridxyz(idx,2),gridxyz(idx,3))
                volint=volint+1
            end do
            avgval=valint/domainsize(intdom)
            valint=valint*dvol
            volint=volint*dvol
            write(*,"(i6,E20.10,f17.6,E20.10)") intdom,valint,volint,avgval
            valinttot=valinttot+valint
            volinttot=volinttot+volint
        end do
        write(*,"(' Integration result of all domains:',E20.10,' a.u.')") valinttot
        write(*,"(' Volume of all domains:',f13.6,' Bohr^3  ',f13.6,' Angstrom^3')") volinttot,volinttot*b2a**3
    else if (isel2==3) then
        idrawdomain=1
        aug3Dold=aug3D
        if (aug3D<3) aug3D=3 !Often we set extension distance to zero, e.g. RDG, in this case the molecule will be truncated, therefore here temporarily augment it
        aug3D=aug3Dold
        idrawdomain=0
    else if (isel2==4) then
        call selfunc_interface(ifuncint)
    else if (isel2==5) then
        write(*,*) "Input the index of the domain to be integrated, e.g. 3"
        read(*,*) intdom
        if (intdom<1.or.intdom>ndomain) then
            write(*,"(a)") " Error: The index of the domain to be integrated is incorrect"
            cycle
        end if
        qatt=0
        qrep=0
        expfac=4D0/3D0
        volneg=0
        volpos=0
        do igrd=1,domainsize(intdom)
            idx=domaingrid(intdom,igrd)
            call signlambda2rho_RDG(gridxyz(idx,1),gridxyz(idx,2),gridxyz(idx,3),sl2r,RDG,rho)
            if (sl2r<0) then
                qatt=qatt+rho**expfac
                volneg=volneg+1
            else
                qrep=qrep+rho**expfac
                volpos=volpos+1
            end if
        end do
        qatt=qatt*dvol
        qrep=qrep*dvol
        qbind=-(qatt-qrep)
        volneg=volneg*dvol
        volpos=volpos*dvol
        write(*,"(' q_att: 'f16.8,' a.u.')") qatt
        write(*,"(' q_rep: 'f16.8,' a.u.')") qrep
        write(*,"(' q_bind:'f16.8,' a.u.')") qbind
        write(*,"(' Volume (lambda2<0):',f13.6,' Bohr^3  ',f13.6,' Angstrom^3')") volneg
        write(*,"(' Volume (lambda2>0):',f13.6,' Bohr^3  ',f13.6,' Angstrom^3')") volpos
        write(*,"(' Volume (Total):    ',f13.6,' Bohr^3  ',f13.6,' Angstrom^3')") volneg+volpos
    end if
end do
end subroutine
