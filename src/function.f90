!! ----------------- Show the list of all supported real space functions
subroutine funclist
use defvar
! write(*,*) "            ----------- Avaliable real space functions -----------"
if (allocated(b)) then
    write(*,*) "1 Electron density"
    write(*,*) "2 Gradient norm of electron density"
    write(*,*) "3 Laplacian of electron density"
    write(*,*) "4 Value of orbital wavefunction"
    if (ipolarpara==0) write(*,*) "5 Electron spin density"
    if (ipolarpara==1) write(*,*) "5 Spin polarization parameter function"
    write(*,*) "6 Hamiltonian kinetic energy density K(r)"
    write(*,*) "7 Lagrangian kinetic energy density G(r)"
    if (ifiletype==4) then
        write(*,*) "8 Electrostatic potential from atomic charges"
    else
        write(*,*) "8 Electrostatic potential from nuclear charges"
    end if
    if (ELFLOL_type==0) write(*,*) "9 Electron Localization Function (ELF)"
    if (ELFLOL_type==1) write(*,*) "9 Electron Localization Function (ELF) defined by Tsirelson" 
    if (ELFLOL_type==2) write(*,*) "9 Electron Localization Function (ELF) defined by Lu, Tian" 
    if (ELFLOL_type==0) write(*,*) "10 Localized orbital locator (LOL)"
    if (ELFLOL_type==1) write(*,*) "10 Localized orbital locator (LOL) defined by Tsirelson" 
    if (ELFLOL_type==2) write(*,*) "10 Localized orbital locator (LOL) defined by Lu, Tian" 
    write(*,*) "11 Local information entropy"
    write(*,*) "12 Total electrostatic potential (ESP)"
    write(*,*) "13 Reduced density gradient (RDG)"
    write(*,*) "14 Reduced density gradient (RDG) with promolecular approximation"
    write(*,*) "15 Sign(lambda2)*rho"
    write(*,*) "16 Sign(lambda2)*rho with promolecular approximation"
    !Fermi hole function only available to single-determinant wavefunction
    if (pairfunctype==1) write(*,"(a,3f10.5)") " 17 Correlation hole for alpha, ref. point:",refx,refy,refz
    if (pairfunctype==2) write(*,"(a,3f10.5)") " 17 Correlation hole for beta, ref. point:",refx,refy,refz
    if (pairfunctype==4) write(*,"(a,3f10.5)") " 17 Correlation factor for alpha, ref. point:",refx,refy,refz
    if (pairfunctype==5) write(*,"(a,3f10.5)") " 17 Correlation factor for beta, ref. point:",refx,refy,refz
    if (pairfunctype==7) write(*,"(a,3f10.5)") " 17 Exc.-corr. density for alpha, ref. point:",refx,refy,refz
    if (pairfunctype==8) write(*,"(a,3f10.5)") " 17 Exc.-corr. density for beta, ref. point:",refx,refy,refz
    if (pairfunctype==10) write(*,"(a,3f10.5)") " 17 Pair density for alpha, ref. point:",refx,refy,refz
    if (pairfunctype==11) write(*,"(a,3f10.5)") " 17 Pair density for beta, ref. point:",refx,refy,refz
    if (pairfunctype==12) write(*,"(a,3f10.5)") " 17 Pair density for all electrons, ref. point:",refx,refy,refz
    write(*,*) "18 Average local ionization energy"
    write(*,"(a,i2,a,3f10.5)") " 19 Source function, mode:",srcfuncmode,", ref. point:",refx,refy,refz
    write(*,"(a,i5)") " 100 User defined real space function, iuserfunc=",iuserfunc
else
    if (ifiletype==4) then
        write(*,*) "8 Electrostatic potential from atomic charges"
    else
        write(*,*) "8 Electrostatic potential from nuclear charges"
    end if
    write(*,*) "14 Reduced density gradient(RDG) with promolecular approximation"
    write(*,*) "16 Sign(lambda2)*rho with promolecular approximation"
    write(*,"(a,i3)") " 100 User defined real space function, iuserfunc=",iuserfunc
end if
end subroutine

!---- Standard interface for selecting real space function
!Note that iorbsel is a global variable
subroutine selfunc_interface(ifunc)
use defvar
integer ifunc
call funclist
read(*,*) ifunc
if (ifunc==4) then
    write(*,"(a,i10)") " Input orbital index, between 1 and",nmo
    read(*,*) iorbsel
end if
end subroutine






! ======================================
! =========== Module function ==========
! ======================================

module function
use defvar
implicit real*8 (a-h,o-z)

contains


!-------- Calculate any supported real space function at a given point
real*8 function calcfuncall(ifunc,x,y,z)
integer ifunc
real*8 x,y,z
if (ifunc==1) then
calcfuncall=fdens(x,y,z)
else if (ifunc==2) then
calcfuncall=fgrad(x,y,z,'t')
else if (ifunc==3) then
calcfuncall=flapl(x,y,z,'t')
else if (ifunc==4) then
calcfuncall=fmo(x,y,z,iorbsel)
else if (ifunc==5) then
calcfuncall=fspindens(x,y,z,'s')
else if (ifunc==6) then
calcfuncall=Hamkin(x,y,z,0)
else if (ifunc==7) then
calcfuncall=lagkin(x,y,z,0)
else if (ifunc==8) then
calcfuncall=nucesp(x,y,z)
else if (ifunc==9) then
calcfuncall=ELF_LOL(x,y,z,"ELF")
else if (ifunc==10) then
calcfuncall=ELF_LOL(x,y,z,"LOL")
else if (ifunc==11) then
calcfuncall=infoentro(1,x,y,z)
else if (ifunc==12) then
calcfuncall=totesp(x,y,z)
else if (ifunc==13) then
calcfuncall=fgrad(x,y,z,'r')
else if (ifunc==14) then
calcfuncall=RDGprodens(x,y,z)
else if (ifunc==15) then
calcfuncall=signlambda2rho(x,y,z)
else if (ifunc==16) then
calcfuncall=signlambda2rho_prodens(x,y,z)
else if (ifunc==17) then
calcfuncall=pairfunc(refx,refy,refz,x,y,z)
else if (ifunc==18) then
calcfuncall=avglocion(x,y,z)
else if (ifunc==19) then
calcfuncall=srcfunc(x,y,z,srcfuncmode)
else if (ifunc==100) then
calcfuncall=userfunc(x,y,z)
end if
end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate wavefunction value of a range of orbitals and their derivatives at a given point, up to third-order
!! istart and iend is the range of the orbitals will be calculated, to calculate all orbitals, use 1,nmo
!! runtype=1: value  =2: value+dx/y/z  =3: value+dxx/yy/zz(diagonal of hess)  =4: value+dx/y/z+Hessian  
!!        =5: value+dx/y/z+hess+3-order derivative tensor 
subroutine orbderv(runtype,istart,iend,x,y,z,wfnval,grad,hess,tens3)
real*8 x,y,z,wfnval(nmo)
real*8,optional :: grad(3,nmo),hess(3,3,nmo),tens3(3,3,3,nmo)
integer runtype,istart,iend

wfnval=0D0
if (present(grad)) grad=0D0
if (present(hess)) hess=0D0
if (present(tens3)) tens3=0D0
lastcen=-1 !arbitrary value

! if the center/exp of current GTF is the same as previous, then needn't recalculate them
do j=1,nprims
    ix=type2ix(b(j)%functype)
    iy=type2iy(b(j)%functype)
    iz=type2iz(b(j)%functype)
    ep=b(j)%exp
    
    if (b(j)%center/=lastcen) then
        sftx=x-a(b(j)%center)%x
        sfty=y-a(b(j)%center)%y
        sftz=z-a(b(j)%center)%z
        sftx2=sftx*sftx
        sfty2=sfty*sfty
        sftz2=sftz*sftz
        rr=sftx2+sfty2+sftz2
    end if
    if (expcutoff>0.or.-ep*rr>expcutoff) then
        expterm=exp(-ep*rr)
    else
        expterm=0D0
    end if
    lastcen=b(j)%center
!     expterm=exp(-ep*dsqrt(rr))
    if (expterm==0D0) cycle
    
    !Calculate value for current GTF
    if (b(j)%functype==1) then !Some functype use manually optimized formula for cutting down computational time
    GTFval=expterm
    else if (b(j)%functype==2) then
    GTFval=sftx*expterm
    else if (b(j)%functype==3) then
    GTFval=sfty*expterm
    else if (b(j)%functype==4) then
    GTFval=sftz*expterm
    else if (b(j)%functype==5) then
    GTFval=sftx2*expterm
    else if (b(j)%functype==6) then
    GTFval=sfty2*expterm
    else if (b(j)%functype==7) then
    GTFval=sftz2*expterm
    else if (b(j)%functype==8) then
    GTFval=sftx*sfty*expterm
    else if (b(j)%functype==9) then
    GTFval=sftx*sftz*expterm
    else if (b(j)%functype==10) then
    GTFval=sfty*sftz*expterm
    else !If above condition is not satisfied(Angular moment higher than f), the function will calculated explicitly
    GTFval=sftx**ix *sfty**iy *sftz**iz *expterm
    end if
    !Calculate orbital wavefunction value
    do imo=istart,iend
        wfnval(imo)=wfnval(imo)+co(imo,j)*GTFval
    end do
    
    if (runtype>=2) then
        !Calculate 1-order derivative for current GTF
        tx=0.0D0
        ty=0.0D0
        tz=0.0D0
        if (ix/=0) tx=ix*sftx**(ix-1)
        if (iy/=0) ty=iy*sfty**(iy-1)
        if (iz/=0) tz=iz*sftz**(iz-1)
        GTFdx=sfty**iy *sftz**iz *expterm*(tx-2*ep*sftx**(ix+1))
        GTFdy=sftx**ix *sftz**iz *expterm*(ty-2*ep*sfty**(iy+1))
        GTFdz=sftx**ix *sfty**iy *expterm*(tz-2*ep*sftz**(iz+1))
        !Calculate 1-order derivative for orbitals
        do imo=istart,iend
            grad(1,imo)=grad(1,imo)+co(imo,j)*GTFdx
            grad(2,imo)=grad(2,imo)+co(imo,j)*GTFdy
            grad(3,imo)=grad(3,imo)+co(imo,j)*GTFdz
        end do

        if (runtype>=3) then
            !Calculate 2-order derivative for current GTF
            txx=0.0D0
            tyy=0.0D0
            tzz=0.0D0
            if (ix>=2) txx=ix*(ix-1)*sftx**(ix-2)
            if (iy>=2) tyy=iy*(iy-1)*sfty**(iy-2)
            if (iz>=2) tzz=iz*(iz-1)*sftz**(iz-2)
            GTFdxx=sfty**iy *sftz**iz *expterm*( txx + 2*ep*sftx**ix*(-2*ix+2*ep*sftx2-1) )
            GTFdyy=sftx**ix *sftz**iz *expterm*( tyy + 2*ep*sfty**iy*(-2*iy+2*ep*sfty2-1) )
            GTFdzz=sftx**ix *sfty**iy *expterm*( tzz + 2*ep*sftz**iz*(-2*iz+2*ep*sftz2-1) )
            ttx=tx-2*ep*sftx**(ix+1)
            tty=ty-2*ep*sfty**(iy+1)
            ttz=tz-2*ep*sftz**(iz+1)
            GTFdxy=sftz**iz *expterm*ttx*tty
            GTFdyz=sftx**ix *expterm*tty*ttz
            GTFdxz=sfty**iy *expterm*ttx*ttz
            !Calculate diagonal Hessian elements for orbitals
            do imo=istart,iend
                hess(1,1,imo)=hess(1,1,imo)+co(imo,j)*GTFdxx !dxx
                hess(2,2,imo)=hess(2,2,imo)+co(imo,j)*GTFdyy !dyy
                hess(3,3,imo)=hess(3,3,imo)+co(imo,j)*GTFdzz !dzz
            end do
            if (runtype>=4) then !Also process nondiagonal elements
                do imo=istart,iend
                    hess(1,2,imo)=hess(1,2,imo)+co(imo,j)*GTFdxy !dxy
                    hess(2,3,imo)=hess(2,3,imo)+co(imo,j)*GTFdyz !dyz
                    hess(1,3,imo)=hess(1,3,imo)+co(imo,j)*GTFdxz !dxz
                end do
                hess(2,1,:)=hess(1,2,:)
                hess(3,2,:)=hess(2,3,:)
                hess(3,1,:)=hess(1,3,:)
            end if
            
            if (runtype>=5) then
                !Calculate 3-order derivative for current GTF
                ep2=ep*2D0
                ep4=ep*4D0
                epep4=ep2*ep2
                epep8=epep4*2D0
                !dxyz
                a1=0D0
                b1=0D0
                c1=0D0
                if (ix>=1) a1=ix*sftx**(ix-1)
                if (iy>=1) b1=iy*sfty**(iy-1)
                if (iz>=1) c1=iz*sftz**(iz-1)
                a2=-ep2*sftx**(ix+1)
                b2=-ep2*sfty**(iy+1)
                c2=-ep2*sftz**(iz+1)
                GTFdxyz=(a1+a2)*(b1+b2)*(c1+c2)*expterm
                !dxyy,dxxy,dxxz,dxzz,dyzz,dyyz
                atmp=0D0
                btmp=0D0
                ctmp=0D0
                if (ix>=2) atmp=ix*(ix-1)*sftx**(ix-2)
                if (iy>=2) btmp=iy*(iy-1)*sfty**(iy-2)
                if (iz>=2) ctmp=iz*(iz-1)*sftz**(iz-2)
                GTFdxyy=(a1+a2)*sftz**iz *expterm*(-ep4*iy*sfty**iy+epep4*sfty**(iy+2)+btmp-ep2*sfty**iy) !ok
                GTFdxxy=(b1+b2)*sftz**iz *expterm*(-ep4*ix*sftx**ix+epep4*sftx**(ix+2)+atmp-ep2*sftx**ix) !=dyxx
                GTFdxxz=(c1+c2)*sfty**iy *expterm*(-ep4*ix*sftx**ix+epep4*sftx**(ix+2)+atmp-ep2*sftx**ix) !=dzxx
                GTFdxzz=(a1+a2)*sfty**iy *expterm*(-ep4*iz*sftz**iz+epep4*sftz**(iz+2)+ctmp-ep2*sftz**iz)
                GTFdyzz=(b1+b2)*sftx**ix *expterm*(-ep4*iz*sftz**iz+epep4*sftz**(iz+2)+ctmp-ep2*sftz**iz)
                GTFdyyz=(c1+c2)*sftx**ix *expterm*(-ep4*iy*sfty**iy+epep4*sfty**(iy+2)+btmp-ep2*sfty**iy) !=dzyy,ok
                !dxxx,dyyy,dzzz
                aatmp1=0D0
                bbtmp1=0D0
                cctmp1=0D0
                if (ix>=1) aatmp1=ep2*ix*sftx**(ix-1)
                if (iy>=1) bbtmp1=ep2*iy*sfty**(iy-1)
                if (iz>=1) cctmp1=ep2*iz*sftz**(iz-1)
                aatmp2=0D0
                bbtmp2=0D0
                cctmp2=0D0
                if (ix>=2) aatmp2=ep2*ix*(ix-1)*sftx**(ix-1)
                if (iy>=2) bbtmp2=ep2*iy*(iy-1)*sfty**(iy-1)
                if (iz>=2) cctmp2=ep2*iz*(iz-1)*sftz**(iz-1)
                aatmp3=0D0
                bbtmp3=0D0
                cctmp3=0D0
                if (ix>=3) aatmp3=ix*(ix-1)*(ix-2)*sftx**(ix-3)
                if (iy>=3) bbtmp3=iy*(iy-1)*(iy-2)*sfty**(iy-3)
                if (iz>=3) cctmp3=iz*(iz-1)*(iz-2)*sftz**(iz-3)
                GTFdxxx=sfty**iy*sftz**iz*expterm*( (-2*ix+ep2*sftx2-1)*(-epep4*sftx**(ix+1) + aatmp1) - aatmp2 + epep8*sftx**(ix+1) + aatmp3 )
                GTFdyyy=sftx**ix*sftz**iz*expterm*( (-2*iy+ep2*sfty2-1)*(-epep4*sfty**(iy+1) + bbtmp1) - bbtmp2 + epep8*sfty**(iy+1) + bbtmp3 )
                GTFdzzz=sfty**iy*sftx**ix*expterm*( (-2*iz+ep2*sftz2-1)*(-epep4*sftz**(iz+1) + cctmp1) - cctmp2 + epep8*sftz**(iz+1) + cctmp3 )
                
                !Calculate 3-order derivative tensor for orbital wavefunction
                do imo=istart,iend
                    tens3(1,1,1,imo)=tens3(1,1,1,imo)+co(imo,j)*GTFdxxx !dxxx
                    tens3(2,2,2,imo)=tens3(2,2,2,imo)+co(imo,j)*GTFdyyy !dyyy
                    tens3(3,3,3,imo)=tens3(3,3,3,imo)+co(imo,j)*GTFdzzz !dzzz
                    tens3(1,2,2,imo)=tens3(1,2,2,imo)+co(imo,j)*GTFdxyy !dxyy*
                    tens3(1,1,2,imo)=tens3(1,1,2,imo)+co(imo,j)*GTFdxxy !dxxy*
                    tens3(1,1,3,imo)=tens3(1,1,3,imo)+co(imo,j)*GTFdxxz !dxxz*
                    tens3(1,3,3,imo)=tens3(1,3,3,imo)+co(imo,j)*GTFdxzz !dxzz*
                    tens3(2,3,3,imo)=tens3(2,3,3,imo)+co(imo,j)*GTFdyzz !dyzz*
                    tens3(2,2,3,imo)=tens3(2,2,3,imo)+co(imo,j)*GTFdyyz !dyyz*
                    tens3(1,2,3,imo)=tens3(1,2,3,imo)+co(imo,j)*GTFdxyz !dxyz
                end do
                tens3(1,2,1,:)=tens3(1,1,2,:) !dxyx=dxxy
                tens3(1,3,1,:)=tens3(1,1,3,:) !dxzx=dxxz
                tens3(1,3,2,:)=tens3(1,2,3,:) !dxzy=dxyz
                tens3(2,1,1,:)=tens3(1,1,2,:) !dyxx=dxxy
                tens3(2,1,2,:)=tens3(1,2,2,:) !dyxy=dxyy
                tens3(2,1,3,:)=tens3(1,2,3,:) !dyxz=dxyz
                tens3(2,2,1,:)=tens3(1,2,2,:) !dyyx=dxyy
                tens3(2,3,1,:)=tens3(1,2,3,:) !dyzx=dxyz
                tens3(2,3,2,:)=tens3(2,2,3,:) !dyzy=dyyz
                tens3(3,1,1,:)=tens3(1,1,3,:) !dzxx=dxxz
                tens3(3,1,2,:)=tens3(1,2,3,:) !dzxy=dxyz
                tens3(3,1,3,:)=tens3(1,3,3,:) !dzxz=dxzz
                tens3(3,2,1,:)=tens3(1,2,3,:) !dzyx=dxyz
                tens3(3,2,2,:)=tens3(2,2,3,:) !dzyy=dyyz
                tens3(3,2,3,:)=tens3(2,3,3,:) !dzyz=dyzz
                tens3(3,3,1,:)=tens3(1,3,3,:) !dzzx=dxzz
                tens3(3,3,2,:)=tens3(2,3,3,:) !dzzy=dyzz
            end if !end runtype>=5
            
        end if !end runtype>=3
    end if !end runtype>=2
end do
end subroutine

!!!----------- Calculate contribution from EDFs (recorded in wfx file) to density and corresponding derivatives (up to third-order)
!Only S-type GTFs are supported
! In wfx files, GTFs are used to expand core density
! runtype=1: Only calculate rho, =2: rho+dx/dy/dz =3: rho+dx/dy/dz+dxx/dyy/dzz
!        =4: rho+dx/dy/dz+full Hessian =5: rho+dx/dy/dz+full Hessian+tens3
subroutine EDFrho(runtype,x,y,z,value,grad,hess,tens3)
integer runtype
real*8 x,y,z,value
real*8,optional :: grad(3),hess(3,3),tens3(3,3,3)
value=0D0
if (present(grad)) grad=0D0
if (present(hess)) hess=0D0
if (present(tens3)) tens3=0D0

do i=1,nEDFprims
    sftx=x-a(b_EDF(i)%center)%x
    sfty=y-a(b_EDF(i)%center)%y
    sftz=z-a(b_EDF(i)%center)%z
    sftx2=sftx*sftx
    sfty2=sfty*sfty
    sftz2=sftz*sftz
    rr=sftx2+sfty2+sftz2
    ep=b_EDF(i)%exp
    expterm=exp(-ep*rr)
    value=value+CO_EDF(i)*expterm
!     write(11,"(i5,3D20.10)") i,value,CO_EDF(i),expterm
    if (runtype>=2) then
        tmp=2*CO_EDF(i)*expterm*ep
        grad(1)=grad(1)-tmp*sftx
        grad(2)=grad(2)-tmp*sfty
        grad(3)=grad(3)-tmp*sftz
        if (runtype>=3) then
            hess(1,1)=hess(1,1)+tmp*(2*ep*sftx2-1)
            hess(2,2)=hess(2,2)+tmp*(2*ep*sfty2-1)
            hess(3,3)=hess(3,3)+tmp*(2*ep*sftz2-1)
            if (runtype>=4) then
                epep4=ep*ep*4
                tmp2=CO_EDF(i)*epep4*expterm
                hess(1,2)=hess(1,2)+tmp2*sftx*sfty
                hess(1,3)=hess(1,3)+tmp2*sftx*sftz
                hess(2,3)=hess(2,3)+tmp2*sfty*sftz
                hess(2,1)=hess(1,2)
                hess(3,1)=hess(1,3)
                hess(3,2)=hess(2,3)
                if (runtype>=5) then
                    tmp3=CO_EDF(i)*epep4*expterm
                    tens3(1,1,1)=tens3(1,1,1)+tmp3*sftx*(3-2*ep*sftx2)
                    tens3(2,2,2)=tens3(2,2,2)+tmp3*sfty*(3-2*ep*sfty2)
                    tens3(3,3,3)=tens3(3,3,3)+tmp3*sftz*(3-2*ep*sftz2)
                    tens3(1,2,2)=tens3(1,2,2)+tmp3*sftx*(1-2*ep*sfty2)
                    tens3(1,1,2)=tens3(1,1,2)+tmp3*sfty*(1-2*ep*sftx2)
                    tens3(1,1,3)=tens3(1,1,3)+tmp3*sftz*(1-2*ep*sftx2)
                    tens3(1,3,3)=tens3(1,3,3)+tmp3*sftx*(1-2*ep*sftz2)
                    tens3(2,3,3)=tens3(2,3,3)+tmp3*sfty*(1-2*ep*sftz2)
                    tens3(2,2,3)=tens3(2,2,3)+tmp3*sftz*(1-2*ep*sfty2)
                    tens3(1,2,3)=tens3(1,2,3)-CO_EDF(i)*8*ep**3*sftx*sfty*sftz*expterm
                    tens3(1,2,1)=tens3(1,1,2) !dxyx=dxxy
                    tens3(1,3,1)=tens3(1,1,3) !dxzx=dxxz
                    tens3(1,3,2)=tens3(1,2,3) !dxzy=dxyz
                    tens3(2,1,1)=tens3(1,1,2) !dyxx=dxxy
                    tens3(2,1,2)=tens3(1,2,2) !dyxy=dxyy
                    tens3(2,1,3)=tens3(1,2,3) !dyxz=dxyz
                    tens3(2,2,1)=tens3(1,2,2) !dyyx=dxyy
                    tens3(2,3,1)=tens3(1,2,3) !dyzx=dxyz
                    tens3(2,3,2)=tens3(2,2,3) !dyzy=dyyz
                    tens3(3,1,1)=tens3(1,1,3) !dzxx=dxxz
                    tens3(3,1,2)=tens3(1,2,3) !dzxy=dxyz
                    tens3(3,1,3)=tens3(1,3,3) !dzxz=dxzz
                    tens3(3,2,1)=tens3(1,2,3) !dzyx=dxyz
                    tens3(3,2,2)=tens3(2,2,3) !dzyy=dyyz
                    tens3(3,2,3)=tens3(2,3,3) !dzyz=dyzz
                    tens3(3,3,1)=tens3(1,3,3) !dzzx=dxzz
                    tens3(3,3,2)=tens3(2,3,3) !dzzy=dyzz
                end if
            end if
        end if
    end if
end do
end subroutine

!!!------ A general routine used to calculate value, gradient and Hessian matrix at a given point for some real space functions
! itype=1 Only calculate value and grad
! itype=2 Calculate value, gradient and Hessian
subroutine gencalchessmat(itype,ifunc,x,y,z,value,grad,hess)
integer ifunc,itype
real*8 x,y,z,value,grad(3),hess(3,3)
real*8 gradaddx(3),gradminx(3),gradaddy(3),gradminy(3),gradaddz(3),gradminz(3)
character selELFLOL*3
diff=8D-4
denom=2D0*diff
if (ifunc==9) selELFLOL="ELF"
if (ifunc==10) selELFLOL="LOL"

!Evaluate both gradient and Hessian analytically then return directly
if (ifunc==1) then
    call calchessmat_dens(itype,x,y,z,value,grad,hess)
    return
else if (ifunc==4) then
    call calchessmat_mo(itype,iorbsel,x,y,z,value,grad,hess)
    return
end if
    
!For other functions, analytical Hessian or even gradient hasn't been realized
!Now calculate gradient
if (ifunc==3) then
    call calchessmat_lapl(1,x,y,z,value,grad,hess)
else if (ifunc==9.or.ifunc==10) then
    if (ELFLOL_type==0) then !Analytical gradient for Becke's ELF/LOL
        call calchessmat_ELF_LOL(1,x,y,z,value,grad,hess,selELFLOL)
    else !Numerical gradient for other definition of ELF/LOL
        value=ELF_LOL(x,y,z,selELFLOL)
        xadd=ELF_LOL(x+diff,y,z,selELFLOL)
        xmin=ELF_LOL(x-diff,y,z,selELFLOL)
        yadd=ELF_LOL(x,y+diff,z,selELFLOL)
        ymin=ELF_LOL(x,y-diff,z,selELFLOL)
        zadd=ELF_LOL(x,y,z+diff,selELFLOL)
        zmin=ELF_LOL(x,y,z-diff,selELFLOL)
        grad(1)=(xadd-xmin)/denom
        grad(2)=(yadd-ymin)/denom
        grad(3)=(zadd-zmin)/denom
    end if
else if (ifunc==12) then
    value=totesp(x,y,z)
    xadd=totesp(x+diff,y,z)
    xmin=totesp(x-diff,y,z)
    yadd=totesp(x,y+diff,z)
    ymin=totesp(x,y-diff,z)
    zadd=totesp(x,y,z+diff)
    zmin=totesp(x,y,z-diff)
    grad(1)=(xadd-xmin)/denom
    grad(2)=(yadd-ymin)/denom
    grad(3)=(zadd-zmin)/denom
else if (ifunc==100) then
    value=userfunc(x,y,z)
    xadd=userfunc(x+diff,y,z)
    xmin=userfunc(x-diff,y,z)
    yadd=userfunc(x,y+diff,z)
    ymin=userfunc(x,y-diff,z)
    zadd=userfunc(x,y,z+diff)
    zmin=userfunc(x,y,z-diff)
    grad(1)=(xadd-xmin)/denom
    grad(2)=(yadd-ymin)/denom
    grad(3)=(zadd-zmin)/denom
end if

!Calculate Hessian semi (namely based on analyical gradient) or pure (based on function value) numerically
if (itype==2) then
    if (ifunc==3) then !Use semi-analytical for Hessian of laplacian
        call calchessmat_lapl(1,x+diff,y,z,tmpval,gradaddx,hess)
        call calchessmat_lapl(1,x-diff,y,z,tmpval,gradminx,hess)
        call calchessmat_lapl(1,x,y+diff,z,tmpval,gradaddy,hess)
        call calchessmat_lapl(1,x,y-diff,z,tmpval,gradminy,hess)
        call calchessmat_lapl(1,x,y,z+diff,tmpval,gradaddz,hess)
        call calchessmat_lapl(1,x,y,z-diff,tmpval,gradminz,hess)
        hess(1,1)=(gradaddx(1)-gradminx(1))/denom
        hess(2,2)=(gradaddy(2)-gradminy(2))/denom
        hess(3,3)=(gradaddz(3)-gradminz(3))/denom
        hess(1,2)=(gradaddy(1)-gradminy(1))/denom
        hess(2,3)=(gradaddz(2)-gradminz(2))/denom
        hess(1,3)=(gradaddz(1)-gradminz(1))/denom
        hess(2,1)=hess(1,2)
        hess(3,2)=hess(2,3)
        hess(3,1)=hess(1,3)
        return !Don't do below procedures for generating pure numerical Hessian
    else if (ifunc==9.or.ifunc==10) then
        if (ELFLOL_type==0) then !Use semi-analytical for Hessian of Becke's ELF/LOL
            call calchessmat_ELF_LOL(1,x+diff,y,z,tmpval,gradaddx,hess,selELFLOL)
            call calchessmat_ELF_LOL(1,x-diff,y,z,tmpval,gradminx,hess,selELFLOL)
            call calchessmat_ELF_LOL(1,x,y+diff,z,tmpval,gradaddy,hess,selELFLOL)
            call calchessmat_ELF_LOL(1,x,y-diff,z,tmpval,gradminy,hess,selELFLOL)
            call calchessmat_ELF_LOL(1,x,y,z+diff,tmpval,gradaddz,hess,selELFLOL)
            call calchessmat_ELF_LOL(1,x,y,z-diff,tmpval,gradminz,hess,selELFLOL)
            hess(1,1)=(gradaddx(1)-gradminx(1))/denom
            hess(2,2)=(gradaddy(2)-gradminy(2))/denom
            hess(3,3)=(gradaddz(3)-gradminz(3))/denom
            hess(1,2)=(gradaddy(1)-gradminy(1))/denom
            hess(2,3)=(gradaddz(2)-gradminz(2))/denom
            hess(1,3)=(gradaddz(1)-gradminz(1))/denom
            hess(2,1)=hess(1,2)
            hess(3,2)=hess(2,3)
            hess(3,1)=hess(1,3)
            return
        else !for other definition of ELF/LOL, use pure numerical Hessian
            xaddxadd=ELF_LOL(x+2D0*diff,y,z,selELFLOL)
            xminxmin=ELF_LOL(x-2D0*diff,y,z,selELFLOL)
            yaddyadd=ELF_LOL(x,y+2D0*diff,z,selELFLOL)
            yminymin=ELF_LOL(x,y-2D0*diff,z,selELFLOL)
            zaddzadd=ELF_LOL(x,y,z+2D0*diff,selELFLOL)
            zminzmin=ELF_LOL(x,y,z-2D0*diff,selELFLOL)
            xaddyadd=ELF_LOL(x+diff,y+diff,z,selELFLOL)
            xminyadd=ELF_LOL(x-diff,y+diff,z,selELFLOL)
            xaddymin=ELF_LOL(x+diff,y-diff,z,selELFLOL)
            xminymin=ELF_LOL(x-diff,y-diff,z,selELFLOL)
            yaddzadd=ELF_LOL(x,y+diff,z+diff,selELFLOL)
            yminzadd=ELF_LOL(x,y-diff,z+diff,selELFLOL)
            yaddzmin=ELF_LOL(x,y+diff,z-diff,selELFLOL)
            yminzmin=ELF_LOL(x,y-diff,z-diff,selELFLOL)
            xaddzadd=ELF_LOL(x+diff,y,z+diff,selELFLOL)
            xminzadd=ELF_LOL(x-diff,y,z+diff,selELFLOL)
            xaddzmin=ELF_LOL(x+diff,y,z-diff,selELFLOL)
            xminzmin=ELF_LOL(x-diff,y,z-diff,selELFLOL)
        end if
    else if (ifunc==12) then !pure numerical Hessian for total ESP
        xaddxadd=totesp(x+2D0*diff,y,z)
        xminxmin=totesp(x-2D0*diff,y,z)
        yaddyadd=totesp(x,y+2D0*diff,z)
        yminymin=totesp(x,y-2D0*diff,z)
        zaddzadd=totesp(x,y,z+2D0*diff)
        zminzmin=totesp(x,y,z-2D0*diff)
        xaddyadd=totesp(x+diff,y+diff,z)
        xminyadd=totesp(x-diff,y+diff,z)
        xaddymin=totesp(x+diff,y-diff,z)
        xminymin=totesp(x-diff,y-diff,z)
        yaddzadd=totesp(x,y+diff,z+diff)
        yminzadd=totesp(x,y-diff,z+diff)
        yaddzmin=totesp(x,y+diff,z-diff)
        yminzmin=totesp(x,y-diff,z-diff)
        xaddzadd=totesp(x+diff,y,z+diff)
        xminzadd=totesp(x-diff,y,z+diff)
        xaddzmin=totesp(x+diff,y,z-diff)
        xminzmin=totesp(x-diff,y,z-diff)
    else if (ifunc==100) then !pure numerical Hessian for user defined function
        xaddxadd=userfunc(x+2D0*diff,y,z)
        xminxmin=userfunc(x-2D0*diff,y,z)
        yaddyadd=userfunc(x,y+2D0*diff,z)
        yminymin=userfunc(x,y-2D0*diff,z)
        zaddzadd=userfunc(x,y,z+2D0*diff)
        zminzmin=userfunc(x,y,z-2D0*diff)
        xaddyadd=userfunc(x+diff,y+diff,z)
        xminyadd=userfunc(x-diff,y+diff,z)
        xaddymin=userfunc(x+diff,y-diff,z)
        xminymin=userfunc(x-diff,y-diff,z)
        yaddzadd=userfunc(x,y+diff,z+diff)
        yminzadd=userfunc(x,y-diff,z+diff)
        yaddzmin=userfunc(x,y+diff,z-diff)
        yminzmin=userfunc(x,y-diff,z-diff)
        xaddzadd=userfunc(x+diff,y,z+diff)
        xminzadd=userfunc(x-diff,y,z+diff)
        xaddzmin=userfunc(x+diff,y,z-diff)
        xminzmin=userfunc(x-diff,y,z-diff)
    end if 
    !Collect above temporary data to evaluate pure numerical Hessian
    gradx_yadd=(xaddyadd-xminyadd)/denom
    gradx_ymin=(xaddymin-xminymin)/denom
    grady_zadd=(yaddzadd-yminzadd)/denom
    grady_zmin=(yaddzmin-yminzmin)/denom
    gradx_zadd=(xaddzadd-xminzadd)/denom
    gradx_zmin=(xaddzmin-xminzmin)/denom
    hess(1,1)=(xaddxadd-2*value+xminxmin)/denom**2
    hess(2,2)=(yaddyadd-2*value+yminymin)/denom**2
    hess(3,3)=(zaddzadd-2*value+zminzmin)/denom**2
    hess(1,2)=(gradx_yadd-gradx_ymin)/denom
    hess(2,3)=(grady_zadd-grady_zmin)/denom
    hess(1,3)=(gradx_zadd-gradx_zmin)/denom
    hess(2,1)=hess(1,2)
    hess(3,2)=hess(2,3)
    hess(3,1)=hess(1,3)
end if
end subroutine



!============================================================
!=================== Real space functions ===================
!============================================================


!!!--------- User defined function, the content is needed to be filled by users or selected by iuserfunc
!The units should be given in a.u.
real*8 function userfunc(x,y,z)
real*8 x,y,z
userfunc=1D0 !Default value. Note: default "iuserfunc" is 0
!Below functions can be selected by "iuserfunc" parameter in settings.ini
if (iuserfunc==-2) userfunc=calcprodens(x,y,z,0) !Promolecular density
if (iuserfunc==-1) userfunc=linintp3d(x,y,z,1) !The function value evaluated by trilinear interpolation from cubmat
if (iuserfunc==1) userfunc=fspindens(x,y,z,'a') !Alpha density
if (iuserfunc==2) userfunc=fspindens(x,y,z,'b') !Beta density
if (iuserfunc==3) userfunc=(x*x+y*y+z*z)*fdens(x,y,z) !Integrand of electronic spatial extent <R**2>
if (iuserfunc==4) userfunc=weizpot(x,y,z) !Weizsacker potential
if (iuserfunc==5) userfunc=weizsacker(x,y,z) !Integrand of weizsacker functional
if (iuserfunc==6) userfunc=4*pi*fdens(x,y,z)*(x*x+y*y+z*z) !Radial distribution function (assume that density is sphericalized)
if (iuserfunc==7) userfunc=2D0/3D0*lagkin(x,y,z,0)/fdens(x,y,z) !Local Temperature(Kelvin), PNAS,81,8028
if (iuserfunc==8) userfunc=totesp(x,y,z)/fdens(x,y,z) !Average local electrostatic potential, useful to exhibit atomic shell structure, see Chapter 8 of Theoretical Aspects of Chemical Reactivity
if (iuserfunc==9) userfunc=fdens(x,y,z)/nelec !Shape function
if (iuserfunc==10) userfunc=-Hamkin(x,y,z,0)-lagkin(x,y,z,0) !Potential energy density, also known as virial field
if (iuserfunc==11) userfunc=-Hamkin(x,y,z,0) !Energy density
if (iuserfunc==12) userfunc=-nucesp(x,y,z)*fdens(x,y,z) !Local nuclear attraction potential energy
if (iuserfunc==13) userfunc=lagkin(x,y,z,0)/fdens(x,y,z) !This quantity at bond critical point is useful to discriminate covalent bonding and close-shell interaction
if (iuserfunc==14) userfunc=eleesp(x,y,z) !Electrostatic potential from electrons
if (iuserfunc==15) userfunc=fdens(x,y,z)/flapl(x,y,z,'t') !Bond metallicity
if (iuserfunc==16) userfunc=36*(3*pi*pi)**(2D0/3D0)/5D0*fdens(x,y,z)**(5D0/3D0)/flapl(x,y,z,'t') !Dimensionless bond metallicity
if (iuserfunc==17) userfunc=-Hamkin(x,y,z,0)/fdens(x,y,z) !Energy density per electron
if (iuserfunc==18) then !Region of Slow Electrons (RoSE), defined in Chem. Phys. Lett., 582, 144 (2013)
    rho=fdens(x,y,z)
    if (wfntype==0.or.wfntype==3) then !close shell cases
        Dh=2.871234000D0*rho**(5.0D0/3.0D0)
    else if (wfntype==1.or.wfntype==2.or.wfntype==4) then !Open shell cases
        rhospin=fspindens(x,y,z,'s') !rhospin=rhoa-rhob, rho=rhoa+rhob
        rhoa=(rhospin+rho)/2D0
        rhob=(rho-rhospin)/2D0
        Dh=4.557799872D0*(rhoa**(5.0D0/3.0D0)+rhob**(5.0D0/3.0D0)) !kinetic energy of HEG
    end if
    G=Lagkin(x,y,z,0)
    userfunc=(Dh-G)/(Dh+G)
end if
if (iuserfunc==19) userfunc=SEDD(x,y,z) !SEDD
if (iuserfunc==20) userfunc=DORI(x,y,z) !DORI
if (iuserfunc==21) userfunc=-x*fdens(x,y,z) !Integrand of X component of electric dipole moment
if (iuserfunc==22) userfunc=-y*fdens(x,y,z) !Integrand of Y component of electric dipole moment
if (iuserfunc==23) userfunc=-z*fdens(x,y,z) !Integrand of Z component of electric dipole moment
if (iuserfunc==24) userfunc=linrespkernel(x,y,z) !Approximate form of DFT linear response kernel for close-shell
if (iuserfunc==25) userfunc=fgrad(x,y,z,'t')/fdens(x,y,z)/2D0 !Magnitude of fluctuation of the electronic momentum
if (iuserfunc==26) userfunc=2.871234D0*rho**(5D0/3D0) !Thomas-Fermi kinetic energy density
if (iuserfunc==27) userfunc=loceleaff(x,y,z) !Local electron affinity
if (iuserfunc==28) userfunc=(avglocion(x,y,z)+loceleaff(x,y,z))/2 !Local Mulliken electronegativity
if (iuserfunc==29) userfunc=(avglocion(x,y,z)-loceleaff(x,y,z))/2 !Local hardness
if (iuserfunc==30) userfunc=densellip(x,y,z,1) !Ellipticity of electron density
if (iuserfunc==31) userfunc=densellip(x,y,z,2) !eta index, Angew. Chem. Int. Ed., 53, 2766-2770 (2014)
if (iuserfunc==32) userfunc=densellip(x,y,z,2)-1 !Modified eta index
if (iuserfunc==33) userfunc=PAEM(x,y,z,1) !PAEM, potential acting on one electron in a molecule, defined by Zhongzhi Yang
if (iuserfunc==34) userfunc=PAEM(x,y,z,2) !The same as 33, but using DFT XC potential directly rather than evaluating the XC potential based on pair density
if (iuserfunc==35) then !|V(r)|/G(r)
    tmpval=lagkin(x,y,z,0)
    userfunc=abs(-Hamkin(x,y,z,0)-tmpval)/tmpval
end if
if (iuserfunc==36) then !On-top pair density, i.e. r1=r2 case of pair density. paircorrtype affects result
    pairfunctypeold=pairfunctype
    pairfunctype=12
    userfunc=pairfunc(x,y,z,x,y,z)
    pairfunctype=pairfunctypeold
end if
if (iuserfunc==38) userfunc=Ang_rhoeigvec_ple(x,y,z,2) !The angle between the second eigenvector of rho and the plane defined by option 4 of main function 1000
if (iuserfunc==39) userfunc=totespskip(x,y,z,iskipnuc) !ESP without contribution of nuclues "iskipnuc"
if (iuserfunc==40) userfunc=weizsacker(x,y,z) !Steric energy
if (iuserfunc==41) userfunc=stericpot(x,y,z) !Steric potential
if (iuserfunc==42) userfunc=stericcharge(x,y,z) !Steric charge
if (iuserfunc==43) userfunc=stericforce(x,y,z) !The magnitude of steric force
if (iuserfunc==44) userfunc=stericpot_damp(x,y,z) !Steric potential with damping function to a given constant value
if (iuserfunc==45) userfunc=stericforce_damp(x,y,z) !Steric force based on damped potential
if (iuserfunc==46) userfunc=stericforce_directdamp(x,y,z) !Steric force directly damped to zero
if (iuserfunc==50) userfunc=infoentro(2,x,y,z) !Shannon entropy density, see JCP,126,191107 for example
if (iuserfunc==51) userfunc=Fisherinfo(1,x,y,z) !Fisher information density, see JCP,126,191107 for example
if (iuserfunc==52) userfunc=Fisherinfo(2,x,y,z) !Second Fisher information density, see JCP,126,191107 for derivation
if (iuserfunc==53) userfunc=Ghoshentro(x,y,z,1) !Ghosh entropy density with G(r) as kinetic energy density, PNAS,81,8028
if (iuserfunc==54) userfunc=Ghoshentro(x,y,z,2) !Ghosh entropy density with G(r)-der2rho/8 as kinetic energy density, exactly corresponds to Eq.22 in PNAS,81,8028
if (iuserfunc==55) userfunc=fdens(x,y,z)**2 !Integrand of quadratic form of Renyi entropy
if (iuserfunc==56) userfunc=fdens(x,y,z)**3 !Integrand of cubic form of Renyi entropy
if (iuserfunc==60) userfunc=paulipot(x,y,z) !Pauli potential, Comp. Theor. Chem., 1006, 92-99
if (iuserfunc==61) userfunc=pauliforce(x,y,z) !The magnitude of Pauli force
if (iuserfunc==62) userfunc=paulicharge(x,y,z) !Pauli charge
if (iuserfunc==63) userfunc=quantumpot(x,y,z) !Quantum potential
if (iuserfunc==64) userfunc=quantumforce(x,y,z) !The magnitude of quantum force
if (iuserfunc==65) userfunc=quantumcharge(x,y,z) !Quantum charge
if (iuserfunc==66) userfunc=elestatforce(x,y,z) !The magnitude of electrostatic force
if (iuserfunc==67) userfunc=elestatcharge(x,y,z) !Electrostatic charge
if (iuserfunc==70) userfunc=4.5D0*fdens(x,y,z)**2/lagkin(x,y,z,0)   !Phase-space-defined Fisher information density
if (iuserfunc==71) userfunc=elemomdens(x,y,z,1) !X component of electron linear momentum density in 3D representation
if (iuserfunc==72) userfunc=elemomdens(x,y,z,2) !Y component of electron linear momentum density in 3D representation
if (iuserfunc==73) userfunc=elemomdens(x,y,z,3) !Z component of electron linear momentum density in 3D representation
if (iuserfunc==74) userfunc=elemomdens(x,y,z,0) !Magnitude of electron linear momentum density in 3D representation
if (iuserfunc==75) userfunc=magmomdens(x,y,z,1) !X component of magnetic dipole moment density
if (iuserfunc==76) userfunc=magmomdens(x,y,z,2) !Y component of magnetic dipole moment density
if (iuserfunc==77) userfunc=magmomdens(x,y,z,3) !Z component of magnetic dipole moment density
if (iuserfunc==78) userfunc=magmomdens(x,y,z,0) !Magnitude of magnetic dipole moment density
if (iuserfunc==81) userfunc=hamkin(x,y,z,1) !X component of Hamiltonian kinetic energy density
if (iuserfunc==82) userfunc=hamkin(x,y,z,2) !Y component of Hamiltonian kinetic energy density
if (iuserfunc==83) userfunc=hamkin(x,y,z,3) !Z component of Hamiltonian kinetic energy density
if (iuserfunc==84) userfunc=Lagkin(x,y,z,1) !X component of Lagrangian kinetic energy density
if (iuserfunc==85) userfunc=Lagkin(x,y,z,2) !Y component of Lagrangian kinetic energy density
if (iuserfunc==86) userfunc=Lagkin(x,y,z,3) !Z component of Lagrangian kinetic energy density
if (iuserfunc==87) userfunc=localcorr(x,y,z,1) !Local total electron correlation function
if (iuserfunc==88) userfunc=localcorr(x,y,z,2) !Local dynamic electron correlation function
if (iuserfunc==89) userfunc=localcorr(x,y,z,3) !Local nondynamic electron correlation function
if (iuserfunc==90) then
    tmpELF=ELF_LOL(x,y,z,"ELF")
    userfunc=tmpELF*tmpELF*(x*x+y*y+z*z)
else if (iuserfunc==91) then
    userfunc=tmpELF*tmpELF
end if
if (iuserfunc==100) userfunc=fdens(x,y,z)**2 !Disequilibrium (also known as semi-similarity), DOI: 10.1002/qua.24510
if (iuserfunc==101) then !Positive part of ESP
    userfunc=totesp(x,y,z)
    if (userfunc<0D0) userfunc=0D0
else if (iuserfunc==102) then !Negative part of ESP
    userfunc=totesp(x,y,z)
    if (userfunc>0D0) userfunc=0D0
end if
if (iuserfunc>=800) userfunc=funcvalLSB(x,y,z,iuserfunc-800)
if (iuserfunc==900) userfunc=x !X coordinate
if (iuserfunc==901) userfunc=y !Y coordinate
if (iuserfunc==902) userfunc=z !Z coordinate
if (iuserfunc==1000) userfunc=DFTxcfunc(x,y,z) !Various kinds of DFT exchange-correlation functions
if (iuserfunc==1100) userfunc=DFTxcpot(x,y,z) !Various kinds of DFT exchange-correlation potentials
!Below are other examples
! userfunc=hamkin(x,y,z,3)-0.5D0*(hamkin(x,y,z,1)+hamkin(x,y,z,2)) !Anisotropy of Hamiltonian kinetic energy in Z, namely K_Z-0.5*(K_X+K_Y)
! userfunc=-x*y*fdens(x,y,z) !Integrand of XY component of electric quadrupole moment
! userfunc=-x*y*z*fdens(x,y,z)*au2debye*b2a**2 !Integrand of XYZ component of electric octapole moment in Debye-Ang**2
end function


!!!----------- Output orbital wavefunction value at a given point (fmo=function for outputting MO)
real*8 function fmo(x,y,z,id)
real*8 x,y,z,orbval(nmo)
integer id
call orbderv(1,id,id,x,y,z,orbval)
fmo=orbval(id)
end function

!!!----- Calculate orbital wavefunction Hessian matrix at x,y,z, store to hess, output value and gradient vector at same time
!itype=1 Only calculate value and gradient, not Hessian
!itype=2 Calculate value, gradient and Hessian
subroutine calchessmat_mo(itype,id,x,y,z,value,grad,hess)
integer id,itype
real*8 x,y,z,grad(3),hess(3,3),value,wfnval(nmo),wfnderv(3,nmo),wfnhess(3,3,nmo)
if (itype==1) call orbderv(2,id,id,x,y,z,wfnval,wfnderv,wfnhess)
if (itype==2) call orbderv(4,id,id,x,y,z,wfnval,wfnderv,wfnhess)
value=wfnval(id)
grad=wfnderv(:,id)
hess=wfnhess(:,:,id)
end subroutine


!--------Output electron density at a point
real*8 function fdens(x,y,z)
real*8 x,y,z,wfnval(nmo)
call orbderv(1,1,nmo,x,y,z,wfnval)
fdens=0D0
do i=1,nmo
    fdens=fdens+MOocc(i)*wfnval(i)**2
end do
! Add in contribution of Electron density function
if (allocated(b_EDF)) then
    call EDFrho(1,x,y,z,EDFdens)
    fdens=fdens+EDFdens
end if
! if (fdens>0.5D0) fdens=0
end function


!!!------------------------- Output spin or Alpha or Beta electron density at a point
!itype='s' output spin density, ='a' output alpha density, ='b' output beta density
real*8 function fspindens(x,y,z,itype)
real*8 :: x,y,z,wfnval(nmo)
character itype
call orbderv(1,1,nmo,x,y,z,wfnval)
adens=0.0D0
bdens=0.0D0
do i=1,nmo
    if (MOtype(i)==1) then
        adens=adens+MOocc(i)*wfnval(i)**2
    else if (MOtype(i)==2) then
        bdens=bdens+MOocc(i)*wfnval(i)**2
    else if (MOtype(i)==0) then
        adens=adens+MOocc(i)/2D0*wfnval(i)**2
        bdens=bdens+MOocc(i)/2D0*wfnval(i)**2
    end if
end do
if (itype=='s') then
    fspindens=adens-bdens
    if (ipolarpara==1) fspindens=fspindens/(adens+bdens)
else if (itype=='a') then
    fspindens=adens
else if (itype=='b') then
    fspindens=bdens
end if
end function


!!!------------------------- Output gradient of rho and RDG(reduced density gradient) at a point
!label=x/y/z output 1-order derivation of x/y/z, =t get norm, =r get RDG, =s get |der_rho|/rho^(4/3)
real*8 function fgrad(x,y,z,label)
real*8 x,y,z,wfnval(nmo),wfnderv(3,nmo),gradrho(3),EDFgrad(3),sumgrad2
character label
call orbderv(2,1,nmo,x,y,z,wfnval,wfnderv)
rho=0D0
gradrho=0D0
do i=1,nmo
    rho=rho+MOocc(i)*wfnval(i)**2
    gradrho(:)=gradrho(:)+MOocc(i)*wfnval(i)*wfnderv(:,i)
end do
gradrho=2*gradrho
! Add in contribution of Electron density function
if (allocated(b_EDF)) then
    call EDFrho(2,x,y,z,EDFdens,EDFgrad)
    rho=rho+EDFdens
    gradrho=gradrho+EDFgrad
end if
if (label=='x') then
    fgrad=gradrho(1)
else if (label=='y') then
    fgrad=gradrho(2)
else if (label=='z') then
    fgrad=gradrho(3)
else if (label=='t') then
    fgrad=dsqrt( sum(gradrho(:)**2) )
else if (label=='r') then
    sumgrad2=sum(gradrho(:)**2)
    if (RDG_maxrho/=0.0D0.and.rho>=RDG_maxrho) then
        fgrad=100D0
!This occurs at distant region when exponent cutoff is used, the actual value should be very large. In order to avoid denominator become zero, we set it artifically to a big value
    else if (sumgrad2==0D0.or.rho==0D0) then
        RDG=999D0
    else
        fgrad=0.161620459673995D0*dsqrt(sumgrad2)/rho**(4.0D0/3.0D0) !0.161620459673995D0=1/(2*(3*pi**2)**(1/3))
    end if
else if (label=='s') then
    if (rho==0D0) then
        fgrad=0
    else
        fgrad=dsqrt(sum(gradrho(:)**2))/rho**(4.0D0/3.0D0)
    end if
end if
end function


!!--- Simultaneously generate electron density, gradient norm for alpha and beta electrons, as well as dot product between grada and gradb
!---- Mainly used to evalute DFT functional. EDF is not taken into account
!adens/bdens/tdens means the density of alpha/beta/total density, similar for *grad
subroutine gendensgradab(x,y,z,adens,bdens,tdens,agrad,bgrad,tgrad,abgrad)
real*8 x,y,z,adens,bdens,agrad,bgrad,abgrad,wfnval(nmo),wfnderv(3,nmo),gradrhoa(3),gradrhob(3),gradrhot(3),tmparr(3)
call orbderv(2,1,nmo,x,y,z,wfnval,wfnderv)
adens=0D0
bdens=0D0
gradrhoa=0D0
gradrhob=0D0
do i=1,nmo
    if (MOtype(i)==1) then
        adens=adens+MOocc(i)*wfnval(i)**2
        gradrhoa(:)=gradrhoa(:)+MOocc(i)*wfnval(i)*wfnderv(:,i)
    else if (MOtype(i)==2) then
        bdens=bdens+MOocc(i)*wfnval(i)**2
        gradrhob(:)=gradrhob(:)+MOocc(i)*wfnval(i)*wfnderv(:,i)
    else if (MOtype(i)==0) then
        tmpval=MOocc(i)/2D0*wfnval(i)**2
        adens=adens+tmpval
        bdens=bdens+tmpval
        tmparr(:)=MOocc(i)/2D0*wfnval(i)*wfnderv(:,i)
        gradrhoa(:)=gradrhoa(:)+tmparr(:)
        gradrhob(:)=gradrhob(:)+tmparr(:)
    end if
end do
tdens=adens+bdens
gradrhoa=gradrhoa*2
gradrhob=gradrhob*2
gradrhot=gradrhoa+gradrhob
agrad=dsqrt(sum(gradrhoa**2))
bgrad=dsqrt(sum(gradrhob**2))
tgrad=dsqrt(sum(gradrhot**2))
abgrad=sum(gradrhoa*gradrhob)
end subroutine


!!!------------------------- Output Laplacian of electron density at a point
!label=x/y/z output 2-order derivative of electron density respect to xx/yy/zz; &
!=t get their summing; =s get der2rho/rho^(5/3), which is used LSB's project
real*8 function flapl(x,y,z,label)
real*8 x,y,z,wfnval(nmo),wfnderv(3,nmo),wfnhess(3,3,nmo),laplx,laply,laplz,EDFgrad(3),EDFhess(3,3)
character label
call orbderv(3,1,nmo,x,y,z,wfnval,wfnderv,wfnhess)
laplx=2*sum( MOocc(1:nmo)*( wfnderv(1,1:nmo)**2 + wfnval(1:nmo)*wfnhess(1,1,1:nmo) ) )
laply=2*sum( MOocc(1:nmo)*( wfnderv(2,1:nmo)**2 + wfnval(1:nmo)*wfnhess(2,2,1:nmo) ) )
laplz=2*sum( MOocc(1:nmo)*( wfnderv(3,1:nmo)**2 + wfnval(1:nmo)*wfnhess(3,3,1:nmo) ) )
!Add in contribution of electron density function, assume EDFs are S type
if (allocated(b_EDF)) then
    call EDFrho(3,x,y,z,EDFdens,EDFgrad,EDFhess)
    laplx=laplx+EDFhess(1,1)
    laply=laply+EDFhess(2,2)
    laplz=laplz+EDFhess(3,3)
end if
if (label=='t') then
    flapl=laplx+laply+laplz
    flapl=flapl*laplfac !laplfac is an external variable
else if (label=='x') then
    flapl=laplx
else if (label=='y') then
    flapl=laply
else if (label=='z') then
    flapl=laplz
else if (label=='s') then
    dens=sum(MOocc(1:nmo)*wfnval(1:nmo)**2)
    if (allocated(b_EDF)) dens=dens+EDFdens
    flapl=(laplx+laply+laplz)/dens**(5D0/3D0)
end if
end function


!!!----- Calculate electron density, its gradient and Hessian matrix
!itype=1 Only calculate value and grad, not Hessian
!itype=2 calculate value, gradient and Hessian
subroutine calchessmat_dens(itype,x,y,z,elerho,elegrad,elehess)
real*8 x,y,z,elerho,elegrad(3),elehess(3,3),wfnval(nmo),wfnderv(3,nmo),wfnhess(3,3,nmo),EDFgrad(3),EDFhess(3,3)
integer itype
call orbderv(4,1,nmo,x,y,z,wfnval,wfnderv,wfnhess)
elerho=sum( MOocc(1:nmo)*wfnval(1:nmo)*wfnval(1:nmo) )
do itmp=1,3
    elegrad(itmp)=2*sum( MOocc(1:nmo)*wfnval(1:nmo)*wfnderv(itmp,1:nmo) )
end do
if (itype==2) then
    elehess(1,1)=2*sum( MOocc(1:nmo)*( wfnderv(1,1:nmo)**2 + wfnval(1:nmo)*wfnhess(1,1,1:nmo) ) )
    elehess(2,2)=2*sum( MOocc(1:nmo)*( wfnderv(2,1:nmo)**2 + wfnval(1:nmo)*wfnhess(2,2,1:nmo) ) )
    elehess(3,3)=2*sum( MOocc(1:nmo)*( wfnderv(3,1:nmo)**2 + wfnval(1:nmo)*wfnhess(3,3,1:nmo) ) )
    elehess(1,2)=2*sum( MOocc(1:nmo)*( wfnderv(1,1:nmo)*wfnderv(2,1:nmo)+wfnhess(1,2,1:nmo)*wfnval(1:nmo) ) )
    elehess(2,3)=2*sum( MOocc(1:nmo)*( wfnderv(2,1:nmo)*wfnderv(3,1:nmo)+wfnhess(2,3,1:nmo)*wfnval(1:nmo) ) )
    elehess(1,3)=2*sum( MOocc(1:nmo)*( wfnderv(1,1:nmo)*wfnderv(3,1:nmo)+wfnhess(1,3,1:nmo)*wfnval(1:nmo) ) )
    elehess(2,1)=elehess(1,2)
    elehess(3,2)=elehess(2,3)
    elehess(3,1)=elehess(1,3)
end if

!Add in contribution of electron density function, assume EDFs are S type
if (allocated(b_EDF)) then
    call EDFrho(4,x,y,z,EDFdens,EDFgrad,EDFhess)
    elerho=elerho+EDFdens
    elegrad=elegrad+EDFgrad
    elehess=elehess+EDFhess
end if
end subroutine


!!------------- Calculate Laplacian of electron density, its gradient and Hessian matrix
!itype=1 calculate value, gradient
!itype=2 calculate value, gradient and Hessian (Not available)
subroutine calchessmat_lapl(itype,x,y,z,value,grad,hess)
use util
real*8 x,y,z,value,grad(3),hess(3,3)
real*8 wfnval(nmo),wfnderv(3,nmo),wfnhess(3,3,nmo),wfntens3(3,3,3,nmo),rhotens3(3,3,3)
real*8 EDFgrad(3),EDFhess(3,3),EDFtens3(3,3,3)
integer itype
!Numerically verify 3-order derivative of orbital wavefunction
! diff=1D-5
! call orbderv(5,1,nmo,x+diff,y,z,wfnval,wfnderv,wfnhess,wfntens3)
! t1=wfnhess(3,3,1)
! call orbderv(5,1,nmo,x-diff,y,z,wfnval,wfnderv,wfnhess,wfntens3)
! t2=wfnhess(3,3,1)
! write(*,*) (t1-t2)/(2*diff)

call orbderv(5,1,nmo,x,y,z,wfnval,wfnderv,wfnhess,wfntens3)
rhotens3=0D0
dxx=0D0
dyy=0D0
dzz=0D0
do i=1,nmo
    dxx=dxx+MOocc(i)*(wfnderv(1,i)**2+wfnval(i)*wfnhess(1,1,i))
    dyy=dyy+MOocc(i)*(wfnderv(2,i)**2+wfnval(i)*wfnhess(2,2,i))
    dzz=dzz+MOocc(i)*(wfnderv(3,i)**2+wfnval(i)*wfnhess(3,3,i))
    rhotens3(1,1,1)=rhotens3(1,1,1)+MOocc(i)*( 3*wfnderv(1,i)*wfnhess(1,1,i)+wfnval(i)*wfntens3(1,1,1,i) )
    rhotens3(2,2,2)=rhotens3(2,2,2)+MOocc(i)*( 3*wfnderv(2,i)*wfnhess(2,2,i)+wfnval(i)*wfntens3(2,2,2,i) )
    rhotens3(3,3,3)=rhotens3(3,3,3)+MOocc(i)*( 3*wfnderv(3,i)*wfnhess(3,3,i)+wfnval(i)*wfntens3(3,3,3,i) )
    rhotens3(1,1,2)=rhotens3(1,1,2)+MOocc(i)*( 2*wfnderv(1,i)*wfnhess(1,2,i)+wfnderv(2,i)*wfnhess(1,1,i)+wfnval(i)*wfntens3(1,1,2,i) )
    rhotens3(1,1,3)=rhotens3(1,1,3)+MOocc(i)*( 2*wfnderv(1,i)*wfnhess(1,3,i)+wfnderv(3,i)*wfnhess(1,1,i)+wfnval(i)*wfntens3(1,1,3,i) )
    rhotens3(2,2,3)=rhotens3(2,2,3)+MOocc(i)*( 2*wfnderv(2,i)*wfnhess(2,3,i)+wfnderv(3,i)*wfnhess(2,2,i)+wfnval(i)*wfntens3(2,2,3,i) )
    rhotens3(1,2,2)=rhotens3(1,2,2)+MOocc(i)*( 2*wfnderv(2,i)*wfnhess(2,1,i)+wfnderv(1,i)*wfnhess(2,2,i)+wfnval(i)*wfntens3(2,2,1,i) ) !=2,2,1 exchange 1<->2 from (1,1,2) to derive this
    rhotens3(1,3,3)=rhotens3(1,3,3)+MOocc(i)*( 2*wfnderv(3,i)*wfnhess(3,1,i)+wfnderv(1,i)*wfnhess(3,3,i)+wfnval(i)*wfntens3(3,3,1,i) ) !2.758568947939382D-002
    rhotens3(2,3,3)=rhotens3(2,3,3)+MOocc(i)*( 2*wfnderv(3,i)*wfnhess(3,2,i)+wfnderv(2,i)*wfnhess(3,3,i)+wfnval(i)*wfntens3(3,3,2,i) )
!     write(*,*) "A",wfnhess(3,1,i),wfntens3(1,3,3,i)
end do
dxx=2D0*dxx
dyy=2D0*dyy
dzz=2D0*dzz
value=laplfac*(dxx+dyy+dzz)
rhotens3=rhotens3*2D0*laplfac
grad(1)=rhotens3(1,1,1)+rhotens3(1,2,2)+rhotens3(1,3,3)
grad(2)=rhotens3(1,1,2)+rhotens3(2,2,2)+rhotens3(2,3,3)
grad(3)=rhotens3(1,1,3)+rhotens3(2,2,3)+rhotens3(3,3,3)

! diff=1D-5
!Check of flapldx,dy,dz is correct!
! write(*,*) grad(:)
! difflapldx=(flapl(x+diff,y,z,'t')-flapl(x-diff,y,z,'t'))/(2*diff)
! difflapldy=(flapl(x,y+diff,z,'t')-flapl(x,y-diff,z,'t'))/(2*diff)
! difflapldz=(flapl(x,y,z+diff,'t')-flapl(x,y,z-diff,'t'))/(2*diff)
! write(*,*) difflapldx,difflapldy,difflapldz

!Check deviation between analytic and numerical solution
! write(*,*) rhotens3(1,1,1),rhotens3(1,2,2),rhotens3(1,3,3)
! diffrhodxxx=(flapl(x+diff,y,z,'x')-flapl(x-diff,y,z,'x'))/(2*diff)
! diffrhodxyy=(flapl(x+diff,y,z,'y')-flapl(x-diff,y,z,'y'))/(2*diff)
! diffrhodxzz=(flapl(x+diff,y,z,'z')-flapl(x-diff,y,z,'z'))/(2*diff)
! write(*,*) diffrhodxxx,diffrhodxyy,diffrhodxzz
! write(*,*)

!Check diagonal term with finite difference
! write(*,*) rhotens3(1,1,1),rhotens3(2,2,2),rhotens3(3,3,3)
! diffrhodxxx=(flapl(x+diff,y,z,'x')-flapl(x-diff,y,z,'x'))/(2*diff)
! diffrhodyyy=(flapl(x,y+diff,z,'y')-flapl(x,y-diff,z,'y'))/(2*diff)
! diffrhodzzz=(flapl(x,y,z+diff,'z')-flapl(x,y,z-diff,'z'))/(2*diff)
! write(*,*) diffrhodxxx,diffrhodyyy,diffrhodzzz

if (allocated(b_EDF)) then
    call EDFrho(5,x,y,z,EDFdens,EDFgrad,EDFhess,EDFtens3)
    grad(1)=grad(1)+EDFtens3(1,1,1)+EDFtens3(1,2,2)+EDFtens3(1,3,3)
    grad(2)=grad(2)+EDFtens3(1,1,2)+EDFtens3(2,2,2)+EDFtens3(2,3,3)
    grad(3)=grad(3)+EDFtens3(1,1,3)+EDFtens3(2,2,3)+EDFtens3(3,3,3)
end if

!Don't consider Laplacian currently
if (itype==2) then !Calculate Hessian of laplacian
end if
end subroutine


!!!------------------------- Output Lagrangian kinetic G(r) at a point. idir=0/1/2/3 means total/x/y/z kinetic energy
real*8 function Lagkin(x,y,z,idir)
real*8 x,y,z,wfnval(nmo),wfnderv(3,nmo)
integer idir
call orbderv(2,1,nmo,x,y,z,wfnval,wfnderv)
lagkin=0D0
if (idir==0) then
    do imo=1,nmo
        lagkin=lagkin+MOocc(imo)*sum(wfnderv(:,imo)**2)
    end do
else
    do imo=1,nmo
        lagkin=lagkin+MOocc(imo)*wfnderv(idir,imo)**2
    end do
end if
lagkin=lagkin/2D0
end function


!!------------- Output Hamiltonian kinetic K(r) at a point. idir=0/1/2/3 means total/X/Y/Z kinetic energy
real*8 function Hamkin(x,y,z,idir)
real*8 x,y,z,wfnval(nmo),wfnderv(3,nmo),wfnhess(3,3,nmo)
integer idir
call orbderv(3,1,nmo,x,y,z,wfnval,wfnderv,wfnhess)
if (idir==0) then
    hamx=sum( MOocc(1:nmo)*wfnhess(1,1,1:nmo)*wfnval(1:nmo) )
    hamy=sum( MOocc(1:nmo)*wfnhess(2,2,1:nmo)*wfnval(1:nmo) )
    hamz=sum( MOocc(1:nmo)*wfnhess(3,3,1:nmo)*wfnval(1:nmo) )
    Hamkin=hamx+hamy+hamz
else if (idir==1) then
    Hamkin=sum( MOocc(1:nmo)*wfnhess(1,1,1:nmo)*wfnval(1:nmo) )
else if (idir==2) then
    Hamkin=sum( MOocc(1:nmo)*wfnhess(2,2,1:nmo)*wfnval(1:nmo) )
else if (idir==3) then
    Hamkin=sum( MOocc(1:nmo)*wfnhess(3,3,1:nmo)*wfnval(1:nmo) )
end if
Hamkin=-Hamkin/2D0
end function


!!!--------- Output electron linear momentum density in 3D representation at a point. idir=0/1/2/3 means magnitude/x/y/z component
real*8 function elemomdens(x,y,z,idir)
real*8 x,y,z,wfnval(nmo),wfnderv(3,nmo),comp(0:3)
integer idir
call orbderv(2,1,nmo,x,y,z,wfnval,wfnderv)
comp=0
do imo=1,nmo
    comp(1:3)=comp(1:3)-MOocc(imo)*wfnval(imo)*wfnderv(1:3,imo)
end do
comp(0)=sum(comp(1:3)**2)
elemomdens=comp(idir)
end function


!!!--------- Output magnetic dipole moment density at a point. idir=0/1/2/3 means magnitude/x/y/z component
real*8 function magmomdens(x,y,z,idir)
real*8 x,y,z,wfnval(nmo),wfnderv(3,nmo),comp(0:3)
integer idir
call orbderv(2,1,nmo,x,y,z,wfnval,wfnderv)
comp=0
do imo=1,nmo
    tmpx=-wfnval(imo)*(y*wfnderv(3,imo)-z*wfnderv(2,imo))
    tmpy=-wfnval(imo)*(z*wfnderv(1,imo)-x*wfnderv(3,imo))
    tmpz=-wfnval(imo)*(x*wfnderv(2,imo)-y*wfnderv(1,imo))
    comp(1)=comp(1)+MOocc(imo)*tmpx
    comp(2)=comp(2)+MOocc(imo)*tmpy
    comp(3)=comp(3)+MOocc(imo)*tmpz
end do
comp(0)=sum(comp(1:3)**2)
magmomdens=comp(idir)
end function


!!!----- Calculate Sign(lambda2(r))*rho(r) function, this is a warpper used to convert subroutine to function form
real*8 function signlambda2rho(x,y,z)
real*8 x,y,z,sl2r,RDG,rho
call signlambda2rho_RDG(x,y,z,sl2r,RDG,rho)
signlambda2rho=sl2r
end function
!!!------ Calculate signlambda2rho and RDG at the same time
subroutine signlambda2rho_RDG(x,y,z,sl2r,RDG,elerho)
use util
real*8 x,y,z,elerho,sl2r,RDG
real*8 eigvecmat(3,3),eigval(3),elehess(3,3),elegrad(3) !Hessian of electron density
call calchessmat_dens(2,x,y,z,elerho,elegrad,elehess)
call diagmat(elehess,eigvecmat,eigval,100,1D-10)
call sort(eigval)
if (eigval(2)==0D0) then !When eigval(2)==0.0D0, eigval(2)/abs(eigval(2)) can't be calculated, elerho generally will be zero, so sign is not important
    sl2r=elerho
else
    sl2r=elerho*eigval(2)/abs(eigval(2))
end if
sumgrad2=sum(elegrad(:)**2)
if (RDG_maxrho/=0D0.and.elerho>=RDG_maxrho) then
    RDG=100D0
!This occurs at distant region when exponent cutoff is used, the actual value should be very large. In order to avoid denominator become zero, we set it artifically to a big value
else if (sumgrad2==0D0.or.elerho==0D0) then
    RDG=999D0
else
    RDG=0.161620459673995D0*dsqrt(sumgrad2)/elerho**(4D0/3D0) !0.161620459673995D0=1/(2*(3*pi**2)**(1/3))
end if
end subroutine



!!!-----Output ELF(Electron Localization Function) or LOL(Localized orbital locator) and similar function value at a point
!label="ELF" or "LOL"
real*8 function ELF_LOL(x,y,z,label)
real*8 x,y,z,wfnval(nmo),wfnderv(3,nmo)
real*8 D,Dh,gradrho(3),gradrhoa(3),gradrhob(3),rho,rhoa,rhob,rhospin,MOoccnow
real*8 :: Fc=2.871234000D0 ! Fermi constant = (3/10)*(3*Pi^2)**(2/3) = 2.871234, 1/2.871234=0.34828
real*8 :: Fc_pol=4.557799872D0 ! Fermi constant for spin polarized = (3/10)*(6*Pi^2)**(2/3) = 4.5578, 1/4.5578=0.2194
character*3 label

!The ELF defined by Tsirelson, CPL,351,142.
!The LOL defined by Tsirelson, Acta Cryst. (2002). B58, 780-785.
!These functions are not important, so I don't write a special code
!for it, since rho, nebla-rho, nebla^2-rho support EDF, this function also support EDF
if (ELFLOL_type==1) then
    rho=fdens(x,y,z)
    if (wfntype==0.or.wfntype==3) then !close shell cases
        Dh=Fc*rho**(5.0D0/3.0D0)
    else if (wfntype==1.or.wfntype==2.or.wfntype==4) then !Open shell cases
        rhospin=fspindens(x,y,z,'s') !rhospin=rhoa-rhob, rho=rhoa+rhob
        rhoa=(rhospin+rho)/2D0
        rhob=(rho-rhospin)/2D0
        Dh=Fc_pol*(rhoa**(5.0D0/3.0D0)+rhob**(5.0D0/3.0D0)) !kinetic energy of HEG
    end if
    if (label=="ELF") then
        !Restrictly speaking, the kinetic energy expansion should be replace by polarized form for open-shell
        D=Dh-(1/9D0)*fgrad(x,y,z,'t')**2/rho+(1/6D0)*flapl(x,y,z,'t')
        ELF_LOL=1/(1+(D/Dh)**2)
    else if (label=="LOL") then
        D=Dh+(1/72D0)*fgrad(x,y,z,'t')**2/rho+(1/6D0)*flapl(x,y,z,'t')
        t=Dh/D
        ELF_LOL=1D0/(1D0/t+1)
    end if
    if (ELF_LOL<ELFLOL_cut) ELF_LOL=0
    return
end if

!For Becke's and Lu Tian's definition below
call orbderv(2,1,nmo,x,y,z,wfnval,wfnderv)

D=0.0D0
rho=0.0D0
rhoa=0D0
rhob=0D0
gradrho=0D0
gradrhoa=0D0
gradrhob=0D0
if (label=="ELF") then
    if (wfntype==0.or.wfntype==3) then !spin-unpolarized case
        do i=1,nmo
            rho=rho+MOocc(i)*wfnval(i)**2
            gradrho(:)=gradrho(:)+2.0D0*MOocc(i)*wfnval(i)*wfnderv(:,i)
            D=D+MOocc(i)*(sum(wfnderv(:,i)**2)) !Calculate actual kinetic term
        end do        
        Dh=Fc*rho**(5.0D0/3.0D0) !Thomas-Fermi uniform electron gas kinetic energy
        D=D/2D0
        if (rho/=0D0) D=D-(sum(gradrho(:)**2))/rho/8D0
    else if (wfntype==1.or.wfntype==2.or.wfntype==4) then !spin-polarized case
        do i=1,nmo
            MOoccnow=MOocc(i)
            if (MOtype(i)==0) MOoccnow=MOocc(i)/2D0 !double occupied, present when wfntype==2 (ROHF), alpha and beta get half part
            if (MOtype(i)==1.or.MOtype(i)==0) then
                rhoa=rhoa+MOoccnow*wfnval(i)**2
                gradrhoa(:)=gradrhoa(:)+2.0D0*MOoccnow*wfnval(i)*wfnderv(:,i)
            end if
            if (MOtype(i)==2.or.MOtype(i)==0) then
                rhob=rhob+MOoccnow*wfnval(i)**2
                gradrhob(:)=gradrhob(:)+2.0D0*MOoccnow*wfnval(i)*wfnderv(:,i)
            end if
            D=D+MOocc(i)*(sum(wfnderv(:,i)**2)) !Calculate actual kinetic term
        end do
        Dh=Fc_pol*(rhoa**(5.0D0/3.0D0)+rhob**(5.0D0/3.0D0))
        D=D/2D0
        if (rhoa/=0D0) D=D-(sum(gradrhoa(:)**2))/rhoa/8
        if (rhob/=0D0) D=D-(sum(gradrhob(:)**2))/rhob/8
    end if
    if (ELFLOL_type==0.or.ELFLOL_type==2) then
        if (ELF_addminimal==1) D=D+1D-5 !add 1D-5 to avoid D become zero, leads to unity in infinite
        if (ELFLOL_type==0) ELF_LOL=1/(1+(D/Dh)**2)
        if (ELFLOL_type==2) ELF_LOL=1/(1+(D/Dh)) !New formalism defined by Tian Lu
    end if
    if (ELFLOL_type==3) ELF_LOL=D/Dh
    if (ELFLOL_type==995) ELF_LOL=Dh !Thomas-Fermi kinetic energy density
    if (ELFLOL_type==996) ELF_LOL=D/Dh !For test
    if (ELFLOL_type==997) ELF_LOL=D !For test
    if (ELFLOL_type==998) ELF_LOL=1/(1+D) !For test
    if (ELFLOL_type==999) ELF_LOL=1/(1+D**2) !For test
!----------------------------
else if (label=="LOL") then
    t=0.0D0
    if (wfntype==0.or.wfntype==3) then !Spin-unpolarized case
        do i=1,nmo !Store actual kinetic energy to t first
            rho=rho+MOocc(i)*wfnval(i)**2
            t=t+MOocc(i)*(sum(wfnderv(:,i)**2))
        end do
        t=t/2D0
        Dh=Fc*rho**(5.0D0/3.0D0)
    else if (wfntype==1.or.wfntype==2.or.wfntype==4) then !spin-polarized case
        do i=1,nmo !Store actual kinetic energy to t first
            MOoccnow=MOocc(i)
            if (MOtype(i)==0) MOoccnow=MOocc(i)/2D0 !double occupied, present when wfntype==2 (ROHF), alpha and beta get half part
            if (MOtype(i)==1.or.MOtype(i)==0) rhoa=rhoa+MOoccnow*wfnval(i)**2
            if (MOtype(i)==2.or.MOtype(i)==0) rhob=rhob+MOoccnow*wfnval(i)**2
            t=t+MOocc(i)*(sum(wfnderv(:,i)**2))
        end do
        t=t/2D0
        Dh=Fc_pol*(rhoa**(5.0D0/3.0D0)+rhob**(5.0D0/3.0D0))
    end if
!--------- A new definition of LOL, however the value range is not as good as LOL
! ELF_LOL=t-Dh
! if (ELF_LOL>0) then
!     ELF_LOL=1D0/(1D0+1D0/ELF_LOL)
! else if (ELF_LOL<0) then
!     ELF_LOL=-1D0/(1D0+1D0/abs(ELF_LOL))
! end if
!-------------
!If there is very long distance between molecule and current point, t (above) is zero,
!and t=Dh/t is also zero (because rho converges faster), but can't be calculate directly, so simply skip
    if (t/=0.0D0) t=Dh/t
    if (ELFLOL_type==0) ELF_LOL=1D0/(1D0/t+1) !namely t/(1+t). This is default case
    if (ELFLOL_type==2) ELF_LOL=1D0/((1D0/t)**2+1) !New form defined by Tian Lu
end if
if (ELF_LOL<ELFLOL_cut) ELF_LOL=0
end function


!!!----- Calculate ELF/LOL, its gradient and Hessian matrix at x,y,z, store to hess 
!!!!!!!!!!!! currently can not calculate Hessian
!funsel="ELF" or "LOL"
!itype=1 Only calculate value and gradient, not Hessian
!itype=2 Calculate value, gradient and Hessian (Not available)
subroutine calchessmat_ELF_LOL(itype,x,y,z,value,grad,hess,funsel)
use util
integer itype
real*8 x,y,z,value,grad(3),hess(3,3),MOoccnow
real*8 wfnval(nmo),wfnderv(3,nmo),wfnhess(3,3,nmo),wfntens3(3,3,3,nmo)
real*8 rho,gradrho(3),hessrho(3,3),Dh,gradDh(3),hessDh(3,3),Ts,gradTs(3),hessTs(3,3),Wei,gradWei(3),hessWei(3,3)
real*8 rhoa,rhob,gradrhoa(3),gradrhob(3),hessrhoa(3,3),hessrhob(3,3),Dha,Dhb,gradDha(3),gradDhb(3),Weia,Weib,gradWeia(3),gradWeib(3)
! real*8 hessDha(3,3),hessDhb(3,3),hessWeia(3,3),hessWeib(3,3)
real*8 :: Fc=2.871234000D0 ! Fermi constant = (3/10)*(3*Pi^2)**(2/3) = 2.871234, 1/2.871234=0.34828
real*8 :: Fc_pol=4.557799872D0 ! Fermi constant for spin polarized = (3/10)*(6*Pi^2)**(2/3) = 4.5578, 1/4.5578=0.2194
real*8 :: corrELF=1D-5
character funsel*3

if (itype==1) call orbderv(4,1,nmo,x,y,z,wfnval,wfnderv,wfnhess) !Get Hessian of GTF, needn't 3-order tensor
if (itype==2) call orbderv(5,1,nmo,x,y,z,wfnval,wfnderv,wfnhess,wfntens3)

!spin-unpolarized case
if (wfntype==0.or.wfntype==3) then
    rho=0D0
    gradrho=0D0
    Ts=0D0
    gradTs=0D0
    do i=1,nmo
        rho=rho+MOocc(i)*wfnval(i)**2
        gradrho(:)=gradrho(:)+MOocc(i)*wfnval(i)*wfnderv(:,i)
        Ts=Ts+MOocc(i)*(sum(wfnderv(:,i)**2))
    end do
    gradrho=2*gradrho
    Ts=Ts/2D0
    Dh=Fc*rho**(5D0/3D0)
    gradDh(:)=5D0/3D0*Fc*rho**(2D0/3D0)*gradrho(:)
    do i=1,nmo
        gradTs(1)=gradTs(1)+MOocc(i)*(wfnderv(1,i)*wfnhess(1,1,i)+wfnderv(2,i)*wfnhess(1,2,i)+wfnderv(3,i)*wfnhess(1,3,i))
        gradTs(2)=gradTs(2)+MOocc(i)*(wfnderv(2,i)*wfnhess(2,2,i)+wfnderv(1,i)*wfnhess(2,1,i)+wfnderv(3,i)*wfnhess(2,3,i))
        gradTs(3)=gradTs(3)+MOocc(i)*(wfnderv(3,i)*wfnhess(3,3,i)+wfnderv(1,i)*wfnhess(3,1,i)+wfnderv(2,i)*wfnhess(3,2,i))
    end do
    if (funsel=="ELF") then
        !Calculate Hessian for rho
        hessrho(1,1)=2*sum( MOocc(1:nmo)*( wfnderv(1,1:nmo)**2 + wfnval(1:nmo)*wfnhess(1,1,1:nmo) ) )
        hessrho(2,2)=2*sum( MOocc(1:nmo)*( wfnderv(2,1:nmo)**2 + wfnval(1:nmo)*wfnhess(2,2,1:nmo) ) )
        hessrho(3,3)=2*sum( MOocc(1:nmo)*( wfnderv(3,1:nmo)**2 + wfnval(1:nmo)*wfnhess(3,3,1:nmo) ) )
        hessrho(1,2)=2*sum( MOocc(1:nmo)*( wfnderv(1,1:nmo)*wfnderv(2,1:nmo)+wfnhess(1,2,1:nmo)*wfnval(1:nmo) ) )
        hessrho(2,3)=2*sum( MOocc(1:nmo)*( wfnderv(2,1:nmo)*wfnderv(3,1:nmo)+wfnhess(2,3,1:nmo)*wfnval(1:nmo) ) )
        hessrho(1,3)=2*sum( MOocc(1:nmo)*( wfnderv(1,1:nmo)*wfnderv(3,1:nmo)+wfnhess(1,3,1:nmo)*wfnval(1:nmo) ) )
        hessrho(2,1)=hessrho(1,2)
        hessrho(3,1)=hessrho(1,3)
        hessrho(3,2)=hessrho(2,3)
        !Calculate Weizsacker functional and its derivatives
        Wei=sum(gradrho(:)**2)/8D0/rho
        D=Ts-Wei+corrELF
        chi=D/Dh
        value=1D0/(1D0+chi**2)
        gradWei(1)= 0.25D0/rho*( gradrho(1)*hessrho(1,1)+gradrho(2)*hessrho(1,2)+gradrho(3)*hessrho(1,3) ) - wei/rho*gradrho(1)
        gradWei(2)= 0.25D0/rho*( gradrho(2)*hessrho(2,2)+gradrho(1)*hessrho(2,1)+gradrho(3)*hessrho(2,3) ) - wei/rho*gradrho(2)
        gradWei(3)= 0.25D0/rho*( gradrho(3)*hessrho(3,3)+gradrho(2)*hessrho(3,2)+gradrho(1)*hessrho(3,1) ) - wei/rho*gradrho(3)
        chidx=(gradTs(1)-gradWei(1))/Dh - gradDh(1)/Dh**2 *D
        chidy=(gradTs(2)-gradWei(2))/Dh - gradDh(2)/Dh**2 *D
        chidz=(gradTs(3)-gradWei(3))/Dh - gradDh(3)/Dh**2 *D
        grad(1)=-2D0*chi/(1+chi**2)**2 * chidx
        grad(2)=-2D0*chi/(1+chi**2)**2 * chidy
        grad(3)=-2D0*chi/(1+chi**2)**2 * chidz
    else if (funsel=="LOL") then
        value=1D0/(1D0+Ts/Dh)
        tmp=-1D0/Dh/(1D0+Ts/Dh)**2
        grad(1)=tmp*(gradTs(1)-Ts/Dh*gradDh(1))
        grad(2)=tmp*(gradTs(2)-Ts/Dh*gradDh(2))
        grad(3)=tmp*(gradTs(3)-Ts/Dh*gradDh(3))
    end if
!spin-polarized case
else if (wfntype==1.or.wfntype==2.or.wfntype==4) then
    rhoa=0D0
    rhob=0D0
    gradrhoa=0D0
    gradrhob=0D0
    Ts=0D0
    gradTs=0D0
    do i=1,nmo
        MOoccnow=MOocc(i)
        if (MOtype(i)==0) MOoccnow=MOocc(i)/2D0
        if (MOtype(i)==1.or.MOtype(i)==0) then
            rhoa=rhoa+MOoccnow*wfnval(i)**2
            gradrhoa(:)=gradrhoa(:)+MOoccnow*wfnval(i)*wfnderv(:,i)
        end if
        if (MOtype(i)==2.or.MOtype(i)==0) then
            rhob=rhob+MOoccnow*wfnval(i)**2
            gradrhob(:)=gradrhob(:)+MOoccnow*wfnval(i)*wfnderv(:,i)            
        end if
        Ts=Ts+MOocc(i)*(sum(wfnderv(:,i)**2))
    end do
    gradrhoa=2*gradrhoa
    gradrhob=2*gradrhob
    Ts=Ts/2D0
    Dha=Fc_pol*rhoa**(5D0/3D0)
    Dhb=Fc_pol*rhob**(5D0/3D0)
    Dh=Dha+Dhb
    gradDha(:)=5D0/3D0*Fc_pol*rhoa**(2D0/3D0)*gradrhoa(:)
    gradDhb(:)=5D0/3D0*Fc_pol*rhob**(2D0/3D0)*gradrhob(:)
    gradDh=gradDha+gradDhb
    do i=1,nmo
        gradTs(1)=gradTs(1)+MOocc(i)*(wfnderv(1,i)*wfnhess(1,1,i)+wfnderv(2,i)*wfnhess(1,2,i)+wfnderv(3,i)*wfnhess(1,3,i))
        gradTs(2)=gradTs(2)+MOocc(i)*(wfnderv(2,i)*wfnhess(2,2,i)+wfnderv(1,i)*wfnhess(2,1,i)+wfnderv(3,i)*wfnhess(2,3,i))
        gradTs(3)=gradTs(3)+MOocc(i)*(wfnderv(3,i)*wfnhess(3,3,i)+wfnderv(1,i)*wfnhess(3,1,i)+wfnderv(2,i)*wfnhess(3,2,i))
    end do
    
    if (funsel=="ELF") then
!         !Calculate Hessian for rho
        hessrhoa=0D0
        hessrhob=0D0
        do i=1,nmo
            MOoccnow=MOocc(i)
            if (MOtype(i)==0) MOoccnow=MOocc(i)/2D0 !double occupied, present when wfntype==2 (ROHF), alpha and beta get half part
            if (MOtype(i)==1.or.MOtype(i)==0) then
                hessrhoa(1,1)=hessrhoa(1,1)+MOoccnow*( wfnderv(1,i)**2 + wfnval(i)*wfnhess(1,1,i) )
                hessrhoa(2,2)=hessrhoa(2,2)+MOoccnow*( wfnderv(2,i)**2 + wfnval(i)*wfnhess(2,2,i) )
                hessrhoa(3,3)=hessrhoa(3,3)+MOoccnow*( wfnderv(3,i)**2 + wfnval(i)*wfnhess(3,3,i) )
                hessrhoa(1,2)=hessrhoa(1,2)+MOoccnow*( wfnderv(1,i)*wfnderv(2,i)+wfnhess(1,2,i)*wfnval(i) )
                hessrhoa(2,3)=hessrhoa(2,3)+MOoccnow*( wfnderv(2,i)*wfnderv(3,i)+wfnhess(2,3,i)*wfnval(i) )
                hessrhoa(1,3)=hessrhoa(1,3)+MOoccnow*( wfnderv(1,i)*wfnderv(3,i)+wfnhess(1,3,i)*wfnval(i) )
            end if
            if (MOtype(i)==2.or.MOtype(i)==0) then
                hessrhob(1,1)=hessrhob(1,1)+MOoccnow*( wfnderv(1,i)**2 + wfnval(i)*wfnhess(1,1,i) )
                hessrhob(2,2)=hessrhob(2,2)+MOoccnow*( wfnderv(2,i)**2 + wfnval(i)*wfnhess(2,2,i) )
                hessrhob(3,3)=hessrhob(3,3)+MOoccnow*( wfnderv(3,i)**2 + wfnval(i)*wfnhess(3,3,i) )
                hessrhob(1,2)=hessrhob(1,2)+MOoccnow*( wfnderv(1,i)*wfnderv(2,i)+wfnhess(1,2,i)*wfnval(i) )
                hessrhob(2,3)=hessrhob(2,3)+MOoccnow*( wfnderv(2,i)*wfnderv(3,i)+wfnhess(2,3,i)*wfnval(i) )
                hessrhob(1,3)=hessrhob(1,3)+MOoccnow*( wfnderv(1,i)*wfnderv(3,i)+wfnhess(1,3,i)*wfnval(i) )
            end if
        end do
        hessrhoa=hessrhoa*2
        hessrhob=hessrhob*2
        hessrhoa(2,1)=hessrhoa(1,2)
        hessrhoa(3,1)=hessrhoa(1,3)
        hessrhoa(3,2)=hessrhoa(2,3)
        hessrhob(2,1)=hessrhob(1,2)
        hessrhob(3,1)=hessrhob(1,3)
        hessrhob(3,2)=hessrhob(2,3)
!         !Calculate Weizsacker functional and its derivatives
        Weia=sum(gradrhoa(:)**2)/8D0/rhoa
        Weib=sum(gradrhob(:)**2)/8D0/rhob
        Wei=Weia+Weib
        D=Ts-Wei+corrELF
        chi=D/Dh
        value=1D0/(1D0+chi**2)
        gradWeia(1)= 0.25D0/rhoa*( gradrhoa(1)*hessrhoa(1,1)+gradrhoa(2)*hessrhoa(1,2)+gradrhoa(3)*hessrhoa(1,3) ) - weia/rhoa*gradrhoa(1)
        gradWeia(2)= 0.25D0/rhoa*( gradrhoa(2)*hessrhoa(2,2)+gradrhoa(1)*hessrhoa(2,1)+gradrhoa(3)*hessrhoa(2,3) ) - weia/rhoa*gradrhoa(2)
        gradWeia(3)= 0.25D0/rhoa*( gradrhoa(3)*hessrhoa(3,3)+gradrhoa(2)*hessrhoa(3,2)+gradrhoa(1)*hessrhoa(3,1) ) - weia/rhoa*gradrhoa(3)
        gradWeib(1)= 0.25D0/rhob*( gradrhob(1)*hessrhob(1,1)+gradrhob(2)*hessrhob(1,2)+gradrhob(3)*hessrhob(1,3) ) - weib/rhob*gradrhob(1)
        gradWeib(2)= 0.25D0/rhob*( gradrhob(2)*hessrhob(2,2)+gradrhob(1)*hessrhob(2,1)+gradrhob(3)*hessrhob(2,3) ) - weib/rhob*gradrhob(2)
        gradWeib(3)= 0.25D0/rhob*( gradrhob(3)*hessrhob(3,3)+gradrhob(2)*hessrhob(3,2)+gradrhob(1)*hessrhob(3,1) ) - weib/rhob*gradrhob(3)
        gradWei=gradWeia+gradWeib
        chidx=(gradTs(1)-gradWei(1))/Dh - gradDh(1)/Dh**2 *D
        chidy=(gradTs(2)-gradWei(2))/Dh - gradDh(2)/Dh**2 *D
        chidz=(gradTs(3)-gradWei(3))/Dh - gradDh(3)/Dh**2 *D
        grad(1)=-2D0*chi/(1+chi**2)**2 * chidx
        grad(2)=-2D0*chi/(1+chi**2)**2 * chidy
        grad(3)=-2D0*chi/(1+chi**2)**2 * chidz
    else if (funsel=="LOL") then
        value=1D0/(1D0+Ts/Dh)
        tmp=-1D0/Dh/(1D0+Ts/Dh)**2
        grad(1)=tmp*(gradTs(1)-Ts/Dh*gradDh(1))
        grad(2)=tmp*(gradTs(2)-Ts/Dh*gradDh(2))
        grad(3)=tmp*(gradTs(3)-Ts/Dh*gradDh(3))
    end if
end if

! if (itype==1) return
! ! Calculate Hessian for LOL, also need Hessian for rho
! hessrho(1,1)=2*sum( MOocc(1:nmo)*( wfnderv(1,1:nmo)**2 + wfnval(1:nmo)*wfnhess(1,1,1:nmo) ) )
! hessrho(2,2)=2*sum( MOocc(1:nmo)*( wfnderv(2,1:nmo)**2 + wfnval(1:nmo)*wfnhess(2,2,1:nmo) ) )
! hessrho(3,3)=2*sum( MOocc(1:nmo)*( wfnderv(3,1:nmo)**2 + wfnval(1:nmo)*wfnhess(3,3,1:nmo) ) )
! hessrho(1,2)=2*sum( MOocc(1:nmo)*( wfnderv(1,1:nmo)*wfnderv(2,1:nmo)+wfnhess(1,2,1:nmo)*wfnval(1:nmo) ) )
! hessrho(2,3)=2*sum( MOocc(1:nmo)*( wfnderv(2,1:nmo)*wfnderv(3,1:nmo)+wfnhess(2,3,1:nmo)*wfnval(1:nmo) ) )
! hessrho(1,3)=2*sum( MOocc(1:nmo)*( wfnderv(1,1:nmo)*wfnderv(3,1:nmo)+wfnhess(1,3,1:nmo)*wfnval(1:nmo) ) )
! hessrho(2,1)=hessrho(1,2)
! hessrho(3,1)=hessrho(1,3)
! hessrho(3,2)=hessrho(2,3)
! 
! hessTs=0D0
! do i=1,nmo
!     hessTs(1,1)=hessTs(1,1)+MOocc(i)*( wfnhess(1,1,i)**2 + wfnderv(1,i)*wfntens3(1,1,1,i) + wfnhess(1,2,i)**2 + wfnderv(2,i)*wfntens3(1,1,2,i) + wfnhess(1,3,i)**2 + wfnderv(3,i)*wfntens3(1,1,3,i) )
!     hessTs(2,2)=hessTs(2,2)+MOocc(i)*( wfnhess(2,2,i)**2 + wfnderv(2,i)*wfntens3(2,2,2,i) + wfnhess(2,1,i)**2 + wfnderv(1,i)*wfntens3(2,2,1,i) + wfnhess(2,3,i)**2 + wfnderv(3,i)*wfntens3(2,2,3,i) )
!     hessTs(3,3)=hessTs(3,3)+MOocc(i)*( wfnhess(3,3,i)**2 + wfnderv(3,i)*wfntens3(3,3,3,i) + wfnhess(3,2,i)**2 + wfnderv(2,i)*wfntens3(3,3,2,i) + wfnhess(3,1,i)**2 + wfnderv(1,i)*wfntens3(3,3,1,i) )
!     hessTs(1,2)=hessTs(1,2)+MOocc(i)*( wfnhess(1,2,i)*wfnhess(1,1,i) + wfnderv(1,i)*wfntens3(1,1,2,i) + wfnhess(2,2,i)*wfnhess(1,2,i) + wfnderv(2,i)*wfntens3(1,2,2,i) + wfnhess(2,3,i)*wfnhess(1,3,i) + wfnderv(3,i)*wfntens3(1,2,3,i) )
!     hessTs(2,3)=hessTs(2,3)+MOocc(i)*( wfnhess(2,3,i)*wfnhess(2,2,i) + wfnderv(2,i)*wfntens3(2,2,3,i) + wfnhess(3,3,i)*wfnhess(2,3,i) + wfnderv(3,i)*wfntens3(2,3,3,i) + wfnhess(3,1,i)*wfnhess(2,1,i) + wfnderv(1,i)*wfntens3(2,3,1,i) )
!     hessTs(1,3)=hessTs(1,3)+MOocc(i)*( wfnhess(1,3,i)*wfnhess(1,1,i) + wfnderv(1,i)*wfntens3(1,1,3,i) + wfnhess(3,3,i)*wfnhess(1,3,i) + wfnderv(3,i)*wfntens3(1,3,3,i) + wfnhess(3,2,i)*wfnhess(1,2,i) + wfnderv(2,i)*wfntens3(1,3,2,i) )
! end do
! hessTs(2,1)=hessTs(1,2)
! hessTs(3,1)=hessTs(1,3)
! hessTs(3,2)=hessTs(2,3)
! 
! tmp1=10D0/9D0*Fc/rho**(1D0/3D0)
! tmp2=5D0/3D0*Fc*rho**(2D0/3D0)
! hessDh(1,1)=tmp1*gradrho(1)**2 + tmp2*hessrho(1,1)
! hessDh(2,2)=tmp1*gradrho(2)**2 + tmp2*hessrho(2,2)
! hessDh(3,3)=tmp1*gradrho(3)**2 + tmp2*hessrho(3,3)
! hessDh(1,2)=tmp1*gradrho(1)*gradrho(2) + tmp2*hessrho(1,2)
! hessDh(2,3)=tmp1*gradrho(2)*gradrho(3) + tmp2*hessrho(2,3)
! hessDh(1,3)=tmp1*gradrho(1)*gradrho(3) + tmp2*hessrho(1,3)
! hessDh(2,1)=hessDh(1,2)
! hessDh(3,1)=hessDh(1,3)
! hessDh(3,2)=hessDh(2,3)
! 
! !Diagonal of Hessian of LOL
! apre=1/Dh**2/(1+Ts/Dh)**3
! bpre=-1/Dh/(1+Ts/Dh)**2
! apartxx=(gradDh(1)+2*gradTs(1)-Ts/Dh*gradDh(1)) * (gradTs(1)-Ts/Dh*gradDh(1))
! bpartxx=hessTs(1,1)-( gradDh(1)*gradTs(1)-Ts/Dh*gradDh(1)**2+Ts*hessDh(1,1) )/Dh
! hess(1,1)=apre*apartxx+bpre*bpartxx
! apartyy=(gradDh(2)+2*gradTs(2)-Ts/Dh*gradDh(2)) * (gradTs(2)-Ts/Dh*gradDh(2))
! bpartyy=hessTs(2,2)-( gradDh(2)*gradTs(2)-Ts/Dh*gradDh(2)**2+Ts*hessDh(2,2) )/Dh
! hess(2,2)=apre*apartyy+bpre*bpartyy
! apartzz=(gradDh(3)+2*gradTs(3)-Ts/Dh*gradDh(3)) * (gradTs(3)-Ts/Dh*gradDh(3))
! bpartzz=hessTs(3,3)-( gradDh(3)*gradTs(3)-Ts/Dh*gradDh(3)**2+Ts*hessDh(3,3) )/Dh
! hess(3,3)=apre*apartzz+bpre*bpartzz
! !Non-diagonal of Hessian of LOL
! bpartxy=hessTs(1,2)-( gradDh(1)*gradTs(2)-Ts/Dh*gradDh(1)*gradDh(2)+Ts*hessDh(1,2) )/Dh
! hess(1,2)=apre*apartxx+bpre*bpartxy
! bpartyz=hessTs(2,3)-( gradDh(2)*gradTs(3)-Ts/Dh*gradDh(2)*gradDh(3)+Ts*hessDh(2,3) )/Dh
! hess(2,3)=apre*apartyy+bpre*bpartyz
! bpartxz=hessTs(1,3)-( gradDh(1)*gradTs(3)-Ts/Dh*gradDh(1)*gradDh(3)+Ts*hessDh(1,3) )/Dh
! hess(1,3)=apre*apartxx+bpre*bpartxz
! hess(2,1)=hess(1,2)
! hess(3,1)=hess(1,3)
! hess(3,2)=hess(2,3)
end subroutine


!-------- Output average local ionization energy at a point
real*8 function avglocion(x,y,z)
real*8 x,y,z,wfnval(nmo)
call orbderv(1,1,nmo,x,y,z,wfnval)
avglocion=0D0
rho=0D0
do i=1,nmo
    avglocion=avglocion+abs(MOene(i))*MOocc(i)*wfnval(i)**2
    rho=rho+MOocc(i)*wfnval(i)**2 !Calculate rho
end do
if (rho==0D0) then
    avglocion=0D0 !Avoid at distant region rho becomes zero when exponent cutoff is used
else
    avglocion=avglocion/rho
end if
end function

!-------- Calculate average local ionization energy at a point and meantime decompose it to occupied orbitals contribution
subroutine avglociondecomp(ifileid,x,y,z)
real*8 x,y,z,wfnval(nmo)
integer ifileid
character orbtype*2
call orbderv(1,1,nmo,x,y,z,wfnval)
totALIE=0D0
rho=0D0
do i=1,nmo
    totALIE=totALIE+abs(MOene(i))*MOocc(i)*wfnval(i)**2
    rho=rho+MOocc(i)*wfnval(i)**2 !Calculate rho
end do
if (rho==0D0) then
    totALIE=0D0 !Avoid at distant region rho becomes zero when exponent cutoff is used
else
    totALIE=totALIE/rho
end if
write(ifileid,"(' Average local ionization energy:',f16.10,' a.u.')") totALIE
write(ifileid,*) "Contribution of each orbital to average local ionization energy (a.u.):"
do i=1,nmo
    if (MOtype(i)==0) orbtype="AB"
    if (MOtype(i)==1) orbtype="A "
    if (MOtype(i)==2) orbtype="B "
    write(ifileid,"(' Orbital',i6,'  Ene:',f12.6,'  Occ:',f5.2,'  Type:',a,'  Contribution:',f12.6)") i,MOene(i),MOocc(i),orbtype,abs(MOene(i))*MOocc(i)*wfnval(i)**2/rho
end do
end subroutine


!-------- Output local electron affinity at a point
!Since virtual orbitals are involved, such as .fch/.molden/.gms must be used
real*8 function loceleaff(x,y,z)
real*8 x,y,z,wfnval(nmo)
call orbderv(1,1,nmo,x,y,z,wfnval)
loceleaff=0D0
rho=0D0
do i=1,nmo
    if (MOocc(i)==0) then !Only cycles unoccupied orbitals 
        loceleaff=loceleaff-MOene(i)*wfnval(i)**2 !Don't need to multiply "occupation number", because ROHF is not allowed, so all orbitals have the same type
        rho=rho+wfnval(i)**2 !Calculate rho
    end if
end do
if (rho==0D0) then
    loceleaff=0D0 !Avoid at distant region rho become zero when exponent cutoff is used
else
    loceleaff=loceleaff/rho
end if
end function

!!!------ Approximate form of DFT linear response kernel for close-shell, X(r1,r2), see Eq.3 of PCCP,14,3960. Only applied to the case wfntype=0
! r1 is taken as reference point and determined by refx,refy,refz. x,y,z in the argument is the coordinate of r2
real*8 function linrespkernel(x,y,z)
real*8 x,y,z,orbvalr1(nmo),orbvalr2(nmo)
call orbderv(1,1,nmo,refx,refy,refz,orbvalr1)
call orbderv(1,1,nmo,x,y,z,orbvalr2)
linrespkernel=0D0
do imo=1,nmo !Cycle occupied MOs
    if (nint(MOocc(imo))==2D0) then
        do jmo=idxHOMO+1,nmo !Cycle unoccupied MOs
            if (nint(MOocc(jmo))==0D0) linrespkernel=linrespkernel+orbvalr1(imo)*orbvalr1(jmo)*orbvalr2(jmo)*orbvalr2(imo)/(MOene(imo)-MOene(jmo))
        end do
    end if
end do
linrespkernel=linrespkernel*4
end function

!!!------------- Output Exchange-correlation density, correlation hole and correlation factor
!rfx,rfy,rfz is reference point (commonly use refx,refy,refz in defvar module), namely r1
!x,y,z in the argument is the coordinate of r2
!Calculate which function is controlled by "pairfunctype" in settings.ini, correlation type is determined by "paircorrtype" in settings.ini
real*8 function pairfunc(rfx,rfy,rfz,x,y,z)
real*8 rfx,rfy,rfz,x,y,z,orbvalr1(nmo),orbvalr2(nmo)
call orbderv(1,1,nmo,rfx,rfy,rfz,orbvalr1)
call orbderv(1,1,nmo,x,y,z,orbvalr2)
!Calculate alpha and beta density at r1 and r2
adensr1=0.0D0
bdensr1=0.0D0
adensr2=0.0D0
bdensr2=0.0D0
do i=1,nmo
    if (MOtype(i)==0) then
        adensr1=adensr1+MOocc(i)/2*orbvalr1(i)**2
        adensr2=adensr2+MOocc(i)/2*orbvalr2(i)**2
        bdensr1=bdensr1+MOocc(i)/2*orbvalr1(i)**2
        bdensr2=bdensr2+MOocc(i)/2*orbvalr2(i)**2
    else if (MOtype(i)==1) then
        adensr1=adensr1+MOocc(i)*orbvalr1(i)**2
        adensr2=adensr2+MOocc(i)*orbvalr2(i)**2
    else if (MOtype(i)==2) then
        bdensr1=bdensr1+MOocc(i)*orbvalr1(i)**2
        bdensr2=bdensr2+MOocc(i)*orbvalr2(i)**2
    end if
end do
totdensr1=adensr1+bdensr1
totdensr2=adensr2+bdensr2

ntime=1
if (pairfunctype==12) ntime=2 !Will need both alpha and beta information (aXCdens and bXCdens), so process twice
do itime=1,ntime
    !Calculate exchange-correlation density first, and then calculate correlation hole and correlation factor
    !For RHF/ROHF, we calculate them as if present system is open-shell
    if (pairfunctype==1.or.pairfunctype==4.or.pairfunctype==7.or.pairfunctype==10.or.(pairfunctype==12.and.itime==1)) then !Set start and end index of alpha orbitals
        !Cycle alpha orbitals first to obtain aXCdens
        istart=1
        if (wfntype==0.or.wfntype==2.or.wfntype==3) then !RHF,ROHF,R-post-HF
            iend=nmo
        else if (wfntype==1.or.wfntype==4) then !UHF, U-post-HF
            do iend=nmo,1,-1
                if (MOtype(iend)==1) exit
            end do
        end if
    else if (pairfunctype==2.or.pairfunctype==5.or.pairfunctype==8.or.pairfunctype==11.or.(pairfunctype==12.and.itime==2)) then !Set start and end index of beta orbitals
        if (wfntype==0.or.wfntype==3) then !RHF,R-post-HF
            istart=1
            iend=nmo
        else if (wfntype==2) then !ROHF
            istart=1
            do iend=1,nmo
                if (MOtype(iend)==1) exit
            end do
            iend=iend-1
        else if (wfntype==1.or.wfntype==4) then !UHF, U-post-HF
            do istart=1,nmo
                if (MOtype(istart)==2) exit
            end do
            iend=nmo
            if (nint(nbelec)==0) iend=0 !less than istart, so below cycle will be skipped
        end if
    end if

    XCtmp=0D0 !Really X+C density
    Xtmp=0D0 !Only X density
    Ctmp=0D0 !Only C density
    do i=istart,iend
        occi=MOocc(i)
        if (MOtype(i)==0) occi=occi/2D0 !Split close-shell orbital to spin orbital
        do j=istart,iend
            occj=MOocc(j)
            if (MOtype(j)==0) occj=occj/2D0
            tmpmul=orbvalr1(i)*orbvalr2(j)*orbvalr1(j)*orbvalr2(i)
            XCtmp=XCtmp-dsqrt(occi*occj)*tmpmul
            Xtmp=Xtmp-occi*occj*tmpmul
            Ctmp=Ctmp+(occi*occj-dsqrt(occi*occj))*tmpmul
        end do
    end do
        
    if (pairfunctype==1.or.pairfunctype==4.or.pairfunctype==7.or.pairfunctype==10.or.(pairfunctype==12.and.itime==1)) then
        if (paircorrtype==1) aXCdens=Xtmp
        if (paircorrtype==2) aXCdens=Ctmp
        if (paircorrtype==3) aXCdens=XCtmp
        acorrhole=aXCdens/adensr1
        acorrfac=acorrhole/adensr2
        if (pairfunctype==1) pairfunc=acorrhole
        if (pairfunctype==4) pairfunc=acorrfac
        if (pairfunctype==7) pairfunc=aXCdens
        if (pairfunctype==10) pairfunc=adensr1*adensr2+aXCdens
    else if (pairfunctype==2.or.pairfunctype==5.or.pairfunctype==8.or.pairfunctype==11.or.(pairfunctype==12.and.itime==2)) then
        if (paircorrtype==1) bXCdens=Xtmp
        if (paircorrtype==2) bXCdens=Ctmp
        if (paircorrtype==3) bXCdens=XCtmp
        bcorrhole=bXCdens/bdensr1
        bcorrfac=bcorrhole/bdensr2
        if (pairfunctype==2) pairfunc=bcorrhole
        if (pairfunctype==5) pairfunc=bcorrfac
        if (pairfunctype==8) pairfunc=bXCdens
        if (pairfunctype==11) pairfunc=bdensr1*bdensr2+bXCdens
    end if
end do
if (pairfunctype==12) pairfunc=adensr1*(adensr2+bdensr2)+aXCdens +bdensr1*(adensr2+bdensr2)+bXCdens
end function



!!!------------- Output source function
real*8 function srcfunc(x,y,z,imode) !Default imode=1
real*8 x,y,z,denomin
integer imode
denomin=4*pi*dsqrt((x-refx)**2+(y-refy)**2+(z-refz)**2)
if (denomin==0D0) denomin=0.001D0
if (imode==1) srcfunc=-flapl(x,y,z,'t')/denomin !Used to study effect of laplacian everywhere on specific point
if (imode==2) srcfunc=-flapl(refx,refy,refz,'t')/denomin !Used to study effect of laplacian at specific point on everywhere
end function


!!!------------------------- Output RDG with promolecular approximation
real*8 function RDGprodens(x,y,z)
real*8 x,y,z,elerho,elegrad(3)
call calchessmat_prodens(x,y,z,elerho,elegrad)
elegradnorm=dsqrt(sum(elegrad**2))
if ((RDGprodens_maxrho/=0.0D0.and.elerho>=RDGprodens_maxrho)) then
    RDGprodens=100D0
else if (elegradnorm==0D0.or.elerho==0D0) then
    RDGprodens=999D0
else
    RDGprodens=0.161620459673995D0*elegradnorm/elerho**(4D0/3D0)
end if
end function
!!!----- Calculate Sign(lambda2(r))*rho(r) with promolecular approximation
!!! this is a shell used to convert subroutine to function form
real*8 function signlambda2rho_prodens(x,y,z)
real*8 x,y,z,sl2r,RDG
call signlambda2rho_RDG_prodens(x,y,z,sl2r,RDG)
signlambda2rho_prodens=sl2r
end function
!!!------ Calculate Sign(lambda2(r))*rho(r) and RDG at the same time with promolecular approximation
subroutine signlambda2rho_RDG_prodens(x,y,z,sl2r,RDG)
use util
real*8 x,y,z,elerho,RDG,sl2r,sumgrad2
real*8 eigvecmat(3,3),eigval(3),elehess(3,3),elegrad(3)
call calchessmat_prodens(x,y,z,elerho,elegrad,elehess)
call diagmat(elehess,eigvecmat,eigval,100,1D-6)
call sort(eigval)
if (eigval(2)/=0.0D0) then
    sl2r=elerho*eigval(2)/abs(eigval(2)) !At nuclei of single atom system, hessian returned may be zero matrix
else
    sl2r=-elerho !Around nuclei, eigval(2)/abs(eigval(2)) always be negative
end if
! sl2r=elerho !Only obtain promolecular density
sumgrad2=sum(elegrad(:)**2)
if ((RDGprodens_maxrho/=0.0D0.and.elerho>=RDGprodens_maxrho).or.elerho==0.0D0) then
    RDG=100D0
else if (sumgrad2==0D0.or.elerho==0D0) then
    RDG=999D0
else
    RDG=0.161620459673995D0*dsqrt(sumgrad2)/elerho**(4D0/3D0)
end if
end subroutine


!!!----- Calculate electron density, its gradient and Hessian matrix at x,y,z with promolecular approximation
!Notice that fragment must be properly defined!!!
subroutine calchessmat_prodens(xin,yin,zin,elerho,elegrad,elehess)
use util
real*8 elerho,xin,yin,zin
real*8,optional :: elegrad(3),elehess(3,3)
real*8 posarr(200),rhoarr(200)
elerho=0D0
derx=0D0
dery=0D0
derz=0D0
dxx=0D0
dyy=0D0
dzz=0D0
dxy=0D0
dyz=0D0
dxz=0D0
idohess=0
if (present(elehess)) idohess=1
do i=1,nfragatmnum
    iatm=fragatm(i)
    ind=a(iatm)%index
    rx=a(iatm)%x-xin !Relative x
    ry=a(iatm)%y-yin
    rz=a(iatm)%z-zin
    rx2=rx*rx
    ry2=ry*ry
    rz2=rz*rz
    r2=rx2+ry2+rz2
    if (atomdenscut==1) then !Tight cutoff, for CHNO corresponding to cutoff at rho=0.00001
        if (ind==1.and.r2>25D0) then !H, 6.63^2=43.9569. But this seems to be unnecessarily large, so I use 5^2=25
            cycle
        else if (ind==6.and.r2>58.6756D0) then !C, 7.66^2=58.6756
            cycle
        else if (ind==7.and.r2>43.917129D0) then !N, 6.627^2=43.917129
            cycle
        else if (ind==8.and.r2>34.9281D0) then !O, 5.91^2=34.9281
            cycle
        else if (r2>(2.5D0*vdwr(ind))**2) then !Other cases, larger than double of its vdw radius will be skipped
            cycle
        end if
    else if (atomdenscut==2) then !Medium cutoff, the result may be not as accurate as atomdenscut==1, but much more cheaper
        if (r2>(2.2D0*vdwr(ind))**2) cycle
    else if (atomdenscut==3) then !Loose cutoff, the most inaccurate
        if (r2>(1.8D0*vdwr(ind))**2) cycle
    else if (atomdenscut==4) then !Foolish cutoff, you need to know what you are doing
        if (r2>(1.5D0*vdwr(ind))**2) cycle
    end if
    r=dsqrt(r2)
    if (ind<=18) then !H~Ar
        r2_1d5=r2**1.5D0
        do j=1,3
            if (YWTatomcoeff(ind,j)==0D0) cycle
            expterm=YWTatomexp(ind,j)
            term=YWTatomcoeff(ind,j)*dexp(-r/expterm)
            elerho=elerho+term
            if (r==0D0) cycle !Derivative of STO at nuclei is pointless
            tmp=term/expterm/r
            derx=derx-tmp*rx !Calculating gradient doesn't cost detectable time, so always calculate it
            dery=dery-tmp*ry
            derz=derz-tmp*rz
            if (idohess==1) then
                tmp1=1/r2_1d5/expterm
                tmp2=1/r2/(expterm*expterm)
                dxx=dxx+term*(tmp1*rx2-1/r/expterm+tmp2*rx2)
                dyy=dyy+term*(tmp1*ry2-1/r/expterm+tmp2*ry2)
                dzz=dzz+term*(tmp1*rz2-1/r/expterm+tmp2*rz2)
                tmp=term*(tmp1+tmp2)
                dxy=dxy+rx*ry*tmp
                dyz=dyz+ry*rz*tmp
                dxz=dxz+rx*rz*tmp
            end if
        end do
    else !Heavier than Ar
!         if (atomdenscut>=1.and.r>3*vdwr(ind)) cycle !Be careful, so use hentai criterion
        call genatmraddens(ind,posarr,rhoarr,npt) !Extract spherically averaged radial density of corresponding element
        call lagintpol(posarr(1:npt),rhoarr(1:npt),npt,r,term,der1r,der2r,3)
        elerho=elerho+term
        der1rdr=der1r/r
        derx=derx+der1rdr*rx
        dery=dery+der1rdr*ry
        derz=derz+der1rdr*rz
        if (idohess==1) then !I don't know how below code works, but it really works. See promolecular_grid routine in props.f90 of NCIPlot
            tmpval=(der2r-der1rdr)/r2
            dxx=dxx+der1rdr+tmpval*rx2
            dyy=dyy+der1rdr+tmpval*ry2
            dzz=dzz+der1rdr+tmpval*rz2
            dxy=dxy+tmpval*rx*ry
            dyz=dyz+tmpval*ry*rz
            dxz=dxz+tmpval*rx*rz
        end if
    end if
end do
if (present(elegrad)) then
    elegrad(1)=derx
    elegrad(2)=dery
    elegrad(3)=derz
end if
if (idohess==1) then
    elehess(1,1)=dxx
    elehess(2,2)=dyy
    elehess(3,3)=dzz
    elehess(1,2)=dxy
    elehess(2,3)=dyz
    elehess(1,3)=dxz
    elehess(2,1)=dxy
    elehess(3,2)=dyz
    elehess(3,1)=dxz
end if
end subroutine


!!---- Calculate atomic density based on STO fitted or radial density
!if indSTO==0, then all atom densities will be evaluated based on interpolation. if indSTO=18, then use STO fitted atomic density for element <18
real*8 function calcatmdens(iatm,x,y,z,indSTO)
use util
real*8 rho,x,y,z,posarr(200),rhoarr(200)
integer iatm,indSTO
calcatmdens=0
r=dsqrt( (a(iatm)%x-x)**2 + (a(iatm)%y-y)**2 + (a(iatm)%z-z)**2 )
ind=a(iatm)%index
! if (r>6*vdwr(ind)) return !Doesn't improve speed evidently but deteriorate result in rare cases
if (ind<=indSTO) then !H~Ar, use STO fitted density. This is faster than using Lagrange interpolation technique, but not normalized to expected electron number
    do j=1,3
        if (YWTatomcoeff(ind,j)==0D0) cycle
        calcatmdens=calcatmdens+YWTatomcoeff(ind,j)*exp(-r/YWTatomexp(ind,j))
    end do
else
    call genatmraddens(ind,posarr,rhoarr,npt) !Extract spherically averaged radial density of corresponding element
    call lagintpol(posarr(1:npt),rhoarr(1:npt),npt,r,calcatmdens,der1r,der2r,1)
end if
end function
!!---- Calculate promolecular density purely based on interpolation of radial density, the promolecular density obtained in this manner is quite accurate
real*8 function calcprodens(x,y,z,indSTO)
real*8 x,y,z
integer indSTO
calcprodens=0
do iatm=1,ncenter
    calcprodens=calcprodens+calcatmdens(iatm,x,y,z,indSTO)
end do
end function



!!!--------------- Output Shannon information entropy function at a point
!itype=1 rho/N*ln(rho/N), this is normal definition
!itype=2 rho*ln(rho), this is Shannon information density, see J. Chem. Phys., 126, 191107
real*8 function infoentro(itype,x,y,z)
real*8 x,y,z,rho
integer itype
if (nelec==0D0) then
    infoentro=0D0
else
    rho=fdens(x,y,z)
    if (itype==1) rho=fdens(x,y,z)/nelec
    if (rho<=1D-100) then
        infoentro=0.0D0
    else
        infoentro=-rho*log(rho)
    end if
end if
end function


!!!------------------------- Output total ESP at a point
real*8 function totesp(x,y,z)
real*8 x,y,z
totesp=eleesp(x,y,z)+nucesp(x,y,z)
end function


!!!--------------- Output total ESP at a point, but skip the nucleus marked by variable "iskipnuc"
real*8 function totespskip(x,y,z,iskip)
real*8 x,y,z
integer iskip
totespskip=0
do i=1,ncenter
    if (i==iskip) cycle
    totespskip=totespskip+a(i)%charge/dsqrt((x-a(i)%x)**2+(y-a(i)%y)**2+(z-a(i)%z)**2)
end do
totespskip=totespskip+eleesp(x,y,z)
end function

!!!------------------------- Output ESP from nuclear or atomic charges at a point
!At nuclear positions, this function returns 1000 instead of infinity to avoid numerical problems
real*8 function nucesp(x,y,z)
nucesp=0D0
do i=1,nfragatmnum
    dist2mpx=(x-a(fragatm(i))%x)**2
    dist2mpy=(y-a(fragatm(i))%y)**2
    dist2mpz=(z-a(fragatm(i))%z)**2
    dist2=dist2mpx+dist2mpy+dist2mpz
    if (dist2==0D0) then
         nucesp=1D3
         return
    end if
    nucesp=nucesp+a(fragatm(i))%charge/dsqrt(dist2)
end do
end function


!------------------------- Calculate ESP from electrons at a point
real*8 function eleesp(Cx,Cy,Cz)
implicit none
integer,parameter :: narrmax=396 !Enough for h-type GTF, 5+5=10. Alri(0:10,0:5,0:5)-->11*6*6=396
real*8 Cx,Cy,Cz,term,ep,Ax,Ay,Az,Bx,By,Bz,Aexp,Bexp,Px,Py,Pz,prefac,tmpval
real*8 sqPC,sqAB,expngaPC,PAx,PAy,PAz,PBx,PBy,PBz,PCx,PCy,PCz,fjtmp,addesp
real*8 Alri(narrmax),Amsj(narrmax),Antk(narrmax),Fn(0:10) !Enough for h-type GTF, 5+5=10
real*8 twoepsqPC,tl,tm,tn,espprivate,espexpcut
integer nu,imo,iprim,jprim,maxFn,maplri(narrmax),mapmsj(narrmax),mapntk(narrmax),tmpnuml,tmpnumm,tmpnumn
integer Aix,Aiy,Aiz,Bix,Biy,Biz,l,r,i,m,s,j,n,t,k,icen,jcen,sumAi,sumBi
eleesp=0.0D0
espexpcut=log10(espprecutoff)*3
nthreads=getNThreads()
!$OMP parallel do private(Aix,Aiy,Aiz,Bix,Biy,Biz,l,r,i,m,s,j,n,t,k,icen,jcen,sumAi,sumBi, &
!$OMP nu,imo,iprim,jprim,maxFn,maplri,mapmsj,mapntk,tmpnuml,tmpnumm,tmpnumn,&
!$OMP twoepsqPC,tl,tm,tn,Alri,Amsj,Antk,Fn,&
!$OMP sqPC,sqAB,expngaPC,PAx,PAy,PAz,PBx,PBy,PBz,PCx,PCy,PCz,fjtmp,addesp,&
!$OMP term,ep,Ax,Ay,Az,Bx,By,Bz,Aexp,Bexp,Px,Py,Pz,prefac,tmpval,espprivate) shared(eleesp) schedule(dynamic) NUM_THREADS(nthreads)
do iprim=1,nprims
    espprivate=0D0
    icen=b(iprim)%center
    Aexp=b(iprim)%exp
    Ax=a(icen)%x
    Ay=a(icen)%y
    Az=a(icen)%z
    Aix=type2ix(b(iprim)%functype)
    Aiy=type2iy(b(iprim)%functype)
    Aiz=type2iz(b(iprim)%functype)
    sumAi=Aix+Aiy+Aiz
    do jprim=iprim,nprims
        jcen=b(jprim)%center
        Bexp=b(jprim)%exp
        Bix=type2ix(b(jprim)%functype)
        Biy=type2iy(b(jprim)%functype)
        Biz=type2iz(b(jprim)%functype)
        Bx=a(jcen)%x
        By=a(jcen)%y
        Bz=a(jcen)%z
        sumBi=Bix+Biy+Biz
        ep=Aexp+Bexp
        Px=(Ax*Aexp+Bx*Bexp)/ep
        Py=(Ay*Aexp+By*Bexp)/ep
        Pz=(Az*Aexp+Bz*Bexp)/ep
        PAx=Px-Ax
        PAy=Py-Ay
        PAz=Pz-Az
        PBx=Px-Bx
        PBy=Py-By
        PBz=Pz-Bz
        sqAB=(Ax-Bx)**2+(Ay-By)**2+(Az-Bz)**2
        PCx=Px-Cx
        PCy=Py-Cy
        PCz=Pz-Cz
        sqPC=PCx*PCx+PCy*PCy+PCz*PCz

        tmpval=-Aexp*Bexp*sqAB/ep
        prefac=2*pi/ep*dexp(tmpval)
        
        if (-ep*sqPC>espexpcut) then
            expngaPC=dexp(-ep*sqPC)
        else
            expngaPC=0
        end if
        maxFn=sumAi+sumBi
        Fn(maxFn)=Fmch(maxFn,ep*sqPC,expngaPC)
        nu=maxFn
        twoepsqPC=2*ep*sqPC
        do while (nu>0)
            Fn(nu-1)=(expngaPC+twoepsqPC*Fn(nu))/(2*nu-1) !cook book p280
            nu=nu-1
        end do

        tmpnuml=0
        do l=0,Aix+Bix
            tl=1.0D0
            if (mod(l,2)==1) tl=-1.0D0
            fjtmp=fj(l,Aix,Bix,PAx,PBx)*tl*fact(l)
            do r=0,l/2.0D0
                do i=0,(l-2*r)/2.0D0
                    tmpnuml=tmpnuml+1
                    Alri(tmpnuml)=Afac(l,r,i,PCx,ep,fjtmp)
                    maplri(tmpnuml)=l-2*r-i
                end do
            end do
        end do

        tmpnumm=0
        do m=0,Aiy+Biy
            tm=1.0D0
            if (mod(m,2)==1) tm=-1.0D0
            fjtmp=fj(m,Aiy,Biy,PAy,PBy)*tm*fact(m)
            do s=0,m/2.0D0
                do j=0,(m-2*s)/2.0D0
                    tmpnumm=tmpnumm+1
                    Amsj(tmpnumm)=Afac(m,s,j,PCy,ep,fjtmp)
                    mapmsj(tmpnumm)=m-2*s-j
                end do
            end do
        end do

        tmpnumn=0
        do n=0,Aiz+Biz
            tn=1.0D0
            if (mod(n,2)==1) tn=-1.0D0
            fjtmp=fj(n,Aiz,Biz,PAz,PBz)*tn*fact(n)
            do t=0,n/2.0D0
                do k=0,(n-2*t)/2.0D0
                    tmpnumn=tmpnumn+1
                    Antk(tmpnumn)=Afac(n,t,k,PCz,ep,fjtmp)
                    mapntk(tmpnumn)=n-2*t-k
                end do
            end do
        end do

        term=0.0D0
        !Now calc "term"=<psi(iprim)|1/r_Z|psi(jprim)>
        do l=1,tmpnuml
            do m=1,tmpnumm
                do n=1,tmpnumn
                    term=term+Alri(l)*Amsj(m)*Antk(n)*Fn(maplri(l)+mapmsj(m)+mapntk(n))
                end do
            end do
        end do

        if (iprim/=jprim) term=2.0*term
        term=term*prefac
        addesp=0.0D0
        do imo=1,nmo
            addesp=addesp+MOocc(imo)*CO(imo,iprim)*CO(imo,jprim)
        end do
        espprivate=espprivate+addesp*term
    end do !end j primitive
!$OMP critical
    eleesp=eleesp+espprivate
!$OMP end critical
end do !end i primitive
!$OMP end parallel do
eleesp=-eleesp
end function


!!!------------------------- Calculate ESP in a plane
!maxnumgrid is the maximum value of ngridnum1 and ngridnum2, this value determine the size of Alrivec,Amsjvec,Antkvec
!We don't allocate Alrivec,Amsjvec,Antkvec dynamically since if we do such thing, this routine will crash in win7-64bit system.
!The reason may be that dynamical arrays are not fully compatiable with private property in OpenMP
subroutine planeesp(maxnumgrid)
implicit none
integer maxnumgrid
integer,parameter :: narrmax=396 !Enough for h-type GTF, 5+5=10. Alri(0:10,0:5,0:5)-->11*6*6=396
real*8 Cx,Cy,Cz,Cxold,Cyold,Czold,term,ep,Ax,Ay,Az,Bx,By,Bz,Aexp,Bexp,Px,Py,Pz,prefac,tmpval
real*8 sqPC,sqAB,expngaPC,PAx,PAy,PAz,PBx,PBy,PBz,PCx,PCy,PCz,fjtmp,addesp,espexpcut
real*8 Alri(narrmax),Amsj(narrmax),Antk(narrmax),Fnmat(0:ngridnum1-1,0:ngridnum2-1),Fnvec(0:10) !Enough for h-type GTF, 5+5=10
real*8:: Alrivec(narrmax,maxnumgrid),Amsjvec(narrmax,maxnumgrid),Antkvec(narrmax,maxnumgrid)
real*8 twoepsqPC,tl,tm,tn,pleprivate(ngridnum1,ngridnum2) !Store plane contribution of GTFs in each thread, then sum up
integer nu,imo,iprim,jprim,maxFn,maplri(narrmax),mapmsj(narrmax),mapntk(narrmax),tmpnuml,tmpnumm,tmpnumn
integer Aix,Aiy,Aiz,Bix,Biy,Biz,l,r,i,m,s,j,n,t,k,icen,jcen,sumAi,sumBi,ii,jj,planetype,numx,numy,numz,ifinish
ifinish=0
espexpcut=log10(espprecutoff)*3
planemat=0D0
!Calc ESP of nuclear contribution
do ii=0,ngridnum1-1
    do jj=0,ngridnum2-1
        Cx=orgx2D+ii*v1x+jj*v2x
        Cy=orgy2D+ii*v1y+jj*v2y
        Cz=orgz2D+ii*v1z+jj*v2z
        planemat(ii+1,jj+1)=nucesp(Cx,Cy,Cz)
    end do
end do

nthreads=getNThreads()
!$OMP PARALLEL DO private(Aix,Aiy,Aiz,Bix,Biy,Biz,l,r,i,m,s,j,n,t,k,icen,jcen,sumAi,sumBi,ii,jj,planetype,numx,numy,numz,&
!$OMP nu,imo,iprim,jprim,maxFn,maplri,mapmsj,mapntk,tmpnuml,tmpnumm,tmpnumn,&
!$OMP twoepsqPC,tl,tm,tn,Alrivec,Amsjvec,Antkvec,Alri,Amsj,Antk,Fnmat,Fnvec,&
!$OMP sqPC,sqAB,expngaPC,PAx,PAy,PAz,PBx,PBy,PBz,PCx,PCy,PCz,fjtmp,addesp,&
!$OMP Cx,Cy,Cz,Cxold,Cyold,Czold,term,ep,Ax,Ay,Az,Bx,By,Bz,Aexp,Bexp,Px,Py,Pz,prefac,tmpval,pleprivate) &
!$OMP shared(planemat,ifinish) schedule(dynamic) NUM_THREADS(nthreads)
do iprim=1,nprims
    pleprivate=0D0
    icen=b(iprim)%center
    Aexp=b(iprim)%exp
    Ax=a(icen)%x
    Ay=a(icen)%y
    Az=a(icen)%z
    Aix=type2ix(b(iprim)%functype)
    Aiy=type2iy(b(iprim)%functype)
    Aiz=type2iz(b(iprim)%functype)
    sumAi=Aix+Aiy+Aiz
    do jprim=iprim,nprims
        jcen=b(jprim)%center
        Bexp=b(jprim)%exp
        Bix=type2ix(b(jprim)%functype)
        Biy=type2iy(b(jprim)%functype)
        Biz=type2iz(b(jprim)%functype)
        Bx=a(jcen)%x
        By=a(jcen)%y
        Bz=a(jcen)%z
        ep=Aexp+Bexp
        Px=(Ax*Aexp+Bx*Bexp)/ep
        Py=(Ay*Aexp+By*Bexp)/ep
        Pz=(Az*Aexp+Bz*Bexp)/ep
        PAx=Px-Ax
        PAy=Py-Ay
        PAz=Pz-Az
        PBx=Px-Bx
        PBy=Py-By
        PBz=Pz-Bz
        sqAB=(Ax-Bx)**2+(Ay-By)**2+(Az-Bz)**2
        tmpval=-Aexp*Bexp*sqAB/ep
        prefac=2*pi/ep*dexp(tmpval)
        if (prefac<espprecutoff) cycle
        
        Cxold=999.99912D0 !An arbitrary number
        Cyold=999.99912D0
        Czold=999.99912D0
        sumBi=Bix+Biy+Biz
        maxFn=sumAi+sumBi

        !! Start cycle grid point
        do ii=0,ngridnum1-1
            do jj=0,ngridnum2-1
                Cx=orgx2D+ii*v1x+jj*v2x
                Cy=orgy2D+ii*v1y+jj*v2y
                Cz=orgz2D+ii*v1z+jj*v2z
                PCx=Px-Cx
                PCy=Py-Cy
                PCz=Pz-Cz
                sqPC=PCx*PCx+PCy*PCy+PCz*PCz
                twoepsqPC=2*ep*sqPC
                term=0.0D0
                if (-ep*sqPC>espexpcut) then
                    expngaPC=dexp(-ep*sqPC)
                else
                    expngaPC=0D0
                end if
                Fnmat(ii,jj)=Fmch(maxFn,ep*sqPC,expngaPC)
                Fnvec(maxFn)=Fnmat(ii,jj)
                do nu=maxFn,1,-1
                    Fnvec(nu-1)=(expngaPC+twoepsqPC*Fnvec(nu))/(2*nu-1) !cook book p280
                end do

                if (Cx/=Cxold) then
                    tmpnuml=0
                    do l=0,Aix+Bix
                        tl=1.0D0
                        if (mod(l,2)==1) tl=-1.0D0
                        fjtmp=fj(l,Aix,Bix,PAx,PBx)*tl*fact(l)
                        do r=0,l/2.0D0
                            do i=0,(l-2*r)/2.0D0
                                tmpnuml=tmpnuml+1
                                Alri(tmpnuml)=Afac(l,r,i,PCx,ep,fjtmp)
                                maplri(tmpnuml)=l-2*r-i
                            end do
                        end do
                    end do
                    Cxold=Cx
                end if
                if (Cy/=Cyold) then
                    tmpnumm=0
                    do m=0,Aiy+Biy
                        tm=1.0D0
                        if (mod(m,2)==1) tm=-1.0D0
                        fjtmp=fj(m,Aiy,Biy,PAy,PBy)*tm*fact(m)
                        do s=0,m/2.0D0
                            do j=0,(m-2*s)/2.0D0
                                tmpnumm=tmpnumm+1
                                Amsj(tmpnumm)=Afac(m,s,j,PCy,ep,fjtmp)
                                mapmsj(tmpnumm)=m-2*s-j
                            end do
                        end do
                    end do
                    Cyold=Cy
                end if
                if (Cz/=Czold) then
                    tmpnumn=0
                    do n=0,Aiz+Biz
                        tn=1.0D0
                        if (mod(n,2)==1) tn=-1.0D0
                        fjtmp=fj(n,Aiz,Biz,PAz,PBz)*tn*fact(n)
                        do t=0,n/2.0D0
                            do k=0,(n-2*t)/2.0D0
                                tmpnumn=tmpnumn+1
                                Antk(tmpnumn)=Afac(n,t,k,PCz,ep,fjtmp)
                                mapntk(tmpnumn)=n-2*t-k
                            end do
                        end do
                    end do
                    Czold=Cz
                end if

                !Now calc "term"=<psi(iprim)|1/r_Z|psi(jprim)>
                do l=1,tmpnuml
                    do m=1,tmpnumm
                        do n=1,tmpnumn
                            term=term+Alri(l)*Amsj(m)*Antk(n)*Fnvec(maplri(l)+mapmsj(m)+mapntk(n))
                        end do
                    end do
                end do

                if (iprim/=jprim) term=2.0*term
                term=term*prefac
                addesp=0.0D0
                do imo=1,nmo
                    addesp=addesp+MOocc(imo)*CO(imo,iprim)*CO(imo,jprim)
                end do
                pleprivate(ii+1,jj+1)=pleprivate(ii+1,jj+1)-addesp*term
            end do !end jj cycle
        end do !end ii cycle

    end do !end j primitive
    ifinish=ifinish+1
    write(*,"(' Finished: ',i6,' /',i6)") ifinish,nprims
!$OMP CRITICAL
    planemat=planemat+pleprivate
!$OMP END CRITICAL
end do !end i primitive
!$OMP END PARALLEL DO
end subroutine


!!!------------------------- Calculate grid data of ESP from electrons and store to cubmat
subroutine espcub
implicit none
integer,parameter :: narrmax=396 !Enough for h-type GTF, 5+5=10. Alri(0:10,0:5,0:5)-->11*6*6=396
real*8 Cx,Cy,Cz,term,ep,Ax,Ay,Az,Bx,By,Bz,Aexp,Bexp,Px,Py,Pz,prefac,tmpval
real*8 sqPC,sqAB,expngaPC,PAx,PAy,PAz,PBx,PBy,PBz,PCx,PCy,PCz,fjtmp,addesp
real*8 Alrivec(narrmax,nx),Amsjvec(narrmax,ny),Antkvec(narrmax,nz),Fnmat(0:nx-1,0:ny-1,0:nz-1),Fnvec(0:10) !Enough for h-type GTF, 5+5=10
real*8 twoepsqPC,tl,tm,tn,time_begin,time_end,espexpcut,cubprivate(nx,ny,nz)
integer nu,imo,iprim,jprim,maxFn,maplri(narrmax),mapmsj(narrmax),mapntk(narrmax),tmpnuml,tmpnumm,tmpnumn
integer Aix,Aiy,Aiz,Bix,Biy,Biz,l,r,i,m,s,j,n,t,k,icen,jcen,sumAi,sumBi,ii,jj,kk,ifinish
CALL CPU_TIME ( time_begin )
ifinish=0
espexpcut=log10(espprecutoff)*3
nthreads=getNThreads()
!$OMP PARALLEL DO private(Aix,Aiy,Aiz,Bix,Biy,Biz,l,r,i,m,s,j,n,t,k,icen,jcen,sumAi,sumBi,ii,jj,kk, &
!$OMP nu,imo,iprim,jprim,maxFn,maplri,mapmsj,mapntk,tmpnuml,tmpnumm,tmpnumn, &
!$OMP twoepsqPC,tl,tm,tn,Alrivec,Amsjvec,Antkvec,Fnmat,Fnvec, &
!$OMP sqPC,sqAB,expngaPC,PAx,PAy,PAz,PBx,PBy,PBz,PCx,PCy,PCz,fjtmp,addesp,cubprivate, &
!$OMP Cx,Cy,Cz,term,ep,Ax,Ay,Az,Bx,By,Bz,Aexp,Bexp,Px,Py,Pz,tmpval,prefac) &
!$OMP shared(cubmat,ifinish) schedule(dynamic) NUM_THREADS(nthreads)
do iprim=1,nprims
    cubprivate=0D0
    icen=b(iprim)%center
    Aexp=b(iprim)%exp
    Ax=a(icen)%x
    Ay=a(icen)%y
    Az=a(icen)%z
    Aix=type2ix(b(iprim)%functype)
    Aiy=type2iy(b(iprim)%functype)
    Aiz=type2iz(b(iprim)%functype)
    sumAi=Aix+Aiy+Aiz
    do jprim=iprim,nprims
        jcen=b(jprim)%center
        Bexp=b(jprim)%exp
        Bix=type2ix(b(jprim)%functype)
        Biy=type2iy(b(jprim)%functype)
        Biz=type2iz(b(jprim)%functype)
        Bx=a(jcen)%x
        By=a(jcen)%y
        Bz=a(jcen)%z
        ep=Aexp+Bexp
        Px=(Ax*Aexp+Bx*Bexp)/ep
        Py=(Ay*Aexp+By*Bexp)/ep
        Pz=(Az*Aexp+Bz*Bexp)/ep
        PAx=Px-Ax
        PAy=Py-Ay
        PAz=Pz-Az
        PBx=Px-Bx
        PBy=Py-By
        PBz=Pz-Bz
        sqAB=(Ax-Bx)**2+(Ay-By)**2+(Az-Bz)**2
        tmpval=-Aexp*Bexp*sqAB/ep
        
        prefac=2*pi/ep*dexp(tmpval)
        if (prefac<espprecutoff) cycle
        sumBi=Bix+Biy+Biz
        maxFn=sumAi+sumBi

        do ii=1,nx
            Cx=orgx+(ii-1)*dx
            PCx=Px-Cx
            tmpnuml=0
            do l=0,Aix+Bix
                tl=1.0D0
                if (mod(l,2)==1) tl=-1.0D0
                fjtmp=fj(l,Aix,Bix,PAx,PBx)*tl*fact(l)
                do r=0,l/2.0D0
                    do i=0,(l-2*r)/2.0D0
                        tmpnuml=tmpnuml+1
                        Alrivec(tmpnuml,ii)=Afac(l,r,i,PCx,ep,fjtmp)
                        maplri(tmpnuml)=l-2*r-i
                    end do
                end do
            end do
        end do
        do ii=1,ny
            Cy=orgy+(ii-1)*dy
            PCy=Py-Cy
            tmpnumm=0
            do m=0,Aiy+Biy
                tm=1.0D0
                if (mod(m,2)==1) tm=-1.0D0
                fjtmp=fj(m,Aiy,Biy,PAy,PBy)*tm*fact(m)
                do s=0,m/2.0D0
                    do j=0,(m-2*s)/2.0D0
                        tmpnumm=tmpnumm+1
                        Amsjvec(tmpnumm,ii)=Afac(m,s,j,PCy,ep,fjtmp)
                        mapmsj(tmpnumm)=m-2*s-j
                    end do
                end do
            end do
        end do
        do ii=1,nz
            Cz=orgz+(ii-1)*dz
            PCz=Pz-Cz
            tmpnumn=0
            do n=0,Aiz+Biz
                tn=1.0D0
                if (mod(n,2)==1) tn=-1.0D0
                fjtmp=fj(n,Aiz,Biz,PAz,PBz)*tn*fact(n)
                do t=0,n/2.0D0
                    do k=0,(n-2*t)/2.0D0
                        tmpnumn=tmpnumn+1
                        Antkvec(tmpnumn,ii)=Afac(n,t,k,PCz,ep,fjtmp)
                        mapntk(tmpnumn)=n-2*t-k
                    end do
                end do
            end do
        end do

        !! Start cycle grid point
        do kk=0,nz-1
            do jj=0,ny-1
                do ii=0,nx-1
                    Cx=orgx+ii*dx
                    Cy=orgy+jj*dy
                    Cz=orgz+kk*dz
                    PCx=Px-Cx
                    PCy=Py-Cy
                    PCz=Pz-Cz
                    sqPC=PCx*PCx+PCy*PCy+PCz*PCz
                    twoepsqPC=2*ep*sqPC
                    term=0.0D0
                    if (-ep*sqPC>espexpcut) then
                        expngaPC=dexp(-ep*sqPC)
                    else
                        expngaPC=0D0
                    end if
                    Fnmat(ii,jj,kk)=Fmch(maxFn,ep*sqPC,expngaPC)
                    Fnvec(maxFn)=Fnmat(ii,jj,kk)
                    do nu=maxFn,1,-1
                        Fnvec(nu-1)=(expngaPC+twoepsqPC*Fnvec(nu))/(2*nu-1) !cook book p280
                    end do
                    do l=1,tmpnuml
                        do m=1,tmpnumm
                            do n=1,tmpnumn
                                term=term+Alrivec(l,ii+1)*Amsjvec(m,jj+1)*Antkvec(n,kk+1)*Fnvec(maplri(l)+mapmsj(m)+mapntk(n))
                            end do
                        end do
                    end do

                    if (iprim/=jprim) term=2D0*term
                    term=term*prefac
                    addesp=0.0D0
                    do imo=1,nmo
                        addesp=addesp+MOocc(imo)*CO(imo,iprim)*CO(imo,jprim)
                    end do
                    cubprivate(ii+1,jj+1,kk+1)=cubprivate(ii+1,jj+1,kk+1)-addesp*term
                end do !end ii cycle 
            end do !end jj cycle
        end do !end kk cycle
        
    end do !end j primitive
    ifinish=ifinish+1
    write(*,"(' Finished: ',i6,'/',i6)") ifinish,nprims
!$OMP CRITICAL
    cubmat=cubmat+cubprivate
!$OMP END CRITICAL
end do !end i primitive
!$OMP END PARALLEL DO
CALL CPU_TIME ( time_end )
if (isys==1) write(*,"(' ESP calculation took up CPU time',f12.2,' seconds')") time_end-time_begin
end subroutine






!========================================================================
!============== Utilities routine for function modeule ==================
!========================================================================

!!!---------- Generate nuclear attraction potential integral matrix between all GTFs at a given point
subroutine genGTFattmat(x,y,z,GTFattmat)
use defvar
integer,parameter :: narrmax=396
real*8 x,y,z,GTFattmat(nprims,nprims),Alri(narrmax),Amsj(narrmax),Antk(narrmax),Fn(0:10) !Enough for h-type GTF, 5+5=10. Alri(0:10,0:5,0:5)-->11*6*6=396
integer maxFn,maplri(narrmax),mapmsj(narrmax),mapntk(narrmax),tmpnuml,tmpnumm,tmpnumn,Aix,Aiy,Aiz,Bix,Biy,Biz,r,s,t,sumAi,sumBi
espexpcut=log10(espprecutoff)*3
nthreads=getNThreads()
!$OMP parallel do private(Aix,Aiy,Aiz,Bix,Biy,Biz,l,r,i,m,s,j,n,t,k,icen,jcen,sumAi,sumBi, &
!$OMP nu,iprim,jprim,maxFn,maplri,mapmsj,mapntk,tmpnuml,tmpnumm,tmpnumn,&
!$OMP twoepsqPC,tl,tm,tn,Alri,Amsj,Antk,Fn,sqPC,sqAB,expngaPC,PAx,PAy,PAz,PBx,PBy,PBz,PCx,PCy,PCz,fjtmp,&
!$OMP term,ep,Ax,Ay,Az,Bx,By,Bz,Aexp,Bexp,Px,Py,Pz,prefac,tmpval) shared(GTFattmat) schedule(dynamic) NUM_THREADS(nthreads)
do iprim=1,nprims
    icen=b(iprim)%center
    Aexp=b(iprim)%exp
    Ax=a(icen)%x
    Ay=a(icen)%y
    Az=a(icen)%z
    Aix=type2ix(b(iprim)%functype)
    Aiy=type2iy(b(iprim)%functype)
    Aiz=type2iz(b(iprim)%functype)
    sumAi=Aix+Aiy+Aiz
    do jprim=iprim,nprims
        jcen=b(jprim)%center
        Bexp=b(jprim)%exp
        Bix=type2ix(b(jprim)%functype)
        Biy=type2iy(b(jprim)%functype)
        Biz=type2iz(b(jprim)%functype)
        Bx=a(jcen)%x
        By=a(jcen)%y
        Bz=a(jcen)%z
        sumBi=Bix+Biy+Biz
        ep=Aexp+Bexp
        Px=(Ax*Aexp+Bx*Bexp)/ep
        Py=(Ay*Aexp+By*Bexp)/ep
        Pz=(Az*Aexp+Bz*Bexp)/ep
        PAx=Px-Ax
        PAy=Py-Ay
        PAz=Pz-Az
        PBx=Px-Bx
        PBy=Py-By
        PBz=Pz-Bz
        sqAB=(Ax-Bx)**2+(Ay-By)**2+(Az-Bz)**2
        PCx=Px-x
        PCy=Py-y
        PCz=Pz-z
        sqPC=PCx*PCx+PCy*PCy+PCz*PCz
        tmpval=-Aexp*Bexp*sqAB/ep
        prefac=2*pi/ep*dexp(tmpval)
        if (-ep*sqPC>espexpcut) then
            expngaPC=dexp(-ep*sqPC)
        else
            expngaPC=0
        end if
        maxFn=sumAi+sumBi
        Fn(maxFn)=Fmch(maxFn,ep*sqPC,expngaPC)
        nu=maxFn
        twoepsqPC=2*ep*sqPC
        do while (nu>0)
            Fn(nu-1)=(expngaPC+twoepsqPC*Fn(nu))/(2*nu-1) !cook book p280
            nu=nu-1
        end do
        
        tmpnuml=0
        do l=0,Aix+Bix
            tl=1.0D0
            if (mod(l,2)==1) tl=-1.0D0
            fjtmp=fj(l,Aix,Bix,PAx,PBx)*tl*fact(l)
            do r=0,l/2.0D0
                do i=0,(l-2*r)/2.0D0
                    tmpnuml=tmpnuml+1
                    Alri(tmpnuml)=Afac(l,r,i,PCx,ep,fjtmp)
                    maplri(tmpnuml)=l-2*r-i
                end do
            end do
        end do
        tmpnumm=0
        do m=0,Aiy+Biy
            tm=1.0D0
            if (mod(m,2)==1) tm=-1.0D0
            fjtmp=fj(m,Aiy,Biy,PAy,PBy)*tm*fact(m)
            do s=0,m/2.0D0
                do j=0,(m-2*s)/2.0D0
                    tmpnumm=tmpnumm+1
                    Amsj(tmpnumm)=Afac(m,s,j,PCy,ep,fjtmp)
                    mapmsj(tmpnumm)=m-2*s-j
                end do
            end do
        end do
        tmpnumn=0
        do n=0,Aiz+Biz
            tn=1.0D0
            if (mod(n,2)==1) tn=-1.0D0
            fjtmp=fj(n,Aiz,Biz,PAz,PBz)*tn*fact(n)
            do t=0,n/2.0D0
                do k=0,(n-2*t)/2.0D0
                    tmpnumn=tmpnumn+1
                    Antk(tmpnumn)=Afac(n,t,k,PCz,ep,fjtmp)
                    mapntk(tmpnumn)=n-2*t-k
                end do
            end do
        end do

        term=0D0
        !Now calc "term"=<psi(iprim)|1/r12|psi(jprim)>, r1 is inputted x,y,z, r2 is the integration variable
        do l=1,tmpnuml
            do m=1,tmpnumm
                do n=1,tmpnumn
                    term=term+Alri(l)*Amsj(m)*Antk(n)*Fn(maplri(l)+mapmsj(m)+mapntk(n))
                end do
            end do
        end do
        term=term*prefac
        GTFattmat(iprim,jprim)=term
        GTFattmat(jprim,iprim)=term
    end do !end j primitive
end do !end i primitive
!$OMP end parallel do
end subroutine

!!!------------------------- Calculate A-factor at Cook book p245
real*8 function Afac(l,r,i,PC,gamma,fjtmp)
integer l,r,i,ti
real*8 gamma,PC,comp,PCterm,fjtmp
ti=1.0D0
if (mod(i,2)==1) ti=-1.0D0 !faster than ti=(-1)**i
PCterm=1.0D0
if (l-2*r-2*i/=0) PCterm=PC**(l-2*r-2*i)
comp=ti*PCterm*(0.25D0/gamma)**(r+i) / ( fact(r)*fact(i)*fact(l-2*r-2*i) )
Afac=fjtmp*comp
! comp=(-1)**i*fact(l)*PC**(l-2*r-2*i)*(1/(4*gamma))**(r+i) / ( fact(r)*fact(i)*fact(l-2*r-2*i) )
! Afac=(-1)**l * fj(l,l1,l2,PA,PB)*comp
end function      

!!!------------------------- Calculate fj at Cook book p237
real*8 function fj(j,l,m,aa,bb)
real*8 aa,bb,pre,at,bt
integer j,l,m,k,imin,imax
imax=min(j,l)
imin=max(0,j-m)
fj=0.0D0
do k=imin,imax
    pre=fact(l)/fact(l-k)/fact(k) * fact(m)/fact(m-j+k)/fact(j-k)
    at=1.0D0
    bt=1.0D0
    if (l-k/=0) at=aa**(l-k)  !This determine helps to improve efficient
    if (m+k-j/=0) bt=bb**(m+k-j)
    fj=fj+pre*at*bt
end do
end function

!!!---- Calculate int('t^(2*m)*exp(-x*t^2)','t',0,1) see Cook book p281 for detail
real*8 function Fmch(m,x,expnx)
! expnx is input parameter, value should be exp(-x), because calculate the value is time-consuming and in
! other place this value also need be calculate, so not recalculate in this subroutine
IMPLICIT none
integer m,i
real*8 x,expnx,a,b,term,partsum,APPROX,xd,FIMULT,NOTRMS,eps,fiprop
eps=1.0D-8  !convergence precision
Fmch=0D0
if (x<=10) then
    if (expnx==0D0) RETURN
    A=m+0.5D0
    term=1.0D0/A
    partsum=term
    DO I=2,50
        A=A+1.0D0
        term=term*X/A
        partsum=partsum+term
        if ( term/partsum < eps) THEN
           Fmch = 0.5D0*partsum*expnx
           RETURN
        END IF
    END DO
    write(*,*) "Error: Fmch didn't converge"
else !x is big, use suitable method for solve this situation
    A=M
    B=A+0.5D0
    A=A-0.5D0
    XD=1.0D0/X
    APPROX=0.88622692D0*(dsqrt(XD)*XD**m)
    DO I=1,m
        B=B-1.0D0
        APPROX=APPROX*B
    END DO
    FIMULT=0.5D0*expnx*XD
    partsum=0.D0
    IF (FIMULT==0.0D0) THEN
        Fmch=APPROX-FIMULT*partsum
        return
    ELSE
        FIPROP=FIMULT/APPROX
        term=1.0d0
        partsum=term
        NOTRMS=X
        NOTRMS=NOTRMS+M
        DO I=2,NOTRMS
           term=term*A*XD
           partsum=partsum+term
           IF (dabs(term*FIPROP/partsum)<eps)  THEN
              Fmch=APPROX-FIMULT*partsum
              RETURN
           END IF
           A=A-1.0D0
        END DO
        write(*,*) "Didn't converge"
    END IF
end if
end function




!!!------------- Generate Becke weight function, only used by Sobereva
!If inp2=0, then return atomic space weight of atom inp1, else return overlap weight of atom inp1 and inp2
real*8 function beckewei(x,y,z,inp1,inp2)
use defvar
real*8 x,y,z
integer inp1,inp2
! integer :: fraglist(13)=(/ 24,12,23,4,11,3,22,1,10,2,18,20,21 /)
real*8 smat(ncenter,ncenter),Pvec(ncenter)
smat=1.0D0
do ii=1,ncenter
    ri=dsqrt( (x-a(ii)%x)**2+(y-a(ii)%y)**2+(z-a(ii)%z)**2 )
    do jj=1,ncenter
        if (ii==jj) cycle
        rj=dsqrt( (x-a(jj)%x)**2+(y-a(jj)%y)**2+(z-a(jj)%z)**2 )
        rmiu=(ri-rj)/distmat(ii,jj)

        chi=covr_tianlu(a(ii)%index)/covr_tianlu(a(jj)%index) !Adjust for heteronuclear
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
Pvec=Pvec/sum(Pvec)
! beckewei=Pvec(2)*Pvec(3)*Pvec(4)
if (inp2==0) then
    beckewei=Pvec(inp1)
else
    beckewei=Pvec(inp1)*Pvec(inp2)
end if
! beckewei=sum(Pvec(fraglist)) !Get fragmental Becke weight
! beckewei=sum(Pvec(1:5))
end function


!!!-------- Ellipticity of electron density (itype=1) or eta (itype=2)
real*8 function densellip(x,y,z,itype)
use util
integer itype
real*8 x,y,z,dens,grad(3),hess(3,3),eigval(3),eigvecmat(3,3)
call calchessmat_dens(2,x,y,z,dens,grad,hess)
call diagmat(hess,eigvecmat,eigval,300,1D-12)
eigmax=maxval(eigval) !At bcp, will be the most positive 
eigmin=minval(eigval) !At bcp, will be the most negative
do itmp=1,3
    tmpval=eigval(itmp)
    if (tmpval/=eigmax.and.tmpval/=eigmin) eigmed=tmpval !At bcp, will be the second most negative
end do
if (itype==1) then
    densellip=eigmin/eigmed-1
else
    densellip=abs(eigmin)/eigmax
end if
end function


!!!-------------- Single exponential decay detector (SEDD)
!Originally proposed in ChemPhysChem, 13, 3462 (2012), current implementation is based on the new definition in DORI paper (10.1021/ct500490b)
real*8 function SEDD(x,y,z)
real*8 x,y,z,rho,grad(3),hess(3,3)
call calchessmat_dens(2,x,y,z,rho,grad,hess)
dersqr=sum(grad**2)
tmp1_1=rho*(grad(1)*hess(1,1)+grad(2)*hess(1,2)+grad(3)*hess(1,3))
tmp1_2=grad(1)*dersqr
tmp2_1=rho*(grad(1)*hess(1,2)+grad(2)*hess(2,2)+grad(3)*hess(2,3))
tmp2_2=grad(2)*dersqr
tmp3_1=rho*(grad(1)*hess(1,3)+grad(2)*hess(2,3)+grad(3)*hess(3,3))
tmp3_2=grad(3)*dersqr
eps=4/rho**8*( (tmp1_1-tmp1_2)**2 + (tmp2_1-tmp2_2)**2 + (tmp3_1-tmp3_2)**2 )
SEDD=dlog(1+eps)
end function


!!!-------------- Density Overlap Regions Indicator (DORI) DOI:10.1021/ct500490b
real*8 function DORI(x,y,z)
real*8 x,y,z,rho,grad(3),hess(3,3)
call calchessmat_dens(2,x,y,z,rho,grad,hess)
dersqr=sum(grad**2)
tmp1_1=rho*(grad(1)*hess(1,1)+grad(2)*hess(1,2)+grad(3)*hess(1,3))
tmp1_2=grad(1)*dersqr
tmp2_1=rho*(grad(1)*hess(1,2)+grad(2)*hess(2,2)+grad(3)*hess(2,3))
tmp2_2=grad(2)*dersqr
tmp3_1=rho*(grad(1)*hess(1,3)+grad(2)*hess(2,3)+grad(3)*hess(3,3))
tmp3_2=grad(3)*dersqr
theta=4/dersqr**3*( (tmp1_1-tmp1_2)**2 + (tmp2_1-tmp2_2)**2 + (tmp3_1-tmp3_2)**2 )
DORI=theta/(1+theta)
end function


!!!----- Integrand of LSDA exchange functional
real*8 function xLSDA(x,y,z)
real*8 x,y,z
call gendensgradab(x,y,z,adens,bdens,tdens,agrad,bgrad,tgrad,abgrad)
xLSDA=-3D0/2D0*(3D0/4D0/pi)**(1D0/3D0)*(adens**(4D0/3D0)+bdens**(4D0/3D0))
end function

!!!----- Integrand of Becke88 exchange functional
real*8 function xBecke88(x,y,z)
real*8 x,y,z
call gendensgradab(x,y,z,adens,bdens,tdens,agrad,bgrad,tgrad,abgrad)
adens4d3=adens**(4D0/3D0)
bdens4d3=bdens**(4D0/3D0)
slatercoeff=-3D0/2D0*(3D0/4D0/pi)**(1D0/3D0)
slaterxa=slatercoeff*adens4d3
slaterxb=slatercoeff*bdens4d3
slaterx=slaterxa+slaterxb
redagrad=agrad/adens4d3 !Reduced density gradient
redbgrad=bgrad/bdens4d3
arshredagrad=log(redagrad+dsqrt(redagrad**2+1))
Beckexa=adens4d3*redagrad**2/(1+6*0.0042D0*redagrad*arshredagrad)
arshredbgrad=log(redbgrad+dsqrt(redbgrad**2+1))
Beckexb=bdens4d3*redbgrad**2/(1+6*0.0042D0*redbgrad*arshredbgrad)
Beckex=-0.0042D0*(Beckexa+Beckexb)
xBecke88=slaterx+Beckex
end function

!!!----- Integrand of LYP corelation functional
real*8 function cLYP(x,y,z)
real*8 x,y,z
call gendensgradab(x,y,z,adens,bdens,tdens,agrad,bgrad,tgrad,abgrad)
parma=0.04918D0
parmb=0.132D0
parmc=0.2533D0
parmd=0.349D0
densn1d3=tdens**(-1D0/3D0)
parmw=exp(-parmc*densn1d3) / (1+parmd*densn1d3) * tdens**(-11D0/3D0)
parmdel=parmc*densn1d3+parmd*densn1d3/(1+parmd*densn1d3)
parmCf=3D0/10D0*(3*pi*pi)**(2D0/3D0)
tmp1=-parma*4D0/(1+parmd*densn1d3)*adens*bdens/tdens
tmp2a=2**(11D0/3D0)*parmCf*(adens**(8D0/3D0)+bdens**(8D0/3D0))
tmp2b=(47D0/18D0-7D0/18D0*parmdel)*tgrad**2
tmp2c=-(2.5D0-parmdel/18D0)*(agrad**2+bgrad**2)
tmp2d=-(parmdel-11D0)/9D0*(adens/tdens*agrad**2+bdens/tdens*bgrad**2)
tmp2=adens*bdens*(tmp2a+tmp2b+tmp2c+tmp2d)
tmp3=-2D0/3D0*tdens**2*tgrad**2+(2D0/3D0*tdens**2-adens**2)*bgrad**2+(2D0/3D0*tdens**2-bdens**2)*agrad**2
cLYP=tmp1-parma*parmb*parmw*(tmp2+tmp3)
end function


!!!--- Integrand of Weizsacker (steric energy)
! DO NOT consider EDF, because the result outputted by quantum chemistry programs is always for only valence electrons!
real*8 function weizsacker(x,y,z)
real*8 x,y,z,wfnval(nmo),wfnderv(3,nmo),gradrho(3) !,EDFgrad(3)
call orbderv(2,1,nmo,x,y,z,wfnval,wfnderv)
rho=0D0
gradrho=0D0
do i=1,nmo
    rho=rho+MOocc(i)*wfnval(i)**2
    gradrho(:)=gradrho(:)+MOocc(i)*wfnval(i)*wfnderv(:,i)
end do
gradrho=2*gradrho
! if (allocated(b_EDF)) then
!     call EDFrho(2,x,y,z,EDFdens,EDFgrad)
!     rho=rho+EDFdens
!     gradrho=gradrho+EDFgrad
! end if
if (rho<1D-30) then
    weizsacker=0
else
    weizsacker=sum(gradrho(:)**2)/8/rho
end if
end function
!!!--- Weizsacker potential
real*8 function weizpot(x,y,z)
real*8 x,y,z,wfnval(nmo),wfnderv(3,nmo),wfnhess(3,3,nmo),gradrho(3),laplx,laply,laplz,lapltot
rho=0D0
gradrho=0D0
laplx=0D0
laply=0D0
laplz=0D0
call orbderv(3,1,nmo,x,y,z,wfnval,wfnderv,wfnhess)
do i=1,nmo
    rho=rho+MOocc(i)*wfnval(i)**2
    gradrho(:)=gradrho(:)+MOocc(i)*wfnval(i)*wfnderv(:,i)
    laplx=laplx+MOocc(i)*( wfnderv(1,i)**2 + wfnval(i)*wfnhess(1,1,i) )
    laply=laply+MOocc(i)*( wfnderv(2,i)**2 + wfnval(i)*wfnhess(2,2,i) )
    laplz=laplz+MOocc(i)*( wfnderv(3,i)**2 + wfnval(i)*wfnhess(3,3,i) )
end do
gradrho=2*gradrho
lapltot=2*(laplx+laply+laplz)
weizpot=sum(gradrho(:)**2)/8D0/rho**2-lapltot/4D0/rho
end function


!!!--- Steric potential (J. Chem. Phys., 126, 244103), which negative value is one-electron potential
real*8 function stericpot(x,y,z)
real*8 x,y,z,gradrho(3),hessrho(3,3),lapltot
call calchessmat_dens(2,x,y,z,rho,gradrho,hessrho)
lapltot=hessrho(1,1)+hessrho(2,2)+hessrho(3,3)
if (rho<steric_potcutrho) then
    stericpot=steric_potcons
    return
end if
stericpot=sum(gradrho**2)/(rho+steric_addminimal)**2/8D0-lapltot/(rho+steric_addminimal)/4D0
end function

!!!----- Calculate the first-order derivative of steric potential
subroutine stericderv(x,y,z,derv)
real*8 x,y,z,derv(3)
real*8 eleval,elegrad(3),elehess(3,3),laplval,laplgrad(3),rhotens3(3,3,3)
real*8 wfnval(nmo),wfnderv(3,nmo),wfnhess(3,3,nmo),wfntens3(3,3,3,nmo)
real*8 EDFval,EDFgrad(3),EDFhess(3,3),EDFtens3(3,3,3)
call orbderv(5,1,nmo,x,y,z,wfnval,wfnderv,wfnhess,wfntens3)
eleval=sum( MOocc(1:nmo)*wfnval(1:nmo)*wfnval(1:nmo) )
elegrad(1)=2*sum( MOocc(1:nmo)*wfnval(1:nmo)*wfnderv(1,1:nmo) )
elegrad(2)=2*sum( MOocc(1:nmo)*wfnval(1:nmo)*wfnderv(2,1:nmo) )
elegrad(3)=2*sum( MOocc(1:nmo)*wfnval(1:nmo)*wfnderv(3,1:nmo) )
elehess(1,1)=2*sum( MOocc(1:nmo)*( wfnderv(1,1:nmo)**2 + wfnval(1:nmo)*wfnhess(1,1,1:nmo) ) )
elehess(2,2)=2*sum( MOocc(1:nmo)*( wfnderv(2,1:nmo)**2 + wfnval(1:nmo)*wfnhess(2,2,1:nmo) ) )
elehess(3,3)=2*sum( MOocc(1:nmo)*( wfnderv(3,1:nmo)**2 + wfnval(1:nmo)*wfnhess(3,3,1:nmo) ) )
elehess(1,2)=2*sum( MOocc(1:nmo)*( wfnderv(1,1:nmo)*wfnderv(2,1:nmo)+wfnhess(1,2,1:nmo)*wfnval(1:nmo) ) )
elehess(2,3)=2*sum( MOocc(1:nmo)*( wfnderv(2,1:nmo)*wfnderv(3,1:nmo)+wfnhess(2,3,1:nmo)*wfnval(1:nmo) ) )
elehess(1,3)=2*sum( MOocc(1:nmo)*( wfnderv(1,1:nmo)*wfnderv(3,1:nmo)+wfnhess(1,3,1:nmo)*wfnval(1:nmo) ) )
elehess(2,1)=elehess(1,2)
elehess(3,2)=elehess(2,3)
elehess(3,1)=elehess(1,3)
laplval=elehess(1,1)+elehess(2,2)+elehess(3,3)
rhotens3=0D0
do i=1,nmo
    rhotens3(1,1,1)=rhotens3(1,1,1)+MOocc(i)*( 3*wfnderv(1,i)*wfnhess(1,1,i)+wfnval(i)*wfntens3(1,1,1,i) )
    rhotens3(2,2,2)=rhotens3(2,2,2)+MOocc(i)*( 3*wfnderv(2,i)*wfnhess(2,2,i)+wfnval(i)*wfntens3(2,2,2,i) )
    rhotens3(3,3,3)=rhotens3(3,3,3)+MOocc(i)*( 3*wfnderv(3,i)*wfnhess(3,3,i)+wfnval(i)*wfntens3(3,3,3,i) )
    rhotens3(1,1,2)=rhotens3(1,1,2)+MOocc(i)*( 2*wfnderv(1,i)*wfnhess(1,2,i)+wfnderv(2,i)*wfnhess(1,1,i)+wfnval(i)*wfntens3(1,1,2,i) )
    rhotens3(1,1,3)=rhotens3(1,1,3)+MOocc(i)*( 2*wfnderv(1,i)*wfnhess(1,3,i)+wfnderv(3,i)*wfnhess(1,1,i)+wfnval(i)*wfntens3(1,1,3,i) )
    rhotens3(2,2,3)=rhotens3(2,2,3)+MOocc(i)*( 2*wfnderv(2,i)*wfnhess(2,3,i)+wfnderv(3,i)*wfnhess(2,2,i)+wfnval(i)*wfntens3(2,2,3,i) )
    rhotens3(1,2,2)=rhotens3(1,2,2)+MOocc(i)*( 2*wfnderv(2,i)*wfnhess(2,1,i)+wfnderv(1,i)*wfnhess(2,2,i)+wfnval(i)*wfntens3(2,2,1,i) )
    rhotens3(1,3,3)=rhotens3(1,3,3)+MOocc(i)*( 2*wfnderv(3,i)*wfnhess(3,1,i)+wfnderv(1,i)*wfnhess(3,3,i)+wfnval(i)*wfntens3(3,3,1,i) )
    rhotens3(2,3,3)=rhotens3(2,3,3)+MOocc(i)*( 2*wfnderv(3,i)*wfnhess(3,2,i)+wfnderv(2,i)*wfnhess(3,3,i)+wfnval(i)*wfntens3(3,3,2,i) )
end do
rhotens3=rhotens3*2D0
laplgrad(1)=rhotens3(1,1,1)+rhotens3(1,2,2)+rhotens3(1,3,3)
laplgrad(2)=rhotens3(1,1,2)+rhotens3(2,2,2)+rhotens3(2,3,3)
laplgrad(3)=rhotens3(1,1,3)+rhotens3(2,2,3)+rhotens3(3,3,3)
if (allocated(b_EDF)) then
    call EDFrho(5,x,y,z,EDFval,EDFgrad,EDFhess,EDFtens3)
    eleval=eleval+EDFval
    elegrad=elegrad+EDFgrad
    elehess=elehess+EDFhess    
    laplgrad(1)=laplgrad(1)+EDFtens3(1,1,1)+EDFtens3(1,2,2)+EDFtens3(1,3,3)
    laplgrad(2)=laplgrad(2)+EDFtens3(1,1,2)+EDFtens3(2,2,2)+EDFtens3(2,3,3)
    laplgrad(3)=laplgrad(3)+EDFtens3(1,1,3)+EDFtens3(2,2,3)+EDFtens3(3,3,3)
end if
! Above codes can be simplified as below two lines, but will consume additional 1/3 time
! call calchessmat_dens(2,x,y,z,eleval,elegrad,elehess)
! call calchessmat_lapl(1,x,y,z,laplval,laplgrad,laplhess)

eleval=eleval+steric_addminimal
elenorm2=sum(elegrad**2) !Square of norm of electron density gradient
do i=1,3 !x,y,z
    tmp1=-elenorm2/eleval**3*elegrad(i)+laplval/eleval**2*elegrad(i)-laplgrad(i)/eleval
    tmp2=(elegrad(1)*elehess(1,i)+elegrad(2)*elehess(2,i)+elegrad(3)*elehess(3,i)) / eleval**2
    derv(i)=(tmp1+tmp2)/4D0
end do
end subroutine

!---- Magnitude of steric force
real*8 function stericforce(x,y,z)
real*8 x,y,z,derv(3)
call stericderv(x,y,z,derv)
stericforce=dsqrt(sum(derv**2))
end function


!memo: Shubin reported that steric potential/force is quite sensitive to steric_addminimal, &
!so 2016-Sep-30 I devised another solution to solve the diverse behavior of steric potential/force via Becke's damping function
!!!--- Steric potential with damping function to zero
real*8 function stericpot_damp(x,y,z)
real*8 x,y,z,gradrho(3),hessrho(3,3),lapltot
call calchessmat_dens(2,x,y,z,rho,gradrho,hessrho)
lapltot=hessrho(1,1)+hessrho(2,2)+hessrho(3,3)
stericpotorg=sum(gradrho**2)/rho**2/8D0-lapltot/rho/4D0
weiwidth=2 !e.g. weiwidth=2 and steric_potcutrho=-13 means the Becke damping function of [1,0] corresponds to 1D-11~1D-15
tmps=-(dlog(rho)-steric_potcutrho)/weiwidth
if (tmps<-1) then
    consorg=1
else if (tmps>1) then
    consorg=0
else
    do iter=1,2
        tmps=1.5D0*(tmps)-0.5D0*(tmps)**3
    end do
    consorg=0.5D0*(1-tmps) !The weight to switch to constant value steric_potcons
end if
stericpot_damp=stericpotorg*consorg+steric_potcons*(1-consorg)
end function
!!!--- Steric force based on damped steric potential
real*8 function stericforce_damp(x,y,z)
real*8 x,y,z,derv(3)
diffstep=2D-5
derv(1)=(stericpot_damp(x+diffstep,y,z)-stericpot_damp(x-diffstep,y,z))/(2*diffstep)
derv(2)=(stericpot_damp(x,y+diffstep,z)-stericpot_damp(x,y-diffstep,z))/(2*diffstep)
derv(3)=(stericpot_damp(x,y,z+diffstep)-stericpot_damp(x,y,z-diffstep))/(2*diffstep)
stericforce_damp=dsqrt(sum(derv**2))
end function
!!!--- Steric force directly damped to zero rather than based on damped steric potential
real*8 function stericforce_directdamp(x,y,z)
real*8 x,y,z
weiwidth=2
tmps=-(dlog(fdens(x,y,z))-steric_potcutrho)/weiwidth
if (tmps<-1) then
    consorg=1
else if (tmps>1) then
    consorg=0
else
    do iter=1,2
        tmps=1.5D0*(tmps)-0.5D0*(tmps)**3
    end do
    consorg=0.5D0*(1-tmps)
end if
steric_addminimalold=steric_addminimal
steric_addminimal=0
stericforce_directdamp=stericforce(x,y,z)*consorg
steric_addminimal=steric_addminimalold
end function



!!!----- Steric charge, =lapl(steric potential)/(-4*pi)
! Based on analytical first-order derivative, using finite difference to obtain d2v/dx2, d2v/dy2 and d2v/dz2
real*8 function stericcharge(x,y,z)
real*8 x,y,z,derv1add(3),derv1min(3)
if (fdens(x,y,z)<steric_potcutrho) then
    stericcharge=0D0
    return
end if
diffstep=2D-5
call stericderv(x+diffstep,y,z,derv1add)
call stericderv(x-diffstep,y,z,derv1min)
derv2x=(derv1add(1)-derv1min(1))/(2*diffstep) !d2v/dx2
call stericderv(x,y+diffstep,z,derv1add)
call stericderv(x,y-diffstep,z,derv1min)
derv2y=(derv1add(2)-derv1min(2))/(2*diffstep) !d2v/dy2
call stericderv(x,y,z+diffstep,derv1add)
call stericderv(x,y,z-diffstep,derv1min)
derv2z=(derv1add(3)-derv1min(3))/(2*diffstep) !d2v/dz2
stericcharge=-(derv2x+derv2y+derv2z)/4D0/pi
end function


!!!----- Calculate Fisher information density, see JCP,126,191107 for example
!itype=1 Normal definition
!itype=2 Second Fisher information density, Eq.5 of JCP,126,191107
real*8 function Fisherinfo(itype,x,y,z)
real*8 x,y,z,eleval,elegrad(3),elehess(3,3)
integer itype
Fisherinfo=0
eleval=fdens(x,y,z)
if (eleval<=1D-30) return
if (itype==1) Fisherinfo=fgrad(x,y,z,'t')**2/eleval
if (itype==2) Fisherinfo=-flapl(x,y,z,'t')*log(eleval)
end function


!!!----- Ghosh entropy density, PNAS, 81, 8028
!If itype==1, G(r) will be used as kinetic energy density
!If itype==2, G(r)-der2rho/8 will be used instead, which is the kinetic energy density exactly corresponding to Eq. 22 of PNAS, 81, 8028.
real*8 function Ghoshentro(x,y,z,itype)
integer itype
real*8 kintot,x,y,z,wfnval(nmo),wfnderv(3,nmo),wfnhess(3,3,nmo)
if (itype==1) call orbderv(2,1,nmo,x,y,z,wfnval,wfnderv)
if (itype==2) call orbderv(3,1,nmo,x,y,z,wfnval,wfnderv,wfnhess) !If K(r) is used, use this!!!
rho=0D0
do i=1,nmo
    rho=rho+MOocc(i)*wfnval(i)**2
end do
ck=2.871234D0
TFkin=ck*rho**(5D0/3D0)
kintot=0D0
!   If we use Hamiltonian kinetic density
! hamx=sum( MOocc(1:nmo)*wfnhess(1,1,1:nmo)*wfnval(1:nmo) )
! hamy=sum( MOocc(1:nmo)*wfnhess(2,2,1:nmo)*wfnval(1:nmo) )
! hamz=sum( MOocc(1:nmo)*wfnhess(3,3,1:nmo)*wfnval(1:nmo) )
! kintot=-(hamx+hamy+hamz)/2
!   If we use Lagrangian kinetic density G(r)
do i=1,nmo
    kintot=kintot+MOocc(i)*sum(wfnderv(:,i)**2)
end do
kintot=kintot/2D0
if (itype==2) then
    xlapl=2*sum( MOocc(1:nmo)*( wfnderv(1,1:nmo)**2 + wfnval(1:nmo)*wfnhess(1,1,1:nmo) ) )
    ylapl=2*sum( MOocc(1:nmo)*( wfnderv(2,1:nmo)**2 + wfnval(1:nmo)*wfnhess(2,2,1:nmo) ) )
    zlapl=2*sum( MOocc(1:nmo)*( wfnderv(3,1:nmo)**2 + wfnval(1:nmo)*wfnhess(3,3,1:nmo) ) )
    kintot=kintot-(xlapl+ylapl+zlapl)/8
end if
! if (kintot<0) then
!     write(*,"(5f16.10)") kintot,TFkin,x,y,z
!     read(*,*)
! end if
rlambda=5D0/3D0+log(4D0*pi*ck/3D0)
if (kintot<0) then
    rlogterm=0
else
    rlogterm=log(kintot/TFkin)
end if
Ghoshentro=1.5D0*rho*(rlambda+rlogterm)
end function



!!!----- Pauli potential, see Shubin's paper: Comp. Theor. Chem., 1006, 92-99
!Only suitable for close-shell cases
real*8 function paulipot(x,y,z)
real*8 x,y,z
paulipot=totesp(x,y,z)-DFTxcpot(x,y,z)-weizpot(x,y,z) !Note that the sign of ESP in shubin's CTC paper is inversed
end function

!!!----- The magnitude of Pauli force, namely the gradient norm of negative Pauli potential
real*8 function pauliforce(x,y,z)
real*8 x,y,z
diff=2D-5
forcex=-(paulipot(x+diff,y,z)-paulipot(x-diff,y,z))/(2*diff)
forcey=-(paulipot(x,y+diff,z)-paulipot(x,y-diff,z))/(2*diff)
forcez=-(paulipot(x,y,z+diff)-paulipot(x,y,z-diff))/(2*diff)
pauliforce=dsqrt(forcex**2+forcey**2+forcez**2)
end function

!!!----- Pauli charge, =lapl(Pauli potential)/(-4*pi)
real*8 function paulicharge(x,y,z)
real*8 x,y,z
diff=2D-4 !Should not be smaller, otherwise some dirty points will be presented
value=paulipot(x,y,z)
valuexaddadd=paulipot(x+2*diff,y,z)
valuexminmin=paulipot(x-2*diff,y,z)
valueyaddadd=paulipot(x,y+2*diff,z)
valueyminmin=paulipot(x,y-2*diff,z)
valuezaddadd=paulipot(x,y,z+2*diff)
valuezminmin=paulipot(x,y,z-2*diff)
xcomp=(valuexaddadd-2*value+valuexminmin)/(2*diff)**2
ycomp=(valueyaddadd-2*value+valueyminmin)/(2*diff)**2
zcomp=(valuezaddadd-2*value+valuezminmin)/(2*diff)**2
paulicharge=(xcomp+ycomp+zcomp)/(-4*pi)
end function

!!!----- Quantum potential
real*8 function quantumpot(x,y,z)
real*8 x,y,z
quantumpot=totesp(x,y,z)-weizpot(x,y,z)
end function

!!!----- The magnitude of quantum force, namely the gradient norm of quantum potential
real*8 function quantumforce(x,y,z)
real*8 x,y,z
diff=2D-5
forcex=-(quantumpot(x+diff,y,z)-quantumpot(x-diff,y,z))/(2*diff)
forcey=-(quantumpot(x,y+diff,z)-quantumpot(x,y-diff,z))/(2*diff)
forcez=-(quantumpot(x,y,z+diff)-quantumpot(x,y,z-diff))/(2*diff)
quantumforce=dsqrt(forcex**2+forcey**2+forcez**2)
end function

!!!----- Quantum charge
real*8 function quantumcharge(x,y,z)
real*8 x,y,z
diff=2D-4 !Should not be smaller, otherwise some dirty points will be presented
value=quantumpot(x,y,z)
valuexaddadd=quantumpot(x+2*diff,y,z)
valuexminmin=quantumpot(x-2*diff,y,z)
valueyaddadd=quantumpot(x,y+2*diff,z)
valueyminmin=quantumpot(x,y-2*diff,z)
valuezaddadd=quantumpot(x,y,z+2*diff)
valuezminmin=quantumpot(x,y,z-2*diff)
xcomp=(valuexaddadd-2*value+valuexminmin)/(2*diff)**2
ycomp=(valueyaddadd-2*value+valueyminmin)/(2*diff)**2
zcomp=(valuezaddadd-2*value+valuezminmin)/(2*diff)**2
quantumcharge=(xcomp+ycomp+zcomp)/(-4*pi)
end function

!!!----- The magnitude of electrostatic force, namely the gradient norm of electrostatic potential
real*8 function elestatforce(x,y,z)
real*8 x,y,z
diff=2D-5
forcex=-(totesp(x+diff,y,z)-totesp(x-diff,y,z))/(2*diff)
forcey=-(totesp(x,y+diff,z)-totesp(x,y-diff,z))/(2*diff)
forcez=-(totesp(x,y,z+diff)-totesp(x,y,z-diff))/(2*diff)
elestatforce=dsqrt(forcex**2+forcey**2+forcez**2)
end function

!!!----- Electrostatic charge
real*8 function elestatcharge(x,y,z)
real*8 x,y,z
diff=2D-4 !Should not be smaller, otherwise some dirty points will be presented
value=totesp(x,y,z)
valuexaddadd=totesp(x+2*diff,y,z)
valuexminmin=totesp(x-2*diff,y,z)
valueyaddadd=totesp(x,y+2*diff,z)
valueyminmin=totesp(x,y-2*diff,z)
valuezaddadd=totesp(x,y,z+2*diff)
valuezminmin=totesp(x,y,z-2*diff)
xcomp=(valuexaddadd-2*value+valuexminmin)/(2*diff)**2
ycomp=(valueyaddadd-2*value+valueyminmin)/(2*diff)**2
zcomp=(valuezaddadd-2*value+valuezminmin)/(2*diff)**2
elestatcharge=(xcomp+ycomp+zcomp)/(-4*pi)
end function



!!!---- Use trilinear interpolation to obtain value at a given point by using cubmat
!itype==1: interpolate from cubmat, =2: from cubmattmp
real*8 function linintp3d(x,y,z,itype)
real*8 x,y,z
integer itype
character*80 c80tmp
do ix=1,nx
    x1=orgx+(ix-1)*dx
    x2=orgx+ix*dx
    if (x>=x1.and.x<x2) exit  !1D-10 is used to avoid numerical uncertainty
end do
do iy=1,ny
    y1=orgy+(iy-1)*dy
    y2=orgy+iy*dy
    if (y>=y1.and.y<y2) exit
end do
do iz=1,nz
    z1=orgz+(iz-1)*dz
    z2=orgz+iz*dz
    if (z>=z1.and.z<z2) exit
end do
if (ix>=nx.or.iy>=ny.or.iz>=nz) then !Out of grid data range
    linintp3d=0D0
else
    if (itype==1) then
        valxy1=( cubmat(ix,iy,iz  )*(x2-x)*(y2-y) + cubmat(ix+1,iy,iz  )*(x-x1)*(y2-y) + &
            cubmat(ix,iy+1,iz  )*(x2-x)*(y-y1) + cubmat(ix+1,iy+1,iz  )*(x-x1)*(y-y1) ) /dx/dy
        valxy2=( cubmat(ix,iy,iz+1)*(x2-x)*(y2-y) + cubmat(ix+1,iy,iz+1)*(x-x1)*(y2-y) + &
            cubmat(ix,iy+1,iz+1)*(x2-x)*(y-y1) + cubmat(ix+1,iy+1,iz+1)*(x-x1)*(y-y1) ) /dx/dy
    else
        valxy1=( cubmattmp(ix,iy,iz  )*(x2-x)*(y2-y) + cubmattmp(ix+1,iy,iz  )*(x-x1)*(y2-y) + &
            cubmattmp(ix,iy+1,iz  )*(x2-x)*(y-y1) + cubmattmp(ix+1,iy+1,iz  )*(x-x1)*(y-y1) ) /dx/dy
        valxy2=( cubmattmp(ix,iy,iz+1)*(x2-x)*(y2-y) + cubmattmp(ix+1,iy,iz+1)*(x-x1)*(y2-y) + &
            cubmattmp(ix,iy+1,iz+1)*(x2-x)*(y-y1) + cubmattmp(ix+1,iy+1,iz+1)*(x-x1)*(y-y1) ) /dx/dy
    end if
    linintp3d=valxy1+(z-z1)*(valxy2-valxy1)/dz
end if
end function


!!!---- Trilinear interpolation of 3D-vector field by using cubmatvec
subroutine linintp3dvec(x,y,z,vecintp)
real*8 x,y,z,vecintp(3),valxy1(3),valxy2(3)
do ix=1,nx
    x1=orgx+(ix-1)*dx
    x2=orgx+ix*dx
    if (x>=x1.and.x<x2) exit
end do
do iy=1,ny
    y1=orgy+(iy-1)*dy
    y2=orgy+iy*dy
    if (y>=y1.and.y<y2) exit
end do
do iz=1,nz
    z1=orgz+(iz-1)*dz
    z1=orgz+iz*dz
    if (z>=z1.and.z<z2) exit
end do
if (ix>nx.or.iy>ny.or.iz>nz) then !Out of grid data range
    vecintp=0D0
else
    valxy1(:)=( cubmatvec(:,ix,iy,iz  )*(x2-x)*(y2-y) + cubmatvec(:,ix+1,iy,iz  )*(x-x1)*(y2-y) + &
        cubmatvec(:,ix,iy+1,iz  )*(x2-x)*(y-y1) + cubmatvec(:,ix+1,iy+1,iz  )*(x-x1)*(y-y1) ) /dx/dy
    valxy2(:)=( cubmatvec(:,ix,iy,iz+1)*(x2-x)*(y2-y) + cubmatvec(:,ix+1,iy,iz+1)*(x-x1)*(y2-y) + &
        cubmatvec(:,ix,iy+1,iz+1)*(x2-x)*(y-y1) + cubmatvec(:,ix+1,iy+1,iz+1)*(x-x1)*(y-y1) ) /dx/dy
    vecintp=valxy1+(z-z1)*(valxy2-valxy1)/dz
end if
end subroutine




!!!---- Calculate various kinds of integrand of DFT exchange-correlation functionals
!The routines are provided by DFT repository (ftp://ftp.dl.ac.uk/qcg/dft_library/index.html)
!The global variable "iDFTxcsel" is used to select the XC functional, see manual
!Note that the inner core density represented by EDF field is not taken into account
real*8 function DFTxcfunc(x,y,z)
real*8 x,y,z
if (wfntype==0.or.wfntype==3) then !Close-shell
    DFTxcfunc=DFTxcfunc_close(x,y,z)
else !Open-shell
    DFTxcfunc=DFTxcfunc_open(x,y,z)
end if
end function
!---- Calculate various kinds of DFT exchange-correlation potentials, see the comment of DFTxcfunc
real*8 function DFTxcpot(x,y,z)
real*8 x,y,z
if (wfntype==0.or.wfntype==3) then !Close-shell
    DFTxcpot=DFTxcpot_close(x,y,z)
else !Open-shell
    DFTxcpot=DFTxcpot_open(x,y,z)
    write(*,*) "XC potential for open-shell has not been supported yet!"
    read(*,*)
end if
end function



!!---Close-shell form of DFTxcfunc routine
real*8 function DFTxcfunc_close(x,y,z)
real*8 x,y,z
call gendensgradab(x,y,z,adens,bdens,tdens,agrad,bgrad,tgrad,abgrad)
call getXCdata_close(0,tdens,tgrad**2,DFTxcfunc_close,rnouse,rnouse,rnouse,rnouse,rnouse)
end function

!!---Close-shell form of DFTxcpot routine
real*8 function DFTxcpot_close(x,y,z)
real*8 x,y,z,wfnval(nmo),wfnderv(3,nmo),wfnhess(3,3,nmo),gradrho(3),hessrho(3,3),tmparr(3,1),tmpval(1,1),lapltot
rho=0D0
gradrho=0D0
call orbderv(4,1,nmo,x,y,z,wfnval,wfnderv,wfnhess)
do i=1,nmo
    rho=rho+MOocc(i)*wfnval(i)**2
    gradrho(:)=gradrho(:)+MOocc(i)*wfnval(i)*wfnderv(:,i)
end do
gradrho=2*gradrho
hessrho(1,1)=2*sum( MOocc(1:nmo)*( wfnderv(1,1:nmo)**2 + wfnval(1:nmo)*wfnhess(1,1,1:nmo) ) )
hessrho(2,2)=2*sum( MOocc(1:nmo)*( wfnderv(2,1:nmo)**2 + wfnval(1:nmo)*wfnhess(2,2,1:nmo) ) )
hessrho(3,3)=2*sum( MOocc(1:nmo)*( wfnderv(3,1:nmo)**2 + wfnval(1:nmo)*wfnhess(3,3,1:nmo) ) )
hessrho(1,2)=2*sum( MOocc(1:nmo)*( wfnderv(1,1:nmo)*wfnderv(2,1:nmo)+wfnhess(1,2,1:nmo)*wfnval(1:nmo) ) )
hessrho(2,3)=2*sum( MOocc(1:nmo)*( wfnderv(2,1:nmo)*wfnderv(3,1:nmo)+wfnhess(2,3,1:nmo)*wfnval(1:nmo) ) )
hessrho(1,3)=2*sum( MOocc(1:nmo)*( wfnderv(1,1:nmo)*wfnderv(3,1:nmo)+wfnhess(1,3,1:nmo)*wfnval(1:nmo) ) )
hessrho(2,1)=hessrho(1,2)
hessrho(3,2)=hessrho(2,3)
hessrho(3,1)=hessrho(1,3)
lapltot=hessrho(1,1)+hessrho(2,2)+hessrho(3,3)
sigma=sum(gradrho(:)**2)
tmparr(:,1)=gradrho
tmpval=matmul(2*matmul(transpose(tmparr),hessrho),tmparr) !dot product between grad(sigma) and grad(rho)
call getXCdata_close(2,rho,sigma,value,d1rho,d1sig,d2rho,d2rhosig,d2sig)
DFTxcpot_close=d1rho-2*(d2rhosig*sigma+d2sig*tmpval(1,1)+d1sig*lapltot)
end function



!!!For close-shell cases. Input rho and gradrho^2, return the value and its derivative of selected XC (or X/C only) functional
!The global variable "iDFTxcsel" is used to select the XC functional, see manual
!ixcderv=0: only get value, =1: also get d1rho and d1sig, =2: also get d2rho, d2rhosig and d2sig
!rho, sigma: inputted rho and gradrho^2
!value: outputted integrand of functional
!d1rho, d1sig: 1st derivative of functional w.r.t. rho and sigma, respectively
!d2rho, d2sig: 2nd derivative of functional w.r.t. rho and sigma, respectively
!d2rhosig: 1st derv w.r.t. rho and 1st derv w.r.t. sigma
subroutine getXCdata_close(ixcderv,rho,sigma,value,d1rho,d1sig,d2rho,d2rhosig,d2sig)
integer ixcderv
real*8 rho,sigma,value,d1rho,d1sig,d2rho,d2rhosig,d2sig,rhoa1(1),sigmaaa1(1)
real*8 XCzk(1),Xzk(1),Czk(1),XCvrhoa(1),Xvrhoa(1),Cvrhoa(1),XCvsigmaaa(1),Xvsigmaaa(1),Cvsigmaaa(1)
real*8 XCv2rhoa2(1),Xv2rhoa2(1),Cv2rhoa2(1),XCv2rhoasigmaaa(1),Xv2rhoasigmaaa(1),Cv2rhoasigmaaa(1)
real*8 XCv2sigmaaa2(1),Xv2sigmaaa2(1),Cv2sigmaaa2(1)
rhoa1(1)=rho
sigmaaa1(1)=sigma
!X part
if (iDFTxcsel==0.or.iDFTxcsel==80) then
    call rks_x_lda(ixcderv,1,rhoa1,sigmaaa1,Xzk,Xvrhoa,Xvsigmaaa,Xv2rhoa2,Xv2rhoasigmaaa,Xv2sigmaaa2)
else if (iDFTxcsel==1.or.iDFTxcsel==81.or.iDFTxcsel==82.or.iDFTxcsel==83) then
    call rks_x_b88(ixcderv,1,rhoa1,sigmaaa1,Xzk,Xvrhoa,Xvsigmaaa,Xv2rhoa2,Xv2rhoasigmaaa,Xv2sigmaaa2)
else if (iDFTxcsel==2.or.iDFTxcsel==84) then
    call rks_x_pbe(ixcderv,1,rhoa1,sigmaaa1,Xzk,Xvrhoa,Xvsigmaaa,Xv2rhoa2,Xv2rhoasigmaaa,Xv2sigmaaa2)
else if (iDFTxcsel==3.or.iDFTxcsel==85) then
    call rks_x_pw91(ixcderv,1,rhoa1,sigmaaa1,Xzk,Xvrhoa,Xvsigmaaa,Xv2rhoa2,Xv2rhoasigmaaa,Xv2sigmaaa2)
end if
!C part
if (iDFTxcsel==30.or.iDFTxcsel==80) then
    call rks_c_vwn5(ixcderv,1,rhoa1,sigmaaa1,Czk,Cvrhoa,Cvsigmaaa,Cv2rhoa2,Cv2rhoasigmaaa,Cv2sigmaaa2)
else if (iDFTxcsel==31.or.iDFTxcsel==81) then
    call rks_c_p86(ixcderv,1,rhoa1,sigmaaa1,Czk,Cvrhoa,Cvsigmaaa,Cv2rhoa2,Cv2rhoasigmaaa,Cv2sigmaaa2)
else if (iDFTxcsel==32.or.iDFTxcsel==82) then
    call rks_c_lyp(ixcderv,1,rhoa1,sigmaaa1,Czk,Cvrhoa,Cvsigmaaa,Cv2rhoa2,Cv2rhoasigmaaa,Cv2sigmaaa2)
else if (iDFTxcsel==33.or.iDFTxcsel==83.or.iDFTxcsel==85) then
    call rks_c_pw91(ixcderv,1,rhoa1,sigmaaa1,Czk,Cvrhoa,Cvsigmaaa,Cv2rhoa2,Cv2rhoasigmaaa,Cv2sigmaaa2)
else if (iDFTxcsel==34.or.iDFTxcsel==84) then
    call rks_c_pbe(ixcderv,1,rhoa1,sigmaaa1,Czk,Cvrhoa,Cvsigmaaa,Cv2rhoa2,Cv2rhoasigmaaa,Cv2sigmaaa2)
end if
!Whole XC
if (iDFTxcsel==70) then
    call rks_xc_b97(ixcderv,1,rhoa1,sigmaaa1,XCzk,XCvrhoa,XCvsigmaaa,XCv2rhoa2,XCv2rhoasigmaaa,XCv2sigmaaa2)
else if (iDFTxcsel==71) then
    call rks_xc_hcth407(ixcderv,1,rhoa1,sigmaaa1,XCzk,XCvrhoa,XCvsigmaaa,XCv2rhoa2,XCv2rhoasigmaaa,XCv2sigmaaa2)
else if (iDFTxcsel>=80.and.iDFTxcsel<99) then
    XCzk=Xzk+Czk
    XCvrhoa=Xvrhoa+Cvrhoa
    XCvsigmaaa=Xvsigmaaa+Cvsigmaaa
    XCv2rhoa2=Xv2rhoa2+Cv2rhoa2
    XCv2rhoasigmaaa=Xv2rhoasigmaaa+Cv2rhoasigmaaa
    XCv2sigmaaa2=Xv2sigmaaa2+Cv2sigmaaa2
end if
!Note that the problem of derivative of the DFT repository is revised by dividing a factor, similarly hereinafter
if (iDFTxcsel<30) then
    value=Xzk(1)
    d1rho=Xvrhoa(1)
    d1sig=Xvsigmaaa(1)/4D0
    d2rho=Xv2rhoa2(1)/2D0
    d2rhosig=Xv2rhoasigmaaa(1)/4D0
    d2sig=Xv2sigmaaa2(1)/16D0
else if (iDFTxcsel<70) then
    value=Czk(1)
    d1rho=Cvrhoa(1)
    d1sig=Cvsigmaaa(1)/4D0
    d2rho=Cv2rhoa2(1)/2D0
    d2rhosig=Cv2rhoasigmaaa(1)/4D0
    d2sig=Cv2sigmaaa2(1)/16D0
else if (iDFTxcsel<100) then
    value=XCzk(1)
    d1rho=XCvrhoa(1)
    d1sig=XCvsigmaaa(1)/4D0
    d2rho=XCv2rhoa2(1)/2D0
    d2rhosig=XCv2rhoasigmaaa(1)/4D0
    d2sig=XCv2sigmaaa2(1)/16D0
end if
end subroutine






!!---Open-shell form of DFTxcfunc routine
real*8 function DFTxcfunc_open(x,y,z)
real*8 x,y,z
call gendensgradab(x,y,z,adens,bdens,tdens,agrad,bgrad,tgrad,abgrad)
call getXCdata_open(0,adens,bdens,agrad**2,bgrad**2,abgrad,DFTxcfunc_open,&
sb,sb,sb,sb,sb,sb,sb,sb,sb,sb,sb,sb,sb,sb,sb,sb,sb,sb,sb,sb)
end function

!---Open-shell form of DFTxcpot routine
real*8 function DFTxcpot_open(x,y,z)
real*8 x,y,z
DFTxcpot_open=0D0
end function

!!!For open-shell cases. Input rho and gradrho^2, return the value and its derivative of selected XC (or X/C only) functional
!The global variable "iDFTxcsel" is used to select the XC functional, see manual
!ixcderv=0: only get value, =1: also get 1st derv., =2: also get 2nd derv.
!Input quantities:
!rhoa=rho_alpha   rhob=rho_beta
!sigaa=|gradrho_alpha|^2  sigbb=|gradrho_beta|^2
!sigab=Dot product between gradrho_alpha and gradrho_beta
subroutine getXCdata_open(ixcderv,rhoa,rhob,sigaa,sigbb,sigab,value,d1rhoa,d1rhob,d1sigaa,d1sigbb,d1sigab,d2rhoaa,d2rhobb,d2rhoab,&
d2rhoasigaa,d2rhoasigab,d2rhoasigbb,d2rhobsigbb,d2rhobsigab,d2rhobsigaa,d2sigaaaa,d2sigaaab,d2sigaabb,d2sigabab,d2sigabbb,d2sigbbbb)
integer ixcderv
!Input arguments
real*8 rhoa,rhob,sigaa,sigbb,sigab,value,d1rhoa,d1rhob,d1sigaa,d1sigbb,d1sigab,d2rhoaa,d2rhobb,d2rhoab,&
d2rhoasigaa,d2rhoasigab,d2rhoasigbb,d2rhobsigbb,d2rhobsigab,d2rhobsigaa,d2sigaaaa,d2sigaaab,d2sigaabb,d2sigabab,d2sigabbb,d2sigbbbb
!Inputted information
real*8 rhoa1(1),rhob1(1),sigmaaa1(1),sigmabb1(1),sigmaab1(1)
!Returned information
real*8 Xzk(1),Xvrhoa(1),Xvrhob(1),Xvsigmaaa(1),Xvsigmabb(1),Xvsigmaab(1),Xv2rhoa2(1),Xv2rhob2(1),Xv2rhoab(1)&
,Xv2rhoasigmaaa(1),Xv2rhoasigmaab(1),Xv2rhoasigmabb(1),Xv2rhobsigmabb(1),Xv2rhobsigmaab(1),Xv2rhobsigmaaa(1)&
,Xv2sigmaaa2(1),Xv2sigmaaaab(1),Xv2sigmaaabb(1),Xv2sigmaab2(1),Xv2sigmaabbb(1),Xv2sigmabb2(1)
real*8 Czk(1),Cvrhoa(1),Cvrhob(1),Cvsigmaaa(1),Cvsigmabb(1),Cvsigmaab(1),Cv2rhoa2(1),Cv2rhob2(1),Cv2rhoab(1)&
,Cv2rhoasigmaaa(1),Cv2rhoasigmaab(1),Cv2rhoasigmabb(1),Cv2rhobsigmabb(1),Cv2rhobsigmaab(1),Cv2rhobsigmaaa(1)&
,Cv2sigmaaa2(1),Cv2sigmaaaab(1),Cv2sigmaaabb(1),Cv2sigmaab2(1),Cv2sigmaabbb(1),Cv2sigmabb2(1)
real*8 XCzk(1),XCvrhoa(1),XCvrhob(1),XCvsigmaaa(1),XCvsigmabb(1),XCvsigmaab(1),XCv2rhoa2(1),XCv2rhob2(1),XCv2rhoab(1)&
,XCv2rhoasigmaaa(1),XCv2rhoasigmaab(1),XCv2rhoasigmabb(1),XCv2rhobsigmabb(1),XCv2rhobsigmaab(1),XCv2rhobsigmaaa(1)&
,XCv2sigmaaa2(1),XCv2sigmaaaab(1),XCv2sigmaaabb(1),XCv2sigmaab2(1),XCv2sigmaabbb(1),XCv2sigmabb2(1)

rhoa1(1)=rhoa
rhob1(1)=rhob
sigmaaa1(1)=sigaa
sigmabb1(1)=sigbb
sigmaab1(1)=sigab
!X part
if (iDFTxcsel==0.or.iDFTxcsel==80) then
    call uks_x_lda(ixcderv,1,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,Xzk,Xvrhoa,Xvrhob,Xvsigmaaa,Xvsigmabb,Xvsigmaab,Xv2rhoa2,Xv2rhob2,Xv2rhoab,&
    Xv2rhoasigmaaa,Xv2rhoasigmaab,Xv2rhoasigmabb,Xv2rhobsigmabb,Xv2rhobsigmaab,Xv2rhobsigmaaa,&
    Xv2sigmaaa2,Xv2sigmaaaab,Xv2sigmaaabb,Xv2sigmaab2,Xv2sigmaabbb,Xv2sigmabb2)
else if (iDFTxcsel==1.or.iDFTxcsel==81.or.iDFTxcsel==82.or.iDFTxcsel==83) then
    call uks_x_b88(ixcderv,1,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,Xzk,Xvrhoa,Xvrhob,Xvsigmaaa,Xvsigmabb,Xvsigmaab,Xv2rhoa2,Xv2rhob2,Xv2rhoab,&
    Xv2rhoasigmaaa,Xv2rhoasigmaab,Xv2rhoasigmabb,Xv2rhobsigmabb,Xv2rhobsigmaab,Xv2rhobsigmaaa,&
    Xv2sigmaaa2,Xv2sigmaaaab,Xv2sigmaaabb,Xv2sigmaab2,Xv2sigmaabbb,Xv2sigmabb2)
else if (iDFTxcsel==2.or.iDFTxcsel==84) then
    call uks_x_pbe(ixcderv,1,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,Xzk,Xvrhoa,Xvrhob,Xvsigmaaa,Xvsigmabb,Xvsigmaab,Xv2rhoa2,Xv2rhob2,Xv2rhoab,&
    Xv2rhoasigmaaa,Xv2rhoasigmaab,Xv2rhoasigmabb,Xv2rhobsigmabb,Xv2rhobsigmaab,Xv2rhobsigmaaa,&
    Xv2sigmaaa2,Xv2sigmaaaab,Xv2sigmaaabb,Xv2sigmaab2,Xv2sigmaabbb,Xv2sigmabb2)
else if (iDFTxcsel==3.or.iDFTxcsel==85) then
    call uks_x_pw91(ixcderv,1,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,Xzk,Xvrhoa,Xvrhob,Xvsigmaaa,Xvsigmabb,Xvsigmaab,Xv2rhoa2,Xv2rhob2,Xv2rhoab,&
    Xv2rhoasigmaaa,Xv2rhoasigmaab,Xv2rhoasigmabb,Xv2rhobsigmabb,Xv2rhobsigmaab,Xv2rhobsigmaaa,&
    Xv2sigmaaa2,Xv2sigmaaaab,Xv2sigmaaabb,Xv2sigmaab2,Xv2sigmaabbb,Xv2sigmabb2)
end if
!C part
if (iDFTxcsel==30.or.iDFTxcsel==80) then
    call uks_c_vwn5(ixcderv,1,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,Czk,Cvrhoa,Cvrhob,Cvsigmaaa,Cvsigmabb,Cvsigmaab,Cv2rhoa2,Cv2rhob2,Cv2rhoab,&
    Cv2rhoasigmaaa,Cv2rhoasigmaab,Cv2rhoasigmabb,Cv2rhobsigmabb,Cv2rhobsigmaab,Cv2rhobsigmaaa,&
    Cv2sigmaaa2,Cv2sigmaaaab,Cv2sigmaaabb,Cv2sigmaab2,Cv2sigmaabbb,Cv2sigmabb2)
else if (iDFTxcsel==31.or.iDFTxcsel==81) then
    call uks_c_p86(ixcderv,1,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,Czk,Cvrhoa,Cvrhob,Cvsigmaaa,Cvsigmabb,Cvsigmaab,Cv2rhoa2,Cv2rhob2,Cv2rhoab,&
    Cv2rhoasigmaaa,Cv2rhoasigmaab,Cv2rhoasigmabb,Cv2rhobsigmabb,Cv2rhobsigmaab,Cv2rhobsigmaaa,&
    Cv2sigmaaa2,Cv2sigmaaaab,Cv2sigmaaabb,Cv2sigmaab2,Cv2sigmaabbb,Cv2sigmabb2)
else if (iDFTxcsel==32.or.iDFTxcsel==82) then
    call uks_c_lyp(ixcderv,1,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,Czk,Cvrhoa,Cvrhob,Cvsigmaaa,Cvsigmabb,Cvsigmaab,Cv2rhoa2,Cv2rhob2,Cv2rhoab,&
    Cv2rhoasigmaaa,Cv2rhoasigmaab,Cv2rhoasigmabb,Cv2rhobsigmabb,Cv2rhobsigmaab,Cv2rhobsigmaaa,&
    Cv2sigmaaa2,Cv2sigmaaaab,Cv2sigmaaabb,Cv2sigmaab2,Cv2sigmaabbb,Cv2sigmabb2)
else if (iDFTxcsel==33.or.iDFTxcsel==83.or.iDFTxcsel==85) then
    call uks_c_pw91(ixcderv,1,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,Czk,Cvrhoa,Cvrhob,Cvsigmaaa,Cvsigmabb,Cvsigmaab,Cv2rhoa2,Cv2rhob2,Cv2rhoab,&
    Cv2rhoasigmaaa,Cv2rhoasigmaab,Cv2rhoasigmabb,Cv2rhobsigmabb,Cv2rhobsigmaab,Cv2rhobsigmaaa,&
    Cv2sigmaaa2,Cv2sigmaaaab,Cv2sigmaaabb,Cv2sigmaab2,Cv2sigmaabbb,Cv2sigmabb2)
else if (iDFTxcsel==34.or.iDFTxcsel==84) then
    call uks_c_pbe(ixcderv,1,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,Czk,Cvrhoa,Cvrhob,Cvsigmaaa,Cvsigmabb,Cvsigmaab,Cv2rhoa2,Cv2rhob2,Cv2rhoab,&
    Cv2rhoasigmaaa,Cv2rhoasigmaab,Cv2rhoasigmabb,Cv2rhobsigmabb,Cv2rhobsigmaab,Cv2rhobsigmaaa,&
    Cv2sigmaaa2,Cv2sigmaaaab,Cv2sigmaaabb,Cv2sigmaab2,Cv2sigmaabbb,Cv2sigmabb2)
end if
!Whole XC
if (iDFTxcsel==70) then
    call uks_xc_b97(ixcderv,1,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,XCzk,XCvrhoa,XCvrhob,XCvsigmaaa,XCvsigmabb,XCvsigmaab,XCv2rhoa2,XCv2rhob2,XCv2rhoab,&
    XCv2rhoasigmaaa,XCv2rhoasigmaab,XCv2rhoasigmabb,XCv2rhobsigmabb,XCv2rhobsigmaab,XCv2rhobsigmaaa,&
    XCv2sigmaaa2,XCv2sigmaaaab,XCv2sigmaaabb,XCv2sigmaab2,XCv2sigmaabbb,XCv2sigmabb2)
else if (iDFTxcsel==71) then
    call uks_xc_hcth407(ixcderv,1,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,XCzk,XCvrhoa,XCvrhob,XCvsigmaaa,XCvsigmabb,XCvsigmaab,XCv2rhoa2,XCv2rhob2,XCv2rhoab,&
    XCv2rhoasigmaaa,XCv2rhoasigmaab,XCv2rhoasigmabb,XCv2rhobsigmabb,XCv2rhobsigmaab,XCv2rhobsigmaaa,&
    XCv2sigmaaa2,XCv2sigmaaaab,XCv2sigmaaabb,XCv2sigmaab2,XCv2sigmaabbb,XCv2sigmabb2)
else if (iDFTxcsel>=80.and.iDFTxcsel<99) then
    XCzk=Xzk+Czk
    XCvrhoa=Xvrhoa+Cvrhoa
    XCvrhob=Xvrhob+Cvrhob
    XCvsigmaaa=Xvsigmaaa+Cvsigmaaa
    XCvsigmabb=Xvsigmabb+Cvsigmabb
    XCvsigmaab=Xvsigmaab+Cvsigmaab
    XCv2rhoa2=Xv2rhoa2+Cv2rhoa2
    XCv2rhob2=Xv2rhob2+Cv2rhob2
    XCv2rhoab=Xv2rhoab+Cv2rhoab
    XCv2rhoasigmaaa=Xv2rhoasigmaaa+Cv2rhoasigmaaa
    XCv2rhoasigmaab=Xv2rhoasigmaab+Cv2rhoasigmaab
    XCv2rhoasigmabb=Xv2rhoasigmabb+Cv2rhoasigmabb
    XCv2rhobsigmabb=Xv2rhobsigmabb+Cv2rhobsigmabb
    XCv2rhobsigmaab=Xv2rhobsigmaab+Cv2rhobsigmaab
    XCv2rhobsigmaaa=Xv2rhobsigmaaa+Cv2rhobsigmaaa
    XCv2sigmaaa2= Xv2sigmaaa2+Cv2sigmaaa2
    XCv2sigmaaaab=Xv2sigmaaaab+Cv2sigmaaaab
    XCv2sigmaaabb=Xv2sigmaaabb+Cv2sigmaaabb
    XCv2sigmaab2= Xv2sigmaab2+Cv2sigmaab2
    XCv2sigmaabbb=Xv2sigmaabbb+Cv2sigmaabbb
    XCv2sigmabb2= Xv2sigmabb2+Cv2sigmabb2
end if

if (iDFTxcsel<30) then
    value=Xzk(1)
    d1rhoa=Xvrhoa(1)
    d1rhob=Xvrhob(1)
    d1sigaa=Xvsigmaaa(1)
    d1sigbb=Xvsigmabb(1)
    d1sigab=Xvsigmaab(1)
    d2rhoaa=Xv2rhoa2(1)
    d2rhobb=Xv2rhob2(1)
    d2rhoab=Xv2rhoab(1)
    d2rhoasigaa=Xv2rhoasigmaaa(1)
    d2rhoasigab=Xv2rhoasigmaab(1)
    d2rhoasigbb=Xv2rhoasigmabb(1)
    d2rhobsigbb=Xv2rhobsigmabb(1)
    d2rhobsigab=Xv2rhobsigmaab(1)
    d2rhobsigaa=Xv2rhobsigmaaa(1)
    d2sigaaaa=Xv2sigmaaa2(1)
    d2sigaaab=Xv2sigmaaaab(1)
    d2sigaabb=Xv2sigmaaabb(1)
    d2sigabab=Xv2sigmaab2(1)
    d2sigabbb=Xv2sigmaabbb(1)
    d2sigbbbb=Xv2sigmabb2(1)
else if (iDFTxcsel<70) then
    value=Czk(1)
    d1rhoa=Cvrhoa(1)
    d1rhob=Cvrhob(1)
    d1sigaa=Cvsigmaaa(1)
    d1sigbb=Cvsigmabb(1)
    d1sigab=Cvsigmaab(1)
    d2rhoaa=Cv2rhoa2(1)
    d2rhobb=Cv2rhob2(1)
    d2rhoab=Cv2rhoab(1)
    d2rhoasigaa=Cv2rhoasigmaaa(1)
    d2rhoasigab=Cv2rhoasigmaab(1)
    d2rhoasigbb=Cv2rhoasigmabb(1)
    d2rhobsigbb=Cv2rhobsigmabb(1)
    d2rhobsigab=Cv2rhobsigmaab(1)
    d2rhobsigaa=Cv2rhobsigmaaa(1)
    d2sigaaaa=Cv2sigmaaa2(1)
    d2sigaaab=Cv2sigmaaaab(1)
    d2sigaabb=Cv2sigmaaabb(1)
    d2sigabab=Cv2sigmaab2(1)
    d2sigabbb=Cv2sigmaabbb(1)
    d2sigbbbb=Cv2sigmabb2(1)
else if (iDFTxcsel<100) then
    value=XCzk(1)
    d1rhoa=XCvrhoa(1)
    d1rhob=XCvrhob(1)
    d1sigaa=XCvsigmaaa(1)
    d1sigbb=XCvsigmabb(1)
    d1sigab=XCvsigmaab(1)
    d2rhoaa=XCv2rhoa2(1)
    d2rhobb=XCv2rhob2(1)
    d2rhoab=XCv2rhoab(1)
    d2rhoasigaa=XCv2rhoasigmaaa(1)
    d2rhoasigab=XCv2rhoasigmaab(1)
    d2rhoasigbb=XCv2rhoasigmabb(1)
    d2rhobsigbb=XCv2rhobsigmabb(1)
    d2rhobsigab=XCv2rhobsigmaab(1)
    d2rhobsigaa=XCv2rhobsigmaaa(1)
    d2sigaaaa=XCv2sigmaaa2(1)
    d2sigaaab=XCv2sigmaaaab(1)
    d2sigaabb=XCv2sigmaaabb(1)
    d2sigabab=XCv2sigmaab2(1)
    d2sigabbb=XCv2sigmaabbb(1)
    d2sigbbbb=XCv2sigmabb2(1)
end if
end subroutine



!!---- The distance from a point (x,y,z) to the nearest atom in the array
real*8 function surfana_di(x,y,z,nlen,atmlist)
real*8 x,y,z
integer nlen,atmlist(nlen)
dist2min=1D100
do iatm=1,ncenter
    if (any(atmlist==iatm)) then !The atom is in the list
        dist2=(a(iatm)%x-x)**2+(a(iatm)%y-y)**2+(a(iatm)%z-z)**2
        if (dist2<dist2min) dist2min=dist2
    end if
end do
surfana_di=dsqrt(dist2min)
end function

!!---- The distance from a point (x,y,z) to the nearest atom not in the array
real*8 function surfana_de(x,y,z,nlen,atmlist)
real*8 x,y,z
integer nlen,atmlist(nlen)
dist2min=1D100
do iatm=1,ncenter
    if (all(atmlist/=iatm)) then !The atom is not in the list
        dist2=(a(iatm)%x-x)**2+(a(iatm)%y-y)**2+(a(iatm)%z-z)**2
        if (dist2<dist2min) dist2min=dist2
    end if
end do
surfana_de=dsqrt(dist2min)
end function

!!---- Normalized contact distance, defined in terms of de, di and the vdW radii of the atoms
real*8 function surfana_norm(x,y,z,nlen,atmlist)
real*8 x,y,z
integer nlen,atmlist(nlen)
dist2minin=1D100 !The nearest distance to atoms inside
dist2minext=1D100 !The nearest distance to atoms outside
do iatm=1,ncenter
    dist2=(a(iatm)%x-x)**2+(a(iatm)%y-y)**2+(a(iatm)%z-z)**2
    if (any(atmlist==iatm)) then !Atoms inside
        if (dist2<dist2minin) then
            dist2minin=dist2
            iminin=iatm
        end if
    else !Atoms outside
        if (dist2<dist2minext) then
            dist2minext=dist2
            iminext=iatm
        end if
    end if
end do
di=dsqrt(dist2minin)
de=dsqrt(dist2minext)
rvdwin=vdwr(a(iminin)%index)
rvdwext=vdwr(a(iminext)%index)
surfana_norm=(di-rvdwin)/rvdwin+(de-rvdwext)/rvdwext
end function




!!-------- PAEM, potential acting on one electron in a molecule, defined by Zhongzhi Yang in JCC,35,965(2014)
!If itype=1, evaluate the XC potential based on pair density; if =2, it will be equivalent to DFT XC potential
real*8 function PAEM(x,y,z,itype)
integer itype
real*8 x,y,z,wfnval(nmo),GTFint(nprims,nprims)
!Evaluate electron contribution to ESP
call genGTFattmat(x,y,z,GTFint)
rhopot=0
do imo=1,nmo
    do iprim=1,nprims
        do jprim=1,nprims
            rhopot=rhopot+MOocc(imo)*CO(imo,iprim)*CO(imo,jprim)*GTFint(iprim,jprim)
        end do
    end do
end do
!Evaluate XC potential
if (itype==1) then !Based on Muller approximation form
    call orbderv(1,1,nmo,x,y,z,wfnval)
    rho=sum(MOocc(1:nmo)*wfnval(1:nmo)**2)
    xcpot=0
    if (wfntype==0.or.wfntype==3) then !Close-shell
        do imo=1,nmo
            if (MOocc(imo)==0D0) cycle
            do jmo=1,nmo
                if (MOocc(jmo)==0D0) cycle
                tmpval=dsqrt(MOocc(imo)*MOocc(jmo))*wfnval(imo)*wfnval(jmo)
                do iprim=1,nprims
                    do jprim=1,nprims
                        xcpot=xcpot-tmpval*CO(imo,iprim)*CO(jmo,jprim)*GTFint(iprim,jprim)
                    end do
                end do
            end do
        end do
    else if (wfntype==1.or.wfntype==4) then !Unrestricted open-shell
        do ialphaend=nmo,1,-1 !Find the ending index of alpha MO
            if (MOtype(ialphaend)==1) exit
        end do
        do imo=1,ialphaend !Alpha part
            do jmo=1,ialphaend
                tmpval=dsqrt(MOocc(imo)*MOocc(jmo))*wfnval(imo)*wfnval(jmo)
                do iprim=1,nprims
                    do jprim=1,nprims
                        xcpot=xcpot-tmpval*CO(imo,iprim)*CO(jmo,jprim)*GTFint(iprim,jprim)
                    end do
                end do
            end do
        end do
        do imo=ialphaend+1,nmo !Beta part
            do jmo=ialphaend+1,nmo
                tmpval=dsqrt(MOocc(imo)*MOocc(jmo))*wfnval(imo)*wfnval(jmo)
                do iprim=1,nprims
                    do jprim=1,nprims
                        xcpot=xcpot-tmpval*CO(imo,iprim)*CO(jmo,jprim)*GTFint(iprim,jprim)
                    end do
                end do
            end do
        end do
    else if (wfntype==2) then !Restricted open-shell
        do imo=1,nmo !Alpha part
            if (MOocc(imo)==0) cycle
            do jmo=1,nmo
                if (MOocc(jmo)==0) cycle
                tmpval=wfnval(imo)*wfnval(jmo) !Every occupied ROHF MOs contributes one alpha electron
                do iprim=1,nprims
                    do jprim=1,nprims
                        xcpot=xcpot-tmpval*CO(imo,iprim)*CO(jmo,jprim)*GTFint(iprim,jprim)
                    end do
                end do
            end do
        end do
        do imo=1,nmo !Beta part
            if (MOocc(imo)/=2D0) cycle
            do jmo=1,nmo
                if (MOocc(jmo)/=2D0) cycle
                tmpval=wfnval(imo)*wfnval(jmo) !Every doubly occupied ROHF MOs contributes one beta electron
                do iprim=1,nprims
                    do jprim=1,nprims
                        xcpot=xcpot-tmpval*CO(imo,iprim)*CO(jmo,jprim)*GTFint(iprim,jprim)
                    end do
                end do
            end do
        end do
    end if
    xcpot=xcpot/rho
    
else if (itype==2) then !Directly using DFT XC potential
    xcpot=DFTxcpot(x,y,z)
end if

PAEM=-nucesp(x,y,z)+rhopot+xcpot
end function



!!----- Angle between the eigenvectors of rho and the plane defined by option 4 of main function 1000
! The plane is represented by global variables pleA,pleB,pleC,pleD
! ivec=1/2/3 means the eigenvector corresponding to the first/second/third highest eigenvalue is calculated
real*8 function Ang_rhoeigvec_ple(x,y,z,ivec)
use util
integer ivec
real*8 x,y,z,eigvecmat(3,3),eigval(3),eigvaltmp(3),elehess(3,3),elegrad(3)
call calchessmat_dens(2,x,y,z,elerho,elegrad,elehess)
call diagmat(elehess,eigvecmat,eigval,100,1D-10)
eigvaltmp=eigval
call sort(eigvaltmp) !From small to large
call invarr(eigvaltmp,1,3) !1/2/3=large to small
do i=1,3
    if (eigval(i)==eigvaltmp(ivec)) exit
end do
Ang_rhoeigvec_ple=vecang(eigvecmat(1,i),eigvecmat(2,i),eigvecmat(3,i),pleA,pleB,pleC)
end function


!!------ Local electron correlation function (DOI: 10.1021/acs.jctc.7b00293)
!itype=1: Local total electron correlation function
!itype=2: Local dynamic electron correlation function
!itype=3: Local nondynamic electron correlation function
real*8 function localcorr(x,y,z,itype)
integer itype
real*8 x,y,z,wfnval(nmo),occ(nmo)
call orbderv(1,1,nmo,x,y,z,wfnval)
localcorr=0D0
if (wfntype==3) then
    occ=MOocc/2
    where(occ>1) occ=1 !Remove unphysical larger than unity occupation number
    where(occ<0) occ=0 !Remove unphysical negative occupation number
    if (itype==1) then
        do i=1,nmo
            localcorr=localcorr+ dsqrt(occ(i)*(1-occ(i)))*wfnval(i)**2
        end do
        localcorr=localcorr/4
    else if (itype==2) then
        do i=1,nmo
            localcorr=localcorr+ ( dsqrt(occ(i)*(1-occ(i))) - 2*occ(i)*(1-occ(i)) ) *wfnval(i)**2
        end do
        localcorr=localcorr/4
    else if (itype==3) then
        do i=1,nmo
            localcorr=localcorr+ occ(i)*(1-occ(i))*wfnval(i)**2
        end do
        localcorr=localcorr/2
    end if
    localcorr=localcorr*2
else if (wfntype==4) then
    occ=MOocc
    where(occ>1) occ=1
    where(occ<0) occ=0
    if (itype==1) then
        do i=1,nmo
            localcorr=localcorr+ dsqrt(occ(i)*(1-occ(i)))*wfnval(i)**2
        end do
        localcorr=localcorr/4
    else if (itype==2) then
        do i=1,nmo
            localcorr=localcorr+ ( dsqrt(occ(i)*(1-occ(i))) - 2*occ(i)*(1-occ(i)) ) *wfnval(i)**2
        end do
        localcorr=localcorr/4
    else if (itype==3) then
        do i=1,nmo
            localcorr=localcorr+ occ(i)*(1-occ(i))*wfnval(i)**2
        end do
        localcorr=localcorr/2
    end if
end if
end function


!For visually examine functions used in DFRT2.0 project
real*8 function funcvalLSB(x,y,z,itype)
integer itype
real*8 x,y,z,valarr(6),rho,gradrho(3)
iexpcutoffold=expcutoff
expcutoff=1
call valaryyLSB(x,y,z,valarr,rho,gradrho)
funcvalLSB=valarr(itype)
expcutoff=iexpcutoffold
end function



!*********************************************************************
!=====================================================================
!The real space functions below are used for shubins' DFRT 2.0 project
!=====================================================================
!*********************************************************************






end module
