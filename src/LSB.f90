!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!==========================================================================
!**************************************************************************
!  The codes in this file are written specific for Shubin Liu's project
!**************************************************************************
!==========================================================================
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fuzzySBL
use defvar
use util
use function
implicit real*8 (a-h,o-z)
integer,parameter :: nfunc=7,nquant=7 !The number of real space function, the number of quantities to be calculated
real*8 smat(ncenter,ncenter),Pvec(ncenter),atmBeckewei(radpot*sphpot),atmHirshwei(radpot*sphpot,ncenter)
real*8 atmcontri(ncenter,0:nfunc,0:nquant)
real*8 promol(radpot*sphpot),atomdens(radpot*sphpot,ncenter),selfdens(radpot*sphpot)
real*8 promolgrad(radpot*sphpot,3),atomgrad(radpot*sphpot,3)
real*8 funcval(radpot*sphpot,0:nfunc),funcgrdn(radpot*sphpot,0:nfunc) !Store function value and gradient norm, respectively
real*8 funcref(radpot*sphpot,0:nfunc) !Store function value of promolecular state (reference state)
real*8 covr_becke(0:nelesupp) !covalent radii used for Becke partition
type(content) gridatm(radpot*sphpot)
real*8 potx(sphpot),poty(sphpot),potz(sphpot),potw(sphpot)
real*8 arrtmp(nfunc),arrtmp2(nfunc),gradrho(3)
character*40 functionname(0:nfunc),quantityname(0:nquant)

if (ifiletype/=2.and.ifiletype/=3) then
    write(*,*) "Note: Using .wfn or .wfx file will make calculation much faster"
end if

nbeckeiter=3
expcutoff=1
covr_becke=covr_tianlu

functionname(0)=" rho"
functionname(1)=" rho/rho0"
functionname(2)=" |der_rho|/rho^(4/3)"
functionname(3)=" der2rho/rho^(5/3)"
functionname(4)=" (tau-t_w)/t_TF"
functionname(5)=" Xi part of SEDD"
functionname(6)=" Theta part of DORI"
functionname(7)=" Spin density"
quantityname(0)=" Function itself" !Special
quantityname(1)=" Shannon entropy"
quantityname(2)=" Fisher information"
quantityname(3)=" Onicescu information energy of order 2"
quantityname(4)=" Onicescu information energy of order 3"
quantityname(5)=" Information gain"
quantityname(6)=" Relative Renyi entropy of orders 2"
quantityname(7)=" Relative Renyi entropy of orders 3"

write(*,*) "Select partition method:"
write(*,"(a)") " 1 Becke partition"
write(*,"(a)") " 2 Hirshfeld partition"
read(*,*) ipartition

atmcontri=0D0
checkacc=0
call setpromol
write(*,"(' Radial points:',i5,'    Angular points:',i5,'   Total:',i10,' per center')") radpot,sphpot,radpot*sphpot
write(*,*) "Please wait..."
write(*,*)
call walltime(nwalltime1)
call Lebedevgen(sphpot,potx,poty,potz,potw)

! call valaryyLSB(1D0,1D0,1D0,arrtmp,rho,rhogrdn)
! write(*,*) arrtmp(4)
! read(*,*)

do iatm=1,ncenter !! Cycle each atom
    write(*,"(/,' Calculating:',i6,'   /',i6)") iatm,ncenter

    !! Prepare grid points on current center
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
    gridatm%x=gridatm%x+a(iatm)%x !Move quadrature point to actual position in molecule
    gridatm%y=gridatm%y+a(iatm)%y
    gridatm%z=gridatm%z+a(iatm)%z
    
    
    !! Calculate function values of reference state, including promolecular density and density of present atom, function index:
    !index=0: rho
    !index=1: rho/rho0, where rho0 is promolecular density
    !index=2: |der_rho|/rho^(4/3)
    !index=3: der2rho/rho^(5/3)
    !index=4: (tau-t_w)/t_TF
    !index=5: Part of SEDD
    !index=6: Part of DORI
    !index=7: Spin density
    funcref=0D0
    promol=0D0
    promolgrad=0D0
    do jatm=1,ncenter_org
        call dealloall
        call readwfn(custommapname(jatm),1)
nthreads=getNThreads()
!$OMP parallel do shared(atomdens,atomgrad,funcref) private(ipt,rnowx,rnowy,rnowz,arrtmp,rho,gradrho) num_threads(nthreads)
        do ipt=1+iradcut*sphpot,radpot*sphpot
            rnowx=gridatm(ipt)%x
            rnowy=gridatm(ipt)%y
            rnowz=gridatm(ipt)%z
            call valaryyLSB(rnowx,rnowy,rnowz,arrtmp,rho,gradrho)
            atomdens(ipt,jatm)=rho
            atomgrad(ipt,:)=gradrho
            funcref(ipt,0)=funcref(ipt,0)+rho
            funcref(ipt,1)=funcref(ipt,1)+1D0 !i.e. rho/rho0, the value is unity if present system is a single atom
            funcref(ipt,2:nfunc)=funcref(ipt,2:nfunc)+arrtmp(2:nfunc)
        end do
!$OMP end parallel do
        promol=promol+atomdens(:,jatm)
        promolgrad=promolgrad+atomgrad(:,:)
        if (jatm==iatm) selfdens=atomdens(:,jatm)
    end do
    call dealloall
    call readinfile(firstfilename,1) !Retrieve the firstly loaded file(whole molecule)
    
    
    !! Calculate function values and gradient norm for present molecule, then store them in funcval/funcgrdn. Function index:
    !index=0: rho
    !index=1: rho/rho0, where rho0 is promolecular density
    !index=2: |der_rho|/rho^(4/3)
    !index=3: der2rho/rho^(5/3)
    !index=4: (tau-t_w)/t_TF
    !index=5: Part of SEDD
    !index=6: Part of DORI
    !index=7: Spin density
nthreads=getNThreads()
!$OMP parallel do shared(funcval,funcgrdn) private(ipt,arrtmp,arrtmp2,rho,gradrho,idir) num_threads(nthreads)
    do ipt=1+iradcut*sphpot,radpot*sphpot
        call valgradarrLSB(gridatm(ipt)%x,gridatm(ipt)%y,gridatm(ipt)%z,arrtmp,arrtmp2,rho,gradrho)
        funcval(ipt,0)=rho
        funcgrdn(ipt,0)=dsqrt(sum(gradrho**2))
        !Note: see word document for explicit expression of |der(rho/rho0)|
        funcval(ipt,1)=rho/promol(ipt)
        funcgrdn(ipt,1)=0D0
        do idir=1,3
            funcgrdn(ipt,1)=funcgrdn(ipt,1)+( gradrho(idir)/promol(ipt)-rho/promol(ipt)**2*promolgrad(ipt,idir) )**2
        end do
        funcgrdn(ipt,1)=dsqrt(funcgrdn(ipt,1))
        funcval(ipt,2:nfunc)=arrtmp(2:nfunc)
        funcgrdn(ipt,2:nfunc)=arrtmp2(2:nfunc)
    end do
!$OMP end parallel do
    
    
    !! Calculate atomic Becke space weight for present atom at all of its points
!$OMP parallel do shared(atmpartwei) private(i,rnowx,rnowy,rnowz,smat,&
!$OMP ii,ri,jj,rj,rmiu,chi,uij,aij,tmps,iter,Pvec) num_threads(nthreads) schedule(dynamic)
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
        atmBeckewei(i)=Pvec(iatm) !Normalized Pvec, Pvec contain partition weight of each atom in current point, namely i
    end do
!$OMP end parallel do
    
    
    !Calculate Hirshfeld weights of all atoms
    if (ipartition==2) then
        do jatm=1,ncenter
            do i=1+iradcut*sphpot,radpot*sphpot
                if (promol(i)/=0D0) then
                    atmHirshwei(i,jatm)=atomdens(i,jatm)/promol(i)
                else
                    atmHirshwei(i,jatm)=0D0
                end if
            end do
        end do
    end if
    
    
    !! Calculate all quantities contributed by this atom
    if (ipartition==1) then !Atomic grid integration for Becke partition
        do ipt=1+iradcut*sphpot,radpot*sphpot
            weitot=atmBeckewei(ipt)*gridatm(ipt)%value
            do ifunc=0,nfunc
                !Value of function itself
                atmcontri(iatm,ifunc,0)=atmcontri(iatm,ifunc,0)+weitot*funcval(ipt,ifunc)
                if (ifunc/=3) then !Because the variable in log must be positive everywhere
                    !Shannon entropy
                    atmcontri(iatm,ifunc,1)=atmcontri(iatm,ifunc,1)+weitot*( -funcval(ipt,ifunc)*log(funcval(ipt,ifunc)) )
                    !Information gain
                    atmcontri(iatm,ifunc,5)=atmcontri(iatm,ifunc,5)+weitot*( funcval(ipt,ifunc)*log(funcval(ipt,ifunc)/funcref(ipt,ifunc)) )
                end if
                !Fisher information
                atmcontri(iatm,ifunc,2)=atmcontri(iatm,ifunc,2)+weitot*( funcgrdn(ipt,ifunc)**2 /funcval(ipt,ifunc) )
                !Onicescu information energy of order 2
                atmcontri(iatm,ifunc,3)=atmcontri(iatm,ifunc,3)+weitot*( funcval(ipt,ifunc)**2 )
                !Onicescu information energy of order 3
                atmcontri(iatm,ifunc,4)=atmcontri(iatm,ifunc,4)+weitot*( funcval(ipt,ifunc)**3 )
                !Relative Renyi entropy of orders 2
                atmcontri(iatm,ifunc,6)=atmcontri(iatm,ifunc,6)+weitot*( funcval(ipt,ifunc)**2/funcref(ipt,ifunc) )
                !Relative Renyi entropy of orders 3
                atmcontri(iatm,ifunc,7)=atmcontri(iatm,ifunc,7)+weitot*( funcval(ipt,ifunc)**3/funcref(ipt,ifunc)**2 )
            end do
        end do
        
        write(*,"(' ========== Contribution of atom',i4,a,':')") iatm,a(iatm)%name
        do ifunc=0,nfunc
            write(*,"(/,' Function:',a)") trim(functionname(ifunc))
            do iquant=0,nquant !0 corresponds to function itself
                if (ifunc==3.and.(iquant==1.or.iquant==5)) cycle
                write(*,"(a,':',1E16.8)") quantityname(iquant),atmcontri(iatm,ifunc,iquant)
            end do
        end do
        
    else if (ipartition==2) then !Molecular grid integration for Hirshfeld partition
        do jatm=1,ncenter
            do ipt=1+iradcut*sphpot,radpot*sphpot
                weitot=atmBeckewei(ipt)*atmHirshwei(ipt,jatm)*gridatm(ipt)%value
                do ifunc=0,nfunc
                    !Value of function itself
                    atmcontri(jatm,ifunc,0)=atmcontri(jatm,ifunc,0)+weitot*funcval(ipt,ifunc)
                    if (ifunc/=3) then !Because the variable in log must be positive everywhere
                        !Shannon entropy
                        atmcontri(jatm,ifunc,1)=atmcontri(jatm,ifunc,1)+weitot*( -funcval(ipt,ifunc)*log(funcval(ipt,ifunc)) )
                        !Information gain
                        atmcontri(jatm,ifunc,5)=atmcontri(jatm,ifunc,5)+weitot*( funcval(ipt,ifunc)*log(funcval(ipt,ifunc)/funcref(ipt,ifunc)) )
                    end if
                    !Fisher information
                    atmcontri(jatm,ifunc,2)=atmcontri(jatm,ifunc,2)+weitot*( funcgrdn(ipt,ifunc)**2 /funcval(ipt,ifunc) )
                    !Onicescu information energy of order 2
                    atmcontri(jatm,ifunc,3)=atmcontri(jatm,ifunc,3)+weitot*( funcval(ipt,ifunc)**2 )
                    !Onicescu information energy of order 3
                    atmcontri(jatm,ifunc,4)=atmcontri(jatm,ifunc,4)+weitot*( funcval(ipt,ifunc)**3 )
                    !Relative Renyi entropy of orders 2
                    atmcontri(jatm,ifunc,6)=atmcontri(jatm,ifunc,6)+weitot*( funcval(ipt,ifunc)**2/funcref(ipt,ifunc) )
                    !Relative Renyi entropy of orders 3
                    atmcontri(jatm,ifunc,7)=atmcontri(jatm,ifunc,7)+weitot*( funcval(ipt,ifunc)**3/funcref(ipt,ifunc)**2 )
                end do
            end do
        end do
    end if
    
end do !End cycling atoms

!Output atomic contribution obtained from molecular grid integration
if (ipartition==2) then
    do iatm=1,ncenter
        write(*,"(/,' ========== Contribution of atom',i4,a,':')") iatm,a(iatm)%name
        do ifunc=0,nfunc
            write(*,"(/,' Function:',a)") trim(functionname(ifunc))
            do iquant=0,nquant !0 corresponds to function itself
                if (ifunc==3.and.(iquant==1.or.iquant==5)) cycle
                write(*,"(a,':',1E16.8)") quantityname(iquant),atmcontri(iatm,ifunc,iquant)
            end do
        end do
    end do
end if

write(*,*)
write(*,*) "========================================="
write(*,*) "             Overall result"
write(*,*) "========================================="
do ifunc=0,nfunc
    write(*,"(/,' Function:',a)") trim(functionname(ifunc))
    do iquant=0,nquant
        if (ifunc==3.and.(iquant==1.or.iquant==5)) cycle
        write(*,"(a,':',1E16.8)") quantityname(iquant),sum(atmcontri(:,ifunc,iquant))
    end do
end do

call walltime(nwalltime2)
write(*,"(/,' Calculation took up',i8,' seconds wall clock time')") nwalltime2-nwalltime1

end subroutine



!Get value for all the functions used in LSB's DFRT 2.0 project
!Slot of "valarr" array is shown below. rho and gradient of rho are individually returned
!slot=2: |der_rho|/rho^(4/3)
!slot=3: der2rho/rho^(5/3)
!slot=4: (tau-t_w)/t_TF
!slot=5: Xi part of SEDD
!slot=6: Theta part of DORI
!slot=7: Spin density
subroutine valaryyLSB(x,y,z,valarr,rho,gradrho)
use defvar
use function
implicit real*8 (a-h,o-z)
real*8 x,y,z,wfnval(nmo),wfnderv(3,nmo),wfnhess(3,3,nmo),gradrho(3),hess(3,3),valarr(6)
real*8 gradrhoa(3),gradrhob(3),MOoccnow
real*8 :: Fc=2.871234000D0 ! Fermi constant = (3/10)*(3*Pi^2)**(2/3) = 2.871234, 1/2.871234=0.34828
real*8 :: Fc_pol=4.557799872D0 ! Fermi constant for spin polarized = (3/10)*(6*Pi^2)**(2/3) = 4.5578, 1/4.5578=0.2194

call orbderv(4,1,nmo,x,y,z,wfnval,wfnderv,wfnhess)
rho=0D0
gradrho=0D0
do i=1,nmo
    rho=rho+MOocc(i)*wfnval(i)**2
    gradrho(:)=gradrho(:)+MOocc(i)*wfnval(i)*wfnderv(:,i)
end do

gradrho=2*gradrho
rhogrdn=dsqrt(sum(gradrho**2))
dersqr=sum(gradrho**2)
hess(1,1)=2*sum( MOocc(1:nmo)*( wfnderv(1,1:nmo)**2 + wfnval(1:nmo)*wfnhess(1,1,1:nmo) ) )
hess(2,2)=2*sum( MOocc(1:nmo)*( wfnderv(2,1:nmo)**2 + wfnval(1:nmo)*wfnhess(2,2,1:nmo) ) )
hess(3,3)=2*sum( MOocc(1:nmo)*( wfnderv(3,1:nmo)**2 + wfnval(1:nmo)*wfnhess(3,3,1:nmo) ) )
hess(1,2)=2*sum( MOocc(1:nmo)*( wfnderv(1,1:nmo)*wfnderv(2,1:nmo)+wfnhess(1,2,1:nmo)*wfnval(1:nmo) ) )
hess(2,3)=2*sum( MOocc(1:nmo)*( wfnderv(2,1:nmo)*wfnderv(3,1:nmo)+wfnhess(2,3,1:nmo)*wfnval(1:nmo) ) )
hess(1,3)=2*sum( MOocc(1:nmo)*( wfnderv(1,1:nmo)*wfnderv(3,1:nmo)+wfnhess(1,3,1:nmo)*wfnval(1:nmo) ) )
hess(2,1)=hess(1,2)
hess(3,2)=hess(2,3)
hess(3,1)=hess(1,3)
der2rho=hess(1,1)+hess(2,2)+hess(3,3)

!|der_rho|/rho^(4/3)
valarr(2)=rhogrdn/rho**(4D0/3D0)

!der2rho/rho^(5/3)
valarr(3)=der2rho/rho**(5D0/3D0)

!(tau-t_w)/t_TF
D=0D0
rhoa=0D0
rhob=0D0
gradrhoa=0D0
gradrhob=0D0
if (wfntype==0.or.wfntype==3) then !spin-unpolarized case
    do i=1,nmo
        D=D+MOocc(i)*(sum(wfnderv(:,i)**2)) !Calculate actual kinetic term
    end do        
    Dh=Fc*rho**(5.0D0/3.0D0) !Thomas-Fermi uniform electron gas kinetic energy
    D=D/2.0D0
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
    D=D/2.0D0
    if (rhoa/=0D0) D=D-(sum(gradrhoa(:)**2))/rhoa/8
    if (rhob/=0D0) D=D-(sum(gradrhob(:)**2))/rhob/8
end if
valarr(4)=D/Dh
    
!Xi part of SEDD
tmp1_1=rho*(gradrho(1)*hess(1,1)+gradrho(2)*hess(1,2)+gradrho(3)*hess(1,3))
tmp1_2=gradrho(1)*dersqr
tmp2_1=rho*(gradrho(1)*hess(1,2)+gradrho(2)*hess(2,2)+gradrho(3)*hess(2,3))
tmp2_2=gradrho(2)*dersqr
tmp3_1=rho*(gradrho(1)*hess(1,3)+gradrho(2)*hess(2,3)+gradrho(3)*hess(3,3))
tmp3_2=gradrho(3)*dersqr
valarr(5)=4/rho**8*( (tmp1_1-tmp1_2)**2 + (tmp2_1-tmp2_2)**2 + (tmp3_1-tmp3_2)**2 )

!Theta part of DORI
valarr(6)=4/dersqr**3*( (tmp1_1-tmp1_2)**2 + (tmp2_1-tmp2_2)**2 + (tmp3_1-tmp3_2)**2 )

!Spin density
rhoa=0.0D0
rhob=0.0D0
do i=1,nmo
    if (MOtype(i)==1) then
        rhoa=rhoa+MOocc(i)*wfnval(i)**2
    else if (MOtype(i)==2) then
        rhob=rhob+MOocc(i)*wfnval(i)**2
    else if (MOtype(i)==0) then
        rhoa=rhoa+MOocc(i)/2D0*wfnval(i)**2
        rhob=rhob+MOocc(i)/2D0*wfnval(i)**2
    end if
end do
valarr(7)=rhoa-rhob

end subroutine


!Get value and norm of numerical gradient for all the functions used in LSB's DFRT 2.0 project
!The definition of slot is identical to valaryyLSB
subroutine valgradarrLSB(x,y,z,valarr,grdnarr,rho,gradrho)
use defvar
implicit real*8 (a-h,o-z)
integer,parameter :: nfunc=7
real*8 x,y,z,rho,gradrho(3),valarr(nfunc),grdnarr(nfunc),tmparr(3)
real*8 xadd(nfunc),xmin(nfunc),yadd(nfunc),ymin(nfunc),zadd(nfunc),zmin(nfunc),grdxarr(nfunc),grdyarr(nfunc),grdzarr(nfunc)
diff=5D-4
denom=2D0*diff
call valaryyLSB(x,y,z,valarr,rho,gradrho)
call valaryyLSB(x+diff,y,z,xadd,tmp,tmparr)
call valaryyLSB(x-diff,y,z,xmin,tmp,tmparr)
call valaryyLSB(x,y+diff,z,yadd,tmp,tmparr)
call valaryyLSB(x,y-diff,z,ymin,tmp,tmparr)
call valaryyLSB(x,y,z+diff,zadd,tmp,tmparr)
call valaryyLSB(x,y,z-diff,zmin,tmp,tmparr)
grdxarr=(xadd-xmin)/denom
grdyarr=(yadd-ymin)/denom
grdzarr=(zadd-zmin)/denom
do ifunc=1,nfunc
    grdnarr(ifunc)=dsqrt(grdxarr(ifunc)**2+grdyarr(ifunc)**2+grdzarr(ifunc)**2)
end do
end subroutine




!!------ Integrate real space function in Hirshfeld space with molecular grid (i.e. the grid is the same as integating over the whole space)
!! The grid used in function 1 of fuzzy space analysis is only the grid centered at the atom to be studied, this lead to inaccurate integration result
!! when the integrand varies fast at the tail region of the Hirshfeld atom (commonly close to the other nuclei). In this case we must use the molecular grid.
!! Because incorporate molecular grid integration into "intatomspace" routine will break the structure of the routine, I decide to write this new routine
!! dedicated to this purpose. This routine is mainly used in shubin's study.
subroutine intHirsh_molgrid
use defvar
use function
use util
implicit real*8 (a-h,o-z)
real*8 funcval(radpot*sphpot),beckeweigrid(radpot*sphpot),atmdens(radpot*sphpot,ncenter),atmintval(ncenter) !Integration value of each atom
type(content) gridatm(radpot*sphpot),gridatmorg(radpot*sphpot)

write(*,*) "Select the real space function"
call selfunc_interface(ifunc)

call setpromol
call gen1cintgrid(gridatmorg,iradcut)

!funcval: real space function at all grid of current center
!weival: Becke * single-center integraton weight at all grid of current center
!atmdens: Free-state density of every atom at all grid of current center
write(*,"(' Radial points:',i5,'    Angular points:',i5,'   Total:',i10,' per center')") radpot,sphpot,radpot*sphpot

call walltime(iwalltime1)
CALL CPU_TIME(time_begin)
atmintval=0

do iatm=1,ncenter !Cycle each atom
    write(*,"(' Processing center',i6,'(',a2,')   /',i6)") iatm,a(iatm)%name,ncenter
    gridatm%x=gridatmorg%x+a(iatm)%x !Move quadrature point to actual position in molecule
    gridatm%y=gridatmorg%y+a(iatm)%y
    gridatm%z=gridatmorg%z+a(iatm)%z
    
    !Calculate real space function value
nthreads=getNThreads()
!$OMP parallel do shared(funcval) private(i) num_threads(nthreads)
    do i=1+iradcut*sphpot,radpot*sphpot
        funcval(i)=calcfuncall(ifunc,gridatm(i)%x,gridatm(i)%y,gridatm(i)%z)
    end do
!$OMP end parallel do
    
    !Calculate Becke weight
    call gen1cbeckewei(iatm,iradcut,gridatm,beckeweigrid)
    
    !Calculate Hirshfeld weight
    do jatm=1,ncenter_org
        call dealloall
        call readwfn(custommapname(jatm),1)
nthreads=getNThreads()
!$OMP parallel do shared(atmdens) private(ipt) num_threads(nthreads)
        do ipt=1+iradcut*sphpot,radpot*sphpot
            atmdens(ipt,jatm)=fdens(gridatm(ipt)%x,gridatm(ipt)%y,gridatm(ipt)%z)
        end do
!$OMP end parallel do
    end do
    call dealloall
    call readinfile(firstfilename,1) !Retrieve to the first loaded file(whole molecule) to calc real rho again
    
    !Generate Hirshfeld weight and integrate the function
    do i=1+iradcut*sphpot,radpot*sphpot
        promol=sum(atmdens(i,:))
        do jatm=1,ncenter
            Hirshwei=atmdens(i,jatm)/promol !Hirshfeld weight of jatm at i point
            atmintval(jatm)=atmintval(jatm)+funcval(i)*Hirshwei*beckeweigrid(i)*gridatmorg(i)%value
        end do
    end do
    
end do

CALL CPU_TIME(time_end)
call walltime(iwalltime2)
write(*,"(' Calculation took up CPU time',f12.2,'s, wall clock time',i10,'s',/)") time_end-time_begin,iwalltime2-iwalltime1

do iatm=1,ncenter
    write(*,"(' Atom',i6,'(',a2,'):',f20.8)") iatm,a(iatm)%name,atmintval(iatm)
end do
write(*,"(' Total:',f20.8,/)") sum(atmintval(:)) 
end subroutine







!!------- Integrate AIM basins using mixed atomic-center and uniform grids for DFRT 2.0 project
subroutine integratebasinmix_LSB
use defvar
use util
use function
use basinintmod
use topo
implicit real*8(a-h,o-z)
real*8 trustrad(numrealatt),grad(3),hess(3,3),k1(3),k2(3),k3(3),k4(3),xarr(nx),yarr(ny),zarr(nz)
real*8,allocatable :: potx(:),poty(:),potz(:),potw(:)
type(content),allocatable :: gridatt(:) !Record x,y,z,weight of grids in trust radius
integer att2atm(numrealatt) !The attractor corresponds to which atom. If =0, means this is a NNA
integer walltime1,walltime2,radpotAIM,sphpotAIM
integer,parameter :: nfunc=7,nquant=7 !The number of real space function, the number of quantities to be calculated
real*8 arrtmp(nfunc),arrtmp2(nfunc),gradrho(3)
character*40 functionname(0:nfunc),quantityname(0:nquant)
real*8 atmcontri(ncenter,0:nfunc,0:nquant)
real*8,allocatable :: promol(:),promolgrad(:,:),atomdens(:),atomgrad(:,:),funcval(:,:),funcgrdn(:,:),funcref(:,:)
integer,allocatable :: grdidx(:),grdidy(:),grdidz(:) !Record x,y,z index of all the grids belonging to an attractor

itype=2 !Exact refinement
radpotAIM=200
nbeckeiter=8
call setpromol
expcutoff=1
functionname(0)=" rho"
functionname(1)=" rho/rho0"
functionname(2)=" |der_rho|/rho^(4/3)"
functionname(3)=" der2rho/rho^(5/3)"
functionname(4)=" (tau-t_w)/t_TF"
functionname(5)=" Xi part of SEDD"
functionname(6)=" Theta part of DORI"
functionname(7)=" Spin density"
quantityname(0)=" Function itself" !Special
quantityname(1)=" Shannon entropy"
quantityname(2)=" Fisher information"
quantityname(3)=" Onicescu information energy of order 2"
quantityname(4)=" Onicescu information energy of order 3"
quantityname(5)=" Information gain"
quantityname(6)=" Relative Renyi entropy of orders 2"
quantityname(7)=" Relative Renyi entropy of orders 3"

call walltime(walltime1)
CALL CPU_TIME(time_begin)

numcp=0
att2atm=0
atmcontri=0
!Determine trust radius and then integrate in the trust sphere
write(*,*) "Integrating in trust sphere..."

do iatt=1,numrealatt !Cycle each attractors

    !Find one-to-one correspondence between NCP and atom, then the NCP position is defined as nuclear position
    do iatm=1,ncenter
        disttest=dsqrt( (realattxyz(iatt,1)-a(iatm)%x)**2+(realattxyz(iatt,2)-a(iatm)%y)**2+(realattxyz(iatt,3)-a(iatm)%z)**2 )
        if (disttest<0.3D0) then
            att2atm(iatt)=iatm
            write(*,"(' Attractor',i6,' corresponds to atom',i6,' (',a,')')") iatt,iatm,a(iatm)%name
            numcpold=numcp
            !Refine the crude position of attractor by exact newton method
            call findcp(a(iatm)%x,a(iatm)%y,a(iatm)%z,1,0)
            if (numcp==numcpold) then
                write(*,*) "Note: Unable to locate exact CP position! Use nuclear position"
                numcp=numcp+1
                CPpos(1,numcp)=a(iatm)%x
                CPpos(2,numcp)=a(iatm)%y
                CPpos(3,numcp)=a(iatm)%z
            end if
            exit
        end if
    end do
    if (att2atm(iatt)==0) then
        write(*,"(a,i6,a)") " Warning: Unable to determine the attractor",iatt," belongs to which atom!"
        write(*,*) "Non-nuclear attractor is not supported."
        write(*,*) "Press ENTER to exit"
        read(*,*)
        return
    end if
    
    !Determine trust radius and set integration points and weight
    parm=1
    isettrustrad=0
    nintgrid=0 !Then number of integration grids within trust radius
    if (allocated(gridatt)) deallocate(gridatt) !Used to record information of grids in trust sphere of this attractor
    allocate(gridatt(radpotAIM*500))
    do ish=1,radpotAIM !Cycle each radial shell. Radius distance is from near to far
        if (isettrustrad==1) exit !The trust radius has been finally determined in last shell cycle
        !Becke, namely the second-kind Gauss-Chebyshev
        itmp=radpotAIM+1-ish !Invert ish to make radr from near to far
        radx=cos(itmp*pi/(radpotAIM+1D0))
        radr=(1+radx)/(1-radx)*parm
        radw=2*pi/(radpotAIM+1)*parm**3 *(1+radx)**2.5D0/(1-radx)**3.5D0 *4*pi
        !Set proper Lebedev grid according to shell radius
        radtmp=covr(a(iatm)%index)
        if (radr<0.2D0*radtmp) then
            sphpotAIM=26
        else if (radr<0.5D0*radtmp) then
            sphpotAIM=74
        else if (radr<0.8D0*radtmp) then
            sphpotAIM=146
        else
            sphpotAIM=194
        end if
        if (allocated(potx)) deallocate(potx,poty,potz,potw)
        allocate(potx(sphpotAIM),poty(sphpotAIM),potz(sphpotAIM),potw(sphpotAIM))
        call Lebedevgen(sphpotAIM,potx,poty,potz,potw)
        !Combine radial point and weights with angular part, and make them centered at current attractor
        gridatt( nintgrid+1:nintgrid+sphpotAIM )%x=radr*potx+CPpos(1,numcp)
        gridatt( nintgrid+1:nintgrid+sphpotAIM )%y=radr*poty+CPpos(2,numcp)
        gridatt( nintgrid+1:nintgrid+sphpotAIM )%z=radr*potz+CPpos(3,numcp)
        gridatt( nintgrid+1:nintgrid+sphpotAIM )%value=radw*potw
        !Find trust radius for present attractor
        angmax=0
        radrinit=0.15D0
        if (a(iatm)%index>2) radrinit=0.5D0
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
                    isettrustrad=1
                    exit
                end if
            end do
            if (isettrustrad==0) trustrad(iatt)=radr !Passed this shell and temporarily set the radius as trust radius. Continue to enlarge the trust radius, until reached angmax>45 degree
        end if
        nintgrid=nintgrid+sphpotAIM
    end do
    if (isettrustrad==0) trustrad(iatt)=20
    write(*,"(' The trust radius of attractor',i6,' is',f10.3,' Bohr',/)") iatt,trustrad(iatt)
    
    allocate(promol(nintgrid),promolgrad(nintgrid,3),atomdens(nintgrid),atomgrad(nintgrid,3),funcval(nintgrid,0:nfunc),funcgrdn(nintgrid,0:nfunc),funcref(nintgrid,0:nfunc))
    
    !! Calculate function value of reference state, including promolecular density and density of present atom, function index:
    !index=0: rho
    !index=1: rho/rho0, where rho0 is promolecular density
    !index=2: |der_rho|/rho^(4/3)
    !index=3: der2rho/rho^(5/3)
    !index=4: (tau-t_w)/t_TF
    !index=5: Part of SEDD
    !index=6: Part of DORI
    !index=7: Spin density
    funcref=0D0
    promol=0D0
    promolgrad=0D0
    do jatm=1,ncenter_org
        call dealloall
        call readwfn(custommapname(jatm),1)
nthreads=getNThreads()
!$OMP parallel do shared(atomdens,atomgrad,funcref) private(ipt,arrtmp,rho,gradrho) num_threads(nthreads)
        do ipt=1,nintgrid
            call valaryyLSB(gridatt(ipt)%x,gridatt(ipt)%y,gridatt(ipt)%z,arrtmp,rho,gradrho)
            atomdens(ipt)=rho
            atomgrad(ipt,:)=gradrho
            funcref(ipt,0)=funcref(ipt,0)+rho
            funcref(ipt,1)=funcref(ipt,1)+1D0 !i.e. rho/rho0, the value is unity if present system is a single atom
            funcref(ipt,2:nfunc)=funcref(ipt,2:nfunc)+arrtmp(2:nfunc)
        end do
!$OMP end parallel do
        promol=promol+atomdens(:)
        promolgrad=promolgrad+atomgrad(:,:)
    end do
    call dealloall
    call readinfile(firstfilename,1) !Retrieve the firstly loaded file(whole molecule)
    
    !! Calculate function value and gradient norm for present molecule, then store them in funcval/funcgrdn
nthreads=getNThreads()
!$OMP parallel do shared(funcval,funcgrdn) private(ipt,arrtmp,arrtmp2,rho,gradrho,idir) num_threads(nthreads)
    do ipt=1,nintgrid
        call valgradarrLSB(gridatt(ipt)%x,gridatt(ipt)%y,gridatt(ipt)%z,arrtmp,arrtmp2,rho,gradrho)
        funcval(ipt,0)=rho
        funcgrdn(ipt,0)=dsqrt(sum(gradrho**2))
        !Note: see word document for explicit expression of |der(rho/rho0)|
        funcval(ipt,1)=rho/promol(ipt)
        funcgrdn(ipt,1)=0D0
        do idir=1,3
            funcgrdn(ipt,1)=funcgrdn(ipt,1)+( gradrho(idir)/promol(ipt)-rho/promol(ipt)**2*promolgrad(ipt,idir) )**2
        end do
        funcgrdn(ipt,1)=dsqrt(funcgrdn(ipt,1))
        funcval(ipt,2:nfunc)=arrtmp(2:nfunc)
        funcgrdn(ipt,2:nfunc)=arrtmp2(2:nfunc)
    end do
!$OMP end parallel do
    
    !! Integrate the region inside trust radius
    do ipt=1,nintgrid
        rx=gridatt(ipt)%x-CPpos(1,numcp) !The relative distance between current point to corresponding attractor
        ry=gridatt(ipt)%y-CPpos(2,numcp)
        rz=gridatt(ipt)%z-CPpos(3,numcp)
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
        if (switchwei<1D-7) cycle !For saving computational time
        
        weitot=switchwei*gridatt(ipt)%value
        do ifunc=0,nfunc
            !Value of function itself
            atmcontri(iatm,ifunc,0)=atmcontri(iatm,ifunc,0)+weitot*funcval(ipt,ifunc)
            if (ifunc/=3) then !Because the variable in log must be positive everywhere
                !Shannon entropy
                atmcontri(iatm,ifunc,1)=atmcontri(iatm,ifunc,1)+weitot*( -funcval(ipt,ifunc)*log(funcval(ipt,ifunc)) )
                !Information gain
                atmcontri(iatm,ifunc,5)=atmcontri(iatm,ifunc,5)+weitot*( funcval(ipt,ifunc)*log(funcval(ipt,ifunc)/funcref(ipt,ifunc)) )
            end if
            !Fisher information
            atmcontri(iatm,ifunc,2)=atmcontri(iatm,ifunc,2)+weitot*( funcgrdn(ipt,ifunc)**2 /funcval(ipt,ifunc) )
            !Onicescu information energy of order 2
            atmcontri(iatm,ifunc,3)=atmcontri(iatm,ifunc,3)+weitot*( funcval(ipt,ifunc)**2 )
            !Onicescu information energy of order 3
            atmcontri(iatm,ifunc,4)=atmcontri(iatm,ifunc,4)+weitot*( funcval(ipt,ifunc)**3 )
            !Relative Renyi entropy of orders 2
            atmcontri(iatm,ifunc,6)=atmcontri(iatm,ifunc,6)+weitot*( funcval(ipt,ifunc)**2/funcref(ipt,ifunc) )
            !Relative Renyi entropy of orders 3
            atmcontri(iatm,ifunc,7)=atmcontri(iatm,ifunc,7)+weitot*( funcval(ipt,ifunc)**3/funcref(ipt,ifunc)**2 )
        end do
    end do
    
!     write(*,"(' ========== Contribution of atom',i4,a,':')") iatm,a(iatm)%name
!     do ifunc=0,nfunc
!         write(*,"(/,' Function:',a)") trim(functionname(ifunc))
!         do iquant=0,nquant !0 corresponds to function itself
!             if (ifunc==3.and.(iquant==1.or.iquant==5)) cycle
!             write(*,"(a,':',1E16.8)") quantityname(iquant),atmcontri(iatm,ifunc,iquant)
!         end do
!     end do
    
    deallocate(promol,promolgrad,atomdens,atomgrad,funcval,funcgrdn,funcref)
end do !End cycle attractors

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

!---------------------------------------------!
!--------- Integrating uniform grids ---------!
!---------------------------------------------!

!--------- Exact refinement boundary grids
write(*,*) "Refining boundary grids..."
if (itype==2.or.itype==3) then
    if (itype==3) then !Calculate grid data of gradient of electron density used to linear interpolation to obtain the value at any point
        write(*,*)
        call gengradmat
    end if
    nrk4lim=100
    nrk4gradswitch=40
    hsizeinit=0.25D0
    ifinish=0
!$OMP PARALLEL do private(ix,iy,iz,iatt,rnowx,rnowy,rnowz,rx,ry,rz,tmpval,tmpval2,tmpval3,&
!$OMP rnowxtmp,rnowytmp,rnowztmp,orgxref,orgyref,orgzref,dxref,dyref,dzref,ixref,iyref,izref,nrefine,ndiv,&
!$OMP k1,k2,k3,k4,dens,denshold,grad,hess,iattref,xtmp,ytmp,ztmp,irk4,hsize,ixtest,iytest,iztest,tmpdist) shared(ifinish) NUM_THREADS(nthreads) schedule(DYNAMIC)
    do iz=2,nz-1
        rnowz=zarr(iz)
        do iy=2,ny-1
            rnowy=yarr(iy)
            do ix=2,nx-1
                rnowx=xarr(ix)
                if (interbasgrid(ix,iy,iz).eqv. .false.) cycle
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
                            end if
                            gridbas(ix,iy,iz)=iattref !Update attribution of boundary grids
                        end do !End refine grid
                    end do
                end do
                
            end do !End cycle ix grid
        end do
        ifinish=ifinish+1
    end do
!$OMP end PARALLEL do
    call detectinterbasgrd(6)
    write(*,*) "Basin boundary has been updated"
end if

!!---------- Calculate part of quantities contributed from uniform grid, and gain final result
write(*,*)
write(*,*) "Calculating information at uniform grid..."
do iatt=1,numrealatt !Cycle each attractors
    write(*,"(' Processing basin',i5,' /',i5)") iatt,numrealatt
    iatm=att2atm(iatt)
    nintgrid=count(gridbas(2:nx-1,2:ny-1,2:nz-1)==iatt)
    allocate(promol(nintgrid),promolgrad(nintgrid,3),atomdens(nintgrid),atomgrad(nintgrid,3),funcval(nintgrid,0:nfunc),funcgrdn(nintgrid,0:nfunc),funcref(nintgrid,0:nfunc))
    allocate(grdidx(nintgrid),grdidy(nintgrid),grdidz(nintgrid))
    itmp=0
    do iz=2,nz-1
        do iy=2,ny-1
            do ix=2,nx-1
                if (gridbas(ix,iy,iz)==iatt) then !Find grids attributing to this attractor
                    itmp=itmp+1
                    grdidx(itmp)=ix
                    grdidy(itmp)=iy
                    grdidz(itmp)=iz
                end if
            end do
        end do
    end do
    
    !Calculate information for promolecular state
    funcref=0D0
    promol=0D0
    promolgrad=0D0
    do jatm=1,ncenter_org
        call dealloall
        call readwfn(custommapname(jatm),1)
nthreads=getNThreads()
!$OMP parallel do shared(atomdens,atomgrad,funcref) private(ipt,ptx,pty,ptz,arrtmp,rho,gradrho) num_threads(nthreads)
        do ipt=1,nintgrid
            ptx=xarr(grdidx(ipt))
            pty=yarr(grdidy(ipt))
            ptz=zarr(grdidz(ipt))
            call valaryyLSB(ptx,pty,ptz,arrtmp,rho,gradrho)
            atomdens(ipt)=rho
            atomgrad(ipt,:)=gradrho
            funcref(ipt,0)=funcref(ipt,0)+rho
            funcref(ipt,1)=funcref(ipt,1)+1D0 !i.e. rho/rho0, the value is unity if present system is a single atom
            funcref(ipt,2:nfunc)=funcref(ipt,2:nfunc)+arrtmp(2:nfunc)
        end do
!$OMP end parallel do
        promol=promol+atomdens(:)
        promolgrad=promolgrad+atomgrad(:,:)
    end do
    call dealloall
    call readinfile(firstfilename,1) !Retrieve the firstly loaded file(whole molecule)
    
    !Calculate information for present actual state
nthreads=getNThreads()
!$OMP PARALLEL do shared(funcval,funcgrdn) private(ipt,ptx,pty,ptz,arrtmp,arrtmp2,rho,gradrho,idir) NUM_THREADS(nthreads) schedule(DYNAMIC)
    do ipt=1,nintgrid
        ptx=xarr(grdidx(ipt))
        pty=yarr(grdidy(ipt))
        ptz=zarr(grdidz(ipt))
        call valgradarrLSB(ptx,pty,ptz,arrtmp,arrtmp2,rho,gradrho)
        funcval(ipt,0)=rho
        funcgrdn(ipt,0)=dsqrt(sum(gradrho**2))
        !Note: see word document for explicit expression of |der(rho/rho0)|
        funcval(ipt,1)=rho/promol(ipt)
        funcgrdn(ipt,1)=0D0
        do idir=1,3
            funcgrdn(ipt,1)=funcgrdn(ipt,1)+( gradrho(idir)/promol(ipt)-rho/promol(ipt)**2*promolgrad(ipt,idir) )**2
        end do
        funcgrdn(ipt,1)=dsqrt(funcgrdn(ipt,1))
        funcval(ipt,2:nfunc)=arrtmp(2:nfunc)
        funcgrdn(ipt,2:nfunc)=arrtmp2(2:nfunc)
    end do
!$OMP end PARALLEL do
    
    !Accumulate uniform grid contribution to quantities
    do ipt=1,nintgrid
        ptx=xarr(grdidx(ipt))
        pty=yarr(grdidy(ipt))
        ptz=zarr(grdidz(ipt))
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
        
        weitot=switchwei*dvol
        do ifunc=0,nfunc
            !Value of function itself
            atmcontri(iatm,ifunc,0)=atmcontri(iatm,ifunc,0)+weitot*funcval(ipt,ifunc)
            if (ifunc/=3) then !Because the variable in log must be positive everywhere
                !Shannon entropy
                atmcontri(iatm,ifunc,1)=atmcontri(iatm,ifunc,1)+weitot*( -funcval(ipt,ifunc)*log(funcval(ipt,ifunc)) )
                !Information gain
                atmcontri(iatm,ifunc,5)=atmcontri(iatm,ifunc,5)+weitot*( funcval(ipt,ifunc)*log(funcval(ipt,ifunc)/funcref(ipt,ifunc)) )
            end if
            !Fisher information
            atmcontri(iatm,ifunc,2)=atmcontri(iatm,ifunc,2)+weitot*( funcgrdn(ipt,ifunc)**2 /funcval(ipt,ifunc) )
            !Onicescu information energy of order 2
            atmcontri(iatm,ifunc,3)=atmcontri(iatm,ifunc,3)+weitot*( funcval(ipt,ifunc)**2 )
            !Onicescu information energy of order 3
            atmcontri(iatm,ifunc,4)=atmcontri(iatm,ifunc,4)+weitot*( funcval(ipt,ifunc)**3 )
            !Relative Renyi entropy of orders 2
            atmcontri(iatm,ifunc,6)=atmcontri(iatm,ifunc,6)+weitot*( funcval(ipt,ifunc)**2/funcref(ipt,ifunc) )
            !Relative Renyi entropy of orders 3
            atmcontri(iatm,ifunc,7)=atmcontri(iatm,ifunc,7)+weitot*( funcval(ipt,ifunc)**3/funcref(ipt,ifunc)**2 )
        end do
    end do
    
    deallocate(promol,promolgrad,atomdens,atomgrad,funcval,funcgrdn,funcref,grdidx,grdidy,grdidz)
end do !End cycling attractor

do iatm=1,ncenter
    write(*,"(/,' ========== Contribution of atom',i4,a,':')") iatm,a(iatm)%name
    do ifunc=0,nfunc
        write(*,"(/,' Function:',a)") trim(functionname(ifunc))
        do iquant=0,nquant !0 corresponds to function itself
            if (ifunc==3.and.(iquant==1.or.iquant==5)) cycle
            write(*,"(a,':',1E16.8)") quantityname(iquant),atmcontri(iatm,ifunc,iquant)
        end do
    end do
end do
write(*,*)
write(*,*) "========================================="
write(*,*) "             Overall result"
write(*,*) "========================================="
do ifunc=0,nfunc
    write(*,"(/,' Function:',a)") trim(functionname(ifunc))
    do iquant=0,nquant
        if (ifunc==3.and.(iquant==1.or.iquant==5)) cycle
        write(*,"(a,':',1E16.8)") quantityname(iquant),sum(atmcontri(:,ifunc,iquant))
    end do
end do

CALL CPU_TIME(time_end)
call walltime(walltime2)
write(*,"(' Integrating basins took up CPU time',f12.2,'s, wall clock time',i10,'s')") time_end-time_begin,walltime2-walltime1

end subroutine
