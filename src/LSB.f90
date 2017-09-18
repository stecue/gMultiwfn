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
integer,parameter :: nfunc=6,nquant=7 !The number of real space function, the number of quantities to be calculated
real*8 smat(ncenter,ncenter),Pvec(ncenter),atmpartwei(radpot*sphpot)
real*8 atmcontri(ncenter,0:nfunc,0:nquant),checkacc(ncenter)
real*8 promol(radpot*sphpot),atomdens(radpot*sphpot),selfdens(radpot*sphpot) !For Hirshfeld partition
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
isupplyEDF=0 !Avoid EDF complicating things
covr_becke=covr_tianlu

functionname(0)=" rho"
functionname(1)=" rho/rho0"
functionname(2)=" |der_rho|/rho^(4/3)"
functionname(3)=" der2rho/rho^(5/3)"
functionname(4)=" (tau-t_w)/t_TF"
functionname(5)=" Xi part of SEDD"
functionname(6)=" Theta part of DORI"
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
    
    
    !! Calculate function value of reference state, including promolecular density and density of present atom, function index:
    !index=0: rho
    !index=1: rho/rho0, where rho0 is promolecular density
    !index=2: |der_rho|/rho^(4/3)
    !index=3: der2rho/rho^(5/3)
    !index=4: (tau-t_w)/t_TF
    !index=5: SEDD
    !index=6: DORI
    funcref=0D0
    promol=0D0
    do jatm=1,ncenter_org
        call dealloall
        call readwfn(custommapname(jatm),1)
nthreads=getNThreads()
!$OMP parallel do shared(atomdens,funcref) private(ipt,rnowx,rnowy,rnowz,arrtmp,rho,gradrho) num_threads(nthreads)
        do ipt=1+iradcut*sphpot,radpot*sphpot
            rnowx=gridatm(ipt)%x
            rnowy=gridatm(ipt)%y
            rnowz=gridatm(ipt)%z
            call valaryyLSB(rnowx,rnowy,rnowz,arrtmp,rho,gradrho)
            atomdens(ipt)=rho
            funcref(ipt,0)=funcref(ipt,0)+rho
            funcref(ipt,1)=funcref(ipt,1)+1D0 !i.e. rho/rho0, the value is unity if present system is a single atom
            funcref(ipt,2:nfunc)=funcref(ipt,2:nfunc)+arrtmp(2:nfunc)
        end do
!$OMP end parallel do
        promol=promol+atomdens
        if (jatm==iatm) selfdens=atomdens
    end do
    call dealloall
    call readinfile(firstfilename,1) !Retrieve the firstly loaded file(whole molecule)
    
    
    !! Calculate function value and gradient norm for present molecule, then store them in funcval/funcgrdn. Function index:
    !index=0: rho
    !index=1: rho/rho0, where rho0 is promolecular density
    !index=2: |der_rho|/rho^(4/3)
    !index=3: der2rho/rho^(5/3)
    !index=4: (tau-t_w)/t_TF
    !index=5: SEDD
    !index=6: DORI
nthreads=getNThreads()
!$OMP parallel do shared(funcval,funcgrdn) private(ipt,rnowx,rnowy,rnowz,arrtmp,arrtmp2,rho,gradrho,idir) num_threads(nthreads)
    do ipt=1+iradcut*sphpot,radpot*sphpot
        rnowx=gridatm(ipt)%x
        rnowy=gridatm(ipt)%y
        rnowz=gridatm(ipt)%z
        call valgradarrLSB(rnowx,rnowy,rnowz,arrtmp,arrtmp2,rho,gradrho)
        funcval(ipt,0)=rho
        funcgrdn(ipt,0)=dsqrt(sum(gradrho**2))
        !Note: see word document for explicit expression of |der(rho/rho0)|
        funcval(ipt,1)=rho/promol(ipt)
        funcgrdn(ipt,1)=0D0
        do idir=1,3
            funcgrdn(ipt,1)=funcgrdn(ipt,1)+( gradrho(idir)*(1/promol(ipt)-rho/promol(ipt)**2) )**2
        end do
        funcgrdn(ipt,1)=dsqrt(funcgrdn(ipt,1))
        funcval(ipt,2:nfunc)=arrtmp(2:nfunc)
        funcgrdn(ipt,2:nfunc)=arrtmp2(2:nfunc)
    end do
!$OMP end parallel do
    
    
    !! Calculate atomic space weight for present atom at all of its points
    if (ipartition==1) then !Becke
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
            atmpartwei(i)=Pvec(iatm) !Normalized Pvec, Pvec contain partition weight of each atom in current point, namely i
        end do
!$OMP end parallel do
    else if (ipartition==2) then !Hirshfeld based on atomic .wfn files
        do i=1+iradcut*sphpot,radpot*sphpot !Get Hirshfeld weight of present atom
            if (promol(i)/=0D0) then
                atmpartwei(i)=selfdens(i)/promol(i)
            else
                atmpartwei(i)=0D0
            end if
        end do
    end if
    
    
    !! Calculate all quantities contributed by this atom
    do ipt=1+iradcut*sphpot,radpot*sphpot
        do ifunc=0,nfunc
            !Value of function itself
            atmcontri(iatm,ifunc,0)=atmcontri(iatm,ifunc,0)+atmpartwei(ipt)*gridatm(ipt)%value*funcval(ipt,ifunc)
            if (ifunc/=3) then !Because the variable in log must be positive everywhere
                !Shannon entropy
                atmcontri(iatm,ifunc,1)=atmcontri(iatm,ifunc,1)+atmpartwei(ipt)*gridatm(ipt)%value*( -funcval(ipt,ifunc)*log(funcval(ipt,ifunc)) )
                !Information gain
                atmcontri(iatm,ifunc,5)=atmcontri(iatm,ifunc,5)+atmpartwei(ipt)*gridatm(ipt)%value*( funcval(ipt,ifunc)*log(funcval(ipt,ifunc)/funcref(ipt,ifunc)) )
            end if
            !Fisher information
            atmcontri(iatm,ifunc,2)=atmcontri(iatm,ifunc,2)+atmpartwei(ipt)*gridatm(ipt)%value*( funcgrdn(ipt,ifunc)**2 /funcval(ipt,ifunc) )
            !Onicescu information energy of order 2
            atmcontri(iatm,ifunc,3)=atmcontri(iatm,ifunc,3)+atmpartwei(ipt)*gridatm(ipt)%value*( funcval(ipt,ifunc)**2 )
            !Onicescu information energy of order 3
            atmcontri(iatm,ifunc,4)=atmcontri(iatm,ifunc,4)+atmpartwei(ipt)*gridatm(ipt)%value*( funcval(ipt,ifunc)**3 )
            !Relative Renyi entropy of orders 2
            atmcontri(iatm,ifunc,6)=atmcontri(iatm,ifunc,6)+atmpartwei(ipt)*gridatm(ipt)%value*( funcval(ipt,ifunc)**2/funcref(ipt,ifunc) )
            !Relative Renyi entropy of orders 3
            atmcontri(iatm,ifunc,7)=atmcontri(iatm,ifunc,7)+atmpartwei(ipt)*gridatm(ipt)%value*( funcval(ipt,ifunc)**3/funcref(ipt,ifunc)**2 )
        end do
!         checkacc(iatm)=checkacc(iatm)+atmpartwei(ipt)*gridatm(ipt)%value*funcval(ipt,0)
    end do
    
    write(*,"(' Contribution of atom',i4,a,':')") iatm,a(iatm)%name
    do ifunc=0,nfunc
        write(*,"(/,' Function:',a)") trim(functionname(ifunc))
        do iquant=0,nquant !0 corresponds to function itself
            if (ifunc==3.and.(iquant==1.or.iquant==5)) cycle
            write(*,"(a,':',1E16.8)") quantityname(iquant),atmcontri(iatm,ifunc,iquant)
        end do
    end do
!     write(*,"(f12.6)") checkacc(iatm)
    
end do !End cycling atoms

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
! write(*,*) sum(checkacc(:))

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
end subroutine


!Get value and norm of numerical gradient for all the functions used in LSB's DFRT 2.0 project
!The definition of slot is identical to valaryyLSB
subroutine valgradarrLSB(x,y,z,valarr,grdnarr,rho,gradrho)
use defvar
implicit real*8 (a-h,o-z)
integer,parameter :: nfunc=6
real*8 x,y,z,rho,gradrho(3),valarr(nfunc),grdnarr(nfunc),tmparr(3)
real*8 xadd(nfunc),xmin(nfunc),yadd(nfunc),ymin(nfunc),zadd(nfunc),zmin(nfunc),grdxarr(nfunc),grdyarr(nfunc),grdzarr(nfunc)
diff=1D-3
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
