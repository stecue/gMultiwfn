!!--------- Use suitable routine to read in the file, infomode=0/1 show/don't show info
subroutine readinfile(thisfilename,infomode)
use defvar
implicit real*8 (a-h,o-z)
character(len=*) thisfilename
integer infomode,inamelen
inamelen=len_trim(thisfilename)
if (infomode==0) write(*,*) "Please wait..."
if (thisfilename(inamelen-2:inamelen)=="fch".or.thisfilename(inamelen-3:inamelen)=="fchk"&
.or.thisfilename(inamelen-3:inamelen)=="FCH".or.thisfilename(inamelen-3:inamelen)=="FChk"&
.or.thisfilename(inamelen-3:inamelen)=="FCHK") then
    call readfch(thisfilename,infomode)
else if (thisfilename(inamelen-2:inamelen)=="wfn".or.thisfilename(inamelen-2:inamelen)=="WFN") then
    call readwfn(thisfilename,infomode)
else if (thisfilename(inamelen-2:inamelen)=="wfx".or.thisfilename(inamelen-2:inamelen)=="WFX") then
    call readwfx(thisfilename,infomode)
else if (thisfilename(inamelen-2:inamelen)=="chg") then
    call readchg(thisfilename,infomode)
else if (thisfilename(inamelen-2:inamelen)=="pdb".or.thisfilename(inamelen-2:inamelen)=="PDB") then
    call readpdb(thisfilename,infomode)
else if (thisfilename(inamelen-2:inamelen)=="xyz".or.thisfilename(inamelen-2:inamelen)=="XYZ") then
    call readxyz(thisfilename,infomode,1)
else if (thisfilename(inamelen-1:inamelen)=="31") then
    call read31(thisfilename,infomode)
else if (thisfilename(inamelen-2:inamelen)=="grd") then
    call readgrd(thisfilename,infomode,0)
else if (thisfilename(inamelen-2:inamelen)=="cub".or.thisfilename(inamelen-3:inamelen)=="cube") then
    call readcube(thisfilename,infomode,0)
else if (len_trim(thisfilename)>=6) then
    if  (index(thisfilename,".molden")/=0) then
        call readmolden(thisfilename,infomode)
    end if
else
    ifiletype=0
end if
!Load EDF information from external atomic .wfx files for .fch/.wfn/.wfx/.molden
if (ireadatmEDF==1.and.(.not.allocated(b_EDF)).and.(ifiletype==1.or.ifiletype==2.or.ifiletype==3.or.ifiletype==9)) then
    call readatmEDF
end if
end subroutine



!!----------- Generate spherical harmonic -> Cartesian basis function conversion table for d,f,g,h.
!The table comes from IJQC,54,83, which is used by Gaussian and also identical to the one used by Molden
subroutine gensphcartab(matd,matf,matg,math)
real*8 matd(6,5),matf(10,7),matg(15,9),math(21,11)
matd=0D0
matf=0D0
matg=0D0
math=0D0
! From 5D: D 0,D+1,D-1,D+2,D-2
! To 6D:  1  2  3  4  5  6
!        XX,YY,ZZ,XY,XZ,YZ
!
! D0=-0.5*XX-0.5*YY+ZZ
matd(1:3,1)=(/ -0.5D0,-0.5D0,1D0 /)
! D+1=XZ
matd(5,2)=1D0
! D-1=YZ
matd(6,3)=1D0
! D+2=SQRT(3)/2*(XX-YY)
matd(1:2,4)=(/ sqrt(3D0)/2D0,-sqrt(3D0)/2D0 /)
! D-2=XY
matd(4,5)=1D0

! From 7F: F 0,F+1,F-1,F+2,F-2,F+3,F-3
! To 10F:  1   2   3   4   5   6   7   8   9  10      
!         XXX,YYY,ZZZ,XYY,XXY,XXZ,XZZ,YZZ,YYZ,XYZ (Gaussian sequence, not identical to Multiwfn)
!
! F 0=-3/(2*¡Ì5)*(XXZ+YYZ)+ZZZ
matf(3,1)=1D0
matf(6,1)=-1.5D0/sqrt(5D0)
matf(9,1)=-1.5D0/sqrt(5D0)
! F+1=-¡Ì(3/8)*XXX-¡Ì(3/40)*XYY+¡Ì(6/5)*XZZ
matf(1,2)=-sqrt(3D0/8D0)
matf(4,2)=-sqrt(3D0/40D0)
matf(7,2)=sqrt(6D0/5D0)
! F-1=-¡Ì(3/40)*XXY-¡Ì(3/8)*YYY+¡Ì(6/5)*YZZ
matf(2,3)=-sqrt(3D0/8D0)
matf(5,3)=-sqrt(3D0/40D0)
matf(8,3)=sqrt(6D0/5D0)
! F+2=¡Ì3/2*(XXZ-YYZ)
matf(6,4)=sqrt(3D0)/2D0
matf(9,4)=-sqrt(3D0)/2D0
! F-2=XYZ
matf(10,5)=1D0
! F+3=¡Ì(5/8)*XXX-3/¡Ì8*XYY
matf(1,6)=sqrt(5D0/8D0)
matf(4,6)=-3D0/sqrt(8D0)
! F-3=3/¡Ì8*XXY-¡Ì(5/8)*YYY
matf(2,7)=-sqrt(5D0/8D0)
matf(5,7)=3D0/sqrt(8D0)

! From 9G: G 0,G+1,G-1,G+2,G-2,G+3,G-3,G+4,G-4
! To 15G:   1    2    3    4    5    6    7    8
!         ZZZZ,YZZZ,YYZZ,YYYZ,YYYY,XZZZ,XYZZ,XYYZ
!           9   10   11   12   13   14   15
!         XYYY,XXZZ,XXYZ,XXYY,XXXZ,XXXY,XXXX
!
!G 0=ZZZZ+3/8*(XXXX+YYYY)-3*¡Ì(3/35)*(XXZZ+YYZZ-1/4*XXYY)
matg(1,1)=1D0
matg(3,1)=-3D0*sqrt(3D0/35D0)
matg(5,1)=3D0/8D0
matg(10,1)=-3D0*sqrt(3D0/35D0)
matg(12,1)=3D0/4D0*sqrt(3D0/35D0)
matg(15,1)=3D0/8D0
!G+1=2*¡Ì(5/14)*XZZZ-3/2*¡Ì(5/14)*XXXZ-3/2/¡Ì14*XYYZ
matg(6,2)=2D0*sqrt(5D0/14D0)
matg(8,2)=-1.5D0/sqrt(14D0)
matg(13,2)=-1.5D0*sqrt(5D0/14D0)
!G-1=2*¡Ì(5/14)*YZZZ-3/2*¡Ì(5/14)*YYYZ-3/2/¡Ì14*XXYZ
matg(2,3)=2D0*sqrt(5D0/14D0)
matg(4,3)=-1.5D0*sqrt(5D0/14D0)
matg(11,3)=-1.5D0/sqrt(14D0)
!G+2=3*¡Ì(3/28)*(XXZZ-YYZZ)-¡Ì5/4*(XXXX-YYYY)
matg(3,4)=-3D0*sqrt(3D0/28D0)
matg(5,4)=sqrt(5D0)/4D0
matg(10,4)=3D0*sqrt(3D0/28D0)
matg(15,4)=-sqrt(5D0)/4D0
!G-2=3/¡Ì7*XYZZ-¡Ì(5/28)*(XXXY+XYYY)
matg(7,5)=3D0/sqrt(7D0)
matg(9,5)=-sqrt(5D0/28D0)
matg(14,5)=-sqrt(5D0/28D0)
!G+3=¡Ì(5/8)*XXXZ-3/¡Ì8*XYYZ
matg(8,6)=-3D0/sqrt(8D0)
matg(13,6)=sqrt(5D0/8D0)
!G-3=-¡Ì(5/8)*YYYZ+3/¡Ì8*XXYZ
matg(4,7)=-sqrt(5D0/8D0)
matg(11,7)=3D0/sqrt(8D0)
!G+4=¡Ì35/8*(XXXX+YYYY)-3/4*¡Ì3*XXYY
matg(5,8)=sqrt(35D0)/8D0
matg(12,8)=-3D0/4D0*sqrt(3D0)
matg(15,8)=sqrt(35D0)/8D0
!G-4=¡Ì5/2*(XXXY-XYYY)
matg(9,9)=-sqrt(5D0)/2D0
matg(14,9)=sqrt(5D0)/2D0

! From 11H: H 0,H+1,H-1,H+2,H-2,H+3,H-3,H+4,H-4,H+5,H-5
! To 21H:   1     2     3     4     5     6     7     8     9    10
!         ZZZZZ YZZZZ YYZZZ YYYZZ YYYYZ YYYYY XZZZZ XYZZZ XYYZZ XYYYZ 
!          11    12    13    14    15    16    17    18    19    20    21
!         XYYYY XXZZZ XXYZZ XXYYZ XXYYY XXXZZ XXXYZ XXXYY XXXXZ XXXXY XXXXX
!
!H 0=ZZZZZ-5/¡Ì21*(XXZZZ+YYZZZ)+5/8*(XXXXZ+YYYYZ)+¡Ì(15/7)/4*XXYYZ
math(1,1)=1D0
math(12,1)=-5D0/sqrt(21D0)
math(3,1)=-5D0/sqrt(21D0)
math(19,1)=5D0/8D0
math(5,1)=5D0/8D0
math(14,1)=sqrt(15D0/7D0)/4D0
!H+1=¡Ì(5/3)*XZZZZ-3*¡Ì(5/28)*XXXZZ-3/¡Ì28*XYYZZ+¡Ì15/8*XXXXX+¡Ì(5/3)/8*XYYYY+¡Ì(5/7)/4*XXXYY
math(7,2)=sqrt(5D0/3D0)
math(16,2)=-3D0*sqrt(5D0/28D0)
math(9,2)=-3D0/sqrt(28D0)
math(21,2)=sqrt(15D0)/8D0
math(11,2)=sqrt(5D0/3D0)/8D0
math(18,2)=sqrt(5D0/7D0)/4D0
!H-1=¡Ì(5/3)*YZZZZ-3*¡Ì(5/28)*YYYZZ-3/¡Ì28*XXYZZ+¡Ì15/8*YYYYY+¡Ì(5/3)/8*XXXXY+¡Ì(5/7)/4*XXYYY
math(2,3)=sqrt(5D0/3D0)
math(4,3)=-3D0*sqrt(5D0/28D0)
math(13,3)=-3D0/sqrt(28D0)
math(6,3)=sqrt(15D0)/8D0
math(20,3)=sqrt(5D0/3D0)/8D0
math(15,3)=sqrt(5D0/7D0)/4D0
!H+2=¡Ì5/2*(XXZZZ-YYZZZ)-¡Ì(35/3)/4*(XXXXZ-YYYYZ)
math(12,4)=sqrt(5D0)/2D0
math(3,4)=-sqrt(5D0)/2D0
math(19,4)=-sqrt(35D0/3D0)/4D0
math(5,4)=sqrt(35D0/3D0)/4D0
!H-2=¡Ì(5/3)*XYZZZ-¡Ì(5/12)*(XXXYZ+XYYYZ)
math(8,5)=sqrt(5D0/3D0)
math(17,5)=-sqrt(5D0/12D0)
math(10,5)=-sqrt(5D0/12D0)
!H+3=¡Ì(5/6)*XXXZZ-¡Ì(3/2)*XYYZZ-¡Ì(35/2)/8*(XXXXX-XYYYY)+¡Ì(5/6)/4*XXXYY
math(16,6)=sqrt(5D0/6D0)
math(9,6)=-sqrt(1.5D0)
math(21,6)=-sqrt(17.5D0)/8D0
math(11,6)=sqrt(17.5D0)/8D0
math(18,6)=sqrt(5D0/6D0)/4D0
!H-3=-¡Ì(5/6)*YYYZZ+¡Ì(3/2)*XXYZZ-¡Ì(35/2)/8*(XXXXY-YYYYY)-¡Ì(5/6)/4*XXYYY
math(4,7)=-sqrt(5D0/6D0)
math(13,7)=sqrt(1.5D0)
math(20,7)=-sqrt(17.5D0)/8D0
math(6,7)=sqrt(17.5D0)/8D0
math(15,7)=-sqrt(5D0/6D0)/4D0
!H+4=¡Ì35/8*(XXXXZ+YYYYZ)-3/4*¡Ì3*XXYYZ
math(19,8)=sqrt(35D0)/8D0
math(5,8)=sqrt(35D0)/8D0
math(14,8)=-0.75D0*sqrt(3D0)
!H-4=¡Ì5/2*(XXXYZ-XYYYZ)
math(17,9)=sqrt(5D0)/2D0
math(10,9)=-sqrt(5D0)/2D0
!H+5=3/8*¡Ì(7/2)*XXXXX+5/8*¡Ì(7/2)*XYYYY-5/4*¡Ì(3/2)*XXXYY
math(21,10)=3D0/8D0*sqrt(3.5D0)
math(11,10)=5D0/8D0*sqrt(3.5D0)
math(18,10)=-1.25D0*sqrt(1.5D0)
!H-5=3/8*¡Ì(7/2)*YYYYY+5/8*¡Ì(7/2)*XXXXY-5/4*¡Ì(3/2)*XXYYY
math(6,11)=3D0/8D0*sqrt(3.5D0)
math(20,11)=5D0/8D0*sqrt(3.5D0)
math(15,11)=-1.25D0*sqrt(1.5D0)
end subroutine



!!-----------------------------------------------------------------
!!------------------------- Read gaussian formatted check file
!I am trying to make this routine compatible with Q-Chem .fch file
!Some temporary arrays, including shelltype, shell2atom, shellcon, primexp, concoeff will be store to corresponding global arrays after some revisions
subroutine readfch(name,infomode) !infomode=0 means output info, =1 silent
use defvar
use util
implicit real*8 (a-h,o-z)
character(len=*) name
character c80*80,fchtitle*79 !c80 for temporary store text
real*8 temp
integer :: i,j,k,l,ifound
integer,allocatable :: shelltype(:),shell2atom(:),shellcon(:) !Degree of shell contraction
real*8,allocatable :: primexp(:),concoeff(:),SPconcoeff(:),amocoeff(:,:),bmocoeff(:,:),tmpmat(:,:),tmparr(:)
integer :: s2f(-5:5,21)=0 !Give shell type & orbital index to get functype
real*8 conv5d6d(6,5),conv7f10f(10,7),conv9g15g(15,9),conv11h21h(21,11)
real*8 conv5d6dtr(5,6),conv7f10ftr(7,10),conv9g15gtr(9,15),conv11h21htr(11,21)
!For backing up spherical basis functions
integer,allocatable :: shelltype5D(:)
real*8,allocatable :: CObasa5D(:,:),CObasb5D(:,:),Sbas5D(:,:),Dbas5D(:,:,:),Magbas5D(:,:,:)
real*8,external :: normgau
ifiletype=1
imodwfn=0
s2f(-5,1:11)=(/ -32,-31,-30,-29,-28,-27,-26,-25,-24,-23,-22 /)
s2f(-4,1:9)=(/ -21,-20,-19,-18,-17,-16,-15,-14,-13 /)
s2f(-3,1:7)=(/ -12,-11,-10,-9,-8,-7,-6 /)
s2f(-2,1:5)=(/ -5,-4,-3,-2,-1 /)
s2f(-1,1:4)=(/ 1,2,3,4 /)
s2f(0,1)=1
s2f(1,1:3)=(/ 2,3,4 /)
s2f(2,1:6)=(/ 5,6,7,8,9,10 /)
s2f(3,1:10)=(/ 11,12,13,17,14,15,18,19,16,20 /) !Note: The sequence of f functions in Multiwfn is not identical to .fch, so convert here. While spdgh are identical
s2f(4,1:15)=(/ 21,22,23,24,25,26,27,28,29,30,31,32,33,34,35 /)
s2f(5,1:21)=(/ 36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56 /)
call gensphcartab(conv5d6d,conv7f10f,conv9g15g,conv11h21h)
conv5d6dtr=transpose(conv5d6d)
conv7f10ftr=transpose(conv7f10f)
conv9g15gtr=transpose(conv9g15g)
conv11h21htr=transpose(conv11h21h)

if (infomode==0.and.ifchprog==2) write(*,*) "Note: This fch file is regarded as produced by Q-Chem!"
open(10,file=name,access="sequential",status="old")
read(10,"(a)") fchtitle
if (infomode==0) write(*,*) "Loading various information of the wavefunction"
isaveNO=0
isaveNBOocc=0
isaveNBOene=0
if (index(fchtitle,'saveNBOocc')/=0.or.index(fchtitle,'SaveNBOocc')/=0) isaveNBOocc=1
if (index(fchtitle,'saveNBOene')/=0.or.index(fchtitle,'SaveNBOene')/=0) isaveNBOene=1
if (index(fchtitle,'saveNO')/=0.or.index(fchtitle,'SaveNO')/=0) isaveNO=1
if ((isaveNBOocc==1.or.isaveNBOene==1).and.infomode==0) write(*,*) "The file contains NBO information"
if (isaveNO==1.and.infomode==0) write(*,*) "The file contains natural orbitals information"
read(10,"(a)") c80
if (c80(11:11)=="R") wfntype=0
if (c80(11:11)=="U") wfntype=1
if (c80(11:12)=="RO") wfntype=2
call loclabel(10,'Number of electrons')
read(10,"(49x,f12.0)") nelec
read(10,"(49x,f12.0)") naelec
read(10,"(49x,f12.0)") nbelec
call loclabel(10,'Number of basis functions')
read(10,"(49x,i12)") nbasis
nindbasis=nbasis
call loclabel(10,'Number of independant functions',ifound) !G03
if (ifound==0) call loclabel(10,'Number of independent functions',ifound) !G09
if (ifound==1) read(10,"(49x,i12)") nindbasis !Number of linear independant functions
ipostHF=0
call loclabel(10,"Post-SCF wavefunction",ipostHF) !If this sentence appears, regardless if "density" is used, the task of the .fch must be post-HF
if (isaveNO==1) ipostHF=1 !When NOs are stored, "Post-SCF wavefunction" may not be present, so set ipostHF=1 explicltly in this case
ibeta=0
call loclabel(10,"Beta MO coefficients",ibeta) !Determine if this is unrestricted system
if (isaveNBOocc==1.or.isaveNBOene==1.or.isaveNO==1.or.ipostHF==1) then
    if (ibeta==0) wfntype=3
    if (ibeta==1) wfntype=4
end if

virialratio=2D0
call loclabel(10,'Virial Ratio',ifound)
if (ifound==1) read(10,"(49x,1PE22.15)") virialratio
totenergy=0.0D0
call loclabel(10,'Total Energy',ifound) !if no this entry, loclabel return ifound=0, else =1
if (ifound==1) read(10,"(49x,1PE22.15)") totenergy
call loclabel(10,'Atomic numbers')
read(10,"(49x,i12)") ncenter
if (allocated(a)) deallocate(a)
allocate(a(ncenter))
read(10,"(6i12)") (a(i)%index,i=1,ncenter)
a%name=ind2name(a%index)
call loclabel(10,'Nuclear charges') !If ECP was used, nuclear charge /= atomic number
read(10,*)
read(10,"(5(1PE16.8))") (a(i)%charge,i=1,ncenter)
call loclabel(10,'Current cartesian coordinates')
read(10,*)
read(10,"(5(1PE16.8))") (a(i)%x,a(i)%y,a(i)%z,i=1,ncenter)
call loclabel(10,'Shell types')
read(10,"(49x,i12)") nshell
allocate(shelltype(nshell))
read(10,"(6i12)") (shelltype(i),i=1,nshell)

!Note that Multiwfn allows cartesian and spherical harmonic basis functions mixed together. If any basis function is spherical harmonic type, then isphergau=1.
!Only the spherical harmonic ones will be treated specially
isphergau=0
if (any(shelltype<=-2)) isphergau=1 
! if (any(shelltype==-2).and.infomode==0) write(*,*) "Note: Spherical d type Gauss functions are detected"
! if (any(shelltype==-3).and.infomode==0) write(*,*) "Note: Spherical f type Gauss functions are detected"
! if (any(shelltype==-4).and.infomode==0) write(*,*) "Note: Spherical g type Gauss functions are detected"
! if (any(shelltype==-5).and.infomode==0) write(*,*) "Note: Spherical h type Gauss functions are detected"
if (any(abs(shelltype)>5).and.infomode==0) then
    write(*,"(a)") " Error: GTFs with angular moment higher h are found, while Multiwfn currently only support up to h. Press Enter to exit"
    read(*,*)
    stop
end if

if (infomode==0) write(*,*) "Loading basis-set definition..."
call loclabel(10,'Number of primitives per shell')
read(10,*)
allocate(shellcon(nshell))
read(10,"(6i12)") (shellcon(i),i=1,nshell)
call loclabel(10,'Shell to atom map')
read(10,*)
allocate(shell2atom(nshell))
read(10,"(6i12)") (shell2atom(i),i=1,nshell)
call loclabel(10,'Primitive exponents')
read(10,"(49x,i12)") nprimshell
allocate(primexp(nprimshell))
read(10,"(5(1PE16.8))") (primexp(i),i=1,nprimshell)
call loclabel(10,'Contraction coefficients')
read(10,*)
allocate(concoeff(nprimshell))
read(10,"(5(1PE16.8))") (concoeff(i),i=1,nprimshell)
read(10,"(a)") c80
if (index(c80,"P(S=P) Contraction coefficients")/=0) then
    backspace(10)
    read(10,*)
    allocate(SPconcoeff(nprimshell))
    read(10,"(5(1PE16.8))") (SPconcoeff(i),i=1,nprimshell)
end if

if (infomode==0) write(*,*) "Loading orbitals..."
call loclabel(10,'Alpha Orbital Energies',ifound)
read(10,*)
!Note: Some basis maybe removed by linear dependence checking, hence the number of orbitals in .fch may less than nbasis(always equals to nmo in Multiwfn)
!Hence when reading information involving the number of orbitals in .fch, use nindbasis instead nmo
!The expansion coefficients, energies in those undefined orbitals are all set to zero
if (wfntype==0.or.wfntype==2.or.wfntype==3) then !Restricted/restricted open-shell, saveNO/NBO
    nmo=nbasis
    allocate(MOene(nmo))
    allocate(MOocc(nmo))
    allocate(MOtype(nmo))
    allocate(amocoeff(nmo,nbasis))
    MOtype=0
    MOocc=0D0
    MOene=0D0
    amocoeff=0D0
    if (wfntype==0.or.wfntype==3) then
        MOocc(1:nint(nelec/2))=2.0D0
    else if (wfntype==2) then
        MOocc(1:nbelec)=2.0D0  !alpha electrons is always more than beta counterpart
        MOocc(nbelec+1:naelec)=1D0
        MOtype(nbelec+1:naelec)=1
    end if
    read(10,"(5(1PE16.8))") (MOene(i),i=1,nindbasis)
    call loclabel(10,'Alpha MO coefficients')
    read(10,*)
    read(10,"(5(1PE16.8))") ((amocoeff(imo,ibasis),ibasis=1,nbasis),imo=1,nindbasis)
else if (wfntype==1.or.wfntype==4) then !unrestricted wavefunction, open-shell post-HF wavefunction
    nmo=2*nbasis
    allocate(MOene(nmo))
    allocate(MOocc(nmo))
    allocate(MOtype(nmo))
    allocate(amocoeff(nbasis,nbasis))
    allocate(bmocoeff(nbasis,nbasis))
    MOocc=0D0
    MOocc(1:naelec)=1D0
    MOocc(nbasis+1:nbasis+nbelec)=1D0
    MOtype(1:nbasis)=1
    MOtype(nbasis+1:nmo)=2
    MOene=0D0
    amocoeff=0D0
    bmocoeff=0D0
    read(10,"(5(1PE16.8))") (MOene(i),i=1,nindbasis)
    call loclabel(10,'Beta Orbital Energies')
    read(10,*)
    read(10,"(5(1PE16.8))") (MOene(i),i=nbasis+1,nbasis+nindbasis)

    call loclabel(10,'Alpha MO coefficients')
    read(10,*)
    read(10,"(5(1PE16.8))") ((amocoeff(imo,ibasis),ibasis=1,nbasis),imo=1,nindbasis)
    call loclabel(10,'Beta MO coefficients')
    read(10,*)
    read(10,"(5(1PE16.8))") ((bmocoeff(imo,ibasis),ibasis=1,nbasis),imo=1,nindbasis)
end if

if (isaveNBOocc==1.or.isaveNO==1) then
    MOocc=MOene
    MOene=0.0D0
end if
if (isaveNBOene==1) MOocc=0.0D0 !For saveNBO, the automatically determined occupation number is meaningless
where (MOocc==1000) MOocc=0.0D0 !When saveNBO is used, the latest several occupation/energy of NBO are 1000, we modify them to zero
where (MOene==1000) MOene=0.0D0

!!!!!! Most reading have finished, now generate basis information

!Backup spherical basis information (some of them may be Cartesian ones) with 5D suffix (of course, may be actually 7f, 9g, 11h...),
!convert them to cartesian type temporarily, at final stage recover them back, namely get Sbas, Ptot... in spherical basis
if (isphergau==1) then
    allocate(shelltype5D(nshell))
    shelltype5D=shelltype
    where (shelltype<=-2) shelltype=-shelltype !Convert to cartesian type
    nbasis5D=nbasis
    nbasis=0
    do i=1,nshell
        nbasis=nbasis+type2norb(shelltype(i))
    end do
end if

!Allocate space for arrays
nprims=0
do i=1,nshell
    nprims=nprims+type2norb(shelltype(i))*shellcon(i)
end do
allocate(b(nprims),co(nmo,nprims),basshell(nbasis),bascen(nbasis),bastype(nbasis),primstart(nbasis),&
primend(nbasis),primconnorm(nprims),basstart(ncenter),basend(ncenter))
!Fill Cobasa and CObasb
if (isphergau==0) then
    allocate(CObasa(nbasis,nbasis))
    CObasa=transpose(amocoeff)
    if (wfntype==1.or.wfntype==4) then
        allocate(CObasb(nbasis,nbasis))
        CObasb=transpose(bmocoeff)
    end if
else if (isphergau==1) then
    allocate(CObasa(nbasis,nbasis),CObasa5D(nbasis5D,nbasis5D))
    CObasa5D=transpose(amocoeff)
    CObasa=0D0
    if (wfntype==1.or.wfntype==4) then
        allocate(CObasb(nbasis,nbasis),CObasb5D(nbasis5D,nbasis5D))
        CObasb5D=transpose(bmocoeff)
        CObasb=0D0
    end if
    !Map 5D coefficient to 6D coefficient
    ipos5D=1
    ipos6D=1
    do ish=1,nshell
        ishtyp5D=shelltype5D(ish)
        ishtyp6D=shelltype(ish)
        numshorb5D=type2norb(ishtyp5D)
        numshorb6D=type2norb(ishtyp6D)
        if (ishtyp5D>=-1) then !S or P or SP or other cartesian shells
            CObasa(ipos6D:ipos6D+numshorb6D-1,1:nbasis5D)=CObasa5D(ipos5D:ipos5D+numshorb5D-1,:)
            if (wfntype==1.or.wfntype==4) CObasb(ipos6D:ipos6D+numshorb6D-1,1:nbasis5D)=CObasb5D(ipos5D:ipos5D+numshorb5D-1,:)            
        else if (ishtyp5D==-2) then
            !5D->6D
            CObasa(ipos6D:ipos6D+5,1:nbasis5D)=matmul(conv5d6d,CObasa5D(ipos5D:ipos5D+4,:))
            if (wfntype==1.or.wfntype==4) CObasb(ipos6D:ipos6D+5,1:nbasis5D)=matmul(conv5d6d,CObasb5D(ipos5D:ipos5D+4,:))
        else if (ishtyp5D==-3) then
            !7F->10F
            CObasa(ipos6D:ipos6D+9,1:nbasis5D)=matmul(conv7f10f,CObasa5D(ipos5D:ipos5D+6,:))
            if (wfntype==1.or.wfntype==4) CObasb(ipos6D:ipos6D+9,1:nbasis5D)=matmul(conv7f10f,CObasb5D(ipos5D:ipos5D+6,:))
        else if (ishtyp5D==-4) then
            !9G->15G
            CObasa(ipos6D:ipos6D+14,1:nbasis5D)=matmul(conv9g15g,CObasa5D(ipos5D:ipos5D+8,:))
            if (wfntype==1.or.wfntype==4) CObasb(ipos6D:ipos6D+14,1:nbasis5D)=matmul(conv9g15g,CObasb5D(ipos5D:ipos5D+8,:))
        else if (ishtyp5D==-5) then
            !11H->21H
            CObasa(ipos6D:ipos6D+20,1:nbasis5D)=matmul(conv11h21h,CObasa5D(ipos5D:ipos5D+10,:))
            if (wfntype==1.or.wfntype==4) CObasb(ipos6D:ipos6D+20,1:nbasis5D)=matmul(conv11h21h,CObasb5D(ipos5D:ipos5D+10,:))
        end if
        ipos5D=ipos5D+numshorb5D
        ipos6D=ipos6D+numshorb6D
    end do
end if

if (infomode==0) write(*,*) "Converting basis function information to GTF information..."
!Distribute exponent, functype to every GTF and generate CO(:,:) from amocoeff/bmocoeff
!Fill: b,basshell,bascen,bastype,co,primstart,primend,primconnorm
k=1 !Current position of GTF
iexp=1
ibasis=1 !Current position of basis
!Note: Below commented with !!! means the line associated to setting basis information
do i=1,nshell !cycle each shell
    b(k:k+shellcon(i)*type2norb(shelltype(i))-1)%center=shell2atom(i)
    basshell(ibasis:ibasis+type2norb(shelltype(i))-1)=i !!! Set basis attributed to which shell
    bascen(ibasis:ibasis+type2norb(shelltype(i))-1)=shell2atom(i) !!! Set basis attributed to which center
    do j=1,type2norb(shelltype(i)) !cycle each basis(orbital) in each shell
        b(k:k+shellcon(i)-1)%functype=s2f(shelltype(i),j)
        bastype(ibasis)=s2f(shelltype(i),j) !!! set basis type
        primstart(ibasis)=k !!! From where the GTFs attributed to ibasis'th basis
        primend(ibasis)=k+shellcon(i)-1 !!! To where the GTFs attributed to ibasis'th basis
        do l=1,shellcon(i) !cycle each GTF in each basis in each shell
            b(k)%exp=primexp(iexp+l-1)
            tnormgau=normgau(b(k)%functype,b(k)%exp)  !!!Normalization coefficient of GTFs
            if (ifchprog==2) tnormgau=1D0 !In the .fch file of Q-chem, normalization factor of GTF has already been multiplied into contraction coefficient of GTFs
            temp=concoeff(iexp+l-1)  !!!Contraction coefficient of GTFs
            if (shelltype(i)==-1.and.j/=1) temp=SPconcoeff(iexp+l-1)
            primconnorm(k)=temp*tnormgau !Combines contraction and normalization coefficient
            do imo=1,nmo
                if (wfntype==0.or.wfntype==2.or.wfntype==3) then
                    co(imo,k)=cobasa(ibasis,imo)*temp*tnormgau
                else if (wfntype==1.or.wfntype==4) then
                    if (isphergau==1) then
                        if (imo<=nbasis5D) co(imo,k)=cobasa(ibasis,imo)*temp*tnormgau
                        if (imo>nbasis5D) co(imo,k)=cobasb(ibasis,imo-nbasis5D)*temp*tnormgau
                    else
                        if (imo<=nbasis) co(imo,k)=cobasa(ibasis,imo)*temp*tnormgau
                        if (imo>nbasis) co(imo,k)=cobasb(ibasis,imo-nbasis)*temp*tnormgau
                    end if
                end if
            end do
            k=k+1
        end do
        ibasis=ibasis+1
    end do
    iexp=iexp+shellcon(i)
end do

!Generate overlap matrix and dipole moment integral matrix for Cartesian Gauss basis functions
if (infomode==0) write(*,*) "Generating overlap matrix..."
allocate(Sbas(nbasis,nbasis))
call genSbas
if (igenDbas==1) then
    if (infomode==0) write(*,*) "Generating electric dipole moment integral matrix..."
    allocate(Dbas(3,nbasis,nbasis))
    call genDbas
end if
if (igenMagbas==1) then
    if (infomode==0) write(*,*) "Generating magnetic dipole moment integral matrix..."
    allocate(Magbas(3,nbasis,nbasis))
    call genMagbas
end if

if (isphergau==1) then
    if (infomode==0) write(*,*) "Converting basis function information from Cartesian to spherical type..."
    !Map cartesian overlap matrix to spherical
    allocate(sbas5D(nbasis5D,nbasis5D))
    if (igenDbas==1) allocate(Dbas5D(3,nbasis5D,nbasis5D))
    if (igenMagbas==1) allocate(Magbas5D(3,nbasis5D,nbasis5D))
    ipos5D=1
    ipos6D=1
    do ish=1,nshell
        ishtyp5D=shelltype5D(ish)
        ishtyp6D=shelltype(ish)
        numshorb5D=type2norb(ishtyp5D)
        numshorb6D=type2norb(ishtyp6D)
        !Now contract columns of Sbas
        if (ishtyp5D>=-1) sbas(:,ipos5D:ipos5D+numshorb5D-1)=sbas(:,ipos6D:ipos6D+numshorb6D-1) !S, P, SP or other Cartesian shells
        if (ishtyp5D==-2) sbas(:,ipos5D:ipos5D+numshorb5D-1)=matmul(sbas(:,ipos6D:ipos6D+numshorb6D-1),conv5d6d)
        if (ishtyp5D==-3) sbas(:,ipos5D:ipos5D+numshorb5D-1)=matmul(sbas(:,ipos6D:ipos6D+numshorb6D-1),conv7f10f)
        if (ishtyp5D==-4) sbas(:,ipos5D:ipos5D+numshorb5D-1)=matmul(sbas(:,ipos6D:ipos6D+numshorb6D-1),conv9g15g)
        if (ishtyp5D==-5) sbas(:,ipos5D:ipos5D+numshorb5D-1)=matmul(sbas(:,ipos6D:ipos6D+numshorb6D-1),conv11h21h)
        !Now contract rows of Sbas
        if (ishtyp5D>=-1) sbas(ipos5D:ipos5D+numshorb5D-1,:)=sbas(ipos6D:ipos6D+numshorb6D-1,:) !S, P, SP or other Cartesian shells
        if (ishtyp5D==-2) sbas(ipos5D:ipos5D+numshorb5D-1,:)=matmul(conv5d6dtr,sbas(ipos6D:ipos6D+numshorb6D-1,:))
        if (ishtyp5D==-3) sbas(ipos5D:ipos5D+numshorb5D-1,:)=matmul(conv7f10ftr,sbas(ipos6D:ipos6D+numshorb6D-1,:))
        if (ishtyp5D==-4) sbas(ipos5D:ipos5D+numshorb5D-1,:)=matmul(conv9g15gtr,sbas(ipos6D:ipos6D+numshorb6D-1,:))
        if (ishtyp5D==-5) sbas(ipos5D:ipos5D+numshorb5D-1,:)=matmul(conv11h21htr,sbas(ipos6D:ipos6D+numshorb6D-1,:))
        
        if (igenDbas==1) then
            do idir=1,3
                !Now contract columns of Dbas
                if (ishtyp5D>=-1) Dbas(idir,:,ipos5D:ipos5D+numshorb5D-1)=Dbas(idir,:,ipos6D:ipos6D+numshorb6D-1) !S, P, SP or other Cartesian shells
                if (ishtyp5D==-2) Dbas(idir,:,ipos5D:ipos5D+numshorb5D-1)=matmul(Dbas(idir,:,ipos6D:ipos6D+numshorb6D-1),conv5d6d)
                if (ishtyp5D==-3) Dbas(idir,:,ipos5D:ipos5D+numshorb5D-1)=matmul(Dbas(idir,:,ipos6D:ipos6D+numshorb6D-1),conv7f10f)
                if (ishtyp5D==-4) Dbas(idir,:,ipos5D:ipos5D+numshorb5D-1)=matmul(Dbas(idir,:,ipos6D:ipos6D+numshorb6D-1),conv9g15g)
                if (ishtyp5D==-5) Dbas(idir,:,ipos5D:ipos5D+numshorb5D-1)=matmul(Dbas(idir,:,ipos6D:ipos6D+numshorb6D-1),conv11h21h)
                !Now contract rows of Dbas
                if (ishtyp5D>=-1) Dbas(idir,ipos5D:ipos5D+numshorb5D-1,:)=Dbas(idir,ipos6D:ipos6D+numshorb6D-1,:) !S, P, SP or other Cartesian shells
                if (ishtyp5D==-2) Dbas(idir,ipos5D:ipos5D+numshorb5D-1,:)=matmul(conv5d6dtr,Dbas(idir,ipos6D:ipos6D+numshorb6D-1,:))
                if (ishtyp5D==-3) Dbas(idir,ipos5D:ipos5D+numshorb5D-1,:)=matmul(conv7f10ftr,Dbas(idir,ipos6D:ipos6D+numshorb6D-1,:))
                if (ishtyp5D==-4) Dbas(idir,ipos5D:ipos5D+numshorb5D-1,:)=matmul(conv9g15gtr,Dbas(idir,ipos6D:ipos6D+numshorb6D-1,:))
                if (ishtyp5D==-5) Dbas(idir,ipos5D:ipos5D+numshorb5D-1,:)=matmul(conv11h21htr,Dbas(idir,ipos6D:ipos6D+numshorb6D-1,:))
            end do
        end if
        if (igenMagbas==1) then
            do idir=1,3
                !Now contract columns of Magbas
                if (ishtyp5D>=-1) Magbas(idir,:,ipos5D:ipos5D+numshorb5D-1)=Magbas(idir,:,ipos6D:ipos6D+numshorb6D-1) !S, P, SP or other Cartesian shells
                if (ishtyp5D==-2) Magbas(idir,:,ipos5D:ipos5D+numshorb5D-1)=matmul(Magbas(idir,:,ipos6D:ipos6D+numshorb6D-1),conv5d6d)
                if (ishtyp5D==-3) Magbas(idir,:,ipos5D:ipos5D+numshorb5D-1)=matmul(Magbas(idir,:,ipos6D:ipos6D+numshorb6D-1),conv7f10f)
                if (ishtyp5D==-4) Magbas(idir,:,ipos5D:ipos5D+numshorb5D-1)=matmul(Magbas(idir,:,ipos6D:ipos6D+numshorb6D-1),conv9g15g)
                if (ishtyp5D==-5) Magbas(idir,:,ipos5D:ipos5D+numshorb5D-1)=matmul(Magbas(idir,:,ipos6D:ipos6D+numshorb6D-1),conv11h21h)
                !Now contract rows of Magbas
                if (ishtyp5D>=-1) Magbas(idir,ipos5D:ipos5D+numshorb5D-1,:)=Magbas(idir,ipos6D:ipos6D+numshorb6D-1,:) !S, P, SP or other Cartesian shells
                if (ishtyp5D==-2) Magbas(idir,ipos5D:ipos5D+numshorb5D-1,:)=matmul(conv5d6dtr,Magbas(idir,ipos6D:ipos6D+numshorb6D-1,:))
                if (ishtyp5D==-3) Magbas(idir,ipos5D:ipos5D+numshorb5D-1,:)=matmul(conv7f10ftr,Magbas(idir,ipos6D:ipos6D+numshorb6D-1,:))
                if (ishtyp5D==-4) Magbas(idir,ipos5D:ipos5D+numshorb5D-1,:)=matmul(conv9g15gtr,Magbas(idir,ipos6D:ipos6D+numshorb6D-1,:))
                if (ishtyp5D==-5) Magbas(idir,ipos5D:ipos5D+numshorb5D-1,:)=matmul(conv11h21htr,Magbas(idir,ipos6D:ipos6D+numshorb6D-1,:))
            end do
        end if
        ipos5D=ipos5D+numshorb5D
        ipos6D=ipos6D+numshorb6D
    end do
    sbas5D=sbas(1:nbasis5D,1:nbasis5D)
    if (igenDbas==1) Dbas5D=Dbas(:,1:nbasis5D,1:nbasis5D)
    if (igenMagbas==1) Magbas5D=Magbas(:,1:nbasis5D,1:nbasis5D)
!     where (abs(sbas5D)<1D-9) sbas5D=0D0 !Ignore too small value to avoid confusing insight
!Test if the sbas generated by Multiwfn is consistent with Gaussian IOP(3/33=1)
!     open(15,file="x\cof3_5d.out",status="old")
!     allocate(tmpmat(nbasis5D,nbasis5D))
!     call loclabel(15,"*** Overlap ***")
!     call readmatgau(15,tmpmat,1)
!     close(15)
!     write(*,*) maxval(abs(tmpmat-sbas5D))

    !Recover spherical Gauss basis function information
    nbasis=nbasis5D
    shelltype=shelltype5D
    ibasis=1
    do i=1,nshell
        basshell(ibasis:ibasis+type2norb(shelltype(i))-1)=i
        bascen(ibasis:ibasis+type2norb(shelltype(i))-1)=shell2atom(i)
        do j=1,type2norb(shelltype(i))
            bastype(ibasis)=s2f(shelltype(i),j)
            ibasis=ibasis+1
        end do
    end do
    deallocate(CObasa)
    allocate(CObasa(nbasis,nbasis))
    CObasa=CObasa5D
    if (wfntype==1.or.wfntype==4) then
        deallocate(CObasb)
        allocate(CObasb(nbasis,nbasis))
        CObasb=CObasb5D
    end if
    deallocate(sbas)
    allocate(sbas(nbasis,nbasis))
    sbas=sbas5D
    if (igenDbas==1) then
        deallocate(Dbas)
        allocate(Dbas(3,nbasis,nbasis))
        Dbas=Dbas5D
    end if
    if (igenMagbas==1) then
        deallocate(Magbas)
        allocate(Magbas(3,nbasis,nbasis))
        Magbas=Magbas5D
    end if
end if

!Split SP shell to S and P shells
noldshell=nshell
noldprimshell=nprimshell
ibasis=1
do i=1,nshell !Count how many basis shells and primitive shells after split SP as S and P, and meantime update basshell
    if (shelltype(i)==-1) then
        nshell=nshell+1
        nprimshell=nprimshell+shellcon(i)
        basshell(ibasis+1:nbasis)=basshell(ibasis+1:nbasis)+1 !The shell index of the basis function after current one should be augmented by 1, since current shell is splitted
    end if
    ibasis=ibasis+type2norb(shelltype(i))
end do
allocate(shtype(nshell),shcen(nshell),shcon(nshell),primshexp(nprimshell),primshcoeff(nprimshell)) !Global array
jsh=1 !New basis shell index
iprsh=1 !Old primitive shell index
jprsh=1 !New primitive shell index
do ish=1,noldshell !Finally determine global arrays shtype,shcen,shcon,primshexp,primshcoeff shell arrays, in which SP shells are not presented
    ncon=shellcon(ish)
    if (shelltype(ish)/=-1) then !Normal shell
        shtype(jsh)=shelltype(ish)
        shcen(jsh)=shell2atom(ish)
        shcon(jsh)=ncon
        primshexp(jprsh:jprsh+ncon-1)=primexp(iprsh:iprsh+ncon-1)
        primshcoeff(jprsh:jprsh+ncon-1)=concoeff(iprsh:iprsh+ncon-1)
        jsh=jsh+1
        jprsh=jprsh+ncon
    else !SP shell
        shtype(jsh)=0 !S
        shtype(jsh+1)=1 !P
        shcen(jsh:jsh+1)=shell2atom(ish)
        shcon(jsh:jsh+1)=ncon
        primshexp(jprsh:jprsh+ncon-1)=primexp(iprsh:iprsh+ncon-1)
        primshexp(jprsh+ncon:jprsh+2*ncon-1)=primexp(iprsh:iprsh+ncon-1)
        primshcoeff(jprsh:jprsh+ncon-1)=concoeff(iprsh:iprsh+ncon-1)
        primshcoeff(jprsh+ncon:jprsh+2*ncon-1)=SPconcoeff(iprsh:iprsh+ncon-1)
        jsh=jsh+2
        jprsh=jprsh+2*ncon
    end if
    iprsh=iprsh+ncon
end do

!Generate basstart and basend
nowcen=0
indcen=0
do ibasis=1,nbasis
    if (bascen(ibasis)/=nowcen) then
        nowcen=bascen(ibasis)
        indcen=indcen+1
        basstart(indcen)=ibasis
        if (indcen/=1) basend(indcen-1)=ibasis-1
    end if
end do
basend(ncenter)=nbasis

!Generate or load one-particle density matrix for basis functions
if (igenP==1) then
    if (ipostHF/=1) then
        if (infomode==0) write(*,*) "Generating SCF density matrix..."
        call genP
    else if (ipostHF==1.and.wfntype==3) then
    !Cannot only use wfntype==3 as criteria directly, when pop=saveNBO, although wfntype==3, but post-HF density matrix is not presented.
    !    For post-HF calculation, if "density" keyword was used, although orbital coefficients are still SCF's,
    !however, post-HF density matrix is presented in .fch file. Read them, then analysis based on basis function
    !will be post-HF's happily. 
    !If not, we must fill MOocc(:) from .wfn or pop=NO result and then generate density matrix, because
    !.fch itself doesn't contain natural orbital occupation information!
        allocate(ptot(nbasis,nbasis))
        Ptot=0D0
        call loclabel(10,"Total SCF Density",ifoundSCFDM) !In rare case, even SCF density is missing!
        if (ifoundSCFDM==1) then !Further find post-HF density matrix section
            read(10,*)
            call loclabel(10,"Total ",ifoundpostDM,0)
        end if
        if (ifoundpostDM==1) then
            if (infomode==0) write(*,*) "Loading post-HF density matrix..."
            read(10,*)
            read(10,"(5(1PE16.8))") ((Ptot(i,j),j=1,i),i=1,nbasis)
            ptot=ptot+transpose(ptot)
            do i=1,nbasis
                ptot(i,i)=ptot(i,i)/2D0
            end do
            if (infomode==0) write(*,"(a)") " Note: Post-HF density matrix is loaded from .fch file, all of the following analyses based on basis information will be for post-HF wavefunction"
        else if (ifoundSCFDM==0.or.ifoundpostDM==0) then
            if (isaveNO==0) then
                if (infomode==0) write(*,"(a)") " Warning: Post-HF density matrix was not found and natural orbitals are seemingly not stored, so all of the following analyses based on basis information will be for HF wavefunction"
                if (infomode==0) write(*,*) "Generating HF density matrix..."
            else if (isaveNO==1) then
                if (infomode==0) write(*,*) "Generating post-HF density matrix based on natural orbitals..."
            end if
            call genP
        end if
    else if (ipostHF==1.and.wfntype==4) then
        allocate(ptot(nbasis,nbasis),palpha(nbasis,nbasis),pbeta(nbasis,nbasis))
        Ptot=0D0
        Palpha=0D0
        Pbeta=0D0
        call loclabel(10,"Total SCF Density",ifoundSCFDM)
        if (ifoundSCFDM==1) then !Further find post-HF density matrix section
            read(10,*)
            call loclabel(10,"Total ",ifoundpostDM,0)
        end if
        if (ifoundpostDM==1) then
            if (infomode==0) write(*,*) "Loading post-HF density matrix..."
            read(10,*)
            read(10,"(5(1PE16.8))") ((Ptot(i,j),j=1,i),i=1,nbasis)
            read(10,*)
            read(10,"(5(1PE16.8))") ((Palpha(i,j),j=1,i),i=1,nbasis) !Read spin density matrix, use Palpha as temporary store
            Palpha=(Palpha+Ptot)/2D0
            Pbeta=Ptot-Palpha
            ptot=ptot+transpose(ptot)
            Palpha=Palpha+transpose(Palpha)
            Pbeta=Pbeta+transpose(Pbeta)
            do i=1,nbasis
                ptot(i,i)=ptot(i,i)/2D0
                Palpha(i,i)=Palpha(i,i)/2D0
                Pbeta(i,i)=Pbeta(i,i)/2D0
            end do
            if (infomode==0) write(*,"(a)") " Note: Post-HF density matrix is loaded from .fch file, all of the following analyses based on basis information will be for post-HF wavefunction"
        else if (ifoundSCFDM==0.or.ifoundpostDM==0) then
            if (isaveNO==0) then
                if (infomode==0) write(*,"(a)") " Warning: Post-HF density matrix was not found and natural orbitals are seemingly not stored, so all of the following analyses based on basis information will be for HF wavefunction"
                if (infomode==0) write(*,*) "Generating HF density matrix..."
            else if (isaveNO==1) then
                if (infomode==0) write(*,*) "Generating post-HF density matrix based on natural orbitals..."
            end if
            call genP
        end if
    end if
end if

!Output summary of present wavefunction
if (infomode==0) then
    write(*,*)
    write(*,"(' Total/Alpha/Beta electrons:',3f12.4)") nelec,naelec,nbelec
    write(*,"(' Net charge:',f12.5,'      Expected multiplicity:',i5)") sum(a(:)%charge)-nelec,nint(naelec-nbelec)+1
    write(*,"(' Atoms:',i7,',  Basis functions:',i7,',  GTFs:',i7)") ncenter,nbasis,nprims
    write(*,"(' Total energy:',f19.12,' Hartree,   Virial ratio:',f12.8)") totenergy,virialratio
    if (wfntype==0) then
        write(*,"(' This is a restricted single-determinant wavefunction')")
        write(*,"(' Orbitals from 1 to',i6,' are occupied')") nint(nelec/2)
    else if (wfntype==1) then
        write(*,"(' This is an unrestricted single-determinant wavefunction')")
        write(*,"(' Orbitals from ',i6,' to',i6,' are alpha, from',i6,' to',i6,' are occupied')") 1,nbasis,1,nint(naelec)
        write(*,"(' Orbitals from ',i6,' to',i6,' are beta,  from',i6,' to',i6,' are occupied')") nbasis+1,nmo,nbasis+1,nbasis+nint(nbelec)
    else if (wfntype==2) then
        write(*,"(' This is a restricted open-shell wavefunction')")
        write(*,"(' Orbitals from',i6,' to',i6,' are doubly occupied')") 1,nint(nbelec)
        write(*,"(' Orbitals from',i6,' to',i6,' are singly occupied')") nint(nbelec)+1,nint(naelec)
    else if (wfntype==3) then
        write(*,"(' This is a restricted post-HF wavefunction')")
    else if (wfntype==4) then
        write(*,"(' This is an unrestricted post-HF wavefunction')")
        write(*,"(' Orbitals from ',i6,' to',i6,' are alpha, from',i6,' to',i6,' are beta')") 1,nbasis,nbasis+1,nmo
    end if
    write(*,"(' Title line of this file: ',a)") trim(fchtitle)
end if

!Find out index of HOMO, will be used in some cases, only for RHF
if (wfntype==0) then
    do idxHOMO=nmo,1,-1
        if (nint(MOocc(idxHOMO))==2D0) exit
    end do
end if

!Clean
close(10)
deallocate(shelltype,shellcon,shell2atom,primexp,concoeff,amocoeff)
if (allocated(SPconcoeff)) deallocate(SPconcoeff)
if (allocated(bmocoeff)) deallocate(bmocoeff)
if (allocated(tmpmat)) deallocate(tmpmat)
if (allocated(tmparr)) deallocate(tmparr)
if (isphergau==1) then
    deallocate(shelltype5D,cobasa5D)
    if (allocated(cobasb5D)) deallocate(cobasb5D)
end if
end subroutine



!!-----------------------------------------------------------------
!!------------------------- Read gaussian formatted check file for AdNDP analysis
!Cobasa is not read from .fch file, which has already written by AdNDP module, which contains AO expansion coefficients of candidate or saved orbitals
subroutine readfchadndp(fchfilename,iusespin,orbocc,adndpCObas,numorb)
use defvar
use util
implicit real*8 (a-h,o-z)
character*80 fchfilename,c80
integer iusespin
real*8 orbocc(numorb),adndpCObas(nbasis,numorb)
integer,allocatable :: shelltype(:),shell2atom(:),shellcon(:) !Degree of shell contraction
real*8,allocatable :: primexp(:),concoeff(:),SPconcoeff(:)
integer :: s2f(-5:5,21)=0 !Give shell type & orbital index to get functype
real*8 conv5d6d(6,5),conv7f10f(10,7),conv9g15g(15,9),conv11h21h(21,11)
 !For backing up spherical basis functions
integer,allocatable :: shelltype5D(:)
real*8,allocatable :: CObasa5D(:,:)
real*8,external :: normgau
nmo=numorb
s2f(-5,1:11)=(/ -32,-31,-30,-29,-28,-27,-26,-25,-24,-23,-22 /)
s2f(-4,1:9)=(/ -21,-20,-19,-18,-17,-16,-15,-14,-13 /)
s2f(-3,1:7)=(/ -12,-11,-10,-9,-8,-7,-6 /)
s2f(-2,1:5)=(/ -5,-4,-3,-2,-1 /)
s2f(-1,1:4)=(/ 1,2,3,4 /)
s2f(0,1)=1
s2f(1,1:3)=(/ 2,3,4 /)
s2f(2,1:6)=(/ 5,6,7,8,9,10 /)
s2f(3,1:10)=(/ 11,12,13,17,14,15,18,19,16,20 /)
s2f(4,1:15)=(/ 21,22,23,24,25,26,27,28,29,30,31,32,33,34,35 /)
s2f(5,1:21)=(/ 36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56 /)
call gensphcartab(conv5d6d,conv7f10f,conv9g15g,conv11h21h)

open(10,file=fchfilename,access="sequential",status="old")

! if (.not.allocated(a)) then !Hasn't loaded some information before
!     call loclabel(10,'Atomic numbers')
!     read(10,"(49x,i12)") ncenter
!     allocate(a(ncenter))
!     read(10,"(6i12)") (a(i)%index,i=1,ncenter)
!     a%name=ind2name(a%index)
    call loclabel(10,'Nuclear charges') !If ECP was used, nuclear charge /= atomic number
    read(10,*)
    read(10,"(5(1PE16.8))") (a(i)%charge,i=1,ncenter)
    call loclabel(10,'Current cartesian coordinates')
    read(10,*)
    read(10,"(5(1PE16.8))") (a(i)%x,a(i)%y,a(i)%z,i=1,ncenter)
! end if

call loclabel(10,'Shell types')
read(10,"(49x,i12)") nshell
allocate(shelltype(nshell))
read(10,"(6i12)") (shelltype(i),i=1,nshell)
isphergau=0
if (any(shelltype<=-2)) isphergau=1
call loclabel(10,'Number of primitives per shell')
read(10,*)
allocate(shellcon(nshell))
read(10,"(6i12)") (shellcon(i),i=1,nshell)
call loclabel(10,'Shell to atom map')
read(10,*)
allocate(shell2atom(nshell))
read(10,"(6i12)") (shell2atom(i),i=1,nshell)
call loclabel(10,'Primitive exponents')
read(10,"(49x,i12)") nprimshell
allocate(primexp(nprimshell))
read(10,"(5(1PE16.8))") (primexp(i),i=1,nprimshell)
call loclabel(10,'Contraction coefficients')
read(10,*)
allocate(concoeff(nprimshell))
read(10,"(5(1PE16.8))") (concoeff(i),i=1,nprimshell)
read(10,"(a)") c80
if (index(c80,"P(S=P) Contraction coefficients")/=0) then
    backspace(10)
    read(10,*)
    allocate(SPconcoeff(nprimshell))
    read(10,"(5(1PE16.8))") (SPconcoeff(i),i=1,nprimshell)
end if

if (allocated(MOene)) deallocate(MOene,MOocc,MOtype)
allocate(MOene(nmo),MOocc(nmo),MOtype(nmo))
MOene=0D0
MOocc=orbocc
MOtype=iusespin

if (isphergau==1) then
    allocate(shelltype5D(nshell))
    shelltype5D=shelltype
    where (shelltype<=-2) shelltype=-shelltype !Convert to cartesian type
    nbasis5D=nbasis
    nbasis=0
    do i=1,nshell
        nbasis=nbasis+type2norb(shelltype(i))
    end do
end if

!Allocate space for arrays
nprims=0
do i=1,nshell
    nprims=nprims+type2norb(shelltype(i))*shellcon(i)
end do
if (.not.allocated(b)) allocate(b(nprims))
if (allocated(CO)) deallocate(CO)
allocate(CO(nmo,nprims))

if (isphergau==0) then
    if (allocated(CObasa)) deallocate(CObasa)
    allocate(CObasa(nbasis,nmo))
    CObasa=adndpCObas
else if (isphergau==1) then
    if (allocated(CObasa)) deallocate(CObasa)
    allocate(CObasa(nbasis,nmo),CObasa5D(nbasis5D,nmo))
    CObasa5D=adndpCObas
    CObasa=0D0
    
    !Map 5D coefficient to 6D coefficient
    ipos5D=1
    ipos6D=1
    do ish=1,nshell
        ishtyp5D=shelltype5D(ish)
        ishtyp6D=shelltype(ish)
        numshorb5D=type2norb(ishtyp5D)
        numshorb6D=type2norb(ishtyp6D)
        if (ishtyp5D==0.or.ishtyp5D==1.or.ishtyp5D==-1) then !S or P or SP
            CObasa(ipos6D:ipos6D+numshorb6D-1,1:nmo)=CObasa5D(ipos5D:ipos5D+numshorb5D-1,1:nmo)
        else if (ishtyp5D==-2) then
            !5D->6D
            CObasa(ipos6D:ipos6D+5,1:nmo)=matmul(conv5d6d,CObasa5D(ipos5D:ipos5D+4,1:nmo))
        else if (ishtyp5D==-3) then
            !7F->10F
            CObasa(ipos6D:ipos6D+9,1:nmo)=matmul(conv7f10f,CObasa5D(ipos5D:ipos5D+6,1:nmo))
        else if (ishtyp5D==-4) then
            !9G->15G
            CObasa(ipos6D:ipos6D+14,1:nmo)=matmul(conv9g15g,CObasa5D(ipos5D:ipos5D+8,1:nmo))
        else if (ishtyp5D==-5) then
            !11H->21H
            CObasa(ipos6D:ipos6D+20,1:nmo)=matmul(conv11h21h,CObasa5D(ipos5D:ipos5D+10,1:nmo))
        end if
        ipos5D=ipos5D+numshorb5D
        ipos6D=ipos6D+numshorb6D
    end do
end if

!Distribute exponent, functype to every GTF and generate CO(:,:) from CObasa
k=1 !current position of GTF
iexp=1
ibasis=1 !current position of basis
!Note: Below commented with !!! means the line associated to setting basis information
do i=1,nshell !cycle each shell
    b(k:k+shellcon(i)*type2norb(shelltype(i))-1)%center=shell2atom(i)
    do j=1,type2norb(shelltype(i)) !cycle each basis(orbital) in each shell
        b(k:k+shellcon(i)-1)%functype=s2f(shelltype(i),j)
        do l=1,shellcon(i) !cycle each GTF in each basis in each shell
            b(k)%exp=primexp(iexp+l-1)
            tnormgau=normgau(b(k)%functype,b(k)%exp)
            temp=concoeff(iexp+l-1)
            if (shelltype(i)==-1.and.j/=1) temp=SPconcoeff(iexp+l-1)
            do imo=1,nmo
                co(imo,k)=cobasa(ibasis,imo)*temp*tnormgau
            end do
            k=k+1
        end do
        ibasis=ibasis+1
    end do
    iexp=iexp+shellcon(i)
end do

close(10)
end subroutine



!!-----------------------------------------------------------------
!!------------- Read .chg file that only contain atomic charge information
subroutine readchg(name,infomode) !infomode=0 means output info, =1 silent
use defvar
use util
implicit real*8 (a-h,o-z)
character(len=*) name
integer i
real*8 dipx,dipy,dipz,tdip
ifiletype=4
open(10,file=name,access="sequential",status="old")
ncenter=totlinenum(10,1)
allocate(a(ncenter))
dipx=0.0D0
dipy=0.0D0
dipz=0.0D0
do i=1,ncenter
    read(10,*) a(i)%name,a(i)%x,a(i)%y,a(i)%z,a(i)%charge
    dipx=dipx+a(i)%x*a(i)%charge
    dipy=dipy+a(i)%y*a(i)%charge
    dipz=dipz+a(i)%z*a(i)%charge
    do j=0,nelesupp
        if (a(i)%name==ind2name(j)) then
            a(i)%index=j
            exit
        end if
        if (j==nelesupp) write(*,*) "Warning: Found unknown element!"
    end do
end do
close(10)
a%x=a%x/b2a
a%y=a%y/b2a
a%z=a%z/b2a
if (infomode==0) then
    write(*,"(' Total',i8,' atoms')") ncenter
    write(*,"(' Summing up charges:',f12.6)") sum(a%charge)
    write(*,"(' Component of dipole in X/Y/Z:',3f12.6,' a.u.')") dipx,dipy,dipz
    tdip=dsqrt(dipx**2+dipy**2+dipz**2)
    write(*,"(' Total dipole:',f12.6,' a.u., equivalent to',f12.6,' Debye')") tdip,tdip*au2debye
end if
end subroutine



!!-----------------------------------------------------------------
!!----- Read .pdb file, for visualization and calculate RDG with promolecular approximation
subroutine readpdb(name,infomode) !mode=0 means output wfn property,=1 not
use defvar
implicit real*8 (a-h,o-z)
character(len=*) name
character test*6,tmpname*4,element*3
integer i,j
ifiletype=5
open(10,file=name,access="sequential",status="old")
ncenter=0
do while(.true.)
    read(10,"(a6)",iostat=ierror) test
    if (ierror/=0) exit
    if (test=="HETATM".OR.test=="ATOM  ") ncenter=ncenter+1
end do
rewind(10)
allocate(a(ncenter))
i=0
do while(.true.)
    read(10,"(a6)",iostat=ierror) test
    if (ierror/=0) exit
    if (test=="HETATM".or.test=="ATOM  ") then
        backspace(10)
        i=i+1
        read(10,"(12x,a4,14x,3f8.3,22x,a3)") tmpname,a(i)%x,a(i)%y,a(i)%z,element
        tmpname=adjustl(tmpname)
        element=adjustl(element)
        ifound=0
        do j=1,nelesupp
            !Use "element" term to determine actual element
            if (ind2name_up(j)==element(1:2).or.ind2name(j)==element(1:2)) then
                a(i)%index=j
                ifound=1
                exit
            end if
        end do
        if (ifound==0) then
            !"Element" term is missing, use first character of atomic name to determine element
            !Check for name such as C5,N11,O3B, also check name such as 1H5* (at this time, the first letter much be digital)
            do j=1,nelesupp
                if (ind2name_up(j)==tmpname(1:1)//' '.or.&
                (ichar(tmpname(1:1))<=57).and.ind2name_up(j)==tmpname(2:2)//' ') then
                    a(i)%index=j
                    ifound=1
                    exit
                end if
            end do
        end if
        if (ifound==0) then
            write(*,"(3a)") "Warning: Found unknown element """,tmpname,""" , assume it is carbon"
            a(i)%index=12
        end if
        a(i)%name=ind2name(a(i)%index)
    end if
end do
close(10)
a%x=a%x/b2a
a%y=a%y/b2a
a%z=a%z/b2a
a%charge=a%index
if (infomode==0) write(*,"(' Total',i8,' atoms')") ncenter
end subroutine



!!-----------------------------------------------------------------
!!----- Read .xyz file, for visualization and calculate RDG with promolecular approximation
!infomode=0 means output wfn property, =1 not
!iopen=0 means don't open and close file in this routine, used to continuously read trjectory. =1 means do these
subroutine readxyz(name,infomode,iopen) 
use defvar
use util
implicit real*8 (a-h,o-z)
integer infomode,iopen
character(len=*) name
character*79 titleline
ifiletype=5
if (iopen==1) open(10,file=name,access="sequential",status="old")
ncenter=0
read(10,*) ncenter
read(10,"(a)") titleline
if (allocated(a)) deallocate(a)
allocate(a(ncenter))
do i=1,ncenter
    read(10,*) a(i)%name,a(i)%x,a(i)%y,a(i)%z
    call lc2uc(a(i)%name(1:1)) !Convert to upper case
    call uc2lc(a(i)%name(2:2)) !Convert to lower case
    do j=1,nelesupp
        if ( a(i)%name==ind2name(j) ) then
            a(i)%index=j
            exit
        end if
    end do
    if (j==nelesupp+1) then !Only use the first letter of atom name to try to assign. For example OW, HW may be in .xyz file
        do j=1,nelesupp
            if ( ind2name(j)(2:2)==" ".and.a(i)%name(1:1)==ind2name(j)(1:1) ) then
                a(i)%index=j
                exit
            end if
        end do
    end if
    if (j==nelesupp+1) then
        write(*,"(3a)") "Warning: Found unknown element """,a(i)%name,""" , assume it is carbon"
        a(i)%index=12
        write(*,*) "Press ENTER to continue"
        read (*,*)
    end if
end do
if (iopen==1) close(10)
a%x=a%x/b2a
a%y=a%y/b2a
a%z=a%z/b2a
a%charge=a%index
if (infomode==0) then
    write(*,"(a)") titleline
    write(*,"(' Total',i8,' atoms')") ncenter
end if
end subroutine



!!-----------------------------------------------------------------
!!------------------------- Read .31 and one of .32 to .40 file generated by NBO program
subroutine read31(name,infomode) !mode=0 means output wfn property,=1 not
use defvar
use util
implicit real*8 (a-h,o-z)
character(len=*) name
character :: name2*200=" ",chartemp*80,tmpc2*2
integer infomode
integer bastype2func(500) !Convert basis type in NBO(NBO5.0 Program Manual P103) to function type in .wfn
integer,allocatable :: shellcon(:),shell2atom(:),shellnumbas(:),shell2prmshl(:),bastype31(:),shellnumbas5D(:)
real*8,allocatable :: orbcoeff(:,:),prmshlexp(:),cs(:),cp(:),cd(:),cf(:),cg(:),orbcoeff5D(:,:)
real*8 :: conv5d6d(6,5)=0D0,conv7f10f(10,7)=0D0,conv9g15g(15,9)=0D0
real*8,external :: normgau
ifiletype=6
bastype2func(1)=1 !s
bastype2func(101)=2 !x
bastype2func(102)=3 !y
bastype2func(103)=4 !z
bastype2func(201)=5 !xx
bastype2func(202)=8 !xy
bastype2func(203)=9 !xz
bastype2func(204)=6 !yy
bastype2func(205)=10 !yz
bastype2func(206)=7 !zz
bastype2func(301)=11 !xxx
bastype2func(302)=14 !xxy
bastype2func(303)=15 !xxz
bastype2func(304)=17 !xyy
bastype2func(305)=20 !xyz
bastype2func(306)=18 !xzz
bastype2func(307)=12 !yyy
bastype2func(308)=16 !yyz
bastype2func(309)=19 !yzz
bastype2func(310)=13 !zzz
!Below sequence comes from line 47384 in NBO_5 src
bastype2func(401)=35 !XXXX
bastype2func(402)=34 !XXXY
bastype2func(403)=33 !XXXZ
bastype2func(404)=32 !XXYY
bastype2func(405)=31 !XXYZ
bastype2func(406)=30 !XXZZ
bastype2func(407)=29 !XYYY
bastype2func(408)=28 !XYYZ
bastype2func(409)=27 !XYZZ
bastype2func(410)=26 !XZZZ
bastype2func(411)=25 !YYYY
bastype2func(412)=24 !YYYZ
bastype2func(413)=23 !YYZZ
bastype2func(414)=22 !YZZZ
bastype2func(415)=21 !ZZZZ
!Used to convert coefficient matrix from 5D to 6D
!5D sequence in .31: 255   252   253   254   251
!namely -0.5*XX-0.5*YY+ZZ, XZ, YZ, ¡Ì3/2*(XX-YY), XY
!to 6D: XX,XY,XZ,YY,YZ,ZZ, namely 201~206, the indexes are consecutive
conv5d6d(1,1)=-0.5D0
conv5d6d(4,1)=-0.5D0
conv5d6d(6,1)=1D0
conv5d6d(3,2)=1D0
conv5d6d(5,3)=1D0
conv5d6d(1,4)=sqrt(3D0)/2D0
conv5d6d(4,4)=-sqrt(3D0)/2D0
conv5d6d(2,5)=1D0
!Used to convert coefficient matrix from 7F to 10F, Standard f set
!7F sequence in .31: 351   352   353   354   355   356   357
!to 10F: XXX,XXY,XXZ,XYY,XYZ,XZZ,YYY,YYZ,YZZ,ZZZ, namely 301~310, the indexes are consecutive
conv7f10f(3,1)=-0.721962098225322D0
conv7f10f(8,1)=-0.721962098225322D0
conv7f10f(10,1)=0.481308065483548D0
conv7f10f(1,2)=-0.281160203343101D0
conv7f10f(4,2)=-0.281160203343101D0
conv7f10f(6,2)=1.1246408133724D0
conv7f10f(2,3)=-0.281160203343101D0
conv7f10f(7,3)=-0.281160203343101D0
conv7f10f(9,3)=1.1246408133724D0
conv7f10f(3,4)=sqrt(3D0)/2D0
conv7f10f(8,4)=-sqrt(3D0)/2D0
conv7f10f(5,5)=1D0
conv7f10f(1,6)=0.369693511996758D0
conv7f10f(4,6)=-1.10908053599027D0
conv7f10f(2,7)=1.10908053599027D0
conv7f10f(7,7)=-0.369693511996758D0

write(*,*) "Input filename with suffix ranging from .32 to .40 (e.g. ltwd.35)"
write(*,*) "Note: 32=PNAO 33=NAO 34=PNHO 35=NHO 36=PNBO 37=NBO 38=PNLMO 39=NLMO 40=MO"
do while(.true.)
    read(*,"(a)") name2
    if (name2(1:2)=="32".or.name2(1:2)=="33".or.name2(1:2)=="34".or.name2(1:2)=="35"&
    .or.name2(1:2)=="36".or.name2(1:2)=="37".or.name2(1:2)=="38".or.name2(1:2)=="39".or.name2(1:2)=="40") then
        tmpc2=name2(1:2)
        name2(1:len(name))=name
        name2(len_trim(name2)-1:len_trim(name2))=tmpc2
    end if
    inquire(file=name2,exist=alive)
    if (alive.eqv..true.) exit
    write(*,*) "File not found, input again"
end do

itmplen=len_trim(name2)
if (name2(itmplen-1:itmplen)=="32") write(*,*) "Loading .32 file(PNAO)"
if (name2(itmplen-1:itmplen)=="33") write(*,*) "Loading .33 file(NAO)"
if (name2(itmplen-1:itmplen)=="34") write(*,*) "Loading .34 file(PNHO)"
if (name2(itmplen-1:itmplen)=="35") write(*,*) "Loading .35 file(NHO)"
if (name2(itmplen-1:itmplen)=="36") write(*,*) "Loading .36 file(PNBO)"
if (name2(itmplen-1:itmplen)=="37") write(*,*) "Loading .37 file(NBO)"
if (name2(itmplen-1:itmplen)=="38") write(*,*) "Loading .38 file(PNLMO)"
if (name2(itmplen-1:itmplen)=="39") write(*,*) "Loading .39 file(NLMO)"
if (name2(itmplen-1:itmplen)=="40") write(*,*) "Loading .40 file(MO)"
open(10,file=name,access="sequential",status="old")
read(10,*)
read(10,*)
read(10,*)
read(10,*) ncenter,nshell,nprimshell
allocate(a(ncenter),shellcon(nshell),shell2atom(nshell),shellnumbas(nshell),shell2prmshl(nshell))
allocate(prmshlexp(nprimshell),cs(nprimshell),cp(nprimshell),cd(nprimshell),cf(nprimshell),cg(nprimshell))
allocate(bastype31(nshell*15)) !We don't know how many basis before read the file, so use the maximum value(up to g function)
read(10,*)
do i=1,ncenter
    read(10,*) a(i)%index,a(i)%x,a(i)%y,a(i)%z
    a(i)%charge=a(i)%index !.31 doesn't record charge and index separately, so we have to make them indentical
end do
a%name=ind2name(a%index)
a%x=a%x/b2a
a%y=a%y/b2a
a%z=a%z/b2a
read(10,*)
j=1
do i=1,nshell
    read(10,*) shell2atom(i),shellnumbas(i),shell2prmshl(i),shellcon(i)
    read(10,*) bastype31(j:j+shellnumbas(i)-1)
    j=j+shellnumbas(i)
end do
isphergau=0
if (any(bastype31==251).or.any(bastype31==351).or.any(bastype31==451)) isphergau=1
!The conversion relationship between Cartesian and pure functions used by NBO is different to mainstream quantum chemistry code, the relationship for f is documented
!in the NBO manual, however there is no way to find out that for g-type. So if g-type is involved, Cartesian type must be used.
if (any(bastype31==451)) then
    write(*,"(a)") "Error: Multiwfn doesn't support spherical harmonic Gauss functions with g or higher angular moment for NBO plot files. For Gaussian, &
    you should add ""6d 10f"" keywords and regenerate these files. Press ENTER to exit"
    read(*,*)
    stop
end if
read(10,*)
read(10,*) prmshlexp
read(10,*)
read(10,*) cs
read(10,*)
read(10,*) cp
read(10,*)
read(10,*) cd
read(10,*)
read(10,*) cf
read(10,"(a)",iostat=ierror) chartemp
if (ierror==0) then
    backspace(10)
    read(10,*) cg
end if
close(10)

totenergy=0D0
virialratio=0D0
nbasis=sum(shellnumbas) !Calculate the number of basis
open(10,file=name2,access="sequential",status="old")
read(10,*)
read(10,*)
read(10,*)
read(10,"(a)") chartemp
naelec=0D0
nbelec=0D0
if (chartemp(1:11)==" ALPHA SPIN".or.chartemp(1:11)==" alpha spin") then
    wfntype=4
    nmo=2*nbasis
    allocate(orbcoeff(nbasis,nmo),MOocc(nmo),MOene(nmo),MOtype(nmo))
    MOtype(1:nbasis)=1 !alpha
    MOtype(nbasis+1:nmo)=2 !beta
    do iorb=1,nbasis
        read(10,*) orbcoeff(1:nbasis,iorb)
    end do
    if (name2(itmplen-1:itmplen)=="37".or.name2(itmplen-1:itmplen)=="39") read(10,*) MOocc(1:nbasis)
    
    if (chartemp(1:11)==" ALPHA SPIN") then
        call loclabel(10," BETA  SPIN",ifound)
    else
        call loclabel(10," beta  spin",ifound)
    end if
    read(10,*)
    do iorb=nbasis+1,nmo
        read(10,*) orbcoeff(1:nbasis,iorb)
    end do
    if (name2(itmplen-1:itmplen)=="37".or.name2(itmplen-1:itmplen)=="39") read(10,*) MOocc(nbasis+1:nmo)
    
    if (name2(itmplen-1:itmplen)/="37".and.name2(itmplen-1:itmplen)/="39") MOocc=1D0 !This setting has no physical meaning
    naelec=sum(MOocc(1:nbasis))
    nbelec=sum(MOocc(nbasis+1:nmo))
    nelec=naelec+nbelec
else !Close shell system
    wfntype=3
    nmo=nbasis
    allocate(orbcoeff(nbasis,nmo),MOocc(nmo),MOene(nmo),MOtype(nmo))
    MOtype=0
    orbcoeff=0
    MOocc=1D0 !This setting has no physical meaning
    backspace(10)
    do iorb=1,nmo
        read(10,*,iostat=ierror) orbcoeff(1:nbasis,iorb)
        if (ierror/=0) then !Note that the number of orbitals may be smaller than nbasis, therefore if we don't do this determination, we may encounter end-of-file!
            write(*,"(' Warning: The number of orbitals is smaller than basis functions!',/)")
            goto 1300 !I am unable to determine how many actually the number of orbitals it is, if this condition really happen, &
            ! the occupation will be incorrectly (but unavoidable) loaded as coefficients, and thus we skip occupation loading
        end if
    end do
    if (name2(itmplen-1:itmplen)=="37".or.name2(itmplen-1:itmplen)=="39") read(10,*) MOocc
1300    nelec=sum(MOocc)
    naelec=nelec/2
    nbelec=naelec
end if
close(10)
MOene=0D0

!Temporarily convert spherical harmonic gauss functions' information to cartesian type
if (isphergau==1) then
    !Calculate how many cartesian basis after conversion
    nbasis5D=nbasis
    nbasis=0
    do i=1,nshell
        if (shellnumbas(i)==1.or.shellnumbas(i)==3.or.shellnumbas(i)==4) nbasis=nbasis+shellnumbas(i) !S,P,SP
        if (shellnumbas(i)==5) nbasis=nbasis+6 !D
        if (shellnumbas(i)==7) nbasis=nbasis+10 !F
        if (shellnumbas(i)==9) nbasis=nbasis+15 !G
    end do
    
    allocate(shellnumbas5D(nshell))
    shellnumbas5D=shellnumbas !Backup
    where(shellnumbas==5) shellnumbas=6 !Convert number of orbitals in each shell from 5D to cartesian type
    where(shellnumbas==7) shellnumbas=10
    where(shellnumbas==9) shellnumbas=15

    allocate(orbcoeff5D(nbasis5D,nmo))
    orbcoeff5D=orbcoeff !Backup
    deallocate(orbcoeff)
    allocate(orbcoeff(nbasis,nmo)) !Enlarge size from spherical type to cartesian type
    orbcoeff=0D0

    deallocate(bastype31)
    allocate(bastype31(nbasis)) !Enlarge size from spherical type to cartesian type

    !Generate cartesian .31 basis type, the indexes are consecutive, in line with conv5d6d and conv7f10f
    i=1
    do ish=1,nshell
        if (shellnumbas(ish)==1) then
            bastype31(i)=1
        else if (shellnumbas(ish)==3) then
            bastype31(i)=101
            bastype31(i+1)=102
            bastype31(i+2)=103
        else if (shellnumbas(ish)==4) then
            bastype31(i)=1
            bastype31(i+1)=101
            bastype31(i+2)=102
            bastype31(i+3)=103
        else if (shellnumbas(ish)==6) then !d
            do j=1,6
                bastype31(i+j-1)=200+j
            end do
        else if (shellnumbas(ish)==10) then !f
            do j=1,10
                bastype31(i+j-1)=300+j
            end do
        else if (shellnumbas(ish)==15) then !g
            do j=1,15
                bastype31(i+j-1)=400+j
            end do
        end if
        i=i+shellnumbas(ish)
    end do
    
    !Map 5D coefficient to 6D coefficient
    ipos5D=1
    ipos6D=1
    do ish=1,nshell
        n5D=shellnumbas5D(ish)
        n6D=shellnumbas(ish)
        if (n5D==1.or.n5D==3.or.n5D==4) then !S or P or SP
            if (wfntype==3) orbcoeff(ipos6D:ipos6D+n6D-1,1:nbasis5D)=orbcoeff5D(ipos5D:ipos5D+n5D-1,:)
            if (wfntype==4) orbcoeff(ipos6D:ipos6D+n6D-1,1:2*nbasis5D)=orbcoeff5D(ipos5D:ipos5D+n5D-1,:)
        else if (n5D==5) then
            !5D->6D
            if (wfntype==3) orbcoeff(ipos6D:ipos6D+5,1:nbasis5D)=matmul(conv5d6d,orbcoeff5D(ipos5D:ipos5D+4,:))
            if (wfntype==4) orbcoeff(ipos6D:ipos6D+5,1:2*nbasis5D)=matmul(conv5d6d,orbcoeff5D(ipos5D:ipos5D+4,:))
        else if (n5D==7) then
            !7F->10F
            if (wfntype==3) orbcoeff(ipos6D:ipos6D+9,1:nbasis5D)=matmul(conv7f10f,orbcoeff5D(ipos5D:ipos5D+6,:))
            if (wfntype==4) orbcoeff(ipos6D:ipos6D+9,1:2*nbasis5D)=matmul(conv7f10f,orbcoeff5D(ipos5D:ipos5D+6,:))
        else if (n5D==9) then
            !9G->15G
            if (wfntype==3) orbcoeff(ipos6D:ipos6D+14,1:nbasis5D)=matmul(conv9g15g,orbcoeff5D(ipos5D:ipos5D+8,:))
            if (wfntype==4) orbcoeff(ipos6D:ipos6D+14,1:2*nbasis5D)=matmul(conv9g15g,orbcoeff5D(ipos5D:ipos5D+8,:))
        end if
        ipos5D=ipos5D+n5D
        ipos6D=ipos6D+n6D
    end do
end if

nprims=0
do i=1,nshell
    nprims=nprims+shellcon(i)*shellnumbas(i)
end do
allocate(b(nprims),co(nmo,nprims))

iGTF=1 !current GTF index
ibasis=1
do i=1,nshell !cycle each shell
    b(iGTF:iGTF+shellcon(i)*shellnumbas(i)-1)%center=shell2atom(i)
    do j=1,shellnumbas(i) !cycle each basis(orbital) in each shell
        b(iGTF:iGTF+shellcon(i)-1)%functype=bastype2func(bastype31(ibasis))
        do k=1,shellcon(i) !cycle each GTF in each basis in each shell
            iprmshlpos=shell2prmshl(i)+k-1
            b(iGTF)%exp=prmshlexp(iprmshlpos)
            if (bastype31(ibasis)==1) then !s
                contract=cs(iprmshlpos)
            else if (bastype31(ibasis)<=200) then !p
                contract=cp(iprmshlpos)
            else if (bastype31(ibasis)<=300) then !d
            !Contract coefficient in .31 contains normalization coefficient, however for d type, the
            !normalization coefficient is for XX,YY,ZZ, for XY,XZ,YZ, we need refresh normalization coefficient
                contract=cd(iprmshlpos)
                if (bastype31(ibasis)==202.or.bastype31(ibasis)==203.or.bastype31(ibasis)==205) then
                    valnorm31=normgau(5,prmshlexp(iprmshlpos)) !Normalization coefficient for XX,YY,ZZ are identical
                    valnormnew=normgau(8,prmshlexp(iprmshlpos)) !Normalization coefficient for XY,XZ,YZ are identical
                    contract=contract/valnorm31*valnormnew
                end if
            else if (bastype31(ibasis)<=400) then !f
                contract=cf(iprmshlpos)
                !For f shell, in .31 normalization coefficient is for XXX,YYY,ZZZ, now refresh
                if (bastype31(ibasis)/=301.and.bastype31(ibasis)/=307.and.bastype31(ibasis)/=310) then !not XXX,YYY,ZZZ
                    valnorm31=normgau(11,prmshlexp(iprmshlpos))
                    if (bastype31(ibasis)==302.or.bastype31(ibasis)==303.or.bastype31(ibasis)==304&
                    .or.bastype31(ibasis)==306.or.bastype31(ibasis)==308.or.bastype31(ibasis)==309) then
                        valnormnew=normgau(14,prmshlexp(iprmshlpos)) !XXY,XXZ,XYY,XZZ,YYZ,YZZ are identical
                    else if (bastype31(ibasis)==305) then 
                        valnormnew=normgau(20,prmshlexp(iprmshlpos)) !XYZ
                    end if
                    contract=contract/valnorm31*valnormnew
                end if
            else if (bastype31(ibasis)<=500) then !g
                contract=cg(iprmshlpos)
                !For g shell, in .31 normalization coefficient is for XXXX,YYYY,ZZZZ, now refresh
                !Note: I haven't verified that, since .37 outputted by NBO3.1 in Gaussian didn't contains g information
                nt=bastype31(ibasis) !now type
                if (nt/=401.and.nt/=411.and.nt/=415) then !not XXXX,YYYY,ZZZZ (4,0)
                    valnorm31=normgau(21,prmshlexp(iprmshlpos))
                    if (nt==402.or.nt==403.or.nt==407.or.nt==410.or.nt==412.or.nt==414) then
                        valnormnew=normgau(22,prmshlexp(iprmshlpos)) !XXXY,XXXZ,XYYY,XZZZ,YYYZ,YZZZ (3,1)
                    else if (nt==404.or.nt==406.or.nt==413) then
                        valnormnew=normgau(23,prmshlexp(iprmshlpos)) !XXYY,XXZZ,YYZZ (2,2)
                    else if (nt==405.or.nt==408.or.nt==409) then
                        valnormnew=normgau(27,prmshlexp(iprmshlpos)) !XYZZ,XYYZ,XXYZ (2,1,1)
                    end if
                    contract=contract/valnorm31*valnormnew
                end if
            end if
            CO(:,iGTF)=orbcoeff(ibasis,:)*contract
            iGTF=iGTF+1
        end do
        ibasis=ibasis+1
    end do
end do
if (isphergau==1) nbasis=nbasis5D

if (infomode==0) then
    if (name2(itmplen-1:itmplen)=="37".or.name2(itmplen-1:itmplen)=="39") write(*,"(' Total/Alpha/Beta electrons:',3f12.4)") nelec,naelec,nbelec
    write(*,"(' Atoms:',i6,',  Basis functions:',i6,',  Orbitals:',i6,',  GTFs:',i6)") ncenter,nbasis,nmo,nprims
    if (wfntype==3) write(*,*) " This is close-shell system"
    if (wfntype==4) then
        write(*,*) " This is open-shell system"
        write(*,"(' Orbitals from 1 to',i6,' are alpha, from',i6,' to',i6,' are beta')") nbasis,nbasis+1,nmo
    end if
    write(*,*)
end if
end subroutine




!!-----------------------------------------------------------------
!!---------------- Read Gaussian cube file and store in cubmat
!infomode=0 means output cube details, =1 not
!ionlygrid=1 means only read grid data, but do not perturb any other variables, =0 means do all
subroutine readcube(cubname,infomode,ionlygrid)
use defvar
implicit real*8 (a-h,o-z)
character(len=*) cubname
character*79 titleline1,titleline2
integer infomode,ionlygrid
integer,allocatable :: mo_serial(:)
real*8,allocatable :: temp_readdata(:)
type(content) maxv,minv
if (ionlygrid==0) ifiletype=7
open(10,file=cubname,access="sequential",status="old")
read(10,"(a)") titleline1
read(10,"(a)") titleline2
read(10,*) ncentertmp,orgx,orgy,orgz
read(10,*) nx,v1x,v1y,v1z
read(10,*) ny,v2x,v2y,v2z
read(10,*) nz,v3x,v3y,v3z
if (ionlygrid==0) ncenter=ncentertmp
dx=v1x !If this is cubic grid. dx,dy,dz are global array in defvar
dy=v2y
dz=v3z
endx=orgx+v1x*(nx-1) !endx,y,z are global array in defvar
endy=orgy+v2y*(ny-1)
endz=orgz+v3z*(nz-1)

mo_number=0
if (ncenter<0) then
    mo_number=1 !This cube file contains at least one MO data
    ncenter=abs(ncenter)
end if
if (ionlygrid==0) then
    if (allocated(a)) deallocate(a)
    allocate(a(ncenter))
end if
if (allocated(cubmat)) deallocate(cubmat)
allocate(cubmat(nx,ny,nz))

if (infomode==0) then
    write(*,*) "Title line of this file:"
    write(*,"(a)") trim(titleline1)
    write(*,"(a)") trim(titleline2)
    write(*,*)
    write(*,"(' Total number of atoms:',i8)") ncenter
    write(*,"(' Translation vector:        X           Y           Z     (Bohr)')")
    write(*,"(a20,3f12.6)") "Vector 1: ",v1x,v1y,v1z
    write(*,"(a20,3f12.6)") "Vector 2: ",v2x,v2y,v2z
    write(*,"(a20,3f12.6)") "Vector 3: ",v3x,v3y,v3z
    write(*,"(' The range of x is from ',f12.6,' to ',f12.6,' Bohr,' i5,' points')") ,orgx,orgx+(nx-1)*v1x,nx
    write(*,"(' The range of y is from ',f12.6,' to ',f12.6,' Bohr,',i5,' points')") ,orgy,orgy+(ny-1)*v2y,ny
    write(*,"(' The range of z is from ',f12.6,' to ',f12.6,' Bohr,',i5,' points')") ,orgz,orgz+(nz-1)*v3z,nz
    write(*,"(' Total number of grid points:',i10)") nx*ny*nz
    write(*,"(' This grid data will take up at least',i6,' MB memory')") nx*ny*nz*8/1024/1024
end if

if (ionlygrid==0) then
    do i=1,ncenter
        read(10,*) a(i)%index,a(i)%charge,a(i)%x,a(i)%y,a(i)%z !%value is its charge, if ECP was used, it not equal to atomindex
    end do
    a%name=ind2name(a%index)
else if (ionlygrid==1) then
    do i=1,ncentertmp
        read(10,*) !Do not read atomic information, simply skip
    end do
end if
write(*,*)

if (mo_number==1) then
    read(10,"(i5)",advance="no") mo_number !Get actual number of MO
    if (mo_number>1) then
        allocate(mo_serial(mo_number))
        allocate(temp_readdata(nz*mo_number))
        read(10,*) mo_serial
        write(*,"(' There are ',i6,' MOs in this cube file, the serial numbers are: ')") mo_number
        do i=1,mo_number
            write(*,"(' Number ',i6,', corresponds to MO',i6)") i,mo_serial(i)
        end do
        write(*,*) "Which MO do you want to load? Input the serial number"
        do while(.true.)
            read(*,*) mo_select
            if (mo_select>0.and.mo_select<=mo_number) exit
            write(*,*) "Invalid input, input again"
        end do
    else
        read(10,*) !Only one MO, pass the MO serial line
    end if
end if

write(*,*) "Loading grid data, please wait..."
!Load data
ii=0
do i=1,nx   !a(x,y,z)
    do j=1,ny
        if (mo_number==0.or.mo_number==1) then
            read(10,*) cubmat(i,j,:)
        else !Load the specified MO from vast of MOs
            read(10,*) temp_readdata
            cubmat(i,j,:)=temp_readdata(mo_select:size(temp_readdata):mo_number)
        end if
    end do
    progress=dfloat(i)/nx*100
    if (progress>ii) then
        ii=ii+10
        write(*,"(f6.1,'%')") progress
    end if
end do

write(*,"(f6.1,'%')") 100D0
write(*,*) "Done!"
write(*,*)
close(10)

maxv%value=cubmat(1,1,1)
maxv%x=orgx
maxv%y=orgy
maxv%z=orgz
minv%value=cubmat(1,1,1)
minv%x=orgx
minv%y=orgy
minv%z=orgz
sumuppos=0.0D0
sumupneg=0.0D0
do k=1,nz
    do j=1,ny
        do i=1,nx
            if (cubmat(i,j,k)>0) sumuppos=sumuppos+cubmat(i,j,k)
            if (cubmat(i,j,k)<0) sumupneg=sumupneg+cubmat(i,j,k)
            if (cubmat(i,j,k)>maxv%value) then
                maxv%value=cubmat(i,j,k)
                maxv%x=orgx+(i-1)*dx
                maxv%y=orgy+(j-1)*dy
                maxv%z=orgz+(k-1)*dz
            end if
            if (cubmat(i,j,k)<minv%value) then
                minv%value=cubmat(i,j,k)
                minv%x=orgx+(i-1)*dx
                minv%y=orgy+(j-1)*dy
                minv%z=orgz+(k-1)*dz
            end if
        end do
    end do
end do

if (infomode==0) then
    fminivol=v1x*v2y*v3z
    write(*,"(' The minimum value:',D16.8,' at',3f12.6,' Bohr')") minv%value,minv%x,minv%y,minv%z
    write(*,"(' The maximum value:',D16.8,' at',3f12.6,' Bohr')") maxv%value,maxv%x,maxv%y,maxv%z
    write(*,"(' Differential element:',f15.10,' Bohr^3')") fminivol
    write(*,"(' Summing up positive value in grid file:  ',f30.10)") sumuppos
    write(*,"(' After multiplied by differential element:',f30.10)") sumuppos*fminivol
    write(*,"(' Summing up negative value in grid file:  ',f30.10)") sumupneg
    write(*,"(' After multiplied by differential element:',f30.10)") sumupneg*fminivol
    write(*,"(' Summing up all value in grid file:       ',f30.10)") sumuppos+sumupneg
    write(*,"(' After multiplied by differential element:',f30.10)") (sumuppos+sumupneg)*fminivol
end if
end subroutine

!!-----------------------------------------------------------------
!!---- Read Gaussian cube file and save to cubmattmp, this is a simple version of readcube, can only be invoked after cubmat has been loaded
!don't read atomic information, don't modify loaded grid infomation such as nx/y/z,orgx/y/z...
!and don't output statistic information, don't specify coordinate for grid points...
!If inconsis==1, that means the grid setting of this cube file is inconsistent with that of cubmat
subroutine readcubetmp(cubname,inconsis)
use defvar
implicit real*8 (a-h,o-z)
character(len=*) cubname
integer inconsis
integer,allocatable :: mo_serial(:)
real*8,allocatable :: temp_readdata(:)
open(10,file=cubname,access="sequential",status="old")
read(10,*)
read(10,*)
read(10,*) ncentertmp
read(10,*) nxtmp,dxtmp,tmpval,tmpval
read(10,*) nytmp,tmpval,dytmp,tmpval
read(10,*) nztmp,tmpval,tmpval,dztmp
inconsis=0
if (nxtmp/=nx.or.nytmp/=ny.or.nztmp/=nz.or.abs(dxtmp-dx)>0.005D0.or.abs(dytmp-dy)>0.005D0.or.abs(dztmp-dz)>0.005D0) then
    write(*,"(a)") " Error: The grid setting of this cube file is inconsistent with that of the grid data stored in memory!"
    inconsis=1
    return
end if

mo_number=0
if (ncentertmp<0) then
    mo_number=1 !This cube file contains at least one MO data
    ncentertmp=abs(ncentertmp)
end if
if (allocated(cubmattmp)) deallocate(cubmattmp)
allocate(cubmattmp(nx,ny,nz)) !nx,ny,nz is identical to cubmat that already loaded into memory

do i=1,ncentertmp
    read(10,*) !Skip a(i)%index,a(i)%charge,a(i)%x,a(i)%y,a(i)%z !%value is its charge, if ECP was used, it not equal to atomindex
end do

if (mo_number==1) then
    read(10,"(i5)",advance="no") mo_number !Get actual number of MO
    if (mo_number>1) then
        allocate(mo_serial(mo_number))
        allocate(temp_readdata(nz*mo_number))
        read(10,*) mo_serial
        write(*,"(' There are ',i6,' MOs in this grid file, the serial number are: ')") mo_number
        do i=1,mo_number
            write(*,"(' Number ',i6,' : MO= ',i5)") i,mo_serial(i)
        end do
        write(*,*) "Which MO do you want to load? Input the serial number"
        do while(.true.)
            read(*,*) mo_select
            if (mo_select>0.and.mo_select<=mo_number) exit
            write(*,*) "Invalid input, input again"
        end do
    else
        read(10,*) !Only one MO, pass the MO serial line
    end if
end if

write(*,*)
write(*,*) "Loading grid data, please wait..."
!Load data
ii=0
do i=1,nx   !a(x,y,z)
    do j=1,ny
        if (mo_number==0.or.mo_number==1) then
            read(10,*) cubmattmp(i,j,:)
        else !Load the specified MO from vast of MOs
            read(10,*) temp_readdata
            cubmattmp(i,j,:)=temp_readdata(mo_select:size(temp_readdata):mo_number)
        end if
    end do
    progress=dfloat(i)/nx*100
    if (progress>ii) then
        ii=ii+10
        write(*,"(f6.1,'%')") progress
    end if
end do
write(*,*) "Done!"
close(10)
end subroutine



!!-----------------------------------------------------------------
!!----------- Load Dmol3 .grd file, see format description in Material Studio help file
!infomode=0 means output grd details, =1 not
!ionlygrid=1 means only read grid data, but do not perturb any other variables, =0 means do all
subroutine readgrd(grdname,infomode,ionlygrid)
use defvar
implicit real*8 (a-h,o-z)
character(len=*) grdname
character*79 titleline
integer infomode,ionlygrid
type(content) maxv,minv
if (ionlygrid==0) ifiletype=8
open(10,file=grdname,access="sequential",status="old")
read(10,"(a)") titleline
read(10,*)
read(10,*) flenx,fleny,flenz,anga,angb,angc !Notice that the length unit is in Angstrom in .grd file
read(10,*) nx,ny,nz !Here nx,ny,nz are total space (between neighbour grid point) in each direction
read(10,*) ifast,ixback,ixforw,iyback,iyforw,izback,izforw
if (anga/=90.or.angb/=90.or.angc/=90) then
    write(*,*) "Error: Only cubic cell is supported in Multiwfn! Press ENTER to exit"
    read(*,*)
    stop
else if (ifast/=1) then !ifast=1 means x varies fastest
    write(*,*) "Error: The first integer in the fifth line must be 1! Press ENTER to exit"
    read(*,*)
    stop
end if
dx=flenx/nx/b2a
dy=fleny/ny/b2a
dz=flenz/nz/b2a
nx=nx+1 !Convert the number of spacings to the number of points
ny=ny+1
nz=nz+1
allocate(cubmat(nx,ny,nz))
orgx=-dx*abs(ixback)
orgy=-dy*abs(iyback)
orgz=-dz*abs(izback)
endx=orgx+dx*(nx-1) !endx,y,z are global array in defvar
endy=orgy+dy*(ny-1)
endz=orgz+dz*(nz-1)

if (infomode==0) then
    write(*,"(' Title line of this file: ',a)") trim(titleline)
    write(*,*)
    write(*,"(' Translation vectors in X/Y/Z (Bohr):',3f12.6)") dx,dy,dz
    write(*,"(' The range of x is from ',f12.6,' to ',f12.6,' Bohr,' i5,' points')") ,orgx,endx,nx
    write(*,"(' The range of y is from ',f12.6,' to ',f12.6,' Bohr,',i5,' points')") ,orgy,endy,ny
    write(*,"(' The range of z is from ',f12.6,' to ',f12.6,' Bohr,',i5,' points')") ,orgz,endz,nz
    write(*,"(' Total number of grid points:',i10)") nx*ny*nz
    write(*,"(' This grid data will take up at least',i6,' MB memory')") nx*ny*nz*8*4/1024/1024
end if
write(*,*)
write(*,*) "Loading grid data, please wait..."
ii=0
do k=1,nz   !a(x,y,z)
    do j=1,ny
        do i=1,nx
            read(10,*) cubmat(i,j,k)
        end do
    end do
    if (dfloat(k)/nz*100>ii) then
        ii=ii+10
        write(*,"(f6.1,'%')") dfloat(k)/nz*100
    end if
end do
write(*,"(f6.1,'%')") 100D0
write(*,*) "Done!"
write(*,*)
close(10)

!Perform statistic
maxv%value=cubmat(1,1,1)
maxv%x=orgx
maxv%y=orgy
maxv%z=orgz
minv%value=cubmat(1,1,1)
minv%x=orgx
minv%y=orgy
minv%z=orgz
sumuppos=0.0D0
sumupneg=0.0D0
do k=1,nz
    do j=1,ny
        do i=1,nx
            if (cubmat(i,j,k)>0) sumuppos=sumuppos+cubmat(i,j,k)
            if (cubmat(i,j,k)<0) sumupneg=sumupneg+cubmat(i,j,k)
            if (cubmat(i,j,k)>maxv%value) then
                maxv%value=cubmat(i,j,k)
                maxv%x=orgx+(i-1)*dx
                maxv%y=orgy+(j-1)*dy
                maxv%z=orgz+(k-1)*dz
            end if
            if (cubmat(i,j,k)<minv%value) then
                minv%value=cubmat(i,j,k)
                minv%x=orgx+(i-1)*dx
                minv%y=orgy+(j-1)*dy
                minv%z=orgz+(k-1)*dz
            end if
        end do
    end do
end do

if (infomode==0) then
    fminivol=dx*dy*dz
    write(*,"(' The minimum value:',D16.8,' at',3f12.6,' Bohr')") minv%value,minv%x,minv%y,minv%z
    write(*,"(' The maximum value:',D16.8,' at',3f12.6,' Bohr')") maxv%value,maxv%x,maxv%y,maxv%z
    write(*,"(' Differential element:',f15.10,' Bohr^3')") fminivol
    write(*,"(' Summing up positive value in grid file:  ',f30.10)") sumuppos
    write(*,"(' After multiplied by differential element:',f30.10)") sumuppos*fminivol
    write(*,"(' Summing up negative value in grid file:  ',f30.10)") sumupneg
    write(*,"(' After multiplied by differential element:',f30.10)") sumupneg*fminivol
    write(*,"(' Summing up all value in grid file:       ',f30.10)") sumuppos+sumupneg
    write(*,"(' After multiplied by differential element:',f30.10)") (sumuppos+sumupneg)*fminivol
end if
end subroutine

!!-----------------------------------------------------------------
!!----------- Load Dmol3 .grd file and save to cubmattmp
!If inconsis==1, that means the grid setting of this grd file is inconsistent with cubmat stored in memory
subroutine readgrdtmp(grdname,inconsis)
use defvar
implicit real*8 (a-h,o-z)
character(len=*) grdname
open(10,file=grdname,access="sequential",status="old")
read(10,*)
read(10,*)
read(10,*)
read(10,*) nxtmp,nytmp,nztmp
read(10,*)
if (nxtmp+1/=nx.or.nytmp+1/=ny.or.nztmp+1/=nz) then
    write(*,"(a)") " Error: The grid setting of this grd file is inconsistent with that of the grid data stored in memory!"
    inconsis=1
    return
end if

if (allocated(cubmattmp)) deallocate(cubmattmp)
allocate(cubmattmp(nx,ny,nz))
write(*,*)
write(*,*) "Loading data, please wait..."
do k=1,nz   !a(x,y,z)
    do j=1,ny
        do i=1,nx
            read(10,*) cubmattmp(i,j,k)
        end do
    end do
end do
write(*,*) "Grid data loading completed!"
close(10)
end subroutine



!!-----------------------------------------------------------------
!!-------------------- Read .wfn file, infomode=0 means output wfn property, =1 not
subroutine readwfn(name,infomode)
use defvar
use util
implicit real*8 (a-h,o-z)
CHARACTER(LEN=*) name
character*80 wfntitle,lastline,c80tmp*80
real*8,allocatable :: tmpCO(:,:),tmpMOocc(:),tmpMOene(:)
integer,allocatable :: tmpMOtype(:)
integer i,j,infomode
!Original .wfn format doesn't support g, however the .wfn outputted by Multiwfn, Molden2AIM and G09 since B.01 formally supports g
!Below is the g sequence used in Molden2AIM, .wfx, .molden and the .wfn outputted by Multiwfn and G09 since B.01
! 21 XXXX 22 YYYY 23 ZZZZ 24 XXXY 25 XXXZ
! 26 XYYY 27 YYYZ 28 XZZZ 29 YZZZ 30 XXYY
! 31 XXZZ 32 YYZZ 33 XXYZ 34 XYYZ 35 XYZZ
!Below is the g sequence internally used in Multiwfn, identical to .fch
! 21 ZZZZ 22 YZZZ 23 YYZZ 24 YYYZ 25 YYYY
! 26 XZZZ 27 XYZZ 28 XYYZ 29 XYYY 30 XXZZ
! 31 XXYZ 32 XXYY 33 XXXZ 34 XXXY 35 XXXX
! convGseq is used to convert g used in .wfn to the internal sequence of Multiwfn 
! I assume that h sequence is identical to .wfx and Multiwfn
integer :: convGseq(35)=(/ (0,i=1,20), 35,25,21,34,33, 29,24,26,22,32, 30,23,31,28,27 /)
ifiletype=2
imodwfn=0
open(10,file=name,access="sequential",status="old")
read(10,"(a)") wfntitle
read(10,"(a8,i15,13x,i7,11x,i9)") c80tmp,nmo,nprims,ncenter
ibasmode=1
if (index(c80tmp,"SLATER")/=0) ibasmode=2
if (ibasmode==2) then
    write(*,"(a)") " Error: Multiwfn does not support the wfn file recording Slater type orbitals! Press ENTER to exit"
    read(*,*)
    stop
end if
if (allocated(a)) deallocate(a)
allocate(a(ncenter))
allocate(b(nprims))
allocate(co(nmo,nprims))
allocate(MOocc(nmo))
allocate(MOene(nmo))
allocate(MOtype(nmo))
do i=1,ncenter
    read(10,"(a24,3f12.8,10x,f5.1)") c80tmp,a(i)%x,a(i)%y,a(i)%z,a(i)%charge
    read(c80tmp,*) a(i)%name
    call lc2uc(a(i)%name(1:1))
    call uc2lc(a(i)%name(2:2))
    do j=0,nelesupp
        if (a(i)%name==ind2name(j)) then
            a(i)%index=j
            exit
        end if
        if (j==nelesupp) then
            write(*,"(3a,i5)") " Warning: Found unknown element ",a(i)%name," with index of",i
            write(*,*) "The atom will be recognized as hydrogen internally without nuclear charge"
            write(*,*) "Press ENTER to continue"
            read (*,*)
            a(i)%index=1
        end if
    end do
end do
read(10,"(20x,20i3)") (b(i)%center,i=1,nprims)
read(10,"(20x,20i3)") (b(i)%functype,i=1,nprims)
do i=1,nprims    
    if (b(i)%functype>=21.and.b(i)%functype<=35) b(i)%functype=convGseq(b(i)%functype)
end do

read(10,"(a)") c80tmp
backspace(10)
if (c80tmp(11:11)==' ') then !This is a standard .wfn file
    read(10,"(10x,5D14.7)") (b(i)%exp,i=1,nprims)
else
    write(*,"(a)") " Note: This .wfn file is probably generated by ORCA, special treatment is employed since it is non-standard"
    if (c80tmp(25:25)==' ') read(10,"(9x,5D15.7)") (b(i)%exp,i=1,nprims) !Windows
    if (c80tmp(24:24)==' ') read(10,"(9x,5D14.7)") (b(i)%exp,i=1,nprims) !Linux
end if

!From Gaussian09 B.01, if ECP is used, additional CENTER, TYPE, EXPONENTS field present to represent EDF(electron density functions) likewise .wfx file
!However such wfn file is foolish, the number of GTFs used to represent EDF is not explicitly recorded, so we must use trick to guess the number
!The coefficient of EDF is not written, how we use these EDF information? Obviously impossible! So we just skip these mad fields
read(10,"(a)") c80tmp
if (c80tmp(1:6)=="CENTRE") then
    call loclabel(10,"MO",ifound,0)
!     ncenline=1
!     do while(.true.) !Move to the first line of TYPE field, and count we passed how many rows (the number of lines of CENTER field)
!         read(10,"(a)") c80tmp
!         if (c80tmp(1:6)/="CENTRE") exit
!         ncenline=ncenline+1
!     end do
!     backspace(10)
!     backspace(10) !Return to the last line of CENTER field
!     read(10,"(a)") c80tmp
!     nEDFprims=(len_trim(c80tmp)-20)/3+(ncenline-1)*20
!     do itmp=1,ncenline !Move to the first line of CENTER field, make ready to load EDF information
!         backspace(10)
!     end do
!     allocate(b_EDF(nEDFprims),CO_EDF(nEDFprims))
!     read(10,"(20x,20i3)") (b_EDF(i)%center,i=1,nEDFprims)
!     read(10,"(20x,20i3)") (b_EDF(i)%functype,i=1,nEDFprims)
!     read(10,"(10x,5D14.7)") (b_EDF(i)%exp,i=1,nEDFprims)
else
    backspace(10)
end if

!Read orbitals
do i=1,nmo
!     read(10,"(35X,f12.7,15x,f12.6)") MOocc(i),MOene(i) !The life is not so easy, programs such as ADF doesn't follow the standard rule
    read(10,"(a)") c80tmp
    do ichar=1,80
        if (c80tmp(ichar:ichar)=='=') then
            read(c80tmp(ichar+1:),*) MOocc(i)
            exit
        end if
    end do
    do ichar=80,1,-1
        if (c80tmp(ichar:ichar)=='=') then
            read(c80tmp(ichar+1:),*) MOene(i)
            exit
        end if
    end do
!     read(10,"(5D16.8)") (co(i,j),j=1,nprims) ! Note: row/column of CO denote MO/basis function respectively, in contrary to convention
    read(10,*) co(i,:) !Read data in this way has better compatability
end do
read(10,*)
!Use free format to read in energy and virial to ensure compatibility
read(10,"(a)") lastline
iequalsign1=0
iequalsign2=0
do i=1,80
    if (lastline(i:i)=='=') then
        iequalsign1=i
        exit
    end if
end do
do i=80,1,-1
    if (lastline(i:i)=='=') then
        iequalsign2=i
        exit
    end if
end do
totenergy=0
virialratio=2
if (iequalsign1/=0) read(lastline(iequalsign1+1:),*) totenergy
if (iequalsign1==0) write(*,*) "Warning: Unable to find system energy in this file!"
if (iequalsign2/=0) read(lastline(iequalsign2+1:),*) virialratio
if (iequalsign2==0) write(*,*) "Warning: Unable to find Virial in this file!"
call loclabel(10,"$MOSPIN $END",ifoundmospin,0) !Also read spin-type from $MOSPIN $END field outputted by Molden2AIM since v2.0.5, if detected
if (ifoundmospin==1) then !Have defined spin-type explicitly, don't reset spin-type by guessing later
    if (infomode==0) write(*,*) "Note: Found $MOSPIN $END field, orbital spin-types are loaded"
    read(10,*)
    read(10,*) MOtype
    where (MOtype==3) MOtype=0
end if

close (10)

! Determine type of wavefunction
if (sum(MOocc)==2*nmo.and.all(int(MOocc)==MOocc)) then
    wfntype=0 !This is restricted wavefunction
    MOtype=0
else if (sum(MOocc)==nmo.and.all(int(MOocc)==MOocc)) then
    wfntype=1 !This is unrestricted wavefunction
    if (ifoundmospin==0) then
        MOtype=1 !Set all MO is alpha
        do i=2,nmo !if nmo=1, i will be set to 2, and no errors will appear
            if (MOene(i)<MOene(i-1)) exit
        end do
        MOtype(i:nmo)=2 !beta
    end if
else if (any(MOocc/=int(MOocc))) then
    if (nint(maxval(MOocc))==2) then !maximum occupation close to 2, so considered as restricted Post-HF wavefunction
        wfntype=3
        MOtype=0
    else
        wfntype=4 !This is unrestricted Post-HF wavefunction
        if (ifoundmospin==0) then
            MOtype=0
            do i=2,nmo
                if (MOocc(i)>MOocc(i-1)) then
                    MOtype(1:i-1)=1
                    MOtype(i:nmo)=2
                    exit
                end if
            end do
        end if
    end if
else
    wfntype=2 !This is RO wavefunction
    MOtype=0
    do i=1,nmo
        if (MOocc(i)==1) MOtype(i)=1 !alpha
    end do
end if
!Count electrons
call updatenelec

!Sort orbitals so that the orbitals with same spin-type are contiguous, because the wfn file outputted by molden2aim is not always in line with this convention
if (ifoundmospin==1.and.(wfntype==1.or.wfntype==4)) then
    allocate(tmpCO(nmo,nprims),tmpMOtype(nmo),tmpMOocc(nmo),tmpMOene(nmo))
    ipos=0
    do itime=1,2
        do imo=1,nmo
            if ((itime==1.and.MOtype(imo)==1).or.(itime==2.and.MOtype(imo)==2)) then !Move alpha orbitals to tmp arrays
                ipos=ipos+1
                tmpCO(ipos,:)=CO(imo,:)
                tmpMOocc(ipos)=MOocc(imo)
                tmpMOene(ipos)=MOene(imo)
                tmpMOtype(ipos)=MOtype(imo)
            end if
        end do
    end do
    CO=tmpCO
    MOocc=tmpMOocc
    MOene=tmpMOene
    MOtype=tmpMOtype
    write(*,*) "Note: Sequence of orbitals has been sorted according to spin-type"
    deallocate(tmpCO,tmpMOocc,tmpMOene,tmpMOtype)
end if

!Summary
if (infomode==0) then
    write(*,*)
    write(*,"(' System energy:',f21.12,' Hartree,   Virial ratio:',f12.8)") totenergy,virialratio
    write(*,"(' Total/Alpha/Beta electrons:',3f12.4)") nelec,naelec,nbelec
    write(*,"(' Net charge:',f12.5,'    Expected multiplicity:',i5)") sum(a(:)%charge)-nelec,nint(naelec-nbelec)+1
    write(*,"(' The number of orbitals:',i6,',  Atoms:',i7,',  GTFs:',i7)") nmo,ncenter,nprims
    if (wfntype==0) write(*,"(' This is restricted close-shell single-determinant wavefunction')")
    if (wfntype==1) write(*,"(' This is unrestricted single-determinant wavefunction')")
    if (wfntype==2) write(*,"(' This is restricted open-shell wavefunction')")
    if (wfntype==3) write(*,"(' This is close-shell post-HF wavefunction')")
    if (wfntype==4) write(*,"(' This is open-shell post-HF wavefunction')")
    if (wfntype==1.or.wfntype==4) then
        do i=1,nmo
            if (MOtype(i)==2) exit
        end do
        write(*,"(' Orbitals from 1 to',i6,' are alpha type, from',i6,' to',i6,' are beta type')") i-1,i,nmo
    end if
    write(*,"(' Title line of this file: ',a)") trim(wfntitle)
end if
end subroutine



!!-----------------------------------------------------------------
!!------------------------- Read .wfx, mode=0 means output wfn property,=1 not
!For more information, see http://aim.tkgristmill.com/wfxformat.html
subroutine readwfx(name,infomode)
use defvar
use util
implicit real*8 (a-h,o-z)
CHARACTER(LEN=*) name
character spintype*20
integer infomode
!Below is the g sequence used in Molden2AIM, .wfx, .molden and the .wfn outputted by Multiwfn and G09 since B.01
! 21 XXXX 22 YYYY 23 ZZZZ 24 XXXY 25 XXXZ
! 26 XYYY 27 YYYZ 28 XZZZ 29 YZZZ 30 XXYY
! 31 XXZZ 32 YYZZ 33 XXYZ 34 XYYZ 35 XYZZ
!Below is the g sequence internally used in Multiwfn, identical to .fch
! 21 ZZZZ 22 YZZZ 23 YYZZ 24 YYYZ 25 YYYY
! 26 XZZZ 27 XYZZ 28 XYYZ 29 XYYY 30 XXZZ
! 31 XXYZ 32 XXYY 33 XXXZ 34 XXXY 35 XXXX
! convGseq is used to convert g used in .wfx to the internal sequence of Multiwfn
! PS: spdfh sequence in .wfx is identical to Multiwfn
integer :: convGseq(35)=(/ (0,i=1,20), 35,25,21,34,33, 29,24,26,22,32, 30,23,31,28,27 /)
ifiletype=3
imodwfn=0
open(10,file=name,access="sequential",status="old")
call loclabel(10,"<Number of Nuclei>")
read(10,*)
read(10,*) ncenter
if (allocated(a)) deallocate(a)
allocate(a(ncenter))
call loclabel(10,"<Number of Primitives>")
read(10,*)
read(10,*) nprims
allocate(b(nprims))
call loclabel(10,"<Number of Occupied Molecular Orbitals>")
read(10,*)
read(10,*) nmo
allocate(MOocc(nmo),MOene(nmo),MOtype(nmo),co(nmo,nprims))
call loclabel(10,"<Nuclear Names>")
read(10,*)
read(10,*) a%name
do i=1,ncenter !Multiwfn does not allow number included in atom name
    if (iachar(a(i)%name(2:2))<=57.and.iachar(a(i)%name(2:2))>=48) a(i)%name(2:2)=' '
end do
call loclabel(10,"<Atomic Numbers>")
read(10,*)
read(10,*) a%index
call loclabel(10,"<Nuclear Charges>")
read(10,*)
read(10,*) a%charge
call loclabel(10,"<Nuclear Cartesian Coordinates>")
read(10,*)
do i=1,ncenter
    read(10,*) a(i)%x,a(i)%y,a(i)%z
end do
call loclabel(10,"<Number of Electrons>")
read(10,*)
read(10,*) nelec
call loclabel(10,"<Number of Alpha Electrons>")
read(10,*)
read(10,*) naelec
call loclabel(10,"<Number of Beta Electrons>")
read(10,*)
read(10,*) nbelec
call loclabel(10,"<Primitive Centers>")
read(10,*)
read(10,*) b%center
call loclabel(10,"<Primitive Types>")
read(10,*)
read(10,*) b%functype
!The g sequence in .wfx is not identical to Multiwfn, convert them here
do i=1,nprims    
    if (b(i)%functype>=21.and.b(i)%functype<=35) b(i)%functype=convGseq(b(i)%functype)
end do
call loclabel(10,"<Primitive Exponents>")
read(10,*)
read(10,*) b%exp
call loclabel(10,"<Molecular Orbital Occupation Numbers>")
read(10,*)
read(10,*) MOocc
call loclabel(10,"<Molecular Orbital Energies>")
read(10,*)
read(10,*) MOene
call loclabel(10,"<Molecular Orbital Spin Types>")
read(10,*)
do i=1,nmo
    read(10,"(a20)") spintype
    if (adjustl(spintype)=="Alpha and Beta") MOtype(i)=0 !adjustl is needed, because the wfx outputted by ORCA is non-standard
    if (adjustl(spintype)=="Alpha") MOtype(i)=1
    if (adjustl(spintype)=="Beta") MOtype(i)=2
end do
call loclabel(10,"<Molecular Orbital Primitive Coefficients>")
read(10,*)
do i=1,nmo
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*) CO(i,:)
end do
if (any(b%functype>56)) then !Angular moment of GTF should be no higher than h
    write(*,"(' Warning: Angular moment of one or more GTFs exceeds h, Multiwfn is unable to deal with this case! Its/their contributions will be discarded')")
    read (*,*)
    do iGTF=1,nprims
        if (b(iGTF)%functype>56) then
            b(iGTF)%functype=1 !Assume it is S type
            CO(:,iGTF)=0D0
        end if
    end do
end if
call loclabel(10,"<Energy = T + Vne + Vee + Vnn>")
read(10,*)
read(10,*) totenergy
call loclabel(10,"<Virial Ratio (-V/T)>")
read(10,*)
read(10,*) virialratio
!------process EDF information
call loclabel(10,"<Number of EDF Primitives>",ifound)
if (ifound==1.and.readEDF==1) then
    read(10,*)
    read(10,*) nEDFprims
    allocate(b_EDF(nEDFprims),CO_EDF(nEDFprims))
    call loclabel(10,"<EDF Primitive Centers>")
    read(10,*)
    read(10,*) b_EDF%center
    call loclabel(10,"<EDF Primitive Types>")
    read(10,*)
    read(10,*) b_EDF%functype !We assume all the type index is 1 (S type)
    if (maxval(b_EDF%functype)>1) then
        write(*,*) "ERROR: All GTFs of electron density function must be S type! Press ENTER to exit"
        read(*,*)
        stop
    end if
    call loclabel(10,"<EDF Primitive Exponents>")
    read(10,*)
    read(10,*) b_EDF%exp
    call loclabel(10,"<EDF Primitive Coefficients>")
    read(10,*)
    read(10,*) CO_EDF
    call loclabel(10,"<Number of Core Electrons>")
    read(10,*)
    read(10,*) ninnerelec
    if (infomode==0) write(*,"(a,i6,a)") "Note: Electron density functions (EDF) represent",ninnerelec," inner electrons"
end if
close(10)

if ( all(MOocc==int(MOocc)) ) then
    wfntype=2
    if (nmo==nelec) wfntype=1
    if (nmo==nelec/2) wfntype=0
else !post-HF
    if (naelec==nbelec) wfntype=3
    if (naelec/=nbelec) wfntype=4
end if
if (infomode==0) then
    write(*,*)
    write(*,"(' Total energy:',f19.12,' Hartree,   Virial ratio:',f12.8)") totenergy,virialratio
    write(*,"(' Total/Alpha/Beta electrons:',3f12.4)") nelec,naelec,nbelec
    write(*,"(' Number of orbital:',i6,',  Atoms:',i7,',  GTFs:',i7)") nmo,ncenter,nprims
    if (wfntype==0) write(*,"(' This is restricted close-shell single-determinant wavefunction')")
    if (wfntype==1) write(*,"(' This is unrestricted single-determinant wavefunction')")
    if (wfntype==2) write(*,"(' This is restricted open-shell wavefunction')")
    if (wfntype==3) write(*,"(' This is close-shell post-HF wavefunction')")
    if (wfntype==4) write(*,"(' This is open-shell post-HF wavefunction')")
    if (wfntype==1.or.wfntype==4) then
        do i=1,nmo
            if (MOtype(i)==2) exit
        end do
        write(*,"(' Orbitals from 1 to',i6,' are alpha type, from',i6,' to',i6,' are beta type')") i-1,i,nmo
    end if
    write(*,*)
end if
end subroutine



!!------- Load EDF information from external atomic .wfx files
subroutine readatmEDF
use defvar
use util
character elewfxfilename(110)*200,c200tmp*200
real*8,allocatable :: EDFCOtmp(:),EDFexptmp(:)
integer,allocatable :: EDFtypetmp(:)
integer atmsel(nelesupp,ncenter),natmsel(nelesupp)
iwfxtime=1
nEDFprims=0
ninnerele=0
do while(.true.)
    write(*,*) "Load the inner-core density (EDF information) for which element? e.g. Fe"
    write(*,*) "You can also input atomic indices, e.g. 5,8-10,31 means selecting 5,8,9,10,31"
    write(*,*) "Note: If finished, input ""q"""
    read(*,"(a)") c200tmp
    itmp=ichar(c200tmp(1:1))
    if (c200tmp=='q') then
        exit
    else if (itmp>=48.and.itmp<=57) then !Inputted is atomic indices
        call str2arr(c200tmp,natmsel(iwfxtime),atmsel(iwfxtime,:))
    else !Inputted is element name
        call lc2uc(c200tmp(1:1)) !Make the first/second character in upper/lower case
        call uc2lc(c200tmp(2:2))
        natmsel(iwfxtime)=0
        do iatm=1,ncenter
            if (a(iatm)%name==c200tmp(1:2)) then
                natmsel(iwfxtime)=natmsel(iwfxtime)+1
                atmsel(iwfxtime,natmsel(iwfxtime))=iatm
            end if
        end do
    end if
    if (natmsel(iwfxtime)==0) then
        write(*,*) "No atoms are selected, input again"
        write(*,*)
    else if (natmsel(iwfxtime)>0) then
        write(*,"(' The number of atoms selected is',i7,',  including:')") natmsel(iwfxtime)
        write(*,"(12i6)") atmsel(iwfxtime,1:natmsel(iwfxtime))
        write(*,*)
        write(*,*) "Load EDF information from which file? e.g. c:\ltwd\Fe_lanl2.wfx"
        do while(.true.)
            read(*,*) elewfxfilename(iwfxtime)
            inquire(file=elewfxfilename(iwfxtime),exist=alive)
            if (alive) exit
            write(*,*) "Cannot find this file, input again"
        end do
        open(10,file=elewfxfilename(iwfxtime),status="old") !Count how many EDF GTFs in this file
        call loclabel(10,"<Number of EDF Primitives>",ifound)
        if (ifound==0) then
            write(*,*) "Error: Unable to find EDF information from this file!"
            cycle
        end if
        read(10,*)
        read(10,*) nEDFtmp
        call loclabel(10,"<Number of Core Electrons>")
        read(10,*)
        read(10,*) ninnerelectmp
        close(10)
        write(*,"(' The number of EDF primitives in this file is',i5,/)") nEDFtmp
        nEDFprims=nEDFprims+nEDFtmp*natmsel(iwfxtime)
        ninnerele=ninnerele+ninnerelectmp*natmsel(iwfxtime)
        iwfxtime=iwfxtime+1
    end if
end do
nwfxtime=iwfxtime-1
write(*,"(' The total number of EDF primitives is',i7)") nEDFprims
write(*,"(' The total number of inner-core electrons represented by EDF is',i8)") ninnerele
allocate(b_EDF(nEDFprims),CO_EDF(nEDFprims))
ipos=1
do iwfxtime=1,nwfxtime
    open(10,file=elewfxfilename(iwfxtime),status="old")
    call loclabel(10,"<Number of EDF Primitives>",ifound)
    read(10,*)
    read(10,*) nEDFtmp
    allocate(EDFCOtmp(nEDFtmp),EDFexptmp(nEDFtmp),EDFtypetmp(nEDFtmp))
    call loclabel(10,"<EDF Primitive Types>")
    read(10,*)
    read(10,*) EDFtypetmp
    call loclabel(10,"<EDF Primitive Exponents>")
    read(10,*)
    read(10,*) EDFexptmp
    call loclabel(10,"<EDF Primitive Coefficients>")
    read(10,*)
    read(10,*) EDFCOtmp
    do iatm=1,natmsel(iwfxtime)
        b_EDF(ipos:ipos+nEDFtmp-1)%functype=EDFtypetmp
        b_EDF(ipos:ipos+nEDFtmp-1)%exp=EDFexptmp
        b_EDF(ipos:ipos+nEDFtmp-1)%center=atmsel(iwfxtime,iatm)
        CO_EDF(ipos:ipos+nEDFtmp-1)=EDFCOtmp
        ipos=ipos+nEDFtmp
    end do
    deallocate(EDFCOtmp,EDFexptmp,EDFtypetmp)
    close(10)
end do
! do iEDFprim=1,nEDFprims
!     write(*,"(2i6,2f18.8)") b_EDF(iEDFprim)%center,b_EDF(iEDFprim)%functype,b_EDF(iEDFprim)%exp,CO_EDF(iEDFprim)
! end do
write(*,*) "The EDF information have been loaded"
end subroutine



!!-----------------------------------------------------------------
!!--------- Read Molden input file, get coordinate, basis function and GTF information
!Known issue:
!CFour sometimes fail (e.g. benzene)
!ORCA result is inaccurate when g functions present
subroutine readmolden(name,infomode) !infomode=0 means output info, =1 silent
use defvar
use util
implicit real*8 (a-h,o-z)
character(len=*) name
character c80*80,symtmp*4
integer,allocatable :: shelltype(:),shellcon(:),shell2atom(:) !The definition of shelltype is identical to .fch
integer :: s2f(-5:5,21)=0 !Give shell type & orbital index to get functype
real*8,allocatable :: primexp(:),concoeff(:)
real*8,allocatable :: amocoeff(:,:),bmocoeff(:,:)
real*8 conv5d6d(6,5),conv7f10f(10,7),conv9g15g(15,9),conv11h21h(21,11)
real*8 conv5d6dtr(5,6),conv7f10ftr(7,10),conv9g15gtr(9,15),conv11h21htr(11,21)
!For backing up spherical basis functions
integer,allocatable :: shelltype5D(:)
real*8,allocatable :: CObasa5D(:,:),CObasb5D(:,:),Sbas5D(:,:),Dbas5D(:,:,:),Magbas5D(:,:,:)
real*8,external :: normgau
ifiletype=9
imodwfn=0
s2f(-5,1:11)=(/ -32,-31,-30,-29,-28,-27,-26,-25,-24,-23,-22 /)
s2f(-4,1:9)=(/ -21,-20,-19,-18,-17,-16,-15,-14,-13 /)
s2f(-3,1:7)=(/ -12,-11,-10,-9,-8,-7,-6 /)
s2f(-2,1:5)=(/ -5,-4,-3,-2,-1 /)
s2f(-1,1:4)=(/ 1,2,3,4 /)
s2f(0,1)=1
s2f(1,1:3)=(/ 2,3,4 /)
s2f(2,1:6)=(/ 5,6,7,8,9,10 /)
!---------- The sequence of f functions in Multiwfn (=wfn=wfx) is not identical to Molden, so convert here
!11  12  13  14  15  16  17  18  19  20  !Multiwfn sequence
!XXX YYY ZZZ XXY XXZ YYZ XYY XZZ YZZ XYZ
!xxx yyy zzz xyy xxy xxz xzz yzz yyz xyz !Molden sequence
s2f(3,1:10)=(/ 11,12,13,17,14,15,18,19,16,20 /)
!---------- The sequence of g functions in Multiwfn (=fch) is not identical to Molden, so convert here
! 21   22   23   24   25   26   27   28   29   30   31   32   33   34   35  !Multiwfn sequence
!ZZZZ YZZZ YYZZ YYYZ YYYY XZZZ XYZZ XYYZ XYYY XXZZ XXYZ XXYY XXXZ XXXY XXXX
!xxxx yyyy zzzz xxxy xxxz yyyx yyyz zzzx zzzy xxyy xxzz yyzz xxyz yyxz zzxy !Molden sequence
s2f(4,1:15)=(/ 35,25,21,34,33,29,24,26,22,32,30,23,31,28,27 /)
!---------- The sequence of h functions in Multiwfn (=fch=wfx) is not identical to Molden, so convert here
! 36    37    38    39    40    41    42    43    44    45    46  !Multiwfn sequence
!ZZZZZ YZZZZ YYZZZ YYYZZ YYYYZ YYYYY XZZZZ XYZZZ XYYZZ XYYYZ XYYYY
!xxxxx yyyyy zzzzz xxxxy xxxxz xyyyy xzzzz yyyyz yzzzz xxxyy xxxzz !Molden sequence
! 47    48    49    50    51    52    53    54    55    56  !Multiwfn sequence
!XXZZZ XXYZZ XXYYZ XXYYY XXXZZ XXXYZ XXXYY XXXXZ XXXXY XXXXX
!xxyyy xxzzz yyyzz yyzzz xxxyz xyyyz xyzzz xxyyz xxyzz xyyzz !Molden sequence
s2f(5,1:21)=(/ 56,41,36,55,54,46,42,40,37,53,51,50,47,39,38,52,45,43,49,48,44 /)

call gensphcartab(conv5d6d,conv7f10f,conv9g15g,conv11h21h)
conv5d6dtr=transpose(conv5d6d)
conv7f10ftr=transpose(conv7f10f)
conv9g15gtr=transpose(conv9g15g)
conv11h21htr=transpose(conv11h21h)

open(10,file=name,access="sequential",status="old")

if (infomode==0) write(*,*) "Loading various information of the wavefunction"
!!!!! Load atom information
call loclabel(10,"[Atoms]",ifound)
if (ifound==0) call loclabel(10,"[ATOMS]",ifound)
ilenunit=1 !Default length unit is a.u. =2 means Angstrom
read(10,"(a)") c80
if (index(c80,"Angs")/=0) ilenunit=2
ncenter=0
do while(.true.) !Count the number of atoms
    read(10,"(a)") c80
    if (c80==" ".or.index(c80,"[")/=0) then !Have passed atom field
        exit
    end if
    ncenter=ncenter+1
end do
if (allocated(a)) deallocate(a)
allocate(a(ncenter))
call loclabel(10,"[Atoms]",ifound) !Return to [Atoms]
if (ifound==0) call loclabel(10,"[ATOMS]",ifound)
read(10,"(a)") c80
do iatm=1,ncenter
    read(10,*) c80,nouse,a(iatm)%index,a(iatm)%x,a(iatm)%y,a(iatm)%z
end do
a%name=ind2name(a%index)
a%charge=a%index
if (ilenunit==2) then !Angstrom->a.u.
    a%x=a%x/b2a
    a%y=a%y/b2a
    a%z=a%z/b2a
end if

!Detect the Molden input file is produced by which program
rewind(10)
iorca=0
icfour=0
do while(.true.)
    read(10,"(a)") c80
    if (index(c80,"[GTO]")/=0) exit
    call struc2lc(c80)
    if (index(c80,"orca")/=0) then
        iorca=1
        exit
    else if (index(c80,"cfour")/=0) then
        icfour=1
        exit
    end if
end do
if (iorca==1.and.infomode==0) write(*,*) "This file is generated by ORCA! Special treatment is applied..."
if (icfour==1.and.infomode==0) write(*,*) "This file is generated by CFour! Special treatment is applied..."

!!!!! Load basis-set, and build up GTF information
if (infomode==0) write(*,*) "Loading basis-set definition..."
call loclabel(10,"[GTO]",ifound)
if (ifound==0) then
    write(*,*) "Error: Unable to locate [GTO] field! Press ENTER to exit"
    write(*,*) "Note: [STO] is currently not supported"
    read(*,*)
    stop
end if

!First time, we count how many shells are there to allocate proper size of allocatable arrays
nshell=0
nprimshell=0
read(10,*)
do while(.true.)
    read(10,*) iatm
    do while(.true.)
        read(10,*) c80,ncon
        nshell=nshell+1
        iaddnprmsh=0
        do ish=1,ncon
            read(10,*) tmpv1,tmpv2
            if (tmpv2/=0D0) iaddnprmsh=iaddnprmsh+1 !Many GTF shells outputted by Molpro have zero concontraction coefficient, they shouldn't be read in
        end do
        nprimshell=nprimshell+iaddnprmsh
        if (index(c80,"sp")/=0.or.index(c80,"SP")/=0) then !sp shell will be separated as s and p
            nshell=nshell+1
            nprimshell=nprimshell+iaddnprmsh
        end if
        read(10,"(a)") c80
        if (c80==" ") exit
        backspace(10)
    end do
    read(10,"(a)") c80
    if (c80==" ".or.index(c80,"[")/=0) exit !Finished reading [GTO] field
    backspace(10)
end do

!Second time, read basis-set information actually
allocate(shelltype(nshell),shellcon(nshell),shell2atom(nshell))
allocate(primexp(nprimshell),concoeff(nprimshell))
call loclabel(10,"[GTO]",ifound)
shellcon=0
ishell=0
iprimshell=0
read(10,*) 
do while(.true.)
    read(10,*) iatm
    do while(.true.)
        read(10,*) c80,ncon
        ishell=ishell+1
        shell2atom(ishell)=iatm
        !Determine shell type of basis function, here we first assume they are all Cartesian type
        if (index(c80,"sp")/=0.or.index(c80,"SP")/=0) then
            shelltype(ishell)=-1
        else if (index(c80,"s")/=0.or.index(c80,"S")/=0) then
            shelltype(ishell)=0
        else if (index(c80,"p")/=0.or.index(c80,"P")/=0) then
            shelltype(ishell)=1
        else if (index(c80,"d")/=0.or.index(c80,"D")/=0) then
            shelltype(ishell)=2
        else if (index(c80,"f")/=0.or.index(c80,"F")/=0) then
            shelltype(ishell)=3
        else if (index(c80,"g")/=0.or.index(c80,"G")/=0) then
            shelltype(ishell)=4
        else if (index(c80,"h")/=0.or.index(c80,"H")/=0) then
            shelltype(ishell)=5
        end if
        iprimshellold=iprimshell
        do ish=1,ncon !Read exponents and contraction coefficients. For SP, here load the S one
            read(10,*) exptmp,concoefftmp
            if (concoefftmp==0D0) cycle !The shell with zero contraction coefficients will be ripped out
            iprimshell=iprimshell+1
            shellcon(ishell)=shellcon(ishell)+1
            primexp(iprimshell)=exptmp
            if (iorca==1) then !ORCA doesn't present SP shell in Molden input file, so don't worry about -1
                !The normalization coefficients of spherical harmonic GTFs are weirdly multiplied into contraction coefficients, so they should be to retrieved to standard case
                rnorm=rnormgau_ORCA(primexp(iprimshell),shelltype(ishell))
                concoefftmp=concoefftmp/rnorm
            end if
            concoeff(iprimshell)=concoefftmp
!             write(*,"(2i4,3f18.10)") iatm,shelltype(ishell),primexp(iprimshell),concoefftmp,rnorm
        end do
        nprmshadd=iprimshell-iprimshellold
        !For ORCA, d,f,g are not properly normalized, e.g. at current stage d are normalized to 3, so renormalization is required
        !But the orbital coefficients always correspond for normalized ones (but some normalized to -1)
        if (shelltype(ishell)/=-1) call renormmoldengau(nprmshadd,shelltype(ishell),primexp(iprimshell-nprmshadd+1:iprimshell),concoeff(iprimshell-nprmshadd+1:iprimshell))
        
        if (shelltype(ishell)==-1) then !Separate SP shell as S and P shells
            shelltype(ishell)=0 !s
            call renormmoldengau(nprmshadd,shelltype(ishell),primexp(iprimshell-nprmshadd+1:iprimshell),concoeff(iprimshell-nprmshadd+1:iprimshell))
            ishell=ishell+1
            shelltype(ishell)=1 !p
            shellcon(ishell)=shellcon(ishell-1)
            shell2atom(ishell)=shell2atom(ishell-1)
            primexp(iprimshell+1:iprimshell+nprmshadd)=primexp(iprimshellold+1:iprimshell)
            do itmp=1,ncon !Backspace and load P contract coefficient
                backspace(10)
            end do
            do itmp=1,ncon
                read(10,*) exptmp,rnouse,concoefftmp
                if (concoefftmp==0D0) cycle
                iprimshell=iprimshell+1
                concoeff(iprimshell)=concoefftmp
            end do
            call renormmoldengau(nprmshadd,shelltype(ishell),primexp(iprimshell-nprmshadd+1:iprimshell),concoeff(iprimshell-nprmshadd+1:iprimshell))
        end if
        read(10,"(a)") c80
        if (c80==" ") exit
        backspace(10)
    end do
    read(10,"(a)") c80
    if (c80==" ".or.index(c80,"[")/=0) exit !Finished reading [GTO] field
    backspace(10)
end do

!Determine if the basis functions are Cartesian or spherical harmonic type. Admixture cartesian and spherical type are permitted
isphergau=0 !Default is Cartesian type
i5D=0
i10Flabel=0
i9G=0
i11H=0
imaxL=maxval(shelltype)
if (infomode==0) then
    if (imaxL==0) write(*,*) "The highest angular moment basis functions is S"
    if (imaxL==1) write(*,*) "The highest angular moment basis functions is P"
    if (imaxL==2) write(*,*) "The highest angular moment basis functions is D"
    if (imaxL==3) write(*,*) "The highest angular moment basis functions is F"
    if (imaxL==4) write(*,*) "The highest angular moment basis functions is G"
    if (imaxL==5) write(*,*) "The highest angular moment basis functions is H"
end if
if (imaxL>=2) then
    rewind(10)
    do while(.true.)
        read(10,"(a)") c80
        if (index(c80,'[5D')/=0.or.index(c80,'[5d')/=0) i5D=1
        if (index(c80,'10F')/=0) i10Flabel=1
        if (index(c80,'9G')/=0.or.index(c80,'9g')/=0) i9G=1
        if (index(c80,'11H')/=0.or.index(c80,'11h')/=0) i11H=1
        if (index(c80,'[MO]')/=0) exit
    end do
    if (i5D==1) then
        i7F=1 !By default, using 5D also implies 7F is used, unless 10F is explicitly specified
        if (i10Flabel==1) i7F=0
    end if
end if

if (i5D==1.or.i7F==1.or.i9G==1.or.i11H==1) isphergau=1
if (i5D==1) where(shelltype==2) shelltype=-2
if (i7F==1) where(shelltype==3) shelltype=-3
if (i9G==1) where(shelltype==4) shelltype=-4
if (i11H==1) where(shelltype==5) shelltype=-5
if (infomode==0) then
    if (isphergau==0) then
        write(*,*) "All basis functions are Cartesian type"
    else if (isphergau==1) then
        if (i5D==1.and.any(abs(shelltype)==2)) write(*,*) "All D basis functions are spherical harmonic type"
        if (i7F==1.and.any(abs(shelltype)==3)) write(*,*) "All F basis functions are spherical harmonic type"
        if (i9G==1.and.any(abs(shelltype)==4)) write(*,*) "All G basis functions are spherical harmonic type"
        if (i11H==1.and.any(abs(shelltype)==5)) write(*,*) "All H basis functions are spherical harmonic type"
    end if
end if
nbasis=0
do ishell=1,nshell
    nbasis=nbasis+type2norb(shelltype(ishell))
end do
if (infomode==0.and.iorca==1.and.imaxL>=4) then
    write(*,"(a)") " Warning! For the Molden input file generated by ORCA, when g angular moment basis functions present, the result may be inaccurate! &
    If you really want to preceed, press ENTER button"
    read (*,*)
end if

!!!!! Load orbital information. The sequence: Alpha(high occ / low ene) -> Alpha(low occ / high ene) -> Beta(high occ / low ene) -> Beta(low occ / high ene)
!Close shell orbitals are formally marked as "Alpha" spin. For singly occupied orbitals of ROHF, the spin marker are also Alpha
if (infomode==0) write(*,*) "Loading orbitals..."
nmo=nbasis
call loclabel(10,"Beta",ibeta)
! if (ibeta==0) call loclabel(10,"BETA",ibeta) !Some programs are strange, they don't follow standard molden format, use ALPHA/BETA instead of Alpha/Beta
if (ibeta==0) then
    nmo=nbasis
    allocate(amocoeff(nmo,nbasis),MOocc(nmo),MOene(nmo),MOtype(nmo),MOsym(nmo))
    amocoeff=0D0
else if (ibeta==1) then
    nmo=2*nbasis
    allocate(amocoeff(nbasis,nbasis),bmocoeff(nbasis,nbasis),MOocc(nmo),MOene(nmo),MOtype(nmo),MOsym(nmo))
    amocoeff=0D0
    bmocoeff=0D0
end if
MOsym=" "
MOocc=0D0
MOene=0D0
call loclabel(10,"[MO]",ifound)
read(10,*)
iMOa=0
iMOb=0
itmp=0
do while(.true.)
    itmp=itmp+1
    read(10,"(a)") c80 !Test if it is "Sym="
    backspace(10)
    if (index(c80,"Sym")/=0.or.index(c80,"SYM")/=0) then
        read(10,*) c80,symtmp
        !Remove digitals before the IRREP, e.g. 23B1 should be changed to B1
        do jtmp=1,len(symtmp)
            if (ichar(symtmp(jtmp:jtmp))<48.or.ichar(symtmp(jtmp:jtmp))>57) exit !Find the first position of non-digital
        end do
        symtmp=symtmp(jtmp:len(symtmp))
    else
        symtmp="?"
    end if
    read(10,*) c80,enetmp !Read orbital energy
!     write(*,*) itmp,c80,enetmp  !<<------ If encountering problem when loading MOs, using this to locate the problematic MO
    read(10,"(a)") c80 !Read orbital spin
    ispintmp=1 !Alpha
    if (index(c80,"Beta")/=0) ispintmp=2 !Beta
    read(10,*) c80,occtmp !Read orbital occupation number
    if (ispintmp==1) then
        iMOa=iMOa+1
        MOocc(iMOa)=occtmp
        MOene(iMOa)=enetmp
        MOsym(iMOa)=symtmp
        do ibasis=1,nbasis
!             write(*,*) ibasis,nbasis
            read(10,*) nouse,amocoeff(iMOa,ibasis)
        end do
    else
        iMOb=iMOb+1
        MOocc(nbasis+iMOb)=occtmp
        MOene(nbasis+iMOb)=enetmp
        MOsym(nbasis+iMOb)=symtmp
        do ibasis=1,nbasis
            read(10,*) nouse,bmocoeff(iMOb,ibasis)
        end do
    end if
    read(10,"(a)",iostat=ierror) c80 !Test if the ending of [MO] field is reached
    if (ierror/=0.or.c80==" ".or.index(c80,"[")/=0) exit
    backspace(10)
end do

!Fix orbital coefficients for ORCA. ORCA is rather rather frantic, the F(+3,-3) and G(+3,-3,+4,-4) in ORCA are normalized to -1 rather than 1,
!therefore the sign of their coefficients in all orbitals must be inverted! Hey ORCA, why did you do this!????? Totally non-understandable!
if (iorca==1) then
    ibasis=0
    do ishell=1,nshell
        if (shelltype(ishell)==-3) then
            amocoeff(:,ibasis+6:ibasis+7)=-amocoeff(:,ibasis+6:ibasis+7)
            if (ibeta==1) bmocoeff(:,ibasis+6:ibasis+7)=-bmocoeff(:,ibasis+6:ibasis+7)
        else if (shelltype(ishell)==-4) then
            amocoeff(:,ibasis+6:ibasis+9)=-amocoeff(:,ibasis+6:ibasis+9)
            if (ibeta==1) bmocoeff(:,ibasis+6:ibasis+9)=-bmocoeff(:,ibasis+6:ibasis+9)
        end if
        ibasis=ibasis+type2norb(shelltype(ishell))
    end do
end if
!Fix orbital coefficients for CFour according to Molden2aim. CFour only use Cartesian type basis function
!Notice that at current stage the GTO recording sequence has not been reordered, which is still identical to Molden sequence
!Also notice that even we do this, the result is still incorrect for e.g. benzene. But if we don't do this, we can't even obtain correct result for test case of Molden2aim
if (icfour==1) then
    ibasis=0
    do ishell=1,nshell
        if (shelltype(ishell)==2) then
            !xx, yy, zz, xy, xz, yz
            amocoeff(:,ibasis+1:ibasis+3)=amocoeff(:,ibasis+1:ibasis+3)*sqrt(3D0) !d(xx,yy,zz)*sqrt(3)
            if (ibeta==1) bmocoeff(:,ibasis+1:ibasis+3)=bmocoeff(:,ibasis+1:ibasis+3)*sqrt(3D0)
        else if (shelltype(ishell)==3) then
            !xxx yyy zzz xyy xxy xxz xzz yzz yyz xyz !Molden sequence
            amocoeff(:,ibasis+1:ibasis+3)=amocoeff(:,ibasis+1:ibasis+3)*sqrt(15D0) !f(xxx,yyy,zzz)*sqrt(15)
            amocoeff(:,ibasis+4:ibasis+9)=amocoeff(:,ibasis+4:ibasis+9)*sqrt(3D0) !f(xyy,xzz,yxx,yzz,zxx,zyy)*sqrt(3)
            if (ibeta==1) then
                bmocoeff(:,ibasis+1:ibasis+3)=bmocoeff(:,ibasis+1:ibasis+3)*sqrt(15D0)
                bmocoeff(:,ibasis+4:ibasis+9)=bmocoeff(:,ibasis+4:ibasis+9)*sqrt(3D0)
            end if
        else if (shelltype(ishell)==4) then
            !xxxx yyyy zzzz xxxy xxxz yyyx yyyz zzzx zzzy xxyy xxzz yyzz xxyz yyxz zzxy !Molden sequence
            amocoeff(:,ibasis+1:ibasis+3)=amocoeff(:,ibasis+1:ibasis+3)*sqrt(105D0) !g(x4,y4,z4)*sqrt(105)
            amocoeff(:,ibasis+4:ibasis+9)=amocoeff(:,ibasis+4:ibasis+9)*sqrt(15D0) !g(x3y,x3z,y3x,y3z,z3x,z3y)*sqrt(15)
            amocoeff(:,ibasis+10:ibasis+12)=amocoeff(:,ibasis+10:ibasis+12)*3D0 !g(x2y2,x2z2,y2z2)*3.0
            amocoeff(:,ibasis+13:ibasis+15)=amocoeff(:,ibasis+13:ibasis+15)*sqrt(3D0) !g(x2yz,y2xz,z2xy)*sqrt(3)
            if (ibeta==1) then
                bmocoeff(:,ibasis+1:ibasis+3)=bmocoeff(:,ibasis+1:ibasis+3)*sqrt(105D0)
                bmocoeff(:,ibasis+4:ibasis+9)=bmocoeff(:,ibasis+4:ibasis+9)*sqrt(15D0)
                bmocoeff(:,ibasis+10:ibasis+12)=bmocoeff(:,ibasis+10:ibasis+12)*3D0
                bmocoeff(:,ibasis+13:ibasis+15)=bmocoeff(:,ibasis+13:ibasis+15)*sqrt(3D0)
            end if
        end if
        ibasis=ibasis+type2norb(shelltype(ishell))
    end do    
end if

!Determine wavefunction type
if (ibeta==0) then
    MOtype=0 !Close shell orbital
    wfntype=0 !RHF
    if (any(MOocc/=nint(MOocc))) then
        wfntype=3 !R-post-HF
    else if (any(MOocc==1D0)) then
        wfntype=2 !ROHF
        do imo=1,nmo
            if (MOocc(imo)==1D0) MOtype(imo)=1
        end do
    end if
    if (infomode==0) write(*,"( ' The actual number of orbitals read:',i10)") iMOa
else if (ibeta==1) then
    wfntype=1 !UHF
    if (any(MOocc/=nint(MOocc))) wfntype=4 !U-post-HF
    MOtype(1:nbasis)=1
    MOtype(nbasis+1:nmo)=2
    if (infomode==0) write(*,"( ' The actual number of Alpha/Beta orbitals read:',i10,'  /',i10)") iMOa,iMOb
end if
call updatenelec !Cound the number of electrons


close(10)
!!!!!! All reading have finished, now generate basis information
!Below codes are adapted from readfch

!Backup spherical gauss basis information with 5D suffix (of course, may be 7f, 9g... in fact), convert them to cartesian type temporarily, 
!at final stage recover them back, namely get Sbas, Ptot... in spherical basis
if (isphergau==1) then
    allocate(shelltype5D(nshell))
    shelltype5D=shelltype
    where (shelltype<=-2) shelltype=-shelltype !Convert to cartesian type
    nbasis5D=nbasis
    nbasis=0
    do i=1,nshell
        nbasis=nbasis+type2norb(shelltype(i))
    end do
end if

!Allocate space for arrays
nprims=0
do i=1,nshell
    nprims=nprims+type2norb(shelltype(i))*shellcon(i)
end do
allocate(b(nprims),co(nmo,nprims),basshell(nbasis),bascen(nbasis),bastype(nbasis),primstart(nbasis),&
primend(nbasis),primconnorm(nprims),basstart(ncenter),basend(ncenter))
!Fill Cobasa and CObasb
if (isphergau==0) then
    allocate(CObasa(nbasis,nbasis))
    CObasa=transpose(amocoeff)
    if (wfntype==1.or.wfntype==4) then
        allocate(CObasb(nbasis,nbasis))
        CObasb=transpose(bmocoeff)
    end if
else if (isphergau==1) then !Since we have artifically convert spherical shells to cartesian shells, here the orbital coefficients are also correspondingly converted
    allocate(CObasa(nbasis,nbasis),CObasa5D(nbasis5D,nbasis5D))
    CObasa5D=transpose(amocoeff)
    CObasa=0D0
    if (wfntype==1.or.wfntype==4) then
        allocate(CObasb(nbasis,nbasis),CObasb5D(nbasis5D,nbasis5D))
        CObasb5D=transpose(bmocoeff)
        CObasb=0D0
    end if
    !Map 5D coefficient to 6D coefficient. Since the number of spherical basis functions is more than cartesian ones, 
    !therefore Cobasa (6D) will have some orbitals with vacant coefficients, only orbitals (1~nbasis5D) are filled
    ipos5D=1
    ipos6D=1
    do ish=1,nshell
        ishtyp5D=shelltype5D(ish)
        ishtyp6D=shelltype(ish)
        numshorb5D=type2norb(ishtyp5D)
        numshorb6D=type2norb(ishtyp6D)
        if (ishtyp5D>=-1) then !S or P or SP or other cartesian shells, directly copy
            CObasa(ipos6D:ipos6D+numshorb6D-1,1:nbasis5D)=CObasa5D(ipos5D:ipos5D+numshorb5D-1,:)
            if (wfntype==1.or.wfntype==4) CObasb(ipos6D:ipos6D+numshorb6D-1,1:nbasis5D)=CObasb5D(ipos5D:ipos5D+numshorb5D-1,:)            
        else if (ishtyp5D==-2) then
            !5D->6D
            CObasa(ipos6D:ipos6D+5,1:nbasis5D)=matmul(conv5d6d,CObasa5D(ipos5D:ipos5D+4,:))
            if (wfntype==1.or.wfntype==4) CObasb(ipos6D:ipos6D+5,1:nbasis5D)=matmul(conv5d6d,CObasb5D(ipos5D:ipos5D+4,:))
        else if (ishtyp5D==-3) then
            !7F->10F
            CObasa(ipos6D:ipos6D+9,1:nbasis5D)=matmul(conv7f10f,CObasa5D(ipos5D:ipos5D+6,:))
            if (wfntype==1.or.wfntype==4) CObasb(ipos6D:ipos6D+9,1:nbasis5D)=matmul(conv7f10f,CObasb5D(ipos5D:ipos5D+6,:))
        else if (ishtyp5D==-4) then
            !9G->15G
            CObasa(ipos6D:ipos6D+14,1:nbasis5D)=matmul(conv9g15g,CObasa5D(ipos5D:ipos5D+8,:))
            if (wfntype==1.or.wfntype==4) CObasb(ipos6D:ipos6D+14,1:nbasis5D)=matmul(conv9g15g,CObasb5D(ipos5D:ipos5D+8,:))
        else if (ishtyp5D==-5) then
            !11H->21H
            CObasa(ipos6D:ipos6D+20,1:nbasis5D)=matmul(conv11h21h,CObasa5D(ipos5D:ipos5D+10,:))
            if (wfntype==1.or.wfntype==4) CObasb(ipos6D:ipos6D+20,1:nbasis5D)=matmul(conv11h21h,CObasb5D(ipos5D:ipos5D+10,:))
        end if
        ipos5D=ipos5D+numshorb5D
        ipos6D=ipos6D+numshorb6D
    end do
end if

if (infomode==0) write(*,*) "Converting basis function information to GTF information..."
!Distribute exponent, functype to every GTF and generate CO(:,:) from amocoeff/bmocoeff
!Fill: b,basshell,bascen,bastype,co,primstart,primend,primconnorm
k=1 !current position of GTF
iexp=1
ibasis=1 !current position of basis
!Note: Below commented with !!! means the line associated to setting basis information
do i=1,nshell !cycle each basis shell
    b(k:k+shellcon(i)*type2norb(shelltype(i))-1)%center=shell2atom(i)
    basshell(ibasis:ibasis+type2norb(shelltype(i))-1)=i !!! set basis attributed to which shell
    bascen(ibasis:ibasis+type2norb(shelltype(i))-1)=shell2atom(i) !!! set basis attributed to which center
    do j=1,type2norb(shelltype(i)) !cycle each basis function in each basis shell
        b(k:k+shellcon(i)-1)%functype=s2f(shelltype(i),j)
        bastype(ibasis)=s2f(shelltype(i),j) !!! set basis type
        primstart(ibasis)=k !!! From where the GTFs attributed to ibasis'th basis
        primend(ibasis)=k+shellcon(i)-1 !!! To where the GTFs attributed to ibasis'th basis
        do l=1,shellcon(i) !cycle each GTF in each basis function
            b(k)%exp=primexp(iexp+l-1)
            tnormgau=normgau(b(k)%functype,b(k)%exp)  !!!Normalization coefficient of cartesian GTFs
            temp=concoeff(iexp+l-1)  !!!Contraction coefficient of GTFs
            primconnorm(k)=temp*tnormgau !Combines contraction and normalization coefficient
            do imo=1,nmo
                if (wfntype==0.or.wfntype==2.or.wfntype==3) then
                    co(imo,k)=cobasa(ibasis,imo)*temp*tnormgau
                else if (wfntype==1.or.wfntype==4) then
                    if (isphergau==1) then
                        if (imo<=nbasis5D) co(imo,k)=cobasa(ibasis,imo)*temp*tnormgau
                        if (imo>nbasis5D) co(imo,k)=cobasb(ibasis,imo-nbasis5D)*temp*tnormgau
                    else
                        if (imo<=nbasis) co(imo,k)=cobasa(ibasis,imo)*temp*tnormgau
                        if (imo>nbasis) co(imo,k)=cobasb(ibasis,imo-nbasis)*temp*tnormgau
                    end if
                end if
            end do
            k=k+1
        end do
        ibasis=ibasis+1
    end do
    iexp=iexp+shellcon(i)
end do

!Generate overlap matrix and dipole moment integral matrix for Cartesian Gauss basis functions
if (infomode==0) write(*,*) "Generating overlap matrix..."
allocate(Sbas(nbasis,nbasis))
call genSbas
if (igenDbas==1) then
    if (infomode==0) write(*,*) "Generating electric dipole moment integral matrix..."
    allocate(Dbas(3,nbasis,nbasis))
    call genDbas
end if
if (igenMagbas==1) then
    if (infomode==0) write(*,*) "Generating magnetic dipole moment integral matrix..."
    allocate(Magbas(3,nbasis,nbasis))
    call genMagbas
end if

!Check normalizaiton of basis functions
! do i=1,size(sbas,1)
!     write(*,"(i10,f12.6)") i,sbas(i,i)
! end do

if (isphergau==1) then
    if (infomode==0) write(*,*) "Converting basis function information from Cartesian to spherical type..."
    !Map cartesian overlap matrix to spherical harmonic overlap matrix
    allocate(Sbas5D(nbasis5D,nbasis5D))
    if (igenDbas==1) allocate(Dbas5D(3,nbasis5D,nbasis5D))
    if (igenMagbas==1) allocate(Magbas5D(3,nbasis5D,nbasis5D))
    ipos5D=1
    ipos6D=1
    do ish=1,nshell
        ishtyp5D=shelltype5D(ish)
        ishtyp6D=shelltype(ish)
        numshorb5D=type2norb(ishtyp5D)
        numshorb6D=type2norb(ishtyp6D)
        !Now contract columns
        if (ishtyp5D>=-1) sbas(:,ipos5D:ipos5D+numshorb5D-1)=sbas(:,ipos6D:ipos6D+numshorb6D-1) !S, P, SP or other Cartesian shells
        if (ishtyp5D==-2) sbas(:,ipos5D:ipos5D+numshorb5D-1)=matmul(sbas(:,ipos6D:ipos6D+numshorb6D-1),conv5d6d)
        if (ishtyp5D==-3) sbas(:,ipos5D:ipos5D+numshorb5D-1)=matmul(sbas(:,ipos6D:ipos6D+numshorb6D-1),conv7f10f)
        if (ishtyp5D==-4) sbas(:,ipos5D:ipos5D+numshorb5D-1)=matmul(sbas(:,ipos6D:ipos6D+numshorb6D-1),conv9g15g)
        if (ishtyp5D==-5) sbas(:,ipos5D:ipos5D+numshorb5D-1)=matmul(sbas(:,ipos6D:ipos6D+numshorb6D-1),conv11h21h)
        !Now contract rows
        if (ishtyp5D>=-1) sbas(ipos5D:ipos5D+numshorb5D-1,:)=sbas(ipos6D:ipos6D+numshorb6D-1,:) !S, P, SP or other Cartesian shells
        if (ishtyp5D==-2) sbas(ipos5D:ipos5D+numshorb5D-1,:)=matmul(conv5d6dtr,sbas(ipos6D:ipos6D+numshorb6D-1,:))
        if (ishtyp5D==-3) sbas(ipos5D:ipos5D+numshorb5D-1,:)=matmul(conv7f10ftr,sbas(ipos6D:ipos6D+numshorb6D-1,:))
        if (ishtyp5D==-4) sbas(ipos5D:ipos5D+numshorb5D-1,:)=matmul(conv9g15gtr,sbas(ipos6D:ipos6D+numshorb6D-1,:))
        if (ishtyp5D==-5) sbas(ipos5D:ipos5D+numshorb5D-1,:)=matmul(conv11h21htr,sbas(ipos6D:ipos6D+numshorb6D-1,:))
        
        if (igenDbas==1) then
            do idir=1,3
                !Now contract columns of Dbas
                if (ishtyp5D>=-1) Dbas(idir,:,ipos5D:ipos5D+numshorb5D-1)=Dbas(idir,:,ipos6D:ipos6D+numshorb6D-1) !S, P, SP or other Cartesian shells
                if (ishtyp5D==-2) Dbas(idir,:,ipos5D:ipos5D+numshorb5D-1)=matmul(Dbas(idir,:,ipos6D:ipos6D+numshorb6D-1),conv5d6d)
                if (ishtyp5D==-3) Dbas(idir,:,ipos5D:ipos5D+numshorb5D-1)=matmul(Dbas(idir,:,ipos6D:ipos6D+numshorb6D-1),conv7f10f)
                if (ishtyp5D==-4) Dbas(idir,:,ipos5D:ipos5D+numshorb5D-1)=matmul(Dbas(idir,:,ipos6D:ipos6D+numshorb6D-1),conv9g15g)
                if (ishtyp5D==-5) Dbas(idir,:,ipos5D:ipos5D+numshorb5D-1)=matmul(Dbas(idir,:,ipos6D:ipos6D+numshorb6D-1),conv11h21h)
                !Now contract rows of Dbas
                if (ishtyp5D>=-1) Dbas(idir,ipos5D:ipos5D+numshorb5D-1,:)=Dbas(idir,ipos6D:ipos6D+numshorb6D-1,:) !S, P, SP or other Cartesian shells
                if (ishtyp5D==-2) Dbas(idir,ipos5D:ipos5D+numshorb5D-1,:)=matmul(conv5d6dtr,Dbas(idir,ipos6D:ipos6D+numshorb6D-1,:))
                if (ishtyp5D==-3) Dbas(idir,ipos5D:ipos5D+numshorb5D-1,:)=matmul(conv7f10ftr,Dbas(idir,ipos6D:ipos6D+numshorb6D-1,:))
                if (ishtyp5D==-4) Dbas(idir,ipos5D:ipos5D+numshorb5D-1,:)=matmul(conv9g15gtr,Dbas(idir,ipos6D:ipos6D+numshorb6D-1,:))
                if (ishtyp5D==-5) Dbas(idir,ipos5D:ipos5D+numshorb5D-1,:)=matmul(conv11h21htr,Dbas(idir,ipos6D:ipos6D+numshorb6D-1,:))
            end do
        end if
        if (igenMagbas==1) then
            do idir=1,3
                !Now contract columns of Magbas
                if (ishtyp5D>=-1) Magbas(idir,:,ipos5D:ipos5D+numshorb5D-1)=Magbas(idir,:,ipos6D:ipos6D+numshorb6D-1) !S, P, SP or other Cartesian shells
                if (ishtyp5D==-2) Magbas(idir,:,ipos5D:ipos5D+numshorb5D-1)=matmul(Magbas(idir,:,ipos6D:ipos6D+numshorb6D-1),conv5d6d)
                if (ishtyp5D==-3) Magbas(idir,:,ipos5D:ipos5D+numshorb5D-1)=matmul(Magbas(idir,:,ipos6D:ipos6D+numshorb6D-1),conv7f10f)
                if (ishtyp5D==-4) Magbas(idir,:,ipos5D:ipos5D+numshorb5D-1)=matmul(Magbas(idir,:,ipos6D:ipos6D+numshorb6D-1),conv9g15g)
                if (ishtyp5D==-5) Magbas(idir,:,ipos5D:ipos5D+numshorb5D-1)=matmul(Magbas(idir,:,ipos6D:ipos6D+numshorb6D-1),conv11h21h)
                !Now contract rows of Magbas
                if (ishtyp5D>=-1) Magbas(idir,ipos5D:ipos5D+numshorb5D-1,:)=Magbas(idir,ipos6D:ipos6D+numshorb6D-1,:) !S, P, SP or other Cartesian shells
                if (ishtyp5D==-2) Magbas(idir,ipos5D:ipos5D+numshorb5D-1,:)=matmul(conv5d6dtr,Magbas(idir,ipos6D:ipos6D+numshorb6D-1,:))
                if (ishtyp5D==-3) Magbas(idir,ipos5D:ipos5D+numshorb5D-1,:)=matmul(conv7f10ftr,Magbas(idir,ipos6D:ipos6D+numshorb6D-1,:))
                if (ishtyp5D==-4) Magbas(idir,ipos5D:ipos5D+numshorb5D-1,:)=matmul(conv9g15gtr,Magbas(idir,ipos6D:ipos6D+numshorb6D-1,:))
                if (ishtyp5D==-5) Magbas(idir,ipos5D:ipos5D+numshorb5D-1,:)=matmul(conv11h21htr,Magbas(idir,ipos6D:ipos6D+numshorb6D-1,:))
            end do
        end if
        ipos5D=ipos5D+numshorb5D
        ipos6D=ipos6D+numshorb6D
    end do
    Sbas5D=sbas(1:nbasis5D,1:nbasis5D)
    if (igenDbas==1) Dbas5D=Dbas(:,1:nbasis5D,1:nbasis5D)
    if (igenMagbas==1) Magbas5D=Magbas(:,1:nbasis5D,1:nbasis5D)

    !Recover spherical gauss basis function information
    nbasis=nbasis5D
    shelltype=shelltype5D
    ibasis=1
    do i=1,nshell
        basshell(ibasis:ibasis+type2norb(shelltype(i))-1)=i
        bascen(ibasis:ibasis+type2norb(shelltype(i))-1)=shell2atom(i)
        do j=1,type2norb(shelltype(i))
            bastype(ibasis)=s2f(shelltype(i),j)
            ibasis=ibasis+1
        end do
    end do
    deallocate(CObasa)
    allocate(CObasa(nbasis,nbasis))
    CObasa=CObasa5D
    if (wfntype==1.or.wfntype==4) then
        deallocate(CObasb)
        allocate(CObasb(nbasis,nbasis))
        CObasb=CObasb5D
    end if
    deallocate(sbas)
    allocate(sbas(nbasis,nbasis))
    sbas=Sbas5D
    if (igenDbas==1) then
        deallocate(Dbas)
        allocate(Dbas(3,nbasis,nbasis))
        Dbas=Dbas5D
    end if
    if (igenMagbas==1) then
        deallocate(Magbas)
        allocate(Magbas(3,nbasis,nbasis))
        Magbas=Magbas5D
    end if
end if

!The orbitals in Molden input file may be not properly normalized, we normalize the coefficients here
do imo=1,nbasis !For open-shell, only check alpha part
    testnorm=0D0
    do ibas=1,nbasis
        do jbas=1,nbasis
            testnorm=testnorm+sbas(ibas,jbas)*cobasa(ibas,imo)*cobasa(jbas,imo)
        end do
    end do
!     write(*,*) imo,testnorm !If the orbitals are not normalized to 1, that means current Molden input file is problematic
    ! For ORCA, we have warned users because of g functions. Sometimes the actual orbitals are more than nbasis, then that orbitals will be normalized to zero
    if (abs(testnorm-nint(testnorm))>0.01.and.iorca==0) then
        write(*,"(/,a)") " Warning: The normalization check of orbital wavefunctions failed! That means this Molden input file cannot be well supported by Multiwfn. &
        (Because the Molden input files generated by many quantum chemistry programs are problematic)"
        write(*,"(a)") " If you really want to proceed, press ENTER button, but notice that the analysis result may be not correct"
        read (*,*)
        exit
    end if
    !Artificially normalize orbitals
!     if (imo>nbasis) then !beta part
!         cobasb(:,imo-nbasis)=cobasb(:,imo-nbasis)/dsqrt(testnorm)
!     else
!         cobasa(:,imo)=cobasa(:,imo)/dsqrt(testnorm)
!     end if
!     co(imo,:)=co(imo,:)/dsqrt(testnorm)
end do

!Move local shell arrays to the global ones, they will be used in other functions
allocate(shtype(nshell),shcen(nshell),shcon(nshell),primshexp(nprimshell),primshcoeff(nprimshell))
shtype=shelltype
shcen=shell2atom
shcon=shellcon
primshexp=primexp
primshcoeff=concoeff

!Generate basstart and basend
nowcen=0
indcen=0
do ibasis=1,nbasis
    if (bascen(ibasis)/=nowcen) then
        nowcen=bascen(ibasis)
        indcen=indcen+1
        basstart(indcen)=ibasis
        if (indcen/=1) basend(indcen-1)=ibasis-1
    end if
end do
basend(ncenter)=nbasis

!Generate one-particle density matrix for basis functions
if (igenP==1) then
    if (infomode==0) write(*,*) "Generating density matrix..."
    call genP
end if

!Output summary of present wavefunction
if (infomode==0) then
    write(*,*)
    write(*,"(' Total/Alpha/Beta electrons:',3f12.4)") nelec,naelec,nbelec
    write(*,"(' Net charge:',f12.5,'      Expected multiplicity:',i5)") sum(a(:)%charge)-nelec,nint(naelec-nbelec)+1
    write(*,"(' Atoms:',i7,',  Basis functions:',i7,',  GTFs:',i7)") ncenter,nmo,nprims
    if (wfntype==0) then
        write(*,"(' This is restricted single-determinant wavefunction')")
        write(*,"(' Orbitals from 1 to',i6,' are occupied')") nint(nelec/2)
    else if (wfntype==1) then
        write(*,"(' This is unrestricted single-determinant wavefunction')")
        write(*,"(' Orbitals from ',i6,' to',i6,' are alpha, from',i6,' to',i6,' are occupied')") 1,nbasis,1,nint(naelec)
        write(*,"(' Orbitals from ',i6,' to',i6,' are beta,  from',i6,' to',i6,' are occupied')") nbasis+1,nmo,nbasis+1,nbasis+nint(nbelec)
    else if (wfntype==2) then
        write(*,"(' This is restricted open-shell wavefunction')")
        write(*,"(' Orbitals from',i6,' to',i6,' are doubly occupied')") 1,nint(nbelec)
        write(*,"(' Orbitals from',i6,' to',i6,' are singly occupied')") nint(nbelec)+1,nint(naelec)
    else if (wfntype==3) then
        write(*,"(' This is restricted post-HF wavefunction')")
    else if (wfntype==4) then
        write(*,"(' This is unrestricted open-shell wavefunction')")
        write(*,"(' Orbitals from ',i6,' to',i6,' are alpha, from',i6,' to',i6,' are beta')") 1,nbasis,nbasis+1,nmo
    end if
end if

!Find out index of HOMO, will be used in some cases, only for RHF
if (wfntype==0) then
    do idxHOMO=nmo,1,-1
        if (nint(MOocc(idxHOMO))==2D0) exit
    end do
end if
end subroutine







!=======================================================================
!=======================================================================
!!!!!!!!!!!!!!!! Below routines are used to output files !!!!!!!!!!!!!!!
!=======================================================================
!=======================================================================



!!---------- Output current coordinate to pdb file
subroutine outpdb(outpdbname,ifileid)
use defvar
character(len=*) outpdbname
integer i,ifileid
open(ifileid,file=outpdbname,status="replace")
write(ifileid,"('REMARK   Generated by Multiwfn, Totally',i10,' atoms')") ncenter
do i=1,ncenter
    write(ifileid,"(a6,i5,1x,a4,1x,a3, 1x,a1,i4,4x,3f8.3,2f6.2,10x,a2)") &
    "HETATM",i,' '//ind2name_up(a(i)%index)//' ',"MOL",'A',1,a(i)%x*b2a,a(i)%y*b2a,a(i)%z*b2a,1.0,0.0,adjustr(ind2name_up(a(i)%index))
end do
write(ifileid,"('END')")
close(ifileid)
write(*,*) "Exporting pdb file finished!"
end subroutine


!!---------- Output current coordinate to xyz file
subroutine outxyz(outxyzname,ifileid)
use defvar
character(len=*) outxyzname
integer i,ifileid
open(ifileid,file=outxyzname,status="replace")
write(ifileid,"(i6)") ncenter
write(ifileid,*) "Generated by Multiwfn"
do i=1,ncenter
    write(ifileid,"(a,3f16.8)") ind2name_up(a(i)%index),a(i)%x*b2a,a(i)%y*b2a,a(i)%z*b2a
end do
close(ifileid)
write(*,*) "Exporting xyz file finished!"
end subroutine


!!---------- Output current coordinate to Gaussian input file
subroutine outgjf(outgjfname,ifileid)
use defvar
character(len=*) outgjfname
open(ifileid,file=outgjfname,status="replace")
write(ifileid,"(a,/,/,a,/)") "#P B3LYP/6-31G*","Generated by Multiwfn"
netcharge=nint(sum(a%charge)-nelec)
if (nelec==0) netcharge=0 !nelec==0 means no electron informations, e.g. pdb file
write(ifileid,"(2i2)") netcharge,nint(naelec-nbelec)+1
do i=1,ncenter
    write(ifileid,"(a,1x,3f14.8)") a(i)%name,a(i)%x*b2a,a(i)%y*b2a,a(i)%z*b2a
end do
write(ifileid,*)  !Two blank line at the end of the file
close(ifileid)
write(*,*) "Exporting Gaussian input file finished!"
end subroutine


!!---------- Output current coordinate to GAMESS-US input file
subroutine outGAMESSinp(outname,ifileid)
use defvar
character(len=*) outname
character SCFTYPE*3,selectyn
ioutguess=0
if (allocated(CObasa)) then
    write(*,*) "If write initial guess information? (y/n)"
    read(*,*) selectyn
    if (selectyn=='y'.or.selectyn=='Y') ioutguess=1
end if
netcharge=nint(sum(a%charge)-nelec)
if (nelec==0) netcharge=0 !nelec==0 means no electron informations, e.g. pdb file
mult=nint(naelec-nbelec)+1
iopsh=0
if (mult/=1.or.wfntype==1.or.wfntype==4) iopsh=1
SCFTYPE="RHF"
if (iopsh==1) SCFTYPE="UHF"
open(ifileid,file=outname,status="replace")
write(ifileid,"(a,i1,a,i1,a)") " $CONTRL SCFTYP="//SCFTYPE//" MULT=",mult," ICHARG=",netcharge," RUNTYP=ENERGY"
write(ifileid,"(a)") " DFTTYP=B3LYPV3 ISPHER=0 MAXIT=60 NPRINT=-5 $END"
write(ifileid,"(a)") " $BASIS GBASIS=N31 NGAUSS=6 NDFUNC=1 NPFUNC=0 DIFFSP=.F. DIFFS=.F. $END"
write(ifileid,"(a)") " $SYSTEM MWORDS=800 $END"
!  $lmoeda matom(1)=3,4 mcharg(1)=0,0 mmult(1)=1,1 $end
write(ifileid,"(a)") " $SCF DIRSCF=.T. $END"
write(ifileid,"(a)") " $DFT DC=.F. $END"
if (ioutguess==1) write(ifileid,"(a,i5,a)") " $GUESS GUESS=MOREAD NORB=",nbasis," $END"
write(ifileid,"(a)") " $DATA"
write(ifileid,"(a)") "Generated by Multiwfn"
write(ifileid,"(a)") "C1"
do i=1,ncenter
    write(ifileid,"(a,f4.1,1x,3f14.8)") a(i)%name,dfloat(a(i)%index),a(i)%x*b2a,a(i)%y*b2a,a(i)%z*b2a
end do
write(ifileid,"(a)") " $END"
! Write $VEC information
if (ioutguess==1) then
    write(ifileid,*)
    write(ifileid,"(a)") " $VEC"
    ntime=ceiling(nbasis/5D0)
    idxmo=1
    do imo=1,nbasis
        do itime=1,ntime
            if (itime<ntime) then
                write(ifileid,"(i2,i3,5(1PE15.8))") idxmo,itime,CObasa((itime-1)*5+1:itime*5,imo)
            else
                write(ifileid,"(i2,i3,5(1PE15.8))") idxmo,itime,CObasa((itime-1)*5+1:nbasis,imo)
            end if
        end do
        idxmo=idxmo+1
        if (idxmo==100) idxmo=0
    end do
    if (iopsh==1) then
        idxmo=1
        do imo=1,nbasis
            do itime=1,ntime
                if (itime<ntime) then
                    write(ifileid,"(i2,i3,5(1PE15.8))") idxmo,itime,CObasb((itime-1)*5+1:itime*5,imo)
                else
                    write(ifileid,"(i2,i3,5(1PE15.8))") idxmo,itime,CObasb((itime-1)*5+1:nbasis,imo)
                end if
            end do
            idxmo=idxmo+1
            if (idxmo==100) idxmo=0
        end do
    end if
    write(ifileid,"(a)") " $END"
end if
close(ifileid)
write(*,*) "Exporting GAMESS-US input file finished!"
write(*,*) "The task corresponds to single point at B3LYP-D3/6-31G* level"
end subroutine


!!!------------------------- Output 3D matrix with property to a cube file. fileid must be opened before invoking this routine, and close it after that
subroutine outcube(matrix,numx,numy,numz,org_x,org_y,org_z,transx,transy,transz,fileid)
use defvar
implicit real*8 (a-h,o-z)
integer numx,numy,numz,fileid
real*8 org_x,org_y,org_z,transx,transy,transz
real*8 matrix(numx,numy,numz)
write(fileid,"(' Generated by Multiwfn')")
write(fileid,"(' Totally ',i12,' grid points')") numx*numy*numz
if (ncenter>=1) then
    write(fileid,"(i5,3f12.6)") ncenter,org_x,org_y,org_z
else
    write(fileid,"(i5,3f12.6)") 1,org_x,org_y,org_z
end if
write(fileid,"(i5,3f12.6)") numx,transx,0.0,0.0
write(fileid,"(i5,3f12.6)") numy,0.0,transy,0.0
write(fileid,"(i5,3f12.6)") numz,0.0,0.0,transz
if (ncenter>=1) then
    do i=1,ncenter
        write(fileid,"(i5,4f12.6)") a(i)%index,a(i)%charge,a(i)%x,a(i)%y,a(i)%z
    end do
else
    write(*,"(a)") " Note: Current system has no atom, in order to maximize compatibility of the generated .cub file, a hydrogen atom is added to 0,0,0"
    write(*,*)
    write(fileid,"(i5,4f12.6)") 1,1D0,0D0,0D0,0D0
end if
where (abs(matrix)<=1D-99) matrix=0D0 !Diminish too small value, otherwise the symbol "E" cannot be shown by 1PE13.5 format e.g. 9.39376-116, 
write(*,*) "Please wait..."
do i=1,numx
    do j=1,numy
        write(fileid,"(6(1PE13.5))",advance="no") matrix(i,j,1:numz)
        write(fileid,*)
    end do
end do
end subroutine


!!!------------------------- Output current wavefunction to new.wfn
!If isortatmind==1, then any atom without GTF posited on it will not be output, and the index is filled to assure contiguous
!Orbitals with zero occupiation will not be outputted
!If ioutinfo==1, output information
subroutine outwfn(outwfnname,isortatmind,ioutinfo,ifileid)
use defvar
implicit real*8 (a-h,o-z)
character(len=*) outwfnname
integer isortatmind,ioutinfo,ifileid,indconv(ncenter)
!convGmul2wfn converts the g sequence used internally in Multiwfn (input) to commonly used g sequence in .wfn (the one outputted by Molden2AIM and g09 since B01 )
integer :: convGmul2wfn(35)=(/ (0,i=1,20), 23,29,32,27,22, 28,35,34,26,31, 33,30,25,24,21 /)
integer :: convGwfn2mul(35)=(/ (0,i=1,20), 35,25,21,34,33, 29,24,26,22,32, 30,23,31,28,27 /)

open(ifileid,file=outwfnname,status="replace")
write(ifileid,*) "Generated by Multiwfn"

if (isortatmind==1) then !Find real number of centers
    j=0
    do i=1,ncenter
        if (any(b%center==i)) j=j+1
    end do
    write(ifileid,"('GAUSSIAN',i15,' MOL ORBITALS',i7,' PRIMITIVES',i9,' NUCLEI')") count(MOocc(1:nmo)/=0D0),nprims,j
else
    write(ifileid,"('GAUSSIAN',i15,' MOL ORBITALS',i7,' PRIMITIVES',i9,' NUCLEI')") count(MOocc(1:nmo)/=0D0),nprims,ncenter
end if

j=1
do i=1,ncenter
    if (isortatmind==1) then
        if (all(b%center/=i)) cycle
    end if
    indconv(j)=i !The j actual atom corresponds to the i original atom
    write(ifileid,"(2x,a2,i4,4x,'(CENTRE',i3,')',1x,3f12.8,'  CHARGE =',f5.1)") a(i)%name,j,j,a(i)%x,a(i)%y,a(i)%z,a(i)%charge
    j=j+1
end do

if (isortatmind==1) then !Convert center of GTF to reordered index
    do i=1,j-1 !Cycle the centers with GTF
        where(b%center==indconv(i)) b%center=i
    end do
end if

!Convert the g sequence to common sequence in .wfn, output, and then convert back
write(ifileid,"('CENTRE ASSIGNMENTS  ',20i3)") b(1:nprims)%center
do i=1,nprims
    if (b(i)%functype>=21.and.b(i)%functype<=35) b(i)%functype=convGmul2wfn(b(i)%functype)
end do
write(ifileid,"('TYPE ASSIGNMENTS    ',20i3)") b(1:nprims)%functype
do i=1,nprims
    if (b(i)%functype>=21.and.b(i)%functype<=35) b(i)%functype=convGwfn2mul(b(i)%functype)
end do

write(ifileid,"('EXPONENTS ',5D14.7)") b(1:nprims)%exp
imo=0
izeroocc=0
do i=1,nmo
    if (MOocc(i)==0D0) then
        izeroocc=izeroocc+1
        cycle
    end if
    imo=imo+1 !Use i instead of imo, this can make MO index contiguous
    write(ifileid,"('MO',I5,'     MO 0.0        OCC NO = ',f12.7,'  ORB. ENERGY =', f12.6)") imo,MOocc(i),MOene(i)
    write(ifileid,"(5D16.8)") (co(i,j),j=1,nprims)
end do
write(ifileid,"('END DATA',/,' THE  HF ENERGY = ',f19.12,' THE VIRIAL(-V/T)= ',f12.8)") totenergy,virialratio
close(ifileid)
if (izeroocc>0.and.ioutinfo==1) write(*,"(a,i10,a)") " Note: Found",izeroocc," zero occupied orbitals and have discarded them"
end subroutine


!!!------------------------- Output current wavefunction to Molden input file
subroutine outmolden(outname,ifileid)
use defvar
use util
implicit real*8 (a-h,o-z)
character(len=*) outname
integer ifileid
character symbol

open(ifileid,file=outname,status="replace")
write(ifileid,"(a)") "[Molden Format]"
write(ifileid,"(a)") "[Atoms] AU"
do i=1,ncenter
    write(ifileid,"(a,i7,i4,3f14.7)") a(i)%name,i,a(i)%index,a(i)%x,a(i)%y,a(i)%z
end do
write(ifileid,"(a)") "[GTO]"
do iatm=1,ncenter
    write(ifileid,"(2i6)") iatm,0
    do ish=1,nshell
        if (shcen(ish)==iatm) then
            symbol=shelltype2name(shtype(ish))
            call struc2lc(symbol)
            write(ifileid,"(a,i4,' 1.0')") symbol,shcon(ish)
            if (ish==0) then
                istart=0
            else
                istart=sum(shcon(1:ish-1))
            end if
            do ipsh=istart+1,istart+shcon(ish)
                write(ifileid,"(2(1PE16.8))") primshexp(ipsh),primshcoeff(ipsh)
            end do
        end if
    end do
    write(ifileid,*)
end do
write(ifileid,*) 

if (any(shtype==-2).and.any(shtype==-3)) then !Default is 6d10f
    write(ifileid,"('[5D]')") !5d7f
else if (any(shtype==-2)) then
    write(ifileid,"('[5D10F]')") !5d10f
else if (any(shtype==-3)) then
    write(ifileid,"('[7F]')") !=6d7f
end if
if (any(shtype==-4)) write(ifileid,"('[9G]')") !Default is Cartesian G
if (any(shtype==-5)) write(ifileid,"('[11H]')") !Default is Cartesian H
write(ifileid,"(a)") "[MO]"
if (wfntype==0.or.wfntype==2.or.wfntype==3) then !Close shell
!     write(*,*) nmo
    do imo=1,nmo
        write(ifileid,"('Ene=',f16.8)") MOene(imo)
        write(ifileid,"('Spin= Alpha')")
        write(ifileid,"('Occup=',f10.6)") MOocc(imo)
        do ibas=1,nbasis
            write(ifileid,"(i6,f18.10)") ibas,CObasa(ibas,imo)
        end do
    end do
else !Open shell
    do isep=nmo,1,-1
        if (MOtype(isep)==1) exit
    end do
    do imo=1,isep !Alpha part
        write(ifileid,"('Ene=',f16.8)") MOene(imo)
        write(ifileid,"('Spin= Alpha')")
        write(ifileid,"('Occup=',f10.6)") MOocc(imo)
        do ibas=1,nbasis
            write(ifileid,"(i6,f18.10)") ibas,CObasa(ibas,imo)
        end do
    end do
    do imo=isep+1,nmo !Beta part
        write(ifileid,"('Ene=',f16.8)") MOene(imo)
        write(ifileid,"('Spin= Beta')")
        write(ifileid,"('Occup=',f10.6)") MOocc(imo)
        do ibas=1,nbasis
            write(ifileid,"(i6,f18.10)") ibas,CObasb(ibas,imo-isep)
        end do
    end do
end if

close(ifileid)
write(*,*) "Exporting Molden input file finished!"
end subroutine


!!!------------------------- Output current wavefunction to .fch file
subroutine outfch(outname,ifileid)
use defvar
use util
implicit real*8 (a-h,o-z)
character(len=*) outname
integer ifileid

open(ifileid,file=outname,status="replace")
write(ifileid,"(a)") "Generated by Multiwfn"
write(ifileid,"(a10,a30,a30)") "SP        ","RB3LYP                        ","                      6-31G(d)"
write(ifileid,"(A40,3X,A1,5X,I12)") "Number of atoms                         ","I",ncenter
write(ifileid,"(A40,3X,A1,5X,I12)") "Charge                                  ","I",sum(a(:)%charge)-nelec
write(ifileid,"(A40,3X,A1,5X,I12)") "Multiplicity                            ","I",nint(naelec-nbelec)+1
write(ifileid,"(A40,3X,A1,5X,I12)") "Number of electrons                     ","I",nint(nelec)
write(ifileid,"(A40,3X,A1,5X,I12)") "Number of alpha electrons               ","I",nint(naelec)
write(ifileid,"(A40,3X,A1,5X,I12)") "Number of beta electrons                ","I",nint(nbelec)
write(ifileid,"(A40,3X,A1,5X,I12)") "Number of basis functions               ","I",nbasis
write(ifileid,"(A40,3X,A1,5X,I12)") "Number of independent functions         ","I",nbasis
write(ifileid,"(A40,3X,A1,3X,'N=',I12)") "Atomic numbers                          ","I",ncenter
write(ifileid,"(6I12)") a(:)%index
write(ifileid,"(A40,3X,A1,3X,'N=',I12)") "Nuclear charges                         ","R",ncenter
write(ifileid,"(5(1PE16.8))") a(:)%charge
write(ifileid,"(A40,3X,A1,3X,'N=',I12)") "Current cartesian coordinates           ","R",ncenter*3
write(ifileid,"(5(1PE16.8))") (a(i)%x,a(i)%y,a(i)%z,i=1,ncenter)
!Basis function definition
write(ifileid,"(A40,3X,A1,5X,I12)") "Number of contracted shells             ","I",nshell
write(ifileid,"(A40,3X,A1,5X,I12)") "Number of primitive shells              ","I",nprimshell
write(ifileid,"(A40,3X,A1,5X,I12)") "Highest angular momentum                ","I",maxval(abs(shtype))
write(ifileid,"(A40,3X,A1,5X,I12)") "Largest degree of contraction           ","I",maxval(abs(shcon))
write(ifileid,"(A40,3X,A1,3X,'N=',I12)") "Shell types                             ","I",nshell
write(ifileid,"(6I12)") shtype(:)
write(ifileid,"(A40,3X,A1,3X,'N=',I12)") "Number of primitives per shell          ","I",nshell
write(ifileid,"(6I12)") shcon(:)
write(ifileid,"(A40,3X,A1,3X,'N=',I12)") "Shell to atom map                       ","I",nshell
write(ifileid,"(6I12)") shcen(:)
write(ifileid,"(A40,3X,A1,3X,'N=',I12)") "Primitive exponents                     ","R",nprimshell
write(ifileid,"(5(1PE16.8))") primshexp(:)
write(ifileid,"(A40,3X,A1,3X,'N=',I12)") "Contraction coefficients                ","R",nprimshell
write(ifileid,"(5(1PE16.8))") primshcoeff(:)
write(ifileid,"(A40,3X,A1,3X,'N=',I12)") "Coordinates of each shell               ","R",nshell*3
write(ifileid,"(5(1PE16.8))") (a(shcen(i))%x,a(shcen(i))%y,a(shcen(i))%z,i=1,nshell)
write(ifileid,"(A40,3X,A1,5X,1PE22.15)") "Virial Ratio                            ","R",virialratio
write(ifileid,"(A40,3X,A1,5X,1PE22.15)") "Total Energy                            ","R",totenergy
!Orbital informaiton
write(ifileid,"(A40,3X,A1,3X,'N=',I12)") "Alpha Orbital Energies                  ","R",nbasis
write(ifileid,"(5(1PE16.8))") MOene(1:nbasis)
if (wfntype==1.or.wfntype==4) then
    write(ifileid,"(A40,3X,A1,3X,'N=',I12)") "Beta Orbital Energies                   ","R",nbasis
    write(ifileid,"(5(1PE16.8))") MOene(nbasis+1:2*nbasis)
end if
write(ifileid,"(A40,3X,A1,3X,'N=',I12)") "Alpha MO coefficients                   ","R",nbasis*nbasis
write(ifileid,"(5(1PE16.8))") ((CObasa(ibasis,imo),ibasis=1,nbasis),imo=1,nbasis)
if (wfntype==1.or.wfntype==4) then
    write(ifileid,"(A40,3X,A1,3X,'N=',I12)") "Beta MO coefficients                    ","R",nbasis*nbasis
    write(ifileid,"(5(1PE16.8))") ((CObasb(ibasis,imo),ibasis=1,nbasis),imo=1,nbasis)
end if
write(ifileid,"(A40,3X,A1,3X,'N=',I12)") "Total SCF Density                       ","R",nbasis*(nbasis+1)/2
write(ifileid,"(5(1PE16.8))") ((Ptot(i,j),j=1,i),i=1,nbasis)
if (wfntype==1.or.wfntype==4) then
    write(ifileid,"(A40,3X,A1,3X,'N=',I12)") "Spin SCF Density                        ","R",nbasis*(nbasis+1)/2
write(ifileid,"(5(1PE16.8))") ((Palpha(i,j)-Pbeta(i,j),j=1,i),i=1,nbasis)
end if
close(ifileid)
write(*,*) "Exporting .fch file finished!"
end subroutine
