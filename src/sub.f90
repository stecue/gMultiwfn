!------------- Modify & Check wavefunction
subroutine modwfn
use defvar
use util
implicit real*8 (a-h,o-z)
character seltmpc*10,orbtype*5,selectyn,c1000tmp*1000
real*8 eigval(nbasis),eigvec(nbasis,nbasis),tmpmat(nbasis,nbasis)
integer orbarr(nmo)
! write(*,*) wfntype
do while(.true.)
    write(*,*) "          ============ Modify & Check wavefunction ============ "
    write(*,"(' Number of GTFs:',i6,', Orb:',i6,', Atoms:',i5,', A/B elec:',2f8.3)") nprims,nmo,ncenter,naelec,nbelec
    write(*,*) "-1 Return"
    write(*,*) "0 Save the modified wavefunction to a new .wfn file"
    write(*,*) "1 List all primitive function"
    if (ifiletype==1.or.ifiletype==9) write(*,*) "2 List all basis function" !appear only when input file contains basis information
    write(*,*) "3 List all orbitals"
    write(*,*) "4 Print detail information of an orbital"
    if (ifiletype==1.or.ifiletype==9) write(*,*) "5 Print coefficient matrix in basis function"
    if (ifiletype==1.or.ifiletype==9) write(*,*) "6 Print density matrix in basis function"
    if (ifiletype==1.or.ifiletype==9) write(*,*) "7 Print overlap matrix in basis function along with eigenvalues"
    write(*,*) "11 Swap some information of two primitive functions"
    write(*,*) "21 Set center of a primitive function"
    write(*,*) "22 Set type of a primitive function"
    write(*,*) "23 Set exponent of a primitive function"
    write(*,*) "24 Set coefficient of a primitive function"
    write(*,*) "25 Set coefficients of some GTFs satisfied certain conditions"
    write(*,*) "26 Set occupation number of some orbitals"
    write(*,*) "27 Set type of some orbitals"
    write(*,*) "31 Translate the system"
    write(*,*) "32 Translate and duplicate the system"
!     write(*,*) "33 Rotate wavefunction, namely X->Y, Y->Z, Z->X"
    if (imodwfn==0) write(*,*) "34 Set occupation number of inner orbitals to zero" !If occupation has been modified, don't do this to complicate things
    if (allocated(MOsym)) write(*,*) "35 Keep or discard orbital contributions according to irreducible rep."
!     write(*,*) "36 Invert phase of an orbital"
    read(*,*) iselect
    
    write(*,*)
    if (iselect==-1) then
        if ((ifiletype==1.or.ifiletype==9).and.imodwfn==1) then
            write(*,*) "Updating density matrix..."
            call genP
            write(*,*) "Density matrix has been updated"
        end if
        exit
    else if (iselect==0) then
        call outwfn("new.wfn",1,1,10)
        write(*,*) "New .wfn file has been outputted to new.wfn in current folder"
    else if (iselect==1) then
        do i=1,nprims
            write(*,"(i6,' Center:',i5,'(',a2,')','   Type: ',a,' Exponent:',D16.7)") i,b(i)%center,a(b(i)%center)%name,GTFtype2name(b(i)%functype),b(i)%exp
        end do
    else if (iselect==2) then
        do i=1,nbasis
            write(*,"(' Basis:',i5,' Shell:',i5,' Center:',i5,'(',a2,') Type:',a)")&
             i,basshell(i),bascen(i),a(bascen(i))%name,GTFtype2name(bastype(i))
        end do
    else if (iselect==3) then
        do i=1,nmo
            if (MOtype(i)==0) orbtype="A+B"
            if (MOtype(i)==1) orbtype="Alpha"
            if (MOtype(i)==2) orbtype="Beta"
            if (allocated(MOsym)) then
                write(*,"(' Orbital:',i5,' Energy(a.u.):',f12.6,' Occ:',f10.6,' Type: ',a,' Sym: ',a)") i,MOene(i),MOocc(i),orbtype,MOsym(i)
            else
                write(*,"(' Orbital:',i5,' Energy(a.u.):',f12.6,' Occ:',f10.6,' Type: ',a)") i,MOene(i),MOocc(i),orbtype
            end if
        end do
    else if (iselect==4) then
        write(*,*) "Input the orbital index, e.g. 12"
        read(*,*) i
        if (i<1.or.i>nmo) then
            write(*,"('Invalid orbital index, should within range of',i5,' and ',i5)") 1,nmo
        else
            write(*,"(' Occupation number is ',f12.7,'     Energy is',f12.6,' Hartree')") MOocc(i),MOene(i)
            if (MOtype(i)==0) write(*,*) "This is a close-shell orbital"
            if (MOtype(i)==1) write(*,*) "This is an alpha orbital"
            if (MOtype(i)==2) write(*,*) "This is a beta orbital"
            write(*,*)
            do j=1,nprims
                write(*,"(' GTF:',i6,' Cen:',i5,'(',a2,')',' Type: ',a,' Coeff:',1PD15.8,'  Exp: ',1PD13.7)") &
                j,b(j)%center,a(b(j)%center)%name,GTFtype2name(b(j)%functype),co(i,j),b(j)%exp
            end do
            write(*,"(a,/)") " Note: The ""coeff."" are expansion coefficients of orbitals with respect to GTFs, including normalization constant"
            if (allocated(b)) then
                do j=1,nbasis
                    if (MOtype(i)==0.or.MOtype(i)==1) covalue=CObasa(j,i)
                    if (MOtype(i)==2) covalue=CObasb(j,i-nbasis)
                    write(*,"(' Basis func:',i6,'  Cen:',i5,'(',a2,')',' Shell:',i5,' Type: ',a,' Coeff:',f12.8)") &
                    j,bascen(j),a(bascen(j))%name,basshell(j),GTFtype2name(bastype(j)),covalue
                end do
            end if
            write(*,"(a,/)") " Note: The ""coeff."" are expansion coefficients of orbitals with respect to basis functions, which are normalized functions"
        end if
    else if (iselect==5) then
        write(*,*) "(i,j) element means coefficient of ith basis function in jth orbital"
        if (wfntype==0.or.wfntype==2.or.wfntype==3) then
            call showmatgau(cobasa,"Coefficient matrix",0)
        else if (wfntype==1.or.wfntype==4) then
            call showmatgau(cobasa,"Alpha coefficient matrix",0)
            call showmatgau(cobasb,"Beta coefficient matrix",0)
        end if
    else if (iselect==6) then
        call showmatgau(Ptot,"Total density matrix",1)
        sum=0
        do i=1,nbasis
            sum=sum+Ptot(i,i)
        end do
        write(*,"(' The trace of the density matrix:',f12.6)") sum
        if (wfntype==1.or.wfntype==2.or.wfntype==4) then
            suma=0
            sumb=0
            do i=1,nbasis
                suma=suma+Palpha(i,i)
                sumb=sumb+Pbeta(i,i)
            end do
            call showmatgau(Palpha,"Alpha density matrix",1)
            write(*,*)
            call showmatgau(Pbeta,"Beta density matrix",1)
            write(*,*)
            write(*,"(' The trace of the alpha and beta density matrix:',2f12.6)") suma,sumb
        end if
    else if (iselect==7) then
        tmpmat=sbas
        call showmatgau(Sbas,"Overlap matrix",1)
        call diagsymat(tmpmat,eigvec,eigval,ierror)
        write(*,*)
        write(*,*) "Eigenvalues:"
        write(*,"(6f12.8)") eigval
    else if (iselect==11) then
        write(*,*) "Swap information of which two GTFs? Input their indices  e.g. 18,21"
        read(*,*) i,j
        write(*,*) "Swap which information for the two GTFs?"
        write(*,*) "1 Swap all propertie"
        write(*,*) "2 Swap center"
        write(*,*) "3 Swap function type"
        write(*,*) "4 Swap exponent"
        write(*,*) "5 Swap orbital expansion coefficient"
        read(*,*) iswapcontent
        if (iswapcontent==1) call swap(i,j,"all")
        if (iswapcontent==2) call swap(i,j,"cen")
        if (iswapcontent==3) call swap(i,j,"typ")
        if (iswapcontent==4) call swap(i,j,"exp")
        if (iswapcontent==5) call swap(i,j,"MO ")
        write(*,*) "Swapping finished!"
    else if (iselect==21) then
        write(*,*) "Input the index of primitive function"
        read(*,*) i
        write(*,*) "Input the center"
        read(*,*) j
        if (j<=ncenter.and.j>0) then
            b(i)%center=j
        else
            write(*,"('Error: The value should >0 and <=',i7)") ncenter
        end if
    else if (iselect==22) then
        write(*,*) "Input the index of primitive function"
        read(*,*) i
        write(*,*) "Input the type"
        write(*,*) "Valid input: S/X/Y/Z/XX/YY/ZZ/XY/XZ/YZ/XXX/YYY/ZZZ/XXY/XXZ/YYZ/XYY/XZZ/YZZ/XYZ"
        write(*,*) "ZZZZ/YZZZ/YYZZ/YYYZ/YYYY/XZZZ/XYZZ/XYYZ/XYYY/XXZZ/XXYZ/XXYY/XXXZ/XXXY/XXXX"
        write(*,*) "ZZZZZ/YZZZZ/YYZZZ/YYYZZ/YYYYZ/YYYYY/XZZZZ/XYZZZ/XYYZZ/XYYYZ/XYYYY/XXZZZ/XXYZZ/XXYYZ/XXYYY/XXXZZ/XXXYZ/XXXYY/XXXXZ/XXXXY/XXXXX"
        read(*,*) seltmpc
        do j=1,size(GTFtype2name)
            if (seltmpc==GTFtype2name(j)) then
                b(i)%functype=j
                exit
            end if
            if (j==20) write(*,*) "Error: Couldn't recognize this type"
        end do
    else if (iselect==23) then
        write(*,*) "Input the index of primitive function"
        read(*,*) i
        write(*,*) "Input the exponent"
        read(*,*) rexp
        b(i)%exp=rexp
    else if (iselect==24) then
        write(*,*) "Input the index of primitive function"
        read(*,*) iprm
        write(*,*) "Input the orbital index, e.g. 12"
        read(*,*) imonum
        if (iprm<=nprims.and.iprm>0.and.imonum<=nmo.and.imonum>0) then
            write(*,*) "Input the coefficient"
            read(*,*) rcoeff
            CO(imonum,iprm)=rcoeff
        else
            write(*,"('Error: The index of function or orbital exceed valid range')")
        end if
    else if (iselect==25) then
        write(*,*) "The following your inputs are condition for filtering GTFs"
        write(*,*) "Rule of range input: 3,17 means from 3 to 17, 6,6 means only 6, 0,0 means all"
        write(*,*)
        write(*,*) "Input the range of index of GTFs"
        read(*,*) ind1,ind2
        write(*,*) "Input the range of atoms"
        read(*,*) iatm1,iatm2
        write(*,*) "Input the type of GTFs (one of S,X,Y,Z,XX,XY... ALL means all types)"
        read(*,*) seltmpc
        write(*,*) "Input the range of orbitals"
        read(*,*) imo1,imo2
        write(*,*) "Input the expansion coefficient you want to set"
        read(*,*) coval
        do itypfil=1,size(type2ix)
            if (seltmpc==GTFtype2name(itypfil)) exit
        end do
        if (ind1==0) ind1=1
        if (ind2==0) ind2=nprims
        if (iatm1==0) iatm1=1
        if (iatm2==0) iatm2=ncenter
        if (imo1==0) imo1=1
        if (imo2==0) imo2=nmo
        do i=ind1,ind2
            if (seltmpc=="ALL".or.seltmpc=="all") then
                if (b(i)%center>=iatm1.and.b(i)%center<=iatm2) CO(imo1:imo2,i)=coval
            else
                if (b(i)%center>=iatm1.and.b(i)%center<=iatm2.and.b(i)%functype==itypfil) CO(imo1:imo2,i)=coval
            end if
        end do
        write(*,*) "Done!"
    else if (iselect==26) then
        do while(.true.)
            write(*,*) "Select the orbitals for which the occupation numbers are needed to be changed"
            write(*,*) "e.g. 2,4,13-16,20 means select orbital 2,4,13,14,15,16,20"
            write(*,*) "Input 0 can select all orbitals, input q or 00 can return"
            read(*,"(a)") c1000tmp
            if (c1000tmp(1:1)=='q'.or.c1000tmp(1:2)=='00') exit
            if (c1000tmp(1:1)=='0') then
                numorbsel=nmo
                do i=1,nmo
                    orbarr(i)=i
                end do
            else
                call str2arr(c1000tmp,numorbsel,orbarr)
                if ( any(orbarr(1:numorbsel)<1).or.any(orbarr(1:numorbsel)>nmo) ) then
                    write(*,*) "Error: One or more orbital indices exceeded valid range!"
                    cycle
                end if
            end if
            write(*,*) "Set to which value? e.g. 1.2"
            write(*,*) "Note:"
            write(*,"(a)") " You can also input for example ""+1.1"" ""-1.1"" ""*1.1"" ""/1.1"" to add, minus, multiply and divide the occupation numbers by 1.1"
            write(*,"(a)") " To recover the initial occupation numbers, input ""i"""
            write(*,"(a)") " To generate occupation state for calculating odd electron density, input ""odd"""
            read(*,"(a)") c1000tmp
            if (index(c1000tmp,"odd")/=0) then
                do iorb=1,numorbsel
                    MOocc(orbarr(iorb))=min(2-MOocc(orbarr(iorb)),MOocc(orbarr(iorb)))
                end do
                write(*,*) "Done!"
            else if (index(c1000tmp,"i")/=0) then
                MOocc(orbarr(1:numorbsel))=MOocc_org(orbarr(1:numorbsel))
                write(*,*) "The occupation numbers have been recovered"
            else if (c1000tmp(1:1)=='+'.or.c1000tmp(1:1)=='-'.or.c1000tmp(1:1)=='*'.or.c1000tmp(1:1)=='/') then
                read(c1000tmp(2:),*) tmpval
                if (c1000tmp(1:1)=='+') MOocc(orbarr(1:numorbsel))=MOocc(orbarr(1:numorbsel))+tmpval
                if (c1000tmp(1:1)=='-') MOocc(orbarr(1:numorbsel))=MOocc(orbarr(1:numorbsel))-tmpval
                if (c1000tmp(1:1)=='*') MOocc(orbarr(1:numorbsel))=MOocc(orbarr(1:numorbsel))*tmpval
                if (c1000tmp(1:1)=='/') MOocc(orbarr(1:numorbsel))=MOocc(orbarr(1:numorbsel))/tmpval
            else
                read(c1000tmp,*) tmpval
                MOocc(orbarr(1:numorbsel))=tmpval
                write(*,*) "Done!"
            end if
            call updatenelec !Update the number of electrons
            if (occtmp==0D0) write(*,"(/,a)") " Note: When saving present wavefunction to new .wfn file, the orbitals with zero occupation number will be discarded"
            imodwfn=1
            if (any(MOocc/=int(MOocc))) then
                if (wfntype==0) then
                    wfntype=3 !RHF-> Restricted post-HF wavefunction
                    write(*,*) "Note: Now the wavefunction is recognized as restricted post-HF wavefunction"
                else if (wfntype==1.or.wfntype==2) then !UHF/ROHF-> Unrestricted post-HF wavefunction
                    wfntype=4
                    write(*,*) "Note: Now the wavefunction is recognized as unrestricted post-HF wavefunction"
                end if
            end if
            write(*,*)
        end do
    else if (iselect==27) then
        write(*,*) "Set type for which range of orbitals? e.g. 14,17"
        write(*,*) "Hint: Input 0,0 can select all orbitals"
        read(*,*) iorb1,iorb2
        if (iorb1==0.and.iorb2==0) then
            iorb1=1
            iorb2=nmo
        end if
        write(*,*) "Set to which type? 0=Alpha&Beta 1=Alpha 2=Beta"
        read(*,*) isettype
        MOtype(iorb1:iorb2)=isettype
        !Recount alpha and beta electrons
        call updatenelec
        write(*,*) "Done!"
        !Update wavefunction type
        if (all(MOtype==0)) then
            if (all(MOocc==nint(MOocc))) then !All A+B orbital & integer occupation
                wfntype=0
                write(*,"('Note: Now the wavefunction is recognized as restricted close-shell single-determinant wavefunction')")
            else !All A+B orbital & partial occupation
                wfntype=3
                write(*,"('Note: Now the wavefunction is recognized as close-shell post-HF wavefunction')")
            end if
            write(*,*)
        else
            if (any(MOocc/=nint(MOocc)).and.all(MOtype/=0)) then !Either A or B, and partial occupation
                wfntype=4
                write(*,"('Note: Now the wavefunction is recognized as open-shell post-HF wavefunction')")
            else if (all(MOocc==nint(MOocc)).and.all(MOtype/=0)) then !Integer occupation and either A or B
                wfntype=1
                write(*,"('Note: Now the wavefunction is recognized as unrestricted single-determinant wavefunction')")
            else if (all(MOocc==nint(MOocc)).and.any(MOtype==0).and.all(MOtype/=2)) then !Integer occupation and at least one orbital is A, and B is unexisted
                wfntype=2
                write(*,"('Note: Now the wavefunction is recognized as restricted open-shell wavefunction')")
            else
                write(*,"('Warning: The type of present wavefunction cannot be identified! You need to reset orbital types')")
            end if
            write(*,*)
        end if
        imodwfn=1
    else if (iselect==31) then
        write(*,*) "Input X,Y,Z of translation vector (e.g. 3.2,1.0,0)"
        read(*,*) pbctransx,pbctransy,pbctransz
        write(*,*) "You inputted coordinates are in which unit?  1:Bohr  2:Angstrom"
        read(*,*) iunit
        if (iunit==2) then
            pbctransx=pbctransx/b2a
            pbctransy=pbctransy/b2a
            pbctransz=pbctransz/b2a
        end if
        do i=1,ncenter
            a(i)%x=a(i)%x+pbctransx
            a(i)%y=a(i)%y+pbctransy
            a(i)%z=a(i)%z+pbctransz
        end do
        imodwfn=1
    else if (iselect==32) then
        write(*,*) "Input X,Y,Z of translation vector (e.g. 3.2,1.0,0)"
        read(*,*) pbctransx,pbctransy,pbctransz
        write(*,*) "You inputted coordinates are in which unit?  1:Bohr  2:Angstrom"
        read(*,*) iunit
        if (iunit==2) then
            pbctransx=pbctransx/b2a
            pbctransy=pbctransy/b2a
            pbctransz=pbctransz/b2a
        end if
        write(*,*) "Duplicate system how many times?"
        read(*,*) numdup
        !_tmp is for backing up current information
        allocate(a_tmp(ncenter))
        allocate(b_tmp(nprims))
        allocate(CO_tmp(nmo,nprims))
        a_tmp=a
        b_tmp=b
        CO_tmp=CO
        deallocate(a,b,CO)
        nprims_tmp=nprims
        ncenter_tmp=ncenter
        nprims=nprims*(numdup+1)
        ncenter=ncenter*(numdup+1)
        nelec=nelec*(numdup+1)
        naelec=naelec*(numdup+1)
        nbelec=nbelec*(numdup+1)
        allocate(a(ncenter))
        allocate(b(nprims))
        allocate(CO(nmo,nprims))
        do idup=0,numdup
            a(ncenter_tmp*idup+1:ncenter_tmp*(idup+1))=a_tmp(1:center_tmp)
            a(ncenter_tmp*idup+1:ncenter_tmp*(idup+1))%x=a_tmp(1:center_tmp)%x+pbctransx*idup
            a(ncenter_tmp*idup+1:ncenter_tmp*(idup+1))%y=a_tmp(1:center_tmp)%y+pbctransy*idup
            a(ncenter_tmp*idup+1:ncenter_tmp*(idup+1))%z=a_tmp(1:center_tmp)%z+pbctransz*idup
            b(nprims_tmp*idup+1:nprims_tmp*(idup+1))=b_tmp(1:nprims_tmp)
            b(nprims_tmp*idup+1:nprims_tmp*(idup+1))%center=b_tmp(1:nprims_tmp)%center+ncenter_tmp*idup
            CO(:,nprims_tmp*idup+1:nprims_tmp*(idup+1))=CO_tmp(:,1:nprims_tmp) !Notice that the orbitals do not satisify normalization condition any more, and the orbital occupation number will be artifical
        end do
        deallocate(a_tmp,b_tmp,CO_tmp)
        imodwfn=1
        call gendistmat !The number of atoms have changed, so we must update distance matrix
    else if (iselect==33) then
        write(*,*) "Rotate which orbital? (Input 0 to rotate all orbitals)"
        read(*,*) iorb
        if (iorb/=0) then
            call orbcoeffrotate(iorb)
        else if (iorb==0) then
            do imo=1,nmo
                call orbcoeffrotate(imo)
            end do
            write(*,*) "Also rotate atomic coordinates? (y/n)"
            read(*,*) selectyn
            if (selectyn=='y'.or.selectyn=='Y') then
                do iatm=1,ncenter
                    tmpval=a(iatm)%x
                    a(iatm)%x=a(iatm)%z
                    a(iatm)%z=a(iatm)%y
                    a(iatm)%y=tmpval
                end do
            end if
        end if
        write(*,*) "Done!"
    else if (iselect==34) then
        innerel=0
        do i=1,ncenter
            if (int(a(i)%charge)/=a(i)%index) cycle !For the atom used ECP, skip it
            if (a(i)%index>2.and.a(i)%index<=10) innerel=innerel+2
            if (a(i)%index>10.and.a(i)%index<=18) innerel=innerel+10
            if (a(i)%index>18.and.a(i)%index<=36) innerel=innerel+18
            if (a(i)%index>36.and.a(i)%index<=54) innerel=innerel+36
            if (a(i)%index>54.and.a(i)%index<=86) innerel=innerel+54
            if (a(i)%index>86) innerel=innerel+86
        end do
        nelec=nelec-innerel
        naelec=naelec-innerel/2
        nbelec=nbelec-innerel/2
        if (wfntype==1.or.wfntype==4) then !UHF and U-post-HF wfn
            MOocc(1:innerel/2)=0D0
            do j=1,nmo !Where the first beta orbital appear now
                if (motype(j)==2) exit
            end do
            MOocc(1:j+innerel/2-1)=0D0
            write(*,"('The effect of ',i7,' lowest energy orbitals have been discarded')") innerel
        else if (wfntype==0.or.wfntype==2.or.wfntype==3) then !restricted(=0) or RO(=2) or post-R(=3) wavefunction
            MOocc(1:innerel/2)=0D0
            write(*,"('The effect of ',i7,' lowest energy orbitals have been discarded')") innerel/2
        end if
        if (wfntype==3.or.wfntype==4) write(*,"('Warning: Discarding inner orbitals for post-HF wavefunction will lead to unexpected result!')") 
        imodwfn=1
    else if (iselect==35) then
        call SelMO_IRREP
    else if (iselect==36) then
        write(*,*) "Input lower and upper limits of orbital index, e.g. 8,15"
        read(*,*) ilow,ihigh
        do idx=ilow,ihigh
            CO(idx,:)=-CO(idx,:)
            if (allocated(CObasa)) then
                if (idx<=nbasis) then
                    CObasa(:,idx)=-CObasa(:,idx)
                else
                    CObasb(:,idx-nbasis)=-CObasb(:,idx-nbasis)
                end if
            end if
        end do
        write(*,*) "Done!"
        imodwfn=1
    end if
    write(*,*)
end do
end subroutine


!!---------------- Select MOs according to irreducible representation
subroutine SelMO_IRREP
use defvar
use util
character symlab(nmo)*4,c200tmp*200,symstat(nmo)*9 !Allocate the array lengths as upper limit
integer tmparr(nmo),symNorb(nmo) !Allocate the array lengths as upper limit
if (wfntype/=0.and.wfntype/=1) then
    write(*,"(a)") " Error: This function only works for RHF or UHF wavefunctions (or the DFT counterparts)"
    return
end if
nsym=0
do imo=1,nmo
    if (MOocc_org(imo)==0D0) cycle
    if (all(symlab(1:nsym)/=MOsym(imo))) then
        nsym=nsym+1
        symlab(nsym)=MOsym(imo)
    end if
end do
symNorb=0
do imo=1,nmo
    if (MOocc_org(imo)==0D0) cycle
    do isym=1,nsym
        if (MOsym(imo)==symlab(isym)) symNorb(isym)=symNorb(isym)+1
    end do
end do
symstat="Normal"
MOocc=MOocc_org
if (imodwfn==1) write(*,*) "Note: Original occupation status has been recovered"
write(*,*) "Note: Only the orbitals that originally occupied are taken into account here"
do while(.true.)
    write(*,*)
    write(*,*) "Information of various irreducible representations:"
    do isym=1,nsym
        write(*,"(i5,'  Sym: ',a,'  N_orb:',i5,'    Status: ',a)") isym,symlab(isym),symNorb(isym),symstat(isym)
    end do
    write(*,*)
    write(*,*) "0 Save and return"
    write(*,*) "1 Discard specific irreducible representations"
    write(*,*) "2 Recover original status"
    write(*,*) "3 Reverse status"
    read(*,*) isel
    
    if (isel==0) then
        call updatenelec
        imodwfn=1
        write(*,*) "The current orbital occupation status has been saved"
        write(*,*) "Updating density matrix..."
        call genP
        write(*,*) "Density matrix has been updated"
        exit
    else if (isel==2) then
        MOocc=MOocc_org
        symstat="Normal"
    else if (isel==1.or.isel==3) then
        if (isel==1) then
            write(*,*) "Input the index of the irreducible representations to be discarded, e.g. 1,3-5"
            read(*,"(a)") c200tmp
            call str2arr(c200tmp,nsymsel,tmparr)
            do isym=1,nsymsel
                symstat(tmparr(isym))="Discarded"
            end do
        else if (isel==3) then
            do isym=1,nsym
                if (symstat(isym)=="Normal") then
                    symstat(isym)="Discarded"
                else
                    symstat(isym)="Normal"
                end if
            end do
        end if
        do imo=1,nmo
            if (MOocc_org(imo)==0D0) cycle
            do isym=1,nsym
                if (MOsym(imo)==symlab(isym)) then
                    if (symstat(isym)=="Normal") MOocc(imo)=MOocc_org(imo)
                    if (symstat(isym)=="Discarded") MOocc(imo)=0D0
                    exit
                end if
            end do
        end do
    end if
    write(*,*) "Done!"
end do
end subroutine

!!---------- Update the number of electrons
subroutine updatenelec
use defvar
integer imo
nelec=0
naelec=0
nbelec=0
do imo=1,nmo
    if (MOtype(imo)==0) then
        naelec=naelec+MOocc(imo)/2D0
        nbelec=nbelec+MOocc(imo)/2D0
    else if (MOtype(imo)==1) then
        naelec=naelec+MOocc(imo)
    else if (MOtype(imo)==2) then
        nbelec=nbelec+MOocc(imo)
    end if
end do
nelec=naelec+nbelec
end subroutine
            

!!!-------- Check if present wavefunction is sanity, i.e. all orbital satisfies normalization condition
subroutine wfnsanity
use defvar
implicit real*8 (a-h,o-z)
real*8 GTFSmat(nprims*(nprims+1)/2)
call genGTFSmat(GTFSmat,nprims*(nprims+1)/2)
rmaxdev=0
do imo=1,nmo
    tmp=0
    do iGTF=1,nprims
        do jGTF=iGTF+1,nprims
            tmp=tmp+2*co(imo,iGTF)*co(imo,jGTF)*GTFSmat(jGTF*(jGTF-1)/2+iGTF)
        end do
        tmp=tmp+co(imo,iGTF)**2*GTFSmat(iGTF*(iGTF-1)/2+iGTF)
    end do
    write(*,"(' Orbital',i7,', Occ:',f8.4,'   Value:',f16.10)") imo,MOocc(imo),tmp
    if (tmp>rmaxdev) rmaxdev=tmp
end do
write(*,"(' Maximum deviation to 1:',f16.10)") rmaxdev-1
end subroutine



!!---------- Return normalization coefficient for specific type of cartesian type GTF, see Levine 5ed p487
!The meaning of itype is defined in GTFtype2name
real*8 function normgau(itype,exp)
use defvar
use util
implicit real*8 (a-h,o-z)
ix=type2ix(itype)
iy=type2iy(itype)
iz=type2iz(itype)
normgau=(2*exp/pi)**0.75D0*dsqrt( (8*exp)**(ix+iy+iz)*ft(ix)*ft(iy)*ft(iz)/(ft(2*ix)*ft(2*iy)*ft(2*iz)) )
end function
!!!---------- Return normalization coefficient for specific type of spherical harmonic GTF
!See http://en.wikipedia.org/wiki/Gaussian_orbital for the formula
!!!This is useless for resolving the contraction coefficent problem of Molden input file of ORCA , I don't know why
!Lval is the angular moment, 0/1/2/3/4=s/p/d/f/g
! real*8 function rnormgau_sph(Lval,exp)
! use defvar
! use util
! integer Lval
! real*8 exp
! rnormgau_sph=exp**(Lval/2D0+0.75D0) / (dsqrt(pi)*2**(-0.25D0-Lval/2D0)) / dsqrt(gamma_ps(Lval+1))
! end function

subroutine genn1n2nf(Lval,n1,n2,nf)
integer Lval,n1,n2,nf
if (Lval==0) then
    n1=3
    n2=3
    nf=1
else if (Lval==1) then
    n1=7
    n2=5
    nf=1
else if (Lval==2) then
    n1=11
    n2=7
    nf=9
else if (Lval==3) then
    n1=15
    n2=9
    nf=225
else if (Lval==4) then
    n1=19
    n2=11
    nf=11025
end if
end subroutine
!!---- Used to produce normalization factor to counteract the contraction coefficient problem of the Molden input file generated by ORCA
!I don't know where the formula comes from
!Lval=0/1/2/3/4 corresponds to s/p/d/f/g
real*8 function rnormgau_ORCA(exp,Lval)
real*8 exp,pi
integer Lval,n1,n2,nf
pi=acos(-1D0)
call genn1n2nf(Lval,n1,n2,nf)
rnormgau_ORCA=dsqrt(dsqrt(2**n1*exp**n2/(pi**3*nf*nf)))
end function
!----- Renormalization of Gauss basis functions for Molden input file
subroutine renormmoldengau(nlen,Lval,exp,con)
implicit real*8 (a-h,o-z)
integer nlen,Lval
real*8 exp(nlen),con(nlen),ctmp(nlen)
pi=acos(-1D0)
call genn1n2nf(Lval,n1,n2,nf)
fc=(2D0**n1)/(pi**3*nf)
do i=1,nlen
    prmnormfac=sqrt(sqrt(fc*(exp(i)**n2)))
    ctmp(i)=con(i)*prmnormfac
end do
facnorm=0D0
do i=1,nlen
    do j=1,i
      expavg=(exp(i)+exp(j))/2D0
      facadd=ctmp(i)*ctmp(j)/sqrt(fc*(expavg**n2))
      if (i/=j) facadd=facadd*2
      facnorm=facnorm+facadd
    end do
end do
if (facnorm>1D-10) facnorm=1/sqrt(facnorm)
con=con*facnorm
end subroutine


!!----- Use Lowdin orthogonalization method to transform density matrix (Xmat) and coefficient to ortho-normalized basis, meanwhile update Sbas to identity matrix
!see Szabo p143
subroutine symmortho
use defvar
use util
implicit real*8 (a-h,o-z)
real*8 Umat(nbasis,nbasis),svalvec(nbasis),Xmat(nbasis,nbasis) !workvec(3*nbasis-1)
! call DSYEV('V','U',nbasis,sbas,nbasis,svalvec,workvec,3*nbasis-1,ierror) !lapack, the resultant sbas is eigenvector matrix
! call diagmat(Sbas,Umat,svalvec) !My slow diagonalization routine
call diagsymat(Sbas,Umat,svalvec,ierror)
if (ierror/=0) write(*,*) "Error: Diagonalization of overlap matrix failed!"
!Now Sbas is diagonalized matrix
forall (i=1:nbasis) Sbas(i,i)=dsqrt(svalvec(i))
Xmat=matmul(matmul(Umat,Sbas),transpose(Umat)) !Then Xmat is S^0.5
Ptot=matmul(matmul(Xmat,Ptot),Xmat)
if (allocated(Palpha)) then
    Palpha=matmul(matmul(Xmat,Palpha),Xmat)
    Pbeta=Ptot-Palpha
end if
Cobasa=matmul(Xmat,Cobasa)
if (allocated(Cobasb)) Cobasb=matmul(Xmat,Cobasb)
forall(i=1:nbasis) Sbas(i,i)=1D0 !Reconstruct overlap matrix in ortho basis function
end subroutine
!!----- Like symmortho, but input overlap matrix, and only output Lowdin orthogonalization transformation matrix Xmat
! innbas=input number of basis. Smatin is input overlap matrix, which will be modified
! inv=1 get S^-0.5, other get S^0.5
subroutine symmorthomat(innbas,Smatin,Xmat,inv)
use defvar
use util
integer inv
real*8 Umat(innbas,innbas),svalvec(innbas),Smatin(nbasis,nbasis),Smat(innbas,innbas),Xmat(innbas,innbas)
Smat=Smatin
call diagsymat(Smat,Umat,svalvec,ierror)
forall (i=1:nbasis) Smat(i,i)=dsqrt(svalvec(i))
if (inv==1) forall (i=1:nbasis) Smat(i,i)=1D0/Smat(i,i)
if (ierror/=0) write(*,*) "Error: Diagonalization of overlap matrix is fail!"
Xmat=matmul(matmul(Umat,Smat),transpose(Umat))
end subroutine


!!!------------------------- Generate distance matrix
subroutine gendistmat
use defvar
implicit real*8 (a-h,o-z)
if (allocated(distmat)) deallocate(distmat)
allocate(distmat(ncenter,ncenter))
distmat=0.0D0
do i=1,ncenter
    do j=i+1,ncenter
        distmat(i,j)=dsqrt((a(i)%x-a(j)%x)**2+(a(i)%y-a(j)%y)**2+(a(i)%z-a(j)%z)**2)
    end do
end do
distmat=distmat+transpose(distmat)
end subroutine


!!!------------------------- Swap two gaussian primitive functions
subroutine swap(i,j,swaptype)
use defvar
integer n,i,j
character*3 swaptype
type(primtype) tempb !For exchanging basis functions' order
if (swaptype=="all") then
    tempb=b(i)
    b(i)=b(j)
    b(j)=tempb
else if (swaptype=="cen") then
    tempb%center=b(i)%center
    b(i)%center=b(j)%center
    b(j)%center=tempb%center
else if (swaptype=="typ") then
    tempb%functype=b(i)%functype
    b(i)%functype=b(j)%functype
    b(j)%functype=tempb%functype
else if (swaptype=="exp") then
    tempb%exp=b(i)%exp
    b(i)%exp=b(j)%exp
    b(j)%exp=tempb%exp
end if
if (swaptype=="all".or.swaptype=="MO ") then
    do n=1,nmo
        temp=co(n,i)
        co(n,i)=co(n,j)
        co(n,j)=temp
    end do
end if
end subroutine


!!!---- Rotate(exchange) GTF and basis function coefficients within all shell in different direction of specific orbital
! use this three times, namely XYZ->ZXY->YZX->XYZ, the coefficient recovered.
! In detail, for examples, for d-type will lead to such coefficient exchange: XX to YY, YY to ZZ, ZZ to XX, XY to YZ, XZ to XY, YZ to XZ
! exchange only involve the GTFs/basis func. in the same shell
subroutine orbcoeffrotate(orb) !orb=Rotate which orbital
use defvar
implicit real*8 (a-h,o-z)
integer orb
real*8 COorborg(nprims) !For backing up origin CO
real*8 CObasa_tmp(nbasis) !For backing up origin CObasa
real*8 CObasb_tmp(nbasis) !For backing up origin CObasb
COorborg(:)=CO(orb,:)
do i=1,nprims
    ixtmp=type2iz(b(i)%functype)
    iytmp=type2ix(b(i)%functype)
    iztmp=type2iy(b(i)%functype)
    do j=1,nprims
        if (type2ix(b(j)%functype)==ixtmp.and.type2iy(b(j)%functype)==iytmp.and.&
        type2iz(b(j)%functype)==iztmp.and.b(j)%exp==b(i)%exp.and.b(j)%center==b(i)%center) CO(orb,j)=COorborg(i)
    end do
end do
if (allocated(CObasa)) then
    CObasa_tmp=CObasa(:,orb)
    if (wfntype==1.or.wfntype==4) CObasb_tmp=CObasb(:,orb)
    do iatm=1,ncenter
        do ibas=basstart(iatm),basend(iatm)
            ityp=bastype(ibas)
            ixtmp=type2iz(ityp)
            iytmp=type2ix(ityp)
            iztmp=type2iy(ityp)
            do jbas=basstart(iatm),basend(iatm)
                jtyp=bastype(jbas)
                if (type2ix(jtyp)==ixtmp.and.type2iy(jtyp)==iytmp.and.type2iz(jtyp)==iztmp.and.&
                basshell(ibas)==basshell(jbas)) then
                    CObasa(jbas,orb)=CObasa_tmp(ibas)
                    if (wfntype==1.or.wfntype==4) CObasb(jbas,orb)=CObasb_tmp(ibas)
                end if
            end do
        end do
    end do
end if
end subroutine


!!!------ Define property/origin/spacing/grid number and then save to a 3D matrix, infomode=1 means silent
!! iorb is used to choose the orbital for whose wavefunction will be calculated. This can be an arbitrary value if functype/=4
subroutine savecubmat(functype,infomode,iorb)
use defvar
use util
use function
implicit none
integer walltime1,walltime2
integer :: i,j,k,ii,infomode,functype,calcfunc,ifinish,iorb !Calculate which orbital wavefunction for fmo routine
real*8 t,time_begin,time_end,time_endtmp,tmpx,tmpy,tmpz,xarr(nx),yarr(ny),zarr(nz)
character c80tmp*80
iorbsel=iorb
calcfunc=functype
if (functype==12) calcfunc=8 !If calculate total ESP, first calculate nuclear ESP, then invoke espcub
if (infomode==0.and.functype/=12) then
    if (expcutoff<0) write(*,"(' Note: All exponential functions exp(x) with x<',f8.3,' will be ignored ')") expcutoff
end if

ii=10
ifinish=0
!Cut the minimal noise at the end of the coordinate, otherwise the originally symmetry points may become unsymmetry
do k=1,nz
    write(c80tmp,"(D20.13)") orgz+(k-1)*dz
    read(c80tmp,*) zarr(k)
end do
do j=1,ny
    write(c80tmp,"(D20.13)") orgy+(j-1)*dy
    read(c80tmp,*) yarr(j)
end do
do i=1,nx
    write(c80tmp,"(D20.13)") orgx+(i-1)*dx
    read(c80tmp,*) xarr(i)
end do

call walltime(walltime1)
CALL CPU_TIME(time_begin)
!$OMP PARALLEL DO SHARED(cubmat,ifinish) PRIVATE(i,j,k,tmpx,tmpy,tmpz) schedule(dynamic) NUM_THREADS( nthreads  )
do k=1,nz
    tmpz=zarr(k)
    do j=1,ny
        tmpy=yarr(j)
        do i=1,nx
            tmpx=xarr(i)
            if (calcfunc==1513) then !Only involved by funcvsfunc routine, when RDG and signlambda2rho is combined
                call signlambda2rho_RDG(tmpx,tmpy,tmpz,cubmat(i,j,k),cubmattmp(i,j,k))
            else if (calcfunc==1614) then !Only involved by funcvsfunc routine, when RDG and signlambda2rho is combined
                call signlambda2rho_RDG_prodens(tmpx,tmpy,tmpz,cubmat(i,j,k),cubmattmp(i,j,k))
            else
                cubmat(i,j,k)=calcfuncall(calcfunc,tmpx,tmpy,tmpz)
            end if
        end do
    end do
    if (infomode==0.and.functype/=12) then
        if ( nthreads  ==1) then
            CALL CPU_TIME(time_endtmp)
            t=dfloat(k)/nz*100 !Completed percent
            if (k==1) then !Show approximate time at start
                t=(time_endtmp-time_begin)*(nz-1) !How many time to work out remain work
                write(*,"(' Calculation will take up CPU time about',f9.2,' seconds (',f8.2,' minutes)')") t,t/60D0
            else if (t>=ii) then
                ii=ii+10
                write(*,"(f5.1,'% completed,',f8.1,' seconds remain')") t,(time_endtmp-time_begin)/k*(nz-k)
            end if
        else
            ifinish=ifinish+1
            write(*,"(' Finished:',i5,'  /',i5)") ifinish,nz
        end if
    end if
end do
!$OMP END PARALLEL DO
CALL CPU_TIME(time_end)
call walltime(walltime2)
if (infomode==0.and.functype/=12.and.isys==1) write(*,"(' Calculation took up CPU time',f12.2,'s, wall clock time',i10,'s')") time_end-time_begin,walltime2-walltime1
if (functype==12) call espcub
end subroutine


!!!------ Output molecular formula
subroutine showformula
use defvar
implicit real*8 (a-h,o-z)
character*6 tmp
write(*,"(' Formula: ')",advance="no")
do i=0,nelesupp
    n=0
    do iatm=1,ncenter
        if (a(iatm)%index==i) n=n+1
    end do
    write(tmp,"(i6)") n
    if (n/=0) write(*,"(a,a,' ')",advance="no") trim(ind2name(i)),trim(adjustl(tmp))
end do
write(*,*)
end subroutine


!!!----------- Resize number of orbitals of CO, MOene, MOtype, MOocc to "newnmo", also resize number of GTFs of CO to "newnprims"
subroutine resizebynmo(newnmo,newnprims)
use defvar
implicit real*8 (a-h,o-z)
real*8,allocatable :: CO_bk(:,:),MOene_bk(:),MOocc_bk(:)
integer,allocatable :: MOtype_bk(:)
integer newnmo,oldnmo,newnprims,oldnprims
oldnmo=size(CO,1)
oldnprims=size(CO,2)
allocate(CO_bk(oldnmo,oldnprims),MOene_bk(oldnmo),MOocc_bk(oldnmo),MOtype_bk(oldnmo))
CO_bk=CO
MOene_bk=MOene
MOocc_bk=MOocc
MOtype_bk=MOtype
deallocate(CO,MOene,MOocc,MOtype)
allocate(CO(newnmo,newnprims),MOene(newnmo),MOocc(newnmo),MOtype(newnmo))
if (newnmo>=oldnmo) then !Enlarge array size, don't forget to fill the gap afterwards
    if (newnprims>=oldnprims) CO(1:oldnmo,1:oldnprims)=CO_bk(:,:)
    if (newnprims<oldnprims) CO(1:oldnmo,:)=CO_bk(:,1:newnprims)
    MOene(1:oldnmo)=MOene_bk(:)
    MOocc(1:oldnmo)=MOocc_bk(:)
    MOtype(1:oldnmo)=MOtype_bk(:)
else if (newnmo<oldnmo) then !Reduce array size
    if (newnprims>=oldnprims) CO(:,1:oldnprims)=CO_bk(1:newnmo,:)
    if (newnprims<oldnprims) CO(:,:)=CO_bk(1:newnmo,1:newnprims)
    MOene(:)=MOene_bk(1:newnmo)
    MOocc(:)=MOocc_bk(1:newnmo)
    MOtype(:)=MOtype_bk(1:newnmo)
end if
deallocate(CO_bk,MOene_bk,MOocc_bk,MOtype_bk)
end subroutine


!!!------------ Generate gjf of atoms in molecule, and invoke Gaussian to get .wfn, then input them into custom list
subroutine setpromol
use defvar
use util
implicit real*8 (a-h,o-z)
integer :: i,j,k,itype=0
character*2 typename(100),nametmp
character*80 basisset,tmpdir,c80tmp
character*80 outwfnname
logical alivegauout,alivewfntmp,aliveatomwfn

!The only difference between c80tmp and tmpdir is that the latter has \ or / separator at the end
if (iwfntmptype==1) then
    if (isys==1) tmpdir="wfntmp\"
    if (isys==2.or.isys==3) tmpdir="wfntmp/"
    c80tmp="wfntmp"
    inquire(file='./wfntmp/.',exist=alivewfntmp)
    if ((isys==1).and.(alivewfntmp.eqv..true.)) then !delete old wfntmp folder
        write(*,*) "Running: rmdir /S /Q wfntmp"
        call system("rmdir /S /Q wfntmp")
    else if ((isys==2.or.isys==3).and.alivewfntmp.eqv..true.) then
        write(*,*) "Running: rm -rf wfntmp"
        call system("rm -rf wfntmp")
    end if
else if (iwfntmptype==2) then
    do i=1,9999 !Find a proper name of temporary folder
        write(c80tmp,"('wfntmp',i4.4)") i
        inquire(file='./c80tmp/.',exist=alivewfntmp)
        if (alivewfntmp.eqv..false.) exit
    end do
    if (isys==1) write(tmpdir,"('wfntmp',i4.4,'\')") i
    if (isys==2.or.isys==3) write(tmpdir,"('wfntmp',i4.4,'/')") i
end if
write(*,*) "Running: mkdir "//trim(c80tmp) !Build new temporary folder
call system("mkdir "//trim(c80tmp))
inquire(file="./atomwfn/.",exist=aliveatomwfn)
if (isys==1.and.aliveatomwfn.eqv..true.) then
    write(*,*) "Running: copy atomwfn\*.wfn "//trim(tmpdir)
    call system("copy atomwfn\*.wfn "//trim(tmpdir))
else if ((isys==2.or.isys==3).and.aliveatomwfn.eqv..true.) then
    write(*,*) "Running: cp atomwfn/*.wfn "//trim(tmpdir)
    call system("cp atomwfn/*.wfn "//trim(tmpdir))
end if

noatmwfn=0 !Check if the atomic wfn file have pre-stored in atomwfn folder, if not, invoke gaussian to calc it
do i=1,nfragatmnum
    if (isys==1) inquire(file="atomwfn\"//a(fragatm(i))%name//".wfn",exist=alive)
    if (isys==2.or.isys==3) inquire(file="atomwfn/"//a(fragatm(i))%name//".wfn",exist=alive)
    if (.not.alive) then
        noatmwfn=1
        exit
    end if
end do

if (noatmwfn==0) then
    write(*,"(a)") " All atom .wfn files needed have already presented in ""atomwfn"" folder, we will not calculate them"
else if (noatmwfn==1) then !Some or all atomic wfn don't exist, calc them
    !Select calculation level
    write(*,"(a)") " Note: Some or all atom .wfn files needed are not present in ""atomwfn"" folder, they must be calculated now. See Section 3.7.3 of the manual for detail."
    write(*,"(a)") " Now please input the level for calculating atom wfn files, theoretical method is optional."
    write(*,"(a)") " For example: 6-31G* or B3LYP/cc-pVDZ    You can also add other keywords at the same time, e.g. M052X/6-311G(2df,2p) scf=(vshift=500) int=ultrafine"
    read(*,"(a)") basisset !Note: 6d 10f is not required for generating wfn files, since the work has been done in L607 internally
    !Check Gaussian path
    inquire(file=gaupath,exist=alive)
    if (.not.alive) then
        write(*,*) "Couldn't find Gaussian path defined in ""gaupath"" variable in settings.ini"
        if (isys==1) write(*,*) "Input the path of Gaussian executable file, e.g. ""d:\study\g03w\g03.exe"""
        if (isys==2.or.isys==3) write(*,*) "Input the path of Gaussian executable file, e.g. ""/sob/g09/g09"""
        do while(.true.)
            read(*,"(a)") gaupath
            inquire(file=gaupath,exist=alive)
            if (alive) exit
            write(*,*) "Couldn't find Gaussian executable file, input again"
        end do
    end if
end if

!Generate .gjf file for all elements, regardless if their wfn file have already presented, meanwhile count the total number of elements
itype=0
do i=1,nfragatmnum
    inquire(file=trim(tmpdir)//a(fragatm(i))%name//".gjf",exist=alive)
    if (.not.alive) then
        itype=itype+1 !How many different types
        typename(itype)=a(fragatm(i))%name
                
        if (a_org(fragatm(i))%index>36) then
            inquire(file=trim(tmpdir)//a(fragatm(i))%name//".wfn",exist=alive)
            if (.not.alive) then !The wfn file of the heavy element hasn't been provided in "atomwfn" and hence cannot be found in "wfntmp" here
                write(*,"(a,a,a)") " Error: Multiwfn can't invoke Gaussian to generate wavefunction file and sphericalize density for ",a(fragatm(i))%name,", since its &
                index is larger than 36! You should provide corresponding atom .wfn files in ""atomwfn"" folder manually"
                write(*,*) "Press ENTER to continue"
                read (*,*)
                return
            end if
        end if
        
        open(14,file=trim(tmpdir)//a(fragatm(i))%name//".gjf",status="replace")
        !If user inputted including "/" e.g. B3LYP/6-31g*, will replace default theoretical method
        if (index(basisset,'/')==0) then
            if (a(fragatm(i))%index<=20.or.a(fragatm(i))%index>=31) then
                write(14,"(a,/)") "#T out=wfn ROHF/"//trim(basisset) !Main group elements. If not use scf=sp, in g09, RO calculations for IIIA elements are to converge
                write(14,"(a,/)") "Temporary file for promolecule, ROHF"//trim(basisset)
            else
                write(14,"(a,/)") "#T out=wfn UB3LYP/"//trim(basisset) !Transition metals
                write(14,"(a,/)") "Temporary file for promolecule, UB3LYP"//trim(basisset)
            end if
        else
            if (a(fragatm(i))%index<=20.or.a(fragatm(i))%index>=31) then
                write(14,"(a,/)") "#T out=wfn RO"//trim(basisset) !Main group elements
                write(14,"(a,/)") "Temporary file for promolecule, RO"//trim(basisset)
            else
                write(14,"(a,/)") "#T out=wfn U"//trim(basisset) !Transition metals (RO may leads to convergence problem)
                write(14,"(a,/)") "Temporary file for promolecule, U"//trim(basisset)
            end if
        end if

        !Currently support up to the fourth row
        if (a(fragatm(i))%name=="H ".or.a(fragatm(i))%name=="Li".or.a(fragatm(i))%name=="Na".or.a(fragatm(i))%name=="K") write(14,*) "0 2"
        if (a(fragatm(i))%name=="Be".or.a(fragatm(i))%name=="Mg".or.a(fragatm(i))%name=="Ca") write(14,*) "0 1"
        if (a(fragatm(i))%name=="B ".or.a(fragatm(i))%name=="Al".or.a(fragatm(i))%name=="Ga") write(14,*) "0 2"
        if (a(fragatm(i))%name=="C ".or.a(fragatm(i))%name=="Si".or.a(fragatm(i))%name=="Ge") then
            if (SpherIVgroup==0) write(14,*) "0 5"
            if (SpherIVgroup==1) write(14,*) "0 3"
        end if
        if (a(fragatm(i))%name=="N ".or.a(fragatm(i))%name=="P ".or.a(fragatm(i))%name=="As") write(14,*) "0 4"
        if (a(fragatm(i))%name=="O ".or.a(fragatm(i))%name=="S ".or.a(fragatm(i))%name=="Se") write(14,*) "0 3"
        if (a(fragatm(i))%name=="F ".or.a(fragatm(i))%name=="Cl".or.a(fragatm(i))%name=="Br") write(14,*) "0 2"
        if (a(fragatm(i))%name=="He".or.a(fragatm(i))%name=="Ne".or.a(fragatm(i))%name=="Ar".or.a(fragatm(i))%name=="Kr") write(14,*) "0 1"
        if (a(fragatm(i))%name=="Sc") write(14,*) "0 2" !3d1 4s2
        if (a(fragatm(i))%name=="Ti") write(14,*) "0 3" !3d2 4s2
        if (a(fragatm(i))%name=="V ") write(14,*) "0 4" !3d3 4s2
        if (a(fragatm(i))%name=="Cr") write(14,*) "0 7" !3d5 4s1, needn't sphericalization
        if (a(fragatm(i))%name=="Mn") write(14,*) "0 6" !3d5 4s2, needn't sphericalization
        if (a(fragatm(i))%name=="Fe") write(14,*) "0 5" !3d6 4s2
        if (a(fragatm(i))%name=="Co") write(14,*) "0 4" !3d7 4s2
        if (a(fragatm(i))%name=="Ni") write(14,*) "0 3" !3d8 4s2
        if (a(fragatm(i))%name=="Cu") write(14,*) "0 2" !3d10 4s1, needn't sphericalization
        if (a(fragatm(i))%name=="Zn") write(14,*) "0 1" !3d10 4s2, needn't sphericalization
        write(14,*) a(fragatm(i))%name,0.0,0.0,0.0
        write(14,*)
        write(14,*) trim(tmpdir)//a(fragatm(i))%name//".wfn" !The output path of wfn file
        write(14,*)
        write(14,*)
        close(14)
    end if
end do

if (noatmwfn==0) then
    if (isys==1) call system("del "//trim(tmpdir)//"*.gjf /Q") !The .gjf generated have valueless now, delete them for avoiding user's misunderstanding
    if (isys==2.or.isys==3) call system("rm "//trim(tmpdir)//"*.gjf -f")
else if (noatmwfn==1) then !Some wfn needs to be genereated by Gaussian and sphericalized here
    do i=1,nfragatmnum
        nametmp=a_org(fragatm(i))%name
        inquire(file=trim(tmpdir)//nametmp//".wfn",exist=alive)
        if (alive) cycle !If the .wfn file had copied from atomwfn folder, needn't recalculate

        write(*,*) "Running:"
        write(*,*) trim(Gaupath)//' "'//trim(tmpdir)//nametmp//'.gjf" "'//trim(tmpdir)//nametmp//'"'
        call system(trim(Gaupath)//' "'//trim(tmpdir)//nametmp//'.gjf" "'//trim(tmpdir)//nametmp//'"')
        !Check if Gaussian task was successfully finished
        if (isys==1) inquire(file=trim(tmpdir)//trim(nametmp)//".out",exist=alivegauout)
        if (isys==2.or.isys==3) inquire(file=trim(tmpdir)//trim(nametmp)//".log",exist=alivegauout)
        if (alivegauout) then
            if (isys==1) open(10,file=trim(tmpdir)//trim(nametmp)//".out",status="old")
            if (isys==2.or.isys==3) open(10,file=trim(tmpdir)//trim(nametmp)//".log",status="old")
            call loclabel(10,"Normal termination",igaunormal)
            close(10)
            if (igaunormal==0) then
                write(*,"(a)") "Gaussian running may be failed! Please manually check Gaussian input and output files in wfntmp folder. Press ENTER to continue"
                read (*,*)
            else if (igaunormal==1) then
                write(*,*) "Finished successfully!"
            end if
        else
            write(*,"(a)") "Gaussian running may be failed! Please manually check Gaussian input and output files in wfntmp folder"
            read (*,*)
        end if
    
        !Load and sphericalize electron density for the just generated wfn, and then save
        if (ispheratm==1.and.igaunormal==1) then
            call dealloall
            call readwfn(trim(tmpdir)//nametmp//".wfn",1)
            !Main group, restrict open-shell
            if (nametmp=="H ".or.nametmp=="Li".or.nametmp=="Na".or.nametmp=="K") MOocc(nmo)=1.0D0
            if (nametmp=="B ".or.nametmp=="Al".or.nametmp=="Ga") then
                nmo=nmo+2
                call resizebynmo(nmo,nprims) !Enlarge nmo by 2, but don't interfere nprims
                MOene(nmo-1:nmo)=MOene(nmo-2)
                MOtype(nmo-2:nmo)=1 !actually no use, because we only use atomic wfn. files to get total density
                MOocc(nmo-2:nmo)=1D0/3D0
                call orbcoeffrotate(nmo-2) !XYZ->ZXY, note: nmo-2 is original single occupied orbital
                CO(nmo-1,:)=CO(nmo-2,:)
                call orbcoeffrotate(nmo-2) !ZXY->YZX
                CO(nmo,:)=CO(nmo-2,:)
                call orbcoeffrotate(nmo-2) !YZX->XYZ, namely recovered
                !Now nmo-2,nmo-1,nmo correspond XYZ,ZXY,YZX
            end if
            if (nametmp=="C ".or.nametmp=="Si".or.nametmp=="Ge") then
                if (SpherIVgroup==0) then
                    MOocc(nmo-3:nmo)=1.0D0
                else if (SpherIVgroup==1) then
                    nmo=nmo+1
                    call resizebynmo(nmo,nprims)
                    MOene(nmo)=MOene(nmo-2) !MOene(nmo-1) is degenerate to MOene(nmo-2)
                    MOtype(nmo-2:nmo)=1
                    MOocc(nmo-2:nmo)=2D0/3D0
                    !Rotate and copy the first occupied p orbital (nmo-5)
                    call orbcoeffrotate(nmo-2) !XYZ->ZXY
                    CO(nmo-1,:)=CO(nmo-2,:) !Overlap the already occupied orbital
                    call orbcoeffrotate(nmo-2) !ZXY->YZX
                    CO(nmo,:)=CO(nmo-2,:)
                    call orbcoeffrotate(nmo-2) !YZX->XYZ, namely recovered
                end if
            end if
            if (nametmp=="N ".or.nametmp=="P ".or.nametmp=="As") MOocc(nmo-2:nmo)=1.0D0
            if (nametmp=="O ".or.nametmp=="S ".or.nametmp=="Se") MOocc(nmo-2:nmo)=4D0/3D0
            if (nametmp=="F ".or.nametmp=="Cl".or.nametmp=="Br") MOocc(nmo-2:nmo)=5D0/3D0
            !Transition metals, unrestrict open-shell, find boundary of alpha and beta first
            do ibound=2,nmo
                if (MOene(ibound)<MOene(ibound-1)) exit !from ii is beta orbitals
            end do
            !For Sc, Ti and V, rotate and duplicate d orbitals in each diection to get *near* spherical density, as for III main group
            !Note: Don't use Hartree-Fock, because correct energy sequence couldn't be reproduced, so can't be sphericalized correctly!
            if (nametmp=="Sc".or.nametmp=="Ti".or.nametmp=="V ") then !3d1 4s2, 3d2 4s2, 3d3 4s2
                if (nametmp=="Sc") then
                    ibeg=1 !alpha 4s orbital, because this s orbital shows very strong unequlitity
                    iend=2
                else if (nametmp=="Ti") then
                    ibeg=2
                    iend=3
                else if (nametmp=="V") then
                    ibeg=2
                    iend=4
                end if
                ienlarge=(iend-ibeg+1)*2
                call resizebynmo(nmo+ienlarge,nprims)
                ipass=0
                do iavgorb=ibeg,iend
                    call orbcoeffrotate(ibound-iavgorb) !rotate this orbital
                    CO(nmo+1+ipass,:)=CO(ibound-iavgorb,:) !Duplicate this orbital
                    call orbcoeffrotate(ibound-iavgorb)
                    CO(nmo+2+ipass,:)=CO(ibound-iavgorb,:)
                    call orbcoeffrotate(ibound-iavgorb) !recover
                    MOocc(ibound-iavgorb)=1D0/3D0
                    MOene(nmo+1+ipass:nmo+2+ipass)=MOene(ibound-iavgorb)
                    ipass=ipass+2 !next time skip nmo+1 and nmo+2
                end do
                MOocc(nmo+1:nmo+ienlarge)=1D0/3D0
                MOtype(nmo+1:nmo+ienlarge)=1
                nmo=nmo+ienlarge
            else if (nametmp=="Fe") then !3d6 4s2
                MOocc(nmo-1)=0D0 !delete the only d-beta orbital, the "nmo"th orbital is 4s-beta
                MOocc(ibound-6:ibound-2)=1.2D0 !Scatter one electron in beta-d orbital to alpha orbitals evenly. MOocc(ibound) is 4s orbital
            else if (nametmp=="Co") then !3d7 4s2
                MOocc(nmo-2:nmo-1)=0D0
                MOocc(ibound-6:ibound-2)=1.4D0
            else if (nametmp=="Ni") then !3d8 4s2
                MOocc(nmo-3:nmo-1)=0D0
                MOocc(ibound-6:ibound-2)=1.6D0
            end if
            call outwfn(trim(tmpdir)//nametmp//".wfn",0,0,10)
        end if
    end do
end if
write(*,*)

!Setup custom operation list
ncustommap=nfragatmnum
if (allocated(custommapname)) deallocate(custommapname)
if (allocated(customop)) deallocate(customop)
allocate(custommapname(ncustommap))
allocate(customop(ncustommap))
customop='-'
if (ipromol==1) customop='+'  !Do promolecular map

!Generate atomic wfn file from element wfn file, meanwhile take them into custom operation list
do i=1,itype !Scan each atomtype in current system
    call dealloall
    call readwfn(trim(tmpdir)//typename(i)//".wfn",1)
    do j=1,nfragatmnum
        if (a_org(fragatm(j))%name==typename(i)) then !Find atoms attributed to current element
            a(1)%x=a_org(fragatm(j))%x !Modify the atomic .wfn, then output to new .wfn
            a(1)%y=a_org(fragatm(j))%y
            a(1)%z=a_org(fragatm(j))%z
            write(outwfnname,"(a2,i4,a4)") typename(i),fragatm(j),".wfn"
            call outwfn(trim(tmpdir)//outwfnname,0,0,10)
            custommapname(j)=trim(tmpdir)//outwfnname !Sequence is identical to atom in fragment
        end if
    end do
end do

call dealloall
call readinfile(firstfilename,1)
! nmo=nmo_org !Original .wfn may be maniplated by user, retrieve them
! b=b_org
! CO=CO_org
! nprims=nprims_org
end subroutine



!!!------------------------- Generate density matrix, currently only for .fch file
subroutine genP
use defvar
use util
implicit real*8 (a-h,o-z)
integer i,j,imo
if (allocated(Ptot)) deallocate(Ptot)
if (allocated(Palpha)) deallocate(Palpha)
if (allocated(Pbeta)) deallocate(Pbeta)
allocate(Ptot(nbasis,nbasis))
if (wfntype==1.or.wfntype==2.or.wfntype==4) then !open-shell
    allocate(Palpha(nbasis,nbasis))
    allocate(Pbeta(nbasis,nbasis))
end if

if (wfntype==0.or.wfntype==2.or.wfntype==3) then !RHF,ROHF,restricted post-HF
    Ptot=0
    do imo=1,nmo
        if (MOocc(imo)==0D0) cycle
        Ptot=Ptot+MOocc(imo)*matmul(CObasa(:,imo:imo),transpose(CObasa(:,imo:imo)))
    end do
    if (wfntype==2) then !ROHF
        Palpha=0D0
!Add the singly occupied orbital contribution to Ptot, as if all orbital is doubled occupied, divided by 2 yielding Palpha
        do imo=1,nmo
            if (MOtype(imo)==1) then
                Palpha=Palpha+MOocc(imo)*matmul(CObasa(:,imo:imo),transpose(CObasa(:,imo:imo)))
            end if
        end do
        Palpha=(Palpha+Ptot)/2D0
        Pbeta=Ptot-Palpha
    end if
else if (wfntype==1.or.wfntype==4) then
    Palpha=0D0
    Pbeta=0D0
    do imo=1,nbasis
        if (MOocc(imo)==0D0) cycle
        Palpha=Palpha+MOocc(imo)*matmul(CObasa(:,imo:imo),transpose(CObasa(:,imo:imo)))
    end do
    do imo=1,nbasis
        if (MOocc(imo+nbasis)==0D0) cycle
        Pbeta=Pbeta+MOocc(imo+nbasis)*matmul(CObasb(:,imo:imo),transpose(CObasb(:,imo:imo)))
    end do
    Ptot=Palpha+Pbeta
end if
end subroutine



!!!------------------ Evaluate overlap integral for two unnormalized GTFs, a warpper of doSintactual for simplicity
real*8 function doSint(iGTF,jGTF)
integer iGTF,jGTF
real*8,external :: doSintactual
doSint=doSintactual(iGTF,jGTF,0,0,0,0,0,0)
end function
!!!------------------ Evaluate overlap integral for two unnormalized GTFs
!~p arguments are the shifts of GTF index, used by doKint, doveloint but not by doSint
real*8 function doSintactual(iGTF,jGTF,ix1p,iy1p,iz1p,ix2p,iy2p,iz2p)
use util
use defvar
implicit real*8(a-h,o-z)
integer iGTF,jGTF,ix1p,iy1p,iz1p,ix2p,iy2p,iz2p
x1=a(b(iGTF)%center)%x
y1=a(b(iGTF)%center)%y
z1=a(b(iGTF)%center)%z
x2=a(b(jGTF)%center)%x
y2=a(b(jGTF)%center)%y
z2=a(b(jGTF)%center)%z
ee1=b(iGTF)%exp
ee2=b(jGTF)%exp
ep=ee1+ee2
sqrtep=dsqrt(ep)
px=(ee1*x1+ee2*x2)/ep
py=(ee1*y1+ee2*y2)/ep
pz=(ee1*z1+ee2*z2)/ep        
expterm=dexp( -ee1*ee2*((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)/ep )
ix1=type2ix(b(iGTF)%functype)+ix1p
iy1=type2iy(b(iGTF)%functype)+iy1p
iz1=type2iz(b(iGTF)%functype)+iz1p
ix2=type2ix(b(jGTF)%functype)+ix2p
iy2=type2iy(b(jGTF)%functype)+iy2p
iz2=type2iz(b(jGTF)%functype)+iz2p
!chen book,P103
numx=ceiling( (ix1+ix2+1)/2D0 ) !Need to calculate n points
sx=0.0D0
do i=1,numx
    tmp=Rhm(numx,i)/sqrtep+px
    term1=(tmp-x1)**ix1
    term2=(tmp-x2)**ix2
    sx=sx+Whm(numx,i)*term1*term2
end do
sx=sx/sqrtep

numy=ceiling( (iy1+iy2+1)/2D0 )
sy=0.0D0
do i=1,numy
    tmp=Rhm(numy,i)/sqrtep+py
    term1=(tmp-y1)**iy1
    term2=(tmp-y2)**iy2
    sy=sy+Whm(numy,i)*term1*term2
end do
sy=sy/sqrtep

numz=ceiling( (iz1+iz2+1)/2D0 )
sz=0.0D0
do i=1,numz
    tmp=Rhm(numz,i)/sqrtep+pz
    term1=(tmp-z1)**iz1
    term2=(tmp-z2)**iz2
    sz=sz+Whm(numz,i)*term1*term2
end do
sz=sz/sqrtep

doSintactual=sx*sy*sz*expterm
end function
!!!------------------ Generate overlap matrix between all GTFs
!nsize should be nprims*(nprims+1)/2
subroutine genGTFSmat(GTFSmat,nsize)
use defvar
implicit real*8 (a-h,o-z)
integer nsize
real*8 GTFSmat(nsize)
!$OMP PARALLEL DO SHARED(GTFSmat) PRIVATE(ides,iGTF,jGTF) schedule(dynamic) NUM_THREADS( nthreads  )
do iGTF=1,nprims
    do jGTF=iGTF,nprims
        ides=jGTF*(jGTF-1)/2+iGTF
        GTFSmat(ides)=doSint(iGTF,jGTF)
    end do
end do
!$OMP END PARALLEL DO
end subroutine
!!!------------------ Generate overlap matrix between all basis functions
!Sbas should be allocated first. The resultant matrix is for Cartesian basis functions, may be converted to spherical-harmonic later
subroutine genSbas
use defvar
implicit real*8 (a-h,o-z)
Sbas=0D0
!$OMP PARALLEL DO SHARED(Sbas) PRIVATE(i,ii,j,jj) schedule(dynamic) NUM_THREADS( nthreads  )
do i=1,nbasis
    do j=i,nbasis
        do ii=primstart(i),primend(i)
            do jj=primstart(j),primend(j)
                Sbas(i,j)=Sbas(i,j)+primconnorm(ii)*primconnorm(jj)*doSint(ii,jj)
            end do
        end do
        Sbas(j,i)=Sbas(i,j)
    end do
end do
!$OMP END PARALLEL DO
end subroutine



!!!-------- Evaluate dipole moment integral for two unnormalized GTFs. The negative charge of electron has been considered!
!~p arguments are the shifts of GTF index as doSintactual
!xint/yint/zint correspond to dipole moment integral in X/Y/Z
subroutine dodipoleint(iGTF,jGTF,ix1p,iy1p,iz1p,ix2p,iy2p,iz2p,xint,yint,zint)
use util
use defvar
implicit real*8(a-h,o-z)
real*8 xint,yint,zint
integer iGTF,jGTF,ix1p,iy1p,iz1p,ix2p,iy2p,iz2p
x1=a(b(iGTF)%center)%x
y1=a(b(iGTF)%center)%y
z1=a(b(iGTF)%center)%z
x2=a(b(jGTF)%center)%x
y2=a(b(jGTF)%center)%y
z2=a(b(jGTF)%center)%z
ee1=b(iGTF)%exp
ee2=b(jGTF)%exp
ep=ee1+ee2
sqrtep=dsqrt(ep)
px=(ee1*x1+ee2*x2)/ep
py=(ee1*y1+ee2*y2)/ep
pz=(ee1*z1+ee2*z2)/ep        
expterm=dexp( -ee1*ee2*((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)/ep )
ix1=type2ix(b(iGTF)%functype)+ix1p
iy1=type2iy(b(iGTF)%functype)+iy1p
iz1=type2iz(b(iGTF)%functype)+iz1p
ix2=type2ix(b(jGTF)%functype)+ix2p
iy2=type2iy(b(jGTF)%functype)+iy2p
iz2=type2iz(b(jGTF)%functype)+iz2p
!First, calculate sx,sy,sz as usual as doSint
numx=ceiling( (ix1+ix2+1)/2D0 ) !Need to calculate n points
sx=0.0D0
do i=1,numx
    tmp=Rhm(numx,i)/sqrtep+px
    term1=(tmp-x1)**ix1
    term2=(tmp-x2)**ix2
    sx=sx+Whm(numx,i)*term1*term2
end do
sx=sx/sqrtep
numy=ceiling( (iy1+iy2+1)/2D0 )
sy=0.0D0
do i=1,numy
    tmp=Rhm(numy,i)/sqrtep+py
    term1=(tmp-y1)**iy1
    term2=(tmp-y2)**iy2
    sy=sy+Whm(numy,i)*term1*term2
end do
sy=sy/sqrtep
numz=ceiling( (iz1+iz2+1)/2D0 )
sz=0.0D0
do i=1,numz
    tmp=Rhm(numz,i)/sqrtep+pz
    term1=(tmp-z1)**iz1
    term2=(tmp-z2)**iz2
    sz=sz+Whm(numz,i)*term1*term2
end do
sz=sz/sqrtep
!Second, calculate overlap integral in X,Y,Z directions but with X,Y,Z coordinate variables (relative to the original point of the whole system) to produce sxx,syy,szz
numx=ceiling( (ix1+ix2+2)/2D0 ) !Because X variable is introduced, ix1+ix2+2 is used instead of ix1+ix2+1
sxx=0.0D0
do i=1,numx
    tmp=Rhm(numx,i)/sqrtep+px
    term1=(tmp-x1)**ix1
    term2=(tmp-x2)**ix2
    sxx=sxx+Whm(numx,i)*term1*term2*tmp
end do
sxx=sxx/sqrtep
numy=ceiling( (iy1+iy2+2)/2D0 )
syy=0.0D0
do i=1,numy
    tmp=Rhm(numy,i)/sqrtep+py
    term1=(tmp-y1)**iy1
    term2=(tmp-y2)**iy2
    syy=syy+Whm(numy,i)*term1*term2*tmp
end do
syy=syy/sqrtep
numz=ceiling( (iz1+iz2+2)/2D0 )
szz=0.0D0
do i=1,numz
    tmp=Rhm(numz,i)/sqrtep+pz
    term1=(tmp-z1)**iz1
    term2=(tmp-z2)**iz2
    szz=szz+Whm(numz,i)*term1*term2*tmp
end do
szz=szz/sqrtep

xint=-sxx*sy*sz*expterm
yint=-sx*syy*sz*expterm
zint=-sx*sy*szz*expterm
end subroutine
!------ A warpper of dodipoleint, used to directly get a single component of dipole moment integral. icomp=1/2/3 corresponds to X/Y/Z component
real*8 function dipintcomp(icomp,iGTF,jGTF,ix1p,iy1p,iz1p,ix2p,iy2p,iz2p)
integer icomp,iGTF,jGTF,ix1p,iy1p,iz1p,ix2p,iy2p,iz2p
real*8 xcomp,ycomp,zcomp
call dodipoleint(iGTF,jGTF,ix1p,iy1p,iz1p,ix2p,iy2p,iz2p,xcomp,ycomp,zcomp)
if (icomp==1) dipintcomp=xcomp
if (icomp==2) dipintcomp=ycomp
if (icomp==3) dipintcomp=zcomp
end function
!!!--------------- Generate dipole moment integral matrix between all GTFs
!nsize should be nprims*(nprims+1)/2
subroutine genGTFDmat(GTFdipmat,nsize)
use defvar
implicit real*8 (a-h,o-z)
integer nsize
real*8 GTFdipmat(3,nsize)
!$OMP PARALLEL DO SHARED(GTFdipmat) PRIVATE(ides,iGTF,jGTF,xdiptmp,ydiptmp,zdiptmp) schedule(dynamic) NUM_THREADS( nthreads  )
do iGTF=1,nprims
    do jGTF=iGTF,nprims
        ides=jGTF*(jGTF-1)/2+iGTF
        call dodipoleint(iGTF,jGTF,0,0,0,0,0,0,xdiptmp,ydiptmp,zdiptmp)
        GTFdipmat(1,ides)=xdiptmp
        GTFdipmat(2,ides)=ydiptmp
        GTFdipmat(3,ides)=zdiptmp
    end do
end do
!$OMP END PARALLEL DO
end subroutine
!!!------------------ Generate dipole moment integral matrix between all basis functions
!Dbas should be allocated first. The resultant matrix is for Cartesian basis functions, may be converted to spherical-harmonic later
subroutine genDbas
use defvar
implicit real*8 (a-h,o-z)
Dbas=0D0
!$OMP PARALLEL DO SHARED(Dbas) PRIVATE(i,ii,j,jj,xdiptmp,ydiptmp,zdiptmp) schedule(dynamic) NUM_THREADS( nthreads  )
do i=1,nbasis
    do j=i,nbasis
        do ii=primstart(i),primend(i)
            do jj=primstart(j),primend(j)
                call dodipoleint(ii,jj,0,0,0,0,0,0,xdiptmp,ydiptmp,zdiptmp)
                Dbas(1,i,j)=Dbas(1,i,j)+primconnorm(ii)*primconnorm(jj)*xdiptmp
                Dbas(2,i,j)=Dbas(2,i,j)+primconnorm(ii)*primconnorm(jj)*ydiptmp
                Dbas(3,i,j)=Dbas(3,i,j)+primconnorm(ii)*primconnorm(jj)*zdiptmp
            end do
        end do
        Dbas(:,j,i)=Dbas(:,i,j)
    end do
end do
!$OMP END PARALLEL DO
end subroutine



!!!------- Evaluate magnetic integral for two unnormalized GTFs 
!The imaginary sign i is ignored. Note that the negative sign of magnetic operator is not occurred here
!Consult doVelint for the method for evaluation of <a|d/dx|b>, and TCA,6,341 for the formula of magnetic integral
subroutine doMagint(iGTF,jGTF,xcomp,ycomp,zcomp)
use defvar
implicit real*8(a-h,o-z)
integer iGTF,jGTF
real*8 xcomp,ycomp,zcomp,term(4)
ee1=b(iGTF)%exp
ee2=b(jGTF)%exp
ix1=type2ix(b(iGTF)%functype)
iy1=type2iy(b(iGTF)%functype)
iz1=type2iz(b(iGTF)%functype)
ix2=type2ix(b(jGTF)%functype)
iy2=type2iy(b(jGTF)%functype)
iz2=type2iz(b(jGTF)%functype)
term=0
!X component, <a|y*d/dz-z*d/dy|b>. Since <a|d/dz|b>=iz2*<a|b-1z>-2*ee2*<a|b+1z>, the term such as <a|y*d/dz|b> can be evaluated in terms of dipole integrals iz2*<a|y|b-1z>-2*ee2*<a|y|b+1z>
if(iz2>0) term(1)=   iz2*dipintcomp(2,iGTF,jGTF,0,0,0,0,0,-1) !viz. iz2*<a|y|b-1z>
          term(2)=-2*ee2*dipintcomp(2,iGTF,jGTF,0,0,0,0,0, 1)
if(iy2>0) term(3)=  -iy2*dipintcomp(3,iGTF,jGTF,0,0,0,0,-1,0)
          term(4)= 2*ee2*dipintcomp(3,iGTF,jGTF,0,0,0,0, 1,0)
xcomp=-sum(term) !Note that the result of dipintcomp has a negative sign due to the negative charge of electron, so here revise the sign
term=0
!Y component, <a|z*d/dx-x*d/dz|b>
if(ix2>0) term(1)=   ix2*dipintcomp(3,iGTF,jGTF,0,0,0,-1,0,0)
          term(2)=-2*ee2*dipintcomp(3,iGTF,jGTF,0,0,0, 1,0,0)
if(iz2>0) term(3)=  -iz2*dipintcomp(1,iGTF,jGTF,0,0,0,0,0,-1)
          term(4)= 2*ee2*dipintcomp(1,iGTF,jGTF,0,0,0,0,0, 1)
ycomp=-sum(term)
term=0
!Z component, <a|x*d/dy-y*d/dx|b>
if(iy2>0) term(1)=   iy2*dipintcomp(1,iGTF,jGTF,0,0,0,0,-1,0)
          term(2)=-2*ee2*dipintcomp(1,iGTF,jGTF,0,0,0,0, 1,0)
if(ix2>0) term(3)=  -ix2*dipintcomp(2,iGTF,jGTF,0,0,0,-1,0,0)
          term(4)= 2*ee2*dipintcomp(2,iGTF,jGTF,0,0,0, 1,0,0)
zcomp=-sum(term)
end subroutine
!!!--------------- Generate magnetic dipole moment integral matrix between all GTFs
!nsize should be nprims*(nprims+1)/2
!Beware that when using this result, (j,i) element should be set to negative value of (i,j) due to the Hermitean character of this operator!
subroutine genGTFMmat(GTFdipmat,nsize)
use defvar
implicit real*8 (a-h,o-z)
integer nsize
real*8 GTFdipmat(3,nsize)
GTFdipmat=0D0
!$OMP PARALLEL DO SHARED(GTFdipmat) PRIVATE(ides,iGTF,jGTF,xdiptmp,ydiptmp,zdiptmp) schedule(dynamic) NUM_THREADS( nthreads  )
do iGTF=1,nprims
    do jGTF=iGTF+1,nprims !For iGTF=jGTF, the value must exactly zero, so don't calculate
        ides=jGTF*(jGTF-1)/2+iGTF
        call doMagint(iGTF,jGTF,xdiptmp,ydiptmp,zdiptmp)
        GTFdipmat(1,ides)=xdiptmp
        GTFdipmat(2,ides)=ydiptmp
        GTFdipmat(3,ides)=zdiptmp
    end do
end do
!$OMP END PARALLEL DO
end subroutine
!!!------------------ Generate magnetic integral matrix between all basis functions
!Magbas should be allocated first. The resultant matrix is for Cartesian basis functions, may be converted to spherical-harmonic later
!Notice that the diagonal element of magnetic integral matrix is zero, and (i,j)=-(j,i) due to Hermitean character
subroutine genMagbas
use defvar
implicit real*8 (a-h,o-z)
Magbas=0D0
!$OMP PARALLEL DO SHARED(Magbas) PRIVATE(i,ii,j,jj,xcomp,ycomp,zcomp) schedule(dynamic) NUM_THREADS( nthreads  )
do i=1,nbasis
    do j=i+1,nbasis
        do ii=primstart(i),primend(i)
            do jj=primstart(j),primend(j)
                call doMagint(ii,jj,xcomp,ycomp,zcomp)
                Magbas(1,i,j)=Magbas(1,i,j)+primconnorm(ii)*primconnorm(jj)*xcomp
                Magbas(2,i,j)=Magbas(2,i,j)+primconnorm(ii)*primconnorm(jj)*ycomp
                Magbas(3,i,j)=Magbas(3,i,j)+primconnorm(ii)*primconnorm(jj)*zcomp
            end do
        end do
        Magbas(:,j,i)=-Magbas(:,i,j)
    end do
end do
!$OMP END PARALLEL DO
end subroutine



!!!------- Evaluate velocity integral for two unnormalized GTFs 
!There are three components. e.g. X direction: i<a|d/dx|b>. The imaginary sign i is ignored. Note that the negative sign of momentum operator is not occurred here
!One can consult p97 of Chen's book for the derivative of GTF. Namely <a|d/dx|b>=ix2*<a|b-1x>-2*ee2*<a|b+1x>
subroutine doVelint(iGTF,jGTF,xcomp,ycomp,zcomp)
use defvar
implicit real*8(a-h,o-z)
integer iGTF,jGTF
real*8 xcomp,ycomp,zcomp
ee1=b(iGTF)%exp
ee2=b(jGTF)%exp
ix1=type2ix(b(iGTF)%functype)
iy1=type2iy(b(iGTF)%functype)
iz1=type2iz(b(iGTF)%functype)
ix2=type2ix(b(jGTF)%functype)
iy2=type2iy(b(jGTF)%functype)
iz2=type2iz(b(jGTF)%functype)
term1=0
if(ix2>0) term1=   ix2*doSintactual(iGTF,jGTF,0,0,0,-1,0,0)
          term2=-2*ee2*doSintactual(iGTF,jGTF,0,0,0, 1,0,0)
xcomp=term1+term2
term1=0
if(iy2>0) term1=   iy2*doSintactual(iGTF,jGTF,0,0,0,0,-1,0)
          term2=-2*ee2*doSintactual(iGTF,jGTF,0,0,0,0, 1,0)
ycomp=term1+term2
term1=0
if(iz2>0) term1=   iz2*doSintactual(iGTF,jGTF,0,0,0,0,0,-1)
          term2=-2*ee2*doSintactual(iGTF,jGTF,0,0,0,0,0, 1)
zcomp=term1+term2
end subroutine
!!!--------------- Generate velocity integral matrix between all GTFs
!nsize should be nprims*(nprims+1)/2
!Beware that when using this result, (j,i) element should be set to negative value of (i,j) due to the Hermitean character of this operator!
subroutine genGTFVelmat(GTFVelmat,nsize)
use defvar
implicit real*8 (a-h,o-z)
integer nsize
real*8 GTFVelmat(3,nsize)
!$OMP PARALLEL DO SHARED(GTFVelmat) PRIVATE(ides,iGTF,jGTF,xtmp,ytmp,ztmp) schedule(dynamic) NUM_THREADS( nthreads  )
do iGTF=1,nprims
    do jGTF=iGTF,nprims
        ides=jGTF*(jGTF-1)/2+iGTF
        call doVelint(iGTF,jGTF,xtmp,ytmp,ztmp)
        GTFVelmat(1,ides)=xtmp
        GTFVelmat(2,ides)=ytmp
        GTFVelmat(3,ides)=ztmp
    end do
end do
!$OMP END PARALLEL DO
end subroutine
!!!------------------ Generate velocity integral matrix between all basis functions
!Velbas should be allocated first. The resultant matrix is for Cartesian basis functions, may be converted to spherical-harmonic later
!Notice that the diagonal element of velocity integral matrix is zero, and (i,j)=-(j,i) due to Hermitean character
subroutine genVelbas
use defvar
implicit real*8 (a-h,o-z)
Velbas=0D0
!$OMP PARALLEL DO SHARED(Velbas) PRIVATE(i,ii,j,jj,xcomp,ycomp,zcomp) schedule(dynamic) NUM_THREADS( nthreads  )
do i=1,nbasis
    do j=i+1,nbasis
        do ii=primstart(i),primend(i)
            do jj=primstart(j),primend(j)
                call doVelint(ii,jj,xcomp,ycomp,zcomp)
                Velbas(1,i,j)=Velbas(1,i,j)+primconnorm(ii)*primconnorm(jj)*xcomp
                Velbas(2,i,j)=Velbas(2,i,j)+primconnorm(ii)*primconnorm(jj)*ycomp
                Velbas(3,i,j)=Velbas(3,i,j)+primconnorm(ii)*primconnorm(jj)*zcomp
            end do
        end do
        Velbas(:,j,i)=-Velbas(:,i,j)
    end do
end do
!$OMP END PARALLEL DO
end subroutine



!!!------------------- Evaluate kinetic integral (i.e. -(1/2)der2 )for two unnormalized GTFs, see Chen's book, p104
real*8 function doTint(iGTF,jGTF)
use defvar
implicit real*8(a-h,o-z)
integer iGTF,jGTF
real*8 term(4)
ee1=b(iGTF)%exp
ee2=b(jGTF)%exp
ix1=type2ix(b(iGTF)%functype)
iy1=type2iy(b(iGTF)%functype)
iz1=type2iz(b(iGTF)%functype)
ix2=type2ix(b(jGTF)%functype)
iy2=type2iy(b(jGTF)%functype)
iz2=type2iz(b(jGTF)%functype)
term=0
if(ix1>0.and.ix2>0)  term(1)=   ix1*ix2*doSintactual(iGTF,jGTF,-1,0,0,-1,0,0)
if(ix1>0)            term(2)=-2*ee2*ix1*doSintactual(iGTF,jGTF,-1,0,0, 1,0,0)
if(ix2>0)            term(3)=-2*ee1*ix2*doSintactual(iGTF,jGTF, 1,0,0,-1,0,0)
                     term(4)= 4*ee1*ee2*doSintactual(iGTF,jGTF, 1,0,0, 1,0,0)
Tx=sum(term)
term=0
if(iy1>0.and.iy2>0)  term(1)=   iy1*iy2*doSintactual(iGTF,jGTF,0,-1,0,0,-1,0)
if(iy1>0)            term(2)=-2*ee2*iy1*doSintactual(iGTF,jGTF,0,-1,0,0, 1,0)
if(iy2>0)            term(3)=-2*ee1*iy2*doSintactual(iGTF,jGTF,0, 1,0,0,-1,0)
                     term(4)= 4*ee1*ee2*doSintactual(iGTF,jGTF,0, 1,0,0, 1,0)
Ty=sum(term)
term=0
if(iz1>0.and.iz2>0)  term(1)=   iz1*iz2*doSintactual(iGTF,jGTF,0,0,-1,0,0,-1)
if(iz1>0)            term(2)=-2*ee2*iz1*doSintactual(iGTF,jGTF,0,0,-1,0,0, 1)
if(iz2>0)            term(3)=-2*ee1*iz2*doSintactual(iGTF,jGTF,0,0, 1,0,0,-1)
                     term(4)= 4*ee1*ee2*doSintactual(iGTF,jGTF,0,0, 1,0,0, 1)
Tz=sum(term)
doTint=(Tx+Ty+Tz)/2
end function
!!!------------------ Generate kinetic energy matrix between all GTFs
!nsize should be nprims*(nprims+1)/2
subroutine genGTFTmat(GTFTmat,nsize)
use defvar
implicit real*8 (a-h,o-z)
integer nsize
real*8 GTFTmat(nsize)
!$OMP PARALLEL DO SHARED(GTFTmat) PRIVATE(ides,iGTF,jGTF) schedule(dynamic) NUM_THREADS( nthreads  )
do iGTF=1,nprims
    do jGTF=iGTF,nprims
        ides=jGTF*(jGTF-1)/2+iGTF
        GTFTmat(ides)=doTint(iGTF,jGTF)
    end do
end do
!$OMP END PARALLEL DO
end subroutine
!!!------------------ Generate kinetic energy matrix between all basis functions
!Tbas should be allocated first. The resultant matrix is for Cartesian basis functions, may be converted to spherical-harmonic later
subroutine genTbas
use defvar
implicit real*8 (a-h,o-z)
Tbas=0D0
!$OMP PARALLEL DO SHARED(Tbas) PRIVATE(i,ii,j,jj) schedule(dynamic) NUM_THREADS( nthreads  )
do i=1,nbasis
    do j=i,nbasis
        do ii=primstart(i),primend(i)
            do jj=primstart(j),primend(j)
                Tbas(i,j)=Tbas(i,j)+primconnorm(ii)*primconnorm(jj)*doTint(ii,jj)
            end do
        end do
        Tbas(j,i)=Tbas(i,j)
    end do
end do
!$OMP END PARALLEL DO
end subroutine



!!!------ Show system one-electron properties based on density matrix and integral matrix between basis functions
!The results are correct only when Cartesian basis functions are used
subroutine sys1eprop
use defvar
if (.not.allocated(Sbas)) allocate(Sbas(nbasis,nbasis))
call genSbas
write(*,"(' Total number of electrons:',f16.8)") sum(Ptot*Sbas)
if (.not.allocated(Tbas)) allocate(Tbas(nbasis,nbasis))
call genTbas
write(*,"(' Kinetic energy:',f18.9,' a.u.')") sum(Ptot*Tbas)
if (.not.allocated(Dbas)) allocate(Dbas(3,nbasis,nbasis))
call genDbas
write(*,"(' Electric dipole moment in X/Y/Z:',3f13.7,' a.u.')") sum(Ptot*Dbas(1,:,:)),sum(Ptot*Dbas(2,:,:)),sum(Ptot*Dbas(3,:,:))
if (.not.allocated(Magbas)) allocate(Magbas(3,nbasis,nbasis))
call genMagbas
write(*,"(' Magnetic dipole moment in X/Y/Z:',3f13.7,' a.u.')") sum(Ptot*Magbas(1,:,:)),sum(Ptot*Magbas(2,:,:)),sum(Ptot*Magbas(3,:,:))
if (.not.allocated(Velbas)) allocate(Velbas(3,nbasis,nbasis))
call genVelbas
write(*,"(' Linear momentum in X/Y/Z:       ',3f13.7,' a.u.')") sum(Ptot*Velbas(1,:,:)),sum(Ptot*Velbas(2,:,:)),sum(Ptot*Velbas(3,:,:))
end subroutine



!!!------------- Show all properties at a point
!ifuncsel: controls the gradient and Hessian for which function
!ifileid: Output to which file destination, of course 6=screen
subroutine showptprop(inx,iny,inz,ifuncsel,ifileid)
use util
use defvar
use function
implicit real*8(a-h,o-z)
real*8 inx,iny,inz
real*8 elehess(3,3),eigvecmat(3,3),eigval(3),elegrad(3),funchess(3,3),funcgrad(3),tmparr(3,1)
integer ifuncsel,ifileid
if (allocated(b)) then !If loaded file contains wavefuntion information
    call gencalchessmat(2,1,inx,iny,inz,elerho,elegrad,elehess) !Generate electron density, gradient and hessian
    write(ifileid,"(' Density of all electrons:',E18.10)") elerho
    if (ipolarpara==0) then
        tmpval=fspindens(inx,iny,inz,'s')
        write(ifileid,"(' Density of Alpha electrons:',E18.10)") (elerho+tmpval)/2D0
        write(ifileid,"(' Density of Beta electrons:',E18.10)") (elerho-tmpval)/2D0
        write(ifileid,"(' Spin density of electrons:',E18.10)") tmpval
    else if (ipolarpara==1) then
        write(ifileid,"(' Spin polarization parameter function:',E18.10)") fspindens(inx,iny,inz,'s')
    end if
    valG=lagkin(inx,iny,inz,0)
    valGx=lagkin(inx,iny,inz,1)
    valGy=lagkin(inx,iny,inz,2)
    valGz=lagkin(inx,iny,inz,3)
    write(ifileid,"(' Lagrangian kinetic energy G(r):',E18.10)") valG
    write(ifileid,"(' G(r) in X,Y,Z:',3E18.10)") valGx,valGy,valGz
    valK=Hamkin(inx,iny,inz,0)
!     valKx=Hamkin(inx,iny,inz,1)
!     valKy=Hamkin(inx,iny,inz,2)
!     valKz=Hamkin(inx,iny,inz,3)
    write(ifileid,"(' Hamiltonian kinetic energy K(r):',E18.10)") valK
!     write(ifileid,"(' K(r) in X,Y,Z:',3E18.10)") valKx,valKy,valKz
    write(ifileid,"(' Potential energy density V(r):',E18.10)") flapl(inx,iny,inz,'t')/4.0D0-2*valG
    write(ifileid,"(' Energy density E(r) or H(r):',E18.10)") -valK
    write(ifileid,"(' Laplacian of electron density:',E18.10)") laplfac*(elehess(1,1)+elehess(2,2)+elehess(3,3))
    write(ifileid,"(' Electron localization function (ELF):',E18.10)") ELF_LOL(inx,iny,inz,"ELF")
    write(ifileid,"(' Localized orbital locator (LOL):',E18.10)") ELF_LOL(inx,iny,inz,"LOL")
    write(ifileid,"(' Local information entropy:',E18.10)") infoentro(1,inx,iny,inz)
    write(ifileid,"(' Reduced density gradient (RDG):',E18.10)") fgrad(inx,iny,inz,'r')
    write(ifileid,"(' Reduced density gradient with promolecular approximation:',E18.10)") RDGprodens(inx,iny,inz)
    write(ifileid,"(' Sign(lambda2)*rho:',E18.10)") signlambda2rho(inx,iny,inz)
    write(ifileid,"(' Sign(lambda2)*rho with promolecular approximation:',E18.10)") signlambda2rho_prodens(inx,iny,inz)
    if (pairfunctype==1) write(ifileid,"(a,3f10.5,' :',E18.10)") " Corr. hole for alpha, ref.:",refx,refy,refz,pairfunc(inx,iny,inz)
    if (pairfunctype==2) write(ifileid,"(a,3f10.5,' :',E18.10)") " Corr. hole for beta, ref.:",refx,refy,refz,pairfunc(inx,iny,inz)
    if (pairfunctype==4) write(ifileid,"(a,3f10.5,' :',E18.10)") " Corr. fac. for alpha, ref.:",refx,refy,refz,pairfunc(inx,iny,inz)
    if (pairfunctype==5) write(ifileid,"(a,3f10.5,' :',E18.10)") " Corr. fac. for beta, ref.:",refx,refy,refz,pairfunc(inx,iny,inz)
    if (pairfunctype==7) write(ifileid,"(a,3f10.5,' :',E18.10)") " Exc.-corr. dens. for alpha, ref:",refx,refy,refz,pairfunc(inx,iny,inz)
    if (pairfunctype==8) write(ifileid,"(a,3f10.5,' :',E18.10)") " Exc.-corr. dens. for beta, ref:",refx,refy,refz,pairfunc(inx,iny,inz)
    write(ifileid,"(' Source function, ref.:',3f10.5,' :',E18.10)") refx,refy,refz,srcfunc(inx,iny,inz,srcfuncmode)
    write(ifileid,"(' Wavefunction value for orbital',i10,' :',E18.10)") iorbsel,fmo(inx,iny,inz,iorbsel)
    if (iALIEdecomp==0) then
        write(ifileid,"(' Average local ionization energy:',E18.10)") avglocion(inx,iny,inz)
    else if (iALIEdecomp==1) then
        call avglociondecomp(ifileid,inx,iny,inz)
    end if
    write(ifileid,"(' User defined real space function:',E18.10)") userfunc(inx,iny,inz)
    fesptmp=nucesp(inx,iny,inz)
    if (ifiletype==4) then
        write(ifileid,"(' ESP from atomic charges:',E18.10)") fesptmp
    else
        write(ifileid,"(' ESP from nuclear charges:',E18.10)") fesptmp
    end if
    if (ishowptESP==1) then
        fesptmpelec=eleesp(inx,iny,inz)
        write(ifileid,"(' ESP from electrons:',E18.10)") fesptmpelec
        write(ifileid,"(' Total ESP:',E18.10,' a.u. (',E14.7,' J/C,',E14.7,' kcal/mol)')") fesptmpelec+fesptmp,(fesptmpelec+fesptmp)*27.2113838D0,(fesptmpelec+fesptmp)*au2kcal
    end if
    write(ifileid,*)
    if (ifuncsel==1) then
        write(ifileid,*) "Note: Below information are for electron density"
        funchess=elehess
        funcgrad=elegrad
    else
        if (ifuncsel==3) write(ifileid,*) "Note: Below information are for Laplacian of electron density"
        if (ifuncsel==4) write(ifileid,*) "Note: Below information are for value of orbital wavefunction"
        if (ifuncsel==9) write(ifileid,*) "Note: Below information are for electron localization function"
        if (ifuncsel==10) write(ifileid,*) "Note: Below information are for localized orbital locator"
        if (ifuncsel==12) write(ifileid,*) "Note: Below information are for total ESP"
        if (ifuncsel==100) write(ifileid,*) "Note: Below information are for user defined real space function"
        call gencalchessmat(2,ifuncsel,inx,iny,inz,funcvalue,funcgrad,funchess)
    end if
    write(ifileid,*)
    write(ifileid,*) "Components of gradient in x/y/z are:"
    write(ifileid,"(3E18.10)") funcgrad(1),funcgrad(2),funcgrad(3)
    write(ifileid,"(' Norm of gradient is:',E18.10)") dsqrt(sum(funcgrad**2))
    write(ifileid,*)
    write(ifileid,*) "Components of Laplacian in x/y/z are:"
    write(ifileid,"(3E18.10)") funchess(1,1),funchess(2,2),funchess(3,3)
    write(ifileid,"(' Total:',E18.10)") funchess(1,1)+funchess(2,2)+funchess(3,3)
    write(ifileid,*)
    write(ifileid,*) "Hessian matrix:"
    write(ifileid,"(3E18.10)") funchess
    call diagmat(funchess,eigvecmat,eigval,300,1D-12)
    write(ifileid,"(' Eigenvalues of Hessian:',3E18.10)") eigval(1:3)
    write(ifileid,*) "Eigenvectors(columns) of Hessian:"
    write(ifileid,"(3E18.10)") ((eigvecmat(i,j),j=1,3),i=1,3)
    write(ifileid,"(' Determinant of Hessian:',D18.10)") detmat(funchess)
    if (ifuncsel==1) then !Output ellipicity for rho
        eigmax=maxval(eigval) !At bcp, will be the most positive 
        eigmin=minval(eigval) !At bcp, will be the most negative
        do itmp=1,3
            tmpval=eigval(itmp)
            if (tmpval/=eigmax.and.tmpval/=eigmin) eigmed=tmpval !At bcp, will be the second most negative
        end do
        write(ifileid,"(a,f12.6)") " Ellipticity of electron density:",eigmin/eigmed-1
        write(ifileid,"(a,f12.6)") " eta index:",abs(eigmin)/eigmax
    end if
!     diffstep=1D-5

else !Only loaded structure, use YWT promolecule density
    if (ifiletype==4) then
        write(ifileid,"(' ESP from atomic charges:',E18.10)") nucesp(inx,iny,inz)
    else
        write(ifileid,"(' ESP from nuclear charges:',E18.10)") nucesp(inx,iny,inz)
    end if
    write(ifileid,*)
    call calchessmat_prodens(inx,iny,inz,elerho,elegrad,elehess)
    write(ifileid,"(a)") " Note: The loaded file does not contain wavefunction information, below results are evaluated from promolecule density"
    write(ifileid,"(' Density of electrons:',E18.10)") elerho
    write(ifileid,"(' Reduced density gradient:',E18.10)") RDGprodens(inx,iny,inz)
    write(ifileid,"(' Sign(lambda2)*rho:',E18.10)") signlambda2rho_prodens(inx,iny,inz)
    write(ifileid,"(' User defined real space function:',E18.10)") userfunc(inx,iny,inz)
    write(ifileid,*)
    write(ifileid,*) "Components of gradient in x/y/z are:"
    write(ifileid,"(3E18.10)") elegrad(1),elegrad(2),elegrad(3)
    write(ifileid,"(' Norm of gradient is:',E18.10)") dsqrt(sum(elegrad**2))
    write(ifileid,*)
    write(ifileid,*) "Components of Laplacian in x/y/z are:"
    write(ifileid,"(3E18.10)") elehess(1,1),elehess(2,2),elehess(3,3)
    write(ifileid,"(' Total:',E18.10)") elehess(1,1)+elehess(2,2)+elehess(3,3)
    write(ifileid,*)
    write(ifileid,*) "Hessian matrix:"
    write(ifileid,"(3E18.10)") elehess
    call diagmat(elehess,eigvecmat,eigval,300,1D-12)
    write(ifileid,"(' Eigenvalues of Hessian:',3E18.10)") eigval(1:3)
    write(ifileid,*) "Eigenvectors(columns) of Hessian:"
    write(ifileid,"(3E18.10)") ((eigvecmat(i,j),j=1,3),i=1,3)
end if
end subroutine


!!------------------ Set up grid setting
!If ienableloadextpt==1, then show the option used to load external points
!If igridsel==100, that means user didn't set up grid here but choose to load a set of point coordinates from external plain text file
subroutine setgrid(ienableloadextpt,igridsel)
use defvar
implicit real*8 (a-h,o-z)
real*8 molxlen,molylen,molzlen,tmpx,tmpy,tmpz
character*200 cubefilename,pointfilename
character c80*80
integer ienableloadextpt
logical filealive
ntotlow=125000
ntotmed=512000
ntothigh=1728000
do while(.true.)
    write(*,*) "Please select a method to set up grid"
    write(*,"(a,f10.6,a)") " -10 Set extension distance of grid range for mode 1~4, current:",aug3D," Bohr"
    write(*,*) "1 Low quality grid   , covering whole system, about 125000 points in total"
    write(*,*) "2 Medium quality grid, covering whole system, about 512000 points in total"
    write(*,*) "3 High quality grid  , covering whole system, about 1728000 points in total"
    write(*,*) "4 Input the number of points or grid spacing in X,Y,Z, covering whole system"
    write(*,*) "5 Input original point, translation vector and the number of points"
    write(*,*) "6 Input center coordinate, number of points and extension distance"
    write(*,*) "7 The same as 6, but input two atoms, the midpoint will be defined as center"
    write(*,*) "8 Use grid setting of another cube file"
    if (ienableloadextpt==1) write(*,*) "100 Load a set of points from external file"
    read(*,*) igridsel
    if ((igridsel==1.or.igridsel==2.or.igridsel==3).and.(isys==2.or.isys==3)) write(*,"(a)") "Note: Don't view isosurface after generating the grid, the graph will be total incorrect"
    if (igridsel/=-10) exit
    write(*,*) "Input extension distance (Bohr) e.g. 6.5"
    read(*,*) aug3D
end do

if (igridsel==100) then !Load points rather than set up grid
    write(*,*) "Input the path of the file containing points, e.g. c:\ltwd.txt"
    write(*,*) "Note: See program manual for the format of the file"
    do while(.true.)
        read(*,"(a)") pointfilename
        inquire(file=pointfilename,exist=filealive)
        if (filealive) then
            open(10,file=pointfilename,status="old")
            read(10,*) numextpt
            write(*,"(a,i10,a)") ' There are',numextpt,' points'
            if (allocated(extpt)) deallocate(extpt)
            allocate(extpt(numextpt,4))
            do itmp=1,numextpt
                read(10,*) extpt(itmp,1:3)
            end do
            close(10)
            exit
        else
            write(*,*) "Error: File cannot be found, input again"
        end if
    end do
    write(*,*) "Please wait..."
else
    molxlen=(maxval(a%x)-minval(a%x))+2*aug3D
    molylen=(maxval(a%y)-minval(a%y))+2*aug3D
    molzlen=(maxval(a%z)-minval(a%z))+2*aug3D
    if (molxlen==0.0D0) molxlen=2.0D0 !Avoid catastrophe when aug3D=0 and system is plane
    if (molylen==0.0D0) molylen=2.0D0
    if (molzlen==0.0D0) molzlen=2.0D0
    if (igridsel==1.or.igridsel==2.or.igridsel==3) then
        if (igridsel==1) dx=(molxlen*molylen*molzlen/dfloat(ntotlow))**(1.0D0/3.0D0)
        if (igridsel==2) dx=(molxlen*molylen*molzlen/dfloat(ntotmed))**(1.0D0/3.0D0)
        if (igridsel==3) dx=(molxlen*molylen*molzlen/dfloat(ntothigh))**(1.0D0/3.0D0)
        dy=dx
        dz=dx
        nx=nint(molxlen/dx)+1
        ny=nint(molylen/dy)+1
        nz=nint(molzlen/dz)+1
        orgx=minval(a%x)-aug3D
        orgy=minval(a%y)-aug3D
        orgz=minval(a%z)-aug3D
    else if (igridsel==4) then
        write(*,*) "Input the number of grid points in X,Y,Z direction, e.g. 139,59,80"
        write(*,"(a)") " or input the grid spacing (bohr) in X,Y,Z direction, e.g. 0.05,0.08,0.08  (if only input one value, it will be used for all directions)"
        read(*,"(a)") c80
        if (index(c80,'.')/=0) then
            if (index(c80,',')/=0) then
                read(c80,*) dx,dy,dz
            else
                read(c80,*) tmp
                dx=tmp
                dy=tmp
                dz=tmp
            end if
            nx=molxlen/dx+1
            ny=molylen/dy+1
            nz=molzlen/dz+1
        else
            read(c80,*) nx,ny,nz
            dx=molxlen/(nx-1)
            dy=molylen/(ny-1)
            dz=molzlen/(nz-1)
        end if
        orgx=minval(a%x)-aug3D
        orgy=minval(a%y)-aug3D
        orgz=minval(a%z)-aug3D
    else if (igridsel==5) then
        write(*,*) "Input X,Y,Z coordinate of original point (Bohr) e.g. 0.1,4,-1"
        read(*,*) orgx,orgy,orgz
        write(*,*) "Input X,Y,Z component of translation vector (Bohr) e.g. 0.1,0.1,0.15"
        read(*,*) dx,dy,dz
        write(*,*) "Input the number of points in X,Y,Z direction e.g. 139,59,80"
        read(*,*) nx,ny,nz
    else if (igridsel==6.or.igridsel==7) then
        if (igridsel==6) then
            write(*,*) "Input X,Y,Z coordinate of center (Angstrom)"
            read(*,*) cenx,ceny,cenz
            cenx=cenx/b2a
            ceny=ceny/b2a
            cenz=cenz/b2a
        else if (igridsel==7) then
            write(*,*) "Input index of the two atoms e.g. 2,5"
            write(*,*) "If the two indices are identical, box center will be placed at the nucleus"
            read(*,*) indatm1,indatm2
            cenx=(a(indatm1)%x+a(indatm2)%x)/2.0D0
            ceny=(a(indatm1)%y+a(indatm2)%y)/2.0D0
            cenz=(a(indatm1)%z+a(indatm2)%z)/2.0D0
        end if
        write(*,*) "Input the number of points in X,Y,Z direction e.g. 40,40,25"
        read(*,*) nx,ny,nz
        write(*,*) "Input the extended distance in X,Y,Z direction (Bohr) e.g. 4.0,4.0,6.5"
        read(*,*) aug3Dx,aug3Dy,aug3Dz
        orgx=cenx-aug3Dx
        orgy=ceny-aug3Dy
        orgz=cenz-aug3Dz
        dx=aug3Dx*2.0D0/(nx-1)
        dy=aug3Dy*2.0D0/(ny-1)
        dz=aug3Dz*2.0D0/(nz-1)
    else if (igridsel==8) then
        write(*,*) "Input filename of a cube file"
        do while(.true.)
            read(*,"(a)") cubefilename
            inquire(file=cubefilename,exist=filealive)
            if (filealive) then
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
    end if
    endx=orgx+dx*(nx-1)
    endy=orgy+dy*(ny-1)
    endz=orgz+dz*(nz-1)
    write(*,"(' Coordinate of origin in X,Y,Z is',3f12.6,' Bohr')") orgx,orgy,orgz
    write(*,"(' Coordinate of end point in X,Y,Z is',3f12.6,' Bohr')") endx,endy,endz
    write(*,"(' Grid spacing in X,Y,Z is',3f12.6,' Bohr')") dx,dy,dz
    write(*,"(' The number of points in X,Y,Z is',3i5,'   Total:',i12)") nx,ny,nz,nx*ny*nz
end if
end subroutine


!!!------------------------- Delete virtual orbitals higher than LUMO+10 for HF/DFT wavefunctions
subroutine delvirorb(infomode)
use defvar
implicit real*8 (a-h,o-z)
integer :: infomode,nvirsave=10 !Lowest nvirsave virtual orbitals will be reserved
if (ifiletype/=1.and.ifiletype/=9) return !Only works for fch and molden
if (idelvirorb==0) return
if (iuserfunc==24) return !linear response kernel require all orbital information
if (imodwfn==1) return !Don't make things more complicated!
!This routine doesn't work for post-HF instances
if (wfntype==0.or.wfntype==2) then !RHF, ROHF
    if (nmo<=naelec+nvirsave) return
    nmo=naelec+10 !Simply shield those virtual orbitals
else if (wfntype==1) then !Perserve up to LUMO+10 for alpha, and identical number of orbitals for beta
    if (nmo/2<=naelec+nvirsave) return !naelec is always >= nbelec
    nperserve=naelec+nvirsave
    !Cobasa and Cobasb are needn't to be modified
    co(naelec+11:naelec+nvirsave+nperserve,:)=co(nmo/2+1:nmo/2+nperserve,:)
    MOene(naelec+nvirsave+1:naelec+nvirsave+nperserve)=MOene(nmo/2+1:nmo/2+nperserve)
    MOocc(naelec+nvirsave+1:naelec+nvirsave+nperserve)=MOocc(nmo/2+1:nmo/2+nperserve)
    MOtype(naelec+nvirsave+1:naelec+nvirsave+nperserve)=MOtype(nmo/2+1:nmo/2+nperserve)
    nmo=2*nperserve
end if
imodwfn=1 !Will not call this routine again
if (infomode==1) then
    write(*,"(a)") " Note: Virtual orbitals higher than LUMO+10 have been discarded for saving computational time"
    write(*,*)
end if
end subroutine


!!!-------- imode=1: Convert unit of grid/plane parameters from Bohr to Angstrom. =2: Convert them back
subroutine convgridlenunit(imode)
use defvar
implicit none
integer imode
real*8 scll
if (imode==1) scll=b2a
if (imode==2) scll=1/b2a
orgx=orgx*scll
orgy=orgy*scll
orgz=orgz*scll
orgx2D=orgx2D*scll
orgy2D=orgy2D*scll
orgz2D=orgz2D*scll
endx=endx*scll
endy=endy*scll
endz=endz*scll
dx=dx*scll
dy=dy*scll
dz=dz*scll
v1x=v1x*scll
v1y=v1y*scll
v1z=v1z*scll
v2x=v2x*scll
v2y=v2y*scll
v2z=v2z*scll
a1x=a1x*scll
a1y=a1y*scll
a1z=a1z*scll
a2x=a2x*scll
a2y=a2y*scll
a2z=a2z*scll
a3x=a3x*scll
a3y=a3y*scll
a3z=a3z*scll
d1=d1*scll
d2=d2*scll
end subroutine


!!-------- Deallocate all arrays about wavefunction except that with _org suffix, to avoid problem when load another file
subroutine dealloall
use defvar
implicit real*8 (a-h,o-z)
if (allocated(a)) deallocate(a)
if (allocated(b)) deallocate(b)
if (allocated(CO)) deallocate(CO)
if (allocated(MOocc)) deallocate(MOocc)
if (allocated(MOsym)) deallocate(MOsym)
if (allocated(MOene)) deallocate(MOene)
if (allocated(MOtype)) deallocate(MOtype)
if (allocated(b_EDF)) deallocate(CO_EDF,b_EDF)
!Loaded file contains basis information
if (allocated(shtype)) deallocate(shtype,shcen,shcon,primshexp,primshcoeff,&
basshell,bascen,bastype,basstart,basend,primstart,primend,primconnorm)
if (allocated(CObasa)) deallocate(CObasa)
if (allocated(CObasb)) deallocate(CObasb)
if (allocated(Ptot)) deallocate(Ptot)
if (allocated(Palpha)) deallocate(Palpha)
if (allocated(Pbeta)) deallocate(Pbeta)
if (allocated(Sbas)) deallocate(Sbas)
if (allocated(Dbas)) deallocate(Dbas)
end subroutine





!!------- Generate atomic/fragmental Hirshfeld weight and store it to planemat, calculate free-atom/fragmental density and store it to planemattmp
!The atoms in the fragment is inputted as "selatm" array, nselatm is the number of its elements
!if itype=1, use atomic wavefunction to calculate Hirshfeld weight, and setpromol must have been invoked; if =2, use built-in atomic density to generate it
subroutine genHirshplanewei(selatm,nselatm,itype)
use defvar
use function
implicit real*8 (a-h,o-z)
integer selatm(nselatm),nselatm,itype
if (allocated(planemat)) deallocate(planemat)
if (allocated(planemattmp)) deallocate(planemattmp)
allocate(planemat(ngridnum1,ngridnum2),planemattmp(ngridnum1,ngridnum2))
planemat=0D0
planemattmp=0D0
do iatm=1,ncenter_org !Calc free atomic density of each atom, get promolecular density and Hirshfeld weight of present atom
    iyes=0
    if (any(selatm==iatm)) iyes=1
    if (itype==1) then
        call dealloall
        call readwfn(custommapname(iatm),1)
    end if
!$OMP PARALLEL DO private(i,j,rnowx,rnowy,rnowz,tmpval) shared(planemat) schedule(dynamic) NUM_THREADS( nthreads  )
    do i=1,ngridnum1 !First calculate promolecular density and store it to planemat
        do j=1,ngridnum2
            rnowx=orgx2D+(i-1)*v1x+(j-1)*v2x
            rnowy=orgy2D+(i-1)*v1y+(j-1)*v2y
            rnowz=orgz2D+(i-1)*v1z+(j-1)*v2z
            if (itype==1) then
                tmpval=fdens(rnowx,rnowy,rnowz)
            else
                tmpval=calcatmdens(iatm,rnowx,rnowy,rnowz,0)
            end if
            planemat(i,j)=planemat(i,j)+tmpval
            if (iyes==1) planemattmp(i,j)=planemattmp(i,j)+tmpval
        end do
    end do
!$OMP END PARALLEL DO
end do
if (itype==1) then
    call dealloall
    call readinfile(firstfilename,1) !Retrieve the first loaded file(whole molecule)
end if

do i=1,ngridnum1 !Calculate Hirshfeld weighting function
    do j=1,ngridnum2
        if (planemat(i,j)/=0D0) then
            planemat(i,j)=planemattmp(i,j)/planemat(i,j)
        else
            planemat(i,j)=0D0
        end if
    end do
end do
end subroutine

!!------- Calculate some quantities involved in shubin's project in a plane
!itype=1: Calculate the sum of atomic relative Shannon entropy (namely total relative Shannon entropy)
!itype=2: Calculate the sum of x=[rhoA-rho0A]/rhoA
!itype=3: Calculate the difference between total relative Shannon entropy and deformation density
subroutine genentroplane(itype)
use defvar
use function
implicit real*8 (a-h,o-z)
integer itype
real*8 planeprodens(ngridnum1,ngridnum2),planedens(ngridnum1,ngridnum2)
if (allocated(planemat)) deallocate(planemat)
allocate(planemat(ngridnum1,ngridnum2))
planeprodens=0D0
planemat=0D0
!Calculate molecular density in the plane and store it to planedens
!$OMP PARALLEL DO private(i,j,rnowx,rnowy,rnowz) shared(planedens) schedule(dynamic) NUM_THREADS( nthreads  )
do i=1,ngridnum1
    do j=1,ngridnum2
        rnowx=orgx2D+(i-1)*v1x+(j-1)*v2x
        rnowy=orgy2D+(i-1)*v1y+(j-1)*v2y
        rnowz=orgz2D+(i-1)*v1z+(j-1)*v2z
        planedens(i,j)=fdens(rnowx,rnowy,rnowz)
    end do
end do
!$OMP END PARALLEL DO
do jatm=1,ncenter_org !Calculate promolecular density in the plane and store it to planeprodens
    call dealloall
    call readwfn(custommapname(jatm),1)
!$OMP PARALLEL DO private(i,j,rnowx,rnowy,rnowz) shared(planeprodens) schedule(dynamic) NUM_THREADS( nthreads  )
    do i=1,ngridnum1
        do j=1,ngridnum2
            rnowx=orgx2D+(i-1)*v1x+(j-1)*v2x
            rnowy=orgy2D+(i-1)*v1y+(j-1)*v2y
            rnowz=orgz2D+(i-1)*v1z+(j-1)*v2z
            planeprodens(i,j)=planeprodens(i,j)+fdens(rnowx,rnowy,rnowz)
        end do
    end do
!$OMP END PARALLEL DO
end do
!Calculate Hirshfeld weight, relative Shannon entropy and x=[rhoA-rho0A]/rhoA for each atom in the plane and accumulate them to planemat
do jatm=1,ncenter_org !Cycle each atom, calculate its contribution in the plane
    call dealloall
    call readwfn(custommapname(jatm),1)
!$OMP PARALLEL DO private(i,j,rnowx,rnowy,rnowz,rho0A,rhoA,tmpval) shared(planemat) schedule(dynamic) NUM_THREADS( nthreads  )
    do i=1,ngridnum1
        do j=1,ngridnum2
            rnowx=orgx2D+(i-1)*v1x+(j-1)*v2x
            rnowy=orgy2D+(i-1)*v1y+(j-1)*v2y
            rnowz=orgz2D+(i-1)*v1z+(j-1)*v2z
            rho0A=fdens(rnowx,rnowy,rnowz)
            rhoA=planedens(i,j)*rho0A/planeprodens(i,j)
            if (itype==1.or.itype==3) tmpval=rhoA*log(rhoA/rho0A) !Relative Shannon entropy
            if (itype==2) tmpval=(rhoA-rho0A)/rhoA !x=[rhoA-rho0A]/rhoA
            planemat(i,j)=planemat(i,j)+tmpval
        end do
    end do
!$OMP END PARALLEL DO
end do
call dealloall
call readinfile(firstfilename,1) !Retrieve the first loaded file(whole molecule)
if (itype==3) planemat=planemat-(planedens-planeprodens) !Diff between total relative Shannon entropy and deformation density
end subroutine



!!----- Generate atomic Hirshfeld weight and store it to cubmat
!The atoms in the fragment is inputted as "selatm" array, nselatm is the number of its elements
!if itype=1, use atomic wavefunction to calculate Hirshfeld weight, and setpromol must have been invoked; if =2, use built-in atomic density to generate it
subroutine genHirshcubewei(selatm,nselatm,itype)
use defvar
use function
implicit real*8 (a-h,o-z)
integer selatm(nselatm),nselatm,itype
if (allocated(cubmat)) deallocate(cubmat)
if (allocated(cubmattmp)) deallocate(cubmattmp)
allocate(cubmat(nx,ny,nz),cubmattmp(nx,ny,nz))
cubmat=0D0
cubmattmp=0D0
do iatm=1,ncenter_org
    write(*,"(' Finished',i6,'  /',i6)") iatm,ncenter_org
    if (itype==1) then
        call dealloall
        call readwfn(custommapname(iatm),1)
    end if
!$OMP PARALLEL DO SHARED(cubmat,cubmattmp,ifinish) PRIVATE(i,j,k,tmpx,tmpy,tmpz,tmpval) schedule(dynamic) NUM_THREADS( nthreads  )
    do k=1,nz !First calculate promolecular density and store it to cubmat
        tmpz=orgz+(k-1)*dz
        do j=1,ny
            tmpy=orgy+(j-1)*dy
            do i=1,nx
                tmpx=orgx+(i-1)*dx
                if (itype==1) then
                    tmpval=fdens(tmpx,tmpy,tmpz)
                else
                    tmpval=calcatmdens(iatm,tmpx,tmpy,tmpz,0)
                end if
                cubmat(i,j,k)=cubmat(i,j,k)+tmpval
                if (any(selatm==iatm)) cubmattmp(i,j,k)=cubmattmp(i,j,k)+tmpval
            end do
        end do
    end do
!$OMP END PARALLEL DO
end do
if (itype==1) then
    call dealloall
    call readinfile(firstfilename,1) !Retrieve the first loaded file(whole molecule)
end if

do k=1,nz !Calculate Hirshfeld weighting function
    do j=1,ny
        do i=1,nx
            if (cubmat(i,j,k)/=0D0) then
                cubmat(i,j,k)=cubmattmp(i,j,k)/cubmat(i,j,k)
            else
                cubmat(i,j,k)=0D0
            end if
        end do
    end do
end do
end subroutine




!!----------- Output spherically averaged atomic radial density, then used for generating promolecular density for heavy atoms
subroutine sphatmraddens
use defvar
use function
implicit real*8 (a-h,o-z)
real*8,allocatable :: potx(:),poty(:),potz(:),potw(:),radpos(:),sphavgval(:)
truncrho=1D-8
rlow=0D0
rhigh=12
nsphpt=2030
nradpt=200 !Totally 200 radial points, but the number of point is truncated at truncrho
allocate(potx(nsphpt),poty(nsphpt),potz(nsphpt),potw(nsphpt),radpos(nradpt),sphavgval(nradpt))
call Lebedevgen(nsphpt,potx,poty,potz,potw)
ifinish=0
iprogstp=20
iprogcrit=iprogstp
write(*,*) "Calculating..."
!$OMP PARALLEL DO SHARED(sphavgval,radpos,ifinish,iprogcrit) PRIVATE(irad,radx,radr,isph,rnowx,rnowy,rnowz,tmpval) schedule(dynamic) NUM_THREADS( nthreads  )
do irad=1,nradpt
    radx=cos(irad*pi/(nradpt+1))
    radr=(1+radx)/(1-radx) !Becke transform
    radpos(irad)=radr
    tmpval=0
    do isph=1,nsphpt
        rnowx=potx(isph)*radr
        rnowy=poty(isph)*radr
        rnowz=potz(isph)*radr
        tmpval=tmpval+fdens(rnowx,rnowy,rnowz)*potw(isph)
    end do
    sphavgval(irad)=tmpval !Spherically average density
    ifinish=ifinish+1
    if (ifinish==iprogcrit) then
        write(*,"(' Finished:',i6,'  /',i6)") ifinish,nradpt
        iprogcrit=iprogcrit+iprogstp
    end if
end do
!$OMP END PARALLEL DO
open(10,file="sphavgval.txt",status="replace")
itmp=0
do irad=nradpt,1,-1
    if (sphavgval(irad)>truncrho) itmp=itmp+1
end do
write(10,"(a,i3,a)") "else if (iele==",a(1)%index,") then  !"
write(10,"('    npt=',i5)") itmp
itmp=0
do irad=nradpt,1,-1
    if (sphavgval(irad)>truncrho) then
        itmp=itmp+1
        write(10,"('    rhoarr(',i3,')=',f25.10,'D0')") itmp,sphavgval(irad)
    end if
end do
close(10)
write(*,*) "The result has been output to sphavgval.txt in current folder"
write(*,*) "The second column is radial distance (Bohr), the third column is value"
end subroutine




!!--- Generate single-center integration grid for Becke's integration. Not adapted according to element. Return iradcut and gridatm
subroutine gen1cintgrid(gridatm,iradcut)
use defvar
implicit real*8 (a-h,o-z)
integer iradcut
real*8 potx(sphpot),poty(sphpot),potz(sphpot),potw(sphpot)
type(content) gridatm(radpot*sphpot)
call Lebedevgen(sphpot,potx,poty,potz,potw)
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
end subroutine
!!--- Generate Becke weight for a batch of points around iatm, sharpness parameter=3
!!--- Input: iatm, iradcut, gridatm   Return: beckeweigrid
subroutine gen1cbeckewei(iatm,iradcut,gridatm,beckeweigrid)
use defvar
implicit real*8 (a-h,o-z)
integer iatm,iradcut
real*8 beckeweigrid(radpot*sphpot),smat(ncenter,ncenter),Pvec(ncenter)
type(content) gridatm(radpot*sphpot)
!$OMP parallel do shared(beckeweigrid) private(i,rnowx,rnowy,rnowz,smat,ii,ri,jj,rj,rmiu,chi,uij,aij,tmps,Pvec) num_threads( nthreads  ) schedule(dynamic)
do i=1+iradcut*sphpot,radpot*sphpot
    smat=1D0
    rnowx=gridatm(i)%x
    rnowy=gridatm(i)%y
    rnowz=gridatm(i)%z
    do ii=1,ncenter
        ri=dsqrt( (rnowx-a(ii)%x)**2+(rnowy-a(ii)%y)**2+(rnowz-a(ii)%z)**2 )
        do jj=1,ncenter
            if (ii==jj) cycle
            rj=dsqrt( (rnowx-a(jj)%x)**2+(rnowy-a(jj)%y)**2+(rnowz-a(jj)%z)**2 )
            rmiu=(ri-rj)/distmat(ii,jj)
             !Adjust for heteronuclear
            chi=covr_tianlu(a(ii)%index)/covr_tianlu(a(jj)%index)
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
    Pvec=1D0
    do ii=1,ncenter
        Pvec=Pvec*smat(:,ii)
    end do
    beckeweigrid(i)=Pvec(iatm)/sum(Pvec)
end do
!$OMP end parallel do
end subroutine

!!-------- Randomly generate the name of Sobereva's lover
subroutine mylover(outname)
integer,parameter :: nlovers=43
character*80 lovername(nlovers),outname
CALL RANDOM_SEED()
CALL RANDOM_NUMBER(tmp)
lovername(1)="K-ON\Mio_Akiyama"
lovername(2)="K-ON\Azusa_Nakano"
lovername(3)="EVA\Rei_Ayanami"
lovername(4)="Ore_no_Imoto\Black_Cat"
lovername(5)="Touhou_project\Ran_Yakumo"
lovername(6)="Haiyore!Nyaruko-san\Nyaruko"
lovername(7)="Bodacious_Space_Pirates\Kurihara_Chiaki"
lovername(8)="Otoboku\Mariya_Mikado"
lovername(9)="Amagami\Miya_Tachibana"
lovername(10)="Shakugan_no_Shana\Shana"
lovername(11)="Yuru_Yuri\Akari_Akaza"
lovername(12)="Natsuiro_Kiseki\Yuka_Hanaki"
lovername(13)="Love_Live!\Nico_Yazawa"
lovername(14)="Vocaloid\Miku_Hatsune"
lovername(15)="iDOLM@STER\Makoto_Kikuchi"
lovername(16)="Last_Exile\Dio_Eraclea"
lovername(17)="NHK_ni_Youkoso!\Misaki_Nakahara"
lovername(18)="Rio_Rainbow_Gate\Rio_Rollins"
lovername(19)="Blood-C\Saya_Kisaragi"
lovername(20)="Mahou_Shoujo_Madoka-Magica\Homura_Akemi"
lovername(21)="Saki\Hisa_Takei"
lovername(22)="Strawberry_Panic\Chikaru_Minamoto"
lovername(23)="Najica\Najica_Hiiragi"
lovername(24)="Blue_Drop\Hagino_Senkouji"
lovername(25)="Fate_Zero\Saber"
lovername(26)="Baka_to_Test_to_Shoukanjuu\Hideyoshi_Kinoshita"
lovername(27)="Watamote\Tomoko_Kuroki"
lovername(28)="Genshiken_Nidaime\Kenjirou_Hato"
lovername(29)="The_World_God_Only_Knows\Chihiro_Kosaka"
lovername(30)="Kan_Colle\Shimakaze"
lovername(31)="Wake_Up,Girls!\Miyu_Okamoto"
lovername(32)="Gokukoku\Kazumi_Schlierenzauer"
lovername(33)="Love_Live!\Nozomi_Tojo"
lovername(34)="Tokimeki_Memorial\Yuina_Himoo"
lovername(35)="MADLAX\MADLAX"
lovername(36)="Gun_Gale_Online\Kirito"
lovername(37)="Denkigai_No_Honyasan\Sennsei"
lovername(38)="Kan_Colle\Kongou"
lovername(39)="Plastic_Memories\Aira"
lovername(40)="MaiMengJun\QinXue_Chen"
lovername(41)="Sakurako-san_no_Ashimoto_ni_wa_Shitai_ga_Umatteiru\Sakurako"
lovername(42)="Hibike!_Euphonium\Reina_Kousaka"
lovername(43)="Planetarian\Yumemi_Hoshino"
!Dear Kanan,
!
!You are the only one I deeply love forever in the real world,
!although you can't be with me, and I am even unable to know your name and touch your finger.
!I believe I will never love anyone else in the rest of my life.
!
!I love your brilliant dance, your kawaii smile, your lovely double ponytail, and especially, your extremely pure and beautiful heart.
!
!                     ----- The author of Multiwfn, Tian Lu, 2015-May-19
outname=lovername(ceiling(tmp*nlovers))
end subroutine


!!------------ Load parameters in settings.ini when boot up multiwfn
subroutine loadsetting
use defvar
use util
implicit real*8 (a-h,o-z)
character c80tmp*80,c200tmp*200,settingpath*80
!Set default color of atomic spheres
atm3Dclr(:,1)=0.85D0
atm3Dclr(:,2)=0.6D0
atm3Dclr(:,3)=0.5D0
atm3Dclr(0,:)=(/0.0D0,  0.5D0,  0.6D0 /) !Bq
atm3Dclr(1,:)=(/0.95D0, 0.95D0, 0.95D0/) !H
atm3Dclr(5,:)=(/0.95D0, 0.7D0,  0.7D0 /) !B
atm3Dclr(6,:)=(/0.85D0, 0.85D0, 0.55D0/) !C
atm3Dclr(7,:)=(/0.5D0,  0.5D0,  1.0D0 /) !N
atm3Dclr(8,:)=(/1.0D0,  0.2D0,  0.2D0 /) !O
atm3Dclr(9,:)=(/0.6D0,  0.9D0,  0.9D0 /) !F
atm3Dclr(15,:)=(/0.9D0, 0.4D0,  0.0D0 /) !P
atm3Dclr(16,:)=(/0.9D0, 0.7D0,  0.1D0 /) !S
atm3Dclr(17,:)=(/0.1D0, 0.9D0,  0.1D0 /) !Cl

inquire(file="settings.ini",exist=alive)
if (alive.eqv..true.) then
    settingpath="settings.ini"
else if (alive.eqv..false.) then
    call getenv("Multiwfnpath",c80tmp)
    if (isys==1) then
        settingpath=trim(c80tmp)//"\settings.ini"
    else if (isys==2.or.isys==3) then
        settingpath=trim(c80tmp)//"/settings.ini"
    end if
    inquire(file=settingpath,exist=alive)
    if (alive.eqv..false.) then
        write(*,"(a)") " Warning: ""settings.ini"" was found neither in current folder nor in the path defined by ""Multiwfnpath"" &
        environment variable. Now using default settings instead"
        write(*,*)
        return
    end if
end if

open(20,file=settingpath,status="old")
! Below are the parameters can affect calculation results
call loclabel(20,'iuserfunc=')
read(20,*) c80tmp,iuserfunc
call loclabel(20,'refxyz=')
read(20,*) c80tmp,refx,refy,refz
call loclabel(20,'iDFTxcsel=')
read(20,*) c80tmp,iDFTxcsel
call loclabel(20,'paircorrtype=')
read(20,*) c80tmp,paircorrtype
call loclabel(20,'pairfunctype=')
read(20,*) c80tmp,pairfunctype
call loclabel(20,'iautointgrid=')
read(20,*) c80tmp,iautointgrid
call loclabel(20,'radpot=')
read(20,*) c80tmp,radpot
call loclabel(20,'sphpot=')
read(20,*) c80tmp,sphpot
call loclabel(20,'radcut=')
read(20,*) c80tmp,radcut
call loclabel(20,'expcutoff=')
read(20,*) c80tmp,expcutoff
call loclabel(20,'espprecutoff=')
read(20,*) c80tmp,espprecutoff
call loclabel(20,'RDG_maxrho=')
read(20,*) c80tmp,RDG_maxrho
call loclabel(20,'RDGprodens_maxrho=')
read(20,*) c80tmp,RDGprodens_maxrho
call loclabel(20,'ELF_addminimal=')
read(20,*) c80tmp,ELF_addminimal
call loclabel(20,'ELFLOL_type=')
read(20,*) c80tmp,ELFLOL_type
call loclabel(20,'ELFLOL_cut=')
read(20,*) c80tmp,ELFLOL_cut
call loclabel(20,'iALIEdecomp=')
read(20,*) c80tmp,iALIEdecomp
call loclabel(20,'srcfuncmode=')
read(20,*) c80tmp,srcfuncmode
call loclabel(20,'atomdenscut=')
read(20,*) c80tmp,atomdenscut
call loclabel(20,'aug1D=')
read(20,*) c80tmp,aug1D
call loclabel(20,'aug2D=')
read(20,*) c80tmp,aug2D
call loclabel(20,'aug3D=')
read(20,*) c80tmp,aug3D
call loclabel(20,'num1Dpoints=')
read(20,*) c80tmp,num1Dpoints
call loclabel(20,'nprevorbgrid=')
read(20,*) c80tmp,nprevorbgrid
call loclabel(20,'bndordthres=')
read(20,*) c80tmp,bndordthres
call loclabel(20,'compthres=')
read(20,*) c80tmp,compthres
call loclabel(20,'compthresCDA=')
read(20,*) c80tmp,compthresCDA
call loclabel(20,'ispheratm=')
read(20,*) c80tmp,ispheratm
call loclabel(20,'laplfac=')
read(20,*) c80tmp,laplfac
call loclabel(20,'ipolarpara=')
read(20,*) c80tmp,ipolarpara
call loclabel(20,'ADCtransfer=')
read(20,*) c80tmp,ADCtransfer
call loclabel(20,'SpherIVgroup=')
read(20,*) c80tmp,SpherIVgroup
call loclabel(20,'MCvolmethod=')
read(20,*) c80tmp,MCvolmethod
call loclabel(20,'readEDF=')
read(20,*) c80tmp,readEDF
call loclabel(20,'ireadatmEDF=')
read(20,*) c80tmp,ireadatmEDF
call loclabel(20,'idelvirorb=')
read(20,*) c80tmp,idelvirorb
call loclabel(20,'ifchprog=')
read(20,*) c80tmp,ifchprog
call loclabel(20,'ishowptESP=')
read(20,*) c80tmp,ishowptESP
call loclabel(20,'imolsurparmode=')
read(20,*) c80tmp,imolsurparmode
call loclabel(20,'steric_addminimal=')
read(20,*) c80tmp,steric_addminimal
call loclabel(20,'steric_potcutrho=')
read(20,*) c80tmp,steric_potcutrho
call loclabel(20,'steric_potcons=')
read(20,*) c80tmp,steric_potcons
call loclabel(20,'NICSnptlim=')
read(20,*) c80tmp,NICSnptlim
call loclabel(20,'iplaneextdata=')
read(20,*) c80tmp,iplaneextdata
call loclabel(20,'igenP=')
read(20,*) c80tmp,igenP
call loclabel(20,'igenDbas=')
read(20,*) c80tmp,igenDbas
call loclabel(20,'igenMagbas=')
read(20,*) c80tmp,igenMagbas
!Below are the parameters involved in plotting
call loclabel(20,'symbolsize=')
read(20,*) c80tmp,symbolsize
call loclabel(20,'pleatmlabsize=')
read(20,*) c80tmp,pleatmlabsize
call loclabel(20,'disshowlabel=')
read(20,*) c80tmp,disshowlabel
call loclabel(20,'iatom_on_contour_far=')
read(20,*) c80tmp,iatom_on_contour_far
call loclabel(20,'iatmlabtype=')
read(20,*) c80tmp,iatmlabtype
call loclabel(20,'iatmlabtype3D=')
read(20,*) c80tmp,iatmlabtype3D
call loclabel(20,'graphformat=')
read(20,*) c80tmp,graphformat
call loclabel(20,'graph1Dsize=')
read(20,*) c80tmp,graph1Dwidth,graph1Dheight
call loclabel(20,'graph2Dsize=')
read(20,*) c80tmp,graph2Dwidth,graph2Dheight
call loclabel(20,'graph3Dsize=')
read(20,*) c80tmp,graph3Dwidth,graph3Dheight
call loclabel(20,'numdigxyz=')
read(20,*) c80tmp,numdigx,numdigy,numdigz
call loclabel(20,'numdiglinexy=')
read(20,*) c80tmp,numdiglinex,numdigliney
call loclabel(20,'numdigctr=')
read(20,*) c80tmp,numdigctr
call loclabel(20,'fillcoloritpxy=')
read(20,*) c80tmp,fillcoloritpx,fillcoloritpy
call loclabel(20,'inowhiteblack=')
read(20,*) c80tmp,inowhiteblack
call loclabel(20,'isurfstyle=')
read(20,*) c80tmp,isurfstyle
call loclabel(20,'bondRGB=')
read(20,*) c80tmp,bondclrR,bondclrG,bondclrB
call loclabel(20,'atmlabRGB=')
read(20,*) c80tmp,atmlabclrR,atmlabclrG,atmlabclrB
call loclabel(20,'CP_RGB=')
read(20,*) c80tmp,CP3n3RGB,CP3n1RGB,CP3p1RGB,CP3p3RGB
call loclabel(20,'atmcolorfile=') !Set atom 3D color either according to external file or default setting
read(20,*) c80tmp,c200tmp
inquire(file=c200tmp,exist=alive)
if (index(c200tmp,"none")==0.and.alive) then
    write(*,"(' Note: Loading atom color settings from ',a)") trim(c200tmp)
    open(21,file=c200tmp,status="old")
    do iele=0,nelesupp
        read(21,*) inouse,atm3Dclr(iele,:)
    end do
    close(21)
end if
!Below are parameters about system
call loclabel(20,'nthreads=')
read(20,*) c80tmp,nthreads
call loclabel(20,'ompstacksize=')
read(20,*) c80tmp,ompstacksize
call loclabel(20,'gaupath=')
read(20,*) c80tmp,gaupath
call loclabel(20,'isilent=')
read(20,*) c80tmp,isilent
call loclabel(20,'isys=')
read(20,*) c80tmp,isys
call loclabel(20,'imodlayout=')
read(20,*) c80tmp,imodlayout
call loclabel(20,'outmedinfo=')
read(20,*) c80tmp,outmedinfo
call loclabel(20,'iopengl=')
read(20,*) c80tmp,iopengl
call loclabel(20,'iwfntmptype=')
read(20,*) c80tmp,iwfntmptype
call loclabel(20,'ispecial=')
read(20,*) c80tmp,ispecial
!The last opened file name
call loclabel(20,'lastfile=')
read(20,"(10x,a)") lastfile
close(20)
end subroutine
