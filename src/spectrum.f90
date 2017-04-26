subroutine plotspectrum
use defvar
use util
implicit real*8 (a-h,o-z)
real*8,allocatable :: str(:),FWHM(:),linex(:),liney(:),datax(:),indcurve(:,:)
integer,allocatable :: tmparr(:),indband2idx(:),idx2indband(:)
character ctest,c200tmp*200,c80tmp*80
real*8 :: sclfreqfac=1D0
integer :: icurveclr=1,ilineclr=5 !curve color=Red, line color=black
if (allocated(curvex)) deallocate(curvex,curvey,curveytmp) !Note that they are global allocatable arrays
allocate(curvex(num1Dpoints),curvey(num1Dpoints),curveytmp(num1Dpoints))
gauweigh=0.5D0
iusersetY1=0 !If user has set the axes by himself
iusersetY2=0
iusersetX=0
idraw=0
isavepic=0
ishowline=1
ishowgrid=1
iunitliney=1 !Only for IR
shiftx=0D0 !Shift value in X direction
iramantype=1 !=1 Raman activities   =2 Raman intensities
if (index(filename,"tda.dat")/=0) then !sTDA output file
    write(*,*) "Select spectrum type, 3:UV-Vis  4:ECD"
else
    write(*,*) "Select spectrum type, 1:IR  2:Raman  3:UV-Vis  4:ECD  5:VCD"
end if
read(*,*) ispectrum
if (ispectrum==1.or.ispectrum==2.or.ispectrum==5) then
    ibroadfunc=1 !Usually, IR, Raman and VCD use Lorentzian
    iunitx=0 !=0 cm^-1, =1 eV, =2 nm, =3 1000cm^-1
else if (ispectrum==3.or.ispectrum==4) then
    ibroadfunc=2 !Usually, UV-Vis, ECD use Gaussian
    iunitx=2
end if

!For ECD when eV is used, integrating the peak of a unit strength is 1
!For Raman and VCD, integrating the peak of a unit strength is 1
if (ispectrum==2.or.ispectrum==4.or.ispectrum==5) then
    scalecurve=1D0
else if(ispectrum==1) then
!For IR when km/L is used, integrating the peak of a unit strength is 100
    scalecurve=100D0
else if (ispectrum==3) then
!Classically, 1 unit oscillator strength can be broadened to 28700 area(X:eV Y:L/mol/cm)
!Can be derived: 1/(4.32D-9)/8065.5447=28700   Here 1eV=8065.5447/cm, 4.32D-9 comes from Swizard manual (see also Review in C.C. vol.20 p168)
!result is also consistent with Gaussview
    if (iunitx==1.or.iunitx==2) scalecurve=28700
    if (iunitx==3) scalecurve=1D0/4.32D-6
end if

!Note: For nm unit, we still store all data and FWHM in eV, and generate curve as usual in eV. Only at final stage, we scale the curve to get the one in nm
!If we choose 1000cm^-1, we immediately convert all data and FWHM in 1000cm^-1 before generate curve.
!When unit is changed, we reset lower and upper limit to auto, rather than convert them to current unit to avoid problems.

open(10,file=filename,status="old")
if (index(filename,"tda.dat")/=0) then !sTDA output file
    if (ispectrum==4) then
        write(*,*) "Read the rotatory strengths in which representation?"
        write(*,*) "1: Length representation     2: Velocity representation"
        write(*,*) "3: Mixed-form representation (recommended)"
        read(*,*) irotrep
    end if
    call loclabel(10,"DATXY")
    read(10,*)
    numdata=0
    do while(.true.)
        read(10,"(a)",iostat=ierror) c80tmp
        if (c80tmp==" ".or.ierror/=0) exit
        numdata=numdata+1
    end do
    allocate(datax(numdata),str(numdata),FWHM(numdata),linex(3*numdata),liney(3*numdata))
    if (ispectrum==3) FWHM=2D0/3D0
    if (ispectrum==4) FWHM=0.2D0
    call loclabel(10,"DATXY",ifound)
    read(10,*)
    do i=1,numdata
        !Note: The last four columns of tda.dat correspond to f_length, f_velocity, R_length, R_velocity
        read(10,*) inouse,datax(i),fl,fv,Rl,Rv
        if (ispectrum==3) then
            str(i)=fl
        else if (ispectrum==4) then
            if (irotrep==1) str(i)=Rl
            if (irotrep==2) str(i)=Rv
            if (irotrep==3) str(i)=Rv*fl/fv
        end if
    end do
else !Gaussian output file or plain text file
    call loclabel(10,"Gaussian, Inc",igauout)
    rewind(10)
    if (igauout==1) then
        write(*,*) "This is a Gaussian output file, loading..."

        !IR, Raman, VCD
        if (ispectrum==1.or.ispectrum==2.or.ispectrum==5) then !Find how many frequencies in the file
            do while(.true.)
                call loclabel(10,"Frequencies -- ",ifound,0) !HPmodes is also compatible, because in this manner we locate to the tranditional output section
                if (ifound==1) then
                    i1=0
                    i2=0
                    i3=0
                    backspace(10)
                    backspace(10)
                    read(10,*,iostat=ierror) i1,i2,i3
                    if (ierror/=0) then
                        read(10,*,iostat=ierror) i1,i2
                        if (ierror/=0) then
                            read(10,*,iostat=ierror) i1
                        end if
                    end if
                    read(10,*)
                    read(10,*)
                    if (i1==0.or.i2==0.or.i3==0) exit
                else
                    exit
                end if
            end do
            numdata=max(i1,i2,i3)
            rewind(10)

            allocate(datax(numdata),str(numdata),FWHM(numdata),linex(3*numdata),liney(3*numdata))
            FWHM=8D0
            ilackdata=numdata
            inow=1
            do while(.true.)
                if (ilackdata>3) then
                    iread=3
                else
                    iread=ilackdata
                end if
                call loclabel(10,"Frequencies -- ",ifound,0)
                read(10,"(16x)",advance="no")
                if (iread==1) read(10,*) datax(inow)
                if (iread==2) read(10,*) datax(inow),datax(inow+1)
                if (iread==3) read(10,*) datax(inow),datax(inow+1),datax(inow+2)
                if (ispectrum==1) call loclabel(10,"IR Inten    --",ifound,0)
                if (ispectrum==2) call loclabel(10,"Raman Activ --",ifound,0)
                if (ispectrum==5) call loclabel(10,"Rot. str.",ifound,0)            
                read(10,"(16x)",advance="no")
                if (iread==1) read(10,*) str(inow)
                if (iread==2) read(10,*) str(inow),str(inow+1)
                if (iread==3) read(10,*) str(inow),str(inow+1),str(inow+2)
                if (ilackdata<=3) exit
                ilackdata=ilackdata-3
                inow=inow+3
            end do
            
            if (ispectrum==1) then
                call loclabel(10,"Anharmonic Infrared Spectroscopy",ifound,0)
                if (ifound==1) then
                    write(*,"(a)") " Note: Found anharmonic frequencies, if load them instead of harmonic frequencies? (y/n)"
                    read(*,*) ctest
                    if (ctest=='y'.or.ctest=='Y') then
                        call loclabel(10,"Mode(Quanta)",ifound,0)
                        read(10,*)
                        do itmp=1,numdata
                            read(10,*) c200tmp,rnouse,datax(itmp),rnouse,str(itmp)
                        end do
                    end if
                end if
            end if
        
        !UV-Vis, ECD
        else if (ispectrum==3.or.ispectrum==4) then
            !Because this may be an excited state optimization task, we need to determine how many steps are there
            numopt=0
            do while(.true.)
                call loclabel(10,"Excitation energies and oscillator strengths",ifound,0)
                numopt=numopt+ifound
                if (ifound==0) exit
                read(10,*)
            end do
            !Check how many states
            rewind(10)
            do i=1,numopt !Locate to the last time of output
                call loclabel(10,"Excitation energies and oscillator strengths",ifound,0)
                read(10,*)
            end do
            numdata=0
            do while(.true.)
                call loclabel(10," Excited State",ifound,0)
                if (ifound==0) exit
                read(10,*)
                numdata=numdata+1
            end do
            allocate(datax(numdata),str(numdata),FWHM(numdata),linex(3*numdata),liney(3*numdata))
            FWHM=2/3D0
            !Locate to the last time of output
            rewind(10)
            do i=1,numopt
                call loclabel(10,"Excitation energies and oscillator strengths",ifound,0)
                read(10,*)
            end do
            do i=1,numdata !Gaussian output is too flexible to use fixed format to read in
                call loclabel(10," Excited State",ifound,0)
                do while(.true.)
                    read(10,"(a)",advance="no") ctest
                    if (ctest=='-') exit
                end do
                read(10,"(5x)",advance="no")
                read(10,*) datax(i) !Read excitation energy (eV)
                backspace(10)
                do while(.true.)
                    read(10,"(a)",advance="no") ctest
                    if (ctest=='=') exit
                end do
                read(10,*) str(i) !Read oscillator strength
            end do
            if (ispectrum==4) then !Read ECD rotatory strength
                write(*,*) "Read the rotatory strengths in which representation?"
                write(*,*) "1: Length representation     2: Velocity representation (Recommended)"
                read(*,*) itmp
                if (itmp==1) then
                    call loclabel(10,"R(length)",ifound,1)
                    read(10,*)
                    do i=1,numdata
                        read(10,*) inouse,rnouse,rnouse,rnouse,str(i)
                    end do
                else
                    call loclabel(10,"R(velocity)",ifound,1)
                    read(10,*)
                    do i=1,numdata
                        read(10,*) inouse,rnouse,rnouse,rnouse,str(i)
                        do while(.true.)
                            read(10,"(a)") c80tmp
                            !Sometimes Gaussian output some additional info.
                            if (index(c80tmp,"Total R(velocity) tensor")/=0.or.index(c80tmp,"R(velocity) tensor in inp. orien.")/=0) then
                                do itmp=1,4
                                    read(10,*)
                                end do
                            else
                                backspace(10)
                                exit
                            end if
                        end do
                    end do
                end if
            end if
        end if
        
    else !Other files
        iORCAout=0
        call loclabel(10,"O   R   C   A",iORCAout,1) !ORCA output is also supported
        rewind(10)
        if (iORCAout==1) then
            write(*,*) "This is an ORCA output file"
            if (ispectrum==1.or.ispectrum==2) then !IR, Raman
                if (ispectrum==1) call loclabel(10,"IR SPECTRUM")
                if (ispectrum==2) call loclabel(10,"RAMAN SPECTRUM")
                call loclabel(10,"Mode    freq (cm**-1)",ifound,0)
                read(10,*)
                read(10,*)
                numdata=0
                do while(.true.)
                    read(10,"(a)") c80tmp
                    if (c80tmp==" ") exit
                    numdata=numdata+1
                end do
                allocate(datax(numdata),str(numdata),FWHM(numdata),linex(3*numdata),liney(3*numdata))
                if (ispectrum==1) call loclabel(10,"IR SPECTRUM",ifound,1)
                if (ispectrum==2) call loclabel(10,"RAMAN SPECTRUM",ifound,1)
                call loclabel(10,"Mode    freq (cm**-1)",ifound,0)
                read(10,*)
                read(10,*)
                do i=1,numdata
                    read(10,*) c80tmp,datax(i),str(i)
                end do
                FWHM=8D0
            else if (ispectrum==3.or.ispectrum==4) then !UV-Vis, ECD
                call loclabel(10,"Number of roots to be determined")
                read(10,"(50x,i7)") numdata
                allocate(datax(numdata),str(numdata),FWHM(numdata),linex(3*numdata),liney(3*numdata))
                if (ispectrum==3) call loclabel(10,"ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS",ifound,0)
                if (ispectrum==4) call loclabel(10,"CD SPECTRUM",ifound,0)
                read(10,*)
                read(10,*)
                read(10,*)
                read(10,*)
                read(10,*)
                do i=1,numdata
                    read(10,*) rnouse,datax(i),rnouse,str(i)
                end do
                FWHM=2D0/3D0
                datax=datax/8065.5447D0 !Convert from cm-1 to eV
            end if
            
        else !Plain text file
            read(10,*) numdata,inptype
            allocate(datax(numdata),str(numdata),FWHM(numdata),linex(3*numdata),liney(3*numdata))
            if (inptype==1) then !Only x-position and strengths
                if (ispectrum==1.or.ispectrum==2.or.ispectrum==5) FWHM=8D0
                if (ispectrum==3) FWHM=2D0/3D0
                if (ispectrum==4) FWHM=0.2D0
                do i=1,numdata
                    read(10,*) datax(i),str(i)
                end do
            else if (inptype==2) then !also with FWHM
                do i=1,numdata
                    read(10,*) datax(i),str(i),FWHM(i)
                end do
            end if
        end if
    end if
end if    
    
close(10)


do while(.true.)
    write(*,*)
    write(*,*) "-3 Return to main menu"
    write(*,*) "-2 Export transition data to plain text file"
    write(*,*) "-1 Show transition data"
    write(*,*) "0 Plot spectrum"
    write(*,*) "1 Save picture"
    write(*,*) "2 Export X-Y data set of lines and curves to plain text file"
    if (iusersetX==0) write(*,*) "3 Set lower and upper limit of X-axis, current: Auto"
    if (iusersetX==1) write(*,"(a,f12.5,a,f12.5)") " 3 Set lower and upper limit of X-axis, current:",xlow," to",xhigh
    if (iusersetY1==0) write(*,*) "4 Set left Y-axis, current: Auto"
    if (iusersetY1==1) write(*,"(' 4 Set left Y-axis, current low:',f12.3,' high:',f12.3,' step:',f11.3)") orgy1,endy1,stepy1
    if (iusersetY2==0) write(*,*) "5 Set right Y-axis, current: Auto"
    if (iusersetY2==1) write(*,"(' 5 Set right Y-axis, current low:',f10.4,' high:',f10.4,' step:',f10.4)") orgy2,endy2,stepy2
    if (ibroadfunc==1) write(*,*) "6 Select broadening function, current: Lorentzian"
    if (ibroadfunc==2) write(*,*) "6 Select broadening function, current: Gaussian"
    if (ibroadfunc==3) write(*,*) "6 Select broadening function, current: Pseudo-Voigt"
    write(*,"(a,f20.5)") " 7 Set scale ratio for curve, current:",scalecurve
    !For nm unit, FWHM couldn't be defined, since is not linear. Instead, we calculate curve as usually as in eV (use FWHM defined for eV), finally, we scale curve in x to obtain the curve in nm.
    if (maxval(FWHM)==minval(FWHM)) then
        if (iunitx==0) then
            write(*,"(a,f20.5,' cm^-1')") " 8 Input full width at half maximum (FWHM), current:",FWHM(1)
        else if (iunitx==1.or.iunitx==2) then
            write(*,"(a,f20.5,' eV')") " 8 Input full width at half maximum (FWHM), current:",FWHM(1)
        else if (iunitx==3) then
            write(*,"(a,f17.5,' 1000cm^-1')") " 8 Input full width at half maximum (FWHM), current:",FWHM(1)
        end if
    else
        write(*,*) "8 Set FWHM for all transitions, current: Inputted from file"
    end if
    if (ishowline==1) write(*,*) "9 Toggle showing discrete lines, current: ON"
    if (ishowline==0) write(*,*) "9 Toggle showing discrete lines, current: OFF"
    if (ispectrum==1) then
        if (iunitliney==1) write(*,*) "10 Switch the unit of infrared intensity, current: km/mol"
        if (iunitliney==2) write(*,*) "10 Switch the unit of infrared intensity, current: esu^2*cm^2"
    else if (ispectrum==3.or.ispectrum==4) then
        if (iunitx==1) write(*,*) "10 Set the unit of excitation energy, current: eV"
        if (iunitx==2) write(*,*) "10 Set the unit of excitation energy, current: nm"
        if (iunitx==3) write(*,*) "10 Set the unit of excitation energy, current: 1000cm^-1"
    end if
    if (ibroadfunc==3) write(*,"(a,f10.5)") " 11 Set Gaussian-weighting coefficient, current:",gauweigh
    write(*,"(a,f12.6)") " 12 Set shift value in X, current:",shiftx
    write(*,*) "13 Set colors of curve and discrete lines"
    if (ispectrum==1.or.ispectrum==2.or.ispectrum==5) write(*,*) "14 Multiply the vibrational frequencies by a factor"
    if (ispectrum==3.or.ispectrum==4) write(*,*) "14 Multiply the transition energies by a factor"
    write(*,*) "15 Output contribution of individual transition to the spectrum"
    write(*,*) "16 Find the positions of local minimum and maximum"
    if (ishowgrid==1) write(*,*) "17 Toggle showing dashed grid lines, current: ON"
    if (ishowgrid==0) write(*,*) "17 Toggle showing dashed grid lines, current: OFF"
    if (ispectrum==2.and.iramantype==1) write(*,*) "19 Convert Raman activities to intensities"
    if (ispectrum==2.and.iramantype==2) write(*,*) "19 Convert Raman intensities to activities"
    if (ispectrum==1.or.ispectrum==3) write(*,*) "20 Modify oscillator strengths"
    if (ispectrum==2.and.iramantype==1) write(*,*) "20 Modify Raman activities"
    if (ispectrum==2.and.iramantype==2) write(*,*) "20 Modify Raman intensities"
    if (ispectrum==4.or.ispectrum==5) write(*,*) "20 Modify rotatory strengths"
    read(*,*) isel
    
    
    if (isel==-3) then
        return
    else if (isel==-2) then
        open(10,file="transinfo.txt",status="replace")
        write(10,"(2i6)") numdata,2
        do i=1,numdata
            write(10,"(3f15.6)") datax(i),str(i),FWHM(i)
        end do
        close(10)
        write(*,"(a)") "The transition data have been saved to transinfo.txt in current directory, &
        this file can be directly used as input file of Multiwfn."
    else if (isel==-1.or.isel==20) then
        if (ispectrum==1) then
            write(*,*) " Index  Freq.(cm^-1)  Intens.( km/mol   esu^2*cm^2)"
            do i=1,numdata
                write(*,"(i6,1x,f12.5,7x,f12.5,f12.5)") i,datax(i),str(i),str(i)/2.5066D0
            end do
        else if (ispectrum==2) then
            if (iramantype==1) write(*,*) " Index  Freq.(cm^-1)      Activities(A^4/amu)"
            if (iramantype==2) write(*,*) " Index  Freq.(cm^-1)          Intensity"
            do i=1,numdata
                write(*,"(i6,3x,f12.5,7x,f12.5)") i,datax(i),str(i)
            end do
        else if (ispectrum==3.or.ispectrum==4) then
            if (ispectrum==3) write(*,*) " Index  Excit.energy(eV       nm         1000cm^-1)       Oscil.str."
            if (ispectrum==4) write(*,*) " Index  Excit.energy(eV       nm         1000cm^-1)       Rotat.str."
            if (iunitx==3) datax=datax/8.0655447D0
            do i=1,numdata
                write(*,"(i6,1x,4f15.5)") i,datax(i),1240.7011D0/datax(i),8.0655447D0*datax(i),str(i)
            end do
            if (iunitx==3) datax=datax*8.0655447D0
        else if (ispectrum==5) then
            write(*,*) " Index  Freq.(cm^-1)      Rotat.str."
            do i=1,numdata
                write(*,"(i6,3x,f12.5,7x,f12.5)") i,datax(i),str(i)
            end do
        end if
        if (isel==20) then
            write(*,*) "Input index range"
            write(*,*) "e.g. 1,3-6,22 means selecting mode 1,3,4,5,6,22"
            read(*,"(a)") c200tmp
            call str2arr(c200tmp,ntmparr)
            allocate(tmparr(ntmparr))
            call str2arr(c200tmp,ntmparr,tmparr)
            write(*,*) "Input value, e.g. 0.25"
            read(*,*) tmpval
            do i=1,ntmparr
                str(tmparr(i))=tmpval
            end do
            deallocate(tmparr)
        end if
    else if (isel==0) then
        idraw=1
        isavepic=0
    else if (isel==1) then
        idraw=1
        isavepic=1
    else if (isel==3) then
        write(*,*) "Input lower limit, upper limit and the step between labels e.g. 200,1700,150"
        read(*,*) xlow,xhigh,stepx
        iusersetX=1
    else if (isel==4) then
        write(*,*) "Input lower limit, upper limit and the step between labels e.g. 200,1700,150"
        read(*,*) orgy1,endy1,stepy1
        iusersetY1=1
    else if (isel==5) then
        write(*,*) "Input lower limit, upper limit and the step between labels e.g. 200,1700,150"
        read(*,*) orgy2,endy2,stepy2
        iusersetY2=1
    else if (isel==6) then
        write(*,*) "1 Lorentzian"
        write(*,*) "2 Gaussian"
        write(*,*) "3 Pseudo-Voigt"
        read(*,*) ibroadfunc
    else if (isel==7) then
        write(*,*) "Input a value  e.g. 0.5"
        read(*,*) scalecurve
    else if (isel==8) then
        if (ispectrum==1.or.ispectrum==2.or.ispectrum==5) then !Always use cm^-1
            write(*,*) "Input the FWHM in cm^-1, e.g. 4"
        else if (ispectrum==3.or.ispectrum==4) then
            if (iunitx==2) write(*,"(a,/)") " NOTE: nm is not a linear unit of energy, so in principle one cannot define FWHM in nm. Nevertheless, in Multiwfn, when nm &
            is chosen as the unit, the curve will be generated in eV as X-axis first, and then convert to nm. Since current unit is nm, now you have to define the FWHM in eV."
            if (iunitx==1.or.iunitx==2) write(*,*) "Input the FWHM in eV, e.g. 4"
            if (iunitx==3) write(*,*) "Input the FWHM in 1000cm^-1, e.g. 4"
        end if
        read(*,*) tmp
        FWHM=tmp !Note that FWHM is an array
    else if (isel==9) then
        if (ishowline==1) then
            ishowline=0
        else
            ishowline=1
        end if
    else if (isel==10) then
        if (ispectrum==1) then
            if (iunitliney==1) then
                iunitliney=2
            else
                iunitliney=1
            end if
        else if (ispectrum==3.or.ispectrum==4) then
            iusersetx=0 !Ensure auto X-axis is used, else will encounter problems 
            iold=iunitx
            write(*,*) "1: eV  2: nm  3: 1000cm^-1"
            read(*,*) iunitx
            if (iunitx==1.or.iunitx==2) then ! For eV and nm, datax are in eV. If select nm, we rescale linex and curvex to nm after broadening
                if (iold==3) then
                    datax=datax/8.0655447D0
                    FWHM=FWHM/8.0655447D0
                    scalecurve=scalecurve/8.0655447D0 !1 unit oscillator strength can be broadened to 28700 area(X:eV Y:L/mol/cm)
                end if
            else if (iunitx==3) then ! cm^-1
                if (iold==1.or.iold==2) then
                    datax=datax*8.0655447D0 !Convert data from eV to nm
                    FWHM=FWHM*8.0655447D0
                    scalecurve=scalecurve*8.0655447D0
                end if
            end if
        end if
    else if (isel==11) then
        write(*,*) "Input a value, e.g. 0.3"
        read(*,*) gauweigh
    else if (isel==12) then
        write(*,*) "Input a value, e.g. 4.5"
        read(*,*) shiftx
    else if (isel==13) then
        write(*,*) "Use which color for curve?"
        write(*,*) "1/2/3/4/5 = Red/Green/Blue/White/Black"
        write(*,*) "6/7/8/9/10 = Gray/Cyan/Yellow/Orange/Magenta"
        write(*,*) "11/12/13/14 = Crimson/Dark green/Purple/Brown"
        read(*,*) icurveclr
        write(*,*) "Use which color for discrete lines?"
        write(*,*) "1/2/3/4/5 = Red/Green/Blue/White/Black"
        write(*,*) "6/7/8/9/10 = Gray/Cyan/Yellow/Orange/Magenta"
        write(*,*) "11/12/13/14 = Crimson/Dark green/Purple/Brown"
        read(*,*) ilineclr
    else if (isel==14) then
        if (ispectrum==1.or.ispectrum==2.or.ispectrum==5) then
            write(*,*) "Input the index range of the transition modes you want to scaled"
            write(*,*) "e.g. 1,3-6,22 means selecting mode 1,3,4,5,6,22"
            write(*,"(a)") " Note: You can check the modes by option -1. Input ""all"" can select all modes. Input 0 can return."
            read(*,"(a)") c200tmp
            if (c200tmp(1:1)=='0'.or.c200tmp(1:2)=='-1') cycle !Avoid user input -1 here
            if (index(c200tmp,"all")==0) then
                call str2arr(c200tmp,nmode)
                allocate(tmparr(nmode))
                call str2arr(c200tmp,nmode,tmparr)
            else
                nmode=numdata
                allocate(tmparr(numdata))
                forall(itmp=1:numdata) tmparr(itmp)=itmp
            end if
            write(*,"(i6,' modes are selected')") nmode
            write(*,*) "Multiply the modes by which factor?  e.g. 0.9614 (suitable for B3LYP/6-31G*)"
            read(*,*) tmpval
            do idx=1,nmode
                datax(tmparr(idx))=datax(tmparr(idx))*tmpval
            end do
            deallocate(tmparr)
        else
            write(*,*) "Input the scale factor, e.g. 0.92"
            read(*,*) tmpval
            do i=1,numdata
                datax(i)=datax(i)*tmpval
            end do
        end if
        write(*,*) "Done! The excitation energies have been scaled"
    else if (isel==17) then
        if (ishowgrid==1) then
            ishowgrid=0
        else if (ishowgrid==0) then
            ishowgrid=1
        end if
    else if (isel==19) then
        write(*,*) "Input wavenumber of incident light (in cm^-1)"
        write(*,*) "e.g. 15797.79 (corresponding to 633nm He-Ne laser)"
        read(*,*) v0
        write(*,*) "Input temperature in K (Input 0 means ignoring the boltzmann term)"
        read(*,*) temper
        Cfac=1D-12
        do i=1,numdata
            vi=datax(i)
            if (temper==0) then
                Bfac=1D0
            else
!                 write(*,*) v0,lightc,planckc,boltzc,temper,-vi*100*lightc*planckc/(boltzc*temper)
                Bfac=1-exp( -vi*100*lightc*planckc/(boltzc*temper) )
            end if
            if (iramantype==1) str(i)= Cfac*(v0-vi)**4 /Bfac /vi * str(i) !Convert from activity to intensity
            if (iramantype==2) str(i)= str(i)*Bfac*vi /Cfac /(v0-vi)**4 !Convert from intensity to activity
        end do
        write(*,*) "Done!"
        if (iramantype==1) then
            iramantype=2
        else if (iramantype==2) then
            iramantype=1
        end if
    end if

    if (isel==0.or.isel==1.or.isel==2.or.isel==15.or.isel==16) then !These functions need calculation
        !======Construct correspondence array if output individual band
        if (isel==15) then
            write(*,"(a)") " Input criterion of strength, e.g. 0.2, the contribution of the bands having larger strength than this value will be outputted"
            read(*,*) critindband
            numindband=0
            do idata=1,numdata
                if (str(idata)>=critindband) numindband=numindband+1
            end do
            allocate(indband2idx(numindband),idx2indband(numdata),indcurve(num1Dpoints,numindband))
            indcurve=0D0
            itmp=0
            idx2indband=0
            do idata=1,numdata
                if (str(idata)<critindband) cycle
                itmp=itmp+1
                indband2idx(itmp)=idata !The ith band whose contribution will be outputted corresponds to the indband2idx(i) actual index
                idx2indband(idata)=itmp !If not =0, the contribution from the ith band will be stored to the idx2indband(i) slot of indcurve
            end do
        end if
        !======Determine upper and lower limit of X axis=======
        if (iusersetx==0) then !Automatical scale, if user hasn't manually set the range
            if (iunitx==0) then !cm^-1 for IR, Raman, VCD
                xlow=4000D0 !In common spectrum the energy is from high to low
                xhigh=0D0
            else if (iunitx==2) then !nm for UV-Vis, ECD, generate proper range in eV
                xhigh=1240.7011D0/ (1240.7011D0/minval(datax)+40) !Note that when nm is used, in common spectrum the energy is from high to low, so we invert xlow and xhigh
                xlow=1240.7011D0/ (1240.7011D0/maxval(datax)-40)
            else
                rangetmp=maxval(datax)-minval(datax)
                if (rangetmp==0) rangetmp=abs(datax(1))
                xlow=minval(datax)-0.3D0*rangetmp
                xhigh=maxval(datax)+0.3D0*rangetmp
            end if
        else if (iusersetx==1) then !The range was inputted by user
            !For nm, still using eV during broadening, so here we convert user inputted range from nm to eV, after broadening then rescale back
            if (iunitx==2) then
                xlow=1240.7011D0/xlow
                xhigh=1240.7011D0/xhigh
            end if
        end if
        !======Set x position of curves==========
        if (iunitx==2) then !For nm, which is not a linear energy unit, we generate points evenly distributed in X-axis with nm as unit
            xhighnm=1240.7011D0/xhigh
            xlownm=1240.7011D0/xlow
            xptstep=(xhighnm-xlownm)/(num1Dpoints-1) !Spacing in nm
            do ipoint=1,num1Dpoints
                curvex(ipoint)=1240.7011D0/(xlownm+(ipoint-1)*xptstep) !Now curvex is in eV
            end do
        else !Common case
            xptstep=(xhigh-xlow)/(num1Dpoints-1)
            do ipoint=1,num1Dpoints
                curvex(ipoint)=xlow+(ipoint-1)*xptstep
            end do
        end if
        !======Generate energy levels line=======
        liney=0D0
        do idata=1,numdata
            inow=3*(idata-1)
            linex(inow+1:inow+3)=datax(idata)
            liney(inow+2)=str(idata) !Unit of str array is km/mol
        end do
        !======Broadening=======
        !Under the current X and Y axes units, below code guarantee that the integral of the peak broadened by one unit of strength is 1
        curvey=0D0
        if (ibroadfunc==1.or.ibroadfunc==3) then !Lorentzian function, see http://mathworld.wolfram.com/LorentzianFunction.html
            do idata=1,numdata !Cycle each transition
                preterm=str(idata)*0.5D0/pi*FWHM(idata) !Integral of the peak equals to str(idata)
                do ipoint=1,num1Dpoints
                    curveytmp(ipoint)=preterm/( (curvex(ipoint)-datax(idata))**2+0.25D0*FWHM(idata)**2 )
                end do
                curvey=curvey+curveytmp
                if (isel==15) then
                    if (idx2indband(idata)/=0) indcurve(:,idx2indband(idata))=curveytmp !Individual contribution
                end if
            end do
        end if
        if (ibroadfunc==2.or.ibroadfunc==3) then !Gaussian function, see http://en.wikipedia.org/wiki/Gaussian_function
            if (ibroadfunc==3) then
                curvey=(1-gauweigh)*curvey !Scale Lorentzian function part of Pseudo-Voigt
                indcurve(:,:)=(1-gauweigh)*indcurve(:,:)
            end if
            do idata=1,numdata
                gauss_c=FWHM(idata)/2D0/sqrt(2*dlog(2D0))
                gauss_a=str(idata)/(gauss_c*sqrt(2D0*pi))
                do ipoint=1,num1Dpoints
                    curveytmp(ipoint)=gauss_a*dexp( -(curvex(ipoint)-datax(idata))**2/(2*gauss_c**2) )
                end do
                if (ibroadfunc==3) curveytmp=gauweigh*curveytmp !Scale Gaussian function part of Pseudo-Voigt
                curvey=curvey+curveytmp
                if (isel==15) then
                    if (idx2indband(idata)/=0) indcurve(:,idx2indband(idata))=indcurve(:,idx2indband(idata))+curveytmp !Individual contribution
                end if
            end do
        end if

        !Change units, scale curve, set axis
        if (ispectrum==1.and.iunitliney==2) liney=liney/2.5066D0 !!IR spectrum, convert strength unit from km/mol to cgs
        if (iunitx==2) then !Scaling plot in eV unit in X direction yields the plot in nm unit
            linex=1240.7011D0/linex
            curvex=1240.7011D0/curvex
            xlow=1240.7011D0/xlow
            xhigh=1240.7011D0/xhigh
        end if
        curvey=scalecurve*curvey
        if (isel==15) indcurve=scalecurve*indcurve
        curvex=curvex+shiftx
        if (iusersetx==0) stepx=(xhigh-xlow)/10
    end if

    !=======Output data to external text file======
    if (isel==2.or.isel==15) then
        !Output curve
        if (isel==2) then
            open(10,file="spectrum_curve.txt",status="replace")
            do ipt=1,num1Dpoints
                write(10,"(2f13.5)") curvex(ipt),curvey(ipt)
            end do
            close(10)
            write(*,*) "Curve data have been written to spectrum_curve.txt in current folder"
        else if (isel==15) then
            open(10,file="spectrum_curve.txt",status="replace")
            do ipt=1,num1Dpoints
                write(10,"(f13.5)",advance="no") curvex(ipt)
                write(10,"(f13.5)",advance="no") curvey(ipt)
                do iindband=1,numindband
                    write(10,"(f13.5)",advance="no") indcurve(ipt,iindband)
                end do
                write(10,*)
            end do
            close(10)
            write(*,"(a,i5,a)") " The total spectrum and the contributions from",numindband," transitions have been outputted to spectrum_curve.txt in current folder"
        end if
        !Show the units used in curve
        if (ispectrum==1) then
            write(*,*) "Column 1: Wavenumber (cm^-1)"
            write(*,*) "Column 2: Molar absorption coefficient (L/mol/cm)"
        else if (ispectrum==2) then
            write(*,*) "Column 1: Wavenumber (cm^-1)"
            write(*,*) "Column 2: Relative Raman intensity"
        else if (ispectrum==3.or.ispectrum==4) then
            if (iunitx==1) write(*,*) "Column 1: Excitation energy (eV)"
            if (iunitx==2) write(*,*) "Column 1: Wavelength (nm)"
            if (iunitx==3) write(*,*) "Column 1: Wavenumber (1000cm^-1)"
            if (ispectrum==3) write(*,*) "Column 2: Molar absorption coefficient (L/mol/cm)"
            if (ispectrum==4) write(*,*) "Column 2: Delta molar absorption coefficient (L/mol/cm)"
        else if (ispectrum==5) then
            write(*,*) "Column 1: Wavenumber (cm^-1)"
            write(*,*) "Column 2: Delta molar absorption coefficient (L/mol/cm)"
        end if
        if (isel==15) then
            write(*,*) "The correspondence between the columns and individual transitions in the file:"
            write(*,*) "     Column#   Transition#"
            do itmp=1,numindband
                write(*,"(2i12)") itmp+2,indband2idx(itmp)
            end do
            deallocate(indband2idx,idx2indband,indcurve)
        end if
        !Output discrete line
        open(10,file="spectrum_line.txt",status="replace")
        do i=1,3*numdata
            write(10,"(2f16.5)") linex(i)+shiftx,liney(i)
        end do
        close(10)
        write(*,*)
        write(*,*) "Discrete line data have been written to spectrum_line.txt in current folder"
        if (ispectrum==1) then
            write(*,*) "Column 1: Frequency (cm^-1)"
            if (iunitliney==1) write(*,*) "Column 2: IR intensities (km/mol)"
            if (iunitliney==2) write(*,*) "Column 2: IR intensities (esu^2*cm^2)"
        else if (ispectrum==2) then
            write(*,*) "Column 1: Frequency (cm^-1)"
            if (iramantype==1) write(*,*) "Column 2: Raman scattering activities (A^4/amu)"
            if (iramantype==2) write(*,*) "Column 2: Raman scattering intensities"
        else if (ispectrum==3.or.ispectrum==4) then
            if (iunitx==1) write(*,*) "Column 1: Excitation energy (eV)"
            if (iunitx==2) write(*,*) "Column 1: Wavelength (nm)"
            if (iunitx==3) write(*,*) "Column 1: Wavenumber (1000cm^-1)"
            if (ispectrum==3) write(*,*) "Column 2: Oscillator strength"
            if (ispectrum==4) write(*,*) "Column 2: Rotatory strength in cgs (10^-40 erg-esu-cm/Gauss)"
        else if (ispectrum==5) then
            write(*,*) "Column 1: Frequency (cm^-1)"
            write(*,*) "Column 2: Rotatory strength (10^-44 esu^2 cm^2)"
        end if
    else if (isel==16) then
        numlocmax=0
        do ipoint=2,num1Dpoints-1
            gradold=curvey(ipoint)-curvey(ipoint-1)
            gradnew=curvey(ipoint+1)-curvey(ipoint)
            if (gradold*gradnew<0D0.and.gradold>gradnew) then
                numlocmax=numlocmax+1
                write(*,"(' Local maximum X:',f15.4,'      Value:',f15.4)") curvex(ipoint),curvey(ipoint)
            end if
        end do
        write(*,*)
        numlocmin=0
        do ipoint=2,num1Dpoints-1
            gradold=curvey(ipoint)-curvey(ipoint-1)
            gradnew=curvey(ipoint+1)-curvey(ipoint)
            if (gradold*gradnew<0D0.and.gradold<gradnew) then
                numlocmin=numlocmin+1
                write(*,"(' Local minimum X:',f15.4,'      Value:',f15.4)") curvex(ipoint),curvey(ipoint)
            end if
        end do
        write(*,"(/,' Totally found',i5,' local minimum,',i5,' local maximum')") numlocmin,numlocmax
    end if

    !======Draw spectrum now=======
    if (idraw==1) then
        if (isavepic==0) then
        else if (isavepic==1) then
        end if
        ! Name of X-axis
        ! Name of Y-axis
        if (ispectrum==1.or.ispectrum==3) then
        else if (ispectrum==2) then
        else if (ispectrum==4.or.ispectrum==5) then
        end if

        if (iusersetY1==0) then
            endy1=1.1D0*maxval(abs(curvey)) !Default lower and upper limit of Y axis
            orgy1=-endy1/30D0
            if (ispectrum==4.or.ispectrum==5) orgy1=-endy1 !Positive and negative of ECD curve may have similar height
            stepy1=(endy1-orgy1)/10
        end if
        if (iusersetY2==0) then
            endy2=1.1D0*maxval(abs(liney)) !Default lower and upper limit of Y axis
            orgy2=-endy2/30D0
            if (ispectrum==4.or.ispectrum==5) orgy2=-endy2
            stepy2=(endy2-orgy2)/10
        end if
        
        if (ishowline==1) then
            if (ispectrum==1) then
            else if (ispectrum==2) then
            else if (ispectrum==3) then
            else if (ispectrum==4) then
            else if (ispectrum==5) then
            end if
        end if
        if (isavepic==1) write(*,*) "Graph file has been saved to current folder with ""DISLIN"" prefix"
    end if
    idraw=0

end do
end subroutine
