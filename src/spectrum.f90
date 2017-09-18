!!!--------- Plot various kinds of spectra
!multiple.txt can records input file of multiple systems, but the system with the maximum number of transitions must be presented as the first entry
!If all(weight=1), then it is assumed that multiple.txt is not oriented for plotting weighted spectrum along with individual spectrum for each system, &
!but only for the latter, and in this case custom legend is allowed, which can be written in multiple.txt as the second column
subroutine plotspectrum
use defvar
use util
implicit real*8 (a-h,o-z)
real*8,allocatable :: weight(:) !Weight of various system for plotting mixed spectrum
real*8,allocatable :: dataxall(:,:),strall(:,:),FWHMall(:,:) !Transition data loaded from multiple files. The first index corresponds to system index
integer,allocatable :: numdataall(:)
character*80,allocatable :: mollegend(:) 
real*8,allocatable :: linexall(:,:),lineyall(:,:) !Array used to draw discrete lines for all systems. The first index corresponds to system index
real*8,allocatable :: curveyall(:,:) !The first index corresponds to system index. curvey in global array is used to record weighted curve
integer,allocatable :: tmparr(:)
real*8,allocatable :: indcurve(:,:) !Y value of curve of each individual band
integer,allocatable :: indband2idx(:),idx2indband(:) !Used to map individual band index
character c200tmp*200,c200tmp2*200,selectyn
character clegend*2000 !Buffer for showing legends
integer :: icurveclr=1,ilineclr=5 !Default: Red for curve, black for discrete lines

!Spectrum types: 1=IR  2=Raman (or pre-resonance Raman)  3=UV-Vis  4=ECD  5=VCD
!
!Definition of units:  iunitx =0 cm^-1, =1 eV, =2 nm, =3 1000cm^-1
!For vibrational spectrum, cm^-1 is always used; For electronic spectrum, eV, nm, 1000cm^-1 can be used 
!Note: For nm unit, we still store all data and FWHM in eV, and generate curve as usual in eV. Only at final stage, we scale the curve to get the one in nm
!If we choose 1000cm^-1, we immediately convert all data and FWHM into 1000cm^-1 before generating curve.
!When unit is changed, we reset lower and upper limit to auto rather than convert them to current unit to avoid problems.
!
!IR may use esu^2*cm^2 or km/mol as Y-axis unit, the data is always recorded in km/mol
!Unit conversion: 1eV=8.0655447*1000cm^-1    (1240.7011/nm)eV    (1240.7011/eV)nm     1esu^2*cm^2=2.5066km/mol

!! Initialize variables
gauweigh=0.5D0 !Gaussian weight used in Pseudo-Voigt broadening
iusersetY1=0 !User has not set the axes definition by himself
iusersetY2=0
iusersetX=0
idraw=0
isavepic=0
ishowline=1
ishowgrid=1
ishowtotal=1 !If showing weighted spectrum
iweisyscurve=0
iunitliney=1 !Only for IR
shiftx=0D0 !Shift value in X direction
iramantype=1 !=1 Raman activities   =2 Raman intensities
write(*,*) "Select type of the spectrum"
write(*,*) "1:IR  2:Raman (or pre-resonance Raman)  3:UV-Vis  4:ECD  5:VCD"
read(*,*) ispectrum
if (ispectrum==1.or.ispectrum==2.or.ispectrum==5) then !IR, Raman and VCD
    ibroadfunc=1 !Use Lorentzian broadening
    iunitx=0 !cm^-1
else if (ispectrum==3.or.ispectrum==4) then !UV-Vis, ECD
    ibroadfunc=2 !Use Gaussian broadening
    iunitx=2 !nm is default unit. But transition energies are loaded as eV
end if

if (ispectrum==2.or.ispectrum==4.or.ispectrum==5) then
!For ECD when eV is used, integrating the peak of a unit strength is assumed to be 1
!For Raman and VCD (cm^-1 is used), integrating the peak of a unit strength is assumed to be 1
    scalecurve=1D0
else if (ispectrum==1) then
!For IR when km/L is used, integrating the peak of a unit strength is 100
    scalecurve=100D0
else if (ispectrum==3) then
!1 unit oscillator strength can be broadened to 28700 area (X:eV Y:L/mol/cm)
!1 unit oscillator strength can be broadened to 1D0/4.32D-6 area (X:1000cm^-1 Y:L/mol/cm)
!Their relationship: 1/(4.32D-9)/8065.5447=28700  1eV=8.0655447*1000*cm^-1 
!4.32D-9 can be found from Swizard manual (see also Review in C.C. vol.20 p168)
!The result is consistent with Gaussview
    if (iunitx==1.or.iunitx==2) scalecurve=28700 !nm,eV, both are recorded as eV internally
    if (iunitx==3) scalecurve=1D0/4.32D-6 !1000cm^-1
end if


!! Load transition data from external files
nsystem=0
if (index(filename,"multiple.txt")/=0) then !Multiple file list with weights is recorded in multiple.txt
    open(11,file=filename,status="old") !Note that 10 is used by loadtransdata
    !Count total number of entries
    do while(.true.)
        read(11,*,iostat=ierror) c200tmp
        if (ierror/=0.or.c200tmp==" ") exit
        nsystem=nsystem+1
    end do
    write(*,"(' There are',i4,' systems')") nsystem
    allocate(weight(nsystem),mollegend(nsystem))
    mollegend=" "
    rewind(11)
    do i=1,nsystem
        read(11,"(a)") c200tmp
        c200tmp2=c200tmp
        read(c200tmp,*,iostat=ierror) c200tmp,weight(i)
        if (ierror/=0) then
            ispc=index(c200tmp," ")
            read(c200tmp2(:ispc-1),*) c200tmp
            read(c200tmp2(ispc+1:),"(a)") mollegend(i)
            weight(i)=1
        end if
        inquire(file=c200tmp,exist=alive)
        if (alive.eqv. .false.) then
            write(*,"(' Error: Cannot find ',a)") trim(c200tmp)
            write(*,*) "Press ENTER to exit program"
            read(*,*)
            stop
        end if
        if (weight(i)==1) then
            write(*,"(' Loading ',a,'    Legend: ',a)") trim(c200tmp),trim(mollegend(i))
        else
            write(*,"(' Loading ',a,'    Weight:',f7.4)") trim(c200tmp),weight(i)
        end if
        call loadtransdata(ispectrum,trim(c200tmp),numdata) !Data are loaded into datax,str,FWHM in global memory
        if (i==1) allocate(dataxall(nsystem,numdata),strall(nsystem,numdata),FWHMall(nsystem,numdata),numdataall(nsystem))
        if (numdata>size(dataxall,2)) then !Because we need to allocate enough slots, while the number of slots is allocated when loading the first file
            write(*,*)
            write(*,"(' The number of transitions in this file:',i6)") numdata
            write(*,"(' The number of transitions in the first file:',i6)") size(dataxall,2)
            write(*,"(a)") " Error: You should put the system with maximum number of transitions to the first entry of multiple.txt"
            write(*,*) "Press ENTER to exit program"
            read(*,*)
            stop
        end if
        dataxall(i,1:numdata)=datax
        strall(i,1:numdata)=str
        FWHMall(i,1:numdata)=FWHM
        numdataall(i)=numdata
    end do
    close(11)
    if (all(weight==1)) then !When all weights are unity, then no weighted spectrum will be plotted, but simply plotting all systems together
        ishowtotal=0
        iweisyscurve=1
    end if
    !Below, numdata indicates maximum number of numdataall array
    numdata=maxval(numdataall)
else !Only one system
    nsystem=1
    allocate(weight(1))
    weight=1D0
    call loadtransdata(ispectrum,filename,numdata)
    allocate(dataxall(nsystem,numdata),strall(nsystem,numdata),FWHMall(nsystem,numdata),numdataall(nsystem))
    dataxall(1,:)=datax
    strall(1,:)=str
    FWHMall(1,:)=FWHM
    numdataall(1)=numdata
end if


!! Allocate arrays properly
! curvey in global array is used to record weighted curve
! curveytmp is a temporary array used to record curve from each transition during generating curve of each system
if (allocated(curvex)) deallocate(curvex) !Global array
if (allocated(curvey)) deallocate(curvey) !Global array
if (allocated(curveytmp)) deallocate(curveytmp) !Global array
allocate(curvex(num1Dpoints),curvey(num1Dpoints),curveytmp(num1Dpoints),curveyall(nsystem,num1Dpoints))
allocate(linexall(nsystem,3*numdata),lineyall(nsystem,3*numdata))
 

!! Main interface
do while(.true.)
    write(*,*)
    write(*,*) "-3 Return to main menu"
    write(*,*) "-2 Export transition data to plain text file"
    write(*,*) "-1 Show transition data"
    write(*,*) "0 Plot spectrum"
    write(*,*) "1 Save graphical file of the spectrum in current folder"
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
    do imol=1,nsystem
        if (any(FWHMall(imol,1:numdataall(imol))/=FWHMall(1,1))) exit
    end do
    if (imol==nsystem+1) then
        if (iunitx==0) then
            write(*,"(a,f20.5,' cm^-1')") " 8 Input full width at half maximum (FWHM), current:",FWHMall(1,1)
        else if (iunitx==1.or.iunitx==2) then
            !FWHM cannot be defined for nm, since it is not a linear unit, so what inputted is eV
            write(*,"(a,f20.5,' eV')") " 8 Input full width at half maximum (FWHM), current:",FWHMall(1,1)
        else if (iunitx==3) then
            write(*,"(a,f17.5,' 1000cm^-1')") " 8 Input full width at half maximum (FWHM), current:",FWHMall(1,1)
        end if
    else
        write(*,*) "8 Set FWHM for all transitions, current: Loaded from input file"
    end if
    if (ishowline==1) write(*,*) "9 Toggle showing discrete lines, current: ON"
    if (ishowline==0) write(*,*) "9 Toggle showing discrete lines, current: OFF"
    if (ispectrum==1) then !IR allows using different strength units
        if (iunitliney==1) write(*,*) "10 Switch the unit of infrared intensity, current: km/mol"
        if (iunitliney==2) write(*,*) "10 Switch the unit of infrared intensity, current: esu^2*cm^2"
    else if (ispectrum==3.or.ispectrum==4) then !UV-Vis or ECD allows using different energy units
        if (iunitx==1) write(*,*) "10 Set the unit of excitation energy, current: eV"
        if (iunitx==2) write(*,*) "10 Set the unit of excitation energy, current: nm"
        if (iunitx==3) write(*,*) "10 Set the unit of excitation energy, current: 1000cm^-1"
    end if
    if (ibroadfunc==3) write(*,"(a,f10.5)") " 11 Set Gaussian-weighting coefficient, current:",gauweigh
    write(*,"(a,f12.6)") " 12 Set shift value in X, current:",shiftx
    if (nsystem==1.or.any(weight/=1)) then
        write(*,*) "13 Set colors of curve and discrete lines"
    end if
    if (ispectrum==1.or.ispectrum==2.or.ispectrum==5) write(*,*) "14 Multiply the vibrational frequencies by a factor"
    if (ispectrum==3.or.ispectrum==4) write(*,*) "14 Multiply the transition energies by a factor"
    if (nsystem==1) write(*,*) "15 Output contribution of individual transition to the spectrum"
    if (.not.(nsystem>1.and.all(weight==1))) write(*,*) "16 Find the positions of local minimum and maximum"
    if (ishowgrid==1) write(*,*) "17 Toggle showing dashed grid lines, current: ON"
    if (ishowgrid==0) write(*,*) "17 Toggle showing dashed grid lines, current: OFF"
    if (nsystem>1.and.any(weight/=1)) then
        if (iweisyscurve==1) write(*,*) "18 Toggle weighting spectrum of each system, current: ON"
        if (iweisyscurve==0) write(*,*) "18 Toggle weighting spectrum of each system, current: OFF"
    end if
    if (ispectrum==2.and.iramantype==1) write(*,*) "19 Convert Raman activities to intensities"
    if (ispectrum==2.and.iramantype==2) write(*,*) "19 Convert Raman intensities to activities"
    if (ispectrum==1.or.ispectrum==3) write(*,*) "20 Modify oscillator strengths"
    if (ispectrum==2.and.iramantype==1) write(*,*) "20 Modify Raman activities"
    if (ispectrum==2.and.iramantype==2) write(*,*) "20 Modify Raman intensities"
    if (ispectrum==4.or.ispectrum==5) write(*,*) "20 Modify rotatory strengths"
    read(*,*) isel
    
    if (isel==-3) then
        return
    else if (isel==-2) then !Export transition data
        if (nsystem==1) then
            open(10,file="transinfo.txt",status="replace")
            write(10,"(2i6)") numdata,2
            do i=1,numdata
                write(10,"(3f15.6)") dataxall(1,i),strall(1,i),FWHMall(1,i)
            end do
            close(10)
            write(*,"(a)") " The transition data have been exported to transinfo.txt in current directory, &
            this file can be directly used as input file of Multiwfn."
        else
            do imol=1,nsystem
                write(c200tmp,"(a,i3.3,a)") "transinfo",imol,".txt"
                open(10,file=c200tmp,status="replace")
                write(10,"(2i6)") numdataall(imol),2
                do i=1,numdataall(imol)
                    write(10,"(3f15.6)") dataxall(imol,i),strall(imol,i),FWHMall(imol,i)
                end do
                close(10)
            end do
            write(*,"(a)") " The transition data have been exported to .txt with ""transinfo"" as prefix in current directory, &
            these files can be directly used as input file of Multiwfn."
        end if
    else if (isel==-1.or.isel==20) then !Show transition data and modify strengths
        do imol=1,nsystem
            if (nsystem>1) write(*,"(/,' Transition data of system',i5)") imol
            if (ispectrum==1) then !IR
                write(*,*) " Index  Freq.(cm^-1)  Intens.( km/mol   esu^2*cm^2)"
                do i=1,numdataall(imol)
                    write(*,"(i6,1x,f12.5,7x,f12.5,f12.5)") i,dataxall(imol,i),strall(imol,i),strall(imol,i)/2.5066D0
                end do
            else if (ispectrum==2) then !Raman
                if (iramantype==1) write(*,*) " Index  Freq.(cm^-1)      Activities(A^4/amu)"
                if (iramantype==2) write(*,*) " Index  Freq.(cm^-1)          Intensity"
                do i=1,numdataall(imol)
                    write(*,"(i6,3x,f12.5,7x,f12.5)") i,dataxall(imol,i),strall(imol,i)
                end do
            else if (ispectrum==3.or.ispectrum==4) then !UV-Vis, ECD
                if (ispectrum==3) write(*,*) " Index  Excit.energy(eV       nm         1000cm^-1)       Oscil.str."
                if (ispectrum==4) write(*,*) " Index  Excit.energy(eV       nm         1000cm^-1)       Rotat.str."
                if (iunitx==3) dataxall=dataxall/8.0655447D0 !If unit is in 1000cm^-1, temporarily convert to eV
                do i=1,numdataall(imol)
                    write(*,"(i6,1x,4f15.5)") i,dataxall(imol,i),1240.7011D0/dataxall(imol,i),8.0655447D0*dataxall(imol,i),strall(imol,i)
                end do
                if (iunitx==3) dataxall=dataxall*8.0655447D0 !Convert back from eV to 1000cm^-1
            else if (ispectrum==5) then !VCD
                write(*,*) " Index  Freq.(cm^-1)      Rotat.str."
                do i=1,numdataall(imol)
                    write(*,"(i6,3x,f12.5,7x,f12.5)") i,dataxall(imol,i),strall(imol,i)
                end do
            end if
        end do
        if (isel==20) then !Modify strengths
            imol=1
            if (nsystem>1) then
                write(*,*)
                write(*,*) "Input index of system, e.g. 2"
                read(*,*) imol
            end if
            write(*,*) "Input index range of transitions"
            write(*,*) "e.g. 1,3-6,22 means selecting mode 1,3,4,5,6,22"
            read(*,"(a)") c200tmp
            call str2arr(c200tmp,ntmparr)
            allocate(tmparr(ntmparr))
            call str2arr(c200tmp,ntmparr,tmparr)
            write(*,*) "Input value, e.g. 0.25"
            read(*,*) tmpval
            do i=1,ntmparr
                strall(imol,tmparr(i))=tmpval
            end do
            deallocate(tmparr)
        end if
    else if (isel==0) then !Draw curve
        idraw=1
        isavepic=0
    else if (isel==1) then !Save curve picture
        idraw=1
        isavepic=1
    else if (isel==3) then !Change X axis
        write(*,*) "Input lower limit, upper limit and step between ticks e.g. 200,1700,150"
        read(*,*) xlow,xhigh,stepx
        iusersetX=1 !User has modified it
    else if (isel==4) then !Change left Y axis
        orgy1old=orgy1
        endy1old=endy1
        stepy1old=stepy1
        write(*,*) "Input lower limit, upper limit and step between ticks e.g. 0,17000,2000"
        read(*,*) orgy1,endy1,stepy1
        iusersetY1=1
        write(*,"(a)") " Do you want to let program properly scale right Y axis so that its zero position exactly corresponds to left Y-axis? (y/n)"
        read(*,*) selectyn
        if (selectyn=='y'.or.selectyn=='Y') then
            iusersetY2=1
            ratiotmp=(endy1old-orgy1old)/(endy2-orgy2)
            endy2=endy1/ratiotmp
            orgy2=orgy1/ratiotmp
            stepy2=stepy1/ratiotmp
        end if
    else if (isel==5) then !Change right Y axis
        orgy2old=orgy2
        endy2old=endy2
        stepy2old=stepy2
        write(*,*) "Input lower limit, upper limit and step between ticks e.g. -10,40,5"
        read(*,*) orgy2,endy2,stepy2
        iusersetY2=1
        write(*,"(a)") " Do you want to let program properly scale left Y axis so that its zero position exactly corresponds to right Y-axis? (y/n)"
        read(*,*) selectyn
        if (selectyn=='y'.or.selectyn=='Y') then
            iusersetY1=1
            ratiotmp=(endy1-orgy1)/(endy2old-orgy2old)
            endy1=endy2*ratiotmp
            orgy1=orgy2*ratiotmp
            stepy1=stepy2*ratiotmp
        end if
    else if (isel==6) then !Set broadening function
        write(*,*) "1 Lorentzian"
        write(*,*) "2 Gaussian"
        write(*,*) "3 Pseudo-Voigt"
        read(*,*) ibroadfunc
    else if (isel==7) then !Scale factor for curve
        write(*,*) "Input the scale factor, e.g. 0.8"
        read(*,*) scalecurve
    else if (isel==8) then !Set FWHM
        if (ispectrum==1.or.ispectrum==2.or.ispectrum==5) then !Always use cm^-1
            write(*,*) "Input the FWHM in cm^-1, e.g. 4"
        else if (ispectrum==3.or.ispectrum==4) then
            if (iunitx==2) write(*,"(a,/)") " NOTE: nm is not a linear unit of energy, so in principle one cannot define FWHM in nm. Nevertheless, in Multiwfn, when nm &
            is chosen as the unit, the curve will be generated in eV as X-axis first, and then convert to nm. Since current unit is nm, now you have to define the FWHM in eV."
            if (iunitx==1.or.iunitx==2) write(*,*) "Input the FWHM in eV, e.g. 0.5"
            if (iunitx==3) write(*,*) "Input the FWHM in 1000cm^-1, e.g. 4"
        end if
        read(*,*) tmp
        FWHMall=tmp
    else if (isel==9) then !If show discrete lines
        if (ishowline==1) then
            ishowline=0
        else
            ishowline=1
        end if
    else if (isel==10) then !Change unit of X or Y axis. Not applied to VCD/ECD
        if (ispectrum==1) then !IR
            if (iunitliney==1) then
                iunitliney=2
            else
                iunitliney=1
            end if
        else if (ispectrum==3.or.ispectrum==4) then !UV-Vis, ECD
            iusersetx=0 !Ensure auto X-axis is used, else will encounter problems 
            iold=iunitx
            write(*,*) "1: eV  2: nm  3: 1000cm^-1"
            read(*,*) iunitx
            if (iunitx==1.or.iunitx==2) then
                if (iold==3) then !Convert data from 1000cm^-1 to eV (For both eV and nm, transition energies are stored in eV)
                    dataxall=dataxall/8.0655447D0
                    FWHMall=FWHMall/8.0655447D0
                    scalecurve=scalecurve/8.0655447D0 !1 unit oscillator strength can be broadened to 28700 area (X:eV Y:L/mol/cm)
                end if
            else if (iunitx==3) then
                if (iold==1.or.iold==2) then !Convert data from eV to 1000cm^-1
                    dataxall=dataxall*8.0655447D0
                    FWHMall=FWHMall*8.0655447D0
                    scalecurve=scalecurve*8.0655447D0
                end if
            end if
        end if
    else if (isel==11) then !Weight of Gaussian function
        write(*,*) "Input a value, e.g. 0.3"
        read(*,*) gauweigh
    else if (isel==12) then !Shift value in X
        write(*,*) "Input a value, e.g. 4.5"
        read(*,*) shiftx
    else if (isel==13) then !Set line and curve colors
        if (nsystem==1) write(*,*) "Use which color for curve?"
        if (nsystem>1) write(*,*) "Use which color for weighted curve?"
        write(*,*) "1/2/3/4/5 = Red/Green/Blue/White/Black"
        write(*,*) "6/7/8/9/10 = Gray/Cyan/Yellow/Orange/Magenta"
        write(*,*) "11/12/13/14 = Crimson/Dark green/Purple/Brown"
        read(*,*) icurveclr
        write(*,*) "Use which color for discrete lines?"
        write(*,*) "1/2/3/4/5 = Red/Green/Blue/White/Black"
        write(*,*) "6/7/8/9/10 = Gray/Cyan/Yellow/Orange/Magenta"
        write(*,*) "11/12/13/14 = Crimson/Dark green/Purple/Brown"
        read(*,*) ilineclr
    else if (isel==14) then !Set scale factor for transition energies or frequencies
        if (nsystem>1) then
            write(*,*) "Note: This operation will be applied to all systems loaded"
            write(*,*)
        end if
        if (ispectrum==1.or.ispectrum==2.or.ispectrum==5) then !Vibrational spectra, mode selection is viable
            write(*,*) "Input the index range of the transitions you want to scaled"
            write(*,*) "e.g. 1,3-6,22 means selecting transitions 1,3,4,5,6,22"
            write(*,"(a)") " Note: You can check all transitions by option -1. Input ""all"" can select all modes. Input 0 can return."
            read(*,"(a)") c200tmp
            if (c200tmp(1:1)=='0'.or.c200tmp(1:2)=='-1') cycle !Avoid foolish user inputting -1 here
            if (index(c200tmp,"all")==0) then
                call str2arr(c200tmp,nmode)
                allocate(tmparr(nmode))
                call str2arr(c200tmp,nmode,tmparr)
            else !Selected all modes
                nmode=numdata
                allocate(tmparr(numdata))
                forall(itmp=1:numdata) tmparr(itmp)=itmp
            end if
            write(*,"(i6,' frequencies are selected')") nmode
            write(*,"(a)") " Multiplying the frequencies by which factor?  e.g. 0.9614 (suitable for B3LYP/6-31G*)"
            read(*,*) tmpval
            do idx=1,nmode
                dataxall(:,tmparr(idx))=dataxall(:,tmparr(idx))*tmpval
            end do
            deallocate(tmparr)
        else !Electronic spectra, use universal scaling
            write(*,*) "Input the scale factor, e.g. 0.92"
            read(*,*) tmpval
            dataxall=dataxall*tmpval
        end if
        write(*,*) "Done! Transition energies have been scaled"
    else if (isel==17) then !If showing grids on the plot
        if (ishowgrid==1) then
            ishowgrid=0
        else if (ishowgrid==0) then
            ishowgrid=1
        end if
    else if (isel==18) then !If weighting curve of each system
        if (iweisyscurve==1) then
            iweisyscurve=0
        else if (iweisyscurve==0) then
            iweisyscurve=1
        end if
    else if (isel==19) then !Convert between Raman activity and Raman intensity
        if (nsystem>1) then
            write(*,*) "Note: This operation will be applied to all systems loaded"
            write(*,*)
        end if
        write(*,*) "Input wavenumber of incident light (in cm^-1)"
        write(*,*) "e.g. 15797.79 (corresponding to 633nm He-Ne laser)"
        read(*,*) v0
        write(*,*) "Input temperature in K (Input 0 means ignoring the Boltzmann term)"
        read(*,*) temper
        Cfac=1D-12
        write(*,"(' Note: C factor of',1PD16.8,' is used')") Cfac
        do imol=1,nsystem
            do i=1,numdataall(imol)
                vi=dataxall(imol,i)
                if (temper==0) then
                    Bfac=1D0
                else
    !                 write(*,*) v0,lightc,planckc,boltzc,temper,-vi*100*lightc*planckc/(boltzc*temper)
                    Bfac=1-exp( -vi*100*lightc*planckc/(boltzc*temper) )
                end if
                if (iramantype==1) strall(imol,i)= Cfac*(v0-vi)**4 /Bfac /vi * strall(imol,i) !Convert from activity to intensity
                if (iramantype==2) strall(imol,i)= strall(imol,i)*Bfac*vi /Cfac /(v0-vi)**4 !Convert from intensity to activity
            end do
        end do
        write(*,*) "Done!"
        if (iramantype==1) then
            iramantype=2
        else if (iramantype==2) then
            iramantype=1
        end if
    else if (isel==21) then !If weighting curve of each system
        if (ishowtotal==1) then
            ishowtotal=0
        else if (ishowtotal==0) then
            ishowtotal=1
        end if
    end if
    
    if (isel==15.and.nsystem>1) then !Showing individual transition contribution is not possible when multiple files are involved
        write(*,*) "Error: This function is not available when multiple files are involved!"
        write(*,*) "Press ENTER to continue"
        read(*,*)
        cycle
    end if
    

    !!=======================================================================!!
    !!=======================================================================!!
    !!=============== Below functions need calculation of curves ============!!
    !!=======================================================================!!
    !!=======================================================================!!
    if (isel==0.or.isel==1.or.isel==2.or.isel==15.or.isel==16) then
        !====== Construct correspondence array if outputting individual bands. Only available when one file is loaded
        !This function is not available when multiple systems are considered
        if (isel==15) then
            write(*,"(a)") " Input criterion of strength, e.g. 0.2, the contribution of the transitions whose absolute value of strength larger than this value will be outputted"
            read(*,*) critindband
            numindband=0 !The number of individual bands satisfying criterion
            do idata=1,numdata
                if (abs(strall(1,idata))>=critindband) numindband=numindband+1
            end do
            allocate(indband2idx(numindband),idx2indband(numdata),indcurve(num1Dpoints,numindband))
            indcurve=0D0
            itmp=0
            idx2indband=0
            do idata=1,numdata
                if (abs(strall(1,idata))<critindband) cycle
                itmp=itmp+1
                indband2idx(itmp)=idata !Map index of outputted bands (indband) to actual index of all transitions (idx)
                idx2indband(idata)=itmp !If not =0, the contribution from i band will be stored to idx2indband(i) slot of indcurve array
            end do
        end if
        !====== Determine upper and lower limit of X axis =======
        if (iusersetx==0) then !Automatical scale, if user has not manually set the range
            tmpmin=minval(dataxall(1,1:numdataall(1)))
            tmphigh=maxval(dataxall(1,1:numdataall(1)))
            if (nsystem>1) then !Find upper and lower values for all systems
                do imol=2,nsystem
                    tmpa=minval(dataxall(imol,1:numdataall(imol)))
                    tmpb=maxval(dataxall(imol,1:numdataall(imol)))
                    if (tmpa<tmpmin) tmpmin=tmpa
                    if (tmpb>tmpmax) tmpmax=tmpb
                end do
            end if
            if (iunitx==0) then !cm^-1 for IR, Raman, VCD
                xlow=4000D0 !In common spectrum the energy is from high to low
                xhigh=0D0
            else if (iunitx==2) then !nm for UV-Vis, ECD, generate proper range in eV
                xhigh=1240.7011D0/ (1240.7011D0/tmpmin+40) !Note that when nm is used, in common spectrum the energy is from high to low, so we invert xlow and xhigh
                xlow=1240.7011D0/ (1240.7011D0/tmphigh-40)
            else
                rangetmp=tmphigh-tmpmin
                if (rangetmp==0) rangetmp=abs(dataxall(1,1)) !Only one data
                xlow=tmpmin-0.3D0*rangetmp
                xhigh=tmphigh+0.3D0*rangetmp
            end if
        else if (iusersetx==1) then !The range was defined by user
            !nm is selected unit, however we still using eV during broadening, so we convert user inputted range from nm to eV, after broadening then convert back
            if (iunitx==2) then
                xlow=1240.7011D0/xlow
                xhigh=1240.7011D0/xhigh
            end if
        end if
        !====== Set x position of curves ==========
        if (iunitx==2) then !For nm, which is not a linear energy unit, we generate points evenly distributed in X-axis with nm as unit
            xhighnm=1240.7011D0/xhigh !nm->eV
            xlownm=1240.7011D0/xlow
            xptstep=(xhighnm-xlownm)/(num1Dpoints-1) !Get proper spacing in nm
            do ipoint=1,num1Dpoints
                curvex(ipoint)=1240.7011D0/(xlownm+(ipoint-1)*xptstep) !Now curvex is recorded in eV
            end do
        else !Common case, linear energy unit is used
            xptstep=(xhigh-xlow)/(num1Dpoints-1)
            do ipoint=1,num1Dpoints
                curvex(ipoint)=xlow+(ipoint-1)*xptstep
            end do
        end if
        !====== Generate energy levels line =======
        lineyall=0D0
        linexall=1D0 !To garantee that linexall will not be zero, otherwise may crash when converting unit
        do imol=1,nsystem
            do idata=1,numdataall(imol)
                inow=3*(idata-1)
                linexall(imol,inow+1:inow+3)=dataxall(imol,idata)
                lineyall(imol,inow+2)=weight(imol)*strall(imol,idata) !Line height is weighted! Otherwise plotting them is meaningless
            end do
        end do
        !===============================================
        !============== Broaden to curve ===============
        !===============================================
        !Under current X and Y axes units, below code guarantees that the integral of the peak broadened by one unit of strength is 1
        curveyall=0D0
        if (ibroadfunc==1.or.ibroadfunc==3) then !Lorentzian function, see http://mathworld.wolfram.com/LorentzianFunction.html
            do imol=1,nsystem
                do idata=1,numdataall(imol) !Cycle each transition
                    preterm=strall(imol,idata)*0.5D0/pi*FWHMall(imol,idata) !Integral of the peak equals to str(idata)
                    do ipoint=1,num1Dpoints
                        curveytmp(ipoint)=preterm/( (curvex(ipoint)-dataxall(imol,idata))**2+0.25D0*FWHMall(imol,idata)**2 )
                    end do
                    curveyall(imol,:)=curveyall(imol,:)+curveytmp
                    if (isel==15) then !Individual contribution
                        if (idx2indband(idata)/=0) indcurve(:,idx2indband(idata))=curveytmp
                    end if
                end do
            end do
        end if
        if (ibroadfunc==2.or.ibroadfunc==3) then !Gaussian or Pseudo-Voigt function, see http://en.wikipedia.org/wiki/Gaussian_function
            if (ibroadfunc==3) then
                curveyall=(1-gauweigh)*curveyall !Scale Lorentzian function part of Pseudo-Voigt
                indcurve=(1-gauweigh)*indcurve
            end if
            do imol=1,nsystem
                do idata=1,numdataall(imol)
                    gauss_c=FWHMall(imol,idata)/2D0/sqrt(2*dlog(2D0))
                    gauss_a=strall(imol,idata)/(gauss_c*sqrt(2D0*pi))
                    do ipoint=1,num1Dpoints
                        curveytmp(ipoint)=gauss_a*dexp( -(curvex(ipoint)-dataxall(imol,idata))**2/(2*gauss_c**2) )
                    end do
                    if (ibroadfunc==3) curveytmp=gauweigh*curveytmp !Scale Gaussian function part of Pseudo-Voigt
                    curveyall(imol,:)=curveyall(imol,:)+curveytmp
                    if (isel==15) then !Individual contribution
                        if (idx2indband(idata)/=0) indcurve(:,idx2indband(idata))=indcurve(:,idx2indband(idata))+curveytmp
                    end if
                end do
            end do
        end if

        !Change units, scale curve, set axis
        if (ispectrum==1.and.iunitliney==2) lineyall=lineyall/2.5066D0 !For IR spectrum, convert strength unit from km/mol to esu^2*cm^2
        if (iunitx==2) then !eV->nm
            linexall=1240.7011D0/linexall
            curvex=1240.7011D0/curvex
            xlow=1240.7011D0/xlow
            xhigh=1240.7011D0/xhigh
        end if
        curveyall=scalecurve*curveyall
        if (isel==15) indcurve=scalecurve*indcurve
        curvex=curvex+shiftx
        linexall=linexall+shiftx
        if (iusersetx==0) stepx=(xhigh-xlow)/10
        
        !Generate weighted curve from multiple curves. Height of discrete lines have already been weighted when generating them
        curvey=0
        do imol=1,nsystem
            curvey=curvey+weight(imol)*curveyall(imol,:)
        end do
        !Weighting spectrum of each system
        if (iweisyscurve==1) then
            do imol=1,nsystem
                curveyall(imol,:)=weight(imol)*curveyall(imol,:)
            end do
        end if
    end if

    !================================================
    !======= Output data to external text file ======
    !================================================
    if (isel==2.or.isel==15) then !Output curve for total and individual contributions, respectively
        if (isel==2) then !Output regular spectrum
            if (nsystem==1) then
                open(10,file="spectrum_curve.txt",status="replace")
                do ipt=1,num1Dpoints
                    write(10,"(2f13.5)") curvex(ipt),curvey(ipt)
                end do
                close(10)
                write(*,*) "Curve data has been written to spectrum_curve.txt in current folder"
            else !Also output curve for all systems
                if (any(weight/=1)) then !Output weighted spectrum 
                    open(10,file="spectrum_curve.txt",status="replace")
                    do ipt=1,num1Dpoints
                        write(10,"(2f13.5)") curvex(ipt),curvey(ipt)
                    end do
                    close(10)
                    write(*,"(a)") " The curve data corresponding to weighted spectrum has been written to spectrum_curve.txt in current folder"
                end if
                open(10,file="curveall.txt",status="replace")
                do ipt=1,num1Dpoints
                    write(10,"(f13.5)",advance="no") curvex(ipt)
                    do imol=1,nsystem
                        write(10,"(f13.5)",advance="no") curveyall(imol,ipt)
                    end do
                    write(10,*)
                end do
                close(10)
                write(*,"(a)") " Curve data of all systems have been exported to curveall.txt in current folder as different columns"
            end if
        else if (isel==15) then !Output individual band contributions
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
        !Explain meaning of each column
        if (ispectrum==1) then !IR
            write(*,*) "Column 1: Wavenumber (cm^-1)"
            write(*,*) "Column 2: Molar absorption coefficient (L/mol/cm)"
        else if (ispectrum==2) then !Raman
            write(*,*) "Column 1: Wavenumber (cm^-1)"
            write(*,*) "Column 2: Relative Raman intensity"
        else if (ispectrum==3.or.ispectrum==4) then !UV-Vis, ECD
            if (iunitx==1) write(*,*) "Column 1: Excitation energy (eV)"
            if (iunitx==2) write(*,*) "Column 1: Wavelength (nm)"
            if (iunitx==3) write(*,*) "Column 1: Wavenumber (1000cm^-1)"
            if (ispectrum==3) write(*,*) "Column 2: Molar absorption coefficient (L/mol/cm)"
            if (ispectrum==4) write(*,*) "Column 2: Delta molar absorption coefficient (L/mol/cm)"
        else if (ispectrum==5) then !VCD
            write(*,*) "Column 1: Wavenumber (cm^-1)"
            write(*,*) "Column 2: Delta molar absorption coefficient (L/mol/cm)"
        end if
        if (isel==15) then
            write(*,*) "Correspondence between the columns and individual transitions in the file:"
            write(*,*) "     Column#   Transition#"
            do itmp=1,numindband
                write(*,"(2i12)") itmp+2,indband2idx(itmp)
            end do
            deallocate(indband2idx,idx2indband,indcurve)
        end if
        
        !Output discrete lines
        open(10,file="spectrum_line.txt",status="replace")
        do imol=1,nsystem
            do ipt=1,3*numdataall(imol)
                write(10,"(2f16.5)") linexall(imol,ipt),lineyall(imol,ipt)
            end do
            write(10,*)
        end do
        close(10)
        write(*,*)
        if (nsystem>1) then
            write(*,"(a)") " Discrete line data of various systems have been written together to spectrum_line.txt in current folder"
            if (any(weight/=1)) write(*,*) "Note: The height of discrete lines in this file have been weighted"
        else
            write(*,*) "Discrete line data have been written to spectrum_line.txt in current folder"
        end if
        if (ispectrum==1) then !IR
            write(*,*) "Column 1: Frequency (cm^-1)"
            if (iunitliney==1) write(*,*) "Column 2: IR intensities (km/mol)"
            if (iunitliney==2) write(*,*) "Column 2: IR intensities (esu^2*cm^2)"
        else if (ispectrum==2) then !Raman
            write(*,*) "Column 1: Frequency (cm^-1)"
            if (iramantype==1) write(*,*) "Column 2: Raman scattering activities (A^4/amu)"
            if (iramantype==2) write(*,*) "Column 2: Raman scattering intensities"
        else if (ispectrum==3.or.ispectrum==4) then !UV-Vis, ECD
            if (iunitx==1) write(*,*) "Column 1: Excitation energy (eV)"
            if (iunitx==2) write(*,*) "Column 1: Wavelength (nm)"
            if (iunitx==3) write(*,*) "Column 1: Wavenumber (1000cm^-1)"
            if (ispectrum==3) write(*,*) "Column 2: Oscillator strength"
            if (ispectrum==4) write(*,*) "Column 2: Rotatory strength in cgs (10^-40 erg-esu-cm/Gauss)"
        else if (ispectrum==5) then !VCD
            write(*,*) "Column 1: Frequency (cm^-1)"
            write(*,*) "Column 2: Rotatory strength (10^-44 esu^2 cm^2)"
        end if
        
    else if (isel==16) then !Find minimum/maximum positions
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
        if (nsystem>1) write(*,*) "Note: The minimum and maximum reported above correspond to weighted spectrum"
    end if

    !========================================
    !============ Draw spectrum =============
    !========================================
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

        if (iusersetY1==0) then !Set default lower and upper limit of left Y axis
            endy1=1.1D0*max(maxval(abs(curvey)),maxval(abs(curveyall)))
            if (nsystem>1.and.all(weight==1)) endy1=1.1D0*maxval(abs(curveyall))
            orgy1=-endy1/30D0 !Slightly lower it to avoid curve touching bottom
            if (ispectrum==4.or.ispectrum==5) orgy1=-endy1 !Positive and negative of ECD curve may have similar height
            stepy1=(endy1-orgy1)/10
        end if
        if (iusersetY2==0) then !Set default lower and upper limit of right Y axis
            endy2=1.1D0*maxval(abs(lineyall))
            orgy2=-endy2/30D0
            if (ispectrum==4.or.ispectrum==5) orgy2=-endy2
            stepy2=(endy2-orgy2)/10
        end if
        if (ishowline==1) then
        else
        end if
        ileg=0
        numleg=1+nsystem
        
        if (ishowgrid==1) then !Draw shallow gray dashed lines
        end if
        !Draw weighted curve
        if (ishowtotal==1) then
            ileg=ileg+1
        end if
        !Draw curve for each system
        if (nsystem>1) then
            do imol=1,nsystem
                iclrtmp=imol+1 !The 1st color is red, which has been used by weighted curve
                if (iclrtmp==4) iclrtmp=10 !4 corresponds to black, however, due to reverse, it corresponds to white, which is unable to use
                if (iclrtmp==2) then !2 corresponds to green, which is too bright
                    iclrtmp=12 !Change to dark green
                else if (iclrtmp==12) then
                    iclrtmp=2 
                end if
                if (weight(imol)==1D0) then
                    write(c200tmp,"(a)") trim(mollegend(imol))
                else
                    write(c200tmp,"(i3,' (',f5.1,'%)')") imol,weight(imol)*100
                end if
                ileg=ileg+1
            end do
        end if
        
        if (ishowline==1) then
            if (ispectrum==1) then
            else if (ispectrum==2) then
            else if (ispectrum==3) then
            else if (ispectrum==4) then
            else if (ispectrum==5) then
            end if
            do imol=1,nsystem
                if (iweisyscurve==1) then
                    iclrtmp=imol+1 !The 1st color is red, which has been used by weighted curve
                    if (iclrtmp==4) iclrtmp=10 !4 corresponds to black, however, due to reverse, it corresponds to white, which is unable to use
                    if (iclrtmp==2) then !2 corresponds to green, which is too bright
                        iclrtmp=12 !Change to dark green
                    else if (iclrtmp==12) then
                        iclrtmp=2 
                    end if
                end if
            end do
        end if
        if (isavepic==1) write(*,*) "Graphic file has been saved to current folder with ""DISLIN"" prefix"
    end if
    idraw=0

end do
end subroutine




!!----- Load transition data from a file to datax,str,FWHM in global memory
! ispectrum decides spectrum type, "filename" indicate the file to be loaded, numdata is the number of transitions in the file
! For electronic spectra, transition energies are loaded as eV
subroutine loadtransdata(ispectrum,loadspecname,numdata)
use defvar
use util
implicit real*8 (a-h,o-z)
character(len=*) loadspecname
character ctest,c80tmp*80,c200tmp*200
integer ispectrum
integer :: nrdfreq=0 ! >0 means pre-resonance raman, which loads external field frequency
real*8,allocatable :: rdfreq(:)

if (allocated(datax)) deallocate(datax,str,FWHM)
open(10,file=loadspecname,status="old")

!sTDA output file
if (index(filename,"tda.dat")/=0) then
    write(*,*) "Recognized as sTDA program output file"
    if (ispectrum==4) then
        write(*,*)
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
    allocate(datax(numdata),str(numdata),FWHM(numdata))
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
            if (irotrep==3) then
                str(i)=Rv
                if (fv/=0) str(i)=Rv*fl/fv
            end if
        end if
    end do

!Gaussian/ORCA output file or plain text file, all of them may use .out or .log as suffix
else
    call loclabel(10,"Gaussian, Inc",igauout)
    rewind(10)
    if (igauout==1) then
        write(*,*) "Recognized as a Gaussian output file"

        !IR, Raman, VCD
        if (ispectrum==1.or.ispectrum==2.or.ispectrum==5) then
            !Detect if this is a pre-resonance raman task, and how many frequencies are readed
            if (ispectrum==2) then
                call loclabel(10,"NFrqRd=",ifound,0)
                if (ifound==1) then
                    read(10,"(a)") c200tmp
                    itmp=index(c200tmp,"NFrqRd=")
                    read(c200tmp(itmp+7:),*) nrdfreq
                    if (nrdfreq>0) then
                        allocate(rdfreq(nrdfreq))
                        do itmp=0,nrdfreq,5
                            nleft=nrdfreq-itmp
                            read(10,"(a)") c200tmp
                            if (nleft>5) then
                                read(c200tmp(14:),*) rdfreq(itmp+1:itmp+5)
                            else
                                read(c200tmp(14:),*) rdfreq(itmp+1:nrdfreq)
                            end if
                        end do
                        write(*,*) "This is a pre-resonance Raman calculation, external field frequencies (a.u.):"
                        do itmp=1,nrdfreq
                            write(*,"(i5,':',f16.8)") itmp,rdfreq(itmp)
                        end do
                        write(*,*) "Load data for which frequency? Input its index, e.g. 3"
                        read(*,*) irdfreq
                    end if
                end if
                rewind(10)
            end if
            
            !Find how many frequencies in the file
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

            allocate(datax(numdata),str(numdata),FWHM(numdata))
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
                if (ispectrum==1) then
                    call loclabel(10,"IR Inten    --",ifound,0)
                else if (ispectrum==2) then
                    if (nrdfreq==0) then !Normal raman
                        call loclabel(10,"Raman Activ --",ifound,0)
                    else !Pre-resonance raman
                        write(c200tmp,"('RamAct Fr=',i2)") irdfreq
                        call loclabel(10,trim(c200tmp),ifound,0)
                    end if
                else if (ispectrum==5) then
                    call loclabel(10,"Rot. str.",ifound,0)            
                end if
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
                    write(*,"(a)") " Note: Found anharmonic IR information, if load them instead of the harmonic ones? (y/n)"
                    read(*,*) ctest
                    if (ctest=='y'.or.ctest=='Y') then
                        call loclabel(10,"Mode(",ifound,0)
                        read(10,*)
                        do itmp=1,numdata
                            read(10,*) c200tmp,rnouse,datax(itmp),rnouse,str(itmp)
                        end do
                    end if
                end if
            else if (ispectrum==2) then
                call loclabel(10,"Anharmonic Raman Spectroscopy",ifound,0)
                if (ifound==1) then
                    write(*,"(a)") " Note: Found anharmonic Raman information, if load them instead of the harmonic ones? (y/n)"
                    read(*,*) ctest
                    if (ctest=='y'.or.ctest=='Y') then
                        call loclabel(10,"Mode(",ifound,0)
                        read(10,*)
                        do itmp=1,numdata
                            read(10,*) c200tmp,harmfreq,datax(itmp),harmact,str(itmp)
                            str(itmp)=0.059320323D0*harmfreq*str(itmp) !The conversion coefficient can be found in output file
                        end do
                    end if
                end if
            else if (ispectrum==5) then
                call loclabel(10,"Anharmonic VCD Spectroscopy",ifound,0)
                if (ifound==1) then
                    write(*,"(a)") " Note: Found anharmonic VCD information, if load them instead of the harmonic ones? (y/n)"
                    read(*,*) ctest
                    if (ctest=='y'.or.ctest=='Y') then
                        call loclabel(10,"Mode(",ifound,0)
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
            allocate(datax(numdata),str(numdata),FWHM(numdata))
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
        call loclabel(10,"O   R   C   A",iORCAout,maxline=100)
        if (iORCAout==1) then !ORCA output file
            write(*,*) "Recognized as an ORCA output file"
            isTDA=0
            if (ispectrum==3.or.ispectrum==4) call loclabel(10,"ORCA sTD",isTDA) !When plotting UV-Vis or ECD, check if this is a sTDA or sTD-DFT calculation
            if (isTDA==0) then !Regular calculation
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
                    allocate(datax(numdata),str(numdata),FWHM(numdata))
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
                    allocate(datax(numdata),str(numdata),FWHM(numdata))
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
            else if (isTDA==1) then !sTDA or sTD-DFT calculation
                write(*,*) "This is a sTDA or sTD-DFT calculation"
                if (ispectrum==4) then
                    write(*,*)
                    write(*,*) "Read the rotatory strengths in which representation?"
                    write(*,*) "1: Length representation     2: Velocity representation"
                    write(*,*) "3: Mixed-form representation (recommended)"
                    read(*,*) irotrep
                end if
                call loclabel(10,"roots found,",ifound,0)
                read(10,*) numdata
                allocate(datax(numdata),str(numdata),FWHM(numdata))
                call loclabel(10,"state   eV        nm        fL",ifound,0)
                read(10,*)
                do i=1,numdata
                    read(10,*) inouse,datax(i),rnouse,fl,fv,Rl,Rv
                    if (ispectrum==3) then
                        str(i)=fl
                    else if (ispectrum==4) then
                        if (irotrep==1) str(i)=Rl
                        if (irotrep==2) str(i)=Rv
                        if (irotrep==3) then
                            str(i)=Rv
                            if (fv/=0) str(i)=Rv*fl/fv
                        end if
                    end if
                end do
                if (ispectrum==3) FWHM=2D0/3D0
                if (ispectrum==4) FWHM=0.2D0
            end if
        
        !Plain text file    
        else
            write(*,*) "Recognized as a plain text file"
            rewind(10)
            read(10,*) numdata,inptype
            allocate(datax(numdata),str(numdata),FWHM(numdata))
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

end subroutine
