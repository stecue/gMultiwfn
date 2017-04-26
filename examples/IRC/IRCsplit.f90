program IRCsplit
implicit real*8 (a-h,o-z)
logical alive
integer totlinenum
character IRCoutfile*200,stdinputfile*200,c79tmp*79,newgjffile*200,newwfnfile*200,TSfilename*200
character,allocatable :: stdincontent(:)*79
character*2 :: ind2name(0:109)=(/ "Bq","H ","He", &   !X(number O) is ghost atom
"Li","Be","B ","C ","N ","O ","F ","Ne", & !3~10
"Na","Mg","Al","Si","P ","S ","Cl","Ar", & !11~18
"K ","Ca","Sc","Ti","V ","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr", & !19~36
"Rb","Sr","Y ","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I ","Xe", & !37~54
"Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu", & !55~71
"Hf","Ta","W ","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn", & !72~86
"Fr","Ra","Ac","Th","Pa","U ","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr", & !87~103
"Rf","Db","Sg","Bh","Hs","Mt" /) !104~109

write(*,*) "IRCsplit v1.0.3: Generate .wfn/.chk file for each point of Gaussian09 IRC task"
write(*,*) "Programmed by Tian Lu, release date: 2015-Apr-23"
write(*,*) "Beijing Kein Research Center for Natural Sciences"
write(*,*) "E-mail: Sobereva@sina.com"
write(*,*)
write(*,*) "Input the path of output file of Gaussian IRC task, e.g. c:\miku.out"
do while(.true.)
	read(*,"(a)") IRCoutfile
! 	IRCoutfile="examples\trimerization.out"
	inquire(file=IRCoutfile,exist=alive)
	if (alive) exit
	write(*,*) "Cannot find the file, input again"
end do

open(10,file=IRCoutfile,status="old")

write(*,*)
write(*,*) "Input the path of Gaussian input file of standard single point task"
write(*,*) "e.g. c:\Love_live\niconiconi.gjf"
do while(.true.) 
	read(*,"(a)") stdinputfile
! 	stdinputfile="examples\trimerization_SP.gjf"
	inquire(file=stdinputfile,exist=alive)
	if (alive) exit
	write(*,*) "Cannot find the file, input again"
end do
open(11,file=stdinputfile,status="old")
totlinenum=0 !Count total number of lines
do while(.true.)
	read(11,"(a)",iostat=ierror) c79tmp
	if (ierror/=0) exit
	totlinenum=totlinenum+1
end do
rewind(11)
allocate(stdincontent(totlinenum))
nblank=0
iatmend=totlinenum
do iline=1,totlinenum !Load all lines into memory
	read(11,"(a)") stdincontent(iline)
	if (stdincontent(iline)==" ") then
		nblank=nblank+1
		if (nblank==2) iatmbegin=iline+2 !From which line to which line record atomic coordinates
		if (nblank==3) iatmend=iline-1
	end if
end do
close(11)

write(*,*)
write(*,*) "1: Generate .wfn   2: Generate .chk   3: Generate both .wfn and .chk"
read(*,*) itypegen
write(*,*) "Input the path of the files to be generated"
if (itypegen==1) write(*,"(a)") " If you input e.g. c:\sob, then after you run all generated Gaussian input files, c:\sob0001.wfn, c:\sob0002.wfn... will be yielded"
if (itypegen==2) write(*,"(a)") " If you input e.g. c:\sob, then after you run all generated Gaussian input files, c:\sob0001.chk, c:\sob0002.chk... will be yielded"
if (itypegen==3) write(*,"(a)") " If you input e.g. c:\sob, then after you run all generated Gaussian input files,&
 c:\sob0001.wfn, c:\sob0002.wfn... and c:\sob0001.chk, c:\sob0002.chk... will be yielded"
read(*,"(a)") newwfnfile
! newwfnfile="c:\IRC\IRC"

!Count the total number of forward and reverse path points
rewind(10)
nforwardall=0
do while(.true.)
	call loclabel(10,"Path Number:   1",ifound,0)
	if (ifound==1) then
		nforwardall=nforwardall+1
		read(10,*)
	else
		exit
	end if
end do
nforwardall=nforwardall-1 !rip out TS point, because TS points is considered as point 0 in forward path
nreverseall=0
rewind(10)
do while(.true.)
	call loclabel(10,"Path Number:   2",ifound,0)
	if (ifound==1) then
		nreverseall=nreverseall+1
		read(10,*)
	else
		exit
	end if
end do
write(*,"(' The total number of forward and reverse points:',2i6)") nforwardall,nreverseall

write(*,*)
write(*,*) "Extract how many forward and reverse points, respectively? e.g. 10,15"
read(*,*) nforward,nreverse
npttot=nforward+nreverse+1

call loclabel(10,"Input orientation:",icoortype,1) !If icoortype=1, will then load "Input orientation:"; else means that "geom=check" is used and thus load "Z-Matrix orientation:"
call loclabel(10,"NAtoms=",ifound,1)
read(10,*) c79tmp,natm
rewind(10)

do ipt=1,npttot
	if (ipt==1) then
		idxtmp=nreverse+1
	else if (ipt<=nforward+1) then !Load forward points
		idxtmp=nreverse+ipt
	else !Load reverse points
		idxtmp=nreverse+1-(ipt-nforward-1)
	end if
	nnamelen=len_trim(stdinputfile)
	write(newgjffile,"(a,i4.4,'.gjf')") stdinputfile(1:nnamelen-4),idxtmp
	write(*,"(' Generating ',a)") trim(newgjffile)
	if (ipt==1) TSfilename=newgjffile
	open(20,file=newgjffile,status="replace")
	!Locate the position of atomic coordinates in IRC output file
	if (ipt==1) then !The first point, namely TS
		call loclabel(10,"Path Number:   1",ifound,0)
		if (icoortype==1) call loclabelup(10,"Input orientation:",ifound,0) !If not found, it is highly possible that rcfc or readfc is used
		if (icoortype==0.or.ifound==0) then !if geom=check is used, then the input file corresponding to TS geometry will be skipped
			close(20)
			write(*,"(a)") " Warning: The TS geometry cannot be found! Please manually fill the TS geometry into the corresponding .gjf file"
			cycle
		end if
	else if (ipt<=nforward+1) then !FORWARD
		call loclabel(10,"Path Number:   1",ifound,0) !Pass the previous "Path Number:   1" line
		read(10,*)
		call loclabel(10,"Path Number:   1",ifound,0)
		if (icoortype==1) call loclabelup(10,"Input orientation:",ifound,0)
		if (icoortype==0) call loclabelup(10,"Z-Matrix orientation:",ifound,0)
	else !REVERSE
		if (ipt/=nforward+2) then
			call loclabel(10,"Path Number:   2",ifound,0) !Pass the previous "Path Number:   2" line
			read(10,*)
		end if
		call loclabel(10,"Path Number:   2",ifound,0)
		if (icoortype==1) call loclabelup(10,"Input orientation:",ifound,0)
		if (icoortype==0) call loclabelup(10,"Z-Matrix orientation:",ifound,0)
	end if
	read(10,*)
	read(10,*)
	read(10,*)
	read(10,*)
	read(10,*)
	if (itypegen==2.or.itypegen==3)  write(20,"('%chk=',a,i4.4,'.chk')") trim(newwfnfile),idxtmp
	do iline=1,iatmbegin-1 !Write the head part
		if (index(stdincontent(iline),'#')/=0.and.(itypegen==1.or.itypegen==3)) then
			write(20,"(a)") trim(stdincontent(iline))//" out=wfn"
		else
			write(20,"(a)") stdincontent(iline)
		end if
	end do
	do iatm=1,natm !Write Geometry
		read(10,*) inouse,iele,inouse,x,y,z
		write(20,"(a,2x,3f12.8)") ind2name(iele),x,y,z
	end do
	do iline=iatmend+1,totlinenum !Write remaining part
		write(20,"(a)") stdincontent(iline)
	end do
	if (itypegen==1.or.itypegen==3) write(20,"(a,i4.4,a,/)") trim(newwfnfile),idxtmp,".wfn"
	close(20)
end do

close(10)
write(*,*)
write(*,"(a)") " All Gaussian input files have been generated current folder, run them by Gaussian then you will obtain corresponding .wfn/.chk files"
if (icoortype==1) then
	write(*,"(' Note: ',a,' corresponds to TS geometry')") trim(TSfilename)
else if (icoortype==0) then
	write(*,"(' Note: ',a,' is an empty file, which should correspond to TS geometry')") trim(TSfilename)
end if
pause
end program



!!!-------- Find the line where the label first appears in fileid
!Return ifound=1 if found the label, else return 0
!irewind=1 will rewind, =0 will not rewind
!If current line already has the label, calling this subroutine will do nothing
subroutine loclabel(fileid,label,ifound,irewind)
integer fileid,ierror,ifound,irewind
character*200 c200
CHARACTER(LEN=*) label
if (irewind==1) rewind(fileid)
do while(.true.)
	read(fileid,"(a)",iostat=ierror) c200
	if (index(c200,label)/=0) then
		backspace(fileid)
		ifound=1 ! Found result
		return
	end if
	if (ierror/=0) exit
end do
ifound=0
end subroutine

!Loclabel, but in reverse direction
subroutine loclabelup(fileid,label,ifound,irewind)
integer fileid,ierror,ifound,irewind
character*200 c200
CHARACTER(LEN=*) label
if (irewind==1) rewind(fileid)
do while(.true.)
	read(fileid,"(a)",iostat=ierror) c200
	if (index(c200,label)/=0) then
		backspace(fileid)
		ifound=1 ! Found result
		return
	else if (index(c200,"Copyright (c)")/=0) then !Up enough
		ifound=0
		return
	else
		backspace(fileid)
		backspace(fileid)
	end if
	if (ierror/=0) exit
end do
ifound=0
end subroutine