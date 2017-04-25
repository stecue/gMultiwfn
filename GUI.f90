module GUI
implicit real*8(a-h,o-z)

contains

!Select input file by GUI
subroutine selfilegui
CALL dwgfil("Choose a Fortran file",filename,"*.wfn")
end subroutine

!!--------- A GUI for drawing molecular structure and orbital isosurface
subroutine drawmolgui
character ictmp*4,molorblist*50000 !max 9999 orbitals (the 0th is "none"), each one take up 4 characters, adding "|",so 10000*(4+1)=50000
integer :: nmothres=144 !If nmo > this value, the list will be long and WGLIS will crash, so use text box instead
! Set variables for viewing orbital
molorblist(1:4)="None"
molorblist(5:)=" "
do i=1,nmo
	write(ictmp,"(i4)") i
	molorblist(i*5+1:i*5+5)="|"//ictmp
end do

GUI_mode=1
idrawmol=1 !Molecular structure must be drawn
CALL swgtit('Molecular structure')
call swgwth(90)
CALL SWGOPT("CENTER","POSITION") !Main window appear in the center of screen
call SWGPOP("NOOK")  !Don't show OK&QUIT&HELP in upper menu
call SWGPOP("NOQUIT")
call SWGPOP("NOHELP")
CALL WGINI('HORI',idiswindow)
call swgatt(idiswindow,"INACTIVE","CLOSE") !Disable close button
CALL WGPOP(idiswindow,"Orbital info.",idisorbinfomenu)
call wgapp(idisorbinfomenu,"Show all",idisorbinfo)
if ((wfntype==0.or.wfntype==1.or.wfntype==2).and.(ifiletype==1.or.ifiletype==9)) call wgapp(idisorbinfomenu,"Show up to LUMO+10",idisorbinfo2)
CALL WGPOP(idiswindow," Isosur#1 style",idisisosur1style)
call wgapp(idisisosur1style,"Use solid face",idisisosur1solid)
call wgapp(idisisosur1style,"Use mesh",idisisosur1mesh)
call wgapp(idisisosur1style,"Use points",idisisosur1point)
call wgapp(idisisosur1style,"Use solid face+mesh",idisisosur1solidmesh)
call wgapp(idisisosur1style,"Use transparent face",idisisosur1tpr)
call wgapp(idisisosur1style,"Set color for face",idisisosur1solidclr)
call wgapp(idisisosur1style,"Set color for mesh and points",idisisosur1meshptclr)
call wgapp(idisisosur1style,"Set opacity for transparent face",idisisosur1opa)
CALL WGPOP(idiswindow," Isosur#2 style",idisisosur2style)
call wgapp(idisisosur2style,"Use solid face",idisisosur2solid)
call wgapp(idisisosur2style,"Use mesh",idisisosur2mesh)
call wgapp(idisisosur2style,"Use points",idisisosur2point)
call wgapp(idisisosur2style,"Use solid face+mesh",idisisosur2solidmesh)
call wgapp(idisisosur2style,"Use transparent face",idisisosur2tpr)
call wgapp(idisisosur2style,"Set color for face",idisisosur2solidclr)
call wgapp(idisisosur2style,"Set color for mesh and points",idisisosur2meshptclr)
call wgapp(idisisosur2style,"Set opacity for transparent face",idisisosur2opa)
CALL WGPOP(idiswindow," Isosur. quality",idisisosurquality)
call wgapp(idisisosurquality,"Set number of grid points",idisisosurnumpt)
CALL WGPOP(idiswindow," Set lighting",idissetlight)
call wgapp(idissetlight,"Light 1",idissetlight1)
call wgapp(idissetlight,"Light 2",idissetlight2)
call wgapp(idissetlight,"Light 3",idissetlight3)
call wgapp(idissetlight,"Light 4",idissetlight4)
call wgapp(idissetlight,"Light 5",idissetlight5)
call wgapp(idissetlight,"Light 6",idissetlight6)
call wgapp(idissetlight,"Light 7",idissetlight7)
call wgapp(idissetlight,"Light 8",idissetlight8)
call wgapp(idissetlight,"Enable all",idissetlightall1)
call wgapp(idissetlight,"Disable all",idissetlightall0)

CALL WGDRAW(idiswindow,idisgraph) !Draw-widget to display molecular structure
CALL SWGWTH(20) !Set parent widget width
CALL SWGSPC(1.3D0,0.0D0) !Set space between widgets below
CALL WGBAS(idiswindow,"VERT",idisright)
if ((isys==1.and.imodlayout==1).or.isys==2.or.isys==3) CALL WGBAS(idiswindow,"VERT",idisright2) !Provide another frame for linux version
call wgpbut(idisright,"RETURN",idisreturn)
call wgpbut(idisright,"Up",idisrotup)
call wgpbut(idisright,"Down",idisrotdown)
call wgpbut(idisright,"Left",idisrotleft)
call wgpbut(idisright,"Right",idisrotright)
call wgpbut(idisright,"Zoom in",idiszoomin)
call wgpbut(idisright,"Zoom out",idiszoomout)
call wgpbut(idisright,"Reset view",idisreset)
call wgpbut(idisright,"Save picture",idissavepic)
call wgbut(idisright,"Show labels",ishowatmlab,idisshowatmlab)
call wgbut(idisright,"Show axis",ishowaxis,idisshowaxis)
call wgbut(idisright,"Show+Sel. isosur#2",isosursec,idisisosursec)
call swgatt(idisisosursec,"INACTIVE","STATUS") !User cannot select isosurface 2 when just entered this GUI, it must be actived by selecting an orbital
call SWGSTP(0.05D0)
call wgscl(idisright,"Bonding threshold",0.0D0,5.0D0,1.15D0,2,idisbondcrit)
call wgscl(idisright,"Ratio of atomic size",0.0D0,5.0D0,1.0D0,2,idisatmsize)
call SWGSTP(0.02D0)
call wgscl(idisright,"Radius of bonds",0.0D0,2.0D0,0.2D0,2,idisbondradius)
call SWGSTP(2.0D0)
call wgscl(idisright,"Size of atomic labels",0.0D0,200.0D0,30.0D0,0,idislabelsize)
CALL SWGSPC(0.0D0,0.0D0)
sur_value=0.05D0
if (isys==1.and.imodlayout==0) then
	!Set region for orbital viewing
	call WGLAB(idisright,"Orbitals:",iorbseltext)
	CALL WGBAS(idisright,"FORM",idisbotrig)
	!Set orbital selection list
	call swgwin(0,5,70,120)
	call swgtyp("VSCROLL","LIST") !This is neccessary, else some latter orbitals cannot be displayed on the list
	CALL WGLIS(idisbotrig,molorblist,1,iorbplot)
	!Set progress bar
	call swgopt("SMOOTH","PBAR")
	call swgclr(0D0,0D0,1D0,"PBAR")
	call swgtyp("VERT","PBAR")
	call swgwin(78,5,15,115)
	call wgpbar(idisbotrig,0.0D0,1.0D0,0.1D0,iprogbar)
	!Set scale bar of isosurvalue
	call swgwin(100,5,70,115)
	call swgtyp("VERT","SCALE")
	call swgstp(0.005D0)
	call wgscl(idisbotrig,"Isovalue",0D0,0.5D0,sur_value,3,idisisosurscl)
else if ((isys==1.and.imodlayout==1).or.isys==2.or.isys==3) then !Use different layout for linux, since the sizes of widgets relative to Windows version are different
	CALL SWGSPC(0.0D0,0.5D0)
	if ((isys==1.and.imodlayout==1).or.nmo<nmothres) then
		call WGLAB(idisright2,"View orbitals:",iorbseltext)
		call swgtyp("SCROLL","LIST")
		call WGLIS(idisright2,molorblist,1,iorbplot)
	else
		call WGLAB(idisright2,"Orbital index:",iorbseltext)
		call WGTXT(idisright2,"  0 ",iorbplot)
	end if
	call swgopt("SMOOTH","PBAR")
	call swgclr(0D0,0D0,1D0,"PBAR")
	call swgtyp("HORI","PBAR")
	call WGLAB(idisright2,"Progress:",iorbseltext)
	call wgpbar(idisright2,0.0D0,1.0D0,0.1D0,iprogbar)
	call swgtyp("HORI","SCALE")
	call swgstp(0.005D0)
	call wgscl(idisright2,"Isovalue of orbital",0D0,0.5D0,sur_value,3,idisisosurscl)
end if


call SWGCBK(idissetlight1,setlight1)
call SWGCBK(idissetlight2,setlight2)
call SWGCBK(idissetlight3,setlight3)
call SWGCBK(idissetlight4,setlight4)
call SWGCBK(idissetlight5,setlight5)
call SWGCBK(idissetlight6,setlight6)
call SWGCBK(idissetlight7,setlight7)
call SWGCBK(idissetlight8,setlight8)
call SWGCBK(idissetlightall0,setlightall0)
call SWGCBK(idissetlightall1,setlightall1)
call SWGCBK(idisisosur1solid,setisosur1solid)
call SWGCBK(idisisosur1mesh,setisosur1line)
call SWGCBK(idisisosur1point,setisosur1point)
call SWGCBK(idisisosur1solidmesh,setisosur1solidmesh)
call SWGCBK(idisisosur1tpr,setisosur1tpr)
call SWGCBK(idisisosur1solidclr,setisosur1solidclr)
call SWGCBK(idisisosur1meshptclr,setisosur1meshptclr)
call SWGCBK(idisisosur1opa,setisosur1opa)
call SWGCBK(idisisosur2solid,setisosur2solid)
call SWGCBK(idisisosur2mesh,setisosur2line)
call SWGCBK(idisisosur2point,setisosur2point)
call SWGCBK(idisisosur2solidmesh,setisosur2solidmesh)
call SWGCBK(idisisosur2tpr,setisosur2tpr)
call SWGCBK(idisisosur2solidclr,setisosur2solidclr)
call SWGCBK(idisisosur2meshptclr,setisosur2meshptclr)
call SWGCBK(idisisosur2opa,setisosur2opa)
call SWGCBK(idisisosurnumpt,setisosurnumpt)
call SWGCBK(idisorbinfo,showorbinfo)
if ((wfntype==0.or.wfntype==1.or.wfntype==2).and.(ifiletype==1.or.ifiletype==9)) call SWGCBK(idisorbinfo2,showorbinfo2)
call SWGCBK(idisreturn,GUIreturn)
call SWGCBK(idisrotleft,rotleft)
call SWGCBK(idisrotright,rotright)
call SWGCBK(idisrotup,rotup)
call SWGCBK(idisrotdown,rotdown)
call SWGCBK(idiszoomin,zoomin)
call SWGCBK(idiszoomout,zoomout)
call SWGCBK(idisreset,resetview)
call SWGCBK(idissavepic,savepic)
call SWGCBK(idisshowatmlab,ifshowatmlabel)
call SWGCBK(idisshowaxis,ifshowaxis)
call SWGCBK(idisisosursec,ifisosursec)
call SWGCBK(idisbondcrit,setbondcrit)
call SWGCBK(idislabelsize,setlabelsize)
call SWGCBK(idisatmsize,setatmsize)
call SWGCBK(idisbondradius,setbondradius)
!! Click button and calculate cube data for selected orbital
if (isys==1) then
	call SWGCBK(iorbplot,showorbsellist)
else !If nmo is too large, using list widget will cause crash, in this case switch to text box input
	if (nmo<nmothres) then
		call SWGCBK(iorbplot,showorbsellist)
	else
		call SWGCBK(iorbplot,showorbselbox)
	end if
end if
call SWGCBK(idisisosurscl,setisosurscl)
CALL SWGSPC(4.0D0,0.5D0) !Reset the default widget spacing
call swgtyp("HORI","SCALE") !Reset the default mode for list widget
idrawisosur=0 !Don't draw the cubmat in memory at first time go into the GUI
!However, in linux, "draw" widget is available only after WGFIN subroutine so we need a mouse event to active it, before this, the draw widget cannot be used, this is why "if (isys==1)"
CALL WGFIN
idrawisosur=0 !After ending this GUI, recover to initial setting
isosur1style=1 !Recover to solid face
isosur2style=1
isosursec=0 !Don't draw the second isosurface
end subroutine


!!------------ A GUI for drawing relief map
subroutine drawplanegui(init1,end1,init2,end2,init3,end3,idrawtype)
real*8 init1,end1,init2,end2,init3,end3
integer,intent (in) :: idrawtype
GUI_mode=2
dp_init1=init1
dp_end1=end1
dp_init2=init2
dp_end2=end2
dp_init3=init3
dp_end3=end3
if (isavepic==1) then
else if (isavepic==0) then
	if (idrawtype==3) CALL swgtit('Relif map')
	if (idrawtype==4) CALL swgtit('Shaded surface map')
	if (idrawtype==5) CALL swgtit('Shaded surface map')
	CALL swgwth(96)
	CALL SWGOPT("CENTER","POSITION") !Main window appear in the center of screen
	call SWGPOP("NOOK")  !Don't show OK&QUIT in upper menu
	call SWGPOP("NOQUIT")
	call SWGPOP("NOHELP")
	CALL WGINI('HORI',idiswindow)
	call swgatt(idiswindow,"INACTIVE","CLOSE") !Disable close button	
	CALL WGDRAW(idiswindow,idisgraph) !Draw-widget to display molecular structure
	CALL SWGWTH(20) !Set parent widget width
	CALL WGBAS(idiswindow,"VERT",idisright)
	CALL WGBAS(idisright,"VERT",idisOK)
	! CALL WGBAS(idisright,"FORM",idisOK) !seems this function has bug in current dislin version
	! call swgsiz(int(iscrwidth*0.12D0),50)
	call wgpbut(idisOK,"RETURN",idisreturn)
	call wgsep(idisright,idissep2)
	call wgpbut(idisright,"Up",idisrotup)
	call wgpbut(idisright,"Down",idisrotdown)
	call wgpbut(idisright,"Left",idisrotleft)
	call wgpbut(idisright,"Right",idisrotright)
	call wgpbut(idisright,"Zoom in",idiszoomin)
	call wgpbut(idisright,"Zoom out",idiszoomout)
	call wgpbut(idisright,"Reset view",idisreset)
	call SWGCBK(idisreturn,GUIreturn)
	call SWGCBK(idisrotleft,rotleft)
	call SWGCBK(idisrotright,rotright)
	call SWGCBK(idisrotup,rotup)
	call SWGCBK(idisrotdown,rotdown)
	call SWGCBK(idiszoomin,zoomin)
	call SWGCBK(idiszoomout,zoomout)
	call SWGCBK(idisreset,resetview)
	CALL WGFIN
end if
end subroutine


!!--------- A GUI for drawing isosurface
!if iallowsetstyle==1, then isosurface style can be customly controlled for isosurface 1
!if iallowsetstyle==2, then isosurface style can be customly controlled for both isosurface 1 and 2
subroutine drawisosurgui(iallowsetstyle)
integer iallowsetstyle
character*20 temp
idrawisosur=1
ishowattlab=1 !Show attractors, so that one can compare attractors with isosurfaces
ishowatt=1
CALL swgtit('Isosurface graph')
CALL swgwth(96)
CALL SWGOPT("CENTER","POSITION") !Main window appear in the center of screen
call SWGPOP("NOOK")  !Don't show OK&QUIT in upper menu
call SWGPOP("NOQUIT")
call SWGPOP("NOHELP")
! CALL SWGHLP("Green region: Isosurface value|Blue region: Negative of the isosurface value|&
! The min & max value of scale bar is -5 and 5 respectively, if the inputted value exceed this range, scale bar will not change")
CALL WGINI('HORI',idiswindow)
call swgatt(idiswindow,"INACTIVE","CLOSE") !Disable close button
if (iallowsetstyle==1) then
	call WGPOP(idiswindow," Isosurface style",idisisosur1style)
	call wgapp(idisisosur1style,"Use solid face",idisisosur1solid)
	call wgapp(idisisosur1style,"Use mesh",idisisosur1mesh)
	call wgapp(idisisosur1style,"Use points",idisisosur1point)
	call wgapp(idisisosur1style,"Use solid face+mesh",idisisosur1solidmesh)
	call wgapp(idisisosur1style,"Use transparent face",idisisosur1tpr)
	call wgapp(idisisosur1style,"Set color for face",idisisosur1solidclr)
	call wgapp(idisisosur1style,"Set color for mesh and points",idisisosur1meshptclr)
	call wgapp(idisisosur1style,"Set opacity for transparent face",idisisosur1opa)
else if (iallowsetstyle==2) then
	call WGPOP(idiswindow," Isosurface style",idisisosurallstyle)
	call wgapp(idisisosurallstyle,"Use solid face",idisisosurallsolid)
	call wgapp(idisisosurallstyle,"Use mesh",idisisosurallmesh)
	call wgapp(idisisosurallstyle,"Use points",idisisosurallpoint)
	call wgapp(idisisosurallstyle,"Use solid face+mesh",idisisosurallsolidmesh)
	call wgapp(idisisosurallstyle,"Use transparent face",idisisosuralltpr)
end if
CALL WGPOP(idiswindow," Set lighting",idissetlight)
call wgapp(idissetlight,"Light 1",idissetlight1)
call wgapp(idissetlight,"Light 2",idissetlight2)
call wgapp(idissetlight,"Light 3",idissetlight3)
call wgapp(idissetlight,"Light 4",idissetlight4)
call wgapp(idissetlight,"Light 5",idissetlight5)
call wgapp(idissetlight,"Light 6",idissetlight6)
call wgapp(idissetlight,"Light 7",idissetlight7)
call wgapp(idissetlight,"Light 8",idissetlight8)
call wgapp(idissetlight,"Enable all",idissetlightall1)
call wgapp(idissetlight,"Disable all",idissetlightall0)

CALL WGDRAW(idiswindow,idisgraph) !Draw-widget to display molecular structure
CALL SWGWTH(20) !Set parent widget width
CALL WGBAS(idiswindow,"VERT",idisright)
CALL WGBAS(idisright,"VERT",idisOK)
! CALL WGBAS(idisright,"FORM",idisOK) !seems this function has bug in current dislin version
! call swgsiz(int(iscrwidth*0.12D0),50)
CALL SWGSPC(1.3D0,0.0D0) !Set space between widgets below
call wgpbut(idisOK,"RETURN",idisreturn)
call wgpbut(idisright,"Up",idisrotup)
call wgpbut(idisright,"Down",idisrotdown)
call wgpbut(idisright,"Left",idisrotleft)
call wgpbut(idisright,"Right",idisrotright)
call wgpbut(idisright,"Zoom in",idiszoomin)
call wgpbut(idisright,"Zoom out",idiszoomout)
call wgpbut(idisright,"Reset view",idisreset)
call wgpbut(idisright,"Save picture",idissavepic)
write(temp,"(f10.5)") sur_value
! call WGLTXT(idisright,"Min:",temp1,70,idisscrmin) !Sadly, up to now dislin don't have routine can change scale widget min&max value
! call WGLTXT(idisright,"Max:",temp2,70,idisscrmax)
call WGLAB(idisright,"Isosurface value:",idislabel)
call WGTXT(idisright,temp,idisscrval)
if (isosursec==0) call wgbut(idisright,"Show both sign",isosurshowboth,idisshowbothsign) !When showing two grid data, this option is meaningless
call wgbut(idisright,"Show molecule",idrawmol,idisshowmol)
call wgbut(idisright,"Show atomic labels",ishowatmlab,idisshowatmlab)
call wgbut(idisright,"Show axis",ishowaxis,idisshowaxis)
call wgbut(idisright,"Show data range",ishowdatarange,idisshowdatarange)
call wgbut(idisright,"Show isosurface",idrawisosur,idisshowisosur)
call swgstp(0.01D0) !Use smaller step size of scale bar than default
if (sur_value>5) then !Do not let sur_value exceed axis range
	call wgscl(idisright,"Isosurface value",-5D0,5D0,5D0,4,idisisosurscl)
else if (sur_value<-5) then
	call wgscl(idisright,"Isosurface value",-5D0,5D0,-5D0,4,idisisosurscl)
else
	call wgscl(idisright,"Isosurface value",-5D0,5D0,sur_value,4,idisisosurscl)
end if
call SWGSTP(0.05D0)
call wgscl(idisright,"Bonding threshold",0.0D0,5.0D0,1.15D0,2,idisbondcrit)
call wgscl(idisright,"Ratio of atomic size",0.0D0,5.0D0,1.0D0,2,idisatmsize)
! call SWGSTP(0.02D0)
! call wgscl(idisright,"Radius of bonds",0.0D0,2.0D0,0.2D0,2,idisbondradius) !There no enough vertical space for this widget in Linux platform
call SWGSTP(2.0D0)
call wgscl(idisright,"Size of atomic labels",0.0D0,200.0D0,30.0D0,0,idislabelsize)

if (iallowsetstyle==1) then
	call SWGCBK(idisisosur1solid,setisosur1solid)
	call SWGCBK(idisisosur1mesh,setisosur1line)
	call SWGCBK(idisisosur1point,setisosur1point)
	call SWGCBK(idisisosur1solidmesh,setisosur1solidmesh)
	call SWGCBK(idisisosur1tpr,setisosur1tpr)
	call SWGCBK(idisisosur1solidclr,setisosur1solidclr)
	call SWGCBK(idisisosur1meshptclr,setisosur1meshptclr)
	call SWGCBK(idisisosur1opa,setisosur1opa)
else if (iallowsetstyle==2) then
	call SWGCBK(idisisosurallsolid,setisosurallsolid)
	call SWGCBK(idisisosurallmesh,setisosurallline)
	call SWGCBK(idisisosurallpoint,setisosurallpoint)
	call SWGCBK(idisisosurallsolidmesh,setisosurallsolidmesh)
	call SWGCBK(idisisosuralltpr,setisosuralltpr)
end if
call SWGCBK(idissetlight1,setlight1)
call SWGCBK(idissetlight2,setlight2)
call SWGCBK(idissetlight3,setlight3)
call SWGCBK(idissetlight4,setlight4)
call SWGCBK(idissetlight5,setlight5)
call SWGCBK(idissetlight6,setlight6)
call SWGCBK(idissetlight7,setlight7)
call SWGCBK(idissetlight8,setlight8)
call SWGCBK(idissetlightall0,setlightall0)
call SWGCBK(idissetlightall1,setlightall1)
call SWGCBK(idisreturn,GUIreturn)
call SWGCBK(idisrotleft,rotleft)
call SWGCBK(idisrotright,rotright)
call SWGCBK(idisrotup,rotup)
call SWGCBK(idisrotdown,rotdown)
call SWGCBK(idiszoomin,zoomin)
call SWGCBK(idiszoomout,zoomout)
call SWGCBK(idisreset,resetview)
call SWGCBK(idisisosurscl,setisosurscl)
if (isosursec==0) call SWGCBK(idisshowbothsign,setshowbothsign)
call SWGCBK(idisshowmol,setshowmolstruct)
call SWGCBK(idisshowatmlab,setshowatmlab)
call SWGCBK(idisshowdatarange,setshowdatarange)
call SWGCBK(idisshowisosur,ifshowisosur)
call SWGCBK(idisshowaxis,ifshowaxis)
call SWGCBK(idissavepic,savepic)
call SWGCBK(idisbondcrit,setbondcrit)
call SWGCBK(idislabelsize,setlabelsize)
call SWGCBK(idisatmsize,setatmsize)
! call SWGCBK(idisbondradius,setbondradius)
CALL WGFIN
idrawisosur=0
isosur1style=1
isosur2style=1
isosursec=0
ishowattlab=0 !Don't show attractors in other GUIs
ishowatt=0
end subroutine


!!----------------- A GUI for drawing molecule and CPs and paths
subroutine drawmoltopogui
GUI_mode=4
CALL swgtit('Molecular structure, critical points and topology paths')
CALL swgwth(90)
CALL SWGOPT("CENTER","POSITION") !Main window appear in the center of screen
call SWGPOP("NOOK")  !Don't show OK&QUIT&HELP in upper menu
call SWGPOP("NOQUIT")
call SWGPOP("NOHELP")
CALL WGINI('HORI',idiswindow)
call swgatt(idiswindow,"INACTIVE","CLOSE") !Disable close button
CALL WGDRAW(idiswindow,idisgraph) !Draw-widget to display molecular structure
CALL SWGWTH(20) !Set parent widget width
CALL SWGSPC(1.0D0,0.0D0) !Set space between widgets below
CALL WGBAS(idiswindow,"VERT",idisright)
CALL WGBAS(idisright,"VERT",idisOK)
call wgpbut(idisOK,"RETURN",idisreturn)
call wgpbut(idisright,"Up",idisrotup)
call wgpbut(idisright,"Down",idisrotdown)
call wgpbut(idisright,"Left",idisrotleft)
call wgpbut(idisright,"Right",idisrotright)
call wgpbut(idisright,"Zoom in",idiszoomin)
call wgpbut(idisright,"Zoom out",idiszoomout)
call wgpbut(idisright,"Reset view",idisreset)
call wgpbut(idisright,"Save picture",idissavepic)
call wgbut(idisright,"Atom labels",ishowatmlab,idisshowatmlab)
call wgbut(idisright,"CP labels",ishowCPlab,idisshowCPlab)
call wgbut(idisright,"Path labels",ishowpathlab,idisshowpathlab)
call wgbut(idisright,"Paths",idrawpath,idisshowpath)
call wgbut(idisright,"Basin surface",idrawpath,idisshowbassurf)
call wgbut(idisright,"Show molecule",idrawmol,idisshowmol)
call wgbut(idisright,"Show axis",ishowaxis,idisshowaxis)
call wgbut(idisright,"Show (3,-3)",ishow3n3,idisshow3n3)
call wgbut(idisright,"Show (3,-1)",ishow3n1,idisshow3n1)
call wgbut(idisright,"Show (3,+1)",ishow3p1,idisshow3p1)
call wgbut(idisright,"Show (3,+3)",ishow3p3,idisshow3p3)
CALL SWGSPC(1.0D0,0.0D0)
call SWGSTP(0.1D0)
call wgscl(idisright,"Ratio of atomic size",0.0D0,1.5D0,ratioatmsphere,2,idisatmsize)
call SWGSTP(0.02D0)
call wgscl(idisright,"Radius of bonds",0.0D0,0.5D0,bondradius,2,idisbondradius)
call SWGSTP(2.0D0)
call wgscl(idisright,"Size of labels",0.0D0,80.0D0,textheigh,0,idislabelsize)
call SWGCBK(idisreturn,GUIreturn)
call SWGCBK(idisrotleft,rotleft)
call SWGCBK(idisrotright,rotright)
call SWGCBK(idisrotup,rotup)
call SWGCBK(idisrotdown,rotdown)
call SWGCBK(idiszoomin,zoomin)
call SWGCBK(idiszoomout,zoomout)
call SWGCBK(idisreset,resetview)
call SWGCBK(idissavepic,savepic)
call SWGCBK(idisshowatmlab,ifshowatmlabel)
call SWGCBK(idisshowCPlab,ifshowCPlabel)
call SWGCBK(idisshowpathlab,ifshowpathlabel)
call SWGCBK(idisshowaxis,ifshowaxis)
call SWGCBK(idisshowpath,setshowpath)
call SWGCBK(idisshowbassurf,setshowbassurf)
call SWGCBK(idisshowmol,setshowmolstruct)
call SWGCBK(idisshow3n3,ifshow3n3)
call SWGCBK(idisshow3n1,ifshow3n1)
call SWGCBK(idisshow3p1,ifshow3p1)
call SWGCBK(idisshow3p3,ifshow3p3)
call SWGCBK(idisatmsize,setatmsize)
call SWGCBK(idisbondradius,setbondradius)
call SWGCBK(idislabelsize,setlabelsize)
CALL SWGSPC(4.0D0,0.5D0) !Reset the default widget spacing
call swgtyp("HORI","SCALE") !Reset the default mode for list widget
CALL WGFIN
end subroutine


!!----------------- A GUI for drawing molecule and surface minima and maxima for quantitative surface analysis
subroutine drawsurfanalysis
GUI_mode=5
CALL swgtit('Molecular structure, surface minima and maxima')
CALL swgwth(90)
CALL SWGOPT("CENTER","POSITION") !Main window appear in the center of screen
call SWGPOP("NOOK")  !Don't show OK&QUIT&HELP in upper menu
call SWGPOP("NOQUIT")
call SWGPOP("NOHELP")
CALL WGINI('HORI',idiswindow)
call swgatt(idiswindow,"INACTIVE","CLOSE") !Disable close button
CALL WGDRAW(idiswindow,idisgraph) !Draw-widget to display molecular structure
CALL SWGWTH(20) !Set parent widget width
CALL SWGSPC(1.0D0,0.0D0) !Set space between widgets below
CALL WGBAS(idiswindow,"VERT",idisright)
CALL WGBAS(idisright,"VERT",idisOK)
call wgpbut(idisOK,"RETURN",idisreturn)
call wgsep(idisright,idissep1)
call wgpbut(idisright,"Up",idisrotup)
call wgpbut(idisright,"Down",idisrotdown)
call wgpbut(idisright,"Left",idisrotleft)
call wgpbut(idisright,"Right",idisrotright)
call wgpbut(idisright,"Zoom in",idiszoomin)
call wgpbut(idisright,"Zoom out",idiszoomout)
call wgpbut(idisright,"Reset view",idisreset)
call wgpbut(idisright,"Save picture",idissavepic)
call wgbut(idisright,"Atom labels",ishowatmlab,idisshowatmlab)
call wgbut(idisright,"Minimum label",ishowlocminlab,idisshowlocminlab)
call wgbut(idisright,"Maximum label",ishowlocmaxlab,idisshowlocmaxlab)
call wgbut(idisright,"Minimum position",ishowlocminpos,idisshowlocminpos)
call wgbut(idisright,"Maximum position",ishowlocmaxpos,idisshowlocmaxpos)
call wgbut(idisright,"Show axis",ishowaxis,idisshowaxis)
CALL SWGSPC(1.0D0,0.0D0)
call SWGSTP(0.1D0)
call wgscl(idisright,"Ratio of atomic size",0.0D0,6D0,ratioatmsphere,2,idisatmsize)
call SWGSTP(0.02D0)
call wgscl(idisright,"Radius of bonds",0.0D0,0.5D0,bondradius,2,idisbondradius)
call SWGSTP(2.0D0)
call wgscl(idisright,"Size of labels",0.0D0,80.0D0,textheigh,0,idislabelsize)
call SWGCBK(idisreturn,GUIreturn)
call SWGCBK(idisrotleft,rotleft)
call SWGCBK(idisrotright,rotright)
call SWGCBK(idisrotup,rotup)
call SWGCBK(idisrotdown,rotdown)
call SWGCBK(idiszoomin,zoomin)
call SWGCBK(idiszoomout,zoomout)
call SWGCBK(idisreset,resetview)
call SWGCBK(idissavepic,savepic)
call SWGCBK(idisshowatmlab,ifshowatmlabel)
call SWGCBK(idisshowlocminlab,setshowlocminlab)
call SWGCBK(idisshowlocmaxlab,setshowlocmaxlab)
call SWGCBK(idisshowlocminpos,setshowlocminpos)
call SWGCBK(idisshowlocmaxpos,setshowlocmaxpos)
call SWGCBK(idisshowaxis,ifshowaxis)
call SWGCBK(idisatmsize,setatmsize)
call SWGCBK(idisbondradius,setbondradius)
call SWGCBK(idislabelsize,setlabelsize)
CALL SWGSPC(4.0D0,0.5D0) !Reset the default widget spacing
call swgtyp("HORI","SCALE") !Reset the default mode for list widget
CALL WGFIN
end subroutine


!!----------------- A GUI for drawing basin integration visualization
subroutine drawbasinintgui
use basinintmod
character ictmp*4,basinlist*50000 !max 9999 basins (the 0th is "none"), each one take up 4 characters, adding "|",so 10000*(4+1)=50000
GUI_mode=6
! Set variables for viewing orbital
basinlist(1:4)="None"
basinlist(5:)=" "
do irealatt=1,numrealatt
	write(ictmp,"(i4)") irealatt
	basinlist(irealatt*5+1:irealatt*5+5)="|"//ictmp
end do
basinlist((numrealatt+1)*5+1:(numrealatt+1)*5+5)="|"//"Unas" !numrealatt+1 means unassigned
basinlist((numrealatt+2)*5+1:(numrealatt+2)*5+5)="|"//"Boun" !numrealatt+2 means the ones go to boundary

CALL swgtit('Molecular structure, attractors and basins')
CALL swgwth(90)
CALL SWGOPT("CENTER","POSITION") !Main window appear in the center of screen
call SWGPOP("NOOK")  !Don't show OK&QUIT&HELP in upper menu
call SWGPOP("NOQUIT")
call SWGPOP("NOHELP")
CALL WGINI('HORI',idiswindow)
call swgatt(idiswindow,"INACTIVE","CLOSE") !Disable close button
CALL WGDRAW(idiswindow,idisgraph) !Draw-widget to display molecular structure
CALL SWGWTH(20) !Set parent widget width
CALL SWGSPC(1.0D0,0.0D0) !Set space between widgets below
CALL WGBAS(idiswindow,"VERT",idisright)
CALL WGBAS(idiswindow,"VERT",idisright2) !Provide another frame for linux version
CALL WGBAS(idisright,"VERT",idisOK)
call wgpbut(idisOK,"RETURN",idisreturn)
call wgpbut(idisright,"Up",idisrotup)
call wgpbut(idisright,"Down",idisrotdown)
call wgpbut(idisright,"Left",idisrotleft)
call wgpbut(idisright,"Right",idisrotright)
call wgpbut(idisright,"Zoom in",idiszoomin)
call wgpbut(idisright,"Zoom out",idiszoomout)
call wgpbut(idisright,"Reset view",idisreset)
call wgpbut(idisright,"Save picture",idissavepic)
call wgbut(idisright,"Show molecule",idrawmol,idisshowmol)
call wgbut(idisright,"Show axis",ishowaxis,idisshowaxis)
call wgbut(idisright,"Atom labels",ishowatmlab,idisshowatmlab)
call wgbut(idisright,"Attractor labels",ishowattlab,idisshowattlab)
call wgbut(idisright,"Show basin interior",idrawinternalbasin,idisdrawinternalbasin)
CALL SWGSPC(1.0D0,0.0D0)
call SWGSTP(0.1D0)
call wgscl(idisright,"Ratio of atomic size",0.0D0,1.5D0,ratioatmsphere,2,idisatmsize)
call SWGSTP(0.02D0)
call wgscl(idisright,"Radius of bonds",0.0D0,0.5D0,bondradius,2,idisbondradius)
call SWGSTP(2.0D0)
call wgscl(idisright,"Size of labels",0.0D0,80.0D0,textheigh,0,idislabelsize)
call SWGSTP(0.02D0)
call wgscl(idisright,"Size of attractors",0.0D0,0.3D0,attsphsize,2,idisattsize)
if (isys==1) then
	call WGLAB(idisright,"Basins:",ibasinseltext)
	CALL WGBAS(idisright,"FORM",idisbotrig)
	call swgwin(0,5,80,130)
	call swgtyp("VSCROLL","LIST") !This is neccessary, else some latter orbitals cannot be displayed on the list
	CALL WGLIS(idisbotrig,basinlist,1,idisbasinplot)
else if (isys==2.or.isys==3) then
	call WGLAB(idisright2,"Basins:",ibasinseltext)
	call swgtyp("SCROLL","LIST") !This is neccessary, else some latter orbitals cannot be displayed on the list
	CALL WGLIS(idisright2,basinlist,1,idisbasinplot)
end if

!Widget response
call SWGCBK(idisreturn,GUIreturn)
call SWGCBK(idisrotleft,rotleft)
call SWGCBK(idisrotright,rotright)
call SWGCBK(idisrotup,rotup)
call SWGCBK(idisrotdown,rotdown)
call SWGCBK(idiszoomin,zoomin)
call SWGCBK(idiszoomout,zoomout)
call SWGCBK(idisreset,resetview)
call SWGCBK(idissavepic,savepic)
call SWGCBK(idisshowmol,setshowmolstruct)
call SWGCBK(idisshowaxis,ifshowaxis)
call SWGCBK(idisshowatmlab,ifshowatmlabel)
call SWGCBK(idisshowattlab,ifshowattlabel)
call SWGCBK(idisdrawinternalbasin,ifdrawinternalbasin)
call SWGCBK(idisatmsize,setatmsize)
call SWGCBK(idisbondradius,setbondradius)
call SWGCBK(idislabelsize,setlabelsize)
call SWGCBK(idisattsize,setattsize)
call SWGCBK(idisbasinplot,showbasinsel) 
CALL SWGSPC(4.0D0,0.5D0) !Reset the default widget spacing
call swgtyp("HORI","SCALE") !Reset the default mode for list widget
idrawbasinidx=-10 !Don't draw basins by default
CALL WGFIN
end subroutine






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! ----------------- Dislin GUI callback subroutine
subroutine showorbinfo(id)
integer,intent (in) :: id
character*10 orbtype
write(*,*) "Orbital list:"
do i=1,nmo
	if (MOtype(i)==0) orbtype="A+B"
	if (MOtype(i)==1) orbtype="Alpha"
	if (MOtype(i)==2) orbtype="Beta"
	write(*,"(' Orb:',i6,' Ene(a.u./eV):',f13.6,f13.4,' Occ:',f10.6,' Type: ',a)") i,MOene(i),MOene(i)*au2eV,MOocc(i),trim(orbtype)
end do
end subroutine

subroutine showorbinfo2(id) !Only show up to LUMO+10, works for wfntype==0,1,2
integer,intent (in) :: id
character*10 orbtype
if (wfntype==0.or.wfntype==2) then
	write(*,*) "Orbital list:"
	do nmoend=1,nmo
		if (MOocc(nmoend)==0D0) exit
	end do
	nmoend=nmoend+10
	if (nmoend>nmo) nmoend=nmo
	do i=1,nmoend
		if (MOtype(i)==0) orbtype="A+B"
		if (MOtype(i)==1) orbtype="Alpha"
		write(*,"(' Orb:',i6,' Ene(a.u./eV):',f13.6,f13.4,' Occ:',f10.6,' Type: ',a)") i,MOene(i),MOene(i)*au2eV,MOocc(i),trim(orbtype)
	end do
else if (wfntype==1) then
	do iLUMOa=1,nmo
		if (MOocc(iLUMOa)==0) exit
	end do
	do iLUMOb=nmo,1,-1
		if (MOocc(iLUMOb)==1) exit
	end do
	iLUMOb=iLUMOb+1
	do ibeta=1,nmo
		if (MOtype(ibeta)==2) exit
	end do
	iaend=iLUMOa+10
	if (iaend>=ibeta) iaend=ibeta-1
	ibend=iLUMOb+10
	if (ibend>nmo) ibend=nmo
	write(*,*) "Alpha orbital list:"
	orbtype="Alpha"
	do i=1,iaend
		write(*,"(' Orb:',i6,' Ene(a.u./eV):',f13.6,f13.4,' Occ:',f10.6,' Type: ',a)") i,MOene(i),MOene(i)*au2eV,MOocc(i),trim(orbtype)
	end do
	write(*,*) "Beta orbital list:"
	orbtype="Beta"
	do i=ibeta,ibend
		write(*,"(' Orb:',i6,' Ene(a.u./eV):',f13.6,f13.4,' Occ:',f10.6,' Type: ',a)") i,MOene(i),MOene(i)*au2eV,MOocc(i),trim(orbtype)
	end do
end if

end subroutine

subroutine rotleft(id)
integer,intent (in) :: id
XVU=XVU+10
end subroutine

subroutine rotright(id)
integer,intent (in) :: id
XVU=XVU-10
end subroutine

subroutine rotup(id)
integer,intent (in) :: id
if (YVU<90D0) YVU=YVU+10
end subroutine

subroutine rotdown(id)
integer,intent (in) :: id
if (YVU>-90D0) YVU=YVU-10
end subroutine

subroutine zoomin(id)
integer,intent (in) :: id
ZVU=ZVU-1
if (ZVU==2) call SWGATT(idiszoomin,"INACTIVE","STATUS") !Too near the molecule, disable zoom in button
end subroutine

subroutine zoomout(id)
integer,intent (in) :: id
if (ZVU==2) call SWGATT(idiszoomin,"ACTIVE","STATUS")
ZVU=ZVU+1
end subroutine

subroutine resetview(id)
integer,intent (in) :: id
XVU=150.0D0
YVU=30.0D0
if (GUI_mode==1.or.GUI_mode==3) then
	bondcrit=1.15D0
	textheigh=30.0D0
	ratioatmsphere=1.0D0
	bondradius=0.2D0
	ishowatmlab=1
	ishowaxis=1
	call swgbut(idisshowatmlab,ishowatmlab)
	call swgbut(idisshowaxis,ishowaxis)
	if (GUI_mode==1) then
		sur_value=0.05D0
		ZVU=6.0D0
		call swgscl(idisisosurscl,sur_value)
		call swgscl(idisbondradius,bondradius)
		call swgscl(idisatmsize,ratioatmsphere)
		call swgscl(idisbondcrit,bondcrit)
		call swgscl(idislabelsize,textheigh)
	else if (GUI_mode==3) then
		ZVU=7.0D0
		isosurshowboth=1
		ishowdatarange=0
		idrawmol=1
		idrawisosur=1
		call swgbut(idisshowbothsign,isosurshowboth)
		call swgbut(idisshowdatarange,ishowdatarange)
		call swgbut(idisshowmol,idrawmol)
		call swgbut(idisshowisosur,idrawisosur)
	end if
else if (GUI_mode==4) then
	ZVU=5.0D0 !Let the system seems closer
	ishowatmlab=0
	ishowCPlab=0
	ishowpathlab=0
	ishowaxis=1
	ishow3n3=1
	ishow3n1=1
	ishow3p1=1
	ishow3p3=1
	idrawpath=1
	idrawbassurf=1
	bondradius=0.07D0
	ratioatmsphere=0.6D0
	textheigh=36
	call swgbut(idisshowatmlab,ishowatmlab)
	call swgbut(idisshowCPlab,ishowCPlab)
	call swgbut(idisshowpathlab,ishowpathlab)
	call swgbut(idisshowaxis,ishowaxis)
	call swgbut(idisshowmol,idrawmol)
	call swgbut(idisshow3n3,ishow3n3)
	call swgbut(idisshow3n1,ishow3n1)
	call swgbut(idisshow3p1,ishow3p1)
	call swgbut(idisshow3p3,ishow3p3)
	call swgbut(idisshowpath,idrawpath)
	call swgbut(idisshowbassurf,idrawbassurf)
	call swgscl(idisbondradius,bondradius)
	call swgscl(idisatmsize,ratioatmsphere)
	call swgscl(idislabelsize,textheigh)
else if (GUI_mode==5) then
	ZVU=6.0D0
	textheigh=30.0D0
	ratioatmsphere=1.0D0
	bondradius=0.2D0
	ishowatmlab=1
	ishowaxis=1
	ishowlocminlab=0
	ishowlocmaxlab=0
	ishowlocminpos=1
	ishowlocmaxpos=1
	call swgscl(idislabelsize,textheigh)
	call swgscl(idisatmsize,ratioatmsphere)
	call swgscl(idisbondradius,bondradius)
	call swgbut(idisshowatmlab,ishowatmlab)
	call swgbut(idisshowaxis,ishowaxis)
	call swgbut(idisshowlocminlab,ishowlocminlab)
	call swgbut(idisshowlocmaxlab,ishowlocmaxlab)
	call swgbut(idisshowlocminpos,ishowlocminpos)
	call swgbut(idisshowlocmaxpos,ishowlocmaxpos)
else if (GUI_mode==6) then
	ZVU=6.0D0
	idrawmol=1
	ishowaxis=1
	ishowatmlab=0
	ishowattlab=1
	idrawinternalbasin=0
	ratioatmsphere=1.0D0
	bondradius=0.2D0
	textheigh=40.0D0
	attsphsize=0.1D0
	call swgbut(idisshowmol,idrawmol)
	call swgbut(idisshowaxis,ishowaxis)
	call swgbut(idisshowatmlab,ishowatmlab)
	call swgbut(idisshowattlab,ishowattlab)
	call swgbut(idisdrawinternalbasin,idrawinternalbasin)
	call swgscl(idisatmsize,ratioatmsphere)
	call swgscl(idisbondradius,bondradius)
	call swgscl(idislabelsize,textheigh)
	call swgscl(idisattsize,attsphsize)
else if (GUI_mode==2) then
	ZVU=7.0D0
end if
call SWGATT(idiszoomin,"ACTIVE","STATUS")
end subroutine

subroutine savepic(id)
integer,intent (in) :: id
!There is a bug, present at lease in DISLIN 10.2, that is in Windows environment, when saving picture and transparent style is used, the SURISO often crashes
if (isys==1.and.(GUI_mode==1.or.GUI_mode==3).and.isosur1style==5) then !When showing orbital isosurface or cubmat isosurface
	call DWGMSG("Error: For Windows verison, picture cannot be saved when transparent face style is used")
else
	isavepic=1
	call DWGMSG("The graph has been saved to a file with ""DISLIN"" prefix in current folder")
	isavepic=0
end if
end subroutine

subroutine ifshowatmlabel(id)
integer,intent (in) :: id
call GWGBUT(id,istat)
if (istat==0) ishowatmlab=0
if (istat==1) ishowatmlab=1
end subroutine

subroutine ifshowattlabel(id)
integer,intent (in) :: id
call GWGBUT(id,istat)
if (istat==0) ishowattlab=0
if (istat==1) ishowattlab=1
end subroutine

subroutine ifdrawinternalbasin(id)
integer,intent (in) :: id
call GWGBUT(id,istat)
if (istat==0) idrawinternalbasin=0
if (istat==1) idrawinternalbasin=1
end subroutine

subroutine ifshowCPlabel(id)
integer,intent (in) :: id
call GWGBUT(id,istat)
if (istat==0) ishowCPlab=0
if (istat==1) ishowCPlab=1
end subroutine

subroutine ifshowpathlabel(id)
integer,intent (in) :: id
call GWGBUT(id,istat)
if (istat==0) ishowpathlab=0
if (istat==1) ishowpathlab=1
end subroutine

subroutine ifshowaxis(id)
integer,intent (in) :: id
call GWGBUT(id,istat)
if (istat==0) ishowaxis=0
if (istat==1) ishowaxis=1
end subroutine

subroutine ifshowisosur(id)
integer,intent (in) :: id
call GWGBUT(id,istat)
if (istat==0) idrawisosur=0
if (istat==1) idrawisosur=1
end subroutine

!Warpper of showorbsel
subroutine showorbsellist(id)
integer,intent (in) :: id
call showorbsel(id,1)
end subroutine
subroutine showorbselbox(id)
integer,intent (in) :: id
call showorbsel(id,2)
end subroutine
!! Select an orbital and calculate cube data for it
!if itype==1, load index from list; if =2, load index from text box
subroutine showorbsel(id,itype)
use function
integer,intent (in) :: id,itype
real*8 molxlen,molylen,molzlen
character*10 corbtype
! Set grid for calculating cube data
molxlen=(maxval(a%x)-minval(a%x))+2*aug3D
molylen=(maxval(a%y)-minval(a%y))+2*aug3D
molzlen=(maxval(a%z)-minval(a%z))+2*aug3D
orgx=minval(a%x)-aug3D
orgy=minval(a%y)-aug3D
orgz=minval(a%z)-aug3D
if (isys==1) then
	dx=(molxlen*molylen*molzlen/dfloat(nprevorbgrid))**(1.0D0/3.0D0)
	dy=dx
	dz=dx
	nx=nint(molxlen/dx)+1
	ny=nint(molylen/dy)+1
	nz=nint(molzlen/dz)+1
else if (isys==2.or.isys==3) then !There is bug in Linux version of suriso, only when nx=ny=nz can ensure the routine doesn't crash
	nx=nint(nprevorbgrid**(1.0D0/3.0D0))
	ny=nx
	nz=nx
	dx=molxlen/(nx-1)
	dy=molylen/(ny-1)
	dz=molzlen/(nz-1)
end if
if (itype==1) then
	call GWGLIS(id,isel)
else if (itype==2) then
	call GWGINT(id,isel) !Input 0 means "none"
	if (isel<0.or.isel>nmo) isel=0
	isel=isel+1
end if
if (isel==1) then !Namely "None" in the orbital list
	if (isosursec==0) then
		idrawisosur=0
		if (allocated(cubmat)) deallocate(cubmat)
		call swgatt(idisisosursec,"INACTIVE","STATUS")
	else if (isosursec==1) then !When the second isosurface is selected, select "NONE" only clean cubmattmp
		if (allocated(cubmattmp)) deallocate(cubmattmp)
	end if
else
	idrawisosur=1
	call swgatt(idisisosursec,"ACTIVE","STATUS")
	corbtype=" "
	if (motype(isel-1)==0) corbtype="A+B"
	if (motype(isel-1)==1) corbtype="Alpha"
	if (motype(isel-1)==2) corbtype="Beta"
	write(*,"(' Orb:',i6,' Ene(a.u./eV):',f13.6,f13.4,' Occ:',f10.6,' Type: ',a)") isel-1,MOene(isel-1),MOene(isel-1)*27.2113838,MOocc(isel-1),trim(corbtype)
	if (isosursec==0) then !Save cube data for isosurface 1 to cubmat
		if (allocated(cubmat)) deallocate(cubmat)
		allocate(cubmat(nx,ny,nz))
		call swgatt(iprogbar,"ACTIVE","STATUS")
		ifinish=0
		!$OMP parallel do PRIVATE(i,j,k) SHARED(cubmat,ifinish) NUM_THREADS(rtNThreads())
		do k=1,nz
			do j=1,ny
				do i=1,nx
					cubmat(i,j,k)=fmo(orgx+(i-1)*dx,orgy+(j-1)*dy,orgz+(k-1)*dz,isel-1)
				end do
			end do
			ifinish=ifinish+1
			!$OMP CRITICAL
			call swgval(iprogbar,ifinish/dfloat(nz))
			!$OMP end CRITICAL
		end do
		!$OMP end parallel do
		if (ifixorbsign==1.and.sum(cubmat)<0) cubmat=-cubmat
	else if (isosursec==1) then !Save cube data for isosurface 2 to cubmattmp
		if (allocated(cubmattmp)) deallocate(cubmattmp)
		allocate(cubmattmp(nx,ny,nz))
		call swgatt(iprogbar,"ACTIVE","STATUS")
		ifinish=0
		!$OMP parallel do PRIVATE(i,j,k) SHARED(cubmat,ifinish) NUM_THREADS(rtNThreads())
		do k=1,nz
			do j=1,ny
				do i=1,nx
					cubmattmp(i,j,k)=fmo(orgx+(i-1)*dx,orgy+(j-1)*dy,orgz+(k-1)*dz,isel-1)
				end do
			end do
			ifinish=ifinish+1
			!$OMP CRITICAL
			call swgval(iprogbar,ifinish/dfloat(nz))
			!$OMP end CRITICAL
		end do
		!$OMP end parallel do
		if (ifixorbsign==1.and.sum(cubmattmp)<0) cubmattmp=-cubmattmp
	end if
end if
end subroutine

subroutine showbasinsel(id)
use basinintmod
integer,intent (in) :: id
call GWGLIS(id,isel)
idrawbasinidx=isel-1
if (isel==1) idrawbasinidx=-10 !Don't draw basins
if (isel==numrealatt+2) idrawbasinidx=0 !Unassigned
if (isel==numrealatt+3) idrawbasinidx=-1 !Moved to boundary
end subroutine

subroutine setbondcrit(id)
integer,intent (in) :: id
call GWGSCL(id,bondcrit)
end subroutine

subroutine GUIreturn(id)
integer,intent (in) :: id
call sendok
end subroutine

subroutine setlabelsize(id)
integer,intent (in) :: id
call GWGSCL(id,textheigh)
end subroutine

subroutine setatmsize(id)
integer,intent (in) :: id
call GWGSCL(id,ratioatmsphere)
end subroutine

subroutine setattsize(id)
integer,intent (in) :: id
call GWGSCL(id,attsphsize)
end subroutine

subroutine setbondradius(id)
integer,intent (in) :: id
call GWGSCL(id,bondradius)
end subroutine

subroutine setisosurscl(id) !Drag scale bar, change sur_value & text
integer,intent (in) :: id
character*20 temp
call GWGSCL(id,sur_value)
if (GUI_mode==3) then
	write(temp,"(f8.3)") sur_value
	call SWGTXT(idisscrval,temp)
end if
end subroutine

integer,intent (in) :: id
call GWGFLT(id,sur_value)
if (sur_value<5D0.and.sur_value>-5D0) call SWGSCL(idisisosurscl,sur_value)
end subroutine

subroutine setshowbothsign(id)
integer,intent (in) :: id
call GWGBUT(id,istat)
if (istat==0) isosurshowboth=0
if (istat==1) isosurshowboth=1
end subroutine

subroutine setshowmolstruct(id)
integer,intent (in) :: id
call GWGBUT(id,istat)
if (istat==0) idrawmol=0
if (istat==1) idrawmol=1
end subroutine

subroutine setshowpath(id)
integer,intent (in) :: id
call GWGBUT(id,istat)
if (istat==0) idrawpath=0
if (istat==1) idrawpath=1
end subroutine

subroutine setshowbassurf(id)
integer,intent (in) :: id
call GWGBUT(id,istat)
if (istat==0) idrawbassurf=0
if (istat==1) idrawbassurf=1
end subroutine

subroutine setshowatmlab(id)
integer,intent (in) :: id
call GWGBUT(id,istat)
if (istat==0) ishowatmlab=0
if (istat==1) ishowatmlab=1
end subroutine

subroutine setshowdatarange(id)
integer,intent (in) :: id
call GWGBUT(id,istat)
if (istat==0) ishowdatarange=0
if (istat==1) ishowdatarange=1
end subroutine

subroutine ifshow3n3(id)
integer,intent (in) :: id
call GWGBUT(id,istat)
if (istat==0) ishow3n3=0
if (istat==1) ishow3n3=1
end subroutine
subroutine ifshow3n1(id)
integer,intent (in) :: id
call GWGBUT(id,istat)
if (istat==0) ishow3n1=0
if (istat==1) ishow3n1=1
end subroutine
subroutine ifshow3p1(id)
integer,intent (in) :: id
call GWGBUT(id,istat)
if (istat==0) ishow3p1=0
if (istat==1) ishow3p1=1
end subroutine
subroutine ifshow3p3(id)
integer,intent (in) :: id
call GWGBUT(id,istat)
if (istat==0) ishow3p3=0
if (istat==1) ishow3p3=1
end subroutine

!For molecular surface analysis
subroutine setshowlocminlab(id)
integer,intent (in) :: id
call GWGBUT(id,istat)
if (istat==0) ishowlocminlab=0
if (istat==1) ishowlocminlab=1
end subroutine
subroutine setshowlocmaxlab(id)
integer,intent (in) :: id
call GWGBUT(id,istat)
if (istat==0) ishowlocmaxlab=0
if (istat==1) ishowlocmaxlab=1
end subroutine
subroutine setshowlocminpos(id)
integer,intent (in) :: id
call GWGBUT(id,istat)
if (istat==0) ishowlocminpos=0
if (istat==1) ishowlocminpos=1
end subroutine
subroutine setshowlocmaxpos(id)
integer,intent (in) :: id
call GWGBUT(id,istat)
if (istat==0) ishowlocmaxpos=0
if (istat==1) ishowlocmaxpos=1
end subroutine



!Set number of grid points for viewing orbitals (namely control isosurface quality)
subroutine setisosurnumpt(id)
integer,intent (in) :: id
character inpstring*30
nprevorbgridold=nprevorbgrid
CALL SWGWTH(40)
write(inpstring,"(i15)") nprevorbgrid
call dwgtxt("Input the number of grid points|Higher number leads to finer quality",inpstring)
read(inpstring,*) nprevorbgrid
if (nprevorbgrid/=nprevorbgridold) then !Remove current isosurface and corresponding grid data, recover initial state when entering the GUI
	if (allocated(cubmat)) deallocate(cubmat)
	if (allocated(cubmattmp)) deallocate(cubmattmp)
	idrawisosur=0
	isosursec=0
	call swgatt(idisisosursec,"INACTIVE","STATUS")
end if
CALL SWGWTH(20) !Recover default
end subroutine

!If show cubmattmp
subroutine ifisosursec(id)
integer,intent (in) :: id
call GWGBUT(id,istat)
if (istat==0) isosursec=0
if (istat==1) isosursec=1
end subroutine

!Set isosurface style. We absolutely avoid that only one isosurface is transparent. So, when we set one of isosurface to transparent, then another too.
!When we set one of isosurface to a specific style, then we check if the old style of another isosurface is transparent, if yes, we set its style to current style together
subroutine setisosur1solid(id)
integer,intent (in) :: id
isosur1style=1
if (isosur2style==5) isosur2style=1
end subroutine
subroutine setisosur1line(id)
integer,intent (in) :: id
isosur1style=2
if (isosur2style==5) isosur2style=2
end subroutine
subroutine setisosur1point(id)
integer,intent (in) :: id
isosur1style=3
if (isosur2style==5) isosur2style=3
end subroutine
subroutine setisosur1solidmesh(id)
integer,intent (in) :: id
isosur1style=4
if (isosur2style==5) isosur2style=4
end subroutine
subroutine setisosur1tpr(id)
integer,intent (in) :: id
isosur1style=5
isosur2style=5
end subroutine
!----
subroutine setisosur2solid(id)
integer,intent (in) :: id
isosur2style=1
if (isosur1style==5) isosur1style=1
end subroutine
subroutine setisosur2line(id)
integer,intent (in) :: id
isosur2style=2
if (isosur1style==5) isosur1style=2
end subroutine
subroutine setisosur2point(id)
integer,intent (in) :: id
isosur2style=3
if (isosur1style==5) isosur1style=3
end subroutine
subroutine setisosur2solidmesh(id)
integer,intent (in) :: id
isosur2style=4
if (isosur1style==5) isosur1style=4
end subroutine
subroutine setisosur2tpr(id)
integer,intent (in) :: id
isosur1style=5
isosur2style=5
end subroutine
!----
subroutine setisosurallsolid(id)
integer,intent (in) :: id
isosur1style=1
isosur2style=1
end subroutine
subroutine setisosurallline(id)
integer,intent (in) :: id
isosur1style=2
isosur2style=2
end subroutine
subroutine setisosurallpoint(id)
integer,intent (in) :: id
isosur1style=3
isosur2style=3
end subroutine
subroutine setisosurallsolidmesh(id)
integer,intent (in) :: id
isosur1style=4
isosur2style=4
end subroutine
subroutine setisosuralltpr(id)
integer,intent (in) :: id
isosur1style=5
isosur2style=5
end subroutine

!Set color for solid representation of isosurface 1
subroutine setisosur1solidclr(id)
integer,intent (in) :: id
character inpstring*30
CALL SWGWTH(40)
write(inpstring,"(f4.2,',',f4.2,',',f4.2)") clrRcub1same,clrGcub1same,clrBcub1same
call dwgtxt("Input R,G,B value, for the same sign part",inpstring)
read(inpstring,*) clrRcub1same,clrGcub1same,clrBcub1same
write(inpstring,"(f4.2,',',f4.2,',',f4.2)") clrRcub1oppo,clrGcub1oppo,clrBcub1oppo
call dwgtxt("Input R,G,B value, for the opposite sign part",inpstring)
read(inpstring,*) clrRcub1oppo,clrGcub1oppo,clrBcub1oppo
CALL SWGWTH(20) !Recover default
end subroutine

!Set color for mesh and points representation of isosurface 1
subroutine setisosur1meshptclr(id)
integer,intent (in) :: id
character inpstring*30
CALL SWGWTH(40)
write(inpstring,"(f4.2,',',f4.2,',',f4.2)") clrRcub1samemeshpt,clrGcub1samemeshpt,clrBcub1samemeshpt
call dwgtxt("Input R,G,B value, for the same sign part",inpstring)
read(inpstring,*) clrRcub1samemeshpt,clrGcub1samemeshpt,clrBcub1samemeshpt
write(inpstring,"(f4.2,',',f4.2,',',f4.2)") clrRcub1oppomeshpt,clrGcub1oppomeshpt,clrBcub1oppomeshpt
call dwgtxt("Input R,G,B value, for the opposite sign part",inpstring)
read(inpstring,*) clrRcub1oppomeshpt,clrGcub1oppomeshpt,clrBcub1oppomeshpt
CALL SWGWTH(20) !Recover default
end subroutine

!Set opacity for transparent face representation of isosurface 1
subroutine setisosur1opa(id)
integer,intent (in) :: id
character inpstring*30
CALL SWGWTH(40)
write(inpstring,"(f4.2)") opacitycub1
call dwgtxt("Input opacity, between 0.0 and 1.0",inpstring)
read(inpstring,*) opacitycub1
CALL SWGWTH(20) !Recover default
end subroutine

!Set color for solid representation of isosurface 2
subroutine setisosur2solidclr(id)
integer,intent (in) :: id
character inpstring*30
CALL SWGWTH(40)
write(inpstring,"(f4.2,',',f4.2,',',f4.2)") clrRcub2same,clrGcub2same,clrBcub2same
call dwgtxt("Input R,G,B value, for the same sign part",inpstring)
read(inpstring,*) clrRcub2same,clrGcub2same,clrBcub2same
write(inpstring,"(f4.2,',',f4.2,',',f4.2)") clrRcub2oppo,clrGcub2oppo,clrBcub2oppo
call dwgtxt("Input R,G,B value, for the opposite sign part",inpstring)
read(inpstring,*) clrRcub2oppo,clrGcub2oppo,clrBcub2oppo
CALL SWGWTH(20) !Recover default
end subroutine

!Set color for mesh and points representation of isosurface 2
subroutine setisosur2meshptclr(id)
integer,intent (in) :: id
character inpstring*30
CALL SWGWTH(40)
write(inpstring,"(f4.2,',',f4.2,',',f4.2)") clrRcub2samemeshpt,clrGcub2samemeshpt,clrBcub2samemeshpt
call dwgtxt("Input R,G,B value, for the same sign part",inpstring)
read(inpstring,*) clrRcub2samemeshpt,clrGcub2samemeshpt,clrBcub2samemeshpt
write(inpstring,"(f4.2,',',f4.2,',',f4.2)") clrRcub2oppomeshpt,clrGcub2oppomeshpt,clrBcub2oppomeshpt
call dwgtxt("Input R,G,B value, for the opposite sign part",inpstring)
read(inpstring,*) clrRcub2oppomeshpt,clrGcub2oppomeshpt,clrBcub2oppomeshpt
CALL SWGWTH(20) !Recover default
end subroutine

!Set opacity for transparent face representation of isosurface 2
subroutine setisosur2opa(id)
integer,intent (in) :: id
character inpstring*30
CALL SWGWTH(40)
write(inpstring,"(f4.2)") opacitycub2
call dwgtxt("Input opacity, between 0.0 and 1.0",inpstring)
read(inpstring,*) opacitycub2
CALL SWGWTH(20) !Recover default
end subroutine

!Set lighting 1~8
subroutine setlight1(id)
integer,intent (in) :: id
if (ienablelight1==1) then
	ienablelight1=0
else
	ienablelight1=1
end if
call showlightsetting
end subroutine
subroutine setlight2(id)
integer,intent (in) :: id
if (ienablelight2==1) then
	ienablelight2=0
else
	ienablelight2=1
end if
call showlightsetting
end subroutine
subroutine setlight3(id)
integer,intent (in) :: id
if (ienablelight3==1) then
	ienablelight3=0
else
	ienablelight3=1
end if
call showlightsetting
end subroutine
subroutine setlight4(id)
integer,intent (in) :: id
if (ienablelight4==1) then
	ienablelight4=0
else
	ienablelight4=1
end if
call showlightsetting
end subroutine
subroutine setlight5(id)
integer,intent (in) :: id
if (ienablelight5==1) then
	ienablelight5=0
else
	ienablelight5=1
end if
call showlightsetting
end subroutine
subroutine setlight6(id)
integer,intent (in) :: id
if (ienablelight6==1) then
	ienablelight6=0
else
	ienablelight6=1
end if
call showlightsetting
end subroutine
subroutine setlight7(id)
integer,intent (in) :: id
if (ienablelight7==1) then
	ienablelight7=0
else
	ienablelight7=1
end if
call showlightsetting
end subroutine
subroutine setlight8(id)
integer,intent (in) :: id
if (ienablelight8==1) then
	ienablelight8=0
else
	ienablelight8=1
end if
call showlightsetting
end subroutine
subroutine setlightall0(id)
integer,intent (in) :: id
ienablelight1=0
ienablelight2=0
ienablelight3=0
ienablelight4=0
ienablelight5=0
ienablelight6=0
ienablelight7=0
ienablelight8=0
end subroutine
subroutine setlightall1(id)
integer,intent (in) :: id
ienablelight1=1
ienablelight2=1
ienablelight3=1
ienablelight4=1
ienablelight5=1
ienablelight6=1
ienablelight7=1
ienablelight8=1
end subroutine
!Show light setting
subroutine showlightsetting
write(*,*) "Current lighting status:"
if (ienablelight1==1) then
	write(*,*) "Light 1: Enabled"
else
	write(*,*) "Light 1: Disabled"
end if
if (ienablelight2==1) then
	write(*,*) "Light 2: Enabled"
else
	write(*,*) "Light 2: Disabled"
end if
if (ienablelight3==1) then
	write(*,*) "Light 3: Enabled"
else
	write(*,*) "Light 3: Disabled"
end if
if (ienablelight4==1) then
	write(*,*) "Light 4: Enabled"
else
	write(*,*) "Light 4: Disabled"
end if
if (ienablelight5==1) then
	write(*,*) "Light 5: Enabled"
else
	write(*,*) "Light 5: Disabled"
end if
if (ienablelight6==1) then
	write(*,*) "Light 6: Enabled"
else
	write(*,*) "Light 6: Disabled"
end if
if (ienablelight7==1) then
	write(*,*) "Light 7: Enabled"
else
	write(*,*) "Light 7: Disabled"
end if
if (ienablelight8==1) then
	write(*,*) "Light 8: Enabled"
else
	write(*,*) "Light 8: Disabled"
end if
end subroutine


end module
