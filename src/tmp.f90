
!!------- Generate atomic radial density files at different states, used for such as Hirshfeld-I
!"atmrad" in current folder is used as working directory
!-2,-1,0,+1,+2 charge states of each element will be calculated to produce atomic .wfn file by Gaussian, predefined ground state multiplicity is used
!After that, radial density file (.rad) is generated for each state of each element
!If atomic wfn file is already existed, calculation will be skipped
!Radial distance values are the same as built-in atomic density, i.e. those in atmraddens.f90
subroutine atmradfile
use defvar
use function
implicit real*8 (a-h,o-z)
character c80tmp*80
integer :: chgmulti(nelesupp,-3:3)=0 !Ground state multiplicity of each charge state of each element. If value=0, means undefined

!Define chgmulti for elements for possible states
!H,Li,Na,K,Rb,Cs
chgmulti(1,0)=2
chgmulti(1,1)=1
chgmulti(1,-1)=1
chgmulti(3,:)=chgmulti(1,:)
chgmulti(11,:)=chgmulti(1,:)
chgmulti(19,:)=chgmulti(1,:)
chgmulti(37,:)=chgmulti(1,:)
chgmulti(55,:)=chgmulti(1,:)
!He,Ne,Ar,Kr,Xe,Rn
chgmulti(2,0)=1
chgmulti(2,1)=2
chgmulti(2,-1)=2
chgmulti(10,:)=chgmulti(2,:)
chgmulti(18,:)=chgmulti(2,:)
chgmulti(36,:)=chgmulti(2,:)
chgmulti(54,:)=chgmulti(2,:)
chgmulti(86,:)=chgmulti(2,:)
!Be,Mg,Ca,Sr,Ba
chgmulti(4,0)=1
chgmulti(4,1)=2
chgmulti(4,2)=1
chgmulti(4,-1)=2
chgmulti(12,:)=chgmulti(4,:)
chgmulti(20,:)=chgmulti(4,:)
chgmulti(38,:)=chgmulti(4,:)
chgmulti(56,:)=chgmulti(4,:)
!B,Al,Ga,In,Tl
chgmulti(5,0)=2
chgmulti(5,1)=1
chgmulti(5,2)=2
chgmulti(5,-1)=3
chgmulti(5,-2)=4
chgmulti(13,:)=chgmulti(5,:)
chgmulti(31,:)=chgmulti(5,:)
chgmulti(49,:)=chgmulti(5,:)
chgmulti(81,:)=chgmulti(5,:)
!C,Si,Ge,Sn,Pb
chgmulti(6,0)=3
chgmulti(6,1)=2
chgmulti(6,2)=1
chgmulti(6,-1)=4
chgmulti(6,-2)=3
chgmulti(14,:)=chgmulti(6,:)
chgmulti(32,:)=chgmulti(6,:)
chgmulti(50,:)=chgmulti(6,:)
chgmulti(82,:)=chgmulti(6,:)
!N,P,As,Sb,Bi
chgmulti(7,0)=4
chgmulti(7,1)=3
chgmulti(7,2)=2
chgmulti(7,-1)=3
chgmulti(7,-2)=2
chgmulti(15,:)=chgmulti(7,:)
chgmulti(33,:)=chgmulti(7,:)
chgmulti(51,:)=chgmulti(7,:)
chgmulti(83,:)=chgmulti(7,:)
!O,S,Se,Te,Po
chgmulti(8,0)=3
chgmulti(8,1)=4
chgmulti(8,2)=3
chgmulti(8,-1)=2
chgmulti(8,-2)=1
chgmulti(16,:)=chgmulti(8,:)
chgmulti(34,:)=chgmulti(8,:)
chgmulti(52,:)=chgmulti(8,:)
chgmulti(84,:)=chgmulti(8,:)
!F,Cl,Br,I,At
chgmulti(9,0)=2
chgmulti(9,1)=3
chgmulti(9,2)=4
chgmulti(9,-1)=1
chgmulti(17,:)=chgmulti(9,:)
chgmulti(35,:)=chgmulti(9,:)
chgmulti(53,:)=chgmulti(9,:)
chgmulti(85,:)=chgmulti(9,:)
!Spin multiplicity of transition metal for each state is determined by chemical intuition as well as a few single point energy data
!For simplicity, I assume that later elements in each row has identical configuration, of course this is incorrect but not too bad
!Sc (3d1,4s2)
chgmulti(21,0)=2
chgmulti(21,1)=3
chgmulti(21,2)=2
chgmulti(21,-1)=3
chgmulti(39,:)=chgmulti(21,:) !Y
chgmulti(57,:)=chgmulti(21,:) !La
!Ti (3d2,4s2)
chgmulti(22,0)=3
chgmulti(22,1)=4
chgmulti(22,2)=3
chgmulti(22,-1)=4
chgmulti(40,:)=chgmulti(22,:) !Zr
chgmulti(72,:)=chgmulti(22,:) !Hf
!V  (3d3,4s2)
chgmulti(23,0)=4
chgmulti(23,1)=5
chgmulti(23,2)=4
chgmulti(23,-1)=5
chgmulti(41,:)=chgmulti(23,:) !Nb
chgmulti(73,:)=chgmulti(23,:) !Ta
!Cr (3d5,4s1)
chgmulti(24,0)=7
chgmulti(24,1)=6
chgmulti(24,2)=5
chgmulti(24,-1)=6
chgmulti(42,:)=chgmulti(24,:) !Mo
chgmulti(74,:)=chgmulti(24,:) !W
!Mn (3d5,4s2)
chgmulti(25,0)=6
chgmulti(25,1)=7
chgmulti(25,2)=6
chgmulti(25,-1)=5
chgmulti(43,:)=chgmulti(25,:) !Tc
chgmulti(75,:)=chgmulti(25,:) !Re
!Fe (3d6,4s2)
chgmulti(26,0)=5
chgmulti(26,1)=6
chgmulti(26,2)=7
chgmulti(26,-1)=4
chgmulti(44,:)=chgmulti(26,:) !Ru
chgmulti(76,:)=chgmulti(26,:) !Os
!Co (3d7,4s2)
chgmulti(27,0)=4
chgmulti(27,1)=5
chgmulti(27,2)=4
chgmulti(27,-1)=3
chgmulti(45,:)=chgmulti(27,:) !Rh
chgmulti(77,:)=chgmulti(27,:) !Ir
!Ni (3d8,4s2)
chgmulti(28,0)=3
chgmulti(28,1)=4
chgmulti(28,2)=3
chgmulti(28,-1)=2
chgmulti(46,:)=chgmulti(28,:) !Pd
chgmulti(78,:)=chgmulti(28,:) !Pt
!Cu (3d10,4s1)
chgmulti(29,0)=2
chgmulti(29,1)=1
chgmulti(29,2)=2
chgmulti(29,-1)=1
chgmulti(47,:)=chgmulti(29,:) !Ag
chgmulti(79,:)=chgmulti(29,:) !Au
!Zn (3d10,4s2)
chgmulti(30,0)=1
chgmulti(30,1)=2
chgmulti(30,2)=1
chgmulti(30,-1)=2
chgmulti(48,:)=chgmulti(30,:) !Cd
chgmulti(80,:)=chgmulti(30,:) !Hg

end subroutine
