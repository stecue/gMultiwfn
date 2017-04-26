module deftype
implicit none
type atomtype
character*2 name !name of atom
integer index !The index in periodic table, if ECP was used, charge will smaller than this value
real*8 x,y,z,charge !Coordinate(Bohr) and charge of atoms
end type

type primtype
integer center,functype !The number of nuclei that the basis function centered on and its function type
real*8 exp !Exponent
end type

type content !Type for grid data points
real*8 x,y,z,value
end type

end module

!============ Store globally shared information
module defvar
use deftype
implicit none
integer :: ido
real*8,parameter :: pi=3.141592653589793D0,b2a=0.529177249D0 !1 Bohr = 0.529177249 Angstrom
real*8,parameter :: au2kcal=627.51D0,au2KJ=2625.5D0,au2eV=27.2113838D0
real*8,parameter :: masse=9.10938215D-31,chge=1.602176487D0,lightc=2.99792458D8,au2debye=2.5417462D0 !masse/chge: Mass/charge of an electron
real*8,parameter :: planckc=6.62606896D-34,h_bar=1.054571628D-34,amu2kg=1.66053878D-27
real*8,parameter :: boltzc=1.3806488D-23,boltzcau=3.1668114D-6,boltzceV=8.6173324D-5 !in J/K, in Hartree/K and in eV/K, respectively
real*8,parameter :: avogacst=6.02214179D23
integer,parameter :: nelesupp=150 !The number of elements supported, ghost(index=0) is not taken into account
real*8 ctrval(1000) !Value of contour lines

!Store important calculated data
real*8,allocatable :: curvex(:),curvey(:),curveytmp(:) !For line plot
real*8,allocatable :: planemat(:,:),planemattmp(:,:) !planemattmp is mainly used to draw contour line of a function on contour map of another function (e.g. vdw surface on ESP contour map)
real*8,allocatable :: cubmat(:,:,:) !cubmat, store density/laplacian...3D-matrix
real*8,allocatable :: cubmattmp(:,:,:) !For cube data exchanging, position of points must be identical to cubmat, so don't use type(content) for lowering memory consumption
real*8,allocatable :: cubmatvec(:,:,:,:) !Used to store vector field
real*8,allocatable :: gradd1(:,:),gradd2(:,:) !Gradient in direction1/2 for gradient line plot
real*8,allocatable :: distmat(:,:) !Distance matrix, in Bohr
character*200 filename,firstfilename
character*80 firstfileonlyname,extctrsetting,cmdarg2 !cmdarg is the parameter of booting multiwfn
character,allocatable :: custommapname(:)*80,customop(:) !Custom operation for custom map/cube file
logical alive
integer,allocatable :: fragatm(:),fragatmbackup(:) !Store the index of atoms in fragment, has no relationship with frag1/frag2. fragatmbackup is used to backup fragatm during custom operation
integer,allocatable :: frag1(:),frag2(:) !These two fragments are only used for bond order analysis/composition analysis etc., store index of basis functions or atoms. Their size just fit their content
integer :: ncustommap=0,imodwfn=0 !if 1, means original occupation number or orbital type have been modified. Modified GTF information is not taken into account
integer :: iorbsel=1 !Which orbital is selected, and its value will be calculated by fmo and calchessmat_mo
integer :: iorbsel2=0 !Which orbital will be plotted together with iorbsel in plane map

character*10 :: colorname(14)=[character(len=10)::"Red","Green","Blue","White","Black","Gray","Cyan","Yellow","Orange","Magenta","Crimson","Dark green","Purple","Brown"] !Color name of useclrind routine
!The name for superheavy atoms are consistent with Stuttgart PP website: http://www.tc.uni-koeln.de/PP/clickpse.en.html
character*2 :: ind2name(0:nelesupp)=(/ "Bq","H ","He", &   !X(number O) is ghost atom
"Li","Be","B ","C ","N ","O ","F ","Ne", & !3~10
"Na","Mg","Al","Si","P ","S ","Cl","Ar", & !11~18
"K ","Ca","Sc","Ti","V ","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr", & !19~36
"Rb","Sr","Y ","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I ","Xe", & !37~54
"Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu", & !55~71
"Hf","Ta","W ","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn", & !72~86
"Fr","Ra","Ac","Th","Pa","U ","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr", & !87~103
"Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","Ut","Fl","Up","Lv","Us","Uo","Un","Ux",("??",ido=121,nelesupp) /) !104~all  Such as Uuo/Uup is replaced by Uo/Up
character*2 :: ind2name_up(0:nelesupp)=(/ "BQ","H ","HE", & !Same as ind2name, but all characters are upper case, to cater to .pdb file
"LI","BE","B ","C ","N ","O ","F ","NE", & !3~10
"NA","MG","AL","SI","P ","S ","CL","AR", & !11~18
"K ","CA","SC","TI","V ","CR","MN","FE","CO","NI","CU","ZN","GA","GE","AS","SE","BR","KR", & !19~36
"RB","SR","Y ","ZR","NB","MO","TC","RU","RH","PD","AG","CD","IN","SN","SB","TE","I ","XE", & !37~54
"CS","BA","LA","CE","PR","ND","PM","SM","EU","GD","TB","DY","HO","ER","TM","YB","LU", & !55~71
"HF","TA","W ","RE","OS","IR","PT","AU","HG","TL","PB","BI","PO","AT","RN", & !72~86
"FR","RA","AC","TH","PA","U ","NP","PU","AM","CM","BK","CF","ES","FM","MD","NO","LR", & !87~103
"RF","DB","SG","BH","HS","MT","DS","RG","CN","UT","FL","UP","LV","US","UO","UN","UX",("??",ido=121,nelesupp) /) !104~all
!Bondi vdW radius, from J.Phys.Chem.,1964,68(3),441-451, unit is Angstrom, will be convert to Bohr when multiwfn start
real*8 :: vdwr(0:nelesupp)=(/ 0.4D0,1.2D0,1.4D0,& !Ghost,H,He
1.82D0,1.77D0,1.74D0,1.7D0,1.55D0,1.52D0,1.47D0,1.54D0,& !Li~Ne
2.27D0,1.73D0,1.73D0,2.1D0,1.8D0,1.8D0,1.75D0,1.88D0,& !Na~Ar
(2.0D0,ido=19,27),1.63D0,1.4D0,1.39D0,1.87D0,2.0D0,1.85D0,1.9D0,1.85D0,2.02D0,& ! Ni~Kr(28~36)
(2.0D0,ido=37,45),1.63D0,1.72D0,1.58D0,1.93D0,2.17D0,2.0D0,2.06D0,1.98D0,2.16D0,& !Pd~Xe(46~54)
(2.0D0,ido=55,77),1.72D0,1.66D0,1.55D0,1.96D0,2.02D0,(2.0D0,ido=83,nelesupp) /) !Pt~Pb(78~82)
!##No use currently!## Modified Bondi vdW radii, but for all main group (except for H,He), use IVA radius in corresponding row. Specifically used to molecular surface decomposition
real*8 :: vdwr_tianlu(0:nelesupp)=(/ 0.4D0,1.7D0,1.7D0,& !Ghost,H,Ne   H and Ne are augmented to carbon radius
(1.7D0,ido=3,10),& !Li~Ne
(2.1D0,ido=11,18),& !Na~Ar
1.87D0,1.87D0,  (2.0D0,ido=21,27),1.63D0,1.40D0,1.39D0,  (1.87D0,ido=31,36),& !K ,Ca,  Ni~Zn(21~30),  Ga~Kr(31,37)
1.93D0,1.93D0,  (2.0D0,ido=39,45),1.63D0,1.72D0,1.58D0,  (1.93D0,ido=49,54),& !Rb,Sr,  Y ~Cd(39~48),  In~Xe(49~54)
1.96D0,1.96D0,  (2.0D0,ido=57,77),1.72D0,1.66D0,1.55D0,  (1.96D0,ido=81,86),& !Cs,Ba,  La~Hg(57~80),  Tl~Rn(81~86)
(2.0D0,ido=87,nelesupp) /) !Rn~Mt(87~109,~all)

!Covalent radius, from "Dalton Trans., 2008, 2832-2838", unit is Angstrom, will be convert to Bohr when multiwfn start
real*8 :: covr(0:nelesupp)=(/ 0.1D0,0.31D0,0.28D0,& !Ghost,H,Ne(1~2)
1.28D0,0.96D0,0.84D0,0.76D0,0.71D0,0.66D0,0.57D0,0.58D0,& !Li~Ne(3~10)     here C is sp3
1.66D0,1.41D0,1.21D0,1.11D0,1.07D0,1.05D0,1.02D0,1.06D0,& !Na~Ar(11~18)
2.03D0,1.76D0,1.70D0,1.60D0,1.53D0,1.39D0,1.39D0,1.32D0,1.26D0,& !K~Co(19~27)  here MnD0,FeD0,Co is low-spinD0, high spin is 1.61D0,1.52D0,1.50
1.24D0,1.32D0,1.22D0,1.22D0,1.20D0,1.19D0,1.20D0,1.20D0,1.16D0,& !Ni~Kr(28~36)
2.20D0,1.95D0,1.90D0,1.75D0,1.64D0,1.54D0,1.47D0,1.46D0,1.42D0,& !Rb~Rh(37~45)
1.39D0,1.45D0,1.44D0,1.42D0,1.39D0,1.39D0,1.38D0,1.39D0,1.40D0,& !Pd~Xe(46~54)
2.44D0,2.15D0,2.07D0,2.04D0,2.03D0,2.01D0,1.99D0,1.98D0,1.98D0,1.96D0,1.94D0,& !Cs~Tb(55~65)
1.92D0,1.92D0,1.89D0,1.90D0,1.87D0,1.87D0,1.75D0,1.70D0,1.62D0,1.51D0,1.44D0,1.41D0,& !Dy~Ir(66~77)
1.36D0,1.36D0,1.32D0,1.45D0,1.46D0,1.48D0,1.40D0,1.50D0,1.50D0,2.60D0,2.21D0,& !Pt~Ra(78~88)
2.15D0,2.06D0,2.00D0,1.96D0,1.90D0,1.87D0,1.80D0,1.69D0,(1.5D0,ido=97,nelesupp) /) !Ac~Cm(89~96),~all
!(Covalent) radius proposed by Suresh, from J. Phys. Chem. A 2001, 105, 5940-5944. For missing values (including all noble gases and very heavy elements), the ones in covr array are used
!Unit is Angstrom, will be convert to Bohr when multiwfn start
real*8 :: covr_Suresh(0:nelesupp)=(/ 0.1D0,0.327D0,0.28D0,& !Ghost,H,Ne(1~2)
1.219D0,0.911D0,0.793D0,0.766D0,0.699D0,0.658D0,0.633D0,0.58D0,& !Li~Ne(3~10)
1.545D0,1.333D0,1.199D0,1.123D0,1.11D0,1.071D0,1.039D0,1.06D0,& !Na~Ar(11~18)
1.978D0,1.745D0,1.337D0,1.274D0,1.236D0,1.128D0,1.18D0,1.091D0,1.089D0,& !K~Co(19~27)
1.077D0,1.146D0,1.187D0,1.199D0,1.179D0,1.209D0,1.201D0,1.201D0,1.16D0,& !Ni~Kr(28~36)
2.217D0,1.928D0,1.482D0,1.377D0,1.353D0,1.24D0,1.287D0,1.212D0,1.229D0,& !Rb~Rh(37~45)
1.24D0,1.362D0,1.429D0,1.385D0,1.38D0,1.421D0,1.4D0,1.397D0,1.40D0,& !Pd~Xe(46~54)
2.442D0,2.149D0,1.653D0,2.04D0,2.03D0,2.01D0,1.99D0,1.98D0,1.98D0,1.96D0,1.94D0,& !Cs~Tb(55~65)
1.92D0,1.92D0,1.89D0,1.90D0,1.87D0,1.87D0,1.364D0,1.346D0,1.256D0,1.258D0,1.222D0,1.227D0,& !Dy~Ir(66~77)
1.227D0,1.273D0,1.465D0,1.531D0,1.434D0,1.496D0,1.40D0,1.50D0,1.50D0,2.60D0,2.21D0,& !Pt~Ra(78~88)
2.15D0,2.06D0,2.00D0,1.96D0,1.90D0,1.87D0,1.80D0,1.69D0,(1.5D0,ido=97,nelesupp) /) !Ac~Cm(89~96),~all
!Covalent radius, from Pyykko "Chem. Eur. J.,15,186 (2009)", unit is in Angstrom, will be convert to Bohr when multiwfn start
real*8 :: covr_pyy(0:nelesupp)=(/ 0.1D0,0.32D0,0.46D0,& !Ghost,H,Ne(1~2)
1.33D0,1.02D0,0.85D0,0.75D0,0.71D0,0.63D0,0.64D0,0.67D0,& !Li~Ne(3~10)
1.55D0,1.39D0,1.26D0,1.16D0,1.11D0,1.03D0,0.99D0,0.96D0,& !Na~Ar(11~18)
1.96D0,1.71D0,1.48D0,1.36D0,1.34D0,1.22D0,1.19D0,1.16D0,1.11D0,& !K~Co(19~27)
1.10D0,1.12D0,1.18D0,1.24D0,1.21D0,1.21D0,1.16D0,1.14D0,1.17D0,& !Ni~Kr(28~36)
2.10D0,1.85D0,1.63D0,1.54D0,1.47D0,1.38D0,1.28D0,1.25D0,1.25D0,& !Rb~Rh(37~45)
1.20D0,1.28D0,1.36D0,1.42D0,1.40D0,1.40D0,1.36D0,1.33D0,1.31D0,& !Pd~Xe(46~54)
2.32D0,1.96D0,1.80D0,1.63D0,1.76D0,1.74D0,1.73D0,1.72D0,1.68D0,1.69D0,1.68D0,& !Cs~Tb(55~65)
1.67D0,1.66D0,1.65D0,1.64D0,1.70D0,1.62D0,1.52D0,1.46D0,1.37D0,1.31D0,1.29D0,1.22D0,& !Dy~Ir(66~77)
1.23D0,1.24D0,1.34D0,1.44D0,1.44D0,1.51D0,1.45D0,1.47D0,1.42D0,2.23D0,2.01D0,& !Pt~Ra(78~88)
1.86D0,1.75D0,1.69D0,1.70D0,1.71D0,1.72D0,1.66D0,1.66D0,1.68D0,1.68D0,1.65D0,1.67D0,1.73D0,1.76D0,1.61D0,& !Ac~Lr(89~103)
1.57D0,1.49D0,1.43D0,1.41D0,1.34D0,1.29D0,1.28D0,1.21D0,1.22D0,1.36D0,1.43D0,1.62D0,1.75D0,1.65D0,1.57D0,(1.5D0,ido=119,nelesupp)  /) !Rf~118(104~118),~all
real*8 :: covr_tianlu(0:nelesupp)=(/ 0.1D0,0.31D0,0.28D0,& !H,Ne(1~2) !Based on CSD radii, but for all main group (except for H,He), use IVA radius in corresponding row
(0.76D0,ido=3,10),& !Li~Ne(3~10)
(1.11D0,ido=11,18),1.2D0,1.2D0,& !Na~Ar(11~18),K,Ca
1.70D0,1.60D0,1.53D0,1.39D0,1.39D0,1.32D0,1.26D0,1.24D0,1.32D0,1.22D0,& !Sc~Zn(21~30)  here MnD0,FeD0,Co is low-spinD0, high spin is 1.61D0,1.52D0,1.50
(1.2D0,ido=31,36),1.42D0,1.42D0,& !Ga~Kr(31~36),Rb,Sr
1.90D0,1.75D0,1.64D0,1.54D0,1.47D0,1.46D0,1.42D0,1.39D0,1.45D0,1.44D0,& !Y~Cd(39~48)
(1.39D0,ido=49,54),1.46D0,1.46D0,& !In~Xe(49~54),Cs,Ba
2.07D0,2.04D0,2.03D0,2.01D0,1.99D0,1.98D0,1.98D0,1.96D0,1.94D0,1.92D0,1.92D0,1.89D0,1.90D0,1.87D0,1.87D0,& !La~Lu(57~71)
1.75D0,1.70D0,1.62D0,1.51D0,1.44D0,1.41D0,1.36D0,1.36D0,1.32D0,& !Hf~Hg(72~80)
(1.46D0,ido=81,86),1.46D0,1.46D0,&!Tl~Rn(81~86),Fr(still 1.46),Ra(still 1.46)
2.15D0,2.06D0,2.00D0,1.96D0,1.90D0,1.87D0,1.80D0,1.69D0,(1.5D0,ido=97,nelesupp) /) !Ac~Cm(89~96),~all
!Radii proposed in Chem. Phys. Lett., 480 (2009) 127-131, the unit is Bohr!
real*8 :: radii_Hugo(0:nelesupp)=(/ 0.10D0,1.00D0,0.74D0,& !Ghost,H,Ne(1~2)
1.59D0,1.21D0,1.28D0,1.10D0,0.97D0,1.00D0,0.88D0,0.79D0,& !Li~Ne(3~10)
1.63D0,1.33D0,1.51D0,1.29D0,1.14D0,1.15D0,1.02D0,0.93D0,& !Na~Ar(11~18)
1.77D0,1.49D0,1.44D0,1.41D0,1.42D0,1.42D0,1.35D0,1.31D0,1.31D0,1.33D0,1.33D0,1.20D0,1.51D0,1.31D0,1.18D0,1.18D0,1.07D0,0.99D0,& !K~Kr
1.80D0,1.55D0,1.48D0,1.43D0,1.42D0,1.38D0,1.37D0,1.36D0,1.35D0,1.28D0,1.34D0,1.23D0,1.53D0,1.36D0,1.26D0,1.23D0,1.14D0,1.06D0,1.87D0,1.62D0,& !Rb~Xe,Cs,Ba
1.56D0,1.57D0,1.58D0,1.57D0,1.56D0,1.55D0,1.55D0,1.49D0,1.52D0,1.51D0,1.50D0,1.49D0,1.48D0,1.47D0,1.58D0,& !La~Lu
1.41D0,1.34D0,1.31D0,1.32D0,1.27D0,1.23D0,1.23D0,1.21D0,1.14D0,1.49D0,1.35D0,1.37D0,1.27D0,1.21D0,1.12D0,1.83D0,1.16D0,& !Hf~Rn,Fr,Ra
1.62D0,1.47D0,1.52D0,1.48D0,1.47D0,1.50D0,1.51D0,1.51D0,1.48D0,1.47D0,1.46D0,1.45D0,1.44D0,1.43D0,1.67D0,1.51D0,(1.5D0,ido=105,nelesupp) /) !Ac~Rf,~all

real*8 :: YWTatomcoeff(18,3)=reshape((/ & !Coef. of fitting B3LYP/6-31G* density by Weitao Yang group for the first three rows, see supporting info. of JACS,132,6498
0.2815D0,2.437D0,11.84D0,31.34D0,67.82D0,120.2D0,190.9D0,289.5D0,406.3D0,561.3D0,760.8D0,1016.0D0,1319.0D0,1658.0D0,2042.0D0,2501.0D0,3024.0D0,3625.0D0, &
0.0D0,0.0D0,0.06332D0,0.3694D0,0.8527D0,1.172D0,2.247D0,2.879D0,3.049D0,6.984D0,22.42D0,37.17D0,57.95D0,87.16D0,115.7D0,158.0D0,205.5D0,260.0D0, &
0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.06358D0,0.3331D0,0.8878D0,0.7888D0,1.465D0,2.17D0,3.369D0,5.211D0 /),(/18,3/))
real*8 :: YWTatomexp(18,3)=reshape((/ & !Corresponding exponent of YWTatom, the value setted to 1.0 don't have any meaning, only for avoiding divide zero
0.5288D0,0.3379D0,0.1912D0,0.139D0,0.1059D0,0.0884D0,0.0767D0,0.0669D0,0.0608D0,0.0549D0,0.0496D0,0.0449D0,0.0411D0,0.0382D0,0.0358D0,0.0335D0,0.0315D0,0.0296D0, &
1.0D0,1.0D0,0.9992D0,0.6945D0,0.53D0,0.548D0,0.4532D0,0.3974D0,0.3994D0,0.3447D0,0.2511D0,0.215D0,0.1874D0,0.1654D0,0.1509D0,0.1369D0,0.1259D0,0.1168D0, &
1.0D0,1.0D0,1.0D0,1.0D0,1.0D0,1.0D0,1.0D0,1.0D0,1.0D0,1.0D0,1.0236D0,0.7753D0,0.5962D0,0.6995D0,0.5851D0,0.5149D0,0.4974D0,0.4412D0 /),(/18,3/))
!Atomic weights, from http://www.chem.qmul.ac.uk/iupac/AtWt/, the data is mainly based on Pure Appl. Chem., 81, 2131-2156 (2009)
real*8 :: atmwei(0:nelesupp)=(/ 0D0,1.00794D0,4.0026D0,6.941D0,9.01218D0,10.811D0,12.0107D0,14.0067D0,15.9994D0,18.9984D0,20.1797D0,& !1~10
22.98977D0,24.305D0,26.98154D0,28.0855D0,30.97376D0,32.065D0,35.453D0,39.948D0,39.0983D0,40.078D0,& !11~20
44.95591D0,47.867D0,50.9415D0,51.9961D0,54.93805D0,55.845D0,58.93319D0,58.6934D0,63.546D0,65.38D0,& !21~30
69.723D0,72.64D0,74.9216D0,78.96D0,79.904D0,83.798D0,85.4678D0,87.62D0,88.90585D0,91.224D0,& !31~40
92.90638D0,95.96D0,98D0,101.07D0,102.9055D0,106.42D0,107.8682D0,112.411D0,114.818D0,118.71D0,& !41~50
121.76D0,127.6D0,126.90447D0,131.293D0,132.90545D0,137.327D0,138.90547D0,140.116D0,140.90765D0,144.242D0,& !51~60
145D0,150.36D0,151.964D0,157.25D0,158.92535D0,162.5D0,164.93032D0,167.259D0,168.93421D0,173.054D0,& !61~70
174.9668D0,178.49D0,180.94788D0,183.84D0,186.207D0,190.23D0,192.217D0,195.084D0,196.96657D0,200.59D0,& !71~80
204.3833D0,207.2D0,208.9804D0,209D0,210D0,222D0,223D0,226D0,227D0,232.03806D0,& !71~90
231.03588D0,238.02891D0,237D0,244D0,243D0,247D0,247D0,251D0,252D0,257D0,258D0,259D0,262D0,265D0,268D0,271D0,272D0,270D0,276D0,& !91~109
281D0,282D0,285D0,285D0,289D0,289D0,293D0,294D0,294D0,& !110~118
(0D0,ido=119,nelesupp) /) !119~all
 
!Series of Lebedev-Laikov routines
integer :: Lebelist(32)=(/ 6,14,26,38,50,74,86,110,146,170,194,230,266,302,350,434,590,770,974,1202,1454,1730,2030,2354,2702,3074,3470,3890,4334,4802,5294,5810 /)
integer :: fact(0:10)=(/ 1,1,2,6,24,120,720,5040,40320,362880,3628800 /) ! Store factorials from 0~10 
integer :: isphergau=0 !By default, all basis functions are cartesian type, =1 means spherical (but some of them can be still cartesian type)
character*5 :: GTFtype2name(-32:56)=(/ &
"H 0  ","H+1  ","H-1  ","H+2  ","H-2  ","H+3  ","H-3  ","H+4  ","H-4  ","H+5  ","H-5  ", & !-32:-22
"G 0  ","G+1  ","G-1  ","G+2  ","G-2  ","G+3  ","G-3  ","G+4  ","G-4  ", & !-21:-13
"F 0  ","F+1  ","F-1  ","F+2  ","F-2  ","F+3  ","F-3  ","D 0  ","D+1  ","D-1  ","D+2  ","D-2  ", & !-12:-6,-5:-1
"     ","S    ","X    ","Y    ","Z    ","XX   ","YY   ","ZZ   ","XY   ","XZ   ","YZ   ", & !0~10
"XXX  ","YYY  ","ZZZ  ","XXY  ","XXZ  ","YYZ  ","XYY  ","XZZ  ","YZZ  ","XYZ  ", & !f 11~20
"ZZZZ ","YZZZ ","YYZZ ","YYYZ ","YYYY ","XZZZ ","XYZZ ","XYYZ ","XYYY ","XXZZ ","XXYZ ","XXYY ","XXXZ ","XXXY ","XXXX ", & !g 21~35
"ZZZZZ","YZZZZ","YYZZZ","YYYZZ","YYYYZ","YYYYY","XZZZZ","XYZZZ","XYYZZ","XYYYZ","XYYYY","XXZZZ","XXYZZ","XXYYZ","XXYYY","XXXZZ","XXXYZ","XXXYY","XXXXZ","XXXXY","XXXXX" /) !h 36~56
!Here s,p,d sequences are identical to .wfn, .wfx, .fch, .molden  !Note: Sequence in .fch = sequence in Gaussian
!Here f sequence is identical to .wfn, .wfx, but not identical to .fch and .molden
!Here g sequence is identical to .fch, .wfn does not support higher than f function, not identical to .wfx and .molden
!here h sequence is identical to .wfx and .fch, .molden doesn't support h
!Notice: The .wfn produced by G09 B.01 and later supports g and h, the definition is identical to here, and thus can be normally loaded
!Overall, spd: Multiwfn=wfn=wfx=fch=molden   f: Multiwfn=wfn=wfx!=fch=molden   g: Multiwfn=fch!=wfx=molden=molden2aim   h: Multiwfn=wfx=fch
integer :: type2ix(56)=(/ 0, 1,0,0, 2,0,0,1,1,0, 3,0,0,2,2,0,1,1,0,1, 0,0,0,0,0,1,1,1,1,2,2,2,3,3,4, 0,0,0,0,0,0,1,1,1,1,1,2,2,2,2,3,3,3,4,4,5 /)
integer :: type2iy(56)=(/ 0, 0,1,0, 0,2,0,1,0,1, 0,3,0,1,0,2,2,0,1,1, 0,1,2,3,4,0,1,2,3,0,1,2,0,1,0, 0,1,2,3,4,5,0,1,2,3,4,0,1,2,3,0,1,2,0,1,0 /)
integer :: type2iz(56)=(/ 0, 0,0,1, 0,0,2,0,1,1, 0,0,3,0,1,1,0,2,2,1, 4,3,2,1,0,3,2,1,0,2,1,0,1,0,0, 5,4,3,2,1,0,4,3,2,1,0,3,2,1,0,2,1,0,1,0,0 /)
!Negative value means the shell use spherical gauss function. -1=SP, and impossible be used in Multiwfn (when detect it, split it as S and P)
character :: shelltype2name(-5:5)=(/ "H","G","F","D"," ","S","P","D","F","G","H" /)
!Convert shell type to the number of orbitals,0=s,1=p,-1=sp,2=6d,-2=5d,3=10f,-3=7f,4=15g,-4=9g,5=21h,-5=11h
integer :: type2norb(-5:5)=(/ 11,9,7,5,4,1,3,6,10,15,21 /) 

!-------- Variables for wfn information(_org means using for backuping the first loaded molecule)
integer :: ibasmode=0 !0/1 = GTO/STO is used in current wavefunction
integer :: nmo=0,nprims=0,ncenter=0,ncenter_org=0,nmo_org=0,nprims_org=0 !Number of orbitals, primitive functions, nuclei
integer :: idxHOMO=0 !For fch and molden, record the index of original HOMO, this will be used to calculate linear response kernel for pi-electrons
integer :: ifiletype=0 !unknown textfile=0, fch=1, wfn=2, wfx=3, chg=4, pdb/xyz=5, .31=6, cube=7, grd=8, molden=9
integer :: wfntype=0 !0/1/2/3/4 means R/U/ROHF /R/U-Post-HF wavefunction
real*8 :: totenergy=0,virialratio=2,nelec=0,naelec=0,nbelec=0
integer :: ninnerelec=0 !Electrons represented by EDF
!-------- Variables for nuclei & basis function & Molecular orbital
type(atomtype),allocatable :: a(:),a_org(:),a_tmp(:)
type(primtype),allocatable :: b(:),b_org(:),b_tmp(:)
real*8,allocatable :: MOocc(:),MOocc_org(:),MOene(:) !Occupation number & energy of orbital
integer,allocatable :: MOtype(:) !The type of orbitals, (alpha&beta)=0/alpha=1/beta=2, not read from .wfn directly
character*4,allocatable :: MOsym(:) !The symmetry of orbitals, meaningful when .molden is used since it sometimes records irreducible representation
real*8,allocatable :: CO(:,:),CO_org(:,:),CO_tmp(:,:) !Coefficient matrix of primitive basis functions, including both normalization and contraction coefficients
                                                        !Note: Row/column of CO denote MO/GTF respectively, in contrary to convention
!-------- Describe inner electron density in EDF section in .wfx file
type(primtype),allocatable :: b_EDF(:)
real*8,allocatable :: CO_EDF(:)
integer :: nEDFprims=0
!-------- Variables when basis functions are basis rather than primitive function as basis
integer :: nbasis=0,nshell=0,nprimshell=0 !The number of basis, basis shell and primitive shell. SP shell is counted as S and P shell separately
integer,allocatable :: shtype(:),shcon(:),shcen(:) !Type, contraction degree and attributed center of a basis shell
real*8,allocatable :: primshexp(:),primshcoeff(:) !Exponent and contraction coefficient of a primitive shell
integer,allocatable :: basshell(:) !The ith element is the shell index that the ith basis attributed to
integer,allocatable :: bascen(:),bastype(:) !Center/type of basis, definition is the same as GTF
integer,allocatable :: basstart(:),basend(:) !The ith element means the basis from where to where is attributed to the ith atom
integer,allocatable :: primstart(:),primend(:)  !The ith element means the GTF from where to where is attributed to the ith basis function
real*8,allocatable :: primconnorm(:) !element i means the contract. coeff. * normalization coeff. of GTF i, can be used for e.g. constructing basis integral from GTF integral
real*8,allocatable :: Sbas(:,:),Sbas_org(:,:) !Overlap matrix and its backup
real*8,allocatable :: Dbas(:,:,:) !Dipole moment integral matrix, the first index 1,2,3=X,Y,Z
real*8,allocatable :: Tbas(:,:) !Kinetic energy integral matrix
real*8,allocatable :: Vbas(:,:) !Nuclear attraction potential integral matrix
real*8,allocatable :: Velbas(:,:,:) !Velocity integral matrix, the first index 1,2,3=X,Y,Z
real*8,allocatable :: Magbas(:,:,:) !Magnetic integral matrix, the first index 1,2,3=X,Y,Z
!Coefficient matrix for alpha/beta orbital, CObasa(i,j) means the coefficient of ith basis in the jth orbital, differ to CO(:,:)
real*8,allocatable,target :: CObasa(:,:),CObasb(:,:) !wfntype==0,2,3 only allocate CObasa(nbasis,nmo), ==1,4 also allocate CObasb, dimension of both cobasa and cobasb would be (nbasis,nbasis)
real*8,allocatable,target :: CObasa_org(:,:),CObasb_org(:,:)
real*8,allocatable,target :: Ptot(:,:),Palpha(:,:),Pbeta(:,:) !Density matrix of total/alpha/beta, for wfntype==0.or.wfntype==3, only Ptot is filled, for others, all of Ptot,Palpha and Pbeta are filled
real*8,allocatable :: Palpha_org(:,:),Pbeta_org(:,:) !Backup P, e.g. for Wiberg bond order calculation
! real*8,allocatable :: twoPDM(:) !Store two-particle density matrix by one-dimension array, Not use currently

!-------- Trajectory
integer :: nframetraj=0 !The number of frames in the trajectory
real*8,allocatable :: traj(:,:,:) !traj(1/2/3,a,i) corresponds to x/y/z of the ath atom in frame i
!-------- Points loaded from external file
integer :: numextpt=0
real*8,allocatable :: extpt(:,:),extpttmp(:) !extpt(i,1:4) corresponds to X/Y/Z/value of point i, length unit is bohr. extpttmp only records function value



!!!!!!!!!!!!!!!!!!!!!! Parameter !!!!!!!!!!!!!!!!!!!!!!
!For passing Dislin main parent GUI and draw widget identifier
integer idissetlight1,idissetlight2,idissetlight3,idissetlight4,idissetlight5,idissetlight6,idissetlight7,idissetlight8,idissetlightall0,idissetlightall1
integer idisgraph,idiszoomin,idiszoomout,idisisosurscl,idisscrval,idisshowbothsign,idisshowisosur,idisshowdatarange,idisshowmol,iprogbar,idisisosursec
integer idisisosurquality,idisisosurnumpt,idisorbinfo2
integer idisshowatmlab,idisshowaxis,idisbondradius,idislabelsize,idisbondcrit,idisatmsize,idisshowpathlab !In draw mol GUI
integer idisshowattlab,idisdrawinternalbasin,idisattsize !Draw basin GUI
integer idisshow3n3,idisshow3n1,idisshow3p1,idisshow3p3,idisshowCPlab,idisshowpath,idisshowbassurf
integer idisshowlocminlab,idisshowlocmaxlab,idisshowlocminpos,idisshowlocmaxpos !For molecular surface analysis
!For setting isosurface style, colors
integer idisisosur1style,idisisosur1solid,idisisosur1mesh,idisisosur1point,idisisosur1solidmesh,idisisosur1tpr,idisisosur1opa
integer idisisosur2style,idisisosur2solid,idisisosur2mesh,idisisosur2point,idisisosur2solidmesh,idisisosur2tpr,idisisosur2opa
integer idisisosurallstyle,idisisosurallsolid,idisisosurallmesh,idisisosurallpoint,idisisosurallsolidmesh,idisisosuralltpr
integer GUI_mode !=1: Show mol and orbitals =2: show plane =3: show isosurface =4: show mol and CPs =5: show mol and surface analysis extremes =6: Show basin integral

!Plotting external parameter, can be set in settings.ini
character :: graphformat*4="png " ! ps/eps/pdf/wmf/gif/tiff/bmp
integer :: graph1Dwidth=1280,graph1Dheight=800,graph2Dwidth=1280,graph2Dheight=1200,graph3Dwidth=1400,graph3Dheight=1400
integer :: itickreverse=0,iticks=2,symbolsize=8,ilenunit1D=1,ilenunit2D=1,iatmlabtype=1,iatmlabtype3D=2,iplaneextdata=0
integer :: numdigx=2,numdigy=2,numdigz=3,numdiglinex=3,numdigliney=3,numdigctr=3
real*8 :: planestpx=1.5D0,planestpy=1.5D0,planestpz=0.1D0
integer :: fillcoloritpx=5,fillcoloritpy=3,pleatmlabsize=50
real*8 :: disshowlabel=0.5D0
real*8 :: bondclrR=0.1D0,bondclrG=1.0D0,bondclrB=0.1D0,atmlabclrR=0D0,atmlabclrG=0D0,atmlabclrB=0D0
real*8 :: CP3n3RGB(3)=(/0.72D0,0D0,0.72D0/),CP3n1RGB(3)=(/1D0,0.5D0,0D0/),CP3p1RGB(3)=(/1D0,1D0,0D0/),CP3p3RGB(3)=(/0D0,1D0,0D0/)
real*8 :: atm3Dclr(0:nelesupp,3) !Colors of the atom spheres shown in 3D plots, set in "loadsetting" routine

!Plotting Internal parameter
integer :: imodlayout=0
integer :: idrawbasinidx=0,idrawinternalbasin=0 !Draw which basin. If draw interal part of the basin
integer :: ifixorbsign=0 !if 1, during generating orbital isosurface by drawmolgui, most part will always be positive (namely if sum(cubmat)<0 or sum(cubmattmp)<0, the data sign will be inverted)
integer :: iatom_on_contour,iatom_on_contour_far=0,plesel,IGRAD_ARROW=0,ILABEL_ON_CONTOUR,LASTCTRVAL,ictrlabsize=20,ivdwctrlabsize=0,iwidthvdwctr=10,iwidthposctr=1,iwidthnegctr=1,iwidthgradline=1
integer :: iclrindctrpos=5,iclrindctrneg=5,ivdwclrindctr=3,iclrindgradline=6,vdwctrstyle(2)=(/ 1,0 /),ctrposstyle(2)=(/ 1,0 /),ctrnegstyle(2)=(/ 10,15 /)
integer :: isavepic=0,icurve_vertlinex=0,iclrindatmlab=1,imarkrefpos=0,ilog10y=0,iclrcurve=1,inowhiteblack=0
integer :: ifragcontri,nfragatmnum,nfragatmnumbackup,ipromol,inucespplot=0,idrawmol=1,idrawisosur=0,isosursec=0,idrawtype=1,idrawcontour=1
integer :: iinvgradvec=0,icolorvecfield=0,vecclrind=30,idrawplanevdwctr=0,iplaneoutall=0
real*8 :: surcolorzmin,surcolorzmax !fillctr is the contour value will be draw on fillcolor map
real*8 :: curve_vertlinex=0D0,curvexyratio=0.618D0 !Gold partition
real*8 :: gradplotstep=0.002D0,gradplotdis=0.01D0,gradplottest=0.2D0,cutgradvec=0.3D0
real*8 :: clrRcub1same=0.3D0,clrGcub1same=0.75D0,clrBcub1same=0.3D0,clrRcub1oppo=0.3D0,clrGcub1oppo=0.45D0,clrBcub1oppo=0.9D0 !Color for isosurface 1 with solid style
real*8 :: clrRcub2same=0.4D0,clrGcub2same=0.5D0,clrBcub2same=0.0D0,clrRcub2oppo=0.35D0,clrGcub2oppo=0.1D0,clrBcub2oppo=0.9D0 !Color for isosurface 2 with solid style
real*8 :: clrRcub1samemeshpt=0.3D0,clrGcub1samemeshpt=0.75D0,clrBcub1samemeshpt=0.3D0,clrRcub1oppomeshpt=0.3D0,clrGcub1oppomeshpt=0.45D0,clrBcub1oppomeshpt=0.9D0 !Color for isosurface 1 with solid style
real*8 :: clrRcub2samemeshpt=0.4D0,clrGcub2samemeshpt=0.5D0,clrBcub2samemeshpt=0.0D0,clrRcub2oppomeshpt=0.35D0,clrGcub2oppomeshpt=0.1D0,clrBcub2oppomeshpt=0.9D0 !Color for isosurface 2 with solid style
real*8 :: opacitycub1=0.7D0,opacitycub2=0.7D0 !Opacity for isosurface 1 and 2 with transparent style
!About topology information on plane
integer :: imark3n3=1,imark3n1=1,imark3p1=1,imark3p3=1,imarkpath=1,sizemarkcp=30,sizemarkpath=5,sizemark3n1path=5,idrawintbasple=0,isurfstyle=2
real*8 :: clrRpath=0.3D0,clrGpath=0.1D0,clrBpath=0D0,clrR3n1path=0.0D0,clrG3n1path=0D0,clrB3n1path=0.5D0
integer,allocatable :: boldlinelist(:)
character*3 :: drawsurmesh="ON "
!Parameter for drawing molecular structure or 3D map
integer :: ienablelight1=1,ienablelight2=1,ienablelight3=1,ienablelight4=0,ienablelight5=0,ienablelight6=0,ienablelight7=0,ienablelight8=0 !If enable lighting 1~8
integer :: ishowatmlab=1,ishowCPlab=0,ishowpathlab=0,ishowaxis=1,isosurshowboth=1,ishowdatarange=0,idrawpath=1,idrawbassurf=1,ishowattlab=0,ishowatt=0
integer :: isosur1style=1,isosur2style=1 !isosurface style,1/2/3/4/5=solid,mesh,points,solid+mesh,transparent
integer :: ishowlocminlab=0,ishowlocmaxlab=0,ishowlocminpos=0,ishowlocmaxpos=0 !For molecular surface analysis
integer :: ishow3n3=0,ishow3n1=0,ishow3p1=0,ishow3p3=0
real*8 :: bondcrit=1.15D0,textheigh=30.0D0,ratioatmsphere=1.0D0,ratioCPsphere=0.8D0,bondradius=0.2D0,attsphsize=0.1D0
real*8 :: XVU=150.0D0,YVU=30.0D0,ZVU=6.0D0 !3D view angle
!For passing ploting parameter from GUI routine to their call-back routine
!sur_value: The value of isosurface will be plot by drawmol routine when idrawisosur=1
real*8 :: dp_init1,dp_end1,dp_init2,dp_end2,dp_init3,dp_end3,sur_value=0.05D0

!!! Other external parameter !!!
integer :: iautointgrid=1,radpot=75,sphpot=434 !sphpot=230/302/434/590/770, low is 50*434, high is 100*590
integer :: ispecial=0 !=0: Normal, =1 specific for Chunying Rong, =2 for Shubin's 2nd project
integer :: igenDbas=0,igenMagbas=0,igenP=1,iwfntmptype=1,outmedinfo=0,intmolcust=0,isilent=0,isys=2,iopengl=0,idelvirorb=1,ifchprog=1 !isys=1/2/3: Windows/Linux/MacOS X
integer :: iuserfunc=0,iDFTxcsel=84,ispheratm=1,ADCtransfer=0,SpherIVgroup=0,MCvolmethod=2,readEDF=1,ireadatmEDF=0,ishowptESP=1,imolsurparmode=1
integer :: NICSnptlim=8000
real*8 :: bndordthres=0.05D0,compthres=0.5D0,compthresCDA=1D0,expcutoff=-40D0,espprecutoff=0D0
integer :: nthreads=2,ompstacksize=100000000
character :: lastfile*200="",gaupath*80=""
!About function calculation, external or internal parameters
integer :: RDG_addminimal=1,ELF_addminimal=1,num1Dpoints=3000,atomdenscut=1,nprevorbgrid=120000,paircorrtype=3,pairfunctype=1,srcfuncmode=1
integer :: ELFLOL_type=0,ipolarpara=0,iALIEdecomp=0,iskipnuc=0
real*8 :: laplfac=1D0,ELFLOL_cut=0D0
real*8 :: RDG_maxrho=0.05D0,RDGprodens_maxrho=0.1D0,aug1D=1.5D0,aug2D=4.5D0,aug3D=6.0D0,radcut=10.0D0
real*8 :: refx=0D0,refy=0D0,refz=0D0
real*8 :: pleA=0D0,pleB=0D0,pleC=0D0,pleD=0D0 !!ABCD of the plane defined by main function 1000, used for special aims
real*8 :: globaltmp=0 !A variable can be used anywhere and can be set by option 5 of main function 1000, for debugging purpose avoiding re-compile code
!About line/plane/grid calculation, inner parameter
!For 3D grid data
real*8 :: orgx,orgy,orgz,endx,endy,endz,dx,dy,dz !Origin, end point and translation length
integer :: nx=80,ny=80,nz=80 !The number of points in three directions
!For 2D plane map
real*8 :: v1x,v1y,v2x,v2y,v1z,v2z,a1x,a1y,a1z,a2x,a2y,a2z,a3x,a3y,a3z,d1,d2 !Translation vector 1 and 2, three point in self-defined plane for projecting label, d1,d2=Length of v1,v2
real*8 :: orgx2D,orgy2D,orgz2D !Origin
integer :: ngridnum1=100,ngridnum2=100 !The number of points in two directions
!Specific for Shubin's project
real*8 :: steric_addminimal=1D-4,steric_potcutrho=0D0,steric_potcons=0D0
!Other
integer :: ifirstMultiwfn=1 !If 1, means we re-load file via main function -11 and don't need to do some initializations

contains
  integer function rtNThreads()
    IMPLICIT NONE
    INTEGER currNThreads, TID, OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
    print *,'TRUE'
!$OMP PARALLEL
      TID = OMP_GET_THREAD_NUM()
      !Only master thread does this
      IF (TID .EQ. 0) THEN
        currNThreads = OMP_GET_NUM_THREADS()
        PRINT *, 'OMP_NUM_THREADS = ', currNThreads
      END IF
!$OMP END PARALLEL
    IF (nthreads .NE. 0) THEN
      currNThreads = nthreads
    END IF
    rtNThreads = currNThreads
    PRINT *, 'Number of threads = ', rtNThreads
  END FUNCTION
end module


!-------- Module for topology analysis
module topo
integer,parameter :: maxnumcp=100000 !Maximum number of CPs
integer :: CPtype(maxnumcp)=0 !0=none 1=(3,-3) 2=(3,-1) 3=(3,+1) 4=(3,+3)
integer pathnumpt(maxnumcp) !pathnumpt: How many points in each path
integer :: topomaxcyc=120,ishowsearchlevel=0,maxpathpt=451,npathtry=30
integer :: numcp=0,numpath=0,ifunctopo=1 !Selected real space function number
real*8 :: CPpos(3,maxnumcp)
real*8,allocatable :: topopath(:,:,:) !Store topological paths. {x,y,z}, index of points in each path, index of paths
real*8 :: CPstepscale=1D0,gradconv=1D-6,dispconv=1D-7,minicpdis=0.03D0,vdwsumcrit=1.2D0,discritpathfin=0.05D0,pathstepsize=0.03D0,singularcrit=5D-22
real*8 :: CPsearchlow=0D0,CPsearchhigh=0D0
character :: CPtyp2lab(0:4)*6=(/ "  ??  ","(3,-3)","(3,-1)","(3,+1)","(3,+3)" /) 
!Basin surface related:
integer :: nsurfpt=100,nsurfpathpercp=40,numbassurf=0 !Number of points in each surface path, number of paths in each interbasin surface, total number of basin surfaces
integer :: cp2surf(100000)=0 !Convert total index of (3,-1) to surface index, if zero, means no surface corresponds to this CP
real*8 :: surfpathstpsiz=0.02D0 !Step size in interbasin surface path
real*8,allocatable :: bassurpath(:,:,:,:) !Store interbasin paths. {x,y,z}, indices of points, index of path, index of (3,-1) CP
!Interbasin path on plane map (Note: Two direction paths are counted as one path
integer :: nple3n1path=0 !Already generated in-plane path from (3,-1)
integer :: n3n1plept=300 !Number of points in each direction of path 
integer :: cp2ple3n1path(10000)=0 !Convert total index of (3,-1) on the given plane to interbasin path index, if zero, means no path corresponds to this CP
real*8 :: ple3n1pathstpsiz=0.02D0
real*8,allocatable :: ple3n1path(:,:,:,:) !Store path derived form (3,-1) on given plane. {x,y,z}, indices of points, direction path (1 or 2), index of (3,-1) CP on plane
end module


!--------- Module for surface analysis
module surfvertex
use defvar
type triangtype
    integer idx(3) !Consists of which three surface vertices
    real*8 area
    real*8 value !mapped function value at geometry center
end type
type surfcor2vtxtype
!if k=surfcor2vtxtype(i,q)%athcor means the two corner with surface corner index of i and k, interpolated to the surface vertex with index of %itpvtx. this information is stored in slot q
    integer athcor !another corner
    integer itpvtx !interpolated to which surface vertex index
end type
type(surfcor2vtxtype),allocatable :: surfcor2vtx(:,:)
integer,allocatable :: surcor2vtxpos(:) !Will add new interpolation relationship to which slot of surfcor2vtx

type(triangtype),allocatable :: surtriang(:) !Record center of generated surface triangle
integer nsurtri !Temporary accummlated index of triangles in generating process
type(content),allocatable :: survtx(:) !Record x,y,z coordinate of surface vertex, with interpolated function value
integer,allocatable :: abs2suridx(:,:,:) ! Convert absolute indices of corners to surface vertex indices
integer,allocatable :: vtxconn(:,:) !(i,j)=k means the two surface vertices with index of i and j are connected, j is storage slot
integer,allocatable :: vtxconnpos(:) !records current slot range of vtxconn
real*8 surfisoval,tetravol0,tetravol1,tetravol2,tetravol3 !volume of interpolated tetradrons of type 1,2,3
integer nsurvtx,nsurlocmin,nsurlocmax
integer surlocmaxidx(10000),surlocminidx(10000) !Store indices of local minimum and maximum points. If =0, means this slot is empty or has been discarded
integer :: nbisec=3 !Do how many times bisection before linear interpolation
integer :: ifuncintp=1 !Use which real space function to do bisection interpolation
integer,allocatable :: elimvtx(:),elimtri(:)
end module


!---------- Module for basin integral
module basinintmod
integer vec26x(26),vec26y(26),vec26z(26)
real*8 len26(26)
real*8 :: valcritclus=0.005D0
integer ifuncbasin !Which function is used to partition the basin
integer :: mergeattdist=5
integer*2,allocatable :: gridbas(:,:,:) !Each grid belongs to which basin(attractor). -2=Boundary grids -1=Traveled to boundary grid, 0=Unassigned, x=basin index
integer numatt !The number of crude attractors after near-grid method
integer numrealatt !The number of actual attractors (the ones left after clustering)
integer*2,allocatable :: attgrid(:,:) !Crude attractor corresponds to which grid. attgrid(i,1/2/3)=The ix/iy/iz of the ith attractor
real*8,allocatable :: attval(:),attxyz(:,:) !Value and xyz coordinate of crude attractors, attxyz(numatt,1:3)
integer,allocatable :: attconv(:) !Attractor conversion list. If attconv(i)=j, means attractor i is belong to actual attractor j. -1 and 0 is also included
integer,allocatable :: nrealatthas(:) !nrealatthas(i)=m means actual attractor i has m crude attractors
integer,allocatable :: realatttable(:,:) !realatttable(i,j)=k means the jth member of the ith actual attractor is crude attractor k
real*8,allocatable :: realattval(:),realattxyz(:,:) !Value and xyz coordinate of actual attractors. For the ones having multiple crude attractors, these arrays record average value
logical*1,allocatable :: interbasgrid(:,:,:) !.true. means this is a boundary grid, else it is a internal grid
logical*1,allocatable :: grdposneg(:,:,:) !.true. means the value at this grid is positive, .false. means negative
end module
