c	Continuous Runoff and Erosion Hillslope model with DYnamic Surface properties CREHDYS
c	CREHDYS v4.0 : multiple corrections and adjustments were made compared to original  CREHDYS model 
c     The model does the job but:
c	 - some parts of the code are redundant and can be optimized
c	 - some variables remain for historical reasons and are not usefull
c	Developed by Eric Laloy and Charles Bielders, revised and adapted by Timothée Clement
c	For more information, please see:
c	Laloy and Bielders (JoH 2008, EJSS 2009a, EJSS 2009b)
c	Laloy et al. (JoH 2010)
c     Clement et al. (HyP 2025)

	implicit none ! implicit none statement is used to inhibit a very old feature of Fortran that by default treats all variables that start with the letters i, j, k, l, m and n as integers and all other variables as real arguments. Implicit None should always be used.

c	Block 1 - General variables, containing information about the pixel size, type and time step during events. BLOCK DATA statement identifies a subprogram that initializes variables and arrays in labeled common blocks
c
	real dxvert, dxdiag, dx2, cu, tcst, t0, t1 ! All variables must be predefined in these BLOCKS before there value can be set and adjusted
	real dxf
	real dstart
	real dt, dx
	common /gen2/ dt, dx ! Fortran 77 has no global variables, i.e. variables that are shared among several program units (subroutines). The only way to pass information between subroutines is to use the subroutine parameter list. Sometimes this is inconvenient, e.g., when many subroutines share a large set of parameters. In such cases one can use a common block. This is a way to specify that certain variables should be shared among certain subroutines.
	integer m2,ierr,moutlet
	character(20) plotname

c	Block 2 - Hill-slope characteristics, including flow routing informations

	integer flacc(100485), neff ! integer (discrete exact number) ARRAY (1 dimension, so vector) of length 100485 (maximum length (spatial dimension) for a simulation). neff=number of grid elements in the basin
	integer outlet, xycoord(400,400), noutlets ! outlet=flow accumulation number at the outlet. xy=matrix of x-y coordinates of grid elements. noutlets=number of outlets
	integer nsoil,i,j,ncrop,m, ! nsoil and ncrops = number of soils and crops. i,j = indices used in different cases, e.g. loops. n=total number of grid elements (cells)=nrow*ncol. m= cell index (grid element) in input file (1 dimension)
     . x, y, count
	integer north, northw,west,southw,south,southe,east,
     .      northe, fdn,fdne,fde,fdse,fds,fdsw,fdw,fdnw
	real sl(100485), fdir(100485),xcoord(100485)
	real ycoord(100485)
	integer spix(100485),cpix(100485)
	real numpix(100485)
	integer simdur, ncol, nrow, nrowout, n, n2
	common /gen/ simdur, ncol, nrow, nrowout, n, n2


c	Block 3 - Variables for rainfall, rainfall kinetic energy and climat

	integer raites,l
	
	real rcc(90000), rtc(90000), rain,pinst,tempc,soitemp,radi


	real tc(90000), rc(90000), dtm
	integer lline, idate
	common /rainfall/ tc, rc, dtm, lline, idate

	real TOTPRECTOT,TOTPRECGRTOT,EKTOT,EKTOTBEL,EKTOTBELCOV,EKTOTTOT
     .	,EKTOTBELTOT, EKTOTBELCOVTOT,TOTPREC,TOTPRECGR,EKBELCOVSPLASH, 
     .	EKBELCOV, EKBELCOV2,EKBEL,DAY,EK


c	Block 4 - Infiltration and drainage processes

	real tp(100485),fc(100485),wp(100485),a(100485),df(100485)
     .	,ksat(100485),dr(100485)

	real FIL,VERTFLUX,FILTREAL, TOTFILT, MEANFILT, starteventtime

	real thetainst(100485),theta(100485),psi(100485), thetaoutlet,
     .	ke(100485),pond(100485), asmlim(100485), lsttr ! because the dimension of cuminf is declared in the common block, the size of the array should not be specified here
	
	integer rrksovar

	real filt(100485),draintot(100485)
	real gamptun3,drain,keb2(100485)
	real kcover2(100485),cfact(100485)
	real ksend(100485),kinio,kinim,keff(100485)

	real kso2, kcho2, kchm2, hfront, c_kso


c	Block 5 - Evapotranspiration process
	real albed(100485) ! because the dimension of these variables are declared in the common block, the size of the arrays should not be specified here

	real cl(100485),sa(100485),st(100485)
	real wcf(100485), vfs(100485)

	real THETBEFETR,TETPR,PERCOL,TOTWLOSS

	real etpmm
	common /etpot/ etpmm

	real esu(100485),edx(100485),es1p(100485),es2p(100485),
     &     es(100485),pep(100485),cuminf(100485),
     &     alphas(100485),edisp,edispedx
	integer ie
	common /etpact/ esu,edx,es1p,es2p,
     &     es,pep,cuminf,
     &     alphas,edisp,edispedx,ie


c	Block 6 - Random roughness and surface storage

	real RRCM, RRMC

	real dir2(100485),rrv(100485),
     .	dirvon(100485),rr(100485),
     .	rrch,rrexp,c_rr, meandirvon, totdirvon,
     .    MB_stor, prop_nodep, connectresh

	real rrintini(10),rrmaizeini(10)


c	Block 7 - Runoff height

	real HF3,HMM,HINDT, MEANHMM, TOTHMM, TOTS, MEANS
	real hsv(100485),ss(100485),hsv1221(100485),s(100485)


c	Block 8 - Runoff dynamics and balance

	real alpha1221,qp1221,totruntot,totrun,qintot11,qintot12
	real QIDT,lat,QOUTTEST,QMMHR, TOTQMMHR, TOTQOUTLET
	real chwid,mnch,mn,qoutlet,weightr

	real q(100485),wheel(100485)

	real oldqbar(100485),oldqtrans,man2
	real qmin,qi2(100485),lstr,prstr,qi(100485),q21
	real q2t1(100485),q2t2(100485),q1t1(100485),q1t2(100485)
	integer t1step,lastru
	integer timnum,timnum2
	real qp12d,qp21d
	real q2,tol2,h11,h21,h12,h22,alp,testck,ck
	real qtest,hf,hf2,q22,hwp11,hwp12,hwp21
	real qp22(100485),wper11,wper12,wper21,wper22,hwp,hmin
c	common /kin1/ qp12d,qp21d
	real qin21,qin22,qin12,qin11,beta,omega
	common /kin/ qin21,qin22,qin12,qin11,beta,omega
	real qp21(100485), qp12(100485), qp11(100485)
	integer iw
	common /kin2/ qp21, qp12, qp11, iw
	real alp11(100485),alp12(100485),alp21(100485),alp22(100485)
	common /kin3/ alp11,alp12,alp21,alp22
	real slwp,qact,width2(100485)
	integer iw2
	common /hwpr/ slwp,qact,width2,iw2
	real mann(100485)
	common /manning/ mann


c	Block 9 - Crop rotation and cover, LAI and biomass development

	integer prstmaize
 	character(15) curcropper ! current cropping period : maize or intercrop
	real cover(100485)

	real kextint, kextmaize, gcoverint, kcovermaize, biomcover
	real rbiomcover,kccover,kext,prstlai,lstbiomcover
	real FN,FW,KDEG, prstlaimaize, prstlaiint

	real sumtemp, sumtempdeg, checkd
	real prstcover, prstcoverint, prstcovermaize
	
	integer dinterp(10), dmaizep(10),wheelinter(10),wheelmaize(10)  ! no more than 10 maize cropping seasons or intercrop periods
	integer ninterp, nmaizep, rrksvarint(10), rrksvarmaize(10)
	integer destruction, burial

 
c	Block 10 - Maize canpoy and stem flow

	real TRA, HEIGHT, BRANDT, TRF


c	Block 11 - Erosion processes

	real qs22(100485), qs21(100485),qs12(100485),qs11(100485)
	real qs2t1(100485), qs2t2(100485)

	real As,rhop,rhow,dclay,dsilt,DSAND,DVFS,DCF,CLAYFRACT
	real SILTFRACT,SANDFRACT,VFSFRACT,CFFRACT,NUDYNWAT,DIAM
	real VS,d50,CTC,DTC,COH,YDET,TOTEROSTOT, TOTEROS, CCOUTLET
	real TOTQSOUTLETS,QSINTOT11,QSINTOT12, DFL,TRANSC,DEP,EROS
	real Dsplash,VEL,WSTREAM,CONC,DFLMAX,EBAR,Lerosionsplash
	real COHW,YDETW,concmax

	
c	Block 12 - Routines for kinematic wave equation:
c	fully implicit four points finite difference approximation using Newton-Raphson iteration
c	internal: Iter and Iter2
c	external:
	external errq2,errq3,errhwp


c	START

	call cpu_time(t0)		! record the current time at start of simulation

c	open files

	open (unit=1, file='gen.inp', status='old')		!In fortran, each file is associated with a unique unit number (identifier), an integer, typically between 1 and 99. Some units are reserved: 5 is a standard input, 6 is standard output
	open (unit=2, file='weather.inp', status='old')
	open (unit=3, file='plot.inp', status='old')
c	open (unit=5, file='cover.out', status='unknown')
c	open (unit=399, file='filt.out', status='unknown')		! Why weird values for unit ?
	open (unit=375,file='daily.out',status='unknown')
c	open (unit=376,file='events.out',status='unknown')
	open (unit=378,file='dailyRE.out',status='unknown')
	open (unit=379,file='eventR.out',status='unknown') ! event runoff output
	open (unit=380,file='eventE.out',status='unknown')
 	open (unit=381,file='crop.out',status='unknown') ! (inter)crop growth and cover
 	open (unit=382, file='soil.out',status='unknown') ! soil outputs (rain, infiltration, drainage, etp, water content) 
	open (unit=9,file='parameros.inp',status='unknown')
	open (unit=10,file='paramhyd.inp',status='unknown')
c	open (unit=11,file='ksat.txt',status='unknown')
c	open (unit=13,file='rr.txt',status='unknown')
c	open (unit=16,file='man.txt',status='unknown')
c	open (unit=17,file='cover.txt',status='unknown')
	open (unit=18,file='theta.inp',status='unknown')

c	read input data
	read (1,*) plotname
	read (1,*) simdur
	read (1,*) dt        ! read(unit number identifier=1='gen.inp', format identifier=1, défined hereabove, followed by the name of variables to which to assign read values
    1 format (/,1x,i5/,1x,f6.0)		! define a format (# 1 in this case) for inputs and outputs (/= go to next line, i.e. introduce a new input , 1x=there is a space at the beginning, i5/=there is an integer with 5 number followed by a line break, 1x=there is a space again, and un real with 6 numbers and 0 decimal digits)
	read (1,*) nsoil		!read(unit identifier, * = pre-formatted reading format, e.g., values are separated by spaces or commas). Making a new "read" statement (after the preceding one)by going to a new line will also go the a new line in the input file.
	read (1,*) tp(1),fc(1),wp(1),a(1),df(1)		! Only the first value of theses elements is defined because only a single soil type is considered. Code should be adapted if spix>1
	read (1,*) cl(1),sa(1),st(1),wcf(1),vfs(1)	
	read (1,*) rrch, mnch
	read (1,*) albed(1)
	read (1,*) dx		![m]  
	read (1,*) ninterp, nmaizep ! not used anymore ?
	read (1,*) dinterp(1:ninterp)
	read (1,*) destruction, burial
	read (1,*) dmaizep(1:nmaizep)
	read (1,*) wheelinter(1:ninterp) , wheelmaize(1:nmaizep)
	read (1,*) rrintini(1:ninterp) 
	read (1,*) rrmaizeini(1:nmaizep)  
	read (1,*) kinio, kinim 
	read (1,*) rrksvarint(1:ninterp) 
	read (1,*) rrksvarmaize(1:nmaizep) 
	read (1,*) kextint, kextmaize 
	read (1,*) gcoverint, kcovermaize
	read (1,*) rbiomcover
	read (10,*) kso2, kcho2, kchm2, hfront,c_kso,		! hydrology parameters 
     .	c_rr, man2, MB_stor, prop_nodep, connectresh
  
	read (3,*) n,nrow,ncol		! n=total number of grid elements (cells)=nrow*ncol

	do i=1,n		! loop on all grid elements (cells)
	read (3,*) numpix(i),xcoord(i),ycoord(i),sl(i),		! record the plot inputs for all grid cells
     . 	fdir(i), ! flow direction
     .	flacc(i) ! flow accumulation
     .	,spix(i),cpix(i),wheel(i) ! indice for soil, crop, and channel (wheel track)
	end do	
	count=0
	do i=1,nrow		! 
		do j=1,ncol
		count=count+1
		xycoord(i,j)=count		! assign the 1-dimension index to each grid cell in the xycoord matrix
		end do
	end do	
	do i=1,n
c		read(11,*) ksat(i)
		ksat(i)=kso2		! spatially-homogeneous Ksat homogène spatialement, dead variable ? this ksat is no more used, it is now ke that is used, and computed as a function of cover etc.
	end do
	do i=1,n
c		read(13,*) rr(i)
		rr(i)=0		! rr is defined null everywhere initially. The value will be changed afterward according to the period (intercrop or cropping season) and to the presence of wheel tracks
	end do
	do i=1,n
c		read(16,*) mann(i)
		if (wheel(i).eq.0) mann(i)=man2		! spatially-homogeneous Manning's n
          if (wheel(i).eq.1) mann(i)=mnch
	end do
      do i=1,n
c		read(17,*) cover(i)
		cover(i)=0.0		! initial cover is null everywhere, whether the simulation begins in intercropping or cropping season period
	end do
	do i=1,n
		read(18,*) theta(i)		! initial soil moisture content [m³ water/m³soil] for each element, can be spatially heterogenous
      end do
            
	read(9,*) As, d50, Coh, Cohw 	! erodibility parameters 

	close(1)		! CLOSE statement disconnects a file from a unit. CLOSE(unit=1) : disconnect the "1" unit which was open by the "OPEN" statement
	close(3)
!	close(9)
	close(10)
	close(11)
	close(13)	
	close(16)
	close(17)		! Not all file units are closed because we still need some of them : weather, parameros, outputs,...

	dxvert=dx
	dxdiag=dx*2**(0.5)		! diagonal distance (pythagorean theorem)

	neff=0
	df(1)=df(1)*10 ! convert df [cm] to df [mm]
	do i=1,n		! loop on all grid elements (cells)
	  if (sl(i).gt.(-9999)) then
	    sl(i)=sl(i)/100	! slope introduced in %
	    neff=neff+1
	    thetainst(i)=theta(i)   ! copy-past theta initial values to thetainst
	    a(i)=a(spix(i))		! infiltration & drainage variables. Only 1 soil type possible so only the 1st value is defined
	    df(i)=df(spix(i)) ! all df(i) are same as df(1)
	    tp(i)=tp(spix(i))
	    wp(i)=wp(spix(i))
	    fc(i)=fc(spix(i))
	    dr(i)=0.
	    s(i)=0. ! total ponding water of a cell at a given time in event
	    albed(i) = albed(cpix(i))	
	    sa(i) = sa(spix(i))
	    cl(i)= cl(spix(i))
	    st(i)=st(spix(i))
	    ksend(i) = kso2		! no wheel track is considered initially
	    psi(i)=hfront
	    cfact(i)=c_kso
	  end if
	end do
	
	rrexp=c_rr
		
	dx2=dx*dx		! pixel surface [m²]
	cu=dx2/3.6e+6	! coefficient to convert (infiltration or other) flow in [mmhr-1] to [m3s-1]
	tcst=dt/3600	! fraction of dt [s] in 1 hour = coefficient to convert rain [mmhr-1] to [mm per time step]

	noutlets=0
	do i=1,n
	  if (flacc(i).eq.maxval(flacc)) then
	    moutlet=i ! moutlet is the index of (the last) outlet cell
	    noutlets=noutlets+1 ! number of outlets
	  end if
 	end do

	do i=1,n
	  if (sl(i).gt.(-9999)) then
	    wp(i)=wp(i)*tp(i) ! !!! Conversion from wp [m³water/m³pores] to wp [m³water/m³soil]
	    fc(i)=fc(i)*tp(i) ! idem
	  end if
	end do	

c	variables for evapotranspiration, Ritchie's method (1972)
c	upper stage of soil evaporation (mm)
	do i=1,n
	  if (sl(i).gt.(-9999)) then		! compute only if there is slope information for the cell
	    
	    alphas(i) = (4.165 + 0.02456 * (sa(i) + vfs(1))
     &            - 0.01703 * cl(i)
     &            - 0.0004 * (sa(i) + vfs(1))**2)		! compute soil tranmissivity for each cell
	    esu(i)= (9*(alphas(i)-3)**0.42)		! compute upper limit of soil evaporation for each cell [mm]

c	    maximum evaporative depth (mm)
	    edx(i)=90.-0.77*cl(i)+0.006*(sa(i)+vfs(1))**2 ! max evaporative depth [mm]

c	    intializing of the evaporation variables
	    es1p(i)=0.
	    es2p(i)=0.
	    es(i)=0.		!soil evaporation
	    pep(i)=0.
	    asmlim(i) = wp(i) !+ (fc(i) - wp(i))*.25
	  end if
	end do

c	Erosion variables pre-processing

c	Vs = settling velocity - stokes law spherical particle Chow et al. 1988
	rhop=2650 ! density of sediment particles [kg/m³]
	rhow=1000 ! density of water
	dclay=0.002 ! diameter of clay particles in [mm]
	dsilt=0.01
	dsand=0.2
	dvfs=0.03
	dcf=0.5
	clayfract=cl(1)/100			!0.15
	siltfract=st(1)/100			!0.25
	sandfract=sa(1)/100			!0.50
	vfsfract=vfs(1)/100			!0.10 ! very fine sand
	cffract=wcf(1)/100			!0 ! coarse fragments
	nudynwat=1.3e-3 ! mu (and not nu!) dyanmic viscosity of water at 10°C [Nsm-2]

	Diam=(dclay*clayfract+dsilt*siltfract+dsand*sandfract+dvfs
     .	*vfsfract+cffract*dcf)*1e-3  ! weighted mean particle diameter [m]

	Vs=(1./18)*((Diam)**2)*9.81*(rhop-rhow)/nudynwat ! Particle settling velocity [m/s]. Original CREHDYS : Vs=	10*(2./9)*((Diam*0.5)**2)*9.81*(rhop-rhow)/nudynwat. The *10 factor at the beginning has been removed
	
	ctc=((d50+5)/0.32)**(-0.6) ! experimental coefficient "c" in transport capacity calculation
	dtc=((d50+5)/300)**(0.25) ! experimental coefficient "d" in transport capacity calculation

!	Coh: 0.1 kg/cm2 = 9.81 kPa

	Ydet=1./(0.89+0.56*Coh)		! efficiency factor of flow detachment on classic (overland flow) cells (0 to 1)  [/]
	Ydetw=1./(0.89+0.56*Cohw)   ! efficiency factor of flow detachment on channel (whell track) cells (0 to 1) [/]

c	Cover degradation constants

	fN=0.57 ! correction factor accounting for the initial nitrogen content of the residue in Douglas & Rickman residue decomposition equation. fN=0.57+0.126*N where N is initial nitrogen content of residue [g/kg]
	fW=0.2 ! coefficient depending on whether the residue is buried or lies on soil surface (=0.2 if surface)
	kdeg=-0.01 ! cover crop decomposition coefficient [°C days]
	
c	Set a few variables
	
	totruntot =0.
	toterostot=0.
	totprectot=0.
	totprecgrtot=0.
	ektot=0.
	ektotbel=0.
	ektotbelcov=0.
	ektottot=0.			! van Dijk et al. 2002.
	ektotbeltot=0.		! Boline et al. 1984.
	ektotbelcovtot=0.

c	Start simulation
      dstart=0
	do 11 idate=1,simdur		! DO LOOP ON THE WHOLE SIMULATION DURATION (index = idate) (time step=1 day)

	write(*,*) 'working on day ',idate

	totrun=0
	toteros=0
	timnum=0
	timnum2=0
	totprec=0
	totprecgr=0.
	ektot=0.
	ektotbel=0.
	ektotbelcov=0.

      
c	Specific implementation for Laloy et al. (EJSS 2009b) calibration - validation case study: 

c	***** WHEEL / KS / RR *****
c---------------------------------------------------------------
	if (idate.eq.dinterp(1)) then		! first day inter1 ! original CREHDYS : if (idate.eq.1)
	write(*,*) ' Intercrop 1'
      curcropper='intercrop'
!	Kcover & Kext
	! kccover=gcoverint ! We are in the intercropping period so crop growth coefficient=cover crop growth coefficient
	! kext=kextint ! cover extinction coefficient (intercrop)
	prstmaize=0
	dstart=dinterp(1)


	! RR & KS & manning

		!	Ks & RR overland decreasing

		if (rrksvarint(1).eq.0) then

			rrksovar=0 ! If Ksat and RR are constant during the intercropping period

		else
			rrksovar=1 ! If Ksat and RR are varying during the intercropping period

	    end if

		do i=1,n
			if (sl(i).gt.(-9999)) then

			if (wheel(i).eq.0) then ! overland flow cell (not wheel track)

c	random roughness
			
			rr(i)=rrintini(1) ! [mm]
			rrcm = rr(i)/10.	! from rr[mm] to rrcm [cm]

c	surface retention, from OpenLISEM (Kamphorst et al., 2000 ?)

			dir2(i)= 0.243*rr(i)+0.01*rr(i)*rr(i)-0.012*rr(i)*sl(i)*100  ! previous CREHDYS: eq. from Mwandera & Feyen (1992):(0.294*rrcm+0.031*rrcm*rrcm-0.012*rrcm*sl(i)*100)*10.

			if (rrksovar.eq.0) then ! ! If Ksat and RR are constant during the season

				ke(i)=ksend(i)
			else
				ke(i)=kinio
			end if

			else ! wheel track cell
			rr(i)=rrch
			rrmc=rrch/10.

			dir2(i)= 0.243*rr(i)+0.01*rr(i)*rr(i)-0.012*rr(i)*sl(i)*100  ! previous CREHDYS: eq. from Mwandera & Feyen (1992):(0.294*rrcm+0.031*rrcm*rrcm-0.012*rrcm*sl(i)*100)*10.
		   
			ke(i)=kcho2

			end if
			end if
		end do


	! initialize sumtemp

		sumtemp=0
		sumtempdeg=0

	! initialize Ek

	ektot=0.
	ektotbel=0.
	ektotbelcov=0.
	ektottot=0.			! van Dijk et al. (2002)
	ektotbeltot=0.		! Bollinne et al. (1984)
	ektotbelcovtot=0.
      end if
	

c----------------------------------------------------------
	if (idate.eq.dmaizep(1)) then		! maize1 ! original CREHDYS : elseif
	write(*,*) ' Maize 1'
      curcropper='maize'
	!	Kcover & Kext
	! kccover=kcovermaize ! We are in the maize cropping season so the cover growth coefficient is the one from maize
	! kext=kextmaize ! cover extinction coefficient (maize)
	prstmaize=1
	dstart=dmaizep(1)



	! RR & KS & manning

		!	Ks & RR overland decreasing

		if (rrksvarmaize(1).eq.0) then

			rrksovar=0 ! If Ksat and RR are constant during the season

		else
			rrksovar=1 ! If Ksat and RR are varying during the season

	    end if

		do i=1,n
		if (sl(i).gt.(-9999)) then

			if (wheel(i).eq.0) then
			
			rr(i)=rrmaizeini(1)
			rrcm = rr(i)/10.	! rrm in cm

			dir2(i)= 0.243*rr(i)+0.01*rr(i)*rr(i)-0.012*rr(i)*sl(i)*100  ! previous CREHDYS: eq. from Mwandera & Feyen (1992):(0.294*rrcm+0.031*rrcm*rrcm-0.012*rrcm*sl(i)*100)*10.

			if (rrksovar.eq.0) then ! If Ksat and RR are constant during the season

				ke(i)=ksend(i)
			else
				ke(i)=kinim
			end if

			else
			rr(i)=rrch
			rrmc=rrch/10.

			dir2(i)= 0.243*rr(i)+0.01*rr(i)*rr(i)-0.012*rr(i)*sl(i)*100  ! previous CREHDYS: eq. from Mwandera & Feyen (1992):(0.294*rrcm+0.031*rrcm*rrcm-0.012*rrcm*sl(i)*100)*10.

		   
			ke(i)=kchm2
              end if
			end if
		end do



	! initialize sumtemp

		sumtemp=0 ! do not initialize sumtempdeg at zero or it would reinitialize intercrop degradation at maize sowing

	! initialize Ek

	ektot=0.
	ektotbel=0.
	ektotbelcov=0.
	ektottot=0.			
	ektotbeltot=0.	
	ektotbelcovtot=0.
      end if
	

c-------------------------------------------------------------------------------
	if (idate.eq.dinterp(2)) then  !inter2 ! original CREHDYS : elseif
	write(*,*) ' Intercrop 2'
      curcropper='intercrop'
	!	Kcover & Kext
	! kccover=gcoverint
	! kext=kextint ! cover extinction coefficient (intercrop)
	prstmaize=0
	dstart=dinterp(2)


	! RR & KS & manning

		!	Ks & RR overland decreasing

		if (rrksvarint(2).eq.0) then

			rrksovar=0 ! ! If Ksat and RR are varying during the intercropping period

		else
			rrksovar=1 ! Ksat and RR are varying during the intercropping period

	    end if
		do i=1,n
          if (sl(i).gt.(-9999)) then

			if (wheel(i).eq.0) then
			
			rr(i)=rrintini(2)
			rrcm = rr(i)/10.	! rrm in cm

			dir2(i)= 0.243*rr(i)+0.01*rr(i)*rr(i)-0.012*rr(i)*sl(i)*100  ! previous CREHDYS: eq. from Mwandera & Feyen (1992):(0.294*rrcm+0.031*rrcm*rrcm-0.012*rrcm*sl(i)*100)*10.

			if (rrksovar.eq.0) then ! If Ksat and RR are constant during the season

				ke(i)=ksend(i)
			else
				ke(i)=kinio
			end if

			else
			rr(i)=rrch
			rrmc=rrch/10.

			dir2(i)= 0.243*rr(i)+0.01*rr(i)*rr(i)-0.012*rr(i)*sl(i)*100  ! previous CREHDYS: eq. from Mwandera & Feyen (1992):(0.294*rrcm+0.031*rrcm*rrcm-0.012*rrcm*sl(i)*100)*10.

		   
			ke(i)=kcho2

			end if
	 	    end if
		end do



	! initialize sumtemp

		sumtemp=0
		sumtempdeg=0

	! initialize Ek

	ektot=0.
	ektotbel=0.
	ektotbelcov=0.
	ektottot=0.			
	ektotbeltot=0.		
	ektotbelcovtot=0.
      end if
	
c--------------------------------------------------------------------------
	if (idate.eq.dmaizep(2)) then	! maize2 ! original CREHDYS : elseif
	write(*,*) ' Maize 2'
      curcropper='maize'
	!	Kcover & Kext
	! kccover=kcovermaize
	! kext=kextmaize ! cover extinction coefficient (maize)
	prstmaize=1
	dstart=dmaizep(2)

	! RR & KS & manning
		!	Ks & RR overland decreasing

		if (rrksvarmaize(2).eq.0) then

			rrksovar=0 ! ! If Ksat and RR are constant during the season

		else
			rrksovar=1 ! Ksat and RR are varying during the season

	    end if
		do i=1,n
		if (sl(i).gt.(-9999)) then

			if (wheel(i).eq.0) then
		
			rr(i)=rrmaizeini(2)
			rrcm = rr(i)/10.	! rrm in cm

			dir2(i)= 0.243*rr(i)+0.01*rr(i)*rr(i)-0.012*rr(i)*sl(i)*100  ! previous CREHDYS: eq. from Mwandera & Feyen (1992):(0.294*rrcm+0.031*rrcm*rrcm-0.012*rrcm*sl(i)*100)*10.

			if (rrksovar.eq.0) then ! If Ksat and RR are constant during the season

				ke(i)=ksend(i)
			else
				ke(i)=kinim
			end if

			else
			rr(i)=rrch
			rrmc=rrch/10.

			dir2(i)= 0.243*rr(i)+0.01*rr(i)*rr(i)-0.012*rr(i)*sl(i)*100  ! previous CREHDYS: eq. from Mwandera & Feyen (1992):(0.294*rrcm+0.031*rrcm*rrcm-0.012*rrcm*sl(i)*100)*10.
		   
			ke(i)=kchm2

			end if
	        end if
		end do

	! initialize sumtemp

		sumtemp=0
		sumtempdeg=0

	! initialize Ek

	ektot=0.
	ektotbel=0.
	ektotbelcov=0.
	ektottot=0.			
	ektotbeltot=0.		
	ektotbelcovtot=0.
      end if
c-------------------------------------------------------------------------------
	if (idate.eq.dinterp(3)) then  !inter3 ! original CREHDYS : elseif
	write(*,*) ' Intercrop 3'
      curcropper='intercrop'
	!	Kcover & Kext
	! kccover=gcoverint
	! kext=kextint ! cover extinction coefficient (intercrop)
	prstmaize=0
	dstart=dinterp(3)

	! RR & KS & manning

		!	Ks & RR overland decreasing

		if (rrksvarint(3).eq.0) then

			rrksovar=0 ! ! Ksat and RR are constant during the season

		else
			rrksovar=1 

	    end if
		do i=1,n ! do loop on all grid elements (cells)
		if (sl(i).gt.(-9999)) then

			if (wheel(i).eq.0) then
			
			rr(i)=rrintini(3)
			rrcm = rr(i)/10.	! rrm in cm

			dir2(i)= 0.243*rr(i)+0.01*rr(i)*rr(i)-0.012*rr(i)*sl(i)*100  ! previous CREHDYS: eq. from Mwandera & Feyen (1992):(0.294*rrcm+0.031*rrcm*rrcm-0.012*rrcm*sl(i)*100)*10.

			if (rrksovar.eq.0) then ! If Ksat and RR are constant during the season

				ke(i)=ksend(i)
			else
				ke(i)=kinio
			end if

			else
			rr(i)=rrch
			rrmc=rrch/10.

			dir2(i)= 0.243*rr(i)+0.01*rr(i)*rr(i)-0.012*rr(i)*sl(i)*100  ! previous CREHDYS: eq. from Mwandera & Feyen (1992):(0.294*rrcm+0.031*rrcm*rrcm-0.012*rrcm*sl(i)*100)*10.

		   
			ke(i)=kcho2

			end if
	        end if
		end do



	! initialize sumtemp

		sumtemp=0
		sumtempdeg=0

	! initialize Ek

	ektot=0.
	ektotbel=0.
	ektotbelcov=0.
	ektottot=0.		
	ektotbeltot=0.	
	ektotbelcovtot=0.
      end if
	
c--------------------------------------------------------------------------
	if (idate.eq.dmaizep(3)) then	! maize3 ! original CREHDYS : elseif
	write(*,*) ' Maize 3'
      curcropper='maize'
	!	Kcover & Kext
	! kccover=kcovermaize
	! kext=kextmaize ! cover exinction coefficient (maize)
	prstmaize=1
	dstart=dmaizep(3)
	
	! RR & KS & manning
		!	Ks & RR overland decreasing

		if (rrksvarmaize(3).eq.0) then

			rrksovar=0

		else
			rrksovar=1

	    end if
		do i=1,n
          if (sl(i).gt.(-9999)) then

			if (wheel(i).eq.0) then ! Si on n'est pas dans une wheel track
		    
			rr(i)=rrmaizeini(3)
			rrcm = rr(i)/10.	! rrm in cm

			dir2(i)= 0.243*rr(i)+0.01*rr(i)*rr(i)-0.012*rr(i)*sl(i)*100  ! previous CREHDYS: eq. from Mwandera & Feyen (1992):(0.294*rrcm+0.031*rrcm*rrcm-0.012*rrcm*sl(i)*100)*10.

			if (rrksovar.eq.0) then ! if Ksat and RR are constant during the season

				ke(i)=ksend(i)
			else
				ke(i)=kinim
			end if

			else ! Si on est dans une wheel track
			rr(i)=rrch
			rrmc=rrch/10.

			dir2(i)= 0.243*rr(i)+0.01*rr(i)*rr(i)-0.012*rr(i)*sl(i)*100  ! previous CREHDYS: eq. from Mwandera & Feyen (1992):(0.294*rrcm+0.031*rrcm*rrcm-0.012*rrcm*sl(i)*100)*10.
		    
			ke(i)=kchm2 ! 

			end if
	        end if
		end do

	! initialize sumtemp

		sumtemp=0
		sumtempdeg=0

	! initialize Ek

	ektot=0.
	ektotbel=0.
	ektotbelcov=0.
	ektottot=0.			
	ektotbeltot=0.
	ektotbelcovtot=0.
	end if	
c---------------------------------------------------------
	! end if ! original CREHDYS : end if all elseif for End of case study specific implementation (Laloy et al. (EJSS 2009b) calibration
c---------------------------------------------------------

c	End of case study specific implementation (Laloy et al. (EJSS 2009b) calibration)

c	random roughness evolution

	do i=1,n ! do loop on all grid elements (cells)	
	  rrv(i)=rr(i)*exp(-rrexp*ektotbelcovtot) ! RR of overland flow cell (not wheel tracks) calculation (Alberts et al., 1995) [mm]
	  rrcm = rrv(i)/10. ! RR from [mm] to [cm]

c	corresponding retention
	
	  if (rrksovar.eq.1) then ! If Ksat and RR are varying during the season
	    dir2(i)=0.243*rrv(i)+0.01*rrv(i)*rrv(i)-0.012*rrv(i)*sl(i)*100 ! ! Surface retention due to random roughness [mm]. In previous CREHDYS versions: eq. from Mendera & Feyen (1992): (0.294*rrcm+0.031*rrcm*rrcm-0.012*rrcm*sl(i)*100)*10.0   
          dirvon(i) = connectresh*(prop_nodep*dir2(i) + MB_stor)	
        else ! If Ksat and RR are constant during the season
	    rrv(i)=rr(i) ! random roughness = initial random roughness
          dirvon(i) = connectresh*(prop_nodep*dir2(i) + MB_stor)	
	  end if
	end do

	do m=1,n ! do loop on all grid elements (cells). (Re)initialize many variables at the start of each event (bewteen each day of the simulation)

	  qp21(m) = 0.
	  qp12(m) = 0.
	  qp22(m) = 0.
	  qi(m) = 0.
	  qi2(m)=0.
	  q2t1(m)=0.
	  q2t2(m)=0.
	  q1t1(m)=0.
	  q1t2(m)=0.
	  alp11(m)=0.
	  alp12(m)=0.
	  alp21(m)=0.
	  alp22(m)=0.
	  oldqbar(m)=0.    
	  s(m)=0. ! total ponding water of a cell at a given time in event
	  ss(m)=0.
	  hsv(m)=0.
	  hsv1221(m)=0.
	  cuminf(m)=0. ! cumulated infiltration during a rain event
	  q(m)=0
	  draintot(m)=0 ! cumulative drainage

c	erosion variables
	  qs11(m)=0.
	  qs12(m)=0.
	  qs21(m)=0.
	  qs22(m)=0.
	  qs2t1(m)=0.
	  qs2t2(m)=0.
      ! original CREHDYS : if(m.lt.90001) rcc(m)=0. ! Initialize the throughfall rate at zero everywhere, but this was an error because we are in the "m" loop while it should be in the "lline" loop
	end do
	
	wper11=0.
	wper12=0.
	wper21=0.
	wper22=0.
	lastru=0.
	lstr=0.
	lsttr=0.
	qoutlet=0.
	ccoutlet=0.
	totqsoutlets=0.
	totqoutlet=0.
	totqmmhr=0.

c	read weather data
	read(2,*) day,tempc,soitemp,radi,raites,starteventtime ! Because this READ(2,*) statement is called multiple times (in the "rainfa" subroutine within an event, and here in the daily simulation loop), each time this READ(2,*) is executed, FORTRAN automatically reads to a new line in this input, so we always stay in the current simulation day      lsttr=0. ! lsttr is time of the end of event [min]. if = 0, indicate that we are before the event
c	cover developement (percentage cover, biomass, LAI) routines

	if (tempc.lt.0) tempc=0 ! Base for degree day calculation is 0°C
	sumtemp=sumtemp+tempc ! cumulative daily temperature in [°C]
	! no cover (no intercropn nor maize) developping
	if (dstart.eq.0) then
	  prstcover=0
	  prstcovermaize=0
	  prstcoverint=0
	  prstlaiint=0
	  prstlaimaize=0
	  biomcover=0
	  lstbiomcover = 0
	end if
	! intercrop growth or decay : 
	if (dstart.gt.0.and.idate.ge.dinterp(1).and.dinterp(1).ne.0) then ! there is a cover developping or decaying (dstart>0), idate>1st day of intercrop period, 
	   if (destruction.eq.0.or.idate.lt.(dinterp(1)+destruction)) then ! before destruction or no destruction : intercrop is growing
	     sumtempdeg = 0 ! sum of temperature for degradation calculation
	     biomcover = gcoverint*(sumtemp**2) ! cover crop dry biomass [g/m²] as function of sumtemp [°C days]) - expon. curve
	     lstbiomcover = biomcover ! last biomass of cover crop
	     prstcoverint =  -1*(exp(-rbiomcover*biomcover)-1) ! conversion from dry biomass to [%] cover for cover crop
	     if (prstcoverint.gt.0.95) prstcoverint=0.95
	     if (kextint.gt.0) prstlaiint= -log(1-prstcoverint)/kextint ! Andrieu et al. 1997 & Baret et al., 1993. 
	   elseif (destruction.gt.0.and.idate.ge.(dinterp(1)+destruction)
     &           .and.
     &           (burial.eq.0
     &            .or.idate.lt.(dinterp(1)+destruction+burial))
     &          ) then ! after destruction and before burial (or no burial) : decay
	     sumtempdeg=sumtempdeg+tempc ! sum of temperature for degradation calculation
	     biomcover = exp(fN*fW*kdeg*sumtempdeg)*lstbiomcover ! cover crop dry biomass [g/m²] after destruction, following douglas & rickman cover crop decay equation
	     prstcoverint =  0.004*biomcover    ! experimental
	     if (prstcoverint.gt.0.95) prstcoverint=0.95
	     prstlaiint=0 ! No LAI after destruction (evapotranspiration and interception by residues from dead intercrop is neglected)
	   else ! after burial
	     biomcover = 0
	     lstbiomcover = 0
	     prstcoverint=0
	     prstlaiint=0
	   end if
	else ! no cover developping (dstart=0) or idate<1st day of intercrop
	   biomcover = 0
	   lstbiomcover = 0
	   prstcoverint=0
	   prstlaiint=0
	end if
    
	! maize growth : 
	if (prstmaize.gt.0) then
	   prstcovermaize=0.01*(95*1*exp(kcovermaize*sumtemp)/     ! verhulst equation for maize cover as function of sum of T° since sowing
     &              (95+1*(exp(kcovermaize*sumtemp)-1))) ! cover=f(sumtemp)
	   if (kextmaize.gt.0) then
	      prstlaimaize= -log(1-prstcovermaize)/kextmaize
	   end if
	else
	   prstcovermaize=0
	   prstlaimaize=0 
	end if
	    
	if (kextmaize.eq.0) prstlaimaize=0
	if (kextint.eq.0) prstlaiint=0
	    
	prstcover=prstcoverint+prstcovermaize  ! total present cover
	if (prstcover.gt.1) prstcover=1
	prstlai=prstlaiint+prstlaimaize ! total present LAI
	
	
	write(381,3772) idate,curcropper,tempc,biomcover,100*prstcoverint,
     &            (100*prstcovermaize),prstlaiint,prstlaimaize          ! write crop outputs
 3772 format (i5,1x,A10,1x,f6.1,1x,f10.2,1x,f10.1,1x,f10.1,1x,f10.2,
     &        1x,f10.2) ! format to write crop outputs
      
	if (raites.eq.0) goto 12	! no rainfall this day

c	raining day

c	read rainfall data

	call rainfa		! call 'rainfa' subroutine

c	compute ksat evolution

	do i=1,n

	if (wheel(i).eq.0) then ! if whe are not in a wheel track

		if (rrksovar.eq.0) then ! if Ksat and RR are not varying during the season

			keff(i)=ksend(i) ! Keff=final K


		elseif (prstcover.gt.0.01) then ! If % cover >0
	
		kcover2(i) = ke(i)*((ksend(i)/ke(i)) + (1. - (ksend(i)/ke(i)))
     .		*exp(-cfact(i)*ektotbelcovtot*(1.- (rrv(i)/40.)))) ! effective Ksat of covered soil as function of cumulative KE through canopy
	
			keff(i) = kcover2(i)			
		else ! if % cover=(~)0

		keb2(i) = ke(i)*((ksend(i)/ke(i)) + (1. - (ksend(i)/ke(i)))*
     .		exp(-cfact(i)*ektotbeltot*(1.- (rrv(i)/40.)))) ! effective Ksat of bare soil as function of cumulative free KE

				!	Risse et al. 1995.

			keff(i) = keb2(i)

		end if

	else ! If we are in a wheel track

		keff(i) =	ke(i)	
		
	end if
		
	end do
      
c     Compute percolation before event
	do i=1,n
	  if (theta(i).gt.fc(i)) then
	    call perco(theta(i),keff(i),tp(i),fc(i),df(i), ! do percolation calculation, either for the whole  day (if there was no rain event) or for the time remaining after event
     &           dt,a(i),idate,percol,lsttr,raites,starteventtime)
	    theta(i) = ((theta(i)*df(i)) - percol)/df(i)		!+ wp(i)
	    if (theta(i).lt.wp(i)) theta(i) = wp(i) ! no end if needed because statement on same line
	    draintot(i)=draintot(i)+percol
	  end if
	  thetainst(i)=theta(i) ! retrieve soil water content values after daily evapotranspiration and percolation
	end do
      


c	start interception-infiltration-runoff-rainfall routine

	do 13 i=1,lline	! Do event loop on total number of computation time step (e.g. 30s) in the rain event
      rcc(i)=0. ! initialize throughfall rate variable
	write(*,*) 'Time = ',tc(i),' min' ! original CREHDYS : write(*,*) 'Time = ',i*dtm,' min'
	
	if (prstlai.eq.0.or.rc(i).eq.0.) then ! LAI=0 or rain=0 : no interception
      rcc(i)=rain(rc(i),prstlai,i,dt,tc(i),idate,cu,prstcover) ! rcc=net rainfall rate = throughfall rate (after interception), computed as function of vegetation cover. The rain function has to be called iven if there is no vegetation for proper variables following
      rcc(i)=rc(i) ! throughfall rate = net rainfall rate
	else	
	rcc(i)=rain(rc(i),prstlai,i,dt,tc(i),idate,cu,prstcover) ! rcc=net rainfall rate = throughfall rate (after interception), computed as function of vegetation cover
	end if
	
	prstr=rcc(i) ! prstr not used?
	rtc(i)=tc(i)-dtm ! rtc=preceding computation time step = tc (=computation time in event) - dtm

	totprec=totprec + rc(i)*tcst
	pinst= rc(i)*tcst
	totprecgr=totprecgr + rcc(i)*tcst
c	kinetic energy content of the rain. Van Dijk et al. 2002. metric units

	ek = 28.3*(1-0.52*exp(-0.042*rc(i))) ! J m-2 mm-1
	ektot = ektot +	 ek*rc(i)*tcst	! J m-2

c	kinetic energy content of the rain for belgium. Bollinne et al. 1980
c				ekbel = emax*(1-a*exp(-b*rainrate)) J m-2 mm-1

	ekbel = 29*(1-0.6*exp(-0.061*rc(i))) ! J m-2 mm-1


	ektotbel = ektotbel + ekbel*rc(i)*tcst	! J m-2

	ekbelcov = 29*(1-0.6*exp(-0.061*rc(i))) ! free rainfall kinetic energy [J m-2 mm-1]
	ekbelcov2 = 29*(1-0.6*exp(-0.061*rcc(i))) ! throughfall rain kinetic energy [J m-2 mm-1]

	Tra = ekbelcov*rc(i)*tcst*(1.-prstcover) ! Tra=Maize (and intercrop residues) throughfall kinetic energy content that reachs directly the ground (J/m2)

c maize stem and canopy flow

	if (prstmaize.gt.0) then ! maize period

		if (idate.lt.(dstart+31)) then ! arbitratry height values for maize stem depending on day simulation
		height=0.5
		elseif (idate.lt.(dstart+76)) then
		height=1
		else
		height=1.5
		end if

		brandt=15.8*(height**0.5)-5.87

		Trf= rcc(i)*tcst*(prstcovermaize)*brandt*0.5  ! Trf=Kinetic energy content of drops falling off maize leaves  [J/m²] (m² of total soil, i.e. covered AND not covered)
              
	else ! not maize period -> intercrop period

		Trf=0 ! No KE falling off maize leaves because zero maize cover (KE falling from intercrop leaves or residues is neglected) 

      end if
      
      ekbelcovsplash=Trf/(prstcovermaize*rcc(i)*tcst)	! [J m-2 mm-1] (m² of covered soil)

	ektotbelcov = ektotbelcov + Tra + Trf ! ektotbelcov = KE from

	timnum2=timnum2+1
	if (timnum2.gt.1.) lstr = rcc(i)

c	compute inflitration and drainage
	
c	do 140 m2=0,maxval(flacc)

	do 14 m=1,n ! loop on all cells of the grid

	if (sl(m).gt.(-9999)) then ! if the cell is defined (has a value for the slope)
		
		if (keff(m).eq.0.) then ! If Keff=0

				filt(m)=0.		! Null infiltration
		else  ! If Keff>0
	
		filt(m)=gamptun3(hsv(m),rcc(i),keff(m),dt,rtc(i),tc(i), ! function computing (Green-Ampt) infiltration. filt [mm/h].
     .					m,n,n2,idate,cu,df(m),tp(m),theta(m),dx2,
     .					psi(m),dtm)
          
          ! if (filt(m).gt.99999) filt(m)=99999
          
		end if
	
		if (thetainst(m).gt.fc(m)) then ! if soil moisture > FC
			dr(m)=drain(thetainst(m),keff(m),tp(m),fc(m),df(m),dt, ! drain function computes drainage (percolation) at event computation time step dr[mm/h]
     .		a(m),idate)
		else ! Soil moisture <= field capacity

			dr(m) = 0. ! no drainage

		end if

	else ! the cell is not defined (has no value for slope)

	filt(m)=-9999 ! no value for infiltration
	dr(m)=-9999 ! no value for drainage

	end if

c	end if   ! flacc=m2	

   14 continue ! end of spatial do loop 14 (on all cells of the grid)
c  140 continue

c	compute runoff 150 - 15

	outlet=maxval(flacc)		! flow accumulation number of the outlet

	do 150 m2=0,maxval(flacc) ! loop on the flow accumulations

		do 15 m=1,n		!older CREHDYS versions : n2 ! spatial loop on all cells of the grid. This double loop, with the if statement at L1153, enables to perform a do loop on all cells in the flow acc order

	if (sl(m).gt.(-9999)) then ! check if valid slope value
	
	if (flacc(m).eq.m2) then   ! pixel corresponding to flow acc

	iw=int(m) ! cell indice for kinematic wave
	iw2=int(m) ! second cell indice for kinematic wave

c	iw=indice m for kinematic wave	
	
      
c	update thetainst after infiltration and drainage

	vertflux = rcc(i) + s(m)*3600/dt ! total surface vertical influx when ponding [mm/h] = throughfall rate [mm/h] + total ponding water (i.e. including previous ponding, rain excess, runoff coming from upstream)
      
	if (filt(m).gt.vertflux) then ! if infiltration rate > total surface vertical influx

		filtreal = vertflux ! real infiltration rate [mm/h]

	else ! if infiltration rate <= total surface vertical influx

		filtreal = filt(m) 

      end if
      
      fil=filtreal*cu	! filt(m) [mm/h] to fil in [m3/s]. original CREHDYS : fil=filt(m)*cu and this operation was executed BEFORE correction for vertflux so fil was sometimes very high because it is only potential infiltration rate (before correction for actual rainfall vertflux)
	
	if (fil.lt.0.) fil=0.
	
      cuminf(m)=cuminf(m) + filtreal*dt/3600 ! cumulated infiltration of cell(m) during a rainfall event [mm]
      
	thetainst(m)= thetainst(m)+(((filtreal-dr(m))/3600)*dt)/df(m) ! update soil moisture of cell m at that computation time step.
      if (m.eq.moutlet) thetaoutlet=thetainst(m)
	draintot(m) = draintot(m) + dr(m)*dt/3600  ! adjust cumulative drainage during the rainfall event [mm]
	
	if (thetainst(m).lt.wp(m)) thetainst(m)=wp(m)
	if (thetainst(m).gt.tp(m)) thetainst(m)=tp(m)
	  
c***** update boundary condition *****
			
		qintot11=0
		qintot12=0
		qsintot11=0
		qsintot12=0

		if (n.gt.1) then		!n>1
c	Initially written for a rectangular plot, works for any plot configuration
		do x=1,nrow
			do y=1,ncol

			if (xycoord(x,y).eq.m) then !  the dooble loop (L1203-1204) combined with this IF statement enables to work on the right (x,y) coordinates corresponding to current cell "m"

			if (ycoord(m).ne.1.and.ycoord(m).ne.ncol.and.
     .			xcoord(m).ne.1.and.xcoord(m).ne.nrow) then ! inner pixel
			north=xycoord(x-1,y)
			northe=xycoord(x-1,y+1)
			east=xycoord(x,y+1)
			southe=xycoord(x+1,y+1)
			south=xycoord(x+1,y)
			southw=xycoord(x+1,y-1)
			west=xycoord(x,y-1)
			northw=xycoord(x-1,y-1)
			fdn=fdir(north) ! flow direction of north cell
			fdne=fdir(northe)
			fde=fdir(east)
			fdse=fdir(southe)
			fds=fdir(south)
			fdsw=fdir(southw)
			fdw=fdir(west)
			fdnw=fdir(northw)
			
			elseif (m.eq.xycoord(1,1)) then  ! outer pixel (upper left)
			southe=xycoord(x+1,y+1)
			south=xycoord(x+1,y)
			east=xycoord(x,y+1)
			fde=fdir(east)
			fdse=fdir(southe)
			fds=fdir(south)
			fdn=0 ! no north-east to south-west cells, so no flow direction
			fdne=0
			fdsw=0
			fdw=0
			fdnw=0

			elseif (m.eq.xycoord(1,ncol)) then  ! outer pixel (upper right)
			southw=xycoord(x+1,y-1)
			south=xycoord(x+1,y)
			west=xycoord(x,y-1)
			fdsw=fdir(southw)
			fds=fdir(south)
			fdw=fdir(west)
			fdn=0
			fdne=0
			fde=0
			fdse=0
			fdnw=0

			elseif (m.eq.xycoord(nrow,1)) then  ! outer pixel (bottom left)
			north=xycoord(x-1,y)
			northe=xycoord(x-1,y+1)
			east=xycoord(x,y+1)
			fdn=fdir(north)
			fdne=fdir(northe)
			fde=fdir(east)
			fdse=0
			fds=0
			fdsw=0
			fdw=0
			fdnw=0

			elseif (m.eq.xycoord(nrow,ncol)) then  ! outer pixel (bottom right)
			north=xycoord(x-1,y)
			northw=xycoord(x-1,y-1)
			west=xycoord(x,y-1)
			fdn=fdir(north)
			fdnw=fdir(northw)
			fdw=fdir(west)
			fdsw=0
			fds=0
			fdse=0
			fde=0
			fdne=0

			elseif (ycoord(m).eq.1) then ! left
			north=xycoord(x-1,y)
			northe=xycoord(x-1,y+1)
			east=xycoord(x,y+1)
			southe=xycoord(x+1,y+1)
			south=xycoord(x+1,y)

			fdn=fdir(north)
			fdne=fdir(northe)
			fde=fdir(east)
			fdse=fdir(southe)
			fds=fdir(south)
			fdsw=0
			fdw=0
			fdnw=0
		
			elseif (ycoord(m).eq.ncol) then ! right

			north=xycoord(x-1,y)
			northw=xycoord(x-1,y-1)
			west=xycoord(x,y-1)
			southw=xycoord(x+1,y-1)
			south=xycoord(x+1,y)
			fdn=fdir(north)
			fdnw=fdir(northw)
			fdw=fdir(west)
			fdsw=fdir(southw)
			fds=fdir(south)
			fdse=0
			fde=0
			fdne=0
			
			elseif (xcoord(m).eq.1) then ! upper

			west=xycoord(x,y-1)
			southw=xycoord(x+1,y-1)
			south=xycoord(x+1,y)
			southe=xycoord(x+1,y+1)
			east=xycoord(x,y+1)
			fdn=0
			fdnw=0
			fdw=fdir(west)
			fdsw=fdir(southw)
			fds=fdir(south)
			fdse=fdir(southe)
			fde=fdir(east)
			fdne=0

			elseif (xcoord(m).eq.nrow) then ! bottom
			
			north=xycoord(x-1,y)
			northe=xycoord(x-1,y+1)
			east=xycoord(x,y+1)
			northw=xycoord(x-1,y-1)
			west=xycoord(x,y-1)
			fdn=fdir(north)
			fdnw=fdir(northw)
			fdw=fdir(west)
			fdsw=0
			fds=0
			fdse=0
			fde=fdir(east)
			fdne=fdir(northe)
			
			else

			write(*,*) ' they forget me ',m
			read(*,*)

			end if  ! innter/outer pixel
c             Compute total water and sediment inflow in current cell (=sum of outflows from surrounding cells)
			if (fdn.eq.4) then ! north cell flow direction=south (north cell supplies flow to current cell)
			qintot11=qintot11+q2t1(north) ! Total cell inflow at beginning of time step [m³/s]. q2t1(north) is the north cell outflow at beginning of time step
			qintot12=qintot12+q2t2(north) ! Total cell inflow at end of time step [m³/s]. q2t2(north) is the north cell outflow at end of time step
			qsintot11=qsintot11+qs2t1(north) ! Total cell sediment inflow at beginning of time step [m³/s]. qs2t1(north) is the north cell sediment outflow at beginning of time step
			qsintot12=qsintot12+qs2t2(north) ! Total cell sediment inflow at end of time step [m³/s]. qs2t2(north) is the north cell sediment outflow at end of time step
			q(north)=1
			end if

			if (fdne.eq.8) then ! north-east cell flow direction=south-west
			qintot11=qintot11+q2t1(northe) !
			qintot12=qintot12+q2t2(northe)
			qsintot11=qsintot11+qs2t1(northe)
			qsintot12=qsintot12+qs2t2(northe)

			q(northe)=1
			end if

			if (fde.eq.16) then ! east cell flow direction=west
			qintot11=qintot11+q2t1(east)
			qintot12=qintot12+q2t2(east)
			qsintot11=qsintot11+qs2t1(east)
			qsintot12=qsintot12+qs2t2(east)

			q(east)=1
			end if

			if (fdse.eq.32) then
			qintot11=qintot11+q2t1(southe)
			qintot12=qintot12+q2t2(southe)
			qsintot11=qsintot11+qs2t1(southe)
			qsintot12=qsintot12+qs2t2(southe)
			q(southe)=1
			end if

			if (fds.eq.64) then
			qintot11=qintot11+q2t1(south)
			qintot12=qintot12+q2t2(south)
			qsintot11=qsintot11+qs2t1(south)
			qsintot12=qsintot12+qs2t2(south)
			q(south)=1
			end if

			if (fdsw.eq.128) then
			qintot11=qintot11+q2t1(southw)
			qintot12=qintot12+q2t2(southw)
			qsintot11=qsintot11+qs2t1(southw)
			qsintot12=qsintot12+qs2t2(southw)

			q(southw)=1
			end if

			if (fdw.eq.1) then
			qintot11=qintot11+q2t1(west)
			qintot12=qintot12+q2t2(west)
			qsintot11=qsintot11+qs2t1(west)
			qsintot12=qsintot12+qs2t2(west)

			q(west)=1
			end if

			if (fdnw.eq.2) then
			qintot11=qintot11+q2t1(northw)
			qintot12=qintot12+q2t2(northw)
			qsintot11=qsintot11+qs2t1(northw)
			qsintot12=qsintot12+qs2t2(northw)

			q(northw)=1
			end if

			qp11(m)=qintot11 ! Total cell inflow at beginning of time step 	
			qp12(m)=qintot12 ! Total cell inflow at end of time step 
			qi(m)=qintot11
			qi2(m)=qintot12
			qs11(m)=qsintot11 ! Total sediment inflow at beginning of the time step [kg/s]
			qs12(m)=qsintot12 ! Total sediment inflow at end of time step [kg/s]

			end if  
			end do
		end do
		
		else 	!n=1	
		
		qp11(m)=qintot11		
		qp12(m)=qintot12	
		qi(m)=qintot11
		qi2(m)=qintot12
		qs11(m)=qsintot11
		qs12(m)=qsintot12

		end	if		
		
	qidt = (qi(m) + qi2(m))/2  ! average incoming runoff flow rate [m³/s] during time step dt (between beginning and end of time step)
	hindt=(qidt*dt/dx2)*1000 ! average incoming runoff depth of cell(m) [mm] during time step dt

	if (hsv(m).gt.0) hindt=0. ! hsv=average runoff depth of cell at end of time step. If hsv>0 (already runoff on current cell), then no incoming runoff because the ponding is already maximum

	s(m) = ss(m) + (rcc(i) - filtreal)*dt/3600 + hindt ! original CREHDYS : + (rcc(i) - filt(m))*cu*1000*dt/dx2 (unnecessarily complicated, + filt(m) while filtreal must be used). adjust total ponding water of the cell [mm] at that computation time step = ss (previous s) + rain - infiltration + incoming runoff

	if (s(m).lt.0.) s(m) = 0.

	if (i.eq.1) then ! first computation time step of the event
		oldqbar(m) = 0.
		else
		oldqbar(m) =oldqtrans
	end if
	
      if (s(m).gt.dirvon(m)) then ! ponding water > surface retention -> there is runoff
	
	dxf=dxvert
					
		goto 22
	
	end if  

c *******	no runoff ******* (s(m) <= dirvon(m))

	ss(m)=(rcc(i) - filtreal)*dt/3600 + ss(m) + hindt ! adjust previous total ponding water of the cell [mm] ! in previous CREHDYS there was not the "+ hindt"

	if (ss(m).lt.0.) ss(m)=0.
	hsv(m)=0.
	hsv1221(m)=0.
      hmm=0.
      meanhmm=0.
      tothmm=0.
      totfilt=0.
      meanfilt=0.
	q2t1(m)=0.
	q2t2(m)=0.
	qi(m)=0.
	qi2(m)=0.
	qp21(m)=0.
	qp22(m)=0.
	qs21(m)=0.
	qs22(m)=0.
	qs2t1(m)=0.
	qs2t2(m)=0.
	q2=0.
	h21=0.
	h22=0.
	h11=0.
	h12=0.
	Dfl=0.
	Transc=0.
	Dep=0.
	eros=0.
	qs2t1(m)=0.
	qs2t2(m)=0.
	goto 23

   22 continue ! there is runoff
c     ******Kinematic wave******
	
c	*** overland & channel (wheel tracks) flow cell *** 
	
	omega=0.8 ! weight for the weighted 4 point finite difference method (theta in Chow et al., 1988 or theta_w in KINEROS documentation)
	beta=0.6 ! momentum correction factor. combined with Manning equation, β = 3/5 = 0.6
	tol2=1.0e-7 ! tolerance value used for the xmin argument of Newton-Raphson iter  (different from the tolerable error for convergence !)

c	handling of diagonal flow
	if ( fdir(m).eq.128.or.fdir(m).eq.2.or.
     & fdir(m).eq.8.or.fdir(m).eq.32) then
	  dx=dxdiag ! correction of dx for diagonal flow
	else
	  dx=dxvert ! vertical flow (dx=dx)
	end if
	
	weightr=1
	qin11=oldqbar(m)/dx ! previous vertical influx [m³.m^-1.s^-1]=[m²/s]
	qin21=oldqbar(m)/dx	! previous vertical influx [m³.m^-1.s^-1]=[m²/s]
	lat=(rcc(i)*cu-fil)	! current vertical influx (rain excess) [m³/s]. Named as lateral inflow in Chow et al. and KINEROS2 manual
	qin12=(lat*weightr/dx)	! current vertical influx [m³.m^-1.s^-1]=[m²/s]
	qin22=(lat*weightr/dx) ! current vertical influx [m³.m^-1.s^-1]=[m²/s]
	
	if (qin21.lt.0.) qin21=0   
	
	if (qin11.lt.0.) qin11=0   
	
c***** alpha *********
	if (wheel(m).gt.0) then	! channel cell (wheel track)
		
	slwp=(sl(m)) ! slope (of the wetted perimeter)
		
	width2(m)=dxvert		! channel width
	
	hmin= -10*tol2

	hwp11=(((qp11(m)*mann(m)*(width2(m)**0.666667))/(sqrt(slwp)))**0.6 ! initial value for the water height in the channel cell at inflow (i) and at beginning of time step (j), calculated by Manning equation and with P=~dx
     .	)/dx
	hwp12=(((qp12(m)*mann(m)*(width2(m)**0.666667))/(sqrt(slwp)))**0.6 ! same for (i, j+1)
     .	)/dx
	hwp21=(((qp21(m)*mann(m)*(width2(m)**0.666667))/(sqrt(slwp)))**0.6 ! (i+1, j)
     .	)/dx
	

	hwp=hwp11
	qact=qp11(m) ! current flow rate
	if (qact.eq.0.) then ! if actual flow rate is null, then the height is null
		hwp11=0.
		go to 2112
	end if
	 
	call iter(errhwp, hwp, hmin, 10., ierr) ! Newton-Raphson iterations to find the water height (at (i,j)) that match the actual flow rate (calculated by preceding Newton-Raphson iterations on the continuity equation)
	hwp11=hwp

 2112	hwp=hwp12
	qact=qp12(m)
	if (qact.eq.0.) then
		hwp12=0.
		go to 2113
	end if
	call iter(errhwp,hwp,hmin,10.,ierr)
	hwp12=hwp

 2113	hwp=hwp21
	qact=qp21(m)
	if (qact.eq.0.) then
		hwp21=0.
		go to 2114
	end if
	call iter(errhwp,hwp,hmin,10.,ierr)
	hwp21=hwp

	

		if (ierr .gt. 0) then
					write(*,*) ' no convergence for wet perimeter ',M,
     .				' at idate = ',idate
					read(*,*)
          end if

	chwid=dxvert

 2114	wper11= chwid + 2*hwp11 ! wetted perimeter of channel (wheel track) at (i,j) = dx+2*h
	wper12= chwid + 2*hwp12 ! at (i, j+1)
	wper21= chwid + 2*hwp21 ! at (i+1,j)
	wper22=wper21 ! qp22 is unkwon so it's impossile to find hwp22 (even with Newton iterations). wper22 is then approximated as=wper21


	alp11(m) = (mann(m)*(wper11**0.666667)/(sqrt(slwp)))**0.6 ! alpha(i,j) for this channel cell 
	alp12(m) = (mann(m)*(wper12**0.666667)/(sqrt(slwp)))**0.6 ! alpha(i,j+1) 
	alp21(m) = (mann(m)*(wper21**0.666667)/(sqrt(slwp)))**0.6 ! alpha(i+1,j) 
	alp22(m) = (mann(m)*(wper22**0.666667)/(sqrt(slwp)))**0.6 ! alpha(i+1,j+1) 
	
	else	! overland flow cell (no channel/wheel track cell)
	
	alp= (mann(m)*(dxf**0.666667)/sqrt(sl(m)))**0.6 ! alpha for overland flow (no channel) cell. P= just dx. Same alpha for each of the four nodes (i,j),(i+1,j), etc
	
	wper11=dxf ! wetted perimeter of overland flow cell is just cell size
	wper12=dxf
	wper21=dxf
	wper22=dxf

	slwp=sl(m)
	
	width2(m)=dxvert

	alp11(m)=alp
	alp12(m)=alp
	alp21(m)=alp
	alp22(m)=alp

	
	end if

c	test for runoff

	qtest = qp11(m) + qp12(m) + qp21(m)
	                          
            if (qtest .le. 1.e-12 .and. (qin22) .le. 1.e-13) then ! very low outflow~=0

              qp22(m) = 0.
		  			
            else ! there is significant outflow
									    								                              						
				if (wheel(m).gt.0.and.timnum.eq.0.) then ! timnum is an indicator if there is flow at the outlet (>0) or not (0)
					q22=1.0e-9
				else
					q22=(qp12(m) + qp21(m))/2 ! eq. (9.6.4) from Chow et al. ; first value for q22 is linear solution as average Qvalues . original CREHDYS : q22=(qp12(m) + qp21(m))/2 (without + lat)
				end if
			
				if (q22.eq.0.) q22=1.0e-9 ! begin with initial very low but non null value for q22 because at this point it means that qp11>0 but qp12=qp21=0

				!qmin = min(-q22,-tol2)
				qmin = 0 ! lower bracket for iter subroutine. CREHDYS original : qmin=-tol2*10=-1*10^-6(=-0.000001) ! lower bracket (minimum q22) for Newton-Raphson iterations
				 							
				call iter (errq2, q22, qmin, 10., ierr) ! Newton-Raphson iterations to find numerical solution for q22 value. original CREHDYS : qmax=10 m³/s
	!		
				if (ierr .gt. 0) then ! no solution was found in iter
					write(*,*) ' no convergence for cell ',m,
     .				' at idate = ',idate
					read(*,*) ! some input from user in the terminal (e.g. enter "enter") is needed to continue
                  end if ! if no convergence, the last value for q22 is either (qp12(m) + qp21(m))/2 (L1615) or 1.0e-9 (L1618)

c	 to avoid to low outflow
				qouttest = 0.	!1.e-12
		
				if (q22 .lt. qouttest) q22 = 0.
				
				if (q22.gt.0.) then
					lastru =1 ! lastru never used in the code
				else
					lastru =0
				end if

				qp22(m) = q22 ! cell outflow [m³/s] at end of time step (Q(i+1, j+1)) computed by Newton-Raphson iterations just here above
							
	      end if
	
	q2t1(m)=qp21(m) ! cell outflow [m³/s] at beginning of time step
	q2t2(m)=qp22(m) ! cell outflow [m³/s] at end of time step (computed by Newton-Raphson iterations just here above)
	q1t1(m)=qp11(m) ! cell inflow [m³/s] at beginning of time step
	q1t2(m)=qp12(m) ! cell inflow [m³/s] at end of time step

	q2=q2t2(m)

	if (flacc(m).eq.maxval(flacc)) then ! one of outlets (and there is runoff)
	  timnum = timnum+1 ! timnum is an indicator for the number of times that kinematic wave arrived to outlet 
	  qoutlet=q2	! outflow at (last !) outlet
	end if ! outlet

	h22=(((q2t2(m)*mann(m)*(wper22**0.666667))/(sqrt(slwp)))**0.6)/ ! Manning eq. (h=A/dx) to convert outflow of cell(m) into runoff depth at out of the cell and at end of time step (h22)
     .											width2(m)
	if (m.gt.n.and.h22.ne.0.) then ! ? why ? this condition should never happen
	hwp=h22
	qact=q2t2(m)
	call iter(errhwp, hwp, hmin, 10., ierr)
	h22=hwp
	end if

	h12=(((qp12(m)*mann(m)*(wper12**0.666667))/(sqrt(slwp)))**0.6)/ ! Manning eq. to convert outflow [m³/s] of cell m at beginning of time step to runoff depth [m]
     .						width2(m)
	
	h21=(((qp21(m)*mann(m)*(wper21**0.666667))/(sqrt(slwp)))**0.6)/ ! Manning eq. to convert inflow of cell m at beginning of time step to runoff depth
     .						width2(m)

	h11=(((qp11(m)*mann(m)*(wper11**0.666667))/(sqrt(slwp)))**0.6)/ ! Manning eq. to convert inflow of cell m at beginning of time step to runoff depth
     .						width2(m)
	
	hf=(h12+h21)/2 !average runoff depth[m] of cell m between inlet at end of time step and outlet at beginning of time step. (eq. 9.6.4 from Chow et al.) 

	hf2=(h21+h22)/2 !average runoff depth [m] of cell at inlet (inflow)
	hf3=(h12+h22)/2 ! average runoff depth [m] of cell at end of time step

c	 update initial condition for next dt

	qp21(m)=qp22(m) ! adjust cell outflow at beginning of time step for next computation time step

	if (m.le.n) oldqtrans = rcc(i)*cu-fil ! oldqtrans=previous rain excess [m³/s]
	
c	***** end of Kinematic wave *****

	hsv1221(m)=hf*1000 ! average runoff depth [mm]
	hsv(m)=hf3*1000 ! hsv=average runoff depth [mm] at the end of computation time step
	ss(m)=hf3*1000+dirvon(m)	! adjust previous total ponding water [mm] of the cell
 
 3766 format(i5,1x,f7.2,1x,f7.2,1x,f10.4) ! format not used.
 3769 format(i5,1x,f7.2,1x,f7.2,1x,f7.2,1x,f10.4,1x,f10.4,1x,f10.4,1x,
     &       f10.4,1x,f15.10,1x,f15.6,1x,f7.4) ! format to write runoff outputs
 3770 format(i5,1x,f7.2,1x,f7.2,1x,f10.4,1x,f15.10,1x,f15.3,1x, ! format not used  
     .	f15.10,1x,f15.10,1x,f15.10,1x,f8.4)
 3788 format(i5,1x,i5,1x,f7.2,1x,f7.2,1x,f8.4)	! format not used
   

 8200 format(i5,1x,i5,1x,f4.1,1x,f10.8,1x,f10.8,1x,f10.8) ! format not used

!	end if   ! if pixel ~ flacc

c	erosion loop

c	!!! if no erosion
c	goto 99

	if (hsv1221(m).gt.0.01) then	! if there is runoff on that cell
			
	hmm=hsv1221(m) ! average runoff depth [mm] of cell m at that computation time step
	hf=hmm/1000 ! average runoff depth [m] of cell m

	if (rcc(i).gt.0) then ! if net rainfall (throughfall) > 0
	Dsplash = Lerosionsplash(dt,dx2,hmm, ekbel,ekbelcovsplash,rc(i),  ! splash detachment on (covered and uncovered fraction of) soil [kg/s] (eq. (2) and (3) from Laloy & Bielders 2009a)
     . rcc(i), As, prstcover, rrv(m)) !kg/s                                 ! original CREHDYS : ekbelcov2 was used instead of ekbelcovsplash !

	else ! net rainfall=0
	Dsplash =0.
	end if
	
	vel=0.5*(q2t1(m)+q1t2(m))/(hf*dxf)	! average flow velocity [m/s]

	wstream=vel*SIN(ATAN(sl(m)))*100	! unit stream power [cm/s] (eq. (6) from Laloy and Bielders, 2009a). Original CREHDYS model : wstream=vel*sl(m)*100. Rem: for slopes < 20%, this does not change a lot (sin(10%)~=10%))
	
      if ((wstream).gt.0.4) then		! unit stream power > critical unit stream power (=0.4) [cm/s]			

	Transc=ctc*((wstream-0.4)**dtc)*2650		! sediment concentration at transport capacity [kg/m³] (eq. (5) from Laloy & Bielders, 2009a)

	else ! unit stream power < critical unit stream power [cm/s]

	Transc=0

	end if

	conc=(qs21(m)+ qs12(m))/(q2t1(m) + q1t2(m)) ! average sediment concentration = average sediment flow / average water flow [kg/s / m³/s]=[kg/m³]

	if (isnan(conc)) conc=0.

	if (Transc.gt.conc) then		! if sediment concentration < transport capacity : there is flow detachment Dfl
	
	if (wheel(m).gt.0) then ! if on channel (wheel track)

	Dfl=Ydetw*(Transc-conc)*Vs*dxf*dx ! eq. (8) from Laloy & Bielders, 2009a. Note : specific cohesion for wheel tracks

	else ! if on oveland flow cell (not channel/wheel track)

	Dfl=Ydet*(Transc-conc)*Vs*dxf*dx ! same eq. BUT specific cohesion for overland flow cell

	end if

	Dflmax=(Transc-conc)*(q1t2(m)+q2t1(m))*0.5 ! maximum flow detachment (eq. 10 from Laloy & Bielders, 2009a)

	if (Dfl.gt.Dflmax) Dfl=Dflmax ! the amount of detached soil cannot exceed the remaining carrying capacity of the flow

	Dep=0 ! No deposition
	
	else    ! sediment concentration >= transport capacity : there is deposition Dep
	
	Dfl=0

	Dep=(Transc-conc)*Vs*dxf*dx ! deposition rate [kg/s]. This is a negative value


	if (abs(Dep*dt).gt.abs((Transc-conc)*dxf*dx*hf)) then ! Dep>Depmax
 
	
	Dep=((Transc-conc)*dxf*dx*hf)/dt ! Dep is maximum deposition, it is the excess of sediment compared to transport capacity (eq. 11 from Laloy & Bielders 2009a)

	end if

	end if ! end if transport capacity vs. sediment concentration
	
	eros=Dep+Dfl+Dsplash  ! total net erosion rate [kg/s] on that cell

! Lisem model:

	alpha1221=(alp12(m)+alp21(m))*0.5
	qp1221=(q1t2(m)+q2t1(m))*0.5

	ebar=eros/dx ! net erosion rate [kg/(m.s)]

! Lisem model:	.and.q2t2(m).gt.0
	if (q2t1(m).gt.0) then ! This "if" statement ensures no NaN problems because null denominators in the following qs22(m) equation
	   qs22(m)=(ebar*dx
     .  + alpha1221*(qp1221**beta)*(dx/dt)*(qs21(m)/q2t1(m)) ! Numerical solution of Qs22 according to a linear backward-difference scheme
     .  - conc*alpha1221*beta*(qp1221**(beta-1))*(dx/dt)          ! ~ eq. (9.6.7) from Chow et al. (but for sediment flow thus for mass balance equation (and not continuity equation))
     .  * (q2t2(m)-q2t1(m)) + qs12(m))
     .  / (1+alpha1221*(qp1221**beta)*(dx/dt)*(1./q2t2(m)))
	else
	  qs22(m)=0.
	end if
	concmax=Transc
	! conc max = 848 kg m-3 in Morgan et al. (1997)
	
	if (qs22(m).lt.0) qs22(m)=0.
	
	if (flacc(m).eq.maxval(flacc)) then ! outlet
	  ccoutlet = qs21(m) ! why not =qs22(m) ? ! sediment flux [kg/s] at (last !) outlet
	end if ! coutlet

c	 update initial condition for next dt

	qs21(m)=qs22(m)
	qs2t2(m)=qs22(m)
	qs2t1(m)=qs21(m)
	
	else ! average runoff depth < 0.01 (there is no runoff) on that cell -> no erosion

	Dfl=0.
	Transc=0.
	Dep=0.
      Dsplash=0.
	eros=0.
	qs21(m)=0.
	qs22(m)=0.
	qs2t1(m)=0.
	qs2t2(m)=0.

c	here you can print the water balance and erosion information associated with pixel m within the event : theta(m), qp22(m) m3/s, qs22(m) kg/s,...
	
	end if		! end erosion loop

   99 continue
   23 continue
     
	end if   ! end if pixel ~ flacc
	
	else     ! sl(m)==-9999

      end if	! end if sl(m) (slope) has valid value
      
   15 continue ! end runoff loop ; spatial loop on all cells of the grid
  150 continue ! end runoff loop ; flow acc loop on the number of cells to the outlet
	
      do m=1,n ! loop to compute total fluxes and variables on all outlets
          if (flacc(m).eq.maxval(flacc)) then ! (one of the) outlet(s) AND there is flow at at least one of the outlets 
          qmmhr=(q2t2(m)/(dx2*neff))*1000*3600
          totfilt=totfilt+filt(m) ! total infiltration
          tots=tots+s(m)
          totdirvon=totdirvon+dirvon(m)
          totqmmhr=totqmmhr+qmmhr ! total outflow [mm/h] on all outlets during this time step
          totqoutlet=totqoutlet+q2t2(m) ! total outflow [m³/s] on all outlets
          tothmm=tothmm+hsv1221(m) ! total runoff depth [mm] on all outlets
          totqsoutlets=totqsoutlets+qs2t1(m) ! total sediment flux from all outlets [kg/s]
          totrun=totrun+qmmhr*dt/3600 ! total daily runoff qmmhr*dt/3600 transforme qmmhr [mm par heure] en q en [mm pour 1 event time step (dt)]
	    toteros=toteros+(qs2t1(m)*dt) ! total daily erosion
          end if
      end do
	meanhmm=tothmm/noutlets
      meanfilt=totfilt/noutlets
      means=tots/noutlets
      meandirvon=totdirvon/noutlets
      
      
	write(379,3769) idate,rc(i), rcc(i), (timnum2)*(dt/60.), ! write event runoff outputs
     .               meanfilt, means, meandirvon,
     .           	   totqmmhr, totqoutlet, meanhmm, thetaoutlet

	write(380,3771) idate,lstr,(timnum2-1)*(dt/60.), ! write event erosion outputs
     .            totqsoutlets, Transc, totqsoutlets/totqoutlet,
     .            Dep, Dfl, Dsplash, eros
 3771 format(i5,1x,f7.2,1x,f7.2,1x,f15.10,1x,f15.10,1x,f15.3,1x,      ! format to write erosion outputs
     .	f15.10,1x,f15.10,1x,f15.10,1x,f15.10)

	if (totqmmhr.eq.0) timnum=0 ! timnum is an indicator if there is flow at the outlets
	
	qoutlet=0.
	ccoutlet=0.
	totqsoutlets=0.
	totqmmhr=0.
	totqoutlet=0.
	tothmm=0.
	meanhmm=0.
	tots=0.
	means=0.
	totfilt=0.
	meanfilt=0.
	totdirvon=0.
	meandirvon=0.

	if (i.eq.lline) lsttr = tc(i)

   13	continue ! END OF EVENT DO LOOP 13 (on all computation time steps within the event)

   12 continue ! no rainfall

c	new soil moisture content after event
      if (raites.eq.1) then ! if there was a rain event this day
	do i=1,n ! do spatial loop on all cells of the grid
	theta(i)= thetainst(i) ! retrieve soil moisture values after event calculations (drainage and infiltration)
	if (theta(i).lt.wp(i)) theta(i) = wp(i)
	if (theta(i).gt.tp(i)) theta(i) = tp(i)
      end do
      end if

c	percolation, evaporation, and transpiration routine

      do i=1,n		!do spatial loop on all cells of the grid
      ! compute percolation
      percol = 0.
      if (theta(i).gt.fc(i)) then ! there is percolation only if theta > FC
		  call perco(theta(i),keff(i),tp(i),fc(i),df(i), ! do percolation calculation, either for the whole  day (if there was no rain event) or for the time remaining after event
     .           dt,a(i),idate,percol,lsttr,raites,starteventtime)
            percol=min(percol, theta(i)*df(i) - fc(i)*df(i)) ! limit percolation until soil water content reach field capacity
            theta(i) = theta(i) - percol/df(i) ! compute theta after percolation
            draintot(i)=draintot(i)+percol ! adjust total daily percolation
      end if

      edisp =	(theta(i)-wp(i))*df(i)			! Available water in the control volume [mm]
      edispedx = (theta(i)-wp(i))*edx(i)      ! Available water in the evaporative depth [mm]

      call etp (tempc,radi,albed(i),etpmm) 
	ie = i
	call etpr2(prstlai,totprecgr,idate)
      ! compute soil evaporation
      theta(i) = theta(i) - es(i)/df(i)
      
      edisp = (theta(i)-wp(i))*df(i) ! compute new available water content after soil evaporation
      
      ! compute plant transpiration
      pep(i) = min(pep(i),edisp)
      theta(i) = theta(i) - pep(i)/df(i)
      
      end do
      
      totwloss=es(i)+pep(i)+percol
      thetaoutlet=theta(moutlet) ! record theta outlet at end of day
      
c	end of percolation, transpiration & evaporation routines

c	here print the runoff and erosion information associated with pixel of interest after the event

	write(375,3777) idate,totprec,totrun,(toteros/(REAL(neff)*dx2))*10.,thetaoutlet ! write daily runoff and erosion outputs in daily.out
 3777 format(i5,f8.4,f8.4,' mm',1x,f16.4,' t/ha',f16.4)
	write(378,3780) idate,totprec,totrun,(toteros/(REAL(neff)*dx2))*10. ! write daily runoff and erosion outputs in dailyRE.out
 3780 format(i5,f8.4,1x,f8.4,1x,f16.4)
      write(382,3773) idate,totprec,keff(moutlet),rrv(moutlet),
     .                cuminf(moutlet),draintot(moutlet),es(moutlet),
     .                pep(moutlet),totwloss,thetaoutlet ! write soil outputs
 3773 format(i5,1x,f8.2,1x,f8.2,1x,f8.2,1x,f8.2,1x,f8.2,1x,f8.2,1x,
     &       f8.2,1x,f8.2,1x,f8.4)
	totruntot = totrun+totruntot
	toterostot=10*(toteros/(neff*dx2))+toterostot
	totprectot = totprec + totprectot
	totprecgrtot = totprecgr + totprecgrtot
	ektottot=ektottot + ektot
	ektotbeltot=ektotbeltot+ektotbel
	ektotbelcovtot=ektotbelcovtot+ektotbelcov
   11 continue		! HERE, THE TEMPORAL SIMULATION STOPS. The "continue" statement, when it is writed alone, does not make anything special in Fortran. BUT here with a number (11 in this case), this statement is equivalent to a statement ending the do loop 11 ("DO 11 ...")

	write(375,3778) totprectot, totruntot, toterostot
 3778 format('total rainfall = ',f8.4,' mm '/,1x,'total runoff = ',
     .	   f8.4,' mm '/,1x,'total erosion = ',f8.4,' t/ha ')

	call cpu_time(t1)
	write(*,*) 'Time of operation was ', t1 - t0, ' seconds'
	!read(*,*)
      end


c********************Functions and subroutines********************
c---------------------------------------------------------------------------------------------
	subroutine rainfa	! Read rainfall data (within event). This subroutine is called each "day" of the simulation, but read each time the whole weather.inp file. As opposed to a function, a subroutine does not return a single value, but can return many values via its arguments and can only be used as a stand-alone command via a "CALL"
!	implicit real (a-z) ! In original Laloy CREHDYS, "implicit real a-z" was written without parentheses, but since fortran 90, parentheses are needed. This statement indicate that, in this routine, all variables beginning with any letter (from a to z) is considered as a REAL if it is not explicitely considered as an other type of numeric, like for example integers below
	save
      integer k,rline,i,j,l,jbeg
	real rcs(3000),tcs(3000)
	real beta, omega
	integer simdur, ncol, nrow, nrowout, n, n2
	common /gen/ simdur, ncol, nrow, nrowout, n, n2
	real dt, dx
	common /gen2/ dt, dx
	real tc(90000), rc(90000), dtm
	integer lline, idate
	common /rainfall/ tc, rc, dtm, lline, idate

  
	dtm=dt/60. ! dtm=fraction of event time step (e.g. 30s) in 1 minute (e.g.dtm=0.5)
	tc(1)=0 
	k=1
	l=0
	stram=0.

	tcs(k)=0.
	rcs(k)=0.
	k=2
   11 read(2,*) day,tcs(k),rcs(k),zero1,jbeg,zero2 ! # day, event time [number of minute since start of event], rainfall rate [mm/h], dummy unused, dummy (1=end of event, 0=continue), dummy unused
		! Each time the "11 READ" is executed, FORTRAN goes to next line. This way, it reads each rain event record at one, and assign the right rain value to the right tcs(k) and rcs(k). Moreover, even if we are not in the first event, FORTRAN knows how much times the "weather.inp" was already READ, so it indeed go to the next event.
	k=k+1 ! next line (next event rain record)
	if (jbeg.eq.0) goto 11 ! rain event continue : keep reading weather file -> goto "11" tag again and read next (k+1) line of event rain record
      rline=k-1 ! if jbeg=/0 : end of event : rline = number of rain records (lines) within the event 
	do i=2,rline ! loop on all records of the event
	tinc=tcs(i)-tcs(i-1) ! tinc[min]=time interval between two recods in the weather.inp file. This time interval does not have to be constant neither to be the same as event time step (e.g. 30s).
	kt=tinc/dtm ! kt=number of computation time step (e.g. 30s) between 2 rain records. e.g. : kt=1[min]/0.5[]=2
	
   	do j=1,int(kt) ! loop within a same event record, because event time step is lower than time between two rain records (within a rain event), thus the same event record has to be read several times. int(kt) round kt down to the nearest integer.
	l=l+1 ! l=indice of computation time step since beginning of the event
	lline=l ! at the end of each event reading, lline=total number of computation time step (dtm) within the rain event 
	rc(l)=rcs(i) ! read rainfall rate of this event record
	if (i.eq.2) then ! start of event (1st record of the rain event) -> have to start from minute 0
	tc(l)=tcs(2)-dtm+(j-1)*dtm ! time in the event decomposed by computation time step. e.g. : 5 min, 5.5 min, 6 min, etc
	else ! not the 1st record of the event
	tc(l)=tbeg+(j-1)*dtm ! time in the event decomposed by computation time step. e.g. : 5 min, 5.5 min, 6 min, etc
	end if
	end do ! end of loop within a same event record (on all computation time steps within a same event record)
	tbeg=tc(l)+dtm ! tbeg = last time step [minute]
	end do ! end of loop and all records of this event (this day)
	end 

c------------------------------------------------------------------------------
c------------------------------------------------------------------------------
	real function rain(rate,prstlai,i,dt,t,idate,cu,prstcover) ! Calculate net rainfall = throughfall rate (after interception)	! As opposed to a subroutine, a function must return a single value, and can be invoked from within expressions, like a write statement, inside conditions statements (if/then), etc
		! L1045 : rcc(i)=rain(rc(i),prstlai,i,dt,tc(i),idate,cu,prstcover) ! rcc=net rainfall rate
c	implicit real (a-h,o-z)
	save
      real rate,prstlai,dt,t,cu,prstcover
	integer i, idate
	real pcum, cint, lstcint,inrate
	rbmin=rate/3600.		![mm/h] to [mm/s]
	p=1-0.046*prstlai
	if (t.le.dt/60) then ! original CREHDYS : if (t.eq.0.or.t.eq.0.5)
	pcum=0.
	lstcint=0.
	cint=0.
	if (t.eq.0.) goto 23
	end if
	
	pcum=pcum+rbmin*dt ! cumulative rainfall [mm].
	
	stmax=0.935+0.498*prstlai-0.00575*prstlai*prstlai ! Maximum vegetation rainfall interception [mm]
	
	cint=stmax*(1-exp(-(1-p)*(pcum/stmax)))  ! Cumulative interception during rainfall [mm]
	inrate=((cint-lstcint)/dt)*3600
   23 rain=rate-inrate
		
	lstcint=cint
	
	return
	end
c---------------------------------------------------------------------------------------------
c---------------------------------------------------------------------------------------------
	function gamptun3(s,r,ke,dt,rt,t,m,n,n2,idate
     .					,cu,df,thets,theti,dx2,psimat,dtm) ! function computing infiltration and drainage, called each computation time step. See L1113 : filt(m)=gamptun3(hsv(m),rcc(i),keff(m),dt,rtc(i),tc(i), m,n,n2,idate,cu,df(m),tp(m),theta(m),dx2,psi(m),dtm)
c	 Chu, WRR 1978 
	save
c	implicit real (a-h,o-z)
      real s, r, ke, dt, rt, t
	integer m, n, n2, idate
	real cu,df,thets,theti,dx2,psimat,dtm
	real fmax
	real(8) fcerr,dfcerr,cinf ! real(8) enables better precision on these variables, preventing any bug of being stuck in while loop of Newton-Raphson iterations (when fcerr is impossible to reach below error tolerance eps)
	real(8) ns,lf,tex(100485),oldtex(100485),oldp(100485)
	real(8) cumf(100485), oldcumf(100485),oldt(100485),tp(100485)
	real(8) ts(100485),tf(100485),ptp(100485),pr(100485),kf
	real(8) excess(100485),rbithr,lastf(100485)
	integer pond(100485),tests(100485)
	
	ns=(thets-theti)*psimat ! ns = M*S in Chu (1978) = (difference in average soil moisture before and after wetting)*(difference in average capillary potential before and after wetting) = (total porosity (=saturated moisture content) - theta)*psi(=Green-Ampt soil matric potential at wetting front)
	date=200.
	n3=n ! n3 = total number of spatial cells
	eps=1.0e-6 ! original CREHDYS : eps=1.0e-6
	rbit=r/60. ! rbit [mm/min] = r(throughfall rate [mm/h])/60
	rbithr=r ! rbithr [mm/h]= throughfall rate
	kf=ke/60. ! kf [mm/min] = effective Ksat [mm/h]/60
	st=s ! Average runoff depth of a cell [mm] at end of time step
	if (t.eq.0.or.t.eq.dtm) then ! If one of the two first time steps of event : it is impossible to calculate variables for the preceding time step AND variables must be returned to null value between successive events
	pond(m)=0.
	excess(m)=0.
	oldp(m)=0.
	oldtex(m)=0. ! cumulative rainfall excess from preceding computation time step
	tex(m)=0.
	oldcumf(m)=0.
	oldt(m) =0.
	ptp(m)=0.
	tf(m)=0.
	ts(m)=0.
	tp(m)=0.
	cumf(m)=0.
	pr(m) =0.
	tests(m) = 0.
	lastf(m)= 0.
	fmax=r		

	if (t.eq.0.) goto 2		! This should never happen because tc(so t)=0 only when l=0, but the do loop "l" begins with l=1. When IF is written like this with the statement on the same line as the IF, no end if is needed
	endif ! end if t=0 OR t=dtm
	
	pr(m) = rbit*(dt/60.) + pr(m) ! (cumulative absolute precipitation [mm] for this rain event)

	if (rbit.gt.0.) then ! IF net rain [mm/min]> 0
	  if (rbit.gt.kf) then ! and IF rain intensity [mm/min] > effective hydraulic conductivity [mm/min]
		tests(m) = 0.
		if (pond(m).eq.0.) then ! and IF no ponding (at initial (last) time step). Indeed, even when recalled, the function gamptun3 saves the values of pond(m) from preceding call of gamptun3
	
			cun = pr(m) - oldtex(m) - (kf*ns/(rbit-kf)) ! eq. (38) from Chu (1978). Cun=surface condition indicator when there is no surface ponding at the initial time.
			if (cun.le.0.) then ! No surface ponding at current (terminal) time
				tex(m)=oldtex(m) ! no change in tex because no surface ponding at initial time so no more rain excess this time step
				cumf(m) = pr(m)-oldtex(m) ! No ponding -> cum. infiltrated water=cum. precipitation - cum. rain excess from last time step
          		!cumf(m)=oldcumf(m) + rbithr*dt/3600.
				fmax=ke*((ns/cumf(m))+1) ! eq. (11) from Chu (1978) : infiltration capacity = K(1+SM/F)
				pond(m)=0. ! Keep no ponding because cun<0
				
				goto 2

			else ! cun>0 -> surface ponding at current (terminal) time (but not at initial time because always in IF pond(m)=0)
	          
			   tp(m) = (kf*ns/(rbit-kf)-oldp(m)+oldtex(m))/rbit + ! eq. (22) in Chu (1978) : time to ponding
     .			       oldt(m)
			   ptp(m)= oldp(m) + (tp(m) - oldt(m))*rbit ! ptp=cumulative precipitation at ponding time, see eq. (13).
			   ts(m) = ((ptp(m) - oldtex(m))/ns - log(1 + (ptp(m) -
     .				   oldtex(m))/ns))*ns/kf ! ts is the pseudotime (limit constant) from Chu (1978), see eq. (15). F0 (cumulative infiltration at initial ponding time step) is calculated here as ptp-oldtex.

			   tf(m) =  t - tp(m) + ts(m) ! time on shifted time scale

			   if (tf(m).lt.0) then
					tf(m)=0
			   end if			   

			   pond(m) = 1. ! there is now ponding because cun>0
			   
			   	if (cumf(m).gt.0.) then
				cinf = cumf(m)
				else
				cinf = 0.1
				end if

				fcerr= cinf - ns*log(1 + cinf/ns) - kf*tf(m) ! eq. (28) from Chu (1978)
				dfcerr= cinf/(cinf+ns) ! error (difference) of approaching solution of eq. (28)
	
				do while (abs(fcerr).gt.eps) ! approach solution of implicit eq. (28) from Chu (1978) with Newton-Raphson iterations
				fcerr=cinf - ns*log(1 + cinf/ns) - kf*tf(m) ! eq. (28) from Chu (1978)
				dfcerr=cinf/(cinf+ns)
				cinf=cinf-(fcerr/dfcerr)
				end do
				cumf(m)=cinf
				tex(m) = pr(m) - cumf(m)
				fmax=kf*60.*(1 + ns/cumf(m)) ! eq. (11) from Chu (1978) : infiltration capacity = K*(1+S*M/F). rem : kf*60=ke [mm/h]
				if (rbit.gt.kf.and.m.eq.n3) then ! last cell
				excess(m) = excess(m) + ((rbithr-(fmax))/60.) ! excess at last cell (outlet) is not ponding because at outlet ?
     .						*(dt/60.)
				end if
				
				goto 2

			end if ! end IF cun
		
		else	! IF there was ponding at initial time (pond(m)=/0 donc pond(m)=1) (and always in case IF rain intensity [mm/min] > effective hydraulic conductivity [mm/min])

			 tf(m) =  t - tp(m) + ts(m) ! pseudotime

			   if (tf(m).lt.0) then
					tf(m)=0
			   end if
			   			 
		     if (cumf(m).gt.0.) then
			 cinf = cumf(m)
			 else
			 cinf = 0.1
			 end if

			fcerr= cinf - ns*log(1 + cinf/ns) - kf*tf(m) ! eq. (28) from Chu (1978)
			dfcerr= cinf/(cinf+ns)
		
			do while (abs(fcerr).gt.eps) ! approach solution of implicit eq. (28) from Chu (1978) with Newton-Raphson iterations
			fcerr=cinf - ns*log(1 + cinf/ns) - kf*tf(m)
			dfcerr=cinf/(cinf+ns)
			cinf=cinf-(fcerr/dfcerr)
	
			end do
			cumf(m)=cinf
			
			cp = pr(m) - cumf(m) - oldtex(m) ! eq. (45) from Chu (1978) : Cp=surface condition indicator when there is ponding at initial time.
			
			if (cp.gt.0.) then ! surface ponding at current (terminal time)
			
				pond(m) = 1. ! there is still ponding because cp>0
				tex(m) = pr(m) - cumf(m) ! cum. rain excess at current time step
				fmax=kf*60.*(1 + ns/cumf(m)) ! eq. (11) from Chu (1978) : infiltration capacity = K(1+SM/F)
				if (rbit.gt.kf.and.m.eq.n3) then ! last cell
					excess(m) = excess(m) + ((rbithr-(fmax))/60.)
     .							*(dt/60.)
				end if
			   
                 goto 2

			else ! cp<=0, no more surface ponding at current (terminal) time (but well at initial time)

				tex(m)=oldtex(m) ! no more rain excess
          		cumf(m)=pr(m)-oldtex(m)
				!cumf(m)=oldcumf(m) + rbithr*dt/3600.
				fmax=ke*((ns/cumf(m))+1) ! eq. (11) from Chu (1978) : infiltration capacity = K(1+SM/F)
				pond(m)=0. ! no more ponding becasue cp<=0
				tests(m) =0.
				goto 2
			   
		    end if ! end IF cp L2228

		end if ! end IF ponding

        else ! IF rain intensity [mm/min] <= effective hydraulic conductivity [mm/min]

		if (st.gt.1.e-7.and.pond(m).eq.1.) then ! if the soil type is defined and if there ponding at last (initial) time step
			tf(m) =  t - tp(m) + ts(m)

			   if (tf(m).lt.0) then
					tf(m)=0
			   end if			
			 
		     if (cumf(m).gt.0.) then
			 cinf = cumf(m)
			 else
			 cinf = 0.1
			 end if

     			fcerr= cinf - ns*log(1 + cinf/ns) - kf*tf(m)
			dfcerr= cinf/(cinf+ns)
		
			do while (abs(fcerr).gt.eps) ! approach solution of implicit eq. (28) from Chu (1978) with Newton-Raphson iterations
			fcerr=cinf - ns*log(1 + cinf/ns) - kf*tf(m)
			dfcerr=cinf/(cinf+ns)
			cinf=cinf-(fcerr/dfcerr)
	
			end do
			cumf(m)=cinf
			cp = pr(m) - cumf(m) - oldtex(m)
			if (cp.gt.0.) then	! still ponding
					tests(m) = 1.
					pond(m) = 1.
					tex(m) = pr(m) - cumf(m)
					fmax=kf*60.*(1 + ns/cumf(m)) ! eq. (11) from Chu (1978) : infiltration capacity = K(1+SM/F)
					if (rbit.gt.kf.and.m.eq.n3) then
					!if (m.eq.n3) then
					excess(m) = excess(m) + ((rbithr-(fmax))/60.)
     .							*(dt/60.)

					end if
				go to 2
			else ! cp<0  -> no more ponding
				tex(m)=oldtex(m)
          		cumf(m)=pr(m)-oldtex(m)
				fmax=ke*((ns/cumf(m))+1) ! eq. (11) from Chu (1978) : infiltration capacity = K(1+SM/F)
				pond(m)=0.
				tests(m) =0.
				goto 2 		 
			end if ! end if cp L2280
		else ! st=-999999 (no soil type defined) and no ponding at initial (last) time step
						
			tex(m)=oldtex(m) ! no ponding and rain<keff so no more rain excess
			cumf(m)=pr(m)-oldtex(m)
			fmax=ke*((ns/cumf(m))+1) ! eq. (11) from Chu (1978) : infiltration capacity = K(1+SM/F)
			pond(m)=0.
			tests(m) = 0.
			
			goto 2
	    end if ! end if st existing and ponding

	  end if ! end if rain intensity >< Keff

	else ! net rain=0

	  if (st.gt.1.e-7.and.pond(m).eq.1.) then
			tf(m) =  t - tp(m) + ts(m)
			 
			   if (tf(m).lt.0) then
					tf(m)=0
			   end if
			   			 
		     if (cumf(m).gt.0.) then
			 cinf = cumf(m)
			 else
			 cinf = 0.1
			 end if

			fcerr= cinf - ns*log(1 + cinf/ns) - kf*tf(m)
			dfcerr= cinf/(cinf+ns)
	
	
			do while (abs(fcerr).gt.eps) ! approach solution of implicit eq. (28) from Chu (1978) with Newton-Raphson iterations
			fcerr=cinf - ns*log(1 + cinf/ns) - kf*tf(m)
			dfcerr=cinf/(cinf+ns)
			cinf=cinf-(fcerr/dfcerr)
	
			end do
			cumf(m)=cinf
			cp = pr(m) - cumf(m) - oldtex(m)
			if (cp.gt.0.) then	! cp>0 = still ponding at current (terminal) time step	
					tests(m) = 1.
					pond(m) = 1.
					tex(m) = pr(m) - cumf(m)
					fmax=kf*60.*(1 + ns/cumf(m)) ! eq. (11) from Chu (1978) : infiltration capacity = K(1+SM/F)
					if (rbit.gt.kf.and.m.eq.n3) then

					excess(m) = excess(m) + ((rbithr-(fmax))/60.)
     .							*(dt/60.)

					end if
				go to 2
			else ! cp<0 : no more ponding at current (terminal) time step
				
				tex(m)=oldtex(m)
          		cumf(m)=pr(m)-oldtex(m)
				fmax=ke*((ns/cumf(m))+1) ! eq. (11) from Chu (1978) : infiltration capacity = K(1+SM/F)
				pond(m)=0.
				tests(m) =0.
				goto 2 		 
			end if
	  else ! (st unexisting) or there is no ponding at last (initial) time step
			
			
			tex(m)=oldtex(m)
			cumf(m)=pr(m)-oldtex(m)
			fmax=ke*((ns/cumf(m))+1) ! eq. (11) from Chu (1978) : infiltration capacity = K(1+SM/F)
			pond(m)=0.
			tests(m) = 0.
			
			goto 2
	    
	  end if ! end if st or ponding
	end if ! end if net rain
	

    2 gamptun3=fmax ! final infiltration rate [mm/h]
      
      oldtex(m)=tex(m) ! End of gamptun3 function so adjust old values (of last time step) for the next call (next time step) of gamptun3 function
	oldcumf(m)=cumf(m)
	oldp(m)=pr(m)
	oldt(m)=t
      lastf(m)=fmax ! last infiltration rate
	
   12 format (i5,1x,i5,1x,f7.2,1x,f8.4,1x,f8.4,1x,f8.4,1x,f8.4,1x,i2 ! format not used 
     .,1x,i2,1x,f8.4,1x,f8.4,1x)
	
	return

	end
	
c------------------------------------------------------------------------------------------------	
c------------------------------------------------------------------------------------------------	       
	real function drain(tetinst,ks,thets,fc,df,dt,ks2,idate) ! percolation calculation within event (time step dt). For between events, see "perco" subroutine
														! L1120 : dr(m)=drain(thetainst(m),keff(m),tp(m),fc(m),df(m),dt,a(m), idatee)
      save
c	implicit real (a-z)
	real tetinst,ks,thets,fc,df,dt,ks2
	integer idate

	dthr=dt/3600 ! ! transform event time step (e.g. 30s) to hour time step. dthr=fraction of event time step (dt) in 1h
 
 10   if (tetinst.gt.fc) goto 20 ! if soil moisture > field capacity. No "end if" needed because statement on the same line as IF
      dr=0. ! else soil moisture<=field capacity
	drain=dr
      return
 
********************************
 20   continue
	rkfc=ks2*(tetinst/thets)**(-2.655/log10(fc/thets)) ! rkfc=adjusted K (hydraulic conductivity) for subsurface [mm/h] because ks2 (a(m) so a(i) in the main code body) is in [mm/h]
	ti=((tetinst-fc)*df)/rkfc ! ti=travel time trough control volume [h]. ? Why *thets in original CREHDYS ?

	if (ti.lt.(dt/3600.)) ti=dt/3600. ! travel time cannot be lower than event time step

	dr=df*(tetinst-fc)*(1-exp(-(dthr/ti))) !dr= [mm/dthr] ! Original Laloy comment : "dr in mm/h", but this is false because see Savabi & Williams. dt=dthr
	
	drain=dr/dthr ! output drain [mm/h]
	return
      end

c------------------------------------------------------------------------------------------------	
c------------------------------------------------------------------------------------------------	
C   Code type: FORTRAN subroutine ITER

C   Compiler: Fortran77

C   Taken from the KINEROS model

C   Date: 10/93

C   Description:

C     Newton-Raphson iterative procedure to locate the root of
C     "errfct", where "errfct" is an equation written in residual
C     form, i.e. f(x) = 0. Successive estimates are required to
C     bracket the root; an estimate which falls outside the previous
C     interval will be recomputed based on a bisection of the
C     interval. The routine will return when the residual or
C     difference between successive estimates is smaller than the
C     tolerance. These test criteria are scaled by the magnitude of
C     the estimate when it is greater than one, so that the
C     relative accuracy is more or less independent of scale.
C------------------------------------------------------------------------------
C  Arguments:

C    errfct     --           external subroutine which computes the residual
C                            and derivative of the error function at an
C                            assumed root,

C         x     real         when called, contains an initial estimate of the
C                            root, on exit, the final estimate

C      xmin     real         lower bound on x,

C      xmax     real         upper bound for the initial interval,

C      ierr     integer      0 = no error,
C                            1 = did not converge.
C
C      trace    character    location of source of call
C------------------------------------------------------------------------------
      subroutine iter (errfct, xvar, xmin, xmax, ierr)
	
      integer i, ierr
	real xvar,xmin,xmax,xlow,xhigh,x1,x2,
     .fx,dfx,ofx,deriv,xdenom,test1,test2
      data tol, rtol /1.E-7,1.e-4/ ! tolerance value for convergence criterion (error tolerated on the final solution of the residual equation). Original CREHDYS : data tol, rtol /1.E-7,1.e-4/

      xlow = xmin
      xhigh = xmax
      ierr = 1 ! by default this flag means that solution did not converge
      call errfct (xlow, fxlow, deriv) ! compute fxlow which is the solution of the residual equation for xlow

      call errfct (xhigh, fxhigh, deriv) ! compute fxhigh which is the solution of the residual equation for xmin

      if(xlow. gt. xhigh .and. fxlow .gt. fxhigh) then
        write(*, 977) xlow, xhigh
	  read(*,*)
        stop ' Limits contradiction in ITER'
      end if
  977  format(' Lower bound higher than upper bound in call to ITER'/
     &       '  xlow = ',g13.4,',  xhigh = ',g13.4)
      x1 = xvar ! initial value for x

      call errfct (x1, fx, deriv) ! first call of the function with previous values as initial value. fx is the first solution of the residual equation, and deriv is its derivative.
C                                                                iteration loop
C------------------------------------------------------------------------------

      do i = 1, 200 ! number of iterations : in KINER0S2: 50; in CREHDYS, 200

        if (deriv .ne. 0. .and. i .le. 100) then ! maximal number of iterations before trying with bissection : in KINEROS2, 10. In CREHDYS, 100. 
C                                                       compute newton estimate
          x2 = x1 - fx / deriv ! next x value

          if (x2 .le. xlow .or. x2 .ge. xhigh) then ! new estimate falls outside the expected interval
C                                                                        bisect
            x2 = (xlow + xhigh) / 2. ! next x value is the center of the expected interval

          end if

        else ! deriv=0 (if deriv=0, the Newton calculation of next x value is not applicable because denominator is zero) OR i(number of Newton iterations)>100
C                                                                        bisect
          x2 = (xlow + xhigh) / 2.

        end if
C                                            get function value of new estimate
        ofx = fx ! ofx=old fx
  10    call errfct (x2, fx, deriv) ! call again the function and derivative, with new x (x2). Original CREHDYs : at L2554
C
  900 format(' itr:', i2,3g12.4)

        test1 = abs(fx)
C                                             ! original comment in KINEROS2 : "scale the test criteria if x2 < 1". Error, they meant if x2>1 ?
        if (abs(x2) .gt. 100485.*tol)   test1 = abs (fx / x2) ! originally in KINEROS2 : if (abs(x2).gt.1000.*tol)
C!!  modify denominator in case x2 is zero:  12/02
        xdenom = max(abs(x1),abs(x2))
        if(xdenom .le. 1.e-9) then
C          write(77,'(a6,i2,6g12.3)') trace,i,x,x1,x2,xmin,xmax
C          stop ' error stop in iter '
          test2 = 2.*rtol
        else
          test2 = abs ((x2 - x1) / xdenom)
        end if
        dfx = abs(ofx - fx) 
C                                                              convergence test
        if (test1 .lt. tol .or. test2 .lt. rtol .or. (dfx .lt. tol
     &  .and. abs(x2) .lt. tol)) then

          xvar = x2
C          ierr = 0
          ierr = -i ! ierr is negative number indicating that convergence was found and at which iteration it was found.
          return ! end subroutine iter

        end if
C                                                              bracket the root
        if (fx * fxhigh .gt. 0.) then ! at this point of the code subroutine, it means that convergence was not found because the "return" was not executed

          xhigh = x2 

        else ! fx*fxhigh<=0

          xlow = x2

        end if

        x1 = x2 ! adjust x1 preceding x value

      end do
C                                                            end iteration loop
C------------------------------------------------------------------------------
      return

      end

C-------------------------------------------------------------------------------
C-------------------------------------------------------------------------------
	subroutine iter2(funcd,rtsafe,xmin,xmax,ierr)
	INTEGER MAXIT,ierr
	REAL rtsafe,x1,x2,xacc,xmin,xmax,xvar
	EXTERNAL funcd
	PARAMETER (MAXIT=100) !Maximum allowed number of iterations.

c	Using a combination of Newton-Raphson and bisection, find the root of a function bracketed 
c	between x1 and x2. The root, returned as the function value rtsafe, will be refined until 
c	its accuracy is known within �xacc. funcd is a user-supplied subroutine which returns 
c	both the function value and the first derivative of the function. 

	INTEGER j
	!REAL df,dx,dxold,f,fh,fl,temp,xh,xl
	real df,dx,dxold,f,fh,fl,temp,xh,xl
	ierr=0
	x1=xmin
	x2=xmax
	xacc=1.0e-7

	call funcd(x1,fl,df)
	call funcd(x2,fh,df)
	!if(x1. gt. x2 .and. fx1 .gt. fx2) then
	if((fl.gt.0..and.fh.gt.0.).or.(fl.lt.0..and.fh.lt.0.)) then
		write(*,*) 'root must be bracketed in rtsafe'
		read(*,*)
	    write(399,*) rtsafe,x1,x2,fl,fh
		stop 'root must be bracketed in rtsafe'
	end if
	
	if(fl.eq.0.)then 
	rtsafe=x1 
	return 
	else if(fh.eq.0.)then 
	rtsafe=x2 
	return 
	else if(fl.lt.0.)then	!Orient the search so that f(xl) < 0. 
	xl=x1 
	xh=x2 
	else 
	xh=x1 
	xl=x2 
	endif 
!	rtsafe=.5*(x1+x2)		!Initialize the guess for root, 
!	rtsafe=xvar
	
	dxold=abs(x2-x1)		!the 'stepsize before last', 
	dx=dxold				!and the last step. 
	call funcd(rtsafe,f,df) 
	do 11 j=1,MAXIT			!Loop over allowed iterations. 

	if(((rtsafe-xh)*df-f)*((rtsafe-xl)*df-f).gt.0. !Bisect if Newton out of range, 

     . .or. abs(2.*f).gt.abs(dxold*df) ) then		   !or not decreasing fast enough. 

	dxold=dx 
	dx=0.5*(xh-xl) 
	rtsafe=xl+dx 
	if(xl.eq.rtsafe)return 

	else 
	dxold=dx 
	dx=f/df 
	temp=rtsafe 

c	Change in root is negligible. 
c	Newton step acceptable. Take it. 

	rtsafe=rtsafe-dx 

	if(temp.eq.rtsafe)return 
	endif 
	if(abs(dx).lt.xacc) return	!Convergence criterion. 
	call funcd(rtsafe,f,df)		!The one new function evaluation per iteration. 
	if(f.lt.0.) then			!Maintain the bracket on the root. 

	xl=rtsafe 

	else 

	xh=rtsafe 
	
	endif 
   11	enddo  
	!pause 'rtsafe exceeding maximum iterations' 
	ierr=1 ! did not converge
	return 
	END 

c-----------------------------------------------------------------------------------------
c-----------------------------------------------------------------------------------------
	subroutine errq2 (q22, ferr, dferr)
	
C   Computes the residual ferr and derivative dferr at q22 of the error
C   function form of the finite difference approximation of the kinematic
C   equation for overland and channel (wheel tracks) flow.

	real qq
	real qin21,qin22,qin12,qin11,beta,omega
	common /kin/ qin21,qin22,qin12,qin11,beta,omega
	real qp21(100485), qp12(100485), qp11(100485)
	integer iw
	common /kin2/ qp21, qp12, qp11, iw
	real alp11(100485),alp12(100485),alp21(100485),alp22(100485)
	common /kin3/ alp11,alp12,alp21,alp22
	real dt, dx
	common /gen2/ dt, dx


	ferr = alp22(iw)*(q22**beta)+alp12(iw)*(qp12(iw)**beta)-alp21(iw)* ! ~ eq. (29) from KINEROS documentation
     .	   (qp21(iw)**beta)-alp11(iw)*(qp11(iw)**beta)+(2*dt/dx)*
     .	   (omega*(q22-qp12(iw))+(1-omega)*(qp21(iw)-qp11(iw)))-dt*
     .	   (qin12+qin22)

	if (q22.ne.0.) then
	dferr = alp22(iw)*beta*q22**(beta-1) + (2*dt/dx)*omega ! first derivative of ferr as q22 function
	else
	dferr = (2*dt/dx)*omega
	end if
	
      end
c-------------------------------------------------------------------------------------------------
c-------------------------------------------------------------------------------------------------
	subroutine errq3 (q22, ferr, dferr)
	
C   Computes the residual ferr and derivative dferr at q22 of the error
C   function form of the finite difference approximation of the kinematic
C   equation for overland and channel (wheel tracks) flow.

	real qq,rterm
	real qin21,qin22,qin12,alp,beta,omega !,quin11
 	common /kin/ qin21,qin22,qin12,alp,beta,omega   ! alp == qin11 ??? TODO check
	real qp21(100485), qp12(100485), qp11(100485)
	integer iw
	common /kin2/ qp21, qp12, qp11, iw
	real dt, dx
	common /gen2/ dt, dx

      
	ferr = alp*(q22**beta)+alp*(qp12(iw)**beta)-alp*(qp21(iw)**beta)-
     .	   alp*(qp11(iw)**beta)+(2*dt/dx)*(omega*(q22-qp12(iw))+
     .	   (1-omega)*(qp21(iw)-qp11(iw)))-dt*(qin21+qin22)

	dferr = alp*beta*q22**(beta-1) + (2*dt/dx)*omega
	
      end
c-------------------------------------------------------------------------------------------------
c-------------------------------------------------------------------------------------------------
	subroutine errhwp (hwp, ferr, dferr)
	real qj
	real slwp, qact, width2(100485)
	integer iw2
	common /hwpr/ slwp, qact, width2, iw2
	real mann(100485)
	common /manning/ mann
 	
	qj= ((slwp**0.5)/mann(iw2))*(((width2(iw2)*hwp)**1.666667)/ ! compute flow rate from manning equation with A=dx*h and P=dx+2h
     .	((width2(iw2)+2*hwp)**0.666667))

	ferr = 1 - (qact/qj) ! residual form (f=0) of a function for which flow rate computed from hwp by Manning equation is equal to the actual flow rate

	dferr = (5*width2(iw2) + 6*hwp)/(3*hwp*(width2(iw2)+2*hwp)) ! derivative of ferr. see eq. (5.6.15) (and table 5.6.1 for rectangle channel) from Chow et al. 
	
      end

c-----------------------------------------------------------------------------------------------
c-------------------------------------------------------------------------------------------------
	subroutine etp (tempc,radi,albed,etpmm) 
c	implicit real (a-z)  
	real tempc, radi, albed, etpmm
      
       tempk = tempc + 273 ! Temp [°C] to [Kelvin]
       tempr = 5304./tempk
 
c	etp Ritchie's equation 
 
       delta = (tempr/tempk)*exp(21.25-tempr) ! Slope of saturation vapor pressure curvepercol !  eq. (7.24) from WEPP documentation or eq. (14) in Bouraoui et Dillaha
	
	 h0 = (1-albed)*radi*23.9	 ! Net radiation (not reflected by albedo) in [Ly] (Langley) . *23.9 is a conversion factor from [MJ/m²(/day)] to Langley [Ly(/day)]  ! original CREHDYS : h0 = (1-albed)*radi*1./548.11

       etpm = (0.00128/58.3)*h0*delta/(0.68+delta) ! potential evapotranspiration E0 [m/day] eq. from WEPP documentation
       etpmm = etpm*1000.		! in [mm(/day)] 
 
      return
      end

c--------------------------------------------------------------------------------------------
c-------------------------------------------------------------------------------------------------	
	subroutine perco(theti,ks,thets,fc,df,dt,ks2,idate,percol,
     .             lsttr, raites, starteventtime) ! this function calculates percolation bewteen events (time step = daily or ~ hours between rain events). the "drain" function is the equivalent but within events
			! this funciton is called like this: perco(theta(i),keff(i),tp(i),fc(i),df(i),dt,a(i),idate,percol,lsttr,raites, starteventtime)
c	implicit real (a-z)
	real theti,ks,thets,fc,df,dt,ks2,percol
	real lsttr, starteventtime
	integer idate, raites

      if (raites.eq.0) then ! if there was no rain event this day : compute percolation for whole day (24h)
	dthr=24. ! delta time step [h] : should be=24h because this function computes the percolation between two consecutive days. original CREHDYS : dthr=dt/3600.
      else ! if there was a rain event this day
          if (lsttr.eq.0) then ! if we are before the event
          dthr=starteventtime ! time before start of event
          else ! we are after the event
          dthr=24.-starteventtime-(lsttr/60) ! time remaining after end of event
          end if
      end if
      
	if (theti.ge.thets) theti=thets ! thets=total porosity. theti (soil moisture content) cannot be > total porosity. No end if is needed because statement written on the same line
 
********************************

	rkfc=ks2*((theti)/thets)**(-2.655/log10(fc/thets)) ! rkfc=hydraulic conductivity of subsurface [mm/h]
	ti=(((theti)-fc)*df)/rkfc ! ti=travel time through control volume [h]. original CREHDYS : ti=(((theti)-fc)*thets*df)/rkfc

	percol=df*((theti)-fc)*(1-exp(-(dthr/ti))) ! [mm(/dthr)]. original CREHDYS : [mm/h]
	
	return
	end
	
c------------------------------------------------------------------------------------------
c-------------------------------------------------------------------------------------------------
	subroutine etpr2 (prstlai,totprec,idate)
c	implicit real (a-h,o-z)
	real prstlai,totprec
	integer idate
	real cump2(100485)
	integer i,ttime(100485)
	real etpmm
	common /etpot/ etpmm
	real esu(100485),edx(100485),es1p(100485),es2p(100485),
     &     es(100485),pep(100485),cuminf(100485),
     &     alphas(100485),edisp,edispedx
	integer ie
	common /etpact/ esu,edx,es1p,es2p,
     &     es,pep,cuminf,
     &     alphas,edisp,edispedx,ie

	i=ie

	esm = etpmm*exp(-0.4*prstlai)   ! Potential soil evaporation [mm]
      esm=min(esm,edisp,edispedx)     ! Limit potential soil evaporation to available water in the evaporative depth (when edx<df) OR to total available water in df (when edx>df)
c	soil evaporation

	if (idate.eq.1) then
	
	  es1p(i)  = 0.   ! es1p mm accumulated since first stage beginning
	  es2p(i)  = 0.   ! es2p mm accumulated since second stage beginning
	  ttime(i) = 0    ! time since stage 2 beginning
	  cump2(i) = 0.   ! cumulative precipitation since stage 2 beginning
	
	end if
      
	if (es2p(i).gt.0.) cump2(i)=cump2(i)+totprec ! if we are in stage 2 : compute cumulative precipitation since stage 2 beginning
      
	if (cump2(i).gt.es2p(i)) then ! new stage 1
	
	  es1p(i)  = 0.
	  es2p(i)  = 0.
	  ttime(i) = 0
	  cump2(i) = 0.
          
	  es1p(i) = es1p(i) + esm	

	  goto 31

	else ! cump2 < es2p : either stage 1 or stage 2 continue

		if (es2p(i).eq.0.) then ! stage 1 continue

			es1p(i) = es1p(i) + esm

			goto 31

		else ! es2p>0 : we were in stage 2 and stage 2 continue

			goto 32

		end if

	end if

   31	if (es1p(i).gt.esu(i)) then		! when accumulated soil evaporation > upper limit of soil evaporation -> pass from stage 1 (moisture not limiting) to stage 2 (moisture limiting)

		goto 32 ! pass to stage 2

	else		! es1p<=esu : we stay in stage 1	

		es(i) = max(0.,esm)		! soil evaporation = potential soil evaporation (except if <0 (Na), then it is set to 0) ! original CREHDYS : max1 and not max
		es2p(i) = 0.
		ttime(i) = 0

		goto 33

	end if


   32 if (es2p(i).ge.0.) then ! stage 2 (moisture limiting)

	ttime(i)=ttime(i) + 1		! time since stage 2 started [days]

	es2=alphas(i)*(ttime(i)**(0.5))-alphas(i)*((ttime(i)-1.)**(0.5))		! potential soil evaporation at stage 2 (moisture limiting)

	es(i) = min(esm,es2)		! Soil evaporation = minimum between potential soil evaporation when moisture is not limiting and potential soil evaporation when moisture is limiting. The latter case is the more likely. BUT that way, if LAI or Radiation or transmissivity are so low that it is impossible to reach the empirical relation of stage 2, then the equation from stage 1 is kept.
		
	es2p(i) = es2p(i) + es(i) 

	goto 33
      end if
c	plant transpiration

   33	if (prstlai.le.3.) then

		pep(i) = (etpmm)*prstlai/3.		! plant evaporation (transpiration) [mm]. original CREHDYS : pep(i) = (etpmm-es(i))*prstlai/3.
	else

		pep(i) = (etpmm-es(i))

	end if
	
	return 
	end

c------------------------------------------------------------------------------------------
c------------------------------------------------------------------------------------------
	real function Lerosionsplash(dt,dx2,hsvmm, ek,ekt,p,pt, ! L1718 : Dsplash=Lerosionsplash(dt,dx2,hmm, ekbel,ekbelcovsplash,rc(i),rcc(i),As,prstcover, rrv(m))
     .                         As, prstcover, RR)           ! original CREHDYS : ekbelcov2 was unsed instead of ekbelcovsplash
c	implicit real (a-z)
	real dt, dx2, hsvmm, ek, ekt, p, pt, As, prstcover, RR
	
c	Ds = splash detachment gs-1 occurs only if ponded depth h > surface storage
c	Ds_r = rainfall (uncovered fraction), Ds_t = troughfall (vegetation-covered fraction)	
      
	afpa=1.406*RR**(-0.942) ! 'a' coefficient for calculation of the fraction of ponded area. See Jetten (2002)
	fpa=1-exp(-afpa*hsvmm) ! fraction of ponded area []
      
	Ar_ponded=dx2*(1-prstcover)*fpa ! uncovered AND ponded area [m²]
	Ar_dry=dx2*(1-prstcover)*(1-fpa) ! uncovered AND dry area [m²]
	At_ponded=dx2*(prstcover)*fpa ! vegetation-covered AND ponded area [m²]
	At_dry=dx2*(prstcover)*(1-fpa) ! vegetation-covered AND dry area [m²]
      
	Ds_r_ponded=1e-3*((2.82/As)*ek*exp(-1.48*hsvmm)+2.96)*
     &            p*Ar_ponded/3600	![kg/s] (eq.(2) from Laloy & Bielders, 2009a) ! ponded
	Ds_r_dry=1e-3*((2.82/As)*ek*exp(-1.48*0)+2.96)*p*Ar_dry/3600	! dry -> hsvmm (water depth) = 0 mm

c	Ds_r=1e-3*((2.82/As)*ek*dexp(-1.48*hsvmm)+2.96)*(p*dt/3600)*Ar/dt	!kg/s

	Ds_t_ponded=1e-3*((2.82/As)*ekt*exp(-1.48*hsvmm)+2.96)*
     &            pt*At_ponded/3600	!kg/s
	Ds_t_dry=1e-3*((2.82/As)*ekt*exp(-1.48*0)+2.96)*pt*At_dry/3600	!kg/s
      
	Ds=Ds_r_ponded + Ds_r_dry + Ds_t_ponded + Ds_t_dry
	Lerosionsplash=Ds
	return 
	end
c------------------------------------------------------------------------------------------
