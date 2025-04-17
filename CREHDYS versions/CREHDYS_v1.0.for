c	Continuous Runoff and Erosion model CREHDYS
c	Beta version 2007, the model does the job but:
c	 - some parts of the code are redundant and can be optimized
c	 - some variables remain for historical reasons and are not usefull
c	Developed by Eric Laloy and Charles Bielders
c	For more information, please see:
c	Laloy and Bielders (JoH 2008, EJSS 2009a, EJSS 2009b)
c	Laloy et al. (JoH 2010)

	implicit none

c	Block 1 - General variables, containing information about the pixel size, type and time step during events
c
	real dx, dxvert, dxdiag, dx2, cu, tcst,dt,t0,t1
	real dxf
	real dstart
	integer simdur, idate,lline
	common /gen2/ dt,dx
	integer m2,laicheck,lch,ierr


c	Block 2 - Hill-slope characteristics, including flow routing informations

	integer flacc(100485), neff
	integer outlet, xycoord(400,400)
	integer nsoil,i,j,ncrop,n,m
     .	,k,ch_check,n2, x,y,
     .	count,ncol, nrow, nrowout
	integer north, northw,west,southw,south,southe,east,
     .		northe, fdn,fdne,fde,fdse,fds,fdsw,fdw,fdnw
	real sl(100485), fdir(100485),xcoord(100485),
     .				ycoord(100485),
     .				spix(100485),cpix(100485),
     .				numpix(100485)
	common /gen/ simdur,ncol,nrow,nrowout,n,n2


c	Block 3 - Variables for rainfall, rainfall kinetic energy and climat

	integer tempc,soitemp,radi,raites,l
	
	real rcc(90000), rtc(90000), dtm,tc,
     .	 rc,rain,pinst

	common /rainfall/ tc(90000),rc(90000),dtm,lline,
     .	idate

	real TOTPRECTOT,TOTPRECGRTOT,EKTOT,EKTOTBEL,EKTOTBELCOV,EKTOTTOT
     .	,EKTOTBELTOT, EKTOTBELCOVTOT,TOTPREC,TOTPRECGR,EKBELSPLASH, 
     .	EKBELCOV, EKBELCOV2,EKBEL,DAY,EK


c	Block 4 - Infiltration and drainage processes

	real tp(100485),fc(100485),wp(100485),a(100485),df(100485)
     .	,ksat(100485),dr(100485)

	real FIL,VERTFLUX,FILTREAL

	real thetainst(100485),theta(100485),psi(100485),
     .	ke(100485),pond(100485), asmlim(100485)
     .	,cuminf,LSTTR
	
	integer rrksovar

	real filt(100485),draintot(100485)
	real gamptun3,drain,keb2(100485),
     .	kcover2(100485),cfact(100485),
     .	ksend(100485),kinio,kinim,keff(100485)

	real kso2, kcho2, kchm2, hfront, c_kso


c	Block 5 - Evapotranspiration process

	real esu,edx,es1p,es2p,
     .	 es,pep,alphas,edispm,
     .	 edispmfc,albed(100485)

	integer ie

	real cl(100485),sa(100485),st(100485)
     .	,wcf(100485), vfs(100485),ETPMM

	real THETBEFETR,EDISP, TETPR,PERCOL,TOTWLOSS

	common /etpot/ etpmm

	common / etpact / esu(100485),edx(100485),es1p(100485),es2p(100485),
     .				  es(100485),pep(100485),cuminf(100485)
     .					,alphas(100485),edispm,edispmfc,ie


c	Block 6 - Random roughness and surface storage

	real RRCM, RRMC

	real dir2(100485),rrv(100485),
     .	dirvon(100485),rr(100485),
     .	rrch,rrexp,c_rr

	real rrintini(10),rrmaizeini(10)

     				
c	Block 7 - Runoff height

	real HF3,HMM,HINDT
	
	real hsv(100485),ss(100485),hsv1221(100485),s(100485)


c	Block 8 - Runoff dynamics and balance

      real alpha1221, qp1221,TOTRUNTOT,TOTRUN,QINTOT11,QINTOT12,
     .	QIDT,lat,QOUTTEST,QMMHR,chwid,mnch,mn,QOUTLET,WEIGHTR

	real q(100485),wheel(100485)

	real mann,omega,slwp,qact,oldqbar(100485),oldqtrans,man2
	real beta,qmin,qi2(100485),lstr,prstr,qi(100485),q21
	real q2t1(100485),q2t2(100485),q1t1(100485),q1t2(100485)
	integer iw,t1step,iw2,lastru
	integer timnum,timnum2
	real qp12d,qp21d,width2
	real q2,tol2,h11,h21,h12,h22,alp,testck,ck
	real qp21,qp12,qtest,hf,hf2,q22,qp11,hwp11,hwp12,hwp21
	real qp22(100485),wper11,wper12,wper21,wper22,hwp
	real alp11,alp12,alp21,alp22,hmin
	real qin21,qin22,qin12,qin11
	common /kin1/ qp12d,qp21d
	common /kin3/ alp11(100485),alp12(100485),alp21(100485),alp22(100485)
	common /kin/ qin21,qin22,qin12,qin11,beta,omega
	common /kin2/ qp21(100485),qp12(100485),qp11(100485),iw
	common /hwpr/ slwp,qact,width2(100485),iw2
      common /manning/ mann(100485)


c	Block 9 - Crop rotation and cover, LAI and biomass development

	integer prstmaize

	real cover(100485)

	real kextint, kextmaize, gcoverint, kcovermaize, biomcover,
     .	rbiomcover,KCCOVER,KEXT,prstlai,lstbiomcover,FN,FW,KDEG

	real sumtemp, sumtempdeg, checkd,PRSTCOVER
	
	integer dinterp(10), dmaizep(10),wheelinter(10),wheelmaize(10)
     .	,ninterp, nmaizep, rrksvarint(10), rrksvarmaize(10), 
     .	destruction, burial

 
c	Block 10 - Maize canpoy and stem flow

	real TRA, HEIGHT, BRANDT, TRF


c	Block 11 - Erosion processes

	real qs22(100485), qs21(100485),qs12(100485),qs11(100485),
     .	qs2t1(100485), qs2t2(100485)

      real As,rhop,rhow,dclay,dsilt,DSAND,DVFS,DCF,CLAYFRACT,
     .	SILTFRACT,SANDFRACT,VFSFRACT,CFFRACT,NUDYNWAT,DIAM,
     .	VS,d50,CTC,DTC,COH,YDET,TOTEROSTOT, TOTEROS, CCOUTLET,
     .    QSINTOT11,QSINTOT12, DFL,TRANSC,DEP,EROS, 
     .	DSPLASH,LEROSIONSPLASH,VEL,WSTREAM,CONC,DFLMAX,EBAR,
     .    COHW,YDETW, concmax

	
c	Block 12 - Routines for kinematic wave equation:
c	fully implicit four points finite difference approximation using Newton-Raphson iteration
c	internal: Iter and Iter2
c	external:
	external errq2,errq3,errhwp


c	START

	call cpu_time(t0)

c	open files

	open (unit=1, file='gen.inp', status='old')
	open (unit=2, file='weather.inp', status='old')
	open (unit=3, file='plot.inp', status='old')
c	open (unit=5, file='cover.out', status='unknown')
c	open (unit=399, file='filt.out', status='unknown')
	open (unit=375,file='daily.out',status='unknown')
c	open (unit=376,file='events.out',status='unknown')
	open (unit=378,file='dailyRE.out',status='unknown')
	open (unit=379,file='eventR.out',status='unknown')
	open (unit=380,file='eventE.out',status='unknown')
	open (unit=9,file='parameros.inp',status='unknown')
	open (unit=10,file='paramhyd.inp',status='unknown')
c	open (unit=11,file='ksat.txt',status='unknown')
c	open (unit=13,file='rr.txt',status='unknown')
c	open (unit=16,file='man.txt',status='unknown')
c	open (unit=17,file='cover.txt',status='unknown')
	open (unit=18,file='theta.inp',status='unknown')

c	read input data

	read (1,1) simdur,dt
    1 format(/,1x,i5/,1x,f6.0)
	read (1,*) nsoil
	read(1,*) tp(1),fc(1),wp(1),a(1),df(1)
     	read(1,*)cl(1),sa(1),st(1),wcf(1),vfs(1)	
	read(1,*) rrch,  mnch	       
	read(1,*) albed(1)
	read(1,*) dx  
	read(1,*) ninterp, nmaizep
	read(1,*) dinterp(1:ninterp)
	read(1,*) destruction, burial
	read(1,*) dmaizep(1:nmaizep)
	read(1,*) wheelinter(1:ninterp) , wheelmaize(1:nmaizep)
	read(1,*) rrintini(1:ninterp) 
	read(1,*) rrmaizeini(1:nmaizep)  
	read(1,*) kinio, kinim 
	read(1,*) rrksvarint(1:ninterp) 
	read(1,*) rrksvarmaize(1:nmaizep) 
	read(1,*) kextint, kextmaize 
	read(1,*) gcoverint, kcovermaize
	read(1,*) rbiomcover
	     
	read(10,*) kso2, kcho2, kchm2, hfront,c_kso, 
     .	c_rr, man2
  
	read(3,*) n,nrow,ncol

	do i=1,n
	read(3,*) numpix(i),xcoord(i),ycoord(i),sl(i),
     . 	fdir(i),
     .	flacc(i)
     .	,spix(i),cpix(i),wheel(i)	
	end do	
	count=0
	do i=1,nrow
		do j=1,ncol
		count=count+1
		xycoord(i,j)=count
		end do
	end do	
	do i=1,n
c		read(11,*) ksat(i)
		ksat(i)=kso2
	end do
	do i=1,n
c		read(13,*) rr(i)
		rr(i)=0
	end do
	do i=1,n
c		read(16,*) mann(i)
		mann(i)=man2
	end do
      do i=1,n
c		read(17,*) cover(i)
		cover(i)=0.0
	end do
	do i=1,n
		read(18,*) theta(i)
	end do
	
	read(9,*) As, d50, Coh, Cohw 

	close(1)
	close(3)
!	close(9)
	close(10)
	close(11)
	close(13)	
	close(16)
	close(17)

	dxvert=dx
	dxdiag=dx*2**(0.5)

	neff=0
	
	do i=1,n
	if (sl(i).gt.(-9999)) then
	sl(i)=sl(i)/100	! slope introduced in %
	neff=neff+1	
	thetainst(i)=theta(i)
	a(i)=a(spix(i))
	df(i)=df(spix(i))
	tp(i)=tp(spix(i))
	wp(i)=wp(spix(i))
	fc(i)=fc(spix(i))
	dr(i)=0.
	s(i)=0.
	albed(i) = albed(cpix(i))	
	sa(i) = sa(spix(i))
	cl(i)= cl(spix(i))
	st(i)=st(spix(i))
	ksend(i) = kso2		 
	psi(i)=hfront		
	cfact(i)=c_kso		
	end if
	end do
	
	rrexp=c_rr
		
	dx2=dx*dx
	cu=dx2/3.6e+6	! mmhr-1 to m3s-1
	tcst=dt/3600	! mmhr-1 to mm

	do i=1,n
	if (sl(i).gt.(-9999)) then
	wp(i)=wp(i)*tp(i)
	fc(i)=fc(i)*tp(i)
	end if
	end do	

c	variables for evapotranspiration, Ritchie's method (1972)
c	upper stage of soil evaporation (mm)
	do i=1,n
	if (sl(i).gt.(-9999)) then
	alphas(i) = (4.165+0.02456*sa(i)-0.01703*cl(i)-0.0004*sa(i)*sa(i))
	esu(i)= (0.9*(alphas(i)-3)**0.42)*10		! mm

c	maximum evaporative depth (mm)
	edx(i)=90.-0.77*cl(i)+0.006*sa(i)*sa(i)

c	intializing of the evaporation variables
	es1p(i)=0.
	es2p(i)=0.
	es(i)=0.
	pep(i)=0.
	asmlim(i) = wp(i) !+ (fc(i) - wp(i))*.25
	end if
	end do

c	Erosion variables pre-processing

c	Vs = settling velocity - stokes law spherical particle Chow et al. 1988
	rhop=2650
	rhow=1000
	dclay=0.002 !diam in mm
	dsilt=0.01
	dsand=0.2
	dvfs=0.03
	dcf=0.5
	clayfract=cl(1)/100			!0.15
	siltfract=st(1)/100			!0.25
	sandfract=sa(1)/100			!0.50
	vfsfract=vfs(1)/100			!0.10
	cffract=wcf(1)/100			!0
	nudynwat=1.3e-3 !Nsm-2

	Diam=(dclay*clayfract+dsilt*siltfract+dsand*sandfract+dvfs
     .	*vfsfract+cffract*dcf)*1e-3  ! m

	Vs=	10*(2./9)*((Diam*0.5)**2)*9.81*(rhop-rhow)/nudynwat !m/s
	
	ctc=((d50+5)/0.32)**(-0.6)
	dtc=((d50+5)/300)**(0.25)

!	Coh: 0.1 kg/cm2 = 9.81 kPa

	Ydet=1./(0.89+0.56*Coh)		!0-1
	Ydetw=1./(0.89+0.56*Cohw)		!0-1

c	Cover degradation constants

	fN=0.57
	fW=0.2
	kdeg=-0.01
	
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

	do 11 idate=1,simdur

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
	if (idate.eq.1) then		! first day inter1
	write(*,*) ' Intercrop 1'

!	Kcover & Kext
	kccover=gcoverint
	kext=kextint
	prstmaize=0
	dstart=dinterp(1)


	! RR & KS & manning

		!	Ks & RR overland decreasing

		if (rrksvarint(1).eq.0) then

			rrksovar=0

		else
			rrksovar=1

	    end if

		do i=1,n
			if (sl(i).gt.(-9999)) then

			if (wheel(i).eq.0) then

c	random roughness
			
			rr(i)=rrintini(1)
			rrcm = rr(i)/10.	! rrm in cm

c	surface retention, Mwandera & Feyen 1992

			dir2(i)= (0.294*rrcm+0.031*rrcm*rrcm-0.012*rrcm*sl(i)
     .			*100)*10.

			if (rrksovar.eq.0) then

				ke(i)=ksend(i)
			else
				ke(i)=kinio
			end if

			else
			rr(i)=rrch
			rrmc=rrch/10.

			dir2(i)=(0.294*rrmc+0.031*rrmc*rrmc-0.012*rrmc*sl(i)
     .			*100)*10.

		   
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

	

c----------------------------------------------------------
	elseif (idate.eq.dmaizep(1)) then		! maize1
	write(*,*) ' Maize 1'

	!	Kcover & Kext
	kccover=kcovermaize
	kext=kextmaize
	prstmaize=1
	dstart=dmaizep(1)



	! RR & KS & manning

		!	Ks & RR overland decreasing

		if (rrksvarmaize(1).eq.0) then

			rrksovar=0

		else
			rrksovar=1

	    end if

		do i=1,n
		if (sl(i).gt.(-9999)) then

			if (wheel(i).eq.0) then
			
			rr(i)=rrmaizeini(1)
			rrcm = rr(i)/10.	! rrm in cm

			dir2(i)= (0.294*rrcm+0.031*rrcm*rrcm-0.012*rrcm*sl(i)*100) 
     .			*10.

			if (rrksovar.eq.0) then

				ke(i)=ksend(i)
			else
				ke(i)=kinim
			end if

			else
			rr(i)=rrch
			rrmc=rrch/10.

			dir2(i)=(0.294*rrmc+0.031*rrmc*rrmc-0.012*rrmc*sl(i)*
     .			100)*10.

		   
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

	

c-------------------------------------------------------------------------------
	elseif (idate.eq.dinterp(2)) then  !inter2
	write(*,*) ' Intercrop 2'

	!	Kcover & Kext
	kccover=gcoverint
	kext=kextint
	prstmaize=0
	dstart=dinterp(2)


	! RR & KS & manning

		!	Ks & RR overland decreasing

		if (rrksvarint(2).eq.0) then

			rrksovar=0

		else
			rrksovar=1

	    end if
		do i=1,n
          if (sl(i).gt.(-9999)) then

			if (wheel(i).eq.0) then
			
			rr(i)=rrintini(2)
			rrcm = rr(i)/10.	! rrm in cm

			dir2(i)= (0.294*rrcm+0.031*rrcm*rrcm-0.012*rrcm*sl(i)
     .			*100)*10.

			if (rrksovar.eq.0) then

				ke(i)=ksend(i)
			else
				ke(i)=kinio
			end if

			else
			rr(i)=rrch
			rrmc=rrch/10.

			dir2(i)=(0.294*rrmc+0.031*rrmc*rrmc-0.012*rrmc*sl(i)*
     .			100)*10.

		   
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

	
c--------------------------------------------------------------------------
	elseif (idate.eq.dmaizep(2)) then	! maize2
	write(*,*) ' Maize 2'

	!	Kcover & Kext
	kccover=kcovermaize
	kext=kextmaize
	prstmaize=1
	dstart=dmaizep(2)

	! RR & KS & manning
		!	Ks & RR overland decreasing

		if (rrksvarmaize(2).eq.0) then

			rrksovar=0

		else
			rrksovar=1

	    end if
		do i=1,n
		if (sl(i).gt.(-9999)) then

			if (wheel(i).eq.0) then
		
			rr(i)=rrmaizeini(2)
			rrcm = rr(i)/10.	! rrm in cm

			dir2(i)= (0.294*rrcm+0.031*rrcm*rrcm-0.012*rrcm*sl(i)
     .			*100)*10.

			if (rrksovar.eq.0) then

				ke(i)=ksend(i)
			else
				ke(i)=kinim
			end if

			else
			rr(i)=rrch
			rrmc=rrch/10.

			dir2(i)=(0.294*rrmc+0.031*rrmc*rrmc-0.012*rrmc*sl(i)
     .			*100)*10.
		   
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

c-------------------------------------------------------------------------------
	elseif (idate.eq.dinterp(3)) then  !inter2
	write(*,*) ' Intercrop 3'

	!	Kcover & Kext
	kccover=gcoverint
	kext=kextint
	prstmaize=0
	dstart=dinterp(3)

	! RR & KS & manning

		!	Ks & RR overland decreasing

		if (rrksvarint(3).eq.0) then

			rrksovar=0

		else
			rrksovar=1

	    end if
		do i=1,n
		if (sl(i).gt.(-9999)) then

			if (wheel(i).eq.0) then
			
			rr(i)=rrintini(3)
			rrcm = rr(i)/10.	! rrm in cm

			dir2(i)= (0.294*rrcm+0.031*rrcm*rrcm-0.012*rrcm*sl(i)
     .			*100)*10.

			if (rrksovar.eq.0) then

				ke(i)=ksend(i)
			else
				ke(i)=kinio
			end if

			else
			rr(i)=rrch
			rrmc=rrch/10.

			dir2(i)=(0.294*rrmc+0.031*rrmc*rrmc-0.012*rrmc*sl(i)*
     .			100)*10.

		   
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

	
c--------------------------------------------------------------------------
	elseif (idate.eq.dmaizep(3)) then	! maize2
	write(*,*) ' Maize 3'

	!	Kcover & Kext
	kccover=kcovermaize
	kext=kextmaize
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

			if (wheel(i).eq.0) then
		    
			rr(i)=rrmaizeini(3)
			rrcm = rr(i)/10.	! rrm in cm

			dir2(i)= (0.294*rrcm+0.031*rrcm*rrcm-0.012*rrcm*sl(i)
     .			*100)*10.

			if (rrksovar.eq.0) then

				ke(i)=ksend(i)
			else
				ke(i)=kinim
			end if

			else
			rr(i)=rrch
			rrmc=rrch/10.

			dir2(i)=(0.294*rrmc+0.031*rrmc*rrmc-0.012*rrmc*sl(i)
     .			*100)*10.
		    
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
		
c---------------------------------------------------------
	end if
c---------------------------------------------------------

c	End of case study specific implementation

c	random roughness evolution

	do i=1,n
		
	thetainst(i)=theta(i)
	rrv(i)=rr(i)*exp(-rrexp*ektotbelcovtot)
	rrcm = rrv(i)/10.

c	corresponding retention
	
	if (rrksovar.eq.1) then
	dirvon(i) = (0.294*rrcm+0.031*rrcm*rrcm-0.012*rrcm*sl(i)*100)*10. ! mm	
	else
	dirvon(i)=dir2(i)
	end if

	end do

	do m=1,n	
			
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
	 s(m)=0.
	 ss(m)=0.
	 hsv(m)=0.
	 hsv1221(m)=0.
	 cuminf(m)=0.
	 q(m)=0
	 draintot(m)=0

c	erosion variables
	 qs11(m)=0.
	 qs12(m)=0.
	 qs21(m)=0.
	 qs22(m)=0.
	 qs2t1(m)=0.
	 qs2t2(m)=0.	

	if(m.lt.90001) rcc(m)=0.
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

c	read weather data
	
	read(2,*) day,tempc,soitemp,radi,raites

c	cover developement (percentage cover, biomass, LAI) rountines

	if (tempc.lt.0) tempc=0
	sumtemp=sumtemp+tempc

	if (kccover.gt.0) then

		if (prstmaize.gt.0) then ! maize

		prstcover=0.01*(95*1*exp(kccover*sumtemp)/			! verhulst equation 
     .		   (95+1*(exp(kccover*sumtemp)-1)))		    ! cover=f(sumtemp)

		biomcover = 0
		lstbiomcover = 0
		
		else					 ! intercropping	period

		if (destruction.gt.0) then		! there is a destruction date
		
		checkd = dstart+destruction

		if (idate.lt.(dstart+destruction)) then ! before destruction

			sumtempdeg=0

			biomcover =	kccover*(sumtemp**2)	! g/m�= f(sumtemp(0) �C) - expon. curve
			lstbiomcover = biomcover

			prstcover =  -1*(exp(-rbiomcover*biomcover)-1)

		elseif (burial.gt.0) then ! there is a destruction and a burial

			checkd=dstart+destruction+burial

			if (idate.lt.(dstart+destruction+burial)) then	! degradation

				sumtempdeg=sumtempdeg+tempc

				biomcover =	exp(fN*fW*kdeg*sumtempdeg)*lstbiomcover		! douglas & rickman

				prstcover =  0.004*biomcover		! experimental
				if (prstcover.gt.0.95) prstcover=0.95

			else		! buried

				prstcover=0
				biomcover =0
			    
			end if
		
		else		! degradation and no burial

			sumtempdeg=sumtempdeg+tempc

			biomcover =	exp(fN*fW*kdeg*sumtempdeg)*lstbiomcover

			prstcover =  0.004*biomcover 
			if (prstcover.gt.0.95) prstcover=0.95

		end if		! before/after destruction loop

		else		! no destruction

		biomcover =	kccover*(sumtemp**2)	! g/m�= f(sumtemp(0) �C) -  expon. curve	

		prstcover =  -1*(exp(-rbiomcover*biomcover)-1)

		end if ! destruction y/n loop

		end if ! intercropping / maize	period loop

	else

	prstcover=0
	biomcover =0
	lstbiomcover = 0

	end if	
	
	if (kext.gt.0) then		

	prstlai= -log(1-prstcover)/kext ! Andrieu et al. 1997

	else

	prstlai=0

	end if

	if (raites.eq.0) goto 12	! no rainfall

c	raining day

c	read rainfall data

	call rainfa

c	compute ksat evolution

	do i=1,n

	if (wheel(i).eq.0) then

		if (rrksovar.eq.0) then

			keff(i)=ksend(i)


		elseif (prstcover.gt.0.01) then
	
		kcover2(i) = ke(i)*((ksend(i)/ke(i)) + (1. - (ksend(i)/ke(i)))
     .		*exp(-cfact(i)*ektotbelcovtot*(1.- (rrv(i)/40.))))
	
			keff(i) = kcover2(i)			
		else

		keb2(i) = ke(i)*((ksend(i)/ke(i)) + (1. - (ksend(i)/ke(i)))*
     .		exp(-cfact(i)*ektotbeltot*(1.- (rrv(i)/40.))))

				!	Rise et al. 1995.

			keff(i) = keb2(i)

		end if

	else

		keff(i) =	ke(i)	
		
	end if
		
	end do

c	start infiltration-runoff-rainfall routine

	do 13 i=1,lline

	write(*,*) 'Time = ',i*0.5,' min'
	
	if (prstlai.eq.0.or.rc(i).eq.0.) then
	rcc(i)=rc(i)
	else	
	rcc(i)=rain(rc(i),prstlai,i,dt,tc(i),idate,cu,prstcover)
	end if
	
	prstr=rcc(i)
	rtc(i)=tc(i)-dtm

	totprec=totprec + rc(i)*tcst
	pinst= rc(i)*tcst
	totprecgr=totprecgr + rcc(i)*tcst
c	kinetic energy content of the rain. Van Dijk et al. 2002. metric units

	ek = 28.3*(1-0.52*exp(-0.042*rc(i))) ! J m-2 mm-1
	ektot = ektot +	 ek*rc(i)*tcst	! J m-2

c	kinetic energy content of the rain for belgium. Bollinne et al. 1980
c				ekbel = emax*(1-a*exp(-b*rainrate)) J m-2 mm-1

	ekbel = 29*(1-0.6*exp(-0.061*rc(i)))

	ekbelsplash=ekbel*rc(i)*tcst	! J m-2

	ektotbel = ektotbel + ekbel*rc(i)*tcst	! J m-2

	ekbelcov = 29*(1-0.6*exp(-0.061*rc(i)))
	ekbelcov2 = 29*(1-0.6*exp(-0.061*rcc(i))) ! J m-2 mm-1

	Tra = ekbelcov*rc(i)*tcst*(1.-prstcover)

c maize stem and canopy flow

	if (prstmaize.gt.0) then

		if (idate.lt.(dstart+31)) then
		height=0.5
		elseif (idate.lt.(dstart+76)) then
		height=1
		else
		height=1.5
		end if

		brandt=15.8*(height**0.5)-5.87

		Trf= rcc(i)*tcst*(prstcover)*brandt*0.5

	else

		Trf=0

	end if

	ektotbelcov = ektotbelcov + Tra + Trf

	timnum2=timnum2+1
	if (timnum2.gt.1.) lstr = rcc(i)

c	compute inflitration and drainage
	
c	do 140 m2=0,maxval(flacc)

	do 14 m=1,n

	if (sl(m).gt.(-9999)) then
		
		if (keff(m).eq.0.) then

				filt(m)=0.		
		else
	
		filt(m)=gamptun3(hsv(m),rcc(i),keff(m),dt,rtc(i),tc(i),
     .					m,n,n2,idate,cu,df(m),tp(m),theta(m),dx2,
     .					psi(m))

		end if
	
		if (thetainst(m).gt.fc(m)) then
			dr(m)=drain(thetainst(m),keff(m),tp(m),fc(m),df(m),dt,
     .		a(m),idate)
		else

			dr(m) = 0.

		end if

	else

	filt(m)=-9999
	dr(m)=-9999

	end if

c	end if   ! flacc=m2	

   14 continue
c  140 continue

c	compute runoff 150 - 15

	outlet=maxval(flacc)

	do 150 m2=0,maxval(flacc)

		do 15 m=1,n		!avant : n2

	if (sl(m).gt.(-9999)) then
	
	if (flacc(m).eq.m2) then   ! pixel corresponding to flacc

	iw=int(m)
	iw2=int(m)

c	iw=indice m for kinematic wave	
	
      fil=filt(m)*cu		! m3/s
	
	if (fil.lt.0.) fil=0.

	if (s(m).gt.0.or.rcc(i).gt.0.) then
	
	cuminf(m)=cuminf(m) + filt(m)*dt/3600

	end if

c	update thetainst after infiltration and drainage

	vertflux = rcc(i) + s(m)*3600/dt

	if (filt(m).gt.vertflux) then

		filtreal = vertflux

	else

		filtreal = filt(m)

	end if	
	
	thetainst(m)= thetainst(m)+(((filtreal-dr(m))/3600)*dt)/df(m)

	draintot(m) = draintot(m) + dr(m)*dt/3600
	
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

			if (xycoord(x,y).eq.m) then

			if (ycoord(m).ne.1.and.ycoord(m).ne.ncol.and.
     .			xcoord(m).ne.1.and.xcoord(m).ne.nrow) then ! pixel de l'int�rieur
			north=xycoord(x-1,y)
			northe=xycoord(x-1,y+1)
			east=xycoord(x,y+1)
			southe=xycoord(x+1,y+1)
			south=xycoord(x+1,y)
			southw=xycoord(x+1,y-1)
			west=xycoord(x,y-1)
			northw=xycoord(x-1,y-1)
			fdn=fdir(north)
			fdne=fdir(northe)
			fde=fdir(east)
			fdse=fdir(southe)
			fds=fdir(south)
			fdsw=fdir(southw)
			fdw=fdir(west)
			fdnw=fdir(northw)
			
			elseif (m.eq.xycoord(1,1)) then  ! pixel de l'ext�rieur
			southe=xycoord(x+1,y+1)
			south=xycoord(x+1,y)
			east=xycoord(x,y+1)
			fde=fdir(east)
			fdse=fdir(southe)
			fds=fdir(south)
			fdn=0
			fdne=0
			fdsw=0
			fdw=0
			fdnw=0

			elseif (m.eq.xycoord(1,ncol)) then  ! pixel de l'ext�rieur
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

			elseif (m.eq.xycoord(nrow,1)) then  ! pixel de l'ext�rieur
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

			elseif (m.eq.xycoord(nrow,ncol)) then  ! pixel de l'ext�rieur
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

			elseif (ycoord(m).eq.1) then
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
		
			elseif (ycoord(m).eq.ncol) then

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
			
			elseif (xcoord(m).eq.1) then

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

			elseif (xcoord(m).eq.nrow) then
			
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

			end if  ! pixel int/ext

			if (fdn.eq.4) then
			qintot11=qintot11+q2t1(north)
			qintot12=qintot12+q2t2(north)
			qsintot11=qsintot11+qs2t1(north)
			qsintot12=qsintot12+qs2t2(north)
			q(north)=1
			end if

			if (fdne.eq.8) then
			qintot11=qintot11+q2t1(northe)
			qintot12=qintot12+q2t2(northe)
			qsintot11=qsintot11+qs2t1(northe)
			qsintot12=qsintot12+qs2t2(northe)

			q(northe)=1
			end if

			if (fde.eq.16) then
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

			qp11(m)=qintot11		
			qp12(m)=qintot12	
			qi(m)=qintot11
			qi2(m)=qintot12
			qs11(m)=qsintot11
			qs12(m)=qsintot12

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
		
	qidt = (qi(m) + qi2(m))/2 
	hindt=(qidt*dt/dx2)*1000

	if (hsv(m).gt.0) hindt=0.

	s(m) = ss(m) + (rcc(i) - filt(m))*cu*1000*dt/dx2  + hindt

	if (s(m).lt.0.) s(m) = 0.

	if (i.eq.1) then
		oldqbar(m) = 0.
		else
		oldqbar(m) =oldqtrans
	end if
	
      if (s(m).gt.dirvon(m)) then
	
	dxf=dxvert
					
		goto 22
	
	end if  

c *******	no runoff *******

	ss(m)=(rcc(i) - filt(m))*cu*1000*dt/dx2 + ss(m)

	if (ss(m).lt.0.) ss(m)=0.
	hsv(m)=0.
	hsv1221(m)=0.
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
	h21=0
	h22=0
	h11=0
	h12=0
	Dfl=0
	Transc=0
	Dep=0
	eros=0
	qs2t1(m)=0.
	qs2t2(m)=0.
	goto 23

   22 continue
c     ******Kinematic wave******
	
c	*** overland & channel flow cell *** 
	
	omega=0.75
	beta=0.6
	tol2=1.0e-7

c	handling of diagonal flow
	if(fdir(m).eq.128.or.fdir(m).eq.2.or.fdir(m).eq.8.or.fdir(m).eq.32)
     .	 then
			dx=dxdiag
		else
			dx=dxvert
	end if
	
	weightr=1
	qin11=oldqbar(m)/dx
	qin21=oldqbar(m)/dx	
	lat=(rcc(i)*cu-fil)	
	qin12=(lat*weightr/dx)	
	qin22=(lat*weightr/dx)
	
	if (qin21.lt.0.) qin21=0   
	
	if (qin11.lt.0.) qin11=0   
	
c***** alpha *********
	if (wheel(m).gt.0) then	!channel cell
		
	slwp=(sl(m))
		
	width2(m)=dxvert		!chwid
	
	hmin= -10*tol2

	hwp11=(((qp11(m)*mann(m)*(width2(m)**0.666667))/(sqrt(slwp)))**0.6
     .	)/dx
	hwp12=(((qp12(m)*mann(m)*(width2(m)**0.666667))/(sqrt(slwp)))**0.6
     .	)/dx
	hwp21=(((qp21(m)*mann(m)*(width2(m)**0.666667))/(sqrt(slwp)))**0.6
     .	)/dx
	

	hwp=hwp11
	qact=qp11(m)
	if (qact.eq.0.) then
		hwp11=0.
		go to 2112
	end if
	 
	call iter(errhwp, hwp, hmin, 10., ierr)
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

 2114	wper11= chwid + 2*hwp11
	wper12= chwid + 2*hwp12
	wper21= chwid + 2*hwp21
	wper22=wper21


	alp11(m) = (mann(m)*(wper11**0.666667)/(sqrt(slwp)))**0.6
	alp12(m) = (mann(m)*(wper12**0.666667)/(sqrt(slwp)))**0.6
	alp21(m) = (mann(m)*(wper21**0.666667)/(sqrt(slwp)))**0.6
	alp22(m) = (mann(m)*(wper22**0.666667)/(sqrt(slwp)))**0.6
	
	else	! overland flow cell	
	
	alp= (mann(m)*(dxf**0.666667)/sqrt(sl(m)))**0.6
	
	wper11=dxf 
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
	                          
            if (qtest .le. 1.e-12 .and. (qin22) .le. 1.e-13) then

              qp22(m) = 0.
		  			
            else
									    								                              						
				if (wheel(m).gt.0.and.timnum.eq.0.) then
					q22=1.0e-9
				else
					q22=(qp12(m) + qp21(m))/2
				end if
			
				if (q22.eq.0.) q22=1.0e-9

				!qmin = min(-q22,-tol2)
				qmin = -tol2*10
				 							
				call iter (errq2, q22, qmin, 10., ierr)
	!		
				if (ierr .gt. 0) then
					write(*,*) ' no convergence for cell ',m,
     .				' at idate = ',idate
					read(*,*)
                  end if

c	 to avoid to low outflow
				qouttest = 0.	!1.e-12
		
				if (q22 .lt. qouttest) q22 = 0.
				
				if (q22.gt.0.) then
					lastru =1
				else
					lastru =0
				end if

				qp22(m) = q22
							
	      end if
	
	q2t1(m)=qp21(m)
	q2t2(m)=qp22(m)
	q1t1(m)=qp11(m)
	q1t2(m)=qp12(m)		   

	q2=q2t2(m)

	if (flacc(m).eq.maxval(flacc)) then ! outlet
		
	timnum = timnum+1
	
			qoutlet=q2	

	end if ! outlet

	h22=(((q2t2(m)*mann(m)*(wper22**0.666667))/(sqrt(slwp)))**0.6)/
     .											width2(m)
	if (m.gt.n.and.h22.ne.0.) then
	hwp=h22
	qact=q2t2(m)
	call iter(errhwp, hwp, hmin, 10., ierr)
	h22=hwp
	end if

	h12=(((qp12(m)*mann(m)*(wper12**0.666667))/(sqrt(slwp)))**0.6)/
     .						width2(m)
	
	h21=(((qp21(m)*mann(m)*(wper21**0.666667))/(sqrt(slwp)))**0.6)/
     .						width2(m)

	h11=(((qp11(m)*mann(m)*(wper11**0.666667))/(sqrt(slwp)))**0.6)/
     .						width2(m)
	
	hf=(h12+h21)/2

	hf2=(h21+h22)/2

	hf3=(h12+h22)/2

c	 update initial condition for next dt

	qp21(m)=qp22(m)

	if (m.le.n) oldqtrans = rcc(i)*cu-fil
	
c	***** end of Kinematic wave *****

	hsv1221(m)=hf*1000
	hsv(m)=hf3*1000
	ss(m)=hf3*1000+dirvon(m)	
 
 3766 format(1x,i5,1x,f7.2,1x,f7.2,1x,f10.4)
 3769 format(1x,i5,1x,f7.2,1x,f7.2,1x,f10.4,1x,f15.10,1x,f15.6)
 3770 format(1x,i5,1x,f7.2,1x,f7.2,1x,f10.4,1x,f15.10,1x,f15.3,1x,     
     .	f15.10,1x,f15.10,1x,f15.10,1x,f8.4)
 3788 format(1x,i5,1x,i5,1x,f7.2,1x,f7.2,1x,f8.4)	
   

 8200 format(i3,1x,i3,1x,f4.1,1x,f10.8,1x,f10.8,1x,f10.8)

!	end if   ! if pixel ~ flacc

c	erosion loop

c	!!! if no erosion
c	goto 99

	if (hsv1221(m).gt.0.01) then		
			
	hmm=hsv1221(m)
	hf=hmm/1000

	if (rcc(i).gt.0) then
	Dsplash = Lerosionsplash(dt,dx2,hmm, ekbel,ekbelcov2,rc(i),rcc(i),
     .	As, prstcover) !kg/s   

	else
	Dsplash =0.
	end if
	
	vel=0.5*(q2t1(m)+q1t2(m))/(hf*dxf)	

	wstream=vel*sl(m)*100			!cm/s

	if ((wstream).gt.0.4) then		!cm/s			

	Transc=ctc*((wstream-0.4)**dtc)*2650		

	else

	Transc=0

	end if

	conc=(qs21(m)+ qs12(m))/(q2t1(m) + q1t2(m))

	if (isnan(conc)) conc=0.

	if (Transc.gt.conc) then		! Dfl
	
	if (wheel(m).gt.0) then

	Dfl=Ydetw*(Transc-conc)*Vs*dxf*dx

	else

	Dfl=Ydet*(Transc-conc)*Vs*dxf*dx

	end if

	Dflmax=(Transc-conc)*(q1t2(m)+q2t1(m))*0.5

	if (Dfl.gt.Dflmax) Dfl=Dflmax

	Dep=0
	
	else			!Dep
	
	Dfl=0

	Dep=(Transc-conc)*Vs*dxf*dx


	if (abs(Dep*dt).gt.abs((Transc-conc)*dxf*dx*hf)) then
 
	
	Dep=((Transc-conc)*dxf*dx*hf)/dt

	end if

	end if
	
	eros=Dep+Dfl+Dsplash  ! kg/s
				
! Lisem model:

	alpha1221=(alp12(m)+alp21(m))*0.5
	qp1221=(q1t2(m)+q2t1(m))*0.5

	ebar=eros/dx
		
! Lisem model:	
	qs22(m)=(ebar*dx+alpha1221*(qp1221**beta)*(dx/dt)*(qs21(m)/q2t1(m))
     .		- conc*alpha1221*beta*(qp1221**(beta-1))*(dx/dt)*
     .		(q2t2(m)-q2t1(m)) + qs12(m))/(1+alpha1221*(qp1221**beta)*
     .		(dx/dt)*(1./q2t2(m)))

	concmax=Transc
	! conc max = 848 kg m-3 in Morgan et al. (1997)
	
	if (qs22(m).lt.0) qs22(m)=0.
	
	if (flacc(m).eq.maxval(flacc)) then ! outlet

		ccoutlet = qs21(m) 
	
	end if ! coutlet

c	 update initial condition for next dt

	qs21(m)=qs22(m)
	qs2t2(m)=qs22(m)
	qs2t1(m)=qs21(m)
	
	else ! no eros

	Dfl=0
	Transc=0
	Dep=0
	eros=0
	qs21(m)=0
	qs22(m)=0
	qs2t1(m)=0.
	qs2t2(m)=0.

c	here you can print the water balance and erosion information associated with pixel m within the event : theta(m), qp22(m) m3/s, qs22(m) kg/s,...
	
	end if		! end erosion loop

   99 continue
   23 continue
   
	end if   ! if pixel ~ flacc   
	
	else     ! sl(m)==-9999

	end if	

   15 continue	! end runoff loop
	
  150 continue	! end runoff loop
	
	qmmhr = (qoutlet/(dx2*neff))*1000*3600			

	totrun=totrun+ qmmhr*dt/3600

	toteros=toteros+ccoutlet*dt
	
	write(379,3769) idate,rc(i),(timnum2)*(dt/60.),
     .           	   qmmhr, qoutlet, hmm

	write(380,3771) idate,lstr,(timnum2-1)*(dt/60.),
     .           	   ccoutlet, Transc, ccoutlet/qoutlet, Dep, Dfl, 
     .				Dsplash,
     .				eros
 3771 format(1x,i5,1x,f7.2,1x,f7.2,1x,f10.4,1x,f15.10,1x,f15.3,1x,     
     .	f15.10,1x,f15.10,1x,f15.10,1x,f10.4)

	if (qmmhr.eq.0) timnum=0
	
	qoutlet=0.
	ccoutlet=0.

	if (i.eq.lline) lsttr = tc(i)

   13	continue

   12 continue

c	new soil moisture content after event

	do i=1,n
c	theta(i)= thetainst(i)
	if (theta(i).lt.wp(i)) theta(i) = wp(i)
	if (theta(i).gt.tp(i)) theta(i) = tp(i)
	end do

c	evaporation, transpiration & percolation routine

	do i=1,n

	call etp (tempc,radi,albed(i),etpmm) 

	ie = i
	
	call etpr2(prstlai,totprecgr,idate)

	thetbefetr=theta(i)
	
	edisp =	(theta(i))*df(i)			!(theta(i)-wp(i))*df(i)

	edispm = min1((theta(i))*edx(i),edisp)    !(theta(i)-wp(i))*edx(i)

	edispmfc = min1(fc(i)*edx(i),fc(i)*df(i))

	if (es(i).gt.edispm) es(i) = edispm

	if (theta(i).gt.fc(i)) then
c	soil moisture not limiting,drainage

		tetpr = es(i) + pep(i)

		percol = 0.

		call perco(theta(i),keff(i),tp(i),fc(i),df(i),
     .				dt,a(i),idate,percol,lsttr)

		!percol = percol*(24-lsttr/60.)

		percol = percol*24.

		totwloss = tetpr + percol

		if (edisp.gt.totwloss) then

			theta(i) = (edisp - totwloss)/df(i)		!+ wp(i)

		else

			totwloss = edisp

			theta(i) = (edisp - totwloss)/df(i)		!+ wp(i)

		end if

	if (theta(i).lt.wp(i)) theta(i) = wp(i)

	else
c	soil moisture limiting, no drainage
		percol=0.

		if (theta(i).lt.asmlim(i)) then

			pep(i) = pep(i) * theta(i)/asmlim(i)

			tetpr = es(i) + pep(i)

			if (edisp.gt.tetpr) then

				theta(i) = (edisp - tetpr)/df(i)	!+ wp(i)

				if (theta(i).lt.wp(i)) theta(i) = wp(i)

			elseif (theta(i).eq.wp(i)) then	           
						
				es(i)=0.
				pep(i)=0.

		    else
			
				es(i) = (es(i)/(es(i)+pep(i))) * edisp/tetpr

				pep(i) = (pep(i)/(es(i)+pep(i))) * edisp/tetpr

				tetpr = edisp

				theta(i) = (edisp - tetpr)/df(i)	!+ wp(i)
			
				if (theta(i).lt.wp(i)) theta(i) = wp(i)

			end if	

		else

			tetpr = es(i) + pep(i)

			if (edisp.gt.tetpr) then

				theta(i) = (edisp - tetpr)/df(i)	!+ wp(i)

				if (theta(i).lt.wp(i)) theta(i) = wp(i)

	       elseif (theta(i).eq.wp(i)) then
	           						
				es(i)=0.
				pep(i)=0.

		   else
			
				es(i) = (es(i)/(es(i)+pep(i))) * edisp/tetpr

				pep(i) = (pep(i)/(es(i)+pep(i))) * edisp/tetpr

				tetpr = edisp

				theta(i) = (edisp - tetpr)/df(i)	!+ wp(i)

				if (theta(i).lt.wp(i)) theta(i) = wp(i)
			
			end if	

	    end if	

	end if
	
	end do
	
	do i=1,n

	if (theta(i).lt.wp(i)) theta(i)=wp(i)

	end do

c	end of percolation, transpiration & evaporation routines

c	here print the runoff and erosion information associated with pixel of interest after the event

	write(375,3777) idate,totprec,totrun,(toteros/(neff*dx2))*10	
 3777 format(1x,i5,f8.4,f8.4,' mm',1x,f16.4,' t/ha')
	write(378,3780) idate,totprec,totrun,(toteros/(neff*dx2))*10	
 3780 format(1x,i5,f8.4,1x,f8.4,1x,f16.4)
	totruntot = totrun+totruntot
	toterostot=toteros+toterostot
	totprectot = totprec + totprectot
	totprecgrtot = totprecgr + totprecgrtot
	ektottot=ektottot + ektot
	ektotbeltot=ektotbeltot+ektotbel
	ektotbelcovtot=ektotbelcovtot+ektotbelcov
   11 continue

	write(375,3778) totruntot,totprectot
 3778 format(1x,'total runoff = ',f8.4,' mm '/,1x,'total rainfall = ',
     .	   f8.4,' mm ')

	call cpu_time(t1)
	write(*,*) 'Time of operation was ', t1 - t0, ' seconds'
	!read(*,*)
      end


c********************Functions nd subroutines********************
c---------------------------------------------------------------------------------------------
	subroutine rainfa
	implicit real (a-z)
	integer k,rline,i,j,l,lline,jbeg,idate,simdur,n,ncol,nrow,
     .	nrowout,n2
	real rcs(3000),tcs(3000)
	real beta,omega
	common /gen/ simdur,ncol,nrow,nrowout,n,n2
	common /gen2/ dt,dx
	common /rainfall/ tc(90000),rc(90000),dtm,lline,idate
  
	dtm=dt/60.
	tc(1)=0
	k=1
	l=0
	stram=0.

	tcs(k)=0.
	rcs(k)=0.
	k=2
   11 read(2,*) day,tcs(k),rcs(k),zero1,jbeg,zero2

	k=k+1
	if (jbeg.eq.0) goto 11
      rline=k-1
	do i=2,rline
	tinc=tcs(i)-tcs(i-1)
	kt=tinc/dtm
	
   	do j=1,int(kt)
	l=l+1
	lline=l
	rc(l)=rcs(i)
	if (i.eq.2) then
	tc(l)=tcs(2)-dtm+(j-1)*dtm
	else
	tc(l)=tbeg+(j-1)*dtm
	end if
	end do
	tbeg=tc(l)+dtm
	end do
	end 

c------------------------------------------------------------------------------
c------------------------------------------------------------------------------
	function rain(rate,prstlai,i,dt,t,idate,cu,prstcover)
	
	implicit real (a-h,o-z)
	real pcum, cint, lstcint,inrate
	rbmin=rate/3600.		!mm/s
	p=1-0.046*prstlai
	if (t.eq.0.or.t.eq.0.5) then
	pcum=0.
	lstcint=0.
	cint=0.
	if (t.eq.0.) goto 23
	end if
	
	pcum=pcum+rbmin*dt !mm
	
	stmax=0.935+0.498*prstlai-0.00575*prstlai*prstlai !mm
	
	cint=stmax*(1-exp(-(1-p)*(pcum/stmax)))  !mm
	inrate=((cint-lstcint)/dt)*3600
   23	rain=rate-inrate
		
	lstcint=cint
	
	return
	end
c---------------------------------------------------------------------------------------------
c---------------------------------------------------------------------------------------------
	function gamptun3(s,r,ke,dt,rt,t,m,n,n2,idate
     .					,cu,df,thets,theti,dx2,psimat)
c	 Chu, WRR 1978 
	
      implicit real (a-h,o-z)
	real s,r,ke, thets,theti,psimat,fmax

      real(8) ns,lf,tex(100485),oldtex(100485),oldp(100485)
	real(8)  cumf(100485), oldcumf(100485),oldt(100485),tp(100485)
	real(8)  ts(100485),tf(100485),ptp(100485),pr(100485),kf
	real(8)  excess(100485),rbithr,lastf(100485)
	integer pond(100485),m,tests(100485)
	
	ns=(thets-theti)*psimat
	date=200.
	n3=n
	eps=1.0e-6
	rbit=r/60.
	rbithr=r
	kf=ke/60.
	st=s
	if (t.eq.0.or.t.eq.0.5) then
	pond(m)=0.
	excess(m)=0.
	oldp(m)=0.
	oldtex(m)=0.
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

	if (t.eq.0.) goto 2		!never happend
	endif
	
	pr(m) = rbit*(dt/60.) + pr(m)

	if (rbit.gt.0.) then
	  if (rbit.gt.kf) then
		tests(m) = 0.
		if (pond(m).eq.0.) then
	
			cun = pr(m) - oldtex(m) - (kf*ns/(rbit-kf))
			if (cun.le.0.) then
				tex(m)=oldtex(m)
				cumf(m) = pr(m)-oldtex(m)
          		!cumf(m)=oldcumf(m) + rbithr*dt/3600.
				fmax=ke*((ns/cumf(m))+1)
				pond(m)=0.
				
				goto 2

			else
	          
			   tp(m) = (kf*ns/(rbit-kf)-oldp(m)+oldtex(m))/rbit + 
     .			       oldt(m)
			   ptp(m)= oldp(m) + (tp(m) - oldt(m))*rbit
			   ts(m) = ((ptp(m) - oldtex(m))/ns - log(1 + (ptp(m) -
     .				   oldtex(m))/ns))*ns/kf

			   tf(m) =  t - tp(m) + ts(m)

			   if (tf(m).lt.0) then
					tf(m)=0
			   end if			   

			   pond(m) = 1.
			   
			   	if (cumf(m).gt.0.) then
				cinf = cumf(m)
				else
				cinf = 0.1
				end if

				fcerr= cinf - ns*log(1 + cinf/ns) - kf*tf(m)
				dfcerr= cinf/(cinf+ns)
	
				do while (abs(fcerr).gt.eps)
				fcerr=cinf - ns*log(1 + cinf/ns) - kf*tf(m)
				dfcerr=cinf/(cinf+ns)
				cinf=cinf-(fcerr/dfcerr)
				end do
				cumf(m)=cinf
				tex(m) = pr(m) - cumf(m)
				fmax=kf*60.*(1 + ns/cumf(m))
				if (rbit.gt.kf.and.m.eq.n3) then
				excess(m) = excess(m) + ((rbithr-(fmax))/60.)
     .						*(dt/60.)
				end if
				
				goto 2

			end if
		
		else	

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
		
			do while (abs(fcerr).gt.eps)
			fcerr=cinf - ns*log(1 + cinf/ns) - kf*tf(m)
			dfcerr=cinf/(cinf+ns)
			cinf=cinf-(fcerr/dfcerr)
	
			end do
			cumf(m)=cinf
			
			cp = pr(m) - cumf(m) - oldtex(m)
			
			if (cp.gt.0.) then
			
				pond(m) = 1.
				tex(m) = pr(m) - cumf(m)
				fmax=kf*60.*(1 + ns/cumf(m))
				if (rbit.gt.kf.and.m.eq.n3) then
					excess(m) = excess(m) + ((rbithr-(fmax))/60.)
     .							*(dt/60.)
				end if
			   
                 goto 2

			else

				tex(m)=oldtex(m)
          		cumf(m)=pr(m)-oldtex(m)
				!cumf(m)=oldcumf(m) + rbithr*dt/3600.
				fmax=ke*((ns/cumf(m))+1)
				pond(m)=0.
				tests(m) =0.
				goto 2
			   
		    end if

		end if

        else

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
		
			do while (abs(fcerr).gt.eps)
			fcerr=cinf - ns*log(1 + cinf/ns) - kf*tf(m)
			dfcerr=cinf/(cinf+ns)
			cinf=cinf-(fcerr/dfcerr)
	
			end do
			cumf(m)=cinf
			cp = pr(m) - cumf(m) - oldtex(m)
			if (cp.gt.0.) then		
					tests(m) = 1.
					pond(m) = 1.
					tex(m) = pr(m) - cumf(m)
					fmax=kf*60.*(1 + ns/cumf(m))
					if (rbit.gt.kf.and.m.eq.n3) then
					!if (m.eq.n3) then
					excess(m) = excess(m) + ((rbithr-(fmax))/60.)
     .							*(dt/60.)

					end if
				go to 2
			else
				tex(m)=oldtex(m)
          		cumf(m)=pr(m)-oldtex(m)
				fmax=ke*((ns/cumf(m))+1)
				pond(m)=0.
				tests(m) =0.
				goto 2 		 
			end if
		else
						
			tex(m)=oldtex(m)
			cumf(m)=pr(m)-oldtex(m)
			fmax=ke*((ns/cumf(m))+1)
			pond(m)=0.
			tests(m) = 0.
			
			goto 2
	    end if

	  end if

	else

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
	
	
			do while (abs(fcerr).gt.eps)
			fcerr=cinf - ns*log(1 + cinf/ns) - kf*tf(m)
			dfcerr=cinf/(cinf+ns)
			cinf=cinf-(fcerr/dfcerr)
	
			end do
			cumf(m)=cinf
			cp = pr(m) - cumf(m) - oldtex(m)
			if (cp.gt.0.) then		
					tests(m) = 1.
					pond(m) = 1.
					tex(m) = pr(m) - cumf(m)
					fmax=kf*60.*(1 + ns/cumf(m))
					if (rbit.gt.kf.and.m.eq.n3) then

					excess(m) = excess(m) + ((rbithr-(fmax))/60.)
     .							*(dt/60.)

					end if
				go to 2
			else
				
				tex(m)=oldtex(m)
          		cumf(m)=pr(m)-oldtex(m)
				fmax=ke*((ns/cumf(m))+1)
				pond(m)=0.
				tests(m) =0.
				goto 2 		 
			end if
	  else
			
			
			tex(m)=oldtex(m)
			cumf(m)=pr(m)-oldtex(m)
			fmax=ke*((ns/cumf(m))+1)
			pond(m)=0.
			tests(m) = 0.
			
			goto 2
	    
	  end if
	end if
	

    2 gamptun3=fmax
      
      oldtex(m)=tex(m)
	oldcumf(m)=cumf(m)
	oldp(m)=pr(m)
	oldt(m)=t
      lastf(m)=fmax
	
   12 format (1x,i3,1x,i3,1x,f7.2,1x,f8.4,1x,f8.4,1x,f8.4,1x,f8.4,1x,i2
     .,1x,i2,1x,f8.4,1x,f8.4,1x)
	
	return

	end
	
c------------------------------------------------------------------------------------------------	
c------------------------------------------------------------------------------------------------	       
	function drain(tetinst,ks,thets,fc,df,dt,ks2,idate)
	
	implicit real (a-z)
	integer idate	

	dthr=dt/3600
 
 10   if (tetinst.gt.fc) goto 20
      dr=0.
	drain=dr
      return
 
********************************
 20   continue
	rkfc=ks2*(tetinst/thets)**(-2.655/log10(fc/thets))
	ti=((tetinst-fc)*thets*df)/rkfc

	if (ti.lt.(dt/3600.)) ti=dt/3600.

	dr=df*(tetinst-fc)*(1-exp(-(dthr/ti))) ! dr in mm/h
	
	drain=dr
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
      data tol, rtol /1.E-7,1.e-4/

      xlow = xmin
      xhigh = xmax
      ierr = 1
      call errfct (xlow, fxlow, deriv)

      call errfct (xhigh, fxhigh, deriv)

      if(xlow. gt. xhigh .and. fxlow .gt. fxhigh) then
        write(*, 977) xlow, xhigh
	  read(*,*)
        stop ' Limits contradiction in ITER'
      end if
  977  format(' Lower bound higher than upper bound in call to ITER'/
     &       '  xlow = ',g13.4,',  xhigh = ',g13.4)
      x1 = xvar

      call errfct (x1, fx, deriv)
C                                                                iteration loop
C------------------------------------------------------------------------------

      do i = 1, 200

        if (deriv .ne. 0. .and. i .le. 100) then
C                                                       compute newton estimate
          x2 = x1 - fx / deriv

          if (x2 .le. xlow .or. x2 .ge. xhigh) then
C                                                                        bisect
            x2 = (xlow + xhigh) / 2.

          end if

        else
C                                                                        bisect
          x2 = (xlow + xhigh) / 2.

        end if
C                                            get function value of new estimate
        ofx = fx
  10    call errfct (x2, fx, deriv)
C
  900 format(' itr:', i2,3g12.4)

        test1 = abs(fx)
C                                             scale the test criteria if x2 < 1
        if (abs(x2) .gt. 100485.*tol)   test1 = abs (fx / x2)
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
          ierr = -i
          return

        end if
C                                                              bracket the root
        if (fx * fxhigh .gt. 0.) then

          xhigh = x2

        else

          xlow = x2

        end if

        x1 = x2

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
	ierr=1
	return 
	END 

c-----------------------------------------------------------------------------------------
c-----------------------------------------------------------------------------------------
	subroutine errq2 (q22, ferr, dferr)
	
C   Computes the residual ferr and derivative dferr at q22 of the error
C   function form of the finite difference approximation of the kinematic
C   equation for overland and channel flow.

      integer iw
	real qp12,qp21,qp11,beta,omega,qq
	real alp11,alp12,alp21,alp22
	real qin21,qin22,qin12,qin11,dt,dx
	common /kin3/ alp11(100485),alp12(100485),alp21(100485),alp22(100485)
      common /kin/ qin21,qin22,qin12,qin11,beta,omega
	common /kin2/ qp21(100485),qp12(100485),qp11(100485),iw
	common /gen2/ dt,dx


	ferr = alp22(iw)*(q22**beta)+alp12(iw)*(qp12(iw)**beta)-alp21(iw)*
     .	   (qp21(iw)**beta)-alp11(iw)*(qp11(iw)**beta)+(2*dt/dx)*
     .	   (omega*(q22-qp12(iw))+(1-omega)*(qp21(iw)-qp11(iw)))-dt*
     .	   (qin12+qin22)

	if (q22.ne.0.) then
	dferr = alp22(iw)*beta*q22**(beta-1) + (2*dt/dx)*omega
	else
	dferr = (2*dt/dx)*omega
	end if
	
      end
c-------------------------------------------------------------------------------------------------
c-------------------------------------------------------------------------------------------------
	subroutine errq3 (q22, ferr, dferr)
	
C   Computes the residual ferr and derivative dferr at q22 of the error
C   function form of the finite difference approximation of the kinematic
C   equation for overland and channel flow.

      integer iw
	real qp12,qp21,qp11,alp,beta,omega,qq,rterm
	real qin21,qin22,DT,DX
      common /kin/ qin21,qin22,alp,beta,omega
	common /kin2/ qp21(100485),qp12(100485),qp11(100485),iw
	common /gen2/ dt,dx

      
	ferr = alp*(q22**beta)+alp*(qp12(iw)**beta)-alp*(qp21(iw)**beta)-
     .	   alp*(qp11(iw)**beta)+(2*DT/DX)*(omega*(q22-qp12(iw))+
     .	   (1-omega)*(qp21(iw)-qp11(iw)))-DT*(qin21+qin22)

	dferr = alp*beta*q22**(beta-1) + (2*DT/DX)*omega
	
      end
c-------------------------------------------------------------------------------------------------
c-------------------------------------------------------------------------------------------------
	subroutine errhwp (hwp, ferr, dferr)
	integer iw2
	real qj,mann,slwp,width2,qact
	common /hwpr/ slwp,qact,width2(100485),iw2
      common /manning/ mann(100485)
	
	qj= ((slwp**0.5)/mann(iw2))*(((width2(iw2)*hwp)**1.666667)/
     .	((width2(iw2)+2*hwp)**0.666667))

	ferr = 1 - (qact/qj)

	dferr = (5*width2(iw2) + 6*hwp)/(3*hwp*(width2(iw2)+2*hwp))
	
      end

c-----------------------------------------------------------------------------------------------
c-------------------------------------------------------------------------------------------------
	subroutine etp (tempc,radi,albed,etpmm) 
      implicit real (a-z)
      integer tempc,radi     
	 
       tempk = tempc + 273
       tempr = 5304./tempk
 
c	etp Ritchie's equation
 
       delta = (tempr/tempk)*exp(21.25-tempr)
	
	 h0 = (1-albed)*radi*1./548.11	 

       etpm = 0.00211*h0*delta/(0.68+delta)
       etpmm = etpm*1000.	
 
      return
      end

c--------------------------------------------------------------------------------------------
c-------------------------------------------------------------------------------------------------	
	subroutine perco(theti,ks,thets,fc,df,dt,ks2,idate,percol,
     .lsttr)
	
	implicit real (a-z)
	integer idate	

	dthr=dt/3600.
	
	if (theti.ge.thets) theti=thets
 
********************************

	rkfc=ks2*((theti)/thets)**(-2.655/log10(fc/thets))
	ti=(((theti)-fc)*thets*df)/rkfc

	percol=df*((theti)-fc)*(1-exp(-(dthr/ti)))
	
	return
      end
	
c------------------------------------------------------------------------------------------
c-------------------------------------------------------------------------------------------------
	subroutine etpr2 (prstlai,totprec,idate)
      implicit real (a-z)
	integer i,ie,idate,ttime(100485)
	common /etpot/ etpmm
	common / etpact / esu(100485),edx(100485),es1p(100485),es2p(100485),
     .				  es(100485),pep(100485),cuminf(100485)
     .					,alphas(100485),edispm,edispmfc,ie

	i=ie

	esm = etpmm*exp(-0.4*prstlai)

c	soil evaporation

	if (idate.eq.1) then
	
		es1p(i) = 0.		! es1p mm accumulated since first stage beginning
		es2p(i) = 0.        ! es2p mm accumulated since second stage beginning
		ttime(i)=0

	end if

	if (edispm.ge.edispmfc.and.totprec.gt.0.) then !new
	
		es1p(i) = 0.
		es2p(i) = 0.
		ttime(i)=0
		
		esm = etpmm*exp(-0.4*prstlai)

		es1p(i) = es1p(i) + esm	

		goto 31

	elseif (-1.gt.alphas(i)) then

c	si on est toujours dans first stage

		if(es2p(i).eq.0.) then


			esm = etpmm*exp(-0.4*prstlai)

			es1p(i) = es1p(i) + esm

			goto 31

		else

c	si on �tait dans second stage, on le recommence au d�but

			es2p(i) = 0.
			ttime(i)=0

			goto 32

		end if

	else

		if(es2p(i).eq.0.) then

			esm = etpmm*exp(-0.4*prstlai)

			es1p(i) = es1p(i) + esm

			goto 31

		else

			goto 32

		end if

	end if

   31	if (es1p(i).gt.esu(i)) then

		goto 32

	else	


		es(i) = max1(0.0,esm)
		es2p(i) = 0.
		ttime(i) = 0

		goto 33

	end if


   32 espot = max1(0.0,esu(i))

	ttime(i)=ttime(i) + 1


	es2=alphas(i)*(ttime(i)**(0.5))-alphas(i)*((ttime(i)-1.)**(0.5))

	es(i) = min1(esm,es2)
		
	es2p(i) = es2p(i) + es(i) 

	goto 33

c	plant transpiration

   33	if (prstlai.le.3.) then

		pep(i) = (etpmm-es(i))*prstlai/3.
	else

		pep(i) = (etpmm-es(i))

	end if
	if (pep(i).gt.etpmm) pep(i)=etpmm
	
	return 
	end

c------------------------------------------------------------------------------------------
c------------------------------------------------------------------------------------------
	function Lerosionsplash(dt,dx2,hsvmm, ek,ekt,p,pt
     .	,As, prstcover)
      implicit real (a-z)
	
c	Ds = splash detachment gs-1 occurs only if ponded depth h > surface storage
c	Ds_r = rainfall, Ds_t = troughfall	
	Ar=dx2*(1-prstcover)
	At=dx2*(prstcover)
	
	Ds_r=1e-3*((2.82/As)*ek*exp(-1.48*hsvmm)+2.96)*p*Ar/3600	!kg/s

c	Ds_r=1e-3*((2.82/As)*ek*dexp(-1.48*hsvmm)+2.96)*(p*dt/3600)*Ar/dt	!kg/s

	Ds_t=1e-3*((2.82/As)*ekt*exp(-1.48*hsvmm)+2.96)*pt*At/3600	!kg/s

	Ds=Ds_r +Ds_t
	Lerosionsplash=Ds
	return 
	end
c------------------------------------------------------------------------------------------
