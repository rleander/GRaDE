	program seco03	! first and second order stats of a group of simulated record compaired to their historical records 
C	use msflib	! This line ony applies to Digital Visual Fortran
	implicit none 
	! calculates second-order stats for station data 
	! Over all runs 
	!    Over all datyrs years in a run 
	!	Over all months m 
	!	   chain mth month of each year in a run to single array 
	!	   jackknife sub returns nlag autocor.coefs estimates and 
	!	         standard error in these estimates 
	! 		 and the sample variance 
	!	average the results over all the months (in the halfyear)
	!	se in this average from standard errors for each month 
	!    se in this average from se in each run 
	! end result: nlag autocorrelation estimates and their standard errors 
	! over the year or the half year. 

	integer maxblock 		! maximum number of years in a block 
	parameter(maxblock=2000)

	integer datum(maxblock*365)	! stored date corresponding to each day in the current block 
	real values(maxblock*365)	! stored value for each day in the current block 

	integer nlag,nm,nblock
	real numRR			! sum of weights of precipitation (RR) files 
	real numTG			! sum of weights of temperature (TG) files 
	real weight			! dummy weight read from file list 
	integer dayndx 			
	parameter(nlag=20)
	real msum(2)			! monthly sum/average
	real msumsqr(2)			! and stddev of monthly sum

!	character*(35) fmt1,fmt2
	character*(15) my_name  

	!-----New Variables 
	real value 			! value read from series 
	integer date 			! date read from series 
	integer blocksize 		! number of years in block 
	character*(12) strmsel		! string read from arg[3] containing one char for each of the 12 months (optional)
					! if i-th month NOT selected the ith position of this string should 
					! contain a '.'  

	character*(30) line_in 		! single line read from file 
	integer TYPE_OF_VAR		! temperature of rainfall 
	!-----New Variables 

	!---------------------------------------------------------------------------------
					! jackknife arrays 
	real Xm(maxblock)		! month sums/averages 
	real Xd(maxblock*31)		! day sums/averages
	real theta(maxblock)		! pseudo-values dummy var 
	real sdX

	real logpstd_m(maxblock)		! log pseudo monthsum stddev averaged over months 
	real logpstd_d(maxblock)		! log pseudo daysum stddev averaged over months 
	real logpavg_m(maxblock)		! log pseudo monthsum mean averaged over months 

	real pseu(maxblock*nlag)

				! RESULTS FROM A SINGLE BLOCK OF DATA
				! 2-component vectors: the firstelement is supposed to be the estimate 
				! and the second the standard error of this estimate 
	real avg_m(2)		! estimated mean monthly		biased estimates 
	real std_m(2)		! estimated stddev monthly
	real std_d(2)		! estimated stddev daily
	real theta_m(2)		! estimated theta monthly
	real theta_d(2)		! estimated theta daily
	real acf(2,nlag)	! estimated autocorrelation 

	real pstd_m(maxblock)	! pseudo monthsum stddev 		pseudovalues 
	real pstd_d(maxblock)	! pseudo daysum stddev 
	real pavg_m(maxblock)	! pseudo monthsum mean 
	real pacf(maxblock*nlag)! pseudo autocorrelation 

	real ptheta_m(maxblock)	! pseudo monthsum theta 		theta is the log of the stdev  
	real ptheta_d(maxblock)	! pseudo daysum theta 

					! RESULTS FROM THE HISTORICAL BLOCK OF DATA
	real hisavg_m(2)		! estimated mean monthly		statistics estimates 
	real hisstd_m(2)		! estimated stddev monthly
	real hisstd_d(2)		! estimated stddev daily
	real hisacf(2,nlag)		! estimated autocorrelation 
	real histheta_m(2)		! estimated theta for monthly amounts 
	real histheta_d(2)		! estimated theta for daily amounts 

	real hispstd_m(maxblock)	! pseudo monthsum stddev 		pseudovalues 
	real hispstd_d(maxblock)	! pseudo daysum stddev 
	real hispavg_m(maxblock)	! pseudo monthsum mean 
	real hispacf(maxblock*nlag)	! pseudo autocorrelation 

					! RESULTS FROM THE SIMULATED BLOCKs OF DATA
	real simavg_m(2)		! estimated mean monthly		statistics estimates 
	real simstd_m(2)		! estimated stddev monthly
	real simstd_d(2)		! estimated stddev daily
	real simacf(2,nlag)		! estimated autocorrelation 
	real simpstd_m(maxblock)	! pseudo monthsum stddev 		pseudovalues 
	real simpstd_d(maxblock)	! pseudo daysum stddev 
	real simpavg_m(maxblock)	! pseudo monthsum mean 
	real simpacf(maxblock*nlag)	! pseudo autocorrelation 

	real biasavg_m			! sim - his for all stats 
	real biasstd_m
	real biasstd_d
	real biasstd_m_rel
	real biasstd_d_rel
	real biasacf(nlag)

!					! variables refreshed for every new sim-his combination 
	!---------------------------------------------------------------------------------
	real RRbiasavg_m 		! bias averaged over all RR files 
	real RRbiasstd_m
	real RRbiasstd_d
	real RRbiasstd_m_rel
	real RRbiasstd_d_rel
	real RRbiasacf(nlag)

	real TGbiasavg_m		! bias averaged over all TG files 
	real TGbiasstd_m
	real TGbiasstd_d
	real TGbiasstd_m_rel
	real TGbiasstd_d_rel
	real TGbiasacf(nlag)

	real RRhispstd_m(maxblock)	! pseudo monthsum stddev averaged over RR files 
	real RRhispstd_d(maxblock)	! pseudo daysum stddev 
	real RRhispavg_m(maxblock)	! pseudo monthsum mean 
	real RRhisptheta_m(maxblock)	! pseudo monthsum theta ( = log[stddev])
	real RRhisptheta_d(maxblock)	! pseudo daysum theta 
	real RRhispacf(maxblock*nlag)	! pseudo autocorrelation 

	real TGhispstd_m(maxblock)	! pseudo monthsum stddev averaged over TG files 
	real TGhispstd_d(maxblock)	! pseudo daysum stddev 
	real TGhispavg_m(maxblock)	! pseudo monthsum mean 
	real TGhisptheta_m(maxblock)	! pseudo monthsum theta 
	real TGhisptheta_d(maxblock)	! pseudo daysum theta 
	real TGhispacf(maxblock*nlag)	! pseudo autocorrelation 


	real RRavg_m(2),RRstd_m(2),RRstd_d(2),RRacf(2,nlag)	! RR-average stats from pseudo values 
	real RRtheta_m(2),RRtheta_d(2)
	real TGavg_m(2),TGstd_m(2),TGstd_d(2),TGacf(2,nlag)	! TG-average stats from pseudo values 
	real TGtheta_m(2),TGtheta_d(2)
	!---------------------------------------------------------------------------------


					! read from historic data 
	
	real R_j(nlag),r_e(nlag),se_j(nlag)
	real sampvar
	integer i,iy
	integer im,ilag,jlag,id 
	real cvn
	integer mselect(12)		! 1 for month taken into account 
					! 0 for month ignored 
	character*(36)	Mnames		! monthnames abbreviated 
	data Mnames/'JanFebMarAprMayJunJulAugSepOctNovDec'/  
!	data mselect /1,1,1,0,0,0,0,0,0,1,1,1/
	character*(*) CAPITALS		! contains all the Capitals in the alphabet 
	parameter(CAPITALS='ABCDEFGHIJKLMNOPQRSTUVWXYZ')
	integer nd(12) 			! days in each month
	integer dayzero(12) 		! day b4 each month starts
	character*(80) regel 
	character*(50) fhis,fsim
	integer instr 			! function determines position of substring in string 


	data dayzero /0,31,59,90,120,151,181,212,243,273,303,334/	! 0th day of each month 
	data nd /31,28,31,30,31,30,31,31,30,31,30,31/			! number of days in each month 

	call getarg(0,my_name)		! program name 
!	if(nargs().lt.3) then 
!	   write(*,*) my_name(1:instr(my_name,' ')-1) 
!     &	               //' historical_record  simulated_record'
!     &	               //' [msel]' 
!	   write(*,*) '    calculates second order statistics '
!     	   write(*,*) '    (mean and stdd. monthsum and day stdd and'
!     &	                //' autocorrelation)'
!	   write(*,*) ' file given in arg[2] '
!           write(*,*) ' compared to the historical record'
!     &	                //', given in arg[1]' 
!	   write(*,*) ' msel is a 12-string with a CAPTITAL for each' 
!     &                  // ' month IN the season' 
!	   write(*,*) ' eg : JFMamjjasOND or XXXrrrrrrFFF selects the' 
!     &	                // ' winter-months (default)' 
!	   stop '#######'
!	endif 

	if(iargc().ge.1) then 		! interprete arg[1] as a string with selected months 
	   call getarg(1,strmsel) 	! string of 12 positions containing a capital for selected month
	else 				! and anything else for UNselected month
	   strmsel='JFMamjjasOND'	! winter season selected 
	endif 
	do im=1,12
	   if(instr(CAPITALS,strmsel(im:im)).le.0) then 
	      mselect(im)=0
	   else 
	      mselect(im)=1		! selected months have mselect > 0
	   endif 
	enddo ! im 
	nm=0
	do im=1,12			! count selected months 
	   if(mselect(im).gt.0) then 
	      nm=nm+1
	   endif 
	enddo 

	write(*,'(a$)') '# '
	do im=1,12
	   if(mselect(im).gt.0) then 	! this month is selected 
	     write(*,'(a$)') Mnames((im-1)*3+1:im*3)//' ' 
	   endif 
	enddo 
	write(*,*) 

	  
!+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_
! initialisation of all-file averages of biasses and pseudo values 
	  
	RRbiasavg_m = 0.0 		! bias averaged over all RR files 
	RRbiasstd_m = 0.0
	RRbiasstd_d = 0.0
	RRbiasstd_m_rel = 0.0
	RRbiasstd_d_rel = 0.0
	do ilag=1,nlag
	   RRbiasacf(nlag) = 0.0
	enddo 

	TGbiasavg_m = 0.0		! bias averaged over all TG files 
	TGbiasstd_m = 0.0
	TGbiasstd_d = 0.0
	TGbiasstd_m_rel = 0.0
	TGbiasstd_d_rel = 0.0
	do ilag=1,nlag
	   TGbiasacf(nlag) = 0.0
	enddo 

	do iy=1,maxblock
	   RRhispstd_m(iy) = 0.0	! pseudo monthsum stddev averaged over RR files 
	   RRhispstd_d(iy) = 0.0	! pseudo daysum stddev 
	   RRhispavg_m(iy) = 0.0	! pseudo monthsum mean 

	   TGhispstd_m(iy) = 0.0	! pseudo monthsum stddev averaged over TG files 
	   TGhispstd_d(iy) = 0.0	! pseudo daysum stddev 
	   TGhispavg_m(iy) = 0.0	! pseudo monthsum mean 
	enddo ! iy 

	do iy=1,maxblock*nlag
	   RRhispacf(iy) = 0.0	! pseudo autocorrelation 
	   TGhispacf(iy) = 0.0	! pseudo autocorrelation 
	enddo ! iy 

	numRR=0.0		! total of weights for the precipitation 
	numTG=0.0		! total of weights for the temperature 
!+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_



!=========================== LOOP OVER ALL FILES READ FROM STDIN =====================------
	DO  

 888	read(*,'(a80)',end=188) regel
	if(instr(regel,'#').gt.0) goto 888	! commented line 
	weight=1.0 
 	read(regel,*) 	fsim,fhis,weight	! read filename of historical series and simulated series from STDIN 
 911	continue 				! ignore read errors, if reading weight raises error flag, 
						! weight automatically retains its former value of 1.0 

	! open historical series 
	open(13
     &	   ,file=fhis(1:instr(fhis,' ')-1)	! open file as #13 for reading 
     &	   ,status='OLD'			! ascii read 
     &	   ,err=213)
	goto 113
 213	continue 
	write(0,*) fhis(1:instr(fhis,' ')-1)//' non-existent'
	stop 
 113	continue
	write(0,*)'#',fhis(1:instr(fhis,' ')-1)
     &	       //' opened as historical record, weight :',weight  
	if(instr(fhis,'.tg').gt.0) then 
	   TYPE_OF_VAR = 1 			! temperature data 
	else 
	   TYPE_OF_VAR = 0 			! rainfall (?) data 
	endif 

!---------------------historical file openend for reading-------------------------------

	   id=1
	   do 
							! assume each year in the block to have 365 days
 613 	     read(13,'(a30)',end=513) line_in 	! read line from ascii file 
	     if(instr(line_in,'#').gt.0) goto 613
	     read(line_in,*) date, value
	     if(mod(date,10000).eq.0229) goto 613	! skip leap days   		
	     datum(id)=date
	     values(id)=value 
	     id=id+1
	   enddo  ! days in file  
 513	   continue 					
	   close(13)

	   blocksize=((id-1)/365)		! block size equals the length of the 
						! historical record (years)
	   write(0,*) '# blocksize = ', blocksize,' years'

 	   avg_m(1)=0.0				
	   std_m(1)=0.0 
	   std_d(1)=0.0
	   do iy=1,maxblock	
	      pavg_m(iy)=0.0 
	      pstd_m(iy)=0.0 
	      pstd_d(iy)=0.0
	   enddo 

	   do ilag=1,nlag
	      acf(1,ilag)=0.0	
	   enddo 
	   do iy=1,maxblock*nlag			
	      pacf(iy)=0.0	
	   enddo 


	   msum(1)=0.0
	   msumsqr(1)=0.0

	   do im=1,12				! Month Loop
	      if(mselect(im).gt.0) then 	
	         msum(2)=0.0			! average monthsum for this calendar month 
	         msumsqr(2)=0.0			! average squared monthsum for this calendar month 
	         do iy=1,blocksize
	            Xm(iy)=0.0
	            do id=1,nd(im) 	
	               dayndx = (iy-1)*365+dayzero(im)+id 	! refer to day id in month im in the current block 	
	               Xm(iy)=Xm(iy)+values(dayndx)
	            enddo ! id 
	            if(TYPE_OF_VAR.eq.1) then 	! TYPE_OF_VAR = 1 for temperatures 
	               Xm(iy)=Xm(iy)/nd(im)		!               and 0 for rainfall 
	            endif 
	               					! Xm(iy) contains the monthsum for the current month 
							! in the year iy
	            msum(2)=msum(2)+Xm(iy)		! Xm summed over years 
	         enddo !iy
	         call sdIjack(Xm,blocksize,1,theta,sdX) 	! jack monthsum statistics 
	         avg_m(1)=avg_m(1)+msum(2)/blocksize/nm		! monthsum mean averaged over months 
	         std_m(1)=std_m(1)+sdX/nm				! monthsum stddev averaged over months 
	         do iy=1,blocksize				! averaging pseudo values of mean and stddev 
	            pavg_m(iy)=pavg_m(iy)+(msum(2)-Xm(iy))	! pseudo array average monthsum 
     &	                                 /(blocksize-1)/nm
	            pstd_m(iy)=pstd_m(iy)+theta(iy)/nm		! pseudo array stdev monthsum 
	            do id=1,nd(im)				! chain months to calculate autocorrelation
	               dayndx = (iy-1)*365+dayzero(im)+id 
	               Xd(id+(iy-1)*nd(im))=values(dayndx)	
	            enddo ! days in month 
	         enddo ! iy

      	         call acjack(Xd, blocksize, nd(im), nlag, r_e	! acjack sees the array pseu as 
     &	                       ,R_j, se_j, pseu, sampvar)	! a 2D (blocksize,nlag) array 
								! but it is actually a 1D array of 
								! length maxblock*nlag
 
	         do ilag=1,nlag
	            acf(1,ilag)=acf(1,ilag)+R_j(ilag)/nm		! average R_j over months 
	            do iy=1,blocksize
	               pacf((ilag-1)*blocksize+iy)
     &	                  = pacf((ilag-1)*blocksize+iy)	
     &	                  + pseu((ilag-1)*blocksize+iy)/nm	! average pseudo values over months 
	            enddo 
	         enddo 

	         call sdPjack(Xd,blocksize,nd(im),theta,sdX) 	! jack DAYsum statistics 
	         std_d(1)=std_d(1)+sdX/nm				! estimate of day stddev averaged over months 
	         do iy=1,blocksize
	            pstd_d(iy)=pstd_d(iy)+theta(iy)/nm		! pseudo values of the day stddev 
	         enddo 
	      endif ! selected months 
	   enddo ! im 		


	! --------------------------------------------------------------------
	! Output of this loop for a single block 
	!
	!	avg_m		sample monthsum average 
	!	std_m		sample monthsum stddev
	!	std_d		sample daily stddev 
	!
	!	pavg_m		pseudo monthsum average 
	!	pstd_m				stddev
	!	pstd_d			daily stddev
	!
	!	acf		sample autocorrelation coefficients 
	!	pacf		pseudovalues of autocorrelation coefficients 
	! --------------------------------------------------------------------

!	   write(*,'(a,3f8.2,3f8.3)')'# his', avg_m(1), std_m(1), std_d(1) 	! show results for current var 
!     &	       ,acf(1,1),acf(1,2),(acf(1,3)+acf(1,4)+acf(1,5))/3.	


	! jackknife estimates of statistics 
	! blocksize Pseudo-values are calculated here for each statistic 
	do iy=1,blocksize						! bias correction 
	   ptheta_m(iy)=log(std_m(1))					! pseudo values of theta 
     &	        +(blocksize-1)*(log(std_m(1))-log(pstd_m(iy)))		! ...on monthly basis
	   ptheta_d(iy)=log(std_d(1))					! or 
     &	        +(blocksize-1)*(log(std_d(1))-log(pstd_d(iy)))		! ...on daily basis
	   pavg_m(iy)=avg_m(1)+(blocksize-1)*(avg_m(1)-pavg_m(iy))	! pseudo monthly mean 
	   pstd_m(iy)=std_m(1)+(blocksize-1)*(std_m(1)-pstd_m(iy))	! pseudo monthly stddev
	   pstd_d(iy)=std_d(1)+(blocksize-1)*(std_d(1)-pstd_d(iy))	! pseudo daily stddev 
	enddo 
	! note : ACJACK returns the pseudo values themselves, 
	!        no tranformation is required

	
	call pseu2est(pavg_m,blocksize,avg_m)
	call pseu2est(pstd_m,blocksize,std_m)
	call pseu2est(pstd_d,blocksize,std_d)
	call pseu2est(ptheta_m,blocksize,theta_m)
	call pseu2est(ptheta_d,blocksize,theta_d)


	do ilag=1,nlag
	   call pseu2est(pacf((ilag-1)*blocksize+1),
     &	        blocksize,acf(1,ilag))
	enddo ! ilag 
!	   write(*,'(a,3f8.2,3f8.3,a)')'# his'
!     &	       ,avg_m(1), std_m(1), std_d(1) 	! show results for current var 
!     &	       ,acf(1,1),acf(1,2),(acf(1,3)+acf(1,4)+acf(1,5))/3.	
!     &	       ,' jackknife' 



C======================================================================================================
C Store historical stats 

	do i=1,2 
	   hisavg_m(i)=avg_m(i)			
	   hisstd_m(i)=std_m(i)
	   hisstd_d(i)=std_d(i)
	   histheta_m(i)=theta_m(i)
	   histheta_d(i)=theta_d(i)
	   do ilag=1,nlag
	      hisacf(i,ilag)=acf(i,ilag)
	   enddo
	enddo ! i 
	
	do iy=1,blocksize
	   hispavg_m(iy)=pavg_m(iy)
	   hispstd_m(iy)=pstd_m(iy)
	   hispstd_d(iy)=pstd_m(iy)	
	enddo 
	do ilag=1,nlag*blocksize
	   hispacf(ilag)=pacf(ilag)
	enddo


C NB : the his... variables are redundant, but are maintained such to keep 
C this program comparable to seco02.f
C Add pseudo-values to the TG en RR averages 
	
	if(TYPE_OF_VAR.eq.0) then 			! precipitation 
	   do iy=1,blocksize
	      RRhispavg_m(iy)=RRhispavg_m(iy)+hispavg_m(iy)
	      RRhispstd_m(iy)=RRhispstd_m(iy)+hispstd_m(iy)
	      RRhispstd_d(iy)=RRhispstd_d(iy)+hispstd_d(iy)
	   enddo 
	   do ilag=1,nlag*blocksize
	      RRhispacf(ilag)=RRhispacf(ilag)+hispacf(ilag)
	   enddo
	endif 

	
	if(TYPE_OF_VAR.eq.1) then 			! temperature
	   do iy=1,blocksize
	      TGhispavg_m(iy)=TGhispavg_m(iy)+hispavg_m(iy)
	      TGhispstd_m(iy)=TGhispstd_m(iy)+hispstd_m(iy)
	      TGhispstd_d(iy)=TGhispstd_d(iy)+hispstd_d(iy)
	   enddo 
	   do ilag=1,nlag*blocksize
	      TGhispacf(ilag)=TGhispacf(ilag)+hispacf(ilag)
	   enddo
	endif 


C======================================================================================================
C Open simulated file of data 

	! open generated series 
	open(14
     &	   ,file=fsim(1:instr(fsim,' ')-1)	! open file as #13 for reading 
     &	   ,status='OLD'			! ascii read 
     &	   ,err=214)
	goto 114
 214	continue 
	write(0,*) fsim(1:instr(fsim,' ')-1)//' non-existent'
	stop 
 114	continue
	write(*,*)'#',fsim(1:instr(fsim,' ')-1)
     &	       //' opened as simulation record'  

	nblock=0
 	avg_m(1)=0.0			! averages over ALL blocks 	
	std_m(1)=0.0 
	std_d(1)=0.0
	do iy=1,maxblock	
	   pavg_m(iy)=0.0 
	   pstd_m(iy)=0.0 
	   pstd_d(iy)=0.0
	enddo 

 	do ilag=1,nlag
	   acf(1,ilag)=0.0 	
	enddo 
	do iy=1,maxblock*nlag			
	   pacf(iy)=0.0	
	enddo 

	do 						! loop processing subsequent runs in simulated 
	   do id=1,365*blocksize
 614 	     read(14,'(a30)',end=514) line_in 	
	     if(instr(line_in,'#').gt.0) goto 614	! skip comment lines 
	     read(line_in,*) date, value
	     if(mod(date,10000).eq.0229) goto 614	! skip leap days   		
	     datum(id)=date
	     values(id)=value 
	   enddo  ! days in file  
	   nblock=nblock+1				! full block read 
	  
	   !=========================================== processing block 

	   msum(1)=0.0
	   msumsqr(1)=0.0

	   do im=1,12				! Month Loop
	      if(mselect(im).gt.0) then 	
	         msum(2)=0.0			! average monthsum for this calendar month 
	         msumsqr(2)=0.0			! average squared monthsum for this calendar month 
	         do iy=1,blocksize
	            Xm(iy)=0.0
	            do id=1,nd(im) 	
	               dayndx = (iy-1)*365+dayzero(im)+id 	! refer to day id in month im in the current block 	
	               Xm(iy)=Xm(iy)+values(dayndx)
	            enddo ! id 
	            if(TYPE_OF_VAR.eq.1) then 	! TYPE_OF_VAR = 1 for temperatures 
	               Xm(iy)=Xm(iy)/nd(im)		!               and 0 for rainfall 
	            endif 
	               					! Xm(iy) contains the monthsum for the current month 
							! in the year iy
	            msum(2)=msum(2)+Xm(iy)		! Xm summed over years 
	         enddo !iy
	         call sdIjack(Xm,blocksize,1,theta,sdX) 	! jack monthsum statistics 
	         avg_m(1)=avg_m(1)+msum(2)/blocksize/nm		! monthsum mean averaged over months 
	         std_m(1)=std_m(1)+sdX/nm				! monthsum stddev averaged over months 
	         do iy=1,blocksize				! averaging pseudo values of mean and stddev 
	            pavg_m(iy)=pavg_m(iy)+(msum(2)-Xm(iy))	! pseudo array average monthsum 
     &	                                 /(blocksize-1)/nm
	            pstd_m(iy)=pstd_m(iy)+theta(iy)/nm		! pseudo array stdev monthsum 
	            do id=1,nd(im)				! chain months to calculate autocorrelation
	               dayndx = (iy-1)*365+dayzero(im)+id 
	               Xd(id+(iy-1)*nd(im))=values(dayndx)	
	            enddo ! days in month 
	         enddo ! iy

      	         call acjack(Xd, blocksize, nd(im), nlag, r_e	! acjack sees the array pseu as 
     &	                       ,R_j, se_j, pseu, sampvar)	! a 2D (blocksize,nlag) array 
								! but it is actually a 1D array of 
								! length maxblock*nlag
	         do ilag=1,nlag
	            acf(1,ilag)=acf(1,ilag)+R_j(ilag)/nm		! average R_j over months 
	            do iy=1,blocksize
	               pacf((ilag-1)*blocksize+iy)
     &	                  = pacf((ilag-1)*blocksize+iy)	
     &	                  + pseu((ilag-1)*blocksize+iy)/nm	! average pseudo values over months 
	            enddo 
	         enddo 

	         call sdPjack(Xd,blocksize,nd(im),theta,sdX) 	! jack DAYsum statistics 
	         std_d(1)=std_d(1)+sdX/nm				! estimate of day stddev averaged over months 
	         do iy=1,blocksize
	            pstd_d(iy)=pstd_d(iy)+theta(iy)/nm		! pseudo values of the day stddev 
	         enddo 
	      endif ! selected months 
	   enddo ! im 		


	! --------------------------------------------------------------------
	! Output of this loop for a single block 
	!
	!	avg_m		sample monthsum average 
	!	std_m		sample monthsum stddev
	!	std_d		sample daily stddev 
	!
	!	pavg_m		pseudo monthsum average 
	!	pstd_m				stddev
	!	pstd_d			daily stddev
	!
	!	acf		sample autocorrelation coefficients 
	!	pacf		pseudovalues of autocorrelation coefficients 
	! --------------------------------------------------------------------

C	   write(*,'(a2,i3,3f8.2,3f8.3)')'#', nblock
C     &	       , avg_m/nblock, std_m/nblock
C     &	       ,std_d/nblock,acf(1)/nblock,acf(2)/nblock
C     &	       ,(acf(3)+acf(4)+acf(5))/3./nblock	

	   !=========================================== done processing block 
       
	enddo ! runs 
 514	continue 					
	close(14)

	if(nblock.le.0) stop 'no complete runs found in simulation'
 	avg_m(1)=avg_m(1)/nblock
	std_m(1)=std_m(1)/nblock
	std_d(1)=std_d(1)/nblock
	do iy=1,blocksize	
	   pavg_m(iy)=pavg_m(iy)/nblock
	   pstd_m(iy)=pstd_m(iy)/nblock
	   pstd_d(iy)=pstd_d(iy)/nblock
	enddo 

 	do ilag=1,nlag
	   acf(1,ilag)=acf(1,ilag)/nblock
	enddo 
	do iy=1,maxblock*nlag			
	   pacf(iy)=pacf(iy)/nblock
	enddo 


!	   write(*,'(a,3f8.2,3f8.3)')'# sim'
!     &	        , avg_m(1), std_m(1), std_d(1) 			! show results for current var 
!     &	        , acf(1,1)
!     &	        , acf(1,2),(acf(1,3)+acf(1,4)+acf(1,5))/3.	


	! Creating the real Pseudo_values from the X_j 
	
	do iy=1,blocksize						! pseudo-values to estimates (eqn 5)
	  logpavg_m(iy) =log(avg_m(1))
     &	    +(blocksize-1)*(log(avg_m(1))-log(pavg_m(iy)))
	  logpstd_m(iy) =log(std_m(1))
     &	    +(blocksize-1)*(log(std_m(1))-log(pstd_m(iy)))
	  logpstd_d(iy) =log(std_d(1))
     &	    +(blocksize-1)*(log(std_d(1))-log(pstd_d(iy)))
	enddo 

	! jackknife estimates of statistics 
	do iy=1,blocksize						! bias correction 
	   ptheta_m(iy)=log(std_m(1))					! pseudo values of theta 
     &	        +(blocksize-1)*(log(std_m(1))-log(pstd_m(iy)))		! ...on monthly basis
	   ptheta_d(iy)=log(std_d(1))					! or 
     &	        +(blocksize-1)*(log(std_d(1))-log(pstd_d(iy)))		! ...on daily basis
	   pavg_m(iy)=avg_m(1)+(blocksize-1)*(avg_m(1)-pavg_m(iy))	! pseudo values of monthly mean 	
	   pstd_m(iy)=std_m(1)+(blocksize-1)*(std_m(1)-pstd_m(iy))	! pseudo values of monthly stddev 
	   pstd_d(iy)=std_d(1)+(blocksize-1)*(std_d(1)-pstd_d(iy))	! pseudo values of daily stddev
	   do ilag=1,nlag
	      pacf((ilag-1)*blocksize+iy) = acf(1,ilag)+(blocksize-1)
     &	           *(acf(1,ilag)-pacf((ilag-1)*blocksize+iy)) 	 	
	   enddo  	
	enddo 
	

	call pseu2est(pavg_m,blocksize,avg_m)
	call pseu2est(pstd_m,blocksize,std_m)
	call pseu2est(pstd_d,blocksize,std_d)
	do ilag=1,nlag
	   call pseu2est(pacf((ilag-1)*blocksize+1),
     &	        blocksize,acf(1,ilag))
	enddo ! ilag 

	write(*,'(a)') '# ' 
!	   write(*,'(a,3f8.2,3f8.3,a)')'# sim'
!     &	       , avg_m(1), std_m(1), std_d(1) 	! show results for current var 
!     &	       ,acf(1,1),acf(1,2),(acf(1,3)+acf(1,4)+acf(1,5))/3.	
!     &	       ,' jackknife' 


C======================================================================================================
C Store simulation stats 

	do i=1,2 
	   simavg_m(i)=avg_m(i)			
	   simstd_m(i)=std_m(i)
	   simstd_d(i)=std_d(i)
	   do ilag=1,nlag
	      simacf(i,ilag)=acf(i,ilag)
	   enddo
	enddo ! i 
	
	do iy=1,blocksize
	   simpavg_m(iy)=pavg_m(iy)
	   simpstd_m(iy)=pstd_m(iy)
	   simpstd_d(iy)=pstd_m(iy)	
	enddo 
	do ilag=1,nlag*blocksize
	   simpacf(ilag)=pacf(ilag)
	enddo


C======================================================================================================
C Calculate biases simulation with respect to the historical stats 
	biasavg_m=simavg_m(1)-hisavg_m(1)			! month average bias 
	biasstd_m=simstd_m(1)-hisstd_m(1)			! month stddev bias
	biasstd_d=simstd_d(1)-hisstd_d(1)			! day stddev bias 
	biasstd_m_rel=(simstd_m(1)/hisstd_m(1)-1.0)*100.0	! month stddev bias, percentage 
	biasstd_d_rel=(simstd_d(1)/hisstd_d(1)-1.0)*100.0	! day stddev bias, percentage 				
	do ilag=1,nlag
	   biasacf(ilag)=simacf(1,ilag)-hisacf(1,ilag)	! acf bias 
	enddo 

C Add biases to the TG en RR averages 

								! weight is the weight assigned to this variable in the average 
	if(TYPE_OF_VAR.eq.0) then 				! precipitation file 
	   RRbiasavg_m=RRbiasavg_m+biasavg_m
	   RRbiasstd_m=RRbiasstd_m+biasstd_m
	   RRbiasstd_d=RRbiasstd_d+biasstd_d
	   RRbiasstd_m_rel=RRbiasstd_m_rel+biasstd_m_rel
	   RRbiasstd_d_rel=RRbiasstd_d_rel+biasstd_d_rel
	   do ilag=1,nlag
	      RRbiasacf(ilag)=RRbiasacf(ilag)+biasacf(ilag)
	   enddo 
	   numRR=numRR+weight
	endif 

	if(TYPE_OF_VAR.eq.1) then 				! temperature file 
	   TGbiasavg_m=TGbiasavg_m+biasavg_m
	   TGbiasstd_m=TGbiasstd_m+biasstd_m
	   TGbiasstd_d=TGbiasstd_d+biasstd_d
	   TGbiasstd_m_rel=TGbiasstd_m_rel+biasstd_m_rel
	   TGbiasstd_d_rel=TGbiasstd_d_rel+biasstd_d_rel
	   do ilag=1,nlag
	      TGbiasacf(ilag)=TGbiasacf(ilag)+biasacf(ilag)
	   enddo 
	   numTG=numTG+weight
	endif 
	
	ENDDO		! READING FILES 
 188	continue 	! reached end of STDIN with filenames 



C======================================================================================================
C Process averages

	if(numRR.gt.0) then 				! precipitation file 
	   RRbiasavg_m=RRbiasavg_m/numRR
	   RRbiasstd_m=RRbiasstd_m/numRR
	   RRbiasstd_d=RRbiasstd_d/numRR
	   RRbiasstd_m_rel=RRbiasstd_m_rel/numRR
	   RRbiasstd_d_rel=RRbiasstd_d_rel/numRR
	   do ilag=1,nlag
	      RRbiasacf(ilag)=RRbiasacf(ilag)/numRR
	   enddo 
	   do iy=1,blocksize
	      RRhispavg_m(iy)=RRhispavg_m(iy)/numRR
	      RRhispstd_m(iy)=RRhispstd_m(iy)/numRR
	      RRhispstd_d(iy)=RRhispstd_d(iy)/numRR
	      RRhisptheta_m(iy)=RRhisptheta_m(iy)/numRR
	      RRhisptheta_d(iy)=RRhisptheta_d(iy)/numRR
	   enddo 
	   do ilag=1,nlag*blocksize
	      RRhispacf(ilag)=RRhispacf(ilag)/numRR
	   enddo
	endif 


	if(numTG.gt.0) then 				! temperature file 
	   TGbiasavg_m=TGbiasavg_m/numTG
	   TGbiasstd_m=TGbiasstd_m/numTG
	   TGbiasstd_d=TGbiasstd_d/numTG
	   TGbiasstd_m_rel=TGbiasstd_m_rel/numTG
	   TGbiasstd_d_rel=TGbiasstd_d_rel/numTG
	   do ilag=1,nlag
	      TGbiasacf(ilag)=TGbiasacf(ilag)/numTG
	   enddo 
	   do iy=1,blocksize
	      TGhispavg_m(iy)=TGhispavg_m(iy)/numTG
	      TGhispstd_m(iy)=TGhispstd_m(iy)/numTG
	      TGhispstd_d(iy)=TGhispstd_d(iy)/numTG
	      TGhisptheta_m(iy)=TGhisptheta_m(iy)/numTG
	      TGhisptheta_d(iy)=TGhisptheta_d(iy)/numTG
	   enddo 
	   do ilag=1,nlag*blocksize
	      TGhispacf(ilag)=TGhispacf(ilag)/numTG
	   enddo
	endif 
	
	! RRbias?????	now contains biases of several stats, averaged over all precipitation variables 
	! TGbias?????	now contains biases of several stats, averaged over all temperature variables 

	! RRhisp?????	pseudo-values of several stats, averaged over all precipitation variables  
	! TGhisp?????	pseudo-values of several stats, averaged over all temperature variables  
	
C======================================================================================================
C Pseudo vars to estimates 
		
	if(numRR.gt.0) then 				! precipitation file 
	   call pseu2est(RRhispavg_m,blocksize,RRavg_m) 
	   call pseu2est(RRhispstd_m,blocksize,RRstd_m)
	   call pseu2est(RRhispstd_d,blocksize,RRstd_d)
	   call pseu2est(RRhisptheta_m,blocksize,RRtheta_m)
	   call pseu2est(RRhisptheta_d,blocksize,RRtheta_d)
	   do ilag=1,nlag
	      call pseu2est(RRhispacf(ilag-1*blocksize+1),blocksize
     &	           ,RRacf(1,ilag))
	   enddo
	endif 

	if(numTG.gt.0) then 				! precipitation file 
	   call pseu2est(TGhispavg_m,blocksize,TGavg_m) 
	   call pseu2est(TGhispstd_m,blocksize,TGstd_m)
	   call pseu2est(TGhispstd_d,blocksize,TGstd_d)
	   call pseu2est(TGhisptheta_m,blocksize,TGtheta_m)
	   call pseu2est(TGhisptheta_d,blocksize,TGtheta_d)
	   do ilag=1,nlag
	      call pseu2est(TGhispacf(ilag-1*blocksize+1),blocksize
     &	           ,TGacf(1,ilag))
	   enddo
	endif 


C======================================================================================================
C Dump averaged results to STDOUT 


	   write(*,'(a)') '# RR'
	   write(*,'(a,3f8.3,4x,a)') 				! precipitation variables 
     &	     '# ',RRbiasavg_m,RRavg_m(1),RRavg_m(2)*2,'avg_m'
	   write(*,'(a,3f8.3,4x,a,4x,f8.3,a)') '# ',RRbiasstd_m_rel
     &	     ,RRstd_m(1),RRtheta_m(2)*2*100.0,'% = ',RRstd_m(2)
     &        ,' std_m'
	   write(*,'(a,3f8.3,4x,a,4x,f8.3,a)') '# ',RRbiasstd_d_rel
     &	     ,RRstd_d(1),RRtheta_d(2)*2*100.0,'% = ',RRstd_d(2)
     &        ,' std_d'
	   write(*,'(a)') '# '

	   write(*,'(a)') '# TG'
	   write(*,'(a,3f8.3,4x,a)') 				! precipitation variables 
     &	     '# ',TGbiasavg_m,TGavg_m(1),TGavg_m(2)*2,'avg_m'
	   write(*,'(a,3f8.3,4x,a,4x,f8.3,a)') '# ',TGbiasstd_m_rel
     &	     ,TGstd_m(1),TGtheta_m(2)*2*100.0,'% = ',TGstd_m(2)
     &        ,' std_m'
	   write(*,'(a,3f8.3,4x,a,4x,f8.3,a)') '# ',TGbiasstd_d_rel
     &	     ,TGstd_d(1),TGtheta_d(2)*2*100.0,'% = ',TGstd_d(2)
     &        ,' std_d'
	   write(*,'(a)') '# '

	   write(*,'(a6,a4,10a8)') '# ','lag','sim','his','2xse'
     &	      ,'sim-his' 

	   do ilag=1,nlag 
	      cvn=1.0
	      do jlag=1,ilag-1
	         cvn=cvn+2*(ilag-jlag)*hisacf(1,jlag)/ilag
	      enddo ! jlag 
	      write(*,'(a6,i4,10f8.3)') '  ',ilag 
     &	        ,simacf(1,ilag),hisacf(1,ilag)
     &	        ,hisacf(2,ilag)*2,biasacf(ilag),cvn/ilag
	   enddo 





C======================================================================================================
	END	! PROGRAM 





! --------------------- SUBROUTINES ---------------------------------
	subroutine pseu2est(arr_in,ny,res)
	! pseudo array of ny years to estimate 
	! estimate = res(1) +/- res(2),    that is : res(2) is the standard error of res(1)
	integer ny,iy 
	real sum,sumsqr,arr_in(ny),res(2)
	sum=0.0					
	sumsqr=0.0				
	do iy=1,ny
	   sumsqr=sumsqr+arr_in(iy)**2.
	   sum=sum+arr_in(iy)
	enddo 
	res(1)=sum/ny
	res(2)=((sumsqr-sum**2./ny)/(ny-1)/(ny))**0.5 		! (single) standard error in res(1)
								! result=res(1)+/-res(2)
	return 
	end 

      subroutine acjack(x, nj, nd, nl, r, r_j, se_j, pseu, var)
      implicit none
      integer nj, nd, nl
      real x(31*nj), r(nl), r_j(nl), se_j(nl)
      real pseu(nj,nl), var

!--------------------------------------------------------------------
! x(i)      :data array, length = number of days x number of years
! nj        :number of years	nj independend sets of nd samples  
! nd        :number of days (max 31)
! nl        :number of acc's to be calculated (lag 1..nl)
! r(i)      :sample estimate of lag i acc
! r_j(i)    :jackknife estimate of lag i acc
! se_j(i)   :jackknife standard error of lag i acc estimate
! pseu(j,i) :pseudovalues theta_j^* (see Equation (2) in Buishand
!            and Beersma (1996)) with theta = rho(i)
!            index j denotes the year that is omitted
!            index i denotes the lag of the acc
!            the pseudovalues are useful in the multivariate
!            extension (averaging over months and/or grid points)
! -------------------------------------------------------------------

      integer ii, i, j, k, n
      real c0, c(nl)
      real P0, P(nl)
      real S, Sl(nl), Sr(nl)
      real ac(nl), ac_d(nl)
! -------------------------------------------------------------------
!     [BB93 = Buishand and Beersma, 1993]
!     c0, c(i)          :c_k                 BB93 Eq (6)
!     P0, P(i)          :P_k                 BB93 Eq (6)
!     S, Sl(i), Sr(i)   :S, S_{L,k}, S_{R,k} BB93 Eq (6)
!     ac_d(i)           :r_{k(.)}            BB93 Eq (5)
! -------------------------------------------------------------------

      real cS, cP0, cP(nl), cSl(nl), cSr(nl)
      real Spseu

! -------------------------------------------------------------------
! -------------------------------------------------------------------


      S=0.0
      P0=0.0
      do i=1, nl
	 P(i)=0.0
	 Sl(i)=0.0
	 Sr(i)=0.0
	 ac_d(i)=0.0
      enddo
      n=nj*nd
      do i=1, n
	 P0=P0+x(i)**2
	 S=S+x(i)
      enddo
      c0=(P0 - S**2/real(n))/real(n)
      var=n*c0/real(n-1)
!	write(0,*) c0

      do k=1, nl
         do i=1, nj
            do j=(i-1)*nd+1, i*nd-k
               P(k)=P(k) + x(j)*x(j+k)
               Sl(k)=Sl(k) + x(j)
               Sr(k)=Sr(k) + x(j+k)
            enddo
         enddo
         c(k)=(P(k) - S*(Sl(k)+Sr(k)-real(nd-k)*S/real(nd))
     &         /real(n))/ real(n-nj*k)
!	 if (c0 .ne. 0.0) then 
	    r(k)=c(k)/c0
!	 else
!	    r(k)=1.0
!	 endif
      enddo

      n=(nj-1)*nd
      do ii=1, nj
	 cS=0.0
	 cP0=0.0
	 do i=(ii-1)*nd+1, ii*nd
	    cP0=cP0+x(i)**2
	    cS=cS+x(i)
	 enddo
	 cP0=P0 - cP0
	 cS=S - cS
	 c0=(cP0 - cS**2/real(n))/real(n)

	 do k=1, nl
	    cP(k)=0.0
	    cSl(k)=0.0
	    cSr(k)=0.0
	    do j=(ii-1)*nd+1, ii*nd-k
	       cP(k)=cP(k) + x(j)*x(j+k)
	       cSl(k)=cSl(k) + x(j)
	       cSr(k)=cSr(k) + x(j+k)
	    enddo
	    cP(k)=P(k) - cP(k)
	    cSl(k)=Sl(k) - cSl(k)
	    cSr(k)=Sr(k) - cSr(k)
	    c(k)=(cP(k) - cS*(cSl(k)+cSr(k)-real(nd-k)*cS/real(nd))
     &          /real(n))/ real(n-(nj-1)*k)
!	    if (c0 .ne. 0.0) then 
	       ac(k)=c(k)/c0
!	    else
!	       ac(k)=1.0
!	    endif
	    ac_d(k)=ac_d(k) + ac(k)/real(nj)
	    pseu(ii,k) = nj*r(k) - (nj-1)*ac(k)
	 enddo
      enddo

      do k=1, nl
	 r_j(k) = nj*r(k) - (nj-1)*ac_d(k)
!	 r_j(k) = 0.0
	 Spseu = 0.0
	 do ii=1, nj
!	    r_j(k) = r_j(k) + pseu(ii,k)/real(nj)
	    Spseu = Spseu + pseu(ii,k)**2
	 enddo
      	 se_j(k) = sqrt((Spseu - nj*r_j(k)**2) / real(nj*(nj-1)))
      enddo
	
      return
      end



      subroutine sdIjack(x, nj, nd, theta, svar)
      real x(nj*nd), theta(nj), svar
      integer nj, nd
	! theta(j) is the jth pseudo value for stdev of the month averages if 
	! j-th month left out 

      integer i, j
      real somy, somyy, y(100)
      real aj, aj1, aj2

      aj  = real(nj)
      aj1 = real(nj-1)
      aj2 = real(nj-2)

      somy=0.0
      somyy=0.0

      do j=1,nj
         y(j)=0.0
         do i=1, nd
            y(j)=y(j)+x((j-1)*nd+i)/real(nd)
         enddo
         somy=somy + y(j)
         somyy=somyy + y(j)**2
      enddo
      svar=(somyy - somy**2/aj)/aj1
      svar=sqrt(svar)	

      do j=1, nj
         theta(j) = (somyy-y(j)**2 - (somy-y(j))**2/aj1)/aj2
	 theta(j) = sqrt(theta(j))
      enddo

      return
      end

      subroutine sdPjack(x, nj, nd, theta, svar)
      real x(nj*31), theta(nj), svar
      integer nj, nd

      real somxx, somx, cxx, cx
      real an, an1
      integer i, j, l, n

      somxx=0.0
      somx=0.0

      n=nj*nd
      an = real(n)
      an1 = real(n-1)

      do i=1, n
         somxx=somxx+x(i)**2
         somx=somx+x(i)
      enddo
      svar=(somxx - somx**2/an)/an1
      svar=sqrt(svar)

      n=(nj-1)*nd
      an = real(n)
      an1 = real(n-1)

      do j=1, nj
         cxx=0.0
         cx=0.0
         do l=1, nd
            cxx=cxx + x((j-1)*nd+l)**2
            cx=cx + x((j-1)*nd+l)
         enddo
         theta(j) = (somxx-cxx - (somx-cx)**2/an)/an1
         theta(j) = sqrt(theta(j))
      enddo

      return
      end


	real function stdevblk(seq,nseq,nb)
	implicit none 
	! splits sequence seq(1:nseq) in nb blocks,
	! calculates the stdev of each and averages the results 
	integer nseq,nb,nbs,ib,i
	real seq(nseq)
	real avg,std,sum,sumsqr
	
	avg=0.0
	nbs=nseq/nb	! floor(nseq/nb)
	do ib=0,nb-1 
	   sum=0.0
	   sumsqr=0.0
	   do i=1,nbs
	      sum=sum+seq(i+ib*nbs)
	      sumsqr=sumsqr+seq(i+ib*nbs)**2.
	   enddo ! i 
	   std=((sumsqr-sum**2./nbs)/(nbs-1))**0.5	! stdev of current block 
	   avg=avg+std
	enddo ! ib 
	stdevblk=avg/nb
	return 
	end 


!======================================================================OLD
!	USED TO CALCULATE AUTOCORRELATION OVER THE ENTIRE SEQUENCE AT ONCE 
	subroutine autoc(x,rho,nx,nr)
	implicit none 
	! primitive autocorrelation subroutine 
	integer nx,nr,ix,ir
	real x(1:nx),rho(0:nr),sx
	do ir=0,nr
	   rho(ir)=0.0
	enddo
	sx=0
	do ix=1,nx-nr 
	   do ir=0,nr 
	      rho(ir)=rho(ir)+x(ix)*x(ix+ir)
	   enddo
	   sx=sx+x(ix)
	enddo
	do ir=0,nr 
	   rho(ir)=(rho(ir)-sx**2./(nx-nr))/(nx-nr-1)
	enddo
	do ir=nr,0,-1 
	   rho(ir)=rho(ir)/rho(0)
	enddo
	return 
	end 

	integer function instr(s,ssub)
	! same as 'index' 
	! it returns the position of substring ssub in string s 
	! 0 if not found 
	implicit none 
	integer j,ls,lsub
	character*(*) s 
	character*(*) ssub 
	ls=len(s)
	lsub=len(ssub)
	instr=0
	if(lsub.le.ls) then 
	   instr=0
	   do j=1,ls-lsub+1
	      if(s(j:j+lsub-1).eq.ssub) then 
	         instr=j
	         exit
	      endif
	   enddo
	endif 
	return 
	end 

