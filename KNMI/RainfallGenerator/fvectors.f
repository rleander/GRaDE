	Program FVECTORS
C 	Purpose : producing feature-vectors spanning the required historical period 
C 	Syntax  : FVECTORS  firstyr  lastyr 
C	          The arguments firstyr and lastyr defining the historical interval for 
C	          which to calculate the averages are compulsory. 
C	Input   : firstyr (format: YYYY) 
C	          lastyr (format: YYYY) 
C	          STDIN (list of filenames and weights) 
C	          Files containing historical rainfall and temperature data, as specified in 
C	          STDIN (see 'Requires'). 
C	Output  : STDOUT (fvector-file) 
C
C	build with : g77 fvectors.f -o fvectors -Wall
C
C
C	(c) KNMI, De Bilt, January 2006
C	Robert Leander
	
	implicit none 
	logical bnd 
	integer iarg, instr, date2cal
	integer MAXSERIES, MAXYR, YR0
	parameter (MAXSERIES=100)	! maximum number of historical series 
	parameter (MAXYR=200)		! maximum number of historical years  
	parameter (YR0=1900)		! year before the first year used in any simulation  

					! YR0 is assumed to precede any possible historical date,
					! i.e. in the indices for local storage in the program 
					! is based on the historical date from which the year YR0
					! is subtracted, i.e. array-index 1 corresponds to  
					! the 1st of Jan of YR0+1.

					! MAXYR determines the time-dimension of arrays for storage,
					! which is MAXYR*365 days, calculated from the 
					! abovementioned index 1.
	real thresh 
	parameter(thresh=0.3)		! wet-day threshold to determine mean-wet-day rainfall
	character*(15) fnin		! filename of input file 
	character*(15) fn_arr(MAXSERIES)! array of input filenames
	character*(70) regel 		
	character*(15) my_name		! program name 
	integer datum
	integer firstyr, lastyr, ny  	! first year, last year, number of years 
	real wt(MAXSERIES)		! weights array
	real wtsum_P,wtsum_T		! sum of weights for rainfall and temperature 	
	
	integer p_or_t(MAXSERIES)	! keeps track whether a series is rainfall or temperature
					! +1 : rainfall 
					! -1 : temperature 
	integer numP,numT		! number of rainfall stations, 
					! number of temperatures station 
	character*(6) paramstr1		! command line parameter 
	real value			! value read from input file 
	real dayseries(MAXYR*365,MAXSERIES)	! daily series loaded into memory 
	real smooths(365,2,MAXSERIES)	! smoothed calendar day characteristics 
					! 1: shift; mean temperature or zero for rainfall
					! 2: scale; stddev. of temperature or mean-wet-day rainfall
	real mean(365)			! temporary calendar day mean array of 365 days 
	real var(365)			! temporary calendar day variance array of 365 days 
	
					! 365 days a year, July 31 is excluded in leapyears 
	integer iy, id, calday  	! year, day, calendar day 
	integer ifile, nfile		! file-index, number of files 
	integer nwet(365)
	integer datums(MAXYR*365)		! datum array 
					! four feature vector elements 
	real Pavg(MAXYR*365)		! Average standardized daily rainfall 
	real Tavg(MAXYR*365)		! Average standardized temperature 

	real Fwet(MAXYR*365)		! fraction of wet stations  	
					! is not used in this implementation, but can be used 
					! as a 4th element in 4D resampling, compatibility with 
					! older versions. 

	integer spanP,spanT		
					! calendarday 
	parameter(spanP=45)
	parameter(spanT=30)
	real kernlP(-spanP:spanP)	! Epanechnikov kernel rainfall
	real kernlT(-spanT:spanT)	! Epanechnikov kernel temperature
	
! ******PROCESSING COMMAND-LINE ARGUMENTS ****************************
        call getarg(0,my_name)		! What's my name ? 
	if(iarg().lt.1) then 
	   write(0,*) 'syntax : '//my_name(1:instr(my_name,' ')-1)
     &	   // '  firstyear lastyear '
	   write(0,*) 
	   stop '...'
	endif 

	   call getarg(1,paramstr1)	
	   read(paramstr1,*) firstyr 		! arg1 = first historical year 
	   call getarg(2,paramstr1)
	   read(paramstr1,*) lastyr 		! arg2 = last historical year 
						! firstyr and lastyr delimit the historical 
						! period one wants to use for resampling and 
						! what is commonly available for all data

	
	ny=lastyr-firstyr+1
	write(0,*) firstyr,' to ',lastyr,' : '
     &	     ,ny,' years '



! ******PROCESSING INPUT FILES ***************************************
	ifile=0
	numP=0					
	numT=0				

	do while(.TRUE.)	! begin reading filenames 
	ifile=ifile+1
	do id=1,MAXYR*365
	   dayseries(id,ifile)=-999.99      
	enddo ! id  
	
 315	read(*,'(a70)',end=666) regel
	if(instr(regel,'#').gt.0) goto 315		! lines in the filelist containing # 
							! are ignored 
	read(regel,*,err=315,end=415) fnin
	fn_arr(ifile)=fnin
 415    open(15,file=fnin(1:instr(fnin,' ')-1),status='OLD',err=315)
	write(0,'(i5,a$)') ifile,' : '//fnin(1:instr(fnin,' ')-1)//' '
	wt(ifile)=1.0					! wt defaults to 1 if not found 
	read(regel(instr(regel,' ')+1:len(regel)),*,end=121) wt(ifile)
 121	continue 
	if(instr(fnin,'.rr').gt.0 .or. instr(fnin,'.pr').gt.0) then  ! rainfall file 
	   numP=numP+1
	   p_or_t(ifile)=1				! 1 for rainfall station 
	   write(0,*) ' rainfall station'
	endif 
	if(instr(fnin,'.t').gt.0) then 			! temperature file 
	   write(0,*) ' temperature station'
	   numT=numT+1
	   p_or_t(ifile)=-1				! -1 for temperature  
	endif 
	do while (.TRUE.)		! begin reading from the current file 
 215	   read(15,'(a30)',end=116,err=215) regel
	   if(instr(regel,'#').le.0) then 		! ignore '#'-lines 
	      read(regel,*,err=215,end=215) datum,value
	      iy=datum/10000-YR0			! what year ? 
	      if(iy.lt.1 .or. iy.gt.MAXYR) goto 215
	      if(mod(datum,10000).eq.0731 		! skip July 31 in leap years 
     &	         .and. mod(iy,4).eq.0) goto 215
	      calday=date2cal(datum,1)			! obtain calendar day 
	      if(calday.ge.1) then 
	         datums(calday+(iy-1)*365)=datum 	! store date for later use 
	         dayseries(calday+(iy-1)*365,ifile)=value
	      endif 
	   endif 
 	enddo 				! end reading from the current file 
 116	continue 
	close(15)
      enddo  			! end reading filenames
 666  continue 
	
	nfile=ifile-1 					! store number of files 
	write(0,*) nfile,' files found '  
	write(0,*) numP,' rainfall records' 
	write(0,*) numT,' temperature stations' 

! ***** CALCULATE CALENDAR DAY STATS *********************************

	do iy=firstyr-YR0,lastyr-YR0 
	   do calday=1,365 
	      fwet(calday+(iy-1)*365)=0.0	! zero the fraction of wet-stations 
						! before looping over files 			
	   enddo ! calday 
	enddo ! iy 
	call epaknl(kernlP,spanP)		! prepare smoothing kernel for rainfall 	
	call epaknl(kernlT,spanT)		! prepare smoothing kernel for temperature

	do ifile=1,nfile 
	   do calday=1,365 
	      smooths(calday,1,ifile)=0.0 
	      smooths(calday,2,ifile)=0.0 
	   enddo ! calday 


	   !-----------------HANDLE RAINFALL DATA 
	   if(p_or_t(ifile).gt.0) then 		
	      do calday=1,365
	         mean(calday)=0.0 		! mean-wet-day rainfall 
	         nwet(calday)=0
	      enddo ! calday 
	      do iy=firstyr-YR0,lastyr-YR0
	         do calday=1,365 
	            value = dayseries(calday+(iy-1)*365,ifile)
	            if(bnd(value,0.,900.)) then
	               if(value.gt.thresh) then 	
	                  mean(calday)=mean(calday)+value
	                  nwet(calday)=nwet(calday)+1
	                  fwet(calday+(iy-1)*365)= 		
     &	                       fwet(calday+(iy-1)*365)+1./numP
	               endif ! wet day 
	            endif ! valid data 
	         enddo ! caldays 
	      enddo ! iy 
	      do calday=1,365 
	         if (nwet(calday).gt.0) then 
	             mean(calday)=mean(calday)/nwet(calday)	
	         endif 
	      enddo 
	      call perconv(mean,kernlP,smooths(1,2,ifile),365,spanP) 	! smooth with kernel 
	      ! SCALE =mean-wet-day rainfall --> column 2 	
	      ! SHIFT = 0 --> column 1 

	      ! NB. smoothing of the calendar-day mean-wet-day amounts is performed
	      !     as a periodic convolution with the smoothing kernel 
	   endif ! rainfall variable 	


	   !-----------------HANDLE TEMPERATURE DATA 
	   if(p_or_t(ifile).lt.0) then 		
	      do calday=1,365
	         mean(calday)=0.0 			! mean Temp. 
	         var(calday)=0.0 			! variance of Temp. 
	      enddo 
	      do iy=firstyr-YR0,lastyr-YR0 
	         do calday=1,365 
	            value = dayseries(calday+(iy-1)*365,ifile) 
	            if(bnd(value,-100.0,900.0)) then 
	               mean(calday)=mean(calday)
     &	                 + dayseries(calday+(iy-1)*365,ifile)
	               var(calday)=var(calday)
     &	                 + dayseries(calday+(iy-1)*365,ifile)**2. 
	            endif 
	         enddo ! calday 
	      enddo ! iy 
	      do calday=1,365 
	          var(calday)=
     &	           (var(calday)-mean(calday)**2./ny)	! variance 
     &	               /(ny-1)
	          var(calday)=var(calday)**0.5 		! variance --> stddev 
	          mean(calday)=mean(calday)/ny		! mean 
	      enddo ! calday 
	      call perconv(mean,kernlT,smooths(1,1,ifile),365,spanT) ! smooth mean 
	      call perconv(var,kernlT,smooths(1,2,ifile),365,spanT)  ! smooth var
	      ! SCALE = stddev of temperature --> column 2 	
	      ! SHIFT = mean temperature --> column 1 
	   endif ! temperature variable 


	   ! write smoothed seasonal scale and translation to files 

	   fnin=fn_arr(ifile)
	   if(p_or_t(ifile).gt.0) then 
	   open(13,file=fnin(1:instr(fnin,'.')-1)//'.smr'		! *.smr = rainfall 
     &	        ,status='REPLACE')
	   else
	   open(13,file=fnin(1:instr(fnin,'.')-1)//'.smt'		! *.smt = temperature 
     &	        ,status='REPLACE')
	   endif 

	do calday=1,365 
	   write(13,'(i6,100f8.3)') calday		! NOTE : smooths(calday,1,ifile)=0 in
     &	          ,smooths(calday,1,ifile)		! .smr files 
     &	          ,smooths(calday,2,ifile)
	enddo ! calday 
	close(13)

	enddo ! ifile 

! ***** STANDARDIZATION **********************************************
	write(0,*) 'standardizing' 
	! standardize daily series 
	! all variables are standardized by subtracting their 
	! calendar day offset and dividing by their calendar day scale 
	! For a rainfall variable  offset= zero 
	!                          scale = mean-wet-day amount 
	! For a temperature variable offset = mean 
	!                            scale = standard deviation 

	do ifile=1,nfile 
	   do iy=firstyr-YR0,lastyr-YR0
	      do calday=1,365 
	         if(abs(dayseries(calday+(iy-1)*365,ifile)).lt.900) then 
	         dayseries(calday+(iy-1)*365,ifile) = 
     &	           (
     &	              dayseries(calday+(iy-1)*365,ifile)		
     &	                  -smooths(calday,1,ifile)
     &	           )
     &	             /smooths(calday,2,ifile)
	          else 
	          endif 
	      enddo ! calday 
	   enddo ! iy 
	enddo ! ifile 		

	wtsum_P=0.0
	wtsum_T=0.0
	do ifile=1,nfile
	   if(p_or_t(ifile).gt.0) wtsum_P=wtsum_P+wt(ifile)
	   if(p_or_t(ifile).lt.0) wtsum_T=wtsum_T+wt(ifile)	
	enddo 

! ***** AVERAGE STANDARDIZED RAINFALL AND TEMPERATURE ****************
	numP=max(numP,1)		! prevent division-by-zero 
	numT=max(numT,1)
	! calculate averages of standardized rainfall and temperture 
	do iy=firstyr-YR0,lastyr-YR0
	   do calday=1,365 
	      id=calday+(iy-1)*365
	      Pavg(id)=0.0 
	      Tavg(id)=0.0 
	      do ifile=1,nfile					! loop over files 
	         if(abs(dayseries(id,ifile)).gt.900) then 	! bogus value alert
	             Pavg(id)=-999.99
	             Tavg(id)=-999.99
	             goto 112
	         endif 
	         if(p_or_t(ifile).gt.0) then 
	             Pavg(id)=Pavg(id)+dayseries(id,ifile)*wt(ifile)	! weighted rainfall average
	         else 
	             Tavg(id)=Tavg(id)+dayseries(id,ifile)*wt(ifile)	! weighted temperature average
	         endif 
	      enddo ! ifile 
	      Pavg(id)=Pavg(id)/wtsum_P
	      Tavg(id)=Tavg(id)/wtsum_T
 112	      continue
	      write(*,'(i10,i5,4e12.4)') YR0+iy,calday		! year and calendar day 
     &	          ,Pavg(id)					! 1st feature-vector element
     &	          ,Tavg(id)					! 2nd        "
     &	          ,fwet(id)					! 3rd        "
	   enddo ! calday 
	enddo ! iy 


	end



!********************************************************************************

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
	         goto 380
	      endif
	   enddo
 380	   continue
	endif 
	return 
	end 

	subroutine epaknl(knl,span)
	implicit none 
	! fill knl(-span:span) with Epanechnikov smoothing kernel 
	integer span, i 
	real knl(-span:span), som, alpha
	
        som = 1.0
        knl(0) = 1.0 
        do i=1, span
	    alpha=(i*1.0/span)
            knl(i) = 1-alpha**2
            knl(-i) = knl(i)
            som = som + 2*knl(i)
        enddo
        do i=-span, span
            knl(i) = knl(i)/som
        enddo
	return 
	end 

	integer function modmod(a,b)
	implicit none 
	integer a,b 
	modmod=mod(mod(a,b)+b,b)
	return 
	end 

	integer function nmod(x,y)
	! returns x mod y as a number in the range 1..y
	implicit none 
	integer x,y,modmod
	nmod=modmod(x-1,y)+1
	end 

	subroutine perconv(x,y,res,nx,ny)
	implicit none 
	! performs a periodical convolution on x(1:nx) 
	! with array y(-ny:ny) 
	! intended for smoothing x 
	! output in res(1:nx) 
	integer nx, ny, l, i, k, nmod  
	real x(1:nx), y(-ny:ny), res(1:nx)

        do l=1, nx 
	   res(l)=0.0
           do i=-ny, ny
              k = nmod(l+i,nx) 			
	      res(l) = res(l)+x(k)*y(i)
           enddo ! i 
	enddo ! j 
	return
	end 

	integer function date2cal(datum,shift)
	implicit none 
	! convert date to day-of-year 
	! shift=1 : leave out 0731 in leapyears, iow shift 1 day between 0228 and 0801
	! shift=0 : leave out 0229 in leapyears 
	integer dayzero1(12)	! day zero in each month 
	integer yr,day,mnth,datum,shift 
	integer dayzero2(12)	! day zero in each month 
	data dayzero1 /0,31,59,90,120,151,181,212,243,273,304,334/ 	!no leap or feb 29 ignored 
	data dayzero2 /0,31,60,91,121,152,182,212,243,273,304,334/ 	!leap AND jul 31 ignored 
	yr=datum/10000
	mnth=mod(datum,10000)/100
	day=mod(datum,100)
	if(mod(yr,4).eq.0 .and. shift.eq.1) then 
	   date2cal=dayzero2(mnth)+day
	else 
	   date2cal=dayzero1(mnth)+day
	endif 
	return 
	end 

	logical function bnd(x,x0,x1)
	! returns true if x element of [x0,x1], false otherwise
	implicit none 
	real x,x0,x1
	if(x.gt.x1 .or. x.lt.x0) then 
	     bnd=.False.
	else 
	     bnd=.True.
	endif
	return 
	end 


	integer function iarg()	
	! returns the number of non-empty command line arguments 
	integer count 
	character*(10) sarg
	count=0
	call getarg(count+1,sarg) 
	do while (sarg(1:1).gt.' ')
	   count=count+1
	   call getarg(count+1,sarg) 
	enddo 
 	continue 
	iarg=count 
	return
	end
