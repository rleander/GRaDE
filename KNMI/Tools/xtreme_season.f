	program xtreme_season
C	Similar to xtreme_win, returns the ordered n-day max of 
C	the defined season 

	implicit none 
	integer maxyr,maxday 
	parameter(maxyr=60000)		! maximum number of years 
	parameter(maxday=40)		! maximum number of days 
	character*(10)	sarg 
	character*(40) regel 
	integer nday,cnr 		! number of days in total 
	integer argint, iarg
	integer instr 
	integer datum 
	integer iy,jy,ic
	real value,ndaysum,maxsum

	integer maxdate 
	real buffer(0:maxday-1)		! buffer storing old values 
	integer bufptr			! pointer into buffer
	integer cnt
	integer date2cal

	logical inseason, inseason_shift 
	integer date0, date1
	integer seasondays 

	! extreme value arrays 
	integer keys(0:maxyr)
	real ymax(0:maxyr)		! year maxima 
	real ymax2(0:maxyr)		! year maxima 
	integer yrs(0:maxyr)
	integer datemax(maxyr)
	integer year(maxyr)
	

! *****PROSSING ARGUMENTS ******************************************************* 
	cnr=argint('-nc',2)		! column number
	nday=argint('-nd',1)		! n-day 
	if(iarg().lt.2) then 
	   write(0,*) 'xtreme_season date0 date1 [options]'
	   write(0,*) 'Returns seasonal max of n-day totals'
	   write(0,*) 'for the season between date0 and date1'
	   write(0,*) 'Options : ' 
	   write(0,*) '     -nd days	'
	   write(0,*) '     -nc column (from which to read data)'
	   write(0,*)
	   write(0,*) 'Reads from STDIN, writes to STDOUT' 
	   write(0,*)
	   stop '####' 
	endif 

	call getarg(1,sarg)
	read(sarg,*) date0
	call getarg(2,sarg)
	read(sarg,*) date1
	
!	produces an ordered list of annual n-day totals of a series 

	do bufptr=0,maxday-1
	   buffer(bufptr)=0.0 
	enddo 
	bufptr=0
	ndaysum=0.0

	maxsum=0.0 
	cnt=0
	seasondays=0

	do while(.True.)
 413	   read(*,'(a40)',end=313) regel 
	   if(instr(regel,'#').gt.0) goto 413	! skip commented lines 
	   read(regel,*) datum, (value,ic=1,cnr-1)
	   ndaysum = ndaysum
     &	           + value
     &	           - buffer(bufptr) 
	   buffer(bufptr)=value 
!	   bufptr=modulo(bufptr+1,nday)
	   bufptr=mod(mod(bufptr+1,nday)+nday,nday)

C	   if(inseason(date0,date1,datum)) then 	
           if(inseason_shift(date0,date1,datum,nday)
     &        .and. inseason_shift(date0,date1,datum,1)) then 
	      seasondays=seasondays+1
              if(ndaysum.gt.maxsum) then 
                 maxsum=ndaysum
                 maxdate=datum
              endif 
	   else 
	      if(seasondays.gt.0) then 
	         ymax(cnt+1)=maxsum 
	         datemax(cnt+1)=maxdate
	         year(cnt+1)=datum/10000
	         cnt=cnt+1 				! count years 
	      endif 
	      seasondays=0
	      maxsum=0.0	! set minimumvalue to zero 
	   endif 
	enddo 	
 313	continue 
	close(13) 

	! copy unsorted yearmaxima
	do iy=1,cnt
	   ymax2(iy)=ymax(iy)
	   yrs(iy)=mod(datemax(iy)/10000,100)
	enddo 

	! sort yearmaxima 
	ymax(0)=1E+07-1
	keys(0)=0 
	do iy=1,cnt
	   do jy=iy,1,-1 
	      if(ymax(iy).gt.ymax(keys(jy-1))) then 
	         keys(jy)=keys(jy-1)
	      else 
	         keys(jy)=iy
	         exit 
	      endif  
	   enddo 
	enddo 		! ...resulting in a list of sorted indices 
	
	do iy=cnt,1,-1		! dump loop : writing ordered list to stdout 
	   write(*,'(2E12.4,1i6,i10,i5,f8.1,i10)') 
     &	-log(-log(((cnt-iy+1)-0.3)/(cnt+0.4)))	! reduced gumbel variate 
     &	             ,ymax(keys(iy))
     &	             ,mod(year(keys(iy)),10000)
     &	             ,datemax(keys(iy))
     &	             ,date2cal(datemax(keys(iy)),1)
C     &	             ,sum30d(keys(iy)),keys(iy)
	enddo 
	
	end 

	
	include 'routines.f' 	
