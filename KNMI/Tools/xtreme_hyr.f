	program xtreme_hyr
C	Calculate the n-day extremes of hydrological years 
	implicit none 
	integer maxyr,maxday 
	parameter(maxyr=60000)		! maximum number of years 
	parameter(maxday=40)		! maximum number of days 
	character*(10)	strarg 
	character*(40) regel 
	integer nday,cnr 		! number of days in total 
	integer instr 
	integer datum 
	integer nmod 
	integer iy,im,id,jy,oldyr,im0,ic
	real value,ndaysum,maxsum

	real sum30,maxsum30		! 30-day sum 
	real sum30d(maxyr)		! 30-day sum 
	real buf30(0:maxday-1)
	integer bufptr30

	integer maxdate 
	logical init 			! required to skip first halfyear 
	real buffer(0:maxday-1)		! buffer storing old values 
	integer bufptr			! pointer into buffer
	integer cnt
	integer date2cal,cal,cal0,cal2

	! extreme value arrays 
	integer keys(0:maxyr)
	real ymax(0:maxyr)		! year maxima 
	real ymax2(0:maxyr)		! year maxima 
	integer yrs(0:maxyr)
	integer datemax(maxyr)
	integer year(maxyr)
	integer yesterday
	
	cal0=date2cal(19010930,1)
! *****PROSSING ARGUMENTS ******************************************************* 
	cnr=2
	if(iargc().ge.1) then 
	  call getarg(1,strarg)		! arg[1] = column number 
	  read(strarg,*) cnr
	endif 

	nday=1
	if(iargc().ge.2) then 
	  call getarg(2,strarg)		! arg[2] = nday 
	  read(strarg,*) nday
	endif 
     	
!	produces an ordered list of annual n-day totals of a series 

	do bufptr=0,maxday-1
	   buffer(bufptr)=0.0 
	   buf30(bufptr)=0.0 
	enddo 
	bufptr=0
	bufptr30=0
	oldyr=-1
	ndaysum=0.0
	sum30=0.0


	ndaysum=0.0
	maxsum=0.0 
	cnt=1
	do 
 413	   read(*,'(a40)',end=313) regel 
	   if(instr(regel,'#').gt.0) goto 413
	   read(regel,*) datum, (value,ic=1,cnr-1)
	   iy=datum/10000
	   im0=im 
	   im=mod(datum,10000)/100		! month number 
	   id=mod(datum,100)			! day number 
	   ndaysum = ndaysum
     &	           + value
     &	           - buffer(bufptr) 
	   buffer(bufptr)=value 
	   bufptr=mod(mod(bufptr+1,nday)+nday,nday)

	   yesterday=cal2
	   cal=date2cal(datum,1)
	   cal2=mod(cal-cal0+365-1,365)+1
	   if(cal2.ge.nday) then
	      if(ndaysum.gt.maxsum) then 
                 maxsum=ndaysum 
		 maxdate=datum
	      endif 
	   endif 

	   if(cal2.lt.yesterday) then 			! add hydro-annual max to the ranks 
	         ymax(cnt)=maxsum 
	         datemax(cnt)=maxdate
	         sum30d(cnt)=maxsum30
	         year(cnt)=oldyr 
	         cnt=cnt+1 				! count years 
	      maxsum=0.0	
	   endif 
	   oldyr=iy 
	enddo 	
 313	continue 
	close(13) 

	cnt=cnt-1		! cnt pointing to the last added seasonmax 


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
     &	             ,sum30d(keys(iy)),keys(iy)
	enddo 

	
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

	integer function date2cal(datum,shift)
	implicit none 
	! convert date to day-of-year 
	! shift=1 : leave out 0731 in leapyears, iow shift between 0228 and 0801
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

	integer function nmod(ia,ib)
	implicit none 
	integer ia,ib
!	nmod=modulo(ia-1,ib)+1
	nmod=mod(mod(ia-1,ib)+ib,ib)+1
	return 
	end 


	
