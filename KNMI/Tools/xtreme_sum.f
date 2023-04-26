	program xtreme01
	implicit none 
	integer maxyr,maxday 
	parameter(maxyr=60000)		! maximum number of years 
	parameter(maxday=40)		! maximum number of days 
	character*(10)	strarg 
	character*(20) fnin 
	character*(40) regel 
	integer nday 			! number of days in total 
	integer instr 
	integer datum 
	integer nmod
	integer iy,im,id,jy,oldyr,im0,ic,cnr
	integer lastdays(0:12)
	data lastdays 
     &	   /0,31,59,90,120,151,181,212,243,273,304,334,365/
	integer nod(12)		! number-of-days in this month 
	data nod 
     &     /31,28,31,30,31,30,31,31,30,31,30,31/    
	real value,ndaysum,maxsum
	real maxavg,maxstd		! average and standard error of the maxima 

	real sum30,maxsum30		! 30-day sum 
	real sum30d(maxyr)		! 30-day sum 
	real buf30(0:maxday-1)
	integer bufptr30

	integer maxdate 
	logical init 			! required to skip first halfyear 
	real buffer(0:maxday-1)		! buffer storing old values 
	integer bufptr			! pointer into buffer
	integer cnt
	integer date2cal
	integer months(12)
	data months /0,0,0,1,1,1,1,1,1,0,0,0/ 	! select months with a '1' 
						! unselect months with a '0'

	integer firstmnd, lastmnd 
	integer season			! 0 is off-season 1 is on-season  	

	character*3 smonths(12) 
	data smonths /'Jan','Feb','Mar','Apr','May','Jun'
     &	             ,'Jul','Aug','Sep','Oct','Nov','Dec'/ 

	! extreme value arrays 
	integer keys(0:maxyr)
	real ymax(0:maxyr)		! year maxima 
	real ymax2(0:maxyr)		! year maxima 
	integer yrs(0:maxyr)
	integer datemax(maxyr)
	integer year(maxyr)
	integer iarg 
        logical arglogical
	character*(15) my_name		! name of the program 
	

! *****PROSSING ARGUMENTS ******************************************************* 
*
        cnr=2
        if(iarg().ge.1) then 
          call getarg(1,strarg)         ! arg[1] = column number 
          read(strarg,*) cnr
        endif 

        nday=1
        if(iarg().ge.2) then 
          call getarg(2,strarg)         ! arg[2] = nday 
          read(strarg,*) nday
        endif 
        
!       produces an ordered list of annual n-day totals of a series 

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
	sum30=0.0
	maxsum=0.0 
	maxsum30=0.0
	cnt=1
	season=0
	im=-1
	init=.False. 
	do while (.TRUE.)
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
!	   bufptr=modulo(bufptr+1,nday)
	   bufptr=mod(mod(bufptr+1,nday)+nday,nday)

	   sum30 = sum30			! keep thrack of 30-day sum 
     &	           + value
     &	           - buf30(bufptr30) 
	   buf30(bufptr30)=value 
!	   bufptr30=modulo(bufptr30+1,30)
	   bufptr30=mod(mod(bufptr30+1,30)+30,30)

	      if(ndaysum.gt.maxsum .and. months(im).gt.0) then 
	         if(months(im).gt.months(nmod(im-1,12)) 	! if it is the first month of the season, check 
     &	             .and. id.lt.nday) then 			! whether ALL days in the n-day sum are WITHIN the season 
	         else 
	            maxsum=ndaysum 
		    maxdate=datum
	            maxsum30=sum30
	   	 endif 
	      endif 

!	   if(iy.ne.oldyr .and. oldyr.gt.0) then	! data from new year 
	   if(months(im).lt.months(im0)) then 		! marking end of season 
	      if(im0.gt.0 .and. init) then  		           ! end of season 
	         ymax(cnt)=maxsum 
	         datemax(cnt)=maxdate
	         sum30d(cnt)=maxsum30
	         year(cnt)=oldyr 
	         cnt=cnt+1 				! count years 
	      endif 
	      maxsum=0.0	
	   endif 
	   if(months(im).gt.months(im0)) then 		! beginning of season 
	      init=.True. 
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
	         goto 178 
	      endif  
	   enddo 
 178	   continue 
	enddo 		! ...resulting in a list of sorted indices 
	
	maxavg=0.0 
	maxstd=0.0 

	do iy=cnt,1,-1		! dump loop : writing ordered list to stdout 
           if(arglogical('-floats')) then 
	      write(*,'(2f12.4,1i6,i10,i5,f8.1,i10)') 
     &	        -log(-log(((cnt-iy+1)-0.3)/(cnt+0.4)))	! reduced gumbel variate 
     &	             ,ymax(keys(iy))
     &	             ,mod(year(keys(iy)),10000)
     &	             ,datemax(keys(iy))
     &	             ,date2cal(datemax(keys(iy)),1)
     &	             ,sum30d(keys(iy)),keys(iy)
           else 
	      write(*,'(2E12.4,1i6,i10,i5,f8.1,i10)') 
     &	        -log(-log(((cnt-iy+1)-0.3)/(cnt+0.4)))	! reduced gumbel variate 
     &	             ,ymax(keys(iy))
     &	             ,mod(year(keys(iy)),10000)
     &	             ,datemax(keys(iy))
     &	             ,date2cal(datemax(keys(iy)),1)
     &	             ,sum30d(keys(iy)),keys(iy)
           endif 
	   maxavg=maxavg+ymax(keys(iy))
	   maxstd=maxavg+ymax(keys(iy))**2.
	enddo 

	maxstd=(maxstd-(maxavg**2.)/cnt)/(cnt-1)
	maxstd=maxstd**0.5				
	maxavg=maxavg/cnt
	end 

        logical function arglogical(prefix)
C       returns .True. if the prefix is found in ARGV[]
        implicit none
        integer jarg
        integer instr, iarg
        logical result
        character*(*) prefix
        character*50 sarg

        result=.False.
        jarg=1
        do jarg=1,iarg()
           call getarg(jarg,sarg)
           if(sarg(1:instr(sarg,' ')-1).eq.prefix) result=.True.
        enddo
        arglogical=result
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
	         goto 231
	      endif
	   enddo
 231	   continue 
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


	
	
	integer function iarg()	
	! returns the number of non-empty cmdl args 
	integer count 
	integer ifun 
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

