	program scangap
	! opens an hourly data file and reports gaps 
	use msflib
	implicit none 
	character*(30) fnin
	character*(30) readstr 
	integer instr,datum
	real values(10)
	integer date2cal
	integer ndx,ndxold,gap,maxgap
	integer datumold,first,last,missing,total

	call getarg(1,fnin)
	open(15
     &	   ,file=fnin(1:instr(fnin,' ')-1)
     &	   ,status='old')
	write(*,*) fnin(1:instr(fnin,' ')-1)//' opened'

	missing=0
	first=0 
	maxgap=0

	do 
 215	   read(15,'(a30)',end=315,err=215) readstr
	   if(instr(readstr,'#').le.0) then 
	      read(readstr,*) datum,values(1)
	      if(mod(datum,10000).eq.0731 
     &	         .and. mod(datum/10000,4).eq.0) goto 215
     	         ndx=date2cal(datum,1)+(datum/10000-1900)*365	! Juli 31st is ignored in leapyears  
	         if(first.eq.0) then  
	            first=ndx
	         else 
	            if(ndx.ne.ndxold+1) then 
	               gap=ndx-ndxold-1
	               missing=missing+gap
	               if(gap.gt.maxgap) maxgap=gap
	               write(*,*) datumold,datum,gap,missing
	            endif 
	         endif 
	   ndxold=ndx
	   datumold=datum
	   endif 
	enddo 
 315	continue 
	close(15)
	last=ndx
	total=last-first+1
	write(*,*) 
	write(*,'(a,i8,a,i8,a,f8.3,a,i5,a)') 
     &     ' Missing ',missing,' of ',total,' = '
     &	     ,missing*100.0/total,' %   maxgap = '
     &       ,maxgap,' days'  
	write(*,*) 
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

	
