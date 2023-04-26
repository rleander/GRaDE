	program scanlist
	! opens and scans the files in the list of names read from STDIN
	! The program evokes GNUPlot, which produces a Postscript image 
	! with bars representing uninterrupted segments of daily sequences.
	! Through this plot interruptions in the sequences are easily 
	! detectable. 
	! 
	! (c) KNMI, De Bilt January 2006
	! Robert Leander
	!
	! build with: g77 scanlist.f -o scanlist -Wall 	
 
	implicit none 
	character*(*) GNUPLOT
	parameter (GNUPLOT='gnuplot')  ! name of the GNUPlot program on this system 

	character*(30) fnin
	character*(60) readstr 
	integer instr,datum
	real values(10)
	integer date2ndx
	integer ndx0,ndxold0
	integer ndx,ndxold,gap,maxgap,maxgapndx
	integer datumold,first,missing
	integer leapyrs
	parameter (leapyrs=1)
	character*(5) args
	real x1,x0,y1,y0
	real dum 
	integer colnr 

	integer ifile,i
	real lat,lon
	integer iloc,oldloc             ! location number 
	character*(20) namestr 
	real absurd

	open(32,file='databar.gnu')
	write(32,*) 'pl ''scan.gnu'' u 1:2 w l'
	
	write(32,*) 'set grid x'
	write(32,*) 'set noborder'
	write(32,*) 'set nokey'
	write(32,'(a$)') 'set yt ('
	
	open(33,file='scan.gnu',status='REPLACE')! gnuplot output 

	colnr=1							! columnr AFTER the date and hour column
	call getarg(1,args)					! default = first 
	read(args,*,err=911,end=911) colnr
 911	continue 

	oldloc=0
	iloc=0 
	ifile=0 
	do 
 266	   read(*,*,end=666) fnin	! read filename from stdin 
	   if(instr(fnin,'#').gt.0) cycle
	   if(instr(fnin,'&').gt.0) then 
	      ifile=ifile+1
	      cycle
	   endif ! empty line 
	   write(0,*) fnin
!	   read(fnin(1:instr(fnin,'.')-1),*) locnr
!	   locnr=locnr/10		! designation of stations location 
	   open(15,file=fnin(1:instr(fnin,' ')-1),status='old'
     &	    ,err=266)
	   ifile=ifile+1
	   if(ifile.eq.1) then 
	      write(32,'(a,i3$)') 
     &	          ''''//fnin(1:instr(fnin,' ')-1)//'''',-ifile
	   else 
	      write(32,'(a,i3$)')
     &	          ','''//fnin(1:instr(fnin,' ')-1)//'''',-ifile
	   endif 

	   missing=0
	   first=0 
	   maxgap=0
	   maxgapndx=0
	   absurd=999.0		! totally insane values not to be taken seriously 


	   lat=1.00
	   lon=1.00

	   namestr=fnin(1:instr(fnin,'.')-1)
	   y0=-ifile-0.4		! bottom of bar 
	   y1=-ifile+0.4		! top of bar 
	   first=0
	     do 
 215	        read(15,'(a60)',end=315) readstr
 	        if(instr(readstr,'#').le.0) then 	! contains name or coordinates 
                         read(readstr,*,err=215,end=215) 
     &	                     datum,(dum,i=1,colnr-1),values(1)
	            if(values(1).le.-90.0) goto 215	! absurd value -90 : detecting missing value 
	            if(values(1).ge.absurd) goto 215
	            if(leapyrs.le.0) then 
	               if(mod(datum,10000).eq.0731
     &	                   .and. mod((datum/10000),4).eq.0) goto 215
	            endif 
	            ndx=date2ndx(datum,leapyrs)
	            ndx0=date2ndx(datum,0)
	           
	            if(first.eq.0) then  
	               first=ndx
	               x0=ndx0*1.0/365.
	            else 
	               if(ndx.lt.ndxold+1) write(*,*) ndx,datum,'!!!'
	               if(ndx.gt.ndxold+1) then 		! found a gap 
	                  gap=ndx-ndxold-1

	                  write(*,'(a,2i12)') ' ',datum,gap

	                  missing=missing+gap
	                  if(gap.gt.maxgap) then 
	                     maxgap=gap
	                     maxgapndx=ndxold
	                  endif 
	                  x1=ndxold0*1.0/365			! end of bar 
			  write(33,'(2f8.3)') x0,y0		! draw a rectangular bar in gnuplot 
			  write(33,'(2f8.3)') x0,y1
			  write(33,'(2f8.3)') x1,y1
			  write(33,'(2f8.3)') x1,y0
			  write(33,'(2f8.3)') x0,y0
	                  write(33,*)
	                  x0=ndx0*1.0/365			! beginning of new bar 
	               endif 
	            endif 
	        ndxold=ndx
	        ndxold0=ndx0
	        datumold=datum
	        endif 
	     enddo ! days 
 315	     continue 
	   close(15)

	     x1=ndxold0*1.0/365			! end of bar 
	     write(33,'(2f8.3)') x0,y0		! draw a rectangular bar in gnuplot 
	     write(33,'(2f8.3)') x0,y1
	     write(33,'(2f8.3)') x1,y1
	     write(33,'(2f8.3)') x1,y0
	     write(33,'(2f8.3)') x0,y0
	     write(33,*)
	enddo 	! done with this file, on to the next 
 666	continue 
	close(33)

	write(32,'(a1)')')'
	write(32,*) 'rep' 
	write(32,*) 'set term pos eps enhanced ''CourierBold'''
	write(32,*) 'set out ''scanlist.ps'''
	write(32,*) 'rep' 
	write(32,*) 'set out'
	write(32,*) 'set term win'	! Windows
!	write(32,*) 'set term x11'	! linux X11

	close(32)
	call system(GNUPLOT//' databar.gnu ')

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

	integer function date2ndx(datum,leap)
	implicit none 
	! returns an index for the given datum by counting days 
	! leap=1 take leap years into account 
	! leap=0 do NOT count leap years: each year 365 days 
	integer dayzero1(12)	! day zero in each month 
	integer iy,day,mnth,datum
	integer dayzero2(12)	! day zero in each month 
	integer dayzero3(12)	! day zero in each month 
	integer nleap 
	integer leap 
	data dayzero1 /0,31,59,90,120,151,181,212,243,273,304,334/ 	! ordinary year 
	data dayzero2 /0,31,60,91,121,152,182,213,244,274,305,335/ 	! leapyear 
	data dayzero3 /0,31,60,91,121,152,182,212,243,273,304,334/ 	! shifted 
!	iy=datum/10000 - 1900				! 1900 is not a leap year 
	iy=datum/10000 					! 1900 is not a leap year 
	nleap=(iy-1)/4					! number of leap years passed since 1900
	mnth=mod(datum,10000)/100
	day=mod(datum,100)
	if(leap.gt.0) then 
C	   if(mod(iy,4).eq.0 .and. mod(iy,100).ne.0) then
	   if(mod(iy,4).eq.0) then
	      date2ndx=dayzero2(mnth)+day+365*(iy)+nleap	! leap year 
	   else 
	      date2ndx=dayzero1(mnth)+day+365*(iy)+nleap	! ordinary year 
	   endif 
	else 
	   if(mod(iy,4).eq.0 .and. mod(iy,100).ne.0) then
	      date2ndx=dayzero3(mnth)+day+365*(iy)	! leap year 
	   else 
	      date2ndx=dayzero1(mnth)+day+365*(iy)	! ordinary year 
	   endif 
	endif 

	return 
	end 

	
