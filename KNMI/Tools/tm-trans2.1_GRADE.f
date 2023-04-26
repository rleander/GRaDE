	program tmtrans21
	implicit none 

!JB     modified by:
!JB     Jules Beersma (KNMI) 
!JB     to transform 1000-yr parts of 50K-yr GRADE simulations
!JB     See !JB and cJB
!JB     20140924        

!JB     tm-multisite (versions) originally created by Alexander Bakker and Carine Homan in 2010
!JB     basis for this program was: tm-multisite-tmp2.f(.gz) dd 20101012  

	! modification with regard to tmtransmonth: NO MORE SCRATCHING !!              !JB???
	! using static arr (it's still good ol' f77, so no pointers and stuff)

	! can be used to transform an observed time series in a future one
	! with the help of the 10% quantile (Q10),
	! the 50% quantile (Q50)
	! and the 90% quantile (Q90)
	! for each of the calendar month

	! syntax : tm-trans-month <scenario> input(STDIN) output(STDOUT)
	! scenario can be one of 'G','W','G+','W+'
	! For the precipitation scenario 'G+' e.g. a file with the relative changes 
	! named 'tm_G+.txt' (and similar for the other scenarios) should be present, which 
	! is expected to have 7 cols like:
	!
        !1   1.320   2.750   6.240  12.570   5.580  11.160 
	!2   0.690   1.480   5.880  11.860   5.610  11.220 
	!3  -0.990  -1.920   4.940   9.980   5.710  11.400
	! ......
        !12 -4.350  -8.700   3.050   6.200   5.900  11.750 
	!
	! containing monthnr, dQ10(2050,2100), dQ50(2050,2100) and dQ90(2050,2100)
	! i.e. for each of the three stats the relative changes for 2050 and 2100.
	! Changes default to point zero	and are dumped to STDERR for check 
	! Piecewise-linear interpolation of the relative changes is based on the given year.
	! Changes are dumped to STDERR for check 

	! The daily data input should be in the format 
	! date1 value1 
	! date2 value2
	!
	! where the date may be formatted as either 
	! YYYYMMDD
cJB	! YYYYMMDDhh	(hh is ignored) 
cJB	! YYYY MM DD    or 
cJB	! DD MM YYYY
cJB	! Any other format will cause the program to halt

	! Missing values, either absent or invalid on a 'date value'-line, or below -90
	! are ignored in determining the transformation parameters and result in a value of -99.99
	! in the output.

	! The output is always in the format YYYYMMDD value 
	! 
	! syntax to compile: "f77 tm-multisite_GRADE.f" 
	

	integer MAXDAYS,tsmax,ts
	parameter(MAXDAYS=1000*365+250)		! a century 'll do - chgd 1000yrs 
	parameter(tsmax=134)
	integer dd(MAXDAYS)			! date 
	real*8 tm(MAXDAYS,tsmax)		! and precipitation amount
cJB	character*(10000) regel 
	character*(5000) regel 
	character*(5) scenario			! name of the scenario: G,W,G+,W+
cJB	character*(*) fmt
cJB	parameter(fmt='7.1')			! format specifier (string) for output 
						! in the form 'spaces.decimals'
cJB	integer cyear_obs, shift_y
	integer im,ivar,instr,j,date,i10alpha,i50alpha,i90alpha
	integer is
	integer lastchar
	integer iarg, iday 
	integer cyear				! the central year, this determines the 
						! relative changes for the scenario.
						! it is assumed that the subjected series corresponds to the 
						! control year 1990.
						! The applied changes are obtained by piecewise linear 
						! interpolation between 1990, 2050 and 2100.

	real*8 alpha				! weighting factor between scenarios
	real*8 value(tsmax)			! temporary temperature amount
	real*8 chg(12,3)			! relative changes [%], (month,var) 
	real*8 chg2050(12,3),chg2100(12,3)
	real*8 months(0:31*1000,12,tsmax) 
	real*8 Q10_o(12,tsmax),Q10_f(12,tsmax)	! Q10 observed and future
	real*8 Q50_o(12,tsmax),Q50_f(12,tsmax)	! Q50 observed and future
	real*8 Q90_o(12,tsmax),Q90_f(12,tsmax)	! Q90 observed and future
	real*8 a(12,tsmax)			! transformation parameter a
	real*8 b(12,tsmax)			! transformation parameter b

	real*8 x90
	parameter (x90=90.d0)			! quantile to be fitted [%]
	real*8 x50
	parameter (x50=50.d0)			! quantile to be fitted [%]
	real*8 x10
	parameter (x10=10.d0)			! quantile to be fitted [%]
	integer cnt(12,tsmax)			! general daycounter per months
	integer cnt_t
        integer ios
        integer argint
        logical argstring, arglogical


	character*3 acroniem(12) 		! acronyms for months 
	data acroniem/'Jan','Feb','Mar','Apr','May','Jun',
     &	              'Jul','Aug','Sep','Oct','Nov','Dec'/

******* CHECK NUMBER OF ARGUMENTS **************************************************************
        if(iarg().lt.2) then
cJB        write(0,*) 'syntax : rr-trans <scenario> <year> -ts' 
           write(0,*) 'syntax : rr-trans -sc scenario -cy year [-n ts]' !JB 20140923
           write(0,*) 'Transforms daily precipitation with changes'
           write(0,*) 'representative of the given central year.'
           write(0,*) 'scenario should be one of:'
           write(0,*) 'G, W, G+ or W+'
           write(0,*) 'Reads from STDIN, writes to STDOUT'
cJB        write(0,*) '-ts: ts is number of time series'
           write(0,*) '-n ts: ts is number of time series '
     &                  //'(default is ts=1)' !JB 20140923
           write(0,*)
           stop '###'
        endif

******* PROCESS COMMAND LINE ARGUMENTS  *****
        ! Get the scenario name from the command line args
        if (.not.argstring('-sc','',scenario)) then
           stop 'No scenario: -sc <scenario>'
        endif
        ! Get the central year from the command line args
        cyear=argint('-cy',-1)
        if (cyear.lt.0) then
           stop 'No valid central year: -cy <year>'
        endif
        if(cyear.lt.1990 .or. cyear.gt.2100) then
           write(0,*) 'central year should be between 1990 and 2100'
           stop '###'
        endif
        ! Get the optional number of series from the command line args, default=1
        ts=argint('-n',1)
        ts=max(min(ts,1000),1)

        open(22,
     &     file='rr_'//trim(scenario)//'.txt',
     &     status='OLD', iostat=ios)
        if (ios.ne.0) then
           print *, ios
           write(0,*) 'Cannot open rr_'//trim(scenario)//'.txt'
           stop
        endif
        do while(.True.)
 222       read(22,'(a100)',end=122) regel
           read(regel,*,end=222,err=222)
     &          im, (chg2050(im,ivar),chg2100(im,ivar),ivar=1,3)
        enddo
 122    continue
        close(22)

	! DETERMINE THE RELATIVE CHANGE FOR EACH MONTH BY 
	! PIECEWISE-LINEAR INTERPOLATION
	
	! cyear should exceed 1990 and fall below 2100 to obtain
	! sane results. Not checked here, users own responsibility .... 

	if(cyear.le.2050) then 
	   alpha=real(cyear-1990)/real(2050-1990)
	   do im=1,12 
	      do ivar=1,3 
	         chg(im,ivar)=alpha*chg2050(im,ivar)
	      enddo ! ivar 
	   enddo ! im 
	else 
	   alpha=real(cyear-2050)/(2100-2050)
	   do im=1,12 
	      do ivar=1,3 
	         chg(im,ivar)=alpha*chg2100(im,ivar)
     &	                     +(1.-alpha)*chg2050(im,ivar)
	      enddo ! ivar 
	   enddo ! im 
	endif 


	! AND WRITE THOSE CHANGES ON THE SCREEN
	write(0,'(a)') 
     & 	       '# month     deltaQ10    deltaQ50     deltaQ90'
	do im=1,12
	   write(0,'(i4,3f14.3,a8)') im,
     &	                  (chg(im,ivar),ivar=1,3),acroniem(im)
	enddo 
	write(0,*)

	! READ THE DATA TO BE TRANSFORMED FROM STDIN, CALCULATE STATISTICS AND 
	! WRITE TO ARRAY tm
	do im=1,12
	   do is=1,ts
	      Q10_o(im,is)=0.0d0
	      Q50_o(im,is)=0.0d0
	      Q90_o(im,is)=0.0d0
	      months(0,im,is)=1E+07
	      cnt(im,is)=0
	   enddo ! is
	enddo ! is
	cnt_t=0

!JB first line is always a header 
        read(*,'(a5000)') regel
        write(*,'(a)') regel(1:lastchar(regel))
!JB

	do while(.True.)
 112	   read(*,'(a5000)',end=666) regel
!JB check if subsequent lines are comment lines (i.e. starting with '#')
	   if(instr(regel,'#').gt.0) then 
	      write(*,'(a)') regel(1:lastchar(regel))
	      goto 112
	   endif
	   
	   read(regel,*) date
	   if(date.eq.0) then
	      write(*,*) regel(2:lastchar(regel))
	   else
	      do is=1,ts
	         value(is)=-99.9				! default value : missing nr -99.9
	      enddo ! is
C	      call parse(regel,date,value,ts)
cJB	      read(regel,*,end=666,err=666) date,
cJB  &                                     (value(is),is=1,ts)
	      read(regel,*) date, (value(is),is=1,ts)
	      cnt_t=cnt_t+1
	      dd(cnt_t)=date

	      do is=1,ts
	         tm(cnt_t,is)=value(is)
	         if(value(is).gt.-90.) then 		!if valid temperature
	            im=mod(date,10000)/100
	            cnt(im,is)=cnt(im,is)+1
	            j = cnt(im,is)
	            do while(value(is).gt.months(j-1,im,is))  ! separate the input time series per season
	               months(j,im,is)=months(j-1,im,is)	   ! and sort them
	               j=j-1
	            enddo					! if valid temperature
	            months(j,im,is)=value(is)
	         endif 
C		 write(0,*) 'is=',is
	      enddo ! is
c	      if (cnt_t .eq. 18263) goto 666
cJB	      if (cnt_t .eq. 365250) goto 666         !JB superfluous: handeled by line with label 112  
	   endif
	enddo 
 666	continue 



*** determine needed shift such that cyear is in middle of transformed series ***	
cJB	cyear_obs=(dd(1)+dd(cnt_t))/20000
cJB	shift_y=(cyear-cyear_obs)*10000

	! DETERMINE THE OBSERVED QUANTILE FOR EACH MONTH
	do im=1,12
	 do is=1,ts
	   i90alpha=int((1.-x90*0.01)*cnt(im,is))
	   Q90_o(im,is)=(months(i90alpha,im,is)+
     &                          months(i90alpha+1,im,is))
     &	                    /(2.+1E-07)
	   i50alpha=int((1.-x50*0.01)*cnt(im,is))
	   Q50_o(im,is)=(months(i50alpha,im,is)+
     &                          months(i50alpha+1,im,is))
     &	                    /(2.+1E-07)
	   i10alpha=int((1.-x10*0.01)*cnt(im,is))
	   Q10_o(im,is)=(months(i10alpha,im,is)+
     &                          months(i10alpha+1,im,is))
     &	                    /(2.+1E-07)
         enddo ! is
	enddo ! im

C	write(0,'(a)')'#month		Q90'
C	do im=1,12
C	 write(0,'(a8,100f10.4)') acroniem(im),
C     &                            (Q90_o(im,is),is=1,ts)
C	enddo
C	write(0,*)
	
C	write(0,'(a)')'#month		Q50'
C	do im=1,12
C	 write(0,'(a8,100f10.4)') acroniem(im),
C     &                            (Q50_o(im,is),is=1,ts)
C	enddo
C	write(0,*)

C	write(0,'(a)')'#month		Q10'
C	do im=1,12
C	 write(0,'(a8,100f10.4)') acroniem(im),
C     &                            (Q10_o(im,is),is=1,ts)
C	enddo
C	write(0,*)


	! DETERMINE THE QUANTILES FOR THE FUTURE CLIMATE 
	do im=1,12
	 do is=1,ts
	   Q10_f(im,is) = Q10_o(im,is)+(chg(im,1))	! 1 = Q10
	   Q50_f(im,is) = Q50_o(im,is)+(chg(im,2))	! 2 = Q50
	   Q90_f(im,is) = Q90_o(im,is)+(chg(im,3))	! 3 = Q90
	 enddo ! is
	enddo ! im

	! CALCULATE a AND b
	do im=1,12
	 do is=1,ts
	   a(im,is) = (Q90_f(im,is) - Q50_f(im,is))/
     &                   (Q90_o(im,is) - Q50_o(im,is))
	   b(im,is) = (Q10_f(im,is) - Q50_f(im,is))/
     &                   (Q10_o(im,is) - Q50_o(im,is))
         enddo
	enddo
	
C	write(0,'(a)')'#month		a'
C	do im=1,12
C	 write(0,'(a8,100f10.4)') acroniem(im),(a(im,is),is=1,ts)
C	enddo
C	write(0,*)
	
C	write(0,'(a)')'#month		b'
C	do im=1,12
C	 write(0,'(a8,100f10.4)') acroniem(im),(b(im,is),is=1,ts)
C	enddo
C	write(0,*)

	
	!CALCULATE NEW DATA FOR THE FUTURE
	do iday=1,cnt_t
	  date=dd(iday)
	  im = mod(date,10000)/100
	  do is=1,ts
	     value(is)=tm(iday,is)
	     if(value(is).gt.-90.) then 				! do not transform missing value code 
	       if(value(is).ge.Q50_o(im,is))then			! -99.99
	value(is)=a(im,is)*(value(is)-Q50_o(im,is))+Q50_f(im,is)
	       elseif (value(is).lt.Q50_o(im,is))then
	value(is)=b(im,is)*(value(is)-Q50_o(im,is))+Q50_f(im,is)
	       endif 
	     endif
	  enddo ! is
C	  write(*,'(i9.8, 2x, 100f'//fmt//')') (date+shift_y),(value(is),is=1,ts)	! write to STDOUT
	  write(*,'(i9.8, 200f6.1)') date,(value(is),is=1,ts)                    	! write to STDOUT
	enddo ! iday
	end


	integer function ndigit(n)
	implicit none 
	integer n 
	ndigit=int(log(n*1.0)/log(10.0))+1
	return 
	end 

	subroutine parse(regel,date,val,ts)
	implicit none 
	integer ndigit,iy,im,id,date,is,i,ts
	character*(2500) regel
	real*8 val(*)
	read(regel,*) date			! try to read the date 
	if(ndigit(date).ge.8) then 
	         read(regel,*,end=197,err=197) date,
     &                                        (val(i),i=1,ts)
 197	         continue 
C     	         date=date/(10**(ndigit(date)-8))
	else
	         if(ndigit(date).ge.4) then 
	             read(regel,*,end=207,err=207) iy,im,id,
     &                                         (val(i),i=1,ts)
	         else 
		     read(regel,*,end=207,err=207) id,im,iy,
     &                                         (val(i),i=1,ts)
     	         endif 
 207	         continue 
	         date=iy*10000+im*100+id
	endif 
	return 
	end 

        include 'routines.f'
