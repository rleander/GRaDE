	Program backtran
C	Purpose : generates the actual synthetic sequences from the index-file 
C	          produced by RESAMPLE.F. 
C	Syntax  : BACKTRAN simname 
C	Input   : simname.log 
C	          .smr and .smt files created by FVECTORS.F 
C	          STDIN (list of filenames, the same list that was used for FVECTORS.F 
C	          can be used) 
C	Output  : For each filename in "list of filenames" a file is returned, 
C	          with the same filename with the name of the simulation (simname) 
C	          added. E.g. if "list of filenames" contains "debilt.tm" and 
C	          "debilt.rr" and simname is "test0", the files "debilt_test0.tm" 
C	          and "debilt_test0.rr" are created. (Note, the number of files 
C	          created during the process is twice the original number of files.)
C
C
C	build with : g77 backtran.f -o backtran -Wall
C
C
C	(c) KNMI, De Bilt, January 2006
C	Robert Leander

	implicit none 
	integer MAXYR,YR0
	parameter (MAXYR=200)	! dimension of the feature-vector array in years 
	parameter (YR0=1900)	! offset of historical years 
        !---------------------------------------------------------------------------

	character*(80) readstr 
	integer instr,last,datum
	integer iarg			
	character*(50) fn_sim,fn_his,fn_out
	real value
	real dayseries(MAXYR*365)	! daily series loaded into memory 
                              
	real shift(365)		! calendarday shift 
	real scale(365)		! calendarday scale 
				! 365 days a year, July 31 is excluded in leapyears 
	integer date2cal,cal2date 
	integer id, iy, calday  
	integer date 
	integer datums(MAXYR*365)
	real a_h,a_s,b_h,b_s
	integer iy_sample 	! year to sample from 
	integer id_sample 	! calendar day to sample from 
				! four feature vector elements 
	character*(15) my_name	! program name 

        logical rainfall, temperature 

! ***** PROCESSING COMMAND-LINE  *************************************
        call getarg(0,my_name)		! What's my name ? 
	if(iarg().lt.1) then 
	   write(0,*) 'syntax : '//my_name(1:instr(my_name,' ')-1)
     &	          //' <simname> '
	   write(0,*) 'names of historical files are read from STDIN'
	   write(0,*) 
	
	   stop ' '
	else 
	   call getarg(1,fn_sim)	! simulation name 
        endif 

	open(12,file=fn_sim(1:instr(fn_sim,' ')-1)//'.log'
     &	        ,status='OLD')
	write(0,*) fn_sim(1:instr(fn_sim,' ')-1)
     &	         //'.log opened as index file' 
	write(0,*) 

	do while (.TRUE.)			! loop over files read from STDIN 
	read(*,*,end=666) fn_his			

! ***** PROCESSING HISTORICAL DATA  *********************************
	open(11,file=fn_his(1:instr(fn_his,' ')-1),status='OLD')
	write(0,*) fn_his(1:instr(fn_his,' ')-1)
     &	     //' opened as historic file' 

	   do id=1,MAXYR*365
	         dayseries(id)=-999.99	! missing values by default
	   enddo

	   do while (.TRUE.)
 211	      read(11,'(a40)',end=111,err=211) readstr
	      if(instr(readstr,'#').le.0) then 		! no comment line 
	         read(readstr,*) datum,value		
	         iy=datum/10000-YR0			
	         if(iy.lt.1 .or. iy.gt.MAXYR) goto 211
	         calday=date2cal(datum,1)		! ndx 1..365 corresponding to calendarday 
	         datums(calday+(iy-1)*365)=datum 	! store date for later use (optional):
							! it can be printed lateron, if one desires 
							! to inspect the seqeuence of historic dates
							! in the simulation.	
	         dayseries(calday+(iy-1)*365)=value	! store value for later use
	      endif					! no comment line  
	   enddo 
 111	   continue 
	   close(11)		! close historical file 

           rainfall=.False.
           temperature=.False. 
           
	   if(instr(fn_his,'.rr').gt.0 .or. instr(fn_his,'.pr').gt.0) then 
	      open(13,file=fn_his(1:last(fn_his,'.'))//'smr'	! rainfall variable 
     &	          ,status='OLD')				! smoothed calendarday statistics 
	      write(0,*) fn_his(1:last(fn_his,'.'))//'smr'
     &	                  //' opened for '//trim(fn_his) 
              rainfall=.True. 
	   endif  
	   if(instr(fn_his,'.t').gt.0) then 
	      open(13,file=fn_his(1:last(fn_his,'.'))//'smt'	! temperature variable 
     &	          ,status='OLD')				! smoothed calendarday statistics 
	      write(0,*) fn_his(1:last(fn_his,'.'))//'smt'
     &	                  //' opened for '//trim(fn_his) 
              temperature=.True.
	   endif  
	   if (.not.(rainfall .or. temperature)) then 
		write(0,*) fn_his(1:instr(fn_his,' '))//
     &		      'is neither rainfall nor temperature'
		      stop '' 
	   endif 
	   
	   do calday=1,365 					! read calendar-day shift and scale
	      read(13,*) id,shift(calday),scale(calday)		! needed for standardization and 
	   enddo  						! de-standardization
	   close(13)


	! simple output name, including only '_sim'
C	fn_out=fn_his(1:instr(fn_his,'.')-1)//'_sim'
C     &	   //fn_his(instr(fn_his,'.'):instr(fn_his,' ')-1)


	! the simulation name is included in the output filename
	fn_out=fn_his(1:instr(fn_his,'.')-1)
     &	   //'_'//fn_sim(1:instr(fn_sim,' ')-1)
     &	   //fn_his(instr(fn_his,'.'):instr(fn_his,' ')-1)
	open(14,file=fn_out(1:instr(fn_out,' ')-1),status='REPLACE')
	write(0,*) fn_out(1:instr(fn_out,' ')-1)//' opened for output' 
	write(0,*) 

! ################ BEGIN BaCKTRANSFORMATTION #########################
	   do while (.TRUE.)					! reading lines from ndx file 			
 	     read(12,'(a80)',end=112) readstr 
	     if(instr(readstr,'#').le.0) then 	! no comment line
                read(readstr,*) iy,calday,date	! parse line from simulation ndx-file 

	        id_sample=mod(date,1000)	! historical year 
	        iy_sample=date/1000-YR0		! historical calendar day	

		value=dayseries(id_sample+(iy_sample-1)*365)
	        if(value.lt.-900) then 				! missing value 
	           write(14,'(i10,f8.2,i10)') 
     &	               cal2date(iy,calday), value		! do not scale/shift
C     &	               ,datums(id_sample+(iy_sample-1)*365)
	        else 
	        a_h=shift(id_sample)            
	        b_h=scale(id_sample)
	        a_s=shift(mod(calday-1,365)+1)     ! day in the simulation can be a leapday
	        b_s=scale(mod(calday-1,365)+1)     ! NB. if calday=366, choose the value for the first calday, i.e. Jan 1st
	     
	        write(14,'(i10,f8.2,i10)') 
     &	           cal2date(iy,calday),
     &	           a_s+(b_s/b_h)*(value-a_h) 

						! if xs is the simulated value and 
						! xh is the value extracted from the historic 
						! file then 
						!
						!      xh-shift[h]
						! xs = ----------- * scale[s] + shift[s]  ,
						!        scale[h]
						!
						! where h represents the historical calendar day 
						! and s de calendar day of the simulated value. 
						! standardization and backtransformation of 
						! historical values is performed in a single 
						! operation.
C     &	           ,datums(id_sample+(iy_sample-1)*365)
	       endif ! not a missing value 
 	     endif ! no comment line 
	   enddo ! reading lines from ndx file 					

! ################ END BaCKTRANSFORMATTION ##########################
 112	   continue 
	close(14)	! close output file 
	rewind(12)	! close ndx-file 
	enddo		! loop over files 
 666	continue 	
	write(0,*) 'Done ! '
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
	         goto 676
	      endif
	   enddo
 676	   continue 
	endif 
	return 
	end 

	integer function last(s,ssub)
	! same as instr, but returns the last position of substring ssub in string s 
	! 0 if not found 
	implicit none 
	integer j,ls,lsub
	character*(*) s 
	character*(*) ssub 
	ls=len(s)
	lsub=len(ssub)
	last=0
	if(lsub.le.ls) then 
	   last=0
	   do j=ls-lsub+1,1,-1
	      if(s(j:j+lsub-1).eq.ssub) then 
	         last=j
	         goto 213
	      endif
	   enddo
 213	   continue 
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

	integer function cal2date(iy,calday)
	implicit none 
	! converts year and calendarday to date YYYYMMDD,
	! taking into account leap year if mod(iy,4)=0 
	integer iy,calday,im,id   
	integer leap(12),no_leap(12)
	data leap /1,32,61,92,122,153,183,214,245,275,306,336/
	data no_leap /1,32,60,91,121,152,182,213,244,274,305,335/
	if(mod(iy,4).eq.0) then 
	   do im=12,1,-1 
	      if(calday.ge.leap(im)) goto 254
	   enddo 
 254	   continue 
	   id=calday-leap(im)+1 
	else 
	   do im=12,1,-1 
	      if(calday.ge.no_leap(im)) goto 260 
	   enddo 
 260	   continue
	   id=calday-no_leap(im)+1 
	endif 
	cal2date=iy*10000+im*100+id
	return 
	end 
	


	integer function iarg()	
	! returns the number of non-empty cmdl args 
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




