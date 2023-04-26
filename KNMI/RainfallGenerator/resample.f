        program RESAMPLE
C	Purpose : Unconditional Daily NN-resampling in feature space, using a 
C	          calendar day window of 121 days [today-60 days, today+60 days]. 
C	Syntax  : RESAMPLE.F   simname   fvector-file   nyears  [random-number-seed] 
C	          The file fvector-file contains the feature-vector information produced 
C	          by FVECTORS.F. The random number seed is optional and is included to 
C	          produce multiple independend sequences with the same settings. 
C	Input   : fvector-file 
C	          nyears (number of simulated years) 
C	          random-number-seed (integer;optional) 
C	Output  : simname.log (resampled historical dates)
C
C
C	build with : g77 resample.f -o resample Wall
C
C	
C	(c) KNMI, De Bilt, January 2006
C	Robert Leander

C	changes:
C	16022006 - in the calculation of the historical memory element, the index 
C		   runs from firstday to lastday, which should be firstday+1 
C	           This error has a crucial impact on the calculation of the 
C	           historical memory element 	
C	           (originally line 174)



      implicit none 
      integer winwidth,nhist,nbuf,nk,ndim,MAXYR,YR0  

      !	PARAMETERS WHICH CAN BE ALTERED TO MODIFY THE SIMULATION
      !------------------------------------------------------------------
      parameter (winwidth=60)	! half-window width (days)				
      			! window is chosen symmetrically from 		
      			! today-winwidth : today+winwidth 		
      parameter (nhist=1)	! 5-day memory (which seems to work best for  
      			! the Meuse basin. For the Rhine basin no  
      			! memory element was included			
      parameter (nk=10)	! fixed number of nearest neighbours 	

      !	PARAMETERS WHICH SHOULD ONLY BE CHANGED FOR SPECIFIC REASONS 
      !------------------------------------------------------------------
      parameter (nbuf=nhist+2)! maximum length of 'memory buffer' (NEVER CHANGE!) 	
      parameter (ndim=4)	! maximum number of dimensions in the feature 
                                ! vector, change ndim in sub 'SORTNNB' 
      			! accordingly!     
      			! Note that the actual choice of the feature 
      			! vector elements is done within SORTNNB
      			! only change whenever feature vector elements 
      			! are customized
      parameter (MAXYR=200)	! dimension of the feature-vector array in years 
      parameter (YR0=1900)	! offset of historical years 
        !---------------------------------------------------------------------------

      integer zaad, iseed	! initial RNG-seed and actual iseed variable
      integer count,i
      integer buf_ptr		! length of and pointer into cyclic buffer where 
      			! previous indices are stored 
      integer buffer(0:nbuf-1)! storage for previously drawn indices (memory)
      real wt(ndim)		! weight of dimensions in the weighted Euclidean distance 
      real sum,sumsqr		! temporary vars for calculating the weights 		
      integer ny 		! simulation length (command line argument) in years 
      integer iarg			
      
      integer instr 
      integer nmod 		! alternative integer modulo function:
      			! nmod(a,b) takes a value between 1 and b
      integer modmod		! alternative integer modulo function to deal with 
      			! negative numbers
      character*(100) fnsim,fndata	
      character*(10) argstr 
      character*(30) myname 
      character*(5) yearstr 

      integer firstyr,lastyr 		! first and last historical year in the simulation
      integer firstday,lastday	! first and last historical day in the simulation
      common /years/ firstyr,lastyr 

      integer Indices(-winwidth:winwidth,MAXYR)
      integer Ir				! random choices 
      integer ndx				! saved indices 
      integer irng,uniran
      real    Fspace(ndim,MAXYR*365)		! historical vectors in feature-space
      common /FeatureSpace/ Fspace
      real	Current(ndim)			! Feature-vector of the current day
      integer dates(MAXYR*365)
      integer nsum 
      real	runsum
      integer year,calday,my_calday 		! year and calendar day read from file 
      real 	vector(ndim)			! feature vector read from file 


      integer nd(24)
      integer ip, im, id, iy, jd, jy 
    	real distr(nk)
      integer sIperm(nk)
      
      real memsum 

      data nd /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31
     &	        ,31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /
      ! First row is NO LEAPYEAR
      ! Second row is LEAPYEAR
      ! So the the right number of days in month im in year iy is 
      ! nd(max(1-mod(iy,4),0)*12+im) taking leap years into account 

      ! cover for no argument problem 

! ******PROCESS CMDL ARGUMENTS ***************************************
        call getarg(0,myname)		! What's my name ? 

      if(iarg().lt.3) then 			! require at least 3 arguments 
           write(0,*) 
           write(0,*) 'syntax : '//myname(1:instr(myname,' ')-1)
     &	           //' <simname> <fvector-file> <nyears>'
     &	             // ' [rng-seed]       ' 
           write(0,*) 
           stop '###'
      endif 

      call getarg(1,fnsim)		! name of the simulation 
      call getarg(2,fndata)		! name of summary statistics-file 
      call getarg(3,argstr)		! number of years to simulate 	
      read(argstr,*) ny 	

      call getarg(0,myname)		! program name 


! ******READ FEATURE-VECTOR DATA IN ELEMENTS 1:NDIM-1 ****************

      do id=1,365*100
         do i=1,ndim  
            Fspace(i,id)=-999.99		
         enddo ! i 
      enddo ! id 

      open(16,file=trim(fndata),status='OLD')
      lastyr=0
      firstyr=99999

      count=0
      do i=1,ndim 
         vector(i)=0.0 
      enddo 
      do while (.TRUE.)
         read(16,*,end=116) year,calday,(vector(i),i=1,ndim-1)     
         if(year.lt.firstyr) firstyr=year
         if(year.gt.lastyr) lastyr=year

         id=(year-YR0-1)*365+calday
         do i=1,ndim-1
            Fspace(i,id)=vector(i)
         enddo ! i 
         dates(id)=year*1000+calday		! YYYYDDD
         count=count+1			! count days read from file 
      enddo 
 116	continue 
      close(16)
      write(0,*)  

C	Should you desire to resample a sub-period within the period spanned by the 
C	feature vectors, please re-specify firstyr and lastyr below accordingly 
C	firstyr = 
C	lastyr = 
      
      write(0,*) 'first year : ',firstyr 
      write(0,*) 'last year  : ',lastyr  

      firstyr = firstyr-YR0		! an offset is subtracted, corresponding to the 
      lastyr = lastyr-YR0		! declaration of arrays 

      write(0,*) 'Done reading feature-vector components' 

! ******CALCULATE THE MEMORY OF THE 1st ELEMENT AS THE NDIM-th ELEMENT
      firstday=(firstyr-1)*365+1
      lastday=lastyr*365

      Fspace(ndim,firstday)=-999.99999
      memsum=0
C	do id=firstday,lastday		! NOTE This error of the first release was corrected
      				! 16-02-2006 (see also program header) 
      do id=firstday+1,lastday
         if(id-1-nhist.ge.firstday) then		! if day subtracted beyond first day 
            memsum=memsum-Fspace(1,id-1-nhist)+Fspace(1,id-1)
            Fspace(ndim,id)=memsum 
         else 
            memsum=memsum+Fspace(1,id-1)
            Fspace(ndim,id)=-999.99				 
         endif
      enddo ! id 
      

! ******MARK MISSING DAYS IN FEATURE-VECTOR SPACE AS NON-ELECTABLE ***

      ! fill Fspace(ndim,1:nhist) with 99999.9 to avoid them being selected 
      do count=1,nhist
         Fspace(ndim,count)=-99999.9
      enddo 

! ******OPEN OUTPUT FOR RESAMPLED INDICES ***************************
      open(7
     &	   ,file=fnsim(1:instr(fnsim,' ')-1)//'.log'
     &	   ,status='UNKNOWN')
      write(0,*) fnsim(1:instr(fnsim,' ')-1)//'.log opened'

! ******INITIALIZE RANDOM NUMBER SEED *******************************
      if(iarg().ge.4) then 
         call getarg(4,argstr)
         read(argstr,*) zaad
         zaad=-abs(zaad)	! initial seed must be negative
      else 
         zaad = -21933
         write(0,*) ' no seed entered, defaulting to ',zaad 
      endif 

! ******CALCULATE WEIGHTS IN THE EUCLIDEAN DISTANCE (INVERSE VARIANCE) 
      write(0,*) 'globally determined weights: ' 
      do ip=1,ndim 
         sumsqr=0.0
         sum=0.0
         nsum=0 
         do iy=firstyr,lastyr
            do id=1,365
               i=(iy-1)*365+id 
               if(abs(Fspace(ip,i)).lt.900) then  
                  sum=sum+Fspace(ip,i)
                  sumsqr=sumsqr+Fspace(ip,i)**2.0
                  nsum=nsum+1
               endif 
            enddo ! id 
         enddo !iy
         sumsqr=(sumsqr-(sum**2.)/nsum)/(nsum-1)		! variance 
         if(sumsqr.gt.0) then 
            wt(ip)=1./sumsqr
         else 
            wt(ip)=0.0
            write(0,*) ip,':', wt(ip),'  Variance Invalid !! ' 
         endif 
      enddo ! ip 
      write(0,*) 

C	Alternatives which are not implemented here are e.g. locally
C	determined weights, i.e. the weights are determined for each 
C	calendar day separately, or "adaptive" weights by means of the
C	Mahalanobis distance  (see KNMI Publication 186-IV).

      do ip=1,ndim 
         write(0,*) ip,':', wt(ip)
      enddo 

      iseed = zaad		! negative seed 
      
      !############################################# 
      write(0,'(a,i6,a)') ' simulation  : ',ny,' years' 
      write(0,'(a,i6,a)') ' memory      : ',nhist,' days' 
      write(0,*) 
      !############################################# 

! ******SELECT 0th HISTORICAL DAY FROM A WINDOW CENTRED ON DEC 31 ****
 
      count = 0				! initialize counter 
      buf_ptr=0				! initialize buffer pointer 

      iy=uniran(firstyr+1,lastyr,iseed)	! pick random year 
      					! the first year is excluded to
      					! avoid problems with 
      					! undefined memory elements

      calday=nmod(uniran(-winwidth,winwidth,iseed),365)	! pick random calendar day 
      ndx = (iy-1)*365+calday					

      write(0,'(a,i7)') ' seed           = ',zaad
      write(0,'(a,i7,a,i4,a,i5)') 
     &	                  ' starting index = ',ndx,'  day' 
     &	             ,mod(ndx-1,365)+1,' of year',(ndx-1)/365+1
      write(0,*)

! ******INITIALIZE MEMORY BUFFER ************************************
      
      buffer(buf_ptr) = ndx			! initialize memory buffer 
      do id=-nhist-1,-1			
         buffer(modmod(buf_ptr+id,nbuf))=buffer(buf_ptr)+id		
      enddo 

! ***** INITIALIZE DECREASING KERNEL *********************************
        do i=1, nk
           distr(i) = 1./real(i)		! P(x) = 1/x
        enddo
      
      write(7,*) '# program : ',myname	! write additional info to .log-file 
      write(7,*) '# Feature vectors :'
     &	       //fndata(1:instr(fndata,' ')-1)
      write(7,*) '# starting index        :',ndx,'  ( day'
     &	   ,mod(ndx-1,365)+1,' of year',(ndx-1)/365+1
      write(7,*) '# seed                  :',zaad
      write(7,*) '# window                :',2*winwidth+1,' days'
      write(7,*) '# memory                :',nhist,' days'
      write(7,*) '# number of nearest nb  :',nk

! ***** BEGIN RESAMPLING ********************************************
      do iy = 1, ny					! For all years in the simulation ...
         write(0,'(a,i10)') ' '// fnsim(1:instr(fnsim,' ')-1)
     &	                        // ', year', iy
         write(yearstr,'(i5)') iy			! For all months in the year .... 
         my_calday=0					! calendar day being simulated 
         do im = 1, 12		 		! For all days in the month ...
            do id=1, nd(max(1-mod(iy,4),0)*12+im)
      	 count =  count + 1			! counts simulated days
               buf_ptr=modmod(buf_ptr+1,nbuf)		! update buffer pointer 	
      	 my_calday = my_calday+1		! calendar-day 1..366 
      	 calday = nmod(my_calday,365)		! simulation calendar-day 1..365
      						! in leap years an additional 1st of 
      						! January is simulated at the end of 
      						! the year
      	
      	! buf_ptr refers to the day currently to be simulated P(i+1)
      	! buf_ptr-1 refers to most recently simulated day P(i)
      	! buf_ptr-2 refers to the day before the most recently simulated day,
      	!    which is the day to be added to the memory (P(i-1))
              ! buf_ptr-2-nhist is P(i-nhist-1)
              !    which is the day to be removed from the memory,
      	!    since memory is the sum of P(i-1)...P(i-nhist)
      	!    equal to bufptr, since the buffer length=nhist+2

              if(count-2.ge.1) then 				! adjust running n-day sum 
                 runsum=runsum
     &	               +Fspace(1,buffer(modmod(buf_ptr-2,nbuf)))	
                 if(count-nhist-2.ge.1) then 
                   runsum=runsum
     &	               -Fspace(1,buffer(buf_ptr))	! last simulation step (count-1) 
     	               					! note that 
      						! modmod(buf_ptr-nhist-2,nbuf)=buf_ptr	
                   runsum=max(runsum,0.0)
                 endif 
              else 
                 runsum=0.0
              endif 

      	do jy=firstyr,lastyr			! fill window with indices 
         	   do jd=-winwidth,winwidth 
            	      Indices(jd,jy)=nmod((calday-1)+jd,365)
     &	                      +(jy-1)*365
         	   enddo 
      	enddo 
      	! The window must be centred on the last simulated calendar day, 
                ! for which an analog is searched.

              do i=1,ndim-1           
                 Current(i)=Fspace(i,ndx)	
              enddo 
              Current(ndim)=runsum 		! is the memory of the first 
      					! feature-vector element 

                ! feature-vector of the last simulated day is used for
                ! comparison

         	Ir=irng(distr,1,nk,iseed)	! random nnb choice 

    		call sortnnb(sIperm, Indices, Current ,wt, Ir
     &                             ,winwidth, firstyr, lastyr)	

      	ndx = sIperm(Ir) + 1		! new ndx is successor of selected nearest nb

              buffer(buf_ptr)=ndx


! **** WRITE INDEX TO FILE *******************************************
    		write(7,'(i6,i4,i8,i3)') 
     &	                             iy		! year in simulation (1..ny)
     &	                           , my_calday	! calendar day in simulation (1..366) 
     &	                           , dates(ndx)	! date of the selected historical day 
      					! year*1000+my_calday
     &	                           , Ir		! selected nearest neighbour rank 	
 
            enddo ! id 
         enddo ! im 
    	enddo ! iy

! ***** END RESAMPLING **********************************************

        ! write last sum to complete index file 
      runsum=runsum
     &	    +Fspace(1,buffer(modmod(buf_ptr-2,nbuf)))		
     &	    -Fspace(1,buffer(modmod(buf_ptr-nhist-2,nbuf)))	
      runsum=max(runsum,0.0)
      close(7)

      write(0,*) fnsim(1:instr(fnsim,' ')-1) //' done !!' 
      stop 
      end





!============================================================================
    	SUBROUTINE sortnnb( 
     &				si, 		! indices of sorted nearest neighbours (out)
     &				Indices, 	! indices of days within the window
     &				Current,	! feature vector of the 'current' state 
     &				invvar,		! weights of feature-vector elements
      					! (=1/variance)
     &				ichoice, 	! randomly selected nearest neighbour 
     &				ww, 		! half-window width  
     &				y0, y1		! first and last year (in the array fspace) to be searched
     &			)
    	implicit none

        integer ndim 
      !---------------------------------------------------------------------------
        parameter (ndim=4)	! maximum number of dimensions in the feature vector !
      			! NOTE ndim should correspond with ndim in the main 
      			! program.
      !---------------------------------------------------------------------------
        
      integer MAXYR	
      parameter(MAXYR=200)

      integer si(*),ichoice
      integer ww 			! window width 
      integer Indices(MAXYR*(2*ww+1))	! window with indices 
      real    Fspace(ndim,MAXYR*365)	! historical vectors in feature-space
      common /FeatureSpace/ Fspace

      real Current(ndim)
      real invvar(ndim) 	
      integer i, j
      integer firstndx,lastndx  
      real dx, sa(0:50)		! The distances of Max. 50 nearest-neighbours 
      				! can be stored here  
      integer y0,y1
      
      ! Each row in Indices had a length of -ww..0..ww = 2*ww+1
      ! this is the first (fastest) index 
      firstndx = (y0-1)*(ww*2+1)+1 	! first considered entry of Indices 
      lastndx = y1*(ww*2+1)		! last considered entry of Indices 

      do i=1,ichoice 		
         sa(i)=99999.9	! any distance is below 99999.9
      enddo ! i 
      sa(0) = -1		! any distance is beyond sa(0)	
      
      ! check the remaining days in the window for closer 
      do j=firstndx, lastndx		! calculate the Euclidean distances and gather the 
      				! ichoice smallest in sa
               dx=0.0
               dx=dx+invvar(1)				! average standardized rainfall 
     &	               *(Fspace(1,Indices(j))-Current(1))**2.	
               dx=dx+invvar(2)				! average standardized temperature 
     &	               *(Fspace(2,Indices(j))-Current(2))**2.	

C Note that this code resembles the code for the Meuse basin and therefore the fraction 
C of wet stations is not used in the search for nearest neighbours
C	         dx=dx+invvar(3)				! fraction of wet stations
C     &	               *(Fspace(3,Indices(j))-Current(3))**2.	

               dx=dx+invvar(4)				! 5-day memory of average standardized rainfall
     &	               *(Fspace(4,Indices(j))-Current(4))**2.	

      	 if(Indices(j).eq.y1*365) dx=9999.9		! exclude last historical day

               if (dx .lt. sa(ichoice))  then 
                  i=ichoice-1
      	    do while(sa(i).gt.dx)
      	       sa(i+1) = sa(i)
      	       si(i+1) = si(i)
      	       i=i-1
      	    enddo 
      	    sa(i+1) = dx
      	    si(i+1) = Indices(j)  
               endif 	!... if distance is shorter than largest in list 
      enddo  ! loop over indices 

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
               goto 676
            endif
         enddo
 676	   continue 
      endif 
      return 
      end 

C  (C) Copr. 1986-92 Numerical Recipes Software -1$393.


      integer function irng(distrib,min,max,iseed)
      implicit none 
      ! returns a random int from min to max
      ! 
      integer min,max,iseed,i
      real distrib(min:max),sum,uni	
      real ran3

      sum = 0.0
      do i=min,max 
         sum=sum+distrib(i) 
      enddo 
      uni=ran3(iseed)*sum	! draw from uniform distribution 
      			! between 0 and sum(distrib) 
      			! now find appropriate interval 
      i=min
      sum=sum-distrib(i)
      do while (uni.le.sum)
         i=i+1
         sum=sum-distrib(i)
      enddo 
      irng=i
      return 
      end 



      integer function nmod(x,y)
      ! returns x mod y as a number in the range 1..y
      implicit none 
      integer x,y,modmod
      nmod=modmod(x-1,y)+1
      end 


      integer function uniran(i1,i2,ISEED)
      implicit none 
      ! uses the RNG to return an int ranging from i1 to i2 
      ! (that is: not larger than i2 and not smaller than i1
      ! integer ISEED is passed by ref. (input/output)
      ! uniform distribution 
  
      integer ISEED,i1,i2
      real ran3
      real xran 
      xran=ran3(ISEED)		! call RNG with seed
      				! returning a number between 0 and 1 
      uniran=i1+nint(xran*(i2-i1))
      return 
      end 

      integer function modmod(a,b)
      implicit none 
      integer a,b 
      modmod=mod(mod(a,b)+b,b)
      return 
      end 

      
      FUNCTION ran3(idum)
      INTEGER idum
      INTEGER MBIG,MSEED,MZ
C     REAL MBIG,MSEED,MZ
      REAL ran3,FAC
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
C     PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
      INTEGER i,iff,ii,inext,inextp,k
      INTEGER mj,mk,ma(55)
C     REAL mj,mk,ma(55)
      SAVE iff,inext,inextp,ma
      DATA iff /0/
      if(idum.lt.0.or.iff.eq.0)then
        iff=1
        mj=MSEED-iabs(idum)
        mj=mod(mj,MBIG)
        ma(55)=mj
        mk=1
        do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.MZ)mk=mk+MBIG
          mj=ma(ii)
11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
12        continue
13      continue
        inext=0
        inextp=31
        idum=1
      endif
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.MZ)mj=mj+MBIG
      ma(inext)=mj
      ran3=mj*FAC
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software v%1jw#<0(9p#3.
      
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

      
      


