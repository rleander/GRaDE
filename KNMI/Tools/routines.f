C	GLOBAL SUBROUTINES AND FUNCTIONS 

        integer function nextinteger(s,i,delim)
        implicit none
        character*(*) s, delim
        character*(20) ssub
        integer str_from_ndx_to_delim, i
        ! returns next int in string s from pos i upto delimiter 'delim', adjusts pos i
        if(str_from_ndx_to_delim(s,ssub,i,delim).gt.0) then
          read(ssub,*) nextinteger
        endif
        return
        end

        real function nextreal(s,i,delim)
        implicit none
        character*(*) s, delim
        character*(20) ssub
        integer str_from_ndx_to_delim, i
        ! returns next int in string s from pos i upto delimiter 'delim', adjusts pos i
        if(str_from_ndx_to_delim(s,ssub,i,delim).gt.0) then
          read(ssub,*) nextreal
        endif
        return
        end

        integer function str_from_ndx_to_delim(s,ssub,startndx,delim)
        ! Extracts the substring ssub out of s from index startndx
        ! upto the first delimiter or end of s
        ! 0 if not found
        implicit none
        integer startndx,endndx,lastchar,substrndx
        integer lendelim
        character*(*) s,ssub,delim
        lendelim=len(delim)
        endndx=substrndx(s,delim,startndx)
        if(endndx.eq.0) then
          endndx=lastchar(s)
        else
          endndx=endndx-1
        endif
        ssub=s(startndx:endndx)
        str_from_ndx_to_delim=endndx-startndx+1
        startndx=endndx+lendelim+1
        return
        end


        integer function substrndx(s,ssub,startndx)
        ! same as 'index', but start at startndx in the string
        ! it returns the position of substring ssub in string s
        ! 0 if not found
        implicit none
        integer j,ls,lsub,startndx,ndx
        character*(*) s
        character*(*) ssub
        ls=len(s)
        lsub=len(ssub)
        if(lsub.le.ls) then
           ndx=0
           do j=startndx,ls-lsub+1
              if(s(j:j+lsub-1).eq.ssub) then
                 ndx=j
                 goto 022
              endif
           enddo
 022       continue
        endif
        if(ndx.ge.startndx) then
            substrndx=ndx
        else
            substrndx=0
        endif
        return
        end

      integer function splits(delim,subject,strings)
      ! Fill array strings with substrings obtained by splitting
      ! with delimiter delim, return the number of substrings
      ! nb. delim must be a const character*(*), though may have more than one
      implicit none
      integer lastchar, lensubject, lendelim, i
      integer firstndx, j
      character*(*) subject, delim
      character*(*) strings(*)
      lensubject=lastchar(subject)
      lendelim=len(delim)
      firstndx=1
      j=0
      do i=1,lensubject-lendelim+1
         if(subject(i:i+lendelim-1).eq.delim) then
           strings(j+1)=subject(firstndx:i-1)
           j=j+1
           firstndx=i+lendelim
         endif
      enddo
      if(firstndx.le.lensubject) then
        strings(j+1)=subject(firstndx:lensubject)
        j=j+1
      endif
      splits=j
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
 676       continue
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
 213       continue
        endif
        return
        end


	integer function firstchar(s)
	implicit none 
	character*(*) s 
	integer i
	do i=1,len(s)
	   if(s(i:i).ne.' ') then 
	      firstchar=i
	      goto 037 
	   endif 
	enddo ! i
 037	continue 
	return 
	end 

        integer function lastchar(s)
        implicit none
        character*(*) s
        integer i, lst

        do i=len(s),1,-1
           if(s(i:i).ne.' ') then
              lst=i
              goto 047
           endif
        enddo ! i
 047    continue
        do i=0,len(s)-1
           if((ichar(s(i+1:i+1)).lt.32)
     &        .and. (ichar(s(i+1:i+1)).ne.13) 
     &        .and. (ichar(s(i+1:i+1)).ne.10)) then    ! ascii-zero string end
              lst=min(lst,i)
              goto 057
           endif
           if(ichar(s(i+1:i+1)).gt.126) then           ! ascii-zero string end
              lst=min(lst,i)
              goto 057
           endif
        enddo ! i
 057    continue
        lastchar=lst
        return
        end


	integer function firstfind(s,ssub)
	! same as 'index' 
	! it returns the position of substring ssub in string s 
	! 0 if not found 
	implicit none 
	integer j,ls,lsub
	character*(*) s 
	character*(*) ssub 
	ls=len(s)
	lsub=len(ssub)
	firstfind=0
	if(lsub.le.ls) then 
	   firstfind=0
	   do j=1,ls-lsub+1
	      if(s(j:j+lsub-1).eq.ssub) then 
	         firstfind=j
	         goto 070
	      endif
	   enddo
 070	   continue 
	endif 
	return 
	end 

	integer function lastfind(s,ssub)
	! same as 'index' 
	! it returns the position of substring ssub in string s 
	! 0 if not found 
	implicit none 
	integer j,ls,lsub
	character*(*) s 
	character*(*) ssub 
	ls=len(s)
	lsub=len(ssub)
	lastfind=0
	if(lsub.le.ls) then 
	   lastfind=0
	   do j=ls-lsub+1,1,-1
	      if(s(j:j+lsub-1).eq.ssub) then 
	         lastfind=j
	         goto 093
	      endif
	   enddo
 093	   continue 
	endif 
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

        integer function nmod(x,y)
        ! returns x mod y as a number in the range 1..y
        implicit none
        integer x,y,modmod
        nmod=modmod(x-1,y)+1
        end


        integer function modmod(a,b)
        implicit none
        integer a,b
        modmod=mod(mod(a,b)+b,b)
        return
        end


!============================================================================
C retrieving flags or parameters from the command-line:
C argint, argreal, arglogical
!============================================================================

        integer function argint(prefix, default)
        implicit none
        integer lastchar 
        integer jarg
        integer instr, iarg
        integer default, result
        character*(*) prefix
        character*50 sarg

        result=default
        do jarg=1,iarg()
           call getarg(jarg,sarg)
           if(sarg(1:instr(sarg,' ')-1).eq.prefix) then
              call getarg(jarg+1,sarg)
              read(sarg,*,end=233,err=233) result
              write(0,'(2a,i0)') 'ARGINT : ',
     &           prefix(1:lastchar(prefix))//' = ', result ! dump to screen     
           endif
 233       continue
        enddo
        argint=result
        return
        end

        real function argreal(prefix, default)
        implicit none
        integer jarg
        integer instr, iarg
	integer lastchar
        real default, result
        character*(*) prefix
        character*50 sarg

        result=default
        do jarg=1,iarg()
           call getarg(jarg,sarg)
           if(sarg(1:instr(sarg,' ')-1).eq.prefix) then
              call getarg(jarg+1,sarg)
              read(sarg,*,end=323,err=323) result
              write(0,'(2a,e15.5)') 'ARGREAL : ',
     &           prefix(1:lastchar(prefix))//' = ', result ! dump to screen     
           endif
 323       continue
        enddo
        argreal=result
        return
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

        logical function argstring(prefix, default, result) 
        implicit none 
        integer jarg 
        integer instr, iarg  
        character*(*) default 
        character*(*) result            !MAX 20 characters!
        character*(*) prefix  
!       character*50 sarg 
        character*100 sarg 
         
        result=default  
        argstring=.False.
        do jarg=1,iarg() 
           call getarg(jarg,sarg)
           if(sarg(1:instr(sarg,' ')-1).eq.prefix) then 
              call getarg(jarg+1,sarg)
!             read(sarg,*,end=233,err=233) result
              result = trim(sarg)
              argstring=.True.
           endif
 233       continue 
        enddo
        return
        end


!============================================================================

	SUBROUTINE weigths(Space,waits,firstyr,lastyr,NDIM)
C	Calculate weights for the Euclidean distances (inverse variance)
	implicit none 
	integer NDIM 		! dimensionality of the feature-vectors 
	integer firstyr, lastyr ! period 	
	real Space(NDIM,*)	! feature-vectors 
	integer ip,id,iy,i,nsum
	real sum, sumsqr, waits(NDIM)

        do ip=1,NDIM
           sumsqr=0.0
           sum=0.0
           nsum=0
           do iy=firstyr,lastyr
              do id=1,365
                 i=(iy-1)*365+id
                 if(abs(Space(ip,i)).lt.900) then	! values outside [-900,900]
                    sum=sum+Space(ip,i)			! are not taken seriously
                    sumsqr=sumsqr+Space(ip,i)**2.0
                    nsum=nsum+1
                 endif
              enddo ! id
           enddo !iy
           sumsqr=(sumsqr-(sum**2.)/nsum)/(nsum-1)              ! variance
           if(sumsqr.gt.0) then
              waits(ip)=1./sumsqr
           else
              waits(ip)=0.0
           endif
        enddo ! ip

C       Alternatives which are not implemented here are e.g. locally
C       determined weights, i.e. the weights are determined for each
C       calendar day separately, or "adaptive" weights by means of the
C       Mahalanobis distance  (see KNMI Publication 186-IV).

	return 
	end 

!============================================================================
        SUBROUTINE sortnnb_yrs(
     &                          si,             ! indices of sorted nearest neighbours (out)
     &				sa,		! sorted distances to be returned 
     &                          Indices,        ! indices of days within the window
     &                          Fspace,         ! feature vector space
     &                          Current,        ! feature vector of the 'current' state
     &                          invvar,         ! weights of feature-vector elements
                                                ! (=1/variance)
     &                          ichoice,        ! randomly selected nearest neighbour
     &                          ww,             ! half-window width
     &                          ndim,           ! feature-vector dimension
     &                          y0, y1          ! first and last year (in the array fspace) to be searched
     &                  )
        implicit none

        integer si(*),ichoice
        integer ww                      ! window width
        integer ndim                    ! feature-vector dimension
        integer Indices(*)              ! window with indices
                                        ! subdivided in rows with
                                        ! length 2*ww+1
        real Fspace(ndim,*)             ! historical vectors in feature-space

        real Current(ndim)
        real invvar(ndim)
        integer i, j, idim
        integer firstndx,lastndx
        real dx, sa(0:50)               ! The distances of Max. 50 nearest-neighbours
                                        ! can be stored here
        integer y0,y1

        ! Each row in Indices had a length of -ww..0..ww = 2*ww+1
        ! this is the first (fastest) index
        firstndx = (y0-1)*(ww*2+1)+1    ! first considered entry of Indices
        lastndx = y1*(ww*2+1)           ! last considered entry of Indices

        do i=1,ichoice
           sa(i)=1E+14
           si(i)=-1
        enddo ! i
        sa(0) = -1              ! any distance is beyond sa(0)

        ! check the remaining days in the window for closer
        do j=firstndx, lastndx          ! calculate the Euclidean distances and gather the
                                        ! ichoice smallest in sa
                 dx=0.0
                 do idim=1,ndim
                    dx=dx+invvar(idim)*
     &                (Fspace(idim,Indices(j))-Current(idim))**2.
                 enddo

                 if(Indices(j).eq.y1*365) dx=9999.9             ! exclude last historical day
                 if (dx .lt. sa(ichoice))  then
                    i=ichoice-1
                    do while(sa(i).gt.dx)
                       sa(i+1) = sa(i)
                       si(i+1) = si(i)
                       i=i-1
                    enddo
                    sa(i+1) = dx
                    si(i+1) = Indices(j)
                 endif  !... if distance is shorter than largest in list
        enddo  ! loop over indices

        return
        end
!============================================================================

!============================================================================
        SUBROUTINE sortnnb(
     &                          si,             ! indices of sorted nearest neighbours (out)
     &                          sa,             ! sorted distances 
     &                          Indices,        ! indices of days within the window
     &                          Fspace,         ! feature vector space
     &                          Current,        ! feature vector of the 'current' state
     &                          invvar,         ! weights of feature-vector elements
                                                ! (=1/variance)
     &                          ichoice,        ! randomly selected nearest neighbour
     &                          nvar,           ! feature-vector dimension
     &                          nday            ! number of days in Indices
     &                  )
        implicit none
	include 'params.f'

        integer si(*),ichoice
        integer nday                    ! number of days to choose from 
        integer nvar                    ! feature-vector dimension
        integer Indices(*)              ! window with indices
                                        ! subdivided in rows with
                                        ! length 2*ww+1
        real Fspace(nvar,*)             ! historical vectors in feature-space

        real Current(nvar)
        real invvar(nvar)
        integer i, j, idim
        real dx, sa(0:*)      

	real valid
	parameter(valid=-900.)		! lower bound of the first feature-
					! vector element

        do i=1,ichoice
           sa(i)=1E+14
           si(i)=-1
        enddo ! i
        sa(0) = -1      ! any distance is beyond sa(0)

        ! check the days in the window for closest
                        
        do j=1,nday
	   if(Fspace(1,Indices(j)+1).ge.valid) then 
                 dx=0.0
                 do idim=1,nvar
		    dx=dx+invvar(idim)
     &	              *(Fspace(idim,Indices(j))-Current(idim))**2.
                 enddo

                 if (dx .lt. sa(ichoice))  then
                    i=ichoice
                    do while(sa(i-1).gt.dx)
                       sa(i) = sa(i-1)
                       si(i) = si(i-1)
                       i=i-1
                    enddo
                    sa(i) = dx
                    si(i) = Indices(j)
                 endif  !... if distance is shorter than largest in list
	   endif  ! valid successor 	     
        enddo  ! loop over indices j

        return
        end
!============================================================================


!OLD SUB USED (AND MODIFIED) FOR OLDER VERSIONS OF RESAMPLE USED (AND MODIFIED) FOR OLDER VERSIONS OF RESAMPLE USED (AND MODIFIED) FOR OLDER VERSIONS OF RESAMPLE USED (AND MODIFIED) FOR OLDER VERSIONS OF RESAMPLE
!============================================================================
        SUBROUTINE sortnnb_old(
     &                          si,             ! indices of sorted nearest neighbours (out)
     &                          Indices,        ! indices of days within the window
     &                          Fspace,         ! feature vector space
     &                          Current,        ! feature vector of the 'current' state
     &                          invvar,         ! weights of feature-vector elements
                                                ! (=1/variance)
     &                          ichoice,        ! randomly selected nearest neighbour
     &                          ww,             ! half-window width
     &                          nvar,           ! feature-vector dimension
     &                          y0, y1          ! first and last year (in the array fspace) to be searched
     &                  )
        implicit none
	include 'params.f'

        integer si(*),ichoice
        integer ww                      ! window width
        integer nvar                    ! feature-vector dimension
        integer Indices(*)              ! window with indices
                                        ! subdivided in rows with
                                        ! length 2*ww+1
        real Fspace(nvar,*)             ! historical vectors in feature-space

        real Current(nvar)
        real invvar(nvar)
        integer i, j, idim
        integer firstndx,lastndx
        real dx, sa(0:NKMAX)      
        integer y0,y1

        ! Each row in Indices had a length of -ww..0..ww = 2*ww+1
        ! this is the first (fastest) index
        firstndx = (y0-1)*(ww*2+1)+1    ! first considered entry of Indices
        lastndx = y1*(ww*2+1)           ! last considered entry of Indices

        do i=1,ichoice
           sa(i)=1E+14
           si(i)=-1
        enddo ! i
        sa(0) = -1              ! any distance is beyond sa(0)

        ! check the remaining days in the window for closer
C        do j=firstndx, lastndx          ! calculate the Euclidean distances and gather the
                                        ! ichoice smallest in sa
        do j=1,(y1-y0+1)*(2*ww+1)
                 dx=0.0
                 do idim=1,nvar
		    dx=dx+invvar(idim)
     &	              *(Fspace(idim,Indices(j))-Current(idim))**2.
                 enddo

                 if(Indices(j).eq.y1*365) dx=9999.9             ! exclude last historical day
                 if (dx .lt. sa(ichoice))  then
                    i=ichoice-1
                    do while(sa(i).gt.dx)
                       sa(i+1) = sa(i)
                       si(i+1) = si(i)
                       i=i-1
                    enddo
                    sa(i+1) = dx
                    si(i+1) = Indices(j)
                 endif  !... if distance is shorter than largest in list
        enddo  ! loop over indices


        return
        end
!============================================================================


	subroutine ssort(x,n)
	implicit none 
	integer n,i,j
	real x(n),dummy
	do i=2,n
	   dummy=x(i)
	   do j=i,2,-1
	      if(x(j-1).le.dummy) then 
	         goto 410 
	      else 
	         x(j)=x(j-1)
	         x(j-1)=dummy
	      endif 
	   enddo ! j 
 410	   continue 
	enddo 
	return 
	end 

	subroutine ssortxy(x,y,n)
	implicit none 
	! sort y along with x as key 
	integer n,i,j
	real x(n),y(n),dummy,dummy_y
	do i=2,n
	   dummy=x(i)
	   dummy_y=y(i)
	   do j=i,2,-1
	      if(x(j-1).le.dummy) then 
	         goto 430 
	      else 
	         x(j)=x(j-1)
	         x(j-1)=dummy
	         y(j)=y(j-1)
	         y(j-1)=dummy_y
	      endif 
	   enddo ! j 
 430	   continue 
	enddo 
	return 
	end 


	integer function doy_leap(iy) 
	implicit none 
	integer iy 
	doy_leap=366-ior(iand(iy,1),ishft(iand(iy,2),-1))
	return 
	end 


        integer function date2cal(datum,shift)
        implicit none 
        ! convert date to day-of-year 
        ! shift=1 : leave out 0731 in leapyears, iow shift between 0228 and 0801
        ! shift=0 : leave out 0229 in leapyears 
        ! shift=-1 : do not leave out days 
        integer yr,day,mnth,datum,shift 
        integer dayzero1(12)    ! day zero in each month 
        integer dayzero2(12)    ! day zero in each month 
        integer dayzero0(12)    ! day zero in each month 
        data dayzero1 /0,31,59,90,120,151,181,212,243,273,304,334/      !no leap
        data dayzero2 /0,31,60,91,121,152,182,212,243,273,304,334/      !leap AN
        data dayzero0 /0,31,60,91,121,152,182,213,244,274,305,335/      !full leap year 
        yr=datum/10000
        mnth=mod(datum,10000)/100
        day=mod(datum,100)
        if(mod(yr,4).eq.0 .and. shift.eq.1) then 
           date2cal=dayzero2(mnth)+day
        else if(mod(yr,4).eq.0 .and. shift.eq.-1) then 
           date2cal=dayzero0(mnth)+day
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
        if(mod(iy,4).eq.0 .and. 
     &     (mod(iy,100).ne.0 .or. mod(iy,1000).eq.0)) then
           do im=12,1,-1
              if(calday.ge.leap(im)) goto 254
           enddo
 254       continue
           id=calday-leap(im)+1
        else
           do im=12,1,-1
              if(calday.ge.no_leap(im)) goto 260
           enddo
 260       continue
           id=calday-no_leap(im)+1
        endif
        cal2date=iy*10000+im*100+id
        return
        end


        logical function inseason(datum1,datum2,datum)
        implicit none
C       returns .True. if datum is within the season, defined by
C       the first day datum1 and the last day datum2
        integer datum1,datum2,datum
        integer d1,d2,dd12,d
        d=mod(datum,10000)
        d1=mod(datum1,10000)
        d2=mod(datum2,10000)
        dd12=mod(mod(d2-d1,1232)+1232,1232)
        inseason=(mod(d-d1+1232,1232).le.dd12)
        return
        end

        logical function inseason_shift(datum1,datum2,datum,shft)
        implicit none
C       returns .True. if at [datum], the day [shft] days ago is within the season,
C       defined by the first day datum1 and the last day datum2
        integer date2cal
        integer datum1,datum2,datum,shft
        integer d1,d2,dd12,d
        d=mod(date2cal(datum,0)-shft+365,365)+1
        d1=date2cal(datum1,0)
        d2=date2cal(datum2,0)
        dd12=mod(mod(d2-d1,366)+366,366)
        inseason_shift=(mod(d-d1+366,366).le.dd12)
        return
        end


        subroutine dumpargsource
        implicit none
        character*(50) myname
        integer lastchar
        call getarg(0,myname)
        write(0,*)
        call system('cat '//myname(1:lastchar(myname))//'.f '
     &      //' | grep -i arg 1>&2')
        write(0,*)
        return
        end

        logical function GregorianLeap(iyr)
        implicit none
        integer iyr
           GregorianLeap = ((mod(iyr,4).eq.0) .and. 
     &          .not. ((mod(iyr,100).eq.0) .and.
     &            .not.(mod(iyr,400).eq.0)))
        return
        end

        integer function days_since(datum,refyr)
        implicit none
        ! Gives the number of days at the specified date  since the year 2000 started,
        ! according to the Gregorian calendar, so that refyr-01-01 yields 1
        integer datum, refyr
        logical GregorianLeap
        integer days, iyr
        integer iy,calday,im,id
        integer leap(13),no_leap(13), leap_shift(13)
        data leap /1,32,61,92,122,153,183,214,245,275,306,336,367/          ! standard leap year (366 days)
        data no_leap /1,32,60,91,121,152,182,213,244,274,305,335,366/       ! standard non-leap year (365 days)
        days = 0
        iy = int(datum/10000)
        im = int(mod(datum,10000)/100)
        id = mod(datum,100)
        do iyr = refyr, iy-1
           if (GregorianLeap(iyr)) then
              days = days + 366
           else
              days = days + 365
           endif
        enddo
        if (GregorianLeap(iy)) then
           days=days+leap(im)-1+id
           if (id.le.leap(im+1)-leap(im)) then 
               days_since=+days
           else
               days_since=-days
           endif
        else
           days=days+no_leap(im)-1+id
           if (id.le.no_leap(im+1)-no_leap(im)) then
               days_since=+days
           else
               days_since=-days
           endif
        endif
        end function days_since




        integer function date2ndx(datum,leap)
        implicit none
        ! returns the order of the day of datum from 19000101 (the first day)
        ! leap=1 take leap years into account 
        ! leap=0 do NOT count leap years: each year 365 days 
        integer dayzero1(12)    ! day zero in each month 
        integer iy,day,mnth,datum
        integer yy
        integer dayzero2(12)    ! day zero in each month 
        integer dayzero3(12)    ! day zero in each month 
        integer nleap
        integer leap
        data dayzero1 /0,31,59,90,120,151,181,212,243,273,304,334/      ! ordinary year 
        data dayzero2 /0,31,60,91,121,152,182,213,244,274,305,335/      ! leapyear 
        data dayzero3 /0,31,60,91,121,152,182,212,243,273,304,334/      ! shifted 
        yy=datum/10000
        iy=yy - 1900                           ! 1900 is not a leap year 
        nleap=(iy-1)/4                                  ! number of leap years passed since 1900
        mnth=mod(datum,10000)/100
        day=mod(datum,100)
        if(leap.gt.0) then
           if(mod(yy,4).eq.0 .and. 
     &        (mod(yy,100).ne.0 .or. mod(yy,400).eq.0)) then
              date2ndx=dayzero2(mnth)+day+365*(iy)+nleap        ! leap year 
           else
              date2ndx=dayzero1(mnth)+day+365*(iy)+nleap        ! ordinary year 
           endif
        else
           if(mod(yy,4).eq.0 .and.
     &        (mod(yy,100).ne.0 .or. mod(yy,400).eq.0)) then
              date2ndx=dayzero3(mnth)+day+365*(iy)      ! leap year 
           else
              date2ndx=dayzero1(mnth)+day+365*(iy)      ! ordinary year 
           endif
        endif

        return
        end



        integer function ndx2date(dayssince, startdate)
        implicit none
        ! (inverse of date2ndx with leap=1)
        ! returns the date (yyyymmdd) given de number of days since startdate (yyyymmdd) 
        ! USES function JD and subroutine GDATE (see below)
        ! and therefore assumes leapyears
        
        integer dayssince, startdate
        integer startyyyy, startmm, startdd
        integer julday
        integer yyyy, mm, dd
        integer jd

        startyyyy = startdate/10000
        startmm   = mod(startdate, 10000)/100
        startdd   = mod(startdate, 100)

        julday = dayssince + jd(startyyyy, startmm, startdd)
        call gdate(julday, yyyy, mm, dd)

        ndx2date = yyyy*10000 + mm*100 + dd
        
        return
        end

        
        INTEGER FUNCTION JD (YEAR,MONTH,DAY)
C
C-------COMPUTES THE JULIAN DATE (JD) GIVEN A GREGORIAN CALENDAR
C       DATE (YEAR,MONTH,DAY).
C
        INTEGER YEAR,MONTH,DAY,I,J,K
C
        I= YEAR
        J= MONTH
        K= DAY
C
        JD= K-32075+1461*(I+4800+(J-14)/12)/4+367*(J-2-(J-14)/12*12)/12
     &       -3*((I+4900+(J-14)/12)/100)/4
C
        RETURN
        END

        
        SUBROUTINE GDATE (JD, YEAR,MONTH,DAY)
C
C-------COMPUTES THE GREGORIAN CALENDAR DATE (YEAR,MONTH,DAY)
C       GIVEN THE JULIAN DATE (JD).
C
        INTEGER JD,YEAR,MONTH,DAY,I,J,K
C
        L= JD+68569
        N= 4*L/146097
        L= L-(146097*N+3)/4
        I= 4000*(L+1)/1461001
        L= L-1461*I/4+31
        J= 80*L/2447
        K= L-2447*J/80
        L= J/11
        J= J+2-12*L
        I= 100*(N-49)+I+L
C
        YEAR= I
        MONTH= J
        DAY= K
C
        RETURN
        END     




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

        logical function ibnd(x,x0,x1)
        ! returns true if x element of [x0,x1], false otherwise
        implicit none
        integer x,x0,x1
        if(x.gt.x1 .or. x.lt.x0) then
             ibnd=.False.
        else
             ibnd=.True.
        endif
        return
        end


!----------------------------------------------------------------------------
!--------------------------------------------------
! adopted from input.f 
!--------------------------------------------------
      subroutine error(myline)
      implicit none
      character     myline*(*)
      write(*,*)
      write(*,*) '***** Error: ',myline
      call abort
      end 
!--------------------------------------------------

