	program thwt_EOBS
	!use portlib 
	! purpose : derive areal precipitation from point precipitation by 
	! weighting over gridpoints 
	! each gridpoint inside the subcatchments is assigned the nearest RCM
        ! gridpoint
	! available for the considered day. The program adapts weigths by inspecting how 
	! the number of stations changes from one day onto another. 

        ! adapted for EOBS-precipitation (rotated grid .22 degrees resolution )


	! note 290403 : 
	!   currently solving the following problem. 
	!   number of polygons seems incorrect or there is a polygon (nr. 9) to 
	!   which none of the gridpoints are assigned, resulting in a sum of weights
	!   equal to zero  

	implicit none 
	include 'netcdf.inc'

	integer maxstn,maxpt,maxyr,maxarea,maxsave
	character*(10) fieldname
        character*(8) ext                               ! extension of output files, default 'out' 
C	parameter (fieldname='precipitation_flux')
	
	character*(20) fnout 

	parameter (maxstn=300*300)		! maximum number of stations
	parameter (maxyr=100)			! maximum number of years 
	parameter (maxpt=500000)		! maximum number of gridpoints 
	parameter (maxarea=200)			! maximum number of sub-areas 
	parameter (maxsave=100)			! maximum number of saved 'states' 
	real scale_factor
	real lowerbound
        integer FillValue
C	parameter (factor=24*3600.0)		! factor to multiply data with 
	parameter (lowerbound=-1e+7)
	
	! NETCDF RELATED MATERIAL
        integer  ierror, ncid
        integer instr

	logical argstring
        character*(255)   ncfile
        character*(50)    name
	character*(255)   lat_name
	
	integer latdim   ! variable shape
        integer ndim, nvar, ngatt, unlim_id
        integer time_dimid,lat_dimid,lon_dimid
        integer time_varid,lat_varid,lon_varid
	integer lat_type
	integer field_varid
        integer nlon,nlat,ntime
	real*8  tot_loc    
	integer start(3),cnt(3)
        real*4 lats(maxstn)
        real*4 lons(maxstn)
	real*4 latval,lonval
        real*4 wts(maxstn)
        real*4 values(maxstn)
	real*4 series(200*365)
	logical used(maxstn)

	real fx,fy 				! fx = km per degrees LON
	common /scale/fx,fy			! fy = km per degrees LAT 

	integer npt,npoly
	integer nstn(maxarea)		! number of stations for each area
	integer stnnr(600,maxarea)	! station numbers for each area 
	real stnwt(600,maxarea)		! station weigths for each area 
	integer stndata(600*365,600)	! data for each station 

        integer since                   ! date corresponding to the zeroth time level 
                                        ! usually indicated in the attributes with the 
                                        ! word 'since'         
        parameter(since=19510101)
	real times(366)			! time levels  		
	real*8 average(200*366)		! weighted average over the entire area 
	integer npt_tot(200*366)	! total number of interpolation gridpunts over the entire area 
        integer*8 sumnpt                ! totaal number of fine-grid points over the entire area 
	integer dates(0:200*366)	! dates for time levels 

	integer istn,cstn,ipt,ipoly,itime
	real gridpt(2,maxpt)	
	real*8 value 
	real areal(0:maxarea)			! area values for all sub-catchments 
						! for a single day  
						! 0 contains entire area, all catchments together
	character*(100) line_of_data		! line read from the file 


	! POLYGON INFORMATION INTO GRIDPOINTS
	integer maxx,maxy,maxintersect,maxnodes
	parameter(maxx=2000,maxy=2000)
	parameter(maxintersect=200)	! no more than 20 intersections on one line 
	parameter(maxnodes=10000)	! maximum number of nodes in pgon 
	
	real xg					! xlevel (which is LON)
	real yg					! ylevel (which is LAT)
	integer polynr
	real xmin, ymin, xmax, ymax	        ! xmin, ymin, xmax, ymax 
	real x_int(maxintersect)		! intersectionpoints 
	real xx,x1,x2,y1,y2					! temporary intersection point 
	integer n_int				! number of intersections at an y-level
	character*(4) str_year			! year information, included in the netcdf filename 
	real delta_lon,delta_lat
	real delta_x,delta_y			! grid spacing in km 
	real dist,distance,mindist
	parameter (delta_x=2.5,delta_y=2.5)     ! according to belgian data 
	real r_e,r_p,pi,theta0
	parameter(r_e = 6377.8)			! equatorial earth readius 
	parameter(r_p = 6357)			! polar earth radius
	parameter(pi = 3.1415926)	
	integer i,ix,iy,nx,ny 
	real stnwt_sum				! keeps track of the sum of station weights for each sub basin
	real node(2,maxnodes)			! nodes of a pgon 
	integer nnode				! number of node in pgon 	
	integer inode 
	
        integer firstchar, lastchar 
        character*(20) dbuser, dbhost, dbname, dbpass 
        character*(300) datqry 

C	OPEN NETCDF DATA 
	if(.not.argstring('-nc','',ncfile)) then 
	   write(0,*) 'No NETCDF-file provided'
	   stop '----'
	endif
	
	str_year = ncfile(7:10)
	write(0,*) str_year
 
	if(.not.argstring('-var','pr',fieldname)) then 
	   write(0,*) 'No varname provided'
	endif 
	if(.not.argstring('-ext',
     &         fieldname(firstchar(fieldname):lastchar(fieldname)),
     &         ext)) then 
	   write(0,*) 'No extension provided'
	endif 

       ierror = nf_open(ncfile(1:instr(ncfile,' ')-1)
     &        , nf_nowrite, ncid)

       ierror = nf_inq(ncid, ndim, nvar, ngatt, unlim_id)	

       write (0,*) ndim,' dimensions, ', nvar,' variables', ngatt
       write (0,*)
       
       ierror = nf_inq_dimid (ncid, 'time', time_dimid)         ! time ID
       ierror = nf_inq_dim(ncid, time_dimid, name, ntime)       ! inquire t leve
       write (0,*) 'time dimID = ', time_dimid,'  ',ntime,' levels'

       ierror = nf_inq_dimid (ncid, 'y', lat_dimid)          ! time ID
       ierror = nf_inq_dim(ncid, lat_dimid, name, nlat)         ! inquire t leve
       write (0,*) 'Lat dimID  = ', lat_dimid,'  ',nlat,' levels', ncid
       
       ierror = nf_inq_dimid (ncid, 'x', lon_dimid)          ! time ID
       ierror = nf_inq_dim(ncid, lon_dimid, name, nlon)         ! inquire t leve
       write (0,*) 'Lon dimID  = ', lon_dimid,'  ',nlon,' levels'

       !if (nlat .le. 0 .and. nlon .le. 0) then
       !	   nlat = 52800.	
       !	   nlon = 52800.
       !endif		

       write(0,*)
       ierror = nf_inq_varid (ncid, 'time', time_varid)
       write (0,*) 'Time  varID = ', time_varid
       ierror = nf_inq_varid (ncid, 'lat', lat_varid)
       write (0,*) 'Lat  varID = ', lat_varid
       ierror = nf_inq_varid (ncid, 'lon', lon_varid)
       write (0,*) 'Lon  varID = ', lon_varid
	
       !check the name of the latitude variable - 29/02/2012 EW	
       ierror = nf_inq_var(ncid, lat_varid, name, lat_type, ndim, 
     & latdim, ngatt)
       write(0,*) 'latitude info: ', ierror, ngatt

       ierror = nf_inq_varid (ncid, 
     &          fieldname(firstchar(fieldname):lastchar(fieldname)), 
     &          field_varid)	                                        ! variable of interest
       write (0,*) fieldname(firstchar(fieldname):lastchar(fieldname))
     &             //' varID = ', field_varid

       scale_factor=-1.0 
       ierror=nf_get_att_real(ncid,field_varid,'scale_factor',
     &      scale_factor)
       if(scale_factor.gt.0) then 
         write(0,*) 'Found scale factor : ', scale_factor 
       else 
         scale_factor=0.1d0
         write(0,*) 'No scale factor  found, defaulting to ', 
     &      scale_factor 
       endif 
        
       ierror=nf_get_att_int(ncid,field_varid,'_FillValue',
     &      FillValue)
        write(0,*) 'Fillvalue = ',FillValue 

       ierror = nf_get_var_real(ncid, lat_varid, lats)
       ierror = nf_get_var_real(ncid, lon_varid, lons)
       ierror = nf_get_var_real(ncid, time_varid, times)
    
       !write(0,*), times

!-------SAMPLE SHOWING HOW TO EXTRACT DATA FOR A STATION FROM THE NETCDF-FILE
C	istn=300
C	start(1)=mod(istn-1,nlon)+1
C	start(2)=(istn-1)/nlon+1
C	start(3)=1
C	cnt(1)=1
C	cnt(2)=1
C	cnt(3)=1
C           ierror = nf_get_vara_real(ncid, lat_varid,start,cnt,latval)
C           ierror = nf_get_vara_real(ncid, lon_varid,start,cnt,lonval)
C	write(*,*) latval,lonval
C	write(*,*) lats(istn),lons(istn)
!-------SAMPLE SHOWING HOW TO EXTRACT DATA FOR A STATION FROM THE NETCDF-FILE
	
C	********************************************************************
	write(0,*) 'nlat*nlon= ', nlat*nlon
	do istn=1,nlat*nlon		
	   used(istn)=.False.
	enddo 
C	********************************************************************

	write(0,'(a,f6.2,a2,f6.2,a)') ' grid resolution : '
     &	         ,delta_x,' x',delta_y,' km' 

	! READ POLYGONS FROM STANDARD INPUT 
	ipoly=0
	open(34,file='dumpgrid.out')
	open(35,file='weights.out')
        open(55,file='area_npt.out')

	do itime=1,ntime 
          average(itime)=0.d0
          npt_tot(itime)=0
        enddo 
        sumnpt=0

	do while (.True.)
	   read(*,'(a100)',end=666) line_of_data
	   !write(0,*) line_of_data
	   read(line_of_data,*) polynr,xmin,ymin,xmax,ymax,nnode
	   ipoly=ipoly+1
	   write(0,'(a,i4.4)') 'Reading definitions of area ',ipoly
	   theta0=(ymax+ymin)/2.		! average latitude of this box 
           fx = r_e*cos(theta0*pi/180)*pi/180   ! factors between degree
           fy = r_p*pi/180
           delta_lon=delta_x/fx
           delta_lat=delta_y/fy
	   !write(0,*) delta_lon, delta_lat

	   do inode=1,nnode
              read(*,'(a100)') line_of_data
	      read(line_of_data,*) node(1,inode),node(2,inode)
	   enddo 
	   ipt=0
	   ny=int((ymax-ymin)/delta_lat)+2
	   
	   do iy=0,ny-1
	      yg=ymin-0.5*delta_lat+iy*delta_lat
	      n_int=0

	      ! start scanning pgon for intersections 
	      do inode=1,nnode-1
	         x1=node(1,inode)
	         x2=node(1,inode+1)
	         y1=node(2,inode)
	         y2=node(2,inode+1)
	         if(abs(y2-yg).gt.1e-014) then 
	            xx=(y2-yg)*(y1-yg)          
	            if(xx.le.1e-014) then 		! level intersected 
		       
	               n_int=n_int+1	
	               x_int(n_int)=
     &	                 (x2*yg+x1*y1-x2*y1-x1*yg)/(y2-y1)+x1
	            endif 
	         endif 
	      enddo ! inode 

	      if(n_int.gt.0) then 
	         call ssort(x_int,n_int)		! sort first n_int elements ascendingly 
	         nx=int((xmax-xmin)/delta_lon)+2
	         do ix=0,nx-1				! run over gridpoints 
	            xg=xmin-0.5*delta_lon+ix*delta_lon
	            do i=1,n_int
	               if(x_int(i).gt.xg) exit 
	            enddo ! intersections 
	            i=i-1				! i holds the number of intersections passed 
	            if(mod(i,2).ne.0) then 
	                ipt=ipt+1 
	                gridpt(1,ipt)=xg		! store x-coordinate 
	                gridpt(2,ipt)=yg		! store y-coordinate 
		    endif 
	         enddo ! ix 				
	      endif
	   enddo ! iy 
	   npt=ipt			! number of gridpoints in this area 

	   ! CALCULATE WEIGHTS FOR EACH OF THE GRIDBOXES IN THIS SUBBASIN 

	   do istn=1,nlon*nlat
	      wts(istn)=0.0
	   enddo 
	   do ipt=1,npt 
	      mindist=1e+07 
	      do istn=1,nlon*nlat
C		 dist=distance(gridpt(1,ipt),gridpt(2,ipt),
C     &		         lons(istn),lats(istn),r_e) 
	         dist=(fx*(gridpt(1,ipt)-lons(istn)))**2. 
     &	             +(fy*(gridpt(2,ipt)-lats(istn)))**2. 
	         if(dist.lt.mindist) then 
		    mindist=dist
		    cstn=istn
		 endif 
	      enddo  ! istn
	      used(cstn)=.True.
	      wts(cstn)=wts(cstn)+1.0
	      
	      write(34,*)
	      write(34,'(2e12.4)') lons(cstn),lats(cstn)
	      write(34,'(2e12.4)') gridpt(1,ipt),gridpt(2,ipt)
	   enddo 
	   !write(0,*) 'wts ' , wts(1)
	   	
	   ! add stations to the list of stations for this area 
	   nstn(ipoly)=0
	   !write(0,*) 'nlon, nlat ', nlon, nlat
	   do istn=1,nlon*nlat
	      
	      if(wts(istn).ge.1.) then
		 nstn(ipoly)=nstn(ipoly)+1
		 stnnr(nstn(ipoly),ipoly)=istn		! number 
		 !write(0,*) nstn(ipoly)
	         stnwt(nstn(ipoly),ipoly)=wts(istn)/npt ! weight
	      endif 
	   enddo 

	   ! WRITE WEIGHTS TO FILE 
	   write(35,'(a,i4.4)') 'area ',ipoly
	   do istn=1,nstn(ipoly)
	      write(35,'(2i5,f10.7)'), istn,stnnr(istn,ipoly),
     &	                               stnwt(istn,ipoly)
	   enddo 
	   write(35,*)

	   ! WRITE SERIES DATA 
C	   write(fnout,'(a,i3.3,a)') 'area',ipoly,'.rr'
	   write(fnout,'(a,i3.3,a,a)') 'area_h',ipoly,'_',str_year,'.out'
	   write(0,*) fnout
	   open(23,file=fnout(1:instr(fnout,' ')-1))
	   open(24,file='fill_values.dat')
	   
	   do istn=1,nstn(ipoly)
	      start(1)=mod(stnnr(istn,ipoly)-1,nlon)+1
	      start(2)=(stnnr(istn,ipoly)-1)/nlon+1
	      start(3)=1
	      cnt(1)=1
	      cnt(2)=1
      	      cnt(3)=ntime
              ierror = nf_get_vara_int(ncid, field_varid,start,cnt
     &           ,stndata(1,istn))
	      !write(0,*) 'ierror ' , ierror, field_varid
              IF (ierror.ne.NF_NOERR) call handle_error(ierror)        !
	   enddo 
	   
	   do itime=1,ntime 
	      value=0.d0
	      do istn=1,nstn(ipoly)
                 if(stndata(itime,istn).le.FillValue+1.) then           ! Assume FillValue to be NEGATIVE !!!!   
		     goto 384
                 endif
		 stnwt_sum = stnwt_sum + stnwt(istn,ipoly)
	         value=value+stnwt(istn,ipoly)*
     &	         max(stndata(itime,istn)*scale_factor,lowerbound) 
		 
		 if (itime .eq. 1) then 
		     write(24,'(f10.8)') stnwt_sum	  
		 endif
 384		 continue	
	      enddo
	      
	      !EW - 07/03/2012, check if the sum of the weights equals 1 and normalize if not
	      if (stnwt_sum .lt. 1. .and. stnwt_sum .ge. 0.5) then
	          value = value / stnwt_sum	
	      endif				
	      	
	      stnwt_sum = 0.		      	
	      average(itime)=average(itime)+value*npt                   ! average over all sub-areas
              npt_tot(itime)=npt_tot(itime)+npt                         ! count interpolation gridpoints 
 !384          continue 
	      write(23,'(i10,f10.3)') dates(times(itime)),value 
	   enddo 
	   	
	   close(23)
           write(55,'(a25,i5)') fnout(1:lastchar(fnout)), npt 
           sumnpt=sumnpt+npt
	enddo ! reading from polygon definitions 	

 666    continue
	close(34)
	close(35)
        close(55)
	close(24)

	ierror=nf_close(ncid)
	write(0,*) 'NetCDF file closed' 

C tot hier 
	write(0,*) nlon*nlat
	open(66,file='used_stns.out')
	do istn=1,nlon*nlat
	   if(used(istn)) then 
	     write(66,'(2e15.5,i6)') lons(istn),lats(istn),istn
	   endif 
	enddo 
	close(66)

        open(44,file='area000.'//ext(firstchar(ext):lastchar(ext)))
        do itime=1,ntime 
          if(npt_tot(itime)*1.0/sumnpt.ge.0.8) then 
	    write(44,'(i10,f10.3)') dates(times(itime)),
     &                             average(itime)/npt_tot(itime)
          else 
	    write(44,'(i10,f10.3)') dates(times(itime)), FillValue*1.0
        endif 
        enddo 
        close(44)

	end 



	real function norm(x1,x2)
	! returns squared norm of x2-x1 scaled with fx and fy 
	implicit none 
	real x1(2),x2(2) 
	real fx,fy 
	common /scale/fx,fy
	norm=(fx*(x1(1)-x2(1)))**2.+(fy*(x1(2)-x2(2)))**2.
	return 
	end 


	subroutine add(a,b)
	implicit none 
	real a,b
	a=a+b
	return 
	end 	

	subroutine mul(a,b)
	implicit none 
	real a,b
	a=a*b
	return 
	end 	

	subroutine inc(a)
	implicit none 
	real a
	a=a+1.0
	return 
	end 

	subroutine dec(a)
	implicit none 
	real a
	a=a-1.0
	return 
	end 



	subroutine bool2hex(P,s,n)
	! in : integer array, consisting of n zeroes/ones 
	! out: string s representing each 4 elements of P as a hex-digit 
	integer n,P(n)
	character(*) s
	integer nibble,ibit,bit_nr 
	character*(*) hexdigit
	
	parameter (hexdigit='0123456789ABCDEF')
	nibble=0 
	s=''
	do ibit=1,n 
	   bit_nr=mod(ibit-1,4)
	   nibble=nibble+P(ibit)*(2**(4-bit_nr-1))
	   if(bit_nr.eq.3) then 
	      s(ibit/4:ibit/4)=hexdigit(nibble+1:nibble+1)
     	      nibble=0
	   endif 
	enddo ! ibit 
	s(n/4+1:n/4+1)=hexdigit(nibble+1:nibble+1)
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

	subroutine myknl(knl,span)
	implicit none 
	! fill knl(-span:span) with Epanechnikov smoothing kernel 
	integer span, i 
	real knl(-span:span), som, alpha 
	
        som = 1.0
        knl(0) = 1.0 
        do i=1, span
	    alpha=(i*1.0/span)
            knl(i) = (1-alpha**2)**2		! fourth order function 
            knl(-i) = knl(i)
            som = som + 2*knl(i)
        enddo
        do i=-span, span
            knl(i) = knl(i)/som
        enddo
	return 
	end 
	
	subroutine perconv(x,y,res,nx,ny)
	implicit none 
	! performs a periodical convolution on x(1:n) 
	! with array y(-ny:ny) 
	! intended for smoothing x 
	! output in res(1:nx) 
	integer nx, ny, l, i, k  
	real x(1:nx), y(-ny:ny), res(1:nx)

        do l=1, nx 
	   res(l)=0.0
           do i=-ny, ny
              k = mod(l+i-1+nx,nx) + 1			! mod can't handle negative
	      res(l) = res(l)+x(k)*y(i)
           enddo ! i 
	enddo ! j 
	return
	end 

	integer function date2ndx(datum,leap)
	implicit none 
	! returns the order of the day of datum from 19000101 (the first day)
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
	iy=datum/10000 - 1900				! 1900 is not a leap year 
	nleap=(iy-1)/4					! number of leap years passed since 1900
	mnth=mod(datum,10000)/100
	day=mod(datum,100)
	if(leap.gt.0) then 
	   if(mod(iy,4).eq.0 .and. mod(iy,100).ne.0) then
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
        
	real function distance(lat1,lon1,lat2,lon2,r)
	implicit none 
	real lat1,lon1,lat2,lon2,r 
	real phi1,phi2,the1,the2
	real pi 
	pi=4.*atan(1.0)
	phi1=lon1*pi/180.
	phi2=lon2*pi/180.
	the1=lat1*pi/180.
	the2=lat2*pi/180.
	distance = 
     &	  r* acos(cos(the1)*cos(the2)*cos(phi2-phi1)
     &	              + sin(the1)*sin(the2))
	return 
	end 

	
        integer function ceil(a)
        implicit none 
        real a 
        if(a.gt.0) then 
                ceil=a+1 
        else 
                ceil=a
        endif 
        return 
        end 

        integer function flor(a)
        implicit none 
        real a 
        if(a.lt.0) then 
                flor=a-1 
        else 
                flor=a
        endif 
        return 
        end 

	
       Subroutine handle_error (ierror)
         include 'netcdf.inc'
         integer ierror
         if (ierror.ne.nf_noerr ) then
           print *, nf_strerror(ierror)
           stop "stopped in handle_error"
         end if
       end


        
       	include 'routines.f' 
