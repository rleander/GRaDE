       program thwt_EOBS_Rhine
       !use portlib 
       ! purpose : derive areal precipitation from point precipitation by 
       ! weighting over gridpoints 
       ! each gridpoint inside the subcatchments is assigned the nearest RCM
       ! gridpoint
       ! available for the considered day. The program adapts weigths by inspecting how 
       ! the number of stations changes from one day onto another. 

       ! adapted for EOBS-precipitation (rotated grid .22 degrees resolution )

       ! Compiled with:
       ! gfortran -c -I /usr/include thwt-EOBS_v15.0_Rhine.f -Wall 
       ! gfortran -o thwt_EOBS_v15 thwt-EOBS_v15.0_Rhine.o -L/usr/lib/x86_64-linux-gnu -lnetcdff

       ! Requires:
       ! routines,f
       ! params.f

       
       ! converts EOBS precip. and temp to
       ! HBV_Rhine subbasins (134) using rhinebasins134.txt and arc.txt

       !note 290403 :
       !JB Is this 20030429 ???
       !   currently solving the following problem. 
       !   number of polygons seems incorrect or there is a polygon (nr. 9) to 
       !   which none of the gridpoints are assigned, resulting in a sum of weights
       !   equal to zero

       ! 20170613 Jules Beersma (JB)
       !          adapted to:
       !          * new lat and long definition in EOBS_v15
       !          * include the historical date to each precip./temp. value (previously done with separate IDL programme)
       !          * renumber subbasins using ARCGIS- to HBV-numbering in arc.txt 
       !
       !

       use point_in_poly
       implicit none 
       include 'netcdf.inc'     
!JB routines.f (van R. Leander) is also included, see last line of this program
!JB params.f   (van R. Leander) is also included via routines.f and should therefore also be in the main directory
!JB compile with: f77 thwt-EOBS_Rhine.f -lnetcdff
       
       integer maxstn,maxpt,maxyr,maxarea,maxsave
       character*(20) lonname, latname, timname
       character*(10) fieldname
       character*(50) time_units
       character*(8)  ext                               ! extension of output files, default 'out' 
       logical udunit2date
C      parameter (fieldname='precipitation_flux')
       
       character*(20) fnout 

       parameter (maxstn=300*300)       ! maximum number of stations
       parameter (maxyr=100)            ! maximum number of years 
       parameter (maxpt=500000)         ! maximum number of gridpoints 
       parameter (maxarea=200)          ! maximum number of sub-areas 
       parameter (maxsave=100)          ! maximum number of saved 'states' 
       real scale_factor,lowerbound
       integer FillValue
C      parameter (factor=24*3600.0)     ! factor to multiply data with 
       parameter (lowerbound=-1e+7)
       logical use_fine_grid
       parameter (use_fine_grid = .False.) ! True:  intermediate fine grid used in interpolation
                                           ! False: siple averaging over source points 

       
       ! NETCDF RELATED MATERIAL
       integer  ierror, ncid
       integer instr
       character*(12) fdump
       logical parse_udunit
       integer tunit

       logical argstring, arglogical
       character*(255)    ncfile, ncbatch
       character*(50)    name

       integer ndim, nvar, ngatt, unlim_id
       integer idim
       integer time_dimid,lat_dimid,lon_dimid
       integer time_varid,lat_varid,lon_varid
       integer field_varid
       integer field_dimids(3)
       integer field_dimlen(3)
       integer ngd1,ngd2,ntime                 ! size n of Grid Dimension 1 and 2 and time  
       integer start(3),cnt(3)
       real*4 lats(maxstn)
       real*4 lons(maxstn)
       real*4 latval,lonval
       real*4 wts(maxstn)
       real*4 values(maxstn)
       real*4 series(200*365)
       logical used(maxstn)
       character*(8) dumpfilename, polyfilename

       real fx,fy ! fx = km per degrees LON
       common /scale/fx,fy! fy = km per degrees LAT 

       integer npt,npoly
       integer nstn(maxarea)            ! number of stations for each area
       integer basin(maxarea)           ! JB: conversion of ARCGIS numbering to HBV_Rhine numbering 
       integer stnnr(100,maxarea)       ! station numbers for each area 
       real stnwt(100,maxarea)          ! station weigths for each area 
       integer stndata(200*365,100)     ! data for each station 

       integer torigin, torigin_default ! date corresponding to the zeroth time level 
                                        ! usually indicated in the attributes with the 
                                        ! word 'since'         
       parameter(torigin_default=19500101)
       integer times(200*365)           ! time levels  
       real*8 average(200*365)          ! weighted average over the entire area 
       integer npt_tot(200*365)         ! total number of interpolation gridpunts over the entire area 
       integer*8 sumnpt                 ! totaal number of fine-grid points over the entire area 
       integer dates(0:200*365)         ! dates for time levels 
                                        ! JB starting with zero since the first time level (days since) 
                                        ! in the NetCDF file is zero

       integer istn,cstn,ipt,ipoly,itime
       real gridpt(2,maxpt)
       real*8 value 
       real areal(0:maxarea)            ! area values for all sub-catchments 
                                        ! for a single day  
                                        ! 0 contains entire area, all catchments together
       character*(100) line_of_data     ! line read from the file 

       ! POLYGON INFORMATION INTO GRIDPOINTS
       integer maxx,maxy,maxintersect,maxnodes
       parameter(maxx=2000,maxy=2000)
       parameter(maxintersect=200)      ! no more than 20 intersections on one line 
       parameter(maxnodes=10000)        ! maximum number of nodes in pgon 
       
       real xg                          ! xlevel (which is LON)
       real yg                          ! ylevel (which is LAT)
       integer partnr, polynr, oldpolynr
       real xmin, ymin, xmax, ymax      ! xmin, ymin, xmax, ymax 
       real x_int(maxintersect)         ! intersectionpoints 
       real xx,x1,x2,y1,y2              ! temporary intersection point 
       integer n_int                    ! number of intersections at an y-level 
       real delta_lon,delta_lat
       real delta_x,delta_y             ! grid spacing in km 
       real dist,distance,mindist
!      parameter (delta_x=2.5,delta_y=2.5)     ! according to belgian data 
       parameter (delta_x=0.6,delta_y=0.6)     ! according to belgian data 
       real r_e,r_p,pi,theta0
       parameter(r_e = 6377.8)          ! equatorial earth readius 
       parameter(r_p = 6357)            ! polar earth radius
       parameter(pi = 3.1415)
       integer i,ix,iy,nx,ny 
       logical series_append, batchmode, verbose

       real node(2,maxnodes)            ! nodes of a pgon 
!JB>
c      real node_x, node_y              !JB20160810 needed because of different structure HBV_Maas_polygons.txt
!JB<
       integer nnode                    ! number of node in pgon
       real stnwt_sum                   ! keeps track of the sum of station weights for each sub basin 
       integer inode
       
       integer firstchar, lastchar 
       character*(20) dbuser, dbhost, dbname, dbpass 
       character*(300) datqry 

!JB>
       integer index2date                      ! needed for using integer function ndx2date 
!JB<
       
C      OPEN NETCDF DATA 
       if(.not.argstring('-nc','',ncfile)) then                     !JB input NetCDF filename; GEEN default  
          if(.not.argstring('-ncbatch','',ncbatch)) then            !Batching nc-files
             stop 'No NETCDF-file or batch-file provided'
          else
             batchmode = .True.
          endif
       else
          batchmode = .False.
       endif 
       if(.not.argstring('-var','',fieldname)) then                 !JB variable name; 'rr' default 
          stop 'No field variable name (-var) !'
       endif 
       verbose = arglogical('-v')
       series_append = arglogical('-append')
       if(.not.argstring('-lat','',latname)) then
          stop 'No latitude variable name (-lat) !'
       endif 
       if(.not.argstring('-lon','',lonname)) then
          stop 'No longitude variable name (-lon) !'
       endif 
       if(.not.argstring('-tim','',timname)) then
          stop 'No time variable name (-tim) !'
       endif 
       if(.not.argstring('-ext', trim(fieldname),                   !JB extension name; 'variable name' default
     &         ext)) then 
          if (verbose) write(0,*) 'No extension provided'
       endif 

       call handle_error ( nf_open(trim(ncfile),nf_nowrite, ncid),
     &          '','Opening file '//trim(ncfile))
       write(0,'(a)') 'nc-file = '//trim(ncfile)

       call handle_error ( nf_inq_varid (ncid, trim(fieldname), 
     &          field_varid), '',trim(fieldname))                   ! variable of interest
       ierror = nf_inq_vardimid (ncid, field_varid, field_dimids) 

       if (verbose) write (0,*) trim(fieldname) 
     &                        //' varID = ', field_varid
       if (verbose) write (0,*) 'Dimensions:'
       do idim=1,3
          ierror = nf_inq_dim(ncid, field_dimids(idim), name, 
     &                        field_dimlen(idim))
          if (verbose) write(0,'(a,i4,a,i0,a)') '    ',
     &                        field_dimids(idim),
     &              ': '//trim(name)//' [',field_dimlen(idim),']' 
       enddo 
       ngd1 = field_dimlen(1)
       ngd2 = field_dimlen(2)
       ntime = field_dimlen(3)

       ! RL: 
       ! The variables holding the coordinates and mapping should be deductable from the header (attributes)
       ! But for EOBS, this is already not the case, not matching ANY convention
       ! Therefore, force the user to specify the names on the commandline, no defaults.
  
       call handle_error ( nf_inq_varid (ncid, trim(timname), 
     &          time_varid), '',trim(timname))
       if (verbose) write (0,*) 'Time  varID = ', time_varid
       call handle_error ( nf_inq_varid (ncid, trim(latname), 
     &          lat_varid), '', trim(latname))
       if (verbose) write (0,*) 'Lat  varID = ', lat_varid
       call handle_error ( nf_inq_varid (ncid, trim(lonname),
     &          lon_varid), '', trim(lonname))

!RL Extract time unit
       call handle_error ( 
     &      nf_get_att_text (ncid,time_varid,'units',time_units),
     &      '', '')
       if (.not.udunit2date (time_units,torigin)) then
          torigin = torigin_default          
       endif

!JB Extract scale factor from NetCDF file
! should be 0.1 for EOBS precip. and 0.01 for EOBS temp.
       scale_factor=-1.0
       ierror=nf_get_att_real(ncid,field_varid,'scale_factor',
     &      scale_factor)
       if (scale_factor.gt.0) then 
         if (verbose) write(0,*) 'Found scalefactor : ', scale_factor 
       else 
         scale_factor=0.1d0
         if (verbose) 
     &      write(0,*) 'No scalefactor  found, defaulting to ', 
     &      scale_factor 
       endif 
        
       ierror=nf_get_att_int(ncid,field_varid,'_FillValue',
     &      FillValue)
        if (verbose) write(0,*) 'Fillvalue = ',FillValue 

       call handle_error(nf_get_var_real(ncid, lat_varid, lats),'','')
       call handle_error(nf_get_var_real(ncid, lon_varid, lons),'','')
       call handle_error(nf_get_var_int(ncid, time_varid, times),'','')

!-------SAMPLE SHOWING HOW TO EXTRACT DATA FOR A STATION FROM THE NETCDF-FILE
C       istn=300
C       start(1)=mod(istn-1,ngd1)+1
C       start(2)=(istn-1)/ngd1+1
C       start(3)=1
C       cnt(1)=1
C       cnt(2)=1
C       cnt(3)=1
C           ierror = nf_get_vara_real(ncid, lat_varid,start,cnt,latval)
C           ierror = nf_get_vara_real(ncid, lon_varid,start,cnt,lonval)
C       write(*,*) latval,lonval
C       write(*,*) lats(istn),lons(istn)
!-------SAMPLE SHOWING HOW TO EXTRACT DATA FOR A STATION FROM THE NETCDF-FILE
       
C       ********************************************************************
       do istn=1,ngd2*ngd1
          used(istn)=.False.
       enddo 
C       ********************************************************************

       if (verbose) 
     &     write(0,'(a,f6.2,a2,f6.2,a)') ' grid resolution : '
     &         ,delta_x,' x',delta_y,' km' 

!JB>   
       ! READ POLYGONS FROM STANDARD INPUT 
       ipoly=0
       open(34,file='dumpgrid.out')           !JB output file for thtw_EOBS_Rhine?
       open(35,file='weights.out')            !JB output file for thtw_EOBS_Rhine?
       open(55,file='area_npt.out')           !JB output file for thtw_EOBS_Rhine?

       do itime=1,ntime 
          average(itime)=0.d0
          npt_tot(itime)=0
       enddo 
       sumnpt=0

!JB>
!-------Construct time series with dates from NetCDF time levels (days since)
       do itime = 1, ntime 
          dates(times(itime)) = index2date(times(itime), torigin)
       enddo 
!JB<

       polynr = -999
       if (verbose) write(0,*) 'EXTRACTING POLYGONS'

       do while (.True.)
          read(*,'(a100)',end=666) line_of_data
          oldpolynr = polynr
          read(line_of_data,*) polynr,partnr,xmin,ymin,xmax,ymax,nnode  
          !RL: Can do multiple parts of one shape, but not NEGATIVE parts 
          if (polynr.ne.oldpolynr) then
             if (oldpolynr.ge.0) then
                npt=ipt ! number of gridpoints in this area 
      
                ! CALCULATE WEIGHTS FOR EACH OF THE GRIDBOXES IN THIS SUBBASIN 
                if (use_fine_grid) then
                   do istn=1,ngd1*ngd2
                      wts(istn)=0.0
                   enddo 
                   do ipt=1,npt 
                      mindist=1e+07 
                      do istn=1,ngd1*ngd2
                         dist=(fx*(gridpt(1,ipt)-lons(istn)))**2. 
     &                       +(fy*(gridpt(2,ipt)-lats(istn)))**2. 
                         if (dist.lt.mindist) then 
                             mindist=dist
                             cstn=istn
                         endif 
                      enddo
                      used(cstn)=.True.
                      wts(cstn)=wts(cstn)+1.0
                   enddo 
         
                   ! add stations to the list of stations for this area 
                   nstn(ipoly)=0
                   do istn=1,ngd1*ngd2
                      if(wts(istn).ge.1.) then
                         nstn(ipoly)=nstn(ipoly)+1
                         stnnr(nstn(ipoly),ipoly)=istn          ! number 
                         stnwt(nstn(ipoly),ipoly)=wts(istn)/npt ! weight
                      endif 
                   enddo 
                else
                   if (pinpok2D(node(1,1:nnode)*1.d0, 
     &                          node(2,1:nnode)*1.d0,
     &                          lons(istn)*1.d0,lats(istn)*1.d0)) then
                      used(cstn)=.True.
                      nstn(ipoly)=nstn(ipoly)+1
                      stnnr(nstn(ipoly),ipoly)=istn          ! number 
                      stnwt(nstn(ipoly),ipoly)=1.0           ! weight
                   else
                      wts(istn)=0.0
                   endif
                endif
                      
                ! WRITE WEIGHTS TO FILE 
                write(35,'(a,i4.4)') 'area ',ipoly
                do istn=1,nstn(ipoly)
                   write(35,'(2i6,f10.7)') istn,stnnr(istn,ipoly),
     &                                  stnwt(istn,ipoly)
                enddo 
                write(35,*)
      
                ! WRITE SERIES DATA 
                write(fnout,'(a,i3.3,a)') 'area', ipoly, '.'
     &                            //ext(firstchar(ext):lastchar(ext))
                if (series_append) then
                   open(23,file=fnout(1:instr(fnout,' ')-1),
     &                access='APPEND')
                else
                   open(23,file=fnout(1:instr(fnout,' ')-1))
                endif
                 !write(0,*) 'nstn ' , nstn(ipoly)
                do istn=1,nstn(ipoly)
                   start(1)=mod(stnnr(istn,ipoly)-1,ngd1)+1
                   start(2)=(stnnr(istn,ipoly)-1)/ngd1+1
                   start(3)=1
                   cnt(1)=1
                   cnt(2)=1
                   cnt(3)=ntime
                   ierror = nf_get_vara_int(ncid, field_varid,start,cnt
     &                    ,stndata(1,istn))
                   call handle_error(ierror,'','')
                enddo 
                
                do itime=1,ntime 
                   value=0.d0
                   do istn=1,nstn(ipoly)
                      if (stndata(itime,istn).le.FillValue+1.) then    ! Assume FillValue to be NEGATIVE !!!!  
                         goto 384
                      endif
                      stnwt_sum = stnwt_sum + stnwt(istn,ipoly) 
                      value=value+stnwt(istn,ipoly)*
     &                max(stndata(itime,istn)*scale_factor,lowerbound)
                   enddo 
      
                   !EW - 08/03/2012, check if the sum of the weights equals 1 and normalize if not
                   if (stnwt_sum .lt. 1. .and. stnwt_sum .ge. 0.5) then
                       value = value / stnwt_sum
                   endif 
                   stnwt_sum = 0.
      
                   average(itime)=average(itime)+value*npt                   ! average over all sub-areas
                   npt_tot(itime)=npt_tot(itime)+npt                         ! count interpolation gridpoints 
 384               continue 
                   write(23,'(i10,f10.3)') dates(times(itime)),value 
                enddo 
                close(23)
                write(55,'(a25,i5)') fnout(1:lastchar(fnout)), npt 
                sumnpt=sumnpt+npt
             endif

             if (verbose) then 
                 write(0,'(a,i4.4,a,i4,a$)') 
     &                     'Reading definitions of area ',
     &                       polynr,', part ', partnr + 1, '' 
             else
                 write(0,'(i3.3,a$)') polynr, '.'
!                write(0,'(a,i4.4,a)') 
!    &                     'Reading definitions of area ',
!    &                       polynr, ' ...'
             endif

             ipoly=ipoly+1
             ipt=0
             close(65)
             close(66)
             write(dumpfilename,'(a5,i3.3)') 'dump.',polynr 
             write(polyfilename,'(a5,i3.3)') 'poly.',polynr 
             open(65,file=polyfilename) 
             open(66,file=dumpfilename) 
          else
             if (verbose)
     &           write(0,'(a,a4,a,i4,a$)') 
     &                    '                            ',
     &                       '   ',', part ', partnr + 1, '' 
          endif
          if (verbose) 
     &       write(0,'(4f8.4,i5)') xmin, ymin, xmax, ymax, nnode
          theta0=(ymax+ymin)/2. ! average latitude of this box 
          fx = r_e*cos(theta0*pi/180)*pi/180   ! factors between degree
          fy = r_p*pi/180
          delta_lon=delta_x/fx
          delta_lat=delta_y/fy

          do inode=1,nnode
             read(*,'(a100)') line_of_data
             read(line_of_data,*) node(1,inode),node(2,inode)
             write(65,'(2e20.7)') node(1,inode),node(2,inode)
          enddo
          if (use_fine_grid) then
             write(65,*)
             ny=int((ymax-ymin)/delta_lat)+2
             do iy=0,ny-1
                yg=ymin-0.5*delta_lat+iy*delta_lat
                n_int=0
   
                ! start scanning pgon for intersections 
                do inode=1,nnode-1
                   x1=node(1,inode)
                   y1=node(2,inode)
                   x2=node(1,inode+1)
                   y2=node(2,inode+1)
                   if (abs(y2-yg).gt.1e-014) then
                      xx=(y2-yg)*(y1-yg)          
                      if(xx.le.1e-014) then  ! level intersected 
                         n_int=n_int+1
                         x_int(n_int)= (x2*(yg-y1) + x1*(y2-yg))/(y2-y1)
                      endif
                   endif
                enddo ! inode 
   
                if(n_int.gt.0) then 
                   if (mod(n_int,2).ne.0) then     ! odd number of intersections: PROBLEMS !!
                      write(0,'(a,E15.5,a)') 
     &                       '! Odd number of intersections of contour'
     &                  //    ' and gridline (Y = ',yg,')! '
                   endif
                   call ssort(x_int,n_int) ! sort first n_int elements ascendingly 
                   nx=int((xmax-xmin)/delta_lon)+2
                   do ix=0,nx-1            ! run over gridpoints 
                      xg=xmin-0.5*delta_lon+ix*delta_lon
                      do i=1,n_int
                         if(x_int(i).gt.xg) exit 
                      enddo                ! intersections 
                      i=i-1                ! i holds the number of intersections passed 
                      if (mod(i,2).ne.0) then 
                          ipt=ipt+1 
                          gridpt(1,ipt)=xg ! store x-coordinate 
                          gridpt(2,ipt)=yg ! store y-coordinate 
                      endif 
                   enddo ! ix 
                endif
             enddo ! iy 
          endif ! use_fine_grid
       enddo ! reading from polygon definitions 
 666   continue
       close(34)
       close(35)
       close(55)
       close(65)
       close(66)

       ierror=nf_close(ncid)
       if (verbose) write(0,*) 'NetCDF file closed' 

C tot hier 

       open(66,file='used_stns.out')
       do istn=1,ngd1*ngd2
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
       write(0,*)
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
          knl(i) = (1-alpha**2)**2   ! fourth order function 
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
              k = mod(l+i-1+nx,nx) + 1  ! mod can't handle negative
             res(l) = res(l)+x(k)*y(i)
           enddo ! i 
       enddo ! j 
       return
       end 

       integer function date2index(datum,leap)
       implicit none 
       ! returns the order of the day of datum from 19000101 (the first day)
       ! leap=1 take leap years into account 
       ! leap=0 do NOT count leap years: each year 365 days 
       integer dayzero1(12)     ! day zero in each month 
       integer iy,day,mnth,datum
       integer dayzero2(12)     ! day zero in each month 
       integer dayzero3(12)     ! day zero in each month 
       integer nleap 
       integer leap 
       data dayzero1 /0,31,59,90,120,151,181,212,243,273,304,334/ ! ordinary year 
       data dayzero2 /0,31,60,91,121,152,182,213,244,274,305,335/ ! leapyear 
       data dayzero3 /0,31,60,91,121,152,182,212,243,273,304,334/ ! shifted 
       iy=datum/10000 - 1900                                      ! 1900 is not a leap year 
       nleap=(iy-1)/4                                             ! number of leap years passed since 1900
       mnth=mod(datum,10000)/100
       day=mod(datum,100)
       if(leap.gt.0) then 
          if(mod(iy,4).eq.0 .and. mod(iy,100).ne.0) then
             date2index=dayzero2(mnth)+day+365*(iy)+nleap         ! leap year 
          else 
             date2index=dayzero1(mnth)+day+365*(iy)+nleap         ! ordinary year 
          endif 
       else 
          if(mod(iy,4).eq.0 .and. mod(iy,100).ne.0) then
             date2index=dayzero3(mnth)+day+365*(iy)               ! leap year 
          else 
             date2index=dayzero1(mnth)+day+365*(iy)               ! ordinary year 
          endif 
       endif 

       return 
       end 

!JB>
       integer function index2date(dayssince, startdate)
       implicit none 
        ! (inverse of date2index with leap=1)
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

       index2date = yyyy*10000 + mm*100 + dd
       
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
     &  r* acos(cos(the1)*cos(the2)*cos(phi2-phi1)
     &              + sin(the1)*sin(the2))
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

        Subroutine handle_error (ierror,pre,post)
          include 'netcdf.inc'
          integer ierror
          character*(*) pre, post
          if (ierror.ne.nf_noerr ) then
             write(0,*) pre//' '//trim(nf_strerror(ierror))//' : '//post
             stop "stopped in handle_error"
          end if
        end

        logical function udunit2date(udunitstr,origin)
        implicit none
        character*(*) udunitstr
        integer dd, mm, yyyy
        integer origin, i, ios
        character*(30) sorigindate, stunit 
        character*(10) dummy
        logical success
        success = .False.
        read(udunitstr,*,IOSTAT=ios) stunit, dummy, sorigindate
        do i=1,len_trim(sorigindate)
           if (sorigindate(i:i).eq.'-') sorigindate(i:i) = ' '
        enddo
        read(sorigindate,*) yyyy,mm,dd
        origin = yyyy*10000+mm*100+dd
        udunit2date = .True.
        return
        end
         

        include 'routines.f' 
