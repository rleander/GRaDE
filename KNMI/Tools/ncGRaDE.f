      program ncGRaDE

!     201910 Jules Beersma (KNMI):
!     * MAXSTN = 134 (in stead of 70; i.e. large enough for HBV-Rhine)
!     * Reading/writing of evaporation data only with (new) command line option: "-ev"!
!     * New command line option "-shift" added (similar to "-offset" but shifts only the original NG years)
!     * Adds creation date to global attributes in NetCDF file (in stead of "simname")      

      implicit none 
      include 'netcdf.inc' 
      integer ncid, ierror
      integer lastchar, firstchar, instr
      integer dims(2)
      integer time_dimid,station_dimid,cli_dimid,cln_dimid
      integer time_varid,date_varid,lon_varid,lat_varid
      integer station_varid,stnam_varid,staid_varid, stn_varid
      integer ios
      integer days_since
      integer refyear, offset, shift, startdate, enddate
      integer istn, istnmax, ndxmax, ndx, ndx0, ndx1, refndx
      integer MAXSTN, MAXTIME
!     parameter (MAXSTN=70,MAXTIME=366*20000)
      parameter (MAXSTN=134,MAXTIME=366*20000)
      real opp(MAXSTN)
      logical do_netcdf, do_ascii
      double precision times(MAXTIME)
      integer dates (MAXTIME)
      character*(50) simname
      character*(30) naam(MAXSTN)
      character*(30) HBVid(MAXSTN)
      logical file_exists(MAXSTN)
      logical argstring, arglogical
      integer argint
      logical include_meta, verbose, include_evap
      character*(100) fmeta, fnc, fasc, fnin, fstn
      integer from(2), thru(2)
      
      character*(32) time_units
      character*(*) time_calendar,time_standard_name
      character*(*) time_long_name
!     parameter(time_units='days since 0001-01-01 00:00:00.0')
      parameter(time_long_name='time')
      parameter(time_calendar='gregorian')
      parameter(time_standard_name='time')

      character*(*) latitude_standard_name,latitude_long_name
      character*(*) latitude_units
      parameter(latitude_standard_name='latitude')
      parameter(latitude_long_name='Latitude values')
      parameter(latitude_units='degrees_N')

      character*(*) longitude_standard_name,longitude_long_name
      character*(*) longitude_units
      parameter(longitude_standard_name='longitude')
      parameter(longitude_long_name='Longitude values')
      parameter(longitude_units='degrees_E')

      character*(*) p_name,p_standard_name,p_long_name
      character*(*) p_units 
      integer p_varid
      real p_FillValue
      character*(*) p_ext
      parameter(p_ext='pr')
      parameter(p_name='rainfall_rate')
      parameter(p_standard_name='rainfall_rate')
      parameter(p_long_name='rainfall amounts per day')
      parameter(p_units='mm')
      parameter(p_FillValue=-9999.0)

      character*(*) t_name,t_standard_name,t_long_name
      character*(*) t_units 
      integer t_varid
      real t_FillValue
      character*(*) t_ext
      parameter(t_ext='tg')
      parameter(t_name='air_temperature')
      parameter(t_standard_name='air_temperature')
      parameter(t_long_name='air temperature values')
      parameter(t_units='degree_Celsius')
      parameter(t_FillValue=-9999.0)

      character*(*) e_name,e_standard_name,e_long_name
      character*(*) e_units 
      integer e_varid
      real e_FillValue
      character*(*) e_ext
      parameter(e_ext='ev')
      parameter(e_name='reference_evaporation_rate')
      parameter(e_standard_name='water_potential_evaporation_amount')
      parameter(e_long_name=
     &               'potential evapotranspiration amounts per day')
      parameter(e_units='mm')
      parameter(e_FillValue=-9999.0)

      ! character
!     character*(500) title
!     character*(100) institute
!     character*(500) source
!     character*(500) reference
!     character*(500) disclaimer
!     character*(100) version
!     character*(500) features
!     character*(500) reference_data
!     character*(100) Conventions
!     character*(1000) comment

      character*(10000) attbuffer               ! buffer for attributes read from file 
      character*(500) att_name
      character*(1000) att_text, regel

      character*(8)   yyyymmdd
      character*(10)  datum
      integer sep 


      ! HANDLE cmdlargs
      if(.not.argstring('-sim','',simname)) then
         write(0,*) 'No simulation name specified (-sim) !'  
      endif
      do_netcdf=argstring('-nc', '',fnc)
      if (.not.do_netcdf) then
         write(0,*) 'No Netcdf-file specified, so no netcdf output!'  
      endif 
      do_ascii=argstring('-ascii', '',fasc)
      if (.not.do_ascii) then
         write(0,*) 'No Ascii-file specified, so no ascii output!'  
      endif 
      if (.not.(do_ascii .or. do_netcdf)) then
         write(0,*) 'No output specified, closing...'  
         stop ""
      endif 

      verbose = arglogical('-v') 
      include_evap = arglogical('-ev')          !JB20191017
      include_meta=argstring('-meta','',fmeta)  ! include meta-data in the file as 
                                                ! global attributes if availabe 
      if (.not.arglogical('-start')) then
         stop 'No startdate specified (-start) !'  
      endif
      startdate = argint('-start',19500101)     ! start date

      if (.not.arglogical('-end')) then
         stop 'No enddate specified (-end) !'  
      endif

      enddate = argint('-end',99999999)         ! stop date
      refyear = argint('-refyear',int(startdate/10000))    ! days since .... time reference
      offset  = argint('-offset',0)             ! shift of years added to the ALL the dates
      shift   = argint('-shift',0)              ! shift of years added to the dates in the inputfile ONLY (used e.g. to change 00010101 into 20010101)
                                                ! NOTE that offset and shift should be multiples of 4 years, otherwise leapyears and non-leapyears
                                                ! do not correspond to the calender anymore
      write(time_units,'(a,i4,a,i2.2,a,i2.2,a)') 
     & 'days since ', refyear,
     &            '-',01,
     &            '-',01, ' 00:00:00.0'     
      if (.not.argstring('-list','',fstn)) then
         write(0,*) 'No station/area list specified (-list) !'  
         stop 'FATAL!'
      endif
                                                ! global attributes if availabe 
      ! READ LIST OF STATIONS/CATCHMENTS
      do istn=1,MAXSTN
         file_exists(istn) = .False.
         naam(istn) = ''
         HBVid(istn) = ''
         opp(istn) = 0.0
      enddo
      istnmax = 0
      if (len_trim(fstn).gt.0) then
         open(10,file=trim(fstn),status='OLD')
         do while(.True.)
            read(10,'(a)',end=110) regel
            if (regel(1:1).eq.'#') cycle
            read(regel,*) istn, opp(istn),naam(istn), HBVid(istn)
            write(fnin,'(a,i3.3,a)') 'area', istn, '.all'
            inquire(file=fnin,exist=file_exists(istn))
            istnmax = max(istnmax,istn)
         enddo
 110     continue
         close(10)
      endif

      ! OPEN NETCDF 
      call hndl(nf_create(fnc(1:lastchar(fnc)),NF_CLOBBER,ncid))
        
      ! DEFINE DIMENSIONS 
      write(0,*) 'Dimensions ... '
      call hndl(nf_def_dim(ncid,'time',NF_UNLIMITED,time_dimid))
      call hndl(nf_def_dim(ncid,'stations',istnmax,station_dimid))
      call hndl(nf_def_dim(ncid,'char_leng_id',30,cli_dimid))
      call hndl(nf_def_dim(ncid,'char_leng_name',30,cln_dimid))

      ! DEFINE VARIABLES
      write(0,*) 'Variables ... '
      dims(1)=time_dimid
      call hndl(nf_def_var(ncid,'time',NF_DOUBLE,1,dims,time_varid))
      call hndl(nf_put_att_text(ncid,time_varid,'standard_name',
     &                  len(time_standard_name),time_standard_name))
      call hndl(nf_put_att_text(ncid,time_varid,'long_name',
     &                  len(time_long_name),time_long_name))
      call hndl(nf_put_att_text(ncid,time_varid,'units',
     &                  len(time_units),time_units))
      call hndl(nf_put_att_text(ncid,time_varid,'calendar',
     &                  len(time_calendar),time_calendar))

      ! DEFINE VARIABLES
      dims(1)=time_dimid
      call hndl(nf_def_var(ncid,'date',NF_INT,1,dims,date_varid))

!     dims(2)=station_dimid
!     dims(1)=cli_dimid
!     call hndl(nf_def_var(ncid,'HBV_id',NF_CHAR,2,dims,
!    &                                              staid_varid))
!     dims(1)=cln_dimid
!     call hndl(nf_def_var(ncid,'HBV_names',NF_CHAR,2,dims,
!    &                                              stnam_varid))

! ----------------------------------------------------------------
      dims(2)=station_dimid
      dims(1)=cli_dimid
      call hndl(nf_def_var(ncid,'station_id',NF_CHAR,2,dims,
     &                                              staid_varid))
       dims(1)=cln_dimid
      call hndl(nf_def_var(ncid,'station_names',NF_CHAR,2,dims,
     &                                              stnam_varid))
! ----------------------------------------------------------------




      dims(1)=station_dimid
      call hndl(nf_def_var(ncid,'stations',NF_DOUBLE,1,dims,
     &                                              stn_varid))

      dims(1)=station_dimid
      dims(2)=time_dimid

      call hndl(nf_def_var(ncid,p_name,NF_FLOAT,2,dims,p_varid))
      call hndl(nf_put_att_text(ncid,p_varid,'standard_name',
     &                          len(p_standard_name),p_standard_name))
      call hndl(nf_put_att_text(ncid,p_varid,'long_name',
     &                                  len(p_long_name),p_long_name))
      call hndl(nf_put_att_text(ncid,p_varid,'units',
     &                                  len(p_units),p_units))
      call hndl(nf_put_att_real(ncid,p_varid,'FillValue',NF_FLOAT,
     &                                             1,-9999.0))

      call hndl(nf_def_var(ncid,t_name,NF_FLOAT,2,dims,t_varid))
      call hndl(nf_put_att_text(ncid,t_varid,'standard_name',
     &                          len(t_standard_name),t_standard_name))
      call hndl(nf_put_att_text(ncid,t_varid,'long_name',
     &                                  len(t_long_name),t_long_name))
      call hndl(nf_put_att_text(ncid,t_varid,'units',
     &                                  len(t_units),t_units))
      call hndl(nf_put_att_real(ncid,t_varid,'FillValue',NF_FLOAT,
     &                                             1,-9999.0))

      if (include_evap) then !JB20191017
        call hndl(nf_def_var(ncid,e_name,NF_FLOAT,2,dims,e_varid))
        call hndl(nf_put_att_text(ncid,e_varid,'standard_name',
     &                            len(e_standard_name),e_standard_name))
        call hndl(nf_put_att_text(ncid,e_varid,'long_name',
     &                                    len(e_long_name),e_long_name))
        call hndl(nf_put_att_text(ncid,e_varid,'units',
     &                                    len(e_units),e_units))
        call hndl(nf_put_att_real(ncid,e_varid,'FillValue',NF_FLOAT,
     &                                               1,-9999.0))
      endif

      ! OPEN TEXTFILE WITH METADATA (GLOBAL ATTRIBUTES)
      write(0,*) 'Meta data ... '
      if(include_meta) then 
        open(12,file=trim(fmeta))
        do while(.True.)
          read(12,'(a1000)',end=112) regel
          if (len_trim(regel).gt.0) then
             sep=instr(regel,'=')
             att_name=regel(instr(regel,':')+1:sep-1)
             att_text=regel(sep+1:1000)
             if (verbose) write(0,*) trim(att_name)//':'//trim(att_text)
             call hndl(nf_put_att_text(ncid,NF_GLOBAL,
     &         att_name(firstchar(att_name):lastchar(att_name)),
     &         lastchar(att_text)-firstchar(att_text)+1,
     &         att_text(firstchar(att_text):lastchar(att_text))))
           endif
        enddo 
 112    continue 
        close(12)
        ! Add simulation name to global attributes
!       call hndl(nf_put_att_text(ncid,NF_GLOBAL,'simname',
!    &            len_trim(simname),simname))
        ! Add creation date to global attributes
        call date_and_time(yyyymmdd)
        datum = yyyymmdd(7:8)//"-"//yyyymmdd(5:6)//"-"//yyyymmdd(1:4)
        call hndl(nf_put_att_text(ncid,NF_GLOBAL,'NetCDF_created',
     &            len_trim(datum),datum))
      endif 

      ! CLOSE NETCDF HEADER, OPEN DATA
      call hndl(nf_enddef(ncid))

      ! WRITE STATION/CATCHMENT NAMES AND ID's
      write(0,*) 'Stations ... '
      if (len_trim(fstn).gt.0) then
         do istn=1,istnmax
            from(1) = 1
            from(2) = istn
            thru(1) = len_trim(HBVid(istn))
            thru(2) = 1
            call hndl(nf_put_vara_text(ncid, staid_varid,from,thru,
     &               HBVid(istn)))
            thru(1) = len_trim(naam(istn))
            call hndl(nf_put_vara_text(ncid, stnam_varid,from,thru,
     &               naam(istn)))
            call hndl(nf_put_vara_double(ncid, stn_varid,
     &               (/istn/),(/1/),dble(istn)))
            write(0,'(i4,a,a30,a30)') istn,':',HBVid(istn), naam(istn)
         enddo
      endif    

      ! WRITE P,T and optional E data
      if (include_evap) then !JB20191017
         write(0,*) 'P, T, E data ... '
      else
         write(0,*) 'P, T data ... '
      endif
      write(0,*) 'Writing '//trim(p_name)
      ndxmax = 0
      call putdata(ncid,p_varid,p_fillvalue,simname,refyear,startdate,
     &                enddate,offset,shift,p_ext,2,istnmax,ndxmax,dates,times)
      write(0,*) 'Writing '//trim(t_name)
      call putdata(ncid,t_varid,t_fillvalue,simname,refyear,startdate,
     &                enddate,offset,shift,t_ext,2,istnmax,ndxmax,dates,times)
      if (include_evap) then !JB20191017
        write(0,*) 'Writing '//trim(e_name)
        call putdata(ncid,e_varid,e_fillvalue,simname,refyear,startdate,
     &               enddate,offset,shift,e_ext,2,istnmax,ndxmax,dates,times)
      endif
        
      ! WRITE TIMES, just from 0 to ndxmax - 1 (consistent with earlier files)
      write(0,*) 'Writing time'
      ndx0 = days_since(startdate+offset*10000,refyear) 
      ndx1 = days_since(enddate+offset*10000,refyear) 
      do ndx = ndx0, ndx1
         times(ndx-ndx0+1) = ndx-1             ! times(1..ndx1-ndx0+1) = ndx0 .. ndx1
      enddo
      call hndl(nf_put_vara_double(ncid, time_varid,
     &   (/1/),(/ndx1-ndx0+1/),times(1:ndx1-ndx0+1)))
      call hndl(nf_put_vara_int(ncid, date_varid,
     &   (/1/),(/ndx1-ndx0+1/),dates(1:ndx1-ndx0+1)))

      ! CLOSE NETCDF
      write(0,*) 'Closing output'
      call hndl(nf_close(ncid))
      end 

      subroutine hndl(errno)
         include 'netcdf.inc'
         integer errno
         if (errno.ne.nf_noerr ) then
           write(0,*) nf_strerror(errno)
           stop "stopped in handle_error"
         end if
       end

      subroutine putdata(ncid,varid,fillvalue,simname,
     &              refyear,startdate,enddate,offset,shift,ext,ndec,
     &                                  istnmax,ndxmax,dates,times)
      implicit none
      include 'netcdf.inc' 
      integer ncid, ndec, varid
      integer ndxmax
      real fillvalue
      character*(*) ext
      integer refyear,startdate,enddate,offset,shift
      integer dayshift,mmdd
      integer istnmax 
      character*(*) simname
      character*(3) sstn
      integer istn, ios, f
      integer MAXSTN
      integer maxndx
!     parameter (MAXSTN=70)
      parameter (MAXSTN=134)
      logical file_opened(MAXSTN)
      real datavalue(MAXSTN)
      character*(100) regel
!     character*(50) fnamestr
      character*(100) fnamestr
      integer dates(*) 
      real*8 times(*)
     
      integer days_since, ndx, ndx0, refndx, datum, datum0
      integer ndx2date !JB20191011

      f = 10**ndec
      ndx0= days_since(startdate+offset*10000,refyear)
      do istn=1,istnmax
         write(sstn,'(i3.3)') istn
         if (len_trim(simname).gt.0) then
            write(fnamestr,'(a)') 'area'//trim(sstn)//'_'
     &              //trim(simname)//'.'//trim(ext)
         else
            write(fnamestr,'(a)') 'area'//trim(sstn)
     &                             //'.'//trim(ext)
         endif
         open(10+istn,file=trim(fnamestr),status='OLD',iostat=ios)
         if (ios.eq.0) then
            file_opened(istn)=.True.
         else
            file_opened(istn)=.False.
            write(0,*) 'Could not open '//trim(fnamestr)
         endif
      enddo
      maxndx = ndxmax
        dayshift = 0
      do while(.True.)
         do istn=1,istnmax
            datavalue(istn)=fillvalue
         enddo
         datum0 = 0
         do istn=1,istnmax
            if (file_opened(istn)) then
               regel=''
               do while (len_trim(regel).eq.0 .or. regel(1:1).eq.'#')
                  read(10+istn,'(a)',end=616) regel
               enddo
               read(regel,*) datum, datavalue(istn)
!              datum=datum+offset*10000
               datum=datum+(offset+shift)*10000  !JB20191008
               datavalue(istn) = int(datavalue(istn)*f+0.5)/dble(f) !afronden op 2 decimalen
               if (datum0.eq.0) then
                  ndx = days_since(datum,refyear)
                  mmdd = mod(datum,10000)
                  if (ndx.lt.0 .and. mmdd.eq.229) then
                     dayshift = 1
                     ndx=-ndx
                  elseif (dayshift.ne.0 .and. mmdd.ge.731) then
                     dayshift = 0
                     ndx=-ndx
                  else
                     ndx = ndx+dayshift
                  endif
                  datum0 = datum
               else
                  if (datum.ne.datum0) then
                     ! datum mismatch !
                     continue
                  endif
               endif
            else
               write(0,'(a,i3.3,a)') 'File #',istn,' not opened !'
            endif
         enddo
!        datum0 = datum                                              !JB20191011 TEST
         if (ndx .gt. 0) datum = ndx2date(ndx-1,refyear*10000+0101)  !JB20191011 schuift indien nodig datum 1 dag omhoog
!        if (datum.ne.datum0) write(0,*) datum0, datum, ndx          !JB20191011 TEST

         if (datum.ge.startdate+offset*10000 .and. 
     &       datum.le.enddate+offset*10000) then  ! only write within requested time window
            if (ndx.gt.0) then                    ! valid index into the time dimension !JB20191011 ndx only negative for 'the 0731 to be removed' 
!                                                   datum, ndx, ' date' ! RL666
!                                                   datavalue(1) = ndx  ! RL666
               call hndl(nf_put_vara_real(ncid, varid,
!    &             (/1,ndx-ndx0+1/),(/istnmax,1/),datavalue(1:istn)))    !JB20191009 gives run-error when istnmax = MAXSTN (which is logical) 
     &             (/1,ndx-ndx0+1/),(/istnmax,1/),datavalue(1:istnmax))) !JB20191009
               dates(ndx-ndx0+1) = datum
               times(ndx-ndx0+1) = ndx - 1.0
               maxndx = max(ndx,maxndx)
            endif
         endif
       enddo
 616   do istn=1,istnmax
          if (file_opened(istn)) close(10+istn)
       enddo
       ndxmax=maxndx
       return
       end


       include 'routines.f'




