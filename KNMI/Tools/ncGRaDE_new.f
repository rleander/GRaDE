      program ncGRaDE
      implicit none 
      include 'netcdf.inc' 
      integer ncid, ierror 
      integer lastchar, firstchar, instr
      integer dims(2)
      integer time_dimid,station_dimid,cli_dimid,cln_dimid
      integer time_varid,lon_varid,lat_varid
      integer station_varid,stnam_varid,staid_varid
      integer ios
      integer ndx2date, date2ndx  ! RL666
      integer refdate, startdate, enddate
      integer istn, istnmax, ndxmax, ndx, ndx0, ndx1, refndx 
      integer MAXSTN, MAXTIME
      parameter (MAXSTN=70,MAXTIME=366*2000)
      real opp(MAXSTN)
      double precision times(MAXTIME)
      character*(50) simname
      character*(30) naam(MAXSTN)
      character*(30) HBVid(MAXSTN)
      logical file_exists(MAXSTN)
      logical argstring, arglogical
      integer argint 
      logical include_meta, verbose
      character*(50) fmeta, fnc, fnin, fstn 
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
      character*(500) title 
      character*(100) institute 
      character*(500) source 
      character*(500) references 
      character*(500) disclaimer 
      character*(100) version 
      character*(500) features 
      character*(500) reference_data 
      character*(100) Conventions 
      character*(1000) comment 

      character*(10000) attbuffer               ! buffer for attributes read from file 
      character*(500) att_name
      character*(1000) att_text, regel 
      integer sep 


      ! HANDLE cmdlargs
      if(.not.argstring('-sim','',simname)) then
         write(0,*) 'No simulation name specified (-sim) !'  
      endif
      if(.not.argstring('-nc', 'out.nc',fnc)) then 
        write(0,*) 'No Netcdf-file specified !'  
        write(0,*) 'Defaulting to ''out.nc'' ...'
      endif 
      verbose = arglogical('-v') 
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
      refdate = argint('-refdate',startdate)    ! days since .... time reference
      write(time_units,'(a,i4,a,i2.2,a,i2.2,a)') 
     & 'days since ', int(refdate/10000),
     &            '-',int(mod(refdate,10000)/100),
     &            '-',mod(refdate,100), ' 00:00:00.0'     
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
      call hndl(nf_def_dim(ncid,'time',NF_UNLIMITED,time_dimid))
      call hndl(nf_def_dim(ncid,'stations',istnmax,station_dimid))
      call hndl(nf_def_dim(ncid,'char_leng_id',30,cli_dimid))
      call hndl(nf_def_dim(ncid,'char_leng_name',30,cln_dimid))


      ! DEFINE VARIABLES
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

!      dims(1)=station_dimid
!      call hndl(nf_def_var(ncid,'stations',NF_DOUBLE,1,dims,
!     &                                          station_varid))
!      call hndl(nf_def_var(ncid,'latitude',NF_FLOAT,1,dims,
!     &                                              lat_varid))
!      call hndl(nf_put_att_text(ncid,lat_varid,'standard_name',
!     &          len(latitude_standard_name),latitude_standard_name))
!      call hndl(nf_put_att_text(ncid,lat_varid,'long_name',
!     &                  len(latitude_long_name),latitude_long_name))
!      call hndl(nf_put_att_text(ncid,lat_varid,'units',
!     &                  len(latitude_units),latitude_units))
!      call hndl(nf_def_var(ncid,'longitude',NF_FLOAT,1,dims,
!     &                                              lon_varid))
!      call hndl(nf_put_att_text(ncid,lon_varid,'standard_name',
!     &          len(longitude_standard_name),longitude_standard_name))
!      call hndl(nf_put_att_text(ncid,lon_varid,'long_name',
!     &                  len(longitude_long_name),longitude_long_name))
!      call hndl(nf_put_att_text(ncid,lon_varid,'units',
!     &                  len(longitude_units),longitude_units))

      dims(2)=station_dimid
      dims(1)=cli_dimid
      call hndl(nf_def_var(ncid,'HBV_id',NF_CHAR,2,dims,
     &                                              staid_varid))
      dims(1)=cln_dimid
      call hndl(nf_def_var(ncid,'HBV_names',NF_CHAR,2,dims,
     &                                              stnam_varid))

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

      call hndl(nf_def_var(ncid,e_name,NF_FLOAT,2,dims,e_varid))
      call hndl(nf_put_att_text(ncid,e_varid,'standard_name',
     &                          len(e_standard_name),e_standard_name))
      call hndl(nf_put_att_text(ncid,e_varid,'long_name',
     &                                  len(e_long_name),e_long_name))
      call hndl(nf_put_att_text(ncid,e_varid,'units',
     &                                  len(e_units),e_units))
      call hndl(nf_put_att_real(ncid,e_varid,'FillValue',NF_FLOAT,
     &                                             1,-9999.0))

      ! OPEN TEXTFILE WITH METADATA (GLOBAL ATTRIBUTES)
      if(include_meta) then 
        open(12,file=trim(fmeta))
        do while(.True.)
          read(12,'(a1000)',end=112) regel
          if (len_trim(regel).gt.0) then
             sep=instr(regel,'=')
             att_name=regel(instr(regel,':')+1:sep-1)
             att_text=regel(sep+1:1000)
             if (verbose) write(*,*) trim(att_name)//':'//trim(att_text)
             call hndl(nf_put_att_text(ncid,NF_GLOBAL,
     &         att_name(firstchar(att_name):lastchar(att_name)),
     &         lastchar(att_text)-firstchar(att_text)+1,
     &         att_text(firstchar(att_text):lastchar(att_text))))
           endif
        enddo 
 112    continue 
        close(12)
      endif 

      ! CLOSE NETCDF HEADER, OPEN DATA
      call hndl(nf_enddef(ncid))

      ! WRITE STATION/CATCHMENT NAMES AND ID's
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
            write(0,'(i4,a30,a30)') istn, HBVid(istn), naam(istn)
         enddo
      endif    

      ! WRITE P,T,E data
      write(0,*) 'Writing '//trim(p_name)
      ndxmax = 0
      call putdata(ncid,p_varid,p_fillvalue,simname,
     &             refdate,startdate,enddate,p_ext,2,istnmax,ndxmax)
      write(*,*) ndxmax, ndx2date(ndxmax-1,refdate)
      write(0,*) 'Writing '//trim(t_name)
      call putdata(ncid,t_varid,t_fillvalue,simname,
     &             refdate,startdate,enddate,t_ext,1,istnmax,ndxmax)
      write(*,*) ndxmax, ndx2date(ndxmax-1,refdate)
      write(0,*) 'Writing '//trim(e_name)
      call putdata(ncid,e_varid,e_fillvalue,simname,
     &             refdate,startdate,enddate,e_ext,2,istnmax,ndxmax)
      write(*,*) ndxmax, ndx2date(ndxmax-1,refdate)
      

      ! WRITE TIMES, just from 0 to ndxmax - 1 (consistent with earlier files)
      write(0,*) 'Writing time'
      refndx = date2ndx(refdate,1)
      ndx0 = date2ndx(startdate,1) - refndx + 1      ! what goes into the time variable at the first entry
      ndx1 = date2ndx(enddate,1) - refndx + 1        ! what goes into the time variable at the first entry
      do ndx = ndx0, ndx1
         times(ndx-ndx0+1) = ndx*1.0
      enddo
      call hndl(nf_put_vara_double(ncid, time_varid,
     &   (/1/),(/ndx1-ndx0+1/),times(1:ndxmax)))

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
     &                   refdate,startdate,enddate,ext,ndec,
     &                                        istnmax,ndxmax)
       implicit none
       include 'netcdf.inc' 
       integer ncid, ndec, varid
       integer ndxmax
       real fillvalue
       character*(*) ext
       integer refdate,startdate,enddate
       integer istnmax 
       character*(*) simname
       character*(3) sstn
       integer istn, ios, f
       integer MAXSTN
       integer maxndx
       parameter (MAXSTN=70)
       logical file_opened(MAXSTN)
       real datavalue(MAXSTN)
       character*(100) regel
       character*(50) fnamestr
     
       integer date2ndx, ndx, ndx0, refndx, datum, datum0

       f = 10**ndec
       refndx = date2ndx(refdate,1)
       ndx0= date2ndx(startdate,1) - refndx + 1      ! what goes into the time variable at the first entry
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
                datavalue(istn) = int(datavalue(istn)*f+0.5)/dble(f)
                if (datum0.eq.0) then
                   ndx = date2ndx(datum,1) - refndx + 1          ! what goes into the time variable 
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
          if (datum.ge.startdate .and. datum.le.enddate) then           ! only write within requested time window
             if (ndx.gt.0) then                                         ! valid index into the time dimension
                call hndl(nf_put_vara_real(ncid, varid,
     &              (/1,ndx-ndx0+1/),(/istnmax,1/),datavalue(1:istn)))
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

      subroutine putndata(ncid,varids,fillvalue,simname,
     &             refdate,startdate,enddate,fu0,fu1,nvar,ndec,nstn,
     &                                        istnmax,ndxmax)
      ! Writes nvar variables to a netcdf file
      ! varids should be in varids(*)
      ! file units should be running consecutively from 
      ! fu0 ... fu1 in nvar blocks of nstn 
      ! so, fu1-fu0+1 = nvar*nstn
      ! file unit for var i, stn j = fu0+(i-1)*nstn + j-1
      implicit none
      include 'netcdf.inc' 
      integer ncid, ndec(*), nvar, nstn 
      integer varids(*), fu0, fu1
      integer ndxmax
      real fillvalue
      integer refdate,startdate,enddate
      integer istnmax 
      character*(*) simname
      character*(3) sstn
      integer istn, ios, f
      integer MAXSTN
      integer maxndx
      parameter (MAXSTN=70)
      logical file_opened(MAXSTN)
      real datavalue(MAXSTN)
      character*(100) regel
      character*(50) fnamestr
     
      integer date2ndx, ndx, ndx0, refndx, datum, datum0
      integer i         ! index over variables 
      integer j         ! index over stations 

      f = 10**ndec(i)
      refndx = date2ndx(refdate,1)
      ndx0= date2ndx(startdate,1) - refndx + 1      ! what goes into the time variable at the first entry
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
               datavalue(istn) = int(datavalue(istn)*f+0.5)/dble(f)
               if (datum0.eq.0) then
                  ndx = date2ndx(datum,1) - refndx + 1          ! what goes into the time variable 
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
         if (datum.ge.startdate .and. datum.le.enddate) then            ! only write within requested time window
            if (ndx.gt.0) then                                          ! valid index into the time dimension
               call hndl(nf_put_vara_real(ncid, varid,
     &             (/1,ndx-ndx0+1/),(/istnmax,1/),datavalue(1:istn)))
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




