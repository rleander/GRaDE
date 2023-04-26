      program ncGRaDE
      implicit none 
      include 'netcdf.inc' 
      integer ncid, ierror 
      integer lastchar, firstchar, instr
      integer dims(2)
      integer time_dimid,station_dimid,cli_dimid,cln_dimid
      integer time_varid,lon_varid,lat_varid
      integer station_varid,stnam_varid,staid_varid
      integer p_varid, t_varid, e_varid 

      logical argstring
      logical include_meta
      character*(50) fmeta, fnc 
      
      character*(*) time_units,time_calendar,time_standard_name
      character*(*) time_long_name
      parameter(time_units='days since 0001-01-01 00:00:00.0')
      parameter(time_long_name='DaysSince001-01-01_000000.0')
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

      character*(*) p_units 
      real p_FillValue
      parameter(p_units='mm/day')
      parameter(p_FillValue=-9999.0)

      character*(*) t_units 
      real t_FillValue
      parameter(t_units='degrees C')
      parameter(t_FillValue=-9999.0)

      character*(*) e_units 
      real e_FillValue
      parameter(e_units='mm/day')
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
      if(.not.argstring('-nc', 'shit.nc',fnc)) then 
        write(0,*) 'No Netcdf-file specified !'  
        write(0,*) 'Defaulting to ''shit.nc'' ...'
      endif 
      include_meta=argstring('-meta','',fmeta)  ! include meta-data in the file as 
                                                ! global attributes if availabe 

      ! OPEN NETCDF 
      call hndl(nf_create(fnc(1:lastchar(fnc)),NF_CLOBBER,ncid))
        
      ! DEFINE DIMENSIONS 
      call hndl(nf_def_dim(ncid,'time',NF_UNLIMITED,time_dimid))
      call hndl(nf_def_dim(ncid,'stations',15,station_dimid))
      call hndl(nf_def_dim(ncid,'char_leng_id',10,cli_dimid))
      call hndl(nf_def_dim(ncid,'char_leng_name',30,cln_dimid))


      ! DEFINE VARIABLES
      dims(1)=time_dimid
      call hndl(nf_def_var(ncid,'time',NF_DOUBLE,1,dims,time_varid))
      call hndl(nf_put_att_text(ncid,time_varid,'units',
     &                  len(time_units),time_units))
      call hndl(nf_put_att_text(ncid,time_varid,'units',
     &                  len(time_calendar),time_calendar))
      call hndl(nf_put_att_text(ncid,time_varid,'standard_name',
     &                  len(time_standard_name),time_standard_name))
      call hndl(nf_put_att_text(ncid,time_varid,'long_name',
     &                  len(time_long_name),time_long_name))

      dims(1)=station_dimid
      call hndl(nf_def_var(ncid,'stations',NF_DOUBLE,1,dims,
     &                                          station_varid))
      call hndl(nf_def_var(ncid,'latitude',NF_FLOAT,1,dims,
     &                                              lat_varid))
      call hndl(nf_put_att_text(ncid,lat_varid,'standard_name',
     &          len(latitude_standard_name),latitude_standard_name))
      call hndl(nf_put_att_text(ncid,lat_varid,'long_name',
     &                  len(latitude_long_name),latitude_long_name))
      call hndl(nf_put_att_text(ncid,lat_varid,'units',
     &                  len(latitude_units),latitude_units))

      call hndl(nf_def_var(ncid,'longitude',NF_FLOAT,1,dims,
     &                                              lon_varid))
      call hndl(nf_put_att_text(ncid,lon_varid,'standard_name',
     &          len(longitude_standard_name),longitude_standard_name))
      call hndl(nf_put_att_text(ncid,lon_varid,'long_name',
     &                  len(longitude_long_name),longitude_long_name))
      call hndl(nf_put_att_text(ncid,lon_varid,'units',
     &                  len(longitude_units),longitude_units))

      dims(2)=cli_dimid
      call hndl(nf_def_var(ncid,'station_id',NF_CHAR,2,dims,
     &                                              staid_varid))
      dims(2)=cln_dimid
      call hndl(nf_def_var(ncid,'station_names',NF_CHAR,2,dims,
     &                                              stnam_varid))

      dims(1)=station_dimid
      dims(2)=time_dimid
      call hndl(nf_def_var(ncid,'rainfall_rate',NF_FLOAT,2,dims,
     &                                                  p_varid))
      call hndl(nf_put_att_text(ncid,p_varid,'units',
     &                                  len(p_units),p_units))
      call hndl(nf_put_att_real(ncid,p_varid,'FillValue',NF_FLOAT,
     &                                             1,-9999.0))
      call hndl(nf_def_var(ncid,'air_temperature',NF_FLOAT,2,dims,
     &                                                  t_varid))
      call hndl(nf_put_att_text(ncid,t_varid,'units',
     &                                  len(t_units),t_units))
      call hndl(nf_put_att_real(ncid,t_varid,'FillValue',NF_FLOAT,
     &                                             1,-9999.0))
      call hndl(nf_def_var(ncid,'reference_evaporation_rate'
     &                                 ,NF_FLOAT,2,dims,e_varid))
      call hndl(nf_put_att_text(ncid,e_varid,'units',
     &                                  len(e_units),e_units))
      call hndl(nf_put_att_real(ncid,e_varid,'FillValue',NF_FLOAT,
     &                                             1,-9999.0))

      ! GLOBAL ATTRIBUTES 

      ! OPEN TEXTFILE WITH METADATA 
      if(include_meta) then 
        open(33,file='meta')
        do while(.True.)
          read(33,'(a1000)',end=333) regel
          sep=instr(regel,'=')
          att_name=regel(1:sep-1)
          att_text=regel(sep+1:1000)
C          write(*,*) att_name(firstchar(att_name):lastchar(att_name))
C     &        //':'
C     &        // att_text(firstchar(att_text):lastchar(att_text))
          call hndl(nf_put_att_text(ncid,NF_GLOBAL,
     &      att_name(firstchar(att_name):lastchar(att_name)),
     &      lastchar(att_text)-firstchar(att_text)+1,
     &      att_text(firstchar(att_text):lastchar(att_text))))
        enddo 
 333    continue 
        close(33)
      endif 

      ! CLOSE NETCDF
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

        integer function lastchar(s)
        implicit none
        character*(*) s
        integer i
        do i=len(s),1,-1
           if(s(i:i).ne.' ') then
              lastchar=i
              goto 047
           endif
        enddo ! i
 047    continue
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
 037    continue
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

        logical function argstring(prefix, default, result)
        implicit none
        integer jarg
        integer instr
        character*(*) default
        character*(*) result            !MAX 20 characters!
        character*(*) prefix
        character*50 sarg

        result=default
        argstring=.False.
        do jarg=1,iargc()
           call getarg(jarg,sarg)
           if(sarg(1:instr(sarg,' ')-1).eq.prefix) then
              call getarg(jarg+1,sarg)
              read(sarg,*,end=233,err=233) result
              argstring=.True.
           endif
 233       continue
        enddo
        return
        end



