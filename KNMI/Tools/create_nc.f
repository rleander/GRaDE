       program create_nc
       ! purpose : put P,E,T series for catchments in netCDF format 

       implicit none 
       include 'netcdf.inc'     

       integer ncid
       
C      OPEN NETCDF DATA 
       if(.not.argstring('-nc','',ncfile)) then                     !JB input NetCDF filename; GEEN default  
          stop 'No NETCDF-file name provided'
!      verbose = arglogical('-v')

       call handle_error ( nf_create(trim(ncfile),nf_nowrite, ncid),

!      ngd1 = field_dimlen(1)
!      ngd2 = field_dimlen(2)
!      ntime = field_dimlen(3)
       ! Therefore, force the user to specify the names on the commandline, no defaults.

       nf_enddef(ncid)
       nf_close(ncid)
       end

       include 'routines.f' 
