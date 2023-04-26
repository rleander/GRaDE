      program joinasciiseries
      implicit none
      integer MAXSTN
      parameter(MAXSTN=300)
      logical file_opened(MAXSTN)    
      integer istnmax, istn
      integer datum, datum0
      integer ios
      real*8 datavalue, fillvalue
      character*(50) fnamestr
      character*(50) simname
      character*(100) regel
      character*(3) sstn
      character*(8) ext
      integer argint
      logical argstring
      real*4 argreal

      if (.not.argstring('-e','',ext)) then
         stop 'No extension ''-e'' specified'
      endif
      simname = ''
      if (.not.argstring('-s','',simname)) then
         write(0,*)  'No simulation name ''-s'' specified'
      endif
      fillvalue = argreal('-fill',-9999.0)
      istnmax   = argint('-nmax',1)
      
!     simname = 'VechtTest_3dmem_w=61_long'
!     fillvalue = -9999.0

      ! Open all source files
      file_opened = .False.
      do istn=1,istnmax
         write(sstn,'(i3.3)') istn
         if (len_trim(simname).gt.0) then
            write(fnamestr,'(a)') 'area'//trim(sstn)//'_'           ! area034_sim.rr
     &              //trim(simname)//'.'//trim(ext)
         else
            write(fnamestr,'(a)') 'area'//trim(sstn)                ! area034.rr
     &                             //'.'//trim(ext)
         endif
         open(10+istn,file=trim(fnamestr),status='OLD',iostat=ios)
         if (ios.eq.0) then
            file_opened(istn)=.True.
            write(0,*) trim(fnamestr)//' opened' 
         else
            file_opened(istn)=.False.
            write(0,*) 'Could not open '//trim(fnamestr)
         endif
      enddo

      ! Write lines
      do while (.true.)
         datum0 = 0
         do istn=1,istnmax
            regel = ''
            do while (len_trim(regel).eq.0 .or. regel(1:1).eq.'#')
               read(10+istn,'(a)',end=616) regel
            enddo
            datavalue = fillvalue
            read(regel,*,err=516) datum, datavalue
 516        continue
            if (datum0.eq.0) then
               datum0 = datum
               write(*,'(i15$)') datum      
            else
               if (datum.ne.datum0) then ! datum mismatch !
                  continue
               endif
            endif
            write(*,'(f12.6$)') datavalue
         enddo
         write(*,*)
      enddo
 616  continue
      ! Close all source files
      do istn=1,istnmax
         if (file_opened(istn)) close(10+istn)
      enddo
      end program




      include 'routines.f'
