        program petf
        implicit none
        character*(25) f_mpv, f_mtm
        character*(50) regel
        double precision mtm(12), mpv(12), tdaily, etf
        parameter (etf=0.05)
        integer ierr, im, ios, date
        logical argstring
        integer instr
        do im=1,12
           mtm(im) = -999.d0
           mpv(im) = -999.d0
        enddo 
        if (.not.argstring('-mtm', 'mtm.dat', f_mtm)) then
          stop 'No monthly temperature file (-mtm argument)'
        endif

        if (.not.argstring('-mpv', 'mpv.dat', f_mpv)) then
          stop 'No monthly pv file (-mpv argument)'
        endif

        open(22,file=trim(f_mtm),status='OLD') 
        do while(.true.)
           read(22,'(a)',end=222) regel
           regel=regel(1:instr(regel//'#','#')-1)
           if (len_trim(regel).gt.0) then
              read(regel,*,iostat=ios) im, mtm(im) 
           endif
        enddo
 222    continue
        close(22)
        open(33,file=trim(f_mpv),status='OLD') 
        do while(.true.)
           read(33,'(a)',end=333) regel
           regel=regel(1:instr(regel//'#','#')-1)
           if (len_trim(regel).gt.0) then
              read(regel,*,iostat=ios) im, mpv(im) 
           endif
        enddo
 333    continue
        close(33)
!       do im=1,12
!          write(0,'(i4,2f12.5)') im, mtm(im),mpv(im) 
!       enddo

        ! --------------------- process file ------------------------
        do while(.True.)
           read(*,'(a)',end=666) regel
           regel=regel(1:instr(regel//'#','#')-1)
           if (len_trim(regel).gt.0) then
              read(regel,*,iostat=ios) date, tdaily 
              im = int(mod(date,10000)/100)
              write(*,'(i10,2f10.3)') 
     &           date, 
     &             mpv(im)*(1.d0+(tdaily-mtm(im))*etf)
           endif
        enddo
 666    continue
        ! --------------------- process file ------------------------

        end
        
        
        

        
        include 'routines.f'
