        program censor
        implicit none

        ! reads values from collumn cnr of STDIN 
        ! writes to STDOUT : 
        ! standardized gumbel variate, value, value of gumbelfit, sigma, mu 
        ! purpose : Read in samples x_i from a stdin and 
        ! estimate the a and u in their distribution function 

        ! F(x)=exp(-exp(-(x-mu)/sigma))	(GUMBEL: x = mu + alpha*y)  or 

        ! using maximum-likelyhood estimate
        integer ix,maxx,nx
        parameter(maxx=200000)
        real xi(maxx),yi(maxx) ! array with variables
        integer idum 
        integer iarg 
        integer cnr 			! column number 
        real euler			! Euler number 
        real Tthreshold1, Tthreshold2   ! Threshold return priod from which to apply weissman

        real Tr, Tmax			! return period 
        real a, b                       ! linear coefficients from least-squares
        real yy, dy, ymax
        integer nstep

        ! CENSORING 
        real sigma_hat,mean_xk
        real mu_hat
        integer k1, k2, n1, n2 
        logical left_fixed

        character*(6) paramstr 
        parameter(euler=0.577)	

        left_fixed = .True.
        
        if (iargc().lt.2) then 
           write(*,*) 
           write(*,*) 'syntax : censor <column> <k>' 
           write(*,*) 'purpose: censored Gumbel dist. on k largest'
           write(*,*) 'reads from stdin' 
           write(*,*) 'writes to stdout' 
           write(*,*) 
           write(*,*) 
           write(*,*) ' T1 T2 T3 ... are return periods for which'
           write(*,*) ' the estimated return levels are to be'
           write(*,*) ' calculated. This output is written to STDERR'
           write(*,*) 
           stop ''
        else 
            call getarg(1,paramstr)
            read(paramstr,*) cnr 
            call getarg(2,paramstr)
            k1 = -1
            read(paramstr,*,err=506) k1 
 506        if (k1<0) then
               read(paramstr,*) Tthreshold1
            endif
            call getarg(3,paramstr)
            k2 = -1
            read(paramstr,*,err=516) k2
 516        if (k2<0) then
               read(paramstr,*) Tthreshold2
            endif
            call getarg(4,paramstr)
            read(paramstr,*,err=516) Tmax
            call getarg(5,paramstr)
            read(paramstr,*,err=516) nstep
        endif 


        ! read in values 
        ix=1
        do 
 206       continue
           read(*,*,err=206,end=106) (xi(ix),idum=1,cnr)
           ix=ix+1
        enddo 
 106    continue 
        nx=ix-1
        if (k1<0) then
            k1 = nx - floor(0.3 + (nx+0.4)*(1.-1./Tthreshold1))
            write(0,'(a,f6.0,a,i0)') 'T1 = ',Tthreshold1,', k1 = ', k1 
        endif
        k1=max(2,min(k1,nx))

        if (k2<0) then
            k2 = nx - floor(0.3 + (nx+0.4)*(1.-1./Tthreshold2))
            write(0,'(a,f6.0,a,i0)') 'T2 = ',Tthreshold2,', k2 = ', k2 
        endif
        k2=max(2,min(k2,nx))

        ! sort values 
        call ssort(xi,nx)

        ! fill array with gumbel variates
        do ix=1,nx 
           yi(ix)=-log(-log((ix-0.3)/(nx+0.4)))
        enddo ! ix 

        ! ESTIMATE 
        n1 = nx-k1+1
        n2 = nx-k2

        if (left_fixed) then
           call lslinfit(yi(n1:n2),xi(n1:n2),a,b,k1-k2+1,
     &                           yi(n1),xi(n1),.True.)    ! left point fixed     
        else
           call lslinfit(yi(n1:n2),xi(n1:n2),a,b,k1-k2+1,0.,0.,.False.)   ! no fixed point     
        endif

        ! write results to stdout 
        ! gumbel var, x_i, fitted x_i, alpha estimate, u estimate

        write(*,*) '# Using a linear fit on the values between' 
        write(*,*) '# k1=', k1,' and k2=',k2,' from ',nx 
        write(*,*) '# 1 gumbel variate'
        write(*,*) '# 2 quantile'
        write(*,*) '# 3 gumbel fitted quantile'
        !write(*,*) '# 4 stddev of fitted quantile'
        
!	do ix=1,nx
!	   write(*,'(10f10.3)') 
!     &	         yi(ix),xi(ix)
!     &	        ,xi(nx-k+1)+sigma_hat
!     &	             *log(k/(nx*(1.-((ix-0.3)/(nx+0.4)) )))
!	enddo !ix 
        do ix=1,nx-k1
           write(*,'(10f10.3)') 
     &         yi(ix),xi(ix),xi(ix)
        enddo !ix 
        do ix=nx-k1+1,nx               ! replace the upper part by the fit for the second column
           write(*,'(10f10.3)') 
     &         yi(ix), b + a*yi(ix), xi(ix)  ! linear approximation
        enddo !ix 

        ! extrapolation
        ymax = -log(-log(1.-1./Tmax))          ! max std gum var
        dy = (ymax - yi(nx))/nstep
        do ix = 1, nstep
           yy = yi(nx) + dy*ix
           write(*,'(10f10.3)') yy, b + a*yy   ! linear extrapolation
        enddo

        end 

        subroutine ssort(x,n)
        implicit none 
        integer n,i,j
        real x(n),dummy
        do i=2,n
           dummy=x(i)
           do j=i,2,-1
              if(x(j-1).le.dummy) then 
                 exit 
              else 
                 x(j)=x(j-1)
                 x(j-1)=dummy
              endif 
           enddo ! j 
        enddo 
        return 
        end 

        subroutine syminv2d(m_in,m_out)
        implicit none 
        ! throw a non-singular symmetric 2x2 matrix (m_in) in and you get 
        ! the inverse out (m_out)
        real m_in(2,2),m_out(2,2)
        real a,b,c,f 
        a=m_in(1,1)
        b=m_in(2,2)
        c=m_in(1,2)
        m_out(2,2)= 1./(b- c*c/a)
        f=m_out(2,2)
        m_out(1,1)=(1.+c*c/a*f)/a
        m_out(1,2)=-c*f/a
        m_out(2,1)=m_out(1,2)
        return 
        end 


        subroutine lslinfit(x,y,a,b,nx,x0,y0,fixedpoint)        
        implicit none
        integer nx
        real x(nx), y(nx), a, b
        real sx,sy,sxx,sxy
        real x0,y0
        logical fixedpoint
        integer i
        sxx = 0.0
        sxy = 0.0
        sx  = 0.0
        sy  = 0.0
        do i = 1,nx
           sxy = sxy + x(i)*y(i)
           sxx = sxx + x(i)*x(i)
           sx = sx + x(i)
           sy = sy + y(i) 
        enddo

        if (fixedpoint) then
           a = (sxy-x0*sy-y0*sx+x0*y0*nx)/(sxx-2*x0*sx+x0*x0*nx)
           b = y0-a*x0
        else
           a = (sxy - sx*sy/nx) / (sxx - sx*sx/nx)
           b = (sy - a*sx)/nx
        endif

        return
        end


        real function expsum(xi,wi,p,a,n)
        implicit none 
        ! returns sum(xi^p exp(-a*xi)),
        !   xi weighted with wi 
        integer n,p,i 
        real xi(n),wi(n),a,res,x
        res=0.0 
        do i=1,n
           x=xi(i)*wi(i)
           res=res+(x**p)*exp(-a*x)
        enddo ! i 
        expsum=res
        return 
        end 
        

