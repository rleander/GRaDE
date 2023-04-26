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
        real Tthreshold			! Threshold return priod from which to apply weissman

        real Tr, Tmax                   ! return period
        real yy, dy, ymax
        integer nstep


        ! CENSORING 
        real sigma_hat,mean_xk
        real mu_hat, F
        integer k 

        character*(6) paramstr 
        parameter(euler=0.577)	
        
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
            k = -1
            read(paramstr,*,err=506) k 
 506        if (k<0) then
               read(paramstr,*) Tthreshold
            endif
            call getarg(3,paramstr)
            read(paramstr,*) Tmax
            call getarg(4,paramstr)
            read(paramstr,*) nstep
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
        if (k<0) then
            k = nx - floor(0.3 + (nx+0.4)*(1.-1./Tthreshold))
            write(0,'(a,f6.0,a,i0)') 'T = ',Tthreshold,', k = ', k 
        endif
        k=max(2,min(k,nx))

        ! sort values 
        call ssort(xi,nx)

        ! fill array with gumbel variates
        do ix=1,nx 
           yi(ix)=-log(-log((ix-0.3)/(nx+0.4)))
        enddo ! ix 

        ! ESTIMATE 
        mean_xk=0.0			! average the k highest values 
        do ix=nx-k+1,nx 
           mean_xk=mean_xk+xi(ix)
        enddo 
        mean_xk=mean_xk/k

        sigma_hat=mean_xk-xi(nx-k+1)	! Weissman formulas
        mu_hat=xi(nx-k+1)+sigma_hat*log(1.*k)


        ! write results to stdout 
        ! gumbel var, x_i, fitted x_i, alpha estimate, u estimate

        write(*,*) '# Censored gumbel fit, using Weissman formula' 
        write(*,*) '# upper', k,' values of',nx 
        write(*,*) '# 1 gumbel variate'
        write(*,*) '# 2 quantile'
        write(*,*) '# 3 gumbel fitted quantile'
        !write(*,*) '# 4 stddev of fitted quantile'
        
        write(0,'(a,E12.4,a,E12.4,a,E12.4,a)') 
     &   '# logfit(x,',sigma_hat,',',xi(nx-k+1),',',1.*k/nx,')'
!	do ix=1,nx
!	   write(*,'(10f10.3)') 
!     &	         yi(ix),xi(ix)
!     &	        ,xi(nx-k+1)+sigma_hat
!     &	             *log(k/(nx*(1.-((ix-0.3)/(nx+0.4)) )))
!	enddo !ix 
        do ix=1,nx-k
           write(*,'(10f10.3)') 
     &         yi(ix),xi(ix),xi(ix)
        enddo !ix 
        do ix=nx-k+1,nx               ! replace the upper part by the fit for the second column
           write(*,'(10f10.3)') 
     &         yi(ix)
     &        ,xi(nx-k+1)+sigma_hat
     &             *log(k/(nx*(1.-((ix-0.3)/(nx+0.4)) )))
     &        ,xi(ix)
        enddo !ix 

        ! extrapolation
        ymax = -log(-log(1.-1./Tmax))          ! max std gum var
        dy = (ymax - yi(nx))/nstep
        do ix = 1, nstep
           yy = yi(nx) + dy*ix
           F = exp(-exp(-yy))
           write(*,'(10f10.3)') yy, 
     &                xi(nx-k+1)+sigma_hat*log(k/(nx*(1.-F))) ! linear extrapolation
        enddo


        do iarg=3,iargc()
           call getarg(iarg,paramstr)
           read(paramstr,*,err=555,end=555) Tr		! return period 		
           write(0,'(f8.1,E14.6)') Tr,xi(nx-k+1)
     &        +sigma_hat*log(k*Tr/nx)
 555       continue  
        enddo ! iarg  
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
        

