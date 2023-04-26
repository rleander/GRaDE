module weibull
contains
function cdf(x, par) result (WeibullCDF)
   implicit none
   real               :: WeibullCDF
   real               :: x   ! argument
   real, dimension(:) :: par ! distribution parameters
   real :: a,k,c
   a = par(1)
   k = par(2)
   c = par(3)
   WeibullCDF = 1.-exp(-((max(x-c,0.0))/a)**k)
end function cdf

function xdf(x, par) result (Weibull_exceed)
   real               :: WeibullCDF
   real               :: x   ! argument
   real, dimension(:) :: par ! distribution parameters
   Weibull_exceed = 1.0 - cdf(x, par)
end function xdf

end module weibull

module integratie
contains
subroutine getpar(qr,qtbl,partbl,parint)
   implicit none
   real, intent(in) :: qr   !< reference discharge for which to find the parameters
   real, intent(in), dimension(:) :: qtbl
   real, intent(in), dimension(:,:) :: partbl
   real, intent(out), dimension(:) :: parint 
   
   integer :: i,j,n
   real :: w

   n = size(qtbl)
   do i=2,n-1
      if (qr<qtbl(i)) then
         exit
      endif
   enddo
   w = (qtbl(i)-qr)/(qtbl(i)-qtbl(i-1))
   parint = w*partbl(:,i-1) + (1.-w)*partbl(:,i)
end subroutine getpar

function integraal(ql,qref,F,par) result (som)
   use weibull
   implicit none
   real :: som, xdf0, xdf1
   real :: ql
   integer :: istep
   real, dimension(:)   :: qref   ! referentie q
   real, dimension(:)   :: F      ! onderschrijdingskans v(cumulatieve verdeling van qref)
   real, dimension(:,:) :: par    ! parameters van de versoring: onzekerheid
   integer :: nstep

   xdf1 = xdf(ql-qref(1),par(:,1))
   som = 0.5*xdf1*F(1)

   nstep = size(F)
   do istep = 2, size(F)
      xdf0 = xdf1
      xdf1 = xdf(ql-qref(istep),par(:,istep))
      som = som + 0.5*(xdf0+xdf1)*(F(istep)-F(istep-1))
   enddo
!  som = som + 0.5*(1.0-F(nstep))*(xdf1+1.)
end function integraal
end module integratie

program uitintegratie
use weibull
use integratie
implicit none

character(len=*), parameter :: fwerklijn = 'werklijn.txt'
! character(len=*), parameter :: fonzekerheid = 'weibull_parameters_uniform.txt'
! integer, parameter :: npar = 3 ! weibull 3-parameter

character(len=*), parameter :: fonzekerheid = 'beta_parameters_uniform.txt'
integer, parameter :: npar = 3 ! beta 4-parameter

character(len=666) :: regel
real :: sgv, exceedance
real, dimension(:), allocatable :: Tr
real, dimension(:), allocatable :: F
real, dimension(:), allocatable :: Qref
real, dimension(:), allocatable :: Qref_onz
real, dimension(:,:), allocatable :: parameters
real, dimension(:,:), allocatable :: parameters_int
integer :: i, j, i_werklijn, n_werklijn, n_onz, i_level
n_werklijn = 0
n_onz = 0

open(33,file=fwerklijn)
do while (.True.)
   read(33,'(a)',end=333) regel
   if (regel(1:1)/='#') n_werklijn = n_werklijn + 1
enddo
333 continue
allocate(Tr(n_werklijn), Qref(n_werklijn), F(n_werklijn))
i = 0
rewind(33)
do while (.True.)
   read(33,'(a)',end=303) regel
   if (regel(1:1)=='#') cycle
   i = i + 1
   read(regel,*) sgv, Tr(i), Qref(i)
   F(i) = 1. - 1./Tr(i)  ! non-exceedance probability
enddo
303 continue
close(33)

n_onz = 0
open(44,file=fonzekerheid)
do while (.True.)
   read(44,'(a)',end=444) regel
   if (regel(1:1)/='#') n_onz = n_onz + 1
enddo
444 continue
allocate(Qref_onz(n_onz), parameters(npar,n_onz))
rewind(44)
i = 0
do while (.True.)
   read(44,'(a)',end=404) regel
   if (regel(1:1)=='#') cycle
   i = i + 1
   read(regel,*) Qref_onz(i), parameters(:,i)
enddo
404 continue
close(44)

! Interpolatie van tabelwaarden naar de Q-waarden in de werklijn
allocate(parameters_int(npar,n_werklijn))
do i_werklijn = 1, n_werklijn
   call getpar(Qref(i_werklijn),Qref_onz,parameters,parameters_int(:,i_werklijn))
enddo

! uitintegratie voor gekozen set van levels
do i_level = 200, 20000, 200
   exceedance = integraal(real(i_level),Qref,F,parameters_int)
   print *,  i_level, 1./exceedance, -log(-log(1.-exceedance))
enddo
end program uitintegratie
