module point_in_poly
contains
   subroutine intersect(x1,y1,x2,y2,x3,y3,x4,y4,k,l)
      ! intersect (x1,y1)->(x2,y2)  with  (x3,y3)->(x4,y4)
      implicit none
      real(kind=8), intent(in) :: x1,x2,x3,x4
      real(kind=8), intent(in) :: y1,y2,y3,y4
      real(kind=8), intent(out) :: k,l
      l = ((x3-x1)*(y4-y3)-(y3-y1)*(x4-x3)) / ((x2-x1)*(y4-y3)-(y2-y1)*(x4-x3))
      k = ((x1-x3)*(y4-y3)+l*(x2-x1)*(y4-y3))/(x4-x3)*(y4-y3)
   end subroutine intersect

   function pinpok2D(x,y,xc,yc) result (s)
      implicit none
      logical                                  :: s
      real(kind=8), dimension(:), intent(in)   :: x, y
      real(kind=8),               intent(in)   :: xc,yc
      real(kind=8) :: k, l, xo, yo
      real(kind=8), dimension(2) :: xbox, ybox
      integer :: nspt = 0
      integer :: n, i
      xbox=(/minval(x),maxval(x)/)
      ybox=(/minval(y),maxval(y)/)
      xo = 2*xbox(1)-xbox(2)
      yo = 2*ybox(1)-ybox(2)
      n = min(size(x),size(y))
      call intersect(x(n),y(n),x(1),y(1),xc,yc,xo,yo,k,l)
      if (k*(k-1.d0)<0.d0 .and. l*(l-1.d0)<0.d0) nspt = nspt + 1
      do i=1,n-1
         call intersect(x(i),y(i),x(i+1),y(i+1),xc,yc,xo,yo,k,l)
         if (k*(k-1.d0)<0.d0 .and. l*(l-1.d0)<0.d0) nspt = nspt + 1
      enddo
      s = (mod(nspt,2) == 1)            ! odd number of intersections
      return
   end function pinpok2D
end module point_in_poly
