module ddotmod
  contains
      double precision function ddot(n,dx,incx,dy,incy)

!     forms the dot product of two vectors.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.

    double precision :: dx(*),dy(*),dtemp
    integer :: i,incx,incy,ix,iy,n

    ddot = 0.0d0
    dtemp = 0.0d0
    if(n <= 0)return
    if(incx == 1 .AND. incy == 1)go to 20

!        code for unequal increments or equal increments
!          not equal to 1

    ix = 1
    iy = 1
    if(incx < 0)ix = (-n+1)*incx + 1
    if(incy < 0)iy = (-n+1)*incy + 1
    do 10 i = 1,n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
    10 END DO
    ddot = dtemp
    return

!        code for both increments equal to 1


!        clean-up loop

    20 continue
    do i=1,n
        dtemp=dtemp+dx(i)*dy(i)
    enddo
    60 ddot = dtemp
    return
    end function ddot
      double precision function ddot1(n,dx,dy)

!     forms the dot product of two vectors.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.

    double precision :: dx(*),dy(*),dtemp
    integer :: n,i

    ddot1 = 0.0d0
    dtemp = 0.0d0
    if(n <= 0)return
    do i=1,n
        dtemp=dtemp+dx(i)*dy(i)
    enddo
    ddot1 = dtemp
    return
    end function ddot1
  end module ddotmod
  
  
