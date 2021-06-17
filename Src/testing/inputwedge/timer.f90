
    subroutine timer(ttime)
    double precision :: ttime
!     *********

!     Subroutine timer

!     This subroutine is used to determine user time. In a typical
!     application, the user time for a code segment requires calls
!     to subroutine timer to determine the initial and final time.

!     The subroutine statement is

!       subroutine timer(ttime)

!     where

!       ttime is an output variable which specifies the user time.

!     Argonne National Laboratory and University of Minnesota.
!     MINPACK-2 Project.

!     Modified October 1990 by Brett M. Averick.

!     **********
    real :: result
    real :: tarray(2)
!      real etime

!     The first element of the array tarray specifies user time

    call etime(tarray,result)

    ttime = dble(tarray(1))
    return

    end subroutine timer
          
!====================== The end of timer ===============================

! rogram test_etime
!    integer(8) :: i, j
!    real, dimension(2) :: tarray
!    real :: result
!    call ETIME(tarray, result)
!    print *, result
!    print *, tarray(1)
!    print *, tarray(2)
!    do i=1,100000000    ! Just a delay
!        j = i * i - i
!    end do
!    call ETIME(tarray, result)
!    print *, result
!    print *, tarray(1)
!    print *, tarray(2)
! nd program test_etime
