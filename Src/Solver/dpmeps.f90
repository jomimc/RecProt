
    double precision function dpmeps()
!     **********

!     Subroutine dpeps

!     This subroutine computes the machine precision parameter
!     dpmeps as the smallest floating point number such that
!     1 + dpmeps differs from 1.

!     This subroutine is based on the subroutine machar described in

!     W. J. Cody,
!     MACHAR: A subroutine to dynamically determine machine parameters,
!     ACM Trans. Math. Soft., 14, 1988, pages 303-311.

!     The subroutine statement is:

!       subroutine dpeps(dpmeps)

!     where

!       dpmeps is a double precision variable.
!         On entry dpmeps need not be specified.
!         On exit dpmeps is the machine precision.

!     MINPACK-2 Project. February 1991.
!     Argonne National Laboratory and University of Minnesota.
!     Brett M. Averick.

!     *******
    dpmeps=2.2D-16
    return
  end function dpmeps
  
