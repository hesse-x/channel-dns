! ---&---1---------2---------3---------4---------5---------6---------7--
! FFTPACK transforms with persistent line and factor work arrays.
MODULE dns_fft_backend
  USE dns_kinds
  IMPLICIT NONE
  PRIVATE
  CHARACTER(*), PARAMETER, PUBLIC :: fft_backend_name = 'FFTPACK'
  PUBLIC :: fft_engine, init_fft_engine, close_fft_engine, backend_wisdom
  PUBLIC :: backend_x_forward, backend_x_backward, backend_z, backend_y
  TYPE :: fft_engine
    INTEGER :: nx, ny, nz, nyp, nxs, nxsp, nzs
    LOGICAL :: serial
    COMPLEX(dp), ALLOCATABLE :: xhat(:, :, :)
    COMPLEX(dp), POINTER, CONTIGUOUS :: k(:, :, :) => NULL(), y(:, :, :) => NULL()
    REAL(dp), ALLOCATABLE :: xl(:), zl(:), yl(:), tx(:), tz(:), ty(:)
    INTEGER :: ix(15), iz(15), iy(15)
  END TYPE
  INTERFACE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
    SUBROUTINE rffti(n, w, iw)
      IMPORT dp
      INTEGER :: n, iw(*)
      REAL(dp) :: w(*)
    END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
    SUBROUTINE cffti(n, w, iw)
      IMPORT dp
      INTEGER :: n, iw(*)
      REAL(dp) :: w(*)
    END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
    SUBROUTINE costi(n, w, iw)
      IMPORT dp
      INTEGER :: n, iw(*)
      REAL(dp) :: w(*)
    END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
    SUBROUTINE rfftf(n, x, w, iw)
      IMPORT dp
      INTEGER :: n, iw(*)
      REAL(dp) :: w(*), x(*)
    END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
    SUBROUTINE rfftb(n, x, w, iw)
      IMPORT dp
      INTEGER :: n, iw(*)
      REAL(dp) :: w(*), x(*)
    END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
    SUBROUTINE cfftf(n, x, w, iw)
      IMPORT dp
      INTEGER :: n, iw(*)
      REAL(dp) :: w(*), x(*)
    END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
    SUBROUTINE cfftb(n, x, w, iw)
      IMPORT dp
      INTEGER :: n, iw(*)
      REAL(dp) :: w(*), x(*)
    END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
    SUBROUTINE cost(n, x, w, iw)
      IMPORT dp
      INTEGER :: n, iw(*)
      REAL(dp) :: w(*), x(*)
    END SUBROUTINE
  END INTERFACE
CONTAINS

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE init_fft_engine(w, nxg, nyg, nzg, ypart, nxs, nxsp, nzs, serial)
    TYPE(fft_engine), INTENT(out) :: w
    INTEGER, INTENT(in) :: nxg, nyg, nzg, ypart, nxs, nxsp, nzs
    LOGICAL, INTENT(in) :: serial

    w%nx = nxg
    w%ny = nyg
    w%nz = nzg
    w%nyp = ypart
    w%nxs = nxs
    w%nxsp = nxsp
    w%nzs = nzs
    w%serial = serial
    ALLOCATE (w%xhat(0:nxg/2, 0:ypart, 0:nzg-1), w%k(0:nxs, 0:ypart, -nzs:nzs))
    IF (serial) THEN
      w%y => w%k
    ELSE
      ALLOCATE (w%y(0:nxsp, 0:nyg, -nzs:nzs))
    END IF
    ALLOCATE (w%xl(nxg), w%zl(2*nzg), w%yl(nyg+1), w%tx(2*nxg), w%tz(4*nzg), w%ty(3*(nyg+1)))
    CALL rffti(nxg, w%tx, w%ix)
    CALL cffti(nzg, w%tz, w%iz)
    CALL costi(nyg+1, w%ty, w%iy)
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE close_fft_engine(w)
    TYPE(fft_engine), INTENT(inout) :: w

    IF (.NOT. w%serial) DEALLOCATE (w%y)
    DEALLOCATE (w%k)
    NULLIFY (w%y, w%k)
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE backend_wisdom(path, exporting)
    CHARACTER(*), INTENT(in) :: path
    LOGICAL, INTENT(in) :: exporting
    ! FFTPACK initializes its deterministic factor tables without wisdom files.
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE backend_x_forward(w, physical)
    TYPE(fft_engine), INTENT(inout) :: w
    REAL(dp), INTENT(inout) :: physical(w%nx, 0:w%nyp, w%nz)
    INTEGER :: y, z, j

    DO z = 0, w%nz-1
      DO y = 0, w%nyp
        w%xl = physical(:, y, z+1)
        CALL rfftf(w%nx, w%xl, w%tx, w%ix)
        w%xhat(0, y, z) = CMPLX(w%xl(1), 0.0_DP, dp)
        DO j = 1, w%nx/2-1
          w%xhat(j, y, z) = CMPLX(w%xl(2*j), w%xl(2*j+1), dp)
        END DO
        w%xhat(w%nx/2, y, z) = CMPLX(w%xl(w%nx), 0.0_DP, dp)
      END DO
    END DO
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE backend_x_backward(w, physical)
    TYPE(fft_engine), INTENT(inout) :: w
    REAL(dp), INTENT(out) :: physical(w%nx, 0:w%nyp, w%nz)
    INTEGER :: y, z, j

    DO z = 0, w%nz-1
      DO y = 0, w%nyp
        w%xl(1) = REAL(w%xhat(0, y, z), dp)
        DO j = 1, w%nx/2-1
          w%xl(2*j) = REAL(w%xhat(j, y, z), dp)
          w%xl(2*j+1) = AIMAG(w%xhat(j, y, z))
        END DO
        w%xl(w%nx) = REAL(w%xhat(w%nx/2, y, z), dp)
        CALL rfftb(w%nx, w%xl, w%tx, w%ix)
        physical(:, y, z+1) = w%xl
      END DO
    END DO
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE backend_z(w, forward)
    TYPE(fft_engine), INTENT(inout) :: w
    LOGICAL, INTENT(in) :: forward
    INTEGER :: x, y, z

    DO y = 0, w%nyp
      DO x = 0, w%nxs
        DO z = 0, w%nz-1
          w%zl(2*z+1) = REAL(w%xhat(x, y, z), dp)
          w%zl(2*z+2) = AIMAG(w%xhat(x, y, z))
        END DO
        IF (forward) THEN
          CALL cfftf(w%nz, w%zl, w%tz, w%iz)
        ELSE
          CALL cfftb(w%nz, w%zl, w%tz, w%iz)
        END IF
        DO z = 0, w%nz-1
          w%xhat(x, y, z) = CMPLX(w%zl(2*z+1), w%zl(2*z+2), dp)
        END DO
      END DO
    END DO
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE backend_y(w)
    TYPE(fft_engine), INTENT(inout) :: w
    INTEGER :: x, z

    DO z = -w%nzs, w%nzs
      DO x = 0, w%nxsp
        w%yl = REAL(w%y(x, :, z), dp)
        CALL cost(w%ny+1, w%yl, w%ty, w%iy)
        w%y(x, :, z) = CMPLX(w%yl, AIMAG(w%y(x, :, z)), dp)
        w%yl = AIMAG(w%y(x, :, z))
        CALL cost(w%ny+1, w%yl, w%ty, w%iy)
        w%y(x, :, z) = CMPLX(REAL(w%y(x, :, z), dp), w%yl, dp)
      END DO
    END DO
  END SUBROUTINE
END MODULE
! ---&---1---------2---------3---------4---------5---------6---------7--
