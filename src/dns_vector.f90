! ---&---1---------2---------3---------4---------5---------6---------7--
! Spectral vector operators and physical-space products.
MODULE dns_vector
  USE dns_kinds, ONLY: dp, fail
  USE mesh
  USE fft, ONLY: channel_xyzfft2, channel_xyzfft2_backward, CWAVX, CWAVZ
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: DOTTIME, CROSSTIME, Grad, Div, Curl, DiV_Grad, Div_BOUND, DISTANCE_L8, DERIV, FOR_Q
  PUBLIC :: prepare_vector_workspace, close_vector_workspace
  REAL(dp), ALLOCATABLE, SAVE :: product_a(:, :, :, :), product_b(:, :, :, :), PRODUCT(:, :, :, :)
  COMPLEX(dp), ALLOCATABLE, SAVE :: gradient_work(:, :, :, :)
CONTAINS

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE prepare_vector_workspace()
    IF (ALLOCATED(product_a)) RETURN
    ! Component-last physical buffers give FFTW contiguous planes for each component.
    ALLOCATE (product_a(NX2, 0:NY2P, NZ2, 3), product_b(NX2, 0:NY2P, NZ2, 3), PRODUCT(NX2, 0:NY2P, NZ2, 3))
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE close_vector_workspace()
    IF (ALLOCATED(product_a)) DEALLOCATE (product_a, product_b, product)
    IF (ALLOCATED(gradient_work)) DEALLOCATE (gradient_work)
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE DOTTIME(adb, a, b)
    COMPLEX(dp), INTENT(out) :: adb(0:NXSP, 0:NYS, -NZS:NZS)
    COMPLEX(dp), INTENT(in) :: a(3, 0:NXSP, 0:NYS, -NZS:NZS), b(3, 0:NXSP, 0:NYS, -NZS:NZS)
    INTEGER :: c, x, y, z

    CALL prepare_vector_workspace()
    DO c = 1, 3
      CALL channel_xyzfft2_backward(product_a(:, :, :, c), a(c, :, :, :))
      CALL channel_xyzfft2_backward(product_b(:, :, :, c), b(c, :, :, :))
    END DO
    DO z = 1, NZ2
      DO y = 0, NY2P
        DO x = 1, NX2
          PRODUCT(x, y, z, 1) = product_a(x, y, z, 1)*product_b(x, y, z, 1)+ &
                                product_a(x, y, z, 2)*product_b(x, y, z, 2)+product_a(x, y, z, 3)*product_b(x, y, z, 3)
        END DO
      END DO
    END DO
    CALL channel_xyzfft2('F', PRODUCT(:, :, :, 1), adb)
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE CROSSTIME(acb, a, b)
    COMPLEX(dp), INTENT(out) :: acb(3, 0:NXSP, 0:NYS, -NZS:NZS)
    COMPLEX(dp), INTENT(in) :: a(3, 0:NXSP, 0:NYS, -NZS:NZS), b(3, 0:NXSP, 0:NYS, -NZS:NZS)
    REAL(dp) :: a1, a2, a3, b1, b2, b3
    INTEGER :: c, x, y, z

    CALL prepare_vector_workspace()
    DO c = 1, 3
      CALL channel_xyzfft2_backward(product_a(:, :, :, c), a(c, :, :, :))
      CALL channel_xyzfft2_backward(product_b(:, :, :, c), b(c, :, :, :))
    END DO
    DO z = 1, NZ2
      DO y = 0, NY2P
        DO x = 1, NX2
          a1 = product_a(x, y, z, 1)
          a2 = product_a(x, y, z, 2)
          a3 = product_a(x, y, z, 3)
          b1 = product_b(x, y, z, 1)
          b2 = product_b(x, y, z, 2)
          b3 = product_b(x, y, z, 3)
          PRODUCT(x, y, z, 1) = a2*b3-a3*b2
          PRODUCT(x, y, z, 2) = a3*b1-a1*b3
          PRODUCT(x, y, z, 3) = a1*b2-a2*b1
        END DO
      END DO
    END DO
    DO c = 1, 3
      CALL channel_xyzfft2('F', PRODUCT(:, :, :, c), acb(c, :, :, :))
    END DO
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE DERIV(axis, out, input)
    INTEGER, INTENT(in) :: axis
    COMPLEX(dp), INTENT(out) :: out(0:, 0:, -NZS:)
    COMPLEX(dp), INTENT(in) :: input(0:, 0:, -NZS:)
    INTEGER :: x, y, z
    SELECT CASE (axis)
    CASE (1, 3)
      DO z = -NZS, NZS
        DO y = 0, NYS
          DO x = 0, NXSP
            IF (axis == 1) THEN
              out(x, y, z) = CWAVX(x)*input(x, y, z)
            ELSE
              out(x, y, z) = CWAVZ(z)*input(x, y, z)
            END IF
          END DO
        END DO
      END DO
    CASE (2)
      DO z = -NZS, NZS
        out(:, NYS, z) = 0
        IF (NYS == 0) CYCLE
        DO x = 0, NXSP
          out(x, NYS-1, z) = REAL(2*NYS, dp)*input(x, NYS, z)
        END DO
        DO y = NYS-2, 0, -1
          DO x = 0, NXSP
            out(x, y, z) = out(x, y+2, z)+REAL(2*(y+1), dp)*input(x, y+1, z)
          END DO
        END DO
        out(:, 0, z) = 0.5_DP*out(:, 0, z)
      END DO
    CASE default
      CALL fail('derivative axis must be 1, 2, or 3')
    END SELECT
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE Grad(out, input)
    COMPLEX(dp), INTENT(out) :: out(3, 0:NXSP, 0:NYS, -NZS:NZS)
    COMPLEX(dp), INTENT(in) :: input(0:NXSP, 0:NYS, -NZS:NZS)
    INTEGER :: x, y, z

    DO z = -NZS, NZS
      out(2, :, NYS, z) = 0
      IF (NYS > 0) THEN
        DO x = 0, NXSP
          out(2, x, NYS-1, z) = REAL(2*NYS, dp)*input(x, NYS, z)
        END DO
        DO y = NYS-2, 0, -1
          DO x = 0, NXSP
            out(2, x, y, z) = out(2, x, y+2, z)+REAL(2*(y+1), dp)*input(x, y+1, z)
          END DO
        END DO
        out(2, :, 0, z) = 0.5_DP*out(2, :, 0, z)
      END IF
      DO y = 0, NYS
        DO x = 0, NXSP
          out(1, x, y, z) = CWAVX(x)*input(x, y, z)
          out(3, x, y, z) = CWAVZ(z)*input(x, y, z)
        END DO
      END DO
    END DO
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE Div(out, input)
    COMPLEX(dp), INTENT(out) :: out(0:NXSP, 0:NYS, -NZS:NZS)
    COMPLEX(dp), INTENT(in) :: input(3, 0:NXSP, 0:NYS, -NZS:NZS)
    INTEGER :: x, y, z

    DO z = -NZS, NZS
      out(:, NYS, z) = 0
      IF (NYS > 0) THEN
        DO x = 0, NXSP
          out(x, NYS-1, z) = REAL(2*NYS, dp)*input(2, x, NYS, z)
        END DO
        DO y = NYS-2, 0, -1
          DO x = 0, NXSP
            out(x, y, z) = out(x, y+2, z)+REAL(2*(y+1), dp)*input(2, x, y+1, z)
          END DO
        END DO
        out(:, 0, z) = 0.5_DP*out(:, 0, z)
      END IF
      DO y = 0, NYS
        DO x = 0, NXSP
          out(x, y, z) = CWAVX(x)*input(1, x, y, z)+out(x, y, z)+CWAVZ(z)*input(3, x, y, z)
        END DO
      END DO
    END DO
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE Curl(out, input)
    COMPLEX(dp), INTENT(out) :: out(3, 0:NXSP, 0:NYS, -NZS:NZS)
    COMPLEX(dp), INTENT(in) :: input(3, 0:NXSP, 0:NYS, -NZS:NZS)
    INTEGER :: x, y, z

    DO z = -NZS, NZS
      out(1, :, NYS, z) = 0
      out(3, :, NYS, z) = 0
      IF (NYS > 0) THEN
        DO x = 0, NXSP
          out(1, x, NYS-1, z) = REAL(2*NYS, dp)*input(3, x, NYS, z)
          out(3, x, NYS-1, z) = REAL(2*NYS, dp)*input(1, x, NYS, z)
        END DO
        DO y = NYS-2, 0, -1
          DO x = 0, NXSP
            out(1, x, y, z) = out(1, x, y+2, z)+REAL(2*(y+1), dp)*input(3, x, y+1, z)
            out(3, x, y, z) = out(3, x, y+2, z)+REAL(2*(y+1), dp)*input(1, x, y+1, z)
          END DO
        END DO
        out(1, :, 0, z) = 0.5_DP*out(1, :, 0, z)
        out(3, :, 0, z) = 0.5_DP*out(3, :, 0, z)
      END IF
      DO y = 0, NYS
        DO x = 0, NXSP
          out(1, x, y, z) = out(1, x, y, z)-CWAVZ(z)*input(2, x, y, z)
          out(2, x, y, z) = CWAVZ(z)*input(1, x, y, z)-CWAVX(x)*input(3, x, y, z)
          out(3, x, y, z) = CWAVX(x)*input(2, x, y, z)-out(3, x, y, z)
        END DO
      END DO
    END DO
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE DiV_Grad(out, input)
    COMPLEX(dp), INTENT(out) :: out(0:NXSP, 0:NYS, -NZS:NZS)
    COMPLEX(dp), INTENT(in) :: input(0:NXSP, 0:NYS, -NZS:NZS)

    IF (.NOT. ALLOCATED(gradient_work)) ALLOCATE (gradient_work(3, 0:NXSP, 0:NYS, -NZS:NZS))
    CALL Grad(gradient_work, input)
    CALL Div(out, gradient_work)
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE Div_BOUND(out, v)
    COMPLEX(dp), INTENT(out) :: out(0:NXSP, 0:1, -NZS:NZS)
    COMPLEX(dp), INTENT(in) :: v(3, 0:NXSP, 0:NYS, -NZS:NZS)
    COMPLEX(dp) :: tangent, normal
    INTEGER :: x, y, z

    out = 0
    DO z = -NZS, NZS
      DO y = 0, NYS
        DO x = 0, NXSP
          tangent = CWAVX(x)*v(1, x, y, z)+CWAVZ(z)*v(3, x, y, z)
          normal = REAL(y*y, dp)*v(2, x, y, z)
          out(x, 0, z) = out(x, 0, z)+tangent+normal
          IF (MOD(y, 2) == 0) THEN
            out(x, 1, z) = out(x, 1, z)+tangent-normal
          ELSE
            out(x, 1, z) = out(x, 1, z)-tangent+normal
          END IF
        END DO
      END DO
    END DO
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE DISTANCE_L8(distance, u, v)
    REAL(dp), INTENT(out) :: distance
    COMPLEX(dp), INTENT(in) :: u(0:NXSP, 0:NYS, -NZS:NZS), v(0:NXSP, 0:NYS, -NZS:NZS)
    INTEGER :: x, y, z

    distance = 0
    DO z = -NZS, NZS
      DO y = 0, NYS
        DO x = 0, NXSP
          distance = MAX(distance, ABS(u(x, y, z)-v(x, y, z)))
        END DO
      END DO
    END DO
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE FOR_Q(q, v)
    REAL(dp), INTENT(out) :: q
    COMPLEX(dp), INTENT(in) :: v(0:, 0:, -NZS:)
    INTEGER :: y

    q = 0
    DO y = 0, NYS, 2
      q = q-REAL(v(0, y, 0), dp)/REAL(y*y-1, dp)
    END DO
  END SUBROUTINE
END MODULE
! ---&---1---------2---------3---------4---------5---------6---------7--
