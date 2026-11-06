! ---&---1---------2---------3---------4---------5---------6---------7--
! Wall boundary values, initial velocity and perturbations.
MODULE dns_initial
  USE dns_kinds
  USE dns_config
  USE mesh
  USE fft
  IMPLICIT NONE
CONTAINS

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE velocity_boundary(bc)
    COMPLEX(dp), INTENT(out) :: bc(3, 0:NXSP, 0:1, -NZS:NZS)

    bc = 0
    IF (my_id == 0) bc(2, 0, :, 0) = V_Microblow0
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE initialize_velocity(v)
    COMPLEX(dp), INTENT(out) :: v(3, 0:NXSP, 0:NYS, -NZS:NZS)
    REAL(dp), ALLOCATABLE :: physical(:, :, :)
    REAL(dp) :: uc, SCALE(3)
    INTEGER :: c, y
    ALLOCATE (physical(NX, 0:NYP, NZ))
    uc = 1.5_DP
    IF (INDEX_GRAD_PRE == 1) uc = Re/2
    scale = [uc*COS(theta_rad), V_Microblow0, uc*SIN(theta_rad)]
    DO c = 1, 3
      DO y = 0, NYP
        physical(:, y, :) = SCALE(c)
        IF (c /= 2) physical(:, y, :) = SCALE(c)*(1-YC(y)**2)
      END DO
      CALL channel_xyzfft('F', physical, v(c, :, :, :))
    END DO
    CALL disturb_velocity(v, EPS_init)
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE disturb_velocity(v, amplitude)
    COMPLEX(dp), INTENT(inout) :: v(3, 0:NXSP, 0:NYS, -NZS:NZS)
    REAL(dp), INTENT(in) :: amplitude
    REAL(dp), ALLOCATABLE :: physical(:, :, :, :), noise(:, :, :, :)
    REAL(dp) :: uc
    INTEGER :: c, y

    IF (amplitude == 0) RETURN
    ALLOCATE (physical(NX, 0:NYP, NZ, 3), noise(NX, 0:NYP, NZ, 3))
    uc = 1.5_DP
    IF (INDEX_GRAD_PRE == 1) uc = Re/2
    CALL RANDOM_NUMBER(noise)
    noise = noise-0.5_DP
    DO c = 1, 3
      CALL channel_xyzfft('B', physical(:, :, :, c), v(c, :, :, :))
      DO y = 0, NYP
        IF (c == 3) noise(:, y, :, c) = noise(:, y, :, c)*INDEX_disturb_Z
        physical(:, y, :, c) = physical(:, y, :, c)+uc*(1-YC(y)**2)*amplitude*noise(:, y, :, c)
      END DO
      CALL channel_xyzfft('F', physical(:, :, :, c), v(c, :, :, :))
    END DO
  END SUBROUTINE
END MODULE
! ---&---1---------2---------3---------4---------5---------6---------7--
