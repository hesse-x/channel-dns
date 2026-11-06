! ---&---1---------2---------3---------4---------5---------6---------7--
PROGRAM test_numerics
  USE dns_kinds
  USE dns_config
  USE mesh
  USE fft
  USE dns_vector, ONLY: DERIV
  USE dns_equations, ONLY: WN_AD, SOLVE
  IMPLICIT NONE
  REAL(dp), ALLOCATABLE :: up(:, :, :), expected(:, :, :), out(:, :, :), wn(:, :), ad(:, :, :)
  COMPLEX(dp), ALLOCATABLE :: u(:, :, :), du(:, :, :), rhs(:, :, :), bc(:, :, :), guarded(:, :, :, :)
  REAL(dp), ALLOCATABLE, TARGET :: shifted_buffer(:)
  REAL(dp), POINTER, CONTIGUOUS :: shifted(:, :, :)
  INTEGER :: offset, nphysical
  REAL(dp) :: pi, angle, yy, xx, zz, global, w
  INTEGER :: x, y, z, err, ranks
  CHARACTER(4096) :: directory

  CALL MPI_Init(err)
  CALL MPI_Comm_size(MPI_COMM_WORLD, ranks, err)
  NX = 16
  NY = 16
  NZ = 16
  NXS = 8
  NYS = 12
  NZS = 8
  NX2 = 24
  NY2 = 24
  NZ2 = 24
  INDEX_Serial = 0
  IF (ranks == 1) INDEX_Serial = 1
  XL = 10
  ZL = 7
  Re = 2800
  CALL GET_COMMAND_ARGUMENT(1, directory)
  CALL setup_mesh(XL, ZL)
  CALL setup_fft(XL, ZL, Re, TRIM(directory))
  ALLOCATE (up(NX, 0:NYP, NZ), out(NX, 0:NYP, NZ), expected(NX, 0:NYP, NZ))
  ALLOCATE (u(0:NXSP, 0:NYS, -NZS:NZS), du(0:NXSP, 0:NYS, -NZS:NZS), rhs(0:NXSP, 0:NYS, -NZS:NZS))
  ALLOCATE (bc(0:NXSP, 0:1, -NZS:NZS), wn(0:NXSP, -NZS:NZS), ad(0:NXSP, 2:NYS, -NZS:NZS))
  pi = ACOS(-1.0_DP)
  DO z = 1, NZ
    zz = 2*pi*(z-1)/NZ
    DO y = 0, NYP
      yy = YC(y)
      DO x = 1, NX
        xx = 2*pi*(x-1)/NX
        up(x, y, z) = SIN(xx)*COS(zz)*(1-yy*yy)+0.3_DP*yy*yy
        expected(x, y, z) = -2*yy*SIN(xx)*COS(zz)+0.6_DP*yy
      END DO
    END DO
  END DO
  CALL channel_xyzfft('F', up, u)
  CALL channel_xyzfft('B', out, u)
  CALL check(MAXVAL(ABS(out-up)), 2.E-13_DP, 'coarse FFT round trip')
  ! Verify direct new-array execution and the safe fallback across alignment classes.
  nphysical = NX*(NYP+1)*NZ
  ALLOCATE (shifted_buffer(nphysical+3), guarded(3, 0:NXSP, 0:NYS, -NZS:NZS))
  DO offset = 1, 4
    shifted(1:NX, 0:NYP, 1:NZ) => shifted_buffer(offset:offset+nphysical-1)
    guarded = CMPLX(777.0_DP, 11.0_DP, dp)
    guarded(2, :, :, :) = u
    CALL channel_xyzfft('B', shifted, guarded(2, :, :, :))
    CALL check(MAXVAL(ABS(shifted-up)), 2.E-13_DP, 'strided component backward transform / alignment')
    CALL channel_xyzfft('F', shifted, guarded(2, :, :, :))
    CALL check(MAXVAL(ABS(guarded(2, :, :, :)-u)), 2.E-13_DP, 'strided component forward transform / alignment')
    CALL check(MAXVAL(ABS(guarded(1, :, :, :)-CMPLX(777.0_DP, 11.0_DP, dp))), 0.0_DP, 'first component guard')
    CALL check(MAXVAL(ABS(guarded(3, :, :, :)-CMPLX(777.0_DP, 11.0_DP, dp))), 0.0_DP, 'last component guard')
  END DO
  NULLIFY (shifted)
  DEALLOCATE (shifted_buffer, guarded)

  CALL DERIV(2, du, u)
  CALL channel_xyzfft('B', out, du)
  CALL check(MAXVAL(ABS(out-expected)), 2.E-12_DP, 'Chebyshev derivative')
  CALL DERIV(1, du, u)
  CALL channel_xyzfft('B', out, du)
  DO z = 1, NZ
    zz = 2*pi*(z-1)/NZ
    DO y = 0, NYP
      DO x = 1, NX
        expected(x, y, z) = 2*pi/XL*COS(2*pi*(x-1)/NX)*COS(zz)*(1-YC(y)**2)
      END DO
    END DO
  END DO
  CALL check(MAXVAL(ABS(out-expected)), 2.E-12_DP, 'Fourier derivative')
  ! Full-spectrum x/z Nyquist modes must split and merge on the coarse grid.
  DO z = 1, NZ
    DO y = 0, NYP
      DO x = 1, NX
        up(x, y, z) = REAL((-1)**(x+z), dp)*(1-YC(y)**2)
      END DO
    END DO
  END DO
  CALL channel_xyzfft('F', up, u)
  CALL channel_xyzfft('B', out, u)
  CALL check(MAXVAL(ABS(out-up)), 2.E-13_DP, 'Nyquist round trip')
  ! Fine-grid reconstruction of a truncated analytic field.
  DEALLOCATE (up, out, expected)
  ALLOCATE (up(NX2, 0:NY2P, NZ2), out(NX2, 0:NY2P, NZ2))
  DO z = 1, NZ2
    DO y = 0, NY2P
      angle = pi*(LSUBY2(my_id)+y)/NY2
      DO x = 1, NX2
        up(x, y, z) = SIN(2*pi*(x-1)/NX2)*COS(2*pi*(z-1)/NZ2)*(1-COS(angle)**2)
      END DO
    END DO
  END DO
  CALL channel_xyzfft2('F', up, u)
  CALL channel_xyzfft2('B', out, u)
  CALL check(MAXVAL(ABS(out-up)), 3.E-13_DP, 'dealiased FFT round trip')
  ! Manufactured Helmholtz problem: u=1-y^2, homogeneous Dirichlet walls.
  CALL WN_AD(wn, ad, 3.0_DP)
  rhs = 0
  bc = 0
  DO z = -NZS, NZS
    DO x = 0, NXSP
      w = -wn(x, z)
      rhs(x, 0, z) = -2-0.5_DP*w
      rhs(x, 2, z) = 0.5_DP*w
    END DO
  END DO
  CALL SOLVE(u, rhs, 3.0_DP, 1.0_DP, 0.0_DP, bc, wn, ad)
  du = 0
  du(:, 0, :) = 0.5_DP
  du(:, 2, :) = -0.5_DP
  CALL check(MAXVAL(ABS(u-du)), 1.E-12_DP, 'manufactured Helmholtz solve')
  CALL close_fft()
  IF (my_id == 0) PRINT *, 'PASS: Fourier/Chebyshev, Nyquist, transpose, Helmholtz'
  CALL MPI_Finalize(err)
CONTAINS

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE check(VALUE, tolerance, name)
    REAL(dp), INTENT(in) :: VALUE, tolerance
    CHARACTER(*), INTENT(in) :: name

    CALL MPI_Allreduce(VALUE, global, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, err)
    CALL mpi_check(err, 'test error reduction')
    IF (global > tolerance) THEN
      IF (my_id == 0) WRITE (*, *) name, global
      CALL fail(name)
    END IF
  END SUBROUTINE
END PROGRAM
! ---&---1---------2---------3---------4---------5---------6---------7--
