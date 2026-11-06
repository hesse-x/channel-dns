! ---&---1---------2---------3---------4---------5---------6---------7--
! Shared transforms, normalization and MPI transpose for both FFT backends.
MODULE fft
  USE dns_kinds
  USE mesh
  USE dns_fft_backend
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: fft_backend_name
  PUBLIC :: setup_fft, close_fft, channel_xyzfft, channel_xyzfft2, xzfft
  PUBLIC :: channel_xyzfft_backward, channel_xyzfft2_backward
  PUBLIC :: WAVX, WAVZ, CWAVX, CWAVZ, CWAVXRe, CWAVZRe, WL, WD, WR
  REAL(dp), ALLOCATABLE :: WAVX(:), WAVZ(:), WL(:), WD(:), WR(:)
  COMPLEX(dp), ALLOCATABLE :: CWAVX(:), CWAVZ(:), CWAVXRe(:), CWAVZRe(:)
  TYPE :: workspace
    INTEGER :: nx, ny, nz, nyp
    INTEGER, ALLOCATABLE :: ybegin(:), yupper(:), sc(:), sd(:), rc(:), rd(:)
    COMPLEX(dp), ALLOCATABLE :: send(:), recv(:)
    TYPE(fft_engine) :: engine
  END TYPE
  TYPE(workspace), SAVE :: coarse, fine
  CHARACTER(:), ALLOCATABLE :: wisdom_path
CONTAINS

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE init_workspace(w, nxg, nyg, nzg, ypart, ybeg, yupper)
    TYPE(workspace), INTENT(out) :: w
    INTEGER, INTENT(in) :: nxg, nyg, nzg, ypart, ybeg(0:), yupper(0:)
    INTEGER :: r

    w%nx = nxg
    w%ny = nyg
    w%nz = nzg
    w%nyp = ypart
    ALLOCATE (w%ybegin(0:NPART-1), w%yupper(0:NPART-1))
    w%ybegin = ybeg
    w%yupper = yupper
    ! Keep the rank-indexed tables explicitly zero based.
    ALLOCATE (w%sc(0:NPART-1), w%sd(0:NPART-1), w%rc(0:NPART-1), w%rd(0:NPART-1))
    CALL init_fft_engine(w%engine, nxg, nyg, nzg, ypart, NXS, NXSP, NZS, NPART == 1)
    DO r = 0, NPART-1
      w%sc(r) = (NXSPT(r)+1)*(ypart+1)*NZS2P1
      w%rc(r) = (NXSP+1)*(yupper(r)+1)*NZS2P1
    END DO
    w%sd(0) = 0
    w%rd(0) = 0
    DO r = 1, NPART-1
      w%sd(r) = w%sd(r-1)+w%sc(r-1)
      w%rd(r) = w%rd(r-1)+w%rc(r-1)
    END DO
    ALLOCATE (w%send(MAX(SUM(w%sc), SUM(w%rc))), w%recv(MAX(SUM(w%sc), SUM(w%rc))))
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE setup_fft(xl, zl, re, directory)
    REAL(dp), INTENT(in) :: xl, zl, re
    CHARACTER(*), INTENT(in) :: directory
    INTEGER :: i
    REAL(dp) :: pi
    CHARACTER(80) :: suffix
    ALLOCATE (WAVX(0:NXSP), WAVZ(-NZS:NZS), CWAVX(0:NXSP), CWAVZ(-NZS:NZS))
    ALLOCATE (CWAVXRe(0:NXSP), CWAVZRe(-NZS:NZS), WL(0:NYS), WD(0:NYS), WR(0:NYS))
    pi = ACOS(-1.0_DP)
    DO i = 0, NXSP
      WAVX(i) = (LSUBX(my_id)+i)*(2*pi/xl)
    END DO
    DO i = -NZS, NZS
      WAVZ(i) = i*(2*pi/zl)
    END DO
    CWAVX = CMPLX(0.0_DP, WAVX, dp)
    CWAVZ = CMPLX(0.0_DP, WAVZ, dp)
    CWAVXRe = CWAVX/re
    CWAVZRe = CWAVZ/re
    WL = 0
    WD = 0
    WR = 0
    DO i = 2, NYS
      WL(i) = 1.0_DP/(4.0_DP*i*(i-1))
      WD(i) = -1.0_DP/(2.0_DP*(i*i-1))
      WR(i) = 1.0_DP/(4.0_DP*i*(i+1))
    END DO
    WL(2) = 2*WL(2)
    WD(NYS-1:NYS) = 0
    WR(NYS-3:NYS) = 0
    suffix = ''
    IF (NPROC > 1) WRITE (suffix, '(a,i0,a,i0)') '.np', NPROC, '.rank', my_id
    wisdom_path = directory//'/channel_dns.wisdom'//TRIM(suffix)
    CALL backend_wisdom(wisdom_path, .FALSE.)
    CALL init_workspace(coarse, NX, NY, NZ, NYP, LSUBY, NYPT)
    CALL init_workspace(fine, NX2, NY2, NZ2, NY2P, LSUBY2, NY2PT)
    CALL backend_wisdom(wisdom_path, .TRUE.)
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE destroy_workspace(w)
    TYPE(workspace), INTENT(inout) :: w

    CALL close_fft_engine(w%engine)
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE close_fft()
    CALL destroy_workspace(coarse)
    CALL destroy_workspace(fine)
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE TRANSPOSE(w, forward)
    TYPE(workspace), INTENT(inout) :: w
    LOGICAL, INTENT(in) :: forward
    INTEGER :: r, x, y, z, p, err, y0

    IF (NPART == 1) RETURN
    DO r = 0, NPART-1
      IF (forward) THEN
        p = w%sd(r)
        DO z = -NZS, NZS
          DO y = 0, w%nyp
            DO x = LSUBX(r), LSUBX(r)+NXSPT(r)
              p = p+1
              w%send(p) = w%engine%k(x, y, z)
            END DO
          END DO
        END DO
      ELSE
        p = w%rd(r)
        y0 = w%ybegin(r)
        DO z = -NZS, NZS
          DO y = y0, y0+w%yupper(r)
            DO x = 0, NXSP
              p = p+1
              w%send(p) = w%engine%y(x, y, z)
            END DO
          END DO
        END DO
      END IF
    END DO
    IF (forward) THEN
      CALL MPI_Alltoallv(w%send, w%sc, w%sd, MPI_DOUBLE_COMPLEX, w%recv, w%rc, w%rd, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, err)
    ELSE
      CALL MPI_Alltoallv(w%send, w%rc, w%rd, MPI_DOUBLE_COMPLEX, w%recv, w%sc, w%sd, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, err)
    END IF
    CALL mpi_check(err, 'FFT transpose')
    DO r = 0, NPART-1
      IF (forward) THEN
        p = w%rd(r)
        y0 = w%ybegin(r)
        DO z = -NZS, NZS
          DO y = y0, y0+w%yupper(r)
            DO x = 0, NXSP
              p = p+1
              w%engine%y(x, y, z) = w%recv(p)
            END DO
          END DO
        END DO
      ELSE
        p = w%sd(r)
        DO z = -NZS, NZS
          DO y = 0, w%nyp
            DO x = LSUBX(r), LSUBX(r)+NXSPT(r)
              p = p+1
              w%engine%k(x, y, z) = w%recv(p)
            END DO
          END DO
        END DO
      END IF
    END DO
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE forward_transform(w, physical, spectral, is_coarse)
    TYPE(workspace), INTENT(inout) :: w
    REAL(dp), INTENT(inout) :: physical(w%nx, 0:w%nyp, w%nz)
    COMPLEX(dp), INTENT(out) :: spectral(0:, 0:, -NZS:)
    LOGICAL, INTENT(in) :: is_coarse
    INTEGER :: x, y, z, slot
    REAL(dp) :: inv_nx, inv_ny, inv_nz
    LOGICAL :: fullx, fullz

    inv_nx = 1.0_DP/REAL(w%nx, dp)
    inv_ny = 1.0_DP/REAL(w%ny, dp)
    inv_nz = 1.0_DP/REAL(w%nz, dp)
    fullx = is_coarse .AND. NXS == w%nx/2
    fullz = is_coarse .AND. NZS == w%nz/2
    CALL backend_x_forward(w%engine, physical)
    w%engine%xhat(0:NXS, :, :) = w%engine%xhat(0:NXS, :, :)*inv_nx
    IF (fullx) w%engine%xhat(NXS, :, :) = 0.5_DP*w%engine%xhat(NXS, :, :)
    CALL backend_z(w%engine, .TRUE.)
    DO z = -NZS, NZS
      slot = MODULO(z, w%nz)
      w%engine%k(:, :, z) = w%engine%xhat(0:NXS, :, slot)*inv_nz
    END DO
    IF (fullz) THEN
      w%engine%k(:, :, NZS) = 0.5_DP*w%engine%k(:, :, NZS)
      w%engine%k(:, :, -NZS) = w%engine%k(:, :, NZS)
    END IF
    CALL TRANSPOSE(w, .TRUE.)
    CALL backend_y(w%engine)
    DO z = -NZS, NZS
      DO y = 0, NYS
        DO x = 0, NXSP
          spectral(x, y, z) = w%engine%y(x, y, z)*inv_ny
        END DO
      END DO
      spectral(:, 0, z) = 0.5_DP*spectral(:, 0, z)
    END DO
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE backward_transform(w, physical, spectral, is_coarse)
    TYPE(workspace), INTENT(inout) :: w
    REAL(dp), INTENT(out) :: physical(w%nx, 0:w%nyp, w%nz)
    COMPLEX(dp), INTENT(in) :: spectral(0:, 0:, -NZS:)
    LOGICAL, INTENT(in) :: is_coarse
    INTEGER :: x, y, z
    LOGICAL :: fullx, fullz

    fullx = is_coarse .AND. NXS == w%nx/2
    fullz = is_coarse .AND. NZS == w%nz/2
    w%engine%y(:, NYS+1:w%ny, :) = 0
    DO z = -NZS, NZS
      DO y = 0, NYS
        DO x = 0, NXSP
          w%engine%y(x, y, z) = spectral(x, y, z)
        END DO
      END DO
      w%engine%y(:, 0, z) = 2*w%engine%y(:, 0, z)
    END DO
    CALL backend_y(w%engine)
    CALL TRANSPOSE(w, .FALSE.)
    ! Clear only slots not written by expansion; do not write retained modes twice.
    w%engine%xhat(NXS+1:w%nx/2, :, :) = 0
    w%engine%xhat(0:NXS, :, NZS+1:w%nz-NZS-1) = 0
    DO z = 0, NZS
      w%engine%xhat(0:NXS, :, z) = 0.5_DP*w%engine%k(:, :, z)
    END DO
    DO z = -NZS, -1
      w%engine%xhat(0:NXS, :, w%nz+z) = 0.5_DP*w%engine%k(:, :, z)
    END DO
    IF (fullz) w%engine%xhat(0:NXS, :, NZS) = 0.5_DP*(w%engine%k(:, :, NZS)+w%engine%k(:, :, -NZS))
    CALL backend_z(w%engine, .FALSE.)
    IF (fullx) w%engine%xhat(NXS, :, :) = 2*w%engine%xhat(NXS, :, :)
    CALL backend_x_backward(w%engine, physical)
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE channel_xyzfft(c, up, u)
    CHARACTER, INTENT(in) :: c
    REAL(dp), INTENT(inout) :: up(NX, 0:NYP, NZ)
    COMPLEX(dp), INTENT(inout) :: u(0:, 0:, -NZS:)

    IF (c == 'F') THEN
      CALL forward_transform(coarse, up, u, .TRUE.)
    ELSE IF (c == 'B') THEN
      CALL backward_transform(coarse, up, u, .TRUE.)
    ELSE
      CALL fail('FFT direction must be F or B')
    END IF
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE channel_xyzfft2(c, up, u)
    CHARACTER, INTENT(in) :: c
    REAL(dp), INTENT(inout) :: up(NX2, 0:NY2P, NZ2)
    COMPLEX(dp), INTENT(inout) :: u(0:, 0:, -NZS:)

    IF (c == 'F') THEN
      CALL forward_transform(fine, up, u, .FALSE.)
    ELSE IF (c == 'B') THEN
      CALL backward_transform(fine, up, u, .FALSE.)
    ELSE
      CALL fail('FFT direction must be F or B')
    END IF
  END SUBROUTINE
  ! Read-only backward entries also accept strided vector components.

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE channel_xyzfft_backward(up, u)
    REAL(dp), INTENT(out) :: up(NX, 0:NYP, NZ)
    COMPLEX(dp), INTENT(in) :: u(0:, 0:, -NZS:)

    CALL backward_transform(coarse, up, u, .TRUE.)
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE channel_xyzfft2_backward(up, u)
    REAL(dp), INTENT(out) :: up(NX2, 0:NY2P, NZ2)
    COMPLEX(dp), INTENT(in) :: u(0:, 0:, -NZS:)

    CALL backward_transform(fine, up, u, .FALSE.)
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE xzfft(up, u)
    REAL(dp), INTENT(in) :: up(NX, NZ)
    COMPLEX(dp), INTENT(out) :: u(0:NXSP, -NZS:NZS)
    REAL(dp), ALLOCATABLE :: physical(:, :, :)
    COMPLEX(dp), ALLOCATABLE :: spectral(:, :, :)
    INTEGER :: y
    ALLOCATE (physical(NX, 0:NYP, NZ), spectral(0:NXSP, 0:NYS, -NZS:NZS))
    DO y = 0, NYP
      physical(:, y, :) = up
    END DO
    CALL channel_xyzfft('F', physical, spectral)
    u = spectral(:, 0, :)
  END SUBROUTINE
END MODULE
! ---&---1---------2---------3---------4---------5---------6---------7--
