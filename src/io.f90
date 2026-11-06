! ---&---1---------2---------3---------4---------5---------6---------7--
! Spectral checkpoints and physical fields in text and binary form.
MODULE dns_io
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_CHAR, C_NULL_CHAR
  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: INT32
  USE, INTRINSIC :: IEEE_ARITHMETIC
  USE dns_kinds
  USE dns_config, ONLY: output_format
  USE mesh
  USE fft
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: rank_file, write_restart, read_restart, snapshot, write_physical, atomic_finish, make_directory
  CHARACTER(8), PARAMETER :: magic = 'CHDNS'//ACHAR(0)//ACHAR(0)//ACHAR(0)
  INTERFACE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
    FUNCTION c_rename(old, new) BIND(c, name='rename') RESULT(err)
      IMPORT C_INT, C_CHAR
      CHARACTER(C_CHAR), INTENT(in) :: old(*), new(*)
      INTEGER(C_INT) :: err
    END FUNCTION

! ---&---1---------2---------3---------4---------5---------6---------7--
!
    FUNCTION c_mkdir(path, mode) BIND(c, name='mkdir') RESULT(err)
      IMPORT C_INT, C_CHAR
      CHARACTER(C_CHAR), INTENT(in) :: path(*)
      INTEGER(C_INT), VALUE :: mode
      INTEGER(C_INT) :: err
    END FUNCTION
  END INTERFACE
CONTAINS

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  FUNCTION rank_file(prefix, suffix) RESULT(name)
    CHARACTER(*), INTENT(in) :: prefix, suffix
    CHARACTER(:), ALLOCATABLE :: name
    CHARACTER(32) :: rank

    IF (my_id < 100) THEN
      WRITE (rank, '(i2)') my_id
    ELSE
      WRITE (rank, '(i0)') my_id
    END IF
    name = prefix//TRIM(rank)//suffix
  END FUNCTION

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE make_directory(path)
    CHARACTER(*), INTENT(in) :: path
    INTEGER :: err
    LOGICAL :: exists

    INQUIRE (file=path, exist=exists)
    IF (exists) RETURN
    err = c_mkdir(path//C_NULL_CHAR, INT(o'755', C_INT))
    INQUIRE (file=path, exist=exists)
    IF (.NOT. exists) CALL fail('cannot create directory: '//path)
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE atomic_finish(unit, path)
    INTEGER, INTENT(in) :: unit
    CHARACTER(*), INTENT(in) :: path
    INTEGER :: err
    FLUSH (unit, iostat=err)
    IF (err /= 0) CALL fail('cannot flush '//path)
    CLOSE (unit, iostat=err)
    IF (err /= 0) CALL fail('cannot close '//path)
    err = c_rename(path//'.tmp'//C_NULL_CHAR, path//C_NULL_CHAR)
    IF (err /= 0) CALL fail('cannot commit '//path)
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  FUNCTION dimensions(kind) RESULT(d)
    INTEGER, INTENT(in) :: kind
    INTEGER(INT64) :: d(12)

    IF (kind == 1) THEN
      d = [INT(NX, INT64), INT(NYP+1, INT64), INT(NZ, INT64), 3_INT64, INT(NZS2P1, INT64), INT(NXSP+1, INT64), &
           INT(my_id, INT64), INT(NPROC, INT64), INT(LSUBX(my_id), INT64), INT(NXSP+1, INT64), 0_INT64, INT(NYS+1, INT64)]
    ELSE
      d = [INT(NX, INT64), INT(NYP+1, INT64), INT(NZ, INT64), 4_INT64, 0_INT64, 0_INT64, INT(my_id, INT64), INT(NPROC, INT64), &
           0_INT64, INT(NX, INT64), INT(LSUBY(my_id), INT64), INT(NYP+1, INT64)]
    END IF
  END FUNCTION

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE header(unit, kind, step, time)
    INTEGER, INTENT(in) :: unit, kind
    INTEGER(INT64), INTENT(in) :: step
    REAL(dp), INTENT(in) :: time
    INTEGER :: ios

    WRITE (unit, iostat=ios) magic, 1_INT32, INT(kind, INT32), step, time, dimensions(kind)
    IF (ios /= 0) CALL fail('cannot write field header')
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE write_restart(directory, step, time, v)
    CHARACTER(*), INTENT(in) :: directory
    INTEGER(INT64), INTENT(in) :: step
    REAL(dp), INTENT(in) :: time
    COMPLEX(dp), INTENT(in) :: v(3, 0:NXSP, 0:NYS, -NZS:NZS)
    CHARACTER(:), ALLOCATABLE :: path
    INTEGER :: u, ios, x, y, z, c, err

    IF (output_format /= 'binary') THEN
      path = directory//'/'//rank_file('field_S_', '.dat')
      OPEN (newunit=u, file=path//'.tmp', status='replace', iostat=ios)
      IF (ios /= 0) CALL fail('cannot write '//path)
      WRITE (u, '(i0,1x,es26.17e3)') step, time
      DO z = -NZS, NZS
        DO y = 0, NYS
          DO x = 0, NXSP
            DO c = 1, 3
              WRITE (u, '("(",es26.17e3,",",es26.17e3,")")') REAL(v(c, x, y, z), dp), AIMAG(v(c, x, y, z))
            END DO
          END DO
        END DO
      END DO
      CALL atomic_finish(u, path)
    END IF
    IF (output_format /= 'ascii') THEN
      path = directory//'/'//rank_file('field_S_', '.bin')
      OPEN (newunit=u, file=path//'.tmp', status='replace', access='stream', form='unformatted', &
            convert='little_endian', iostat=ios)
      IF (ios /= 0) CALL fail('cannot write '//path)
      CALL header(u, 1, step, time)
      WRITE (u, iostat=ios) v
      IF (ios /= 0) CALL fail('cannot write restart payload')
      CALL atomic_finish(u, path)
    END IF
    CALL MPI_Barrier(MPI_COMM_WORLD, err)
    CALL mpi_check(err, 'restart commit barrier')
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE read_restart(directory, step, time, v)
    CHARACTER(*), INTENT(in) :: directory
    INTEGER(INT64), INTENT(out) :: step
    REAL(dp), INTENT(out) :: time
    COMPLEX(dp), INTENT(out) :: v(3, 0:NXSP, 0:NYS, -NZS:NZS)
    CHARACTER(:), ALLOCATABLE :: path, binary
    CHARACTER(8) :: probe
    COMPLEX(dp) :: extra
    REAL(dp) :: nan
    INTEGER :: u, ios, err
    INTEGER(INT32) :: version, kind
    INTEGER(INT64) :: d(12), length, lo, hi, bits
    LOGICAL :: text_exists, binary_exists

    path = directory//'/'//rank_file('field_S_', '.dat')
    binary = directory//'/'//rank_file('field_S_', '.bin')
    INQUIRE (file=path, exist=text_exists)
    INQUIRE (file=binary, exist=binary_exists)
    IF (.NOT. text_exists .OR. (output_format /= 'ascii' .AND. binary_exists)) path = binary
    OPEN (newunit=u, file=path, status='old', access='stream', form='unformatted', &
          convert='little_endian', iostat=ios)
    IF (ios /= 0) CALL fail('cannot open restart: '//path)
    READ (u, iostat=ios) probe
    IF (ios /= 0) CALL fail('truncated restart')
    IF (probe == magic) THEN
      READ (u, iostat=ios) version, kind, step, time, d
      IF (ios /= 0) CALL fail('truncated restart header')
      IF (version /= 1 .OR. kind /= 1 .OR. ANY(d /= dimensions(1))) CALL fail('restart header does not match mesh/rank')
      INQUIRE (unit=u, size=length)
      IF (length /= 128_INT64+16_INT64*SIZE(v, kind=INT64)) CALL fail('restart payload length mismatch')
      READ (u, iostat=ios) v
      IF (ios /= 0) CALL fail('truncated restart payload')
      CLOSE (u)
    ELSE
      CLOSE (u)
      OPEN (newunit=u, file=path, status='old', iostat=ios)
      IF (ios /= 0) CALL fail('cannot open text restart')
      READ (u, *, iostat=ios) step, time
      IF (ios /= 0) CALL fail('invalid text restart metadata')
      ! List-directed input accepts scalar and array records.
      ! Reading one extra value detects surplus tokens even on the final record.
      ! Unread/null entries remain NaN, so a truncated payload cannot pass at EOF.
      nan = ieee_value(0.0_DP, ieee_quiet_nan)
      v = CMPLX(nan, nan, dp)
      READ (u, *, iostat=ios) v, extra
      IF (ios >= 0) CALL fail('invalid or extra text restart payload')
      CLOSE (u)
    END IF
    IF (step < 0 .OR. .NOT. ieee_is_finite(time) .OR. .NOT. ALL(ieee_is_finite(REAL(v, dp))) .OR. &
        .NOT. ALL(ieee_is_finite(AIMAG(v)))) CALL fail('invalid or nonfinite restart')
    CALL MPI_Allreduce(step, lo, 1, MPI_INTEGER8, MPI_MIN, MPI_COMM_WORLD, err)
    CALL mpi_check(err, 'restart step minimum')
    CALL MPI_Allreduce(step, hi, 1, MPI_INTEGER8, MPI_MAX, MPI_COMM_WORLD, err)
    CALL mpi_check(err, 'restart step maximum')
    IF (lo /= hi) CALL fail('inconsistent restart steps')
    bits = TRANSFER(time, bits)
    CALL MPI_Allreduce(bits, lo, 1, MPI_INTEGER8, MPI_MIN, MPI_COMM_WORLD, err)
    CALL mpi_check(err, 'restart time minimum')
    CALL MPI_Allreduce(bits, hi, 1, MPI_INTEGER8, MPI_MAX, MPI_COMM_WORLD, err)
    CALL mpi_check(err, 'restart time maximum')
    IF (lo /= hi) CALL fail('inconsistent restart times')
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE snapshot(v, p, physical)
    COMPLEX(dp), INTENT(in) :: v(3, 0:NXSP, 0:NYS, -NZS:NZS), p(0:NXSP, 0:NYS, -NZS:NZS)
    REAL(dp), INTENT(out) :: physical(NX, 0:NYP, NZ, 4)
    INTEGER :: c

    DO c = 1, 3
      CALL channel_xyzfft_backward(physical(:, :, :, c), v(c, :, :, :))
    END DO
    CALL channel_xyzfft_backward(physical(:, :, :, 4), p)
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE write_physical(directory, step, time, physical)
    CHARACTER(*), INTENT(in) :: directory
    INTEGER(INT64), INTENT(in) :: step
    REAL(dp), INTENT(in) :: time, physical(NX, 0:NYP, NZ, 4)
    CHARACTER(:), ALLOCATABLE :: path
    INTEGER :: u, ios, x, y, z

    IF (output_format /= 'binary') THEN
      path = directory//'/'//rank_file('field_P_', '.dat')
      OPEN (newunit=u, file=path//'.tmp', status='replace', iostat=ios)
      IF (ios /= 0) CALL fail('cannot write '//path)
      WRITE (u, '(i0,1x,es26.17e3)') step, time
      DO z = 1, NZ
        DO y = 0, NYP
          DO x = 1, NX
            WRITE (u, '(3(es26.17e3,1x))') physical(x, y, z, 1:3)
          END DO
        END DO
      END DO
      WRITE (u, '(es26.17e3)') physical(:, :, :, 4)
      CALL atomic_finish(u, path)
    END IF
    IF (output_format /= 'ascii') THEN
      path = directory//'/'//rank_file('field_P_', '.bin')
      OPEN (newunit=u, file=path//'.tmp', status='replace', access='stream', form='unformatted', &
            convert='little_endian', iostat=ios)
      IF (ios /= 0) CALL fail('cannot write '//path)
      CALL header(u, 2, step, time)
      WRITE (u, iostat=ios) physical
      IF (ios /= 0) CALL fail('cannot write physical payload')
      CALL atomic_finish(u, path)
    END IF
  END SUBROUTINE
END MODULE
! ---&---1---------2---------3---------4---------5---------6---------7--
