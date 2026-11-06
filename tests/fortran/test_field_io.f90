! ---&---1---------2---------3---------4---------5---------6---------7--
PROGRAM test_field_io
  USE dns_kinds
  USE dns_config
  USE mesh
  USE dns_io
  IMPLICIT NONE
  COMPLEX(dp), ALLOCATABLE :: original(:, :, :, :), restored(:, :, :, :)
  INTEGER(INT64), PARAMETER :: expected_step = 9000000000_INT64
  INTEGER(INT64) :: step
  REAL(dp) :: time
  INTEGER :: x, y, z, c, err, ranks, format_index
  CHARACTER(4096) :: directory

  CALL MPI_Init(err)
  CALL MPI_Comm_size(MPI_COMM_WORLD, ranks, err)
  NX = 8
  NY = 8
  NZ = 8
  NXS = 4
  NYS = 6
  NZS = 4
  NX2 = 12
  NY2 = 12
  NZ2 = 12
  INDEX_Serial = 0
  IF (ranks == 1) INDEX_Serial = 1
  CALL setup_mesh(10.0_DP, 10.0_DP)
  CALL GET_COMMAND_ARGUMENT(1, directory)
  ALLOCATE (original(3, 0:NXSP, 0:NYS, -NZS:NZS), restored(3, 0:NXSP, 0:NYS, -NZS:NZS))
  DO z = -NZS, NZS
    DO y = 0, NYS
      DO x = 0, NXSP
        DO c = 1, 3
          original(c, x, y, z) = CMPLX(REAL(c+3*(LSUBX(my_id)+x)+37*y+91*z, dp)/127, &
                                       REAL(c-2*x+11*y-53*z, dp)/113, dp)
        END DO
      END DO
    END DO
  END DO
  DO format_index = 1, 2
    output_format = 'ascii'
    IF (format_index == 2) output_format = 'binary'
    CALL write_restart(TRIM(directory), expected_step, 0.125_DP, original)
    restored = 0
    CALL read_restart(TRIM(directory), step, time, restored)
    IF (step /= expected_step .OR. time /= 0.125_DP) CALL fail('restart 64-bit step/time round trip')
    IF (ANY(original /= restored)) CALL fail('restart spectral complex payload round trip')
  END DO
  IF (my_id == 0) PRINT *, 'PASS: text/binary complex spectral ordering and 64-bit metadata'
  CALL MPI_Finalize(err)
END PROGRAM
! ---&---1---------2---------3---------4---------5---------6---------7--
