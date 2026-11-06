! ---&---1---------2---------3---------4---------5---------6---------7--
PROGRAM test_vector_mesh
  USE dns_kinds
  USE dns_config
  USE mesh
  USE fft
  USE dns_vector
  USE dns_initial
  IMPLICIT NONE
  COMPLEX(dp), ALLOCATABLE :: scalar(:, :, :), RESULT(:, :, :), vector(:, :, :, :), &
                              gradp(:, :, :, :), rotation(:, :, :, :), bc(:, :, :, :)
  REAL(dp) :: local_error, global_error, flow
  INTEGER :: err, ranks, x, y, z, c, r, next
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
  ZL = 10
  Re = 2800
  INDEX_GRAD_PRE = 0
  theta_rad = 0
  EPS_init = 0
  V_Microblow0 = 0.001_DP
  CALL setup_mesh(XL, ZL)
  next = 0
  DO r = 0, ranks-1
    IF (LSUBX(r) /= next .OR. NXSPT(r) < 0) CALL fail('x partition continuity')
    next = next+NXSPT(r)+1
  END DO
  IF (next /= NXS+1) CALL fail('x partition coverage')
  next = 0
  DO r = 0, ranks-1
    IF (LSUBY(r) /= next .OR. NYPT(r) < 0) CALL fail('y partition continuity')
    next = next+NYPT(r)+1
  END DO
  IF (next /= NY+1) CALL fail('y partition coverage')
  CALL GET_COMMAND_ARGUMENT(1, directory)
  CALL setup_fft(XL, ZL, Re, TRIM(directory))
  ALLOCATE (scalar(0:NXSP, 0:NYS, -NZS:NZS), RESULT(0:NXSP, 0:NYS, -NZS:NZS), &
            vector(3, 0:NXSP, 0:NYS, -NZS:NZS), gradp(3, 0:NXSP, 0:NYS, -NZS:NZS), &
            rotation(3, 0:NXSP, 0:NYS, -NZS:NZS), bc(3, 0:NXSP, 0:1, -NZS:NZS))
  scalar = 0
  vector = 0
  DO z = -NZS, NZS
    DO y = 0, 4
      DO x = 0, NXSP
        scalar(x, y, z) = CMPLX(SIN(REAL(LSUBX(my_id)+x+y+z, dp)), COS(REAL(x-y+z, dp)), dp)
        DO c = 1, 3
          vector(c, x, y, z) = scalar(x, y, z)*CMPLX(REAL(c, dp), REAL(c-2, dp), dp)
        END DO
      END DO
    END DO
  END DO
  CALL Grad(gradp, scalar)
  CALL Curl(rotation, gradp)
  local_error = MAXVAL(ABS(rotation))
  CALL Curl(rotation, vector)
  CALL Div(RESULT, rotation)
  local_error = MAX(local_error, MAXVAL(ABS(RESULT)))
  CALL MPI_Allreduce(local_error, global_error, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, err)
  IF (global_error > 2.E-11_DP) CALL fail('curl(grad)=0 or div(curl)=0')
  CALL initialize_velocity(vector)
  CALL velocity_boundary(bc)
  local_error = 0
  IF (my_id == 0) THEN
    scalar = vector(1, :, :, :)
    CALL FOR_Q(flow, scalar)
    local_error = ABS(flow-1)
    local_error = MAX(local_error, ABS(bc(2, 0, 0, 0)-V_Microblow0), ABS(bc(2, 0, 1, 0)-V_Microblow0))
  ELSE
    local_error = MAXVAL(ABS(bc))
  END IF
  CALL MPI_Allreduce(local_error, global_error, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, err)
  IF (global_error > 2.E-13_DP) CALL fail('laminar flow rate / zero-mode wall ownership')
  CALL close_fft()
  IF (my_id == 0) PRINT *, 'PASS: mesh partition, vector identities, laminar flow and wall ownership'
  CALL MPI_Finalize(err)
END PROGRAM
! ---&---1---------2---------3---------4---------5---------6---------7--
