! ---&---1---------2---------3---------4---------5---------6---------7--
PROGRAM test_zero_mode
  USE dns_kinds
  USE dns_config
  USE mesh
  USE fft
  USE eqsolve_pv
  USE dns_equations
  USE dns_timestep
  USE dns_vector, ONLY: DERIV
  IMPLICIT NONE
  COMPLEX(dp), ALLOCATABLE :: velocity(:, :, :, :), forcing(:, :, :, :), bcv(:, :, :, :), &
                              p(:, :, :), bcp(:, :, :), gradient(:, :, :)
  REAL(dp) :: local_error, global_error
  INTEGER :: ierr, ranks
  CHARACTER(4096) :: directory

  CALL MPI_Init(ierr)
  CALL MPI_Comm_size(MPI_COMM_WORLD, ranks, ierr)
  NX = 16
  NY = 16
  NZ = 16
  NXS = 6
  NYS = 12
  NZS = 6
  NX2 = 24
  NY2 = 24
  NZ2 = 24
  XL = 10
  ZL = 10
  Re = 2800
  G_pre = 0
  a_V = 1
  b_V = 0
  INDEX_Serial = 0
  IF (ranks == 1) INDEX_Serial = 1
  CALL GET_COMMAND_ARGUMENT(1, directory)
  CALL setup_mesh(XL, ZL)
  CALL setup_fft(XL, ZL, Re, TRIM(directory))
  CALL allocate_equations()
  G_V = 1500*Re
  a_pre = 1
  b_pre = 0
  CALL EQSOLVE_PV_setup()
  a_pre = 0
  b_pre = 1
  ALLOCATE (velocity(3, 0:NXSP, 0:NYS, -NZS:NZS), forcing(3, 0:NXSP, 0:NYS, -NZS:NZS), &
            bcv(3, 0:NXSP, 0:1, -NZS:NZS), p(0:NXSP, 0:NYS, -NZS:NZS), bcp(0:NXSP, 0:1, -NZS:NZS), &
            gradient(0:NXSP, 0:NYS, -NZS:NZS))
  forcing = 0
  bcv = 0
  bcp = 0
  IF (my_id == 0) THEN
    ! Manufactured pressure p=y+y^3/4 = 19*T1/16 + T3/16.
    ! Normal momentum requires p_y = RHS_y - (G_V/Re)*v_y.
    forcing(2, 0, 0, 0) = 1.375_DP+1500*0.001_DP
    forcing(2, 0, 2, 0) = 0.375_DP
    bcv(2, 0, :, 0) = 0.001_DP
    bcp(0, 0, 0) = 1.75_DP
    bcp(0, 1, 0) = -1.75_DP
  END IF
  CALL STEP23(p, bcp, velocity, forcing, bcv)
  CALL DERIV(2, gradient, p)
  local_error = 0
  IF (my_id == 0) THEN
    local_error = MAXVAL(ABS(velocity(2, 0, 1:NYS, 0)))
    local_error = MAX(local_error, ABS(velocity(2, 0, 0, 0)-0.001_DP))
    gradient(0, 0, 0) = gradient(0, 0, 0)-1.375_DP
    gradient(0, 2, 0) = gradient(0, 2, 0)-0.375_DP
    local_error = MAX(local_error, MAXVAL(ABS(gradient(0, :, 0))))
    WRITE (*, *) 'zero-mode manufactured momentum max error = ', local_error
  END IF
  CALL MPI_Allreduce(local_error, global_error, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  IF (global_error > 1.E-10_DP) CALL fail('zero-mode pressure gradient / normal momentum constraint')
  CALL close_fft()
  IF (my_id == 0) PRINT *, 'PASS: zero-mode manufactured pressure and normal velocity'
  CALL MPI_Finalize(ierr)
END PROGRAM
! ---&---1---------2---------3---------4---------5---------6---------7--
