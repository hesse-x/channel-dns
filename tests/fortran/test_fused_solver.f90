! ---&---1---------2---------3---------4---------5---------6---------7--
! Compare the optimized tile path with the independent staged tau/influence path.
PROGRAM test_fused_solver
  USE dns_kinds
  USE dns_config
  USE mesh
  USE fft
  USE eqsolve_pv
  USE dns_equations
  USE dns_timestep
  USE dns_vector, ONLY: Div
  IMPLICIT NONE
  COMPLEX(dp), ALLOCATABLE :: forcing(:, :, :, :), bcv(:, :, :, :), bcp(:, :, :), rhs_pressure(:, :, :)
  COMPLEX(dp), ALLOCATABLE :: reference_p(:, :, :), actual_p(:, :, :), reference_v(:, :, :, :), actual_v(:, :, :, :), sigma(:, :, :)
  REAL(dp) :: local_error, global_error, scale
  INTEGER :: err, ranks, x, y, z, c, scenario, q
  CHARACTER(4096) :: directory

  CALL MPI_Init(err)
  CALL MPI_Comm_size(MPI_COMM_WORLD, ranks, err)
  NX = 32
  NY = 16
  NZ = 16
  NXS = 16
  NYS = 12
  NZS = 6
  NX2 = 48
  NY2 = 24
  NZ2 = 24
  XL = 10
  ZL = 7
  Re = 2800
  INDEX_Serial = 0
  IF (ranks == 1) INDEX_Serial = 1
  CALL GET_COMMAND_ARGUMENT(1, directory)
  CALL setup_mesh(XL, ZL)
  CALL setup_fft(XL, ZL, Re, TRIM(directory))
  CALL allocate_equations()
  ALLOCATE (forcing(3, 0:NXSP, 0:NYS, -NZS:NZS), reference_v(3, 0:NXSP, 0:NYS, -NZS:NZS), &
            actual_v(3, 0:NXSP, 0:NYS, -NZS:NZS), bcv(3, 0:NXSP, 0:1, -NZS:NZS), &
            bcp(0:NXSP, 0:1, -NZS:NZS), rhs_pressure(0:NXSP, 0:NYS, -NZS:NZS), &
            reference_p(0:NXSP, 0:NYS, -NZS:NZS), actual_p(0:NXSP, 0:NYS, -NZS:NZS), sigma(0:NXSP, 0:1, -NZS:NZS))
  DO z = -NZS, NZS
    DO y = 0, NYS
      DO x = 0, NXSP
        DO c = 1, 3
          forcing(c, x, y, z) = CMPLX(SIN(REAL(3*(LSUBX(my_id)+x)+7*y+z+c, dp)), &
                                      COS(REAL(LSUBX(my_id)+x-2*y+5*z-c, dp)), dp)*0.001_DP
        END DO
      END DO
    END DO
  END DO
  bcv = 0
  bcp = 0
  IF (my_id == 0) bcv(2, 0, :, 0) = 0.001_DP
  DO scenario = 1, 3
    G_pre = 0
    G_V = 1500*Re
    a_V = 1
    b_V = 0
    IF (scenario == 2) b_V = [0.1_DP, 0.0_DP, 0.2_DP]
    IF (scenario == 3) G_pre = 0.25_DP
    a_pre = 1
    b_pre = 0
    CALL EQSOLVE_PV_setup()
    a_pre = 0
    b_pre = 1
    IF (scenario == 3) a_pre = 0.25_DP
    CALL Div(rhs_pressure, forcing)
    CALL EQSOLVE_PV_2(reference_p, rhs_pressure, bcp, reference_v, forcing, bcv)
    sigma(:, 0, :) = -WN_V*reference_v(2, :, NYS, :)-Re*forcing(2, :, NYS, :)
    sigma(:, 1, :) = -WN_V*reference_v(2, :, NYS-1, :)+Re*(REAL(2*NYS, dp)*reference_p(:, NYS, :)-forcing(2, :, NYS-1, :))
    sigma = sigma/SIGMA0
    DO y = 0, NYS
      q = MOD(y+1, 2)
      reference_p(:, y, :) = reference_p(:, y, :)+sigma(:, q, :)*PRE0(:, y, :)
      reference_v(1, :, y, :) = reference_v(1, :, y, :)+sigma(:, q, :)*V0(1, :, y, :)
      reference_v(2, :, y, :) = reference_v(2, :, y, :)+sigma(:, 1-q, :)*V0(2, :, y, :)
      reference_v(3, :, y, :) = reference_v(3, :, y, :)+sigma(:, q, :)*V0(3, :, y, :)
    END DO
    IF (my_id == 0) THEN
      reference_v(2, 0, 0, 0) = 0.5_DP*(bcv(2, 0, 0, 0)+bcv(2, 0, 1, 0))
      reference_v(2, 0, 1:NYS, 0) = 0
    END IF
    CALL STEP23(actual_p, bcp, actual_v, forcing, bcv)
    scale = 1+MAX(MAXVAL(ABS(reference_p)), MAXVAL(ABS(reference_v)))
    local_error = MAX(MAXVAL(ABS(actual_p-reference_p)), MAXVAL(ABS(actual_v-reference_v)))/scale
    CALL MPI_Allreduce(local_error, global_error, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, err)
    IF (global_error > 2.E-11_DP) THEN
      PRINT *, 'scenario, normalized error: ', scenario, global_error
      CALL fail('fused/staged pressure-velocity solve mismatch')
    END IF
  END DO
  CALL close_timestep_workspace()
  CALL close_fft()
  IF (my_id == 0) PRINT *, 'PASS: fused/staged solve, full tiles, tails, zero mode and Robin boundaries'
  CALL MPI_Finalize(err)
END PROGRAM
! ---&---1---------2---------3---------4---------5---------6---------7--
