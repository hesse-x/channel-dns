! ---&---1---------2---------3---------4---------5---------6---------7--
! MPI setup, time integration, diagnostics and output.
PROGRAM channel_dns
  USE dns_kinds
  USE dns_config
  USE mesh
  USE fft
  USE eqsolve_pv
  USE dns_equations, ONLY: EQSOLVE_PV_setup
  USE dns_timestep, ONLY: NLTERM, STEP23, prepare_timestep_workspace, close_timestep_workspace
  USE dns_initial
  USE dns_io
  USE dns_statistics
  IMPLICIT NONE
  COMPLEX(dp), ALLOCATABLE :: v(:, :, :, :), previous(:, :, :, :), nl(:, :, :, :), oldnl(:, :, :, :), bcv(:, :, :, :)
  COMPLEX(dp), ALLOCATABLE :: p(:, :, :), bcp(:, :, :), oldbcp(:, :, :)
  REAL(dp), ALLOCATABLE :: physical(:, :, :, :)
  REAL(dp) :: compute_elapsed, step_start, compute_max
  REAL(dp) :: time, a0, a1, b0, b1, dpdx, dpdz, qx, qz, gradients(2), start, elapsed, maximum
  INTEGER(INT64) :: nt, step, completed
  INTEGER :: err, i, it, effective, remaining, argn, threads, ios, physical_cadence, statistics_cadence, warmup
  CHARACTER(4096) :: arg
  CHARACTER(:), ALLOCATABLE :: directory
  LOGICAL :: benchmark, timing, print_fft, diagnostic_step, restart_step, physical_step, check_step, sample

  CALL MPI_Init(err)
  CALL mpi_check(err, 'MPI initialization')
  CALL MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN, err)
  CALL mpi_check(err, 'MPI error handler')
  directory = ''
  warmup = 0
  benchmark = .FALSE.
  timing = .FALSE.
  print_fft = .FALSE.
  argn = COMMAND_ARGUMENT_COUNT()
  i = 1
  DO WHILE (i <= argn)
    CALL GET_COMMAND_ARGUMENT(i, arg)
    SELECT CASE (TRIM(arg))
    CASE ('--benchmark-warmup')
      benchmark = .TRUE.
      i = i+1
      IF (i > argn) CALL fail('--benchmark-warmup requires a nonnegative integer')
      CALL GET_COMMAND_ARGUMENT(i, arg)
      READ (arg, *, iostat=ios) warmup
      IF (ios /= 0) CALL fail('invalid benchmark warmup')
      IF (warmup < 0) CALL fail('invalid benchmark warmup')
    CASE ('--timing')
      timing = .TRUE.
    CASE ('--print-fft')
      print_fft = .TRUE.
    CASE ('--threads', '-t')
      ! Each MPI rank computes on its main thread.
      i = i+1
      IF (i > argn) CALL fail('--threads requires a positive integer')
      CALL GET_COMMAND_ARGUMENT(i, arg)
      READ (arg, *, iostat=ios) threads
      IF (ios /= 0) CALL fail('--threads requires a positive integer')
      IF (threads < 1) CALL fail('--threads requires a positive integer')
    CASE ('--help', '-h')
      WRITE (*, '(a)') 'usage: channel_dns [--timing] [--print-fft] <run-directory>'
      CALL MPI_Finalize(err)
      STOP
    CASE default
      IF (arg(1:1) == '-' .OR. LEN(directory) > 0) CALL fail('usage: channel_dns [--timing] <run-directory>')
      directory = TRIM(arg)
    END SELECT
    i = i+1
  END DO
  IF (LEN(directory) == 0) CALL fail('run-directory is required')
  CALL load_config(directory)
  CALL setup_mesh(XL, ZL)
  IF (my_id == 0) WRITE (*, '(a,i0,a)') 'channel-dns: Fortran, ', NPROC, ' MPI rank(s), single-thread '//fft_backend_name
  CALL setup_fft(XL, ZL, Re, directory)
  IF (print_fft .AND. my_id == 0) THEN
    WRITE (*, *) 'WAVX ', WAVX
    WRITE (*, *) 'WAVZ ', WAVZ
    WRITE (*, *) 'WL ', WL
    WRITE (*, *) 'WD ', WD
    WRITE (*, *) 'WR ', WR
  END IF
  CALL allocate_equations()
  CALL setup_statistics()
  CALL prepare_timestep_workspace()
  ALLOCATE (v(3, 0:NXSP, 0:NYS, -NZS:NZS), previous(3, 0:NXSP, 0:NYS, -NZS:NZS), &
            nl(3, 0:NXSP, 0:NYS, -NZS:NZS), oldnl(3, 0:NXSP, 0:NYS, -NZS:NZS), bcv(3, 0:NXSP, 0:1, -NZS:NZS))
  ALLOCATE (p(0:NXSP, 0:NYS, -NZS:NZS), bcp(0:NXSP, 0:1, -NZS:NZS), oldbcp(0:NXSP, 0:1, -NZS:NZS))
  ALLOCATE (physical(NX, 0:NYP, NZ, 4))
  v = 0
  previous = 0
  nl = 0
  oldnl = 0
  bcv = 0
  p = 0
  bcp = 0
  oldbcp = 0
  nt = 0
  time = 0
  dpdx = 0
  dpdz = 0
  remaining = NCOUNT_disturb
  IF (INDEX_CONT == 0) THEN
    CALL initialize_velocity(v)
  ELSE
    CALL read_restart(directory, nt, time, v)
    IF (remaining > 0) THEN
      completed = MIN(INT(NCOUNT_disturb, INT64), nt/NDISTURBFRQ)
      remaining = remaining-INT(completed)
      IF (nt > 0 .AND. MOD(nt, INT(NDISTURBFRQ, INT64)) == 0 .AND. nt/NDISTURBFRQ <= NCOUNT_disturb) THEN
        CALL disturb_velocity(v, EPS_disturb)
      END IF
    END IF
  END IF
  CALL velocity_boundary(bcv)
  physical_cadence = physical_interval
  IF (physical_cadence < 0) physical_cadence = NPRNFRQ
  statistics_cadence = statistics_interval
  IF (statistics_cadence < 0) statistics_cadence = NPRNFRQ
  IF (timing) CALL MPI_Barrier(MPI_COMM_WORLD, err)
  start = MPI_Wtime()
  compute_elapsed = 0
  DO it = 1, NSTEPS
    IF (benchmark .AND. it == warmup+1) THEN
      CALL MPI_Barrier(MPI_COMM_WORLD, err)
      CALL mpi_check(err, 'benchmark warmup barrier')
    END IF
    IF (benchmark) step_start = MPI_Wtime()
    effective = MIN(it, 2)
    a0 = REAL(effective, dp)/DT
    a1 = REAL(1-effective, dp)/(2*DT)
    b0 = REAL(effective, dp)
    b1 = REAL(1-effective, dp)
    IF (it < 3) THEN
      G_V = (a0+a1)*Re
      a_pre = 1
      b_pre = 0
      CALL EQSOLVE_PV_setup()
      a_pre = a_pre0
      b_pre = b_pre0
    END IF
    step = nt+it
    time = time+DT
    CALL NLTERM(nl, v, bcp)
    oldbcp = b0*bcp+b1*oldbcp
    oldnl = a0*v+a1*previous+b0*nl+b1*oldnl
    CALL swap_vectors(v, previous)
    IF (INDEX_UNSTEADY_BC == 1) CALL velocity_boundary(bcv)
    CALL STEP23(p, oldbcp, v, oldnl, bcv)
    CALL swap_vectors(nl, oldnl)
    CALL swap_scalars(bcp, oldbcp)
    IF (my_id == 0) THEN
      qx = 0
      qz = 0
      DO i = 0, NYS, 2
        qx = qx-REAL(v(1, 0, i, 0), dp)/REAL(i*i-1, dp)
        qz = qz-REAL(v(3, 0, i, 0), dp)/REAL(i*i-1, dp)
      END DO
      IF (INDEX_GRAD_PRE == 1) THEN
        dpdx = grad_pre_x
        dpdz = grad_pre_z
      ELSE
        IF (QX_DPDX_1 == 0 .OR. QZ_DPDZ_1 == 0) CALL fail('zero unit pressure-gradient flow response')
        dpdx = (qxmean-qx)/QX_DPDX_1
        dpdz = (qzmean-qz)/QZ_DPDZ_1
      END IF
      v(1, 0, :, 0) = v(1, 0, :, 0)+dpdx*U_DPDX_1
      v(3, 0, :, 0) = v(3, 0, :, 0)+dpdz*W_DPDZ_1
    END IF
    gradients = [dpdx, dpdz]
    CALL MPI_Bcast(gradients, 2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, err)
    CALL mpi_check(err, 'pressure gradient broadcast')
    dpdx = gradients(1)
    dpdz = gradients(2)
    IF (benchmark .AND. it > warmup) compute_elapsed = compute_elapsed+MPI_Wtime()-step_start
    diagnostic_step = MOD(step, INT(NPRNFRQ, INT64)) == 0
    restart_step = MOD(step, INT(NDUMFRQ, INT64)) == 0
    physical_step = .FALSE.
    IF (physical_cadence > 0) physical_step = MOD(step, INT(physical_cadence, INT64)) == 0
    check_step = diagnostic_step .OR. restart_step .OR. it == NSTEPS
    sample = step > Nstat
    IF (check_step .OR. physical_step .OR. sample) THEN
      CALL snapshot(v, p, physical)
      CALL evaluate_fields(directory, step, time, physical, v, p, check_step, diagnostic_step)
      IF (physical_step) CALL write_physical(directory, step, time, physical)
      IF (sample) THEN
        CALL accumulate_statistics(dpdx, dpdz)
        IF (MOD(step, INT(statistics_cadence, INT64)) == 0 .OR. it == NSTEPS) CALL write_statistics(directory)
      END IF
    END IF
    IF (restart_step) CALL write_restart(directory, step, time, v)
    IF (remaining > 0) THEN
      IF (MOD(step, INT(NDISTURBFRQ, INT64)) == 0) THEN
        CALL disturb_velocity(v, EPS_disturb)
        remaining = remaining-1
        IF (my_id == 0) WRITE (*, '(a,i0,a,i0,a)') 'disturbance at step ', step, ' (', remaining, ' remaining)'
      END IF
    END IF
  END DO
  elapsed = MPI_Wtime()-start
  IF (benchmark) THEN
    CALL MPI_Reduce(compute_elapsed, compute_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, err)
    CALL mpi_check(err, 'benchmark compute time')
    CALL MPI_Reduce(elapsed, maximum, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, err)
    CALL mpi_check(err, 'benchmark loop time')
    IF (my_id == 0) THEN
      WRITE (*, '(a,es24.16,a)') 'benchmark compute wall time = ', compute_max, ' s'
      WRITE (*, '(a,es24.16,a)') 'benchmark loop wall time = ', maximum, ' s'
      WRITE (*, '(a,i0)') 'benchmark measured steps = ', MAX(0, NSTEPS-warmup)
    END IF
  END IF
  IF (timing) THEN
    CALL MPI_Reduce(elapsed, maximum, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, err)
    CALL mpi_check(err, 'timestep timing')
    IF (my_id == 0) WRITE (*, '(a,es16.8,a)') 'channel-dns: timestep wall time = ', maximum, ' s'
  END IF
  CALL close_statistics()
  CALL close_timestep_workspace()
  CALL close_fft()
  CALL MPI_Finalize(err)
CONTAINS

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE swap_vectors(a, b)
    COMPLEX(dp), ALLOCATABLE, INTENT(inout) :: a(:, :, :, :), b(:, :, :, :)
    COMPLEX(dp), ALLOCATABLE :: tmp(:, :, :, :)

    CALL MOVE_ALLOC(a, tmp)
    CALL MOVE_ALLOC(b, a)
    CALL MOVE_ALLOC(tmp, b)
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE swap_scalars(a, b)
    COMPLEX(dp), ALLOCATABLE, INTENT(inout) :: a(:, :, :), b(:, :, :)
    COMPLEX(dp), ALLOCATABLE :: tmp(:, :, :)

    CALL MOVE_ALLOC(a, tmp)
    CALL MOVE_ALLOC(b, a)
    CALL MOVE_ALLOC(tmp, b)
  END SUBROUTINE
END PROGRAM
! ---&---1---------2---------3---------4---------5---------6---------7--
