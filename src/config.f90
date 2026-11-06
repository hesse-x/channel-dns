! ---&---1---------2---------3---------4---------5---------6---------7--
! TOML and text parameters, validation and derived constants.
MODULE dns_config
  USE dns_kinds
  USE dns_toml
  USE, INTRINSIC :: IEEE_ARITHMETIC
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: load_config, validate_config, lower
  INTEGER :: NX, NY, NZ, INDEX_disturb_Z, NXS, NYS, NZS, NX2, NY2, NZ2, INDEX_Serial
  PUBLIC :: NX, NY, NZ, INDEX_disturb_Z, NXS, NYS, NZS, NX2, NY2, NZ2, INDEX_Serial
  INTEGER :: INDEX_GRAD_PRE, INDEX_UNSTEADY_BC
  PUBLIC :: INDEX_GRAD_PRE, INDEX_UNSTEADY_BC
  REAL(dp) :: XL, ZL, THETA, Re, G_pre, a_pre0, b_pre0
  PUBLIC :: XL, ZL, THETA, Re, G_pre, a_pre0, b_pre0
  REAL(dp), DIMENSION(3) :: a_V, b_V
  PUBLIC :: a_V, b_V
  INTEGER :: INDEX_CONT, NSTEPS, NDUMFRQ, NPRNFRQ, NCOUNT_disturb, NDISTURBFRQ, Nstat
  PUBLIC :: INDEX_CONT, NSTEPS, NDUMFRQ, NPRNFRQ, NCOUNT_disturb, NDISTURBFRQ, Nstat
  REAL(dp) :: EPS_init, DT, EPS_disturb, V_Microblow0
  PUBLIC :: EPS_init, DT, EPS_disturb, V_Microblow0
  INTEGER, PUBLIC :: physical_interval = -1, statistics_interval = -1
  CHARACTER(6), PUBLIC :: output_format = 'ascii'
  REAL(dp), PUBLIC :: theta_rad, grad_pre_x, grad_pre_z, qxmean, qzmean
  INTEGER, PARAMETER, PUBLIC :: E_ijk(3, 2) = RESHAPE([2, 3, 1, 3, 1, 2], [3, 2])
CONTAINS

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  FUNCTION lower(s) RESULT(t)
    CHARACTER(*), INTENT(in) :: s
    CHARACTER(LEN(s)) :: t
    INTEGER :: i, k

    t = s
    DO i = 1, LEN(s)
      k = IACHAR(s(i:i))
      IF (k >= 65 .AND. k <= 90) t(i:i) = ACHAR(k+32)
    END DO
  END FUNCTION

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE load_config(directory)
    CHARACTER(*), INTENT(in) :: directory
    TYPE(toml_config), TARGET :: document
    TYPE(toml_table), POINTER :: child, root
    CHARACTER(:), ALLOCATABLE :: text
    INTEGER :: unit, ios, file_index, err
    INTEGER(INT64) :: fingerprint, minimum, maximum
    CHARACTER :: ch
    LOGICAL :: exists
    CHARACTER(48), PARAMETER :: files(4) = [CHARACTER(48) :: &
                                            'config.toml', 'meshparm_channel_with_Microblow.txt', &
                                            'parm_channel_with_Microblow.txt', 'control_channel_with_Microblow.txt']

    INQUIRE (file=directory//'/config.toml', exist=exists)
    fingerprint = 0
    DO file_index = 1, 4
      IF (exists .NEQV. (file_index == 1)) CYCLE
      OPEN (newunit=unit, file=directory//'/'//TRIM(files(file_index)), status='old', &
            access='stream', form='unformatted', iostat=ios)
      IF (ios /= 0) CALL fail('cannot open '//TRIM(files(file_index)))
      DO
        READ (unit, iostat=ios) ch
        IF (ios /= 0) EXIT
        fingerprint = IEOR(ISHFTC(fingerprint, 5), INT(IACHAR(ch), INT64))
      END DO
      CLOSE (unit)
    END DO
    CALL MPI_Allreduce(fingerprint, minimum, 1, MPI_INTEGER8, MPI_MIN, MPI_COMM_WORLD, err)
    CALL mpi_check(err, 'parameter fingerprint minimum')
    CALL MPI_Allreduce(fingerprint, maximum, 1, MPI_INTEGER8, MPI_MAX, MPI_COMM_WORLD, err)
    CALL mpi_check(err, 'parameter fingerprint maximum')
    IF (minimum /= maximum) CALL fail('parameter files differ between MPI ranks')
    IF (.NOT. exists) THEN
      CALL load_legacy(directory)
    ELSE
      CALL document%load(directory//'/config.toml')
      root => document%handle
      CALL check_keys(root, 'mesh|parm|control|output', 'root')
      CALL section(root, 'mesh', child)
      CALL check_keys(child, 'nx|ny|nz|index_disturb_z|nxs|nys|nzs|nx2|ny2|nz2|index_serial', 'mesh')
      CALL read_integer(child, 'nx', NX)
      CALL read_integer(child, 'ny', NY)
      CALL read_integer(child, 'nz', NZ)
      CALL read_integer(child, 'index_disturb_z', INDEX_disturb_Z)
      CALL read_integer(child, 'nxs', NXS)
      CALL read_integer(child, 'nys', NYS)
      CALL read_integer(child, 'nzs', NZS)
      CALL read_integer(child, 'nx2', NX2)
      CALL read_integer(child, 'ny2', NY2)
      CALL read_integer(child, 'nz2', NZ2)
      CALL read_integer(child, 'index_serial', INDEX_Serial)
      CALL section(root, 'parm', child)
      CALL check_keys(child, 'index_grad_pre|index_unsteady_bc|xl|zl|theta|re|g_pre|a_pre0|b_pre0|a_v|b_v', 'parm')
      CALL read_integer(child, 'index_grad_pre', INDEX_GRAD_PRE)
      CALL read_integer(child, 'index_unsteady_bc', INDEX_UNSTEADY_BC)
      CALL read_real(child, 'xl', XL)
      CALL read_real(child, 'zl', ZL)
      CALL read_real(child, 'theta', THETA)
      CALL read_real(child, 're', Re)
      CALL read_real(child, 'g_pre', G_pre)
      CALL read_real(child, 'a_pre0', a_pre0)
      CALL read_real(child, 'b_pre0', b_pre0)
      CALL read_vector(child, 'a_v', a_V)
      CALL read_vector(child, 'b_v', b_V)
      CALL section(root, 'control', child)
      CALL check_keys(child, &
                      'index_cont|nsteps|ndumfrq|nprnfrq|ncount_disturb|ndisturbfrq|nstat|eps_init|dt|eps_disturb|v_microblow0', &
                      'control')
      CALL read_integer(child, 'index_cont', INDEX_CONT)
      CALL read_integer(child, 'nsteps', NSTEPS)
      CALL read_integer(child, 'ndumfrq', NDUMFRQ)
      CALL read_integer(child, 'nprnfrq', NPRNFRQ)
      CALL read_integer(child, 'ncount_disturb', NCOUNT_disturb)
      CALL read_integer(child, 'ndisturbfrq', NDISTURBFRQ)
      CALL read_integer(child, 'nstat', Nstat)
      CALL read_real(child, 'eps_init', EPS_init)
      CALL read_real(child, 'dt', DT)
      CALL read_real(child, 'eps_disturb', EPS_disturb)
      CALL read_real(child, 'v_microblow0', V_Microblow0)
      IF (has_key(root, 'output')) THEN
        CALL section(root, 'output', child)
        CALL check_keys(child, 'format|physical_interval|statistics_interval', 'output')
        CALL read_string(child, 'format', text)
        text = lower(text)
        IF (text /= 'ascii' .AND. text /= 'binary' .AND. text /= 'both') CALL fail('invalid output.format')
        output_format = text
        IF (has_key(child, 'physical_interval')) THEN
          CALL read_integer(child, 'physical_interval', physical_interval)
          IF (physical_interval < 0) CALL fail('physical_interval must be nonnegative')
        END IF
        IF (has_key(child, 'statistics_interval')) THEN
          CALL read_integer(child, 'statistics_interval', statistics_interval)
          IF (statistics_interval <= 0) CALL fail('statistics_interval must be positive')
        END IF
      END IF
      CALL document%CLOSE()
    END IF
    theta_rad = THETA*ACOS(-1.0_DP)/180.0_DP
    grad_pre_x = 0
    grad_pre_z = 0
    qxmean = 0
    qzmean = 0
    IF (INDEX_GRAD_PRE == 1) THEN
      grad_pre_x = -COS(theta_rad)
      grad_pre_z = -SIN(theta_rad)
    ELSE
      qxmean = COS(theta_rad)
      qzmean = SIN(theta_rad)
    END IF
    CALL validate_config()
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE section(root, name, child)
    TYPE(toml_table), INTENT(inout) :: root
    CHARACTER(*), INTENT(in) :: name
    TYPE(toml_table), POINTER, INTENT(out) :: child

    child => table(root, name)
    IF (.NOT. ASSOCIATED(child)) CALL fail('missing TOML table: '//name)
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE load_legacy(directory)
    CHARACTER(*), INTENT(in) :: directory
    INTEGER :: u, ios, i

    OPEN (newunit=u, file=directory//'/meshparm_channel_with_Microblow.txt', status='old', iostat=ios)
    IF (ios /= 0) CALL fail('missing legacy mesh parameters')
    READ (u, *, iostat=ios)
    IF (ios /= 0) CALL fail('invalid legacy parameter record')
    READ (u, *, iostat=ios) NX, NY, NZ, INDEX_disturb_Z
    IF (ios /= 0) CALL fail('invalid legacy parameter record')
    READ (u, *, iostat=ios)
    IF (ios /= 0) CALL fail('invalid legacy parameter record')
    READ (u, *, iostat=ios) NXS, NYS, NZS
    IF (ios /= 0) CALL fail('invalid legacy parameter record')
    READ (u, *, iostat=ios)
    IF (ios /= 0) CALL fail('invalid legacy parameter record')
    READ (u, *, iostat=ios) NX2, NY2, NZ2
    IF (ios /= 0) CALL fail('invalid legacy parameter record')
    READ (u, *, iostat=ios)
    IF (ios /= 0) CALL fail('invalid legacy parameter record')
    READ (u, *, iostat=ios) INDEX_Serial
    IF (ios /= 0) CALL fail('invalid legacy parameter record')
    IF (ios /= 0) CALL fail('invalid legacy mesh parameters')
    CLOSE (u)
    OPEN (newunit=u, file=directory//'/parm_channel_with_Microblow.txt', status='old', iostat=ios)
    IF (ios /= 0) CALL fail('missing legacy physical parameters')
    READ (u, *, iostat=ios)
    IF (ios /= 0) CALL fail('invalid legacy parameter record')
    READ (u, *, iostat=ios) XL, ZL, THETA
    IF (ios /= 0) CALL fail('invalid legacy parameter record')
    READ (u, *, iostat=ios)
    IF (ios /= 0) CALL fail('invalid legacy parameter record')
    READ (u, *, iostat=ios) INDEX_GRAD_PRE, Re
    IF (ios /= 0) CALL fail('invalid legacy parameter record')
    READ (u, *, iostat=ios)
    IF (ios /= 0) CALL fail('invalid legacy parameter record')
    READ (u, *, iostat=ios) INDEX_UNSTEADY_BC
    IF (ios /= 0) CALL fail('invalid legacy parameter record')
    READ (u, *, iostat=ios)
    IF (ios /= 0) CALL fail('invalid legacy parameter record')
    READ (u, *, iostat=ios) G_pre, a_pre0, b_pre0
    IF (ios /= 0) CALL fail('invalid legacy parameter record')
    READ (u, *, iostat=ios)
    IF (ios /= 0) CALL fail('invalid legacy parameter record')
    DO i = 1, 3
      READ (u, *, iostat=ios) a_V(i), b_V(i)
      IF (ios /= 0) CALL fail('invalid legacy boundary parameters')
    END DO
    CLOSE (u)
    OPEN (newunit=u, file=directory//'/control_channel_with_Microblow.txt', status='old', iostat=ios)
    IF (ios /= 0) CALL fail('missing legacy control parameters')
    READ (u, *, iostat=ios)
    IF (ios /= 0) CALL fail('invalid legacy parameter record')
    READ (u, *, iostat=ios) INDEX_CONT, EPS_init
    IF (ios /= 0) CALL fail('invalid legacy parameter record')
    READ (u, *, iostat=ios)
    IF (ios /= 0) CALL fail('invalid legacy parameter record')
    READ (u, *, iostat=ios) NSTEPS, NDUMFRQ, NPRNFRQ, DT
    IF (ios /= 0) CALL fail('invalid legacy parameter record')
    READ (u, *, iostat=ios)
    IF (ios /= 0) CALL fail('invalid legacy parameter record')
    READ (u, *, iostat=ios) NCOUNT_disturb, NDISTURBFRQ, EPS_disturb
    IF (ios /= 0) CALL fail('invalid legacy parameter record')
    READ (u, *, iostat=ios)
    IF (ios /= 0) CALL fail('invalid legacy parameter record')
    READ (u, *, iostat=ios) V_Microblow0
    IF (ios /= 0) CALL fail('invalid legacy parameter record')
    READ (u, *, iostat=ios)
    IF (ios /= 0) CALL fail('invalid legacy parameter record')
    READ (u, *, iostat=ios) Nstat
    IF (ios /= 0) CALL fail('invalid legacy parameter record')
    IF (ios /= 0) CALL fail('invalid legacy control parameters')
    CLOSE (u)
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE validate_config()
    INTEGER :: ranks, err

    CALL MPI_Comm_size(MPI_COMM_WORLD, ranks, err)
    CALL mpi_check(err, 'configuration MPI size')
    IF (INDEX_Serial /= 0 .AND. INDEX_Serial /= 1) CALL fail('INDEX_Serial must be 0 or 1')
    IF ((INDEX_Serial == 1) .NEQV. (ranks == 1)) CALL fail('INDEX_Serial must match MPI rank count')
    IF (MIN(NX, NZ, NX2, NZ2) < 2 .OR. MIN(NY, NY2) < 4) CALL fail('grid dimensions too small')
    IF (MOD(NX, 2) /= 0 .OR. MOD(NZ, 2) /= 0 .OR. MOD(NX2, 2) /= 0 .OR. MOD(NZ2, 2) /= 0 .OR. MOD(NY, 2) /= 0) &
      CALL fail('Fourier grids and coarse NY must be even')
    IF (NXS < 0 .OR. NXS > MIN(NX, NX2)/2 .OR. NZS < 0 .OR. NZS > MIN(NZ, NZ2)/2) CALL fail('invalid Fourier truncation')
    IF (NYS < 4 .OR. NYS > MIN(NY, NY2)) CALL fail('invalid Chebyshev truncation')
    IF (ranks > MIN(NXS+1, NY+1, NY2+1)) CALL fail('MPI size exceeds decomposition extent')
    IF (NDUMFRQ <= 0 .OR. NPRNFRQ <= 0 .OR. DT <= 0 .OR. Re <= 0 .OR. XL <= 0 .OR. ZL <= 0 .OR. NSTEPS < 0) &
      CALL fail('invalid run parameters')
    IF (NCOUNT_disturb < 0) CALL fail('NCOUNT_disturb must be nonnegative')
    IF (NCOUNT_disturb > 0 .AND. NDISTURBFRQ <= 0) CALL fail('NDISTURBFRQ must be positive')
    IF (INDEX_CONT /= 0 .AND. INDEX_CONT /= 1) CALL fail('INDEX_CONT must be 0 or 1')
    IF (INDEX_UNSTEADY_BC /= 0 .AND. INDEX_UNSTEADY_BC /= 1) CALL fail('INDEX_UNSTEADY_BC must be 0 or 1')
    IF (output_format /= 'ascii' .AND. output_format /= 'binary' .AND. output_format /= 'both') CALL fail('invalid output.format')
    IF (.NOT. ieee_is_finite(XL)) CALL fail('nonfinite parameter: XL')
    IF (.NOT. ieee_is_finite(ZL)) CALL fail('nonfinite parameter: ZL')
    IF (.NOT. ieee_is_finite(THETA)) CALL fail('nonfinite parameter: THETA')
    IF (.NOT. ieee_is_finite(Re)) CALL fail('nonfinite parameter: Re')
    IF (.NOT. ieee_is_finite(G_pre)) CALL fail('nonfinite parameter: G_pre')
    IF (.NOT. ieee_is_finite(a_pre0)) CALL fail('nonfinite parameter: a_pre0')
    IF (.NOT. ieee_is_finite(b_pre0)) CALL fail('nonfinite parameter: b_pre0')
    IF (.NOT. ALL(ieee_is_finite(a_V))) CALL fail('nonfinite parameter: a_V')
    IF (.NOT. ALL(ieee_is_finite(b_V))) CALL fail('nonfinite parameter: b_V')
    IF (.NOT. ieee_is_finite(EPS_init)) CALL fail('nonfinite parameter: EPS_init')
    IF (.NOT. ieee_is_finite(DT)) CALL fail('nonfinite parameter: DT')
    IF (.NOT. ieee_is_finite(EPS_disturb)) CALL fail('nonfinite parameter: EPS_disturb')
    IF (.NOT. ieee_is_finite(V_Microblow0)) CALL fail('nonfinite parameter: V_Microblow0')
  END SUBROUTINE
END MODULE
! ---&---1---------2---------3---------4---------5---------6---------7--
