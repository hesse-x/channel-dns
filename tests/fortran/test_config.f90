! ---&---1---------2---------3---------4---------5---------6---------7--
PROGRAM test_config
  USE dns_kinds
  USE dns_config
  USE dns_toml, ONLY: toml_config, toml_table, table, has_key, read_integer, read_real, read_string, read_vector
  IMPLICIT NONE
  CHARACTER(4096) :: directory, source
  CHARACTER(:), ALLOCATABLE :: text
  REAL(dp) :: vector(3), real_value
  INTEGER :: unit
  TYPE(toml_config), TARGET :: document
  TYPE(toml_table), POINTER :: child
  INTEGER :: err, VALUE

  CALL MPI_Init(err)
  CALL GET_COMMAND_ARGUMENT(1, directory)
  CALL load_config(TRIM(directory))
  IF (NX /= 16 .OR. NZ /= 16 .OR. NX2 /= 24 .OR. Nstat /= 50000) CALL fail('decimal integer conversion')
  IF (Re /= 2800.0_DP .OR. XL /= 10.0_DP .OR. DT /= 0.003_DP) CALL fail('TOML integer/real conversion')
  IF (ANY(a_V /= 1.0_DP) .OR. ANY(b_V /= 0.0_DP)) CALL fail('inline numeric array conversion')
  IF (output_format /= 'binary') CALL fail('quoted string conversion')
  IF (physical_interval /= 0 .OR. statistics_interval /= 2) CALL fail('TOML output configuration')
  IF (ABS(theta_rad) > EPSILON(0.0_DP) .OR. qxmean /= 1.0_DP) CALL fail('derived physical parameters')
  ! Reloading and repeated close must release the owned tree safely.
  CALL document%load(TRIM(directory)//'/config.toml')
  CALL document%load(TRIM(directory)//'/config.toml')
  child => table(document%handle, 'MESH')
  IF (.NOT. ASSOCIATED(child)) CALL fail('case-insensitive TOML table lookup')
  IF (.NOT. has_key(child, 'nx')) CALL fail('case-insensitive TOML key lookup')
  CALL read_integer(child, 'Nx', VALUE)
  IF (VALUE /= 16) CALL fail('typed TOML integer access')
  child => table(document%handle, 'missing')
  IF (ASSOCIATED(child)) CALL fail('missing TOML table lookup')
  IF (has_key(document%handle, 'missing')) CALL fail('TOML lookup must not insert missing tables')
  CALL document%CLOSE()
  CALL document%CLOSE()
  IF (ALLOCATED(document%handle)) CALL fail('TOML document close')
  CALL GET_COMMAND_ARGUMENT(2, source)
  CALL document%load(TRIM(source)//'/examples/config.toml')
  child => table(document%handle, 'mesh')
  CALL read_integer(child, 'nx', VALUE)
  IF (VALUE /= 32) CALL fail('example config.toml')
  CALL document%load(TRIM(source)//'/examples/Re_tau180/Re_tau180.toml')
  child => table(document%handle, 'parm')
  CALL read_real(child, 're', real_value)
  IF (real_value /= 180.0_DP) CALL fail('Re_tau180 configuration')
  OPEN (newunit=unit, file='config-parser-test.toml', status='replace')
  WRITE (unit, '(a)') '[MiXeD] # section comment', &
    'label = "a#b=c" # outside comment', "single = 'literal'", &
    ACHAR(9)//'number'//ACHAR(9)//'='//ACHAR(9)//'+17', &
    'vector = [+1, -2.5, 3E+1]'
  CLOSE (unit)
  CALL document%load('config-parser-test.toml')
  child => table(document%handle, 'mixed')
  CALL read_string(child, 'label', text)
  IF (text /= 'a#b=c') CALL fail('quoted comment marker and equals sign')
  CALL read_string(child, 'single', text)
  IF (text /= 'literal') CALL fail('literal string')
  CALL read_integer(child, 'number', VALUE)
  IF (VALUE /= 17) CALL fail('tabs and signed integer')
  CALL read_vector(child, 'vector', vector)
  IF (ANY(vector /= [1.0_DP, -2.5_DP, 30.0_DP])) CALL fail('signed numeric array')
  CALL document%CLOSE()
  OPEN (newunit=unit, file='config-parser-test.toml', status='old')
  CLOSE (unit, status='delete')
  PRINT *, 'PASS: configuration subset and typed Fortran access'
  CALL MPI_Finalize(err)
END PROGRAM
! ---&---1---------2---------3---------4---------5---------6---------7--
