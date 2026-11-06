! ---&---1---------2---------3---------4---------5---------6---------7--
MODULE dns_fft_backend
  USE ISO_C_BINDING
  USE dns_kinds
  USE mesh, ONLY: fail
  IMPLICIT NONE
  PRIVATE
  INCLUDE 'fftw3.f03'
  CHARACTER(*), PARAMETER, PUBLIC :: fft_backend_name = 'FFTW'
  PUBLIC :: fft_engine, init_fft_engine, close_fft_engine, backend_wisdom
  PUBLIC :: backend_x_forward, backend_x_backward, backend_z, backend_y
  TYPE :: fft_engine
    INTEGER :: nx, ny, nz, nyp, physical_alignment
    REAL(dp), ALLOCATABLE :: physical(:, :, :)
    REAL(dp), POINTER, CONTIGUOUS :: yr(:) => NULL()
    COMPLEX(dp), ALLOCATABLE :: xhat(:, :, :)
    COMPLEX(dp), POINTER, CONTIGUOUS :: k(:, :, :) => NULL(), y(:, :, :) => NULL()
    TYPE(C_PTR) :: k_memory = C_NULL_PTR, y_memory = C_NULL_PTR
    TYPE(C_PTR) :: xp = C_NULL_PTR, xm = C_NULL_PTR, zp = C_NULL_PTR, zm = C_NULL_PTR, yp = C_NULL_PTR
  END TYPE
CONTAINS

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE init_fft_engine(w, nxg, nyg, nzg, ypart, nxs, nxsp, nzs, serial)
    TYPE(fft_engine), INTENT(out) :: w
    INTEGER, INTENT(in) :: nxg, nyg, nzg, ypart, nxs, nxsp, nzs
    LOGICAL, INTENT(in) :: serial
    INTEGER(C_INT) :: n(1), embed(1)
    INTEGER(C_SIZE_T) :: nk, nywork
    COMPLEX(dp), POINTER, CONTIGUOUS :: flat(:)
    TYPE(fftw_iodim) :: ydim(1), batches(2), zdim(1), zbatches(2)

    w%nx = nxg
    w%ny = nyg
    w%nz = nzg
    w%nyp = ypart
    ALLOCATE (w%physical(nxg, 0:ypart, nzg), w%xhat(0:nxg/2, 0:ypart, 0:nzg-1))
    nk = INT(nxs+1, C_SIZE_T)*(ypart+1)*(2*nzs+1)
    w%k_memory = fftw_alloc_complex(nk)
    IF (.NOT. C_ASSOCIATED(w%k_memory)) CALL fail('FFTW workspace allocation failed')
    CALL C_F_POINTER(w%k_memory, flat, [INT(nk)])
    w%k(0:nxs, 0:ypart, -nzs:nzs) => flat
    nywork = INT(nxsp+1, C_SIZE_T)*(nyg+1)*(2*nzs+1)
    IF (serial) THEN
      ! The serial transpose is a layout identity: use one allocation for both views.
      w%y => w%k
      CALL C_F_POINTER(w%k_memory, w%yr, [INT(2*nywork)])
    ELSE
      w%y_memory = fftw_alloc_complex(nywork)
      IF (.NOT. C_ASSOCIATED(w%y_memory)) CALL fail('FFTW y-workspace allocation failed')
      CALL C_F_POINTER(w%y_memory, flat, [INT(nywork)])
      w%y(0:nxsp, 0:nyg, -nzs:nzs) => flat
      CALL C_F_POINTER(w%y_memory, w%yr, [INT(2*nywork)])
    END IF
    n = nxg
    embed = nxg
    w%xp = fftw_plan_many_dft_r2c(1, n, (ypart+1)*nzg, w%physical, embed, 1, nxg, &
                                  w%xhat, [nxg/2+1], 1, nxg/2+1, FFTW_MEASURE)
    w%xm = fftw_plan_many_dft_c2r(1, n, (ypart+1)*nzg, w%xhat, [nxg/2+1], 1, nxg/2+1, &
                                  w%physical, embed, 1, nxg, FFTW_MEASURE)
    w%physical_alignment = fftw_alignment_of(w%physical)
    ! Include all y planes in the z plan so alignment and vectorization are
    ! planned for the complete layout, rather than a new-array call per plane.
    zdim(1) = fftw_iodim(nzg, (nxg/2+1)*(ypart+1), (nxg/2+1)*(ypart+1))
    zbatches(1) = fftw_iodim(nxs+1, 1, 1)
    zbatches(2) = fftw_iodim(ypart+1, nxg/2+1, nxg/2+1)
    w%zp = fftw_plan_guru_dft(1, zdim, 2, zbatches, w%xhat, w%xhat, FFTW_FORWARD, FFTW_MEASURE)
    w%zm = fftw_plan_guru_dft(1, zdim, 2, zbatches, w%xhat, w%xhat, FFTW_BACKWARD, FFTW_MEASURE)
    ! Batch all real/imaginary x lanes and z planes without extracting y lines.
    ydim(1) = fftw_iodim(nyg+1, 2*(nxsp+1), 2*(nxsp+1))
    batches(1) = fftw_iodim(2*(nxsp+1), 1, 1)
    batches(2) = fftw_iodim((2*nzs+1), 2*(nxsp+1)*(nyg+1), 2*(nxsp+1)*(nyg+1))
    w%yp = fftw_plan_guru_r2r(1, ydim, 2, batches, w%yr, w%yr, [FFTW_REDFT00], FFTW_MEASURE)
    IF (.NOT. C_ASSOCIATED(w%xp) .OR. .NOT. C_ASSOCIATED(w%xm) .OR. .NOT. C_ASSOCIATED(w%zp) .OR. &
        .NOT. C_ASSOCIATED(w%zm) .OR. .NOT. C_ASSOCIATED(w%yp)) CALL fail('FFTW plan creation failed')
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE close_fft_engine(w)
    TYPE(fft_engine), INTENT(inout) :: w

    CALL fftw_destroy_plan(w%xp)
    CALL fftw_destroy_plan(w%xm)
    CALL fftw_destroy_plan(w%zp)
    CALL fftw_destroy_plan(w%zm)
    CALL fftw_destroy_plan(w%yp)
    CALL fftw_free(w%k_memory)
    IF (C_ASSOCIATED(w%y_memory)) CALL fftw_free(w%y_memory)
    NULLIFY (w%k, w%y, w%yr)
    w%k_memory = C_NULL_PTR
    w%y_memory = C_NULL_PTR
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE backend_wisdom(path, exporting)
    CHARACTER(*), INTENT(in) :: path
    LOGICAL, INTENT(in) :: exporting
    INTEGER :: ignored

    IF (exporting) THEN
      ignored = fftw_export_wisdom_to_filename(path//C_NULL_CHAR)
    ELSE
      ignored = fftw_import_wisdom_from_filename(path//C_NULL_CHAR)
    END IF
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE backend_z(w, forward)
    TYPE(fft_engine), INTENT(inout) :: w
    LOGICAL, INTENT(in) :: forward

    IF (forward) THEN
      CALL fftw_execute_dft(w%zp, w%xhat, w%xhat)
    ELSE
      CALL fftw_execute_dft(w%zm, w%xhat, w%xhat)
    END IF
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE backend_y(w)
    TYPE(fft_engine), INTENT(inout) :: w

    CALL fftw_execute_r2r(w%yp, w%yr, w%yr)
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE backend_x_forward(w, physical)
    TYPE(fft_engine), INTENT(inout) :: w
    REAL(dp), INTENT(inout) :: physical(w%nx, 0:w%nyp, w%nz)

    IF (fftw_alignment_of(physical) == w%physical_alignment) THEN
      CALL fftw_execute_dft_r2c(w%xp, physical, w%xhat)
    ELSE
      ! External buffers may have another alignment class; retain a safe fallback.
      w%physical = physical
      CALL fftw_execute_dft_r2c(w%xp, w%physical, w%xhat)
    END IF
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE backend_x_backward(w, physical)
    TYPE(fft_engine), INTENT(inout) :: w
    REAL(dp), INTENT(out) :: physical(w%nx, 0:w%nyp, w%nz)

    IF (fftw_alignment_of(physical) == w%physical_alignment) THEN
      CALL fftw_execute_dft_c2r(w%xm, w%xhat, physical)
    ELSE
      CALL fftw_execute_dft_c2r(w%xm, w%xhat, w%physical)
      physical = w%physical
    END IF
  END SUBROUTINE
END MODULE
! ---&---1---------2---------3---------4---------5---------6---------7--
