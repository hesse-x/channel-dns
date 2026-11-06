! ---&---1---------2---------3---------4---------5---------6---------7--
! Numerical kinds, MPI errors and collective failure handling.
MODULE dns_kinds
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: dp => C_DOUBLE
  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: INT64, ERROR_UNIT
  USE mpi_f08
  IMPLICIT NONE
CONTAINS

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE fail(message)
    CHARACTER(*), INTENT(in) :: message
    INTEGER :: ierr, rank

    CALL MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    WRITE (ERROR_UNIT, '(a,i0,2a)') 'rank ', rank, ': ', message
    CALL MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    error STOP 1
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE mpi_check(ierr, operation)
    INTEGER, INTENT(in) :: ierr
    CHARACTER(*), INTENT(in) :: operation
    CHARACTER(MPI_MAX_ERROR_STRING) :: message
    INTEGER :: n, err

    IF (ierr == MPI_SUCCESS) RETURN
    CALL MPI_Error_string(ierr, message, n, err)
    CALL fail(operation//': '//message(:n))
  END SUBROUTINE
END MODULE
! ---&---1---------2---------3---------4---------5---------6---------7--
