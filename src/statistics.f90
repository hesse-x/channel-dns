! ---&---1---------2---------3---------4---------5---------6---------7--
! Volume diagnostics and plane-averaged velocity moments.
MODULE dns_statistics
  USE, INTRINSIC :: IEEE_ARITHMETIC
  USE dns_kinds
  USE dns_config, ONLY: Nstat
  USE mesh
  USE fft
  USE dns_vector, ONLY: Div
  USE dns_io, ONLY: rank_file, atomic_finish, make_directory
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: close_statistics, setup_statistics, evaluate_fields, accumulate_statistics, write_statistics
  COMPLEX(dp), ALLOCATABLE :: divergence(:, :, :)
  REAL(dp), ALLOCATABLE :: divp(:, :, :)
  REAL(dp), ALLOCATABLE :: weights(:), moments(:, :), average(:, :)
  INTEGER(INT64) :: samples = 0
  REAL(dp) :: dpdx_mean = 0, dpdz_mean = 0
CONTAINS

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE setup_statistics()
    INTEGER :: j, g, k
    REAL(dp) :: theta, VALUE, pi
    ALLOCATE (weights(0:NYP), moments(20, 0:NYP), average(20, 0:NYP))
    ALLOCATE (divergence(0:NXSP, 0:NYS, -NZS:NZS), divp(NX, 0:NYP, NZ))
    average = 0
    pi = ACOS(-1.0_DP)
    DO j = 0, NYP
      g = LSUBY(my_id)+j
      theta = pi*g/NY
      IF (g == 0 .OR. g == NY) THEN
        VALUE = 1.0_DP/(REAL(NY, dp)**2-1)
      ELSE
        VALUE = 1
        DO k = 1, (NY-1)/2
          VALUE = VALUE-2*COS(2*k*theta)/(4.0_DP*k*k-1)
        END DO
        VALUE = 2*(VALUE-COS(NY*theta)/(REAL(NY, dp)**2-1))/NY
      END IF
      weights(j) = VALUE/2
    END DO
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE close_statistics()
    IF (ALLOCATED(weights)) DEALLOCATE (weights, moments, average, divergence, divp)
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE evaluate_fields(directory, step, time, physical, v, p, check, append)
    CHARACTER(*), INTENT(in) :: directory
    INTEGER(INT64), INTENT(in) :: step
    REAL(dp), INTENT(in) :: time, physical(NX, 0:NYP, NZ, 4)
    COMPLEX(dp), INTENT(in) :: v(3, 0:NXSP, 0:NYS, -NZS:NZS), p(0:NXSP, 0:NYS, -NZS:NZS)
    LOGICAL, INTENT(in) :: check, append
    REAL(dp) :: u1, v1, w1, p1, uu, vv, ww, total(3), local(3), variance
    INTEGER :: x, y, z, c, err, unit, ios
    LOGICAL :: finite, allfinite, exists

    moments = 0
    DO z = 1, NZ
      DO y = 0, NYP
        DO x = 1, NX
          u1 = physical(x, y, z, 1)
          v1 = physical(x, y, z, 2)
          w1 = physical(x, y, z, 3)
          uu = u1*u1
          vv = v1*v1
          ww = w1*w1
          p1 = physical(x, y, z, 4)-0.5_DP*(uu+vv+ww)
          moments(:, y) = moments(:, y)+[u1, v1, w1, p1, uu, vv, ww, u1*v1, v1*w1, w1*u1, &
                                         p1*u1, p1*v1, p1*w1, p1*p1, uu*u1, vv*v1, ww*w1, uu*uu, vv*vv, ww*ww]
        END DO
      END DO
    END DO
    IF (.NOT. check) RETURN
    CALL Div(divergence, v)
    CALL channel_xyzfft('B', divp, divergence)
    finite = ALL(ieee_is_finite(physical)) .AND. ALL(ieee_is_finite(moments)) .AND. ALL(ieee_is_finite(divp)) &
             .AND. ALL(ieee_is_finite(REAL(v, dp))) .AND. ALL(ieee_is_finite(AIMAG(v))) &
             .AND. ALL(ieee_is_finite(REAL(p, dp))) .AND. ALL(ieee_is_finite(AIMAG(p)))
    local = 0
    local(1) = SUM(moments(6, :))
    DO y = 0, NYP
      variance = 0
      DO c = 1, 3
        variance = variance+MAX(0.0_DP, moments(4+c, y)/NXTNZ-(moments(c, y)/NXTNZ)**2)
      END DO
      local(2) = local(2)+0.5_DP*weights(y)*variance
      local(3) = local(3)+weights(y)*SUM(divp(:, y, :)**2)/NXTNZ
    END DO
    CALL MPI_Reduce(local, total, 3, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, err)
    CALL mpi_check(err, 'diagnostic reduction')
    IF (my_id == 0) finite = finite .AND. ALL(ieee_is_finite(total))
    CALL MPI_Allreduce(finite, allfinite, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, err)
    CALL mpi_check(err, 'diagnostic finite check')
    IF (.NOT. allfinite) CALL fail('non-finite field or diagnostic')
    IF (my_id /= 0) RETURN
    total(3) = SQRT(total(3))
    WRITE (*, '(a,i0,4(a,es24.16e3),a)') 'step = ', step, ' TIME = ', time, ' EKT = ', total(1), &
      ' TKE = ', total(2), ' div_L2 = ', total(3), ' finite = true'
    IF (.NOT. append) RETURN
    INQUIRE (file=directory//'/diagnostics.csv', exist=exists)
    OPEN (newunit=unit, file=directory//'/diagnostics.csv', position='append', status='unknown', iostat=ios)
    IF (ios /= 0) CALL fail('cannot append diagnostics.csv')
    IF (.NOT. exists) WRITE (unit, '(a)') 'step,time,ekt_sum_v2,tke_plane_fluctuation_volume_mean,divergence_l2,finite'
    WRITE (unit, '(i0,4(",",es26.17e3),",1")') step, time, total
    CLOSE (unit, iostat=ios)
    IF (ios /= 0) CALL fail('cannot close diagnostics.csv')
    OPEN (newunit=unit, file=directory//'/'//rank_file('EKT', '.dat'), position='append', status='unknown', iostat=ios)
    IF (ios /= 0) CALL fail('cannot append EKT')
    WRITE (unit, '(2(es26.17e3,1x))') time, total(1)
    CLOSE (unit, iostat=ios)
    IF (ios /= 0) CALL fail('cannot close EKT')
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE accumulate_statistics(dpdx, dpdz)
    REAL(dp), INTENT(in) :: dpdx, dpdz

    average = (average*REAL(samples, dp)+moments/REAL(NXTNZ, dp))/REAL(samples+1, dp)
    dpdx_mean = (dpdx_mean*REAL(samples, dp)+dpdx)/REAL(samples+1, dp)
    dpdz_mean = (dpdz_mean*REAL(samples, dp)+dpdz)/REAL(samples+1, dp)
    samples = samples+1
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE write_statistics(directory)
    CHARACTER(*), INTENT(in) :: directory
    CHARACTER(:), ALLOCATABLE :: path
    CHARACTER :: number
    INTEGER :: d, u, ios, y

    CALL make_directory(directory//'/Result')
    DO d = 0, 4
      WRITE (number, '(i1)') d
      CALL make_directory(directory//'/Result/'//number)
    END DO
    IF (my_id == 0) THEN
      path = directory//'/Result/0/DPDX_DPDZ.dat'
      OPEN (newunit=u, file=path//'.tmp', status='replace', iostat=ios)
      IF (ios /= 0) CALL fail('cannot write statistics')
      WRITE (u, '(i0)') samples
      WRITE (u, '(2(es26.17e3,1x))') dpdx_mean, dpdz_mean
      CALL atomic_finish(u, path)
    END IF
    DO d = 1, 4
      WRITE (number, '(i1)') d
      IF (d == 1) THEN
        path = directory//'/Result/1/'//rank_file('V_Pre_mean_Field_', '.dat')
      ELSE
        path = directory//'/Result/'//number//'/'//rank_file('V_Pre_moment'//number//'_Field_', '.dat')
      END IF
      OPEN (newunit=u, file=path//'.tmp', status='replace', iostat=ios)
      IF (ios /= 0) CALL fail('cannot write statistics')
      DO y = 0, NYP
        SELECT CASE (d)
        CASE (1)
          WRITE (u, '(*(f19.8))') YC(y), average(1:4, y)
        CASE (2)
          WRITE (u, '(*(f19.8))') YC(y), average(5:10, y), average(14, y), average(11:13, y)
        CASE (3)
          WRITE (u, '(*(f19.8))') YC(y), average(15:17, y)
        CASE (4)
          WRITE (u, '(*(f19.8))') YC(y), average(18:20, y)
        END SELECT
      END DO
      CALL atomic_finish(u, path)
    END DO
  END SUBROUTINE
END MODULE
! ---&---1---------2---------3---------4---------5---------6---------7--
