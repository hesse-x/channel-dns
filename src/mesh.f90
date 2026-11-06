! ---&---1---------2---------3---------4---------5---------6---------7--
! Fourier mode partitions and physical-space collocation points.
MODULE mesh
  USE dns_kinds
  USE dns_config, ONLY: NX, NY, NZ, NXS, NYS, NZS, NX2, NY2, NZ2, INDEX_Serial
  IMPLICIT NONE
  INTEGER :: NXH, NZH, NXS2P1, NZS2P1, NXSP, NXSP2P1, NYP, NY2P, NXTNZ
  INTEGER :: my_id, NPROC, INDEX_MPI0, NPART
  INTEGER, ALLOCATABLE :: LSUBX(:), NXSPT(:), LSUBX2(:), NXSP2P1T(:), LSUBY(:), NYPT(:), LSUBY2(:), NY2PT(:)
  REAL(dp), ALLOCATABLE :: XC(:), YC(:), ZC(:)
CONTAINS

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE partition(n, beg, upper)
    INTEGER, INTENT(in) :: n
    INTEGER, ALLOCATABLE, INTENT(out) :: beg(:), upper(:)
    INTEGER :: r
    ALLOCATE (beg(0:NPART-1), upper(0:NPART-1))
    upper = n/NPART
    IF (MOD(n, NPART) > 0) upper(NPART-MOD(n, NPART):) = upper(NPART-MOD(n, NPART):)+1
    beg(0) = 0
    DO r = 1, NPART-1
      beg(r) = beg(r-1)+upper(r-1)
    END DO
    upper = upper-1
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE setup_mesh(xl, zl)
    REAL(dp), INTENT(in) :: xl, zl
    INTEGER :: err, j, g
    REAL(dp) :: pi

    CALL MPI_Comm_rank(MPI_COMM_WORLD, my_id, err)
    CALL mpi_check(err, 'mesh rank')
    CALL MPI_Comm_size(MPI_COMM_WORLD, NPROC, err)
    CALL mpi_check(err, 'mesh size')
    NPART = NPROC
    INDEX_MPI0 = my_id
    NXH = NX/2
    NZH = NZ/2
    NXS2P1 = 2*NXS+1
    NZS2P1 = 2*NZS+1
    NXTNZ = NX*NZ
    CALL partition(NXS+1, LSUBX, NXSPT)
    CALL partition(NY+1, LSUBY, NYPT)
    CALL partition(NY2+1, LSUBY2, NY2PT)
    LSUBX2 = 2*LSUBX
    NXSP2P1T = 2*NXSPT+1
    NXSP = NXSPT(my_id)
    NXSP2P1 = 2*NXSP+1
    NYP = NYPT(my_id)
    NY2P = NY2PT(my_id)
    ALLOCATE (XC(NX), YC(0:NYP), ZC(NZ))
    pi = ACOS(-1.0_DP)
    DO j = 1, NX
      XC(j) = j*(xl/NX)
    END DO
    DO j = 1, NZ
      ZC(j) = j*(zl/NZ)
    END DO
    DO j = 0, NYP
      g = LSUBY(my_id)+j
      IF (g < NY/2) THEN
        YC(j) = COS(pi*g/NY)
      ELSE IF (g == NY/2) THEN
        YC(j) = 0
      ELSE
        YC(j) = -COS(pi*(NY-g)/NY)
      END IF
    END DO
  END SUBROUTINE
END MODULE
! ---&---1---------2---------3---------4---------5---------6---------7--
