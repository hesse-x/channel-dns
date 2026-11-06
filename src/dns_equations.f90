! ---&---1---------2---------3---------4---------5---------6---------7--
! Fourier-Chebyshev pressure and velocity solvers.
MODULE dns_equations
  USE dns_kinds, ONLY: dp
  USE dns_vector
  IMPLICIT NONE
  INTEGER, PARAMETER :: mode_tile = 8
  COMPLEX(dp), ALLOCATABLE, SAVE :: tile_p(:, :), tile_dy(:, :), tile_force(:, :), tile_v(:, :, :)
  PRIVATE
  PUBLIC :: prepare_fused_workspace, close_fused_workspace, solve_modes_fused
  PUBLIC :: EQSOLVE_PV_setup, WN_AD, Infulence_coefficient, B_PROBLEM, EQSOLVE_UW_DPDXZ_1, EQSOLVE_PV_1, EQSOLVE_PV_2, SOLVE, SOLVE2
CONTAINS

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE prepare_fused_workspace()
    USE mesh, ONLY: NYS

    IF (ALLOCATED(tile_p)) RETURN
    ALLOCATE (tile_p(mode_tile, 0:NYS), tile_dy(mode_tile, 0:NYS), &
              tile_force(mode_tile, 0:NYS), tile_v(mode_tile, 0:NYS, 3))
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE close_fused_workspace()
    IF (ALLOCATED(tile_p)) DEALLOCATE (tile_p, tile_dy, tile_force, tile_v)
  END SUBROUTINE
! Eight adjacent x modes share one ordered tau recurrence.
! The lane index is contiguous; work arrays are local to the tile.

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE solve_tau_tile(out, rhs, g, a, b, bc, wn, ad, xb, z, lanes, zero_lane)
    USE mesh
    USE fft, ONLY: WL, WD, WR
    COMPLEX(dp), INTENT(out) :: out(mode_tile, 0:NYS)
    COMPLEX(dp), INTENT(in) :: rhs(mode_tile, 0:NYS), bc(0:, 0:, -NZS:)
    REAL(dp), INTENT(in) :: g, a, b, wn(0:, -NZS:), ad(0:, 2:, -NZS:)
    INTEGER, INTENT(in) :: xb, z, lanes, zero_lane
    REAL(dp) :: weights(mode_tile), ad01(mode_tile, 0:1), func
    COMPLEX(dp) :: rhs0(mode_tile)
    INTEGER :: l, x, y, q

    DO l = 1, lanes
      x = xb+l-1
      weights(l) = wn(x, z)
      rhs0(l) = rhs(l, 0)
    END DO
    IF (g == 0 .AND. a == 0 .AND. zero_lane > 0) THEN
      rhs0(zero_lane) = 0
      DO y = 2, NYS, 2
        rhs0(zero_lane) = rhs0(zero_lane)+rhs(zero_lane, y)/REAL(y*y-1, dp)
      END DO
      x = xb+zero_lane-1
      IF (b /= 0) rhs0(zero_lane) = rhs0(zero_lane)+(bc(x, 0, z)+bc(x, 1, z))/(2*b)
    END IF
    DO l = 1, lanes
      x = xb+l-1
      out(l, 0) = 0.5_DP*(bc(x, 0, z)+bc(x, 1, z))
      out(l, 1) = 0.5_DP*(bc(x, 0, z)-bc(x, 1, z))
      out(l, 2) = WL(2)*rhs0(l)+WD(2)*rhs(l, 2)+WR(2)*rhs(l, 4)
    END DO
    DO y = 3, NYS-2
      DO l = 1, lanes
        out(l, y) = WL(y)*rhs(l, y-2)+WD(y)*rhs(l, y)+WR(y)*rhs(l, y+2)
      END DO
    END DO
    DO l = 1, lanes
      out(l, NYS-1) = WL(NYS-1)*rhs(l, NYS-3)+WD(NYS-1)*rhs(l, NYS-1)
      out(l, NYS) = WL(NYS)*rhs(l, NYS-2)+WD(NYS)*rhs(l, NYS)
    END DO
    DO y = NYS-2, 2, -1
      DO l = 1, lanes
        out(l, y) = out(l, y)-WR(y)*weights(l)*ad(xb+l-1, y+2, z)*out(l, y+2)
      END DO
    END DO
    DO y = NYS, NYS-1, -1
      q = MOD(y, 2)
      func = a+b*REAL(y*y, dp)
      IF (a == 0 .AND. b == 0) func = REAL((1-q)*y*y+q, dp)
      ad01(1:lanes, q) = func
    END DO
    DO y = NYS-2, 0, -1
      q = MOD(y, 2)
      func = a+b*REAL(y*y, dp)
      IF (a == 0 .AND. b == 0) func = REAL((1-q)*y*y+q, dp)
      DO l = 1, lanes
        out(l, q) = out(l, q)-ad01(l, q)*ad(xb+l-1, y+2, z)*out(l, y+2)
        ad01(l, q) = func-WL(y+2)*weights(l)*ad01(l, q)*ad(xb+l-1, y+2, z)
      END DO
    END DO
    IF (g == 0 .AND. a == 0 .AND. zero_lane > 0) THEN
      out(zero_lane, 0) = 0
      ad01(zero_lane, 0) = 1
    END IF
    DO l = 1, lanes
      out(l, 0) = out(l, 0)/ad01(l, 0)
      out(l, 1) = out(l, 1)/ad01(l, 1)
    END DO
    DO y = 2, NYS
      DO l = 1, lanes
        out(l, y) = (out(l, y)-WL(y)*weights(l)*out(l, y-2))*ad(xb+l-1, y, z)
      END DO
    END DO
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE solve_modes_fused(pressure, pressure_bc, velocity, rhs_velocity, velocity_bc)
    USE mesh
    USE eqsolve_pv
    USE dns_config, ONLY: Re
    USE fft, ONLY: CWAVX, CWAVZ
    COMPLEX(dp), INTENT(out) :: pressure(0:NXSP, 0:NYS, -NZS:NZS), velocity(3, 0:NXSP, 0:NYS, -NZS:NZS)
    COMPLEX(dp), INTENT(in) :: pressure_bc(0:NXSP, 0:1, -NZS:NZS), &
                               rhs_velocity(3, 0:NXSP, 0:NYS, -NZS:NZS), velocity_bc(3, 0:NXSP, 0:1, -NZS:NZS)
    COMPLEX(dp) :: divwall(mode_tile, 0:1), cc(mode_tile, 2), sigma(mode_tile, 0:1), tail(mode_tile, 2)
    COMPLEX(dp) :: tangent, normal, ptop, vtop, vprev
    INTEGER :: z, xb, lanes, l, x, y, c, q, zero_lane

    CALL prepare_fused_workspace()
    DO z = -NZS, NZS
      DO xb = 0, NXSP, mode_tile
        lanes = MIN(mode_tile, NXSP+1-xb)
        zero_lane = 0
        IF (INDEX_MPI0 == 0 .AND. xb == 0 .AND. z == 0) zero_lane = 1
        DO c = 1, 3
          DO y = 0, NYS
            DO l = 1, lanes
              tile_v(l, y, c) = rhs_velocity(c, xb+l-1, y, z)
            END DO
          END DO
        END DO
        tile_dy(1:lanes, NYS) = 0
        DO l = 1, lanes
          tile_dy(l, NYS-1) = REAL(2*NYS, dp)*tile_v(l, NYS, 2)
        END DO
        DO y = NYS-2, 0, -1
          DO l = 1, lanes
            tile_dy(l, y) = tile_dy(l, y+2)+REAL(2*(y+1), dp)*tile_v(l, y+1, 2)
          END DO
        END DO
        tile_dy(1:lanes, 0) = 0.5_DP*tile_dy(1:lanes, 0)
        DO y = 0, NYS
          DO l = 1, lanes
            tile_dy(l, y) = CWAVX(xb+l-1)*tile_v(l, y, 1)+tile_dy(l, y)+CWAVZ(z)*tile_v(l, y, 3)
          END DO
        END DO
        CALL solve_tau_tile(tile_p, tile_dy, G_pre, a_pre, b_pre, pressure_bc, WN_pre, AD_pre, xb, z, lanes, zero_lane)
        tile_dy(1:lanes, NYS) = 0
        DO l = 1, lanes
          tile_dy(l, NYS-1) = REAL(2*NYS, dp)*tile_p(l, NYS)
        END DO
        DO y = NYS-2, 0, -1
          DO l = 1, lanes
            tile_dy(l, y) = tile_dy(l, y+2)+REAL(2*(y+1), dp)*tile_p(l, y+1)
          END DO
        END DO
        tile_dy(1:lanes, 0) = 0.5_DP*tile_dy(1:lanes, 0)
        DO l = 1, lanes
          tail(l, 1) = tile_v(l, NYS, 2)
          tail(l, 2) = tile_v(l, NYS-1, 2)
        END DO
        DO c = 1, 3
          DO y = 0, NYS
            DO l = 1, lanes
              SELECT CASE (c)
              CASE (1)
                tile_force(l, y) = Re*(CWAVX(xb+l-1)*tile_p(l, y)-tile_v(l, y, c))
              CASE (2)
                tile_force(l, y) = Re*(tile_dy(l, y)-tile_v(l, y, c))
              CASE (3)
                tile_force(l, y) = Re*(CWAVZ(z)*tile_p(l, y)-tile_v(l, y, c))
              END SELECT
            END DO
          END DO
          CALL solve_tau_tile(tile_v(:, :, c), tile_force, G_V, a_V(c), b_V(c), velocity_bc(c, :, :, :), &
                              WN_V, AD_V, xb, z, lanes, zero_lane)
        END DO
        divwall = 0
        DO y = 0, NYS
          DO l = 1, lanes
            tangent = CWAVX(xb+l-1)*tile_v(l, y, 1)+CWAVZ(z)*tile_v(l, y, 3)
            normal = REAL(y*y, dp)*tile_v(l, y, 2)
            divwall(l, 0) = divwall(l, 0)+tangent+normal
            IF (MOD(y, 2) == 0) THEN
              divwall(l, 1) = divwall(l, 1)+tangent-normal
            ELSE
              divwall(l, 1) = divwall(l, 1)-tangent+normal
            END IF
          END DO
        END DO
        DO l = 1, lanes
          x = xb+l-1
          cc(l, 1) = divwall(l, 1)*Div_Vc_BOUND(2, x, 0, z)-divwall(l, 0)*Div_Vc_BOUND(2, x, 1, z)
          cc(l, 2) = divwall(l, 0)*Div_Vc_BOUND(1, x, 1, z)-divwall(l, 1)*Div_Vc_BOUND(1, x, 0, z)
          ptop = tile_p(l, NYS)+(cc(l, 1)*PREc(1, x, NYS, z)+cc(l, 2)*PREc(2, x, NYS, z))
          vtop = tile_v(l, NYS, 2)+(cc(l, 1)*Vc(1, 2, x, NYS, z)+cc(l, 2)*Vc(2, 2, x, NYS, z))
          vprev = tile_v(l, NYS-1, 2)+(cc(l, 1)*Vc(1, 2, x, NYS-1, z)+cc(l, 2)*Vc(2, 2, x, NYS-1, z))
          sigma(l, 0) = (-WN_V(x, z)*vtop-Re*tail(l, 1))/SIGMA0(x, 0, z)
          sigma(l, 1) = (-WN_V(x, z)*vprev+Re*(REAL(2*NYS, dp)*ptop-tail(l, 2)))/SIGMA0(x, 1, z)
        END DO
        DO y = 0, NYS
          q = MOD(y+1, 2)
          DO l = 1, lanes
            x = xb+l-1
            pressure(x, y, z) = tile_p(l, y)+(cc(l, 1)*PREc(1, x, y, z)+cc(l, 2)*PREc(2, x, y, z)) &
                                +sigma(l, q)*PRE0(x, y, z)
            velocity(1, x, y, z) = tile_v(l, y, 1)+(cc(l, 1)*Vc(1, 1, x, y, z)+cc(l, 2)*Vc(2, 1, x, y, z)) &
                                   +sigma(l, q)*V0(1, x, y, z)
            velocity(2, x, y, z) = tile_v(l, y, 2)+(cc(l, 1)*Vc(1, 2, x, y, z)+cc(l, 2)*Vc(2, 2, x, y, z)) &
                                   +sigma(l, 1-q)*V0(2, x, y, z)
            velocity(3, x, y, z) = tile_v(l, y, 3)+(cc(l, 1)*Vc(1, 3, x, y, z)+cc(l, 2)*Vc(2, 3, x, y, z)) &
                                   +sigma(l, q)*V0(3, x, y, z)
          END DO
        END DO
      END DO
    END DO
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE EQSOLVE_PV_setup
    USE MESH
    USE EQSOLVE_PV

    IMPLICIT NONE

    CALL WN_AD(WN_pre, AD_pre, G_pre)
    CALL WN_AD(WN_V, AD_V, G_V)
    CALL Infulence_coefficient
    CALL B_PROBLEM
    IF (INDEX_MPI0 == 0) CALL EQSOLVE_UW_DPDXZ_1
    RETURN
  END

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE WN_AD(WN, AD, G)
    USE FFT
    USE MESH

    IMPLICIT NONE
    INTEGER :: ix, iy, iz
    REAL(dp) :: G, FC

    REAL(dp) WN(0:NXSP, -NZS:NZS), AD(0:NXSP, 2:NYS, -NZS:NZS)

    DO iz = -NZS, NZS
      DO ix = 0, NXSP
        WN(ix, iz) = -(WAVX(ix)**2+WAVZ(iz)**2+G)
      END DO
    END DO

    AD(0:NXSP, NYS, -NZS:NZS) = 1.0D0/(1.0D0+WD(NYS)*WN)
    AD(0:NXSP, NYS-1, -NZS:NZS) = 1.0D0/(1.0D0+WD(NYS-1)*WN)

    DO iy = NYS-2, 2, -1
      FC = WR(iy)*WL(iy+2)
      AD(0:NXSP, iy, -NZS:NZS) = 1.0D0/(1.0D0+(WD(iy)-FC*WN*AD(0:NXSP, iy+2, -NZS:NZS))*WN)
    END DO

    RETURN
  END

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE Infulence_coefficient
    USE MESH
    USE EQSOLVE_PV

    IMPLICIT NONE
    INTEGER :: I, J

    COMPLEX(dp), ALLOCATABLE :: RHS_PRE(:, :, :), BC_PRE(:, :, :)
    COMPLEX(dp), ALLOCATABLE :: RHS_V(:, :, :, :), BC_V(:, :, :, :)

    ALLOCATE (RHS_PRE(0:NXSP, 0:NYS, -NZS:NZS), BC_PRE(0:NXSP, 0:1, -NZS:NZS))
    ALLOCATE (RHS_V(3, 0:NXSP, 0:NYS, -NZS:NZS), BC_V(3, 0:NXSP, 0:1, -NZS:NZS))

    DO I = 1, 2
      RHS_PRE = CMPLX(0.0D0, 0.0D0, kind=dp)
      BC_PRE(0:NXSP, 0, -NZS:NZS) = CMPLX(REAL(2-I, kind=dp), 0.0D0, kind=dp)
      BC_PRE(0:NXSP, 1, -NZS:NZS) = CMPLX(REAL(I-1, kind=dp), 0.0D0, kind=dp)

      RHS_V = CMPLX(0.0D0, 0.0D0, kind=dp)
      BC_V = CMPLX(0.0D0, 0.0D0, kind=dp)

      CALL EQSOLVE_PV_1(PREc(I, 0:NXSP, 0:NYS, -NZS:NZS), RHS_PRE, BC_PRE, Vc(I, 1:3, 0:NXSP, 0:NYS, -NZS:NZS), RHS_V, BC_V)

      CALL Div_BOUND(Div_Vc_BOUND(I, 0:NXSP, 0:1, -NZS:NZS), Vc(I, 1:3, 0:NXSP, 0:NYS, -NZS:NZS))
    END DO

    BC_PRE(0:NXSP, 0, -NZS:NZS) = (Div_Vc_BOUND(1, 0:NXSP, 0, -NZS:NZS)*Div_Vc_BOUND(2, 0:NXSP, 1, -NZS:NZS) &
                                   -Div_Vc_BOUND(1, 0:NXSP, 1, -NZS:NZS)*Div_Vc_BOUND(2, 0:NXSP, 0, -NZS:NZS))
    IF (INDEX_MPI0 == 0) THEN
      IF (a_V(2) == 0.0D0) THEN
        BC_PRE(0, 0, 0) = 1.0D0
      ELSE
        BC_PRE(0, 0, 0) = 0.25D0/Div_Vc_BOUND(1, 0, 0, 0)
      END IF
    END IF

    DO J = 0, 1
      DO I = 1, 2
        Div_Vc_BOUND(I, 0:NXSP, J, -NZS:NZS) = Div_Vc_BOUND(I, 0:NXSP, J, -NZS:NZS)/BC_PRE(0:NXSP, 0, -NZS:NZS)
        ! At kx=kz=0, the two wall divergence constraints are dependent.
        ! Apply one independent wall condition to the zero mode.
        IF (INDEX_MPI0 == 0) THEN
          Div_Vc_BOUND(I, 0, J, 0) = 0
          IF (J == 0) Div_Vc_BOUND(I, 0, J, 0) = BC_PRE(0, 0, 0)
        END IF
      END DO
    END DO

    DEALLOCATE (RHS_PRE, BC_PRE, RHS_V, BC_V)

    RETURN
  END

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE B_PROBLEM
    USE MESH
    USE EQSOLVE_PV
    USE dns_config, ONLY: Re

    IMPLICIT NONE
    INTEGER :: iy, INDEX_iy

    COMPLEX(dp), ALLOCATABLE :: RHS_PRE(:, :, :), BC_PRE(:, :, :)
    COMPLEX(dp), ALLOCATABLE :: RHS_V(:, :, :, :), BC_V(:, :, :, :)

    ALLOCATE (RHS_PRE(0:NXSP, 0:NYS, -NZS:NZS), BC_PRE(0:NXSP, 0:1, -NZS:NZS))
    ALLOCATE (RHS_V(3, 0:NXSP, 0:NYS, -NZS:NZS), BC_V(3, 0:NXSP, 0:1, -NZS:NZS))

    RHS_PRE(0:NXSP, 0, -NZS:NZS) = CMPLX(REAL(NYS-1, kind=dp), 0.0D0, kind=dp)
    DO iy = 1, NYS-1
      INDEX_iy = MOD(iy+1, 2)
      RHS_PRE(0:NXSP, iy, -NZS:NZS) = CMPLX(REAL(2*(NYS-INDEX_iy), kind=dp), 0.0D0, kind=dp)
    END DO
    RHS_PRE(0:NXSP, NYS, -NZS:NZS) = CMPLX(0.0D0, 0.0D0, kind=dp)

    BC_PRE = CMPLX(0.0D0, 0.0D0, kind=dp)

    RHS_V = CMPLX(0.0D0, 0.0D0, kind=dp)
    BC_V = CMPLX(0.0D0, 0.0D0, kind=dp)

    CALL EQSOLVE_PV_2(PRE0, RHS_PRE, BC_PRE, V0, RHS_V, BC_V)

    SIGMA0(0:NXSP, 0, -NZS:NZS) = -WN_V*V0(2, 0:NXSP, NYS, -NZS:NZS)
    SIGMA0(0:NXSP, 1, -NZS:NZS) = -WN_V*V0(2, 0:NXSP, NYS-1, -NZS:NZS)+Re*REAL(NYS+NYS, kind=dp)*PRE0(0:NXSP, NYS, -NZS:NZS)

    SIGMA0 = Re-SIGMA0

    IF (INDEX_MPI0 == 0) SIGMA0(0, 1, 0) = 1.0D0

    DEALLOCATE (RHS_PRE, BC_PRE, RHS_V, BC_V)

    RETURN
  END

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE EQSOLVE_UW_DPDXZ_1

    USE MESH
    USE EQSOLVE_PV
    USE dns_config, ONLY: Re

    IMPLICIT NONE
    INTEGER :: iy

    COMPLEX(dp), ALLOCATABLE :: V(:, :, :), RHS_V(:, :, :), BC_V(:, :, :)
    ALLOCATE (V(0:NXSP, 0:NYS, -NZS:NZS), RHS_V(0:NXSP, 0:NYS, -NZS:NZS), BC_V(0:NXSP, 0:1, -NZS:NZS))
    RHS_V = CMPLX(0.0D0, 0.0D0, kind=dp)
    RHS_V(0, 0, 0) = CMPLX(Re, 0.0D0, kind=dp)

    BC_V = CMPLX(0.0D0, 0.0D0, kind=dp)

    CALL SOLVE(V, RHS_V, G_V, a_V(1), b_V(1), BC_V, WN_V, AD_V)

    CALL FOR_Q(QX_DPDX_1, V)

    DO iy = 0, NYS
      U_DPDX_1(iy) = REAL(V(0, iy, 0), kind=dp)
    END DO
    RHS_V = CMPLX(0.0D0, 0.0D0, kind=dp)
    RHS_V(0, 0, 0) = CMPLX(Re, 0.0D0, kind=dp)

    BC_V = CMPLX(0.0D0, 0.0D0, kind=dp)

    CALL SOLVE(V, RHS_V, G_V, a_V(3), b_V(3), BC_V, WN_V, AD_V)

    CALL FOR_Q(QZ_DPDZ_1, V)

    DO iy = 0, NYS
      W_DPDZ_1(iy) = REAL(V(0, iy, 0), kind=dp)
    END DO

    DEALLOCATE (V, RHS_V, BC_V)

    RETURN
  END

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE EQSOLVE_PV_1(PRE, RHS_PRE, BC_PRE, V, RHS_V, BC_V)

    USE MESH
    USE EQSOLVE_PV
    USE dns_config, ONLY: Re

    IMPLICIT NONE
    INTEGER :: I

    COMPLEX(dp) PRE(0:NXSP, 0:NYS, -NZS:NZS), RHS_PRE(0:NXSP, 0:NYS, -NZS:NZS), BC_PRE(0:NXSP, 0:1, -NZS:NZS)
    COMPLEX(dp) V(3, 0:NXSP, 0:NYS, -NZS:NZS), RHS_V(3, 0:NXSP, 0:NYS, -NZS:NZS), BC_V(3, 0:NXSP, 0:1, -NZS:NZS)

    COMPLEX(dp), ALLOCATABLE :: TEMP(:, :, :, :)

    ALLOCATE (TEMP(3, 0:NXSP, 0:NYS, -NZS:NZS))

    CALL SOLVE(PRE, RHS_PRE, G_pre, a_pre, b_pre, BC_PRE, WN_pre, AD_pre)

    CALL Grad(TEMP, PRE)

    TEMP = Re*(TEMP-RHS_V)

    DO I = 1, 3
      CALL SOLVE(V(I, 0:NXSP, 0:NYS, -NZS:NZS), TEMP(I, 0:NXSP, 0:NYS, -NZS:NZS), &
                 G_V, a_V(I), b_V(I), BC_V(I, 0:NXSP, 0:1, -NZS:NZS), WN_V, AD_V)
    END DO

    DEALLOCATE (TEMP)

    RETURN
  END

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE EQSOLVE_PV_2(PRE, RHS_PRE, BC_PRE, V, RHS_V, BC_V)

    USE MESH
    USE EQSOLVE_PV
    USE dns_config, ONLY: Re

    IMPLICIT NONE
    INTEGER :: I, iy

    COMPLEX(dp) PRE(0:NXSP, 0:NYS, -NZS:NZS), RHS_PRE(0:NXSP, 0:NYS, -NZS:NZS), BC_PRE(0:NXSP, 0:1, -NZS:NZS)
    COMPLEX(dp) V(3, 0:NXSP, 0:NYS, -NZS:NZS), RHS_V(3, 0:NXSP, 0:NYS, -NZS:NZS), BC_V(3, 0:NXSP, 0:1, -NZS:NZS)

    COMPLEX(dp), ALLOCATABLE :: TEMP(:, :, :, :)
    COMPLEX(dp), ALLOCATABLE :: Div_V_BOUND(:, :, :)

    ALLOCATE (TEMP(3, 0:NXSP, 0:NYS, -NZS:NZS), Div_V_BOUND(0:NXSP, 0:1, -NZS:NZS))

    CALL SOLVE(PRE, RHS_PRE, G_pre, a_pre, b_pre, BC_PRE, WN_pre, AD_pre)

    CALL Grad(TEMP, PRE)

    TEMP = Re*(TEMP-RHS_V)

    DO I = 1, 3
      CALL SOLVE(V(I, 0:NXSP, 0:NYS, -NZS:NZS), TEMP(I, 0:NXSP, 0:NYS, -NZS:NZS), &
                 G_V, a_V(I), b_V(I), BC_V(I, 0:NXSP, 0:1, -NZS:NZS), WN_V, AD_V)
    END DO

    CALL Div_BOUND(Div_V_BOUND, V)

    TEMP(1, 0:NXSP, 0, -NZS:NZS) = Div_V_BOUND(0:NXSP, 1, -NZS:NZS)*Div_Vc_BOUND(2, 0:NXSP, 0, -NZS:NZS) &
                                   -Div_V_BOUND(0:NXSP, 0, -NZS:NZS)*Div_Vc_BOUND(2, 0:NXSP, 1, -NZS:NZS)
    TEMP(2, 0:NXSP, 0, -NZS:NZS) = Div_V_BOUND(0:NXSP, 0, -NZS:NZS)*Div_Vc_BOUND(1, 0:NXSP, 1, -NZS:NZS) &
                                   -Div_V_BOUND(0:NXSP, 1, -NZS:NZS)*Div_Vc_BOUND(1, 0:NXSP, 0, -NZS:NZS)

    DO iy = 0, NYS
      DO I = 1, 3
        V(I, 0:NXSP, iy, -NZS:NZS) = V(I, 0:NXSP, iy, -NZS:NZS)+TEMP(1, 0:NXSP, 0, -NZS:NZS)*Vc(1, I, 0:NXSP, iy, -NZS:NZS) &
                                     +TEMP(2, 0:NXSP, 0, -NZS:NZS)*Vc(2, I, 0:NXSP, iy, -NZS:NZS)
      END DO
      PRE(0:NXSP, iy, -NZS:NZS) = PRE(0:NXSP, iy, -NZS:NZS)+TEMP(1, 0:NXSP, 0, -NZS:NZS)*PREc(1, 0:NXSP, iy, -NZS:NZS) &
                                  +TEMP(2, 0:NXSP, 0, -NZS:NZS)*PREc(2, 0:NXSP, iy, -NZS:NZS)
    END DO

    DEALLOCATE (TEMP, Div_V_BOUND)

    RETURN
  END

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE SOLVE(U, RHS, G, a, b, BC, WN, AD)

    USE FFT
    USE MESH

    IMPLICIT NONE
    INTEGER :: iy, index
    REAL(dp) :: G, a, b, FUNC

    COMPLEX(dp), INTENT(out) :: U(0:NXSP, 0:NYS, -NZS:NZS)
    COMPLEX(dp), INTENT(in) :: RHS(0:NXSP, 0:NYS, -NZS:NZS)
    COMPLEX(dp) :: F(0:NXSP, 0:NYS, -NZS:NZS)
    COMPLEX(dp), DIMENSION(0:NXSP, 0:1, -NZS:NZS) :: BC

    REAL(dp) WN(0:NXSP, -NZS:NZS), AD(0:NXSP, 2:NYS, -NZS:NZS)

    COMPLEX(dp) COR

    REAL(dp), ALLOCATABLE :: AD01(:, :, :)
    ALLOCATE (AD01(0:NXSP, 0:1, -NZS:NZS))

    ! The singular pressure compatibility condition changes only a local RHS copy.
    F = RHS
    IF (G == 0.0D0 .AND. a == 0.0D0 .AND. INDEX_MPI0 == 0) THEN
      COR = CMPLX(0.0D0, 0.0D0, kind=dp)
      DO iy = 2, NYS, 2
        COR = COR+F(0, iy, 0)/REAL(iy**2-1, kind=dp)
      END DO
      IF (b == 0.0D0) THEN
        F(0, 0, 0) = COR
      ELSE
        COR = COR+(BC(0, 0, 0)+BC(0, 1, 0))/2.0D0/b
        F(0, 0, 0) = COR
      END IF
    END IF

    U(0:NXSP, 0, -NZS:NZS) = (BC(0:NXSP, 0, -NZS:NZS)+BC(0:NXSP, 1, -NZS:NZS))/2.0D0
    U(0:NXSP, 1, -NZS:NZS) = (BC(0:NXSP, 0, -NZS:NZS)-BC(0:NXSP, 1, -NZS:NZS))/2.0D0
    DO iy = 2, NYS-2
      U(0:NXSP, iy, -NZS:NZS) = WL(iy)*F(0:NXSP, iy-2, -NZS:NZS) &
                                +WD(iy)*F(0:NXSP, iy, -NZS:NZS) &
                                +WR(iy)*F(0:NXSP, iy+2, -NZS:NZS)
    END DO
    U(0:NXSP, NYS-1, -NZS:NZS) = WL(NYS-1)*F(0:NXSP, NYS-3, -NZS:NZS) &
                                 +WD(NYS-1)*F(0:NXSP, NYS-1, -NZS:NZS)
    U(0:NXSP, NYS, -NZS:NZS) = WL(NYS)*F(0:NXSP, NYS-2, -NZS:NZS) &
                               +WD(NYS)*F(0:NXSP, NYS, -NZS:NZS)

    DO iy = NYS-2, 2, -1
      U(0:NXSP, iy, -NZS:NZS) = U(0:NXSP, iy, -NZS:NZS)-WR(iy)*WN &
                                *AD(0:NXSP, iy+2, -NZS:NZS)*U(0:NXSP, iy+2, -NZS:NZS)
    END DO

    DO iy = NYS, NYS-1, -1
      index = MOD(iy, 2)
      IF (a == 0.0D0 .AND. b == 0.0D0) THEN
        FUNC = REAL((1-index)*iy**2+index, kind=dp)
      ELSE
        FUNC = a+b*REAL(iy**2, kind=dp)
      END IF
      AD01(0:NXSP, index, -NZS:NZS) = FUNC
    END DO

    DO iy = NYS-2, 0, -1
      index = MOD(iy, 2)
      IF (a == 0.0D0 .AND. b == 0.0D0) THEN
        FUNC = REAL((1-index)*iy**2+index, kind=dp)
      ELSE
        FUNC = a+b*REAL(iy**2, kind=dp)
      END IF
      U(0:NXSP, index, -NZS:NZS) = U(0:NXSP, index, -NZS:NZS)-AD01(0:NXSP, index, -NZS:NZS) &
                                   *AD(0:NXSP, iy+2, -NZS:NZS)*U(0:NXSP, iy+2, -NZS:NZS)
      AD01(0:NXSP, index, -NZS:NZS) = FUNC-WL(iy+2)*WN &
                                      *AD01(0:NXSP, index, -NZS:NZS)*AD(0:NXSP, iy+2, -NZS:NZS)
    END DO

    IF (G == 0.0D0 .AND. a == 0.0D0 .AND. INDEX_MPI0 == 0) THEN
      U(0, 0, 0) = CMPLX(0.0D0, 0.0D0, kind=dp)
      AD01(0, 0, 0) = 1.0D0
    END IF

    U(0:NXSP, 0, -NZS:NZS) = U(0:NXSP, 0, -NZS:NZS)/AD01(0:NXSP, 0, -NZS:NZS)
    U(0:NXSP, 1, -NZS:NZS) = U(0:NXSP, 1, -NZS:NZS)/AD01(0:NXSP, 1, -NZS:NZS)

    DO iy = 2, NYS
      U(0:NXSP, iy, -NZS:NZS) = (U(0:NXSP, iy, -NZS:NZS)-WL(iy)*WN &
                                 *U(0:NXSP, iy-2, -NZS:NZS))*AD(0:NXSP, iy, -NZS:NZS)
    END DO

    DEALLOCATE (AD01)

    RETURN
  END

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE SOLVE2(U, F, G, a, b, BC, WN, AD)

    USE FFT
    USE MESH

    IMPLICIT NONE
    INTEGER :: iy, index, ikx, ikz
    REAL(dp) :: G, DFLOAT_iy2, FUNC_AD, FUNC_BD, CTEMP

    COMPLEX(dp), DIMENSION(0:NXSP, 0:NYS, -NZS:NZS) :: U, F
    COMPLEX(dp), DIMENSION(0:NXSP, 0:1, -NZS:NZS) :: BC

    REAL(dp) WN(0:NXSP, -NZS:NZS), AD(0:NXSP, 2:NYS, -NZS:NZS), a(0:1), b(0:1)

    COMPLEX(dp) C0, C1

    REAL(dp), ALLOCATABLE :: AD01(:, :, :), BD01(:, :, :)
    ALLOCATE (AD01(0:NXSP, 0:1, -NZS:NZS), BD01(0:NXSP, 0:1, -NZS:NZS))

    U(0:NXSP, 0:1, -NZS:NZS) = BC(0:NXSP, 0:1, -NZS:NZS)
    DO iy = 2, NYS-2
      U(0:NXSP, iy, -NZS:NZS) = WL(iy)*F(0:NXSP, iy-2, -NZS:NZS) &
                                +WD(iy)*F(0:NXSP, iy, -NZS:NZS) &
                                +WR(iy)*F(0:NXSP, iy+2, -NZS:NZS)
    END DO
    U(0:NXSP, NYS-1, -NZS:NZS) = WL(NYS-1)*F(0:NXSP, NYS-3, -NZS:NZS) &
                                 +WD(NYS-1)*F(0:NXSP, NYS-1, -NZS:NZS)
    U(0:NXSP, NYS, -NZS:NZS) = WL(NYS)*F(0:NXSP, NYS-2, -NZS:NZS) &
                               +WD(NYS)*F(0:NXSP, NYS, -NZS:NZS)

    DO iy = NYS-2, 2, -1
      U(0:NXSP, iy, -NZS:NZS) = U(0:NXSP, iy, -NZS:NZS)-WR(iy)*WN &
                                *AD(0:NXSP, iy+2, -NZS:NZS)*U(0:NXSP, iy+2, -NZS:NZS)
    END DO

    DO iy = NYS, NYS-1, -1
      index = MOD(iy, 2)
      DFLOAT_iy2 = REAL(iy**2, kind=dp)
      FUNC_AD = REAL(1-index-index, kind=dp)*(a(index)+b(index)*DFLOAT_iy2)
      FUNC_BD = a(1-index)+b(1-index)*DFLOAT_iy2
      AD01(0:NXSP, index, -NZS:NZS) = FUNC_AD
      BD01(0:NXSP, index, -NZS:NZS) = FUNC_BD
    END DO

    DO iy = NYS-2, 0, -1
      index = MOD(iy, 2)
      DFLOAT_iy2 = REAL(iy**2, kind=dp)
      FUNC_AD = REAL(1-index-index, kind=dp)*(a(index)+b(index)*DFLOAT_iy2)
      FUNC_BD = a(1-index)+b(1-index)*DFLOAT_iy2
      U(0:NXSP, index, -NZS:NZS) = U(0:NXSP, index, -NZS:NZS)-AD01(0:NXSP, index, -NZS:NZS) &
                                   *AD(0:NXSP, iy+2, -NZS:NZS)*U(0:NXSP, iy+2, -NZS:NZS)
      U(0:NXSP, 1-index, -NZS:NZS) = U(0:NXSP, 1-index, -NZS:NZS)-BD01(0:NXSP, index, -NZS:NZS) &
                                     *AD(0:NXSP, iy+2, -NZS:NZS)*U(0:NXSP, iy+2, -NZS:NZS)
      AD01(0:NXSP, index, -NZS:NZS) = FUNC_AD-WL(iy+2)*WN &
                                      *AD01(0:NXSP, index, -NZS:NZS)*AD(0:NXSP, iy+2, -NZS:NZS)
      BD01(0:NXSP, index, -NZS:NZS) = FUNC_BD-WL(iy+2)*WN &
                                      *BD01(0:NXSP, index, -NZS:NZS)*AD(0:NXSP, iy+2, -NZS:NZS)
    END DO

    DO ikx = 0, NXSP
      DO ikz = -NZS, NZS
        CTEMP = AD01(ikx, 0, ikz)*AD01(ikx, 1, ikz)-BD01(ikx, 0, ikz)*BD01(ikx, 1, ikz)
        C0 = (U(ikx, 0, ikz)*AD01(ikx, 1, ikz)-U(ikx, 1, ikz)*BD01(ikx, 1, ikz))/CTEMP
        C1 = (AD01(ikx, 0, ikz)*U(ikx, 1, ikz)-BD01(ikx, 0, ikz)*U(ikx, 0, ikz))/CTEMP
        U(ikx, 0, ikz) = C0
        U(ikx, 1, ikz) = C1
      END DO
    END DO

    DO iy = 2, NYS
      U(0:NXSP, iy, -NZS:NZS) = (U(0:NXSP, iy, -NZS:NZS)-WL(iy)*WN &
                                 *U(0:NXSP, iy-2, -NZS:NZS))*AD(0:NXSP, iy, -NZS:NZS)
    END DO

    DEALLOCATE (AD01, BD01)

    RETURN
  END
END MODULE
! ---&---1---------2---------3---------4---------5---------6---------7--
