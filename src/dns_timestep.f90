! ---&---1---------2---------3---------4---------5---------6---------7--
! Rotational nonlinear term and coupled pressure-velocity step.
MODULE dns_timestep
  USE dns_kinds, ONLY: dp
  USE mesh
  USE fft, ONLY: CWAVXRe, CWAVZRe
  USE dns_vector, ONLY: Curl, CROSSTIME, prepare_vector_workspace, close_vector_workspace
  USE dns_equations, ONLY: solve_modes_fused, prepare_fused_workspace, close_fused_workspace
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: NLTERM, STEP23, prepare_timestep_workspace, close_timestep_workspace
  COMPLEX(dp), ALLOCATABLE, SAVE :: curl_work(:, :, :, :)
CONTAINS

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE prepare_timestep_workspace()
    IF (.NOT. ALLOCATED(curl_work)) ALLOCATE (curl_work(3, 0:NXSP, 0:NYS, -NZS:NZS))
    CALL prepare_vector_workspace()
    CALL prepare_fused_workspace()
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE close_timestep_workspace()
    IF (ALLOCATED(curl_work)) DEALLOCATE (curl_work)
    CALL close_vector_workspace()
    CALL close_fused_workspace()
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE NLTERM(nl_v, v, bc_pre)
    COMPLEX(dp), INTENT(out) :: nl_v(3, 0:NXSP, 0:NYS, -NZS:NZS), bc_pre(0:NXSP, 0:1, -NZS:NZS)
    COMPLEX(dp), INTENT(in) :: v(3, 0:NXSP, 0:NYS, -NZS:NZS)
    COMPLEX(dp) :: wall_term
    INTEGER :: x, y, z

    CALL prepare_timestep_workspace()
! W = rot(V); NL_V = V x W.
    CALL Curl(curl_work, v)
    CALL CROSSTIME(nl_v, v, curl_work)
! Pressure Neumann data at y = +1 and y = -1.
    bc_pre = 0
    DO z = -NZS, NZS
      DO y = 0, NYS
        DO x = 0, NXSP
          wall_term = CWAVXRe(x)*curl_work(3, x, y, z)-CWAVZRe(z)*curl_work(1, x, y, z)+nl_v(2, x, y, z)
          bc_pre(x, 0, z) = bc_pre(x, 0, z)+wall_term
          IF (MOD(y, 2) == 0) THEN
            bc_pre(x, 1, z) = bc_pre(x, 1, z)-wall_term
          ELSE
            bc_pre(x, 1, z) = bc_pre(x, 1, z)+wall_term
          END IF
        END DO
      END DO
    END DO
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE STEP23(pre, bc_pre, v, rhs_v, bc_v)
    COMPLEX(dp), INTENT(out) :: pre(0:NXSP, 0:NYS, -NZS:NZS), v(3, 0:NXSP, 0:NYS, -NZS:NZS)
    COMPLEX(dp), INTENT(in) :: bc_pre(0:NXSP, 0:1, -NZS:NZS), &
                               rhs_v(3, 0:NXSP, 0:NYS, -NZS:NZS), bc_v(3, 0:NXSP, 0:1, -NZS:NZS)

    CALL solve_modes_fused(pre, bc_pre, v, rhs_v, bc_v)
! The zero-mode normal velocity is constant across the channel.
    IF (INDEX_MPI0 == 0) THEN
      v(2, 0, 0, 0) = 0.5_DP*(bc_v(2, 0, 0, 0)+bc_v(2, 0, 1, 0))
      v(2, 0, 1:NYS, 0) = 0
    END IF
  END SUBROUTINE
END MODULE
! ---&---1---------2---------3---------4---------5---------6---------7--
