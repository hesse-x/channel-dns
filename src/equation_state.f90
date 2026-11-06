! ---&---1---------2---------3---------4---------5---------6---------7--
! Helmholtz factors, influence functions and zero-mode responses.
MODULE eqsolve_pv
  USE dns_kinds
  USE mesh
  USE dns_config, ONLY: G_pre, a_V, b_V
  IMPLICIT NONE
  REAL(dp) :: a_pre, b_pre, G_V, QX_DPDX_1, QZ_DPDZ_1
  REAL(dp), ALLOCATABLE :: WN_pre(:, :), AD_pre(:, :, :), WN_V(:, :), AD_V(:, :, :), U_DPDX_1(:), W_DPDZ_1(:)
  COMPLEX(dp), ALLOCATABLE :: PREc(:, :, :, :), Vc(:, :, :, :, :), Div_Vc_BOUND(:, :, :, :)
  COMPLEX(dp), ALLOCATABLE :: PRE0(:, :, :), V0(:, :, :, :), SIGMA0(:, :, :)
CONTAINS

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE allocate_equations()
    ALLOCATE (WN_pre(0:NXSP, -NZS:NZS), WN_V(0:NXSP, -NZS:NZS))
    ALLOCATE (AD_pre(0:NXSP, 2:NYS, -NZS:NZS), AD_V(0:NXSP, 2:NYS, -NZS:NZS))
    ALLOCATE (PREc(2, 0:NXSP, 0:NYS, -NZS:NZS), Vc(2, 3, 0:NXSP, 0:NYS, -NZS:NZS))
    ALLOCATE (Div_Vc_BOUND(2, 0:NXSP, 0:1, -NZS:NZS), U_DPDX_1(0:NYS), W_DPDZ_1(0:NYS))
    ALLOCATE (PRE0(0:NXSP, 0:NYS, -NZS:NZS), V0(3, 0:NXSP, 0:NYS, -NZS:NZS), SIGMA0(0:NXSP, 0:1, -NZS:NZS))
  END SUBROUTINE
END MODULE
! ---&---1---------2---------3---------4---------5---------6---------7--
