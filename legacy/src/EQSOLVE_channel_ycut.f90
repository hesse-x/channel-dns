! ---&---1---------2---------3---------4---------5---------6---------7--
!
SUBROUTINE EQSOLVE_PV_setup
  USE MESH
  USE EQSOLVE_PV
  IMPLICIT REAL*8 (A-H, O-Z)
  CALL WN_AD(WN_pre,AD_pre,G_pre)
  CALL WN_AD(WN_V,AD_V,G_V)
  CALL Infulence_coefficient
  CALL B_PROBLEM
  IF(INDEX_MPI0==0) CALL EQSOLVE_UW_DPDXZ_1
  RETURN
END
! ---&---1---------2---------3---------4---------5---------6---------7--
!
SUBROUTINE WN_AD(WN,AD,G)
  USE FFT
  USE MESH
  IMPLICIT REAL*8 (A-H, O-Z)

  REAL*8 WN(0:NXSP,-NZS:NZS),AD(0:NXSP,2:NYS,-NZS:NZS)

  DO iz = -NZS,NZS
    DO ix = 0,NXSP
      WN(ix,iz) = -(WAVX(ix)**2 + WAVZ(iz)**2 + G)
    ENDDO
  ENDDO

  AD(0:NXSP,NYS,  -NZS:NZS) = 1.0D0/(1.0D0+ WD(NYS  )*WN)
  AD(0:NXSP,NYS-1,-NZS:NZS) = 1.0D0/(1.0D0+ WD(NYS-1)*WN)

  DO iy = NYS-2,2,-1
    FC = WR(iy)*WL(iy+2)
    AD(0:NXSP,iy,-NZS:NZS) = 1.0D0/(1.0D0+(WD(iy)-FC*WN*AD(0:NXSP,iy+2,-NZS:NZS))*WN)
  ENDDO

  RETURN
END
! ---&---1---------2---------3---------4---------5---------6---------7--
!
SUBROUTINE Infulence_coefficient
  USE MESH
  USE EQSOLVE_PV
  IMPLICIT REAL*8 (A-H,O-Z)

  COMPLEX*16, ALLOCATABLE :: RHS_PRE(:,:,:),BC_PRE(:,:,:)
  COMPLEX*16, ALLOCATABLE :: RHS_V(:,:,:,:),BC_V(:,:,:,:)

  ALLOCATE( RHS_PRE(0:NXSP,0:NYS,-NZS:NZS),BC_PRE(0:NXSP,0:1,-NZS:NZS) )
  ALLOCATE( RHS_V(3,0:NXSP,0:NYS,-NZS:NZS),BC_V(3,0:NXSP,0:1,-NZS:NZS) )

  DO I=1,2
    RHS_PRE = DCMPLX(0.0D0, 0.0D0)
    BC_PRE(0:NXSP,0,-NZS:NZS) = DCMPLX(DFLOAT(2-I), 0.0D0)
    BC_PRE(0:NXSP,1,-NZS:NZS) = DCMPLX(DFLOAT(I-1), 0.0D0)

    RHS_V = DCMPLX(0.0D0, 0.0D0)
    BC_V = DCMPLX(0.0D0, 0.0D0)

    CALL EQSOLVE_PV_1 (PREc(I,0:NXSP,0:NYS,-NZS:NZS),RHS_PRE,BC_PRE, Vc(I,1:3,0:NXSP,0:NYS,-NZS:NZS),RHS_V,BC_V)

    CALL Div_BOUND(Div_Vc_BOUND(I,0:NXSP,0:1,-NZS:NZS),Vc(I,1:3,0:NXSP,0:NYS,-NZS:NZS))
  ENDDO

  BC_PRE(0:NXSP,0,-NZS:NZS) = (Div_Vc_BOUND(1,0:NXSP,0,-NZS:NZS)*Div_Vc_BOUND(2,0:NXSP,1,-NZS:NZS) &
       - Div_Vc_BOUND(1,0:NXSP,1,-NZS:NZS)*Div_Vc_BOUND(2,0:NXSP,0,-NZS:NZS))
  IF(INDEX_MPI0==0) THEN
    IF(a_V(2)==0.0D0) THEN
      BC_PRE(0,0,0) = 1.0D0
    ELSE
      BC_PRE(0,0,0) = 0.25D0/Div_Vc_BOUND(1,0,0,0)
    ENDIF
  ENDIF

  DO J=0,1
    DO I=1,2
      Div_Vc_BOUND(I,0:NXSP,J,-NZS:NZS)=Div_Vc_BOUND(I,0:NXSP,J,-NZS:NZS)/BC_PRE(0:NXSP,0,-NZS:NZS)
      IF(INDEX_MPI0==0) Div_Vc_BOUND(I,0,J,0)=BC_PRE(0,0,0)
    ENDDO
  ENDDO

  DEALLOCATE(RHS_PRE,BC_PRE, RHS_V,BC_V)

  RETURN
END
! ---&---1---------2---------3---------4---------5---------6---------7--
!
SUBROUTINE B_PROBLEM
USE MESH
USE EQSOLVE_PV
INCLUDE 'com_channel.h'

COMPLEX*16, ALLOCATABLE :: RHS_PRE(:,:,:),BC_PRE(:,:,:)
COMPLEX*16, ALLOCATABLE :: RHS_V(:,:,:,:),BC_V(:,:,:,:)

ALLOCATE(RHS_PRE(0:NXSP,0:NYS,-NZS:NZS),BC_PRE(0:NXSP,0:1,-NZS:NZS))
ALLOCATE(RHS_V(3,0:NXSP,0:NYS,-NZS:NZS),BC_V(3,0:NXSP,0:1,-NZS:NZS))

RHS_PRE(0:NXSP,0,-NZS:NZS) = DCMPLX(DFLOAT(NYS-1), 0.0D0)
DO iy=1,NYS-1
INDEX_iy=MOD(iy+1,2)
RHS_PRE(0:NXSP,iy,-NZS:NZS) = DCMPLX(DFLOAT(2*(NYS-INDEX_iy)), 0.0D0)
ENDDO
RHS_PRE(0:NXSP,NYS,-NZS:NZS) = DCMPLX(0.0D0, 0.0D0)

BC_PRE = DCMPLX(0.0D0, 0.0D0)

RHS_V = DCMPLX(0.0D0, 0.0D0)
BC_V = DCMPLX(0.0D0, 0.0D0)

CALL EQSOLVE_PV_2 (PRE0,RHS_PRE,BC_PRE, V0,RHS_V,BC_V)

SIGMA0(0:NXSP,0,-NZS:NZS)=-WN_V*V0(2,0:NXSP,NYS  ,-NZS:NZS)
SIGMA0(0:NXSP,1,-NZS:NZS)=-WN_V*V0(2,0:NXSP,NYS-1,-NZS:NZS) + Re*DFLOAT(NYS+NYS)*PRE0(0:NXSP,NYS,-NZS:NZS)

SIGMA0=Re-SIGMA0

IF(INDEX_MPI0==0) SIGMA0(0,1,0)=1.0D0

DEALLOCATE(RHS_PRE,BC_PRE, RHS_V,BC_V)

RETURN
END
! ---&---1---------2---------3---------4---------5---------6---------7--
!
 SUBROUTINE EQSOLVE_UW_DPDXZ_1

 USE MESH
  USE EQSOLVE_PV
 INCLUDE 'com_channel.h'

 COMPLEX*16, ALLOCATABLE :: V(:,:,:),RHS_V(:,:,:),BC_V(:,:,:)
 ALLOCATE(V(0:NXSP,0:NYS,-NZS:NZS),RHS_V(0:NXSP,0:NYS,-NZS:NZS),BC_V(0:NXSP,0:1,-NZS:NZS))
!FOR U_DPDX_1
 RHS_V = DCMPLX(0.0D0, 0.0D0)
 RHS_V(0,0,0) = DCMPLX(Re, 0.0D0)

 BC_V = DCMPLX(0.0D0, 0.0D0)

 CALL SOLVE (V,RHS_V,G_V,a_V(1),b_V(1),BC_V,WN_V,AD_V)

 CALL FOR_Q (QX_DPDX_1,V)

 DO iy = 0,NYS
  U_DPDX_1(iy) = DREAL(V(0,iy,0))
 ENDDO
!FOR W_DPDZ_1
 RHS_V = DCMPLX(0.0D0, 0.0D0)
 RHS_V(0,0,0) = DCMPLX(Re, 0.0D0)

 BC_V = DCMPLX(0.0D0, 0.0D0)

 CALL SOLVE (V,RHS_V,G_V,a_V(3),b_V(3),BC_V,WN_V,AD_V)

 CALL FOR_Q (QZ_DPDZ_1,V)

 DO iy = 0,NYS
  W_DPDZ_1(iy) = DREAL(V(0,iy,0))
 ENDDO

 DEALLOCATE(V,RHS_V,BC_V)

 RETURN
 END
! ---&---1---------2---------3---------4---------5---------6---------7--
!
 SUBROUTINE EQSOLVE_PV_1(PRE,RHS_PRE,BC_PRE, V,RHS_V,BC_V)

  USE MESH
 USE EQSOLVE_PV
 INCLUDE 'com_channel.h'

 COMPLEX*16 PRE(0:NXSP,0:NYS,-NZS:NZS),RHS_PRE(0:NXSP,0:NYS,-NZS:NZS),BC_PRE(0:NXSP,0:1,-NZS:NZS)
 COMPLEX*16 V(3,0:NXSP,0:NYS,-NZS:NZS),RHS_V(3,0:NXSP,0:NYS,-NZS:NZS),BC_V(3,0:NXSP,0:1,-NZS:NZS)

 COMPLEX*16, ALLOCATABLE :: TEMP(:,:,:,:)

 ALLOCATE( TEMP(3,0:NXSP,0:NYS,-NZS:NZS) )

 CALL SOLVE (PRE,RHS_PRE,G_pre,a_pre,b_pre,BC_PRE,WN_pre,AD_pre)

 CALL Grad(TEMP,PRE)

 TEMP = Re*(TEMP - RHS_V)

 DO I=1,3
  CALL SOLVE (V(I,0:NXSP,0:NYS,-NZS:NZS),TEMP(I,0:NXSP,0:NYS,-NZS:NZS),G_V,a_V(I),b_V(I),BC_V(I,0:NXSP,0:1,-NZS:NZS),WN_V,AD_V)
 ENDDO

 DEALLOCATE(TEMP)

 RETURN
 END
! ---&---1---------2---------3---------4---------5---------6---------7--
!
 SUBROUTINE EQSOLVE_PV_2(PRE,RHS_PRE,BC_PRE, V,RHS_V,BC_V)

  USE MESH
 USE EQSOLVE_PV
 INCLUDE 'com_channel.h'

 COMPLEX*16 PRE(0:NXSP,0:NYS,-NZS:NZS),RHS_PRE(0:NXSP,0:NYS,-NZS:NZS),BC_PRE(0:NXSP,0:1,-NZS:NZS)
 COMPLEX*16 V(3,0:NXSP,0:NYS,-NZS:NZS),RHS_V(3,0:NXSP,0:NYS,-NZS:NZS),BC_V(3,0:NXSP,0:1,-NZS:NZS)

 COMPLEX*16, ALLOCATABLE :: TEMP(:,:,:,:)
 COMPLEX*16, ALLOCATABLE :: Div_V_BOUND(:,:,:)

 ALLOCATE( TEMP(3,0:NXSP,0:NYS,-NZS:NZS),Div_V_BOUND(0:NXSP,0:1,-NZS:NZS) )

 CALL SOLVE (PRE,RHS_PRE,G_pre,a_pre,b_pre,BC_PRE,WN_pre,AD_pre)

 CALL Grad(TEMP,PRE)

 TEMP = Re*(TEMP - RHS_V)

 DO I=1,3
  CALL SOLVE (V(I,0:NXSP,0:NYS,-NZS:NZS),TEMP(I,0:NXSP,0:NYS,-NZS:NZS),G_V,a_V(I),b_V(I),BC_V(I,0:NXSP,0:1,-NZS:NZS),WN_V,AD_V)
 ENDDO

 CALL Div_BOUND(Div_V_BOUND,V)

 TEMP(1,0:NXSP,0,-NZS:NZS) = Div_V_BOUND(0:NXSP,1,-NZS:NZS)*Div_Vc_BOUND(2,0:NXSP,0,-NZS:NZS) - Div_V_BOUND(0:NXSP,0,-NZS:NZS)*Div_Vc_BOUND(2,0:NXSP,1,-NZS:NZS)
 TEMP(2,0:NXSP,0,-NZS:NZS) = Div_V_BOUND(0:NXSP,0,-NZS:NZS)*Div_Vc_BOUND(1,0:NXSP,1,-NZS:NZS) - Div_V_BOUND(0:NXSP,1,-NZS:NZS)*Div_Vc_BOUND(1,0:NXSP,0,-NZS:NZS)

 DO iy = 0,NYS
  DO I=1,3
   V(I,0:NXSP,iy,-NZS:NZS) = V(I,0:NXSP,iy,-NZS:NZS) + TEMP(1,0:NXSP,0,-NZS:NZS)*Vc(1,I,0:NXSP,iy,-NZS:NZS) &
               + TEMP(2,0:NXSP,0,-NZS:NZS)*Vc(2,I,0:NXSP,iy,-NZS:NZS)
  ENDDO
  PRE(0:NXSP,iy,-NZS:NZS) = PRE(0:NXSP,iy,-NZS:NZS) + TEMP(1,0:NXSP,0,-NZS:NZS)*PREc(1,0:NXSP,iy,-NZS:NZS) &
              + TEMP(2,0:NXSP,0,-NZS:NZS)*PREc(2,0:NXSP,iy,-NZS:NZS)
 ENDDO

 DEALLOCATE(TEMP,Div_V_BOUND)

 RETURN
 END
! ---&---1---------2---------3---------4---------5---------6---------7--
!
 SUBROUTINE SOLVE (U,F,G,a,b,BC,WN,AD)

 USE FFT
  USE MESH
 IMPLICIT REAL*8 (A-H, O-Z)

 COMPLEX*16, DIMENSION(0:NXSP,0:NYS,-NZS:NZS) :: U,F
 COMPLEX*16, DIMENSION(0:NXSP,0:1,-NZS:NZS) :: BC

 REAL*8 WN(0:NXSP,-NZS:NZS),AD(0:NXSP,2:NYS,-NZS:NZS)

 COMPLEX*16 COR

 REAL*8, ALLOCATABLE :: AD01(:,:,:)
 ALLOCATE( AD01(0:NXSP,0:1,-NZS:NZS) )

 IF(G==0.0D0 .AND. a==0.0D0 .AND. INDEX_MPI0==0) THEN
  COR = DCMPLX(0.0D0,0.0D0)
  DO iy=2,NYS,2
   COR = COR + F(0,iy,0)/DFLOAT(iy**2-1)
  ENDDO
  IF(b==0.0D0) THEN
   F(0,0,0) = COR
  ELSE
   COR = COR + (BC(0,0,0)+BC(0,1,0))/2.0D0/b
!---if to adjust BC---
!   BC(0,0:1,0) = BC(0,0:1,0) - b*( COR - F(0,0,0) )
!---if to adjust F---
   F(0,0,0) = COR
  ENDIF
 ENDIF
   
 U(0:NXSP,0,-NZS:NZS)=(BC(0:NXSP,0,-NZS:NZS)+BC(0:NXSP,1,-NZS:NZS) )/2.0D0
 U(0:NXSP,1,-NZS:NZS)=(BC(0:NXSP,0,-NZS:NZS)-BC(0:NXSP,1,-NZS:NZS) )/2.0D0
 DO iy = 2,NYS-2
  U(0:NXSP,iy,-NZS:NZS) = WL(iy)*F(0:NXSP,iy-2,-NZS:NZS) &
        + WD(iy)*F(0:NXSP, iy ,-NZS:NZS) &
        + WR(iy)*F(0:NXSP,iy+2,-NZS:NZS)
 ENDDO
 U(0:NXSP,NYS-1,-NZS:NZS) =  WL(NYS-1)*F(0:NXSP,NYS-3,-NZS:NZS) &
       + WD(NYS-1)*F(0:NXSP,NYS-1,-NZS:NZS)
 U(0:NXSP, NYS ,-NZS:NZS) =  WL( NYS )*F(0:NXSP,NYS-2,-NZS:NZS) &
       + WD( NYS )*F(0:NXSP, NYS ,-NZS:NZS)

 DO iy = NYS-2,2,-1
  U(0:NXSP,iy,-NZS:NZS) = U(0:NXSP,iy,-NZS:NZS) - WR(iy)*WN &
          *AD(0:NXSP,iy+2,-NZS:NZS)*U(0:NXSP,iy+2,-NZS:NZS)
 ENDDO

 DO iy = NYS,NYS-1,-1
  index = MOD(iy,2)
  IF(a==0.0D0 .AND. b==0.0D0) THEN
   FUNC = DFLOAT((1-index)*iy**2+index)
  ELSE
   FUNC = a+b*DFLOAT(iy**2)
  ENDIF
  AD01(0:NXSP,index,-NZS:NZS) = FUNC
 ENDDO

 DO iy = NYS-2,0,-1
  index = MOD(iy,2)
  IF(a==0.0D0 .AND. b==0.0D0) THEN
   FUNC = DFLOAT((1-index)*iy**2+index)
  ELSE
   FUNC = a+b*DFLOAT(iy**2)
  ENDIF
  U(0:NXSP,index,-NZS:NZS)=U(0:NXSP,index,-NZS:NZS)-AD01(0:NXSP,index,-NZS:NZS) &
           *AD(0:NXSP,iy+2,-NZS:NZS)*U(0:NXSP,iy+2,-NZS:NZS)
  AD01(0:NXSP,index,-NZS:NZS) = FUNC - WL(iy+2)*WN &
           *AD01(0:NXSP,index,-NZS:NZS)*AD(0:NXSP,iy+2,-NZS:NZS)
 ENDDO
    
 IF(G==0.0D0 .AND. a==0.0D0 .AND. INDEX_MPI0==0) THEN
  U(0,0,0) = DCMPLX(0.0D0,0.0D0)
  AD01(0,0,0) = 1.0D0
 ENDIF

 U(0:NXSP,0,-NZS:NZS) = U(0:NXSP,0,-NZS:NZS)/AD01(0:NXSP,0,-NZS:NZS)
 U(0:NXSP,1,-NZS:NZS) = U(0:NXSP,1,-NZS:NZS)/AD01(0:NXSP,1,-NZS:NZS)

 DO iy = 2,NYS
  U(0:NXSP,iy,-NZS:NZS) = (U(0:NXSP,iy,-NZS:NZS) - WL(iy)*WN &
        *U(0:NXSP,iy-2,-NZS:NZS))*AD(0:NXSP,iy,-NZS:NZS)
 ENDDO

 DEALLOCATE( AD01 )

 RETURN
 END
! ---&---1---------2---------3---------4---------5---------6---------7--
!
 SUBROUTINE SOLVE2 (U,F,G,a,b,BC,WN,AD)

 USE FFT
  USE MESH
 IMPLICIT REAL*8 (A-H, O-Z)

 COMPLEX*16, DIMENSION(0:NXSP,0:NYS,-NZS:NZS) :: U,F
 COMPLEX*16, DIMENSION(0:NXSP,0:1,-NZS:NZS) :: BC

 REAL*8  WN(0:NXSP,-NZS:NZS),AD(0:NXSP,2:NYS,-NZS:NZS),a(0:1),b(0:1)

 COMPLEX*16 C0,C1

 REAL*8, ALLOCATABLE :: AD01(:,:,:),BD01(:,:,:)
 ALLOCATE( AD01(0:NXSP,0:1,-NZS:NZS),BD01(0:NXSP,0:1,-NZS:NZS) )


 U(0:NXSP,0:1,-NZS:NZS) = BC(0:NXSP,0:1,-NZS:NZS)
 DO iy = 2,NYS-2
  U(0:NXSP,iy,-NZS:NZS) = WL(iy)*F(0:NXSP,iy-2,-NZS:NZS) &
        + WD(iy)*F(0:NXSP, iy ,-NZS:NZS) &
        + WR(iy)*F(0:NXSP,iy+2,-NZS:NZS)
 ENDDO
 U(0:NXSP,NYS-1,-NZS:NZS) =  WL(NYS-1)*F(0:NXSP,NYS-3,-NZS:NZS) &
       + WD(NYS-1)*F(0:NXSP,NYS-1,-NZS:NZS)
 U(0:NXSP, NYS ,-NZS:NZS) =  WL( NYS )*F(0:NXSP,NYS-2,-NZS:NZS) &
       + WD( NYS )*F(0:NXSP, NYS ,-NZS:NZS)

 DO iy = NYS-2,2,-1
  U(0:NXSP,iy,-NZS:NZS) = U(0:NXSP,iy,-NZS:NZS) - WR(iy)*WN &
          *AD(0:NXSP,iy+2,-NZS:NZS)*U(0:NXSP,iy+2,-NZS:NZS)
 ENDDO

 DO iy = NYS,NYS-1,-1
  index = MOD(iy,2)
  DFLOAT_iy2 = DFLOAT(iy**2)
  FUNC_AD = DFLOAT(1-index-index)*(a(index)+b(index)*DFLOAT_iy2)
  FUNC_BD = a(1-index)+b(1-index)*DFLOAT_iy2
  AD01(0:NXSP,index,-NZS:NZS) = FUNC_AD
  BD01(0:NXSP,index,-NZS:NZS) = FUNC_BD
 ENDDO

 DO iy = NYS-2,0,-1
  index = MOD(iy,2)
  DFLOAT_iy2 = DFLOAT(iy**2)
  FUNC_AD = DFLOAT(1-index-index)*(a(index)+b(index)*DFLOAT_iy2)
  FUNC_BD = a(1-index)+b(1-index)*DFLOAT_iy2
  U(0:NXSP,index,-NZS:NZS)=U(0:NXSP,index,-NZS:NZS)-AD01(0:NXSP,index,-NZS:NZS) &
           *AD(0:NXSP,iy+2,-NZS:NZS)*U(0:NXSP,iy+2,-NZS:NZS)
  U(0:NXSP,1-index,-NZS:NZS)=U(0:NXSP,1-index,-NZS:NZS)-BD01(0:NXSP,index,-NZS:NZS) &
           *AD(0:NXSP,iy+2,-NZS:NZS)*U(0:NXSP,iy+2,-NZS:NZS)
  AD01(0:NXSP,index,-NZS:NZS) = FUNC_AD - WL(iy+2)*WN &
           *AD01(0:NXSP,index,-NZS:NZS)*AD(0:NXSP,iy+2,-NZS:NZS)
  BD01(0:NXSP,index,-NZS:NZS) = FUNC_BD - WL(iy+2)*WN &
           *BD01(0:NXSP,index,-NZS:NZS)*AD(0:NXSP,iy+2,-NZS:NZS)
 ENDDO

 DO ikx=0,NXSP
  DO ikz=-NZS,NZS
   CTEMP=AD01(ikx,0,ikz)*AD01(ikx,1,ikz)-BD01(ikx,0,ikz)*BD01(ikx,1,ikz)
   C0=(U(ikx,0,ikz)*AD01(ikx,1,ikz)-U(ikx,1,ikz)*BD01(ikx,1,ikz))/CTEMP
   C1=(AD01(ikx,0,ikz)*U(ikx,1,ikz)-BD01(ikx,0,ikz)*U(ikx,0,ikz))/CTEMP
   U(ikx,0,ikz)=C0
   U(ikx,1,ikz)=C1
  END DO
 END DO

 DO iy = 2,NYS
  U(0:NXSP,iy,-NZS:NZS) = (U(0:NXSP,iy,-NZS:NZS) - WL(iy)*WN &
        *U(0:NXSP,iy-2,-NZS:NZS))*AD(0:NXSP,iy,-NZS:NZS)
 ENDDO

 DEALLOCATE( AD01,BD01 )

 RETURN
 END
