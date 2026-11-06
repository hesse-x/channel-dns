! ---&---1---------2---------3---------4---------5---------6---------7--
!
  CHARACTER*2 FUNCTION itochar2(iy)

  IMPLICIT REAL*8 (A-H, O-Z)

  WRITE(itochar2,'(I2)') iy

  END
! ---&---1---------2---------3---------4---------5---------6---------7--
!
  CHARACTER*3 FUNCTION itochar3(iy)

  IMPLICIT REAL*8 (A-H, O-Z)

  WRITE(itochar3,'(I3)') iy

  END
! ---&---1---------2---------3---------4---------5---------6---------7--
!
  CHARACTER*4 FUNCTION itochar4(iy)

  IMPLICIT REAL*8 (A-H, O-Z)

  WRITE(itochar4,'(I4)') iy

  END
! ---&---1---------2---------3---------4---------5---------6---------7--
!
  CHARACTER*8 FUNCTION Ftochar(Q)

  IMPLICIT REAL*8 (A-H, O-Z)

  WRITE(Ftochar,'(F8.2)') Q

  RETURN
  END
! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE DOTTIME(ADB,A,B)

   USE MESH
  IMPLICIT REAL*8 (A-H,O-Z)

  COMPLEX*16 ADB(0:NXSP,0:NYS,-NZS:NZS)
  COMPLEX*16, DIMENSION(3,0:NXSP,0:NYS,-NZS:NZS) :: A,B

  REAL*8, ALLOCATABLE :: A2(:,:,:,:),B2(:,:,:,:),ADB2(:,:,:)

  ALLOCATE (A2(3,NX2,0:NY2P,NZ2),B2(3,NX2,0:NY2P,NZ2),ADB2(NX2,0:NY2P,NZ2))

  DO I = 1,3
    CALL channel_xyzfft2 ('B',A2(I,1:NX2,0:NY2P,1:NZ2),A(I,0:NXSP,0:NYS,-NZS:NZS))
    CALL channel_xyzfft2 ('B',B2(I,1:NX2,0:NY2P,1:NZ2),B(I,0:NXSP,0:NYS,-NZS:NZS))
  ENDDO

  ADB2 =  A2(1,1:NX2,0:NY2P,1:NZ2)*B2(1,1:NX2,0:NY2P,1:NZ2) + &
      A2(2,1:NX2,0:NY2P,1:NZ2)*B2(2,1:NX2,0:NY2P,1:NZ2) + &
      A2(3,1:NX2,0:NY2P,1:NZ2)*B2(3,1:NX2,0:NY2P,1:NZ2)


  CALL channel_xyzfft2 ('F',ADB2,ADB)

  DEALLOCATE (A2,B2,ADB2)

  RETURN
  END
! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE CROSSTIME(ACB,A,B)

   USE MESH
  INCLUDE 'com_channel.h'

  COMPLEX*16, DIMENSION(3,0:NXSP,0:NYS,-NZS:NZS) :: ACB,A,B

  REAL*8, ALLOCATABLE :: A2(:,:,:,:),B2(:,:,:,:),ACB2(:,:,:,:)

  ALLOCATE(A2(3,NX2,0:NY2P,NZ2),B2(3,NX2,0:NY2P,NZ2),ACB2(3,NX2,0:NY2P,NZ2))

  DO I=1,3
    CALL channel_xyzfft2 ('B',A2(I,1:NX2,0:NY2P,1:NZ2),A(I,0:NXSP,0:NYS,-NZS:NZS))
    CALL channel_xyzfft2 ('B',B2(I,1:NX2,0:NY2P,1:NZ2),B(I,0:NXSP,0:NYS,-NZS:NZS))
  ENDDO

  DO I=1,3
    K1=E_ijk(I,1);  K2=E_ijk(I,2)
    ACB2(I,1:NX2,0:NY2P,1:NZ2) = A2(K1,1:NX2,0:NY2P,1:NZ2)*B2(K2,1:NX2,0:NY2P,1:NZ2) - &
                    A2(K2,1:NX2,0:NY2P,1:NZ2)*B2(K1,1:NX2,0:NY2P,1:NZ2)
  ENDDO

  DO I=1,3
    CALL channel_xyzfft2 ('F',ACB2(I,1:NX2,0:NY2P,1:NZ2),ACB(I,0:NXSP,0:NYS,-NZS:NZS))
  ENDDO

  DEALLOCATE(A2,B2,ACB2)

  RETURN
  END
! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE Grad(Grad_U,U)

   USE MESH
  IMPLICIT REAL*8 (A-H,O-Z)

  COMPLEX*16 Grad_U(3,0:NXSP,0:NYS,-NZS:NZS),U(0:NXSP,0:NYS,-NZS:NZS)

! [FIX 2026-08-14] 原函数体为损坏代码（DU 未声明、U 被按 4 维索引、混入 Div_U 行），
! 任何编译器均无法通过；按 Curl/DiV_Grad 惯例重建为标量梯度：Grad_U(I)=dU/dx_I。
  DO I=1,3
    CALL DERIV (I,Grad_U(I,0:NXSP,0:NYS,-NZS:NZS),U(0:NXSP,0:NYS,-NZS:NZS))
  ENDDO

  RETURN
  END
! ---&---1---------2---------3---------4---------5---------6---------7--
!
! [FIX 2026-08-14] 补回缺失的 Div 子程序（STEP23:55 与本文件 DiV_Grad:157 均调用，
! 但全源码无定义）。函数体即原 Grad 内混入的残片，语义为矢量场散度。
  SUBROUTINE Div(Div_U,U)

   USE MESH
  IMPLICIT REAL*8 (A-H,O-Z)

  COMPLEX*16 Div_U(0:NXSP,0:NYS,-NZS:NZS),U(3,0:NXSP,0:NYS,-NZS:NZS)

  COMPLEX*16, ALLOCATABLE :: DU(:,:,:,:)

  ALLOCATE(DU(3,0:NXSP,0:NYS,-NZS:NZS))

  DO I=1,3
    CALL DERIV (I,DU(I,0:NXSP,0:NYS,-NZS:NZS),U(I,0:NXSP,0:NYS,-NZS:NZS))
  ENDDO

  Div_U = DU(1,0:NXSP,0:NYS,-NZS:NZS) + DU(2,0:NXSP,0:NYS,-NZS:NZS) + DU(3,0:NXSP,0:NYS,-NZS:NZS)

  DEALLOCATE(DU)

  RETURN
  END
! ---&---1---------2---------3---------4---------5---------6---------7--
! rot 旋度
  SUBROUTINE Curl(Curl_U,U)

   USE MESH
  INCLUDE 'com_channel.h'

  COMPLEX*16 Curl_U(3,0:NXSP,0:NYS,-NZS:NZS),U(3,0:NXSP,0:NYS,-NZS:NZS)

  COMPLEX*16, ALLOCATABLE :: DU(:,:,:,:)

  ALLOCATE(DU(2,0:NXSP,0:NYS,-NZS:NZS))

  DO I=1,3
    K1=E_ijk(I,1);  K2=E_ijk(I,2)
    CALL DERIV (K1,DU(1,0:NXSP,0:NYS,-NZS:NZS),U(K2,0:NXSP,0:NYS,-NZS:NZS))
    CALL DERIV (K2,DU(2,0:NXSP,0:NYS,-NZS:NZS),U(K1,0:NXSP,0:NYS,-NZS:NZS))
    Curl_U(I,0:NXSP,0:NYS,-NZS:NZS) = DU(1,0:NXSP,0:NYS,-NZS:NZS) - DU(2,0:NXSP,0:NYS,-NZS:NZS)
  ENDDO

  DEALLOCATE(DU)

  RETURN
  END
! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE DiV_Grad(DiV_Grad_U,U)

   USE MESH
  IMPLICIT REAL*8 (A-H,O-Z)

  COMPLEX*16 DiV_Grad_U(0:NXSP,0:NYS,-NZS:NZS),U(0:NXSP,0:NYS,-NZS:NZS)

  COMPLEX*16, ALLOCATABLE :: Grad_U(:,:,:,:)

  ALLOCATE(Grad_U(3,0:NXSP,0:NYS,-NZS:NZS))

  CALL Grad(Grad_U,U)
  CALL Div(DiV_Grad_U,Grad_U)

  DEALLOCATE(Grad_U)

  RETURN
  END
! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE Div_BOUND(Div_V_BOUND,V)

   USE MESH
  IMPLICIT REAL*8 (A-H,O-Z)

  COMPLEX*16 Div_V_BOUND(0:NXSP,0:1,-NZS:NZS),V(3,0:NXSP,0:NYS,-NZS:NZS)

  COMPLEX*16, ALLOCATABLE :: DUDX(:,:,:),DWDZ(:,:,:),TEMP1(:,:),TEMP2(:,:)

  ALLOCATE(DUDX(0:NXSP,0:NYS,-NZS:NZS),DWDZ(0:NXSP,0:NYS,-NZS:NZS),TEMP1(0:NXSP,-NZS:NZS),TEMP2(0:NXSP,-NZS:NZS))

  CALL DERIV (1,DUDX,V(1,0:NXSP,0:NYS,-NZS:NZS))
  CALL DERIV (3,DWDZ,V(3,0:NXSP,0:NYS,-NZS:NZS))

  Div_V_BOUND=DCMPLX(0.0D0,0.0D0)
  DO iy=0,NYS
    TEMP1=DUDX(0:NXSP,iy,-NZS:NZS)+DWDZ(0:NXSP,iy,-NZS:NZS)
    TEMP2=DFLOAT(iy**2)*V(2,0:NXSP,iy,-NZS:NZS)
    Div_V_BOUND(0:NXSP,0,-NZS:NZS) = Div_V_BOUND(0:NXSP,0,-NZS:NZS) + TEMP1 + TEMP2
    IF(MOD(iy,2)==0) THEN
      Div_V_BOUND(0:NXSP,1,-NZS:NZS) = Div_V_BOUND(0:NXSP,1,-NZS:NZS) + TEMP1 - TEMP2
    ELSE
      Div_V_BOUND(0:NXSP,1,-NZS:NZS) = Div_V_BOUND(0:NXSP,1,-NZS:NZS) - TEMP1 + TEMP2
    ENDIF
  ENDDO

  DEALLOCATE(DUDX,DWDZ,TEMP1,TEMP2)

  RETURN
  END
! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE DISTANCE_L8 (DISTANCE,U,V)

   USE MESH
  IMPLICIT REAL*8 (A-H,O-Z)

  COMPLEX*16, DIMENSION(0:NXSP,0:NYS,-NZS:NZS) :: U,V

  DISTANCE = 0.0D0
  DO ix = 0,NXSP
    DO iy = 0,NYS
      DO iz = -NZS,NZS
        DISTANCE = DMAX1(DISTANCE,CDABS(U(ix,iy,iz)-V(ix,iy,iz)))
      ENDDO
    ENDDO
  ENDDO

  RETURN
  END
! ---&---1---------2---------3---------4---------5---------6---------7--
! deriv 散度
  SUBROUTINE DERIV (IC,DF,F)

  USE FFT
   USE MESH
  IMPLICIT REAL*8 (A-H, O-Z)

  COMPLEX*16, DIMENSION(0:NXSP,0:NYS,-NZS:NZS) :: F,DF

  IF(IC==1) GO TO 10
  IF(IC==2) GO TO 20
  IF(IC==3) GO TO 30

10  DO ix = 0,NXSP
    DF(ix,0:NYS,-NZS:NZS) = CWAVX(ix)*F(ix,0:NYS,-NZS:NZS)
  ENDDO

  RETURN

20  DF(0:NXSP, NYS ,-NZS:NZS) = DCMPLX(0.0D0, 0.0D0)
  DF(0:NXSP,NYS-1,-NZS:NZS) = DFLOAT(NYS+NYS)*F(0:NXSP,NYS,-NZS:NZS)

  DO iy = NYS-2,0,-1
    iy1=iy+1
    DF(0:NXSP,iy,-NZS:NZS) = DF(0:NXSP,iy+2,-NZS:NZS) + DFLOAT(iy1+iy1)*F(0:NXSP,iy1,-NZS:NZS)
  ENDDO

  DF(0:NXSP,0,-NZS:NZS) = DF(0:NXSP,0,-NZS:NZS)/2.0D0

  RETURN

30  DO iz = -NZS,NZS
    DF(0:NXSP,0:NYS,iz) = CWAVZ(iz)*F(0:NXSP,0:NYS,iz)
  ENDDO

  RETURN
  END
! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE FOR_Q (Q,V)

  USE MESH
  IMPLICIT REAL*8 (A-H, O-Z)

  COMPLEX*16 V(0:NXSP,0:NYS,-NZS:NZS)

  Q=0.0D0
  DO iy=0,NYS,2
    Q = Q - DREAL(V(0,iy,0))/DFLOAT(iy**2-1)
  ENDDO

  RETURN
  END
