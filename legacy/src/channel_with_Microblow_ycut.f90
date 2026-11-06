! ---&---1---------2---------3---------4---------5---------6---------7--
!
  MODULE MESH
  INTEGER,SAVE :: NX,  NXS, NX2,  NXH,NXS2P1,    NXSP,NXSP2P1
  INTEGER,SAVE :: NY,  NYS, NY2,          NYP,NY2P
  INTEGER,SAVE :: NZ,  NZS, NZ2,  NZH,NZS2P1
  INTEGER,SAVE :: NXTNZ
  INTEGER,SAVE :: INDEX_Serial,my_id,NPROC,INDEX_MPI0
  INTEGER,ALLOCATABLE,SAVE :: LSUBX(:),NXSPT(:),LSUBX2(:),NXSP2P1T(:)
  INTEGER,ALLOCATABLE,SAVE :: LSUBY(:),NYPT(:), LSUBY2(:),NY2PT(:)
  REAL*8, ALLOCATABLE,SAVE :: XC(:),YC(:),ZC(:)
  END MODULE
! ---&---1---------2---------3---------4---------5---------6---------7--
!
  MODULE FFT
  INTEGER,ALLOCATABLE,SAVE :: itrigx(:),itrigx2(:)
  REAL*8,ALLOCATABLE,SAVE :: trigx(:),trigx2(:),WAVX(:)
  INTEGER,ALLOCATABLE,SAVE :: itrigy(:),itrigy2(:)
  REAL*8,ALLOCATABLE,SAVE :: trigy(:),trigy2(:),WL(:),WD(:),WR(:)
  INTEGER,ALLOCATABLE,SAVE :: itrigz(:),itrigz2(:)
  REAL*8,ALLOCATABLE,SAVE :: trigz(:),trigz2(:),WAVZ(:)
  COMPLEX*16,ALLOCATABLE,SAVE :: CWAVX(:),CWAVZ(:),CWAVXRe(:),CWAVZRe(:)
  END MODULE
! ---&---1---------2---------3---------4---------5---------6---------7--
!
  MODULE EQSOLVE_PV
  REAL*8,SAVE :: G_pre,a_pre,b_pre, G_V,a_V(3),b_V(3), QX_DPDX_1,QZ_DPDZ_1
  REAL*8,ALLOCATABLE,SAVE :: WN_pre(:,:),AD_pre(:,:,:),WN_V(:,:),AD_V(:,:,:), U_DPDX_1(:),W_DPDZ_1(:)
  COMPLEX*16,ALLOCATABLE,SAVE :: PREc(:,:,:,:),Vc(:,:,:,:,:),Div_Vc_BOUND(:,:,:,:)
  COMPLEX*16,ALLOCATABLE,SAVE :: PRE0(:,:,:),V0(:,:,:,:),SIGMA0(:,:,:)
  END MODULE
! ---&---1---------2---------3---------4---------5---------6---------7--
!
  MODULE STATISTICS
  INTEGER,SAVE :: ICOUNT
  REAL*8,SAVE :: DPDX,DPDZ,DPDX_mean,DPDZ_mean
!----1
  REAL*8,ALLOCATABLE,SAVE :: V_mean(:,:),Pre_mean(:)
  REAL*8,ALLOCATABLE,SAVE :: VV_mean(:,:),Pre_V_mean(:,:),Pre_Pre_mean(:)
  REAL*8,ALLOCATABLE,SAVE :: VVV_mean(:,:)
  REAL*8,ALLOCATABLE,SAVE :: VVVV_mean(:,:)
!----2
!   REAL*8,ALLOCATABLE,SAVE :: V_mean(:,:,:),Pre_mean(:,:)
!   REAL*8,ALLOCATABLE,SAVE :: VV_mean(:,:,:),Pre_V_mean(:,:,:),Pre_Pre_mean(:,:)
!   REAL*8,ALLOCATABLE,SAVE :: VVV_mean(:,:,:)
!   REAL*8,ALLOCATABLE,SAVE :: VVVV_mean(:,:,:)
!----3
!   REAL*8,ALLOCATABLE,SAVE :: V_mean(:,:,:,:),Pre_mean(:,:,:)
!   REAL*8,ALLOCATABLE,SAVE :: VV_mean(:,:,:,:),Pre_V_mean(:,:,:,:),Pre_Pre_mean(:,:,:)
!   REAL*8,ALLOCATABLE,SAVE :: VVV_mean(:,:,:,:)
!   REAL*8,ALLOCATABLE,SAVE :: VVVV_mean(:,:,:,:)
  END MODULE
! ---&---1---------2---------3---------4---------5---------6---------7--
!
  PROGRAM channel_3d_with_LBC_Vreduce

  USE FFT
  USE MESH
  USE EQSOLVE_PV
  USE STATISTICS
  INCLUDE 'com_channel.h'

  INCLUDE 'mpi_head.h'

  character*2,external :: itochar2
  character*3,external :: itochar3
  character*4,external :: itochar4

  REAL*8  a(0:1),b(0:1)

  COMPLEX*16, ALLOCATABLE :: PRE(:,:,:),BC_PRE(:,:,:),BC_PREm(:,:,:)
  COMPLEX*16, ALLOCATABLE :: V(:,:,:,:),Vm(:,:,:,:),BC_V(:,:,:,:)
  COMPLEX*16, ALLOCATABLE :: NL_V(:,:,:,:),NL_Vm(:,:,:,:)
!   COMPLEX*16, ALLOCATABLE :: Vn(:,:,:,:),Vnm(:,:,:,:)

!---
!   my_id = 0;  NPROC = 1
  call MPI_Init(ierr)
  call MPI_Comm_Rank(MPI_Comm_World,my_id,ierr)
  call MPI_Comm_Size(MPI_Comm_World,NPROC,ierr)

  CALL CPU_TIME(time1)

  OPEN(10,FILE='meshparm_channel_with_Microblow.txt',STATUS="UNKNOWN")
    READ(10,*)
    READ(10,*) NX,NY,NZ,INDEX_disturb_Z
    READ(10,*)
    READ(10,*) NXS,NYS,NZS
    READ(10,*)
    READ(10,*) NX2,NY2,NZ2
    READ(10,*)
    READ(10,*) INDEX_Serial
  CLOSE(10)
  NXH=NX/2;  NZH=NZ/2
  NXS2P1=2*NXS+1;  NZS2P1=2*NZS+1

  INDEX_MPI0 = (1-INDEX_Serial)*my_id

  NPART=INDEX_Serial+(1-INDEX_Serial)*NPROC

  ALLOCATE( LSUBX(0:NPART-1),NXSPT(0:NPART-1),LSUBX2(0:NPART-1),NXSP2P1T(0:NPART-1) )

  NXSP = (NXS+1)/NPART
  NXSPT=NXSP
  K=MOD(NXS+1,NPART)
  IF(K>0) THEN
    DO I=NPART-K,NPART-1
      NXSPT(I)=NXSP+1
    ENDDO
  ENDIF
  LSUBX(0)=0
  DO I=1,NPART-1
    LSUBX(I)=LSUBX(I-1)+NXSPT(I-1)
  ENDDO
  LSUBX2 = LSUBX + LSUBX
  NXSPT=NXSPT-1;  NXSP2P1T=NXSPT+NXSPT+1
  NXSP=NXSPT(INDEX_MPI0);  NXSP2P1=NXSP2P1T(INDEX_MPI0)

  ALLOCATE( LSUBY(0:NPART-1),NYPT(0:NPART-1) )

  NYP = (NY+1)/NPART
  NYPT=NYP
  K=MOD(NY+1,NPART)
  IF(K>0) THEN
    DO I=NPART-K,NPART-1
      NYPT(I)=NYP+1
    ENDDO
  ENDIF
  LSUBY(0)=0
  DO I=1,NPART-1
    LSUBY(I)=LSUBY(I-1)+NYPT(I-1)
  ENDDO
  NYPT=NYPT-1
  NYP=NYPT(INDEX_MPI0)

  ALLOCATE( LSUBY2(0:NPART-1),NY2PT(0:NPART-1) )

  NY2P = (NY2+1)/NPART
  NY2PT=NY2P
  K=MOD(NY2+1,NPART)
  IF(K>0) THEN
    DO I=NPART-K,NPART-1
      NY2PT(I)=NY2P+1
    ENDDO
  ENDIF
  LSUBY2(0)=0
  DO I=1,NPART-1
    LSUBY2(I)=LSUBY2(I-1)+NY2PT(I-1)
  ENDDO
  NY2PT=NY2PT-1
  NY2P=NY2PT(INDEX_MPI0)

  NXTNZ=NX*NZ

  ALLOCATE( XC(NX),YC(0:NYP),ZC(NZ) )

  ALLOCATE( itrigx(15),itrigx2(15) )
  ALLOCATE( trigx(2*NX),trigx2(2*NX2),WAVX(0:NXSP) )
  ALLOCATE( itrigy(15),itrigy2(15) )
  ALLOCATE( trigy(3*NY+3),trigy2(3*NY2+3),WL(0:NYS),WD(0:NYS),WR(0:NYS) )
  ALLOCATE( itrigz(15),itrigz2(15) )
  ALLOCATE( trigz(4*NZ),trigz2(4*NZ2),WAVZ(-NZS:NZS) )

  ALLOCATE( CWAVX(0:NXSP),CWAVZ(-NZS:NZS),CWAVXRE(0:NXSP),CWAVZRE(-NZS:NZS) )

  ALLOCATE( WN_pre(0:NXSP,-NZS:NZS),AD_pre(0:NXSP,2:NYS,-NZS:NZS), WN_V(0:NXSP,-NZS:NZS),AD_V(0:NXSP,2:NYS,-NZS:NZS) )
  ALLOCATE( PREc(2,0:NXSP,0:NYS,-NZS:NZS),Vc(2,3,0:NXSP,0:NYS,-NZS:NZS),&
           &Div_Vc_BOUND(2,0:NXSP,0:1,-NZS:NZS),U_DPDX_1(0:NYS),W_DPDZ_1(0:NYS) )
  ALLOCATE( PRE0(0:NXSP,0:NYS,-NZS:NZS),V0(3,0:NXSP,0:NYS,-NZS:NZS),SIGMA0(0:NXSP,0:1,-NZS:NZS) )

   ALLOCATE( V_mean(3,0:NYP),Pre_mean(0:NYP) )
   ALLOCATE( VV_mean(6,0:NYP),Pre_V_mean(3,0:NYP),Pre_Pre_mean(0:NYP) )
   ALLOCATE( VVV_mean(3,0:NYP) )
   ALLOCATE( VVVV_mean(3,0:NYP) )

  ALLOCATE( PRE(0:NXSP,0:NYS,-NZS:NZS),BC_PRE(0:NXSP,0:1,-NZS:NZS),BC_PREm(0:NXSP,0:1,-NZS:NZS) )
  ALLOCATE( V(3,0:NXSP,0:NYS,-NZS:NZS),Vm(3,0:NXSP,0:NYS,-NZS:NZS),BC_V(3,0:NXSP,0:1,-NZS:NZS) )
  ALLOCATE( NL_V(3,0:NXSP,0:NYS,-NZS:NZS),NL_Vm(3,0:NXSP,0:NYS,-NZS:NZS) )
!   ALLOCATE( Vn(3,0:NXSP,0:NYS,-NZS:NZS),Vnm(3,0:NXSP,0:NYS,-NZS:NZS) )

! [FIX 2026-08-14] 阶段0排查发现：it=1 时 a(1)=b(1)=0，但 NL_Vm=...+0*NL_Vm、
! BC_PREm=...+0*BC_PREm、以及 Vm 仍被引用；未初始化内存若含 NaN 位型则 0*NaN=NaN
! 污染全场（MPI 下 rank1 堆复用脏内存触发，串行恰好拿到零页所以"能跑"）。
! 显式置零，数学上严格等价（系数为 0）。
  Vm=DCMPLX(0.0D0,0.0D0);  NL_Vm=DCMPLX(0.0D0,0.0D0);  BC_PREm=DCMPLX(0.0D0,0.0D0)

  OPEN(10,FILE='parm_channel_with_Microblow.txt',STATUS="UNKNOWN")
    READ(10,*)
    READ(10,*) XL,ZL,THETA
    READ(10,*)
    READ(10,*) INDEX_GRAD_PRE,Re
    READ(10,*)
    READ(10,*) INDEX_UNSTEADY_BC
    READ(10,*)
    READ(10,*) G_pre,a_pre0,b_pre0
    READ(10,*)
    DO I=1,3
      READ(10,*) a_V(I),b_V(I)
    ENDDO
  CLOSE(10)
  THETA_RAD=THETA*DACOS(-1.0D0)/180.0D0
  IF(INDEX_GRAD_PRE==1) THEN
    GRAD_PRE_X=-DCOS(THETA_RAD)
    GRAD_PRE_Z=-DSIN(THETA_RAD)
  ELSE
    QXmean=DCOS(THETA_RAD)
    QZmean=DSIN(THETA_RAD)
  ENDIF

!---6 READ DT AND NSTEPS
  OPEN(10,FILE='control_channel_with_Microblow.txt',STATUS="UNKNOWN")
    READ(10,*)
    READ(10,*) INDEX_CONT,EPS_init
    READ(10,*)
    READ(10,*) NSTEPS,NDUMFRQ,NPRNFRQ,DT
    READ(10,*)
    READ(10,*) NCOUNT_disturb,NDISTURBFRQ,EPS_disturb
    READ(10,*)
    READ(10,*) V_Microblow0
    READ(10,*)
    READ(10,*) Nstat
    CLOSE(10)
!---6 READ DT AND NSTEPS

  CALL SETUP

! [TEMP 2026-08-14] 阶段1对照补丁：打印参数解析值 + 分块表 + 坐标，
! 与 C++ 骨架程序的 DUMP 行逐值对照（tools/dump_compare.py）。对拍前可移除。
  WRITE(*,*) 'DUMP rank', my_id, ' meshparm', NX,NY,NZ,INDEX_disturb_Z,NXS,NYS,NZS,NX2,NY2,NZ2,INDEX_Serial
  WRITE(*,*) 'DUMP rank', my_id, ' parm1', XL,ZL,THETA,INDEX_GRAD_PRE,Re,INDEX_UNSTEADY_BC
  WRITE(*,*) 'DUMP rank', my_id, ' parm2', G_pre,a_pre0,b_pre0
  WRITE(*,*) 'DUMP rank', my_id, ' parm3', a_V(1),b_V(1),a_V(2),b_V(2),a_V(3),b_V(3)
  WRITE(*,*) 'DUMP rank', my_id, ' ctrl', INDEX_CONT,EPS_init,NSTEPS,NDUMFRQ,NPRNFRQ,DT, &
             NCOUNT_disturb,NDISTURBFRQ,EPS_disturb,V_Microblow0,Nstat
  IF(INDEX_GRAD_PRE==1) THEN
    WRITE(*,*) 'DUMP rank', my_id, ' derived', THETA_RAD,GRAD_PRE_X,GRAD_PRE_Z
  ELSE
    WRITE(*,*) 'DUMP rank', my_id, ' derived', THETA_RAD,QXmean,QZmean
  ENDIF
  WRITE(*,*) 'DUMP rank', my_id, ' sizes', NXH,NZH,NXS2P1,NZS2P1,INDEX_MPI0,NPART,NXSP,NXSP2P1,NYP,NY2P,NXTNZ
  WRITE(*,*) 'DUMP rank', my_id, ' LSUBX', LSUBX
  WRITE(*,*) 'DUMP rank', my_id, ' NXSPT', NXSPT
  WRITE(*,*) 'DUMP rank', my_id, ' LSUBX2', LSUBX2
  WRITE(*,*) 'DUMP rank', my_id, ' NXSP2P1T', NXSP2P1T
  WRITE(*,*) 'DUMP rank', my_id, ' LSUBY', LSUBY
  WRITE(*,*) 'DUMP rank', my_id, ' NYPT', NYPT
  WRITE(*,*) 'DUMP rank', my_id, ' LSUBY2', LSUBY2
  WRITE(*,*) 'DUMP rank', my_id, ' NY2PT', NY2PT
  WRITE(*,*) 'DUMP rank', my_id, ' XC', XC
  WRITE(*,*) 'DUMP rank', my_id, ' YC', YC
  WRITE(*,*) 'DUMP rank', my_id, ' ZC', ZC

! [TEMP 2026-08-15] 阶段2对照补丁：波数与 Chebyshev 权重（SETUP 36-67 行计算结果），
! 与 C++ 侧 --print-fft 的 DUMP 行逐值对照（<1e-15）。对拍前可移除。
  WRITE(*,*) 'DUMP rank', my_id, ' WAVX', WAVX
  WRITE(*,*) 'DUMP rank', my_id, ' WAVZ', WAVZ
  WRITE(*,*) 'DUMP rank', my_id, ' WL', WL(0:NYS)
  WRITE(*,*) 'DUMP rank', my_id, ' WD', WD(0:NYS)
  WRITE(*,*) 'DUMP rank', my_id, ' WR', WR(0:NYS)

  IF(INDEX_CONT==0) THEN
    NT=0;   TIME=0.0D0
    CALL INITIAL(V,EPS_init)
  ELSE
    OPEN(50,FILE='field_S_'//itochar2(my_id)//'.dat',STATUS='UNKNOWN')
      READ(50,*) NT,TIME
      READ(50,*) V
    CLOSE(50)
  ENDIF

  IF(INDEX_UNSTEADY_BC==0) CALL VEL_BOUND (BC_V,TIME)

  ICOUNT=0
!---
  DO 29 it=1,NSTEPS

    IF(it<3) THEN
      a(0) = DFLOAT(it);  a(1) = DFLOAT(1-it)/2.0D0
      b(0) = DFLOAT(it);  b(1) = DFLOAT(1-it)
      a = a/DT;  G_V = (a(0)+a(1))*Re
      a_pre=1.0D0;  b_pre=0.0D0
       CALL EQSOLVE_PV_setup
      a_pre=a_pre0;  b_pre=b_pre0
    ENDIF

    NTPit = NT + it;  TIME = TIME + DT
    IF(INDEX_MPI0==0) WRITE(*,*) 'it = ', it,'    TIME = ', TIME

    CALL NLTERM(NL_V,V, BC_PRE)
    BC_PREm = b(0)*BC_PRE+b(1)*BC_PREm
    NL_Vm = a(0)*V+a(1)*Vm + (b(0)*NL_V+b(1)*NL_Vm)
    Vm=V
    IF(INDEX_UNSTEADY_BC==1) CALL VEL_BOUND (BC_V,TIME)

    CALL STEP23 (PRE,BC_PREm,V,NL_Vm,BC_V)

    IF(INDEX_MPI0==0) THEN
      V(2,0,0,0) = (BC_V(2,0,0,0)+BC_V(2,0,1,0))/2.0D0
      V(2,0,1:NYS,0) = DCMPLX(0.0D0,0.0D0)
    ENDIF

    NL_Vm=NL_V;  BC_PREm=BC_PRE

    IF(INDEX_MPI0==0) THEN
      CALL FOR_Q(QX,V(1,0:NXSP,0:NYS,-NZS:NZS))
      CALL FOR_Q(QZ,V(3,0:NXSP,0:NYS,-NZS:NZS))
      IF(INDEX_GRAD_PRE==1) THEN
        QX = QX + GRAD_PRE_X*QX_DPDX_1
        QZ = QZ + GRAD_PRE_Z*QZ_DPDZ_1
        V(1,0,0:NYS,0) = V(1,0,0:NYS,0) + GRAD_PRE_X*U_DPDX_1
        V(3,0,0:NYS,0) = V(3,0,0:NYS,0) + GRAD_PRE_Z*W_DPDZ_1
      ELSE
        DPDX = (QXmean-QX)/QX_DPDX_1
         DPDZ = (QZmean-QZ)/QZ_DPDZ_1
        V(1,0,0:NYS,0) = V(1,0,0:NYS,0) + DPDX*U_DPDX_1
         V(3,0,0:NYS,0) = V(3,0,0:NYS,0) + DPDZ*W_DPDZ_1
      ENDIF
    ENDIF

    IF(MOD(it,NDUMFRQ)==0) THEN
      OPEN(60,FILE='field_S_'//itochar2(my_id)//'.dat',STATUS='UNKNOWN')
        WRITE(60,*) NTPit,TIME
        WRITE(60,*) V
      CLOSE(60)
    ENDIF

    IF(MOD(it,NPRNFRQ)==0) CALL OUTPUT(NTPit,TIME,V,PRE, Nstat)

    IF(NCOUNT_disturb>0 .AND. MOD(it,NDISTURBFRQ)==0) THEN
      CALL VDISTURB (V,EPS_disturb)
      NCOUNT_disturb = NCOUNT_disturb - 1
    ENDIF

29  CONTINUE

  DEALLOCATE( LSUBX,NXSPT,LSUBX2,NXSP2P1T )
  DEALLOCATE( LSUBY,NYPT,LSUBY2,NY2PT )
  DEALLOCATE( XC,YC,ZC )

  DEALLOCATE( itrigx,itrigx2 )
  DEALLOCATE( trigx,trigx2,WAVX )
  DEALLOCATE( itrigy,itrigy2 )
  DEALLOCATE( trigy,trigy2,WL,WD,WR )
  DEALLOCATE( itrigz,itrigz2 )
  DEALLOCATE( trigz,trigz2,WAVZ )

  DEALLOCATE( CWAVX,CWAVZ,CWAVXRE,CWAVZRE )

  DEALLOCATE( WN_pre,AD_pre, WN_V,AD_V )
  DEALLOCATE( PREc,Vc,Div_Vc_BOUND,U_DPDX_1,W_DPDZ_1 )
  DEALLOCATE( PRE0,V0,SIGMA0 )

   DEALLOCATE( V_mean,Pre_mean )
   DEALLOCATE( VV_mean,Pre_V_mean,Pre_Pre_mean )
   DEALLOCATE( VVV_mean )
   DEALLOCATE( VVVV_mean )

  DEALLOCATE( PRE,BC_PRE,BC_PREm )
  DEALLOCATE( V,Vm,BC_V )
  DEALLOCATE( NL_V,NL_Vm )
!   DEALLOCATE( Vn,Vnm )
  
  CALL CPU_TIME(time2)
  WRITE(*,*) 'Time = ', time2-time1

  call MPI_finalize(ierr)

  END
