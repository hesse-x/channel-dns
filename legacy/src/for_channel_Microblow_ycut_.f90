!----------------
! Initial section
!----------------
! ---&---1---------2---------3---------4---------5---------6---------7--
!
	SUBROUTINE INITIAL (V,EPS)   !给定初场

	USE MESH
	INCLUDE 'com_channel.h'

	COMPLEX*16 V(3,0:NXSP,0:NYS,-NZS:NZS)

	REAL*8, ALLOCATABLE :: Vp(:,:,:,:)

	ALLOCATE( Vp(3,NX,0:NYP,NZ) )

	Vp=0.0D0

	IF(INDEX_GRAD_PRE==1) THEN
		uc=Re/2.0D0
	ELSE
		uc=1.5D0
	ENDIF
	ucx=uc*DCOS(THETA_RAD)
	ucz=uc*DSIN(THETA_RAD)

	DO iy = 0,NYP
		uy = ucx*(1.0D0-yc(iy)**2)
		Vp(1,1:NX,iy,1:NZ) = uy
		Vp(2,1:NX,iy,1:NZ) = V_Microblow0
		uy = ucz*(1.0D0-yc(iy)**2)
		Vp(3,1:NX,iy,1:NZ) = uy
	ENDDO

	DO I=1,3
		CALL channel_xyzfft ('F',Vp(I,1:NX,0:NYP,1:NZ),V(I,0:NXSP,0:NYS,-NZS:NZS))
	ENDDO

! [FIX 2026-08-15] EPS=0 is a strict no-op; avoid a roundtrip that seeds numerical noise.
	IF(EPS/=0.0D0) CALL VDISTURB (V,EPS)

	DEALLOCATE(Vp)

	RETURN
	END
! ---&---1---------2---------3---------4---------5---------6---------7--
!
	SUBROUTINE VDISTURB (V,EPS)   !施加扰动

	USE MESH
	INCLUDE 'com_channel.h'

	COMPLEX*16 V(3,0:NXSP,0:NYS,-NZS:NZS)

	REAL*8, ALLOCATABLE :: Vp(:,:,:,:)

	INTEGER DATE_TIME (8)
	CHARACTER (LEN = 12) REAL_CLOCK (3)

	CALL DATE_AND_TIME (REAL_CLOCK (1), REAL_CLOCK (2), REAL_CLOCK (3), DATE_TIME)

 	XRAND=RAND(1000*my_id*DATE_TIME (8)+16*DATE_TIME (7))

	ALLOCATE( Vp(3,NX,0:NYP,NZ) )

	DO I=1,3
		CALL channel_xyzfft ('B',Vp(I,1:NX,0:NYP,1:NZ),V(I,0:NXSP,0:NYS,-NZS:NZS))
	ENDDO

	IF(INDEX_GRAD_PRE==0) THEN
		uc=1.5D0
	ELSE
		uc=Re/2.0D0
	ENDIF

	DO iy = 0,NYP
		uy = uc*(1.0D0-yc(iy)**2)
		DO iz = 1,NZ
			DO ix = 1,NX
				Vp(1,ix,iy,iz) = Vp(1,ix,iy,iz)+uy*EPS*(RAND(0)-0.5D0)
				Vp(2,ix,iy,iz) = Vp(2,ix,iy,iz)+uy*EPS*(RAND(0)-0.5D0)
        WRITE(*,*) INDEX_disturb_Z
				Vp(3,ix,iy,iz) = Vp(3,ix,iy,iz)+uy*EPS*(RAND(0)-0.5D0)*DFLOAT(INDEX_disturb_Z)
			ENDDO
		ENDDO
	ENDDO

	DO I=1,3
		CALL channel_xyzfft ('F',Vp(I,1:NX,0:NYP,1:NZ),V(I,0:NXSP,0:NYS,-NZS:NZS))
	ENDDO

	DEALLOCATE(Vp)

	RETURN
	END
!---------------------------
! End of the Initial section
!---------------------------





!-----------------
! Boundary section
!-----------------
! ---&---1---------2---------3---------4---------5---------6---------7--
!
	SUBROUTINE VEL_BOUND (BC_V,TIME)

  	USE MESH
	INCLUDE 'com_channel.h'

	COMPLEX*16 BC_V(3,0:NXSP,0:1,-NZS:NZS)

! 	REAL*8, ALLOCATABLE :: BC_Vp(:,:,:,:)
! 	COMPLEX*16, ALLOCATABLE ::BC_V_TEMP(:,:)
! 
! 	ALLOCATE( BC_Vp(3,NX,0:1,NZ) )
! 	ALLOCATE( BC_V_TEMP(0:NXSP,-NZS:NZS) )

	BC_V = DCMPLX(0.0D0, 0.0D0)
 	IF(INDEX_MPI0==0) BC_V(2,0,0:1,0) = DCMPLX(V_Microblow0, 0.0D0)

! 	DEALLOCATE(BC_Vp,BC_V_TEMP)
	
	RETURN
	END
!----------------------------
! End of the Boundary section
!----------------------------



!----------------
! output section
!----------------
! ---&---1---------2---------3---------4---------5---------6---------7--
!
	SUBROUTINE OUTPUT (NTPit,TIME,V,PRE, Nstat)

	USE MESH
	IMPLICIT REAL*8 (A-H,O-Z)

	INCLUDE 'mpi_head.h'

	character*2,external :: itochar2
	character*3,external :: itochar3
	character*4,external :: itochar4

	COMPLEX*16 V(3,0:NXSP,0:NYS,-NZS:NZS)
	COMPLEX*16 PRE(0:NXSP,0:NYS,-NZS:NZS)

	REAL*8, ALLOCATABLE :: Vp(:,:,:,:),PREp(:,:,:)

	ALLOCATE( Vp(3,NX,0:NYP,NZ),PREp(NX,0:NYP,NZ) )

	DO I=1,3
		CALL channel_xyzfft ('B',Vp(I,1:NX,0:NYP,1:NZ),V(I,0:NXSP,0:NYS,-NZS:NZS))
	ENDDO
	CALL channel_xyzfft ('B',PREp,PRE)

 	EKT_myid=0.0D0
	DO iy=0,NYP
		DO ix=1,NX
			DO iz=1,NZ
				EKT_myid=EKT_myid+Vp(2,ix,iy,iz)**2
			ENDDO
		ENDDO
	ENDDO

	IF(INDEX_Serial==1) THEN
		EKT=EKT_myid
	ELSE
		CALL MPI_REDUCE(EKT_myid,EKT,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
	ENDIF

	IF(INDEX_MPI0==0) THEN
		OPEN(61,FILE='EKT'//itochar2(my_id)//'.dat',STATUS='UNKNOWN',ACCESS='APPEND')
			WRITE(61,*) TIME,EKT
		CLOSE(61)
	ENDIF

	OPEN(60,FILE='field_P_'//itochar2(my_id)//'.dat',STATUS='UNKNOWN')
		WRITE(60,*) NTPit,TIME
		WRITE(60,*) Vp,PREp
	CLOSE(60)

 	IF(NTPit>Nstat) THEN
		CALL FOR_STATISTICS (Vp,Prep)
	ENDIF

	DEALLOCATE(Vp,PREp)

	END
!---------------------------
! End of the output section
!---------------------------
