! ---&---1---------2---------3---------4---------5---------6---------7--
!
	SUBROUTINE SETUP

	USE FFT
	USE MESH
	INCLUDE 'com_channel.h'

	REAL*8, ALLOCATABLE :: YC0TONY(:)
	ALLOCATE(YC0TONY(0:NY))

	E_ijk(1,1)=2;  E_ijk(1,2)=3
	E_ijk(2,1)=3;  E_ijk(2,2)=1
	E_ijk(3,1)=1;  E_ijk(3,2)=2

	PI = DACOS(-1.0D0)

	DX=XL/DFLOAT(NX)
	DO I = 1,NX
		XC(I) = DFLOAT(I)*DX
	ENDDO

	DO J = 0,NY/2-1
		YC0TONY(J) = DCOS(PI*DFLOAT(J)/DFLOAT(NY))
		YC0TONY(NY-J) = -YC0TONY(J)
	ENDDO
	YC0TONY(NY/2) = 0.0D0

	YC(0:NYP)=YC0TONY( LSUBY(INDEX_MPI0):LSUBY(INDEX_MPI0)+NYP)

	DZ=ZL/DFLOAT(NZ)
	DO K = 1,NZ
		ZC(K) = DFLOAT(K)*DZ
	ENDDO

	ALPHAX = 2.0D0*PI/XL
	DO k = 0,NXSP
		ikx = LSUBX(INDEX_MPI0)+k
		WAVX(k) = DFLOAT(ikx)*ALPHAX
	ENDDO
	CWAVX = DCMPLX(0.0D0,1.0D0)*WAVX
	CWAVXRE = CWAVX/RE

	ALPHAZ = 2.0D0*PI/ZL
	DO k = -NZS,NZS
		WAVZ(K) = DFLOAT(k)*ALPHAZ
	ENDDO
	CWAVZ = DCMPLX(0.0D0,1.0D0)*WAVZ
	CWAVZRE = CWAVZ/RE

	CALL rffti (NX,trigx,itrigx)
	CALL rffti (NX2,trigx2,itrigx2)
	CALL costi (NY+1,trigy,itrigy)
	CALL costi (NY2+1,trigy2,itrigy2)
	IF(NZS>0) THEN
		CALL cffti (NZ,trigz,itrigz)
		CALL cffti (Nz2,trigz2,itrigz2)
	ENDIF

	DO iy = NYS,2,-1
		WL(iy) =  1.0D0/DFLOAT(4*iy*(iy-1))
		WD(iy) = -1.0D0/DFLOAT(2*(iy**2-1))
		WR(iy) =  1.0D0/DFLOAT(4*iy*(iy+1))
	ENDDO
	WL(2) = 2.0D0*WL(2)
	WD(NYS-1:NYS) = 0.0D0
	WR(NYS-3:NYS) = 0.0D0

	DEALLOCATE(YC0TONY)

	RETURN
	END
! ---&---1---------2---------3---------4---------5---------6---------7--
!
	SUBROUTINE xzfft (UP,U)

	USE FFT
	USE MESH
	IMPLICIT REAL*8 (A-H, O-Z)

	REAL*8 UP(NX,NZ)
	COMPLEX*16 U(0:NXSP,-NZS:NZS)

	REAL*8, ALLOCATABLE :: UMX(:),UMZ(:),UK(:,:),UX_TEMP(:)

	ALLOCATE(UMX(NX),UMZ(2*NZ),UK(0:NXSP2P1,NZ),UX_TEMP(0:NXS2P1))
	  
	DO iz = 1,NZ
		UMX = UP(1:NX,iz)
		CALL rfftf(NX,UMX,trigx,itrigx)
		UM_TEMP=0.0D0
		UX_TEMP(0)=UMX(1)/DFLOAT(NX)
		UX_TEMP(1)=0.0D0
		IF(NXS<NXH) THEN
			UX_TEMP(2:NXS2P1)=UMX(2:NXS2P1)/DFLOAT(NX)
		ELSE
			UX_TEMP(2:NX-1) = UMX(2:NX-1)/DFLOAT(NX)
			UX_TEMP(NX) = UMX(NX)/DFLOAT(NX+NX)
		ENDIF
		UK(0:NXSP2P1,iz)=UX_TEMP( LSUBX2(INDEX_MPI0): LSUBX2(INDEX_MPI0)+NXSP2P1)
	ENDDO

	U = DCMPLX( 0.0D0, 0.0D0 )
	IF(NZS>0) THEN
		DO ix = 0,NXSP
			DO iz = 1,NZ
				UMZ(2*iz-1) = UK(2*ix,  iz)
				UMZ(2*iz  ) = UK(2*ix+1,iz)
			ENDDO
			CALL cfftf(NZ,UMZ,trigz,itrigz)
			IF(NZS<NZH) THEN
				DO iz = 0,NZS
					U(ix,iz) = DCMPLX( UMZ(2*iz+1)/DFLOAT(NZ), UMZ(2*iz+2)/DFLOAT(NZ) )
				ENDDO
				DO iz = -NZS,-1
					iz1 = NZ + iz
					U(ix,iz) = DCMPLX( UMZ(2*iz1+1)/DFLOAT(NZ), UMZ(2*iz1+2)/DFLOAT(NZ) )
				ENDDO
			ELSE
				DO iz = 0,NZH-1
					U(ix,iz) = DCMPLX( UMZ(2*iz+1)/DFLOAT(NZ), UMZ(2*iz+2)/DFLOAT(NZ) )
				ENDDO
				U(ix, NZH) = DCMPLX( UMZ(NZ+1)/DFLOAT(NZ+NZ), UMZ(NZ+2)/DFLOAT(NZ+NZ) )
				U(ix,-NZH) = U(ix,NZH) 
				DO iz = -NZH+1,-1
					iz1 = NZ + iz
					U(ix,iz) = DCMPLX( UMZ(2*iz1+1)/DFLOAT(NZ), UMZ(2*iz1+2)/DFLOAT(NZ) )
				ENDDO
			ENDIF
		ENDDO
	ELSE
		DO iz = -NZS,NZS
			DO ix = 0,NXSP
				U(ix,iz) = DCMPLX(UK(2*ix,iz),UK(2*ix+1,iz))
			ENDDO
		ENDDO
	ENDIF

	DEALLOCATE(UMX,UMZ,UK)

	RETURN
	END
! ---&---1---------2---------3---------4---------5---------6---------7--
!
	SUBROUTINE channel_xyzfft (C,UP,U)

	USE FFT
	USE MESH
	IMPLICIT REAL*8 (A-H, O-Z)

	INCLUDE 'mpi_head.h'

	CHARACTER C
	REAL*8 UP(NX,0:NYP,NZ)
	COMPLEX*16 U(0:NXSP,0:NYS,-NZS:NZS)

	REAL*8, ALLOCATABLE :: UKX(:,:,:),UKZ(:,:,:),UY(:,:,:),UMX(:),UMY(:),UMZ(:)

	ALLOCATE(UKX(0:NXS2P1,0:NYP,NZ),UKZ(0:NXS2P1,0:NYP,-NZS:NZS),UY(0:NXSP2P1,0:NY,-NZS:NZS),UMX(NX),UMY(0:NY),UMZ(2*NZ))

	IF(C .EQ. 'F') GO TO 10
	IF(C .EQ. 'B') GO TO 20

10	UKX=0.0D0
	DO iy = 0,NYP
		DO iz = 1,NZ
			UMX = UP(1:NX,iy,iz)
			CALL rfftf(NX,UMX,trigx,itrigx)
			UKX(0,iy,iz) = UMX(1)/DFLOAT(NX)
			UKX(1,iy,iz) = 0.0D0
			IF(NXS<NXH) THEN
				UKX(2:NXS2P1,iy,iz)=UMX(2:NXS2P1)/DFLOAT(NX)
			ELSE
				UKX(2:NX-1,iy,iz) = UMX(2:NX-1)/DFLOAT(NX)
				UKX(NX,iy,iz) = UMX(NX)/DFLOAT(NX+NX)
			ENDIF
		ENDDO
	ENDDO

	UKZ=0.0D0
	IF(NZS>0) THEN
		DO ix = 0,NXS
			DO iy = 0,NYP
				DO iz = 1,NZ
					UMZ(2*iz-1) = UKX(2*ix  ,iy,iz)
					UMZ(2*iz  ) = UKX(2*ix+1,iy,iz)
				ENDDO
				CALL cfftf(NZ,UMZ,trigz,itrigz)
				IF(NZS<NZH) THEN
					DO iz = 0,NZS
						UKZ(2*ix  ,iy,iz) = UMZ(2*iz+1)/DFLOAT(NZ)
						UKZ(2*ix+1,iy,iz) = UMZ(2*iz+2)/DFLOAT(NZ)
					ENDDO
					DO iz = -NZS,-1
						iz1 = NZ + iz
						UKZ(2*ix  ,iy,iz) = UMZ(2*iz1+1)/DFLOAT(NZ)
						UKZ(2*ix+1,iy,iz) = UMZ(2*iz1+2)/DFLOAT(NZ)
					ENDDO
				ELSE
					DO iz = 0,NZH-1
						UKZ(2*ix  ,iy,iz) = UMZ(2*iz+1)/DFLOAT(NZ)
						UKZ(2*ix+1,iy,iz) = UMZ(2*iz+2)/DFLOAT(NZ)
					ENDDO
					UKZ(2*ix  ,iy, NZH) = UMZ(NZ+1)/DFLOAT(NZ+NZ)
					UKZ(2*ix+1,iy, NZH) = UMZ(NZ+2)/DFLOAT(NZ+NZ)
					UKZ(2*ix  ,iy,-NZH) = UKZ(2*ix  ,iy,NZH)
					UKZ(2*ix+1,iy,-NZH) = UKZ(2*ix+1,iy,NZH)
					DO iz = -NZH+1,-1
						iz1 = NZ + iz
						UKZ(2*ix  ,iy,iz) = UMZ(2*iz1+1)/DFLOAT(NZ)
						UKZ(2*ix+1,iy,iz) = UMZ(2*iz1+2)/DFLOAT(NZ)
					ENDDO
				ENDIF
			ENDDO
		ENDDO
	ELSE
		UKZ(0:NXS2P1,0:NYP,-NZS:NZS)=UKX(0:NXS2P1,0:NYP,1:NZ)
	ENDIF

	IF(INDEX_Serial==1) THEN
		UY = UKZ
	ELSE
		DO K = 0,NPROC-1
			CALL MPI_SENDRECV(UKZ(LSUBX2(K):LSUBX2(K)+NXSP2P1T(K),0:NYP,-NZS:NZS),(NXSP2P1T(K)+1)*(NYP+1)*NZS2P1,MPI_REAL8,K,10, &
							UY(0:NXSP2P1,LSUBY(K):LSUBY(K)+NYPT(K),-NZS:NZS),(NXSP2P1+1)*(NYPT(K)+1)*NZS2P1,MPI_REAL8,K,10,MPI_COMM_WORLD,STATUS,ierr)
		ENDDO

		CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	ENDIF

	DO iz = -NZS,NZS
		DO ix = 0,NXSP2P1
			UMY = UY(ix,0:NY,iz)
			CALL cost (NY+1,UMY,trigy,itrigy)
			UY(ix,0,iz) = UMY(0)/DFLOAT(NY+NY)
			UY(ix,1:NYS,iz) = UMY(1:NYS)/DFLOAT(NY)
		ENDDO
	ENDDO

	DO iy = 0,NYS
		DO iz = -NZS,NZS
			DO ix = 0,NXSP
				U(ix,iy,iz) = DCMPLX(UY(2*ix,iy,iz),UY(2*ix+1,iy,iz))
			ENDDO
		ENDDO
	ENDDO

	RETURN


20	DO iy = 0,NYS
		DO iz = -NZS,NZS
			DO ix = 0,NXSP
				UY(2*ix  ,iy,iz) = DREAL(U(ix,iy,iz))
				UY(2*ix+1,iy,iz) = DIMAG(U(ix,iy,iz))
			ENDDO
		ENDDO
	ENDDO

	DO iz = -NZS,NZS
		DO ix = 0,NXSP2P1
			UMY = 0.0D0
			UMY(0:NYS) = UY(ix,0:NYS,iz)
			UMY(0) = UMY(0)+UMY(0)
			CALL cost (NY+1,UMY,trigy,itrigy)
			UY(ix,0:NY,iz) = UMY/2.0D0
		ENDDO
	ENDDO	  

	IF(INDEX_Serial==1) THEN
		UKZ = UY
	ELSE
		DO K = 0,NPROC-1
			CALL MPI_SENDRECV(UY(0:NXSP2P1,LSUBY(K):LSUBY(K)+NYPT(K),-NZS:NZS),(NXSP2P1+1)*(NYPT(K)+1)*NZS2P1,MPI_REAL8,K,10, &
							UKZ(LSUBX2(K):LSUBX2(K)+NXSP2P1T(K),0:NYP,-NZS:NZS),(NXSP2P1T(K)+1)*(NYP+1)*NZS2P1,MPI_REAL8,K,10,MPI_COMM_WORLD,STATUS,ierr)
		ENDDO

		CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	ENDIF

	IF(NZS>0) THEN
		DO ix = 0,NXS
			DO iy = 0,NYP
				UMZ = 0.0D0
				IF(NZS<NZH) THEN
					DO iz = 0,NZS
						UMZ(2*iz+1) = UKZ(2*ix  ,iy,iz)
						UMZ(2*iz+2) = UKZ(2*ix+1,iy,iz)
					ENDDO
					DO iz = -NZS,-1
						iz1 = NZ + iz
						UMZ(2*iz1+1) = UKZ(2*ix  ,iy,iz)
						UMZ(2*iz1+2) = UKZ(2*ix+1,iy,iz)
					ENDDO
				ELSE
					DO iz = 0,NZH-1
						UMZ(2*iz+1) = UKZ(2*ix  ,iy,iz)
						UMZ(2*iz+2) = UKZ(2*ix+1,iy,iz)
					ENDDO
					UMZ(NZ+1) = UKZ(2*ix  ,iy,NZH)+UKZ(2*ix  ,iy,-NZH)
					UMZ(NZ+2) = UKZ(2*ix+1,iy,NZH)+UKZ(2*ix+1,iy,-NZH)
					DO iz = -NZH+1,-1
						iz1 = NZ + iz
						UMZ(2*iz1+1) = UKZ(2*ix  ,iy,iz)
						UMZ(2*iz1+2) = UKZ(2*ix+1,iy,iz)
					ENDDO
				ENDIF
				CALL cfftb(NZ,UMZ,trigz,itrigz)
				DO iz = 1,NZ
					UKX(2*ix  ,iy,iz) = UMZ(2*iz-1)
					UKX(2*ix+1,iy,iz) = UMZ(2*iz  )
				ENDDO
			ENDDO
		ENDDO
	ELSE
 		UKX(0:NXS2P1,0:NYP,1:NZ)=UKZ(0:NXS2P1,0:NYP,-NZS:NZS)
	ENDIF
   
	DO iy = 0,NYP
		DO iz = 1,NZ
			UMX = 0.0D0
			UMX(1) = UKX(0,iy,iz)
			IF(NXS<NXH) THEN
				UMX(2:NXS2P1)=UKX(2:NXS2P1,iy,iz)
			ELSE
				UMX(2:NX-1)=UKX(2:NX-1,iy,iz)
				UMX(NX)=UKX(NX,iy,iz)+UKX(NX,iy,iz)
			ENDIF
			CALL rfftb(NX,UMX,trigx,itrigx)
			UP(1:NX,iy,iz) = UMX
		ENDDO
	ENDDO

	DEALLOCATE(UKX,UKZ,UY,UMX,UMY,UMZ)

	RETURN
	END
! ---&---1---------2---------3---------4---------5---------6---------7--
!
	SUBROUTINE channel_xyzfft2 (C,UP2,U)

	USE FFT
	USE MESH
	IMPLICIT REAL*8 (A-H, O-Z)

	INCLUDE 'mpi_head.h'

	CHARACTER C
	REAL*8 UP2(NX2,0:NY2P,NZ2)
	COMPLEX*16 U(0:NXSP,0:NYS,-NZS:NZS)

	REAL*8, ALLOCATABLE :: UKX(:,:,:),UKZ(:,:,:),UY(:,:,:),UMX(:),UMY(:),UMZ(:)

	ALLOCATE(UKX(0:NXS2P1,0:NY2P,NZ2),UKZ(0:NXS2P1,0:NY2P,-NZS:NZS),UY(0:NXSP2P1,0:NY2,-NZS:NZS),UMX(NX2),UMY(0:NY2),UMZ(2*NZ2))

	IF(C .EQ. 'F') GO TO 10
	IF(C .EQ. 'B') GO TO 20

10	UKX=0.0D0
	DO iy = 0,NY2P
		DO iz = 1,NZ2
			UMX = UP2(1:NX2,iy,iz)
			CALL rfftf(NX2,UMX,trigx2,itrigx2)
			UKX(0,iy,iz) = UMX(1)/DFLOAT(NX2)
			UKX(1,iy,iz) = 0.0D0
			IF(NXS<NXH) THEN
				UKX(2:NXS2P1,iy,iz)=UMX(2:NXS2P1)/DFLOAT(NX2)
			ELSE
				UKX(2:NX+1,iy,iz) = UMX(2:NX+1)/DFLOAT(NX2)
			ENDIF
		ENDDO
	ENDDO

	UKZ=0.0D0
	IF(NZS>0) THEN
		DO ix = 0,NXS
			DO iy = 0,NY2P
				DO iz = 1,NZ2
					UMZ(2*iz-1) = UKX(2*ix  ,iy,iz)
					UMZ(2*iz  ) = UKX(2*ix+1,iy,iz)
				ENDDO
				CALL cfftf(NZ2,UMZ,trigz2,itrigz2)
				IF(NZS<NZH) THEN
					DO iz = 0,NZS
						UKZ(2*ix  ,iy,iz) = UMZ(2*iz+1)/DFLOAT(NZ2)
						UKZ(2*ix+1,iy,iz) = UMZ(2*iz+2)/DFLOAT(NZ2)
					ENDDO
					DO iz = -NZS,-1
						iz1 = NZ2 + iz
						UKZ(2*ix  ,iy,iz) = UMZ(2*iz1+1)/DFLOAT(NZ2)
						UKZ(2*ix+1,iy,iz) = UMZ(2*iz1+2)/DFLOAT(NZ2)
					ENDDO
				ELSE
					DO iz = 0,NZH
						UKZ(2*ix  ,iy,iz) = UMZ(2*iz+1)/DFLOAT(NZ2)
						UKZ(2*ix+1,iy,iz) = UMZ(2*iz+2)/DFLOAT(NZ2)
					ENDDO
					DO iz = -NZH,-1
						iz1 = NZ2 + iz
						UKZ(2*ix  ,iy,iz) = UMZ(2*iz1+1)/DFLOAT(NZ2)
						UKZ(2*ix+1,iy,iz) = UMZ(2*iz1+2)/DFLOAT(NZ2)
					ENDDO
				ENDIF
			ENDDO
		ENDDO
	ELSE
		UKZ(0:NXS2P1,0:NY2P,-NZS:NZS)=UKX(0:NXS2P1,0:NY2P,1:NZ2)
	ENDIF

	IF(INDEX_Serial==1) THEN
		UY = UKZ
	ELSE
		DO K = 0,NPROC-1
			CALL MPI_SENDRECV(UKZ(LSUBX2(K):LSUBX2(K)+NXSP2P1T(K),0:NY2P,-NZS:NZS),(NXSP2P1T(K)+1)*(NY2P+1)*NZS2P1,MPI_REAL8,K,10, &
							UY(0:NXSP2P1,LSUBY2(K):LSUBY2(K)+NY2PT(K),-NZS:NZS),(NXSP2P1+1 )*&
                           &(NY2PT(K)+1)*NZS2P1,MPI_REAL8,K,10,MPI_COMM_WORLD,STATUS,ierr)
		ENDDO

		CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	ENDIF

	DO iz = -NZS,NZS
		DO ix = 0,NXSP2P1
			UMY = UY(ix,0:NY2,iz)
			CALL cost (NY2+1,UMY,trigy2,itrigy2)
			UY(ix,0,iz) = UMY(0)/DFLOAT(NY2+NY2)
			UY(ix,1:NYS,iz) = UMY(1:NYS)/DFLOAT(NY2)
		ENDDO
	ENDDO

	DO iy = 0,NYS
		DO iz = -NZS,NZS
			DO ix = 0,NXSP
				U(ix,iy,iz) = DCMPLX(UY(2*ix,iy,iz),UY(2*ix+1,iy,iz))
			ENDDO
		ENDDO
	ENDDO

	RETURN


20	DO iy = 0,NYS
		DO iz = -NZS,NZS
			DO ix = 0,NXSP
				UY(2*ix  ,iy,iz) = DREAL(U(ix,iy,iz))
				UY(2*ix+1,iy,iz) = DIMAG(U(ix,iy,iz))
			ENDDO
		ENDDO
	ENDDO

	DO iz = -NZS,NZS
		DO ix = 0,NXSP2P1
			UMY = 0.0D0
			UMY(0:NYS) = UY(ix,0:NYS,iz)
			UMY(0) = UMY(0)+UMY(0)
			CALL cost (NY2+1,UMY,trigy2,itrigy2)
			UY(ix,0:NY2,iz) = UMY/2.0D0
		ENDDO
	ENDDO

	IF(INDEX_Serial==1) THEN
		UKZ = UY
	ELSE
		DO K = 0,NPROC-1
			CALL MPI_SENDRECV(UY(0:NXSP2P1,LSUBY2(K):LSUBY2(K)+NY2PT(K),-NZS:NZS),(NXSP2P1+1)*(NY2PT(K)+1)*NZS2P1,MPI_REAL8,K,10, &
							UKZ(LSUBX2(K):LSUBX2(K)+NXSP2P1T(K),0:NY2P,-NZS:NZS),(NXSP2P1T(K)+1&
                               & )*(NY2P+1)*NZS2P1,MPI_REAL8,K,10,MPI_COMM_WORLD,STATUS,ierr)
		ENDDO

		CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	ENDIF

	IF(NZS>0) THEN
		DO ix = 0,NXS
			DO iy = 0,NY2P
				UMZ = 0.0D0
				IF(NZS<NZH) THEN
					DO iz = 0,NZS
						UMZ(2*iz+1) = UKZ(2*ix  ,iy,iz)
						UMZ(2*iz+2) = UKZ(2*ix+1,iy,iz)
					ENDDO
					DO iz = -NZS,-1
						iz1 = NZ2 + iz
						UMZ(2*iz1+1) = UKZ(2*ix  ,iy,iz)
						UMZ(2*iz1+2) = UKZ(2*ix+1,iy,iz)
					ENDDO
				ELSE
					DO iz = 0,NZH
						UMZ(2*iz+1) = UKZ(2*ix  ,iy,iz)
						UMZ(2*iz+2) = UKZ(2*ix+1,iy,iz)
					ENDDO
					DO iz = -NZH,-1
						iz1 = NZ2 + iz
						UMZ(2*iz1+1) = UKZ(2*ix  ,iy,iz)
						UMZ(2*iz1+2) = UKZ(2*ix+1,iy,iz)
					ENDDO
				ENDIF
				CALL cfftb(NZ2,UMZ,trigz2,itrigz2)
				DO iz = 1,NZ2
					UKX(2*ix  ,iy,iz) = UMZ(2*iz-1)
					UKX(2*ix+1,iy,iz) = UMZ(2*iz  )
				ENDDO
			ENDDO
		ENDDO
	ELSE
		UKX(0:NXS2P1,0:NY2P,1:NZ2)=UKZ(0:NXS2P1,0:NY2P,-NZS:NZS)
	ENDIF

	DO iy = 0,NY2P
		DO iz = 1,NZ2
			UMX = 0.0D0
			UMX(1) = UKX(0,iy,iz)
			IF(NXS<NXH) THEN
				UMX(2:NXS2P1)=UKX(2:NXS2P1,iy,iz)
			ELSE
				UMX(2:NX+1)=UKX(2:NX+1,iy,iz)
			ENDIF
			CALL rfftb(NX2,UMX,trigx2,itrigx2)
			UP2(1:NX2,iy,iz) = UMX
		ENDDO
	ENDDO

	DEALLOCATE(UKX,UKZ,UY,UMX,UMY,UMZ)

	RETURN
	END
