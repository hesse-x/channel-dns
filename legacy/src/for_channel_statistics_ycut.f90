! ---&---1---------2---------3---------4---------5---------6---------7--
!
	SUBROUTINE FOR_STATISTICS (Vp,Prep)

	USE MESH
	USE STATISTICS
	IMPLICIT REAL*8 (A-H,O-Z)

	character*2,external :: itochar2

	REAL*8 Vp(3,NX,0:NYP,NZ),PREp(NX,0:NYP,NZ),Utemp(3)

	REAL*8, ALLOCATABLE :: VVp(:,:,:,:)

	ALLOCATE( VVP(3,NX,0:NYP,NZ) )

    VVp(1:3,1:NX,0:NYP,1:NZ) = Vp(1:3,1:NX,0:NYP,1:NZ)*Vp(1:3,1:NX,0:NYP,1:NZ)
    PREp = PREp-(VVp(1,1:NX,0:NYP,1:NZ)+VVp(2,1:NX,0:NYP,1:NZ)+VVp(3,1:NX,0:NYP,1:NZ))/2.0D0

	IF(INDEX_MPI0==0) THEN
		OPEN(10,FILE='./Result/0/DPDX_DPDZ.dat',STATUS='UNKNOWN')
			DPDX_mean=(DPDX_mean*DFLOAT(ICOUNT)+DPDX)/DFLOAT(ICOUNT+1)
			DPDZ_mean=(DPDZ_mean*DFLOAT(ICOUNT)+DPDZ)/DFLOAT(ICOUNT+1)
			WRITE(10,*) ICOUNT+1
			WRITE(10,*) DPDX_mean,DPDZ_mean
		CLOSE(10)
	ENDIF

    OPEN(11,FILE='./Result/1/V_Pre_mean_Field_'//itochar2(my_id)//'.dat',STATUS='UNKNOWN')
		DO iy=0,NYP
			Utemp=0.0D0;	SPre=0.0D0
			DO iz=1,NZ
				DO ix=1,NX
					DO I=1,3
						Utemp(I)=Utemp(I)+Vp(I,ix,iy,iz)
					ENDDO
					SPre=SPre+PREp(ix,iy,iz)
				ENDDO
			ENDDO
			V_mean(1:3,iy)=(V_mean(1:3,iy)*DFLOAT(ICOUNT)+Utemp/DFLOAT(NXTNZ))/DFLOAT(ICOUNT+1)
			Pre_mean(iy)=(Pre_mean(iy)*DFLOAT(ICOUNT)+SPre/DFLOAT(NXTNZ))/DFLOAT(ICOUNT+1)
			WRITE(11,'(5(1X,F18.8))') YC(iy),(V_mean(L,iy),L=1,3),Pre_mean(iy)
! 			DO iz=0,NZ-1
! 				DO ix=0,NX-1
! 					Vp(1:3,ix,iy,iz)=Vp(1:3,ix,iy,iz)-V_mean(1:3,iy)
! 					Prep(ix,iy,iz)=Prep(ix,iy,iz)-Pre_mean(iy)
! 				ENDDO
! 			ENDDO
		ENDDO
	CLOSE(11)

!	Vp = Vp - V_mean
!	PREp = PREp - Pre_mean

!	V_mean = (V_mean*DFLOAT(ICOUNT)+Vp)/DFLOAT(ICOUNT+1)
!	Pre_mean = (Pre_mean*DFLOAT(ICOUNT)+PREp)/DFLOAT(ICOUNT+1)

!	Vp = Vp - V_mean
!	PREp = PREp - Pre_mean

! 	VVp(1:3,0:NX-1,LY0:NYP,0:NZ-1) = Vp(1:3,0:NX-1,LY0:NYP,0:NZ-1)*Vp(1:3,0:NX-1,LY0:NYP,0:NZ-1)

	OPEN(12,FILE='./Result/2/V_Pre_moment2_Field_'//itochar2(my_id)//'.dat',STATUS='UNKNOWN')
	OPEN(13,FILE='./Result/3/V_Pre_moment3_Field_'//itochar2(my_id)//'.dat',STATUS='UNKNOWN')
	OPEN(14,FILE='./Result/4/V_Pre_moment4_Field_'//itochar2(my_id)//'.dat',STATUS='UNKNOWN')
		DO iy=0,NYP

!-----uu,vv,ww,pp
			Utemp=0.0D0;	SPre=0.0D0
			DO iz=1,NZ
				DO ix=1,NX
					DO I=1,3
						Utemp(I)=Utemp(I)+VVp(I,ix,iy,iz)
					ENDDO
					SPre=SPre+PREp(ix,iy,iz)*PREp(ix,iy,iz)
				ENDDO
			ENDDO
			VV_mean(1:3,iy)=(VV_mean(1:3,iy)*DFLOAT(ICOUNT)+Utemp/DFLOAT(NXTNZ))/DFLOAT(ICOUNT+1)
! 			VV_mean(1:3,iy)=VV_mean(1:3,iy)-V_mean(1:3,iy)*V_mean(1:3,iy)
			Pre_Pre_mean(iy)=(Pre_Pre_mean(iy)*DFLOAT(ICOUNT)+SPre/DFLOAT(NXTNZ))/DFLOAT(ICOUNT+1)
! 			Pre_Pre_mean(iy)=Pre_Pre_mean(iy)-Pre_mean(iy)*Pre_mean(iy)

!-----uv,vw,wv
			Utemp=0.0D0;	SPre=0.0D0
			DO I=1,3
				I1=I+1-I/3*3
				DO iz=1,NZ
					DO ix=1,NX
						Utemp(I)=Utemp(I)+Vp(I,ix,iy,iz)*Vp(I1,ix,iy,iz)
					ENDDO
				ENDDO
				VV_mean(I+3,iy)=(VV_mean(I+3,iy)*DFLOAT(ICOUNT)+Utemp(I)/DFLOAT(NXTNZ))/DFLOAT(ICOUNT+1)
! 				VV_mean(I+3,iy)=VV_mean(I+3,iy)-V_mean(I,iy)*V_mean(I1,iy)
			ENDDO

!-----pu,pv,pw
			Utemp=0.0D0;	SPre=0.0D0
			DO iz=1,NZ
				DO ix=1,NX
					DO I=1,3
						Utemp(I)=Utemp(I)+Prep(ix,iy,iz)*Vp(I,ix,iy,iz)
					ENDDO
				ENDDO
			ENDDO
			Pre_V_mean(1:3,iy) = (Pre_V_mean(1:3,iy)*DFLOAT(ICOUNT)+Utemp/DFLOAT(NXTNZ))/DFLOAT(ICOUNT+1)
! 			Pre_V_mean(1:3,iy) = Pre_V_mean(1:3,iy)-Pre_mean(iy)*V_mean(1:3,iy)

			WRITE(12,'(11(1X,F18.8))') YC(iy),(VV_mean(L1,iy),L1=1,6),Pre_Pre_mean(iy),(Pre_V_mean(L2,iy),L2=1,3)

!-----uuu,vvv,www
			Utemp=0.0D0;	SPre=0.0D0
			DO iz=1,NZ
				DO ix=1,NX
					DO I=1,3
						Utemp(I)=Utemp(I)+VVp(I,ix,iy,iz)*Vp(I,ix,iy,iz)
					ENDDO
				ENDDO
			ENDDO
			VVV_mean(1:3,iy) = (VVV_mean(1:3,iy)*DFLOAT(ICOUNT)+Utemp/DFLOAT(NXTNZ))/DFLOAT(ICOUNT+1)
! 			VVV_mean(1:3,iy) = VVV_mean(1:3,iy)-3.0D0*V_mean(1:3,iy)*VV_mean(1:3,iy)-V_mean(1:3,iy)**3

			WRITE(13,'(4(1X,F18.8))') YC(iy),(VVV_mean(L,iy),L=1,3)

!-----uuuu,vvvv,wwww
			Utemp=0.0D0;	SPre=0.0D0
			DO iz=1,NZ
				DO ix=1,NX
					DO I=1,3
						Utemp(I)=Utemp(I)+VVp(I,ix,iy,iz)*VVp(I,ix,iy,iz)
					ENDDO
				ENDDO
			ENDDO
			VVVV_mean(1:3,iy) = (VVVV_mean(1:3,iy)*DFLOAT(ICOUNT)+Utemp/DFLOAT(NXTNZ))/DFLOAT(ICOUNT+1)
! 			VVVV_mean(1:3,iy) = VVVV_mean(1:3,iy)-6.0D0*(V_mean(1:3,iy)**2)*VV_mean(1:3,iy)-4.0D0*V_mean(1:3,iy)*VVV_mean(1:3,iy)-V_mean(1:3,iy)**4

			WRITE(14,'(4(1X,F18.8))') YC(iy),(VVVV_mean(L,iy),L=1,3)

		ENDDO
	CLOSE(12)
	CLOSE(13)
	CLOSE(14)

	ICOUNT=ICOUNT+1

	DEALLOCATE(VVp)

	END
! ---&---1---------2---------3---------4---------5---------6---------7--
!
	  SUBROUTINE TEST (uavq,Ruv)

! 	  INCLUDE 'dim.h'

! 	  REAL*8 uavq(0:NY),Ruv(0:NY)
! 
! 	  REAL*8, ALLOCATABLE :: Duavq(:)
! 
!       ALLOCATE(Duavq(0:NY))
! 
! 	  CALL cost (NY+1,uavq,trigy)
! 	  uavq(0) = uavq(0)/2.0D0
! 	  uavq(ny) = uavq(ny)/2.0D0
! 	  uavq = uavq/DFLOAT(NY)
!  
! 	  Duavq(ny) = 0.0D0
! 	  Duavq(ny-1) = 2.0D0*DFLOAT(ny)*uavq(ny)
! 	  DO iky = NY-2,0,-1
! 		 Duavq(iky) = Duavq(iky+2) + 2.0D0*DFLOAT(iky+1)*uavq(iky+1)
! 	  ENDDO
! 	  Duavq(0) = Duavq(0)/2.0D0
! 
! 	  Duavq(0) = 2.0D0*Duavq(0)
! 	  Duavq(ny) = 2.0D0*Duavq(ny)
! 	  CALL cost (NY+1,Duavq,trigy)
! 	  Duavq = Duavq/2.0D0
! 
! 	  OPEN(66,file="test.dat",status="unknown")
! 	  DO iy=NY,NY/2,-1
! 		 WRITE(66,100) 1.0D0+YC(iy),Duavq(iy)/Re,ABS(Ruv(iy)),Duavq(iy)/Re-Ruv(iy)
! 	  ENDDO
! 	  CLOSE(66)
! 
! 	  DEALLOCATE(Duavq)
! 
!  100  format(2x,4(f14.8,2x))

	  END