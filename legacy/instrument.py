#!/usr/bin/env python3
"""Add wall-clock instrumentation to the copied driver, without changing kernels."""
import pathlib
import sys

source = pathlib.Path(sys.argv[1]).read_bytes().decode('latin1').replace('\r\n', '\n')
source = source.replace('  REAL*8  a(0:1),b(0:1)', '''  REAL*8  a(0:1),b(0:1)
  REAL*8 BM_START,BM_COMPUTE,BM_MAX,BM_LOOP,BM_TOTAL
  INTEGER BM_WARMUP,BM_ARGI,BM_IOS
  CHARACTER(128) BM_ARG''', 1)
source = source.replace('  CALL CPU_TIME(time1)', '''  BM_WARMUP=0
  BM_ARGI=1
  DO WHILE(BM_ARGI<=COMMAND_ARGUMENT_COUNT())
    CALL GET_COMMAND_ARGUMENT(BM_ARGI,BM_ARG)
    IF(TRIM(BM_ARG)/='--benchmark-warmup') CALL MPI_ABORT(MPI_COMM_WORLD,2,ierr)
    BM_ARGI=BM_ARGI+1
    CALL GET_COMMAND_ARGUMENT(BM_ARGI,BM_ARG)
    READ(BM_ARG,*,IOSTAT=BM_IOS) BM_WARMUP
    IF(BM_IOS/=0) CALL MPI_ABORT(MPI_COMM_WORLD,2,ierr)
    IF(BM_WARMUP<0) CALL MPI_ABORT(MPI_COMM_WORLD,2,ierr)
    BM_ARGI=BM_ARGI+1
  ENDDO
  CALL CPU_TIME(time1)''', 1)
# The imported diagnostics are temporary dumps, outside the physics calculation.
a=source.index('! [TEMP 2026-08-14]')
b=source.index('  IF(INDEX_CONT==0) THEN',a)
source=source[:a]+source[b:]
source=source.replace("    IF(INDEX_MPI0==0) WRITE(*,*) 'it = ', it,'    TIME = ', TIME",'')
source=source.replace('  DO 29 it=1,NSTEPS', '''  BM_COMPUTE=0.0D0
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  BM_LOOP=MPI_WTIME()
  DO 29 it=1,NSTEPS
    IF(it==BM_WARMUP+1) CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    BM_START=MPI_WTIME()''', 1)
source=source.replace('    IF(MOD(it,NDUMFRQ)==0) THEN', '''    IF(it>BM_WARMUP) BM_COMPUTE=BM_COMPUTE+MPI_WTIME()-BM_START
    IF(MOD(it,NDUMFRQ)==0) THEN''', 1)
source=source.replace('29  CONTINUE','''29  CONTINUE
  BM_TOTAL=MPI_WTIME()-BM_LOOP
  CALL MPI_REDUCE(BM_COMPUTE,BM_MAX,1,MPI_REAL8,MPI_MAX,0,MPI_COMM_WORLD,ierr)
  IF(my_id==0) WRITE(*,'(a,es24.16,a)') 'benchmark compute wall time = ',BM_MAX,' s'
  CALL MPI_REDUCE(BM_TOTAL,BM_MAX,1,MPI_REAL8,MPI_MAX,0,MPI_COMM_WORLD,ierr)
  IF(my_id==0) THEN
    WRITE(*,'(a,es24.16,a)') 'benchmark loop wall time = ',BM_MAX,' s'
    WRITE(*,'(a,i0)') 'benchmark measured steps = ',MAX(0,NSTEPS-BM_WARMUP)
  ENDIF''',1)
pathlib.Path(sys.argv[2]).write_bytes(source.encode('latin1'))
