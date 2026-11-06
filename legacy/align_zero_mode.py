#!/usr/bin/env python3
"""Isolate the one known zero-mode change; leave the preserved source untouched."""
import pathlib
import sys
source=pathlib.Path(sys.argv[1]).read_bytes().decode('latin1').replace('\r\n','\n')
old='IF(INDEX_MPI0==0) Div_Vc_BOUND(I,0,J,0)=BC_PRE(0,0,0)'
new='''IF(INDEX_MPI0==0) THEN
        Div_Vc_BOUND(I,0,J,0)=DCMPLX(0.0D0,0.0D0)
        IF(J==0) Div_Vc_BOUND(I,0,J,0)=BC_PRE(0,0,0)
      ENDIF'''
assert source.count(old)==1
pathlib.Path(sys.argv[2]).write_bytes(source.replace(old,new).encode('latin1'))
