from outer_model import *
import numpy as np
for j in [2,5,6,8,9,10]:
 print('\nPARAM',j,flush=True)
 amps=np.linspace(-1.5,1.5,9) if j<8 else np.linspace(-2,2,9)
 for a in amps:
  x=np.zeros(12);x[j]=a
  try:
   vals,meta,*_=mins(x,[6,42,1000])
   print(f'{a:+.3f}',*(f'{v:+.6g}' for v in vals),*(f'u{m[2]:.6f}' for m in meta),flush=True)
  except Exception as e: print('ERR',a,e,flush=True)
