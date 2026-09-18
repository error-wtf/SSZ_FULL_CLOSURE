from outer_model import *
import numpy as np
for j in [8,9,10,11]:
 print('\nPARAM',j,flush=True)
 for a in [-2,-1,-.5,-.25,0,.1,.25,.5,1,2]:
  x=np.zeros(12);x[j]=a
  vals,meta,*_=mins(x,[6,42,1000])
  print(f'{a:+.3f}',*(f'{v:+.6g}' for v in vals),*(f'u{m[2]:.6f}' for m in meta),flush=True)
