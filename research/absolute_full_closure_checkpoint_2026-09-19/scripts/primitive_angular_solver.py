from pathlib import Path
import importlib.util, numpy as np, pandas as pd
from scipy.interpolate import CubicSpline
ROOT=Path('/mnt/data/ssz_final_work/SSZ_FULL_CLOSURE')
base=pd.read_csv(ROOT/'data/generated/absolute_attempt_2026-09-19/ELECTRIC_HYBRID_ONSHELL_CENTRAL_41.csv').reset_index(drop=True)
r=base.x.to_numpy(float); u=base.u.to_numpy(float); n=len(base); order=np.argsort(r); rs=r[order]
def dr(y):
 y=np.asarray(y,float); sp=CubicSpline(rs,y[order]); z=np.empty_like(y); z[order]=sp(rs,1); return z
# slots ordering same as mh_general_primitives.SLOTS
slots=[*[f'a{i}' for i in range(1,10)],*[f'b{i}' for i in range(1,6)],*[f'c{i}' for i in range(1,7)],*[f'd{i}' for i in range(1,5)],*[f'e{i}' for i in range(1,5)],*[f'v{i}' for i in range(1,14)]]
slot_index={s:i for i,s in enumerate(slots)}
z=np.load('/mnt/data/primitive_basis_responses.npz',allow_pickle=True); R=z['R']; B=z['B']; centers=z['centers']; widths=z['widths']
prims=['a1','c2','c4','F_tensor','G_tensor','H_tensor']; nb=len(centers); nv=R.shape[0]
baseS=np.stack([base[s].to_numpy(float) for s in slots],axis=0)
# geometry fixed
f=base.f.to_numpy(float); h=base.h.to_numpy(float); ph=base.phiprime.to_numpy(float)
ep=np.array([1.,0,0]); ef=np.array([0,1.,0]); ev=np.array([0,0,1.]); ephi=np.r_[ef,np.zeros(3)]; ephip=np.r_[np.zeros(3),ef]; eV=np.r_[ev,np.zeros(3)]; eVp=np.r_[np.zeros(3),ev]
TRUST=(u>=.622)&(u<.70)

def candidate(p):
 p=np.asarray(p,float); S=baseS+np.tensordot(p,R,axes=(0,0)); d=base[['u','x','phi','f','h','phiprime','A0prime','X']].copy()
 for i,s in enumerate(slots): d[s]=S[i]
 return d,S

def sym(a,b): return .5*(np.einsum('ni,nj->nij',a,b)+np.einsum('ni,nj->nij',b,a))

def direct_mass(S,L):
 q={s:S[slot_index[s]] for s in slots}
 # algebraic H0 map
 if np.any(np.abs(q['a3'])<1e-12) or np.any(np.abs(q['v1'])<1e-12): return None,None
 pp=q['a1']/q['a3']; qq=L*q['a4']/q['a3']; ppr=dr(pp); qqr=dr(qq)
 Braw=q['a2']-q['v2']*q['v4']/(2*q['v1'])
 Cphi=q['a5']+L*q['a6']-q['v2']*q['v5']/(2*q['v1'])
 CH2=q['a7']+L*q['a8']-q['v2']*q['v3']/(2*q['v1'])
 Ch1=L*(q['a9']-q['v2']*q['v6']/(2*q['v1']))
 Beff=Braw-q['a3']*ppr-CH2*pp; D=Ch1-q['a3']*qqr-CH2*qq
 if np.any(~np.isfinite(D)) or np.min(np.abs(D[TRUST]))<1e-9: return None,D
 h0v=np.c_[-CH2/D,-Cphi/D,-q['v2']/D]; h1v=np.c_[-q['a3']/D,-Beff/D,np.zeros(n)]
 H0=ep[None,:]-qq[:,None]*h0v; H1=-qq[:,None]*h1v-pp[:,None]*ef[None,:]
 H=np.concatenate([H0,H1],axis=1); hh=np.concatenate([h0v,h1v],axis=1)
 EPHI=np.tile(ephi,(n,1)); EPHIP=np.tile(ephip,(n,1)); EV=np.tile(eV,(n,1)); EVP=np.tile(eVp,(n,1))
 C=np.zeros((n,6,6))
 def add(co,a,b):
  nonlocal C; C += np.asarray(co)[:,None,None]*sym(a,b)
 add(q['c2'],H,EPHIP); add(q['c3']+L*q['c4'],H,EPHI); add(L*q['c5'],H,hh); add(q['c6'],H,H)
 add(L*q['d2'],hh,EPHIP); add(L*q['d3'],hh,EPHI); add(L*q['d4'],hh,hh); add(q['e2'],EPHIP,EPHIP); add(q['e3']+L*q['e4'],EPHI,EPHI)
 S0=q['v3'][:,None]*H+q['v4'][:,None]*EPHIP+q['v5'][:,None]*EPHI+L*q['v6'][:,None]*hh
 add(-q['v1'],EV,EV); add(np.ones(n),EV,S0); add(-1/(4*q['v1']),S0,S0)
 v1p=dr(q['v1']); J0=-2*v1p[:,None]*EV-2*q['v1'][:,None]*EVP+L*(q['v8'][:,None]*hh+q['v12'][:,None]*H+q['v13'][:,None]*EPHI)
 add(-1/(4*L*q['v9']),J0,J0)
 A=C[:,:3,:3]; BB=C[:,:3,3:]; Bs=.5*(BB+np.swapaxes(BB,1,2)); Bp=np.empty_like(Bs)
 for i in range(3):
  for j in range(3): Bp[:,i,j]=dr(Bs[:,i,j])
 return A-Bp,D

def angular(p,L1=1e7,L2=2e7):
 d,S=candidate(p); A,D1=direct_mass(S,L1); Bm,D2=direct_mass(S,L2)
 bad=np.full(n,-1e30)
 if A is None or Bm is None: return {k:bad.copy() for k in ['cminus','cplus','c55','c56','disc','cV']},d,S
 q={s:S[slot_index[s]] for s in slots}
 M11_0=-r*r*q['v6']**2/(4*q['v1']); M13_0=-r*q['v6']/2; M33_0=-q['v1']
 M11_1=2*L1*(A[:,0,0]-Bm[:,0,0]); M13_1=2*L1*(A[:,0,2]-Bm[:,0,2]); M33_1=2*L1*(A[:,2,2]-Bm[:,2,2])
 M12_0=2*Bm[:,0,1]-A[:,0,1]; M23_0=2*Bm[:,1,2]-A[:,1,2]; M22_0=2*Bm[:,1,1]/L2-A[:,1,1]/L1
 # Reconstruct alpha7 from the current emitted 41-slot candidate.  This is
 # the common Appendix-A identity inverted from a6 and avoids freezing the
 # odd-sector tensor coefficient while Horndeski primitives are varied.
 fp=dr(f); hp=dr(h); a4=q['a4']; a4p=dr(a4); Ap=base.A0prime.to_numpy(float)
 geom=(r*(fp*h+f*hp)-2*f*h)/(2*r*f*h*h*ph)
 alpha7=ph/(2*r*f)*(q['a6']-geom*a4+a4p/(h*ph)+Ap*q['v6']/(4*h*ph))
 if np.any(~np.isfinite(alpha7)) or np.min(alpha7[TRUST])<=1e-12:
  return {k:bad.copy() for k in ['cminus','cplus','c55','c56','disc','cV']},d,S
 Mc1=(M11_0*M33_1-2*M13_0*M13_1+M33_0*M11_1)/(4*r*r*f*h*h*alpha7*M33_0)
 Mc2=ph*(r*h*ph*M22_0*M33_0+2*M13_0*M23_0-2*M12_0*M33_0)/(4*r*f*h*q['v1']*alpha7*M33_0)
 K22=q['e1']; a4=q['a4']; v10=q['v10']; K1=-(2*r*r*a4/f)/(1-f*q['v6']**2/(8*a4*v10)); K33=f*q['v1']**2*K1/(2*r*r*a4*v10); cV=-r*r*M33_0/(f*K33)
 c55=-(Mc1+Mc2+r*r*M22_0/(f*K22)); cross=r*r*(4*f*alpha7*Mc2+ph*ph*M22_0)**2/(4*f*f*ph*ph*alpha7*K22); c56=r*r*(Mc1+Mc2)*M22_0/(f*K22)-cross
 disc=c55*c55-4*c56; sq=np.sqrt(np.maximum(disc,0)); cm=.5*(c55-sq); cp=.5*(c55+sq)
 return {'cminus':cm,'cplus':cp,'c55':c55,'c56':c56,'disc':disc,'cV':cV,'alpha7':alpha7,'D1':D1,'D2':D2},d,S

if __name__=='__main__':
 a,_,_=angular(np.zeros(nv));
 for k in ['cminus','cplus','c55','c56','disc','cV']:
  j=np.flatnonzero(TRUST)[np.nanargmin(a[k][TRUST])]; print(k,float(a[k][j]),'u',float(u[j]))


def angular_from_S(S,L1=1e7,L2=2e7):
 """Evaluate the action-derived high-L angular gates for an explicit 41-slot array."""
 S=np.asarray(S,float)
 d=base[['u','x','phi','f','h','phiprime','A0prime','X']].copy()
 for i,sn in enumerate(slots): d[sn]=S[i]
 A,D1=direct_mass(S,L1); Bm,D2=direct_mass(S,L2)
 bad=np.full(n,-1e30)
 if A is None or Bm is None: return {k:bad.copy() for k in ['cminus','cplus','c55','c56','disc','cV','alpha7']},d
 q={sn:S[slot_index[sn]] for sn in slots}
 M11_0=-r*r*q['v6']**2/(4*q['v1']); M13_0=-r*q['v6']/2; M33_0=-q['v1']
 M11_1=2*L1*(A[:,0,0]-Bm[:,0,0]); M13_1=2*L1*(A[:,0,2]-Bm[:,0,2]); M33_1=2*L1*(A[:,2,2]-Bm[:,2,2])
 M12_0=2*Bm[:,0,1]-A[:,0,1]; M23_0=2*Bm[:,1,2]-A[:,1,2]; M22_0=2*Bm[:,1,1]/L2-A[:,1,1]/L1
 fp=dr(f); hp=dr(h); a4=q['a4']; a4p=dr(a4); Ap=base.A0prime.to_numpy(float)
 geom=(r*(fp*h+f*hp)-2*f*h)/(2*r*f*h*h*ph)
 alpha7=ph/(2*r*f)*(q['a6']-geom*a4+a4p/(h*ph)+Ap*q['v6']/(4*h*ph))
 if np.any(~np.isfinite(alpha7)) or np.min(alpha7[TRUST])<=1e-12:
  return {k:bad.copy() for k in ['cminus','cplus','c55','c56','disc','cV','alpha7']},d
 Mc1=(M11_0*M33_1-2*M13_0*M13_1+M33_0*M11_1)/(4*r*r*f*h*h*alpha7*M33_0)
 Mc2=ph*(r*h*ph*M22_0*M33_0+2*M13_0*M23_0-2*M12_0*M33_0)/(4*r*f*h*q['v1']*alpha7*M33_0)
 K22=q['e1']; v10=q['v10']; K1=-(2*r*r*a4/f)/(1-f*q['v6']**2/(8*a4*v10)); K33=f*q['v1']**2*K1/(2*r*r*a4*v10); cV=-r*r*M33_0/(f*K33)
 c55=-(Mc1+Mc2+r*r*M22_0/(f*K22)); cross=r*r*(4*f*alpha7*Mc2+ph*ph*M22_0)**2/(4*f*f*ph*ph*alpha7*K22); c56=r*r*(Mc1+Mc2)*M22_0/(f*K22)-cross
 disc=c55*c55-4*c56; sq=np.sqrt(np.maximum(disc,0)); cm=.5*(c55-sq); cp=.5*(c55+sq)
 return {'cminus':cm,'cplus':cp,'c55':c55,'c56':c56,'disc':disc,'cV':cV,'alpha7':alpha7,'D1':D1,'D2':D2},d
