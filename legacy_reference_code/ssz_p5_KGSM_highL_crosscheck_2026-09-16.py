from pathlib import Path
import importlib.util, numpy as np, pandas as pd
B=Path('/mnt/data')
# reducer
spec=importlib.util.spec_from_file_location('op',B/'ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py')
op=importlib.util.module_from_spec(spec); spec.loader.exec_module(op)
d=op.prepare_central(); r=d.x.to_numpy(float); u=d.u.to_numpy(float)
trusted=(u>=0.612)&(u<=0.708)
D=lambda y: op.deriv(r,np.asarray(y,float),1,9,8)
h=d.h.to_numpy(float); Ap=d.A0prime.to_numpy(float); ph=d.phiprime.to_numpy(float)
a4=d.a4.to_numpy(float); a6=d.a6.to_numpy(float); a9=d.a9.to_numpy(float)
v1=d.v1.to_numpy(float); v4=d.v4.to_numpy(float); v5=d.v5.to_numpy(float); v6=d.v6.to_numpy(float); v8=d.v8.to_numpy(float); v9=d.v9.to_numpy(float); v13=d.v13.to_numpy(float)
c4=d.c4.to_numpy(float); c5=d.c5.to_numpy(float); d2=d.d2.to_numpy(float); d3=d.d3.to_numpy(float); e4=d.e4.to_numpy(float)
z1=(a4+r*a9-.5*r*Ap*v6)*v6
m1minus=-r*a4*v8+z1
m2=2*c5*v1+(Ap*v1+.5*ph*v4)*v6
m3plus=a4*D(v1)+.5*Ap*v1*v6
m4=2*Ap*v1+ph*v4
m5minus=a4*v13-a6*v6
M11=-r*r*v6*v6/(4*v1)
M13=-r*v6/2
M33=-v1
M22=e4+2*h*c4*a6/a4-m5minus*m5minus/(4*a4*a4*v9)
M12=(-r*d3/2-h*c4*z1/(a4*v6)-r*h*a6*m2/(2*a4*v1)+r*v5*v6/(4*v1)-m1minus*m5minus/(4*a4*a4*v9)+.25*D(2*r*h*c4+r*(2*d2*v1-v4*v6)/(2*v1)+r*v6*m5minus/(2*a4*v9)))
M23=(h*Ap*c4*v1/a4+v5/2+m3plus*m5minus/(2*a4*a4*v9)-h*a6*m4/(2*a4)-.25*D(v4+v1*m5minus/(a4*v9)))
T={'11':M11,'12':M12,'13':M13,'22/L':M22,'23':M23,'33':M33}
def rel(a,b): return np.abs(a-b)/np.maximum(1.0,np.abs(b))
rows=[]
for L in (1e3,1e4,1e5,1e6,1e7):
    a=op.canonical_audit(d,float(L)); M=a['M']; S=a['S']
    candidates={
      '11':M[:,0,0], '12':M[:,0,1], '13':M[:,0,2],
      '22/L':M[:,1,1]/L, '23':M[:,1,2], '33':M[:,2,2]
    }
    row={'L':L,'S_asym_res':np.max(np.abs(S+np.swapaxes(S,1,2)))}
    for k,x in candidates.items():
      # compare project M and opposite sign, retain better for diagnosis
      eplus=rel(x,T[k]); eminus=rel(-x,T[k])
      row[k+'_proj_med']=float(np.median(eplus[trusted])); row[k+'_proj_max']=float(np.max(eplus[trusted]))
      row[k+'_tilde_med']=float(np.median(eminus[trusted])); row[k+'_tilde_max']=float(np.max(eminus[trusted]))
    rows.append(row)
out=pd.DataFrame(rows)
out.to_csv(B/'ssz_p5_KGSM_highL_crosscheck_2026-09-16.csv',index=False)
print(out.to_string(index=False))
