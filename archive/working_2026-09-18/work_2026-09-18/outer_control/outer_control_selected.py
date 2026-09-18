exec(open('/mnt/data/ssz_work/outer_control_grid.py').read().split('rows=[]')[0])
for Hfac,a1fac in [(1,1),(1.2,.2),(1.2,.5),(1.2,1),(1.5,.1),(1.5,.2),(1.5,.5),(2,.05),(2,.1),(2,.2),(3,.05),(3,.1),(4,.05)]:
 o=build(Hfac=Hfac,a1fac=a1fac)
 m=metrics(o,Ls=(6,))
 print(Hfac,a1fac,m,flush=True)
