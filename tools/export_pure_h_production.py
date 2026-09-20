#!/usr/bin/env python3
from pathlib import Path
from ssz_p5.config import repo_root
from ssz_p5.production.pure_h_export import export_pure_h

def main():
    root=repo_root(); out=root/'data/generated/pure_h'
    w,c=export_pure_h(root,out)
    print('WEAK_DIRECT_41 =',w['WEAK_DIRECT_41'])
    print('CORE_DIRECT_41 =',c['CORE_DIRECT_41'])
    print('weak rows =',w['rows'],'max scaled regression =',w['max_regression_scaled_error'])
    print('core rows =',c['rows'],'max scaled regression =',c['max_regression_scaled_error'])
    return 0 if w['WEAK_DIRECT_41']=='PASS' and c['CORE_DIRECT_41']=='PASS' else 2
if __name__=='__main__': raise SystemExit(main())
