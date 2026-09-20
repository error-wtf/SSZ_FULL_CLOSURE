#!/usr/bin/env python3
from ssz_p5.config import repo_root
from ssz_p5.production.outer_export import export_outer

def main():
 c=export_outer(repo_root(),repo_root()/'data/generated/outer')
 print('OUTER_DIRECT_41 =',c['OUTER_DIRECT_41'])
 print('rows =',c['rows'],'max scaled regression =',c['max_scaled_regression_error'])
 return 0 if c['OUTER_DIRECT_41']=='PASS' else 2
if __name__=='__main__': raise SystemExit(main())
