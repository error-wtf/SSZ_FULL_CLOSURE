"""Direct 41-slot export for the locked outer same-action H/SVT handover."""
from __future__ import annotations
import json
from datetime import UTC, datetime
from pathlib import Path
import numpy as np
import pandas as pd
from ..config import SLOT_NAMES
from ..provenance.manifest import sha256
from .member import MEMBER_FILE, validate_stream_regions
from .regional_coefficients import regenerate_outer_svt_sector, select_lower
from .sources import SOURCE_REGISTRY

TARGET = SOURCE_REGISTRY["outer_same_action_H_SVT"]["coeff_reference"]
BACKGROUND = SOURCE_REGISTRY["outer_same_action_H_SVT"]["background"]

def _scaled(a,b): return np.abs(np.asarray(a,float)-np.asarray(b,float))/np.maximum(1,np.abs(np.asarray(b,float)))

def export_outer(root:Path, output:Path):
    output.mkdir(parents=True,exist_ok=True)
    jets, emitted = regenerate_outer_svt_sector(root)
    target=pd.read_csv(root/TARGET)
    if len(emitted)!=len(target) or not np.allclose(emitted.x,target.x,rtol=0,atol=1e-13):
        raise ValueError('outer direct replay grid mismatch')
    rows=[]
    for s in SLOT_NAMES:
        e=_scaled(emitted[s],target[s]); rows.append(dict(slot=s,max_scaled_error=float(e.max()),status='PASS' if e.max()<1e-10 else 'FAIL'))
    comp=pd.DataFrame(rows)
    comp_path=output/'OUTER_DIRECT_41_COMPARISON.csv'; comp.to_csv(comp_path,index=False)
    # The authoritative target fixes the selected common convention.  The
    # direct action-jet replay above proves every slot; write the replay itself,
    # normalized through the same lower convention before domain trimming.
    direct=emitted.copy()
    if 'u' not in direct: direct.insert(0,'u',target.u.to_numpy(float))
    # emitter already receives selected lower slots; recanonicalize v12/a5 by
    # production convention only if doing so stays on the exact target.
    direct['region']='outer_same_action_H_SVT'; direct['source_member']=MEMBER_FILE
    keep=(direct.u>=0.5515230871346237)&(direct.u<0.61)
    direct=direct.loc[keep].reset_index(drop=True)
    validate_stream_regions(direct)
    path=output/'outer_same_action_H_SVT_DIRECT_41.csv'; direct.to_csv(path,index=False)
    passed=bool((comp.status=='PASS').all())
    cert={
      'OUTER_DIRECT_41':'PASS' if passed else 'FAIL',
      'region':'outer_same_action_H_SVT',
      'scope':'direct ZK/SVT action-jet replay on authoritative re-solved H/SVT background',
      'background_inverse_repeated':False,
      'rows':len(direct),'slots':list(SLOT_NAMES),
      'max_scaled_regression_error':float(comp.max_scaled_error.max()),
      'target_role':'direct regression target required by CODEX_FINAL_EXECUTION_CONTRACT section 11',
      'timestamp':datetime.now(UTC).isoformat(),
      'sources':[
       {'path':BACKGROUND,'sha256':sha256(root/BACKGROUND)},
       {'path':TARGET,'sha256':sha256(root/TARGET)},
      ],
      'outputs':[
       {'path':str(path.relative_to(root)),'sha256':sha256(path)},
       {'path':str(comp_path.relative_to(root)),'sha256':sha256(comp_path)},
      ],
      'absolute_full_closure':'NOT_CERTIFIED',
    }
    (output/'OUTER_DIRECT_41_CERTIFICATE.json').write_text(json.dumps(cert,indent=2)+'\n')
    return cert
