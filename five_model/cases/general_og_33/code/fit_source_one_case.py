#!/usr/bin/env python3
"""Reuse the validated source-branch five-model implementation for one gate case."""
from __future__ import annotations
import argparse,csv,json,sys
from pathlib import Path
ap=argparse.ArgumentParser();ap.add_argument('--case',type=Path,required=True);ap.add_argument('--code-dir',type=Path,required=True);a=ap.parse_args()
sys.path.insert(0,str(a.code_dir.resolve()))
from run_source_branch_length_validation import MODELS,fit_task,ard_task,write_tsv
case=a.case.resolve();fits=[fit_task((str(case),*m)) for m in MODELS]
minimum=min(float(r['AIC']) for r in fits)
for r in fits:r['delta_AIC']=float(r['AIC'])-minimum
fits.sort(key=lambda r:r['model']);out=case/'results/source_branch_length_validation';out.mkdir(parents=True,exist_ok=True)
write_tsv(out/'five_model_fit_GTDB_source_branch_lengths.tsv',fits)
ard=next(r for r in fits if r['model']=='ARD_unrestricted');ardrow=ard_task((str(case),ard));write_tsv(out/'ARD_contribution_GTDB_source_branch_lengths.tsv',[ardrow])
(out/'COMPLETE.marker').write_text('status=complete\ncases=1\nbranch_scale=source_branch_lengths\n')
print(json.dumps({'case_id':ard['case_id'],'status':'complete'}))
