import ast,sys
from pathlib import Path
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from src.pruning_policy import pruning_ratio,has_observed_tree_content
assert pruning_ratio(1,0)==1e8
assert pruning_ratio(0,0)==0
assert pruning_ratio(3,1)==2
assert not has_observed_tree_content(pd.DataFrame({'ROOT':[1]}))
assert not has_observed_tree_content(pd.DataFrame({'M1':[0]}))
assert has_observed_tree_content(pd.DataFrame({'M1':[1]}))
s=ast.parse((ROOT/'src/run_phylosilid_fullTree_scRNA.py').read_text())
f=next(x for x in s.body if isinstance(x,ast.FunctionDef) and x.name=='select_best_cv')
ns={'pd':pd,'compute_pruning_ratio':pruning_ratio,'PRUNING_CONFIDENCE_THRESHOLD':10}
exec(compile(ast.Module(body=[f],type_ignores=[]),'<selection>','exec'),ns)
def row(cv,pre,final,stage='final'):
 return dict(cv_value=cv,omega_pre_qc=pre,omega_final=final,scaffold_count=6,retention_rate=.9,final_cells=105,final_muts=24,output_stage=stage)
for rows,expected in [([row(.9,1,0),row(.3,1,0)],.3),([row(.3,3,1),row(.8,2,1)],.8),([row(.3,0,0),row(.9,1,0)],.3),([row(.3,0,0,'pre_fallback'),row(.9,1,0)],.9),([row(.3,2,1,'pre_fallback'),row(.9,1,2,'pre_fallback')],.9),([row(.3,1,2,'pre_fallback'),row(.9,1,0,'pre_fallback')],.3)]:
 result=ns['select_best_cv'](pd.DataFrame(rows));assert result[1]==expected,(result[1],expected)
for name in ['run_phylosilid_fullTree_scRNA.py','full_tree_builder.py','pruning_policy.py']:
 ast.parse((ROOT/'src'/name).read_text())
print('PASS: zero/zero, positive/zero, ordinary ratio, empty content, CV ties and final priority')
# Execute the actual orchestration function with controlled QC outputs. This
# exercises snapshot restoration after in-place mutation, rather than a copy
# of the fallback implementation.
import copy,json,logging,os,tempfile
mod=ast.parse((ROOT/'src/full_tree_builder.py').read_text())
f=next(x for x in mod.body if isinstance(x,ast.FunctionDef) and x.name=='build_fully_resolved_tree')
f.returns=None
for a in list(f.args.args)+list(f.args.kwonlyargs):a.annotation=None
fns={'copy':copy,'json':json,'logging':logging,'os':os,'has_observed_tree_content':has_observed_tree_content}
initial=pd.DataFrame({'M1':[1]},index=['cell1'])
def init(**kw):return {'name':'ROOT'},initial.copy(),[],[]
def refined(**kw):return kw['T_current'],kw['M_current'],[],[],0.0,[]
def empty_qc(**kw):
 kw['M_current'].drop(index=kw['M_current'].index,inplace=True)
 return kw['T_current'],kw['M_current'],['cell1'],[],['M1'],[],[],[]
fns.update(build_initial_tree=init,build_refined_tree=refined,tree_QC_and_filter=empty_qc)
exec(compile(ast.Module(body=[f],type_ignores=[]),'<builder>','exec'),fns)
with tempfile.TemporaryDirectory() as temp:
 params={}
 result=fns['build_fully_resolved_tree'](None,initial,['M1'],initial,initial,initial,params,temp,'sample',[])
 assert result[1].equals(initial) and result[4]==0 and params['_output_stage']=='pre_fallback'
 assert params['_pre_qc_snapshot'][1].equals(initial)
 assert json.loads(Path(temp,'output_stage.json').read_text())['output_stage']=='pre_fallback'
 assert result[5]==[] and result[7]==[]
 fns['build_initial_tree']=lambda **kw: ({'name':'ROOT'},pd.DataFrame(),[],[])
 try:fns['build_fully_resolved_tree'](None,initial,[],initial,initial,initial,{},temp,'sample',[])
 except ValueError as e:assert 'No usable' in str(e)
 else:raise AssertionError('Both empty must reject')
print('PASS: real builder restores nonempty zero-omega pre snapshot; rejects both empty')
