"""Recheck all retained PASS certificates and all three exact witnesses."""
from pathlib import Path
from fractions import Fraction as F
import argparse,json,time
import verify as v
HERE=Path(__file__).resolve().parent
def main():
 parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=HERE/'V17_MPFR_RECHECK.json');args=parser.parse_args()
 start=time.perf_counter();certificates=[]
 for path in sorted((HERE.parent/'certificates').glob('*.json')):
  record=json.loads(path.read_text())
  if record['status']!='PASS_COMPUTATIONAL':continue
  result=v.check_record(record);result['file']=path.name;result['source_sha256']=v.sha(path)
  certificates.append(result)
  print(path.name,result['status'],result['floor_display'],flush=True)
 path=HERE.parent/'witness_search/FIXED_RUNG_WITNESSES.json';search=json.loads(path.read_text());witnesses=[]
 v.require(search['immutable_input_sha256']['model_inputs.json']==v.INPUT_SHA,'witness model hash')
 for record in search['witnesses']:
  v.require(record['status']=='PASS_WITNESS','incomplete witness')
  result=v.witness(record['raw_baseline_normalized_inputs'],record['q_exact'])
  v.require(F(result['actual_linf_distance'])==F(record['actual_linf_change'])<=F(record['radius_limit']),'wrong witness distance')
  v.require(result['new_winner']==record['winner'] and record['rival_is_unique_winner'] is True,'wrong witness winner')
  raw=list(map(F,record['raw_baseline_normalized_inputs']));pairmax=max(raw[i]**4+raw[i+1]**4 for i in range(0,40,2))
  v.require(record['admissible'] is True and pairmax==F(record['maximum_canonical_pair_sum'])<=1,'wrong admissibility summary')
  # The integer intervals are much narrower than outward binary64 conversion.
  # Verify their reported values are compatible, while the independent integer
  # replay remains the authority for exact equality of those endpoint records.
  for stored,computed in zip(record['scores_exact'],result['scores']):
   sl,su=map(F,stored);cl,cu=map(F,computed)
   v.require(cl<=sl<=su<=cu,'integer/MPFR score enclosure incompatibility')
  witnesses.append(result);print('witness q='+result['q'],result['status'],result['actual_linf_distance'],flush=True)
 v.require({w['q'] for w in witnesses}=={'4','8','16'} and len(witnesses)==3,'expected three distinct fixed-rung witnesses')
 result={'status':'PASS','mpfr_bits':256,'certificate_count':len(certificates),
         'certificates':certificates,'witness_count':len(witnesses),'witnesses':witnesses,
         'witness_source_sha256':v.sha(path),'input_model_sha256':v.INPUT_SHA,
         'elapsed_seconds':time.perf_counter()-start,
         'scope':'Recomputed source-box certificate margins and rational point witnesses; no optimal-radius assertion.'}
 args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
 print('FINAL',result['status'],'certificates',len(certificates),'witnesses',len(witnesses))
if __name__=='__main__':main()
