"""Directed Decimal80 fixed-source winner cover; no sampling uncertainty."""
from decimal import Decimal as D, localcontext, ROUND_FLOOR, ROUND_CEILING
from fractions import Fraction as F
from pathlib import Path
from datetime import datetime, timezone
import hashlib, json

BASE=Path(__file__).resolve().parent
PREC=80

def calc(fn, rounding):
    with localcontext() as c:
        c.prec=PREC; c.rounding=rounding
        return fn()

def exact_fraction(f):
    f=F(f)
    a,b=D(f.numerator),D(f.denominator)
    return (calc(lambda:a/b,ROUND_FLOOR),calc(lambda:a/b,ROUND_CEILING))

ZERO=(D(0),D(0)); ONE=(D(1),D(1))

def add(a,b):
    return calc(lambda:a[0]+b[0],ROUND_FLOOR),calc(lambda:a[1]+b[1],ROUND_CEILING)

def neg(a):
    return a[1].copy_negate(),a[0].copy_negate()

def sub(a,b):
    return add(a,neg(b))

def mul(a,b):
    if a[0]>=0 and b[0]>=0:
        return calc(lambda:a[0]*b[0],ROUND_FLOOR),calc(lambda:a[1]*b[1],ROUND_CEILING)
    if a[0]>=0 and b[1]<=0:
        return calc(lambda:a[1]*b[0],ROUND_FLOOR),calc(lambda:a[0]*b[1],ROUND_CEILING)
    if a[1]<=0 and b[0]>=0:
        return mul(b,a)
    return (min(calc(lambda:x*y,ROUND_FLOOR) for x in a for y in b),
            max(calc(lambda:x*y,ROUND_CEILING) for x in a for y in b))

def divide(a,b):
    assert b[0]>0
    inv=(calc(lambda:D(1)/b[1],ROUND_FLOOR),calc(lambda:D(1)/b[0],ROUND_CEILING))
    return mul(a,inv)

def transc(a,kind):
    # Decimal ln/exp are correctly rounded to nearest, even in directed contexts.
    # Widen each independently evaluated monotone endpoint by one representable.
    with localcontext() as c:
        c.prec=PREC
        lo=getattr(a[0],kind)().next_minus(c)
        hi=getattr(a[1],kind)().next_plus(c)
    return lo,hi

def total(items):
    items=list(items)
    return (calc(lambda:sum((a[0] for a in items),D(0)),ROUND_FLOOR),
            calc(lambda:sum((a[1] for a in items),D(0)),ROUND_CEILING))

def encode(a):
    return [str(a[0]),str(a[1])]

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def prepare(rows,keys):
    selected=[r for r in rows if all('excluded' not in r['criteria'][k] for k in keys)]
    assert len(selected)==83
    atoms=[[] for _ in keys]
    scorelogs=[]
    for row in selected:
        one_minus_mu=[]; nus=[]
        for j,k in enumerate(keys):
            item=row['criteria'][k]
            mu,nu=F(item['mu']),F(item['nu'])
            assert 0<mu<1 and 0<nu<1 and mu+nu<=1
            lm=transc(exact_fraction(mu),'ln'); ln=transc(exact_fraction(nu),'ln')
            atoms[j].extend([neg(lm),neg(ln)])
            one_minus_mu.append(transc(exact_fraction(1-mu),'ln'))
            nus.append(ln)
        scorelogs.append((one_minus_mu,nus))
    winner=next(i for i,r in enumerate(selected) if r['id']=='3102')
    return selected,atoms,scorelogs,winner

def evaluate(lo,hi,prepared):
    rows,atoms,scorelogs,winner=prepared
    p=(exact_fraction(lo)[0],exact_fraction(hi)[1])
    e=[]
    for col in atoms:
        vals=[]
        for rate in col:
            y=mul(rate,p)
            vals.append(mul(y,transc(neg(y),'exp')))
        e.append(sub(ONE,divide(total(vals),exact_fraction(len(rows)))))
    if any(x[0]<=0 for x in e):
        return None
    denom=total(e)
    weights=[divide(x,denom) for x in e]
    scores=[]
    for logs_a,logs_b in scorelogs:
        a=transc(total(mul(w,l) for w,l in zip(weights,logs_a)),'exp')
        b=transc(total(mul(w,l) for w,l in zip(weights,logs_b)),'exp')
        scores.append(sub(sub(ONE,a),b))
    rivals=[(i,sub(scores[winner],s)[0]) for i,s in enumerate(scores) if i!=winner]
    worst,gap=min(rivals,key=lambda t:t[1])
    return {'p_lo':str(lo),'p_hi':str(hi),'gap_lower':str(gap),
            'worst_rival_id':rows[worst]['id'],'weights':[encode(w) for w in weights],
            'winner_score':encode(scores[winner]),
            'rival_score_upper':{rows[i]['id']:str(scores[i][1]) for i,_ in rivals}}

def run():
    ip=BASE/'root_results_R1/EXACT_INPUTS_R1.json'
    proto=BASE/'PROTOCOL_FREEZE_R1.json'
    freeze=BASE/'CONTINUUM_VERIFICATION_FREEZE_R1.md'
    data=json.loads(ip.read_text(encoding='utf-8'))
    spec=json.loads(proto.read_text(encoding='utf-8'))
    dest=BASE/'continuum_decimal_R1.json'
    assert not dest.exists()
    result={'started_utc':datetime.now(timezone.utc).isoformat(),'code_sha256':sha(__file__),
            'source_input_sha256':sha(ip),'protocol_sha256':sha(proto),'verification_freeze_sha256':sha(freeze),
            'precision_decimal_digits':PREC,'parameter':'p=1/q in[1/16,1], so q in[1,16]',
            'claim':'fixed-source winner only, no sampling or input-uncertainty certificate','groups':{}}
    for group,keys in [('primary',spec['primary_criteria']),('secondary',spec['secondary_criteria'])]:
        prepared=prepare(data['rows'],keys)
        queue=[(F(1,16),F(1),0)]; leaves=[]; unresolved=[]; processed=0
        while queue:
            lo,hi,depth=queue.pop()
            record=evaluate(lo,hi,prepared); processed+=1
            if record is not None and D(record['gap_lower'])>0:
                record['depth']=depth; leaves.append(record)
            elif depth>=16 or len(queue)+len(leaves)+len(unresolved)+2>4096:
                unresolved.append({'p_lo':str(lo),'p_hi':str(hi),'depth':depth,'bounds':record})
            else:
                mid=(lo+hi)/2
                queue.extend([(mid,hi,depth+1),(lo,mid,depth+1)])
        leaves.sort(key=lambda x:F(x['p_lo']))
        intervals=sorted([(F(x['p_lo']),F(x['p_hi'])) for x in leaves+unresolved])
        assert intervals[0][0]==F(1,16) and intervals[-1][1]==1
        assert all(a[1]==b[0] for a,b in zip(intervals,intervals[1:]))
        result['groups'][group]={'status':'PASS' if not unresolved else 'UNRESOLVED',
            'criteria':keys,'winner_id':'3102','number_of_rivals':82,'processed_nodes':processed,
            'accepted_leaves':len(leaves),'unresolved_leaves':len(unresolved),
            'minimum_certified_gap_lower':str(min(D(x['gap_lower']) for x in leaves)) if leaves else None,
            'leaves':leaves,'unresolved':unresolved}
        print(group,result['groups'][group]['status'],len(leaves),'leaves',flush=True)
    result['completed_utc']=datetime.now(timezone.utc).isoformat()
    with dest.open('x',encoding='utf-8') as h:
        json.dump(result,h,ensure_ascii=False,indent=2); h.write('\n')
    print('receipt',sha(dest),flush=True)

if __name__=='__main__':
    run()
