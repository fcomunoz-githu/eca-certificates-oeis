"""Reproduce every result of the paper from the 260 OEIS entries in data/oeis_seq.zip.
Usage (from the repository root):  python src/run_all.py
Outputs: results/certificates.json, results/claims.json, results/summary.json,
         results/sanity_annihilators.txt, results/tables/*.tex
"""
import os,sys,re,json,collections
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0,os.path.join(ROOT,'src'))
import sympy as sp
from verify import *            # annihilator, statements, strip_sig, verify_statement, parse_rec, check_rec, check_closed, S
from seqs import classify,gen
from castruct import certify
from ca import rows
from data_setup import ensure_data
ensure_data()

SEQ=os.path.join(ROOT,'data','oeis_seq'); OUT=os.path.join(ROOT,'results')
os.makedirs(os.path.join(OUT,'tables'),exist_ok=True)
ids=[l.strip() for l in open(os.path.join(ROOT,'data','eca_entries.txt')) if l.strip()]
w2s=lambda w:''.join(map(str,w))

def read(i):
    L=open(os.path.join(SEQ,i+'.seq'),encoding='utf-8').read().splitlines()
    d={}
    for l in L: d.setdefault(l[:2],[]).append(l[11:])
    data=[int(x) for x in ''.join(d.get('%S',[])+d.get('%T',[])+d.get('%U',[])).replace(' ','').strip(',').split(',')]
    return L,d['%N'][0],int(d['%O'][0].split(',')[0]),data

# 1. simulation reproduces OEIS data
mism=[]
for i in ids:
    L,name,off,data=read(i); rule,typ=classify(name)
    if gen(rule,typ,off+len(data))[off:off+len(data)]!=data: mism.append(i)
print('entries:',len(ids),' data mismatches:',len(mism))

# 2. certificates
certs={}; claims=[]
for i in ids:
    L,name,off,data=read(i); rule,typ=classify(name)
    if rule not in certs: certs[rule]=certify(rule)
    c=certs[rule]
    for t,hdr in statements(L):
        auth=re.findall(r'_([^_]+)_,\s*([A-Z][a-z]{2} \d{1,2} \d{4})',t) or re.findall(r'_([^_]+)_,\s*([A-Z][a-z]{2} \d{1,2} \d{4})',hdr)
        rec=dict(id=i,rule=rule,type=typ,offset=off,statement=t,author=auth[0][0] if auth else '',date=auth[0][1] if auth else '')
        if c is None: rec['status']='no-certificate'; claims.append(rec); continue
        p,t0,F=c; P,N=annihilator(typ,p,F,t0)
        a=gen(rule,typ,N+300)
        st=[]
        for q in re.split(r';\s*(?=a\(n\))',strip_sig(t)):
            try: kind,res=verify_statement(q,a,off,P,N)
            except Exception: kind,res='error',None
            st.append(dict(part=q,kind=kind,result=None if res is None else res[0],index=None if res is None else res[1]))
        if all(s['result']=='ok' for s in st): rec['status']='proved'
        elif any(s['result']=='fail' for s in st): rec['status']='false-as-stated'
        else: rec['status']='not-covered'
        rec.update(parts=st,annihilator=str(P.as_expr()),N=N)
        claims.append(rec)

# 3. corrections for false statements
for rec in claims:
    if rec['status']!='false-as-stated': continue
    c=certs[rec['rule']]; P,N=annihilator(rec['type'],c[0],c[2],c[1]); a=gen(rec['rule'],rec['type'],N+300)
    s=strip_sig(rec['statement']); pr=parse_rec(s); out=None
    if pr:
        coef,_=pr
        for lo in range(0,80):
            r=check_rec(a,coef,lo,P,N)
            if r and r[0]=='ok': out=f'holds exactly for n>={lo}'; break
        if out is None: out='does not hold for any range n>=n0'
    else:
        for lo in range(0,30):
            ok=True
            for q in re.split(r';\s*(?=a\(n\))',s):
                m=re.match(r'^a\(n\)\s*=\s*(.*?)(?:\s*(?:,|;)?\s*for\s+(.*))?$',q)
                cond=re.sub(r'n\s*>=?\s*-?\d+','',m.group(2) or '')+f' n>={lo}'
                r=check_closed(a,m.group(1),cond,P,N)
                if not r or r[0]!='ok': ok=False;break
            if ok: out=f'holds for n>={lo}'; break
    rec['correction']=out

certjson={str(k):(None if v is None else dict(p=v[0],t0=v[1],lB=len(v[2][0][1]),lD=len(v[2][0][3]),
          forms=[dict(r=r,A=w2s(f[0]),B=w2s(f[1]),C=w2s(f[2]),D=w2s(f[3]),E=w2s(f[4])) for r,f in enumerate(v[2])])) for k,v in sorted(certs.items())}
json.dump(certjson,open(os.path.join(OUT,'certificates.json'),'w'),indent=1)
json.dump(claims,open(os.path.join(OUT,'claims.json'),'w'),indent=1)
st=collections.Counter(c['status'] for c in claims)
summ=dict(entries=len(ids),data_mismatches=len(mism),rules=len(certs),certified_rules=sum(v is not None for v in certs.values()),
          statements=len(claims),status=st,entries_with_a_proof=len({c['id'] for c in claims if c['status']=='proved'}),
          entries_fully_proved=len([i for i in ids if all(c['status']=='proved' for c in claims if c['id']==i)]),
          by_type={t:dict(collections.Counter(c['status'] for c in claims if c['type']==t)) for t in sorted({c['type'] for c in claims})})
json.dump(summ,open(os.path.join(OUT,'summary.json'),'w'),indent=1)
print(json.dumps(summ,indent=1))

# 4. independent numerical sanity check of Proposition 5 (not part of the proofs)
viol=0; tested=0
for ru,c in certs.items():
    if c is None: continue
    p,t0,F=c
    for typ in ['dec','bin','on','off','ton','toff','mid','mbin','mdec']:
        P,N=annihilator(typ,p,F,t0); cs=P.all_coeffs()[::-1]; a=gen(ru,typ,N+len(cs)+300)
        for t in range(N,N+300):
            tested+=1; viol+= sum(x*a[t+i] for i,x in enumerate(cs))!=0
    R=rows(ru,400)
    for t in range(t0+2*p,401):
        r=(t-t0)%p; m=(t-t0-r)//p; A,B,C,D,E=F[r]
        if R[t][0]!=A+B*m+C+D*m+E: viol+=1
open(os.path.join(OUT,'sanity_annihilators.txt'),'w').write(f'annihilator tests: {tested}, row-shape tests up to t=400, violations: {viol}\n')
print('sanity violations:',viol)

# 5. Remark 6: with beta^{l_B} instead of beta^{l_D} the annihilator fails for most asymmetric certificates
asym=0; fails=0
for ru,c in certs.items():
    if c is None: continue
    p,t0,F=c; lb,ld=len(F[0][1]),len(F[0][3])
    if lb==ld: continue
    for typ,beta in (('dec',2),('bin',10)):
        asym+=1; _,N=annihilator(typ,p,F,t0)
        Q=sp.Poly(sp.expand((S**p-1)*(S**p-beta**lb)*(S**p-beta**(2*p))),S); qs=Q.all_coeffs()[::-1]
        a=gen(ru,typ,N+len(qs)+60)
        fails+=any(sum(x*a[t+i] for i,x in enumerate(qs))!=0 for t in range(N,N+50))
with open(os.path.join(OUT,'sanity_annihilators.txt'),'a') as f:
    f.write(f'asymmetric certificates (rule,base) pairs: {asym}; beta^lB variant fails in: {fails}\n')
print('remark check:',asym,fails)
