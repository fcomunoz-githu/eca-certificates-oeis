"""Statistics on the certificate search (Section 'The certificate search'): timing, alternative certificates,
behaviour with enlarged limits. Usage: python src/search_stats.py  (about 45 minutes)"""
import os,sys,json,time
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from ca import rows,find_structure,prove_structure
from castruct import certify
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
K=json.load(open(os.path.join(ROOT,'results','certificates.json')))
out={}
# 1. time of the default search per rule
tt={}
for ru in K:
    t=time.time(); certify(int(ru)); tt[ru]=time.time()-t
out['time_default_total_s']=sum(tt.values()); out['time_default_max_s']=max(tt.values())
out['time_uncertified_s']={ru:tt[ru] for ru,v in K.items() if v is None}
# 2. all valid certificates with p<=12, t0<40 for certified rules
alt={}
for ru,v in K.items():
    if v is None: continue
    R=rows(int(ru),40+12*8+5); found=[]
    for p in range(1,13):
        for t0 in range(0,40):
            for lB in range(1,2*p):
                F=find_structure(int(ru),R,p,t0,lB=lB)
                if F and prove_structure(int(ru),R,p,t0,F)[0]: found.append((p,t0,lB))
    alt[ru]=found
out['n_valid_certificates_per_rule']={ru:len(f) for ru,f in alt.items()}
out['min_p_per_rule']={ru:min(f)[0] for ru,f in alt.items()}
out['valid_p_values_per_rule']={ru:sorted({x[0] for x in f}) for ru,f in alt.items()}
out['max_minimal_t0']=max(min(x[1] for x in f if x[0]==min(f)[0]) for f in alt.values())
# 3. enlarged search for the 5 uncertified rules
big={}
for ru,v in K.items():
    if v is not None: continue
    t=time.time(); c=certify(int(ru),maxp=24,maxt0=80); big[ru]=dict(found=c is not None,seconds=time.time()-t)
out['enlarged_p24_t080']=big
json.dump(out,open(os.path.join(ROOT,'results','search_stats.json'),'w'),indent=1)
print(json.dumps({k:v for k,v in out.items() if k not in('n_valid_certificates_per_rule','min_p_per_rule','valid_p_values_per_rule')},indent=1))
