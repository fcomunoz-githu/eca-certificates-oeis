from ca import *
def certify(rule,maxp=12,maxt0=40):
    R=rows(rule,maxt0+maxp*8+5)
    for p in range(1,maxp+1):
        for t0 in range(0,maxt0):
            for lB in [p]+[l for l in range(1,2*p) if l!=p]:
                F=find_structure(rule,R,p,t0,lB=lB)
                if F:
                    ok,why=prove_structure(rule,R,p,t0,F)
                    if ok: return p,t0,F
    return None
