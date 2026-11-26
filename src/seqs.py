import re
from ca import rows
TYPES=[('Decimal representation of the n-th iteration','dec'),('Binary representation of the n-th iteration','bin'),
('Number of ON (black) cells in the n-th iteration','on'),('Number of OFF (white) cells in the n-th iteration','off'),
('Total number of ON (black) cells after n iterations','ton'),('Total number of OFF (white) cells after n iterations','toff'),
('Binary representation of the middle column','mbin'),('Decimal representation of the middle column','mdec'),('Middle column','mid')]
def classify(name):
    m=re.search(r'"Rule (\d+)"',name); rule=int(m.group(1))
    for k,v in TYPES:
        if name.startswith(k): return rule,v
    return rule,None
def gen(rule,typ,N):
    R=rows(rule,N+1)
    out=[];acc=0;mc=0
    for t in range(N+1):
        W=R[t][0][1:-1]   # window [-t,t]
        if typ=='dec': out.append(int(''.join(map(str,W)),2))
        elif typ=='bin': out.append(int(''.join(map(str,W))))
        elif typ=='on': out.append(sum(W))
        elif typ=='off': out.append(len(W)-sum(W))
        elif typ=='ton': acc+=sum(W); out.append(acc)
        elif typ=='toff': acc+=len(W)-sum(W); out.append(acc)
        elif typ in('mbin','mdec','mid'):
            c=W[t]
            if typ=='mid': out.append(c)
            else: mc=mc*(10 if typ=='mbin' else 2)+c; out.append(mc)
    return out
