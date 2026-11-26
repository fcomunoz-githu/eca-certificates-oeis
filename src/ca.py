# Elementary CA from a single ON cell, exact infinite-background evolution.
# Row t is represented by W_t = c_t restricted to [-t-1, t+1]  (length 2t+3).
def f(rule,a,b,c): return (rule>>(4*a+2*b+c))&1
def step(rule,W,b):
    """W covers [-t-1,t+1], background b outside. returns (W', b')"""
    nb=f(rule,b,b,b)
    P=[b]+W+[b]
    inner=[f(rule,P[i-1],P[i],P[i+1]) for i in range(1,len(P)-1)]
    return [nb]+inner+[nb], nb
def rows(rule,T):
    W=[0,1,0]; b=0; out=[(W,b)]
    for t in range(T):
        W,b=step(rule,W,b); out.append((W,b))
    return out

def find_structure(rule,Wlist,p,t0,mlist=range(2,7),lB=None):
    """For each residue r in 0..p-1 find (A,B,C,D,E) with W_{t0+r+pm}=A B^m C D^m E for m in mlist,
       |B|=|D|=p, |A|,|E|>=1. Return list of forms or None."""
    forms=[]
    for r in range(p):
        Ws={m:Wlist[t0+r+p*m][0] for m in mlist}
        m2=min(mlist); L2=len(Ws[m2]); ok=None
        lb=p if lB is None else lB; ld=2*p-lb
        for a in range(1,L2):
            for e in range(1,L2-a-m2*(lb+ld)+1):
                good=True
                W=Ws[m2]
                A=W[:a];B=W[a:a+lb];E=W[len(W)-e:];D=W[len(W)-e-ld:len(W)-e]
                C=W[a+m2*lb:len(W)-e-m2*ld]
                for m in mlist:
                    if A+B*m+C+D*m+E!=Ws[m]: good=False;break
                if good: ok=(A,B,C,D,E);break
            if ok: break
        if not ok: return None
        forms.append(ok)
    return forms

def insert(W,pos,X): return W[:pos]+X+W[pos:]

def prove_structure(rule,Wlist,p,t0,forms):
    """Rigorous inductive check (see paper, Lemma 2/Theorem 3).
    Hypothesis H(r,m): W_{t0+r+pm} = A_r B_r^m C_r D_r^m E_r, for all m>=2.
    Base: H(0,2) checked directly. Step: H(r,m) => H(r+1,m) (r<p-1) or H(0,m+1) (r=p-1), for all m>=2."""
    # base
    A,B,C,D,E=forms[0]
    if Wlist[t0+2*p][0]!=A+B*2+C+D*2+E: return False,'base'
    bgs=[Wlist[t0+r][1] for r in range(p)]
    # background must be periodic with period p (b_t determined by rule; check b at t0+r+pm equals for all m)
    # b_t is eventually periodic with period <=2; require consistency:
    for r in range(p):
        for m in range(2,6):
            if Wlist[t0+r+p*m][1]!=Wlist[t0+r+p*2][1]: return False,'bg'
    for r in range(p):
        A,B,C,D,E=forms[r]
        b=Wlist[t0+r+2*p][1]
        x2=A+B*2+C+D*2+E
        Y,nb=step(rule,x2,b)
        # lemma: image of x_m = Y with (m-2) copies of Bt inserted at position P1 and Dt at position P2 (from right)
        # where the insertion in x is between copy1 and copy2 of B: index a+p in x -> index a+p+1 in Y (shift by 1)
        lb,ld=len(B),len(D)
        BBB=B*3; Bt=[f(rule,BBB[lb+i-1],BBB[lb+i],BBB[lb+i+1]) for i in range(lb)]
        DDD=D*3; Dt=[f(rule,DDD[ld+i-1],DDD[ld+i],DDD[ld+i+1]) for i in range(ld)]
        P1=len(A)+lb+1            # insertion index in Y
        Q2=len(E)+ld+1            # insertion index counted from right end of Y
        if r<p-1: r2,delta=r+1,0
        else: r2,delta=0,1
        A2,B2,C2,D2,E2=forms[r2]
        T2=A2+B2*(2+delta)+C2+D2*(2+delta)+E2
        if Y!=T2: return False,('img',r)
        # check background consistency of next row
        if nb!=Wlist[t0+r2+p*(2+delta)][1]: return False,'nb'
        # inserted word must be the rotation of B2 at the insertion offset, inside the run
        if len(B2)!=lb or len(D2)!=ld: return False,'len'
        o=P1-len(A2)
        if not (0<=o<=(2+delta)*lb): return False,('posB',r)
        if Bt!=(B2[o%lb:]+B2[:o%lb]): return False,('rotB',r)
        o2=Q2-len(E2)    # offset from right end of D-run
        if not (0<=o2<=(2+delta)*ld): return False,('posD',r)
        # inserting at distance o2 from the right end of run D2^k: word inserted must be rotation so that
        # suffix structure preserved: D2^k = ... ; inserting X at position (k p - o2) from left of run -> X = rot(D2,(-o2)%p)
        s_=(-o2)%ld
        if Dt!=(D2[s_:]+D2[:s_]): return False,('rotD',r)
        # left and right insertion points must not cross: P1 <= len(Y)-Q2
        if P1>len(Y)-Q2: return False,'cross'
    return True,'ok'
