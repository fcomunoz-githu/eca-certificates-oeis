import re,sys,json,csv
import sympy as sp
from seqs import classify,gen
from castruct import certify
n,x,k,S=sp.symbols('n x k S')

def annihilator(typ,p,forms,t0):
    """return (list of coefficient polys in forward shift as sympy poly in S, N) such that P(S)a(t)=0 for t>=N"""
    T=t0+2*p
    L=max(sum(len(w) for w in f) for f in forms)
    if typ in('dec','bin'):
        beta=2 if typ=='dec' else 10
        P=sp.Integer(1)
        # partial sums of run lengths from the right (runs at odd indices of forms)
        sums=set([0]); 
        for f in forms:
            s=0; runs=f[1::2]
            for B in reversed(runs): s+=len(B); sums.add(s)
        for s in sums: P*= (S**p-beta**s)
        return sp.Poly(sp.expand(P),S),T
    if typ in('on','off'): return sp.Poly(sp.expand((S**p-1)**2),S),T
    if typ in('ton','toff'): return sp.Poly(sp.expand((S-1)*(S**p-1)**2),S),T
    lb,ld=len(forms[0][1]),len(forms[0][3])
    per=p if lb==p else p*lb*ld
    Nmid=T+p*(L+t0+p+2)
    if typ=='mid': return sp.Poly(S**per-1,S),Nmid
    beta=10 if typ=='mbin' else 2
    return sp.Poly(sp.expand((S**per-1)*(S-beta)),S),Nmid

def clean(e):
    e=e.strip().rstrip('.').strip()
    e=e.replace('^','**')
    return e

def parse_cond(c):
    c=c.strip().rstrip('.').strip()
    lo=0;par=None
    m=re.search(r'n\s*>\s*(-?\d+)',c); 
    if m: lo=int(m.group(1))+1
    m=re.search(r'n\s*>=\s*(-?\d+)',c)
    if m: lo=int(m.group(1))
    if re.search(r'\bodd\b',c): par=1
    if re.search(r'\beven\b',c): par=0
    return lo,par

def check_rec(a,coef,lo,P,N):
    d=max(coef) if coef else 0
    dP=P.degree()
    hi=max(lo,N+d)+dP+d+2
    if hi>=len(a): return None
    for t in range(max(lo,d),hi+1):
        if a[t]!=sum(c*a[t-j] for j,c in coef.items()): return ('fail',t)
    return ('ok',hi)

def parse_rec(s):
    m=re.match(r'^a\(n\)\s*=\s*(.*?)(?:\s+for\s+(.*))?$',s)
    if not m or 'a(n-' not in m.group(1): return None
    rhs=m.group(1); cond=m.group(2) or ''
    A=sp.Function('A')
    try: ex=sp.sympify(clean(rhs).replace('a(n-','A(n-'),locals={'A':A})
    except Exception: return None
    coef={}
    for term in sp.Add.make_args(sp.expand(ex)):
        fs=term.atoms(sp.Function)
        if len(fs)!=1: return None
        f=list(fs)[0]
        if f.func.__name__!='A': return None
        j=-(f.args[0]-n)
        c=term/f
        if not (c.is_Rational and j.is_Integer): return None
        coef[int(j)]=coef.get(int(j),0)+c
    return coef,cond

def check_gf(a,off,num,den,P,N):
    num=sp.Poly(num,x); den=sp.Poly(den,x)
    assert den.eval(0)!=0 and N>=off
    dn=num.degree(); dd=den.degree()
    hi=max(dn+1,N+dd)+P.degree()+dd+2
    if hi>=len(a): return None
    ser=sp.series(num.as_expr()/den.as_expr(),x,0,hi+1).removeO()
    ser=sp.Poly(ser,x)
    for t in range(0,hi+1):
        v=ser.coeff_monomial(x**t)
        want=a[t] if t>=off else 0
        if v!=want: return ('fail',t)
    return ('ok',hi)

def exp_annihilator(expr,var):
    """expr: exponential polynomial in var. return Poly in S annihilating it (forward shift), or None"""
    expr=sp.expand(sp.powsimp(sp.expand_power_exp(sp.expand(expr)),force=True))
    roots={}
    for term in sp.Add.make_args(expr):
        lam=sp.Integer(1); rest=term
        for pw in term.atoms(sp.Pow):
            if pw.exp.has(var):
                e=sp.expand(pw.exp); a1=e.coeff(var); 
                if not (e-a1*var).free_symbols.isdisjoint({var}): return None
                lam*=pw.base**a1
                rest=rest/pw
        rest=sp.simplify(rest)
        if not rest.is_polynomial(var): return None
        deg=sp.degree(rest,var) if rest.has(var) else 0
        lam=sp.nsimplify(lam)
        roots[lam]=max(roots.get(lam,0),deg+1)
    Pq=sp.Integer(1)
    for lam,mult in roots.items(): Pq*=(S-lam)**mult
    return sp.Poly(sp.expand(Pq),S)

def subseq_annihilator(P,q):
    y=sp.Symbol('y')
    R=sp.resultant(P.as_expr(),y-S**q,S)
    return sp.Poly(sp.expand(R).subs(y,S),S)

def check_closed(a,expr_s,cond,P,N):
    lo,par=parse_cond(cond)
    e=clean(expr_s)
    try: ex=sp.sympify(e,locals={'n':n,'floor':sp.floor,'I':sp.I})
    except Exception: return None
    q=4 if ex.has(sp.I) else (2 if (par is not None or ex.has(sp.floor) or ex.has((-1)**n)) else 1)
    res=[]
    for s in range(q):
        if par is not None and (s%2)!=par: continue
        exk=sp.expand(ex.subs(n,q*k+s))
        exk=exk.replace(sp.floor,lambda z: sp.floor(sp.expand(z)))
        # floor(k + c) with c rational in [0,1): k + floor(c)
        def fl(z):
            z=sp.expand(z); c=z-z.coeff(k)*k
            if z.coeff(k).is_Integer and c.is_Rational: return z.coeff(k)*k+sp.floor(c)
            return sp.floor(z)
        exk=exk.replace(sp.floor,fl)
        if exk.has(sp.floor):
            # floor of non-integer quantity: fall back (cannot certify)
            return None
        exk=sp.expand(sp.powsimp(sp.expand(exk.subs((-1)**(2*k),1)),force=True))
        exk=sp.nsimplify(sp.expand(sp.simplify(exk)))
        Pe=exp_annihilator(exk,k)
        if Pe is None: return None
        cs=Pe.all_coeffs()[::-1]
        chk=sp.simplify(sp.expand(sp.powsimp(sp.expand(sum(c*exk.subs(k,k+i) for i,c in enumerate(cs))),force=True)))
        if chk!=0: return None
        Pq=subseq_annihilator(P,q) if q>1 else P
        # u_k=a(qk+s) satisfies Pq for qk+s>=N ; check enough k
        kmin=0
        while q*kmin+s<lo: kmin+=1
        kN=0
        while q*kN+s<N: kN+=1
        hi=max(kmin,kN)+Pq.degree()+Pe.degree()+2
        if q*hi+s>=len(a): return None
        f=sp.lambdify(k,exk,'sympy')
        for kk in range(kmin,hi+1):
            if sp.nsimplify(exk.subs(k,kk))!=a[q*kk+s]: return ('fail',q*kk+s)
        res.append(hi)
    return ('ok',max(res) if res else 0)

def statements(L):
    out=[];inblock=False;blockhdr=''
    for l in L:
        if l[:2] not in ('%F','%C'): continue
        t=l[11:]
        if re.search(r'onjecture',t) and '(Start)' in t: inblock=True; blockhdr=t; continue
        if inblock:
            if '(End)' in t: inblock=False; continue
            out.append((t,blockhdr))
        elif re.search(r'onjectur',t): out.append((t,t))
    return out

def strip_sig(t):
    t=re.sub(r'\s*-\s*_[^_]+_,.*$','',t)
    t=re.sub(r'\s*=\s*A\d{6}\([^)]*\)','',t)   # drop trailing '= Axxxxxx(..)' cross-references
    t=re.sub(r',?\s*conjectured\.?\s*$','',t)
    if re.search(r'where i = sqrt\(-1\)',t):
        t=re.sub(r'\s*where i = sqrt\(-1\)\.?','',t); t=re.sub(r'\bi\b','I',t)
    t=re.sub(r'^(Empirical\s*|Conjectured?\s*(recurrence)?:?\s*|Conjecture:\s*)','',t.strip(),flags=re.I)
    return t.strip()

def verify_statement(t,a,off,P,N):
    s=strip_sig(t)
    m=re.match(r'^(?:Empirical\s+)?G\.f\.:\s*(.*)$',s,re.I)
    if m:
        g=clean(m.group(1))
        try:
            ge=sp.sympify(g,locals={'x':x}); num,den=sp.fraction(sp.together(ge))
        except Exception: return 'gf',None
        return 'gf',check_gf(a,off,sp.expand(num),sp.expand(den),P,N)
    r=parse_rec(s)
    if r:
        coef,cond=r; lo,par=parse_cond(cond)
        if par is not None: return 'rec',None
        return 'rec',check_rec(a,coef,max(lo,off),P,N)
    m=re.match(r'^a\(n\)\s*=\s*(.*?)(?:\s*(?:,|;)?\s*for\s+(.*))?$',s)
    if m and 'a(' not in m.group(1) and 'A' not in m.group(1):
        cond=m.group(2) or ''
        lo,par=parse_cond(cond)
        if lo<off: cond=cond+f' n>={off}'
        return 'closed',check_closed(a,m.group(1),cond if lo>=off else f'n>={off} '+cond,P,N)
    return 'other',None
