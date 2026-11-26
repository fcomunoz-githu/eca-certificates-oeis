"""Assemble paper/main_full.tex from paper/main.tex and results/tables/*.tex, then run pdflatex twice."""
import os,subprocess
H=os.path.dirname(os.path.abspath(__file__)); T=os.path.join(H,'..','results','tables')
s=open(os.path.join(H,'main.tex')).read()
import json
st=os.path.join(H,'..','results','search_stats.json')
S=json.load(open(st)) if os.path.exists(st) else None
txt={}
if S:
    txt['%%S_TIME%%']=f"the default search over all 91 rules takes {S['time_default_total_s']:.0f} seconds on a single core of a standard machine, at most {S['time_default_max_s']:.0f} seconds for a single rule, about two thirds of it spent on the five rules for which none is found."
    mp=S['min_p_per_rule']; import collections; c=collections.Counter(mp.values())
    txt['%%S_MINP%%']="the smallest period of a valid certificate is "+', '.join(f"$p={k}$ for {v} rules" for k,v in sorted(c.items()))+f", and the smallest $t_0$ for that period is at most {S['max_minimal_t0']}."
    b=S['enlarged_p24_t080']
    txt['%%S_BIG%%']="with $p^{\\max}=24$ and $t_0^{\\max}=80$ the search, run for "+f"{sum(v['seconds'] for v in b.values())/60:.0f} minutes in total, finds no certificate for any of them."
    n=S['n_valid_certificates_per_rule']
    txt['%%S_ALT%%']=f"Within the limits $p\\le12$, $t_0<40$, every certified rule has many valid certificates (between {min(n.values())} and {max(n.values())}, counting each choice of $(p,t_0,\\ell_B)$ once); for instance, the valid periods of Rule~37 are $p\\in\\{{{','.join(map(str,S['valid_p_values_per_rule']['37']))}\\}}$."
for k,v in txt.items(): s=s.replace(k,v)
for k,f in [('%%SPACETIME%%','spacetime'),('%%TAB_CERT2%%','tab_cert2'),('%%TAB_CERT%%','tab_cert'),('%%TAB_FALSE%%','tab_false'),('%%TAB_PROVED%%','tab_proved'),('%%TAB_OTHER%%','tab_other'),('%%BARS%%','bars')]:
    s=s.replace(k,open(os.path.join(T,f+'.tex')).read())
open(os.path.join(H,'main_full.tex'),'w').write(s)
for _ in range(2): subprocess.run(['pdflatex','-interaction=nonstopmode','main_full.tex'],cwd=H,stdout=subprocess.DEVNULL)
