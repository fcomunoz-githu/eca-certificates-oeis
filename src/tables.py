"""Build the LaTeX tables of the paper from results/*.json.  Usage: python src/tables.py"""
import os,re,json,collections
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); OUT=os.path.join(ROOT,'results','tables')
os.makedirs(OUT,exist_ok=True)
C=json.load(open(os.path.join(ROOT,'results','claims.json'))); K=json.load(open(os.path.join(ROOT,'results','certificates.json')))
def esc(s):
    s=s.replace('\\','\\textbackslash{}')
    for c in '&%$#_{}': s=s.replace(c,'\\'+c)
    return s.replace('^','\\^{}').replace('~','\\~{}')
def bt(s): return ''.join(esc(ch)+('\\allowbreak{}' if ch in '+-*,)=/ ' else '') for ch in s)
sig=lambda t: re.sub(r'\s*-\s*_[^_]+_,.*$','',t).strip()
TYP={'dec':'row, base 2','bin':'row, binary digits','on':'ON in row','off':'OFF in row','ton':'total ON','toff':'total OFF','mbin':'centre col., digits','mdec':'centre col., base 2','mid':'centre column'}
rows=[]
for k,v in K.items():
    if v is None: continue
    L=max(len(f['A'])+len(f['B'])+len(f['C'])+len(f['D'])+len(f['E']) for f in v['forms'])
    rows.append(f"{k} & {v['p']} & {v['t0']} & {v['lB']} & {v['lD']} & {L} & {sum(1 for c in C if c['rule']==int(k) and c['status']=='proved')} \\\\")
open(os.path.join(OUT,'tab_cert.tex'),'w').write('\n'.join(rows))
open(os.path.join(OUT,'tab_proved.tex'),'w').write('\n'.join(f"{c['id']} & {c['rule']} & {TYP[c['type']]} & {esc(c['author'] or '--')} & \\texttt{{{bt(sig(c['statement']))}}} \\\\" for c in C if c['status']=='proved'))
open(os.path.join(OUT,'tab_false.tex'),'w').write('\n'.join(f"{c['id']} & {c['rule']} & \\texttt{{{bt(sig(c['statement']))}}} & $n={[p['index'] for p in c['parts'] if p['result']=='fail'][0]}$ & {esc(c['correction'])} \\\\" for c in C if c['status']=='false-as-stated'))
open(os.path.join(OUT,'tab_other.tex'),'w').write('\n'.join(f"{c['id']} & {c['rule']} & \\texttt{{{bt(sig(c['statement'])[:160])}}} & {'no certificate' if c['status']=='no-certificate' else 'format not covered'} \\\\" for c in C if c['status'] in('no-certificate','not-covered')))
# per-type summary for the bar chart
order=['dec','bin','on','off','ton','toff','mid','mbin','mdec']
cnt={t:collections.Counter(c['status'] for c in C if c['type']==t) for t in order}
with open(os.path.join(OUT,'bars.tex'),'w') as f:
    for st in ['proved','false-as-stated','not-covered','no-certificate']:
        f.write('\\addplot coordinates {'+' '.join(f'({t},{cnt[t][st]})' for t in order)+'};\n')
print('tables written to',OUT)
# compact two-block certificate table (Table 3)
half=(len(rows)+1)//2; L1=rows[:half]; L2=rows[half:]+['& & & & & & \\\\']*(half-len(rows[half:]))
hdr='rule & $p$ & $t_0$ & $\\ell_B$ & $\\ell_D$ & $L$ & pr.'
body='\n'.join(a.replace('\\\\','').rstrip()+' & '+b for a,b in zip(L1,L2))
open(os.path.join(OUT,'tab_cert2.tex'),'w').write('\\begin{tabular}{rrrrrrr|rrrrrrr}\\toprule\n'+hdr+' & '+hdr+'\\\\\\midrule\n'+body+'\n\\bottomrule\\end{tabular}')
