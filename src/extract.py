import os,re,sys,json,csv
root=sys.argv[1]; out=sys.argv[2]  # usage: python extract.py <oeisdata/seq> <out.csv>
rows=[]
auth=re.compile(r'_([^_]{2,60})_,\s*([A-Z][a-z]{2} \d{1,2} \d{4})')
for d,_,fs in os.walk(root):
  for f in fs:
    if not f.endswith('.seq'): continue
    try: L=open(os.path.join(d,f),encoding='utf-8',errors='replace').read().splitlines()
    except: continue
    name='';kw='';data=''
    for l in L:
      if l.startswith('%N'): name=l[11:]
      if l.startswith('%K'): kw=l[11:]
      if l[:2] in('%S','%T','%U'): data+=l[11:]
    if 'dead' in kw: continue
    full='\n'.join(L).lower()
    proved_hint=bool(re.search(r'\b(proved|proof|proven|confirmed|established|theorem)\b',full))
    for i,l in enumerate(L):
      if re.search(r'conjectur',l,re.I) and l[:2] in ('%C','%F','%E'):
        txt=l[11:]; m=auth.findall(txt); j=i+1
        while not m and j<len(L) and j<i+15 and L[j][:2]==l[:2]:
          txt+=' '+L[j][11:]; m=auth.findall(L[j][11:]); j+=1
        if not m: txt=l[11:]
        rows.append(dict(id=f[:-4],field=l[1],name=name,keywords=kw,text=txt,author=m[-1][0] if m else '',date=m[-1][1] if m else '',proved_hint=proved_hint,data=data[:300]))
rows.sort(key=lambda r:r['id'])
with open(out,'w',newline='',encoding='utf-8') as fo:
  w=csv.DictWriter(fo,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print(len(rows),len({r['id'] for r in rows}))
