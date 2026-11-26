"""TikZ space-time diagram with the certificate regions highlighted (Figure 3). Usage: python src/figure_spacetime.py RULE TMAX"""
import os,sys
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from ca import rows
from castruct import certify
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
rule=int(sys.argv[1]) if len(sys.argv)>1 else 25; TM=int(sys.argv[2]) if len(sys.argv)>2 else 18
p,t0,F=certify(rule); R=rows(rule,TM); T=t0+2*p
o=[r'\begin{tikzpicture}[x=2.6mm,y=-2.6mm]']
for t in range(TM+1):
    W=R[t][0]; x0=-(t+1)
    if t>=T:
        r=(t-t0)%p; m=(t-t0-r)//p; A,B,C,D,E=F[r]
        xb=x0+len(A); xd=x0+len(A)+m*len(B)+len(C)
        o.append(rf'\fill[blue!25] ({xb},{t}) rectangle ({xb+m*len(B)},{t+1});')
        o.append(rf'\fill[red!25] ({xd},{t}) rectangle ({xd+m*len(D)},{t+1});')
        o.append(rf'\fill[green!20] ({xb+m*len(B)},{t}) rectangle ({xd},{t+1});')
    for i,c in enumerate(W):
        if c: o.append(rf'\fill[black!80] ({x0+i+0.18},{t+0.18}) rectangle ({x0+i+0.82},{t+0.82});')
    o.append(rf'\node[anchor=east,font=\tiny] at ({-TM-2.5},{t+0.5}) {{{t}}};')
o.append(rf'\draw[dashed] ({-TM-2},{T}) -- ({TM+2},{T}) node[right,font=\scriptsize]{{$T={T}$}};')
o.append(r'\end{tikzpicture}')
open(os.path.join(ROOT,'results','tables','spacetime.tex'),'w').write('\n'.join(o))
print(p,t0,[[''.join(map(str,w)) for w in f] for f in F])
