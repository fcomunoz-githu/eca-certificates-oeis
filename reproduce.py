"""One command to reproduce everything:  python reproduce.py   (add --full to also run the search statistics, ~20 min)"""
import subprocess,sys,os
H=os.path.dirname(os.path.abspath(__file__))
def run(*a): print('>>',' '.join(a)); subprocess.run([sys.executable,*a],cwd=H,check=True)
run('src/run_all.py'); run('src/tables.py'); run('src/figure_spacetime.py','25','18')
if '--full' in sys.argv: run('src/search_stats.py')
run('-m','unittest','discover','-s','tests','-v')
try: run('paper/build.py'); print('PDF: paper/main_full.pdf')
except Exception as e: print('LaTeX build skipped:',e)
