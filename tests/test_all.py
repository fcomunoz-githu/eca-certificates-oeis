"""Quick test suite (about 1 minute).  Usage:  python -m unittest discover tests"""
import os,sys,json,unittest,re
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0,os.path.join(ROOT,'src'))
from data_setup import ensure_data
ensure_data()
from ca import rows,step,prove_structure,find_structure
from castruct import certify
from seqs import gen,classify
from verify import annihilator,verify_statement,check_rec,S
import sympy as sp

class TestSimulation(unittest.TestCase):
    def test_background_odd_rule(self):
        # Rule 67 has f(0,0,0)=1: background alternates, rows 0..2 as in A266839 (1, 4, 10)
        self.assertEqual(gen(67,'dec',2),[1,4,10])
    def test_oeis_data_sample(self):
        for i in ['A262779','A266592','A266839','A267684','A265382']:
            L=open(os.path.join(ROOT,'data','oeis_seq',i+'.seq')).read().splitlines()
            d={}; [d.setdefault(l[:2],[]).append(l[11:]) for l in L]
            data=[int(x) for x in ''.join(d.get('%S',[])+d.get('%T',[])+d.get('%U',[])).replace(' ','').strip(',').split(',')]
            off=int(d['%O'][0].split(',')[0]); rule,typ=classify(d['%N'][0])
            self.assertEqual(gen(rule,typ,off+len(data))[off:off+len(data)],data,i)

class TestCertificates(unittest.TestCase):
    def test_stored_certificates_are_valid(self):
        K=json.load(open(os.path.join(ROOT,'results','certificates.json')))
        for ru,v in K.items():
            if v is None: continue
            F=[tuple([int(c) for c in f[k]] for k in 'ABCDE') for f in v['forms']]
            R=rows(int(ru),v['t0']+6*v['p']+8)
            self.assertTrue(prove_structure(int(ru),R,v['p'],v['t0'],F)[0],ru)
    def test_corrupted_certificate_rejected(self):
        p,t0,F=certify(37); R=rows(37,t0+6*p+8)
        A,B,C,D,E=F[0]; bad=[(A,B,C[:-1]+[1-C[-1]],D,E)]+list(F[1:])
        self.assertFalse(prove_structure(37,R,p,t0,bad)[0])

class TestVerification(unittest.TestCase):
    def setUp(self):
        p,t0,F=certify(175); self.P,self.N=annihilator('bin',p,F,t0); self.a=gen(175,'bin',self.N+100)
    def test_true_statements(self):
        self.assertEqual(verify_statement('a(n) = 101*a(n-1)-100*a(n-2) for n>2.',self.a,0,self.P,self.N)[1][0],'ok')
        self.assertEqual(verify_statement('G.f.: (1+9*x+x^2) / ((1-x)*(1-100*x)).',self.a,0,self.P,self.N)[1][0],'ok')
    def test_false_statements(self):
        self.assertEqual(verify_statement('a(n) = 101*a(n-1)-101*a(n-2) for n>2.',self.a,0,self.P,self.N)[1][0],'fail')
        self.assertEqual(verify_statement('G.f.: (1+9*x+x^2) / ((1-x)*(1-101*x)).',self.a,0,self.P,self.N)[1][0],'fail')
    def test_annihilator_numerically(self):
        cs=self.P.all_coeffs()[::-1]
        for t in range(self.N,self.N+50): self.assertEqual(sum(c*self.a[t+i] for i,c in enumerate(cs)),0)

class TestResults(unittest.TestCase):
    def test_summary_numbers(self):
        s=json.load(open(os.path.join(ROOT,'results','summary.json')))
        self.assertEqual((s['statements'],s['status']['proved'],s['status']['false-as-stated'],s['certified_rules']),(640,575,15,86))

if __name__=='__main__': unittest.main()
