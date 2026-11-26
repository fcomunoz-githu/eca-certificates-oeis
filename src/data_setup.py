"""Unpack data/oeis_seq.zip (the 260 OEIS entries) into data/oeis_seq/ if needed."""
import os,zipfile
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
def ensure_data():
    d=os.path.join(ROOT,'data','oeis_seq')
    if not os.path.isdir(d) or len(os.listdir(d))<260:
        with zipfile.ZipFile(os.path.join(ROOT,'data','oeis_seq.zip')) as z: z.extractall(os.path.join(ROOT,'data'))
    return d
if __name__=='__main__': print(ensure_data())
