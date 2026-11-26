# Uniform proofs of OEIS conjectures on elementary cellular automata

Supplementary material for

> F. J. Muñoz, *Finite certificates for eventually regular elementary cellular automata, with proofs of
> 575 conjectures from the OEIS* (preprint, 2026).

`python reproduce.py` reproduces every number, table and figure of the paper from the 260 OEIS entries
stored in `data/oeis_seq.zip`, runs the tests and rebuilds the PDF (about 3 minutes).

## Contents

| path | content |
|---|---|
| `paper/main.tex`, `paper/main_full.tex`, `paper/main_full.pdf` | the paper (source with table placeholders, assembled source, PDF) |
| `src/ca.py` | ECA simulation with infinite background; certificate search; checks (C1)–(C4) of Definition 2 (`prove_structure`) |
| `src/castruct.py` | certificate search driver (`certify`) |
| `src/seqs.py` | the nine sequence types and their generation |
| `src/verify.py` | annihilating polynomials of Proposition 5 (`annihilator`) and the finite checks of Lemma 7 (`check_rec`, `check_gf`, `check_closed`) |
| `src/run_all.py` | runs everything and writes `results/` |
| `src/tables.py` | builds the LaTeX tables in `results/tables/` |
| `src/figure_spacetime.py` | the space-time figure with the certificate regions (Rule 25) |
| `src/search_stats.py` | statistics of Section 4: timing, alternative certificates, enlarged limits (`results/search_stats.json`, about 45 min) |
| `tests/test_all.py` | test suite (`python -m unittest discover -s tests`) |
| `reproduce.py` | runs everything in order |
| `src/extract.py` | extracts all "conjecture" lines from a full OEIS data export (used to build `data/conjectures_all_oeis.csv`) |
| `data/oeis_seq.zip` | the 260 OEIS entries used, in OEIS internal format, unchanged (unpacked automatically into `data/oeis_seq/` on first run) |
| `data/eca_entries.txt` | their A-numbers |
| `data/conjectures_all_oeis.csv` | all 27,776 "conjecture" lines of the OEIS (21,086 entries), for future work |
| `data/OEIS_SNAPSHOT.txt` | exact version of the OEIS data |
| `results/certificates.json` | the certificates of the 86 rules (p, t0, ℓB, ℓD and the words A_r, B_r, C_r, D_r, E_r) |
| `results/claims.json` | the 640 statements, each with A-number, rule, type, author, date, outcome, parts checked, annihilator, threshold N, and correction when false |
| `results/summary.json` | the counts reported in the paper |
| `results/sanity_annihilators.txt` | independent numerical checks (not needed for the proofs) |

## How to reproduce

Requirements: Python ≥ 3.10 and SymPy (tested with Python 3.13 and SymPy 1.14.0); a TeX Live
installation with `pgfplots` to rebuild the PDF.

```
pip install -r requirements.txt
python reproduce.py          # results, tables, figure, tests, PDF (about 3 minutes)
python reproduce.py --full   # also the search statistics of Section 4 (about 45 more minutes)
```

Step by step: `python src/run_all.py` (results/*.json), `python src/tables.py` (results/tables),
`python src/figure_spacetime.py 25 18`, `python -m unittest discover -s tests`, `python paper/build.py`.

To rebuild `data/conjectures_all_oeis.csv` from scratch:

```
git clone --filter=blob:none --sparse https://github.com/oeis/oeisdata.git
cd oeisdata && git checkout be2f360122f661bddebad22b15440935d21ac598 && git sparse-checkout set seq && cd ..
python src/extract.py oeisdata/seq data/conjectures_all_oeis.csv
```

## Status values in `results/claims.json`

- `proved`: every part of the statement verified by Lemma 7 with a valid certificate.
- `false-as-stated`: a counterexample was found; `correction` gives the exact valid range when it exists.
- `not-covered`: format outside Lemma 7 (floor of exponentials, e.g.f., identities with other A-numbers, prose).
- `no-certificate`: rules 129, 137, 167, 169, 182.

## Licences

Code: MIT licence (see `LICENSE`). The files in `data/oeis_seq.zip` and the text of the statements are
from The On-Line Encyclopedia of Integer Sequences (https://oeis.org), licensed under CC BY-SA 4.0.
