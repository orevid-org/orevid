# Hit rate by evidence (record)

Generated with `py -3 -m evaluation.evidence` (14-10-2026, after P40: hypergeometric population = the 19,735 genes with at least one disease, p < 1e-5; `docs/HYPERGEOMETRIC_N.md`). Hits among the first 25 candidates of each disease. Disease-level confidence intervals, threshold sensitivity and the permutation null: `docs/EVIDENCE_CI.md` (P38: the shared genes are not evidence of a shared mechanism). Decisions: `docs/DECISIONS.md` P27-P31, P38, P40.

```
complex_disease bench, seeds 1-5 (where the levels were fixed)
  all       A: 45.3% ( 2100)   B: 22.1% ( 1561)   C: 10.0% ( 1486)   -:  1.5% ( 3228)   neighbour: 49.9% ( 1691)   table: 23.0% ( 1970)   hidden:  4.2% ( 4714)
  hard      A: 29.2% (  514)   B: 18.9% (  623)   C:  4.8% (  421)   -:  2.4% ( 1667)   neighbour: 13.3% (  309)   table: 27.4% (  828)   hidden:  2.9% ( 2088)

Confirmation: 9 hidden areas, seed 1
  all       A: 37.4% ( 2098)   B: 24.6% ( 1151)   C:  4.3% (  558)   -:  1.4% ( 4818)   neighbour: 59.4% (  641)   table: 26.3% ( 2608)   hidden:  1.7% ( 5376)
  hard      A: 19.5% ( 1034)   B: 26.3% (  911)   C:  1.1% (  364)   -:  1.6% ( 3516)   neighbour: 24.7% (  219)   table: 22.5% ( 1726)   hidden:  1.6% ( 3880)

External exam (645 diseases), hits among the top 25 not known:
  neighbour                12.6% (3025)
  whole table               4.3% (3936)
  table, specific genes     4.8% (2875)
  table, generic genes      2.8% (1061)
  hidden                    1.5% (9164)
  all                       4.2% (16125)
```
