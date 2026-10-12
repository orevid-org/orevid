# Hit rate by evidence (record)

Generated with `py -3 -m evaluation.evidence` (12-10-2026). Hits among the first 25 candidates of each disease. Decisions: `docs/DECISIONS.md` P27-P31.

```

Benchmark complex_disease, seeds 1-5 (where the levels were fixed)
  all       A: 44.9% ( 2115)   B: 22.1% ( 1568)   C: 10.1% ( 1476)   -:  1.6% ( 3216)   neighbour: 49.9% ( 1691)   table: 22.7% ( 1992)   hidden:  4.2% ( 4692)
  hard      A: 28.7% (  523)   B: 18.8% (  627)   C:  4.8% (  416)   -:  2.4% ( 1659)   neighbour: 13.3% (  309)   table: 27.0% (  841)   hidden:  2.9% ( 2075)

Confirmation: 9 held-out areas, seed 1
  all       A: 36.2% ( 2167)   B: 23.6% ( 1199)   C:  4.3% (  552)   -:  1.4% ( 4707)   neighbour: 59.4% (  641)   table: 25.2% ( 2725)   hidden:  1.7% ( 5259)
  hard      A: 18.5% ( 1094)   B: 25.3% (  948)   C:  1.1% (  358)   -:  1.7% ( 3425)   neighbour: 24.7% (  219)   table: 21.3% ( 1823)   hidden:  1.6% ( 3783)

External exam (645 diseases), hits among the first 25 not already known:
  neighbour                12.6% (3025)
  full table                4.2% (4119)
  table, specific genes     4.8% (3007)
  table, generic genes      2.7% (1112)
  hidden                    1.4% (8981)
  all                       4.2% (16125)
```
