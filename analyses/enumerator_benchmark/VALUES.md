# VALUES: benchmark of the topological-definition enumerator on planted modules

Lines count the header as line 1. Tables are written by `benchmark.py`.
`benchmark_summary.csv` gives, per measure (column measure), the number of replicates and the median, min, max
and mean over replicates. Below, each value is the median (min–max) over replicates, with its line in that file.
Per-replicate values are in `benchmark_per_replicate.csv` (one row per replicate, columns as the measures).

## 1. the topological definition topologies (`planted_topologies.csv`, column D1_D4_topologies_up_to_isomorphism)

- line 2: n = 3: 1 topologies
- line 3: n = 4: 7 topologies
- line 4: n = 5: 79 topologies
- line 5: n = 6: 1,759 topologies

## 2. Graphs

- background edges (the contracted analysed network, rewired): 5,606 (5,606–5,606) [line 3]
- noise arcs: 561 (561–561) [line 82]
- arcs of the final graph: 6,408 (6,387–6,425) [line 2]
- replicates: 200

## 3. The topological-definition enumerator

| n | planted | intact | reported | recall of intact planted | destroyed planted reported | precision (sample of 1,000) |
|---|---|---|---|---|---|---|
| 3 | 10 (10–10) [line 19] | 10 (10–10) [line 18] | 5,950.5 (5,655–6,149) [line 22] | 1 (1–1) [line 21] | 0 (0–0) [line 14] | 1 (1–1) [line 20] |
| 4 | 10 (10–10) [line 40] | 10 (10–10) [line 39] | 93,719 (87,834–99,188) [line 43] | 1 (1–1) [line 42] | 0 (0–0) [line 35] | 1 (1–1) [line 41] |
| 5 | 10 (10–10) [line 58] | 10 (10–10) [line 57] | 1,307,937 (1,200,566–1,434,029) [line 61] | 1 (1–1) [line 60] | 0 (0–0) [line 56] | 1 (1–1) [line 59] |
| 6 | 10 (10–10) [line 76] | 10 (10–10) [line 75] | 16,980,879.5 (14,967,550–19,485,262) [line 79] | 1 (1–1) [line 78] | 0 (0–0) [line 74] | 1 (1–1) [line 77] |

Exhaustive cross-check against the independent Python enumeration:

- n = 3 (200 replicates): modules 5,950.5 (5,655–6,149) [line 17]; recall 1 (1–1) [line 16]; precision 1 (1–1) [line 15]
- n = 4 (5 replicates): modules 95,291 (93,095–96,488) [line 38]; recall 1 (1–1) [line 37]; precision 1 (1–1) [line 36]

## 4. Relaxations

ratio = modules reported under the relaxed test / modules reported under the topological definition; contains = the relaxed set contains
every topological-definition module (1 = true in every replicate); share planted = planted modules reported / all reported.

| n | test | reported | ratio | contains the topological definition set | recall of intact planted | share planted |
|---|---|---|---|---|---|---|
| 3 | the topological definition | 5,950.5 (5,655–6,149) [line 22] | 1 | – | 1 (1–1) [line 21] | 0.001681 (0.001626–0.001768) [line 24] |
| 3 | D3 relaxed | 5,950.5 (5,655–6,149) [line 7] | 1 (1–1) [line 5] | 1 (1–1) [line 4] | 1 (1–1) [line 6] | 0.001681 (0.001626–0.001768) [line 8] |
| 3 | D4 relaxed | 69,627.5 (65,919–73,888) [line 12] | 11.72 (10.97–12.56) [line 10] | 1 (1–1) [line 9] | 1 (1–1) [line 11] | 0.0001436 (0.0001353–0.0001517) [line 13] |
| 4 | the topological definition | 93,719 (87,834–99,188) [line 43] | 1 | – | 1 (1–1) [line 42] | 0.0001067 (0.0001008–0.0001139) [line 45] |
| 4 | D3 relaxed | 93,719 (87,834–99,188) [line 28] | 1 (1–1) [line 26] | 1 (1–1) [line 25] | 1 (1–1) [line 27] | 0.0001067 (0.0001008–0.0001139) [line 29] |
| 4 | D4 relaxed | 608,173.5 (553,538–678,680) [line 33] | 6.477 (5.989–7.085) [line 31] | 1 (1–1) [line 30] | 1 (1–1) [line 32] | 1.644e-05 (1.473e-05–1.807e-05) [line 34] |
| 5 | the topological definition | 1,307,937 (1,200,566–1,434,029) [line 61] | 1 | – | 1 (1–1) [line 60] | 7.646e-06 (6.973e-06–8.329e-06) [line 63] |
| 5 | D3 relaxed | 1,307,937 (1,200,566–1,434,029) [line 49] | 1 (1–1) [line 47] | 1 (1–1) [line 46] | 1 (1–1) [line 48] | 7.646e-06 (6.973e-06–8.329e-06) [line 50] |
| 5 | D4 relaxed | 5,062,192.5 (4,548,622–5,877,839) [line 54] | 3.871 (3.526–4.314) [line 52] | 1 (1–1) [line 51] | 1 (1–1) [line 53] | 1.975e-06 (1.701e-06–2.198e-06) [line 55] |
| 6 | the topological definition | 16,980,879.5 (14,967,550–19,485,262) [line 79] | 1 | – | 1 (1–1) [line 78] | 5.889e-07 (5.132e-07–6.681e-07) [line 81] |
| 6 | D3 relaxed | 16,980,879.5 (14,967,550–19,485,262) [line 67] | 1 (1–1) [line 65] | 1 (1–1) [line 64] | 1 (1–1) [line 66] | 5.889e-07 (5.132e-07–6.681e-07) [line 68] |
| 6 | D4 relaxed | 45,257,474.5 (39,502,662–53,845,966) [line 72] | 2.66 (2.401–2.983) [line 70] | 1 (1–1) [line 69] | 1 (1–1) [line 71] | 2.21e-07 (1.857e-07–2.531e-07) [line 73] |

D3 relaxed, ratio 1 in every replicate at every n: True. When true, the D3-relaxed test reports the same modules as the topological definition. topological-definition modules are
acyclic, and in an acyclic subgraph with one source and one sink every vertex lies on a source–sink path, so dropping
that clause of D3 removes nothing.

