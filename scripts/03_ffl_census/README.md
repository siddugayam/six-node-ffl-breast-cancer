# Feed-forward loops: three-node cores and the n-node census (Methods 2.2, Fig. S2, Supplementary Note S2)

- Definitions and three-node cores: `01a_counting_convention.py`, `01b_coherence_corrected.py` (Mangan–Alon coherence
  sub-types, Table S4); `ffl_def.py`, `ffl_graph.py`, `ffl_cores.py` → `02_ffl_census.py` and the follow-up analyses
  `02b` … `02y` (edge-class coverage, candidate space, the deposited exemplar circuits, cores per source network).
- n-node census: `03_ffl_census.py`, `03_ffl_enumerate.py`, the C enumerators (`ffl_enum*.c`, `v2_census*.c`,
  `R5_census.c`) with their run scripts (`run_*.sh`), the `R*` and `v2_*` graph and class-diversity scripts, and
  `03_ffl_census_nolegacymirna.py` (the census graph without the 30 exemplar edges).
- `11_ffl_module_membership.R`: node sets of the 3- to 6-node modules. `60_census_*`: STRING sensitivity.
The census of the paper was re-run in `analyses/census_and_motif_nulls/`.
