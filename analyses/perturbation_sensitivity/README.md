# Sensitivity analyses of the perturbation tests (Supplementary Notes S1 and S2)

Three sensitivity analyses of the primary results. They read the result tables of `analyses/perturbation_tests/`.

- `C1/c1_p1b_series_level.py`: the miR-29 fibroblast–carcinoma difference with the cell lines of each series averaged.
- `C2/c2_p3a_without_gse115646.py`: the strong-tier shift without the GSE115646 series.
- `C3/c3_p2_mouse_leave_one_series_out.py`: the mouse results, leaving out one series at a time.

`analyses/perturbation_tests/_rma.R` pools the estimates (random-effects model, metafor).
