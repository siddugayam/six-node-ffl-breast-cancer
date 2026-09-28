#!/bin/bash
# 99c_dynamics_driver2.sh -- runs the analyses that were written but never executed,
# plus the additions from this pass, on the existing 1,024-set 3-node ensemble.
REV=/path/to/revision
cd "$REV" || exit 1
L="$REV/logs/v3"
run () {
  echo "--- $1 start $(date -Is)"
  python3 -u "$REV/scripts/v3/$1" > "$L/$2" 2>&1
  echo "--- $1 exit=$? $(date -Is)"
}
echo "driver2 start $(date -Is)"
run 02b_summary_table.py                        02b_summary_table.log
run 08_mir_vs_txn_paired.py                     08_mir_vs_txn_paired.log
run 02e_paired_tests_and_textbook_validation.py 02e_paired_textbook.log
run 02g_fcd_recomputed.py                       02g_fcd_recomputed.log
run 02f_pulse_longwindow.py                     02f_pulse_longwindow.log
run 06b_composite_showcase.py                   06b_composite_showcase.log
run 03e_higher_order_convergence_audit.py       03e_convergence_audit.log
echo "driver2 done $(date -Is)"
