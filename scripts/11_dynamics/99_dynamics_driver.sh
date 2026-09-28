#!/bin/bash
# 99_dynamics_driver.sh
# Runs the remaining dynamical-analysis steps in order once the enlarged
# (N = 16,384) 3-node ensemble has been written by 02c_characterise_3node_N16384.py.
# Each step logs to logs/v3/ and its exit status is appended to the driver log.
REV=/path/to/revision
cd "$REV" || exit 1
L="$REV/logs/v3"

echo "driver start $(date -Is)"

# --- wait for the enlarged ensemble -----------------------------------------
while ! grep -q "^DONE in" "$L/02c_characterise_3node_N16384.log" 2>/dev/null; do
  if grep -q "^EXIT=" "$L/02c_characterise_3node_N16384.log" 2>/dev/null; then
    if ! grep -q "^DONE in" "$L/02c_characterise_3node_N16384.log"; then
      echo "02c FAILED - aborting driver"; exit 1
    fi
  fi
  sleep 30
done
echo "02c complete $(date -Is)"

run () {   # run <script> <logname>
  echo "--- $1 start $(date -Is)"
  python3 -u "$REV/scripts/11_dynamics/$1" > "$L/$2" 2>&1
  echo "--- $1 exit=$? $(date -Is)"
}

run 02b_summary_table.py                        02b_summary_table.log
run 08_mir_vs_txn_paired.py                     08_mir_vs_txn_paired.log
run 02e_paired_tests_and_textbook_validation.py 02e_paired_textbook.log
run 06b_composite_showcase.py                   06b_composite_showcase.log
run 02f_pulse_longwindow.py                     02f_pulse_longwindow.log
run 03d_higher_order_compI1.py                  03d_higher_order_compI1.log
run 09b_composite_oscillation.py                09b_oscillation.log
run 03e_higher_order_convergence_audit.py       03e_convergence_audit.log
run 05_figures.py                               05_figures.log
run 07_consolidate.py                           07_consolidate.log

echo "driver done $(date -Is)"
