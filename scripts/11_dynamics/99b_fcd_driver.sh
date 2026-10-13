#!/bin/bash
# waits for the enlarged 3-node ensemble, then recomputes fold-change detection
REV=/path/to/six-node-ffl-breast-cancer
L="$REV/logs/systems_and_dynamics"
while ! grep -q "^DONE in" "$L/02c_characterise_3node_N16384.log" 2>/dev/null; do sleep 30; done
sleep 60   # let the driver's fast steps go first
python3 -u "$REV/scripts/11_dynamics/02g_fcd_peak_amplitude.py" > "$L/02g_fcd_peak_amplitude.log" 2>&1
echo "02g exit=$?"
