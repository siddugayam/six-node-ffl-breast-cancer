#!/usr/bin/env python3
"""Table S1 with the analysed-network flag: every row of results/v5/tables/TableS1_all_interactions.csv
(6,859 deposited edges) plus 'in_analysed_network' (FALSE for the 30 legacy GeneMANIA miRNA_miRNA edges
of the original circuits, TRUE for the 6,829 analysed edges). The project file is only read."""
import csv
SRC = '/path/to/revision/results/v5/tables/TableS1_all_interactions.csv'
R = list(csv.DictReader(open(SRC))); cols = list(R[0].keys())
n = 0
with open('TableS1_all_interactions_flagged.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=cols + ['in_analysed_network']); w.writeheader()
    for r in R:
        leg = r['edge_type'] == 'miRNA_miRNA'; n += leg
        w.writerow({**r, 'in_analysed_network': 'FALSE' if leg else 'TRUE'})
assert len(R) == 6859 and n == 30, (len(R), n)
print(f'{len(R)} rows, {n} flagged FALSE, {len(R) - n} analysed')
