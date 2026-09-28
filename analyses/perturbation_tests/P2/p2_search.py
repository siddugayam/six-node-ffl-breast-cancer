#!/usr/bin/env python3
"""P2 search: GEO for perturbations of co-transcribed miRNA clusters (the analysis plan's five primary and four secondary
clusters): cluster knockout/deletion/knockdown, cluster over-expression, or co-transfection of >= 2 members.
Homo sapiens (primary) and Mus musculus (secondary), expression profiling; then the first-pass screen."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import _screen
ACT = '(knockout OR deletion OR deleted OR knockdown OR overexpression OR "over-expression" OR transgenic OR mimic OR mimics OR inhibitor OR inhibitors OR sponge OR antagomir OR "co-transfection" OR transfection OR cluster)'
CL = {
 'miR17_92': '("miR-17-92" OR "miR-17~92" OR "miR-17/92" OR "mir-17-92" OR "miR17-92" OR "MIR17HG" OR "Mir17hg" OR "oncomiR-1")',
 'miR106a_363': '("miR-106a-363" OR "miR-106a~363" OR "miR-106a/363" OR "mir-106a-363")',
 'DLK1_DIO3': '("DLK1-DIO3" OR "Dlk1-Dio3" OR "miR-379/410" OR "miR-379-410" OR "miR-379~410" OR "Mirg" OR "14q32 microRNA" OR "14q32 miRNA" OR "miR-379/656")',
 'miR29ab1': '("miR-29a/b-1" OR "miR-29a/b1" OR "miR-29ab1" OR "miR-29a~b-1" OR "Mir29ab1" OR "miR-29 cluster" OR "miR-29a/b")',
 'miR183_96_182': '("miR-183-96-182" OR "miR-183/96/182" OR "miR-183~96~182" OR "miR-183 cluster" OR "Mir183" OR "miR-183/182/96")',
 'miR106b_25': '("miR-106b-25" OR "miR-106b~25" OR "miR-106b/25" OR "miR-106b/93/25" OR "Mir106b-25")',
 'miR200b_200a_429': '("miR-200b-200a-429" OR "miR-200b/200a/429" OR "miR-200b~200a~429" OR "miR-200 family" OR "miR-200 cluster")',
 'miR200c_141': '("miR-200c-141" OR "miR-200c/141" OR "miR-200c~141")',
 'miR23a_27a_24': '("miR-23a-27a-24" OR "miR-23a/27a/24" OR "miR-23a~27a~24" OR "miR-23a cluster" OR "miR-23~27~24")',
}
hits = {}
for org, tag in (('Homo sapiens', 'hs'), ('Mus musculus', 'mm')):
    h = _screen.search(f'P2_{tag}', [(f'P2_{tag}_{k}', f'{v} AND {ACT}') for k, v in CL.items()], organism=org)
    for g, r in h.items(): r['organism_query'] = org; hits[g] = r
MIR = r'(mi-?R-?(17|18a|19a|19b|20a|92|106a|18b|20b|363|29|183|96|182|106b|93|25|200|141|429|23a|27a|24|379|410|127|134|136|154|299|323|329|337|369|370|376|377|381|382|409|411|412|431|432|433|485|487|493|494|495|496|539|541|543|544|654|655|656)|mir-?1[07]|Mirg|cluster|MIR17HG|Mir17hg)'
PERT = r'(knock|\bKO\b|-/-|del|null|flox|over-?express|\bOE\b|transgen|mimic|inhibit|sponge|antago|transfect|\bTg\b)'
rows = _screen.screen('P2', hits, r'(' + MIR + r'.{0,40}' + PERT + r'|' + PERT + r'.{0,40}' + MIR + r')',
                      r'(control|ctrl|scrambl|non[- ]?target|\bNT\b|\bNC\b|negative|mock|\bWT\b|wild[- ]?type|littermate|empty|vector|GFP|untreated|\+/\+)')
print(len(hits), 'GSE hits; >=2 auto-perturbed and >=2 auto-control:', sum(1 for r in rows if r['auto_perturbed'] >= 2 and r['auto_control'] >= 2))
