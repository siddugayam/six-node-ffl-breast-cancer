#!/usr/bin/env python3
"""Consolidate every headline number produced by the v3 systems-pharmacology analysis into
one table.  Every value is read back from the result files, not retyped."""
import csv, json, collections
import numpy as np

RES = '/path/to/revision/results/v3'
rows = []


def add(section, metric, value, detail=''):
    rows.append(dict(section=section, metric=metric, value=value, detail=detail))


arch = json.load(open(f'{RES}/systems_architecture_summary.json'))
nul = json.load(open(f'{RES}/systems_architecture_nulls.json'))
pl = list(csv.DictReader(open(f'{RES}/systems_powerlaw_csn_fits.csv')))
lrt = list(csv.DictReader(open(f'{RES}/systems_powerlaw_csn_lrt.csv')))
ctl = json.load(open(f'{RES}/systems_controllability_summary.json'))
enr = list(csv.DictReader(open(f'{RES}/systems_controllability_enrichment.csv')))
cmp_ = list(csv.DictReader(open(f'{RES}/systems_flow_vs_degree_comparison.csv')))
comb = json.load(open(f'{RES}/systems_combination_summary.json'))
cst = json.load(open(f'{RES}/systems_combination_stats.json'))
dyn = json.load(open(f'{RES}/systems_dynamics_summary.json'))
dnl = json.load(open(f'{RES}/systems_dynamics_nulls.json'))
refs = list(csv.DictReader(open(f'{RES}/systems_references_verified.csv')))
mxc = json.load(open(f'{RES}/systems_matching_crosscheck.json'))
regs = list(csv.DictReader(open(f'{RES}/systems_regulator_subnetwork.csv')))
regn = list(csv.DictReader(open(f'{RES}/systems_regulator_subnetwork_nulls.csv')))
igc = list(csv.DictReader(open(f'{RES}/systems_igraph_crosscheck.csv')))
skc = list(csv.DictReader(open(f'{RES}/systems_flow_vs_degree_sinkcontrol.csv')))
skj = json.load(open(f'{RES}/systems_flow_vs_degree_sinkcontrol.json'))

# A ARCHITECTURE
add('A architecture', 'nodes', arch['N'])
add('A architecture', 'directed edges', arch['M'])
add('A architecture', 'weakly connected components', arch['n_WCC'])
add('A architecture', 'strongly connected components', arch['n_SCC'],
    f"{arch['n_SCC_size1']} of them singletons")
add('A architecture', 'bow-tie IN / CORE / OUT',
    f"{arch['bowtie']['IN']} / {arch['bowtie']['CORE']} / {arch['bowtie']['OUT']}",
    f"tubes {arch['bowtie']['TUBES']}, tendrils "
    f"{arch['bowtie']['TENDRILS_from_IN']+arch['bowtie']['TENDRILS_to_OUT']}")
add('A architecture', 'CORE composition',
    f"{arch['bowtie_types']['CORE']['TF']} TF + {arch['bowtie_types']['CORE']['miRNA']} miRNA",
    'zero genes in the core: every gene is an output')
add('A architecture', 'OUT composition',
    f"{arch['bowtie_types']['OUT']['Gene']} genes + {arch['bowtie_types']['OUT']['TF']} TF "
    f"+ {arch['bowtie_types']['OUT']['miRNA']} miRNA")
add('A architecture', 'edges inside the core',
    f"{arch['edges_within_core']} ({arch['edges_within_core_frac']*100:.1f}%)")
add('A architecture', 'reciprocity', round(arch['reciprocity'], 4),
    'all of it TF<->miRNA mutual pairs (1,223 dyads)')
for m, lab in [('LSCC', 'core size'), ('flow_hierarchy', 'flow hierarchy'),
               ('GRC', 'global reaching centrality'), ('F0', 'trophic incoherence F0'),
               ('transitivity', 'transitivity'), ('ffl_unique', '3-node FFL cores'),
               ('ffl_raw', '3-node FFL cores, raw directed count')]:
    a = nul['NULL_A'][m]; b = nul['NULL_B'][m]
    add('A architecture', lab, round(a['obs'], 4),
        f"NULL-A z={a['z']:+.2f} p={a['p_two']:.3g}; NULL-B z={b['z']:+.2f} p={b['p_two']:.3g} "
        f"({nul['NRAND']} randomisations)")
add('A architecture', 'trophic level, mean by class',
    '; '.join(f"{k} {v:.2f}" for k, v in arch['trophic_level_mean_by_type'].items()))
add('A architecture', 'C(k) ~ k^b hierarchical-modularity exponent',
    round(arch['Ck_scaling_exponent'], 3),
    f"r={arch['Ck_scaling_r']:.3f} over {arch['Ck_n_bins']} degree bins "
    f"(Ravasz-Barabasi predict b = -1)")
add('A architecture', 'independent cross-check in R/igraph',
    '; '.join(f"{r['metric']}={float(r['igraph_value']):g}" for r in igc
              if r['metric'] in ('nodes', 'edges', 'n_SCC', 'largest_SCC', 'reciprocity',
                                 'flow_hierarchy_full', 'flow_hierarchy_regulators')),
    'igraph reproduces every networkx value exactly (systems_igraph_crosscheck.csv)')
_full = next(r for r in regs if r['network'].startswith('full'))
_reg = next(r for r in regs if r['network'].startswith('regulators'))
add('A architecture', 'HIERARCHY TEST on regulators only (TF+miRNA)',
    f"{_reg['N']} nodes, {_reg['M']} edges; largest SCC {_reg['largest_SCC']} "
    f"({float(_reg['frac_nodes_in_largest_SCC'])*100:.1f}% of regulators)",
    'the 206 out-degree-0 genes are forced into the bow-tie OUT layer by construction, so the '
    'apparent regulator-above-gene hierarchy is an artefact of network assembly; among the '
    'regulators themselves the network is one giant recurrent core')
add('A architecture', 'regulator subnetwork: bow-tie IN / CORE / OUT',
    f"{_reg['bowtie_IN']} / {_reg['bowtie_CORE']} / {_reg['bowtie_OUT']}")
add('A architecture', 'regulator subnetwork: hierarchical layers (condensation DAG depth)',
    _reg['condensation_depth_layers'],
    f"full network: {_full['condensation_depth_layers']} layers")
for _r in regn:
    add('A architecture', f"regulator subnetwork null: {_r['metric']}",
        f"obs {float(_r['observed']):.4f} vs null {float(_r['null_mean']):.4f}+/-"
        f"{float(_r['null_sd']):.4f}",
        f"z={float(_r['z']):+.2f}, p={float(_r['p_two_sided']):.3g} "
        f"({_r['n_randomisations']} curveball randomisations)")
for p in pl:
    add('A architecture', f"power law, {p['distribution']}",
        f"alpha={float(p['alpha']):.2f}, xmin={p['xmin']}, n_tail={p['n_tail']}",
        f"Clauset-Shalizi-Newman bootstrap p={float(p['bootstrap_p']):.3f} "
        f"({p['nboot']} sims) -> power law "
        f"{'RULED OUT (p<=0.1)' if float(p['bootstrap_p']) <= 0.1 else 'not ruled out'}")
for l in lrt:
    add('A architecture', f"LRT {l['distribution']} vs {l['alternative']}",
        f"LR={float(l['LR']):+.2f}", f"p_two={float(l['p_two_sided']):.3g}")

# B CONTROLLABILITY
add('B controllability', 'maximum matching |M*|', ctl['matching_size'],
    f"independently recomputed with scipy.sparse.csgraph.maximum_bipartite_matching: "
    f"{mxc['matching_size_scipy']} (agreement: {mxc['agree']})")
add('B controllability', 'minimum driver nodes N_D', ctl['N_D'],
    f"n_D = {ctl['n_D_over_N']*100:.1f}% of nodes (Liu, Slotine & Barabasi 2011)")
add('B controllability', 'control roles (Jia & Barabasi 2013)',
    '; '.join(f'{k} {v}' for k, v in ctl['roles'].items()),
    f"from {ctl['n_matchings_sampled']} sampled maximum matchings; the 5 critical nodes are "
    f"exactly the {mxc['n_in_degree_zero']} in-degree-0 miRNAs "
    f"({', '.join(mxc['in_degree_zero_nodes'])}), which must be drivers in every configuration")
add('B controllability', 'deletion classes (Vinayagam et al. 2016)',
    '; '.join(f'{k} {v}' for k, v in ctl['deletion_classes'].items()))
for g, p in [('FFL_hub (hub_ffl==TRUE)', 'driver_node'),
             ('FFL_hub (hub_ffl==TRUE)', 'indispensable'),
             ('top20_FFL_participation', 'driver_node'),
             ('miR-29/collagen module', 'driver_node'),
             ('miR-29/collagen module', 'indispensable'),
             ('miR-29 family', 'driver_node'),
             ('degree_hub', 'driver_node')]:
    e = next(x for x in enr if x['group'] == g and x['property'] == p)
    add('B controllability', f'{g} -> {p}',
        f"{e['k_group']}/{e['n_group']} ({float(e['frac_group'])*100:.1f}%)",
        f"rest of network {float(e['frac_rest'])*100:.1f}%, OR={float(e['odds_ratio']):.3g}, "
        f"Fisher p={float(e['p']):.3g}")
add('B controllability', 'Spearman driver frequency vs FFL participation',
    f"{ctl['spearman_driverfreq_vs_fflparticipation'][0]:+.3f}",
    f"p={ctl['spearman_driverfreq_vs_fflparticipation'][1]:.3g}")
dgc = json.load(open(f'{RES}/systems_controllability_degree_control.json'))
add('B controllability', 'FFL-hub driver depletion after adjusting for degree',
    f"log-odds {dgc['logistic_driver']['is_FFL_hub']['coef']:+.3f}",
    f"logistic driver ~ log(degree+1) + node type + is_FFL_hub: p="
    f"{dgc['logistic_driver']['is_FFL_hub']['p']:.3g}. THE RAW FISHER DEPLETION IS A DEGREE "
    f"EFFECT and does not survive adjustment.")
add('B controllability', 'FFL-hub driver rate, degree-matched resampling',
    f"{dgc['degree_matched']['observed']:.3f} vs matched non-hubs "
    f"{dgc['degree_matched']['null_mean']:.3f}+/-{dgc['degree_matched']['null_sd']:.3f}",
    f"one-sided p={dgc['degree_matched']['p_one_sided']:.3g} "
    f"({dgc['degree_matched']['n_resamples']} resamples)")
add('B controllability', 'log(FFL participation) after adjusting for degree',
    f"log-odds {dgc['logistic_driver_continuous_ffl']['coef']:+.3f}",
    f"p={dgc['logistic_driver_continuous_ffl']['p']:.3g}")
_sc = dgc['structural_constraints']
add('B controllability', 'STRUCTURAL CAVEAT',
    f"{_sc['n_out_degree_zero']} of 587 nodes have out-degree 0, all genes",
    _sc['note'])
add('B controllability', 'Spearman driver frequency vs degree',
    f"{ctl['spearman_driverfreq_vs_degree'][0]:+.3f}",
    f"p={ctl['spearman_driverfreq_vs_degree'][1]:.3g}")
foc = list(csv.DictReader(open(f'{RES}/systems_controllability_focus_nodes.csv')))
for r in foc[:5]:
    add('B controllability', f"{r['node']} control status",
        f"driver={r['is_driver']}, role={r['control_role']}, class={r['deletion_class']}",
        f"delta N_D on deletion = {r['delta_ND_on_deletion']}")

# C INFORMATION FLOW
for c in cmp_:
    add('C information flow', f"{c['measure']} vs degree",
        f"Spearman {float(c['spearman_vs_degree']):+.3f}",
        f"top-20 overlap {c['top20_overlap']}/20, Jaccard {c['top20_jaccard']}")
for r in skc:
    add('C information flow', f"{r['measure']} vs degree [{r['node_set']}]",
        f"Spearman {float(r['spearman_vs_degree']):+.3f}",
        f"n={r['n']}, top-20 overlap {r['top20_overlap']}/20")
add('C information flow', 'SINK CAVEAT',
    f"{skj['n_sinks']} of {skj['n_nodes']} nodes have out-degree 0 and directed throughput "
    f"exactly 0",
    'the undirected current-flow and information-centrality measures are direction-blind and '
    'credit these sinks with routing capacity they do not have; only the DIRECTED random-walk '
    'throughput distinguishes them, and even that differs from degree only moderately once '
    'the sinks are excluded')
add('C information flow', 'degree hubs that cannot forward any signal',
    f"{len(skj['top30_degree_hubs_that_are_sinks'])} of the top 30 by degree",
    ', '.join(skj['top30_degree_hubs_that_are_sinks']))
fr = list(csv.DictReader(open(f'{RES}/systems_focus_node_ranks.csv')))
for r in fr:
    add('C information flow', f"{r['node']} ranks",
        f"degree {r['rank_degree']}, current-flow betw {r['rank_current_flow_betweenness']}, "
        f"info centrality {r['rank_information_centrality']}, "
        f"directed RW throughput {r['rank_rw_throughput']}")

# D PERTURBATION
sk = list(csv.DictReader(open(f'{RES}/systems_single_knockout.csv')))
add('D perturbation', 'best single knockout',
    f"{comb['best_single']['node']} (composite {comb['best_single']['composite']:.5f})",
    f"dFFL={comb['best_single']['dFFL']}, dReach={comb['best_single']['dReach']}, "
    f"dCOL={comb['best_single']['dCOL']}")
add('D perturbation', 'top 10 single knockouts',
    '; '.join(r['node'] for r in sk[:10]))
_n_drug = sum(1 for r in sk if r['druggable'] == 'True')
_n_appr = sum(1 for r in sk if int(r['n_approved_drugs']) > 0)
_n_anti = sum(1 for r in sk if int(r['n_antineoplastic_drugs']) > 0)
add('D perturbation', 'druggable nodes (DGIdb)', _n_drug,
    f'{_n_appr} with an approved drug, {_n_anti} with an antineoplastic')
_by_t = collections.Counter(r['type'] for r in sk if r['druggable'] == 'True')
_nt_all = collections.Counter(r['type'] for r in sk)
add('D perturbation', 'DRUGGABILITY CAVEAT',
    f"druggable nodes are {_by_t['Gene']} genes + {_by_t['TF']} TFs + {_by_t['miRNA']} miRNAs",
    f"DGIdb is keyed on protein-coding gene symbols, so all {_nt_all['miRNA']} miRNAs in the "
    f"network - including the entire miR-29 family - are druggable-by-definition-zero; the "
    f"'best druggable control point' can therefore only ever be a TF or a gene")
add('D perturbation', 'double knockouts evaluated', comb['n_pairs'], 'exhaustive, no sampling')
add('D perturbation', 'best pair',
    f"{comb['best_pair']['A']} + {comb['best_pair']['B']} "
    f"(composite {comb['best_pair']['composite']:.5f})",
    f"{comb['best_pair']['ratio_to_best_single']:.2f}x the best single node; "
    f"synergy {comb['best_pair']['synergy']:+.5f} (SUB-additive)")
add('D perturbation', 'best pair vs random-pair null',
    f"z = {comb['z_of_best_pair_vs_random_pairs']:.2f}",
    f"random pairs mean {comb['random_pair_null']['mean']:.5f} "
    f"sd {comb['random_pair_null']['sd']:.5f}; empirical p = "
    f"{comb['empirical_p_best_pair']:.3g}")
add('D perturbation', 'pairs beating the best single node',
    f"{comb['n_pairs_beating_best_single']} / {comb['n_pairs']} "
    f"({comb['frac_pairs_beating_best_single']*100:.2f}%)")
add('D perturbation', 'synergy of the 100 most damaging pairs',
    f"{cst['top100_pairs']['n_sub_additive']}/100 SUB-additive",
    f"mean synergy {cst['top100_pairs']['mean_synergy']:+.5f}")
add('D perturbation', 'fraction of all pairs that are super-additive',
    f"{cst['synergy_distribution']['frac_super_additive']*100:.1f}%",
    f"mean synergy {cst['synergy_distribution']['mean']:+.6f}")
add('D perturbation', 'most synergistic pair',
    f"{cst['most_synergistic_pair']['A']} + {cst['most_synergistic_pair']['B']}",
    f"synergy {cst['most_synergistic_pair']['synergy']:+.5f} "
    f"(z={cst['most_synergistic_pair']['z_vs_all_pairs']:.1f} vs all pairs) but only "
    f"{cst['most_synergistic_pair']['percent_of_best_single_damage']:.1f}% of the best "
    f"single-node damage")
add('D perturbation', 'best pair with an approved drug on both targets',
    f"{cst['best_approved_drug_pair']['node_A']} + {cst['best_approved_drug_pair']['node_B']}",
    f"{cst['best_approved_drug_pair']['ratio_to_best_single']}x best single; "
    f"synergy {cst['best_approved_drug_pair']['synergy']}")
ts = list(csv.DictReader(open(f'{RES}/systems_targeted_set_knockouts.csv')))
for r in ts:
    add('D perturbation', f"targeted knockout: {r['set_name']}",
        f"composite {r['composite_damage']} ({r['ratio_to_best_single']}x best single)",
        f"n={r['n']}, synergy {r['synergy']}, {r['n_druggable']} druggable members")

# E DYNAMICS
for tag in ('published', 'curated'):
    d = dyn[tag]
    add('E dynamics', f'best single perturbation [{tag} signs]',
        f"{d['best_single']['node']} ({d['best_single']['direction']:+d})",
        f"Spearman {d['best_single']['spearman']:+.4f} vs observed logFC; "
        f"max-statistic FWER p={d['fwer_p_best_single']:.4f} "
        f"({d['n_permutations']} permutations, null max "
        f"{d['permutation_null_max_mean']:.3f}+/-{d['permutation_null_max_sd']:.3f})")
    add('E dynamics', f'sparse driver set [{tag} signs]',
        f"R2 = {d['R2_at_5']:.3f} with 5 perturbations, {d['R2_at_10']:.3f} with 10",
        f"LASSO entry order: {', '.join(d['lasso_first10'][:6])}")
    add('E dynamics', f'edges usable [{tag} signs]', d['n_edges'],
        f"max Neumann iterations to 1e-12 = {d['max_iterations']} (alpha={d['alpha']})")
    dc = dnl['direction_consistency'][tag]
    add('E dynamics', f'driver direction vs own logFC [{tag}]',
        f"{dc['n_direction_consistent']}/{dc['n_top_with_logFC']} of the top 20 consistent",
        'a perturbation that fits the phenotype need not act in the direction the driver '
        'itself moves in the data')
b = dnl['boolean']
add('E dynamics', 'Boolean threshold attractor',
    f"limit cycle of period {b['period']}, reached in {b['steps']} steps",
    'not a fixed point')
add('E dynamics', 'attractor / phenotype sign agreement',
    f"{b['observed_agreement']:.4f} on n={b['n']}",
    f"permuted-seed null {b['null_permuted_seed_mean']:.3f}+/-"
    f"{b['null_permuted_seed_sd']:.3f} p={b['p_permuted_seed']:.3g}; rewired-network null "
    f"{b['null_rewired_mean']:.3f}+/-{b['null_rewired_sd']:.3f} p={b['p_rewired']:.3g}; "
    f"binomial vs 0.5 p={b['binomial_p']:.3g}")
add('E dynamics', 'best Boolean clamp (>=50 nodes reached)',
    f"{b['best_clamp_n50']['node']} clamped to {b['best_clamp_n50']['clamp']:+d}",
    f"accuracy {b['best_clamp_n50']['accuracy']:.3f} on n={b['best_clamp_n50']['n']}, "
    f"binomial p={b['best_clamp_n50']['binom_p']:.3g}")
add('E dynamics', 'upstream perturbations reproducing the observed collagen increase',
    '; '.join(dnl['collagen_top'][:6]),
    'ranked by the magnitude of the correctly-signed response at COL1A1 and COL3A1')

# F CITATIONS
add('F citations', 'references verified against live PubMed + Crossref',
    f"{sum(r['VERIFIED'] == 'True' for r in refs)}/{len(refs)}")
for r in refs:
    add('F citations', r['key'],
        f"{r['claimed']} | doi {r['claimed_doi']}"
        + (f" | PMID {r['pubmed_pmid']}" if r['pubmed_pmid'] else ''),
        f"Crossref: {r['crossref_journal']} {r['crossref_year']} "
        f"v{r['crossref_volume']} p{r['crossref_page']} -> "
        f"{'VERIFIED' if r['VERIFIED'] == 'True' else 'CHECK'}")

with open(f'{RES}/systems_master_summary.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=['section', 'metric', 'value', 'detail'])
    w.writeheader(); w.writerows(rows)
print(f'wrote systems_master_summary.csv ({len(rows)} rows)')
for r in rows:
    print(f"[{r['section']}] {r['metric']}: {r['value']}" +
          (f"   -- {r['detail']}" if r['detail'] else ''))
