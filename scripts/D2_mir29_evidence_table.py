#!/usr/bin/env python3
"""D) Curated PubMed literature table: miR-29 perturbation -> collagen genes.
Every row was read from the PubMed record (title + abstract) retrieved by
scripts/D_mir29_literature.py. This is LITERATURE evidence, not our own analysis."""
import csv, json, urllib.parse, urllib.request, os
OUT = "/path/to/revision/results/multiomics"
E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

REC = [
 # pmid, mirna, perturbation, system, readout, direction, direct_target_assay, note
 ("18723672","miR-29a/b/c","anti-miR (inhibition) and over-expression","mouse heart / cardiac fibroblasts",
  "COL1A1, COL1A2, COL3A1 (and elastin, fibrillins)","inhibition -> collagen UP; over-expression -> collagen DOWN","yes (3'UTR reporter)",
  "Foundational paper: 'down-regulation of miR-29 would ... induce the expression of collagens, whereas over-expression of miR-29 in fibroblasts reduces collagen expression'"),
 ("18390668","miR-29c","re-expression in nasopharyngeal carcinoma cells","NPC cell lines / patient tumours",
  "multiple collagens and laminin gamma1","miR-29c LOW in tumour -> ECM mRNAs UP; restoring miR-29c represses them","yes (reporter)",
  "miR-29c down-regulated in NPC, up-regulating mRNAs encoding extracellular matrix proteins"),
 ("30472058","miR-29 mimic (remlarsen, MRG-201)","intradermal miR-29 mimic, human phase-2 skin study",
  "human skin, incisional wounds","collagen genes / fibroplasia","mimic -> collagen expression DOWN","no (clinical/expression)",
  "First-in-class miR-29 mimic in humans; remlarsen repressed collagen expression and development of fibroplasia in skin"),
 ("31034464","miR-29a-3p","mimic and inhibitor in nasopharyngeal carcinoma cells","NPC cell lines, radiotherapy",
  "COL1A1","miR-29a UP -> COL1A1 DOWN","yes (COL1A1 3'UTR luciferase)",
  "COL1A1 confirmed a direct target of miR-29a; COL1A1 knockdown phenocopies miR-29a for radiosensitivity"),
 ("30518744","miR-29a/b/c family","over-expression in methotrexate-resistant osteosarcoma lines","MG-63/MTX, U2OS/MTX",
  "COL3A1 (and MCL1)","miR-29 over-expression -> COL3A1 mRNA and protein DOWN","yes (target-relationship assay)",
  "Over-expression of miR-29 family decreased COL3A1 and MCL1 mRNA and protein; COL3A1 over-expression rescued the phenotype"),
 ("35256956","miR-29a-3p","exosome-/liposome-based nanovesicle delivery, in vivo","lung metastatic niche, tumour-bearing mice",
  "collagen I","miR-29a-3p delivery -> collagen I secretion DOWN, lung colonisation reduced","no (in vivo delivery)",
  "Clinical lung tumour data show miR-29a-3p expression negatively correlates with collagen I; delivery down-regulated collagen I secretion by lung fibroblasts in vivo"),
 ("41900803","miR-29a-3p, miR-29b-3p","miRNA mimic transfection of osteoblasts; miRNAs carried by MDA-MB-231 breast-cancer EVs",
  "mouse osteoblasts; MDA-MB-231 and 4T1 breast cancer EVs","Col1a1, Col1a2",
  "miR-29a/29b mimic -> collagen 1a1 and 1a2 mRNA DOWN","no (mimic transfection)",
  "Breast-cancer-derived extracellular vesicles shuttle miR-29a/b to bone cells; mimics reduced collagen1a1 and 1a2 mRNA and mineralisation"),
 ("27233758","miR-29a/b/c family","over-expression and knockdown in leiomyoma / myometrial cells","human uterine leiomyoma vs myometrium",
  "fibrillar collagens type I, II, III","miR-29 knockdown -> collagen III UP; miR-29 over-expression -> fibrillar collagens DOWN","no",
  "Down-regulation of miR-29 species in myometrium increases collagen type III; when miR-29 members are over-expressed the major fibrillar collagens decrease"),
 ("34491311","miR-29 family","lncRNA MIAT sponging of miR-29 family; TGF-beta3 induction","human leiomyoma smooth-muscle cell spheroids",
  "COL1A1, COL3A1","MIAT sequesters miR-29 -> COL1A1/COL3A1 UP; miR-29 restoration blocks induction","no",
  "TGF-beta3 induced COL1A1, COL3A1 and MIAT; miR-29 mediates the effect on COL1A1 and COL3A1 induction"),
 ("32125507","miR-29","pharmacological reduction of scleral miR-29 by genipin","guinea-pig myopia model, sclera",
  "COL1A1","miR-29 DOWN -> COL1A1 UP","no",
  "Genipin inhibits scleral expression of miR-29 and MMP2 and promotes COL1A1 expression - reverse-direction confirmation"),
 ("27490049","miR-29b","over-expression in human corneal endothelial cells","primary human corneal endothelial cells",
  "COL1A1, COL4A1, LAMC1","miR-29b over-expression -> COL1A1 mRNA to 1.9% of control; protein reduced","no",
  "Quantitative ECM reduction after miR-29b over-expression"),
 ("36206860","miR-29a","mimic and inhibitor","human fetal scleral fibroblasts",
  "COL1A1 (via HSP47/Smad3)","mimic -> COL1A1 DOWN; inhibitor -> collagen production UP","no",
  "Both directions tested in the same system"),
 ("29230533","miR-29","expression modulation after ionizing radiation","irradiated fibroblasts",
  "type I collagen genes","miR-29 targets type I collagen genes via 3'UTR binding","yes (3'UTR binding)",
  "Radiation-induced type I collagen increase is miR-29 dependent"),
 ("23467423","miR-29c","restoration by HIF-alpha activation","rat remnant kidney, human IgA nephropathy, HK2 cells",
  "extracellular-matrix genes, COL2A1","miR-29c DOWN -> ECM/collagen UP; restoration attenuates fibrosis","yes (COL2A1 direct target)",
  "The miR-29 family directly targets a large number of extracellular matrix genes"),
]

def summaries(pmids):
    u = E+"esummary.fcgi?"+urllib.parse.urlencode({"db":"pubmed","id":",".join(pmids),"retmode":"json"})
    return json.loads(urllib.request.urlopen(u, timeout=120).read().decode())["result"]

s = summaries([r[0] for r in REC])
rows = []
for r in REC:
    it = s.get(r[0], {})
    rows.append(dict(pmid=r[0],
        year=(it.get("pubdate") or "")[:4], journal=it.get("source",""),
        title=it.get("title","").strip(),
        doi=next((a["value"] for a in it.get("articleids",[]) if a["idtype"]=="doi"), ""),
        mirna=r[1], perturbation=r[2], model_system=r[3], readout_genes=r[4],
        reported_direction=r[5], direct_target_assay=r[6], evidence_note=r[7],
        evidence_class="published literature (not our analysis)",
        supports_miR29_represses_collagen=True))
with open(os.path.join(OUT,"mir29_collagen_literature.csv"),"w",newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader()
    for r in rows: w.writerow(r)
print("wrote mir29_collagen_literature.csv rows:", len(rows))
for r in rows:
    print(f"  PMID {r['pmid']} ({r['year']}) {r['journal']}: {r['title'][:90]}")
    print(f"      {r['mirna']} | {r['perturbation'][:60]} | {r['readout_genes'][:55]} | {r['reported_direction'][:70]}")
