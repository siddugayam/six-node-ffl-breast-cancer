#!/usr/bin/env python3
"""Part A deliverable: the full search + triage inventory of breast-cancer miRNA
cohorts, with the outcome fields actually present in GEO and the usability verdict."""
import json, os, re, csv, collections
BASE="/path/to/revision"
CACHE=BASE+"/cache/v6/cohorts"; OUT=BASE+"/results/v6"
soft=CACHE+"/soft"
allr={r["accession"]:r for r in json.load(open(CACHE+"/triage_all.json"))}
cand=json.load(open(CACHE+"/triage_candidates.json"))

# verdicts I established by reading each series' deposited sample characteristics
AUTO=json.load(open(CACHE+"/autoclass.json")) if os.path.exists(CACHE+"/autoclass.json") else {}
VERDICT = {
 "GSE22220":("SUPERSERIES of GSE22216 + GSE22219","the miRNA arm (GSE22216) is used for Cox; the mRNA arm (GSE22219) is already in the paper","do not count twice"),
 "GSE59595":("NOT USABLE","the only numeric field is 'age (years)'; no survival time deposited","Chinese/Italian series"),
 "GSE59594":("NOT USABLE","miRNA arm of GSE59595; age only, no survival time",""),
 "GSE59590":("NOT USABLE","mRNA arm of GSE59595; age only, no survival time",""),
 "GSE59246":("NOT USABLE","subseries of GSE59248; follow-up time in a minority of samples, <=12 events",""),
 "GSE59247":("NOT USABLE","subseries of GSE59248; follow-up time in a minority of samples, <=12 events",""),
 "GSE31309":("NOT USABLE","only 'year of birth' is numeric; no outcome",""),
 "GSE75669":("NOT USABLE","subset of GSE75685; dates only, no baseline date",""),
 "GSE59248":("NOT USABLE","follow-up time present for only 49 of 153 samples, 4 DMFS / 10 recurrence / 12 death events","DCIS-to-IBC superseries (GSE59246+GSE59247); too few events for Cox"),
 "GSE75685":("NOT USABLE","date-at-relapse / date-at-death only, no baseline date to build a time from",""),
 "GSE68373":("LIQUID BIOPSY - binary only","metastasis class (metastasis / no metastasis) in 43 plasmas","plasma, n=64"),
 "GSE37407":("NOT USABLE","time to recurrence for only 15 samples; primary/metastasis pair design",""),
 "GSE22216":("USED - Cox","DRFS time (yr) + distant-relapse event, complete for 210/210","Buffa/Oxford early primary BC"),
 "GSE37405":("USED - Cox","time.to.recurrence / time.to.death / time.to.ok (yr)","DBCG high-risk ER+ on tamoxifen; 3 array generations, platform as stratum"),
 "GSE19783":("USED - Cox (already in the paper)","death status + disease-free-survival time (months)","Oslo MicMa; the paper's existing n=101 miRNA cohort"),
 "GSE19536":("DUPLICATE of GSE19783","same 101 GPL8227 samples","superseries/duplicate deposit"),
 "GSE59829":("USED - Cox (NEW in v6)","DFS (months) + distant_metastasis event, complete for 123/123","Illumina v2 miRNA beadchip"),
 "GSE78870":("USED - Cox (NEW in v6)","time-to-progression + progression event; disease-free interval + relapse event","ER+ advanced BC, first-line aromatase inhibitor; TLDA cards"),
 "GSE57897":("USED - logistic only","binary 'survival' flag (117/422 events); NO follow-up time deposited","largest tissue miRNA series with any outcome"),
 "GSE103161":("USED - logistic only","metastasis yes/no + dead/alive; NO follow-up time","systemically untreated LN- ER+, matched pairs"),
 "GSE126125":("USED - logistic only (subset design)","metastasis yes/no; NO time","same design/platform as GSE103161, likely overlapping patients"),
 "GSE97811":("USED - logistic only","5-year recurrence Yes/No; NO time","women <35"),
 "GSE28321":("USED - logistic only","rfs status flag; NO time","Ambion bioarray"),
 "GSE40267":("USED - logistic only","death_dd_after_diagnosis + cause_of_death, but NO censoring time for the 87 survivors","TNBC; cannot be Cox-modelled as deposited"),
 "GSE26659":("USED - logistic only (already in the paper)","72-month relapse flag; NO time","= GSE26666, identical 94 GSMs"),
 "GSE26666":("DUPLICATE of GSE26659","identical 94 GSMs verified","companion deposit"),
 "GSE45498":("NOT USABLE","'status' is a tissue label (Metastatic/Primary), no survival fields deposited","TNBC; survival is in the paper, not in GEO"),
 "GSE41970":("NOT USABLE","same as GSE45498","TNBC"),
 "GSE81002":("NOT USABLE","only ER/HER2/grade/PAM50 deposited; no survival fields","Aure et al. luminal-A split; outcome not in GEO"),
 "GSE81000":("NOT USABLE","no survival fields",""),
 "GSE80999":("NOT USABLE","no survival fields",""),
 "GSE40525":("NOT USABLE","no outcome fields (tumour/normal design)","already assembled in v4 for DE only"),
 "GSE68085":("NOT USABLE","no outcome fields deposited (tissue type/subtype/batch only)","named in the task as 'SCAN-B miRNA' - it is not; 114 NGS samples, Norway"),
 "GSE267543":("NOT USABLE","only tissue/gender deposited",""),
 "GSE131599":("NOT USABLE","ER/PR/HER2 only",""),
 "GSE100769":("NOT USABLE","ER/HER2 only",""),
 "GSE225292":("NOT USABLE","no outcome fields",""),
 "GSE86277":("NOT USABLE","no outcome fields","TNBC subtype series named in the task"),
 "GSE86278":("NOT USABLE","no outcome fields","TNBC subtype series named in the task"),
 "GSE86281":("NOT USABLE","no outcome fields","TNBC subtype series named in the task"),
 "GSE70754":("NOT USABLE","no outcome fields deposited",""),
 "GSE39543":("NOT USABLE","'OS time' present for 46 but NO vital-status/event field",""),
 "GSE162670":("NOT USABLE","LCM morphology design, no outcome",""),
 "GSE22049":("NOT USABLE","6 platforms, no outcome fields",""),
 "GSE37963":("NOT USABLE","no outcome fields",""),
 "GSE93624":("NOT BREAST","pediatric Crohn disease ileal transcriptome (GPL11154, n=245)","named in the task; the accession is not a breast miRNA series"),
 "GSE59993":("LIQUID BIOPSY - no case/control or outcome usable","'time from surgery' only; circulating miRNA","named in the task"),
 "GSE73002":("USED - liquid biopsy","diagnosis: breast cancer 1280 / non-cancer 2686 / benign breast 54 / prostate 93","serum, 3D-Gene; no survival"),
 "GSE44281":("USED - liquid biopsy","case/control status, prospectively collected (Sister Study)","serum"),
 "GSE110651":("USED - liquid biopsy","new distant metastasis Yes/No","serum, metastatic BC on eribulin"),
 "GSE118782":("USED - liquid biopsy","dfs_time + dfs_status for the cancer plasmas","plasma, n=40"),
}
rows=[]
OUTCOME_RE=re.compile(r"surviv|dfs|rfs|\bos\b|relapse|recur|event|follow.?up|death|dead|alive|"
  r"dmfs|drfs|\bmfs\b|\bpcr\b|response|outcome|metasta|censor|prognos|time|status",re.I)
seen=set()
for r in sorted(cand, key=lambda x:-x["n_samples"]) + [allr[a] for a in VERDICT if a in allr]:
    acc=r["accession"]
    if acc in seen: continue
    seen.add(acc)
    p=os.path.join(soft,acc+".soft.txt")
    keys=[]; plats=[]; nprobed=""
    if os.path.exists(p):
        t=open(p,encoding="utf-8",errors="replace").read()
        c=collections.Counter()
        for line in t.splitlines():
            if line.startswith("!Sample_characteristics_ch1"):
                v=line.split("=",1)[1].strip()
                c[v.split(":")[0].strip() if ":" in v else v[:40]]+=1
        keys=[k for k in c if OUTCOME_RE.search(k)]
        plats=sorted(set(re.findall(r"!Sample_platform_id = (GPL\d+)",t)))
        nprobed=len(re.findall(r"^\^SAMPLE",t,re.M))
    a=AUTO.get(acc)
    lab = {"no outcome field deposited":"NOT USABLE - no outcome field deposited in GEO",
           "event only (logistic possible)":"NOT USED - binary outcome only, below the priority cut (n or design)",
           "time+event (Cox possible)":"REVIEWED - see outcome_detail"}
    default = (lab.get(a["auto_class"], a["auto_class"]), "auto-classified from the deposited sample characteristics", "") if a else \
              ("not probed (n<40 or not a breast miRNA series)","","")
    v=VERDICT.get(acc, default)
    rows.append(dict(accession=acc, n_geo=r["n_samples"], n_probed=nprobed,
        title=r["title"], platforms=";".join(plats) if plats else r["gpl"],
        outcome_fields_in_GEO="; ".join(keys),
        auto_class=(AUTO.get(acc) or {}).get("auto_class",""),
        auto_time_fields=(AUTO.get(acc) or {}).get("time_fields",""),
        auto_event_fields=(AUTO.get(acc) or {}).get("event_fields",""),
        verdict=v[0], outcome_detail=v[1], note=v[2],
        search_queries=r.get("queries","")))
with open(OUT+"/cohorts_mirna_inventory.csv","w",newline="") as f:
    w=csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader()
    for r in rows: w.writerow(r)
print("wrote", len(rows), "rows")
vc=collections.Counter(r["verdict"].split(" - ")[0].split(" (")[0] for r in rows)
for k,v in vc.most_common(): print("%5d  %s"%(v,k))
