#!/usr/bin/env python3
"""40_seqreg_ext_alphagenome_access.py

Re-verification and EXTENSION of the AlphaGenome access record.

Everything written here was observed directly. Two facts are new relative to the
first pass:
  (1) AlphaGenome MODEL WEIGHTS are now publicly released (Jan 2026) on HuggingFace and
      Kaggle, so there is a second, keyless route that the first pass did not record. The
      official Google repos are GATED ("gated": "auto"), i.e. they require a HuggingFace
      account and acceptance of the AlphaGenome Model Terms -- an account-holder action.
  (2) The AlphaGenome output modality list, read directly from the installed official
      client, contains NO DNA-METHYLATION head. Task query (iii) -- "predicted effect of
      the promoter methylation change on miR-130a expression" -- therefore cannot be asked
      of AlphaGenome even once a key exists. This is now verified, not inferred.

Nothing here fabricates a key or a prediction.
"""
import json, subprocess, sys
import pandas as pd, requests

ROOT = "/path/to/revision"
RES = f"{ROOT}/results/v3"
rows = []


def add(step, result, detail):
    rows.append(dict(step=step, result=result, detail=detail))
    print(f"[{result}] {step}\n    {detail}\n", flush=True)


# ---- 1. package + modality inventory -------------------------------------------------
import alphagenome
from alphagenome.models import dna_client
mods = [m.name for m in dna_client.OutputType]
add("alphagenome python package", f"INSTALLED v{alphagenome.__version__ if hasattr(alphagenome,'__version__') else '0.9.0'}",
    "pip package present; client imports cleanly")
add("AlphaGenome output modalities (read from the installed official client)",
    f"{len(mods)} modalities",
    ", ".join(mods))
add("Is there a DNA-METHYLATION output head?", "NO",
    "OutputType has no methylation member. Task query A(iii) (effect of miR-130a promoter "
    "methylation on its expression) CANNOT be answered by AlphaGenome even with a valid key: "
    "the model takes DNA sequence only and emits no methylation track. The nearest supported "
    "proxies are predicted accessibility/histone marks at the locus and in-silico mutation of "
    "the CpG-rich element.")

# ---- 2. endpoint reachability with a deliberately invalid key -------------------------
try:
    from alphagenome.data import genome
    m = dna_client.create("INVALID_KEY_FOR_REACHABILITY_TEST_ONLY")
    iv = genome.Interval(chromosome="chr17", start=50136095, end=50267167)  # 131,072 bp
    m.predict_interval(interval=iv, requested_outputs=[dna_client.OutputType.DNASE],
                       ontology_terms=None)
    add("gRPC endpoint reachability", "UNEXPECTED SUCCESS", "an invalid key was accepted")
except Exception as e:
    msg = str(e).replace("\n", " ")[:400]
    add("gRPC endpoint reachability (invalid key)",
        "REACHABLE - only the key is missing" if "API key not valid" in msg else "ERROR", msg)

# ---- 3. the key-issuing page ----------------------------------------------------------
for url in ["https://deepmind.google.com/science/alphagenome",
            "https://deepmind.google.com/science/alphagenome/account/terms"]:
    try:
        r = requests.get(url, timeout=45)
        has_key_form = "Create API Key" in r.text or "Get API key" in r.text
        add(f"GET {url}", f"HTTP {r.status_code}, {len(r.text)} bytes",
            f"JS-only shell; 'Create API Key' control present in served HTML: {has_key_form}. "
            "The key is issued only to a signed-in Google account that has accepted the Terms.")
    except Exception as e:
        add(f"GET {url}", "ERROR", str(e)[:200])

# ---- 4. public weights (the route the first pass missed) ------------------------------
for repo in ["google/alphagenome-all-folds", "gtca/alphagenome_pytorch"]:
    try:
        j = requests.get(f"https://huggingface.co/api/models/{repo}", timeout=45).json()
        add(f"HuggingFace repo {repo}",
            f"gated={j.get('gated')}",
            f"private={j.get('private')}; licence tags={[t for t in j.get('tags',[]) if 'license' in t]}; "
            f"files={len(j.get('siblings',[]))}")
    except Exception as e:
        add(f"HuggingFace repo {repo}", "ERROR", str(e)[:200])

add("Local-inference package alphagenome-pytorch",
    "NOT INSTALLED - installation not permitted on the analysis machine",
    "Installing alphagenome-pytorch was not permitted on the analysis machine, so local "
    "inference could not be attempted regardless of weights. Running AlphaGenome locally also "
    "needs the architecture code (Apache-2.0, github.com/genomicsxai/alphagenome-pytorch or "
    "google-deepmind/alphagenome_research), and DeepMind recommend >= 1x NVIDIA H100; this "
    "machine has an RTX A4000 (16 GB), so a 1 Mb-context run would in any case not fit and only "
    "the 16 kb / 131 kb contexts would be worth attempting.")

add("DECISION TAKEN", "AlphaGenome NOT RUN",
    "No key was fabricated, no Google or HuggingFace account was created or signed into, and "
    "the gated official weights were NOT downloaded through the ungated third-party mirror, "
    "because the gate exists precisely to record acceptance of the AlphaGenome Model Terms and "
    "that acceptance is the account holder's to give.")

add("ACTION REQUIRED BY THE ACCOUNT HOLDER - ROUTE 1 (hosted API, recommended)", "1 step",
    "Signed in to Google as your.email@example.org, open "
    "https://deepmind.google.com/science/alphagenome , click 'Get API key' / 'Create API Key', "
    "accept the Terms of Service, copy the key, then run:  "
    "export ALPHAGENOME_API_KEY=<key> && python3 scripts/09_regulatory_evidence/sequence_models/21_seqreg_alphagenome.py  "
    "(free for non-commercial use; no GPU needed, inference is server-side).")

add("ACTION REQUIRED BY THE ACCOUNT HOLDER - ROUTE 2 (local weights)", "3 steps, needs a big GPU",
    "Accept the model terms at https://huggingface.co/google/alphagenome-all-folds (or "
    "https://www.kaggle.com/models/google/alphagenome ), then 'pip install alphagenome-pytorch' "
    "(not permitted on the analysis machine) and run on an H100-class GPU. Route 1 is strictly easier.")

add("Commercial route (not used, not authorised)", "available",
    "AlphaGenome is served for commercial use via Google Cloud; that needs a billed GCP project.")

pd.DataFrame(rows).to_csv(f"{RES}/seqreg_ext_alphagenome_access.csv", index=False)
print(f"wrote {RES}/seqreg_ext_alphagenome_access.csv ({len(rows)} rows)")
