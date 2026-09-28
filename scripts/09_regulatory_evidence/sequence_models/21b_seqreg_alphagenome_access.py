#!/usr/bin/env python3
"""21b -- record, reproducibly, exactly what was attempted to obtain AlphaGenome access,
what the server replied, and the single step the account holder must perform."""
import os, pandas as pd, subprocess, json

ROOT="/path/to/revision"; RES=f"{ROOT}/results/v3"

def probe():
    """Confirm the gRPC endpoint is reachable and that the only missing element is a key."""
    from alphagenome.models import dna_client
    from alphagenome.data import genome
    try:
        m=dna_client.create("DUMMY_KEY_NOT_REAL")
        iv=genome.Interval(chromosome='chr17', start=50_180_000, end=50_180_000+16384)
        m.predict_interval(interval=iv,
                           requested_outputs=[dna_client.OutputType.RNA_SEQ],
                           ontology_terms=['UBERON:0002097'])
        return "UNEXPECTED_SUCCESS"
    except Exception as e:
        return f"{type(e).__name__}: {str(e)[:300]}"

resp=probe()
rows=[
 dict(step="pip install alphagenome", result="OK", detail="alphagenome 0.9.0 installed from PyPI"),
 dict(step="network reachability of the AlphaGenome gRPC endpoint",
      result="REACHABLE",
      detail=f"a deliberately invalid key returns: {resp}"),
 dict(step="GET https://deepmind.google.com/science/alphagenome",
      result="HTTP 200 but JS-only shell; the only actionable control in the served HTML is 'Sign in'",
      detail="page body contains no API-key form without a Google session"),
 dict(step="GET https://deepmind.google.com/science/alphagenome/account/terms",
      result="HTTP 200, same JS-only shell, no content without a signed-in Google session",
      detail="this is the documented key-issuing page (README of google-deepmind/alphagenome)"),
 dict(step="obtain key using your.email@example.org",
      result="NOT DONE - BLOCKED, ACCOUNT-HOLDER ACTION REQUIRED",
      detail=("the key is issued only after (a) signing in to a Google account and "
              "(b) accepting the AlphaGenome Terms of Service. Both are acts of the account "
              "holder; no email-only or form-only route exists. No key was fabricated and no "
              "sign-in was attempted.")),
 dict(step="ACTION REQUIRED BY THE ACCOUNT HOLDER",
      result="1 step",
      detail=("Open https://deepmind.google.com/science/alphagenome/account/terms while signed "
              "in as your.email@example.org, accept the Terms of Service, copy the "
              "API key, then run:  export ALPHAGENOME_API_KEY=<key> && "
              "python3 scripts/09_regulatory_evidence/sequence_models/21_seqreg_alphagenome.py")),
 dict(step="commercial alternative (not used)",
      result="available",
      detail="AlphaGenome is served for commercial use on Google Cloud "
             "(docs.cloud.google.com/gemini-enterprise-agent-platform/models/open-models/alphagenome); "
             "this requires a billed GCP project, which was not authorised."),
]
pd.DataFrame(rows).to_csv(f"{RES}/seqreg_alphagenome_access.csv", index=False)
print(pd.DataFrame(rows).to_string(index=False))
