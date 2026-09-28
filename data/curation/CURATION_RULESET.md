# Gene-curation rule set

The rules below are those stated in the manuscript (Methods, Section 2.1). They describe the
procedure; the per-gene decisions that reduced the DisGeNET and GeneCards gene lists to the final
1,057 breast-cancer genes were not recorded in a form that can be re-audited, and inter-curator
agreement was not recorded (both are stated as limitations in the paper).

## Source thresholds (applied before curation)

- DisGeNET v7.0: gene–disease association score ≥ 0.4 (the score combines the number and curation
  type of the reporting sources, animal-model data and the number of supporting publications).
  1,612 genes.
- GeneCards v5.19: relevance score ≥ 0.4 after min–max rescaling of the returned set. 5,247 genes.

## Manual curation (two curators, applied independently)

1. A gene was retained if its association was supported by at least one primary experimental report
   in breast tissue or a breast cancer cell line.
2. A gene was excluded if it was solely a germline risk locus with no reported somatic or expression
   phenotype.

Disagreements were resolved by the senior author. Curators: G.P.K.R. and J.M.A.; senior author F.S.M.

Result: 1,057 genes. Transcription-factor status was not curated by hand: it was assigned by
intersection with the Lambert et al. (2018) human TF census (665 of the 1,057 genes).
