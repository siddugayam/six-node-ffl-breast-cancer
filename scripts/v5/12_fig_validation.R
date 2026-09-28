# Figure: coherence typology, expression validation by evidence tier, and the mediation result.
suppressMessages({library(data.table); library(ggplot2); library(patchwork); library(scales)})
REV <- "/path/to/revision"
source(file.path(REV,"scripts/v5/10_theme.R")); FIG <- file.path(REV,"figures/v5")

# ---- A: coherence typology ------------------------------------------------------------
co <- fread(file.path(REV,"results/v2/ffl_cores_coherence_corrected.csv"))
cc <- co[coherence!="unresolved", .N, by=.(class, coherence)]
cc[, kind := ifelse(substr(coherence,1,1)=="C","coherent","incoherent")]
tot <- cc[, .(N=sum(N)), by=.(coherence, kind)][order(-N)]
tot[, coherence := factor(coherence, levels=rev(tot$coherence))]
pA <- ggplot(tot, aes(N, coherence, fill=kind)) +
  geom_col(width=0.66) +
  geom_text(aes(label=N), hjust=-0.18, size=2.6, colour=INK) +
  scale_fill_manual(values=c(coherent=unname(SEM["pos"]), incoherent=unname(SEM["neg"])), name=NULL) +
  scale_x_continuous(expand=expansion(mult=c(0,0.16))) +
  labs(x="3-node cores", y=NULL, title="a   Coherence typology of the 3-node cores",
       subtitle="898 of 1,649 cores resolved\nI1, the classical pulse generator, is commonest") +
  theme_pub() + theme(legend.position="top", legend.justification="left",
                      panel.grid.major.y=element_blank())

cls <- cc[, .(coherent=sum(N[kind=="coherent"]), total=sum(N)), by=class]
cls[, pct := 100*coherent/total]
cls[, lab := sprintf("%s\n(n = %d)", class, total)]
cls[, lab := factor(lab, levels=cls[order(pct)]$lab)]
pB <- ggplot(cls, aes(pct, lab)) +
  geom_col(fill=SEM["accent"], width=0.6, alpha=0.9) +
  geom_vline(xintercept=50, linetype="dashed", colour=MUTED, linewidth=0.4) +
  geom_text(aes(label=sprintf("%.1f%%", pct)), hjust=-0.25, size=2.7, colour=INK) +
  scale_x_continuous(limits=c(0,105), breaks=seq(0,90,30), expand=expansion(mult=c(0,0.10))) +
  labs(x="% of resolved cores that are coherent", y=NULL,
       title="b   The motif classes differ functionally",
       subtitle="miRNA-led circuits are mostly coherent\ncomposite circuits mostly incoherent") +
  theme_pub() + theme(panel.grid.major.y=element_blank())

# ---- C: concordance by class and evidence tier ----------------------------------------
sc <- fread(file.path(REV,"results/sign_concordance_summary.csv"))
sc <- sc[stratum!="ALL_sign_annotated"]
nice <- c(miRNA_target_ALL="miRNA -| target, all",
          miRNA_target_strong="miRNA -| target, strong",
          miRNA_target_weak="miRNA -| target, weak (CLIP)",
          miRNA_target_predonly="miRNA -| target, predicted only",
          TF_target_Activation="TF -> target, activation",
          TF_target_Repression="TF -> target, repression",
          TF_miRNA_Activation="TF -> miRNA, activation",
          TF_miRNA_Repression="TF -> miRNA, repression")
sc[, lab := nice[stratum]]
sc <- sc[!is.na(lab)]
sc[, sig := binom_p_vs_null < 0.05]
sc[, lab := factor(lab, levels=sc[order(concordance_rate - null_concordance_rate)]$lab)]
long <- melt(sc[, .(lab, observed=100*concordance_rate, null=100*null_concordance_rate, sig)],
             id.vars=c("lab","sig"), variable.name="what", value.name="pct")
pC <- ggplot(long, aes(pct, lab)) +
  geom_line(aes(group=lab), colour=RULE, linewidth=1.6) +
  geom_point(aes(colour=what, shape=what), size=2.4) +
  geom_text(data=sc, aes(x=100*pmax(concordance_rate,null_concordance_rate)+1.4, y=lab,
                         label=ifelse(sig, sprintf("p = %s", signif(binom_p_vs_null,2)), "n.s.")),
            hjust=0, size=2.3, colour=INK) +
  scale_colour_manual(values=c(observed=unname(SEM["accent"]), null=unname(SEM["null"])),
                      labels=c("observed","expression-matched null"), name=NULL) +
  scale_shape_manual(values=c(observed=16, null=1),
                     labels=c("observed","expression-matched null"), name=NULL) +
  scale_x_continuous(limits=c(35,82), expand=expansion(mult=c(0.02,0.16))) +
  labs(x="sign concordance (%)", y=NULL,
       title="c   Only transcriptional activation predictions validate",
       subtitle="Predicted edges vs random node pairs matched on the expression decile of both partners, 1,066 tumours.\nPrediction-only miRNA edges fall BELOW their null - they correlate in the direction opposite to the model.") +
  theme_pub() + theme(legend.position="top", legend.justification="left",
                      panel.grid.major.y=element_blank())

# ---- D: mediation ----------------------------------------------------------------------
med <- data.table(
  path = c("ETS1 -> COL1A1","NFKB1 -> COL1A1","SP1 -> COL1A1","RELA -> COL1A1",
           "miR-29a -| COL3A1","miR-29a -| COL1A1","miR-29b -| COL3A1"),
  total= c(0.301, 0.106, 0.039, 0.051, -0.237, -0.230, -0.286),
  acme = c(0.724, 0.242, 0.094, 0.006, 0.009, 0.008, -0.132),
  lo   = c(0.662, 0.195, 0.054,-0.033,-0.037,-0.034,-0.178),
  hi   = c(0.794, 0.288, 0.136, 0.048, 0.057, 0.053,-0.087),
  ade  = c(-0.423,-0.136,-0.055, 0.045,-0.246,-0.239,-0.154))
med[, kind := ifelse(grepl("miR", path), "miRNA arm", "TF arm")]
med[, path := factor(path, levels=rev(med$path))]
pD <- ggplot(med, aes(y=path)) +
  geom_vline(xintercept=0, colour=RULE, linewidth=0.5) +
  geom_errorbarh(aes(xmin=lo, xmax=hi), height=0.16, colour=MUTED, linewidth=0.5) +
  geom_point(aes(x=acme, colour=kind), size=2.8) +
  geom_point(aes(x=ade), shape=4, size=2.2, colour=INK, stroke=0.7) +
  scale_colour_manual(values=c(`TF arm`=unname(SEM["neg"]), `miRNA arm`=unname(SEM["pos"])), name=NULL) +
  labs(x="effect size", y=NULL,
       title="d   Fibroblast content mediates the TF arm, not the miR-29a arm",
       subtitle="Filled circle = mediated effect (ACME) with 95% bootstrap CI;  x = direct effect (ADE). 1,000 resamples.\nFor ETS1 and NFKB1 the mediated effect exceeds the total and the direct effect flips sign - inconsistent mediation.\nFor miR-29a the mediated effect is indistinguishable from zero: the relationship is direct.",
       caption="CAF covariate: 123-gene stromal signature with all collagens removed; neither it nor the marker-panel alternative contains a network node.") +
  theme_pub() + theme(legend.position="top", legend.justification="left",
                      panel.grid.major.y=element_blank())

out <- (pA | pB) / pC / pD + plot_layout(heights=c(0.95, 1, 1.2))
ggsave(file.path(FIG,"Fig4_validation_and_mediation.png"), out, width=11, height=12.6, dpi=400, bg="white")
ggsave(file.path(FIG,"Fig4_validation_and_mediation.pdf"), out, width=11, height=12.6, bg="white")
cat("wrote Fig4_validation_and_mediation\n")
