# Figure 2: coherence typology, expression validation by evidence tier, and the mediation result.
# Revised 2026-09-25: mediation summarised over the 42 distinct stromal estimates (three runs that
# reproduced other estimates exactly are excluded), Unicode arrows, italic gene symbols, p values
# formatted as in Table 5, in-panel text at the tick-label size, value labels on a white ground so
# that reference and grid lines do not cross them, and a signed effect-size axis in panel d. Run in a UTF-8 locale (LC_ALL=en_US.UTF-8).
suppressMessages({library(data.table); library(ggplot2); library(patchwork); library(scales); library(ragg); library(ggtext)})
REV <- "/path/to/revision"
source(file.path(REV,"scripts/15_figures_and_tables/10_theme.R")); FIG <- file.path(REV,"figures/final")
co <- fread(file.path(REV,"results/v2/ffl_cores_coherence_corrected.csv"))
TXT <- 8 / .pt  # drawn at print width (174 mm): tick labels, legends and in-panel text 8 pt, axis titles 9 pt
cc <- co[coherence!="unresolved", .N, by=.(class, coherence)]
cc[, kind := ifelse(substr(coherence,1,1)=="C","coherent","incoherent")]
tot <- cc[, .(N=sum(N)), by=.(coherence, kind)][order(-N)]
tot[, coherence := factor(coherence, levels=rev(tot$coherence))]
pA <- ggplot(tot, aes(N, coherence, fill=kind)) + geom_col(width=0.66) +
  geom_text(aes(label=N), hjust=-0.18, size=TXT, colour=INK) +
  scale_fill_manual(values=c(coherent=unname(SEM["pos"]), incoherent=unname(SEM["neg"])), name=NULL) +
  scale_x_continuous(expand=expansion(mult=c(0,0.16))) + labs(x="3-node cores", y=NULL, tag="a") +
  theme_pub(9.5) + theme(legend.position="top", legend.justification="left", panel.grid.major.y=element_blank())
cls <- cc[, .(coherent=sum(N[kind=="coherent"]), total=sum(N)), by=class]
cls[, pct := 100*coherent/total]; cls[, lab := sprintf("%s\n(n = %d)", class, total)]
cls[, lab := factor(lab, levels=cls[order(pct)]$lab)]
pB <- ggplot(cls, aes(pct, lab)) + geom_vline(xintercept=50, linetype="dashed", colour=MUTED, linewidth=0.4) +
  geom_col(fill=SEM["accent"], width=0.6, alpha=0.9) +
  geom_label(aes(label=sprintf("%.1f%%", pct)), hjust=0, nudge_x=1.5, size=TXT, colour=INK, fill="white",
             border.colour=NA, label.padding=grid::unit(1.5,"pt")) +
  scale_x_continuous(limits=c(0,105), breaks=seq(0,90,30), expand=expansion(mult=c(0,0.10))) +
  labs(x="% of resolved cores that are coherent", y=NULL, tag="b") + theme_pub(9.5) + theme(panel.grid.major.y=element_blank())
sc <- fread(file.path(REV,"results/sign_concordance_summary.csv")); sc <- sc[stratum!="ALL_sign_annotated"]
nice <- c(miRNA_target_ALL="miRNA ⊣ target, all", miRNA_target_strong="miRNA ⊣ target, strong",
          miRNA_target_weak="miRNA ⊣ target, weak (CLIP)", miRNA_target_predonly="miRNA ⊣ target, predicted only",
          TF_target_Activation="TF → target, activation", TF_target_Repression="TF → target, repression",
          TF_miRNA_Activation="TF → miRNA, activation", TF_miRNA_Repression="TF → miRNA, repression")
sc[, lab := nice[stratum]]; sc <- sc[!is.na(lab)]; sc[, sig := binom_p_vs_null < 0.05]
sc[, lab := factor(lab, levels=sc[order(concordance_rate - null_concordance_rate)]$lab)]
fmtp <- function(p){ if (p < 0.01) { e <- floor(log10(p)); m <- p/10^e; sprintf("*p* = %.1f \u00d7 10<sup>\u2212%d</sup>", m, -e) } else sprintf("*p* = %s", format(signif(p,2), scientific=FALSE)) }
sc[, plab := ifelse(sig, sapply(binom_p_vs_null, fmtp), "n.s.")]
long <- melt(sc[, .(lab, observed=100*concordance_rate, null=100*null_concordance_rate, sig)], id.vars=c("lab","sig"), variable.name="what", value.name="pct")
pC <- ggplot(long, aes(pct, lab)) + geom_line(aes(group=lab), colour=RULE, linewidth=1.6) +
  geom_point(aes(colour=what, shape=what), size=2.4) +
  ggtext::geom_richtext(data=sc, aes(x=100*pmax(concordance_rate,null_concordance_rate)+1.4, y=lab, label=plab), hjust=0, size=TXT, colour=INK, fill="white", label.colour=NA, label.padding=grid::unit(1,"pt")) +
  scale_colour_manual(values=c(observed=unname(SEM["accent"]), null=unname(SEM["null"])), labels=c("observed","expression-matched null"), name=NULL) +
  scale_shape_manual(values=c(observed=16, null=1), labels=c("observed","expression-matched null"), name=NULL) +
  scale_x_continuous(limits=c(35,82), expand=expansion(mult=c(0.02,0.16))) +
  labs(x="sign concordance (%)", y=NULL, tag="c") + theme_pub(9.5) + theme(legend.position="top", legend.justification="left", panel.grid.major.y=element_blank())
mm <- fread(file.path(REV,"results/v3/deconv_mediation_extended.csv")); mm <- mm[!is.na(ACME) & !is.na(ADE)]
mm <- mm[!(mediator %in% c("MCPcounter_Fibroblasts_recomputed","NNLS_Wu641_full_CAFs","QPROG_Wu641_full_CAFs"))]
want <- data.table(x=c("ETS1","NFKB1","SP1","RELA","hsa-miR-29a","hsa-miR-29a","hsa-miR-29b"),
                   y=c("COL1A1","COL1A1","COL1A1","COL1A1","COL3A1","COL1A1","COL3A1"),
                   path=c("ETS1>COL1A1","NFKB1>COL1A1","SP1>COL1A1","RELA>COL1A1","miR-29a|COL3A1","miR-29a|COL1A1","miR-29b|COL3A1"),
                   expr=c("italic('ETS1')~' \u2192 '~italic('COL1A1')","italic('NFKB1')~' \u2192 '~italic('COL1A1')","italic('SP1')~' \u2192 '~italic('COL1A1')","italic('RELA')~' \u2192 '~italic('COL1A1')",
                          "'miR-29a'~' \u22a3 '~italic('COL3A1')","'miR-29a'~' \u22a3 '~italic('COL1A1')","'miR-29b'~' \u22a3 '~italic('COL3A1')"))
med <- merge(want, mm, by=c("x","y"), all.x=TRUE)[, .(acme=median(ACME), lo=min(ACME), hi=max(ACME), ade=median(ADE),
      n_inc=sum(toupper(as.character(inconsistent)) %in% c("TRUE","1")), n=.N), by=path]
med <- med[match(want$path, path)]; print(med)
med[, kind := ifelse(grepl("miR", path), "miRNA arm", "TF arm")]; med[, lab := paste0(n_inc, "/", n)]
med[, path := factor(path, levels=rev(want$path))]
labmap <- setNames(want$expr, want$path)
pD <- ggplot(med, aes(y=path)) + geom_vline(xintercept=0, colour=RULE, linewidth=0.5) +
  geom_errorbarh(aes(xmin=lo, xmax=hi), height=0.16, colour=MUTED, linewidth=0.5) +
  geom_point(aes(x=acme, colour=kind), size=2.8) + geom_point(aes(x=ade), shape=4, size=2.2, colour=INK, stroke=0.7) +
  geom_text(aes(x=0.80, label=lab), size=TXT, colour=MUTED, hjust=0) +
  scale_colour_manual(values=c(`TF arm`=unname(SEM["neg"]), `miRNA arm`=unname(SEM["pos"])), name=NULL) +
  scale_y_discrete(labels=function(v) parse(text=unname(labmap[v]))) +
  scale_x_continuous(breaks=seq(-0.25,0.75,0.25), labels=label_number(accuracy=0.01, style_negative="minus"),
                     limits=c(-0.35,0.90), expand=expansion(0)) + labs(x="effect size", y=NULL, tag="d") +
  theme_pub(9.5) + theme(legend.position="top", legend.justification="left", panel.grid.major.y=element_blank())
out <- (pA | pB) / pC / pD + plot_layout(heights=c(0.95, 1, 1.2))
ggsave(file.path(FIG,"Fig4_validation_and_mediation.png"), out, width=6.85, height=8.6, dpi=600, bg="white", device=ragg::agg_png)
ggsave(file.path(FIG,"Fig4_validation_and_mediation.pdf"), out, width=6.85, height=8.6, bg="white", device=cairo_pdf)
cat("wrote\n")
