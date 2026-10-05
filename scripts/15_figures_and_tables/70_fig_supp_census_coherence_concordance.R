# Supplementary Figs S2, S3 and S6 of the research-article version (1 Oct 2026): panels that moved out of the main
# figures. S2 = census and edge-class saturation (formerly Fig. 1a, b); S3 = coherence typology (formerly Fig. 2a, b);
# S6 = sign concordance by edge class and evidence tier (formerly Fig. 2c). Code and data as in 11_fig_census_motif.R
# and 12_fig_validation.R; drawn at print width (174 mm), 8 pt text, 600 dpi. Run in a UTF-8 locale (LC_ALL=en_US.UTF-8).
suppressMessages({library(data.table); library(ggplot2); library(patchwork); library(scales); library(ragg); library(ggtext)})
REV <- Sys.getenv("FFL_REV", unset = "/path/to/revision")
source(file.path(REV,"scripts/15_figures_and_tables/10_theme.R"))
FIG <- file.path(REV,"figures/final"); dir.create(FIG, showWarnings=FALSE, recursive=TRUE)
BASE <- 9.5; TXT <- 8 / .pt

# ---- S2: census (a) and edge-class saturation (b) --------------------------------------
RR <- file.path(REV,"results/v2/census_rerun")
cc <- fread(file.path(RR,"census_all_graphs.csv"))[graph=="nolegacy"]
cen <- data.table(n=cc$n, ffl=as.numeric(cc$modules), sd=suppressWarnings(as.numeric(cc$modules_sd)),
                  exact=cc$method=="exhaustive", maxclass=cc$max_edge_classes)
cen[is.na(sd), sd := 0]
sci1 <- function(x) { e <- floor(log10(x)); sprintf("%.1f %%*%% 10^%d", x/10^e, e) }
cen[, lab := ifelse(n==3, sprintf("'%s'", format(ffl, big.mark=",")), sapply(ffl, sci1))]
pA <- ggplot(cen, aes(factor(n), ffl)) +
  geom_col(aes(fill=exact), width=0.62) +
  geom_errorbar(data=cen[!(exact)], aes(ymin=pmax(ffl-sd,1), ymax=ffl+sd), width=0.18, colour=INK, linewidth=0.35) +
  geom_text(aes(label=lab), parse=TRUE, vjust=-0.5, size=TXT, colour=INK) +
  scale_y_log10(breaks=10^c(0,3,6,9), labels=label_math(10^.x, format=log10), expand=expansion(mult=c(0,0.18))) +
  scale_fill_manual(values=c(`TRUE`=INK, `FALSE`=MUTED),
                    labels=c(`TRUE`="exhaustive", `FALSE`="RAND-ESU estimate (eight seeds)"), name=NULL) +
  labs(x="module size (nodes)", y=expression(italic(n)*"-node FFLs (log scale)"), tag="a") +
  theme_pub(BASE) + theme(panel.grid.major.x=element_blank(), legend.position="top", legend.justification="left",
                          legend.background=element_rect(fill="white", colour=NA))
pB <- ggplot(cen, aes(factor(n), maxclass)) +
  geom_col(fill=SEM["accent"], width=0.62, alpha=0.85) +
  geom_hline(yintercept=6, linetype="dashed", colour=SEM["neg"], linewidth=0.4) +
  annotate("text", x=0.6, y=6.32, label="6 of 7 classes", size=TXT, colour=SEM["neg"], hjust=0) +
  annotate("segment", x=3.55, xend=3.1, y=7.0, yend=6.2, arrow=arrow(length=unit(1.5,"mm"), type="closed"), colour=INK, linewidth=0.35) +
  annotate("text", x=4.1, y=7.3, label="first~reached~at~italic(n)==5", parse=TRUE, size=TXT, colour=INK, hjust=0.5) +
  scale_y_continuous(limits=c(0,7.6), breaks=0:7, expand=expansion(mult=c(0,0.02))) +
  labs(x="module size (nodes)", y="max. edge classes per module", tag="b") +
  theme_pub(BASE) + theme(panel.grid.major.x=element_blank())
if (packageVersion("patchwork") >= "1.2.0") pA <- free(pA)
S2 <- pA | pB
ggsave(file.path(FIG,"FigS2_census.png"), S2, width=6.85, height=3.2, dpi=600, bg="white", device=ragg::agg_png)
ggsave(file.path(FIG,"FigS2_census.pdf"), S2, width=6.85, height=3.2, bg="white", device=cairo_pdf)

# ---- S3: coherence typology (a) and coherent share by class (b) -------------------------
co <- fread(file.path(REV,"results/v2/ffl_cores_coherence_corrected.csv"))
cc3 <- co[coherence!="unresolved", .N, by=.(class, coherence)]
cc3[, kind := ifelse(substr(coherence,1,1)=="C","coherent","incoherent")]
tot <- cc3[, .(N=sum(N)), by=.(coherence, kind)][order(-N)]
tot[, coherence := factor(coherence, levels=rev(tot$coherence))]
qA <- ggplot(tot, aes(N, coherence, fill=kind)) + geom_col(width=0.66) +
  geom_text(aes(label=N), hjust=-0.18, size=TXT, colour=INK) +
  scale_fill_manual(values=c(coherent=unname(SEM["pos"]), incoherent=unname(SEM["neg"])), name=NULL) +
  scale_x_continuous(expand=expansion(mult=c(0,0.16))) + labs(x="3-node cores", y=NULL, tag="a") +
  theme_pub(BASE) + theme(legend.position="top", legend.justification="left", panel.grid.major.y=element_blank())
cls <- cc3[, .(coherent=sum(N[kind=="coherent"]), total=sum(N)), by=class]
cls[, pct := 100*coherent/total]; cls[, lab := sprintf("%s\n(n = %d)", class, total)]
cls[, lab := factor(lab, levels=cls[order(pct)]$lab)]
qB <- ggplot(cls, aes(pct, lab)) + geom_vline(xintercept=50, linetype="dashed", colour=MUTED, linewidth=0.4) +
  geom_col(fill=SEM["accent"], width=0.6, alpha=0.9) +
  geom_label(aes(label=sprintf("%.1f%%", pct)), hjust=0, nudge_x=1.5, size=TXT, colour=INK, fill="white",
             border.colour=NA, label.padding=grid::unit(1.5,"pt")) +
  scale_x_continuous(limits=c(0,105), breaks=seq(0,90,30), expand=expansion(mult=c(0,0.10))) +
  labs(x="% of resolved cores that are coherent", y=NULL, tag="b") + theme_pub(BASE) + theme(panel.grid.major.y=element_blank())
S3 <- qA | qB
ggsave(file.path(FIG,"FigS3_coherence.png"), S3, width=6.85, height=3.0, dpi=600, bg="white", device=ragg::agg_png)
ggsave(file.path(FIG,"FigS3_coherence.pdf"), S3, width=6.85, height=3.0, bg="white", device=cairo_pdf)

# ---- S6: sign concordance of predicted edges ------------------------------------------
sc <- fread(file.path(REV,"results/sign_concordance_summary.csv")); sc <- sc[stratum!="ALL_sign_annotated"]
nice <- c(miRNA_target_ALL="miRNA ⊣ target, all", miRNA_target_strong="miRNA ⊣ target, strong",
          miRNA_target_weak="miRNA ⊣ target, weak (CLIP)", miRNA_target_predonly="miRNA ⊣ target, predicted only",
          TF_target_Activation="TF → target, activation", TF_target_Repression="TF → target, repression",
          TF_miRNA_Activation="TF → miRNA, activation", TF_miRNA_Repression="TF → miRNA, repression")
sc[, lab := nice[stratum]]; sc <- sc[!is.na(lab)]; sc[, sig := binom_p_vs_null < 0.05]
sc[, lab := factor(lab, levels=sc[order(concordance_rate - null_concordance_rate)]$lab)]
fmtp <- function(p){ if (p < 0.01) { e <- floor(log10(p)); m <- p/10^e; sprintf("*p* = %.1f × 10<sup>−%d</sup>", m, -e) } else sprintf("*p* = %s", format(signif(p,2), scientific=FALSE)) }
sc[, plab := ifelse(sig, sapply(binom_p_vs_null, fmtp), "n.s.")]
long <- melt(sc[, .(lab, observed=100*concordance_rate, null=100*null_concordance_rate, sig)], id.vars=c("lab","sig"), variable.name="what", value.name="pct")
rC <- ggplot(long, aes(pct, lab)) + geom_line(aes(group=lab), colour=RULE, linewidth=1.6) +
  geom_point(aes(colour=what, shape=what), size=2.4) +
  ggtext::geom_richtext(data=sc, aes(x=100*pmax(concordance_rate,null_concordance_rate)+1.4, y=lab, label=plab), hjust=0, size=TXT, colour=INK, fill="white", label.colour=NA, label.padding=grid::unit(1,"pt")) +
  scale_colour_manual(values=c(observed=unname(SEM["accent"]), null=unname(SEM["null"])), labels=c("observed","expression-matched null"), name=NULL) +
  scale_shape_manual(values=c(observed=16, null=1), labels=c("observed","expression-matched null"), name=NULL) +
  scale_x_continuous(limits=c(35,82), expand=expansion(mult=c(0.02,0.16))) +
  labs(x="sign concordance (%)", y=NULL) + theme_pub(BASE) +
  theme(legend.position="top", legend.justification="left", panel.grid.major.y=element_blank())
ggsave(file.path(FIG,"FigS6_sign_concordance.png"), rC, width=6.85, height=3.5, dpi=600, bg="white", device=ragg::agg_png)
ggsave(file.path(FIG,"FigS6_sign_concordance.pdf"), rC, width=6.85, height=3.5, bg="white", device=cairo_pdf)
cat("wrote FigS2_census, FigS3_coherence, FigS6_sign_concordance\n")
