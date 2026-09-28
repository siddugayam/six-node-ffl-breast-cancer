# Figure: FFL census across module sizes, edge-class saturation, and motif significance
# against three null models.
# Revised 2026-09-25: drawn at print width (174 mm) with 8 pt text (9 pt axis titles), counts in the
# number format of Table 3, italic n. Run in a UTF-8 locale (LC_ALL=en_US.UTF-8).
# Revised 2026-09-27: every number is read from the re-run tables of the analysis machine (INBOX_2026-09-26b,
# copied to results/v2/legacy_rerun_2026-09-26b/): the census graph and the network used for the motif tests
# exclude the 30 legacy miRNA-miRNA edges of the exemplar circuits (graphs `nolegacy` and `dep_nolegacy`);
# counts at n = 3-6 are exhaustive and n = 7 is the mean +- SD of eight RAND-ESU seeds.
suppressMessages({library(data.table); library(ggplot2); library(patchwork); library(scales)})
REV <- "/path/to/revision"
source(file.path(REV,"scripts/v6/10_theme.R"))
FIG <- file.path(REV,"figures/final"); dir.create(FIG, showWarnings=FALSE, recursive=TRUE)
BASE <- 9.5; TXT <- 8 / .pt   # theme_pub(9.5): tick labels and legends 8 pt; in-panel text 8 pt

# ---- panel A: census -----------------------------------------------------------------
RR <- file.path(REV,"results/v2/legacy_rerun_2026-09-26b")
cc <- fread(file.path(RR,"census_all_graphs.csv"))[graph=="nolegacy"]
cen <- data.table(n=cc$n, ffl=as.numeric(cc$modules), sd=suppressWarnings(as.numeric(cc$modules_sd)),
                  exact=cc$method=="exhaustive", maxclass=cc$max_edge_classes)
cen[is.na(sd), sd := 0]
sci <- function(x) { e <- floor(log10(x)); sprintf("%.1f %%*%% 10^%d", x/10^e, e) }
cen[, lab := ifelse(n==3, sprintf("'%s'", format(ffl, big.mark=",")), sapply(ffl, sci))]
pA <- ggplot(cen, aes(factor(n), ffl)) +
  geom_col(aes(fill=exact), width=0.62) +
  geom_errorbar(data=cen[!(exact)], aes(ymin=pmax(ffl-sd,1), ymax=ffl+sd), width=0.18, colour=INK, linewidth=0.35) +
  geom_text(aes(label=lab), parse=TRUE, vjust=-0.5, size=TXT, colour=INK) +
  scale_y_log10(breaks=10^c(0,3,6,9), labels=label_math(10^.x, format=log10),
                expand=expansion(mult=c(0,0.18))) +
  scale_fill_manual(values=c(`TRUE`=INK, `FALSE`=MUTED),
                    labels=c(`TRUE`="exhaustive", `FALSE`="RAND-ESU estimate (eight seeds)"), name=NULL) +
  labs(x="module size (nodes)", y=expression(italic(n)*"-node FFLs (log scale)"), tag = "a") +
  theme_pub(BASE) + theme(panel.grid.major.x=element_blank(), legend.position="top", legend.justification="left",
                      legend.background=element_rect(fill="white", colour=NA))

# ---- panel B: edge-class saturation --------------------------------------------------
pB <- ggplot(cen, aes(factor(n), maxclass)) +
  geom_col(fill=SEM["accent"], width=0.62, alpha=0.85) +
  geom_hline(yintercept=6, linetype="dashed", colour=SEM["neg"], linewidth=0.4) +
  annotate("text", x=0.6, y=6.32, label="6 of 7 classes",
           size=TXT, colour=SEM["neg"], hjust=0) +
  annotate("segment", x=3.55, xend=3.1, y=7.0, yend=6.2,
           arrow=arrow(length=unit(1.5,"mm"), type="closed"), colour=INK, linewidth=0.35) +
  annotate("text", x=4.1, y=7.3, label="first~reached~at~italic(n)==5", parse=TRUE, size=TXT, colour=INK, hjust=0.5) +
  scale_y_continuous(limits=c(0,7.6), breaks=0:7, expand=expansion(mult=c(0,0.02))) +
  labs(x="module size (nodes)", y="max. edge classes per module", tag = "b") +
  theme_pub(BASE) + theme(panel.grid.major.x=element_blank())

# ---- panel C: motif significance -----------------------------------------------------
m <- fread(file.path(RR,"motif_nulls_all_graphs.csv"))
m <- m[graph=="dep_nolegacy" & metric %in% c("Composite-FFL (once per reciprocal pair)","TF-FFL",
                                    "miRNA-FFL","All 3-node FFL modules","reciprocal TF<->miRNA pairs")]
m[, null := fcase(grepl("NULL-A", null_desc), "A: degree-preserving",
                  grepl("NULL-B", null_desc), "B: reciprocity preserved",
                  grepl("NULL-C", null_desc), "C: B + bipartite curveball")]
m[, lab := fcase(metric=="Composite-FFL (once per reciprocal pair)","Composite-FFL",
                 metric=="All 3-node FFL modules","all 3-node FFLs",
                 metric=="reciprocal TF<->miRNA pairs","reciprocal TF\u2013miRNA pairs",
                 default=metric)]
m[, lab := factor(lab, levels=c("reciprocal TF\u2013miRNA pairs","Composite-FFL","all 3-node FFLs",
                                "miRNA-FFL","TF-FFL"))]
m[, dir := fifelse(fold>1.02, "enriched", fifelse(fold<0.98, "depleted", "not different"))]
pC <- ggplot(m, aes(x=log2(fold), y=lab, fill=dir)) +
  geom_vline(xintercept=0, colour=RULE, linewidth=0.5) +
  geom_col(width=0.6) +
  geom_text(aes(label=sprintf("%.2f×", fold),
                hjust=ifelse(log2(fold)>0, -0.15, 1.15)), size=TXT, colour=INK) +
  facet_wrap(~null, nrow=1) +
  scale_fill_manual(values=c(enriched=unname(SEM["pos"]), depleted=unname(SEM["neg"]),
                             `not different`=unname(SEM["null"])), name=NULL) +
  scale_x_continuous(labels=label_number(accuracy=1, style_negative="minus"), expand=expansion(mult=c(0.36,0.36))) +
  labs(x=expression(log[2]*" fold change vs randomised networks"), y=NULL, tag = "c") +
  theme_pub(BASE) + theme(legend.position="top", panel.grid.major.y=element_blank())

if (packageVersion("patchwork") >= "1.2.0") pA <- free(pA)   # keep the y-axis title of panel a next to its axis
out <- (pA | pB) / pC + plot_layout(heights=c(1, 1.15))
ggsave(file.path(FIG,"Fig3_census_and_motif_significance.png"), out, width=6.85, height=5.9, dpi = 600, bg="white", device=ragg::agg_png)
ggsave(file.path(FIG,"Fig3_census_and_motif_significance.pdf"), out, width=6.85, height=5.9, bg="white", device=cairo_pdf)
cat("wrote Fig3_census_and_motif_significance\n")
