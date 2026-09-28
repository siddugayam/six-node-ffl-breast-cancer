# Census and motif-significance figure, rebuilt from the INBOX_2026-09-26b re-runs.
# Same layout, theme and panel logic as scripts/v6/11_fig_census_motif.R, but every number is read
# from the re-run tables instead of being typed in, and the output goes to this folder only.
# usage: Rscript fig_census_motif.R <census graph variant> <null graph> <output stem>
#   census variant: table3 | nolegacy | table3_nostring | nolegacy_nostring
#   null graph:     dep | dep_nolegacy | pub | pub_nolegacy | pub_nostring | pub_nolegacy_nostring
suppressMessages({library(data.table); library(ggplot2); library(patchwork); library(scales)})
args <- commandArgs(trailingOnly=TRUE)
VAR <- args[1]; NULLG <- args[2]; STEM <- args[3]
REV <- "/path/to/revision"
HERE <- file.path(REV, "INBOX_2026-09-26b")
source(file.path(REV, "scripts/v6/10_theme.R"))

# ---- panel a: census -----------------------------------------------------------------
cc <- fread(file.path(HERE, "census/census_all_graphs.csv"))[graph == VAR]
cen <- data.table(n = cc$n, ffl = as.numeric(cc$modules),
                  sd = suppressWarnings(as.numeric(cc$modules_sd)),
                  exact = cc$method == "exhaustive", maxclass = cc$max_edge_classes)
cen[is.na(sd), sd := 0]
pA <- ggplot(cen, aes(factor(n), ffl)) +
  geom_col(aes(fill = exact), width = 0.62) +
  geom_errorbar(data = cen[!(exact)], aes(ymin = pmax(ffl - sd, 1), ymax = ffl + sd), width = 0.18,
                colour = INK, linewidth = 0.35) +
  geom_text(aes(label = label_number(scale_cut = cut_short_scale(), accuracy = 0.1)(ffl)),
            vjust = -0.6, size = 2.7, colour = INK) +
  scale_y_log10(labels = label_number(scale_cut = cut_short_scale()), expand = expansion(mult = c(0, 0.18))) +
  scale_fill_manual(values = c(`TRUE` = INK, `FALSE` = MUTED),
                    labels = c(`TRUE` = "exhaustive", `FALSE` = "RAND-ESU estimate (8 seeds)"), name = NULL) +
  labs(x = "module size (nodes)", y = expression(italic(n)*"-node FFLs  (log scale)"), tag = "a") +
  theme_pub() + theme(panel.grid.major.x = element_blank(), legend.position = "top", legend.justification = "left",
                      legend.background = element_rect(fill = "white", colour = NA))

# ---- panel b: edge classes -----------------------------------------------------------
m5 <- cen[n == 5, maxclass]; mbig <- max(cen[n >= 6, maxclass])
first_n <- min(cen[maxclass == max(cen$maxclass), n])
pB <- ggplot(cen, aes(factor(n), maxclass)) +
  geom_col(fill = SEM["accent"], width = 0.62, alpha = 0.85) +
  geom_hline(yintercept = max(cen$maxclass), linetype = "dashed", colour = SEM["neg"], linewidth = 0.4) +
  annotate("text", x = 0.6, y = max(cen$maxclass) + 0.32,
           label = sprintf("maximum: %d of 7 classes, first reached at n = %d", max(cen$maxclass), first_n),
           size = 2.4, colour = SEM["neg"], hjust = 0) +
  scale_y_continuous(limits = c(0, 7.6), breaks = 0:7, expand = expansion(mult = c(0, 0.02))) +
  labs(x = "module size (nodes)", y = "max. distinct edge classes per module", tag = "b") +
  theme_pub() + theme(panel.grid.major.x = element_blank())

# ---- panel c: motif significance -----------------------------------------------------
m <- fread(file.path(HERE, "nulls/motif_nulls_all_graphs.csv"))
m <- m[graph == NULLG & metric %in% c("Composite-FFL (once per reciprocal pair)", "TF-FFL", "miRNA-FFL",
                                      "All 3-node FFL modules", "reciprocal TF<->miRNA pairs")]
m[, null := fcase(null == "A", "A: full randomisation", null == "B", "B: reciprocity preserved",
                  null == "C", "C: B + bipartite curveball")]
m[, lab := fcase(metric == "Composite-FFL (once per reciprocal pair)", "Composite-FFL",
                 metric == "All 3-node FFL modules", "all 3-node FFLs",
                 metric == "reciprocal TF<->miRNA pairs", "reciprocal TF-miRNA pairs", default = metric)]
m[, lab := factor(lab, levels = c("reciprocal TF-miRNA pairs", "Composite-FFL", "all 3-node FFLs", "miRNA-FFL", "TF-FFL"))]
m[, dir := fifelse(fold > 1.02, "enriched", fifelse(fold < 0.98, "depleted", "not different"))]
pC <- ggplot(m, aes(x = log2(fold), y = lab, fill = dir)) +
  geom_vline(xintercept = 0, colour = RULE, linewidth = 0.5) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.2f×", fold), hjust = ifelse(log2(fold) > 0, -0.15, 1.15)), size = 2.4, colour = INK) +
  facet_wrap(~null, nrow = 1) +
  scale_fill_manual(values = c(enriched = unname(SEM["pos"]), depleted = unname(SEM["neg"]),
                               `not different` = unname(SEM["null"])), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.22, 0.22))) +
  labs(x = expression(log[2]*" fold change vs randomised networks"), y = NULL, tag = "c") +
  theme_pub() + theme(legend.position = "top", panel.grid.major.y = element_blank())

out <- (pA | pB) / pC + plot_layout(heights = c(1, 1.15))
ggsave(file.path(HERE, "figures", paste0(STEM, ".png")), out, width = 11.2, height = 8.2, dpi = 600, bg = "white")
ggsave(file.path(HERE, "figures", paste0(STEM, ".pdf")), out, width = 11.2, height = 8.2, bg = "white")
cat("wrote", STEM, "\n")
