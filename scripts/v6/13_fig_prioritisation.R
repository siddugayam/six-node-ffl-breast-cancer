# Figure: the 30 prioritised nodes, their evidence-domain profiles, and the divergence between
# ranking by FFL topology and ranking by evidence.
# Revised 2026-09-25: drawn at print width (174 mm) with 8 pt text (9 pt axis titles); gene symbols
# italic, miRNA names roman. Run in a UTF-8 locale (LC_ALL=en_US.UTF-8).
suppressMessages({library(data.table); library(ggplot2); library(patchwork); library(ggrepel)})
REV <- "/path/to/revision"
source(file.path(REV,"scripts/v6/10_theme.R")); FIG <- file.path(REV,"figures/final")
BASE <- 9.5; TXT <- 8 / .pt   # theme_pub(9.5): tick labels and legends 8 pt; in-panel text 8 pt
P <- fread(file.path(REV,"results/v5/node_prioritisation_full.csv"))
P <- P[!is.na(priority)]

# ---- A: rank divergence, TFs -----------------------------------------------------------
D <- P[type=="TF"][order(-ffl_cores)][, rank_ffl := seq_len(.N)]
D <- D[order(-priority)][, rank_comp := seq_len(.N)]
D[, gap := rank_comp - rank_ffl]
lab <- D[rank_ffl<=10 | rank_comp<=10]
pA <- ggplot(D, aes(rank_ffl, rank_comp)) +
  geom_abline(slope=1, intercept=0, colour=RULE, linewidth=0.5) +
  geom_point(aes(size=ffl_cores, fill=abs(gap)>40), shape=21, colour="grey30", stroke=0.25, alpha=0.9) +
  geom_text_repel(data=lab, aes(label=name), size=TXT, fontface="italic", max.overlaps=30, seed=42,
                  segment.size=0.15, segment.colour=MUTED, min.segment.length=0.1) +
  scale_fill_manual(values=c(`TRUE`=unname(SEM["neg"]), `FALSE`=unname(SEM["null"])),
                    labels=c("agree","diverge > 40 ranks"), name=NULL) +
  scale_size_continuous(range=c(0.6,4.5), name="FFL cores") +
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  labs(x="rank by FFL participation (1 = most central)",
       y="rank by seven-domain evidence composite", tag = "a") +
  theme_pub(BASE) + theme(legend.position="right")

# ---- B: domain profile heatmap of the 30 ----------------------------------------------
dom <- c(T_topology="topology", E_evidence="edge evidence", D_expression="expression",
         R_replication="replication", M_multiomic="multi-omic", C_clinical="clinical",
         F_functional="essentiality")
top <- P[, head(.SD, 10), by=type]
H <- melt(top[, c("name","type",names(dom)), with=FALSE], id.vars=c("name","type"),
          variable.name="domain", value.name="score")
H[, domain := factor(dom[as.character(domain)], levels=unname(dom))]
ord <- top[order(type, -priority), name]
H[, name := factor(name, levels=rev(ord))]
H[, type := factor(type, levels=c("TF","Gene","miRNA"))]
ylab <- function(v) parse(text=ifelse(grepl("^hsa-", v), paste0("'", v, "'"), paste0("italic('", v, "')")))
pB <- ggplot(H, aes(domain, name, fill=score)) +
  geom_tile(colour="white", linewidth=0.4) +
  facet_grid(type ~ ., scales="free_y", space="free_y", switch="y") +
  scale_fill_gradient2(low="#F2F6F8", mid="#7FB3C4", high="#0B4A5C", midpoint=0.5,
                       na.value="#EEF1F3", limits=c(0,1), name="rank-normalised\nevidence") +
  scale_y_discrete(labels=ylab) +
  labs(x=NULL, y=NULL, tag = "b") +
  theme_pub(BASE) +
  theme(axis.text.x=element_text(angle=38, hjust=1), panel.grid=element_blank(),
        strip.placement="outside", strip.text.y.left=element_text(angle=0),
        legend.position="right", legend.key.height=unit(24,"pt"))

out <- pA / pB + plot_layout(heights=c(1, 1.5))
ggsave(file.path(FIG,"Fig5_prioritisation.png"), out, width=6.85, height=8.6, dpi = 600, bg="white", device=ragg::agg_png)
ggsave(file.path(FIG,"Fig5_prioritisation.pdf"), out, width=6.85, height=8.6, bg="white", device=cairo_pdf)
cat("wrote Fig5_prioritisation\n")
