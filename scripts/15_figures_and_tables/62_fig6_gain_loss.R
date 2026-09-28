# Figure 6: behaviours gained and lost by higher-order modules relative to their embedded three-node
# cores, for the composite circuit (a) and the incoherent type-1 circuit I1 (b).
# R version of the Fig_dyn5 panel of 05_figures_springer.py (same data, same encoding: gains right of
# zero in red, losses left of zero in blue, colour intensity by module size), drawn at print width
# (174 mm) with 8 pt text in the style of the other main figures. Written 2026-09-25; it supersedes
# the Fig_dyn5 output of the Python script. Run in a UTF-8 locale (LC_ALL=en_US.UTF-8).
suppressMessages({library(data.table); library(ggplot2); library(patchwork)})
REV <- "/path/to/revision"
source(file.path(REV,"scripts/15_figures_and_tables/10_theme.R")); FIG <- file.path(REV,"figures/final")
BASE <- 9.5   # theme_pub(9.5): tick labels and legend 8 pt, axis titles 9 pt

GL <- fread(file.path(REV,"results/v3/dynamics_higher_order_gain_loss.csv"))
beh <- c("is_bistable","is_sustained_osc","is_damped_osc","is_ultrasensitive","is_memory",
         "is_noise_rejecting","is_pulse")
blab <- c("bistability", "sustained\noscillation", "damped\noscillation", "ultrasensitivity",
          "memory", "noise\nrejection", "pulse")   # thresholds are defined in the figure legend
g <- GL[metric_type == "binary" & module %in% c("n4","n5","n6") & behaviour %in% beh]
g[, y := match(behaviour, beh)]
g[, k := match(module, c("n4","n5","n6"))]          # n4 above n6 within each behaviour
w <- 0.24
g[, `:=`(ymin = y + (k - 2) * w - 0.45 * w, ymax = y + (k - 2) * w + 0.45 * w)]
bars <- rbind(g[, .(family, module, y, ymin, ymax, xmin = 0, xmax = pct_gain, dir = "gain")],
              g[, .(family, module, y, ymin, ymax, xmin = -pct_loss, xmax = 0, dir = "loss")])
bars[, key := paste(dir, module, sep = ".")]
# colours of 05_figures_springer.py (DIV_HI, DIV_LO) at alpha 0.45, 0.65 and 0.85 over white
tint <- function(hex, a) { v <- col2rgb(hex)[, 1]; rgb(t(a * v + (1 - a) * 255), maxColorValue = 255) }
cols <- c(gain.n4 = tint("#e34948", 0.45), gain.n5 = tint("#e34948", 0.65), gain.n6 = tint("#e34948", 0.85),
          loss.n4 = tint("#2a78d6", 0.45), loss.n5 = tint("#2a78d6", 0.65), loss.n6 = tint("#2a78d6", 0.85))

panel <- function(fam, title, tag, show_y) {
  ggplot(bars[family == fam]) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = key)) +
    geom_vline(xintercept = 0, colour = INK, linewidth = 0.5) +
    scale_fill_manual(values = cols, breaks = c("gain.n4","gain.n5","gain.n6"),
                      labels = c("4-node","5-node","6-node"), name = "module") +
    scale_y_reverse(breaks = seq_along(beh), labels = if (show_y) blab else NULL,
                    expand = expansion(add = 0.3)) +
    scale_x_continuous(limits = c(-21, 38), breaks = seq(-20, 30, 10),
                       labels = function(x) sub("-", "−", x)) +
    labs(x = "% of parameter sets\n← lost vs core   |   gained vs core →", y = NULL,
         title = title, tag = tag) +
    theme_pub(BASE) +
    theme(panel.grid.major.y = element_blank(), axis.title.x = element_text(lineheight = 1.1),
          plot.title = element_text(size = BASE - 0.5, colour = INK, hjust = 0.5))
}
pa <- panel("COMP_C2_toggle", "composite circuit", "a", TRUE)
pb <- panel("I1_miRNA_FFL", "I1 circuit", "b", FALSE)
fig <- (pa | pb) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
ggsave(file.path(FIG,"Fig_dyn5_higher_order_gain_loss.png"), fig, width = 6.85, height = 4.6, dpi = 600,
       bg = "white", device = ragg::agg_png)
ggsave(file.path(FIG,"Fig_dyn5_higher_order_gain_loss.pdf"), fig, width = 6.85, height = 4.6, bg = "white",
       device = cairo_pdf)
cat("wrote Fig_dyn5_higher_order_gain_loss\n")
