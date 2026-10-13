# Shared publication theme for the figures.
suppressMessages({library(ggplot2); library(grid)})
PAL   <- c(TF="#C2410C", miRNA="#1D4ED8", Gene="#047857")
EPAL  <- c(TF_target="#EA580C", TF_miRNA="#F59E0B", miRNA_target="#2563EB",
           miRNA_miRNA="#7C3AED", gene_gene="#059669")
SEM   <- c(pos="#1B6E52", neg="#A2372B", weak="#8A6212", null="#64748B", accent="#0B6E7F")
INK   <- "#0F171C"; MUTED <- "#63767F"; RULE <- "#DDE4E9"

theme_pub <- function(base=10) {
  theme_minimal(base_size=base, base_family="sans") +
    theme(
      plot.title      = element_blank(),
      plot.subtitle   = element_blank(),
      plot.tag        = element_text(size=base+2, face="bold", colour=INK),
      plot.tag.position= c(0.005, 0.995),
      plot.caption    = element_text(size=base-1.5, colour=MUTED, hjust=0, margin=margin(t=8)),
      axis.title      = element_text(size=base-0.5, colour=INK),
      axis.text       = element_text(size=base-1.5, colour=INK),
      panel.grid.minor= element_blank(),
      panel.grid.major= element_line(colour=RULE, linewidth=0.3),
      legend.title    = element_text(size=base-1, colour=INK),
      legend.text     = element_text(size=base-1.5, colour=INK),
      legend.key.size = unit(9,"pt"),
      strip.text      = element_text(size=base-0.5, face="bold", colour=INK),
      plot.margin     = margin(6,8,6,6)
    )
}
sci <- function(x, d=2) formatC(x, format="g", digits=d)
