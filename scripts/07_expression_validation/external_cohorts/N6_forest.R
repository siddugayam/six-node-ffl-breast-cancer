## N6_forest.R -- forest plots from newcohorts_meta.csv (base graphics, no deps)
setwd("/path/to/revision")
OUT <- "results/v2"
M <- read.csv(file.path(OUT,"newcohorts_meta.csv"), stringsAsFactors=FALSE)
panels <- list(
  c("COL1A1~COL3A1","raw"),  c("COL1A1~COL3A1","cafA"),
  c("ETS1->COL1A1","raw"),   c("ETS1->COL1A1","cafA"),
  c("NFKB1->COL1A1","raw"),  c("NFKB1->COL1A1","cafA"))
draw <- function(ed, aj){
  d <- M[M$edge==ed & M$adjustment==aj & M$pool=="tumour_only", ]
  if(!nrow(d)) return(invisible(NULL))
  d <- rbind(d[d$row_type=="study", ], d[d$row_type=="RE_pooled", ])
  k <- nrow(d); y <- rev(seq_len(k))
  xr <- range(c(d$ci_lo_rho, d$ci_hi_rho, 0), na.rm=TRUE)
  xr <- xr + c(-0.05, 0.05)*diff(xr)
  par(mar=c(6.5,8.5,3,4.5), mgp=c(2.2,0.7,0))
  plot(NA, xlim=xr, ylim=c(0.4, k+0.6), yaxt="n", ylab="", bty="n",
       xlab="Spearman rho (95% CI)",
       main=paste0(ed, "  [", ifelse(aj=="raw","unadjusted","CAF-adjusted"), "]"),
       cex.main=1)
  abline(v=0, col="grey60", lty=2)
  for(i in seq_len(k)){
    pooled <- d$row_type[i]=="RE_pooled"
    segments(d$ci_lo_rho[i], y[i], d$ci_hi_rho[i], y[i],
             lwd=ifelse(pooled,2.5,1.5), col=ifelse(pooled,"firebrick","black"))
    points(d$rho[i], y[i], pch=ifelse(pooled,18,15),
           cex=ifelse(pooled, 2.0, 0.7 + 1.4*sqrt(d$weight_pct[i]/100)),
           col=ifelse(pooled,"firebrick","black"))
  }
  axis(2, at=y, labels=paste0(sub("POOLED .*","POOLED (RE)", d$cohort), " (n=", d$n, ")"), las=1, tick=FALSE, cex.axis=0.75)
  lab <- sprintf("%.2f", d$rho)
  axis(4, at=y, labels=lab, las=1, tick=FALSE, cex.axis=0.75)
  pr <- d[d$row_type=="RE_pooled", ]
  mtext(sprintf("RE pooled rho=%.3f [%.3f, %.3f], p=%.3g, I2=%.1f%% (k=%d, tau2=%.3f)",
                pr$rho, pr$ci_lo_rho, pr$ci_hi_rho, pr$p_pooled, pr$I2, pr$k, pr$tau2),
        side=1, line=4.6, cex=0.62)
}
pdf(file.path(OUT,"newcohorts_forest.pdf"), width=8.5, height=12)
par(mfrow=c(3,2))
for(p in panels) draw(p[1], p[2])
dev.off()
png(file.path(OUT,"newcohorts_forest.png"), width=1700, height=2400, res=150)
par(mfrow=c(3,2))
for(p in panels) draw(p[1], p[2])
dev.off()
cat("wrote", file.path(OUT,"newcohorts_forest.pdf"), "and .png\n")
