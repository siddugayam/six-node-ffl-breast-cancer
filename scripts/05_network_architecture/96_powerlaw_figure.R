suppressPackageStartupMessages(library(poweRlaw)); set.seed(20260908)
RES <- "/path/to/revision/results/v3"
FIG <- "/path/to/revision/figures/v3"
x <- read.csv(file.path(RES,"systems_degree_table.csv"))$degree_total; x <- x[x>0]
draw <- function() {
  m33 <- displ$new(x); m33$setXmin(33); m33$setPars(estimate_pars(m33))
  m24 <- displ$new(x); m24$setXmin(24); m24$setPars(estimate_pars(m24))
  ln  <- dislnorm$new(x); ln$setXmin(33); ln$setPars(estimate_pars(ln))
  plot(m33, pch=16, cex=.55, col="grey35", xlab="total degree k", ylab="P(K >= k)",
       main="Total degree is not scale-free")
  lines(m33, col="#c0392b", lwd=2.5); lines(m24, col="#e67e22", lwd=2.5, lty=3)
  lines(ln,  col="#2471a3", lwd=2.5, lty=2)
  legend("bottomleft", bty="n", cex=.8, lwd=2.5, lty=c(1,3,2),
         col=c("#c0392b","#e67e22","#2471a3"),
         legend=c("power law, R/poweRlaw KS-optimal xmin=33 (a=3.46, CSN p=0.065)",
                  "power law, python powerlaw xmin=24 (a=2.94, CSN p=0.008)",
                  "lognormal at xmin=33 (LR vs power law -1.62, p=0.106)"))
  mtext("CSN bootstrap p <= 0.1 at BOTH x_min choices: the power law is rejected either way",
        side=3, line=0.1, cex=.72, col="grey25")
}
pdf(file.path(FIG,"Fig_v3_A2_powerlaw_xmin_sensitivity.pdf"), width=6.5, height=5.4); draw(); dev.off()
png(file.path(FIG,"Fig_v3_A2_powerlaw_xmin_sensitivity.png"), width=1300, height=1080, res=190); draw(); dev.off()
cat("figure written\n")
