#!/usr/bin/env Rscript
## jam_04_permutation_control.R
## (a) the component structure of the final annealed state (why ranks 2-5 of a
##     single run are not modules), and
## (b) the control the corrected score does not provide: re-run the whole search
##     on node scores permuted over the network.  The permutation destroys every
##     relation between differential expression and network position while
##     keeping the z distribution, so the best score it reaches is the score a
##     module carrying no network information would obtain.
source("/path/to/revision/scripts/14_module_detection/jam_common.R")
sourceCpp(file.path(SCR, "jam_engine.cpp"))
R <- readRDS(file.path(BASE, "cache/v7/jam/jam_runs.rds"))
g <- R$g; z <- R$scores$primary$z; cal <- R$cals$sets; nodes <- V(g)$name

r <- jam_run(g, z, cal, seed = 1, iters = 2e5)
fs <- data.table(rank = seq_along(r$final_modules), size = lengths(r$final_modules),
                 score = unname(r$final_scores))
fwrite(fs, file.path(OUT, "jam_final_state_components.csv"))
jam_log("final state: ", nrow(fs), " components, sizes ",
        paste(head(fs$size, 8), collapse = ","), " ; scores ",
        paste(round(head(fs$score, 8), 2), collapse = ","))

NPERM <- 20
perm <- rbindlist(lapply(seq_len(NPERM), function(b) {
  set.seed(9000 + b)
  zp <- z[sample(length(z))]; names(zp) <- nodes
  rr <- jam_run(g, zp, cal, seed = b, iters = 2e5)
  data.table(perm = b, size = length(rr$modules[[1]]), score = unname(rr$scores[1]))
}))
obs <- max(R$modules[variant == "primary" & rank == 1]$score)
perm[, observed := obs]
perm[, p_emp := (sum(score >= obs) + 1) / (.N + 1)]
fwrite(perm, file.path(OUT, "jam_permutation_control.csv"))
jam_log("observed best corrected score ", round(obs, 2),
        " ; node-label-permuted searches ", round(min(perm$score), 2), "-",
        round(max(perm$score), 2), " (median ", round(median(perm$score), 2),
        ", sizes ", min(perm$size), "-", max(perm$size), ") ; empirical p = ",
        signif(perm$p_emp[1], 3))
