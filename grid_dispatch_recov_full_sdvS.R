rm(list=ls()) 
source("dmc/dmc.R")
load_model("LBA", "lbaN_B.R")

# Load
load("samples/recov_sTPPM_full_sdvS.RData")

recov_samples <- h.run.dmc(recov_samples,
                     cores = 24,
                     report = 10,
                     p.migrate = 0.05)

save(recov_samples, file = "samples/recov_sTPPM_full_sdvS.RData")

# # Check for convergence
# gelman.diag.dmc(recov_samples)
# check.recovery.dmc(recov_samples)

recov_samples1 <- h.RUN.dmc(h.samples.dmc(samples = recov_samples,
                                    nmc = 100,
                                    thin = 10,
                                    n.chains = 159),
                      cores = 24,
                      report = 10, 
                      verbose = TRUE)

# Save samples
save(recov_samples, recov_samples1, file = "samples/recov_sTPPM_full_sdvS.RData")

# recov_samples2 <- h.RUN.dmc(h.samples.dmc(samples = recov_samples1,
#                                     nmc = 100,
#                                     thin = 10,
#                                     n.chains = 159),
#                       cores = 16,
#                       report = 10)
# # Save samples
# save(recov_samples, recov_samples1, recov_samples2, file = "samples/recov_sTPPM_full_sdvS.RData")


# Sampling plots ----------------------------------------------------------

# Load samples
load("samples/recov_sTPPM_full_sdvS.RData")

# samples1 <- samples2

# Chains
pdf("plots/chains_recov_full_sdvS.pdf", height = 6, width = 8)
par(mfrow = c(3, 3))
for(i in 1:length(recov_samples1)) {
  plot.dmc(recov_samples1, pll.chain = TRUE, subject = i, start = 1)
  title(paste0('Subject ', names(recov_samples1)[i]))
}
dev.off()

# Chains by subject
pdf("plots/chains_by_sub_recov_full_sdvS.pdf", height = 6, width = 8)
for(sub in names(recov_samples1)) {
  par(mfrow = c(1, 1))
  plot.dmc(recov_samples1, hyper = FALSE, subject = sub, pll.chain = TRUE)
  mtext(paste0('Participant ', sub, ' posterior log likelihoods of all chains'),
        outer = TRUE, line = -1.5)
  plot.dmc(recov_samples1, hyper = FALSE, subject = sub, layout = c(2, 2), density = TRUE)
  mtext(paste0('Participant ', sub), outer = TRUE, line = -1.5)
}
dev.off()


# Posterior predictive fits -----------------------------------------------

# Sample posterior predictives
post_preds <- h.post.predict.dmc(recov_samples1)
save(post_preds, file = "deriv/post_preds_recov_full_sdvS.RData")

# Load posterior predictives
load("deriv/post_preds_recov_full_sdvS.RData")

# Plot fits to CDFs
pdf("plots/cdfs_recov_full_sdvS.pdf", height = 6, width = 8)
plot.pp.dmc(post_preds, "cdf", layout = c(2, 2), model.legend = FALSE)
dev.off()

# Subject-level fits
pdf("plots/cdfs_by_sub_recov_full_sdvS.pdf", height = 6, width = 8)
for (i in seq_len(length(post_preds))) {
  plot.pp.dmc(post_preds[[i]], "cdf", layout = c(2, 2), model.legend = FALSE)
}
dev.off()


# Fit indices
h.IC.dmc(recov_samples1, DIC = TRUE)

# Save posterior predictive sims
post_pred_sims <- h.post.predict.dmc(recov_samples1, save.simulation = TRUE)
save(post_pred_sims, file = "deriv/post_pred_sims_recov_full_sdvS.RData")



# Parameters --------------------------------------------------------------

# Summary
parms <- summary.dmc(recov_samples1)
parms <- do.call(rbind, lapply(parms, function(x) x$statistics[, 1]))
# head(parms)
save(parms, file = "deriv/map_parms_recov_full_sdvS.RData")

# Recovery
h.check.recovery.dmc(recov_samples1, digits = 2)

# Effective sample size
ess <- effectiveSize.dmc(recov_samples1)
do.call(rbind, ess)
min(do.call(rbind, ess))

# # Correlations
# pdf("plots/pairs_recov_full_sdvS.pdf", height = 6, width = 8)
# pairs.dmc(recov_samples1, start = 50)
# dev.off()