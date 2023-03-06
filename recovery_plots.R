# Clear workspace
rm(list = ls())

# Set working directory to top-level folder containing DMC
getwd()

# Load libraries and DMC functions
source("dmc/dmc.R")
load_model("LBA", "lbaN_B.R")
# source("dmc/dmc_ATC.R")
# source("dmc/dmc_PMDC.R")


# -------------------------------------------------------------------------

# Functions ---------------------------------------------------------------
# Get data frame for recovery plots
get.ggdf.recov <- function(post_summaries, top_msds, grepchr = "B", n) {
  tmp <- list()
  j = 0
  for (i in grep(grepchr, colnames(post_summaries[[1]]))) {
    j = j + 1
    tmp[[j]] <- data.frame(cbind(post_summaries[[1]][,i],
                                 post_summaries[[2]][,i], 
                                 post_summaries[[3]][,i]))
    name <- rownames(top_msds)[i]
    colnames(tmp[[j]]) <- c("M", "LCI", "HCI")
    tmp[[j]]$reps <- 1:n
    tmp[[j]]$param <- name
    tmp[[j]]$top <- top_msds$M[i]
  }
  do.call(rbind, tmp)
}

# Make recovery plot
ggplot.recov <- function(ggdf, ncol = 2) {
  ggplot(ggdf, aes(reps, M)) +
    geom_ribbon(aes(ymin = LCI, ymax = HCI), 
                col = "lightblue",
                fill = "lightblue", 
                alpha = 0.8) + 
    geom_point(size = 0.5, 
               alpha = 0.8) +
    geom_hline(aes(yintercept = top), 
               linetype = 1, 
               size = 1,
               alpha = 0.5) + 
    ylab("") +
    facet_wrap(~ param, ncol = ncol)
}

# Get data frame for mimicry plots
get.ggdf.mimic <- function(post_summaries, top_msds, recov_msds, sim_pop_mu, grepchr = "B", n) {
  tmp <- list()
  j = 0
  for (i in grep(grepchr, colnames(post_summaries[[1]]))){
    j = j + 1
    tmp[[j]] <- data.frame(cbind(post_summaries[[1]][,i],
                                 post_summaries[[2]][,i],
                                 post_summaries[[3]][,i]))
    name <- rownames(recov_msds)[i]
    colnames(tmp[[j]]) <- c("M", "LCI", "HCI")
    tmp[[j]]$reps <- 1:n
    tmp[[j]]$param <- name
    tmp[[j]]$top <- top_msds$M[i]
    tmp[[j]]$recov <- recov_msds$M[i]
    tmp[[j]]$true <- as.numeric(sim_pop_mu[i])
  }
  do.call(rbind, tmp)
}

# Make mimicry plot (plot recovered values against true values and top model)
ggplot.mimic <- function(ggdf, ncol = 2) {
  ggplot(ggdf, aes(reps, M)) +
    geom_ribbon(aes(ymin = LCI, ymax = HCI), 
                col = "lightblue",
                fill = "lightblue",
                alpha = 0.8) + ylab("") +
    geom_point(size = 0.5,
               alpha = 0.8) + 
    geom_hline(aes(yintercept = top, linetype = "Top"), 
               color = "blue", size = 1, alpha = 0.5) +
    geom_hline(aes(yintercept = recov, linetype = "Recovered"), 
               color = "red", size = 1, alpha = 0.5) +
    geom_hline(aes(yintercept = true, linetype = "True"), 
               color = "black", size = 1, alpha = 0.5, show.legend = TRUE) +
    scale_linetype_manual(name = "Mean", values = c(1, 1, 1),
                          guide = guide_legend(override.aes = list(color = c("red", "blue", "black")))) +
    facet_wrap( ~ param, ncol = ncol)
}


# -------------------------------------------------------------------------

# Load top samples object
print(load("samples/sTPPM_full_sdvS.RData"))
top_samples <- samples1
colnames(top_samples[[1]]$theta)

# Check your samples are all the same length for every participant
for(i in 1:length(top_samples)) {
  print(top_samples[[i]]$nmc)
}

# Bring longer thetas down to min nmc by sampling
nmcs <- sapply(top_samples, function(x) x$nmc)
nmc <- min(nmcs)
for (i in 1:length(top_samples)) if (nmcs[i] > nmc) top_samples[[i]]$theta <-
  top_samples[[i]]$theta[,,sample(1:dim(top_samples[[i]]$theta)[3], nmc)]
top_samps <- lapply(top_samples, function(x) x["theta"])

# # Get summary
# top_post_summ <- get.participant.median.CIs(top_samps)
# top_post_summ
# 
# # Save
# save(top_post_summ, file = "deriv/top_post_summ_full_sdvS.RData")
print(load("deriv/top_post_summ_full_sdvS.RData"))

# Get the posterior mean values that we simulated from (see recovery_setup.R)
top_msds <- data.frame(get.msds(top_samples))
top_msds
sim_pop_mu <- top_msds$M
names(sim_pop_mu) <- rownames(top_msds)
sim_pop_mu


# Load recovery samples
print(load("samples/recov_sTPPM_full_sdvS.RData"))
recov_samples <- recov_samples1

# Check your samples are all the same length for every participant
for(i in 1:length(recov_samples)) {
  print(recov_samples[[i]]$nmc)
}

# Bring longer thetas down to min nmc by sampling
nmcs <- sapply(recov_samples, function(x) x$nmc)
nmc <- min(nmcs)
for (i in 1:length(recov_samples)) if (nmcs[i] > nmc) recov_samples[[i]]$theta <-
  recov_samples[[i]]$theta[,,sample(1:dim(recov_samples[[i]]$theta)[3], nmc)]
recov_samps <- lapply(recov_samples, function(x) x["theta"])

# # Get summary
# recov_post_summ <- get.participant.median.CIs(recov_samps)
# recov_post_summ
# 
# # Save
# save(recov_post_summ, file = "deriv/recov_post_summ_full_sdvS.RData")
print(load("deriv/recov_post_summ_full_sdvS.RData"))




# -------------------------------------------------------------------------


# Parameter correlations --------------------------------------------------

# Plot correlation between true and recovered parameters

# Get true pars from top model
top_msds <- data.frame(get.msds(top_samples))
top_msds
top_pars <- top_msds$M; names(top_pars) <- rownames(top_msds)
top_pars
top_pars == sim_pop_mu

# Get recovered pars 
recov_msds <- data.frame(get.msds(recov_samples))
recov_msds
recov_pars <- recov_msds$M; names(recov_pars) <- rownames(recov_msds)
recov_pars


pars <- data.frame(cbind(top_pars, recov_pars))
pars$type <- NA
pars$type[ grep("^A", rownames(pars)) ] <- "A"
pars$type[ grep("^B", rownames(pars)) ] <- "B"
pars$type[ grep("^t0", rownames(pars)) ] <- "t0"
pars$type[ grep("^mean_v", rownames(pars)) ] <- "mean_v"
pars$type[ grep("^sd_v", rownames(pars)) ] <- "sd_v"
pars

pars_correlation <- cor.test(pars$top_pars, pars$recov_pars)
pars_correlation$estimate

# Plot true vs. recovered values
ggplot(data = pars, mapping = aes(x = top_pars, y = recov_pars)) +
  geom_point(aes(col = type), size = 2) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  labs(title = "True vs. recovered parameters",
       x = "Theta true", 
       y = "Theta recovered", 
       col = "Par. type") +
  theme_minimal()


# -------------------------------------------------------------------------

# Check recovery ----------------------------------------------------------

# # Top model
# sim_pop_mu <- top_pars
# check_recov <- h.check.recovery.dmc(ps = sim_pop_mu, samples = top_samples)
# check_recov
# save(check_recov, file = "deriv/check_recov_full_sdvS.RData")
# print(load("deriv/check_recov_full_sdvS.RData"))
# 
# # Top recovery model
# check_recov <- h.check.recovery.dmc(ps = sim_pop_mu, samples = recov_samples)
# check_recov
# save(check_recov, file = "deriv/check_recov_recov_full_sdvS.RData")
# print(load("deriv/check_recov_recov_full_sdvS.RData"))



# -------------------------------------------------------------------------

# Plot recovery against top means -----------------------------------------

# Get data frames for ggplot
recov_A <- get.ggdf.recov(recov_post_summ, top_msds, grepchr = "^A", n = 100)
recov_B <- get.ggdf.recov(recov_post_summ, top_msds, grepchr = "^B", n = 100)
recov_t0 <- get.ggdf.recov(recov_post_summ, top_msds, grepchr = "^t0", n = 100)
recov_mean_v <- get.ggdf.recov(recov_post_summ, top_msds, grepchr = "^mean_v", n = 100)
recov_sd_v <- get.ggdf.recov(recov_post_summ, top_msds, grepchr = "^sd_v", n = 100)


# Plot recovery
A_plot <- ggplot.recov(recov_A, ncol = 2) + xlab("") + ggtitle("A") +
  theme_classic()

B_plot <- ggplot.recov(recov_B, ncol = 4) + xlab("") + ggtitle("Threshold") +
  theme_classic()
ggsave("plots/recovery_recov_full_sdvS_B.png", plot = B_plot, 
       width = 2000, height = 1400, units = "px")

t0_plot <- ggplot.recov(recov_t0, ncol = 1) + xlab("") + ggtitle("Non-decision time") +
  theme_classic()

mean_v_plot <- ggplot.recov(recov_mean_v, ncol = 5) + xlab("") + ggtitle("Accumulation rate") +
  theme_classic()
ggsave("plots/recovery_recov_full_sdvS_mean_v.png", plot = mean_v_plot, 
       width = 2000, height = 2600, units = "px")

sd_v_plot <- ggplot.recov(recov_sd_v, ncol = 2) + xlab("") + ggtitle("sd_v") +
  theme_classic()


library("grid")
pdf("plots/recovery_recov_full_sdvS.pdf", height = 4, width = 8)
grid.arrange(t0_plot, A_plot, sd_v_plot,
             ncol = 3, nrow = 1,
             top = textGrob("", gp = gpar(fontsize = 15, font = 8)),
             layout_matrix = rbind(c(3, 3, 4, 4, 5, 5, 5, 5)))
dev.off()


# -------------------------------------------------------------------------

# Mimicry -----------------------------------------------------------------

# If these all recover then we are in good shape to start taking 
# the effects in the top model seriously

# No threshold effects ----------------------------------------------------

# Load samples object from reduced model
print(load("samples/recov_sTPPM_v_sdvS.RData"))
recov_samples <- recov_samples1

# Bring longer thetas down to min nmc by sampling
nmcs <- sapply(recov_samples, function(x) x$nmc)
nmc <- min(nmcs)
for (i in 1:length(recov_samples)) if (nmcs[i] > nmc) recov_samples[[i]]$theta <-
  recov_samples[[i]]$theta[,,sample(1:dim(recov_samples[[i]]$theta)[3], nmc)]
recov_samps <- lapply(recov_samples, function(x) x["theta"])

# # Get summary
# recov_post_summ <- get.participant.median.CIs(recov_samps)
# recov_post_summ
# 
# # Save
# save(recov_post_summ, file = "deriv/recov_post_summ_v_sdvS.RData")
print(load("deriv/recov_post_summ_v_sdvS.RData"))

# Load true parameter vector
print(load("deriv/recov_p_vector_v_sdvS.RData"))
sim_pop_mu

recov_msds <- data.frame(get.msds(recov_samples))
recov_msds
top_msds

# Get data frames for ggplot
mimic_A <- get.ggdf.mimic(recov_post_summ, top_msds, recov_msds, sim_pop_mu, grepchr = "^A", n = 100)
mimic_B <- get.ggdf.mimic(recov_post_summ, top_msds, recov_msds, sim_pop_mu, grepchr = "^B", n = 100)
mimic_t0 <- get.ggdf.mimic(recov_post_summ, top_msds, recov_msds, sim_pop_mu, grepchr = "^t0", n = 100)
mimic_mean_v <- get.ggdf.mimic(recov_post_summ, top_msds, recov_msds, sim_pop_mu, grepchr = "^mean_v", n = 100)
mimic_sd_v <- get.ggdf.mimic(recov_post_summ, top_msds, recov_msds, sim_pop_mu, grepchr = "^sd_v", n = 100)

# Plot recovery
ggplot.mimic(mimic_A, ncol = 1) + xlab("") + 
  theme_classic()

ggplot.mimic(mimic_B, ncol = 4) + xlab("") + 
  theme_classic()

ggplot.mimic(mimic_t0, ncol = 1) + xlab("") + 
  theme_classic()

ggplot.mimic(mimic_mean_v, ncol = 5) + xlab("") + 
  theme_classic()

ggplot.mimic(mimic_sd_v, ncol = 2) + xlab("") + 
  theme_classic()


# -------------------------------------------------------------------------

# No drift rate effects ---------------------------------------------------

# Load samples object from reduced model
print(load("samples/recov_sTPPM_B_sdvS.RData"))
recov_samples <- recov_samples1

# Bring longer thetas down to min nmc by sampling
nmcs <- sapply(recov_samples, function(x) x$nmc)
nmc <- min(nmcs)
for (i in 1:length(recov_samples)) if (nmcs[i] > nmc) recov_samples[[i]]$theta <-
  recov_samples[[i]]$theta[,,sample(1:dim(recov_samples[[i]]$theta)[3], nmc)]
recov_samps <- lapply(recov_samples, function(x) x["theta"])

# # Get summary
# recov_post_summ <- get.participant.median.CIs(recov_samps)
# recov_post_summ
# 
# # Save
# save(recov_post_summ, file = "deriv/recov_post_summ_B_sdvS.RData")
print(load("deriv/recov_post_summ_B_sdvS.RData"))

# Load true parameter vector
print(load("deriv/recov_p_vector_B_sdvS.RData"))
sim_pop_mu

recov_msds <- data.frame(get.msds(recov_samples))
recov_msds
top_msds

# Get data frames for ggplot
mimic_A <- get.ggdf.mimic(recov_post_summ, top_msds, recov_msds, sim_pop_mu, grepchr = "^A", n = 100)
mimic_B <- get.ggdf.mimic(recov_post_summ, top_msds, recov_msds, sim_pop_mu, grepchr = "^B", n = 100)
mimic_t0 <- get.ggdf.mimic(recov_post_summ, top_msds, recov_msds, sim_pop_mu, grepchr = "^t0", n = 100)
mimic_mean_v <- get.ggdf.mimic(recov_post_summ, top_msds, recov_msds, sim_pop_mu, grepchr = "^mean_v", n = 100)
mimic_sd_v <- get.ggdf.mimic(recov_post_summ, top_msds, recov_msds, sim_pop_mu, grepchr = "^sd_v", n = 100)

# Plot recovery
ggplot.mimic(mimic_A, ncol = 1) + xlab("") + 
  theme_classic()

ggplot.mimic(mimic_B, ncol = 4) + xlab("") + 
  theme_classic()

ggplot.mimic(mimic_t0, ncol = 1) + xlab("") + 
  theme_classic()

ggplot.mimic(mimic_mean_v, ncol = 5) + xlab("") + 
  theme_classic()

ggplot.mimic(mimic_sd_v, ncol = 2) + xlab("") + 
  theme_classic()

