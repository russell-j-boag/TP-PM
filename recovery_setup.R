# Clear workspace
rm(list = ls())

# Set working directory to top-level folder containing DMC
getwd()

# Load libraries and DMC functions
source("dmc/dmc.R")
load_model("LBA", "lbaN_B.R")
# source("dmc/dmc_ATC.R")
# source("dmc/dmc_PMDC.R")

# Load samples object
print(load("samples/sTPPM_full_sdvS.RData"))
samples <- samples1
colnames(samples[[1]]$theta)


# Set up parameters and priors --------------------------------------------

# Mean (here put your actual priors that you used to fit the model)
p.vector <- c(t0 = 0.3, A = 1.5,
              
              B.3s10C = 2,  B.3s30C = 2,  B.6s10C = 2,  B.6s30C = 2,  
              B.3s10N = 2,  B.3s30N = 2,  B.6s10N = 2,  B.6s30N = 2, 
              B.3s10P = 2,  B.3s30P = 2,  B.6s10P = 2,  B.6s30P = 2,
              
              mean_v.cc3s10C = 1, mean_v.nn3s10C = 0, mean_v.pc3s10C = 0, mean_v.pn3s10C = 0, 
              mean_v.cc6s10C = 1, mean_v.nn6s10C = 0, mean_v.pc6s10C = 0, mean_v.pn6s10C = 0,
              mean_v.cc3s30C = 1, mean_v.nn3s30C = 0, mean_v.pc3s30C = 0, mean_v.pn3s30C = 0, 
              mean_v.cc6s30C = 1, mean_v.nn6s30C = 0, mean_v.pc6s30C = 0, mean_v.pn6s30C = 0,
              mean_v.cc3s10N = 0, mean_v.nn3s10N = 1, mean_v.pc3s10N = 0, mean_v.pn3s10N = 0,
              mean_v.cc6s10N = 0, mean_v.nn6s10N = 1, mean_v.pc6s10N = 0, mean_v.pn6s10N = 0,
              mean_v.cc3s30N = 0, mean_v.nn3s30N = 1, mean_v.pc3s30N = 0, mean_v.pn3s30N = 0,
              mean_v.cc6s30N = 0, mean_v.nn6s30N = 1, mean_v.pc6s30N = 0, mean_v.pn6s30N = 0,
              mean_v.pp3s10P = 1, mean_v.pp6s10P = 1, mean_v.pp3s30P = 1, mean_v.pp6s30P = 1,
              mean_v.PMFA = 0,
              
              sd_v.ccSDV = 0.5, sd_v.nnSDV = 0.5)


# Set priors
p.prior <- prior.p.dmc(
  dists = c("beta", rep("tnorm", length(p.vector)-1)),
  p1 = c(p.vector),                           
  p2 = c(0.2, 0.1, rep(1, 12), rep(1, 37), rep(0.5, 2)), 
  lower = c(0.1, 0, rep(0, 12), rep(NA, 37), rep(0, 2)),
  upper = c(1, 10, rep(Inf, 12), rep(Inf, 37), rep(Inf, 2))
)
length(p.prior)
length(p.prior) * 3

# Get the average amount of data per participant
Ns <- lapply(samples, get.ns.dmc)
Ns

# Get the average total number of trials per participant observed in the experiment
# (this accounts for reduction in power from removed RT outlier trials etc.)
ns <- ceiling(Reduce("+", Ns)/length(Ns))
ns

# Get model from samples object
model <- attr(samples[[1]]$data, "model")
attr(model, "match.map")

# Get representative parameters from the top model
# Get posterior mean (or median) parameters of the subject-average distribution
msds <- get.msds(samples)
msds <- data.frame(msds)
msds
sim_pop_mu <- msds$M; names(sim_pop_mu) <- rownames(msds)
sim_pop_mu


# Simulate synthetic data -------------------------------------------------

# Loop through and make 100 synthetic participants with the representative data frame size 
for (i in 1:100) {
  data <- simulate.dmc(sim_pop_mu, model, n = ns)
  data <- cbind(i, data)
  if (i == 1) okdats <- data else okdats <- rbind(okdats, data)
}

# Label subjects column 's'
head(okdats)
names(okdats)[1] <- "s"
okdats$s <-factor(okdats$s)

# Make data model
dm <- data.model.dmc(okdats, model)


# Sampling ----------------------------------------------------------------

# Settings
n.chains <- length(p.prior) * 3

# Generate start points
recov_samples <- h.samples.dmc(nmc = 100, 
                               p.prior = p.prior, 
                               data = dm, 
                               thin = 10, 
                               n.chains = n.chains)

# Save samples and run on servers with "batch" file (batch_<model_name>.R)
# with the command: nohup R CMD BATCH batch_<model_name>.R &
save(recov_samples, file = "samples/recov_sTPPM_full_sdvS.RData")

# Run fitting on grid with batch script



# -------------------------------------------------------------------------

# Model mimicry recovery study --------------------------------------------
#
# 1. Generate data from a model with certain effects 
#    deleted (by changing the p.vector you simulate from)
# 2. Fit the top model to that and check whether the
#    deleted effects spontaneously emerge from nothing. 
#    If so, we should be suspicious of them.

# Setup environment -------------------------------------------------------

# Clear workspace
rm(list = ls())

# Set working directory
getwd()

# Load libraries and DMC functions
source("dmc/dmc.R")
load_model("LBA", "lbaN_B.R")

# Load top samples object
print(load("samples/sTPPM_full_sdvS.RData"))
top_samples <- samples1
colnames(top_samples[[1]]$theta)


# Set up parameters and priors --------------------------------------------

# Mean (here put your actual priors that you used to fit the model)
p.vector <- c(t0 = 0.3, A = 1.5,
              
              B.3s10C = 2,  B.3s30C = 2,  B.6s10C = 2,  B.6s30C = 2,  
              B.3s10N = 2,  B.3s30N = 2,  B.6s10N = 2,  B.6s30N = 2, 
              B.3s10P = 2,  B.3s30P = 2,  B.6s10P = 2,  B.6s30P = 2,
              
              mean_v.cc3s10C = 1, mean_v.nn3s10C = 0, mean_v.pc3s10C = 0, mean_v.pn3s10C = 0, 
              mean_v.cc6s10C = 1, mean_v.nn6s10C = 0, mean_v.pc6s10C = 0, mean_v.pn6s10C = 0,
              mean_v.cc3s30C = 1, mean_v.nn3s30C = 0, mean_v.pc3s30C = 0, mean_v.pn3s30C = 0, 
              mean_v.cc6s30C = 1, mean_v.nn6s30C = 0, mean_v.pc6s30C = 0, mean_v.pn6s30C = 0,
              mean_v.cc3s10N = 0, mean_v.nn3s10N = 1, mean_v.pc3s10N = 0, mean_v.pn3s10N = 0,
              mean_v.cc6s10N = 0, mean_v.nn6s10N = 1, mean_v.pc6s10N = 0, mean_v.pn6s10N = 0,
              mean_v.cc3s30N = 0, mean_v.nn3s30N = 1, mean_v.pc3s30N = 0, mean_v.pn3s30N = 0,
              mean_v.cc6s30N = 0, mean_v.nn6s30N = 1, mean_v.pc6s30N = 0, mean_v.pn6s30N = 0,
              mean_v.pp3s10P = 1, mean_v.pp6s10P = 1, mean_v.pp3s30P = 1, mean_v.pp6s30P = 1,
              mean_v.PMFA = 0,
              
              sd_v.ccSDV = 0.5, sd_v.nnSDV = 0.5)


# Set priors
p.prior <- prior.p.dmc(
  dists = c("beta", rep("tnorm", length(p.vector)-1)),
  p1 = c(p.vector),                           
  p2 = c(0.2, 0.1, rep(1, 12), rep(1, 37), rep(0.5, 2)), 
  lower = c(0.1, 0, rep(0, 12), rep(NA, 37), rep(0, 2)),
  upper = c(1, 10, rep(Inf, 12), rep(Inf, 37), rep(Inf, 2))
)
length(p.prior)
length(p.prior) * 3

# Get the average amount of data per participant
Ns <- lapply(top_samples, get.ns.dmc)
Ns

# Get the average total number of trials per participant observed in the experiment
# (this accounts for reduction in power from removed RT outlier trials etc.)
ns <- ceiling(Reduce("+", Ns)/length(Ns))
ns

# Get model from your samples object
model <- attr(top_samples[[1]]$data, "model")
attr(model, "match.map")

# Get representative parameters from the top model
# Get posterior mean (or median) parameters of the subject-average distribution
top_msds <- get.msds(top_samples)
top_msds <- data.frame(top_msds)
top_msds
sim_pop_mu <- top_msds$M; names(sim_pop_mu) <- rownames(top_msds)
sim_pop_mu



# -------------------------------------------------------------------------

# No threshold effects recovery model -------------------------------------
#
# Here we use the model fitted without threshold effects as the source 
# of the parameters to get parameter values that actually could have fit 
# our study without threshold effects
# 

# Load samples object from reduced model
print(load("samples/sTPPM_v_sdvS.RData"))
replacement_samples <- samples1
colnames(replacement_samples[[1]]$theta)

# Get replacement parameters
replacement_msds <- get.msds(replacement_samples) # Reduced model
replacement_pop_mu <- data.frame(replacement_msds)$M
names(replacement_pop_mu) <- rownames(replacement_msds)
replacement_pop_mu
sim_pop_mu
length(replacement_pop_mu)
length(sim_pop_mu)

# Replace full model parms with reduced model parms
sim_pop_mu[grep(names(sim_pop_mu), pattern = "^A")] <- 
  replacement_pop_mu[grep(names(replacement_pop_mu), pattern = "^A")]  # A

# Conflict thresholds
sim_pop_mu[grep(names(sim_pop_mu), pattern = "^B.*C")] <- 
  replacement_pop_mu[grep(names(replacement_pop_mu), pattern = "^B.*C")]

# Non-conflict thresholds
sim_pop_mu[grep(names(sim_pop_mu), pattern = "^B.*N")] <- 
  replacement_pop_mu[grep(names(replacement_pop_mu), pattern = "^B.*N")]

# PM thresholds
sim_pop_mu[grep(names(sim_pop_mu), pattern = "^B.*P")] <- 
  replacement_pop_mu[grep(names(replacement_pop_mu), pattern = "^B.*P")]

sim_pop_mu[grep(names(sim_pop_mu), pattern = "^t0")] <- 
  replacement_pop_mu[grep(names(replacement_pop_mu), pattern = "^t0")]  # t0

# Rates
# names(sim_pop_mu[grep(names(sim_pop_mu), pattern = "^mean_v")]) ==
#   names(replacement_pop_mu[grep(names(replacement_pop_mu), pattern = "^mean_v")])

sim_pop_mu[grep(names(sim_pop_mu), pattern = "^mean_v")] <- 
  replacement_pop_mu[grep(names(replacement_pop_mu), pattern = "^mean_v")]

sim_pop_mu[grep(names(sim_pop_mu), pattern = "^sd_v")] <- 
  replacement_pop_mu[grep(names(replacement_pop_mu), pattern = "^sd_v")]  # sdv

save(sim_pop_mu, file = "deriv/recov_p_vector_v_sdvS.RData")

# Simulate synthetic data -------------------------------------------------

# Loop through and make 100 synthetic participants with the representative data frame size 
for (i in 1:100) {
  data <- simulate.dmc(sim_pop_mu, model, n = ns)
  data <- cbind(i, data)
  if (i == 1) okdats <- data else okdats <- rbind(okdats, data)
}

# Label subjects column 's'
head(okdats)
names(okdats)[1] <- "s"
okdats$s <-factor(okdats$s)

# Make data model
dm <- data.model.dmc(okdats, model)


# Sampling ----------------------------------------------------------------

# Settings
n.chains <- length(p.prior) * 3

# Set up a samples (start points) object to fit to the synthetic data
recov_samples <- h.samples.dmc(nmc = 100, 
                               p.prior = p.prior, 
                               data = dm, 
                               thin = 10, 
                               n.chains = n.chains
)

save(recov_samples,  file = "samples/recov_sTPPM_v_sdvS.RData")

# Run fitting on grid with batch script



# -------------------------------------------------------------------------

# No drift rate effects recovery model ------------------------------------
#
# Use the model fitted without drift rate effects as the source of the 
# parameters to get parameter values that actually could have fit our 
# study without drift rate effects
# 

# Load samples object from reduced model
print(load("samples/sTPPM_B_sdvS.RData"))
replacement_samples <- samples1
colnames(replacement_samples[[1]]$theta)

# Get replacement parameters
replacement_msds <- get.msds(replacement_samples) # Reduced model
replacement_pop_mu <- data.frame(replacement_msds)$M
names(replacement_pop_mu) <- rownames(replacement_msds)
replacement_pop_mu
sim_pop_mu
length(replacement_pop_mu)
length(sim_pop_mu)

# Replace full model parameters with reduced model parameters
sim_pop_mu[grep(names(sim_pop_mu), pattern = "^A")] <- 
  replacement_pop_mu[grep(names(replacement_pop_mu), pattern = "^A")]  # A

sim_pop_mu[grep(names(sim_pop_mu), pattern = "^B")] <- 
  replacement_pop_mu[grep(names(replacement_pop_mu), pattern = "^B")]  # B

sim_pop_mu[grep(names(sim_pop_mu), pattern = "^t0")] <- 
  replacement_pop_mu[grep(names(replacement_pop_mu), pattern = "^t0")]  # t0

# Conflict rates
sim_pop_mu[grep(names(sim_pop_mu), pattern = "^mean_v.cc.*C")] <- 
  replacement_pop_mu[grep(names(replacement_pop_mu), pattern = "^mean_v.cc.*C")]
sim_pop_mu[grep(names(sim_pop_mu), pattern = "^mean_v.nn.*C")] <- 
  replacement_pop_mu[grep(names(replacement_pop_mu), pattern = "^mean_v.nn.*C")]
sim_pop_mu[grep(names(sim_pop_mu), pattern = "^mean_v.pc.*C")] <- 
  replacement_pop_mu[grep(names(replacement_pop_mu), pattern = "^mean_v.pc.*C")]
sim_pop_mu[grep(names(sim_pop_mu), pattern = "^mean_v.pn.*C")] <- 
  replacement_pop_mu[grep(names(replacement_pop_mu), pattern = "^mean_v.pn.*C")]

# Non-conflict rates
sim_pop_mu[grep(names(sim_pop_mu), pattern = "^mean_v.cc.*N")] <- 
  replacement_pop_mu[grep(names(replacement_pop_mu), pattern = "^mean_v.cc.*N")]
sim_pop_mu[grep(names(sim_pop_mu), pattern = "^mean_v.nn.*N")] <- 
  replacement_pop_mu[grep(names(replacement_pop_mu), pattern = "^mean_v.nn.*N")]
sim_pop_mu[grep(names(sim_pop_mu), pattern = "^mean_v.pc.*N")] <- 
  replacement_pop_mu[grep(names(replacement_pop_mu), pattern = "^mean_v.pc.*N")]
sim_pop_mu[grep(names(sim_pop_mu), pattern = "^mean_v.pn.*N")] <- 
  replacement_pop_mu[grep(names(replacement_pop_mu), pattern = "^mean_v.pn.*N")]

# PM rates
sim_pop_mu[grep(names(sim_pop_mu), pattern = "^mean_v.pp.*P")] <- 
  replacement_pop_mu[grep(names(replacement_pop_mu), pattern = "^mean_v.pp.*P")]
sim_pop_mu[grep(names(sim_pop_mu), pattern = "^mean_v.PMFA")] <- 
  replacement_pop_mu[grep(names(replacement_pop_mu), pattern = "^mean_v.PMFA")]

sim_pop_mu[grep(names(sim_pop_mu), pattern = "^sd_v")] <- 
  replacement_pop_mu[grep(names(replacement_pop_mu), pattern = "^sd_v")]  # sdv

save(sim_pop_mu, file = "deriv/recov_p_vector_B_sdvS.RData")

# Simulate synthetic data -------------------------------------------------

# Loop through and make 100 synthetic participants with the representative data frame size 
for (i in 1:100) {
  data <- simulate.dmc(sim_pop_mu, model, n = ns)
  data <- cbind(i, data)
  if (i == 1) okdats <- data else okdats <- rbind(okdats, data)
}

# Label subjects column 's'
head(okdats)
names(okdats)[1] <- "s"
okdats$s <-factor(okdats$s)

# Make data model
dm <- data.model.dmc(okdats, model)


# Sampling ----------------------------------------------------------------

# Settings
n.chains <- length(p.prior) * 3

# Set up a samples (start points) object to fit to the synthetic data
recov_samples <- h.samples.dmc(nmc = 100, 
                               p.prior = p.prior, 
                               data = dm, 
                               thin = 10, 
                               n.chains = n.chains
)

save(recov_samples,  file = "samples/recov_sTPPM_B_sdvS.RData")

# Run fitting on grid with batch script
