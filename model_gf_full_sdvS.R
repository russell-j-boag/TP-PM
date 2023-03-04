# Clear workspace
rm(list=ls())

# Set working directory to top-level folder containing DMC
# setwd("~/DMC")
# setwd("D:/ATC experiments/Exp5 TP-PM/Analysis/DMC")

# Load libraries and DMC functions
source("dmc/dmc.R")
load_model("LBA", "lbaN_B_gf.R")


# Load data
print(load("Data/ATC_TP-PM_clean.RData"))
dat <- cleandats[,c("s","S","TP","PM","R","RT")]
# names(dat)[c(1)] <- c("subjects")
levels(dat$R)
levels(dat$S)
head(dat)
str(dat)
# NB: RT has been truncated below 0.20s

# Rename PM factor levels
levels(dat$PM)
dat$PM <- factor(dat$PM, levels = c("10%", "30%"), labels = c("10", "30"))

# Check factor levels
lapply(dat, levels)

# 3200 trials per subject
table(dat$s)


# Create match maps -------------------------------------------------------

# expand.grid(list(S = c("cc","nn","pc","pn"),
#                  TP = c("3s","6s"),
#                  PM = c("10","30"),
#                  R = c("C","N","P")))y

# Bottom-up inputs map
map_I <- empty.map(
  
  list(S = c("cc","nn","pc","pn"),
       TP = c("3s","6s"),
       PM = c("10","30"),
       R = c("C","N","P")),
  
  levels = c("cc3s10C","nn3s10C","pc3s10C","pn3s10C",
             "cc6s10C","nn6s10C","pc6s10C","pn6s10C",
             "cc3s30C","nn3s30C","pc3s30C","pn3s30C",
             "cc6s30C","nn6s30C","pc6s30C","pn6s30C",
             
             "cc3s10N","nn3s10N","pc3s10N","pn3s10N",
             "cc6s10N","nn6s10N","pc6s10N","pn6s10N",
             "cc3s30N","nn3s30N","pc3s30N","pn3s30N",
             "cc6s30N","nn6s30N","pc6s30N","pn6s30N",
             
             "pp3s10P",
             "pp6s10P",
             "pp3s30P",
             "pp6s30P",
             
             "PMFA"))
map_I
length(map_I)
length(levels(map_I))

map_I[1:48] <- c("cc3s10C","nn3s10C","pc3s10C","pn3s10C",
                 "cc6s10C","nn6s10C","pc6s10C","pn6s10C",
                 "cc3s30C","nn3s30C","pc3s30C","pn3s30C",
                 "cc6s30C","nn6s30C","pc6s30C","pn6s30C",
                 
                 "cc3s10N","nn3s10N","pc3s10N","pn3s10N",
                 "cc6s10N","nn6s10N","pc6s10N","pn6s10N",
                 "cc3s30N","nn3s30N","pc3s30N","pn3s30N",
                 "cc6s30N","nn6s30N","pc6s30N","pn6s30N",
                 
                 "PMFA","PMFA","pp3s10P","pp3s10P",
                 "PMFA","PMFA","pp6s10P","pp6s10P",
                 "PMFA","PMFA","pp3s30P","pp3s30P",
                 "PMFA","PMFA","pp6s30P","pp6s30P")
map_I


# # Gain map
# map_g <- empty.map(
#   
#   list(S = c("cc","nn","pc","pn"),
#        TP = c("3s","6s"),
#        PM = c("10","30"),
#        R = c("C","N","P")),
#   
#   levels = c("3s10",
#              "3s30",
#              "6s10",
#              "6s30"))
# map_g
# length(map_g)
# length(levels(map_g))
# 
# map_g[1:48] <- c("3s10","3s10","3s10","3s10",
#                  "6s10","6s10","6s10","6s10",
#                  "3s30","3s30","3s30","3s30",
#                  "6s30","6s30","6s30","6s30",
#                  
#                  "3s10","3s10","3s10","3s10",
#                  "6s10","6s10","6s10","6s10",
#                  "3s30","3s30","3s30","3s30",
#                  "6s30","6s30","6s30","6s30",
#                  
#                  "3s10","3s10","3s10","3s10",
#                  "6s10","6s10","6s10","6s10",
#                  "3s30","3s30","3s30","3s30",
#                  "6s30","6s30","6s30","6s30")
# map_g


# # Focus map
# map_f <- empty.map(
#   
#   list(S = c("cc","nn","pc","pn"),
#        TP = c("3s","6s"),
#        PM = c("10","30"),
#        R = c("C","N","P")),
#   
#   levels = c("3s10",
#              "3s30",
#              "6s10",
#              "6s30"))
# map_f
# length(map_f)
# length(levels(map_f))
# 
# map_f[1:48] <- c("3s10","3s10","3s10","3s10",
#                  "6s10","6s10","6s10","6s10",
#                  "3s30","3s30","3s30","3s30",
#                  "6s30","6s30","6s30","6s30",
#                  
#                  "3s10","3s10","3s10","3s10",
#                  "6s10","6s10","6s10","6s10",
#                  "3s30","3s30","3s30","3s30",
#                  "6s30","6s30","6s30","6s30",
#                  
#                  "3s10","3s10","3s10","3s10",
#                  "6s10","6s10","6s10","6s10",
#                  "3s30","3s30","3s30","3s30",
#                  "6s30","6s30","6s30","6s30")
# map_f


# sdv map
map_sdv <- empty.map(
  
  list(S = c("cc","nn","pc","pn"),
       TP = c("3s","6s"),
       PM = c("10","30"),
       R = c("C","N","P")),
  
  levels = c("ccSDV","nnSDV","ppSDV"))
map_sdv
length(map_sdv)
length(levels(map_sdv))

map_sdv[1:48] <- c("ccSDV","nnSDV","ppSDV","ppSDV",
                   "ccSDV","nnSDV","ppSDV","ppSDV",
                   "ccSDV","nnSDV","ppSDV","ppSDV",
                   "ccSDV","nnSDV","ppSDV","ppSDV",
                   
                   "ccSDV","nnSDV","ppSDV","ppSDV",
                   "ccSDV","nnSDV","ppSDV","ppSDV",
                   "ccSDV","nnSDV","ppSDV","ppSDV",
                   "ccSDV","nnSDV","ppSDV","ppSDV",
                   
                   "ccSDV","nnSDV","ppSDV","ppSDV",
                   "ccSDV","nnSDV","ppSDV","ppSDV",
                   "ccSDV","nnSDV","ppSDV","ppSDV",
                   "ccSDV","nnSDV","ppSDV","ppSDV")
map_sdv


# Threshold map
map_B <- empty.map(
  
  list(S = c("cc","nn","pc","pn"),
       TP = c("3s","6s"),
       PM = c("10","30"),
       R = c("C","N","P")),
  
  levels = c("3s10C",
             "3s30C",
             "6s10C",
             "6s30C",
             
             "3s10N",
             "3s30N",
             "6s10N",
             "6s30N",
             
             "3s10P",
             "3s30P",
             "6s10P",
             "6s30P"))
map_B
length(map_B)
length(levels(map_B))

map_B[1:48] <- c("3s10C","3s10C","3s10C","3s10C",
                 "6s10C","6s10C","6s10C","6s10C",
                 "3s30C","3s30C","3s30C","3s30C",
                 "6s30C","6s30C","6s30C","6s30C",
                 
                 "3s10N","3s10N","3s10N","3s10N",
                 "6s10N","6s10N","6s10N","6s10N",
                 "3s30N","3s30N","3s30N","3s30N",
                 "6s30N","6s30N","6s30N","6s30N",
                 
                 "3s10P","3s10P","3s10P","3s10P",
                 "6s10P","6s10P","6s10P","6s10P",
                 "3s30P","3s30P","3s30P","3s30P",
                 "6s30P","6s30P","6s30P","6s30P")
map_B



# Build model -------------------------------------------------------------

model <- model.dmc(
  p.map = list(
    A = "1",
    B = "MAPB",
    t0 = "1",
    I = "MAPI",
    g = c("TP", "PM"),
    f = c("TP", "PM"),
    sd_v = "MAPSDV",
    st0 = "1",
    N = "PM"), 
  match.map = list(
    M = list(
      cc = "C", 
      nn = "N", 
      pc = "P", 
      pn = "P"),
    MAPB = map_B,
    MAPSDV = map_sdv,
    MAPI = map_I),
  factors = list(
    S = c("cc", "nn", "pc", "pn"),
    TP = c("3s", "6s"),
    PM = c("10", "30")),
  constants = c(N.10 = 3, N.30 = 3, 
                st0 = 0,
                sd_v.ppSDV = 0.5), 
  responses = c("C","N","P"),
  type = "normN")

length(attr(model, "p.vector"))

# Create parameter vector
p.vector <- c(t0 = 0.3, A = 1.5,
              
              B.3s10C = 2,  B.3s30C = 2,  B.6s10C = 2,  B.6s30C = 2,  
              B.3s10N = 2,  B.3s30N = 2,  B.6s10N = 2,  B.6s30N = 2, 
              B.3s10P = 2,  B.3s30P = 2,  B.6s10P = 2,  B.6s30P = 2,
              
              I.cc3s10C = 1, I.nn3s10C = 0, I.pc3s10C = 0, I.pn3s10C = 0, 
              I.cc6s10C = 1, I.nn6s10C = 0, I.pc6s10C = 0, I.pn6s10C = 0,
              I.cc3s30C = 1, I.nn3s30C = 0, I.pc3s30C = 0, I.pn3s30C = 0, 
              I.cc6s30C = 1, I.nn6s30C = 0, I.pc6s30C = 0, I.pn6s30C = 0,
              I.cc3s10N = 0, I.nn3s10N = 1, I.pc3s10N = 0, I.pn3s10N = 0,
              I.cc6s10N = 0, I.nn6s10N = 1, I.pc6s10N = 0, I.pn6s10N = 0,
              I.cc3s30N = 0, I.nn3s30N = 1, I.pc3s30N = 0, I.pn3s30N = 0,
              I.cc6s30N = 0, I.nn6s30N = 1, I.pc6s30N = 0, I.pn6s30N = 0,
              I.pp3s10P = 1, I.pp6s10P = 1, I.pp3s30P = 1, I.pp6s30P = 1,
              I.PMFA = 0,
              
              g.3s.10 = 1.5, g.6s.10 = 1.5, g.3s.30 = 1.5, g.6s.30 = 1.5,    
              f.3s.10 = 1.5, f.6s.10 = 1.5, f.3s.30 = 1.5, f.6s.30 = 1.5,
              
              sd_v.ccSDV = 0.5, sd_v.nnSDV = 0.5)

length(p.vector)

# Check parameter vector matches model
check.p.vector(p.vector, model)

# Check model simulates
# simulate.dmc(p.vector, model)

# Set priors
p.prior <- prior.p.dmc(
  dists = c("beta", rep("tnorm", length(p.vector)-1)),
  p1 = c(p.vector),                           
  p2 = c(0.2, 0.1, rep(1, 12), rep(1, 37), rep(0.5, 4), rep(0.5, 4), rep(0.5, 2)), 
  lower = c(0.1, 0, rep(0, 12), rep(NA, 37), rep(0, 4), rep(0, 4), rep(0, 2)),
  upper = c(1, 10, rep(Inf, 12), rep(Inf, 37), rep(Inf, 4), rep(Inf, 4), rep(Inf, 2))
)
length(p.prior)
length(p.vector)

# Plot priors
# par(mfcol = c(2, 4)); for (i in names(p.prior)) plot.prior(i, p.prior)


# Make data model
dm <- data.model.dmc(dat, model)
save(dm, file = "samples/dmTPPM_gf_full_sdvS.RData")


# Initialize samples object
# 183 chains
n.chains <- length(p.prior) * 3

# Generate start points for fixed effects model
samples <- h.samples.dmc(nmc = 100, 
                         p.prior = p.prior, 
                         data = dm, 
                         thin = 10, 
                         n.chains = n.chains)

# Save
save(samples, file = "samples/sTPPM_gf_full_sdvS.RData")

# Load
print(load("samples/sTPPM_gf_full_sdvS.RData"))

# -------------------------------------------------------------------------
# 
# # Generate start points for hierarchical model
# # Hyper-level
# hstart <- make.hstart(samples)
# 
# # Subject-level
# theta <- make.theta1(samples)
# 
# 
# # Hyper-level priors
# # Mu
# mu.prior <- prior.p.dmc(
#   dists = rep("tnorm", 51),
#   p1 = c(t0 = 1, p.vector[-1]),
#   p2 = c(1, 1, rep(1, 12), rep(1.5, 37)), 
#   lower = c(0.1, 0, rep(0, 12), rep(NA, 37)),
#   upper = c(1, 10, rep(Inf, length(p.vector)-2))
# )
# 
# # Sigma
# sigma.prior <- prior.p.dmc(
#   dists = c(rep("gamma", 51)), 
#   p1 = c(1, 1, rep(1, 12), rep(1.5, 37)),                          
#   p2 = c(rep(1, 51))
# )
# 
# # Create hyper prior object
# pp.prior <- list(mu.prior, sigma.prior)
# 
# # Initialize samples object
# n.chains <- length(p.prior) * 3
# 
# # Generate start points for hierarchical model
# hsamples <- h.samples.dmc(nmc = 100, 
#                           p.prior = p.prior, 
#                           data = dm, 
#                           pp.prior = pp.prior, 
#                           thin = 10,
#                           hstart.prior = hstart, 
#                           theta1 = theta, 
#                           n.chains = n.chains
# )
# 
# # Save
# save(hsamples, file = "dmc/samples/hsTPPM_gf_full_sdvS.RData")