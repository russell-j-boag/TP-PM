# Clear workspace
rm(list = ls())

# Set working directory to top-level folder containing DMC
getwd()

source("dmc/dmc.R")
source("exploration_functions.R")
load_model("LBA", "lbaN_B.R")
# require("stringr")
# require("cowplot")
# source("dmc/dmc_ATC.R")
# source("dmc/dmc_PMDC.R")

# Load top samples object
print(load("samples/sTPPM_full_sdvS.RData"))
top_samples <- samples1
colnames(top_samples[[1]]$theta)

# Load top model posterior predictives
# post_pred_sims <- h.post.predict.dmc(top_samples, save.simulation = TRUE)
# save(post_pred_sims, file = "deriv/post_pred_sims_full_sdvS.RData")
print(load("deriv/post_pred_sims_full_sdvS.RData"))

sim <- do.call(rbind, post_pred_sims)
currentsim <- sim

# Test summary functions
RP_effects_ongoing(currentsim)
RT_effects_ongoing(currentsim)
RP_effects_PM(currentsim)
RT_effects_PM(currentsim)

RP_effects_ongoing_top <- get.effects.dmc(post_pred_sims, fun = RP_effects_ongoing)
RT_effects_ongoing_top <- get.effects.dmc(post_pred_sims, fun = RT_effects_ongoing)
RP_effects_PM_top <- get.effects.dmc(post_pred_sims, fun = RP_effects_PM)
RT_effects_PM_top <- get.effects.dmc(post_pred_sims, fun = RT_effects_PM)

# Get p.names
top_samples[[1]]$p.names
#  [1] "A"              "B.3s10C"        "B.3s30C"        "B.6s10C"        "B.6s30C"        "B.3s10N"       
#  [7] "B.3s30N"        "B.6s10N"        "B.6s30N"        "B.3s10P"        "B.3s30P"        "B.6s10P"       
# [13] "B.6s30P"        "t0"             "mean_v.cc3s10C" "mean_v.nn3s10C" "mean_v.pc3s10C" "mean_v.pn3s10C"
# [19] "mean_v.cc6s10C" "mean_v.nn6s10C" "mean_v.pc6s10C" "mean_v.pn6s10C" "mean_v.cc3s30C" "mean_v.nn3s30C"
# [25] "mean_v.pc3s30C" "mean_v.pn3s30C" "mean_v.cc6s30C" "mean_v.nn6s30C" "mean_v.pc6s30C" "mean_v.pn6s30C"
# [31] "mean_v.cc3s10N" "mean_v.nn3s10N" "mean_v.pc3s10N" "mean_v.pn3s10N" "mean_v.cc6s10N" "mean_v.nn6s10N"
# [37] "mean_v.pc6s10N" "mean_v.pn6s10N" "mean_v.cc3s30N" "mean_v.nn3s30N" "mean_v.pc3s30N" "mean_v.pn3s30N"
# [43] "mean_v.cc6s30N" "mean_v.nn6s30N" "mean_v.pc6s30N" "mean_v.pn6s30N" "mean_v.pp3s10P" "mean_v.pp6s10P"
# [49] "mean_v.pp3s30P" "mean_v.pp6s30P" "mean_v.PMFA"    "sd_v.ccSDV"     "sd_v.nnSDV"

# First generate posterior predictives with certain parameters averaged across conditions.
# avps.post.predict.dmc & avps.h.post.predict.dmc
# grep match all parameters to the av.posts vector,
# and ones that match each other are averaged.


# -------------------------------------------------------------------------

# Average drift rates over all conditions
av_posts_v <- c(
  "^mean_v.cc*C",
  "^mean_v.nn*C",
  "^mean_v.pc*C",
  "^mean_v.pn*C",
  
  "^mean_v.cc*N",
  "^mean_v.nn*N",
  "^mean_v.pc*N",
  "^mean_v.pn*N",
  
  "mean_v.pp*P"
)

# This converts your av_posts to regex using glob2rx
av_posts_v <- glob2rx(av_posts_v)
av_posts_v

post_preds_fixed_v <- avps.h.post.predict.dmc(top_samples, 
                                              n.post = 100, 
                                              av.posts = av_posts_v,
                                              save.simulation = TRUE)

save(post_preds_fixed_v, file = "deriv/av_posts_full_sdvS_fixed_v.RData")
print(load("deriv/av_posts_full_sdvS_fixed_v.RData"))

# get.effects.dmc gets whatever effects are specified by fun for the data, and for each rep of the sim. 
# then it calculates posterior mean and quantiles
# the output is a data frame with post.mean, credible intervals of the effects plus data effect.

RP_effects_ongoing_fixed_v <- get.effects.dmc(post_preds_fixed_v, fun = RP_effects_ongoing)
RT_effects_ongoing_fixed_v <- get.effects.dmc(post_preds_fixed_v, fun = RT_effects_ongoing)
RP_effects_PM_fixed_v <- get.effects.dmc(post_preds_fixed_v, fun = RP_effects_PM)
RT_effects_PM_fixed_v <- get.effects.dmc(post_preds_fixed_v, fun = RT_effects_PM)


# -------------------------------------------------------------------------

# Average thresholds over all conditions
av_posts_B <- c(
  "^B.*C",
  "^B.*N",
  "^B.*P"
)

# This converts your av_posts to regex using glob2rx
av_posts_B <- glob2rx(av_posts_B)
av_posts_B

post_preds_fixed_B <- avps.h.post.predict.dmc(top_samples, 
                                              n.post = 100, 
                                              av.posts = av_posts_B,
                                              save.simulation = TRUE)

save(post_preds_fixed_B, file = "deriv/av_posts_full_sdvS_fixed_B.RData")
print(load("deriv/av_posts_full_sdvS_fixed_B.RData"))

# get.effects.dmc gets whatever effects are specified by fun for the data, and for each rep of the sim. 
# then it calculates posterior mean and quantiles
# the output is a data frame with post.mean, credible intervals of the effects plus data effect.

RP_effects_ongoing_fixed_B <- get.effects.dmc(post_preds_fixed_B, fun = RP_effects_ongoing)
RT_effects_ongoing_fixed_B <- get.effects.dmc(post_preds_fixed_B, fun = RT_effects_ongoing)
RP_effects_PM_fixed_B <- get.effects.dmc(post_preds_fixed_B, fun = RP_effects_PM)
RT_effects_PM_fixed_B <- get.effects.dmc(post_preds_fixed_B, fun = RT_effects_PM)


# Save
save(
  RP_effects_ongoing_top, 
  RT_effects_ongoing_top,
  RP_effects_PM_top, 
  RT_effects_PM_top,
  RP_effects_ongoing_fixed_v,
  RT_effects_ongoing_fixed_v,
  RP_effects_PM_fixed_v,
  RT_effects_PM_fixed_v,
  RP_effects_ongoing_fixed_B,
  RT_effects_ongoing_fixed_B,
  RP_effects_PM_fixed_B,
  RT_effects_PM_fixed_B,
  file = "deriv/av_posts_full_sdvS_effects.RData"
)
print(load("deriv/av_posts_full_sdvS_effects.RData"))

# RT_effects_ongoing_top
# RT_effects_ongoing_fixed_v


# -------------------------------------------------------------------------

# Here we don't average but instead remove various model mechanisms

# Turn off reactive control

# "mean_v.cc3s10C" "mean_v.nn3s10C" "mean_v.pc3s10C" "mean_v.pn3s10C"
# "mean_v.cc6s10C" "mean_v.nn6s10C" "mean_v.pc6s10C" "mean_v.pn6s10C" 
# "mean_v.cc3s30C" "mean_v.nn3s30C" "mean_v.pc3s30C" "mean_v.pn3s30C" 
# "mean_v.cc6s30C" "mean_v.nn6s30C" "mean_v.pc6s30C" "mean_v.pn6s30C"
# "mean_v.cc3s10N" "mean_v.nn3s10N" "mean_v.pc3s10N" "mean_v.pn3s10N" 
# "mean_v.cc6s10N" "mean_v.nn6s10N" "mean_v.pc6s10N" "mean_v.pn6s10N" 
# "mean_v.cc3s30N" "mean_v.nn3s30N" "mean_v.pc3s30N" "mean_v.pn3s30N"
# "mean_v.cc6s30N" "mean_v.nn6s30N" "mean_v.pc6s30N" "mean_v.pn6s30N" 
# "mean_v.pp3s10P" "mean_v.pp6s10P" "mean_v.pp3s30P" "mean_v.pp6s30P"

# Original parameters to replace
reactive <- c("mean_v.cc3s10C", "mean_v.nn3s10C", "mean_v.pc3s10C", "mean_v.pn3s10C",
                "mean_v.cc6s10C", "mean_v.nn6s10C", "mean_v.pc6s10C", "mean_v.pn6s10C",
                "mean_v.cc3s30C", "mean_v.nn3s30C", "mean_v.pc3s30C", "mean_v.pn3s30C", 
                "mean_v.cc6s30C", "mean_v.nn6s30C", "mean_v.pc6s30C", "mean_v.pn6s30C",
                "mean_v.cc3s10N", "mean_v.nn3s10N", "mean_v.pc3s10N", "mean_v.pn3s10N",
                "mean_v.cc6s10N", "mean_v.nn6s10N", "mean_v.pc6s10N", "mean_v.pn6s10N", 
                "mean_v.cc3s30N", "mean_v.nn3s30N", "mean_v.pc3s30N", "mean_v.pn3s30N",
                "mean_v.cc6s30N", "mean_v.nn6s30N", "mean_v.pc6s30N", "mean_v.pn6s30N", 
                "mean_v.pp3s10P", "mean_v.pp6s10P", "mean_v.pp3s30P", "mean_v.pp6s30P")

# Replacement parameters
no_reactive <- c("mean_v.cc3s10C", "mean_v.nn3s10C", "mean_v.cc3s10C", "mean_v.nn3s10C",
                   "mean_v.cc6s10C", "mean_v.nn6s10C", "mean_v.cc6s10C", "mean_v.nn6s10C",
                   "mean_v.cc3s30C", "mean_v.nn3s30C", "mean_v.cc3s30C", "mean_v.nn3s30C", 
                   "mean_v.cc6s30C", "mean_v.nn6s30C", "mean_v.cc6s30C", "mean_v.nn6s30C",
                   "mean_v.cc3s10N", "mean_v.nn3s10N", "mean_v.cc3s10N", "mean_v.nn3s10N",
                   "mean_v.cc6s10N", "mean_v.nn6s10N", "mean_v.cc6s10N", "mean_v.nn6s10N", 
                   "mean_v.cc3s30N", "mean_v.nn3s30N", "mean_v.cc3s30N", "mean_v.nn3s30N",
                   "mean_v.cc6s30N", "mean_v.nn6s30N", "mean_v.cc6s30N", "mean_v.nn6s30N", 
                   "mean_v.pp3s10P", "mean_v.pp6s10P", "mean_v.pp3s30P", "mean_v.pp6s30P")

# Get posterior predictives
post_preds_no_reactive <- pickps.h.post.predict.dmc(top_samples, 
                                                      n.post = 100,
                                                      pickps_others = reactive, 
                                                      pickps_set = no_reactive, 
                                                      save.simulation = TRUE)
save(post_preds_no_reactive, file = "deriv/mechanism_test_full_sdvS_no_reactive.RData")

# Load posterior predictives
print(load("deriv/mechanism_test_full_sdvS_no_reactive.RData"))

RP_effects_ongoing_no_reactive <- get.effects.dmc(post_preds_no_reactive, fun = RP_effects_ongoing)
RT_effects_ongoing_no_reactive <- get.effects.dmc(post_preds_no_reactive, fun = RT_effects_ongoing)
RP_effects_PM_no_reactive <- get.effects.dmc(post_preds_no_reactive, fun = RP_effects_PM)
RT_effects_PM_no_reactive <- get.effects.dmc(post_preds_no_reactive, fun = RT_effects_PM)


# -------------------------------------------------------------------------

# Turn off proactive control (across TP levels)

# "B.3s10C" "B.3s30C" "B.6s10C" "B.6s30C"        
# "B.3s10N" "B.3s30N" "B.6s10N" "B.6s30N"        
# "B.3s10P" "B.3s30P" "B.6s10P" "B.6s30P"

# Original parameters to replace
proactive <- c("B.3s10C", "B.3s30C", "B.6s10C", "B.6s30C",        
                 "B.3s10N", "B.3s30N", "B.6s10N", "B.6s30N",       
                 "B.3s10P", "B.3s30P", "B.6s10P", "B.6s30P")

# Replacement parameters
no_proactive_TP <- c("B.6s10C", "B.6s30C", "B.6s10C", "B.6s30C",        
                       "B.6s10N", "B.6s30N", "B.6s10N", "B.6s30N",       
                       "B.6s10P", "B.6s30P", "B.6s10P", "B.6s30P")

# Get posterior predictives
post_preds_no_proactive_TP <- pickps.h.post.predict.dmc(top_samples, 
                                                          n.post = 100,
                                                          pickps_others = proactive, 
                                                          pickps_set = no_proactive_TP, 
                                                          save.simulation = TRUE)
save(post_preds_no_proactive_TP, file = "deriv/mechanism_test_full_sdvS_no_proactive_TP.RData")

# Load posterior predictives
print(load("deriv/mechanism_test_full_sdvS_no_proactive_TP.RData"))

RP_effects_ongoing_no_proactive_TP <- get.effects.dmc(post_preds_no_proactive_TP, fun = RP_effects_ongoing)
RT_effects_ongoing_no_proactive_TP <- get.effects.dmc(post_preds_no_proactive_TP, fun = RT_effects_ongoing)
RP_effects_PM_no_proactive_TP <- get.effects.dmc(post_preds_no_proactive_TP, fun = RP_effects_PM)
RT_effects_PM_no_proactive_TP <- get.effects.dmc(post_preds_no_proactive_TP, fun = RT_effects_PM)


# -------------------------------------------------------------------------

# Turn off proactive control (across PM frequency levels)

# "B.3s10C" "B.3s30C" "B.6s10C" "B.6s30C"        
# "B.3s10N" "B.3s30N" "B.6s10N" "B.6s30N"        
# "B.3s10P" "B.3s30P" "B.6s10P" "B.6s30P"

# Original parameters to replace
proactive <- c("B.3s10C", "B.3s30C", "B.6s10C", "B.6s30C",        
                 "B.3s10N", "B.3s30N", "B.6s10N", "B.6s30N",       
                 "B.3s10P", "B.3s30P", "B.6s10P", "B.6s30P")

# Replacement parameters
no_proactive_PM <- c("B.3s10C", "B.3s10C", "B.6s10C", "B.6s10C",        
                       "B.3s10N", "B.3s10N", "B.6s10N", "B.6s10N",       
                       "B.3s10P", "B.3s10P", "B.6s10P", "B.6s10P")

# Get posterior predictives
post_preds_no_proactive_PM <- pickps.h.post.predict.dmc(top_samples, 
                                                          n.post = 100,
                                                          pickps_others = proactive, 
                                                          pickps_set = no_proactive_PM, 
                                                          save.simulation = TRUE)
save(post_preds_no_proactive_PM, file = "deriv/mechanism_test_full_sdvS_no_proactive_PM.RData")

# Load posterior predictives
print(load("deriv/mechanism_test_full_sdvS_no_proactive_PM.RData"))

RP_effects_ongoing_no_proactive_PM <- get.effects.dmc(post_preds_no_proactive_PM, fun = RP_effects_ongoing)
RT_effects_ongoing_no_proactive_PM <- get.effects.dmc(post_preds_no_proactive_PM, fun = RT_effects_ongoing)
RP_effects_PM_no_proactive_PM <- get.effects.dmc(post_preds_no_proactive_PM, fun = RP_effects_PM)
RT_effects_PM_no_proactive_PM <- get.effects.dmc(post_preds_no_proactive_PM, fun = RT_effects_PM)


# Save
save(
  RP_effects_ongoing_top, 
  RT_effects_ongoing_top,
  RP_effects_PM_top,
  RT_effects_PM_top,
  RP_effects_ongoing_no_reactive,
  RT_effects_ongoing_no_reactive,
  RP_effects_PM_no_reactive,
  RT_effects_PM_no_reactive,
  RP_effects_ongoing_no_proactive_TP,
  RT_effects_ongoing_no_proactive_TP,
  RP_effects_PM_no_proactive_TP,
  RT_effects_PM_no_proactive_TP,
  RP_effects_ongoing_no_proactive_PM,
  RT_effects_ongoing_no_proactive_PM,
  RP_effects_PM_no_proactive_PM,
  RT_effects_PM_no_proactive_PM,
  file = "deriv/mechanism_test_full_sdvS_effects.RData"
)
print(load("deriv/mechanism_test_full_sdvS_effects.RData"))



# -------------------------------------------------------------------------

# factors=NA
# n.post=100
# random=TRUE
# 
# samples <- top_samples[[1]]
# model <- attributes(samples$data)$model
# facs <- names(attr(model,"factors"))
# if (any(is.na(factors))) factors <- facs
# if (!all(factors %in% facs))
#   stop(paste("Factors argument must contain one or more of:",paste(facs,collapse=",")))
# resp <- names(attr(model,"responses"))
# ns <- table(samples$data[,facs],dnn=facs)
# n.par <- dim(samples$theta)[2]
# thetas <- matrix(aperm(samples$theta,c(3,1,2)),ncol=n.par)
# colnames(thetas) <- dimnames(samples$theta)[[2]]
# if (is.na(n.post)) use <- c(1:dim(thetas)[1]) else {
#   if (random) use <- sample(c(1:dim(thetas)[1]),n.post,replace=F) else
#     use <- round(seq(1,dim(thetas)[1],length.out=n.post))
# }
# n.post <- length(use)
# posts <- thetas[use,]
# str(posts)
# 
# # Substitution
# pickps_others=t0.ref.swi
# pickps_others
# pickps_set=t0.ref.nos
# pickps_set
# colnames(posts)
# colnames(posts) %in% pickps_others
# colnames(posts) %in% pickps_set
# 
# head(posts[,colnames(posts) %in% pickps_others][,pickps_others],2)
# head(posts[,colnames(posts) %in% pickps_set][,pickps_set],2)
# 
# 
# # Replace some parameter vlaues with others.
# posts[,colnames(posts) %in% pickps_others][,pickps_others] <- 
#   posts[,colnames(posts) %in% pickps_set][,pickps_set] 
# 
# colnames(posts)
# 
