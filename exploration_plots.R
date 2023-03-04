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


# Load posterior predictive objects for reduced models
print(load("deriv/av_posts_full_sdvS_fixed_v.RData"))
print(load("deriv/av_posts_full_sdvS_fixed_B.RData"))
print(load("deriv/av_posts_full_sdvS_effects.RData"))

print(load("deriv/mechanism_test_full_sdvS_no_reactive.RData"))
print(load("deriv/mechanism_test_full_sdvS_no_proactive_TP.RData"))
print(load("deriv/mechanism_test_full_sdvS_no_proactive_PM.RData"))
print(load("deriv/mechanism_test_full_sdvS_effects.RData"))


# Plots and tables --------------------------------------------------------

# -------------------------------------------------------------------------

# Ongoing task accuracy effects with drift rates averaged over all conditions
RP_effects_ongoing_fixed_v <- data.frame(RP_effects_ongoing_fixed_v)
RP_effects_ongoing_fixed_v$effect <- NA
RP_effects_ongoing_fixed_v$effect <- factor(rownames(RP_effects_ongoing_fixed_v))
str(RP_effects_ongoing_fixed_v)
RP_effects_ongoing_fixed_v

ggplot(RP_effects_ongoing_fixed_v, aes(effect, mean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25) +
  geom_point(aes(effect, data), pch = 21, size = 4, colour = "black") +
  xlab("Effect") + ylab("Accuracy") + 
  theme_minimal() +
  ggtitle("Drift rates averaged over all conditions")


# % of data effect explained by averaged mechanism
RP_percent_ongoing_fixed_v <- (RP_effects_ongoing_fixed_v$mean/RP_effects_ongoing_fixed_v$data)*100

RP_percent_ongoing_fixed_v <- data.frame(RP_percent_ongoing_fixed_v)
rownames(RP_percent_ongoing_fixed_v) <- rownames(RP_effects_ongoing_fixed_v)
colnames(RP_percent_ongoing_fixed_v) <- "Effect %"
RP_percent_ongoing_fixed_v

RP_percent_ongoing_fixed_v$`Diff. from 100%` <- NA
RP_percent_ongoing_fixed_v$`Diff. from 100%` <- RP_percent_ongoing_fixed_v$`Effect %` - rep(100,length(RP_percent_ongoing_fixed_v$`Diff. from 100%`))
RP_percent_ongoing_fixed_v$`Effect %` <- round(RP_percent_ongoing_fixed_v$`Effect %`,2)
RP_percent_ongoing_fixed_v$`Diff. from 100%` <- round(RP_percent_ongoing_fixed_v$`Diff. from 100%`,2)
RP_percent_ongoing_fixed_v


# -------------------------------------------------------------------------

# Ongoing task RT effects with drift rates averaged over all conditions
RT_effects_ongoing_fixed_v <- data.frame(RT_effects_ongoing_fixed_v)
RT_effects_ongoing_fixed_v$effect <- NA
RT_effects_ongoing_fixed_v$effect <- factor(rownames(RT_effects_ongoing_fixed_v))
str(RT_effects_ongoing_fixed_v)
RT_effects_ongoing_fixed_v

ggplot(RT_effects_ongoing_fixed_v, aes(effect, mean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25) +
  geom_point(aes(effect, data), pch = 21, size = 4, colour = "black") +
  xlab("Effect") + ylab("RT (s)") + 
  theme_minimal() +
  ggtitle("Drift rates averaged over all conditions")


# % of data effect explained by averaged mechanism
RT_percent_ongoing_fixed_v <- (RT_effects_ongoing_fixed_v$mean/RT_effects_ongoing_fixed_v$data)*100

RT_percent_ongoing_fixed_v <- data.frame(RT_percent_ongoing_fixed_v)
rownames(RT_percent_ongoing_fixed_v) <- rownames(RT_effects_ongoing_fixed_v)
colnames(RT_percent_ongoing_fixed_v) <- "Effect %"
RT_percent_ongoing_fixed_v

RT_percent_ongoing_fixed_v$`Diff. from 100%` <- NA
RT_percent_ongoing_fixed_v$`Diff. from 100%` <- RT_percent_ongoing_fixed_v$`Effect %` - rep(100,length(RT_percent_ongoing_fixed_v$`Diff. from 100%`))
RT_percent_ongoing_fixed_v$`Effect %` <- round(RT_percent_ongoing_fixed_v$`Effect %`,2)
RT_percent_ongoing_fixed_v$`Diff. from 100%` <- round(RT_percent_ongoing_fixed_v$`Diff. from 100%`,2)
RT_percent_ongoing_fixed_v


# -------------------------------------------------------------------------

# PM task accuracy effects with drift rates averaged over all conditions
RP_effects_PM_fixed_v <- data.frame(RP_effects_PM_fixed_v)
RP_effects_PM_fixed_v$effect <- NA
RP_effects_PM_fixed_v$effect <- factor(rownames(RP_effects_PM_fixed_v))
str(RP_effects_PM_fixed_v)
RP_effects_PM_fixed_v

ggplot(RP_effects_PM_fixed_v, aes(effect, mean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25) +
  geom_point(aes(effect, data), pch = 21, size = 4, colour = "black") +
  xlab("Effect") + ylab("Accuracy") + 
  theme_minimal() +
  ggtitle("Drift rates averaged over all conditions")


# % of data effect explained by averaged mechanism
RP_percent_PM_fixed_v <- (RP_effects_PM_fixed_v$mean/RP_effects_PM_fixed_v$data)*100

RP_percent_PM_fixed_v <- data.frame(RP_percent_PM_fixed_v)
rownames(RP_percent_PM_fixed_v) <- rownames(RP_effects_PM_fixed_v)
colnames(RP_percent_PM_fixed_v) <- "Effect %"
RP_percent_PM_fixed_v

RP_percent_PM_fixed_v$`Diff. from 100%` <- NA
RP_percent_PM_fixed_v$`Diff. from 100%` <- RP_percent_PM_fixed_v$`Effect %` - rep(100,length(RP_percent_PM_fixed_v$`Diff. from 100%`))
RP_percent_PM_fixed_v$`Effect %` <- round(RP_percent_PM_fixed_v$`Effect %`,2)
RP_percent_PM_fixed_v$`Diff. from 100%` <- round(RP_percent_PM_fixed_v$`Diff. from 100%`,2)
RP_percent_PM_fixed_v


# -------------------------------------------------------------------------

# PM task RT effects with drift rates averaged over all conditions
RT_effects_PM_fixed_v <- data.frame(RT_effects_PM_fixed_v)
RT_effects_PM_fixed_v$effect <- NA
RT_effects_PM_fixed_v$effect <- factor(rownames(RT_effects_PM_fixed_v))
str(RT_effects_PM_fixed_v)
RT_effects_PM_fixed_v

ggplot(RT_effects_PM_fixed_v, aes(effect, mean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25) +
  geom_point(aes(effect, data), pch = 21, size = 4, colour = "black") +
  xlab("Effect") + ylab("RT (s)") + 
  theme_minimal() +
  ggtitle("Drift rates averaged over all conditions")


# % of data effect explained by averaged mechanism
RT_percent_PM_fixed_v <- (RT_effects_PM_fixed_v$mean/RT_effects_PM_fixed_v$data)*100

RT_percent_PM_fixed_v <- data.frame(RT_percent_PM_fixed_v)
rownames(RT_percent_PM_fixed_v) <- rownames(RT_effects_PM_fixed_v)
colnames(RT_percent_PM_fixed_v) <- "Effect %"
RT_percent_PM_fixed_v

RT_percent_PM_fixed_v$`Diff. from 100%` <- NA
RT_percent_PM_fixed_v$`Diff. from 100%` <- RT_percent_PM_fixed_v$`Effect %` - rep(100,length(RT_percent_PM_fixed_v$`Diff. from 100%`))
RT_percent_PM_fixed_v$`Effect %` <- round(RT_percent_PM_fixed_v$`Effect %`,2)
RT_percent_PM_fixed_v$`Diff. from 100%` <- round(RT_percent_PM_fixed_v$`Diff. from 100%`,2)
RT_percent_PM_fixed_v


# -------------------------------------------------------------------------

# -------------------------------------------------------------------------

# Ongoing task accuracy effects with thresholds averaged over all conditions
RP_effects_ongoing_fixed_B <- data.frame(RP_effects_ongoing_fixed_B)
RP_effects_ongoing_fixed_B$effect <- NA
RP_effects_ongoing_fixed_B$effect <- factor(rownames(RP_effects_ongoing_fixed_B))
str(RP_effects_ongoing_fixed_B)
RP_effects_ongoing_fixed_B

ggplot(RP_effects_ongoing_fixed_B, aes(effect, mean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25) +
  geom_point(aes(effect, data), pch = 21, size = 4, colour = "black") +
  xlab("Effect") + ylab("Accuracy") + 
  theme_minimal() +
  ggtitle("Thresholds averaged over all conditions")


# % of data effect explained by averaged mechanism
RP_percent_ongoing_fixed_B <- (RP_effects_ongoing_fixed_B$mean/RP_effects_ongoing_fixed_B$data)*100

RP_percent_ongoing_fixed_B <- data.frame(RP_percent_ongoing_fixed_B)
rownames(RP_percent_ongoing_fixed_B) <- rownames(RP_effects_ongoing_fixed_B)
colnames(RP_percent_ongoing_fixed_B) <- "Effect %"
RP_percent_ongoing_fixed_B

RP_percent_ongoing_fixed_B$`Diff. from 100%` <- NA
RP_percent_ongoing_fixed_B$`Diff. from 100%` <- RP_percent_ongoing_fixed_B$`Effect %` - rep(100,length(RP_percent_ongoing_fixed_B$`Diff. from 100%`))
RP_percent_ongoing_fixed_B$`Effect %` <- round(RP_percent_ongoing_fixed_B$`Effect %`,2)
RP_percent_ongoing_fixed_B$`Diff. from 100%` <- round(RP_percent_ongoing_fixed_B$`Diff. from 100%`,2)
RP_percent_ongoing_fixed_B


# -------------------------------------------------------------------------

# Ongoing task RT effects with thresholds averaged over all conditions
RT_effects_ongoing_fixed_B <- data.frame(RT_effects_ongoing_fixed_B)
RT_effects_ongoing_fixed_B$effect <- NA
RT_effects_ongoing_fixed_B$effect <- factor(rownames(RT_effects_ongoing_fixed_B))
str(RT_effects_ongoing_fixed_B)
RT_effects_ongoing_fixed_B

ggplot(RT_effects_ongoing_fixed_B, aes(effect, mean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25) +
  geom_point(aes(effect, data), pch = 21, size = 4, colour = "black") +
  xlab("Effect") + ylab("RT (s)") + 
  theme_minimal() +
  ggtitle("Thresholds averaged over all conditions")


# % of data effect explained by averaged mechanism
RT_percent_ongoing_fixed_B <- (RT_effects_ongoing_fixed_B$mean/RT_effects_ongoing_fixed_B$data)*100

RT_percent_ongoing_fixed_B <- data.frame(RT_percent_ongoing_fixed_B)
rownames(RT_percent_ongoing_fixed_B) <- rownames(RT_effects_ongoing_fixed_B)
colnames(RT_percent_ongoing_fixed_B) <- "Effect %"
RT_percent_ongoing_fixed_B

RT_percent_ongoing_fixed_B$`Diff. from 100%` <- NA
RT_percent_ongoing_fixed_B$`Diff. from 100%` <- RT_percent_ongoing_fixed_B$`Effect %` - rep(100,length(RT_percent_ongoing_fixed_B$`Diff. from 100%`))
RT_percent_ongoing_fixed_B$`Effect %` <- round(RT_percent_ongoing_fixed_B$`Effect %`,2)
RT_percent_ongoing_fixed_B$`Diff. from 100%` <- round(RT_percent_ongoing_fixed_B$`Diff. from 100%`,2)
RT_percent_ongoing_fixed_B


# -------------------------------------------------------------------------

# PM task accuracy effects with thresholds averaged over all conditions
RP_effects_PM_fixed_B <- data.frame(RP_effects_PM_fixed_B)
RP_effects_PM_fixed_B$effect <- NA
RP_effects_PM_fixed_B$effect <- factor(rownames(RP_effects_PM_fixed_B))
str(RP_effects_PM_fixed_B)
RP_effects_PM_fixed_B

ggplot(RP_effects_PM_fixed_B, aes(effect, mean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25) +
  geom_point(aes(effect, data), pch = 21, size = 4, colour = "black") +
  xlab("Effect") + ylab("Accuracy") + 
  theme_minimal() +
  ggtitle("Thresholds averaged over all conditions")


# % of data effect explained by averaged mechanism
RP_percent_PM_fixed_B <- (RP_effects_PM_fixed_B$mean/RP_effects_PM_fixed_B$data)*100

RP_percent_PM_fixed_B <- data.frame(RP_percent_PM_fixed_B)
rownames(RP_percent_PM_fixed_B) <- rownames(RP_effects_PM_fixed_B)
colnames(RP_percent_PM_fixed_B) <- "Effect %"
RP_percent_PM_fixed_B

RP_percent_PM_fixed_B$`Diff. from 100%` <- NA
RP_percent_PM_fixed_B$`Diff. from 100%` <- RP_percent_PM_fixed_B$`Effect %` - rep(100,length(RP_percent_PM_fixed_B$`Diff. from 100%`))
RP_percent_PM_fixed_B$`Effect %` <- round(RP_percent_PM_fixed_B$`Effect %`,2)
RP_percent_PM_fixed_B$`Diff. from 100%` <- round(RP_percent_PM_fixed_B$`Diff. from 100%`,2)
RP_percent_PM_fixed_B


# -------------------------------------------------------------------------

# PM task RT effects with thresholds averaged over all conditions
RT_effects_PM_fixed_B <- data.frame(RT_effects_PM_fixed_B)
RT_effects_PM_fixed_B$effect <- NA
RT_effects_PM_fixed_B$effect <- factor(rownames(RT_effects_PM_fixed_B))
str(RT_effects_PM_fixed_B)
RT_effects_PM_fixed_B

ggplot(RT_effects_PM_fixed_B, aes(effect, mean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25) +
  geom_point(aes(effect, data), pch = 21, size = 4, colour = "black") +
  xlab("Effect") + ylab("RT (s)") + 
  theme_minimal() +
  ggtitle("Thresholds averaged over all conditions")


# % of data effect explained by averaged mechanism
RT_percent_PM_fixed_B <- (RT_effects_PM_fixed_B$mean/RT_effects_PM_fixed_B$data)*100

RT_percent_PM_fixed_B <- data.frame(RT_percent_PM_fixed_B)
rownames(RT_percent_PM_fixed_B) <- rownames(RT_effects_PM_fixed_B)
colnames(RT_percent_PM_fixed_B) <- "Effect %"
RT_percent_PM_fixed_B

RT_percent_PM_fixed_B$`Diff. from 100%` <- NA
RT_percent_PM_fixed_B$`Diff. from 100%` <- RT_percent_PM_fixed_B$`Effect %` - rep(100,length(RT_percent_PM_fixed_B$`Diff. from 100%`))
RT_percent_PM_fixed_B$`Effect %` <- round(RT_percent_PM_fixed_B$`Effect %`,2)
RT_percent_PM_fixed_B$`Diff. from 100%` <- round(RT_percent_PM_fixed_B$`Diff. from 100%`,2)
RT_percent_PM_fixed_B



# -------------------------------------------------------------------------

# -------------------------------------------------------------------------

# Ongoing task accuracy effects with reactive control turned off

RP_effects_ongoing_no_reactive <- data.frame(RP_effects_ongoing_no_reactive)
RP_effects_ongoing_no_reactive$effect <- NA
RP_effects_ongoing_no_reactive$effect <- factor(rownames(RP_effects_ongoing_no_reactive))
str(RP_effects_ongoing_no_reactive)
RP_effects_ongoing_no_reactive

ggplot(RP_effects_ongoing_no_reactive, aes(effect, mean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25) +
  geom_point(aes(effect, data), pch = 21, size = 4, colour = "black") +
  xlab("Effect") + ylab("Accuracy") + 
  theme_minimal() +
  ggtitle("No reactive control")

# % of data effect explained by removed mechanism
RP_percent_ongoing_no_reactive <- (RP_effects_ongoing_no_reactive$mean/RP_effects_ongoing_no_reactive$data)*100

RP_percent_ongoing_no_reactive <- data.frame(RP_percent_ongoing_no_reactive)
rownames(RP_percent_ongoing_no_reactive) <- rownames(RP_effects_ongoing_no_reactive)
colnames(RP_percent_ongoing_no_reactive) <- "Effect %"
RP_percent_ongoing_no_reactive

RP_percent_ongoing_no_reactive$`Diff. from 100%` <- NA
RP_percent_ongoing_no_reactive$`Diff. from 100%` <- RP_percent_ongoing_no_reactive$`Effect %` - 
  rep(100,length(RP_percent_ongoing_no_reactive))
RP_percent_ongoing_no_reactive$`Effect %` <- round(RP_percent_ongoing_no_reactive$`Effect %`,2)
RP_percent_ongoing_no_reactive$`Diff. from 100%` <- round(RP_percent_ongoing_no_reactive$`Diff. from 100%`,2)
RP_percent_ongoing_no_reactive


# -------------------------------------------------------------------------

# Ongoing task RT effects with reactive control turned off

RT_effects_ongoing_no_reactive <- data.frame(RT_effects_ongoing_no_reactive)
RT_effects_ongoing_no_reactive$effect <- NA
RT_effects_ongoing_no_reactive$effect <- factor(rownames(RT_effects_ongoing_no_reactive))
str(RT_effects_ongoing_no_reactive)
RT_effects_ongoing_no_reactive

ggplot(RT_effects_ongoing_no_reactive, aes(effect, mean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25) +
  geom_point(aes(effect, data), pch = 21, size = 4, colour = "black") +
  xlab("Effect") + ylab("RT (s)") + 
  theme_minimal() +
  ggtitle("No reactive control")

# % of data effect explained by removed mechanism
RT_percent_ongoing_no_reactive <- (RT_effects_ongoing_no_reactive$mean/RT_effects_ongoing_no_reactive$data)*100

RT_percent_ongoing_no_reactive <- data.frame(RT_percent_ongoing_no_reactive)
rownames(RT_percent_ongoing_no_reactive) <- rownames(RT_effects_ongoing_no_reactive)
colnames(RT_percent_ongoing_no_reactive) <- "Effect %"
RT_percent_ongoing_no_reactive

RT_percent_ongoing_no_reactive$`Diff. from 100%` <- NA
RT_percent_ongoing_no_reactive$`Diff. from 100%` <- RT_percent_ongoing_no_reactive$`Effect %` - 
  rep(100,length(RT_percent_ongoing_no_reactive))
RT_percent_ongoing_no_reactive$`Effect %` <- round(RT_percent_ongoing_no_reactive$`Effect %`,2)
RT_percent_ongoing_no_reactive$`Diff. from 100%` <- round(RT_percent_ongoing_no_reactive$`Diff. from 100%`,2)
RT_percent_ongoing_no_reactive


# -------------------------------------------------------------------------

# PM task accuracy effects with reactive control turned off

RP_effects_PM_no_reactive <- data.frame(RP_effects_PM_no_reactive)
RP_effects_PM_no_reactive$effect <- NA
RP_effects_PM_no_reactive$effect <- factor(rownames(RP_effects_PM_no_reactive))
str(RP_effects_PM_no_reactive)
RP_effects_PM_no_reactive

ggplot(RP_effects_PM_no_reactive, aes(effect, mean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25) +
  geom_point(aes(effect, data), pch = 21, size = 4, colour = "black") +
  xlab("Effect") + ylab("Accuracy") + 
  theme_minimal() +
  ggtitle("No reactive control")

# % of data effect explained by removed mechanism
RP_percent_PM_no_reactive <- (RP_effects_PM_no_reactive$mean/RP_effects_PM_no_reactive$data)*100

RP_percent_PM_no_reactive <- data.frame(RP_percent_PM_no_reactive)
rownames(RP_percent_PM_no_reactive) <- rownames(RP_effects_PM_no_reactive)
colnames(RP_percent_PM_no_reactive) <- "Effect %"
RP_percent_PM_no_reactive

RP_percent_PM_no_reactive$`Diff. from 100%` <- NA
RP_percent_PM_no_reactive$`Diff. from 100%` <- RP_percent_PM_no_reactive$`Effect %` - 
  rep(100,length(RP_percent_PM_no_reactive))
RP_percent_PM_no_reactive$`Effect %` <- round(RP_percent_PM_no_reactive$`Effect %`,2)
RP_percent_PM_no_reactive$`Diff. from 100%` <- round(RP_percent_PM_no_reactive$`Diff. from 100%`,2)
RP_percent_PM_no_reactive


# -------------------------------------------------------------------------

# PM task RT effects with reactive control turned off

RT_effects_PM_no_reactive <- data.frame(RT_effects_PM_no_reactive)
RT_effects_PM_no_reactive$effect <- NA
RT_effects_PM_no_reactive$effect <- factor(rownames(RT_effects_PM_no_reactive))
str(RT_effects_PM_no_reactive)
RT_effects_PM_no_reactive

ggplot(RT_effects_PM_no_reactive, aes(effect, mean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25) +
  geom_point(aes(effect, data), pch = 21, size = 4, colour = "black") +
  xlab("Effect") + ylab("RT (s)") + 
  theme_minimal() +
  ggtitle("No reactive control")

# % of data effect explained by removed mechanism
RT_percent_PM_no_reactive <- (RT_effects_PM_no_reactive$mean/RT_effects_PM_no_reactive$data)*100

RT_percent_PM_no_reactive <- data.frame(RT_percent_PM_no_reactive)
rownames(RT_percent_PM_no_reactive) <- rownames(RT_effects_PM_no_reactive)
colnames(RT_percent_PM_no_reactive) <- "Effect %"
RT_percent_PM_no_reactive

RT_percent_PM_no_reactive$`Diff. from 100%` <- NA
RT_percent_PM_no_reactive$`Diff. from 100%` <- RT_percent_PM_no_reactive$`Effect %` - 
  rep(100,length(RT_percent_PM_no_reactive))
RT_percent_PM_no_reactive$`Effect %` <- round(RT_percent_PM_no_reactive$`Effect %`,2)
RT_percent_PM_no_reactive$`Diff. from 100%` <- round(RT_percent_PM_no_reactive$`Diff. from 100%`,2)
RT_percent_PM_no_reactive



# -------------------------------------------------------------------------

# -------------------------------------------------------------------------

# Ongoing task accuracy effects with proactive control turned off across TP levels

RP_effects_ongoing_no_proactive_TP <- data.frame(RP_effects_ongoing_no_proactive_TP)
RP_effects_ongoing_no_proactive_TP$effect <- NA
RP_effects_ongoing_no_proactive_TP$effect <- factor(rownames(RP_effects_ongoing_no_proactive_TP))
str(RP_effects_ongoing_no_proactive_TP)
RP_effects_ongoing_no_proactive_TP

ggplot(RP_effects_ongoing_no_proactive_TP, aes(effect, mean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25) +
  geom_point(aes(effect, data), pch = 21, size = 4, colour = "black") +
  xlab("Effect") + ylab("Accuracy") + 
  theme_minimal() +
  ggtitle("No proactive control by time pressure")

# % of data effect explained by removed mechanism
RP_percent_ongoing_no_proactive_TP <- (RP_effects_ongoing_no_proactive_TP$mean/RP_effects_ongoing_no_proactive_TP$data)*100

RP_percent_ongoing_no_proactive_TP <- data.frame(RP_percent_ongoing_no_proactive_TP)
rownames(RP_percent_ongoing_no_proactive_TP) <- rownames(RP_effects_ongoing_no_proactive_TP)
colnames(RP_percent_ongoing_no_proactive_TP) <- "Effect %"
RP_percent_ongoing_no_proactive_TP

RP_percent_ongoing_no_proactive_TP$`Diff. from 100%` <- NA
RP_percent_ongoing_no_proactive_TP$`Diff. from 100%` <- RP_percent_ongoing_no_proactive_TP$`Effect %` - 
  rep(100,length(RP_percent_ongoing_no_proactive_TP))
RP_percent_ongoing_no_proactive_TP$`Effect %` <- round(RP_percent_ongoing_no_proactive_TP$`Effect %`,2)
RP_percent_ongoing_no_proactive_TP$`Diff. from 100%` <- round(RP_percent_ongoing_no_proactive_TP$`Diff. from 100%`,2)
RP_percent_ongoing_no_proactive_TP


# -------------------------------------------------------------------------

# Ongoing task RT effects with proactive control turned off across TP levels

RT_effects_ongoing_no_proactive_TP <- data.frame(RT_effects_ongoing_no_proactive_TP)
RT_effects_ongoing_no_proactive_TP$effect <- NA
RT_effects_ongoing_no_proactive_TP$effect <- factor(rownames(RT_effects_ongoing_no_proactive_TP))
str(RT_effects_ongoing_no_proactive_TP)
RT_effects_ongoing_no_proactive_TP

ggplot(RT_effects_ongoing_no_proactive_TP, aes(effect, mean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25) +
  geom_point(aes(effect, data), pch = 21, size = 4, colour = "black") +
  xlab("Effect") + ylab("RT (s)") + 
  theme_minimal() +
  ggtitle("No proactive control by time pressure")

# % of data effect explained by removed mechanism
RT_percent_ongoing_no_proactive_TP <- (RT_effects_ongoing_no_proactive_TP$mean/RT_effects_ongoing_no_proactive_TP$data)*100

RT_percent_ongoing_no_proactive_TP <- data.frame(RT_percent_ongoing_no_proactive_TP)
rownames(RT_percent_ongoing_no_proactive_TP) <- rownames(RT_effects_ongoing_no_proactive_TP)
colnames(RT_percent_ongoing_no_proactive_TP) <- "Effect %"
RT_percent_ongoing_no_proactive_TP

RT_percent_ongoing_no_proactive_TP$`Diff. from 100%` <- NA
RT_percent_ongoing_no_proactive_TP$`Diff. from 100%` <- RT_percent_ongoing_no_proactive_TP$`Effect %` - 
  rep(100,length(RT_percent_ongoing_no_proactive_TP))
RT_percent_ongoing_no_proactive_TP$`Effect %` <- round(RT_percent_ongoing_no_proactive_TP$`Effect %`,2)
RT_percent_ongoing_no_proactive_TP$`Diff. from 100%` <- round(RT_percent_ongoing_no_proactive_TP$`Diff. from 100%`,2)
RT_percent_ongoing_no_proactive_TP


# -------------------------------------------------------------------------

# PM task accuracy effects with proactive control turned off across TP levels

RP_effects_PM_no_proactive_TP <- data.frame(RP_effects_PM_no_proactive_TP)
RP_effects_PM_no_proactive_TP$effect <- NA
RP_effects_PM_no_proactive_TP$effect <- factor(rownames(RP_effects_PM_no_proactive_TP))
str(RP_effects_PM_no_proactive_TP)
RP_effects_PM_no_proactive_TP

ggplot(RP_effects_PM_no_proactive_TP, aes(effect, mean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25) +
  geom_point(aes(effect, data), pch = 21, size = 4, colour = "black") +
  xlab("Effect") + ylab("Accuracy") + 
  theme_minimal() +
  ggtitle("No proactive control by time pressure")

# % of data effect explained by removed mechanism
RP_percent_PM_no_proactive_TP <- (RP_effects_PM_no_proactive_TP$mean/RP_effects_PM_no_proactive_TP$data)*100

RP_percent_PM_no_proactive_TP <- data.frame(RP_percent_PM_no_proactive_TP)
rownames(RP_percent_PM_no_proactive_TP) <- rownames(RP_effects_PM_no_proactive_TP)
colnames(RP_percent_PM_no_proactive_TP) <- "Effect %"
RP_percent_PM_no_proactive_TP

RP_percent_PM_no_proactive_TP$`Diff. from 100%` <- NA
RP_percent_PM_no_proactive_TP$`Diff. from 100%` <- RP_percent_PM_no_proactive_TP$`Effect %` - 
  rep(100,length(RP_percent_PM_no_proactive_TP))
RP_percent_PM_no_proactive_TP$`Effect %` <- round(RP_percent_PM_no_proactive_TP$`Effect %`,2)
RP_percent_PM_no_proactive_TP$`Diff. from 100%` <- round(RP_percent_PM_no_proactive_TP$`Diff. from 100%`,2)
RP_percent_PM_no_proactive_TP


# -------------------------------------------------------------------------

# PM task RT effects with proactive control turned off across TP levels

RT_effects_PM_no_proactive_TP <- data.frame(RT_effects_PM_no_proactive_TP)
RT_effects_PM_no_proactive_TP$effect <- NA
RT_effects_PM_no_proactive_TP$effect <- factor(rownames(RT_effects_PM_no_proactive_TP))
str(RT_effects_PM_no_proactive_TP)
RT_effects_PM_no_proactive_TP

ggplot(RT_effects_PM_no_proactive_TP, aes(effect, mean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25) +
  geom_point(aes(effect, data), pch = 21, size = 4, colour = "black") +
  xlab("Effect") + ylab("RT (s)") + 
  theme_minimal() +
  ggtitle("No proactive control by time pressure")

# % of data effect explained by removed mechanism
RT_percent_PM_no_proactive_TP <- (RT_effects_PM_no_proactive_TP$mean/RT_effects_PM_no_proactive_TP$data)*100

RT_percent_PM_no_proactive_TP <- data.frame(RT_percent_PM_no_proactive_TP)
rownames(RT_percent_PM_no_proactive_TP) <- rownames(RT_effects_PM_no_proactive_TP)
colnames(RT_percent_PM_no_proactive_TP) <- "Effect %"
RT_percent_PM_no_proactive_TP

RT_percent_PM_no_proactive_TP$`Diff. from 100%` <- NA
RT_percent_PM_no_proactive_TP$`Diff. from 100%` <- RT_percent_PM_no_proactive_TP$`Effect %` - 
  rep(100,length(RT_percent_PM_no_proactive_TP))
RT_percent_PM_no_proactive_TP$`Effect %` <- round(RT_percent_PM_no_proactive_TP$`Effect %`,2)
RT_percent_PM_no_proactive_TP$`Diff. from 100%` <- round(RT_percent_PM_no_proactive_TP$`Diff. from 100%`,2)
RT_percent_PM_no_proactive_TP


# -------------------------------------------------------------------------

# -------------------------------------------------------------------------

# Ongoing task accuracy effects with proactive control turned off across PM levels

RP_effects_ongoing_no_proactive_PM <- data.frame(RP_effects_ongoing_no_proactive_PM)
RP_effects_ongoing_no_proactive_PM$effect <- NA
RP_effects_ongoing_no_proactive_PM$effect <- factor(rownames(RP_effects_ongoing_no_proactive_PM))
str(RP_effects_ongoing_no_proactive_PM)
RP_effects_ongoing_no_proactive_PM

ggplot(RP_effects_ongoing_no_proactive_PM, aes(effect, mean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25) +
  geom_point(aes(effect, data), pch = 21, size = 4, colour = "black") +
  xlab("Effect") + ylab("Accuracy") + 
  theme_minimal() +
  ggtitle("No proactive control by PM frequency")

# % of data effect explained by removed mechanism
RP_percent_ongoing_no_proactive_PM <- (RP_effects_ongoing_no_proactive_PM$mean/RP_effects_ongoing_no_proactive_PM$data)*100

RP_percent_ongoing_no_proactive_PM <- data.frame(RP_percent_ongoing_no_proactive_PM)
rownames(RP_percent_ongoing_no_proactive_PM) <- rownames(RP_effects_ongoing_no_proactive_PM)
colnames(RP_percent_ongoing_no_proactive_PM) <- "Effect %"
RP_percent_ongoing_no_proactive_PM

RP_percent_ongoing_no_proactive_PM$`Diff. from 100%` <- NA
RP_percent_ongoing_no_proactive_PM$`Diff. from 100%` <- RP_percent_ongoing_no_proactive_PM$`Effect %` - 
  rep(100,length(RP_percent_ongoing_no_proactive_PM))
RP_percent_ongoing_no_proactive_PM$`Effect %` <- round(RP_percent_ongoing_no_proactive_PM$`Effect %`,2)
RP_percent_ongoing_no_proactive_PM$`Diff. from 100%` <- round(RP_percent_ongoing_no_proactive_PM$`Diff. from 100%`,2)
RP_percent_ongoing_no_proactive_PM


# -------------------------------------------------------------------------

# Ongoing task RT effects with proactive control turned off across PM levels

RT_effects_ongoing_no_proactive_PM <- data.frame(RT_effects_ongoing_no_proactive_PM)
RT_effects_ongoing_no_proactive_PM$effect <- NA
RT_effects_ongoing_no_proactive_PM$effect <- factor(rownames(RT_effects_ongoing_no_proactive_PM))
str(RT_effects_ongoing_no_proactive_PM)
RT_effects_ongoing_no_proactive_PM

ggplot(RT_effects_ongoing_no_proactive_PM, aes(effect, mean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25) +
  geom_point(aes(effect, data), pch = 21, size = 4, colour = "black") +
  xlab("Effect") + ylab("RT (s)") + 
  theme_minimal() +
  ggtitle("No proactive control by PM frequency")

# % of data effect explained by removed mechanism
RT_percent_ongoing_no_proactive_PM <- (RT_effects_ongoing_no_proactive_PM$mean/RT_effects_ongoing_no_proactive_PM$data)*100

RT_percent_ongoing_no_proactive_PM <- data.frame(RT_percent_ongoing_no_proactive_PM)
rownames(RT_percent_ongoing_no_proactive_PM) <- rownames(RT_effects_ongoing_no_proactive_PM)
colnames(RT_percent_ongoing_no_proactive_PM) <- "Effect %"
RT_percent_ongoing_no_proactive_PM

RT_percent_ongoing_no_proactive_PM$`Diff. from 100%` <- NA
RT_percent_ongoing_no_proactive_PM$`Diff. from 100%` <- RT_percent_ongoing_no_proactive_PM$`Effect %` - 
  rep(100,length(RT_percent_ongoing_no_proactive_PM))
RT_percent_ongoing_no_proactive_PM$`Effect %` <- round(RT_percent_ongoing_no_proactive_PM$`Effect %`,2)
RT_percent_ongoing_no_proactive_PM$`Diff. from 100%` <- round(RT_percent_ongoing_no_proactive_PM$`Diff. from 100%`,2)
RT_percent_ongoing_no_proactive_PM


# -------------------------------------------------------------------------

# PM task accuracy effects with proactive control turned off across PM levels

RP_effects_PM_no_proactive_PM <- data.frame(RP_effects_PM_no_proactive_PM)
RP_effects_PM_no_proactive_PM$effect <- NA
RP_effects_PM_no_proactive_PM$effect <- factor(rownames(RP_effects_PM_no_proactive_PM))
str(RP_effects_PM_no_proactive_PM)
RP_effects_PM_no_proactive_PM

ggplot(RP_effects_PM_no_proactive_PM, aes(effect, mean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25) +
  geom_point(aes(effect, data), pch = 21, size = 4, colour = "black") +
  xlab("Effect") + ylab("Accuracy") + 
  theme_minimal() +
  ggtitle("No proactive control by PM frequency")

# % of data effect explained by removed mechanism
RP_percent_PM_no_proactive_PM <- (RP_effects_PM_no_proactive_PM$mean/RP_effects_PM_no_proactive_PM$data)*100

RP_percent_PM_no_proactive_PM <- data.frame(RP_percent_PM_no_proactive_PM)
rownames(RP_percent_PM_no_proactive_PM) <- rownames(RP_effects_PM_no_proactive_PM)
colnames(RP_percent_PM_no_proactive_PM) <- "Effect %"
RP_percent_PM_no_proactive_PM

RP_percent_PM_no_proactive_PM$`Diff. from 100%` <- NA
RP_percent_PM_no_proactive_PM$`Diff. from 100%` <- RP_percent_PM_no_proactive_PM$`Effect %` - 
  rep(100,length(RP_percent_PM_no_proactive_PM))
RP_percent_PM_no_proactive_PM$`Effect %` <- round(RP_percent_PM_no_proactive_PM$`Effect %`,2)
RP_percent_PM_no_proactive_PM$`Diff. from 100%` <- round(RP_percent_PM_no_proactive_PM$`Diff. from 100%`,2)
RP_percent_PM_no_proactive_PM


# -------------------------------------------------------------------------

# PM task RT effects with proactive control turned off across PM levels

RT_effects_PM_no_proactive_PM <- data.frame(RT_effects_PM_no_proactive_PM)
RT_effects_PM_no_proactive_PM$effect <- NA
RT_effects_PM_no_proactive_PM$effect <- factor(rownames(RT_effects_PM_no_proactive_PM))
str(RT_effects_PM_no_proactive_PM)
RT_effects_PM_no_proactive_PM

ggplot(RT_effects_PM_no_proactive_PM, aes(effect, mean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25) +
  geom_point(aes(effect, data), pch = 21, size = 4, colour = "black") +
  xlab("Effect") + ylab("RT (s)") + 
  theme_minimal() +
  ggtitle("No proactive control by PM frequency")

# % of data effect explained by removed mechanism
RT_percent_PM_no_proactive_PM <- (RT_effects_PM_no_proactive_PM$mean/RT_effects_PM_no_proactive_PM$data)*100

RT_percent_PM_no_proactive_PM <- data.frame(RT_percent_PM_no_proactive_PM)
rownames(RT_percent_PM_no_proactive_PM) <- rownames(RT_effects_PM_no_proactive_PM)
colnames(RT_percent_PM_no_proactive_PM) <- "Effect %"
RT_percent_PM_no_proactive_PM

RT_percent_PM_no_proactive_PM$`Diff. from 100%` <- NA
RT_percent_PM_no_proactive_PM$`Diff. from 100%` <- RT_percent_PM_no_proactive_PM$`Effect %` - 
  rep(100,length(RT_percent_PM_no_proactive_PM))
RT_percent_PM_no_proactive_PM$`Effect %` <- round(RT_percent_PM_no_proactive_PM$`Effect %`,2)
RT_percent_PM_no_proactive_PM$`Diff. from 100%` <- round(RT_percent_PM_no_proactive_PM$`Diff. from 100%`,2)
RT_percent_PM_no_proactive_PM

