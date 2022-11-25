# Clear workspace
rm(list=ls())

# Set working directory
getwd()

# Load packages
library(tidyverse)
library(gridExtra)

# Source model functions
source("dmc/dmc.R")
load_model("LBA", "lba_B.R")


# Summary functions -------------------------------------------------------

fixedeffects.meanthetas <- function(samples){
  ## bring longer thetas down to min nmc by sampling
  nmcs <- sapply(samples, function(x) x$nmc)
  nmc <- min(nmcs)
  for (i in 1:length(samples)) if (nmcs[i] > nmc) samples[[i]]$theta <-
    samples[[i]]$theta[,,sample(1:dim(samples[[i]]$theta)[3], nmc)]
  samps <- lapply(samples, function(x) x["theta"])
  
  ## thetas into big array for apply
  samps2 <- unlist(samps)
  dim3 <- c(dim(samps[[1]]$theta), length(samps2)/prod(dim(samps[[1]]$theta)))
  dim(samps2) <- dim3
  samps3 <- apply(samps2, c(1,2,3), mean)
  
  ## back to a theta list after applied
  colnames(samps3) <- colnames(samps[[1]]$theta)
  samps5 <- list(samps3)
  attributes(samps5) <- attributes(samps[[1]])
  samps5
}


# -------------------------------------------------------------------------

# Load samples
print(load("samples/sTPPM_full_sdvS.RData"))
samples <- samples1
# print(load("samples/sTPPM_full.RData"))
# samples <- samples2
# print(load("samples/sTPPM_B.RData"))
# samples <- samples1
samples[[1]]$p.names


# Get parameter summary ---------------------------------------------------
# summary.dmc()
# summary.dmc(samples[[1]])$statistics

# parms <- summary.dmc(samples)
# parms <- do.call(rbind, lapply(parms, function(x) x$statistics[,1]))
# parms

# Load parameter summary
print(load("deriv/map_parms_full_sdvS.RData"))
# print(load("deriv/map_parms_full.RData"))
str(parms)
head(parms)
nrow(parms)

# Mean over participants
colMeans(parms)

# Summarize thetas
mean_thetas <- fixedeffects.meanthetas(samples)[[1]]

# Save
save(mean_thetas, file = "deriv/mean_thetas_full_sdvS.RData")
# save(mean_thetas, file = "deriv/mean_thetas_full.RData")

# Load
print(load("deriv/mean_thetas_full_sdvS.RData"))
# print(load("deriv/mean_thetas_full.RData"))

# Explore mean thetas
str(mean_thetas)
dim(mean_thetas)
mean_thetas

# Create parameter data frame ---------------------------------------------

msds <- cbind(apply(mean_thetas, 2, median), apply(mean_thetas, 2, sd))
colnames(msds) <- c("M", "SD")
msds <- data.frame(msds)
msds
colMeans(parms)


# Add factors for plotting ------------------------------------------------
ps <- data.frame(msds)
ps$TP <- NA; ps$PM <- NA; ps$S <- NA; ps$R <- NA

ps$TP[grep("3s", rownames(ps))] <- "3s"
ps$TP[grep("6s", rownames(ps))] <- "6s"

ps$PM[grep("10", rownames(ps))] <- "10% PM"
ps$PM[grep("30", rownames(ps))] <- "30% PM"

ps$S[grep("cc", rownames(ps))] <- "Conflict"
ps$S[grep("nn", rownames(ps))] <- "Non-conflict"
ps$S[grep("pc", rownames(ps))] <- "PM conflict"
ps$S[grep("pn", rownames(ps))] <- "PM non-conflict"
ps$S[grep("pp", rownames(ps))] <- "PM"

ps$R[grep("C", rownames(ps))] <- "Conflict"
ps$R[grep("N", rownames(ps))] <- "Non-conflict"
ps$R[grep("P", rownames(ps))] <- "PM"

ps$TP <- factor(ps$TP)
ps$PM <- factor(ps$PM)
ps$S <- factor(ps$S)
ps$R <- factor(ps$R)
str(ps)
ps

# Get A
A <- ps[ grep("^A", rownames(ps)), c("M", "SD") ]
A

# Get B
B <- ps[ grep("B.", rownames(ps)), c("M", "SD", "TP", "PM", "R") ]
B

# Get v
v <- ps[ grep("mean_v.", rownames(ps)), ]
v <- v[ -grep("PMFA", rownames(v)), ]
v

# Exclude PM miss rates
v <- v[!(v$S == "PM conflict" & v$R != "PM") & !(v$S == "PM non-conflict" & v$R != "PM"),]
v

# Get t0
t0 <- ps[ grep("t0",rownames(ps)), c("M", "SD") ]
t0


# Make plots --------------------------------------------------------------

# Plot thresholds
B_plot <- B %>%
  ggplot(aes(x = factor(TP), y = M)) +
  geom_point(stat = "identity", aes(color = R, shape = R), size = 3) +
  geom_line(aes(y = M, group = R, color = R), 
            linetype = "dashed", 
            size = 0.8) +
  geom_errorbar(aes(ymin = M - SD, 
                    ymax = M + SD, 
                    width = 0.3, color = R)) +
  ylim(1.1, 2.7) +
  facet_grid(. ~ PM) +
  labs(title = "Threshold", 
       x = "Time pressure", 
       y = "B", 
       color = "Response",
       shape = "Response") +
  theme_minimal()
B_plot

ggsave("plots/B_plot.png", plot = B_plot, 
       width = 2000, height = 1400, units = "px")


# Plot rates
v_plot <- v %>%
  ggplot(aes(x = factor(TP), y = M, shape = R, color = R)) +
  geom_point(stat = "identity", aes(), size = 3) +
  geom_line(aes(y = M, group = R), 
            linetype = "dashed", 
            size = 0.8) +
  geom_errorbar(aes(ymin = M - SD, 
                    ymax = M + SD, 
                    width = 0.3)) +
  ylim(0.5, 2.1) +
  facet_grid(S ~ PM) +
  labs(title = "Accumulation rate", 
       x = "Time pressure", 
       y = "v", 
       color = "Response",
       shape = "Response") +
  theme_minimal()
v_plot

ggsave("plots/v_plot.png", plot = v_plot, 
       width = 2000, height = 1400, units = "px")


# Rates arranged slightly differently
v_plot2 <- v %>%
  ggplot(aes(x = factor(PM), y = M, shape = R, color = R)) +
  geom_point(stat = "identity", aes(), size = 3) +
  geom_line(aes(y = M, group = R), 
            linetype = "dashed", 
            size = 0.8) +
  geom_errorbar(aes(ymin = M - SD, 
                    ymax = M + SD, 
                    width = 0.3)) +
  ylim(0, 2.8) +
  facet_grid(S ~ TP) +
  labs(title = "Accumulation rate", 
       x = "Time pressure", 
       y = "v", 
       color = "Response",
       shape = "Response") +
  theme_minimal()
v_plot2

