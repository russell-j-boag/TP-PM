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

group.inference.dist <- function (hsamples, fun) {
  # Bring longer thetas down to minimum nmc by sampling
  nmcs <- sapply(hsamples, function(x) x$nmc)
  nmc  <- min(nmcs)
  for (i in 1:length(hsamples)) if (nmcs[i] > nmc) {
    hsamples[[i]]$theta <- 
      hsamples[[i]]$theta[, , sample(1:dim(hsamples[[i]]$theta)[3], nmc)]
  }
  # inference <- lapply(hsamples, function(x) x["theta"])
  inference <- list()
  for (i in 1:length(hsamples)) {
    thetas <- hsamples[[i]]$theta
    inference [[i]] <- fun (thetas)
  }
  inf2 <- unlist(inference)
  dim3 <- c(dim(inference[[1]]), length(inf2)/prod(dim(inference[[1]])))
  dim(inf2) <- dim3
  apply(inf2, c(3, 1), mean)
}

subject.inference.dist <- function (hsamples, fun) {
  # Bring longer thetas down to minimum nmc by sampling
  nmcs <- sapply(hsamples, function(x) x$nmc)
  nmc  <- min(nmcs)
  for (i in 1:length(hsamples)) if (nmcs[i] > nmc) {
    hsamples[[i]]$theta <- 
      hsamples[[i]]$theta[, , sample(1:dim(hsamples[[i]]$theta)[3], nmc)]
  }
  # inference <- lapply(hsamples, function(x) x["theta"])
  inference <- list()
  for (i in 1:length(hsamples)) {
    thetas <- hsamples[[i]]$theta
    inference [[i]] <- fun (thetas)
  }
  inf2 <- unlist(inference)
  dim3 <- c(dim(inference[[1]]), length(inf2)/prod(dim(inference[[1]])))
  dim(inf2) <- dim3
  effect <- apply(inf2, c(4), mean)
  names(effect) <- names(samples)
  data.frame(effect)
}

minp <- function (effect) min(ecdf(effect)(0), 1 - ecdf(effect)(0))

zandp <- function (samples, fun) {
  effect <- group.inference.dist(samples, fun)
  Z <- mean(effect) / sd(effect)
  p <- minp(effect)
  round(data.frame(Z, p), 3)
}

mean.sd <- function (samples, fun) {
  effect <- group.inference.dist(samples, fun)
  M <- median(effect)
  SD <- sd(effect)
  round(data.frame(M, SD), 3)
}


# -------------------------------------------------------------------------

# Load samples
print(load("samples/sTPPM_full.RData"))
samples <- samples1
length(samples)
samples[[1]]$p.names

# -------------------------------------------------------------------------
# Threshold effects

# Time pressure contrast (mean 6s minus mean 3s)
fun <- function (thetas) {
  (thetas[,"B.6s10C",, drop=F] + thetas[,"B.6s10N",, drop=F] + 
     thetas[,"B.6s30C",, drop=F] + thetas[,"B.6s30N",, drop=F])/4 -
    (thetas[,"B.3s10C",, drop=F] + thetas[,"B.3s10N",, drop=F] + 
       thetas[,"B.3s30C",, drop=F] + thetas[,"B.3s30N",, drop=F])/4 
}

mean.sd(samples, fun)
zandp(samples, fun)

# Get raw samples to plot effect density
effect <- as.numeric(group.inference.dist(samples, fun))
effect <- data.frame(effect = effect)
head(effect) 
mean(effect$effect)

# Density plot
B_TP_plot <- ggplot(effect, aes(x = effect)) + 
  geom_histogram(aes(y = stat(density)), binwidth = 0.01,
                 alpha = 0.1,
                 col = "black",
                 fill = "steelblue",
                 size = 0.3) +
  geom_density(alpha = 0.1,
               linetype = "dashed",
               size = 0.3) +
  geom_vline(xintercept = quantile(as.numeric(group.inference.dist(samples, fun)), probs = c(0.025, 0.975)),
             alpha = 0.5, linetype = "dashed",
             size = 0.5) +
  geom_vline(aes(xintercept = 0),
             alpha = 0.8,
             col = "red",
             size = 1) +
  labs(title = "Threshold time pressure contrast", 
       subtitle = paste("M =", mean.sd(samples, fun)$M,
                        " Z =", zandp(samples, fun)$Z,
                        " p =", zandp(samples, fun)$p),
       x = "Contrast", 
       y = "Density") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0),
    plot.subtitle = element_text(hjust = 0.5)
  )
B_TP_plot


# PM contrast (mean 30% minus mean 10%)
fun <- function (thetas) {
  (thetas[,"B.3s30C",, drop=F] + thetas[,"B.3s30N",, drop=F] + 
     thetas[,"B.6s30C",, drop=F] + thetas[,"B.6s30N",, drop=F])/4 -
    (thetas[,"B.3s10C",, drop=F] + thetas[,"B.3s10N",, drop=F] + 
       thetas[,"B.6s10C",, drop=F] + thetas[,"B.6s10N",, drop=F])/4 
}

mean.sd(samples, fun)
zandp(samples, fun)

# Get raw samples to plot effect density
effect <- as.numeric(group.inference.dist(samples, fun))
effect <- data.frame(effect = effect)
head(effect) 
mean(effect$effect)

# Density plot
B_PM_plot <- ggplot(effect, aes(x = effect)) + 
  geom_histogram(aes(y = stat(density)), binwidth = 0.01,
                 alpha = 0.1,
                 col = "black",
                 fill = "steelblue",
                 size = 0.3) +
  geom_density(alpha = 0.1,
               linetype = "dashed",
               size = 0.3) +
  geom_vline(xintercept = quantile(as.numeric(group.inference.dist(samples, fun)), probs = c(0.025, 0.975)),
             alpha = 0.5, linetype = "dashed",
             size = 0.5) +
  geom_vline(aes(xintercept = 0),
             alpha = 0.8,
             col = "red",
             size = 1) +
  labs(title = "Threshold PM contrast", 
       subtitle = paste("M =", mean.sd(samples, fun)$M,
                        " Z =", zandp(samples, fun)$Z,
                        " p =", zandp(samples, fun)$p),
       x = "Contrast", 
       y = "Density") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0),
    plot.subtitle = element_text(hjust = 0.5)
  )
B_PM_plot


# -------------------------------------------------------------------------
# Rate quality effects

# TP contrast (mean 6s minus mean 3s)
fun <- function (thetas) {
  ((thetas[,"mean_v.cc6s10C",, drop=F] + 
      thetas[,"mean_v.cc6s30C",, drop=F] +
      thetas[,"mean_v.nn6s10N",, drop=F] + 
      thetas[,"mean_v.nn6s30N",, drop=F])/4 -
     (thetas[,"mean_v.cc6s10N",, drop=F] + 
        thetas[,"mean_v.cc6s30N",, drop=F] +
        thetas[,"mean_v.nn6s10C",, drop=F] + 
        thetas[,"mean_v.nn6s30C",, drop=F])/4) -
    ((thetas[,"mean_v.cc3s10C",, drop=F] + 
        thetas[,"mean_v.cc3s30C",, drop=F] +
        thetas[,"mean_v.nn3s10N",, drop=F] + 
        thetas[,"mean_v.nn3s30N",, drop=F])/4 -
       (thetas[,"mean_v.cc3s10N",, drop=F] + 
          thetas[,"mean_v.cc3s30N",, drop=F] +
          thetas[,"mean_v.nn3s10C",, drop=F] + 
          thetas[,"mean_v.nn3s30C",, drop=F])/4) 
}

mean.sd(samples, fun)
zandp(samples, fun)

# Get raw samples to plot effect density
effect <- as.numeric(group.inference.dist(samples, fun))
effect <- data.frame(effect = effect)
head(effect) 
mean(effect$effect)

# Density plot
quality_TP_plot <- ggplot(effect, aes(x = effect)) + 
  geom_histogram(aes(y = stat(density)), binwidth = 0.005,
                 alpha = 0.1,
                 col = "black",
                 fill = "steelblue",
                 size = 0.3) +
  geom_density(alpha = 0.1,
               linetype = "dashed",
               size = 0.3) +
  geom_vline(xintercept = quantile(as.numeric(group.inference.dist(samples, fun)), probs = c(0.025, 0.975)),
             alpha = 0.5, linetype = "dashed",
             size = 0.5) +
  geom_vline(aes(xintercept = 0),
             alpha = 0.8,
             col = "red",
             size = 1) +
  labs(title = "Rate quality time pressure contrast", 
       subtitle = paste("M =", mean.sd(samples, fun)$M,
                        " Z =", zandp(samples, fun)$Z,
                        " p =", zandp(samples, fun)$p),
       x = "Contrast", 
       y = "Density") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0),
    plot.subtitle = element_text(hjust = 0.5)
  )
quality_TP_plot


# PM contrast (mean 10% minus mean 30%)
fun <- function (thetas) {
  ((thetas[,"mean_v.cc3s10C",, drop=F] + 
      thetas[,"mean_v.cc6s10C",, drop=F] +
      thetas[,"mean_v.nn3s10N",, drop=F] + 
      thetas[,"mean_v.nn6s10N",, drop=F])/4 -
     (thetas[,"mean_v.cc3s10N",, drop=F] + 
        thetas[,"mean_v.cc6s10N",, drop=F] +
        thetas[,"mean_v.nn3s10C",, drop=F] + 
        thetas[,"mean_v.nn6s10C",, drop=F])/4) -
    ((thetas[,"mean_v.cc3s30C",, drop=F] + 
        thetas[,"mean_v.cc6s30C",, drop=F] +
        thetas[,"mean_v.nn3s30N",, drop=F] + 
        thetas[,"mean_v.nn6s30N",, drop=F])/4 -
       (thetas[,"mean_v.cc3s30N",, drop=F] + 
          thetas[,"mean_v.cc6s30N",, drop=F] +
          thetas[,"mean_v.nn3s30C",, drop=F] + 
          thetas[,"mean_v.nn6s30C",, drop=F])/4) 
}

mean.sd(samples, fun)
zandp(samples, fun)

# Get raw samples to plot effect density
effect <- as.numeric(group.inference.dist(samples, fun))
effect <- data.frame(effect = effect)
head(effect) 
mean(effect$effect)

# Density plot
quality_PM_plot <- ggplot(effect, aes(x = effect)) + 
  geom_histogram(aes(y = stat(density)), binwidth = 0.005,
                 alpha = 0.1,
                 col = "black",
                 fill = "steelblue",
                 size = 0.3) +
  geom_density(alpha = 0.1,
               linetype = "dashed",
               size = 0.3) +
  geom_vline(xintercept = quantile(as.numeric(group.inference.dist(samples, fun)), probs = c(0.025, 0.975)),
             alpha = 0.5, linetype = "dashed",
             size = 0.5) +
  geom_vline(aes(xintercept = 0),
             alpha = 0.8,
             col = "red",
             size = 1) +
  labs(title = "Rate quality PM contrast", 
       subtitle = paste("M =", mean.sd(samples, fun)$M,
                        " Z =", zandp(samples, fun)$Z,
                        " p =", zandp(samples, fun)$p),
       x = "Contrast", 
       y = "Density") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0),
    plot.subtitle = element_text(hjust = 0.5)
  )
quality_PM_plot


# TP x PM interaction (10%-30% at low TP minus 10%-30% at high TP)
fun <- function (thetas) {
  (((thetas[,"mean_v.cc6s10C",, drop=F] + 
       thetas[,"mean_v.nn6s10N",, drop=F])/2 -
      (thetas[,"mean_v.cc6s10N",, drop=F] + 
         thetas[,"mean_v.nn6s10C",, drop=F])/2) -
     ((thetas[,"mean_v.cc6s30C",, drop=F] + 
         thetas[,"mean_v.nn6s30N",, drop=F])/2 -
        (thetas[,"mean_v.cc6s30N",, drop=F] + 
           thetas[,"mean_v.nn6s30C",, drop=F])/2))-
    (((thetas[,"mean_v.cc3s10C",, drop=F] + 
         thetas[,"mean_v.nn3s10N",, drop=F])/2 -
        (thetas[,"mean_v.cc3s10N",, drop=F] + 
           thetas[,"mean_v.nn3s10C",, drop=F])/2) -
       ((thetas[,"mean_v.cc3s30C",, drop=F] + 
           thetas[,"mean_v.nn3s30N",, drop=F])/2 -
          (thetas[,"mean_v.cc3s30N",, drop=F] + 
             thetas[,"mean_v.nn3s30C",, drop=F])/2))
}

mean.sd(samples, fun)
zandp(samples, fun)

# Get raw samples to plot effect density
effect <- as.numeric(group.inference.dist(samples, fun))
effect <- data.frame(effect = effect)
head(effect) 
mean(effect$effect)

# Density plot
quality_TPxPM_plot <- ggplot(effect, aes(x = effect)) + 
  geom_histogram(aes(y = stat(density)), binwidth = 0.005,
                 alpha = 0.1,
                 col = "black",
                 fill = "steelblue",
                 size = 0.3) +
  geom_density(alpha = 0.1,
               linetype = "dashed",
               size = 0.3) +
  geom_vline(xintercept = quantile(as.numeric(group.inference.dist(samples, fun)), probs = c(0.025, 0.975)),
             alpha = 0.5, linetype = "dashed",
             size = 0.5) +
  geom_vline(aes(xintercept = 0),
             alpha = 0.8,
             col = "red",
             size = 1) +
  labs(title = "Rate quality TP x PM interaction", 
       subtitle = paste("M =", mean.sd(samples, fun)$M,
                        " Z =", zandp(samples, fun)$Z,
                        " p =", zandp(samples, fun)$p),
       x = "Contrast", 
       y = "Density") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0),
    plot.subtitle = element_text(hjust = 0.5)
  )
quality_TPxPM_plot


# -------------------------------------------------------------------------
# Rate quantity effects

# TP contrast (mean 6s minus mean 3s)
fun <- function (thetas) {
  ((thetas[,"mean_v.cc6s10C",, drop=F] + 
      thetas[,"mean_v.cc6s30C",, drop=F] +
      thetas[,"mean_v.nn6s10N",, drop=F] + 
      thetas[,"mean_v.nn6s30N",, drop=F] +
      thetas[,"mean_v.cc6s10N",, drop=F] + 
      thetas[,"mean_v.cc6s30N",, drop=F] +
      thetas[,"mean_v.nn6s10C",, drop=F] + 
      thetas[,"mean_v.nn6s30C",, drop=F])/8) -
    ((thetas[,"mean_v.cc3s10C",, drop=F] + 
        thetas[,"mean_v.cc3s30C",, drop=F] +
        thetas[,"mean_v.nn3s10N",, drop=F] + 
        thetas[,"mean_v.nn3s30N",, drop=F] +
        thetas[,"mean_v.cc3s10N",, drop=F] + 
        thetas[,"mean_v.cc3s30N",, drop=F] +
        thetas[,"mean_v.nn3s10C",, drop=F] + 
        thetas[,"mean_v.nn3s30C",, drop=F])/8) 
}

mean.sd(samples, fun)
zandp(samples, fun)

# Get raw samples to plot effect density
effect <- as.numeric(group.inference.dist(samples, fun))
effect <- data.frame(effect = effect)
head(effect) 
mean(effect$effect)

# Density plot
quantity_TP_plot <- ggplot(effect, aes(x = effect)) + 
  geom_histogram(aes(y = stat(density)), binwidth = 0.005,
                 alpha = 0.1,
                 col = "black",
                 fill = "steelblue",
                 size = 0.3) +
  geom_density(alpha = 0.1,
               linetype = "dashed",
               size = 0.3) +
  geom_vline(xintercept = quantile(as.numeric(group.inference.dist(samples, fun)), probs = c(0.025, 0.975)),
             alpha = 0.5, linetype = "dashed",
             size = 0.5) +
  geom_vline(aes(xintercept = 0),
             alpha = 0.8,
             col = "red",
             size = 1) +
  labs(title = "Rate quantity time pressure contrast", 
       subtitle = paste("M =", mean.sd(samples, fun)$M,
                        " Z =", zandp(samples, fun)$Z,
                        " p =", zandp(samples, fun)$p),
       x = "Contrast", 
       y = "Density") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0),
    plot.subtitle = element_text(hjust = 0.5)
  )
quantity_TP_plot


# PM contrast (mean 10% minus mean 30%)
fun <- function (thetas) {
  ((thetas[,"mean_v.cc3s10C",, drop=F] + 
      thetas[,"mean_v.cc6s10C",, drop=F] +
      thetas[,"mean_v.nn3s10N",, drop=F] + 
      thetas[,"mean_v.nn6s10N",, drop=F] +
      thetas[,"mean_v.cc3s10N",, drop=F] + 
      thetas[,"mean_v.cc6s10N",, drop=F] +
      thetas[,"mean_v.nn3s10C",, drop=F] + 
      thetas[,"mean_v.nn6s10C",, drop=F])/8) -
    ((thetas[,"mean_v.cc3s30C",, drop=F] + 
        thetas[,"mean_v.cc6s30C",, drop=F] +
        thetas[,"mean_v.nn3s30N",, drop=F] + 
        thetas[,"mean_v.nn6s30N",, drop=F] +
        thetas[,"mean_v.cc3s30N",, drop=F] + 
        thetas[,"mean_v.cc6s30N",, drop=F] +
        thetas[,"mean_v.nn3s30C",, drop=F] + 
        thetas[,"mean_v.nn6s30C",, drop=F])/8) 
}

mean.sd(samples, fun)
zandp(samples, fun)

# Get raw samples to plot effect density
effect <- as.numeric(group.inference.dist(samples, fun))
effect <- data.frame(effect = effect)
head(effect) 
mean(effect$effect)

# Density plot
quantity_PM_plot <- ggplot(effect, aes(x = effect)) + 
  geom_histogram(aes(y = stat(density)), binwidth = 0.005,
                 alpha = 0.1,
                 col = "black",
                 fill = "steelblue",
                 size = 0.3) +
  geom_density(alpha = 0.1,
               linetype = "dashed",
               size = 0.3) +
  geom_vline(xintercept = quantile(as.numeric(group.inference.dist(samples, fun)), probs = c(0.025, 0.975)),
             alpha = 0.5, linetype = "dashed",
             size = 0.5) +
  geom_vline(aes(xintercept = 0),
             alpha = 0.8,
             col = "red",
             size = 1) +
  labs(title = "Rate quantity PM contrast", 
       subtitle = paste("M =", mean.sd(samples, fun)$M,
                        " Z =", zandp(samples, fun)$Z,
                        " p =", zandp(samples, fun)$p),
       x = "Contrast", 
       y = "Density") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0),
    plot.subtitle = element_text(hjust = 0.5)
  )
quantity_PM_plot


# -------------------------------------------------------------------------

# Reactive control

# Overall reactive control
# Ongoing task rates lower with PM target present versus absent
fun <- function (thetas) {
  ((thetas[,"mean_v.cc3s10C",, drop=F] + 
      thetas[,"mean_v.cc6s10C",, drop=F] +
      thetas[,"mean_v.nn3s10N",, drop=F] + 
      thetas[,"mean_v.nn6s10N",, drop=F] +    
     thetas[,"mean_v.cc3s30C",, drop=F] + 
        thetas[,"mean_v.cc6s30C",, drop=F] +
        thetas[,"mean_v.nn3s30N",, drop=F] + 
        thetas[,"mean_v.nn6s30N",, drop=F])/8) -
     ((thetas[,"mean_v.pc3s10C",, drop=F] + 
         thetas[,"mean_v.pc6s10C",, drop=F] +
         thetas[,"mean_v.pn3s10N",, drop=F] + 
         thetas[,"mean_v.pn6s10N",, drop=F] +    
         thetas[,"mean_v.pc3s30C",, drop=F] + 
         thetas[,"mean_v.pc6s30C",, drop=F] +
         thetas[,"mean_v.pn3s30N",, drop=F] + 
         thetas[,"mean_v.pn6s30N",, drop=F])/8)
}

mean.sd(samples, fun)
zandp(samples, fun)

# Get raw samples to plot effect density
effect <- as.numeric(group.inference.dist(samples, fun))
effect <- data.frame(effect = effect)
head(effect) 
mean(effect$effect)

# Density plot
reactive_total_plot <- ggplot(effect, aes(x = effect)) + 
  geom_histogram(aes(y = stat(density)), binwidth = 0.005,
                 alpha = 0.1,
                 col = "black",
                 fill = "steelblue",
                 size = 0.3) +
  geom_density(alpha = 0.1,
               linetype = "dashed",
               size = 0.3) +
  geom_vline(xintercept = quantile(as.numeric(group.inference.dist(samples, fun)), probs = c(0.025, 0.975)),
             alpha = 0.5, linetype = "dashed",
             size = 0.5) +
  geom_vline(aes(xintercept = 0),
             alpha = 0.8,
             col = "red",
             size = 1) +
  labs(title = "Rate quantity PM contrast", 
       subtitle = paste("M =", mean.sd(samples, fun)$M,
                        " Z =", zandp(samples, fun)$Z,
                        " p =", zandp(samples, fun)$p),
       x = "Contrast", 
       y = "Density") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0),
    plot.subtitle = element_text(hjust = 0.5)
  )
reactive_total_plot


# PM contrast (mean 10% minus mean 30%)
# Reactive control not stronger for more/less frequent PM
fun <- function (thetas) {
  ((thetas[,"mean_v.cc3s10C",, drop=F] + 
      thetas[,"mean_v.cc6s10C",, drop=F] +
      thetas[,"mean_v.nn3s10N",, drop=F] + 
      thetas[,"mean_v.nn6s10N",, drop=F])/4 -
     (thetas[,"mean_v.pc3s10C",, drop=F] + 
        thetas[,"mean_v.pc6s10C",, drop=F] +
        thetas[,"mean_v.pn3s10N",, drop=F] + 
        thetas[,"mean_v.pn6s10N",, drop=F])/4) -
    ((thetas[,"mean_v.cc3s30C",, drop=F] + 
        thetas[,"mean_v.cc6s30C",, drop=F] +
        thetas[,"mean_v.nn3s30N",, drop=F] + 
        thetas[,"mean_v.nn6s30N",, drop=F])/4 -
       (thetas[,"mean_v.pc3s30C",, drop=F] + 
          thetas[,"mean_v.pc6s30C",, drop=F] +
          thetas[,"mean_v.pn3s30N",, drop=F] + 
          thetas[,"mean_v.pn6s30N",, drop=F])/4) 
}

mean.sd(samples, fun)
zandp(samples, fun)

# Get raw samples to plot effect density
effect <- as.numeric(group.inference.dist(samples, fun))
effect <- data.frame(effect = effect)
head(effect) 
mean(effect$effect)

# Density plot
reactive_PM_plot <- ggplot(effect, aes(x = effect)) + 
  geom_histogram(aes(y = stat(density)), binwidth = 0.005,
                 alpha = 0.1,
                 col = "black",
                 fill = "steelblue",
                 size = 0.3) +
  geom_density(alpha = 0.1,
               linetype = "dashed",
               size = 0.3) +
  geom_vline(xintercept = quantile(as.numeric(group.inference.dist(samples, fun)), probs = c(0.025, 0.975)),
             alpha = 0.5, linetype = "dashed",
             size = 0.5) +
  geom_vline(aes(xintercept = 0),
             alpha = 0.8,
             col = "red",
             size = 1) +
  labs(title = "Reactive control PM contrast", 
       subtitle = paste("M =", mean.sd(samples, fun)$M,
                        " Z =", zandp(samples, fun)$Z,
                        " p =", zandp(samples, fun)$p),
       x = "Contrast", 
       y = "Density") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0),
    plot.subtitle = element_text(hjust = 0.5)
  )
reactive_PM_plot


# TP contrast (mean 3s minus mean 6s)
# Reactive control weaker at high time pressure
fun <- function (thetas) {
  ((thetas[,"mean_v.cc3s10C",, drop=F] + 
      thetas[,"mean_v.cc3s30C",, drop=F] +
      thetas[,"mean_v.nn3s10N",, drop=F] + 
      thetas[,"mean_v.nn3s30N",, drop=F])/4 -
     (thetas[,"mean_v.pc3s10C",, drop=F] + 
        thetas[,"mean_v.pc3s30C",, drop=F] +
        thetas[,"mean_v.pn3s10N",, drop=F] + 
        thetas[,"mean_v.pn3s30N",, drop=F])/4) -
    ((thetas[,"mean_v.cc6s10C",, drop=F] + 
        thetas[,"mean_v.cc6s30C",, drop=F] +
        thetas[,"mean_v.nn6s10N",, drop=F] + 
        thetas[,"mean_v.nn6s30N",, drop=F])/4 -
       (thetas[,"mean_v.pc6s10C",, drop=F] + 
          thetas[,"mean_v.pc6s30C",, drop=F] +
          thetas[,"mean_v.pn6s10N",, drop=F] + 
          thetas[,"mean_v.pn6s30N",, drop=F])/4) 
}

mean.sd(samples, fun)
zandp(samples, fun)

# Get raw samples to plot effect density
effect <- as.numeric(group.inference.dist(samples, fun))
effect <- data.frame(effect = effect)
head(effect) 
mean(effect$effect)

# Density plot
reactive_TP_plot <- ggplot(effect, aes(x = effect)) + 
  geom_histogram(aes(y = stat(density)), binwidth = 0.005,
                 alpha = 0.1,
                 col = "black",
                 fill = "steelblue",
                 size = 0.3) +
  geom_density(alpha = 0.1,
               linetype = "dashed",
               size = 0.3) +
  geom_vline(xintercept = quantile(as.numeric(group.inference.dist(samples, fun)), probs = c(0.025, 0.975)),
             alpha = 0.5, linetype = "dashed",
             size = 0.5) +
  geom_vline(aes(xintercept = 0),
             alpha = 0.8,
             col = "red",
             size = 1) +
  labs(title = "Reactive control PM contrast", 
       subtitle = paste("M =", mean.sd(samples, fun)$M,
                        " Z =", zandp(samples, fun)$Z,
                        " p =", zandp(samples, fun)$p),
       x = "Contrast", 
       y = "Density") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0),
    plot.subtitle = element_text(hjust = 0.5)
  )
reactive_TP_plot

