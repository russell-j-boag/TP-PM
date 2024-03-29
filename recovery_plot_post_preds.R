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
print(load("samples/recov_sTPPM_full_sdvS.RData"))
samples <- recov_samples1

# Bring longer thetas down to min nmc by sampling
nmcs <- sapply(samples, function(x) x$nmc)
nmc  <- min(nmcs)
for (i in 1:length(samples)) if (nmcs[i] > nmc) samples[[i]]$theta <-
  samples[[i]]$theta[,,sample(1:dim(samples[[i]]$theta)[3], nmc)]

ci <- subject.average.ci(samples)
round(ci,3)


# Check how many runs it took to converge
# If any say "NA" then it didn't converge
for(i in 1:length(samples)) {
  print(attr(samples[[i]], "auto"))
}

# Check samples are all the same length for every participant
for(i in 1:length(samples)) {
  print(samples[[i]]$nmc)
}

# How many parameters? How many chains?
# str(samples[[1]])
# names(samples[[1]])
# samples[[1]]$theta

# Checkout some trace plots for participants and confirm that the samples look
# like mcmc samples should look.
# plot.dmc(samples[[1]])

# Plot the posterior parameter distributions against the priors to make sure
# that you don't have huge prior influence
# plot.dmc(samples[[1]], p.prior = samples[[1]]$p.prior)


# Plot mean accuracy and RT -----------------------------------------------

# Generate posterior predictives
# post_pred_sims <- h.post.predict.dmc(samples, save.simulation = TRUE)
# save(post_pred_sims, file = "deriv/post_pred_sims_recov_full_sdvS.RData")

# Load posterior predictives
print(load("deriv/post_pred_sims_recov_full_sdvS.RData"))

# Stack into data frame
sims <- do.call(rbind, post_pred_sims)
# str(sims)
sims$correct <- NULL
dim(sims)
head(sims, 10)
sims <- sims[,c("reps","TP","PM","S","R","RT")]

# Do the same for the data
data <- lapply(post_pred_sims, function(x) attr(x, "data"))
data <- do.call(rbind, data)
data$correct <- NULL
head(data, 10)
data <- data[,c("TP","PM","S","R","RT")]

# Create match function to code accuracy
matchfun = function (d) {
  (d$S=="cc" & d$R=="C")|(d$S=="nn" & d$R=="N")|
    (d$S=="pc" & d$R=="P")|(d$S=="pn" & d$R=="P")
}

table(matchfun(data))
table(matchfun(sims))

# Get a gglist for plotting (as stored in post_pred_sims)
GGLIST <- get.fitgglist.dmc(sims, data, acc.fun = matchfun)
GGLIST

# Remove NAs
GGLIST <- lapply(GGLIST, function(x)  x[is.finite(x$data),])
names(GGLIST)


# Correct response proportion ---------------------------------------------

# Re-order columns for desired panel order
rp <- GGLIST$pps
dim(rp)
head(rp, 10)
# rp <- rp[,c(1,3,2,4,5,6,7,8)]

matchfun(rp)

# Take only the correct responses and drop the R column
rp_correct <- rp[ matchfun(rp),
                  !(names(rp) %in% c("R")) ]
rp_correct

levels(rp_correct$TP) 
levels(rp_correct$PM) 
levels(rp_correct$PM) <- c("10% PM", "30% PM")
levels(rp_correct$S)
levels(rp_correct$S) <- c("Conflict", "Non-conflict", "PM conflict", "PM non-conflict")
str(rp_correct)
head(rp_correct, 10)


# Ongoing task accuracy
rp_correct_ongoing_plot <- rp_correct[rp_correct$S == "Conflict"|rp_correct$S == "Non-conflict",] %>%
  ggplot(aes(x = TP)) +
  geom_point(stat = "identity", aes(y = data), size = 3, shape = 1) +
  geom_point(stat = "identity", aes(y = median), size = 2) +
  geom_line(aes(y = median, group = S), 
            linetype = "dashed", 
            size = 0.8) +
  geom_errorbar(aes(ymin = lower, 
                    ymax = upper, 
                    width = 0.3)) +
  ylim(0.6, 1) +
  facet_grid(S ~ PM) +
  labs(title = "Ongoing task accuracy", 
       x = "Time pressure", 
       y = "Response proportion", 
       color = "Stimulus",
       shape = "Stimulus") +
  theme_minimal()
rp_correct_ongoing_plot 


# PM accuracy
rp_correct_PM_plot <- rp_correct[rp_correct$S == "PM conflict"|rp_correct$S == "PM non-conflict",] %>%
  ggplot(aes(x = TP)) +
  geom_point(stat = "identity", aes(y = data), size = 3, shape = 1) +
  geom_point(stat = "identity", aes(y = median), size = 2) +
  geom_line(aes(y = median, group = S), 
            linetype = "dashed", 
            size = 0.8) +
  geom_errorbar(aes(ymin = lower, 
                    ymax = upper, 
                    width = 0.3)) +
  ylim(0.6, 1) +
  facet_grid(S ~ PM) +
  labs(title = "PM accuracy", 
       x = "Time pressure", 
       y = "Response proportion", 
       color = "Stimulus",
       shape = "Stimulus") +
  theme_minimal()
rp_correct_PM_plot 



# Correct response time ---------------------------------------------------

# Re-order columns for desired panel order
rt <- GGLIST$RTs
dim(rt)
head(rt, 10)
# rt <- rt[,c(1,3,2,4,5,6,7,8)]


# Take only the correct responses and drop the R column
rt_correct <- rt[ matchfun(rt),
                  !(names(rt) %in% c("R")) ]
rt_correct

levels(rt_correct$TP) 
levels(rt_correct$PM) 
levels(rt_correct$PM) <- c("10% PM", "30% PM")
levels(rt_correct$S)
levels(rt_correct$S) <- c("Conflict", "Non-conflict", "PM conflict", "PM non-conflict")
str(rt_correct)
head(rt_correct, 10)


# Correct ongoing task RT
rt_correct_ongoing_plot <- rt_correct[rt_correct$S == "Conflict"|rt_correct$S == "Non-conflict",] %>%
  ggplot(aes(x = TP)) +
  geom_point(stat = "identity", aes(y = data), size = 3, shape = 1) +
  geom_point(stat = "identity", aes(y = median), size = 2) +
  geom_line(aes(y = median, group = quantile), 
            linetype = "dashed", 
            size = 0.8) +
  geom_errorbar(aes(ymin = lower, 
                    ymax = upper, 
                    width = 0.3)) +
  ylim(0.5, 4.3) +
  facet_grid(S ~ PM) +
  labs(title = "Ongoing task RT (correct responses)", 
       x = "Time pressure", 
       y = "RT (s)", 
       color = "Stimulus",
       shape = "Stimulus") +
  theme_minimal()
rt_correct_ongoing_plot 


# Correct PM RT
rt_correct_PM_plot <- rt_correct[rt_correct$S == "PM conflict"|rt_correct$S == "PM non-conflict",] %>%
  ggplot(aes(x = TP)) +
  geom_point(stat = "identity", aes(y = data), size = 3, shape = 1) +
  geom_point(stat = "identity", aes(y = median), size = 2) +
  geom_line(aes(y = median, group = quantile), 
            linetype = "dashed", 
            size = 0.8) +
  geom_errorbar(aes(ymin = lower, 
                    ymax = upper, 
                    width = 0.3)) +
  ylim(0.5, 3) +
  facet_grid(S ~ PM) +
  labs(title = "PM RT (correct responses)", 
       x = "Time pressure", 
       y = "RT (s)", 
       color = "Stimulus",
       shape = "Stimulus") +
  theme_minimal()
rt_correct_PM_plot 


# Error response time -----------------------------------------------------

# Pull out RTs
rt <- GGLIST$RTs
dim(rt)
head(rt, 10)

# Take only the incorrect responses and drop the R column
rt_error <- rt[ !matchfun(rt), ]
rt_error

levels(rt_error$TP) 
levels(rt_error$PM) 
levels(rt_error$PM) <- c("10% PM", "30% PM")
levels(rt_error$S)
levels(rt_error$S) <- c("Conflict", "Non-conflict", "PM conflict", "PM non-conflict")
levels(rt_error$R) 
levels(rt_error$R) <- c("Conflict", "Non-conflict", "PM")
str(rt_error)
head(rt_error, 10)


# Incorrect ongoing task RT
rt_error_ongoing_plot <- rt_error[(rt_error$S == "Conflict" & rt_error$R == "Non-conflict")|
                                    (rt_error$S == "Non-conflict" & rt_error$R == "Conflict"),] %>%
  ggplot(aes(x = TP)) +
  geom_point(stat = "identity", aes(y = data), size = 3, shape = 1) +
  geom_point(stat = "identity", aes(y = median), size = 2) +
  geom_line(aes(y = median, group = quantile), 
            linetype = "dashed", 
            size = 0.8) +
  geom_errorbar(aes(ymin = lower, 
                    ymax = upper, 
                    width = 0.3)) +
  ylim(0.5, 4.8) +
  facet_grid(S ~ PM) +
  labs(title = "Ongoing task RT (incorrect responses)", 
       x = "Time pressure", 
       y = "RT (s)", 
       color = "Stimulus",
       shape = "Stimulus") +
  theme_minimal()
rt_error_ongoing_plot 


# Incorrect PM RT
rt_error_PM_plot <- rt_error[(rt_error$S != "PM conflict" & rt_error$R == "PM")|
                               (rt_error$S != "PM non-conflict" & rt_error$R == "PM"),] %>%
  ggplot(aes(x = TP)) +
  geom_point(stat = "identity", aes(y = data), size = 3, shape = 1) +
  geom_point(stat = "identity", aes(y = median), size = 2) +
  geom_line(aes(y = median, group = quantile), 
            linetype = "dashed", 
            size = 0.8) +
  geom_errorbar(aes(ymin = lower, 
                    ymax = upper, 
                    width = 0.3)) +
  ylim(0, 10) +
  facet_grid(S ~ PM) +
  labs(title = "PM RT (incorrect responses)", 
       x = "Time pressure", 
       y = "RT (s)", 
       color = "Stimulus",
       shape = "Stimulus") +
  theme_minimal()
rt_error_PM_plot 


# Combine plots -----------------------------------------------------------

library("gridExtra")
library("ggpubr")

ggsave("plots/fits_accuracy_ongoing.png", plot = rp_correct_ongoing_plot, 
       width = 2000, height = 1400, units = "px")

ggsave("plots/fits_accuracy_PM.png", plot = rp_correct_PM_plot, 
       width = 2000, height = 1400, units = "px")

ggsave("plots/fits_correct_RT_PM.png", plot = rt_correct_PM_plot, 
       width = 2000, height = 1400, units = "px")


# Ongoing task RT
# Remove plot legends
rt_correct_ongoing_plot <- rt_correct_ongoing_plot + theme(legend.position = "none",
                                                           axis.title.x = element_blank(),
                                                           axis.text.x = element_blank())

rt_error_ongoing_plot <- rt_error_ongoing_plot + theme(legend.position = "none",
                                                       strip.text.x = element_blank())

ggarrange(rt_correct_ongoing_plot, rt_error_ongoing_plot, nrow = 2,
          common.legend = TRUE, legend = "right")

ggsave("plots/fits_RT_ongoing.png", plot = last_plot(), 
       width = 2200, height = 2000, units = "px")


# PM RT
# Remove plot legends
rt_correct_PM_plot <- rt_correct_PM_plot + theme(legend.position = "none",
                                                 axis.title.x = element_blank(),
                                                 axis.text.x = element_blank())

rt_error_PM_plot <- rt_error_PM_plot + theme(legend.position = "none",
                                             strip.text.x = element_blank())

ggarrange(rt_correct_PM_plot, rt_error_PM_plot, nrow = 2,
          common.legend = TRUE, legend = "right")

ggsave("plots/fits_RT_PM.png", plot = last_plot(), 
       width = 2200, height = 2000, units = "px")
