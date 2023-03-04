
### Below function averages parameters across conditions by grepping out
# from av.posts. So if you set av.posts to match anything containing mean_v,
# it would average all rates together and replace all values with the avg before
# simulating. We use it to parse the effects of rates/thresholds on the manifests
# in terms of experimental factors.

# samples=samples[[3]]
# probs=c(1:99)/100;random=TRUE
# bw="nrd0";report=10;save.simulation=TRUE;factors=NA
# av.posts<-av.posts.threscond
avps.post.predict.dmc = function(samples,n.post=100,probs=c(1:99)/100,random=TRUE,
                                 bw="nrd0",report=10,save.simulation=TRUE,factors=NA, av.posts=c())
  # make list of posterior preditive density, quantiles and response p(robability)
{
  
  get.dqp <- function(sim,facs,probs,n.post=NA) {
    
    quantile.names <- function(x,probs=seq(0, 1, 0.25),na.rm=FALSE,type=7, ...) {
      out <- quantile(x,probs=probs,na.rm=na.rm,type=type,names=FALSE,...)
      names(out) <- probs*100
      if (all(is.na(out))) NULL else out
    }
    
    qs <- tapply(sim$RT,sim[,c(facs,"R")],quantile.names,probs=probs,na.rm=TRUE)
    #     qs <- apply(qs,1:length(dim(qs)),function(x){
    #       if ( is.null(x[[1]]) || all(is.na(x[[1]]))) NULL else x[[1]]})
    
    # get probabilities given not na
    simOK <- sim[!is.na(sim$RT),]
    p <- tapply(simOK$RT,simOK[,c(facs,"R")],length)
    p[is.na(p)] <- 0 # In case some cells are empty
    np <- rep(apply(p,1:length(facs),sum),times=length(levels(simOK$R)))
    p <- p/np
    # get p NA
    pNA <- tapply(sim$RT,sim[,c(facs,"R")],length)
    pNA[is.na(pNA)] <- 0 # In case some cells are empty
    npNA <- rep(apply(pNA,1:length(facs),sum),times=length(levels(sim$R)))
    pNA <- tapply(is.na(sim$RT),sim[,c(facs,"R")],sum)
    pNA[is.na(pNA)] <- 0 # In case some cells are empty
    pNA <- pNA/npNA
    
    # For a simulation get proability replicates
    if (!is.na(n.post)) {
      repfac <- rep(1:n.post,each=sum(ns))
      repfac <- repfac[!is.na(sim$RT)]
      ps <- tapply(simOK$RT,cbind(simOK[,c(facs,"R")],rep=repfac),length)
      ps[is.na(ps)] <- 0 # In case some cells are empty
      ps <- n.post*ps/np
      # and NA replicates
      pNAs <- tapply(is.na(sim$RT),cbind(sim[,c(facs,"R")],rep=rep(1:n.post,each=sum(ns))),sum)
      pNAs[is.na(pNAs)] <- 0 # In case some cells are empty
      pNAs <- n.post*pNAs/npNA
    } else {
      ps=NULL
      pNAs=NULL
    }
    
    # cell names
    cell.names <- dimnames(qs)[[1]]
    n.cell <- length(facs)+1
    for ( i in 2:n.cell )
      cell.names <- outer(cell.names,dimnames(qs)[[i]],"paste",sep=".")
    cell.names <- as.vector(cell.names)
    # Get density and make defective
    dens <- tapply(sim$RT,sim[,c(facs,"R")],function(x){
      if (all(is.na(x))) NULL else {
        x <- x[x>=quantile(x,.01,na.rm=TRUE) & x<=quantile(x,.99,na.rm=TRUE)]
        if (length(x[!is.na(x)])<2) NULL else density(x[!is.na(x)],bw=bw)
      }
    })
    for (i in 1:length(p)) if ( is.finite(p[i]) && !(p[i]==0) )
    {
      if (!is.null(qs[i][[1]])) {
        names(qs[i][[1]]) <- as.numeric(names(qs[i][[1]]))*p[i]/100
        attr(qs[i][[1]],"cell.name") <- cell.names[i]
      }
      if (!is.null(dens[i][[1]]) ) {
        dens[i][[1]]$y <- dens[i][[1]]$y*p[i]
        attr(dens[i][[1]],"cell.name") <- cell.names[i]
      }
    }
    dnd <- dimnames(dens)
    dnq <- dimnames(qs)
    dens <- apply(dens,1:length(facs),function(x){x})
    qs <- apply(qs,1:length(facs),function(x){x})
    if ( is.null(dim(dens)) ) {
      dens <- array(dens,dim=c(length(dens)))
      dimnames(dens) <- dnd[-length(dnd)]
      qs <- array(qs,dim=c(length(qs)))
      dimnames(qs) <- dnq[-length(dnq)]
    }
    list(pdf=dens,cdf=qs,p=p,ps=ps,pNA=pNA,pNAs=pNAs)
  }
  
  model <- attributes(samples$data)$model
  facs <- names(attr(model,"factors"))
  if (any(is.na(factors))) factors <- facs
  if (!all(factors %in% facs))
    stop(paste("Factors argument must contain one or more of:",paste(facs,collapse=",")))
  resp <- names(attr(model,"responses"))
  ns <- table(samples$data[,facs],dnn=facs)
  n.par <- dim(samples$theta)[2]
  thetas <- matrix(aperm(samples$theta,c(3,1,2)),ncol=n.par)
  colnames(thetas) <- dimnames(samples$theta)[[2]]
  if (is.na(n.post)) use <- c(1:dim(thetas)[1]) else {
    if (random) use <- sample(c(1:dim(thetas)[1]),n.post,replace=F) else
      use <- round(seq(1,dim(thetas)[1],length.out=n.post))
  }
  n.post <- length(use)
  posts <- thetas[use,]
  
  
  cat("Below is how I'm averaging (each row is averaged). If this is wrong, adjust your
      av.posts to grep correctly.")
  for (q in 1:length(av.posts)) print(colnames(posts[, grep(av.posts[q], colnames(posts))]))
  ### tweak to average av.posts
  q=1
  
  if(length(av.posts)!= 0) {
    ### loop through all the averaged posts
    for (q in 1:length(av.posts)) {
      
      num.params <- dim(posts[, grep(av.posts[q], colnames(posts))])[2]
      average.params <- rowMeans(posts[, grep(av.posts[q], colnames(posts))])
      posts[, grep(av.posts[q], colnames(posts))] <- matrix(average.params,nrow=length(average.params),ncol=num.params,byrow=F)
      
    }
  }
  
  ###
  n.rep <- sum(ns)
  sim <- data.frame(matrix(nrow=n.post*n.rep,ncol=dim(samples$data)[2]))
  names(sim) <- names(samples$data)
  # Tweaks for Stop Signal
  if ( !any(names(samples$data)=="SSD") ) {
    SSD <- rep(Inf,sum(ns))
    leave.out <- -c(1:dim(samples$data)[2])[names(samples$data)=="RT"]
  } else {
    # Assumes last two are SSD and RT! FIX ME.
    if (is.null(facs)) SSD <- samples$data$SSD else
      SSD <- unlist(tapply(samples$data$SSD,samples$data[,facs],identity))
    leave.out <- -c((dim(samples$data)[2]-1):dim(samples$data)[2])
  }
  cat("\n")
  cat(paste("Simulating (\'.\'=",report,"): ",sep=""))
  for (i in names(samples$data)[leave.out])
    sim[[i]] <- factor(rep(NA,n.post*n.rep),levels=levels(samples$data[[i]]))
  for (i in 1:n.post) {
    tmp <- simulate.dmc(posts[i,],model,ns,SSD=SSD)
    sim[(1+(i-1)*n.rep):(i*n.rep),names(tmp)] <- tmp
    if ( (i %% report) == 0) cat(".")
  }
  if ( any(names(sim)=="R2") ) { # MTR model
    sim$R <- paste(as.character(sim$R),as.character(sim$R2),sep="")
    sim$R[sim$R2=="DK"] <- "DK"
    sim$R <- factor(sim$R)
    samples$data$R <- paste(as.character(samples$data$R),as.character(samples$data$R2),sep="")
    samples$data$R[samples$data$R2=="DK"] <- "DK"
    samples$data$R <- factor(samples$data$R)
  }
  if (save.simulation) {
    sim <- cbind(reps=rep(1:n.post,each=dim(samples$data)[1]),sim)
    attr(sim,"data") <- samples$data
    sim
  } else {
    sim.dqp <- get.dqp(sim,factors,probs,n.post)
    dat.dqp <- get.dqp(sim=samples$data,factors,probs)
    names(dat.dqp) <- paste("data",names(dat.dqp),sep=".")
    c(sim.dqp,dat.dqp)
  }
}


# lapplys the above function on everybody
avps.h.post.predict.dmc <- function(samples,n.post=100,probs=c(1:99)/100,
                                    bw="nrd0",save.simulation=FALSE,av.posts=c()){
  # lapply post.predict to each subject
  lapply(samples,avps.post.predict.dmc,n.post=n.post,probs=probs,bw=bw,
         save.simulation=save.simulation, av.posts=av.posts)
  
}


get.effects.dmc <- function (PPs, fun = function (x) {mean (x)}, lower=.025, upper=.975) {
  
  simdata <- do.call(rbind, PPs)
  data <- lapply(PPs, function(x) attr(x, "data"))
  data <- do.call(rbind, data)
  nreps=max(PPs[[1]]$reps)
  
  data.effects <- fun(data)
  noutput <- length(data.effects)
  
  sim.effects <- matrix(NA, nrow = nreps, ncol = noutput + 1)
  sim.effects[ , noutput + 1] <- 1:nreps
  
  colnames(sim.effects) <- c(names(data.effects), "n.rep")
  
  # Calculate effects separately for each rep
  for (j in 1:nreps) {
    
    currentsim.effects <- simdata[simdata$reps == j, ]
    sim.effects[j, 1:noutput] <- as.numeric(fun(currentsim.effects))
    
  }
  
  # Get a ggplot df with posterior mean, lower, and upper.
  effects.ggdf <-  t(apply(sim.effects, c(2), 
                           function(x) c(mean(x), quantile(x, probs = c(lower, upper)))))
  effects.ggdf <- data.frame(effects.ggdf)
  effects.ggdf <- effects.ggdf[(!rownames(effects.ggdf) %in% "n.rep"),]
  colnames(effects.ggdf) <- c("mean", "lower", "upper")
  contrast <- rownames(effects.ggdf)
  effects.ggdf$data <- as.numeric(data.effects)
  attr(effects.ggdf, "post.effects.samples") <- sim.effects
  effects.ggdf
  
}

# The below function picks certain parameters from a list called pickps_set,
# and replaces them with pickps_other, before performing posterior prediciton.
# We use it to turn control mechanisms off in the model. To turn off proactive
# control, we set the ongoing task thresholds equal in the PM block to the control
# thresholds. To turn off reactive, we set the ongiong rates on PM trials (PM block)
# to the ongoing rates on non-PM trials (PM block)

# samples=samples[[3]]
# probs=c(1:99)/100;random=TRUE
# bw="nrd0";report=10;save.simulation=TRUE;factors=NA; n.post=100
# pickps_set <- c("B.A2C",        "B.B2C",        "B.C2C" ,
#                 "B.D2C"  ,             "B.A2N"  ,      "B.B2N"  ,      "B.C2N"  ,
#                 "B.D2N")
#
# pickps_others <- c("B.A3C"   ,     "B.B3C",        "B.C3C" ,
#                   "B.D3C",        "B.A3N"   ,     "B.B3N" ,       "B.C3N" ,
#                   "B.D3N" )
#
#

pickps.post.predict.dmc = function(samples,n.post=100,probs=c(1:99)/100,random=TRUE,
                                   bw="nrd0",report=10,save.simulation=TRUE,factors=NA, pickps_others, pickps_set)
  # make list of posterior predictive density, quantiles and response probability
{
  
  get.dqp <- function(sim,facs,probs,n.post=NA) {
    
    quantile.names <- function(x,probs=seq(0, 1, 0.25),na.rm=FALSE,type=7, ...) {
      out <- quantile(x,probs=probs,na.rm=na.rm,type=type,names=FALSE,...)
      names(out) <- probs*100
      if (all(is.na(out))) NULL else out
    }
    
    qs <- tapply(sim$RT,sim[,c(facs,"R")],quantile.names,probs=probs,na.rm=TRUE)
    #     qs <- apply(qs,1:length(dim(qs)),function(x){
    #       if ( is.null(x[[1]]) || all(is.na(x[[1]]))) NULL else x[[1]]})
    
    # get probabilities given not na
    simOK <- sim[!is.na(sim$RT),]
    p <- tapply(simOK$RT,simOK[,c(facs,"R")],length)
    p[is.na(p)] <- 0 # In case some cells are empty
    np <- rep(apply(p,1:length(facs),sum),times=length(levels(simOK$R)))
    p <- p/np
    # get p NA
    pNA <- tapply(sim$RT,sim[,c(facs,"R")],length)
    pNA[is.na(pNA)] <- 0 # In case some cells are empty
    npNA <- rep(apply(pNA,1:length(facs),sum),times=length(levels(sim$R)))
    pNA <- tapply(is.na(sim$RT),sim[,c(facs,"R")],sum)
    pNA[is.na(pNA)] <- 0 # In case some cells are empty
    pNA <- pNA/npNA
    
    # For a simulation get proability replicates
    if (!is.na(n.post)) {
      repfac <- rep(1:n.post,each=sum(ns))
      repfac <- repfac[!is.na(sim$RT)]
      ps <- tapply(simOK$RT,cbind(simOK[,c(facs,"R")],rep=repfac),length)
      ps[is.na(ps)] <- 0 # In case some cells are empty
      ps <- n.post*ps/np
      # and NA replicates
      pNAs <- tapply(is.na(sim$RT),cbind(sim[,c(facs,"R")],rep=rep(1:n.post,each=sum(ns))),sum)
      pNAs[is.na(pNAs)] <- 0 # In case some cells are empty
      pNAs <- n.post*pNAs/npNA
    } else {
      ps=NULL
      pNAs=NULL
    }
    
    # cell names
    cell.names <- dimnames(qs)[[1]]
    n.cell <- length(facs)+1
    for ( i in 2:n.cell )
      cell.names <- outer(cell.names,dimnames(qs)[[i]],"paste",sep=".")
    cell.names <- as.vector(cell.names)
    # Get density and make defective
    dens <- tapply(sim$RT,sim[,c(facs,"R")],function(x){
      if (all(is.na(x))) NULL else {
        x <- x[x>=quantile(x,.01,na.rm=TRUE) & x<=quantile(x,.99,na.rm=TRUE)]
        if (length(x[!is.na(x)])<2) NULL else density(x[!is.na(x)],bw=bw)
      }
    })
    for (i in 1:length(p)) if ( is.finite(p[i]) && !(p[i]==0) )
    {
      if (!is.null(qs[i][[1]])) {
        names(qs[i][[1]]) <- as.numeric(names(qs[i][[1]]))*p[i]/100
        attr(qs[i][[1]],"cell.name") <- cell.names[i]
      }
      if (!is.null(dens[i][[1]]) ) {
        dens[i][[1]]$y <- dens[i][[1]]$y*p[i]
        attr(dens[i][[1]],"cell.name") <- cell.names[i]
      }
    }
    dnd <- dimnames(dens)
    dnq <- dimnames(qs)
    dens <- apply(dens,1:length(facs),function(x){x})
    qs <- apply(qs,1:length(facs),function(x){x})
    if ( is.null(dim(dens)) ) {
      dens <- array(dens,dim=c(length(dens)))
      dimnames(dens) <- dnd[-length(dnd)]
      qs <- array(qs,dim=c(length(qs)))
      dimnames(qs) <- dnq[-length(dnq)]
    }
    list(pdf=dens,cdf=qs,p=p,ps=ps,pNA=pNA,pNAs=pNAs)
  }
  
  model <- attributes(samples$data)$model
  facs <- names(attr(model,"factors"))
  if (any(is.na(factors))) factors <- facs
  if (!all(factors %in% facs))
    stop(paste("Factors argument must contain one or more of:",paste(facs,collapse=",")))
  resp <- names(attr(model,"responses"))
  ns <- table(samples$data[,facs],dnn=facs)
  n.par <- dim(samples$theta)[2]
  thetas <- matrix(aperm(samples$theta,c(3,1,2)),ncol=n.par)
  colnames(thetas) <- dimnames(samples$theta)[[2]]
  if (is.na(n.post)) use <- c(1:dim(thetas)[1]) else {
    if (random) use <- sample(c(1:dim(thetas)[1]),n.post,replace=F) else
      use <- round(seq(1,dim(thetas)[1],length.out=n.post))
  }
  n.post <- length(use)
  posts <- thetas[use,]
  
  ### Replace some parameter values with others
  posts[,colnames(posts) %in% pickps_others][,pickps_others] <- 
    posts[,colnames(posts) %in% pickps_set][,pickps_set] 
  
  ###
  n.rep <- sum(ns)
  sim <- data.frame(matrix(nrow=n.post*n.rep,ncol=dim(samples$data)[2]))
  names(sim) <- names(samples$data)
  # Tweaks for Stop Signal
  if ( !any(names(samples$data)=="SSD") ) {
    SSD <- rep(Inf,sum(ns))
    leave.out <- -c(1:dim(samples$data)[2])[names(samples$data)=="RT"]
  } else {
    # Assumes last two are SSD and RT! FIX ME.
    if (is.null(facs)) SSD <- samples$data$SSD else
      SSD <- unlist(tapply(samples$data$SSD,samples$data[,facs],identity))
    leave.out <- -c((dim(samples$data)[2]-1):dim(samples$data)[2])
  }
  cat("\n")
  cat(paste("Simulating (\'.\'=",report,"): ",sep=""))
  for (i in names(samples$data)[leave.out])
    sim[[i]] <- factor(rep(NA,n.post*n.rep),levels=levels(samples$data[[i]]))
  for (i in 1:n.post) {
    tmp <- simulate.dmc(posts[i,],model,ns,SSD=SSD)
    sim[(1+(i-1)*n.rep):(i*n.rep),names(tmp)] <- tmp
    if ( (i %% report) == 0) cat(".")
  }
  if ( any(names(sim)=="R2") ) { # MTR model
    sim$R <- paste(as.character(sim$R),as.character(sim$R2),sep="")
    sim$R[sim$R2=="DK"] <- "DK"
    sim$R <- factor(sim$R)
    samples$data$R <- paste(as.character(samples$data$R),as.character(samples$data$R2),sep="")
    samples$data$R[samples$data$R2=="DK"] <- "DK"
    samples$data$R <- factor(samples$data$R)
  }
  if (save.simulation) {
    sim <- cbind(reps=rep(1:n.post,each=dim(samples$data)[1]),sim)
    attr(sim,"data") <- samples$data
    sim
  } else {
    sim.dqp <- get.dqp(sim,factors,probs,n.post)
    dat.dqp <- get.dqp(sim=samples$data,factors,probs)
    names(dat.dqp) <- paste("data",names(dat.dqp),sep=".")
    c(sim.dqp,dat.dqp)
  }
}

# lapply the above to the whole samples object.
pickps.h.post.predict.dmc <- function(samples, n.post=100, probs=c(1:99)/100,
                                      bw="nrd0",
                                      save.simulation=FALSE, pickps_set, pickps_others)
  # apply lost.predict to each subject
{
  lapply(samples, pickps.post.predict.dmc, n.post=n.post, probs=probs, bw=bw,
         save.simulation=save.simulation, pickps_set=pickps_set, pickps_others=pickps_others)
}


# Accuracy effects (ongoing task)
RP_effects_ongoing <- function (data) {
  # In: simulated or observed choice-RT data
  # Out: response proportions associated with various data effects
  
  # Overall accuracy
  total <- length(data$R[ ((data$S=="cc" & data$R=="C")|(data$S=="nn" & data$R=="N")|
                             (data$S=="pc" & data$R=="P")|(data$S=="pn" & data$R=="P")) ]) / length(data$R)
  
  # Ongoing task accuracy
  ongoing <- length(data$R[ ((data$S=="cc" & data$R=="C")|(data$S=="nn" & data$R=="N")) ]) /
    length(data$R[ ((data$S=="cc")|(data$S=="nn")) ])
  
  # Stimulus
  nonconflict <- length(data$R[ (data$S == "nn" & data$R == "N") ]) /
    length(data$R[ (data$S == "nn") ])

  conflict <- length(data$R[ (data$S == "cc" & data$R == "C") ]) /
    length(data$R[ (data$S == "cc") ])

  stimulus <- conflict - nonconflict
  
  # PM (high TP)
  PM10_TP3 <- length(data$R[ ((data$S=="cc" & data$R=="C")|(data$S=="nn" & data$R=="N")) &
                               (data$PM == "10" & data$TP == "3s")]) /
    length(data$R[ ((data$S=="cc")|(data$S=="nn")) & (data$PM == "10" & data$TP == "3s")])

  PM30_TP3 <- length(data$R[ ((data$S=="cc" & data$R=="C")|(data$S=="nn" & data$R=="N")) &
                               (data$PM == "30" & data$TP == "3s")]) /
    length(data$R[ ((data$S=="cc")|(data$S=="nn")) & (data$PM == "30" & data$TP == "3s")])

  PM_TP3 <- PM10_TP3 - PM30_TP3
  
  # PM (low TP)
  PM10_TP6 <- length(data$R[ ((data$S=="cc" & data$R=="C")|(data$S=="nn" & data$R=="N")) &
                               (data$PM == "10" & data$TP == "6s")]) /
    length(data$R[ ((data$S=="cc")|(data$S=="nn")) & (data$PM == "10" & data$TP == "6s")])

  PM30_TP6 <- length(data$R[ ((data$S=="cc" & data$R=="C")|(data$S=="nn" & data$R=="N")) &
                               (data$PM == "30" & data$TP == "6s")]) /
    length(data$R[ ((data$S=="cc")|(data$S=="nn")) & (data$PM == "30" & data$TP == "6s")])

  PM_TP6 <- PM10_TP6 - PM30_TP6
  
  # TP (PM 10%)
  TP6_PM10 <- length(data$R[ ((data$S=="cc" & data$R=="C")|(data$S=="nn" & data$R=="N")) &
                                (data$TP == "6s" & data$PM == "10")]) /
    length(data$R[ ((data$S=="cc")|(data$S=="nn")) & (data$TP == "6s" & data$PM == "10")])

  TP3_PM10 <- length(data$R[ ((data$S=="cc" & data$R=="C")|(data$S=="nn" & data$R=="N")) &
                                (data$TP == "3s" & data$PM == "10")]) /
    length(data$R[ ((data$S=="cc")|(data$S=="nn")) & (data$TP == "3s" & data$PM == "10")])

  TP_PM10 <- TP6_PM10 - TP3_PM10
  
  # TP (PM 30%)
  TP6_PM30 <- length(data$R[ ((data$S=="cc" & data$R=="C")|(data$S=="nn" & data$R=="N")) &
                               (data$TP == "6s" & data$PM == "30")]) /
    length(data$R[ ((data$S=="cc")|(data$S=="nn")) & (data$TP == "6s" & data$PM == "30")])
  
  TP3_PM30 <- length(data$R[ ((data$S=="cc" & data$R=="C")|(data$S=="nn" & data$R=="N")) &
                               (data$TP == "3s" & data$PM == "30")]) /
    length(data$R[ ((data$S=="cc")|(data$S=="nn")) & (data$TP == "3s" & data$PM == "30")])
  
  TP_PM30 <- TP6_PM30 - TP3_PM30
  
  # Output data frame  
  out <- data.frame(
    total = total,
    ongoing = ongoing,
    stimulus = stimulus,
    PM_TP3 = PM_TP3,
    PM_TP6 = PM_TP6,
    TP_PM10 = TP_PM10,
    TP_PM30 = TP_PM30
  )
  names(out) <- c(
    "Total", 
    "Ongoing", 
    "Stimulus", 
    "PM_TP3", 
    "PM_TP6",
    "TP_PM10",
    "TP_PM30"
  )
  out
  
}


# Accuracy effects (PM task)
RP_effects_PM <- function (data) {
  # In: simulated or observed choice-RT data
  # Out: response proportions associated with various data effects
  
  # PM task accuracy
  PM <- length(data$R[ ((data$S == "pc" & data$R == "P")|(data$S == "pn" & data$R == "P")) ]) /
    length(data$R[ ((data$S == "pc")|(data$S == "pn")) ])
  
  # Stimulus
  nonconflict <- length(data$R[ (data$S == "pn" & data$R == "P") ]) /
    length(data$R[ (data$S == "pn") ])
  
  conflict <- length(data$R[ (data$S == "pc" & data$R == "P") ]) /
    length(data$R[ (data$S == "pc") ])
  
  stimulus <- conflict - nonconflict
  
  # PM (high TP)
  PM10_TP3 <- length(data$R[ ((data$S == "pc" & data$R == "P")|(data$S == "pn" & data$R == "P")) &
                               (data$PM == "10" & data$TP == "3s") ]) /
    length(data$R[ ((data$S == "pc")|(data$S == "pn")) & (data$PM == "10" & data$TP == "3s") ])
  
  PM30_TP3 <- length(data$R[ ((data$S == "pc" & data$R == "P")|(data$S == "pn" & data$R == "P")) &
                               (data$PM == "30" & data$TP == "3s") ]) /
    length(data$R[ ((data$S == "pc")|(data$S == "pn")) & (data$PM == "30" & data$TP == "3s") ])
  
  PM_TP3 <- PM10_TP3 - PM30_TP3
  
  # PM (low TP)
  PM10_TP6 <- length(data$R[ ((data$S == "pc" & data$R == "P")|(data$S == "pn" & data$R == "P")) &
                               (data$PM == "10" & data$TP == "6s") ]) /
    length(data$R[ ((data$S == "pc")|(data$S == "pn")) & (data$PM == "10" & data$TP == "6s") ])
  
  PM30_TP6 <- length(data$R[ ((data$S == "pc" & data$R == "P")|(data$S == "pn" & data$R == "P")) &
                               (data$PM == "30" & data$TP == "6s") ]) /
    length(data$R[ ((data$S == "pc")|(data$S == "pn")) & (data$PM == "30" & data$TP == "6s") ])
  
  PM_TP6 <- PM10_TP6 - PM30_TP6
  
  # TP (PM 10%)
  TP6_PM10 <- length(data$R[ ((data$S == "pc" & data$R == "P")|(data$S == "pn" & data$R == "P")) &
                               (data$TP == "6s" & data$PM == "10") ]) /
    length(data$R[ ((data$S == "pc")|(data$S == "pn")) & (data$TP == "6s" & data$PM == "10") ])
  
  TP3_PM10 <- length(data$R[ ((data$S == "pc" & data$R == "P")|(data$S == "pn" & data$R == "P")) &
                               (data$TP == "3s" & data$PM == "10") ]) /
    length(data$R[ ((data$S == "pc")|(data$S == "pn")) & (data$TP == "3s" & data$PM == "10") ])
  
  TP_PM10 <- TP6_PM10 - TP3_PM10
  
  # TP (PM 30%)
  TP6_PM30 <- length(data$R[ ((data$S == "pc" & data$R == "P")|(data$S == "pn" & data$R == "P")) &
                               (data$TP == "6s" & data$PM == "30") ]) /
    length(data$R[ ((data$S == "pc")|(data$S == "pn")) & (data$TP == "6s" & data$PM == "30") ])
  
  TP3_PM30 <- length(data$R[ ((data$S == "pc" & data$R == "P")|(data$S == "pn" & data$R == "P")) &
                               (data$TP == "3s" & data$PM == "30") ]) /
    length(data$R[ ((data$S == "pc")|(data$S == "pn")) & (data$TP == "3s" & data$PM == "30") ])
  
  TP_PM30 <- TP6_PM30 - TP3_PM30
  
  # Output data frame  
  out <- data.frame(
    PM = PM,
    stimulus = stimulus,
    PM_TP3 = PM_TP3,
    PM_TP6 = PM_TP6,
    TP_PM10 = TP_PM10,
    TP_PM30 = TP_PM30
  )
  names(out) <- c(
    "PM",
    "Stimulus", 
    "PM_TP3", 
    "PM_TP6",
    "TP_PM10", 
    "TP_PM30"
  )
  out
  
}


# RT effects (ongoing task)
RT_effects_ongoing <- function (data) {
  # In: simulated or observed choice-RT data
  # Out: mean response times associated with various data effects
  
  # Correct trial RT
  correct <- mean(data$RT[ ((data$S == "cc" & data$R == "C")|
                              (data$S == "nn" & data$R == "N")) ])
  
  # Error trial RT
  error <- mean(data$RT[ ((data$S == "cc" & data$R == "N")|
                            (data$S == "nn" & data$R == "C")) ])
  
  # Stimulus
  conflict <- mean(data$RT[ (data$S == "cc" & data$R == "C") ])
  
  nonconflict <- mean(data$RT[ (data$S == "nn" & data$R == "N") ])
  
  stimulus <- conflict - nonconflict
  
  # PM (high TP)
  PM10_TP3 <- mean(data$RT[ ((data$S == "cc" & data$R == "C")|
                               (data$S == "nn" & data$R == "N")) &
                              (data$PM == "10" & data$TP == "3s") ])
  
  PM30_TP3 <- mean(data$RT[ ((data$S == "cc" & data$R == "C")|
                               (data$S == "nn" & data$R == "N")) &
                              (data$PM == "30" & data$TP == "3s") ])
  
  PM_TP3 <- PM30_TP3 - PM10_TP3
  
  # PM (low TP)
  PM10_TP6 <- mean(data$RT[ ((data$S == "cc" & data$R == "C")|
                               (data$S == "nn" & data$R == "N")) &
                              (data$PM == "10" & data$TP == "6s") ])
  
  PM30_TP6 <- mean(data$RT[ ((data$S == "cc" & data$R == "C")|
                               (data$S == "nn" & data$R == "N")) &
                              (data$PM == "30" & data$TP == "6s") ])
  
  PM_TP6 <- PM30_TP6 - PM10_TP6
  
  # TP (PM 10%)
  TP3_PM10 <- mean(data$RT[ ((data$S == "cc" & data$R == "C")|
                               (data$S == "nn" & data$R == "N")) &
                              (data$PM == "10" & data$TP == "3s") ])
  
  TP6_PM10 <- mean(data$RT[ ((data$S == "cc" & data$R == "C")|
                               (data$S == "nn" & data$R == "N")) &
                              (data$PM == "10" & data$TP == "6s") ])
  
  TP_PM10 <- TP6_PM10 - TP3_PM10
  
  # TP (PM 30%)
  TP3_PM30 <- mean(data$RT[ ((data$S == "cc" & data$R == "C")|
                               (data$S == "nn" & data$R == "N")) &
                              (data$PM == "30" & data$TP == "3s") ])
  
  TP6_PM30 <- mean(data$RT[ ((data$S == "cc" & data$R == "C")|
                               (data$S == "nn" & data$R == "N")) &
                              (data$PM == "30" & data$TP == "6s") ])
  
  TP_PM30 <- TP6_PM30 - TP3_PM30
  
  # Output data frame  
  out <- data.frame(
    correct = correct,
    error = error,
    stimulus = stimulus,
    PM_TP3 = PM_TP3,
    PM_TP6 = PM_TP6,
    TP_PM10 = TP_PM10,
    TP_PM30 = TP_PM30
  )
  names(out) <- c(
    "Correct", 
    "Error", 
    "Stimulus", 
    "PM_TP3", 
    "PM_TP6",
    "TP_PM10", 
    "TP_PM30"
  )
  out
  
}


# RT effects (PM task)
RT_effects_PM <- function (data) {
  # In: simulated or observed choice-RT data
  # Out: mean response times associated with various data effects
  
  # Correct trial RT
  correct <- mean(data$RT[ ((data$S == "pc" & data$R == "P")|
                              (data$S == "pn" & data$R == "P")) ])
  
  # PM miss
  PM_miss <- mean(data$RT[ ((data$S == "pc" & data$R == "C")|
                            (data$S == "pn" & data$R == "N")) ])
  
  # PM false alarm
  PM_FA <- mean(data$RT[ ((data$S == "cc" & data$R == "P")|
                              (data$S == "nn" & data$R == "P")) ])
  
  # Stimulus
  conflict <- mean(data$RT[ (data$S == "pc" & data$R == "P") ])
  
  nonconflict <- mean(data$RT[ (data$S == "pn" & data$R == "P") ])
  
  stimulus <- conflict - nonconflict
  
  # PM (high TP)
  PM10_TP3 <- mean(data$RT[ ((data$S == "pc" & data$R == "P")|
                               (data$S == "pn" & data$R == "P")) &
                              (data$PM == "10" & data$TP == "3s") ])
  
  PM30_TP3 <- mean(data$RT[ ((data$S == "pc" & data$R == "P")|
                               (data$S == "pn" & data$R == "P")) &
                              (data$PM == "30" & data$TP == "3s") ])
  
  PM_TP3 <- PM30_TP3 - PM10_TP3
  
  # PM (low TP)
  PM10_TP6 <- mean(data$RT[ ((data$S == "pc" & data$R == "P")|
                               (data$S == "pn" & data$R == "P")) &
                              (data$PM == "10" & data$TP == "6s") ])
  
  PM30_TP6 <- mean(data$RT[ ((data$S == "pc" & data$R == "P")|
                               (data$S == "pn" & data$R == "P")) &
                              (data$PM == "30" & data$TP == "6s") ])
  
  PM_TP6 <- PM30_TP6 - PM10_TP6
  
  # TP (PM 10%)
  TP3_PM10 <- mean(data$RT[ ((data$S == "pc" & data$R == "P")|
                               (data$S == "pn" & data$R == "P")) &
                              (data$PM == "10" & data$TP == "3s") ])
  
  TP6_PM10 <- mean(data$RT[ ((data$S == "pc" & data$R == "P")|
                               (data$S == "pn" & data$R == "P")) &
                              (data$PM == "10" & data$TP == "6s") ])
  
  TP_PM10 <- TP6_PM10 - TP3_PM10
  
  # TP (PM 30%)
  TP3_PM30 <- mean(data$RT[ ((data$S == "pc" & data$R == "P")|
                               (data$S == "pn" & data$R == "P")) &
                              (data$PM == "30" & data$TP == "3s") ])
  
  TP6_PM30 <- mean(data$RT[ ((data$S == "pc" & data$R == "P")|
                               (data$S == "pn" & data$R == "P")) &
                              (data$PM == "30" & data$TP == "6s") ])
  
  TP_PM30 <- TP6_PM30 - TP3_PM30
  
  # Output data frame  
  out <- data.frame(
    correct = correct,
    PM_miss = PM_miss,
    PM_FA = PM_FA,
    stimulus = stimulus,
    PM_TP3 = PM_TP3,
    PM_TP6 = PM_TP6,
    TP_PM10 = TP_PM10,
    TP_PM30 = TP_PM30
  )
  names(out) <- c(
    "Correct", 
    "PM miss", 
    "PM FA",
    "Stimulus", 
    "PM_TP3", 
    "PM_TP6",
    "TP_PM10", 
    "TP_PM30"
  )
  out
  
}


# # Function to calculate desired RT effects
# RT.median <- function (data) {
#   # In: simulated or observed choice-RT data
#   # Out: median response times associated with various data effects
#   
#   # Correct trial RT
#   correct <- median(data$RT[ ((data$C == "sam" & data$R == "SAM")|
#                                 (data$C == "dif" & data$R == "DIF" )) ])
#   
#   # Error trial RT
#   error <- median(data$RT[ ((data$C == "sam" & data$R == "DIF")|
#                               (data$C == "dif" & data$R == "SAM" )) ])
#   
#   
#   # Updating
#   # Reference trial RT
#   ref_nos <- median(data$RT[ ((data$U == "ref") & (data$S == "nos")) & 
#                                ((data$C == "sam" & data$R == "SAM")|
#                                   (data$C == "dif" & data$R == "DIF" )) ]) 
#   # Comparison trial RT
#   cmp_nos <- median(data$RT[ ((data$U == "cmp") & (data$S == "nos")) & 
#                                ((data$C == "sam" & data$R == "SAM")|
#                                   (data$C == "dif" & data$R == "DIF" )) ]) 
#   # Updating cost
#   updating <- ref_nos - cmp_nos
#   
#   
#   # Mode switching
#   # Switch trial RT
#   swi <- median(data$RT[ (data$S == "swi") & ((data$C == "sam" & data$R == "SAM")|
#                                                 (data$C == "dif" & data$R == "DIF" )) ])
#   
#   # No-switch trial RT
#   nos <- median(data$RT[ (data$S == "nos") & ((data$C == "sam" & data$R == "SAM")|
#                                                 (data$C == "dif" & data$R == "DIF" )) ]) 
#   
#   # Mode-switching cost
#   switching <- swi - nos
#   
#   
#   # Comparison
#   # Different trial RT
#   dif_nos <- median(data$RT[ (data$S == "nos") & (data$C == "dif" & data$R == "DIF") ]) 
#   
#   # Same trial RT
#   sam_nos <- median(data$RT[ (data$S == "nos") & (data$C == "sam" & data$R == "SAM") ]) 
#   
#   # Comparison cost
#   comparison <- dif_nos - sam_nos
#   
#   
#   # Gate opening
#   # Reference/switch trial RT
#   ref_swi <- median(data$RT[ ((data$U == "ref") & (data$S == "swi")) & 
#                                ((data$C == "sam" & data$R == "SAM")|
#                                   (data$C == "dif" & data$R == "DIF" )) ]) 
#   
#   # Reference/no-switch trial RT
#   ref_nos <- median(data$RT[ ((data$U == "ref") & (data$S == "nos")) & 
#                                ((data$C == "sam" & data$R == "SAM")|
#                                   (data$C == "dif" & data$R == "DIF" )) ]) 
#   
#   # Gate-opening cost
#   gate_opening <- ref_swi - ref_nos
#   
#   
#   # Gate closing
#   # Comparison/switch trial RT
#   cmp_swi <- median(data$RT[ ((data$U == "cmp") & (data$S == "swi")) & 
#                                ((data$C == "sam" & data$R == "SAM")|
#                                   (data$C == "dif" & data$R == "DIF" )) ]) 
#   
#   # Comparison/no-switch trial RT
#   cmp_nos <- median(data$RT[ ((data$U == "cmp") & (data$S == "nos")) & 
#                                ((data$C == "sam" & data$R == "SAM")|
#                                   (data$C == "dif" & data$R == "DIF" )) ]) 
#   
#   # Gate-closing cost
#   gate_closing <- cmp_swi - cmp_nos
#   
#   
#   # Substitution
#   # Comparison/different trial RT
#   cmp_dif <- median(data$RT[ ((data$U == "cmp") & (data$C == "dif" & data$R == "DIF")) ]) 
#   
#   # Comparison/same trial RT
#   cmp_sam <- median(data$RT[ ((data$U == "cmp") & (data$C == "sam" & data$R == "SAM")) ]) 
#   
#   # Reference/different trial RT
#   ref_dif <- median(data$RT[ ((data$U == "ref") & (data$C == "dif" & data$R == "DIF")) ]) 
#   
#   # Reference/same trial RT
#   ref_sam <- median(data$RT[ ((data$U == "ref") & (data$C == "sam" & data$R == "SAM")) ]) 
#   
#   # Substitution cost
#   substitution <- (ref_dif - ref_sam) - (cmp_dif - cmp_sam) 
#   
#   
#   # Output data frame  
#   out <- data.frame(correct = correct,
#                     error = error,
#                     updating = updating,
#                     switching = switching,
#                     comparison = comparison,
#                     gate_opening = gate_opening,
#                     gate_closing = gate_closing,
#                     substitution = substitution
#   )
#   names(out) <- c("Correct", "Error", "Updating", "Switching",
#                   "Comparison", "Gate opening", "Gate closing",
#                   "Substitution")
#   out
#   
# }

