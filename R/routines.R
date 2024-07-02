KM_resolve <- function(S, t.risk, n.risk, ndeath=NULL,ticks=NULL,strict_tick=FALSE,strict_dec=TRUE,totaltime=NULL,totaltime_power=1, c.event=NULL, t.event=NULL, optmethod="approx",cen_penalty=1e-3,constr_tol=1e-8,nprobe=0,epagap=1e-6,epgap=1e-4,tilim=1000,trace=1,ttol=0.01,cen_max=NULL,ceventinc=FALSE) {
  #Wrapper function that just uses the CIF method (should be equivalent in the case of one risk)
  #S: a list containing $time: times KM extracted, $surv: KM values at those times
  #t.risk = vector of times at which number at risk is known
  #n.risk = vector of numbers at risk
  #ndeath = total number of events/deaths (can be omitted)
  #ticks = optional vector of locations of tick marks representing censoring times
  #strict_tick = whether to enforce a constraint that censoring only occurs at the identified tick points
  #strict_dec = whether to enforce a constraint that any identified decrement in the curve must be associated with at least 1 observed event
  #totaltime = observed/inferred total time at risk or total transformed times at risk (sum t_{i}^pow)
  #totaltime_power = the power by which each of the individual times has been raised in the calculation of totaltime, defaults to 1, potentially useful if the MLEs of the Weibull rate and shape parameters are known.
  #c.event = optional vector of cumulative number of events at the t.event values
  #t.event = optional vector of times at which number of cumulative events is supplied, if NULL makes t.event=t.risk
  #optmethod = Method of finding the solution: "approx" - solve the continuous QP and then do integer rounding or "miqp" - use full mixed integer QP optimization using Rcplex
  #cen_penalty = Penalty on the square of the objective matrix (penalizes solutions with lots of ties in censoring)
  #constr_tol = Tolerance value on constraints to avoid solve.QP erroneously claiming no solutions (increase the value if get this error)
  #nprobe, epgap, tilim = control parameters for MIQP method (tilim is maximum computation time for optimizer in seconds)
  #cen_max = Maximum number of censoring times at a particular point in the solution (another way of preventing solutions with over clumpy censoring distributions)
  #ttol = percentage tolerance by which the totaltime value is to be recovered
  #ceventinc = if c.event supplied, whether the cumulative total is events at times <=t (rather than <t if FALSE).
  if (is.data.frame(S)) S <- list(time=S$time,surv=S$surv)
  if (!is.null(ndeath) & !is.null(c.event)) {
    if (max(c.event) > ndeath) {
       stop("Cumulative events in c.event contradict total events in ndeath.")
    }
  }
  S$cif <- cbind(1- S$surv)
  output <- CIF_resolve(S, t.risk, n.risk, nevent=NULL,ndeath,ticks,strict_tick,strict_dec, totaltime,totaltime_power,c.event, t.event, optmethod,cen_penalty,constr_tol,nprobe,epagap,epgap,tilim,trace,ttol,cen_max, ceventinc)
  #Translate this output into the standard output (need to retain attributes)
  int_obj <- attr(output,"int_obj")
  status <- attr(output,"object")$status
  names(output)[names(output)=="nevent.1"]<-"nevent"
  #output <- output[,-which(names(output)=="tevent")]
  attr(output,"int_obj") <- int_obj
  attr(output,"status") <- status
  class(output) <- "KMresolve"
  return(output)
}

CIF_resolve <- function(S, t.risk, n.risk, nevent=NULL,ndeath=NULL,ticks=NULL,strict_tick=FALSE,strict_dec=TRUE, totaltime=NULL, totaltime_power=1, c.event=NULL,t.event=NULL,optmethod="approx",cen_penalty=1e-3,constr_tol=1e-8,nprobe=0,epagap=1e-6,epgap=1e-4,tilim=1000,trace=1,ttol=0.01,cen_max=NULL,ceventinc=FALSE) {
  #S$cif should be a list containing $time: times at which the CIF extracted, $cif: a matrix of CIF values at those times.
  #t.risk = vector of times at which number at risk is known
  #n.risk = vector of numbers at risk
  #nevent = vector of number of events of each type (can be omitted)
  #ndeath = total number of events/deaths across all types (can be omitted)
  #ticks = optional vector of locations of tick marks representing censoring times
  #strict_tick = whether to enforce a constraint that censoring only occurs at the tick points
  #strict_dec = whether to enforce a constraint that any identified decrement (increment in estimated CIF) in the curve must be associated with at least 1 observed event
  #totaltime = observed/inferred total time at risk or total transformed times at risk (sum t_{i}^pow)
  #totaltime_power = the power by which each of the individual times has been raised in the calculation of totaltime, defaults to 1, potentially useful if the MLEs of the Weibull rate and shape parameters are known.
  #c.event = optional vector or matrix of cumulative number of events; should either be a vector (in which case is interpreted as the cumulative total events) or a matrix where each column represents a given risk
  #t.event = optional vector of times at which number of cumulative events is supplied, if NULL makes t.event=t.risk
  #optmethod = Method of finding the solution: "approx" - solve the continuous QP and then do integer rounding or "miqp" - use full mixed integer QP optimization using Rcplex
  #constr_tol = Tolerance value on constraints to avoid solve.QP erroneously claiming no solutions (increase the value if get this error)
  #nprobe, epgap, tilim = control parameters for MIQP method.
  #cen_max = Maximum number of censoring times at a particular point in the solution (another way of preventing solutions with over clumpy censoring distributions)
  #ttol = percentage tolerance by which the totaltime value is to be recovered
  #ceventinc = if c.event supplied, whether the cumulative total is events at times <=t (rather than <t if FALSE).
  if (!is.null(totaltime) & optmethod!="miqp") warning("Total time constraint may not be met in the integer solution when the approximate integer rounding method is used. Use MIQP, if required.")
  if (is.data.frame(S)) stop("S should be a list not a data frame")
  if (!"time"%in%names(S)) stop("S needs to include a component named time")
  if (!"cif"%in%names(S)) stop("S needs to include a component named cif")
  origS <- S
  if (!optmethod%in%c("approx","miqp")) stop("Optimization method must either be approx or miqp")
  #if (optmethod=="miqp") require("Rcplex")
  #if (optmethod=="approx") require("quadprog")
  if (length(unique(n.risk))!=length(n.risk)) {
    #Retain only the later time with the same numbers at risk.
    g <- sapply(unique(n.risk), function(x) max(which(n.risk==x)))
    n.risk <- n.risk[g]
    t.risk <- t.risk[g]
  }
  myapproxfun <- function(x,y, method) { suppressWarnings(approxfun(x,y,method=method))}
  #nevent: vector of length ncomp giving the total number of events of each type
  #ndeath: alternatively give just the total number of events of all types

  if (!is.null(nevent) & !is.null(ndeath)) {
    if (sum(nevent)!=ndeath) stop("nevent disagrees with ndeath: ndeath should give total number of events of all types.")
  }
  if (min(n.risk)>0) {
    #Impute an extra time with zero at risk
    cmax <- max(t.risk)
    t.risk <- c(t.risk,max(S$time,t.risk,ticks)+1)
    n.risk <- c(n.risk,0)
    if (!is.null(ticks)) {
      if (max(ticks) < cmax) {
        #Add an extra tick point
        ticks <- c(ticks,max(S$time,t.risk,ticks)+1)
      }
    }
  }
  last_time <- max(S$time,t.risk,ticks)

  ncomp <- dim(S$cif)[2]
  if (!is.null(nevent)) {
    if (min(nevent)<0) stop("Number of events cannot be negative!")
    if (max(abs(nevent - floor(nevent)))>1e-6) stop("nevent should be integers")
    if (length(nevent)!=ncomp) stop("nevent should give number of events of each type")
  }

  extraS <- list()
  for (i in 1:ncomp) {
    extraS[[i]]<-myapproxfun(c(0,S$time,last_time),c(0,S$cif[,i],max(S$cif[,i])),method="constant")(t.risk)
  }
  miss <- which(!t.risk%in%S$time)
  if (length(miss)>0) {
    S$time <-c(S$time,t.risk[miss])
    CIFmiss <- array(0,c(length(miss),ncomp))
    for (i in 1:ncomp) CIFmiss[,i] <- extraS[[i]][miss]
    S$cif <- rbind(S$cif, CIFmiss)
    S$cif <- S$cif[order(S$time),,drop=FALSE]
    S$time <- sort(S$time)
  }
  if (!is.null(ticks)) {
    extraS2 <- list()
    for (i in 1:ncomp) {
      extraS2[[i]]<-myapproxfun(c(0,S$time,max(ticks,S$time)),c(0,S$cif[,i],max(S$cif[,i])),method="constant")(ticks)
    }
    miss2 <- which(!ticks%in%S$time)
    if (length(miss2)>0) {
      S$time <-c(S$time,ticks[miss2])
      CIFmiss2 <- array(0,c(length(miss2),ncomp))
      for (i in 1:ncomp) CIFmiss2[,i] <- extraS2[[i]][miss2]
      S$cif <-rbind(S$cif, CIFmiss2)
      S$cif <- S$cif[order(S$time),,drop=FALSE]
      S$time <- sort(S$time)
    }
  }

  #For number at risk assume is number at risk at time t (i.e. before events occurred)
  S$chunk <- sapply(S$time,function(x) sum(t.risk <= x))
  #For cumulative events assume number includes those that happened at time t
  if (last_time%in%ticks) {
    if (S$chunk[length(S$chunk)] > S$chunk[(length(S$chunk)-1)]) S$chunk[length(S$chunk)] <- S$chunk[(length(S$chunk)-1)]
  }
  chunk_total <- -diff(c(n.risk,0))
  if (!is.null(c.event)) {
    if (is.null(t.event)) t.event <- t.risk

    #Set up the chunks for cumulative events
    if (identical(t.event, t.risk) & !ceventinc) {
      S$chunk2 <- S$chunk
    }else{
      if (ceventinc) {
        S$chunk2 <- sapply(S$time,function(x) sum(t.event < x))
      }else{
        S$chunk2 <- sapply(S$time,function(x) sum(t.event <= x))
      }
    }
    if (is.vector(c.event)) {
       chunk_event_totals <- NULL
       chunk_event_total <- diff(c(c.event,c.event[length(c.event)]))
    }else{
       chunk_event_total <- NULL
       chunk_event_totals <- apply(c.event, 2, function(x) diff(c(x,x[length(x)])))
    }
  }else{
   chunk_event_total <- chunk_event_totals <- NULL
  }

  #Need to do something about chunks of zero length
  if (max(S$chunk) < length(chunk_total)) {
    chunk_total <- chunk_total[1:max(S$chunk)]
  }

  if (!is.null(chunk_event_total)) {
    if (max(S$chunk2) < length(chunk_event_total)) {
       chunk_event_total <- chunk_event_total[1:max(S$chunk2)]
    }
  }
  if (!is.null(chunk_event_totals)) {
    if (max(S$chunk2) < dim(chunk_event_totals)[1]) {
       chunk_event_totals <- chunk_event_totals[1:max(S$chunk2),]
    }
  }

  S$surv <- 1 - apply(S$cif,1,sum)

  Obs <- (S$cif - rbind(rep(0,ncomp),S$cif[1:c(dim(S$cif)[1]-1),,drop=FALSE]))/array(c(1,S$surv[1:(length(S$surv)-1)]),dim=dim(S$cif))
  Obs <- replace(Obs,which(is.nan(Obs)),0) #Remove NaN's which should be due to 0s going to 0s.
  decs <- which(Obs >0)
  #Set up the objective matrix
  B <- array(0,c(ncomp*length(S$time),(ncomp+1)*length(S$time)))
  for (j in 1:ncomp) {
    B[cbind((1:length(S$time) + (j-1)*length(S$time)),(1:length(S$time) + (j-1)*length(S$time)))]<-1
    for (i in 2:length(S$time)) {
      for (k in 1:(ncomp+1)) {
        B[i + (j-1)*length(S$time),(1:(i-1) + (k-1)*length(S$time))]<- Obs[i,j]
      }
    }
  }
  Dmat <- t(B)%*%B
  dvec <- n.risk[1]* t(c(Obs))%*%B

  event_totals <- nevent
  event_total <- ndeath

  #Five possibilities
  #i) No event info: nchunk
  #ii) Just total deaths: nchunk + 1
  #iii) Total of each event: nchunk + ncomp
  #iv) Chunked total deaths: nchunk + nchunk2
  #v) Chunked total for each event: nchunk + (ncomp)*nchunk2
  chunkscenario <- 1*(is.null(chunk_event_total) & is.null(chunk_event_totals) & is.null(event_totals) & is.null(event_total)) +2*(is.null(chunk_event_total) & is.null(chunk_event_totals) & is.null(event_totals) & !is.null(event_total))+3*(is.null(chunk_event_total) & is.null(chunk_event_totals) & !is.null(event_totals))+4*(!is.null(chunk_event_total) & is.null(chunk_event_totals)) + 5*(!is.null(chunk_event_totals))

  #Set up the constraint matrix:

  #Each chunk needs two constraints
  nchunk <- length(unique(S$chunk))
  if (chunkscenario%in%c(4,5)) {
    nchunk2 <- max(S$chunk2)
  }else{
    nchunk2 <- 0
  }
  nconstrE <- c(nchunk, nchunk+1,nchunk+ncomp,nchunk+nchunk2,nchunk + (ncomp)*nchunk2)[chunkscenario]

  Aeq <- array(0,c(nconstrE,(ncomp+1)*length(S$time)))
  Beq <-rep(0,nconstrE)
  for (i in 1:nchunk) {
    mm <- which(S$chunk==i)
    for (j in 1:ncomp) {
      mm <- c(mm,which(S$chunk==i)+j*length(S$time))
    }
    Aeq[i,mm] <- 1
    Beq[i] <- chunk_total[i]
  }

  if (chunkscenario==2) {

    for (j in 1:ncomp) {
      Aeq[nchunk+1,((1:length(S$time))+(j-1)*length(S$time))]<-1
    }
    Beq[nchunk+1] <- event_total
  }
  if (chunkscenario==3) {
    for (j in 1:ncomp) {
      Aeq[nchunk+j,((1:length(S$time))+(j-1)*length(S$time))]<-1
      Beq[nchunk+j]<-event_totals[j]
    }
  }
  if (chunkscenario==4) {

    for (i in 1:nchunk2) {
      mm <- which(S$chunk2==i)
      if (ncomp >1) {
       for (j in 1:(ncomp-1)) {
         mm <- c(mm,which(S$chunk2==i)+j*length(S$time))
       }
      }
      Aeq[nchunk+i,mm] <- 1
      Beq[nchunk+i] <- chunk_event_total[i]
    }

  }
  if (chunkscenario==5) {

   for (j in 1:ncomp) {
    for (i in 1:nchunk2) {
      mm <- which(S$chunk2==i) + (j-1)*length(S$time)
      Aeq[(nchunk+i+(j-1)*nchunk2),mm] <- 1
      Beq[(nchunk+i+(j-1)*nchunk2)] <- chunk_event_totals[i,j]
    }
   }
  }

   #Only penalize the censoring times
  Dmat2 <- Dmat + diag(c(rep(c(0,cen_penalty),c(ncomp*length(S$time),length(S$time)))))

  excl <- NULL
  for (j in 1:ncomp) {
    wvalj <- which(Obs[,j] == 0)
    if (j==1) {
      Aeq4 <- array(0,c(length(wvalj),(ncomp+1)*length(S$time)))
      Aeq4[cbind(1:length(wvalj),wvalj)]<-1
    }else{
      Aj <- array(0,c(length(wvalj),(ncomp+1)*length(S$time)))
      Aj[cbind(1:length(wvalj),wvalj + (j-1)*length(S$time))]<-1
      Aeq4 <- rbind(Aeq4,Aj)
    }
    excl <- c(excl, wvalj + (j-1)*length(S$time))
  }
  ncon1 <- dim(Aeq4)[1]




  if (!is.null(ticks)) {
    wval2 <- which(!S$time %in% ticks)
    Aeq5 <- array(0,c(length(wval2),(ncomp+1)*length(S$time)))
    Aeq5[cbind(1:length(wval2),(wval2 + ncomp*length(S$time)))]<-1

    #Is it really necessary to add in the 0 constraints??
    Aeq2 <- diag((ncomp+1)*length(S$time))
    Beq2 <- rep(0,dim(Aeq2)[1])
    ####################################################
    if (strict_tick) {
      wval3 <- which(S$time %in% ticks)
      Beq2[(wval3 + (ncomp)*length(S$time))] <- 1-constr_tol #avoid method thinking no solution
    }
    excl <- c(excl, wval2 + ncomp*length(S$time))
    decI <- 1*(strict_dec)*(1:((ncomp+1)*length(S$time)) %in% decs)[-excl]
    decS <- pmax(decI*(1-constr_tol), Beq2[-excl])
    Aeq2 <- Aeq2[-excl,]
    Aeq3 <- rbind(Aeq,Aeq2)
    Aeq6 <- rbind(Aeq4,Aeq5,Aeq3)
    #Beq3 <- c(rep(0,ncon1+length(wval2)),Beq,Beq2[-excl])
    Beq3 <- c(rep(0,ncon1+length(wval2)),Beq,decS)
  }else{
    Aeq2 <- diag((ncomp+1)*length(S$time))
    decI <- 1*(strict_dec)*(1:((ncomp+1)*length(S$time)) %in% decs)[-excl]
    Aeq2 <- Aeq2[-excl,]
    wval2 <- NULL
    Aeq3 <- rbind(Aeq,Aeq2)
    Aeq6 <-   rbind(Aeq4,Aeq3)
    Beq3 <- c(rep(0,ncon1),Beq,decI*(1-constr_tol)) #Force any decrement to be associated with at least one event
  }


  if (!is.null(totaltime)) {
    if (!is.numeric(totaltime)) stop("totaltime needs to be a numeric scalar corresponding to the total patient time at risk")
   #Need to add two more inequality constraints
    #If ticks supplied assume censoring is at the times themselves.
    #Otherwise take midpoint between the time and the next time
    if (is.null(ticks)) {
      cenpoints <- (S$time + c(S$time,S$time[length(S$time)])[2:(length(S$time)+1)])*0.5
    }else{
      cenpoints <- S$time
    }
    Aeq7 <- rbind(c(rep(S$time^totaltime_power, ncomp),cenpoints^totaltime_power),c(-rep(S$time^totaltime_power, ncomp),-cenpoints^totaltime_power))
    Aeq6 <- rbind(Aeq6,Aeq7)
    Beq3 <- c(Beq3, totaltime*(1-ttol), -totaltime*(1+ttol))
  }
  if (!is.null(cen_max)) {
    Aeq8 <- array(0,c(length(S$time)*(ncomp + 1), length(S$time)))
    Aeq8[cbind((ncomp*length(S$time) +1):((ncomp+1)*length(S$time)),1:length(S$time))] <- 1
    Aeq6 <- rbind(Aeq6,Aeq8)
    Beq3 <- c(Beq3, rep(cen_max, length(S$time)))
  }


  mats <- list(dvec=dvec, Q=Dmat2, A=Aeq6, bvec=Beq3)

  if (optmethod=="miqp") {
    qpobj <- Rcplex::Rcplex(cvec=-c(dvec),Amat=-(Aeq6),bvec=-Beq3,Qmat=Dmat2,sense=rep(c("E","L"),c(ncon1 + length(wval2) + dim(Aeq)[1], length(Beq3) - ncon1 - length(wval2) - dim(Aeq)[1])),lb=0,vtype="I",control=list(round=1,probe=as.integer(nprobe),epagap=epagap,epgap=epgap,tilim=tilim,trace=trace))
    #Need to add something that checks whether optimization took place i.e. status code.

    solution <- qpobj$xopt
    nevent <- array(0,c(length(S$time),ncomp))
    for (j in 1:ncomp) {
      nevent[,j] <- solution[(1:length(S$time) + (j-1)*length(S$time))]
    }
    intcenevent <- solution[(1+ncomp*length(S$time)):((ncomp+1)*length(S$time))]
    int_obj <- qpobj$obj
  }
  if (optmethod=="approx") {
    qpobj <- tryCatch(quadprog::solve.QP(Dmat=Dmat2,dvec=t(dvec),Amat = t(Aeq6),bvec = Beq3,meq=ncon1 + length(wval2) + dim(Aeq)[1]),error=function(e) return(NA))
    if (!is.list(qpobj)) {
      #Attempt a scaled version
      sc <- norm(Dmat2,"2")
      qpobj <- tryCatch(quadprog::solve.QP(Dmat=Dmat2/sc,dvec=t(dvec)/sc,Amat = t(Aeq6),bvec = Beq3,meq=ncon1 + length(wval2) + dim(Aeq)[1]),error=function(e) return(NA))
      if (!is.list(qpobj)) stop("Unable to find a feasible solution of QP. Check input data or try increasing constr_tol value.")
    }

    solution <- qpobj$solution
    nevent <- array(0,c(length(S$time),ncomp))
    for (j in 1:ncomp) {
      nevent[,j] <- diff(c(0,floor(round(cumsum(solution[(1:length(S$time) + (j-1)*length(S$time))]),5)+0.5)))
    }
    tevent <- apply(nevent,1,sum)
    cenevent <- pmax(0,solution[(1+ncomp*length(S$time)):((ncomp+1)*length(S$time))])
    cen_needed <- chunk_total - tapply(tevent, factor(S$chunk),sum)
    cur_cen <-  round(tapply(cenevent, S$chunk,sum),8)
    newcenevent <- cenevent * rep(cen_needed,table(S$chunk))/rep(cur_cen + 1*(cur_cen==0),table(S$chunk))
    intcenevent <- diff(c(0,floor(round(cumsum(newcenevent),5)+0.5)))
    intsol <- c(nevent,intcenevent)
    int_obj <- -dvec%*%intsol + 0.5*t(intsol)%*%Dmat2%*%intsol
  }
  pen <- 0.5*cen_penalty*sum(intcenevent^2)
  data <- data.frame(time=S$time, nevent = nevent, ncen = intcenevent, tevent=apply(nevent,1,sum))
  data$nrisk <- sum(data$tevent+data$ncen) - cumsum(c(0,data$tevent+data$ncen))[1:dim(data)[1]]
  attr(data,"use_ticks") <- 1*(!is.null(ticks))
  attr(data,"S") <- origS
  attr(data,"c.event") <- !is.null(c.event)
  attr(data,"events") <- data.frame(t.event=t.event, c.event=c.event)
  attr(data,"risks") <- data.frame(t.risk=t.risk,n.risk=n.risk)
  attr(data,"Obs")<-Obs
  attr(data,"object") <- qpobj
  attr(data,"int_obj") <- int_obj
  attr(data,"int_obj_wpen") <- int_obj - pen
  attr(data,"mats") <- mats
  class(data) <- "CIFresolve"
  return(data)
}

print.CIFresolve <- function(x,...) {
  print(as.data.frame(lapply(x,as.vector)))
}

print.KMresolve <- function(x,...) {
  print(as.data.frame(lapply(x,as.vector)))
}

make_data <- function(cif, cen_method=NULL) {
  #cif: CIFresolve or KMresolve object
  #cen_method: Convention for censoring events if ticks not supplied;
  # "start": Censor at the start of the interval
  # "mid": Censor at the middle of the potential interval
  if (is.null(cen_method)) {
    if (attr(cif,"use_ticks")==1) {
      cen_method <-"start"
    }else{
      cen_method <- "mid"
    }
  }
  if (attr(cif,"use_ticks")==1 & cen_method=="mid") {
    warning("Mid-interval censoring only makes sense for pseudo-IPD without supplied tick marks.")
    cen_method <- "start"
  }
  if (!inherits(cif, "CIFresolve") & !inherits(cif,"KMresolve")) stop("Should supply an object made by CIF_resolve or KM_resolve")
  times<-event<-NULL
  if (inherits(cif,"CIFresolve") & length(cif)>5) {
    ncomp <- length(cif)-4
    for (i in 1:ncomp) {
      times <- c(times,rep(cif$time,unlist(cif[paste("nevent",i,sep=".")])))
      event <- c(event,rep(i,sum(unlist(cif[paste("nevent",i,sep=".")]))))
    }
  }else{
    times <- rep(cif$time,cif$nevent)
    event <- rep(1,sum(cif$nevent))
  }
  if (cen_method=="start") {
    times <- c(times,rep(cif$time,cif$ncen))
  }else{
    mids <- (cif$time + c(cif$time[2:length(cif$time)],cif$time[length(cif$time)]))/2
    times <- c(times,rep(mids,cif$ncen))
  }
  event <- c(event,rep(0,sum(cif$ncen)))
  data.frame(time=times,event=event)
}


plot.CIFresolve <- function(x,...) {
  cif <- x
  S <- attr(cif,"S")
  ncomp <- dim(S$cif)[2]
  opar <- par(no.readonly =TRUE)
  on.exit(par(opar))
  if (attr(cif,"c.event")) {
    par(mfrow=c(1,3))
  }else{
    par(mfrow=c(1,2))
  }
  plot(S$time,S$cif[,1],xlab="Time",ylab="CIF",lty=2,type="s",ylim=c(0,1))
  if (ncomp > 1) {
    for (j in 2:ncomp) {
      lines(S$time,S$cif[,j],type="s",col=j,lty=2)
    }
  }
  KM_l <- c(1,cumprod(1-cif$tevent/cif$nrisk))[1:length(cif$nrisk)]
  CIF <- list()
  if (ncomp>1) {
    for (j in 1:ncomp) {
      CIF[[j]] <- cumsum(cif[[paste("nevent",j,sep=".")]]/cif$nrisk * KM_l)
      lines(cif$time,CIF[[j]],type="s",col=j)
    }
  }else{
    CIF <- cumsum(cif$nevent/cif$nrisk * KM_l)
    lines(cif$time,CIF,type="s")
  }
  legend("topright",lty=c(1,2),legend=c("Reconstruction","Original"),bty="n")

  plot(cif$time,cif$nrisk,type="s",xlab="Time",ylab="Number at risk")
  points(attr(cif,"risks")$t.risk,attr(cif,"risks")$n.risk,pch=16)

  if (attr(cif,"c.event")) {
    #Additional panel for cumulative events.
    j <- dim(attr(cif,"events"))[2]
    if (j==2) {
      #Plot cumulative total events
      plot(cif$time, cumsum(cif$tevent),type="s",xlab="Time",ylab="Cumulative total events")
      points(attr(cif,"events")$t.event, attr(cif,"events")$c.event,pch=16)
    }else{
      plot(cif$time, cumsum(cif$nevent.1),type="s",xlab="Time",ylab="Cumulative events")
      points(attr(cif,"events")$t.event,attr(cif,"events")[,2],pch=16)
      for (j in 2:ncomp) {
        lines(cif$time, cumsum(cif[[paste("nevent",j,sep=".")]]),type="s",col=j)
        points(attr(cif,"events")$t.event,attr(cif,"events")[,(j+1)],col=j,pch=16)
      }
      #Plot cumulative events for each risk
    }
  }

}


plot.KMresolve <- function(x,...) {
  km <- x
  S <- attr(km,"S")
  opar <- par(no.readonly =TRUE)
  on.exit(par(opar))
  if (attr(km,"c.event")) {
     par(mfrow=c(1,3))
  }else{
    par(mfrow=c(1,2))
  }
  plot(S$time,S$surv,xlab="Time",ylab="S(t)",lty=2,type="s",ylim=c(0,1))
  KM <- cumprod(1-km$tevent/km$nrisk)
  lines(km$time,KM,type="s")
  legend("topright",lty=c(1,2),legend=c("Reconstruction","Original"),bty="n")

  plot(km$time,km$nrisk,type="s",xlab="Time",ylab="Number at risk")
  points(attr(km,"risks")$t.risk,attr(km,"risks")$n.risk,pch=16)

  if (attr(km,"c.event")) {
    plot(km$time,cumsum(km$nevent),type="s",xlab="Time",ylab="Cumulative events")
    points(attr(km,"events")$t.event,attr(km,"events")$c.event,pch=16)
  }
}

