#' Clusterwise cluster stability - Weighted 
#' Modification \code{\link{fpc::clusterboot}} to weight channels in each permutation, to sva
#' @param data n x p dataframe or dissimilarity matrix (if distances = TRUE) or list of dissimilarity matrices (if multichannel = TRUE)  
#' @param indvec optional vector of indices if data is a dissimilarity matrix or  \code{'dist'} object repressenting only unique observations.
#' @param dt.m optional list of dataframe with data needed to weight channels (if multichannel = TRUE)  
#' 
#' 


disthclustCBI <-  function (dmatrix, k, cut = "number", method, weights=NULL, noisecut = 0, 
          ...) 
{
  n <- nrow(as.matrix(dmatrix))
  c1 <- fastcluster::hclust(as.dist(dmatrix), method = method, members=weights)
  noise <- FALSE
  if (cut == "number") 
    partition <- cutree(c1, k = k)
  else partition <- cutree(c1, h = k)
  nc <- max(partition)
  clsizes <- numeric(0)
  for (i in 1:nc) clsizes[i] <- sum(partition == i)
  ncn <- sum(clsizes > noisecut)
  if (ncn < nc) {
    noise <- TRUE
    newcln <- (1:nc)[clsizes > noisecut]
    nc <- ncn + 1
    newpart <- rep(nc, n)
    for (i in 1:ncn) newpart[partition == newcln[i]] <- i
    partition <- newpart
  }
  cl <- list()
  for (i in 1:nc) cl[[i]] <- partition == i
  out <- list(result = c1, noise = noise, nc = nc, nccl = ncn, 
              clusterlist = cl, partition = partition, clustermethod = "hclust")
  out
}



mdshclustCBI <-  function (dmatrix, k, k.mds, cut = "number", method, weights=NULL, noisecut = 0, 
                           ...) 
{
  n <- nrow(as.matrix(dmatrix))
  mds <- cmdscale(dmatrix,k = k.mds)
  c1 <- fastcluster::hclust(dist(mds), method = method, members=weights)
  noise <- FALSE
  if (cut == "number") 
    partition <- cutree(c1, k = k)
  else partition <- cutree(c1, h = k)
  nc <- max(partition)
  clsizes <- numeric(0)
  for (i in 1:nc) clsizes[i] <- sum(partition == i)
  ncn <- sum(clsizes > noisecut)
  if (ncn < nc) {
    noise <- TRUE
    newcln <- (1:nc)[clsizes > noisecut]
    nc <- ncn + 1
    newpart <- rep(nc, n)
    for (i in 1:ncn) newpart[partition == newcln[i]] <- i
    partition <- newpart
  }
  cl <- list()
  for (i in 1:nc) cl[[i]] <- partition == i
  out <- list(result = c1, noise = noise, nc = nc, nccl = ncn, 
              clusterlist = cl, partition = partition, clustermethod = "hclust")
  out
}



clusterbootw <- function (data, indvec=NULL, dt.m=NULL, multichannel=F, B = 100, distances = (class(data) == "dist"), 
                          bootmethod = "boot", bscompare = TRUE, multipleboot = FALSE, 
                          jittertuning = 0.05, noisetuning = c(0.05, 4), subtuning = floor(nrow(as.matrix(data))/2), 
                          clustermethod, members= NULL, noisemethod = FALSE, count = TRUE, showplots = FALSE, 
                          dissolution = 0.5, recover = 0.75, mc.cores=1, seed = NULL, ...) 
{
  sumlogic <- function(x, y, relation = "eq") switch(relation, 
                                                     eq = sum(x == y, na.rm = TRUE), s = sum(x < y, na.rm = TRUE), 
                                                     l = sum(x > y, na.rm = TRUE), se = sum(x <= y, na.rm = TRUE), 
                                                     le = sum(x >= y, na.rm = TRUE))
  if (!is.null(seed)) 
    set.seed(seed)
  invisible(distances)
  if(multichannel) data <- lapply(data, as.matrix) else data <- as.matrix(data)
  if (distances & showplots) 
    dpoints <- cmdscale(data)
  if(multichannel) {n <- nrow(data[[1]]) 
                    p <- ncol(data[[1]])
                    }else {n <- nrow(data)
                           p <- ncol(data[[1]])}
  if(!is.null(indvec)) {n <- length(indvec)
                        subtuning <- floor(n/2)}
  p <- ncol(data[[1]])
  #cod <- cov(data)
  #md <- colMeans(data)
  lb <- length(bootmethod)
  if(multichannel){
    dt1 <- dt.m
    we<- sapply(1:length(dt1), function(x) 1/(ceiling(max(dt1[[x]][,time]))-mean(dt1[[x]][event!="missing", min(time), by= pid]$V1)))
    cdata <- Reduce('+', lapply(1:length(data), function(x) {data[[x]]*we[x] }))
    c1 <- clustermethod(as.dist(cdata), ...) }else 
    {cdata<- as.dist(data)
     c1 <- clustermethod(as.dist(cdata), weights=members, ...)}
  if (noisemethod) {
    if (c1$nccl == 0) 
      stop("No clusters, only noise estimated!")
  }
  else c1$nccl <- c1$nc
  bootresult <- jitterresult <- noiseresult <- bojitresult <- subsetresult <- matrix(0, 
                                                                                     nrow = c1$nc, ncol = B)
  lbootresult <- rep(0,c1$nc)
  if(is.null(indvec))  {
    lbootpartition <- rep(0,length(c1$partition)) 
    lbootpartition_uni <-rep(0,length(c1$partition)) }else {
      lbootpartition <- rep(0,length(indvec))
      lbootpartition_uni <-rep(0,length(c1$partition)) }
  bootpartition <-subsetpartition <- matrix(0, 
                          nrow = length(lbootpartition), ncol = B)
  bootpartition_uni <-subsetpartition_uni <- matrix(0, 
                                            nrow = length(lbootpartition_uni), ncol = B)
  
  if (("jitter" %in% bootmethod) | ("bojit" %in% bootmethod)) {
    jsd <- numeric(0)
    ecd <- eigen(cod, symmetric = TRUE)
    ecd$values[ecd$values < 0] <- 0
    ecd$values[is.na(ecd$values)] <- 0
    rotdata <- data %*% solve(t(ecd$vectors))
    for (i in 1:p) {
      sx <- sort(rotdata[, i])
      dx <- sx[2:n] - sx[1:(n - 1)]
      dx <- dx[dx > 0]
      jsd[i] <- quantile(dx, jittertuning)
    }
  }
  if ("noise" %in% bootmethod) {
    ecd <- eigen(cod, symmetric = TRUE)
    ecd$values[ecd$values < 0] <- 0
  }
  if (showplots) {
    if (distances) 
      plot(dpoints, pch = sapply(c1$partition, toString), 
           col = c1$partition)
    else plot(data, pch = sapply(c1$partition, toString), 
              col = c1$partition)
  }
  for (l in 1:lb) {
    bootlist <-mclapply(1:B, mc.cores= mc.cores, 
                       FUN= function(i) {
      if (count) 
        cat(bootmethod[l], i, "\n")
      if (bootmethod[l] == "boot") {
        bsamp <- bsamp_ind <-  sample(n, n, replace = TRUE)
        if(!is.null(indvec)) bsamp <- indvec[bsamp]
        if (!multipleboot) 
          bsamp <- unique(bsamp)
        if (distances) {
          if(multichannel){
          dt1 <- lapply(dt.m, function(x) x[pid %in% unique(x[,pid])[bsamp]])
          weigth <- sapply(1:length(dt1), function(x) 1/(ceiling(max(dt1[[x]][,time]))-mean(dt1[[x]][event!="missing", min(time), by= pid]$V1)))
          mdata <- Reduce('+', lapply(1:length(data), function(x) { data[[x]][bsamp,bsamp]*weigth[x] }))
        }else {mdata <- as.dist(as.matrix(cdata)[bsamp,bsamp]) }}}
      if (bootmethod[l] == "subset") {
        samp <- sample(n, subtuning, replace = FALSE)
        print(head(samp))
        if(!is.null(indvec)) {bsamp_ind <-  indvec[samp]                      
        bsamp_unique <- unique(bsamp_ind)
        bsamp_corr <- match(bsamp_ind, bsamp_unique)
        bsamp <- bsamp_unique
        print(head(bsamp))
        if(!is.null(members)) { 
          wt<- data.frame("w"=members[bsamp_ind], "corr" = bsamp_corr)
    
          b_members <- aggregate(w~corr,data = wt,sum)
          print(head(b_members))}
        }
        if (distances) {
          if(multichannel){
            dt1 <- lapply(dt.m, function(x) x[pid %in% unique(x[,pid])[bsamp]])
            weigth <- sapply(1:length(dt1), function(x) 1/(ceiling(max(dt1[[x]][,time]))-mean(dt1[[x]][event!="missing", min(time), by= pid]$V1)))
            mdata <- Reduce('+', lapply(1:length(data), function(x) { data[[x]][bsamp,bsamp]*weigth[x] }))
          }else {mdata <- as.dist(as.matrix(cdata)[bsamp,bsamp]) }}
      }
      print(dim(as.matrix(mdata)))
      if (bootmethod[l] == "jitter") {
        jnoise <- matrix(0, ncol = p, nrow = n)
        for (j in 1:p) jnoise[, j] <- rnorm(n, sd = jsd[j])
        jnoise <- jnoise %*% t(ecd$vectors)
        mdata <- data + jnoise
        bsamp <- 1:n
      }
      if (bootmethod[l] == "bojit") {
        bsamp <- sample(n, n, replace = TRUE)
        jnoise <- matrix(0, ncol = p, nrow = n)
        for (j in 1:p) jnoise[, j] <- rnorm(n, sd = jsd[j])
        jnoise <- jnoise %*% t(ecd$vectors)
        mdata <- data[bsamp, ] + jnoise
      }
      if (bootmethod[l] == "noise") {
        noiseind <- as.logical(rbinom(n, 1, noisetuning[1]))
        nn <- sum(noiseind)
        jnoise <- matrix(0, ncol = p, nrow = nn)
        for (j in 1:p) jnoise[, j] <- runif(nn, min = -noisetuning[2] * 
                                              sqrt(ecd$values[j]), max = noisetuning[2] * 
                                              sqrt(ecd$values[j]))
        jnoise <- t(t(jnoise %*% t(ecd$vectors)) + md)
        mdata <- data
        mdata[noiseind, ] <- jnoise
        bsamp <- (1:n)[!noiseind]
      }
      if(!is.null(members)){
      bc1 <- clustermethod(as.dist(mdata), weights=b_members$w, ...)}else{    bc1 <- clustermethod(as.dist(mdata), ...)}
      print(head(bc1$partition))
      if (showplots) {
        if (distances) 
          plot(dpoints[bsamp, ], pch = sapply(bc1$partition, 
                                              toString), col = bc1$partition)
        else plot(mdata, pch = sapply(bc1$partition, 
                                      toString), col = bc1$partition)
      }
      if (noisemethod) {
        effnc1 <- c1$nccl
        effnb1 <- bc1$nccl
      }
      else {
        effnc1 <- c1$nc
        effnb1 <- bc1$nc
      }
      for (j in 1:effnc1) {
        maxgamma <- 0
        if (effnb1 > 0) {
          for (k in 1:effnb1) {
            if (multipleboot) {
              if (bscompare) 
                ncases <- 1:n
              else {
                ncases <- 1
                m <- 2
                if (m <= n) {
                  if (!(bsamp[m] %in% bsamp[1:(m - 1)])) 
                    ncases <- c(ncases, m)
                  m <- m + 1
                }
              }
            }
            else ncases <- 1:length(bsamp)
            cg <- switch(bootmethod[l], boot = clujaccard(c1$clusterlist[[j]][bsamp][ncases], 
                                                          bc1$clusterlist[[k]][ncases], zerobyzero = 0), 
                         bojit = clujaccard(c1$clusterlist[[j]][bsamp][ncases], 
                                            bc1$clusterlist[[k]][ncases], zerobyzero = 0), 
                         subset = clujaccard(c1$clusterlist[[j]][bsamp_ind], 
                                             bc1$clusterlist[[k]][bsamp_corr], zerobyzero = 0), 
                         jitter = clujaccard(c1$clusterlist[[j]], 
                                             bc1$clusterlist[[k]], zerobyzero = 0), 
                         noise = clujaccard(c1$clusterlist[[j]][!noiseind], 
                                            bc1$clusterlist[[k]][!noiseind], zerobyzero = 0))
            if (cg > maxgamma) 
              maxgamma <- cg
          }
        }
        print(head(maxgamma))
        if (bootmethod[l] == "boot") 
        {cat("c1: ")
         cat(paste(length(c1$partition)))
         cat("\n")
         cat("bc1: ")
         cat(paste(length(bc1$partition)))
         cat("\n")
         cat("bsamp: ")
         cat(paste(length(bsamp_ind)))
         cat("\n")
         cat("ncases: ")
         cat(paste(length(ncases)))
         lbootresult[j] <- maxgamma
         if(!is.null(indvec)) {
         lbootpartition[bsamp_ind] <- bc1$partition[bsamp_corr]
         lbootpartition_uni[bsamp_unique] <- bc1$partition}else{
           lbootpartition[bsamp] <- bc1$partition
         }}
        if (bootmethod[l] == "subset") {
          lbootresult[j] <- maxgamma
          lbootpartition[bsamp_ind] <- bc1$partition[bsamp_corr]
          lbootpartition_uni[bsamp_unique] <- bc1$partition}
        if (bootmethod[l] == "bojit") 
          bojitresult[j, i] <- maxgamma
        if (bootmethod[l] == "jitter") 
          jitterresult[j, i] <- maxgamma
        if (bootmethod[l] == "noise") 
          noiseresult[j, i] <- maxgamma
      }
      if (noisemethod) {
        if (c1$nc > c1$nccl) {
          j <- c1$nc
          if (bc1$nc > bc1$nccl) 
            maxgamma <- switch(bootmethod[l], boot = clujaccard(c1$clusterlist[[c1$nc]][bsamp][ncases], 
                                                                bc1$clusterlist[[bc1$nc]][ncases], zerobyzero = 0), 
                               bojit = clujaccard(c1$clusterlist[[c1$nc]][bsamp][ncases], 
                                                  bc1$clusterlist[[bc1$nc]][ncases], zerobyzero = 0), 
                               subset = clujaccard(c1$clusterlist[[c1$nc]][bsamp], 
                                                   bc1$clusterlist[[bc1$nc]], zerobyzero = 0), 
                               jitter = clujaccard(c1$clusterlist[[c1$nc]], 
                                                   bc1$clusterlist[[bc1$nc]], zerobyzero = 0), 
                               noise = clujaccard(c1$clusterlist[[c1$nc]][!noiseind], 
                                                  bc1$clusterlist[[bc1$nc]][!noiseind], 
                                                  zerobyzero = 0))
          else maxgamma <- 0
          if (bootmethod[l] == "boot") 
            lbootresult[j] <- maxgamma
          if (bootmethod[l] == "subset") 
            subsetresult[j, i] <- maxgamma
          if (bootmethod[l] == "bojit") 
            bojitresult[j, i] <- maxgamma
          if (bootmethod[l] == "jitter") 
            jitterresult[j, i] <- maxgamma
          if (bootmethod[l] == "noise") 
            noiseresult[j, i] <- maxgamma
        }
      }
   bootout <- list(result = lbootresult, partition = lbootpartition, partition_uni = lbootpartition_uni) 
   bootout})
}
  if (!("boot" %in% bootmethod)) 
    bootresult <- bootmean <- bootbrd <- bootrecover <- NULL
  else {
    for(i in 1:B){
    bootresult[,i] <- bootlist[[i]]$result
    bootpartition[,i] <- bootlist[[i]]$partition
    }
    bootmean = apply(bootresult, 1, mean, na.rm = TRUE)
    bootbrd = apply(bootresult, 1, sumlogic, y = dissolution, 
                    relation = "se")
    bootrecover = apply(bootresult, 1, sumlogic, y = recover, 
                        relation = "l")
  }
  if (!("jitter" %in% bootmethod)) 
    jitterresult <- jittermean <- jitterbrd <- jitterrecover <- NULL
  else {
    jittermean = apply(jitterresult, 1, mean, na.rm = TRUE)
    jitterbrd = apply(jitterresult, 1, sumlogic, y = dissolution, 
                      relation = "se")
    jitterrecover = apply(jitterresult, 1, sumlogic, y = recover, 
                          relation = "l")
  }
  if (!("subset" %in% bootmethod)) 
    subsetresult <- subsetmean <- subsetbrd <- subsetrecover <- NULL
  else {
    for(i in 1:B){
      subsetresult[,i] <- bootlist[[i]]$result
      subsetpartition[,i] <- bootlist[[i]]$partition
      subsetpartition_uni[,i] <- bootlist[[i]]$partition_uni
    }
    subsetmean = apply(subsetresult, 1, mean, na.rm = TRUE)
    subsetbrd = apply(subsetresult, 1, sumlogic, y = dissolution, 
                    relation = "se")
    subsetrecover = apply(subsetresult, 1, sumlogic, y = recover, 
                        relation = "l")
  }
  if (!("noise" %in% bootmethod)) 
    noiseresult <- noisemean <- noisebrd <- noiserecover <- NULL
  else {
    noisemean = apply(noiseresult, 1, mean, na.rm = TRUE)
    noisebrd = apply(noiseresult, 1, sumlogic, y = dissolution, 
                     relation = "se")
    noiserecover = apply(noiseresult, 1, sumlogic, y = recover, 
                         relation = "l")
  }
  if (!("bojit" %in% bootmethod)) 
    bojitresult <- bojitmean <- bojitbrd <- bojitrecover <- NULL
  else {
    bojitmean = apply(bojitresult, 1, mean, na.rm = TRUE)
    bojitbrd = apply(bojitresult, 1, sumlogic, y = dissolution, 
                     relation = "se")
    bojitrecover = apply(bojitresult, 1, sumlogic, y = recover, 
                         relation = "l")
  }
  if (showplots) {
    if (distances) 
      plot(dpoints, pch = sapply(c1$partition, toString), 
           col = c1$partition)
    else plot(data, pch = sapply(c1$partition, toString), 
              col = c1$partition)
  }
  out <- list(result = c1, partition = c1$partition, nc = c1$nc, 
              nccl = c1$nccl, clustermethod = c1$clustermethod, B = B, 
              noisemethod = noisemethod, bootpartition = bootpartition, bootmethod = bootmethod, 
              multipleboot = multipleboot, dissolution = dissolution, 
              recover = recover, bootresult = bootresult, bootmean = bootmean, 
              bootbrd = bootbrd, bootrecover = bootrecover, jitterresult = jitterresult, 
              jittermean = jittermean, jitterbrd = jitterbrd, jitterrecover = jitterrecover, 
              subsetresult = subsetresult, subsetmean = subsetmean, 
              subsetbrd = subsetbrd, subsetrecover = subsetrecover, subsetpartition=subsetpartition,
              subsetpartition_uni=subsetpartition_uni,
              bojitresult = bojitresult, bojitmean = bojitmean, bojitbrd = bojitbrd, 
              bojitrecover = bojitrecover, noiseresult = noiseresult, 
              noisemean = noisemean, noisebrd = noisebrd, noiserecover = noiserecover)
  class(out) <- "clboot"
  out
}
