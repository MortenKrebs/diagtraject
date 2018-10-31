#' Modification of \code{\link{wcCmpCluster\{WeightedCluster\}}} that distinguishes hclust(... ,method="ward.D") and hclust(... ,method="ward.D2").
#' 
#' See \code{\link{wcCmpCluster\{WeightedCluster\}}}

wcCmpCluster2 <- function (diss, weights = NULL, maxcluster, method = "all", 
                           pam.combine = TRUE) {
  if (maxcluster < 2) {
    stop(" [!] maxcluster should be greater than 2")
  }
  hclustmethods <- c("ward.D", "ward.D2", "single", "complete", "average", 
                     "mcquitty", "median", "centroid")
  all.methods <- c(hclustmethods, "pam", "diana", "beta.flexible")
  noweights.methods <- c("diana", "beta.flexible")
  if (any(method == "all")) {
    if (is.null(weights)) {
      method <- all.methods
    }
    else {
      method <- c(hclustmethods, "pam")
    }
  }
  if (!is.null(weights) && any(method %in% noweights.methods)) {
    stop(" [!] methods ", paste(noweights.methods, collapse = ", "), 
         " cannot be used with weights.")
  }
  ret <- list()
  if (any(method %in% hclustmethods)) {
    dd <- as.dist(diss)
  }
  pamrange <- 2:maxcluster
  range <- lapply(1:length(method), function(meth){
    #range <- foreach(meth= 1:length(method), .packages = c("WeightedCluster","TraMineR")) %dopar% {
    if (method[meth] != "pam") {
      if (method[meth] %in% hclustmethods) {
        hc <- hclust(dd, method = method[meth], members = weights)
      } else if (method[meth] == "diana") {
        hc <- diana(diss, diss = TRUE)
      } else if (method[meth] == "beta.flexible") {
        hc <- agnes(diss, diss = TRUE, method = "flexible", 
                    par.method = 0.625)
      }
      ret <- as.clustrange(hc, diss = diss, weights = weights, 
                           ncluster = maxcluster)
      l <- list(ret)
      names(l) <- c(method[meth])
      if (pam.combine) {
        pam_result <- wcKMedRange(diss, kvals = pamrange, weights = weights, initialclust = hc)
        l <- list(ret, pam_result)
        names(l) <- c(method[meth], paste0(method[meth],".pam"))}
    }
    else {
      ret <- wcKMedRange(diss, kvals = pamrange, 
                         weights = weights)
      l <- list(ret)
      names(l) <- c(method[meth])}
    l})
  range <- unlist(range, recursive = F)
  allstats <- list()
  for (meth in names(range)) {
    allstats[[meth]] <- cbind(range[[meth]]$stats, method = meth, 
                              ngroup = pamrange)
  }
  range$param <- list(method = method, pam.combine = pam.combine, 
                      all.methods = names(range), kvals = pamrange)
  range$allstats <- do.call(rbind, allstats)
  class(range) <- "clustrangefamily"
  return(range)
}