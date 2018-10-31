#' Modification of \code{\link{prediction.strength\{fpc\}}} that allows input format to be a dissimilarity matrix and clustering method to be hierachical.
#'
#' @param dist distance matrix or object of class \code{'dist'}.
#' @return object of class predstr

pred.strength.dist <- function(dist,Gmin= 2, Gmax= 10, M=50, cutoff = .8, clustermethod=c("ward.D2","average"), class.method=c("ward.D2","average","pam"),
                               indx= NULL, weights=NULL, mc.cores=NULL){

xdata <- as.matrix(dist)
#xdata[upper.tri(xdata)] <- NA 
n <- nrow(xdata)
if(!is.null(indx)) n <-length(indx) else indx <-1:n 
nf <- c(floor(n/2), n - floor(n/2))
indvec <-mcor <-members<- samp <- clcenters <- clusterings <- jclusterings <- classifications <- list()
corrpred <- list()
pred <- list()
k=Gmin:Gmax
if(class.method=="ward.D2")    {clust_full <- hclust(as.dist(xdata), method = "ward.D2",members=weights)}
corrpred_func <- function(k){
  
  
  for (l in 1:M) {
    pred[[l]] <- numeric(0)
    nperm <- sample(n, n)
    indvec[[l]] <-samp[[l]]<-mcor[[l]]<-members[[l]]<- list()
    samp[[l]][[1]] <- indx[nperm[1:nf[1]]]
    samp[[l]][[2]] <-indx[nperm[(nf[1] + 1):n]]
    indvec[[l]][[1]] <- unique(samp[[l]][[1]])
    indvec[[l]][[2]] <- unique(samp[[l]][[2]])
    mcor[[l]][[1]] <- match(samp[[l]][[1]],indvec[[l]][[1]])
    mcor[[l]][[2]] <- match(samp[[l]][[2]],indvec[[l]][[2]])
    if(!is.null( weights)) { 
      wt<- data.frame("w"=weights[samp[[l]][[1]]], "corr" = mcor[[l]][[1]]) 
      members[[l]][[1]] <-  aggregate(w~corr,data = wt,sum)
      wt<- data.frame("w"=weights[samp[[l]][[2]]], "corr" = mcor[[l]][[2]]) 
      members[[l]][[2]] <-  aggregate(w~corr,data = wt,sum)
    }
    for (i in 1:2) {
     
      clusterings[[i]] <- disthclustCBI(as.dist(xdata[indvec[[l]][[i]], 
                                                      indvec[[l]][[i]]]), 
                                        method = clustermethod, members=members[[l]][[i]], k)
  if(class.method=="ward.D2")    {
#     w <- rep(1, n)
#       w[indvec[[l]][[i]]] <- 1e-10
#       classifications[[i]] <- disthclustCBI(as.dist(xdata), members=w, method = "ward.D2", k)$partition[indvec[[l]][[i]]]
    classifications[[i]] <- cutree(clust_full,k)[samp[[l]][[i]]]}
#   if(class.method=="average")    {
#     w <- rep(1, n)
#     w[indvec[[l]][[i]]] <- 0
#     classifications[[i]] <- disthclustCBI(as.dist(xdata), members=w, method = "average", k)$partition[indvec[[l]][[i]]]}
#   
  if(class.method=="pam"){
    jclusterings[[i]] <- rep(-1, n)
    jclusterings[[i]][indvec[[l]][[i]]] <- clusterings[[i]]$partition
    medoids <- sapply(1:k, function(x) {
     m <- as.data.frame(as.matrix(as.dist((xdata[indvec[[l]][[i]], indvec[[l]][[i]]])[clusterings[[i]]$partition==x,clusterings[[i]]$partition==x] )))
    med <- which(names(which.min(colMeans(m)))== rownames(xdata)) 
    if(length(med)==0) med <-  indvec[[l]][[i]][clusterings[[i]]$partition==x]
     med[sample(1:length(med),1)]
     })
  
    j <- 3 - i
    classifications[[j]] <- classifdist(as.dist(xdata), 
                                        jclusterings[[i]], method = "centroid", 
                                        centroids = medoids, nnk = nnk)[indvec[[l]][[j]]]
    }
    }
    ps <- matrix(0, nrow = 2, ncol = k)
    for (i in 1:2) {
      ctable <- table(clusterings[[i]]$partition[mcor[[l]][[i]]], 
                      classifications[[i]][mcor[[l]][[1]]])
      for (kk in 1:k) {
        ps[i, kk] <- sum(ctable[kk, ]^2 - ctable[kk, 
                                                 ])
        cpik <- clusterings[[i]]$partition[mcor[[l]][[i]]] == kk
        print(length(clusterings[[i]]$partition))
        print(length(unique(mcor[[l]][[i]])))
        
        nik <- sum(cpik)
        if (nik > 1) 
          ps[i, kk] <- ps[i, kk]/(nik * (nik - 1))
        else ps[i, kk] <- 1
      }
    }
    pred[[l]] <- mean(c(min(ps[1, ]), min(ps[2, ])))
  }
  unlist(pred)}

if(!is.null(mc.cores)){ 
  corrpred <- mclapply(k, corrpred_func, mc.cores=mc.cores) } else {
  corrpred <- lapply(k, corrpred_func)}

mean.pred <- numeric(0)
mean.pred <- c(1)
cp <- list()
cp[Gmin:Gmax] <- corrpred 
for (k in Gmin:Gmax) 
  mean.pred <- c(mean.pred, mean(cp[[k]]))
optimalk <- max(which(mean.pred > cutoff))
out <- list(predcorr = cp, mean.pred = mean.pred, 
            optimalk = optimalk, cutoff = cutoff, method = clustermethod, 
            Gmax = Gmax, M = M)
class(out) <- "predstr"
out}
