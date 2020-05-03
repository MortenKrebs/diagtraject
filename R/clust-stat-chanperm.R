#' Clustering statistics with channel-wise permutations
#' 
#' Calculates clustering statistcs for the clustering of the (weighted) sum of a list of dissimilarity matrices. 
#' Performs (stratified) permutations and calculates clustering statistics for permuted  datasets.
#'   
#' @param dist_list list of dissimilarities in either format \code{'dist'} or as a symmetric \code{n x n} matrix.    
#' @param B number of permutations
#' @param maxclust maximum number of clusters
#' @param nclust number of clusters. Ignored if maxclust is specified.
#' @param strat.var optional factor of length n by which the permutation is stratified.  
#' @param stat vector of strings with cluster statistics to be calculated. Possible stats: 
#' \code{"ASW"} Average silhouette width.
#' \code{"CH"} Calinsky Harabasz index on squared distances.
#' @param c_weight Optional vector of numerics with channel weights.
#' @details ...
#' @return Object of class "ch_perm.stats".
#' @examples p <- channel_perm(d_imp_UD2,B = 50, maxclust = 15,strat.var = as.vector(TraMineR::seqlength(seq_mis[[1]])),c_weight = weigth, stat="ASW")
#' @export


channel_perm <-  function(dist_list, B, maxclust=NULL, nclust=NULL, strat.var=NULL, stat=c("ASW","CH"),c_weight=rep(1,length(dist_list))) {

if(!is.null(maxclust)){nclust <- 2:maxclust}  
l <- length(dist_list)
n <- nrow(dist_list[[1]])
if(is.null(strat.var)) strat.var <- rep(1,n)
dist_list <- lapply(dist_list, as.matrix)
dt <- data.table(id=1:n, strat=strat.var)
CH <- ASW <- perm_ASW <- perm_CH <- NULL
d_c <- Reduce('+',lapply(1:length(d), function(x) { d[[x]]*weigth[x] }))
c <- hclust(as.dist(d_c),method = "ward.D2")
if("CH" %in% stat){
  CH <- sapply(nclust, function(x) fpc::cluster.stats(d=d_c, cutree(c, x), sepindex=F, silhouette=F, sepwithnoise = F, wgap=F, aggregateonly = T)$ch)}
if("ASW" %in% stat){
  ASW <- sapply(nclust, function(x) summary(cluster::silhouette(dist=d_c, x= cutree(c, x)))$avg.width)
}
stats <- list(CH=CH, ASW=ASW)


  
perm_stats  <- lapply(1:B, function(x){
                      d <- list()
                      for(i in 1:l){
                        perm <- dt[,id[sample.int(.N,.N,F)], by= strat]$V1
                        d[[i]] <- dist_list[[i]][perm, perm]
                      }
                      d_c <- Reduce('+',lapply(1:length(d), function(x) { d[[x]]*weigth[x] }))
                      c <- hclust(as.dist(d_c),method = "ward.D2")
                      if("CH" %in% stat)
                      CH = sapply(nclust, function(x) fpc::cluster.stats(d=d_c, cutree(c, x), sepindex=F, silhouette=F, sepwithnoise = F, wgap=F, aggregateonly = T)$ch)
                      if("ASW" %in% stat)
                      ASW = sapply(nclust, function(x) summary(cluster::silhouette(dist=d_c, x= cutree(c, x)))$avg.width)
                      stats <- list(CH=CH, ASW=ASW)
                      return(stats)
                    })
if("ASW" %in% stat) {
  perm_ASW <- lapply(perm_stats, function(x) x$ASW) }
if("CH" %in% stat){
  perm_CH <- lapply(perm_stats, function(x) x$CH) }

out <- list(stats=stats, perm_stats= list(CH=perm_CH, ASW=perm_ASW))
class(out) <- "ch_perm.stats"
return(out)}

