#' Imputation of missing states 
#' 
#' imputes states given state at last observation before censoring and timespecific transition rates and calculate probability weighted substitution costs. 
#' @param seq_mis object of class \code{'stslist'} as created by \code{\link[TraMineR:seqdef]{seqdef}}.   
#' @param cens.type character indicating the type of censoring. Must be either \code{"rigth"}, \code{"left"} or \code{"both"}.     
#' @param last optional caracter string containing state levels at last observation before censoring. 
#' @param trans.rates object of class \code{'array'} with transition rates in format created by  \code{\link[TraMineR:seqtrate]{seqtrate}}. If \code{NULL} (default) transition rates will be calculated using  \code{\link[TraMineR:seqtrate]{seqtrate}} function.
#' @param smooth \code{'logical'} indicating if transition rates should be smoothed. 
#' @param sm  \code{'character'} indicating substitution cost setting. Must be \code{"CONSTANT"} or \code{"TRATE"} for sm calculated by \code{\link[TraMineR::seqsubm]{TraMineR::seqsubm}} or object of class \code{'matrix'} containing substitution costs. 
#' @param method Currently only \code{"prob"} (default) is available. See \code{Details}.
#' @param prob.out logical indicating if imputed probabilities should be included in output. Defaults to \code{FALSE}.
#' @param diag logical indicating if diagonal should be printed in dist object. Defaults to \code{FALSE}.
#' @param resol.ratio optional numeric specified if increaments differ between calculations of dissimmilarities in imputed and complete sequences. Defaults to \code{1}.
#' @param resol.comp optional vector of integers. If increaments differ between calculations of dissimmilarities in complete and imputed sequences the differences can be specified for compensation. 
#' @param mc.cores optional integer specifying the number of cores for parallel computation.
#' @details Calculates dissimilarities for right and left censored state sequence objects using probability weighted substitution costs
#'  \deqn{ d_{inf}(i,j) = \sum\limits_{t=1}^{t_{max}} \sum Pr(i)_t Pr(j)_t^T \circ  SC }
#' @return Object of class \code{'dist'} containing dissimilarities.
#' @examples 
#' 
#' ## Creating a sequence object
#' data(mvad)
#' mvad.alphabet <- c("employment", "FE", "HE", "joblessness", "school", 
#'                    "training")
#' mvad.labels <- c("employment", "further education", "higher education", 
#'                  "joblessness", "school", "training")
#' mvad.scodes <- c("EM", "FE", "HE", "JL", "SC", "TR")
#' mvad.seq <- seqdef(mvad, 17:86, alphabet = mvad.alphabet, states = mvad.scodes, 
#'                    labels = mvad.labels, xtstep = 6)
#' 
#' ## Introducing right-censoring
#' addMissing <- function(x){
#' if(is.factor(x)) return(factor(x, levels=c(levels(x), "missing")))
#' return(x)}
#' mvad.perm <- mvad
#' mvad.perm <- as.data.frame(lapply(mvad.perm, addMissing))
#' row.perm.r <- sample(1:nrow(mvad))[1:floor(nrow(mvad)*.5)]
#' row.perm <- 1:nrow(mvad) %in% row.perm.r 
#' col.perm.r <- sample(floor(ncol(mvad[,17:86])*.8):ncol(mvad[,17:86]),size = length(row.perm.r),replace = T)     
#' for(i in 1:length(row.perm.r)){
#'   mvad.perm[row.perm.r[i],(col.perm.r[i]+16):ncol(mvad)] <- "missing"}
#' perm.seq <- seqdef(mvad.perm, 17:86, alphabet = mvad.alphabet, states = mvad.scodes, missing = "missing", labels = mvad.labels, xtstep = 6)
#' 
#' ## Computing Hamming distance in observed states
#' perm.seq2 <- seqdef(mvad.perm, 17:86, xtstep = 6)
#' sub.cost2 <- seqsubm(seqdata = perm.seq2, method = "CONSTANT")
#' sub.cost2["missing->",] <- sub.cost2[,"missing->"] <- 0
#' dist.obs <- seqdist(perm.seq2, method = "HAM", sm = sub.cost2)
#' 
#' ## Computing Probability weighted Hamming distance in censored states:
#' dist.mis <- mis.cost(perm.seq, cens.type="right",sum_to_1 = F, 
#' method = "prob",sm = "CONSTANT",smooth = F)
#' dist <- dist.obs + as.matrix(dist.mis$dist)
#' 
#' ## Obtaining imputed probabilities
#' prob <- mis.cost(perm.seq, cens.type="right",sum_to_1 = F, 
#' method = "prob",smooth = F, prob.out=T)
#' 
#' @export







mis.cost <- function(seq_mis, cens.type = c("right","left","both"),
                     last=NULL, trans.rates=NULL, smooth=F,
                     sum_to_1=T,MM=T, sm="CONSTANT", method="prob", prob.out=F, diag = F, resol.comp=NULL, resol.ratio=1, mc.cores=NULL){
  debut <- proc.time()
  if(prob.out==T) {message(" [>] NB: Returning imputed probabilities only.")
                   if(cens.type =="both") stop("prob.out=T can only be used for right of left censored, not both - do them seperately!")}
  if(sum(sm %in% c("CONSTANT","TRATE"))) { 
    message(" [>] creating substitution cost matrix using ", sm, " costs.")
    sm <- seqsubm(seq_mis, method = sm)} else message(" [>] using user defined substitution cost matrix.") 
  if(sum(grepl(rownames(sm), pattern = "miss"))){ 
    message(" [>] removing missing state(s) \"",grep(rownames(sm), pattern = "miss", value=T)  ,"\" from substitution cost-matrix.")
    sm <- sm[-grep(rownames(subst.cost), pattern = "miss"),-grep(colnames(sm), pattern = "miss")]}
  tri <- tri.ineq2(sm)
  if(length(tri$ineq)>1) message("Triangular inequality violated: \n substituting \"", 
                                     rownames(sm)[tri$arr.indices[2,1]], "\" with \"",
                                     rownames(sm)[tri$arr.indices[2,2]], "\" costs ", sm[tri$arr.indices[2,1],tri$arr.indices[2,2]], "\n substituting \"", 
                                     rownames(sm)[tri$arr.indices[2,1]], "\" with \"",
                                     rownames(sm)[tri$arr.indices[2,3]], "\" costs ", sm[tri$arr.indices[2,1],tri$arr.indices[2,3]], "\n substituting \"", 
                                     rownames(sm)[tri$arr.indices[2,3]], "\" with \"",
                                     rownames(sm)[tri$arr.indices[2,2]], "\" costs ", sm[tri$arr.indices[2,3],tri$arr.indices[2,2]]) else(message(" [>] Substitution matrix respects triangular inequality"))

  
   max_imp <-function(i) {
    message("     [>] comparing nr ", i, " shortest sequence to")
    li <- seqlength(seqd_o[i,])
    seqd_lend <- seqd_o[i:nrow(seqd_o),c(li,li:gmax)]
    seqd2_lend <- seqd2_o[i:nrow(seqd_o),c(li,li:(gmax+1))]
    seqd2_lend[,c(1,2)] <- "None"
    tr <- trans.rates[,,c(li,li:gmax)]
    #tr <- trans.rates[,,c(li,li:gmax,gmax)]
    seqd2_l <- unique(seqd2_lend) 
    if(method=="bg.samp") seqd2_l <- seqd2_lend 
    seqd_l <- seqd_lend[rownames(seqd_lend) %in% rownames(seqd2_l),]
    #print(seqd_l)    
    message("     [>] " ,nrow(seqd_l), " distinct equal or longer sequences")
    mcorr_l <-match(seqconc(seqd2_lend),seqconc(seqd2_l))
    mis_l <- seqlength(seqd_l) #length of j 
    lmax <- max(mis_l)
    seqd_l <- seqd_l[,c(1:max(mis_l),ncol(seqd_l))] #extending sequences
    stmin <- seqd2_o$last[i]
    #t_min <- c((1:lmax)[-c(1:2)],lmax) # seq i imputed from li+1 to lmax 
    t_min <- c((1:lmax)[-c(1:2)]) # positions in tr to be imputed for seq i
    if(prob.out==T){
      # Using markov algebra:
      if(length(t_min)==0) p_mint <-0 else {
      p_mint <- sapply(1:length(t_min), function(z) as.numeric(paste0("[",stmin," ->]")== rownames(trans.rates)) %*% Reduce("%*%",lapply(1:z, function(z2) tr[,,t_min[z2]])))
      rownames(p_mint) <- colnames(tr) }
    p_mint
    
    }else{
      
    d_l <- lapply(1:nrow(seqd_l), function(j){
      if(li==gmax) {0}  else{
      stmax <- seqd2_l$last[j]
#       t_max <- c((mis_l[j]:lmax)[-c(1)],lmax)#,lmax)  # seq j imputed lj+1 to lmax
#       if(mis_l[j]==lmax) t_max <- c(t_max,lmax)
      t_max <- c((mis_l[j]:lmax)[-c(1)])#,lmax)  #positions in tr to be imputed for seq j
 if(mis_l[j]==lmax) t_max2 <- t_min else if(min(t_max)==3) t_max2 <- NULL else t_max2 <- seq(2,min(t_max)-1)[-c(1)]
if (!is.null( t_max2)) stjt <- unlist(seqd_l[j, t_max2]) else stjt <- NULL# vector of observed states

  #t_max2 <- t_min[seq(2,length(t_min)-length(t_max),l=length(t_min)-length(t_max))]
# 
# cat("i impute:") 
# print(t_min)
# cat("j impute:") 
# print(t_max)
# cat("j obs:") 
# print(t_max2)
# print(stjt)
      #    (1:(length(t_min)-length(t_max)))]  # seq j still observ        
     

      # create alphabet x t probability-matrix for shorter sequence  
      if(MM==F){
        p_mint <- t(apply(tr[paste0("[",stmin," ->]"), ,t_min], 1, function(z1) sapply(
          1:length(t_min), function(z2) 1-prod(1-z1[1:z2]))))
        if(min(p_mint<0)) stop("negative probabilities")
        p_mint[paste0("[-> ",stmin,"]"),] <- 0
        p_mint[paste0("[-> ",stmin,"]"),] <- 1-colSums(p_mint)}else{
          # Using markov algebra:
          p_mint <- sapply(1:length(t_min), function(z) as.numeric(paste0("[",stmin," ->]")== rownames(trans.rates)) %*% Reduce("%*%",lapply(1:z, function(z2) tr[,,t_min[z2]])))
          rownames(p_mint) <- colnames(tr)}
      
      # create alphabet x t probability-matrix for longer sequence  
      p_max1 <- sapply(stjt, function(z) as.numeric(paste0("[",z," ->]")== rownames(trans.rates))) 
      if(MM==F){  
        p_max2 <- t(apply(tr[paste0("[",stmax," ->]"), ,t_max], 1, function(z1) sapply(
          1:length(t_max), function(z2) 1-prod(1-z1[1:z2]))))
        p_max2[paste0("[-> ",stmax,"]"),] <- 0
        p_max2[paste0("[-> ",stmax,"]"),] <- 1-colSums(p_max2)}else{
          # Using markov algebra:
          if(mis_l[j]==lmax) p_max2 <- NULL else {p_max2 <- sapply(1:length(t_max), function(z) as.numeric(paste0("[",stmax," ->]")== rownames(trans.rates)) %*% Reduce("%*%",lapply(1:z, function(z2) tr[,,t_max[z2]])))
          rownames(p_max2) <- colnames(tr)}}
      p_max <- cbind(p_max1, p_max2)
    
      if(method=="bg" | method=="bg.samp"){ 
        p_mint_tmp <- p_mint
        p_max_tmp <- p_max
        p_mint[1:nrow(p_mint),] <- p_max[1:nrow(p_max),]<- 0
        for(x1 in 1:ncol(p_mint)){
          if(method=="bg"){
        p_mint[which.max(p_mint_tmp[,x1]),x1] <-1
        p_max[which.max(p_max_tmp[,x1]),x1] <-1} else {
        p_mint[sample(1:nrow(p_mint),1,prob=p_mint_tmp[,x1]),x1] <-1
        p_max[sample(1:nrow(p_max),1,prob=p_max_tmp[,x1]),x1] <-1}
        }}
#       if(method=="bg.samp"){ 
#         #cat(apply(p_array, c(1,3),FUN=x))
#         p_mint <- apply(p_mint, 2, FUN=function(x) as.numeric(sample(x,1,prob = x+1e-10)==x))
#         p_max <- apply(p_max, 2, FUN=function(x) as.numeric(sample(x,1,prob = x+1e-10)==x))
#       } 
      
      
      p_array <- array(sapply(1:ncol(p_mint), function(z1) 
        sapply(p_mint[,z1], function(z2) z2*unlist(p_max[,z1]))), 
        dim = c(nrow(p_mint),nrow(p_max),ncol(p_mint)))

      #if(sum(p_array[,,dim(p_array)[3]])==0) {cat(stmax)}
      dt  <- apply(p_array, 3, function(z) sum(sm*z)) 
      #dt  <- apply(p_array, 3, function(z) sum(z)) 

      #if(max(dt[-c((ncol(p_mint)-1):ncol(p_mint))])>2.1) print(p_mint)
      #      cat(length(dt[-c((ncol(p_mint)-1):ncol(p_mint))]))
      if(!is.null(resol.comp)) { dt[1] <- dt[1]*(1-resol.comp2[i]*resol.ratio)}   

sum(dt)     #  <----  #[-c((ncol(p_mint)-1):ncol(p_mint))])
      #length(dt)
      #dt[-c((ncol(p_mint)-1):ncol(p_mint))] -resol.comp2[i]*resol.ratio
    }}
    )
    d_o <- d_l[mcorr_l]
    d_o}
    #lapply(d_o, function(x) length(x)- resol.comp2[i]*resol.ratio)
    #[order(as.numeric(rownames(seqd_lend)))]
  }
reps <- 1
if(cens.type=="both") {reps <- 2}  
  
for(k in 1:reps) {

  if(cens.type=="right"|cens.type=="both") { message(" [>] imputing right-censored sequences")   
  l <- seqlength(seq_mis)
  message(" [>] ", nrow(seq_mis), " sequences.")
  message(" [>] sequence length from ", min(l), " to ", max(l), "." ) 
  
  if(is.null(trans.rates)){
    message(" [>] computing transition rates" )
    
    trans.rates <- TraMineR::seqtrate(seq_mis, time.varying = T)
  }}

  if(cens.type=="left") { 
    message(" [>] imputing left-censored sequences")
    #seq_mis <- seqdef(seq_mis[,ncol(seq_mis):diff(range(seqlength(seq_mis)))],missing="*",void = "#")
    seq_mis_r <- seq_mis[,ncol(seq_mis):1]
    trans.rates <- seqtrate(seq_mis_r, time.varying = T)
    seq_mis_r[seq_mis_r=="%"]<- alphabet(seq_mis)[1]
    seq_mis <- seqdef(seq_mis_r, missing="*")
    #seq_mis <-seq_mis[,ncol(seq_mis):(diff(range(seqlength(seq_mis)))+1)]
    l <- seqlength(seq_mis)
    message(" [>] ", nrow(seq_mis), " sequences.")
    message(" [>] sequence length from ", min(l), " to ", max(l), "." )}
  
  if(smooth){
    message(" [>] smoothing transition rates" )
    trans.rates <- aperm(apply(trans.rates, 1:2, smooth),c(2,3,1))}
  if(sum_to_1){
    message(" [>] making all transition probabilities sum to 1 by assuming same state if they don't")
    tr2 <- trans.rates
    not1 <- which(!apply(tr2[,,],c(1,3),sum)==1,arr.ind = T)
    not1s <- apply(tr2[,,],c(1,3),sum)
    for(f in 1:length(not1[,2]))
    diag(tr2[,,not1[,2][f]])[not1[,1][f]] <-1 -not1s[not1[,1][f],not1[,2][f]]
    trans.rates <- tr2
  }
  
  seqend <- seq_mis[,min(l):max(l)] #cut sequence object to area with missing states
  trans.rates <- trans.rates[,,(min(l)-1):(max(l)-1)] #cut trans.rates object to area with missing states
  
  seqend2 <- seqend
  
  if(!is.null(last)){ if(!length(last)==nrow(seq_mis)) stop("length of last must be the the same as as number of observations in sequence object") 
  #last_s <- last[1:nrow(seqend)] #temp
  seqend2$last <- last_s } else {seqend2$last <- sapply(1:nrow(seqend2), function(x) seqend2[x,][seqlength(seqend2[x,])])}
  if(!is.null(resol.comp)){
  resol.comp <- resol.comp[1:nrow(seqend)]  #temp
  seqend2$resol <- resol.comp }
  seqd2 <- unique(seqend2)
  if(method=="bg.samp") seqd2 <- seqend2
  seqd <-seqend[rownames(seqend2) %in% rownames(seqd2),]
  #resol.comp <-resol.comp[rownames(seqend2) %in% rownames(seqd2)]
  mcorr <-match(seqconc(seqend2),seqconc(seqd2)) 
  if(!is.null(resol.comp)){ resol.comp2 <-seqd2$resol}
  #seqd2 <- seqd2[,-ncol(seqend2)]
  len <- nrow(seqd)
  rownames(seqd) <- 1:len
  rownames(seqd2) <- 1:len
  message(" [>] ", len, " distinct sequence endings after removing first ", min(l)-1, " states."  )
  mis<-seqlength(seqd)
  if(!is.null(resol.comp)){
    seqd_o <- seqd[order(mis,resol.comp2),]
    seqd2_o <- seqd2[order(mis,resol.comp2),]} else {
      seqd_o <- seqd[order(mis),]
      seqd2_o <- seqd2[order(mis),]
    }
  if(!is.null(resol.comp)){resol.comp2 <- seqd2_o$resol}
  gmax <-max(mis)
  if(!sum(c(dim(sm),nrow(trans.rates),ncol(trans.rates))==length(c(alphabet(seqd),  unique(last)[which(!unique(last) %in% alphabet(seqd))])))==4) stop("substition cost matrix must have dimentions alph.length X alpha.length")
  pi <- rep(list(list()),len)
  

  
  if(is.null(mc.cores)){
    message(" [>] computing substitution cost using imputed states uptill length of longest sequence")
    pi <-  lapply(1:len, max_imp)} else{
    library(parallel)
    message(" [>] computing (in parallel) substitution cost using imputed states uptill length of longest sequence")
    pi <-  mclapply(1:len, max_imp, mc.cores=mc.cores)}
  if(prob.out==T){
    
    out <- list(prob.out = pi[order(as.numeric(rownames(seqd_o)))][mcorr])}else{
  mf <- matrix(NA,length(pi), length(pi))
  mf[lower.tri(mf,diag = T)] <- unlist(pi)
  mf[upper.tri(mf)] <-  t(mf)[upper.tri(mf)]
  m <- mf[order(as.numeric(rownames(seqd_o))),order(as.numeric(rownames(seqd_o)))][mcorr,mcorr]
  rownames(m) <-colnames(m) <- rownames(seq_mis)
  if(!diag) m <-  as.dist(m) 
  if(cens.type=="right"|cens.type=="both") {  
  right <- list(dist = m/resol.ratio, 
             order= order(as.numeric(rownames(seqd_o))),
             pi=pi, mcorr= mcorr, mf = mf, last = seqend2$last) }
  if(cens.type=="left") {  
    left <- list(dist = m/resol.ratio, 
                  order= order(as.numeric(rownames(seqd_o))),
                  pi=pi, mcorr= mcorr, mf = mf, last = seqend2$last) }
  if(cens.type=="both") { cens.type <- "left" }

  
}}

fin <- proc.time()
totaltime <- format(round(difftime(as.POSIXct(sum(fin[1:2]), 
                                              origin = "1960-01-01"), as.POSIXct(sum(debut[1:2]), 
                                                                                 origin = "1960-01-01")), 3))
message(" [>] total time: ", totaltime)

if(reps==2) return(list(dist =right$dist + left$dist)) else {
  if(prob.out==T){ return(out) } else {
    if(cens.type=="right") return(right)
    if(cens.type=="left") return(left) }}
}
