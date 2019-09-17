bootmds <- function(dist, k,  w, nrep = 100, mc=NULL, verbose = FALSE, ...) 
{
library(vegan)  
library(smacof)
library(parallel)
    
  n <- nrow(dist)          ## number of objects
  p <- k        ## number of dimensions
 # val <- object$stress     
  #smacall <- object$call
  
  #N <- dim(data)[1]
  coord <- list()
  object <- wcmdscale(d = dist, k = k,w = w)
  
coord <- lapply(1:nrep,#mc.cores = mc, 
                function(x) {
   samp <- sample(1:n, size = n, replace = TRUE,prob = w)
   

   key <- cbind(as.numeric(names(table(samp))),table(samp))
   st <- dist[key[,1],key[,1]]     
   ## bootstrap sample data
    message("Replication: ", x, "\n")
    
    ## compute input dissimilarities
    
   if(!is.null(w)){ 
     mds <-wcmdscale(d = st, k = k,w = key[,2])
   }else{ mds <-wcmdscale(d = st, k = k)}
   
   message(dim(object), " ", dim(mds[match(samp,key[,1]),]))
   
   fit <- Procrustes(object, mds[match(samp,key[,1]),])
   fit$Yhat
  }  )
  
  ## stability measure
print(lapply(coord,dim))
  y0 <- Reduce("+", coord)/length(coord)
  stab.num <- sum(sapply(coord, function(yy) (norm(yy-y0))^2))
  stab.denom <- sum(sapply(coord, function(yy) (norm(yy))^2))
  stab <- 1 - stab.num/stab.denom
  stab}
