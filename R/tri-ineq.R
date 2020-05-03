#' Modification of \code{\link{tri.ineq\{fossil\}}} that gives magnitude and indices.
#' 
#'   
#' @param d distance matrix
#' @return list containing \code{ineq} and \code{arr.indices}.
#' @export

#tri.ineq{fossil} Rip:

tri.ineq2 <- function(d){
mat <- d
n <- dim(mat)[1]
ineq <- 0
arr <- c(0,0,0)
for (i in 1:(n - 2)) {
  for (j in (i + 1):(n - 1)) {
    for (k in (j + 1):n) {
      sds <- c(mat[j, i], mat[k, i], mat[k, j])
      lng <- max(sds)
      if (lng > (sum(sds) - lng)) {
        ineq <- c(ineq, lng - (sum(sds) - lng))        
        arr <- rbind(arr,c(i,j,k))
    }}
  }
}
if(length(ineq)>1) message("Triangular inequality violated")
list(ineq=ineq,arr.indices=arr)}
