# diagtraject

Here we provide the R-code for imputation of right-censored sequence data used in the paper "Patterns in Comorbid Diagnostic Trajectories of Individuals with Schizophrenia Associate with Etiological Factors". 

The r-package can be installed by running the following lines of code in R: 

``` 
library(devtools)
install_github("MortenKrebs/diagtraject")
```

We here provide an example using publicly available life course data:

```
## Creating a sequence object
data(mvad)
mvad.alphabet <- c("employment", "FE", "HE", "joblessness", "school", 
                   "training")
mvad.labels <- c("employment", "further education", "higher education", 
                 "joblessness", "school", "training")
mvad.scodes <- c("EM", "FE", "HE", "JL", "SC", "TR")
mvad.seq <- seqdef(mvad, 17:86, alphabet = mvad.alphabet, states = mvad.scodes, 
                   labels = mvad.labels, xtstep = 6)

## Introducing right-censoring
addMissing <- function(x){
if(is.factor(x)) return(factor(x, levels=c(levels(x), "missing")))
return(x)}
mvad.perm <- mvad
mvad.perm <- as.data.frame(lapply(mvad.perm, addMissing))
row.perm.r <- sample(1:nrow(mvad))[1:floor(nrow(mvad)*.5)]
row.perm <- 1:nrow(mvad) %in% row.perm.r 
col.perm.r <- sample(floor(ncol(mvad[,17:86])*.8):ncol(mvad[,17:86]),size = length(row.perm.r),replace = T)     
for(i in 1:length(row.perm.r)){
  mvad.perm[row.perm.r[i],(col.perm.r[i]+16):ncol(mvad)] <- "missing"}
perm.seq <- seqdef(mvad.perm, 17:86, alphabet = mvad.alphabet, states = mvad.scodes, missing = "missing", labels = mvad.labels, xtstep = 6)

## Computing Hamming distance in observed states
perm.seq2 <- seqdef(mvad.perm, 17:86, xtstep = 6)
sub.cost2 <- seqsubm(seqdata = perm.seq2, method = "CONSTANT")
sub.cost2["missing->",] <- sub.cost2[,"missing->"] <- 0
dist.obs <- seqdist(perm.seq2, method = "HAM", sm = sub.cost2)
 
## Computing Probability weighted Hamming distance in censored states:
dist.mis <- mis.cost(perm.seq, cens.type="right",sum_to_1 = F, 
method = "prob",sm = "CONSTANT",smooth = F)
dist <- dist.obs + as.matrix(dist.mis$dist)
```
