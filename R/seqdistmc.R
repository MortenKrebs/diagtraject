#' Dissimilarity between multichannel sequence object - save combined sub.matrix
#' Modification \code{\link{TraMineR::seqdistmc}} to save combined sub.matrix
#' @param channels list of stslist objects
#' @param method
#' @param norm
#' @param indel
#' @param sm 
#' @param with.missing
#' @param full.matrix
#' @param link
#' @param cval
#' @param miss.cost
#' @param cweight
#' 

seqdistmc2 <- function (channels, method, norm = FALSE, indel = 1, sm = NULL, 
          with.missing = FALSE, full.matrix = TRUE, link = "sum", 
          cval = 2, miss.cost = 2, cweight = NULL) {
  nchannels <- length(channels)
  if (nchannels < 2) {
    stop("[!] please specify at least two channels")
  }
  if (is.null(cweight)) {
    cweight <- rep(1, nchannels)
  }
  numseq <- sapply(channels, nrow)
  if (any(numseq != numseq[1])) {
    stop(" [!] sequence objects have different numbers of rows")
  }
  numseq <- numseq[1]
  message(" [>] ", nchannels, " channels with ", numseq, " sequences")
  metlist <- c("OM", "LCS", "DHD", "HAM")
  if (!method %in% metlist) {
    stop(" [!] method must be one of: ", paste(metlist, 
                                               collapse = " "), call. = FALSE)
  }
  if (method == "LCS") {
    method <- "OM"
    sm <- "CONSTANT"
    indel <- 1
    cval <- 2
    miss.cost <- 2
  }
  timeVarying <- method %in% c("DHD")
  if (is.null(sm)) {
    costmethod <- "CONSTANT"
    if (method == "DHD") {
      costmethod <- "TRATE"
    }
    sm <- rep(costmethod, nchannels)
  }
  else if (length(sm) == 1 && sm %in% c("CONSTANT", "TRATE")) {
    sm <- rep(sm, nchannels)
  }
  if (length(indel) == 1) {
    indel <- rep(indel, nchannels)
  }
  if ((length(indel) != nchannels) || (length(sm) != nchannels) || 
        (length(cweight) != nchannels)) {
    stop(" [!] you should supply one weight, one substitution matrix and one indel per channel")
  }
  indel_list <- numeric(length = nchannels)
  substmat_list <- list()
  alphabet_list <- list()
  alphsize_list <- list()
  maxlength_list <- numeric(length = nchannels)
  for (i in 1:nchannels) {
    if (!inherits(channels[[i]], "stslist")) {
      stop(" [!] channel ", i, " is not a state sequence object, use 'seqdef' function to create one", 
           call. = FALSE)
    }
    alphabet_list[[i]] <- attr(channels[[i]], "alphabet")
    if (with.missing) {
      alphabet_list[[i]] <- c(alphabet_list[[i]], attr(channels[[i]], 
                                                       "nr"))
      message(" [>] including missing value as an additional state")
    }
    else {
      if (any(channels[[i]] == attr(channels[[i]], "nr"))) {
        stop(" [!] found missing values in channel ", 
             i, ", please set 'with.missing=T' to nevertheless compute distances")
      }
    }
    alphsize_list[[i]] <- length(alphabet_list[[i]])
    maxlength_list[i] <- ncol(channels[[i]])
    indel_list[i] <- indel[i]
    if (is.character(sm[[i]])) {
      message(" [>] computing substitution cost matrix for channel ", 
              i)
      substmat_list[[i]] <- seqsubm(channels[[i]], sm[[i]], 
                                    with.missing = with.missing, time.varying = timeVarying, 
                                    cval = cval, miss.cost = miss.cost)
    }
    else {
#       if (method == "OM") {
#         TraMineR.checkcost(sm[[i]], channels[[i]], with.missing = with.missing, 
#                            indel = indel[i])
#       }
#       else {
#         TraMineR.checkcost(sm[[i]], channels[[i]], with.missing = with.missing)
#       }
      substmat_list[[i]] <- sm[[i]]
    }
    substmat_list[[i]] <- cweight[i] * substmat_list[[i]]
  }
  slength1 <- seqlength(channels[[1]])
  for (i in 2:nchannels) {
    if (sum(slength1 != seqlength(channels[[i]])) > 0) {
      if (!with.missing) {
        stop(" [!] some channels have sequences of different length for the same individual. Please set 'with.missing=TRUE' to nevertheless compute distances")
      }
      else {
        warning(" [!] some channels have sequences of different length for the same individual. Shorter sequences will be filled with missing values.")
        break
      }
    }
  }
  message(" [>] building combined sequences...", appendLF = F)
  sep <- "@@@@TraMineRSep@@@@"
  maxlength = max(maxlength_list)
  newseqdata <- matrix("", nrow = numseq, ncol = maxlength)
  newseqdataNA <- matrix(TRUE, nrow = numseq, ncol = maxlength)
  for (i in 1:nchannels) {
    seqchan <- channels[[i]]
    void <- attr(seqchan, "void")
    nr <- attr(seqchan, "nr")
    for (j in 1:maxlength) {
      if (j > maxlength_list[i]) {
        newCol <- as.character(rep(void, numseq))
      }
      else {
        newCol <- as.character(seqchan[, j])
      }
      newseqdataNA[, j] <- newseqdataNA[, j] & newCol == 
        void
      newCol[newCol == void] <- nr
      if (i > 1) {
        newseqdata[, j] <- paste(newseqdata[, j], newCol, 
                                 sep = sep)
      }
      else {
        newseqdata[, j] <- newCol
      }
    }
  }
  newseqdata[newseqdataNA] <- NA
  alphabet_size <- length(unique(as.character(newseqdata))) - 
    as.integer(sum(is.na(newseqdata)) > 0)
  suppressMessages(newseqdata <- seqdef(newseqdata, cpal = rep("blue", 
                                                               alphabet_size)))
  message(" OK")
  message(" [>] computing combined substitution and indel costs...", 
          appendLF = FALSE)
  alphabet <- attr(newseqdata, "alphabet")
  alphabet_size <- length(alphabet)
  if (!timeVarying) {
    newsm <- matrix(0, nrow = alphabet_size, ncol = alphabet_size)
    for (i in 1:(alphabet_size - 1)) {
      statelisti <- strsplit(alphabet[i], sep)[[1]]
      for (j in (i + 1):alphabet_size) {
        cost <- 0
        statelistj <- strsplit(alphabet[j], sep)[[1]]
        for (chan in 1:nchannels) {
          ipos <- match(statelisti[chan], alphabet_list[[chan]])
          jpos <- match(statelistj[chan], alphabet_list[[chan]])
          cost <- cost + substmat_list[[chan]][ipos, 
                                               jpos]
        }
        newsm[i, j] <- cost
        newsm[j, i] <- cost
      }
    }
  }
  else {
    newsm <- array(0, dim = c(alphabet_size, alphabet_size, 
                              maxlength))
    for (t in 1:maxlength) {
      for (i in 1:(alphabet_size - 1)) {
        statelisti <- strsplit(alphabet[i], sep)[[1]]
        for (j in (i + 1):alphabet_size) {
          cost <- 0
          statelistj <- strsplit(alphabet[j], sep)[[1]]
          for (chan in 1:nchannels) {
            ipos <- match(statelisti[chan], alphabet_list[[chan]])
            jpos <- match(statelistj[chan], alphabet_list[[chan]])
            cost <- cost + substmat_list[[chan]][ipos, 
                                                 jpos, t]
          }
          newsm[i, j, t] <- cost
          newsm[j, i, t] <- cost
        }
      }
    }
  }
  message(" OK")
rownames(newsm) <- colnames(newsm) <- gsub(alphabet, pattern = "@@@@TraMineRSep@@@@",replacement = "")

  newindel <- sum(indel_list * cweight)
  if (link == "mean") {
    newindel <- newindel/sum(cweight)
    newsm <- newsm/sum(cweight)
  }
  message(" [>] computing distances ...")
  return(list(dist=seqdist(newseqdata, method = method, norm = norm, 
                 indel = newindel, sm = newsm, with.missing = FALSE, 
                 full.matrix = full.matrix),sm=newsm))
}
