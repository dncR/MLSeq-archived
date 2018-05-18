####### 1. NBLDA FUNCTIONS #######
GetDnb <- function(ns, x, y,beta){
  uniq <- sort(unique(y))
  ds <- matrix(1, nrow = length(uniq), ncol = ncol(x))

  for (k in 1:length(uniq)){
    a <- colSums(x[y == uniq[k], ]) + beta
    b <- colSums(ns[y == uniq[k], ]) + beta
    ds[k, ] <- a/b
  }

  return(ds)
}


predictNBLDA <- function(x, y, xte = NULL, beta = 1, type = c("mle", "deseq", "quantile", "none", "TMM"),
                         prior = NULL, truephi = NULL, null.out = NULL, ds = NULL, disperhat = NULL, ...){
  ## Args:
  ##    x: raw count matrix which is used while training NBLDA classifier. Samples in the rows and genes in the columns
  ##    y: a vector of classes for train set samples.
  ##    xte: a raw count matrix for test samples.
  ##    beta: a numeric value. This is the beta parameter from trained model.
  ##    type: normalization type. Should be same as in trained model.
  ##    prior: prior probabilities of each class. Should be same as in trained model.
  ##    truephi:
  ##    null.out: train set parameters.
  ##    ds: offset parameters for each gene. gene expression parameters.
  ##    disperhat: dispersion estimates from trained model.

  if (is.null(prior)){
    prior <- rep(1/length(unique(y)), length(unique(y)))
    ## prior vectorunun sinif sayisi ile esit uzunlukta olmasi lazım.
    ## Bu durumu kontrol edip gerekirse uyari verilecek.
  }

  # null.out <- NullModel(x, type = type)  ## Trained model'den alınacak.
  # ns <- null.out$n
  nste <- NullModelTest(null.out, x, xte, type = type)$nste

  uniq <- sort(unique(y))
  # ds <- GetDnb(ns, x, y, beta)  ### Trained modelden alınacak.

  ###### ÖNEMLİ:
  ###### GENLERE AİT DISPERSION KESTİRİMLERİ TRAIN EDILEN MODELDEN ALINIYOR.
  ###### BU DOĞRU BİR YAKLAŞIM MIDIR İNCELENECEK!!!!
  # if (!is.null(truephi)){
  #   if (length(truephi) >= 2){
  #     truephi <- truephi[1]
  #     warning("\"truephi\" should be given as a single value. Only first element is used.")
  #   }
  #   disperhat <- rep(truephi, ncol(nste))
  #
  # } else {
  #   tt = getT(t(x), sizeFactors = rep(1, nrow(x)), verbose = FALSE)$target
  #
  #   ### Moment estimation of gene-wise dispersions.
  #   rM = rowMeans(t(x))
  #   rV = rowVars(t(x))
  #
  #   disp = (rV - rM)/rM^2
  #   disp0 = numeric()
  #
  #   for (i in 1:length(disp)){
  #     a1 = 0
  #     a2 = disp[i]
  #     disp0[i] = max(a1, a2)  ## Negative dispersions are set to 0.
  #   }
  #   names(disp0) <- names(disp)
  #
  #   disperhat = getAdjustDisp(disp0, shrinkTarget = tt, verbose = FALSE)$adj
  # }

  phihat <- as.numeric(disperhat)
  discriminant <- matrix(NA, nrow = nrow(xte), ncol = length(uniq))

  # Replace Infinity with zero dispersions.
  inv.phihat <- (1/phihat)
  if (any(inv.phihat == Inf | inv.phihat == -Inf)){
    id.inf <- which(inv.phihat == Inf | inv.phihat == -Inf)
    inv.phihat[id.inf] <- 0
  } else {
    id.inf <- NULL
  }

  for (k in 1:length(uniq)){
    for (l in 1:nrow(xte)){
      dstar = ds[k, ]
      part2 = 1 + nste[l, ] * dstar * phihat
      part1 = dstar/part2

      part3 <- (1/phihat) * log(part2)

      if (!is.null(id.inf)){
        part3.limits <- nste[l, ] * dstar
        part3[id.inf] <- part3.limits[id.inf]
      }

      discriminant[l, k]<- sum(xte[l, ] * log(part1)) - sum(part3) + log(prior[k])
    }
  }

  preds <- uniq[apply(discriminant, 1, which.max)]  ## Predicted class labels

  return(preds)
}

######## 1. PLDA FUNCTIONS ##########

# x[tr, ], y[tr], x[te, ], rhos = rhos, beta = beta,
# type = "quantile", prior = prior, transform = FALSE
## type seçenekleri "TMM", "deseq" ve "none" olacak.
Classify <- function(x, y, xte = NULL, rho = 0, beta = 1, rhos = NULL, type = c("mle", "deseq", "quantile", "TMM"),
                     prior = NULL, transform = TRUE, alpha = NULL, return.selected.gene.names = FALSE){
  if (is.null(xte)){
    xte <- x
    warning("Since no xte was provided, testing was performed on training data set.")
  }
  if (!is.null(rho) && length(rho) > 1) stop("Can only enter 1 value of rho. If you would like to enter multiple values, use rhos argument.")
  type <- match.arg(type)
  if (!transform && !is.null(alpha)) stop("You have asked for NO transformation but have entered alpha.")

  if (transform && is.null(alpha)){
    alpha <- FindBestTransform(x)  ### Train ve test seti buradaki "alpha" değerine göre transform edilmeli.
  }

  if (transform){
    if (alpha <= 0 || alpha > 1) stop("alpha must be between 0 and 1")
    x <- x^alpha
    xte <- xte^alpha
  }

  if (is.null(prior)) prior <- rep(1/length(unique(y)), length(unique(y)))
  if (is.null(rho) && is.null(rhos)) stop("Must enter rho or rhos.")

  null.out <- NullModel(x, type = type)
  ns <- null.out$n
  nste <- NullModelTest(null.out, x, xte, type = type)$nste

  uniq <- sort(unique(y))
  if (is.null(rhos)){
    ds <- GetD(ns, x, y, rho, beta)
    discriminant <- matrix(NA, nrow = nrow(xte), ncol = length(uniq))
    for (k in 1:length(uniq)){
      discriminant[ ,k] <- rowSums(scale(xte, center = FALSE, scale = (1/log(ds[k, ])))) - rowSums(scale(nste, center = FALSE, scale = (1/ds[k, ]))) + log(prior[k])
    }
    save <- list(ns = ns, nste = nste, ds = ds, discriminant = discriminant, ytehat = uniq[apply(discriminant, 1, which.max)],
                 alpha = alpha, rho = rho, x = x, y = y, xte = xte, alpha = alpha, type = type, prior = prior,
                 transform = transform, trainParameters = null.out)
    # class(save) <- "poicla"
    return(save)
  } else {
    save <- list()
    ds.list <- GetD(ns, x, y, rho = NULL, rhos = rhos, beta)
    for (rho in rhos){
      ds <- ds.list[[which(rhos == rho)]]
      discriminant <- matrix(NA, nrow = nrow(xte), ncol = length(uniq))
      for (k in 1:length(uniq)){
        discriminant[,k] <- rowSums(scale(xte, center = FALSE, scale = (1/log(ds[k, ])))) - rowSums(scale(nste, center = FALSE, scale = (1/ds[k, ]))) + log(prior[k])
      }
      save[[which(rhos == rho)]] <- (list(ns = ns, nste = nste, ds = ds, discriminant = discriminant,
                                          ytehat = uniq[apply(discriminant, 1, which.max)], alpha = alpha, prior = prior,
                                          rho = rho, x = x, y = y, xte = xte, alpha = alpha, type = type, transform = transform,
                                          trainParameters = null.out))
    }
    # class(save) <- "poicla"
    return(save)
  }
}


predictPLDA <- function(x, y = NULL, xte = NULL, rho = 0, beta = 1, type = c("mle", "deseq", "quantile", "TMM"),
                        prior = NULL, transform = TRUE, alpha = NULL, null.out = NULL, ds = NULL){
  ## Args:
  ##  x: raw count matrix which is used while training PLDA classifier. Samples in the rows and genes in the columns
  ##  y: a vector of classes for train set samples.
  ##  xte: a raw count matrix for test samples.
  ##  rho: a numeric value. This is the best value of tuning parameter.
  ##  beta: a numeric value. This is the beta parameter from trained model.
  ##  type: normalization type. Should be same as in trained model.
  ##  prior: prior probabilities of each class. Should be same as in trained model.
  ##  transform: TRUE/FALSE.
  ##  alpha: a numeri value between 0 and 1. It is used to perform power transformation on raw data.
  ##  null.out: train set parameters.
  ##  ds: offset parameters for each gene. gene expression parameters.

  type <- match.arg(type)

  if (transform){
    if (alpha <= 0 || alpha > 1){
      stop("alpha must be between 0 and 1")
    }
    x <- x^alpha
    xte <- xte^alpha
  }

  if (is.null(prior)){ ### prior trained model'den alınacak.
    prior <- rep(1/length(unique(y)), length(unique(y)))
  }

  # null.out <- NullModel(x, type = type)  ### trained modelden alınacak
  # ns <- null.out$n
  nste <- NullModelTest(null.out, x, xte, type = type)$nste

  uniq <- sort(unique(y))
  # ds <- GetD(ns, x, y, rho, beta)   ### trained modelden alınacak
  discriminant <- matrix(NA, nrow = nrow(xte), ncol = length(uniq))
  for (k in 1:length(uniq)){
    discriminant[ ,k] <- rowSums(scale(xte, center = FALSE, scale = (1/log(ds[k, ])))) - rowSums(scale(nste, center = FALSE, scale = (1/ds[k, ]))) + log(prior[k])
  }
  pred <- uniq[apply(discriminant, 1, which.max)]
  return(pred)
}

## Generate Count Data using Poisson Model.
CountDataSet <- function(n, p, K, param, sdsignal){
  if (n < 4*K) stop("We require n to be at least 4*K.")
  q0 <- rexp(p, rate = 1/25)
  isDE <- runif(p) < 0.3
  classk <- matrix(NA, nrow=K, ncol=p)
  for (k in 1:K){
    lfc <- rnorm(p, sd = sdsignal)
    classk[k,] <- ifelse(isDE, q0*exp(lfc), q0)
  }
  truesf <- runif(n)*2+.2 #size factors for training observations
  truesfte <- runif(n)*2+.2 #size factors for test observations
  conds <- sample(c(rep(1:K, 4), sample(1:K, n-4*K, replace=TRUE))) # class labels for training observations
  condste <- sample(c(rep(1:K, 4), sample(1:K, n-4*K, replace=TRUE))) # class labels for test observations
  x <- xte <- matrix(NA, nrow=n, ncol=p)
  for (i in 1:n){
    for (k in 1:K){
      if (conds[i]==k) x[i,] <- rnbinom(p, mu = truesf[i]*classk[k,], size=param)
      if (condste[i]==k) xte[i,] <- rnbinom(p, mu = truesfte[i]*classk[k,], size=param)
    }
  }
  rm <- apply(x,2,sum)==0
  return(list(x=x[,!rm],xte=xte[,!rm],y=conds,yte=condste, truesf=truesf, truesfte=truesfte))
}



## Find best value of trainsformation parameter "alpha" for power transformation.
FindBestTransform <- function(x){
  alphas <- seq(.01, 1, len = 50)
  gof <- rep(NA, length(alphas))
  for (alpha in alphas){
    gof[alphas == alpha] <- GoodnessOfFit(x^alpha, type = "mle")
  }
  return(alphas[which.min(abs(gof - (nrow(x) - 1)*(ncol(x) - 1)))])
}


NullModel <- function(x, type = c("mle", "deseq", "quantile", "none", "TMM")){
  # x MUST be a n times p matrix - i.e. observations on the rows and features on the columns
  type <- match.arg(type)
  rowsum <- rowSums(x)
  colsum <- colSums(x)

  if (type == "mle"){
    sizes <- rowSums(x)/sum(x)  # size factors for each sample
    mle <- outer(sizes, colsum, "*")
    return(list(n = mle, sizes = sizes))

  } else if (type == "quantile"){
    # This is quantile normalization idea of Bullard et al 2010 -- quantile-normalize using 75th quantile of observed counts for each sample, excluding zero-counts
    sample.qts <- apply(x, 1, quantile, .75)
    sample.qts <- pmax(sample.qts, 1) # Don't wait to standardize by 0... min allowed is 1
    sample.qts <- sample.qts/sum(sample.qts)   # size factors for each sample
    fit <- outer(sample.qts, colsum, "*")
    return(list(n = fit, sizes = sample.qts))

  } else if (type == "deseq"){ #Trying to implement idea from Anders and Huber Genome Biology 2010 paper.
    # I stole the next 3 lines from the DESeq bioconductor package.... it was in the function estimateSizeFactorsForMatrix
    counts <- t(x)
    geomeans <- exp(rowMeans(log(counts)))   # Geometric means for each features.
    sizes <- apply(counts, 2, function(cnts) median((cnts/geomeans)[geomeans > 0]))
    rawsizestr <- sizes  ## This part will be used to normalize test samples using train sample information.
    sizes <- sizes/sum(sizes)  # size factors for each sample
    fit <- outer(sizes, colsum, "*")
    return(list(n = fit, sizes = sizes, geomeans = geomeans, rawsizestr = rawsizestr))

  } else if (type == "TMM"){
    countsDGE <- DGEList(counts = t(x), genes = colnames(x))
    countsDGE.normalized <- calcNormFactors(countsDGE, method = "TMM")   ## RLE: DESeq mantigi ile normalize ediyor.

    # This function is used to find reference sample.
    # Codes are copied from edgeR and modified here.
    # rawCounts are count matrix with genes in the row and samples in the column
    findRefSample <- function (rawCounts, lib.size, p = 0.75){
      y <- t(t(rawCounts) / lib.size)
      f75 <- apply(y, 2, function(x){
        quantile(x, p = p)
      })
      refSample <- which.min(abs(f75 - mean(f75)))
      return(refSample)
    }

    nf <- countsDGE.normalized$samples$norm.factors   ## Normalization factors
    lsize <- countsDGE.normalized$samples$lib.size  ## Library sizes

    rawCounts <- t(x)
    refSampleID <- findRefSample(rawCounts, lib.size = colSums(rawCounts), p = 0.75)

    # TMM Normalized Counts
    # From: http://grokbase.com/t/r/bioconductor/127r11kp23/bioc-edger-and-libsize-normalization
    # n.normalized <- 1e6 * (x / (nf * lsize))   ## same results with cpm(..., log = FALSE)
    n.normalized <- t(cpm(countsDGE.normalized, log = FALSE))  # genes in the column, samples in the rows.

    return(list(n = n.normalized, refSampleID = refSampleID, refSample = x[refSampleID, ],
                normalization.factors = nf, lib.size = lsize))

  } else if (type == "none"){


  }
}


# type olarak "none" eklenebilir mi? bu durum incelenecek.
NullModelTest <- function(null.out, x, xte = NULL, type = c("mle", "deseq", "quantile", "none", "TMM")){
  type <- match.arg(type)

  if (is.null(xte)){
    xte <- x
  }

  if (type == "mle"){
    sizeste <- rowSums(xte)/sum(x)
    nste <- outer(sizeste, colSums(x), "*")

  } else if (type == "quantile"){
    sizeste <- pmax(1, apply(xte, 1, quantile, .75)) # don't want to scale by anything less than 1...
    sizeste <- sizeste/sum(apply(x, 1, quantile, .75))
    nste <- outer(sizeste, colSums(x), "*")

  } else if (type == "deseq"){
    countste <- t(xte)
    ##### sizeste İLE HESAPLANAN SONUÇLAR MLSEQ'DE YAZDIĞIMIZ sizeF.ts İLE AYNI SONUCU vermiyor.
    ##### BURADA hesaplanan sizefactor değerleri kendi içinde tekrar orantılanarak kullanılıyor.
    ##### normalized counts değerleri ise her genetotal değerinin sizefactorler ile çarpımı olarak elde ediliyor.
    sizeste <- apply(countste, 2, function(cnts) median((cnts/null.out$geomeans)[null.out$geomeans > 0], na.rm=TRUE))/sum(null.out$rawsizestr)
    nste <- outer(sizeste, colSums(x), "*")

  } else if (type == "TMM"){
    referenceSample <- null.out$refSample

    ## Check if the feature names of reference sample are in the same order with rawCounts.
    if (identical(colnames(xte), names(referenceSample))){
      xte <- rbind(xte, referenceSample)
    } else {
      stop(warning("Reference sample either does not have same features or the features are not in the same order as test set. Calculation stops.."))
    }

    countsDGE <- DGEList(counts = t(xte), genes = colnames(xte))

    ## Reference column is selected from train set.
    ## Train Set'de kullanılan ref sample'a göre normalization factor hesaplanıyor.
    countsDGE.nf <- calcNormFactors(countsDGE, method = "TMM", refColumn = ncol(countsDGE))

    nf <- countsDGE.nf$samples$norm.factors   ## Normalization factors
    lsize <- countsDGE.nf$samples$lib.size  ## Library sizes

    # TMM Normalized Counts
    # From: http://grokbase.com/t/r/bioconductor/127r11kp23/bioc-edger-and-libsize-normalization
    # n.normalized <- 1e6 * (xte / (nf * lsize))
    n.normalized <- t(cpm(countsDGE.nf, log = FALSE))  # genes in the column, samples in the rows.

    nste <- n.normalized   ## Input data from transformed expression data.
    nste <- nste[-nrow(nste), ]   ## Remove reference sample from test data.
    sizeste <- nf[-length(nf)]

  } else if (type == "none"){


  }

  return(list(nste = nste, sizeste = sizeste))
}


# Poisson Dissimilarities for Clustering.
PoissonDistance <- function(x, beta = 1, type = c("mle", "deseq", "quantile"), transform = TRUE, alpha = NULL, perfeature = FALSE){
  type <- match.arg(type)
  if (!transform && !is.null(alpha)) stop("You have asked for NO transformation but have entered alpha.")
  if (transform && !is.null(alpha)){
    if (alpha > 0 && alpha <= 1) x <- x^alpha
    if (alpha <= 0 || alpha>1) stop("alpha must be between 0 and 1")
  }

  if (transform && is.null(alpha)){
    alpha <- FindBestTransform(x)
    x <- x^alpha
  }
  dd <- matrix(0, nrow=nrow(x), ncol=nrow(x))
  ddd <- NULL

  if (perfeature) ddd <- array(0, dim=c(nrow(x), nrow(x), ncol(x)))

  for (i in 2:nrow(dd)){
    xi <- x[i,]

    for (j in 1:(i-1)){
      xj <- x[j,]
      n <- NullModel(x[c(i,j),],type=type)$n
      ni <- n[1,]
      nj <- n[2,]
      di <- (xi+beta)/(ni+beta)
      dj <- (xj+beta)/(nj+beta)
      dd[i,j] <- sum(ni+nj-ni*di-nj*dj+xi*log(di)+xj*log(dj))

      if (perfeature) ddd[j,i,] <- ddd[i,j,] <- ni+nj-ni*di-nj*dj+xi*log(di)+xj*log(dj)
    }
  }

  save <- list(dd=as.dist(dd+t(dd)), alpha=alpha, x=x, ddd=ddd, alpha=alpha, type=type)
  class(save) <- "poidist"
  return(save)
}

print.poidist <- function(x,...){
  if(!is.null(x$alpha)) cat("Value of alpha used to transform data: ", x$alpha, fill = TRUE)
  if(is.null(x$alpha)) cat("No transformation performed.", fill = TRUE)
  cat("This type of normalization was performed:", x$type, fill = TRUE)
  cat("Dissimilarity computed for ", nrow(x$x), " observations.", fill = TRUE)
}




######### 2. HELPER FUNCTIONS   ##########
# 02/07/11 -- Helper functions.  (Copied from PoiClaClu source files.)

GoodnessOfFit <- function(x, type){
  ns <- NullModel(x, type = type)$n
  return(sum(((x - ns)^2)/ns, na.rm=TRUE))
}



ColorDendrogram <- function(hc, y, main = "", branchlength = 0.7, labels = NULL, xlab = NULL, sub = NULL, ylab = "", cex.main = NULL){
  if (is.null(labels))
    labels <- rep("", length(y))
  plot(hc, hang = 0.2, main = main, labels = labels, xlab = xlab,sub = sub, ylab = ylab, cex.main = cex.main)
  y <- y[hc$ord]
  if (is.numeric(y)) {
    y <- y + 1
    y[y == 7] <- "orange"
  }
  for (i in 1:length(hc$ord)) {
    o = hc$merge[, 1] == -hc$ord[i] | hc$merge[, 2] == -hc$ord[i]
    segments(i, hc$he[o] - branchlength, i, hc$he[o], col = y[i])
  }
}

permute.rows <- function(x){
  dd <- dim(x)
  n <- dd[1]
  p <- dd[2]
  mm <- runif(length(x)) + rep(seq(n) * 10, rep(p, n))
  matrix(t(x)[order(mm)], n, p, byrow = T)
}

#' @importFrom caret createFolds
foldIndex <- function(data = NULL, n = NULL, nFolds = 5, repeats = 2){

  if (!is.null(data)){
    n = nrow(data)
  }

  indIn <- indOut <- list()

  for (j in 1:repeats){
    tmpIn = createFolds(y = 1:n, k = nFolds, list = TRUE, returnTrain = TRUE)
    tmpOut = lapply(tmpIn, function(x)c(1:n)[-x])

    indIn = c(indIn, tmpIn)
    indOut = c(indOut, tmpOut)
  }

  nms = paste(rep(paste("Fold", 1:nFolds, sep = ""), repeats),
              rep(paste(".Rep", 1:repeats, sep = ""), c(rep(nFolds, repeats))), sep = "")
  names(indIn) <- names(indOut) <- nms
  return(list(indexIn = indIn, indexOut = indOut))
}

# Create folds and determine the indices of samples in each fold.
balanced.folds <- function(y, nfolds = min(min(table(y)), 10)){
  totals <- table(y)
  fmax <- max(totals)
  nfolds <- min(nfolds, fmax)
  # makes no sense to have more folds than the max class size
  folds <- as.list(seq(nfolds))
  yids <- split(seq(y), y)
  # nice way to get the ids in a list, split by class
  ###Make a big matrix, with enough rows to get in all the folds per class
  bigmat <- matrix(NA, ceiling(fmax/nfolds) * nfolds, length(totals))
  for(i in seq(totals)) {
    bigmat[seq(totals[i]), i] <- sample(yids[[i]])
  }
  smallmat <- matrix(bigmat, nrow = nfolds) # reshape the matrix
  ### Now do a clever sort to mix up the NAs
  smallmat <- permute.rows(t(smallmat)) ### Now a clever unlisting
  x <- apply(smallmat, 2, function(x) x[!is.na(x)])
  if(is.matrix(x)){
    xlist <- list()
    for(i in 1:ncol(x)){
      xlist[[i]] <- x[,i]
    }
    return(xlist)
  }
  return(x)
}

## Standart error of a given vector of parameter estimation.
se <- function(vec){
  return(sd(vec)/sqrt(length(vec)))
}

Soft <- function(x, a){
  return(sign(x)*pmax(abs(x) - a, 0))
}

# Find d.kj estimates
GetD <- function(ns, x, y, rho, beta, rhos = NULL){
  if (!is.null(rho) && !is.null(rhos)) stop("do you want to use rho or rhos in GetD function???")
  if (is.null(rhos)){
    uniq <- sort(unique(y))
    ds <- matrix(1, nrow = length(uniq), ncol = ncol(x))
    for (k in 1:length(uniq)){
      a <- colSums(x[y == uniq[k], ]) + beta
      b <- colSums(ns[y == uniq[k], ]) + beta
      ds[k, ] <- 1 + Soft(a/b - 1, rho/sqrt(b))
    }
    return(ds)
  } else {
    uniq <- sort(unique(y))
    ds.list <- list()
    for(rho in rhos){
      ds <- matrix(1, nrow = length(uniq), ncol = ncol(x))
      for(k in 1:length(uniq)){
        a <- colSums(x[y == uniq[k], ]) + beta
        b <- colSums(ns[y == uniq[k],]) + beta
        ds[k,] <- 1 + Soft(a/b - 1, rho/sqrt(b))
      }
      ds.list[[which(rhos == rho)]] <- ds
    }
    return(ds.list)
  }
}


# print.poicla <- function(x, ...){
#   if (is.null(x$rhos) && is.null(x$rho)){
#     if (!is.null(x[[1]]$alpha)) cat("Value of alpha used to transform data: ", x[[1]]$alpha, fill=TRUE)
#     if (is.null(x[[1]]$alpha)) cat("No transformation performed.", fill=TRUE)
#     cat("Type of normalization performed: ", x[[1]]$type, fill=TRUE)
#     cat("Number of training observations: ", nrow(x[[1]]$x), fill=TRUE)
#     if(!is.null(x[[1]]$xte)) cat("Number of test observations: ", nrow(x[[1]]$xte), fill=TRUE)
#     cat("Number of features: ", ncol(x[[1]]$x), fill=TRUE)
#     cat("-------------------------", fill=TRUE)
#     cat("-------------------------", fill=TRUE)
#     for (i in 1:length(x)){
#       cat("-------------------------", fill=TRUE)
#       cat("Value of rho used: ", x[[i]]$rho, fill=TRUE)
#       cat("Number of features used in classifier: ", sum(apply(x[[i]]$ds!=1, 2, sum)!=0), fill=TRUE)
#       if (sum(apply(x[[i]]$ds!=1, 2, sum)!=0)<100) cat("Indices of features used in classifier: ", which(apply(x[[i]]$ds!=1, 2, sum)!=0), fill=TRUE)
#     }
#   } else {
#     if (!is.null(x$alpha)) cat("Value of alpha used to transform data: ", x$alpha, fill=TRUE)
#     if (is.null(x$alpha)) cat("No transformation performed.", fill=TRUE)
#     cat("Type of normalization performed: ", x$type, fill=TRUE)
#     cat("Number of training observations: ", nrow(x$x), fill=TRUE)
#     if (!is.null(x$xte)) cat("Number of test observations: ", nrow(x$xte), fill=TRUE)
#     cat("Number of features: ", ncol(x$x), fill=TRUE)
#     cat("Value of rho used: ", x$rho, fill=TRUE)
#     cat("Number of features used in classifier: ", sum(apply(x$ds!=1, 2, sum)!=0), fill=TRUE)
#     if (sum(apply(x$ds!=1, 2, sum)!=0)<100) cat("Indices of features used in classifier: ", which(apply(x$ds!=1, 2, sum)!=0), fill=TRUE)
#   }
# }



