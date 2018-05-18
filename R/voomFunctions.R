##### HELPER FUNCTIONS ######
# Weighted statistics:
weighted.stats <- function(x, w, conditions){
  n = ncol(x) #number of samples
  p = nrow(x) #number of genes
  nclass = length(unique(conditions)) #number of class

  if (is.factor(conditions)){
    cNames = sort(levels(conditions))
  }
  if (is.numeric(conditions)){
    cNames = as.character(sort(unique(conditions)))
  }

  WM = WS = wSum = se.scale = matrix(0, p, nclass)
  rownames(WS) = rownames(WM) = rownames(wSum) = rownames(se.scale) = rownames(x)
  colnames(WS) = colnames(WM) = colnames(wSum) = colnames(se.scale) = cNames

  c.ind = as.numeric(conditions)

  w.mean00 = function (x, w){
    wm = NULL
    for (i in 1:p){
      wm0 = sum(w[i,]*x[i,]) / sum(w[i,])
      wm = c(wm, wm0)
    }
    return(wm)
  }

  w.mean = function (x, w, conditions){
    for (j in 1:nclass){
      WM[,j] = w.mean00(x[,c.ind == j], w[,c.ind == j])
    }
    return(WM)
  }

  w.sd = function (x, w, conditions){
    w.sd00 = function (x, w){
      ws = NULL
      w.sd0 = function (x, w){
        sumw = sum(w)
        sumw.sq = sum(w)^2
        w.sq = sum(w^2)
        denom = sum(w * ((x - mean(x))^2))
        sqrt((sumw * denom) / (sumw.sq - w.sq))
      }
      for (i in 1:p){
        ws0 = w.sd0(x[i,], w[i,])
        ws = c(ws, ws0)
      }
      return(ws)
    }
    for (j in 1:nclass){
      WS[,j] = w.sd00(x[,c.ind == j], w[,c.ind == j])
    }
    return(WS)
  }

  weightedMean = w.mean00(x, w) #Overall weighted mean
  weightedMean.C = w.mean(x, w, conditions) #Weighted means for each group
  weightedSD.C = w.sd(x, w, conditions) #Weighted standard deviations for each group
  #weightedSD.pooled = weightedSD.C


  for (i in 1:nclass){
    tmp = w[,which(c.ind == i)]
    rSum.tmp = rowSums(tmp)

    wSum[,i] = rSum.tmp
  }

  weightedSD.pooled = sqrt(rowSums((wSum - 1) * (weightedSD.C^2)) / (rowSums(wSum) - nclass))
  se.scale = sqrt(1 / wSum - 1 / rowSums(wSum))

  s0 = median(weightedSD.pooled)

  delta = (weightedMean.C - weightedMean)/(se.scale*(weightedSD.pooled + s0))


  weightedSD.pooled = sqrt(rowSums(as.data.frame(weightedSD.pooled)) / (n - nclass))
  stats = list(n = n, p = p, nclass = nclass, se.scale = se.scale, weightedMean = weightedMean, weightedMean.C = weightedMean.C, weightedSD.C = weightedSD.C, weightedSD.pooled = weightedSD.pooled, delta = delta)

  return(stats)
}

# soft.shrink(...)
soft.shrink = function(delta, threshold){
  dif = abs(delta) - threshold
  delta = sign(delta) * dif * (dif > 0)
  nonzero = sum(drop((dif > 0) %*% rep(1, ncol(delta))) > 0)
  attr(delta, "nonzero") = nonzero
  delta
}

# diag.disc(...)
diag.disc = function(x, centroids, prior, weight) {
  if (!missing(weight)) {
    posid = (weight > 0)
    if (any(posid)) {
      weight = sqrt(weight[posid])
      centroids = centroids[posid, , drop = FALSE] * weight
      x = x[posid, , drop = FALSE] * weight
    } else {
      mat = outer(rep(1, ncol(x)), log(prior), "*")
      dimnames(mat) = list(NULL, dimnames(centroids)[[2]])
      return(mat)
    }
  }
  dd = t(x) %*% centroids
  dd0 = drop(rep(1, nrow(centroids)) %*% (centroids^2))/2 - log(prior)
  names(dd0) = NULL
  scale(dd, dd0, FALSE)
}

# safe.exp(...)
safe.exp = function(x){
  xx = sign(x) * pmin(abs(x), 500)
  return(exp(xx))
}

# softmax(...)
softmax = function(x, gap = FALSE){
  d = dim(x)
  maxdist = x[, 1]
  pclass = rep(1, d[1])
  for (i in seq(2, d[2])) {
    l = x[, i] > maxdist
    pclass[l] = i
    maxdist[l] = x[l, i]
  }
  dd = dimnames(x)[[2]]
  if (gap) {
    x = abs(maxdist - x)
    x[cbind(seq(d[1]), pclass)] = drop(x %*% rep(1, d[2]))
    gaps = do.call("pmin", data.frame(x))
  }
  pclass <- if (is.null(dd) || !length(dd)){
    pclass
  } else {
    factor(pclass, levels = seq(d[2]), labels = dd)
  }
  if (gap){
    list(class = pclass, gaps = gaps)
  } else {
    pclass
  }
}



######### 1. voomDLDA & voomDQDA functions:   ###########
# control arguements for voom classifier.

#' Define controlling parameters for voom-based classifiers
#'
#' This function sets the control parameters for voom based classifiers while training the model.
#'
#' @param method validation method. Support repeated cross validation only ("repeatedcv").
#' @param number a positive integer. Number of folds.
#' @param repeats a positive integer. Number of repeats.
#' @param tuneLength a positive integer. If there is a tuning parameter in the classifier, this value
#' is used to define total number of tuning parameter to be searched.
#'
#'
#' @author Gokmen Zararsiz, Dincer Goksuluk, Selcuk Korkmaz, Vahap Eldem, Bernd Klaus, Ahmet Ozturk and Ahmet Ergun Karaagaoglu
#'
#' @keywords RNA-seq classification
#'
#' @seealso \code{\link{classify}}, \code{\link[caret]{trainControl}}, \code{\link{discreteControl}}
#'
#' @examples
#' 1L
#'
#' @name voomControl
#' @rdname voomControl
#'
#' @export
voomControl <- function(method = "repeatedcv", number = 5, repeats = 10, tuneLength = 10){
  list(method = method, number = number,
       repeats = repeats, tuneLength = tuneLength, controlClass = "voom.control")
}

#Calculation of quantile normalization factors necessary for calcNormFactorsGSD function. Same for both train and test sets.
calcFactorQuantileGSD <- function(data, lib.size, p = 0.75){
  y <- t(t(data)/lib.size)
  f <- apply(y, 2, function(x) quantile(x, p = p))
}

#Calculation of weighted normalization factors necessary for calcNormFactorsGSD function.
calcFactorWeightedGSD <- function(obs, ref, libsize.obs = NULL, libsize.ref = NULL, logratioTrim = 0.3,
                                   sumTrim = 0.05, doWeighting = TRUE, Acutoff = -1e+10){
  if (all(obs == ref)){
    return(1)
  }

  obs <- as.numeric(obs)
  ref <- as.numeric(ref)

  if (is.null(libsize.obs)){
    nO <- sum(obs)
  } else {
    nO <- libsize.obs
  }

  if (is.null(libsize.ref)){
    nR <- sum(ref)
  } else {
    nR <- libsize.ref
  }

  logR <- log2((obs/nO)/(ref/nR))
  absE <- (log2(obs/nO) + log2(ref/nR))/2
  v <- (nO - obs)/nO/obs + (nR - ref)/nR/ref
  fin <- is.finite(logR) & is.finite(absE) & (absE > Acutoff)
  logR <- logR[fin]
  absE <- absE[fin]
  v <- v[fin]
  n <- sum(fin)
  loL <- floor(n * logratioTrim) + 1
  hiL <- n + 1 - loL
  loS <- floor(n * sumTrim) + 1
  hiS <- n + 1 - loS
  keep <- (rank(logR) >= loL & rank(logR) <= hiL) & (rank(absE) >= loS & rank(absE) <= hiS)

  if (doWeighting){
    2^(sum(logR[keep]/v[keep], na.rm = TRUE)/sum(1/v[keep],na.rm = TRUE))
  } else {
    2^(mean(logR[keep], na.rm = TRUE))
  }
}

#Calculation of deseq normalization factors necessary for calcNormFactorsGSD function.
calcFactorRLEGSD <- function(data.train, data.test, lib.size, lib.size.test){
  gm <- exp(rowMeans(log(data.train)))
  f = apply(data.train, 2, function(u) median((u/gm)[gm > 0]))
  f.test = apply(data.test, 2, function(u) median((u/gm)[gm > 0]))
  f = f / lib.size
  f.test = f.test / lib.size.test
  deseqsizefactors = list(f, f.test)
  return(deseqsizefactors)
}

#Calculation of normalization factors using calcNormFactorsGSD function.
calcNormFactorsGSD <- function(data.train, data.test, lib.size = NULL, method = c("TMM", "deseq", "none"), refColumn = NULL, logratioTrim = 0.3, sumTrim = 0.05,
                                doWeighting = TRUE, Acutoff = -1e+10, p = 0.75, ...){
  x <- as.matrix(data.train)
  xtest <- as.matrix(data.test)
  if (any(is.na(x)||is.na(xtest)))
    stop("NAs not permitted")
  if (is.null(lib.size))
    lib.size <- colSums(x)

  lib.size.test <- colSums(xtest)
  method <- match.arg(method)
  allzero <- rowSums(x > 0) == 0
  if (any(allzero)){
    x <- x[!allzero, , drop = FALSE]
  }

  xtest <- xtest[!allzero, , drop = FALSE]
  if (nrow(x) == 0 || ncol(x) == 1){
    method = "none"
  }

  if (method == "TMM"){
    f75 <- calcFactorQuantileGSD(data = x, lib.size = lib.size, p = 0.75)
    f75.test <- calcFactorQuantileGSD(data = xtest, lib.size = lib.size.test, p = 0.75)

    refColumn <- which.min(abs(f75 - mean(f75)))
    f <- rep(NA, ncol(x))
    f.test <- rep(NA, ncol(xtest))
    for (i in 1:ncol(x)) f[i] <- calcFactorWeightedGSD(obs = x[,i], ref = x[, refColumn], libsize.obs = lib.size[i],
                                                       libsize.ref = lib.size[refColumn], logratioTrim = logratioTrim,
                                                       sumTrim = sumTrim, doWeighting = doWeighting, Acutoff = Acutoff)
    for (i in 1:ncol(xtest)) f.test[i] <- calcFactorWeightedGSD(obs = xtest[,i], ref = x[, refColumn], libsize.obs = lib.size.test[i],
                                                                libsize.ref = lib.size[refColumn], logratioTrim = logratioTrim,
                                                                sumTrim = sumTrim, doWeighting = doWeighting, Acutoff = Acutoff)
    normf = list(f,f.test)
  } else if (method == "deseq"){
    normf = calcFactorRLEGSD(data.train = x, data.test = xtest, lib.size = lib.size, lib.size.test = lib.size.test)#/lib.size
  } else {
    normf = list(rep(1, ncol(x)), rep(1, ncol(xtest)))
  }

  names(normf) = c("train", "test")

  f = as.numeric(normf[[1]]) / (exp(mean(log(normf[[1]]))))
  f.test = as.numeric(normf[[2]]) / (exp(mean(log(normf[[1]]))))
  normf2 = list(f, f.test, lib.size, lib.size.test)
  names(normf2) = c("TrainNormFactor","TestNormFactor","TrainLibSize","TestLibSize")
  return(normf2)
}

#' @importFrom stats lowess approxfun
#' @importFrom limma lmFit
voomGSD <- function(data.train, data.test, group, norm.method = c("TMM", "deseq", "none"), design = NULL, lib.size = NULL, span = 0.5){
  #traindatayÄ±, testdatayÄ± ve trainclassÄ± alacak, method = TMM, RLE, none olacak.
  #train ve test icin voom gibi expression ve weights dÃ¶ndÃ¼recek.
  out <- list()
  NormFactors = calcNormFactorsGSD(data.train = data.train, data.test = data.test, method = norm.method)
  TrainNormFactor = NormFactors$TrainNormFactor
  TestNormFactor = NormFactors$TestNormFactor
  TrainLibSize = NormFactors$TrainLibSize
  TestLibSize = NormFactors$TestLibSize
  lib.size.tr = TrainNormFactor * TrainLibSize
  lib.size.ts = TestNormFactor * TestLibSize

  design.tr = model.matrix(~group)
  rownames(design.tr) = colnames(data.train)

  design.ts <- matrix(1, ncol(data.test), 1)
  rownames(design.ts) <- colnames(data.test)
  colnames(design.ts) <- "GrandMean"

  y.tr <- t(log2(t(data.train + 0.5)/(lib.size.tr + 1) * 1e+06))
  y.ts <- t(log2(t(data.test + 0.5)/(lib.size.ts + 1) * 1e+06))
  fit.tr <- lmFit(y.tr, design.tr)
  fit.ts <- lmFit(y.ts, design.ts)


  if (is.null(fit.tr$Amean))
    fit$Amean <- rowMeans(y.tr, na.rm = TRUE)

  fit.ts$Amean = fit.tr$Amean
  fit.ts$sigma = fit.tr$sigma
  fit.ts$coefficients = fit.tr$coefficients[,1]

  sx <- fit.tr$Amean + mean(log2(lib.size.tr + 1)) - log2(1e+06)
  sy <- sqrt(fit.tr$sigma)
  l <- lowess(sx, sy, f = span)
  f <- approxfun(l, rule = 2)

  fitted.values.tr <- fit.tr$coefficients %*% t(fit.tr$design)
  fitted.values.ts <- fit.ts$coefficients %*% t(fit.ts$design)
  fitted.cpm.tr <- 2^fitted.values.tr
  fitted.cpm.ts <- 2^fitted.values.ts
  fitted.count.tr <- 1e-06 * t(t(fitted.cpm.tr) * (lib.size.tr + 1))
  fitted.count.ts <- 1e-06 * t(t(fitted.cpm.ts) * (lib.size.ts + 1))
  fitted.logcount.tr <- log2(fitted.count.tr)
  fitted.logcount.ts <- log2(fitted.count.ts)
  w.tr <- 1/f(fitted.logcount.tr)^4
  w.ts <- 1/f(fitted.logcount.ts)^4
  dim(w.tr) <- dim(fitted.logcount.tr)
  dim(w.ts) <- dim(fitted.logcount.ts)
  dimnames(w.tr) = dimnames(y.tr)
  dimnames(w.ts) = dimnames(y.ts)
  out$TrainExp <- y.tr
  out$TestExp <- y.ts
  out$TrainWeights <- w.tr
  out$TestWeights <- w.ts
  new("EList", out)
}


# input <- counts(data)
# output <- colData(data)[ ,class.labels]
# pooled.var = ifelse(method == "voomDLDA", TRUE, FALSE)


# input: raw counts extracted from 'DESeqDataSet' object.
# output: sınıf labellarını içeren bir vector.
voomDDA.train = function(input, output, pooled.var = ifelse(method == "voomDLDA", TRUE, FALSE),
                         normalize = normalize, class.labels = NULL){

  rawCounts <- input
  countsDGE <- DGEList(counts = input, genes = rownames(input))

  normalizedCounts <- function(x, type){
    switch(type,
           deseq = calcNormFactors(x, method = "RLE"),
           TMM = calcNormFactors(x, method = "TMM"),
           none = calcNormFactors(x, method = "none"))
  }

  countsDGE.normalized <- normalizedCounts(x = countsDGE, type = normalize)   ## RLE: DESeq mantigi ile normalize ediyor.
  rawData <- countsDGE
  rawData[[class.labels]] <- output

  ### Design matrisi daha sonra düzenlenecek.
  # des <- model.matrix(~ rep(1, length(rawData[[class.labels]])) + colData(data)[ ,"treat"])[ ,1, drop = FALSE]

  voomRes = voom(counts = countsDGE.normalized, plot = F)
  x = voomRes$E # pxn dim gene expression matrix  - logCPM values
  w = voomRes$weights #pxn dim weight matrix - Precision Weigts
  dimnames(w) = dimnames(x)

  transformedData <- list(normalizedData = countsDGE.normalized, transformedData = voomRes)
  transformedData[[class.labels]] <- output

  # input <- list(x = x, w = w)   ## Input data from transformed expression data.
  # output <- output   ## Output: class labels of samples.

  trainParameters <- list(NULL)

  p = nrow(x)
  n = ncol(x)

  nclass = length(unique(output))
  n.class = table(output)

  if (min(n.class) == 1) {
    stop(warning("Warning: a class contains only 1 sample"))
  }

  if (is.factor(output)){
    classNames = levels(output)
  }

  wStats = weighted.stats(x = x, w = w, output)

  results = structure(list(counts = rawCounts, input = list(x = x, w = w), output = factor(output),
                           weightedStats = wStats, normalization = normalize,
                           nclass = nclass, classNames = classNames, PooledVar = pooled.var,
                           rawData = countsDGE, transformedData = transformedData), class = "voomDDA")

  return(results)
}


#' @importFrom stats napredict
predict.voomDDA = function(object, newdata){
  ## Object: "voomDDA.train" fonksiyonundan dönen bir obje olacak.
  ## newdata: features in the rows, samples in the columns. count datalar.
  ##
  n = ncol(newdata)
  p = nrow(newdata)

  disc = matrix(0, n, object$nclass)  ## Discriminant scores for each class
  dimnames(disc) = list(colnames(newdata), object$classNames)
  vm = voomGSD(data.train = object$counts, data.test = newdata, group = object$output,
               norm.method = object$normalization)
  x = vm$TestExp

  x2 = t(x)

  ### Bu kısım ayrıca bir fonksiyon olarak düzenlenecek. calculateDiscriminantScores(...)
  {
    if (object$PooledVar) {
      vp = (object$weightedStats$weightedSD.pooled)^2
      if (any(i0 <- vp == 0))
        vp[i0] <- 1e-07 * min(vp[!i0])
      ivp <- rep(1/vp, each = n)
      for (k in 1:(object$nclass)) {
        y = x2 - rep(object$weightedStats$weightedMean.C[, k], each = n)
        disc[, k] = rowSums(y * y * ivp)
      }
    } else {
      if (FALSE) {
        for (k in 1:(object$nclass)) {
          x2 = x2 - rep(object$weightedStats$weightedMean.C[, k], each = n)
          vsd = (object$weightedStats$weightedSD.C)^2
          disc[, k] = rowSums((x2 * x2)/rep(vsd[, k], each = n)) + sum(log(vsd[, k]))
        }
      } else {
        vsd = (object$weightedStats$weightedSD.C)^2
        for (k in 1:(object$nclass)) {
          disc[, k] = apply(x2, 1, function(z) sum((z - object$weightedStats$weightedMean.C[, k])^2/vsd[, k])) + sum(log(vsd[, k]))
        }
      }
    }
  }

  idx = apply(disc, 1, which.min)
  pred = colnames(disc)[idx]

  if (inherits(attr(x2, "na.action"), "exclude"))
    pred = napredict(omit = attr(x2, "na.action"), pred)
  pred

}


######### 2. voomNSC functions:   ###########
### voom NSC train

# counts = counts(data)
# conditions = colData(data)[ ,class.labels]
# normalization = "TMM"
# n.threshold = control$tuneLength
# thresholds = NULL
# getInitialThresholds = FALSE
# remove.zeros = TRUE
# offset.percent = 50
# idxIn = folds$indexIn$Fold1.Rep1

voomNSC.train = function(counts, conditions, n.threshold = 30, thresholds = NULL, offset.percent = 50,
                         remove.zeros = TRUE, normalization = c("TMM", "deseq", "none"), idxIn = NULL,
                         getInitialThresholds = FALSE, ...){

  y = as.factor(conditions)
  n.class = table(y)
  prior = n.class/length(y)

  #   weights = w
  ## wnsc(....)

  # y = conditions
  # offset.percent = offset.percent
  # n.threshold = n.threshold
  # prior = prior
  # remove.zeros = remove.zeros
  # weights = w
  # idx = idxIn

  wnsc = function(x, y, n.threshold = 30, offset.percent = 50, prior = NULL, remove.zeros = TRUE,
                  weights = NULL, idx = idxIn){
    selected.genes <- selected.genesIndex <- list()
    this.call = match.call()
    Y = model.matrix(~factor(y) - 1, data = list(y = y))   #### BURADAKİ DESIGN MATRİSİNİ DE KONTROL EDELİM

    if (!is.null(idx)){
      xtest = x[ ,-idx]
      ytest = y[-idx]
    } else {
      xtest <- x
      ytest <- y
    }

    wStats = weighted.stats(x, weights, y)

    if (min(n.class) == 1) {
      stop(warning("Warning: a class contains only 1 sample"))
    }

    # n = sum(n.class)
    ntest = ncol(xtest)
    K = length(prior)
    # p = nrow(x)
    dimnames(Y) = list(NULL, names(n.class))
    centroids = wStats$weightedMean.C
    sd = wStats$weightedSD.pooled
    offset = quantile(sd, offset.percent/100)
    sd = sd + offset
    centroid.overall = wStats$weightedMean
    se.scale = wStats$se.scale
    delta = wStats$delta

    if (is.null(thresholds)){
      threshold = seq(0, max(abs(delta)), length = n.threshold)
    } else {
      threshold <- thresholds
    }

    nonzero = numeric(length(threshold))
    errors = numeric(length(threshold))
    predictions = as.list(seq(n.threshold))
    prob = array(0, c(ntest, K, n.threshold))

    dshrunkAll <- list()
    for (ii in 1:n.threshold){
      delta.shrunk = soft.shrink(delta, threshold[ii])
      delta.shrunk = t(t(delta.shrunk) * as.numeric(se.scale))
      dshrunkAll[[ii]] <- delta.shrunk

      nonzero[ii] = attr(delta.shrunk, "nonzero")
      posid = drop(abs(delta.shrunk) %*% rep(1, K)) > 0
      selected.genes[[ii]] = rownames(wStats$weightedMean.C)[posid]
      selected.genesIndex[[ii]] = which(posid == TRUE)

      ##### BURADA KULLANILAN "xtest" MATRISI RAW COUNTS DEGERLERI ICERIYOR. "centroid.overall" DEGERI ISE
      ##### VOOM DONUSUMU SONUCUNDA ELDE EDILEN w ve x MATRISLERI YARDIMI ILE HESAPLANIYOR. BU SEBEPLE TEST
      ##### MATRISININ DE VOOM DONUSUMU ELDE EDILMIS HALI ILE KULLANILMASI GEREKIR. SONUCLAR BU SEBEPLE
      ##### YANLIS CIKIYOR OLABILIR.
      dd = diag.disc((xtest - centroid.overall)/sd, delta.shrunk,
                     prior, weight = posid)
      predictions[[ii]] = softmax(dd)
      dd = safe.exp(dd)
      prob[, , ii] = dd/drop(dd %*% rep(1, K))
      if (!is.null(ytest)) {
        errors[ii] = sum(predictions[[ii]] != ytest)
      }
    }

    thresh.names = format(round(threshold, 3))
    names(predictions) = names(selected.genes) = thresh.names
    attr(predictions, "row.names") = paste(seq(ntest))
    class(predictions) = "data.frame"

    if (remove.zeros){
      n.threshold = match(0, nonzero, n.threshold)
    }

    dimnames(prob) = list(paste(seq(ntest)), names(n.class), thresh.names)

    object = list(counts = counts, conditions = factor(conditions), predictions = predictions, prob = prob[, , seq(n.threshold)],
                  weightedMean.C = centroids, weightedMean = centroid.overall, delta = delta, normalization = normalization,
                  weightedSD.pooled = sd, threshold = threshold[seq(n.threshold)], nonzero = nonzero[seq(n.threshold)],
                  se.scale = se.scale, call = this.call, prior = prior, offset = offset, SelectedGenes = selected.genes,
                  SelectedGenesIndex = selected.genesIndex, weightedStats = wStats)

    # object$errors = errors
    object$accuracy = 1 - errors / ntest

    # opt.threshold = max(object$threshold[which(object$errors == min(object$errors))])
    opt.threshold = max(object$threshold[which(object$accuracy == max(object$accuracy))])

    model.res = as.data.frame(cbind(object$threshold, object$nonzero, object$accuracy))
    colnames(model.res) = c("threshold", "nonzero", "accuracy")
    object$modelRes = model.res

    object$delta.shrunk <- dshrunkAll[[which(object$threshold == opt.threshold)]]
    object$opt.threshold = opt.threshold
    object
  }

  ### Initial Calculations:
  if (!is.null(thresholds)){
    n.threshold <- length(thresholds)
  }

  normalization = match.arg(normalization)
  this.call = match.call()

  #### voom transformation steps
  if (normalization == "TMM"){
    design = model.matrix(~ conditions)   ##### Design sadece "conditions"a göre olmamali.
    rownames(design) = colnames(counts)
    dge = DGEList(counts = counts)
    dge = calcNormFactors(dge, method = "TMM")
    vm = voom(dge, design, plot = F)
    x = vm $ E #pxn dim gene expression matrix
    w = vm $ weights #pxn dim weight matrix
    dimnames(w) = dimnames(x)
  } else if (normalization == "deseq"){
    design = model.matrix(~ conditions)
    rownames(design) = colnames(counts)
    dge = DGEList(counts = counts)
    dge = calcNormFactors(dge, method = "RLE")
    vm = voom(dge, design, plot=F)
    x = vm $ E #pxn dim gene expression matrix
    w = vm $ weights #pxn dim weight matrix
    dimnames(w) = dimnames(x)
  } else {
    design = model.matrix(~ conditions)
    rownames(design) = colnames(counts)
    vm = voom(counts, design, plot = F)
    x = vm $ E #pxn dim gene expression matrix
    w = vm $ weights #pxn dim weight matrix
    dimnames(w) = dimnames(x)
  }

  if (getInitialThresholds){
    wStats = weighted.stats(x, w, y)
    delta = wStats$delta
    threshold = seq(0, max(abs(delta)), length = n.threshold)
    return(threshold)
  } else {
    junk = wnsc(x, y = conditions, offset.percent = offset.percent, n.threshold = n.threshold, prior = prior,
                remove.zeros = remove.zeros, weights = w)

    junk$call = this.call
    # class(junk) = "pamrtrained"   #### voom.train sınıfında bir obje olarak dönecek.
    return(junk)
  }
}


## voom NSC prediction
predict.voomNSC = function(fit, newdata, threshold = NULL, prior = NULL){
  # Args:
  #   fit: fitted model from voomNSC.train(...)
  #   newdata: (p x n) raw counts matrix for test set.
  #   threshold: a numeric value for threshold (tuning parameter).
  #              If NULL, it is obtained from fitted model.
  #   prior: a numeric value in the interval [0, 1] for prior probabilities of classes.
  #          If NULL, it is obtained from fitted model.

  if (is.null(prior)) prior = fit$prior
  if (is.null(threshold)) threshold = fit$opt.threshold
  vm = voomGSD(data.train = fit$counts, data.test = newdata, group = fit$conditions,
               norm.method = fit$normalization)
  x = vm$TestExp

  sd = fit$weightedSD.pooled
  centroid.overall = fit$weightedMean
  centroids = fit$weightedMean.C
  se.scale = fit$se.scale
  delta = fit$delta

  delta.shrunk = soft.shrink(delta, threshold)
  delta.shrunk = t(t(delta.shrunk) * as.numeric(se.scale))
  posid = drop(abs(delta.shrunk) %*% rep(1, length(prior))) > 0
  dd = diag.disc((x - fit$weightedMean)/(fit$weightedSD.pooled), delta.shrunk, prior, posid)
  softmax(dd)
}


