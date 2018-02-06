#' @importFrom quadprog solve.QP
projectCellType = function (Y, coefCellType, contrastCellType = NULL, nonnegative = TRUE,
          lessThanOne = FALSE)
{
  if (is.null(contrastCellType))
    Xmat <- coefCellType
  else Xmat <- tcrossprod(coefCellType, contrastCellType)
  nCol <- dim(Xmat)[2]
  if (nCol == 2) {
    Dmat <- crossprod(Xmat)
    mixCoef <- t(apply(Y, 2, function(x) {
      solve(Dmat, crossprod(Xmat, x))
    }))
    colnames(mixCoef) <- colnames(Xmat)
    return(mixCoef)
  }
  else {
    nSubj <- dim(Y)[2]
    mixCoef <- matrix(0, nSubj, nCol)
    rownames(mixCoef) <- colnames(Y)
    colnames(mixCoef) <- colnames(Xmat)
    if (nonnegative) {
      if (lessThanOne) {
        Amat <- cbind(rep(-1, nCol), diag(nCol))
        b0vec <- c(-1, rep(0, nCol))
      }
      else {
        Amat <- diag(nCol)
        b0vec <- rep(0, nCol)
      }
      for (i in 1:nSubj) {
        obs <- which(!is.na(Y[, i]))
        Dmat <- crossprod(Xmat[obs, ])
        mixCoef[i, ] <- solve.QP(Dmat, crossprod(Xmat[obs,
                                                      ], Y[obs, i]), Amat, b0vec)$sol
      }
    }
    else {
      for (i in 1:nSubj) {
        obs <- which(!is.na(Y[, i]))
        Dmat <- crossprod(Xmat[obs, ])
        mixCoef[i, ] <- solve(Dmat, t(Xmat[obs, ]) %*%
                                Y[obs, i])
      }
    }
    return(mixCoef)
  }
}

#' @importFrom stats model.matrix lm pf vcov
#' @importFrom nlme lme getVarCov
#'
validationCellType = function (Y, pheno, modelFix, modelBatch = NULL, L.forFstat = NULL,
          verbose = FALSE)
{
  N <- dim(pheno)[1]
  pheno$y <- rep(0, N)
  xTest <- model.matrix(modelFix, pheno)
  sizeModel <- dim(xTest)[2]
  M <- dim(Y)[1]
  if (is.null(L.forFstat)) {
    L.forFstat <- diag(sizeModel)[-1, ]
    colnames(L.forFstat) <- colnames(xTest)
    rownames(L.forFstat) <- colnames(xTest)[-1]
  }
  sigmaResid <- sigmaIcept <- nObserved <- nClusters <- Fstat <- rep(NA,
                                                                     M)
  coefEsts <- matrix(NA, M, sizeModel)
  coefVcovs <- list()
  if (verbose)
    cat("[validationCellType] ")
  for (j in 1:M) {
    ii <- !is.na(Y[j, ])
    nObserved[j] <- sum(ii)
    pheno$y <- Y[j, ]
    if (j%%round(M/10) == 0 && verbose)
      cat(".")
    try({
      if (!is.null(modelBatch)) {
        fit <- try(lme(modelFix, random = modelBatch,
                       data = pheno[ii, ]))
        OLS <- inherits(fit, "try-error")
      }
      else OLS <- TRUE
      if (OLS) {
        fit <- lm(modelFix, data = pheno[ii, ])
        fitCoef <- fit$coef
        sigmaResid[j] <- summary(fit)$sigma
        sigmaIcept[j] <- 0
        nClusters[j] <- 0
      }
      else {
        fitCoef <- fit$coef$fixed
        sigmaResid[j] <- fit$sigma
        sigmaIcept[j] <- sqrt(getVarCov(fit)[1])
        nClusters[j] <- length(fit$coef$random[[1]])
      }
      coefEsts[j, ] <- fitCoef
      coefVcovs[[j]] <- vcov(fit)
      useCoef <- L.forFstat %*% fitCoef
      useV <- L.forFstat %*% coefVcovs[[j]] %*% t(L.forFstat)
      Fstat[j] <- (t(useCoef) %*% solve(useV, useCoef))/sizeModel
    })
  }
  if (verbose)
    cat(" done\n")
  rownames(coefEsts) <- rownames(Y)
  colnames(coefEsts) <- names(fitCoef)
  degFree <- nObserved - nClusters - sizeModel + 1
  Pval <- 1 - pf(Fstat, sizeModel, degFree)
  out <- list(coefEsts = coefEsts, coefVcovs = coefVcovs, modelFix = modelFix,
              modelBatch = modelBatch, sigmaIcept = sigmaIcept, sigmaResid = sigmaResid,
              L.forFstat = L.forFstat, Pval = Pval, orderFstat = order(-Fstat),
              Fstat = Fstat, nClusters = nClusters, nObserved = nObserved,
              degFree = degFree)
  out
}
