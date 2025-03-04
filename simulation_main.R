library(fdapace)
##' @title generate sample path X for each subject based on observation times
##' @param eigenVals a vector with length k
##' @param t observation times, a value or a vector from [0, 1]
##' @param eigenFuns default - Fourier basis functions
##' @param Gau if scores are normal distributed. True: N(0, lambda); False: exp(1/sqrt(lambda_i)) - sqrt(lambda)
##' @return a sample path with the same length as input t on [0, 1]
getX <- function(eigenVals, t = NULL, eigenFuns = NULL, Gau = TRUE){
  stopifnot(length(eigenVals) > 0, length(t) >= 1)
  mat <- sapply(1:length(eigenVals), function(k.){
    theScore <- ifelse(Gau, rnorm(n=1, mean=0, sd=sqrt(eigenVals[k.])),
                       rexp(n=1, rate=1/sqrt(eigenVals[k.])) - sqrt(eigenVals[k.]))
    # fourier basis
    if(k. %% 2 == 1){
      if(k. == 1){
        rep(theScore, times = length(t))
      } else{
        theScore*cos((k. - 1)*t*pi)*sqrt(2)
      }
    } else{
      theScore*sin(k.*t*pi)*sqrt(2)
    }
  })
  if(length(t) == 1){
    return(sum(mat))
  } else{
    if(length(eigenVals) == 1){
      return(mat)
    } else{
      return(rowSums(mat))
    }
  }
}


##' @title Generate discrete observations for EACH subject
##' @param N is observation times
##' @param obsGrid is the location of observations.
##' @param errSig2ma the variance of errors
##' @return Return a N by 6 matrix, columns are t, mu, ranX, trueY, errors, y.
datGP <- function(obsGrid, N, eigenVals, errSig2ma, ...){
  muTrue <- obsGrid + 10*exp(-(obsGrid-0.5)^2)
  relizationX <- getX(eigenVals = eigenVals, t = obsGrid, ...)
  trueY <- muTrue + relizationX
  errors <- rnorm(N, mean = 0, sd = sqrt(errSig2ma))
  y <- trueY + errors
  return(data.frame(t = obsGrid, mu = muTrue, ranX = relizationX, trueY = trueY, errors = errors, y = y))
}

##' @title select optimal bandwidth by the generalized cross-validation score.
##' @param const the constant depends on the kernel
##' @param t.in/y.in observation grid/responses, must be sorted with respect to t.in
##' @param bw.min minimal bandwidth
##' @param t.out where the function should be fitted, only for 2d case to reduce computational cost
##' @return Return a optimal bandwidth
bwSelection.GCV <- function(const, t.in, y.in, t.range, bw.min, t.out,
                            d = c(1, 2))
{
  stopifnot(d %in% c(1, 2))
  totObs <- length(y.in)
  w.in <- rep(1, totObs)
  if(d == 1L){
    stopifnot(is.atomic(t.in)) # t must be a vector
    q <- (0.25*t.range/bw.min)^(1/9)
    bwCandidates = q^(0:9)*bw.min
    gcvScores <- sapply(bwCandidates, function(h.) {
      newmu = Lwls1D(h., kernel_type='gauss', npoly = 1, nder = 0, xin = t.in,
                     yin = y.in, win = w.in, xout = t.in)
      stopifnot(sum(is.na(newmu)) == 0)
      cvsum = sum((newmu - y.in)^2)
      return(cvsum/(1-(t.range*const)/(totObs*h.))^2)
    })
    return(bwCandidates[which.min(gcvScores)])
  } else{
    stopifnot(!is.null(dim(t.in))) # require a 2d domain input t.in
    stopifnot(ncol(t.in) == 2) # only support 2d input
    t.range <- t.range*sqrt(2) # adjust t.range for 2d
    q <- (0.25*t.range/bw.min)^(1/9)
    bwCandidates = q^(0:9)*bw.min
    if(!is.matrix(t.in)) t.in <- as.matrix(t.in)
    gcvScores <- rep(Inf, length(bwCandidates))
    counter <- 0
    # get GCV from largest bandwidth
    for (i. in rev(seq_along(bwCandidates))) {
      h. <- bwCandidates[i.]
      # 2d fitting can not have more than 201 observations per axis. Otherwise it will dead
      fit <- Lwls2D(h., 'gauss', xin = t.in, yin = y.in, win = w.in,
                    xout1 = t.out, xout2 = t.out)
      stopifnot(sum(is.na(fit)) == 0)
      obsFit <- fdapace:::interp2lin(t.out, t.out, fit, t.in[, 1], t.in[, 2])
      stopifnot(sum(is.na(obsFit)) == 0)
      rss <- sum((y.in - obsFit)^2)
      bottom <- max(1 - 3*(1/totObs)*(t.range*const/h.)^2, 0)
      temp <- rss/bottom^2
      cond1 <- is.infinite(temp)
      cond2 <- temp > min(gcvScores) & (temp - min(gcvScores))/min(gcvScores) > 0.05
      if(cond1 | cond2){
        break
      }
      gcvScores[i.] <- temp
      counter <- counter + 1
    }
    if (counter == 0) {
      stop('data is too sparse. Need larger bandwidth candidates')
    } else{
      return(bwCandidates[which.min(gcvScores)])
    }
  }
}

##' @title conduct functional PCA
##' @param df data input, in data.frame
##' @param max_K the maximum searching range
##' @param bw.mu/bw.err/bw.cov bandwidth for estimating mean/error/covariance
##' @param estMethod smooth for independent design and sample average for dense common design
fPCA.rev <- function(df, fPCAinput, max_K = 5, bw.mu = NA, bw.err = NA,
                     bw.cov = NA, estMethod = NULL, need.score = TRUE)
{
  df.sort <- df[order(df$t), ]
  idx <- unlist(fPCAinput$Lid)
  N.seq <- sapply(idx, function(n.) nrow(subset(df.sort, id == n.)))
  obsGrid <- subset(df.sort, select = 't', drop = TRUE)
  y.sort <- subset(df.sort, select = 'y', drop = TRUE)
  obs.range <- range(obsGrid)
  buffer <- .Machine$double.eps*10 # ~ 1e-15 scale
  
  if(estMethod == 'cross-sectional'){
    regGrid <- subset(df.sort, id == min(idx), select = t, drop = TRUE)
    if (abs(obsGrid[1] - regGrid[1]) < buffer) 
      obsGrid[1] <- regGrid[1]
    if (abs(obsGrid[length(obsGrid)] - regGrid[length(regGrid)]) < buffer) 
      obsGrid[length(obsGrid)] <- regGrid[length(regGrid)]
    
    ymat.t <- sapply(idx, function(n.) subset(df, id == n., select = y, drop = TRUE)) # N by n matrix
    ymat <- t(ymat.t) # n by N matrix
    muhat.reg <- colMeans(ymat) # cross-sectional mean estimation
    rawCov <- cov(ymat, use = 'pairwise.complete.obs') # cross-sectional covariance estimation
    rawCov <- 0.5*(rawCov + t(rawCov))
    sigma2hat <- mean(diff(ymat.t, differences=2)^2, na.rm=TRUE)/choose(4, 2)
    covhat.reg <- rawCov
    diag(covhat.reg) <- diag(covhat.reg) - sigma2hat
  } else{
    rational_cut <- c(0.1, 0.9)
    k0 <- 1/sqrt(2*pi) # constant for Gaussian kernel
    maxOut <- 51 # the number is to reduce computational cost in 2d fitting
    regGrid <- seq(obs.range[1], obs.range[2], length.out = maxOut)
    
    if (regGrid[1] < obsGrid[1] & abs(obsGrid[1] - regGrid[1]) < buffer){
      obsGrid[1] <- regGrid[1]
    }
    if (regGrid[maxOut] > obsGrid[length(obsGrid)] &
        abs(obsGrid[length(obsGrid)] - regGrid[maxOut]) < buffer){
      obsGrid[length(obsGrid)] <- regGrid[maxOut]
    }
    
    if(is.na(bw.mu)){
      # start GCV bandwidth selection for mean estimation
      k.mu <- 3 # minimal number of neighbors needed to decide minimum bandwidth for mu estimation
      d.1d <- dist(obsGrid)
      # obtain the minimum bandwidth by GCV to approximate mu
      bwmu.min <- max(dbscan::kNNdist(d.1d, k.mu, approx = 0.05)) # add approx to speed up
      bw.mu <- bwSelection.GCV(const = k0, t.in = obsGrid, y.in = y.sort,
                               t.range = diff(obs.range), bw.min = bwmu.min,
                               t.out = regGrid, d = 1)
    }
    muhat.obs <- Lwls1D(bw.mu, kernel_type = 'gauss', npoly = 1, nder = 0,
                        xin = obsGrid, yin= y.sort, win = rep(1, length(y.sort)),
                        xout = obsGrid)
    muhat.reg <- Lwls1D(bw.mu, kernel_type = 'gauss', npoly = 1, nder = 0,
                        xin = obsGrid, yin= y.sort, win = rep(1, length(y.sort)),
                        xout = regGrid)
    df.sort$muhat <- muhat.obs
    df.sort$resihat <- df.sort$y - muhat.obs
    
    # this will break the ordering of t
    rcov <- lapply(idx, function(n.){
      temp <- subset(df.sort, id == n.)
      if(nrow(temp) > 1){
        t2d <- expand.grid(x=temp$t, y=temp$t)
        resiprod <- matrix(temp$resihat) %*% temp$resihat
        return(cbind(sub.id = n., t = t2d, resiprod = c(resiprod)))
      }
    })
    rcov <- do.call(rbind, rcov)
    stopifnot(nrow(rcov[which(is.na(rcov$resiprod)), ]) == 0)
    cov.diag <- subset(rcov, t.x == t.y)
    cov.offd <- subset(rcov, t.x != t.y)
    
    tmat.neq <- as.matrix(subset(cov.offd, select = c('t.x', 't.y')))
    resimat.neq <- subset(cov.offd, select = 'resiprod', drop = TRUE)
    if(!is.vector(resimat.neq)) resimat.neq <- unlist(resimat.neq)
    
    if(is.na(bw.cov)){
      k.cov <- 6 # minimal number of neighbors needed to decide minimum bandwidth for cov estimation
      bwcov.min <- max(dbscan::kNNdist(tmat.neq, k.cov, approx = 0.05))
      bw.cov <- bwSelection.GCV(const = k0, t.in = tmat.neq, y.in = resimat.neq,
                                t.range = diff(obs.range), bw.min = bwcov.min,
                                t.out = regGrid, d=2)
    }
    
    covhat.reg <- Lwls2D(bw.cov, 'gauss', xin = tmat.neq, yin = resimat.neq,
                         win = rep(1, length(resimat.neq)), xout1 = regGrid,
                         xout2 = regGrid)
    
    estCov.diag <- spline(x = regGrid, y = diag(covhat.reg), xout = obsGrid)$y
    # ordering the diagonal entries for 1d smoothing
    cov.diag <- cov.diag[order(cov.diag$t.x), ]
    # stopifnot(all.equal(cov.diag$t.x, obsGrid))
    nDiag <- nrow(cov.diag)
    # 1d smoother for the diagonals of cov(t,s)
    # use GCV to select bandwidth for estimating variance of errors.
    
    if(is.na(bw.err)){
      # use bwmu.min to reduce computational cost
      bw.err <- bwSelection.GCV(const = k0, t.in = obsGrid, y.in = cov.diag$resiprod,
                                t.range = diff(obs.range), bw.min = bwmu.min,
                                t.out = regGrid, d = 1)
    }
    estRaw.diag <- Lwls1D(bw.err, kernel_type = 'gauss', npoly = 1, nder = 0,
                          xin = obsGrid, yin = cov.diag$resiprod,
                          win = rep(1, nDiag), xout = obsGrid)
    
    t.idx.in <- which(obsGrid <= obs.range[2]*rational_cut[2] & obsGrid >= obs.range[2]*rational_cut[1])
    sigma2hat <- mean(estRaw.diag[t.idx.in] - estCov.diag[t.idx.in])
    sigma2hat <- ifelse(sigma2hat < 0, 1e-6, sigma2hat)
  }
  
  # perform eigen-analysis
  eig <- eigen(covhat.reg)
  positiveInd <- eig[['values']] >= 0
  stopifnot(sum(positiveInd) > 0)
  # padding zero to estimated eigenvalues if we estimate eigen-pairs too well.
  if(sum(positiveInd) < max_K){
    d <- c(eig[["values"]][positiveInd], rep(0, times = max_K - sum(positiveInd)))
  } else{
    d <- eig[["values"]][positiveInd]
  }
  # adjust the maximum number used in selection if necessary
  eigenV <- eig[["vectors"]][, 1:length(d), drop = FALSE]
  
  # normalization
  phi <- apply(eigenV, 2, function(x) {
    x <- x/sqrt(trapzRcpp(regGrid, x^2))
    if (0 <= sum(x*muhat.reg)) 
      return(x)
    else return(-x)
  })
  
  # obtain eigenvalues
  step.size <- diff(regGrid)[1]
  lambda.all <- step.size*d
  
  # re-fit covariance surface by discarding all negative eigenvalues
  covhat.reg <- phi[, 1:sum(positiveInd)] %*% diag(x = lambda.all, nrow = sum(positiveInd)) %*% t(phi[, 1:sum(positiveInd)])
  
  nlambda <- min(max_K, sum(positiveInd))
  # re-refine the max number of fPCs under consideration
  lambda <- lambda.all[1:max_K]
  phi <- phi[, 1:max_K, drop = FALSE]
  # phi must be a N by k (k can be 1) matrix
  if (is.vector(phi)) {
    phi = matrix(as.numeric(phi), nrow = length(phi), ncol = 1)
  }
  stopifnot(length(lambda) == ncol(phi))
  if(need.score){
    # perform score estimation - only needs max_K of scores
    if(estMethod == 'cross-sectional'){
      scoreEst <- mapply(function(yvec, tvec, mu, lambda, phi, nlambda){
        if(is.vector(phi)){
          phi = matrix(as.numeric(phi), nrow = length(phi), ncol = 1)
        }
        stopifnot(length(lambda) == ncol(phi))
        xiEst = matrix(0, length(lambda)) 
        cy = yvec - mu
        
        # Get Scores xiEst
        for(i in 1:length(lambda)){
          temp = cy*phi[, i]
          xiEst[i, 1] = trapzRcpp(X = tvec[!is.na(temp)], Y = temp[!is.na(temp)])
        }
        
        fittedY = mu + phi %*% xiEst
        res <- list('xiEst' = xiEst, 'xiVar' = NA, 'fittedY' = fittedY)
        return(res)
      }, fPCAinput$Ly, fPCAinput$Lt,
      MoreArgs = list(mu = muhat.reg, lambda = lambda, phi = phi, nlambda = nlambda))
    } else{
      scoreEst <- mapply(function(sub.idx, yVec, tVec, df, regGrid, lambda, phi, cov){
        muVec <- subset(df.sort, id == sub.idx, select = 'muhat', drop = TRUE)
        phiVec <- sapply(1:ncol(phi), function(c.) {
          spline(regGrid, phi[, c.], xout = tVec)$y
        })
        #phiVec <- t(phiVec)
        t2d <- expand.grid(t.x = tVec, t.y = tVec)
        covVec <- matrix(fdapace:::interp2lin(regGrid, regGrid, cov, t2d$t.x, t2d$t.y),
                         nrow = length(tVec))
        covVec <- 0.5*(covVec + t(covVec))
        sigmaY.vec <- covVec + diag(sigma2hat, length(tVec))
        res <- fdapace:::GetIndCEScoresCPP(yVec, muVec, lambda, phiVec, sigmaY.vec)
        return(res)
      }, fPCAinput$Lid, fPCAinput$Ly, fPCAinput$Lt,
      MoreArgs = list(df = df.sort, regGrid, lambda = lambda, phi = phi, cov = covhat.reg))
    }
    
    # Get fitted X
    scores <- t(do.call(cbind, scoreEst[1, ])) 
    fittedX <- t(do.call(cbind, scoreEst[3, ])) # n by N matrix
    
    return(list(sig2ma = sigma2hat, regGrid = regGrid, muhat = muhat.reg,
                covhat = covhat.reg, eigVals = lambda, eigVals.all = lambda.all,
                eigFuncs = phi, scores = scores, yhat = fittedX,
                bwmu = bw.mu, bwcov = bw.cov))
    
  } else{
    # no scores needed. No fitted y will be provided
    return(list(sig2ma = sigma2hat, regGrid = regGrid, muhat = muhat.reg,
                covhat = covhat.reg, eigVals = lambda, eigVals.all = lambda.all,
                eigFuncs = phi, bwmu = bw.mu, bwcov = bw.cov))
  }
}

fPCAselection <- function(df, n, numObs, fPCAinput, max_K = 5, estMethod = NULL, ...)
{
  
  # determine whether to use non-para smoothing in estimation
  idx <- unlist(fPCAinput$Lid)
  N.seq <- sapply(idx, function(n.) nrow(subset(df, id == n.)))
  obsGrid <- subset(df, select = 't', drop = TRUE)
  obs.range <- range(obsGrid)
  unique.t <- unique(obsGrid)
  if(is.null(estMethod)){
    cond1 <- length(unique(N.seq)) == 1
    cond2 <- length(unique.t)/diff(obs.range) >= 5
    if(cond1 & cond2 & nrow(df)/(length(unique.t)*length(N.seq)) == 1){
      estMethod <- 'cross-sectional'
    } else{
      estMethod <- 'smooth'
    }
  } else{
    estMethod <- match.arg(estMethod)
  }
  
  # use full data to obtain mu/eigenvalues/scores/yhat/sigma2hat
  # fix bandwidth for cov function estimation if neccessary to reduce computational time
  bw.cov.user = ifelse(numObs > 26, n^{-1/5}/6, NA)
  # bw.cov.user = NA
  #FPCA.mod <- fPCA.rev(df, fPCAinput, max_K = max_K, estMethod = estMethod, bw.cov = 0.5)
  FPCA.mod <- fPCA.rev(df, fPCAinput, max_K = max_K, estMethod = estMethod,
                       bw.cov = bw.cov.user)
  bwcov.opt <- FPCA.mod$bwcov
  # Get fitted values - all functions are evaluated on regGrid - full data
  muhat <- FPCA.mod$mu
  covhat <- FPCA.mod$covhat
  lambda <- FPCA.mod$eigVals
  lambda.all <- FPCA.mod$eigVals.all
  phi <- FPCA.mod$eigFuncs
  yhat <- FPCA.mod$yhat
  scores <- FPCA.mod$scores
  regGrid.full <- FPCA.mod$regGrid
  step.size <- diff(regGrid.full)[1]
  sig2mahat <- FPCA.mod$sig2ma
  sig2ma_LYH <- sum(diag(covhat) + sig2mahat)*step.size
  # data-splitting to obtain eigenfunctions on each subset
  block.size <- as.integer(n/2)
  stopifnot(all.equal(block.size, n/2))
  outIndex <- seq(1, block.size, by = 1)
  df.1 <- subset(df, id %in% outIndex)
  df.2 <- subset(df, !id %in% outIndex)
  
  fPCAinput.1 <- MakeFPCAInputs(IDs = df.1$id, tVec=df.1$t, yVec=df.1$y)
  fPCAinput.2 <- MakeFPCAInputs(IDs = df.2$id, tVec=df.2$t, yVec=df.2$y)
  
  # to reduce the computational cost, use the same bandwidth obtained from full data
  # if sparse, we could also re-calculate the bandwidth for covariance
  bw.cov.user.sub = ifelse(numObs > 26, bwcov.opt, NA)
  FPCA.mod.1 <- fPCA.rev(df.1, fPCAinput.1, max_K = max_K, bw.cov = bw.cov.user.sub,
                         estMethod = estMethod, need.score = FALSE)
  FPCA.mod.2 <- fPCA.rev(df.2, fPCAinput.2, max_K = max_K, bw.cov = bw.cov.user.sub,
                         estMethod = estMethod, need.score = FALSE)
  
  # if two regular grid are different
  # smooth two functions onto regGrid from full data
  B1 <- FPCA.mod.1$eigFuncs # eigenfunctions for 1st sub-sample
  B2 <- FPCA.mod.2$eigFuncs # eigenfunctions for 2nd sub-sample
  regGrid.1 <- FPCA.mod.1$regGrid
  regGrid.2 <- FPCA.mod.2$regGrid
  
  # re-define maximum number of fPC allowed in subsample
  max_K <- min(length(lambda), ncol(B1), ncol(B2))
  if(is.character(all.equal(regGrid.1, regGrid.2))){
    B1 <- sapply(1:max_K, function(c.) {
      spline(x=regGrid.1, y=B1[,c.], xout = regGrid.full, method = 'natural')$y
    })
    
    B2 <- sapply(1:max_K, function(c.) {
      spline(x=regGrid.2, y=B2[,c.], xout = regGrid.full, method = 'natural')$y
    })
  }
  
  # initializing f0 value
  f0Mat.cv <- matrix(rep(NA, max_K), ncol = max_K, nrow = 1)
  IC.df <- data.frame(matrix(NA, nrow = max_K, ncol = 4))
  names(IC.df) <- c('FY.AIC', 'FY.BIC', 'LYH.BIC', 'LYH.AIC')
  for(k. in 1:max_K){
    eigVals.K <- lambda[1:k.]
    eigFuncs.K <- phi[, 1:k., drop = FALSE]
    scores.K <- scores[, 1:k., drop = FALSE]
    theAngle <- abs(det(step.size*t(B1[, 1:k., drop = FALSE])%*%B2[, 1:k., drop = FALSE]))
    f0Mat.cv[k.] <- 1 - ifelse(theAngle > 1, 1, theAngle)
    
    fittedY.k <-  t(muhat + eigFuncs.K %*% t(scores.K)) # n by N matrix
    if(estMethod == 'smooth'){
      fittedY.list <- lapply(idx, function(n.) {
        tSub <- subset(df, id == n., select = 't', drop = TRUE)
        spline(regGrid.full, fittedY.k[n.,], xout = tSub)$y
      }) 
    }
    resi.k <- lapply(idx, function(n.) {
      if(estMethod == 'smooth'){
        fittedY.list[[n.]] - subset(df, id == n., select = 'y', drop = TRUE)
      } else{
        fittedY.k[n.,] - subset(df, id == n., select = 'y', drop = TRUE)
      }
    })
    
    rss <- sapply(1:n, function(n.){
      matrix(resi.k[[n.]], nrow = 1)%*%matrix(resi.k[[n.]], ncol = 1)
    })
    rss <- sum(rss)
    loglik <- rss/sig2mahat + sum(N.seq)*log(2*pi*sig2mahat)
    IC.df[k., 1] <- loglik + k. # FY AIC
    IC.df[k., 2] <- loglik + log(n)*k. # FY BIC
    # LYH BIC
    IC.df[k., 3] <- sig2ma_LYH - sum(lambda.all[1:k.])
    if(IC.df[k., 3] <= 0) {
      IC.df[k., 3] <- Inf}
    else{
      if(k.+1 <= length(lambda.all)){
        idx.out.lambda <- seq(k.+1, length(lambda.all), by = 1)
        IC.df[k., 3] <- log(sum(N.seq))*k.*
          sqrt(sum(lambda.all[idx.out.lambda]^2))/sig2mahat
      } else{
        # degenerate case
        IC.df[k., 3] <- Inf
      }
    }
    IC.df[k., 4] <- (sum(N.seq))*(log(rss) - log(sum(N.seq)) + 1) + 2*n*k. # LYH AIC
  }
  return(list(g0 = lambda, f0 = f0Mat.cv, IC = IC.df))
}

# one-simulation run
simuFunc <- function(n, numObs, eigenVals, err_sig2ma, maxK = 5,
                     FVE.threshold = 0.99, is.Gaussian = TRUE, is.regular = FALSE)
{
  # obtain data
  dat.list <- lapply(1:n, function(n.) {
    if(!is.regular){
      # use sample.int instead of runif to avoid duplicates with accuracy at 1e-10
      cur.t <- sort(sample.int(1e10, size = numObs, replace = F)/1e10)
    } else{
      cur.t <- seq(0, 1, length.out = numObs)
    }
    subjData <- datGP(obsGrid = cur.t, N = numObs, eigenVals = eigenVals,
                      errSig2ma = err_sig2ma, Gau = is.Gaussian)
    return(data.frame(id = n., subjData))
  })
  
  df <- do.call(rbind, dat.list)
  # make FPCA inputs for further analysis
  fPCAobj <- MakeFPCAInputs(IDs = df$id, tVec=df$t, yVec=df$y)
  kSelection <- fPCAselection(df = df, n = n, numObs = numObs,
                              fPCAinput = fPCAobj, max_K = maxK)
  FVE.lambda <- cumsum(kSelection$g0)/sum(kSelection$g0)
  k.FVE = match(TRUE, FVE.lambda >= FVE.threshold)
  k.prop <- sapply(1:length(kSelection$g0), function(c.){
    gk <- kSelection$g0[c.]/sum(kSelection$g0)
    if(abs(sum(kSelection$f0)) <= 1e-10){
      fk1 <- Inf
      fk3 <- Inf
    } else{
      fk1 <- kSelection$f0[c.]/sum(kSelection$f0)
      fk3 <- kSelection$f0[c.]/mean(kSelection$f0)
    }
    
    fk2 <- kSelection$f0[c.]/(1 + sum(kSelection$f0))
    fk4 <- kSelection$f0[c.]/(1 + mean(kSelection$f0))
    k.prop.wn1 <- gk + fk1
    k.prop.wn2 <- gk + fk2
    k.prop.wn3 <- gk + fk3
    k.prop.wn4 <- gk + fk4
    k.prop.won <- gk + kSelection$f0[c.]
    return(c(k.prop.wn1, k.prop.wn2, k.prop.wn3, k.prop.wn4, k.prop.won))
  })
  k.prop <- sapply(1:nrow(k.prop), function(r.) {
    ifelse(all(k.prop[r.,] == Inf), maxK, which.min(k.prop[r.,]))
  })
  names(k.prop) <- c(paste0('k.prop.wn', seq(1, (length(k.prop) - 1), by = 1)),
                     "k.prop.won")
  k.IC <- sapply(1:ncol(kSelection$IC), function(c.) {
    return(which.min(kSelection$IC[, c.]))
  })
  names(k.IC) <- names(kSelection$IC)
  final_res <- c(k.prop, 'FVE' = k.FVE, k.IC)
  return(final_res)
}

## working path
# wkPath <- ifelse(Sys.info()[['sysname']] == "Windows",
#                    "C:/Users/Chi Zhang/OneDrive - University of Waterloo/Desktop/Research/Data",
#                    "/Users/czhang/Library/CloudStorage/OneDrive-UniversityofWaterloo/Desktop/Research/Projects/order-determination"
# )
# 
# if(getwd() != wkPath) setwd(wkPath)
##' fixed parameters set-up - free to modify
eigenVals <- list(lam.1 = c(9, 4, 1), lam.2 = c(36, 25, 16, 9, 4, 1))
maxK <- c(6, 12) # maximum number of fPCs under consideration
sig2ma <- c(0.1, 0.5, 1, 4) # variance for random errors
n <- list(n.1 = c(100, 200), n.2 = c(200, 300)) # sample size
numObs <- c(11, 26, 51) # number of observations per curve
M <- 500 # simulation size
sparse.lv <- c('sparse', 'neither', 'dense')
scenarios <- lapply(1:length(eigenVals), function(i.) {
  expand.grid(n = n[[i.]], sparse.lv = sparse.lv, eigenvals = eigenVals[i.],
              sig2ma = sig2ma, isGau = c(TRUE, FALSE), isRegular = c(TRUE, FALSE),
              stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
})

scenarios <- do.call(rbind, scenarios)
# if the data is sparse/neither, we do not consider regular collected data
scenarios <- subset(scenarios, !(sparse.lv %in% c('sparse', 'neither') & isRegular == TRUE))
#scenarios$isRegular <- ifelse(scenarios$sparse.lv == 'dense', TRUE, FALSE)
scenarios$maxK <- ifelse(names(scenarios$eigenvals) == 'lam.1', maxK[1], maxK[2])
scenarios <- cbind(scenario = 1:nrow(scenarios), scenarios)
row.names(scenarios) <- NULL
# Complete simulation runs - parallel computing
library(parallel)
RNGkind("L'Ecuyer-CMRG")
set.seed(1)
scenarios.rep <- scenarios[rep(seq_len(nrow(scenarios)), each = M), ]
rownames(scenarios.rep) <- NULL
scenarios.rep$seed <- sample.int(1e8L, nrow(scenarios.rep))
num.cores <- detectCores()
num.workers <- min(150L, num.cores, nrow(scenarios.rep))
start_time <- format(Sys.time(), '%Y_%m_%d_%H%M')
message(paste('Program starts at', start_time))
fileName <- paste('simulationRes', start_time, sep = '_')
allRes <- lapply(1:nrow(scenarios), function(s.){
  theScenarios <- subset(scenarios.rep, scenario == s.)
  setting <- scenarios[s., ]
  message(paste0('Working on scenario ', s.))
  numObs <- ifelse(setting$sparse.lv == 'dense', numObs[3],
                   ifelse(setting$sparse.lv == 'sparse', numObs[1], numObs[2]))
  eigenvals <- unlist(setting$eigenvals)
  selected.K <- mclapply(1:nrow(theScenarios), function(case){
    set.seed(theScenarios[case, 'seed'])
    Kselection <-
      simuFunc(n = setting$n, numObs = numObs, eigenVals = eigenvals,
               err_sig2ma = setting$sig2ma, maxK = setting$maxK,
               is.Gaussian = setting$isGau, is.regular = setting$isRegular)
    return(Kselection)
  }, mc.cores = num.workers, mc.preschedule = FALSE) # some cases need much more time than others
  selected.K <- do.call(rbind, selected.K)
  accuracy.mat <- apply(selected.K, 2, function(c.) sum(c. == length(eigenvals)))
  accuracy.mat <- accuracy.mat/M
  count.lv <- seq(1, max(selected.K), by = 1)
  occurence_table <- apply(selected.K, 2, function(c.) {
    tabulate(factor(c., levels = count.lv), nbins = max(count.lv))
  })
  occurence_table <- data.frame(t(occurence_table))
  names(occurence_table) <- as.character(count.lv)
  final_res <- list(raw.df = data.frame(selected.K),
                    accuracy.tb = cbind(setting, t(accuracy.mat)),
                    occurence.tb = occurence_table)
  saveRDS(final_res, file = paste0(fileName, '_scenario_', s., '.rds'))
  end_time <- format(Sys.time(), '%Y_%m_%d_%H%M')
  message(paste0('Scenario ', s., ' is done at ', end_time))
  return(final_res)
})

saveRDS(allRes, file = paste0(fileName, '_summary.rds'))
if(paste0(fileName, '_summary.rds') %in% list.files()){
  file.remove(dir(pattern = paste0(fileName, '_scenario_\\d+.rds'), full.names = TRUE))
}
end_time <- format(Sys.time(), '%Y_%m_%d_%H%M')
message(paste('Program ended at', end_time))

