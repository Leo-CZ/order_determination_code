library(dplyr)
library(fdapace)
library(gam)
library(lubridate)
library(tidyr)
library(httr)
library(jsonlite)
# functions
## data pre-process
## Note: original data already considered day light saving on Mar 12/Nov 5, 2017.
data_preprocess <- function(df, startDate, pub_holidays, time.gap,
                            maxDuration.h = 5, 
                            membership.status = c('Casual', 'Member'))
{
  df <- subset(df, Member.type == membership.status)
  df$asDate <- as.Date(df$Start.date, tz = "EST")
  # if we want to consider how the total number of rentals affect
  # the total number of duration in that day, remove all these over-night records
  df$pass2day <- as.integer(as.Date(df$End.date, tz = "EST") - df$asDate)
  # table(df$pass2day)
  # sum(df$pass2day)/nrow(df)*100 # 880 of them, account for 0.745% of data
  df <- df[df$pass2day == 0, ] # all the same day records
  # remove all records that duration exceeds the max_duration (in hour)
  df <- df[df$Duration/3600 < maxDuration.h, ]
  
  # manipulate data
  df$day <- day(df$Start.date)
  df$start.hour <- hour(df$Start.date)
  df$min <- minute(df$Start.date)
  if(time.gap == 15){
    df <- df %>% mutate(hour.num = case_when(
      min < 15 ~ 1,
      min >= 15 & min < 30 ~ 2,
      min >= 30 & min < 45 ~ 3,
      min >= 45 & min < 60 ~ 4) + start.hour*4)
  } else if(time.gap == 20){
    df <- df %>% mutate(hour.num = case_when(
      min < 20 ~ 1,
      min >= 20 & min < 40 ~ 2,
      min >= 40 & min < 60 ~ 3) + start.hour*3)
  } else if(time.gap == 30){
    df <- df %>% mutate(hour.num = case_when(
      min < 30 ~ 1,
      min >= 30 & min < 60 ~ 2) + start.hour*2)
  } else {
    df <- df %>% mutate(hour.num = start.hour)
  }
  df$daysElapsed <- as.integer(df$asDate - as.Date(startDate, tz = 'EST'))
  # if the rental happens in holidays or weekends, assign false
  df$wday <- ifelse(
    wday(df$Start.date) %in% c(1, 7) | as.Date(df$Start.date) %in% pub_holidays$date,
    FALSE, TRUE)
  
  rental.counts <- df %>% group_by(wday, daysElapsed, hour.num) %>% summarise(n=n())
  rentalX.wday <- rental.counts[rental.counts$wday == TRUE,]
  rentalX.weekend <- rental.counts[rental.counts$wday == FALSE,]
  X.wday <- pivot_wider(rentalX.wday, id_cols = daysElapsed,
                        names_from = hour.num, values_from = n,
                        values_fill = 0, names_sort = TRUE)
  X.weekend <- pivot_wider(rentalX.weekend, id_cols = daysElapsed,
                           names_from = hour.num, values_from = n,
                           values_fill = 0, names_sort = TRUE)
  rental.time <- df %>%
    group_by(wday, daysElapsed) %>%
    summarise(rental_totalhours = sum(Duration)/3600)
  y.wday <- rental.time[rental.time$wday == TRUE, c('daysElapsed', 'rental_totalhours')]
  y.weekend <- rental.time[rental.time$wday == FALSE, c('daysElapsed', 'rental_totalhours')]
  dat.wday <- inner_join(y = y.wday, X.wday, by = join_by(daysElapsed))
  dat.weekend <- inner_join(y = y.weekend, X.weekend, by = join_by(daysElapsed))
  return(list(workday = dat.wday, weekend = dat.weekend))
}
# df is a n by L matrix
# each row represents discrete observations taken from the underlying random process.
fPCA.rev <- function(df, obsGrid, need.score = TRUE, FVE.threshold = 0.99)
{
  muhat.reg <- colMeans(df, na.rm = TRUE) # cross-sectional mean estimation
  rawCov <- cov(df, use = 'pairwise.complete.obs') # cross-sectional covariance estimation
  rawCov <- 0.5*(rawCov + t(rawCov))
  sigma2hat <- mean(diff(t(df), differences=2)^2, na.rm=TRUE)/choose(4, 2)
  sigma2hat <- ifelse(sigma2hat >= 0, sigma2hat, 0)
  covhat.reg <- rawCov
  diag(covhat.reg) <- diag(covhat.reg) - sigma2hat
  
  eig <- eigen(covhat.reg)
  positiveInd <- eig[['values']] >= 0
  stopifnot(sum(positiveInd) > 0)
  # padding zero if we estimate eigen-pairs too well.
  d <- eig[["values"]][positiveInd]
  # if(sum(positiveInd) < max_K){
  #   d <- c(eig[["values"]][positiveInd], rep(0, times = max_K - sum(positiveInd)))
  # } else{
  #   d <- eig[["values"]][positiveInd]
  # }
  # adjust the maximum number used in selection if necessary
  eigenV <- eig[["vectors"]][, 1:length(d), drop = FALSE]
  
  # normalization
  phi.all <- apply(eigenV, 2, function(x) {
    x <- x/sqrt(trapzRcpp(obsGrid, x^2))
    if (0 <= sum(x*muhat.reg)) 
      return(x)
    else return(-x)
  })
  
  # obtain eigenvalues
  step.size <- diff(obsGrid)[1]
  lambda.all <- step.size*d
  max_K <- match(TRUE, cumsum(lambda.all)/sum(lambda.all) >= 0.999)
  #lambda.all <- lambda.all[lambda.all >= 1e-10] # increase numerical stability.
  # in case we estimate eigenpairs too well. i.e. sum(positiveInd) == p_0
  nlambda <- min(max_K, sum(positiveInd))
  # re-fit covariance surface by discarding all (near-)zero and negative eigenvalues
  # re-refine the max number of fPCs under consideration
  lambda <- lambda.all[1:max_K]
  phi <- phi.all[, 1:max_K, drop = FALSE]
  covhat.reg <- phi[, 1:nlambda] %*% diag(x = lambda.all, nrow = nlambda) %*% t(phi[, 1:nlambda])
  # phi must be a N by k (k can be 1) matrix
  if (is.vector(phi)) {
    phi = matrix(as.numeric(phi), nrow = length(phi), ncol = 1)
  }
  stopifnot(length(lambda) == ncol(phi))
  
  if(need.score){
    y.list <- lapply(1:nrow(df), function(r.) df[r., ])
    # perform score estimation - only needs max_K of scores
    scoreEst <- mapply(function(yvec, tvec, mu, lambda, phi, nlambda){
      if(is.vector(phi)){
        phi = matrix(as.numeric(phi), nrow = length(phi), ncol = 1)
      }
      stopifnot(length(lambda) == ncol(phi))
      xiEst = matrix(0, length(lambda)) 
      cy = yvec - mu
      
      # Get Scores xiEst - integral method
      for(i in 1:length(lambda)){
        temp = cy*phi[, i]
        xiEst[i, 1] = trapzRcpp(X = tvec[!is.na(temp)], Y = temp[!is.na(temp)])
      }
      
      fittedY = mu + phi %*% xiEst
      res <- list('xiEst' = xiEst, 'xiVar' = NA, 'fittedY' = fittedY)
      return(res)
    },
    y.list,
    MoreArgs = list(tvec = obsGrid, mu = muhat.reg, lambda = lambda, phi = phi,
                    nlambda = nlambda))
    scores <- t(do.call(cbind, scoreEst[1, ]))  # n by d matrix
    fittedX <- t(do.call(cbind, scoreEst[3, ])) # n by N matrix
    return(list(sig2ma = sigma2hat, muhat = muhat.reg, covhat = covhat.reg,
                eigVals = lambda, eigVals.all = lambda.all, eigFuncs = phi,
                eigFuncs.all = phi.all, scores = scores, yhat = fittedX))
  } else{
    return(list(sig2ma = sigma2hat, muhat = muhat.reg, covhat = covhat.reg,
                eigVals = lambda, eigVals.all = lambda.all, eigFuncs = phi,
                eigFuncs.all = phi.all))
  }
}


fPCAselection <- function(df, obsGrid, FVE.threshold = 0.99)
{
  n <- nrow(df)
  Nobs <- length(obsGrid)
  FPCA.mod <- fPCA.rev(df, obsGrid, need.score = TRUE, FVE.threshold = FVE.threshold)
  # Get fitted values - all functions are evaluated on regGrid - full data
  muhat <- FPCA.mod$mu
  covhat <- FPCA.mod$covhat
  lambda <- FPCA.mod$eigVals
  lambda.all <- FPCA.mod$eigVals.all
  max_K <- match(TRUE, cumsum(lambda.all)/sum(lambda.all) >= 0.999)
  phi <- FPCA.mod$eigFuncs
  phi.all <- FPCA.mod$eigFuncs.all
  yhat <- FPCA.mod$yhat
  scores <- FPCA.mod$scores
  step.size <- diff(obsGrid)[1]
  sig2mahat <- FPCA.mod$sig2ma
  sig2ma_LYH <- sum(diag(covhat) + sig2mahat)*step.size
  # data-splitting to obtain eigenfunctions on each subset
  # 2-fold data-splitting
  split.idx <- sample.int(n, size = floor(n/2))
  df.1 <- df[split.idx, ] 
  df.2 <- df[-split.idx, ]
  
  # to reduce the computational cost, use the same bandwidth obtained from full data
  # if sparse, we could also re-calculate the bandwidth for covariance
  FPCA.mod.1 <- fPCA.rev(df.1, obsGrid, need.score = FALSE)
  FPCA.mod.2 <- fPCA.rev(df.2, obsGrid, need.score = FALSE)
  
  # if two regular grid are different
  # smooth two functions onto regGrid from full data
  B1 <- FPCA.mod.1$eigFuncs # eigenfunctions for 1st sub-sample
  B2 <- FPCA.mod.2$eigFuncs # eigenfunctions for 2nd sub-sample
  
  # re-define maximum number of fPC allowed in subsample
  if(min(length(lambda.all), ncol(B1), ncol(B2)) < max_K){
    max_K <- min(length(lambda), ncol(B1), ncol(B2))
    lambda.all <- lambda.all[1:max_K]
    phi.all <- phi.all[, 1:max_K]
    if(length(lambda) > max_K) {
      lambda <- lambda[1:max_K]
      phi <- phi[,1:max_K]
    }
    
  }
  
  # initializing f0 value
  f0Mat.cv <- matrix(rep(NA, max_K), ncol = max_K, nrow = 1)
  IC.df <- data.frame(matrix(NA, nrow = max_K, ncol = 4))
  names(IC.df) <- c('FY.AIC', 'FY.BIC', 'LYH.BIC', 'LYH.AIC')
  for(k. in 1:max_K){
    eigVals.K <- lambda[1:k.]
    eigFuncs.K <- phi[, 1:k., drop = FALSE]
    scores.K <- scores[, 1:k., drop = FALSE]
    f0Mat.cv[k.] <- 1 - abs(det(step.size*t(B1[, 1:k., drop = FALSE])%*%B2[, 1:k., drop = FALSE]))
    
    fittedY.k <-  t(muhat + eigFuncs.K %*% t(scores.K)) # n by N matrix
    
    resi.k <- lapply(1:n, function(n.) {
      fittedY.k[n.,] - df[n.,]
    })
    
    rss <- sapply(1:n, function(n.){
      matrix(unlist(resi.k[[n.]]), nrow = 1)%*%matrix(unlist(resi.k[[n.]]), ncol = 1)
    })
    rss <- sum(rss)
    loglik <- rss/sig2mahat + n*Nobs*log(2*pi*sig2mahat)
    IC.df[k., 1] <- loglik + k. # FY AIC
    IC.df[k., 2] <- loglik + log(n)*k. # FY BIC
    # LYH BIC
    IC.df[k., 3] <- sig2ma_LYH - sum(lambda[1:k.])
    if(IC.df[k., 3] <= 0) {
      IC.df[k., 3] <- Inf}
    else{
      if(k.+1 <= length(lambda)){
        idx.out.lambda <- seq(k.+1, length(lambda), by = 1)
        IC.df[k., 3] <- log(n*Nobs)*k.*
          sqrt(sum(lambda[idx.out.lambda]^2))/sig2mahat
      } else{
        # degenerate case
        IC.df[k., 3] <- Inf
      }
    }
    IC.df[k., 4] <- (n*Nobs)*(log(rss) - log(n*Nobs) + 1) + 2*n*k. # LYH AIC
  }
  
  FVE.lambda <- cumsum(lambda.all)/sum(lambda.all)
  k.FVE = match(TRUE, FVE.lambda >= FVE.threshold)
  k.prop <- sapply(1:length(lambda), function(c.){
    gk <- lambda[c.]/sum(lambda)
    fk <- f0Mat.cv[c.]/(1+sum(f0Mat.cv))
    k.prop.wn <- gk + fk
    k.prop.won <- gk + f0Mat.cv[c.]
    return(c(k.prop.wn, k.prop.won))
  })
  k.prop <- sapply(1:nrow(k.prop), function(r.) which.min(k.prop[r.,]))
  names(k.prop) <- c("k.cv.wn", "k.cv.won")
  k.IC <- sapply(1:ncol(IC.df), function(c.) {
    return(which.min(IC.df[, c.]))
  })
  names(k.IC) <- names(IC.df)
  k.res <- c(k.prop, 'FVE' = k.FVE, k.IC)
  return(list(muhat = muhat, lambda = lambda.all, eigenfuncs = phi,
              eigenfuncs.all = phi.all, f0Mat = f0Mat.cv, IC = IC.df,
              k.selction = k.res))
}
mod.sofr <- function(df, startDate, gap, split.ratio = 0.1,
                     dayType = c('workday', 'weekend'))
{
  df$daysElapsed <- NULL
  # train-test split
  n.train <- floor(nrow(df)*(1-split.ratio))
  n.test <- nrow(df) - n.train
  train.idx <- sample.int(nrow(df), size = n.train, replace = FALSE)
  df.train <- df[train.idx,] 
  df.test <- df[-train.idx, ]
  # modelling
  X.mat <- as.matrix(df.train[, -ncol(df)])
  X.mat.test <- as.matrix(df.test[, -ncol(df)])
  yVec <- unlist(df.train[, 'rental_totalhours'])
  yVec.test <- unlist(df.test[, 'rental_totalhours'])
  t.obs <- seq(1, ncol(X.mat), by = 1)*gap/60
  fPCA.res <- fPCAselection(X.mat, obsGrid = t.obs, FVE.threshold = 0.99)
  muhat <- fPCA.res$muhat
  k.selection <- fPCA.res$k.selction
  #eigen.basis <- fPCA.res$eigenfuncs # p by max_K matrix 
  eigen.basis <- fPCA.res$eigenfuncs.all # p by all eigenfuncs matrix 
  deMeanX.train <- lapply(1:nrow(X.mat), function(n.) X.mat[n., ] - muhat)
  deMeanX.test <- lapply(1:nrow(X.mat.test), function(n.) X.mat.test[n., ] - muhat)
  df.FPCA <- lapply(k.selection, function(k.){
    #if(k. == 1) browser()
    temp <- mapply(function(xVec, tVec, phi){
      res <- sapply(1:ncol(phi), function(c.){
        return(trapzRcpp(tVec, xVec*phi[, c.]))
      })
      return(res)
    }, deMeanX.train, MoreArgs = list(tVec = t.obs, phi = eigen.basis[, 1:k., drop = FALSE]))
    # temp is d by n matrix
    if(is.null(ncol(temp))) temp <- matrix(temp, nrow = 1)
    res.df <- data.frame(y = yVec, X = t(temp))
    return(res.df)
  })
  
  df.FPCA.test <- lapply(k.selection, function(k.){
    temp <- mapply(function(xVec, tVec, phi){
      estScores <- sapply(1:ncol(phi), function(c.){
        return(trapzRcpp(tVec, xVec*phi[, c.]))
      })
      return(estScores)
    }, deMeanX.test, MoreArgs = list(tVec = t.obs, phi = eigen.basis[, 1:k., drop = FALSE]))
    # temp is d by n
    # fitted X is p by n
    if(is.null(ncol(temp))) temp <- matrix(temp, nrow = 1)
    fittedX <- sapply(1:ncol(temp), function(c.){
      muhat + eigen.basis[, 1:k., drop = FALSE] %*% matrix(temp[,c.], ncol = 1) 
    })
    res.df <- list(data = data.frame(y = yVec.test, X = t(temp)),
                   fittedX = t(fittedX))
    return(res.df)
  })
  
  mod.sofr <- lapply(1:length(df.FPCA), function(l.){
    fit <- lm(y ~., data=df.FPCA[[l.]])
    beta.scores <- fit$coefficients[-1]
    beta <- rowSums(sapply(1:length(beta.scores), function(l.) beta.scores[l.]*eigen.basis[, l.]))
    pred <- predict(fit, newdata = df.FPCA.test[[l.]]$data)
    MSE <- mean((pred - yVec.test)^2)
    IMSE <- mean(sapply(1:n.test, function(n.) {
      trapzRcpp(t.obs, (X.mat.test[n., ] - df.FPCA.test[[l.]]$fittedX[n.,])^2)
    }))
    return(list(eigenfunctions = eigen.basis, mod.fit = fit, beta = beta,
                mse = MSE, imse = IMSE))
  })
  theNames <- names(df.FPCA)
  names(mod.sofr) <- theNames
  mseVec <- sapply(mod.sofr, function(l.) l.$mse)
  imseVec <- sapply(mod.sofr, function(l.) l.$imse)
  quanSummary <- rbind(k.selection, mseVec, imseVec)
  return(list(mod.fit = mod.sofr, quanSummary = quanSummary))
}
# determine data path
is.remote <- FALSE
dataPath <- ifelse(
  is.remote, "/u/c378zhan/projects/order_determination/data",
  ifelse(Sys.info()[['sysname']] == "Windows",
         "C:/Users/Chi Zhang/OneDrive - University of Waterloo/Desktop/Research/Data",
         "/Users/czhang/Library/CloudStorage/OneDrive-UniversityofWaterloo/Desktop/Research/Data"
  )
)

if(!is.remote){
  bikePath = paste0(dataPath, '/Bike sharing')
  setwd(bikePath)
} else{
  setwd(dataPath)
}
 
# import data
raw.df <- lapply(list.files(pattern = "\\d+Q\\d-[a-z]+\\-[a-z]+.csv"),
                 function(l.){
                   read.csv(file = l., header = TRUE)
                   })
raw.df <- do.call(rbind, raw.df)
startDate <- as.character('2017-01-01', tz = "EST") # the minimum date of the data
# obtain holidays online
holidays.online <- GET('https://date.nager.at/api/v3/publicholidays/2017/US')
holidays <- fromJSON(rawToChar(holidays.online$content))
# need manual adjustment
washington.DC.holidays <- subset(holidays, global == TRUE, select = c(localName, date))
addon.holidays <- data.frame(localName = c('Emancipation Day', 'Juneteenth'),
                             date = c('2017-04-17', '2017-06-19'))
washington.DC.holidays <- rbind(washington.DC.holidays, addon.holidays)
washington.DC.holidays$date <- as.Date(washington.DC.holidays$date)
washington.DC.holidays <- washington.DC.holidays[order(washington.DC.holidays$date), ]
row.names(washington.DC.holidays) <- NULL

gap <- 20
dayType <- 'weekend'
membership.status <- 'Member'
scenarios <- expand.grid(membership = membership.status, gap = gap, stringsAsFactors = FALSE)
scenarios <- cbind(scenario = seq(1, nrow(scenarios)), scenarios)
rownames(scenarios) <- NULL
##' parallel computing for train-test split
library(parallel)
M <- 500
RNGkind("L'Ecuyer-CMRG")
set.seed(1)
scenarios.rep <- scenarios[rep(seq_len(nrow(scenarios)), each = M), ]
rownames(scenarios.rep) <- NULL
scenarios.rep$seed <- sample.int(1e8L, nrow(scenarios.rep))
num.cores <- detectCores()
num.workers <- min(150L, num.cores, nrow(scenarios.rep))
start_time <- format(Sys.time(), '%Y_%m_%d_%H%M')
message(paste('Program starts at', start_time))
fileName <- paste('realData_BikeShare', start_time, sep = '_')

allRes <- lapply(1:nrow(scenarios), function(r.){
  theScenarios <- subset(scenarios.rep, scenario == r.)
  message(paste0('Working on scenario ', r.))
  scenario <- scenarios[r.,]
  df <- data_preprocess(raw.df, startDate = startDate, maxDuration.h = 5,
                        pub_holidays = washington.DC.holidays,
                        time.gap = scenario$gap, membership.status = scenario$membership)
  df.workday <- df[[dayType[1]]]
  df.weekend <- df[[dayType[2]]]
  res.list <- mclapply(1:nrow(theScenarios), function(case){
    set.seed(theScenarios[case, 'seed'])
    mod.workday <- mod.sofr(df = df.workday, startDate = startDate, split.ratio = 0.1,
                            dayType = dayType[1], gap = scenario$gap)
    mod.weekend <- mod.sofr(df = df.weekend, startDate = startDate, split.ratio = 0.1,
                            dayType = dayType[2], gap = scenario$gap)
    theName <- rownames(mod.workday$quanSummary)
    quan.workday <- data.frame(rep = case, time.gap = scenario$gap, 
                               metric = theName, mod.workday$quanSummary)
    quan.weekend <- data.frame(rep = case, time.gap = scenario$gap,
                               metric = theName, mod.weekend$quanSummary)
    return(list(workday = quan.workday, weekend = quan.weekend))
  }, mc.cores = num.workers, mc.preschedule = TRUE)
  temp.workday <- lapply(res.list, function(i.) i.[['workday']])
  temp.weekend <- lapply(res.list, function(i.) i.[['weekend']])
  workday.res <- do.call(rbind, temp.workday)
  weekend.res <- do.call(rbind, temp.weekend)
  return(list(workday = workday.res, weekend = weekend.res))
})

saveRDS(allRes, file = paste0(fileName, '.rds'))
end_time <- format(Sys.time(), '%Y_%m_%d_%H%M')
message(paste('Program ended at', end_time))

final_summary <- function(df, dayType){
  temp <- do.call(rbind, temp)
  temp$rep <- NULL
  theSummary <- temp %>%
    group_by(membership.status, metric, time.gap) %>%
    summarise(across(everything(), mean))

  kTable <- temp %>%
    filter(metric == 'k.selection') %>%
    select(!metric) %>%
    group_by(membership.status, time.gap) %>%
    reframe(across(everything(), ~ mlv1(.x, method='mfv')))

  row.idx <- which(theSummary$metric == 'k.selection')
  for (i. in 1:length(row.idx)) {
    newTable <- inner_join(theSummary[row.idx[i.], 1:3], kTable)
    theSummary[row.idx[i.], ] <- newTable
  }
  return(theSummary)
}

## read saved results
allRes <- readRDS('realData_BikeShare_2025_02_12_1309.rds')
library(modeest)
theSummary <- allRes %>%
  group_by(metric) %>%
  summarise(across(everything(), mean)) %>% 
  select(!c(rep, time.gap))
kTable <- allRes %>%
  filter(metric == 'k.selection') %>%
  select(!c(metric, rep)) %>%
  reframe(across(everything(), ~ mlv1(.x, method='mfv'))) %>% 
  select(!time.gap)
row.idx <- which(theSummary$metric == 'k.selection')
theSummary[row.idx, 2:8] <- kTable
theSummary

##' Illustration example ----------------------------------------
df <- data_preprocess(raw.df, startDate = startDate, maxDuration.h = 5,
                      pub_holidays = washington.DC.holidays,
                      time.gap = 20, membership.status = 'Member') # slow
#df.workday <- df[['workday']]
df.weekend <- df[['weekend']]
pdf('./plots/theNumOfRentals.pdf', width = 12, height = 7)
opar <- par(mar = c(4.2, 4.8, 1, 1))
matplot(0:23, t(df.weekend[, 2:25]), type = 'l', xaxt= "n", yaxt = 'n',
        xlab = substitute(paste(bold('hours'))),
        ylab = substitute(paste(bold('the number of rentals'))),
        cex.lab = 1.8)
axis(1, at = seq(0, 23, by = 2), cex.axis = 1.5)
axis(1, at = setdiff(seq(0, 23, 1), seq(0, 23, 2)), labels = NA, tck=-0.01)
axis(2, at = seq(0, 300, by = 100), cex.axis = 1.5)
par(opar)
dev.off()

pdf('./plots/totalRentalHours.pdf', width = 12, height = 8)
opar <- par(mar = c(4.5,4.5,2,0))
hist(df.weekend$rental_totalhours, main = '', xlab = 'total rental time (in hours)',
     cex.lab = 1.8, cex.axis = 1.5)
par(opar)
dev.off()

# fit the model by the entire set
df <- df.weekend
df$daysElapsed <- NULL
gap <- 20
# modelling
X.mat <- as.matrix(df[, -ncol(df)])
yVec <- unlist(df[, 'rental_totalhours'])
t.obs <- seq(1, ncol(X.mat), by = 1)*gap/60
fPCA.res <- fPCAselection(X.mat, obsGrid = t.obs, FVE.threshold = 0.99)
muhat <- fPCA.res$muhat
k.selection <- 4
#eigen.basis <- fPCA.res$eigenfuncs # p by max_K matrix
eigen.basis <- fPCA.res$eigenfuncs.all # p by all eigenfuncs matrix
deMeanX.train <- lapply(1:nrow(X.mat), function(n.) X.mat[n., ] - muhat)
df.FPCA <- lapply(k.selection, function(k.){
  #if(k. == 1) browser()
  temp <- mapply(function(xVec, tVec, phi){
    res <- sapply(1:ncol(phi), function(c.){
      return(trapzRcpp(tVec, xVec*phi[, c.]))
    })
    return(res)
  }, deMeanX.train, MoreArgs = list(tVec = t.obs, phi = eigen.basis[, 1:k., drop = FALSE]))
  # temp is d by n matrix
  if(is.null(ncol(temp))) temp <- matrix(temp, nrow = 1)
  res.df <- data.frame(y = yVec, X = t(temp))
  return(res.df)
})

mod.fit <- lapply(1:length(df.FPCA), function(l.){
  fit <- lm(y ~., data=df.FPCA[[l.]])
  beta.scores <- fit$coefficients[-1]
  beta <- rowSums(sapply(1:length(beta.scores), function(l.) beta.scores[l.]*eigen.basis[, l.]))
  return(list(eigenfunctions = eigen.basis, mod.fit = fit, beta = beta))
})

## plot - based on \hat{d} = 5
pdf('./plots/bikeSlope.pdf', width = 12, height = 7)
opar <- par(mar = c(4.4, 4.2, 2, 0.8))
plot(0:71, mod.fit[[1]]$beta, type = 'l', lwd = 2,
     xlab = substitute(paste(bold('hours'))), ylab = '', cex.lab = 1.5,
     xaxt = 'n', yaxt = 'n')
title(ylab=expression(beta(t)), line=2.0, cex.lab=1.5)
axis(1, at = seq(0, 71, by = 2), label = seq(0, 71, by = 2)/3, outer = FALSE, cex.axis = 1.5, tck=-0.02)
axis(1, at = setdiff(seq(0, 71, 1), seq(0, 23, 4)), labels = NA, tck=-0.01)
#axis(1, at=c(0,48), labels=c("",""), lwd.ticks=0)
axis(2, at = seq(0, 1.2, by = 0.2), cex.axis = 1.5)
par(opar)
dev.off()

## eigen-function plots
# first three eigenfunctions
library("colorspace")
pdf('./plots/bikeEigenFuncs.pdf', width = 12, height = 7)
opar <- par(mar = c(4.4, 4.2, 2, 0.8))
cols <- c('slateblue1', 'orangered4', 'black')
matplot(0:71, eigen.basis[, 1:3], type = 'l', lwd = 2,
        xlab = substitute(paste(bold('hours'))), ylab = '', cex.lab = 1.5,
        xaxt = 'n', yaxt = 'n', col = cols)
title(ylab=expression(phi(t)), line=2.0, cex.lab=1.5)
axis(1, at = seq(0, 71, by = 2), label = seq(0, 71, by = 2)/3, outer = FALSE, cex.axis = 1.5, tck=-0.02)
axis(1, at = setdiff(seq(0, 71, 1), seq(0, 23, 4)), labels = NA, tck=-0.01)
axis(2, at = seq(-0.5, 0.5, by = 0.2), cex.axis = 1.5)
par(opar)
dev.off()
# eigenfunction 6-7
pdf('./plots/bikeEigenFuncs_Extra.pdf', width = 12, height = 7)
opar <- par(mar = c(4.4, 4.2, 2, 0.8))
matplot(0:71, eigen.basis[, 5:6], type = 'l', lwd = 2,
        xlab = substitute(paste(bold('hours'))), ylab = '', cex.lab = 1.5,
        xaxt = 'n', yaxt = 'n', col = cols[1:2])
title(ylab=expression(phi(t)), line=2.0, cex.lab=1.5)
axis(1, at = seq(0, 71, by = 2), label = seq(0, 71, by = 2)/3, outer = FALSE, cex.axis = 1.5, tck=-0.02)
axis(1, at = setdiff(seq(0, 71, 1), seq(0, 23, 4)), labels = NA, tck=-0.01)
axis(2, at = seq(-0.5, 0.5, by = 0.2), cex.axis = 1.5)
par(opar)
dev.off()


