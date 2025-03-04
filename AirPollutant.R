library(dplyr)
library(tidyverse)
library(httr)
library(jsonlite)
library(fdapace)
library(nnet)
### --------------- Data Import and Load Functions -----------------------------
# functions
fPCA.rev <- function(df, obsGrid, need.score = TRUE)
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
  
  d <- eig[["values"]][positiveInd]
  eigenV <- eig[["vectors"]][, 1:length(d), drop = FALSE]
  
  # normalization
  phi <- apply(eigenV, 2, function(x) {
    x <- x/sqrt(trapzRcpp(obsGrid, x^2))
    if (0 <= sum(x*muhat.reg)) 
      return(x)
    else return(-x)
  })
  
  # obtain eigenvalues
  step.size <- diff(obsGrid)[1]
  lambda <- step.size*d
  # re-fit covariance surface by discarding all negative eigenvalues
  covhat.reg <- phi[, 1:sum(positiveInd)] %*% diag(x = lambda, nrow = sum(positiveInd)) %*% t(phi[, 1:sum(positiveInd)])
  # phi must be a N by k (k can be 1) matrix
  if (is.vector(phi)) {
    phi = matrix(as.numeric(phi), nrow = length(phi), ncol = 1)
  }
  stopifnot(length(lambda) == ncol(phi))
  
  if(need.score){
    y.list <- lapply(1:nrow(df), function(r.) df[r., ])
    scoreEst <- mapply(function(yvec, tvec, mu, lambda, phi){
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
    MoreArgs = list(tvec = obsGrid, mu = muhat.reg, lambda = lambda, phi = phi))
    scores <- t(do.call(cbind, scoreEst[1, ]))  # n by d matrix
    fittedX <- t(do.call(cbind, scoreEst[3, ])) # n by N matrix
    return(list(sig2ma = sigma2hat, muhat = muhat.reg, covhat = covhat.reg,
                eigVals = lambda, eigFuncs = phi, scores = scores, yhat = fittedX))
  } else{
    return(list(sig2ma = sigma2hat, muhat = muhat.reg, covhat = covhat.reg,
                eigVals = lambda, eigFuncs = phi))
  }
}


fPCAselection <- function(df, obsGrid, FVE.threshold = 0.99)
{
  n <- nrow(df)
  Nobs <- length(obsGrid)
  FPCA.mod <- fPCA.rev(df, obsGrid, need.score = TRUE)
  # Get fitted values - all functions are evaluated on regGrid - full data
  muhat <- FPCA.mod$mu
  covhat <- FPCA.mod$covhat
  lambda <- FPCA.mod$eigVals
  phi <- FPCA.mod$eigFuncs
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
  max_K <- min(ncol(B1), ncol(B2), ncol(phi))
  
  # initializing f0 value
  f0Mat.cv <- matrix(rep(NA, max_K), ncol = max_K, nrow = 1)
  IC.df <- data.frame(matrix(NA, nrow = max_K, ncol = 4))
  names(IC.df) <- c('FY.AIC', 'FY.BIC', 'LYH.BIC', 'LYH.AIC')
  for(k. in 1:max_K){
    eigVals.K <- lambda[1:k.]
    eigFuncs.K <- phi[, 1:k., drop = FALSE]
    scores.K <- scores[, 1:k., drop = FALSE]
    temp.innerproduct <- abs(det(step.size*t(B1[, 1:k., drop = FALSE])%*%B2[, 1:k., drop = FALSE]))
    f0Mat.cv[k.] <- 1 - ifelse(temp.innerproduct > 1, 1, temp.innerproduct)
    
    fittedY.k <-  t(muhat + eigFuncs.K %*% t(scores.K)) # n by N matrix
    
    resi.k <- lapply(1:n, function(n.) {
      fittedY.k[n.,] - df[n.,]
    })
    
    rss <- sapply(1:n, function(n.){
      sum(unlist(resi.k[[n.]])^2, na.rm = TRUE)
    })
    rss <- sum(rss, na.rm = TRUE)
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
  
  FVE.lambda <- cumsum(lambda)/sum(lambda)
  k.FVE = match(TRUE, FVE.lambda >= FVE.threshold)
  k.prop <- sapply(1:max_K, function(c.){
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
  return(list(muhat = muhat, lambda = lambda, eigenfuncs = phi,
              f0Mat = f0Mat.cv, IC = IC.df, k.selction = k.res))
}

mod.res <- function(df, obsGrid, workday.status, season.status,
                    split.ratio = 0.1){
  df <- subset(df, workday == workday.status & season == season.status)
  # train-test split
  n.train <- floor(nrow(df)*(1-split.ratio))
  n.test <- nrow(df) - n.train
  test.idx <- sample(nrow(df), size = n.test, replace = FALSE)
  df.train <- df[-test.idx,] 
  df.test <- df[test.idx, ]
  # modelling
  X.mat <- as.matrix(df.train[, paste(seq(0, 23, by = 1))])
  X.mat.test <- as.matrix(df.test[, paste(seq(0, 23, by = 1))])
  yVec <- unlist(df.train[, 'PM2.5lv'])
  yVec.test <- unlist(df.test[, 'PM2.5lv'])
  
  fPCA.res <- fPCAselection(X.mat, obsGrid = obsGrid, FVE.threshold = 0.99)
  muhat <- fPCA.res$muhat
  k.selection <- fPCA.res$k.selction
  eigenvalues <- fPCA.res$lambda
  eigen.basis <- fPCA.res$eigenfuncs # p by all eigenfuncs matrix 
  deMeanX.train <- lapply(1:nrow(X.mat), function(n.) X.mat[n., ] - muhat)
  deMeanX.test <- lapply(1:nrow(X.mat.test), function(n.) X.mat.test[n., ] - muhat)
  
  df.FPCA <- lapply(k.selection, function(k.){
    temp <- mapply(function(xVec, tVec, phi){
      res <- sapply(1:ncol(phi), function(c.){
        return(trapzRcpp(tVec, xVec*phi[, c.]))
      })
      return(res)
    }, deMeanX.train, MoreArgs = list(tVec = t.obs, phi = eigen.basis[, 1:k.]))
    res.df <- data.frame(y = yVec, X = t(temp))
    return(res.df)
  })
  
  df.FPCA.test <- lapply(k.selection, function(k.){
    temp <- mapply(function(xVec, tVec, phi){
      estScores <- sapply(1:ncol(phi), function(c.){
        return(trapzRcpp(tVec, xVec*phi[, c.]))
      })
      return(estScores)
    }, deMeanX.test, MoreArgs = list(tVec = obsGrid, phi = eigen.basis[, 1:k.]))
    # temp is d by n
    # fitted X is p by n
    fittedX <- sapply(1:ncol(temp), function(c.){
      muhat + eigen.basis[, 1:k.] %*% matrix(temp[,c.], ncol = 1) 
    })
    IMSE <- sapply(1:n.test, function(n.) {
      trapzRcpp(obsGrid, (X.mat.test[n., ] - fittedX[,n.])^2)
    })
    IMSE <- mean(IMSE, na.rm = TRUE)
    temp <-t(temp)
    colnames(temp) <- paste0('X.', seq(1, k., by = 1))
    res.df <- list(scores = temp, fittedX = t(fittedX), IMSE = IMSE)
    return(res.df)
  })
  
  fit.res <- lapply(1:length(df.FPCA), function(l.){
    fit <- multinom(y ~ ., data = df.FPCA[[l.]])
    pred.prob <- predict(fit, newdata = df.FPCA.test[[l.]]$scores, type = "class")
    accuracy <- sum(pred.prob == yVec.test)/length(yVec.test)
    return(list(mod = fit, accuracy = accuracy))
  })
  
  imseVec <- sapply(df.FPCA.test, function(l.){
    l.$IMSE
  })
  accuracyVec <- sapply(fit.res, function(l.) l.$accuracy)
  quanSummary <- rbind(k.selection, imseVec, accuracyVec)
  return(list(muhat = muhat, eigen.val = eigenvalues, eigenfunctions = eigen.basis,
              quanSummary = quanSummary))
}
##'------------------- Data import and preprocessing ---------------------------
is.remote <- FALSE # FAlSE when on local machine
dataPath <- ifelse(
  is.remote, "/u/c378zhan/projects/order_determination/data",
  ifelse(Sys.info()[['sysname']] == "Windows",
         "C:/Users/Chi Zhang/OneDrive - University of Waterloo/Desktop/Research/Data/Air pollution",
         "/Users/czhang/Library/CloudStorage/OneDrive-UniversityofWaterloo/Desktop/Research/Data"
         )
)
Huairoufile <- 'PRSA_Data_Huairou_20130301-20170228.csv'
HuairouData <- read.csv(paste0(dataPath, '/', Huairoufile), header = TRUE)
s.idx <- c(6, 7, 8) # summer month, June, July, Aug
w.idx <- c(12, 1, 2) # winter month, Dec, Jan, Feb
select.col <- c("PM2.5", "PM10", "SO2", "NO2", "CO", "O3", "TEMP", "PRES",
                "DEWP", "RAIN", "wd", "WSPM")
df <- subset(HuairouData, month %in% c(s.idx, w.idx),
             select = c('year', 'month', 'day', 'hour', select.col))
sum(is.na(df$PM2.5))/nrow(df) # 2.34% missing value
df$season <- ifelse(df$month %in% c(6, 7, 8), 1, 2) # 1 = summer; 2 = winter
df$date <- paste(df$year, df$month, df$day, sep = '-')
df$date <- as.Date(df$date)

# compensated days off, between 2013/03/01 and 2017/02/28
year.range <- range(df$year)
holiday.list <- lapply(seq(year.range[1], year.range[2], by = 1), function(y.){
  url.temp <- paste0('https://date.nager.at/api/v3/publicholidays/', y., '/CN')
  holidays.online <- GET(url.temp)
  holidays <- fromJSON(rawToChar(holidays.online$content))
  return(holidays)
})
holidays <- do.call(rbind, holiday.list)
holidays$date <- as.Date(holidays$date)
holidays <- subset(holidays, date <= '2017-02-28' & date >= '2013-03-01',
                   select=c('date', 'name'))


y2013.addon <- c('2013-04-04', '2013-04-06', '2013-04-29', '2013-04-30',
                 '2013-06-10', '2013-06-11', '2013-09-20', '2013-09-21',
                 paste('2013', '10', paste0('0', seq(2, 7)), sep = '-')
)

y2013.compensated <-
  c('2013-04-07', '2013-04-27', '2013-04-28', '2013-06-08',
    '2013-06-09', '2013-09-22', '2013-09-29', '2013-10-12')

y2014.addon <- c(paste('2014', '02', paste0('0', seq(1, 6)), sep = '-'),
                 '2014-04-07', '2014-05-02', '2014-05-03',
                 paste('2014', '10', paste0('0', seq(2, 7)), sep = '-'))

y2014.compensated <- 
  c('2014-01-26', '2014-02-08', '2014-05-04', '2014-09-28', '2014-10-11')

y2015.addon <- c('2015-01-02', '2015-01-03', '2015-02-18',
                 paste('2015', '02', paste0(seq(20, 24)), sep = '-'),
                 '2015-04-06', '2015-06-22',
                 paste('2015', '10', paste0('0', seq(2, 7)), sep = '-'))

y2015.compensated <- c('2015-01-04', '2015-02-15', '2015-02-28', '2015-10-10')


y2016.addon <- 
  c('2016-02-07', paste('2016', '02', paste0(seq(9, 13)), sep = '-'),
    '2016-04-04', '2016-05-02', '2016-06-10', '2016-06-11',
    '2016-09-16', '2016-09-17',
    paste('2016', '10', paste0('0', seq(2, 7)), sep = '-'))

y2016.compensated <- 
  c('2016-02-06', '2016-02-14', '2016-06-12', '2016-09-18',
    '2016-10-08', '2016-10-09')

y2017.addon <-
  c('2017-01-02', '2017-01-27', '2017-01-29', '2017-01-30', '2017-01-31',
    '2017-02-01', '2017-02-02')

y2017.compensated <- c('2017-01-22', '2017-02-04')


addon.holidays <- c(y2013.addon, y2014.addon, y2015.addon, y2016.addon, y2017.addon)
addon.workdays <- c(y2013.compensated, y2014.compensated, y2015.compensated,
                    y2016.compensated, y2017.compensated)

addon.workdays <- as.Date(addon.workdays)

holidays.all <- data.frame(date = sort(as.Date(c(holidays$date, addon.holidays))))
holidays.final <- left_join(holidays.all, holidays)
compensated.days <- data.frame(date = as.Date(addon.workdays))
# year 2016 Qingming jie is on Apr 4 instead of 5.
holidays.final[holidays.final$date == '2016-04-04', ]$name <-
  holidays.final[holidays.final$date == '2016-04-05', ]$name
holidays.final <- holidays.final[-which(holidays.final$date == '2016-04-05'), ]
compensated.days <- data.frame(date = sort(as.Date(c(addon.workdays))))

df.sub <- subset(df, month %in% c(s.idx, w.idx),
                 select = c('year', 'month', 'day', 'hour', 'DEWP', "season", "date"))
df.wider <- pivot_wider(data = df.sub, names_from = hour,
                        values_from = DEWP, values_fill = NA)
df.wider$which.day <- weekdays(df.wider$date)

df.wider$workday <- ifelse((df.wider$which.day %in% c('Sunday', 'Saturday') &
                              !df.wider$date %in% addon.workdays) |
                             df.wider$date %in% holidays.final$date, 0, 1) # 0 = weekend, 1 = workday
df.wider$week <- week(df.wider$date)
# multi-level
PM2.5DayAvg <- df %>% group_by(date) %>% 
  summarise(avgPM2.5 = mean(PM2.5, na.rm = TRUE)) %>% 
  mutate(level = case_when(
    avgPM2.5 <= 35 ~ 0,
    avgPM2.5 > 35 & avgPM2.5 <= 75 ~ 1,
    avgPM2.5 > 75 & avgPM2.5 <= 150 ~ 2,
    avgPM2.5 > 150 ~ 3))

df.wider$PM2.5lv <- PM2.5DayAvg$level
df.final <- df.wider[, c('date', 'season', 'workday', 'PM2.5lv', paste(seq(0, 23)))]
df.final <- subset(df.final, !is.na(PM2.5lv)) # rm days that NO PM2.5 records the whole day
NA.idx <- which(sapply(1:nrow(df.final), function(r.) sum(is.na(df.final[r.,])) > 0))
NA.df <- df.final[which(sapply(1:nrow(df.final), function(r.) sum(is.na(df.final[r.,])) > 0)), ]
df.final <- df.final[-NA.idx, ]


# workday has ~200 obs but weekend only has ~100
df.workday.winter <- subset(df.final, workday == 1 & season == 2)
# FPCA-selection
t.obs <- seq(0, 23, by = 1)/23
X.mat.winter <- as.matrix(subset(df.workday.winter, select = c(paste(seq(0, 23))), drop = TRUE), dimnames = NULL)
library(parallel)
RNGkind("L'Ecuyer-CMRG")
set.seed(1)
num.rep <- 500
seeds.case <- sample.int(1e7, num.rep, replace = FALSE)
num.cores <- detectCores()
num.workers <- min(150L, num.cores)
start_time <- format(Sys.time(), '%Y_%m_%d_%H%M')
message(paste('Program starts at', start_time))
fileName <- paste('realData_AirQuality', start_time, sep = '_')
allRes <- mclapply(1:num.rep, function(case){
  set.seed(seeds.case[case])
  res.winter <- mod.res(df.final, obsGrid = t.obs, workday.status = 1,
                        season.status = 2, split.ratio = 0.1)
  res.summer <- mod.res(df.final, obsGrid = t.obs, workday.status = 1,
                        season.status = 1, split.ratio = 0.1)
  theName <- rownames(res.winter$quanSummary)
  winter.df <- data.frame(rep = case, metric = theName, res.winter$quanSummary)
  summer.df <- data.frame(rep = case, metric = theName, res.summer$quanSummary)
  return(list(winter = winter.df, summer = summer.df))
}, mc.cores = num.workers, mc.preschedule = TRUE)

saveRDS(allRes, file = paste0(fileName, '.rds'))
end_time <- format(Sys.time(), '%Y_%m_%d_%H%M')
message(paste('Program ended at', end_time))
## ------------------- plot and summary ---------------------------------------
list.files(pattern = '.rds$')
allRes <- readRDS("realData_AirQuality_2025_02_11_1616.rds")

final_summary <- function(res.list, season){
  temp <- lapply(res.list, function(i.) i.[[season]])
  temp <- do.call(rbind, temp)
  theSummary <- temp %>% group_by(metric) %>% summarise(across(2:8, mean))
  row.idx <- temp$metric == 'k.selection'
  kTable <- temp[row.idx, 3:ncol(temp)]
  kVec <- sapply(1:ncol(kTable), function(c.) names(which.max(table(kTable[,c.]))))
  theSummary[3, 2:ncol(theSummary)] <- matrix(as.numeric(kVec), nrow = 1)
  return(theSummary)
}

finalRes.winter <- final_summary(allRes, 'winter')
finalRes.winter



