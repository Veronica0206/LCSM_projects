## Main function to get parameters for data generation 
## and calculate true values for parameters would be estimated
###############################################################
getPara <- function(Y_eta0.mean, Y_eta1.mean, Y_eta.diff, Y_gamma.mean, 
                    Y_eta0.var, Y_eta1.var, Y_eta2.var, Y_rho, beta, 
                    X_eta0.mean, X_eta1.mean, X_eta0.var, X_eta1.var, X_rho){
  ### Mean vectors and var-cov matrix of growth factors of Y, TIC, and growth factors of TVC
  #### Mean vector
  Y_eta2.mean <- Y_eta1.mean + Y_eta.diff
  mean0 <- c(Y_eta0.mean, Y_eta1.mean, Y_eta2.mean, Y_gamma.mean, 0, X_eta0.mean, X_eta1.mean)
  #### Var-cov matrix
  ##### var-cov matrix of Y
  Y_eta0eta1 <- Y_rho * sqrt(Y_eta0.var * Y_eta1.var) 
  Y_eta0eta2 <- Y_rho * sqrt(Y_eta0.var * Y_eta2.var)
  Y_eta1eta2 <- Y_rho * sqrt(Y_eta1.var * Y_eta2.var)
  sigma.yy <- matrix(c(Y_eta0.var, Y_eta0eta1, Y_eta0eta2, 
                       Y_eta0eta1, Y_eta1.var, Y_eta1eta2, 
                       Y_eta0eta2, Y_eta1eta2, Y_eta2.var), nrow = 3, ncol = 3)
  ##### var-cov matrix of TIC and growth factors of TVC
  covBL <- 0.3 * sqrt(1.0 * X_eta0.var)
  X_eta0eta1 <- X_rho * sqrt(X_eta0.var * X_eta1.var) 
  phi <- matrix(c(1.0, covBL, 0, 
                  covBL, X_eta0.var, X_eta0eta1,
                  0, X_eta0eta1, X_eta1.var), nrow = 3, ncol = 3)
  ##### var-cov matrix of Y and covariates at baseline
  sigma.xy <-  beta %*% phi[1:2, 1:2]
  sigma0 <- rbind(cbind(sigma.yy, sigma.xy), cbind(t(sigma.xy), phi[1:2, 1:2]))
  ##### Unexplained var-cov matrix of growth factors of Y
  psi <- sigma.yy - beta %*% phi[1:2, 1:2] %*% t(beta) 
  
  func0 <- grad0 <- matrix(c(1, mean0[4], 0, 
                             0, 0.5, 0.5, 
                             0, -0.5, 0.5), nrow = 3, byrow = T)
  
  true <- c(mean0[1:4], psi[row(psi) >= col(psi)], 
            mean0[6:7], phi[2:3, 2:3][row(phi[2:3, 2:3]) >= col(phi[2:3, 2:3])], 
            0, 1.0, covBL, c(beta))
  mean0.s <- func0 %*% mean0[1:3]
  psi.s <- grad0 %*% psi %*% t(grad0)
  beta.s <- grad0 %*% beta
  true.s <- c(mean0.s[1:3], mean0[4], psi.s[row(psi.s) >= col(psi.s)], 
              mean0[6:7], phi[2:3, 2:3][row(phi[2:3, 2:3]) >= col(phi[2:3, 2:3])], c(beta.s))
  meanYX <- mean0
  sigmaYX <- rbind(cbind(sigma0, matrix(c(rep(0, 4), X_eta0eta1), ncol = 1)), 
                   cbind(matrix(c(rep(0, 4), X_eta0eta1), nrow = 1), phi[3, 3]))
  return(list(meanYX, sigmaYX, true, true.s))
}

## Main function to generate the data
######################################
getDat <- function(n, p, meanYX, sigmaYX, rateX, kappa, sdX, sdY, rhoXY, delta = 0.25, seed = NA){
  ### Set seed if not an NA
  if (is.na(seed) == 0){
    set.seed(seed)
  }
  
  ### Shell for long format data
  long_dat <- expand.grid(id = 1:n, wave = 0:(p - 1))
  long_dat <- long_dat[order(long_dat$id, long_dat$wave), ]
  
  ### Generate individual intercept, slope1, slope2, knot and covariates
  eta0.eta1.eta2.gamma.x <- matrix(0, nrow = n, ncol = length(meanYX))
  eta0.eta1.eta2.gamma.x[, -4] <- MASS::mvrnorm(n, mu = meanYX[-4], Sigma = sigmaYX)
  eta0.eta1.eta2.gamma.x[, 4] <- meanYX[4]
  
  ### Obtain individually varying time-points
  record <- rep(0, nrow(long_dat))
  for (i in 1:nrow(long_dat)){
    record[i] <- ifelse(long_dat$wave[i] != 0 & long_dat$wave[i] != p - 1,
                        runif(1, long_dat$wave[i] - delta, long_dat$wave[i] + delta), 
                        long_dat$wave[i])
  }
  long_dat$time <- record
  
  ### Calculate the time-varying covariate for each individual at each time-point
  true_val <- true_chgBL <- true_chg <- true_rate <- matrix(0, nrow = nrow(long_dat), ncol = 1)
  scoreXY <- matrix(0, nrow = nrow(long_dat), ncol = 2)
  for (i in 1:n){
    for (j in 1:p){
      df <- long_dat[long_dat$id == i, ][j, ] ### subset the j th observation of the ith individual
      if (j == 1){ # At t = 0, the measurement is the initial status, or the intercept
        true_rate[(i - 1) * p + j, 1] <- NA
        true_chg[(i - 1) * p + j, 1] <- 0
        true_chgBL[(i - 1) * p + j, 1] <- 0
        true_val[(i - 1) * p + j, 1] <- eta0.eta1.eta2.gamma.x[df$id, 6]
      }
      else{
        lag <- long_dat$time[long_dat$id == i][j] - long_dat$time[long_dat$id == i][j - 1]
        true_rate[(i - 1) * p + j, 1] <- eta0.eta1.eta2.gamma.x[df$id, 7] * rateX[j]
        true_chg[(i - 1) * p + j, 1] <- eta0.eta1.eta2.gamma.x[df$id, 7] * rateX[j] * lag
        true_chgBL[(i - 1) * p + j, 1] <- true_chgBL[(i - 1) * p + j - 1, 1] + eta0.eta1.eta2.gamma.x[df$id, 7] * rateX[j] * lag
        true_val[(i - 1) * p + j, 1] <- true_val[(i - 1) * p + j - 1, 1] + eta0.eta1.eta2.gamma.x[df$id, 7] * rateX[j] * lag
      }
    }
    scoreXY[, 1] <- true_val[, 1]
    long_dat$Rate <- true_rate[, 1]
    long_dat$CHG <- true_chg[, 1]
    long_dat$CHG_BL <- true_chgBL[, 1]
  }
  
  ### Add (individual) knot as a variable in long-format dataset
  long_dat$gamma <- rep(eta0.eta1.eta2.gamma.x[, 4], each = p)
  
  ### Calculate the loadings for slope1 and slope2 at each time-point
  long_dat$time1 <- ifelse(long_dat$time <= long_dat$gamma, long_dat$time, long_dat$gamma)
  long_dat$time2 <- ifelse(long_dat$time <= long_dat$gamma, 0, long_dat$time - long_dat$gamma)
  
  ### Calculate the repeated measurement for each individual at each time-point
  for (i in 1:n){
    for (j in 1:p){
      df <- long_dat[long_dat$id == i, ][j, ]
      if (j == 1){
        scoreXY[(i - 1) * p + j, 2] <- eta0.eta1.eta2.gamma.x[df$id, 1] + 
          df$time1 * eta0.eta1.eta2.gamma.x[df$id, 2] +
          df$time2 * eta0.eta1.eta2.gamma.x[df$id, 3] 
      }
      else{
        scoreXY[(i - 1) * p + j, 2] <- eta0.eta1.eta2.gamma.x[df$id, 1] + 
          df$time1 * eta0.eta1.eta2.gamma.x[df$id, 2] +
          df$time2 * eta0.eta1.eta2.gamma.x[df$id, 3] + 
          kappa * long_dat$Rate[long_dat$id == i][j] 
      }
      scoreXY[(i - 1) * p + j, ] <- scoreXY[(i - 1) * p + j, ] + 
        MASS::mvrnorm(1, mu = rep(0, 2), Sigma = matrix(c(sdX^2, rep(rhoXY * sdX * sdY, 2), sdY^2), nrow = 2, ncol = 2))
    }
  }
  long_dat$X <- scoreXY[, 1]
  long_dat$Y <- scoreXY[, 2]
  
  ### Convert long format to wide format
  wide_Y <- reshape2::dcast(long_dat, id ~ wave, value.var = "Y")
  colnames(wide_Y) <- c("id", paste0("Y", 1:p))
  wide_X <- reshape2::dcast(long_dat, id ~ wave, value.var = "X")
  colnames(wide_X) <- c("id", paste0("TVC", 1:p))
  #wide_CHG <- reshape2::dcast(long_dat, id ~ wave, value.var = "obs_CHG")
  #colnames(wide_CHG) <- c("id", paste0("CHG", 1:p))
  wide_time <- reshape2::dcast(long_dat, id ~ wave, value.var = "time")
  colnames(wide_time) <- c("id", paste0("T", 1:p))
  dat <- merge(merge(wide_Y, wide_time, by = "id"), wide_X, by = "id")
  ### Add time-invariant covariate
  dat$TIC <- eta0.eta1.eta2.gamma.x[, 5]
  return(list(long_dat, dat))
}

getBLSGM_TVCslp <- function(dat, p, init, manifests, rateX, extratry = extratry, loop = loop){
  ### Define latent variables
  latents <- c("eta0Ys", "eta1Ys", "eta2Ys", "eta0x", "eta1x", paste0("lx", 1:p), paste0("dx", 2:p))
  outDef <- outLag <- outLoads1 <- outLoads2 <- list()
  for (j in 1:p){
    outDef[[j]] <- mxMatrix("Full", 1, 1, free = F, labels = paste0("data.T", j), name = paste0("t", j))
    outLag[[j]] <- mxAlgebraFromString(paste0("t", j , " -  t", j - 1), name = paste0("lag", j))
    outLoads1[[j]] <- mxAlgebraFromString(paste0("t", j, " -  mug"), name = paste0("L1", j))
    outLoads2[[j]] <- mxAlgebraFromString(paste0("abs(t", j, " - mug)"), name = paste0("L2", j))
  }
  outLag[[1]] <- NULL
  abs_rate <- list()
  abs_rate[[2]] <- mxAlgebra(X_mueta1, name = paste0("abs_rate", 2))
  for (j in 3:p){
    abs_rate[[j]] <- mxAlgebraFromString(paste0("X_mueta1 * ", "rel_rate", j), name = paste0("abs_rate", j))
  }
  ### Create a mxModel object
  model_mx <- mxModel(name = "BLSGM with TVC slopes", type = "RAM", 
                      manifestVars = manifests, latentVars = latents,
                      mxData(observed = dat, type = "raw"),
                      #### Define factor loadings from latent variables to manifests related to Y
                      mxPath(from = "eta0Ys", to = manifests[(p + 1):(p * 2)], arrows = 1, free = F, values = 1),
                      mxPath(from = "eta1Ys", to = manifests[(p + 1):(p * 2)], arrows = 1, free = F, values = 0,
                             labels = paste0("L1", 1:p, "[1,1]")),
                      mxPath(from = "eta2Ys", to = manifests[(p + 1):(p * 2)], arrows = 1, free = F, values = 0, 
                             labels = paste0("L2", 1:p, "[1,1]")),
                      #### Define the variances of residuals
                      mxPath(from = manifests[(p + 1):(p * 2)], to = manifests[(p + 1):(p * 2)], arrows = 2, free = T,
                             values = init[23], labels = paste0("Y_residual")),
                      #### Define means of latent variables
                      mxPath(from = "one", to = latents[1:3], arrows = 1, free = T, values = init[1:3],
                             labels = paste0(c("Y_mueta0s", "Y_mueta1s", "Y_mueta2s"))),
                      #### Define var-cov matrix of latent variables
                      mxPath(from = latents[1:3], to = latents[1:3], arrows = 2,
                             connect = "unique.pairs", free = T,
                             values = init[c(5:10)],
                             labels = paste0(c("Y_psi00s", "Y_psi01s", "Y_psi02s", 
                                               "Y_psi11s", "Y_psi12s", "Y_psi22s"))),
                      #### Add additional parameter and constraints
                      mxMatrix("Full", 1, 1, free = T, values = init[4], 
                               labels = "muknot", name = "mug"),
                      ############################################################################################
                      #### Define factor loadings from latent variables to manifests related to time-varying covariate
                      mxPath(from = "eta0x", to = "lx1", arrows = 1, free = F, values = 1),
                      mxPath(from = "eta1x", to = paste0("dx", 2:p), arrows = 1, free = c(F, rep(T, p - 2)), 
                             values = rateX[-1], labels = paste0("rel_rate", 2:p)),
                      #### Define means of TIC and latent variables of TVC
                      mxPath(from = "one", to = c("TIC", latents[4:5]), arrows = 1, free = T, values = c(0, init[11:12]),
                             labels = c("muTIC", "X_mueta0", "X_mueta1")),
                      #### Define var-cov matrix of TIC and latent variables of TVC
                      mxPath(from = c("TIC", latents[4:5]), to = c("TIC", latents[4:5]), arrows = 2,
                             connect = "unique.pairs", free = c(rep(T, 2), F, rep(T, 3)),
                             values = c(1, 1.5, 0, init[13:15]),
                             labels = c("phiTIC", "covBL", NA, "X_psi00", "X_psi01", "X_psi11")),
                      mxPath(from = manifests[1:p], to = manifests[1:p], arrows = 2, free = TRUE,
                             values = init[24], labels = paste0("X_residual")),
                      mxPath(from = manifests[1:p], to = manifests[(p + 1):(p * 2)], arrows = 2, free = TRUE,
                             values = init[25], labels = paste0("Cov_XYres")),
                      #### Define latent true scores
                      mxPath(from = paste0("lx", 1:p), to = manifests[1:p], arrows = 1, free = F, values = 1),
                      #### Define path from latent instantaneous rate of change at each measurement to true scores
                      mxPath(from = paste0("dx", 2:p), to = paste0("lx", 2:p), arrows = 1, free = F, values = 0,
                             labels = paste0("lag", 2:p, "[1,1]")),
                      #### Define autoregressive paths
                      mxPath(from = paste0("lx", 1:(p - 1)), to = paste0("lx", 2:p),
                             arrows = 1, free = F, values = 1),
                      ############################################################################################
                      ##### Regression coefficients from TIC to growth factors of Y
                      mxPath(from = "TIC", to = latents[1:3], arrows = 1, free = T, 
                             values = init[16:18], labels = paste0("betaTIC", 0:2)),
                      ############################################################################################
                      #### Include time-varying covariate standardized baseline
                      mxPath(from = "lx1", to = latents[1:3], arrows = 1, free = T, 
                             values = init[19:21], labels = paste0("betaTVC", 0:2)),
                      mxPath(from = paste0("dx", 2:p), to = manifests[(p + 2):(p * 2)], arrows = 1, free = TRUE, 
                             values = init[22], labels = "kappa"),
                      ############################################################################################
                      #### Add mx objects defined out of mxModel
                      outDef, outLag, outLoads1, outLoads2, abs_rate,
                      ############################################################################################
                      #### Inverse transformation of mean vector and var-cov of growth factors of Y
                      ##### Define transformation function and matrix
                      mxAlgebra(rbind(cbind(1, -mug, mug),
                                      cbind(0, 1, -1),
                                      cbind(0, 1, 1)), name = "func"),
                      mxAlgebra(rbind(cbind(1, -mug, mug), 
                                      cbind(0, 1, -1),  
                                      cbind(0, 1, 1)), name = "grad"),
                      ##### Mean vector and var-cov matrix of growth factors of Y
                      mxAlgebra(rbind(Y_mueta0s, Y_mueta1s, Y_mueta2s), name = "Y_mean_s"),
                      mxAlgebra(rbind(cbind(Y_psi00s, Y_psi01s, Y_psi02s),
                                      cbind(Y_psi01s, Y_psi11s, Y_psi12s),
                                      cbind(Y_psi02s, Y_psi12s, Y_psi22s)), name = "Y_psi_s"),
                      ##### Mean vector and var-cov matrix of growth factors of TVC
                      mxAlgebra(rbind(X_mueta0, X_mueta1), name = "X_mean"),
                      mxAlgebra(rbind(cbind(X_psi00, X_psi01),
                                      cbind(X_psi01, X_psi11)), name = "X_var"),
                      ##### Mean vector and var-cov matrix of growth factors of time_varying covariate
                      mxAlgebra(rbind(muTIC, X_mueta0), name = "BL_mean"),
                      mxAlgebra(rbind(cbind(phiTIC, covBL),
                                      cbind(covBL, X_psi00)), name = "BL_var"),
                      ##### Coefficients from baseline covariates to  
                      mxAlgebra(rbind(cbind(betaTIC0, betaTVC0),
                                      cbind(betaTIC1, betaTVC1),
                                      cbind(betaTIC2, betaTVC2)), name = "beta_s"),
                      mxAlgebra(func %*% (Y_mean_s + beta_s %*% BL_mean), name = "Y_mean"),
                      mxAlgebra(grad %*% Y_psi_s %*% t(grad), name = "Y_var_unexp"),
                      mxAlgebra(grad %*% beta_s, name = "beta"))
  model0 <- mxTryHard(model_mx, extraTries = 9, 
                      initialGradientIterations = 20, OKstatuscodes = 0)
  model <- mxRun(model0, intervals = T)
  return(model)
}

## Main function to get point estimate (se) of parameters of interest
######################################################################
getEstimation <- function(model, paraNames, p){
  model.para <- summary(model)$parameters[, c(1, 5, 6)]
  abs_rate_est <- abs_rate_se <- rep(0, p)
  for (j in 3:p){
    abs_rate_est[j] <- mxEvalByName(paste0("abs_rate", j), model)
    abs_rate_se[j] <- mxSE(paste0("abs_rate", j), model, forceName = T)
  }
  model.est <- c(model$Y_mean$result, model.para[model.para$name == "muknot", 2],
                 model$Y_var_unexp$result[row(model$Y_var_unexp$result) >= col(model$Y_var_unexp$result)],
                 model$X_mean$result, model$X_var$result[row(model$X_var$result) >= col(model$X_var$result)],
                 model.para[grep("rel_rate", model.para$name), 2], abs_rate_est[3:p],
                 model.para[model.para$name == "muTIC", 2], model.para[model.para$name == "phiTIC", 2], 
                 model.para[model.para$name == "covBL", 2], model.para[model.para$name == "kappa", 2],
                 model$beta$result[, 1], model$beta$result[, 2], 
                 model.para[model.para$name == "Y_residual", 2], model.para[model.para$name == "X_residual", 2], 
                 model.para[model.para$name == "Cov_XYres", 2])
  Y_var_se <- mxSE(Y_var_unexp, model)
  X_var_se <- mxSE(X_var, model)
  model.se <- c(mxSE(Y_mean, model), model.para[model.para$name == "muknot", 3],
                Y_var_se[row(Y_var_se) >= col(Y_var_se)],
                mxSE(X_mean, model), X_var_se[row(X_var_se) >= col(X_var_se)],
                model.para[grep("rel_rate", model.para$name), 3], abs_rate_se[3:p],
                model.para[model.para$name == "muTIC", 3], model.para[model.para$name == "phiTIC", 3], 
                model.para[model.para$name == "covBL", 3], model.para[model.para$name == "kappa", 3],
                mxSE(beta, model)[, 1], mxSE(beta, model)[, 2],
                model.para[model.para$name == "Y_residual", 3], model.para[model.para$name == "X_residual", 3], 
                model.para[model.para$name == "Cov_XYres", 3])
  est <- data.frame(Name = paraNames, Estimate = model.est, SE = model.se)
  return(est)
}

## Set the number of individuals and of repeated measures and seed
n <- 200; p <- 10

## Set "true" values to parameters
###################################
### Set true values of parameters of one TIC and the baseline of TVC
# mean.x1 <- mean.x2 <- 0; sd.x1 <- sd.x2 <- 1; rho12 <- 0.3

### Set true values of parameters for latent variables of time-varying covariate
#### eta0: the intercept (the measurement at initial status)
X_eta0.mean <- 10; X_eta0.var <- 16
#### eta1: the shape factor 
X_eta1.mean <- 5; X_eta1.var <- 1
#### rate loadings
X_rate <- c(0, seq(from = 1, by = -0.1, length.out = p - 1))
#### correlation between latent variables of X
X_rho <- 0.3
#### Residuals of X
X_sd <- 1

### Set true values of parameters for latent variables of longitudinal outcome variable
#### eta0: the intercept (the measurement at initial status)
Y_eta0.mean <- 100; Y_eta0.var <- 25
#### eta1: the slope of the first stage
Y_eta1.mean <- 5; Y_eta1.var <- 1
#### eta2: the slope of the second stage
Y_eta.diff <- -3.2; Y_eta2.var <- 1
#### gamma: the transition time between two stages
Y_gamma.mean <- 4.5; 
#### correlation between latent variables of X
Y_rho <- 0.3
#### Residuals of Y
Y_sd <- 1
#### Correlation between residuals X and Y
rhoXY <- 0.3

### Coefficients from baseline of X to growth factors of Y
beta <- matrix(c(1.251505, 1.877258/sqrt(X_eta0.var),
                 0.250301, 0.3754515/sqrt(X_eta0.var),
                 0.250301, 0.3754515/sqrt(X_eta0.var)), nrow = 3, byrow = T)

### Coefficient from change-from-baseline of X to growth factors of Y
kappa <- 0.5

paraNames <- c("Y_mueta0", "Y_mueta1", "Y_mueta2", "Y_mug", paste0("Y_psi", c("00", "01", "02", "11", "12", "22"), "_unexp"), 
               "X_mueta0", "X_mueta1", paste0("X_psi", c("00", "01", "11")), paste0("rel_rate", 3:p), paste0("abs_rate", 3:p),
               "muTIC", "varTIC", "covBL", "kappa", paste0("betaTIC", 0:2), paste0("betaTVC", 0:2), "Y_residual", "X_residual", "Cov_XYres")

parameters <- getPara(Y_eta0.mean = Y_eta0.mean, Y_eta1.mean = Y_eta1.mean, Y_eta.diff = Y_eta.diff, Y_gamma.mean = Y_gamma.mean, 
                      Y_eta0.var = Y_eta0.var, Y_eta1.var = Y_eta1.var, Y_eta2.var = Y_eta2.var, Y_rho = Y_rho, beta = beta, 
                      X_eta0.mean = X_eta0.mean, X_eta1.mean = X_eta1.mean, X_eta0.var = X_eta0.var, X_eta1.var = X_eta1.var, X_rho = X_rho)
dat <- getDat(n = n, p = p, meanYX = parameters[[1]], sigmaYX = parameters[[2]], rateX = X_rate, kappa = kappa, sdX = X_sd, sdY = Y_sd, rhoXY = rhoXY)
init <- c(parameters[[4]], kappa, Y_sd^2, X_sd^2, rhoXY * X_sd * Y_sd)

### Fit the proposed model
###########################
tmpTVCslp <- try(getBLSGM_TVCslp(dat = dat[[2]], p = p, init = init, manifests = c(paste0("TVC", 1:p), paste0("Y", 1:p), "TIC"), 
                                 rateX = X_rate, extratry = extratry, loop = loop))

out <- getEstimation(model = tmpTVCslp, paraNames = paraNames, p = p)

out$true <- c(parameters[[3]][1:15], X_rate[3:p], X_rate[3:p] * parameters[[3]][12], 
              parameters[[3]][16:18], kappa, parameters[[3]][19:24], Y_sd^2, X_sd^2, rhoXY * X_sd * Y_sd)