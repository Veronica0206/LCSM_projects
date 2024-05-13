## Main function to obtain matrix of covariates
### Input: Parameteres setting of covariates & logistic coefficients
### Output: Mean-vector and var-cov matrix of covariates
#####################################################################
getBetaX <- function(mean.x1 = 0, mean.x2 = 0, sd.x1 = 1, sd.x2 = 1, rho12 = 0.3, beta){
  #### Mean vector of covariates
  mean.x <- c(mean.x1, mean.x2)
  #### var-cov matrix of covariates
  var.x1 <- sd.x1^2; var.x2 <- sd.x2^2; x1x2 <- rho12 * sd.x1 * sd.x2
  phi <- matrix(c(var.x1, x1x2, x1x2, var.x2), nrow = 2, ncol = 2)
  true.x <- c(mean.x, phi[row(phi) >= col(phi)])
  ### Parameters setting of Beta
  true.B <- c(t(beta))
  return(list(mean.x, phi, true.x, true.B))
}

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

## Main function to generate data
### Input: Population values of parameters
### Output: Data (long-format and wide-format)
###############################################
getDat <- function(n, p, c = 2, mean.x, var.x, betaX, MeanYX, SigmaYX, X_Rate, Kappa, X_sd, Y_sd, XY_rho, delta = 0.25, seed = NA){
  ### Set seed if not an NA
  if (is.na(seed) == 0){
    set.seed(seed)
  }
  #### n is the sample size; p is the # of repated measures; c is the number of latent classes
  #### alpha & zeta are the mean vector and residuals of growth factors (when X are standardized)
  #### mean.x & var.x are for generating covariates
  #### Step 1: Generate component labels for each individual from corresponding covariates;
  #### Here, we suppose for each individual, we have two continuous covariates (as we set in parameters function)
  #### Generate covariates matrix for all individuals
  X <- MASS::mvrnorm(n = n, mu = mean.x, Sigma = var.x)
  #### Combine to the one vector
  oneX <- cbind(rep(1, n), X)
  #### Calculate exp(oneXbeta^T) as mixing proportions
  Prop <- exp(oneX %*% t(betaX))
  #### Calculate the probabilities for each individual to be labeled as each class
  Prob <- Prop/apply(Prop, 1, sum)
  mChoice <- t(apply(Prob, 1, rmultinom, n = 1, size = 1))
  dat0 <- data.frame(id = 1:n, gx1 = X[, 1], gx2 = X[, 2], 
                     z = apply(mChoice, 1, function(choose) {which(choose == 1)}))
  N <- rep(0, length(unique(dat0$z)))
  for (k in 1:length(unique(dat0$z))){
    N[k] <- sum(dat0$z == k)
  }
  ### Generate individual intercept, slope1, slope2, knot and covariates
  eta0.eta1.eta2.gamma.x <- as.data.frame(matrix(0, nrow = n, ncol = length(MeanYX[[1]])))
  colnames(eta0.eta1.eta2.gamma.x) <- c(paste0("muetaY", 0:2), "gamma", "ex", paste0("muetaX", 0:1))
  eta0.eta1.eta2.gamma.x$id <- 1:n
  KAPPA <- rep(0, n)
  SD_X <- SD_Y <- rho_XY <- rep(0, n)
  RateX <- matrix(0, nrow = n, ncol = p)
  for (k in 1:c){
    eta0.eta1.eta2.gamma.x[dat0$z == k, c(1:3, 5:7)] <- MASS::mvrnorm(N[k], ##### # of individuals in kth cluster
                                                                      mu = MeanYX[[k]][-4], 
                                                                      Sigma = SigmaYX[[k]])
    eta0.eta1.eta2.gamma.x[dat0$z == k, 4] <- MeanYX[[k]][4]
    KAPPA[dat0$z == k] <- rep(Kappa[k], N[k])
    SD_X[dat0$z == k] <- rep(X_sd[k], N[k])
    SD_Y[dat0$z == k] <- rep(Y_sd[k], N[k])
    rho_XY[dat0$z == k] <- rep(XY_rho[k], N[k])
    RateX[dat0$z == k, ] <- rep(X_Rate[[k]], each = N[k])
  }
  ### Long-format data framework
  long_dat <- expand.grid(id = 1:n, wave = 0:(p - 1))
  long_dat <- long_dat[order(long_dat$id, long_dat$wave), ]
  
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
        true_rate[(i - 1) * p + j, 1] <- eta0.eta1.eta2.gamma.x[df$id, 7] * RateX[df$id, j]
        true_chg[(i - 1) * p + j, 1] <- eta0.eta1.eta2.gamma.x[df$id, 7] * RateX[df$id, j] * lag
        true_chgBL[(i - 1) * p + j, 1] <- true_chgBL[(i - 1) * p + j - 1, 1] + eta0.eta1.eta2.gamma.x[df$id, 7] * RateX[df$id, j] * lag
        true_val[(i - 1) * p + j, 1] <- true_val[(i - 1) * p + j - 1, 1] + eta0.eta1.eta2.gamma.x[df$id, 7] * RateX[df$id, j] * lag
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
          KAPPA[df$id] * long_dat$CHG[long_dat$id == i][j] 
      }
      scoreXY[(i - 1) * p + j, ] <- scoreXY[(i - 1) * p + j, ] + 
        MASS::mvrnorm(1, mu = rep(0, 2), Sigma = matrix(c(SD_X[df$id]^2, rep(rho_XY[df$id] * SD_X[df$id] * SD_Y[df$id], 2), SD_Y[df$id]^2), nrow = 2, ncol = 2))
    }
  }
  long_dat$X <- scoreXY[, 1]
  long_dat$Y <- scoreXY[, 2]
  long_dat$class <- rep(dat0$z, each = p)
  
  ### Convert long format to wide format
  wide_Y <- reshape2::dcast(long_dat, id ~ wave, value.var = "Y")
  colnames(wide_Y) <- c("id", paste0("Y", 1:p))
  wide_X <- reshape2::dcast(long_dat, id ~ wave, value.var = "X")
  colnames(wide_X) <- c("id", paste0("TVC", 1:p))
  #wide_CHG <- reshape2::dcast(long_dat, id ~ wave, value.var = "obs_CHG")
  #colnames(wide_CHG) <- c("id", paste0("CHG", 1:p))
  wide_time <- reshape2::dcast(long_dat, id ~ wave, value.var = "time")
  colnames(wide_time) <- c("id", paste0("T", 1:p))
  names(dat0) <- c("id", paste0("gx", 1:2), "class")
  dat <- merge(merge(merge(wide_Y, wide_time, by = "id"), wide_X, by = "id"), dat0, by = "id")
  ### Add time-invariant covariate
  dat$TIC <- eta0.eta1.eta2.gamma.x[, 5]
  return(list(long_dat, dat))
}

getBLSGMM_TVCchg <- function(dat, p, c, init, beta, manifests, rateX, extratry = 10, loop = 20){
  ### Define latent variables
  latents <- c("eta0Ys", "eta1Ys", "eta2Ys", "eta0x", "eta1x", paste0("lx", 1:p), paste0("dx", 2:p), paste0("deltax", 2:p))
  outDef <- outLag <- outLoads1 <- outLoads2 <- list()
  class.list <- list()
  for (k in 1:c){
    for (j in 1:p){
      outDef[[j]] <- mxMatrix("Full", 1, 1, free = F, labels = paste0("data.T", j), name = paste0("t", j))
      outLag[[j]] <- mxAlgebraFromString(paste0("t", j , " -  t", j - 1), name = paste0("lag", j))
      outLoads1[[j]] <- mxAlgebraFromString(paste0("t", j, " -c", k, "mug"), name = paste0("c", k, "L1", j))
      outLoads2[[j]] <- mxAlgebraFromString(paste0("abs(t", j, " -c", k, "mug)"), name = paste0("c", k, "L2", j))
    }
    outLag[[1]] <- NULL
    abs_rate <- list()
    abs_rate[[2]] <- mxAlgebraFromString(paste0("c", k, "X_mueta1"), name = paste0("c", k, "abs_rate", 2))
    for (j in 3:p){
      abs_rate[[j]] <- mxAlgebraFromString(paste0("c", k, "X_mueta1 * c", k, "rel_rate", j), name = paste0("c", k, "abs_rate", j))
    }
    ### Create a mxModel object
    class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM", 
                               manifestVars = manifests, latentVars = latents,
                               mxData(observed = dat, type = "raw"),
                               #### Define factor loadings from latent variables to manifests related to Y
                               mxPath(from = "eta0Ys", to = manifests[(p + 1):(p * 2)], arrows = 1, free = F, values = 1),
                               mxPath(from = "eta1Ys", to = manifests[(p + 1):(p * 2)], arrows = 1, free = F, values = 0,
                                      labels = paste0("c", k, "L1", 1:p, "[1,1]")),
                               mxPath(from = "eta2Ys", to = manifests[(p + 1):(p * 2)], arrows = 1, free = F, values = 0, 
                                      labels = paste0("c", k, "L2", 1:p, "[1,1]")),
                               #### Define the variances of residuals
                               mxPath(from = manifests[(p + 1):(p * 2)], to = manifests[(p + 1):(p * 2)], arrows = 2, free = T,
                                      values = init[[k]][23], labels = paste0("c", k, "Y_residual")),
                               #### Define means of latent variables
                               mxPath(from = "one", to = latents[1:3], arrows = 1, free = T, values = init[[k]][1:3],
                                      labels = paste0("c", k, c("Y_mueta0s", "Y_mueta1s", "Y_mueta2s"))),
                               #### Define var-cov matrix of latent variables
                               mxPath(from = latents[1:3], to = latents[1:3], arrows = 2,
                                      connect = "unique.pairs", free = T,
                                      values = init[[k]][c(5:10)],
                                      labels = paste0("c", k, c("Y_psi00s", "Y_psi01s", "Y_psi02s", 
                                                                "Y_psi11s", "Y_psi12s", "Y_psi22s"))),
                               #### Add additional parameter and constraints
                               mxMatrix("Full", 1, 1, free = T, values = init[[k]][4], 
                                        labels = paste0("c", k, "muknot"), name = paste0("c", k, "mug")),
                               ############################################################################################
                               #### Define factor loadings from latent variables to manifests related to time-varying covariate
                               mxPath(from = "eta0x", to = "lx1", arrows = 1, free = F, values = 1),
                               mxPath(from = "eta1x", to = paste0("dx", 2:p), arrows = 1, free = c(F, rep(T, p - 2)), 
                                      values = rateX[[k]][-1], labels = paste0("c", k, "rel_rate", 2:p)),
                               #### Define means of TIC and latent variables of TVC
                               mxPath(from = "one", to = c("TIC", latents[4:5]), arrows = 1, free = T, values = c(0, init[[k]][11:12]),
                                      labels = paste0("c", k, c("muTIC", "X_mueta0", "X_mueta1"))),
                               #### Define var-cov matrix of TIC and latent variables of TVC
                               mxPath(from = c("TIC", latents[4:5]), to = c("TIC", latents[4:5]), arrows = 2,
                                      connect = "unique.pairs", free = c(rep(T, 2), F, rep(T, 3)),
                                      values = c(1, 0.3 * sqrt(init[[k]][13]), 0, init[[k]][13:15]),
                                      labels = c(paste0("c", k, c("phiTIC", "covBL")), NA, 
                                                 paste0("c", k, c("X_psi00", "X_psi01", "X_psi11")))),
                               mxPath(from = manifests[1:p], to = manifests[1:p], arrows = 2, free = TRUE,
                                      values = init[[k]][24], labels = paste0("c", k, "X_residual")),
                               mxPath(from = manifests[1:p], to = manifests[(p + 1):(p * 2)], arrows = 2, free = TRUE,
                                      values = init[[k]][25], labels = paste0("c", k, "Cov_XYres")),
                               #### Define latent true scores
                               mxPath(from = paste0("lx", 1:p), to = manifests[1:p], arrows = 1, free = F, values = 1),
                               #### Define path from latent instantaneous rate of change at each measurement to true scores
                               mxPath(from = paste0("dx", 2:p), to = paste0("deltax", 2:p), arrows = 1, free = F, values = 0,
                                      labels = paste0("lag", 2:p, "[1,1]")),
                               mxPath(from = paste0("deltax", 2:p), to = paste0("lx", 2:p), arrows = 1, free = F, values = 1),
                               #### Define autoregressive paths
                               mxPath(from = paste0("lx", 1:(p - 1)), to = paste0("lx", 2:p),
                                      arrows = 1, free = F, values = 1),
                               ############################################################################################
                               ##### Regression coefficients from TIC to growth factors of Y
                               mxPath(from = "TIC", to = latents[1:3], arrows = 1, free = T, 
                                      values = init[[k]][16:18], labels = paste0("c", k, "betaTIC", 0:2)),
                               ############################################################################################
                               #### Include time-varying covariate standardized baseline
                               mxPath(from = "lx1", to = latents[1:3], arrows = 1, free = T, 
                                      values = init[[k]][19:21], labels = paste0("c", k, "betaTVC", 0:2)),
                               mxPath(from = paste0("deltax", 2:p), to = manifests[(p + 2):(p * 2)], arrows = 1, free = TRUE, 
                                      values = init[[k]][22], labels = paste0("c", k, "kappa")),
                               ############################################################################################
                               #### Add mx objects defined out of mxModel
                               outDef, outLag, outLoads1, outLoads2, abs_rate,
                               ############################################################################################
                               #### Inverse transformation of mean vector and var-cov of growth factors of Y
                               ##### Define transformation function and matrix
                               mxAlgebraFromString(paste0("rbind(cbind(", "1,", "-c", k, "muknot,", "c", k, "muknot),",
                                                          "cbind(0, 1, -1), cbind(0, 1, 1))"),
                                                   name = paste0("c", k, "func")),
                               mxAlgebraFromString(paste0("rbind(cbind(", "1,", "-c", k, "muknot,", "c", k, "muknot),",
                                                          "cbind(0, 1, -1), cbind(0, 1, 1))"),
                                                   name = paste0("c", k, "grad")),
                               ##### Mean vector and var-cov matrix of growth factors of Y
                               mxAlgebraFromString(paste0("rbind(c", k, "Y_mueta0s, c", k, "Y_mueta1s, c", k, "Y_mueta2s)"), 
                                                   name = paste0("c", k, "Y_mean_s")),
                               mxAlgebraFromString(paste0("rbind(cbind(c", k, "Y_psi00s, c", k, "Y_psi01s, c", k, "Y_psi02s),",
                                                          "cbind(c", k, "Y_psi01s, c", k, "Y_psi11s, c", k, "Y_psi12s),",
                                                          "cbind(c", k, "Y_psi02s, c", k, "Y_psi12s, c", k, "Y_psi22s))"), 
                                                   name = paste0("c", k, "Y_psi_s")),
                               ##### Mean vector and var-cov matrix of growth factors of TVC
                               mxAlgebraFromString(paste0("rbind(c", k, "X_mueta0, c", k, "X_mueta1)"), name = paste0("c", k, "X_mean")),
                               mxAlgebraFromString(paste0("rbind(cbind(c", k, "X_psi00, c", k, "X_psi01),",
                                                          "cbind(c", k, "X_psi01, c", k, "X_psi11))"), name = paste0("c", k, "X_var")),
                               ##### Mean vector and var-cov matrix of growth factors of time_varying covariate
                               mxAlgebraFromString(paste0("rbind(c", k, "muTIC, c", k, "X_mueta0)"), name = paste0("c", k, "BL_mean")),
                               mxAlgebraFromString(paste0("rbind(cbind(c", k, "phiTIC, c", k, "covBL),",
                                                          "cbind(c", k, "covBL, c", k, "X_psi00))"), name = paste0("c", k, "BL_var")),
                               ##### Coefficients from baseline covariates to Y growth factors
                               mxAlgebraFromString(paste0("rbind(cbind(c", k, "betaTIC0, c", k, "betaTVC0),",
                                                          "cbind(c", k, "betaTIC1, c", k, "betaTVC1),",
                                                          "cbind(c", k, "betaTIC2, c", k, "betaTVC2))"),
                                                   name = paste0("c", k, "beta_s")),
                               mxAlgebraFromString(paste0("c", k, "func %*% (c", k, "Y_mean_s + c", k, "beta_s %*% c", k, "BL_mean)"), 
                                                   name = paste0("c", k, "Y_mean")),
                               
                               mxAlgebraFromString(paste0("c", k, "grad %*% c", k, "Y_psi_s %*% t(c", k, "grad)"), 
                                                   name = paste0("c", k, "Y_var_unexp")),
                               mxAlgebraFromString(paste0("c", k, "grad %*% c", k, "beta_s"), 
                                                   name = paste0("c", k, "beta")),
                               mxFitFunctionML(vector = T))
  }
  ### Make the class proportion matrix, fixing one parameter at a non-zero constant (one)
  classBeta <- mxMatrix(type = "Full", nrow = c, ncol = dim(beta)[2],
                        free = rep(c(F, rep(T, c - 1)), 3), values = beta,
                        labels = paste0("beta", rep(1:c), rep(0:2, each = c)), 
                        name = "classbeta")
  classPV <- mxMatrix(nrow = 3, ncol = 1, labels = c(NA, "data.gx1", "data.gx2"), 
                      values = 1, name = "weightsV")
  classP <- mxAlgebra(classbeta %*% weightsV, name = "weights")
  algebraObjective <- mxExpectationMixture(paste0("Class", 1:c), 
                                           weights = "weights", scale = "softmax")
  objective <- mxFitFunctionML()
  GMM_mx <- mxModel("BLSGMM with TVC interval specific change", 
                    mxData(observed = dat, type = "raw"), class.list, classBeta,
                    classPV, classP, algebraObjective, objective)
  model0 <- mxTryHard(GMM_mx, extraTries = extratry, 
                      initialGradientIterations = loop, OKstatuscodes = 0)
  model <- mxRun(model0, intervals = T)
  return(model)
}

## Main function to get point estimate (se) of parameters of interest
######################################################################
getEstimation <- function(model, paraNames, p, c){
  model.para <- summary(model)$parameters[, c(1, 5, 6)]
  model.est <- model.se <- est <- list()
  for (k in 1:c){
    abs_rate_est <- abs_rate_se <- rep(0, p)
    for (j in 3:p){
      abs_rate_est[j] <- mxEvalByName(paste0("c", k, "abs_rate", j), model@submodels[[k]])
      abs_rate_se[j] <- mxSE(paste0("Class", k, ".c", k, "abs_rate", j), model, forceName = T)
    }
    model.est[[k]] <- c(mxEvalByName(paste0("c", k, "Y_mean"), model@submodels[[k]]),
                        model.para[model.para$name == paste0("c", k, "muknot"), 2],
                        mxEvalByName(paste0("c", k, "Y_var_unexp"), model@submodels[[k]])[row(mxEvalByName(paste0("c", k, "Y_var_unexp"), model@submodels[[k]])) >= 
                                                                                            col(mxEvalByName(paste0("c", k, "Y_var_unexp"), model@submodels[[k]]))],
                        mxEvalByName(paste0("c", k, "X_mean"), model@submodels[[k]]),
                        mxEvalByName(paste0("c", k, "X_var"), model@submodels[[k]])[row(mxEvalByName(paste0("c", k, "X_var"), model@submodels[[k]])) >= 
                                                                                          col(mxEvalByName(paste0("c", k, "X_var"), model@submodels[[k]]))],
                        model.para[grep(paste0("c", k, "rel_rate"), model.para$name), 2], abs_rate_est[3:p],
                        model.para[model.para$name == paste0("c", k, "muTIC"), 2], model.para[model.para$name == paste0("c", k, "phiTIC"), 2], 
                        model.para[model.para$name == paste0("c", k, "covBL"), 2], model.para[model.para$name == paste0("c", k, "kappa"), 2],
                        mxEvalByName(paste0("c", k, "beta"), model@submodels[[k]])[, 1],
                        mxEvalByName(paste0("c", k, "beta"), model@submodels[[k]])[, 2], 
                        model.para[model.para$name == paste0("c", k, "Y_residual"), 2], model.para[model.para$name == paste0("c", k, "X_residual"), 2], 
                        model.para[model.para$name == paste0("c", k, "Cov_XYres"), 2])
    Y_mean_se <- mxSE(paste0("Class", k, ".c", k, "Y_mean"), model, forceName = T)
    Y_var_se <- mxSE(paste0("Class", k, ".c", k, "Y_var_unexp"), model, forceName = T)
    X_mean_se <- mxSE(paste0("Class", k, ".c", k, "X_mean"), model, forceName = T)
    X_var_se <- mxSE(paste0("Class", k, ".c", k, "X_var"), model, forceName = T)
    beta_se <- mxSE(paste0("Class", k, ".c", k, "beta"), model, forceName = T)
    model.se[[k]] <- c(Y_mean_se, model.para[model.para$name == paste0("c", k, "muknot"), 3], Y_var_se[row(Y_var_se) >= col(Y_var_se)],
                       X_mean_se, X_var_se[row(X_var_se) >= col(X_var_se)],
                       model.para[grep(paste0("c", k, "rel_rate"), model.para$name), 3], abs_rate_est[3:p],
                       model.para[model.para$name == paste0("c", k, "muTIC"), 3], model.para[model.para$name == paste0("c", k, "phiTIC"), 3], 
                       model.para[model.para$name == paste0("c", k, "covBL"), 3], model.para[model.para$name == paste0("c", k, "kappa"), 3],
                       beta_se[, 1], beta_se[, 2], model.para[model.para$name == paste0("c", k, "Y_residual"), 3], 
                       model.para[model.para$name == paste0("c", k, "X_residual"), 3], model.para[model.para$name == paste0("c", k, "Cov_XYres"), 3])
    est[[k]] <- data.frame(Name = paste0("c", k, paraNames), Estimate = model.est[[k]], SE = model.se[[k]])
  }
  est.beta <- data.frame(Name = paste0("beta", rep(2:c, 3), rep(0:2, each = c - 1)), 
                         Estimate = c(mxEval(classbeta, model)[-1, ]),
                         SE = c(mxSE(classbeta, model)[-1, ]))
  estimates <- rbind(do.call(rbind.data.frame, est), est.beta)
  return(estimates)
}

## Set the number of individuals and of repeated measures and seed
n <- 500; p <- 10; c <- 2

## Set "true" values to parameters
###################################
## Population values of covariates parameters
mean.x1 <- 0; mean.x2 <- 0; sd.x1 <- 1; sd.x2 <- 1; rho12 <- 0.3

## Population values of log coefficients
log_Beta <- matrix(c(0, 0, 0, 0, log(1.5), log(1.7)), byrow = T, nrow = 2)
BetaX <- getBetaX(beta = log_Beta)

### First latent class
#### Set true values of parameters for latent variables of time-varying covariate
##### eta0: the intercept (the measurement at initial status)
X_eta0.mean1 <- 0; X_eta0.var1 <- 1
##### eta1: the shape factor 
X_eta1.mean1 <- 5; X_eta1.var1 <- 1
##### rate loadings
X_rate1 <- c(0, seq(from = 1, by = -0.1, length.out = p - 1))
##### correlation between latent variables of X
X_rho1 <- 0.3
##### Residuals of X
X_sd1 <- 1

#### Set true values of parameters for latent variables of longitudinal outcome variable
##### eta0: the intercept (the measurement at initial status)
Y_eta0.mean1 <- 48; Y_eta0.var1 <- 25
##### eta1: the slope of the first stage
Y_eta1.mean1 <- 4.50; Y_eta1.var1 <- 1
##### eta2: the slope of the second stage
Y_eta.diff1 <- -2.85; Y_eta2.var1 <- 1
##### gamma: the transition time between two stages
Y_gamma.mean1 <- 5.0; 
##### correlation between latent variables of Y
Y_rho1 <- 0.3
##### Residuals of Y
Y_sd1 <- 1
##### Correlation between residuals X and Y
XY_rho1 <- 0.3

##### Coefficients from baseline of X to growth factors of Y
beta1 <- matrix(c(1.251505, 1.877258/sqrt(X_eta0.var1),
                  0.250301, 0.3754515/sqrt(X_eta0.var1),
                  0.250301, 0.3754515/sqrt(X_eta0.var1)), nrow = 3, byrow = T)

##### Coefficient from change-from-baseline of X to growth factors of Y
kappa1 <- 0.3

#### Second latent class
#### Set true values of parameters for latent variables of time-varying covariate
##### eta0: the intercept (the measurement at initial status)
X_eta0.mean2 <- 0; X_eta0.var2 <- 1
##### eta1: the shape factor 
X_eta1.mean2 <- 5; X_eta1.var2 <- 1
##### rate loadings
X_rate2 <- c(0, seq(from = 1, by = -0.1, length.out = p - 1))
##### correlation between latent variables of X
X_rho2 <- 0.3
##### Residuals of X
X_sd2 <- 1

#### Set true values of parameters for latent variables of longitudinal outcome variable
##### eta0: the intercept (the measurement at initial status)
Y_eta0.mean2 <- 52; Y_eta0.var2 <- 25
##### eta1: the slope of the first stage
Y_eta1.mean2 <- 5.00; Y_eta1.var2 <- 1
##### eta2: the slope of the second stage
Y_eta.diff2 <- -3.2; Y_eta2.var2 <- 1
##### gamma: the transition time between two stages
Y_gamma.mean2 <- 4.0; 
##### correlation between latent variables of Y
Y_rho2 <- 0.3
##### Residuals of Y
Y_sd2 <- 1
##### Correlation between residuals X and Y
XY_rho2 <- 0.3

##### Coefficients from baseline of X to growth factors of Y
beta2 <- matrix(c(1.251505, 1.877258/sqrt(X_eta0.var2),
                  0.250301, 0.3754515/sqrt(X_eta0.var2),
                  0.250301, 0.3754515/sqrt(X_eta0.var2)), nrow = 3, byrow = T)/sqrt(2)

##### Coefficient from change-from-baseline of X to growth factors of Y
kappa2 <- 0.6

parameters1 <- getPara(Y_eta0.mean = Y_eta0.mean1, Y_eta1.mean = Y_eta1.mean1, Y_eta.diff = Y_eta.diff1, Y_gamma.mean = Y_gamma.mean1, 
                       Y_eta0.var = Y_eta0.var1, Y_eta1.var = Y_eta1.var1, Y_eta2.var = Y_eta2.var1, Y_rho = Y_rho1, beta = beta1, 
                       X_eta0.mean = X_eta0.mean1, X_eta1.mean = X_eta1.mean1, X_eta0.var = X_eta0.var1, X_eta1.var = X_eta1.var1, X_rho = X_rho1)
parameters2 <- getPara(Y_eta0.mean = Y_eta0.mean2, Y_eta1.mean = Y_eta1.mean2, Y_eta.diff = Y_eta.diff2, Y_gamma.mean = Y_gamma.mean2, 
                       Y_eta0.var = Y_eta0.var2, Y_eta1.var = Y_eta1.var2, Y_eta2.var = Y_eta2.var2, Y_rho = Y_rho2, beta = beta2, 
                       X_eta0.mean = X_eta0.mean2, X_eta1.mean = X_eta1.mean2, X_eta0.var = X_eta0.var2, X_eta1.var = X_eta1.var2, X_rho = X_rho2)
parameters <- list(parameters1, parameters2)

X_Rate <- list(X_rate1, X_rate2)
Kappa <- c(kappa1, kappa2)
X_sd <- c(X_sd1, X_sd2)
Y_sd <- c(Y_sd1, Y_sd2)
XY_rho <- c(XY_rho1, XY_rho2)

mean.x <- BetaX[[1]]; var.x <- BetaX[[2]]
MeanYX <- SigmaYX <- true <- list(); initial <- list()
for (k in 1:c){
  MeanYX[[k]] <- parameters[[k]][[1]]
  SigmaYX[[k]] <- parameters[[k]][[2]]
  true[[k]] <- c(parameters[[k]][[3]][1:15], X_Rate[[k]][3:p], X_Rate[[k]][3:p] * parameters[[k]][[3]][12], 
                 parameters[[k]][[3]][16:18], Kappa[[k]], parameters[[k]][[3]][19:24], Y_sd[k]^2, X_sd[k]^2, XY_rho[k] * X_sd[k] * Y_sd[k])
  initial[[k]] <- c(parameters[[k]][[4]], Kappa[[k]], Y_sd[k]^2, X_sd[k]^2, XY_rho[k] * X_sd[k] * Y_sd[k])
}

paraNames <- c("Y_mueta0", "Y_mueta1", "Y_mueta2", "Y_mug", paste0("Y_psi", c("00", "01", "02", "11", "12", "22"), "_unexp"), 
               "X_mueta0", "X_mueta1", paste0("X_psi", c("00", "01", "11")), paste0("rel_rate", 3:p), paste0("abs_rate", 3:p),
               "muTIC", "varTIC", "covBL", "kappa", paste0("betaTIC", 0:2), paste0("betaTVC", 0:2), "Y_residual", "X_residual", "Cov_XYres")

dat <- getDat(n = n, p = p, c = c, mean.x = mean.x, var.x = var.x, betaX = log_Beta, 
              MeanYX = MeanYX, SigmaYX = SigmaYX, X_Rate = X_Rate, Kappa = Kappa, X_sd = X_sd, Y_sd = Y_sd, XY_rho = XY_rho)
### Fit the proposed model
###########################
tmpGMMTVCchg <- try(getBLSGMM_TVCchg(dat = dat[[2]], p = p, c = c, init = initial, beta = log_Beta,
                                     manifests = c(paste0("TVC", 1:p), paste0("Y", 1:p), "TIC"), 
                                     rateX = X_Rate, extratry = 10, loop = 20))

outTVC <- getEstimation(model = tmpGMMTVCchg, paraNames = paraNames, p = p, c = c)
outTVC$true <- c(true[[1]], true[[2]], log_Beta[2, ])
