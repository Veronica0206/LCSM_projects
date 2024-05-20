Decomposing Impact on Longitudinal Outcome of Time-varying Covariate
into Baseline Effect and Temporal Effect
================
Jin Liu

## OS, R version and OpenMx Version

``` r
library(nlpsem)
```

    ## Loading required package: OpenMx

    ## OpenMx may run faster if it is compiled to take advantage of multiple cores.

``` r
OpenMx::mxOption(model = NULL, key = "Default optimizer", "CSOLNP", reset = FALSE)
OpenMx::mxVersion()
```

    ## OpenMx version: 2.21.8 [GIT v2.21.8]
    ## R version: R version 4.2.2 (2022-10-31)
    ## Platform: aarch64-apple-darwin20 
    ## MacOS: 14.4.1
    ## Default optimizer: CSOLNP
    ## NPSOL-enabled?: No
    ## OpenMP-enabled?: No

## Specify parameters need to be print out

``` r
paraBLS_TVC.r <- c(
  "Y_alpha0", "Y_alpha1", "Y_alpha2", "Y_mug",
  paste0("Y_psi", c("00", "01", "02", "11", "12", "22")), "Y_residuals",
  "X_mueta0", "X_mueta1", paste0("X_psi", c("00", "01", "11")), 
  paste0("X_rel_rate", 2:9), paste0("X_abs_rate", 1:9), "X_residuals",
  paste0("betaTIC1", 0:2), paste0("betaTIC2", 0:2), paste0("betaTVC", 0:2),
  "muTIC1", "muTIC2", "phiTIC11", "phiTIC12", "phiTIC22",
  "Y_mueta0", "Y_mueta1", "Y_mueta2", 
  "covBL1", "covBL2", "kappa", "Cov_XYres")
```

## First Model: Latent Growth Curve Model with Decomposed TVC (Baseline and Interval-specific Slopes)

### Read in dataset for analyses (wide-format data)

``` r
load("TVC1_BLS_dat.RData")
```

### Fit the model

``` r
BLS_TVCslp_out <- getTVCmodel(dat = TVC1_BLS_dat, t_var = "T", records = 1:10, y_var = "Y", curveFun = "BLS", intrinsic = FALSE, y_model = "LGCM", TVC = "TVC",
                              decompose = 1, growth_TIC = c("ex1", "ex2"), res_scale = c(0.25, 0.25), res_cor = 0.3, tries = 10, paramOut = TRUE, names = "paraBLS_TVC.r")
```

### Visulize longitudinal outcomes

``` r
xstarts <- mean(TVC1_BLS_dat$T1)
Figure1 <- getFigure(
  model = BLS_TVCslp_out@mxOutput, y_var = "Y", curveFun = "BLS", sub_Model = "TVC", y_model = "LGCM", t_var = "T", records = 1:10, m_var = NULL, x_var = NULL,
  x_type = NULL, xstarts = xstarts, xlab = "Time", outcome = "Outcome"
)
```

    ## Treating first argument as an object that stores a character
    ## Treating first argument as an object that stores a character
    ## Treating first argument as an object that stores a character

``` r
show(Figure1)
```

    ## figOutput Object
    ## --------------------
    ## Trajectories: 1 
    ## Figure 1:

    ## `geom_smooth()` using method = 'gam' and formula = 'y ~ s(x, bs = "cs")'

![](OpenMx_demo9_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

## Second Model: Latent Growth Curve Model with Decomposed TVC (Baseline and Interval-specific Changes)

### Read in dataset for analyses (wide-format data)

``` r
load("TVC2_BLS_dat.RData")
```

### Fit the model

``` r
BLS_TVCchg_out <- getTVCmodel(dat = TVC2_BLS_dat, t_var = "T", records = 1:10, y_var = "Y", curveFun = "BLS", intrinsic = FALSE, y_model = "LGCM", TVC = "TVC",
                              decompose = 2, growth_TIC = c("ex1", "ex2"), res_scale = c(0.25, 0.25), res_cor = 0.3, tries = 10, paramOut = TRUE, names = "paraBLS_TVC.r")
```

### Visulize longitudinal outcomes

``` r
xstarts <- mean(TVC2_BLS_dat$T1)
Figure2 <- getFigure(
  model = BLS_TVCchg_out@mxOutput, y_var = "Y", curveFun = "BLS", sub_Model = "TVC", y_model = "LGCM", t_var = "T", records = 1:10, m_var = NULL, x_var = NULL,
  x_type = NULL, xstarts = xstarts, xlab = "Time", outcome = "Outcome"
)
```

    ## Treating first argument as an object that stores a character
    ## Treating first argument as an object that stores a character
    ## Treating first argument as an object that stores a character

``` r
show(Figure2)
```

    ## figOutput Object
    ## --------------------
    ## Trajectories: 1 
    ## Figure 1:

    ## `geom_smooth()` using method = 'gam' and formula = 'y ~ s(x, bs = "cs")'

![](OpenMx_demo9_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

## Third Model: Latent Growth Curve Model with Decomposed TVC (Baseline and Change from Baseline)

### Read in dataset for analyses (wide-format data)

``` r
load("TVC3_BLS_dat.RData")
```

### Fit the model

``` r
BLS_TVCchgBL_out <- getTVCmodel(dat = TVC2_BLS_dat, t_var = "T", records = 1:10, y_var = "Y", curveFun = "BLS", intrinsic = FALSE, y_model = "LGCM", TVC = "TVC",
                                decompose = 3, growth_TIC = c("ex1", "ex2"), res_scale = c(0.25, 0.25), res_cor = 0.3, tries = 10, paramOut = TRUE, names = "paraBLS_TVC.r")
```

### Visulize longitudinal outcomes

``` r
xstarts <- mean(TVC3_BLS_dat$T1)
Figure3 <- getFigure(
  model = BLS_TVCchgBL_out@mxOutput, y_var = "Y", curveFun = "BLS", sub_Model = "TVC", y_model = "LGCM", t_var = "T", records = 1:10, m_var = NULL, x_var = NULL,
  x_type = NULL, xstarts = xstarts, xlab = "Time", outcome = "Outcome"
)
```

    ## Treating first argument as an object that stores a character
    ## Treating first argument as an object that stores a character
    ## Treating first argument as an object that stores a character

``` r
show(Figure3)
```

    ## figOutput Object
    ## --------------------
    ## Trajectories: 1 
    ## Figure 1:

    ## `geom_smooth()` using method = 'gam' and formula = 'y ~ s(x, bs = "cs")'

![](OpenMx_demo9_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->
