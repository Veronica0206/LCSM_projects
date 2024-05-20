Jenss-Bayley Latent Change Score Model with Individual Ratio of Growth
Acceleration in the Framework of Individual Measurement Occasions
================
Jin Liu

## Require package would be used

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

## OS, R version and OpenMx Version

``` r
mxOption(model = NULL, key = "Default optimizer", "CSOLNP", reset = FALSE)
mxVersion()
```

    ## OpenMx version: 2.21.8 [GIT v2.21.8]
    ## R version: R version 4.2.2 (2022-10-31)
    ## Platform: aarch64-apple-darwin20 
    ## MacOS: 14.4.1
    ## Default optimizer: CSOLNP
    ## NPSOL-enabled?: No
    ## OpenMP-enabled?: No

## Read in dataset for analyses (wide-format data)

``` r
load("JB_dat.RData")
```

## First Model: Latent Change Score Model with Quadratic Functional Form w/ Random Ratio

### Specify parameters need to be print out

``` r
paraLCSM_JB.f <- c("mueta0", "mueta1", "mueta2", "ratio", paste0("psi", c("00", "01", "02", "0g", "11", "12", "1g", "22", "2g", "gg")), "residuals", 
                   paste0("abs_rate", 1:9), paste0("abs_rate_se", 1:9), paste0("change_in_interval", 1:9), paste0("change_in_interval_se", 1:9),
                   paste0("change_from_baseline", 1:9), paste0("change_in_interval_se", 1:9))
```

### Fit the model

``` r
JB.f_out <- getLCSM(dat = JB_dat, t_var = "T", y_var = "Y", curveFun = "JB", intrinsic = TRUE, records = 1:10, res_scale = 0.1, paramOut = TRUE, 
                    names = "paraLCSM_JB.f")
```

### Visulize longitudinal outcomes

``` r
xstarts <- mean(JB_dat$T1)
Figure1 <- getFigure(
  model = JB.f_out@mxOutput, y_var = "Y", curveFun = "JB", sub_Model = "LCSM", t_var = "T", records = 1:10, m_var = NULL, x_var = NULL,
  x_type = NULL, xstarts = xstarts, xlab = "Time", outcome = "Outcome"
)
```

    ## Treating first argument as an object that stores a character

``` r
show(Figure1)
```

    ## figOutput Object
    ## --------------------
    ## Trajectories: 1 
    ## Figure 1:

    ## `geom_smooth()` using method = 'gam' and formula = 'y ~ s(x, bs = "cs")'

![](OpenMx_demo4_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

    ## Figure 2:

![](OpenMx_demo4_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->

## Second Model: Latent Change Score Model with Jenss-Bayley Functional Form w/ Fixed Ratio

### Specify parameters need to be print out

``` r
paraLCSM_JB.r <- c("mueta0", "mueta1", "mueta2", "ratio", paste0("psi", c("00", "01", "02", "11", "12", "22")), "residuals", paste0("abs_rate", 1:9), 
                   paste0("abs_rate_se", 1:9), paste0("change_in_interval", 1:9), paste0("change_in_interval_se", 1:9), paste0("change_from_baseline", 1:9),
                   paste0("change_in_interval_se", 1:9))
```

### Fit the model

``` r
JB.r_out <- getLCSM(dat = JB_dat, t_var = "T", y_var = "Y", curveFun = "JB", intrinsic = FALSE, records = 1:10, res_scale = 0.1, paramOut = TRUE, 
                    names = "paraLCSM_JB.r")
```

### Visulize longitudinal outcomes

``` r
xstarts <- mean(JB_dat$T1)
Figure2 <- getFigure(
  model = JB.r_out@mxOutput, y_var = "Y", curveFun = "JB", sub_Model = "LCSM", t_var = "T", records = 1:10, m_var = NULL, x_var = NULL,
  x_type = NULL, xstarts = xstarts, xlab = "Time", outcome = "Outcome"
)
```

    ## Treating first argument as an object that stores a character

``` r
show(Figure2)
```

    ## figOutput Object
    ## --------------------
    ## Trajectories: 1 
    ## Figure 1:

    ## `geom_smooth()` using method = 'gam' and formula = 'y ~ s(x, bs = "cs")'

![](OpenMx_demo4_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

    ## Figure 2:

![](OpenMx_demo4_files/figure-gfm/unnamed-chunk-9-2.png)<!-- -->
