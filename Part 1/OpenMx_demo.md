Jenss-Bayley Latent Change Score Model with Individual Ratio of Growth
Acceleration in the Framework of Individual Measurement Occasions
================
Jin Liu

## Require package would be used

``` r
library(OpenMx)
```

    ## To take full advantage of multiple cores, use:
    ##   mxOption(key='Number of Threads', value=parallel::detectCores()) #now
    ##   Sys.setenv(OMP_NUM_THREADS=parallel::detectCores()) #before library(OpenMx)

``` r
library(tidyr)
library(ggplot2)
```

## OS, R version and OpenMx Version

``` r
mxOption(model = NULL, key = "Default optimizer", "CSOLNP", reset = FALSE)
mxVersion()
```

    ## OpenMx version: 2.19.8 [GIT v2.19.8-dirty]
    ## R version: R version 4.1.2 (2021-11-01)
    ## Platform: x86_64-apple-darwin17.0 
    ## MacOS: 12.3.1
    ## Default optimizer: CSOLNP
    ## NPSOL-enabled?: Yes
    ## OpenMP-enabled?: Yes

## “True” values of parameters

``` r
## Set "true" values to parameters
###################################
#### eta0: the intercept (the measurement at initial status)
# eta0.mean <- 50; eta0.var <- 16
#### eta1: the slope of the linear asymptote
# eta1.mean <- 2.5; eta1.var <- 1
#### eta2: the vertical distance between the actual intercept and the intercept of linear asymptote
# eta2.mean <- -30; eta2.var <- 36
#### exp(gamma): the ratio of acceleration of (t+1)/t 
# gamma.mean <- -0.7; gamma.var <- 0.10;
# rho <- 0.3
```

## Define Parameter lists

``` r
paraJB_LCSM_f <- c("mueta0", "mueta1", "mueta2", "gamma", 
                   paste0("psi", c("00", "01", "02", "11", "12", "22")),
                   paste0("instant_rate_est", 1:9), paste0("instant_rate_var", 1:9), paste0("change_in_interval", 1:9),
                   paste0("change_from_baseline", 1:9), "residuals")

paraJB_LCSM_r <- c("mueta0", "mueta1", "mueta2", "gamma", 
                   paste0("psi", c("00", "01", "02", "0g", "11", "12", "1g", "22", "2g", "gg")),
                   paste0("instant_rate_est", 1:9), paste0("instant_rate_var", 1:9), paste0("change_in_interval", 1:9),
                   paste0("change_from_baseline", 1:9), "residuals")
```

## Read in dataset for analyses (wide-format data)

``` r
load("JB_random_dat.RData")
```

## Summarize data

``` r
summary(JB_random_dat)
```

    ##        id               Y1              Y2              Y3       
    ##  Min.   :  1.00   Min.   :37.43   Min.   :54.10   Min.   :60.96  
    ##  1st Qu.: 50.75   1st Qu.:47.61   1st Qu.:63.52   1st Qu.:73.76  
    ##  Median :100.50   Median :49.77   Median :67.09   Median :77.48  
    ##  Mean   :100.50   Mean   :49.98   Mean   :67.23   Mean   :77.27  
    ##  3rd Qu.:150.25   3rd Qu.:52.45   3rd Qu.:70.54   3rd Qu.:80.20  
    ##  Max.   :200.00   Max.   :59.78   Max.   :81.33   Max.   :94.11  
    ##        Y4               Y5               Y6               Y7        
    ##  Min.   : 65.49   Min.   : 69.20   Min.   : 72.02   Min.   : 75.14  
    ##  1st Qu.: 78.97   1st Qu.: 82.90   1st Qu.: 86.04   1st Qu.: 88.33  
    ##  Median : 83.48   Median : 87.26   Median : 90.85   Median : 93.97  
    ##  Mean   : 83.37   Mean   : 87.64   Mean   : 90.99   Mean   : 93.79  
    ##  3rd Qu.: 87.05   3rd Qu.: 92.03   3rd Qu.: 95.58   3rd Qu.: 99.08  
    ##  Max.   :100.05   Max.   :106.39   Max.   :110.67   Max.   :114.85  
    ##        Y8               Y9              Y10               T1   
    ##  Min.   : 76.40   Min.   : 75.84   Min.   : 76.92   Min.   :0  
    ##  1st Qu.: 91.18   1st Qu.: 93.00   1st Qu.: 95.13   1st Qu.:0  
    ##  Median : 96.90   Median : 99.15   Median :101.23   Median :0  
    ##  Mean   : 96.67   Mean   : 99.13   Mean   :101.57   Mean   :0  
    ##  3rd Qu.:102.24   3rd Qu.:104.33   3rd Qu.:107.84   3rd Qu.:0  
    ##  Max.   :119.29   Max.   :123.88   Max.   :128.53   Max.   :0  
    ##        T2               T3              T4              T5       
    ##  Min.   :0.7503   Min.   :1.751   Min.   :2.753   Min.   :3.751  
    ##  1st Qu.:0.8653   1st Qu.:1.890   1st Qu.:2.853   1st Qu.:3.890  
    ##  Median :1.0085   Median :2.007   Median :3.011   Median :3.991  
    ##  Mean   :0.9945   Mean   :2.005   Mean   :2.998   Mean   :4.002  
    ##  3rd Qu.:1.1166   3rd Qu.:2.133   3rd Qu.:3.126   3rd Qu.:4.109  
    ##  Max.   :1.2488   Max.   :2.249   Max.   :3.246   Max.   :4.248  
    ##        T6              T7              T8              T9             T10   
    ##  Min.   :4.753   Min.   :5.751   Min.   :6.750   Min.   :7.756   Min.   :9  
    ##  1st Qu.:4.890   1st Qu.:5.897   1st Qu.:6.905   1st Qu.:7.892   1st Qu.:9  
    ##  Median :5.013   Median :6.006   Median :7.020   Median :7.992   Median :9  
    ##  Mean   :5.012   Mean   :6.010   Mean   :7.013   Mean   :8.000   Mean   :9  
    ##  3rd Qu.:5.134   3rd Qu.:6.128   3rd Qu.:7.147   3rd Qu.:8.124   3rd Qu.:9  
    ##  Max.   :5.249   Max.   :6.247   Max.   :7.247   Max.   :8.247   Max.   :9

## Visualize data

``` r
long_dat_T <- gather(JB_random_dat, var.T, time, T1:T10)
long_dat_Y <- gather(JB_random_dat, var.Y, measures, Y1:Y10)
long_dat <- data.frame(id = long_dat_T[, 1], time = long_dat_T[, 13],
                       measures = long_dat_Y[, 13])
ggplot(aes(x = time, y = measures), data = long_dat) +
  geom_line(aes(group = id), color = "lightgrey") +
  geom_point(aes(group = id), color = "darkgrey", size = 0.5) +
  geom_smooth(aes(group = 1), size = 1.8, col = "lightblue", se = F) + 
  labs(title = "Nonlinear Pattern with Individually Varying Measurement Time",
       x ="Time", y = "Measurement") + 
  theme(plot.title = element_text(hjust = 0.5))
```

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

![](OpenMx_demo_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

## Jenss-Bayley Latent Change Model with Random Ratio of Growth Acceleration

``` r
source("LCSM_JB_random.R")
JB_LCSM_r <- getLCSM_JB_random(dat = JB_random_dat, T_records = 1:10, traj_var = "Y", t_var = "T",
                               paraNames = paraJB_LCSM_r)
JB_LCSM_r[[2]]
```

    ##                     Name Estimate     SE
    ## 1                 mueta0  49.9774 0.2705
    ## 2                 mueta1   2.4459 0.0732
    ## 3                 mueta2 -30.3284 0.4669
    ## 4                  gamma  -0.7018 0.0080
    ## 5                  psi00  13.6655 1.4637
    ## 6                  psi01   1.3354 0.2944
    ## 7                  psi02   5.2062 1.8296
    ## 8                  psi0g  -0.0113 0.0309
    ## 9                  psi11   0.9467 0.1071
    ## 10                 psi12   1.9853 0.5218
    ## 11                 psi1g  -0.0040 0.0088
    ## 12                 psi22  36.7981 4.3597
    ## 13                 psi2g   0.0133 0.0561
    ## 14                 psigg   0.0000 0.0014
    ## 15     instant_rate_est1  17.4604 0.2200
    ## 16     instant_rate_est2   9.8743 0.1077
    ## 17     instant_rate_est3   6.1234 0.0749
    ## 18     instant_rate_est4   4.2713 0.0662
    ## 19     instant_rate_est5   3.3462 0.0657
    ## 20     instant_rate_est6   2.8909 0.0675
    ## 21     instant_rate_est7   2.6664 0.0695
    ## 22     instant_rate_est8   2.5556 0.0710
    ## 23     instant_rate_est9   2.5005 0.0719
    ## 24     instant_rate_var1   7.9998 0.9043
    ## 25     instant_rate_var2   2.1817 0.2303
    ## 26     instant_rate_var3   1.0062 0.1028
    ## 27     instant_rate_var4   0.8410 0.0889
    ## 28     instant_rate_var5   0.8612 0.0942
    ## 29     instant_rate_var6   0.8963 0.0997
    ## 30     instant_rate_var7   0.9197 0.1032
    ## 31     instant_rate_var8   0.9328 0.1051
    ## 32     instant_rate_var9   0.9396 0.1061
    ## 33   change_in_interval1  17.3646 0.2188
    ## 34   change_in_interval2   9.9819 0.1089
    ## 35   change_in_interval3   6.0788 0.0744
    ## 36   change_in_interval4   4.2857 0.0664
    ## 37   change_in_interval5   3.3823 0.0664
    ## 38   change_in_interval6   2.8844 0.0674
    ## 39   change_in_interval7   2.6750 0.0698
    ## 40   change_in_interval8   2.5210 0.0700
    ## 41   change_in_interval9   2.5012 0.0720
    ## 42 change_from_baseline1  17.3646 0.2188
    ## 43 change_from_baseline2  27.3465 0.3166
    ## 44 change_from_baseline3  33.4253 0.3693
    ## 45 change_from_baseline4  37.7110 0.4079
    ## 46 change_from_baseline5  41.0933 0.4442
    ## 47 change_from_baseline6  43.9777 0.4828
    ## 48 change_from_baseline7  46.6527 0.5262
    ## 49 change_from_baseline8  49.1736 0.5737
    ## 50 change_from_baseline9  51.6749 0.6263
    ## 51             residuals   0.9937 0.0405

## Jenss-Bayley Latent Change Model with Fixed Ratio of Growth Acceleration

``` r
source("LCSM_JB_fixed.R")
JB_LCSM_f <- getLCSM_JB_fixed(dat = JB_random_dat, T_records = 1:10, traj_var = "Y", t_var = "T",
                              paraNames = paraJB_LCSM_f)
JB_LCSM_f[[2]]
```

    ##                     Name Estimate     SE
    ## 1                 mueta0  49.9761 0.2710
    ## 2                 mueta1   2.4469 0.0724
    ## 3                 mueta2 -30.3222 0.4685
    ## 4                  gamma  -0.7022 0.0079
    ## 5                  psi00  13.7324 1.4585
    ## 6                  psi01   1.3163 0.2813
    ## 7                  psi02   4.9723 1.7685
    ## 8                  psi11   0.9252 0.0958
    ## 9                  psi12   1.9520 0.4673
    ## 10                 psi22  37.2419 4.0425
    ## 11     instant_rate_est1  17.4635 0.2177
    ## 12     instant_rate_est2   9.8735 0.1083
    ## 13     instant_rate_est3   6.1221 0.0760
    ## 14     instant_rate_est4   4.2705 0.0669
    ## 15     instant_rate_est5   3.3460 0.0658
    ## 16     instant_rate_est6   2.8911 0.0673
    ## 17     instant_rate_est7   2.6669 0.0691
    ## 18     instant_rate_est8   2.5563 0.0704
    ## 19     instant_rate_est9   2.5014 0.0712
    ## 20     instant_rate_var1   8.1256 0.8685
    ## 21     instant_rate_var2   2.2031 0.2297
    ## 22     instant_rate_var3   0.9992 0.1017
    ## 23     instant_rate_var4   0.8251 0.0839
    ## 24     instant_rate_var5   0.8422 0.0863
    ## 25     instant_rate_var6   0.8760 0.0902
    ## 26     instant_rate_var7   0.8988 0.0928
    ## 27     instant_rate_var8   0.9116 0.0943
    ## 28     instant_rate_var9   0.9183 0.0950
    ## 29   change_in_interval1  17.3676 0.2165
    ## 30   change_in_interval2   9.9812 0.1095
    ## 31   change_in_interval3   6.0776 0.0754
    ## 32   change_in_interval4   4.2849 0.0671
    ## 33   change_in_interval5   3.3821 0.0666
    ## 34   change_in_interval6   2.8846 0.0672
    ## 35   change_in_interval7   2.6755 0.0693
    ## 36   change_in_interval8   2.5217 0.0694
    ## 37   change_in_interval9   2.5021 0.0712
    ## 38 change_from_baseline1  17.3676 0.2165
    ## 39 change_from_baseline2  27.3488 0.3148
    ## 40 change_from_baseline3  33.4264 0.3684
    ## 41 change_from_baseline4  37.7113 0.4078
    ## 42 change_from_baseline5  41.0933 0.4443
    ## 43 change_from_baseline6  43.9779 0.4826
    ## 44 change_from_baseline7  46.6534 0.5254
    ## 45 change_from_baseline8  49.1751 0.5722
    ## 46 change_from_baseline9  51.6771 0.6238
    ## 47             residuals   0.9938 0.0376
