# LCSM_projects

## Part I: Jenss-Bayley Latent Change Score Model with Individual Ratio of Growth Acceleration in the Framework of Individual Measurement Occasions
**Manuscript Title:** <br>
Jenss-Bayley Latent Change Score Model with Individual Ratio of Growth Acceleration in the Framework of Individual
Measurement Occasions (accepted for publication in Journal of Educational and Behavioral Statistics)

**Description:** <br>
In this part, we extend an existing LCSM with the Jenss-Bayley growth curve and propose a novel expression of change scores that allows for (1) unequally-spaced study waves and (2) individual measurement occasions around each wave. We also extend the existing model to estimate the individual ratio of growth acceleration (that largely determines the trajectory shape and is viewed as the most important parameter in the Jenss-Bayley model). 
- Jenss-Bayley Latent Change Score Model with Random Ratio of Growth Acceleration
- Jenss-Bayley Latent Change Score Model with Fixed Ratio of Growth Acceleration

**Demo:** 
- [*R* package: *OpenMx*](https://github.com/Veronica0206/LCSM_projects/blob/main/Part%201/OpenMx_P1/OpenMx_demo.md)
(For OS, R version, and OpenMx version, see the demo)

**Example data:**
- [Data](https://github.com/Veronica0206/NonLinearCurve/blob/main/data/JB_random_dat.RData)
(For OS, R version, and OpenMx version, see the demo)

**Source Code:** <br>
***R package: OpenMx*** <br>
**The models developed in this project are now part of *R* package *NonLinearCurve* (dependency: *OpenMx*), where we provide functions capable of 'calculating' starting values from the input and generate the estimates described in the manuscript.**
- [Jenss-Bayley Latent Change Score Model with Random Ratio of Growth Acceleration](https://github.com/Veronica0206/NonLinearCurve/blob/main/R/LCSM_JB_random.R)
- [Jenss-Bayley Latent Change Score Model with Fixed Ratio of Growth Acceleration](https://github.com/Veronica0206/NonLinearCurve/blob/main/R/LCSM_JB_fixed.R)

***MPlus 8*** <br>
- [Jenss-Bayley Latent Change Score Model with Random Ratio of Growth Acceleration](https://github.com/Veronica0206/LCSM_projects/blob/main/Part%201/MPlus8_P1/JB_LCSM_random.inp)
- [Jenss-Bayley Latent Change Score Model with Fixed Ratio of Growth Acceleration](https://github.com/Veronica0206/LCSM_projects/blob/main/Part%201/MPlus8_P1/JB_LCSM_fixed.inp)

## Part II: Parallel Latent Basis Growth Model in the Framework of Individual Measurement Occasions
**Manuscript Title:** <br>
Extending Latent Basis Growth Model to Explore Joint Development in the Framework of Individual Measurement Occasions

**Description:** <br>
In this part, we propose a novel specification for LBGMs that allows for (1) unequally-spaced study waves and (2) individual measurement occasions around each wave. We then extend LBGMs to explore multiple repeated outcomes because longitudinal processes rarely unfold in isolation. 
- Parallel Latent Basis Growth Model

**Demo:** 
- [*R* package: *OpenMx*](https://github.com/Veronica0206/LCSM_projects/blob/main/Part%202/OpenMx_demo2.md)
(For OS, R version, and OpenMx version, see the demo)

**Source Code:** <br>
Will upload later.

## Part IV: Decomposing Impact on Longitudinal Outcome of Time-varying Covariate into Baseline Effect and Temporal Effect
**Manuscript Title:** <br>
Decomposing Impact on Longitudinal Outcome of Time-varying Covariate into Baseline Effect and Temporal Effect (accepted for publication in Journal of Educational and Behavioral Statistics)

**Description:** <br>
In this part, we propose decomposing the TVC’s effect into initial trait and temporal states using three methods to address this limitation. In each method, the baseline of the TVC is viewed as an initial trait, and the corresponding effects are obtained by regressing random intercepts and slopes on the baseline value. Temporal states are characterized as 
- interval-specific slopes,
- interval-specific changes,
- changes from the baseline at each measurement occasion.

**Demo:** 
- [*R* package: *OpenMx*](https://github.com/Veronica0206/LCSM_projects/blob/main/Part%204/OpenMx_E9/OpenMx_demo9.md)
(For OS, R version, and OpenMx version, see the demo)

**Example data:**
- [Data for Model 1](https://github.com/Veronica0206/LCSM_projects/blob/main/Part%204/OpenMx_E9/TVC1_BLS_dat.csv)
- [Data for Model 2](https://github.com/Veronica0206/LCSM_projects/blob/main/Part%204/OpenMx_E9/TVC2_BLS_dat.csv)
- [Data for Model 3](https://github.com/Veronica0206/LCSM_projects/blob/main/Part%204/OpenMx_E9/TVC3_BLS_dat.csv)
(For OS, R version, and OpenMx version, see the demo)

**Source Code:** <br>
***R package: OpenMx*** <br>
- [Decomposing TVC into baseline and interval-specific slopes](https://github.com/Veronica0206/LCSM_projects/blob/main/Part%204/OpenMx_E9/TVC1.R)
- [Decomposing TVC into baseline and interval-specific changes](https://github.com/Veronica0206/LCSM_projects/blob/main/Part%204/OpenMx_E9/TVC2.R)
- [Decomposing TVC into baseline and changes from the baseline at each measurement occasion](https://github.com/Veronica0206/LCSM_projects/blob/main/Part%204/OpenMx_E9/TVC3.R)

***MPlus 8*** <br>
- [Decomposing TVC into baseline and interval-specific slopes](https://github.com/Veronica0206/LCSM_projects/blob/main/Part%204/Mplus8_E9/BLSGM_TVCslp.inp)
- [Decomposing TVC into baseline and interval-specific changes](https://github.com/Veronica0206/LCSM_projects/blob/main/Part%204/Mplus8_E9/BLSGM_TVCslp.inp)
- [Decomposing TVC into baseline and changes from the baseline at each measurement occasion](https://github.com/Veronica0206/LCSM_projects/blob/main/Part%204/Mplus8_E9/BLSGM_TVCslp.inp)

## Part V: Further Exploration of the Effects of Time-varying Covariate in Growth Mixture Models with Nonlinear Trajectories
**Manuscript Title:** <br>
Further Exploration of the Effects of Time-varying Covariate in Growth Mixture Models with Nonlinear Trajectories (Behavior Research Methods)

**Description:** <br>
In this part, we extended Part IV to mixture modeling framework and proposed 3 mixture models with TVC. In each model, the baseline of the TVC is viewed as an initial trait, and the corresponding effects are obtained by regressing random intercepts and slopes on the baseline value. Temporal states are characterized as 
- interval-specific slopes,
- interval-specific changes,
- changes from the baseline at each measurement occasion.

**Demo:** 
- [*R* package: *OpenMx*](https://github.com/Veronica0206/LCSM_projects/blob/main/Part%205/OpenMx_E10/OpenMx_demo10.md)
(For OS, R version, and OpenMx version, see the demo)

**Example data:**
- [Data for Model 1](https://github.com/Veronica0206/LCSM_projects/blob/main/Part%205/OpenMx_E10/GMMTVC1_BLS_dat.csv)
- [Data for Model 2](https://github.com/Veronica0206/LCSM_projects/blob/main/Part%205/OpenMx_E10/GMMTVC2_BLS_dat.csv)
- [Data for Model 3](https://github.com/Veronica0206/LCSM_projects/blob/main/Part%205/OpenMx_E10/GMMTVC3_BLS_dat.csv)
(For OS, R version, and OpenMx version, see the demo)

**Source Code:** <br>
***R package: OpenMx*** <br>
- [Decomposing TVC into baseline and interval-specific slopes](https://github.com/Veronica0206/LCSM_projects/blob/main/Part%205/OpenMx_E10/GMM_TVC1.R)
- [Decomposing TVC into baseline and interval-specific changes](https://github.com/Veronica0206/LCSM_projects/blob/main/Part%205/OpenMx_E10/GMM_TVC2.R)
- [Decomposing TVC into baseline and changes from the baseline at each measurement occasion](https://github.com/Veronica0206/LCSM_projects/blob/main/Part%205/OpenMx_E10/GMM_TVC3.R)

***MPlus 8*** <br>
- [Decomposing TVC into baseline and interval-specific slopes](https://github.com/Veronica0206/LCSM_projects/blob/main/Part%205/Mplus8_E10/BLSGMM_TVCslp.inp)
- [Decomposing TVC into baseline and interval-specific changes](https://github.com/Veronica0206/LCSM_projects/blob/main/Part%205/Mplus8_E10/BLSGMM_TVCslp.inp)
- [Decomposing TVC into baseline and changes from the baseline at each measurement occasion](https://github.com/Veronica0206/LCSM_projects/blob/main/Part%205/Mplus8_E10/BLSGMM_TVCslp.inp)



