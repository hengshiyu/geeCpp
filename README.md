# geeCpp: a c++ package for implementation of generalized estimating equations (GEE)
Hengshi Yu, Tahmeed Tureen, Weiqi Zhou


## Introduction
cvcrand is an R package for the design and analysis of cluster randomized trials (CRTs). 

A cluster is the unit of randomization for a cluster randomized trial. Thus, when the number of clusters is small, there might be some baseline imbalance from the randomization between the arms. Constrained randomization constrained the randomization space. Given the baseline values of some cluster-level covariates, users can perform a constrained randomization on the clusters into two arms, with an optional input of user-defined weights on the covariates. 

In addition to the covariate-constrained randomization, covariate-by-covariate constrained randomization is also applicable for cluster randomized trails. User could directly perform constraints on each covariate of interest and randomize from the constrained space. 

At the end of the study, the individual outcome is collected. The `cvcrand` package then performs clustered permutation test on either continuous outcome or binary outcome adjusted for some individual-level covariates, producing p-value of the intervention effect.

## Functions and references
The cvcrand package constains three main functions. In the design of CRTs with two arms, users can use the `cvrall()` function to perform constrained randomization or use the `cvrcov` to perform covariate-by-covariate constrained randomization. And for the analysis part, user will use the `cptest()` function for clustered permutation test. 

1. cvrall function: covariate-constrained randomization for two-arm cluster randomized trials
    * Raab, G.M. and Butcher, I., 2001. Balance in cluster randomized trials. Statistics in medicine, 20(3), pp.351-365.
    * Li, F., Lokhnygina, Y., Murray, D.M., Heagerty, P.J. and DeLong, E.R., 2016. An evaluation of constrained randomization for the design and analysis of group‐randomized trials. Statistics in medicine, 35(10), pp.1565-1579.
    * Li, F., Turner, E. L., Heagerty, P. J., Murray, D. M., Vollmer, W. M., & DeLong, E. R. (2017). An evaluation of constrained randomization for the design and analysis of group‐randomized trials with binary outcomes. Statistics in medicine, 36(24), 3791-3806.
    * Dickinson, L. M., Beaty, B., Fox, C., Pace, W., Dickinson, W. P., Emsermann, C., & Kempe, A. (2015). Pragmatic cluster randomized trials using covariate constrained randomization: A method for practice-based research networks (PBRNs). The Journal of the American Board of Family Medicine, 28(5), 663-672.
    * Bailey, R. A., & Rowley, C. A. (1987). Valid randomization. Proceedings of the Royal Society of London. A. Mathematical and Physical Sciences, 410(1838), 105-124.
    
2. cvrcov function: covariate-by-covariate constrained randomization for two-arm cluster randomized trials

    * Greene, E. J. (2017). A SAS macro for covariate-constrained randomization of general cluster-randomized and unstratified designs. Journal of statistical software, 77(CS1).


2. cptest function: clustered permutation test for two-arm cluster randomized trial
    * Gail, M.H., Mark, S.D., Carroll, R.J., Green, S.B. and Pee, D., 1996. On design considerations and randomization‐based inference for community intervention trials. Statistics in medicine, 15(11), pp.1069-1092.
    * Li, F., Lokhnygina, Y., Murray, D.M., Heagerty, P.J. and DeLong, E.R., 2016. An evaluation of constrained randomization for the design and analysis of group‐randomized trials. Statistics in medicine, 35(10), pp.1565-1579.
    * Li, F., Turner, E. L., Heagerty, P. J., Murray, D. M., Vollmer, W. M., & DeLong, E. R. (2017). An evaluation of constrained randomization for the design and analysis of group‐randomized trials with binary outcomes. Statistics in medicine, 36(24), 3791-3806.
    * Dickinson, L. M., Beaty, B., Fox, C., Pace, W., Dickinson, W. P., Emsermann, C., & Kempe, A. (2015). Pragmatic cluster randomized trials using covariate constrained randomization: A method for practice-based research networks (PBRNs). The Journal of the American Board of Family Medicine, 28(5), 663-672.
    * Eldridge, S. M., Ukoumunne, O. C., & Carlin, J. B. (2009). The Intra‐Cluster Correlation Coefficient in Cluster Randomized Trials: A Review of Definitions. International Statistical Review, 77(3), 378-394.
    * Hannan, P. J., Murray, D. M., Jacobs Jr, D. R., & McGovern, P. G. (1994). Parameters to aid in the design and analysis of community trials: intraclass correlations from the Minnesota Heart Health Program. Epidemiology, 88-95.


### `cvrall()` example: covariate-constrained randomization

The balance score for constrained randomization in the program is developed from Raab and Butcher (2001).  

A study presented by Dickinson et al (2015) is about two approaches (interventions) for increasing the "up-to-date" immunization rate in 19- to 35-month-old children. They planned to randomize 16 counties in Colorado 1:1 to either a population-based approach or a practice-based approach. There are several county-level variables. The program will randomize on a subset of these variables. The continuous variable of average income is categorized to illustrate the use of the `cvrall()` on multi-category variables. And the percentage in Colorado Immunization Information System (CIIS) variable is truncated at 100%.

For the constrained randomization, we used the `cvrall()` function to randomize 8 out of the 16 counties into the practice-based. For the definition of the whole randomization space, if the total number of all possible schemes is smaller than `50,000`, we enumerate all the schemes as the whole randomization space. Otherwise, we simulate `50,000` schemes and choose the unique shemes among them as the whole randomization space. We calculate the balance scores of `"l2"` metric on three continuous covariates as well as two categorical covariates of location and income category. Location has `"Rural"` and `"Urban"`. The level of `"Rural"` was then dropped in `cvrall()`. As income category has three levels of `"low"`, `"med"`, and `"high"`,  the level of `"high"` was dropped to create dummy variables according to the alphanumerical order as well. Then we constrained the randomization space to the schemes with `"l2"` balance scores less than the `0.1` quantile of that in the whole randomization space. Finally, a randomization scheme is sampled from the constrained space.  

We saved the constrained randomization space in a CSV file in `"dickinson_constrained.csv"`, the first column of which is an indicator variable of the finally selected scheme (`1`) or not (`0`). We also have the output of a histogram displaying the distribution of all balance scores with a red line indicating our selected cutoff (the `0.1` quantile).


```r

Design_result <- cvrall(clustername = Dickinson_design$county,
                        balancemetric = "l2",
                        x = data.frame(Dickinson_design[ , c("location", "inciis",
                            "uptodateonimmunizations", "hispanic", "incomecat")]),
                        ntotal_cluster = 16,
                        ntrt_cluster = 8,
                        categorical = c("location", "incomecat"),
                        savedata = "dickinson_constrained.csv",
                        bhist = TRUE,
                        cutoff = 0.1,
                        seed = 12345, 
                        check_validity = TRUE)

```


### `cvrall()` example: stratified constrained randomization

User-defined weights can be used to induce stratification on one or more categorical variables. In the study presented by Dickinson et al (2015), there are 8 `"Urban"` and 8 `"Rural"` counties. A user-defined weight of `1,000` is added to the covariate of `location`, while these weights for other covariates are all `1`. Intuitively, a large weight assigned to a covariate sharply penalizes any imbalance of that covariates, therefore including schemes that are optimally balanced with respect to that covariate in the constrained randomization space. In practice, the resulting constrained space approximates the stratified randomization space on that covariate. In our illustrative data example, since half of the counties are located in rural areas, perfect balance is achieved by considering constrained randomization with the large weight for `location` variable. Alternatively, the option of `stratify` is able to perform the equivalent stratification on the stratifying variables specified.


```r
# Stratification on location

Design_stratified_result1 <- cvrall(clustername = Dickinson_design$county,
                                     balancemetric = "l2",
                                     x = data.frame(Dickinson_design[ , c("location", "inciis",
                                     "uptodateonimmunizations", 
                                     "hispanic", "incomecat")]),
                                     ntotal_cluster = 16,
                                     ntrt_cluster = 8,
                                     categorical = c("location", "incomecat"),
                                     weights = c(1000, 1, 1, 1, 1),
                                     cutoff = 0.1,
                                     seed = 12345) 
                                  
                                  
# An alternative and equivalent way to stratify on location

Design_stratified_result2 <- cvrall(clustername = Dickinson_design$county,
                                     balancemetric = "l2",
                                     x = data.frame(Dickinson_design[ , c("location", "inciis",
                                     "uptodateonimmunizations", 
                                     "hispanic", "incomecat")]),
                                     ntotal_cluster = 16,
                                     ntrt_cluster = 8,
                                     categorical = c("location", "incomecat"),
                                     stratify = "location",
                                     cutoff = 0.1,
                                     seed = 12345)

```

### `cvrcov()` example: covariate-by-covariate constrained randomization




For the constrained randomization, we used the `cvrcov()` function to randomize 8 out of the 16 counties into the practice-based. For the definition of the whole randomization space, if the total number of all possible schemes is smaller than `100,000`, we enumerate all the schemes as the whole randomization space. Otherwise, we simulate `100,000` unique schemes. Location has `"Rural"` and `"Urban"`. The level of `"Rural"` was then kept as 1 in `cvrcov()` and `"Urban"` is 0. Then we constrained the randomization space to have the schemes with absolute total difference of location be smaller than or equal to `5`, absolute mean difference of percentages of children ages 19-35 months in the CIIS less than or equal to `0.5` fraction of the overall mean, and absolute mean difference of income to be less than or equal to the `0.4` fraction of the overall mean. Finally, a randomization scheme is sampled from the constrained space. 

We saved the constrained randomization space in a CSV file in `"dickinson_cov_constrained.csv"`, the first column of which is an indicator variable of the finally selected scheme (`1`) or not (`0`). 


```r

Design_cov_result <- cvrcov(clustername = Dickinson_design_numeric$county,
                            x = data.frame(Dickinson_design_numeric[ , c("location", "inciis", 
                                                                          "uptodateonimmunizations", 
                                                                          "hispanic", "income")]),
                            ntotal_cluster = 16,
                            ntrt_cluster = 8,
                            constraints = c("s5", "mf.5", "any", "any", "mf0.4"), 
                            categorical = c("location"),
                            savedata = "dickinson_cov_constrained.csv",
                            seed = 12345, 
                            check_validity = TRUE)

```

### `cptest()` example: Clustered Permutation Test

At the end of cluster randomized trials, individual outcomes are collected. Permutation test based on Gail et al (1996) and Li et al (2016) is then applied to the continuous or binary outcome with some individual-level covariates. 

Suppose that the researchers were able to assess 300 children in each cluster in a study presented by Dickinson et al (2015), and the cluster randomized trial is processed with the selected randomization scheme from the example above of the `cvrall()` function. We expanded the values of the cluster-level covariates on the covariates' values of the individuals, according to which cluster they belong to. The correlated individual outcome of up-to-date on immunizations (`1`) or not (`0`) is then simulated using a generalized linear mixed model (GLMM) with a logistic link to induce correlation by including a random effect at the county level. The intracluster correlation (ICC) was set to be 0.01, using the latent response definition provided in Eldridge et al (2009). This is a reasonable value for population health studies Hannan et al (1994). We simulated one data set, with the outcome data dependent on the county-level covariates used in the constrained randomization design and a positive treatment effect so that the practice-based intervention increases up-to-date immunization rates more than the community-based intervention. For each individual child, the outcome is equal to `1` if he or she is up-to-date on immunizations and `0` otherwise. 

We used the `cptest()` function to process the clustered permutation test on the binary outcome of the status of up-to-date on immunizations. We input the file about the constrained space with the first column indicating the final scheme. The permutation test is on the continuous covariates of `"inciis"`, `"uptodateonimmunizations"`, `"hispanic"`, as well as categorical variables of `"location"` and `"incomecat"`. Location has `"Rural"` and `"Urban"`. The level of `"Rural"` was then dropped in `cptest()`. As income category has three levels of `"low"`, `"med"`, and `"high"`,  the level of `"high"` was dropped to create dummy variables according to the alphanumerical order as well.

```r
 Analysis_result <- cptest(outcome = Dickinson_outcome$outcome,
                           clustername = Dickinson_outcome$county,
                           z = data.frame(Dickinson_outcome[ , c("location", "inciis",
                               "uptodateonimmunizations", "hispanic", "incomecat")]), 
                           cspacedatname = "dickinson_constrained.csv",                                 
                           outcometype = "binary",                                                      
                           categorical = c("location","incomecat"))
```


