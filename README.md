---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "README-"
)
```


## Casecohort-pooling
  
*casecohort-pooling.R* is an R file for case-cohort pooling analysis as described in "Min Shi, David M. Umbach, Clarice R. Weinberg, Pooling biospecimens for exposure assessment for case-cohort analyses in cohort studies".

### Simulation study
**Sections 1-7 contain code for the simulation study.**

Section 1 sets up the necessary R libraries.

Section 2 is the code used for simulating case-cohort studies. Exposure has a log normal distribution. Age at enrollment can have four distributions as described in the paper: 
  1) uniform between 35 and 74.
  2) Gaussian with mean 55.2 and standard deviation 9.
  3) same as 2) but exposure is correlted with age at enrollment.
  4) Left skewed.  

Section 3 gives the pool forming strategy. It can be used to map single samples to pools using case-cohort data (ccoh.dat).

The data frame ccoh.dat contains case-cohort data. The required columns are:

ID:           sample ID variable
age.start:    age at enrollment
age.end:      age at the diagnosis for all cases (in or not in subcohort) and age at the end of followup for non-cases in the subcohort
status:       case (value=1) or non-case (value=0) status
subcohort:    subcohort status (TRUE/FALSE)        
 
The pool mapping will be saved in a text file "pool_mapping.txt". There are three vairables in this file:

ID:           sample ID as in ccoh.dat
pools:        pool ID 
p.size:       pool size

The file gives mapping of individual samples to pools. For example, as shown below pool 1 is a pool of size 4 with four samples in the pool: 1472, 6916, 10262 and 28880 while pool 3 is a singleton with only one sample 30940. 

   ID pools p.size
   
 1472     1      4
 
 6916     1      4
 
10262     1      4

28880     1      4

 5318     2      4
 
12968     2      4

15212     2      4

24344     2      4

30940     3      1

If pooled exposure was measured the values can be entered in the file.

Section 4 creates pools and partial pools. In addition to variables in ccoh.dat listed above there is a new variable "frac" which is the weight of each sample. For pools with reduced number of individuals due to drop out the weight will be reduced accordingly. This variables will be fed to the "weights" option in logistic regression. Note that pools of different sizes were given different ages so that they can have different stratification parameters in Section 6.

Section 5 creates per-year format data based on the start and end of followup time.

Section 6 calls logistic regression to do the analysis. The parameter of interest is "E".

Section 7 runs bootstrapping to get variable estimate.

### Functions to form pools and call logistic regression to analyze data
For actual study one can use function *create.pools* to form pools and function *cch.logistic* to call logistic regression to analyze the data. Please check *casecohort-pooling.R* for example code.
