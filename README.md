# UBVSR: Update Bayesian Variable Selection Regression

Bayesian variable selection regression (BVSR) is able to jointly analyze genome-wide genetic data to produce the posterior probability of association for each covariate and estimate hyperparameters such as heritability and number of covariates that have nonzero effects. But the slow computation via MCMC hampered its wide-spread usage. In this package, we incorporate the Multiple-try importance tempering(MT-IT) method into BVSR and combined with Cholesky decomposition updating method when evaluating the posteriors, which is much faster than computing from scratch in certain cases.

## Installation
```{r eval = FALSE}
devtools::install_github("yangym58/UBVSR")
```

## Usage
```{r eval = False}
library(UBVSR)
```

### Input and output Specification
Function updatechol_add
- rtri is the right triangular matrix $R$ satisfies $R^tR=X^tX+\sigma^{-2}I$ for current $X$.
- sigma2 is $\sigma^2$.
- X is the design matrix(each column is a variable).
- indvec is a vector indicates which variable is contained in current X.
- newv is the variable to be added.
- output: The right triangular matrix satisfies $R'^tR'=X'^tX'+\sigma^{-2}I$ after adding new variable.

Function updatechol_remove
- rtri is the right triangular matrix $R$ satisfies $R^tR=X^tX+\sigma^{-2}I$ for current $X$.
- indvec is a vector indicates which variable is contained in current X.
- newv is the variable to be removed.
- output: The right triangular matrix satisfies $R'^tR'=X'^tX'+\sigma^{-2}I$ after removing a variable.

Function MT_IT 
- output:
- gammamat: Each row is a vector indicating whether the variable has an effect in an iteration of MT_IT.
- sigmavec: Each element is the hyperparameter in prior of beta in an iteration of MT_IT.
wmat: Each element is the weight in an iteration of MT_IT.

### Example
Using the function generate_data to generate data: the probability to have effect for each variable is 0.05. Heritability is 0.4. 1000 observations with 200 total variable.
```{r eval = False}
# Generate data
data <- generate_data(0.05, 0.4, 1000, 200)
```
Run the model:
```{r eval = False}
Y <- data$Y
X <- data$X
m <- 20
a <- 2
delta <- 0.5
niter <- 1500
res <- MT_IT(Y, X, m, a, delta, niter)
```
Evaluate the occuracy:
```{r eval = False}
calculate_accu(data$beta, res$gammamat)
```

## Learn more
For more information about BVSR:[BVSR](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6788783/).

For more information about MT_IT:[MT_IT](https://arxiv.org/abs/2304.06251).




