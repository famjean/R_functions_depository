Table of contents
-----------------
1/ Function to convert tau and rho in r and to compute Fisher Z score.   
2/ Function to compute power   
3/ Function to perform meta-analysis with studies with missing data (mean difference)   
4/ Function to perform meta-analysis with studies with missing data (generalization)   
5/ Function to test variance between studies (between means)   
6/ Function to test variance between studies (generalization)  
7/ Function to perform meta-imputation (between means)     
8/ Function to perform meta-imputation (generalization)   
9/ Function to prepare data for multivariate meta-analysis

---------------------------------------------------------------------------------------------------

Functions
---------
1/ Function to convert tau and rho in r and to compute Fisher Z score.

For maths, see [Walker, 2003, Journal of Modern Applied Statistical Methods](https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&cad=rja&uact=8&ved=2ahUKEwiLr6L1847iAhUC8hoKHdm2BYUQFjAAegQIARAC&url=http%3A%2F%2Fwww.cedu.niu.edu%2F~walker%2Fpersonal%2FWalker%2520Kendall%2527s%2520Tau.pdf&usg=AOvVaw0qglHyyIwZV-so3y07CyCM) and [Zimmerman, Zumbo, & Williams, 2003, Psicología ](https://www.researchgate.net/publication/26421626_Bias_in_Estimation_and_Hypothesis_Testing_of_Correlation)

+ rho.tau.r: a numerical vector with the Spearman rho, Kendall tau, Pearson r correlation coefficients   
+ n: a numerical vector with the sample size associated with the correlation coefficients   
+ type: a string vector with the used methods to compute the correlation coefficients, must be "Kendall" or  "Pearson" or "Spearman"   

```r
corToFisherZ <- function( rho.tau.r, n, type )
{
  # Checking
  if ( !any( type %in% c( "Pearson", "Spearman", "Kendall" ) ) ) 
    { stop( "type should be Pearson or Spearman or Kendall" ) }
  
  # convert rho and tau
  rho.tau.r -> transformed.rho.tau.And.r ; 
  sin( rho.tau.r[ type == "Kendall" ] * pi * 1/2 ) -> transformed.rho.tau.And.r[ type == "Kendall" ] ;
  (6 / ( pi * (n[ type == "Spearman"] + 1))) * (asin( rho.tau.r[ type == "Spearman" ] ) + ((n[ type == "Spearman"] - 2) * asin( rho.tau.r[ type == "Spearman" ] / 2)) ) -> 
    transformed.rho.tau.And.r[ type == "Spearman" ] ;
  
  # Fisher Z
  atanh( transformed.rho.tau.And.r ) -> fisher.Z  ; 
  sqrt( 1/ (n - 3) ) ->  fisher.Z.se ;
  
  # Fisher Z special Method
  cat( "Special method: Fisher Z transformation applied on raw Spearman rho (not on r transformation) and the Fisher Z se for Spearman original values is sqrt( 1.060 / (n - 3) )\n\n" )
  
  atanh( transformed.rho.tau.And.r ) -> fisher.Z.specialmethod ;
  atanh( rho.tau.r[ type == "Spearman" ] ) -> fisher.Z.specialmethod[ type == "Spearman" ] ;
  
  as.numeric( rep( NA, length( rho.tau.r ) ) ) -> fisher.Z.se.specialmethod ;
  sqrt( 1/ (n[ type == "Pearson"] - 3) ) ->  fisher.Z.se.specialmethod[ type == "Pearson"] ;
  sqrt( 1/ (n[ type == "Kendall"] - 3) ) ->  fisher.Z.se.specialmethod[ type == "Kendall"] ;
  sqrt( 1.060 / (n[ type == "Spearman"] - 3) ) ->  fisher.Z.se.specialmethod[ type == "Spearman"] ;
  
  # output
  data.frame( n = n,
              type = type,
              rho.tau.r = rho.tau.r,
              r = transformed.rho.tau.And.r,
              fisher.Z = fisher.Z,
              fisher.Z.se = fisher.Z.se,
              fisher.Z.specialmethod = fisher.Z.specialmethod, 
              fisher.Z.se = fisher.Z.se.specialmethod
              )
}

# Test 
corToFisherZ( c( -0.5, -0.5, -0.5, 0, 0, 0, 0.5, 0.5, 0.5 ), 
              rep(20,9), 
              c( "Pearson", "Spearman", "Kendall", "Pearson", "Spearman",  "Kendall", "Pearson", "Spearman", "Kendall" ) )
```

NB: perform a metaregression on type to test whether transformations bias the meta-analysis.

-------------------------------------------------------------------

2/ Function to compute power   
   
For math an source see [Harrer, Cuijpers, and Ebert, 2019, Online book, Doing Meta-Analysis in R: A Hand-on Guiden chap 14](https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/power-analysis.html) 
I had generalized the formula by let the user put the ESaverafeSE.  

+ theta: wanting result of the meta-analysis, common effect of the studies  
+ k: number of studies in the meta-analysis  
+ ES.average.SE: the average standard error of the analyzed Effect Size. Be carrefull with the transformation in log, etc. 
+ Heterogeneity: hetérogeneity of the studies. "high", "moderate", "low". Default is "high".
+ alpha: level of risk of type II. 

```r
power.meta.fixed.test <-function( 
  theta, k, ES.average.SE, heterogeneity = "high", alpha = 0.05 )
{
  # Checking
  if ( !any( heterogeneity %in% c( "low", "moderate", "high" ) ) ) 
  { stop( "heterogeneity should be low or moderate or high" ) }
  
  # Prepare computations
  zval <- qnorm( 1 - alpha / 2 , 0, 1 )
  theta.f.var <- 1 / sum( (1 / ES.average.SE^2 ) * k )
  H <- switch( heterogeneity, "low" = 1.33, "moderate" = 1.67, "high" =  2 )
  theta.r.var <- theta.f.var * H
  lambda.f <- theta / sqrt( theta.f.var )
  lambda.r <- theta / sqrt( theta.r.var )
  
  # compute powers
  power.fixed <- 1 - pnorm( zval - lambda.f ) + pnorm( -zval - lambda.f ) 
  power.random <- 1 - pnorm( zval - lambda.r ) + pnorm( -zval - lambda.r ) 
  
  list( power.fixed = power.fixed,
        power.random = power.random )
}
```
    
-------------------------------------------------------------------
   
3/ Function to perform meta-analysis with studies with missing data (mean difference)   

For maths and details on functions, see [Schartzer, Carpenter, and Rücker, 2014, book: Leta-Analysis with R, Springer](https://www.springer.com/gb/book/9783319214153)

+ For n.e, mean.e, sd.e, n.c, mean.c, sd.c, see ?meta::metacont.   
+ p.miss.e: percentage of missing data in experimental arm If there is no missing data in study, put 0.      
+ p.miss.c: percentage of missing data in control arm If there is no missing data in study, put 0.      
+ data: An optional data frame containing the study information.   
+ studlab: An optional vector with study labels.   
+ mu.e: mean difference between missing group and observed group for the studies in the experimental arm.   
+ nu.e: standard error of differences between missing group and observed group for the studies in the experimental arm.  
+ mu.c: mean difference between missing group and observed group for the studies in the control arm.   
+ nu.c.c: standard error of differences between missing group and observed group for the studies in the control arm.   
+ rho: correlation between 1/ the average difference between missing group and observed group for the studies in the experimental arm, and 2/ between missing group and observed group for the studies in the control arm.   
+ ...: arguments to pass to metacont. See ?meta::metacont.   

Details:   
  There are four strategies to deal with missing data. 
  + Fixed equal: mu.e = mu.c = a determined value (e.g. 8 or 20); nu.e = nu.c = 0 ; rho = 0   
  + Fixed opposite: mu.e =a determined value (e.g. 8 or 20) ; mu.c = -mu.e ; nu.e = nu.c = 0 ; rho = 0   
  + Random equal: mu.e = mu.c = a determined value (e.g. 8 or 20) ; nu.e = nu.c = a determined value (e.g. sqrt(5) or sqrt(15)) ; rho = 1  
  + Random opposite (uncorrelated):  mu.e =a determined value (e.g. 8 or 20) ; mu.c = -mu.e ; nu.e = nu.c = a determined value (e.g. sqrt(5) or sqrt(15)) ; rho = 0   
     
```r
metacont.miss <- function( 
  n.e, mean.e, sd.e, p.miss.e, 
  n.c, mean.c, sd.c, p.miss.c,
  data = NULL, studlab,
  mu.e, nu.e, mu.c, nu.c, rho,
  ...
  )
{
  # require package
  require( meta )
  
  # internal functions
  chknull <- function(x, name = NULL) {
    ##
    ## Check whether argument is NULL
    ##
    if (is.null(name))
      name <- deparse(substitute(x))
    ##
    if (is.null(x))
      stop("Argument '", name, "' is NULL.", call. = FALSE)
    ##
    invisible(NULL)
  }
  
  setstudlab <- function(x, k) {
    ##
    ## Set study labels
    ##
    if (is.null(x))
      x <- seq(along = rep(1, k))
    if (is.factor(x))
      x <- as.character(x)
    ##
    x
  }
  
  # check input
  nulldata <- is.null(data)
  if (nulldata) 
    data <- sys.frame(sys.parent())
  mf <- match.call()
  n.e <- eval(mf[[match("n.e", names(mf))]], data, enclos = sys.frame(sys.parent()))
  chknull(n.e)
  k.All <- length(n.e)
  mean.e <- eval(mf[[match("mean.e", names(mf))]], data, enclos = sys.frame(sys.parent()))
  chknull(mean.e)
  sd.e <- eval(mf[[match("sd.e", names(mf))]], data, enclos = sys.frame(sys.parent()))
  chknull(sd.e)
  n.c <- eval(mf[[match("n.c", names(mf))]], data, enclos = sys.frame(sys.parent()))
  chknull(n.c)
  mean.c <- eval(mf[[match("mean.c", names(mf))]], data, enclos = sys.frame(sys.parent()))
  chknull(mean.c)
  sd.c <- eval(mf[[match("sd.c", names(mf))]], data, enclos = sys.frame(sys.parent()))
  chknull(sd.c)
  studlab <- eval(mf[[match("studlab", names(mf))]], data, 
                  enclos = sys.frame(sys.parent()))
  studlab <- setstudlab(studlab, k.All)
  
  # function for se of studies with missing data
  semiss <- function(seTE, n.e, p.miss.e, n.c, p.miss.c,
                     mu.e, nu.e, mu.c, nu.c, rho){
    V1 <- function(p.miss.e, p.miss.c, nu.e, nu.c, rho)
      nu.e^2*p.miss.e^2 + nu.c^2*p.miss.c^2 - 2*rho*nu.e*nu.c*p.miss.e*p.miss.c
    V2 <- function(n.e, p.miss.e, n.c, p.miss.c, mu.e, nu.e, mu.c, nu.c)
      (mu.e^2+nu.e^2)*p.miss.e*(1-p.miss.e)/n.e +
      (mu.c^2+nu.c^2)*p.miss.c*(1-p.miss.c)/n.c
    sqrt(seTE^2 +
           V1(p.miss.e, p.miss.c, nu.e, nu.c, rho) +
           V2(n.e, p.miss.e, n.e, p.miss.c,
              mu.e, nu.e, mu.c, nu.c))
  }
  
  # Classic meta 
  mm1 <- metacont(  
             n.e = n.e, mean.e = mean.e, sd.e = sd.e, 
             n.c = n.c, mean.c = mean.c, sd.c = sd.c,
             data = data, ... )
  data$TE <- mm1$TE
  data$seTE <- mm1$seTE
  
  # Adjustment for missing data
  data$TE.miss <- with( data, TE + mu.e * p.miss.e - mu.c * p.miss.c )
  data$seTE.miss <- with( data,
                            semiss( seTE = seTE, 
                                    n.e = n.e, p.miss.e = p.miss.e, 
                                    n.c = n.c, p.miss.c = p.miss.c,
                                    mu.e = mu.e, nu.e = nu.e, 
                                    mu.c = mu.c, nu.c = nu.c, 
                                    rho = rho ) )
  # Relaunch meta
  metagen( TE = TE.miss, seTE = seTE.miss, 
           data = data, ... )
}
```
   
-------------------------------------------------------------------
   
4/ Function to perform meta-analysis with studies with missing data (generalization)  

For maths and details on functions, see [Schartzer, Carpenter, and Rücker, 2014, book: Leta-Analysis with R, Springer](https://www.springer.com/gb/book/9783319214153)

+ For TE, seTE, n.e, n.c, see ?meta::metagen. For OR, put logOR and SE of logOR.  
+ p.miss.e: percentage of missing data in experimental arm If there is no missing data in study, put 0.      
+ p.miss.c: percentage of missing data in control arm If there is no missing data in study, put 0.      
+ data: An optional data frame containing the study information.   
+ studlab: An optional vector with study labels.   
+ mu.e: mean difference between missing group and observed group for the studies in the experimental arm.   
+ nu.e: standard error of differences between missing group and observed group for the studies in the experimental arm.  
+ mu.c: mean difference between missing group and observed group for the studies in the control arm.   
+ nu.c.c: standard error of differences between missing group and observed group for the studies in the control arm.   
+ rho: correlation between 1/ the average difference between missing group and observed group for the studies in the experimental arm, and 2/ between missing group and observed group for the studies in the control arm.   
+ ...: arguments to pass to metacont. See ?meta::metagen.   

Details:   
  There are four strategies to deal with missing data. 
  + Fixed equal: mu.e = mu.c = a determined value (e.g. 8 or 20); nu.e = nu.c = 0 ; rho = 0   
  + Fixed opposite: mu.e =a determined value (e.g. 8 or 20) ; mu.c = -mu.e ; nu.e = nu.c = 0 ; rho = 0   
  + Random equal: mu.e = mu.c = a determined value (e.g. 8 or 20) ; nu.e = nu.c = a determined value (e.g. sqrt(5) or sqrt(15)) ; rho = 1  
  + Random opposite (uncorrelated):  mu.e =a determined value (e.g. 8 or 20) ; mu.c = -mu.e ; nu.e = nu.c = a determined value (e.g. sqrt(5) or sqrt(15)) ; rho = 0   
     
```r
metagen.miss <- function( 
  n.e, mean.e, sd.e, p.miss.e, 
  n.c, mean.c, sd.c, p.miss.c,
  data = NULL, studlab,
  mu.e, nu.e, mu.c, nu.c, rho,
  ...
  )
{
  # require package
  require( meta )
  
  # internal functions
  chknull <- function(x, name = NULL) {
    ##
    ## Check whether argument is NULL
    ##
    if (is.null(name))
      name <- deparse(substitute(x))
    ##
    if (is.null(x))
      stop("Argument '", name, "' is NULL.", call. = FALSE)
    ##
    invisible(NULL)
  }
  
  setstudlab <- function(x, k) {
    ##
    ## Set study labels
    ##
    if (is.null(x))
      x <- seq(along = rep(1, k))
    if (is.factor(x))
      x <- as.character(x)
    ##
    x
  }
  
  # check input
  nulldata <- is.null(data)
  if (nulldata) 
    data <- sys.frame(sys.parent())
  mf <- match.call()
  TE <- eval(mf[[match("TE", names(mf))]], data, enclos = sys.frame(sys.parent()))
  seTE <- eval(mf[[match("seTE", names(mf))]], data, enclos = sys.frame(sys.parent()))
  n.e <- eval(mf[[match("n.e", names(mf))]], data, enclos = sys.frame(sys.parent()))
  chknull(n.e)
  k.All <- length(n.e)
  n.c <- eval(mf[[match("n.c", names(mf))]], data, enclos = sys.frame(sys.parent()))
  chknull(n.c)
  studlab <- eval(mf[[match("studlab", names(mf))]], data, 
                  enclos = sys.frame(sys.parent()))
  studlab <- setstudlab(studlab, k.All)
  
  # function for se of studies with missing data
  semiss <- function(seTE, n.e, p.miss.e, n.c, p.miss.c,
                     mu.e, nu.e, mu.c, nu.c, rho){
    V1 <- function(p.miss.e, p.miss.c, nu.e, nu.c, rho)
      nu.e^2*p.miss.e^2 + nu.c^2*p.miss.c^2 - 2*rho*nu.e*nu.c*p.miss.e*p.miss.c
    V2 <- function(n.e, p.miss.e, n.c, p.miss.c, mu.e, nu.e, mu.c, nu.c)
      (mu.e^2+nu.e^2)*p.miss.e*(1-p.miss.e)/n.e +
      (mu.c^2+nu.c^2)*p.miss.c*(1-p.miss.c)/n.c
    sqrt(seTE^2 +
           V1(p.miss.e, p.miss.c, nu.e, nu.c, rho) +
           V2(n.e, p.miss.e, n.e, p.miss.c,
              mu.e, nu.e, mu.c, nu.c))
  }
  
  # Adjustment for missing data
  data$TE.miss <- with( data, TE + mu.e * p.miss.e - mu.c * p.miss.c )
  data$seTE.miss <- with( data,
                            semiss( seTE = seTE, 
                                    n.e = n.e, p.miss.e = p.miss.e, 
                                    n.c = n.c, p.miss.c = p.miss.c,
                                    mu.e = mu.e, nu.e = nu.e, 
                                    mu.c = mu.c, nu.c = nu.c, 
                                    rho = rho ) )
  # Launch meta
  metagen( TE = TE.miss, seTE = seTE.miss, 
           data = data, ... )
}
```
NB: The TE must be approximately normally distribued. e.g. for logOR, such analyse could not be applied for rare events.   

-------------------------------------------------------------------
   
5/ Function to test variance between studies (between means)   
   
+ Se: Standard deviation in experimental arm.    
+ Sc: Standard deviation in control arm.    
+ Ne: Number of observations in experimental arm.   
+ Nc: Number of observations in control arm.   
+ data: An optional data frame containing the study information.   
+ studlab: An optional vector with study labels.   
+ round: number of decimal to round.   
   
```r
meta.variance.test.means <- function( 
  Se, Sc, Ne, Nc, 
  data = NULL, studlab,
  round = 3
)
{
  # internal functions
  setstudlab <- function(x, k) {
    ##
    ## Set study labels
    ##
    if (is.null(x))
      x <- seq(along = rep(1, k))
    if (is.factor(x))
      x <- as.character(x)
    ##
    x
  }
  
  # Check input
  nulldata <- is.null(data)
  if (nulldata) 
    data <- sys.frame(sys.parent())
  mf <- match.call()
  Ne <- eval(mf[[match("Nc", names(mf))]], data, enclos = sys.frame(sys.parent()))
  Nc <- eval(mf[[match("Ne", names(mf))]], data, enclos = sys.frame(sys.parent()))
  Se <- eval(mf[[match("Se", names(mf))]], data, enclos = sys.frame(sys.parent()))
  Sc <- eval(mf[[match("Sc", names(mf))]], data, enclos = sys.frame(sys.parent()))
  k.All <- length(Ne)
  studlab <- eval(mf[[match("studlab", names(mf))]], data, 
                  enclos = sys.frame(sys.parent()))
  studlab <- setstudlab(studlab, k.All)
  
  # Test variances
  Tests.e <- matrix( rep( NA, length(Ne)^2 ), ncol = length(Ne) )
  Tests.c <- matrix( rep( NA, length(Ne)^2 ), ncol = length(Ne) )
  
  for (a in 1:length(Ne) )
  {
    f.tests <- c(pf(Se^2/Se[a]^2, Ne-1, Ne[a]-1),
                 pf(Sc^2/Sc[a]^2, Nc-1, Nc[a]-1))
    f.tests <- matrix(f.tests, ncol=2)
    Tests.e[,a] <- f.tests[,1]
    Tests.c[,a] <- f.tests[,2]
  }
  
  round( Tests.e, round ) -> Tests.e
  round( Tests.c, round ) -> Tests.c
  rownames( Tests.e ) <- colnames( Tests.e ) <- studlab
  rownames( Tests.c ) <- colnames( Tests.c ) <- studlab
  
  list( "F Var test on experimental arm" = Tests.e,
        "F Var test on control arm" = Tests.c )
}
```
   
-------------------------------------------------------------------
   
6/ Function to test variance between studies (generalization)  
   
+ seES: Standard error of the effect size.      
+ N: Number of observations in the study.   
+ data: An optional data frame containing the study information.   
+ studlab: An optional vector with study labels.   
+ round: number of decimal to round.   
   
```r
meta.variance.test <- function( 
  seES, N, 
  data = NULL, studlab,
  round = 3
)
{
  # internal functions
  setstudlab <- function(x, k) {
    ##
    ## Set study labels
    ##
    if (is.null(x))
      x <- seq(along = rep(1, k))
    if (is.factor(x))
      x <- as.character(x)
    ##
    x
  }
  
  # Check input
  nulldata <- is.null(data)
  if (nulldata) 
    data <- sys.frame(sys.parent())
  mf <- match.call()
  seES <- eval(mf[[match("seES", names(mf))]], data, enclos = sys.frame(sys.parent()))
  N <- eval(mf[[match("N", names(mf))]], data, enclos = sys.frame(sys.parent()))
  k.All <- length(N)
  studlab <- eval(mf[[match("studlab", names(mf))]], data, 
                  enclos = sys.frame(sys.parent()))
  studlab <- setstudlab(studlab, k.All)
  
  # Test variances
  Tests <- matrix( rep( NA, length(N)^2 ), ncol = length(N) )
  
  for (a in 1:length(N) )
  {
    Tests[,a] <- pf( seES^2 / seES[a]^2, N-1, N[a]-1 )
  }
  
  round( Tests, round ) -> Tests
  rownames( Tests ) <- colnames( Tests ) <- studlab
  
  list( "F Var test" = Tests )
}
```
   
-------------------------------------------------------------------
   
  7/ Function to perform meta-imputation (between means)   
   
For maths and details on functions, see Schartzer, Carpenter, and Rücker, 2014, book: Leta-Analysis with R, Springer   
   
Perform a test on variances before applying imputation and exclude study very different variances.   
   
+ missing: a string vector with the names of study with missing seES.
+ excluded: a string vector with the names of excluded study for imputation. If not, let "".  
+ Me: Mean in experimental arm.    
+ Mc: Mean in control arm.   
+ Se: Standard deviation in experimental arm.    
+ Sc: Standard deviation in control arm.    
+ Ne: Number of observations in experimental arm.   
+ Nc: Number of observations in control arm. 
+ data: An optional data frame containing the study information.   
+ B: number of simulation.    
+ seed: definie number to fix the random.  
  
```r
meta.imputation.means <- function(
  studies, missing, excluded = "", 
  Me, Mc, Se, Sc, Ne, Nc,
  data = NULL, 
  B = 10000,  
  seed = 1234
  )
{
  # load package
  require( meta ) 
  
  # Set seed
  set.seed( seed )
  
  # check input
  nulldata <- is.null(data)
  if (nulldata) 
    data <- sys.frame(sys.parent())
  mf <- match.call()
  Me <- eval(mf[[match("Me", names(mf))]], data, enclos = sys.frame(sys.parent()))
  Mc <- eval(mf[[match("Mc", names(mf))]], data, enclos = sys.frame(sys.parent()))
  Ne <- eval(mf[[match("Nc", names(mf))]], data, enclos = sys.frame(sys.parent()))
  Nc <- eval(mf[[match("Ne", names(mf))]], data, enclos = sys.frame(sys.parent()))
  Se <- eval(mf[[match("Se", names(mf))]], data, enclos = sys.frame(sys.parent()))
  Sc <- eval(mf[[match("Sc", names(mf))]], data, enclos = sys.frame(sys.parent()))
  
  # select studies
  toimpute <- !(  studies %in% c( missing, excluded ) )
  tobeimputed <- studies %in% missing 
  
  # Form pooled estimate of variability:
  S2.e <- sum( ( Ne[toimpute] - 1 ) * Se[toimpute]^2 )
  S2.c <- sum( ( Nc[toimpute] - 1 ) * Sc[toimpute]^2 )
  
  # Calculate degrees of freedom:
  df.e <- sum( Ne[toimpute] - 1 )
  df.c <- sum( Nc[toimpute] - 1 )
  
  # Stock sd for computations
  Se.miss <- Se
  Sc.miss <- Sc
  
  # Prepare object to stock results
  stocker <- data.frame( TE.fixed = rep( NA, B), seTE.fixed = NA,
                         TE.random = NA, seTE.random = NA,
                         tau = NA )
  
  # Run multiple imputations
  for (b in 1:B) 
  {
    # Draw sigma2.e, sigma2.c
    sigma2.e <- S2.e / rchisq( 1, df = df.e)
    sigma2.c <- S2.c / rchisq( 1, df = df.c)
    
    # Draw standard deviations for missiing
    sd.e.miss <- sigma2.e *
      rchisq( 1, df = Ne[ tobeimputed ] - 1 ) / ( Ne[ tobeimputed ] - 1 )
    sd.c.miss <- sigma2.c *
      rchisq( 1, df = Nc[ tobeimputed ] - 1 ) / ( Nc[ tobeimputed ] - 1 )
    
    # Meta analyses of current imputed dataset
    Se.miss[ tobeimputed ] <- sqrt( sd.e.miss )
    Sc.miss[ tobeimputed ] <- sqrt( sd.e.miss )
    
    # Run meta
    m.miss <- metacont( n.e = Ne, mean.e = Me, sd.e = Se.miss,
                        n.c = Nc, mean.c = Mc, sd.c = Sc.miss,
                        data = data )
    # Store results
    stocker$TE.fixed[b]    <- m.miss$TE.fixed 
    stocker$seTE.fixed[b]  <- m.miss$seTE.fixed
    stocker$TE.random[b]   <- m.miss$TE.random
    stocker$seTE.random[b] <- m.miss$seTE.random
    stocker$tau[b]         <- m.miss$tau
  }
  
  # Calculate between and within variances
  s2.b.fixed  <- var( stocker$TE.fixed )
  s2.b.random <- var( stocker$TE.random )
  s2.w.fixed  <- mean( stocker$seTE.fixed^2 )
  s2.w.random <- mean( stocker$seTE.random^2  )
 
  # Fixed effect estimate using multiple imputation
  TE.fixed.imp   <- mean( stocker$TE.fixed )
  seTE.fixed.imp <- sqrt( var( stocker$TE.fixed ) * (1 + 1/B) +
                           mean( stocker$seTE.fixed^2 ) )
  
  # Random effects estimate using multiple imputation
  TE.random.imp   <- mean( stocker$TE.random )
  seTE.random.imp <- sqrt( var( stocker$TE.random ) * (1 + 1/B) +
                            mean( stocker$seTE.random^2 ) )
  
  # Calculate degrees of freedom
  df.fixed  <- (B - 1) * (1 + s2.w.fixed / ( (1 + 1 / B) * s2.b.fixed ) )^2
  df.random <- (B - 1) * (1 + s2.w.random / ( (1 + 1 / B) * s2.b.random ) )^2
  
  cbind( Type = c( "Fixed", "Random" ),
    rbind(
    data.frame( ci( TE.fixed.imp, seTE.fixed.imp, df =  df.fixed) )[c(1:5,8,6)],
    data.frame( ci( TE.random.imp, seTE.random.imp, df =  df.random) )[c(1:5,8,6)]
  ) )
}
```
   
---------------------------------------------------------------------------------------------------
  
8/ Function to perform meta-imputation (generalization)   
   
For maths and details on functions, see Schartzer, Carpenter, and Rücker, 2014, book: Leta-Analysis with R, Springer   

Perform a test on variances before applying imputation and exclude study very different variances.   

+ missing: a string vector with the names of study with missing seES.
+ excluded: a string vector with the names of excluded study for imputation. If not, let "".  
+ ES: Effect size. 
+ seES: Standard error of the effect size.   
+ N: Number of observations.   
+ data:An optional data frame containing the study information.   
+ B: number of simulation.    
+ seed: definie number to fix the random.  

```r
meta.imputation <- function(
  studies, missing, excluded = "", 
  ES, seES, N, 
  data = NULL, 
  B = 10000,  
  seed = 1234
)
{
  # load package
  require( meta ) 
  
  # Set seed
  set.seed( seed )
  
  # check input
  nulldata <- is.null(data)
  if (nulldata) 
    data <- sys.frame(sys.parent())
  mf <- match.call()
  ES <- eval(mf[[match("ES", names(mf))]], data, enclos = sys.frame(sys.parent()))
  seES <- eval(mf[[match("seES", names(mf))]], data, enclos = sys.frame(sys.parent()))
  N <- eval(mf[[match("N", names(mf))]], data, enclos = sys.frame(sys.parent()))
  
  # select studies
  toimpute <- !(  studies %in% c( missing, excluded ) )
  tobeimputed <- studies %in% missing 
  
  # Form pooled estimate of variability:
  S2 <- sum( ( N[toimpute] - 1 ) * seES[toimpute]^2 )
  
  # Calculate degrees of freedom:
  df <- sum( N[toimpute] - 1 )
  
  # Stock sd for computations
  S.miss <- seES
  
  # Prepare object to stock results
  stocker <- data.frame( TE.fixed = rep( NA, B), seTE.fixed = NA,
                         TE.random = NA, seTE.random = NA,
                         tau = NA )
  
  # Run multiple imputations
  for (b in 1:B) 
  {
    # Draw sigma2.
    sigma2 <- S2 / rchisq( 1, df = df)
    
    # Draw standard deviations for missing
    sd.miss <- sigma2 *
      rchisq( 1, df = N[ tobeimputed ] - 1 ) / ( N[ tobeimputed ] - 1 )
    
    # Meta analyses of current imputed dataset
    S.miss[ tobeimputed ] <- sqrt( sd.miss )
    
    # Run meta
    m.miss <- metagen( TE = ES, seTE = S.miss,
                       data = data )
                        
    # Store results
    stocker$TE.fixed[b]    <- m.miss$TE.fixed 
    stocker$seTE.fixed[b]  <- m.miss$seTE.fixed
    stocker$TE.random[b]   <- m.miss$TE.random
    stocker$seTE.random[b] <- m.miss$seTE.random
    stocker$tau[b]         <- m.miss$tau
  }
  
  # Calculate between and within variances
  s2.b.fixed  <- var( stocker$TE.fixed )
  s2.b.random <- var( stocker$TE.random )
  s2.w.fixed  <- mean( stocker$seTE.fixed^2 )
  s2.w.random <- mean( stocker$seTE.random^2  )
  
  # Fixed effect estimate using multiple imputation
  TE.fixed.imp   <- mean( stocker$TE.fixed )
  seTE.fixed.imp <- sqrt( var( stocker$TE.fixed ) * (1 + 1/B) +
                            mean( stocker$seTE.fixed^2 ) )
  
  # Random effects estimate using multiple imputation
  TE.random.imp   <- mean( stocker$TE.random )
  seTE.random.imp <- sqrt( var( stocker$TE.random ) * (1 + 1/B) +
                             mean( stocker$seTE.random^2 ) )
  
  # Calculate degrees of freedom
  df.fixed  <- (B - 1) * (1 + s2.w.fixed / ( (1 + 1 / B) * s2.b.fixed ) )^2
  df.random <- (B - 1) * (1 + s2.w.random / ( (1 + 1 / B) * s2.b.random ) )^2
  
  cbind( Type = c( "Fixed", "Random" ),
         rbind(
           data.frame( ci( TE.fixed.imp, seTE.fixed.imp, df =  df.fixed) )[c(1:5,8,6)],
           data.frame( ci( TE.random.imp, seTE.random.imp, df =  df.random) )[c(1:5,8,6)]
         ) )
}
```
   
---------------------------------------------------------------------------------------------------
  

9/ Function to prepare data for multivariate meta-analysis
   
For maths and details on functions, see Schartzer, Carpenter, and Rücker, 2014, book: Leta-Analysis with R, Springer   

+ ES: A data.frame with effect size for the evaluations. E.g cbind( score1, score2, ..., scoreN).   
+ Var: A data.frame of variance of the evaluations. E.g cbind( var(score1), var(score2), ..., var(score) )    
+ Covar:  A data.frame of covariance of the evaluations. E.g cbind( covar(score1, score2), covar(score1, score3), ..., covar(score1, scoreN), covar(score2, score3), ..., covar(score2, scoreN), ..., covar(scoreN-1, scoreN) ). NULL if you impute covar with rho.   
+ rho: NULL if covar are know. If not, a data.frame with correlation between the evaluations. If details are not know, impute a specific value. E.g data.frame( cor(score1,score2) = c(0.3,...,0.4), cor(score1, score3) = 0.7, ..., cor(score1, scoreN) = 0.8, cor(score2, score3) = 0.3, ... cor(score2, scoreN), cor(score3,scor4), ...,cor( scoreN-1, scoreN) = 0.3)      
   
```r
preparation_meta.multivariate <- function(
 ES, Var, Covar = NULL, rho = NULL
)
{
   # Internal function 
   cor.sdTOcov <- function(sd1, sd2, cor) {sd1*sd2*cor}
   
   # Check input
    # Compute covar
   if ( is.null( rho ) & is.null( Covar ) ) { stop( "rho or Covar must be filled" ) }
   if ( is.null( Covar ) )
   {
      if ( ncol(rho ) !=  sum( ncol(Var) - 1:ncol(Var) ) ) { stop( paste( "There must be",  sum( ncol(Var) - 1:ncol(Var) ), "covars or rho" ) ) }
      
      Covar <- data.frame( matrix( rep( NA, sum( ncol(Var) - 1:ncol(Var) ) * nrow(Var) ), nrow = nrow(Var) ) )
      for (a in 1:(ncol(Var) - 1) )
      {
         for (b in (a + 1):ncol(Var) )
         {
            c <- sum( ncol(Var) - 1:a ) - (ncol(Var) - a) + (b - a)
            Covar[c] <- cor.sdTOcov( Var[a], Var[b], rho[c] ) 
         }
      }
   }
   
   if ( ncol( Covar) !=  sum( ncol(Var) - 1:ncol(Var) ) ) { stop( paste( "There must be",  sum( ncol(Var) - 1:ncol(Var) ), "covars or rho" ) ) }
   
   # Make VarCovar data.frame
   VarCovar <- NULL
  
   for (a in 1:(ncol(Var) - 1) )
   {
      VarCovar <- cbind( VarCovar, Var[,a] )
      for (b in (a + 1):ncol(Var) )
      {
         c <- sum( ncol(Var) - 1:a ) - (ncol(Var) - a) + (b - a)
         VarCovar <- cbind( VarCovar, Covar[,c] )
      }
   }
   VarCovar <- cbind( VarCovar, Var[,ncol(Var)] )

   list( ES = as.matrix(ES), VarCovar = as.matrix(VarCovar) )
}
```
   
---------------------------------------------------------------------------------------------------
  
      
