Table of contents
-----------------
1/ Function to convert tau and rho in r and to compute Fisher Z score.   
2/ Function to compute power   
3/ Function to perform meta-analysis with studies with missing data   
   
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
   
3/ Function to perform meta-analysis with studies with missing data   

For maths and details on functions, see [Schartzer, Carpenter, and Rücker, 2014, book: Leta-Analysis with R, Springer](https://www.springer.com/gb/book/9783319214153)

+ For n.e, mean.e, sd.e, n.c, mean.c, sd.c, data, see ?meta::metacont.   
+ p.miss.e: percentage of missing data in experimental arm If there is no missing data in study, put 0.      
+ p.miss.c: percentage of missing data in control arm If there is no missing data in study, put 0.      
+ data: An optional data frame containing the study information.   
+ studlab: An optional vector with study labels.   
+ mu.e: mean difference between missing group and observed group for the studies in the experimental arm.   
+ nu.e: variance of differences between missing group and observed group for the studies in the experimental arm.  
+ mu.c: mean difference between missing group and observed group for the studies in the control arm.   
+ nu.c.c: variance of differences between missing group and observed group for the studies in the control arm.   
+ rho: correlation between 1/ the average difference between missing group and observed group for the studies in the experimental arm, and 2/ between missing group and observed group for the studies in the control arm.   
+ ...: arguments to pass to metacont. See ?meta::metacont.   

Details:   
  There are four strategies to deal with missing data. 
  + Fixed equal: mu.e = mu.c = a determined value (e.g. 8 or 20); nu.e = nu.c = 0 ; rho = 0   
  +
  +
  +
     
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
   
