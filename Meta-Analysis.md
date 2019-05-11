Table of contents
-----------------
1/ Function to convert tau and rho in r and to compute Fisher Z score.
2/ Function to compute power   
   
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

