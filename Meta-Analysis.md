Table of contents
-----------------
1/ Function to convert tau and rho in r and to compute Fisher Z score.

---------------------------------------------------------------------------------------------------

Functions
---------
1/ Function to convert tau and rho in r and to compute Fisher Z score.

For maths, see [Walker, 2003, Journal of Modern Applied Statistical Methods](https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&cad=rja&uact=8&ved=2ahUKEwiLr6L1847iAhUC8hoKHdm2BYUQFjAAegQIARAC&url=http%3A%2F%2Fwww.cedu.niu.edu%2F~walker%2Fpersonal%2FWalker%2520Kendall%2527s%2520Tau.pdf&usg=AOvVaw0qglHyyIwZV-so3y07CyCM) and [Zimmerman, Zumbo, & Williams, 2003, Psicología ](https://www.researchgate.net/publication/26421626_Bias_in_Estimation_and_Hypothesis_Testing_of_Correlation)

rho.tau.r: a numerical vector with the Spearman rho, Kendall tau, Pearson r correlation coefficients   
n: an numerical vector with the sample size associated with the correlation coefficients   
type: a string vector with the used methods to compute the correlation coefficients, must be "Kendall" or  "Pearson" or "Spearman"   

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


