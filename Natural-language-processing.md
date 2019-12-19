Table of contents
-----------------
1/ Function to lemmatize corpus

---------------------------------------------------------------------------------------------------

Functions  
---------  
1/ Function to lemmatize corpus  

!! treetagger had to be installed. https://www.cis.uni-muenchen.de/~schmid/tools/TreeTagger/   

!! to know available langages : `r koRpus::available.koRpus.lang()` keep the last two letters : fr for frenche, en for english

+ corpus : a corpus of tm package (! different of qunateda) 
+ treetaggerfilepath : the filepath of treetagger library, e.g. "~/Programs/Treetagger/"
+ unknown.words : if TRUE keep the words which are unknown of the lexic. If FALSE, there are deleted. 

```r
lemmatization.corpus <- function( corpus, lang = "fr",
                                  treetaggerfilepath = "~/Programs/Treetagger/",
                                  unknown.words = TRUE )
{
  # require package
  require( quanteda ) ; require( koRpus )
  
  require( paste0( "koRpus.lang.", lang ), character.only = TRUE )
  
  # function
  String.Concatener <- function( vector, separator )
  {   
    # exclude NA
    na.omit( vector ) -> vector ;
    
    # selection of length for loop
    if ( length( vector ) >= 2) 
    { 
      # stock fisrt element
      tmp.val1 <- vector[1] ;
      
      # loop to concatenate the nect elements
      for ( a in 2:length( vector ) ) 
      { 
        tmp.val1 <- paste( tmp.val1, vector[a], sep = separator ) ; 
      } ;
    } 
    else
    {
      tmp.val1 <- vector ;
    } ;
    
    return( tmp.val1 ) ;
  } ;
  
  # Job
  if ( unknown.words )
  {
    lapply( 1:length(corpus), 
            function( a, corpus) 
            {
              
              vector1 <- NULL 
              
              lapply( 1:length( corpus[[a]][[1]] ),  
                      function(b, corpus){
                        corpus[[a]][[1]][b] %>%
                          quanteda::tokens(.) %>% 
                          as.character(.) %>%
                          koRpus::treetag( ., 
                                           treetagger = "manual", 
                                           format = "obj", TT.tknz = FALSE, 
                                           lang = lang, 
                                           debug = TRUE,
                                           TT.options = list(path = treetaggerfilepath, 
                                                             preset = lang ) ) ->
                          lemma ;
                        
                        lemma[ lemma[, "lemma"] == "<unknown>","token"] ->
                          lemma[ lemma[, "lemma"] == "<unknown>","lemma"] ;
                        
                        lemma[ is.na( lemma[, "lemma"] ) ,"token"] ->
                          lemma[is.na( lemma[, "lemma"] ) ,"lemma"] ;
                        
                        lemma[, "lemma"] %>%
                          String.Concatener( ., " " ) ->
                          lemma ; 
                        
                        return( lemma )
                      }, 
                      corpus) %>%
                c( vector1, .) ->
                vector1 ;
              
              return( unlist( vector1 ) )
            },
            corpus ) ->
      corpus[[a]][[1]] ;
  } else
  {
    lapply( 1:length(corpus), 
            function( a, corpus) 
            {
              
              vector1 <- NULL 
              
              lapply( 1:length( corpus[[a]][[1]] ),  
                      function(b, corpus){
                        corpus[[a]][[1]][b] %>%
                          quanteda::tokens(.) %>% 
                          as.character(.) %>%
                          koRpus::treetag( ., 
                                           treetagger = "manual", 
                                           format = "obj", TT.tknz = FALSE, 
                                           lang = lang, 
                                           debug = TRUE,
                                           TT.options = list(path = treetaggerfilepath, 
                                                             preset = lang ) ) ->
                          lemma ;
                        
                        lemma[ lemma[, "lemma"] == "<unknown>","token"] ->
                          lemma[ lemma[, "lemma"] == "<unknown>","lemma"] ;
                        
                        lemma[ is.na( lemma[, "lemma"] ) ,"token"] ->
                          lemma[is.na( lemma[, "lemma"] ) ,"lemma"] ;
                        
                        lemma[, "lemma"] %>%
                          String.Concatener( ., " " ) ->
                          lemma ; 
                        
                        return( lemma )
                      }, 
                      corpus) %>%
                c( vector1, .) ->
                vector1 ;
              
              return( unlist( vector1 ) )
            },
            corpus ) ->
      corpus[[a]][[1]] ;
  }
  
  return( corpus )
}
```