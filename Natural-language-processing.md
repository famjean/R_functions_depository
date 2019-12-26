Table of contents
-----------------
1/ Function to lemmatize corpus
2/ Function to clean text (remove special characters, numbers, uppercases, punctuation and stripWhitespace, +- remove stopwords, +- remove defined words)
3/ Function to describe the corpus (number of words, sentences, paragraphs, TTR and spartisity)
4/ Function to class terms by tf*idf

---------------------------------------------------------------------------------------------------

Functions  
---------  
1/ Function to lemmatize corpus  

!! treetagger had to be installed. https://www.cis.uni-muenchen.de/~schmid/tools/TreeTagger/   

!! to know available langages : `r koRpus::available.koRpus.lang()` keep the last two letters : fr for frenche, en for english

+ corpus : a corpus of tm package (! different of quanteda) 
+ treetaggerfilepath : the filepath of treetagger library, e.g. "~/Programs/Treetagger/"
+ unknown.words : if TRUE keep the words which are unknown of the lexic. If FALSE, there are deleted. 

```r
lemmatization.corpus <- function( corpus, lang = "fr",
                                  treetaggerfilepath = "~/Programs/Treetagger/",
                                  unknown.words = TRUE )
{
  # install packages
  load.packages <- function( PackagesNames, ..., 
                             install = TRUE, 
                             load = TRUE )
  {
    # Get a vector with packages names 
    packages1 <- c( PackagesNames, ... ) ;
    
    # Require or install and require packages
    lapply( 1:length( packages1 ),
            function(a)
            {
              if ( install ) 
              { if ( !packages1[a] %in% installed.packages()[,1] ) 
              { install.packages( packages1[a] ) } }
              if ( load ) { require( packages1[a], character.only = TRUE  )  }
            } ) -> tmp
  } 
  
  load.packages( c( "quanteda", "tm") , install = TRUE, load = FALSE )
  
  # load packages
  load.packages( c( "koRpus", paste0( "koRpus.lang.", lang ), "dplyr" ) , install = TRUE, load = TRUE )
  
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
   
   for ( a in  1:length(corpus) )
   {
     for ( b in 1:length( corpus[[a]][[1]] ) )
     {
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
                        lemma -> 
                          corpus[[a]][[1]][b] ; 
     }
   }
   
  } else
  {
    for ( a in  1:length(corpus) )
   {
     for ( b in 1:length( corpus[[a]][[1]] ) )
     {
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
              
                        lemma[, "lemma"] %>%
                          String.Concatener( ., " " ) ->
                              lemma ; 
                        lemma -> 
                          corpus[[a]][[1]][b] ; 
     }
   }
  }
  
  return( corpus )
}
```

-------------------------------------------------------------------

2/ Function to clean text (remove special characters, numbers, uppercases, punctuation and stripWhitespace, +- remove stopwords, +- remove defined words)

+ corpus : a corpus of tm package (! different of quanteda) 
+ langage : langage of the corpus, put two letters inside "", see tm::stopwords to explore available langages.
+ RMstopwords : if TRUE, remove the stopwords, if FALSE, let its. Stop words are generally the most common words in a language.
+ RMpunctuation : if TRUE, remove the punctuation, if FALSE, let it.
+ WordsToRM : a string with the words to delete, if NULL, it doesn't delete words.

```r
TextCleaning <- function( corpus, 
                          language = "en", 
                          RMstopwords = TRUE,
                          RMpunctuation = TRUE,
                          WordsToRM = NULL
                         )
{
  # install/load packages
  load.packages <- function( PackagesNames, ..., 
                             install = TRUE, 
                             load = TRUE )
  {
    # Get a vector with packages names 
    packages1 <- c( PackagesNames, ... ) ;
    
    # Require or install and require packages
    lapply( 1:length( packages1 ),
            function(a)
            {
              if ( install ) 
              { if ( !packages1[a] %in% installed.packages()[,1] ) 
              { install.packages( packages1[a] ) } }
              if ( load ) { require( packages1[a], character.only = TRUE  )  }
            } ) -> tmp
  } 
  
  load.packages( c( "tm") , install = TRUE, load = FALSE )
  load.packages( c( "dplyr" ) , install = TRUE, load = TRUE )
  
  # domestic functions 
  toSpace <- tm::content_transformer( function(x, pattern) { return (gsub(pattern, " ", x) ) } ) 
  
  # job
  if ( RMstopwords ) 
  {
    corpus %>%
      tm::tm_map(., tm::removeWords, tm::stopwords(language) ) -> 
      corpus ;
  }
  
  if ( RMpunctuation ) 
  {
    corpus %>%
      tm::tm_map(., tm::removePunctuation) ->
      corpus ;
  }
  
  
  if ( !is.null(WordsToRM) ) 
    { corpus %>% 
      tm::tm_map(., tm::removeWords, WordsToRM ) -> 
      corpus ; } ;
  
  
  corpus %>% 
    tm::tm_map(., toSpace, "/") %>%
    tm::tm_map(., toSpace, "@") %>% 
    tm::tm_map(., toSpace, "\\|") %>%
    tm::tm_map(., toSpace, "#") %>% 
    tm::tm_map(., toSpace, "$") %>% 
    tm::tm_map(., toSpace, "€" ) %>%
    tm::tm_map(., toSpace, ":") %>% 
    tm::tm_map(., toSpace, "-") %>% 
    tm::tm_map(., toSpace, "'") %>%
    tm::tm_map(., toSpace, "’") %>%
    tm::tm_map(., toSpace, "\\(") %>%
    tm::tm_map(., toSpace, "\\)") %>%
    tm::tm_map(., toSpace, "\\{") %>%
    tm::tm_map(., toSpace, "\\}") %>%
    tm::tm_map(., toSpace, "\\[") %>%
    tm::tm_map(., toSpace, "\\]") %>%
    tm::tm_map(., tm::removeNumbers) %>% 
    tm::tm_map(., tm::content_transformer(tolower) ) %>% 
    tm::tm_map(., tm::stripWhitespace) -> 
    corpus ;
  
  return( corpus ) ;
} 
```
-------------------------------------------------------------------

3/ Function to describe the corpus (number of words, sentences, paragraphs, TTR and spartisity)

+ corpus : a corpus of tm package (! different of quanteda) 
+ doc.name : definie a string with the names of documents, put NULL to do nothing
+ round : number of decimals for indicators, put NULL to do nothing
```r
Corpus.Description <- function( corpus, 
                                doc.name = NULL, 
                                round = 2 )
{
  # install packages
  load.packages <- function( PackagesNames, ..., 
                             install = TRUE, 
                             load = TRUE )
  {
    # Get a vector with packages names 
    packages1 <- c( PackagesNames, ... ) ;
    
    # Require or install and require packages
    lapply( 1:length( packages1 ),
            function(a)
            {
              if ( install ) 
              { if ( !packages1[a] %in% installed.packages()[,1] ) 
              { install.packages( packages1[a] ) } }
              if ( load ) { require( packages1[a], character.only = TRUE  )  }
            } ) -> tmp
  } 
  
  load.packages( c( "quanteda", "tm") , install = TRUE, load = FALSE )
  
  # load packages
  load.packages( c( "dplyr" ) , install =   TRUE, load = TRUE )
  
  # domestic functions 
  as.numeric.data.frame <- function( dataframe )
  {
    for ( a in 1:dim( dataframe )[2] )
    {
      if ( sum( is.na( as.numeric( as.character( dataframe[,a] ) ) ) ) == 
           sum( is.na( dataframe[,a] ) ) )
      { 
        dataframe[,a] <- as.numeric( as.character( dataframe[,a] ) )
      } ;
    } ;
    return( dataframe ) ;
  } ;

  
  # job
  corpus %>% 
    lapply( ., as.character) %>% 
    lapply( ., length) %>% 
    unlist(.) %>% 
    cbind( document = names(.) , paragraphs = .) -> 
    paragraphs ;
  
  corpus %>%
    lapply( ., as.character ) %>%
    lapply( ., function(x) gsub( ",|’|-|—|_|(|)|[|]|/", "", x ) ) %>%
    lapply( ., function(x) unlist( strsplit( x, "\\.|\\!|\\?|\\;|:" ) ) ) %>%
    lapply( ., length ) %>%
    unlist(.) %>% 
    cbind( document = names(.) , sentences = .) ->
    sentences ;
  
  corpus %>%
    tm::TermDocumentMatrix(.) %>%
    as.matrix(.) %>% 
    colSums(.) %>%  
    cbind( document = names(.), words = .) ->
    words ;
  
  corpus %>%
    quanteda::corpus(.) %>%
    quanteda::dfm(.) %>%
    quanteda::textstat_lexdiv(.) %>%
    mutate( ., document = gsub("text", "", .$document ) ) ->
    lexdiv ;
  
  merge( words, sentences, by = "document" ) %>%
    merge( ., paragraphs, by = "document" ) %>%
    merge( ., lexdiv, by = "document" ) ->
    output ;
  
  output %>%
    as.numeric.data.frame(.) ->
    output ;
  
  output$document %>%
    as.character(.) ->
    output$document ;
  
  output[,2:4] %>% 
    colSums(.) -> 
    tmp ;
  
  corpus %>%
    quanteda::corpus(.) %>%
    quanteda::dfm(.) %>% 
    as.matrix()%>% 
    apply( ., 2, sum) %>% 
    matrix( ., nrow = 1, dimnames = list( c( "All documents"), names(.) ) ) %>% 
    quanteda::as.dfm(.) %>% 
    quanteda::textstat_lexdiv(.) %>%
    dplyr::select( -1 ) %>%
    as.numeric(.) %>%
    c( tmp, . ) ->
    tmp ;
  
  rbind( output,
         c( "Total", tmp) 
         ) %>%
    as.numeric.data.frame(.) ->
    output ;
  
  corpus %>%
    tm::TermDocumentMatrix(.) %>% 
    as.matrix(.) %>%
    apply( .,
           2,
           function(x) { length( x[ x == 0 ] ) / length( x ) } ) %>%
    c( .,
       quanteda::sparsity( quanteda::dfm( quanteda::corpus( corpus ) ) ) ) ->
    sparsity ;
  
  output %>%
    cbind.data.frame( .,
                      sparsity = sparsity ) ->
    output ;
    
  if ( !is.null( round ) ) 
  {
    output[,-1] %>%
      round( ., round ) ->
      output[,-1] ;
  }
  
    if ( !is.null( doc.name ) ) 
  {
    output[,1] <- c( doc.name, "Total" ) ;
  }
  
  NULL -> rownames( output ) ; 
  
  return( output )
}
```

-------------------------------------------------------------------

4/ Function to class terms by tf*idf

+ tdm : a term document matrix of tm package (! different of quanteda) 
```r
tf.idf.rank <- function( tdm )
{
  # install packages
  load.packages <- function( PackagesNames, ..., 
                             install = TRUE, 
                             load = TRUE )
  {
    # Get a vector with packages names 
    packages1 <- c( PackagesNames, ... ) ;
    
    # Require or install and require packages
    lapply( 1:length( packages1 ),
            function(a)
            {
              if ( install ) 
              { if ( !packages1[a] %in% installed.packages()[,1] ) 
              { install.packages( packages1[a] ) } }
              if ( load ) { require( packages1[a], character.only = TRUE  )  }
            } ) -> tmp
  } 
  
  # load packages
  load.packages( c( "dplyr" ) , install =   TRUE, load = TRUE )
  
  tdm %>% 
    as.matrix(.) ->
    tdm ;
  
  # Total count 
  tdm %>% 
    data.frame( ., total = rowSums(.) ) %>%
    dplyr::select( ., total) ->
      total ;
  
  # tf
  tdm %>% 
    colSums() ->
    Total.words ;
  
  tdm -> tdm.tfprop ;
  for (a in 1:dim(tdm)[2] )
  {
    tdm[,a] / Total.words[a] ->
      tdm.tfprop[,a] ;
  }
  
  paste( colnames(tdm), 
         "tfprop", 
         sep = "." ) ->
    colnames( tdm.tfprop ) ;
  
  # idf 
  dim( tdm )[2] -> nb.doc ;
  
  NULL -> idf ;
  for (a in 1:dim(tdm)[1] )
  {
    tdm[a,] > 0 ->
      nb.doc.with.term ;
    
    sum( nb.doc.with.term ) ->
      nb.doc.with.term ;
    
    c( idf, log( nb.doc / nb.doc.with.term ) ) ->
      idf
  }
    
  # tf * idf
  tdm -> tdm.tf.x.idf ;
  for (a in 1:dim(tdm)[2] )
  {
    tdm[,a] * idf ->
      tdm.tf.x.idf[,a] ;
  }
  
  paste( colnames(tdm), 
         "tf.x.idf", 
         sep = "." ) ->
    colnames( tdm.tf.x.idf ) ;
  
  # tfprop * idf
  tdm -> tdm.tfprop.x.idf ;
  for (a in 1:dim(tdm)[2] )
  {
    tdm.tfprop[,a] * idf ->
      tdm.tfprop.x.idf[,a] ;
  }
  
  paste( colnames(tdm), 
         "tfprop.x.idf", 
         sep = "." ) ->
    colnames( tdm.tfprop.x.idf ) ;
  
  # rank
  tdm -> tdm.rank ;
  for (a in 1:dim(tdm)[2] )
  {
    rank( tdm[,a] ) ->
      tdm.rank[,a] ;
  }
  
  paste( colnames(tdm), 
         "rank", 
         sep = "." ) ->
    colnames( tdm.rank ) ;
  
  # fusion
  cbind.data.frame( terms = rownames( tdm ),
                    total,
                    idf = idf,
                    tdm,
                    tdm.tfprop,
                    tdm.tf.x.idf,
                    tdm.tfprop.x.idf,
                    tdm.rank ) ->
    output ;

  rownames( tdm ) ->
    rownames( output ) ;
  
  return( output )
}
```

-------------------------------------------------------------------