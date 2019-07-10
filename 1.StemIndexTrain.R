library(gelnet)
library(dplyr)
library(gdata)
library(biomaRt)
library(synapseClient)
#mRNAsi
synapseLogin()
genes2hugo <- function( v, srcType = "ensembl_gene_id" )
{
    ## Retrieve the EMSEMBL -> HUGO mapping
    ensembl <- biomaRt::useMart( "ENSEMBL_MART_ENSEMBL", host="www.ensembl.org", dataset="hsapiens_gene_ensembl" )
    ID <- biomaRt::getBM( attributes=c(srcType, "hgnc_symbol"), filters=srcType, values=v, mart=ensembl )

    ## Make sure there was at least one mapping
    if( nrow(ID) < 1 ) stop( "No IDs mapped successfully" )
    
    ## Drop empty duds
    j <- which( ID[,2] == "" )
    if( length(j) > 0 ) ID <- ID[-j,]
    stopifnot( all( ID[,1] %in% v ) )

    ID
}

main.train <- function( fnOut = "pcbc-stemsig.tsv", fnGenes = NULL )
{
  ## Load RNAseq data
#  synRNA <- synGet( "syn2701943", downloadLocation = "/data/PCBC" )
#  X <- read.delim( synRNA@filePath ) %>%
  X <- read.delim( "PCBC/rnaseq_norm.tsv" ) %>%
    tibble::column_to_rownames( "tracking_id" ) %>%
    as.matrix()
  
  ## Retrieve metadata
  synMeta <- synTableQuery( "SELECT UID, Diffname_short FROM syn3156503" )
  Y <- synMeta@values %>%
    mutate( UID = gsub("-", ".", UID) ) %>%
    tibble::column_to_rownames( "UID" )
  
  ## Retrieve the labels from the metadata
  y <- Y[colnames(X),]
  names(y) <- colnames(X)
  
  ## Fix the missing labels by hand
  y["SC11.014BEB.133.5.6.11"] <- "EB"
  y["SC12.039ECTO.420.436.92.16"] <- "ECTO"
  
  ## Drop the splice form ID from the gene names
  v <- strsplit( rownames(X), "\\." ) %>% lapply( "[[", 1 ) %>% unlist()
  rownames(X) <- v
  
  ## Map Ensembl IDs to HUGO
  V <- genes2hugo( rownames(X) )
  X <- X[V[,1],]
  rownames(X) <- V[,2]
  
  ## Reduce the gene set to the provided list (if applicable)
  if( is.null( fnGenes ) == FALSE )
  {
    vGenes <- read.delim( fnGenes, header=FALSE ) %>% as.matrix() %>% drop()
    VE <- genes2hugo( vGenes, "entrezgene" )
    X <- X[intersect( rownames(X), VE[,2] ),]
  }
  
  ## Mean-center the data
  m <- apply( X, 1, mean )
  X <- X - m
  
  ## Identify stem cell samples
  j <- which( y == "SC" )
  X.tr <- X[,j]
  X.bk <- X[,-j]
  
  ## Train a one-class model
  mm <- gelnet( t(X.tr), NULL, 0, 1 )
  
  ## Store the signature to a file
  write.table(mm$w, file = fnOut, sep = "\t", quote = FALSE, col.names = FALSE)
  
  ## Perform leave-one-out cross-validation
  auc <- c()
  for( i in 1:ncol(X.tr) )
  {
    ## Train a model on non-left-out data
    X1 <- X.tr[,-i]
    m1 <- gelnet( t(X1), NULL, 0, 1 )
    
    ## Score the left-out sample against the background
    s.bk <- apply( X.bk, 2, function(z) {cor( m1$w, z, method="sp" )} )
    s1 <- cor( m1$w, X.tr[,i], method="sp" )
    
    ## AUC = P( left-out sample is scored above the background )
    auc[i] <- sum( s1 > s.bk ) / length(s.bk)
    cat( "Current AUC: ", auc[i], "\n" )
    cat( "Average AUC: ", mean(auc), "\n" )
  }
  
  return(auc)
}

main.train()


#mDNAsi
load("pcbc.pd.f.Rda")
load("pcbc.data.Rda")
m <- apply(pcbc.data, 1, mean )
pcbc.data.2 <- pcbc.data - m
M1_smp <- pcbc.pd.f[pcbc.pd.f$Diffname_short %in% "SC",] #SC
M2_smp <- pcbc.pd.f[!(pcbc.pd.f$Diffname_short %in% "SC"),] #non-SC
X.tr <- pcbc.data.2[, as.character(M1_smp$UID)] # 44 samples
X.bk <- pcbc.data.2[, as.character(M2_smp$UID)] # 55 samples
mm <- gelnet(t(X.tr), NULL, 0, 1) # NULL for a one-class task 
save( mm, file = "pcbc-stemsig.p219.Rda")
