library(dplyr)
#mRNAsi
main.predict <- function( fnSig = "pcbc-stemsig.tsv", fnOut = "mRNA_StemScore.tsv" )
{
  ## Load the signature
  w <- read.delim( fnSig, header=FALSE, row.names=1 ) %>% as.matrix() %>% drop()
  
  ## Reduces HUGO|POSITION gene IDs to just HUGO
  f <- function( v ) unlist( lapply( strsplit( v, "\\|" ), "[[", 1 ) )
  
  s <- "Medulloblastoma_data/GSE85217_stemsig_rna_exp_matrix.txt"
  X <- read.delim( s, as.is=TRUE, check.names=FALSE ) %>%  ## Read the raw values
    filter( !grepl( "\\?", HGNC_symbol_from_ensemblv77 ) ) %>%      ## Drop genes with no mapping to HUGO
    mutate( HGNC_symbol_from_ensemblv77 = f( HGNC_symbol_from_ensemblv77 ) ) %>%        ## Clip gene ids to HUGO
    filter( HGNC_symbol_from_ensemblv77 %in% names(w) )         ## Reduce to the signature's gene set
  
  ## SLC35E2 has multiple entries with the same HUGO id
  ## Keep the first entry only
  j <- grep( "SLC35E2", X[,1] )
  if( length(j) > 1 )
    X <- X[-j[-1],]
  
  ## Convert to a matrix
  rownames(X) <- NULL
  X <- X %>% tibble::column_to_rownames( "HGNC_symbol_from_ensemblv77" ) %>% as.matrix()
  
  ## Reduce the signature to the common set of genes
  stopifnot( all( rownames(X) %in% names(w) ) )
  w <- w[ rownames(X) ]
  
  ####### Score via Spearman correlation
  s <- apply( X, 2, function(z) {cor( z, w, method = "sp", use = "complete.obs" )} )
  
  ## Scale the scores to be between 0 and 1
  s <- s - min(s)
  s <- s / max(s)
  
  write.table(cbind(s), file = fnOut, sep = "\t", quote = FALSE, col.names = FALSE)
}


#mDNAsi
load("pcbc-stemsig.p219.Rda")
w <- mm$w

da <- read.table("Medulloblastoma_data/GSE85212_stemsig_methylation.txt",sep="\t",header=T,row.names=1) 
idx <- intersect(rownames(da),names(w))
w <- w[as.character(idx)]
X <- da[as.character(idx),]
X <- as.matrix(X)

## Score via linear model
ss <- t(w) %*% X
## Scale the scores to be between 0 and 1
ss <- ss - min(ss)
ss <- ss / max(ss)
ss <- as.data.frame(t(ss))

#colnames(ss) <- "mDNAsi"
write.table(ss,file="mDNA_StemScore.tsv",sep="\t",quote=F,col.names=F) 
save(s,ss, file = "Medulloblastoma_StemScore.Rda")


