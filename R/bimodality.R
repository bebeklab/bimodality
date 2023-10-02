library(Biobase)
library(assertthat)
library(dplyr)
library(impute)
library(GEOquery)

setwd("/Users/alpun/Desktop/Bimodality/")

assayData_function <- function(txt_file){
  read.delim(txt_file, sep="\t", row.names = 1, header=T)
} #reads assay data text file (expression set) and transforms it into a data frame

featureData_function <- function(txt_file) {
  read.delim(txt_file, sep=c("\t"), row.names = 1, header=T)
} #reads feature data text file (fdata) and transforms it into a data frame

phenoData_function <- function(txt_file) {
  read.delim(txt_file, sep="\t", row.names = 1, header=T)
} #reads pheno data text file (pdata) and transforms it into a data frame

collapseProbes <- function (x, method = 'mean', topTable = NULL, colNameOfStat = NULL, colNameOfGeneSymbol = 'SYMBOL'){

  # Validate parameters ----------

  if ( !(method %in% list('variance',
                          'logFC','P.Value',
                          'adj.P.Value', 'mean')) ) {
    stop("Method ", paste(method, "not supported"))
  }

  # Make sure that a topTable is provided
  if (is.data.frame(topTable) ){

    # Check that the metric is one of the columns
    if(!(method %in% list("variance",
                          "logFC", "P.Value",
                          "adj.P.Value")) ) {
      stop("Method ", paste(method, "not supported") )
    }

    # Make sure the column name that has the gene symbol is provided
    if (!(colNameOfGeneSymbol %in% colnames(topTable)) ) {
      stop(paste(colNameOfGeneSymbol, 'must be a column in topTable'))
    }
  }


  # Main Workflow ----------
  if (!is.matrix(x))  {

    # Get Matrix from eset
    mat_exprs <- exprs(x)

  }


  if (method == 'variance') {
    cat("Processing eset using:",method,"\n")

    # Get the Variance Vector
    vec_variance  <- apply(mat_exprs,1,var)

    # get the symbols
    vec_symbol    <- fData(x)[rownames(mat_exprs),
                              colNameOfGeneSymbol]

    # Create a new data frame to cross reference
    df_xref      <- as.data.frame(cbind(mat_exprs, "variance"= vec_variance))

    # Attach the Symbols and the ProbIDs
    df_xref <- cbind(df_xref,"SYMBOL"=vec_symbol,"probeID"=rownames(mat_exprs))

    # order by variance and find the unique SYMBOL list
    df_xref      <- df_xref[order(df_xref[[method]], decreasing=T),]
    uniqueGeneList  <- unique(df_xref[["SYMBOL"]])

    # get first occurance, (highest variance)
    new_matrix    <- df_xref[match(uniqueGeneList,
                                   table = df_xref[["SYMBOL"]]), ]

    # Update the eset
    x2 <- x[as.vector(new_matrix$probeID),]

    # Return the new eset

    esetToReturn <- x2

  }else if (method == 'logFC') {

    cat("Processing eset using:",method,"\n")

    # order by logFC
    topTable <- topTable[order(topTable[[method]], decreasing = TRUE),]

    # subset the df_xref by the rows in the topTable, make sure that:
    #   nrows of df_xref match nrows of topTable
    df_xref <- df_xref[rownames(topTable),]

    # Check that they match
    if(!assertthat::are_equal(rownames(df_xref),rownames(topTable)) ) {
      stop("Matrix rows do not match topTable rows in method=logFC")
    }

    # Combine
    df_xref <- cbind(df_xref,topTable)
    df_xref[[colNameOfGeneSymbol]] <- df_xref[[colNameOfGeneSymbol]]

    # Order by method

    df_xref      <- df_xref[order(df_xref[[method]], decreasing=T),]
    uniqueGeneList  <- unique(df_xref[[colNameOfGeneSymbol]])


    # get first occurance, (highest variance) and create the new matrix
    new_matrix    <- df_xref[match(uniqueGeneList,
                                   table = df_xref[[colNameOfGeneSymbol]]), ]

    mat_exprs <- new_matrix[,colnames(mat_exprs)]
    rownames(mat_exprs) <- as.vector(new_matrix[[colNameOfGeneSymbol]])

    esetToReturn <- x[as.vector(new_matrix$probeID),]
    exprs(esetToReturn) <- as.matrix(mat_exprs)

    new_fData <- new_matrix[,c("SYMBOL","variance","probeID")]
    rownames(new_fData)  <- new_fData[[colNameOfGeneSymbol]]

    fData(esetToReturn) <- new_fData

    assertthat::are_equal(rownames(exprs(esetToReturn)),
                          rownames(fData(esetToReturn)))
  }

  else if (method == 'mean') {

    cat("Processing eset using:",method,"\n")

    # Get the mean Vector
    vec_mean  <- apply(mat_exprs,1,mean)

    # get the symbols
    vec_symbol    <- fData(x)[rownames(mat_exprs),
                              colNameOfGeneSymbol]

    # Create a new data frame to cross reference
    df_xref      <- as.data.frame(cbind(mat_exprs, "mean"= vec_mean))

    # Attach the Symbols and the ProbIDs
    df_xref <- cbind(df_xref,"SYMBOL"=vec_symbol,"probeID"=rownames(mat_exprs))

    # order by mean and find the unique SYMBOL list
    df_xref      <- df_xref[order(df_xref[[method]], decreasing=T),]
    uniqueGeneList  <- unique(df_xref[["SYMBOL"]])

    # get first occurance, (highest mean)
    new_matrix    <- df_xref[match(uniqueGeneList,
                                   table = df_xref[["SYMBOL"]]), ]

    # Update the eset
    x2 <- x[as.vector(new_matrix$probeID),]

    # Return the new eset

    esetToReturn <- x2
  }

  esetToReturn
}

mapGmt2NamedList <- function(gmtFileName=NULL, type='simple') {
  # This function  takes in a GMT file and creates a named list
  if(type == 'simple') {
    # means name\tgene1\tgene2\t...
    # if default type is simple set to 2
    fieldStart = 2
  } else if (type == 'broad'){
    # name\tURL\t\gene1\t\gene2\
    # if the default type is broad, set to 3
    fieldStart = 3
  } else {
    # if there default is missing
    # Same as broad but 2nd position is not necessarily URL
    fieldStart = 3
  }
  geneSet <- list()
  gmtFile <- readLines(con=gmtFileName)
  gmtFile <- strsplit(gmtFile, '\t')
  for (entry in gmtFile) {
    entry <- unlist(entry)
    geneSet[[unlist(entry[1])[1]]] <- entry[fieldStart:length(entry)]
  }
  rm(gmtFile,entry, gmtFileName,type)
  geneSet
}

CalculateTorque <- function(rho_t, ApcLocs, pow){
  # Takes in rho, the correlation matrix for a network to all genes
  # And DIGElocs, the position of the DIGE targets
  # Calculate the torque of the difference in CDFs
  # First, calculate center of mass, c
  # Then, multiply c*A

  corr_vec <- as.vector(rho_t);
  vals <- sort(corr_vec, decreasing = FALSE)
  pos <- order(corr_vec) ; # default: least to greatest

  samp_vec <- rho_t[,ApcLocs]
  samp_vec <- as.vector(samp_vec)

  ###
  #[C,ia,ib] = intersect(___) also returns index vectors ia and ib using any of the previous syntaxes.

  a <- intersect(vals, samp_vec);
  b <- match(a,vals)
  c <- match(a,samp_vec)


  missed <- matrix(1, nrow = length(vals), ncol = 1);
  hits <- matrix(0, nrow = length(vals), ncol = 1)

  missed[b] <- 0; # Arrange by order of 'vals'
  hits[b] <- 1; #abs(vals(b)).^pow;

  Pmiss = cumsum(missed/sum(missed))
  Phit = cumsum(hits/sum(hits))
  diff_cdf = Phit - Pmiss;


  cut <- min(which(vals>=0)) # Evaluate area on each side of rho=0
  F.neg = diff_cdf[1:cut];
  F.pos = diff_cdf[(cut+1):(length(diff_cdf))];

  c_of_m.neg = sum(F.neg*vals[1:cut])/sum(F.neg);
  c_of_m.pos = sum(F.pos*vals[(cut+1):(length(diff_cdf))])/sum(F.pos);

  # Weight = area = Riemann sum
  weight.neg = sum(F.neg * diff(vals[1:(cut+1)]));
  weight.pos = sum(F.pos * diff(vals[(cut): (length(diff_cdf))]));
  torque = (c_of_m.neg*weight.neg + c_of_m.pos*weight.pos);
  # Should be negative value and negative value, for   greatest bimodality
  return(torque)
}

apcProts_function <- function(txt_file) {
  read.delim(txt_file, sep="\t", header=T)
} #reads protein names text file and transforms it into a data frame

eset.raw_function <- function(Apc_dChip_matrix, Apc_dChip, p_data, Apc_p21_ids){
  if (assertthat::are_equal(colnames(Apc_dChip), rownames(p_data)) &
      assertthat::are_equal(rownames(Apc_dChip), rownames(Apc_p21_ids))){
    eset.raw <-  Biobase::ExpressionSet(assayData = Apc_dChip_matrix,
                                        phenoData = Biobase::AnnotatedDataFrame(data = p_data),
                                        featureData = Biobase::AnnotatedDataFrame(data = Apc_p21_ids))
  }

  else{
    Apc_dChip_matrix <-  Apc_dChip_matrix[rownames(Apc_dChip) %in% rownames(Apc_p21_ids), ]

    Apc_p21_ids <- Apc_p21_ids[rownames(Apc_dChip), ]
    Apc_p21_ids1 <- as.data.frame(Apc_p21_ids)
    colnames(Apc_p21_ids1) <- "SYMBOL"
    rownames(Apc_p21_ids1) <- rownames(Apc_dChip)
    Apc_p21_ids <- Apc_p21_ids1

    eset.raw <-  Biobase::ExpressionSet(assayData = Apc_dChip_matrix,
                                        phenoData = Biobase::AnnotatedDataFrame(data = p_data),
                                        featureData = Biobase::AnnotatedDataFrame(data = Apc_p21_ids))
  }
  eset.raw
}

eset.collapsed_function <- function(eset.raw){
  x <- eset.raw
  exprs(x) <-  limma::normalizeBetweenArrays(object = exprs(x), method = 'quantile' )

  eset.collapsed <- collapseProbes(x,
                                   method = "mean",
                                   topTable = NULL,
                                   colNameOfStat = NULL,
                                   colNameOfGeneSymbol = 'SYMBOL')
  eset.collapsed
}

mRNA_function <- function(eset.collapsed, Apc_p21_ids){
  mRNA <- exprs(eset.collapsed)
  mRNA2 <- merge(x=mRNA, y=Apc_p21_ids, by="row.names")

  row.names(mRNA2)<- mRNA2$SYMBOL
  mRNA2$SYMBOL<-NULL
  mRNA2[,1]<-NULL

  mRNA <- as.matrix(mRNA2)
  mRNA
}

tscore_function <- function(mRNA){

  tscore <- (rowMeans((mRNA[, 1:8]))-rowMeans((mRNA[, 9:16]))) /
    sqrt((apply((mRNA[,1:8]),1,var)/8) + (apply((mRNA[,9:16]),1,var)/8))

  tscore
}

ts_function <- function(tscore){

  ts <- abs(tscore)/max(abs(tscore))
  ts_bin <- matrix(0, nrow = length(ts),ncol =  1)

  x <- sort(ts)
  y <- order(ts)

  cut <- round(.75*length(ts))
  ts_bin[which(ts>=x[cut])] = 1

  ts
}


main_function <- function(APC_RMA_txt, Apc_p21_txt, pdata_txt, Apc_prots_txt, petal_networks_txt){
  Apc_dChip <- assayData_function(APC_RMA_txt)
  Apc_p21_ids <- featureData_function(Apc_p21_txt)
  p_data <- phenoData_function(pdata_txt)
  Apc_dChip_matrix <- data.matrix(Apc_dChip)
  ApcProts_raw <- apcProts_function(Apc_prots_txt)

  eset.raw <- eset.raw_function(Apc_dChip_matrix, Apc_dChip, p_data, Apc_p21_ids)

  eset.collapsed <- eset.collapsed_function(eset.raw)

  mRNA <- mRNA_function(eset.collapsed, Apc_p21_ids)

  #density plot
  limma::plotDensities(eset.collapsed, legend=FALSE, main="All raw data" )

  #boxplot
  graphics::boxplot(exprs(eset.collapsed), las=2 )

  petal_networks_list <- mapGmt2NamedList(petal_networks_txt, type = "broad")


  z <- dplyr::intersect(toupper(row.names(mRNA)), ApcProts_raw$ALDH1L1)
  ApcProts<- which(toupper(row.names(mRNA)) %in% ApcProts_raw$ALDH1L1)
  ApcLocs<-ApcProts; # needs the location of protiens in microarray data
  ApcProts<-z

  i<-1
  NetLocs<-as.list(NULL);
  for(i in 1:length(petal_networks_list)){
    NetLocs[[i]] <- which(toupper(row.names(mRNA)) %in% petal_networks_list[[i]])
  }

  tscore <- tscore_function(mRNA)

  ts <- ts_function(tscore)

  cores = parallel::detectCores()
  cl <- parallel::makeCluster(cores[1]-1) #not to overload your computer
  doParallel::registerDoParallel(cl)

  ########
  rho.all<- list((matrix(0,2,2)));
  rho.t <- list((matrix(0,2,2)));

  torque.net<-rep(-1,24)
  torque.pval<-0
  p<-1;

  p_value <- 0
  bimodality_score <- 0
  petal_names <- 0


  for(p in 1:length(petal_networks_list)){
    rho.all[[p]] =  cor(x = t(mRNA[NetLocs[[p]],]), y = t(mRNA),
                        method = 'pearson',use="pairwise.complete.obs")  ;

    rho.t[[p]] =   rho.all[[p]] *  (ts[NetLocs[[p]] ]) ;  # rho, the correlation matrix for a network to all genes

    torque.net[p] = CalculateTorque(rho.t[[p]], ApcLocs, 0);
    maxperm = 1000;

    D = matrix(0,nrow = maxperm, ncol=1);
    for( m in 1:maxperm){
      r_group=sample(x = 1:ncol(rho.t[[p]]),size= length(ApcLocs),replace=F)
      D[m] = CalculateTorque(rho.t[[p]], r_group, 0);
    }

    torque.pval[p] = length(which(D<torque.net[p]))/length(which(D!=0));
    print(paste(p,torque.net[p],names(petal_networks_list)[p],torque.pval[p]));


    petal_names[p] <- names(petal_networks_list)[p]
    bimodality_score[p] <- torque.net[p]
    p_value[p] <- torque.pval[p]
  }
  bimodality_df <- data.frame(petal_names, bimodality_score, p_value)
  bimodality_df
}


bimodality_df <- main_function("Apc RMA.txt", "Apc p21 - ids.txt", "pdata.txt", "Apc prots.txt", "petal_network.txt")

bimodality_df <- bimodality_df[order(bimodality_df$p_value, decreasing = TRUE), ]
#bimodality_df[order(bimodality_df$bimodality_score, decreasing = TRUE), ]

write.csv(bimodality_df[order(bimodality_df$p_value, decreasing = TRUE), ], file='/Users/alpun/Downloads/bimodality/vignettes/Bimodality_with_variance.csv', row.names = FALSE, quote = FALSE)
