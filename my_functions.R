my_transpose <- function(df){
  library(data.table,quietly = 1)
  t_df=data.table::transpose(df)
  colnames(t_df) <- rownames(df)
  rownames(t_df) <- colnames(df)
  return(t_df)
}

get_abundance_and_tax_from_table <- function (exp_file, # required
                                              species_are_rows=TRUE,
                                              sep="\t", # extra internal options you can change
                                              row.names=1,
                                              skip=1,
                                              taxa_fields=NULL,
                                              tax_col_name="taxonomy",
                                              tax_sep=";",
                                              NA_option="___",
                                              check.names=FALSE) { # 4 chars or longer "breaks" renamer (will include this as if if were a real name)
  # This function reads a .csv/.tsv speciesXsamples file, where the last "sample" 
  # column (or row) is the taxonomy. Then returns the abundances and taxa data in
  # separate data.frames
  require(gsubfn)
  require(tidyverse)
  if (is.null(taxa_fields)) {taxa_fields = c("Kingdom","Phylum","Class","Order","Family","Genus","Species")}
  
  exp=read.csv(exp_file,sep = sep,skip = skip,row.names=row.names,check.names = check.names)
  if (!species_are_rows) {exp<-my_transpose(exp)}
  
  tax<-exp["taxonomy"]; exp<-exp[1:dim(exp)[2]-1]
  tax <- tax %>% separate("taxonomy",sep = tax_sep,taxa_fields)
  tax[is.na(tax)]<- NA_option # avoid na-related errors
  
  return(list(exp,tax))
}


# 24 sep 2021 gamma_dist_model.R
create_gamma_distr_simuls <- function(data, n) {
  ## data --> real (or not) abundance datasets from where we will obtain
  ##          a mean abundance and its variance for every present OTU
  ## n -----> number of simulations wanted
  data["mean"] = apply(data,MARGIN = 1,FUN = mean) %>% as.data.frame #x_
  data["var"]  = apply(data,MARGIN = 1,FUN = var)  %>% as.data.frame #o2
  beta = data["mean"]**2/data["var"]
  beta[is.na(beta)]<-0
  data["shape"]=beta #a
  data["scale"]=1/(beta/data["mean"]) #s
  data["scale"][is.na(data["scale"])] <-0
  
  
  simul <-apply(data[,c("shape","scale")], MARGIN = 1, FUN=function(row){
    rgamma(n=n, shape=row[1] , scale=row[2])})
}


# 30-sep-2021 2021-09-22_gamma_distr_model n2_correlations.R
# método de Cordero de usar simulaciones con distribuciones gamma para
# dar significancia a correlaciones entre OTUs y hacer buenos grafos
# ++++++++++++++++++++++++++++
# flattenCorrMatrix http://www.sthda.com/english/wiki/correlation-matrix-a-quick-start-guide-to-analyze-format-and-visualize-a-correlation-matrix-using-r-software
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat=NULL) { # incluyo opcion de quitar pmat
  ut <- upper.tri(cormat)
  if (is.null(pmat)) {
    data.frame(
      row = rownames(cormat)[row(cormat)[ut]],
      column = rownames(cormat)[col(cormat)[ut]],
      cor  =cormat[ut]
    )
  } else {
    data.frame(
      row = rownames(cormat)[row(cormat)[ut]],
      column = rownames(cormat)[col(cormat)[ut]],
      cor  =cormat[ut],
      p = pmat[ut]
    )
  }
}


# aplicable a una pmat tmb
flattenCorrMatrixCor <- function(cormat_list) {
  first <- flattenCorrMatrix(cormat_list[[1]])
  next_ones <- lapply(cormat_list[-1], FUN = function(cormat){
    ut <- upper.tri(cormat)
    data.frame(cor=cormat[ut])
  })
  cbind(first,next_ones)%>%setNames(c("row","column",c(names(cormat_list))))
}


unflatten <- function(orig,colnam=c('otu1','otu2'),values="weight") {
  cn1 = colnam[1]
  cn2 = colnam[2]
  
  df=data.frame(row.names = c(orig[[cn1]],orig[[cn2]])%>%unique)
  for (o1 in unique(orig[[cn1]])) {
    for (o2 in unique(orig[[cn2]])) {
      tryCatch(expr={
        df[o1,o2] = orig[orig[[cn1]]==o1 & orig[[cn2]]==o2,values]
        df[o2,o1] = orig[orig[[cn1]]==o1 & orig[[cn2]]==o2,values]
      },error=function(cond){
        # message(cond)
        # message(o1)
        df[o1,o2] = 0
        df[o2,o1] = 0
      }
      )
    }
  }
  df[is.na(df)] <- 0
  df <-df[colnames(df),colnames(df)]
  return(df)
}


pcgcodemaker <- function(tax,f_=FALSE) {
  # Entero == 1; Pseudo == 2; others == 3
  if (f_==TRUE){
    e=" f__Enterobacteriaceae"
    p=" f__Pseudomonadaceae"
  } else {
    e="Enterobacteriaceae"
    p="Pseudomonadaceae"
  }
  entero000 = as.numeric(tax["Family"]!=e)*2+1
  pseudo000 = as.numeric(tax["Family"]==p)
  pcgcode= factor(entero000-pseudo000, ordered=TRUE); rm(entero000); rm(pseudo000)
  return(pcgcode)
}



prettyrenamer <- function(tax, taxtoinclude = 1){
  levels <- c("(K)", "(P)", "(C)", "(O)", "(F)", "(G)", "(S)")
  final_name <- list()
  for (otu in rownames(tax)) {
    lvl <- 0
    tax_list <- c()
    for (tax_level in colnames(tax)){
      taxon <- tax[otu,tax_level]
      if (!is.na(taxon)) {
        if (nchar(taxon) > 4) {
          lvl <- lvl + 1
          tr<-str_split(taxon, "__")[[1]][2]
          tr<-str_trunc(tr, ellipsis = "", width = nchar(taxon)-3)
          tax_list <- c(tax_list, tr)
        }
      }
    }
    if (lvl==7) {tti <- taxtoinclude + 1} else {tti <- taxtoinclude}# if species, I should include the genus too
    final_name[[otu]] <- paste(tax_list[(lvl-(tti-1)):lvl], collapse = "_")
    final_name[[otu]] <- paste(as.character(otu), final_name[[otu]], levels[lvl], sep = "_")
  }
  return(final_name)
}

