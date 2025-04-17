library(readxl)
library(tidyverse)
library(Seurat)
library(org.Hs.eg.db)
library(Biobase)
library(BisqueRNA)
library(ggpubr)
setwd("/data/Warner_Lab/bulk_sting_analysis")
metadata <- read.csv("metadata.csv")

# Gland
sc <- readRDS("/data/Warner_Lab/19_014Warner/Joe/MSG/labels.RDS")
bulk <- read.delim("/data/Warner_Lab/STING_BulkRNA_SjD_HV_nonSjD/MSG_with_corrected_sample_names/report/RawCountFile_filtered.txt",row.names=1,check.names=FALSE)%>%
  dplyr::select(metadata$MSG_Samples)
colnames(bulk) <- gsub("^MSG_", "", colnames(bulk))
ens_to_symbol <- function(matrix) {
  ens <- rownames(matrix)
  symbols <- mapIds(org.Hs.eg.db, keys = ens,
                    column = c('SYMBOL'), keytype = 'ENSEMBL')
  symbols <- symbols[!is.na(symbols)]
  symbols <- symbols[!duplicated(symbols)]
  symbols <- symbols[match(rownames(matrix), names(symbols))]
  symbols <- symbols[!is.na(symbols)]
  matrix <- matrix[names(symbols), ]
  rownames(matrix) <- symbols
  return(matrix)
}
bulk <- ens_to_symbol(bulk)
rm(cell_type_mapping)

# Reformat data to ExpressionSet objects
Idents(sc) <- "celltypes"
sc.data <- base::as.matrix(GetAssayData(sc))
sc.pheno <- base::data.frame(check.names = F, check.rows = F,
                             stringsAsFactors = F, row.names = base::names(Idents(sc)), SubjectName = sc@meta.data$scRNA.ID,
                             cellType = Idents(sc))
sc.meta <- base::data.frame(labelDescription = base::c("SubjectName",
                                                       "cellType"), row.names = base::c("SubjectName", "cellType"))
sc.pdata <- methods::new("AnnotatedDataFrame", data = sc.pheno,
                         varMetadata = sc.meta)
sc.eset <- Biobase::ExpressionSet(assayData = sc.data, phenoData = sc.pdata)
rm(sc.data, sc.meta, sc.pdata, sc.pheno)
bulk.eset <- Biobase::ExpressionSet(assayData = as.matrix(bulk))

# Deconvolution via BisqueRNA
res <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset, markers=NULL, use.overlap=F)
deconv <- as.data.frame(res[["bulk.props"]])
saveRDS(deconv, "deconv_celltypes_gland.RDS")

# PBMC/Blood
sc <- readRDS("/vf/users/Warner_Lab/19_014Warner/Joe/PBMC/labels.RDS")
sc$celltypes <- sc$generaltypes
bulk <- read.delim("/data/Warner_Lab/STING_BulkRNA_SjD_HV_nonSjD/Blood_with_corrected_sample_names/report/RawCountFile_filtered.txt",row.names=1,check.names=FALSE) %>%
  dplyr::select(metadata$Blood_Samples)
colnames(bulk) <- gsub("^Blood_", "", colnames(bulk))
bulk <- ens_to_symbol(bulk)
rm(cell_type_mapping)

# Reformat data to ExpressionSet objects
Idents(sc) <- "celltypes"
sc.data <- base::as.matrix(GetAssayData(sc))
sc.pheno <- base::data.frame(check.names = F, check.rows = F,
                             stringsAsFactors = F, row.names = base::names(Idents(sc)), SubjectName = sc@meta.data$scRNA.ID,
                             cellType = Idents(sc))
sc.meta <- base::data.frame(labelDescription = base::c("SubjectName",
                                                       "cellType"), row.names = base::c("SubjectName", "cellType"))
sc.pdata <- methods::new("AnnotatedDataFrame", data = sc.pheno,
                         varMetadata = sc.meta)
sc.eset <- Biobase::ExpressionSet(assayData = sc.data, phenoData = sc.pdata)
rm(sc.data, sc.meta, sc.pdata, sc.pheno)
bulk.eset <- Biobase::ExpressionSet(assayData = as.matrix(bulk))

# Deconvolution via BisqueRNA
res <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset, markers=NULL, use.overlap=F)
deconv <- as.data.frame(res[["bulk.props"]])
saveRDS(deconv, "deconv_celltypes_blood.RDS")