#!/usr/bin/env Rscript
rm(list=ls())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages(require(pacman))
if(!suppressPackageStartupMessages(require(rhdf5))){
  source("https://bioconductor.org/biocLite.R")
  biocLite("rhdf5")
}
if(!suppressPackageStartupMessages(require(cmapR)))
  devtools::install_github("cmap/cmapR")

suppressPackageStartupMessages(require(cmapR))
if (!require("pacman")) install.packages ("pacman")
require('pacman')

library(stringr)
library(devtools)
library(S4Vectors)
library(data.table)
library(rmimp)
library(plyr)
p_load(magrittr)
p_load(seqinr)


#Needs to be changed according to the working directory
setwd('/Users/amavi/Downloads/Abhijeet_Mavi/data/prospective_breast_cancer') 


#Input files
ds_path<-'/Users/amavi/Downloads/Abhijeet_Mavi/MIMP/p_original.gct'
ids_path <- "/Users/amavi/Downloads/Abhijeet_Mavi/ids.RData" 
mut_data_path1 <- "~/Downloads/Abhijeet_Mavi/MIMP/mut_original.maf.txt"
mut_data_path <- "/Users/amavi/Downloads/Abhijeet_Mavi/data/prospective_breast_cancer/mutations/germline/CPTAC2_BRCA_prospective.v1.1.germline.variants.122017.maf"
fasta_path <- '/Users/amavi/Downloads/Abhijeet_Mavi/MIMP/RefSeq.20160914_Human_ucsc_hg19_customProDBnr_mito_150contams.fasta'


#Function for creating patient-wise mutation file using terminal
slice.maf <- function(maf, id){
  system(paste("grep '\t", id,"_T\t' ", maf, " > ", id, ".maf", sep=''))
}

load(ids_path) #loading NP data, different for different datatypes


#Preprocessing the fasta file in RefSeq format
fasta <-  read.fasta(fasta_path, seqtype = 'AA')
names(fasta) <- names(fasta) %>% sub('\\..*', '', .) %>% sub('^>','', .)
write.fasta(fasta, names(fasta) , file.out = 'annotated_fasta.fasta')


#Extracting .gct file from directory and the respective metadata and removing Xs in front of patient IDs
my_ds <- parse.gctx(ds_path)
col_metdata <- my_ds@cdesc
rownames(col_metdata ) <- sub('^X','', rownames(col_metdata))
row_metdata <- my_ds@rdesc
mat_metdata <- my_ds@mat
colnames(mat_metdata) <- sub('^X','', colnames(mat_metdata))
cid_metadata <- my_ds@cid %>% sub('^X', '', .)


#Preprocessing the phospho data for equal number of localised and actual sites
row_metdata <- row_metdata[
  which(row_metdata$Best_numActualVMSites_sty==row_metdata$Best_numLocalizedVMsites_sty),
  ]
mat_metdata <- mat_metdata[which(row_metdata$Best_numActualVMSites_sty==row_metdata$Best_numLocalizedVMsites_sty),]
row_metdata$NP_revised<-sub('\\..*','',row_metdata$id)#need to account for the localized data


unique_NP <- unique(row_metdata$NP_revised) #unique NP ids
mut_header <- read.delim(mut_data_path1, nrows = 1, header=F) #mutation file headers


#Algorithm to run MIMP on whole dataset, patient-wise
for (i in cid_metadata ){
  #Creating directory for each patient
  dir.create(i)
  setwd(i)
  
  #Considering only missense mutation and patients with non-zero data
  slice.maf(mut_data_path, i)
  if(file.info(paste(i, 'maf', sep='.'))$size > 0){
    mut_data <- read.delim(paste(i, 'maf', sep='.'), header=TRUE)
    colnames(mut_data) <- mut_header
    
    mut_data <- data.frame(subset(mut_data,Variant_Classification=='Missense_Mutation'))
    mut_file_NM<-sub('\\..*','',mut_data$RefSeq)
    mut_data$NM_revised<-mut_file_NM
    
    keep.idx <- nchar(mut_file_NM) > 0
    
    mut_file_NM <- mut_file_NM[keep.idx] 
    mut_data <- mut_data[keep.idx, ]
    
    mut_data <- data.frame(mut_data, NP_ids=ids[match(mut_file_NM, ids$tx_name), 'pro_name'])
    
    #Creating the dataframe for mutation file in MIMP format from respective columns of mutation file
    df_mut<- data.frame(
      acc=mut_data[ ,'NP_ids'],
      pos=sub('p\\.','',mut_data[ ,'HGVSp_Short'])
    )
    
    #Validating if df_mut is in correct format 
    mut_rows <- which(grepl('^[A-Z][0-9]+[A-Z]$',df_mut$pos) &  
                        (!is.na(df_mut$acc)) & (df_mut$acc!=''))
    df_mut <- df_mut[ mut_rows,]
    
    #Creating the phospho data in MIMP format
    phos_pat_id_match<-which((is.na(mat_metdata[,i])==FALSE))
    
    phos_tmp <- (row_metdata$variableSites[phos_pat_id_match] %>% sub(
      ' $','', .) %>% strsplit(., ' ') %>% lapply(., function(x) x[grep('S|T|Y', x)]))
    
    names(phos_tmp) <- row_metdata[phos_pat_id_match, 'NP_revised']
    
    phos_tmp2 <- lapply(1:length(phos_tmp), function(x) paste(names(phos_tmp)[x], phos_tmp[[x]])) %>% unlist 
    
    df_phospho <- data.frame(
      acc=sub(' .*$', '',phos_tmp2),
      pos=sub('.* ', '',phos_tmp2) %>% sub('^.(.*).$','\\1', .)
    )
    
    #Get the path to example FASTA sequence file 
    data2.file <- "../annotated_fasta.fasta"
    
    # Run rewiring analysis in MIMP
    results_missense_i = mimp(df_mut,data2.file,df_phospho, display.results=TRUE,cores=4)
    write.csv(results_missense_i,file = "result_file.csv")
    
  }
  setwd('..') #resetting the directory to parent directory
}
