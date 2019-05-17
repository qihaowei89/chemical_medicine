#!/usr/bin/env Rscript
# by: qihao
# date : 2018-07-02,2018-08-02

options(warn = -1)
package_list <- c("optparse","readxl","magrittr","stringr")
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos="http://cran.r-project.org")
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

# clean enviroment object
rm(list=ls()) 

# Load essential packages
suppressMessages({
  library(optparse)
  library(readxl)
  library(magrittr)
  library(stringr)
})


if (T){
  option_list <- list(make_option(c("-v", "--vcf"), type="character",help = "input file, eg, xx.vcf"), 
                      make_option(c("-t", "--targetdir"), type="character",help = "dir of xx.RS.target file"),
                      make_option(c("-r","--ref"),type="character",help="sample reorder referece file"),
                      make_option(c("-o", "--outdir"), type="character",help = "output dir"))
  opts <- parse_args(OptionParser(option_list=option_list,usage = "usage: %prog [options]",add_help_option = T))
}

# get sample name and target type 
split = opts$vcf %>% str_split(pattern = "/",simplify = T)  
sample_name = split[length(split)] %>% str_split(pattern = "\\.",simplify = T)  %>% '['(1)
type = sample_name %>% str_split(pattern = "_",simplify = T)  %>% '['(2) %>% str_split(pattern = "",simplify = T) %>% '['(3)
if(!is.na(type)){
  if(type == 1) target = "all.RS.target"
  if(type == 2) target = "lung.RS.target"
  if(type == 3) target = "jzc.RS.target"
  if(type == 4) target = "BC.RS.target"
}else{
  target = "all.RS.target"
}

# get sample genotypes from vcf file 
vcf <- file(opts$vcf) 
open(vcf)
sample_db =　NULL
while(T){
  line <- readLines(vcf,n = 1)
  if (length(line) == 0) break
  if(!grepl(line,pattern = "^#")){
    tmp <- str_split(line,"\t") %>% '[['(1)
    cut <- tmp[10] %>% str_split(":") %>% "[["(1) %>% "["(3) %>% as.numeric()
    if ( cut >= 0.8 ) geno = 1
    if ( cut >= 0.3 & cut < 0.8 ) geno = 0.5
    if ( cut <  0.3 ) geno = 0
    Row <- c(tmp[1],tmp[2],tmp[4],tmp[5],geno,cut)
    sample_db <- rbind(sample_db,Row)
  }
}
close(vcf)
colnames(sample_db) <- c("chr","pos","ref","aln","geno","H")
sample_db  <- data.frame(sample_db,row.names = NULL,stringsAsFactors = F)


# get target type infomations 
target_db  <- read.table(paste(opts$targetdir,target,sep = "/"),header = T,stringsAsFactors = F)
for (line in 1:dim(target_db)[1]) {
  tmp = target_db[line,]
  ref=tmp$ref %>% str_split(pattern = "",simplify = T)
  alt=tmp$alt %>% str_split(pattern = "",simplify = T)
  n = length(ref) - length(alt)
  if(n > 0 ) {
    tmp$pos = tmp$pos+1
    tmp$ref=ref[2:(n+1)] %>% paste0(collapse = "")
    tmp$alt="-" 
    target_db[line,] = tmp
  }
  if( n < 0 ) {
    tmp$pos = tmp$pos+1
    tmp$alt=alt[2:(1-n)] %>% paste0(collapse = "")
    tmp$ref="-" 
    target_db[line,] = tmp
  }
}



CHR_sample = sample_db$chr %>% unique()
CHR_target = target_db$chr %>% unique()
db_vcf     = sapply(CHR_sample, function(i) list(sample_db[which(sample_db$chr %in% i),]))
db_target  = sapply(CHR_target, function(i) list(target_db[which(target_db$chr %in% i),]))
output = NULL
for (i in 1:length(db_target)){
  chr = CHR_target[i]
  tmp_target = db_target[[i]]
  if(chr %in% CHR_sample){
    tmp_target_pos = tmp_target$pos %>% unique()
    ID_vcf = which(CHR_sample %in% chr)
    tmp_vcf_pos = db_vcf[[ID_vcf]]$pos
    for (pos in tmp_target_pos){
      if(pos %in% tmp_vcf_pos){
        tmp_vcf = db_vcf[[ID_vcf]][tmp_vcf_pos %in% pos,]
        tmp = cbind(tmp_target[tmp_target$pos == pos & tmp_target$geno == tmp_vcf$geno,][,-c(7:11)],H=tmp_vcf$H)
        if(!dim(tmp)[1]==0) output = rbind(output,tmp)
      }else{
        tmp =cbind(tmp_target[tmp_target$pos == pos & tmp_target$geno == 0,][,-c(7:11)],H="")
        if(!dim(tmp)[1]==0) output = rbind(output,tmp)
      }
    }
  }else{
    tmp = cbind(tmp_target[tmp_target$geno == 0,][,-c(7:11)],H="")
    if(!dim(tmp)[1]==0) output = rbind(output,tmp)
  }
}

# reorder output file 
ref <- read.table(opts$ref,header = F, stringsAsFactors = F)[,1] %>% unique() 
out=NULL
for(i in ref){
  tmp=output[which(output$药物 %in% i),]
  if(dim(tmp)[1] > 1)
    tmp <- tmp[order(tmp$基因,tmp$检测位点),]
  out = rbind(out,tmp)
}

colnames(out)[5] <- "疗效或毒副作用预测（仅供参考）"
write.table(out,file=paste0(opts$outdir,"/",sample_name,".CHEMICAL.medicine.xls"),quote=F,sep="\t",row.names=F)



