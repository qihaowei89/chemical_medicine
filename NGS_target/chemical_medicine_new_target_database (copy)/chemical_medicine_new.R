#!/usr/bin/env Rscript
# by: qihao
# date : 2018-07-02,
 


options(warn = -1)
package_list <- c("optparse","readxl","magrittr","stringr","tidyr","reshape","dplyr")
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos="http://cran.r-project.org")
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

# clean enviroment object
rm(list=ls()) 

# Load essential packages
suppressMessages({library(optparse);library(readxl);library(magrittr);library(stringr);library(tidyr);library(reshape);library(dplyr)})


if (T){
  option_list <- list(make_option(c("-x", "--xls"), type="character",help = "input file, eg, xx.annotate.filter.xls"), 
                      make_option(c("-t", "--targetdir"), type="character",help = "dir of xx.RS.target file"),
                      make_option(c("-o", "--outdir"), type="character",help = "output dir"))
  opts <- parse_args(OptionParser(option_list=option_list,usage = "usage: %prog [options]",add_help_option = T))
}



# get sample name and target type 
split = opts$xls %>% str_split(pattern = "/") %>% '[['(1)  
sample_name = split[length(split)] %>% str_split(pattern = "\\.",simplify = T)  %>% '['(1)
type = sample_name %>% str_split(pattern = "_",simplify = T)  %>% '['(2) %>% str_split(pattern = "",simplify = T) %>% '['(3)
if(!is.na(type)){
  if(type == 1) target = "ALL.rs.target"
  if(type == 2) target = "FA.rs.target"
  if(type == 3) target = "JZC.rs.target"
  if(type == 4) target = "RXA.rs.target"
  if(type == 5) target = "XHD.rs.target"
}else{
  target = "ALL.rs.target"
}


# get sample genotypes from xls file 
xls <- file(opts$xls) 
open(xls)
sample_db =　NULL
h=readLines(xls,n = 1) 

while(T){
  line <- readLines(xls,n = 1)
  line 
  if (length(line) == 0) break
  tmp <- str_split(line,"\t") %>% '[['(1)
  cut <- tmp[18]  %>% as.numeric()
  if ( cut >= 0.8 ) geno = 1
  if ( cut >= 0.3 & cut < 0.8 ) geno = 0.5
  if ( cut <  0.3 ) geno = 0
  Row <- c(tmp[5],tmp[6],tmp[8],tmp[9],geno,cut,tmp[2])
  sample_db <- rbind(sample_db,Row)
}
close(xls)

colnames(sample_db) <- c("chr","pos","ref","alt","geno","H","gene")
sample_db  <- data.frame(sample_db,row.names = NULL,stringsAsFactors = F)
sample_db$genotype = NA

for (line in 1:dim(sample_db)[1]) {
  tmp = sample_db[line,]
  tmp$pos %<>% as.numeric()
  n = (str_length(tmp$ref) == 1 & str_length(tmp$alt) == 1 & !any(grepl(pattern = "-",x = c(tmp$ref,tmp$alt))))
  if(!n) {
    if(tmp$geno == 0  ) tmp$genotype = paste(tmp$ref,tmp$ref,sep = "/")
    if(tmp$geno == 0.5) tmp$genotype = paste(tmp$ref,tmp$alt,sep = "/")
    if(tmp$geno == 1  ) tmp$genotype = paste(tmp$alt,tmp$alt,sep = "/")
    sample_db[line,] = tmp
  }else{
    if(tmp$geno == 0  ) tmp$genotype = paste(tmp$ref,tmp$ref,sep = "")
    if(tmp$geno == 0.5) tmp$genotype = paste(tmp$ref,tmp$alt,sep = "")
    if(tmp$geno == 1  ) tmp$genotype = paste(tmp$alt,tmp$alt,sep = "")
    sample_db[line,] = tmp
  }
}


sample_db %<>% unite(col = "id", c("chr","pos","geno"),remove = F,sep = ":")
target_db <- read.table(paste(opts$targetdir,target,sep = "/"),header = T,stringsAsFactors = F,sep = "\t")
target_db %<>% unite(col = "id",c("chr","pos","geno"),sep = ":",remove = F)

oo = merge(target_db[,-c(8:11)],sample_db[,c(1,7)],all.x = T) %$% .[order(药物,基因,检测位点,疗效或毒副作用),] %>% '['(-c(1))  %$% unite(. ,col = "id",c(药物,基因,检测位点,疗效或毒副作用))
score_out = NULL
for (ID in unique(oo$id)) {
  condition = !(subset(oo, id == ID)$H %>% is.na())
  if (any(condition)) {
    score_tmp = subset(oo, id == ID)[condition,][,-c(4,5)] %$% separate(., col = "id",into = c("药物","基因","检测位点","疗效或毒副作用"),sep = "_")
  }else{
    score_tmp = subset(oo, id == ID & geno == 0)[,-c(4,5)] %$% separate(., col = "id",into = c("药物","基因","检测位点","疗效或毒副作用"),sep = "_")
  }
  score_out = rbind(score_out,score_tmp)
}


ref = score_out$药物 %>% unique()
Grade_score=NULL
for(i in ref){
  tmp = score_out[which(score_out$药物 %in% i),] %>% cast(药物 + 基因 + 检测位点 + 基因型 ~ 疗效或毒副作用,fill = NA)
  
  if(is.null(try(tmp$毒副作用))) tmp$毒副作用 = NA 
  if(is.null(try(tmp$疗效))) tmp$疗效 = NA 
  
  # tmp$毒副作用 = ifelse(is.null(try(tmp$毒副作用)), rep(0,length(tmp$药物)), tmp$毒副作用)
  # tmp$疗效 = ifelse(is.null(try(tmp$疗效)), rep(0,length(tmp$药物)), tmp$疗效)
  
  tmp_sum = summarise(tmp, 疗效评分 = mean(疗效,na.rm=T)*100,毒副作用评分 = mean(毒副作用,na.rm=T)*100) %>% round(digits = 1)
  
  if(tmp_sum$疗效评分 > 40  & tmp_sum$疗效评分 < 60 ) tmp_sum$疗效 = "一般"
  if(tmp_sum$疗效评分 <= 40 ) tmp_sum$疗效 = "较差"
  if(tmp_sum$疗效评分 >= 60 ) tmp_sum$疗效 = "较好"
  
  if(tmp_sum$毒副作用评分 > 40 & tmp_sum$毒副作用评分< 60 ) tmp_sum$毒副作用 = "一般"
  if(tmp_sum$毒副作用评分<= 40 ) tmp_sum$毒副作用 = "较弱"
  if(tmp_sum$毒副作用评分>= 60 ) tmp_sum$毒副作用 = "较强"
  
  if(tmp_sum$疗效 == "较好" & tmp_sum$毒副作用 == "较强") tmp_sum$综合评价 = "备用"
  if(tmp_sum$疗效 == "一般" & tmp_sum$毒副作用 == "较强") tmp_sum$综合评价 = "慎用"
  if(tmp_sum$疗效 == "较差" & tmp_sum$毒副作用 == "较强") tmp_sum$综合评价 = "慎用"
  
  if(tmp_sum$疗效 == "较好" & tmp_sum$毒副作用 == "一般") tmp_sum$综合评价 = "优选"
  if(tmp_sum$疗效 == "一般" & tmp_sum$毒副作用 == "一般") tmp_sum$综合评价 = "备选"
  if(tmp_sum$疗效 == "较差" & tmp_sum$毒副作用 == "一般") tmp_sum$综合评价 = "慎用"
  
  if(tmp_sum$疗效 == "较好" & tmp_sum$毒副作用 == "较弱") tmp_sum$综合评价 = "优选"
  if(tmp_sum$疗效 == "一般" & tmp_sum$毒副作用 == "较弱") tmp_sum$综合评价 = "备选"
  if(tmp_sum$疗效 == "较差" & tmp_sum$毒副作用 == "较弱") tmp_sum$综合评价 = "备选"
  
  row_num = dim(tmp)[1] -1
  Grade_score_tmp = cbind(tmp[,-c(5,6)],rbind(tmp_sum,data.frame(疗效评分 = rep("",row_num),毒副作用评分 = rep("",row_num),疗效 = rep("",row_num),毒副作用 = rep("",row_num),综合评价 = rep("",row_num)))) 
  Grade_score = rbind(Grade_score,Grade_score_tmp)
}

colnames(Grade_score)[1] <- "#药物" 
write.table(rbind(colnames(Grade_score)),file=paste0(opts$outdir,"/",sample_name,".CHEMICAL.medicine.xls"),quote=F,sep="\t",row.names=F,col.names = F)
write.table(Grade_score,file=paste0(opts$outdir,"/",sample_name,".CHEMICAL.medicine.xls"),quote=F,sep="\t",col.names = F,append = T,row.names=F)

