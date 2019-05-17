#!/usr/bin/env Rscript
# by: qihao
# date : 2018-07-02,2018-09-11

#check required packages 
options(warn = -1)
package_list <- c("optparse","readxl","magrittr","stringr","tidyr","reshape","dplyr","openxlsx")
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos="http://cran.r-project.org")
    if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
      source("https://bioconductor.org/biocLite.R");biocLite(p)
    }
  }
}

# clean enviroment object
rm(list=ls()) 

# Load essential packages
suppressMessages({library(optparse);library(readxl);library(magrittr);library(stringr);library(tidyr);library(reshape);library(dplyr);library(openxlsx)})

# get options 
if (T){
  option_list <- list(make_option(c("-x", "--xls"), type="character",help = "input file, eg, xx.annotate.filter.xls"), 
                      make_option(c("-t", "--targetdir"), type="character",help = "dir of xx.RS.target file"),
                      make_option(c("-o", "--outdir"), type="character",help = "output dir"))
  opts <- parse_args(OptionParser(option_list=option_list,usage = "Usage: %prog [options]",add_help_option = T))
}

if(F){
  opts = list() 
  opts$xls = "/home/wqh/workdir/chemical_medicine/chemical_medicine_project/test/input_file/Ct1807_1113_20022156.annotate.filter.xls"
  opts$targetdir = "/home/wqh/workdir/chemical_medicine/chemical_medicine_project/target_database_NGS"
  opts$outdir = ""
}



sort_genotype = function(data){
  data[-grep(pattern = "/",data)] %<>% str_split(pattern = "",   simplify = T) %>% data.frame(stringsAsFactors = F) %>% 
    apply( 1, sort) %>% t()%>% data.frame(stringsAsFactors = F) %>%  unite(col = "基因型",c("X1","X2"),sep = "" )  %>% '[['(1)
  data[grep(pattern = "/", data)] %<>% str_split(pattern = "/",simplify = T) %>% data.frame(stringsAsFactors = F) %>% 
    apply( 1, sort) %>% t()%>% data.frame(stringsAsFactors = F) %>%  unite(col = "基因型",c("X1","X2"),sep = "/")  %>% '[['(1)
  return(data)
}

# get sample name and target type 
sample_name = basename(opts$xls) %>% str_split(pattern = "\\.",simplify = T)  %>% '['(1)
type = sample_name %>% str_split(pattern = "_",simplify = T)  %>% '['(2) %>% str_split(pattern = "",simplify = T) %>% '['(3)
if(!is.na(type)){
  target = switch (type,
                   "1" = "ALL.rs.target",
                   "2" = "FA.rs.target",
                   "3" = "JZC.rs.target",
                   "4" = "RXA.rs.target",
                   "5" = "XHD.rs.target",
                   "ALL.rs.target")
}else
  target = "ALL.rs.target"


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
  if ( cut >= 0.2 & cut < 0.8 ) geno = 0.5
  if ( cut <  0.2 ) geno = 0
  Row <- c(tmp[5],tmp[6],tmp[8],tmp[9],geno,cut,tmp[2])
  sample_db <- rbind(sample_db,Row)
}
close(xls)

colnames(sample_db) <- c("chr","pos","ref","alt","geno","H","gene")
sample_db  <- data.frame(sample_db,row.names = NULL,stringsAsFactors = F)
sample_db <- sample_db[sample_db$geno!=0,]
sample_db$genotype = rep(0,nrow(sample_db))


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
sample_db$genotype %<>% sort_genotype()
target_db <- read.table(paste(opts$targetdir,target,sep = "/"),header = T,stringsAsFactors = F,sep = "\t")
target_db$基因型  %<>% sort_genotype()
target_db %<>% unite(col = "id",c("chr","pos","geno"),sep = ":",remove = F)

oo = merge(target_db[,-c(8:11)],sample_db[,c(1,7)],by="id",all.x = T) %$% .[order(药物,基因,检测位点,疗效或毒副作用),] %>% '['(-c(1))  %$% unite(. ,col = "id",c(药物,基因,检测位点,疗效或毒副作用))
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


a = score_out %>% cast(药物 + 基因 + 检测位点 + 基因型 ~ 疗效或毒副作用,fill = NA,value = "等级")
# a[a$药物=="长春瑞滨",]$毒副作用 = NA
A = aggregate(疗效 ~ 药物  ,data=a, FUN=function(n) round(mean(n),digits = 3)*100) %$% .[order(药物),]
B = aggregate(毒副作用 ~ 药物+基因 ,data=a,max)  %$% aggregate(毒副作用 ~ 药物  ,data=., FUN=function(n) round(mean(n),digits = 3)*100 ) %$% .[order(药物),]

Totole = merge(A,B,all.x=T,all.y=T) 
colnames(Totole)[2:3] <- c("疗效评分","毒副作用评分")
Totole$毒副作用评分[is.na(Totole$毒副作用评分)] = "-"

Totole$疗效 = "适中"
Totole$疗效[Totole$疗效评分 <= 40]  = "较差"
Totole$疗效[Totole$疗效评分 >= 60]  = "较好"

Totole$毒副作用 = "适中"
Totole$毒副作用[Totole$毒副作用评分 <= 40 ] = "较弱"
Totole$毒副作用[Totole$毒副作用评分 >= 60 ] = "较强"
Totole$毒副作用[Totole$毒副作用评分 == "-"] = "-"

Totole$综合评价[Totole$疗效 == "较好" & Totole$毒副作用 == "较强" ] = "备选"
Totole$综合评价[Totole$疗效 == "适中" & Totole$毒副作用 == "较强" ] = "慎用"
Totole$综合评价[Totole$疗效 == "较差" & Totole$毒副作用 == "较强" ] = "慎用"

Totole$综合评价[Totole$疗效 == "较好" & Totole$毒副作用 == "适中" ] = "优选"
Totole$综合评价[Totole$疗效 == "适中" & Totole$毒副作用 == "适中" ] = "备选"

Totole$综合评价[Totole$疗效 == "较好" & Totole$毒副作用 == "较弱" ] = "优选"
Totole$综合评价[Totole$疗效 == "适中" & Totole$毒副作用 == "较弱" ] = "备选"
Totole$综合评价[Totole$疗效 == "较差" & Totole$毒副作用 == "较弱" ] = "备选"

Totole$综合评价[is.na(Totole$综合评价)] = "备选"


Grade_score = NULL
for (gene in Totole$药物) {
  out_tmp = subset(a[,1:4],药物 == gene)
  Totole_tmp = subset(Totole,药物 == gene)
  tmp = merge(out_tmp,Totole_tmp,all.x=T)
  tmp[2:nrow(tmp),5:9] = ""
  Grade_score = rbind(Grade_score,tmp)
}
Grade_score = Grade_score[!duplicated(Grade_score[,1:3]),]

# write output into Excel file
colnames(Grade_score)[1] <- "#药物" 
write.table(rbind(colnames(Grade_score)),file=paste0(opts$outdir,"/",sample_name,".CHEMICAL.medicine.xls"),quote=F,sep="\t",row.names=F,col.names = F)
write.table(Grade_score,file=paste0(opts$outdir,"/",sample_name,".CHEMICAL.medicine.xls"),quote=F,sep="\t",col.names = F,append = T,row.names=F)

t_row = target_db[,1:3] %>% unique() %>% nrow()
g_row = Grade_score %>% nrow()
print(t_row == g_row)

rm(list = ls())