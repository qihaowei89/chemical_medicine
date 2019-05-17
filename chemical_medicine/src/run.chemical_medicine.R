#!/usr/bin/env Rscript

plist <-  c("tidyr","magrittr","readxl","stringr","reshape")
# library(tidyr)
# library(magrittr)
# library(readxl)
# library(stringr)
# library(reshape)
options(warn = -1)
for (package in plist){
  library(package,character.only = T,quietly = T,warn.conflicts = F)
}
# 
arg <- commandArgs(T)
if(interactive()){
  arg <- c("../NGS_target/化疗.xlsx","化疗")
  pwd <- "/home/wqh/work/chemical_medicine/NGS_target"
}


target <- read.delim("/home/wqh/work/chemical_medicine/chemical_medicine_project/target_database_NGS/ALL.rs.target",stringsAsFactors = F)  %>% 
            unite(col="ID",c("基因","检测位点"),sep = ":",remove = F)
use <- target  %$% ID %>% unique()

pwd <- system("pwd",intern = T) #

 
  
genotype_file <- file.path(pwd,arg[1])

if ((str_detect(arg[1],pattern = ".xlsx$"))) {
  sample <- read_xlsx(genotype_file)[,1:3] 
}

if (str_detect(arg[1],pattern = ".xls$")) {
  sample <- read.delim(genotype_file,stringsAsFactors = F)[,1:3] 
}

sample %<>% unite(col = "IID",c("基因","检测位点"),sep = ":",remove = F)

sample_sub <- subset(sample,IID %in% use) 

sample_sub$检测结果 %<>%  str_replace_all(pattern = "\\(|\\)",replacement = "")  


target_sub <- subset(target,ID %in% sample_sub$IID)

sample_IID <- sample_sub[,-1] %>% unite(col = "IID",c("基因","检测位点","检测结果"),sep = ":") 
target_IID <- target_sub[,-2] %>% unite(col = "IID",c("基因","检测位点","基因型"),sep = ":")

a = merge(sample_IID,target_IID,all.x = T) %$% .[order(IID),]  %$% separate(., col = "IID",into = c("基因","检测位点","基因型"),sep = ":") %>% '['(-7:-11) %>% '['(c(4,1:3,5,6)) %>% 
  cast(药物 + 基因 + 检测位点 + 基因型 ~ 疗效或毒副作用,fill = NA,value = "等级")



A = aggregate(疗效 ~ 药物  ,data=a, FUN=function(n) round(mean(n),digits = 3)*100) %$% .[order(药物),]
B = aggregate(毒副作用 ~ 药物+基因 ,data=a,max)  %$% aggregate(毒副作用 ~ 药物  ,data=., FUN=function(n) round(mean(n),digits = 3)*100 ) %$% .[order(药物),]

Totole = merge(A,B,all.x=T,all.y=T) 

colnames(Totole)[2:3] <- c("疗效评分","毒副作用评分")
Totole$毒副作用评分[is.na(Totole$毒副作用评分)] = "-"
Totole$疗效评分[is.na(Totole$疗效评分)] = "-"

Totole$疗效[Totole$疗效评分!="-"] = "一般"
Totole$疗效[Totole$疗效评分 <= 40]  = "较差"
Totole$疗效[Totole$疗效评分 >= 60]  = "较好"
Totole$疗效[Totole$疗效评分=="-"] = "-"

Totole$毒副作用 = "一般"
Totole$毒副作用[Totole$毒副作用评分 <= 40] = "较弱"
Totole$毒副作用[Totole$毒副作用评分 >= 60] = "较强"
Totole$毒副作用[Totole$毒副作用评分 == "-"] = "-"

Totole$综合评价[Totole$疗效 == "较好" & Totole$毒副作用 == "较强" ]= "备选"
Totole$综合评价[Totole$疗效 == "一般" & Totole$毒副作用 == "较强" ]= "慎用"
Totole$综合评价[Totole$疗效 == "较差" & Totole$毒副作用 == "较强" ]= "慎用"

Totole$综合评价[Totole$疗效 == "较好" & Totole$毒副作用 == "一般" ]= "优选"
Totole$综合评价[Totole$疗效 == "一般" & Totole$毒副作用 == "一般" ]= "备选"
Totole$综合评价[Totole$疗效 == "较差" & Totole$毒副作用 == "一般" ]= "慎用"

Totole$综合评价[Totole$疗效 == "较好" & Totole$毒副作用 == "较弱" ]= "优选"
Totole$综合评价[Totole$疗效 == "一般" & Totole$毒副作用 == "较弱" ]= "备选"
Totole$综合评价[Totole$疗效 == "较差" & Totole$毒副作用 == "较弱" ]= "备选"

Totole$综合评价[is.na(Totole$综合评价)&Totole$毒副作用%in%c("较强","一般")] = "慎用"
Totole$综合评价[is.na(Totole$综合评价)&Totole$疗效%in%c("较好","一般")] = "备选"

Grade_score = NULL
for (gene in Totole$药物) {
  out_tmp = subset(a[,1:4],药物 == gene)
  Totole_tmp = subset(Totole,药物 == gene)
  tmp = merge(out_tmp,Totole_tmp,all.x=T)
  if(nrow(tmp) == 1) Grade_score = rbind(Grade_score,tmp)
  if(nrow(tmp) >  1) {tmp[2:nrow(tmp),5:9] = "";Grade_score = rbind(Grade_score,tmp)}

}

colnames(Grade_score)[1] <- "#药物" 

outpath <- file.path(pwd,Sys.Date() %>% as.character())

if(!dir.exists(outpath)) dir.create(outpath)

name <- arg[2]
if(nrow(sample_sub) == (Grade_score[,2:3] %>% unique() %>% nrow())) {
  cat("finish!\n")
  write.table(Grade_score,file=sprintf("%s/%s.%s.CHEMICAL.medicine.xls",outpath,Sys.Date(),name),quote=F,sep="\t",col.names = T,row.names=F)
}else{
  stop("There were something wrong ! Go check!! and re-run the script.\n")
}




