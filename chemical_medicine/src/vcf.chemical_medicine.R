#!/usr/bin/env Rscript
# by: qihao
# date : 2018-08-16,2018-09-11,2018-10-9

#check required packages 
options(warn = -1)
rm(list=ls()) 
package_list <- c("optparse","readxl","magrittr","stringr","tidyr","reshape","dplyr","openxlsx","progress")
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos="http://cran.r-project.org")
    if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
      source("https://bioconductor.org/biocLite.R")
      biocLite(p)
    }
  }
}

# clean enviroment object


# Load essential packages
sapply(package_list, function(n) require(n,character.only = T,quietly = T,warn.conflicts = FALSE)) 

if (T){
  option_list <- list(make_option(c("-v", "--vcf"), type="character",help = "input file, eg, xx.vcf"), 
                      make_option(c("-t", "--targetdir"), type="character",help = "dir of xx.RS.target file"),
                      make_option(c("-o", "--outdir"), type="character",help = "output dir"))
  opts <- parse_args(OptionParser(option_list=option_list,usage = "usage: %prog [options]",add_help_option = T))
}

if (F){
  opts = list()
  opts$vcf = "/home/wqh/work/chemical_medicine/chemical_medicine_project/test/vcf/Ct18100347_1011_10024073.MuTect2.vcf"
  opts$targetdir = "/home/wqh/work/chemical_medicine/chemical_medicine_project/target_database_NGS/"
  opts$outdir = "/home/wqh/work/chemical_medicine/chemical_medicine_project/test"
}





sort_genotype = function(data){
  data[-grep(pattern = "/",data)] %<>% str_split(pattern = "",   simplify = T) %>% data.frame(stringsAsFactors = F) %>% 
    apply( 1, sort) %>% t()%>% data.frame(stringsAsFactors = F) %>%  unite(col = "基因型",c("X1","X2"),sep = "" )  %>% '[['(1)
  data[grep(pattern = "/", data)] %<>% str_split(pattern = "/",simplify = T) %>% data.frame(stringsAsFactors = F) %>% 
    apply( 1, sort) %>% t()%>% data.frame(stringsAsFactors = F) %>%  unite(col = "基因型",c("X1","X2"),sep = "/")  %>% '[['(1)
  return(data)
}

main <- function(){
  pb <- progress_bar$new(
    format = "  :what | :bar | :percent , runtime: :elapsed",
    clear = FALSE, total = 4, width = 60)
  # get sample name and target type 
  # step1
  pb$tick(tokens = list(what = "step1"))
  {
    # split = opts$vcf %>% str_split(pattern = "/",simplify = T)  
    sample_name = basename(opts$vcf) %>% str_split(pattern = "\\.",simplify = T)  %>% '['(1)
    type = sample_name %>% str_split(pattern = "_",simplify = T)  %>% '['(2) %>% str_split(pattern = "",simplify = T) %>% '['(3)
    if(!is.na(type)){
      target = switch (type,
                       "1" = "ALL.rs.target",
                       "2" = "FA.rs.target" ,
                       "3" = "JZC.rs.target",
                       "4" = "RXA.rs.target",
                       "5" = "XHD.rs.target",
                       "ALL.rs.target")
    }else{
      target = "ALL.rs.target"
    }
  }
  
  # target_db  <- read.table(paste(opts$targetdir,target,sep = "/"),header = T,stringsAsFactors = F,sep = "\t")
  # 
  # tmp_target.bed <- target_db[,c(3,7:8)] %>% unique() %>% .[,c(2,3,1)]  %$% .[order(chr,pos),] 
  # tmp_target.bed  %>% write.table(file.path(opts$outdir,"tmp_target.bed"),sep = "\t",quote = F,row.names = F,col.names = F)
  # 
  # 
  # 
  # target_region_depth <- system(sprintf("/home/wqh/miniconda3/bin/samtools depth -a -b %s -d 100000 -q 30 -Q 30 %s.target.realigned.removemultiple.bam",
  #                                       file.path(opts$outdir,"tmp_target.bed"), paste(dirname(opts$vcf),"../bam",sample_name,sep = "/")),intern = T)
  # 
  # target_region_depth %<>% cbind %>% str_split(pattern = "\t",simplify = T) %>% as.data.frame()
  # 
  # cbind(tmp_target.bed,depth =target_region_depth[,3]) %>% as.data.frame() %$% subset(.,depth !=0)
  # out_of_pannle_site <- cbind(tmp_target.bed,depth =target_region_depth[,3]) %>% as.data.frame() %$% subset(.,depth ==0,select = 检测位点)  
  # 
  # save(out_of_pannle_site,file = "out_of_pannle_site.RData")
  
  
  # get sample genotypes from vcf file 
  # step2
  pb$tick(tokens = list(what = "step2"))
  {
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
        if ( cut >= 0.2 & cut < 0.8 ) geno = 0.5
        if ( cut <  0.2 ) geno = 0
        Row <- c(tmp[1],tmp[2],tmp[4],tmp[5],geno,cut)
        sample_db <- rbind(sample_db,Row)
      }
    }
    close(vcf)
    colnames(sample_db) <- c("chr","pos","ref","alt","geno","H")
    sample_db  <- data.frame(sample_db,row.names = NULL,stringsAsFactors = F)
    sample_db <- sample_db[sample_db$geno!=0,]
    
    sample_db$genotype = rep(0,nrow(sample_db))
    for (line in 1:dim(sample_db)[1]) {
      tmp = sample_db[line,]
      tmp$pos %<>% as.numeric()
      n = str_length(tmp$ref) - str_length(tmp$alt)
      if(n != 0) {
        if(n > 0 ) { #缺失
          tmp$pos = tmp$pos+1
          tmp$ref = sub(pattern = tmp$alt,replacement = "",x=tmp$ref)
          tmp$alt = "-"
        }
        if( n < 0 ) { #插入
          tmp$pos = tmp$pos+1
          tmp$alt = sub(pattern = tmp$ref,replacement = "",x=tmp$alt)
          tmp$ref = "-"
        }
        if(tmp$geno == 0 &  (tmp$ref %in% c("A","T","C","G"))) tmp$genotype = paste(tmp$ref,tmp$ref,sep = "" )
        if(tmp$geno == 0 & !(tmp$ref %in% c("A","T","C","G"))) tmp$genotype = paste(tmp$ref,tmp$ref,sep = "/")
        
        if(tmp$geno == 0.5) tmp$genotype = paste(tmp$ref,tmp$alt,sep = "/")
        
        if(tmp$geno == 1 &  (tmp$alt %in% c("A","T","C","G"))) tmp$genotype = paste(tmp$alt,tmp$alt,sep = "" )
        if(tmp$geno == 1 & !(tmp$alt %in% c("A","T","C","G"))) tmp$genotype = paste(tmp$alt,tmp$alt,sep = "/")
        
        
        sample_db[line,] = tmp
      }else{
        if(tmp$geno == 0 & tmp$ref == "-" ) tmp$genotype = paste(tmp$ref,tmp$ref,sep = "/")
        if(tmp$geno == 0 & tmp$ref != "-" ) tmp$genotype = paste(tmp$ref,tmp$ref,sep = "" )
        if(tmp$geno == 0.5 &  (tmp$ref != "-" & tmp$ref != "-")) tmp$genotype = paste(tmp$ref,tmp$alt,sep = "")
        if(tmp$geno == 0.5 &  (tmp$ref == "-" | tmp$ref == "-")) tmp$genotype = paste(tmp$ref,tmp$alt,sep = "")
        if(tmp$geno == 1 & tmp$alt == "-" ) tmp$genotype = paste(tmp$alt,tmp$alt,sep = "/")
        if(tmp$geno == 1 & tmp$alt != "-" ) tmp$genotype = paste(tmp$alt,tmp$alt,sep = "" )
        
        sample_db[line,] = tmp
      }
    }
  }
  
  # step3
  pb$tick(tokens = list(what = "step3"))
  {
    sample_db %<>% unite(col = "id", c("chr","pos","geno"),remove = F,sep = ":")
    target_db  <- read.table(paste(opts$targetdir,target,sep = "/"),header = T,stringsAsFactors = F,sep = "\t")
    target_db  %<>% unite(col = "id",c("chr","pos","geno"),sep = ":",remove = F)
    target_db$基因型  %<>% sort_genotype()
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
    A = aggregate(疗效 ~ 药物  ,data=a, FUN=function(n) round(mean(n),digits = 3)*100) %$% .[order(药物),]
    B = aggregate(毒副作用 ~ 药物+基因 ,data=a,max)  %$% aggregate(毒副作用 ~ 药物  ,data=., FUN=function(n) round(mean(n),digits = 3)*100 ) %$% .[order(药物),]
    
    Totole = merge(A,B,all.x=T,all.y=T) 
    colnames(Totole)[2:3] <- c("疗效评分","毒副作用评分")
    Totole$毒副作用评分[is.na(Totole$毒副作用评分)] = "-"
    
    Totole$疗效 = "适中"
    Totole$疗效[Totole$疗效评分 <= 40]  = "较差"
    Totole$疗效[Totole$疗效评分 >= 60]  = "较好"
    Totole$疗效[Totole$疗效评分=="-"] = "-"
    
    Totole$毒副作用 = "适中"
    Totole$毒副作用[Totole$毒副作用评分 <= 40] = "较弱"
    Totole$毒副作用[Totole$毒副作用评分 >= 60] = "较强"
    Totole$毒副作用[Totole$毒副作用评分 == "-"] = "-"
    
    Totole$综合评价[Totole$疗效 == "较好" & Totole$毒副作用 == "较强" ]= "备选"
    Totole$综合评价[Totole$疗效 == "适中" & Totole$毒副作用 == "较强" ]= "慎用"
    Totole$综合评价[Totole$疗效 == "较差" & Totole$毒副作用 == "较强" ]= "慎用"
    
    Totole$综合评价[Totole$疗效 == "较好" & Totole$毒副作用 == "适中" ]= "优选"
    Totole$综合评价[Totole$疗效 == "适中" & Totole$毒副作用 == "适中" ]= "备选"
    Totole$综合评价[Totole$疗效 == "较差" & Totole$毒副作用 == "适中" ]= "慎用"
    
    Totole$综合评价[Totole$疗效 == "较好" & Totole$毒副作用 == "较弱" ]= "优选"
    Totole$综合评价[Totole$疗效 == "适中" & Totole$毒副作用 == "较弱" ]= "备选"
    Totole$综合评价[Totole$疗效 == "较差" & Totole$毒副作用 == "较弱" ]= "备选"
    
    Totole$综合评价[is.na(Totole$综合评价)] = "备选"
    

    Grade_score = NULL
    for (gene in Totole$药物) {
      out_tmp = subset(a[,1:4],药物 == gene)
      Totole_tmp = subset(Totole,药物 == gene)
      tmp = merge(out_tmp,Totole_tmp,all.x=T)
      if(nrow(tmp) == 1) Grade_score = rbind(Grade_score,tmp)
      if(nrow(tmp) >  1) {tmp[2:nrow(tmp),5:9] = "";Grade_score = rbind(Grade_score,tmp)}
    }
    Grade_score = Grade_score[!duplicated(Grade_score[,1:3]),]
  } 
  
  # write output into Excel file
  #step4
  pb$tick(tokens = list(what = "Complate!"))
  {
    colnames(Grade_score)[1] <- "#药物" 
    write.table(rbind(colnames(Grade_score)),file=paste0(opts$outdir,"/",sample_name,".CHEMICAL.medicine.xls"),quote=F,sep="\t",row.names=F,col.names = F)
    write.table(Grade_score,file=paste0(opts$outdir,"/",sample_name,".CHEMICAL.medicine.xls"),quote=F,sep="\t",col.names = F,append = T,row.names=F)
    t_row = target_db[,1:3] %>% unique() %>% nrow()
    g_row = Grade_score %>% nrow()
    print(t_row == g_row)
  }
  rm(list = ls())
}

#----------------------------- main ----------------------------# 
main()

