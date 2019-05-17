#!/usr/bin/env Rscript
# by: qihao
# date : 2018-08-20,2018-08-23
options(warn = -1)
package_list <- c("tidyr","readxl","magrittr","stringr","reshape","dplyr","optparse",'Biostrings')
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos="http://cran.r-project.org")
    if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
      source("https://bioconductor.org/biocLite.R");biocLite(p)
    }
  }
}

# Load essential packages
sapply(package_list, function(n) require(n,character.only = T,quietly = T,warn.conflicts = FALSE))  



if (!interactive()){
  option_list <- list(make_option(c("-c", "--csv"), type="character",help = "input csv file"), 
                      make_option(c("-s", "--samplelist"), type="character",help = "sample ID list"), 
                      make_option(c("-t", "--targetdir"), type="character",default = "",help = " database dir"),
                      make_option(c("-o", "--outdir"), type="character",help = "output dir"))
  opts <- parse_args(OptionParser(option_list=option_list,usage = "usage: %prog [options]",add_help_option = T))
}

if(interactive()){
  setwd("/home/wqh/work/chemical_medicine/Massarray_target/")
  opts = list()
  opts$csv = "/home/wqh/下载/20190213.csv"
  opts$samplelist = "/home/wqh/20190213sample.list"
  opts$targetdir  = "/home/wqh/work/chemical_medicine/chemical_medicine_project/target_database_Massarray"
  opts$outdir = "20190213/"
}

main <- function(){
  sort_genotype = function(data){
    data[-grep(pattern = "\\.",data)] %<>% str_split(pattern = "",   simplify = T) %>% data.frame(stringsAsFactors = F) %>% 
      apply( 1, sort) %>% t()%>% data.frame(stringsAsFactors = F) %>%  unite(col = "基因型",c("X1","X2"),sep = "" )  %>% '[['(1)
    data[grep(pattern = "\\.", data)] %<>% str_split(pattern = "\\.",simplify = T) %>% data.frame(stringsAsFactors = F) %>% 
      apply( 1, sort) %>% t()%>% data.frame(stringsAsFactors = F) %>%  unite(col = "基因型",c("X1","X2"),sep = ".")  %>% '[['(1)
    
    return(data)
  }
  
  filt_genotype = function(data){
    data$基因型 %<>% gsub(data$基因型,pattern = "(^[ATCG]$)",replacement = "\\1\\1",perl = T) %>% 
      gsub(pattern = "^DEL$",replacement = "del.del",perl = T)  %>% 
      gsub(pattern = "DEL",replacement = "del",perl = T) %>%
      gsub(pattern = "^([ATCG]{3,}$)",replacement = "\\1\\.\\1",perl = T) %>% 
      gsub(pattern = "-",replacement = "del") %>%  
      gsub(pattern = "\\(TA\\)",replacement = "TA") %>% 
      gsub(pattern = "\\/",replacement = "\\.") %>%
      gsub(pattern = "^(TA\\d$)",replacement = "\\1.\\1") %>%
      sort_genotype()
    if("ref" %in% colnames(data)) data$ref %<>% gsub(pattern = "-",replacement = "del")
    return(data)
  }
  
  # read Massrray raw genotype data
  raw_data = read.csv(opts$csv,skip = 2,stringsAsFactors = F,na.strings = "")[,c(1,3,4)] 
  colnames(raw_data) =c("Sample.ID","基因型","ID")
  sample_list = read.csv(opts$samplelist,header = F,stringsAsFactors = F,col.names = c("Sample.ID","type"))
  
  if(!dir.exists(opts$outdir)) system(sprintf("mkdir %s",opts$outdir))
  
  for (l in 1:nrow(sample_list)) { 
    type = sample_list[l,]$type
    id = sample_list[l,]$Sample.ID
    database = read.table(sprintf("%s/%s.rs.target",opts$targetdir,type),sep = "\t",header = T,stringsAsFactors = F)[,-c(7:8,10)] %>% filt_genotype() 
    sample_data = raw_data[raw_data$Sample.ID == id & !is.na(raw_data$基因型),]
    # sample_data = raw_data[raw_data$Sample.ID == id,]
    # sample_data_raw = raw_data[raw_data$Sample.ID == id & !is.na(raw_data$基因型),]
    both_id = intersect(sample_data$ID,unique(database$ID)) %>% sort()
    id_2_rs_list <- database[,c(3,9)] %>% unique()
    
    sample_data = sample_data[sample_data$ID %in% both_id,]  %>% filt_genotype() %$% .[order(ID),]
    # sample_data_raw = sample_data_raw[sample_data_raw$ID %in% both_id,] %$% .[order(ID),]
    sample_data  %<>% merge( y = id_2_rs_list,all.x = T)
    
    # sample_data_raw  %<>% merge( y = id_2_rs_list,all.x = T)
    # c_bind <- data.frame(sample_data, sample_data_raw )
    
    
    b = unite(database[,-6],col = "u_id",c("药物","基因","检测位点","疗效或毒副作用","ID"),sep = ":") %>% unique() %$% .[order(u_id),]
    m_sites = data.frame(x = cbind((b$u_id %>% table())[(b$u_id %>% table()) != 3] %>% names()),stringsAsFactors = F) %>% 
      separate(col="x",c("药物","基因","检测位点","疗效或毒副作用","ID"),sep = ":") %>% '['(c(3,5)) %>% unique()
    
    if(length(m_sites$ID) != 0) {
      m_sites_database = database[database$检测位点==m_sites$检测位点,] %>% unite(col = "IID",c("ID","基因型"),sep = ":")
      score_out_2 = m_sites_database[m_sites_database$IID %in% unite(data = sample_data[sample_data$ID %in% m_sites$ID,],col = "IID",c("ID","基因型"),sep = ":")$IID,] %>% separate(col = "IID",c("ID","基因型"),sep = ":") %>% '['(-c(8,9))
      s_sites_database = database[database$检测位点!=m_sites$检测位点,]
    }else{
      score_out_2 = NULL
      s_sites_database = database
    }
    
    s_sites_database$geno[s_sites_database$geno == 1  ] = 2
    s_sites_database$geno[s_sites_database$geno == 0.5] = 1
    s_sites_database$geno[s_sites_database$geno == 0  ] = 0
    
    sample_data_1 = sample_data[!(sample_data$ID %in% m_sites$ID),]
    sample_data_1 = merge(sample_data_1,(s_sites_database %$% .[order(ID),] %>%'['(c(7,9)) %>% unique()),all.x = T)
    
    raw_1 = sample_data_1[ grepl(pattern = "\\.",x = sample_data_1$基因型),] %>% separate(col = "基因型",c("alt1","alt2"),sep = "\\.")
    raw_2 = sample_data_1[!grepl(pattern = "\\.",x = sample_data_1$基因型),] 
    raw_2 = data.frame(raw_2,(raw_2$基因型 %>% str_split(pattern = "",simplify = T)),stringsAsFactors = F)[,c(1,2,6,7,4,5)]
    colnames(raw_2)[3:4] = c("alt1","alt2")
    sample_data_1 = rbind(raw_1,raw_2)
    
    sample_data_1$ref1[str_length(sample_data_1$ref)==1] = sapply(sample_data_1$ref[str_length(sample_data_1$ref)==1],function(n) DNAString(n) %>% reverseComplement() %>% as.character()) 
    sample_data_1$ref1[str_length(sample_data_1$ref)!=1] = sample_data_1$ref[str_length(sample_data_1$ref)!=1]
    
    
    sample_data_1$geno[sample_data_1$alt1 != sample_data_1$alt2] = 1
    sample_data_1$geno[(sample_data_1$alt1 == sample_data_1$alt2) & (sample_data_1$alt1 == sample_data_1$ref |sample_data_1$alt1 == sample_data_1$ref1)] = 0
    sample_data_1$geno[(sample_data_1$alt1 == sample_data_1$alt2) & (sample_data_1$alt1 != sample_data_1$ref & sample_data_1$alt1 != sample_data_1$ref1)] = 2
    
    
    s_sites_database %<>% unite(col = "IID",c(ID,geno),sep = ":")
    sample_data_1 %<>% unite(col = "IID",c(ID,geno),sep = ":")
    # score_out_1  = merge(sample_data_1,s_sites_database,all.x=T) %>% separate(col = "IID",into = c("ID","x"),sep = ":",remove = T) %>% '['(c(-2,-4:-8)) %>% '['(c(3,4,2,1,5,6,7))
    
    # n <- merge(sample_data_1,s_sites_database,all.x=T) %>% separate(col = "IID",into = c("ID","x"),sep = ":",remove = T) %>% '['(c(-2,-4,-5,-8,-11))  
    conbine <- function(n){
      need <- apply(n[,3:4] %>% sapply(str_length) != 1 ,1,any)
      a <- n[need, ] %>%  unite(col="基因型",c("alt1","alt2"),sep = "/") %>% '['(c(4,5,2,1,3,6,7))
      b <- n[!need,] %>%  unite(col="基因型",c("alt1","alt2"),sep = "") %>% '['(c(4,5,2,1,3,6,7))
      return(rbind(a,b))
    }
    
    
    # score_out_1  = merge(sample_data_1,s_sites_database,all.x=T) %>% separate(col = "IID",into = c("ID","x"),sep = ":",remove = T) %>% '['(c(-2,-4,-5,-8,-11)) %>% unite(col="基因型",c("alt1","alt2"),sep = "") %>% '['(c(4,5,2,1,3,6,7))
    score_out_1  = merge(sample_data_1,s_sites_database,all.x=T) %>% separate(col = "IID",into = c("ID","x"),sep = ":",remove = T) %>% '['(c(-2,-4,-5,-8,-11))  %>% conbine()
    
    ##   
    score_out = rbind(score_out_1,score_out_2) %$% .[order(药物,基因,检测位点,ID,疗效或毒副作用),]
    score_out$基因[is.na(score_out$基因)] = "intergenic" 
    
    # ref = score_out$药物 %>% unique()
    
    a = score_out %>% cast(药物 + 基因 + 检测位点 + ID + 基因型 ~ 疗效或毒副作用,fill = NA,value = "等级")
    # a[a$药物=="长春瑞滨",]$毒副作用 = NA
    A = aggregate(疗效 ~ 药物  ,data=a, FUN=function(n) round(mean(n),digits = 3)*100) %$% .[order(药物),]
    B = aggregate(毒副作用 ~ 药物+基因 ,data=a,max)  %$% aggregate(毒副作用 ~ 药物  ,data=., FUN=function(n) round(mean(n),digits = 3)*100 ) %$% .[order(药物),]
    Totole = merge(A,B,all.x=T,all.y=T) 
    colnames(Totole)[2:3] <- c("疗效评分","毒副作用评分")
    Totole$毒副作用评分[is.na(Totole$毒副作用评分)] = "-"
    Totole$疗效评分[is.na(Totole$疗效评分)] = "-"
    
    Totole$疗效 = "一般"
    Totole$疗效[Totole$疗效评分 <= 40]  = "较差"
    Totole$疗效[Totole$疗效评分 >= 60]  = "较好"
    Totole$疗效[Totole$疗效评分 =="-"]  = "-"
    
    Totole$毒副作用 = "一般"
    Totole$毒副作用[Totole$毒副作用评分 <= 40] = "较弱"
    Totole$毒副作用[Totole$毒副作用评分 >= 60] = "较强"
    Totole$毒副作用[Totole$毒副作用评分 =="-"] = "-"
    
    Totole$综合评价[Totole$疗效 == "较好" & Totole$毒副作用 == "较强" ]= "备选"
    Totole$综合评价[Totole$疗效 == "一般" & Totole$毒副作用 == "较强" ]= "慎用"
    Totole$综合评价[Totole$疗效 == "较差" & Totole$毒副作用 == "较强" ]= "慎用"
    
    Totole$综合评价[Totole$疗效 == "较好" & Totole$毒副作用 == "一般" ]= "优选"
    Totole$综合评价[Totole$疗效 == "一般" & Totole$毒副作用 == "一般" ]= "备选"
    Totole$综合评价[Totole$疗效 == "较差" & Totole$毒副作用 == "一般" ]= "慎用"
    
    Totole$综合评价[Totole$疗效 == "较好" & Totole$毒副作用 == "较弱" ]= "优选"
    Totole$综合评价[Totole$疗效 == "一般" & Totole$毒副作用 == "较弱" ]= "备选"
    Totole$综合评价[Totole$疗效 == "较差" & Totole$毒副作用 == "较弱" ]= "备选"
    
    # Totole$综合评价[is.na(Totole$综合评价)] = "备选"
    Totole$综合评价[Totole$疗效%in%c("一般","较好") & Totole$毒副作用=="-"] = "备选"
    Totole$综合评价[Totole$疗效=="-"] = "慎用"
    ##,"巯嘌呤","柔红霉素","替加氟"
    sort = switch (type,
      "ALL" = c("奥沙利铂", "卡铂", "奈达铂", "顺铂", "环磷酰胺", "氟尿嘧啶", "吉西他滨", "甲氨蝶呤", "甲酰四氢叶酸", "卡培他滨", "培美曲塞", "多西他赛", "伊立替康", "依托泊苷", "长春瑞滨", "长春新碱", "紫杉醇", "阿那曲唑", "来曲唑", "他莫昔芬", "表柔比星", "多柔比星", "替加氟"),
      "FA"  = c("奥沙利铂", "卡铂", "奈达铂", "顺铂", "环磷酰胺", "吉西他滨", "培美曲塞", "多西他赛", "伊立替康", "依托泊苷", "长春瑞滨", "紫杉醇"),
      "SGA" = c("奥沙利铂", "卡铂", "奈达铂", "顺铂", "氟尿嘧啶", "吉西他滨", "卡培他滨", "多西他赛", "伊立替康", "紫杉醇", "表柔比星"),
      "WA"  = c("奥沙利铂", "卡铂", "奈达铂", "顺铂", "环磷酰胺", "氟尿嘧啶", "卡培他滨", "多西他赛", "伊立替康", "紫杉醇", "表柔比星"),
      "JZC" = c("奥沙利铂", "氟尿嘧啶", "甲酰四氢叶酸", "卡培他滨", "伊立替康"),
      "GJA" = c("卡铂", "顺铂", "氟尿嘧啶", "吉西他滨", "培美曲塞", "多西他赛", "伊立替康", "长春瑞滨", "长春新碱", "紫杉醇"),
      "RXA" = c("卡铂", "顺铂", "环磷酰胺", "氟尿嘧啶", "吉西他滨", "甲氨蝶呤", "多西他赛", "紫杉醇", "阿那曲唑", "来曲唑", "他莫昔芬", "表柔比星")
    )
    
    
  
    Totole = Totole[match(sort,Totole$药物),]
    
    
    Grade_score = NULL
    for (gene in Totole$药物) {
      out_tmp = subset(a[,1:5],药物 == gene)
      Totole_tmp = subset(Totole,药物 == gene)
      tmp = merge(out_tmp,Totole_tmp,all.x=T)
      if(nrow(tmp) == 1) Grade_score = rbind(Grade_score,tmp)
      if(nrow(tmp) >  1) {tmp[2:nrow(tmp),6:10] = "";Grade_score = rbind(Grade_score,tmp)}
    }
    colnames(Grade_score)[1] <- "#药物" 
    Grade_score %<>% '['(-4)      # For check results
    Grade_score$基因型 %<>% str_replace_all(pattern = "del",replacement = "Del") %>% str_replace_all(pattern = "TGGCGCGTCC",replacement = "Ins")  %>% str_replace_all(pattern = "ATTTGTTCATGCCT",replacement = "Ins")
    # print(Grade_score)
    
    write.table(rbind(colnames(Grade_score)),file=paste0(opts$outdir,"/",id,"_",type,".CHEMICAL.medicine.xls"),quote=F,sep="\t",row.names=F,col.names = F)
    write.table(Grade_score,file=paste0(opts$outdir,"/",id,"_",type,".CHEMICAL.medicine.xls"),quote=F,sep="\t",col.names = F,append = T,row.names=F)
    cat(sprintf("%s : finished!\n",id))
  }
}

#--------------------------------- main -----------------------------#
main()



