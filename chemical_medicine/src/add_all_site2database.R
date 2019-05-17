# Date : 2018-11-05 
# For : add all.rs.target database 


#----------------------------------- read raw data from xlsx file ---------------------------------------#
library(readxl)
file <- read_xlsx("../多柔比星_数据库_grade_20181106.xlsx",sheet = 3)[,1:10][,c(-2,-6,-8)]
file 


#----------------------------------- get snp infomation from 1000G(Grch37) web site ---------------------------------------#

library(foreach)
library(doParallel) #install.packages("doParallel")
library(magrittr)
library(stringr)
get_items = function(snp){
  library(xml2)
  library(stringr)
  library(rvest) 
  site <- sprintf("http://grch37.ensembl.org/Homo_sapiens/Variation/Explore?db=core;v=%s;vdb=variation;",snp)
  webpage <- read_html(site)
  tmp <- html_nodes(webpage,css="p span")[4]  %>% html_text(trim = T) %>% str_replace_all(pattern = "\\s",replacement = ":") %>% str_split(pattern = "::",simplify = T)
  return(tmp)
} 

#get CPU core NO.
no_cores <- detectCores() 
cl<-makeCluster(no_cores)
registerDoParallel(cl)
s0 <- Sys.time()
b <- foreach(snp = unique(file$variants),.combine = rbind) %dopar% get_items(snp )
t <- Sys.time()-s0
# a <-  system.time(foreach(snp = tt$variants,.combine = rbind) %dopar% get_items(snp ))
stopImplicitCluster()
b
#----------------------------------- merge snp infomation to raw data ---------------------------------------#
colnames(b) <- c("chr","pos","snp","ref","alt")
c <- data.frame(b)
file_ann <- merge(file,c,by.x = "variants",by.y = "snp",all.x = T)

file_ann %$% .[order(genes),] %>% write.table("file_ann_20181106.xls",quote = F,sep = "\t",row.names = F)

file_ann <- read.table("file_ann_20181106.xls",sep = "\t",header = T,stringsAsFactors = F)

head(file_ann)

file_ann$geno <- str_count(file_ann$genotype,file_ann$alt)/2 
n$geno <- str_count(file_ann$genotype,file_ann$alt) 
file_ann <- unique(file_ann)

  
file_ann  %<>% '['(c(2,3,1,5:12,4)) 

a <- c("药物","基因","检测位点","基因型","疗效或毒副作用","等级","chr","pos","ref","alt","geno","ID")
colnames(file_ann) <- a

head(file_ann)
# add to exist file
file_ann %>% write.table("../chemical_medicine_project/target_database_Massarray/ALL.rs.target",sep = "\t",append = T,quote = F,row.names = F,col.names = F)

# add to new   file 
file_ann %>% write.table("../chemical_medicine_project/target_database_Massarray/ALL.rs.target",sep = "\t",append = F,quote = F,row.names = F,col.names = T)