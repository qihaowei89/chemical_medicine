---
output:
  word_document: default
  html_document: default
---
# 1. R packages 依赖

`"optparse","readxl","magrittr","stringr","tidyr","reshape","dplyr","openxlsx","progress"`


通过如下命令自动安装依赖的R packages(首次运行时)

```
Rscript <Massarray.chemical_medicine.R> -h 

```
# 2. 文件夹介绍

## 2.1 target_database_Massarray/ 存放质谱化疗的数据库

- ALL.rs.target **综合数据库**
- FA.rs.target  **肺癌数据库**
- GJA.rs.target **宫颈癌数据库**
- JZC.rs.target **结直肠癌数据库**
- RXA.rs.target **乳腺癌数据库**
- SGA.rs.target **食管癌数据库**
- WA.rs.target  **胃癌数据库**


## 2.2 target_database_NGS/ 存放NGS化疗的数据库

- ALL.rs.target **综合数据库**
- FA.rs.target  **肺癌数据库**
- JZC.rs.target **结直肠癌数据库**
- RXA.rs.target **乳腺癌数据库**
- XHD.rs.target **消化道癌数据库**

## 2.2.1 样本命名与数据库对应的规则：   
- **"1" = "ALL.rs.target"**   
- **"2" = "FA.rs.target"**   
- **"3" = "JZC.rs.target"**   
- **"4" = "RXA.rs.target"**   
- **"5" = "XHD.rs.target"**

例如：   
sample1 = 'Ct19020074_1011_10032024'   
sample2 = 'Ct19020076_1031_BB19020006' 

sample1: __1011__中第三位数字为1，则对应数据库为ALL.rs.target   
sample2: __1031__中第三位数字为3，则对应数据库为JZC.rs.target

# 3. src/ 存放分析脚本文件
## **3.1 ctDNA标签法流程脚本为: `chemical_medicine.R`**


## 使用方法： 
**`Rscript chemical_medicine.R`    `[options]`**

```    
    Options:
        -x XLS, --xls=XLS
                input file, eg, xx.annotate.filter.xls
                
        -t TARGETDIR, --targetdir=TARGETDIR
                dir of xx.RS.target file

        -o OUTDIR, --outdir=OUTDIR
                output dir

        -h, --help
                Show this help message and exit
```             
**示例：**
    
`Rscript src/chemical_medicine.R \
    -x test/input_file/Ct1807_1113_20022156.annotate.filter.xls \
    -t  target_database_NGS \
    -o test/output`
    
    
    
    
## **3.2 ctDNA旧流程(刘蓓)脚本为: `vcf.chemical_medicine.R`**    

## 使用方法
**`Rscript vcf.chemical_medicine.R`    `[options]`**
```    
    Options:
        -v VCF, --vcf=VCF
                input file, eg, xx.vcf
                
        -t TARGETDIR, --targetdir=TARGETDIR
                dir of xx.RS.target file

        -o OUTDIR, --outdir=OUTDIR
                output dir

        -h, --help
                Show this help message and exit
```         
**示例：**

`Rscript src/vcf.chemical_medicine.R \
-v test/input_file/Ct18090036_1011_FA18080044.MuTect2.vcf \
-t target_database_NGS \
-o test/output`
    
    
# **3.3 Massarray数据分析脚本为: `Massarray.chemical_medicine.R`**

## 使用方法
**`Rscipt Massarray.chemical_medicine.R`    `[options]`**
```
    Options:
        -c CSV, --csv=CSV
                input csv file

        -s SAMPLELIST, --samplelist=SAMPLELIST
                sample ID list

        -t TARGETDIR, --targetdir=TARGETDIR
                 database dir

        -o OUTDIR, --outdir=OUTDIR
                output dir

        -h, --help
                Show this help message and exit
```

**示例：**

`Rscript src/Massarray.chemical_medicine.R \
-c test/input_file/20180905_recall.csv \
-s test/input_file/sample.list \  
-t target_database_Massarray \   
-o test/output`
    

### **sample.list** 格式示例：

```
#ID,type
10024562,SGA
10024562,WA
10024562,FA
10024562,JZC
```


# **3.4 泛生子外包化疗数据分析脚本: `run.chemical_medicine.R`**   
## 使用方法
**`Rscipt run.chemical_medicine.R`    `<input_file.xlsx>` `<output_file_name>`**


**示例：** 

`Rscipt run.chemical_medicine.R 杨好人.xlsx 杨好人`



# 结果文件示例

```
#药物         基因     检测位点    基因型   疗效评分  毒副作用评分  疗效  毒副作用  综合评价
奥沙利铂      ERCC1    rs11615     GG       56.5      50        一般  一般      备选
奥沙利铂      ERCC2    rs13181     TT
奥沙利铂      GSTP1    rs1695      AA
奥沙利铂      MTHFR    rs1801131   TT
奥沙利铂      MTHFR    rs1801133   AG
奥沙利铂      XRCC1    rs25487     CC
表柔比星      ERCC1    rs11615     GG       46.7      50        一般  一般      备选
表柔比星      GSTP1    rs1695      AA
表柔比星      XRCC1    rs25487     CC
多西他赛      ABCB1    rs1045642   AG       46.7      55        一般  一般      备选
多西他赛      ERCC1    rs11615     GG
多西他赛      ERCC2    rs13181     TT
氟尿嘧啶      DPYD     rs3918290   CC       56.7      62.2      一般  较强      慎用
氟尿嘧啶      DPYD     rs55886062  AA
氟尿嘧啶      DPYD     rs67376798  TT
氟尿嘧啶      GSTP1    rs1695      AA
氟尿嘧啶      TP53     rs1042522   CC
氟尿嘧啶      UMPS     rs1801019   GG
...


```


