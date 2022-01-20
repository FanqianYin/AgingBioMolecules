#Prepare data from article supplementary
#Lung

source("D:/Rpackages/AgingBioMolecules/data/Mus musculus/Proteome/addData2list.R")

library(tidyverse)
load("D:/Rpackages/AgingBioMolecules/data/Mus musculus/Proteome/Mice_aging_proteome.Rdata")

#1. Title:Sample multiplexing for targeted pathway proteomics in aging mice
## 10.1073/pnas.1919410117; n=10; months:
raw_data <- read_csv("D:/Rpackages/AgingBioMolecules/data/Mus musculus/Proteome/Article_sup/Yu_et_al_2020_PNAS_Lung_10samples.csv",skip = 1)
raw_data <- separate(raw_data,1,into = c("type","Uniprot.id","Symbol.2"))
fdata <- raw_data[,c(1:6,(ncol(raw_data)-1):(ncol(raw_data)))]
expr <- raw_data[,7:16] %>% as.matrix()
pdata <- data.frame("samples"=colnames(expr),
                    "group"=factor(c(rep("Yung",5),rep("Old",5)),levels = c("Yung","Old")),
                    "month"=c(rep("4",5),rep("20",5)))
FC_temp <- apply(expr, 1, function(x) mean(x[6:10])/mean(x[1:5]))
fdata <- fdata %>% mutate("p.value"=t.pvalue,"q.value"=t.qvalue,"log2FC"=log2(FC_temp),"Significance"=ifelse(p.value<0.05,"sig","nonsig"))
Mice_aging_proteome <- addData2list(list= Mice_aging_proteome, expr = expr, tissue = "Lung",fdata = fdata,pdata = pdata,title = "Yu_2020_PNAS", value.type="Scaled (rowsum=100)",
                                    Title = "Sample multiplexing for targeted pathway proteomics in aging mice",doi = "10.1073/pnas.1919410117",
                                    n = 10,tissue.name = "Lung",age = "4m-20m", quant.method = "TMT",specie.sex="male",specie.strain = "C57BL/6J")

#2. Title: An atlas of the aging lung mapped by single cell transcriptomics and deep tissue proteomics
## 10.1038/s41467-019-08831-9; n=8; months: 3 vs 24 months
raw_data <- read_csv("D:/Rpackages/AgingBioMolecules/data/Mus musculus/Proteome/Article_sup/Angelidis_et_al_2019_Nat Commun_Lung_8samples.csv")
colnames(raw_data)[c(1,3,4,20,32,33)] <- c("Uniprot.id","Gene.Symbol","log2FC","Significance","-logP","q.value")
fdata <- raw_data[,-5:-12]
fdata$Significance <- ifelse(is.na(fdata$Significance),"nonsig","sig")
fdata <- fdata[-5:-11] %>% mutate("p.value"=10^(-(`-logP`)))
expr <- raw_data[,5:12] %>% as.matrix()
pdata <- data.frame("samples"=colnames(expr),
                    "group"=factor(c(rep("Old",4),rep("Yung",4)),levels = c("Yung","Old")),
                    "month"=c(rep("3",4),rep("24",4)))

Mice_aging_proteome <- addData2list(list= Mice_aging_proteome, expr = expr, tissue = "Lung",fdata = fdata,pdata = pdata,title = "Angelidis_2019_Nat Commun", value.type="exp.mat",
                                    Title = "An atlas of the aging lung mapped by single cell transcriptomics and deep tissue proteomics",doi = "10.1038/s41467-019-08831-9",
                                    stat.method = "t.test",age = "3m-24m",n=8,specie.strain = "C57BL/6")




save(Mice_aging_proteome,file = "D:/Rpackages/AgingBioMolecules/data/Mus musculus/Proteome/Mice_aging_proteome.Rdata")
