#Prepare data from article supplementary
#Heart

source("D:/Rpackages/AgingBioMolecules/data/Mus musculus/Proteome/addData2list.R")

library(tidyverse)
load("D:/Rpackages/AgingBioMolecules/data/Mus musculus/Proteome/Mice_aging_proteome.Rdata")

#1. Title:Sample multiplexing for targeted pathway proteomics in aging mice
## 10.1073/pnas.1919410117; n=10; months:
raw_data <- read_csv("D:/Rpackages/AgingBioMolecules/data/Mus musculus/Proteome/Article_sup/Yu_et_al_2020_PNAS_Heart_10samples.csv",skip = 1)
raw_data <- separate(raw_data,1,into = c("type","Uniprot.id","Symbol.2"))
fdata <- raw_data[,c(1:6,(ncol(raw_data)-1):(ncol(raw_data)))]
expr <- raw_data[,7:16] %>% as.matrix()
pdata <- data.frame("samples"=colnames(expr),
                    "group"=factor(c(rep("Yung",5),rep("Old",5)),levels = c("Yung","Old")),
                    "month"=c(rep("4",5),rep("20",5)))
FC_temp <- apply(expr, 1, function(x) mean(x[6:10])/mean(x[1:5]))
fdata <- fdata %>% mutate("p.value"=t.pvalue,"q.value"=t.qvalue,"log2FC"=log2(FC_temp),"Significance"=ifelse(p.value<0.05,"sig","nonsig"))
Mice_aging_proteome <- addData2list(list= Mice_aging_proteome, expr = expr, tissue = "Heart",fdata = fdata,pdata = pdata,title = "Yu_2020_PNAS", value.type="Scaled (rowsum=100)",
                                    Title = "Sample multiplexing for targeted pathway proteomics in aging mice",doi = "10.1073/pnas.1919410117",
                                    n = 10,tissue.name = "Heart",age = "4m-20m", quant.method = "TMT",specie.sex="male",specie.strain = "C57BL/6J")


#2. Title: Accurate Quantification of More Than 4000 Mouse Tissue Proteins Reveals Minimal Proteome Changes During Aging
## 10.1074/mcp.M110.004523; n=2; months:5 vs 26
raw_data <- read_csv("D:/Rpackages/AgingBioMolecules/data/Mus musculus/Proteome/Article_sup/Walther_et_al_2011_MOL CELL PROTEOMICS_CerebHeartKidney_samples.csv")
fdata <- raw_data[-c(2,3,4,6)]
colnames(fdata)[c(2,12,13)] <- c("log2FC","Gene.Symbol","Uniprot.id")
Mice_aging_proteome <- addData2list(list= Mice_aging_proteome, tissue = "Heart",fdata = fdata, title = "Walther_2011_MOL CELL PROTEOMICS",
                                    Title = "Accurate Quantification of More Than 4000 Mouse Tissue Proteins Reveals Minimal Proteome Changes During Aging",value.type = "log2FC",
                                    doi = "10.1074/mcp.M110.004523",age = "5m-26m",n=2,specie.sex = "female",specie.strain = "C57BL/6JN")


save(Mice_aging_proteome,file = "D:/Rpackages/AgingBioMolecules/data/Mus musculus/Proteome/Mice_aging_proteome.Rdata")
