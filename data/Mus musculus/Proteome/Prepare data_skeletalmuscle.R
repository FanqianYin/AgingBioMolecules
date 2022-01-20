#Prepare data from article supplementary
#SkeletalMuscle

source("D:/Rpackages/AgingBioMolecules/data/Mus musculus/Proteome/addData2list.R")

library(tidyverse)
load("D:/Rpackages/AgingBioMolecules/data/Mus musculus/Proteome/Mice_aging_proteome.Rdata")

#1. Title:Sample multiplexing for targeted pathway proteomics in aging mice
## 10.1073/pnas.1919410117; n=10; months:
raw_data <- read_csv("D:/Rpackages/AgingBioMolecules/data/Mus musculus/Proteome/Article_sup/Yu_et_al_2020_PNAS_SkeletalMuscle_10samples.csv",skip = 1)
raw_data <- separate(raw_data,1,into = c("type","Uniprot.id","Symbol.2"))
fdata <- raw_data[,c(1:6,(ncol(raw_data)-1):(ncol(raw_data)))]
expr <- raw_data[,7:16] %>% as.matrix()
pdata <- data.frame("samples"=colnames(expr),
                    "group"=factor(c(rep("Yung",5),rep("Old",5)),levels = c("Yung","Old")),
                    "month"=c(rep("4",5),rep("20",5)))
FC_temp <- apply(expr, 1, function(x) mean(x[6:10])/mean(x[1:5]))
fdata <- fdata %>% mutate("p.value"=t.pvalue,"q.value"=t.qvalue,"log2FC"=log2(FC_temp),"Significance"=ifelse(p.value<0.05,"sig","nonsig"))
Mice_aging_proteome <- addData2list(list= Mice_aging_proteome, expr = expr, tissue = "SkeletalMuscle",fdata = fdata,pdata = pdata,title = "Yu_2020_PNAS", value.type="Scaled (rowsum=100)",
                                    Title = "Sample multiplexing for targeted pathway proteomics in aging mice",doi = "10.1073/pnas.1919410117",
                                    n = 10,tissue.name = "Kidney",age = "4m-20m", quant.method = "TMT",specie.sex="male",specie.strain = "C57BL/6J")

#2. Title: Comparative proteomic profiling reveals a role for Cisd2 in skeletal muscle aging
## 10.1111/acel.12705; n=; months: 26 vs 3 months
raw_data <- read_csv("D:/Rpackages/AgingBioMolecules/data/Mus musculus/Proteome/Article_sup/Huang_et_al_2018_Aging Cell_Gastrocnemius_samples.csv")
raw_data$p.value <- 10^(-raw_data$`Significance(-10 log p)`/10)
fdata <- raw_data %>% mutate("log2FC"=log2(FC),"Significance"=ifelse(p.value<0.05,"sig","nonsig"));colnames(fdata)[2] <- "Gene.Symbol"
Mice_aging_proteome <- addData2list(list= Mice_aging_proteome, tissue = "Muscle.Gastrocnemius",tissue.name = "gastrocnemius muscle",fdata = fdata,title = "Huang_2018_Aging Cell", value.type="stat",
                                    Title = "Comparative proteomic profiling reveals a role for Cisd2 in skeletal muscle aging",
                                    age = "3m-26m", doi = "10.1111/acel.12705",specie.strain = "C57BL/6",specie.sex = "male")

#3. Title: Advanced aging causes diaphragm functional abnormalities, global proteome remodeling, and loss of mitochondrial cysteine redox flexibility in mice
## 10.1016/j.exger.2017.12.017; n=12; months: 30 vs 6 months
raw_data <- read_csv("D:/Rpackages/AgingBioMolecules/data/Mus musculus/Proteome/Article_sup/Kelley_et_al_2018_Experimental Gerontology_Diaphragm _12samples.csv",na = "-")
raw_data <- raw_data[-1:-2]
colnames(raw_data)
fdata <- raw_data[,-7:-18];fdata <- fdata[-13];colnames(fdata)[12] <- "log2FC"
fdata <- separate(fdata,1,into = c("Uniprot.id","Symbol.2"),sep = "\\|")
fdata$Gene.Symbol <- get_genesymbol(fdata$Description)
fdata$Significance <- ifelse(fdata$`Significance -10logP`>20,"sig","nonsig")
fdata$p.value <- 10^(-fdata$`Significance -10logP`/10)
expr <- raw_data[,7:18] %>% as.matrix();expr[expr==0] <- NA
pdata <- data.frame("samples"=colnames(expr),
                    "group"=factor(c(rep("Yung",6),rep("Old",6)),levels = c("Yung","Old")),
                    "month"=c(rep("6",6),rep("30",6)))

Mice_aging_proteome <- addData2list(list= Mice_aging_proteome, expr = expr, tissue = "Muscle.Diaphragm",tissue.name = "Inspiratory muscle (diaphragm)",
                                    fdata = fdata,pdata = pdata,title = "Kelley_2018_Experimental Gerontology", value.type="intensity",
                                    Title = "Advanced aging causes diaphragm functional abnormalities, global proteome remodeling, and loss of mitochondrial cysteine redox flexibility in mice",
                                    doi = "10.1016/j.exger.2017.12.017",age = "6m-30m",specie.strain = "C57BL/6",specie.sex = "male",n=12)

#4. Title: Biochemical isolation of myonuclei employed to define changes to the myonuclear proteome that occur with aging
## 10.1111/acel.12604; n=10; months:3 vs 24
#Muscle.nuclei
raw_data <- read_csv("D:/Rpackages/AgingBioMolecules/data/Mus musculus/Proteome/Article_sup/Cutler_et_al_2017_Aging Cell_Muscle nuclei_samples.csv",skip = 1)
colnames(raw_data)[4] <- c("log2FC")
raw_data$p.value <- 10^(-raw_data$`-log 10 t.test`)
fdata <- raw_data
fdata$Significance <- ifelse(raw_data$p.value<0.05,"sig","nonsig");colnames(fdata)[2] <- "Gene.Symbol"
Mice_aging_proteome <- addData2list(list= Mice_aging_proteome, tissue = "Muscle.nuclei",fdata = fdata,title = "Cutler_2017_Aging Cell",
                                    Title = "Biochemical isolation of myonuclei employed to define changes to the myonuclear proteome that occur with aging",
                                    doi = "10.1111/acel.12604",stat.method = "t-test",age = "3m-24m",n=10,specie.strain = "C57BL/6",specie.sex = "male")


save(Mice_aging_proteome,file = "D:/Rpackages/AgingBioMolecules/data/Mus musculus/Proteome/Mice_aging_proteome.Rdata")
