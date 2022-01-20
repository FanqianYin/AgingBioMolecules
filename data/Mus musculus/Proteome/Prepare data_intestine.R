#Prepare data from article supplementary
#Intestine

source("D:/Rpackages/AgingBioMolecules/data/Mus musculus/Proteome/addData2list.R")

library(tidyverse)
load("D:/Rpackages/AgingBioMolecules/data/Mus musculus/Proteome/Mice_aging_proteome.Rdata")

#1. Title:Region-Specific Proteome Changes of the Intestinal Epithelium during Aging and Dietary Restriction
## 10.1016/j.celrep.2020.107565; n=; months:
### Small intestine: 26m vs 3m
raw_data <- read_csv("D:/Rpackages/AgingBioMolecules/data/Mus musculus/Proteome/Article_sup/Gebert_et_al_2020_Cell Reports_Small Intestine_samples.csv",skip = 2)
colnames(raw_data)[c(2,3,10,11)] <- c("Uniprot.id","Gene.Symbol","FC","p.adj")
fdata <- raw_data %>% mutate("p.value"=p.adj,"log2FC"=log2(FC),"Significance"=ifelse(p.value<0.05,"sig","nonsig"))
Mice_aging_proteome <- addData2list(list= Mice_aging_proteome, tissue = "Small intestine",fdata = fdata,title = "Gebert_2020_Cell Reports.1", value.type="stat",
                                    Title = "Region-Specific Proteome Changes of the Intestinal Epithelium during Aging and Dietary Restriction",doi = "10.1016/j.celrep.2020.107565",
                                    age = "3m-26m",tissue.name = "Isolated crypts from small intestine",specie.strain = "C57BL/6J",specie.sex = "male")
### Small intestine: 18m vs 3m
raw_data <- read_csv("D:/Rpackages/AgingBioMolecules/data/Mus musculus/Proteome/Article_sup/Gebert_et_al_2020_Cell Reports_Small Intestine_samples2.csv",skip = 2)
colnames(raw_data)[c(2,3,10,11)] <- c("Uniprot.id","Gene.Symbol","FC","p.adj")
fdata <- raw_data %>% mutate("p.value"=p.adj,"log2FC"=log2(FC),"Significance"=ifelse(p.value<0.05,"sig","nonsig"))
Mice_aging_proteome <- addData2list(list= Mice_aging_proteome, tissue = "Small intestine",fdata = fdata,title = "Gebert_2020_Cell Reports.2", value.type="stat",
                                    Title = "Region-Specific Proteome Changes of the Intestinal Epithelium during Aging and Dietary Restriction",doi = "10.1016/j.celrep.2020.107565",
                                    age = "3m-18m",tissue.name = "Isolated crypts from small intestine",specie.strain = "C57BL/6J",specie.sex = "male")
### Duodenum: 26m vs 3m
raw_data <- read_csv("D:/Rpackages/AgingBioMolecules/data/Mus musculus/Proteome/Article_sup/Gebert_et_al_2020_Cell Reports_Duodenum_samples.csv",skip = 2)
colnames(raw_data)[c(2,3,10,11)] <- c("Uniprot.id","Gene.Symbol","FC","p.adj")
fdata <- raw_data %>% mutate("p.value"=p.adj,"log2FC"=log2(FC),"Significance"=ifelse(p.value<0.05,"sig","nonsig"))
Mice_aging_proteome <- addData2list(list= Mice_aging_proteome, tissue = "Duodenum",fdata = fdata,title = "Gebert_2020_Cell Reports", value.type="stat",
                                    Title = "Region-Specific Proteome Changes of the Intestinal Epithelium during Aging and Dietary Restriction",doi = "10.1016/j.celrep.2020.107565",
                                    age = "3m-26m",specie.strain = "C57BL/6J",specie.sex = "male")
### Ileum: 26m vs 3m
raw_data <- read_csv("D:/Rpackages/AgingBioMolecules/data/Mus musculus/Proteome/Article_sup/Gebert_et_al_2020_Cell Reports_Ileum_samples.csv",skip = 2)
colnames(raw_data)[c(2,3,10,11)] <- c("Uniprot.id","Gene.Symbol","FC","p.adj")
fdata <- raw_data %>% mutate("p.value"=p.adj,"log2FC"=log2(FC),"Significance"=ifelse(p.value<0.05,"sig","nonsig"))
Mice_aging_proteome <- addData2list(list= Mice_aging_proteome, tissue = "Ileum",fdata = fdata,title = "Gebert_2020_Cell Reports", value.type="stat",
                                    Title = "Region-Specific Proteome Changes of the Intestinal Epithelium during Aging and Dietary Restriction",doi = "10.1016/j.celrep.2020.107565",
                                    age = "3m-26m",specie.strain = "C57BL/6J",specie.sex = "male")
### Jejunum: 26m vs 3m
raw_data <- read_csv("D:/Rpackages/AgingBioMolecules/data/Mus musculus/Proteome/Article_sup/Gebert_et_al_2020_Cell Reports_Jejunum_samples.csv",skip = 2)
colnames(raw_data)[c(2,3,10,11)] <- c("Uniprot.id","Gene.Symbol","FC","p.adj")
fdata <- raw_data %>% mutate("p.value"=p.adj,"log2FC"=log2(FC),"Significance"=ifelse(p.value<0.05,"sig","nonsig"))
Mice_aging_proteome <- addData2list(list= Mice_aging_proteome, tissue = "Jejunum",fdata = fdata,title = "Gebert_2020_Cell Reports", value.type="stat",
                                    Title = "Region-Specific Proteome Changes of the Intestinal Epithelium during Aging and Dietary Restriction",doi = "10.1016/j.celrep.2020.107565",
                                    age = "3m-26m",specie.strain = "C57BL/6J",specie.sex = "male")


save(Mice_aging_proteome,file = "D:/Rpackages/AgingBioMolecules/data/Mus musculus/Proteome/Mice_aging_proteome.Rdata")
