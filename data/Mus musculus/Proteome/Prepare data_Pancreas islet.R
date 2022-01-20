#Prepare data from article supplementary
#Pancreas islet

source("D:/Rpackages/AgingBioMolecules/data/Mus musculus/Proteome/addData2list.R")

library(tidyverse)
load("D:/Rpackages/AgingBioMolecules/data/Mus musculus/Proteome/Mice_aging_proteome.Rdata")

#1. Title: Integrated In Vivo Quantitative Proteomics and Nutrient Tracing Reveals Age-Related Metabolic Rewiring of Pancreatic β Cell Function
## 10.1016/j.celrep.2018.11.031; n=8; months: 1 vs 12 months
raw_data <- read_csv("D:/Rpackages/AgingBioMolecules/data/Mus musculus/Proteome/Article_sup/Wortham_et_al_2018_Cell Reports_Pancreas islet_samples.csv")
colnames(raw_data)[c(2,3)] <- c("log2FC","p")
fdata <- raw_data
fdata$Significance <- "sig"
colnames(fdata)[c(1,3)] <- c("Gene.Symbol","p.value")
Mice_aging_proteome <- addData2list(list= Mice_aging_proteome, tissue = "Pancreas islet",fdata = fdata,title = "Wortham_2018_Cell Reports", value.type="stat",
                                    Title = "Integrated In Vivo Quantitative Proteomics and Nutrient Tracing Reveals Age-Related Metabolic Rewiring of Pancreatic β Cell Function",doi = "10.1016/j.celrep.2018.11.031",
                                    age.record = "4-weeks-old vs 1 years",age = "1m-12m",specie.sex = "male-female",specie.strain = "C57BL/6N",n=8)
#could add nonsig protein data later


save(Mice_aging_proteome,file = "D:/Rpackages/AgingBioMolecules/data/Mus musculus/Proteome/Mice_aging_proteome.Rdata")
