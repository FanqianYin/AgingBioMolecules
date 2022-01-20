
#Prepare data from article supplementary
#HSC

source("D:/Rpackages/AgingBioMolecules/data/Mus musculus/Proteome/addData2list.R")
library(tidyverse)
load("D:/Rpackages/AgingBioMolecules/data/Mus musculus/Proteome/Mice_aging_proteome.Rdata")

#1. Title:Proteomic analysis of young and old mouse hematopoietic stem cells and their progenitors reveals post-transcriptional regulation in stem cells
## 10.7554/eLife.62210; n=10; months: 8–14 weeks vs (24–27 months)
raw_data <- read_csv("D:/Rpackages/AgingBioMolecules/data/Mus musculus/Proteome/Article_sup/Zaro_et_al_2020_elife_HSC_20samples.csv")
fdata <- raw_data[,1];colnames(fdata) <- "Gene.Symbol"
#HSC
expr <- raw_data[,c(2:7,50:55)] %>% as.matrix();expr[expr==0] <- NA
keep.index <- apply(expr, 1, function(x) (sum(is.na(x))/length(x)) <0.5 )
expr <- expr[keep.index,];fdata <- fdata[keep.index,]
expr_imputed <- impute::impute.knn(expr) %>% .$data
rownames(expr_imputed) <- fdata$Gene.Symbol
pdata <- data.frame("samples"=colnames(expr),
                    "group"=factor(c(rep("Young",6),rep("Old",6)),levels = c("Young","Old")),
                    "month"=c(rep("3",6),rep("24",6)))
res <- omicstoolkits::DEA_t.test(Exp = expr_imputed %>% t,groups = pdata$group,log2.trans = T,fearture.name = "Gene.Symbol")
res <- res %>% mutate("p.value"=p,"Significance"=ifelse(p.value<0.05,"sig","nonsig"),.keep="unused")
Mice_aging_proteome <- addData2list(list= Mice_aging_proteome, expr = expr, tissue = "HSC",tissue.name = "Hematopoietic stem cell",
                                    fdata = res,pdata = pdata,title = "Zaro_2020_elife", value.type="Scaled intensity (normalized to sum to 1,000,000)",
                                    Title = "Proteomic analysis of young and old mouse hematopoietic stem cells and their progenitors reveals post-transcriptional regulation in stem cells",
                                    doi = "10.7554/eLife.62210",age.record = "8–14 weeks vs (24–27 months)",age = "3m-25m",n = 12,specie.sex = "male-female",specie.strain = "C57BL/6")
#MPPa
expr <- raw_data[,c(8:13,56:59)] %>% as.matrix();expr[expr==0] <- NA
fdata <- raw_data[,1];colnames(fdata) <- "Gene.Symbol"
keep.index <- apply(expr, 1, function(x) (sum(is.na(x))/length(x)) <0.5 )
expr <- expr[keep.index,];fdata <- fdata[keep.index,]
expr_imputed <- impute::impute.knn(expr) %>% .$data
rownames(expr_imputed) <- fdata$Gene.Symbol
pdata <- data.frame("samples"=colnames(expr),
                    "group"=factor(c(rep("Young",6),rep("Old",4)),levels = c("Young","Old")),
                    "month"=c(rep("3",6),rep("24",4)))
res <- omicstoolkits::DEA_t.test(Exp = expr_imputed %>% t,groups = pdata$group,log2.trans = T,fearture.name = "Gene.Symbol")
res <- res %>% mutate("p.value"=p,"Significance"=ifelse(p.value<0.05,"sig","nonsig"),.keep="unused")
Mice_aging_proteome <- addData2list(list= Mice_aging_proteome, expr = expr, tissue = "MPPa",tissue.name = "Multipotent progenitor a",
                                    fdata = res,pdata = pdata,title = "Zaro_2020_elife", value.type="Scaled intensity (normalized to sum to 1,000,000)",
                                    Title = "Proteomic analysis of young and old mouse hematopoietic stem cells and their progenitors reveals post-transcriptional regulation in stem cells",
                                    doi = "10.7554/eLife.62210",age.record = "8–14 weeks vs (24–27 months)",age = "3m-25m",n = 10,specie.sex = "male-female",specie.strain = "C57BL/6")
#MPPb
expr <- raw_data[,c(14:19,60:63)] %>% as.matrix();expr[expr==0] <- NA
fdata <- raw_data[,1];colnames(fdata) <- "Gene.Symbol"
keep.index <- apply(expr, 1, function(x) (sum(is.na(x))/length(x)) <0.5 )
expr <- expr[keep.index,];fdata <- fdata[keep.index,]
expr_imputed <- impute::impute.knn(expr) %>% .$data
rownames(expr_imputed) <- fdata$Gene.Symbol
pdata <- data.frame("samples"=colnames(expr),
                    "group"=factor(c(rep("Young",6),rep("Old",4)),levels = c("Young","Old")),
                    "month"=c(rep("3",6),rep("24",4)))
res <- omicstoolkits::DEA_t.test(Exp = expr_imputed %>% t,groups = pdata$group,log2.trans = T,fearture.name = "Gene.Symbol")
res <- res %>% mutate("p.value"=p,"Significance"=ifelse(p.value<0.05,"sig","nonsig"),.keep="unused")
Mice_aging_proteome <- addData2list(list= Mice_aging_proteome, expr = expr, tissue = "MPPb",tissue.name = "Multipotent progenitor b",
                                    fdata = res,pdata = pdata,title = "Zaro_2020_elife", value.type="Scaled intensity (normalized to sum to 1,000,000)",
                                    Title = "Proteomic analysis of young and old mouse hematopoietic stem cells and their progenitors reveals post-transcriptional regulation in stem cells",
                                    doi = "10.7554/eLife.62210",age.record = "8–14 weeks vs (24–27 months)",age = "3m-25m",n = 10,specie.sex = "male-female",specie.strain = "C57BL/6")
#MPPc
expr <- raw_data[,c(20:25,64:67)] %>% as.matrix();expr[expr==0] <- NA
fdata <- raw_data[,1];colnames(fdata) <- "Gene.Symbol"
keep.index <- apply(expr, 1, function(x) (sum(is.na(x))/length(x)) <0.5 )
expr <- expr[keep.index,];fdata <- fdata[keep.index,]
expr_imputed <- impute::impute.knn(expr) %>% .$data
rownames(expr_imputed) <- fdata$Gene.Symbol
pdata <- data.frame("samples"=colnames(expr),
                    "group"=factor(c(rep("Young",6),rep("Old",4)),levels = c("Young","Old")),
                    "month"=c(rep("3",6),rep("24",4)))
res <- omicstoolkits::DEA_t.test(Exp = expr_imputed %>% t,groups = pdata$group,log2.trans = T,fearture.name = "Gene.Symbol")
res <- res %>% mutate("p.value"=p,"Significance"=ifelse(p.value<0.05,"sig","nonsig"),.keep="unused")
Mice_aging_proteome <- addData2list(list= Mice_aging_proteome, expr = expr, tissue = "MPPc",tissue.name = "Multipotent progenitor c",
                                    fdata = res,pdata = pdata,title = "Zaro_2020_elife", value.type="Scaled intensity (normalized to sum to 1,000,000)",
                                    Title = "Proteomic analysis of young and old mouse hematopoietic stem cells and their progenitors reveals post-transcriptional regulation in stem cells",
                                    doi = "10.7554/eLife.62210",age.record = "8–14 weeks vs (24–27 months)",age = "3m-25m",n = 10,specie.sex = "male-female",specie.strain = "C57BL/6")


save(Mice_aging_proteome,file = "D:/Rpackages/AgingBioMolecules/data/Mus musculus/Proteome/Mice_aging_proteome.Rdata")
