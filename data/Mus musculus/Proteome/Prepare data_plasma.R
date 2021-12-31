#Prepare data from article supplementary

library(tidyverse)
Mice_aging_proteome <- list()

##plasma
Mice_aging_proteome[["plasma"]] <- list()
###1. Title:Plasma proteomic profiling of young and old mice reveals cadherin-13 prevents age-related bone loss
  ###doi:10.18632/aging.103184; n=12; months:
Mice_aging_proteome[["plasma"]][[1]] <- list()
expr_data <- read_csv("D:/Rpackages/AgingBioMolecules/data/Mus musculus/Proteome/Article_sup/Yang_et_al_2020_aging_plasma_12samples.csv")
expr_data <- expr_data[1:3280,]
expr_data$`Gene names`[expr_data$`Gene names`=="0"] <- NA
expr_data1 <- expr_data[,c(5:10,17:22)] %>% as.matrix();rownames(expr_data1) <- expr_data$Accession
expr_data2 <- expr_data[,c(11:16,23:28)] %>% as.matrix();rownames(expr_data2) <- expr_data$Accession

fdata <- expr_data[,1:4]
sdata1 <- expr_data[,c(1:4,29:31)] %>% mutate("Significance"=ifelse(`G_value OG`>= 3.841,"sig","nonsig"))
sdata2 <- expr_data[,c(1:4,32:34)] %>% mutate("Significance"=ifelse(`G_value HP`>= 3.841 & `SAM (D) HP` >= 1.96,"sig","nonsig"))

pdata <- data.frame("samples"=colnames(expr_data1),
                    "group"=factor(c(rep("Yung",6),rep("Old",6)),levels = c("Yung","Old")),
                    "month"=c(rep("2",6),rep("22",6)))

Mice_aging_proteome[["plasma"]][[1]][["expr1"]] <- expr_data1;Mice_aging_proteome[["plasma"]][[1]][["expr2"]] <- expr_data2
Mice_aging_proteome[["plasma"]][[1]][["expr_zscore"]] <- scale((expr_data1/rep(apply(expr_data1, 2, sum,na.rm=T),each=nrow(expr_data1))) %>% t) %>% t
#Mice_aging_proteome[["plasma"]][[1]][["expr_zscore"]] <- scale(log2(expr_data1+1) %>% t) %>% t
Mice_aging_proteome[["plasma"]][[1]][["fdata"]] <- fdata
Mice_aging_proteome[["plasma"]][[1]][["sdata1"]] <- sdata1;Mice_aging_proteome[["plasma"]][[1]][["sdata2"]] <- sdata2
Mice_aging_proteome[["plasma"]][[1]][["pdata"]] <- pdata
Mice_aging_proteome[["plasma"]][[1]][["title"]] <- "Yang_2020_Aging"


save(Mice_aging_proteome,file = "D:/Rpackages/AgingBioMolecules/data/Mus musculus/Proteome/Mice_aging_proteome.Rdata")
