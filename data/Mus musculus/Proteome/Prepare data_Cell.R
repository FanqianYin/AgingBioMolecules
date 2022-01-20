#Prepare data from article supplementary
#Cell (separate from living mice)

source("D:/Rpackages/AgingBioMolecules/data/Mus musculus/Proteome/addData2list.R")
library(tidyverse)
load("D:/Rpackages/AgingBioMolecules/data/Mus musculus/Proteome/Mice_aging_proteome.Rdata")


#1. Title: Ageing-inducedchangesintheredoxstatusofperipheralmotornerves implyaneffectonredoxsignallingratherthanoxidativedamage
## 10.1016/j.freeradbiomed.2016.02.008; n=8; months:6-8 vs 26-28
#Peripheral motor nerve
raw_data <- read_csv("D:/Rpackages/AgingBioMolecules/data/Mus musculus/Proteome/Article_sup/McDonagh_et_al_2016_Free Radical Biology and Medicine_peripheral motor nerve_8samples.csv",na = "-")
raw_data <- separate(raw_data,1,into = c("Uniprot.id","Symbol.2"),"\\|")
expr <- as.matrix(raw_data[8:15]);expr[expr==0] <- NA
pdata <- data.frame("samples"=colnames(expr),
                    "group"=factor(c(rep("Yung",4),rep("Old",4)),levels = c("Yung","Old")),
                    "month"=c(rep("7",4),rep("27",4)))
fdata <- raw_data[-8:-15]
fdata <- fdata %>% mutate("log2FC"=(log2(`Old Area`/`Adult Area`)),"Gene.Symbol"=get_genesymbol(Description),
                          "p.value"=10^(-(`Significance (-10lgP)`/10)),"Significance"=ifelse(p.value<0.05,"sig","nonsig"))
Mice_aging_proteome <- addData2list(list= Mice_aging_proteome, tissue = "Peripheral motor nerve",fdata = fdata,expr = expr,pdata = pdata,
                                    title = "McDonagh_2016_Free Radical Biology and Medicine",
                                    Title = "Ageing-inducedchangesintheredoxstatusofperipheralmotornerves implyaneffectonredoxsignallingratherthanoxidativedamage",
                                    doi = "10.1016/j.freeradbiomed.2016.02.008",age = "7m-27m",n=8,specie.strain = "C57BL/6",specie.sex = "male")


save(Mice_aging_proteome,file = "D:/Rpackages/AgingBioMolecules/data/Mus musculus/Proteome/Mice_aging_proteome.Rdata")
