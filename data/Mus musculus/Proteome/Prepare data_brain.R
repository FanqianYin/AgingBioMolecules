#Prepare data from article supplementary
#Brain

source("D:/Rpackages/AgingBioMolecules/data/Mus musculus/Proteome/addData2list.R")
library(tidyverse)
load("D:/Rpackages/AgingBioMolecules/data/Mus musculus/Proteome/Mice_aging_proteome.Rdata")

#1. Title:Sample multiplexing for targeted pathway proteomics in aging mice
## 10.1073/pnas.1919410117; n=10; months:
raw_data <- read_csv("D:/Rpackages/AgingBioMolecules/data/Mus musculus/Proteome/Article_sup/Yu_et_al_2020_PNAS_Brain_10samples.csv",skip = 1)
raw_data <- separate(raw_data,1,into = c("type","Uniprot.id","Symbol.2"))
fdata <- raw_data[,c(1:6,(ncol(raw_data)-1):(ncol(raw_data)))]
expr <- raw_data[,7:16] %>% as.matrix()
pdata <- data.frame("samples"=colnames(expr),
                    "group"=factor(c(rep("Young",5),rep("Old",5)),levels = c("Young","Old")),
                    "month"=c(rep("4",5),rep("20",5)))
FC_temp <- apply(expr, 1, function(x) mean(x[6:10])/mean(x[1:5]))
fdata <- fdata %>% mutate("p.value"=t.pvalue,"q.value"=t.qvalue,"log2FC"=log2(FC_temp),"Significance"=ifelse(p.value<0.05,"sig","nonsig"))
Mice_aging_proteome <- addData2list(list= Mice_aging_proteome, expr = expr, tissue = "Brain",fdata = fdata,pdata = pdata,title = "Yu_2020_PNAS", value.type="Scaled (rowsum=100)",
             Title = "Sample multiplexing for targeted pathway proteomics in aging mice",doi = "10.1073/pnas.1919410117",
             n = 10,tissue.name = "Whole Brain",age = "4m-20m", quant.method = "TMT",specie.sex="male",specie.strain = "C57BL/6J")

#2. Title:Quantitative Proteomics Reveals Significant Differences between Mouse Brain Formations in Expression of Proteins Involved in Neuronal Plasticity during Aging
## 10.3390/cells10082021; n=NA; months:1 vs 22
raw_data <- read_csv("D:/Rpackages/AgingBioMolecules/data/Mus musculus/Proteome/Article_sup/Drulis-Fajdasz_et_al_2021_Cells_brain_samples.csv")
#Hippocampus
fdata <- raw_data[,1:11]
colnames(fdata) <- c("Multi.id","Uniprot.id","Protein names","Gene.Symbol","function_annotation","Average.Old(pmol/mg)","Old.SD","Average.Young(pmol/mg)","Young.SD","p.value","FC")
fdata <- fdata %>% mutate("Significance"=ifelse(p.value<0.05,"sig","nonsig"),"log2FC"=log2(FC))
Mice_aging_proteome <- addData2list(list= Mice_aging_proteome, tissue = "Hippocampus",fdata = fdata,title = "Drulis-Fajdasz_2021_Cells", value.type="mean and SD (Conc. (pmol/mg))",
                                    Title = "Quantitative Proteomics Reveals Significant Differences between Mouse Brain Formations in Expression of Proteins Involved in Neuronal Plasticity during Aging",
                                    doi = "10.3390/cells10082021",n = 10,tissue.name = "Hippocampus",age = "1m-22m",specie.strain = "C57BL/10J",specie.sex = "female")
#Cortex
fdata <- raw_data[,c(1:5,12:17)]
colnames(fdata) <- c("Multi.id","Uniprot.id","Protein names","Gene.Symbol","function_annotation","Average.Old(pmol/mg)","Old.SD","Average.Young(pmol/mg)","Young.SD","p.value","FC")
fdata <- fdata %>% mutate("Significance"=ifelse(p.value<0.05,"sig","nonsig"),"log2FC"=log2(FC))
Mice_aging_proteome <- addData2list(list= Mice_aging_proteome, tissue = "Cortex",fdata = fdata,title = "Drulis-Fajdasz_2021_Cells", value.type="mean and SD (Conc. (pmol/mg))",
                                    Title = "Quantitative Proteomics Reveals Significant Differences between Mouse Brain Formations in Expression of Proteins Involved in Neuronal Plasticity during Aging",
                                    doi = "10.3390/cells10082021",n = 10,tissue.name = "Cortex",age = "1m-22m",specie.strain = "C57BL/10J",specie.sex = "female")
#Cerebellum
fdata <- raw_data[,c(1:5,18:23)]
colnames(fdata) <- c("Multi.id","Uniprot.id","Protein names","Gene.Symbol","function_annotation","Average.Old(pmol/mg)","Old.SD","Average.Young(pmol/mg)","Young.SD","p.value","FC")
fdata <- fdata %>% mutate("Significance"=ifelse(p.value<0.05,"sig","nonsig"),"log2FC"=log2(FC))
Mice_aging_proteome <- addData2list(list= Mice_aging_proteome, tissue = "Cerebellum",fdata = fdata,title = "Drulis-Fajdasz_2021_Cells", value.type="mean and SD (Conc. (pmol/mg))",
                                    Title = "Quantitative Proteomics Reveals Significant Differences between Mouse Brain Formations in Expression of Proteins Involved in Neuronal Plasticity during Aging",
                                    doi = "10.3390/cells10082021",n = 10,tissue.name = "Cerebellum",age = "1m-22m",specie.strain = "C57BL/10J",specie.sex = "female")

#3. Title: Proteomic Profile of Mouse Brain Aging Contributions to Mitochondrial Dysfunction, DNA Oxidative Damage, Loss of Neurotrophic Factor, and Synaptic and Ribosomal Proteins
## 10.1155/2020/5408452; n=20; months:4 vs 16
#Hippocampus
raw_data <- read_csv("D:/Rpackages/AgingBioMolecules/data/Mus musculus/Proteome/Article_sup/Li_et_al_2020_Oxidative Medicine and Cellular Longevity_Hippocampus_samples.csv")
raw_data$Accession <- str_replace_all(raw_data$Accession,"tr\\|","")
raw_data <- raw_data %>% separate(1,into = c("Uniprot.id","Uniprot.id2"),sep = "\\|")
raw_data$Gene.Symbol <- lapply(str_split(raw_data$Description,"="),function(x) x[3]) %>% unlist %>% str_replace_all(" PE","")
raw_data <- raw_data %>% mutate("FC"=get("Ratio Aged"),log2FC=log2(FC),"Significance-p.value"=Significance,"Significance"=ifelse(Significance>=5,"sig","nonsig"))
fdata <- raw_data
Mice_aging_proteome <- addData2list(list= Mice_aging_proteome, tissue = "Hippocampus",fdata = fdata,title = "Li_2020_OXID MED CELL LONGEV", value.type="stat",
                                    Title = "Proteomic Profile of Mouse Brain Aging Contributions to Mitochondrial Dysfunction, DNA Oxidative Damage, Loss of Neurotrophic Factor, and Synaptic and Ribosomal Proteins",
                                    doi = "10.1155/2020/5408452",stat.method = "Peaks Q algorithm",age = "4m-16m",n = 10,specie.strain = "B6129SF2/J",specie.sex = "male")

#Cortex
raw_data <- read_csv("D:/Rpackages/AgingBioMolecules/data/Mus musculus/Proteome/Article_sup/Li_et_al_2020_Oxidative Medicine and Cellular Longevity_Cortex_samples.csv")
raw_data$Accession <- str_replace_all(raw_data$Accession,"tr\\|","")
raw_data <- raw_data %>% separate(1,into = c("Uniprot.id","Uniprot.id2"),sep = "\\|")
raw_data$Gene.Symbol <- lapply(str_split(raw_data$Description,"="),function(x) x[3]) %>% unlist %>% str_replace_all(" PE","")
raw_data <- raw_data %>% mutate("FC"=get("Ratio Aged"),log2FC=log2(FC),"Significance-p.value"=Significance,"Significance"=ifelse(Significance>=5,"sig","nonsig"))

fdata <- raw_data
Mice_aging_proteome <- addData2list(list= Mice_aging_proteome, tissue = "Cortex",fdata = fdata,title = "Li_2020_OXID MED CELL LONGEV", value.type="mean",
                                    Title = "Proteomic Profile of Mouse Brain Aging Contributions to Mitochondrial Dysfunction, DNA Oxidative Damage, Loss of Neurotrophic Factor, and Synaptic and Ribosomal Proteins",
                                    doi = "10.1155/2020/5408452",stat.method = "Peaks Q algorithm",age = "4m-16m",n=12,specie.strain = "B6129SF2/J",specie.sex = "male")

#4. Title: Proteomic analysis of aged microglia: shifts in transcription, bioenergetics, and nutrient response
## 10.1186/s12974-017-0840-7; n=; months:3-5 vs 20-24
#Microglia
raw_data <- read_csv("D:/Rpackages/AgingBioMolecules/data/Mus musculus/Proteome/Article_sup/Flowers_et_al_2017_J Neuroinflammation_Microglia_samples.csv",skip = 2)
colnames(raw_data) <- c("Gene.Symbol","FC","p.value","Zscore")
raw_data$Gene.Symbol <- str_replace_all(raw_data$Gene.Symbol,"\\*","")
fdata <- raw_data %>% mutate("log2FC"=log2(FC),"Significance"=ifelse(p.value<0.05,"sig","nonsig"))
Mice_aging_proteome <- addData2list(list= Mice_aging_proteome, tissue = "Microglia",fdata = fdata,title = "Flowers_2017_J Neuroinflammation",
                                    Title = "Proteomic analysis of aged microglia: shifts in transcription, bioenergetics, and nutrient response",
                                    doi = "10.1186/s12974-017-0840-7",stat.method = "Welch's t-test",age = "4m-22m",age.record = "3–5 months old vs 20–24 months old",
                                    specie.strain = "C57BL/6N")

#5. Title: Biochemical isolation of myonuclei employed to define changes to the myonuclear proteome that occur with aging
## 10.1111/acel.12604; n=10; months:3 vs 24
#brain.nuclei
raw_data <- read_csv("D:/Rpackages/AgingBioMolecules/data/Mus musculus/Proteome/Article_sup/Cutler_et_al_2017_Aging Cell_Brain nuclei_samples.csv",skip = 1)
colnames(raw_data)[4] <- c("log2FC")
raw_data$p.value <- 10^(-raw_data$`-log 10 t.test`)
fdata <- raw_data;colnames(fdata)[1:2] <- c("Uniprot.id","Gene.Symbol")
fdata$Significance <- ifelse(raw_data$p.value<0.05,"sig","nonsig")
Mice_aging_proteome <- addData2list(list= Mice_aging_proteome, tissue = "Brain.Nuclei",fdata = fdata,title = "Cutler_2017_Aging Cell",
                                    Title = "Biochemical isolation of myonuclei employed to define changes to the myonuclear proteome that occur with aging",
                                    doi = "10.1111/acel.12604",stat.method = "t-test",age = "3m-24m",n=10,specie.strain = "C57BL/6",specie.sex = "male")

#6. Title: Proteomic analysis and functional characterization of mouse brain mitochondria during aging reveal alterations in energy metabolism
## 10.1002/pmic.201400277; n=6;
#months:5 vs 12
#Brain.mitochondria
raw_data <- read_csv("D:/Rpackages/AgingBioMolecules/data/Mus musculus/Proteome/Article_sup/Stauch_et_al_2015_Proteomics_Brain mitochondria_6samples_5vs12m.csv")
expr <- as.matrix(raw_data[3:8])
pdata <- data.frame("samples"=colnames(expr),
                    "group"=factor(c(rep("Young",3),rep("Old",3)),levels = c("Young","Old")),
                    "month"=c(rep("5",3),rep("12",3)))
log2FC_temp <- apply(expr, 1, function(x) log2(mean(x[4:6]/x[1:3])))
fdata <- raw_data[-3:-8] %>% mutate("log2FC"=log2FC_temp,"p.value"=`P-value`,"Significance"=ifelse(p.value<0.05,"sig","nonsig"),"Gene.Symbol"=get_genesymbol(`Protein Name`))
colnames(fdata)[1] <- "Uniprot.id"
Mice_aging_proteome <- addData2list(list= Mice_aging_proteome, tissue = "Brain.mitochondria",fdata = fdata,expr = expr,pdata = pdata,title = "Stauch_2015_Proteomics.1",
                                    Title = "Proteomic analysis and functional characterization of mouse brain mitochondria during aging reveal alterations in energy metabolism",
                                    doi = "10.1002/pmic.201400277",age = "5m-12m",n=6,specie.strain = "C57BL/6",specie.sex = "male")
#months:5 vs 24
#Brain.mitochondria
raw_data <- read_csv("D:/Rpackages/AgingBioMolecules/data/Mus musculus/Proteome/Article_sup/Stauch_et_al_2015_Proteomics_Brain mitochondria_6samples_5vs24m.csv")
expr <- as.matrix(raw_data[3:8])
pdata <- data.frame("samples"=colnames(expr),
                    "group"=factor(c(rep("Young",3),rep("Old",3)),levels = c("Young","Old")),
                    "month"=c(rep("5",3),rep("24",3)))
log2FC_temp <- apply(expr, 1, function(x) log2(mean(x[4:6]/x[1:3])))
fdata <- raw_data[-3:-8] %>% mutate("log2FC"=log2FC_temp,"p.value"=`P-value`,"Significance"=ifelse(p.value<0.05,"sig","nonsig"),"Gene.Symbol"=get_genesymbol(`Protein Name`))
Mice_aging_proteome <- addData2list(list= Mice_aging_proteome, tissue = "Brain.mitochondria",fdata = fdata,expr = expr,pdata = pdata,title = "Stauch_2015_Proteomics.2",
                                    Title = "Proteomic analysis and functional characterization of mouse brain mitochondria during aging reveal alterations in energy metabolism",
                                    doi = "10.1002/pmic.201400277",age = "5m-24m",n=6,specie.strain = "C57BL/6",specie.sex = "male")

#months:12 vs 24
#Brain.mitochondria
raw_data <- read_csv("D:/Rpackages/AgingBioMolecules/data/Mus musculus/Proteome/Article_sup/Stauch_et_al_2015_Proteomics_Brain mitochondria_6samples_12vs24m.csv")
expr <- as.matrix(raw_data[3:8])
pdata <- data.frame("samples"=colnames(expr),
                    "group"=factor(c(rep("Young",3),rep("Old",3)),levels = c("Young","Old")),
                    "month"=c(rep("12",3),rep("24",3)))
log2FC_temp <- apply(expr, 1, function(x) log2(mean(x[4:6]/x[1:3])))
fdata <- raw_data[-3:-8] %>% mutate("log2FC"=log2FC_temp,"p.value"=`P-value`,"Significance"=ifelse(p.value<0.05,"sig","nonsig"),"Gene.Symbol"=get_genesymbol(`Protein Name`))
Mice_aging_proteome <- addData2list(list= Mice_aging_proteome, tissue = "Brain.mitochondria",fdata = fdata,expr = expr,pdata = pdata,title = "Stauch_2015_Proteomics.3",
                                    Title = "Proteomic analysis and functional characterization of mouse brain mitochondria during aging reveal alterations in energy metabolism",
                                    doi = "10.1002/pmic.201400277",age = "12m-24m",n=6,specie.strain = "C57BL/6",specie.sex = "male")

#7. Title: Accurate Quantification of More Than 4000 Mouse Tissue Proteins Reveals Minimal Proteome Changes During Aging
## 10.1074/mcp.M110.004523; n=8; months:5 vs 26
#Hippocampus
#Expression data is provided as log2 SILAC ratios between biological samples and the corresponding metabolically labeled protein standard
raw_data <- read_csv("D:/Rpackages/AgingBioMolecules/data/Mus musculus/Proteome/Article_sup/Walther_et_al_2011_MOL CELL PROTEOMICS_Hippocampus_8samples.csv")
expr <- as.matrix(raw_data[2:9])
pdata <- data.frame("samples"=colnames(expr),
                    "group"=factor(c(rep("Young",4),rep("Old",4)),levels = c("Young","Old")),
                    "month"=c(rep("5",4),rep("26",4)))
fdata <- raw_data[-2:-9] %>% mutate("Difference"=`T-test Difference`,"p.value"=10^(-`T-test p value (-log10)`),
                                    "Significance"=ifelse(p.value<0.05,"sig","nonsig"),"p.adj"=p.adjust(p.value,"BH"))
colnames(fdata)[c(4,14,15)] <- c("log2FC","Gene.Symbol","Uniprot.id")
Mice_aging_proteome <- addData2list(list= Mice_aging_proteome, tissue = "Hippocampus",fdata = fdata,expr = expr,pdata = pdata,title = "Walther_2011_MOL CELL PROTEOMICS",
                                    Title = "Accurate Quantification of More Than 4000 Mouse Tissue Proteins Reveals Minimal Proteome Changes During Aging",value.type = "log2 SILAC ratios",
                                    doi = "10.1074/mcp.M110.004523",age = "5m-26m",n=8,specie.sex = "female",specie.strain = "C57BL/6JN")
#Frontal Cortex
#Expression data is provided as log2 SILAC ratios between biological samples and the corresponding metabolically labeled protein standard
raw_data <- read_csv("D:/Rpackages/AgingBioMolecules/data/Mus musculus/Proteome/Article_sup/Walther_et_al_2011_MOL CELL PROTEOMICS_Frontal Cortex_8samples.csv")
expr <- as.matrix(raw_data[2:9])
pdata <- data.frame("samples"=colnames(expr),
                    "group"=factor(c(rep("Young",4),rep("Old",4)),levels = c("Young","Old")),
                    "month"=c(rep("5",4),rep("26",4)))
fdata <- raw_data[-2:-9] %>% mutate("Difference"=`T-test Difference`,"p.value"=10^(-`T-test p value (-log10)`),
                                    "Significance"=ifelse(p.value<0.05,"sig","nonsig"),"p.adj"=p.adjust(p.value,"BH"))
colnames(fdata)[c(4,14,15)] <- c("log2FC","Gene.Symbol","Uniprot.id")
Mice_aging_proteome <- addData2list(list= Mice_aging_proteome, tissue = "Frontal Cortex",fdata = fdata,expr = expr,pdata = pdata,title = "Walther_2011_MOL CELL PROTEOMICS",
                                    Title = "Accurate Quantification of More Than 4000 Mouse Tissue Proteins Reveals Minimal Proteome Changes During Aging",value.type = "log2 SILAC ratios",
                                    doi = "10.1074/mcp.M110.004523",age = "5m-26m",n=8,specie.sex = "female",specie.strain = "C57BL/6JN")
#Cerebellum
raw_data <- read_csv("D:/Rpackages/AgingBioMolecules/data/Mus musculus/Proteome/Article_sup/Walther_et_al_2011_MOL CELL PROTEOMICS_CerebHeartKidney_samples.csv")
fdata <- raw_data[-c(2,3,5,6)]
colnames(fdata)[c(2,12,13)] <- c("log2FC","Gene.Symbol","Uniprot.id")
Mice_aging_proteome <- addData2list(list= Mice_aging_proteome, tissue = "Cerebellum",fdata = fdata, title = "Walther_2011_MOL CELL PROTEOMICS",
                                    Title = "Accurate Quantification of More Than 4000 Mouse Tissue Proteins Reveals Minimal Proteome Changes During Aging",value.type = "log2FC",
                                    doi = "10.1074/mcp.M110.004523",age = "5m-26m",n=2,specie.sex = "female",specie.strain = "C57BL/6JN")

#8. Title: Aging in Mouse Brain Is a Cell/Tissue-Level Phenomenon Exacerbated by Proteasome Loss
## 10.1021/pr100059j; n=; weeks:embryonic day (ED) 10 to 100 weeks ()
#Brain
#Expression data is provided as log2 SILAC ratios between biological samples and the corresponding metabolically labeled protein standard
FC_data <- read_csv("D:/Rpackages/AgingBioMolecules/data/Mus musculus/Proteome/Article_sup/Mao_et_al_2010_J Proteome Res_brain_samples_15TP_FC.csv")
p_data <- read_csv("D:/Rpackages/AgingBioMolecules/data/Mus musculus/Proteome/Article_sup/Mao_et_al_2010_J Proteome Res_brain_samples_15TP_p.csv")
f_data <- read_csv("D:/Rpackages/AgingBioMolecules/data/Mus musculus/Proteome/Article_sup/Mao_et_al_2010_J Proteome Res_brain_samples_15TP_fdata.csv",skip = 1)
f_data <- f_data %>% arrange(desc(`Mowse Score`));f_data <- f_data[!duplicated(f_data$`Spot ID`),]
FC_data$`Spot ID` <- paste0("S",FC_data$`Spot ID`);p_data$`Spot ID` <- paste0("S",p_data$`Spot ID`);f_data$`Spot ID` <- paste0("S",f_data$`Spot ID`)
f_data <- as.data.frame(f_data);rownames(f_data) <- f_data$`Spot ID`;
FC_data <- FC_data[FC_data$`Spot ID` %in% f_data$`Spot ID`,];p_data <- p_data[p_data$`Spot ID` %in% f_data$`Spot ID`,]
f_data <- f_data[p_data$`Spot ID`,]
#months: 5(22weeks) vs 25(100weeks)
fdata <- cbind(f_data,"FC"=1/FC_data$`22w vs 100w`,"p.value"=p_data$`22w vs 100w`);fdata <- fdata[!is.na(fdata$FC),];fdata$Significance <- ifelse(fdata$p.value<0.05,"sig","nonsig")
fdata <- fdata %>% mutate("log2FC"=log2(FC));colnames(fdata)[2] <- "Gene.Symbol"
Mice_aging_proteome <- addData2list(list= Mice_aging_proteome, tissue = "Brain",fdata = fdata,title = "Mao_2010_J Proteome Res.1",
                                      Title = "Aging in Mouse Brain Is a Cell/Tissue-Level Phenomenon Exacerbated by Proteasome Loss",value.type = "stat",
                                      doi = "10.1021/pr100059j",age.record = "22weeks vs 100weeks",age = "5m-25m",specie.sex = "male",specie.strain = "C57BL/6",n=12)
#weeks:embryonic day (ED) 10 to 100 weeks: ET10 ET12 ET14 ET16 ET18 NB 1w 2w 4w 8w 14w 22w 42w 75w 100w
expr_fc <- as.matrix(FC_data[,str_detect(colnames(FC_data),".*100w")]);expr_fc <- cbind(expr_fc,"100w"=rep(1,196))
colnames(expr_fc) <- str_replace_all(colnames(expr_fc)," vs 100w","")
expr_p <- as.matrix(p_data[,str_detect(colnames(p_data),".*100w")]);expr_p <- cbind(expr_p,"100w"=rep(NA,196))
colnames(expr_p) <- str_replace_all(colnames(expr_fc)," vs 100w","")
fdata <- f_data;colnames(fdata)[2] <- "Gene.Symbol"
Mice_aging_proteome <- addData2list(list= Mice_aging_proteome, tissue = "Brain",fdata = fdata,expr_fc = expr_fc,expr_p = expr_p,pdata = pdata,title = "Mao_2010_J Proteome Res.2",
                                    Title = "Aging in Mouse Brain Is a Cell/Tissue-Level Phenomenon Exacerbated by Proteasome Loss",value.type = "stat",
                                    doi = "10.1021/pr100059j",age = "ET10-ET12-ET14-ET16-ET18-NB-1w-2w-4w-8w-14w-22w-42w-75w-100w",group.n = 15,group.type = "timepoint",
                                    specie.sex = "male",specie.strain = "C57BL/6",n=64)


save(Mice_aging_proteome,file = "D:/Rpackages/AgingBioMolecules/data/Mus musculus/Proteome/Mice_aging_proteome.Rdata")

