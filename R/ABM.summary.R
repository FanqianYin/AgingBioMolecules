

ABM.summary <- function(query,plot=T,keep.table=F){
  if(query=="summary.articles"){


  }


}

load("D:/Rpackages/AgingBioMolecules/data/Summary/Mice_aging_proteome.SummaryTable.Rdata")
library(RColorBrewer)
library(tidyverse)
col.red <- brewer.pal(9,"Set1")[1];col.blue <- brewer.pal(9,"Set1")[2];col.grey <- brewer.pal(9,"Set1")[9];col.green <- brewer.pal(9,"Set1")[3]
#plotting of summary table information

#article.pub years
df.plot <- data.summary.df %>% filter(!duplicated(id2));df.plot <-  table(df.plot$article.pub.year)  %>% as.data.frame()
colnames(df.plot) <- c("Year","numbers")
ggplot(df.plot,aes(Year,numbers,fill=numbers))+geom_bar(stat = "identity")+
  theme_classic()+scale_y_continuous(breaks = seq(1,10) )+ggtitle("Publication numbers")+scale_fill_gradient(low = "lightblue", high = "darkblue", na.value = NA)+
  theme(axis.text.x = element_text(size=15))
#article.journal
df.plot <- data.summary.df %>% filter(!duplicated(id2));df.plot <-  table(df.plot$article.journal) %>% sort(decreasing = T)  %>% as.data.frame()
colnames(df.plot) <- c("Journals","numbers");df.plot$Journals <- factor(df.plot$Journals,levels = df.plot$Journals)
ggplot(df.plot,aes(Journals,numbers,fill=numbers))+geom_bar(stat = "identity")+
  theme_classic()+scale_y_continuous(breaks = seq(1,10) )+ggtitle("Journal numbers")+scale_fill_gradient(low = "lightblue", high = col.red, na.value = NA)+
  theme(axis.text.x = element_text(angle=30,hjust = 0.9,size=12))
#sample sex
df.plot <- data.summary.df %>% filter(!duplicated(id2));df.plot <-  table(df.plot$sample.sex) %>% sort(decreasing = T)  %>% as.data.frame()
colnames(df.plot) <- c("Sex","numbers")
ggplot(df.plot,aes(Sex,numbers,fill=numbers))+geom_bar(stat = "identity")+coord_polar()+ggtitle("Journal numbers")+theme_bw()
ggplot(df.plot, aes(x="", y=numbers, fill=Sex)) +geom_col(color = "black")+
  geom_bar(stat="identity", width=1) +
  coord_polar(theta = "y")+geom_text(aes(label = numbers),
                                      position = position_stack(vjust = 0.5))+theme_classic()
#tissue
df.plot <- data.summary.df %>% filter(!duplicated(id2));df.plot <-  table(df.plot$tissue) %>% sort(decreasing = T)  %>% as.data.frame()
colnames(df.plot) <- c("Tissues","numbers");df.plot$Tissues <- factor(df.plot$Tissues,levels = df.plot$Tissues)
ggplot(df.plot,aes(Tissues,numbers,fill=numbers))+geom_bar(stat = "identity")+
  theme_classic()+scale_y_continuous(breaks = seq(1,10) )+ggtitle("Tissue numbers")+scale_fill_gradient(low = col.green, high = col.red, na.value = NA)+
  theme(axis.text.x = element_text(angle=30,hjust = 0.9,size=12))
