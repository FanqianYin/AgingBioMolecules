#plot.AgingBioMolecules
#keep.data c("all","sig")

plot.AgingBioMolecules <- function(query, type="Gene_name", specie="Mus musculus", omic="proteome", plot.type="dotplot",keep.data="sig"){
  #query <- "Psap"
  #plot single molecule
  require(tidyverse)
  require(RColorBrewer)
  #create col
  col.red <- brewer.pal(9,"Set1")[1];col.blue <- brewer.pal(9,"Set1")[2];col.grey <- brewer.pal(9,"Set1")[9]
  load("D:/Rpackages/AgingBioMolecules/data/Summary/Mice_aging_proteome.SummaryTable.Rdata")
  if(length(query)==1){
    if(specie=="Mus musculus"){
      if(omic=="proteome"){
        if(query %in% DEA.information$Gene.Symbol){
          gp_df <- data.frame(data.summary.df[c("id","id2","tissue","age.timepoint")],
                              "Significance"= DEA.information$Significance[query,],
                              "p.value"= DEA.information$p.value[query,],
                              "log2FC"= DEA.information$log2FC[query,]) %>%
            mutate("direction"=ifelse(log2FC>0,"up","down"),
                   "direction"=ifelse(Significance=="sig",direction,"nonsig"))
          if(keep.data=="sig") gp_df <- gp_df[!is.na(gp_df$Significance),]
          p <- ggplot(gp_df,aes(id2,tissue,size=abs(log2FC),col=direction))+
            geom_point()+
            scale_color_manual(values = c("down"=col.blue,"nonsig"=col.grey,"up"=col.red))+
            theme_bw()+
            labs(x="Publications",y="Tissues",title = query)+
            theme(axis.text.x = element_text(angle = 30,hjust = 1), plot.title = element_text(hjust = 0.5))


        }else print(paste0("Protein ",query," does not included in database."))


      }


    }



  }
  #plot multi molecule
  if(length(query)>1){




  }



  #print(gp_df)
  print(p)
}

#next to do
#1.add sec x axis with age.group informaiton
