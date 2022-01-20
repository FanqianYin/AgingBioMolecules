#plot.AgingBioMolecules
#keep.data c("all","sig")

plot.AgingBioMolecules <- function(query, type="Gene_name", specie="Mus musculus", omic="proteome", plot.type="dotplot",keep.data="sig"){
  #query <- "Ctsd"
  #plot single molecule
  require(tidyverse)
  require(RColorBrewer)
  #create col
  col.red <- brewer.pal(9,"Set1")[1];col.blue <- brewer.pal(9,"Set1")[2];col.grey <- brewer.pal(9,"Set1")[9];col.green <- brewer.pal(9,"Set1")[3]
  load("./data/Mice_aging_proteome.SummaryTable.Rdata")
  if(length(query)==1){
    if(specie=="Mus musculus"){
      if(omic=="proteome"){
        if(query %in% DEA.information$Gene.Symbol){
          if(plot.type=="dotplot"){
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
            print(p)
          }else if(plot.type=="boxplot"){
            load("D:/Rpackages/AgingBioMolecules/data/Mus musculus/Proteome/Mice_aging_proteome.Rdata")
            if(query %in% DEA.information$Gene.Symbol.expr){
              gp_df <- data.frame()
              p.list <- list()
              for (t in names(Mice_aging_proteome)) {
                for (d in names(Mice_aging_proteome[[t]])) {
                  id=paste0(t,".",d)
                  if(any(names(Mice_aging_proteome[[t]][[d]])=="expr") &
                     (query %in% Mice_aging_proteome[[t]][[d]][["fdata"]][["Gene.Symbol"]])){
                    query.boolean <- query==Mice_aging_proteome[[t]][[d]][["fdata"]][["Gene.Symbol"]]
                    query.boolean[is.na(query.boolean)] <- FALSE
                    query.exp <- Mice_aging_proteome[[t]][[d]][["expr"]][query.boolean,]
                    if(is.matrix(query.exp)) query.exp <- query.exp[1,]
                    query.exp <- scale(query.exp)
                    samp.n <- nrow(Mice_aging_proteome[[t]][[d]][["pdata"]])
                    query.df <- data.frame("id"=rep(id,samp.n),"tissue"=rep(t,samp.n),
                                           Mice_aging_proteome[[t]][[d]][["pdata"]][c("samples","group","month")],"expr"=query.exp)
                    query.df$group <- as.character(query.df$group)
                    query.df$group[query.df$group=="Yung"] <- "Young"
                    gp_df <- rbind(gp_df,query.df)
                    p.list[[id]] <- ggplot(gp_df,aes(x = group,y=expr,fill=group))+geom_boxplot()+
                      scale_color_manual(values = c("Young"=col.green,"Old"=col.red))+
                      theme_bw()+ggtitle(id)
                  }

                }
              }
             # p <- ggplot(gp_df,aes(x = group,y=expr,fill=group))+geom_boxplot()+
             #   scale_color_manual(values = c("Young"=col.green,"Old"=col.red))+
             #   theme_bw()+
             #   facet_grid(vars(id))
              gridExtra::grid.arrange(grobs=p.list)


            }else print(paste0("Protein ",query," does not have expression data in database."))
          }



        }else print(paste0("Protein ",query," does not included in database."))


      }


    }



  }
  #plot multi molecule
  if(length(query)>1){




  }



  #print(gp_df)
  #print(p)
  #if(exists("p.list")) gridExtra::grid.arrange(grobs=p.list)
}

#next to do
#1.add sec x axis with age.group informaiton
#2.reorder dotplot
