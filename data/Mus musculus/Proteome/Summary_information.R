#summary information of database
library(tidyverse)
load("D:/Rpackages/AgingBioMolecules/data/Mus musculus/Proteome/Mice_aging_proteome.Rdata")
#article information
#check data
for (t in names(Mice_aging_proteome)) {
  for (d in names(Mice_aging_proteome[[t]])) {
    id=paste0(t,".",d)
    if(is.null(Mice_aging_proteome[[t]][[d]][["metadata"]][["Title"]])) print(paste0(id,": Lack of Title"))
    if(is.null(Mice_aging_proteome[[t]][[d]][["metadata"]][["doi"]])) print(paste0(id,": Lack of doi"))
    if(is.null(Mice_aging_proteome[[t]][[d]][["metadata"]][["age"]])) print(paste0(id,": Lack of age"))
    if(is.null(Mice_aging_proteome[[t]][[d]][["metadata"]][["n"]])) print(paste0(id,": Lack of sample number n"))
    if(is.null(Mice_aging_proteome[[t]][[d]][["metadata"]][["specie.strain"]])) print(paste0(id,": Lack of specie.strain"))
    if(is.null(Mice_aging_proteome[[t]][[d]][["metadata"]][["specie.sex"]])) print(paste0(id,": Lack of specie.sex"))
    #if(is.null(Mice_aging_proteome[[t]][[d]][["metadata"]][["group.n"]])) print(paste0(id,": Lack of group.n"))
    if(is.null(Mice_aging_proteome[[t]][[d]][["metadata"]][["group.type"]])) print(paste0(id,": Lack of group.type"))
    if(is.null(Mice_aging_proteome[[t]][[d]][["fdata"]])) print(paste0(id,": Lack of fdata"))
    if(is.null(Mice_aging_proteome[[t]][[d]][["fdata"]][["Significance"]])) print(paste0(id,": Lack of Significance"))
    if(is.null(Mice_aging_proteome[[t]][[d]][["fdata"]][["log2FC"]])) print(paste0(id,": Lack of log2FC"))
    if(is.null(Mice_aging_proteome[[t]][[d]][["fdata"]][["p.value"]])) print(paste0(id,": Lack of p.value"))
    # if(is.null(Mice_aging_proteome[[t]][[d]][["fdata"]])) print(paste0(id,": Lack of fdata"))
    # if(is.null(Mice_aging_proteome[[t]][[d]][["fdata"]])) print(paste0(id,": Lack of fdata"))

  }
}
#data.summary.table <- data.frame()
timer <- 1
for (t in names(Mice_aging_proteome)) {
  for (d in names(Mice_aging_proteome[[t]])) {
    d.vect <- list("id"=paste0(t,".",d),
                   "id2"=d,
                   "article.author"=unlist(str_split(string = d,pattern = "_"))[1],
                   "article.pub.year"=unlist(str_split(d,"_"))[2],
                   "article.journal"=str_replace_all(unlist(str_split(d,"_"))[3],"\\..+",""),
                   "article.title"=Mice_aging_proteome[[t]][[d]][["metadata"]][["Title"]],
                   "article.doi"=Mice_aging_proteome[[t]][[d]][["metadata"]][["doi"]],

                   "tissue.name"=Mice_aging_proteome[[t]][[d]][["metadata"]][["tissue.name"]],
                   "tissue"=t,
                   "tissue.super"=NA,

                   "age.record"=Mice_aging_proteome[[t]][[d]][["metadata"]][["age.record"]],
                   "age.young"=unlist(str_split(Mice_aging_proteome[[t]][[d]][["metadata"]][["age"]],"-"))[1],
                   "age.old"=unlist(str_split(Mice_aging_proteome[[t]][[d]][["metadata"]][["age"]],"-"))[2],
                   "age.timepoint"=Mice_aging_proteome[[t]][[d]][["metadata"]][["age"]],

                   "sample.numbers"=Mice_aging_proteome[[t]][[d]][["metadata"]][["n"]],
                   "sample.strain"=Mice_aging_proteome[[t]][[d]][["metadata"]][["specie.strain"]],
                   "sample.sex"=Mice_aging_proteome[[t]][[d]][["metadata"]][["specie.sex"]], # male, female, male-female

                   "group.numbers"=Mice_aging_proteome[[t]][[d]][["metadata"]][["group.n"]],
                   "group.type"=Mice_aging_proteome[[t]][[d]][["metadata"]][["group.type"]], # Old-young, timepoints

                   "data.type"=ifelse("expr" %in% names(Mice_aging_proteome[[t]][[d]]),"expr.mat","stat"), #expr.mat, stat

                   "Protein.n"=nrow(Mice_aging_proteome[[t]][[d]][["fdata"]]),
                   "DEP.n"=ifelse(!is.null(Mice_aging_proteome[[t]][[d]][["fdata"]][["Significance"]]),
                                  sum(Mice_aging_proteome[[t]][[d]][["fdata"]][["Significance"]]=="sig",na.rm = T),
                                  NA),
                   "DEP.up.n"=ifelse(!is.null(Mice_aging_proteome[[t]][[d]][["fdata"]][["Significance"]]) & !is.null(Mice_aging_proteome[[t]][[d]][["fdata"]][["log2FC"]]),
                                     sum(Mice_aging_proteome[[t]][[d]][["fdata"]][["Significance"]]=="sig" & Mice_aging_proteome[[t]][[d]][["fdata"]][["log2FC"]] > 0,na.rm = T),
                                     NA),
                   "DEP.down.n"=ifelse(!is.null(Mice_aging_proteome[[t]][[d]][["fdata"]][["Significance"]]) & !is.null(Mice_aging_proteome[[t]][[d]][["fdata"]][["log2FC"]]),
                                       sum(Mice_aging_proteome[[t]][[d]][["fdata"]][["Significance"]]=="sig" & Mice_aging_proteome[[t]][[d]][["fdata"]][["log2FC"]] < 0,na.rm = T),
                                       NA),
                   "stat.method"=Mice_aging_proteome[[t]][[d]][["metadata"]][["stat.method"]],

                   "protein.separate"=NA, #SDS-PAGE, LC
                   "quant.platform"=NA, #LC/MS
                   "quant.method"=NA, #TMT
                   "equipment"=NA,
                   "Data.acquisition.mode"=NA #DDA, DIA

                   )
    if(timer==1) data.summary.table <- d.vect
    if(timer>1) data.summary.table <- rbind(data.summary.table,d.vect)
    timer <- timer+1
  }
}
data.summary.table <- as.data.frame(data.summary.table);rownames(data.summary.table) <- data.summary.table$id
list2df <- function(dfwithList,first.colname="id"){
  new.df <- data.frame("id"=unlist(dfwithList[["id"]]))
  for (n in names(dfwithList)) {
    new.df[[n]] <- unlist(dfwithList[[n]])
    #print(paste0("Finish: ", n))
  }
  return(new.df)
}
data.summary.df <- list2df(data.summary.table)

#protein information
#check fdata
for (t in names(Mice_aging_proteome)) {
  for (d in names(Mice_aging_proteome[[t]])) {
    id=paste0(t,".",d)
    if(is.null(Mice_aging_proteome[[t]][[d]][["fdata"]][["Gene.Symbol"]])) print(paste0(id,": Lack of fdata-Gene.Symbol"))
    if(is.null(Mice_aging_proteome[[t]][[d]][["fdata"]][["Significance"]])) print(paste0(id,": Lack of fdata-Significance"))
    if(is.null(Mice_aging_proteome[[t]][[d]][["fdata"]][["p.value"]])) print(paste0(id,": Lack of fdata-p.value"))
    if(is.null(Mice_aging_proteome[[t]][[d]][["fdata"]][["log2FC"]])) print(paste0(id,": Lack of fdata-log2FC"))
  }
}
#extract single gene.symbol from gene.symbols & create p.value,sig,log2FC when there is not.
for (t in names(Mice_aging_proteome)) {
  for (d in names(Mice_aging_proteome[[t]])) {
    if(str_detect(Mice_aging_proteome[[t]][[d]][["fdata"]][["Gene.Symbol"]],";")){
      Mice_aging_proteome[[t]][[d]][["fdata"]][["Gene.Symbols"]] <- Mice_aging_proteome[[t]][[d]][["fdata"]][["Gene.Symbol"]]
      Mice_aging_proteome[[t]][[d]][["fdata"]][["Gene.Symbol"]] <- sapply(str_split(Mice_aging_proteome[[t]][[d]][["fdata"]][["Gene.Symbol"]],";"), function(x) x[[1]]) #keep only the firt gene.symbol when occur multiple names.
    }
    if(is.null(Mice_aging_proteome[[t]][[d]][["fdata"]][["Significance"]])) Mice_aging_proteome[[t]][[d]][["fdata"]][["Significance"]] <- rep(NA,nrow(Mice_aging_proteome[[t]][[d]][["fdata"]]))
    if(is.null(Mice_aging_proteome[[t]][[d]][["fdata"]][["p.value"]])) Mice_aging_proteome[[t]][[d]][["fdata"]][["p.value"]] <- rep(NA,nrow(Mice_aging_proteome[[t]][[d]][["fdata"]]))
    if(is.null(Mice_aging_proteome[[t]][[d]][["fdata"]][["log2FC"]])) Mice_aging_proteome[[t]][[d]][["fdata"]][["log2FC"]] <- rep(NA,nrow(Mice_aging_proteome[[t]][[d]][["fdata"]]))
  }
}
save(Mice_aging_proteome,file = "D:/Rpackages/AgingBioMolecules/data/Mus musculus/Proteome/Mice_aging_proteome.Rdata")
#create differentially expressed protein information
Gene.Symbols <- c()
for (t in names(Mice_aging_proteome)) {
  for (d in names(Mice_aging_proteome[[t]])) {
    id=paste0(t,".",d)
    Gene.Symbol <- Mice_aging_proteome[[t]][[d]][["fdata"]][["Gene.Symbol"]]
    Gene.Symbols <- c(Gene.Symbols,Gene.Symbol)

  }
}
Gene.Symbols <- unique(Gene.Symbols);Gene.Symbols <- Gene.Symbols[!is.na(Gene.Symbols)]
#gene with expr data
Gene.Symbols.expr <- c()
for (t in names(Mice_aging_proteome)) {
  for (d in names(Mice_aging_proteome[[t]])) {
    id=paste0(t,".",d)
    if(id %in% data.summary.df$id[data.summary.df$data.type=="expr.mat"]){
      Gene.Symbol <- Mice_aging_proteome[[t]][[d]][["fdata"]][["Gene.Symbol"]]
      Gene.Symbols.expr <- c(Gene.Symbols.expr,Gene.Symbol)
    }
  }
}
Gene.Symbols.expr <- unique(Gene.Symbols.expr);Gene.Symbols.expr <- Gene.Symbols.expr[!is.na(Gene.Symbols.expr)]

Gene.Symbol.df <- data.frame("Gene.Symbol"=Gene.Symbols)
DEA.information <- list()
DEA.information[["Gene.Symbol"]] <- Gene.Symbols
DEA.information[["Gene.Symbol.expr"]] <- Gene.Symbols.expr

DEA.information[["Significance"]] <- Gene.Symbol.df;DEA.information[["p.value"]] <- Gene.Symbol.df;DEA.information[["log2FC"]] <- Gene.Symbol.df
rownames(DEA.information[["Significance"]]) <- Gene.Symbols;rownames(DEA.information[["p.value"]]) <- Gene.Symbols;rownames(DEA.information[["log2FC"]]) <- Gene.Symbols
for (t in names(Mice_aging_proteome)) {
  for (d in names(Mice_aging_proteome[[t]])) {
    id=paste0(t,".",d)
    print(paste0("Processing: ",id))
    fdata <- Mice_aging_proteome[[t]][[d]][["fdata"]]
    fdata <- fdata[order(fdata[["p.value"]]),];fdata <- fdata[!duplicated(fdata[["Gene.Symbol"]]),]
    sig.df <- fdata[c("Gene.Symbol","Significance")];colnames(sig.df)[2] <- id
    p.df <- fdata[c("Gene.Symbol","p.value")];colnames(p.df)[2] <- id
    log2fc.df <- fdata[c("Gene.Symbol","log2FC")];colnames(log2fc.df)[2] <- id
    DEA.information[["Significance"]] <- left_join(DEA.information[["Significance"]],sig.df)
    DEA.information[["p.value"]] <- left_join(DEA.information[["p.value"]],p.df)
    DEA.information[["log2FC"]] <- left_join(DEA.information[["log2FC"]],log2fc.df)

  }
}
DEA.information[["Significance"]] <- as.matrix(DEA.information[["Significance"]][-1])
DEA.information[["p.value"]] <- as.matrix(DEA.information[["p.value"]][-1])
DEA.information[["log2FC"]] <- as.matrix(DEA.information[["log2FC"]][-1])

DEA.information[["AgingProteins_mouse"]] <- Gene.Symbol.df
DEA.information[["AgingProteins_mouse"]]$Detected.n <- apply(DEA.information[["p.value"]][-1],1,function(x) sum(!is.na(x)))
DEA.information[["AgingProteins_mouse"]]$DiffExpr.n <- apply(DEA.information[["Significance"]][-1],1,function(x) sum(x=="sig",na.rm = T))
DEA.information[["AgingProteins_mouse"]]$Expr_up.n <- apply(DEA.information[["log2FC"]][-1],1,function(x) sum(x>0,na.rm = T))
DEA.information[["AgingProteins_mouse"]]$Expr_down.n <- apply(DEA.information[["log2FC"]][-1],1,function(x) sum(x<0,na.rm = T))
# temp_d <- as.matrix(DEA.information[["log2FC"]][-1]);temp_d <- temp_d[as.matrix(DEA.information[["Significance"]][-1])=="sig"]
# DEA.information[["AgingProteins_mouse"]]$DiffExpr_up.n <- apply(DEA.information[["Significance"]][-1],1,function(x) sum(x=="sig",na.rm = T))


dir.create("D:/Rpackages/AgingBioMolecules/data/Summary")
save(data.summary.df,DEA.information,file = "D:/Rpackages/AgingBioMolecules/data/Summary/Mice_aging_proteome.SummaryTable.Rdata")
save(data.summary.df,DEA.information,file = "D:/Rpackages/AgingBioMolecules/data/Mice_aging_proteome.SummaryTable.Rdata")
