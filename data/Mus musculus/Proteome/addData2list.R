#n sample numbers
addData2list <- function(list=NULL,expr=NULL,expr_fc=NULL,expr_p=NULL,expr_mean=NULL,n=NA,tissue=NULL,fdata=NULL,pdata=NULL,sdata=NULL,title=NA,Title=NA,doi=NA,
                         value.type=NA,tissue.name=NA,stat.method=NA,age=NA,age.record=NA,
                         MS.equipment=NA,MS.platform=NA,specie.strain=NA,specie.sex=NA,group.n=2,group.type="Young-Old",quant.method=NA){
  if(!is.null(doi)) list[[tissue]][[title]][["metadata"]][["doi"]] <- doi
  if(!is.null(Title)) list[[tissue]][[title]][["metadata"]][["Title"]] <- Title
  if(!is.null(value.type)) list[[tissue]][[title]][["metadata"]][["value.type"]] <- value.type
  if(!is.null(tissue.name)) list[[tissue]][[title]][["metadata"]][["tissue.name"]] <- tissue.name
  if(!is.null(stat.method)) list[[tissue]][[title]][["metadata"]][["stat.method"]] <- stat.method
  if(!is.null(MS.equipment)) list[[tissue]][[title]][["metadata"]][["MS.equipment"]] <- MS.equipment
  if(!is.null(MS.platform)) list[[tissue]][[title]][["metadata"]][["MS.platform"]] <- MS.platform
  if(!is.null(specie.strain)) list[[tissue]][[title]][["metadata"]][["specie.strain"]] <- specie.strain
  if(!is.null(specie.sex)) list[[tissue]][[title]][["metadata"]][["specie.sex"]] <- specie.sex
  if(!is.null(group.type)) list[[tissue]][[title]][["metadata"]][["group.type"]] <- group.type
  if(!is.null(age)) list[[tissue]][[title]][["metadata"]][["age"]] <- age
  if(!is.null(age.record)) list[[tissue]][[title]][["metadata"]][["age.record"]] <- age.record
  if(!is.null(n)) list[[tissue]][[title]][["metadata"]][["n"]] <- n
  if(!is.null(quant.method)) list[[tissue]][[title]][["metadata"]][["quant.method"]] <- quant.method

  if(!is.null(expr)) list[[tissue]][[title]][["expr"]] <- expr
  if(!is.null(expr)) list[[tissue]][[title]][["expr_fc"]] <- expr_fc
  if(!is.null(expr)) list[[tissue]][[title]][["expr_p"]] <- expr_p
  if(!is.null(expr)) list[[tissue]][[title]][["expr_mean"]] <- expr_mean
  #list[[tissue]][[title]][["expr_zscore"]] <- scale(expr %>% t) %>% t
  #list[[tissue]][[title]][["expr_zscore"]] <- scale(log2(expr_data1+1) %>% t) %>% t
  if(!is.null(fdata)) list[[tissue]][[title]][["fdata"]] <- fdata
  if(!is.null(sdata)) list[[tissue]][[title]][["sdata"]] <- sdata
  if(!is.null(pdata)) list[[tissue]][[title]][["pdata"]] <- pdata
  if(!is.null(group.n)) list[[tissue]][[title]][["group.n"]] <- group.n
  return(list)
}

#get genesymbol from description value like: Ferritin heavy chain OS=Mus musculus GN=Fth1 PE=1 SV=2
get_genesymbol <- function(vect){
  vect.r <- str_replace(vect,pattern = ".+GN=",replacement = "")
  vect.r <- str_replace(vect.r,pattern = " PE=.+",replacement = "")
  vect.r[str_detect(vect.r,".* OS=.*")] <- NA
  return(vect.r)
}

