#get genesymbol from description value like: Ferritin heavy chain OS=Mus musculus GN=Fth1 PE=1 SV=2
get_genesymbol <- function(vect){
  vect.r <- str_replace(vect,pattern = ".+GN=",replacement = "")
  vect.r <- str_replace(vect.r,pattern = " PE=.+",replacement = "")
  vect.r[str_detect(vect.r,".* OS=.*")] <- NA
  return(vect.r)
}
