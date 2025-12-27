files <- list.files('data/cache', pattern='\\.rds$', full.names=TRUE)
for(f in files){
  cat('---', f, '---\n')
  obj <- tryCatch(readRDS(f), error=function(e) {cat('READ ERROR:', e$message,'\n'); NULL})
  if(!is.null(obj)){
    cat('Class:',class(obj),'\n')
    if(inherits(obj,'xts')||inherits(obj,'zoo')){
      print(head(obj))
    } else if(is.data.frame(obj)||is.matrix(obj)){
      print(head(obj))
    } else {
      print(str(obj, max.level=1))
    }
  }
}
