print.fem <- function(x,...){
 cat('* Model: ')
 cat('The chosen model is ',as.character(x$model),' with K = ',x$K,' (',x$crit,'=',x$critValue,')\n',sep='')
 if (nrow(x$U)<10){ cat('* Loading matrix:\n'); print(x$U)}
 else if(x$call[[1]]=='sfem') cat('* Loading matrix:',sum(rowSums(abs(x$U))>1e-5),
                                  'variables are active over the',nrow(x$U),'original ones\n')
}