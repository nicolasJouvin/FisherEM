print.fem <- function(x,...){
  cat('* Model: ')
  cat('The chosen model is ',as.character(x$model),' with K = ',x$K,' (',x$crit,'=',x$critValue,')\n',sep='')
  if (nrow(x$U)<20){ 
    cat('* Loading matrix (zero values are denoted by .):\n')
    U = x$U
    U[U==0] = NA
    U[U<0.001] = 0.001
    print.default(U, digits=3, na.print='.')
  } else if(x$call[[1]]=='sfem') {
    cat('* Loading matrix:',sum(rowSums(abs(x$U))>1e-5),'variables are active over the',nrow(x$U),'original ones\n')
  }
}