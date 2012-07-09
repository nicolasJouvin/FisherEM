fem <-
function(Y,K,init='random',maxit=100,eps=1e-6,Tinit=c(),model='AkjBk',kernel='',graph=F,Hess=F,method='SVD',l1=0.3,l2=0,nbit=2){

MOD = c('DkBk','DkB','DBk','DB','AkjBk','AkjB','AkBk','AkB','AjBk','AjB','ABk','AB','all')
MET = c('SVD','REG','GS','sparse')
KER = c('','sigmoid','linear','rbf')
if (any(!model%in%MOD)) stop("Invalid model name\n",call.=FALSE)
if (any(!method%in%MET)) stop("Invalid method name\n",call.=FALSE)  
if (any(!kernel%in%KER)) stop("Invalid kernel name\n",call.=FALSE)
if (K>=ncol(Y)) stop("K must be strictly less than the number of variables",call.=FALSE)

      if ( model=='all' & method!='sparse') {
	    resultat = pairlist()
	    bic = rep(NA,12)
		  for(i in 1:12){
		    try(resultat[[i]]<-fem_main(Y,K,init=init,maxit=maxit,eps=eps,Tinit=Tinit,model=MOD[i],kernel=kernel,graph=F,Hess=F,method=method))
		    try(bic[i]<-resultat[[i]]$bic) 
		    cat('\n','model:',MOD[i],'   bic:',resultat[[i]]$bic)
		    }
          
	     max_bic = max(bic,na.rm=T)
       
	    for(i in 1:12){ if(max_bic==bic[i]) res=resultat[[i]] }
      
	    cat('\n')
	    cat('\n','The best model is:',res$prms$model, 'with a BIC equal to:',max_bic,'\n') 
	    cat('\n')
      }

      else if (model=='all' & method=='sparse'){cat('\n','The user needs to choose in first the model before using the sparse option','\n')
						break 
      }

      else if (model!='all' & method=='sparse'){ res0 = fem_main(Y,K,init=init,maxit=50,eps=eps,Tinit=Tinit,model=model,kernel=kernel,graph=F,Hess=F,method='REG') 
						 res  = fem_sparse(Y,K,maxit=15,eps=eps,Tinit=res0$P,model=model,method='sparse',l1=l1,l2=l2,nbit=nbit)
						 cat('\n','The model is:',res$prms$model,'with a l1 term equal to:',l1, 'and with a BIC equal to:',res$bic,'\n')
						 cat('\n')
      }
      
      else { res = fem_main(Y,K,init=init,maxit=maxit,eps=eps,kernel=kernel,graph=graph,Hess=Hess,Tinit=Tinit,model=model,method=method) 
             cat('\n','The model is:',res$prms$model, 'with a BIC equal to:',res$bic,'\n')
	     cat('\n')
      }

class(res)='fem'
res
}

