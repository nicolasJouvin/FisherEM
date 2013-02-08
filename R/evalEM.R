evalEM <-
function(cls,obj,disp=1){

Klbl = length(unique(as.numeric(cls)))
if (Klbl!=obj$prms$K) stop("K is different between the partitions \n",call.=FALSE)

K = Klbl
res = obj$cls

	# Initialization
	t = c()
	lbl = c()

	# Main loop
	P = permutations(K)
	for (i in 1:nrow(P)){
		for (j in 1:K){
			lbl[res==P[i,j]] = j
		}
		t[i] = sum(lbl==cls)/length(cls)
	}
	tx = max(t)
	ind = which.max(t)
	 
	# Compute correct labels
	perm = P[ind,]
	for (j in 1:K){
		lbl[res==perm[j]] = j
	}
	
	# Display
	if (disp){
		cat('- Clustering accuracy rate: ',tx,'\n')
		cat('- Confusion matrix:\n')
		print(table(cls,lbl))
	}
	
	# Return results
	prms=list(lbl=lbl,tx=tx)
}

