## eval ##
#
# Author: 	C.Bouveyron & C.Brunet
# Institute:	SAMM - Paris 1 / Evry
# Date: 	February 2011
# License:	GPL v2
# 
# ***************************************************************
# function which computes clustering accuracy when labels are known.
# Use: evalEM(cls,res,K)
# ***************************************************************
# Call of libraries needed
library('e1071')
# ***************************************************************
#
eval.fem <- function(obj,cls,disp=1){

K = length(unique(cls))
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