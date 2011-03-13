############ demos iris
library('FisherEM')

data(iris)
Y  = iris[,-5]
lbl= as.numeric(iris[,5])
K  = 3
Tinit = t(rmultinom(150,1,c(rep(1/K,K))))

fstep.GramSc <- function(Y,T,kernel){
	n = nrow(Y)
	p = ncol(Y)
	K = ncol(T)
	m = colMeans(Y)
	d = min(p-1,(K-1))

	# Compute S
# browser()
	XX = as.matrix(Y - t(m*t(matrix(1,n,p))))
	TT = t(apply(T,1,"/",sqrt(colSums(T))))
	
	if (n>p & kernel==''){
		S = t(XX) %*% XX /n
		B = t(XX)%*%TT%*%t(TT)%*%XX / n
		
		# Eigendecomposition of S^-1 * B
		eig = svd(ginv(S)%*%B,nu=d,nv=0)
		U = eig$u[,1:d]
	}
	else{
		cat('Kernel mode!\n')
		if (n<p | kernel=='linear') G = XX %*% t(XX)
		if (kernel=='rbf') {sigma=1; G = as.matrix(exp(dist(XX,diag=T)^2/(2*sigma^2)))}
		if (kernel=='sigmoid') {a=1;r=0.1;G = tanh(a * X %*% t(X) + r)}
		lambda = 0
		S = G %*% G + lambda*diag(n)
		B = G %*% TT %*% t(TT) %*% G
		H = svd(ginv(S)%*%B,nu=d,nv=0)$u[,1:d]
		V = svd(t(Y) %*% H,nu=d,nv=0)$u[,1:d]
	}
}


V  = fstep.GramSc(Y,Tinit,'')
cls = max.col(Tinit)

for (i in 1:15){
# x11()
# def.par <- par(no.readonly = TRUE)
x = as.matrix(Y)%*%V
min1= round(min(x[,1]),1)-1
max1= round(max(x[,1],1))+1
min2= round(min(x[,2]),1)-1
max2= round(max(x[,2]),1)+1
x1 = split(x[,1],cls)
x2 = split(x[,2],cls)

xhist1 <- hist(x1$'1',breaks=seq(min1,max1,0.1), plot=FALSE) 
xhist2 <- hist(x1$'2',breaks=seq(min1,max1,0.1), plot=FALSE) 
xhist3 <- hist(x1$'3',breaks=seq(min1,max1,0.1), plot=FALSE) 
yhist1 <- hist(x2$'1',breaks=seq(min2,max2,0.1), plot=FALSE) 
yhist2 <- hist(x2$'2',breaks=seq(min2,max2,0.1), plot=FALSE) 
yhist3 <- hist(x2$'3',breaks=seq(min2,max2,0.1), plot=FALSE) 
topX <- max(c(xhist1$counts,xhist2$counts,xhist3$counts))
topY <- max(c(yhist3$counts,yhist3$counts,yhist3$counts))
nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE) 
xrange <- c(min1,max1)
yrange <- c(min2,max2)

par(mar=c(3,3,1,1)) 
plot(x,col=cls,pch=cls,xlim=xrange, ylim=yrange,)
# , xlim=xrange, ylim=yrange, xlab="", ylab="") 
par(mar=c(0,3,1,1)) 
barplot(xhist1$counts, axes=FALSE,col='gray',ylim=c(0,topX), space=0) 
barplot(xhist2$counts, axes=FALSE,col='red', ylim=c(0,topX), space=0,add=T)
barplot(xhist3$counts, axes=FALSE,col='lightgreen', ylim=c(0,topX), space=0,add=T)
par(mar=c(3,0,1,1)) 
barplot(yhist1$counts, axes=FALSE,col='gray', xlim=c(0,topY), space=0, horiz=TRUE) 
barplot(yhist2$counts, axes=FALSE,col='red', xlim=c(0,topY), space=0, horiz=TRUE,add=T) 
barplot(yhist3$counts, axes=FALSE,col='lightgreen', xlim=c(0,topY), space=0, horiz=TRUE,add=T) 
res = fem(Y,3,init='user',Tinit=Tinit,model='AkB',maxit=1)
Tinit = res$P
cls = res$cls
V = res$V
}

