
# An implementation of the randomization test for periodicity of a time serie, as described in Mainly, B. "Randomization, Bootstrapp and Monte Carlo Method in Biology" 
# 2nd edition - CAP 11
# fev 2024  v4  msantos@est.ufmg.br

# For more details, see a complete description of this test in Mainly, B. "Randomization, Bootstrapp and Monte Carlo Method in Biology" 


# vpk function 
# Return pk values given vector serie x (n even)
vpk<- function(x){
	n=length(x)
	m=n/2
	# A and B coefficients
	i=1:n
	A=rep(0,m-1)
	A0=mean(x)
	Am=(1/n)*sum((x*(-1)^i))
	for(k in 1:(m-1)){
		wk=2*pi*k/n
		A[k]=(2/n)*sum(x*cos(wk*i))
}
	B=rep(0,m-1)
	for(k in 1:(m-1)){
		wk=2*pi*k/n
		B[k]=(2/n)*sum(x*sin(wk*i))
	}
	SK=A^2 + B^2
	y=n*SK
	k=1:length(y)
	w=2*pi*k/n
	pk0=A0^2/((n-1)*var(x))
	pk=SK/((n-1)*var(x))
	pkm=Am^2/((n-1)*var(x))
	pk=c(pk,pkm)
	pk=round(pk/(sum(pk)),4)   # sum(p(k))=1
	return(pk)
}

TestCycle=function(serie, plim, N){
	# input:  vector serie, plim (default=0.05),  N(number of permutations)
	# return data frame with    k, w(k), p(k), sig.level 
	# serie = vector 
	# plim = only estimated p-value < plim will be printed
	# N = number of permutations
	
	pk=vpk(serie)
	n=length(pk)
	# permutacoes
	PSIM=NULL
	for (i in 1:N){
		serieP=sample(serie)
		p=vpk(serieP)
		PSIM=rbind(PSIM,p)
	}
	sig.level=rep(0,ncol(PSIM))
	for (k in 1:ncol(PSIM)){
		sig.level[k]=sum(PSIM[,k]>pk[k])/nrow(PSIM)
	}
	n=length(serie)
	j=1:(n/2)  # n even 
	w=2*pi*j/n
	w=round(w,3)
	cycle=round(n/j,1)
	k=j
	# data.frame with results
	V=data.frame(k,w,cycle,pk,sig.level)
	return(subset(V, sig.level<plim))

}


# An example with the dataset extracted from the book
vserie=read.table("serie.grain.txt", header=T)
head(vserie)

# Usage
TestCycle(serie=vserie$Yield, plim=0.05, N=3000) # plim = show results only if p-value <= plim , N is the number of permutations





