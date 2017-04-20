# TODO: Add comment
# 
# Author: ursus
###############################################################################

Kozak_88<-function(DBH,HT,TOTHT,a0,a1,a2,b1,b2,b3,b5){
	
	Z <- HT/TOTHT
	P <- 1.3/TOTHT
	Q <- 1 - (Z)^(1/3)
	X <- Q / (1 - (P)^(1/3))
	
	result<-(a0*DBH^a1*a2^DBH)*X^(b1*Z^2+b2*log(Z+0.0001)+b3*sqrt(Z)+b4*exp(Z)+b5*(DBH/TOTHT))
	
}
 
Kozak_01<-function(DBH,HT,TOTHT,a0,a1,a2,b1,b2,b3){
	
	Z <- HT/TOTHT
	DH<-DBH/TOTHT
	X1<-(1-(Z)^1/4)/(1-0.01^1/4) 
	
	result<-((a0*DBH^a1)*X1^(b0+b1*(1/exp(DH))+b2*DBH^X1+b3*(X1^DH)))
	
}

Kozak_02a<-function(DBH,HT,TOTHT,a0=0,a1=0,a2=0,b1=0,b2=0,b3=0,b4=0,b5=0,b6=0,gradient=TRUE){
	
	if(length(a0)==1){a0<-rep(a0,length(DBH))}
	if(length(a1)==1){a1<-rep(a1,length(DBH))}
	if(length(a2)==1){a2<-rep(a2,length(DBH))}
	if(length(b1)==1){b1<-rep(b1,length(DBH))}
	if(length(b2)==1){b2<-rep(b2,length(DBH))}
	if(length(b3)==1){b3<-rep(b3,length(DBH))}
	if(length(b4)==1){b4<-rep(b4,length(DBH))}
	if(length(b5)==1){b5<-rep(b5,length(DBH))}
	if(length(b6)==1){b6<-rep(b6,length(DBH))}
	
	Z <- HT/TOTHT
	P <- 1.3/TOTHT
	Q <- 1-(Z^1/3)
	X <- Q / (1 - (P)^(1/3))
	
	A0<-a0
	A1<-DBH^a1
	A2<-TOTHT^a2
	A3<-(X^(Z^4))^b1
	A4<-(X^(1/exp(DBH/TOTHT)))^(b2)
	A5<-(X^(X^0.1))^b3
	A6<-(X^(1/DBH))^b4
	A7<-(X^(TOTHT^Q))^b5
	A8<-(X^(X))^b6

	A<-matrix(c(A0,A1,A2,A3,A4,A5,A6,A7,A8),nrow=length(A1))
	
	res<-apply(A,1,prod)
	
	if(gradient) {
		
		deriv_A0<-rep(1,length(DBH))
		deriv_A1<-ifelse(a1==0,1,a1*DBH^(a1-1))
		deriv_A2<-ifelse(a2==0,1,a2*TOTHT^(a2-1))
		deriv_A3<-ifelse(b1==0,1,b1*(X^(Z^4))^(b1-1))
		deriv_A4<-ifelse(b2==0,1,(b2)*(X^(1/exp(DBH/TOTHT)))^(b2-1))
		deriv_A5<-ifelse(b3==0,1,(b3)*(X^(X^0.1))^(b3-1))
		deriv_A6<-ifelse(b4==0,1,(b4)*(X^(X^Q))^(b4-1))
		deriv_A7<-ifelse(b5==0,1,(b5)*(X^(TOTHT^Q))^(b5-1))
		deriv_A8<-ifelse(b6==0,1,(b6)*(X^(X))^(b6-1))
		
		deriv_A<-matrix(c(deriv_A0,deriv_A1,deriv_A2,deriv_A3,deriv_A4,deriv_A5,
						deriv_A6,deriv_A7,deriv_A8),nrow=length(A1))
		
		grad<-matrix(rep(res,times=dim(A)[2]),nrow=length(A1))
		grad<-grad/A
		grad<-grad*deriv_A
		colnames(grad)<-c("a0","a1","a2","b1","b2","b3","b4","b5","b6")
		attr(res, "gradient") <- grad
		
	}
	
	return(res)
	
}

Kozak_02b<-function(DBH,HT,TOTHT,CR,a0=0,a1=0,a2=0,b1=0,b2=0,b3=0,b4=0,b5=0,b6=0,b7=0,gradient=TRUE){
	
	if(length(a0)==1){a0<-rep(a0,length(DBH))}
	if(length(a1)==1){a1<-rep(a1,length(DBH))}
	if(length(a2)==1){a2<-rep(a2,length(DBH))}
	if(length(b1)==1){b1<-rep(b1,length(DBH))}
	if(length(b2)==1){b2<-rep(b2,length(DBH))}
	if(length(b3)==1){b3<-rep(b3,length(DBH))}
	if(length(b4)==1){b4<-rep(b4,length(DBH))}
	if(length(b5)==1){b5<-rep(b5,length(DBH))}
	if(length(b6)==1){b6<-rep(b6,length(DBH))}
	if(length(b7)==1){b7<-rep(b7,length(DBH))}
	
	Z <- HT/TOTHT
	P <- 1.3/TOTHT
	Q <- 1-(Z^1/3)
	X <- Q / (1 - (P)^(1/3))
	
	A0<-a0
	A1<-DBH^a1
	A2<-TOTHT^a2
	A3<-(X^(Z^4))^b1
	A4<-(X^(1/exp(DBH/TOTHT)))^(b2)
	A5<-(X^(X^0.1))^b3
	A6<-(X^(1/DBH))^b4
	A7<-(X^(TOTHT^Q))^b5
	A8<-(X^(X))^b6
	A9<-(X^CR)^b7
	A<-matrix(c(A0,A1,A2,A3,A4,A5,A6,A7,A8,A9),nrow=length(A1))
	
	res<-apply(A,1,prod)
	if(gradient) {
		
		deriv_A0<-rep(1,length(DBH))
		deriv_A1<-ifelse(a1==0,1,a1*DBH^(a1-1))
		deriv_A2<-ifelse(a2==0,1,a2*TOTHT^(a2-1))
		deriv_A3<-ifelse(b1==0,1,b1*(X^(Z^4))^(b1-1))
		deriv_A4<-ifelse(b2==0,1,(b2)*(X^(1/exp(DBH/TOTHT)))^(b2-1))
		deriv_A5<-ifelse(b3==0,1,(b3)*(X^(X^0.1))^(b3-1))
		deriv_A6<-ifelse(b4==0,1,(b4)*(X^(X^Q))^(b4-1))
		deriv_A7<-ifelse(b5==0,1,(b5)*(X^(TOTHT^Q))^(b5-1))
		deriv_A8<-ifelse(b6==0,1,(b6)*(X^(X))^(b6-1))
		deriv_A9<-ifelse(b7==0,1,b7*(X^CR)^(b7-1))
		
		deriv_A<-matrix(c(deriv_A0,deriv_A1,deriv_A2,deriv_A3,deriv_A4,deriv_A5,
						deriv_A6,deriv_A7,deriv_A8,deriv_A9),nrow=length(A1))
		
		grad<-matrix(rep(res,times=dim(A)[2]),nrow=length(A1))
		grad<-grad/A
		grad<-grad*deriv_A
		colnames(grad)<-c("a0","a1","a2","b1","b2","b3","b4","b5","b6","b7")
		attr(res, "gradient") <- grad
		
	}
	
	
	

	return(res)
	
}
