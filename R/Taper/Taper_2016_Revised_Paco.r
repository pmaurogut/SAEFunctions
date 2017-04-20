# ****   Variables definitions used in analysis    ****
# Study from which the data came from (OSU biomass sampling and modeling project)
# LOC			Location identification  
# TREEID	Tree number identifier
# DBH			Tree diameter at breast height (cm)
# TOTHT			Tree total height (m)
# HT			Height of sectioned disk from the ground (m)
# DOB			Diameter outside bark at a given height, HT (cm)
# DIB			Diameter inside bark at a given height, HT(cm)
# HCB			Tree height to crown base (m; lowest live branch)
# CL			Tree crown length (m)

Taper <- read.csv("Taper/Data/Sample_Trees.csv",stringsAsFactors =FALSE)

Taper <- Taper[which(Taper$DOB > 0),]  #Excludes tree tip.

#attach(Taper); names(Taper)

Taper <- Taper[which(Taper$SPP=="DF" & Taper$DOB > 0),]
what<-Taper$TREEID%in%unique(Taper$TREEID)[1:2]
new_tree<-Taper[what,]
Taper<-Taper[!what,]#Excludes tree tip.

# Definitions of the variables used in the nonlinear mixed-effects model
Taper$Z <- Taper$HT/Taper$TOTHT
Taper$P <- 1.3/Taper$TOTHT
Taper$Q <- 1 - (Taper$HT/Taper$TOTHT)^(1/3)
Taper$X <- Taper$Q / (1 - (Taper$P)^(1/3))

Taper$DH<-Taper$DBH/Taper$TOTHT
Taper$HD<-Taper$TOTHT/Taper$DBH
Taper$CL<-Taper$TOTHT-Taper$HCB
Taper$CR<-Taper$CL/Taper$TOTHT

#Taper$X<-(1-sqrt(Taper$Z))/(1-sqrt(Taper$p))    #X ranges from 0 at the tip to 1 at the reference point p.
Taper$X1<-(1-(Taper$Z)^1/4)/(1-0.01^1/4) 
Taper$X2<-(1-(Taper$Z)^1/3)/(1-Taper$P^1/3)  
Taper$Q <- 1-(Taper$Z^1/3)

library(nlme)
# fit nls for Kozak 1988;
nls1 <- nls(DIB~(a0*DBH^a1*a2^DBH)*X^(b1*Z^2+b2*log(Z+0.0001)+b3*sqrt(Z)+b4*exp(Z)+b5*(DBH/TOTHT)),
            data = Taper, 
            start=c(a0=1.02,a1=0.888, a2=1.01, b1=0.95,b2=-0.18,b3=0.61,b4=-0.35,b5=0.057))

summary(nls1)
AIC(nls1)

#Kozak 2001 model  
kozak01<-nls(DIB~((a0*DBH^a1)*X1^(b0+b1*(1/exp(DH))+b2*DBH^X1+b3*(X1^DH))),data=Taper,
                     start=c(a0=1.02,a1=1.01,b0=0.95,b1=1.18,b2=0.61,b3=0.35))
summary(kozak01)

n=length(Taper)
m=5
p01=predict(kozak01)

res01=p01-DIB
E01=mean(res01)
E01
RMSE01=sqrt(sum((res01)^2)/(n-1))
RMSE01
#DIB.mean=mean(DIB)
R201=1-((sum((res01)^2))/(sum((DIB-mean(DIB))^2)))
R201
AIC(kozak01)
BIC(kozak01) 
logLik(kozak01)

# fit nls for Kozak 2002;

kozak02 <- nls(DIB ~ (a0*(DBH^a1)*(TOTHT^a2)*X^(b1*(Z^4)+b2*(1/exp(DBH/TOTHT))+b3*(X^0.1)+
                                                  b4*(1/DBH)+b5*(TOTHT^Q)+b6*X)), data=Taper,
               start=c(a0=1.37, a1=0.89, a2=0.2, b1=0.21, b2=-0.53, b3=0.02, b4=-0.01, b5=-0.04, b6=-0.1),
               nls.control(minFactor=1e-010, maxiter = 2000))

summary(kozak02)
AIC(kozak02)
BIC(kozak02) 


#Kozak 2002 model  with CR
kozak02cr <- nls(DIB ~ (a0*(DBH^a1)*(TOTHT^a2)*X^(b1*(Z^4)+b2*(1/exp(DBH/TOTHT))+b3*(X^0.1)+
                                                    b4*(1/DBH)+b5*(TOTHT^Q)+b6*X+b7*CR)), data=Taper,
                 start=c(a0=1.37, a1=0.89, a2=0.2, b1=0.21, b2=-0.53, b3=0.02, b4=-0.01, b5=-0.04, b6=-0.1, b7=0.1),
                 nls.control(minFactor=1e-010, maxiter = 2000))
AIC(kozak02cr)
BIC(kozak02cr) 


##Arias-Rodil et al 2015 CJFR 45:647-658 paper.
	
f.nls=nls(DIB~2*((b0*DBH)/(1-exp(b2*(1.3-TOTHT)))+(DBH/2-b0*DBH)*(1-(1/(1-exp(b1*(1.3-TOTHT)))))+
  (exp(-b1*HT))*(((DBH/2-b0*DBH)*exp(1.3*b1))/(1-exp(b1*(1.3-TOTHT))))-exp(b2*HT)*((b0*DBH*exp(-b2*TOTHT))/(1-exp(b2*(1.3-TOTHT))))),
  data=Taper,start=list(b0=0.3435,b1=0.65,b2=0.09)) 
summary(f.nls)

n=length(Taper)
m=2
p=predict(f.nls)
#write.table((p),file="E://workfile//Sabbatical//Taper//MGY.txt")
res=p-DIB
E=mean(res)
E
RMSE=sqrt(sum((res)^2)/(n-1))
RMSE
#DIB.mean=mean(DIB)
R2=1-((sum((res)^2))/(sum((DIB-mean(DIB))^2)))
R2
AIC(f.nls)
BIC(f.nls)
logLik(f.nls)

fit_stat <- cbind(AIC(nls1),AIC(kozak01),AIC(kozak02),AIC(kozak02cr),AIC(f.nls))


                      
                      

	