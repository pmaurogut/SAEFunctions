# TODO: Add comment
# 
# Author: ursus
###############################################################################
#set working directory in the folder where you put Taper
#e.g. if taper is in c:\myfolder\work\Taper, then do sewd("c:/myfolder/work")
source("Taper/Scripts/R_Scripts/Calibration/Models.R")
source("Taper/Scripts/R_Scripts/Calibration/Calibration.R")

Taper <- read.csv("Taper/Data/Sample_Trees.csv",stringsAsFactors =FALSE)

Taper <- Taper[which(Taper$DOB > 0),]  #Excludes tree tip.

Taper$DH<-Taper$DBH/Taper$TOTHT
Taper$HD<-Taper$TOTHT/Taper$DBH
Taper$CL<-Taper$TOTHT-Taper$HCB
Taper$CR<-Taper$CL/Taper$TOTHT

Taper <- Taper[which(Taper$SPP=="DF" & Taper$DOB > 0 & !is.na(Taper$DOB)),]
what<-Taper$TREEID%in%unique(Taper$TREEID)[1:2]|Taper$LOC=="MD"
new_tree<-Taper[what,]
Taper<-Taper[!what,]
Taper<-Taper[!is.na(Taper$DOB),]

model_lin0<-gls(DIB ~DBH+HT,data=Taper,
		weights=varPower(0.2,form=~ HT),
		correlation=corAR1(value=0.2))

model0<-nls(DIB ~ Kozak_02b(DBH,HT,TOTHT,CR,a0,a1,a2,b1,b2,b3,b4,b5,b6,b7,gradient=FALSE),
		data=Taper,
		start=c(a0=1.37, a1=0.89, a2=0.2, b1=0.21, b2=-0.53, b3=0.02,
				b4=-0.01, b5=-0.04, b6=-0.1, b7=0.1),
		nls.control(minFactor=1e-010, maxiter = 2000))

model1<-nlme(DIB ~ Kozak_02b(DBH,HT,TOTHT,CR,a0,a1,a2,b1,b2,b3,b4,b5,b6,b7,gradient=FALSE),
		data=Taper,fixed=a0+a1+a2+b1+b2+b3+b4+b5+b6+b7~1,
		random=list(LOC=a1~1,TREEID=a0+a1+b2~1),
		weights=varPower(0.2,form=~ HT),
		correlation=corCAR1(0.2,form = ~ HT|LOC/TREEID),
		start=model0$m$getPars(),
		control=nlmeControl(pnlsTol=100,maxIter = 2000,msVerbose=TRUE))
summary(model1)

b<-CalibrateRandomEffects(model1,new_tree)


model_lin<-lme(DIB~DOB+HT,data=Taper,
		random=list(LOC=~1+HT,TREEID=~1+DOB+HT),
		correlation=corCAR1(0.2,form = ~ DOB|LOC/TREEID),
		weights=varPower(0.2,form=~ DOB),
		control=lmeControl(opt="optim",maxIter = 20000,msVerbose=TRUE))


model2<-nlme(DIB ~  Kozak_02b(DBH,HT,TOTHT,CR,a0,a1,a2,b1,b2,b3,b4,b5,b6,b7,gradient=FALSE),
		data=Taper,fixed=a0+a1+a2+b1+b2+b3+b4+b5+b6+b7~1,
		random=list(LOC=a2+a1~1,TREEID=a2+a3~1),
		start=model0$m$getPars(),
		)
summary(model2)

b<-CalibrateRandomEffects(model2,new_tree)


