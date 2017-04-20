# TODO: Add comment
# 
# Author: Paco
###############################################################################

library(nlme)
library(sae)
library(Matrix)

create_sample_XZ<-function(sample.data,fixed_pred,variances,groupslme){
	
	nlevels<-length(levels(sample.data[,groupslme]))
	glevels<-levels(sample.data[,groupslme])
	
	
#	Cambiar para incluir efectos mas complicados
	Z<-matrix(rep(0,times=length(sample.data[,1])*nlevels*(length(variances)-1)),nrow=length(sample.data[,1]),
			ncol=length(levels(sample.data[,groupslme]))*(length(variances)-1))
	
	colnames(Z)<-paste(names(variances)[-length(variances)],rep(levels(sample.data[,groupslme]),each=length(variances)-1),sep="_")
	
	
	for(row in 1:length(sample.data[,1])){
		
		col<-levels(sample.data[row,groupslme])[sample.data[row,groupslme]]
		col<-which(col==levels(sample.data[row,groupslme]))
		
		for(i in 1:(length(variances)-1)){
			
			if(names(variances)[i]=="(Intercept)"){
				
				Z[row,(col-1)*(length(variances)-1)+i]<-1
				
			}else{
				
				Z[row,(col-1)*(length(variances)-1)+i]<-sample.data[row,names(variances)[i]]
				
			}
			
			
			
		}		
		
	}
	
	Z<-Matrix(Z, sparse = TRUE)
	
	X<-sample.data[,fixed_pred]
	namesX<-c("(Intercept)",fixed_pred)
	X<-as.matrix(cbind(1,X))
	colnames(X)<-namesX
	X<-Matrix(X)
	c(list(Z=Z),list(X=X))
	
}

create_Gs<-function(variances,levels){
	
#	Cambiar esta funcion para incluir mas random effects usar la funcion create
#	gblocks y revisar la ordenacion de los efectos
	
	
	G<-variances[1]*diag(length(levels))
	dGv<-diag(length(levels))
	dGe<-0*diag(length(levels))
	
	list(G=G,dGv=dGv,dGe=dGe)
	
}

create_var_Matrices<-function(variances,levels,Z,lme.obj){
	
#	esto habrá que cambiarlo si se incluyen mas efectos aleatorios

	Gs<-create_Gs(variances,levels)	
	
	if(is.null(lme.obj$modelStruct$varStruct)){
		
		R<-Matrix(variances[length(variances)]*diag(nrow(Z)), sparse = TRUE)
		dRv<-Matrix(0*diag(nrow(Z)), sparse = TRUE)
		dRe<-Matrix(diag(nrow(Z)), sparse = TRUE)
		
		
	}else{
		
		if(is(lme.obj$modelStruct$varStruct,"varPower")){
			
			
			formula<-attributes(lme.obj$modelStruct$varStruct)$formula
			
			cov_name<-all.vars(formula)
			
			cov<-lme.obj$data[,cov_name]
			exp<-attributes(lme.obj$modelStruct$varStruct)$fixed
			cov<-cov^(2*exp)
			
			R<-Matrix(variances[length(variances)]*diag(cov), sparse = TRUE)
			dRv<-Matrix(0*diag(nrow(Z)), sparse = TRUE)
			dRe<-Matrix(diag(cov), sparse = TRUE)
			
			
		}else{
			
			
			formula<-attributes(lme.obj$modelStruct$varStruct)$formula
			
			cov_name<-all.vars(formula)
			
			cov<-lme.obj$data[,cov_name]
			
			R<-Matrix(variances[length(variances)]*diag(cov), sparse = TRUE)
			dRv<-0*diag(nrow(Z))
			dRe<-Matrix(diag(cov), sparse = TRUE)
			
		}
		
	}
	
	V<-R+Z%*%Gs$G%*%t(Z)
	dVv<-Z%*%Gs$dGv%*%t(Z)+dRv
	dVe<-Z%*%Gs$dGe%*%t(Z)+dRe
	Vinv<-solve(V)
	
	c(list(V=V,dVv=dVv,dVe=dVe),Gs,list(R=R,dRv=dRv,dRe=dRe),list(Vinv=Vinv))
	
	
}

#matrices con los coeficientes m,l,q para los que se quiere estimar
create_MLQ<-function(groupslme,predict.data,fixed_pred,random_pred,lme.obj){
	
	grouplist<-list(predict.data[,groupslme])
	names(grouplist)<-groupslme
	
	nlevels<-length(levels(predict.data[,groupslme]))
	glevels<-levels(predict.data[,groupslme])
	
	vars<-names(predict.data)
	predict.data<-cbind(predict.data,1)
	names(predict.data)<-c(vars,"(Intercept)")
	
	rows<-dim(predict.data)[1]
	
	Lp<-as.matrix(cbind(1,predict.data[,fixed_pred]))
	colnames(Lp)<-c("(Intercept)",fixed_pred)
	

	
	Lg1<-predict.data[,c("Area",fixed_pred)]
	Lg1[,-1]<-Lg1[,-1]*Lg1[,1]
	
	Lt<-apply(Lg1,2,sum)
	Lt<-Lt/Lt[1]
	Lt<-matrix(Lt,nrow=1)
	colnames(Lt)<-c("(Intercept)",fixed_pred)
	
	Lg1<-aggregate(Lg1,by=grouplist,sum)
	present_groups<-levels(Lg1[,groupslme])[Lg1[,groupslme]]
	colnames(Lg1)<-c(groupslme,"(Intercept)",fixed_pred)
	
	Area_t<-0

	cont<-0
	for(group in 1:nlevels){
		
		lev_group<-glevels[group]				
		if(lev_group%in%present_groups){
			
#				Aqui hay que cambiar por un apply y subirlo fuera del bucle,
#				Si no va a ser muy lento
			cont<-cont+1
			if(cont==1){
				
				Lgrupo<-Lg1[Lg1[,groupslme]==lev_group,-1]
				Lgdeni<-Lg1[Lg1[,groupslme]==lev_group,"(Intercept)"]
				Lg<-Lgrupo/Lgdeni
				
				
			}else{
				
				Lgrupo<-Lg1[Lg1[,groupslme]==lev_group,-1]
				Lgdeni<-Lg1[Lg1[,groupslme]==lev_group,"(Intercept)"]
				Lgrupo<-Lgrupo/Lgdeni
				Lg<-rbind(Lg,Lgrupo)
				
				
			}

		}

	}	
	
	L<-as.matrix(rbind(Lp,Lg,Lt))
	L<-t(L)
	
#	Creación de M
	
#	Creacion de Mp
	
	Mp<-matrix(rep(0,nlevels*length(random_pred)*length(Lp[,1])),
			nrow=length(Lp[,1]))
	
	new_names<-names(random_pred)
	ind_groups<-which(new_names==groupslme)
	random_pred[ind_groups]<-"(Intercept)"
	new_names[ind_groups]<-("(Intercept)")
	names(random_pred)<-new_names
	
	colnames(Mp)<-paste(new_names,rep(glevels,
					each=length(new_names)),sep="_")
	
	Mg<-matrix(rep(0,nlevels*length(random_pred)*length(present_groups)),
			nrow=length(present_groups))
	
	colnames(Mg)<-colnames(Mp)
	
	for(i in 1:length(random_pred)){
		
		v<-predict.data[,random_pred[i]]*predict.data[,"Area"]
		
		for(group in 1:nlevels){
			
			lev_group<-glevels[group]
			Mp[,(length(random_pred)*(group-1)+i)]<-
					ifelse(lev_group==glevels[predict.data[,groupslme]],
							v,0)
			
		}
		
		
	}
	
#	Relleno Mg 
	vg<-apply(Mp,2,sum)
	Area_t<-c(0)
	row_Mg<-0
	row_names_Mg<-c()
	for(group in 1:nlevels){
		
		lev_group<-glevels[group]				
		if(lev_group%in%present_groups){
			
			row_Mg<-row_Mg+1
			row_names_Mg<-c(row_names_Mg,lev_group)
			for(i in 1:length(random_pred)){
				
				if(i==1){
					
					
#				Aqui hay que cambiar por un apply y subirlo fuera del bucle,
#				Si no va a ser muy lento
					Area_g<-sum(predict.data[glevels[predict.data[,groupslme]]==lev_group,
									"Area"])
					Area_t<-Area_t+Area_g
					
				}
				
				
				
				Mg[row_Mg,(length(random_pred)*(group-1)+i)]<-vg[(length(random_pred)*(group-1)+i)]/Area_g
				
				
			}
			
			
		}
			
		
	}
	
#	pongo nombres de filas en Mg
	rownames(Mg)<-row_names_Mg
	
#	Creo Mt
	Mt<-vg/Area_t
#	Escalo Mp
	Mp<-Mp/predict.data[,"Area"]
	
	M<-as.matrix(rbind(Mp,Mg,Mt))
	M<-t(M)
	
#	Ojo!!, No es exactamente Q, si no que incluye los pesos de la seigma_e

	if(is.null(lme.obj$modelStruct$varStruct)){
		
		
		Q<-rep(1,length(predict.data[,1]))
		
		Qg<-aggregate((predict.data[,"Area"]^2),by=grouplist,sum)
		Qgden<-aggregate(predict.data[,"Area"],by=grouplist,sum)
		
		
		
	}else{
		
		if(is(lme.obj$modelStruct$varStruct,"varPower")){
			
			
			formula<-attributes(lme.obj$modelStruct$varStruct)$formula
			exp<-attributes(lme.obj$modelStruct$varStruct)$fixed
			
			cov_name<-all.vars(formula)
			cov<-predict.data[,cov_name]
			
			Q<-cov^(2*exp)
			Qg<-aggregate(((predict.data[,cov_name]^(2*exp))*(predict.data[,"Area"]^2)),by=grouplist,sum)
			Qgden<-aggregate(predict.data[,"Area"],by=grouplist,sum)
			

			
		}else{
			
			
			formula<-attributes(lme.obj$modelStruct$varStruct)$formula
			
			cov_name<-all.vars(formula)
			
			cov<-predict.data[,cov_name]
			Q<-cov
			
			Qg<-aggregate((predict.data[,cov_name]*(predict.data[,"Area"]^2)),by=grouplist,sum)
			Qgden<-aggregate(predict.data[,"Area"],by=grouplist,sum)
			
			
		}
		
		
		
	}

	Area_t<-0
	for(group in 1:nlevels){
		
		lev_group<-glevels[group]				
		if(lev_group%in%present_groups){
			
#				Aqui hay que cambiar por un apply y subirlo fuera del bucle,
#				Si no va a ser muy lento

				Qgi<-Qg[Qg[,groupslme]==lev_group,2]
				Qgdeni<-Qgden[Qg[,groupslme]==lev_group,2]
				Area_t<-Area_t+Qgdeni
				Q<-c(Q,Qgi/(Qgdeni^2))
				
				
			}
			
			
		}
	
	Qt<-sum(Q[1:length(predict.data[,"Area"])]*(predict.data[,"Area"]^2))
	Qt<-Qt/((sum(predict.data[,"Area"]))^2)
	Q<-c(Q,Qt)
	Q<-matrix(Q,nrow=1)
	
	M<-Matrix(M, sparse = TRUE)
	L<-Matrix(L, sparse = TRUE)
	Q<-Matrix(Q)
	
	indexes<-list(ip=1:length(Lp[,1]),ig=length(Lp[,1])+c(1:length(Lg[,1])),it=dim(L)[2])
	list(M=M,L=L,Q=Q,indexes=indexes)
	
}

create_BD<-function(M,L,var_Matrices,X,Z,random_pred){
	
	B<-Matrix(t(t(M)%*%var_Matrices$G%*%t(Z)%*%var_Matrices$Vinv))
	
#	esta parte hay que cambiarla y hacerla como una lista con los nombres de las
#	variables
	dBv<-Matrix(t(t(M)%*%var_Matrices$dGv%*%t(Z)%*%var_Matrices$Vinv-t(M)%*%var_Matrices$G%*%t(Z)%*%var_Matrices$Vinv%*%var_Matrices$dVv%*%var_Matrices$Vinv))
	dBe<-Matrix(t(t(M)%*%var_Matrices$dGe%*%t(Z)%*%var_Matrices$Vinv-t(M)%*%var_Matrices$G%*%t(Z)%*%var_Matrices$Vinv%*%var_Matrices$dVe%*%var_Matrices$Vinv))

	D<-t(t(L)-t(B)%*%X)
#	D<-t(D)
	
	list(D=D,B=B,dBv=dBv,dBe=dBe)
	
}

create_predictions<-function(lme.obj,predict.data,groupslme,B,L){
	
	#	Calculo de los EBLUPS
	
	nlevels<-length(levels(predict.data[,groupslme]))
	glevels<-levels(predict.data[,groupslme])
	
	grouplist<-list(predict.data[,groupslme])
	names(grouplist)<-groupslme
	
	mat_groups_aux<-predict.data[,"Area"]
	mat_groups_aux<-aggregate(mat_groups_aux,by=grouplist,sum)
	present_groups<-levels(mat_groups_aux[,groupslme])[mat_groups_aux[,groupslme]]
	
	if(dim(predict.data)[1]==1){
		
		predict.data2<-rbind(predict.data,predict.data)
		eblups<-predict(lme.obj,predict.data2,0:1)
		eblups<-eblups[1,]
		eblups[,groupslme]<-factor(levels(eblups[,groupslme])[eblups[,groupslme]],levels=levels(predict.data[,groupslme]))
		
	}else{
		
		eblups<-predict(lme.obj,predict.data,0:1)
		eblups[,groupslme]<-factor(levels(eblups[,groupslme])[eblups[,groupslme]],levels=levels(predict.data[,groupslme]))
				
	}
	
	predict.data<-cbind(predict.data,eblups[,-1])
	lastcol<-dim(predict.data)[2]
	predict.data[,c(lastcol-1,lastcol)]<-predict.data[,c(lastcol-1,lastcol)]*predict.data[,"Area"]
	
	lastcols<-colnames(predict.data)[c(lastcol-1,lastcol)]
	eblups_gauxb<-aggregate(predict.data[,c("Area",lastcols)],by=grouplist,sum)
	eblups_gauxb[,c(3,4)]<-eblups_gauxb[,c(3,4)]/eblups_gauxb[,2]
	
	Area_t<-0
	row_eblup_gaux<-0
	for(group in 1:nlevels){
		
		lev_group<-glevels[group]				
		if(lev_group%in%present_groups){
			
#				Aqui hay que cambiar por un apply y subirlo fuera del bucle,
#				Si no va a ser muy lento
			row_eblup_gaux<-row_eblup_gaux+1
			
			if(row_eblup_gaux==1){
				
				eblups_gaux<-eblups_gauxb[eblups_gauxb[,groupslme]==lev_group,c(1,3,4)]

			}else{
				
				
				eblups_gauxc<-eblups_gauxb[eblups_gauxb[,groupslme]==lev_group,c(1,3,4)]
				eblups_gaux<-rbind(eblups_gaux,eblups_gauxc)
				
				
			}

			
			
		}
		
		
	}
	
	colnames(eblups_gaux)<-colnames(eblups)
	eblups<-rbind(eblups,eblups_gaux)
	
	b_tot<-B[,dim(B)[2]]
	l_tot<-L[,dim(L)[2]]
	
	fixed_tot<-sum(lme.obj$coefficients$fixed*l_tot)
	eblup_tot<-sum(fixed_tot+t(b_tot)%*%(lme.obj$residuals[,"fixed"]))
	eblups<-rbind(eblups,c(0,fixed_tot,eblup_tot))
	rownames(eblups)<-c(1:length(eblups[,1]))
	
	data.frame(eblups)
	
}

create_Fisher<-function(var_Matrices){
	
	V_minus<-var_Matrices$Vinv
#	K<-solve(Xt%*%V_minus%*%X)
#	P<-V_minus-V_minus%*%X%*%K%*%Xt%*%V_minus	
#	Aqui habrá que cambiar para introducir más efectos aleatorios Rao 2003 pp139
	V_v<-V_minus%*%var_Matrices$dVv
	V_e<-V_minus%*%var_Matrices$dVe
	Fvv<-1/2*sum(V_v*t(V_v))
	Fve<-1/2*sum(V_v*t(V_e))
	Fev<-1/2*sum(V_e*t(V_v))
	Fee<-1/2*sum(V_e*t(V_e))
	Fisher<-matrix(c(Fvv,Fve,Fev,Fee),nrow=2,byrow=TRUE)
	
	V_hat<-solve(Fisher)
	
	list(Fisher=Fisher,V_hat=V_hat)
	
}

create_G1G2G3G4<-function(M,Q,X,Z,G,V,V_inv,aux_Matrices,var_e,V_hat){
	
	D<-aux_Matrices$D
	g1<-t(M)%*%(G-G%*%t(Z)%*%V_inv%*%Z%*%G)

	g2<-t(D)%*%(solve(t(X)%*%(solve(V))%*%X))
	
	g3_v1<-t(aux_Matrices$dBv)%*%V
	g3_v1a<-rowSums(g3_v1*t(aux_Matrices$dBv))
	g3_v2b<-rowSums(g3_v1*t(aux_Matrices$dBe))
	
	g3_e1<-t(aux_Matrices$dBe)%*%V
	g3_e1a<-rowSums(g3_e1*t(aux_Matrices$dBv))
	g3_e2b<-rowSums(g3_e1*t(aux_Matrices$dBe))
	
	G1<-rowSums(g1*t(M))
	G2<-rowSums(g2*t(D))
	G3<-g3_v1a*V_hat[1,1]+g3_v2b*V_hat[2,1]+g3_e1a*V_hat[1,2]+g3_e2b*V_hat[2,2]
	G4<-Q[1,]*var_e
	
	res<-data.frame(G1=G1,G2=G2,G3=G3,G4=G4)
	
}

merge_pred_errors<-function(lme.obj,G1G2G3G4,eblups,IDs,
		groupslme,response,predict.data){
	
	grouplist<-list(predict.data[,groupslme])
	names(grouplist)<-groupslme
	Lg_MLQ<-predict.data[,"Area"]
	Lg_MLQ<-aggregate(Lg_MLQ,by=grouplist,sum)
	present_groups<-levels(Lg_MLQ[,groupslme])[Lg_MLQ[,groupslme]]
	
	sample<-lme.obj$data[,c(groupslme,response)]
	
	Ids<-c(predict.data[,IDs],paste("grupo_",present_groups,sep=""),"total")
	eblups<-cbind(eblups,g1=G1G2G3G4$G1,g2=G1G2G3G4$G2,g3=G1G2G3G4$G3,g4=G1G2G3G4$G4,Ids=Ids)
	
	eblups[is.na(eblups[,3]),3]<-eblups[is.na(eblups[,3]),2]
	
	
	eblups$MSEN<-eblups$g1+eblups$g2
	eblups$MSE1<-eblups$g1+eblups$g2+eblups$g3
	eblups$MSE2<-eblups$g1+eblups$g2+2*eblups$g3
	
	eblups$MSENp<-eblups$g1+eblups$g2+eblups$g4
	eblups$MSEp1<-eblups$g1+eblups$g2+eblups$g3+eblups$g4
	eblups$MSEp2<-eblups$g1+eblups$g2+2*eblups$g3+eblups$g4
	
	eblups$RMSE<-sqrt(eblups$MSE2)
	eblups$SE<-1.96*sqrt(eblups$MSE2)
	eblups$CV<-sqrt(eblups$MSE2)/eblups[,3]
	eblups$L95<-eblups[,3]-eblups$SE
	eblups$U95<-eblups[,3]+eblups$SE
	
	eblups$RMSEp<-sqrt(eblups$MSEp2)
	eblups$SEp<-1.96*sqrt(eblups$MSEp2)
	eblups$CVp<-sqrt(eblups$MSEp2)/eblups[,3]
	eblups$L95p<-eblups[,3]-eblups$SEp
	eblups$U95p<-eblups[,3]+eblups$SEp
	
	
	
	eblups<-list(pixels=eblups[1:length(predict.data[,1]),],stands=eblups[(length(predict.data[,1])+1):(length(eblups[,1])-1),],
			total=eblups[length(eblups[,1]),])
	
	
	
	eblups$total$mu_direct<-mean(sample[,2])
	eblups$total$SD_samp<-sd(sample[,2])/sqrt((length(sample[,2])))
	eblups$total$SE_direct<-eblups$total$SD
	

	eblups$stands$mu_direct<-NA
	eblups$stands$SD_samp<-NA
	eblups$stands$SE_direct<-NA
	eblups$stands$ni<-NA
	eblups$stands$CVd<-NA

	for(gr in present_groups){
		
		subsampl<-sample[which(sample[,groupslme]==gr),]
		mui<-mean(subsampl[,2])
		ni<-length(subsampl[,2])
		sdi<-sd(subsampl[,2])/sqrt((length(subsampl[,2])))
		sei<-1.96*sdi
		eblups$stands[which(eblups$stands[,groupslme]==gr),"mu_direct"]<-mui
		eblups$stands[which(eblups$stands[,groupslme]==gr),"SD_samp"]<-sdi
		eblups$stands[which(eblups$stands[,groupslme]==gr),"SE_direct"]<-sei
		eblups$stands[which(eblups$stands[,groupslme]==gr),"ni"]<-ni
		eblups$stands[which(eblups$stands[,groupslme]==gr),"CVd"]<-sdi/mui

				
	}
	
	eblups$stands$delta_eff<-(eblups$stands$CVd-eblups$stands$CV)/eblups$stands$CVd
	eblups$stands$delta_effp<-(eblups$stands$CVd-eblups$stands$CVp)/eblups$stands$CVd
	
	
	eblups
	
	
}

eblups<-function(lme.obj,predict.data,IDs){
	
	response<-all.vars(attr(terms(lme.obj),"variables"))[
			attr(terms(lme.obj),"response")]

#	Obtencion de los grupos de toda la población
	
	groupslme<-colnames(lme.obj$groups)
	
	predict.data[,groupslme]<-as.factor(predict.data[,groupslme])
	
	totallevels<-levels(factor(c(levels(lme.obj$groups[,groupslme]),levels(predict.data[,groupslme]))))
	
	sample.data<-lme.obj$data
	sample.data[,groupslme]<-as.factor(factor(levels(lme.obj$groups[,groupslme])[lme.obj$groups[,groupslme]],levels=totallevels))
	predict.data[,groupslme]<-as.factor(factor(levels(predict.data[,groupslme])[predict.data[,groupslme]],levels=totallevels))
	
	variances<-VarCorr(lme.obj)[,"Variance"]
	names_variances<-names(variances)
	variances<-as.numeric(variances)
	names(variances)<-names_variances
	
#	Aquí hay que retocar
	fixed_pred<-names(lme.obj$coefficients$fixed)[-which(names(lme.obj$coefficients$fixed)=="(Intercept)")]
	random_pred<-names(lme.obj$coefficients$random)
	names(random_pred)<-random_pred
	
	if(all(random_pred==groupslme)){
		
		random_pred<-NULL
		vars<-c(IDs,fixed_pred,"Area")
		
	}else{
		
		random_pred<-random_pred[-which(random_pred==groupslme)]
		vars<-c(IDs,fixed_pred,random_pred,"Area")
	}
	
#	Datos predicciones por pixel
	predict.data<-predict.data[,c(groupslme,vars)]
	
#	Cambio otra vez random_pred para que incluya los indices de grupo ID_SMA
	random_pred<-as.character(names(lme.obj$coefficients$random))
	names(random_pred)<-random_pred
	
	sample_XZ<-create_sample_XZ(sample.data,fixed_pred,variances,groupslme)
	
	var_Matrices<-create_var_Matrices(variances,totallevels,sample_XZ$Z,lme.obj)
	
	est_Matrices<-create_MLQ(groupslme,predict.data,fixed_pred,random_pred,lme.obj)
	
	aux_Matrices<-create_BD(est_Matrices$M,est_Matrices$L,
			var_Matrices,sample_XZ$X,sample_XZ$Z,random_pred)
	
	Fisher<-create_Fisher(var_Matrices)
	
	eblups<-create_predictions(lme.obj,predict.data,groupslme,aux_Matrices$B,est_Matrices$L)
	
	G1G2G3G4<-create_G1G2G3G4(est_Matrices$M,est_Matrices$Q,sample_XZ$X,sample_XZ$Z,
			var_Matrices$G,var_Matrices$V,var_Matrices$Vinv,aux_Matrices,variances[length(variances)],Fisher$V_hat)
	
	eblups<-list(eblups=merge_pred_errors(lme.obj,G1G2G3G4,eblups,IDs,
					groupslme,response,predict.data),Fisher=Fisher)
	
}

prepare_data<-function(lme.obj,predict.data,ext=0.05){
	
	cols<-names(lme.obj$coefficients$fixed)
	if(any(cols=="(Intercept)")){
	
		cols<-cols[-which(cols=="(Intercept)")]
		
	}
	
	for(i in cols){
		
		rang<-range(lme.obj$data[,i],na.rm=TRUE)
		width<-rang[2]-rang[1]
		predict.data<-predict.data[predict.data[,i]>rang[1]-width*ext,]
		predict.data<-predict.data[predict.data[,i]<rang[2]+width*ext,]
		
	}
	
	
	predict.data
	
}

eblups_gls<-function(gls.obj,predict.data,groups){
	
	response<-all.vars(formula(gls.obj))[1]
	predictors<-all.vars(formula(gls.obj))[-1]
	
	lt<-data.frame(cbind(1,predict.data[,predictors]))
	colnames(lt)<-c("(Intercept)",predictors)
	
	lt_total<-as.matrix(colSums(lt*predict.data[,"Area"])/sum(predict.data["Area"]))
	
	lt_groups<-split(lt*predict.data[,"Area"],predict.data[,groups],drop=TRUE)
	lt_groups<-lapply(lt_groups,function(x){
				
				res<-as.matrix(colSums(x)/sum(x[,1]))

			})
	
	lt<-as.matrix(lt)
	
	R<-make_R(gls.obj$modelStruct$corStruct,
			gls.obj$modelStruct$varStruct,
			gls.obj$sigma,predict.data)

	VarBeta<-gls.obj$varBeta
	
	index_groups<-split(c(1:length(predict.data[,1])),predict.data[,groups],drop=TRUE)
	
	pixels<-data.frame(eblup=predict(gls.obj,predict.data),
			g1=diag(R$R),g2=rowSums((lt%*%VarBeta)*lt))
	
	stands<-aggregate(x=data.frame(eblup=pixels$eblup*predict.data[,"Area"],Area=predict.data[,"Area"]),
			by=list(stands=predict.data[,groups]),FUN=function(x){
				sum(x)
				
			})
	
	stands$eblup<-stands$eblup/stands$Area
	
	g1_groups<-lapply(index_groups,function(x,R){
				
				a<-matrix(predict.data[x,"Area"]/sum(predict.data[x,"Area"]))
				res<-sum(a*diag(R)[x]*a)
				res
				
			},R=R$R)
	
	g1_groups<-data.frame(g1=unlist(g1_groups),stands=names(g1_groups))
	
	g2_groups<-lapply(lt_groups,function(x,VarBeta){
				
				res<-t(x)%*%VarBeta%*%x
				res[1,1]
				
			},VarBeta=VarBeta)
	
	g2_groups<-data.frame(g2=unlist(g2_groups),stands=names(g2_groups))
	
	stands<-merge(stands,g1_groups,by="stands")
	stands<-merge(stands,g2_groups,by="stands")
	
	total<-data.frame(eblup=sum(pixels$eblup*predict.data[,"Area"])/sum(predict.data[,"Area"]),
			Area=sum(predict.data[,"Area"]))
	
	
	a<-matrix(predict.data[,"Area"]/sum(predict.data[,"Area"]))
	g1_total<-sum(a*R$R*a)
	total$g1<-g1_total
	
	g2_total<-t(lt_total)%*%VarBeta%*%lt_total
	total$g2<-g2_total
	
	eblups<-list(pixels=pixels,stands=stands,total=total)
	
	eblups<-lapply(eblups,function(x){
				
				x$MSEp<-x$g1+x$g2
				x$RMSEp<-sqrt(x$MSEp)
				x$CVp<-x$RMSEp/x$eblup
				x$L95p<-x$eblup-1.96*x$RMSEp
				x$U95p<-x$eblup+1.96*x$RMSEp
				x
			})
	
}

calculate_h<-function(formula,sampled_data,vardir,sigma2_v){
	
	X_sampled<-model.matrix(formula,sampled_data)

	gamma<-1/(sigma2_v+sampled_data[,vardir])
	hc<-0
	for(i in 1:length(X_sampled[,1])){
		
		hc<-hc+((X_sampled[i,]%*%t(X_sampled[i,]))*gamma[i])
		
	}
	hc<-solve(hc)
	
	
}

unsampled_mseFH<-function(formula,unsampled_data,sampled_data,vardir,coefs,sigma2_v,ni){
	
	ni<-deparse(substitute(ni))
	if(!is.character(vardir)){
		
		vardir<-deparse(substitute(vardir))	
		
	}
	if(ni%in%colnames(unsampled_data)){
		
		ni<-unsampled_data[,ni]
		
	}else{
		
		ni<-NA
		
	}
	res<-data.frame(zone=row.names(unsampled_data),ni=ni,est_dir=unsampled_data[,variable],
			vardir=NA,var_field=NA,variable=all.vars(formula)[1])
	row.names(res)<-row.names(unsampled_data)
	hc<-calculate_h(formula,sampled_data,vardir,sigma2_v)
	
	unsampled_data<-unsampled_data[,all.vars(formula)]
	unsampled_data[,1]<-1
	X_unsampled<-model.matrix(formula,unsampled_data)
	synthetic<-data.frame(eblup=X_unsampled%*%coefs)
	res<-merge(res,synthetic,by=0,all.x=TRUE)
	
	row.names(res)<-res[,"Row.names"]
	res<-res[,!colnames(res)=="Row.names"]
	
	mse_unsampled<-sigma2_v+rowSums((X_unsampled%*%hc)*X_unsampled)
	mse_unsampled<-data.frame(mse=mse_unsampled)
	res<-merge(res,mse_unsampled,by=0,all.x=TRUE)
	row.names(res)<-res[,"Row.names"]
	res<-res[,!colnames(res)=="Row.names"]
	
	res$rmse<-sqrt(res$mse)
	res$rrmse<-res$rmse/res$eblup
	res
	
}

