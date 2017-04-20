# TODO: Add comment
# 
# Author: ursus
###############################################################################
library(nlme)
library(reshape2)
library(Matrix)

cor2var<-function(x){
	
	std<-attr(x,"stdDev")
	res<-std%*%t(std)
	res<-res*x
	
}

addstdDev<-function(x,y){
	
	names<-names(x)
	
	for(i in names){
		
		attr(x[[i]],"stdDev")<-y[[i]]
		
	}
	
	x
	
}

arrange<-function(to_split,groups){
	
	index<-1:dim(as.data.frame(to_split))[1]
	times<-dim(as.data.frame(to_split))[2]
	index<-split(index,groups,drop=TRUE)
	index<-lapply(index,function(x)diag(x))
	index<-bdiag(index)
	index<-diag(index)
	index_pos<-order(index)
	
	res<-split(to_split,groups)
	res<-lapply(res,function(x){matrix(x,ncol=times,byrow=FALSE)})
	groups_data<-names(res)
	
	res<-bdiag(res)

	res<-res[index_pos,]
	attr(res,"groups_data")<-groups_data
	res
	
}

clean_Z<-function(coefficients_random,listZ){
	
	to_remove<-c()
	last<-0
	for(i in 1:length(listZ$listZ)){
		
		coefs_names<-dimnames(coefficients_random[[i]])[[2]]
		groups_model_i<-dimnames(coefficients_random[[i]])[[1]]
		
		groups_data_i<-attr(listZ$listZ[[i]],"groups_data")
		
		groups_data_i<-rep(groups_data_i,each=length(coefs_names))
		groups_model_i<-rep(groups_model_i,each=length(coefs_names))
		
		to_remove_i<-which(groups_data_i%in%groups_model_i)
		to_remove_i<-to_remove_i+last
		to_remove<-c(to_remove,to_remove_i)
		
		last<-last+length(groups_data_i)
		
	}
	
	if(length(to_remove)>0){
		
		column_names<-attr(listZ$Z,"column_names")[-to_remove]
		group_names<-attr(listZ$Z,"group_names")[-to_remove]
		group_levels<-attr(listZ$Z,"group_levels")[-to_remove]
		group_level_names<-attr(listZ$Z,"group_level_names")[-to_remove]
		
		Z_clean<-listZ$Z[,-to_remove]
		
		attr(Z_clean,"column_names")<-column_names
		attr(Z_clean,"group_names")<-group_names
		attr(Z_clean,"group_levels")<-group_levels
		attr(Z_clean,"group_level_names")<-group_level_names
		
		listZ$Z@Dimnames[[1]]<-as.character(1:dim(listZ$Z)[1])
		listZ$Z@Dimnames[[2]]<-as.character(column_names)
		res<-list(cleanZ=Z_clean,to_remove=to_remove)
		
	}else{
		
		res<-list(cleanZ=listZ$Z,to_remove=to_remove)
		
	}
	
	res
	
}

break_estimates<-function(estimates){
	
	data<-data.frame(values=c(estimates),
			grouping=attr(estimates,"group_level_names"),
			coefficient_name=attr(estimates,"column_names"),
			level=attr(estimates,"group_names"),
			stringsAsFactors=FALSE)
	
	random_eff<-split(data,data$grouping)
}

make_Z<-function(gradient,coefficients_random,new_data,clean=TRUE){
	
	grps<-names(coefficients_random)
	listZ<-list()
	Z<-c()
	prev_group<-c()
	cont<-1
	column_names<-c()
	group_names<-c()
	group_levels<-c()
	group_level_names<-c()
	for(i in 1:length(grps)){
		
		coefs_names<-dimnames(coefficients_random[[i]])[[2]]
		
		if(is.null(gradient)){
			
			if("(Intercept)"%in%coefs_names){
				
				coefs_names<-coefs_names[-which(coefs_names=="(Intercept)")]
				if(length(coefs_names)>0){
					
					Zi<-new_data[,coefs_names]
					Zi<-cbind(1,Zi)
					colnames(Zi)<-c("(Intercept)",coefs_names)
					
				}else{
					
					Zi<-data.frame(v1=rep(1,dim(new_data[,1])))
					colnames(Zi)<-"(Intercept)"
					
				}
				
			}else{
				
				Zi<-new_data[,coefs_names]
				
			}

			Zi<-as.matrix(Zi)
			
		}else{
			
			Zi<-gradient[,coefs_names]
			
		}
		
		cols<-grps[1:i]
		to_interact<-new_data[,cols]
		level_groups<-interaction(to_interact,sep="/",drop=TRUE)
		Zi<-arrange(Zi,level_groups)
		
		groups_data_i<-attr(Zi,"groups_data")
		column_names<-c(column_names,rep(coefs_names,times=length(groups_data_i)))
		groups_data_i<-rep(groups_data_i,each=length(coefs_names))
		group_names<-c(group_names,groups_data_i)
		group_levels<-c(group_levels,rep(i,length(groups_data_i)))
		group_level_names<-c(group_level_names,rep(grps[i],length(groups_data_i)))
		
		if(cont==1){
			
			Z<-Zi
			cont<-cont+1
			
		}else{
			
			Z<-cbind2(Z,Zi)
			
		}
		
		listZ<-c(listZ,Zi)
		
	}
	
	names(listZ)<-grps
	
	listZ<-list(listZ=listZ,Z=Z,to_remove=c())
	attr(listZ$Z,"column_names")<-column_names
	attr(listZ$Z,"group_names")<-group_names
	attr(listZ$Z,"group_levels")<-group_levels
	attr(listZ$Z,"group_level_names")<-group_level_names

	listZ$Z@Dimnames[[1]]<-as.character(1:dim(listZ$Z)[1])
	listZ$Z@Dimnames[[2]]<-as.character(column_names)
	
	if(clean){
		
		cleanZ<-clean_Z(coefficients_random,listZ)
		listZ$Z<-cleanZ$cleanZ
		listZ$to_remove<-cleanZ$to_remove
	}
	
	listZ
}

make_G<-function(reStruct,coefficients_random,new_data,clean=TRUE){
	
#	Random effects matrices
	G<-reStruct(reStruct,
			formula(reStruct),data=new_data)
	G<-corMatrix(G)
	G<-lapply(G,function(x)cor2var(x))
	
	group_names<-names(G)
	
	for(i in 1:length(group_names)){
		
		level_groups<-interaction(new_data[,group_names[c(1:i)]],sep="/",drop=TRUE)
		
		if(clean){
			
			levels_model<-dimnames(coefficients_random[[i]])[[1]]
			level_groups<-level_groups[!level_groups%in%levels_model]
			
		}
		
		list_names<-c()
		if(length(level_groups)==0){
				
			G<-G[-i]
					
		}else{
			
			G[[group_names[i]]]<-bdiag(rep(G[group_names[i]],length(unique(level_groups))))
			list_names<-c(list_names,group_names[i])
		}
		
	}
	
	listG<-G
	SparseG<-bdiag(G)
	names(listG)<-list_names
	list(G=SparseG,listG=listG)
	
}

make_R<-function(corStruct,varStruct,sigma,new_data){

	if(is.null(varStruct)){
		
#	Variance of residuals
		new_data$weightsmodel<-sigma
		weightsmodel<-new_data$weightsmodel
		
	}else{
		
#	Variance of residuals
		weightsmodel<-varStruct
		weightsmodel<-Initialize(weightsmodel,new_data)
		weightsmodel<-1/varWeights(weightsmodel)
#		sigma<-sigma
		new_data$weightsmodel<-weightsmodel
		
	}
	
#	Covariances of residuals
	if(is.null(corStruct)){
		
		cor<-diag(dim(new_data)[1])
		weightmodel<-weightsmodel*weightsmodel
		R<-diag(weightsmodel)
		listR<-list(R=R)
		sparseR<-bdiag(R)
		index<-c(1:length(new_data[,1]))
		index_pos<-index

	}else{
		
		corformula<-formula(corStruct)
		grps<-getGroups(new_data,corformula)
		new_data$cor<-getCovariate(new_data,form=corformula)
		covariates<-split(new_data$cor,grps,drop=TRUE)
		cor<-corMatrix(corStruct,covariates)
		
		#	Indexes
		index<-split(c(1:length(new_data[,1])),grps,drop=TRUE)
		index<-lapply(index,function(x)diag(x))
		index<-bdiag(index)
		index<-diag(index)
		index_pos<-order(index)
		
		weightsmodel<-split(new_data$weightsmodel,grps,drop=TRUE)
		
		cor<-addstdDev(cor,weightsmodel)
		listR<-lapply(cor,cor2var)
		
#	Variance covariance matrix combination
		sparseR<-bdiag(listR)
		sparseR<-sparseR[index_pos,]
		sparseR<-sparseR[,index_pos]
		
	}
	
	list(R=sparseR,listR=listR,cor=cor,
			index=index,index_pos=index_pos,sigma=sigma)
	
}

init_new_data<-function(coefficients_random,reStructmodel,new_data){
	
	new_data[,"grlevel"]<-0
	group_names<-rev(names(reStruct(reStructmodel,
			formula(reStructmodel),data=new_data)))
	
	pars.r.i<-c()
	labels_all<-c()
	groups_ind<-c()
	for(i in 1:length(group_names)){
		
		level_groups<-interaction(new_data[,group_names[c(1:i)]],sep="/",drop=TRUE)
		new_data[,paste("Lev_",group_names[i],sep="")]<-level_groups

		groups_model_i<-dimnames(coefficients_random[[group_names[i]]])[[1]]
		to_sum<-as.character(level_groups)%in%groups_model_i
		new_data[,"grlevel"]<-new_data[,"grlevel"]+to_sum

	}
	
	new_data
	
}

init_pars<-function(coefficients_random,new_data){
	
	grps<-names(coefficients_random)
	prev_group<-c()
	cont<-1
	column_names<-c()
	group_names<-c()
	group_levels<-c()
	group_level_names<-c()
	pars.r.i<-c()
	to_remove<-c()
	last<-0
	for(i in 1:length(grps)){

		coefs_names<-dimnames(coefficients_random[[i]])[[2]]
		cols<-grps[1:i]
		to_interact<-new_data[,cols]
		level_groups<-interaction(to_interact,sep="/",drop=TRUE)
		pars.r.ib<-matrix(rep(1,2*length(new_data[,1])),ncol=2)
		pars.r.ib<-arrange(pars.r.ib,level_groups)
		
		groups_data_i<-attr(pars.r.ib,"groups_data")
		groups_model_i<-dimnames(coefficients_random[[i]])[[1]]
		groups_model_i<-rep(groups_model_i,each=length(coefs_names))
		
		column_names<-c(column_names,rep(coefs_names,times=length(groups_data_i)))
		groups_data_i<-rep(groups_data_i,each=length(coefs_names))
		
		group_names<-c(group_names,groups_data_i)
		group_levels<-c(group_levels,rep(i,length(groups_data_i)))
		group_level_names<-c(group_level_names,rep(grps[i],length(groups_data_i)))
		pars.r.i<-c(pars.r.i,rep(0,length(groups_data_i)))
		
		to_remove_i<-which(groups_data_i%in%groups_model_i)
		to_remove_i<-to_remove_i+last
		to_remove<-c(to_remove,to_remove_i)
		
		last<-last+length(groups_data_i)
		
	}
	
	if(length(to_remove)>0){

		pars.r.i<-pars.r.i[-to_remove]
		
		attr(pars.r.i,"column_names")<-column_names[-to_remove]
		attr(pars.r.i,"group_names")<-group_names[-to_remove]
		attr(pars.r.i,"group_levels")<-group_levels[-to_remove]
		attr(pars.r.i,"group_level_names")<-group_level_names[-to_remove]
		
	}else{
		
		attr(pars.r.i,"column_names")<-column_names
		attr(pars.r.i,"group_names")<-group_names
		attr(pars.r.i,"group_levels")<-group_levels
		attr(pars.r.i,"group_level_names")<-group_level_names
		
	}
	
	pars.r.i

}

make_to_eval<-function(coefs,pars.r.i,to_eval){
	
	list_random_eff<-break_estimates(pars.r.i)
	previous<-c()
	for(i in names(list_random_eff)){
		
		previous<-c(previous,i)
		
		to_eval[,i]<-as.character(interaction(to_eval[,previous],sep="/",drop=TRUE))
		
		df_i<-list_random_eff[[i]]
		
		cols_to_add<-as.character(unique(df_i$coefficient_name))
		
		cols_new_names<-paste("est_",unique(df_i$coefficient_name),sep="")
		
		df_i<-reshape2::dcast(df_i,...~coefficient_name,value.var="values")
		colnames(df_i)[which(colnames(df_i)%in%cols_to_add)]<-cols_new_names
		
		new_to_eval<-merge(to_eval,df_i,by.x=paste("Lev_",i,sep=""),by.y="level",all.x=TRUE)
		for(i in cols_new_names){
			
			new_to_eval[is.na(new_to_eval[,i]),i]<-0
			
		}
		new_to_eval[,cols_to_add]<-new_to_eval[,cols_to_add]+new_to_eval[,cols_new_names]
		
		to_eval[,cols_to_add]<-new_to_eval[,cols_to_add]
		
	}
	
	to_eval
	
}

makeCovMatrices<-function(modelStruct,sigma,coefficients_random,new_data,clean=TRUE){
	
	G<-make_G(modelStruct$reStruct,coefficients_random,new_data)
	
	R<-make_R(modelStruct$corStruct,modelStruct$varStruct,sigma,new_data)
	res<-c(G,R)
	res

}
