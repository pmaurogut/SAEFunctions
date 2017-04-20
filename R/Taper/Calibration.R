# TODO: Add comment
# 
# Author: ursus
###############################################################################
source("Taper/Scripts/R_Scripts/Calibration/Auxiliary_New.R")
DoAnIterationCalibration <- function(model_call,coefs,new_data,Z,G,R,pars.r.i,modelStruct,coefficients_random){
	
	call_1<-model_call$model[[3]]
	response<-as.character(model_call$model[[2]])
	response<-new_data[,response]
	
	to_eval<-cbind(coefs,new_data)
	to_eval<-make_to_eval(coefs,pars.r.i,to_eval)
	fxBb <-eval(call_1,to_eval)
	
	grad <- attr(fxBb,"gradient")
	
	Zlist<-make_Z(grad,coefficients_random,new_data)
	Z<-Zlist$Z
	
	b <- G %*% t(Z) %*% solve (R + Z %*% G %*% t(Z)) %*% ((response - fxBb) + Z %*% pars.r.i)
	b<-as.vector(b)
	attributes(b)<-attributes(pars.r.i)
	return(b)
	
}
# Function for calibrating by expanding three fixed parameters
CalibrateRandomEffects <- function(model,new_data,tolerance = 1e-06, maxit = 2000){
	
	matrices<-makeCovMatrices(model,new_data)
	group_names<-names(matrices$listG)
	
	model_call<-model$call
	model_call$model[[3]]$gradient<-TRUE
	
	parms.f <- model$coefficients$fixed
	coefsfixed<-c(as.list(parms.f))
	new_data<-init_new_data(model$coefficients$random,model$modelStruct$reStruct,new_data)

	pars.r.i<-init_pars(model$coefficients$random,new_data)
	to_eval<-cbind(as.data.frame(coefsfixed),new_data)
	to_eval<-make_to_eval(coefsfixed,pars.r.i,to_eval)
	
	pred<-eval(model_call$model[[3]],to_eval)
	grad<-attr(pred,"gradient")
	
	Zlist<-make_Z(grad,model$coefficients$random,new_data)
	Z<-Zlist$Z
	R<-matrices$R
	G<-matrices$G
	
	tol <- rep(1, length(pars.r.i))
	it <- 1
	while (all(tol > tolerance) | !all(is.finite(pars.r.i))){
		
		
		if(!all(is.finite(pars.r.i))){
		
			pars.r.i[!is.finite(pars.r.i)]<-0
			
		}
		
		parms.r.f <- pars.r.i
		
		pars.r.i <- DoAnIterationCalibration(model_call,coefsfixed,
				new_data,Z,G,R,pars.r.i,model$modelStruct,model$coefficients$random)
		
		it <- it + 1
		
		tol <- abs((pars.r.i - parms.r.f) / pars.r.i)
		
		if(it > maxit){
			
			stop("Reached maximum number of iterations without convergence")
			
		}

	}
	
	pars.r.i.df<-cbind(as.data.frame(coefsfixed),new_data)
	pars.r.i.df<-make_to_eval(coefs_fixed,pars.r.i,pars.r.i.df)
	new_col_names<-colnames(pars.r.i.df)
	
	call_1<-model_call$model[[3]]
	prediction <-eval(call_1,to_eval)
	
	response<-as.character(model_call$model[[2]])
	response<-c(paste("hat_",response,sep=""))
	
	new_col_names<-c(new_col_names,response)
	
	pars.r.i.df<-cbind(pars.r.i.df,prediction)
	colnames(pars.r.i.df)<-new_col_names
	
	pars.r.i<-break_estimates(pars.r.i)
	list(pars.r.i=pars.r.i,pars.r.i.df=pars.r.i.df)
	
}

