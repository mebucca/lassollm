
# Object that will save Deviance Scores, MSE and Poisson Deviances 

DS  <- NULL
MSE <- NULL
PD  <- NULL

DoF <- NULL # Degrees of freedom

# Object that will contain all predictions

Pred <- NULL

iter = 1

for (run in 1:10) {

	# Take Training set and creates k-training.data

	n <- nrow(mydata) # length of training set
	k <- 10 # number of folds

	# Creates folds
	fold <- round(runif(n,1,k+1),0)  # create random training.data(k)
	fold[fold==k+1] <- 1

	# Add folds to dataset
	mydata$fold <- fold


	# Loops over all folds 

	for (fld in 1:k) {

		print(paste0("===================================  ", fld, "  ======================================="))

		training  <- mydata[mydata$fold!=fld,]  # training set of Training data 
		testing   <- mydata[mydata$fold==fld,]  # tests set of Training data

		# Creates training contingency table from trainin microdata 

		Frequencies.training <- cont.table(training)
		Frequencies.training <- llm.pars(Frequencies.training) 

		#  =================================== Estimates different LogLinear models on training data ============================== 

		# ::::::::::::::::::::::::::::::::: Null Model ::::::::::::::::::::::::::::::::


		model.train.null <- glm(F ~ 1, family=poisson, data=Frequencies.training)

		# ::::::::::::::::::::::::::::::: Saturante Model :::::::::::::::::::::::::::::


		model.train.sat <- glm(F ~ Year*Educ_H*Educ_W, family=poisson, data=Frequencies.training)

		# :::::::::::::::::::::::::::: Model of Independence :::::::::::::::::::::::::::

		model.train.indep <- glm(F ~ Year + Educ_H + Educ_W, family=poisson, data=Frequencies.training)


		# ::::::::::::::::::::: Model of *Conditional Independence :::::::::::::::::::::

		model.train.condindep <- glm(F ~ Year + Educ_H + Educ_W + Educ_H*Year + Educ_W*Year + Educ_W*Educ_H, family=poisson, data=Frequencies.training)


		# ::::::::::::::: Model of Quasi Perfect Mobility Constrained ::::::::::::::::::

		model.train.qpmc <- glm(F ~ Year + Educ_H + Educ_W + Educ_H*Year + Educ_W*Year + Educ_W*Educ_H + diagC + diagC*Year, family=poisson, data=Frequencies.training)


		# :::::::::::::: Model of Quasi Perfect Mobility Unconstrained :::::::::::::::::

		model.train.qpm <- glm(F ~ Year + Educ_H + Educ_W + Educ_H*Year + Educ_W*Year + Educ_W*Educ_H + diag + diag*Year, family=poisson, data=Frequencies.training)


		# ::: Model of  Symmetric movement across the minor diagonals + Main diagonal constrained :::

		model.train.sdc <- glm(F ~ Year + Educ_H + Educ_W + Educ_H*Year + Educ_W*Year + Educ_W*Educ_H + diagC + diagC*Year + move_sym + move_sym*Year, family=poisson, data=Frequencies.training)


		# ::: Model of Assymmetric movement across the minor diagonals + Main diagonal constrained :::

		model.train.hyp <- glm(F ~ Year + Educ_H + Educ_W + Educ_H*Year + Educ_W*Year + Educ_W*Educ_H + diagC + diagC*Year + move_asym + move_asym*Year, family=poisson, data=Frequencies.training)

		# ::: Model of  Symmetric movement across the minor diagonals + Main diagonal unconstrained :::

		model.train.sduc <- glm(F ~ Year + Educ_H + Educ_W + Educ_H*Year + Educ_W*Year + Educ_W*Educ_H + diag + diag*Year + move_sym + move_sym*Year, family=poisson, data=Frequencies.training)

		# ::: Model of Assymmetric movement across the minor diagonals + Main diagonal unconstrained :::

		model.train.hypuc <- glm(F ~ Year + Educ_H + Educ_W + Educ_H*Year + Educ_W*Year + Educ_W*Educ_H + diag + diag*Year + move_asym + move_asym*Year, family=poisson, data=Frequencies.training)

		# ::::::::::::::::::::::::::::::::::::::: Crossing Models :::::::::::::::::::::::::::::::::::::::::::::

		model.train.crossing <- glm(F ~ Year + Educ_H + Educ_W + Educ_H*Year + Educ_W*Year + Educ_W*Educ_H + 
		crossing_1 + crossing_2 + crossing_3 + crossing_4 + crossing_1*Year + crossing_2*Year + crossing_3*Year + crossing_4*Year,
		family=poisson, data=Frequencies.training)


		# ::::::::::::::::::::::::::::: Crossing Models + Constrained diagonal ::::::::::::::::::::::::::::::::

		model.train.crossingH <- glm(F ~ Year + Educ_H + Educ_W + Educ_H*Year + Educ_W*Year + Educ_W*Educ_H + diagC + diagC*Year +
		crossing_1 + crossing_2 + crossing_3 + crossing_4 + crossing_1*Year + crossing_2*Year + crossing_3*Year + crossing_4*Year,
		family=poisson, data=Frequencies.training)

		# ::::::::::::::::::::::::::::: Crossing Models + Unconstrained diagonal ::::::::::::::::::::::::::::::::

		model.train.crossingMD <- glm(F ~ Year + Educ_H + Educ_W + Educ_H*Year + Educ_W*Year + Educ_W*Educ_H + diag + diag*Year +
		crossing_1 + crossing_2 + crossing_3 + crossing_4 + crossing_1*Year + crossing_2*Year + crossing_3*Year + crossing_4*Year,
		family=poisson, data=Frequencies.training)


		# ::::::::::::::::::::::::::::::::::::: ORDINAL MODELS ::::::::::::::::::::::::::::::::

		# ::::::::::::::::::::::::::::::::::::: Linear by Linear Association ::::::::::::::::::::::::::::::::

		model.train.linear <- glm(F ~ Year + Educ_W + Educ_H + Educ_W*Year + Educ_H*Year + Educ_W*Educ_H + LA*Year , family=poisson, 
			data=Frequencies.training)


		# ::::::::::::::::::::::::::::::::::::: UNIDIFF model ::::::::::::::::::::::::::::::::

		model.train.unidiff <- gnm(F ~ Year*Educ_W + Year*Educ_H + Mult(Year, Educ_W:Educ_H),
		                     family = poisson,
		                     data = Frequencies.training)
		#gof(Frequencies$F,model.unidiff)
		#tabla de contingencia solo importa diagonal. Parametros se van a 0. Ano a otro, pais a otro. 
		#patrom AM mismo todos los amos cambia intensidad. 


		#  ====================================== Estimates different Lasso model on training data ================================ 


		# Formula that includes all variables and interactions in a given dataframe

		f <- as.formula(~ .*.*.)

		y.training <- Frequencies.training$F

		# Dataframe with Year, Educ_H and Educ_W

		df <- Frequencies.training[,c(1,2,3)]

		# Creates all three-way interactions

		x.training <- model.matrix(f, df) # design matrix with interaction between Year, Educ_H and Educ_W


		# Set penalty factors

		# Defines length of variables

		delta.years  <- length(names(table(mydata$Year)))
		delta.educ.h <- length(names(table(mydata$Educ_H)))
		delta.educ.w <- length(names(table(mydata$Educ_W)))

		######### Marginals are not penalized ##############
		
		# Non regularized
		position.ini = 1
		position.fin = 1 + (delta.years - 1) + (delta.educ.w - 1) + (delta.educ.h - 1) 
		p.fac.nonreg = rep(0,length(seq(position.ini,position.fin)))  

		# Regularized

		position.ini = position.fin + 1
		position.fin = length(y.training)
		p.fac.reg = rep(1,length(seq(position.ini,position.fin)))  

		# penalty factors

		p.fac.indep = c(p.fac.nonreg,p.fac.reg) 


		######### Marginals and Marginals*Year are not penalized ##############
		
		# Non regularized
		position.ini = 1
		position.fin = 1 + (delta.years - 1) + (delta.educ.w - 1) + (delta.educ.h - 1) + (delta.years - 1)*(delta.educ.h - 1) + (delta.years - 1)*(delta.educ.w - 1) 
		p.fac.nonreg = rep(0,length(seq(position.ini,position.fin)))  

		# Regularized

		position.ini = position.fin + 1
		position.fin = length(y.training)
		p.fac.reg = rep(1,length(seq(position.ini,position.fin)))  

		# penalty factors

		p.fac.hywy = c(p.fac.nonreg,p.fac.reg) 


		######### Adaptive Lasso ##############

	    coefs.sat <- as.matrix(model.train.sat$coefficients)
	    p.fac.adaptive <- 1/pmax(abs(coefs.sat), .Machine$double.eps)^(1/2)


		######### Marginals and Two-way interactions are not penalized ##############
		
		# Non regularized
		position.ini = 1
		position.fin = 1 + (delta.years - 1) + (delta.educ.w - 1) + (delta.educ.h - 1) + (delta.years - 1)*(delta.educ.h - 1) + (delta.years - 1)*(delta.educ.w - 1) + (delta.educ.w - 1)*(delta.educ.h - 1)
		p.fac.nonreg = rep(0,length(seq(position.ini,position.fin)))  

		# Regularized

		position.ini = position.fin + 1
		position.fin = length(y.training)
		p.fac.reg = rep(1,length(seq(position.ini,position.fin)))  

		# penalty factors

		p.fac.condindep = c(p.fac.nonreg,p.fac.reg) 



		# Estimates lasso regression and determines best lambda with cross-validations


		# lasso free
		model.train.lasso.free <- cv.glmnet(x.training,y.training,alpha=1,family='poisson')
		opt.lam.free           <- c(model.train.lasso.free$lambda.1se)

		df.model.train.lasso.free <-  length(y.training) - length(which(as.matrix(coef(model.train.lasso.free, s = opt.lam.free))!=0))
		print(df.model.train.lasso.free)
		print("Hello 0")


		# lasso adaptive
		model.train.lasso.adaptive <- cv.glmnet(x.training,y.training,alpha=1,family='poisson', penalty.factor = p.fac.adaptive)
		opt.lam.adaptive           <- c(model.train.lasso.adaptive$lambda.1se)

		df.model.train.lasso.adaptive <-  length(y.training) - length(which(as.matrix(coef(model.train.lasso.adaptive, s = opt.lam.adaptive))!=0))
		print(df.model.train.lasso.adaptive)
		print("Hello 1")


		# lasso with constained penalty factors, margins 
		model.train.lasso.indep <- cv.glmnet(x.training,y.training,alpha=1,family='poisson', penalty.factor = p.fac.indep)
		opt.lam.indep           <- c(model.train.lasso.indep$lambda.1se)

		df.model.train.lasso.indep <-  length(y.training) - length(which(as.matrix(coef(model.train.lasso.indep, s = opt.lam.indep ))!=0))
		print(df.model.train.lasso.indep)

		print("Hello 2")


		# lasso with constained penalty factors, margins + margins*year
		model.train.lasso.hywy <- cv.glmnet(x.training,y.training,alpha=1,family='poisson',penalty.factor = p.fac.hywy)
		opt.lam.hywy           <- c(model.train.lasso.hywy$lambda.1se)

		df.model.train.lasso.hywy <-  length(y.training) - length(which(as.matrix(coef(model.train.lasso.hywy, s = opt.lam.hywy))!=0))
		print(df.model.train.lasso.hywy)
		print("Hello 3")


	#	# lasso with constained penalty factors, margins + two-way interections
	#	
	#	model.train.lasso.condindep <- try(cv.glmnet(x.training,y.training,alpha=1,family='poisson',penalty.factor = p.fac.condindep))
	#	
	#	# In case of error
	#	while (inherits(model.train.lasso.condindep, "try-error")==TRUE) {
	#		model.train.lasso.condindep <- try(cv.glmnet(x.training,y.training,alpha=1,family='poisson',penalty.factor = p.fac.condindep))
	#	}
	#
	#	opt.lam.condindep <- c(model.train.lasso.condindep$lambda.1se)
	#
	#	df.model.train.lasso.condindep <-  length(y.training) - length(which(as.matrix(coef(model.train.lasso.condindep, s = opt.lam.condindep))!=0))
	#	print(df.model.train.lasso.condindep)
	#	print("Hello 4")


		#  =========================== Predicts of Testing data  using models fit in training data ============================== 


		# Creates testing contingency table from testing microdata 

		Frequencies.testing <- cont.table(testing)
		Frequencies.testing <- llm.pars(Frequencies.testing) 

		# Creates Design  matrix for lasso

		# Formula that includes all variables and interactions in a given dataframe

		f <- as.formula(~ .*.*.)

		# Dataframe with Year, Educ_H and Educ_W

		df <- Frequencies.testing[,c(1,2,3)]

		# Creates all three-way interactions

		x.testing <- model.matrix(f, df) # design matrix with interaction between Year, Educ_H and Educ_W


		# Creates predictions from all models

		models <- c("model.train.indep","model.train.condindep","model.train.qpmc","model.train.qpm","model.train.sdc",
			"model.train.hyp","model.train.sduc","model.train.hypuc","model.train.crossing","model.train.crossingH",
			"model.train.crossingMD","model.train.linear","model.train.unidiff",
			"model.train.lasso.free","model.train.lasso.adaptive","model.train.lasso.indep","model.train.lasso.hywy")
		predictions <- log(Frequencies.testing$F) # Starts with vector of observed frequencies

		# Objects that will contain deviance score and poisson deviance of each model at current iteration
		
		Deviance.Score   <- NULL
		Mse              <- NULL
		Poisson.Deviance <- NULL
		Dof              <- NULL


		for (i in models) {

			if (i %in% c("model.train.lasso.free","model.train.lasso.adaptive","model.train.lasso.indep","model.train.lasso.hywy")) {

				# Creates predicted counts for lasso models

				lambda <- eval(parse(text = paste0(i,"$lambda.1se")))

				prediction.logcounts <- predict( eval(parse(text = i)), newx=x.testing, type = "link", s = lambda)
				assign(paste0("logFhat_", i), prediction.logcounts)
				predictions <- cbind(predictions, eval(parse(text = paste0("logFhat_", i)))) # dataframe with predictions from all models

				# Deviance Score
				ds <- deviance.score(exp(prediction.logcounts),Frequencies.testing$F)
				Deviance.Score <- c(Deviance.Score,ds)

				# MSE
				mse <- Mean.SE(log(Frequencies.testing$F),prediction.logcounts)
				Mse <- c(Mse,mse)

				# Poisson Deviance
				pd <- poisson.dev(Frequencies.testing$F, exp(prediction.logcounts))
				Poisson.Deviance <- c(Poisson.Deviance,pd)

				# Degrees of freedom
				dof <- eval(parse(text = paste0("df.",i) ))
				Dof <- c(Dof,dof)

			}

			else {
				# Creates predicted counts for each loglinear model
				prediction.logcounts <- predict( eval(parse(text = i)), newdata=Frequencies.testing)
				assign(paste0("logFhat_", i), prediction.logcounts)
				predictions <- cbind(predictions, eval(parse(text = paste0("logFhat_", i)))) # dataframe with predictions from all models

				# Deviance Score
				ds <- deviance.score(exp(prediction.logcounts),Frequencies.testing$F)
				Deviance.Score <- c(Deviance.Score,ds)

				# MSE
				mse <- Mean.SE(log(Frequencies.testing$F),prediction.logcounts)
				Mse <- c(Mse,mse)

				# Poisson Deviance
				pd <- poisson.dev(Frequencies.testing$F, exp(prediction.logcounts))
				Poisson.Deviance <- c(Poisson.Deviance,pd)
			
				# Degrees of freedom
				dof <- eval(parse(text = paste0(i,"$df.residual") ))
				Dof <- c(Dof,dof)
			
			}
		}

	# Add predictions from current iteration

	Pred <- rbind(Pred,predictions)
	print(paste0("rows = ",nrow(Pred)))


	# Add deviance scores, MSE and poisson deviances from current iteration
	DS  <-  cbind(DS,Deviance.Score)
	MSE <-  cbind(MSE,Mse)
	PD  <-  cbind(PD,Poisson.Deviance)
	DoF <-  cbind(DoF,Dof)

	}

}


# Compares Deviance Score and Possion Deviances of all models

DS.mean   <- as.matrix(apply(DS,1,mean))
MSE.mean  <- as.matrix(apply(MSE,1,mean))
PD.mean   <- as.matrix(apply(PD,1,mean))
DoF.mean  <- as.matrix(apply(DoF,1,mean))



CV.metrics <- cbind(PD.mean,DoF.mean) # combine matrices

colnames(CV.metrics) <- c("Poisson Deviance","Average df")


names.models <- c(" Indep = [W][H][Y]"," CI = [W][H][Y][HW][WY][HY]","CI+[CDiag][CDiagY]",
						 "CI+[Diag][DiagY]", "CI+[CDiag][CDiagY][Sym][SymY]",
						 "CI+[CDiag][CDiagY][Asym][AsymY]",
						 "CI+[Diag][DiagY][Sym][SymY]",
						 "CI+[Diag][DiagY][Asym][AsymY]",
						 "CI+[Cross][CrossY]",
						 "CI+[Cross][CrossY][Sym][SymY]",
						 "CI+[Cross][CrossY][Asym][AsymY]","Linear by Linear","Unidiff",
						 "Lasso Free","Lasso Adaptive","Lasso Indep", "Lasso Indep+[WY][HY]")

rownames(CV.metrics) <- names.models


# Saves tables as a latex table

title.table <- paste0("Error Metrics for competing models from ",iter," iterations of ",k,"-fold Cross-Validation")

lab <- paste0("tab:","CV_errormetrics_chile")

print(xtable(CV.metrics, caption = title.table, align = "lcc", label=lab), type="latex", floating=TRUE, table.placement="h!", hline.after=c(-1,0,17),
size = "footnotesize", latex.environments = "center", caption.placement = "top",
sanitize.colnames.function = identity, file= paste0(path.output.tables,"CV_errormetrics.tex"))


