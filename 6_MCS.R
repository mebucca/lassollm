#  ============== Creates function that reports goodness of fit statistics after log-linear models ============


################################################################## Open Dataset ########################################################## 


Frequencies.sim <- cont.table(mydata)
Frequencies.sim <- llm.pars(Frequencies) 


################################################################ Fit LogLinear Models ############################################# 


## Creates X-matrix with predictors for saturated model

df.full <- Frequencies.sim  %>% mutate_at(-4,funs(factor(.))) %>% as_tibble()


sim.bic <- NULL  # Simulated BICs


# Selected sample sizes for reduced simulation
sz <- c(log(2),log(12),log(85))


#sample.size <- log(85)

#for (sample.size in seq(0.5,6.4,0.1)) {
for (sample.size in sz) {
	# ::::::::::::::::::::::: Creates Simulates Frequencies  :::::::::::::::::::::::

	# DGP include a contrained diagonal + decrease in homogamy among least educated 
	# + decrease in homogamy among most educated + constant assymetric movement across minor diagonals

	n = length(df.full$F)  # Size of dataset

	# Creates coefficients 

	intercept <- sample.size # This number will determine the size of the dataset (sum of all frequencies)

	Beta_Educ_H <- c(1,1.5,.8,.5)   # Simulated coefficients marginal distribution men's education
	Beta_Educ_W <- c(1,1.5,.8,.5)  # Simulated coefficients marginal distribution women's education
	Beta_Year   <- rep(0,length(table(df.full$Year))-1) #  Simulated coefficients marginal distribution years

	Beta_Educ_WH <- rep(0, (length(table(df.full$Educ_W))-1) * (length(table(df.full$Educ_H))-1) ) #  Simulated coefficients for null assortative mating
	Beta_Educ_WYear <- rep(0, (length(table(df.full$Educ_W))-1) * (length(table(df.full$Year))-1) ) #  Simulated coefficients for null educational expansion women
	Beta_Educ_HYear <- rep(0, (length(table(df.full$Educ_H))-1) * (length(table(df.full$Year))-1) ) #  Simulated coefficients for null educational expansion men

	Beta_diag <- c(0,2,0,0,2) # Simulated coefficients for Assortative Mating in unconstrained Diagonal 
	
	# Simulated coefficients for decreasing homogamy among least educated 
	Beta_diag2_Year <- NULL
	for(i in 1:(length(table(df.full$Year)))) { Beta_diag2_Year <- c(Beta_diag2_Year, - 0.1*i) } 

	# Simulated coefficients for increasing homogamy among most educated 
	Beta_diag5_Year <- NULL
	for(i in 1:(length(table(df.full$Year)))) { Beta_diag5_Year <- c(Beta_diag5_Year, + 0.1*i) } 

	# Simulated coefficients for assymetric movements along the minor diagonals 

	Beta_move_asym <- c(0.5,0.5)


	# Simulates Logged Frequencies 
	logmu <- intercept + 
	dummy(df.full$Educ_H)%*% Beta_Educ_H + 
	dummy(df.full$Educ_W)%*% Beta_Educ_W + 
	dummy(df.full$diag) %*%  Beta_diag + 
	dummy(df.full$diag:df.full$Year)[,24:(24+11)]%*% Beta_diag2_Year +
	dummy(df.full$diag:df.full$Year)[,60:(60+11)]%*% Beta_diag5_Year +
	dummy(df.full$move_asym) %*% Beta_move_asym

	# Simulates Frequencies 

	F_sim <- round(exp(logmu),0)

	# Random frequencies simulation using a poisson dgp. -> Makes it count data and adds randomness to the process

	F_sim <- rpois(n, F_sim) 

	# Add Simulated Frequencies to dataset 

	df.full$F_sim <- F_sim
	
	df.full$F_sim

	tag <- df.full %>% summarise(n = mean(F_sim)) %>% round(0)
	
	# Summary of distribution of frequencies
	
	name <-paste("Summary_Freqs_",tag, sep='')
	summary <- df.full %>% summarise(N = sum(F_sim), Nmean= mean(F_sim), min=min(F_sim), max=max(F_sim), ncell=n(), small.5 = mean(as.numeric(F_sim<5)), small.10 = mean(as.numeric(F_sim<10)) ) 
	assign(name, summary)

	# Plot distribution of frequencies
	
	name <-paste("Distrib_Freqs_",tag, sep='')
 	plot.distrib.freq <- ggplot(df.full, aes(x=F_sim)) + geom_histogram(color="red", fill="red", bins=200, size=0.05) + labs(x="cell frequencies", y="count") + xlim(0,2200) + ylim(0,50)
	assign(name, plot.distrib.freq)

	# ::::::::::::::::::::::::::::::::: Null Model ::::::::::::::::::::::::::::::::


	model.null <- glm(F_sim ~ 1, family=poisson, data=df.full)
	gof(df.full$F_sim,model.null)


	# ::::::::::::::::::::::::::::::: Saturante Model :::::::::::::::::::::::::::::


	model.sat <- glm(F_sim ~ Educ_W*Educ_H*Year, family=poisson, data=df.full)
	gof(df.full$F_sim,model.sat)

	# :::::::::::::::::::::::::::: Model of Independence :::::::::::::::::::::::::::

	model.indep <- glm(F_sim ~ Educ_W + Educ_H + Year, family=poisson, data=df.full)
	gof(df.full$F_sim,model.indep)


	# ::::::::::::::::::::: Model of *Conditional Independence :::::::::::::::::::::

	model.condindep <- glm(F_sim ~ Educ_W + Educ_H + Year + Educ_W*Educ_H + Educ_W*Year + Educ_H*Year, family=poisson, data=df.full)
	gof(df.full$F_sim,model.condindep)


	# ::::::::::::::: Model of Quasi Perfect Mobility Constrained ::::::::::::::::::

	model.qpmc <- glm(F_sim ~ Educ_W + Educ_H + Year + Educ_W*Educ_H + Educ_W*Year + Educ_H*Year + diagC + diagC*Year, family=poisson, data=df.full)
	gof(df.full$F_sim,model.qpmc)


	# :::::::::::::: Model of Quasi Perfect Mobility Unconstrained :::::::::::::::::

	model.qpm <- glm(F_sim ~ Educ_W + Educ_H + Year + Educ_W*Educ_H + Educ_W*Year + Educ_H*Year + diag + diag*Year, family=poisson, data=df.full)
	gof(df.full$F_sim,model.qpm)


	# ::: Model of  Symmetric movement across the minor diagonals + Main diagonal constrained :::

	model.sdc <- glm(F_sim ~ Educ_W + Educ_H + Year + Educ_W*Educ_H + Educ_W*Year + Educ_H*Year + diagC + diagC*Year + move_sym + move_sym*Year, family=poisson, data=df.full)
	gof(df.full$F_sim,model.sdc)


	# ::: Model of Assymmetric movement across the minor diagonals + Main diagonal constrained :::

	model.hyp <- glm(F_sim ~ Educ_W + Educ_H + Year + Educ_W*Educ_H + Educ_W*Year + Educ_H*Year + diagC + diagC*Year + move_asym + move_asym*Year, family=poisson, data=df.full)
	gof(df.full$F_sim,model.hyp)

	# ::: Model of  Symmetric movement across the minor diagonals + Main diagonal unconstrained :::

	model.sduc <- glm(F_sim ~ Educ_W + Educ_H + Year + Educ_W*Educ_H + Educ_W*Year + Educ_H*Year + diag + diag*Year + move_sym + move_sym*Year, family=poisson, data=df.full)
	gof(df.full$F_sim,model.sduc)


	# ::: Model of Assymmetric movement across the minor diagonals + Main diagonal unconstrained :::

	model.hypuc <- glm(F_sim ~ Educ_W + Educ_H + Year + Educ_W*Educ_H + Educ_W*Year + Educ_H*Year + diag + diag*Year + move_asym + move_asym*Year, family=poisson, data=df.full)
	gof(df.full$F_sim,model.hypuc)

	# ::::::::::::::::::::::::::::::::::::::: Crossing Models :::::::::::::::::::::::::::::::::::::::::::::

	model.crossing <- glm(F_sim ~ Educ_W + Educ_H + Year + Educ_W*Educ_H + Educ_W*Year + Educ_H*Year + 
	crossing_1 + crossing_2 + crossing_3 + crossing_4 + crossing_1*Year + crossing_2*Year + crossing_3*Year + crossing_4*Year,
	family=poisson, data=df.full)

	gof(df.full$F_sim,model.crossing)

	# ::::::::::::::::::::::::::::: Crossing Models + Constrained diagonal ::::::::::::::::::::::::::::::::

	model.crossingH <- glm(F_sim ~ Educ_W + Educ_H + Year + Educ_W*Educ_H + Educ_W*Year + Educ_H*Year + diagC + diagC*Year +
	crossing_1 + crossing_2 + crossing_3 + crossing_4 + crossing_1*Year + crossing_2*Year + crossing_3*Year + crossing_4*Year,
	family=poisson, data=df.full)

	gof(df.full$F_sim,model.crossingH)

	# ::::::::::::::::::::::::::::: Crossing Models + Unconstrained diagonal ::::::::::::::::::::::::::::::::

	model.crossingMD <- glm(F_sim ~ Educ_W + Educ_H + Year + Educ_W*Educ_H + Educ_W*Year + Educ_H*Year + diag + diag*Year +
	crossing_1 + crossing_2 + crossing_3 + crossing_4 + crossing_1*Year + crossing_2*Year + crossing_3*Year + crossing_4*Year,
	family=poisson, data=df.full)

	gof(df.full$F_sim,model.crossingMD)

	# ::::::::::::::::::::::::::::::::::::: ORDINAL MODELS ::::::::::::::::::::::::::::::::

	# ::::::::::::::::::::::::::::::::::::: Linear by Linear Association ::::::::::::::::::::::::::::::::

	model.linear <- glm(F_sim ~ Year + Educ_W + Educ_H + Educ_W*Year + Educ_H*Year + Educ_W*Educ_H + LA*Year , family=poisson, data=df.full)
	gof(df.full$F_sim,model.linear)


	# ::::::::::::::::::::::::::::::::::::: UNIDIFF model ::::::::::::::::::::::::::::::::

	model.unidiff <- gnm(F_sim ~ Year*Educ_W + Year*Educ_H + Mult(Year, Educ_W:Educ_H),
	                     family = poisson,
	                     data = df.full)
	gof.unidiff(df.full$F_sim,model.unidiff)
	#tabla de contingencia solo importa diagonal. Parametros se van a 0. Ano a otro, pais a otro. 
	#patrom AM mismo todos los amos cambia intensidad. 


	# ::::::::::::::::::::::::::::: Table with statistcis of Goodness of fit ::::::::::::::::::::::::::::::::


	models <- c("model.indep","model.condindep","model.qpmc","model.qpm","model.sdc","model.hyp","model.sduc","model.hypuc",
		        "model.crossing","model.crossingH","model.crossingMD","model.linear","model.unidiff")


	table.gof   <- NULL # Empty table

	for (i in models) {

		# Computes Goodness of fit statistics for each model
		
		if (i=="model.unidiff") {
			tgof <- gof.unidiff(df.full$F_sim, eval(parse(text = i))  )

		}
		else {
			tgof <- gof(df.full$F_sim, eval(parse(text = i))  )
		}
		

		table.gof <- rbind(table.gof,tgof)  # Piles stats for each model

	}

	rownames(table.gof) <- c(" Indep = [W][H][Y]"," CI = [W][H][Y][HW][WY][HY]","CI+[CDiag][CDiagY]",
							 "CI+[Diag][DiagY]", "CI+[CDiag][CDiagY][Sym][SymY]",
							 "CI+[CDiag][CDiagY][Asym][AsymY]",
							 "CI+[Diag][DiagY][Sym][SymY]",
							 "CI+[Diag][DiagY][Asym][AsymY]",
							 "CI+[Cross][CrossY]",
							 "CI+[Cross][CrossY][Sym][SymY]",
							 "CI+[Cross][CrossY][Asym][AsymY]","Linear by Linear","Unidiff")

	print(table.gof)



	title.table <- "Log-Linear Models of the Association Between Husband’s and Wife’s Educational Attainment, simulated data"
	lab <- paste0("tab:","Gofllm_",tag)

	print(xtable(table.gof, caption = title.table, align = "lccccc", label = lab), type="latex", floating=TRUE, table.placement="h!", hline.after=c(-1,0,13),
	size = "footnotesize", latex.environments = "center", caption.placement = "top", 
	sanitize.colnames.function = identity, file= paste0(path.output.tables,"Gofllm_",tag,".tex"))



	# ::::::::::::::::::::::  Evaluates Performance of BIC with different sample size::::::::::::::::::::::::

	sim.bic <- cbind(sim.bic,table.gof[,5]) # Extracts the column with BIC for all specifications for simulated Ns


	rownames(sim.bic) <- c(" Indep = [W][H][Y]","Indep-Year = [WY][HY]"," CI = [W][H][Y][HW][WY][HY]","CI+[CDiag][CDiagY]",
							 "CI+[Diag][DiagY]", "CI+[CDiag][CDiagY][Sym][SymY]",
							 "CI+[CDiag][CDiagY][Asym][AsymY]",
							 "CI+[Diag][DiagY][Sym][SymY]",
							 "CI+[Diag][DiagY][Asym][AsymY]",
							 "CI+[Cross][CrossY]",
							 "CI+[Cross][CrossY][Sym][SymY]",
							 "CI+[Cross][CrossY][Asym][AsymY]","N")

	print(round(sim.bic))


	# Plot evolution of BIC for some specifications

	#rows   <- c(2,3,4,5,6,9)  # index of loglinear models specifications
	#colors <- wes_palette("Darjeeling", 5) # colors
	#colors <- c(colors,"black")
	#
	#sim.bic.selected <- sim.bic[rows, ] # only keep models of
	#
	#sim.bic.select <- sim.bic.selected
	#for (row in 1:length(rows)) {
	#	sim.bic.select[row,] <- sim.bic.selected[row,]/sim.bic.selected[1,]
	#}
	#
	#
	#min = min(sim.bic.select) - .01
	#max = max(sim.bic.select) + .1
	#
	#
	#key<-paste0(path.output, "MC_BIC.pdf")
	#pdf(key)
	#
	#par(mar=c(4,4,1,1))
	#
	#plot(sim.bic[12,],sim.bic.select[1,],ylim=c(min,max),type="b", col=colors[1], cex=.6, main="", xlab="Average n by cell", ylab="BIC.m / BIC.ci")
	#
	#legend(900, max , c("CI","DiagC","Diag","Diag+Sym","Diag+Asym","Crossing"), lty=c(1,1), lwd=c(2.5,2.5), col=colors, cex=0.6, ncol=2) 
	#
	#for (i in 2:length(rows)) {
	#  lines(sim.bic[12,],sim.bic.select[i,],ylim=c(min,max),type="b", col=colors[i], ylab="", xlab="", cex=.6)
	#}
	#
	#abline(v=sum(df.full$F)/n, lty=2)
	#
	#dev.off()



################################## Implements Lasso regression for selection  of Log-Linear Models ######################################## 


	# ================================= Creates Design Matrix for Saturate Model ==================================


	# Formula that includes all variables and interactions in a given dataframe

	f1 <- as.formula(y ~ .*.*.)

	y <- df.full$F_sim

	# Dataframe with Year, Educ_H and Educ_W

	df1 <- df.full[,c(1,2,3)]

	# Creates all three-way interactions

	x1 <- model.matrix(f1, df1) # design matrix with interaction between Year, Educ_H and Educ_W
	x  <- cbind(x1) 

	# display names
	head(x,0)


	# ======================= Creates list with names corresponding to vars in saturated model ====================

	# Names for coefficients

	years <- c(seq(1992,2000,2),seq(2003,2009,3),seq(2011,2015,2))
	h.educ <- c("(H) E","(H) HS", "(H) C-", "(H) C")
	w.educ <- c("(W) E","(W) HS", "(W) C-", "(W) C")

	main.eff <- c("Intercept",years,h.educ,w.educ) # main effects


	h.educ.y <-  NULL  # husband education by year 

	for (i in h.educ ) {
		for (j in years) {
			hy <- paste0(i," * ",j)
			h.educ.y <- c(h.educ.y,hy)
		}
	}

	w.educ.y <-  NULL  # wife education by year 

	for (i in w.educ ) {
		for (j in years) {
			wy <- paste0(i," * ",j)
			w.educ.y <- c(w.educ.y,wy)
		}
	}


	w.educ.h.educ <- NULL # husband education * wife education


	for (i in w.educ ) {
		for (j in h.educ ) {
				wh <- paste0(i," * ",j)
				w.educ.h.educ  <- c(w.educ.h.educ,wh)
		}
	}


	w.educ.h.educ.y <-  NULL  # wife education * husband education, by year 


	for (i in w.educ.h.educ ) {
		for (j in years ) {
				why <- paste0(i," * ",j)
				w.educ.h.educ.y  <- c(w.educ.h.educ.y,why)
		}
	}


	labels <- c(main.eff,h.educ.y,w.educ.y,w.educ.h.educ,w.educ.h.educ.y); labels



	# ======================================== Estimates Lasso regression =========================================

	glmmod <- glmnet(x,y,alpha=1,family='poisson')
	cv.glmmod <- cv.glmnet(x,y,alpha=1,family='poisson')


	# Coefficients from 'best' lambda value 

	opt.lam     <- c(cv.glmmod$lambda.1se)
	lasso.coefs <- coef(cv.glmmod, s = opt.lam)


	# ================================ Extracts Coefficients from 'best' lambda value  ============================

	 
	opt.lam     <- c(cv.glmmod$lambda.1se)
	lasso.coefs <- coef(cv.glmmod, s = opt.lam)



	# Create predictions for margins and assortative mating. Used in plot 
	
 	new_x <- df.full %>% data_grid(Year,Educ_W,Educ_H,.model=model.hypuc) %>% 
 		select(Year,Educ_W,Educ_H) %>% model.matrix(f1, .)
	
 	intercept = lasso.coefs[1]


	# full prediction
	predictions <- cbind(df.full %>% data_grid(Year,Educ_W,Educ_H,.model=model.hypuc), predict(cv.glmmod, new_x, s=opt.lam)) %>%
		as_tibble() %>% rename(pred = `1`) %>% mutate(pred = pred - intercept)

	# margins
	predictions_year  <- predictions %>% filter(Educ_H==1, Educ_W==1)  %>% rename(margin_year=pred) %>% select(Year,margin_year)
	predictions_educh <- predictions %>% filter(Year==1990, Educ_W==1) %>% rename(margin_educh=pred) %>% select(Educ_H,margin_educh)
	predictions_educw <- predictions %>% filter(Year==1990, Educ_H==1) %>% rename(margin_educw=pred) %>% select(Educ_W,margin_educw)


	# match
	
	predictions <- predictions %>% left_join(predictions_year, by="Year")
	predictions <- predictions %>% left_join(predictions_educh, by="Educ_H")
	predictions <- predictions %>% left_join(predictions_educw, by="Educ_W")
	

	# margins by year
	predictions_year_educh <- predictions %>% filter(Educ_H!=1,Year!=1990,Educ_W==1) %>%
		rename(margin_year_educh=pred) %>% mutate(margin_year_educh = margin_year_educh - (margin_year + margin_educh )) %>%
		select(Year,Educ_H,margin_year_educh)

	predictions_year_educw <- predictions %>% filter(Educ_H==1,Year!=1990,Educ_W!=1) %>%
		rename(margin_year_educw=pred) %>% mutate(margin_year_educw = margin_year_educw - (margin_year + margin_educw )) %>% 
		select(Year,Educ_W,margin_year_educw)

	predictions <- predictions %>% left_join(predictions_year_educh, by=c("Year","Educ_H")) %>%  replace_na(list(margin_year_educh = 0))
	predictions <- predictions %>% left_join(predictions_year_educw, by=c("Year","Educ_W")) %>%  replace_na(list(margin_year_educw = 0))


	# computes assortative mating
	predictions <- predictions %>% 
		mutate(assort = pred - (margin_year + margin_educh + margin_educw + margin_year_educh + margin_year_educw) ) %>% select(-diag,-move_asym)


	# =========================================== Plots ================================================

	# Plot Cross-Validation error

	name <-paste("MC_CV_error_",tag, sep='')

	pdf(NULL)
	dev.control(displaylist="enable")
	par(mar=c(8,8,8,8))
	plot(cv.glmmod)
	abline(v=log(opt.lam), lty=2)
	text(log(opt.lam)+0.22,280,expression(paste(lambda,"*")))
	plot.mc.cv.error <- recordPlot()
	invisible(dev.off())

	assign(name, plot.mc.cv.error)


	# Plot coefficients path 

	name <-paste("MC_coef_path_",tag, sep='')
	pdf(NULL)
	dev.control(displaylist="enable")
	par(mar=c(6,6,6,6))
	plot(glmmod,xvar="lambda")
	abline(v=log(opt.lam), lty=2)
	text(log(opt.lam)+0.15,3.2,expression(paste(lambda,"*")))
	plot.mc.coef.path <- recordPlot()
	invisible(dev.off())
	
	assign(name, plot.mc.coef.path)


	# Plots predicted log-Odds


	name <-paste("MC_Lasso_free_logOdds_",tag, sep='')

	plot.assort.lasso <- 
		predictions  %>%
		ggplot(aes(x=factor(Year), y=assort, group=factor(Educ_H) )) +
		facet_grid( . ~ Educ_W, labeller=as_labeller(educ.labeller) ) +
		geom_line(aes(colour=factor(Educ_H, labels = educ.labels)), alpha=0.5) + 
		geom_point( aes(colour=factor(Educ_H, labels = educ.labels), shape=factor(Educ_H, labels = educ.labels)), size=1.5) + 
		theme(axis.text.x = element_text(angle=90, vjust=0.5, size=8)) + 		
		labs(subtitle= "Wive's Educ",  y="Log-Odds", x="", shape="Husband's Educ", colour= "Husband's Educ", title = "") +
		ylim(-2.5, 7.5) +
		scale_color_d3() + scale_fill_d3()


	assign(name, plot.assort.lasso)

	# ====================================== Simulated Confidence intervals ===========================================


	# Iterated cross-validation to choose the "best" lamda value. (Only works with more than 200 simulations)


	sims = n_sims
	iter = 1
	Betas = matrix(NA,dim(x)[2],sims)
	Betas.CI = matrix(NA,dim(x)[2],3)

	lambda.set <- NULL # Will store all lambda values chosen 

	while (iter <= sims) {

		cv.glmmod   <- cv.glmnet(x,y,alpha=1,family='poisson')
		opt.lam <- c(cv.glmmod$lambda.1se)
		lasso.coefs <- coef(cv.glmmod, s = opt.lam)
		lasso.coefs <- as.matrix(lasso.coefs[-2])

		lambda.set <- c(lambda.set,opt.lam) # saves lambda value

		Betas[,iter]<- lasso.coefs	

		iter = iter + 1
	}

#	# Plot stability of coefficients across iterations 


#	key<-paste(path.output.images,"MC_Lasso_simulation_",tag,".pdf", sep='')
#	pdf(key)
#
#	lambdas <- data.frame(lambda.set)
#	min = min(Betas)
#	max = max(Betas) + 4
#
#	plot(Betas[1,], ylim=c(min,max), col="white", ylab = "Coefficient", xlab="Iteration")
#
#
	for (i in c(1, 2:dim(x)[2])) {
		#lines(Betas[i,], col= alpha("blue", 0.2) )
		Betas.CI[i,1] = mean(Betas[i, ]) # Estimator: mean across different solutions in cross-validation
		Betas.CI[i,2] = quantile(Betas[i,],0.25)
		Betas.CI[i,3] = quantile(Betas[i,],0.975)
	}
#		
#		# Add distribution of lamda to the same plot. 
#
#
#	plot.lambdas <- ggplot(lambdas, aes(x=lambda.set)) + geom_histogram(binwidth=.05, colour="red", fill="red") + xlab("Lambda") 
#
#	plot.lambdas <- plot.lambdas + theme(text = element_text(size=7))
#
#	print(plot.lambdas, vp=viewport(.75, .75, .2, .2))
#
#	dev.off()


	# Creates table with lasso coefficients, simulated mean lasso coefficients and simulated confidence intervals

	rownames(Betas.CI) <- labels

	cv.lasso.coefs <- round(Betas.CI,4)

	# Combines the two tables

	lasso.coefs <- cbind(lasso.coefs,cv.lasso.coefs)

	colnames(lasso.coefs) <- c("Coeff.", "CV Mean Coeff","Q(0.025)","Q(.975)")
	row.names(lasso.coefs) <- labels

#	#Save Table with coefficientsplot coefficients
#
	title.table <- paste0("Lasso Estimates, Saturated Model for Association Between Husband’s and Wife’s Educational Attainment 1990-2015, Simulated Data (N=",sum(df.full$F_sim),")")

	#Trick for long table 

	addtorow          <- list()
	addtorow$pos      <- list()
	addtorow$pos[[1]] <- c(0)
	addtorow$command  <- c(paste("\\hline \n",
	                            "\\endhead \n",
	                            "\\hline \n",
	                            "{\\footnotesize Continued on next page} \n",
	                            "\\endfoot \n",
	                            "\\endlastfoot \n",sep=""))

	lab <- paste0("tab:","MC_Lasso_free_",tag)

	print(xtable(lasso.coefs, caption = title.table, align = "lcccc", label = lab), type="latex", floating=FALSE, table.placement="H", 
	size = "scriptsize", tabular.environment="longtable", latex.environments = "center", caption.placement = "top",
	include.rownames = TRUE, add.to.row = addtorow, hline.after=c(-1,0,dim(lasso.coefs)[1]-1),
	sanitize.colnames.function = identity, file= paste0(path.output.tables,"MC_Lasso_free_",tag,".tex"))


}