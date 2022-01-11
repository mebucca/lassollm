set.seed(5506)

# ============================================== SATURATED MODEL ============================================== 

# ================================= Creates Design Matrix for Saturate Model ==================================

options(contrasts=c("contr.treatment","contr.poly"))

# Formula that includes all variables and interactions in a given dataframe

f <- as.formula(~ .*.*.)

y <- Frequencies$F

# Dataframe with Year, Educ_H and Educ_W

df <- Frequencies[,c(1,3,2)]

# Creates all three-way interactions

x <- model.matrix(f, df) # design matrix with interaction between Year, Educ_H and Educ_W


# ======================= Creates list with names corresponding to vars in saturated model ====================

# Names for coefficients

years <- names(table(mydata$Year))[-1]
h.educ <- c("(H) E","(H) HS", "(H) C-", "(H) C")
w.educ <- c("(W) E","(W) HS", "(W) C-", "(W) C")

main.eff <- c("Intercept",years,w.educ,h.educ) # main effects


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


for (i in h.educ ) {
	for (j in w.educ ) {
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


labels <- c(main.eff,w.educ.y,h.educ.y,w.educ.h.educ,w.educ.h.educ.y); labels




# Defines length of variables

delta.years     <- length(names(table(mydata$Year)))
delta.educ.h    <- length(names(table(mydata$Educ_H)))
delta.educ.w    <- length(names(table(mydata$Educ_W)))

# ===============================================  Lasso regression ==========================================================


# Set penalty factors: Here we define the penality factor of two models. The first one doesnt impose 
# restrictions of the parameters to be penalized (lasso free), while the second one leave the marginal 
# distributions unpenalized (lasso indep)


#########  lasso free ######

p.fac.free = rep(1,length(y))


######### Marginals are not penalized ##############

# Non regularized
position.ini = 1
position.fin = 1 + (delta.years - 1) + (delta.educ.w - 1) + (delta.educ.h - 1) 
p.fac.nonreg = rep(0,length(seq(position.ini,position.fin)))  

# Regularized

position.ini = position.fin + 1
position.fin = length(y)
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
position.fin = length(y)
p.fac.reg = rep(1,length(seq(position.ini,position.fin)))  

# penalty factors

p.fac.hywy = c(p.fac.nonreg,p.fac.reg) 


######### Marginals and Two-way interactions are not penalized ##############

# Non regularized
position.ini = 1
position.fin = 1 + (delta.years - 1) + (delta.educ.w - 1) + (delta.educ.h - 1) + (delta.years - 1)*(delta.educ.h - 1) + (delta.years - 1)*(delta.educ.w - 1) + (delta.educ.w - 1)*(delta.educ.h - 1)
p.fac.nonreg = rep(0,length(seq(position.ini,position.fin)))  

# Regularized

position.ini = position.fin + 1
position.fin = length(y)
p.fac.reg = rep(1,length(seq(position.ini,position.fin)))  

# penalty factors

p.fac.condindep = c(p.fac.nonreg,p.fac.reg) 


######### Adaptive Lasso ##############


coefs.sat <- as.matrix(model.sat$coefficients)
p.fac.adaptive <- 1/pmax(abs(coefs.sat), .Machine$double.eps)^(1/4)


# ================================ Estimates models and Extracts Coefficients from 'best' lambda value  ============================


penalty      <- c("free","adaptive","indep","hywy")
titles.lasso <- c("Lasso Free", "Adaptive Lasso", "Lasso Indep.", "Lasso Indep. + [WY][HY]")

Lassotype  <- 1

for (Lassotype in 1:length(penalty)) {


	# Corresponding penalty factor
	penal.fact <- eval(parse(text = paste0("p.fac.",penalty[Lassotype])))

	# Estimates model 
	 
	glmmod<-glmnet(x,y,alpha=1,family='poisson', penalty.factor = penal.fact)


	# Determines best lambda with cross-validations

	cv.glmmod <- cv.glmnet(x,y,alpha=1,family='poisson', penalty.factor = penal.fact)
	
	 
	opt.lam     <- c(cv.glmmod$lambda.1se)
	lasso.coefs <- coef(cv.glmmod, s = opt.lam)


	# Create predictions for margins and assortative mating. Used in plot 
	
 	new_x <- Frequencies %>% data_grid(Year,Educ_W,Educ_H,.model=model.hypuc) %>% 
 		select(Year,Educ_W,Educ_H) %>% model.matrix(f, .)

	
 	intercept = lasso.coefs[1]

 	# full prediction
	predictions <- cbind(Frequencies %>% data_grid(Year,Educ_W,Educ_H,.model=model.hypuc), predict(cv.glmmod, new_x, s=opt.lam)) %>%
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

	name <- paste0("plot_CVerror_Lasso",penalty[Lassotype])

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

	name <- paste0("plot_Coefpath_Lasso",penalty[Lassotype])

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

	# Calculates log-Odds of assortative mating 


	name <- paste0("plot_logOdds_Lasso",penalty[Lassotype])

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

	#save_plot(key, plot.assort.lasso, base_height =5, base_width=12)

	# ====================================== Simulated Confidence intervals ===========================================


	# Iterated cross-validation to choose the "best" lamda value. (Only works with more than 200 simulations)

	sims = n_sims
	iter = 1
	Betas = matrix(NA,dim(x)[2],sims)
	Betas.CI = matrix(NA,dim(x)[2],3)

	lambda.set <- NULL # Will store all lambda values chosen 

	while (iter <= sims) {

		cv.glmmod <- cv.glmnet(x,y,alpha=1,family='poisson', penalty.factor = penal.fact)
		opt.lam <- c(cv.glmmod$lambda.1se)
		lasso.coefs <- coef(cv.glmmod, s = opt.lam)
		lasso.coefs <- as.matrix(lasso.coefs[-2])

		lambda.set <- c(lambda.set,opt.lam) # saves lambda value

		Betas[,iter]<- lasso.coefs	

		iter = iter + 1
	}

#	# Plot stability of coefficients across iterations 
#
#
#	key<-paste0(path.output.tables, "SimCoefs_Lasso",penalty[Lassotype],".pdf")
#	pdf(key)
#
#	lambdas <- data.frame(lambda.set)
#	min = min(Betas)
#	max = max(Betas) + 4
#
#	plot(Betas[1,], ylim=c(min,max), col="white", ylab = "Coefficient", xlab="Iteration")
#
	for (i in c(1, 2:dim(x)[2])) {
#		lines(Betas[i,], col= alpha("blue", 0.2) )
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

	title.table <- " Lasso Estimates, Saturated Model for Association Between Husband’s and Wife’s Educational Attainment (Husbands Aged 30-35): Chile, 1990-2015"

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

	lab <- paste0("tab:","Lasso_",penalty[Lassotype],"_chile")

	print(xtable(lasso.coefs, caption = title.table, align = "lcccc", label = lab), type="latex", floating=FALSE, table.placement="H", 
	size = "scriptsize", tabular.environment="longtable", latex.environments = "center", caption.placement = "top",
	include.rownames = TRUE, add.to.row = addtorow, hline.after=c(-1,0,dim(lasso.coefs)[1]-1),
	sanitize.colnames.function = identity, file= paste0(path.output.tables,"TabCoefs_Lasso",penalty[Lassotype],".tex"))

}



