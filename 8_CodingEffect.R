


####### Effect/Deviation coding ########

options(contrasts=c("contr.sum","contr.poly"))

# Formula that includes all variables and interactions in a given dataframe

seed <- round(runif(1,1,100000),0)
set.seed(seed)

f <- as.formula(~ .*.*.)

y <- Frequencies$F

# Dataframe with Year, Educ_H and Educ_W

df <- Frequencies[,c(1,3,2)]

# Creates all three-way interactions

x.effect <- model.matrix(f, df) # design matrix with interaction between Year, Educ_H and Educ_W


# Corresponding penalty factor

# Estimates model 
 
glmmod.effect <-glmnet(x.effect,y,alpha=1,family='poisson')

# Determines best lambda with cross-validations

cv.glmmod.effect <- cv.glmnet(x.effect,y,alpha=1,family='poisson')

opt.lam.effect      <- c(cv.glmmod.effect$lambda.1se)

lasso.coefs.effect  <- coef(cv.glmmod.effect, s = opt.lam.effect)


new_x <- Frequencies %>% data_grid(Year,Educ_W,Educ_H,.model=model.hypuc) %>% 
	select(Year,Educ_W,Educ_H) %>% model.matrix(f, .)

predictions.effect <- cbind(Frequencies %>% data_grid(Year,Educ_W,Educ_H,.model=model.hypuc), predict(cv.glmmod.effect, new_x, s=opt.lam.effect)) %>%
	as_tibble() %>% rename(pred = `1`) %>% mutate(contrast = "effect")


predictions.effect <- predictions.effect %>% left_join(Frequencies %>% as_tibble() %>% select(F,Year,Educ_H,Educ_W), by=c("Year","Educ_H","Educ_W"))


predictions.effect %>% summarise(cor1 = cor(log(F),pred), cor2=cor(F,exp(pred)))

plot.res.effect <- predictions.effect %>% ggplot(aes(x=log(F),y= log(F) - pred)) + geom_point(alpha=0.3) + ylim(-4,1) + 
geom_hline(aes(yintercept=0), colour="red") + labs(subtitle="Effect coding", x="Log(counts)", y="residual")


# Plot for effect coding 

intercept <- predictions.effect %>% filter(Year==1990,Educ_H==1, Educ_W==1) %>% with(pred)

 	# full prediction
	predictions.effect <- cbind(Frequencies %>% data_grid(Year,Educ_W,Educ_H,.model=model.hypuc), predict(cv.glmmod.effect, new_x, s=opt.lam.effect)) %>%
		as_tibble() %>% rename(pred = `1`) %>% mutate(pred = pred - intercept)

	# margins
	predictions.effect_year  <- predictions.effect %>% filter(Educ_H==1, Educ_W==1)  %>% rename(margin_year=pred) %>% select(Year,margin_year)
	predictions.effect_educh <- predictions.effect %>% filter(Year==1990, Educ_W==1) %>% rename(margin_educh=pred) %>% select(Educ_H,margin_educh)
	predictions.effect_educw <- predictions.effect %>% filter(Year==1990, Educ_H==1) %>% rename(margin_educw=pred) %>% select(Educ_W,margin_educw)


	# match
	
	predictions.effect <- predictions.effect %>% left_join(predictions.effect_year, by="Year")
	predictions.effect <- predictions.effect %>% left_join(predictions.effect_educh, by="Educ_H")
	predictions.effect <- predictions.effect %>% left_join(predictions.effect_educw, by="Educ_W")
	

	# margins by year
	predictions.effect_year_educh <- predictions.effect %>% filter(Educ_H!=1,Year!=1990,Educ_W==1) %>%
		rename(margin_year_educh=pred) %>% mutate(margin_year_educh = margin_year_educh - (margin_year + margin_educh )) %>%
		select(Year,Educ_H,margin_year_educh)

	predictions.effect_year_educw <- predictions.effect %>% filter(Educ_H==1,Year!=1990,Educ_W!=1) %>%
		rename(margin_year_educw=pred) %>% mutate(margin_year_educw = margin_year_educw - (margin_year + margin_educw )) %>% 
		select(Year,Educ_W,margin_year_educw)

	predictions.effect <- predictions.effect %>% left_join(predictions.effect_year_educh, by=c("Year","Educ_H")) %>%  replace_na(list(margin_year_educh = 0))
	predictions.effect <- predictions.effect %>% left_join(predictions.effect_year_educw, by=c("Year","Educ_W")) %>%  replace_na(list(margin_year_educw = 0))


	# computes assortative mating
	predictions.effect <- predictions.effect %>% 
		mutate(assort = pred - (margin_year + margin_educh + margin_educw + margin_year_educh + margin_year_educw) ) %>% select(-diag,-move_asym)

	# plot
	predictions.effect   %>%
	ggplot(aes(x=factor(Year), y=assort, colour=factor(Educ_H, labels = educ.labels), group=factor(Educ_H) )) + geom_line(alpha=0.5)  + geom_point(size=1.5) + 
	facet_grid( . ~ Educ_W, labeller=as_labeller(educ.labeller) ) +
	theme(axis.text.x = element_text(angle=90, vjust=0.5)) + 
	labs(subtitle= "Wive's Educ",  y="Log-Odds", x="", colour= "Husband's Educ", title = "")  + ylim(-2.5, 7.5) +
	scale_color_d3() + scale_fill_d3()





####### Dummy coding ########

options(contrasts=c("contr.treatment","contr.poly"))

# Formula that includes all variables and interactions in a given dataframe

f <- as.formula(~ .*.*.)

y <- Frequencies$F

# Dataframe with Year, Educ_H and Educ_W

df <- Frequencies[,c(1,3,2)]

# Creates all three-way interactions

x.dummy <- model.matrix(f, df) # design matrix with interaction between Year, Educ_H and Educ_W


# Corresponding penalty factor

# Estimates model 
 
glmmod.dummy <- glmnet(x.dummy,y,alpha=1,family='poisson')

# Determines best lambda with cross-validations

cv.glmmod.dummy <- cv.glmnet(x.dummy,y,alpha=1,family='poisson')

opt.lam.dummy      <- c(cv.glmmod.dummy$lambda.1se)
lasso.coefs.dummy  <- coef(cv.glmmod.dummy, s = opt.lam.dummy)


new_x <- Frequencies %>% data_grid(Year,Educ_W,Educ_H,.model=model.hypuc) %>% 
	select(Year,Educ_W,Educ_H) %>% model.matrix(f, .)

predictions.dummy <- cbind(Frequencies %>% data_grid(Year,Educ_W,Educ_H,.model=model.hypuc), predict(cv.glmmod.dummy, new_x, s=opt.lam.dummy)) %>%
	as_tibble() %>% rename(pred = `1`)  %>% mutate(contrast = "dummy")

predictions.dummy <- predictions.dummy %>% left_join(Frequencies %>% as_tibble() %>% select(F,Year,Educ_H,Educ_W), by=c("Year","Educ_H","Educ_W"))


predictions.dummy  %>% summarise(cor1 = cor(log(F),pred), cor2=cor(F,exp(pred)))

plot.res.dummy  <- predictions.dummy %>% ggplot(aes(x=log(F),y= log(F) - pred)) + geom_point(alpha=0.3) + ylim(-4,1) + 
geom_hline(aes(yintercept=0), colour="red") + labs(subtitle="Dummy coding", x="Log(counts)", y="residual")




# Plot for dummy coding 

intercept <- predictions.dummy %>% filter(Year==1990,Educ_H==1, Educ_W==1) %>% with(pred)

 	# full prediction
	predictions.dummy <- cbind(Frequencies %>% data_grid(Year,Educ_W,Educ_H,.model=model.hypuc), predict(cv.glmmod.dummy, new_x, s=opt.lam.dummy)) %>%
		as_tibble() %>% rename(pred = `1`) %>% mutate(pred = pred - intercept)

	# margins
	predictions.dummy_year  <- predictions.dummy %>% filter(Educ_H==1, Educ_W==1)  %>% rename(margin_year=pred) %>% select(Year,margin_year)
	predictions.dummy_educh <- predictions.dummy %>% filter(Year==1990, Educ_W==1) %>% rename(margin_educh=pred) %>% select(Educ_H,margin_educh)
	predictions.dummy_educw <- predictions.dummy %>% filter(Year==1990, Educ_H==1) %>% rename(margin_educw=pred) %>% select(Educ_W,margin_educw)


	# match
	
	predictions.dummy <- predictions.dummy %>% left_join(predictions.dummy_year, by="Year")
	predictions.dummy <- predictions.dummy %>% left_join(predictions.dummy_educh, by="Educ_H")
	predictions.dummy <- predictions.dummy %>% left_join(predictions.dummy_educw, by="Educ_W")
	

	# margins by year
	predictions.dummy_year_educh <- predictions.dummy %>% filter(Educ_H!=1,Year!=1990,Educ_W==1) %>%
		rename(margin_year_educh=pred) %>% mutate(margin_year_educh = margin_year_educh - (margin_year + margin_educh )) %>%
		select(Year,Educ_H,margin_year_educh)

	predictions.dummy_year_educw <- predictions.dummy %>% filter(Educ_H==1,Year!=1990,Educ_W!=1) %>%
		rename(margin_year_educw=pred) %>% mutate(margin_year_educw = margin_year_educw - (margin_year + margin_educw )) %>% 
		select(Year,Educ_W,margin_year_educw)

	predictions.dummy <- predictions.dummy %>% left_join(predictions.dummy_year_educh, by=c("Year","Educ_H")) %>%  replace_na(list(margin_year_educh = 0))
	predictions.dummy <- predictions.dummy %>% left_join(predictions.dummy_year_educw, by=c("Year","Educ_W")) %>%  replace_na(list(margin_year_educw = 0))


	# computes assortative mating
	predictions.dummy <- predictions.dummy %>% 
		mutate(assort = pred - (margin_year + margin_educh + margin_educw + margin_year_educh + margin_year_educw) ) %>% select(-diag,-move_asym)

	# plot
	predictions.dummy   %>%
	ggplot(aes(x=factor(Year), y=assort, colour=factor(Educ_H, labels = educ.labels), group=factor(Educ_H) )) + geom_line(alpha=0.5)  + geom_point(size=1.5) + 
	facet_grid( . ~ Educ_W, labeller=as_labeller(educ.labeller) ) +
	theme(axis.text.x = element_text(angle=90, vjust=0.5)) + 
	labs(subtitle= "Wive's Educ",  y="Log-Odds", x="", colour= "Husband's Educ", title = "")  + ylim(-2.5, 7.5) +
	scale_color_d3() + scale_fill_d3()





