

#  ====================================== Selects variables for analysis ====================================== 


# ::::::::::::: Creates contingency table from microdata :::::::::::::::::::


Frequencies <- cont.table(mydata)
Frequencies <- llm.pars(Frequencies) 


#  ======================================== Estimates different models ======================================== 


# ::::::::::::::::::::::::::::::::: Null Model ::::::::::::::::::::::::::::::::


model.null <- glm(F ~ 1, family=poisson, data=Frequencies)
gof(Frequencies$F,model.null)

# ::::::::::::::::::::::::::::::: Saturante Model :::::::::::::::::::::::::::::


model.sat <- glm(F ~ Year*Educ_W*Educ_H, family=poisson, data=Frequencies)
gof(Frequencies$F,model.sat)

# :::::::::::::::::::::::::::: Model of Independence :::::::::::::::::::::::::::

model.indep <- glm(F ~ Year + Educ_W + Educ_H, family=poisson, data=Frequencies)
gof(Frequencies$F,model.indep)


# ::::::::::::::::::::: Model of *Conditional Independence :::::::::::::::::::::

model.condindep <- glm(F ~ Year + Educ_W + Educ_H +  Educ_W*Year + Educ_H*Year + Educ_W*Educ_H , family=poisson, data=Frequencies)
gof(Frequencies$F,model.condindep)


# ::::::::::::::: Model of Quasi Perfect Mobility Constrained ::::::::::::::::::

model.qpmc <- glm(F ~ Year + Educ_W + Educ_H +  Educ_W*Year + Educ_H*Year + Educ_W*Educ_H  + diagC + diagC*Year, family=poisson, data=Frequencies)
gof(Frequencies$F,model.qpmc)


# :::::::::::::: Model of Quasi Perfect Mobility Unconstrained :::::::::::::::::

model.qpm <- glm(F ~ Year + Educ_W + Educ_H +  Educ_W*Year + Educ_H*Year + Educ_W*Educ_H  + diag + diag*Year, family=poisson, data=Frequencies)
gof(Frequencies$F,model.qpm)


# ::: Model of  Symmetric movement across the minor diagonals + Main diagonal constrained :::

model.sdc <- glm(F ~ Year + Educ_W + Educ_H +  Educ_W*Year + Educ_H*Year + Educ_W*Educ_H  + diagC + diagC*Year + move_sym + move_sym*Year, family=poisson, data=Frequencies)
gof(Frequencies$F,model.sdc)


# ::: Model of Assymmetric movement across the minor diagonals + Main diagonal constrained :::

model.hyp <- glm(F ~ Year + Educ_W + Educ_H +  Educ_W*Year + Educ_H*Year + Educ_W*Educ_H  + diagC + diagC*Year + move_asym + move_asym*Year, family=poisson, data=Frequencies)
gof(Frequencies$F,model.hyp)

# ::: Model of  Symmetric movement across the minor diagonals + Main diagonal unconstrained :::

model.sduc <- glm(F ~ Year + Educ_W + Educ_H +  Educ_W*Year + Educ_H*Year + Educ_W*Educ_H  + diag + diag*Year + move_sym + move_sym*Year, family=poisson, data=Frequencies)
gof(Frequencies$F,model.sduc)


# ::: Model of Assymmetric movement across the minor diagonals + Main diagonal unconstrained :::

model.hypuc <- glm(F ~ Year + Educ_W + Educ_H +  Educ_W*Year + Educ_H*Year + Educ_W*Educ_H  + diag + diag*Year + move_asym + move_asym*Year, family=poisson, data=Frequencies)
gof(Frequencies$F,model.hypuc)

# ::::::::::::::::::::::::::::::::::::::: Crossing Models :::::::::::::::::::::::::::::::::::::::::::::

model.crossing <- glm(F ~ Year + Educ_W + Educ_H +  Educ_W*Year + Educ_H*Year + Educ_W*Educ_H  + 
crossing_1 + crossing_2 + crossing_3 + crossing_4 + crossing_1*Year + crossing_2*Year + crossing_3*Year + crossing_4*Year,
family=poisson, data=Frequencies)

gof(Frequencies$F,model.crossing)

# ::::::::::::::::::::::::::::: Crossing Models + Constrained diagonal ::::::::::::::::::::::::::::::::

model.crossingH <- glm(F ~ Year + Educ_W + Educ_H +  Educ_W*Year + Educ_H*Year + Educ_W*Educ_H  + diagC + diagC*Year +
crossing_1 + crossing_2 + crossing_3 + crossing_4 + crossing_1*Year + crossing_2*Year + crossing_3*Year + crossing_4*Year,
family=poisson, data=Frequencies)

gof(Frequencies$F,model.crossingH)

# ::::::::::::::::::::::::::::: Crossing Models + Unconstrained diagonal ::::::::::::::::::::::::::::::::

model.crossingMD <- glm(F ~ Year + Educ_W + Educ_H +  Educ_W*Year + Educ_H*Year + Educ_W*Educ_H  + diag + diag*Year +
crossing_1 + crossing_2 + crossing_3 + crossing_4 + crossing_1*Year + crossing_2*Year + crossing_3*Year + crossing_4*Year,
family=poisson, data=Frequencies)

gof(Frequencies$F,model.crossingMD)

# ::::::::::::::::::::::::::::::::::::: ORDINAL MODELS ::::::::::::::::::::::::::::::::

# ::::::::::::::::::::::::::::::::::::: Linear by Linear Association ::::::::::::::::::::::::::::::::

model.linear <- glm(F ~ Year + Educ_W + Educ_H + Educ_W*Year + Educ_H*Year + Educ_W*Educ_H + as.numeric(LA)*Year , family=poisson, data=Frequencies)
gof(Frequencies$F,model.linear)


# ::::::::::::::::::::::::::::::::::::: UNIDIFF model ::::::::::::::::::::::::::::::::

model.unidiff <- gnm(F ~ Year*Educ_W + Year*Educ_H + Mult(Year, Educ_W:Educ_H),
                     family = poisson,
                     data = Frequencies)
gof.unidiff(Frequencies$F,model.unidiff)



# ================================== Table with statistcis of Goodness of fit ==================================


models <- c("model.indep","model.condindep","model.qpmc","model.qpm","model.sdc","model.hyp","model.sduc","model.hypuc",
	        "model.crossing","model.crossingH","model.crossingMD","model.linear","model.unidiff")


table.gof   <- NULL # Empty table

for (i in models) {

	# Computes Goodness of fit statistics for each model
	
	if (i=="model.unidiff") {
		tgof <- gof.unidiff(Frequencies$F, eval(parse(text = i))  )

	}
	else {
		tgof <- gof(Frequencies$F, eval(parse(text = i))  )
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


# Exports to latex
		 
		title.table <- "Log-Linear Models of the Association Between Husband’s and Wife’s Educational Attainment (Husbands Aged 30-35): Chile, 1990-2015"
		lab <- paste0("tab:","Gofllm_chile")

		print(xtable(table.gof, caption = title.table, align = "lccccc", label = lab), type="latex", floating=TRUE, table.placement="h!", hline.after=c(-1,0,13),
		size = "footnotesize", latex.environments = "center", caption.placement = "top",
		sanitize.colnames.function = identity, file= paste0(path.output.tables,"Gofllm_chile",".tex"))




# ::::::::::::::::::::::: Plots predicted log Odds of Assortative Mating :::::::::::::::::::::::


educ.labels = c("E-","E","H","C-","C")
educ.labeller = c(`1`="E-",`2`="E",`3`="H",`4`="C-",`5`="C")


#### Predictions with LLM parameters set to zero

intercept = model.hypuc$coefficients[1]

predictions <- Frequencies %>% data_grid(Year,Educ_W,Educ_H,.model=model.hypuc) %>% mutate_all(funs(factor(.))) %>%
add_predictions(model.hypuc) %>% mutate(pred = pred - intercept)  %>% select(Year,Educ_H,Educ_W,pred)

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


# assortative mating baseline 
 
predictions_assort_base <- predictions %>% filter(Year==1990) %>% rename(assort_base= pred) %>% mutate(assort_base = assort_base - margin_educh - margin_educw) %>% select(Educ_H,Educ_W,assort_base)

predictions <- predictions %>% left_join(predictions_assort_base, by=c("Educ_H","Educ_W")) %>%  replace_na(list(margin_year_educw = 0)) 


# Prediction adding LLM parameters
predictions_full <- Frequencies %>% data_grid(Year,Educ_W,Educ_H,.model=model.hypuc) %>%
	mutate(diag=seq(1,5,1)*(Educ_W==Educ_H),  move_asym=((as.numeric(Educ_H)-as.numeric(Educ_W))==-1)*1 + ((as.numeric(Educ_H)-as.numeric(Educ_W)==1)*2)) %>%
	mutate_all(funs(factor(.))) %>%
	add_predictions(model.hypuc) %>% rename(pred_full = pred) %>% mutate(pred_full = pred_full - intercept) 

predictions <- predictions %>% left_join(predictions_full, by=c("Year","Educ_H","Educ_W")) %>%  replace_na(list(margin_year_educw = 0)) 

# create assortative mating parameter

predictions <- predictions %>% mutate(assort = pred_full - (margin_year +  margin_educh + margin_educw + margin_year_educh + margin_year_educw ) )


plot.assort.ll <- 
	predictions %>% 
	ggplot(aes(x=Year %>% as.character(), y=assort, group=factor(Educ_H))) + 
	facet_grid( . ~ Educ_W, labeller=as_labeller(educ.labeller) ) +
	geom_line(aes(colour=factor(Educ_H, labels = educ.labels)), alpha=0.5) + 
	geom_point( aes(colour=factor(Educ_H, labels = educ.labels), shape=factor(Educ_H, labels = educ.labels)), size=1.5) + 
	theme(axis.text.x = element_text(angle=90, vjust=0.5, size=8)) + 
	labs(subtitle= "Wive's Educ",  y="Log-Odds", x="", shape="Husband's Educ", colour= "Husband's Educ", title = "")  +
	scale_color_d3() + scale_fill_d3()



key<-paste0(path.output.images, "logOdds_hypuc",".pdf")

save_plot(key, plot.assort.ll, base_height = 4 , base_width = 12)




