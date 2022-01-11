
# ::::::::::::: Creates a random frequencie with mean 1 (from poisson(1)) and higher that zero

small.freq <- function(Frequencies, mean) {
	
	n = length(Frequencies[which(Frequencies==0)])

	if (n>0) {

	freqs <- NULL # stores freqs
	
	for (i in 1:n) {
		f=0
		while (f==0) {f=rpois(1,mean)}
		freqs <- c(freqs,f)
		}
	}

	return(freqs=freqs)
} 


# ::::::::::::: Creates contingency table from microdata :::::::::::::::::::


cont.table <- function(data) {

	Frequencies <- table(factor(data$Educ_H),factor(data$Educ_W),factor(data$Year))

	# makes table into dataframe
	Frequencies <- data.frame(Frequencies) 
	names(Frequencies) <- c("Educ_H","Educ_W","Year","F") # name variables

	# re-arrange variables
	Frequencies <- data.frame(Frequencies$Year,Frequencies$Educ_H,Frequencies$Educ_W,Frequencies$F) 
	names(Frequencies) <- c("Year","Educ_H","Educ_W","F")

	# Eliminates rows with zeros
	#Frequencies <- Frequencies[Frequencies$F!=0,]   
	
	# An alternative: add small constant (1) to sampling-zeros
	Frequencies$F[Frequencies$F==0] <- small.freq(Frequencies$F,1)

	return(Frequencies = Frequencies)

}



# :::::::::::::::::::::::::: Parameters Loglinear Models ::::::::::::::::::::::::::::::

# Functions that creates special parameters for loglinear models


llm.pars <- function(Frequencies) {

	# constrained homogamy
	Frequencies$diagC <- (Frequencies$Educ_H==Frequencies$Educ_W)*1 


	# unconstrained homogamy
	Frequencies$diag  <- (Frequencies$Educ_H==Frequencies$Educ_W)*1

	for (i in 1:length(levels(Frequencies$Educ_H))) {
		Frequencies$diag[Frequencies$diag==1 & Frequencies$Educ_H==i] <- i
	}

	# Symmetric movements minor diagonals 

	Frequencies$move_sym <- (abs(as.numeric(Frequencies$Educ_H)-as.numeric(Frequencies$Educ_W))==1)*1

	# Asymmetric movements minor diagonals 


	Frequencies$move_asym <- ((as.numeric(Frequencies$Educ_H)-as.numeric(Frequencies$Educ_W))==-1)*1 + ((as.numeric(Frequencies$Educ_H)-as.numeric(Frequencies$Educ_W)==1)*2)
	
	# Linear-by-Linear Association 
	
	Frequencies$HL[Frequencies$Educ_H == 1] <- 8
	Frequencies$HL[Frequencies$Educ_H == 2] <- 10.5
	Frequencies$HL[Frequencies$Educ_H == 3] <- 12
	Frequencies$HL[Frequencies$Educ_H == 4] <- 14
	Frequencies$HL[Frequencies$Educ_H == 5] <- 17
	
	Frequencies$WL[Frequencies$Educ_W == 1] <- 8
	Frequencies$WL[Frequencies$Educ_W == 2] <- 10.5
	Frequencies$WL[Frequencies$Educ_W == 3] <- 12
	Frequencies$WL[Frequencies$Educ_W == 4] <- 14
	Frequencies$WL[Frequencies$Educ_W == 5] <- 17
	
	Frequencies$LA <- (as.numeric(Frequencies$HL)*as.numeric(Frequencies$WL))

	# Crossing parameters


	for (i in 1:(length(levels(Frequencies$Educ_H))-1)) {
		name <-paste("crossing_",i,sep="")

		Frequencies[[name]] <- 0
		Frequencies[[name]] <- (as.numeric(Frequencies$Educ_H)<=i &  as.numeric(Frequencies$Educ_W)>i)*1
		Frequencies[[name]] <- (as.numeric(Frequencies$Educ_W)<=i &  as.numeric(Frequencies$Educ_H)>i)*1

	}


	# Transform new vars to factors

	Frequencies$diagC      <- factor(Frequencies$diagC )
	Frequencies$diag       <- factor(Frequencies$diag)
	Frequencies$move_sym   <- factor(Frequencies$move_sym)
	Frequencies$move_asym  <- factor(Frequencies$move_asym)
	Frequencies$LA         <- factor(Frequencies$LA)
	Frequencies$crossing_1 <- factor(Frequencies$crossing_1)
	Frequencies$crossing_2 <- factor(Frequencies$crossing_2)
	Frequencies$crossing_3 <- factor(Frequencies$crossing_3)
	Frequencies$crossing_4 <- factor(Frequencies$crossing_4)


	return(Frequencies = Frequencies)
}


# ::::::::::::::::::::::::::: Goodness of Fit Stats ::::::::::::::::::::::::::::::

# Creates function that reports goodness of fit statistics after log-linear models 

gof <- function(freq,model) {
	n    <- sum(freq); 
	df   <- model$df.residual
	q    <- model$df.null - model$df.residual;
	G2   <- model$deviance
	D 	 <- 100*(sum(abs(freq-exp(predict(model))))/ (2*n))
	AIC  <- G2 - 2*df
	BIC  <- G2 + log(n)*q
	stats <- c(df,G2,D,AIC,BIC)
	names(stats) <- c("df","G2","100*D","AIC","BIC")
	return(stats)
}



gof.unidiff <- function(freq,model) {
	n    <- sum(freq); 
	df   <- model$df.residual
	q    <- length(model$coefficients)
	G2   <- model$deviance
	D 	 <- 100*(sum(abs(freq-exp(predict(model))))/ (2*n))
	AIC  <- G2 - 2*df
	BIC  <- G2 + log(n)*q
	stats <- c(df,G2,D,AIC,BIC)
	names(stats) <- c("df","G2","100*D","AIC","BIC")
	return(stats)
}

# ::::::::::::::::::::::::::::::: Poisson Deviance ::::::::::::::::::::::::::::::


poisson.dev <- function(F.observed,F.predicted) {
    pd = (F.observed* (log(F.observed/round(F.predicted,0)))) - (F.observed-round(F.predicted,0))
	pd[pd==Inf] <- NA
	pd[pd==-Inf] <- NA
    PD = 2*(mean(pd, na.rm=TRUE)*length(F.observed))  # equivalent to 2*sum(pd) but avoinding NAs
    return(PD = PD)
}


# ::::::::::::::::::::::::::::::: MSE ::::::::::::::::::::::::::::::

Mean.SE <- function(F.predicted,F.observed) {
	se =  (F.predicted-F.observed)^2
	se[se==Inf] <- NA
	se[se==-Inf] <- NA
	MSE  = mean(se, na.rm=TRUE)
	return(MSE = MSE)
}

# ::::::::::::::::::::::::::::::: Deviance Score ::::::::::::::::::::::::::::::

deviance.score <- function(F.predicted,F.observed) {
	ds = -2*log(dpois(round(F.predicted,0),F.observed))
	ds[ds==Inf] <- NA
	ds[ds==-Inf] <- NA
	DS  = mean(ds, na.rm=TRUE)
	return(DS = DS)
}



