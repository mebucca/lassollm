
# =======================================================================================================================================================
# ==================================================================== PREAMBLE =========================================================================
# =======================================================================================================================================================


# Clear Screen 

cat("\014")
rm(list = ls())


# Load Packages 

library("gridGraphics")
library("pryr")
library("foreign")
library("lme4")
library("pander")
library("stargazer")
library("xtable")
library("ggplot2")
library("splitstackshape")
library("plyr")
library("lars")
library("coda")
library("hglm")
library("glmnet")
library("ade4")
library("wesanderson")
library("grplasso")
library("caret")
library("grid")
library("scales")
library("scales")
library("tidyverse")
library("broom")
library("modelr")
library("ggsci")
library("ggthemes")
library("gnm")
library("cowplot")
library("glmc")


options(xtable.floating = FALSE)
options(xtable.timestamp = "")


# Set preferences ggplot theme



theme_set(theme_few() +  
	theme(strip.text.x = element_text(colour = "black", face = "bold"),
	text = element_text(colour = "black", size = 12),	
	strip.text.y = element_text(colour = "black", face = "bold"),
	#panel.grid.minor = element_blank(),
	panel.grid.major = element_blank(),
	panel.background = element_rect(colour = "black", fill = "white", size=1.05),
	axis.text.x = element_text(angle=90, vjust=0.5), 
	strip.background = element_blank()
	)
)


# Open Dataset 

users <- c("Mauricio","danielaurbina","mbucca")
user  <- users[1] 

path.output          <- paste("/Users/",user,"/Dropbox/Apps/ShareLaTeX/SMR_RR/",sep="")
path.output.images   <- paste("/Users/",user,"/Dropbox/Apps/ShareLaTeX/SMR_RR/images/",sep="")
path.output.tables   <- paste("/Users/",user,"/Dropbox/Apps/ShareLaTeX/SMR_RR/tables/",sep="")
path.analysis        <- paste("/Users/",user,"/Google Drive/Research/Assortative Mating and Inequality Chile/Analysis/",sep="")
path.data            <- paste("/Users/",user,"/Google Drive/Research/Assortative Mating and Inequality Chile/Data/", sep="")
setwd(path.data)


# Inputs microdata from CASEN
name <- paste("casen_un","csv", sep=".")
mydata <- read.csv(name)

# Set number of iterations for Lasso repeated CV lambda selection
n_sims <- 10

# =======================================================================================================================================================
# =================================================================== FUNCTIONS =========================================================================
# =======================================================================================================================================================

setwd(path.analysis)
source("Functions.R")

# =======================================================================================================================================================
# =============================================================== LOG-LINEAR MODELS =====================================================================
# =======================================================================================================================================================


setwd(path.analysis)
source("4_LogLinearModeling.R")

# =======================================================================================================================================================
# ==================================================================== LASSO ============================================================================
# =======================================================================================================================================================


setwd(path.analysis)
source("5_Lasso.R")

# Combine plots

key  <- paste0(path.output.images, "logOdds_Lasso.pdf")

plot.assort.lasso <- 
	plot_grid(plot_logOdds_Lassofree + theme(legend.position="none"), get_legend(plot_logOdds_Lassofree),
	plot_logOdds_Lassoadaptive + theme(legend.position="none"), NULL,
	plot_logOdds_Lassoindep + theme(legend.position="none"), NULL,
	plot_logOdds_Lassohywy + theme(legend.position="none"), NULL,
	labels = c("A", NA,"B", NA, "C", NA, "D", NA), nrow=4, ncol=2, rel_widths = c(4, 1.5))


save_plot(key, plot.assort.lasso, base_height =12, base_width=10)


# PLot Coefficients path

key  <- paste0(path.output.images, "CoefPath_Lasso.pdf")

plot.coeffpath.lasso <- 
	plot_grid(plot_Coefpath_Lassofree, 
	plot_Coefpath_Lassoadaptive, 
	plot_Coefpath_Lassoindep, 
	plot_Coefpath_Lassohywy, 
	labels = c("A","B", "C", "D"), nrow=2, ncol=2, scale = 1)

save_plot(key, plot.coeffpath.lasso, base_height =11, base_width=11)


# PLot CV error

key  <- paste0(path.output.images, "CVerror_Lasso.pdf")

plot.cverror.lasso <- 
	plot_grid(plot_CVerror_Lassofree, 
	plot_CVerror_Lassoadaptive,
	plot_CVerror_Lassoindep,
	plot_CVerror_Lassohywy,
	labels = c("A","B", "C", "D"), nrow=2, ncol=2, scale = 1)


save_plot(key, plot.cverror.lasso, base_height =11, base_width=11)



# =======================================================================================================================================================
# ============================================================= MONTE CARLO SIMULATION ==================================================================
# =======================================================================================================================================================

set.seed(8610)

setwd(path.analysis)
source("6_MCS.R")


# Distribution of Frequencies

key  <- paste0(path.output.images, "Distrib_Freqs_MC.pdf")

plot.distrib.freqs <- 
	plot_grid(Distrib_Freqs_19, Distrib_Freqs_117, Distrib_Freqs_828, labels = c("A", "B","C", NA), nrow=2, ncol=2)

save_plot(key, plot.distrib.freqs , base_height =11, base_width=11)

# PLot Coefficients path

key  <- paste0(path.output.images, "CoefPath_Lasso_MC.pdf")

plot.mc.coef.path <- 
	plot_grid(MC_coef_path_19, MC_coef_path_117, MC_coef_path_828, labels = c("A", "B","C", NA), nrow=2, ncol=2, scale = 1)

save_plot(key, plot.mc.coef.path , base_height =11, base_width=11)


# PLot CV error

key  <- paste0(path.output.images, "CVerror_Lasso_MC.pdf")

plot.mc.cv.error <- 
	plot_grid(MC_CV_error_19, MC_CV_error_117, MC_CV_error_828, labels = c("A", "B","C", NA), nrow=2, ncol=2, scale = 1)

save_plot(key, plot.mc.cv.error, base_height =11, base_width=11)


# PLot Assortative Mating Coefficients
 
key  <- paste0(path.output.images, "logOdds_Lasso_MC.pdf")

plot.mc.assort.lasso <- 
	plot_grid(MC_Lasso_free_logOdds_19 + theme(legend.position="none"), get_legend(MC_Lasso_free_logOdds_19),
	MC_Lasso_free_logOdds_117 + theme(legend.position="none"), NULL,
	MC_Lasso_free_logOdds_828 + theme(legend.position="none"), NULL,
	labels = c("A", NA,"B", NA, "C", NA), nrow=3, ncol=2, rel_widths = c(4, 1.5))

save_plot(key, plot.mc.assort.lasso, base_height =9, base_width=10)



# =======================================================================================================================================================
# ============================================================== CROSS VALIDATION ==================================================================
# =======================================================================================================================================================

set.seed(6621)

setwd(path.analysis)
source("7_Cross_Validation.R")


print(paste0("================================================= DONE! ================================================="))





