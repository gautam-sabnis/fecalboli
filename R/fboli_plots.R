#' Exploratory plot 
#' 
#' Plots the raw data (both equally and unequally spaced time points) 
#' @param data data frame containing equally/unequally spaced time 
#' points
#' @return plot
#' @examples 
#' plot_data(data = df_equal)
#' @export 

plot_data <- function(data){

	p <- ggplot(df_unequal, aes(x = variable, y = value, fill = MouseID)) + geom_point() + geom_line() + theme(legend.position = "none") + facet_wrap(. ~ Strain, scales = "free") + labs(x = "Time", y = "Count")
	return(p)

}

#' Exploratory plot 
#' 
#' Plots the raw data (both equally and unequally spaced time points)
#' for a single animal/strain
#' @param data data frame containing equally/unequally spaced time 
#' points
#' @param MouseID mouseID
#' @param Strain strain
#' @return plot
#' @examples 
#' plot_strain_animal(data = df_equal, MouseID = "LL1-1_BalbcJ",
#' Strain = NULL)
#' plot_strain_animal(data = df_equal, MouseID = NULL, Strain = 
#' Strain = "BALB/cJ")
#' @export 

plot_strain_animal <- function(data, MouseID, Strain){

	if (is.null(Strain)){
		df_tmp <- data[data$MouseID %in% MouseID, ]
		p <- ggplot(df_tmp, aes(x = variable, y = value)) + geom_point() + geom_line() + theme_bw(base_size = 16) + ggtitle(paste0("Strain:", unique(df_tmp$Strain), ", MouseID:", paste0(MouseID))) + labs(x = "Time", y = "Count")
	} else {
		df_tmp <- data[data$Strain %in% Strain, ]
		p <- ggplot(df_tmp, aes(x = variable, y = value)) + geom_point() + geom_line() + theme_bw(base_size = 16) + facet_wrap(. ~ MouseID, scales = "free") + ggtitle(paste0("Strain:", unique(df_tmp$Strain))) + labs(x = "Time", y = "Count")
	}

	return(p)

}


#' Plot fitted growth curves on raw data 
#' 
#' Plots the fitted growth curve for all animals belonging to a strain
#' @param data data frame containing equally/unequally spaced time 
#' points
#' @param Strain strain
#' @return a list containing the plot and a data frame with estimated #' growth parameters 
#' @examples 
#' plot_growth_model(data = df_unequal, Strain = "C57BL/6NJ")
#' @export 

plot_growth_model <- function(data, Strain, shrinkage){

	df_Strain <- data[data$Strain %in% Strain, ]
	initVals <- getInitial(value ~ SSgompertz(variable, Asym, b2, b3), data = df_Strain)
    fit <- suppressWarnings(nlme(value ~ SSgompertz(variable, Asym, b2, b3), data = df_Strain, fixed = Asym+b2+b3 ~ 1, random = pdDiag(Asym+b2+b3 ~ 1), groups = ~MouseID, control = nlmeControl(maxIter = 5000, msMaxIter = 1500)))
    fit_fe <- data.frame(Asym = rep(fixef(fit)[1], length(unique(df_Strain$MouseID))), b2 = rep(fixef(fit)[2], length(unique(df_Strain$MouseID))), b3 = rep(fixef(fit)[3], length(unique(df_Strain$MouseID))))
    fit_beta <- fit_fe + ranef(fit)
    df_params <- data.frame(Asym = fit_beta$Asym, k = -log(fit_beta$b3), t_m = -log(fit_beta$b2)/log(fit_beta$b3))
    df_params$MouseID <- unique(df_Strain$MouseID)
    df_params <- df_params[,c("MouseID", "Asym", "k", "t_m")]
    tmp <- ranef(fit, augFrame = T)
                
    Asym_covariate <- plot(tmp, form = Asym ~ Sex + CoatColor)
    b2_covariate <- plot(tmp, form = b2 ~ Sex + CoatColor)
    b3_covariate <- plot(tmp, form = b3 ~ Sex + CoatColor)

    df_fit <- augPred(fit,level=0:1)
	for (x in seq(unique(df_fit$.groups))){
  		df_fit[df_fit$.group == unique(df_fit$.group)[x], "MouseID"] <- (rep(x,table(df_fit$.groups)[x]))
	}
	df_fit$MouseID <- df_fit$.groups
	p <- ggplot(df_fit, aes(x = variable, y=value)) + geom_point(data = df_fit[df_fit$.type == "original",], aes(x = variable, y=value), alpha = 0.5) + geom_line(data = df_fit[df_fit$.type == "predict.MouseID",], aes(x = variable, y = value), color = 'magenta') + geom_line(data = df_fit[df_fit$.type == "predict.fixed",], aes(x=variable, y = value), color = 'blue') + facet_wrap(.~MouseID, scales="free") + labs(x = 'Time (in mins)', y = 'Value') + ggtitle(paste0(Strain)) + theme_bw()

	return(list(p, df_params, Asym_covariate, b2_covariate, b3_covariate))
	

}

#' Plot simulated growth curves 
#' 
#' Plots the simulated growth curves using estimated growth 
#' parameters for comparison across all strains
#' @param data an aggregated data frame containing estimated growth #' parameters for each strain 
#' @param type type of plot
#' @return a plot containing simulated growth curves for all strains 
#' @examples 
#' plot_growth_comparison(data = df_growth)
#' @export 

plot_growth_comparison <- function(data, type){

	if (type == "Strain"){

		df_growth <- aggregate(data[, which(names(data) %in% c("Asym","k","t_m"))], list(data$Strain), mean)
		names(df_growth)[1] <- "Strain"
		df_growth <- df_growth[!apply(df_growth[,-1], 1, function(r) all(r == 0.00)), ]
		Strains <- df_growth$Strain
		df_sim <- list()
		x <- seq(0,60,length.out = 60)
		invisible(lapply(seq_along(Strains), function(s) {
  			Asym <- df_growth[s,"Asym"]
  			slope <- df_growth[s,"k"]
  			inflection <- df_growth[s,"t_m"]  
  			df_sim[[s]] <<- data.frame(Strain = paste0(Strains[s]), t = x, f = Asym*exp(-exp(-slope*(x - inflection))))
  		}))
		df_sim <- do.call(rbind, df_sim)

		p <- ggplot(df_sim, aes(x = t, y = f, color = Strain)) + geom_line() + theme_bw(base_size = 16) + labs(x = "Time", y = "FBoli Count") + theme(legend.position = "none")
		return(p)
	} else if (type == "Sex"){

		df_growth <- aggregate(data[, which(names(data) %in% c("Asym","k","t_m"))], list(data$Strain, data$Sex), mean)
		names(df_growth)[1:2] <- c("Strain", "Sex")
		df_growth <- df_growth[!apply(df_growth[,-1], 1, function(r) all(r == 0.00)), ]
		df_sim <- list()
		x <- seq(0,60,length.out = 60)
		invisible(lapply(seq(nrow(df_growth)), function(s) {
  			Asym <- df_growth[s,"Asym"]
  			slope <- df_growth[s,"k"]
  			inflection <- df_growth[s,"t_m"]  
  			df_sim[[s]] <<- data.frame(Strain = df_growth[s, "Strain"], Sex = df_growth[s,"Sex"], t = x, f = Asym*exp(-exp(-slope*(x - inflection))))
  		}))
		df_sim <- do.call(rbind, df_sim)
		p <- ggplot(df_sim, aes(x = t, y = f, color = Sex)) + geom_line(lwd = 1.1) + theme_bw(base_size = 16) + labs(x = "Time", y = "FBoli Count") + scale_color_brewer(palette = "Set1") + facet_wrap(. ~ Strain, scales = "free")
		return(p)	
	} else if (type == "CoatColor"){

		df_growth <- aggregate(data[, which(names(data) %in% c("Asym","k","t_m"))], list(data$CoatColor), mean)
		names(df_growth)[1] <- c("CoatColor")
		df_growth <- df_growth[!apply(df_growth[,-1], 1, function(r) all(r == 0.00)), ]
		df_sim <- list()
		x <- seq(0,60,length.out = 60)
		invisible(lapply(seq(nrow(df_growth)), function(s) {
  			Asym <- df_growth[s,"Asym"]
  			slope <- df_growth[s,"k"]
  			inflection <- df_growth[s,"t_m"]  
  			df_sim[[s]] <<- data.frame(CoatColor = df_growth[s,"CoatColor"], t = x, f = Asym*exp(-exp(-slope*(x - inflection))))
  		}))
		df_sim <- do.call(rbind, df_sim)
		p <- ggplot(df_sim, aes(x = t, y = f, color = CoatColor)) + geom_line(lwd = 1.1) + theme_bw(base_size = 16) + labs(x = "Time", y = "FBoli Count")
		return(p)	

	}


}


#' Plot grooming paper like plots 
#' 
#' Compare growth parameters with open field traits across strains
#' and highlight strains for which correlations are signficant 
#' @param data data frame containing estimated growth parameters 
#' and other open field and grooming traits  
#' @param of_trait open field trait
#' @param growth_param estimated growth curve parameter 
#' @return a plot highlighting Strains with significant correlations 
#' @examples 
#' plot_signif_corr(data, of_trait = "bin_avg_55.periphery_time_secs"
#' ,growth_param = "Asym")
#' @export 

plot_signif_corr <- function(data, of_trait, growth_param){

	r <- by(data, data$Strain, FUN = function(x) cor(x[,paste0(growth_param)], x[,paste0(of_trait)], method = "pearson"));
	rp <- by(data, data$Strain, FUN = function(x) ifelse(cor.test(x[,paste0(growth_param)], x[,paste0(of_trait)], method = "pearson", exact = FALSE)$p.value < 0.05, 1, 0));
	df_tmp <- data.frame(Strain = dimnames(r)[[1]], Corr = as.vector(r));
	df_tmp2 <- data.frame(Strain = dimnames(rp)[[1]], Corr = as.vector(rp));
	df_tmp[,2] <- df_tmp[,2]*df_tmp2[,2];
	df_tmp <- df_tmp[complete.cases(df_tmp),]

	df1 <- data[, names(data) %in% c("Strain", paste0(growth_param), paste0(of_trait))]
	df1[,paste0(of_trait)] <- df1[,paste0(of_trait)]/33
	df1_agg <- aggregate(. ~ Strain, data =  df1, FUN = function(x) c(m = mean(x), s = sd(x)))
	df_agg <- data.frame(Strain = df1_agg$Strain, Asym_m = df1_agg[,paste0(growth_param)][,1], Asym_s = df1_agg[,paste0(growth_param)][,2], corner_m = df1_agg[,paste0(of_trait)][,1], corner_s = df1_agg[,paste0(of_trait)][,2])

	df_agg$signif <- as.factor(ifelse(df_agg$Strain %in% df_tmp[!df_tmp$Corr == 0, "Strain"], 1, 0))

	p <- ggplot(df_agg, aes(x = Asym_m, y = corner_m)) + geom_point(size = 2, stroke = 1, aes(color = signif), alpha = 0.5) + geom_errorbar(aes(ymin = corner_m - corner_s, ymax = corner_m + corner_s, color = signif), alpha = 0.3) + geom_errorbarh(aes(xmin = Asym_m - Asym_s,xmax = Asym_m + Asym_s, color = signif), alpha = 0.3) + ggrepel::geom_text_repel(data = df_agg[df_agg$signif == 1,], aes(label = Strain, color = signif, size = 5)) + scale_colour_manual(values = c("black", "red")) + theme_bw(base_size = 16) + theme(legend.position = "none") + labs(x = paste0(growth_param), y = "Periphery Time %")

	return(p)

}