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
    fit <- suppressWarnings(nlme(value ~ SSgompertz(variable, Asym, b2, b3), data = df_Strain, fixed = Asym+b2+b3 ~ Sex, random = pdDiag(Asym+b2+b3 ~ 1), groups = ~MouseID, control = nlmeControl(maxIter = 5000, msMaxIter = 1500)))
    fit_fe <- data.frame(Asym = rep(fixef(fit)[1], length(unique(df_Strain$MouseID))), b2 = rep(fixef(fit)[2], length(unique(df_Strain$MouseID))), b3 = rep(fixef(fit)[3], length(unique(df_Strain$MouseID))))
    fit_beta <- fit_fe + ranef(fit)
    df_params <- data.frame(Asym = fit_beta$Asym, k = -log(fit_beta$b3), t_m = -log(fit_beta$b2)/log(fit_beta$b3))
    df_params$MouseID <- unique(df_Strain$MouseID)
    df_params <- df_params[,c("MouseID", "Asym", "k", "t_m")]

    df_fit <- augPred(fit,level=0:1)
	for (x in seq(unique(df_fit$.groups))){
  		df_fit[df_fit$.group == unique(df_fit$.group)[x], "MouseID"] <- (rep(x,table(df_fit$.groups)[x]))
	}
	df_fit$MouseID <- df_fit$.groups
	p <- ggplot(df_fit, aes(x = variable, y=value)) + geom_point(data = df_fit[df_fit$.type == "original",], aes(x = variable, y=value), alpha = 0.5) + geom_line(data = df_fit[df_fit$.type == "predict.MouseID",], aes(x = variable, y = value), color = 'magenta') + geom_line(data = df_fit[df_fit$.type == "predict.fixed",], aes(x=variable, y = value), color = 'blue') + facet_wrap(.~MouseID, scales="free") + labs(x = 'Time (in mins)', y = 'Value') + ggtitle(paste0(Strain)) + theme_bw()

	return(list(p, df_params))
	

}