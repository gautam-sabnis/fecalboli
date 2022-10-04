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