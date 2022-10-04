#' Filter data 
#' 
#' Remove mice with too few data points  
#' @param data data frame containing unequally spaced time 
#' points
#' @param threshold minimum cutoff for unequally spaced time points
#' @return filtered data set 
#' @examples 
#' filter_data(data = df_unequal, threshold = 5)
#' @export 

filter_data <- function(data, threshold){

	mIDs <- names(which(table(data$MouseID) < threshold))
	data <- data[!data$MouseID %in% mIDs, ]
	data$MouseID <- droplevels(data$MouseID)
	return(data)

}