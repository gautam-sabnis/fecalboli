#' Load and preprocess data 
#' 
#' Loads the data and preprocess the features. 
#' @param csv a csv file containing data per stride
#' @return a data frame containing the relevant fecal boli data 
#' @examples 
#' load_data(data = "fecalboli.csv")
#' @export 

load_data <- function(csv){

	data <- read.csv(paste0(csv), stringsAsFactors = TRUE)
	data <- cbind(data, TestAge=as.Date(FBdata$TestDate, format='%m/%d/%Y') - as.Date(FBdata$DOB,format='%m/%d/%Y')) #days
	data <- data[!duplicated(data$MouseID),]

}