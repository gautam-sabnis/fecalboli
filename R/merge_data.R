#' Load and preprocess data 
#' 
#' Loads the data and preprocess the features. 
#' @param df1 data frame containing growth parameters 
#' @param df2 data frame containing grooming and open field data
#' @return a data frame containing growth parameters and other traits #' merged on common MouseIDs 
#' @examples 
#' merge_data(df1, df2)
#' @export 

merge_data <- function(df1, df2){

	mIDs <- (intersect(unique(df1$MouseID),unique(df2$MouseID))) 
	df2 <- df2[df2$MouseID %in% mIDs, ]
	df1 <- df1[!df1$MouseID %in% setdiff(unique(df1$MouseID), unique(df2$MouseID)),]

	df1$MouseID <- droplevels(df1$MouseID)
	df2$MouseID <- droplevels(df2$MouseID)
	df1 <- df1[!duplicated(df1$MouseID),]
	df2 <- df2[!duplicated(df2$MouseID),]

	df <- merge(df1, df2)
	return(df)
}