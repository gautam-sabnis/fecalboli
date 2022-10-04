#' Beta growth model 
#' 
#' Fit a beta growth model (with four parameters) to the data 
#' @param data data frame containing equally/unequally spaced time 
#' points
#' @param MouseID mouseID
#' @param Strain strain
#' @return fitted model?
#' @examples 
#' beta(data = df_equal, MouseID = "LL1-1_BalbcJ", Strain = "BALB/cJ")
#' @export 

beta <- function(data, Strain){

	df_tmp <- data[data$Strain %in% Strain, ]
	#df_tmp$MouseID <- droplevels(df_tmp$MouseID)

	fit.lis <- nlsList(value ~ SSbgf4(variable, w.max, t.e, t.m, t.b), data = df_tmp)
	fit.me <- nlme(fit.lis, control = list(maxIter = 100, msMaxIter = 300, pnlsMaxIter = 20))
	return(fit.me)
}