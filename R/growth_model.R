#' growth model 
#' 
#' Fit a growth model to the data 
#' @param data data frame containing equally/unequally spaced time 
#' points
#' @param type type of growth curve
#' @return Estimated growth curve parameters 
#' @examples 
#' growth_model(data = df_equal, type = "Gompertz")
#' @export 

growth_model <- function(data, type){

	Strains <- unique(data$Strain)
	#Strains <- setdiff(Strains, c("A/J","C3H/HeJ","C57L/J","CAST/EiJ"))
	param_list <- list()
	for (s in 1:length(Strains)){
		cat("Strain: ", paste0(Strains[s]), "\n")
		df_tmp <- data[data$Strain %in% Strains[s], ]
		#df_tmp <- filter_data(data = df_tmp, threshold = 5)
		df_tmp$MouseID <- droplevels(df_tmp$MouseID)
		if (type == "Gompertz"){

			initVals <- getInitial(log(value) ~ SSgompertz(variable, Asym, b2, b3), data = df_tmp[df_tmp$MouseID %in% unique(df_tmp$MouseID)[2],])
			gomp <- nlme(value ~ (Y_asym)*exp(-exp(-k*(variable - t_m))), data = df_tmp, random = Y_asym + k + t_m ~ 1, fixed = Y_asym + k + t_m ~ 1, start = c(Y_asym = round(initVals[[1]]), k = 1, t_m = 30), groups = ~MouseID, control = nlmeControl(msMaxIter = 1500))
			tmp <- nls2(log(value) ~ (Y_asym)*exp(-exp(-k*(variable - t_m))), data = df_tmp[df_tmp$MouseID %in% unique(df_tmp$MouseID)[2],], start = c(Y_asym = 18, k = 1, t_m = 30), algo = "brute")
			param_est <-  data.frame(matrix(rep(gomp$coefficients$fixed, length(unique(df_tmp$MouseID))), nrow = length(unique(df_tmp$MouseID)), ncol = 3, byrow = TRUE) + gomp$coefficients$random$MouseID)
			param_est$MouseID <- as.factor(rownames(param_est))
			param_est$Strain <- rep(Strains[s],nrow(param_est))
			param_list[[s]] <- param_est  


		} else if (type == "Logistic"){

			#initVals <- getInitial(value ~ SSgompertz(variable, Asym, b2, b3), data = df_tmp)
			logist <- nlme(value ~ (Y_asym)/(1 + exp(-k*(variable - t_m))), data = df_tmp, random = Y_asym + k + t_m ~ 1, fixed = Y_asym + k + t_m ~ 1, start = c(Y_asym = 12, k = 1, t_m = 30), groups = ~MouseID)
			param_est <-  data.frame(logist$coefficients$fixed + logist$coefficients$random$MouseID)
			param_est$MouseID <- as.factor(rownames(param_est))
			param_est$Strain <- rep(Strains[s],nrow(param_est))
			param_list[[s]] <- param_est 

		}

	}
	return(param_list)

}