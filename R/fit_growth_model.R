#' Fit growth model 
#' 
#' Fit a growth model to the data 
#' @param data data frame containing equally/unequally spaced time 
#' points
#' @param Strain Strain of interest
#' @param type type of growth curve
#' @return Estimated growth curve parameters 
#' @examples 
#' fit_growth_model(data = df_equal, Strain = "C57BL/6NJ", type = "Gompertz")
#' @export 

fit_growth_model <- function(data, Strain, type, shrinkage, plot){
    setTimeLimit(60)
	df_Strain <- data[data$Strain %in% Strain, ]
	mIDs <- unique(df_Strain$MouseID)
    df_params <- data.frame(matrix(0,length(mIDs),3))
    df_params_strain <- data.frame(matrix(0,1,7))
    names(df_params) <- c("Asym", "k", "t_m")
    if (Strain == "MOLF/EiJ"){
      mIDs <- names(which(table(df_Strain$MouseID) > 3));
      df_Strain <- df_Strain[df_Strain$MouseID %in% mIDs, ]
    } else if (Strain %in% c("CZECHII/EiJ", "NZB/BlNJ")) {
      mIDs <- names(which(table(df_Strain$MouseID) > 5));
      df_Strain <- df_Strain[df_Strain$MouseID %in% mIDs, ]
    } else {
        df_Strain <- df_Strain
    }
    df_Strain$Sex <- as.factor(df_Strain$Sex)
    if (shrinkage == TRUE){

        tryCatch(
            expr = {
                initVals <- getInitial(value ~ SSgompertz(variable, Asym, b2, b3), data = df_Strain)
                fit <- suppressWarnings(nlme(value ~ SSgompertz(variable, Asym, b2, b3), data = df_Strain, fixed = Asym + b2 + b3 ~ Sex, random = pdDiag(Asym+b2+b3 ~ 1), groups = ~MouseID, control = nlmeControl(maxIter = 5000, msMaxIter = 1500)))
                fit_fe <- data.frame(Asym = rep(fixef(fit)[1], length(unique(df_Strain$MouseID))), b2 = rep(fixef(fit)[2], length(unique(df_Strain$MouseID))), b3 = rep(fixef(fit)[3], length(unique(df_Strain$MouseID))))
                fit_beta <- fit_fe + ranef(fit)
                fit_fe_se <- sqrt(diag(vcov(fit)))
                df_params <- data.frame(Asym = fit_beta$Asym, k = -log(fit_beta$b3), t_m = -log(fit_beta$b2)/log(fit_beta$b3))
                df_params_strain <- data.frame(Strain = Strain, Asym = fixef(fit)[1], k = -log(fixef(fit)[3]), t_m = -log(fixef(fit)[2])/log(fixef(fit)[3]), Asym_se = fit_fe_se[1], k_se = (fit_fe_se[3]/fixef(fit)[3]), t_m_se = (fit_fe_se[2]/abs(fixef(fit)[2]*log(fixef(fit)[3]))))
                    },
            error = function(e){cat("ERROR :",conditionMessage(e), "\n")}
            )
    } else {

        for (m in 1:length(mIDs)){
            skip_to_next <- FALSE
            df_tmpp <- df_Strain[df_Strain$MouseID %in% mIDs[m], ]
            df_tmpp$MouseID <- droplevels(df_tmpp$MouseID)
            tryCatch(
                expr = {
                    if (type == "Gompertz"){
                        initVals <- getInitial(value ~ SSgompertz(variable, Asym, b2, b3), data = df_tmpp)
                        fit <- nls(value ~ SSgompertz(variable, Asym, b2, b3), data = df_tmpp, control = list(maxiter = 1000))
                        df_params[m,] <- c(summary(fit)$coefficients[1], -log(summary(fit)$coefficients[3]), -log(summary(fit)$coefficients[2])/log(summary(fit)$coefficients[3]))
                    }
                    else if (type == "Logistic"){

                    }
                    }, 
                error = function(e) {
                    df_params[m,] <- c(rep(NA,3))
                    skip_to_next <- TRUE
                })
            if (skip_to_next) {next}
        }

    }
    
    df_params$Strain <- as.factor(rep(paste0(Strain), nrow(df_params)))
    df_params$MouseID <- mIDs
    df_params$Sex <- df_Strain[!duplicated(df_Strain$MouseID), "Sex"]
    #df_params$CoatColor <- df_Strain[!duplicated(df_Strain$MouseID), "CoatColor"]
    df_params <- df_params[,c("Strain", "MouseID", "Sex", "Asym", "k", "t_m")]
    #df_params_strain <- df_params_strain[,c("Strain", "Asym", "k", "t_m", "Asym_se", "k_se", "t_m_se")]


    return(list(df_params, df_params_strain))
}


