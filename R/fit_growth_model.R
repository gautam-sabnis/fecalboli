#' Fit growth model 
#' 
#' Fit a growth model to the data 
#' @param data data frame containing equally/unequally spaced time 
#' points
#' @param Strain Strain of interest
#' @param type type of growth curve
#' @return Estimated growth curve parameters 
#' @examples 
#' growth_model(data = df_equal, Strain = "C57BL/6NJ", type = "Gompertz")
#' @export 

fit_growth_model <- function(data, Strain, type, shrinkage, plot){
    setTimeLimit(100)
	df_Strain <- data[data$Strain %in% Strain, ]
	mIDs <- unique(df_Strain$MouseID)
    df_params <- data.frame(matrix(0,length(mIDs),3))
    names(df_params) <- c("Asym", "k", "t_m")
    if (shrinkage == TRUE){

        tryCatch(
            expr = {
                initVals <<- getInitial(value ~ SSgompertz(variable, Asym, b2, b3), data = df_Strain)
                fit <<- suppressWarnings(nlme(value ~ SSgompertz(variable, Asym, b2, b3), data = df_Strain, fixed = list(Asym+b2+b3 ~ 1), random = pdDiag(Asym+b2+b3 ~ 1), groups = ~MouseID, control = nlmeControl(maxIter = 5000, msMaxIter = 1500)))
                fit_fe <<- data.frame(Asym = rep(fixef(fit)[1], length(unique(df_Strain$MouseID))), b2 = rep(fixef(fit)[2], length(unique(df_Strain$MouseID))), b3 = rep(fixef(fit)[3], length(unique(df_Strain$MouseID))))
                fit_beta <<- fit_fe + ranef(fit)
                df_params <<- data.frame(Asym = fit_beta$Asym, k = -log(fit_beta$b3), t_m = -log(fit_beta$b2)/log(fit_beta$b3))
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
    df_params$CoatColor <- df_Strain[!duplicated(df_Strain$MouseID), "CoatColor"]
    df_params <- df_params[,c("Strain", "MouseID", "Sex", "CoatColor","Asym", "k", "t_m")]


    if (plot == FALSE){
        return(list(df_params,Asym_covariate))
    } else {
        if (shrinkage == TRUE){
                df_fit <- augPred(fit,level=0:1)
                for (x in seq(unique(df_fit$.groups))){
                    df_fit[df_fit$.group == unique(df_fit$.group)[x], "MouseID"] <- (rep(x,table(df_fit$.groups)[x]))
                }
                df_fit$MouseID <- df_fit$.groups
                p <- ggplot(df_fit, aes(x = variable, y=value)) + geom_point(data = df_fit[df_fit$.type == "original",], aes(x = variable, y=value), alpha = 0.5) + geom_line(data = df_fit[df_fit$.type == "predict.MouseID",], aes(x = variable, y = value), color = 'magenta') + geom_line(data = df_fit[df_fit$.type == "predict.fixed",], aes(x=variable, y = value), color = 'blue') + facet_wrap(.~MouseID, scales="free") + labs(x = 'Time (in mins)', y = 'Value') + ggtitle(paste0(Strain)) + theme_bw()
                } else {
                    df_fit <- df_Strain[,c("variable", "MouseID", "value")]
                    df_fit$.type <- as.factor("original")
                    #df_pred <- data.frame()
                    

                }
            }
        }