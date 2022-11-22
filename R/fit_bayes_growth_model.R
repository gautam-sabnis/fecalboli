#' Fit growth model in a Bayesian framework
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

fit_bayes_growth_model <- function(data, Strain, type){
	df_Strain <- data[data$Strain %in% Strain, ]
    mIDs <- names(table(df_Strain$MouseID)[which(as.numeric(table(df_Strain$MouseID)) < 4)])
    df_Strain <- df_Strain[!df_Strain$MouseID %in% mIDs, ]
    df_Strain$MouseID <- droplevels(df_Strain$MouseID)
    df_Strain$Sex <- as.factor(as.numeric(df_Strain$Sex))
    if (type == "Gompertz"){
    	fit <- brm(bf(value ~ theta1*exp(-exp(-theta2*(variable - theta3))), theta1 ~ Sex + (1|MouseID), theta2 ~ Sex + (1|MouseID), theta3 ~ Sex + (1|MouseID), nl = TRUE), prior = c(prior(normal(0,4), class = "b", nlpar = "theta1"), prior(normal(0,4), class = "b", nlpar = "theta2"), prior(normal(0,4), class = "b", nlpar = "theta3")), data = df_Strain, family = gaussian(), iter = 2000, cores = 4, backend = "cmdstanr", threads = threading(2), silent = 2, control = list(max_treedepth = 10))
    } else if (type == "Logistic"){
    	fit <- brm(bf(value ~ theta1/(1 + exp(-theta2*(variable - theta3))), theta1 ~ Sex + (1|MouseID), theta2 ~ Sex + (1|MouseID), theta3 ~ Sex + (1|MouseID), nl = TRUE), prior = c(prior(normal(0,4), class = "b", nlpar = "theta1"), prior(normal(0,4), class = "b", nlpar = "theta2"), prior(normal(0,4), class = "b", nlpar = "theta3")), data = df_Strain, family = gaussian(), iter = 2000, cores = 4, backend = "cmdstanr", threads = threading(2), silent = 2, control = list(max_treedepth = 10))
    } else {
    	fit <- brm(bf(value ~ theta1*(1 + psi*exp(-theta2*(variable - theta3)))^(-1/psi), theta1 ~ Sex + (1|MouseID), theta2 ~ (1|MouseID), theta3 ~ (1|MouseID), psi ~ 1, nl = TRUE), prior = c(prior(normal(0,4), class = "b", nlpar = "theta1"), prior(normal(0,4), class = "b", nlpar = "theta2"), prior(normal(0,4), class = "b", nlpar = "theta3"), prior(normal(0,4), class = "b", nlpar = "psi")), data = df_Strain, family = gaussian(), iter = 3000, cores = 4, backend = "cmdstanr", threads = threading(2), refresh = 0, silent = 2)
    }

    return(fit)
}