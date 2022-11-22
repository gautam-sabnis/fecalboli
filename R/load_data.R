#' Load and preprocess data 
#' 
#' Loads the data and preprocess the features. 
#' @param csv a csv file containing data per stride or grooming data
#' @param type load either grooming or fecal boli data
#' @return a list containing two data frames
#' 1) df_equal: equally spaced time points
#' 2) df_unequal: un-equally spaced time points
#' @examples 
#' load_data(data = "fecalboli.csv", type = "fboli")
#' load_data(data = "grooming.tsv", type = "grooming")
#' @export 

load_data <- function(csv, type){

	if (type == "grooming"){
		tsv <- csv
		grooming_data <- read.delim(paste0(tsv),stringsAsFactors = TRUE)
		grooming_data <- cbind(grooming_data, bin_avg_5 = apply(grooming_data[,names(grooming_data) %in% c(sapply(seq(0,4), function(x) paste0('bin_',x)))],1,mean,na.rm=TRUE), bin_avg_20 = apply(grooming_data[,names(grooming_data) %in% c(sapply(seq(0,19), function(x) paste0('bin_',x)))],1,mean,na.rm=TRUE),bin_avg_55 = apply(grooming_data[,names(grooming_data) %in% c(sapply(seq(0,54), function(x) paste0('bin_',x)))],1,mean,na.rm=TRUE))
		grooming_data <- reshape(grooming_data[,names(grooming_data) %in% c('Strain','MouseID','Sex','CoatColor','measure',sapply(c(5,20,55), function(x) paste0('bin_avg_',x)), 'bin_avg')], idvar = c('Strain','MouseID','Sex','CoatColor'), timevar = c('measure'), direction = 'wide')


		} else {
			FBdata <- read.csv(paste0(csv), stringsAsFactors = TRUE)
			data <- FBdata[, which(names(FBdata) %in% c("NetworkFilename", "MouseID", "Strain", "Sex", "Weight", sapply(seq(60), function(x) paste0("minuteMax", x))))]
#			data <- data[, c("Strain", "MouseID", "Sex", "Weight", sapply(seq(60), function(x) paste0("minuteMax", x)))]
#			names(data)[names(data) %in% c(sapply(seq(60), function(x) paste0("minuteMax", x)))] <- as.factor(sapply(seq(60), function(x) x))
			data <- cbind(data, TestAge=as.Date(FBdata$TestDate, format='%m/%d/%Y') - as.Date(FBdata$DOB,format='%m/%d/%Y'))#days
			data <- data[!duplicated(data$MouseID),]
			data$MouseID <- as.factor(data$MouseID)
			Strains <- unique(data$Strain)

			df_equal <- data.frame()
			df_unequal <- data.frame()
			pb <- txtProgressBar(min = 0, max = length(Strains), style = 3)

			for (s in seq(Strains)){

  				df <- data[data$Strain %in% Strains[s],]
  				tmp <- list()
  				invisible(suppressMessages({(lapply(seq(nrow(df)), function(x) {
  				tmp[[x]] <<- reshape2::melt(df[df$NetworkFilename == unique(df$NetworkFilename)[x], which(names(df) %in% c('MouseID','Sex','Strain','CoatColor', sapply(seq(60), function(x) paste0("minuteMax",x))))])}))}))
  				df_melt <- do.call(rbind,tmp)
  				df_melt$value <- as.numeric(df_melt$value)
  				df_melt$variable <- as.integer(df_melt$variable)
  				df_melt$Sex <- as.factor(rep(df$Sex,each = 60)) 
  				df_melt$Age <- as.numeric(rep(df$TestAge, each = 60))
  				df_melt$Weight <- rep(df$Weight, each = 60)
  				df_equal <- rbind(df_equal, df_melt)

  				tmp <- list()
  				invisible(lapply(seq(length(unique(df_melt$MouseID))), function(x) {
  					tmpp <- df_melt[df_melt$MouseID == unique(df_melt$MouseID)[x],]
    				tmp[[x]] <<- tmpp[!duplicated(tmpp$value),]}))
  					df_melt2 <- do.call(rbind,tmp)
  					df_melt2$variable <- as.integer(df_melt2$variable)
  					df_unequal <- rbind(df_unequal, df_melt2)

  				setTxtProgressBar(pb, s)
  				close(pb)
  			}


			df_equal <- groupedData(value ~ variable|MouseID, data = df_equal, order.groups=FALSE)
			df_unequal <- groupedData(value ~ variable|MouseID, data = df_unequal, order.groups=FALSE)
			df_equal <- df_equal[complete.cases(df_equal),]
			df_equal <- df_equal[!df_equal$value == 0, ]
			df_unequal <- df_unequal[complete.cases(df_unequal),]
			df_unequal <- df_unequal[!df_unequal$value == 0, ]

			return(list(df_equal, df_unequal))
		}
	

}