make_sort_column <- function(data, sort_col, sort_metric, FUN = median){
    temp <- aggregate(formula(paste0(sort_metric,'~',sort_col)), data=data, FUN=FUN)
    return(factor(unlist(subset(data, select = sort_col)), levels=subset(temp, select = sort_col)[order(subset(temp,
     select = sort_metric)),]))
}