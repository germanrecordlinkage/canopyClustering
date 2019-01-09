canopyClustering.statistics <- function() {
	options("scipen"=16)
	result <- .Call(ccStatisticsCall)
	return(data.frame(result))
}
