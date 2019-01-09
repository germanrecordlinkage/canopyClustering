canopyClustering.matchFiles <- function(filenameA, filenameB, minSimilarity, looseThreshold, tightThreshold, method = "tanimoto") {
        if (method == "tanimoto") {
                methodId <- 1;
        } else if (method == "tanimotoXOR") {
                methodId <- 2;
        } else {
                stop("Unknown method");
        }
	result <- .Call(ccMatchFilesCall, filenameA, filenameB, minSimilarity, looseThreshold, tightThreshold, methodId)
	return(data.frame(result))
}
