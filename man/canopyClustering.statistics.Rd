\name{canopyClustering.statistics}
\alias{canopyClustering.statistics}
\title{
Print Statistics of last Search
}
\description{
This function returns the statistic results of the last call to matchFiles.
}
\usage{
canopyClustering.statistics()
}
\value{
The function returns a data.frame with three columns:
\item{Checkpoint}{
  this column contains the statistics parameter
}
\item{Count}{
  this column contains the corresponding value
}
\item{Percentage}{
  this column contains the value Count in relation to the total sum of comparisons
}
}
\seealso{
\code{\link{canopyClustering.matchFiles}}
}
\examples{
## get name of example files with fingerprints in package directory

fileA <- file.path(path.package("canopyClustering"), "extdata/A.csv")
fileB <- file.path(path.package("canopyClustering"), "extdata/B.csv")

## search all prints from file A in file B and print result

print(canopyClustering.matchFiles(fileA, fileB, 0.8, 0.7, 0.8))

## print statistics

canopyClustering.statistics()

}
\keyword{misc}
