\name{canopyClustering.matchFiles}
\alias{canopyClustering.matchFiles}
\title{
Match multiple Objects 
}
\description{
This function matches objects of two input files according to their similarity
using Canopy Clustering algorithm.
}
\usage{
canopyClustering.matchFiles(filenameA, filenameB, minSimilarity, looseThreshold, tightThreshold, method = "tanimoto")
}
\arguments{
\item{filenameA}{
  a character string containing the filename of the first input file
}
\item{filenameB}{
  a character string containing the filename of the second input file
}
\item{minSimilarity}{
  a numeric value giving the lower bound of similarity coefficient to search for
}
\item{looseThreshold}{
  a numeric value describing the loose threshold for canopy clustering
}
\item{tightThreshold}{
  a numeric value describing the tight threshold for canopy clustering
}
\item{method}{
  an optional string parameter, that defines the method to work with; currently
  two methods are implemented: "tanimto" (default, reading objects as fingerprints,
  comparing objects by tanimoto coefficient) and "tanimotoXOR" (the same as
  "tanimoto" extended by a 128-bit XOR-hash-filter strategy for fast tanimoto
  estimation and comparison)
}
}
\value{
The function returns a data.frame with three columns:
\item{query}{
  this column contains the query object-ids
}
\item{fingerprint}{
  this column contains the matching object-ids
}
\item{tanimoto}{
  this column contains the corresponding similarity coefficients
}
}
\examples{
## get name of example files with fingerprints in package directory

fileA <- file.path(path.package("canopyClustering"), "extdata/A.csv")
fileB <- file.path(path.package("canopyClustering"), "extdata/B.csv")

## search all prints from file A in file B and print result

print(canopyClustering.matchFiles(fileA, fileB, 0.8, 0.7, 0.8))

}
\keyword{misc}
