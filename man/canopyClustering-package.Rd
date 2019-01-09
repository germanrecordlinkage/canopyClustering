\name{canopyClustering-package} 
\alias{canopyClustering-package}
\alias{canopyClustering}
\docType{package}
\title{
Similarity Matching by Sorted Neighbourhood Method
}
\description{
This package contains functions for searching similar objects
in two files of data.

Currently the package expects the data to consist of binary vectors
(fingerprints). The measure of similarity is calculated by the
Tanimoto-Coefficient.

The package's data structure and algorithms are based on
Canopy Clustering and XOR-folding
as discribed in the article referenced below.
}
\details{
\tabular{ll}{
Package: \tab canopyClustering\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2014-04-01\cr
License: \tab GPL-3\cr
}
}
\author{
Stefani Kamphausen, Markus Koetter, iQbis consulting GmbH

Maintainer: Prof. Dr. Rainer Schnell <rainer.schnell@uni-due.de>
}
\references{
Andrew McCallum, Kamal Nigam, and Lyle H. Ungar:
Efficient clustering of high-dimensional data sets
with application to reference matching
The Sixth ACM SIGKDD International Conference on
Knowledge Discovery and Data Mining KDD-2000, page 169--178. (2000)

Baldi P, Hirschberg DS, Nasr RJ:
Speeding Up Chemical Database Searches
Using a Proximity Filter Based on the Logical Exclusive OR.
Journal of Chemical Information and Modeling 2008
}
\keyword{misc}
