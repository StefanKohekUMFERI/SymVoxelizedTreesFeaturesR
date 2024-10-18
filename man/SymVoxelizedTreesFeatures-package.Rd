\name{SymVoxelizedTreesFeatures-package}
\alias{SymVoxelizedTreesFeatures-package}
\alias{SymVoxelizedTreesFeatures}
\docType{package}
\title{
  Calculate symmetry aware features of voxelized XYZ trees for tree species classification
}
\description{
  This package computes features from voxelized trees in XYZ data format. The features measure symmetry through the crown by using IOUs of symmetric trees. 
}
\details{
  
}
\author{
Štefan Kohek,
Damjan Strnad

Maintainer: Štefan Kohek <stefan.kohek@um.si>
}
\references{
  
}
\keyword{ package }
\seealso{

}
\examples{
  GP = new(UGridAndPoints);
  GP$ReadFileToGridR(fileName, swapYZ=TRUE, cellSize=0.5)
  pts=GP$GetPointsR();
  f=GP$CalculateFeaturesR()
}
