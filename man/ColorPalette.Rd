% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IMCMapperFunctions.R
\name{ColorPalette}
\alias{ColorPalette}
\title{ColorPalette}
\usage{
ColorPalette(clusters, color.pal = NULL)
}
\arguments{
\item{clusters}{names of the clusters in the data set}

\item{color.pal}{(optional) name of a palette in the "RColorBrewer" package or in the "viridis" package, or a custom palette provided by the user}
}
\value{
cluster.palette a palette containing the same number of colors as the number of clusters
}
\description{
Generate a color palette depending on user input and on the number of clusters to display
}
\keyword{cluster}
\keyword{color}
\keyword{palette}
