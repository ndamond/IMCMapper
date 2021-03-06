% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IMCMapperFunctions.R
\name{DisplayImages}
\alias{DisplayImages}
\title{DisplayImages}
\usage{
DisplayImages(im.1, im.2 = NULL, im.3 = NULL, images, bcg1 = c(0, 1, 1),
  bcg2 = c(0, 1, 1), bcg3 = c(0, 1, 1), masks = NULL, ...)
}
\arguments{
\item{im.1}{1st image to display (green)}

\item{im.2}{2nd image to display (red)}

\item{im.3}{3rd image to display (blue)}

\item{images}{list of images in data format}

\item{bcg1}{brightness / contrast / gamma for the 1st image (default = c(0,1,1))}

\item{bcg2}{brightness / contrast / gamma for the 2nd image (default = c(0,1,1))}

\item{bcg3}{brightness / contrast / gamma for the 3rd image (default = c(0,1,1))}

\item{masks}{list of image data with cell borders}

\item{...}{additional parameters for the label and legend}
}
\value{
nothing is explicitly returned
}
\description{
Display 1 to 3 combined channel images as RGB
}
\keyword{RGB}
\keyword{display}
\keyword{image}
