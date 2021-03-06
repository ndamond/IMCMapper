% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IMCMapperFunctions.R
\name{ClusterDisplay}
\alias{ClusterDisplay}
\title{ClusterDisplay}
\usage{
ClusterDisplay(x, im.list, id.col, image.col, cluster.col, color.pal = NULL,
  meta, path.col, file.col, slide.col, borders = FALSE, bcg = c(0, 1, 1),
  leg = TRUE, ...)
}
\arguments{
\item{x}{your dataset in the data table format}

\item{im.list}{vector containing the names of the images to be displayed (in the same format as in the image.col in your dataset)}

\item{id.col}{name of the column in your dataset that contains the cell IDs}

\item{image.col}{name of the column in your dataset that contains the image names}

\item{cluster.col}{name of the column in your dataset that contains the clusters. The column can contain discrete or continuous values}

\item{color.pal}{color palette: if nothing is provided , a default palette will be used.\cr
- option 1: name of a color palette in the RColorBrewer ("Set1", "Dark2", "Accent", "Paired", "RdYlBu", "Greys", etc.).\cr
For the full list of palettes see: display.brewer.all() (library(RColorBrewer))\cr
- option 2: name of a color palette in the viridis package ("viridis", "magma", "inferno", "plasma", "cividis")\cr
- option 3: a custom palette in the form of a vector containing Hex color codes (e.g., c("#DC050C", "#FB8072", ...) )}

\item{meta}{path to the "Image.csv" file extracted from CellProfiler}

\item{path.col}{name of the column in the "Image.csv" file that contains the path to the masks (usually starts with "ObjectsPathName_")}

\item{file.col}{name of the column in the "Image.csv" file that contains the names of the masks (usually starts with "ObjectsFileName_")}

\item{slide.col}{name of the column in the "Image.csv" file that contains the image names (usually starts with "Metadata_")}

\item{borders}{(optional) TRUE/FALSE: should cell borders be displayed?}

\item{bcg}{(optional) background/contrast/gamma to apply to the image}

\item{leg}{TRUE/FALSE: should a color legend be displayed. The legend will be displayed as a separate image (default = TRUE)\cr
if "save" is set to TRUE, the legend will be saved together with the pictures as a .png image}

\item{...}{any other argument taken by the "legend" function (see ?legend)}

\item{method}{'raster' or 'browser' option from the EBImage package to display the images in a separate (new tab) JavaScript viewer ('browser')\cr
or directly in R ('raster') (default = 'raster')}

\item{save}{TRUE/FALSE: should the images be saved (default = FALSE)}

\item{path}{path to the folder in which the images should be saved. If the folder does not exist it will be created (default = working directory)}

\item{type}{'png', 'jpeg' or 'tiff': format of the saved images (default = 'png')}

\item{legend}{vector containing the names of the elements displayed in the legend}

\item{fill}{vector containing the colors to be displayed in the legend}

\item{border}{vector containing the color of the legend borders}

\item{cex}{size of the text of the legend (character expansion)}

\item{label}{TRUE/FALSE: should the name of the image be displayed on the image (default = TRUE)\cr
Only works with the 'raster' method. The label is also not displayed on saved images}

\item{label.txt}{text to be display on the image label. The same text will be displayed on all images (default = name of the image)}

\item{label.color}{color of the image label displayed on the image}

\item{label.cex}{size of the text of the image name label (character expansion)}

\item{label.x}{x position of the image name label}

\item{label.y}{y position of the image name label}
}
\value{
currently, this function either prints or saves images but nothing is explicitly returned
}
\description{
Displays cell clusters or cell types
}
\keyword{cell}
\keyword{cluster}
\keyword{type}
