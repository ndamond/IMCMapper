## MAIN FUNCTIONS


#' CellMap
#'
#' Display 1 to 3 cell populations
#'
#' @param x your dataset in the data table format
#' @param id.col name of the column in your dataset that contains the cell IDs
#' @param cell1 vector containing the cell IDs of the 1st cell population to be displayed
#' @param cell2 (optional) vector containing the cell IDs of the 2nd cell population to be displayed
#' @param cell3 (optional) vector containing the cell IDs of the 3rd cell population to be displayed
#' @param meta path to the "Image.csv" file extracted from CellProfiler
#' @param path.col name of the column in the "Image.csv" file that contains the path to the masks (usually starts with "ObjectsPathName_")
#' @param file.col name of the column in the "Image.csv" file that contains the names of the masks (usually starts with "ObjectsFileName_")
#' @param slide.col name of the column in the "Image.csv" file that contains the image names (usually starts with "Metadata_")
#' @param borders (optional) TRUE/FALSE: should cell borders be displayed?
#' @param method 'raster' or 'browser' option from the EBImage package to display the images in a separate (new tab) JavaScript viewer ('browser')\cr
#' or directly in R ('raster') (default = 'raster')
#' @param save TRUE/FALSE: should the images be saved (default = FALSE)
#' @param path path to the folder in which the images should be saved. If the folder does not exist it will be created (default = working directory)
#' @param type 'png', 'jpeg' or 'tiff': format of the saved images (default = 'png')
#' @param leg TRUE/FALSE: should a color legend be displayed. The legend will be displayed as a separate image (default = TRUE)\cr
#' if "save" is set to TRUE, the legend will be saved together with the pictures as a .png image
#' @param legend vector containing the names of the elements displayed in the legend
#' @param fill vector containing the colors to be displayed in the legend
#' @param border vector containing the color of the legend borders
#' @param cex size of the text of the legend (character expansion)
#' @param ... any other argument taken by the "legend" function (see ?legend)
#' @param label TRUE/FALSE: should the name of the image be displayed on the image (default = TRUE)\cr
#' Only works with the 'raster' method. The label is also not displayed on saved images
#' @param label.txt text to be display on the image label. The same text will be displayed on all images (default = name of the image)
#' @param label.color color of the image label displayed on the image
#' @param label.cex size of the text of the image name label (character expansion)
#' @param label.x x position of the image name label
#' @param label.y y position of the image name label
#' @keywords population
#' @return currently, this function either prints or saves images but nothing is explicitly returned
#' @import EBImage
#' @import data.table
#' @export
#'
CellMap <- function (x, id.col,
                     cell1, cell2=NULL, cell3=NULL,
                     meta, path.col, file.col, slide.col,
                     borders = TRUE, leg = TRUE, ...
) {

  # Get the masks
  masks <- getMasks(image.metadata = meta, path.col = path.col, file.col = file.col, slide.col = slide.col)

  # Recover the names of the images from cell ids and read the masks
  im.list <- unique(c(sub("_[^_]+$", "", unique(cell1)), sub("_[^_]+$", "", unique(cell2)), sub("_[^_]+$", "", unique(cell3))))
  my.images <- readMasks(masks = masks, im.list = im.list)

  # Get the cells
  cell.sel.1 <- getCells(images = my.images, pop = unique(cell1))
  cell.sel.2 <- getCells(images = my.images, pop = unique(cell2))
  cell.sel.3 <- getCells(images = my.images, pop = unique(cell3))

  # Get the cell borders
  if (borders==TRUE)
    mask.borders <- getBorders(images = my.images)
  else
    mask.borders = NULL

  # Print out legend
  if (leg == TRUE){
    nb.items = length(Filter(Negate(is.null), list(cell1, cell2, cell3)))
    makeLegend(leg.names = c("Population 1", "Population 2", "Population 3"), title = "Cell populations",
               fun="CellMap", nb.items=nb.items, ...)
  }

  # Display the images
  DisplayImages(im.1=cell.sel.1, im.2=cell.sel.2, im.3=cell.sel.3, images=im.list, masks=mask.borders, fun="CellMap", ...)
}




#' ChannelHM
#'
#' Heatmap 1 to 3 channels
#'
#' @param x your dataset in the data table format. The dataset should have one column that contains the channel names (channel.col) and one column that contains the values for all channels (counts.col)
#' @param id.col name of the column in your dataset that contains the cell IDs
#' @param counts.col name of the column in your dataset that contains the values to display (e.g. counts)
#' @param channel.col name of the column in your dataset that contains the channel names
#' @param image.col name of the column in your dataset that contains the image names
#' @param im.list vector containing the names of the images to be displayed (in the same format as in the image.col in your dataset)
#' @param ch1 name of the 1st channel to be displayed (same name as in the "channel.col" column)
#' @param ch2 (optional) name of the 2nd channel to be displayed (same name as in the "channel.col" column)
#' @param ch3 (optional) name of the 3rd channel to be displayed (same name as in the "channel.col" column)
#' @param bcg1 (optional) background/contrast/gamma to apply to the 1st channel (ch1)
#' @param bcg2 (optional) background/contrast/gamma to apply to the 2nd channel (ch2)
#' @param bcg3 (optional) background/contrast/gamma to apply to the 3rd channel (ch3)
#' @param meta path to the "Image.csv" file extracted from CellProfiler
#' @param path.col name of the column in the "Image.csv" file that contains the path to the masks (usually starts with "ObjectsPathName_")
#' @param file.col name of the column in the "Image.csv" file that contains the names of the masks (usually starts with "ObjectsFileName_")
#' @param slide.col name of the column in the "Image.csv" file that contains the image names (usually starts with "Metadata_")
#' @param borders (optional) TRUE/FALSE: should cell borders be displayed?
#' @param method 'raster' or 'browser' option from the EBImage package to display the images in a separate (new tab) JavaScript viewer ('browser')\cr
#' or directly in R ('raster') (default = 'raster')
#' @param save TRUE/FALSE: should the images be saved (default = FALSE)
#' @param path path to the folder in which the images should be saved. If the folder does not exist it will be created (default = working directory)
#' @param type 'png', 'jpeg' or 'tiff': format of the saved images (default = 'png')
#' @param leg TRUE/FALSE: should a color legend be displayed. The legend will be displayed as a separate image (default = TRUE)\cr
#' if "save" is set to TRUE, the legend will be saved together with the pictures as a .png image
#' @param legend vector containing the names of the elements displayed in the legend
#' @param fill vector containing the colors to be displayed in the legend
#' @param border vector containing the color of the legend borders
#' @param cex size of the text of the legend (character expansion)
#' @param ... any other argument taken by the "legend" function (see ?legend)
#' @param label TRUE/FALSE: should the name of the image be displayed on the image (default = TRUE)\cr
#' Only works with the 'raster' method. The label is also not displayed on saved images
#' @param label.txt text to be display on the image label. The same text will be displayed on all images (default = name of the image)
#' @param label.color color of the image label displayed on the image
#' @param label.cex size of the text of the image name label (character expansion)
#' @param label.x x position of the image name label
#' @param label.y y position of the image name label
#' @keywords channel heatmap
#' @return currently, this function either prints or saves images but nothing is explicitly returned
#' @import EBImage
#' @import data.table
#' @export

ChannelHM <- function (x, id.col, counts.col, channel.col,
                       image.col, im.list,
                       ch1, ch2=NULL, ch3=NULL,
                       bcg1 = c(0,4,2), bcg2 = c(0,4,2), bcg3 = c(0,4,2),
                       meta, path.col, file.col, slide.col,
                       borders = FALSE, leg = TRUE, ...
) {

  # Get the masks
  masks <- getMasks(image.metadata = meta, path.col = path.col, file.col = file.col, slide.col = slide.col)

  # Get the images
  my.images <- readMasks(masks = masks, im.list = im.list)

  # Get the number of cells per image
  nb.cells <- getCellNumber(images = my.images)

  # Subset relevant images and columns
  setkeyv(x, image.col)
  x <- unique(x[as.character(im.list), .(image = get(image.col), id = get(id.col), channel = get(channel.col), counts = get(counts.col))])

  # Get the counts and ids, and create the heatmaps for each channel
  counts.ids.1 <- getCounts(x = x[channel==ch1], ch = ch1, im.list = im.list)
  value.image.1 <- CreateHeatmap(images = my.images, nb.cells = nb.cells, counts.ids = counts.ids.1)

  if (!is.null(ch2)) {
    counts.ids.2 <- getCounts(x = x[channel==ch2], ch = ch2, im.list = im.list)
    value.image.2 <- CreateHeatmap(images = my.images, nb.cells = nb.cells, counts.ids = counts.ids.2)
  }
  else
    value.image.2 <- NULL

  if (!is.null(ch3)) {
    counts.ids.3 <- getCounts(x = x[channel==ch3], ch = ch3, im.list = im.list)
    value.image.3 <- CreateHeatmap(images = my.images, nb.cells = nb.cells, counts.ids = counts.ids.3)
  }
  else
    value.image.3 <- NULL

  # Get the cell borders
  if (borders==TRUE)
    mask.borders <- getBorders(images = my.images)
  else
    mask.borders = NULL

  # Print out legend
  if (leg == TRUE){
    nb.items = length(Filter(Negate(is.null), list(ch1, ch2, ch3)))
    makeLegend(leg.names = c(ch1, ch2, ch3),
               leg.fill = switch(as.character(is.null(ch2) & is.null(ch3) & borders==F), 'TRUE' = c('white'),  'FALSE' = c('green', 'red', 'blue')),
               title="Channels", fun="ChannelHM", nb.items=nb.items, ...)
  }

  # Display the images
  DisplayImages(im.1=value.image.1, im.2=value.image.2, im.3=value.image.3, images=im.list,
                bcg1=bcg1, bcg2=bcg2, bcg3=bcg3, masks=mask.borders, fun="ChannelHM", ...)
}



#' CellsAndCounts
#'
#' Heatmap 1 channel + highlight the borders of 1 to 3 cell populations
#'
#' @param x your dataset in the data table format. \cr
#' The dataset should have one column that contains the channel names (channel.col) and one column that contains the values for all channels (counts.col)
#' @param id.col name of the column in your dataset that contains the cell IDs
#' @param image.col name of the column in your dataset that contains the image names
#' @param counts.col name of the column in your dataset that contains the values to display (e.g. counts)
#' @param channel.col name of the column in your dataset that contains the channel names
#' @param ch name of the channel to be displayed (same name as in the "channel.col" column)
#' @param cell1 vector containing the cell IDs of the 1st cell population to be displayed
#' @param cell2 (optional) vector containing the cell IDs of the 2nd cell population to be displayed
#' @param cell3 (optional) vector containing the cell IDs of the 3rd cell population to be displayed
#' @param meta path to the "Image.csv" file extracted from CellProfiler
#' @param path.col name of the column in the "Image.csv" file that contains the path to the masks (usually starts with "ObjectsPathName_")
#' @param file.col name of the column in the "Image.csv" file that contains the names of the masks (usually starts with "ObjectsFileName_")
#' @param slide.col name of the column in the "Image.csv" file that contains the image names (usually starts with "Metadata_")
#' @param bcg (optional) background/contrast/gamma to apply to the displayed channel (ch)
#' @param method 'raster' or 'browser' option from the EBImage package to display the images in a separate (new tab) JavaScript viewer ('browser')\cr
#' or directly in R ('raster') (default = 'raster')
#' @param save TRUE/FALSE: should the images be saved (default = FALSE)
#' @param path path to the folder in which the images should be saved. If the folder does not exist it will be created (default = working directory)
#' @param type 'png', 'jpeg' or 'tiff': format of the saved images (default = 'png')
#' @param leg TRUE/FALSE: should a color legend be displayed. The legend will be displayed as a separate image (default = TRUE)\cr
#' if "save" is set to TRUE, the legend will be saved together with the pictures as a .png image
#' @param legend vector containing the names of the elements displayed in the legend
#' @param fill vector containing the colors to be displayed in the legend
#' @param border vector containing the color of the legend borders
#' @param cex size of the text of the legend (character expansion)
#' @param ... any other argument taken by the "legend" function (see ?legend)
#' @param label TRUE/FALSE: should the name of the image be displayed on the image (default = TRUE)\cr
#' Only works with the 'raster' method. The label is also not displayed on saved images
#' @param label.txt text to be display on the image label. The same text will be displayed on all images (default = name of the image)
#' @param label.color color of the image label displayed on the image
#' @param label.cex size of the text of the image name label (character expansion)
#' @param label.x x position of the image name label
#' @param label.y y position of the image name label
#' @keywords channel heatmap population
#' @return currently, this function either prints or saves images but nothing is explicitly returned
#' @import EBImage
#' @import data.table
#' @export

CellsAndCounts <- function (x,
                            id.col, image.col, counts.col, channel.col,
                            ch, cell1, cell2=NULL, cell3=NULL,
                            meta, path.col, file.col, slide.col,
                            bcg = c(0,4,0.5), leg = TRUE, ...
) {

  # Get the masks
  masks <- getMasks(image.metadata = meta, path.col = path.col, file.col = file.col, slide.col = slide.col)

  # Recover the names of the images from cell ids and read the masks
  im.list <- unique(c(sub("_[^_]+$", "", unique(cell1)), sub("_[^_]+$", "", unique(cell2)), sub("_[^_]+$", "", unique(cell3))))
  my.images <- readMasks(masks = masks, im.list = im.list)

  # Get the number of cells per image
  nb.cells <- getCellNumber(images = my.images)

  # Get the cells
  cell.sel.1 <- getCells(images = my.images, pop = unique(cell1))
  cell.sel.2 <- getCells(images = my.images, pop = unique(cell2))
  cell.sel.3 <- getCells(images = my.images, pop = unique(cell3))

  # Subset relevant images and columns
  setkeyv(x, image.col)
  x <- unique(x[as.character(im.list), .(image = get(image.col), id = get(id.col), channel = get(channel.col), counts=get(counts.col))])

  # Get the counts and ids, and create the heatmaps
  counts.ids <- getCounts(x = x[channel==ch], ch = ch, im.list = im.list)
  value.image <- CreateHeatmap(images = my.images, nb.cells = nb.cells, counts.ids = counts.ids)

  # Get cells borders
  mask.borders <- getBorders(images = my.images)

  # Print out legend
  if (leg == TRUE){
    nb.items = length(Filter(Negate(is.null), list(cell1, cell2, cell3)))
    makeLegend(leg.names = c(ch, "Population 1", "Population 2", "Population 3"),
               leg.fill = c('white', 'grey', 'grey', 'grey'), leg.border = c('black', 'green', 'red', 'blue'),
               title="Channel / Populations", fun="CellsAndCounts", nb.items=nb.items+1, ...)
  }

  # Display the images
  DisplayCellsAndCounts(value.image=value.image, cell.sel.1=cell.sel.1, cell.sel.2=cell.sel.2, cell.sel.3=cell.sel.3,
                        im.list=im.list, cb=mask.borders, bcg=bcg, fun="CellsAndCounts", ...)
}




#' ClusterDisplay
#'
#' Displays cell clusters or cell types
#'
#' @param x your dataset in the data table format
#' @param im.list vector containing the names of the images to be displayed (in the same format as in the image.col in your dataset)
#' @param id.col name of the column in your dataset that contains the cell IDs
#' @param image.col name of the column in your dataset that contains the image names
#' @param cluster.col name of the column in your dataset that contains the clusters. The column can contain discrete or continuous values
#' @param color.pal   color palette: if nothing is provided , a default palette will be used.\cr
#' - option 1: name of a color palette in the RColorBrewer ("Set1", "Dark2", "Accent", "Paired", "RdYlBu", "Greys", etc.).\cr
#' For the full list of palettes see: display.brewer.all() (library(RColorBrewer))\cr
#' - option 2: name of a color palette in the viridis package ("viridis", "magma", "inferno", "plasma", "cividis")\cr
#' - option 3: a custom palette in the form of a vector containing Hex color codes (e.g., c("#DC050C", "#FB8072", ...) )
#' @param meta path to the "Image.csv" file extracted from CellProfiler
#' @param path.col name of the column in the "Image.csv" file that contains the path to the masks (usually starts with "ObjectsPathName_")
#' @param file.col name of the column in the "Image.csv" file that contains the names of the masks (usually starts with "ObjectsFileName_")
#' @param slide.col name of the column in the "Image.csv" file that contains the image names (usually starts with "Metadata_")
#' @param borders (optional) TRUE/FALSE: should cell borders be displayed?
#' @param bcg (optional) background/contrast/gamma to apply to the image
#' @param method 'raster' or 'browser' option from the EBImage package to display the images in a separate (new tab) JavaScript viewer ('browser')\cr
#' or directly in R ('raster') (default = 'raster')
#' @param save TRUE/FALSE: should the images be saved (default = FALSE)
#' @param path path to the folder in which the images should be saved. If the folder does not exist it will be created (default = working directory)
#' @param type 'png', 'jpeg' or 'tiff': format of the saved images (default = 'png')
#' @param leg TRUE/FALSE: should a color legend be displayed. The legend will be displayed as a separate image (default = TRUE)\cr
#' if "save" is set to TRUE, the legend will be saved together with the pictures as a .png image
#' @param legend vector containing the names of the elements displayed in the legend
#' @param fill vector containing the colors to be displayed in the legend
#' @param border vector containing the color of the legend borders
#' @param cex size of the text of the legend (character expansion)
#' @param ... any other argument taken by the "legend" function (see ?legend)
#' @param label TRUE/FALSE: should the name of the image be displayed on the image (default = TRUE)\cr
#' Only works with the 'raster' method. The label is also not displayed on saved images
#' @param label.txt text to be display on the image label. The same text will be displayed on all images (default = name of the image)
#' @param label.color color of the image label displayed on the image
#' @param label.cex size of the text of the image name label (character expansion)
#' @param label.x x position of the image name label
#' @param label.y y position of the image name label
#' @keywords cluster cell type
#' @return currently, this function either prints or saves images but nothing is explicitly returned
#' @import EBImage
#' @import data.table
#' @export

ClusterDisplay <- function (x, im.list,
                            id.col, image.col, cluster.col,
                            color.pal = NULL,
                            meta, path.col, file.col, slide.col,
                            borders = FALSE, bcg = c(0,1,1), leg = TRUE, ...
) {

  # Get the masks
  masks <- getMasks(image.metadata = meta, path.col = path.col, file.col = file.col, slide.col = slide.col)

  # Get cluster list (looks in the entire dataset so that colors remain consistent even if not all clusters are present in the selected images)
  clusters <- as.factor(unique(x[, get(cluster.col)]))

  # Subset relevant images and columns
  setkeyv(x, image.col)
  x <- unique(x[as.character(im.list), .(image = get(image.col), id = get(id.col), cluster = get(cluster.col))])

  # Get the images
  my.images <- readMasks(masks = masks, im.list = im.list)

  # Get the number of cells per image
  nb.cells <- getCellNumber(images = my.images)

  # Get cluster list
  # clusters <- as.factor(unique(x[, cluster]))

  # Match each cluster with a specific color
  cluster.palette <- ColorPalette(clusters = clusters, color.pal = color.pal)
  x <- MatchColorCluster(x = x, clusters = clusters, cluster.palette = cluster.palette, return.x = TRUE)

  # Attribute a color to each cell based on the cluster the cell belongs to
  cluster.map <- getClusters(x = x, im.list = im.list, ids = FALSE)
  cluster.ids <- getClusters(x = x, im.list = im.list, ids = TRUE)

  # Create images
  cluster.image <- CreateClusterHeatmap(images = my.images, nb.cells = nb.cells, cluster.map = cluster.map, cluster.ids = cluster.ids)

  # Print out legend
  if (leg == TRUE){
    if (length(clusters) < 75){
      color.table <- MatchColorCluster(x = x, clusters = clusters, cluster.palette = cluster.palette, return.x = FALSE)
      makeLegend(leg.names = color.table$cluster, leg.fill = color.table$cluster.palette,
                 ncol=ceiling(length(clusters)/12), title=cluster.col, fun="ClusterDisplay", nb.items=length(clusters), ...)
    }
    else
      makeLegend(leg.names = c(signif(max(x$cluster), 3), rep('', 8), signif(min(x$cluster), 3)),
                 leg.fill = cluster.palette[round(rev(seq(1, length(cluster.palette), length.out=10)))],
                 y.intersp=0.5, title=paste0(cluster.col, "\n"), fun="ClusterDisplay", nb.items=10, ...)
  }

  # Display color images
  DisplayColorImages(cluster.image = cluster.image, images = my.images, borders = borders, bcg = bcg, fun = "ClusterDisplay", ...)
}




#' IMCDisplay
#'
#' Displays raw IMC images
#'
#' @param x your dataset in the data table format
#' @param image.col name of the column in your dataset that contains the image names
#' @param im.list vector containing the names of the images to be displayed (in the same format as in the image.col in your dataset)
#' @param ch1 name of the 1st channel to be displayed (same name as in the "channel.col" column)
#' @param ch2 (optional) name of the 2nd channel to be displayed (same name as in the "channel.col" column)
#' @param ch3 (optional) name of the 3rd channel to be displayed (same name as in the "channel.col" column)
#' @param bcg1 (optional) background/contrast/gamma to apply to the 1st channel (ch1)
#' @param bcg2 (optional) background/contrast/gamma to apply to the 2nd channel (ch2)
#' @param bcg3 (optional) background/contrast/gamma to apply to the 3rd channel (ch3)
#' @param meta path to the "Image.csv" file extracted from CellProfiler
#' @param path.col name of the column in the "Image.csv" file that contains the path to the masks (usually starts with "ObjectsPathName_")
#' @param file.col name of the column in the "Image.csv" file that contains the names of the masks (usually starts with "ObjectsFileName_")
#' @param slide.col name of the column in the "Image.csv" file that contains the image names (usually starts with "Metadata_")
#' @param borders (optional) TRUE/FALSE: should cell borders be displayed?
#' @param method 'raster' or 'browser' option from the EBImage package to display the images in a separate (new tab) JavaScript viewer ('browser')\cr
#' or directly in R ('raster') (default = 'raster')
#' @param save TRUE/FALSE: should the images be saved (default = FALSE)
#' @param path path to the folder in which the images should be saved. If the folder does not exist it will be created (default = working directory)
#' @param type 'png', 'jpeg' or 'tiff': format of the saved images (default = 'png')
#' @param leg TRUE/FALSE: should a color legend be displayed. The legend will be displayed as a separate image (default = TRUE)\cr
#' if "save" is set to TRUE, the legend will be saved together with the pictures as a .png image
#' @param legend vector containing the names of the elements displayed in the legend
#' @param fill vector containing the colors to be displayed in the legend
#' @param border vector containing the color of the legend borders
#' @param cex size of the text of the legend (character expansion)
#' @param ... any other argument taken by the "legend" function (see ?legend)
#' @param label TRUE/FALSE: should the name of the image be displayed on the image (default = TRUE)\cr
#' Only works with the 'raster' method. The label is also not displayed on saved images
#' @param label.txt text to be display on the image label. The same text will be displayed on all images (default = name of the image)
#' @param label.color color of the image label displayed on the image
#' @param label.cex size of the text of the image name label (character expansion)
#' @param label.x x position of the image name label
#' @param label.y y position of the image name label
#' @keywords IMC channel raw
#' @return currently, this function either prints or saves images but nothing is explicitly returned
#' @import EBImage
#' @import data.table
#' @export

IMCDisplay <- function (x,
                        im.list, image.col,
                        ch1,  ch2 = NULL, ch3 = NULL,
                        bcg1 = c(0,0.5,1), bcg2 = c(0,0.5,1), bcg3 = c(0,0.5,1),
                        meta, path.col, file.col = NULL, slide.col,
                        borders=FALSE, leg = TRUE, ...
) {

  # Get the image paths
  pics <- getIMCImages(image.metadata = meta, ch1 = ch1, ch2 = ch2, ch3 = ch3, path.col = path.col, file.col = file.col, slide.col = slide.col, borders = borders)

  # Get the images
  my.image.1 <- readImages(pics = pics, im.list = im.list, chan = 'pic1')

  my.image.2 <- ch2
  if(!is.null(ch2)) my.image.2 <- readImages(pics = pics, im.list = im.list, chan = 'pic2')

  my.image.3 <- ch3
  if(!is.null(ch3)) my.image.3 <- readImages(pics = pics, im.list = im.list, chan = 'pic3')

  # Get the cell borders
  if (borders==TRUE)
    mask.borders <- getBorders(images = readMasks(masks = pics, im.list = im.list))
  else
    mask.borders = NULL

  # Print out legend
  if (leg == TRUE){
    nb.items = length(Filter(Negate(is.null), list(ch1, ch2, ch3)))
    makeLegend(leg.names = c(ch1, ch2, ch3),
               leg.fill = switch(as.character(is.null(ch2) & is.null(ch3) & borders==F), 'TRUE' = c('white'),  'FALSE' = c('green', 'red', 'blue')),
               title="Channel", fun="IMCDisplay", nb.items=nb.items, ...)
  }

  # Display the images
  DisplayImages(im.1=my.image.1, im.2=my.image.2, im.3=my.image.3, images=im.list,
                bcg1=bcg1, bcg2=bcg2, bcg3=bcg3, masks=mask.borders, fun="IMCDisplay", ...)
}




#' CellsandIMC
#'
#' Display 1 IMC image + highlight the borders of 1 to 3 cell populations
#'
#' @param x your dataset in the data table format.\cr
#' The dataset should have one column that contains the channel names (channel.col) and one column that contains the values for all channels (counts.col)
#' @param id.col name of the column in your dataset that contains the cell IDs
#' @param image.col name of the column in your dataset that contains the image names
#' @param counts.col name of the column in your dataset that contains the values to display (e.g. counts)
#' @param channel.col name of the column in your dataset that contains the channel names
#' @param ch channel image to display. The name must match a single channel image in the image folder
#' @param cell1 vector containing the cell IDs of the 1st cell population to be displayed
#' @param cell2 (optional) vector containing the cell IDs of the 2nd cell population to be displayed
#' @param cell3 (optional) vector containing the cell IDs of the 3rd cell population to be displayed
#' @param meta path to the "Image.csv" file extracted from CellProfiler
#' @param path.col name of the column in the "Image.csv" file that contains the path to the masks (usually starts with "ObjectsPathName_")
#' @param file.col name of the column in the "Image.csv" file that contains the names of the masks (usually starts with "ObjectsFileName_")
#' @param slide.col name of the column in the "Image.csv" file that contains the image names (usually starts with "Metadata_")
#' @param bcg (optional) background/contrast/gamma to apply to the displayed channel (ch)
#' @param method 'raster' or 'browser' option from the EBImage package to display the images in a separate (new tab) JavaScript viewer ('browser')\cr
#' or directly in R ('raster') (default = 'raster')
#' @param save TRUE/FALSE: should the images be saved (default = FALSE)
#' @param path path to the folder in which the images should be saved. If the folder does not exist it will be created (default = working directory)
#' @param type 'png', 'jpeg' or 'tiff': format of the saved images (default = 'png')
#' @param leg TRUE/FALSE: should a color legend be displayed. The legend will be displayed as a separate image (default = TRUE)\cr
#' if "save" is set to TRUE, the legend will be saved together with the pictures as a .png image
#' @param legend vector containing the names of the elements displayed in the legend
#' @param fill vector containing the colors to be displayed in the legend
#' @param border vector containing the color of the legend borders
#' @param cex size of the text of the legend (character expansion)
#' @param ... any other argument taken by the "legend" function (see ?legend)
#' @param label TRUE/FALSE: should the name of the image be displayed on the image (default = TRUE)\cr
#' Only works with the 'raster' method. The label is also not displayed on saved images
#' @param label.txt text to be display on the image label. The same text will be displayed on all images (default = name of the image)
#' @param label.color color of the image label displayed on the image
#' @param label.cex size of the text of the image name label (character expansion)
#' @param label.x x position of the image name label
#' @param label.y y position of the image name label
#' @keywords IMC population channel
#' @return currently, this function either prints or saves images but nothing is explicitly returned
#' @import EBImage
#' @import data.table
#' @export

CellsandIMC <- function (x,
                         id.col, image.col, counts.col, channel.col,
                         ch, cell1, cell2=NULL, cell3=NULL,
                         meta, path.col, file.col, slide.col,
                         bcg = c(0,4,2), leg = TRUE, ...
) {

  # Get the image and mask paths
  pics <- getIMCImages(image.metadata = meta, ch1 = ch,
                       path.col = path.col, file.col = file.col, slide.col = slide.col, borders=TRUE)


  # Recover the names of the images from cell ids and read the masks
  im.list <- unique(c(sub("_[^_]+$", "", unique(cell1)), sub("_[^_]+$", "", unique(cell2)), sub("_[^_]+$", "", unique(cell3))))
  my.images <- readMasks(masks = pics, im.list = im.list)

  # Get the number of cells per image
  nb.cells <- getCellNumber(images = my.images)

  # Get the IMC images
  IMC.images.1 <- readImages(pics = pics, im.list = im.list, chan = 'pic1')

  # Get the cells
  cell.sel.1 <- getCells(images = my.images, pop = unique(cell1))
  cell.sel.2 <- getCells(images = my.images, pop = unique(cell2))
  cell.sel.3 <- getCells(images = my.images, pop = unique(cell3))

  # Get cells borders
  mask.borders <- getBorders(images = my.images)

  # Print out legend
  if (leg == TRUE){
    nb.items = length(Filter(Negate(is.null), list(cell1, cell2, cell3)))
    makeLegend(leg.names = c(ch, "Population 1", "Population 2", "Population 3"),
               leg.fill = c('white', 'grey', 'grey', 'grey'), leg.border = c('black', 'green', 'red', 'blue'),
               title="Channel / Populations", fun="CellsandIMC", nb.items=nb.items+1, ...)
  }

  # Display the images
  DisplayCellsAndCounts(value.image=IMC.images.1, cell.sel.1=cell.sel.1, cell.sel.2=cell.sel.2, cell.sel.3=cell.sel.3,
                        im.list=im.list, cb=mask.borders, bcg=bcg, fun="CellsandIMC", ...)
}




#' IMCandCells
#'
#' Display 1 to 3 IMC images + highlight the borders of 1 cell population
#'
#' @param x your dataset in the data table format.\cr
#' The dataset should have one column that contains the channel names (channel.col) and one column that contains the values for all channels (counts.col)
#' @param id.col name of the column in your dataset that contains the cell IDs
#' @param image.col name of the column in your dataset that contains the image names
#' @param cell1 vector containing the cell IDs of the cell population to be displayed
#' @param ch1 1st channel image to display. The name must match a single channel image in the image folder
#' @param ch2 (optional) 2nd channel image to display. The name must match a single channel image in the image folder
#' @param ch3 (optional) 3rd channel image to display. The name must match a single channel image in the image folder
#' @param bcg1 background/contrast/gamma to apply to the 1st channel (ch1)
#' @param bcg2 (optional) background/contrast/gamma to apply to the 2nd channel (ch2)
#' @param bcg3 (optional) background/contrast/gamma to apply to the 3rd channel (ch3)
#' @param meta path to the "Image.csv" file extracted from CellProfiler
#' @param path.col name of the column in the "Image.csv" file that contains the path to the masks (usually starts with "ObjectsPathName_")
#' @param file.col name of the column in the "Image.csv" file that contains the names of the masks (usually starts with "ObjectsFileName_")
#' @param slide.col name of the column in the "Image.csv" file that contains the image names (usually starts with "Metadata_")
#' @param method 'raster' or 'browser' option from the EBImage package to display the images in a separate (new tab) JavaScript viewer ('browser')\cr
#' or directly in R ('raster') (default = 'raster')
#' @param save TRUE/FALSE: should the images be saved (default = FALSE)
#' @param path path to the folder in which the images should be saved. If the folder does not exist it will be created (default = working directory)
#' @param type 'png', 'jpeg' or 'tiff': format of the saved images (default = 'png')
#' @param leg TRUE/FALSE: should a color legend be displayed. The legend will be displayed as a separate image (default = TRUE)\cr
#' if "save" is set to TRUE, the legend will be saved together with the pictures as a .png image
#' @param legend vector containing the names of the elements displayed in the legend
#' @param fill vector containing the colors to be displayed in the legend
#' @param border vector containing the color of the legend borders
#' @param cex size of the text of the legend (character expansion)
#' @param ... any other argument taken by the "legend" function (see ?legend)
#' @param label TRUE/FALSE: should the name of the image be displayed on the image (default = TRUE)\cr
#' Only works with the 'raster' method. The label is also not displayed on saved images
#' @param label.txt text to be display on the image label. The same text will be displayed on all images (default = name of the image)
#' @param label.color color of the image label displayed on the image
#' @param label.cex size of the text of the image name label (character expansion)
#' @param label.x x position of the image name label
#' @param label.y y position of the image name label
#' @keywords IMC population channel
#' @return currently, this function either prints or saves images but nothing is explicitly returned
#' @import EBImage
#' @import data.table
#' @export

IMCandCells <- function (x,
                         id.col, image.col,
                         cell1, ch1, ch2=NULL, ch3=NULL,
                         bcg1 = c(0,0.5,1), bcg2 = c(0,0.5,1), bcg3 = c(0,0.5,1),
                         meta, path.col, file.col, slide.col, leg = TRUE, ...
) {

  # Get the image and mask paths
  pics <- getIMCImages(image.metadata = meta, ch1 = ch1, ch2 = ch2, ch3 = ch3, path.col = path.col, file.col = file.col, slide.col = slide.col, borders=TRUE)

  # Recover the names of the images from cell ids and read the masks
  im.list <- unique(c(sub("_[^_]+$", "", unique(cell1))))
  my.images <- readMasks(masks = pics, im.list = im.list)

  # Get the IMC images
  IMC.images.1 <- readImages(pics = pics, im.list = im.list, chan = 'pic1')

  IMC.images.2 <- ch2
  if(!is.null(ch2)) IMC.images.2 <- readImages(pics = pics, im.list = im.list, chan = 'pic2')

  IMC.images.3 <- ch3
  if(!is.null(ch3)) IMC.images.3 <- readImages(pics = pics, im.list = im.list, chan='pic3')

  # Get the cells
  cell.sel.1 <- getCells(images = my.images, pop = unique(cell1))

  # Get the cell borders
  mask.borders <- getBorders(images = my.images)

  # Print out legend
  if (leg == TRUE){
    nb.items = length(Filter(Negate(is.null), list(ch1, ch2, ch3)))
    makeLegend(leg.names = c("Population 1", ch1, ch2, ch3),
               leg.fill = c('white', 'green', 'red', 'blue'),
               title="Population / Channels", fun="IMCandCells", nb.items=nb.items+1, ...)
  }

  # Display the images
  DisplayIMCandCells(value.image.1=IMC.images.1, value.image.2=IMC.images.2, value.image.3=IMC.images.3, cell.sel.1=cell.sel.1,
                     im.list=im.list, cb=mask.borders, bcg1=bcg1, bcg2=bcg2, bcg3=bcg3, fun="IMCandCells", ...)
}



## INTERNAL FUNCTIONS


#' getMasks
#'
#' Recover the paths to the mask tiff files from the "Image.csv" file
#'
#' @param image.metadata path to the "Image.csv" file extracted from CellProfiler
#' @param path.col name of the column in the "Image.csv" file that contains the path to the masks
#' @param file.col name of the column in the "Image.csv" file that contains the names of the masks
#' @param slide.col name of the column in the "Image.csv" file that contains the image names
#' @keywords mask path
#' @return masks a data table with links to the masks
#' @import data.table
#' @export

getMasks <- function (image.metadata, path.col, file.col, slide.col) {

  if (is.character(image.metadata))
    image.metadata <- data.table::fread(image.metadata, header = T, stringsAsFactors = F)

  masks <- image.metadata[, .(path = get(path.col), file = get(file.col), slide = get(slide.col))]
  masks <- masks[, files := file.path(path, file)]

  masks <- masks[dir.exists(path)]
}



#' getIMCImages
#'
#' Recover the paths to the IMC tiff files from the "Image.csv" file and the image names provided by the user
#'
#' @param image.metadata path to the "Image.csv" file extracted from CellProfiler
#' @param ch1 name of the 1st channel to display
#' @param ch2 name of the 2nd channel to display
#' @param ch3 name of the 3rd channel to display
#' @param path.col name of the column in the "Image.csv" file that contains the path to the masks
#' @param file.col name of the column in the "Image.csv" file that contains the names of the masks
#' @param slide.col name of the column in the "Image.csv" file that contains the image names
#' @param borders TRUE/FALSE should cell borders be displayed?
#' @keywords image mask path
#' @return a data table with links to the images and masks to be displayed
#' @import data.table
#' @export

getIMCImages <- function (image.metadata, ch1, ch2=NULL, ch3=NULL, path.col, file.col=NULL, slide.col, borders=FALSE) {

  if (is.character(image.metadata))
    image.metadata <- data.table::fread(image.metadata, header = T, stringsAsFactors = F)

  if(borders==F){
    pics <- image.metadata[, .(path = get(path.col), slide = get(slide.col))]
  } else if(!is.null(file.col)){
    pics <- image.metadata[, .(path = get(path.col), file = get(file.col), slide = get(slide.col))]
  }

  pics <- pics[dir.exists(path)]
  if (nrow(pics) == 0) print("Directories listed in the 'Image.csv' file do not exist")

  pics[, pic1 := file.path(path, list.files(path, pattern = ch1)[1])]

  if(!is.null(ch2)) pics[, pic2 := file.path(path, list.files(path, pattern = ch2)[1])]

  if(!is.null(ch3)) pics[, pic3 := file.path(path, list.files(path, pattern = ch3)[1])]

  if(borders==T) pics[, files := file.path(path, file)]

  return(pics)
}



#' readMasks
#'
#' Read mask images (16-bit tiff) and store the image data in a list
#'
#' @param masks data table with links to the masks
#' @param im.list vector with names of images to be displayed
#' @keywords read masks images
#' @return list of images in data format
#' @import EBImage
#' @import data.table
#' @export

readMasks <- function(masks, im.list) {

  my.images <- list()

  for (i in seq_along(im.list)) {
    image.path <- masks[slide == im.list[i], files]
    my.images[[i]] <- EBImage::readImage(image.path)
    my.images[[i]] <- EBImage::imageData(my.images[[i]])*65535
  }

  names(my.images) <- im.list
  return(my.images)
}



#' readImages
#'
#' Read IMC images and store them in a list
#'
#' @param pics data table with links to the images and masks to display
#' @param im.list vector with names of images to be displayed
#' @param chan name of the column in pics containing the image to read
#' @keywords read image mask
#' @return list of images in data format
#' @import EBImage
#' @import data.table
#' @export

readImages <- function(pics, im.list, chan='pic1') {

  my.images <- list()
  for (i in seq_along(im.list)) {
    image.path <- pics[slide == im.list[i], get(chan)]
    my.images[[i]] <- EBImage::readImage(image.path)
  }

  names(my.images) <- im.list
  return(my.images)
}



#' getCells
#'
#' From a matrix of pixels corresponding to the mask image: pixels belonging to the cell population to be displayed = 1; other pixels = 0
#'
#' @param images list of images in data format
#' @param pop ids of cells to be displayed
#' @keywords image cells
#' @return list of images in data format where cells to be displayed are indicated
#' @export

getCells <- function (images, pop) {

  cell.sel <- list()

  # Generate a matrix with image names in the first column and cell ids in the second column
  ids <- matrix(nrow=length(pop), ncol=2)

  for (i in seq_along(pop)){
    ids[i,1] <- sub("_[^_]+$", "", pop[i])
    ids[i,2] <- sub(".*_", "", pop[i])
  }

  # Identify the pixels corresponding to cells to be displayed
  for (i in seq_along(images)) {
    image.id <- names(images[i])
    cell.id <- ids[ids[,1] == image.id, 2]

    cell.sel[[i]] <- as.numeric(images[[i]] %in% as.numeric(cell.id))
    cell.sel[[i]] <- matrix(cell.sel[[i]], nrow = nrow(images[[i]]))
  }

  names(cell.sel) <- names(images)
  return(cell.sel)
}



#' getBorders
#'
#' Determine cell borders from the mask using a filter function
#'
#' @param images list of images in data format
#' @keywords cell border mask
#' @return list of images in data format where cell borders are indicated
#' @import EBImage
#' @export

getBorders <- function (images) {

  la <- matrix(1, nc=3, nr=3)
  la[2,2] <- -8

  borders <- list()
  for (i in seq_along(images)){
    borders[[i]] <- EBImage::filter2(images[[i]], la)
  }

  names(borders) <- names(images)
  return(borders)
}



#' getCounts
#'
#' Get the counts value for all cells in a given channel
#'
#' @param x data table with cell ids, channels and counts
#' @param ch channel to be displayed
#' @param im.list vector with names of images to be displayed
#' @keywords channel counts
#' @return data table with cell ids and respective counts for the channel to be displayed
#' @import data.table
#' @export

getCounts <- function (x, ch, im.list) {

  channel.counts <- list()

  for (i in seq_along(im.list)){
    channel.counts[[i]] <- x[image==im.list[i] & channel==ch, .(counts, id)]
  }

  names(channel.counts) <- im.list
  return(channel.counts)
}



#' ColorPalette
#'
#' Generate a color palette depending on user input and on the number of clusters to display
#'
#' @param clusters names of the clusters in the data set
#' @param color.pal (optional) name of a palette in the "RColorBrewer" package or in the "viridis" package, or a custom palette provided by the user
#' @keywords color palette cluster
#' @return cluster.palette a palette containing the same number of colors as the number of clusters
#' @export

ColorPalette <- function(clusters, color.pal=NULL){

  # List of palettes in the RColorBrewer and in the viridis packages
  brewer.pal.palettes <- c("BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", "RdYlGn", "Spectral", "Accent", "Dark2", "Paired",
                           "Pastel1", "Pastel2", "Set1", "Set2", "Set3", "Blues", "BuGn", "BuPu", "GnBu", "Greens", "Greys", "Oranges",
                           "OrRd", "PuBu", "PuBuGn", "PuRd", "Purples", "RdPu", "Reds", "YlGn", "YlGnBu", "YlOrBr", "YlOrRd")

  viridis.palettes <- c("viridis", "magma", "inferno", "plasma", "cividis")

  if (is.null(color.pal)){

    # default palette for displaying clusters
    if (length(clusters) < 75){
      default.palette <- c("#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#BF5B17", "#666666", "#1B9E77", "#D95F02", "#7570B3",
                           "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#999999", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C",
                           "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4", "#FED9A6",
                           "#FFFFCC", "#E5D8BD", "#FDDAEC", "#F2F2F2", "#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4", "#E6F5C9", "#FFF2AE", "#F1E2CC",
                           "#CCCCCC", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5",
                           "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072",
                           "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F"
      )
      cluster.palette <- default.palette[1:length(clusters)]

    } else {

      # default palette for displaying continuous values
      require(RColorBrewer)
      cols <- rev(brewer.pal(11, "Spectral"))
      cluster.palette <- colorRampPalette(cols)(length(clusters))
    }

  } else if (color.pal[1] %in% brewer.pal.palettes){

    # if user input matches one of the RColorBrewer palettes
    require(RColorBrewer)
    color.nb <- as.data.frame(brewer.pal.info)[color.pal,]$maxcolors

    cols <- brewer.pal(color.nb, color.pal)
    cluster.palette <- colorRampPalette(cols)(length(clusters))

  } else if (color.pal[1] %in% viridis.palettes){

    # if user input matches one of the viridis palettes
    require(viridisLite)
    cluster.palette <- get(color.pal)(length(clusters), alpha = 1, begin = 0, end = 1, direction = 1)

  } else if (length(color.pal) < length(clusters)) {

    # if user inputs a palette that does not have enough colors, the provided palette is resampled
    cat(paste("Warning: not enough colors in the color vector provided, your custom palette was resampled\n",
              "Number of clusters:", length(clusters), "\n",
              "Number of colors provided:", length(color.pal)))
    cluster.palette <- colorRampPalette(color.pal)(length(clusters))

  } else {

    # custom palette provided by the user
    cluster.palette <- color.pal[1:length(clusters)]
  }

  return(cluster.palette)
}



#' getClusters
#'
#' Match each cell id to a color based on the cluster to which the cell belong
#'
#' @param x data table with cell ids, image names, colors and clusters
#' @param im.list vector with names of images to be displayed
#' @param ids should ids (TRUE) or colors (FALSE) be returned
#' @keywords cluster cell id
#' @return a data table with either the cell ids or the colors to be displayed
#' @import data.table
#' @export

getClusters <- function (x, im.list, ids=FALSE) {

  cluster.map <- list()

  if (ids == FALSE){
    for (i in seq_along(im.list))
      cluster.map[[i]] <- x[image==im.list[i], color]
  } else if (ids == TRUE){
    for (i in seq_along(im.list))
      cluster.map[[i]] <- x[image==im.list[i], id]
  }

  names(cluster.map) <- im.list
  return(cluster.map)
}



#' getCellNumber
#'
#' Returns a list of cell numbers (one list element per image)
#'
#' @param images list of images in data format
#' @keywords cell number
#' @return list with number of cells per image
#' @export

getCellNumber <- function(images) {
  nb.cells <- c()

  for (i in seq_along(images))
    nb.cells[i] <- max(images[[i]]) +1

  names(nb.cells) <- names(images)
  return(nb.cells)
}



#' MatchColorCluster
#'
#' Match colors and clusters
#'
#' @param x data table with colors and clusters to display
#' @param clusters vector with cluster names
#' @param cluster.palette color palette to use, with one color per cluster
#' @param return.x should the data table (TRUE) or the lookup table (FALSE) be returned
#' @keywords color cluster
#' @return data table with clusters and matched color for each cluster
#' @import data.table
#' @export

MatchColorCluster <- function(x, clusters, cluster.palette, return.x = TRUE){

  lookup <- data.table::as.data.table(cbind(cluster = levels(clusters), cluster.palette))

  if(return.x == FALSE)
    return(lookup)
  else {
    x[, color := lookup$cluster.palette[match(unlist(x$cluster), lookup$cluster)]]
    return(x)
  }
}



#' CreateHeatmap
#'
#' Create a heatmap based on the counts
#'
#' @param images list of images in data format
#' @param nb.cells list with number of cells per image
#' @param counts.ids data table with number of counts for all cells
#' @keywords cell counts
#' @return list of images in data format where number of counts is indicated
#' @export

CreateHeatmap <- function (images, nb.cells, counts.ids) {

  value.image <- list()

  for (i in seq_along(images)) {
    channel.counts <- counts.ids[[i]]$counts
    channel.ids <- counts.ids[[i]]$id

    pop <- as.numeric(sub(".*_", "", channel.ids))
    counts.vector <- integer(nb.cells[[i]])
    counts.vector[pop + 1] <- channel.counts

    value.image[[i]] <- counts.vector[images[[i]]+1]
    value.image[[i]] <- matrix(value.image[[i]], nrow = nrow(images[[i]]))
  }

  names(value.image) <- names(images)
  return(value.image)
}



#' CreateClusterHeatmap
#'
#' Create a heatmap based on the clusters
#'
#' @param images list of images in data format
#' @param nb.cells list with number of cells per image
#' @param cluster.map colors
#' @param cluster.ids clusters
#' @keywords color clusters
#' @return list of images in data format where cell clusters have been replaced by colors
#' @export

CreateClusterHeatmap <- function (images, nb.cells, cluster.map, cluster.ids) {

  cluster.image <- list()

  for (i in seq_along(images)) {
    id.list <- as.numeric(sub(".*_", "", cluster.ids[[names(images[i])]]))
    cluster.vector <- integer(nb.cells[[names(images[i])]])
    cluster.vector[id.list + 1] <- cluster.map[[names(images[i])]]

    cluster.image[[i]] <- cluster.vector[images[[i]]+1]
    cluster.image[[i]] <- replace(cluster.image[[i]], cluster.image[[i]] == "0", "#000000")
    cluster.image[[i]] <- matrix(cluster.image[[i]], nrow = nrow(images[[i]]))
  }

  names(cluster.image) <- names(images)
  return(cluster.image)
}



#' DisplayImages
#'
#' Display 1 to 3 combined channel images as RGB
#'
#' @param im.1 1st image to display (green)
#' @param im.2 2nd image to display (red)
#' @param im.3 3rd image to display (blue)
#' @param images list of images in data format
#' @param bcg1 brightness / contrast / gamma for the 1st image (default = c(0,1,1))
#' @param bcg2 brightness / contrast / gamma for the 2nd image (default = c(0,1,1))
#' @param bcg3 brightness / contrast / gamma for the 3rd image (default = c(0,1,1))
#' @param masks list of image data with cell borders
#' @param ... additional parameters for the label and legend
#' @keywords display image RGB
#' @return nothing is explicitly returned
#' @import EBImage
#' @export

DisplayImages <- function(im.1, im.2=NULL, im.3=NULL, images, bcg1 = c(0,1,1), bcg2 = c(0,1,1), bcg3 = c(0,1,1), masks=NULL, ...) {

  borderImg <- list()

  for (i in seq_along(images)) {

    if (!is.null(masks)){
      borderImg[[i]] <- EBImage::Image(EBImage::combine(masks[[i]]))
      border.comb <- EBImage::combine(borderImg[[i]])
    }

    # If only one channel is to be displayed without cell borders, it will be displayed in white
    if(is.null(im.2) & is.null(im.3)){
      img.comb <- EBImage::combine((((im.1[[i]] + bcg1[1]) * bcg1[2]) ^ bcg1[3]))

      # If cell borders are displayed (white color), the channel image will be displayed in green
      if (!is.null(masks)){
        img.comb <-   EBImage::combine(border.comb, img.comb + border.comb, border.comb)
        img.comb <- EBImage::Image(img.comb)
        colorMode(img.comb) <- Color
      }
    }

    # Two channels: red and green
    if(!is.null(im.2) & is.null(im.3)){
      img.comb <- EBImage::combine((((im.2[[i]] + bcg2[1]) * bcg2[2]) ^ bcg2[3]),
                                   (((im.1[[i]] + bcg1[1]) * bcg1[2]) ^ bcg1[3]))

      if (!is.null(masks)){
        img.comb <- EBImage::combine((((im.2[[i]] + bcg2[1]) * bcg2[2]) ^ bcg2[3]) + border.comb,
                                     (((im.1[[i]] + bcg1[1]) * bcg1[2]) ^ bcg1[3]) + border.comb, border.comb)
      }
      img.comb <- EBImage::Image(img.comb)
      colorMode(img.comb) <- Color
    }

    # Three channels: red, green and blue
    if(!is.null(im.2) & !is.null(im.3)){
      img.comb <- EBImage::combine((((im.2[[i]] + bcg2[1]) * bcg2[2]) ^ bcg2[3]),
                                   (((im.1[[i]] + bcg1[1]) * bcg1[2]) ^ bcg1[3]),
                                   (((im.3[[i]] + bcg3[1]) * bcg3[2]) ^ bcg3[3]))

      if (!is.null(masks)){
        img.comb <- EBImage::combine((((im.2[[i]] + bcg2[1]) * bcg2[2]) ^ bcg2[3]) + border.comb,
                                     (((im.1[[i]] + bcg1[1]) * bcg1[2]) ^ bcg1[3]) + border.comb,
                                     (((im.3[[i]] + bcg3[1]) * bcg3[2]) ^ bcg3[3]) + border.comb)
      }
      img.comb <- EBImage::Image(img.comb)
      colorMode(img.comb) <- Color
    }

    DisplayorSave(image=img.comb, filename = images[i], ...)
  }
}



#' DisplayColorImages
#'
#' Display a color image when custom colors have already been defined
#'
#' @param cluster.image list of images in data format where cell ids have been replaced by cluster-dependent colors
#' @param images list of images in data format
#' @param borders TRUE/FALSE should cell borders be displayed?
#' @param bcg brightness / contrast / gamma (default = c(0,1,1))
#' @param ... additional parameters for the label and legend
#' @keywords display color image counts
#' @return nothing is explicitly returned
#' @import EBImage
#' @export

DisplayColorImages <- function(cluster.image, images, borders, bcg=c(0,1,1), ...) {

  if(borders==TRUE)
    cb <- getBorders(images = images)

  for(i in seq_along(cluster.image)){
    colorImg <- EBImage::Image(cluster.image[[i]])
    colorImg <- EBImage::combine(((colorImg + bcg[1]) * bcg[2]) ^ bcg[3])

    if(borders==TRUE){
      borderImg <- EBImage::Image(combine(cb[[i]], cb[[i]], cb[[i]]))
      colorMode(borderImg) <- Color
      colorImg <- combine(colorImg + borderImg)
    }
    DisplayorSave(image=colorImg, filename = names(images[i]), ...)
  }
}



#' DisplayCellsAndCounts
#'
#' Display a channel heatmap or an IMC image (white) together with the borders of one to three cell populations (RGB)
#'
#' @param value.image channel heatmap to display
#' @param cell.sel.1 list of images in data format containing the 1st cell population to display
#' @param cell.sel.2 list of images in data format containing the 2nd cell population to display
#' @param cell.sel.3 list of images in data format containing the 3rd cell population to display
#' @param im.list vector with names of images to be displayed
#' @param cb list of images in data format with cell borders indicated
#' @param bcg brightness / contrast / gamma (e.g. c(0,4,2))
#' @param ... additional parameters for the label and legend
#' @keywords display cell population channel count
#' @return nothing is explicitly returned
#' @import EBImage
#' @export

DisplayCellsAndCounts <- function(value.image, cell.sel.1, cell.sel.2, cell.sel.3, im.list, cb, bcg, ...) {

  for (i in seq_along(im.list)) {
    img.comb <- EBImage::combine(
      ((((value.image[[i]] + bcg[1]) * bcg[2]) ^ bcg[3]) + (cell.sel.2[[i]]) * cb[[i]] ^ 2),
      ((((value.image[[i]] + bcg[1]) * bcg[2]) ^ bcg[3]) + (cell.sel.1[[i]]) * cb[[i]] ^ 2),
      ((((value.image[[i]] + bcg[1]) * bcg[2]) ^ bcg[3]) + (cell.sel.3[[i]]) * cb[[i]] ^ 2)
    )

    img.comb <- EBImage::Image(img.comb)
    colorMode(img.comb) <- Color
    DisplayorSave(image=img.comb, filename = im.list[i], ...)
  }
}



#' DisplayIMCandCells
#'
#' Display one to three IMC images (RGB) together with the borders of one cell population (white)
#'
#' @param value.image.1 heatmaps for the 1st channel to display
#' @param value.image.2 heatmaps for the 2nd channel to display
#' @param value.image.3 heatmaps for the 3rd channel to display
#' @param cell.sel.1 list of images in data format containing the cell population to display
#' @param im.list vector with names of images to be displayed
#' @param cb list of images in data format with cell borders indicated
#' @param bcg1 brightness / contrast / gamma for the 1st image (e.g. = c(0,4,2))
#' @param bcg2 brightness / contrast / gamma for the 2nd image (e.g. = c(0,4,2))
#' @param bcg3 brightness / contrast / gamma for the 3rd image (e.g. = c(0,4,2))
#' @param ... additional parameters for the label and legend
#' @keywords display image borders IMC
#' @return nothing is explicitly returned
#' @import EBImage
#' @export

DisplayIMCandCells <- function(value.image.1, value.image.2, value.image.3, cell.sel.1, im.list, cb, bcg1, bcg2, bcg3, ...) {

  for (i in seq_along(im.list)) {

    if(is.null(value.image.2) & is.null(value.image.3)){
      img.comb <- EBImage::combine(
        ((cell.sel.1[[i]]) * cb[[i]] ^ 2),
        ((((value.image.1[[i]] + bcg1[1]) * bcg1[2]) ^ bcg1[3]) + (cell.sel.1[[i]]) * cb[[i]] ^ 2),
        ((cell.sel.1[[i]]) * cb[[i]] ^ 2)
      )
    }

    if(!is.null(value.image.2) & is.null(value.image.3)){
      img.comb <- EBImage::combine(
        ((((value.image.2[[i]] + bcg2[1]) * bcg2[2]) ^ bcg2[3]) + (cell.sel.1[[i]]) * cb[[i]] ^ 2),
        ((((value.image.1[[i]] + bcg1[1]) * bcg1[2]) ^ bcg1[3]) + (cell.sel.1[[i]]) * cb[[i]] ^ 2),
        ((cell.sel.1[[i]]) * cb[[i]] ^ 2)
      )
    }

    if(!is.null(value.image.2) & !is.null(value.image.3)){
      img.comb <- EBImage::combine(
        ((((value.image.2[[i]] + bcg2[1]) * bcg2[2]) ^ bcg2[3]) + (cell.sel.1[[i]]) * cb[[i]] ^ 2),
        ((((value.image.1[[i]] + bcg1[1]) * bcg1[2]) ^ bcg1[3]) + (cell.sel.1[[i]]) * cb[[i]] ^ 2),
        ((((value.image.3[[i]] + bcg3[1]) * bcg3[2]) ^ bcg3[3]) + (cell.sel.1[[i]]) * cb[[i]] ^ 2)
      )
    }

    img.comb <- EBImage::Image(img.comb)
    colorMode(img.comb) <- Color

    DisplayorSave(image=img.comb, filename = im.list[i], ...)
  }
}



#' DisplayorSave
#'
#' Display or save an image
#'
#' @param ... additional parameters for the legend (see ?legend)
#' @param image image to be displayed
#' @param method 'raster' or 'browser'
#' @param save TRUE/FALSE should the image be saved?
#' @param path path to the folder where images should be saved
#' @param filename optional filename for the images (note: should be modified and used as a suffix or prefix)
#' @param type format in which the image should be saved ("png", "tiff" or "jpeg")
#' @param fun name of the parent function
#' @param label TRUE/FALSE: should the name of the image be displayed on the image (default = TRUE)\cr
#' Only works with the 'raster' method. The label is also not displayed on saved images
#' @param label.txt text to be display on the image label. The same text will be displayed on all images (default = name of the image)
#' @param label.color color of the image label displayed on the image
#' @param label.cex size of the text of the image name label (character expansion)
#' @param label.x x position of the image name label
#' @keywords display save image label
#' @return the images are displayed or save but nothing is explicitly returned
#' @import EBImage
#' @export

DisplayorSave <- function(..., image, method='raster', save=FALSE, path=getwd(), filename="image", type="png", fun="",
                          label=TRUE, label.txt=NULL, label.color="red", label.x=50, label.y=50, label.cex=3){

  if(save == FALSE){
    display(image, method=method)
    if (label != FALSE)
      text(x=label.x, y=label.y, label=ifelse(!is.null(label.txt), label.txt, filename), col=label.color, cex=label.cex)

  } else {
    if (!dir.exists(path)) dir.create(path)
    fn <- file.path(path, paste0(fun, "_", filename, ".", type))
    writeImage(image, files=fn, type=type)
  }
}



#' makeLegend
#'
#' Create a legend as a separate image
#'
#' @param ... additional parameters for the legend (see ?legend)
#' @param leg.names vector containing the names of the elements displayed in the legend
#' @param leg.fill vector containing the colors to be displayed in the legend
#' @param leg.border vector containing the color of the legend borders
#' @param save TRUE/FALSE should the legend be saved?
#' @param path path to the folder where legend should be saved
#' @param filename optional filename for the legend (note: should be modified and used as a suffix or prefix)
#' @param fun name of the parent function
#' @param nb.items number of items to be displayed in the legend
#' @keywords legend
#' @return the legend is displayed or saved, but nothing is explicitly returned
#' @export

makeLegend <- function(..., leg.names = c("1", "2", "3"), leg.fill = c('green', 'red', 'blue'), leg.border = c('black'),
                       save=FALSE, path=getwd(), filename="legend", fun="", nb.items=3){

  dots <- list(...)
  dots <- dots[unlist(names(dots) %in% names(formals(legend)))]

  if (!'x' %in% names(dots))          dots$x <- 1
  if (!'y' %in% names(dots))          dots$y <- 1
  if (!'fill' %in% names(dots))       dots$fill   <- leg.fill[1:nb.items]       else dots$fill   <- dots$fill[1:nb.items]
  if (!'border' %in% names(dots))     dots$border <- rep(leg.border, nb.items)  else dots$border <- dots$border[1:nb.items]
  if (!'legend' %in% names(dots))     dots$legend <- leg.names[1:nb.items]      else dots$legend <- dots$legend[1:nb.items]
  if (!'cex' %in% names(dots))        dots$cex <- 2
  if (!'xjust' %in% names(dots))      dots$xjust = 0.5
  if (!'yjust' %in% names(dots))      dots$yjust = 0.5
  if (!'title' %in% names(dots))      dots$title = fun

  if (save == FALSE) {
    plot(1, type="n", axes=FALSE, xlab="", ylab="")
    do.call(legend, args=dots)

  } else {
    if (!dir.exists(path)) dir.create(path)
    fn <- file.path(path, paste0(fun, "_", filename, ".png"))
    png(filename=fn)
    plot(1, type="n", axes=FALSE, xlab="", ylab="")
    do.call(legend, args=dots)
  }
}







