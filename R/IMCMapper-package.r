#' IMCMapper.
#'
#' @name IMCMapper
#' @docType package
#'
#' Functions for displaying IMC images and data in R
#'
#' - CellMap: Displays 1 to 3 cell populations
#'
#' - ChannelHM: Heatmap 1 to 3 channels
#'
#' - CellsAndCounts: Heatmap 1 channel + highlight the borders of 1 to 3 cell populations
#'
#' - ClusterDisplay: Display cell clusters (and unique continuous values)
#'
#' - IMCDisplay: Display 1 to 3 raw IMC images
#'
#' - CellsandIMC: Display 1 IMC image + highlight the borders of 1 to 3 cell populations
#'
#' - IMCandCells: Display 1 to 3 IMC images + highlight the borders of 1 cell population
#'
#'
#'
#' GENERAL REQUIREMENTS
#'
#' R packages
#' - "data.table"
#' - "EBImage"
#'
#' CellProfiler data
#' - Cell masks (16-bit tiff files).
#' - When exporting CSV from CellProfiler, in "Select measurements" > "Image", you have to select FileName, PathName and Metadata.
#' - This will generate an "Image.csv" file, which is required to run this script.
#'
#' A data set in the data.table format
#' - The cell objects must be the same as those used to create the masks in CellProfiler.
#' - Must contain an "image" column where the values correspond to the image names in the "Metadata" column of your "Image.csv" file.
#' - Must contain a "cell id" column in the format "###_nnn", where ### is the name of the image in the "image" column and nnn is the number of the cell.
#'
#'
#'
NULL
