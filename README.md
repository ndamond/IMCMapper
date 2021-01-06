
**IMCMapper** is now deprecated, please use **cytomapper** instead
- Source code: [Github](https://github.com/BodenmillerGroup/cytomapper/)
- R package: [Bioconductor](https://www.bioconductor.org/packages/release/bioc/html/cytomapper.html)
- Publication: [Bioinformatics](https://doi.org/10.1093/bioinformatics/btaa1061)


# MAIN FUNCTIONS

**- CellMap:** Displays 1 to 3 cell populations

**- ChannelHM:** Heatmap 1 to 3 channels

**- CellsAndCounts:** Heatmap 1 channel + highlight the borders of 1 to 3 cell populations

**- ClusterDisplay:** Display cell clusters (and unique continuous values)

**- IMCDisplay:** Display 1 to 3 raw IMC images

**- CellsandIMC:** Display 1 IMC image + highlight the borders of 1 to 3 cell populations

**- IMCandCells:** Display 1 to 3 IMC images + highlight the borders of 1 cell population  
  
***

# GENERAL REQUIREMENTS

### R packages
- "data.table"  
- "EBImage"  

### CellProfiler data
- Cell masks (16-bit tiff files).
- When exporting CSV from CellProfiler, in "Select measurements" > "Image", you have to select FileName, PathName and Metadata.
- This will generate an "Image.csv" file, which is required to run this script.

### A data set in the data.table format
- The cell objects must be the same as those used to create the masks in CellProfiler.
- Must contain an "image" column where the values correspond to the image names in the "Metadata" column of your "Image.csv" file.
- Must contain a "cell id" column in the format "###_nnn", where ### is the name of the image in the "image" column and nnn is the number of the cell.

