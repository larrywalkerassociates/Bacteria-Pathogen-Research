# pathogen bacteria analaysis
scripts developed to support investigation of bacteria and pathogens

## Getting started

Configure environmental variables. The easiest way to do this is from within RStudio with the `{usethis}` package. Enter `usethis::edit_r_environ(scope = "project")` to edit the project-level `.Renviron` file, then edit the file and save with:

```R
DATA_PATH = "<local path to Dropbox data>"
```

The `DATA_PATH` is where project data dependencies are stored on a synced Dropbox folder. On my computer it's DATA_PATH = "C:/Users/RebeccaSmith/LWA Dropbox/Rebecca Smith/5 - Projects/Bacteria Study (ULAR and SMB)/14_Rec Zone Characterization (2025)/03_Data Analysis/Regression" 

To determine if the DATA_PATH variable has been read, enter the following into the console: 
  readRenviron(".Renviron")
  Sys.getenv("DATA_PATH")
  
  and add the following to the script below:
  #import data:
  data_path <- Sys.getenv("DATA_PATH")

The `EPSG` is the projection used in this study that all spatial data are standardized to.  


Restart `R` for changes to take effect.

If you notice issues with the location of your libraries (ggplot2, tidyverse, etc.), try running .libPaths() and see where they are being stored
  on my system they are stored locally in the following: 
      [1] "C:/Users/RebeccaSmith/AppData/Local/R/win-library/4.2" "C:/Program Files/R/R-4.2.2/library"
  you can specify the path to your libraries by writing `.libPaths("file-path-here')`

***

Last updated by *Rebecca Smith* on 2026-05-20