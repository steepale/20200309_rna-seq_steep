Examination of Rat RNA-Seq Data: Examine Batch Effects
================
Alec Steep
20200309

# Goals of Analysis

  - 
## Setup the Environment

``` r
################################################################################
##### Resources and Dependencies ###############################################
################################################################################

# Set the working directory
WD <- '/Volumes/Frishman_4TB/motrpac/20200309_rna-seq_steep'
setwd(WD)

# Load the dependencies
#source("https://bioconductor.org/biocLite.R")
#BiocManager::install("ff")
#install.packages("tidyverse")

# Load dependencies
pacs...man <- c("tidyverse","data.table","R.utils")
lapply(pacs...man, FUN = function(X) {
        do.call("library", list(X)) 
})
```

    ## [[1]]
    ##  [1] "forcats"   "stringr"   "dplyr"     "purrr"     "readr"     "tidyr"    
    ##  [7] "tibble"    "ggplot2"   "tidyverse" "stats"     "graphics"  "grDevices"
    ## [13] "utils"     "datasets"  "methods"   "base"     
    ## 
    ## [[2]]
    ##  [1] "data.table" "forcats"    "stringr"    "dplyr"      "purrr"     
    ##  [6] "readr"      "tidyr"      "tibble"     "ggplot2"    "tidyverse" 
    ## [11] "stats"      "graphics"   "grDevices"  "utils"      "datasets"  
    ## [16] "methods"    "base"      
    ## 
    ## [[3]]
    ##  [1] "R.utils"     "R.oo"        "R.methodsS3" "data.table"  "forcats"    
    ##  [6] "stringr"     "dplyr"       "purrr"       "readr"       "tidyr"      
    ## [11] "tibble"      "ggplot2"     "tidyverse"   "stats"       "graphics"   
    ## [16] "grDevices"   "utils"       "datasets"    "methods"     "base"

``` r
############################################################
##### Functions ############################################
############################################################

# Make the 'not in' operator
################################################################################
'%!in%' <- function(x,y) {
        !('%in%'(x,y))
}
################################################################################

# Capture the Date and AUthor
################################################################################
date <- format.Date( Sys.Date(), '%Y%m%d' )
auth <- "steep"
################################################################################

## explicit gc, then execute `expr` `n` times w/o explicit gc, return timings
################################################################################
benchmark <- function(n = 1, expr, envir = parent.frame()) {
        expr <- substitute(expr)
        gc()
        map(seq_len(n), ~ system.time(eval(expr, envir), gcFirst = FALSE))
}
################################################################################

# Function to speed up making rows into lists for interation with lapply
################################################################################
f_pmap_aslist <- function(df) {
        purrr::pmap(as.list(df), list)
}
################################################################################
```

## Load & Clean Data

##### Data Files to Load:

  - GATK 600K-Targetted Raw SNPs: unfiltered calls from 22 F1 control
    samples as g vcf (collective
genotypes)

<!-- end list -->

``` r
################################################################################
#####     Load & Clean Data      ###############################################
################################################################################
```
