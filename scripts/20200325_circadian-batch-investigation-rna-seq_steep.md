PASS1A Rat Liver (Stanford Batch 1) RNA-Seq: Circadian vs. Exercise:
What Drives Variance?
================
Alec Steep and Jiayu Zhang
20200331

## Goals of Analysis

  - TODO: Continue Step-by-step goals of analysis
  - Perform a PCA of Liver
      - If Liver clusters by sex, investigate clustering by autosomal
        and sex genes
          - Histogram of p values by autosomal genes and sex genes (BOTH
            X and Y chromosomes\!)
          - Match with volcano plots
      - Either correct for sex batch, remove sex genes, or subset one
        sex.
  - Perform another PCA with exercise groups, do they cluster together
    by:
      - Time of day?
      - Exercise cohort?
      - Other metadata? – EDA (dunno yet)
  - Build a correlation matrix – EDA (dunno yet)
  - Subset genes into 2 groups:
      - Genes associated with circadian rhythm in Rat Liver (We don’t
        have this set, we use Mouse Liver)
      - Genes associated with exercise effects (TODO: Gene sets need to
        be collected).
  - Investiagte batch
      - Examine overlap of 2 gene sets
      - Rank genes by p-value (dunno how to do this yet, unless by
        differential expression)
          - TODO: Determien ranking strategy
      - Statistical test involving rank – dunno: wilcoxen rank sum test?
          - TODO: Determine proper test
      - Cluster: PCAs – EDA
  - To be continued
…

## Setup the Environment

``` r
################################################################################
##### Resources and Dependencies ###############################################
################################################################################

# Set the working directory
WD <- '/Volumes/Frishman_4TB/motrpac/20200309_rna-seq_steep'
#setwd(WD)

# Load the dependencies
#source("https://bioconductor.org/biocLite.R")
#BiocManager::install("qvalue")
#install.packages("tidyverse")

# Load dependencies
pacs...man <- c("tidyverse","GenomicRanges", "DESeq2","devtools","rafalib","GO.db","vsn","hexbin","ggplot2", "GenomicFeatures","Biostrings","BSgenome","AnnotationHub","plyr","dplyr", "org.Rn.eg.db","pheatmap","sva","formula.tools","pathview","biomaRt", "PROPER","SeqGSEA",'purrr','BioInstaller','RColorBrewer','lubridate', "hms","ggpubr", "ggrepel","genefilter","qvalue")
lapply(pacs...man, FUN = function(X) {
        do.call("library", list(X)) })
```

    ## [[1]]
    ##  [1] "forcats"   "stringr"   "dplyr"     "purrr"     "readr"     "tidyr"    
    ##  [7] "tibble"    "ggplot2"   "tidyverse" "stats"     "graphics"  "grDevices"
    ## [13] "utils"     "datasets"  "methods"   "base"     
    ## 
    ## [[2]]
    ##  [1] "GenomicRanges" "GenomeInfoDb"  "IRanges"       "S4Vectors"    
    ##  [5] "BiocGenerics"  "parallel"      "stats4"        "forcats"      
    ##  [9] "stringr"       "dplyr"         "purrr"         "readr"        
    ## [13] "tidyr"         "tibble"        "ggplot2"       "tidyverse"    
    ## [17] "stats"         "graphics"      "grDevices"     "utils"        
    ## [21] "datasets"      "methods"       "base"         
    ## 
    ## [[3]]
    ##  [1] "DESeq2"               "SummarizedExperiment" "DelayedArray"        
    ##  [4] "BiocParallel"         "matrixStats"          "Biobase"             
    ##  [7] "GenomicRanges"        "GenomeInfoDb"         "IRanges"             
    ## [10] "S4Vectors"            "BiocGenerics"         "parallel"            
    ## [13] "stats4"               "forcats"              "stringr"             
    ## [16] "dplyr"                "purrr"                "readr"               
    ## [19] "tidyr"                "tibble"               "ggplot2"             
    ## [22] "tidyverse"            "stats"                "graphics"            
    ## [25] "grDevices"            "utils"                "datasets"            
    ## [28] "methods"              "base"                
    ## 
    ## [[4]]
    ##  [1] "devtools"             "usethis"              "DESeq2"              
    ##  [4] "SummarizedExperiment" "DelayedArray"         "BiocParallel"        
    ##  [7] "matrixStats"          "Biobase"              "GenomicRanges"       
    ## [10] "GenomeInfoDb"         "IRanges"              "S4Vectors"           
    ## [13] "BiocGenerics"         "parallel"             "stats4"              
    ## [16] "forcats"              "stringr"              "dplyr"               
    ## [19] "purrr"                "readr"                "tidyr"               
    ## [22] "tibble"               "ggplot2"              "tidyverse"           
    ## [25] "stats"                "graphics"             "grDevices"           
    ## [28] "utils"                "datasets"             "methods"             
    ## [31] "base"                
    ## 
    ## [[5]]
    ##  [1] "rafalib"              "devtools"             "usethis"             
    ##  [4] "DESeq2"               "SummarizedExperiment" "DelayedArray"        
    ##  [7] "BiocParallel"         "matrixStats"          "Biobase"             
    ## [10] "GenomicRanges"        "GenomeInfoDb"         "IRanges"             
    ## [13] "S4Vectors"            "BiocGenerics"         "parallel"            
    ## [16] "stats4"               "forcats"              "stringr"             
    ## [19] "dplyr"                "purrr"                "readr"               
    ## [22] "tidyr"                "tibble"               "ggplot2"             
    ## [25] "tidyverse"            "stats"                "graphics"            
    ## [28] "grDevices"            "utils"                "datasets"            
    ## [31] "methods"              "base"                
    ## 
    ## [[6]]
    ##  [1] "GO.db"                "AnnotationDbi"        "rafalib"             
    ##  [4] "devtools"             "usethis"              "DESeq2"              
    ##  [7] "SummarizedExperiment" "DelayedArray"         "BiocParallel"        
    ## [10] "matrixStats"          "Biobase"              "GenomicRanges"       
    ## [13] "GenomeInfoDb"         "IRanges"              "S4Vectors"           
    ## [16] "BiocGenerics"         "parallel"             "stats4"              
    ## [19] "forcats"              "stringr"              "dplyr"               
    ## [22] "purrr"                "readr"                "tidyr"               
    ## [25] "tibble"               "ggplot2"              "tidyverse"           
    ## [28] "stats"                "graphics"             "grDevices"           
    ## [31] "utils"                "datasets"             "methods"             
    ## [34] "base"                
    ## 
    ## [[7]]
    ##  [1] "vsn"                  "GO.db"                "AnnotationDbi"       
    ##  [4] "rafalib"              "devtools"             "usethis"             
    ##  [7] "DESeq2"               "SummarizedExperiment" "DelayedArray"        
    ## [10] "BiocParallel"         "matrixStats"          "Biobase"             
    ## [13] "GenomicRanges"        "GenomeInfoDb"         "IRanges"             
    ## [16] "S4Vectors"            "BiocGenerics"         "parallel"            
    ## [19] "stats4"               "forcats"              "stringr"             
    ## [22] "dplyr"                "purrr"                "readr"               
    ## [25] "tidyr"                "tibble"               "ggplot2"             
    ## [28] "tidyverse"            "stats"                "graphics"            
    ## [31] "grDevices"            "utils"                "datasets"            
    ## [34] "methods"              "base"                
    ## 
    ## [[8]]
    ##  [1] "hexbin"               "vsn"                  "GO.db"               
    ##  [4] "AnnotationDbi"        "rafalib"              "devtools"            
    ##  [7] "usethis"              "DESeq2"               "SummarizedExperiment"
    ## [10] "DelayedArray"         "BiocParallel"         "matrixStats"         
    ## [13] "Biobase"              "GenomicRanges"        "GenomeInfoDb"        
    ## [16] "IRanges"              "S4Vectors"            "BiocGenerics"        
    ## [19] "parallel"             "stats4"               "forcats"             
    ## [22] "stringr"              "dplyr"                "purrr"               
    ## [25] "readr"                "tidyr"                "tibble"              
    ## [28] "ggplot2"              "tidyverse"            "stats"               
    ## [31] "graphics"             "grDevices"            "utils"               
    ## [34] "datasets"             "methods"              "base"                
    ## 
    ## [[9]]
    ##  [1] "hexbin"               "vsn"                  "GO.db"               
    ##  [4] "AnnotationDbi"        "rafalib"              "devtools"            
    ##  [7] "usethis"              "DESeq2"               "SummarizedExperiment"
    ## [10] "DelayedArray"         "BiocParallel"         "matrixStats"         
    ## [13] "Biobase"              "GenomicRanges"        "GenomeInfoDb"        
    ## [16] "IRanges"              "S4Vectors"            "BiocGenerics"        
    ## [19] "parallel"             "stats4"               "forcats"             
    ## [22] "stringr"              "dplyr"                "purrr"               
    ## [25] "readr"                "tidyr"                "tibble"              
    ## [28] "ggplot2"              "tidyverse"            "stats"               
    ## [31] "graphics"             "grDevices"            "utils"               
    ## [34] "datasets"             "methods"              "base"                
    ## 
    ## [[10]]
    ##  [1] "GenomicFeatures"      "hexbin"               "vsn"                 
    ##  [4] "GO.db"                "AnnotationDbi"        "rafalib"             
    ##  [7] "devtools"             "usethis"              "DESeq2"              
    ## [10] "SummarizedExperiment" "DelayedArray"         "BiocParallel"        
    ## [13] "matrixStats"          "Biobase"              "GenomicRanges"       
    ## [16] "GenomeInfoDb"         "IRanges"              "S4Vectors"           
    ## [19] "BiocGenerics"         "parallel"             "stats4"              
    ## [22] "forcats"              "stringr"              "dplyr"               
    ## [25] "purrr"                "readr"                "tidyr"               
    ## [28] "tibble"               "ggplot2"              "tidyverse"           
    ## [31] "stats"                "graphics"             "grDevices"           
    ## [34] "utils"                "datasets"             "methods"             
    ## [37] "base"                
    ## 
    ## [[11]]
    ##  [1] "Biostrings"           "XVector"              "GenomicFeatures"     
    ##  [4] "hexbin"               "vsn"                  "GO.db"               
    ##  [7] "AnnotationDbi"        "rafalib"              "devtools"            
    ## [10] "usethis"              "DESeq2"               "SummarizedExperiment"
    ## [13] "DelayedArray"         "BiocParallel"         "matrixStats"         
    ## [16] "Biobase"              "GenomicRanges"        "GenomeInfoDb"        
    ## [19] "IRanges"              "S4Vectors"            "BiocGenerics"        
    ## [22] "parallel"             "stats4"               "forcats"             
    ## [25] "stringr"              "dplyr"                "purrr"               
    ## [28] "readr"                "tidyr"                "tibble"              
    ## [31] "ggplot2"              "tidyverse"            "stats"               
    ## [34] "graphics"             "grDevices"            "utils"               
    ## [37] "datasets"             "methods"              "base"                
    ## 
    ## [[12]]
    ##  [1] "BSgenome"             "rtracklayer"          "Biostrings"          
    ##  [4] "XVector"              "GenomicFeatures"      "hexbin"              
    ##  [7] "vsn"                  "GO.db"                "AnnotationDbi"       
    ## [10] "rafalib"              "devtools"             "usethis"             
    ## [13] "DESeq2"               "SummarizedExperiment" "DelayedArray"        
    ## [16] "BiocParallel"         "matrixStats"          "Biobase"             
    ## [19] "GenomicRanges"        "GenomeInfoDb"         "IRanges"             
    ## [22] "S4Vectors"            "BiocGenerics"         "parallel"            
    ## [25] "stats4"               "forcats"              "stringr"             
    ## [28] "dplyr"                "purrr"                "readr"               
    ## [31] "tidyr"                "tibble"               "ggplot2"             
    ## [34] "tidyverse"            "stats"                "graphics"            
    ## [37] "grDevices"            "utils"                "datasets"            
    ## [40] "methods"              "base"                
    ## 
    ## [[13]]
    ##  [1] "AnnotationHub"        "BiocFileCache"        "dbplyr"              
    ##  [4] "BSgenome"             "rtracklayer"          "Biostrings"          
    ##  [7] "XVector"              "GenomicFeatures"      "hexbin"              
    ## [10] "vsn"                  "GO.db"                "AnnotationDbi"       
    ## [13] "rafalib"              "devtools"             "usethis"             
    ## [16] "DESeq2"               "SummarizedExperiment" "DelayedArray"        
    ## [19] "BiocParallel"         "matrixStats"          "Biobase"             
    ## [22] "GenomicRanges"        "GenomeInfoDb"         "IRanges"             
    ## [25] "S4Vectors"            "BiocGenerics"         "parallel"            
    ## [28] "stats4"               "forcats"              "stringr"             
    ## [31] "dplyr"                "purrr"                "readr"               
    ## [34] "tidyr"                "tibble"               "ggplot2"             
    ## [37] "tidyverse"            "stats"                "graphics"            
    ## [40] "grDevices"            "utils"                "datasets"            
    ## [43] "methods"              "base"                
    ## 
    ## [[14]]
    ##  [1] "plyr"                 "AnnotationHub"        "BiocFileCache"       
    ##  [4] "dbplyr"               "BSgenome"             "rtracklayer"         
    ##  [7] "Biostrings"           "XVector"              "GenomicFeatures"     
    ## [10] "hexbin"               "vsn"                  "GO.db"               
    ## [13] "AnnotationDbi"        "rafalib"              "devtools"            
    ## [16] "usethis"              "DESeq2"               "SummarizedExperiment"
    ## [19] "DelayedArray"         "BiocParallel"         "matrixStats"         
    ## [22] "Biobase"              "GenomicRanges"        "GenomeInfoDb"        
    ## [25] "IRanges"              "S4Vectors"            "BiocGenerics"        
    ## [28] "parallel"             "stats4"               "forcats"             
    ## [31] "stringr"              "dplyr"                "purrr"               
    ## [34] "readr"                "tidyr"                "tibble"              
    ## [37] "ggplot2"              "tidyverse"            "stats"               
    ## [40] "graphics"             "grDevices"            "utils"               
    ## [43] "datasets"             "methods"              "base"                
    ## 
    ## [[15]]
    ##  [1] "plyr"                 "AnnotationHub"        "BiocFileCache"       
    ##  [4] "dbplyr"               "BSgenome"             "rtracklayer"         
    ##  [7] "Biostrings"           "XVector"              "GenomicFeatures"     
    ## [10] "hexbin"               "vsn"                  "GO.db"               
    ## [13] "AnnotationDbi"        "rafalib"              "devtools"            
    ## [16] "usethis"              "DESeq2"               "SummarizedExperiment"
    ## [19] "DelayedArray"         "BiocParallel"         "matrixStats"         
    ## [22] "Biobase"              "GenomicRanges"        "GenomeInfoDb"        
    ## [25] "IRanges"              "S4Vectors"            "BiocGenerics"        
    ## [28] "parallel"             "stats4"               "forcats"             
    ## [31] "stringr"              "dplyr"                "purrr"               
    ## [34] "readr"                "tidyr"                "tibble"              
    ## [37] "ggplot2"              "tidyverse"            "stats"               
    ## [40] "graphics"             "grDevices"            "utils"               
    ## [43] "datasets"             "methods"              "base"                
    ## 
    ## [[16]]
    ##  [1] "org.Rn.eg.db"         "plyr"                 "AnnotationHub"       
    ##  [4] "BiocFileCache"        "dbplyr"               "BSgenome"            
    ##  [7] "rtracklayer"          "Biostrings"           "XVector"             
    ## [10] "GenomicFeatures"      "hexbin"               "vsn"                 
    ## [13] "GO.db"                "AnnotationDbi"        "rafalib"             
    ## [16] "devtools"             "usethis"              "DESeq2"              
    ## [19] "SummarizedExperiment" "DelayedArray"         "BiocParallel"        
    ## [22] "matrixStats"          "Biobase"              "GenomicRanges"       
    ## [25] "GenomeInfoDb"         "IRanges"              "S4Vectors"           
    ## [28] "BiocGenerics"         "parallel"             "stats4"              
    ## [31] "forcats"              "stringr"              "dplyr"               
    ## [34] "purrr"                "readr"                "tidyr"               
    ## [37] "tibble"               "ggplot2"              "tidyverse"           
    ## [40] "stats"                "graphics"             "grDevices"           
    ## [43] "utils"                "datasets"             "methods"             
    ## [46] "base"                
    ## 
    ## [[17]]
    ##  [1] "pheatmap"             "org.Rn.eg.db"         "plyr"                
    ##  [4] "AnnotationHub"        "BiocFileCache"        "dbplyr"              
    ##  [7] "BSgenome"             "rtracklayer"          "Biostrings"          
    ## [10] "XVector"              "GenomicFeatures"      "hexbin"              
    ## [13] "vsn"                  "GO.db"                "AnnotationDbi"       
    ## [16] "rafalib"              "devtools"             "usethis"             
    ## [19] "DESeq2"               "SummarizedExperiment" "DelayedArray"        
    ## [22] "BiocParallel"         "matrixStats"          "Biobase"             
    ## [25] "GenomicRanges"        "GenomeInfoDb"         "IRanges"             
    ## [28] "S4Vectors"            "BiocGenerics"         "parallel"            
    ## [31] "stats4"               "forcats"              "stringr"             
    ## [34] "dplyr"                "purrr"                "readr"               
    ## [37] "tidyr"                "tibble"               "ggplot2"             
    ## [40] "tidyverse"            "stats"                "graphics"            
    ## [43] "grDevices"            "utils"                "datasets"            
    ## [46] "methods"              "base"                
    ## 
    ## [[18]]
    ##  [1] "sva"                  "genefilter"           "mgcv"                
    ##  [4] "nlme"                 "pheatmap"             "org.Rn.eg.db"        
    ##  [7] "plyr"                 "AnnotationHub"        "BiocFileCache"       
    ## [10] "dbplyr"               "BSgenome"             "rtracklayer"         
    ## [13] "Biostrings"           "XVector"              "GenomicFeatures"     
    ## [16] "hexbin"               "vsn"                  "GO.db"               
    ## [19] "AnnotationDbi"        "rafalib"              "devtools"            
    ## [22] "usethis"              "DESeq2"               "SummarizedExperiment"
    ## [25] "DelayedArray"         "BiocParallel"         "matrixStats"         
    ## [28] "Biobase"              "GenomicRanges"        "GenomeInfoDb"        
    ## [31] "IRanges"              "S4Vectors"            "BiocGenerics"        
    ## [34] "parallel"             "stats4"               "forcats"             
    ## [37] "stringr"              "dplyr"                "purrr"               
    ## [40] "readr"                "tidyr"                "tibble"              
    ## [43] "ggplot2"              "tidyverse"            "stats"               
    ## [46] "graphics"             "grDevices"            "utils"               
    ## [49] "datasets"             "methods"              "base"                
    ## 
    ## [[19]]
    ##  [1] "formula.tools"        "sva"                  "genefilter"          
    ##  [4] "mgcv"                 "nlme"                 "pheatmap"            
    ##  [7] "org.Rn.eg.db"         "plyr"                 "AnnotationHub"       
    ## [10] "BiocFileCache"        "dbplyr"               "BSgenome"            
    ## [13] "rtracklayer"          "Biostrings"           "XVector"             
    ## [16] "GenomicFeatures"      "hexbin"               "vsn"                 
    ## [19] "GO.db"                "AnnotationDbi"        "rafalib"             
    ## [22] "devtools"             "usethis"              "DESeq2"              
    ## [25] "SummarizedExperiment" "DelayedArray"         "BiocParallel"        
    ## [28] "matrixStats"          "Biobase"              "GenomicRanges"       
    ## [31] "GenomeInfoDb"         "IRanges"              "S4Vectors"           
    ## [34] "BiocGenerics"         "parallel"             "stats4"              
    ## [37] "forcats"              "stringr"              "dplyr"               
    ## [40] "purrr"                "readr"                "tidyr"               
    ## [43] "tibble"               "ggplot2"              "tidyverse"           
    ## [46] "stats"                "graphics"             "grDevices"           
    ## [49] "utils"                "datasets"             "methods"             
    ## [52] "base"                
    ## 
    ## [[20]]
    ##  [1] "pathview"             "org.Hs.eg.db"         "formula.tools"       
    ##  [4] "sva"                  "genefilter"           "mgcv"                
    ##  [7] "nlme"                 "pheatmap"             "org.Rn.eg.db"        
    ## [10] "plyr"                 "AnnotationHub"        "BiocFileCache"       
    ## [13] "dbplyr"               "BSgenome"             "rtracklayer"         
    ## [16] "Biostrings"           "XVector"              "GenomicFeatures"     
    ## [19] "hexbin"               "vsn"                  "GO.db"               
    ## [22] "AnnotationDbi"        "rafalib"              "devtools"            
    ## [25] "usethis"              "DESeq2"               "SummarizedExperiment"
    ## [28] "DelayedArray"         "BiocParallel"         "matrixStats"         
    ## [31] "Biobase"              "GenomicRanges"        "GenomeInfoDb"        
    ## [34] "IRanges"              "S4Vectors"            "BiocGenerics"        
    ## [37] "parallel"             "stats4"               "forcats"             
    ## [40] "stringr"              "dplyr"                "purrr"               
    ## [43] "readr"                "tidyr"                "tibble"              
    ## [46] "ggplot2"              "tidyverse"            "stats"               
    ## [49] "graphics"             "grDevices"            "utils"               
    ## [52] "datasets"             "methods"              "base"                
    ## 
    ## [[21]]
    ##  [1] "biomaRt"              "pathview"             "org.Hs.eg.db"        
    ##  [4] "formula.tools"        "sva"                  "genefilter"          
    ##  [7] "mgcv"                 "nlme"                 "pheatmap"            
    ## [10] "org.Rn.eg.db"         "plyr"                 "AnnotationHub"       
    ## [13] "BiocFileCache"        "dbplyr"               "BSgenome"            
    ## [16] "rtracklayer"          "Biostrings"           "XVector"             
    ## [19] "GenomicFeatures"      "hexbin"               "vsn"                 
    ## [22] "GO.db"                "AnnotationDbi"        "rafalib"             
    ## [25] "devtools"             "usethis"              "DESeq2"              
    ## [28] "SummarizedExperiment" "DelayedArray"         "BiocParallel"        
    ## [31] "matrixStats"          "Biobase"              "GenomicRanges"       
    ## [34] "GenomeInfoDb"         "IRanges"              "S4Vectors"           
    ## [37] "BiocGenerics"         "parallel"             "stats4"              
    ## [40] "forcats"              "stringr"              "dplyr"               
    ## [43] "purrr"                "readr"                "tidyr"               
    ## [46] "tibble"               "ggplot2"              "tidyverse"           
    ## [49] "stats"                "graphics"             "grDevices"           
    ## [52] "utils"                "datasets"             "methods"             
    ## [55] "base"                
    ## 
    ## [[22]]
    ##  [1] "PROPER"               "biomaRt"              "pathview"            
    ##  [4] "org.Hs.eg.db"         "formula.tools"        "sva"                 
    ##  [7] "genefilter"           "mgcv"                 "nlme"                
    ## [10] "pheatmap"             "org.Rn.eg.db"         "plyr"                
    ## [13] "AnnotationHub"        "BiocFileCache"        "dbplyr"              
    ## [16] "BSgenome"             "rtracklayer"          "Biostrings"          
    ## [19] "XVector"              "GenomicFeatures"      "hexbin"              
    ## [22] "vsn"                  "GO.db"                "AnnotationDbi"       
    ## [25] "rafalib"              "devtools"             "usethis"             
    ## [28] "DESeq2"               "SummarizedExperiment" "DelayedArray"        
    ## [31] "BiocParallel"         "matrixStats"          "Biobase"             
    ## [34] "GenomicRanges"        "GenomeInfoDb"         "IRanges"             
    ## [37] "S4Vectors"            "BiocGenerics"         "parallel"            
    ## [40] "stats4"               "forcats"              "stringr"             
    ## [43] "dplyr"                "purrr"                "readr"               
    ## [46] "tidyr"                "tibble"               "ggplot2"             
    ## [49] "tidyverse"            "stats"                "graphics"            
    ## [52] "grDevices"            "utils"                "datasets"            
    ## [55] "methods"              "base"                
    ## 
    ## [[23]]
    ##  [1] "SeqGSEA"              "DESeq"                "lattice"             
    ##  [4] "locfit"               "doParallel"           "iterators"           
    ##  [7] "foreach"              "PROPER"               "biomaRt"             
    ## [10] "pathview"             "org.Hs.eg.db"         "formula.tools"       
    ## [13] "sva"                  "genefilter"           "mgcv"                
    ## [16] "nlme"                 "pheatmap"             "org.Rn.eg.db"        
    ## [19] "plyr"                 "AnnotationHub"        "BiocFileCache"       
    ## [22] "dbplyr"               "BSgenome"             "rtracklayer"         
    ## [25] "Biostrings"           "XVector"              "GenomicFeatures"     
    ## [28] "hexbin"               "vsn"                  "GO.db"               
    ## [31] "AnnotationDbi"        "rafalib"              "devtools"            
    ## [34] "usethis"              "DESeq2"               "SummarizedExperiment"
    ## [37] "DelayedArray"         "BiocParallel"         "matrixStats"         
    ## [40] "Biobase"              "GenomicRanges"        "GenomeInfoDb"        
    ## [43] "IRanges"              "S4Vectors"            "BiocGenerics"        
    ## [46] "parallel"             "stats4"               "forcats"             
    ## [49] "stringr"              "dplyr"                "purrr"               
    ## [52] "readr"                "tidyr"                "tibble"              
    ## [55] "ggplot2"              "tidyverse"            "stats"               
    ## [58] "graphics"             "grDevices"            "utils"               
    ## [61] "datasets"             "methods"              "base"                
    ## 
    ## [[24]]
    ##  [1] "SeqGSEA"              "DESeq"                "lattice"             
    ##  [4] "locfit"               "doParallel"           "iterators"           
    ##  [7] "foreach"              "PROPER"               "biomaRt"             
    ## [10] "pathview"             "org.Hs.eg.db"         "formula.tools"       
    ## [13] "sva"                  "genefilter"           "mgcv"                
    ## [16] "nlme"                 "pheatmap"             "org.Rn.eg.db"        
    ## [19] "plyr"                 "AnnotationHub"        "BiocFileCache"       
    ## [22] "dbplyr"               "BSgenome"             "rtracklayer"         
    ## [25] "Biostrings"           "XVector"              "GenomicFeatures"     
    ## [28] "hexbin"               "vsn"                  "GO.db"               
    ## [31] "AnnotationDbi"        "rafalib"              "devtools"            
    ## [34] "usethis"              "DESeq2"               "SummarizedExperiment"
    ## [37] "DelayedArray"         "BiocParallel"         "matrixStats"         
    ## [40] "Biobase"              "GenomicRanges"        "GenomeInfoDb"        
    ## [43] "IRanges"              "S4Vectors"            "BiocGenerics"        
    ## [46] "parallel"             "stats4"               "forcats"             
    ## [49] "stringr"              "dplyr"                "purrr"               
    ## [52] "readr"                "tidyr"                "tibble"              
    ## [55] "ggplot2"              "tidyverse"            "stats"               
    ## [58] "graphics"             "grDevices"            "utils"               
    ## [61] "datasets"             "methods"              "base"                
    ## 
    ## [[25]]
    ##  [1] "BioInstaller"         "SeqGSEA"              "DESeq"               
    ##  [4] "lattice"              "locfit"               "doParallel"          
    ##  [7] "iterators"            "foreach"              "PROPER"              
    ## [10] "biomaRt"              "pathview"             "org.Hs.eg.db"        
    ## [13] "formula.tools"        "sva"                  "genefilter"          
    ## [16] "mgcv"                 "nlme"                 "pheatmap"            
    ## [19] "org.Rn.eg.db"         "plyr"                 "AnnotationHub"       
    ## [22] "BiocFileCache"        "dbplyr"               "BSgenome"            
    ## [25] "rtracklayer"          "Biostrings"           "XVector"             
    ## [28] "GenomicFeatures"      "hexbin"               "vsn"                 
    ## [31] "GO.db"                "AnnotationDbi"        "rafalib"             
    ## [34] "devtools"             "usethis"              "DESeq2"              
    ## [37] "SummarizedExperiment" "DelayedArray"         "BiocParallel"        
    ## [40] "matrixStats"          "Biobase"              "GenomicRanges"       
    ## [43] "GenomeInfoDb"         "IRanges"              "S4Vectors"           
    ## [46] "BiocGenerics"         "parallel"             "stats4"              
    ## [49] "forcats"              "stringr"              "dplyr"               
    ## [52] "purrr"                "readr"                "tidyr"               
    ## [55] "tibble"               "ggplot2"              "tidyverse"           
    ## [58] "stats"                "graphics"             "grDevices"           
    ## [61] "utils"                "datasets"             "methods"             
    ## [64] "base"                
    ## 
    ## [[26]]
    ##  [1] "RColorBrewer"         "BioInstaller"         "SeqGSEA"             
    ##  [4] "DESeq"                "lattice"              "locfit"              
    ##  [7] "doParallel"           "iterators"            "foreach"             
    ## [10] "PROPER"               "biomaRt"              "pathview"            
    ## [13] "org.Hs.eg.db"         "formula.tools"        "sva"                 
    ## [16] "genefilter"           "mgcv"                 "nlme"                
    ## [19] "pheatmap"             "org.Rn.eg.db"         "plyr"                
    ## [22] "AnnotationHub"        "BiocFileCache"        "dbplyr"              
    ## [25] "BSgenome"             "rtracklayer"          "Biostrings"          
    ## [28] "XVector"              "GenomicFeatures"      "hexbin"              
    ## [31] "vsn"                  "GO.db"                "AnnotationDbi"       
    ## [34] "rafalib"              "devtools"             "usethis"             
    ## [37] "DESeq2"               "SummarizedExperiment" "DelayedArray"        
    ## [40] "BiocParallel"         "matrixStats"          "Biobase"             
    ## [43] "GenomicRanges"        "GenomeInfoDb"         "IRanges"             
    ## [46] "S4Vectors"            "BiocGenerics"         "parallel"            
    ## [49] "stats4"               "forcats"              "stringr"             
    ## [52] "dplyr"                "purrr"                "readr"               
    ## [55] "tidyr"                "tibble"               "ggplot2"             
    ## [58] "tidyverse"            "stats"                "graphics"            
    ## [61] "grDevices"            "utils"                "datasets"            
    ## [64] "methods"              "base"                
    ## 
    ## [[27]]
    ##  [1] "lubridate"            "RColorBrewer"         "BioInstaller"        
    ##  [4] "SeqGSEA"              "DESeq"                "lattice"             
    ##  [7] "locfit"               "doParallel"           "iterators"           
    ## [10] "foreach"              "PROPER"               "biomaRt"             
    ## [13] "pathview"             "org.Hs.eg.db"         "formula.tools"       
    ## [16] "sva"                  "genefilter"           "mgcv"                
    ## [19] "nlme"                 "pheatmap"             "org.Rn.eg.db"        
    ## [22] "plyr"                 "AnnotationHub"        "BiocFileCache"       
    ## [25] "dbplyr"               "BSgenome"             "rtracklayer"         
    ## [28] "Biostrings"           "XVector"              "GenomicFeatures"     
    ## [31] "hexbin"               "vsn"                  "GO.db"               
    ## [34] "AnnotationDbi"        "rafalib"              "devtools"            
    ## [37] "usethis"              "DESeq2"               "SummarizedExperiment"
    ## [40] "DelayedArray"         "BiocParallel"         "matrixStats"         
    ## [43] "Biobase"              "GenomicRanges"        "GenomeInfoDb"        
    ## [46] "IRanges"              "S4Vectors"            "BiocGenerics"        
    ## [49] "parallel"             "stats4"               "forcats"             
    ## [52] "stringr"              "dplyr"                "purrr"               
    ## [55] "readr"                "tidyr"                "tibble"              
    ## [58] "ggplot2"              "tidyverse"            "stats"               
    ## [61] "graphics"             "grDevices"            "utils"               
    ## [64] "datasets"             "methods"              "base"                
    ## 
    ## [[28]]
    ##  [1] "hms"                  "lubridate"            "RColorBrewer"        
    ##  [4] "BioInstaller"         "SeqGSEA"              "DESeq"               
    ##  [7] "lattice"              "locfit"               "doParallel"          
    ## [10] "iterators"            "foreach"              "PROPER"              
    ## [13] "biomaRt"              "pathview"             "org.Hs.eg.db"        
    ## [16] "formula.tools"        "sva"                  "genefilter"          
    ## [19] "mgcv"                 "nlme"                 "pheatmap"            
    ## [22] "org.Rn.eg.db"         "plyr"                 "AnnotationHub"       
    ## [25] "BiocFileCache"        "dbplyr"               "BSgenome"            
    ## [28] "rtracklayer"          "Biostrings"           "XVector"             
    ## [31] "GenomicFeatures"      "hexbin"               "vsn"                 
    ## [34] "GO.db"                "AnnotationDbi"        "rafalib"             
    ## [37] "devtools"             "usethis"              "DESeq2"              
    ## [40] "SummarizedExperiment" "DelayedArray"         "BiocParallel"        
    ## [43] "matrixStats"          "Biobase"              "GenomicRanges"       
    ## [46] "GenomeInfoDb"         "IRanges"              "S4Vectors"           
    ## [49] "BiocGenerics"         "parallel"             "stats4"              
    ## [52] "forcats"              "stringr"              "dplyr"               
    ## [55] "purrr"                "readr"                "tidyr"               
    ## [58] "tibble"               "ggplot2"              "tidyverse"           
    ## [61] "stats"                "graphics"             "grDevices"           
    ## [64] "utils"                "datasets"             "methods"             
    ## [67] "base"                
    ## 
    ## [[29]]
    ##  [1] "ggpubr"               "magrittr"             "hms"                 
    ##  [4] "lubridate"            "RColorBrewer"         "BioInstaller"        
    ##  [7] "SeqGSEA"              "DESeq"                "lattice"             
    ## [10] "locfit"               "doParallel"           "iterators"           
    ## [13] "foreach"              "PROPER"               "biomaRt"             
    ## [16] "pathview"             "org.Hs.eg.db"         "formula.tools"       
    ## [19] "sva"                  "genefilter"           "mgcv"                
    ## [22] "nlme"                 "pheatmap"             "org.Rn.eg.db"        
    ## [25] "plyr"                 "AnnotationHub"        "BiocFileCache"       
    ## [28] "dbplyr"               "BSgenome"             "rtracklayer"         
    ## [31] "Biostrings"           "XVector"              "GenomicFeatures"     
    ## [34] "hexbin"               "vsn"                  "GO.db"               
    ## [37] "AnnotationDbi"        "rafalib"              "devtools"            
    ## [40] "usethis"              "DESeq2"               "SummarizedExperiment"
    ## [43] "DelayedArray"         "BiocParallel"         "matrixStats"         
    ## [46] "Biobase"              "GenomicRanges"        "GenomeInfoDb"        
    ## [49] "IRanges"              "S4Vectors"            "BiocGenerics"        
    ## [52] "parallel"             "stats4"               "forcats"             
    ## [55] "stringr"              "dplyr"                "purrr"               
    ## [58] "readr"                "tidyr"                "tibble"              
    ## [61] "ggplot2"              "tidyverse"            "stats"               
    ## [64] "graphics"             "grDevices"            "utils"               
    ## [67] "datasets"             "methods"              "base"                
    ## 
    ## [[30]]
    ##  [1] "ggrepel"              "ggpubr"               "magrittr"            
    ##  [4] "hms"                  "lubridate"            "RColorBrewer"        
    ##  [7] "BioInstaller"         "SeqGSEA"              "DESeq"               
    ## [10] "lattice"              "locfit"               "doParallel"          
    ## [13] "iterators"            "foreach"              "PROPER"              
    ## [16] "biomaRt"              "pathview"             "org.Hs.eg.db"        
    ## [19] "formula.tools"        "sva"                  "genefilter"          
    ## [22] "mgcv"                 "nlme"                 "pheatmap"            
    ## [25] "org.Rn.eg.db"         "plyr"                 "AnnotationHub"       
    ## [28] "BiocFileCache"        "dbplyr"               "BSgenome"            
    ## [31] "rtracklayer"          "Biostrings"           "XVector"             
    ## [34] "GenomicFeatures"      "hexbin"               "vsn"                 
    ## [37] "GO.db"                "AnnotationDbi"        "rafalib"             
    ## [40] "devtools"             "usethis"              "DESeq2"              
    ## [43] "SummarizedExperiment" "DelayedArray"         "BiocParallel"        
    ## [46] "matrixStats"          "Biobase"              "GenomicRanges"       
    ## [49] "GenomeInfoDb"         "IRanges"              "S4Vectors"           
    ## [52] "BiocGenerics"         "parallel"             "stats4"              
    ## [55] "forcats"              "stringr"              "dplyr"               
    ## [58] "purrr"                "readr"                "tidyr"               
    ## [61] "tibble"               "ggplot2"              "tidyverse"           
    ## [64] "stats"                "graphics"             "grDevices"           
    ## [67] "utils"                "datasets"             "methods"             
    ## [70] "base"                
    ## 
    ## [[31]]
    ##  [1] "ggrepel"              "ggpubr"               "magrittr"            
    ##  [4] "hms"                  "lubridate"            "RColorBrewer"        
    ##  [7] "BioInstaller"         "SeqGSEA"              "DESeq"               
    ## [10] "lattice"              "locfit"               "doParallel"          
    ## [13] "iterators"            "foreach"              "PROPER"              
    ## [16] "biomaRt"              "pathview"             "org.Hs.eg.db"        
    ## [19] "formula.tools"        "sva"                  "genefilter"          
    ## [22] "mgcv"                 "nlme"                 "pheatmap"            
    ## [25] "org.Rn.eg.db"         "plyr"                 "AnnotationHub"       
    ## [28] "BiocFileCache"        "dbplyr"               "BSgenome"            
    ## [31] "rtracklayer"          "Biostrings"           "XVector"             
    ## [34] "GenomicFeatures"      "hexbin"               "vsn"                 
    ## [37] "GO.db"                "AnnotationDbi"        "rafalib"             
    ## [40] "devtools"             "usethis"              "DESeq2"              
    ## [43] "SummarizedExperiment" "DelayedArray"         "BiocParallel"        
    ## [46] "matrixStats"          "Biobase"              "GenomicRanges"       
    ## [49] "GenomeInfoDb"         "IRanges"              "S4Vectors"           
    ## [52] "BiocGenerics"         "parallel"             "stats4"              
    ## [55] "forcats"              "stringr"              "dplyr"               
    ## [58] "purrr"                "readr"                "tidyr"               
    ## [61] "tibble"               "ggplot2"              "tidyverse"           
    ## [64] "stats"                "graphics"             "grDevices"           
    ## [67] "utils"                "datasets"             "methods"             
    ## [70] "base"                
    ## 
    ## [[32]]
    ##  [1] "qvalue"               "ggrepel"              "ggpubr"              
    ##  [4] "magrittr"             "hms"                  "lubridate"           
    ##  [7] "RColorBrewer"         "BioInstaller"         "SeqGSEA"             
    ## [10] "DESeq"                "lattice"              "locfit"              
    ## [13] "doParallel"           "iterators"            "foreach"             
    ## [16] "PROPER"               "biomaRt"              "pathview"            
    ## [19] "org.Hs.eg.db"         "formula.tools"        "sva"                 
    ## [22] "genefilter"           "mgcv"                 "nlme"                
    ## [25] "pheatmap"             "org.Rn.eg.db"         "plyr"                
    ## [28] "AnnotationHub"        "BiocFileCache"        "dbplyr"              
    ## [31] "BSgenome"             "rtracklayer"          "Biostrings"          
    ## [34] "XVector"              "GenomicFeatures"      "hexbin"              
    ## [37] "vsn"                  "GO.db"                "AnnotationDbi"       
    ## [40] "rafalib"              "devtools"             "usethis"             
    ## [43] "DESeq2"               "SummarizedExperiment" "DelayedArray"        
    ## [46] "BiocParallel"         "matrixStats"          "Biobase"             
    ## [49] "GenomicRanges"        "GenomeInfoDb"         "IRanges"             
    ## [52] "S4Vectors"            "BiocGenerics"         "parallel"            
    ## [55] "stats4"               "forcats"              "stringr"             
    ## [58] "dplyr"                "purrr"                "readr"               
    ## [61] "tidyr"                "tibble"               "ggplot2"             
    ## [64] "tidyverse"            "stats"                "graphics"            
    ## [67] "grDevices"            "utils"                "datasets"            
    ## [70] "methods"              "base"

``` r
############################################################
##### Functions ############################################
############################################################

# Set select
select <- dplyr::select

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

##### Data files to load:

  - Count Matrix and Metadata Table from:
      - RNA-Seq from Mt. Sinai
          - 3 sequencing batches & metadata
      - RNA-Seq from Stanford
          - 2 sequencing batches &
metadata

<!-- end list -->

``` r
################################################################################
#####     Load & Clean Data      ###############################################
################################################################################

# Files last saved in: 20200309_exploration-rna-seq-phase1_steep.R

# Count matrix
in_file <- paste0(WD,'/data/20200326_rnaseq-countmatrix-pass1a-stanford-sinai_steep.csv')
count_data <- read.table(in_file,sep = ',', header = TRUE,row.names = 1,check.names = FALSE)

# Meatdata table
in_file <- paste0(WD,'/data/20200326_rnaseq-meta-pass1a-stanford-sinai_steep.txt')
col_data <- read.table(in_file, header = TRUE, check.names = FALSE, sep = '\t')
row.names(col_data) <- col_data$sample_key
```

## Place Genes in Genomic Ranges

#### Reference Genome and Annotation: Rnor\_6.0 (GCA\_000001895.4) assembly from Ensembl database (Release 96)

Found at: <http://uswest.ensembl.org/Rattus_norvegicus/Info/Index>.
FASTA: Rattus\_norvegicus.Rnor\_6.0.dna.toplevel.fa.gz
<ftp://ftp.ensembl.org/pub/release-96/fasta/rattus_norvegicus/dna/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz>
GTF: Rattus\_norvegicus.Rnor\_6.0.96.gtf.gz
<ftp://ftp.ensembl.org/pub/release-96/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.96.gtf.gz>

``` r
################################################################################
#####     Annotate Genes by Chromosome       ###################################
################################################################################

# TODO: Load Reference genome and annotations

### Determine which control samples are male and female
# Get the list of genes on the W chromosome

# Construct your own personal galgal5 reference genome annotation
# Construct from gtf file from Ensembl (same file used in mapping)
ens_gtf <- paste0(WD,'/data/Rattus_norvegicus.Rnor_6.0.96.gtf')
Rn_TxDb <- makeTxDbFromGFF(ens_gtf,
                           format=c("gtf"),
                           dataSource="Ensembl_Rattus6_gtf",
                           organism="Rattus norvegicus",
                           taxonomyId=NA,
                           circ_seqs=DEFAULT_CIRC_SEQS,
                           chrominfo=NULL,
                           miRBaseBuild=NA,
                           metadata=NULL)

# Define Female specific sex genes (X chromosome)
# To examine chromosome names
seqlevels(Rn_TxDb)[1:35]
```

    ##  [1] "1"              "2"              "3"              "4"             
    ##  [5] "5"              "6"              "7"              "8"             
    ##  [9] "9"              "10"             "11"             "12"            
    ## [13] "13"             "14"             "15"             "16"            
    ## [17] "17"             "18"             "19"             "20"            
    ## [21] "X"              "Y"              "MT"             "AABR07022258.1"
    ## [25] "AABR07022620.1" "AABR07022926.1" "AABR07024031.1" "AABR07024032.1"
    ## [29] "AABR07024040.1" "AABR07024041.1" "AABR07024044.1" "AABR07024102.1"
    ## [33] "AABR07024104.1" "AABR07024106.1" "AABR07024115.1"

``` r
# Extract genes as GRanges object, then names
X_genes_gr <- genes(Rn_TxDb, columns = "TXCHROM", filter = list(tx_chrom=c("X")))
# Collect ensembl gene ids for female specific genes
X_ens_id <- names(X_genes_gr)
# Examine the gene symbols
X_sym <- mapIds(org.Rn.eg.db, names(X_genes_gr), "SYMBOL", "ENSEMBL")
# Extract genes as GRanges object, then names
Y_genes_gr <- genes(Rn_TxDb, columns = "TXCHROM", filter = list(tx_chrom=c("Y")))
# Collect ensembl gene ids for female specific genes
Y_ens_id <- names(Y_genes_gr)
sex_ens_id <- c(X_ens_id,Y_ens_id)
# Examine the gene symbols
Y_sym <- mapIds(org.Rn.eg.db, names(Y_genes_gr), "SYMBOL", "ENSEMBL")

################################################################################
```

## PCA of Liver Samples

Samples from Stanford Batch 1, which we suspect demonstrates a batch
effect

``` r
################################################################################
#####     PCA of Liver Samples       ###########################################
################################################################################

design = ~ 1 # Primary variable needs to be last
title = paste0('Design: ',as.character(design))
# Create a DESeqDataSet Object
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = col_data,
                              design = design)

# Perform pre-filtering.
# Filter genes with average count of 10 or less.
# Reasoning from:
#citation("PROPER")
#dds
keep <- rowSums(counts(dds))/ncol(dds) >= 10
dds <- dds[keep,]
```

#### Summary of counts and annotation data in a DESeqDataSet

``` r
dds
```

    ## class: DESeqDataSet 
    ## dim: 24387 1360 
    ## metadata(1): version
    ## assays(1): counts
    ## rownames(24387): ENSRNOG00000055717 ENSRNOG00000061350 ...
    ##   ENSRNOG00000056235 ENSRNOG00000058021
    ## rowData names(0):
    ## colnames(1360): 90110015502_SN1 90112015502_SN1 ... 80001995527_SF2
    ##   80001995535_SF2
    ## colData names(198): vial_label X2D_barcode ...
    ##   calculated.variables.time_to_freeze Seq_batch

``` r
dds <- estimateSizeFactors(dds)
rs <- rowSums(counts(dds))
# Normalize the counts
rld <- vst(dds) #vst and rlog comparable with all samples
#rld <- rlog(dds, blind=FALSE)

# Extract matrix of normalized counts
norm_counts <- assay(rld)
counts <- as.data.frame(norm_counts)

# Annotate normalized counts
counts$ensembl <- rownames(counts)
counts$symbol <- mapIds(org.Rn.eg.db, counts$ensembl, "SYMBOL", "ENSEMBL")
counts$entrez <- mapIds(org.Rn.eg.db, counts$ensembl, "ENTREZID", "ENSEMBL")
counts$genename <- mapIds(org.Rn.eg.db, counts$ensembl, "GENENAME", "ENSEMBL")
counts$go <- mapIds(org.Rn.eg.db, counts$ensembl, "GO", "ENSEMBL")
counts$path <- mapIds(org.Rn.eg.db, counts$ensembl, "PATH", "ENSEMBL")
```

#### Liver samples cluster by sex.

Grey samples represent “reference samples”.

``` r
rld.sub <- rld[ , (rld$Tissue == "Liver") ]
DESeq2::plotPCA(rld.sub, intgroup ="animal.registration.sex") +
        guides(color=guide_legend(title="Sex"))
```

![](20200325_circadian-batch-investigation-rna-seq_steep_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
# Variables of interest
male_livers <- (col_data %>% 
                        filter(Tissue == 'Liver') %>% 
                        filter(animal.registration.sex == 'Male'))$sample_key %>% 
        as.character()
female_livers <- (col_data %>% 
                          filter(Tissue == 'Liver') %>% 
                          filter(animal.registration.sex == 'Female'))$sample_key %>% 
        as.character()
ref_livers <- (col_data %>% 
                       filter(Tissue == 'Liver') %>% 
                       filter(is.na(animal.registration.sex)))$sample_key %>% 
        as.character()
livers <- c(male_livers,female_livers,ref_livers)
Y_genes <- Y_ens_id[Y_ens_id %in% row.names(norm_counts)]
X_genes <- X_ens_id[X_ens_id %in% row.names(norm_counts)]
sex <- col_data[livers,"animal.registration.sex"]
group <- col_data[livers,"animal.key.anirandgroup"]
liver_counts <- norm_counts[,livers]
```

#### Predict the sex of reference samples (all samples for that matter) by calculating the median expression of genes on the Y chromosome. We should expect a bimodal distribution with males demonstrating significantly higher median expression.

``` r
chryexp <- colMeans(norm_counts[Y_genes,livers])
```

##### If we create a histogram of the median gene expression values on chromosome Y, we should expect to see a bimodal distribution. However, distinct peaks are not detected. This was a surprising result.

``` r
mypar()
hist(chryexp, breaks = 200)
```

![](20200325_circadian-batch-investigation-rna-seq_steep_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
summary(chryexp)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##   2.235   2.474   2.544   2.587   2.681   3.148

We will not use this common strategy to determine sex of unknown
samples, rather we will use clustering from PCA. The distribution of sex
by group

``` r
table(group, sex) # <- Bad idea.
```

    ##                    sex
    ## group               Female Male
    ##   Control - 7 hr         5    5
    ##   Control - IPE          5    5
    ##   Exercise - 0.5 hr      5    5
    ##   Exercise - 1 hr        5    5
    ##   Exercise - 24 hr       3    3
    ##   Exercise - 4 hr        5    5
    ##   Exercise - 48 hr       3    3
    ##   Exercise - 7 hr        3    3
    ##   Exercise - IPE         5    5

#### Generate a heatmap of expression from 3 sets of genes:

  - Genes from the Y chromosome
  - The top and bottom 25 genes (50 total) associated with sex
  - Randomly selected genes \#\#\#\#\# Males and Females demonstrate
    distinctly different gene expression profiles.
  - Genes on the Y chromosome are not a good predictor of sex in Liver
    mRNA measures (Figure 1; suprising result)
  - Male and female samples show distinct correlation to one another
    (Figure 2)

<!-- end list -->

``` r
# T-test of expression associated with sex
tt <- rowttests(liver_counts,sex)

# Take genes from the Y chromosome
# Y_genes
# Take the top and bottom 25 genes associated with variable of interest (remove any genes in Y chromosome)
top <- row.names(tt[order(-tt$dm),][1:25,])
bot <- row.names(tt[order(tt$dm),][1:25,])
top_n_bot <- setdiff(c(top,bot), Y_genes)

# Randomly select 50 genes not in prior sets 
set.seed(123)
randos <- setdiff(row.names(tt[sample(seq(along=tt$dm),50),]), c(Y_genes,top_n_bot))
geneindex <- c(randos,top_n_bot,Y_genes)
# Generate the heatmap and support with a plot of a correlation matrix
mat <- liver_counts[geneindex,]
mat <- mat -rowMeans(mat)
icolors <- colorRampPalette(rev(brewer.pal(11,"RdYlBu")))(100)
mypar(1,2)
image(t(mat),xaxt="n",yaxt="n",col=icolors)
y <- liver_counts - rowMeans(liver_counts)
image(1:ncol(y),1:ncol(y),cor(y),col=icolors,zlim=c(-1,1),
      xaxt="n",xlab="",yaxt="n",ylab="")
axis(2,1:ncol(y),sex,las=2)
axis(1,1:ncol(y),sex,las=2)
```

![](20200325_circadian-batch-investigation-rna-seq_steep_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

#### A naive t-test and genes with q values less than or equal to 0.05

##### The left figure represents a histogram of p values from a naive t-test (all genes). We see that a significant proportion of genes correlate with sex. To investigate if these genes are located on sex chromosomes or autosomal chromosomes, we contruct a volcano plot on the right. To our suprise, again, genes on the Y chromosome do not show signifcant correlation to sex. Rather some genes on the X chromosome demonstrate significance, however, these genes do not make up the majority of significantly assocaited genes.

## The variance in gene expression is not dicated by differential expression of genes on sex chromosomes.

``` r
mypar(1,2)
# Histogram of p values associated with ttest
hist(tt$p.value,main="",ylim=c(0,1300), breaks = 100)
plot(tt$dm,-log10(tt$p.value))
tt[X_genes,]$dm
```

    ##   [1]  0.0000000000  0.0000000000 -0.1039153165  0.0000000000  0.0379719804
    ##   [6]  0.0000000000  0.2729022505  0.0000000000  0.0090663243  0.0000000000
    ##  [11]  0.0000000000 -0.0468708207  0.0957258983  0.0000000000 -0.1326385623
    ##  [16] -0.0766205772 -0.0438796928  0.0000000000  0.0000000000  0.0000000000
    ##  [21] -0.0375551405  0.0000000000  0.0000000000  0.0000000000 -0.6870205409
    ##  [26] -0.1902717589  0.0000000000  0.0000000000  0.0000000000  0.0000000000
    ##  [31]  0.0000000000  0.0000000000  0.0000000000 -0.8836673599  0.0000000000
    ##  [36]  0.0000000000  0.0000000000  0.0000000000  0.0000000000 -0.0574059067
    ##  [41]  0.0000000000  0.0585261001  0.0000000000  0.0000000000  0.0000000000
    ##  [46]  0.0000000000 -0.0548102222  0.0000000000  0.0000000000  0.0000000000
    ##  [51]  0.0000000000  0.0000000000  0.0000000000  0.0000000000 -0.0611312200
    ##  [56]  0.0028500494  0.0000000000  0.0000000000  0.0000000000  0.0000000000
    ##  [61]  0.0345970569 -0.3436007611 -0.3550732746  0.0000000000  0.0000000000
    ##  [66]  0.0000000000 -0.2223737321  0.0151143412  0.0000000000 -5.2669228431
    ##  [71] -0.0507812022  0.2549965751  0.0465547209  0.0000000000  0.0000000000
    ##  [76]  0.0000000000  0.0000000000  0.0146976078  0.1031009651  0.3001476135
    ##  [81] -0.1669416501  0.2753667831  0.0000000000  0.0000000000 -0.0715319398
    ##  [86] -0.0244025976 -0.4383774046 -0.0337727911  0.0000000000  0.0000000000
    ##  [91]  0.0000000000  0.0000000000  0.1082666741  0.0000000000  0.0000000000
    ##  [96] -0.3056478448  0.0000000000  0.1921325012  0.0000000000  0.0000000000
    ## [101]  0.0000000000  0.5552954575  0.0000000000  0.0000000000 -0.0509280814
    ## [106]  0.0000000000  0.0000000000  0.0000000000  0.0000000000  0.0000000000
    ## [111] -0.3497410644  0.0000000000  0.0000000000  0.0000000000  0.0000000000
    ## [116]  0.0000000000 -0.0360358193  0.0000000000 -0.0600618579  0.0000000000
    ## [121]  0.0000000000  0.0262442528  0.0576143582  0.0000000000 -0.1171137013
    ## [126]  0.0000000000  0.1171643895 -0.0413555398  0.0000000000  0.0000000000
    ## [131]  0.0000000000  0.0000000000 -0.0198449923  0.0000000000  0.0000000000
    ## [136] -0.1305840299  0.0000000000 -0.2906095136  0.0256309800  0.0000000000
    ## [141]  0.0000000000 -0.0115171693  0.0000000000  0.0000000000  0.0000000000
    ## [146] -0.4202809285  0.0758108208  0.1306404605  0.0000000000  0.0000000000
    ## [151]  0.0000000000  0.0000000000  0.0000000000  0.0000000000  0.0000000000
    ## [156]  0.0000000000 -0.1375695529  0.0000000000 -0.1119697103  0.0000000000
    ## [161]  0.0000000000 -0.0619582361  0.0000000000  0.0000000000  0.1182901538
    ## [166]  0.0000000000 -0.0294378361  0.0000000000  0.0095563282 -0.0032299686
    ## [171]  0.0000000000 -0.0101399520  0.0000000000  0.0104990652 -1.2680471915
    ## [176]  0.0224760024  0.0000000000  0.0000000000  0.0000000000  0.0000000000
    ## [181]  0.0989494797  0.0000000000  0.0000000000  0.0000000000 -0.0079668698
    ## [186]  0.0000000000  0.3445441619 -0.3248552857  0.0000000000  0.0000000000
    ## [191] -0.6602305244  0.0000000000  0.5340939953  0.0000000000  0.0133974655
    ## [196]  0.0000000000 -0.0259748178  0.0000000000 -0.2039862472  0.0000000000
    ## [201]  0.1282535986  0.0000000000  0.0000000000  0.0000000000  0.0000000000
    ## [206]  0.0000000000 -0.1010450036  0.0000000000  0.0000000000  0.0000000000
    ## [211]  0.0000000000  0.0000000000  0.1909344831  0.0000000000  0.0000000000
    ## [216]  0.0000000000  0.0000000000  0.2663201212  0.0000000000  0.0000000000
    ## [221]  0.0000000000 -0.3167976859  0.0578249562  0.0837854466  0.0052218998
    ## [226]  0.0000000000  0.0000000000  0.0000000000 -0.1588294361  0.0000000000
    ## [231]  0.0000000000  0.0000000000  0.0000000000  0.0000000000 -0.2186313887
    ## [236]  0.0000000000  0.0000000000  0.1232354517  0.0000000000  0.0000000000
    ## [241] -0.0566479373 -0.0321718774  0.0000000000  0.0000000000 -0.0277042566
    ## [246]  0.0000000000  0.0000000000  0.0000000000  0.0590360466 -0.2739508947
    ## [251]  0.0000000000  0.0000000000 -0.8066811795  0.0302125774 -0.0103756884
    ## [256]  0.0000000000  0.0000000000  0.1056613938  0.0000000000  0.0000000000
    ## [261]  0.0000000000 -0.0002412377 -0.0452656791  0.0000000000  0.0000000000
    ## [266]  0.0000000000  0.0000000000  0.0599569583  0.0000000000  0.0000000000
    ## [271] -3.0690100211  0.0000000000  0.0000000000 -0.6175598867  0.0000000000
    ## [276]  0.0000000000  0.0000000000  1.6601184640  0.0000000000  0.0000000000
    ## [281]  0.0000000000  0.0000000000  0.0000000000  0.1359360778  0.0000000000
    ## [286]  0.3964260760  0.0000000000 -0.0397608099 -0.0504196740  0.1354443927
    ## [291]  0.0000000000  0.0270900018  0.0000000000  0.0000000000  0.8930077453
    ## [296] -0.0967430864  0.0268500545  0.0000000000 -0.2329754269  0.0000000000
    ## [301]  0.0301124080 -0.1135872374  0.0325935588  0.0000000000 -0.4514452227
    ## [306]  0.0000000000  0.1521377665  0.0000000000  0.0000000000  0.0000000000
    ## [311]  0.0000000000 -0.1902691580  0.1137400744  0.0000000000  0.0974242473
    ## [316]  0.0346436696  0.0000000000  0.0000000000  0.0251012626  0.0000000000
    ## [321] -0.6126456274  0.0000000000  0.0000000000  0.0023178493  0.0000000000
    ## [326]  0.5258208808  0.0000000000  0.2990772774 -0.2628921592  0.0000000000
    ## [331]  0.0000000000  0.0000000000  0.0000000000  0.0000000000 -0.3953918393
    ## [336]  0.0309106364  0.0000000000  0.0000000000  0.0000000000  0.0000000000
    ## [341]  0.0000000000  0.0000000000  0.4322811951 -0.0884338200  0.0000000000
    ## [346]  0.0000000000  0.0000000000  0.0000000000  0.0000000000  0.2211070967
    ## [351] -0.0594022598 -0.1749966387  0.0000000000  0.0000000000  0.0000000000
    ## [356]  0.0877774780 -0.0313675380  0.1918161246  0.0000000000 -0.0655144519
    ## [361]  0.0000000000  0.0000000000  0.0522275739  0.0000000000  0.0000000000
    ## [366]  0.0568206621  0.0104662118  0.0000000000  0.0000000000  0.2099298387
    ## [371]  0.0000000000  0.0000000000 -2.0061147640  0.0000000000  0.0000000000
    ## [376]  0.0000000000  0.0000000000 -0.1681827070  0.0000000000  0.0000000000
    ## [381] -5.0461271149  0.1352099630  0.0002441168  0.0000000000  0.1252666787
    ## [386]  0.0000000000  0.0000000000  0.0000000000  0.0000000000 -0.0430604364
    ## [391]  0.0000000000  0.0000000000 -0.6778237121  0.1335646748  1.6402011639
    ## [396] -0.0865031320  0.0124401126 -0.1140577914  0.0000000000  0.0000000000
    ## [401]  0.0000000000  0.0000000000  0.0000000000  0.0000000000  0.0000000000
    ## [406]  0.0000000000  0.0000000000  0.0000000000  0.0000000000  0.0000000000
    ## [411] -0.0280570317  0.0000000000  0.0000000000  0.0000000000  0.0000000000
    ## [416]  0.3590681975  0.0000000000 -0.6985076363 -0.2208423608  0.0000000000
    ## [421] -0.0389761711  0.2168411561 -0.1192214790  0.0000000000 -0.0260805715
    ## [426]  0.0311504696  0.0000000000 -0.0446014448  2.6289715621 -0.0598705828
    ## [431]  0.0000000000  0.0000000000  0.0000000000 -0.0258131014  0.0000000000
    ## [436] -0.0195542750 -0.0861681057  0.0000000000  0.0000000000  0.0000000000
    ## [441]  0.0000000000  1.1494220269  0.0000000000 -0.1279601628  0.0000000000
    ## [446]  0.0000000000  0.0000000000  0.0000000000  0.0000000000  0.0000000000
    ## [451]  0.0000000000  0.0000000000  0.0000000000 -0.0259986691 -0.0294937037
    ## [456]  0.0078147530  0.0000000000  0.0000000000  0.0000000000  0.0000000000
    ## [461]  0.1626438351  0.7995038076  0.1955225775  0.0725636457  0.0000000000
    ## [466] -0.0826045596  0.2000120072  0.0000000000  0.0387269036  0.3437195497
    ## [471]  0.0000000000  0.0000000000  0.2739109910 -0.0446722106  0.2193097025
    ## [476]  0.0000000000  0.0803624796  0.0000000000  0.0224760024 -0.0287026578
    ## [481]  0.0000000000  0.0656331196 -0.1559284816 -0.3880366414  0.1637542421
    ## [486]  0.0000000000  0.0000000000  0.0000000000  0.0000000000  0.0000000000
    ## [491] -0.0111487867  0.0000000000 -0.0304837547 -0.0243943064  0.0000000000
    ## [496]  0.0052394762  0.0000000000  0.0000000000  0.0000000000  0.0000000000
    ## [501]  0.0269479537 -0.0231664079  0.0000000000  0.0000000000  0.0404614295
    ## [506]  0.0000000000  0.0000000000  0.0391923665  0.3260843792  0.0000000000
    ## [511]  0.0000000000  0.0578584849  0.5949132414  0.4174785924  0.0000000000
    ## [516] -0.0841065580  0.0264003000  0.0000000000  0.0000000000  0.0513715099
    ## [521]  0.0000000000  0.0000000000  0.0000000000  0.0742358668 -0.0084239776
    ## [526]  0.0000000000  0.0000000000 -0.2839721667  0.0000000000  0.0000000000
    ## [531]  0.1475254323 -0.0368810501  0.0028320337  0.0000000000  0.0015908303
    ## [536]  0.0000000000  0.0000000000  0.0000000000  0.0000000000  0.0000000000
    ## [541]  0.2077909885 -0.5205517071  0.0000000000  0.0289302115  0.3293579158
    ## [546]  0.0000000000 -0.0054485849  0.0000000000  0.1331392713  0.0000000000
    ## [551]  0.0000000000 -0.3252253671 -0.3866067363 -0.0014444230 -8.8261882249
    ## [556] -0.2983164439 -0.0980596460 -0.0192369087  0.0000000000 -0.0165705526
    ## [561]  0.0874161657  0.2130711048  0.0000000000  0.3797547594 -0.3334440331
    ## [566]  0.0175309918  0.0182787576 -0.0743266689 -0.4965639791 -0.2326825720
    ## [571]  0.0000000000 -0.0252913320  0.0000000000 -0.2868867891  0.4048264991
    ## [576]  0.2672402473  0.0119924476 -0.1993889830  0.1854912117  0.2365519954
    ## [581]  0.0000000000 -1.5528460025  0.0000000000  0.0000000000 -0.0197574664
    ## [586]  0.0000000000  0.0000000000  0.0000000000 -0.4762158042 -0.0662516501
    ## [591]  0.0000000000  0.0000000000  0.0000000000  0.0000000000  0.0000000000
    ## [596]  0.0000000000 -0.1089151576 -0.0302072893 -0.0123553647  0.0050235031
    ## [601]  0.0000000000 -6.8707043844  0.0000000000  0.0000000000  0.0000000000
    ## [606] -0.1225708710 -0.1520939415  0.0000000000 -0.0268681072  0.0000000000
    ## [611] -0.0101624618  0.0034259090 -0.0083199542  0.0000000000  0.0444507916
    ## [616] -0.0358079757  0.0000000000  0.0786168035  0.0000000000  0.0000000000
    ## [621]  0.0000000000  0.0000000000 -0.0812863133  0.0478602545  0.1505321295
    ## [626]  0.0000000000 -0.8443772926  0.0000000000  0.0000000000  0.0000000000
    ## [631]  0.1597518804  0.1172773401  0.0000000000 -0.1212250782 -1.3990687293
    ## [636]  0.2668496337  0.0000000000 -0.3746240733  0.0531541124  0.0000000000
    ## [641]  0.0046441942  0.0000000000  0.2785890063  0.0000000000 -0.2572221645
    ## [646]  0.0000000000  0.0577546559  0.0001224378  0.0000000000  0.0000000000
    ## [651] -0.0903877376 -0.4901621827  0.0000000000  0.0000000000  0.1039930478
    ## [656]  0.3046228318  0.0000000000 -0.2542418962 -0.1152146748  0.1875926182
    ## [661]  0.0000000000  0.0000000000  0.0000000000  0.0000000000  0.1503422270
    ## [666]  0.0000000000  0.0000000000 -0.3635000460  0.0000000000  0.0812992222
    ## [671]  0.0076733475  0.0000000000  0.0000000000 -0.0941329838 -0.0173064124
    ## [676] -0.0057993331  0.0448362854 -0.0323703937  0.0000000000  0.2349788295
    ## [681]  0.0443367650  0.2369123025  0.0000000000  0.0463356908  0.2653225669
    ## [686]  0.0000000000 -0.1178015013  0.0000000000  0.0000000000 -0.3206173722
    ## [691]  0.0000000000  0.1002294694  0.1689969461 -0.6303033916  0.0000000000
    ## [696]  0.0000000000  0.0085955954  0.1661565733  0.0000000000  0.0000000000
    ## [701]  0.0590401323  0.0000000000 -0.0221827288 -0.0197470290 -0.4716907314
    ## [706]  0.0278861129 -0.0259010773  1.0661726537 -0.0735606355  0.0570613353
    ## [711]  0.0000000000  0.0000000000  0.0000000000  0.2317934116  0.0000000000
    ## [716]  0.0000000000  0.0000000000  0.0000000000 -1.2505711251  0.0361432436
    ## [721] -0.1350238685 -0.0358079757  0.0000000000 -0.2046959189  0.0000000000
    ## [726]  0.0000000000 -0.5893155868  0.0094806779  0.0000000000  0.0000000000
    ## [731] -0.0071278335 -0.1388081955  0.0000000000  0.0000000000 -0.0564622915
    ## [736] -0.1642395445  0.0000000000  0.0000000000 -0.1509411340 -0.6806953305
    ## [741]  0.0000000000  0.0000000000  0.0000000000 -0.1573815215  0.0000000000
    ## [746] -0.1324435183  0.0449994340  0.0005582181  0.0000000000  0.0639876621
    ## [751]  0.0000000000  0.4287640356  0.0346436696  0.0000000000  0.0000000000
    ## [756]  0.0000000000  0.0000000000  0.1393766273 -0.1955146850  0.0000000000
    ## [761]  0.0000000000  0.0183504569  0.0000000000  0.0000000000  0.0000000000
    ## [766]  0.0000000000 -0.2397521406  0.0000000000 -0.3731514300  0.0000000000
    ## [771]  0.0000000000  0.0000000000  0.0503773470  0.8451323688 -0.5425294637

``` r
tt[Y_genes,]
```

    ##                     statistic          dm      p.value
    ## ENSRNOG00000052020 -6.8779504 -0.81074302 1.489674e-09
    ## ENSRNOG00000052326        NaN  0.00000000          NaN
    ## ENSRNOG00000056639 -0.4465435 -0.08239261 6.564745e-01
    ## ENSRNOG00000057231        NaN  0.00000000          NaN
    ## ENSRNOG00000058664  1.1647418  0.19777679 2.477665e-01
    ## ENSRNOG00000060048 -1.7496886 -0.10319542 8.420856e-02
    ## ENSRNOG00000060437  0.8293016  0.14304904 4.095309e-01
    ## ENSRNOG00000060496  0.6435822  0.06007540 5.217837e-01
    ## ENSRNOG00000060617 -0.4817986 -0.03505260 6.313336e-01
    ## ENSRNOG00000062104 -1.5910895 -0.46900729 1.157407e-01
    ## ENSRNOG00000062115        NaN  0.00000000          NaN
    ## ENSRNOG00000062152 -1.0000000 -0.05692271 3.204838e-01
    ## ENSRNOG00000062172        NaN  0.00000000          NaN

``` r
points(tt[X_genes,]$dm,-log10(tt[X_genes,]$p.value),col=1,pch=16)
points(tt[Y_genes,]$dm,-log10(tt[Y_genes,]$p.value),col=2,pch=16, xlab="Effect size",ylab="-log10(p-value)")
legend("bottomright",c("X","Y"),col=1:2,pch=16)
as.vector(tt$p.value)
```

    ##     [1] 5.771059e-01 3.149417e-02 5.606650e-03 2.636152e-01 1.808643e-03
    ##     [6] 2.443304e-02 2.421932e-04 3.204838e-01 1.046794e-02 9.565892e-01
    ##    [11] 8.891084e-02 4.243865e-02 9.394027e-01 5.623663e-01 1.264133e-01
    ##    [16] 2.396761e-01 2.030320e-04 7.326915e-05 9.019484e-01 7.643274e-01
    ##    [21] 1.807860e-01 3.545321e-01 7.511736e-01 5.806583e-04 3.985102e-02
    ##    [26] 1.502629e-04 5.968553e-01 7.086574e-12 2.694884e-16 2.536247e-01
    ##    [31] 5.851233e-03 1.079261e-17 1.066090e-02 4.286805e-01 4.682130e-09
    ##    [36] 1.003321e-02 1.953870e-03          NaN 3.289540e-01 1.256899e-03
    ##    [41] 3.289808e-01 8.509993e-02 2.485876e-11 9.217547e-01 4.924328e-01
    ##    [46] 7.668798e-01 1.937902e-01 9.510309e-05 6.340490e-01 6.819651e-01
    ##    [51] 9.436754e-02 7.418622e-07 3.255839e-01 2.169765e-03 4.346458e-02
    ##    [56] 7.967984e-01 7.766967e-02 1.247151e-01 2.326771e-02 5.025178e-06
    ##    [61] 6.426388e-16 4.237496e-02 3.557298e-01 3.966763e-01 4.421414e-01
    ##    [66] 1.627181e-11 1.098841e-03 2.478843e-03 4.806788e-05 1.303020e-17
    ##    [71] 6.238199e-15 2.872690e-02 6.533508e-01 2.793094e-01 2.851058e-02
    ##    [76] 3.836132e-02 3.204838e-01 2.647359e-06 1.066920e-03 2.389572e-11
    ##    [81] 5.883383e-01 3.113547e-01 3.498899e-01          NaN 7.259475e-05
    ##    [86] 9.654992e-01 1.535403e-02 6.931891e-01 1.932020e-01 1.679456e-04
    ##    [91] 4.598207e-01 1.509612e-04 2.982918e-05 6.898622e-01 7.287055e-01
    ##    [96] 3.465839e-10 2.378898e-07 1.076582e-04 6.723886e-06 4.273880e-01
    ##   [101] 5.210676e-01 3.381439e-11 4.182180e-01 3.204838e-01 8.227591e-01
    ##   [106] 8.837380e-01 6.122873e-01 2.830701e-13 8.730050e-01 1.319474e-02
    ##   [111] 3.563094e-05 2.730983e-01 4.395277e-01 8.645382e-02 8.617932e-06
    ##   [116] 6.817691e-01 8.507743e-02 5.651756e-03 8.373380e-06 1.301273e-02
    ##   [121] 1.666220e-07 2.343831e-01 1.171542e-02 9.428038e-01 6.193954e-01
    ##   [126] 1.384798e-01 1.526049e-01 7.673199e-01 6.500021e-01 6.825910e-07
    ##   [131]          NaN 1.775553e-04 1.320055e-01 7.437835e-01 1.480906e-08
    ##   [136] 6.140041e-02 9.579819e-01 2.757974e-01 6.960543e-02 6.573083e-01
    ##   [141]          NaN 2.016346e-04 6.266847e-02 7.630051e-01 1.483594e-01
    ##   [146] 6.189442e-06 8.220411e-02 9.648012e-01 3.182952e-01 9.140341e-02
    ##   [151] 9.403091e-03 2.737065e-05 5.902934e-01 1.262154e-05 6.213953e-01
    ##   [156] 3.102108e-02 3.044454e-01 2.611186e-01 3.672536e-01 7.543220e-05
    ##   [161] 2.643866e-01 8.952871e-01 4.771981e-14 2.599566e-04 3.261379e-14
    ##   [166] 3.141590e-10 1.249302e-01 7.919993e-12 2.503291e-02 1.698268e-01
    ##   [171] 4.602869e-11 5.019198e-01 6.356240e-04 1.573793e-01 1.466237e-01
    ##   [176] 6.483346e-02 2.406468e-07 3.604139e-01 2.152811e-01 1.547236e-02
    ##   [181] 7.290784e-04 1.294303e-01 2.036231e-03 7.238001e-01 8.027244e-03
    ##   [186] 4.768035e-04 1.108068e-01 2.596755e-03 2.174326e-01 3.432180e-05
    ##   [191] 1.081932e-02 1.863543e-01 6.189663e-01 1.882681e-01 1.240051e-01
    ##   [196] 3.275927e-01 6.437097e-02 6.080148e-02 8.573211e-12 2.673696e-01
    ##   [201] 7.930214e-01 1.091675e-03 3.705643e-02 7.433295e-01 7.532626e-01
    ##   [206] 2.331084e-04 9.206281e-01 1.758148e-19 3.589152e-05 1.596673e-01
    ##   [211]          NaN 2.113795e-01 5.151665e-01 4.628786e-15 1.648544e-09
    ##   [216] 5.538685e-01 6.804242e-01 1.575357e-01 4.257991e-02 4.464676e-01
    ##   [221] 3.440914e-11 3.594547e-01 6.546715e-01 6.712478e-01 8.197425e-01
    ##   [226] 5.231942e-01 4.056736e-02 2.495342e-08 6.273934e-05 8.419540e-01
    ##   [231] 2.587411e-03 3.006290e-09 5.859945e-04 5.163355e-11 3.252500e-01
    ##   [236] 8.026084e-01 5.344311e-01 2.353774e-01 5.028216e-05 9.533881e-01
    ##   [241] 3.644028e-05 4.729979e-06 6.632889e-02 4.195384e-02 6.243574e-02
    ##   [246] 4.294053e-02 8.852388e-01 6.335414e-01 5.209046e-02 3.190038e-04
    ##   [251] 6.746386e-01 2.126244e-01 9.407957e-01 4.107241e-08 4.861985e-01
    ##   [256] 1.042913e-02 2.820667e-13 7.064744e-02 4.645645e-03 2.040084e-09
    ##   [261] 2.378555e-05 7.827422e-01 1.650577e-05 3.701557e-02 7.929366e-01
    ##   [266] 4.819561e-01 2.094688e-19 1.937539e-02 2.841758e-01 5.456654e-02
    ##   [271] 1.067668e-02 7.207443e-01 3.175855e-01 9.955556e-02 8.840296e-14
    ##   [276] 1.962231e-01 7.346780e-01 3.387985e-01 9.848625e-04 9.645690e-05
    ##   [281] 1.378084e-13 2.485257e-03 6.250992e-01 5.962994e-02 2.057906e-10
    ##   [286] 4.077672e-06 1.950124e-02 1.246451e-01 4.730869e-01 3.036819e-01
    ##   [291] 1.516223e-06 1.102058e-04 1.137856e-06 7.017712e-01 8.857898e-02
    ##   [296] 1.263972e-03 1.754407e-19 5.092168e-02 4.070909e-01 1.088989e-01
    ##   [301] 3.953800e-04 9.441545e-02 8.893069e-02 6.645940e-06          NaN
    ##   [306] 1.649451e-06 3.721826e-03          NaN 4.814582e-01 6.406751e-01
    ##   [311] 1.190782e-07 6.810528e-01 1.682042e-03 5.913631e-02 1.493977e-04
    ##   [316] 1.331727e-03 6.291055e-04 7.426830e-01 5.777227e-02 1.633914e-03
    ##   [321] 9.853798e-04 4.477175e-02 9.520614e-32 2.995148e-01 5.938592e-05
    ##   [326] 3.096900e-01 7.139602e-01 7.884845e-06 6.697609e-01 2.209877e-01
    ##   [331] 6.259011e-03 2.130729e-10 8.143592e-01 4.260555e-01 2.686072e-05
    ##   [336] 9.875379e-01 4.756598e-01 3.253228e-05 3.761652e-01 6.628748e-01
    ##   [341] 1.449498e-01 3.450872e-02 2.926501e-06 5.025538e-03 7.263001e-01
    ##   [346] 2.636758e-01 7.189870e-02 9.695149e-05 1.634644e-01 1.067829e-19
    ##   [351] 7.839471e-01 9.176040e-07 4.384484e-02 3.903748e-01 3.562395e-01
    ##   [356] 4.803197e-01 6.453763e-01 1.866602e-01 3.623274e-03 3.223441e-02
    ##   [361] 2.474170e-01          NaN 4.925050e-03 3.152435e-14 5.418270e-01
    ##   [366] 2.153737e-11 5.412349e-01 3.980538e-01 1.704798e-01 8.891430e-07
    ##   [371] 2.134667e-04 3.423933e-09 1.568965e-03 2.123639e-01          NaN
    ##   [376] 5.131680e-02 6.895678e-10 6.571254e-12 1.188871e-01 6.809000e-15
    ##   [381] 8.806522e-01 9.252588e-01 1.825502e-01 3.441366e-01 3.962848e-01
    ##   [386] 8.264622e-23 2.059118e-01 1.384850e-02 4.262269e-01 8.643987e-01
    ##   [391] 6.936253e-03 2.865919e-02 3.109177e-06 1.667759e-01 2.912064e-03
    ##   [396] 1.002831e-01 1.912531e-03 9.317964e-01 2.074259e-03 2.231710e-01
    ##   [401] 7.614123e-01 3.063983e-02 3.039569e-01          NaN 1.054137e-15
    ##   [406] 3.662958e-01 1.092958e-01 5.005915e-06 9.878859e-03 5.079054e-01
    ##   [411] 6.150683e-01 3.304048e-01 7.965391e-14 3.180931e-02 2.827599e-63
    ##   [416] 5.430519e-06 3.094403e-02 1.184934e-01 2.463685e-59 9.147100e-09
    ##   [421] 8.828222e-03 1.689794e-01 1.820634e-03 1.111215e-02 9.399875e-01
    ##   [426] 6.986485e-01 1.074015e-01 1.055127e-02 8.138732e-03 2.359507e-01
    ##   [431] 1.585416e-06 5.934564e-01 8.058888e-01 5.085459e-02 7.698776e-02
    ##   [436] 4.248097e-01 4.483661e-18 5.632158e-01 4.724122e-01 4.446770e-04
    ##   [441] 9.382756e-01          NaN 2.389696e-05 7.210189e-01 5.949623e-06
    ##   [446] 1.801388e-01 1.204499e-01 7.688752e-01 3.023945e-02 8.435849e-04
    ##   [451] 5.175567e-02 7.713887e-06 1.604943e-03 3.121599e-05 5.888550e-03
    ##   [456] 3.304486e-19 5.755685e-04 9.431525e-04 8.773071e-01 1.149670e-01
    ##   [461] 7.322585e-01 7.564394e-01 9.774709e-04 7.057792e-01 5.792529e-05
    ##   [466] 5.062664e-03 7.717383e-06 3.129582e-02 1.443607e-01 3.378315e-10
    ##   [471] 2.336213e-08 1.184877e-02 3.785470e-04 2.062232e-01 7.788328e-01
    ##   [476] 8.340885e-01 2.686665e-05 1.632801e-04 4.264599e-01 1.152598e-01
    ##   [481] 5.541959e-12 6.412872e-02 1.113748e-02 3.215695e-03 8.763632e-03
    ##   [486] 2.132538e-01 2.241073e-02 7.738961e-16 2.050701e-02 9.738514e-01
    ##   [491] 4.471164e-01 1.199198e-05 4.142937e-01 7.427578e-01 1.745397e-01
    ##   [496] 1.258725e-01 1.414338e-28 5.219736e-01 7.097136e-23 3.653123e-01
    ##   [501] 6.221341e-01 1.629980e-02 2.161721e-01 5.790151e-06 4.510654e-01
    ##   [506] 2.882804e-06 1.887891e-02 7.038327e-01 3.725984e-01 1.048041e-09
    ##   [511] 3.559332e-01 1.569798e-02 1.755915e-01 6.068163e-09 2.048312e-01
    ##   [516] 3.778267e-21 1.584154e-02 4.693637e-05 4.287844e-02 3.631451e-01
    ##   [521] 5.647257e-01 3.936127e-04 1.930065e-01 7.772929e-01          NaN
    ##   [526] 1.067346e-01 2.483127e-01 4.918310e-01 3.654192e-02 2.018691e-02
    ##   [531] 1.809993e-03 9.853306e-01 8.538007e-06 2.062464e-02 1.710548e-03
    ##   [536] 2.669856e-06 7.761616e-01 8.807439e-01 2.503422e-03 1.139915e-04
    ##   [541] 5.491035e-01 4.163606e-03 5.959570e-02 8.008794e-01 8.924806e-02
    ##   [546] 8.478073e-08 4.563726e-02 3.032715e-02 2.693302e-08 7.456607e-01
    ##   [551] 1.340303e-07 1.331054e-25 1.786567e-10 4.966992e-07 8.122489e-01
    ##   [556] 5.211633e-02 3.204838e-01 2.597896e-07 9.555922e-04 3.879538e-07
    ##   [561] 9.335455e-02 2.397657e-01 5.679485e-01 1.078029e-08 1.118939e-05
    ##   [566] 8.740339e-01 9.212035e-07 2.694905e-07 1.707528e-01 5.844444e-03
    ##   [571] 5.003962e-09 1.713519e-01 2.028795e-01 5.574320e-03 1.111333e-05
    ##   [576] 1.984786e-01 2.140794e-13 3.733571e-02 6.594436e-01 2.622692e-15
    ##   [581] 8.923594e-01 1.958656e-03 3.846527e-01 7.188025e-02 1.344590e-05
    ##   [586] 1.620048e-02 2.823029e-12 1.309952e-01 1.967409e-04 9.095696e-47
    ##   [591] 3.841698e-01 5.428817e-02 5.183269e-02 4.219248e-01 2.276767e-06
    ##   [596] 5.314495e-01 3.961717e-01 2.381797e-02 2.870207e-01 4.336246e-02
    ##   [601] 2.313013e-01 6.577629e-02 2.169913e-01 6.630649e-02 1.495759e-02
    ##   [606] 2.763775e-02 9.968054e-02 3.781178e-02 5.244355e-02 1.887607e-01
    ##   [611] 3.425485e-03 1.444398e-08 2.641428e-01 4.067985e-02 9.346081e-02
    ##   [616] 1.320837e-01 7.885365e-01 6.247631e-01 1.534528e-03 6.406772e-05
    ##   [621] 6.063540e-01 3.404732e-01 7.868335e-04 7.538224e-02 1.165736e-01
    ##   [626] 3.673866e-02 1.072987e-19 2.741892e-01 1.199843e-02 8.462704e-06
    ##   [631] 1.512701e-05 1.260772e-13 2.332874e-20 4.793612e-03 1.398407e-05
    ##   [636] 8.660314e-01 1.000060e-02 5.225654e-01 7.366786e-03 4.544670e-01
    ##   [641] 1.347810e-01 3.204687e-01 1.485681e-05 1.160633e-03 1.259869e-01
    ##   [646] 2.549333e-03 8.497259e-01 4.284773e-01 5.432020e-01 5.658789e-01
    ##   [651] 1.568131e-01 3.200254e-03 6.264866e-01 5.702846e-01 4.229943e-14
    ##   [656] 8.707590e-02 2.154642e-03 2.398724e-18 5.978079e-04 2.747791e-02
    ##   [661] 1.878043e-02 2.912910e-01 3.586647e-02 3.318852e-08 7.919942e-01
    ##   [666] 5.884874e-03 3.204838e-01 3.488479e-01 3.593900e-01 8.792506e-22
    ##   [671] 5.203944e-01 6.810853e-02 9.594793e-02 1.734393e-01 3.766981e-04
    ##   [676] 8.473806e-01 9.429754e-02 1.401308e-01 4.861059e-02 2.295595e-06
    ##   [681] 2.272634e-03 6.385616e-01 2.170708e-02 8.767049e-03 8.090861e-01
    ##   [686] 1.160460e-03 7.429152e-01 4.349276e-02 1.735385e-01 2.397707e-01
    ##   [691] 6.876748e-01 4.974517e-01 3.362648e-02 8.637714e-02 9.024678e-01
    ##   [696] 3.558405e-01 3.019645e-24 2.280959e-01 1.741080e-05 2.297528e-01
    ##   [701] 1.569154e-01 1.001552e-01 1.252444e-02 1.275481e-01 8.367601e-01
    ##   [706] 1.753394e-03 3.150457e-02 3.452286e-32 2.730604e-01 1.745141e-01
    ##   [711]          NaN 6.889499e-01 6.348511e-04 1.863199e-03 2.440809e-01
    ##   [716] 1.219009e-02 3.920467e-02 6.324735e-01 1.260443e-01 9.178491e-01
    ##   [721] 8.784474e-02 2.231630e-02 2.489216e-01 2.274054e-03 5.334489e-01
    ##   [726] 2.705894e-01 4.796060e-01 2.831716e-01 3.204838e-01          NaN
    ##   [731] 1.922183e-08 6.728197e-02 1.131694e-10 1.699535e-12 1.629044e-21
    ##   [736] 3.204838e-01 4.416443e-01 7.540740e-01 6.825391e-07 2.765304e-01
    ##   [741] 3.927592e-01 5.900148e-06 1.593879e-01 4.484292e-06 2.162119e-05
    ##   [746] 1.418814e-03 8.728795e-02 2.083948e-02 1.531504e-01 2.897559e-02
    ##   [751]          NaN          NaN 7.329948e-03 3.204838e-01 6.997272e-01
    ##   [756]          NaN 2.323866e-01 5.310342e-01 5.007240e-01 2.196544e-01
    ##   [761] 8.544824e-01 4.280209e-01          NaN          NaN          NaN
    ##   [766] 9.658259e-05 8.953313e-01 5.200617e-01 3.998278e-01 5.031783e-03
    ##   [771] 7.474720e-01          NaN 2.885934e-01 4.854472e-01 1.853461e-01
    ##   [776] 3.187318e-01 1.766765e-01 1.922121e-16 8.930907e-01 5.873008e-01
    ##   [781] 1.727448e-02 3.012741e-01 1.478446e-02 7.621686e-01 6.646172e-01
    ##   [786] 8.167490e-04 3.669691e-01 1.171643e-03 7.495298e-06 4.552191e-01
    ##   [791] 1.925953e-01 5.558227e-11 6.013694e-02 4.469098e-06 1.409248e-01
    ##   [796] 3.463976e-01 3.058438e-01 6.807339e-01 5.838826e-13 1.545362e-02
    ##   [801] 1.527464e-01 1.610493e-01 1.993039e-01 3.289554e-04 7.492671e-05
    ##   [806] 1.378694e-01 6.565835e-02 2.510021e-02 8.805118e-02 1.593127e-04
    ##   [811] 3.233451e-01 1.412842e-01 9.354760e-01 4.047813e-02 9.737194e-05
    ##   [816] 8.908638e-08 7.237640e-06 7.758508e-01 4.627601e-02 6.686168e-01
    ##   [821] 6.812388e-27 6.927897e-01 5.742222e-01 6.186119e-02 1.669907e-01
    ##   [826] 7.206186e-02 7.812957e-01 1.654602e-01 8.895810e-01 4.926202e-03
    ##   [831] 9.148003e-01 1.029851e-02 4.897965e-05 3.280165e-01 1.313211e-01
    ##   [836] 6.894779e-01 3.052327e-01 2.431185e-02 5.739915e-05 7.080455e-01
    ##   [841] 6.767875e-01 6.377194e-01 1.524796e-01 3.873528e-03 2.671197e-01
    ##   [846] 1.175146e-10 1.726364e-02 4.416611e-01 6.820910e-03 4.055730e-01
    ##   [851] 6.917788e-01 4.186814e-03 6.237050e-02 2.838151e-01 5.296207e-03
    ##   [856] 1.744087e-01 7.538769e-01 1.341141e-05 1.356571e-09 7.365558e-05
    ##   [861] 8.586225e-01 6.835418e-11 3.683952e-03 5.184429e-06 6.679797e-05
    ##   [866] 5.661315e-04 9.602866e-01 8.681567e-11 1.269305e-01 4.889129e-36
    ##   [871] 6.581041e-01 1.197165e-01 8.162444e-03 5.651462e-12 2.661462e-01
    ##   [876] 1.769387e-04 5.697801e-08 6.399629e-03 3.584508e-02 2.620305e-06
    ##   [881] 1.309810e-13 9.724939e-01 9.555420e-04 5.999613e-01 7.154308e-01
    ##   [886] 2.105744e-02 5.619863e-01 1.289438e-02 2.185814e-01 7.517929e-01
    ##   [891] 7.808277e-01 7.185066e-02 2.322341e-12 4.247664e-01 6.975869e-06
    ##   [896] 6.419346e-01 5.980537e-01 1.986118e-09 1.752521e-02 2.068494e-04
    ##   [901] 2.199277e-09 1.931241e-03 3.809635e-03 1.038227e-02 1.178623e-06
    ##   [906] 8.496955e-01 1.477106e-01 7.269050e-01 1.198192e-01 7.287381e-01
    ##   [911] 4.480561e-02 6.776733e-04 8.299246e-01 6.828103e-02 5.425612e-04
    ##   [916] 8.592500e-01 1.793270e-01 7.270817e-02 2.857531e-01 1.396707e-01
    ##   [921] 1.315201e-56 1.589992e-02 7.782032e-01 1.146652e-28 1.784426e-01
    ##   [926] 3.657167e-01 4.175810e-08 1.227394e-05 1.748187e-01 3.535694e-17
    ##   [931] 6.525244e-01 2.573265e-02 4.037674e-02 8.259414e-03 3.204838e-01
    ##   [936] 3.204838e-01 5.105165e-29 2.682698e-14 1.802364e-01 8.555596e-08
    ##   [941] 1.249971e-01 1.392692e-01 7.400432e-01 4.353094e-02 5.720394e-04
    ##   [946] 4.834799e-01 2.753280e-02 9.047129e-01 7.203086e-19 4.860891e-03
    ##   [951] 9.448469e-01 4.741859e-03 9.870549e-01 7.344076e-01 8.265288e-01
    ##   [956] 7.731645e-06 1.330578e-06 1.285650e-01 7.608918e-02 3.099427e-49
    ##   [961] 6.381339e-02 5.933927e-01 7.123633e-01 4.338188e-03 5.219296e-01
    ##   [966] 3.091842e-01 7.715834e-01 9.148716e-01 8.581411e-03 1.145500e-03
    ##   [971] 6.228033e-03 8.721708e-09 1.029355e-08 2.918453e-01 4.947333e-01
    ##   [976] 3.332075e-01 1.060143e-01 1.295110e-15 7.010138e-06 2.568290e-06
    ##   [981] 1.407571e-04 8.045948e-03 1.609031e-01 1.016918e-09 4.031036e-01
    ##   [986] 9.001420e-04 1.069558e-02 2.019328e-04 2.095077e-03 3.672900e-02
    ##   [991] 7.789868e-01 6.485313e-26 1.027695e-03 4.331245e-01          NaN
    ##   [996] 1.857365e-13 5.973545e-01 4.878964e-04 3.778706e-01 3.009890e-01
    ##  [1001] 3.433962e-03 2.450706e-02 5.403830e-01          NaN 7.283962e-01
    ##  [1006] 4.137476e-01 3.869429e-02 4.053713e-12 6.862428e-01 7.962360e-01
    ##  [1011] 8.825296e-01 7.130842e-01 6.167910e-02 1.426835e-01 2.609336e-05
    ##  [1016] 1.057203e-02 1.801566e-01 9.022955e-01 5.228410e-06 6.427512e-05
    ##  [1021] 7.656915e-03 4.144277e-02 6.091732e-01 2.438046e-02 7.636065e-04
    ##  [1026] 4.143693e-01 1.748046e-02 8.669467e-01 8.286538e-08 5.663047e-11
    ##  [1031] 1.226986e-01 3.996300e-01 2.197960e-04 8.364876e-01 3.368058e-01
    ##  [1036] 9.416902e-07          NaN 3.901502e-02 3.116931e-05 4.152352e-03
    ##  [1041] 3.982103e-01 5.191781e-01 4.419857e-02 4.359142e-13 1.471187e-03
    ##  [1046] 2.034468e-02 6.540364e-08 4.642309e-02 3.354267e-01 7.318388e-02
    ##  [1051] 3.592259e-01 7.385258e-01 1.056075e-01 2.859088e-02 1.118977e-33
    ##  [1056] 1.955126e-03 7.417948e-01 9.549793e-06 1.262804e-02 4.923748e-04
    ##  [1061] 4.474354e-01 3.617713e-01 3.031948e-01 9.169426e-02 2.103856e-07
    ##  [1066] 7.976200e-01 4.148216e-01 1.654964e-01 3.647060e-01 8.408742e-01
    ##  [1071] 9.274977e-02 2.876622e-02 3.352460e-01 9.129875e-01 3.684985e-01
    ##  [1076] 6.908855e-01 9.377701e-01 1.323233e-01 6.922224e-01 1.061618e-01
    ##  [1081] 8.750218e-01 3.162584e-04 1.529823e-01 5.928313e-01          NaN
    ##  [1086] 5.060688e-01 1.190456e-01 1.095482e-01 1.110959e-01 4.671313e-01
    ##  [1091] 1.062252e-04 3.204838e-01 4.364296e-02 8.778173e-01 6.466254e-01
    ##  [1096] 6.070963e-30          NaN 8.154461e-19 3.232559e-03 3.579417e-02
    ##  [1101] 4.910846e-01 8.154469e-04 9.804953e-01 5.741014e-03 3.204838e-01
    ##  [1106] 7.489662e-01 4.492202e-02 8.060469e-01 7.811962e-01 8.651328e-01
    ##  [1111] 2.777464e-01          NaN 2.184669e-05 9.264634e-04 8.892554e-26
    ##  [1116] 2.228253e-01 4.332749e-05 8.010078e-05 1.211940e-03 3.358768e-01
    ##  [1121] 6.458698e-01 1.637074e-01 7.512552e-01 3.491433e-03 8.254441e-01
    ##  [1126] 2.038650e-02 2.591403e-07 4.266253e-01 3.204838e-01 1.434937e-01
    ##  [1131] 6.882658e-02 5.291281e-01 8.444374e-01 5.025244e-01 6.145317e-01
    ##  [1136] 4.649053e-03 4.169992e-01 5.802466e-08 1.766272e-01 3.601274e-04
    ##  [1141] 8.035958e-04 1.497809e-05 8.331047e-01 8.295358e-02 1.875245e-08
    ##  [1146] 3.048385e-29 2.011027e-03 3.580210e-01 5.958930e-05 8.357719e-01
    ##  [1151] 8.684533e-02 6.601767e-05 1.080279e-03 5.229842e-04 3.204838e-01
    ##  [1156] 6.487550e-28 6.376754e-01 2.984852e-01 6.520203e-01 2.137254e-01
    ##  [1161] 2.055487e-02 8.596231e-01 9.244924e-03 3.783996e-01 3.498196e-01
    ##  [1166] 4.442277e-01 3.306687e-02 3.692197e-02          NaN 9.291909e-03
    ##  [1171] 9.513677e-01 5.519070e-03 9.019267e-01 5.677899e-02 3.158099e-01
    ##  [1176] 1.627271e-01 1.549175e-04 9.466139e-01 1.610588e-03 4.015983e-01
    ##  [1181] 3.873184e-25 3.354308e-02 3.664564e-02 2.730425e-01          NaN
    ##  [1186] 3.124463e-01 4.423183e-01 2.791914e-02 1.371349e-02 1.341799e-01
    ##  [1191] 3.979120e-14 8.722149e-01 1.379594e-10 7.740254e-05 3.044936e-01
    ##  [1196] 1.446636e-04 4.229449e-05 7.111956e-12 1.066384e-02 6.748134e-01
    ##  [1201] 2.937864e-01 3.530499e-12 3.961692e-01 3.461375e-03 7.547365e-08
    ##  [1206] 7.110901e-03 8.508790e-01 4.574482e-03 5.074422e-08 1.689097e-12
    ##  [1211] 3.242249e-11 2.484974e-02 3.376296e-02 1.446461e-01 3.996352e-01
    ##  [1216] 4.071341e-01 9.613316e-02 7.163588e-02 6.157061e-02 1.077989e-02
    ##  [1221] 3.743095e-03 1.058264e-07 5.656473e-08 3.704838e-01 1.379248e-05
    ##  [1226] 1.132274e-06 6.582377e-08 9.577903e-01 6.263992e-01 3.204838e-01
    ##  [1231] 4.925954e-01 3.707976e-01 6.435041e-04 1.587150e-05 3.989409e-04
    ##  [1236] 2.882521e-02 2.869047e-12 2.919280e-01 1.137150e-01 7.457675e-01
    ##  [1241] 1.290781e-03 8.014606e-01 4.969339e-02 7.262982e-03 1.092795e-02
    ##  [1246] 8.202635e-01 2.454423e-05 8.522477e-03 3.461149e-01 8.250668e-02
    ##  [1251] 1.473006e-01 8.099853e-03 6.781694e-02 1.525234e-26 4.159003e-02
    ##  [1256] 3.399322e-05 5.137496e-03 4.146830e-01 3.417425e-01 1.551185e-02
    ##  [1261] 3.175862e-02          NaN 6.276636e-01 4.863535e-01 2.957241e-02
    ##  [1266] 6.534161e-01 9.224832e-04 1.817571e-04          NaN 5.230463e-05
    ##  [1271] 9.825003e-01 2.950443e-01 8.062818e-03 8.621013e-01 9.779869e-02
    ##  [1276] 1.032326e-02 9.510475e-01 1.720941e-04 7.124551e-01 2.904785e-01
    ##  [1281] 6.665881e-02          NaN 9.169538e-01 1.012053e-03 9.306417e-01
    ##  [1286] 1.203728e-05 9.147345e-02 1.335265e-01 4.422423e-01 1.531472e-01
    ##  [1291] 1.502250e-14 1.210780e-03 1.758171e-01 1.617848e-25 5.894892e-05
    ##  [1296] 6.299842e-03 5.208798e-10 3.986342e-02 2.100330e-02 5.677920e-02
    ##  [1301] 2.377648e-06 9.823209e-01 7.678876e-01 6.915318e-01          NaN
    ##  [1306] 2.897305e-02 2.939959e-01 4.782416e-03 2.016211e-02          NaN
    ##  [1311] 6.766563e-03 1.765622e-04 9.727165e-01 8.277838e-01 8.030096e-01
    ##  [1316] 6.895734e-01 3.204838e-01 3.119043e-01 1.661862e-01 1.235614e-03
    ##  [1321] 1.488874e-22 1.153384e-01 7.492760e-01 7.536894e-01 4.433672e-01
    ##  [1326] 1.026065e-01 7.871390e-01 3.065720e-11 1.793986e-01 5.977498e-04
    ##  [1331] 7.712677e-01 2.486399e-01 8.537507e-01 9.557207e-01 2.407743e-05
    ##  [1336] 1.516480e-01 2.182669e-01 7.134885e-04 5.597736e-02 4.280588e-07
    ##  [1341] 1.403613e-01 1.784922e-01 9.185180e-06 2.103183e-01 4.542293e-01
    ##  [1346] 5.871999e-01 4.345471e-01 1.658861e-02 6.788506e-01 8.777764e-01
    ##  [1351] 5.034477e-12 2.953310e-04 1.618179e-01 1.308843e-15 4.699747e-02
    ##  [1356] 3.099526e-02 8.318581e-01 1.906859e-07 3.978943e-18 1.274874e-01
    ##  [1361] 1.769616e-01 1.037605e-01 7.226395e-06 4.959066e-01 3.290665e-06
    ##  [1366] 2.181426e-44 1.989190e-01 5.686359e-01 1.545307e-04 1.640143e-06
    ##  [1371] 8.670648e-01 3.155184e-01 6.191830e-02 3.391805e-01 4.166014e-10
    ##  [1376] 9.090757e-01 8.909237e-01 6.810488e-03 4.992071e-02 9.206734e-01
    ##  [1381] 1.061657e-04 1.142048e-17 4.149723e-02 1.196451e-01 1.218570e-01
    ##  [1386] 1.216445e-06 2.263390e-08 9.084540e-01 1.571776e-01 1.576783e-02
    ##  [1391] 5.347208e-07 4.948729e-04 5.616781e-02 1.825786e-02 1.308604e-03
    ##  [1396] 8.646846e-12 3.317804e-05 3.978253e-10 1.721055e-01          NaN
    ##  [1401] 1.873724e-01 4.614082e-01 7.286430e-01 6.347293e-01          NaN
    ##  [1406] 1.559186e-01 1.525832e-02 4.404713e-01 1.274998e-05 1.654470e-27
    ##  [1411] 7.306867e-01 1.662694e-01 3.014338e-04 1.008313e-02 9.488132e-01
    ##  [1416] 2.356469e-04 1.959354e-02 2.436918e-02 7.764483e-02 5.465107e-01
    ##  [1421] 9.074031e-01 7.011336e-02 3.262364e-10 1.650382e-01 2.281646e-12
    ##  [1426] 1.345281e-11 1.354317e-02 4.246064e-08 3.393208e-01 8.511577e-01
    ##  [1431] 7.547020e-01 1.166057e-45 2.825045e-01 1.401485e-02 3.917903e-01
    ##  [1436] 8.992552e-01 5.103235e-01 8.054420e-02          NaN 3.722277e-01
    ##  [1441] 3.580672e-02 6.449864e-05 7.541643e-01 1.555085e-13 1.378952e-02
    ##  [1446] 7.535161e-09 3.101471e-01 6.723592e-05 5.631152e-38 8.591073e-01
    ##  [1451] 3.586705e-06 8.665260e-01 1.004862e-21 1.769490e-01 4.386073e-01
    ##  [1456]          NaN 1.371398e-11 3.195538e-05 8.847408e-01 1.642794e-01
    ##  [1461] 1.332947e-01 3.113749e-01 4.783345e-06 1.235089e-11 4.004436e-02
    ##  [1466] 4.607581e-01 9.900572e-01 7.311669e-01 1.026623e-01 2.969253e-01
    ##  [1471] 8.991891e-01 6.868285e-02 1.988520e-05 7.798088e-02 2.566637e-01
    ##  [1476] 3.740089e-01 7.225831e-05 1.273168e-22 1.298196e-01 3.871723e-01
    ##  [1481] 7.161475e-02 9.773935e-01 9.597047e-01 3.927987e-30          NaN
    ##  [1486] 3.223383e-01 2.216253e-14 1.859953e-01 9.117487e-01 1.590337e-01
    ##  [1491] 4.434122e-02 1.867847e-03 1.785051e-05 6.225024e-04 6.902559e-04
    ##  [1496] 3.549257e-02 6.646693e-01          NaN 1.965869e-07 2.840478e-18
    ##  [1501] 7.356873e-04 8.063564e-01 1.335295e-10 3.347500e-01 2.511503e-01
    ##  [1506] 6.531378e-01 2.523547e-03 6.206741e-02 2.195280e-30 3.413612e-04
    ##  [1511] 1.541711e-01 8.177178e-02 1.481383e-01 9.639822e-10 1.447150e-10
    ##  [1516] 8.897869e-02 3.080301e-02 2.339093e-11 4.384964e-01 2.096549e-02
    ##  [1521] 9.417349e-02 2.826130e-01 9.761221e-01 2.836071e-01 5.023332e-02
    ##  [1526] 9.518914e-21 4.162316e-01 4.567152e-01 2.495859e-01 2.068021e-01
    ##  [1531] 7.468599e-01 6.993467e-01 2.629521e-01 9.855066e-01 6.344889e-02
    ##  [1536] 3.997368e-01 3.099456e-01 1.318499e-03 2.716605e-01 2.262538e-02
    ##  [1541] 5.925741e-06 1.125544e-07 8.832321e-01 1.468118e-03 1.407081e-08
    ##  [1546] 5.214358e-01 9.954069e-04 2.707704e-01 4.255418e-01 7.457021e-33
    ##  [1551] 1.714753e-02 1.802185e-03 5.592096e-01 9.830471e-01          NaN
    ##  [1556] 4.758511e-03 8.412044e-01          NaN 4.670704e-05 1.643756e-01
    ##  [1561] 3.501014e-06 5.587581e-02 1.574992e-02 8.511937e-03 3.109077e-02
    ##  [1566] 4.092779e-01 1.139668e-01 1.973079e-01          NaN 7.205203e-01
    ##  [1571] 3.811160e-01 4.403927e-02 6.449819e-01 4.288903e-01 4.019299e-04
    ##  [1576] 1.025195e-04 6.241208e-01 7.322761e-03 2.597417e-14 6.328152e-01
    ##  [1581] 1.156969e-01 3.822342e-02 3.536019e-01 2.001934e-01 1.223811e-01
    ##  [1586] 4.517514e-03 8.334119e-06 1.034164e-01 3.078893e-02 3.108397e-02
    ##  [1591] 6.914139e-08 4.541141e-02 6.654684e-04          NaN 2.198318e-01
    ##  [1596] 7.355164e-03 1.082132e-02 4.051250e-01 3.266200e-01 3.159610e-01
    ##  [1601] 1.422386e-05 7.206040e-04 9.944136e-01 6.759355e-01 6.700075e-01
    ##  [1606] 4.174766e-03 4.329351e-27 6.760356e-02 1.742209e-01 1.207073e-09
    ##  [1611] 5.643181e-01 4.968035e-11 4.950508e-06 1.058327e-06 1.727715e-01
    ##  [1616] 6.221068e-04 2.983484e-04 1.591223e-10 1.988983e-01 9.212569e-01
    ##  [1621] 8.794523e-01 5.297216e-01 9.531143e-01 2.359614e-01 1.246157e-01
    ##  [1626]          NaN 4.771335e-03 2.023630e-01 1.357606e-01 4.266434e-01
    ##  [1631] 2.416184e-04 7.170010e-01 7.159228e-02 5.775986e-02 7.196146e-02
    ##  [1636]          NaN 3.182826e-07 1.291965e-14 8.589239e-07 6.644492e-01
    ##  [1641] 9.656884e-01 8.564068e-05 8.979617e-01 1.430631e-04 3.328911e-01
    ##  [1646] 2.400810e-01          NaN          NaN 2.582907e-03 2.540956e-07
    ##  [1651] 1.190935e-01 6.169703e-01 6.345948e-14 2.395106e-08 3.097047e-01
    ##  [1656] 3.204838e-01 5.658851e-02 5.586896e-02 1.306309e-01 1.920262e-14
    ##  [1661] 1.491922e-02 4.547059e-02 1.934155e-02 1.441590e-01 1.158167e-04
    ##  [1666] 7.635434e-13 3.964450e-08 2.401657e-23 3.150813e-01 1.043466e-23
    ##  [1671] 4.647046e-01 6.727033e-02 6.591177e-08 9.821435e-01 7.743663e-02
    ##  [1676]          NaN 2.735410e-05 8.878923e-43 9.017030e-01 7.516922e-01
    ##  [1681] 4.384797e-03 1.213742e-36 1.182627e-01 1.169979e-06 1.245128e-01
    ##  [1686] 5.650740e-01 1.964212e-02 6.557010e-07          NaN 1.767304e-01
    ##  [1691] 1.294179e-01          NaN 6.801072e-01 3.334131e-01 4.269568e-07
    ##  [1696] 8.587557e-02 2.158385e-01          NaN 9.867367e-01 6.177597e-04
    ##  [1701] 9.931731e-01 1.843687e-01 9.819439e-01          NaN 1.528244e-16
    ##  [1706] 7.142832e-02 9.477094e-04          NaN 5.904306e-01 2.345891e-02
    ##  [1711] 7.760187e-02 1.634219e-04 7.042455e-04 9.553827e-02 1.635382e-02
    ##  [1716] 8.215248e-01 1.668370e-01 2.211809e-51 4.570940e-01 1.216245e-01
    ##  [1721] 1.771248e-01 1.154519e-02          NaN 7.832040e-02 9.643343e-02
    ##  [1726] 6.665612e-29 8.201336e-02 1.309216e-01 5.802118e-06 9.452219e-01
    ##  [1731] 1.637034e-05 1.736993e-01 5.048516e-01 3.673931e-01 8.557324e-01
    ##  [1736] 8.525965e-03 1.342799e-18 3.693009e-02 1.634215e-01 9.349859e-03
    ##  [1741] 4.949050e-04 6.892758e-01 1.264511e-01 7.040499e-01 6.716155e-39
    ##  [1746] 3.875669e-01 3.501595e-01 4.438323e-01 1.414250e-04          NaN
    ##  [1751] 4.324499e-03 1.327224e-03 3.167013e-03 1.552475e-01 6.856032e-05
    ##  [1756] 2.445806e-01 4.414390e-06 9.375820e-01 7.124611e-01 7.143781e-02
    ##  [1761] 1.534621e-03 1.548268e-15 3.118727e-01 8.268390e-01 1.231958e-04
    ##  [1766] 3.677105e-01 4.050565e-03 5.895372e-05 4.031648e-01 1.666138e-01
    ##  [1771] 4.020862e-01 7.920485e-03 3.695425e-01 3.159833e-01 9.130100e-04
    ##  [1776] 1.363412e-04 2.014224e-01 1.925657e-02 1.266523e-02 7.688642e-03
    ##  [1781] 7.309746e-08 1.713706e-01 6.108806e-01 2.206178e-01 2.570742e-04
    ##  [1786] 9.596533e-02 1.879333e-01          NaN 9.115682e-01 1.087256e-02
    ##  [1791] 3.313317e-25 7.444385e-03 6.298231e-01 4.963599e-04 2.896704e-04
    ##  [1796] 8.216300e-27 7.747659e-02 3.192403e-03 1.093721e-02 2.233405e-01
    ##  [1801]          NaN 2.353777e-03 4.230480e-01 1.771228e-03 1.727757e-01
    ##  [1806] 9.151402e-06 3.075362e-01 1.654598e-01 2.685096e-24 3.777167e-02
    ##  [1811] 3.682281e-01 4.210529e-02 1.490347e-08 8.289670e-02 3.964904e-01
    ##  [1816] 4.214345e-04 4.682242e-02 2.096612e-02 8.221560e-01 3.416685e-02
    ##  [1821] 8.672406e-01 3.360315e-05 7.799595e-01 2.168859e-01 8.491801e-01
    ##  [1826] 1.405866e-07          NaN 1.346642e-08 2.886710e-02 6.221347e-02
    ##  [1831] 1.934346e-09 3.689095e-01 1.602366e-01 6.883546e-01 6.246110e-06
    ##  [1836] 2.252173e-02 5.134738e-01 1.573173e-08 9.630103e-01 4.935602e-04
    ##  [1841] 3.601249e-04 3.177242e-04 4.832185e-01 3.694928e-01 1.850933e-01
    ##  [1846] 2.333295e-03 7.690947e-01 1.512547e-01 9.232324e-01 1.524754e-02
    ##  [1851] 9.698821e-01 2.556841e-02 8.749908e-01 1.668251e-14 8.489769e-02
    ##  [1856] 3.206282e-02 9.773039e-01 2.652401e-03 1.611244e-24 1.369026e-01
    ##  [1861] 5.216410e-01 7.488771e-01 5.346467e-01 2.581564e-06 2.606962e-02
    ##  [1866] 6.168425e-02 5.356477e-01 1.754502e-01 5.935464e-01 8.809622e-02
    ##  [1871] 5.363271e-17 1.238557e-02 2.280510e-10 5.799938e-01 3.631881e-01
    ##  [1876] 2.460472e-01 5.156994e-01 5.304004e-03 1.710185e-02 8.016617e-04
    ##  [1881] 2.919437e-01 2.427931e-01 8.098593e-06 2.012720e-02 2.758759e-01
    ##  [1886] 1.375226e-01 1.503577e-04 3.204838e-01 6.654168e-02 2.619018e-01
    ##  [1891] 5.971769e-04 9.965044e-01 4.061982e-01          NaN 4.068446e-01
    ##  [1896] 6.904930e-02 4.838877e-01 1.146175e-36 3.378124e-01 2.038214e-02
    ##  [1901] 1.640672e-07 9.693170e-01 7.384882e-02 2.406045e-02 1.355526e-05
    ##  [1906] 6.465647e-01 3.482422e-01 2.181593e-02 3.099007e-04 7.361355e-22
    ##  [1911] 1.897790e-07 7.772836e-13 3.204838e-01 8.099975e-01 5.247533e-01
    ##  [1916] 1.349633e-04 6.600067e-02 5.836084e-02 1.218566e-04 1.420902e-23
    ##  [1921] 3.381231e-01 1.877650e-02 1.208085e-01 3.986832e-12 3.302264e-01
    ##  [1926] 4.216167e-02 5.681220e-05 1.107817e-22 2.055192e-03 7.066541e-01
    ##  [1931] 3.905518e-24          NaN 8.482620e-01 4.428909e-01 6.653330e-02
    ##  [1936]          NaN 4.383206e-02 7.416545e-04 1.663879e-03 8.197375e-04
    ##  [1941] 7.841235e-02 8.286882e-01 1.565977e-05 3.358195e-01 9.477035e-01
    ##  [1946] 4.304120e-02 3.507332e-04 1.430393e-07 4.193198e-01 4.734744e-01
    ##  [1951] 8.434378e-01 2.567153e-02 1.361203e-22 1.483079e-01 3.095021e-02
    ##  [1956] 2.224317e-01 8.781470e-01 1.544206e-11 7.331761e-06 2.448720e-01
    ##  [1961] 3.406855e-02 2.764592e-02 8.367743e-15 1.318237e-01 6.565527e-48
    ##  [1966] 1.866178e-02 2.805012e-01 2.096226e-01 2.829052e-03 4.232385e-01
    ##  [1971] 2.130762e-01 8.355320e-11 3.311181e-01 8.257628e-01 4.605226e-01
    ##  [1976] 2.740197e-01 4.972988e-09 1.765655e-02 2.753238e-01 1.332539e-04
    ##  [1981] 7.174079e-01 6.022351e-01 6.319071e-03 2.185113e-02 2.821400e-08
    ##  [1986] 5.041465e-01 2.247634e-01 1.803298e-13 5.459866e-01 6.903166e-03
    ##  [1991] 6.834492e-01 2.062556e-01 8.126731e-02 2.942346e-01 1.416580e-02
    ##  [1996] 8.876498e-04 1.373792e-14 4.406223e-01 3.006505e-01 2.453686e-03
    ##  [2001] 5.758563e-01 2.769964e-14 4.614383e-01 5.562814e-04 1.409110e-03
    ##  [2006] 6.101008e-01 6.175670e-01 2.996749e-03 2.382166e-03 6.834557e-02
    ##  [2011] 4.337937e-18 2.747422e-02 8.782737e-06 3.930297e-03 9.269746e-02
    ##  [2016] 6.542012e-07 7.487706e-01 7.430619e-01 7.311660e-01 1.876735e-01
    ##  [2021] 5.939941e-54 1.139103e-02 8.033122e-01 5.150254e-02 1.277827e-01
    ##  [2026] 6.595463e-02 4.728281e-14 1.377618e-03 4.263703e-01 5.725659e-02
    ##  [2031] 5.494204e-08 9.873881e-01 4.036002e-01 1.601415e-05 5.576142e-08
    ##  [2036] 3.329184e-01 4.264847e-01 8.643101e-03 5.546933e-01 2.352360e-04
    ##  [2041] 1.976511e-04 1.404373e-21 5.705875e-01 1.800507e-05 6.342272e-04
    ##  [2046] 1.372865e-01 1.908023e-01 2.047470e-03 2.580770e-02 6.693245e-05
    ##  [2051] 1.871297e-03 3.867416e-01 2.236864e-01 4.079390e-03 1.158425e-04
    ##  [2056] 8.312357e-01 5.316120e-06 9.366963e-04 2.895212e-11 9.770321e-01
    ##  [2061] 9.615995e-01 3.938761e-02 7.944319e-04 8.793124e-02 1.799333e-02
    ##  [2066] 6.921698e-03 8.638796e-01 1.507566e-01 2.296774e-02 9.513715e-02
    ##  [2071] 5.069903e-07 3.132303e-01 3.754354e-10 3.618396e-01 1.032663e-01
    ##  [2076] 2.945434e-07 6.876908e-01 8.263878e-01 4.181193e-01 6.060612e-01
    ##  [2081] 1.977518e-01 2.370528e-01 3.267937e-12 7.568575e-03 1.341925e-01
    ##  [2086] 8.999113e-01 1.407555e-01 2.263573e-10 3.563289e-01 3.467607e-36
    ##  [2091] 5.127809e-01 1.570410e-01 3.415713e-02 1.614646e-03 8.436877e-01
    ##  [2096] 3.389627e-01 1.611887e-10 6.360616e-02 3.187751e-01 9.775005e-12
    ##  [2101] 4.091268e-17 5.579341e-02 6.737902e-01 4.388191e-02 7.393658e-01
    ##  [2106] 1.576063e-01 2.370315e-04 2.079124e-02 3.089508e-05 5.104627e-02
    ##  [2111] 6.243239e-01 2.390103e-20 6.979611e-01 7.078211e-01 1.941415e-01
    ##  [2116] 2.146224e-01 8.854821e-05 6.906442e-05 1.882738e-04 2.723723e-03
    ##  [2121] 2.840119e-01 6.954248e-02 5.814890e-01          NaN 4.858128e-14
    ##  [2126] 3.277904e-06 4.602748e-02 2.802807e-01 8.959127e-01 7.383992e-01
    ##  [2131] 8.109990e-03 4.865365e-02 7.663666e-02 9.969379e-04 1.267071e-04
    ##  [2136] 5.750735e-01 7.020417e-02 3.566189e-01 6.437162e-07 3.118777e-15
    ##  [2141] 4.455031e-04 8.528515e-05 7.474042e-09 8.544609e-04 7.928671e-01
    ##  [2146] 1.285300e-01 3.137404e-02 6.480186e-01 4.040498e-03 1.177620e-01
    ##  [2151] 9.488232e-05 2.189333e-01 3.663360e-02 1.321361e-02 4.685239e-01
    ##  [2156] 1.235356e-01 3.549718e-01 2.453651e-03 1.145676e-01 2.533748e-01
    ##  [2161] 6.511928e-01 3.786491e-04 1.325969e-04 1.432864e-03 8.479916e-01
    ##  [2166] 2.327785e-01 3.279261e-03 6.471712e-01 1.523658e-32 1.473610e-01
    ##  [2171] 4.281029e-01 1.452370e-07 3.536368e-05 7.923025e-01 2.082050e-02
    ##  [2176] 1.142790e-03 8.814255e-01 2.543095e-18 7.623754e-09 1.761925e-04
    ##  [2181] 6.852241e-03 2.306545e-01 9.284481e-04 3.883258e-01 5.406722e-08
    ##  [2186] 9.944198e-01 2.119723e-01 9.068921e-01 4.271721e-03          NaN
    ##  [2191] 5.227905e-01 5.504676e-01 4.374968e-02 3.139365e-04 8.071406e-01
    ##  [2196] 1.157188e-04          NaN 2.985008e-01 5.656801e-11 2.533750e-06
    ##  [2201] 4.256607e-04 5.697216e-01 1.706608e-01 5.886776e-01          NaN
    ##  [2206] 3.157504e-01 9.248721e-01 4.463729e-04 4.891178e-03 1.410756e-04
    ##  [2211] 7.834263e-01 2.192265e-01 1.265403e-46 2.228805e-01 2.033578e-01
    ##  [2216] 8.896204e-01 1.472146e-01 1.359320e-02 2.783903e-04 2.013116e-01
    ##  [2221] 2.352145e-01 6.660569e-01 4.954999e-01 1.207687e-07 2.437819e-02
    ##  [2226] 3.423041e-01 2.046899e-01 5.695005e-01 1.105080e-01 5.462809e-02
    ##  [2231] 8.908184e-02 2.357656e-16 1.456619e-01 4.190738e-01 3.884373e-01
    ##  [2236] 3.204838e-01 4.782375e-02 3.088805e-02 1.823306e-04 5.808884e-01
    ##  [2241] 5.650735e-01 9.191635e-02 1.967973e-04 8.409578e-01 4.895938e-19
    ##  [2246] 9.716348e-01 5.035635e-07 2.538277e-01 1.525564e-04 1.117571e-03
    ##  [2251] 4.963594e-01 1.166346e-03 7.635737e-05 2.326406e-04 3.056441e-25
    ##  [2256] 4.844145e-01 1.252209e-08 3.176001e-07 1.087500e-04          NaN
    ##  [2261] 7.640850e-01 5.834757e-03 5.093366e-08 7.440982e-01 3.352613e-02
    ##  [2266] 2.171360e-02 1.157958e-18 1.478229e-01 5.766556e-01 4.068516e-03
    ##  [2271] 4.662518e-02 2.072077e-03          NaN          NaN 5.447528e-09
    ##  [2276] 5.519308e-01 4.687836e-01 6.577387e-07 8.324260e-16 3.098517e-05
    ##  [2281] 3.322946e-35 2.176991e-03 1.464035e-02 4.189107e-02 1.668134e-01
    ##  [2286] 3.406104e-01 5.642780e-02 9.605209e-02 1.542209e-07 1.890592e-05
    ##  [2291] 6.687198e-12 1.049459e-02 3.208722e-02 1.754911e-01 4.282926e-02
    ##  [2296]          NaN 3.271272e-01 5.226374e-07 5.181313e-01 9.938987e-02
    ##  [2301] 2.736500e-01 4.865293e-01          NaN 6.209821e-02 3.165302e-02
    ##  [2306]          NaN 2.693454e-02 9.219652e-01 7.647603e-01 3.548912e-01
    ##  [2311] 7.471546e-01 1.559512e-01 8.407505e-01 2.989977e-04          NaN
    ##  [2316] 9.642087e-01 1.559302e-01 8.384490e-04 1.734173e-20 2.751436e-14
    ##  [2321]          NaN 2.795939e-02 4.980491e-03 3.486819e-01 1.413114e-04
    ##  [2326]          NaN 9.202616e-01          NaN          NaN 9.125562e-01
    ##  [2331] 5.208300e-11 2.396475e-01 1.060432e-04 7.342685e-02 9.048313e-01
    ##  [2336] 4.021769e-01 6.387264e-01 3.305890e-05          NaN 1.608370e-01
    ##  [2341] 8.622356e-06 9.853688e-01          NaN 4.697849e-01 1.116453e-03
    ##  [2346] 8.638259e-01 5.404784e-03 1.455303e-02 1.396961e-03 6.671927e-01
    ##  [2351] 3.184605e-02 1.757311e-01 8.862572e-01 6.957590e-02 1.408155e-07
    ##  [2356] 2.358083e-01 6.217582e-03 1.626088e-01 9.248938e-03 3.967682e-01
    ##  [2361]          NaN 1.586955e-01          NaN          NaN 9.173237e-02
    ##  [2366] 5.164000e-01 5.503693e-01 5.642998e-02 6.389302e-02 2.892868e-01
    ##  [2371] 3.314741e-01 1.187313e-06 2.585251e-01 8.464006e-02 7.380651e-01
    ##  [2376] 8.971836e-01 7.199796e-02 4.868308e-01 8.899731e-01 2.045918e-02
    ##  [2381] 5.223921e-01 3.413798e-01 9.518809e-01 1.266900e-03 1.163293e-03
    ##  [2386] 8.129840e-01 7.166313e-01 6.716532e-32 1.433948e-01 9.769392e-01
    ##  [2391] 7.026043e-01 5.332789e-01 1.762310e-04 7.444501e-01 3.507324e-01
    ##  [2396] 4.683900e-01 5.100029e-05 7.784353e-01 1.865345e-04 1.000213e-01
    ##  [2401] 4.384021e-01 4.020387e-01 4.629363e-01 7.026774e-04 4.467536e-01
    ##  [2406] 1.296601e-01 2.758052e-05 1.583617e-02 7.232995e-04 2.038325e-01
    ##  [2411] 6.793025e-02 4.332261e-11 2.823594e-02          NaN 9.578347e-19
    ##  [2416] 9.190564e-04          NaN 6.242947e-08 6.084030e-01 6.841572e-02
    ##  [2421] 5.852665e-01 8.030629e-02 2.463169e-25 6.365417e-01 2.730058e-01
    ##  [2426] 9.823189e-06 2.766044e-01 2.587679e-01 6.050868e-01 4.534430e-02
    ##  [2431] 6.588694e-04 4.759969e-02 2.094053e-01 1.567535e-16 9.515144e-01
    ##  [2436] 9.215126e-03 1.169576e-04 9.611217e-01 9.964883e-01 8.272920e-01
    ##  [2441] 6.070260e-01 5.135069e-01 3.204838e-01          NaN 8.865743e-02
    ##  [2446] 3.578209e-01 8.568442e-02 4.107265e-01 4.468304e-04 3.728290e-02
    ##  [2451] 6.792294e-15 3.160537e-02 1.239486e-05 4.923103e-01 9.407443e-01
    ##  [2456] 1.985654e-02 2.439206e-04 3.122899e-02 9.701964e-01 7.623099e-04
    ##  [2461] 5.972094e-07 1.264380e-01 4.622747e-01 7.873695e-02 1.758316e-03
    ##  [2466] 1.388027e-03 5.622818e-09 5.756480e-06          NaN 3.629624e-03
    ##  [2471] 5.199339e-01 1.361904e-09 2.601708e-01 1.232074e-05 1.355446e-02
    ##  [2476] 4.188553e-02 2.623800e-06 2.216075e-01 4.640808e-01 5.663137e-01
    ##  [2481] 1.094518e-02 4.694593e-13 2.207868e-04 7.216660e-02 3.999520e-03
    ##  [2486] 1.004991e-04 3.970048e-05 9.328852e-01 6.925574e-01 5.821569e-02
    ##  [2491] 5.803775e-01 2.028354e-05 7.593317e-01 5.019551e-01 2.993885e-03
    ##  [2496] 3.103586e-02 6.468636e-01 1.010419e-09 2.287037e-01 2.190862e-02
    ##  [2501] 3.339825e-01 1.224686e-01 1.169333e-07 2.135948e-06 3.512453e-02
    ##  [2506] 9.437255e-01 1.543446e-01 3.525848e-02 6.650015e-01 8.176192e-01
    ##  [2511] 9.773699e-08 2.567181e-01 8.054263e-07 7.512605e-01 4.716984e-01
    ##  [2516] 7.043507e-01 3.587111e-02 4.382875e-02 2.456606e-01 3.204838e-01
    ##  [2521]          NaN          NaN 5.269036e-03 6.906719e-07 3.204838e-01
    ##  [2526] 1.118171e-01 2.169369e-01 3.445298e-07 7.138072e-01 1.814431e-03
    ##  [2531] 9.754951e-07 5.004215e-01 2.844140e-02 5.893008e-07 2.081509e-02
    ##  [2536] 7.827293e-02 4.370059e-01 8.549077e-01 1.876273e-10 5.749099e-01
    ##  [2541] 4.193834e-01 1.156330e-07 3.406822e-01 1.222156e-01 1.916985e-01
    ##  [2546] 5.482439e-01 3.796391e-03 1.187053e-16 2.467431e-01          NaN
    ##  [2551] 7.814737e-01 6.694141e-02 4.900047e-02 4.908715e-03 8.880082e-01
    ##  [2556] 4.868294e-02 2.727886e-01 1.531622e-04 8.539144e-01 8.054776e-04
    ##  [2561]          NaN          NaN 7.926214e-02 1.832760e-01 3.204838e-01
    ##  [2566] 7.955985e-18 9.859997e-01 8.035118e-01 6.266521e-13          NaN
    ##  [2571] 1.303264e-01 1.129899e-01 3.165297e-03 1.634453e-01 4.497072e-02
    ##  [2576] 9.079985e-03 4.885535e-01 3.204838e-01 1.364126e-02 3.849320e-02
    ##  [2581] 3.226297e-05 5.276130e-03 3.068585e-01 5.330678e-01 1.287884e-01
    ##  [2586] 4.948083e-01          NaN 8.916300e-02 1.052952e-02 6.937616e-03
    ##  [2591] 6.924328e-01 1.712802e-12 2.846085e-05 4.715172e-10 4.435864e-01
    ##  [2596] 6.608010e-02 1.967539e-03 1.330116e-02 1.295241e-07 1.768852e-01
    ##  [2601] 8.743984e-01 5.296746e-01 3.425493e-03 1.223124e-02 2.452184e-01
    ##  [2606] 1.530954e-04 2.221517e-01 2.804367e-01 7.433573e-02 2.338429e-12
    ##  [2611] 5.272240e-01 2.796205e-01 4.988815e-01 6.344708e-01 9.701655e-01
    ##  [2616] 3.255597e-07 3.101845e-05 9.085536e-02 3.107942e-02          NaN
    ##  [2621] 2.779000e-02 7.295962e-01 3.250438e-08 1.798206e-06 2.354087e-01
    ##  [2626] 3.894276e-02 1.410918e-01 5.614078e-01 5.384884e-02 2.141264e-10
    ##  [2631] 2.485426e-01 1.423070e-01 2.391412e-01 1.972772e-02 3.564613e-06
    ##  [2636] 1.550352e-05 3.033325e-02 5.421009e-01 2.201802e-05 2.729895e-07
    ##  [2641] 3.550912e-01 4.094630e-01 2.300339e-14 4.741720e-02 6.194711e-01
    ##  [2646] 2.620994e-01 9.974665e-06 1.877294e-02 9.648575e-03 1.922217e-01
    ##  [2651] 4.702957e-01 5.387905e-01 5.636814e-01 7.295281e-03 4.821257e-02
    ##  [2656] 2.028602e-03 8.491685e-03 1.332341e-44 2.468779e-01 3.648323e-01
    ##  [2661]          NaN 6.393134e-01 7.809166e-01          NaN 2.996168e-01
    ##  [2666] 3.271638e-02 9.916461e-01 4.400709e-01 4.253918e-05 1.932284e-01
    ##  [2671] 1.277096e-01 3.106605e-01 3.052299e-03          NaN 9.478398e-01
    ##  [2676] 3.531538e-01 1.040138e-09 1.503046e-02 6.557289e-01          NaN
    ##  [2681] 3.204838e-01 6.007045e-02 5.183752e-01 7.174211e-01 3.204838e-01
    ##  [2686] 3.495730e-07 1.634698e-04 9.251840e-01 3.466226e-01 1.116032e-01
    ##  [2691] 2.839987e-13 5.615449e-01 1.355344e-03 3.204838e-01 4.105432e-01
    ##  [2696] 1.003146e-02          NaN 3.698568e-01 1.623901e-01 4.873192e-02
    ##  [2701] 2.665178e-02 5.128866e-06 2.021170e-01          NaN 2.637856e-07
    ##  [2706] 2.013600e-09          NaN 6.240914e-01 3.399441e-02 1.348569e-01
    ##  [2711] 1.629913e-20 1.392368e-07 8.495712e-02 8.215774e-03 4.228526e-08
    ##  [2716] 5.446066e-04 2.355309e-01 1.892737e-01 4.542685e-01 2.836260e-01
    ##  [2721] 8.729939e-01 9.290895e-01 2.086073e-01 1.260466e-01 7.367048e-02
    ##  [2726] 1.043560e-01 3.251962e-01          NaN 2.347082e-01 2.778519e-05
    ##  [2731] 7.913588e-01 1.313188e-03 2.843488e-02 3.755572e-01 6.340855e-01
    ##  [2736] 1.288005e-08 5.822048e-01 1.891251e-01 3.835840e-02 1.353007e-02
    ##  [2741] 8.955132e-01 8.518255e-01 2.566670e-01 1.243225e-02 6.366115e-05
    ##  [2746] 5.291343e-04 1.602094e-13 9.965986e-01 1.439774e-05 1.274518e-13
    ##  [2751] 6.072352e-12 9.929741e-01 4.562287e-01 8.088831e-01 2.818408e-02
    ##  [2756] 2.643979e-50          NaN 1.438102e-01          NaN 6.444488e-08
    ##  [2761] 6.293124e-02 1.848466e-01 1.131287e-02 8.346498e-01 4.782306e-02
    ##  [2766] 8.052737e-01 1.040633e-07 1.538729e-01 2.504848e-03 1.629130e-01
    ##  [2771] 2.134806e-04          NaN 3.536165e-01 4.421458e-02 1.832399e-03
    ##  [2776] 6.117423e-01 3.835147e-05 2.408701e-22 1.377940e-02 3.632451e-01
    ##  [2781] 2.279051e-03 6.184163e-03 1.370964e-06 6.000159e-01 2.393931e-18
    ##  [2786] 1.526117e-03 7.190330e-15 6.996946e-02 8.068933e-21 3.951510e-14
    ##  [2791] 9.798714e-02 6.926355e-04 5.420780e-01 4.006626e-01 8.046744e-01
    ##  [2796] 1.502937e-01 4.698870e-12 4.078873e-01 7.369101e-02 5.021159e-01
    ##  [2801] 4.148672e-04 3.852397e-01 3.028742e-01 1.100113e-04 3.737411e-06
    ##  [2806] 1.867522e-02 1.676370e-05 1.853463e-02 4.674701e-02 3.553538e-02
    ##  [2811] 1.377824e-01 1.213291e-02 3.232294e-04 9.571369e-04 4.091443e-01
    ##  [2816] 9.399404e-01 2.546500e-02 1.535673e-01 1.285924e-03 5.721159e-07
    ##  [2821] 1.541922e-05 4.166473e-13 4.284270e-01 5.909630e-01          NaN
    ##  [2826] 3.436767e-01 5.621936e-01 6.511589e-05 3.115659e-02 1.237309e-04
    ##  [2831] 3.107505e-03 9.769932e-02 3.479456e-02 4.475458e-01 9.762852e-12
    ##  [2836] 9.804528e-02 7.087706e-02 2.839112e-02 9.306790e-01 1.313682e-07
    ##  [2841] 6.779778e-14 1.907882e-18 1.562381e-01 4.882599e-19 3.728011e-01
    ##  [2846] 3.204838e-01 3.204838e-01 8.692389e-01 6.963754e-02 6.269698e-01
    ##  [2851] 9.062380e-06 4.110401e-06 2.108912e-01 5.315239e-02 1.136161e-01
    ##  [2856] 1.746483e-01 1.527165e-06 9.123265e-01 8.527452e-01 5.586678e-01
    ##  [2861] 1.998209e-01 8.539690e-03 4.114438e-02 3.204838e-01 2.674347e-01
    ##  [2866] 8.878757e-01 1.319993e-03 2.667949e-06 7.920801e-04 4.015269e-25
    ##  [2871] 3.416646e-05 1.771791e-04 3.204838e-01 3.665102e-06 9.499970e-02
    ##  [2876] 3.336246e-01 1.139080e-04 7.654297e-02 6.477332e-01 1.533219e-05
    ##  [2881] 2.086115e-13 9.253332e-01 1.091877e-02 7.486346e-05 8.228609e-04
    ##  [2886]          NaN 1.330033e-02 1.250913e-06 1.205502e-01 9.099331e-01
    ##  [2891] 3.812182e-01 5.700917e-01 3.102448e-01 4.244924e-01 6.025495e-01
    ##  [2896] 6.618966e-08 2.242912e-01 1.655245e-01 6.937831e-01 2.600522e-03
    ##  [2901] 1.742682e-01 8.246933e-06 1.162898e-13 4.100236e-02          NaN
    ##  [2906] 2.364593e-05 4.861053e-05 3.576366e-10 6.816992e-01 6.589063e-01
    ##  [2911] 2.805729e-01 2.901068e-04 9.383846e-01 1.763533e-01 1.496678e-09
    ##  [2916] 7.647031e-01 4.136069e-01 5.068809e-02 9.420741e-01 2.459490e-06
    ##  [2921] 1.210383e-15 3.590253e-04 7.631426e-01 1.174086e-01 5.301037e-08
    ##  [2926] 2.320598e-01 5.790558e-02 9.487718e-02 5.273527e-04          NaN
    ##  [2931] 9.573790e-10 2.393339e-01 1.657801e-01 3.846968e-01 4.507263e-01
    ##  [2936] 5.421864e-02 1.000150e-01 9.604920e-02 2.523165e-12 2.443913e-14
    ##  [2941] 1.666843e-03 7.369197e-01 7.025669e-01 8.302785e-02 5.106714e-01
    ##  [2946] 1.352730e-01          NaN 4.509236e-07 1.651938e-05 6.886528e-05
    ##  [2951] 7.789065e-01 3.930412e-01 6.871435e-08 1.735722e-01 3.204838e-01
    ##  [2956] 3.530866e-04          NaN 3.798033e-02          NaN 9.070473e-09
    ##  [2961]          NaN 7.626228e-04 4.801040e-01 5.877992e-01          NaN
    ##  [2966] 4.877834e-01 2.946695e-01 5.433919e-01 7.355935e-02 3.493598e-05
    ##  [2971] 4.504591e-01 1.192933e-03 8.204008e-01 4.286328e-01 1.438520e-01
    ##  [2976] 2.391273e-02 2.599972e-01 9.323955e-01 7.497440e-07 1.256869e-02
    ##  [2981] 2.227723e-01 5.975571e-01 1.497462e-01 1.494595e-04 6.591792e-01
    ##  [2986] 3.573012e-01 1.424501e-03 7.255567e-02 1.050410e-01 2.540503e-12
    ##  [2991] 8.577027e-18 1.448405e-02 1.966389e-03 4.788132e-01 7.213438e-01
    ##  [2996] 1.483925e-02 5.797966e-01 1.511861e-01 9.435404e-01 7.992292e-01
    ##  [3001] 2.743306e-01 1.218522e-02 2.397419e-03 6.252614e-03 2.126867e-02
    ##  [3006] 4.080282e-01 3.497397e-01 9.854015e-01 3.898992e-01 8.849402e-02
    ##  [3011]          NaN 3.204838e-01 8.772743e-02 3.073847e-02 4.412026e-03
    ##  [3016]          NaN 2.318077e-04 2.201126e-03 8.199752e-03 3.559725e-01
    ##  [3021] 8.546472e-09 8.356275e-01 4.434461e-01 5.595542e-01 9.772458e-03
    ##  [3026] 3.910300e-16 4.400387e-01 1.530582e-02 3.527839e-03 5.916004e-01
    ##  [3031] 6.339278e-01 1.050097e-14          NaN 4.034113e-01 8.015657e-02
    ##  [3036] 4.828183e-01 1.433788e-01          NaN 5.900717e-01 6.254081e-01
    ##  [3041] 5.905565e-01 2.997269e-05 1.993565e-01 7.072045e-03 2.404884e-01
    ##  [3046] 4.627884e-01 5.455884e-01 1.093518e-01 4.583216e-01          NaN
    ##  [3051] 3.239635e-01 6.190270e-01 6.291078e-07 4.743043e-01 1.572970e-01
    ##  [3056] 3.847454e-03 8.731277e-01 3.086291e-09 1.298189e-01 4.430768e-08
    ##  [3061] 2.034733e-03 6.276237e-01 5.023167e-02 4.271076e-02 6.872341e-01
    ##  [3066] 3.204838e-01 3.361427e-04 1.130826e-01 8.666637e-01 3.912816e-02
    ##  [3071] 8.614393e-01 3.009416e-01 8.191511e-01 3.585661e-08 2.391864e-02
    ##  [3076] 1.088617e-02 8.728050e-01 2.984003e-01 9.698386e-02          NaN
    ##  [3081] 1.677486e-01 3.207052e-01 1.428375e-02 8.729344e-02 2.223538e-01
    ##  [3086] 5.319480e-08          NaN 1.262382e-01 2.966770e-01 3.518987e-01
    ##  [3091] 8.298979e-01 2.583053e-01 1.207739e-01 2.752832e-02 4.821741e-02
    ##  [3096] 2.706817e-02 1.093545e-15 5.367491e-01 4.102565e-04          NaN
    ##  [3101] 8.087796e-03 3.498139e-06 7.210975e-03 7.291808e-04 2.014385e-01
    ##  [3106] 3.820372e-01 3.235510e-05 1.440125e-03 9.901062e-01 1.354402e-05
    ##  [3111] 1.170246e-07 1.234851e-03 1.048159e-01 4.215431e-01 5.242406e-01
    ##  [3116] 8.982558e-08 6.123622e-01 6.103238e-01 7.825559e-02 4.829665e-04
    ##  [3121] 1.081599e-08 1.111202e-01 5.791211e-02 1.992031e-03 3.569622e-02
    ##  [3126] 8.746523e-01 2.298437e-04 1.008752e-02 7.080066e-04 1.379912e-02
    ##  [3131] 3.364783e-01 2.231187e-01 3.204838e-01 1.702108e-01 1.913831e-21
    ##  [3136] 6.599878e-01 6.912072e-01 1.469366e-01 1.108190e-01 7.472988e-02
    ##  [3141] 8.975657e-01 5.764950e-01 1.982842e-05 8.177157e-01 3.975726e-16
    ##  [3146] 3.082459e-01 7.221353e-01 5.331920e-03 6.195725e-03 3.911571e-01
    ##  [3151] 8.542803e-01          NaN 1.739825e-02 1.166985e-07 1.177162e-01
    ##  [3156] 5.070115e-04 5.447132e-01 2.846693e-03 9.736641e-01 3.953498e-08
    ##  [3161] 2.823896e-01          NaN 3.567192e-03 4.945584e-01 1.115454e-01
    ##  [3166] 2.788806e-01          NaN 8.905017e-03 1.748641e-03 1.269705e-01
    ##  [3171] 4.218531e-01 1.866317e-04 5.435494e-02 2.891755e-02 2.479176e-06
    ##  [3176] 7.709673e-04 5.481605e-03 4.474195e-01 6.820810e-01 6.689178e-04
    ##  [3181] 5.283580e-03 9.517927e-09 3.634864e-01 8.009765e-02 8.998649e-01
    ##  [3186] 2.390199e-01 2.007952e-04 3.351636e-01 9.593279e-04 9.070697e-01
    ##  [3191] 1.102292e-04 8.233867e-01 2.879012e-02 3.204838e-01 1.527859e-11
    ##  [3196] 7.917960e-01 1.037687e-01 2.982442e-04 1.959032e-02 5.689087e-01
    ##  [3201] 3.310097e-04 2.150643e-05 7.670336e-02 6.763530e-01 2.524643e-01
    ##  [3206] 7.424151e-09 2.928474e-02 1.168368e-01 8.529230e-26          NaN
    ##  [3211] 1.002839e-02 1.200179e-01 2.838281e-01 7.713719e-02 3.274865e-09
    ##  [3216] 2.434579e-04 2.442673e-01 2.055342e-08 3.616144e-09 6.535886e-01
    ##  [3221]          NaN 2.199818e-01 1.129800e-02 7.731106e-01 2.821133e-02
    ##  [3226] 2.267666e-02 3.120537e-01 1.044571e-04 1.337159e-02 8.707051e-03
    ##  [3231] 2.807644e-02 4.543751e-01 3.159465e-03 2.408582e-01 4.413594e-15
    ##  [3236] 1.613958e-01 5.340843e-02 2.316637e-01 1.179723e-03 3.780143e-02
    ##  [3241] 5.932536e-01 7.671704e-01 9.186390e-01 3.759703e-07 6.191345e-01
    ##  [3246] 7.834858e-03 1.418451e-01 2.873451e-01 2.999040e-01 7.819511e-01
    ##  [3251] 3.204838e-01 9.145385e-01 4.348887e-05 1.011191e-02 3.721638e-07
    ##  [3256] 9.928077e-02 1.270327e-08 1.095074e-02 2.937283e-02 4.390119e-01
    ##  [3261] 5.623028e-01 2.554459e-01 1.793527e-01 2.065082e-01 5.849710e-02
    ##  [3266] 9.758588e-01 2.873151e-01 5.077946e-02 8.499212e-02 1.201537e-08
    ##  [3271] 6.991716e-02 9.527396e-01 3.112304e-04 1.388944e-01 8.327994e-01
    ##  [3276] 9.616996e-01 9.463957e-08 7.483473e-01 1.401858e-04 3.760496e-01
    ##  [3281] 5.837668e-05 3.620332e-04 5.001213e-01 6.866192e-05 6.799909e-03
    ##  [3286] 2.116848e-01 4.075987e-01 9.985740e-01 4.031601e-02 1.828489e-01
    ##  [3291] 8.363957e-01 4.564943e-02 7.325390e-01 8.885551e-01 8.121006e-06
    ##  [3296] 2.037145e-03 3.352684e-03 7.169886e-07 7.821208e-01 1.025055e-01
    ##  [3301] 3.462990e-01 1.947416e-05 6.988568e-02 9.573320e-01 3.590131e-01
    ##  [3306] 5.531247e-01 1.121212e-01 3.398912e-02 1.664480e-02 3.588163e-01
    ##  [3311] 9.127968e-01 3.809125e-01 1.269900e-03 7.675594e-01 4.906520e-01
    ##  [3316] 5.380110e-03 7.491747e-01 5.031406e-03 7.455994e-01 9.495687e-01
    ##  [3321] 3.204838e-01          NaN 4.999455e-22 1.025521e-01 8.988736e-03
    ##  [3326] 2.031378e-04 2.231008e-01 2.882479e-01 7.629330e-01 2.840700e-02
    ##  [3331] 2.272076e-01 7.727876e-02 2.285535e-02 3.614252e-01 8.261433e-01
    ##  [3336] 4.304831e-01 5.915522e-01 1.556982e-13 4.697792e-01 4.621309e-07
    ##  [3341] 9.164820e-04 2.926343e-04 3.640244e-03 5.340703e-07 9.862993e-01
    ##  [3346] 3.656521e-01 2.463976e-04 4.637067e-04 6.330131e-07 4.237633e-02
    ##  [3351] 8.636188e-01          NaN 3.726827e-04 8.242288e-02 1.798116e-01
    ##  [3356] 5.279812e-06 1.007444e-02 1.726564e-01 8.630366e-01 6.530274e-01
    ##  [3361] 2.124293e-04 1.928007e-18 4.373853e-01 3.136621e-02 4.259828e-01
    ##  [3366] 7.113587e-03 2.692276e-01 1.250597e-04 7.047949e-03 2.687834e-02
    ##  [3371] 3.119207e-08 9.560950e-01 3.730945e-01 7.602759e-05 1.475125e-03
    ##  [3376] 2.001982e-01 2.037724e-03 3.795762e-01 5.987829e-01 3.506909e-01
    ##  [3381] 7.839653e-01 9.167905e-01 6.704147e-01 1.070078e-04 9.872491e-01
    ##  [3386]          NaN 7.294245e-02 2.481894e-01 8.063676e-01          NaN
    ##  [3391] 3.879108e-01 4.547586e-09 2.087405e-04 3.180623e-01 3.204838e-01
    ##  [3396] 8.711519e-04 9.754571e-03 1.344293e-03 7.643553e-01 6.560807e-01
    ##  [3401]          NaN 1.026209e-01 2.126380e-08 3.264084e-01 3.156913e-02
    ##  [3406] 6.322643e-02 2.324122e-09 8.024670e-02 3.982321e-17 7.854554e-01
    ##  [3411] 6.400773e-03 4.643245e-01 1.438973e-01 4.737150e-06 1.270891e-04
    ##  [3416] 2.150603e-01 4.655937e-11 6.848792e-01 7.953885e-04 6.718458e-01
    ##  [3421] 8.307912e-01 3.131808e-02 1.178671e-02 3.614054e-43 4.792457e-07
    ##  [3426] 7.036147e-01 6.272634e-06 2.911120e-03 3.251097e-17 7.909475e-01
    ##  [3431] 1.093148e-02 2.988721e-02 5.527245e-01 6.111277e-01 2.599947e-01
    ##  [3436] 6.010884e-02 2.054682e-01 4.183366e-07 7.804780e-05 1.771269e-03
    ##  [3441] 2.005717e-03 3.296928e-01 2.300490e-01 8.966412e-02 2.439098e-03
    ##  [3446] 1.046898e-03 9.415790e-01 1.057731e-02 4.628623e-01 5.261680e-01
    ##  [3451] 6.471700e-02 2.535547e-02 1.123387e-02 9.872113e-21 8.297494e-01
    ##  [3456] 8.739996e-23 5.031628e-03 1.307603e-01 3.055335e-01          NaN
    ##  [3461] 1.458501e-01 5.009460e-01 1.034666e-03 2.737850e-01 1.076582e-01
    ##  [3466] 3.011603e-01 1.270117e-01 4.794598e-28 7.629993e-36 1.631482e-01
    ##  [3471] 8.381191e-08 8.398448e-02 6.471553e-04 5.642810e-05 8.740142e-05
    ##  [3476] 2.718798e-01 5.127338e-05 7.077134e-01 1.284455e-01 3.925569e-01
    ##  [3481] 6.378232e-03 2.306886e-02 7.775368e-01 7.694202e-01 3.419627e-01
    ##  [3486] 2.907850e-01 1.033800e-03 6.975475e-04 2.270341e-01 1.292793e-02
    ##  [3491] 7.499584e-09 7.541783e-12 1.184403e-01 8.085189e-01 1.240067e-02
    ##  [3496] 7.790680e-02 1.268088e-20 8.491049e-01 7.700100e-03 1.214526e-01
    ##  [3501] 1.230925e-01 1.995612e-01 9.112860e-01 1.687198e-02 5.565575e-10
    ##  [3506] 1.271803e-01 5.232564e-01 1.645316e-01 2.208386e-12 9.462864e-01
    ##  [3511] 6.595220e-01 3.627927e-04 3.840785e-02 3.028897e-01 2.676690e-01
    ##  [3516] 1.087076e-01 5.255773e-02 8.005754e-01 3.244227e-03 4.219092e-01
    ##  [3521] 7.794626e-02 1.841520e-06 1.155887e-01 8.353205e-01 2.113633e-02
    ##  [3526] 2.267393e-08 9.507920e-01 6.662495e-01 2.273551e-02 4.324773e-01
    ##  [3531] 1.228759e-04 5.441899e-04 1.226244e-10 1.024744e-29 3.204838e-01
    ##  [3536] 4.309795e-04 3.234020e-05 3.559370e-13          NaN 2.577743e-01
    ##  [3541] 3.696193e-01 2.766092e-05 2.752312e-02 1.752862e-01 5.060247e-01
    ##  [3546] 7.728282e-01 1.115858e-02 8.908522e-01 1.708797e-01 3.407385e-02
    ##  [3551] 7.402399e-01 1.348062e-01 2.640432e-06 1.426902e-21 2.748360e-01
    ##  [3556] 6.962710e-01 4.611989e-01 2.223654e-03 1.294535e-04 3.951520e-10
    ##  [3561] 8.942238e-01 6.035620e-02 1.278185e-02 9.877412e-01 1.370812e-01
    ##  [3566] 5.505589e-05 5.701872e-12 1.911382e-03 2.993814e-04 9.431454e-03
    ##  [3571] 2.063152e-02 6.444321e-01 2.282030e-03 2.587631e-01 1.560966e-01
    ##  [3576] 9.790332e-02 6.701588e-09 2.660019e-02 7.173484e-06 4.062601e-03
    ##  [3581] 6.647512e-02 4.518655e-01 5.281167e-01 7.756137e-02 2.468563e-03
    ##  [3586] 2.405861e-01 8.992124e-03 4.848681e-01 4.034727e-10 3.094515e-05
    ##  [3591] 4.740813e-02 7.580862e-02 5.941811e-03 7.244937e-01 6.226412e-01
    ##  [3596] 9.399162e-08 6.042932e-01          NaN 5.797232e-02 1.050033e-01
    ##  [3601] 7.072657e-03 5.804869e-02 9.143256e-01 5.721291e-01 9.882670e-01
    ##  [3606] 6.683006e-01 3.576385e-01 6.545800e-01 2.768951e-01 5.052760e-05
    ##  [3611] 8.533198e-01 5.066079e-01 5.697508e-02 1.343275e-01 1.381529e-03
    ##  [3616] 8.603739e-01 2.727675e-01 2.558867e-02 3.419545e-04 5.817411e-02
    ##  [3621] 1.133071e-06 6.183723e-01 2.975816e-02 5.856616e-02 4.159616e-02
    ##  [3626] 4.286484e-01 1.535000e-03 7.781651e-03 3.306696e-16 6.396224e-02
    ##  [3631] 9.969290e-01 1.090053e-09 7.180928e-18 5.981802e-03 2.132751e-01
    ##  [3636] 5.036307e-07 1.660596e-05 1.942480e-01 7.822767e-05 7.022951e-01
    ##  [3641] 6.784397e-01 3.623168e-01 3.254630e-01 1.403066e-10 4.523465e-03
    ##  [3646] 4.251077e-04 5.081135e-07          NaN 3.125995e-01 1.615320e-02
    ##  [3651] 2.493205e-08 3.322623e-05 6.586769e-03 3.292883e-01 1.173138e-05
    ##  [3656] 7.520414e-07 2.076518e-08 9.934361e-01 4.203166e-01          NaN
    ##  [3661] 3.668131e-15 8.628540e-27 1.895362e-01 8.832619e-03 1.827366e-01
    ##  [3666] 6.280912e-01 4.408378e-01 9.862499e-01 1.776671e-02 2.250554e-01
    ##  [3671] 1.063666e-18 3.793662e-01 1.611545e-01 1.290262e-01 1.776263e-19
    ##  [3676] 1.259569e-01 3.168429e-01 1.703095e-02 2.950297e-01 5.689171e-03
    ##  [3681] 5.721105e-09 7.778615e-01          NaN 1.025224e-01 6.116280e-03
    ##  [3686] 1.929387e-06          NaN 1.373626e-01 4.919355e-02 8.712099e-01
    ##  [3691] 2.378551e-01 4.432762e-06 5.278344e-01 1.527668e-03 3.362762e-01
    ##  [3696] 7.911020e-01 1.716590e-04 2.940452e-03          NaN 1.627038e-02
    ##  [3701] 6.372930e-01 2.554178e-05 8.813713e-03 9.732924e-01 6.579997e-03
    ##  [3706] 3.843523e-03 6.486300e-01 9.935844e-01 4.303030e-03 5.539793e-02
    ##  [3711] 3.456321e-02 2.097568e-02 2.130654e-01 1.670698e-01 1.585455e-15
    ##  [3716] 9.032367e-01 2.312835e-01 1.040094e-03 5.172858e-01 5.644083e-18
    ##  [3721] 7.145343e-01 8.891329e-03 1.111467e-03 1.291405e-01 6.248639e-03
    ##  [3726] 4.390954e-02 5.413918e-02 7.874017e-31 3.945180e-02 5.104253e-01
    ##  [3731] 1.218680e-09 9.352248e-33 8.966230e-08 6.708503e-01 1.603168e-23
    ##  [3736] 6.522681e-02 2.434892e-06 1.191381e-01 2.441861e-04 5.911916e-01
    ##  [3741] 4.844626e-09 1.699135e-02 1.154281e-01 5.907236e-01 9.470472e-01
    ##  [3746] 4.478125e-01 3.513509e-03 9.545956e-03 2.879467e-02 1.035454e-01
    ##  [3751] 6.072517e-02 8.415586e-08 1.129186e-13 3.414277e-02 7.586063e-02
    ##  [3756] 2.985336e-01 1.698014e-03 4.652896e-01 6.366034e-06 4.647858e-01
    ##  [3761] 1.365989e-09 1.269376e-02 1.136188e-02 6.055858e-01 2.525499e-01
    ##  [3766] 9.042398e-15 7.895681e-12 1.255962e-02 1.505929e-01          NaN
    ##  [3771] 5.244318e-02 1.189390e-01 8.724685e-01 2.943671e-01 5.201130e-01
    ##  [3776]          NaN 6.920837e-01          NaN 1.155617e-22 9.314591e-01
    ##  [3781] 8.187497e-01 9.768098e-01 1.989163e-01 9.079234e-01 7.241355e-02
    ##  [3786] 2.982112e-01 1.736550e-01 1.390416e-14 4.808597e-03 1.952875e-01
    ##  [3791] 1.211464e-01 5.443992e-01 1.770202e-08 1.261384e-03 9.029080e-02
    ##  [3796] 4.569338e-05 4.322510e-01 5.125852e-03 5.431774e-01 9.825623e-08
    ##  [3801] 9.317200e-01 6.666447e-01 1.427614e-01 8.448995e-01 1.136669e-06
    ##  [3806] 8.973471e-01 7.619871e-35 5.407598e-07 3.327033e-01 1.746313e-01
    ##  [3811]          NaN 3.201175e-13 1.470399e-39 1.329803e-04 3.067946e-09
    ##  [3816] 1.416589e-03 8.628900e-01 1.633674e-15 3.204838e-01 1.269489e-01
    ##  [3821] 4.266011e-01 6.520949e-04 6.550621e-01 9.129951e-03          NaN
    ##  [3826] 8.346945e-02 2.095149e-02 5.150832e-05 3.764965e-15 7.353348e-01
    ##  [3831] 1.569339e-02 5.907020e-01 4.120972e-01 3.623847e-01 1.100775e-02
    ##  [3836] 6.341008e-02          NaN 1.352681e-01 4.685463e-02 2.812946e-03
    ##  [3841] 2.330295e-03 4.585255e-01 2.564440e-01 3.906518e-07 5.594911e-01
    ##  [3846] 8.413566e-03 1.024658e-01 3.384653e-02 3.110138e-01 3.403583e-03
    ##  [3851] 2.177467e-04 2.998797e-02 3.253616e-01 3.308304e-03 2.173827e-02
    ##  [3856] 2.583859e-01 1.897258e-01 9.705775e-01 9.648956e-02 2.999991e-08
    ##  [3861] 1.427055e-21 4.010253e-08 1.307566e-01 5.943368e-01 8.201919e-01
    ##  [3866] 7.387896e-01 1.895773e-01 6.359202e-03 6.042423e-02 8.641712e-01
    ##  [3871] 7.498191e-01 7.894379e-02 1.978616e-03 7.271441e-01 1.276577e-01
    ##  [3876] 3.373530e-01          NaN 2.037235e-10 6.830999e-01 5.608157e-11
    ##  [3881] 9.723862e-01 7.862003e-01 3.831724e-01 2.569607e-03 5.825255e-03
    ##  [3886] 5.082619e-01 3.379071e-01 3.204838e-01 1.129625e-18 4.846853e-01
    ##  [3891]          NaN 6.203791e-01 9.133370e-25 4.632498e-01 5.787595e-03
    ##  [3896] 2.787405e-01 1.164721e-04 3.816871e-02 8.428246e-01 6.024108e-03
    ##  [3901] 2.088610e-02 5.354141e-01 6.288434e-17 7.731286e-01 3.204838e-01
    ##  [3906] 3.670707e-06 4.093228e-03 2.287873e-06 6.810488e-02 7.138863e-01
    ##  [3911] 8.100042e-01 9.530621e-05 5.208960e-02 5.872609e-09 3.612399e-07
    ##  [3916] 5.134882e-01 1.288340e-01 2.939755e-01 9.942387e-01 3.777423e-01
    ##  [3921] 1.769723e-06          NaN 6.406204e-01 1.379276e-15          NaN
    ##  [3926] 4.645978e-03 8.635915e-10 7.574623e-02 2.796910e-06 6.453366e-01
    ##  [3931] 7.562831e-01 5.547934e-01 4.372143e-11 1.523846e-01 4.734812e-01
    ##  [3936] 9.428716e-01 8.385909e-02 2.307447e-02 8.199027e-01 7.315395e-17
    ##  [3941] 5.203071e-01 1.427441e-02 5.276287e-01 1.149114e-01 2.339895e-10
    ##  [3946] 3.948533e-01          NaN 1.935850e-02 5.088921e-01 2.108770e-03
    ##  [3951] 1.928019e-06 3.466212e-02 7.258661e-10 3.664881e-01 3.702831e-01
    ##  [3956] 7.869683e-03          NaN          NaN 4.655466e-03 3.202165e-02
    ##  [3961] 3.512102e-01 3.712909e-01          NaN 9.477982e-02 8.227096e-01
    ##  [3966] 1.965940e-02 9.951241e-02 3.070452e-01 3.249237e-01 2.953834e-01
    ##  [3971] 8.989059e-02          NaN 8.919540e-01 1.323539e-02 1.940508e-05
    ##  [3976] 3.877166e-08 9.057128e-38 2.228177e-02 1.834338e-01 5.301911e-01
    ##  [3981] 8.957487e-01 8.557588e-01 3.849570e-19 1.485284e-01 3.302766e-01
    ##  [3986] 9.222559e-03 1.959989e-05 3.254286e-01 8.065417e-01 1.559026e-01
    ##  [3991] 1.844613e-01 9.102217e-01 7.567485e-03 9.412641e-01 9.367712e-01
    ##  [3996] 2.175513e-55 5.300828e-02 8.234648e-03 5.624769e-03 2.307739e-01
    ##  [4001] 3.616760e-01 2.297172e-04 1.104250e-05 5.692752e-06 4.127314e-01
    ##  [4006] 3.124552e-01 9.126893e-04 3.583584e-02 9.426696e-01 5.453792e-06
    ##  [4011] 9.416349e-01 7.428946e-02 1.079731e-04 1.103198e-01 6.488521e-01
    ##  [4016] 2.102480e-01 1.082524e-02 4.545333e-01 4.892355e-02 7.548503e-02
    ##  [4021] 3.607557e-01 9.201638e-01 6.533707e-05 7.349900e-06 1.318791e-03
    ##  [4026] 7.282243e-04 1.129172e-01 5.951074e-04 3.069460e-02 2.190587e-01
    ##  [4031] 5.167560e-09 1.126820e-04 1.262913e-11 3.277510e-05 2.242649e-10
    ##  [4036] 1.969868e-03 9.167262e-01 1.776547e-18 4.441928e-20 2.251359e-03
    ##  [4041] 1.527004e-12 9.393271e-01 3.204838e-01 7.154363e-03 3.250285e-03
    ##  [4046] 6.595365e-01 3.823929e-01 3.406567e-01 6.497776e-04 3.030130e-05
    ##  [4051] 3.489187e-01          NaN 3.943290e-01 3.204838e-01 7.409283e-08
    ##  [4056] 2.098146e-04 9.889918e-01 6.487713e-03 3.817555e-01 9.623522e-03
    ##  [4061] 4.867133e-01 9.077917e-01 4.252849e-01 4.907139e-03 3.445084e-13
    ##  [4066] 3.523848e-01 3.678888e-02 4.874774e-01 6.716317e-05 1.601989e-04
    ##  [4071] 3.868597e-03 1.528183e-13 9.221991e-01 9.497760e-01 1.375992e-04
    ##  [4076] 3.253113e-01 1.389922e-07 2.838393e-04 3.188706e-01 1.410726e-01
    ##  [4081] 1.460385e-02 3.628181e-01 1.578495e-01 3.416939e-01 3.010549e-02
    ##  [4086] 6.231927e-04 1.694819e-01 1.961733e-01 3.535620e-03 2.371889e-01
    ##  [4091] 2.096369e-01 2.086842e-21 5.540019e-03 8.133889e-03 2.837235e-11
    ##  [4096] 1.298650e-04 5.217873e-01 7.975003e-01 3.671036e-04 4.071521e-02
    ##  [4101] 4.748320e-01 2.331131e-01 1.918533e-15 1.734155e-01 5.235286e-01
    ##  [4106] 2.463021e-01 8.371240e-02 8.698278e-01 2.186131e-01 5.156952e-05
    ##  [4111]          NaN 2.021457e-01 5.146586e-01 5.689541e-20 1.489975e-03
    ##  [4116] 2.682953e-01 4.494726e-01 3.003966e-01 1.159257e-03 2.675544e-01
    ##  [4121] 6.815695e-01 1.531685e-04 2.540341e-01 3.744374e-01 1.126779e-18
    ##  [4126] 7.287309e-02 3.446861e-01 3.103992e-03 5.001554e-01 1.844097e-17
    ##  [4131] 1.326302e-01 3.473420e-43 7.510820e-01 2.110125e-01 1.669556e-01
    ##  [4136] 6.671096e-01 7.969407e-02 4.947137e-01 9.929351e-18 1.983779e-01
    ##  [4141] 9.046079e-04 2.985717e-01 3.288925e-01 9.390711e-01 4.792718e-01
    ##  [4146] 2.867484e-03 7.893148e-03 8.260271e-01 1.463406e-01 6.296816e-01
    ##  [4151] 1.772556e-01          NaN 3.737345e-01 4.571752e-08 6.241706e-05
    ##  [4156] 1.324725e-01 2.975914e-10 2.092408e-01 2.911154e-01 2.928510e-08
    ##  [4161] 1.071429e-05 9.449571e-01 1.299626e-10 4.675775e-01 2.103465e-02
    ##  [4166] 3.990039e-02 3.204838e-01 2.802838e-01 8.509343e-02 1.984149e-02
    ##  [4171] 1.293287e-01 7.824462e-01 3.977217e-04 1.723876e-02 1.906707e-04
    ##  [4176] 8.774803e-01 6.402117e-20 9.191901e-01 1.474230e-28 8.189233e-01
    ##  [4181] 6.938697e-01 1.134039e-14 3.834877e-01 9.891731e-08 7.553190e-01
    ##  [4186] 2.228172e-05 4.633830e-01 6.657680e-04 1.064063e-11 2.489912e-04
    ##  [4191] 6.340487e-03 2.079329e-01 7.240870e-02 2.295495e-01 9.600228e-05
    ##  [4196] 4.395313e-01 1.407668e-06 8.646587e-03 1.705649e-01 3.189194e-02
    ##  [4201]          NaN 6.098426e-01 1.685549e-01 3.880351e-15 1.083975e-08
    ##  [4206] 2.412493e-01 9.031072e-01 2.767913e-01 4.837642e-07 1.339174e-12
    ##  [4211] 6.633702e-01 2.015046e-03 5.408002e-05 8.200446e-20 4.962930e-01
    ##  [4216] 1.847077e-02          NaN 1.085248e-01 5.089596e-01 2.602140e-01
    ##  [4221] 1.764836e-02 5.037055e-01 8.551511e-12 1.164610e-02 3.494360e-01
    ##  [4226] 9.351954e-01 6.667258e-01 3.146111e-01 3.011372e-03 6.776404e-01
    ##  [4231] 8.427581e-01 9.254615e-01 1.501752e-11 8.384734e-02 3.050795e-04
    ##  [4236] 3.934764e-01 7.429156e-02 5.749468e-01 9.871562e-01 1.514997e-03
    ##  [4241] 7.117376e-03 1.261295e-01 2.797159e-01 4.780034e-08 5.091435e-08
    ##  [4246] 7.064314e-01          NaN 2.557269e-04 1.973229e-08 2.441385e-17
    ##  [4251] 1.978939e-07 2.140900e-02 4.370933e-05 2.970104e-32 1.702072e-04
    ##  [4256] 4.009791e-04 2.496531e-01 5.182121e-01 1.638238e-08 9.035480e-05
    ##  [4261] 1.386274e-13 1.497572e-04 4.857694e-01 2.175738e-13 3.225684e-15
    ##  [4266] 3.404690e-01 1.130602e-07 4.818074e-02 6.047646e-01 7.117050e-02
    ##  [4271] 3.204838e-01 1.089039e-03 2.194015e-01 3.788178e-01 2.663642e-01
    ##  [4276] 4.821971e-01 4.123211e-01 3.678829e-01 4.582996e-02 3.694433e-05
    ##  [4281] 1.559899e-01 2.564914e-04 6.548336e-01 3.388099e-04          NaN
    ##  [4286] 1.651745e-02 1.264561e-01 7.592707e-01 4.218335e-01 8.913682e-01
    ##  [4291] 3.477434e-01 1.490685e-35 4.960876e-01 4.995619e-05 6.200956e-03
    ##  [4296] 4.673708e-03 7.844487e-03 1.185114e-01 5.922280e-01 1.309653e-07
    ##  [4301] 4.997940e-04 4.365759e-01 3.204838e-01 1.712510e-29 5.412885e-01
    ##  [4306] 1.797772e-01 4.015599e-02 8.236195e-01 5.775464e-01 8.501932e-02
    ##  [4311] 7.069117e-03 7.088974e-02 1.089630e-01 7.948316e-03 9.518895e-01
    ##  [4316] 9.946202e-01 1.746036e-01 9.619596e-03 4.595427e-01 4.451783e-05
    ##  [4321]          NaN 1.139582e-06 8.041399e-02 4.807165e-01 4.183498e-05
    ##  [4326] 2.504209e-01 1.485129e-03 1.176611e-01 8.312174e-01 1.326923e-18
    ##  [4331] 5.801062e-01 6.810281e-08 1.949205e-08 1.645110e-03 5.774435e-04
    ##  [4336] 1.641124e-06 7.804267e-06 5.028894e-03 7.804112e-01 1.325160e-01
    ##  [4341] 8.709719e-02 2.783381e-03 8.903296e-01 3.075146e-04 1.118368e-02
    ##  [4346] 4.158525e-03 4.828236e-04          NaN 5.990082e-01 1.303548e-53
    ##  [4351] 3.888421e-01 4.877463e-07 1.045731e-01 6.778584e-02 4.324326e-03
    ##  [4356] 5.213319e-01 1.864047e-01 7.102441e-12 4.608310e-01 1.915236e-02
    ##  [4361]          NaN 2.356407e-06 2.352505e-01 3.815672e-01 4.697005e-01
    ##  [4366] 1.370777e-02 3.152087e-01 8.291449e-01 4.046013e-02 6.100481e-01
    ##  [4371] 9.243075e-01 1.025541e-01 2.602366e-11 4.599534e-04 3.456813e-03
    ##  [4376] 4.055567e-02 3.605721e-12 7.375269e-01 9.561134e-04 2.502575e-01
    ##  [4381] 2.337582e-13 9.043079e-01 1.462184e-03 5.624530e-03 6.588574e-01
    ##  [4386] 4.367971e-04 7.794022e-01 8.583265e-01 8.274582e-02 7.009981e-01
    ##  [4391] 5.433657e-01 3.131203e-08 7.261841e-01 3.130611e-04 5.407356e-03
    ##  [4396] 1.072366e-01 7.189313e-01 4.463544e-01 2.229913e-01 2.975221e-01
    ##  [4401] 3.277128e-03 2.287675e-02 5.012593e-03          NaN 2.398216e-03
    ##  [4406] 7.764118e-16 9.494618e-01 1.172366e-01          NaN 2.438217e-01
    ##  [4411] 9.547270e-01 3.329284e-01 2.319317e-03 9.367831e-01 6.882277e-01
    ##  [4416] 7.244930e-01 1.187757e-02 1.800843e-08 2.125703e-01 1.221754e-04
    ##  [4421] 1.061866e-02 3.777900e-07 5.757703e-09 5.080064e-07 4.836933e-06
    ##  [4426] 1.436793e-01 6.615481e-02 6.296329e-03 1.023847e-04 8.143433e-01
    ##  [4431] 6.896741e-01 2.413082e-15 6.712757e-01 2.817626e-19 4.054497e-23
    ##  [4436] 5.642175e-04          NaN 2.545778e-01 8.233845e-04 9.756290e-03
    ##  [4441] 4.256759e-01 8.863413e-02 4.719198e-01 8.806784e-03 1.621591e-12
    ##  [4446] 2.579514e-02 5.449735e-01 7.386689e-03 2.958673e-03 1.442013e-01
    ##  [4451] 5.926599e-01 7.851204e-06 1.279144e-01 7.048474e-01 3.509646e-22
    ##  [4456] 7.159774e-01 4.200704e-01 9.423899e-03 3.535551e-02 1.353249e-01
    ##  [4461] 3.832410e-05 5.194734e-02 2.030610e-02 5.993301e-01 8.939518e-01
    ##  [4466] 8.819053e-01 9.504652e-05 1.660003e-02 7.527769e-02 8.594018e-17
    ##  [4471] 6.669500e-04 5.830622e-01 4.059356e-01 1.675240e-03          NaN
    ##  [4476] 1.498407e-02 9.269511e-01 4.056302e-01 8.250234e-01 5.471239e-01
    ##  [4481] 4.300869e-03 4.212085e-02 4.949888e-03 1.104758e-01 9.321762e-01
    ##  [4486] 3.559441e-01 1.038256e-01 9.311898e-09 3.656933e-02 5.877975e-01
    ##  [4491] 5.087894e-01 1.532483e-01 6.760039e-02 1.788595e-04 1.431608e-02
    ##  [4496] 8.533642e-05 3.917078e-01 1.313379e-02 3.660624e-01 3.059728e-18
    ##  [4501]          NaN 2.576915e-09 3.270059e-05          NaN 9.720607e-01
    ##  [4506]          NaN 5.466712e-01 2.730664e-01 3.395952e-01 7.157912e-02
    ##  [4511] 4.306189e-01 3.181003e-01 7.384536e-06 5.231682e-01 9.378469e-12
    ##  [4516] 2.104494e-05 9.724673e-01          NaN 5.450554e-01 9.982762e-02
    ##  [4521] 1.609137e-02 2.492097e-01 1.303425e-01 5.673124e-04 1.519940e-01
    ##  [4526] 4.478110e-01 9.101595e-01          NaN 1.230662e-01 1.508185e-22
    ##  [4531] 2.909449e-07 2.008987e-03 6.503818e-01 1.100192e-03 2.845965e-09
    ##  [4536] 3.542491e-02 9.481566e-03 8.781836e-17 6.489093e-03 6.582500e-06
    ##  [4541] 2.714137e-01 4.993475e-03 3.204838e-01 9.053954e-03 4.005091e-01
    ##  [4546] 9.959746e-01 8.472440e-01 8.927474e-09 6.074190e-10 3.405749e-02
    ##  [4551] 9.712724e-05 8.193296e-06 3.789973e-05 1.793633e-01 1.069622e-12
    ##  [4556] 2.684673e-01 8.849610e-01 2.012173e-01 3.450538e-05 9.014328e-01
    ##  [4561] 7.198627e-01 7.765862e-06 5.002671e-09 1.655578e-13 3.411715e-03
    ##  [4566] 1.169878e-04 8.137229e-04 7.142151e-03 8.675827e-01 5.373506e-01
    ##  [4571] 2.015207e-01 3.041753e-03 2.832505e-23 6.549995e-11 1.455127e-03
    ##  [4576] 9.497387e-01 1.584691e-01 6.025677e-01 2.011166e-01 4.093755e-02
    ##  [4581] 1.049711e-01 7.889123e-07 6.087072e-01 4.770858e-09 2.658481e-06
    ##  [4586] 5.200636e-23 6.075333e-03 4.016783e-08 8.589333e-02 1.321174e-05
    ##  [4591] 2.679931e-03 8.200175e-01 1.807089e-12 1.263139e-02 7.049659e-01
    ##  [4596] 9.822401e-01 1.257455e-01 1.385762e-04 3.150140e-02 2.499723e-02
    ##  [4601] 2.572180e-01 4.104870e-01 8.248754e-13 1.249112e-01 1.350550e-01
    ##  [4606] 2.177670e-02 8.239161e-01 2.989035e-10          NaN 3.578526e-01
    ##  [4611] 3.265842e-02 2.338153e-01 8.063768e-19 2.453559e-01 7.193786e-13
    ##  [4616] 9.698840e-01 9.704322e-01 2.228992e-02 1.996463e-21 1.005774e-01
    ##  [4621] 3.855389e-01 5.095347e-01 7.899940e-02 3.242945e-03          NaN
    ##  [4626] 5.639221e-04 1.070071e-05 4.382074e-02 6.899459e-01 3.487130e-01
    ##  [4631] 1.206118e-01 4.547451e-01 7.878320e-01 8.931075e-01 6.819604e-58
    ##  [4636] 3.179087e-01 1.946342e-07 1.712882e-01 1.299153e-01 3.804003e-01
    ##  [4641] 5.741304e-08 4.066109e-02 1.124754e-02 4.908217e-01 2.117168e-04
    ##  [4646] 3.389217e-01 2.219724e-02 4.917247e-01          NaN 2.083816e-06
    ##  [4651] 1.076076e-02 2.329726e-01          NaN 2.747527e-01 2.647219e-01
    ##  [4656] 4.241210e-03 3.516895e-04 5.329085e-01 9.895907e-01 1.747675e-30
    ##  [4661] 5.783283e-02 1.536955e-06 9.803490e-01 1.335905e-03 1.817663e-03
    ##  [4666] 5.307820e-01 5.336615e-01          NaN 2.195117e-14 7.482332e-03
    ##  [4671] 1.094550e-01 2.298067e-10 1.444239e-03 9.696466e-01          NaN
    ##  [4676] 5.337783e-01 2.339338e-03 7.502436e-01 4.546576e-06 2.719514e-01
    ##  [4681] 3.840248e-08 9.745388e-01 2.229793e-01 2.046917e-03 8.360374e-01
    ##  [4686] 3.263991e-01          NaN 2.580721e-01 1.524781e-01 7.966420e-05
    ##  [4691] 4.732292e-04 2.002942e-05 3.064200e-01 3.628679e-02 9.633203e-26
    ##  [4696] 6.289987e-15 1.732836e-02 5.792119e-01 2.707511e-18 5.237536e-01
    ##  [4701] 2.173747e-10 2.758682e-01 3.835740e-12 6.562570e-01          NaN
    ##  [4706] 1.172302e-01 7.347497e-03 1.585282e-02 2.027873e-42 2.370992e-01
    ##  [4711] 6.208883e-02 4.113987e-01 5.710998e-02 2.994997e-03 1.036734e-01
    ##  [4716] 8.417261e-02 3.457764e-03 3.161454e-01 1.322026e-09 3.098596e-06
    ##  [4721]          NaN 3.404076e-02 3.987605e-05 5.611912e-01 1.093020e-01
    ##  [4726] 4.005269e-01 5.376520e-02 1.011449e-06 3.614638e-03          NaN
    ##  [4731]          NaN 9.767266e-02 1.107528e-04 6.329343e-01 1.583472e-01
    ##  [4736] 7.363107e-01 6.822824e-01 4.246598e-01 9.806684e-01 5.872431e-01
    ##  [4741] 2.637175e-07 1.785536e-47 8.626654e-26 8.479910e-12 1.468633e-01
    ##  [4746] 5.800848e-01 3.348598e-07 3.379259e-01 8.890833e-01 9.979929e-05
    ##  [4751] 5.421421e-01 4.837037e-01 2.274067e-01 5.796623e-03 3.778250e-01
    ##  [4756] 3.445323e-01 6.113167e-01 9.440831e-16 7.789690e-04 1.207718e-02
    ##  [4761] 3.930647e-01 7.292149e-01 1.025260e-04 3.893943e-03 4.045863e-01
    ##  [4766] 3.990585e-10 2.570237e-01 5.816960e-01 2.949910e-10 8.183670e-02
    ##  [4771] 8.539103e-05 1.786963e-01 9.370116e-16 7.454116e-01 3.491086e-05
    ##  [4776] 4.636609e-04 2.882615e-03 3.711610e-01 1.424827e-02 3.064550e-03
    ##  [4781] 2.196237e-03 7.787990e-03 7.548439e-01 9.506964e-01 1.393240e-05
    ##  [4786] 3.044463e-01 3.204838e-01 8.821920e-01 5.700337e-01 1.383672e-01
    ##  [4791] 9.087948e-04 3.775388e-04 4.123929e-02 8.236489e-02 7.805697e-01
    ##  [4796] 1.477951e-03 9.771981e-08          NaN 1.738663e-01 9.251003e-01
    ##  [4801] 6.619841e-02 5.528264e-02 5.709560e-03 1.637416e-02 9.932738e-01
    ##  [4806] 1.110603e-04 1.619529e-19 2.308626e-03 1.889863e-01 1.559114e-01
    ##  [4811] 5.776169e-01 2.660092e-07 4.206960e-02 5.619607e-01 7.655504e-02
    ##  [4816] 3.260573e-01 2.528775e-06 2.637057e-01 8.827299e-01 3.697381e-01
    ##  [4821] 1.268107e-08 5.882069e-02 3.862141e-01 3.572379e-14 1.220968e-05
    ##  [4826] 3.764278e-01 6.282062e-03 1.661083e-12 3.096669e-23 5.183080e-01
    ##  [4831] 1.947817e-05 4.700222e-19 3.913092e-01 1.435822e-01 9.405552e-05
    ##  [4836] 1.121055e-05          NaN 5.872523e-02 2.328221e-01 6.503029e-01
    ##  [4841] 1.842214e-02 4.452258e-02 2.373651e-01 9.943699e-02 8.363530e-01
    ##  [4846] 2.295125e-01 8.135774e-10 6.743637e-01 6.350073e-01 3.623684e-01
    ##  [4851] 9.421796e-03 8.715442e-01 6.479698e-03 9.902429e-05 6.092965e-01
    ##  [4856] 1.390108e-03 8.741253e-01 5.287084e-07 8.483672e-24 4.738615e-08
    ##  [4861] 3.510003e-02 9.880494e-02 9.579057e-07 6.436661e-01 2.588350e-09
    ##  [4866] 6.544912e-04 3.426785e-03          NaN 5.849509e-01 7.875659e-01
    ##  [4871] 1.413952e-02 5.387507e-11 9.170103e-03 9.260485e-16 1.080262e-01
    ##  [4876] 6.718540e-22 3.307565e-04 4.939125e-02 6.732082e-01 6.534091e-01
    ##  [4881] 3.643252e-01 4.457358e-01 8.376656e-02 3.474264e-03 3.400517e-01
    ##  [4886] 3.834636e-02 1.282726e-01 8.133751e-01 8.456908e-01 3.311783e-02
    ##  [4891] 1.597698e-01 2.587384e-01 1.423352e-01 2.268539e-01 5.048793e-02
    ##  [4896] 2.447287e-01 2.256541e-01 7.915592e-12 2.334418e-06 2.431898e-01
    ##  [4901] 4.472025e-02 2.308709e-02 5.061895e-03 6.650203e-01 8.650394e-02
    ##  [4906] 1.219974e-01 7.316837e-01 6.420765e-11 2.574944e-01 1.505779e-12
    ##  [4911] 1.309643e-04 4.783550e-04 8.434814e-01 8.662824e-02 2.501473e-01
    ##  [4916] 4.551138e-07 8.022912e-08 6.321883e-02 2.868724e-03 1.261606e-06
    ##  [4921] 2.718593e-02 8.138164e-01 8.286254e-02          NaN 5.532583e-01
    ##  [4926] 2.327655e-05          NaN 2.849800e-01          NaN 1.583393e-07
    ##  [4931] 2.462160e-03 4.782131e-01 8.124246e-01 1.520387e-01 3.121182e-01
    ##  [4936] 2.521147e-01 1.172832e-02 1.892724e-17 1.159934e-01 1.263013e-02
    ##  [4941] 3.152513e-04 2.993955e-16 2.219846e-04 2.071261e-11 8.723729e-10
    ##  [4946] 3.078253e-01 2.481964e-02 5.276134e-02          NaN 5.811813e-01
    ##  [4951] 9.997863e-17 1.036869e-09 1.262204e-01 2.812701e-01 8.626513e-01
    ##  [4956] 3.326002e-15 6.192292e-01 3.960033e-14 3.913894e-02 1.190568e-02
    ##  [4961] 2.386170e-01 5.290392e-02 2.558675e-01 4.059018e-02 1.386265e-05
    ##  [4966] 3.787798e-01 9.464132e-01 3.559364e-02 8.844130e-01 4.791073e-01
    ##  [4971] 3.201223e-09 3.132883e-01 1.973061e-05 2.075244e-03 3.679407e-01
    ##  [4976] 8.950049e-02 6.146495e-01 5.668250e-10 2.552457e-02 1.269053e-02
    ##  [4981] 5.750694e-02 5.006164e-02 9.124149e-01 3.254102e-01 2.224065e-01
    ##  [4986] 1.844224e-03 9.880117e-03 5.676597e-02 1.529613e-02          NaN
    ##  [4991] 2.974682e-01 7.051439e-01 8.249013e-01 5.888240e-13 1.826462e-01
    ##  [4996] 1.458045e-03 3.935041e-03 1.011502e-01 2.497794e-05 4.501294e-01
    ##  [5001] 1.810184e-20 8.320999e-02 2.521319e-02 7.937514e-01 2.286439e-05
    ##  [5006] 2.079269e-06 3.204838e-01 2.622577e-01 7.610174e-01 3.587703e-02
    ##  [5011] 2.330765e-12          NaN 5.505435e-01 5.079028e-49 3.919730e-11
    ##  [5016] 5.572782e-22 3.100903e-02 2.636855e-03 6.701036e-06 2.324718e-01
    ##  [5021] 2.146678e-24 9.524875e-03 4.678965e-13 5.531588e-04 7.914413e-02
    ##  [5026] 3.204838e-01 8.401408e-01 5.144738e-01 1.102742e-01 3.853001e-01
    ##  [5031] 6.415240e-01 1.069207e-01 3.352294e-02 4.905838e-02 1.672866e-28
    ##  [5036] 8.027923e-03 4.100320e-01 5.105860e-03 8.107392e-02 2.208297e-03
    ##  [5041] 4.635107e-04 6.592807e-02 2.390386e-02 1.731119e-04 1.169051e-02
    ##  [5046] 6.504571e-02 1.505989e-11 6.618206e-01 9.732691e-01 8.073811e-01
    ##  [5051] 8.155458e-01 1.421619e-02 6.299248e-01 1.083311e-29 2.779056e-16
    ##  [5056] 2.362001e-43 2.042830e-02 2.859536e-02 4.060151e-03 9.215121e-01
    ##  [5061] 2.455766e-04 6.473925e-03 8.096753e-17 2.968171e-01 8.133155e-05
    ##  [5066] 3.958923e-01 2.917311e-55 6.160937e-04 2.280089e-19 3.068821e-02
    ##  [5071] 3.431421e-08 1.931437e-04 4.691603e-03 2.380644e-06 6.832292e-01
    ##  [5076] 1.069721e-04 4.784189e-01 2.301601e-01 5.968144e-01 4.248862e-16
    ##  [5081] 3.116767e-01 1.007168e-02 1.706687e-02 4.379437e-04 9.764834e-03
    ##  [5086] 3.779536e-01 6.670305e-04 9.396741e-01 3.170440e-02 7.878142e-02
    ##  [5091] 9.040112e-05 4.361829e-01 3.521951e-01 2.199925e-04 4.102109e-06
    ##  [5096] 6.822983e-01 3.920010e-01 1.085546e-01 9.206564e-01 3.908483e-02
    ##  [5101]          NaN 1.996283e-04 3.457875e-02 6.075243e-02 1.102275e-03
    ##  [5106] 5.176544e-03 6.409057e-05 1.546409e-02 3.460269e-01 3.141382e-01
    ##  [5111] 3.204838e-01          NaN 3.827091e-04          NaN 4.705643e-19
    ##  [5116] 2.075655e-01 3.611535e-02 2.815963e-01 9.434063e-01          NaN
    ##  [5121] 2.011779e-05 8.340342e-02 9.409304e-01 1.710971e-08 4.574789e-01
    ##  [5126] 1.479020e-02 6.076065e-01 1.387800e-03 6.112894e-02 5.351616e-01
    ##  [5131] 9.616681e-01 5.888627e-01 9.160716e-05 1.503857e-01 9.753947e-05
    ##  [5136] 4.655690e-08 3.980840e-03 3.497544e-02 8.204790e-01 8.361452e-08
    ##  [5141] 3.024737e-01 8.172131e-01 3.423099e-01 3.204838e-01 9.866962e-01
    ##  [5146] 8.138794e-01 4.779772e-03 3.391673e-01 9.078227e-01 5.966796e-01
    ##  [5151] 9.932769e-01 4.834140e-03 2.679444e-02 1.129312e-01 1.030914e-02
    ##  [5156] 6.600645e-07 3.159261e-08 2.183917e-06 3.112873e-01 9.989510e-01
    ##  [5161] 1.403612e-01 9.248150e-01 4.891957e-06 9.749623e-01 2.281840e-02
    ##  [5166] 8.195699e-01 1.968325e-03 6.963647e-01 5.500661e-01 8.771353e-06
    ##  [5171] 5.323979e-02 8.368053e-01 1.707451e-03 1.603035e-01 1.321265e-15
    ##  [5176] 9.694078e-01 2.706216e-01 4.512553e-06 1.181681e-13 2.416359e-26
    ##  [5181] 4.930210e-01 2.273493e-34 3.204838e-01 9.295343e-20 4.698902e-06
    ##  [5186] 1.365440e-01 2.404971e-06 1.723939e-02 3.607649e-02 5.881240e-01
    ##  [5191] 7.809034e-14 8.108931e-07 6.937759e-07 1.766196e-11 4.153918e-01
    ##  [5196] 7.489131e-04 4.603955e-05 6.297567e-09 4.043225e-03 4.136847e-29
    ##  [5201] 3.632112e-03 6.131087e-01 9.201819e-02 5.162257e-01 3.204838e-01
    ##  [5206] 5.855431e-01 7.572212e-01 2.842309e-02 2.500933e-01          NaN
    ##  [5211] 5.757942e-12 8.165935e-01 2.101924e-02 6.150129e-02 4.240803e-02
    ##  [5216] 6.375109e-01 1.557513e-01 2.650287e-02 9.115447e-04          NaN
    ##  [5221] 1.863734e-10 8.636125e-11 5.479063e-03 2.562176e-02 1.367639e-04
    ##  [5226] 5.634130e-05 4.018266e-01 3.324749e-01 3.097445e-01 3.069887e-01
    ##  [5231] 8.597497e-06 1.327506e-10 8.381671e-01 1.946426e-07 2.647044e-01
    ##  [5236] 9.509673e-01 7.980978e-09 3.204838e-01 2.751135e-05 7.695017e-07
    ##  [5241] 1.861047e-04 1.918358e-02 2.477683e-10 1.109334e-09 9.442945e-11
    ##  [5246] 2.232241e-01 8.203452e-04 1.406934e-02 6.413735e-03 5.935165e-02
    ##  [5251] 8.897485e-01 1.563334e-01 3.066267e-05 7.659156e-02 8.130066e-05
    ##  [5256] 2.396700e-01 2.297010e-01 9.400239e-06 2.737502e-03 2.894622e-01
    ##  [5261] 7.278958e-02 8.873366e-19 3.204838e-01 9.275139e-01 7.318317e-01
    ##  [5266] 2.299451e-06 7.712814e-04 3.602315e-01 1.749198e-16 2.109324e-02
    ##  [5271] 3.697691e-03          NaN 9.407185e-01 7.597640e-01 1.151897e-02
    ##  [5276] 7.140946e-04 1.702736e-13 4.598205e-01 1.999036e-06 4.897219e-01
    ##  [5281] 2.600357e-01 2.876752e-01 6.222385e-07 6.840032e-01 8.543360e-01
    ##  [5286] 5.777867e-11 1.261809e-07 2.828779e-10 1.888975e-09 1.831924e-03
    ##  [5291] 5.910952e-02 2.127661e-03 6.378921e-02 1.850004e-02 1.535648e-01
    ##  [5296] 3.945935e-01 1.272703e-01 4.814438e-02 3.522621e-01 2.753695e-02
    ##  [5301] 6.943693e-09 2.322909e-01 1.803627e-01 2.008040e-03 9.647547e-07
    ##  [5306]          NaN 8.277128e-05 1.636342e-01 8.405559e-03 9.345928e-02
    ##  [5311] 5.348187e-01 2.734676e-01 1.127642e-01 4.472541e-18 9.113568e-01
    ##  [5316] 1.279701e-01 2.705453e-17 1.575541e-01          NaN 1.663703e-01
    ##  [5321] 5.971456e-16 9.152045e-02 1.238696e-02 3.349431e-04 7.431096e-11
    ##  [5326] 7.556651e-01 8.151307e-02 6.513943e-03 5.785953e-01 4.159600e-01
    ##  [5331] 3.902748e-02 3.862374e-01 5.609007e-01 4.474751e-01 1.430140e-04
    ##  [5336] 2.069792e-05 5.146643e-02 2.900788e-02 2.090545e-01 5.997452e-01
    ##  [5341] 9.180511e-03          NaN 5.401979e-06          NaN 2.690501e-10
    ##  [5346] 2.918468e-14 4.621409e-01 3.912034e-02 2.982478e-02 1.317553e-03
    ##  [5351] 3.229245e-03 1.913755e-03 2.600422e-02 2.147508e-35 5.158429e-01
    ##  [5356] 3.244383e-02          NaN 9.433980e-01 7.818353e-01 1.135506e-15
    ##  [5361] 1.458648e-07 4.471580e-04 1.617111e-03 6.431854e-05 4.477995e-01
    ##  [5366] 7.614640e-01 1.666459e-01 2.205636e-01 3.257882e-05 1.877876e-01
    ##  [5371] 5.585889e-10 4.318269e-04 2.733422e-01 9.461661e-01 7.741188e-03
    ##  [5376] 4.046929e-01 6.974413e-02 2.240684e-01 6.008690e-01 6.529481e-02
    ##  [5381] 7.564911e-04 7.726757e-03 1.114059e-01 3.335916e-01 5.935481e-01
    ##  [5386] 7.640198e-01 1.796703e-26 7.228863e-01 7.127609e-01 3.642913e-01
    ##  [5391] 9.923521e-01 5.022977e-01 5.447674e-01 8.527535e-01 8.794484e-06
    ##  [5396] 3.374781e-01 7.459782e-05 6.358916e-02 6.182612e-01 2.247181e-02
    ##  [5401] 3.511973e-01 3.308521e-02 3.768857e-02 1.172912e-03 5.660652e-16
    ##  [5406] 6.172245e-01 1.339009e-10          NaN 3.645768e-01 3.744924e-02
    ##  [5411] 7.589774e-01 1.173234e-03 5.921676e-05 9.614445e-03 7.449555e-05
    ##  [5416] 4.290812e-03 2.089075e-07 3.876033e-01 4.812754e-01 6.457460e-01
    ##  [5421] 6.206823e-02 5.944055e-01 2.222016e-04 1.425855e-01 2.993169e-03
    ##  [5426] 9.268132e-01 1.071856e-01 4.062076e-01 3.001696e-03 6.618020e-01
    ##  [5431] 1.201302e-04 1.176798e-06 2.249478e-02 4.514305e-01 1.977224e-15
    ##  [5436] 5.675840e-01 4.696552e-04 6.613012e-02 1.496283e-07 5.044433e-01
    ##  [5441] 3.204838e-01 9.895943e-06 5.488900e-01 4.558178e-01 8.661726e-01
    ##  [5446] 3.804838e-01 8.716135e-04 2.383923e-03 6.786100e-07 2.521011e-01
    ##  [5451] 6.284223e-04 9.940456e-09 3.500322e-06 3.390190e-02 3.149658e-08
    ##  [5456] 3.666686e-01 9.255524e-01 1.489495e-10 5.151746e-01 8.404857e-01
    ##  [5461] 5.457649e-03 1.999160e-07 1.909434e-02 2.485927e-01 6.394591e-01
    ##  [5466] 8.505792e-04 5.232417e-03 1.774048e-01 7.282406e-01 4.178421e-01
    ##  [5471] 4.843136e-01 2.542273e-01 1.265896e-15 1.575370e-22 2.300068e-01
    ##  [5476] 6.521611e-01 4.687691e-08 7.322857e-02 2.115967e-02 9.863499e-07
    ##  [5481] 4.486688e-01 4.062738e-01 7.969479e-01 4.354424e-02 8.058022e-02
    ##  [5486] 2.105402e-02 6.041126e-01 7.924533e-01 1.145370e-01 2.199808e-01
    ##  [5491] 5.061180e-25 2.111910e-03 1.274072e-01 6.789041e-01 1.186287e-02
    ##  [5496] 4.017169e-03 1.370374e-02 2.574302e-05 1.108370e-01 3.035524e-01
    ##  [5501] 8.189315e-01 1.190793e-01 1.398913e-01 1.402908e-04 2.239321e-02
    ##  [5506] 2.732547e-07 4.468582e-01 1.472294e-06 4.461138e-01 1.667081e-04
    ##  [5511] 9.028564e-02          NaN 8.985567e-02 9.920658e-03 4.902678e-02
    ##  [5516] 4.127375e-01 6.679160e-31 3.204838e-01 1.368738e-04 8.390482e-01
    ##  [5521] 5.954503e-02 6.815765e-01 4.466477e-01 6.718348e-01 3.737871e-01
    ##  [5526] 8.825512e-01 8.253516e-01 1.653943e-01 2.119023e-02 1.773070e-02
    ##  [5531] 1.455730e-03 9.289707e-03 1.167349e-05 3.374446e-03 5.777427e-05
    ##  [5536] 1.849828e-01 5.186727e-01 4.855553e-12 2.384253e-06 6.948863e-01
    ##  [5541] 3.507167e-15 1.336052e-01 2.948146e-01 1.898527e-01 1.035479e-02
    ##  [5546] 9.984591e-03 6.845315e-05 8.160002e-01 6.334724e-13 7.693939e-04
    ##  [5551] 3.461770e-01 5.797669e-01          NaN 5.786212e-02 2.513186e-02
    ##  [5556] 7.906161e-01 9.273460e-02 6.328907e-02 6.962925e-05 1.968776e-04
    ##  [5561] 3.633210e-04 7.375988e-02 2.926003e-02 3.204838e-01 3.603971e-02
    ##  [5566] 1.213838e-03 4.809144e-01 2.720020e-01 3.284310e-03 5.746093e-03
    ##  [5571] 5.610111e-02 9.586752e-01 1.106027e-02 1.094705e-06 3.752961e-02
    ##  [5576] 9.918263e-01 1.077459e-02 2.392865e-01 3.021061e-02 1.043945e-01
    ##  [5581]          NaN 1.655549e-38 6.434763e-03 1.687628e-01 6.173596e-01
    ##  [5586] 8.877673e-02 1.548716e-01 4.078602e-04 2.793755e-02 4.992483e-01
    ##  [5591] 1.774324e-08 5.121619e-01 3.690779e-01 7.828441e-02 2.924164e-01
    ##  [5596] 7.735359e-01 2.218667e-15 9.786975e-01 8.843207e-05 6.527726e-01
    ##  [5601] 4.094564e-02 5.581713e-06 3.779871e-04 3.264956e-01 6.859793e-05
    ##  [5606] 6.362899e-01 3.667503e-02 4.474150e-17 1.619967e-01 2.200986e-03
    ##  [5611] 8.842031e-01 7.191470e-01 1.979787e-02 2.241009e-06 3.848468e-02
    ##  [5616] 2.473334e-04 8.032949e-01 4.687141e-01 1.021454e-01 2.286662e-17
    ##  [5621] 1.373895e-07 7.270462e-02 3.247120e-01 1.461767e-03 3.989702e-05
    ##  [5626] 1.115353e-09 1.125747e-08 5.483716e-01 1.268742e-01 9.823179e-05
    ##  [5631] 2.117803e-07 6.224140e-02 1.656062e-23 2.484526e-02 8.424455e-01
    ##  [5636] 1.126231e-01 3.458073e-02 9.492676e-01 3.532422e-01 9.853766e-13
    ##  [5641] 2.813278e-01 8.866003e-01 8.452659e-01 1.038577e-01 6.345345e-07
    ##  [5646] 6.558838e-02 8.695113e-02 2.614312e-04 2.657396e-01 2.917393e-01
    ##  [5651] 2.153825e-05 4.785732e-02 1.749357e-11 4.857086e-07 8.961340e-01
    ##  [5656] 1.813896e-01 1.473017e-12 2.458267e-06 5.160299e-01 4.206889e-06
    ##  [5661] 9.726166e-01 3.289236e-04 5.701771e-04 3.163664e-08 3.084630e-01
    ##  [5666] 8.637850e-01 2.323726e-02 4.426917e-03 3.130400e-09 3.954455e-02
    ##  [5671] 1.113596e-01 1.700913e-06 6.110284e-01 1.467617e-19 5.190432e-01
    ##  [5676] 3.518600e-01 4.649873e-04 5.343744e-01 1.544030e-21 9.992556e-01
    ##  [5681] 1.138994e-03 7.238766e-01 7.539601e-02 1.145378e-08 7.182950e-13
    ##  [5686] 4.061887e-01 6.852297e-08 5.942377e-01 1.998410e-03 4.242002e-02
    ##  [5691] 1.954511e-05 6.978726e-03 4.237406e-09 2.579119e-03 6.061843e-06
    ##  [5696] 2.805383e-01 2.692865e-03 1.260356e-03 2.390475e-03 2.106878e-08
    ##  [5701] 2.450008e-01 4.010107e-02 6.529697e-01 5.348795e-06 1.800685e-06
    ##  [5706] 4.251222e-01 9.807631e-03 2.149140e-01          NaN 4.169294e-01
    ##  [5711] 1.682571e-01 1.042236e-05 5.306216e-04 9.331310e-02 2.451784e-18
    ##  [5716] 7.062851e-03 6.695554e-01 5.118308e-03 7.828099e-05 3.204838e-01
    ##  [5721] 8.358632e-02 2.027742e-33 7.632030e-01 4.625382e-02 4.175816e-01
    ##  [5726] 8.922580e-01 9.944995e-04 7.617804e-01 1.141949e-13 6.473628e-01
    ##  [5731] 1.405165e-05 1.158471e-01 4.184504e-05 3.740099e-04 7.760115e-01
    ##  [5736] 6.167983e-05 7.832837e-01 1.808196e-01 4.718917e-09 7.914120e-06
    ##  [5741] 1.315273e-01 6.464903e-09 2.138201e-01 9.809249e-01 6.785089e-01
    ##  [5746] 9.163526e-01 7.110695e-01 3.552050e-22 1.546685e-10 1.069212e-02
    ##  [5751] 5.615600e-01 3.204838e-01 1.341456e-10 2.118048e-02 4.039154e-01
    ##  [5756] 6.657793e-01 1.742742e-01 3.448393e-04 5.989995e-09 4.170204e-01
    ##  [5761] 3.532862e-01 5.957911e-01 8.271142e-03 1.310966e-22 8.241321e-01
    ##  [5766] 6.380655e-03 4.448944e-03 1.672952e-05 4.480039e-01 9.479408e-07
    ##  [5771] 1.969280e-04 5.758259e-07 1.322404e-10 1.412123e-04 3.203360e-04
    ##  [5776] 8.482649e-01 3.181317e-03 5.378257e-03 2.113181e-02 8.073650e-05
    ##  [5781] 2.850066e-01 4.744262e-01 3.176267e-01 3.870667e-03          NaN
    ##  [5786] 5.869695e-01 4.943717e-03 7.389859e-01 3.225394e-19 6.736648e-02
    ##  [5791]          NaN 1.392737e-04 6.083912e-01 4.270157e-02 1.167967e-05
    ##  [5796] 9.822390e-01 3.233528e-03 6.120207e-01 9.102010e-02 1.624080e-02
    ##  [5801] 4.605931e-01 1.749126e-01 9.805957e-02 1.354017e-04 4.994234e-02
    ##  [5806] 4.089903e-02 2.189819e-01 5.933309e-01 7.076186e-01 1.579095e-01
    ##  [5811] 1.014156e-01 1.347981e-01 9.741079e-01 9.518622e-01 2.091604e-02
    ##  [5816] 3.159299e-01 1.725943e-02 1.160599e-01 3.866065e-02 9.963880e-04
    ##  [5821] 7.107334e-07 3.664981e-01 9.078832e-19 1.079447e-03 6.813083e-01
    ##  [5826]          NaN 5.383500e-01          NaN 1.283334e-01 6.090322e-01
    ##  [5831] 5.174376e-07 3.626289e-01 9.037045e-07 3.701616e-02 9.660926e-02
    ##  [5836] 1.511016e-04 2.375446e-02 9.010389e-01 8.735721e-12 7.936357e-02
    ##  [5841] 8.335615e-01 3.203654e-01          NaN 1.080987e-01 1.455131e-08
    ##  [5846] 4.252884e-03 2.199474e-01 5.690665e-01 7.027260e-02 8.621621e-02
    ##  [5851] 1.455209e-02 8.353771e-01 3.999036e-07 6.877742e-05 1.734314e-01
    ##  [5856] 9.560075e-01 2.336629e-13 4.806127e-05 5.973135e-01 1.721199e-01
    ##  [5861]          NaN 1.064246e-02 6.645761e-11 4.539208e-01 3.335569e-01
    ##  [5866] 1.690275e-01 3.644741e-06 8.137065e-01 5.283534e-02 3.395021e-03
    ##  [5871] 5.035884e-01          NaN 1.693511e-03 1.986995e-01 9.060736e-02
    ##  [5876] 5.557110e-01 7.629998e-04 1.022258e-02 2.301075e-06 1.658711e-02
    ##  [5881] 1.043167e-20 8.364117e-01 4.001719e-01 2.029598e-04 5.521852e-01
    ##  [5886] 1.780005e-01 1.764848e-01 2.254656e-03 2.002916e-04 1.965355e-01
    ##  [5891] 5.960725e-02 1.887910e-04 3.056213e-04 8.656512e-01 1.794215e-01
    ##  [5896] 3.588048e-01 1.580269e-02 6.826176e-03 1.277467e-05 5.607122e-08
    ##  [5901] 1.364839e-01 9.081807e-01 2.671944e-01 3.204838e-01 9.347202e-01
    ##  [5906] 5.405877e-10 5.493545e-05 5.718749e-01 1.911899e-02 3.094124e-01
    ##  [5911] 9.406183e-04 1.707276e-03 3.509880e-02 3.750915e-01 6.890528e-01
    ##  [5916] 8.103608e-03 5.909553e-04 5.365711e-01 1.629021e-01 3.925126e-02
    ##  [5921] 3.204838e-01 1.406252e-03 3.897991e-08 5.538366e-02 2.927630e-01
    ##  [5926] 2.038756e-07 1.234608e-04 1.932364e-05 1.686459e-02 1.730515e-06
    ##  [5931] 3.298005e-01 1.997930e-06 7.572685e-12 3.900946e-20 2.060232e-01
    ##  [5936] 4.778784e-09 6.047729e-01 8.070922e-01 1.039808e-01 4.976039e-06
    ##  [5941] 6.173535e-01 5.890463e-02 6.377886e-01 5.676976e-01 5.143292e-29
    ##  [5946] 1.029823e-03 4.142392e-01 4.558477e-02 7.582147e-01 1.364178e-02
    ##  [5951] 1.833454e-01 8.290476e-02 6.891079e-01 1.586499e-04 1.052307e-01
    ##  [5956] 9.526185e-01 4.372182e-01 1.620014e-01 7.729878e-01 4.923249e-01
    ##  [5961] 6.172660e-04 2.719840e-05          NaN 1.863154e-02 1.704253e-01
    ##  [5966] 9.745318e-01          NaN 1.361275e-03          NaN 1.777447e-01
    ##  [5971] 8.872031e-03 4.108733e-01 2.849556e-02          NaN 6.916060e-10
    ##  [5976] 6.953013e-01          NaN 1.747732e-01 3.070858e-03 7.653345e-02
    ##  [5981] 1.156077e-01 2.906967e-01 2.165253e-06 1.098455e-05          NaN
    ##  [5986] 8.169115e-01 3.258694e-02 6.294270e-02 8.743617e-01 2.995323e-03
    ##  [5991] 3.362575e-01 9.019571e-01 3.949859e-07 6.046085e-01 1.360586e-03
    ##  [5996] 2.445657e-06 1.044864e-02          NaN 1.706403e-05 6.385900e-01
    ##  [6001] 1.148780e-01 8.096111e-01 1.374852e-15 2.996827e-05 2.621580e-02
    ##  [6006] 1.952975e-08 1.215174e-02 1.759418e-04 7.779779e-01 5.688452e-04
    ##  [6011] 2.627766e-01 9.874116e-01 6.254054e-01 4.437835e-01 9.015464e-04
    ##  [6016] 1.569003e-04 3.489888e-08 3.575524e-17 3.268411e-01 4.443505e-02
    ##  [6021] 5.862856e-10 7.523157e-01 1.965151e-01 2.185339e-02 4.571383e-01
    ##  [6026] 8.031326e-01 2.009474e-02 8.066346e-01 2.966456e-01 2.254971e-03
    ##  [6031] 1.386957e-01 6.535719e-02 1.799810e-01 6.564784e-01 3.489576e-03
    ##  [6036] 1.046785e-04 1.145434e-11 3.204838e-01 7.156957e-01 1.652779e-01
    ##  [6041] 4.649615e-01 7.280646e-01 7.493922e-05 7.425506e-02 3.551243e-01
    ##  [6046] 1.970413e-01          NaN 7.752290e-01 3.464168e-01 1.687511e-04
    ##  [6051] 2.776595e-01 3.543523e-01 1.677221e-23 9.458112e-01 1.654433e-01
    ##  [6056] 3.634102e-04 4.932142e-06 7.578562e-01 9.460401e-01 3.090599e-03
    ##  [6061] 4.592382e-03 2.319014e-01 4.440189e-16 1.749507e-05 3.657521e-04
    ##  [6066] 2.659478e-02 4.148315e-01 8.054038e-27 8.082645e-02 1.234135e-01
    ##  [6071] 4.601246e-03 4.063133e-01 7.067523e-02 1.500314e-08 8.227522e-01
    ##  [6076] 9.912261e-01 2.236888e-01 2.821788e-04 3.258965e-04 9.433417e-01
    ##  [6081] 6.184870e-04 4.832958e-01 1.125895e-02 2.748492e-03 6.656857e-01
    ##  [6086] 1.718997e-09 2.479816e-04 5.917576e-01 4.004088e-16 2.257437e-01
    ##  [6091] 2.961469e-02 2.587458e-01 3.963582e-02 1.093052e-16 7.709576e-18
    ##  [6096] 8.058294e-16 3.066941e-02 1.953190e-08 9.158664e-01 1.116514e-03
    ##  [6101] 7.651251e-03 2.880124e-04 9.660469e-01 4.812209e-12 2.203234e-04
    ##  [6106] 3.269329e-05 1.270592e-01 1.015061e-04 3.210891e-04 2.390585e-01
    ##  [6111] 8.151316e-01 2.591118e-02 2.630918e-01 2.725532e-01 3.607594e-02
    ##  [6116] 5.947817e-01 1.028704e-09 5.088930e-05 1.391566e-04 6.445567e-01
    ##  [6121] 2.756861e-09 7.741641e-01 9.880448e-02 4.133043e-01          NaN
    ##  [6126] 1.110717e-01 1.014355e-02 1.715901e-09 2.128246e-02 3.576266e-02
    ##  [6131] 2.849269e-01 8.236675e-03 3.005052e-02 3.312112e-01 3.892196e-13
    ##  [6136] 5.481237e-21 9.886926e-03 2.502357e-06 2.284148e-01 1.918008e-01
    ##  [6141] 6.568942e-02 2.425163e-03 2.079239e-09 2.594568e-01 7.058803e-01
    ##  [6146] 2.506036e-04 8.151086e-01 9.526575e-04 1.122189e-01 5.314471e-06
    ##  [6151] 1.108839e-03 1.592162e-01 4.699109e-02 9.339082e-02 1.515842e-01
    ##  [6156] 2.635233e-01 1.544971e-03 5.513447e-07 7.376147e-06 4.534776e-07
    ##  [6161] 1.569804e-04 1.108846e-07 4.970476e-01 4.480651e-05 2.542188e-03
    ##  [6166] 2.546632e-07 1.013268e-04 1.996706e-01 4.709457e-05 1.200424e-01
    ##  [6171] 2.460191e-01 1.647133e-01 1.469962e-02 1.296108e-01 9.645625e-02
    ##  [6176] 2.804125e-10 3.103627e-06 4.650144e-27 3.414542e-05 1.020962e-01
    ##  [6181] 4.094434e-01 5.766716e-01 3.427926e-01          NaN 8.005782e-01
    ##  [6186] 1.411540e-02 7.909471e-01 6.191265e-02 2.540884e-05 5.826269e-01
    ##  [6191] 2.369505e-01          NaN 4.331644e-11 2.761635e-01 4.329286e-23
    ##  [6196] 9.959357e-03 9.826165e-01 3.864586e-03 7.857576e-28 7.094784e-01
    ##  [6201]          NaN 4.076917e-03 3.012640e-01 2.893721e-22 3.312052e-02
    ##  [6206] 3.019081e-11 1.913869e-14 7.067247e-13 9.668268e-01 3.585692e-03
    ##  [6211] 4.319089e-01 7.779618e-03 2.609222e-02 1.868906e-07 1.202927e-19
    ##  [6216] 2.586971e-02 3.830628e-01 8.295895e-01 9.862211e-02 9.703928e-09
    ##  [6221] 2.058089e-01 9.787911e-02 4.265735e-01 4.824798e-05 2.673570e-04
    ##  [6226] 3.204838e-01 7.566007e-01 3.583667e-01 5.780514e-10          NaN
    ##  [6231] 7.552367e-02 9.018496e-03          NaN 1.542078e-04          NaN
    ##  [6236] 6.282210e-04 8.590311e-02 1.583786e-03 3.719940e-01 8.100144e-05
    ##  [6241] 1.515883e-01 5.971625e-02 3.509468e-02 7.117838e-01 8.310368e-26
    ##  [6246] 5.937616e-01 1.423157e-01 2.923664e-07          NaN 1.965791e-01
    ##  [6251] 1.597762e-01 9.119519e-01 6.735229e-07 4.989150e-01 8.523594e-01
    ##  [6256] 5.439237e-03 9.843764e-03 7.036785e-01 4.843034e-01 3.969165e-01
    ##  [6261] 6.852111e-01          NaN 2.012636e-01 6.780765e-04 3.163087e-01
    ##  [6266] 5.686139e-01 4.159770e-25 2.443487e-02 1.184904e-01 3.876951e-02
    ##  [6271] 3.200534e-01 3.196244e-03 2.659517e-04 2.079313e-02 1.383225e-05
    ##  [6276] 4.428587e-04 1.171511e-03 8.670226e-02          NaN 5.931606e-01
    ##  [6281] 2.086235e-13 9.576103e-04 6.792862e-01          NaN 1.743621e-09
    ##  [6286] 4.027590e-08 9.043991e-01 5.038296e-02 2.935009e-06 1.860995e-02
    ##  [6291] 2.420213e-01 7.293103e-02          NaN 8.978834e-06 9.835727e-04
    ##  [6296] 8.286856e-02 3.514740e-01 6.758442e-01 9.454056e-02 6.960006e-02
    ##  [6301] 8.551695e-01 3.569232e-27          NaN 1.315891e-10 2.164364e-01
    ##  [6306] 6.348671e-44 1.084281e-02 7.418184e-05 1.789321e-01 4.267640e-10
    ##  [6311] 1.953592e-07 9.367966e-01 9.728884e-04 1.076921e-05 2.464401e-04
    ##  [6316] 7.350743e-02 8.574433e-02 9.785109e-02 6.904734e-03 1.181902e-01
    ##  [6321] 5.353926e-01 8.243642e-05 1.950994e-01 2.248216e-10 2.621571e-01
    ##  [6326] 7.152438e-01 6.467614e-01 2.844353e-16 5.651373e-01 1.719916e-02
    ##  [6331] 2.908532e-01 1.441264e-01 4.505489e-03 6.094229e-04 1.230701e-08
    ##  [6336] 4.559780e-01 1.594321e-07          NaN 3.006536e-01 3.463812e-01
    ##  [6341] 2.220093e-01 1.387834e-02 3.340489e-01 2.106875e-01 7.766620e-03
    ##  [6346] 4.233058e-02 8.466645e-01 1.168561e-01 2.844060e-02 2.860076e-01
    ##  [6351] 2.018924e-07 2.925289e-02 1.133806e-01 2.209864e-01 2.345176e-02
    ##  [6356] 7.004141e-01 1.048379e-02 3.021381e-02 1.679912e-06 5.561966e-03
    ##  [6361] 1.102744e-06 5.579410e-05 8.235828e-01 3.204838e-01 2.561125e-03
    ##  [6366]          NaN 5.774691e-02          NaN 3.204838e-01 8.664211e-01
    ##  [6371] 4.299040e-08 2.946672e-02 1.351620e-01 5.524222e-06 1.858112e-01
    ##  [6376] 2.989939e-01 7.721914e-01 3.399497e-04 4.807060e-01 1.556353e-03
    ##  [6381] 2.126004e-08 3.500609e-12 2.953794e-01 5.010307e-01 2.931365e-02
    ##  [6386] 2.406438e-17 1.588825e-59 1.005556e-01 3.339241e-01 6.561589e-01
    ##  [6391] 5.799208e-02 8.803443e-02 2.275635e-01 4.291608e-01 2.761901e-01
    ##  [6396] 8.129209e-01 1.298279e-01 1.898995e-01 4.031788e-03 2.439145e-03
    ##  [6401] 5.530525e-16 3.790456e-03 1.440250e-06 2.025584e-01 6.439435e-04
    ##  [6406] 1.531077e-02 5.732021e-01 7.492348e-01 5.227597e-01 2.553794e-01
    ##  [6411] 2.516647e-01 8.258076e-01 8.757094e-06 6.798764e-02 5.875350e-08
    ##  [6416] 6.998409e-13 1.566457e-01 4.297393e-02 9.049378e-01 3.422312e-01
    ##  [6421] 4.775604e-16 1.735257e-01 8.582499e-03 7.914602e-03 2.735569e-04
    ##  [6426] 9.260349e-06 4.026888e-02 1.142368e-01 5.461016e-10 1.234660e-04
    ##  [6431] 2.222680e-16 9.609062e-08 2.104991e-04 9.760850e-01 6.202046e-01
    ##  [6436] 3.848580e-15 1.982982e-01 2.985206e-01 4.601948e-01 9.775742e-01
    ##  [6441] 1.460163e-03 6.554422e-04 4.442337e-02 4.325815e-01 2.693290e-02
    ##  [6446] 4.729190e-01          NaN 1.434830e-02 9.152331e-01 3.981438e-01
    ##  [6451] 5.567319e-03 1.123914e-01 4.419348e-02 8.076708e-02          NaN
    ##  [6456] 2.915214e-01 7.657639e-01          NaN 1.372991e-05 4.823486e-01
    ##  [6461] 3.487886e-01 1.977209e-01 9.635503e-01 2.001618e-01 6.701759e-01
    ##  [6466] 2.168632e-03 2.499593e-07 1.261850e-06 3.711970e-14 1.046659e-01
    ##  [6471]          NaN 1.640287e-03 1.743213e-02          NaN 5.034459e-07
    ##  [6476] 1.839885e-01 1.084769e-01 3.988212e-01 3.639515e-05 1.000170e-03
    ##  [6481] 3.274450e-13 5.235203e-01          NaN 9.251256e-02 3.610488e-02
    ##  [6486] 5.326057e-01 1.120359e-01 1.090905e-01 3.139302e-05 5.035319e-01
    ##  [6491] 6.550140e-01 2.285284e-01 4.580723e-07 8.799402e-04 3.191785e-03
    ##  [6496] 4.750533e-01 1.461316e-06 2.603895e-06 5.166060e-02 5.532480e-01
    ##  [6501] 4.143248e-01 2.741426e-01 3.357227e-05 5.307033e-01          NaN
    ##  [6506] 4.655098e-07 1.870803e-01 8.052878e-04 6.144476e-01 9.498227e-25
    ##  [6511] 6.335598e-03 2.647295e-03 2.455965e-09 2.085088e-03 6.513871e-01
    ##  [6516] 1.580545e-01 1.769047e-06 1.135853e-02 2.424261e-07 5.646989e-02
    ##  [6521] 5.878248e-02 3.338232e-01 5.823038e-02 6.119261e-02 3.902007e-02
    ##  [6526] 2.674654e-09 7.343617e-01 8.868577e-03 4.028697e-03 2.423413e-03
    ##  [6531] 4.395911e-01 6.520375e-11 2.047700e-01 2.401505e-03          NaN
    ##  [6536] 4.792247e-02 4.022615e-01 6.771905e-02 7.802373e-01 1.783354e-01
    ##  [6541] 7.265193e-01 9.989642e-01 4.319493e-01 1.056330e-06 8.623854e-01
    ##  [6546] 3.716717e-01          NaN 4.610983e-02          NaN 1.301610e-07
    ##  [6551] 3.444467e-03 2.250011e-10 1.163421e-03 1.691078e-01 6.685393e-01
    ##  [6556] 2.844027e-01 5.898338e-01 2.590379e-05 8.309472e-01 3.259571e-04
    ##  [6561]          NaN 6.163316e-06 3.045613e-02 6.073141e-02 1.937922e-07
    ##  [6566] 8.124182e-03          NaN 2.833409e-01 1.002460e-16 6.811839e-01
    ##  [6571]          NaN 5.383521e-03 7.758138e-01 1.145242e-01 1.005253e-01
    ##  [6576] 2.418968e-04 2.279739e-12 3.981900e-04 6.114699e-02 1.751878e-03
    ##  [6581] 1.212232e-02 2.444384e-01 1.715121e-11 7.508173e-06 4.106901e-01
    ##  [6586] 1.605357e-01 2.549626e-02          NaN 3.510540e-01 2.753603e-02
    ##  [6591] 8.187979e-01 2.159168e-01          NaN 4.667031e-02 5.176690e-01
    ##  [6596] 3.665607e-01 8.524915e-02 2.657538e-01 9.569124e-01 3.538979e-03
    ##  [6601] 2.881278e-01 4.789810e-01 3.203346e-01 7.469233e-02 4.133228e-02
    ##  [6606] 1.823574e-10 8.504014e-01 1.116165e-01 2.554884e-18 4.108502e-02
    ##  [6611] 4.261830e-04 5.181662e-02 5.629151e-01 5.048842e-05 2.160436e-03
    ##  [6616] 6.546838e-01 1.155355e-12 2.056953e-04 3.864540e-01 2.152718e-01
    ##  [6621] 1.713634e-01 4.524671e-25 1.658099e-01 2.156852e-10 6.443411e-02
    ##  [6626] 1.963965e-02 1.799089e-01 2.379710e-04 4.407684e-07 2.881516e-02
    ##  [6631] 9.736431e-53 7.237312e-03 6.814219e-01 1.013565e-01 9.588398e-02
    ##  [6636] 3.714412e-02 3.218561e-01 3.579495e-02 3.204838e-01 6.850124e-02
    ##  [6641] 6.755785e-02 2.950769e-01 7.028356e-03 5.455176e-03 5.906606e-07
    ##  [6646] 2.522520e-13 1.158098e-29 9.203600e-04 6.204553e-01          NaN
    ##  [6651] 5.115305e-03 1.056579e-05 9.359806e-01 4.691480e-01 1.759973e-01
    ##  [6656] 2.757538e-02 1.317218e-05 3.059281e-02 1.905442e-01 6.966547e-02
    ##  [6661] 5.878651e-01 2.104724e-01 2.639382e-12 5.116371e-06 2.200658e-02
    ##  [6666] 1.213677e-01 2.318499e-01 6.210927e-03 3.267790e-01 8.250828e-02
    ##  [6671] 1.424599e-04 9.592252e-03 6.884280e-01          NaN 9.211946e-01
    ##  [6676] 7.395268e-02 7.412332e-02 8.630792e-01 3.273359e-01 1.184074e-06
    ##  [6681] 8.630133e-01 4.958770e-03 1.929694e-04 1.313059e-01 4.703005e-01
    ##  [6686] 5.396097e-01 1.252422e-01 2.402327e-01 9.202033e-01 5.712960e-01
    ##  [6691] 8.344039e-01 2.938074e-06 4.918436e-01 6.192837e-03 3.578800e-03
    ##  [6696] 9.610277e-01 1.424539e-02 9.640253e-03 7.465833e-01 3.775540e-02
    ##  [6701] 5.519605e-01 3.204838e-01 1.068006e-13 7.868366e-01 2.020775e-02
    ##  [6706] 1.185223e-09 1.075186e-02 2.048917e-07 4.980997e-01 4.929535e-04
    ##  [6711] 3.346052e-01 4.442030e-09 6.569179e-02 1.604569e-20 1.814139e-03
    ##  [6716] 2.026722e-01 3.437158e-06          NaN 4.907824e-01 8.418390e-01
    ##  [6721] 2.181373e-03 2.105836e-02 6.434516e-01 1.548441e-07 1.079103e-01
    ##  [6726] 1.083714e-07 2.365256e-01 1.670047e-01 6.415267e-01 3.204838e-01
    ##  [6731] 5.185646e-01 1.364941e-02 1.448977e-02 9.860621e-07 7.084921e-03
    ##  [6736] 3.204838e-01 6.495186e-01 9.268881e-01 3.140011e-02 1.991928e-01
    ##  [6741] 9.835153e-01 1.172940e-01          NaN 5.443539e-12 3.273622e-02
    ##  [6746] 3.981916e-03 3.739857e-01 1.063653e-01 1.651590e-01 9.893137e-01
    ##  [6751] 3.067375e-02 1.685234e-03 3.337969e-19 6.195416e-05 1.706609e-02
    ##  [6756] 7.119853e-20 7.527488e-02 5.728060e-01 7.365913e-02 1.138742e-01
    ##  [6761] 1.308992e-04 5.261704e-04 1.016615e-26          NaN 3.531743e-01
    ##  [6766] 2.009368e-11 3.242563e-25 8.803459e-03 1.406339e-01 4.246216e-01
    ##  [6771]          NaN 7.517821e-02 7.574159e-03 2.578389e-02 6.803033e-02
    ##  [6776] 1.340226e-03 3.159259e-01 1.118502e-01 3.579597e-02 1.531161e-03
    ##  [6781] 4.236366e-01 8.443471e-01 1.451494e-16 5.943272e-01 7.057575e-01
    ##  [6786] 1.423701e-03 2.752208e-03 6.118087e-18 2.639703e-01 7.679135e-02
    ##  [6791] 7.699659e-01 5.552643e-01 4.103937e-06 7.797334e-02 1.052198e-03
    ##  [6796] 9.343875e-06 1.950538e-07 9.176360e-01 3.384443e-01          NaN
    ##  [6801] 6.829631e-01 3.409739e-01 1.717955e-01 4.298304e-01 8.089156e-01
    ##  [6806] 2.387661e-03 9.244317e-02 1.521732e-07 2.337766e-01 2.184639e-05
    ##  [6811] 3.808301e-01 4.051253e-03 1.368001e-02 1.922326e-01 1.891856e-01
    ##  [6816] 8.831907e-02 3.558936e-01 1.527139e-15 4.630290e-05 4.520620e-01
    ##  [6821] 1.077853e-14 4.381956e-02 8.703845e-01 9.404067e-01 2.974241e-04
    ##  [6826] 9.356661e-01 4.352817e-01 8.949217e-04 2.113128e-11 5.629655e-01
    ##  [6831] 6.282678e-01 3.934948e-01 6.130350e-02 3.204838e-01 4.116615e-01
    ##  [6836] 5.512929e-01 5.086045e-12 8.021917e-02 5.256080e-18 8.428911e-01
    ##  [6841] 6.962586e-01 5.529641e-14 3.878003e-04 3.652940e-04          NaN
    ##  [6846] 2.411580e-03 1.676072e-04 1.969482e-01 9.940360e-01 1.289259e-01
    ##  [6851] 7.384776e-01 9.902424e-03 3.503698e-02 2.967064e-01 5.327738e-01
    ##  [6856] 9.577753e-01 8.537180e-01 2.939850e-18 3.569469e-02 9.649134e-01
    ##  [6861] 1.834500e-01 1.097826e-03 9.882607e-06 4.016363e-01 1.429886e-01
    ##  [6866] 9.145031e-01 1.885908e-01 3.366831e-03 1.726794e-01 3.385368e-29
    ##  [6871] 2.906319e-08 2.001418e-04 1.093083e-01 5.452849e-04 7.157038e-06
    ##  [6876] 7.506759e-01 4.511172e-01 2.621257e-02 2.091707e-02 5.115546e-01
    ##  [6881] 1.565492e-04 1.352160e-04          NaN 8.597190e-01 4.753168e-01
    ##  [6886] 6.722743e-01 2.983000e-02 4.988657e-03 1.020922e-01 2.603127e-10
    ##  [6891] 2.245165e-03 9.601610e-01 8.044916e-03 1.200335e-02          NaN
    ##  [6896] 9.595066e-01 5.672657e-01 3.040800e-02 8.098493e-01 4.253911e-01
    ##  [6901] 2.431106e-02 1.340293e-13 1.308202e-05 2.707049e-02 7.309511e-01
    ##  [6906] 8.797795e-01 5.967217e-01 1.815532e-04 4.948383e-02 6.345258e-01
    ##  [6911] 2.553123e-01 5.063078e-01 1.120982e-02 4.032960e-10 1.100540e-03
    ##  [6916] 2.082786e-02 9.199705e-01 1.721770e-01 6.091688e-01 8.301164e-01
    ##  [6921] 5.264978e-01 7.730448e-02 9.078962e-01 5.068705e-01 1.218395e-04
    ##  [6926] 2.979868e-06 7.864523e-05 1.321815e-04 1.089671e-02 2.221820e-27
    ##  [6931] 2.991769e-03 2.040199e-07 1.676104e-05          NaN 5.398370e-02
    ##  [6936] 4.422358e-01 9.878810e-02 5.151329e-01 5.300929e-03 8.732743e-04
    ##  [6941] 9.676835e-01 5.147979e-03 1.854043e-01 9.575772e-01 1.289161e-02
    ##  [6946] 9.855777e-03 6.853703e-01 9.591096e-02 1.200015e-02 7.354366e-04
    ##  [6951] 4.826679e-01 2.394140e-03 7.475216e-02 5.920178e-01 4.172281e-03
    ##  [6956] 4.054692e-01 2.412997e-01 5.402171e-04 1.234296e-02 3.789380e-01
    ##  [6961] 5.113620e-02 8.563050e-07 7.975014e-03 1.463359e-05 5.249673e-04
    ##  [6966] 1.493437e-02 8.453848e-02 1.949270e-02 4.858932e-01 9.225981e-08
    ##  [6971] 1.335620e-01 7.288498e-01 4.291342e-01 1.299059e-01 9.407937e-04
    ##  [6976] 1.275781e-01 5.042431e-01 8.928162e-01 2.265155e-04 2.308672e-07
    ##  [6981] 6.368756e-01 4.640699e-03 6.574279e-03 9.199114e-01 1.502540e-08
    ##  [6986] 1.213024e-22 4.876015e-06 4.240362e-01 1.371310e-04 4.970718e-01
    ##  [6991] 4.225767e-05 3.864524e-01 3.650789e-01 4.760426e-10 2.074466e-02
    ##  [6996] 8.267373e-01 3.482227e-02 1.156756e-05 1.091705e-01 1.016317e-01
    ##  [7001] 5.007542e-02 6.666397e-05 4.634329e-03 1.926602e-01 1.036178e-02
    ##  [7006] 2.307222e-01 1.565292e-01 8.112797e-01 2.587752e-02 2.655635e-02
    ##  [7011] 1.120518e-03 1.106307e-01 9.560194e-01 2.985133e-01 3.014853e-03
    ##  [7016] 8.057466e-01 1.535383e-01 9.324914e-03 9.059061e-01 1.961506e-01
    ##  [7021] 1.259220e-16 5.314079e-01 2.399927e-08 8.360921e-02 4.304336e-09
    ##  [7026] 1.325284e-05 1.112851e-07 9.495855e-01 1.279440e-01 1.425452e-03
    ##  [7031] 6.740987e-02 1.891488e-01 4.943501e-01 2.039173e-01 2.601111e-10
    ##  [7036] 2.689857e-01 2.886850e-02 5.319064e-01 3.370534e-09 8.690035e-27
    ##  [7041] 1.940941e-01 3.559091e-02          NaN 1.389311e-02 6.009570e-39
    ##  [7046] 2.803732e-03 1.141870e-01 1.477657e-01 5.221345e-09 2.334256e-01
    ##  [7051] 7.668413e-01 9.961529e-01 6.463007e-01 1.503168e-01 1.780624e-12
    ##  [7056] 3.827737e-34          NaN 5.604419e-02 1.662325e-23 5.178320e-02
    ##  [7061] 7.463210e-15 2.239464e-01 9.542106e-06 3.414386e-01 9.292629e-09
    ##  [7066] 1.406142e-02 1.400126e-24 4.987390e-02 2.942290e-01 3.982215e-01
    ##  [7071] 9.851486e-01 3.204838e-01 4.721199e-01 2.561520e-02 1.903016e-03
    ##  [7076] 1.424638e-03 2.376054e-02 9.495693e-01 6.070257e-26          NaN
    ##  [7081] 2.913729e-02 6.733878e-06 1.665982e-01 3.915771e-01 4.410913e-01
    ##  [7086] 4.127289e-02          NaN 2.736772e-01 7.497859e-03 5.504492e-04
    ##  [7091] 2.749668e-01 3.032096e-01 9.243460e-15 1.363021e-04 3.204838e-01
    ##  [7096] 1.778282e-06 7.038420e-03 1.453177e-01 3.415555e-02 9.173633e-01
    ##  [7101] 6.483608e-02 9.238726e-02 4.242295e-03 7.351715e-12 3.204838e-01
    ##  [7106] 8.935737e-02 3.730397e-02 1.761128e-02 8.206405e-01 2.619512e-01
    ##  [7111] 4.373644e-04 9.107822e-01 1.254018e-02 2.732230e-01 9.823233e-02
    ##  [7116] 1.938869e-01 9.412698e-01 2.111425e-04 7.475094e-02 6.134445e-03
    ##  [7121] 1.589283e-15 9.563103e-01 5.720474e-01 4.582936e-01 7.987042e-02
    ##  [7126] 4.602550e-12 5.151772e-06 3.012798e-01 6.854461e-01 2.079264e-01
    ##  [7131] 3.186386e-05 7.259676e-06 8.657729e-01 2.689430e-01 7.001442e-01
    ##  [7136] 6.161689e-01 4.907799e-02 4.901096e-01 4.292992e-02 1.074900e-01
    ##  [7141] 5.778073e-01 3.777167e-01 1.061005e-09 4.584819e-07 1.537921e-03
    ##  [7146] 6.868648e-16 3.204838e-01 4.855339e-18 6.304237e-03 9.925129e-01
    ##  [7151] 9.233160e-12 9.919007e-02 6.051044e-01 9.530745e-01 8.506287e-09
    ##  [7156] 2.506540e-24 6.291064e-02 1.560424e-01 1.132489e-01 2.949423e-01
    ##  [7161] 9.993841e-01 2.319361e-01 1.287247e-19 2.317545e-01 4.054942e-06
    ##  [7166] 1.207171e-09 7.316659e-02 2.244017e-03 1.681945e-03 8.973840e-01
    ##  [7171] 9.828815e-01 4.710243e-02 1.659753e-01 5.963906e-01 5.455519e-14
    ##  [7176] 2.895260e-01 9.462961e-02 1.614194e-08 8.309869e-01 1.159193e-01
    ##  [7181] 2.455938e-23 8.301161e-01 2.977871e-02 6.398182e-07 3.308322e-03
    ##  [7186] 1.296492e-01 2.821862e-06 1.155579e-02 6.159779e-01 1.989618e-03
    ##  [7191] 4.921778e-01 1.566408e-01 6.714606e-02 7.884450e-01 1.876506e-01
    ##  [7196] 9.824874e-07 4.916349e-01 1.729883e-13 2.925159e-01 4.735981e-01
    ##  [7201] 3.149885e-02 2.168728e-07 7.743364e-01 8.444666e-01 1.111801e-08
    ##  [7206] 1.672023e-14 2.171033e-12 7.727172e-07 5.760091e-03 1.806177e-15
    ##  [7211] 1.829166e-02 1.049219e-07 3.784391e-06 1.654673e-04 8.642981e-03
    ##  [7216] 1.742239e-09 6.238958e-04 1.806195e-07 4.574286e-01 8.005302e-01
    ##  [7221] 2.858457e-01 1.457178e-01 1.752292e-10 1.004116e-07 2.495522e-03
    ##  [7226] 5.671459e-03 5.383433e-01 3.085136e-03 6.615820e-12 8.830812e-03
    ##  [7231] 1.049550e-01 3.935778e-02 2.800793e-01 1.195246e-02 1.073656e-02
    ##  [7236] 4.153092e-10 1.002464e-01 5.222794e-01 1.468693e-04 5.609591e-02
    ##  [7241] 1.107544e-01 2.991043e-09 1.683287e-17 1.519571e-01 1.263358e-02
    ##  [7246] 5.529727e-01 5.665799e-03 3.204838e-01 3.714834e-01 3.010585e-01
    ##  [7251] 3.204838e-01 5.096956e-01 4.933491e-01 2.470244e-09 1.675417e-01
    ##  [7256] 3.628962e-03 3.700638e-01 2.195092e-05 1.088840e-03 3.662041e-02
    ##  [7261] 9.083444e-02 3.699809e-02 1.242556e-01 8.719601e-01 1.990497e-13
    ##  [7266] 7.470574e-07 3.263527e-04 2.864616e-01 6.538213e-01 3.547484e-02
    ##  [7271] 3.635746e-01 5.825864e-02 1.541840e-01 9.000325e-01 8.214143e-01
    ##  [7276] 5.957559e-02 6.732838e-04 3.862050e-01 5.797322e-02 2.113328e-06
    ##  [7281] 7.274388e-01 4.747687e-01 2.682934e-01 3.755223e-02 8.242164e-01
    ##  [7286] 2.719676e-02 1.555091e-01 1.160890e-02 1.688861e-02 2.040696e-06
    ##  [7291] 7.172572e-01 6.424549e-02 2.928639e-07 6.821459e-05 2.418186e-01
    ##  [7296] 1.064562e-20          NaN 9.561268e-01 6.163594e-01 1.571688e-06
    ##  [7301] 7.422672e-01 7.638006e-01 1.215716e-25 7.758922e-05 2.344221e-01
    ##  [7306] 2.807966e-01 2.697845e-04 3.427687e-06 1.624641e-01 1.838451e-02
    ##  [7311] 2.225976e-01 8.145669e-06 1.002214e-03 3.541911e-03 9.816305e-21
    ##  [7316] 4.084375e-01 7.063452e-04 6.815546e-01 1.578407e-01 5.939002e-01
    ##  [7321] 5.482331e-02 9.588689e-02 1.721219e-10          NaN 8.024123e-01
    ##  [7326] 4.097384e-01 8.504729e-05 1.051383e-03 3.204838e-01 5.640124e-01
    ##  [7331] 9.553299e-01 5.721294e-09 5.661706e-02 7.547652e-01 1.771735e-01
    ##  [7336] 1.407743e-15 1.818693e-01 3.830035e-02 9.747545e-01 5.290284e-03
    ##  [7341] 2.459916e-01 8.127682e-06 9.459863e-01          NaN 4.927682e-01
    ##  [7346] 3.882321e-03 6.074789e-05 4.818071e-01 3.543264e-05 1.409066e-14
    ##  [7351] 9.461929e-01 3.466222e-03 1.560038e-09 2.804743e-22 3.204838e-01
    ##  [7356] 7.581089e-01 2.809937e-01 8.551013e-01 1.052584e-03 3.066914e-02
    ##  [7361] 5.514960e-04 3.320034e-01 1.893847e-04 9.325085e-07 1.802488e-01
    ##  [7366] 1.157590e-01 2.257633e-01 5.496752e-01 1.703808e-01 8.872639e-03
    ##  [7371] 7.402343e-01          NaN 5.523271e-02 2.162530e-01 1.323976e-03
    ##  [7376] 6.365282e-01 2.217433e-01 6.911497e-01 7.469122e-11 6.862958e-12
    ##  [7381] 6.557661e-01 1.829499e-02 4.388724e-01          NaN 2.342350e-03
    ##  [7386] 8.503353e-01 9.349307e-01 1.260642e-03 3.587559e-05 9.667579e-01
    ##  [7391] 7.084195e-02 2.697538e-03 3.114623e-01 4.008968e-02 9.992962e-02
    ##  [7396] 7.789991e-03          NaN 3.383813e-01 2.298359e-03 4.817415e-01
    ##  [7401]          NaN 3.204838e-01 1.307633e-07 5.222841e-03 1.378047e-01
    ##  [7406] 3.249144e-01 8.044807e-01 4.132045e-01 1.602017e-02 8.693628e-03
    ##  [7411] 2.060013e-02 3.962092e-03 1.134344e-02 6.169170e-01 1.498385e-01
    ##  [7416] 9.386580e-03 1.656700e-01 2.192633e-01 7.689448e-02 3.686318e-03
    ##  [7421] 8.072012e-12 1.850733e-23 4.094927e-01 1.593678e-02 1.171305e-02
    ##  [7426] 8.258330e-02 8.780129e-02 1.767315e-36 5.916160e-04 7.577158e-01
    ##  [7431] 4.653203e-01 5.648166e-01 2.059936e-12 8.047566e-01 4.608658e-02
    ##  [7436] 8.945733e-01 1.791779e-03 5.473066e-02 1.523249e-54 1.515268e-05
    ##  [7441] 9.156092e-01 9.646250e-01 1.095165e-01 4.514947e-01 9.196833e-01
    ##  [7446] 1.158180e-02 6.945007e-01 7.672605e-02 3.683386e-02 2.577383e-01
    ##  [7451] 2.690170e-03 2.069072e-08 6.825079e-04 1.190412e-04 2.528052e-01
    ##  [7456] 8.233580e-02 2.961577e-01 2.035818e-04 5.859602e-01 8.551357e-01
    ##  [7461] 4.478615e-01 3.288718e-24 3.217085e-01 7.839719e-02 1.769357e-02
    ##  [7466] 9.042968e-01 1.644891e-07 2.317331e-01 4.740258e-01 8.875521e-06
    ##  [7471] 2.082151e-04 5.444796e-01 4.628849e-03 8.521334e-02 3.365659e-02
    ##  [7476] 5.196266e-01 2.111387e-30 5.776757e-04 4.228862e-01          NaN
    ##  [7481] 2.124095e-02 9.564044e-01          NaN 7.923762e-12 2.786116e-02
    ##  [7486] 7.284591e-02 7.550800e-04 1.227158e-02 4.190549e-01 5.987792e-01
    ##  [7491] 7.556261e-02 7.178393e-01 5.684776e-01 3.748970e-01 5.037289e-01
    ##  [7496] 3.829683e-01 1.760825e-02 5.555696e-03          NaN 6.412078e-01
    ##  [7501] 7.265368e-02 4.546399e-02 1.730944e-03 2.882707e-01 7.016276e-01
    ##  [7506] 9.293657e-04 1.452771e-03 1.865425e-02 1.718655e-01 3.300970e-01
    ##  [7511] 8.317043e-02 7.361341e-02 2.359600e-01 1.237236e-09 3.629287e-01
    ##  [7516] 8.198910e-18 1.350358e-01 2.607711e-03 4.758524e-06 2.689647e-01
    ##  [7521] 8.355705e-01 4.665485e-01 5.770694e-01 1.026577e-03 5.957016e-02
    ##  [7526]          NaN 2.974012e-03 4.400509e-14          NaN 3.746539e-02
    ##  [7531]          NaN 2.161334e-02 2.397797e-01 9.237932e-04 1.211941e-01
    ##  [7536] 5.026859e-15 8.080108e-01 5.556310e-02          NaN 5.955480e-01
    ##  [7541] 7.083462e-01 8.040553e-02 1.268500e-01 5.826860e-03 2.580735e-01
    ##  [7546] 1.241968e-01 8.638032e-10 8.387553e-01 9.853695e-13 6.103158e-01
    ##  [7551]          NaN 1.594970e-02 2.549045e-01 2.214719e-10 2.894830e-01
    ##  [7556] 1.362410e-01 2.127011e-01 7.675980e-01 4.192799e-01 6.624687e-01
    ##  [7561]          NaN 2.330651e-02 1.809941e-01 3.937284e-01 1.719498e-01
    ##  [7566] 7.777005e-01 5.457796e-06 8.120790e-01 7.032713e-02 4.458911e-01
    ##  [7571] 7.417585e-01 6.555621e-02 8.663738e-02 2.810967e-12 1.403412e-16
    ##  [7576] 8.950560e-01 9.361059e-04 9.479154e-01 8.603493e-02 7.344841e-01
    ##  [7581] 6.353798e-01          NaN          NaN 9.633177e-01 8.016750e-01
    ##  [7586] 2.347865e-05 5.961360e-27          NaN 9.447358e-01 3.651161e-02
    ##  [7591] 9.043338e-06 7.194166e-01 4.653703e-01 4.463803e-03 3.412494e-01
    ##  [7596] 9.476079e-06 2.644137e-01 8.882488e-01 3.087672e-01 2.505403e-02
    ##  [7601]          NaN 2.857687e-01          NaN 1.133369e-01 2.010986e-02
    ##  [7606] 6.574464e-08 1.823696e-07 4.435235e-02 2.960450e-02 6.965467e-01
    ##  [7611] 2.992827e-06 3.289552e-02 9.771050e-02 3.744247e-04 7.119410e-02
    ##  [7616] 2.993480e-01 6.279627e-01 5.932346e-04 6.755538e-01 8.187368e-04
    ##  [7621] 2.049189e-01 9.952262e-01 1.046502e-02 8.164252e-02 2.481099e-01
    ##  [7626] 5.778676e-01 8.579787e-03 7.391258e-07 1.418244e-03 3.834519e-01
    ##  [7631] 5.313894e-01 5.477275e-01 3.176697e-01 9.152629e-02 6.083137e-02
    ##  [7636] 2.934505e-02 3.147431e-01 8.954245e-01 6.895008e-01 5.613466e-04
    ##  [7641]          NaN 4.529084e-01 2.962377e-02          NaN 2.816254e-12
    ##  [7646] 1.487043e-02 1.775098e-01 4.359461e-01 2.711098e-04 1.654641e-01
    ##  [7651] 1.143547e-01 6.146642e-24 1.212785e-02 1.181483e-04 2.908702e-25
    ##  [7656] 1.597842e-01 2.910715e-01 1.226050e-02          NaN 5.640982e-02
    ##  [7661] 2.315040e-02 5.866185e-01 6.543181e-03 1.976232e-19 1.504722e-09
    ##  [7666] 6.889549e-03 6.197910e-01 4.223017e-01 3.071493e-01 4.914817e-01
    ##  [7671] 1.627084e-01          NaN 4.511306e-04 5.051935e-03 1.079652e-01
    ##  [7676] 7.017424e-02 2.630641e-02 1.598344e-19 3.344466e-02 4.167592e-01
    ##  [7681] 2.299408e-01 7.376128e-01 9.997439e-09 2.879262e-19 5.396503e-02
    ##  [7686] 8.228810e-12          NaN 2.557927e-01 1.493289e-05 9.302723e-01
    ##  [7691] 1.306779e-02 2.223646e-04 7.295968e-05 4.730006e-05 5.576749e-01
    ##  [7696] 5.072215e-06 3.335791e-10          NaN 1.188630e-03 9.746891e-01
    ##  [7701]          NaN 3.460954e-05          NaN 9.161704e-02 1.519841e-06
    ##  [7706] 1.046467e-01 3.854662e-01 2.305380e-04 5.843799e-01 2.891061e-01
    ##  [7711] 1.613095e-11 4.864078e-05 8.215583e-03 3.047131e-17 2.719966e-01
    ##  [7716] 1.472146e-01 2.341441e-01 4.874784e-02 9.034455e-02 7.247837e-15
    ##  [7721] 8.322488e-01 8.500510e-01 8.233178e-03 3.264486e-01 1.297657e-01
    ##  [7726] 3.206711e-01 1.424419e-02 6.933660e-01 3.746543e-01 2.430098e-01
    ##  [7731] 3.412597e-01 3.196997e-01 2.469407e-01 6.810329e-03 2.247002e-06
    ##  [7736] 8.200909e-01 5.572482e-01 2.827670e-01 9.653821e-01 2.629118e-01
    ##  [7741] 1.613716e-01 2.118323e-13 8.161425e-21 3.682104e-02 9.451961e-03
    ##  [7746] 7.893881e-01 4.493974e-03 1.325919e-03 4.309920e-01 7.953698e-01
    ##  [7751] 1.040466e-03 2.147680e-03 6.300950e-01 1.145747e-02 2.521810e-01
    ##  [7756] 9.126963e-02 1.111886e-04 5.506916e-02 4.492812e-01 3.090410e-01
    ##  [7761] 1.781500e-01 1.832994e-01 3.076449e-01          NaN 9.622792e-01
    ##  [7766] 2.787460e-01 3.609953e-04 7.924835e-01 6.704757e-14 2.523965e-01
    ##  [7771] 2.519543e-29 9.271050e-31 7.078259e-03 7.979750e-01 6.030345e-03
    ##  [7776] 8.129640e-12 3.799641e-01 4.073745e-13 9.274865e-01 1.169005e-15
    ##  [7781] 7.210628e-01 3.996916e-01 2.004603e-05 2.442167e-08          NaN
    ##  [7786]          NaN 7.580797e-04 2.528535e-01 4.936487e-01 1.253944e-11
    ##  [7791] 4.714569e-01 5.277236e-06 1.602166e-01 2.183343e-02 5.220748e-01
    ##  [7796] 2.750437e-02 3.943075e-02 5.791942e-01 5.783326e-11          NaN
    ##  [7801] 1.495948e-01 2.528769e-01 5.286159e-01 8.704318e-01 9.736032e-05
    ##  [7806] 1.183453e-03 3.585339e-01 1.944600e-26 5.936046e-09 2.731176e-02
    ##  [7811] 2.650398e-01 6.264705e-15 3.766113e-01 4.914900e-02 4.228256e-02
    ##  [7816] 5.057545e-01 1.545603e-13 3.561148e-02 5.509177e-03 6.228233e-01
    ##  [7821] 8.850647e-03 1.426414e-04 5.696440e-01 6.330271e-11 9.625952e-09
    ##  [7826] 3.236994e-01 4.934679e-01 9.979533e-01 1.059112e-03 3.605501e-01
    ##  [7831] 1.561432e-01 1.017952e-04 1.600102e-04 5.377611e-01 8.653646e-02
    ##  [7836] 4.429848e-01 1.011787e-03 7.019521e-01 4.262570e-01 2.728203e-01
    ##  [7841] 4.418675e-01 1.979876e-07 1.957327e-06 3.319275e-14 3.088242e-02
    ##  [7846] 8.553020e-01 8.744429e-01 5.052950e-06 1.993468e-02 1.130001e-18
    ##  [7851] 3.108934e-04 2.382126e-01 5.329807e-02 3.901807e-01          NaN
    ##  [7856] 9.879157e-01 2.031261e-17 9.768891e-01 7.597996e-01 1.904290e-03
    ##  [7861] 6.993978e-01 1.852436e-01 5.358979e-01 5.436163e-01 4.649776e-02
    ##  [7866] 9.478128e-01 7.681504e-15 3.185562e-01 8.103520e-06 7.744830e-03
    ##  [7871] 3.362505e-07 7.710960e-01 1.449808e-01 2.694548e-06 5.884762e-02
    ##  [7876] 1.210566e-14 7.000215e-01 2.860716e-01 1.235107e-26 6.584029e-04
    ##  [7881] 6.980092e-17 1.005555e-01          NaN 3.573693e-02 8.147668e-01
    ##  [7886] 4.735375e-01 1.324394e-06          NaN 9.220263e-04 1.082047e-03
    ##  [7891] 4.682165e-26 8.015535e-01 7.154271e-01 7.749709e-12 1.320109e-02
    ##  [7896] 9.344305e-06 6.074655e-02 9.995837e-03 5.802589e-18 2.437158e-01
    ##  [7901] 1.235099e-02 9.073495e-02 2.615821e-01 2.218520e-03          NaN
    ##  [7906] 3.037175e-03 4.215567e-01 1.519710e-05 1.249669e-02 6.858296e-09
    ##  [7911] 9.494436e-01 2.014169e-11 4.168649e-01 3.183104e-01 3.923658e-02
    ##  [7916] 7.495954e-04 4.361597e-01 9.243811e-01 2.746505e-03 8.927413e-01
    ##  [7921] 5.231952e-01 3.473562e-05 6.829047e-01 4.850020e-02 1.635974e-01
    ##  [7926] 5.853433e-01 8.705232e-04 1.909062e-01 3.875042e-01 1.434265e-03
    ##  [7931] 1.971416e-01 3.448559e-01 1.563788e-01 3.193305e-21 1.501130e-01
    ##  [7936] 9.614305e-02 3.810039e-01 3.126607e-01 3.196077e-10 7.223510e-02
    ##  [7941] 1.431926e-01 9.331287e-01 3.017450e-02 2.053128e-01 4.136746e-05
    ##  [7946] 2.041419e-03 2.072324e-01 5.063163e-10 3.632047e-02 4.507224e-01
    ##  [7951] 3.186424e-02 5.361697e-04 7.596475e-01 5.512089e-01 1.014979e-04
    ##  [7956] 2.254800e-01 1.210707e-13 7.561364e-01 4.119472e-01 5.383658e-03
    ##  [7961] 1.530595e-04 8.931130e-03 2.630700e-02 3.842933e-03 3.204838e-01
    ##  [7966] 1.405349e-01 1.978541e-04 4.722830e-01 1.897141e-03 4.649335e-01
    ##  [7971] 2.998976e-01 2.231278e-02 1.877661e-04 2.008240e-01 3.260444e-01
    ##  [7976] 4.370682e-04 5.129404e-02 1.352683e-09 8.631181e-01 3.204838e-01
    ##  [7981] 4.928235e-01 2.969437e-21 3.247642e-01 1.278590e-01 5.441638e-01
    ##  [7986] 7.046740e-05 7.511689e-01 2.967437e-05 1.183301e-02 7.019153e-01
    ##  [7991] 6.996011e-01 3.204838e-01 8.785791e-01 2.465153e-02 7.514005e-01
    ##  [7996] 1.829089e-04 1.727161e-01 2.603657e-01 5.143997e-01 4.484722e-14
    ##  [8001] 8.503069e-04 6.980005e-08 3.129527e-01 8.751811e-01 1.305145e-02
    ##  [8006] 6.974895e-02 1.302227e-02 2.483262e-01 6.286664e-01 3.696413e-03
    ##  [8011] 2.906670e-01 9.498130e-01 5.919261e-01 3.204838e-01 1.780675e-04
    ##  [8016] 4.880213e-01 5.655284e-04 1.265029e-01 2.720467e-02 2.343153e-01
    ##  [8021] 2.101595e-03 6.335652e-07 7.165922e-07 3.628112e-01 2.698229e-03
    ##  [8026] 3.061271e-01 1.173106e-06 3.204838e-01 3.864087e-02          NaN
    ##  [8031] 2.759796e-01 4.943339e-03 4.486706e-18 3.123599e-01 4.183864e-01
    ##  [8036] 8.243881e-02 5.858313e-01 8.211996e-03 2.516868e-01 8.514534e-04
    ##  [8041] 5.290641e-01          NaN 1.003879e-04 6.258899e-05 9.200194e-02
    ##  [8046] 6.234510e-02 6.975106e-03 2.223593e-09 3.047933e-09 8.764070e-07
    ##  [8051] 4.031651e-01 5.801790e-03 2.295922e-01 9.385529e-01 4.541301e-03
    ##  [8056] 9.087916e-01 8.497298e-02 6.089031e-11 4.596595e-01 5.561532e-02
    ##  [8061] 2.006987e-01          NaN 1.745560e-03 3.606473e-11 7.504847e-03
    ##  [8066] 4.556285e-01 3.052555e-13 9.838000e-03 1.086956e-02 3.864154e-02
    ##  [8071] 2.192706e-01 2.297479e-09 6.676099e-01 2.709845e-06 2.380434e-01
    ##  [8076] 6.919674e-01 1.893265e-01          NaN 3.850508e-04 7.527491e-01
    ##  [8081]          NaN 1.189026e-01 7.861200e-03 6.849173e-06 6.657283e-04
    ##  [8086] 5.797692e-01 8.913435e-07 1.319667e-01 1.695239e-04 2.272243e-01
    ##  [8091] 1.789882e-01 1.777372e-01 2.847897e-02 7.873915e-01 2.811990e-01
    ##  [8096]          NaN 2.239018e-12 7.032853e-01 4.179206e-01 1.085042e-02
    ##  [8101] 4.969015e-03 7.650270e-02 5.440633e-18 2.879884e-06 5.852205e-03
    ##  [8106] 7.795626e-01 2.613194e-01 8.910000e-03 8.368983e-04 5.256323e-01
    ##  [8111] 7.462820e-02 4.645479e-04 1.280000e-03 4.759362e-09 3.204838e-01
    ##  [8116] 9.948769e-01 5.993338e-01 4.506218e-05 6.428891e-11 2.415906e-06
    ##  [8121] 3.746340e-05 3.204838e-01 5.971245e-04 7.162113e-02 4.599367e-01
    ##  [8126] 2.870348e-02 3.907751e-02 4.163673e-03 3.101538e-07 1.706983e-01
    ##  [8131] 7.973063e-01 1.905310e-09          NaN 2.474405e-08 5.243632e-01
    ##  [8136] 5.427348e-01 1.395034e-01 1.015594e-15 5.949682e-01 3.204838e-01
    ##  [8141] 1.786028e-01          NaN 1.365309e-02 4.696605e-01 1.773350e-01
    ##  [8146] 2.755994e-02 1.526688e-02 1.917451e-01 8.356106e-41 3.270333e-07
    ##  [8151] 2.217161e-10 4.769816e-12 3.371783e-05          NaN 3.280901e-01
    ##  [8156] 4.812247e-20 1.169706e-02 5.914690e-01 4.573361e-01 6.325436e-01
    ##  [8161] 2.986837e-01 4.679456e-01 1.329923e-01 1.114710e-01 1.609423e-04
    ##  [8166] 1.935982e-01          NaN 4.702597e-03 1.392406e-01 2.802204e-05
    ##  [8171] 3.204838e-01 4.774505e-13 9.187914e-01 1.285105e-01 2.934206e-01
    ##  [8176] 6.063024e-02 1.804169e-08 6.948510e-01 1.869378e-05 7.083478e-01
    ##  [8181] 1.774870e-01 2.615577e-17 2.459750e-01 4.888369e-01 2.792779e-02
    ##  [8186] 3.892158e-01 3.969345e-02 8.901491e-01 1.102560e-01 5.906814e-13
    ##  [8191] 7.880851e-03          NaN 8.611878e-06 1.711570e-01 7.124759e-03
    ##  [8196] 7.899517e-01 9.988780e-03 5.264812e-13 6.602611e-01 2.221514e-01
    ##  [8201] 1.907528e-34 4.124789e-06 5.773328e-01 5.979822e-01 9.349322e-01
    ##  [8206] 3.432664e-01 3.824156e-01 1.496703e-01 1.448863e-18 1.599181e-01
    ##  [8211] 3.529072e-03 9.143076e-02 1.421905e-05 2.422756e-01 2.121351e-01
    ##  [8216] 3.493427e-02 7.172295e-07 4.743660e-01 5.442251e-01 1.137626e-01
    ##  [8221] 8.002184e-02 6.252410e-02 6.729105e-02 1.744386e-03 2.413214e-04
    ##  [8226] 3.168681e-01 7.019209e-10 8.767168e-01 8.951998e-03          NaN
    ##  [8231] 3.473796e-03 1.680862e-07 9.806454e-01 9.297170e-01 5.307867e-01
    ##  [8236] 1.013192e-15 1.247988e-02 7.561778e-01 6.633360e-01 6.011252e-01
    ##  [8241] 4.222401e-22 6.458185e-11 4.440189e-12 9.957122e-05 1.216817e-02
    ##  [8246] 2.690864e-08 8.823888e-02 1.166780e-03 6.555993e-01          NaN
    ##  [8251] 6.446607e-04 2.669528e-02 3.231781e-04 6.313641e-02 3.919410e-07
    ##  [8256] 1.892489e-01 5.643540e-01 6.442368e-04 1.242562e-10 6.862797e-01
    ##  [8261] 6.182301e-01 4.711868e-01 2.854647e-01 6.777450e-02 4.287351e-08
    ##  [8266] 7.353665e-01 8.974945e-04 1.014551e-01 1.005464e-01 4.388022e-02
    ##  [8271] 5.685368e-02 3.129877e-06          NaN 4.043246e-01 8.995081e-09
    ##  [8276] 8.802109e-09 9.106167e-12 5.012511e-01 9.185725e-06 4.047301e-06
    ##  [8281] 8.907045e-02 8.497378e-01 2.935959e-01          NaN 3.295870e-07
    ##  [8286] 8.474002e-01 1.187978e-01 9.292440e-08 6.814794e-04 6.426987e-05
    ##  [8291] 1.150778e-01 3.632926e-01 6.356535e-08 1.261937e-03 2.739009e-02
    ##  [8296] 3.620335e-03 3.429093e-03 7.843776e-03 2.832170e-07 2.988596e-01
    ##  [8301] 1.696015e-01 1.829041e-01 5.668768e-01 3.296280e-01 2.764961e-05
    ##  [8306] 6.516001e-02 6.748533e-02 2.019131e-02 7.591716e-01 7.689178e-01
    ##  [8311] 2.012175e-09 1.371333e-04 1.484976e-02 2.881308e-01 1.480970e-29
    ##  [8316] 2.026418e-01 3.063686e-01 3.376303e-07 3.769384e-04 2.337940e-01
    ##  [8321] 1.079839e-01 2.144515e-05 6.939684e-04 2.048732e-01 4.205665e-01
    ##  [8326] 4.712578e-04 2.249743e-03 1.681469e-02 5.820338e-05 5.619036e-01
    ##  [8331] 2.554699e-02 1.077872e-13 1.116550e-01 2.087849e-04 4.855949e-04
    ##  [8336] 8.997420e-01 2.311555e-11 1.137116e-06 9.110373e-01 8.735793e-01
    ##  [8341] 9.386708e-01 1.183864e-02 6.496238e-01 2.199455e-02          NaN
    ##  [8346] 8.496289e-08 8.235551e-05          NaN 1.806355e-01          NaN
    ##  [8351] 3.204838e-01 6.330096e-01          NaN 3.093081e-02 2.282004e-02
    ##  [8356]          NaN 1.246634e-01          NaN 2.234291e-01          NaN
    ##  [8361] 2.083921e-02 2.585881e-04 3.752382e-01          NaN 1.733597e-03
    ##  [8366] 3.431434e-23 3.666785e-07 7.677049e-05 5.321418e-01 3.204838e-01
    ##  [8371] 1.058874e-01 7.832472e-01 9.195904e-01 1.307477e-01 6.633157e-02
    ##  [8376] 4.330770e-05 5.475304e-02 5.323302e-03 5.753929e-11 2.405374e-09
    ##  [8381] 6.609475e-01 1.465860e-01 2.335144e-02 3.728859e-02 2.601286e-01
    ##  [8386] 1.183974e-02 6.938149e-15 3.217553e-01 5.249212e-04 2.278476e-01
    ##  [8391] 9.020678e-01 8.566921e-01 1.329952e-01 2.619798e-01 4.030443e-19
    ##  [8396] 3.143327e-01 5.117230e-02 6.883031e-03 1.226967e-60 9.382519e-01
    ##  [8401] 3.756434e-01 1.181943e-02 1.164638e-02 6.671103e-02 3.456827e-03
    ##  [8406] 9.484732e-01 3.932159e-01 1.659864e-02 6.982875e-01 1.317992e-04
    ##  [8411] 2.928626e-01 1.015947e-17 2.314370e-01 6.305765e-01 3.285932e-01
    ##  [8416] 9.364248e-02 3.957992e-01 1.412888e-01 8.209872e-05 8.745994e-02
    ##  [8421] 6.128015e-01 5.860045e-01 1.575321e-04 4.762948e-02 6.483247e-01
    ##  [8426] 1.154808e-01 2.424433e-01 1.471562e-03 5.921662e-01 9.416039e-02
    ##  [8431] 1.954537e-01 1.300300e-02 4.278374e-01 2.313479e-01 2.432352e-01
    ##  [8436] 1.637080e-06 6.585425e-01 1.369633e-03 7.234109e-03 8.817580e-01
    ##  [8441] 9.782347e-01          NaN 7.094739e-04 3.204838e-01 2.682163e-06
    ##  [8446] 1.164366e-01 5.583578e-58 2.925948e-03 7.505821e-01 2.164860e-01
    ##  [8451] 5.789368e-08 1.444966e-03 4.187213e-04 5.524682e-02 2.529229e-02
    ##  [8456] 4.259631e-04          NaN 8.240711e-01 6.075485e-20 5.690871e-01
    ##  [8461] 2.911822e-01 7.221111e-01 1.270643e-01 8.801591e-01 9.173879e-05
    ##  [8466] 3.204838e-01 7.794816e-06 1.648032e-16 1.316090e-01 1.018439e-03
    ##  [8471] 3.699834e-02 4.173217e-01 6.250640e-16 5.902618e-02 7.793009e-02
    ##  [8476] 9.922759e-01 1.348567e-09 3.212137e-38 3.996044e-01 9.588873e-01
    ##  [8481] 8.260369e-01 8.629495e-02 4.034123e-04 1.445374e-01 6.061450e-16
    ##  [8486]          NaN 1.859780e-02 8.416039e-07 1.776906e-12          NaN
    ##  [8491] 6.671039e-01 2.652419e-02 7.845752e-01 1.364408e-01 3.941932e-01
    ##  [8496] 1.952773e-03 1.379271e-03 9.719389e-01 6.676518e-06 6.071502e-02
    ##  [8501] 1.134086e-01 1.129549e-03 2.771542e-02 2.855247e-01          NaN
    ##  [8506] 8.587617e-03 4.422912e-01 5.873655e-01 5.372388e-03 4.486828e-05
    ##  [8511]          NaN 1.267982e-21 8.487548e-04 5.043750e-16 7.344300e-01
    ##  [8516] 4.200283e-03 3.810938e-02 1.356495e-02          NaN 3.203365e-06
    ##  [8521] 7.876169e-01 2.290013e-01 6.744354e-02 9.605398e-08 1.384219e-01
    ##  [8526] 4.988662e-03 4.961631e-02 1.496309e-02 1.221394e-01 3.204838e-01
    ##  [8531] 1.092500e-18 3.397561e-01 3.178639e-09 3.578282e-02 2.473548e-04
    ##  [8536] 9.297934e-02          NaN 1.299776e-01 1.213743e-01 5.014563e-01
    ##  [8541] 2.433876e-05 2.903557e-03 8.496337e-06 1.790475e-01 1.063570e-02
    ##  [8546] 2.037002e-01 5.010748e-01 2.741534e-01 2.401511e-07 1.166828e-04
    ##  [8551] 2.615488e-01 5.658932e-01 4.060374e-14 3.925444e-01 8.945077e-02
    ##  [8556] 6.461111e-29 7.261835e-01 8.344476e-01 1.919191e-01 4.647396e-01
    ##  [8561] 3.399311e-01 5.385247e-03 2.272755e-36 4.867450e-20 9.043114e-01
    ##  [8566] 4.110911e-01 3.052491e-06 2.702288e-02 4.746639e-03 3.884379e-07
    ##  [8571] 1.840013e-02 4.988682e-03 1.903429e-04 5.578039e-01 1.829980e-27
    ##  [8576] 7.394122e-01 5.880681e-04 5.554975e-01 2.862989e-02 5.973730e-01
    ##  [8581] 1.483914e-01 6.505077e-02          NaN 2.724381e-01 4.178290e-02
    ##  [8586] 1.014017e-04 1.411731e-01 9.267477e-04 4.464898e-02 1.392361e-02
    ##  [8591] 7.479471e-06 2.471226e-02 6.166683e-01 8.024303e-07 2.062058e-01
    ##  [8596] 4.894513e-02 2.710374e-02 1.625659e-01 5.622116e-14 4.324850e-02
    ##  [8601] 1.696503e-11 1.334592e-01 8.691765e-03 9.977194e-03 6.556247e-01
    ##  [8606] 6.645708e-01 5.304199e-01 3.618706e-02 1.136809e-01 8.477051e-01
    ##  [8611] 7.137564e-01 2.321631e-01 1.170922e-01 1.464991e-01 2.628665e-01
    ##  [8616] 7.681354e-01 1.869956e-01 4.630718e-06 9.523070e-01 8.907182e-01
    ##  [8621] 2.363616e-03 1.074714e-01 9.423053e-16 4.851901e-02 2.191863e-01
    ##  [8626] 2.288979e-02          NaN 3.334570e-01          NaN 2.973621e-01
    ##  [8631] 1.332579e-11 2.588373e-18 5.730135e-08 5.050984e-01 4.707528e-30
    ##  [8636] 6.805080e-03 1.583526e-01          NaN 3.052277e-01 2.222423e-01
    ##  [8641] 6.330011e-02 5.092549e-03 4.374354e-02 4.294860e-20 3.505839e-01
    ##  [8646] 2.330298e-01 5.395513e-01 8.458914e-01 3.535554e-01 3.338199e-01
    ##  [8651] 6.088747e-05 2.262790e-01 1.169189e-01 2.963265e-01 2.790229e-02
    ##  [8656] 7.576999e-02 2.246407e-01 1.380449e-01 1.199158e-19 4.297716e-01
    ##  [8661] 8.339975e-05 2.575897e-04 3.026364e-07 8.155109e-01 1.529570e-04
    ##  [8666] 6.025858e-05 3.031397e-01 1.975811e-02 9.807647e-05 2.552863e-18
    ##  [8671] 4.039992e-01 2.924104e-01 3.856443e-02 3.193422e-01 1.640805e-01
    ##  [8676] 7.248996e-09 6.920386e-03 3.985112e-01 3.990700e-02 1.124459e-01
    ##  [8681] 5.739894e-01 1.342421e-01 2.912850e-03          NaN 5.114610e-01
    ##  [8686]          NaN 2.898782e-07 9.212039e-02 3.846026e-01          NaN
    ##  [8691] 8.149041e-06 1.208985e-02 3.042707e-01 5.918624e-01 5.589169e-03
    ##  [8696] 5.339802e-01 1.943724e-01 5.266623e-01 9.025423e-01 3.647822e-07
    ##  [8701] 6.314843e-01 4.888814e-01 1.171571e-06 3.204838e-01 5.343653e-01
    ##  [8706] 7.652865e-03 2.666065e-08 2.229334e-02 3.747971e-09 7.584773e-01
    ##  [8711] 7.246065e-02 3.050906e-04 3.742893e-07 7.901645e-17 6.773374e-01
    ##  [8716] 3.928408e-01 9.492925e-04 9.966541e-01 5.407346e-01 6.907106e-01
    ##  [8721] 5.286305e-02 1.888888e-03 2.826882e-02 1.211523e-04 1.010391e-21
    ##  [8726] 2.103038e-23 7.041962e-01 1.847895e-18 1.875905e-01 8.439834e-01
    ##  [8731] 3.577030e-62 6.257002e-01 4.936620e-02 4.513033e-02 6.384891e-01
    ##  [8736] 4.006802e-01 5.976908e-01 2.619172e-01 6.949016e-01 2.968991e-01
    ##  [8741] 3.486008e-02 9.155054e-01 8.272234e-01 4.213643e-01 4.479811e-01
    ##  [8746] 5.640265e-02 1.349044e-01          NaN 3.976534e-01 4.181574e-01
    ##  [8751] 2.662344e-03 2.979926e-07 1.634680e-01 6.384852e-01 1.940915e-01
    ##  [8756] 9.240571e-08 3.159888e-01 2.158941e-03 5.669943e-02 2.541634e-15
    ##  [8761] 3.729297e-09 1.171764e-03 7.036320e-04 9.984169e-01 1.007768e-01
    ##  [8766] 2.067071e-17 3.544178e-05 1.145234e-03 3.095006e-03 2.539760e-01
    ##  [8771] 2.747632e-09 9.964516e-04 1.360815e-11 3.820717e-05 4.713480e-01
    ##  [8776] 5.248533e-01 9.256371e-01 6.271890e-01 5.505185e-01 3.881688e-03
    ##  [8781] 5.961134e-01 7.305728e-06 5.978498e-02 1.607037e-02 5.281565e-02
    ##  [8786] 2.664179e-04 9.450844e-06          NaN          NaN 6.866499e-03
    ##  [8791] 8.198227e-01 5.123546e-04 6.023873e-04 3.856701e-02 6.171716e-01
    ##  [8796] 9.041578e-01 2.527099e-06 1.381305e-10 2.134700e-12 1.747163e-02
    ##  [8801] 1.912635e-01 2.657055e-01 6.950510e-01 4.568185e-01 7.606768e-22
    ##  [8806] 7.327725e-01 1.012126e-11 6.792580e-01 9.746848e-01 7.468653e-01
    ##  [8811] 2.506788e-05 6.956568e-07 2.895815e-01 9.104649e-04 7.667411e-01
    ##  [8816] 1.442616e-06 1.573382e-03 3.404854e-01 2.662621e-01 2.163075e-02
    ##  [8821] 1.649009e-01 5.513504e-01 2.129365e-06 4.357019e-02 1.589784e-03
    ##  [8826] 3.261276e-11 3.553949e-01 1.223033e-01 1.559233e-01 9.073688e-02
    ##  [8831] 5.309242e-01 1.075319e-08          NaN 1.644139e-01 1.492309e-01
    ##  [8836] 7.988070e-01 1.994096e-03 3.004719e-01 1.520753e-10 6.994828e-01
    ##  [8841] 7.019009e-01 1.525828e-10 8.827272e-01 5.364698e-02 1.368370e-02
    ##  [8846] 2.884991e-01 2.897198e-02 3.764298e-01 2.542418e-12 3.431762e-02
    ##  [8851] 5.372819e-05 7.215395e-10 3.278157e-02 1.026734e-03 4.000337e-03
    ##  [8856] 3.336840e-09 9.872479e-02 5.433306e-04 1.783042e-26 1.719634e-01
    ##  [8861] 2.500269e-02 3.045451e-03 1.213059e-10 3.869812e-04 8.588499e-02
    ##  [8866] 2.195652e-03 6.335384e-03 2.205163e-01 2.228889e-04 2.712743e-01
    ##  [8871] 1.518983e-01 2.311175e-02 9.154711e-01 7.534292e-10 2.261466e-01
    ##  [8876] 5.966181e-01 1.820184e-01 7.806805e-04 7.570329e-01 1.188087e-02
    ##  [8881] 2.778770e-03 2.838943e-01 1.698741e-05 2.815305e-01 3.168512e-01
    ##  [8886] 1.657124e-05 1.083366e-01 2.506125e-01 1.951629e-13 7.451579e-04
    ##  [8891] 3.866710e-02 3.968623e-01 6.322534e-05 1.594568e-01 2.119402e-01
    ##  [8896] 2.311624e-16 2.360038e-01 7.354028e-01 6.344290e-01 4.303573e-01
    ##  [8901] 9.970597e-07 3.569763e-04 6.763919e-01 3.341132e-05 3.478237e-03
    ##  [8906] 4.775466e-03 6.727007e-01 6.952533e-03 6.819523e-03 3.933688e-01
    ##  [8911] 1.698691e-12          NaN 9.158498e-01 2.625487e-02 8.648341e-19
    ##  [8916] 9.647518e-01 2.906501e-01 7.788255e-02 1.809372e-01 7.108050e-15
    ##  [8921] 5.499015e-03 3.767503e-02 2.937071e-05 5.193451e-03 9.739411e-03
    ##  [8926] 5.428230e-03 2.570462e-04          NaN 9.178167e-01 2.804061e-01
    ##  [8931] 9.284057e-02 5.532688e-01 7.057358e-01 3.876832e-03 3.703892e-05
    ##  [8936] 4.982262e-05 4.328713e-01 8.021624e-05 8.117927e-01 9.267499e-02
    ##  [8941] 5.061827e-01 4.164010e-02          NaN 2.492370e-12          NaN
    ##  [8946] 2.845700e-01 7.707697e-03 1.685752e-02 6.090670e-01 3.506023e-04
    ##  [8951] 4.399726e-01 1.597844e-01 3.247598e-02 6.968926e-03 4.642903e-01
    ##  [8956] 6.244629e-02 4.591390e-04 7.974641e-03 1.361116e-03 4.025246e-01
    ##  [8961] 4.707636e-01 4.383299e-06 2.177166e-01 2.065377e-06 7.340980e-01
    ##  [8966]          NaN 4.783996e-02          NaN          NaN 4.657225e-01
    ##  [8971] 2.992217e-02 9.891947e-01 7.269477e-01 2.109774e-02 9.453069e-02
    ##  [8976] 1.817092e-02 5.110774e-01 4.278731e-01          NaN 6.361781e-04
    ##  [8981] 8.234477e-01 9.668020e-04 4.585456e-04 4.362015e-02 3.236230e-06
    ##  [8986] 6.368416e-14 5.731473e-04 5.793350e-02 1.370892e-02 1.763733e-01
    ##  [8991] 9.404176e-14 9.493688e-01 6.301746e-05 6.447029e-01 3.439724e-01
    ##  [8996] 4.674002e-01 8.323160e-03 6.608564e-01 1.652938e-01 2.505749e-08
    ##  [9001] 2.116342e-22 2.516385e-01 8.531020e-01 1.585859e-01 3.260847e-02
    ##  [9006] 4.923761e-01 5.026043e-03 7.007713e-05 6.894600e-05 1.946051e-02
    ##  [9011] 7.320455e-01 2.437190e-01 5.297950e-02 3.931980e-17 1.201760e-01
    ##  [9016] 1.178736e-01 5.329623e-02 6.382585e-04 7.273123e-07 3.268851e-01
    ##  [9021] 6.334264e-02 7.107189e-01 4.852458e-03 5.326452e-02 7.013312e-01
    ##  [9026] 6.968120e-02          NaN 1.454103e-01 5.800296e-01 1.699216e-08
    ##  [9031] 2.647701e-09 2.941630e-01 3.205437e-08 6.865458e-01 1.121871e-13
    ##  [9036] 3.911100e-19 4.644943e-12 1.954323e-02 6.874566e-02 7.527581e-07
    ##  [9041] 9.684400e-01 9.313059e-01 3.958042e-01 1.085410e-01 3.204838e-01
    ##  [9046] 9.972702e-03 1.309036e-02 3.559350e-05 4.543036e-02 3.204838e-01
    ##  [9051] 7.136111e-04 2.909347e-02 4.188590e-01 5.300691e-01 4.356188e-01
    ##  [9056] 9.550153e-01 1.884406e-06 9.433260e-02 5.832746e-01 4.755350e-01
    ##  [9061] 3.012042e-01 9.513656e-01 8.097926e-01 4.811570e-01 7.951476e-09
    ##  [9066] 8.330463e-01 9.938103e-02 8.819230e-02 8.699509e-02 2.414174e-02
    ##  [9071] 1.579729e-04 4.267901e-06          NaN 2.924797e-02 3.557462e-03
    ##  [9076] 1.337779e-01 1.719197e-14          NaN 2.309540e-01 8.920674e-15
    ##  [9081] 7.579191e-01 3.259740e-01 6.667140e-01 2.377665e-01 4.860653e-03
    ##  [9086] 9.154098e-01 2.004544e-09 1.062832e-03 8.928519e-04 4.510736e-03
    ##  [9091] 3.626081e-03 1.664932e-08 2.063204e-01 4.493310e-10 4.597382e-01
    ##  [9096] 5.708963e-01 5.912931e-01 1.211718e-08 6.631877e-01 1.819025e-11
    ##  [9101] 5.138979e-05 3.591936e-02 3.749566e-16 6.008690e-01 5.295890e-01
    ##  [9106] 5.289842e-01 2.691206e-02 8.218093e-01 6.367674e-02 6.343728e-01
    ##  [9111] 3.164440e-03 1.940438e-02 3.061654e-01 5.943544e-02 1.770871e-01
    ##  [9116] 4.902456e-01 2.202667e-03 4.216443e-02 2.848340e-09 2.632439e-02
    ##  [9121] 2.274114e-03 5.098108e-01          NaN 4.915197e-02          NaN
    ##  [9126] 9.561130e-01 6.378562e-22 1.729218e-10 8.527744e-02 1.448970e-01
    ##  [9131] 2.655449e-01          NaN 4.183180e-01 2.824505e-03 5.171715e-03
    ##  [9136] 2.091526e-02 3.895173e-11 7.892604e-01 8.492136e-01 8.408072e-04
    ##  [9141] 5.221828e-01 3.219975e-01          NaN 2.325782e-11 2.773727e-02
    ##  [9146] 8.329721e-02 2.553517e-01 2.241023e-06 5.768141e-01 6.192753e-02
    ##  [9151] 9.020783e-01 2.377323e-09 1.132505e-02 4.806551e-01 2.529580e-04
    ##  [9156] 3.204838e-01 5.418891e-01 5.542837e-02 8.925868e-04 4.575082e-02
    ##  [9161] 5.741087e-03 1.246362e-01 3.204838e-01 2.543554e-03 3.699732e-02
    ##  [9166] 5.231111e-01 4.857540e-01 1.015639e-01 8.010582e-03 8.380559e-01
    ##  [9171]          NaN 1.363348e-01 2.528849e-02 5.960384e-01 5.828112e-01
    ##  [9176] 3.757553e-04 2.369665e-01 1.460766e-01 3.186931e-05 6.775640e-01
    ##  [9181] 6.568184e-01 5.943126e-29 4.246076e-01 1.445818e-05 2.026927e-07
    ##  [9186] 1.540291e-04 3.160660e-14 8.018544e-04 1.964672e-01 8.064731e-01
    ##  [9191] 2.125758e-03 9.549348e-02 1.068955e-03 4.047779e-02          NaN
    ##  [9196]          NaN 1.410760e-15 2.517515e-05 1.850904e-03 1.019245e-05
    ##  [9201] 2.740399e-01 2.132100e-02 2.120140e-01 5.584037e-01 2.755130e-04
    ##  [9206] 4.413776e-01 4.910511e-03 7.548101e-01 4.148287e-05 6.408296e-01
    ##  [9211] 5.758591e-02 1.544385e-04 1.897355e-03 3.158392e-01 1.028322e-01
    ##  [9216] 6.022901e-01 6.872857e-01 1.123972e-03 1.273950e-01 9.957344e-08
    ##  [9221] 1.019052e-09          NaN 7.172272e-03 3.548804e-11 2.413708e-01
    ##  [9226] 1.792518e-12 1.690748e-04 1.515462e-02 7.188499e-02 3.047926e-09
    ##  [9231] 2.334868e-29 9.102999e-01 3.204838e-01 1.064310e-01 3.603445e-08
    ##  [9236] 1.162167e-02 4.805539e-40 9.127081e-01 1.224839e-01 4.633624e-01
    ##  [9241] 9.153463e-03 7.077800e-03 5.362903e-01 1.466157e-01 3.546227e-03
    ##  [9246] 5.016592e-11 4.602762e-02 6.044437e-02 1.660443e-01 1.761984e-01
    ##  [9251] 5.934554e-01 3.310732e-05 3.101532e-02 6.364491e-01 5.901700e-27
    ##  [9256] 9.755614e-12 7.547209e-23 7.319384e-01 3.432816e-10 7.213783e-01
    ##  [9261] 7.621753e-09 2.087599e-03 3.799287e-01 7.851048e-01 8.983660e-01
    ##  [9266] 3.917816e-03 9.059050e-01 7.637221e-01 5.479237e-04 9.888727e-01
    ##  [9271] 1.244099e-08 4.027110e-01 3.578128e-02 8.874402e-03 9.695176e-02
    ##  [9276] 4.629212e-01 6.705420e-01 3.654355e-05 4.477938e-08 3.300146e-01
    ##  [9281] 5.369611e-01 5.576043e-01 3.784744e-01 3.979463e-03 1.210683e-02
    ##  [9286] 5.664055e-01 8.618040e-01 2.782699e-04 6.233154e-04 4.782698e-01
    ##  [9291] 8.783218e-02 7.169581e-01 4.695243e-01 3.160884e-06 7.046452e-02
    ##  [9296] 1.674722e-03 1.948950e-05 8.087703e-02 5.806014e-07 2.152720e-01
    ##  [9301] 6.775600e-01 3.578159e-04 7.967513e-02 2.474775e-01 2.303432e-19
    ##  [9306] 9.484500e-01 5.389772e-02 8.140898e-01 3.861210e-05 5.658225e-01
    ##  [9311] 3.885852e-14 3.657767e-02 6.217708e-03 7.451262e-01 4.612230e-01
    ##  [9316] 1.713721e-02 2.799184e-01 1.992456e-01 7.909214e-01 2.411721e-16
    ##  [9321] 4.785733e-01 5.971091e-01 5.122407e-02 3.681310e-09 5.624125e-02
    ##  [9326] 4.914731e-01 2.730931e-05 4.439200e-01 1.282002e-01 3.695006e-06
    ##  [9331] 1.638125e-01 6.843534e-01 4.133448e-02 2.444375e-04 1.924628e-04
    ##  [9336] 1.663087e-03 1.469070e-01 8.469636e-12 4.932269e-03 3.884304e-01
    ##  [9341] 1.760223e-05 8.889427e-01 1.752939e-02 5.424680e-01          NaN
    ##  [9346] 2.863907e-03 3.629403e-01 3.028926e-04 6.565307e-01 1.786924e-18
    ##  [9351] 6.777454e-04 6.727102e-05 1.626391e-02 1.662381e-01 5.497108e-01
    ##  [9356] 1.224457e-02 1.425528e-02 8.095375e-01 6.353877e-01 4.674279e-01
    ##  [9361] 2.894665e-04 4.102661e-17 8.116413e-04 7.211261e-24 1.300405e-01
    ##  [9366] 8.387194e-01 5.235444e-02 7.980770e-09 8.088494e-01 3.774194e-21
    ##  [9371] 2.913186e-03 1.749772e-04 2.913906e-01 2.423445e-02 5.370682e-35
    ##  [9376] 3.884053e-02 5.667444e-02 1.978378e-01 8.567377e-14 1.093452e-02
    ##  [9381] 3.899076e-06 6.695470e-02 4.544754e-01 2.019208e-03 1.469594e-01
    ##  [9386] 5.583262e-05 2.581151e-09 1.098103e-02 2.203315e-02 8.892828e-03
    ##  [9391] 7.653945e-01 1.444532e-01 3.366095e-02 1.405082e-02 2.891455e-02
    ##  [9396] 1.186043e-09 6.734486e-14 6.741662e-02 7.611521e-01 2.501765e-01
    ##  [9401] 2.385765e-03 9.556673e-01 2.199146e-03 1.835621e-06 9.141156e-01
    ##  [9406] 3.530634e-14 8.311693e-08 3.959632e-15 3.242554e-05 9.527178e-02
    ##  [9411] 5.938394e-01 1.699312e-01 7.359362e-04 4.033957e-02 8.595441e-01
    ##  [9416] 1.867981e-03 4.370439e-01 9.572282e-01 6.890593e-04 1.408127e-03
    ##  [9421] 8.561110e-04 2.594701e-01 1.422692e-03 2.911034e-01 3.204838e-01
    ##  [9426] 8.496015e-23 1.168001e-05 2.536726e-01 3.278726e-01 3.718601e-56
    ##  [9431] 2.729953e-01 2.271646e-02 2.765566e-01 2.031412e-05 7.818778e-01
    ##  [9436] 6.961485e-10 1.006254e-03 5.193397e-01 3.291066e-03 1.517273e-01
    ##  [9441] 4.163865e-01 7.310012e-01 9.947648e-01 3.861290e-03 1.441717e-01
    ##  [9446] 1.627127e-02 5.172622e-01 9.228461e-01          NaN 1.287014e-01
    ##  [9451] 4.325305e-01 1.915978e-07 2.258078e-03 5.964875e-06 2.734260e-04
    ##  [9456] 3.037891e-01 2.815203e-01 4.500375e-02 2.523701e-01 1.657881e-06
    ##  [9461] 3.380548e-02 3.204838e-01 1.046502e-07 3.343157e-44 9.160675e-03
    ##  [9466] 4.032638e-02 7.589277e-01 1.356740e-01 6.163321e-01 6.626295e-02
    ##  [9471]          NaN 1.169577e-07 5.481586e-01 7.388917e-04 2.654232e-09
    ##  [9476] 6.289333e-01 2.736856e-06 5.327486e-06 8.761893e-01 3.629487e-01
    ##  [9481] 7.738885e-06 6.816108e-01 6.650710e-03 3.924612e-01 9.288193e-01
    ##  [9486] 2.611529e-04 2.035381e-02 5.894493e-03 1.703810e-05 4.241948e-01
    ##  [9491] 3.177274e-03 5.439180e-03 1.360940e-01 1.334513e-07 2.037456e-03
    ##  [9496] 2.037223e-03 2.740129e-04 4.892409e-04 4.430962e-01 3.195794e-01
    ##  [9501] 5.159624e-01 7.183291e-01 2.918142e-01 8.540093e-01 6.365230e-01
    ##  [9506] 9.669402e-01 3.039661e-01 4.129212e-04 2.882929e-01 5.977579e-01
    ##  [9511] 1.175541e-06 1.837395e-01 4.296310e-03 1.130092e-02 1.616268e-01
    ##  [9516] 7.683892e-03 4.801646e-02 7.729180e-01 1.089760e-10 9.105560e-03
    ##  [9521] 8.529213e-05          NaN 4.511357e-01 4.904406e-03 3.852848e-02
    ##  [9526] 4.251975e-04 7.093631e-02          NaN          NaN 7.100604e-02
    ##  [9531] 7.592245e-01 5.079065e-01 7.188652e-01 7.187825e-07 1.684627e-01
    ##  [9536] 2.433461e-01 2.699402e-20 2.754061e-01 2.438821e-01 3.256549e-01
    ##  [9541] 2.450098e-01 7.395643e-01 4.830293e-01 7.626645e-04 3.090607e-12
    ##  [9546] 8.104277e-03 1.263595e-01 2.873076e-11 2.166521e-05 8.195564e-01
    ##  [9551] 5.900317e-01 5.193485e-01 7.252655e-03 9.622054e-02 7.366705e-02
    ##  [9556] 3.842609e-02 1.597907e-01 1.263875e-12 6.851673e-01 2.582901e-09
    ##  [9561] 1.978557e-01 2.497796e-01 2.858227e-01 3.818639e-04 2.055365e-03
    ##  [9566] 1.965662e-02 2.253322e-02 6.822584e-02 7.992677e-01 2.741894e-08
    ##  [9571] 9.171461e-01 9.483731e-01 1.775736e-02 1.469271e-02 6.505125e-05
    ##  [9576] 1.061407e-03 6.202457e-11 7.117684e-01 1.208384e-12 1.398763e-03
    ##  [9581] 4.363782e-01 2.995696e-02 2.969437e-05 5.844875e-03 9.811747e-01
    ##  [9586] 2.760589e-31 4.745666e-20 9.136224e-01 4.191910e-01 7.071462e-05
    ##  [9591] 4.320848e-02 1.108360e-02 1.812562e-02 3.425799e-01 5.263793e-03
    ##  [9596] 7.640392e-02 7.406228e-01 2.815052e-02 2.717309e-02 1.261698e-01
    ##  [9601] 5.035873e-07 4.767358e-12 1.284501e-01 1.110597e-02 6.412933e-01
    ##  [9606] 2.032908e-31 2.667556e-03 8.360550e-16 9.340235e-05 6.593241e-01
    ##  [9611] 2.410699e-06 6.300431e-04 1.685151e-01 7.987512e-01 1.020631e-01
    ##  [9616] 3.578496e-03 3.529804e-01 6.823360e-01 1.757715e-07 9.468172e-01
    ##  [9621] 1.944701e-03 8.470228e-04 1.653239e-01 8.983007e-05 1.116958e-17
    ##  [9626] 4.285480e-03 2.460593e-04 5.087677e-23 8.584941e-01 2.659227e-03
    ##  [9631] 2.389293e-01 1.600060e-17 5.503618e-01 8.835206e-01 9.793764e-01
    ##  [9636] 2.894630e-01 4.213149e-02 8.020957e-01 2.467600e-01 6.733531e-02
    ##  [9641] 2.608843e-02 6.973217e-02 7.198959e-04 5.615754e-01 7.694018e-04
    ##  [9646] 1.767282e-01 5.613487e-01 6.597245e-01 3.770655e-02 6.180815e-02
    ##  [9651] 1.597150e-06 1.255510e-02 3.476709e-03 3.267524e-04 6.819451e-01
    ##  [9656] 1.237268e-01 8.832109e-01 4.552788e-04 1.915642e-04 9.130559e-02
    ##  [9661] 2.937525e-01 8.580518e-02 4.015104e-03 2.146798e-05 5.739354e-01
    ##  [9666] 9.098734e-02 9.346521e-02 1.821298e-02 8.172197e-01 2.066729e-10
    ##  [9671] 7.192856e-06 3.931830e-02 4.205594e-02 8.671273e-01 7.231373e-03
    ##  [9676] 2.850937e-01 2.532539e-03 2.048546e-08 4.876406e-07 2.873747e-01
    ##  [9681] 6.002023e-01 4.452587e-01 5.787035e-04 2.295957e-04 1.728631e-02
    ##  [9686] 8.362032e-01 9.098633e-01 1.367331e-01 2.685418e-01 4.087904e-04
    ##  [9691] 8.481486e-01 7.432409e-03 6.540469e-08 6.746620e-01 9.408666e-02
    ##  [9696] 1.271498e-02 2.112754e-12 8.956053e-02 5.664683e-01 5.295891e-03
    ##  [9701] 7.463515e-04 7.102867e-01 9.302263e-02 9.843729e-03 1.475183e-01
    ##  [9706] 2.264896e-03 1.810040e-01 2.854639e-07 4.911389e-01 5.130666e-02
    ##  [9711] 2.507699e-08 7.724916e-01 3.933219e-07 8.133795e-01 7.997060e-02
    ##  [9716] 9.305450e-02 6.696223e-02 5.396406e-01 1.579131e-01 2.048666e-01
    ##  [9721] 7.654116e-05 9.393075e-02 2.860149e-04 1.920373e-17 7.532742e-02
    ##  [9726] 2.694728e-01 1.961306e-01 3.532759e-01 4.332595e-02 1.988410e-01
    ##  [9731] 5.861957e-01 5.673518e-01 5.094154e-09 1.077956e-02 5.928140e-10
    ##  [9736] 3.195293e-01 6.614875e-04 2.549937e-11 9.618387e-01 1.511960e-01
    ##  [9741] 3.589932e-02 7.872728e-01 1.173623e-01 4.643563e-22 6.812331e-03
    ##  [9746] 4.731867e-02 8.818680e-01 8.298679e-02 1.935373e-04 9.113597e-01
    ##  [9751] 6.008392e-03 5.276759e-01 2.930963e-01 7.117508e-04 1.079838e-01
    ##  [9756] 8.790853e-01 8.919604e-01 6.771081e-04 8.136947e-07 6.108753e-01
    ##  [9761] 1.619534e-02 1.052682e-19 2.283014e-07          NaN 4.219185e-02
    ##  [9766] 2.331949e-01 4.044367e-03 1.943833e-02 2.126825e-01 1.894773e-01
    ##  [9771] 9.446236e-01 1.349306e-01 1.954838e-05 6.113224e-01 4.186441e-01
    ##  [9776]          NaN 2.238021e-06 3.358431e-02 1.735419e-01 1.217808e-02
    ##  [9781] 5.981432e-02 6.720353e-01 5.266597e-01 1.051911e-01 1.062246e-04
    ##  [9786] 1.646780e-08 6.226126e-01 2.548190e-01 3.223215e-01 1.385635e-04
    ##  [9791] 5.045006e-01 8.505857e-02 7.021528e-01 7.309355e-01 2.340455e-02
    ##  [9796] 4.909455e-01 1.742939e-07 1.684751e-01 8.900047e-09 5.610728e-09
    ##  [9801] 2.407257e-08 1.257688e-01 6.825991e-07 7.628487e-01 3.702088e-01
    ##  [9806] 2.186954e-05 4.430503e-01 5.223945e-02 7.642001e-01 1.768313e-02
    ##  [9811] 1.448057e-09 2.074143e-01 9.523553e-02 1.544632e-13 2.784101e-02
    ##  [9816] 4.088626e-03 1.733739e-01 3.699646e-01 6.023536e-01 1.780987e-03
    ##  [9821] 2.697789e-01 3.204838e-01 5.459673e-03 5.278902e-03 6.811691e-12
    ##  [9826] 5.378724e-01 6.770027e-06 4.235547e-13 6.183784e-03 2.897726e-02
    ##  [9831] 4.752297e-01 3.622959e-02 2.904903e-03 6.137888e-05 2.020292e-02
    ##  [9836] 5.271501e-01 1.568538e-02 9.887055e-03 2.407187e-03 8.089359e-01
    ##  [9841] 1.198618e-01 4.700292e-01 4.650009e-01 5.290463e-01 6.558043e-05
    ##  [9846] 3.328961e-06 4.213729e-01 7.674002e-06 7.013972e-08 1.765686e-02
    ##  [9851] 2.890796e-05 4.956285e-25 4.538350e-03 3.236118e-04 3.270646e-07
    ##  [9856] 9.334523e-01 7.967940e-01 1.862052e-03 4.816892e-01 3.155472e-01
    ##  [9861] 6.502395e-01 9.475406e-02 6.038947e-01 3.746605e-01 3.823975e-02
    ##  [9866] 5.832614e-01          NaN 5.804369e-01 5.344360e-01 6.811039e-01
    ##  [9871] 7.176459e-03 1.548830e-01 1.623246e-01          NaN 2.457268e-04
    ##  [9876] 1.603088e-06          NaN 3.592998e-01 3.129401e-01 4.247735e-04
    ##  [9881] 5.662843e-02 1.726135e-02 5.237871e-01 1.023736e-01 4.483548e-02
    ##  [9886] 9.549502e-01 5.554118e-01 3.253700e-01 6.004437e-04 3.713074e-04
    ##  [9891] 2.567063e-01 8.924140e-02 1.568776e-02 8.296111e-02 2.746817e-01
    ##  [9896] 1.365528e-02 2.668482e-10 7.575548e-05 2.565545e-01 3.095346e-02
    ##  [9901] 2.721628e-02 5.437182e-01 3.552925e-01 1.249165e-01 2.475602e-01
    ##  [9906] 3.556151e-02 3.811022e-01 1.830005e-02 9.240175e-06 7.493849e-02
    ##  [9911] 9.537133e-01 1.531965e-01 1.966776e-01 9.805801e-02 4.092460e-01
    ##  [9916] 1.538422e-01 6.890259e-01 8.159271e-02 4.846144e-04 4.693508e-03
    ##  [9921] 4.541159e-02 8.733018e-01 9.389873e-04 3.602911e-01 2.261464e-01
    ##  [9926] 1.204902e-06 3.647580e-06 7.517116e-01 8.321860e-01 2.818054e-02
    ##  [9931] 2.081860e-03 7.303192e-02 8.032014e-02 1.465881e-02 3.204838e-01
    ##  [9936] 3.629226e-01 1.256211e-01 7.177051e-01 6.607413e-01 2.439748e-01
    ##  [9941] 6.333890e-01 5.093889e-01 5.365136e-01 8.913877e-02 9.789668e-01
    ##  [9946] 7.549456e-01 6.519251e-01 4.802345e-01 2.805548e-02 2.496844e-04
    ##  [9951] 1.289518e-01 6.800113e-01 1.009699e-01 3.983829e-01 8.645322e-01
    ##  [9956] 2.360926e-02 1.192027e-01 3.772007e-06 9.194401e-02 1.403487e-01
    ##  [9961] 4.564967e-02 7.220788e-10 3.731313e-01 4.614821e-02 8.164799e-30
    ##  [9966] 4.502405e-01 2.384264e-01 8.713539e-02 4.493808e-02 6.403810e-01
    ##  [9971] 2.873144e-06 8.753173e-14 4.755701e-04 9.835506e-01 1.433025e-10
    ##  [9976] 2.225426e-04 9.556335e-01 6.589405e-01 4.813850e-01 7.097258e-11
    ##  [9981] 4.370427e-01 9.865315e-01 7.195395e-01 2.947751e-05 6.163081e-03
    ##  [9986] 1.422049e-05 8.773151e-02 7.762063e-03 2.409925e-03 1.008128e-03
    ##  [9991] 2.103665e-02 4.665863e-01 1.467176e-01 6.030680e-01 5.553898e-01
    ##  [9996] 4.863453e-06 5.657780e-06 2.938601e-01 8.136026e-17 6.984626e-01
    ## [10001] 5.767091e-03 8.800564e-01 1.026315e-01 3.373512e-01 8.428432e-01
    ## [10006]          NaN 4.955518e-01 7.126094e-01 4.728036e-02 3.366655e-08
    ## [10011] 4.325460e-02 7.308777e-01 5.734079e-02 6.259924e-01 1.064557e-02
    ## [10016] 8.443199e-19 1.889724e-02 9.096802e-03 3.929456e-10 4.456689e-01
    ## [10021] 4.909561e-01 5.479730e-01 5.448220e-06 9.899264e-02 1.340531e-01
    ## [10026] 6.101304e-01 1.485090e-01 1.003434e-06 3.764933e-02 3.048040e-01
    ## [10031] 5.315260e-01 7.287194e-01 1.193707e-01 5.102767e-01 3.592631e-10
    ## [10036] 1.031745e-04 9.561943e-01 2.094850e-01 7.779947e-03 7.293069e-02
    ## [10041] 5.613859e-01 6.109336e-02 8.210769e-01 9.424052e-01 9.523706e-01
    ## [10046] 2.932133e-03 3.092869e-05 7.452441e-01 6.579217e-01 2.961561e-01
    ## [10051] 4.492876e-01 1.389249e-08 2.047574e-01 3.629332e-04 1.804711e-01
    ## [10056] 8.198183e-03 4.221830e-01 1.033933e-01 8.647466e-11 2.150945e-01
    ## [10061] 2.980651e-01 1.050587e-01 3.094497e-05 2.219197e-10 1.559036e-01
    ## [10066] 5.382621e-01 7.473223e-01 1.577770e-03 3.413393e-02 7.750659e-01
    ## [10071] 5.264155e-01 9.426469e-01 5.608951e-08 3.045865e-01 4.250179e-02
    ## [10076] 2.441721e-05 1.922880e-03 4.549356e-02 5.095212e-03 2.158646e-02
    ## [10081] 4.935968e-02 8.571586e-01 8.405632e-01 1.570254e-01 5.144578e-01
    ## [10086] 1.664785e-14 7.834931e-01 1.278271e-07 4.387724e-01 2.324558e-10
    ## [10091] 3.431667e-03 1.070870e-05 8.780285e-01 1.864578e-01 5.508344e-02
    ## [10096] 1.349032e-03 3.054707e-05 1.735403e-05 2.635518e-10 4.842002e-03
    ## [10101] 4.555721e-18 9.707930e-06 2.914502e-01 2.597385e-01 2.238697e-02
    ## [10106] 5.869448e-06 7.368798e-01 1.739341e-03 8.329185e-12 3.698507e-19
    ## [10111] 3.630064e-05 5.836207e-01 2.191491e-04 3.738214e-01 1.685670e-06
    ## [10116] 5.094975e-01 1.125145e-01          NaN 7.080592e-01 7.846633e-09
    ## [10121] 7.117814e-02 7.660791e-01 5.375369e-04 3.400591e-07 1.852038e-04
    ## [10126] 4.708827e-01 7.730041e-04 2.849171e-01 1.573672e-01 3.204838e-01
    ## [10131] 6.334065e-01 1.143871e-01 5.114477e-02 5.712492e-01 3.159607e-01
    ## [10136] 7.006942e-01 9.823143e-08 3.200088e-01 4.484409e-01 1.102696e-04
    ## [10141] 6.493874e-01 3.123662e-10 6.486561e-15 5.012676e-02 6.552658e-04
    ## [10146] 7.462407e-01 9.118800e-02 4.981397e-01 2.410897e-01 3.217985e-01
    ## [10151] 8.532239e-02 7.115902e-01 8.185914e-01          NaN 5.178073e-04
    ## [10156] 1.804880e-01 1.932815e-02 1.098651e-01 3.558306e-01 2.384389e-07
    ## [10161] 8.090864e-01 8.638006e-01 3.711664e-01 1.352454e-24 3.464314e-02
    ## [10166] 1.039485e-01 1.706737e-01 3.921146e-03 3.580605e-02 8.346912e-01
    ## [10171] 4.925176e-01 5.039381e-01 7.803387e-01 5.684230e-01 1.393725e-07
    ## [10176] 1.647320e-12 3.914621e-02 3.424739e-04 2.800749e-01 6.447942e-01
    ## [10181] 1.220632e-03 1.196118e-01 3.802955e-11 8.085590e-01 6.612371e-01
    ## [10186] 8.411866e-01 3.204838e-01 3.572641e-04          NaN 8.640301e-04
    ## [10191] 3.817717e-04 7.222793e-01 1.350276e-05 4.503323e-03 5.890091e-01
    ## [10196] 8.974627e-01          NaN 3.512892e-11 1.518006e-06 3.580405e-03
    ## [10201] 3.773883e-02 8.686456e-01 2.339933e-02 3.115464e-01 8.793256e-01
    ## [10206] 1.010950e-01 1.164167e-01 3.209992e-02 1.432376e-06 2.250176e-05
    ## [10211] 6.228059e-01 1.974271e-02 8.326823e-01 6.529078e-01 4.805074e-01
    ## [10216] 7.960261e-02 2.842175e-02 1.750319e-02 6.292238e-01 2.048501e-02
    ## [10221] 3.667803e-01 2.357430e-04 1.044599e-14 3.579225e-01 6.745184e-01
    ## [10226] 1.343237e-14 6.342520e-01 1.266155e-02 1.315287e-08 1.944060e-02
    ## [10231] 7.113776e-01 4.406202e-04 5.082963e-03 3.181143e-01 4.456859e-01
    ## [10236] 9.993824e-01 1.072968e-01 1.143301e-02 3.356536e-36 1.994193e-01
    ## [10241] 2.685109e-01 6.427201e-01 1.655528e-03 1.230984e-04 9.039408e-01
    ## [10246] 8.555985e-02 2.109077e-01 2.064845e-11 8.789214e-01 1.758542e-08
    ## [10251] 1.015082e-04 2.984935e-03 5.140125e-01 3.204838e-01 7.221167e-01
    ## [10256] 7.517497e-01 2.906293e-01 7.248815e-01 1.762662e-01 6.211872e-01
    ## [10261] 4.011107e-02 3.186862e-01 9.678461e-03 3.984448e-01 3.010318e-01
    ## [10266] 2.193433e-15 2.462815e-04 4.930617e-02 2.721514e-09 3.067316e-01
    ## [10271] 9.279486e-01 6.544767e-01 4.189865e-01 1.770526e-01 1.611787e-01
    ## [10276] 3.043350e-01 2.529250e-01 2.250237e-03 1.135870e-03 1.209000e-04
    ## [10281] 5.156265e-02 6.414623e-01 9.051052e-01 2.134189e-01 1.540917e-02
    ## [10286] 8.536736e-01          NaN          NaN 6.475754e-01 8.279481e-02
    ## [10291]          NaN          NaN 3.204838e-01 5.307815e-01 5.695411e-01
    ## [10296] 1.419740e-07 5.278819e-01 3.132363e-03 4.255076e-10 2.511472e-02
    ## [10301] 6.167353e-02 8.493465e-01 2.077940e-06 3.926413e-04 5.612685e-02
    ## [10306] 1.471721e-02 4.772137e-01 2.504713e-08 1.005794e-02 3.972028e-01
    ## [10311] 5.882558e-01 8.912767e-01 2.238500e-03 1.108592e-04 6.802894e-01
    ## [10316] 2.593512e-01 5.183909e-01 3.809406e-01 3.464988e-01 8.527205e-01
    ## [10321] 5.102611e-02 5.861146e-01 4.334742e-04 8.267902e-02 9.550703e-02
    ## [10326] 1.099143e-03 3.376689e-14 1.844236e-01 8.344806e-01 2.161091e-03
    ## [10331] 1.152192e-02 2.844753e-01 1.398964e-03 7.290335e-01 5.006926e-02
    ## [10336] 8.366708e-05 7.989923e-02 9.925963e-01 3.168525e-01 1.832906e-03
    ## [10341] 5.069955e-01 2.420923e-01 3.633758e-02 2.472673e-02 2.060837e-02
    ## [10346] 4.410253e-02 3.172347e-14 4.074587e-05 2.167723e-04 8.930557e-01
    ## [10351] 6.832310e-01 8.818575e-01 4.235523e-01 2.613657e-02 1.021604e-01
    ## [10356] 3.413024e-01 2.076796e-01 6.627253e-03 2.150169e-05 3.641103e-09
    ## [10361] 5.396053e-01 5.122844e-01 2.094511e-03 1.193827e-01 8.196642e-01
    ## [10366] 8.696687e-03 1.927966e-03 9.920836e-04 2.535733e-05 3.236501e-04
    ## [10371] 6.947286e-04 3.739782e-10 1.647131e-02 2.864245e-02 1.356907e-38
    ## [10376] 3.183177e-01 9.944163e-16 3.412006e-12 3.052087e-01 1.630224e-01
    ## [10381] 7.968026e-01 3.055771e-04 8.305793e-01 9.277060e-02 1.271940e-01
    ## [10386] 6.034296e-04 4.605616e-01          NaN 2.447136e-01 7.456624e-01
    ## [10391] 3.281108e-05 1.373903e-02 7.624307e-01 1.050763e-01 6.160029e-09
    ## [10396] 3.699183e-01 8.170967e-01 2.923745e-04 8.078423e-01 8.143703e-01
    ## [10401] 9.988654e-01 2.662223e-01 2.104773e-06 8.697283e-02 1.387882e-10
    ## [10406] 7.760098e-01 1.077139e-03 3.736158e-01 1.144629e-02 9.077449e-07
    ## [10411] 7.880025e-01 1.483289e-02 4.783475e-01 2.622213e-02 4.877897e-01
    ## [10416] 7.136578e-01 4.885595e-02 1.154581e-01 1.397260e-15 6.480879e-01
    ## [10421] 6.571313e-01 5.942301e-01 1.891711e-02 1.642180e-01          NaN
    ## [10426] 1.983884e-09 1.871423e-02 3.838615e-01 2.935565e-02 4.912357e-01
    ## [10431] 4.391996e-03 4.546323e-02          NaN 4.067554e-02 6.637044e-08
    ## [10436] 2.064460e-05 6.805549e-06 8.975693e-01 1.324342e-01 4.697567e-01
    ## [10441] 6.680243e-03          NaN 2.158224e-01 2.124612e-01 9.790560e-09
    ## [10446] 8.826339e-01 8.863899e-01 5.006436e-01 1.146376e-01          NaN
    ## [10451] 7.412177e-01          NaN 5.320432e-02 2.012478e-04 1.967940e-05
    ## [10456] 1.424372e-12 7.547534e-01 2.097846e-01 2.810613e-09          NaN
    ## [10461] 6.239241e-01 9.496988e-01 3.186030e-01 7.102526e-01 2.519664e-06
    ## [10466] 2.173002e-41 6.443779e-01 5.320667e-02 8.837358e-01 5.909569e-02
    ## [10471] 3.856126e-01 6.407285e-14 7.855530e-01          NaN 1.260731e-04
    ## [10476] 1.846219e-03 2.173444e-02 8.913984e-01 5.216358e-01 8.035321e-06
    ## [10481] 1.944140e-01 8.146913e-04          NaN 5.054698e-02          NaN
    ## [10486] 4.036328e-01          NaN 2.263985e-01          NaN 3.204838e-01
    ## [10491] 9.645474e-02 3.820222e-01 8.205742e-11 7.018359e-01 3.204838e-01
    ## [10496] 1.300646e-10 1.346567e-02 2.794789e-02 3.445940e-01 9.016474e-01
    ## [10501] 3.149264e-04 2.085262e-09 8.231181e-01 3.591936e-01 2.584666e-01
    ## [10506] 6.672839e-01 9.938806e-06          NaN 2.764882e-01 2.117300e-01
    ## [10511] 4.755439e-01          NaN 2.348139e-01 2.996023e-05 1.740867e-12
    ## [10516]          NaN 5.260537e-01 7.178247e-01 3.035435e-01 1.133391e-02
    ## [10521] 8.060451e-01 8.963491e-04 4.076242e-02 3.365981e-01 3.635514e-19
    ## [10526] 7.028091e-01 1.317843e-04 6.815401e-01 8.218814e-01 3.445911e-05
    ## [10531] 4.397497e-01 4.943170e-01 4.513138e-01 9.771951e-18 5.559818e-01
    ## [10536] 3.253194e-01 2.664126e-02 2.754551e-02 7.524212e-22 7.722720e-15
    ## [10541] 3.228213e-01 1.882976e-02 7.209544e-01 4.251035e-02 4.533772e-01
    ## [10546] 8.905510e-01 2.292989e-13 3.790669e-01 7.230070e-03 4.962329e-01
    ## [10551] 4.868378e-02 6.186421e-01 6.287852e-02 9.624452e-01 5.394726e-01
    ## [10556] 3.993498e-06 1.010986e-01 1.933625e-01 2.852392e-01 3.227925e-01
    ## [10561] 1.945975e-11          NaN 3.301345e-03 5.031036e-02 6.313242e-04
    ## [10566] 2.512136e-06 3.196097e-01 1.977178e-01 6.510628e-09          NaN
    ## [10571] 2.070052e-01 7.232512e-02 9.326675e-03 7.304018e-01 2.202227e-01
    ## [10576] 6.340276e-02 6.816162e-01 5.815959e-04 1.893706e-02 1.069850e-08
    ## [10581] 4.180163e-02 4.962183e-02 2.559011e-02 9.689929e-01          NaN
    ## [10586] 2.564146e-04 9.768812e-03 7.138587e-01 9.756555e-01          NaN
    ## [10591] 4.291810e-01 4.043449e-01 1.473272e-01 2.139888e-02 6.122967e-03
    ## [10596] 3.454807e-02 4.950531e-01 9.180054e-01 1.269438e-01 3.930210e-10
    ## [10601] 2.272723e-01 9.834157e-02 4.850353e-01 9.404045e-01 8.141634e-01
    ## [10606] 1.102922e-02 5.075693e-04 3.204838e-01 1.131040e-01 1.111478e-01
    ## [10611] 1.235924e-13 3.730416e-10 9.855273e-01 4.441634e-02 7.924723e-01
    ## [10616] 6.015524e-04          NaN 8.869965e-03 6.902600e-01 1.865403e-01
    ## [10621] 1.538541e-01 9.675238e-24 2.887173e-06 3.708650e-01          NaN
    ## [10626]          NaN 6.977638e-03 4.050615e-01 2.756456e-01 5.106715e-01
    ## [10631] 1.545235e-14 3.204838e-01 6.423724e-05 5.238107e-01 3.204838e-01
    ## [10636]          NaN          NaN 4.601446e-01 8.351796e-01 5.104570e-01
    ## [10641] 3.334624e-02 4.998002e-01 2.198703e-01 9.195545e-01 1.977392e-05
    ## [10646] 8.506287e-01 1.844991e-05 3.424623e-14 4.340613e-01          NaN
    ## [10651] 1.339394e-02 1.316481e-02 3.462729e-01 9.289609e-01          NaN
    ## [10656] 2.899579e-05 1.372183e-02 2.107352e-01          NaN 3.625321e-01
    ## [10661] 2.569115e-01 4.535411e-01 1.731025e-01 1.588789e-01 1.568995e-01
    ## [10666] 3.981691e-02 4.420641e-01 1.138555e-01 1.097767e-04 6.117990e-01
    ## [10671] 1.559749e-01 1.425315e-01 3.661782e-11 7.569584e-01 3.097415e-02
    ## [10676] 9.608956e-01 7.179884e-01 2.283985e-01 3.631521e-04          NaN
    ## [10681] 5.777953e-01 6.603429e-01 7.964365e-03 1.792830e-03 5.344258e-01
    ## [10686] 7.834330e-03 3.983247e-03 1.174145e-03 1.618147e-02          NaN
    ## [10691] 5.589487e-01 3.859385e-05 1.704560e-03 6.155219e-01 7.204844e-01
    ## [10696] 3.131612e-02 3.281611e-02 3.770728e-02 4.048067e-01 1.405192e-01
    ## [10701] 4.051008e-01 7.766578e-04 2.356877e-01 3.058553e-01 9.513985e-06
    ## [10706] 7.641427e-01 9.297795e-01 9.612375e-01 2.971272e-01 1.852997e-01
    ## [10711] 6.867288e-02 1.293040e-03 1.482041e-02 3.068187e-02 9.286822e-01
    ## [10716]          NaN 1.986528e-03 9.712249e-01 7.210943e-02 7.729672e-01
    ## [10721] 9.258483e-01          NaN 9.570303e-01 1.527354e-01 2.594415e-02
    ## [10726] 1.895144e-01          NaN          NaN 1.722271e-03 1.679658e-01
    ## [10731] 3.022804e-01 4.368094e-01 9.782821e-01 2.892164e-01 7.304373e-04
    ## [10736] 3.174566e-01 3.469897e-03          NaN 6.893214e-06 3.495070e-01
    ## [10741] 1.308672e-01 1.579710e-01 2.163061e-04 5.515165e-01 2.549507e-01
    ## [10746] 6.494982e-01 9.001094e-06          NaN 1.706828e-02 2.044048e-01
    ## [10751] 6.265765e-01 1.934284e-06 2.500146e-01 1.384585e-11 7.358307e-01
    ## [10756] 1.435242e-08 1.009544e-01 1.916989e-05 5.546279e-01 3.371878e-06
    ## [10761] 1.883347e-01 7.161659e-01 6.716354e-01 4.939366e-01 7.441894e-01
    ## [10766] 1.188238e-02 1.246622e-01 4.553910e-01 6.064835e-02 2.174570e-04
    ## [10771] 7.653958e-01 9.790128e-01 8.476389e-02 8.509832e-02 5.046168e-01
    ## [10776] 3.204838e-01 2.724493e-01 6.751996e-01 6.280944e-03 1.342282e-01
    ## [10781] 3.068494e-01 1.127983e-05 9.407650e-01 1.012695e-05          NaN
    ## [10786] 2.328213e-01 1.888427e-06 3.054089e-04 6.908732e-01 6.351654e-01
    ## [10791] 1.209943e-04 6.000573e-01 1.157688e-05 8.742670e-01 3.992765e-01
    ## [10796] 9.767606e-01          NaN 8.156921e-01 8.348345e-09 9.859798e-02
    ## [10801] 5.973298e-01 1.044866e-02 1.031827e-01 8.162013e-01 7.459952e-01
    ## [10806] 2.277804e-09          NaN 9.969388e-02          NaN 2.630575e-01
    ## [10811] 4.318457e-04 1.677237e-02 3.204838e-01 7.886719e-01 1.052613e-01
    ## [10816]          NaN 3.204838e-01 7.918023e-02 1.328336e-04 2.067902e-04
    ## [10821]          NaN 3.398906e-05 5.533205e-01 3.180240e-05          NaN
    ## [10826] 5.628090e-04 7.470257e-05 6.904255e-01 1.679379e-01 9.272057e-01
    ## [10831] 6.810138e-01 3.291410e-03 1.050880e-01 8.814934e-02 1.789178e-01
    ## [10836] 6.647218e-01 1.797995e-01 8.838654e-01 3.355264e-01 2.925257e-18
    ## [10841] 3.773932e-03 1.160675e-03 4.728138e-05 5.871287e-02 9.006011e-01
    ## [10846]          NaN 1.800770e-03 8.952588e-01 7.701599e-01          NaN
    ## [10851] 7.672365e-01 5.466679e-03 1.609248e-04 6.106511e-08 3.091601e-03
    ## [10856] 6.021535e-04 1.095924e-01 2.992665e-01          NaN 8.014098e-04
    ## [10861] 2.692509e-01 7.283985e-01 1.325591e-01 5.146376e-04 5.181693e-02
    ## [10866] 1.640447e-01 4.272204e-17 6.289901e-05 7.365631e-02 3.204838e-01
    ## [10871] 6.741829e-01 7.311165e-01 7.877045e-04          NaN 8.933306e-01
    ## [10876] 8.013649e-02 7.450532e-37 3.204838e-01 5.197441e-01 2.357543e-08
    ## [10881] 4.579590e-04 3.072724e-01 1.367381e-10 7.841368e-01 5.173575e-02
    ## [10886] 2.562711e-05 3.811119e-01 2.299635e-02 2.068340e-10 1.700193e-01
    ## [10891] 4.186623e-04 5.223696e-12 3.470284e-01 4.325245e-01 1.731999e-01
    ## [10896] 7.857321e-01 7.591109e-01 2.057673e-04 4.500275e-02 3.799875e-04
    ## [10901] 1.994413e-05 5.120679e-01 4.353629e-06 1.032369e-06 3.883932e-02
    ## [10906] 6.627028e-01 1.305394e-02 5.357999e-01 4.282755e-02 2.051423e-10
    ## [10911] 4.506973e-02 6.466354e-01          NaN 1.586154e-01 8.051908e-01
    ## [10916] 1.601395e-01 3.511279e-01 5.574018e-01 4.348479e-01 2.700193e-17
    ## [10921] 1.492361e-31 1.643812e-02 3.820770e-01 1.722381e-02 1.547296e-01
    ## [10926] 1.403599e-01 5.201873e-02 3.153078e-01 8.129089e-09 3.204838e-01
    ## [10931] 1.231199e-03 7.879271e-01 2.778891e-01 8.301450e-01 6.053672e-01
    ## [10936] 1.843898e-02          NaN          NaN 8.380785e-04 9.450626e-04
    ## [10941] 5.172724e-13 5.335472e-01 1.269957e-01 1.785527e-01 1.093449e-03
    ## [10946] 2.340953e-01 3.029618e-01 5.728974e-05 9.259212e-01 1.746661e-02
    ## [10951] 1.492006e-01          NaN 6.742902e-02 9.086045e-01 3.824263e-04
    ## [10956] 2.800551e-01 6.568893e-04 6.769602e-03 2.298582e-01 3.238097e-02
    ## [10961]          NaN 8.805879e-01 3.204838e-01 2.298643e-01 9.552213e-01
    ## [10966] 8.720213e-01 1.543401e-01          NaN 7.474722e-09 3.109381e-01
    ## [10971] 2.310767e-05          NaN 1.131626e-01 2.825427e-05 3.022137e-01
    ## [10976]          NaN 2.570362e-20 1.472394e-02 3.185199e-07 1.954149e-01
    ## [10981] 1.568929e-01 5.147305e-02 3.525062e-01 3.204838e-01 9.158168e-01
    ## [10986] 3.204838e-01 4.461164e-01 9.098559e-01 1.405109e-02 2.449985e-01
    ## [10991] 2.695427e-07 4.647211e-14 9.988777e-01 6.364714e-01 8.188041e-01
    ## [10996] 9.849837e-03 2.458586e-01 9.819565e-05 3.899451e-01 7.717540e-01
    ## [11001] 8.520608e-04          NaN 8.318141e-01 3.841329e-01 2.496247e-06
    ## [11006] 2.192357e-01 1.213284e-02 2.881709e-01 9.294411e-01 1.454782e-07
    ## [11011] 3.353228e-10 6.896296e-15 3.921784e-10 2.497003e-01 2.913187e-01
    ## [11016]          NaN 8.102150e-01          NaN 3.265881e-01 5.701962e-01
    ## [11021] 1.643094e-01 3.810052e-07 4.567421e-01 9.974406e-01 4.141177e-02
    ## [11026] 4.465835e-01 6.395833e-02          NaN 6.474012e-02 9.348470e-03
    ## [11031] 2.077699e-21          NaN 5.511204e-01 1.534463e-06 2.770482e-01
    ## [11036] 5.487118e-12 7.878869e-12 1.270682e-05 4.131348e-02 8.427571e-01
    ## [11041] 7.461742e-01 1.128933e-01 1.983353e-01 7.822144e-04 3.517836e-02
    ## [11046] 3.848220e-01 2.849426e-02          NaN 4.641665e-01 3.684751e-06
    ## [11051] 6.096432e-04 2.923881e-01 2.641437e-01 5.058194e-01          NaN
    ## [11056] 3.204838e-01 4.113543e-02          NaN 4.092277e-01 5.566964e-02
    ## [11061] 5.900131e-01 8.370215e-01 8.802603e-01 1.562133e-01 1.171163e-01
    ## [11066] 8.151305e-01 1.711390e-03 9.412988e-08 2.340242e-02 2.903179e-01
    ## [11071] 2.975230e-17 9.752128e-36 2.877588e-01 1.440422e-01 5.430051e-01
    ## [11076] 2.828746e-01 1.454830e-15 3.125150e-01 9.401855e-01 9.888641e-01
    ## [11081] 4.175794e-12 8.708184e-01 7.384977e-01 1.162722e-02 3.834521e-02
    ## [11086] 7.782756e-01          NaN 1.201758e-05 4.090503e-02 7.892820e-08
    ## [11091]          NaN 7.285498e-01 1.437555e-01          NaN 1.244407e-01
    ## [11096] 3.924799e-01 9.774552e-02 6.535390e-01 8.579359e-01 3.624847e-01
    ## [11101] 3.557573e-01 9.896821e-07 6.475315e-01 8.997333e-01          NaN
    ## [11106] 1.811596e-05 1.091220e-01 2.650628e-07 3.750295e-01 4.446581e-01
    ## [11111] 1.330812e-05          NaN 1.435809e-03          NaN 2.099696e-20
    ## [11116] 9.297375e-01 8.141816e-02 4.136713e-01 2.987980e-01 6.846845e-01
    ## [11121] 3.844451e-01 2.829622e-01 8.045801e-01 4.101028e-01 4.697863e-04
    ## [11126] 4.346854e-01 3.204838e-01          NaN 4.505855e-01 7.553488e-05
    ## [11131] 3.474746e-05 1.513660e-04 6.752137e-05 9.314374e-15 4.949536e-01
    ## [11136] 1.779683e-05 1.762938e-04 2.318550e-01 3.204838e-01 1.759892e-01
    ## [11141] 6.841860e-12 1.357980e-01 4.504847e-15 1.789999e-01 2.109850e-04
    ## [11146] 6.665229e-01 6.659477e-01 3.123532e-01 5.350534e-07 1.743706e-01
    ## [11151] 4.550360e-01 7.319987e-21 6.698744e-01 5.657975e-03          NaN
    ## [11156] 4.862647e-01 6.587683e-01 2.850399e-01 7.023470e-01 9.172198e-01
    ## [11161] 4.527022e-01 1.463039e-10 1.253866e-01 2.827932e-02 7.557283e-06
    ## [11166] 9.032519e-01 3.504880e-04 2.591930e-15 7.056052e-01 5.145897e-02
    ## [11171] 1.098409e-17 6.312484e-01 8.379227e-01 2.450085e-02 6.226149e-02
    ## [11176] 1.797268e-01 4.408820e-01 2.881323e-05 1.926691e-05 7.874806e-02
    ## [11181] 2.407137e-01          NaN          NaN 4.213567e-01 1.779583e-02
    ## [11186] 5.812799e-02 1.055566e-07 7.837494e-01 9.890134e-01          NaN
    ## [11191] 1.412526e-02 3.825403e-01 5.426662e-06 4.520680e-01 4.763541e-03
    ## [11196] 6.889658e-09          NaN 9.388336e-03 6.020237e-01 2.011118e-01
    ## [11201] 2.867382e-01 8.987349e-01 5.726908e-11 7.594462e-01 9.484867e-01
    ## [11206] 1.743804e-02 1.083791e-01 2.333738e-01 1.794297e-01 1.856350e-01
    ## [11211] 8.108340e-01 8.911689e-01 2.557325e-01 9.503321e-01 4.901700e-01
    ## [11216] 1.268184e-04 1.425784e-04 9.455707e-01 3.217856e-01 4.370637e-04
    ## [11221] 1.866828e-03 2.013074e-12 1.475187e-01 2.626624e-01 6.836717e-04
    ## [11226] 9.732872e-02 7.584522e-01 3.562914e-01 8.377140e-03 1.197789e-01
    ## [11231] 9.283536e-09 6.494807e-01 2.923724e-03 2.410573e-02 5.819166e-01
    ## [11236]          NaN 6.276786e-03 1.643423e-02 8.875787e-01 2.824966e-03
    ## [11241] 5.088388e-01 2.637014e-02 1.636222e-01 9.787885e-01 8.032544e-01
    ## [11246]          NaN 2.876673e-04 1.586265e-04 5.562053e-06 3.830034e-02
    ## [11251] 4.875163e-13 7.200178e-02 1.421608e-01          NaN 5.512755e-02
    ## [11256] 1.566039e-01 1.334770e-02 8.060992e-01 1.011324e-05 9.591596e-01
    ## [11261] 3.117636e-01 1.755616e-02 2.066555e-06          NaN 5.509863e-02
    ## [11266] 7.425111e-01 7.857125e-01 5.878093e-01 3.204838e-01          NaN
    ## [11271] 1.148377e-01 2.242891e-11 3.992438e-01 4.103662e-03 2.817649e-08
    ## [11276] 1.950807e-01 1.242019e-01 9.493293e-01 6.701847e-01 9.569496e-01
    ## [11281]          NaN 9.028124e-06 2.930285e-02 4.078207e-01 8.207522e-05
    ## [11286] 7.019552e-04 3.630738e-01 4.559157e-01 4.650441e-01 5.132131e-01
    ## [11291] 5.869117e-01 2.845619e-03 3.421916e-41 8.656771e-02 5.215308e-05
    ## [11296] 1.741592e-14          NaN 1.286736e-02 1.997590e-01 7.954210e-02
    ## [11301] 2.474904e-01 1.623244e-01 5.692111e-01 3.204838e-01 8.840560e-02
    ## [11306] 4.747710e-01 3.618744e-03 5.686457e-01 5.653889e-01 7.682588e-01
    ## [11311] 1.565238e-04 1.709384e-01 1.219893e-01 1.120419e-06 2.232056e-03
    ## [11316] 1.401626e-10          NaN 4.393834e-01 4.311168e-01          NaN
    ## [11321] 8.303539e-01 1.706918e-15 6.314813e-10 3.786777e-01 7.122051e-02
    ## [11326] 2.992112e-01 9.471268e-02 2.566986e-01 2.858925e-01 1.366017e-02
    ## [11331] 2.125923e-01 1.327247e-02 3.204838e-01 3.206358e-09 2.277158e-02
    ## [11336]          NaN 5.635676e-01 7.529882e-02 5.968702e-07 2.435374e-06
    ## [11341] 9.177213e-04 9.588949e-01 2.308400e-04          NaN 8.025358e-01
    ## [11346] 8.095264e-02 2.500641e-01 1.354916e-02 3.693674e-01 4.656384e-01
    ## [11351] 8.542370e-02 5.282606e-01 1.011315e-01          NaN 6.117381e-01
    ## [11356] 1.365706e-03 5.409561e-03 8.119020e-02 4.773632e-02 7.341441e-03
    ## [11361] 1.506191e-01 9.955944e-01 1.518286e-01 3.056283e-07 3.958168e-01
    ## [11366] 1.557395e-01 1.488540e-14 7.985939e-02 4.232725e-03 1.093708e-09
    ## [11371] 3.109650e-04 2.124862e-01 4.923977e-01 5.922710e-01 1.113647e-01
    ## [11376] 4.469425e-04 1.368459e-03 6.757787e-01 5.183295e-02 3.679615e-02
    ## [11381] 1.788770e-01 2.050564e-14 7.603440e-01 7.486550e-01 2.313573e-01
    ## [11386] 4.399051e-01 2.987729e-01 7.881262e-01 7.451749e-05 4.065865e-01
    ## [11391] 3.251188e-01 6.674303e-03 1.395613e-02 4.645105e-01          NaN
    ## [11396] 7.203817e-01 8.574897e-01 3.128878e-01 8.724782e-01 1.448090e-01
    ## [11401] 3.191779e-05 3.530680e-02 7.923146e-03 9.589443e-01 2.188890e-10
    ## [11406] 2.560631e-12 1.872851e-01 3.204838e-01          NaN 5.012424e-01
    ## [11411] 3.830973e-03 2.710750e-04 6.383470e-03 7.560599e-02 7.709525e-02
    ## [11416] 8.615942e-01 8.539277e-01 6.055765e-01 8.528756e-01 1.205969e-01
    ## [11421] 8.949341e-01 4.382384e-03 1.898274e-01 5.224713e-03 1.167850e-01
    ## [11426] 4.128731e-02 8.811985e-01 6.626052e-02 4.999852e-01          NaN
    ## [11431] 9.743102e-02 6.900727e-01 9.893448e-01 2.775125e-02 7.058325e-03
    ## [11436] 2.039433e-08 9.142704e-02 1.438201e-03 1.994234e-25 2.361538e-04
    ## [11441] 4.590099e-01 6.778871e-02          NaN 1.155995e-05 3.509311e-01
    ## [11446] 1.551562e-02 3.611981e-01 2.034367e-11 3.504761e-01 8.700783e-04
    ## [11451] 8.693542e-01 9.352784e-01 8.106159e-01 8.258181e-02 1.227098e-04
    ## [11456] 3.049973e-01 6.828800e-01          NaN 2.695948e-01 1.943900e-01
    ## [11461] 4.737716e-03 5.578265e-02 7.712269e-02 1.729420e-04 2.132952e-36
    ## [11466] 6.291012e-01 4.255156e-02 4.260718e-01 8.076032e-04 6.753238e-02
    ## [11471] 6.797235e-01 8.745406e-01 5.636016e-01 3.185298e-03 7.064497e-02
    ## [11476] 5.909931e-04 4.605993e-01          NaN 2.751679e-02 2.660920e-01
    ## [11481] 1.094761e-01          NaN 9.835457e-01 7.808810e-02 3.815925e-19
    ## [11486] 6.926611e-01 6.028817e-22 4.305435e-01 3.464067e-02 9.253819e-05
    ## [11491] 2.495500e-01 2.620079e-01 3.074478e-01          NaN 6.038743e-01
    ## [11496] 4.315123e-15 7.192110e-01 8.711416e-02 2.466673e-01 8.402055e-09
    ## [11501] 3.226108e-17 2.550339e-02 2.893494e-05 2.584062e-03 7.962838e-01
    ## [11506] 8.534203e-01          NaN 5.292255e-01 2.588412e-02          NaN
    ## [11511] 2.333899e-01 3.921866e-03 1.380112e-06 3.597849e-03 5.488699e-01
    ## [11516] 9.535789e-01 1.397078e-06 2.025261e-01 1.280913e-03 1.982416e-01
    ## [11521] 8.147641e-03 3.204838e-01          NaN 6.696855e-01          NaN
    ## [11526] 2.593490e-01 7.454097e-01 4.633715e-03 6.732281e-02 5.879055e-03
    ## [11531] 2.690369e-01 2.601511e-09          NaN 5.290480e-01 9.130032e-01
    ## [11536] 1.917894e-10 8.135520e-01 8.897913e-03 4.247179e-01 3.369665e-03
    ## [11541] 6.788787e-01 2.594776e-10 5.948889e-03 4.061376e-01          NaN
    ## [11546] 1.401352e-05 1.315610e-03 2.352340e-01 4.895363e-01 3.821736e-01
    ## [11551] 4.738523e-07 1.349800e-01 9.239826e-01 1.366841e-04 1.529468e-02
    ## [11556] 8.104737e-02          NaN 2.911447e-04 4.014003e-01 1.559026e-01
    ## [11561] 9.710082e-01 4.451426e-04 1.899730e-60 5.432627e-01 8.254445e-01
    ## [11566] 1.822685e-18 1.572434e-01 4.190352e-15 6.351103e-01 2.310774e-02
    ## [11571] 7.245885e-05 3.905196e-11 9.617957e-02 1.904083e-01 1.764976e-01
    ## [11576] 8.161957e-01          NaN 3.458239e-02          NaN 1.157680e-02
    ## [11581] 8.463242e-01 5.412166e-01 1.351000e-02 7.676998e-01 4.642621e-01
    ## [11586] 7.946028e-04 5.789958e-09 4.283271e-01 2.582364e-01 6.712414e-01
    ## [11591] 5.340112e-05 2.308787e-04 1.599123e-01 3.655982e-01 8.564116e-01
    ## [11596] 4.220293e-01 9.842426e-04 1.024594e-03 3.204838e-01 1.755054e-01
    ## [11601] 9.508954e-01 1.074716e-01 1.933605e-02 4.791878e-03 1.944251e-04
    ## [11606] 1.000973e-22 9.380875e-01 3.204077e-02          NaN 1.006573e-01
    ## [11611] 6.744223e-01 2.596284e-01 3.463892e-03          NaN          NaN
    ## [11616] 5.032331e-01 1.862224e-02 7.204186e-01 2.741126e-08 1.216231e-02
    ## [11621] 6.992309e-02 2.331436e-01 7.507361e-01 2.378011e-01 1.192247e-02
    ## [11626] 1.909457e-02 1.948218e-04 1.293883e-02 1.519477e-02 3.119468e-02
    ## [11631] 1.111147e-01 1.057812e-02          NaN 1.655502e-16 4.989612e-02
    ## [11636] 1.667753e-01 2.685230e-02 6.666173e-04 8.529339e-04 7.394021e-01
    ## [11641] 5.051787e-01 6.271950e-02 1.264886e-01 2.384330e-02 1.975640e-01
    ## [11646] 4.820793e-02 9.003983e-03          NaN          NaN          NaN
    ## [11651] 5.325401e-01 2.171428e-01 2.793768e-09 7.483367e-01 5.623563e-02
    ## [11656] 1.131694e-03 4.087325e-01 1.713284e-13 1.361478e-13 4.274642e-02
    ## [11661] 7.689654e-02 8.536507e-01 8.577791e-09 7.595247e-01 9.623710e-01
    ## [11666] 1.726475e-08 1.821239e-02 9.745174e-04 3.174218e-28 7.677272e-01
    ## [11671] 1.499926e-01 1.794130e-02 5.330870e-02 6.146237e-07 7.528381e-01
    ## [11676] 7.704148e-07 2.539210e-03 4.164495e-01          NaN 3.884910e-03
    ## [11681] 1.892772e-01 4.749626e-03 6.081542e-06 1.639841e-01 2.963251e-02
    ## [11686] 4.654595e-01 7.996729e-02 3.728814e-04          NaN 1.140680e-17
    ## [11691]          NaN 3.375364e-01 3.337194e-01 6.966035e-01 1.635392e-05
    ## [11696] 3.745684e-15 7.169192e-01 6.389017e-04          NaN 3.072633e-03
    ## [11701] 1.148504e-01 3.495843e-09          NaN 3.204838e-01 9.312748e-16
    ## [11706] 1.612178e-02 4.413686e-01 9.801380e-01 3.418515e-01 6.273270e-04
    ## [11711] 9.760993e-01 2.581558e-01 4.634941e-01 3.106207e-02 7.496096e-01
    ## [11716] 7.087705e-01 3.204838e-01 4.000022e-01 4.017254e-05 1.623938e-02
    ## [11721] 5.591627e-01          NaN 5.519167e-01 1.846831e-13          NaN
    ## [11726] 4.754246e-01 8.075064e-01 5.691315e-02 1.596504e-10 4.758189e-07
    ## [11731] 4.438633e-02 7.769461e-02 8.206337e-06 6.135608e-02 7.553639e-01
    ## [11736] 2.525422e-01 3.978686e-06          NaN          NaN 7.636536e-03
    ## [11741] 2.161720e-04 1.364990e-01 1.729312e-03 6.258228e-04 2.024129e-25
    ## [11746] 2.601245e-02 4.709061e-19 1.872015e-01 1.331947e-01 4.121840e-01
    ## [11751] 9.575935e-07 9.385996e-01 3.710228e-01 3.928600e-01 6.122562e-01
    ## [11756] 1.339602e-01 2.178206e-01 4.601329e-03 5.792494e-05 3.329248e-04
    ## [11761] 1.013934e-01 1.958757e-01 2.064021e-06          NaN 3.204838e-01
    ## [11766]          NaN 2.566210e-06 3.204838e-01 1.194075e-01 2.944971e-01
    ## [11771] 7.775889e-04 4.294887e-03 2.406306e-03          NaN 3.935050e-03
    ## [11776] 3.204838e-01 3.204838e-01 2.635318e-03 4.278744e-01 2.159878e-02
    ## [11781] 4.514130e-01 9.225911e-01 1.840503e-02 1.535238e-03 4.292746e-01
    ## [11786]          NaN 2.291990e-01 7.841719e-01 6.674659e-01          NaN
    ## [11791] 3.892458e-03 1.559125e-01 3.435704e-01 3.448412e-01          NaN
    ## [11796]          NaN 1.896853e-01 8.473356e-01 6.450077e-01          NaN
    ## [11801]          NaN 2.433462e-11 3.204838e-01 1.961800e-08 1.200816e-03
    ## [11806] 7.263547e-02 2.950633e-02 3.838856e-01 2.849612e-01 2.801546e-24
    ## [11811] 5.769175e-01 8.796880e-11          NaN 5.782399e-02          NaN
    ## [11816] 8.780502e-01 1.746359e-03          NaN 3.876270e-01 4.224312e-01
    ## [11821]          NaN 1.916651e-03 1.584917e-01 1.478786e-02 3.204838e-01
    ## [11826] 2.696635e-05 4.968759e-02 6.025166e-01 4.393917e-01 4.235012e-01
    ## [11831] 3.336016e-04 4.892706e-01 1.247277e-08 9.221286e-02 4.832405e-01
    ## [11836] 8.724937e-04 6.609129e-01 4.677072e-09 8.680092e-01          NaN
    ## [11841]          NaN 5.937057e-02 7.440878e-01 8.983354e-01 6.257021e-09
    ## [11846] 3.175605e-01 8.777349e-01 6.740630e-02 8.331439e-02          NaN
    ## [11851] 5.549388e-02 2.835715e-01 6.865106e-01 6.921166e-01 1.871641e-02
    ## [11856] 1.299934e-01 5.869120e-02 4.941167e-01 2.197216e-01 8.762557e-01
    ## [11861]          NaN 2.032743e-01 8.484030e-03 4.182751e-03 1.899057e-02
    ## [11866] 6.230140e-05          NaN 1.508712e-03 7.235813e-03 2.105291e-03
    ## [11871] 5.204010e-02 1.419231e-01 3.462184e-08 2.934402e-02          NaN
    ## [11876] 3.185218e-01 2.296785e-03 4.229577e-01 3.858833e-01 3.149580e-02
    ## [11881] 6.939684e-04 9.138354e-14 2.959869e-01 3.342331e-10 4.293060e-04
    ## [11886] 9.001802e-01 1.214394e-01 1.478789e-02 2.797421e-03 3.893110e-23
    ## [11891] 2.367839e-02          NaN          NaN 4.257748e-02 7.159659e-03
    ## [11896]          NaN 6.437409e-01 2.775908e-01 6.023553e-04 1.550155e-02
    ## [11901] 2.902585e-01 8.369069e-01 6.522366e-01 7.243657e-01 1.766189e-01
    ## [11906] 2.274090e-01 1.123899e-22 1.098653e-02 6.767387e-02 6.136328e-02
    ## [11911]          NaN 7.082070e-01 1.958365e-03 2.822176e-01 7.518636e-01
    ## [11916] 4.686230e-01 2.014992e-02 3.185006e-01 6.470855e-03 3.183204e-02
    ## [11921] 1.887968e-04 1.697374e-07 2.110677e-03 1.837047e-01          NaN
    ## [11926] 8.674427e-02 2.967628e-01          NaN 6.120957e-01 7.464820e-16
    ## [11931] 5.394623e-11 3.020547e-02 9.104729e-02 6.988725e-01 2.760503e-05
    ## [11936] 8.490063e-01 5.817082e-01 8.492150e-01 3.799475e-03 2.065116e-01
    ## [11941] 9.787881e-03 7.480454e-01 5.392324e-06 5.749837e-01 4.588104e-04
    ## [11946] 4.222851e-01 2.120445e-02 1.830607e-01 9.515496e-17 5.261345e-01
    ## [11951] 2.055373e-02          NaN          NaN 2.129411e-01 1.537834e-01
    ## [11956] 8.658498e-01 1.358837e-02 1.886358e-12 7.959104e-01 9.879467e-09
    ## [11961] 1.312749e-03 8.544796e-01 3.007883e-01 1.164836e-02 2.067615e-01
    ## [11966] 1.390461e-09 3.204838e-01          NaN 5.694853e-01 4.496432e-01
    ## [11971] 3.204838e-01          NaN 5.479881e-02 5.739837e-05 6.085132e-08
    ## [11976] 1.979807e-05          NaN 1.555138e-07 2.606159e-06 1.688730e-01
    ## [11981] 8.643256e-02 3.575175e-01 1.995866e-01 4.947053e-23 8.277924e-09
    ## [11986] 9.640524e-02 1.006515e-01 4.025918e-05 1.823855e-10 1.045741e-01
    ## [11991] 9.446855e-03 3.093019e-01 4.308167e-08          NaN 3.327311e-02
    ## [11996] 3.023800e-02 8.594322e-02 2.210500e-10 2.113257e-01 3.335061e-01
    ## [12001] 3.198643e-01 8.206916e-01 3.369420e-01 1.197678e-02 2.308548e-01
    ## [12006] 1.175448e-01 5.404138e-05 2.632157e-01 3.204838e-01 5.423865e-01
    ## [12011] 1.885535e-02 7.123335e-01 9.910668e-01          NaN 4.619159e-01
    ## [12016] 2.414555e-02 9.682593e-01 8.555132e-02 2.795989e-05 6.223413e-04
    ## [12021] 7.560829e-03 7.836685e-01          NaN          NaN 5.216688e-01
    ## [12026] 3.043230e-07 3.204838e-01          NaN 5.140788e-11 5.774563e-01
    ## [12031] 1.628517e-08 4.852866e-02 3.233147e-01 3.393706e-31 3.877401e-01
    ## [12036] 5.304458e-03 9.742407e-02 9.791970e-01 4.822863e-01 3.644051e-01
    ## [12041] 3.560381e-01 1.569999e-14 4.387169e-04 7.127840e-01 3.461507e-01
    ## [12046] 8.254379e-01 2.156034e-11 2.063652e-04 6.307434e-01 5.680421e-01
    ## [12051] 1.335015e-02 9.706972e-01 3.542574e-01 3.399410e-01 5.084110e-03
    ## [12056]          NaN 8.208999e-02 7.390411e-01          NaN 1.364239e-01
    ## [12061] 4.287573e-01 8.227655e-01 1.890496e-01 5.337801e-01 6.126898e-02
    ## [12066] 4.826119e-02 3.204838e-01 3.204838e-01 7.992933e-01 5.036701e-01
    ## [12071] 1.680500e-01 1.122917e-03 9.288967e-01 4.042876e-04 6.928632e-02
    ## [12076] 8.840895e-01          NaN 1.787656e-05 6.213899e-02 4.334759e-01
    ## [12081] 7.038990e-05 4.200380e-01          NaN 3.204838e-01 7.732792e-01
    ## [12086] 9.629983e-01          NaN 3.824964e-01 2.906364e-01 2.263554e-03
    ## [12091] 1.631687e-01 9.904187e-01          NaN 6.110093e-01 4.151021e-01
    ## [12096] 6.055848e-01 3.318476e-01 9.831032e-02 6.747687e-01 5.259290e-03
    ## [12101] 5.700184e-02 4.851573e-01 3.204838e-01 3.478338e-24 2.353154e-01
    ## [12106] 3.181732e-01 8.962010e-05 3.560159e-03 7.434361e-01 1.114467e-01
    ## [12111] 2.971867e-04 4.956926e-01 7.902597e-05 6.207089e-01 9.254639e-01
    ## [12116]          NaN 1.582853e-01 1.649178e-01 3.204838e-01          NaN
    ## [12121] 3.328427e-07 1.835993e-01 3.792474e-01 7.188943e-07 1.671777e-11
    ## [12126] 8.140266e-01 8.128964e-02 5.367149e-22 7.230953e-06          NaN
    ## [12131] 6.094321e-06          NaN 1.876787e-01 1.014771e-01 4.022979e-08
    ## [12136] 4.189989e-03 9.412578e-01 2.322737e-02 2.811499e-43 3.186561e-02
    ## [12141] 9.249792e-01 9.319619e-01 4.463493e-01 7.434011e-02 1.148854e-01
    ## [12146] 7.771285e-03 5.024445e-02 1.177490e-05 7.996630e-01 2.736547e-02
    ## [12151]          NaN 5.098979e-01 2.340104e-02 1.087155e-15 1.529880e-02
    ## [12156] 9.260120e-01 1.970062e-01 5.655949e-01 4.729404e-03 2.661615e-05
    ## [12161]          NaN          NaN          NaN 7.991061e-01 1.694568e-01
    ## [12166]          NaN 1.172390e-02 3.061276e-01 1.913368e-09 8.917984e-01
    ## [12171]          NaN 1.061074e-03 3.678988e-01 2.410933e-02 6.365132e-08
    ## [12176] 8.740220e-02 2.671512e-01 4.376968e-02 1.133679e-08 9.758168e-03
    ## [12181] 3.344402e-01 8.410289e-01 8.639040e-02 6.606928e-03 9.991372e-01
    ## [12186] 1.652080e-06 5.143042e-01 3.204838e-01 1.299330e-01 7.837526e-01
    ## [12191] 7.546831e-02 1.529659e-01 5.053413e-02 3.204838e-01 3.616979e-01
    ## [12196] 1.624044e-01 3.854321e-01 5.217602e-02 6.462026e-02 8.255966e-08
    ## [12201] 1.173347e-01 6.174173e-01 9.714781e-01 2.518513e-02 5.999574e-07
    ## [12206] 1.257604e-02 4.897663e-01 8.068700e-06 9.135414e-01 2.966201e-01
    ## [12211] 2.906877e-01 3.156102e-03 3.290039e-01 6.974672e-02 6.276524e-01
    ## [12216] 9.179794e-01 8.000813e-01 8.104073e-06 3.204838e-01 1.398176e-07
    ## [12221] 2.210985e-04 9.532972e-01 1.175316e-03 3.709509e-01 2.976966e-01
    ## [12226]          NaN 1.119810e-01 1.568602e-01 3.339384e-01 4.685201e-01
    ## [12231] 6.199071e-01 3.280237e-03 3.204838e-01 1.934087e-01 3.270648e-01
    ## [12236] 2.232869e-01 5.325936e-11 7.025100e-05 1.166363e-01 5.173900e-04
    ## [12241]          NaN 4.621096e-01 3.327659e-02 1.684641e-01 4.617287e-10
    ## [12246] 6.652815e-01 4.512839e-18 1.992235e-02 8.545329e-01          NaN
    ## [12251] 3.204838e-01 9.141932e-01 7.188971e-02 8.234166e-01 1.206717e-01
    ## [12256]          NaN 6.314362e-03 1.881149e-09          NaN 5.677986e-01
    ## [12261] 3.150183e-06 4.229610e-02 1.092387e-20 2.482972e-02 9.288070e-01
    ## [12266] 3.520504e-05 2.270742e-03 1.403771e-09 1.132639e-01          NaN
    ## [12271] 1.078330e-04 2.280305e-01 3.561427e-01 9.688378e-01 3.455676e-02
    ## [12276] 1.576274e-03 7.732874e-01 2.478977e-01 2.257305e-25 7.289986e-18
    ## [12281] 9.725850e-01 8.917379e-01 9.759514e-06 1.986297e-01 2.152061e-03
    ## [12286] 5.314789e-05 2.247545e-02 8.805020e-02 5.011424e-11 4.891972e-02
    ## [12291] 5.925434e-03 1.448930e-27 2.412257e-20 2.316840e-02 4.615404e-02
    ## [12296] 3.037790e-05 1.021994e-01          NaN 4.960826e-01 4.407251e-03
    ## [12301] 4.509529e-11 2.350776e-06 9.256606e-05 1.074901e-04 1.821709e-06
    ## [12306] 7.364269e-04          NaN 3.715331e-01 1.309150e-02 9.239978e-02
    ## [12311] 2.919166e-02 4.465104e-01 7.794578e-01          NaN 1.730662e-01
    ## [12316] 1.551112e-03 3.394840e-02          NaN 9.685949e-01 2.253987e-07
    ## [12321] 9.228090e-01 7.429244e-01 6.679509e-03 3.783946e-03 7.407436e-02
    ## [12326] 2.608927e-02 1.691136e-04          NaN 6.914669e-01 9.427844e-01
    ## [12331] 1.829611e-01 1.073213e-04 7.046113e-05 2.055539e-04 9.994412e-01
    ## [12336] 8.235600e-01          NaN 2.805053e-03 2.339561e-02          NaN
    ## [12341] 3.060576e-02          NaN          NaN 7.897569e-04 4.540203e-01
    ## [12346] 9.703975e-15 1.728203e-01 9.494806e-03 9.664436e-02 3.474366e-01
    ## [12351]          NaN 1.490742e-10 4.039198e-01          NaN 8.227610e-02
    ## [12356]          NaN 2.214492e-02 7.148650e-05 1.332245e-01 1.402997e-04
    ## [12361] 9.880139e-01 8.343571e-04 2.249663e-01          NaN 3.128821e-01
    ## [12366] 3.993110e-01 9.818869e-01 1.047169e-02 9.331656e-01          NaN
    ## [12371] 5.229436e-01 5.498192e-01          NaN          NaN          NaN
    ## [12376] 4.112016e-08 7.455140e-01 7.686535e-06          NaN          NaN
    ## [12381] 3.364806e-02 2.500596e-01          NaN 1.098646e-02 5.206670e-03
    ## [12386] 6.388355e-11 6.732320e-01 2.197039e-23 1.415422e-08          NaN
    ## [12391] 4.658480e-01 6.233139e-01 4.934719e-01 4.086745e-01 6.510131e-07
    ## [12396]          NaN 5.443927e-04 8.384948e-11 3.773724e-03 3.346831e-01
    ## [12401] 4.334337e-05 3.204838e-01          NaN 2.869037e-02          NaN
    ## [12406] 5.891429e-02 2.518665e-01 5.734175e-08          NaN          NaN
    ## [12411] 1.388105e-01 7.746280e-01 8.098504e-01 5.811127e-04 8.200873e-03
    ## [12416]          NaN 1.273135e-01 7.961489e-02 1.333203e-03 4.435287e-04
    ## [12421] 3.204838e-01          NaN 7.918846e-01 7.935890e-01          NaN
    ## [12426] 2.538298e-03 3.102641e-05          NaN 3.154517e-02 5.757751e-12
    ## [12431] 2.934237e-02 2.458117e-03 2.756392e-02 2.514335e-02 5.923310e-15
    ## [12436] 2.321974e-02 5.361009e-01 6.240540e-01 7.777573e-33 2.690853e-03
    ## [12441] 4.636473e-01 5.589652e-01 1.107884e-01          NaN 9.493427e-02
    ## [12446] 8.666584e-02 1.994518e-02 8.626079e-07 3.630213e-01 3.204838e-01
    ## [12451] 9.408366e-02 9.541866e-01          NaN          NaN          NaN
    ## [12456]          NaN 2.502367e-01 8.992693e-02 3.614764e-10 7.805868e-02
    ## [12461] 5.049141e-20 4.181004e-07 8.953744e-04 5.763246e-01 5.597365e-01
    ## [12466] 3.216662e-01 4.874143e-07          NaN 3.764475e-03 5.344261e-01
    ## [12471] 1.007753e-01 1.874136e-02 1.582300e-01          NaN 1.357949e-01
    ## [12476] 7.707268e-03 7.268319e-01          NaN 1.565069e-01 1.620701e-01
    ## [12481] 7.924360e-01 1.473438e-13 9.524963e-01 3.503246e-02 1.899477e-05
    ## [12486]          NaN 2.822234e-01 4.121543e-01 8.596221e-02 6.892284e-02
    ## [12491] 1.095574e-01          NaN          NaN 7.931132e-03          NaN
    ## [12496]          NaN 1.779727e-07 1.590993e-01 1.318773e-06 1.859213e-12
    ## [12501] 1.970709e-02 9.187372e-01 2.172873e-01 2.992186e-02 3.219913e-01
    ## [12506]          NaN 2.156704e-01 9.283760e-10 2.410687e-09          NaN
    ## [12511] 1.936799e-01 5.238193e-01 2.121433e-20 7.302969e-02 3.593393e-17
    ## [12516]          NaN 1.103927e-24 6.978322e-01 7.011945e-01 5.742381e-06
    ## [12521] 5.371469e-02 9.978388e-01          NaN 3.204838e-01          NaN
    ## [12526] 8.192259e-01          NaN 2.291850e-09 9.370568e-01 1.423491e-07
    ## [12531]          NaN 3.326166e-01          NaN 7.928019e-29          NaN
    ## [12536] 2.806230e-01 9.940166e-01 4.632566e-01 3.565295e-04 3.018088e-01
    ## [12541] 6.713026e-06          NaN          NaN 9.713287e-02 3.173986e-02
    ## [12546] 2.273938e-03 2.752361e-02 3.204838e-01 1.172238e-04 1.182479e-26
    ## [12551]          NaN          NaN 1.242129e-11          NaN 4.527489e-01
    ## [12556]          NaN 1.457809e-08 3.929633e-24 1.239058e-01 3.108573e-01
    ## [12561] 1.585819e-01          NaN 2.904668e-01 3.204838e-01          NaN
    ## [12566] 2.646262e-02 5.485277e-01 3.847627e-18 8.749863e-01 4.871190e-01
    ## [12571] 4.480829e-01 8.909185e-01 4.176119e-01          NaN 9.312086e-01
    ## [12576] 2.373718e-09 2.081713e-01 1.051268e-04 1.602011e-03 2.619124e-01
    ## [12581] 8.050172e-01 5.275036e-02 4.388678e-12          NaN          NaN
    ## [12586] 1.911309e-03 2.929703e-03 1.034535e-21 5.220091e-01 2.808072e-01
    ## [12591] 9.553145e-01 4.077057e-02 1.678433e-03 8.362986e-02          NaN
    ## [12596] 8.694797e-01 9.491097e-02 1.926631e-01 7.836990e-01          NaN
    ## [12601] 2.079994e-03 5.237216e-01 1.726786e-01 3.862809e-04          NaN
    ## [12606] 9.379924e-01 2.272672e-03 1.480229e-14 3.781579e-01 1.912471e-01
    ## [12611] 2.215695e-03 2.549034e-07 5.719372e-01 2.775219e-10 5.154684e-04
    ## [12616] 6.830440e-03 4.471502e-06 4.852039e-10          NaN          NaN
    ## [12621] 4.346361e-01 3.937994e-03 6.489796e-01 1.029438e-05 1.627832e-01
    ## [12626] 7.080411e-05 1.290318e-03          NaN 6.188938e-01 7.177762e-11
    ## [12631] 6.331731e-02 2.429346e-08 2.698263e-01 1.461060e-01 1.157420e-02
    ## [12636] 4.376610e-01 9.940983e-01 1.841588e-02 9.205895e-02 8.778823e-01
    ## [12641] 3.204838e-01 7.434890e-01 9.774603e-04 5.210947e-05 1.028070e-02
    ## [12646] 9.747314e-04 1.284458e-01          NaN          NaN 1.803993e-01
    ## [12651] 7.535144e-01          NaN 5.380441e-53 4.583246e-01          NaN
    ## [12656] 3.374506e-09 3.204838e-01 7.945462e-01          NaN          NaN
    ## [12661] 8.085720e-08 2.233471e-01 1.652570e-01 7.646490e-01 3.204838e-01
    ## [12666] 1.998415e-18 1.102700e-17          NaN 1.447576e-03 3.204838e-01
    ## [12671] 3.839918e-05 8.146418e-01          NaN 6.610882e-01 1.833801e-01
    ## [12676] 4.188703e-01 2.279113e-01 5.835316e-06 3.204838e-01 8.166900e-02
    ## [12681] 5.496051e-01 5.967487e-01 1.602913e-01 9.252893e-02 2.733120e-01
    ## [12686]          NaN          NaN 1.302613e-02 5.708000e-01 1.624260e-03
    ## [12691] 6.082693e-01          NaN 2.341330e-01          NaN 7.232566e-01
    ## [12696] 1.659288e-02 1.705572e-05 3.404287e-02 1.519953e-01 1.568901e-01
    ## [12701] 2.343354e-02 6.408317e-02 1.030469e-53 4.232941e-01 4.245075e-01
    ## [12706] 8.418715e-02          NaN 4.150907e-02 7.995271e-01 1.584141e-01
    ## [12711] 4.688162e-10          NaN          NaN 5.426460e-01 1.728431e-05
    ## [12716] 5.844026e-01 3.064848e-03 2.848531e-07 2.202200e-01 9.104825e-01
    ## [12721] 6.053046e-01 1.179686e-16 5.473557e-09 3.031438e-01 1.922830e-01
    ## [12726] 4.357662e-12 5.562577e-05          NaN 5.698554e-01 3.204838e-01
    ## [12731] 5.829202e-20          NaN          NaN 1.078070e-02 8.740564e-02
    ## [12736] 2.052866e-02 6.193880e-01 1.239081e-05          NaN          NaN
    ## [12741] 1.485186e-02 1.106426e-01 3.340232e-05 2.424489e-01 9.179006e-16
    ## [12746] 8.779015e-09 1.335205e-31 1.694116e-01 2.949321e-02 8.205554e-02
    ## [12751] 4.566205e-01 5.030447e-01 3.204838e-01 2.913709e-01 3.236408e-01
    ## [12756] 1.184228e-05          NaN 6.910138e-01 2.614243e-05          NaN
    ## [12761] 1.938088e-02 6.599614e-01 6.742930e-01 1.106500e-02          NaN
    ## [12766] 4.878763e-05          NaN          NaN 1.483206e-03 2.023943e-14
    ## [12771] 3.204838e-01 1.143627e-01 5.087934e-03 3.204838e-01 3.208757e-04
    ## [12776] 6.453778e-01 6.710908e-02 1.938525e-01          NaN 4.028267e-02
    ## [12781] 1.601976e-01          NaN 4.302575e-03          NaN          NaN
    ## [12786] 8.139115e-01 8.724821e-04 3.761155e-01          NaN 5.345816e-02
    ## [12791] 6.407442e-04 1.746107e-01 6.735594e-01 4.748236e-14 3.535216e-01
    ## [12796] 8.282249e-01 1.196275e-01 6.495780e-01 1.872960e-02 3.530915e-07
    ## [12801] 1.110044e-02          NaN 3.986679e-05 2.327969e-11 5.174542e-19
    ## [12806] 7.235591e-01 4.529484e-02 1.416953e-04 5.901954e-01 9.343726e-01
    ## [12811]          NaN          NaN 9.242945e-01          NaN 7.843527e-15
    ## [12816] 9.560084e-02          NaN 6.665134e-05          NaN          NaN
    ## [12821] 3.847409e-01 7.029077e-01 4.769126e-06 4.170390e-03 1.632156e-01
    ## [12826] 7.644510e-01 6.842846e-01 3.204838e-01 8.952409e-01          NaN
    ## [12831]          NaN 3.794549e-02 8.707096e-01          NaN 9.806625e-01
    ## [12836] 5.021548e-03          NaN          NaN          NaN 9.257282e-03
    ## [12841]          NaN 9.656846e-01 3.203910e-64 3.204838e-01          NaN
    ## [12846] 2.551761e-03 3.204838e-01          NaN 1.048854e-03          NaN
    ## [12851] 6.134557e-03 3.298203e-01          NaN 9.805227e-02 4.761957e-02
    ## [12856] 2.748056e-01 5.496446e-04 1.926994e-01 2.405073e-01 7.332581e-01
    ## [12861] 1.016280e-02 4.030450e-02 1.226694e-10          NaN 2.493793e-03
    ## [12866] 3.204838e-01 1.136478e-01 5.165436e-05          NaN 3.204838e-01
    ## [12871] 4.445917e-01 7.678797e-01 6.948134e-01          NaN          NaN
    ## [12876] 3.023880e-01 5.638949e-03 7.085058e-01 3.204838e-01 8.261069e-01
    ## [12881] 2.018410e-03 3.351589e-14 4.403466e-01 9.738559e-01 2.580580e-03
    ## [12886] 3.409032e-05 3.744756e-01 3.643499e-03          NaN 3.657792e-02
    ## [12891] 9.839906e-01          NaN          NaN          NaN          NaN
    ## [12896] 2.652631e-01 7.823284e-07 3.472097e-01 3.204838e-01 2.769804e-02
    ## [12901] 5.324663e-03 5.858856e-01          NaN          NaN          NaN
    ## [12906]          NaN 9.551743e-01          NaN          NaN 3.447079e-01
    ## [12911] 3.004377e-01 1.197616e-08 9.553644e-01 1.026833e-01 6.037180e-02
    ## [12916] 4.458499e-06 8.857326e-01          NaN          NaN 8.179940e-01
    ## [12921] 1.157956e-01 5.669846e-02 1.645734e-01 7.333986e-04 5.299548e-01
    ## [12926] 1.558098e-01          NaN          NaN          NaN 9.819281e-01
    ## [12931]          NaN 3.344918e-01 6.011777e-01 8.760135e-01 1.144280e-01
    ## [12936] 4.774591e-12 2.038094e-26 6.594889e-01 6.762520e-02 3.204838e-01
    ## [12941] 3.576471e-32 6.861049e-01 1.156900e-15          NaN 4.917101e-03
    ## [12946] 1.607531e-02 3.494684e-08 9.331230e-16 1.120540e-02 2.689235e-01
    ## [12951] 4.512375e-01 3.204838e-01          NaN 3.740289e-01 2.412049e-02
    ## [12956] 6.489062e-02          NaN 6.971732e-01 4.489404e-02 6.846579e-01
    ## [12961] 2.986600e-08          NaN          NaN 9.500595e-01 7.898843e-02
    ## [12966] 2.091040e-02 9.574270e-03 2.865972e-06          NaN 2.221601e-01
    ## [12971] 3.042455e-06 3.204838e-01 4.993360e-03 7.522650e-02 9.154819e-01
    ## [12976] 7.466782e-01 5.575685e-07          NaN 3.960496e-01          NaN
    ## [12981] 3.204838e-01 1.060935e-01 6.088914e-01          NaN 6.076301e-01
    ## [12986] 1.759388e-02 3.204838e-01 6.218008e-05 3.216469e-01 1.499696e-01
    ## [12991] 6.327176e-02          NaN          NaN 9.077747e-01 5.855268e-02
    ## [12996]          NaN 2.888092e-01 6.774433e-02 6.061088e-01 6.214809e-05
    ## [13001] 5.772319e-01          NaN 1.118689e-05 2.365450e-01 5.597274e-01
    ## [13006] 2.470038e-01 5.813631e-01 3.837202e-03 1.063389e-03 6.210553e-01
    ## [13011] 7.310669e-01 5.083232e-03          NaN 4.057526e-01          NaN
    ## [13016] 4.129682e-01 4.151768e-03 3.204838e-01 1.369974e-01          NaN
    ## [13021] 3.204838e-01 1.515603e-01 4.394010e-01          NaN          NaN
    ## [13026] 7.700620e-01 9.205701e-01 4.661356e-01 4.126992e-04          NaN
    ## [13031]          NaN          NaN 5.241821e-01 5.268932e-10 7.277365e-01
    ## [13036] 1.957657e-01          NaN          NaN 3.691313e-03 1.686364e-01
    ## [13041] 2.166908e-01 1.970672e-02          NaN 7.536805e-01          NaN
    ## [13046]          NaN          NaN 7.766050e-01 3.204838e-01 3.617408e-01
    ## [13051]          NaN 1.319658e-02 4.534517e-02 8.287089e-01 2.342015e-01
    ## [13056] 1.003099e-02          NaN 2.801584e-06 5.762176e-05          NaN
    ## [13061]          NaN 4.732545e-01 7.865103e-02 3.732431e-01 3.073894e-44
    ## [13066] 6.969613e-12 5.128050e-01 7.618909e-03          NaN          NaN
    ## [13071] 2.452998e-02          NaN 5.809576e-04 2.524120e-01 3.204838e-01
    ## [13076] 2.754024e-02          NaN 5.711137e-01 1.350344e-16          NaN
    ## [13081] 1.051311e-16 1.381282e-02 3.204838e-01 1.949816e-02 1.323882e-02
    ## [13086]          NaN          NaN 8.772701e-01 1.397928e-15 1.277130e-04
    ## [13091]          NaN 3.466838e-01 5.723958e-02          NaN          NaN
    ## [13096]          NaN 2.929970e-05 8.996545e-03 3.593027e-03 2.590928e-01
    ## [13101] 7.199614e-01 1.279283e-01 8.886739e-01 1.213695e-03 1.093994e-01
    ## [13106] 8.686743e-01 3.713616e-07 4.357878e-12 7.874234e-01 9.123314e-06
    ## [13111] 3.425552e-01 2.587454e-01 4.243300e-03 4.853670e-03 3.204838e-01
    ## [13116] 1.508161e-01 2.324432e-01 7.625327e-01          NaN 7.867028e-05
    ## [13121] 3.648086e-01 2.418436e-02          NaN 1.002308e-19          NaN
    ## [13126] 7.113519e-01 4.372721e-14 8.351324e-02 9.222084e-01 9.831590e-01
    ## [13131] 9.851495e-01 2.499486e-01 9.263271e-02 5.193965e-14 5.067624e-01
    ## [13136] 9.815569e-01 2.072730e-01 3.145964e-01 5.597411e-02 3.361472e-01
    ## [13141] 3.580088e-01 3.161211e-05          NaN          NaN 7.012461e-01
    ## [13146]          NaN          NaN 5.089927e-04          NaN          NaN
    ## [13151]          NaN          NaN 1.120534e-03          NaN 3.078265e-01
    ## [13156]          NaN 4.268968e-01          NaN 1.796134e-01 2.740959e-06
    ## [13161] 4.984116e-02 6.651771e-01 2.201482e-01 1.790922e-01 5.324590e-01
    ## [13166] 1.512889e-04 9.860014e-01          NaN 7.577143e-01 6.103376e-06
    ## [13171]          NaN 2.607747e-01          NaN 4.443851e-03          NaN
    ## [13176] 8.447508e-01 1.845839e-02          NaN 1.309166e-02 4.124736e-01
    ## [13181] 5.597669e-01 2.055626e-01          NaN 3.285715e-14 2.158897e-02
    ## [13186] 3.042090e-01 3.606184e-03 1.075607e-21 6.737345e-01 7.716288e-04
    ## [13191] 5.675897e-01 2.063165e-02          NaN 3.121129e-05 3.085292e-04
    ## [13196] 4.740081e-01 1.318546e-01 8.821771e-02 7.063786e-01          NaN
    ## [13201] 8.453493e-03          NaN 4.617811e-01 1.088538e-06 9.739309e-01
    ## [13206] 7.049446e-01 1.489110e-12 6.707773e-11 2.416119e-22 6.776866e-01
    ## [13211] 2.411392e-02 2.444267e-05          NaN          NaN 4.348444e-02
    ## [13216] 7.177705e-01          NaN 3.204838e-01 1.391764e-01 2.169451e-01
    ## [13221]          NaN          NaN 7.098770e-01          NaN 6.730852e-01
    ## [13226] 9.172754e-03 1.078542e-03 2.555785e-01 4.044686e-04 4.285551e-01
    ## [13231] 8.694312e-01          NaN 3.204838e-01 1.427788e-02 2.751516e-03
    ## [13236] 5.933214e-03 7.232273e-01 3.101717e-09 9.891258e-01 1.139913e-01
    ## [13241] 4.627412e-02 2.196145e-22 4.545839e-01          NaN 7.692887e-07
    ## [13246] 9.862571e-02          NaN 9.523756e-01          NaN 2.875269e-01
    ## [13251] 2.141534e-01 3.372837e-03 1.696425e-18 4.154868e-01 4.829018e-01
    ## [13256] 2.720893e-02 5.538162e-02 9.817744e-02          NaN          NaN
    ## [13261] 8.037960e-01 2.537119e-07 3.638712e-06          NaN          NaN
    ## [13266] 1.217921e-09 9.798360e-01 2.175789e-26 7.190919e-01 2.120464e-01
    ## [13271] 1.412619e-11 1.971672e-01 1.561914e-01          NaN 9.207413e-32
    ## [13276]          NaN 4.168469e-01 4.093930e-01 6.228702e-01 5.541288e-07
    ## [13281]          NaN 3.883410e-01 5.198455e-01 3.204838e-01 5.231171e-02
    ## [13286] 7.806940e-01 2.892105e-02 5.422416e-01 6.476773e-02 4.089622e-03
    ## [13291] 4.593533e-01          NaN 2.750965e-02 6.120350e-01 1.560481e-03
    ## [13296] 1.467655e-33 2.037269e-01 1.864528e-23 1.684787e-07 1.371455e-13
    ## [13301] 4.021694e-01 6.342419e-05 3.167857e-03          NaN 4.330482e-02
    ## [13306]          NaN 3.503531e-06 1.076860e-02 5.675517e-01 8.102569e-02
    ## [13311]          NaN 5.589400e-01 7.449836e-08          NaN 9.399975e-01
    ## [13316] 2.364846e-04          NaN 3.593011e-01 6.487795e-01 1.276071e-01
    ## [13321] 2.551354e-01          NaN 3.204838e-01 1.569886e-01          NaN
    ## [13326] 5.205804e-09          NaN 8.393944e-03 3.174157e-03 1.059541e-03
    ## [13331]          NaN          NaN 3.761106e-03 2.714229e-01 4.753588e-05
    ## [13336] 6.994074e-01 2.216446e-01 5.650767e-04 9.157766e-01 2.069548e-02
    ## [13341]          NaN          NaN 8.590992e-01 2.361938e-09 3.204838e-01
    ## [13346] 4.677970e-01 1.739334e-13 3.204838e-01 3.486932e-01 1.566136e-01
    ## [13351] 4.629002e-01 9.327781e-04 2.885597e-09 1.364266e-01 2.947112e-06
    ## [13356] 8.379785e-01          NaN 6.583448e-01 8.099927e-01 3.866661e-04
    ## [13361]          NaN          NaN 5.534569e-14 7.294567e-02          NaN
    ## [13366]          NaN 8.056937e-07 7.913212e-01 5.219365e-01 6.158631e-06
    ## [13371] 9.494106e-14 1.875013e-03 3.208920e-02 9.383955e-02 3.321832e-02
    ## [13376] 7.239014e-01          NaN 2.330888e-44 6.208902e-01 8.761471e-17
    ## [13381] 2.473247e-01          NaN 1.899256e-02 7.007817e-01 5.767992e-04
    ## [13386] 9.462767e-20 6.278415e-01 7.537036e-03 6.470341e-01 2.121686e-60
    ## [13391] 8.888517e-01 4.954425e-02 4.893721e-01 3.796514e-05 2.796055e-01
    ## [13396]          NaN 9.150547e-01 1.228416e-22 2.011238e-04          NaN
    ## [13401] 2.749475e-01 2.067113e-12 4.765595e-01          NaN          NaN
    ## [13406] 3.272978e-01 3.607758e-01 5.521945e-03 1.459878e-03 2.834654e-02
    ## [13411] 3.646052e-01 9.280260e-01 3.204838e-01 7.671292e-01 8.309974e-10
    ## [13416] 6.227891e-02 1.548848e-02 8.982646e-11          NaN 1.283825e-01
    ## [13421]          NaN 8.315658e-01 9.516349e-05 3.204838e-01 7.098674e-03
    ## [13426] 6.677588e-09 1.955441e-23 8.511148e-01 9.077435e-01 5.354671e-02
    ## [13431] 3.976945e-02 1.901127e-01 1.209579e-04          NaN          NaN
    ## [13436] 7.312764e-01 5.256370e-03 5.597769e-01 1.234378e-31 2.088487e-06
    ## [13441] 5.767296e-01 8.340949e-04 7.950773e-10          NaN          NaN
    ## [13446] 1.817668e-03 2.277561e-01 3.204838e-01 3.328495e-06          NaN
    ## [13451] 6.745808e-01 2.363634e-73 9.725892e-01 8.501811e-01 1.129322e-03
    ## [13456] 3.032308e-03 4.118461e-19 1.640134e-01          NaN 4.212566e-01
    ## [13461] 8.566000e-31 5.727633e-11 3.013009e-04 3.204838e-01 7.170393e-01
    ## [13466] 9.136190e-11 5.437191e-02 1.089259e-10 3.740178e-01 2.301129e-01
    ## [13471] 8.786828e-01 1.551713e-15 3.654130e-02 7.200330e-01 2.328560e-01
    ## [13476] 8.937585e-03 2.814841e-07 7.154536e-01 3.816493e-01 3.336643e-03
    ## [13481] 2.168794e-03 3.865717e-01 7.671318e-01 5.328639e-02 1.378760e-07
    ## [13486] 4.743278e-01 1.233334e-06 8.924130e-01 1.081766e-04 1.331012e-07
    ## [13491] 3.236860e-03 5.921300e-06 8.924305e-01 2.078302e-02 6.983929e-02
    ## [13496]          NaN 9.101858e-01 5.390553e-01          NaN 9.205910e-01
    ## [13501] 5.639206e-03 3.204838e-01          NaN 7.109109e-02 1.717436e-31
    ## [13506] 6.772928e-01 1.065330e-01 9.097289e-02 6.698685e-01 3.556475e-22
    ## [13511] 4.507864e-12          NaN          NaN          NaN 2.135159e-01
    ## [13516] 2.936533e-01 2.782201e-02 6.680733e-02 4.429378e-01 1.920911e-01
    ## [13521] 2.443175e-01 1.632000e-03 1.045288e-01          NaN 6.573329e-02
    ## [13526] 4.792012e-05 3.818801e-01          NaN 3.639621e-01 2.377020e-01
    ## [13531] 1.768002e-04 2.331125e-01 2.424665e-01 3.876763e-03          NaN
    ## [13536]          NaN 5.427361e-03 1.044188e-01          NaN 2.591927e-01
    ## [13541] 1.985740e-05 2.931180e-21 3.204838e-01 5.886590e-10 5.882173e-01
    ## [13546] 9.182089e-02 3.751801e-01 2.481624e-19 6.495325e-02 2.958922e-12
    ## [13551]          NaN 3.204838e-01 5.522029e-01          NaN 6.604594e-01
    ## [13556] 2.328898e-08 1.663371e-04 1.610338e-01 1.245575e-08 3.204838e-01
    ## [13561]          NaN          NaN 2.140943e-02 1.591788e-01 6.343959e-01
    ## [13566] 3.204838e-01          NaN          NaN 4.222291e-01          NaN
    ## [13571] 8.658751e-01 5.096703e-01 5.551662e-01 3.188198e-01 5.289854e-01
    ## [13576] 1.370782e-01 8.149400e-01 3.384624e-09 2.075546e-04 2.711992e-05
    ## [13581] 5.276549e-06 2.173785e-04 3.164171e-11 3.156401e-02          NaN
    ## [13586] 2.971070e-03 6.769817e-01          NaN          NaN 9.064356e-01
    ## [13591] 1.889979e-02 6.186046e-01 7.227369e-01 1.583483e-08 2.307353e-14
    ## [13596] 3.204838e-01          NaN 2.536962e-06 9.029907e-01          NaN
    ## [13601] 9.667827e-01 1.491817e-01 1.387008e-01          NaN 4.522667e-06
    ## [13606]          NaN 9.853434e-01 6.010327e-01          NaN 4.718311e-01
    ## [13611] 1.649572e-05 3.596967e-38 7.068626e-01 2.836125e-01 2.498233e-01
    ## [13616]          NaN 6.697248e-01 7.416538e-04 1.607501e-10 2.762452e-02
    ## [13621]          NaN 6.177611e-01 1.681588e-06 1.742266e-26 1.945285e-04
    ## [13626] 3.171687e-10 2.639170e-20          NaN 8.241298e-01          NaN
    ## [13631] 1.707239e-04 3.860045e-03 7.969782e-01 6.049207e-01 1.429115e-01
    ## [13636] 1.875872e-05 2.311135e-02 5.882095e-01 1.234305e-03 1.339855e-11
    ## [13641] 2.460541e-11 1.756748e-05 3.204838e-01 9.824499e-02 2.804370e-09
    ## [13646]          NaN          NaN 8.680498e-01          NaN 1.985163e-01
    ## [13651] 3.204838e-01          NaN          NaN          NaN          NaN
    ## [13656]          NaN          NaN          NaN          NaN          NaN
    ## [13661]          NaN          NaN          NaN          NaN          NaN
    ## [13666]          NaN          NaN          NaN          NaN          NaN
    ## [13671]          NaN          NaN          NaN          NaN          NaN
    ## [13676]          NaN          NaN 3.204838e-01          NaN          NaN
    ## [13681]          NaN          NaN          NaN          NaN          NaN
    ## [13686]          NaN          NaN          NaN          NaN          NaN
    ## [13691] 3.204838e-01          NaN          NaN          NaN          NaN
    ## [13696]          NaN          NaN 3.204838e-01          NaN          NaN
    ## [13701]          NaN          NaN          NaN          NaN          NaN
    ## [13706]          NaN          NaN          NaN          NaN          NaN
    ## [13711]          NaN          NaN          NaN          NaN          NaN
    ## [13716]          NaN 1.577219e-01          NaN          NaN          NaN
    ## [13721]          NaN          NaN          NaN 3.204838e-01          NaN
    ## [13726]          NaN          NaN          NaN          NaN          NaN
    ## [13731]          NaN          NaN 3.204838e-01 3.204838e-01          NaN
    ## [13736]          NaN          NaN          NaN          NaN          NaN
    ## [13741]          NaN          NaN          NaN          NaN          NaN
    ## [13746]          NaN          NaN          NaN          NaN          NaN
    ## [13751]          NaN          NaN          NaN          NaN          NaN
    ## [13756]          NaN          NaN          NaN          NaN          NaN
    ## [13761]          NaN          NaN          NaN          NaN          NaN
    ## [13766]          NaN          NaN          NaN          NaN          NaN
    ## [13771]          NaN          NaN          NaN          NaN          NaN
    ## [13776] 3.204838e-01 3.204838e-01          NaN 3.204838e-01          NaN
    ## [13781]          NaN          NaN          NaN          NaN          NaN
    ## [13786]          NaN          NaN          NaN          NaN          NaN
    ## [13791]          NaN          NaN          NaN          NaN          NaN
    ## [13796]          NaN          NaN          NaN          NaN 5.669176e-01
    ## [13801] 3.204838e-01          NaN          NaN          NaN          NaN
    ## [13806]          NaN 3.204838e-01          NaN          NaN          NaN
    ## [13811]          NaN          NaN          NaN          NaN          NaN
    ## [13816]          NaN          NaN          NaN          NaN          NaN
    ## [13821]          NaN          NaN          NaN          NaN          NaN
    ## [13826]          NaN          NaN 3.204838e-01          NaN          NaN
    ## [13831]          NaN          NaN          NaN          NaN          NaN
    ## [13836]          NaN          NaN          NaN          NaN          NaN
    ## [13841]          NaN          NaN          NaN          NaN          NaN
    ## [13846]          NaN          NaN          NaN          NaN          NaN
    ## [13851]          NaN          NaN          NaN          NaN          NaN
    ## [13856]          NaN          NaN          NaN          NaN          NaN
    ## [13861]          NaN          NaN          NaN          NaN          NaN
    ## [13866]          NaN          NaN          NaN 9.394044e-11 4.718243e-05
    ## [13871] 1.287013e-02 1.894854e-01 3.177988e-02 2.714328e-02 2.014461e-02
    ## [13876] 6.928040e-01 1.429340e-01 3.204838e-01 1.502758e-15 1.437338e-02
    ## [13881]          NaN          NaN 8.313581e-01          NaN 1.104647e-03
    ## [13886] 1.013321e-11 1.585853e-01 2.069154e-01          NaN 1.645631e-05
    ## [13891] 5.412566e-03 8.737484e-02 8.335893e-02 1.662719e-01 1.062418e-01
    ## [13896] 2.987286e-01 9.273683e-01 4.802576e-01 7.843532e-01 1.324264e-03
    ## [13901] 4.134871e-01 2.862520e-01 3.589930e-01 2.962479e-06 3.417390e-02
    ## [13906] 6.761387e-01 7.940271e-01 1.202037e-03 2.265057e-01 3.711627e-01
    ## [13911] 3.066548e-04 1.960067e-09 3.955508e-01 5.344484e-01 7.178955e-01
    ## [13916] 3.687788e-02 2.740406e-01 1.153478e-02 5.080347e-04 2.184687e-04
    ## [13921] 2.715630e-03 2.580127e-01 3.045435e-01 1.038798e-01 7.351997e-01
    ## [13926] 7.251565e-01 1.529803e-01 1.343331e-01 7.708153e-02 1.628341e-01
    ## [13931] 3.204838e-01 2.354509e-02 4.000156e-02 3.186839e-06 7.807480e-01
    ## [13936]          NaN 1.573410e-02          NaN 6.476628e-01 1.309685e-01
    ## [13941] 3.228360e-01 3.204838e-01 2.068606e-02 9.362672e-01 3.822878e-02
    ## [13946] 1.631652e-20 9.979225e-01 2.267221e-01          NaN 1.594044e-01
    ## [13951] 3.494274e-03 8.418067e-01 2.888525e-02          NaN          NaN
    ## [13956] 1.280714e-02 5.267692e-01          NaN 9.380449e-05 5.770381e-02
    ## [13961] 2.406623e-01 2.074866e-10 9.760309e-01 8.359201e-01 8.248641e-27
    ## [13966]          NaN          NaN 2.728826e-20 2.179509e-01 2.743859e-01
    ## [13971] 6.216510e-01 6.497992e-14 9.791700e-01 1.079865e-20 9.599144e-01
    ## [13976] 9.144149e-01          NaN          NaN 6.201478e-01 9.621560e-01
    ## [13981] 3.302550e-09 1.625689e-11 9.047669e-07 4.011354e-01 4.635827e-02
    ## [13986] 7.782889e-02 5.156014e-25          NaN 4.862072e-01 4.720373e-01
    ## [13991] 2.860034e-12 1.447860e-05 1.929215e-07 4.403271e-01          NaN
    ## [13996] 3.273737e-04 3.204838e-01          NaN          NaN 1.355003e-19
    ## [14001] 1.180032e-02 7.171858e-04 2.701789e-02 8.887590e-01          NaN
    ## [14006] 1.198756e-02 1.981906e-02 6.651346e-01 1.324001e-03 4.177413e-01
    ## [14011] 3.039295e-02          NaN 3.113445e-10          NaN 8.144936e-03
    ## [14016] 4.833518e-01 7.701237e-01          NaN 5.391596e-01 1.677579e-02
    ## [14021] 3.296133e-05 2.115806e-01          NaN 9.451639e-01 1.516631e-12
    ## [14026] 1.188592e-01          NaN 1.586228e-01 2.638005e-05 3.682863e-01
    ## [14031] 9.146997e-01 7.234099e-03 4.670881e-01 4.557610e-01 3.348502e-02
    ## [14036] 5.975733e-01          NaN 1.542005e-03 2.749042e-01 1.063475e-05
    ## [14041] 2.436599e-01 8.075913e-05 5.584091e-01 8.229848e-03          NaN
    ## [14046] 3.204838e-01 9.680885e-02 4.342218e-02          NaN 2.872735e-01
    ## [14051] 2.780259e-01 2.274787e-13 2.958804e-01 3.555355e-02 6.482858e-01
    ## [14056] 5.105095e-01 8.543654e-01 7.520779e-01 3.583130e-01          NaN
    ## [14061] 4.107837e-03 2.482411e-01 1.674988e-03 9.592679e-03          NaN
    ## [14066] 3.204838e-01 9.357832e-13          NaN 1.891005e-01 9.135651e-01
    ## [14071] 4.705515e-05 2.537054e-01 3.779787e-01 5.239803e-02 3.024617e-01
    ## [14076] 9.539392e-01 1.059828e-01 6.728069e-01 2.016569e-03 2.633869e-04
    ## [14081]          NaN 8.897594e-01 3.138902e-01          NaN 7.547700e-01
    ## [14086] 4.170107e-02 2.334891e-01 5.965826e-01 4.764911e-01 1.234508e-01
    ## [14091] 6.486165e-02 3.094580e-03 3.419018e-19          NaN 3.480503e-01
    ## [14096] 3.239302e-04 2.700417e-08 1.048539e-03 1.758826e-01 1.356455e-01
    ## [14101] 7.612821e-34 5.649566e-01 2.111377e-01 4.225584e-02 2.159719e-03
    ## [14106] 2.598328e-01 7.572833e-01 3.547415e-03 1.269355e-01 1.823013e-01
    ## [14111] 1.407571e-04 3.204838e-01          NaN 2.852446e-03 5.578091e-01
    ## [14116] 7.418142e-12 1.221522e-08 8.665283e-11 5.925211e-03 3.041389e-01
    ## [14121] 1.512885e-02 6.603261e-01 1.158929e-04 1.907359e-01 1.963090e-01
    ## [14126] 4.627514e-09 2.749970e-03 8.476843e-05 2.886117e-01 3.204838e-01
    ## [14131] 3.149899e-01 5.217428e-01 2.240234e-01          NaN 9.277266e-01
    ## [14136] 7.860895e-02 1.697938e-02 7.175138e-01 1.648299e-02 5.392928e-23
    ## [14141] 4.778927e-01 9.696747e-01 6.937225e-01 3.507584e-01          NaN
    ## [14146] 7.847912e-03 7.178758e-01 9.699322e-01 6.258387e-01          NaN
    ## [14151]          NaN          NaN 3.468724e-01 7.528680e-01 6.342097e-01
    ## [14156] 7.449304e-01 1.642946e-06 5.372702e-09 4.806983e-04          NaN
    ## [14161] 2.019040e-04 5.106824e-02 7.901965e-01          NaN 5.359617e-02
    ## [14166] 1.462489e-28 1.571788e-01 6.683093e-01 9.204435e-02          NaN
    ## [14171] 5.525747e-01 8.529666e-01 3.389627e-01 2.842247e-01 3.198259e-01
    ## [14176]          NaN 2.412268e-07 3.877196e-01 4.051545e-01          NaN
    ## [14181] 9.180933e-01          NaN 9.271053e-01 1.666503e-07 4.760224e-01
    ## [14186] 1.650089e-05          NaN 9.361314e-01          NaN 2.510991e-01
    ## [14191] 3.673081e-06 3.204838e-01 1.598487e-04 2.015940e-01 3.204838e-01
    ## [14196] 9.721622e-01 2.149761e-01 1.722763e-01          NaN          NaN
    ## [14201] 3.204838e-01          NaN          NaN 3.204838e-01          NaN
    ## [14206] 3.610621e-06          NaN 3.428024e-01 1.413671e-02          NaN
    ## [14211]          NaN          NaN 8.408001e-01 7.639459e-11 8.276471e-01
    ## [14216]          NaN 1.258893e-11 9.631613e-03 1.672106e-13 1.881594e-02
    ## [14221] 5.205453e-02 6.430584e-01 1.025463e-01 3.428523e-01 1.924457e-01
    ## [14226] 1.562738e-04 9.517091e-01          NaN 6.838028e-01          NaN
    ## [14231] 4.177710e-04          NaN          NaN 8.582268e-01          NaN
    ## [14236] 3.465422e-01 3.498034e-04 6.141314e-01 5.404385e-01 1.837418e-01
    ## [14241] 1.440613e-02 4.329054e-02 3.764298e-03          NaN 1.123035e-01
    ## [14246] 9.826512e-03 5.180277e-01 6.359011e-01 2.892014e-07 5.908210e-02
    ## [14251] 9.614700e-04          NaN          NaN          NaN          NaN
    ## [14256]          NaN 8.660273e-02 4.189075e-02 9.270838e-01 8.427713e-01
    ## [14261] 1.214618e-01 9.954228e-03 3.786763e-01 3.204838e-01          NaN
    ## [14266] 6.926288e-01 8.034278e-02 2.629289e-07 3.127886e-01          NaN
    ## [14271] 3.371092e-05          NaN          NaN 5.961838e-01 3.204838e-01
    ## [14276]          NaN 9.658123e-01 2.323280e-01 9.661387e-01 7.619410e-10
    ## [14281] 1.179787e-03 4.702961e-03 5.755931e-01          NaN 4.290724e-01
    ## [14286] 8.326659e-03          NaN 2.360868e-01          NaN          NaN
    ## [14291] 3.739206e-01          NaN 2.918747e-02          NaN 1.125237e-01
    ## [14296] 5.763086e-01 1.994889e-01 1.561208e-01 5.942595e-01 7.506324e-09
    ## [14301] 4.255832e-12 1.755391e-07 2.928332e-02 2.890901e-07 3.705073e-04
    ## [14306] 7.755673e-01 2.689980e-01 2.329722e-15 2.480837e-02 8.549803e-02
    ## [14311] 5.831322e-01 2.313201e-06 7.364135e-04 7.281240e-01 7.180621e-03
    ## [14316] 4.933810e-01          NaN 2.817720e-01 3.403065e-15 7.938239e-03
    ## [14321] 6.241771e-05 5.006703e-08 3.032997e-01 3.204838e-01 7.222327e-01
    ## [14326] 2.135739e-01 6.457957e-05          NaN 9.041125e-01 1.577219e-01
    ## [14331] 3.380887e-01 3.684008e-02 9.909257e-01          NaN 1.136784e-01
    ## [14336] 4.848189e-01 3.295380e-01          NaN          NaN 1.793991e-01
    ## [14341] 1.616576e-11 1.154648e-02 6.705808e-01 2.777921e-02 3.122864e-01
    ## [14346] 3.204838e-01          NaN 8.412821e-01 2.216492e-01          NaN
    ## [14351] 2.293117e-01 2.673589e-02 2.997686e-01          NaN 2.017965e-01
    ## [14356] 5.738475e-06          NaN          NaN 2.029930e-01          NaN
    ## [14361] 3.759106e-09          NaN 3.879437e-03 2.683460e-03 3.879756e-04
    ## [14366] 1.138387e-02 1.029550e-02 6.004136e-02          NaN 6.709904e-01
    ## [14371] 7.909126e-01          NaN          NaN          NaN 3.348357e-01
    ## [14376] 8.987110e-02 4.898802e-25          NaN 9.124901e-02 6.791090e-01
    ## [14381]          NaN 4.377590e-01 4.718843e-03 3.204838e-01 1.083990e-01
    ## [14386] 5.661289e-02 6.459870e-01 6.375712e-01 1.497804e-01 4.852438e-01
    ## [14391] 5.420187e-01 5.520615e-03 2.505498e-01 3.204838e-01 5.409076e-01
    ## [14396] 1.077465e-03 9.380106e-01 1.070678e-05 1.598146e-02 3.513884e-03
    ## [14401] 4.883957e-01 3.204838e-01 2.615728e-04 1.307261e-02 1.075092e-48
    ## [14406] 1.112409e-02 9.046220e-01 4.006716e-42          NaN 1.151318e-02
    ## [14411] 1.280848e-03 3.483111e-04 7.978616e-02 4.606496e-01 2.005284e-01
    ## [14416] 8.887675e-01 1.256037e-01          NaN          NaN 1.559039e-01
    ## [14421]          NaN 5.117381e-01          NaN 3.635728e-22 4.514683e-06
    ## [14426] 4.307145e-01 3.204838e-01          NaN 2.694983e-01          NaN
    ## [14431]          NaN 4.574028e-01 2.447389e-01          NaN 6.322099e-03
    ## [14436]          NaN 7.193604e-01 1.215570e-04 2.155485e-01 1.230020e-01
    ## [14441] 1.568656e-01 3.885015e-01 1.600811e-04 8.688846e-01          NaN
    ## [14446] 1.108365e-05          NaN          NaN 1.812365e-01 6.259316e-01
    ## [14451]          NaN 8.922693e-01 6.122406e-03 2.247772e-03 3.891984e-01
    ## [14456] 1.773148e-01          NaN 6.095135e-04 1.246610e-09 5.129900e-01
    ## [14461] 9.000189e-01 3.177935e-07 5.038425e-01 2.992463e-02 5.543792e-01
    ## [14466] 3.204838e-01          NaN 1.236697e-01          NaN 4.579643e-09
    ## [14471] 2.042371e-02 5.194987e-05 1.259530e-01 1.595718e-01 8.611349e-02
    ## [14476] 5.426553e-01 1.615214e-09 8.976395e-01 2.061939e-01 1.746708e-08
    ## [14481] 1.346892e-12 2.454758e-01          NaN 5.200659e-01 4.242756e-02
    ## [14486] 2.517803e-01          NaN          NaN          NaN          NaN
    ## [14491]          NaN          NaN 1.669974e-39 1.646849e-02 3.204838e-01
    ## [14496] 2.140489e-01 3.633112e-04 1.772964e-06 7.743179e-01 4.055943e-02
    ## [14501] 3.204838e-01 5.698385e-01 1.721504e-02 8.105893e-01 3.625801e-02
    ## [14506] 8.390420e-02 8.022711e-18 6.859996e-01 3.165610e-02 1.711319e-03
    ## [14511]          NaN 1.243114e-02 1.176334e-02 1.879139e-01 1.435600e-03
    ## [14516] 4.861753e-01          NaN 6.733942e-01          NaN          NaN
    ## [14521]          NaN          NaN          NaN          NaN          NaN
    ## [14526]          NaN          NaN          NaN          NaN 3.204838e-01
    ## [14531] 3.204838e-01          NaN          NaN          NaN          NaN
    ## [14536]          NaN          NaN          NaN          NaN          NaN
    ## [14541]          NaN          NaN          NaN          NaN          NaN
    ## [14546] 3.204838e-01          NaN          NaN          NaN          NaN
    ## [14551]          NaN          NaN          NaN          NaN          NaN
    ## [14556]          NaN          NaN          NaN          NaN          NaN
    ## [14561]          NaN          NaN          NaN          NaN          NaN
    ## [14566]          NaN          NaN          NaN          NaN          NaN
    ## [14571]          NaN          NaN          NaN          NaN          NaN
    ## [14576] 6.220815e-01          NaN          NaN 8.338780e-01          NaN
    ## [14581]          NaN          NaN          NaN 3.204838e-01          NaN
    ## [14586]          NaN          NaN          NaN          NaN          NaN
    ## [14591]          NaN          NaN 1.163934e-01          NaN          NaN
    ## [14596]          NaN          NaN          NaN 3.204838e-01          NaN
    ## [14601]          NaN          NaN          NaN          NaN          NaN
    ## [14606]          NaN          NaN          NaN          NaN          NaN
    ## [14611]          NaN          NaN          NaN          NaN          NaN
    ## [14616]          NaN          NaN          NaN          NaN          NaN
    ## [14621]          NaN          NaN          NaN          NaN          NaN
    ## [14626] 1.593236e-01          NaN          NaN          NaN          NaN
    ## [14631]          NaN          NaN          NaN          NaN          NaN
    ## [14636]          NaN          NaN          NaN          NaN          NaN
    ## [14641]          NaN          NaN          NaN          NaN          NaN
    ## [14646]          NaN          NaN          NaN          NaN          NaN
    ## [14651]          NaN          NaN          NaN          NaN          NaN
    ## [14656]          NaN 6.547443e-02 1.244780e-01          NaN 3.204838e-01
    ## [14661] 1.658431e-05 9.019321e-01 3.616033e-04 7.567165e-06          NaN
    ## [14666]          NaN 1.400713e-01 4.182004e-02 6.975144e-01 9.228576e-01
    ## [14671] 2.426451e-01 1.619968e-25          NaN 1.880052e-02          NaN
    ## [14676] 9.451471e-01 9.881740e-02 2.474955e-03 5.881402e-01 2.292887e-02
    ## [14681] 7.084889e-02 2.433566e-01 2.111445e-02 6.493215e-01 1.935979e-01
    ## [14686] 9.443296e-01 8.004746e-01 4.099248e-03 4.492594e-35          NaN
    ## [14691] 7.760180e-03 5.664620e-01 2.590519e-01 1.896333e-02          NaN
    ## [14696] 3.204838e-01 5.178407e-01 5.407809e-02          NaN 1.644203e-02
    ## [14701] 1.541121e-06          NaN          NaN 1.516623e-02 2.483769e-01
    ## [14706] 7.591516e-01 7.106147e-04 7.117086e-02 3.204838e-01 1.811751e-13
    ## [14711] 8.331644e-02          NaN 2.844215e-01 1.432056e-04 9.556583e-01
    ## [14716] 4.309106e-06          NaN 5.393257e-01          NaN          NaN
    ## [14721]          NaN 5.269819e-02 1.909327e-02 8.003426e-01 2.822872e-01
    ## [14726] 3.204838e-01 3.204838e-01 7.660146e-01 4.175736e-01          NaN
    ## [14731] 1.083501e-03 1.492839e-01 3.204838e-01 1.864774e-13 1.304094e-01
    ## [14736] 6.085024e-01 6.152867e-01 9.760495e-01 6.036224e-01 1.070029e-06
    ## [14741]          NaN 1.031749e-16 9.959902e-01 6.036489e-01 4.156608e-01
    ## [14746] 2.549607e-02 1.482341e-16 1.268156e-01 1.768100e-01 3.759083e-02
    ## [14751] 8.091026e-03 7.621542e-01 9.859112e-01 3.371811e-07 3.027217e-02
    ## [14756] 1.262002e-06 3.992821e-04          NaN 9.414512e-01 2.355055e-01
    ## [14761] 7.404103e-02          NaN          NaN 2.129080e-04 2.197306e-01
    ## [14766] 1.059393e-01 1.956746e-01 6.362131e-01 2.363482e-04 1.510627e-01
    ## [14771] 2.272519e-10 1.944878e-06 1.912180e-03 2.392978e-01 5.712535e-01
    ## [14776] 7.761612e-02 1.643753e-04 9.750525e-01 3.755833e-01          NaN
    ## [14781] 4.547168e-33 4.794019e-01 3.744485e-02 3.538985e-03 2.320369e-09
    ## [14786]          NaN 3.218718e-02 3.050717e-01          NaN 6.654780e-05
    ## [14791] 6.343549e-02 1.079624e-01 1.601565e-01 1.335054e-01 3.017576e-01
    ## [14796] 1.339392e-02 8.718628e-01 3.204838e-01 2.450341e-01 1.234245e-05
    ## [14801] 2.056171e-55          NaN          NaN 1.894618e-05 8.018055e-04
    ## [14806] 4.564845e-02 4.907798e-03          NaN 3.894100e-06 4.782651e-01
    ## [14811] 6.265122e-01 4.144373e-01 6.536622e-01 8.872456e-01          NaN
    ## [14816] 3.204838e-01 1.313557e-13 2.723363e-15          NaN 7.866573e-01
    ## [14821] 1.538833e-50 4.526613e-01 1.008404e-11 5.700523e-07 3.204838e-01
    ## [14826]          NaN 3.052307e-01 2.463895e-01 6.190125e-01 1.264584e-02
    ## [14831] 7.069241e-10 4.622528e-01          NaN 2.238159e-09 1.168731e-07
    ## [14836]          NaN 5.238943e-01 1.271461e-01 2.050750e-01 3.220001e-03
    ## [14841] 1.881496e-02          NaN 3.294494e-36 7.848400e-01 2.281306e-01
    ## [14846]          NaN 3.285405e-01 3.125573e-47 3.907280e-02 2.962738e-01
    ## [14851] 3.846858e-02 1.639587e-01          NaN 4.720726e-05 5.775583e-01
    ## [14856] 6.657451e-01          NaN 2.042188e-01          NaN 1.852858e-01
    ## [14861] 6.914378e-03 1.440414e-01 7.898988e-01 3.512928e-01 1.574619e-01
    ## [14866] 3.957698e-04          NaN 1.534913e-01 7.369648e-01          NaN
    ## [14871] 4.127754e-01          NaN 3.204838e-01 3.204838e-01 3.204838e-01
    ## [14876] 2.445987e-01          NaN 2.711187e-05          NaN 7.697891e-02
    ## [14881] 4.534473e-01 7.661235e-01 1.786775e-01          NaN          NaN
    ## [14886] 1.583101e-01 2.856296e-01 5.427149e-01          NaN 5.218169e-01
    ## [14891] 6.085287e-01          NaN 6.562980e-01 1.955583e-04 3.503027e-02
    ## [14896] 4.267231e-12 6.366175e-11 1.047713e-11 4.294098e-02 5.705806e-02
    ## [14901] 1.788108e-01 8.197817e-04 5.444385e-02 3.255815e-01 1.165591e-02
    ## [14906] 2.347379e-01 4.685362e-02          NaN 1.622876e-01 4.069460e-02
    ## [14911] 6.780482e-03 2.819415e-01 4.468896e-01 5.173769e-04 3.935609e-01
    ## [14916] 7.162991e-01 2.595695e-01 7.415798e-02 6.659036e-01 6.699695e-02
    ## [14921] 5.062014e-01 1.121843e-18 2.511627e-03 3.562303e-02 8.837071e-01
    ## [14926]          NaN          NaN 3.204838e-01 9.089741e-02 1.513227e-01
    ## [14931] 2.538229e-01          NaN 5.452602e-01 7.147568e-01 1.553393e-02
    ## [14936] 2.174562e-03          NaN          NaN 1.564744e-01 6.415326e-01
    ## [14941] 1.176377e-08 9.311844e-02 7.339745e-01 1.085565e-05 2.001009e-01
    ## [14946] 1.172015e-01 5.568197e-01 1.166175e-03 3.204838e-01 1.613411e-01
    ## [14951] 7.946242e-03 1.696153e-01 3.726286e-04 2.499780e-05 8.852297e-01
    ## [14956] 3.177495e-02 8.533960e-01 3.571057e-02 3.204838e-01 8.380694e-01
    ## [14961] 9.777177e-04 2.551007e-03 7.128021e-59 1.742314e-02 8.333764e-01
    ## [14966]          NaN          NaN 1.231028e-03 9.561056e-03 3.629493e-03
    ## [14971] 4.814159e-01          NaN 3.570035e-03 3.204838e-01 9.049206e-02
    ## [14976] 7.585702e-02 6.047141e-01 9.158470e-01          NaN 1.670619e-03
    ## [14981] 1.639535e-01 4.119948e-01 7.469542e-08 9.314346e-01          NaN
    ## [14986] 6.999143e-01 2.432683e-01 6.689505e-01 2.758901e-06 1.244443e-04
    ## [14991] 7.880284e-01 1.604101e-01 8.568059e-01 5.541430e-13          NaN
    ## [14996]          NaN 3.204838e-01          NaN 3.691786e-02 3.872729e-01
    ## [15001] 7.782501e-37 3.204838e-01          NaN 1.750437e-01 4.766361e-29
    ## [15006]          NaN 5.301578e-01          NaN 9.069958e-02          NaN
    ## [15011]          NaN 1.676445e-03 1.593207e-01          NaN 1.074501e-05
    ## [15016] 1.749474e-01 4.299720e-01 2.932667e-01          NaN 2.574129e-01
    ## [15021] 2.364040e-01 9.380303e-01 2.571744e-15 2.220385e-02 5.328492e-01
    ## [15026] 5.157689e-01 2.865020e-01 3.282464e-01 1.220923e-04 9.384255e-01
    ## [15031] 1.283811e-17 4.196966e-01 3.841422e-06          NaN 1.614536e-01
    ## [15036] 7.272117e-01 2.634271e-01 9.385838e-01 1.014758e-06 5.504435e-02
    ## [15041] 8.146116e-02 3.141258e-01 8.701232e-01 3.204838e-01 3.347907e-01
    ## [15046]          NaN 7.506597e-06 5.308010e-08 4.115744e-01 4.305968e-01
    ## [15051] 8.993379e-01 3.016856e-01 3.376950e-02          NaN 3.676306e-01
    ## [15056] 3.636445e-01 1.622887e-01 2.347885e-01          NaN 2.640687e-22
    ## [15061]          NaN 6.037241e-02          NaN 1.565655e-01          NaN
    ## [15066]          NaN 9.935528e-01          NaN          NaN          NaN
    ## [15071]          NaN          NaN          NaN          NaN          NaN
    ## [15076]          NaN          NaN          NaN          NaN 3.204838e-01
    ## [15081]          NaN          NaN          NaN          NaN          NaN
    ## [15086] 2.664592e-15          NaN          NaN          NaN 9.791599e-01
    ## [15091]          NaN          NaN          NaN          NaN          NaN
    ## [15096]          NaN          NaN          NaN          NaN          NaN
    ## [15101]          NaN          NaN          NaN          NaN          NaN
    ## [15106]          NaN          NaN          NaN          NaN          NaN
    ## [15111]          NaN          NaN          NaN          NaN          NaN
    ## [15116]          NaN          NaN          NaN          NaN          NaN
    ## [15121]          NaN          NaN          NaN          NaN          NaN
    ## [15126]          NaN          NaN          NaN          NaN          NaN
    ## [15131]          NaN          NaN          NaN          NaN          NaN
    ## [15136]          NaN          NaN          NaN          NaN          NaN
    ## [15141]          NaN          NaN 3.204838e-01          NaN          NaN
    ## [15146]          NaN          NaN          NaN          NaN          NaN
    ## [15151]          NaN 3.204838e-01          NaN          NaN          NaN
    ## [15156]          NaN          NaN          NaN          NaN          NaN
    ## [15161]          NaN          NaN          NaN          NaN          NaN
    ## [15166]          NaN          NaN          NaN          NaN          NaN
    ## [15171]          NaN          NaN 1.041700e-01          NaN          NaN
    ## [15176]          NaN 6.741102e-02          NaN          NaN 4.707310e-01
    ## [15181]          NaN 7.373575e-01          NaN 2.755066e-03          NaN
    ## [15186] 1.247489e-01          NaN 9.261560e-01          NaN          NaN
    ## [15191] 2.366411e-04 4.397730e-01 1.709836e-02          NaN          NaN
    ## [15196] 1.771308e-01 1.405849e-02 8.150797e-01          NaN          NaN
    ## [15201] 1.984894e-01 1.850520e-03 2.545162e-02          NaN 9.974059e-01
    ## [15206]          NaN          NaN          NaN 1.591270e-01 4.023510e-01
    ## [15211]          NaN          NaN 3.793468e-01 3.204838e-01          NaN
    ## [15216]          NaN          NaN 2.906886e-02 5.052258e-01          NaN
    ## [15221] 1.621321e-01 4.409923e-02          NaN 9.387595e-02 6.149507e-01
    ## [15226]          NaN          NaN 3.414373e-02 2.700654e-01 4.739576e-02
    ## [15231] 1.706115e-01          NaN          NaN          NaN          NaN
    ## [15236] 4.243377e-01 2.741308e-01 4.786534e-01 3.346754e-01          NaN
    ## [15241]          NaN 9.980244e-01          NaN          NaN 2.668350e-07
    ## [15246]          NaN          NaN          NaN          NaN          NaN
    ## [15251] 8.387105e-07 4.933351e-05          NaN 1.938649e-01 5.726430e-03
    ## [15256]          NaN 1.133704e-01 1.482387e-01 3.386144e-03 1.569837e-01
    ## [15261] 9.884892e-02 7.580368e-07 2.752881e-03 1.862706e-01 2.245051e-09
    ## [15266] 3.924621e-03 3.204838e-01 6.413854e-02 6.547002e-02 3.449754e-01
    ## [15271] 3.621256e-01 1.984751e-01 9.639552e-01 1.976375e-01          NaN
    ## [15276]          NaN          NaN          NaN          NaN 3.204838e-01
    ## [15281] 2.691924e-03          NaN 6.459117e-03          NaN          NaN
    ## [15286] 1.797231e-03 1.949588e-03 1.578568e-01 3.786704e-02 1.586034e-06
    ## [15291]          NaN 3.281276e-22 4.804350e-01          NaN 9.973974e-02
    ## [15296] 5.087176e-01          NaN 3.976701e-01 5.327993e-01          NaN
    ## [15301]          NaN 2.784244e-01 1.325654e-02 5.604740e-01          NaN
    ## [15306] 1.414409e-01 2.039816e-06 1.101004e-01 2.077394e-11 4.737166e-02
    ## [15311]          NaN          NaN 5.593685e-01          NaN 4.237618e-01
    ## [15316] 3.747407e-02          NaN 5.596811e-01 9.131558e-01 3.011820e-04
    ## [15321]          NaN          NaN 9.297278e-01 4.640671e-04          NaN
    ## [15326]          NaN          NaN 9.407257e-02 6.969377e-04          NaN
    ## [15331]          NaN 1.781945e-01 4.375694e-04 3.780268e-03 1.477603e-03
    ## [15336] 1.236314e-03 5.385808e-02 4.798133e-01 1.730948e-01 4.103972e-01
    ## [15341] 1.127389e-03 8.257393e-01 1.879678e-02 9.334189e-01 1.548556e-01
    ## [15346] 1.856679e-05          NaN 7.472519e-01          NaN 4.645169e-01
    ## [15351] 1.956741e-01 4.352045e-01 3.113676e-01          NaN          NaN
    ## [15356]          NaN 6.353953e-01 2.291031e-02          NaN 1.037703e-09
    ## [15361] 2.450193e-04          NaN          NaN 1.094811e-01 9.074027e-15
    ## [15366] 3.204838e-01 3.813431e-01 3.134111e-02          NaN 4.296143e-02
    ## [15371]          NaN 4.105380e-01          NaN 1.091138e-01 7.223233e-01
    ## [15376] 3.740506e-01 2.286893e-03          NaN 5.236858e-01 4.899982e-01
    ## [15381]          NaN          NaN 3.721831e-52 7.139011e-02          NaN
    ## [15386]          NaN 4.708367e-01 2.310921e-01 9.783731e-02 7.848727e-01
    ## [15391] 6.467678e-01 7.882235e-02 7.870763e-24 5.722978e-02 8.102762e-01
    ## [15396] 6.079249e-01 1.558447e-01 3.204838e-01          NaN          NaN
    ## [15401] 3.204838e-01 2.362385e-01 5.295068e-01 9.722232e-01          NaN
    ## [15406]          NaN 6.962632e-01 3.430076e-01 1.602305e-01 5.054465e-01
    ## [15411]          NaN 6.112678e-01 9.532610e-01          NaN          NaN
    ## [15416] 4.497684e-01          NaN          NaN          NaN 3.583954e-01
    ## [15421] 8.576358e-02 6.369636e-01 9.179179e-01          NaN 8.817785e-25
    ## [15426] 8.303526e-01          NaN          NaN          NaN 5.363620e-01
    ## [15431]          NaN          NaN 6.150894e-01          NaN          NaN
    ## [15436]          NaN          NaN 3.204838e-01          NaN 1.075307e-03
    ## [15441]          NaN 5.091492e-02 4.335298e-02          NaN          NaN
    ## [15446]          NaN 1.029143e-12          NaN          NaN          NaN
    ## [15451]          NaN          NaN          NaN 5.977545e-07          NaN
    ## [15456] 4.175978e-01          NaN 9.379307e-01          NaN          NaN
    ## [15461] 7.791467e-15          NaN 1.396264e-01 3.953052e-01          NaN
    ## [15466]          NaN          NaN 1.908545e-01          NaN 5.660259e-02
    ## [15471] 3.891516e-01          NaN 1.351581e-02 7.487118e-01 4.605947e-01
    ## [15476] 1.335302e-01 2.376228e-03 7.173411e-01 8.846957e-01 1.894455e-10
    ## [15481]          NaN          NaN          NaN 3.594507e-02 6.538913e-04
    ## [15486]          NaN          NaN 2.087223e-02          NaN 2.862720e-03
    ## [15491]          NaN 8.235413e-01 1.703841e-01 4.536222e-01 3.204838e-01
    ## [15496]          NaN 1.032770e-04 1.391265e-02          NaN 4.588837e-11
    ## [15501]          NaN          NaN 7.090144e-01 5.209171e-06          NaN
    ## [15506] 1.195308e-01 1.607384e-01 3.484723e-01 4.214813e-01          NaN
    ## [15511] 4.115556e-04          NaN          NaN 4.406034e-05          NaN
    ## [15516]          NaN 1.122687e-03          NaN 4.943375e-02          NaN
    ## [15521] 3.204838e-01          NaN 4.424778e-01 8.649284e-01 2.780491e-01
    ## [15526] 3.563717e-01 3.085914e-01 9.707338e-01 5.393430e-04 3.952631e-03
    ## [15531] 1.427163e-13 3.128888e-01 3.204838e-01 1.215659e-01          NaN
    ## [15536]          NaN          NaN 7.909891e-02 3.702814e-01 9.029183e-01
    ## [15541] 9.457262e-03 1.421525e-02 3.170624e-01          NaN 1.943017e-02
    ## [15546] 1.187399e-01 6.028866e-01          NaN 5.863195e-01 1.002300e-05
    ## [15551] 5.578172e-01 3.176697e-01 7.409651e-01          NaN 3.541002e-01
    ## [15556]          NaN 9.024107e-01 4.013610e-01 5.326733e-01          NaN
    ## [15561]          NaN          NaN 1.645460e-12          NaN          NaN
    ## [15566] 1.827494e-03          NaN 4.779265e-03          NaN          NaN
    ## [15571] 3.041106e-01          NaN 2.178434e-01          NaN 1.197067e-17
    ## [15576] 1.962145e-13 3.204838e-01          NaN          NaN 2.963249e-03
    ## [15581] 3.204838e-01 9.833419e-01 3.888473e-02          NaN          NaN
    ## [15586] 6.432213e-01          NaN 5.006111e-04          NaN 8.193021e-06
    ## [15591]          NaN 2.839793e-01          NaN          NaN 7.270255e-01
    ## [15596] 3.887104e-01 3.262217e-01          NaN          NaN          NaN
    ## [15601] 5.920692e-02 5.989192e-01          NaN 7.622831e-01          NaN
    ## [15606]          NaN 2.006856e-16 7.554642e-01 2.274360e-01          NaN
    ## [15611] 1.562197e-01 6.930759e-01          NaN          NaN 1.078009e-09
    ## [15616]          NaN 2.383708e-07          NaN 6.993939e-03 3.204838e-01
    ## [15621] 1.310604e-02 8.089101e-02 8.533995e-02          NaN          NaN
    ## [15626]          NaN 4.392552e-01          NaN 2.120024e-03 2.871425e-02
    ## [15631]          NaN 9.013710e-01          NaN 6.049311e-01 1.174315e-02
    ## [15636] 4.871994e-01 1.253459e-01          NaN          NaN 3.204838e-01
    ## [15641]          NaN 3.383025e-01 7.582696e-01          NaN 6.240747e-02
    ## [15646] 8.010186e-04 2.993438e-02          NaN 3.263445e-01 3.655187e-02
    ## [15651] 3.270907e-02 1.630574e-01 8.903331e-02 4.493401e-01 2.154437e-01
    ## [15656] 4.308745e-06 5.122385e-01          NaN 5.081700e-03          NaN
    ## [15661] 3.323188e-01          NaN 8.258019e-01          NaN          NaN
    ## [15666] 8.043910e-02 1.760925e-01 3.879099e-56 7.574595e-02          NaN
    ## [15671] 5.419494e-02 3.534675e-02          NaN 1.256959e-19 9.964203e-03
    ## [15676] 3.204838e-01          NaN          NaN 4.078325e-01 5.250675e-01
    ## [15681]          NaN          NaN          NaN          NaN          NaN
    ## [15686]          NaN 7.538336e-01 5.457753e-01          NaN          NaN
    ## [15691] 4.365151e-03 1.468362e-02          NaN          NaN 9.812964e-01
    ## [15696] 4.392603e-02          NaN 3.133169e-01          NaN          NaN
    ## [15701]          NaN 1.640541e-03          NaN 1.040659e-34 8.782480e-02
    ## [15706]          NaN 4.934058e-02          NaN 3.204838e-01          NaN
    ## [15711]          NaN 1.781695e-02 1.132191e-01 7.003934e-02 5.821302e-01
    ## [15716] 9.829182e-02 1.608180e-01          NaN          NaN          NaN
    ## [15721] 1.170723e-01          NaN 6.260282e-02 1.705518e-01          NaN
    ## [15726] 6.581252e-01          NaN 3.204838e-01          NaN          NaN
    ## [15731]          NaN 6.405368e-01 4.116566e-04 1.551006e-02 1.562680e-02
    ## [15736]          NaN          NaN 2.489835e-01 9.279215e-01          NaN
    ## [15741]          NaN          NaN 2.816946e-07          NaN 3.204838e-01
    ## [15746]          NaN          NaN          NaN          NaN 6.468299e-01
    ## [15751] 7.274277e-01 3.204838e-01          NaN          NaN 2.571983e-01
    ## [15756]          NaN 1.443201e-01          NaN 6.543053e-01          NaN
    ## [15761] 9.135285e-01          NaN 1.567805e-01          NaN          NaN
    ## [15766] 3.593055e-01          NaN          NaN 8.243829e-01 7.913047e-01
    ## [15771] 3.178949e-02 1.087391e-02 4.536778e-01          NaN 9.249522e-04
    ## [15776] 3.088253e-02 3.427432e-02          NaN          NaN 1.540623e-02
    ## [15781]          NaN          NaN          NaN          NaN          NaN
    ## [15786] 2.361887e-01 1.727700e-03 3.392588e-01          NaN 3.204838e-01
    ## [15791]          NaN 3.204838e-01          NaN 4.090871e-01 1.112653e-01
    ## [15796]          NaN 5.133660e-01          NaN 1.151457e-02 7.694028e-01
    ## [15801] 2.736670e-01 9.921426e-01          NaN 2.146170e-01 9.502811e-01
    ## [15806]          NaN 9.207473e-06          NaN          NaN 2.584463e-01
    ## [15811] 4.286512e-12          NaN 9.936149e-02          NaN 1.265793e-04
    ## [15816] 1.997993e-12 5.515964e-02          NaN          NaN 2.980582e-01
    ## [15821]          NaN          NaN 5.227062e-01 4.658738e-02          NaN
    ## [15826]          NaN 9.183717e-01          NaN          NaN          NaN
    ## [15831]          NaN 7.290240e-01 5.708740e-01          NaN 3.883897e-01
    ## [15836]          NaN          NaN          NaN 8.793020e-02          NaN
    ## [15841]          NaN          NaN 4.365446e-02 5.813094e-01          NaN
    ## [15846]          NaN 2.282362e-03          NaN 1.948150e-03 3.460481e-01
    ## [15851] 1.082238e-03          NaN 8.097677e-01          NaN 1.614687e-01
    ## [15856]          NaN          NaN          NaN          NaN 1.931662e-01
    ## [15861] 3.705914e-01 2.002848e-01 6.574231e-01 1.700882e-02 2.752000e-02
    ## [15866]          NaN          NaN 9.202838e-02 7.776558e-01 2.121495e-08
    ## [15871] 1.241825e-01 1.407214e-01          NaN 5.009304e-01          NaN
    ## [15876] 1.423831e-02 1.663517e-01 9.847860e-01 7.960196e-01 8.728858e-01
    ## [15881]          NaN 6.650616e-02 4.013424e-01 3.204838e-01 4.016477e-02
    ## [15886] 1.080726e-01          NaN          NaN 1.074543e-02          NaN
    ## [15891] 1.494392e-03          NaN          NaN          NaN 6.423956e-02
    ## [15896] 6.078817e-02          NaN 1.248501e-03 3.474057e-12          NaN
    ## [15901]          NaN          NaN 1.525667e-03          NaN          NaN
    ## [15906]          NaN 1.043959e-07          NaN          NaN          NaN
    ## [15911] 6.786659e-01 1.928756e-02          NaN 5.713078e-01 4.103094e-01
    ## [15916]          NaN 9.780088e-02          NaN 4.974803e-01 3.079063e-05
    ## [15921]          NaN 6.462883e-01 3.204838e-01          NaN 3.753861e-01
    ## [15926]          NaN 3.774940e-01 7.501831e-02          NaN 1.186186e-02
    ## [15931] 5.427990e-02 7.527789e-01 1.008750e-06          NaN 6.927547e-01
    ## [15936] 5.623745e-31 2.715216e-04 3.204838e-01          NaN 1.101551e-01
    ## [15941] 3.552692e-04 9.538280e-01 2.219945e-09          NaN 2.595626e-04
    ## [15946]          NaN 1.049637e-01          NaN 3.204838e-01 2.839132e-14
    ## [15951] 3.204838e-01          NaN 2.069255e-02          NaN 3.369712e-01
    ## [15956] 1.559089e-01 9.245962e-01 2.589913e-01 3.203708e-08 7.014853e-02
    ## [15961] 8.215701e-04 8.568442e-02          NaN          NaN 1.709877e-01
    ## [15966]          NaN 9.647069e-04 7.246500e-01          NaN 1.417742e-01
    ## [15971]          NaN 2.485402e-04 4.082586e-13 2.429531e-03 1.229905e-02
    ## [15976] 5.874443e-10 4.639827e-01          NaN 7.997829e-01 2.407145e-06
    ## [15981] 1.743561e-04 2.429128e-02          NaN          NaN          NaN
    ## [15986] 3.988188e-03 7.993074e-02          NaN 5.777418e-06 3.895659e-01
    ## [15991] 3.204838e-01          NaN          NaN          NaN          NaN
    ## [15996] 3.204838e-01          NaN          NaN 8.156676e-01          NaN
    ## [16001] 5.876959e-01 7.429208e-01 1.524294e-07 6.933004e-01 8.099083e-01
    ## [16006] 1.115530e-01 3.480835e-01 4.468705e-02 1.699086e-04 3.509321e-04
    ## [16011] 1.633444e-01 2.450847e-03          NaN          NaN 3.204838e-01
    ## [16016] 4.508904e-08 5.803048e-02          NaN 1.092900e-01 1.825693e-09
    ## [16021]          NaN          NaN 4.850832e-01          NaN 3.696710e-02
    ## [16026]          NaN          NaN          NaN          NaN          NaN
    ## [16031] 3.401799e-03 2.306501e-01          NaN 2.084985e-03          NaN
    ## [16036]          NaN          NaN          NaN 7.661799e-01 9.467171e-01
    ## [16041] 1.594985e-01 4.130264e-01 1.627006e-02 5.492786e-07          NaN
    ## [16046]          NaN 7.002344e-18          NaN 8.319358e-01          NaN
    ## [16051]          NaN 9.780088e-02 3.072176e-01 5.766339e-01 3.204838e-01
    ## [16056]          NaN          NaN 7.543543e-01 3.909544e-01 8.802140e-03
    ## [16061]          NaN 3.619572e-01 9.921426e-01          NaN          NaN
    ## [16066] 1.621778e-01 6.509445e-01          NaN 8.267461e-01          NaN
    ## [16071] 5.521185e-03          NaN 1.265335e-02          NaN          NaN
    ## [16076] 5.405341e-01 9.949090e-01          NaN 9.619554e-01 1.123101e-05
    ## [16081] 5.181095e-02 4.163474e-02 4.893966e-01          NaN 3.175362e-01
    ## [16086]          NaN 3.316102e-01 9.258699e-03          NaN 4.247045e-02
    ## [16091] 6.210956e-01          NaN 8.001329e-01          NaN          NaN
    ## [16096]          NaN 1.588314e-01 3.094844e-03 2.203284e-05 8.390249e-01
    ## [16101] 1.261037e-13 2.897375e-02 1.694556e-01 3.204838e-01 8.847134e-03
    ## [16106]          NaN          NaN 2.058674e-01 8.382879e-04 9.160828e-01
    ## [16111]          NaN 6.143966e-01          NaN          NaN 2.510005e-01
    ## [16116] 8.036279e-01 3.204838e-01 2.943059e-02          NaN          NaN
    ## [16121] 3.103844e-02 9.251633e-01 1.148087e-11 9.510539e-01          NaN
    ## [16126] 3.291933e-02 9.085935e-03 4.959569e-04          NaN          NaN
    ## [16131]          NaN 9.496565e-02 2.505490e-01 7.205555e-02 7.664580e-01
    ## [16136] 3.310805e-01          NaN          NaN          NaN 5.401565e-03
    ## [16141] 5.507384e-05 9.445046e-01          NaN          NaN 8.858044e-02
    ## [16146] 6.244921e-14 9.029135e-01          NaN          NaN          NaN
    ## [16151] 5.578519e-01 1.219678e-01 4.757032e-03          NaN          NaN
    ## [16156]          NaN 9.628660e-02 1.986798e-03 7.849893e-10 1.302999e-01
    ## [16161] 8.489112e-01          NaN 4.632569e-14          NaN          NaN
    ## [16166]          NaN 3.204838e-01 3.600331e-01 3.204838e-01          NaN
    ## [16171]          NaN 7.405532e-02 3.578062e-01          NaN 3.204838e-01
    ## [16176]          NaN          NaN 2.611335e-15 3.994540e-10 9.921426e-01
    ## [16181]          NaN          NaN 9.174173e-01          NaN 2.954453e-03
    ## [16186]          NaN 2.623880e-05 1.112993e-03          NaN 1.683798e-04
    ## [16191]          NaN          NaN 3.237381e-32          NaN 3.204838e-01
    ## [16196] 3.467535e-20 5.081737e-01          NaN          NaN          NaN
    ## [16201]          NaN 8.157151e-04 6.460501e-01          NaN          NaN
    ## [16206] 3.572765e-01 2.490354e-13          NaN          NaN 5.056782e-02
    ## [16211] 9.594541e-27          NaN          NaN 8.846957e-01          NaN
    ## [16216] 9.120412e-02          NaN          NaN 4.672250e-01 2.553150e-26
    ## [16221]          NaN 4.912732e-01 6.048888e-01 1.078388e-09          NaN
    ## [16226] 2.188783e-65 1.767784e-02          NaN          NaN 6.861579e-01
    ## [16231]          NaN 1.451468e-02          NaN          NaN 1.643400e-03
    ## [16236]          NaN 8.417535e-01 8.679665e-01          NaN          NaN
    ## [16241]          NaN 2.160509e-01 1.321630e-08 7.082377e-01 1.735818e-01
    ## [16246]          NaN 2.991626e-01          NaN          NaN 5.009758e-01
    ## [16251] 3.204838e-01 9.025406e-01 2.699814e-16 4.826038e-01          NaN
    ## [16256]          NaN          NaN 1.197580e-01 3.761361e-05 1.938767e-04
    ## [16261]          NaN          NaN          NaN 3.455260e-01          NaN
    ## [16266]          NaN 6.727787e-01 2.356985e-03          NaN          NaN
    ## [16271] 4.414220e-01 5.675549e-03          NaN 3.637936e-01 4.678454e-01
    ## [16276] 5.736972e-01 1.900178e-02          NaN 8.708459e-01 2.970294e-02
    ## [16281]          NaN 5.257328e-01 1.496123e-02          NaN 5.516926e-01
    ## [16286] 5.184194e-01 3.890978e-01 3.288047e-04 6.316141e-01 3.854681e-20
    ## [16291] 3.204838e-01          NaN 8.519418e-01          NaN 2.026491e-02
    ## [16296]          NaN 2.238129e-01 3.210066e-03 3.254910e-01          NaN
    ## [16301] 5.213370e-02          NaN 7.707878e-05 1.595536e-08          NaN
    ## [16306]          NaN 3.037435e-01 4.703617e-05 5.593725e-02 7.640947e-02
    ## [16311] 9.675438e-01          NaN 4.103594e-03 5.344078e-02 5.721649e-01
    ## [16316] 8.207897e-01 5.160973e-60 2.841062e-01 6.021129e-02          NaN
    ## [16321] 3.204838e-01          NaN 5.023231e-01          NaN 4.553422e-01
    ## [16326]          NaN 3.873004e-02 1.550654e-03          NaN 7.084002e-03
    ## [16331] 1.685813e-01          NaN          NaN 1.369713e-01          NaN
    ## [16336] 6.520000e-08          NaN 6.067313e-04 5.003272e-56 3.237140e-03
    ## [16341]          NaN          NaN 8.453256e-01          NaN          NaN
    ## [16346] 2.021401e-01 7.059345e-01 2.990886e-01 9.648626e-01          NaN
    ## [16351] 3.786339e-03 6.384277e-01 1.793117e-01 5.755847e-01 1.641477e-08
    ## [16356] 8.497280e-02          NaN 1.165849e-01 5.727378e-01 5.389882e-01
    ## [16361] 7.642272e-03 8.287973e-01 9.544446e-02 8.722299e-01          NaN
    ## [16366]          NaN          NaN 4.566770e-02 7.462117e-01          NaN
    ## [16371] 7.714374e-01          NaN          NaN          NaN 8.409550e-10
    ## [16376] 1.477904e-05 1.314173e-02          NaN 3.118821e-01          NaN
    ## [16381] 9.788291e-01 2.889015e-01 2.233961e-01          NaN          NaN
    ## [16386] 3.204838e-01 7.730575e-01          NaN          NaN 8.095264e-02
    ## [16391] 3.646816e-04 5.850552e-01 5.487732e-01          NaN          NaN
    ## [16396] 1.508125e-02 1.103629e-01 1.609430e-02 3.812535e-02 9.086001e-02
    ## [16401] 4.705515e-05 1.937498e-01          NaN 5.227574e-08 4.726200e-02
    ## [16406] 3.599632e-01 2.895287e-02 8.412007e-05 9.919006e-01 8.786135e-05
    ## [16411] 3.281276e-22 8.220015e-01 6.258510e-04          NaN 3.078454e-03
    ## [16416]          NaN 8.479721e-01 2.049202e-01 3.316871e-01          NaN
    ## [16421]          NaN 6.128094e-01          NaN          NaN 1.436820e-01
    ## [16426]          NaN 2.745404e-01 8.043910e-02 3.201043e-01 1.053990e-03
    ## [16431] 1.673190e-01          NaN 1.341199e-03          NaN          NaN
    ## [16436]          NaN 4.282090e-02 4.667735e-12          NaN 2.194729e-01
    ## [16441]          NaN 2.781873e-01          NaN          NaN 2.176217e-01
    ## [16446] 3.204838e-01 6.322387e-01 4.747685e-04 7.074328e-01          NaN
    ## [16451] 6.592257e-01 4.092824e-01 4.589349e-05 3.204838e-01 3.258115e-01
    ## [16456]          NaN          NaN          NaN          NaN          NaN
    ## [16461] 4.619458e-01 6.387504e-06 4.362976e-01          NaN 4.286763e-01
    ## [16466] 2.353528e-02 9.984747e-01          NaN          NaN          NaN
    ## [16471] 3.204838e-01 2.118405e-02 5.577403e-01 5.532260e-04          NaN
    ## [16476] 2.685630e-04 3.658083e-04          NaN          NaN          NaN
    ## [16481] 8.046119e-03 3.910734e-01          NaN 9.147499e-01          NaN
    ## [16486]          NaN          NaN          NaN          NaN 4.547008e-01
    ## [16491]          NaN          NaN 9.258498e-01 8.333912e-01          NaN
    ## [16496] 9.210060e-02 1.523663e-01 5.812181e-01          NaN          NaN
    ## [16501] 3.204838e-01          NaN 3.204838e-01          NaN 5.515440e-01
    ## [16506] 1.478791e-03          NaN 1.407214e-01          NaN          NaN
    ## [16511]          NaN 1.985799e-01 3.204838e-01 8.797445e-01 6.834392e-01
    ## [16516] 3.970236e-01 7.581625e-01 5.542465e-02 2.863673e-01 1.313362e-09
    ## [16521]          NaN 4.906331e-01 6.206022e-02          NaN          NaN
    ## [16526]          NaN 2.710897e-03 2.391705e-01 4.901439e-02          NaN
    ## [16531] 9.289400e-17 9.336477e-01          NaN          NaN 5.137850e-02
    ## [16536]          NaN 7.242454e-01 6.609061e-01 5.736559e-02 3.646487e-02
    ## [16541]          NaN 9.134256e-01          NaN          NaN 3.204838e-01
    ## [16546] 2.296396e-03 8.774854e-01          NaN 8.476911e-01          NaN
    ## [16551] 4.604446e-01          NaN          NaN 2.825147e-12          NaN
    ## [16556] 4.263382e-01          NaN          NaN 9.289400e-17          NaN
    ## [16561] 2.195166e-01 2.808132e-01          NaN 9.066307e-01          NaN
    ## [16566] 3.307406e-03          NaN          NaN 1.266455e-02 1.844303e-01
    ## [16571] 5.524751e-01          NaN 8.111983e-01          NaN          NaN
    ## [16576]          NaN          NaN 5.514853e-43 1.835926e-01 1.277074e-05
    ## [16581] 4.202676e-01 1.739391e-01 2.986068e-05          NaN          NaN
    ## [16586] 4.700843e-01          NaN 5.419069e-01          NaN          NaN
    ## [16591]          NaN          NaN 5.685774e-01 3.777334e-01          NaN
    ## [16596] 4.590099e-01 1.064444e-01          NaN 1.132105e-02          NaN
    ## [16601] 5.008602e-03 8.486446e-01          NaN          NaN          NaN
    ## [16606]          NaN 1.080494e-08          NaN 4.463973e-03          NaN
    ## [16611] 2.549494e-04 1.686282e-01 9.331625e-02          NaN 7.245302e-04
    ## [16616]          NaN 3.204838e-01 3.553899e-01          NaN 1.264431e-01
    ## [16621]          NaN          NaN 5.031006e-03          NaN          NaN
    ## [16626]          NaN          NaN          NaN 1.560753e-01          NaN
    ## [16631] 4.768594e-01 4.524628e-01 1.891802e-01 4.865517e-01          NaN
    ## [16636] 1.975664e-03 5.801463e-01          NaN 3.875464e-04          NaN
    ## [16641]          NaN 2.222075e-01 1.318163e-02 4.090871e-01 5.681111e-03
    ## [16646] 4.884594e-01 3.527161e-01 8.935356e-04 1.861181e-13 7.825339e-01
    ## [16651] 7.589010e-04 5.217984e-01 5.991236e-01 9.867352e-02 1.110064e-15
    ## [16656] 3.097027e-02 2.030648e-21 9.577891e-01 2.920830e-01 8.831984e-01
    ## [16661] 1.619713e-02          NaN          NaN          NaN 2.056564e-01
    ## [16666] 7.782501e-37 9.806545e-31 1.307607e-02 5.684627e-08 5.813094e-01
    ## [16671] 5.637960e-01 3.672999e-01 1.170581e-03 1.035141e-01 9.017953e-02
    ## [16676] 2.457568e-01 7.068718e-02 1.500807e-07 6.145401e-01 6.785709e-31
    ## [16681] 9.198576e-01 3.698806e-01 9.773031e-05 9.817781e-01 6.138731e-04
    ## [16686] 2.064565e-01 6.187453e-01 1.476259e-01 1.058243e-02 3.091594e-01
    ## [16691] 5.836515e-04 9.031019e-02 2.439914e-01 5.603463e-19 6.328893e-01
    ## [16696] 6.946922e-01 3.711562e-01 5.094705e-03 5.213522e-01 2.401901e-03
    ## [16701] 2.110348e-02 8.931419e-01 1.035649e-17 9.516105e-03 4.473886e-01
    ## [16706] 2.553150e-26 1.219804e-01 1.216509e-01 9.836480e-01 1.816945e-02
    ## [16711] 1.196737e-03 6.127283e-02 3.069705e-02 1.056636e-01 4.785263e-01
    ## [16716] 1.821323e-02          NaN 2.245323e-02 3.084814e-07 8.989327e-13
    ## [16721]          NaN 9.495829e-01          NaN 4.936525e-01          NaN
    ## [16726]          NaN 2.694972e-04          NaN          NaN          NaN
    ## [16731] 4.549167e-01          NaN          NaN 4.526499e-01          NaN
    ## [16736]          NaN 2.462846e-03          NaN 9.830180e-01          NaN
    ## [16741] 1.221857e-03          NaN 2.745221e-02 6.174397e-06          NaN
    ## [16746]          NaN 4.163147e-02 1.001206e-04          NaN 2.017469e-01
    ## [16751] 9.841295e-01 6.542064e-01 2.981896e-02 9.350802e-01 3.835127e-02
    ## [16756] 6.451496e-02 3.204838e-01          NaN 3.971203e-01 3.204838e-01
    ## [16761]          NaN 7.638269e-01          NaN          NaN          NaN
    ## [16766]          NaN 3.466707e-04 1.002996e-01 3.588970e-01          NaN
    ## [16771]          NaN          NaN 1.059842e-01 8.073480e-01          NaN
    ## [16776]          NaN 2.246803e-01 2.673972e-01          NaN 3.204838e-01
    ## [16781] 8.208413e-02 1.100543e-01 1.103393e-02 9.887065e-04          NaN
    ## [16786] 8.132841e-01          NaN          NaN 3.924163e-07          NaN
    ## [16791] 7.217223e-01 1.437617e-03 5.040180e-01          NaN 8.944636e-01
    ## [16796] 1.484600e-01 3.914321e-03          NaN          NaN          NaN
    ## [16801]          NaN          NaN 4.989410e-01 1.471942e-01          NaN
    ## [16806] 3.984578e-01 1.932485e-02          NaN          NaN 5.396378e-01
    ## [16811]          NaN 3.204838e-01          NaN          NaN 1.772579e-02
    ## [16816] 1.145563e-02 3.204838e-01 8.760440e-02 5.546690e-01          NaN
    ## [16821] 9.241027e-01 6.364714e-01          NaN 8.400506e-01          NaN
    ## [16826] 3.314843e-01 6.042795e-01          NaN 3.263360e-01 2.866314e-01
    ## [16831]          NaN 3.204838e-01          NaN          NaN 2.419398e-01
    ## [16836] 3.330211e-01 5.733628e-02          NaN 7.664783e-04          NaN
    ## [16841]          NaN 4.526444e-01 8.500760e-01 1.272808e-02          NaN
    ## [16846]          NaN 7.340082e-01          NaN 1.087632e-03          NaN
    ## [16851]          NaN 5.413787e-06          NaN          NaN          NaN
    ## [16856] 3.204838e-01          NaN 2.294783e-63          NaN 3.098659e-01
    ## [16861] 2.325756e-04          NaN          NaN 8.689707e-01          NaN
    ## [16866] 5.191802e-05          NaN 4.298810e-05 5.398431e-01          NaN
    ## [16871] 4.206714e-01 9.659745e-01 2.231344e-01 9.986713e-08 8.595510e-02
    ## [16876]          NaN 5.613049e-01 1.797023e-01 4.504032e-02 2.637307e-02
    ## [16881] 3.517767e-01 3.742150e-02 3.204838e-01          NaN 3.204838e-01
    ## [16886]          NaN          NaN 6.361392e-05 2.752634e-01          NaN
    ## [16891]          NaN 1.148087e-11 3.130737e-14          NaN          NaN
    ## [16896] 5.087894e-01 8.176606e-01          NaN          NaN          NaN
    ## [16901]          NaN          NaN 5.367669e-09 1.234683e-01          NaN
    ## [16906]          NaN          NaN 8.668850e-01          NaN 6.861642e-01
    ## [16911] 1.645424e-02 5.333841e-03 3.281200e-03 1.504040e-01 6.061645e-02
    ## [16916] 2.715533e-01          NaN 1.707249e-01          NaN          NaN
    ## [16921]          NaN          NaN 2.021544e-01 3.566152e-01          NaN
    ## [16926]          NaN          NaN 3.140052e-02          NaN          NaN
    ## [16931]          NaN 3.022509e-01          NaN 1.403741e-03          NaN
    ## [16936] 2.817580e-04 1.051624e-09          NaN 8.966627e-02          NaN
    ## [16941]          NaN 1.443764e-03          NaN 3.204838e-01 8.773732e-12
    ## [16946]          NaN 1.817227e-02 9.305829e-01 1.911368e-01 1.433392e-03
    ## [16951]          NaN 3.079184e-01 2.286620e-08          NaN 8.804821e-01
    ## [16956]          NaN 3.319667e-62          NaN          NaN 7.339458e-01
    ## [16961] 7.113185e-01          NaN 4.728233e-01 9.155054e-01 3.311026e-01
    ## [16966] 3.204838e-01          NaN 5.177305e-03 1.184948e-02 5.129764e-16
    ## [16971] 5.486142e-01          NaN 4.376132e-16 1.190233e-01 4.717145e-01
    ## [16976] 1.101579e-02 2.211817e-01          NaN          NaN          NaN
    ## [16981]          NaN          NaN          NaN          NaN          NaN
    ## [16986] 8.203720e-02 1.042377e-01 2.620927e-24 6.943443e-01 8.288812e-02
    ## [16991] 2.515336e-08          NaN          NaN          NaN 6.910715e-03
    ## [16996]          NaN          NaN          NaN 5.583993e-01 4.131309e-02
    ## [17001] 8.622136e-01 9.145151e-02 7.425334e-01 4.507516e-01          NaN
    ## [17006]          NaN 3.910218e-01          NaN 2.375646e-01 7.660940e-07
    ## [17011]          NaN 6.454525e-05          NaN 3.749382e-01 6.533462e-02
    ## [17016]          NaN 8.712569e-01          NaN          NaN          NaN
    ## [17021] 3.204838e-01 3.637936e-01          NaN 8.941851e-01          NaN
    ## [17026] 1.604768e-02          NaN 9.313057e-01 3.356810e-01 2.439415e-06
    ## [17031] 3.546146e-06 3.998953e-01          NaN 8.735685e-01 1.575480e-01
    ## [17036] 1.327779e-01 9.648180e-01 1.354615e-01 2.487110e-01 9.720802e-01
    ## [17041] 7.360759e-02 5.122179e-11          NaN 8.608552e-04 2.520638e-03
    ## [17046] 3.163052e-02 1.409573e-01 2.669984e-01          NaN 2.410024e-10
    ## [17051] 3.633201e-01          NaN          NaN 5.143015e-01 9.541946e-01
    ## [17056]          NaN 1.641628e-01          NaN 2.681069e-01 9.731126e-02
    ## [17061]          NaN          NaN 1.498196e-05 1.420578e-01 1.685257e-01
    ## [17066]          NaN 1.478465e-01 5.460155e-01 3.204838e-01 2.302592e-01
    ## [17071] 7.038899e-01 5.693233e-02          NaN 3.846432e-01          NaN
    ## [17076] 1.761201e-01 8.130198e-02          NaN 1.315460e-01 7.245302e-04
    ## [17081]          NaN 1.222790e-03          NaN          NaN          NaN
    ## [17086]          NaN 3.204838e-01 2.298474e-01 2.297567e-01          NaN
    ## [17091]          NaN 2.974155e-01          NaN          NaN 4.205764e-01
    ## [17096] 4.439478e-03 2.125043e-13          NaN 1.320803e-01 1.507602e-16
    ## [17101] 6.501494e-02          NaN 5.106545e-01          NaN 5.441524e-01
    ## [17106] 2.153433e-01 3.204838e-01 2.439757e-02 9.512767e-02          NaN
    ## [17111] 8.411726e-01 1.642775e-03 3.160142e-01          NaN          NaN
    ## [17116] 6.833813e-01 5.290444e-01 3.317979e-03          NaN 3.283612e-01
    ## [17121]          NaN 1.656626e-02          NaN 7.314692e-01          NaN
    ## [17126] 8.892496e-01 8.971179e-01 4.865392e-03 3.363186e-01 4.370124e-72
    ## [17131] 9.004471e-01          NaN          NaN          NaN          NaN
    ## [17136]          NaN 7.224968e-01 3.204838e-01          NaN 2.700870e-01
    ## [17141]          NaN 3.001544e-01 7.041737e-02          NaN          NaN
    ## [17146] 2.949712e-01 7.060669e-01 5.497005e-01 3.204838e-01          NaN
    ## [17151] 7.388318e-02 3.204838e-01          NaN          NaN          NaN
    ## [17156] 3.678258e-01          NaN          NaN 6.792197e-04 5.424292e-01
    ## [17161] 3.204838e-01 3.948691e-03 8.990163e-02          NaN 3.280140e-02
    ## [17166]          NaN 9.417189e-01          NaN 2.103793e-05 1.798193e-03
    ## [17171] 3.424483e-05 8.389997e-01 8.335781e-01 1.151398e-01          NaN
    ## [17176] 1.272310e-03 8.924045e-25 6.472136e-02          NaN 1.534259e-01
    ## [17181]          NaN 3.134179e-08 7.188284e-01          NaN          NaN
    ## [17186] 1.353869e-25 1.745927e-01 6.559951e-20          NaN 5.878313e-01
    ## [17191] 5.947099e-08 7.135630e-02 8.089101e-02          NaN          NaN
    ## [17196] 1.069259e-01          NaN 9.528209e-01 1.895279e-01 1.206060e-02
    ## [17201] 3.204838e-01          NaN 1.374212e-01 7.154799e-05 1.052986e-01
    ## [17206]          NaN 5.141308e-08          NaN 1.365533e-03 5.970612e-11
    ## [17211]          NaN          NaN          NaN          NaN 3.690800e-01
    ## [17216]          NaN 3.125720e-01 1.890845e-01          NaN          NaN
    ## [17221] 7.781506e-07 4.609409e-03 1.420795e-01 1.443750e-02 2.703322e-04
    ## [17226] 5.798448e-02          NaN          NaN 3.204838e-01          NaN
    ## [17231]          NaN          NaN 2.533092e-02 1.839944e-04 8.714961e-01
    ## [17236] 1.196737e-03          NaN          NaN          NaN          NaN
    ## [17241] 3.890341e-01          NaN          NaN          NaN          NaN
    ## [17246]          NaN 1.042830e-01          NaN 9.804184e-01          NaN
    ## [17251]          NaN          NaN 9.387057e-01 1.081334e-03          NaN
    ## [17256]          NaN 1.028924e-03          NaN          NaN          NaN
    ## [17261] 2.181337e-02 7.005300e-01          NaN 1.950041e-02 7.831054e-14
    ## [17266]          NaN          NaN 2.211864e-01 9.920836e-01 3.832170e-01
    ## [17271]          NaN          NaN          NaN          NaN 4.964638e-01
    ## [17276] 8.563664e-01 8.670093e-02 1.167749e-16 8.998469e-03          NaN
    ## [17281]          NaN 1.374027e-01 6.635348e-01          NaN          NaN
    ## [17286] 6.030881e-02          NaN 3.811347e-01          NaN 1.095592e-04
    ## [17291]          NaN          NaN 6.717201e-01 8.685448e-01 5.150527e-23
    ## [17296]          NaN          NaN 1.491785e-01          NaN 3.450133e-07
    ## [17301]          NaN 2.527569e-15          NaN          NaN          NaN
    ## [17306]          NaN          NaN          NaN 4.629703e-01 1.727017e-59
    ## [17311] 6.163383e-03          NaN          NaN 6.241062e-06 2.675473e-01
    ## [17316] 1.987550e-01 9.807605e-02 4.896682e-01 5.056143e-21 1.555488e-01
    ## [17321]          NaN          NaN 3.895292e-01 5.669962e-01          NaN
    ## [17326] 1.954261e-01          NaN          NaN          NaN 5.069108e-03
    ## [17331]          NaN          NaN 9.403895e-01 2.543392e-02 8.779525e-02
    ## [17336]          NaN          NaN 5.445901e-01 1.319771e-01 8.270081e-07
    ## [17341]          NaN          NaN          NaN 2.289933e-01 1.362165e-02
    ## [17346]          NaN 4.689121e-01          NaN 7.884173e-01          NaN
    ## [17351]          NaN 1.289717e-01          NaN          NaN          NaN
    ## [17356] 4.706449e-59          NaN 4.391607e-02 2.699016e-04          NaN
    ## [17361]          NaN          NaN          NaN          NaN 3.958988e-01
    ## [17366]          NaN          NaN 1.550078e-09 8.238609e-03 4.683796e-05
    ## [17371] 9.618291e-04          NaN          NaN 7.301899e-01 3.523215e-03
    ## [17376] 1.018460e-01 2.766204e-09 8.017697e-01 4.092680e-01          NaN
    ## [17381] 9.214724e-01          NaN 6.494607e-01 5.275335e-01 5.865205e-01
    ## [17386] 1.244964e-01          NaN          NaN          NaN 5.825486e-03
    ## [17391] 5.021095e-05 9.370512e-01 8.881906e-01          NaN          NaN
    ## [17396]          NaN          NaN 2.201313e-69 3.204838e-01          NaN
    ## [17401] 7.560747e-01          NaN 2.434731e-02 1.320251e-01 3.671626e-34
    ## [17406]          NaN 3.834891e-01          NaN 8.514257e-02          NaN
    ## [17411] 3.861285e-09 9.697196e-01 4.595837e-01          NaN          NaN
    ## [17416]          NaN          NaN          NaN 1.759249e-04          NaN
    ## [17421]          NaN          NaN 4.705038e-07          NaN          NaN
    ## [17426]          NaN 1.301481e-01 4.592472e-01 3.204838e-01 2.055690e-01
    ## [17431]          NaN 3.769944e-07          NaN          NaN 2.045549e-01
    ## [17436] 4.348034e-02          NaN 1.142025e-01 9.648385e-02          NaN
    ## [17441] 7.330508e-01 2.487471e-18          NaN 3.321178e-03          NaN
    ## [17446] 2.868383e-02 9.862002e-03 6.176177e-01          NaN          NaN
    ## [17451] 7.430123e-01 9.054332e-07 8.975376e-03          NaN 5.628687e-02
    ## [17456]          NaN 5.271410e-02          NaN 3.204838e-01          NaN
    ## [17461]          NaN          NaN 7.056387e-01 4.857413e-28 8.169524e-02
    ## [17466] 6.251817e-01          NaN 1.935214e-01          NaN 3.989431e-05
    ## [17471]          NaN 8.831556e-01 1.727619e-01 1.188173e-01 8.713923e-01
    ## [17476] 7.634951e-01 5.962355e-01          NaN 3.790876e-08          NaN
    ## [17481] 3.130346e-01          NaN          NaN 7.529871e-01          NaN
    ## [17486]          NaN          NaN 9.616366e-01          NaN 1.439150e-01
    ## [17491]          NaN          NaN 6.750194e-02 1.713734e-01 3.196677e-07
    ## [17496]          NaN 5.803023e-01          NaN          NaN 7.723841e-01
    ## [17501] 7.371070e-01 3.204838e-01          NaN          NaN          NaN
    ## [17506] 2.295271e-01 9.350802e-01          NaN 9.452296e-01          NaN
    ## [17511]          NaN 9.109014e-01 9.886074e-02          NaN 8.288220e-01
    ## [17516]          NaN          NaN 4.486289e-02 7.206147e-03 9.032093e-03
    ## [17521] 1.684001e-01 1.699137e-01 2.987500e-01 3.909734e-01 4.822384e-01
    ## [17526] 5.098504e-02 4.624700e-01 7.425426e-01 7.055011e-02          NaN
    ## [17531] 9.700384e-01          NaN          NaN 4.327033e-01 4.473560e-01
    ## [17536] 3.621483e-03          NaN          NaN 1.589175e-01          NaN
    ## [17541] 5.420893e-01          NaN          NaN 4.934524e-02          NaN
    ## [17546]          NaN 3.754420e-02 7.617830e-01          NaN 8.614569e-01
    ## [17551] 6.942751e-01          NaN 4.880422e-01 2.417122e-04 9.338025e-01
    ## [17556]          NaN          NaN          NaN 3.111507e-12          NaN
    ## [17561]          NaN 5.339758e-01 9.260620e-01 1.764462e-02 2.867584e-01
    ## [17566] 1.270156e-02          NaN          NaN          NaN          NaN
    ## [17571]          NaN          NaN 3.204838e-01          NaN          NaN
    ## [17576] 1.552761e-01 2.267716e-01 3.204838e-01 3.568963e-01          NaN
    ## [17581] 6.212630e-01 1.935243e-68          NaN 5.231313e-01 5.643786e-04
    ## [17586] 3.204838e-01 3.852995e-01 4.070279e-01          NaN          NaN
    ## [17591]          NaN          NaN 6.472396e-29 3.204838e-01 3.404730e-02
    ## [17596] 7.376424e-01 6.954486e-01 1.578817e-01 3.685324e-01 3.636444e-01
    ## [17601] 8.491426e-54 2.562838e-02 6.398063e-04 4.445570e-01 3.204838e-01
    ## [17606]          NaN          NaN          NaN 1.418642e-04 9.006655e-01
    ## [17611] 3.786832e-01          NaN          NaN 4.045166e-04          NaN
    ## [17616] 8.419990e-01 3.673201e-04 3.137890e-01 3.204838e-01          NaN
    ## [17621] 6.185633e-01          NaN          NaN 3.204838e-01 1.185443e-01
    ## [17626]          NaN 6.838785e-01 5.374654e-20 7.395791e-01 2.586380e-01
    ## [17631] 1.560800e-01 2.111596e-01 3.204838e-01          NaN          NaN
    ## [17636]          NaN 9.193304e-01 5.068076e-02 4.558723e-76 3.204838e-01
    ## [17641] 1.091189e-15          NaN 9.674461e-01 2.619282e-01 1.322431e-02
    ## [17646] 9.405119e-03          NaN 5.665524e-01          NaN 1.535171e-01
    ## [17651]          NaN 1.246581e-01 3.156343e-01          NaN          NaN
    ## [17656] 9.220097e-01 9.241365e-60 2.391327e-01          NaN          NaN
    ## [17661] 6.822247e-01 3.204838e-01 3.204838e-01          NaN 3.204838e-01
    ## [17666] 5.999583e-01          NaN          NaN 2.302570e-01          NaN
    ## [17671] 4.388790e-01 4.490501e-01 1.333960e-02          NaN 1.569847e-01
    ## [17676] 8.066107e-03          NaN 6.215461e-01          NaN          NaN
    ## [17681] 5.191797e-01          NaN          NaN 6.852403e-01 6.647077e-01
    ## [17686] 2.418565e-01 4.142149e-01 9.802785e-01          NaN          NaN
    ## [17691] 3.416479e-01          NaN          NaN          NaN          NaN
    ## [17696]          NaN 1.713927e-01 3.036556e-03          NaN 5.876751e-04
    ## [17701]          NaN 8.222183e-01          NaN 3.204838e-01          NaN
    ## [17706]          NaN 4.109654e-01 4.634659e-01          NaN 4.611572e-01
    ## [17711]          NaN 8.837391e-02          NaN 7.648501e-02 3.253993e-01
    ## [17716]          NaN          NaN 3.852995e-01          NaN 9.199722e-01
    ## [17721]          NaN 8.495766e-01          NaN          NaN          NaN
    ## [17726]          NaN 3.167422e-02          NaN          NaN          NaN
    ## [17731]          NaN          NaN          NaN          NaN 5.767336e-01
    ## [17736]          NaN 4.107263e-05 5.176983e-01          NaN          NaN
    ## [17741] 8.185185e-02          NaN 5.767254e-01          NaN          NaN
    ## [17746] 6.456194e-01 1.620129e-01 8.406732e-01 7.964030e-03 1.561956e-01
    ## [17751]          NaN 3.204838e-01          NaN 3.815118e-05          NaN
    ## [17756]          NaN 5.280639e-01          NaN          NaN 9.441593e-01
    ## [17761] 3.204838e-01          NaN          NaN          NaN 3.204838e-01
    ## [17766] 7.342135e-01 5.045239e-01          NaN          NaN          NaN
    ## [17771] 6.764634e-01          NaN 7.900888e-02          NaN          NaN
    ## [17776]          NaN          NaN 3.204838e-01 5.264463e-15          NaN
    ## [17781]          NaN          NaN          NaN          NaN 1.565038e-07
    ## [17786] 9.059230e-01 3.204838e-01          NaN 3.204838e-01          NaN
    ## [17791]          NaN          NaN          NaN 8.848180e-02 6.307063e-01
    ## [17796] 3.734670e-01          NaN          NaN 4.200976e-02 9.794758e-34
    ## [17801] 1.854485e-01          NaN 1.141684e-02          NaN          NaN
    ## [17806] 8.634868e-01 1.606509e-05          NaN 9.022256e-01 4.469309e-01
    ## [17811]          NaN          NaN          NaN          NaN          NaN
    ## [17816] 3.359823e-01          NaN          NaN 6.702505e-01 3.204838e-01
    ## [17821]          NaN          NaN          NaN          NaN          NaN
    ## [17826]          NaN          NaN          NaN 7.013404e-01          NaN
    ## [17831]          NaN          NaN 9.899165e-02 2.688166e-02 9.140745e-01
    ## [17836] 6.275028e-01          NaN 1.129547e-01          NaN          NaN
    ## [17841]          NaN          NaN          NaN          NaN          NaN
    ## [17846]          NaN 8.681423e-01          NaN          NaN 5.945833e-01
    ## [17851]          NaN          NaN          NaN 1.705218e-01          NaN
    ## [17856] 1.340825e-04 1.026815e-01          NaN          NaN          NaN
    ## [17861] 9.540708e-02 4.480395e-19 3.204838e-01 2.235064e-01          NaN
    ## [17866]          NaN          NaN 7.578832e-01          NaN          NaN
    ## [17871]          NaN          NaN          NaN          NaN 3.804752e-02
    ## [17876]          NaN          NaN          NaN          NaN 4.258184e-01
    ## [17881] 3.828233e-02          NaN 5.712777e-01          NaN 5.398428e-01
    ## [17886]          NaN          NaN 9.925189e-01          NaN          NaN
    ## [17891]          NaN          NaN 1.253628e-19          NaN 3.689025e-01
    ## [17896] 4.001486e-08          NaN          NaN 1.539832e-04 7.764550e-03
    ## [17901] 5.694236e-02          NaN 5.005106e-02 3.975847e-03          NaN
    ## [17906]          NaN          NaN          NaN 1.560564e-01 5.204436e-01
    ## [17911]          NaN          NaN          NaN 5.572373e-01          NaN
    ## [17916]          NaN          NaN          NaN          NaN          NaN
    ## [17921] 3.204838e-01 3.204838e-01          NaN          NaN          NaN
    ## [17926] 2.496400e-01 1.752057e-01          NaN          NaN 7.941216e-01
    ## [17931]          NaN          NaN 3.204838e-01          NaN 1.684195e-01
    ## [17936]          NaN 2.281958e-01 7.577114e-01          NaN          NaN
    ## [17941]          NaN 1.178284e-02 1.001775e-01 1.663401e-01          NaN
    ## [17946]          NaN 1.727422e-02 3.204838e-01 9.132244e-01 7.035201e-01
    ## [17951]          NaN          NaN          NaN 3.724183e-01 6.505591e-01
    ## [17956] 8.754475e-01          NaN 2.719817e-02          NaN          NaN
    ## [17961] 3.057774e-02          NaN          NaN 8.088908e-01 9.927796e-01
    ## [17966]          NaN 3.204838e-01 1.978496e-01          NaN          NaN
    ## [17971]          NaN          NaN 2.818222e-01          NaN 5.590670e-01
    ## [17976]          NaN 9.571508e-02          NaN          NaN          NaN
    ## [17981]          NaN          NaN 3.204838e-01 9.575376e-05 3.057656e-01
    ## [17986]          NaN          NaN          NaN 2.101443e-01          NaN
    ## [17991]          NaN 3.204838e-01          NaN 3.204838e-01          NaN
    ## [17996] 2.499380e-08          NaN          NaN          NaN          NaN
    ## [18001] 3.204838e-01          NaN 6.549668e-02          NaN          NaN
    ## [18006] 8.997653e-01          NaN          NaN          NaN          NaN
    ## [18011]          NaN 3.204838e-01 3.842170e-01          NaN          NaN
    ## [18016]          NaN 3.204838e-01          NaN 5.774808e-07          NaN
    ## [18021]          NaN 8.125871e-01          NaN          NaN 2.967793e-01
    ## [18026]          NaN          NaN 5.171431e-01          NaN          NaN
    ## [18031]          NaN          NaN 3.776822e-01          NaN 2.136399e-01
    ## [18036] 1.129049e-01 2.317527e-03 5.053509e-03 3.204838e-01          NaN
    ## [18041]          NaN          NaN 3.501410e-02          NaN          NaN
    ## [18046] 4.128150e-27          NaN 3.204838e-01          NaN 9.945937e-01
    ## [18051]          NaN          NaN 2.767141e-04 2.739748e-01          NaN
    ## [18056]          NaN 4.745032e-01 9.790647e-01          NaN 6.506530e-21
    ## [18061]          NaN          NaN 5.926761e-05 3.253823e-03          NaN
    ## [18066]          NaN          NaN          NaN          NaN 3.680662e-01
    ## [18071]          NaN          NaN          NaN          NaN          NaN
    ## [18076] 6.266149e-01 8.161125e-01          NaN          NaN 5.834086e-01
    ## [18081] 4.920794e-56          NaN          NaN 5.133639e-01          NaN
    ## [18086] 1.721680e-01 2.139296e-01          NaN          NaN          NaN
    ## [18091] 1.054708e-04 7.498921e-02          NaN          NaN          NaN
    ## [18096] 1.583781e-01 1.591750e-04 5.924637e-01 4.404285e-05          NaN
    ## [18101] 4.917865e-10          NaN          NaN 7.090923e-04 3.403108e-06
    ## [18106]          NaN          NaN          NaN 1.944583e-02          NaN
    ## [18111]          NaN 5.690112e-01          NaN 3.204838e-01          NaN
    ## [18116]          NaN 8.406030e-01 1.923491e-01          NaN 3.579755e-02
    ## [18121]          NaN          NaN 8.951987e-01          NaN          NaN
    ## [18126]          NaN 1.689030e-01          NaN          NaN 4.669244e-03
    ## [18131] 1.927526e-03 1.954477e-05          NaN 3.204838e-01 7.450867e-01
    ## [18136]          NaN 5.028819e-04 3.773726e-01 1.561702e-01 3.000220e-01
    ## [18141]          NaN          NaN          NaN 3.204838e-01          NaN
    ## [18146] 9.447648e-01 2.183628e-01          NaN 4.183491e-01          NaN
    ## [18151] 7.258885e-01          NaN          NaN          NaN 9.425724e-01
    ## [18156] 1.285460e-01          NaN          NaN          NaN          NaN
    ## [18161]          NaN 2.301659e-01          NaN          NaN 8.110888e-04
    ## [18166] 5.915463e-02          NaN 1.277392e-07          NaN 6.271118e-04
    ## [18171]          NaN          NaN          NaN          NaN          NaN
    ## [18176]          NaN          NaN 3.204838e-01          NaN 5.003926e-03
    ## [18181]          NaN          NaN 3.119763e-02 1.573691e-01          NaN
    ## [18186] 3.396267e-01          NaN          NaN 4.421694e-01 2.729980e-02
    ## [18191]          NaN          NaN 3.204838e-01 3.204838e-01 4.948956e-02
    ## [18196]          NaN 3.110801e-10 1.687775e-07          NaN          NaN
    ## [18201] 6.484872e-01 5.428898e-01          NaN          NaN          NaN
    ## [18206] 7.346638e-02 8.492621e-02          NaN 1.523875e-01          NaN
    ## [18211]          NaN          NaN 1.680949e-01 9.025207e-01 1.526980e-08
    ## [18216] 7.524397e-01          NaN          NaN 1.937821e-09          NaN
    ## [18221] 3.832937e-01          NaN          NaN 8.046161e-06          NaN
    ## [18226]          NaN 2.669271e-01 7.213369e-10 6.138984e-01          NaN
    ## [18231]          NaN          NaN          NaN 8.973215e-01          NaN
    ## [18236] 3.890178e-04 2.600769e-02 5.099701e-03          NaN          NaN
    ## [18241]          NaN 4.427146e-10 1.561952e-01          NaN 5.448018e-02
    ## [18246] 7.153194e-01          NaN 8.945962e-02 5.574697e-01          NaN
    ## [18251] 5.667060e-01          NaN          NaN          NaN 3.204838e-01
    ## [18256] 2.868684e-01 8.761683e-05          NaN 6.161448e-02          NaN
    ## [18261] 3.991827e-01          NaN          NaN          NaN 3.567851e-03
    ## [18266]          NaN          NaN 1.679853e-01          NaN 5.612415e-01
    ## [18271] 3.407453e-01          NaN          NaN 4.928413e-07          NaN
    ## [18276]          NaN          NaN 1.288848e-03          NaN          NaN
    ## [18281]          NaN 8.655218e-05          NaN          NaN          NaN
    ## [18286] 9.537330e-01          NaN 7.594816e-01          NaN 1.783028e-01
    ## [18291]          NaN          NaN 9.741027e-01          NaN          NaN
    ## [18296]          NaN          NaN 7.757121e-01          NaN 9.192439e-01
    ## [18301]          NaN          NaN          NaN          NaN          NaN
    ## [18306] 1.704674e-01 7.648664e-01 3.204838e-01 3.204838e-01          NaN
    ## [18311]          NaN          NaN          NaN          NaN          NaN
    ## [18316] 1.250374e-02 2.212973e-02 2.549530e-02          NaN          NaN
    ## [18321] 9.550205e-01 7.928488e-01          NaN 2.641071e-01          NaN
    ## [18326] 1.843742e-04          NaN          NaN          NaN          NaN
    ## [18331]          NaN          NaN          NaN          NaN 1.651154e-01
    ## [18336] 3.804085e-02          NaN          NaN          NaN          NaN
    ## [18341] 3.204838e-01          NaN 3.723829e-02          NaN 3.204838e-01
    ## [18346] 4.395932e-04          NaN 2.582770e-01 5.597096e-06          NaN
    ## [18351]          NaN 5.133330e-01 5.525447e-01          NaN          NaN
    ## [18356]          NaN 6.268917e-08          NaN 4.073083e-01 5.912253e-01
    ## [18361] 9.216254e-01          NaN          NaN 5.698584e-02 1.568198e-01
    ## [18366]          NaN 1.040694e-14          NaN          NaN          NaN
    ## [18371]          NaN          NaN          NaN          NaN          NaN
    ## [18376]          NaN 8.995252e-02          NaN 4.419509e-01 2.749040e-02
    ## [18381]          NaN 3.204838e-01          NaN          NaN          NaN
    ## [18386]          NaN          NaN 9.996101e-01          NaN          NaN
    ## [18391] 9.725920e-01          NaN 3.794796e-01          NaN          NaN
    ## [18396]          NaN          NaN 3.637166e-03          NaN 4.937767e-01
    ## [18401]          NaN 1.565649e-01          NaN          NaN 1.123370e-02
    ## [18406]          NaN          NaN          NaN 1.649390e-01          NaN
    ## [18411]          NaN 1.683382e-01          NaN 5.412696e-01          NaN
    ## [18416]          NaN 3.193021e-02 5.749340e-06 3.204838e-01          NaN
    ## [18421] 2.066847e-01          NaN          NaN          NaN 1.282456e-02
    ## [18426] 3.084300e-01          NaN          NaN 3.571313e-01 3.504408e-03
    ## [18431] 4.149761e-01          NaN 1.517730e-01          NaN          NaN
    ## [18436]          NaN 3.438466e-02 5.369129e-01 9.884207e-01          NaN
    ## [18441]          NaN 5.612598e-01          NaN          NaN          NaN
    ## [18446]          NaN 9.905209e-01 3.408416e-01          NaN 7.711363e-01
    ## [18451]          NaN          NaN          NaN 7.518019e-01 3.742314e-04
    ## [18456]          NaN          NaN 4.336180e-02 3.204838e-01 7.068473e-01
    ## [18461] 8.971453e-04          NaN          NaN          NaN          NaN
    ## [18466] 3.204838e-01          NaN 1.277455e-05          NaN          NaN
    ## [18471] 8.327692e-02 2.598492e-01          NaN          NaN          NaN
    ## [18476]          NaN 1.923423e-02          NaN          NaN          NaN
    ## [18481]          NaN          NaN 3.882277e-03          NaN          NaN
    ## [18486] 3.204838e-01 2.471178e-01          NaN 1.491788e-01          NaN
    ## [18491]          NaN          NaN 8.799339e-01          NaN          NaN
    ## [18496] 1.818651e-01          NaN 8.690836e-01 1.269905e-13          NaN
    ## [18501]          NaN          NaN 3.204838e-01 3.204838e-01          NaN
    ## [18506]          NaN          NaN 3.821789e-01          NaN 8.120677e-01
    ## [18511]          NaN          NaN 1.738238e-05          NaN 3.204838e-01
    ## [18516] 3.204838e-01          NaN          NaN          NaN 3.204838e-01
    ## [18521]          NaN 3.064643e-01          NaN          NaN          NaN
    ## [18526] 1.033917e-02          NaN          NaN          NaN          NaN
    ## [18531]          NaN          NaN 2.543812e-03 3.120793e-02          NaN
    ## [18536] 8.993026e-04 1.373740e-16 1.405618e-01          NaN          NaN
    ## [18541]          NaN 7.759384e-01 3.204838e-01          NaN          NaN
    ## [18546]          NaN 4.809833e-01          NaN          NaN 1.052986e-01
    ## [18551] 5.519518e-01 3.754752e-01 2.069552e-01 3.296045e-02          NaN
    ## [18556] 4.217442e-01 7.241830e-02          NaN 5.411110e-01          NaN
    ## [18561] 2.091463e-16          NaN 9.824320e-01 7.839100e-04 3.204838e-01
    ## [18566]          NaN 1.560211e-01          NaN          NaN          NaN
    ## [18571] 9.629987e-01          NaN          NaN          NaN 8.133375e-02
    ## [18576]          NaN 3.714002e-29          NaN          NaN 1.404823e-01
    ## [18581] 3.204838e-01          NaN          NaN 3.118276e-01          NaN
    ## [18586]          NaN 6.992672e-02 7.929318e-02          NaN          NaN
    ## [18591]          NaN 3.204838e-01          NaN          NaN 7.360457e-01
    ## [18596] 1.509334e-01          NaN 3.572099e-06          NaN          NaN
    ## [18601]          NaN          NaN          NaN 1.281749e-02          NaN
    ## [18606] 1.163997e-15          NaN 8.481398e-01          NaN          NaN
    ## [18611] 4.990067e-01 3.204838e-01          NaN          NaN          NaN
    ## [18616] 5.430274e-02 5.999240e-01          NaN          NaN          NaN
    ## [18621]          NaN 4.954462e-01          NaN 3.454243e-02          NaN
    ## [18626] 4.841219e-06          NaN          NaN 6.341271e-03          NaN
    ## [18631]          NaN          NaN          NaN 3.660543e-05 7.355935e-02
    ## [18636]          NaN          NaN          NaN          NaN 5.586232e-54
    ## [18641]          NaN 7.363024e-01 4.901675e-02          NaN          NaN
    ## [18646]          NaN 6.622707e-01          NaN 4.062912e-01          NaN
    ## [18651]          NaN 3.967426e-01 2.605520e-01          NaN          NaN
    ## [18656]          NaN          NaN          NaN          NaN 6.036153e-01
    ## [18661] 7.680573e-01 6.383099e-01 2.434888e-03          NaN          NaN
    ## [18666] 2.056893e-01          NaN 3.269554e-01          NaN          NaN
    ## [18671] 3.204838e-01          NaN 8.962272e-01          NaN          NaN
    ## [18676] 3.384073e-01 5.912817e-01          NaN          NaN          NaN
    ## [18681] 1.008332e-01          NaN          NaN 2.108578e-03 7.649487e-08
    ## [18686] 4.510338e-02          NaN          NaN 2.719651e-29          NaN
    ## [18691]          NaN          NaN          NaN 3.204838e-01 3.134720e-02
    ## [18696]          NaN          NaN          NaN          NaN 6.769819e-01
    ## [18701] 1.779162e-01          NaN          NaN          NaN          NaN
    ## [18706] 1.594595e-02          NaN          NaN          NaN          NaN
    ## [18711]          NaN 3.204838e-01          NaN          NaN          NaN
    ## [18716] 3.192259e-01 3.204838e-01 3.204838e-01          NaN 8.811108e-01
    ## [18721]          NaN 3.650552e-01 3.204838e-01          NaN          NaN
    ## [18726]          NaN 4.463313e-01 7.269719e-01          NaN          NaN
    ## [18731]          NaN          NaN          NaN          NaN          NaN
    ## [18736] 3.136068e-03          NaN          NaN          NaN          NaN
    ## [18741] 3.264143e-01          NaN          NaN 3.533149e-12          NaN
    ## [18746] 9.050555e-01          NaN          NaN          NaN 7.432913e-01
    ## [18751] 7.599116e-01 1.536529e-01          NaN 7.948060e-02          NaN
    ## [18756]          NaN 2.656033e-01 1.423835e-04          NaN          NaN
    ## [18761]          NaN 1.348863e-01 2.216838e-01          NaN          NaN
    ## [18766]          NaN          NaN          NaN          NaN          NaN
    ## [18771] 1.385478e-01 8.423690e-03 1.247579e-03          NaN          NaN
    ## [18776]          NaN          NaN          NaN          NaN 7.460849e-03
    ## [18781]          NaN          NaN          NaN          NaN 3.204838e-01
    ## [18786]          NaN          NaN          NaN          NaN 9.306306e-02
    ## [18791]          NaN 7.459290e-02          NaN          NaN          NaN
    ## [18796]          NaN 8.389986e-01          NaN 1.219194e-02          NaN
    ## [18801]          NaN          NaN          NaN          NaN          NaN
    ## [18806]          NaN 7.017291e-01          NaN 4.581410e-02 6.279423e-01
    ## [18811]          NaN          NaN          NaN 8.958487e-01          NaN
    ## [18816] 1.020191e-01 8.507725e-01 4.283381e-01 1.920405e-01 7.979252e-02
    ## [18821] 2.039382e-01          NaN          NaN 4.838516e-02          NaN
    ## [18826] 4.349515e-01 6.226850e-01          NaN 7.360095e-01          NaN
    ## [18831]          NaN          NaN 8.852429e-01 3.204838e-01 5.450256e-01
    ## [18836] 9.111772e-01          NaN 6.927077e-01 2.522370e-03          NaN
    ## [18841] 1.469360e-02          NaN          NaN 3.639132e-01          NaN
    ## [18846] 1.560382e-01          NaN 2.012594e-01          NaN 2.685586e-07
    ## [18851]          NaN          NaN          NaN          NaN          NaN
    ## [18856] 4.733269e-02          NaN 5.880605e-01          NaN          NaN
    ## [18861]          NaN          NaN          NaN          NaN          NaN
    ## [18866] 4.567888e-01          NaN          NaN          NaN          NaN
    ## [18871] 4.380094e-02 5.984850e-01          NaN 2.831527e-02          NaN
    ## [18876]          NaN          NaN 2.008940e-01          NaN          NaN
    ## [18881]          NaN 3.204838e-01          NaN 9.003135e-01          NaN
    ## [18886]          NaN          NaN          NaN          NaN 3.204838e-01
    ## [18891]          NaN 7.158362e-01          NaN          NaN 9.663652e-01
    ## [18896] 2.573469e-02          NaN 4.422553e-01          NaN 3.698407e-05
    ## [18901] 4.853160e-01          NaN          NaN          NaN          NaN
    ## [18906] 6.397073e-01          NaN          NaN          NaN 6.870381e-01
    ## [18911]          NaN          NaN          NaN 3.204838e-01 2.153531e-06
    ## [18916]          NaN          NaN          NaN          NaN          NaN
    ## [18921] 4.280163e-01 7.616606e-01          NaN          NaN 4.476800e-01
    ## [18926]          NaN          NaN 3.204838e-01 3.705350e-03          NaN
    ## [18931] 3.748429e-01          NaN 3.133051e-02 1.618026e-07 9.487801e-01
    ## [18936] 8.351548e-01 3.181217e-01          NaN          NaN          NaN
    ## [18941] 8.489443e-01          NaN          NaN 3.204838e-01          NaN
    ## [18946]          NaN          NaN          NaN          NaN 3.204838e-01
    ## [18951]          NaN 3.377839e-06          NaN          NaN          NaN
    ## [18956] 1.324390e-01 6.681955e-01 2.791753e-01 5.317832e-01 1.017480e-04
    ## [18961]          NaN 2.448101e-01 3.204838e-01          NaN          NaN
    ## [18966]          NaN 8.821074e-01 1.006405e-01          NaN 1.577563e-01
    ## [18971]          NaN 3.002439e-07          NaN          NaN 5.105474e-01
    ## [18976] 4.975917e-01          NaN 3.204838e-01 3.204838e-01 1.187662e-01
    ## [18981]          NaN          NaN          NaN 1.673892e-01 5.004597e-01
    ## [18986] 9.081026e-01          NaN          NaN          NaN          NaN
    ## [18991] 3.204838e-01          NaN          NaN 6.175162e-01 1.109886e-02
    ## [18996]          NaN 3.685925e-02 1.560495e-01          NaN          NaN
    ## [19001]          NaN          NaN          NaN 6.979132e-01          NaN
    ## [19006] 8.737559e-02          NaN 3.204838e-01 2.193587e-02 5.613328e-01
    ## [19011]          NaN          NaN          NaN          NaN 5.541667e-01
    ## [19016]          NaN          NaN 8.652613e-05          NaN 1.680658e-01
    ## [19021] 8.294152e-01          NaN          NaN          NaN          NaN
    ## [19026] 2.441876e-02          NaN 7.235423e-01 3.271850e-01 7.380392e-01
    ## [19031]          NaN          NaN 4.508725e-03          NaN 6.827980e-01
    ## [19036]          NaN          NaN          NaN          NaN 1.231105e-01
    ## [19041] 3.931800e-01          NaN 3.204838e-01          NaN          NaN
    ## [19046]          NaN          NaN          NaN          NaN 1.336411e-01
    ## [19051] 5.542697e-02 1.801169e-01 2.894135e-01 8.059656e-05 9.808781e-01
    ## [19056]          NaN 4.189070e-01          NaN          NaN          NaN
    ## [19061]          NaN 3.204838e-01 2.130159e-17 9.774915e-01          NaN
    ## [19066]          NaN 1.053238e-01          NaN 4.198777e-01          NaN
    ## [19071]          NaN 9.141130e-07          NaN 4.230474e-01 5.203833e-01
    ## [19076] 1.511909e-04          NaN          NaN 9.225764e-01 8.676220e-16
    ## [19081]          NaN          NaN          NaN          NaN          NaN
    ## [19086]          NaN          NaN          NaN          NaN 6.042205e-01
    ## [19091] 8.768731e-01          NaN 3.204838e-01 6.759589e-01          NaN
    ## [19096]          NaN          NaN 6.747589e-01 1.946266e-04 9.862801e-02
    ## [19101] 1.560551e-01          NaN          NaN 4.208547e-05 8.045035e-03
    ## [19106]          NaN          NaN 8.034712e-01 5.612799e-01          NaN
    ## [19111]          NaN          NaN          NaN          NaN 8.741798e-01
    ## [19116]          NaN          NaN          NaN          NaN          NaN
    ## [19121] 3.203985e-01          NaN 1.196760e-03          NaN          NaN
    ## [19126] 2.777399e-01          NaN          NaN          NaN 5.822501e-01
    ## [19131] 3.843303e-04          NaN          NaN          NaN          NaN
    ## [19136] 2.630516e-01 2.435988e-01          NaN          NaN          NaN
    ## [19141]          NaN 5.628012e-01 3.204838e-01 3.204838e-01          NaN
    ## [19146] 7.748348e-01          NaN          NaN          NaN          NaN
    ## [19151]          NaN 3.385867e-02          NaN 1.908762e-01          NaN
    ## [19156]          NaN 3.204838e-01 3.094364e-01          NaN 1.199379e-01
    ## [19161] 3.204838e-01 3.684657e-01          NaN          NaN 3.204838e-01
    ## [19166]          NaN 3.204838e-01 1.086453e-01 4.465728e-01 1.800305e-02
    ## [19171]          NaN          NaN 5.973162e-01          NaN 9.453201e-01
    ## [19176] 5.600487e-04          NaN          NaN          NaN 1.464777e-03
    ## [19181]          NaN          NaN          NaN 5.034216e-01 3.232640e-01
    ## [19186]          NaN          NaN          NaN          NaN 1.243801e-02
    ## [19191] 6.344134e-01 4.286580e-01          NaN 3.204838e-01 7.670388e-02
    ## [19196] 6.058711e-01 1.878870e-02 3.935523e-01          NaN          NaN
    ## [19201]          NaN          NaN 1.120300e-05          NaN          NaN
    ## [19206] 3.204838e-01          NaN          NaN          NaN 9.342919e-01
    ## [19211] 1.595650e-01          NaN 7.396899e-01 3.204838e-01          NaN
    ## [19216]          NaN          NaN 1.605970e-01 3.204838e-01 3.347358e-01
    ## [19221]          NaN 5.165049e-01 1.577826e-02          NaN          NaN
    ## [19226]          NaN 3.204838e-01 3.204838e-01 3.305139e-01 2.937829e-01
    ## [19231]          NaN 8.990570e-01 4.492186e-02          NaN          NaN
    ## [19236]          NaN 3.636738e-01 1.601202e-01          NaN 1.562204e-01
    ## [19241]          NaN 1.585354e-01 6.607618e-01          NaN 1.355914e-02
    ## [19246]          NaN          NaN          NaN 3.942371e-01 1.994726e-02
    ## [19251]          NaN          NaN          NaN          NaN          NaN
    ## [19256]          NaN 2.305499e-01          NaN          NaN 2.466155e-16
    ## [19261] 3.204838e-01          NaN          NaN          NaN          NaN
    ## [19266] 2.294904e-01 3.204838e-01 1.559913e-01          NaN 3.204838e-01
    ## [19271]          NaN 4.283567e-07          NaN 1.211029e-02          NaN
    ## [19276]          NaN 3.101174e-01 3.485726e-02          NaN 4.396458e-01
    ## [19281] 1.449181e-01 9.862707e-01 1.661790e-01 2.588338e-01          NaN
    ## [19286] 4.973065e-15 7.564160e-01 4.117151e-01          NaN          NaN
    ## [19291] 3.204838e-01          NaN          NaN 3.204838e-01 6.947074e-01
    ## [19296]          NaN          NaN          NaN          NaN 9.510120e-01
    ## [19301]          NaN 8.869890e-01 8.773332e-04 3.204838e-01 9.841446e-01
    ## [19306]          NaN          NaN 9.448429e-01          NaN 8.331218e-20
    ## [19311] 2.244658e-02          NaN          NaN          NaN          NaN
    ## [19316] 3.204838e-01          NaN          NaN          NaN 1.690174e-01
    ## [19321] 3.204838e-01          NaN          NaN          NaN          NaN
    ## [19326] 4.061993e-01 8.012360e-01 3.140180e-01 4.688176e-08          NaN
    ## [19331]          NaN 5.390150e-01          NaN          NaN          NaN
    ## [19336]          NaN          NaN 7.300212e-01          NaN 1.569301e-01
    ## [19341]          NaN          NaN 6.328398e-01          NaN          NaN
    ## [19346] 3.978869e-01 4.197217e-01          NaN 1.422535e-01          NaN
    ## [19351]          NaN 8.689112e-01 2.881279e-01 9.699863e-07          NaN
    ## [19356] 9.985814e-01 3.204838e-01          NaN 3.204838e-01          NaN
    ## [19361]          NaN 2.986903e-01 3.204838e-01          NaN 6.887495e-01
    ## [19366]          NaN 1.314326e-01          NaN          NaN 2.690486e-01
    ## [19371] 9.023677e-01          NaN          NaN 4.109564e-01          NaN
    ## [19376]          NaN 1.922841e-05          NaN 7.233940e-01 9.393697e-02
    ## [19381] 4.056047e-01 6.431531e-01 1.274335e-01          NaN 3.204838e-01
    ## [19386]          NaN          NaN          NaN          NaN 1.245551e-27
    ## [19391]          NaN          NaN          NaN          NaN          NaN
    ## [19396]          NaN 2.320427e-02          NaN          NaN          NaN
    ## [19401] 3.204838e-01          NaN 3.654458e-01          NaN 1.615557e-01
    ## [19406] 7.922916e-01          NaN 8.562008e-01          NaN 2.166656e-01
    ## [19411] 3.204838e-01          NaN 9.368727e-01 8.156945e-01 3.204838e-01
    ## [19416]          NaN 1.190841e-01 8.664498e-01 1.027303e-01 8.280987e-01
    ## [19421]          NaN          NaN          NaN          NaN 5.249931e-03
    ## [19426] 6.561858e-01          NaN          NaN          NaN 8.955163e-01
    ## [19431]          NaN          NaN          NaN          NaN 1.831940e-04
    ## [19436]          NaN 4.485884e-02          NaN          NaN          NaN
    ## [19441]          NaN          NaN          NaN 4.834218e-06          NaN
    ## [19446]          NaN          NaN          NaN 3.100829e-12          NaN
    ## [19451]          NaN 1.403398e-07 8.640675e-02 5.359500e-08 1.922334e-21
    ## [19456]          NaN          NaN 1.887315e-03 3.204838e-01 2.815957e-11
    ## [19461]          NaN 7.630227e-01          NaN          NaN 2.783684e-01
    ## [19466] 1.930690e-11          NaN 1.844556e-01          NaN          NaN
    ## [19471] 4.549868e-12          NaN          NaN          NaN 9.451787e-14
    ## [19476]          NaN 5.382423e-01          NaN 1.714598e-01          NaN
    ## [19481]          NaN          NaN          NaN          NaN 3.498616e-03
    ## [19486]          NaN          NaN 1.015931e-02          NaN 3.305030e-03
    ## [19491]          NaN 3.204838e-01          NaN 6.263681e-01 5.882036e-01
    ## [19496] 7.018753e-01 5.408361e-01          NaN          NaN          NaN
    ## [19501]          NaN          NaN          NaN          NaN          NaN
    ## [19506]          NaN          NaN 3.204838e-01          NaN          NaN
    ## [19511]          NaN          NaN 2.081221e-06          NaN 7.578810e-01
    ## [19516]          NaN          NaN          NaN 9.352354e-03 6.463227e-01
    ## [19521] 2.756594e-01          NaN          NaN 3.204838e-01 3.243322e-01
    ## [19526]          NaN 1.559710e-01          NaN 7.703412e-09 2.716538e-01
    ## [19531] 1.162971e-23          NaN          NaN 9.554845e-01          NaN
    ## [19536]          NaN          NaN          NaN 2.852717e-01          NaN
    ## [19541]          NaN 3.204838e-01 8.965938e-01 3.690907e-04          NaN
    ## [19546]          NaN          NaN          NaN          NaN          NaN
    ## [19551]          NaN          NaN          NaN          NaN 4.533497e-23
    ## [19556] 6.952693e-01 1.345092e-01 9.999183e-06          NaN          NaN
    ## [19561] 2.557767e-06          NaN 3.204838e-01          NaN          NaN
    ## [19566]          NaN          NaN 3.204838e-01          NaN          NaN
    ## [19571]          NaN          NaN          NaN 3.141627e-02          NaN
    ## [19576]          NaN          NaN          NaN 1.786378e-01 1.263273e-03
    ## [19581]          NaN          NaN          NaN 1.514311e-04          NaN
    ## [19586] 7.685900e-01 3.823426e-01 9.499290e-01 3.204838e-01 1.214089e-01
    ## [19591]          NaN          NaN          NaN          NaN          NaN
    ## [19596]          NaN          NaN          NaN          NaN          NaN
    ## [19601] 6.068310e-01          NaN          NaN          NaN          NaN
    ## [19606] 3.358869e-02 3.150927e-02          NaN 4.102926e-03 9.921115e-01
    ## [19611] 8.205698e-01          NaN          NaN 4.053483e-01 3.204838e-01
    ## [19616] 9.403238e-01 4.710820e-02          NaN 1.401891e-02          NaN
    ## [19621]          NaN 2.154497e-01 1.365193e-02 4.379400e-02 3.899073e-02
    ## [19626]          NaN          NaN 1.628894e-01          NaN          NaN
    ## [19631]          NaN 1.614779e-01          NaN 3.204838e-01 2.331361e-11
    ## [19636]          NaN          NaN 1.575845e-01 6.342717e-02          NaN
    ## [19641]          NaN 2.427576e-02          NaN 9.823275e-01 2.100724e-01
    ## [19646] 5.720336e-01          NaN 1.167666e-01          NaN 9.093908e-01
    ## [19651]          NaN          NaN          NaN 4.522969e-01          NaN
    ## [19656]          NaN 2.099758e-04          NaN          NaN 9.505968e-03
    ## [19661]          NaN          NaN          NaN 1.831947e-01 5.132606e-02
    ## [19666] 4.304954e-01          NaN          NaN          NaN          NaN
    ## [19671] 1.567835e-02 9.475347e-01          NaN          NaN 3.204838e-01
    ## [19676] 4.297146e-02          NaN 5.024533e-01          NaN 2.940480e-02
    ## [19681] 1.677746e-03          NaN 3.204838e-01 3.204838e-01 4.628386e-01
    ## [19686] 8.256713e-01          NaN          NaN 9.069959e-01          NaN
    ## [19691]          NaN          NaN 1.565203e-01          NaN 8.092636e-01
    ## [19696] 1.428577e-12 3.242283e-06          NaN          NaN 8.849977e-01
    ## [19701]          NaN          NaN          NaN          NaN 3.225388e-01
    ## [19706] 7.009165e-01          NaN          NaN 1.369731e-01          NaN
    ## [19711]          NaN 1.602332e-02          NaN          NaN          NaN
    ## [19716] 8.918283e-01 9.459029e-05          NaN          NaN          NaN
    ## [19721] 2.528662e-01 1.649926e-01 1.846767e-03 1.124587e-02 7.562879e-01
    ## [19726] 3.947667e-01 1.562186e-01          NaN 6.530460e-05 4.933947e-08
    ## [19731]          NaN          NaN 2.002249e-01 1.850323e-02          NaN
    ## [19736] 2.357618e-07          NaN          NaN          NaN          NaN
    ## [19741] 4.804980e-01 1.430779e-26 2.807582e-01 9.034778e-07 1.553811e-01
    ## [19746] 1.317366e-01 2.395955e-08          NaN 1.583444e-01          NaN
    ## [19751] 5.480268e-03 8.065481e-02 6.613053e-01 8.044463e-01 3.204838e-01
    ## [19756]          NaN          NaN          NaN          NaN          NaN
    ## [19761]          NaN 8.393939e-01 4.078997e-01          NaN          NaN
    ## [19766] 2.294471e-02          NaN          NaN 9.483737e-01          NaN
    ## [19771] 1.335966e-01 5.607375e-01          NaN 6.228463e-01 1.778028e-01
    ## [19776] 3.204838e-01          NaN          NaN          NaN 1.625181e-07
    ## [19781]          NaN 9.837240e-01          NaN          NaN          NaN
    ## [19786]          NaN          NaN          NaN          NaN          NaN
    ## [19791] 1.170884e-01          NaN 9.414991e-01          NaN          NaN
    ## [19796] 3.669542e-01 2.994943e-01          NaN          NaN          NaN
    ## [19801]          NaN 8.915893e-01          NaN          NaN          NaN
    ## [19806] 1.771894e-01 4.489843e-01 1.994847e-16          NaN          NaN
    ## [19811]          NaN          NaN 6.564657e-01          NaN          NaN
    ## [19816]          NaN          NaN 3.204838e-01          NaN          NaN
    ## [19821]          NaN          NaN 7.224192e-09          NaN          NaN
    ## [19826] 6.937000e-01          NaN          NaN          NaN          NaN
    ## [19831]          NaN 8.937920e-02          NaN 5.424450e-01 1.048514e-03
    ## [19836]          NaN 1.979598e-02          NaN 8.887724e-01          NaN
    ## [19841]          NaN          NaN          NaN 3.204838e-01          NaN
    ## [19846]          NaN          NaN          NaN          NaN 1.076160e-03
    ## [19851]          NaN          NaN 8.575038e-02 3.204838e-01          NaN
    ## [19856]          NaN 1.576198e-01 3.484717e-01 8.531849e-01          NaN
    ## [19861]          NaN 1.673877e-01          NaN          NaN          NaN
    ## [19866]          NaN 2.429995e-04 2.057364e-06 8.950228e-02          NaN
    ## [19871]          NaN          NaN          NaN          NaN 9.007501e-01
    ## [19876]          NaN 5.138727e-01          NaN 3.204838e-01          NaN
    ## [19881]          NaN 3.925137e-05 2.114717e-01          NaN          NaN
    ## [19886] 9.266755e-01          NaN          NaN 1.823428e-01 4.813624e-01
    ## [19891] 4.646347e-02 3.948669e-16          NaN          NaN 4.951558e-10
    ## [19896]          NaN 3.204838e-01          NaN 2.510560e-01 1.284064e-02
    ## [19901] 7.097750e-01 2.085276e-01          NaN          NaN          NaN
    ## [19906]          NaN          NaN          NaN 1.067519e-03          NaN
    ## [19911]          NaN          NaN 6.672375e-01          NaN          NaN
    ## [19916]          NaN 3.528846e-02          NaN          NaN 8.749549e-01
    ## [19921] 1.845312e-01 3.200257e-01          NaN 8.085512e-02          NaN
    ## [19926]          NaN          NaN 1.707999e-01          NaN          NaN
    ## [19931]          NaN 3.650954e-01 8.747125e-01          NaN 1.566164e-01
    ## [19936]          NaN 3.204838e-01          NaN 7.150171e-01          NaN
    ## [19941] 3.380653e-24          NaN          NaN          NaN 4.263614e-02
    ## [19946]          NaN 6.092589e-01 3.816409e-01          NaN          NaN
    ## [19951]          NaN 7.822450e-01 2.359893e-02          NaN 4.101428e-02
    ## [19956]          NaN          NaN 1.752247e-01          NaN          NaN
    ## [19961]          NaN          NaN 9.113940e-01          NaN          NaN
    ## [19966]          NaN 9.542923e-01          NaN          NaN          NaN
    ## [19971]          NaN          NaN          NaN          NaN 1.559120e-01
    ## [19976]          NaN          NaN          NaN 2.807540e-13 8.673674e-01
    ## [19981] 3.427773e-02          NaN          NaN          NaN          NaN
    ## [19986] 9.378905e-05          NaN          NaN 2.656685e-04 2.995564e-01
    ## [19991]          NaN 9.108366e-01          NaN 4.081700e-03 4.332024e-02
    ## [19996] 3.591138e-04          NaN          NaN          NaN          NaN
    ## [20001]          NaN          NaN          NaN 1.118504e-06          NaN
    ## [20006]          NaN          NaN 3.204838e-01          NaN          NaN
    ## [20011] 4.837616e-40          NaN 7.487520e-01          NaN 1.667423e-14
    ## [20016]          NaN          NaN          NaN 8.111939e-01 1.559035e-01
    ## [20021] 3.188863e-04 1.295158e-03          NaN          NaN 1.559899e-01
    ## [20026] 9.350709e-01 3.369787e-01 4.866348e-01 3.392023e-01 1.598918e-02
    ## [20031]          NaN 2.349305e-01          NaN 4.326191e-01 1.079497e-09
    ## [20036] 1.319575e-01 4.195097e-01 3.204838e-01          NaN          NaN
    ## [20041]          NaN          NaN          NaN 1.933567e-06          NaN
    ## [20046]          NaN          NaN 1.563949e-01          NaN          NaN
    ## [20051] 6.567619e-01 8.919268e-09          NaN 7.624364e-02          NaN
    ## [20056]          NaN 3.204838e-01 6.981261e-01          NaN 9.597574e-01
    ## [20061] 2.261682e-01          NaN 6.771752e-04 3.204838e-01 3.655236e-01
    ## [20066] 3.204838e-01          NaN 2.113910e-01          NaN 5.130712e-02
    ## [20071]          NaN          NaN          NaN 7.968529e-02 6.909060e-01
    ## [20076]          NaN          NaN          NaN          NaN          NaN
    ## [20081] 6.602076e-01          NaN 2.028809e-02 8.342578e-02          NaN
    ## [20086]          NaN          NaN          NaN 8.206459e-01          NaN
    ## [20091] 6.455260e-02 9.436798e-01          NaN 7.847929e-01          NaN
    ## [20096]          NaN          NaN 8.696442e-01          NaN 1.255216e-02
    ## [20101] 1.176511e-02 3.180719e-01          NaN          NaN          NaN
    ## [20106]          NaN 2.598711e-24          NaN 4.920092e-01 9.571639e-01
    ## [20111]          NaN          NaN          NaN 9.492116e-03          NaN
    ## [20116] 7.703022e-02 6.933989e-01          NaN          NaN          NaN
    ## [20121]          NaN          NaN 1.560285e-01 8.628900e-01 4.129580e-01
    ## [20126]          NaN          NaN          NaN          NaN          NaN
    ## [20131]          NaN          NaN 9.378255e-01          NaN          NaN
    ## [20136] 3.204838e-01          NaN 3.204838e-01          NaN 7.192389e-01
    ## [20141] 1.426751e-16          NaN          NaN 9.791802e-01          NaN
    ## [20146]          NaN          NaN 4.179141e-04          NaN          NaN
    ## [20151]          NaN          NaN 9.825085e-04          NaN          NaN
    ## [20156]          NaN          NaN          NaN 2.780531e-21 4.993552e-06
    ## [20161]          NaN          NaN          NaN 8.471249e-01          NaN
    ## [20166] 8.267122e-01 4.286981e-02          NaN 5.451862e-03          NaN
    ## [20171]          NaN          NaN          NaN          NaN          NaN
    ## [20176] 3.204838e-01 1.646483e-02          NaN          NaN 6.997554e-01
    ## [20181] 7.839666e-01          NaN          NaN 1.495442e-04 3.082362e-02
    ## [20186]          NaN          NaN 1.493484e-01 9.476789e-04 4.148795e-01
    ## [20191]          NaN          NaN          NaN 9.759421e-01 1.798193e-05
    ## [20196]          NaN 8.425788e-01          NaN          NaN          NaN
    ## [20201]          NaN          NaN 5.184723e-02          NaN 2.738459e-03
    ## [20206] 5.550727e-01          NaN 8.802412e-01 4.121386e-01          NaN
    ## [20211]          NaN          NaN 6.245452e-01 6.682563e-01          NaN
    ## [20216] 5.181393e-09          NaN 3.228959e-01 9.780940e-01          NaN
    ## [20221] 8.857276e-01          NaN 9.154001e-01          NaN          NaN
    ## [20226] 1.451916e-11 7.072129e-01 1.503421e-01 3.669651e-01          NaN
    ## [20231]          NaN          NaN 1.967250e-25          NaN 8.031226e-01
    ## [20236] 3.204838e-01 1.490132e-04 4.028640e-01          NaN          NaN
    ## [20241]          NaN 2.022761e-01 6.761460e-01 3.947260e-01          NaN
    ## [20246]          NaN          NaN          NaN 6.430520e-11          NaN
    ## [20251]          NaN 5.096042e-01          NaN          NaN          NaN
    ## [20256]          NaN          NaN          NaN 3.204838e-01 1.487607e-01
    ## [20261] 6.783270e-01          NaN 2.417853e-01 7.675899e-01          NaN
    ## [20266]          NaN 4.541077e-01 1.533673e-01          NaN          NaN
    ## [20271]          NaN          NaN          NaN 5.530206e-13          NaN
    ## [20276]          NaN          NaN 1.816340e-01          NaN 2.761480e-01
    ## [20281]          NaN 3.266138e-01          NaN          NaN          NaN
    ## [20286] 3.204838e-01          NaN          NaN 3.534144e-02 7.080403e-01
    ## [20291]          NaN          NaN 3.204838e-01 5.734414e-01          NaN
    ## [20296]          NaN          NaN          NaN 1.620761e-01 9.770215e-04
    ## [20301]          NaN          NaN 1.559069e-01          NaN 1.466784e-15
    ## [20306] 8.707812e-04 1.585039e-01          NaN          NaN 2.653746e-01
    ## [20311] 6.862718e-02          NaN          NaN          NaN          NaN
    ## [20316]          NaN          NaN          NaN 2.056840e-04          NaN
    ## [20321] 6.069145e-01          NaN 2.593761e-03          NaN          NaN
    ## [20326]          NaN          NaN 5.355515e-01          NaN 4.803994e-01
    ## [20331] 5.721839e-01          NaN 8.740821e-01 5.495049e-01 3.204838e-01
    ## [20336]          NaN          NaN          NaN 6.944367e-03          NaN
    ## [20341] 4.159836e-02          NaN          NaN 9.107455e-01          NaN
    ## [20346] 8.544377e-01          NaN          NaN          NaN 6.667442e-01
    ## [20351]          NaN          NaN          NaN 5.925402e-01 9.225769e-01
    ## [20356]          NaN 3.204838e-01          NaN 3.239133e-02 1.494624e-03
    ## [20361] 2.171336e-01 8.026466e-01          NaN          NaN 3.752869e-04
    ## [20366] 3.942909e-03          NaN 9.306461e-01          NaN          NaN
    ## [20371]          NaN 6.270245e-01          NaN          NaN          NaN
    ## [20376]          NaN          NaN          NaN          NaN          NaN
    ## [20381] 1.282173e-02          NaN 3.658862e-02 2.187597e-01          NaN
    ## [20386] 1.001771e-02          NaN 1.431071e-03 9.739189e-01          NaN
    ## [20391]          NaN 2.676262e-07 6.015779e-01 5.823356e-02 7.355266e-01
    ## [20396]          NaN          NaN          NaN          NaN 5.818370e-01
    ## [20401] 3.204838e-01 3.304724e-01 3.192259e-01 2.098633e-01 6.187823e-01
    ## [20406]          NaN          NaN 1.641597e-01          NaN          NaN
    ## [20411]          NaN          NaN 2.232418e-04          NaN          NaN
    ## [20416] 6.844858e-03          NaN          NaN          NaN 4.417338e-03
    ## [20421] 2.282478e-01          NaN 6.938625e-07 8.294048e-01          NaN
    ## [20426] 6.286665e-02          NaN 3.204838e-01          NaN          NaN
    ## [20431]          NaN 1.733232e-01 2.483706e-06 3.204838e-01          NaN
    ## [20436] 2.365761e-01          NaN          NaN          NaN 5.078269e-10
    ## [20441] 7.771941e-02          NaN 2.196383e-01          NaN          NaN
    ## [20446]          NaN          NaN 3.204838e-01 8.037432e-02          NaN
    ## [20451] 3.744982e-03          NaN          NaN 4.259939e-02          NaN
    ## [20456]          NaN          NaN          NaN          NaN          NaN
    ## [20461] 2.534119e-01          NaN          NaN          NaN          NaN
    ## [20466]          NaN          NaN 9.336908e-01 2.904278e-03          NaN
    ## [20471]          NaN 4.336208e-02 4.031118e-01 1.582053e-01          NaN
    ## [20476] 2.893400e-16          NaN          NaN          NaN          NaN
    ## [20481]          NaN          NaN          NaN          NaN          NaN
    ## [20486]          NaN 3.204838e-01 8.216099e-01 3.382672e-01 2.840814e-01
    ## [20491]          NaN          NaN 1.190841e-01          NaN 5.761040e-01
    ## [20496]          NaN          NaN 6.917986e-01 8.451823e-02          NaN
    ## [20501]          NaN          NaN 2.153615e-02 1.470000e-02          NaN
    ## [20506]          NaN 1.627699e-01 9.451729e-06 8.045947e-01          NaN
    ## [20511]          NaN 1.146943e-42 1.786657e-01 8.021242e-02          NaN
    ## [20516] 1.973682e-03 6.305650e-02          NaN          NaN          NaN
    ## [20521]          NaN          NaN 5.562707e-01          NaN 3.961666e-02
    ## [20526]          NaN          NaN 3.170041e-01          NaN          NaN
    ## [20531]          NaN 1.123807e-01          NaN          NaN 8.310780e-01
    ## [20536]          NaN          NaN 9.650186e-01          NaN 5.443897e-08
    ## [20541] 3.204838e-01 1.706949e-01 3.204838e-01          NaN 3.204838e-01
    ## [20546]          NaN 8.141370e-01 1.163346e-02          NaN 9.425724e-01
    ## [20551]          NaN 2.408216e-01 1.032668e-01          NaN          NaN
    ## [20556] 2.386112e-02          NaN          NaN          NaN          NaN
    ## [20561]          NaN 4.893247e-01 3.204838e-01          NaN          NaN
    ## [20566]          NaN          NaN          NaN          NaN          NaN
    ## [20571]          NaN          NaN          NaN 9.209140e-01          NaN
    ## [20576] 9.183225e-01          NaN          NaN          NaN          NaN
    ## [20581] 8.197080e-02 3.204838e-01          NaN 7.534787e-01          NaN
    ## [20586]          NaN          NaN 6.129961e-01 6.794714e-01 8.613410e-01
    ## [20591]          NaN          NaN 2.110981e-09 1.398143e-03          NaN
    ## [20596]          NaN          NaN 3.204838e-01          NaN 1.775747e-01
    ## [20601]          NaN          NaN          NaN          NaN          NaN
    ## [20606]          NaN          NaN          NaN          NaN          NaN
    ## [20611]          NaN          NaN          NaN          NaN          NaN
    ## [20616] 3.204838e-01          NaN          NaN          NaN          NaN
    ## [20621]          NaN          NaN 3.204838e-01 1.364457e-02          NaN
    ## [20626] 5.084900e-01 6.933568e-01 1.601924e-01 1.719619e-01          NaN
    ## [20631]          NaN          NaN          NaN 5.651548e-01 3.438650e-09
    ## [20636] 3.204838e-01 5.218862e-01          NaN          NaN          NaN
    ## [20641] 1.440619e-01          NaN 4.074275e-02          NaN          NaN
    ## [20646]          NaN          NaN          NaN          NaN 2.654573e-03
    ## [20651]          NaN 9.494713e-01          NaN 4.539166e-06 6.851550e-01
    ## [20656] 5.779548e-01 3.295648e-05 3.204838e-01          NaN 3.204838e-01
    ## [20661] 3.204838e-01          NaN 9.394503e-02          NaN          NaN
    ## [20666]          NaN          NaN 2.231681e-01          NaN          NaN
    ## [20671]          NaN 1.784519e-06          NaN          NaN          NaN
    ## [20676] 5.363107e-01          NaN          NaN          NaN          NaN
    ## [20681]          NaN 2.670766e-08          NaN 2.194424e-05          NaN
    ## [20686]          NaN          NaN 6.173201e-02          NaN          NaN
    ## [20691] 1.840467e-01          NaN 3.204838e-01          NaN          NaN
    ## [20696]          NaN          NaN 4.777183e-01          NaN 3.204838e-01
    ## [20701] 2.453307e-01          NaN          NaN          NaN          NaN
    ## [20706]          NaN 5.464035e-01          NaN 7.990226e-02          NaN
    ## [20711]          NaN          NaN          NaN 1.665562e-06 7.249254e-02
    ## [20716] 3.520293e-04 1.740737e-01 1.220056e-01          NaN          NaN
    ## [20721] 6.669561e-01 5.654165e-01          NaN          NaN 6.755411e-01
    ## [20726]          NaN 3.972062e-01 6.250540e-01          NaN 5.186138e-01
    ## [20731] 8.784064e-03          NaN 4.109455e-11 4.045892e-01          NaN
    ## [20736]          NaN 5.361324e-01 4.538122e-01 1.160781e-02          NaN
    ## [20741] 1.270139e-02 3.862950e-15          NaN 9.336477e-01          NaN
    ## [20746]          NaN 6.481733e-01          NaN          NaN          NaN
    ## [20751] 9.516256e-01          NaN 8.914285e-01          NaN          NaN
    ## [20756] 3.204838e-01          NaN 1.053947e-01 3.204838e-01 5.522191e-05
    ## [20761]          NaN          NaN          NaN          NaN          NaN
    ## [20766] 1.099783e-06          NaN 3.204838e-01 9.624301e-01 8.199752e-03
    ## [20771]          NaN 3.970223e-01          NaN 6.520528e-02          NaN
    ## [20776]          NaN          NaN 2.393839e-07          NaN 1.807453e-01
    ## [20781] 6.688852e-01          NaN          NaN 6.842824e-03 8.191155e-01
    ## [20786]          NaN          NaN          NaN          NaN 8.798217e-06
    ## [20791]          NaN 5.558586e-01 3.204838e-01          NaN          NaN
    ## [20796]          NaN          NaN 1.437256e-11 3.204838e-01          NaN
    ## [20801]          NaN          NaN          NaN 3.070732e-01 1.201344e-01
    ## [20806] 4.668839e-01          NaN          NaN 5.294541e-03          NaN
    ## [20811] 7.417499e-01          NaN          NaN 3.204838e-01          NaN
    ## [20816]          NaN 7.228028e-01          NaN          NaN          NaN
    ## [20821]          NaN 5.427185e-01          NaN 2.266108e-01 1.480122e-02
    ## [20826] 4.557814e-01 4.514558e-02          NaN 4.677377e-01          NaN
    ## [20831] 1.465114e-05          NaN          NaN          NaN          NaN
    ## [20836] 2.753265e-01          NaN          NaN          NaN          NaN
    ## [20841]          NaN          NaN          NaN 3.776822e-01          NaN
    ## [20846]          NaN          NaN 3.204838e-01          NaN          NaN
    ## [20851]          NaN          NaN 2.423756e-01          NaN 4.695656e-01
    ## [20856]          NaN          NaN 4.614343e-01 2.686784e-58          NaN
    ## [20861]          NaN          NaN 8.144527e-01          NaN          NaN
    ## [20866] 1.796279e-01 1.318001e-01          NaN          NaN          NaN
    ## [20871]          NaN 1.590102e-02 9.435884e-01 2.660725e-04 3.204838e-01
    ## [20876] 4.105277e-01          NaN          NaN          NaN 1.145293e-01
    ## [20881] 7.661717e-01 5.610814e-02 7.893008e-04 1.733747e-05 1.683603e-01
    ## [20886]          NaN          NaN 5.352753e-01 3.204838e-01 8.726903e-01
    ## [20891]          NaN          NaN          NaN          NaN 8.130129e-01
    ## [20896]          NaN          NaN 2.871775e-01 1.634264e-01          NaN
    ## [20901]          NaN          NaN          NaN          NaN          NaN
    ## [20906] 1.636764e-15          NaN 4.078514e-01 9.250427e-01          NaN
    ## [20911]          NaN 2.592642e-01 4.700178e-01          NaN          NaN
    ## [20916]          NaN          NaN          NaN          NaN          NaN
    ## [20921] 8.810482e-01          NaN 3.204838e-01          NaN          NaN
    ## [20926]          NaN 5.200602e-01          NaN          NaN          NaN
    ## [20931]          NaN          NaN          NaN 1.593386e-16          NaN
    ## [20936] 2.001535e-05 2.125663e-01          NaN 4.122888e-01          NaN
    ## [20941]          NaN 6.096862e-01          NaN          NaN          NaN
    ## [20946] 1.174821e-02 9.879438e-01          NaN          NaN          NaN
    ## [20951] 9.250912e-01 8.955768e-01 4.512626e-01          NaN 3.204838e-01
    ## [20956]          NaN 3.870700e-01          NaN          NaN          NaN
    ## [20961]          NaN          NaN          NaN          NaN 7.018417e-02
    ## [20966]          NaN 1.516960e-04          NaN          NaN          NaN
    ## [20971] 8.164605e-01          NaN 1.516885e-02 4.176040e-01          NaN
    ## [20976] 3.284534e-04 3.817489e-02 6.549717e-03          NaN 3.204838e-01
    ## [20981]          NaN          NaN 7.182390e-01          NaN          NaN
    ## [20986]          NaN 1.266802e-01          NaN 2.296420e-01          NaN
    ## [20991] 6.292326e-02 7.145834e-02          NaN 5.804884e-01          NaN
    ## [20996] 8.922192e-02          NaN          NaN          NaN          NaN
    ## [21001]          NaN          NaN          NaN 3.204838e-01          NaN
    ## [21006]          NaN          NaN          NaN 9.380847e-01          NaN
    ## [21011]          NaN          NaN 9.247367e-01 3.204838e-01 1.017407e-02
    ## [21016]          NaN 2.793190e-02          NaN          NaN          NaN
    ## [21021]          NaN 5.770172e-01          NaN          NaN          NaN
    ## [21026]          NaN          NaN 9.173550e-01 9.754480e-02          NaN
    ## [21031]          NaN 2.255301e-03          NaN          NaN 4.137745e-02
    ## [21036] 1.486532e-03 8.621923e-02          NaN 8.519489e-01          NaN
    ## [21041]          NaN 3.204838e-01          NaN 3.985374e-07 5.967487e-01
    ## [21046] 3.204838e-01          NaN 8.612483e-01 7.996365e-01 7.581429e-02
    ## [21051]          NaN 9.860367e-01 1.506064e-03          NaN 1.249513e-01
    ## [21056]          NaN          NaN          NaN          NaN          NaN
    ## [21061]          NaN 3.204838e-01          NaN          NaN          NaN
    ## [21066]          NaN          NaN          NaN 3.204838e-01          NaN
    ## [21071] 7.692892e-07 4.541634e-01          NaN 3.204838e-01          NaN
    ## [21076]          NaN          NaN 2.828920e-06 3.204838e-01          NaN
    ## [21081] 9.422719e-01 1.763826e-01          NaN          NaN 8.795121e-01
    ## [21086] 8.857344e-02 2.427559e-01 5.863632e-01          NaN          NaN
    ## [21091]          NaN 5.444757e-01 5.759016e-06 8.855567e-01          NaN
    ## [21096] 4.700917e-01          NaN 3.247590e-01          NaN          NaN
    ## [21101] 3.204838e-01          NaN 1.688158e-02          NaN          NaN
    ## [21106]          NaN 1.864001e-01          NaN          NaN          NaN
    ## [21111] 2.053724e-01 4.069254e-05 3.204838e-01          NaN          NaN
    ## [21116]          NaN 3.342226e-01          NaN 4.665900e-01 8.012803e-01
    ## [21121]          NaN 8.957075e-01 4.349515e-01          NaN          NaN
    ## [21126]          NaN 2.735287e-08 9.964920e-03          NaN          NaN
    ## [21131]          NaN 7.560370e-01 1.444487e-04          NaN 3.204838e-01
    ## [21136]          NaN          NaN 1.724888e-03          NaN          NaN
    ## [21141]          NaN 3.569724e-08          NaN          NaN 3.742867e-01
    ## [21146] 3.204838e-01 1.169672e-01          NaN 5.798168e-01 1.648174e-01
    ## [21151]          NaN          NaN          NaN          NaN          NaN
    ## [21156]          NaN 8.421449e-01          NaN 3.204838e-01          NaN
    ## [21161] 1.474731e-01 8.100520e-01          NaN          NaN          NaN
    ## [21166]          NaN          NaN 2.073119e-56 9.251497e-01          NaN
    ## [21171]          NaN          NaN          NaN 8.160563e-01 2.335159e-03
    ## [21176]          NaN          NaN          NaN 4.297820e-01          NaN
    ## [21181] 4.328403e-05          NaN          NaN          NaN 2.013843e-01
    ## [21186]          NaN          NaN          NaN          NaN 6.591851e-01
    ## [21191] 4.735536e-02 3.204838e-01          NaN          NaN 6.933886e-01
    ## [21196] 2.032624e-03 1.507001e-01          NaN 5.312296e-01          NaN
    ## [21201] 8.707474e-06 4.698270e-01          NaN 5.326502e-01          NaN
    ## [21206] 1.739035e-01          NaN          NaN          NaN          NaN
    ## [21211] 8.452112e-02 5.678763e-03          NaN          NaN 3.249332e-02
    ## [21216] 3.357127e-07 9.779618e-01 6.386251e-02 3.204838e-01 3.204838e-01
    ## [21221]          NaN 9.336477e-01          NaN 4.501946e-28 2.175284e-01
    ## [21226] 4.506760e-01          NaN 1.225055e-01          NaN          NaN
    ## [21231]          NaN          NaN 3.204838e-01          NaN 5.458688e-01
    ## [21236] 4.125497e-13 9.802275e-02          NaN 1.342540e-01          NaN
    ## [21241] 2.322270e-01          NaN          NaN 8.228071e-01 8.664498e-01
    ## [21246]          NaN 4.653996e-01          NaN 7.961282e-01 4.162950e-09
    ## [21251]          NaN 1.103712e-04          NaN          NaN 3.204838e-01
    ## [21256] 1.032670e-04          NaN 9.718136e-01 9.717079e-01 9.603054e-01
    ## [21261] 4.182925e-01 1.136489e-03 3.204838e-01          NaN 1.128527e-36
    ## [21266] 5.309550e-01          NaN          NaN 2.341983e-02 8.558843e-01
    ## [21271]          NaN          NaN          NaN 4.253996e-01 9.093622e-02
    ## [21276]          NaN          NaN 1.605839e-01          NaN 3.204838e-01
    ## [21281]          NaN          NaN          NaN 4.861643e-01          NaN
    ## [21286] 4.581687e-01 7.460531e-02          NaN          NaN          NaN
    ## [21291]          NaN 1.116943e-02          NaN 1.017629e-02          NaN
    ## [21296]          NaN          NaN          NaN 7.940139e-02 1.722024e-01
    ## [21301]          NaN          NaN          NaN          NaN 1.561691e-01
    ## [21306]          NaN          NaN 3.662588e-01          NaN 1.455828e-05
    ## [21311] 7.882651e-01 1.420384e-01 2.929533e-01          NaN          NaN
    ## [21316]          NaN          NaN 7.863198e-01 5.120521e-07 5.988786e-02
    ## [21321]          NaN 3.204838e-01 3.204838e-01          NaN          NaN
    ## [21326] 9.055255e-01 2.816565e-01          NaN          NaN 3.890278e-01
    ## [21331]          NaN 3.204838e-01          NaN          NaN 8.360771e-01
    ## [21336] 9.777578e-01 8.243692e-02          NaN          NaN          NaN
    ## [21341] 4.394256e-01          NaN 8.206976e-01          NaN          NaN
    ## [21346] 8.211725e-01          NaN          NaN 8.093934e-02 3.166921e-04
    ## [21351] 8.912189e-01          NaN          NaN 3.204838e-01          NaN
    ## [21356] 5.334695e-01          NaN          NaN          NaN          NaN
    ## [21361] 2.932674e-04 3.739130e-15 3.047470e-02          NaN          NaN
    ## [21366]          NaN 6.915410e-02 7.410691e-01 6.903993e-03 8.122100e-02
    ## [21371]          NaN 3.204838e-01          NaN          NaN          NaN
    ## [21376]          NaN          NaN 4.397920e-04 3.499295e-01 7.644208e-01
    ## [21381] 3.204838e-01          NaN 3.445742e-01          NaN 8.424423e-01
    ## [21386] 4.712004e-01 5.277450e-01          NaN 2.296457e-01 2.450193e-04
    ## [21391]          NaN          NaN 5.046541e-01 4.849812e-01 6.991557e-01
    ## [21396] 6.983565e-01 2.067615e-01 1.699402e-13 6.328549e-01          NaN
    ## [21401] 9.118693e-01          NaN 2.093218e-10 1.204077e-02          NaN
    ## [21406]          NaN          NaN          NaN 4.544304e-01          NaN
    ## [21411]          NaN          NaN          NaN 1.627747e-01          NaN
    ## [21416]          NaN          NaN          NaN 7.437735e-01 8.626793e-01
    ## [21421] 3.204838e-01          NaN          NaN 1.265014e-02 1.037752e-02
    ## [21426] 2.550068e-01          NaN 6.539934e-03          NaN 3.204838e-01
    ## [21431] 3.633949e-01          NaN 1.664510e-01          NaN 1.564672e-23
    ## [21436]          NaN          NaN          NaN 4.325050e-01 4.283896e-01
    ## [21441] 4.229806e-01          NaN          NaN          NaN 1.312042e-05
    ## [21446]          NaN          NaN          NaN          NaN          NaN
    ## [21451]          NaN          NaN 6.970834e-01 8.153343e-02          NaN
    ## [21456]          NaN 3.204838e-01 1.677848e-07 3.204838e-01          NaN
    ## [21461]          NaN          NaN          NaN 1.265857e-10          NaN
    ## [21466]          NaN          NaN 8.442201e-01 8.139592e-01 1.323886e-71
    ## [21471]          NaN 9.862027e-01 1.737044e-01          NaN          NaN
    ## [21476]          NaN 8.211476e-06          NaN 2.824459e-09          NaN
    ## [21481]          NaN          NaN 5.535835e-01          NaN          NaN
    ## [21486] 3.204838e-01          NaN          NaN 9.627688e-02          NaN
    ## [21491] 9.146662e-01 3.038477e-03          NaN 6.093881e-01          NaN
    ## [21496] 4.001326e-01 4.538880e-01          NaN          NaN          NaN
    ## [21501] 5.122729e-01          NaN          NaN 7.463420e-01 5.532171e-01
    ## [21506]          NaN 3.321722e-02 2.330463e-02          NaN          NaN
    ## [21511] 9.387330e-01 1.594655e-01 2.270288e-01          NaN          NaN
    ## [21516] 2.465164e-02          NaN 6.704384e-03          NaN          NaN
    ## [21521]          NaN 4.014789e-01          NaN          NaN 1.101620e-01
    ## [21526]          NaN          NaN          NaN 2.911746e-04          NaN
    ## [21531]          NaN          NaN          NaN          NaN 9.795948e-01
    ## [21536]          NaN          NaN 8.111669e-01          NaN          NaN
    ## [21541] 3.218546e-08          NaN 5.739091e-03          NaN 7.790138e-02
    ## [21546]          NaN          NaN 2.720262e-01          NaN 8.322035e-02
    ## [21551]          NaN 4.880223e-02 1.614683e-02          NaN          NaN
    ## [21556] 1.688781e-01          NaN          NaN          NaN          NaN
    ## [21561] 1.226034e-01          NaN          NaN 8.525865e-01          NaN
    ## [21566]          NaN          NaN 7.840748e-01          NaN          NaN
    ## [21571] 3.926238e-01          NaN          NaN          NaN 6.437314e-01
    ## [21576]          NaN 3.204838e-01          NaN          NaN 1.560496e-01
    ## [21581] 2.422402e-01          NaN 3.204838e-01 1.375817e-02          NaN
    ## [21586]          NaN          NaN          NaN 5.286855e-01 2.845172e-01
    ## [21591] 3.204838e-01          NaN          NaN          NaN          NaN
    ## [21596]          NaN 1.673985e-01 4.268416e-02          NaN 9.066630e-05
    ## [21601] 6.741434e-01          NaN 3.822283e-01          NaN          NaN
    ## [21606]          NaN 3.366048e-01 8.088516e-01 3.287144e-02          NaN
    ## [21611] 8.396172e-01          NaN          NaN          NaN          NaN
    ## [21616] 3.204838e-01          NaN 3.470954e-01          NaN 3.204838e-01
    ## [21621] 8.385914e-04 1.707484e-01 3.204838e-01          NaN          NaN
    ## [21626]          NaN          NaN          NaN 3.204838e-01          NaN
    ## [21631]          NaN 4.137507e-01 3.204838e-01 1.534922e-02          NaN
    ## [21636]          NaN 1.559032e-01          NaN          NaN          NaN
    ## [21641] 3.204838e-01          NaN 4.550838e-01          NaN 5.364205e-01
    ## [21646]          NaN 3.204838e-01 6.052582e-01 5.961490e-01          NaN
    ## [21651] 8.070471e-01          NaN          NaN          NaN 2.190669e-02
    ## [21656] 6.107295e-42          NaN          NaN          NaN          NaN
    ## [21661]          NaN          NaN 3.204838e-01 4.351167e-01          NaN
    ## [21666]          NaN 1.049397e-20          NaN          NaN          NaN
    ## [21671] 3.204838e-01          NaN 2.439891e-01          NaN 2.714012e-10
    ## [21676] 3.204838e-01          NaN          NaN          NaN          NaN
    ## [21681]          NaN          NaN 3.569675e-01          NaN          NaN
    ## [21686]          NaN          NaN          NaN 2.925860e-01          NaN
    ## [21691]          NaN          NaN 1.178343e-01          NaN 3.352200e-01
    ## [21696] 6.069647e-01 4.041025e-11 2.852252e-01 1.289972e-11          NaN
    ## [21701]          NaN 4.348975e-10 3.204838e-01 1.154420e-01 6.044780e-01
    ## [21706]          NaN          NaN 2.393363e-01          NaN          NaN
    ## [21711]          NaN          NaN 5.304086e-01          NaN          NaN
    ## [21716] 4.215204e-02          NaN 8.297104e-02          NaN 5.086655e-01
    ## [21721]          NaN          NaN          NaN 3.333189e-02 1.265152e-01
    ## [21726]          NaN          NaN          NaN 1.566671e-01          NaN
    ## [21731] 3.920568e-01 6.194069e-01 3.204838e-01 6.556935e-02          NaN
    ## [21736]          NaN 5.173958e-01          NaN          NaN          NaN
    ## [21741]          NaN          NaN 4.172625e-01          NaN 3.204838e-01
    ## [21746]          NaN 1.056883e-01 5.033103e-01          NaN          NaN
    ## [21751]          NaN 4.528094e-01          NaN 8.128930e-01          NaN
    ## [21756]          NaN          NaN 3.204838e-01 5.071531e-01 5.223867e-01
    ## [21761]          NaN          NaN 3.826689e-02 1.379917e-02          NaN
    ## [21766] 3.204838e-01          NaN 3.103743e-01 1.886380e-02 1.806848e-13
    ## [21771] 9.377455e-01          NaN 1.287815e-10          NaN          NaN
    ## [21776] 3.204838e-01 6.527187e-01          NaN          NaN 4.226353e-05
    ## [21781]          NaN          NaN          NaN          NaN          NaN
    ## [21786]          NaN          NaN          NaN 1.304304e-02          NaN
    ## [21791]          NaN 9.606863e-01          NaN          NaN 3.204838e-01
    ## [21796] 1.470757e-01          NaN          NaN          NaN 1.557135e-01
    ## [21801]          NaN          NaN 4.139944e-01          NaN          NaN
    ## [21806]          NaN          NaN 6.114475e-01 6.941183e-01          NaN
    ## [21811] 1.550621e-04 3.828831e-02 2.601676e-01 1.309574e-01 4.865924e-03
    ## [21816] 8.357679e-01 8.882051e-01          NaN          NaN          NaN
    ## [21821] 2.711504e-01          NaN          NaN          NaN 4.286755e-04
    ## [21826] 3.204838e-01 3.134133e-01          NaN 3.204838e-01          NaN
    ## [21831]          NaN          NaN 1.593962e-01          NaN          NaN
    ## [21836]          NaN          NaN 1.949298e-01          NaN          NaN
    ## [21841]          NaN          NaN          NaN 4.864987e-05 1.001828e-03
    ## [21846]          NaN          NaN          NaN 5.539699e-03          NaN
    ## [21851]          NaN 2.608531e-03          NaN          NaN          NaN
    ## [21856]          NaN          NaN 2.532662e-02 2.973744e-05 1.274215e-01
    ## [21861] 6.111855e-08 1.835795e-01 5.267976e-02 2.221719e-04          NaN
    ## [21866] 8.572047e-01 3.204838e-01 3.032679e-01          NaN          NaN
    ## [21871] 7.343700e-01 8.067687e-01          NaN          NaN 1.925642e-02
    ## [21876] 1.638085e-01          NaN 1.258143e-01          NaN          NaN
    ## [21881]          NaN 7.842266e-01 7.238846e-01          NaN 2.377704e-11
    ## [21886]          NaN          NaN 6.082548e-01 1.387014e-02 1.057499e-02
    ## [21891] 5.539501e-02 3.204838e-01 4.872243e-01          NaN 8.279494e-01
    ## [21896] 4.079477e-01 5.194493e-02          NaN          NaN          NaN
    ## [21901]          NaN 2.228454e-02          NaN          NaN          NaN
    ## [21906]          NaN          NaN          NaN 8.693563e-02          NaN
    ## [21911]          NaN          NaN          NaN 4.180651e-01          NaN
    ## [21916] 4.897506e-01 6.885229e-01          NaN          NaN          NaN
    ## [21921]          NaN          NaN          NaN          NaN          NaN
    ## [21926]          NaN 3.204838e-01          NaN          NaN 3.427140e-01
    ## [21931]          NaN          NaN          NaN          NaN          NaN
    ## [21936] 1.134359e-01 9.698851e-01          NaN          NaN 4.209706e-01
    ## [21941]          NaN 6.611517e-02          NaN 1.269549e-01          NaN
    ## [21946]          NaN          NaN          NaN          NaN          NaN
    ## [21951] 1.711158e-01          NaN          NaN 2.526074e-01          NaN
    ## [21956]          NaN          NaN 6.145204e-01 2.529754e-23          NaN
    ## [21961] 1.561702e-01          NaN          NaN          NaN          NaN
    ## [21966] 4.767940e-01 5.848671e-01          NaN          NaN          NaN
    ## [21971]          NaN 1.716447e-11          NaN          NaN          NaN
    ## [21976]          NaN          NaN 5.407091e-02          NaN          NaN
    ## [21981] 2.699589e-01 3.204838e-01 3.204838e-01          NaN 4.125889e-01
    ## [21986] 8.092547e-01 3.204838e-01 5.625938e-01 2.020474e-01          NaN
    ## [21991] 3.204838e-01          NaN 3.405826e-01          NaN          NaN
    ## [21996]          NaN          NaN          NaN          NaN 2.558413e-01
    ## [22001] 4.858905e-03          NaN 9.511008e-01 9.975198e-01          NaN
    ## [22006]          NaN          NaN 8.942974e-01          NaN 6.096483e-01
    ## [22011]          NaN          NaN          NaN          NaN 2.100403e-01
    ## [22016]          NaN          NaN 7.526236e-03          NaN 3.204838e-01
    ## [22021] 4.955701e-01          NaN 1.798965e-04          NaN          NaN
    ## [22026] 2.868581e-01 1.347570e-17 4.062286e-01 1.862100e-02 8.244867e-01
    ## [22031] 3.929969e-01          NaN          NaN 3.204838e-01 5.332496e-01
    ## [22036]          NaN          NaN          NaN          NaN          NaN
    ## [22041]          NaN          NaN          NaN          NaN 3.388325e-01
    ## [22046]          NaN          NaN 3.141661e-02 6.564066e-07 5.452455e-01
    ## [22051]          NaN          NaN          NaN          NaN 1.524560e-01
    ## [22056]          NaN          NaN 3.204838e-01          NaN 3.204838e-01
    ## [22061]          NaN 3.204838e-01          NaN          NaN          NaN
    ## [22066] 3.204838e-01          NaN 2.708639e-01          NaN          NaN
    ## [22071] 4.600542e-01          NaN 3.204838e-01 2.097620e-02          NaN
    ## [22076] 3.204838e-01          NaN          NaN 7.100883e-01 7.760001e-01
    ## [22081]          NaN          NaN          NaN 8.816244e-02          NaN
    ## [22086]          NaN          NaN          NaN 2.968413e-01 2.601486e-01
    ## [22091]          NaN 1.655574e-01          NaN          NaN 4.583494e-03
    ## [22096]          NaN          NaN          NaN 3.204838e-01          NaN
    ## [22101]          NaN 3.204838e-01 2.005641e-01          NaN          NaN
    ## [22106]          NaN          NaN          NaN          NaN          NaN
    ## [22111] 9.284075e-01          NaN          NaN          NaN          NaN
    ## [22116]          NaN          NaN 6.727656e-03 9.962204e-01 6.274614e-01
    ## [22121] 6.209926e-01 5.041084e-01 1.018916e-01 3.204838e-01 2.012286e-15
    ## [22126]          NaN          NaN          NaN          NaN          NaN
    ## [22131]          NaN          NaN 2.942654e-01          NaN 1.981815e-01
    ## [22136]          NaN          NaN 1.611733e-01          NaN 1.228305e-04
    ## [22141]          NaN 4.373498e-01 2.771004e-01          NaN          NaN
    ## [22146] 3.204838e-01          NaN 7.882275e-01 8.294012e-01 3.578485e-01
    ## [22151]          NaN          NaN          NaN 4.560298e-02          NaN
    ## [22156]          NaN          NaN          NaN          NaN          NaN
    ## [22161]          NaN 9.122532e-01 7.272247e-01 1.254880e-01 5.252539e-01
    ## [22166]          NaN 3.995639e-03 3.204838e-01 5.190831e-05 9.573398e-01
    ## [22171] 5.733049e-01 4.767977e-01          NaN 9.544415e-01          NaN
    ## [22176]          NaN          NaN          NaN 3.097677e-01          NaN
    ## [22181]          NaN          NaN          NaN 7.972982e-12          NaN
    ## [22186] 8.238609e-03          NaN          NaN          NaN          NaN
    ## [22191] 3.862429e-01 2.259066e-02 9.613978e-02          NaN 8.722299e-01
    ## [22196] 9.780128e-01          NaN 1.722521e-02 1.667707e-01 9.102068e-01
    ## [22201] 4.893687e-02 1.926821e-01          NaN 9.781823e-03 1.085164e-02
    ## [22206] 4.159101e-02 2.315578e-01 2.775501e-12 9.149124e-02 9.417619e-01
    ## [22211]          NaN          NaN 1.864176e-03          NaN 3.204838e-01
    ## [22216] 1.401801e-01 7.157169e-01          NaN 2.249639e-01          NaN
    ## [22221]          NaN          NaN          NaN          NaN          NaN
    ## [22226] 3.204838e-01 5.396694e-02 1.610215e-01 3.204838e-01 3.204838e-01
    ## [22231] 3.789631e-01 3.204838e-01          NaN          NaN          NaN
    ## [22236]          NaN          NaN          NaN 8.008254e-02 2.269202e-01
    ## [22241]          NaN          NaN          NaN 3.204838e-01 1.079839e-01
    ## [22246]          NaN          NaN          NaN          NaN 3.204838e-01
    ## [22251] 8.645602e-01          NaN 2.301710e-02          NaN          NaN
    ## [22256]          NaN 9.943826e-01 1.113179e-05 7.236589e-01 1.623861e-01
    ## [22261] 3.204838e-01          NaN 5.037361e-01 3.204838e-01 4.912451e-02
    ## [22266]          NaN          NaN          NaN          NaN 6.534501e-01
    ## [22271] 6.001007e-01 7.292417e-01          NaN          NaN          NaN
    ## [22276] 3.204838e-01 4.301557e-01          NaN 1.159139e-04 4.381950e-01
    ## [22281]          NaN          NaN 2.451429e-01          NaN          NaN
    ## [22286] 5.216930e-01 3.204838e-01          NaN 5.619436e-01          NaN
    ## [22291]          NaN          NaN          NaN          NaN 1.577781e-01
    ## [22296]          NaN 3.204838e-01 1.593795e-01          NaN 3.525634e-03
    ## [22301]          NaN          NaN          NaN          NaN 9.757962e-01
    ## [22306]          NaN          NaN 8.095320e-01          NaN          NaN
    ## [22311]          NaN          NaN          NaN          NaN          NaN
    ## [22316] 8.295716e-01          NaN          NaN          NaN          NaN
    ## [22321] 2.204916e-01          NaN 1.136800e-01 1.620693e-02          NaN
    ## [22326] 8.512269e-01          NaN 6.688386e-01          NaN          NaN
    ## [22331]          NaN 2.163747e-05          NaN 1.559060e-01          NaN
    ## [22336] 2.047825e-04 7.099708e-01 3.204838e-01 3.204838e-01 1.710568e-02
    ## [22341]          NaN          NaN          NaN 1.870398e-04 3.204838e-01
    ## [22346]          NaN 8.246929e-12 8.823796e-04 3.204838e-01 7.488265e-01
    ## [22351] 1.365087e-01          NaN 7.371070e-01 1.543672e-01 2.417415e-02
    ## [22356] 4.480055e-24          NaN          NaN          NaN          NaN
    ## [22361] 8.717367e-01          NaN 8.788132e-01          NaN          NaN
    ## [22366]          NaN          NaN          NaN          NaN          NaN
    ## [22371] 2.428650e-01 3.776822e-01 2.470846e-01 1.570187e-01          NaN
    ## [22376]          NaN          NaN 3.016923e-03          NaN          NaN
    ## [22381]          NaN          NaN 3.204838e-01 4.206650e-01          NaN
    ## [22386] 8.796629e-01 4.534887e-15          NaN 3.204838e-01          NaN
    ## [22391]          NaN          NaN          NaN          NaN 6.385951e-02
    ## [22396] 4.064428e-01 2.244781e-13 6.603745e-01          NaN          NaN
    ## [22401]          NaN 4.673608e-01          NaN          NaN          NaN
    ## [22406]          NaN          NaN 5.122684e-01          NaN          NaN
    ## [22411]          NaN          NaN          NaN 3.094171e-01 3.204838e-01
    ## [22416] 7.568619e-01          NaN 1.243129e-02          NaN 7.774801e-01
    ## [22421] 2.388962e-01          NaN          NaN 1.559139e-01 2.662986e-03
    ## [22426]          NaN          NaN 1.297555e-03 1.151220e-01 8.456156e-01
    ## [22431] 3.204838e-01          NaN          NaN 7.599926e-01          NaN
    ## [22436] 1.322588e-02          NaN          NaN          NaN 4.959595e-01
    ## [22441]          NaN          NaN          NaN          NaN          NaN
    ## [22446]          NaN 9.420697e-01 3.204838e-01          NaN 1.874134e-01
    ## [22451]          NaN          NaN 7.013105e-01          NaN 4.132019e-01
    ## [22456] 3.204838e-01          NaN          NaN 5.150186e-01 1.559081e-01
    ## [22461]          NaN          NaN          NaN          NaN 8.580639e-02
    ## [22466] 2.336383e-04          NaN          NaN 3.463607e-07          NaN
    ## [22471] 5.896988e-01          NaN          NaN 8.303430e-02          NaN
    ## [22476] 7.741777e-01          NaN 7.758284e-01          NaN 1.734960e-01
    ## [22481]          NaN          NaN 4.997455e-01          NaN          NaN
    ## [22486]          NaN          NaN          NaN          NaN          NaN
    ## [22491] 1.559369e-01 2.342434e-02 3.204838e-01          NaN          NaN
    ## [22496]          NaN          NaN          NaN 8.008145e-01          NaN
    ## [22501] 2.245783e-08 3.204838e-01          NaN 3.204838e-01          NaN
    ## [22506]          NaN 9.975603e-07 1.352150e-09 1.282168e-01 1.037447e-05
    ## [22511]          NaN 3.749079e-01          NaN          NaN          NaN
    ## [22516] 1.087575e-01          NaN          NaN          NaN          NaN
    ## [22521] 3.218810e-02 4.344587e-01          NaN          NaN          NaN
    ## [22526] 1.008072e-06 8.479765e-01          NaN 9.772988e-01          NaN
    ## [22531] 7.454992e-01 3.204838e-01 4.042216e-01 8.159118e-03          NaN
    ## [22536]          NaN          NaN          NaN          NaN          NaN
    ## [22541]          NaN 1.600075e-01 3.204838e-01          NaN          NaN
    ## [22546]          NaN 8.901073e-01          NaN          NaN          NaN
    ## [22551]          NaN          NaN          NaN 6.924410e-01          NaN
    ## [22556] 5.659078e-03          NaN          NaN 4.569169e-01 9.995932e-01
    ## [22561] 3.204838e-01 3.005510e-01          NaN          NaN 8.154117e-02
    ## [22566]          NaN          NaN          NaN          NaN 6.466669e-01
    ## [22571] 1.096682e-02          NaN          NaN          NaN          NaN
    ## [22576] 1.491243e-01 6.066628e-01 9.889190e-03          NaN 3.204838e-01
    ## [22581] 4.300734e-01          NaN          NaN 7.960888e-01          NaN
    ## [22586]          NaN          NaN          NaN          NaN          NaN
    ## [22591]          NaN          NaN 3.008415e-01 1.566292e-01          NaN
    ## [22596] 1.928704e-01 2.833259e-01 2.349921e-02 2.760110e-03          NaN
    ## [22601] 1.581878e-02 1.334343e-01          NaN 5.738586e-01          NaN
    ## [22606] 1.141440e-02          NaN 6.627533e-01          NaN          NaN
    ## [22611] 5.474482e-01 3.204838e-01          NaN 4.605924e-01 1.240006e-06
    ## [22616] 1.359088e-02          NaN          NaN          NaN          NaN
    ## [22621]          NaN          NaN 9.483124e-01          NaN          NaN
    ## [22626] 1.191180e-01          NaN          NaN          NaN          NaN
    ## [22631] 1.699870e-01          NaN 3.755902e-01 8.450506e-03 1.266798e-02
    ## [22636]          NaN 7.068541e-01          NaN 2.079148e-07 3.204838e-01
    ## [22641]          NaN          NaN 8.453981e-01 2.020242e-01 3.365167e-01
    ## [22646] 3.204838e-01 2.216510e-01          NaN          NaN          NaN
    ## [22651]          NaN 1.797792e-01          NaN 6.492289e-01 1.652271e-01
    ## [22656]          NaN          NaN          NaN          NaN 5.978905e-01
    ## [22661]          NaN          NaN 1.670234e-02 4.613070e-02 6.524420e-02
    ## [22666] 2.080982e-01 4.801333e-01          NaN          NaN          NaN
    ## [22671] 9.789398e-01 1.187507e-01 3.204838e-01          NaN          NaN
    ## [22676] 7.266011e-01          NaN 1.748269e-01          NaN          NaN
    ## [22681] 4.668142e-01 3.204838e-01 6.603745e-01          NaN 1.231937e-01
    ## [22686] 4.574561e-02          NaN          NaN 2.959507e-01          NaN
    ## [22691]          NaN 4.278437e-03          NaN 1.373101e-01 3.120908e-01
    ## [22696] 2.747510e-01          NaN          NaN 7.095001e-01 1.740763e-01
    ## [22701] 1.778365e-01 2.602731e-01          NaN          NaN          NaN
    ## [22706]          NaN          NaN 3.944745e-01          NaN          NaN
    ## [22711]          NaN          NaN 3.071279e-01          NaN          NaN
    ## [22716]          NaN          NaN 1.363833e-04 1.275180e-02 6.852666e-01
    ## [22721]          NaN 2.343484e-01 3.151713e-01          NaN          NaN
    ## [22726] 2.027541e-01          NaN 9.358552e-01          NaN          NaN
    ## [22731] 3.382066e-01          NaN          NaN          NaN          NaN
    ## [22736] 3.204838e-01          NaN 4.684734e-02          NaN          NaN
    ## [22741] 6.550287e-01          NaN 5.563286e-01          NaN          NaN
    ## [22746] 5.294960e-12 7.262359e-07          NaN          NaN          NaN
    ## [22751] 3.204838e-01          NaN          NaN 5.176034e-02          NaN
    ## [22756]          NaN 3.524495e-01          NaN 3.590162e-39 7.075464e-01
    ## [22761]          NaN 7.004511e-01 1.618941e-01          NaN          NaN
    ## [22766] 8.174823e-02          NaN          NaN 1.141276e-03          NaN
    ## [22771] 6.914945e-03 1.538079e-01          NaN 1.589661e-03          NaN
    ## [22776]          NaN          NaN          NaN          NaN 1.277525e-01
    ## [22781]          NaN          NaN          NaN          NaN          NaN
    ## [22786] 3.204838e-01 8.828686e-05 3.014436e-02 8.399644e-01          NaN
    ## [22791]          NaN          NaN 2.461662e-03          NaN          NaN
    ## [22796] 9.943849e-04          NaN 5.921416e-01          NaN          NaN
    ## [22801]          NaN 3.776822e-01          NaN          NaN 6.114145e-01
    ## [22806] 5.322276e-01          NaN 6.237836e-02          NaN 1.365728e-02
    ## [22811]          NaN          NaN          NaN          NaN          NaN
    ## [22816] 8.342613e-01 3.204838e-01          NaN          NaN          NaN
    ## [22821]          NaN 6.487325e-01          NaN          NaN          NaN
    ## [22826]          NaN          NaN          NaN          NaN          NaN
    ## [22831] 5.921461e-01          NaN 1.109310e-01 4.632512e-01 5.826210e-15
    ## [22836] 1.147872e-01 5.698894e-02          NaN 7.677608e-01          NaN
    ## [22841]          NaN 5.693558e-01          NaN 3.446758e-01 6.154976e-03
    ## [22846] 1.059402e-05          NaN 6.567466e-04          NaN          NaN
    ## [22851] 4.683239e-01 2.194574e-10          NaN 2.276469e-01 3.773567e-02
    ## [22856] 3.204838e-01          NaN 6.426003e-07 8.838170e-01          NaN
    ## [22861]          NaN          NaN 1.878944e-01 1.564033e-01 1.590278e-01
    ## [22866]          NaN 3.393195e-01          NaN 8.020743e-01 4.561040e-02
    ## [22871] 5.827764e-01          NaN 4.997391e-01          NaN          NaN
    ## [22876] 3.204838e-01 1.449441e-03          NaN 3.639532e-01 1.504004e-01
    ## [22881]          NaN 3.204838e-01 6.688904e-01          NaN 3.204838e-01
    ## [22886]          NaN          NaN 9.497456e-01 1.866501e-04          NaN
    ## [22891]          NaN          NaN          NaN          NaN          NaN
    ## [22896]          NaN 3.204838e-01          NaN          NaN 8.827131e-01
    ## [22901]          NaN 4.105127e-05 6.471775e-01          NaN 2.301342e-01
    ## [22906]          NaN          NaN          NaN 1.558921e-02 3.204838e-01
    ## [22911] 6.472372e-01 3.618578e-02 2.166486e-03          NaN          NaN
    ## [22916] 7.393915e-01          NaN          NaN          NaN          NaN
    ## [22921] 6.266128e-03          NaN          NaN          NaN          NaN
    ## [22926]          NaN          NaN          NaN          NaN 1.203391e-03
    ## [22931]          NaN          NaN          NaN 1.333673e-01 4.017481e-03
    ## [22936]          NaN          NaN 2.265426e-01          NaN          NaN
    ## [22941]          NaN 1.562167e-01 2.711755e-01          NaN          NaN
    ## [22946] 7.957103e-01 4.723021e-01          NaN 8.484252e-01          NaN
    ## [22951] 3.711159e-01 9.545585e-03 1.261475e-01          NaN 1.484744e-01
    ## [22956] 4.785212e-01          NaN          NaN          NaN          NaN
    ## [22961]          NaN 3.504430e-05 1.582728e-01 4.652771e-02 6.539326e-01
    ## [22966] 5.158994e-01 6.945896e-01          NaN          NaN 7.971879e-01
    ## [22971]          NaN          NaN          NaN 3.204838e-01          NaN
    ## [22976]          NaN          NaN          NaN          NaN 5.127771e-02
    ## [22981] 7.080040e-04          NaN          NaN          NaN          NaN
    ## [22986]          NaN          NaN 1.602198e-01 1.466700e-02 9.758964e-01
    ## [22991]          NaN 7.391863e-01 7.437835e-01 3.204838e-01 3.204838e-01
    ## [22996]          NaN 9.283608e-01 2.433093e-01          NaN          NaN
    ## [23001] 6.673932e-01          NaN 3.360102e-09          NaN          NaN
    ## [23006] 9.450605e-01 9.069959e-01 1.191012e-02 9.980836e-01          NaN
    ## [23011] 5.234536e-02 1.253389e-01          NaN          NaN          NaN
    ## [23016] 1.969691e-01          NaN 3.686671e-02          NaN 9.999721e-01
    ## [23021]          NaN          NaN          NaN 4.009472e-01          NaN
    ## [23026] 9.888970e-01 5.562396e-01          NaN 1.153236e-02 3.204838e-01
    ## [23031] 1.108660e-04          NaN 8.341235e-03          NaN 1.167063e-01
    ## [23036] 1.339793e-04          NaN 3.249418e-01          NaN          NaN
    ## [23041] 3.164681e-01          NaN          NaN 6.144998e-02 5.183523e-01
    ## [23046] 7.878781e-02 4.931248e-01          NaN 3.976316e-02 3.417296e-03
    ## [23051]          NaN 6.731620e-03          NaN          NaN          NaN
    ## [23056] 5.605355e-01          NaN          NaN 1.536251e-01          NaN
    ## [23061] 3.204838e-01          NaN          NaN          NaN 2.886167e-01
    ## [23066] 6.701475e-02 4.239076e-01 9.966063e-01          NaN          NaN
    ## [23071]          NaN 6.831604e-01 9.121340e-27          NaN 1.614148e-01
    ## [23076] 4.985939e-16 4.110252e-01          NaN          NaN          NaN
    ## [23081]          NaN          NaN          NaN          NaN 1.082722e-01
    ## [23086] 9.372427e-06 3.776822e-01          NaN          NaN          NaN
    ## [23091]          NaN 3.534853e-03 6.235023e-02          NaN          NaN
    ## [23096] 2.569531e-01 3.772809e-01          NaN          NaN          NaN
    ## [23101] 2.433390e-01 2.862991e-02          NaN 1.645605e-02          NaN
    ## [23106] 1.636538e-04          NaN          NaN 1.845940e-04          NaN
    ## [23111] 5.254580e-04          NaN 5.134668e-02 8.063124e-01          NaN
    ## [23116] 7.147109e-02          NaN          NaN 4.432657e-01 9.788746e-01
    ## [23121] 9.519841e-01 7.313532e-23          NaN 4.865781e-01          NaN
    ## [23126]          NaN 3.236345e-05          NaN 2.272640e-02 8.124769e-02
    ## [23131] 1.559796e-01          NaN 1.340019e-04 3.204838e-01          NaN
    ## [23136]          NaN 5.466119e-01          NaN          NaN          NaN
    ## [23141]          NaN 3.204838e-01 3.204838e-01 3.988777e-06          NaN
    ## [23146]          NaN          NaN 9.976555e-04          NaN 1.975382e-01
    ## [23151]          NaN          NaN 9.425724e-01          NaN          NaN
    ## [23156] 5.891189e-13          NaN          NaN          NaN 2.169725e-01
    ## [23161] 1.680269e-01          NaN 2.666327e-01          NaN          NaN
    ## [23166]          NaN 9.487194e-01 8.312598e-01          NaN          NaN
    ## [23171]          NaN          NaN 1.843236e-03          NaN 5.340763e-01
    ## [23176]          NaN          NaN          NaN          NaN          NaN
    ## [23181] 4.139978e-01          NaN          NaN 1.559340e-01 4.036664e-01
    ## [23186] 3.204838e-01          NaN 5.728489e-01          NaN 1.500723e-08
    ## [23191] 8.537215e-02 7.884943e-01 7.491314e-01          NaN 2.867730e-01
    ## [23196] 1.017827e-02 5.117933e-03          NaN          NaN          NaN
    ## [23201]          NaN          NaN          NaN          NaN 2.705179e-01
    ## [23206]          NaN          NaN 1.806237e-01 6.513182e-01 8.170538e-02
    ## [23211]          NaN 2.444428e-02          NaN 4.663770e-01 1.007170e-78
    ## [23216] 3.204838e-01 3.204838e-01 1.315295e-06          NaN          NaN
    ## [23221]          NaN          NaN          NaN          NaN 9.286003e-04
    ## [23226] 4.932710e-02 8.123042e-01          NaN          NaN          NaN
    ## [23231] 1.565658e-01          NaN 9.011346e-01          NaN          NaN
    ## [23236]          NaN          NaN          NaN          NaN 2.639113e-01
    ## [23241]          NaN 3.596337e-01 3.204838e-01 4.351231e-01          NaN
    ## [23246]          NaN          NaN 1.555866e-03 4.649477e-03          NaN
    ## [23251] 1.358890e-01          NaN          NaN          NaN 6.765414e-01
    ## [23256]          NaN          NaN          NaN          NaN 1.601126e-05
    ## [23261]          NaN 1.669098e-06          NaN          NaN          NaN
    ## [23266] 7.720515e-06 2.060273e-01          NaN          NaN          NaN
    ## [23271] 3.729237e-01          NaN 1.143182e-02 8.866835e-05          NaN
    ## [23276] 9.283000e-01          NaN          NaN 1.229061e-07          NaN
    ## [23281]          NaN          NaN 3.012185e-01          NaN          NaN
    ## [23286]          NaN 9.679122e-01 5.530207e-01 4.953694e-01 9.146461e-01
    ## [23291]          NaN 2.813092e-01          NaN          NaN          NaN
    ## [23296] 5.044725e-01          NaN 5.114560e-01          NaN          NaN
    ## [23301] 7.632294e-01 3.204838e-01          NaN 9.588360e-01 9.501696e-02
    ## [23306]          NaN          NaN 1.117013e-06          NaN          NaN
    ## [23311]          NaN          NaN          NaN 3.442817e-02          NaN
    ## [23316]          NaN          NaN          NaN          NaN          NaN
    ## [23321]          NaN          NaN 8.773601e-04          NaN 8.128894e-01
    ## [23326] 3.684246e-05          NaN          NaN          NaN          NaN
    ## [23331] 9.886515e-01          NaN 6.200947e-02          NaN          NaN
    ## [23336] 3.204838e-01 9.515649e-01          NaN          NaN          NaN
    ## [23341] 3.204838e-01 3.186034e-03          NaN 1.231881e-15          NaN
    ## [23346]          NaN          NaN          NaN          NaN          NaN
    ## [23351] 1.231100e-01 5.035225e-10 4.428504e-01          NaN          NaN
    ## [23356] 4.056211e-01          NaN          NaN          NaN          NaN
    ## [23361]          NaN          NaN 4.087872e-01          NaN 9.437586e-01
    ## [23366] 9.252256e-01          NaN          NaN          NaN 3.100225e-01
    ## [23371]          NaN 4.507516e-01 1.041717e-04          NaN          NaN
    ## [23376]          NaN          NaN 6.324096e-01          NaN          NaN
    ## [23381] 5.652436e-12 7.071513e-01 1.531493e-01          NaN          NaN
    ## [23386]          NaN 1.623348e-02          NaN          NaN 3.062966e-01
    ## [23391]          NaN          NaN          NaN 5.312606e-01          NaN
    ## [23396]          NaN          NaN 5.055232e-01 2.271469e-01          NaN
    ## [23401] 5.367094e-01 9.878710e-01          NaN 7.642291e-02 3.572522e-01
    ## [23406]          NaN 8.968456e-02 4.006181e-02          NaN          NaN
    ## [23411] 1.510450e-29 6.165187e-01 1.419523e-01          NaN          NaN
    ## [23416] 1.183406e-03 4.892330e-01          NaN          NaN 3.255817e-01
    ## [23421]          NaN          NaN 1.213383e-01 1.860846e-02 1.695024e-01
    ## [23426]          NaN          NaN          NaN 3.336843e-01          NaN
    ## [23431]          NaN 6.875297e-02          NaN 4.171954e-03 4.966575e-01
    ## [23436] 8.061140e-01 3.204838e-01          NaN          NaN          NaN
    ## [23441]          NaN          NaN 1.639726e-01          NaN          NaN
    ## [23446]          NaN 6.124431e-05          NaN          NaN          NaN
    ## [23451]          NaN 5.839179e-04          NaN 1.023650e-06          NaN
    ## [23456] 1.753510e-02 9.319982e-01 5.327054e-02          NaN          NaN
    ## [23461] 3.204838e-01          NaN 1.810169e-01          NaN 9.580025e-01
    ## [23466]          NaN          NaN          NaN          NaN 1.705143e-01
    ## [23471] 4.994925e-01 3.388502e-01          NaN          NaN 3.199944e-01
    ## [23476] 1.765327e-02 9.711963e-81 4.258455e-06 9.677746e-01          NaN
    ## [23481] 3.133151e-03 5.604693e-01          NaN          NaN          NaN
    ## [23486]          NaN 3.204838e-01          NaN 7.526132e-01 6.651163e-01
    ## [23491]          NaN          NaN 6.162713e-01 1.464747e-03          NaN
    ## [23496] 1.190841e-01          NaN          NaN          NaN          NaN
    ## [23501] 9.675003e-01          NaN          NaN          NaN          NaN
    ## [23506] 6.824446e-01 1.635352e-02          NaN          NaN 4.293948e-02
    ## [23511] 3.204838e-01          NaN          NaN 9.725451e-01 4.485244e-01
    ## [23516]          NaN          NaN 8.981591e-74          NaN          NaN
    ## [23521]          NaN 3.204838e-01          NaN          NaN 2.212044e-01
    ## [23526] 9.214441e-02          NaN 1.918661e-01 3.204838e-01 3.485263e-01
    ## [23531] 1.325168e-03          NaN 8.115875e-01 9.526911e-01          NaN
    ## [23536] 3.242961e-07          NaN 8.299534e-01 5.487048e-08          NaN
    ## [23541] 1.511656e-02          NaN          NaN          NaN          NaN
    ## [23546]          NaN 7.790086e-02          NaN          NaN 7.361802e-01
    ## [23551] 4.125353e-06          NaN 1.797161e-12 1.318087e-01          NaN
    ## [23556]          NaN 8.064320e-01 1.354197e-06          NaN          NaN
    ## [23561]          NaN 3.204838e-01 5.454102e-01 1.560382e-01 2.292796e-10
    ## [23566] 1.111088e-01          NaN 8.877873e-02          NaN          NaN
    ## [23571]          NaN          NaN 2.704719e-01          NaN          NaN
    ## [23576] 7.139593e-01 4.076019e-42          NaN          NaN          NaN
    ## [23581] 1.207124e-01          NaN          NaN 1.555184e-01 1.635396e-01
    ## [23586]          NaN 9.895700e-01 5.741266e-03 2.874651e-02          NaN
    ## [23591] 5.589407e-01 1.487950e-02          NaN 3.964393e-01          NaN
    ## [23596]          NaN 6.151098e-02 2.493761e-01 8.579689e-01 8.214552e-68
    ## [23601] 8.818792e-01 1.751417e-01          NaN          NaN 8.115219e-08
    ## [23606]          NaN          NaN 5.327928e-17          NaN 7.601175e-01
    ## [23611]          NaN          NaN 9.805910e-01 5.996496e-01          NaN
    ## [23616]          NaN          NaN 5.360167e-01          NaN          NaN
    ## [23621]          NaN          NaN 1.535545e-01          NaN 8.155634e-02
    ## [23626]          NaN          NaN 3.490351e-01          NaN 4.316294e-02
    ## [23631] 3.204838e-01 2.006147e-01          NaN          NaN          NaN
    ## [23636] 2.122713e-03          NaN          NaN 8.261869e-10          NaN
    ## [23641] 1.585642e-02          NaN 6.840628e-01 6.685672e-01          NaN
    ## [23646]          NaN          NaN 3.357696e-10 4.739993e-01 1.508279e-01
    ## [23651] 2.814707e-02 8.936024e-02          NaN          NaN          NaN
    ## [23656]          NaN          NaN          NaN 9.932693e-01 7.649728e-01
    ## [23661] 5.290375e-01          NaN 8.767627e-01          NaN          NaN
    ## [23666]          NaN          NaN          NaN 2.413355e-10 8.289838e-01
    ## [23671] 8.445613e-01          NaN 5.459249e-02          NaN          NaN
    ## [23676]          NaN          NaN          NaN          NaN          NaN
    ## [23681]          NaN          NaN          NaN 1.942336e-02          NaN
    ## [23686] 5.486922e-01 8.262070e-01 8.718612e-01          NaN          NaN
    ## [23691]          NaN          NaN 8.772076e-47          NaN          NaN
    ## [23696]          NaN          NaN          NaN          NaN          NaN
    ## [23701]          NaN          NaN          NaN 2.651922e-01 4.190597e-03
    ## [23706]          NaN 6.282016e-05          NaN          NaN          NaN
    ## [23711] 1.448845e-01 2.946208e-07 9.632266e-01 5.874289e-01          NaN
    ## [23716] 3.204838e-01 1.563038e-01 3.204838e-01 3.204838e-01 8.386532e-02
    ## [23721] 8.975422e-02          NaN 3.666305e-01          NaN 5.558678e-01
    ## [23726] 2.312661e-43          NaN          NaN          NaN          NaN
    ## [23731] 9.007833e-02 8.905344e-01          NaN          NaN 9.067600e-01
    ## [23736] 7.504952e-03          NaN          NaN          NaN 2.507284e-02
    ## [23741]          NaN          NaN          NaN          NaN          NaN
    ## [23746] 1.752346e-01          NaN          NaN          NaN 3.204838e-01
    ## [23751]          NaN          NaN          NaN          NaN          NaN
    ## [23756] 2.147542e-01          NaN 4.840859e-02          NaN          NaN
    ## [23761]          NaN          NaN 5.164654e-04 1.824033e-02          NaN
    ## [23766] 3.204838e-01          NaN          NaN 8.653253e-01          NaN
    ## [23771]          NaN 2.382706e-02          NaN          NaN 7.873128e-04
    ## [23776]          NaN 4.724791e-01          NaN 9.323137e-01 8.737458e-01
    ## [23781]          NaN          NaN 3.204838e-01 9.639973e-05 2.025259e-08
    ## [23786]          NaN          NaN 3.102481e-03          NaN          NaN
    ## [23791]          NaN          NaN          NaN 7.958102e-01          NaN
    ## [23796]          NaN 5.920558e-06 8.199455e-01 9.199926e-01 3.204838e-01
    ## [23801]          NaN          NaN 1.696246e-06 3.204838e-01 4.373509e-02
    ## [23806] 3.373081e-01 3.204838e-01          NaN 4.057496e-14 7.354446e-01
    ## [23811]          NaN          NaN 3.925308e-02          NaN 5.419346e-01
    ## [23816]          NaN          NaN          NaN          NaN 1.572690e-01
    ## [23821]          NaN 2.602731e-01          NaN          NaN 1.384759e-01
    ## [23826]          NaN          NaN 1.212812e-01 4.499823e-02          NaN
    ## [23831]          NaN 3.068780e-01          NaN          NaN 3.591510e-01
    ## [23836] 5.819257e-01          NaN          NaN 6.333555e-02          NaN
    ## [23841]          NaN 1.123654e-08          NaN          NaN 6.612884e-01
    ## [23846]          NaN          NaN 1.379809e-02          NaN          NaN
    ## [23851]          NaN          NaN 9.522157e-01          NaN 7.757822e-01
    ## [23856] 1.032300e-01          NaN          NaN 6.271943e-02          NaN
    ## [23861] 3.204838e-01          NaN 4.340274e-01          NaN 6.219519e-01
    ## [23866] 6.335408e-19 8.208556e-01          NaN          NaN 3.204838e-01
    ## [23871]          NaN          NaN 5.197448e-01          NaN          NaN
    ## [23876]          NaN          NaN          NaN 9.742266e-01          NaN
    ## [23881]          NaN 3.541759e-01          NaN 1.781017e-04 7.427845e-01
    ## [23886]          NaN          NaN 3.204838e-01          NaN 2.092352e-02
    ## [23891] 9.116343e-01          NaN          NaN 8.067101e-01          NaN
    ## [23896] 3.382316e-10          NaN          NaN 5.857029e-01 5.566038e-05
    ## [23901] 3.168389e-01          NaN 9.535713e-01 3.181668e-01 1.757069e-02
    ## [23906] 4.112146e-02          NaN          NaN          NaN 5.390864e-01
    ## [23911]          NaN          NaN 4.576697e-01 2.332593e-05          NaN
    ## [23916] 9.764178e-01          NaN 2.017959e-01          NaN          NaN
    ## [23921]          NaN 9.380847e-01          NaN          NaN 3.204838e-01
    ## [23926] 3.614736e-01          NaN          NaN          NaN          NaN
    ## [23931]          NaN          NaN          NaN          NaN 9.202489e-01
    ## [23936] 3.681117e-04          NaN 2.781992e-43 9.022252e-03 7.160700e-02
    ## [23941]          NaN          NaN          NaN          NaN          NaN
    ## [23946] 8.290733e-01 5.255889e-01          NaN 2.988662e-02 5.232631e-03
    ## [23951]          NaN          NaN 3.777334e-01 5.525102e-18 3.204838e-01
    ## [23956]          NaN 3.789732e-02 2.054367e-02          NaN 7.676674e-01
    ## [23961]          NaN          NaN          NaN 2.369735e-01 8.555749e-01
    ## [23966] 4.358528e-01          NaN 3.204838e-01          NaN 2.965362e-02
    ## [23971] 2.138927e-01 1.822156e-01          NaN          NaN          NaN
    ## [23976]          NaN          NaN          NaN          NaN          NaN
    ## [23981]          NaN 9.517807e-01 1.354096e-01 3.493509e-01          NaN
    ## [23986] 5.374605e-01          NaN 2.182153e-02          NaN          NaN
    ## [23991] 9.416882e-01          NaN 9.973778e-01          NaN 8.564100e-01
    ## [23996] 7.288892e-03 3.565274e-01          NaN          NaN 9.600571e-01
    ## [24001]          NaN          NaN          NaN          NaN          NaN
    ## [24006] 3.896464e-09          NaN          NaN          NaN          NaN
    ## [24011]          NaN 4.092503e-07          NaN          NaN          NaN
    ## [24016]          NaN          NaN 1.383712e-02 6.704684e-02          NaN
    ## [24021]          NaN          NaN 1.071561e-01 4.921515e-02          NaN
    ## [24026] 8.621643e-47 6.292543e-01 8.648276e-01 2.335211e-01 9.513484e-01
    ## [24031] 3.204838e-01 7.234953e-03 4.663159e-01          NaN          NaN
    ## [24036]          NaN          NaN 6.987886e-01          NaN          NaN
    ## [24041] 3.204838e-01 2.625538e-03          NaN          NaN          NaN
    ## [24046]          NaN 1.233518e-01 2.085801e-02          NaN 1.567495e-01
    ## [24051] 4.207100e-02 6.358661e-15 1.088547e-03          NaN 4.735887e-02
    ## [24056] 2.186769e-08 3.594018e-01          NaN 1.560992e-01 3.575208e-02
    ## [24061] 4.282861e-02          NaN          NaN          NaN          NaN
    ## [24066] 3.967913e-18 6.898909e-01          NaN 1.870354e-01          NaN
    ## [24071]          NaN          NaN          NaN 5.564562e-01          NaN
    ## [24076]          NaN          NaN 3.806877e-01          NaN 9.402364e-12
    ## [24081]          NaN 2.453655e-01          NaN 3.204838e-01 1.581604e-01
    ## [24086]          NaN          NaN 1.188675e-03          NaN 7.716360e-01
    ## [24091] 2.513050e-01          NaN 8.390318e-01 2.519094e-01          NaN
    ## [24096]          NaN          NaN          NaN          NaN 1.812173e-10
    ## [24101]          NaN 7.722000e-01          NaN 3.204838e-01 9.984977e-01
    ## [24106] 7.182834e-01 6.504810e-01 8.309490e-02          NaN 1.408528e-04
    ## [24111] 4.067812e-01          NaN          NaN 6.961685e-01          NaN
    ## [24116]          NaN          NaN 5.320840e-01          NaN          NaN
    ## [24121]          NaN          NaN 1.855826e-03          NaN 3.204838e-01
    ## [24126]          NaN          NaN 2.337728e-05          NaN          NaN
    ## [24131] 9.948242e-01          NaN          NaN 7.607842e-02          NaN
    ## [24136] 7.589581e-02 4.068881e-04          NaN          NaN 6.838953e-01
    ## [24141] 8.959638e-01          NaN 1.835939e-01 4.889662e-03 3.204838e-01
    ## [24146]          NaN          NaN          NaN          NaN 6.463456e-01
    ## [24151]          NaN          NaN 2.586102e-02          NaN          NaN
    ## [24156]          NaN          NaN          NaN 1.790632e-01          NaN
    ## [24161] 3.204838e-01          NaN          NaN          NaN 3.204838e-01
    ## [24166]          NaN          NaN          NaN 5.867464e-01          NaN
    ## [24171] 3.862594e-01 3.204838e-01          NaN          NaN          NaN
    ## [24176]          NaN          NaN          NaN          NaN          NaN
    ## [24181] 9.311736e-01          NaN          NaN          NaN 8.298863e-03
    ## [24186]          NaN          NaN 3.204838e-01 3.204838e-01          NaN
    ## [24191]          NaN 6.470294e-47          NaN          NaN          NaN
    ## [24196] 6.837977e-01          NaN 1.260685e-03 4.793666e-01          NaN
    ## [24201]          NaN 3.830158e-01          NaN          NaN          NaN
    ## [24206] 3.087594e-01          NaN          NaN          NaN          NaN
    ## [24211]          NaN 2.211050e-01          NaN 3.204838e-01 3.204838e-01
    ## [24216] 5.266760e-04 3.204838e-01          NaN 7.514309e-02 6.899673e-07
    ## [24221] 8.457220e-01 1.466199e-09 1.760980e-07          NaN          NaN
    ## [24226] 4.331922e-01          NaN          NaN 8.391606e-02          NaN
    ## [24231]          NaN          NaN          NaN          NaN 1.382463e-08
    ## [24236]          NaN          NaN          NaN          NaN          NaN
    ## [24241]          NaN          NaN 5.263197e-01 6.475325e-01          NaN
    ## [24246]          NaN          NaN 1.495096e-01          NaN          NaN
    ## [24251] 1.838032e-03 1.475304e-02          NaN 9.276227e-01          NaN
    ## [24256] 1.559125e-01          NaN 4.101147e-01          NaN 8.620420e-01
    ## [24261] 3.961354e-01 3.204838e-01          NaN 9.680031e-01          NaN
    ## [24266]          NaN 2.441175e-15 6.280252e-22          NaN          NaN
    ## [24271] 9.310493e-01          NaN 3.204838e-01 1.405217e-04 8.513289e-01
    ## [24276]          NaN 4.925049e-02 1.528208e-06 7.527805e-13 4.811999e-57
    ## [24281] 9.306922e-02 8.538927e-01 1.005702e-02 5.344078e-02          NaN
    ## [24286] 6.617305e-01 1.194200e-06 2.958616e-01 6.905641e-01 1.779360e-05
    ## [24291]          NaN 6.933957e-01 5.622309e-01          NaN 9.711464e-01
    ## [24296]          NaN 3.204838e-01 5.000394e-03 9.485650e-01 1.047254e-02
    ## [24301] 9.311871e-01 7.397166e-01          NaN 6.364543e-01 8.931256e-01
    ## [24306] 5.770583e-01          NaN 1.265359e-01          NaN          NaN
    ## [24311]          NaN 9.280088e-01 5.569007e-01 7.927403e-01 6.480404e-03
    ## [24316] 2.714194e-01          NaN 1.715335e-23 8.337119e-01 9.250071e-01
    ## [24321]          NaN 6.313336e-01 2.477665e-01          NaN          NaN
    ## [24326] 5.217837e-01 8.420856e-02 1.157407e-01          NaN 4.095309e-01
    ## [24331] 3.204838e-01 1.489674e-09 6.564745e-01          NaN 8.512029e-02
    ## [24336] 5.753890e-01 1.466057e-01          NaN 7.889396e-01          NaN
    ## [24341]          NaN 5.737403e-01 6.574245e-01 5.120757e-01 7.132679e-01
    ## [24346] 3.204838e-01          NaN 5.009663e-06          NaN 1.937357e-24
    ## [24351] 9.191210e-01 1.673940e-01          NaN 7.773176e-02 1.559566e-01
    ## [24356]          NaN 5.301175e-03 5.712896e-01 9.958863e-28 3.714146e-19
    ## [24361] 5.326147e-01 3.675097e-03          NaN 2.297654e-04          NaN
    ## [24366] 9.323086e-01 9.218676e-03 2.366749e-02 2.965194e-01 3.204838e-01
    ## [24371] 9.608386e-01          NaN          NaN 7.860814e-01          NaN
    ## [24376]          NaN 9.708462e-01 2.390526e-10 3.989451e-01 1.994936e-01
    ## [24381] 2.383908e-02 3.204838e-01 8.551838e-01 3.037125e-02 5.996047e-01
    ## [24386] 1.589102e-10 4.126175e-01

``` r
p <- tt$p.value
qvals <- qvalue(tt$p.value)$qvalue
index <- which(qvals<=0.05)
abline(h=-log10(max(tt$p.value[index])))
```

![](20200325_circadian-batch-investigation-rna-seq_steep_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
load(file = paste0(WD,"/data/GSE5859Subset.rda"))
library(rafalib)
library(RColorBrewer)
library(genefilter)

load(file = paste0(WD,"/data/GSE5859.rda"))








pcaData <- DESeq2::plotPCA(rld.sub, intgroup=c("animal.registration.sex","Seq_batch"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=animal.registration.sex, shape=Seq_batch)) +
        geom_point(size=3) +
        #geom_text(aes(label=Tissue),hjust=0, vjust=0) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
        #coord_fixed() +
        ggtitle("Gastrocnemius Samples Sequenced Across Batches")
```

![](20200325_circadian-batch-investigation-rna-seq_steep_files/figure-gfm/unnamed-chunk-7-2.png)<!-- -->
