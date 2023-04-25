# Geographically weighted accuracy for hard and soft land cover classifications: 4 approaches with coded illustrations

Alexis Comber<sup>1*</sup> and Naru Tsutsumida<sup>2</sup> 


<sup>1</sup> School of Geography, University of Leeds, Leeds, UK.\
<sup>2</sup> Department of Information and Computer Sciences, Saitama University, Japan\
<sup>*</sup> contact author: a.comber@leeds.ac.uk

## Abstract
This paper proposes a novel spatially varying coefficient (SVC) regression through a Geographical Gaussian Process GAM (GGP-GAM): a Generalized Additive Model (GAM) with Gaussian Process (GP) splines parameterised at observation locations. As with Multiscale Geographically Weighted Regression (MGWR), the proposed GGP-GAM is multiscale but has fewer theoretical and technical limitations. A GGP-GAM was applied to simulated coefficient datasets exhibiting varying degrees of spatial heterogeneity and was shown to out-perform MGWR under a range of fit metrics. The proposed SVC model was then applied to a Brexit vote case study with covariates describing a range of socio-economic, life-stage and voting factors. Resultant coefficient estimates were mapped, showing regional scales of relationship non-stationarity in the UK's Brexit voting process. For the proposed GGP-GAM, a number of areas of further work are identified includeing the creation of a user-friendly package to support model creation and coefficient mapping, and to facilitate ease of comparison with alternate SVC models. Calibration issues are discussed such as GAM tuning (particularly of knots and spline smoothing parameters) and the difficulty in linking the GGP-GAM spline smoothing parameters to intuitive user understandings of process spatial heterogeneity, where corresponding parameters for MGWR (i.e., its kernel bandwidths) hold an advantage.

This paper will be submitted to IJRS.

## Code
To run the analysis in this paper you should download the the R script `gw_acc.R`, the data files (`validation_data.csv`) and install the packages. Package and other info is below. The data files and supporting scripts will need will need to be locally available . The code recreates the results as the same sequence in the paper. 

If you have any problems with data / code / versions etc please contact Lex Comber at the email above.
```{r}
sessionInfo()
R version 4.2.0 (2022-04-22)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Monterey 12.6.1

Matrix products: default
BLAS:   /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRblas.0.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib

locale:
[1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] robCompositions_2.3.1 data.table_1.14.2     pls_2.8-1             gwxtab_0.1.0          scales_1.2.1         
 [6] cols4all_0.6          cowplot_1.1.1         ggspatial_1.1.6       spgwr_0.6-35          spData_2.2.0         
[11] forcats_0.5.2         stringr_1.5.0         dplyr_1.1.0           purrr_1.0.1           readr_2.1.2          
[16] tidyr_1.3.0           tibble_3.1.8          ggplot2_3.4.1         tidyverse_1.3.2       sp_1.6-0             
[21] sf_1.0-9             

loaded via a namespace (and not attached):
  [1] googledrive_2.0.0     colorspace_2.1-0      modeltools_0.2-23     ellipsis_0.3.2        class_7.3-20         
  [6] rgdal_1.5-32          mclust_5.4.10         fs_1.6.1              proxy_0.4-27          farver_2.1.1         
 [11] cvTools_0.3.2         flexmix_2.3-19        fansi_1.0.4           mvtnorm_1.1-3         lubridate_1.8.0      
 [16] ranger_0.14.1         xml2_1.3.3            codetools_0.2-18      splines_4.2.0         robustbase_0.95-0    
 [21] knitr_1.42            jsonlite_1.8.4        broom_1.0.1           kernlab_0.9-31        cluster_2.1.4        
 [26] dbplyr_2.2.1          png_0.1-8             rrcov_1.7-2           compiler_4.2.0        httr_1.4.4           
 [31] rainbow_3.7           backports_1.4.1       assertthat_0.2.1      Matrix_1.5-1          gargle_1.2.1         
 [36] cli_3.6.0             tools_4.2.0           gtable_0.3.1          glue_1.6.2            RANN_2.6.1           
 [41] Rcpp_1.0.10           carData_3.0-5         cellranger_1.1.0      raster_3.6-14         zCompositions_1.4.0-1
 [46] vctrs_0.5.2           fpc_2.2-10            robustHD_0.7.4        lmtest_0.9-40         xfun_0.37            
 [51] laeken_0.5.2          rvest_1.0.3           lifecycle_1.0.3       googlesheets4_1.0.1   terra_1.7-3          
 [56] DEoptimR_1.0-11       MASS_7.3-58.1         zoo_1.8-11            VIM_6.2.2             hms_1.1.2            
 [61] parallel_4.2.0        RColorBrewer_1.1-3    NADA_1.6-1.1          gridExtra_2.3         reshape_0.8.9        
 [66] stringi_1.7.12        pcaPP_2.0-2           e1071_1.7-13          boot_1.3-28           truncnorm_1.0-9      
 [71] hdrcde_3.4            prabclus_2.3-2        rlang_1.0.6           pkgconfig_2.0.3       bitops_1.0-7         
 [76] pracma_2.4.2          fda_6.0.5             lattice_0.20-45       labeling_0.4.2        ks_1.14.0            
 [81] tidyselect_1.2.0      deSolve_1.35          GGally_2.1.2          plyr_1.8.7            magrittr_2.0.3       
 [86] R6_2.5.1              generics_0.1.3        DBI_1.1.3             pillar_1.8.1          haven_2.5.1          
 [91] withr_2.5.0           units_0.8-1           survival_3.4-0        abind_1.4-5           RCurl_1.98-1.8       
 [96] nnet_7.3-17           perry_0.3.1           ggfortify_0.4.16      modelr_0.1.9          crayon_1.5.2         
[101] car_3.1-0             KernSmooth_2.23-20    utf8_1.2.3            tzdb_0.3.0            fds_1.8              
[106] grid_4.2.0            readxl_1.4.1          diptest_0.76-0        vcd_1.4-11            reprex_2.0.2         
[111] classInt_0.4-8        stats4_4.2.0          munsell_0.5.0        
```
