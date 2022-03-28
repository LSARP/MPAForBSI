# Metabolic Preference Assay for Rapid Diagnosis of Bloodstream Infections #
Authors: Thomas Rydzak, Ryan A. Groves, Ruichuan Zhang, Raied Aburashed, Rajnigandha Pushpker, Maryam Mapar, and Ian A. Lewis


MET is an R statistics-based utility tool that supports untargeted metabolomics analysis, visualization, and statistical analyses of LC-MS datasets. MET reads .csv format lists, in a format compatible with the MAVEN LC-MS data analysis software tool [ref].

## METRead ##
  Read in a Maven peak picked CSV file 
## METWrite ##
  Write out a CSV file of the dataset
## METStats ##
* This function takes a MAVEN pick picked format list and performances various statistical techinques 
* inDat - METen format table
* rCst - Logical - TRUE clusters by row
* cCst - Logical - TRUE clusters by column
* pCst - Logical - TRUE ranks rows by P value 
* nCst - Logical - TRUE sorts samples by name
* grid - Logical - TRUE plots the a sample grid over the heatmap
* pca - Logical - TRUE makes the PCA plots
* linePlot - Logical - TRUE makes a line plot
* stackPlot - Logical - TRUE makes a staked line plot
* heatMap - Logical - TRUE makes a classic heatmap
* barPlot - Logical - Plots bar figure, note this needs replicates see "by"
* zPlot - Logical - this plots data relative to a reference sample set
* scale - controls how data are normalized must be one of the following: "log", "durbin", "row", "column", or "none"
   - 'log' applies base 2 log transform
   - 'durbin' applies durbin transform
   - 'row' normalizes by row
   - 'column' noralizes by column
   - 'zScore' normalizes relative to a control, sets equal variance
   - 'none' does not apply any normalization

 
