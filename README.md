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
* c - integer - constant for durbin trasform (0 is ln) 
* excludeName - character, names of columns to be excluded from analysis   
* na.rm - Logical - TRUE converts NA to 0, INF to positive max, -INF to negative max
* avgRep - logical, returns the mean intensity for each group 
* sdRep - logical, returns the standard deviation for each group note, avgRep and sdRep cannont be both set to T
* inverse - logical, inverts the incoming data (1/data)
* pThresh - logical - excludes data above the pValue threshold
* returnMeta - Logical, if T, will return all input meta data note: this cannot be used with column/p clustering or thresholding
* rmZero - logical - removes rows with all zero intensities

## METPlot ##

* This function takes a MAVEN pick picked format list and makes several types of plots
* inDat - METen format table
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
* c - integer - constant for durbin trasform (0 is ln) 
* excludeName - character, names of columns to be excluded from analysis   
* n - number of levels used in heatmap
* zlim - numeric vector of length 2 for the c(min, max) plotted value, if zlim is missing or NULL, then all data are plotted
* logScale - logical - TRUE uses log2 intervals and disable n
* byClass - Logical - TRUE agregates data according to values in the 'notes' column 
* excludeName - character, names of columns to be excluded from analysis
* na.rm - Logical - TRUE converts NA to 0, INF to positive max, -INF to negative max
* returnPCA - logical - TRUE returnes PCA scores
* zJitter - Logical, if true will offset the x scale for zPlots
* returns the dataframe used to plot the data


# METCluster
* Function for clustering untargeted mz data
* mzTable - dataframe - METen format dataframe from untargeted csv METen output
* mzThresh - numeric - Mass tolerance in ppm
* rtThresh - numeric - Renention time tolerance in minutes
* pW - numeric vector, 0 to 1 - weight given to mass,retention time, and correlation
* pThresh - numeric, 0 to 1 - overall probablility threshold
* adList - numeric vector - list of known adducts
* collapse - logical - T collapses groups to highest intensity signal
* clusterOnly - logical - T returns collapsed groups with multiple signals 
