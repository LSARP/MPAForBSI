##################################################################
##################################################################
# Part 1: untargeted Analysis 
##################################################################
##################################################################

##################################################################
# Initialization and read input file
##################################################################

#Set your working directory - your files will be saved here too
setwd(choose.dir())

#Load in the source file MET_LIBRARY.R 
#contains in-house functions
source(file.choose())

# load in the input file Raw_Input_File_S1.csv
dat <- METRead()


##################################################################
# Filtering biomarkers based on p-val, intensity, and fold change
##################################################################
# Run ANOVA on all markers and filter insignificant
# markers w/ bonferroni correction
alpha_val <- 0.05
pDat <- METStats(dat, pCalc = T)                  # Total # of markers = 4362
num_samples = nrow(pDat)
idx_p <- pDat$pVal < alpha_val/(num_samples)     # bonferroni Correction and filter
pDat <- pDat[idx_p,]

# Remove NA data from frame
pDat <- na.omit(pDat)                            # 1864 markers

# Write 
METWrite(pDat)

# Calculate the average threshold per species 
mDat <- METStats(pDat, avgRep = T)      
# Plot average threshold heatmap
METPlot(mDat, heat = T, scale = 'row', rCst = T, grid = F, cCst = F)

# Remove markers if the max average across all groups < 20000
thresh_num <- 20000
idx_num <- apply(abs(metaSep(mDat)$data), 1, max) > thresh_num
length(idx_num[idx_num== TRUE])                  # removed 296 markers 

# filter markers from data set if the threshold is less than 
# 4-fold change in comparison to MHB
thresh_fold <- 4
# select MHB or 1 as reference when calling this function 
fDat <- METStats(mDat, scale = 'fold')        
idx_fold <- (apply(abs(metaSep(fDat)$data), 1, max) > thresh_fold)

## Find data with correct intensity, pValue and fold change (533 markers)
idx <- which((idx_num + idx_fold) == 2)
out <- pDat[idx,]

METWrite(out)

##################################################################
# Clustering Markers
##################################################################
# resolution of masses to delta 5ppm between mass range of 100-500 m/z
mz_ad <- c(0,
           0.0005,        #(e-)
           1.0078,        #(H+)         https://www.chem.ualberta.ca/~massspec/atomic_mass_abund.pdf
           1.0062,        #(H neutron), https://www.chem.ualberta.ca/~massspec/atomic_mass_abund.pdf
           1.0034,        #(C neutron), https://www.chem.ualberta.ca/~massspec/atomic_mass_abund.pdf
           0.9694,        #(N neutron), https://www.chem.ualberta.ca/~massspec/atomic_mass_abund.pdf
           18.0106,       #(H2O)        https://hmdb.ca/metabolites/HMDB0002111
           34.9689,       #(Cl)         https://www.chem.ualberta.ca/~massspec/atomic_mass_abund.pdf
           38.9637,       # K           https://www.chem.ualberta.ca/~massspec/atomic_mass_abund.pdf
           22.9898,       #(Na)         https://www.chem.ualberta.ca/~massspec/atomic_mass_abund.pdf
           44.9977,       #(Formate),   https://hmdb.ca/metabolites/HMDB0304356
           1.0042,        #(O neutrons) https://www.chem.ualberta.ca/~massspec/atomic_mass_abund.pdf
           0.9994,        #(S neutrons) https://www.chem.ualberta.ca/~massspec/atomic_mass_abund.pdf
           1.9958,        #(S neutrons) https://www.chem.ualberta.ca/~massspec/atomic_mass_abund.pdf
           3.9950)        #(S neutrons) https://www.chem.ualberta.ca/~massspec/atomic_mass_abund.pdf


cDat <- METCluster(out, pThresh = 0.65, collapse = FALSE, clusterOnly = FALSE, adList = mz_ad, pW = c(1,1,1))
cDatCleaned <- cDat[!duplicated(cDat$grp),]
METWrite(cDatCleaned)				#210 markers
METPlot(cDatCleaned, heat = T, scale = 'row', rCst = T, grid = F, cCst = F)

##################################################################
##################################################################
# Part 2: Validation of Biomarkers
##################################################################
##################################################################

# Read file
# load in the input file Raw_Input_File_S2.csv
dat <- METRead()
# store the names of each col
markers_name <- names(dat)
# Remove the column not associated with markers
markers_name <- markers_name[6:length(markers_name)]

################################################################
# Statistics: 
################################################################
confi_inteval = 0.95

dfANOVA <- data.frame()
dfANOVA_Species <- data.frame()
dfTukey <- data.frame()

################################################################
# One way ANOVA with Post ADHOC TUKEY
################################################################
for(marker in markers_name){
  
  ################################################################
  # ANOVA Species with MHB
  ################################################################
  aov_species <- aov(dat[,marker] ~ Species,data=dat)
  aov_species <- unclass(summary(aov_species))
  
  df <- data.frame(aov_species)
  df['terms']=row.names(df)
  df <- cbind(markerID=paste(gsub("X", "Marker_", marker)),df)
  row.names(df)=NULL
  
  dfANOVA_Species <- rbind(dfANOVA_Species, df)
  
  ################################################################
  # ANOVA 
  ################################################################
  
  # Some of the data has no Sex data removing data with NA
  dat_no_NA <- na.omit(dat)
  
  ANOVA <- aov(dat_no_NA[,marker] ~ Species*Sex*AgeGroup,data=dat_no_NA)
  ANOVA<-unclass(summary(ANOVA))
  
  df <- data.frame(ANOVA)
  df['terms']=row.names(df)
  df <- cbind(markerID=paste(gsub("X", "Marker_", marker)),df)
  row.names(df)=NULL
  
  dfANOVA <- rbind(dfANOVA, df)
  
  ################################################################
  # Post ADHOC TUKEY
  ################################################################
  
  ############################ Species ############################
  oneWay<- aov(dat[,marker] ~ Species, data = dat)
  Tukey <- TukeyHSD(oneWay,conf.level = confi_inteval)
  Tukey <- data.frame(Tukey$Species)
  Tukey<-cbind(row.names(Tukey),Tukey)
  row.names(Tukey) <- NULL
  Tukey <- cbind(markerID=paste(gsub("X", "Marker_", marker)),Tukey)
  dfTukey <- rbind(dfTukey,Tukey)
  
  ############################## Sex ##############################
  # Some of the data has no Sex data removing data with NA
  dat_no_NA <- na.omit(dat)
  oneWay<- aov(dat_no_NA[,marker] ~ Sex, data = dat_no_NA)
  Tukey <- TukeyHSD(oneWay,conf.level = confi_inteval)
  Tukey <- data.frame(Tukey$Sex)
  Tukey<-cbind(row.names(Tukey),Tukey)
  row.names(Tukey) <- NULL
  Tukey <- cbind(markerID=paste(gsub("X", "Marker_", marker)),Tukey)
  dfTukey <- rbind(dfTukey,Tukey)
  
  ############################## Age ##############################
  oneWay<- aov(dat_no_NA[,marker] ~ AgeGroup, data = dat_no_NA)
  Tukey <- TukeyHSD(oneWay,conf.level = confi_inteval)
  Tukey <- data.frame(Tukey$AgeGroup)
  Tukey<-cbind(row.names(Tukey),Tukey)
  row.names(Tukey) <- NULL
  Tukey <- cbind(markerID=paste(gsub("X", "Marker_", marker)),Tukey)
  dfTukey <- rbind(dfTukey,Tukey)
}

# Store only values without residuals
dfANOVA_Species <- dfANOVA_Species[ which(dfANOVA_Species$terms=='Species    '), ]
dfANOVA <- dfANOVA[which(dfANOVA$terms!='Residuals           '), ]


# Write data to CSV
METWrite(dfANOVA_Species)			#Just marker-species associations
METWrite(dfTukey)				#Tukey-Kramer for marker-species
METWrite(dfANOVA)				#Age and sex linked associations (no MHB)
