################################################################################
################################################################################
##                                                                            ##
##                                                                            ##
##    MET version 1.0.0, Tools for processing and visualization of            ##
##    mass spectroscpy files                                                  ##
##    Copyright (C) 2022 Ian A. Lewis                                         ##
##                                                                            ##
##    This program is free software: you can redistribute it and/or modify    ##
##    it under the terms of the GNU General Public License as published by    ##
##    the Free Software Foundation, either version 3 of the License, or       ##
##    any later version.                                                      ##
##                                                                            ##
##    This program is distributed in the hope that it will be useful,         ##
##    but WITHOUT ANY WARRANTY; without even the implied warranty of          ##
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           ## 
##    GNU General Public License for more details.                            ##
##                                                                            ##
##    A copy of the GNU General Public License can be found at:               ##
##    www.r-project.org/Licenses/GPL-3                                        ##
##                                                                            ##
##                                                                            ##
################################################################################
################################################################################

##################
# Functions List
##################
# METRead
# METWrite
# METStats
# METPlot
# METCluster
##################

################################################################################
##                                                                            ##
##                  Simplifed Reading and writing files                       ##
##                                                                            ##
################################################################################

# Read Files
METRead <- function( inPath = file.choose()){
	return(read.csv( inPath, head = TRUE, stringsAsFactors = FALSE))
} 

# Write Files 
METWrite <- function( table ){
	write.csv( table, file.choose(new=TRUE), quote = FALSE, row.names = FALSE )
}
################################################################################

################################################################################
##                                                                            ##
##                           Internal functions                               ##
##                                                                            ##
################################################################################

myNaN <- function( x ){
	if( any(is.nan(x)) )
		x[ is.nan(x) ] <- 0
	return(x)
} 
## Function for dividing by zero
## x - dataframe with numerator data
## y - dataframe with denominator data
## num - numeric, the value to be assigned for non-zero numerators dived by 0
## den - numeric, the value to be assigned for 0/0
## returns dataframe of dim(x)
myDiv <- function(x, y, num = 0, den = 1){
	dX <- dim(x)
	dY <- dim(y)
	if (all(dX != dY))
		stop("x and y must have the same number of rows and columns")
	vX <- as.numeric(as.vector(unlist(x)))
	vY <- as.numeric(as.vector(unlist(y)))
	if(any(vY == 0) )
		vY[vY == 0 ] <- den
	out <- data.frame(matrix(vX/vY, ncol = ncol(x)))
	names(out) <- names(x)
	return(out)
}


findPval <- function( sLab, mzData, type = 'p'  ){
	nr <- nrow( mzData )
	pVal <- list()
	sLab <- factor(sLab)
	
	for( i in 1:nr ){	
		
		tVal <- tryCatch( pVal[[i]] <- aov(as.vector(unlist(mzData[ i, ])) ~ sLab), 
				error=function(er) return(NA))
		
		if( is.na(tVal[1]) ){
			pVal[[i]] <- NA
			next
		}
		
		if(type == 'p')
			pVal[[i]] <- unclass(summary.aov(tVal))[[1]]$Pr[1]
		
		if(type =='f')
			pVal[[i]] <- unclass(summary.aov(tVal))[[1]]$"F value"[1]

		if(type == 'df')
		pVal[[i]] <- unclass(summary.aov(tVal))[[1]]$"rank"

	}
	return(as.vector(unlist(pVal)))
}

## Internal function heatmap function adapted from base heatmap function in r 
## minor modification to handle errors
## see: https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/heatmap
## for documentation
myheatmap <- function (x, Rowv = NULL, Colv = if (symm) "Rowv" else NULL, 
		reorderfun = function(d,w) reorder(d, w), add.expr, 
		symm = FALSE, revC = identical(Colv, "Rowv"), na.rm = TRUE, 
		grid = TRUE, breaks = NULL, margins = c(5, 5), 
		ColSideColors, RowSideColors, cexRow = 0.2 + 
				1/log10(nr), cexCol = 0.2 + 1/log10(nc), labRow = NULL, 
		labCol = NULL, main = NULL, xlab = NULL, ylab = NULL, keep.dendro = TRUE, 
		verbose = getOption("verbose"), ...) {
	if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
		stop("'x' must be a numeric matrix")
	
	x[ x == Inf ] <- max(x[is.finite(unlist(x))], na.rm = T) * 100 	
	x[ x == -Inf ] <- min(x[is.finite(unlist(x))], na.rm = T) / 100 
	
	nr <- di[1L]
	nc <- di[2L]
	if (nr <= 1 || nc <= 1) 
		stop("'x' must have at least 2 rows and 2 columns")
	if (!is.numeric(margins) || length(margins) != 2L) 
		stop("'margins' must be a numeric vector of length 2")
	doRdend <- !identical(Rowv, NA)
	doCdend <- !identical(Colv, NA)
	if (!doRdend && identical(Colv, "Rowv")) 
		doCdend <- FALSE
	if (is.null(Rowv)) 
		Rowv <- rowMeans(x, na.rm = na.rm)
	if (is.null(Colv)) 
		Colv <- colMeans(x, na.rm = na.rm)
	if (doRdend) {
		if (inherits(Rowv, "dendrogram")) 
			ddr <- Rowv
		else {
			hcr <- hclust(dist(x))
			ddr <- as.dendrogram(hcr)
			if (!is.logical(Rowv) || Rowv) 
				ddr <- reorderfun(ddr, Rowv)
		}
		if (nr != length(rowInd <- order.dendrogram(ddr))) 
			stop("row dendrogram ordering gave index of wrong length")
	}
	else rowInd <- 1L:nr
	if (doCdend) {
		if (inherits(Colv, "dendrogram")) 
			ddc <- Colv
		else if (identical(Colv, "Rowv")) {
			if (nr != nc) 
				stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
			ddc <- ddr
		}
		else {
			hcc <- hclust(dist(if (symm) 
										x
									else t(x)))
			ddc <- as.dendrogram(hcc)
			if (!is.logical(Colv) || Colv) 
				ddc <- reorderfun(ddc, Colv)
		}
		if (nc != length(colInd <- order.dendrogram(ddc))) 
			stop("column dendrogram ordering gave index of wrong length")
	}
	else colInd <- 1L:nc
	x <- x[rowInd, colInd]
	labRow <- if (is.null(labRow)) 
				if (is.null(rownames(x))) 
					(1L:nr)[rowInd]
				else rownames(x)
			else labRow[rowInd]
	labCol <- if (is.null(labCol)) 
				if (is.null(colnames(x))) 
					(1L:nc)[colInd]
				else colnames(x)
			else labCol[colInd]
	lmat <- rbind(c(NA, 3), 2:1)
	lwid <- c(if (doRdend) 1 else 0.05, 4)
	lhei <- c((if (doCdend) 1 else 0.05) + if (!is.null(main)) 0.2 else 0, 
			4)
	if (!missing(ColSideColors)) {
		if (!is.character(ColSideColors) || length(ColSideColors) != 
				nc) 
			stop("'ColSideColors' must be a character vector of length ncol(x)")
		lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
		lhei <- c(lhei[1L], 0.2, lhei[2L])
	}
	if (!missing(RowSideColors)) {
		if (!is.character(RowSideColors) || length(RowSideColors) != 
				nr) 
			stop("'RowSideColors' must be a character vector of length nrow(x)")
		lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 
						1), lmat[, 2] + 1)
		lwid <- c(lwid[1L], 0.2, lwid[2L])
	}
	lmat[is.na(lmat)] <- 0
	if (verbose) {
		cat("layout: widths = ", lwid, ", heights = ", lhei, 
				"; lmat=\n")
		print(lmat)
	}
	op <- par(no.readonly = TRUE)
	on.exit(par(op))
	
	layout(lmat, widths = lwid, heights = lhei, respect = TRUE)
	if (!missing(RowSideColors)) {
		par(mar = c(margins[1L], 0, 0, 0.5))
		image(rbind(1L:nr), col = RowSideColors[rowInd], axes = FALSE)
	}
	if (!missing(ColSideColors)) {
		par(mar = c(0.5, 0, 0, margins[2L]))
		image(cbind(1L:nc), col = ColSideColors[colInd], axes = FALSE)
	}
	par(mar = c(margins[1L], 0, 0, margins[2L]))
	x <- t(x)
	if (revC) {
		iy <- nr:1
		if (doRdend) 
			ddr <- rev(ddr)
		x <- x[, iy]
	}
	else iy <- 1L:nr
	
	
	if( !is.null(breaks) && length(breaks) > 1 ){
		x[ x > max(breaks) ] <- max(breaks)
		x[ x < min(breaks) ] <- min(breaks)
	}
	
	image(1L:nc, 1L:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
					c(0, nr), axes = FALSE, xlab = "", ylab = "", breaks = breaks, ...)
	axis(1, 1L:nc, labels = labCol, las = 2, line = -0.5, tick = 0, 
			cex.axis = cexCol, col.axis = par('fg'))
	box()
	if( grid )
		abline(v = 1L:nc + .5, h = iy + .5)
	if (!is.null(xlab)) 
		mtext(xlab, side = 1, line = margins[1L] - 1.25)
	axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0, 
			cex.axis = cexRow, col.axis = par('fg'))
	if (!is.null(ylab)) 
		mtext(ylab, side = 4, line = margins[2L] - 1.25)
	if (!missing(add.expr)) 
		eval(substitute(add.expr))
	par(mar = c(margins[1L], 0, 0, 0))
	if (doRdend) 
		plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
	else frame()
	par(mar = c(0, 0, if (!is.null(main)) 1 else 0, margins[2L]))
	if (doCdend) 
		plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
	else if (!is.null(main)) 
		frame()
	if (!is.null(main)) {
		par(xpd = NA)
		title(main, cex.main = 1.5 * op[["cex.main"]])
	}
	invisible(list(rowInd = rowInd, colInd = colInd, Rowv = if (keep.dendro && 
							doRdend) ddr, Colv = if (keep.dendro && doCdend) ddc))
}

foldChange <- function( inTable, normVec, zeroOffset = F, normMin = .01,scale = 'fold' ){
	## Calculate standard deviation for normalization vector
	if( scale == "zScore"){
		if( length(normVec) < 1  )
			stop('Z scales require a standard deviation from multiple replicates')
		sdTab <- apply(inTable[, normVec], 1, mySd )
		normTab <- apply(inTable[, normVec], 1, myMean )
	}
    else{
		if( length(normVec) > 1 )
			normTab <- apply(inTable[, normVec], 1, myMean)
		else
			normTab <- inTable[, normVec ]		
	}
	
	## Make norm matrix
	normTab <- as.vector(unlist(normTab))
	normTab <- data.frame(matrix(rep(normTab, ncol(inTable)), 
					ncol = ncol(inTable)))
	
	## Make sd matrix
	if( scale == "zScore"){
		sdTab <- as.vector(unlist(sdTab))
		sdTab <- data.frame(matrix(rep(sdTab, ncol(inTable)),ncol = ncol(inTable)))
    }

	## Set all 0 values to 1% of min observed value
	if( zeroOffset ){
		zSet <- inTable == 0 
		if( mySum(zSet) ){
			tTab <- inTable 
			tTab[ zSet ] <- NA
			mSet <- as.vector(unlist(apply(tTab, 1, myMin )))
			for( i in 1:length(mSet) ){
				if( any(inTable[i,] == 0, na.rm =T) ){
					idx <- which( inTable[i,] == 0 )
					inTable[i, idx] <- mSet[i]				
				}
			}
		}
	}



	## Calculate z score
	if( scale == 'zScore' ){
		out <- (inTable - normTab)/sdTab 
		#out <- apply( out, 2, myNaN )		
		return(out)		
	}
	

	## Normalize	
	out <- myDiv(inTable,  normTab)

	if( scale == 'ratio' ){
		out <- apply( out, 2, myNaN )		
		return(out)		
	}
	
	if( any(out < 1 ) ){
		foldFun <- function(x) { 
			tIdx <- which( x < 1)
			x[ tIdx ] <- abs(1/ x[tIdx]) * -1
			return(x)
		}
		out <- apply(out, 2, foldFun)	
		out <- apply( out, 2, myNaN )
	}
	
	return(out)
	
}


## Function for plotting metabolite heatmaps
metImage <- function( x, col, breaks, cCst = FALSE, rCst = FALSE, grid = TRUE ){
	## Generate levels and colors
	if( missing(col) || missing(breaks) ){
		lDat <- levelFun(x)
		col <- lDat$col
		breaks <- lDat$breaks
	}
	## Set col/row dendrogram
	Colv <- Rowv <- NA
	if( cCst )
		Colv <- NULL
	if( rCst )
		Rowv <- NULL
	windows(width = 2)
	w <- (3 + par('mar')[2L]) * par("csi") * 2.54
	layout(matrix(c(2, 1), ncol = 2L), widths = c(1, lcm(w)))
	par( mar = c( 5.1, 1, 4.1, 4.1), las = 1)
	plot.new()
	plot.window(xlim = c(0, 1), ylim = range(breaks, na.rm = T), xaxs = "i", yaxs = "i")
	rect(0, breaks[-length(breaks)], 1, breaks[-1L], col = col)
	axis(4, at = breaks)
	windows()
	out <- myheatmap(as.matrix(x), Colv = Colv, Rowv = Rowv, na.rm = T,
			col = col, breaks = breaks, grid = grid )
	return(invisible(out))	
}

levelFun <- function(x, n = 10, r = c(0, 1), g = c(0,0), b = c(1,0), zlim, 
		logScale = T, returnInterval = FALSE){
	
	x <- suppressWarnings(as.numeric(as.vector(unlist(x))))
	
	
	## Remove NA and Inf data
	if( any(is.na(x)) )
		x[is.na(x)] <- 0
	if( any(abs(x) == Inf))
		x[ abs(x) == Inf ] <- 0
	
	## Threshold data
	if( !missing(zlim) && length(zlim) == 2 && all(is.numeric(zlim)) ){
		zlim <- sort(abs(zlim) )
	}else
		zlim <- c(  min(abs(x)),  max(abs( x )) )
	
	## Set breaks
	if( logScale ){
		
		scale <- 2^(0:50)
		xSeq <- scale[which( scale >= zlim[1] )[1] : which( scale >= zlim[2] )[1]]
		n <- length(xSeq) 		
		xSeq <- unique(c( -rev(xSeq), xSeq))
	
	}else{
		xSeq <- (c( 
							-rev(seq( zlim[1], zlim[2], length.out = n)), 
							seq( zlim[1], zlim[2], length.out = n) ))								
	}
	
	
	## Set colors	 
	col <- c(
			rgb(rep(r[1],(n-1)),
					rep(g[1],(n-1)), 
					rep(b[1],(n-1)), alpha = rev((1:n)[-1] / n)), 
			'transparent',
		  rgb(rep(r[2],(n-1)),
					rep(g[2],(n-1)), 
					rep(b[2],(n-1)),  alpha = (1:n)[-1] / n) 		
	)		 
	
	if( !returnInterval )
		return( list(breaks = xSeq, col = col))
	return(list(breaks = xSeq, 
							col = col, 
							interval = findInterval(x, xSeq))) 
}


## Internal function for seperating METEN metadata from datasets
## excludeName - character, names of columns to be excluded from analysis
## METNam - character string - default names of METen format metadata
metaSep <- function(inDat, 
		METNam = c('label', 'metaGroupId', 'groupId', 'goodPeakCount', 'medMz', 'medRt', 
				'maxQuality', 'note', 'compound', 'compoundId', 'expectedRtDiff', 
				'ppmDiff', 'parent', 'adductName', 'isotopeLabel', 'formula', 'ctName', 'pVal','grp','nGrp','gPar'), excludeName){
	## Check input dataset
	if( missing(inDat) || !is.data.frame(inDat) || (nrow(inDat) < 1) )
		stop('A METEN format dataset is required')
	## Combine default and user-provided names
	if( !missing(excludeName) )
		METNam <- c(METNam, excludeName)	
	## Remove all standard METen meta data
	mIdx <- na.omit(match(METNam, names(inDat)))
	## Find metadaFta
	if( length(mIdx) == 1 ){
		metaDat <- data.frame(compound = inDat[,mIdx], stringsAsFactors = F)
	}else
		metaDat <- inDat[,mIdx]
	## remove the metaData
	outDat <- inDat[,-na.omit(match(METNam, names(inDat)))]
	return(list( metaData = metaDat, data = outDat))
}



################################################################################
##                                                                            ##
##                              USER FUNCTIONS                                ##
##                                                                            ##
################################################################################
## Function for clustering untargeted mz data
## mzTable - dataframe - METen format dataframe from untargeted csv METen output
## mzThresh - numeric - Mass tolerance in ppm
## rtThresh - numeric - Renention time tolerance in minutes
## pW - numeric vector, 0 to 1 - weight given to mass, 
##                      retention time, and correlation
## pThresh - numeric, 0 to 1 - overall probablility threshold
## adList - numeric vector - list of known adducts
## collapse - logical - T collapses groups to highest intensity signal
## clusterOnly - logical - T returns collapsed groups with multiple signals 
METStats <- function(inDat, rCst = F, cCst = F, pCst = F, nCst = F, 
		usrSelect = F, sdRep = F, by = "\\.", excludeName = NULL, 
		scale = c("log", "durbin", "row", "column", "fold", "ratio", "zScore", "none"), 
		inverse = F, c = 1, n = 10, maxDup = F, avgRep = F, 
		na.rm  = T, sLab, dThresh, pCalc = F, pThresh = 1, bonferroni = F,
		returnMeta = T, rmZero = F,  ...){
	
	## check format of peaklist
	if( !is.data.frame(inDat) || nrow(inDat) < 1 || ncol(inDat) < 2 
			|| !any(names(inDat) == 'compound'))
		stop("A dataframe with a 'compound' and at least one data column is required")
	
	## Set Scale argument
	scale <- if (missing(scale)) 
				"none"
			else match.arg(scale)
	
	## Check pThresh
	if( avgRep && pThresh < 1)
		stop("avgRep and pThresh can not be used in combination")
	if( pThresh < 1 && bonferroni )
		stop("pThresh and bonferroni can not be used in combination")
	if( avgRep && sdRep)
		stop("avgRep and sdRep can not be used in combination")
	
	## Make all compound names unique
	com <- inDat$compound
	if(any(is.na(com)) &&  any(names(inDat) == 'medMz') 
			&& any(names(inDat) == 'medRt')){
		com[which(is.na(com))] <- paste(inDat$medMz[which(is.na(com))], 
				inDat$medRt[which(is.na(com))], sep = "@")
	}
	
	for( i in unique(com) ){
		tCpd <- com[ com == i ]
		if( length(tCpd) == 1 )
			next
		com[com == i] <- paste(tCpd, 1:length(tCpd), sep = '_')
	}
	inDat$compound <- com
	
	
	## Remove duplicate metabolite records and find max for each
	if(maxDup)
		inDat <- findMaxRecord( inDat, excludeName = excludeName )
	
	## Combine any metadata fields to be excluded
	tDat <- metaSep(inDat, excludeName = excludeName)
	metaDat <- tDat[[1]]
	inDat <- tDat[[2]]
	
	## Order samples by name
	if (nCst){
		inDat <- inDat[,order(names(inDat))]
	}
	
	## hold sample names
	sam <- names(inDat) ## read tab header to get sample names
	rownames(inDat) <- metaDat$compound
	
	## Zero any data below the data threshold
	if( !missing(dThresh) )
		inDat[ inDat < dThresh ] <- 0
	
	## Remove all zero data
	if( any(!apply(inDat, 1, function(x) all(x == 0))) && rmZero ){
		inDat <- inDat[- which(apply(inDat, 1, function(x) all(x == 0))),]		
	}
	
	## Get sample labels
	if( (scale == 'zScore' || pCst || pCalc || pThresh < 1) && missing (sLab))
		sLab <- sapply( strsplit(sam, split = by), function (x) x[1])		
	
	
	## Average sample replicates
	if( avgRep || sdRep ){
		
		tLab <- list(factor(
						sapply(strsplit(names(inDat), split = by), function (x) x[1])))
		
		if( avgRep )
			inDat <- aggregate(t(inDat), FUN = function(x) mean(x, na.rm = T), by =tLab)
		else
			inDat <- aggregate(t(inDat), FUN = function(x) sd(x, na.rm = T), by =tLab)
		
		
		uName <- levels(inDat[,1]) 
		inDat <- data.frame(t(inDat[,-1]))
		names(inDat) <- uName
		
		## Reorder colums to match input list
		sam <- unique(sapply(strsplit(sam, split = by), function (x) x[1]))
		inDat <- inDat[,c(which(names(inDat) == 'compound'), 
						match(sam, names(inDat)))]
	}
	
	
	## Scale data
	if( scale == 'log'){
		inDat <- log( inDat )	
	}else if( scale == 'durbin' ){
		inDat <- durbin(inDat, c = c )	
	}else if (scale == "row") {
		inDat <- sweep(inDat, 1L, rowMeans(inDat, na.rm = TRUE), check.margin = FALSE)
		sx <- apply(inDat, 1L, sd, na.rm = TRUE)
		inDat <- sweep(inDat, 1L, sx, "/", check.margin = FALSE)
	}else if (scale == "column") {
		inDat <- sweep(inDat, 2L, colMeans(inDat, na.rm = TRUE), check.margin = FALSE)
		sx <- apply(inDat, 2L, sd, na.rm = TRUE)
		inDat <- sweep(inDat, 2L, sx, "/", check.margin = FALSE)
	}else if (scale == "fold" || scale == "ratio" || scale == 'zScore'){
		userSel <- select.list( names( inDat), multiple = T, 
				title = "Select reference:" )
		if( !length(userSel) || !nzchar(userSel) )
			stop('No file selected')
		normCol <- match(userSel, names(inDat))
		inDat <- foldChange( inTable = inDat, normVec = normCol, scale = scale)	
	}
	
	## Remove NA and Inf data
	if( na.rm  ){
		printList <- newList <- NULL
		
		if( any(is.na(inDat)) ){
			inDat[is.na(inDat)] <- 0
			printList <- c(printList, "NA")
			newList <- c(newList, 0)
		}
		
		if( any(inDat == Inf)){
			inDat[ inDat == Inf ] <- range(unlist(inDat), finite = T)[2]	
			newList <- c(printList, "Inf")
			newList <- c(newList, max(inDat, na.rm = T))
		}
		
		if( any(inDat == -Inf)){
			inDat[ inDat == -Inf ] <- range(unlist(inDat), finite = T)[1]	
			printList <- c(printList, "-Inf")
			newList <- c(newList, min(inDat, na.rm = T))
		}
		if(!is.null(printList))
			print( paste(c(paste(printList, collapse = ", "), 
									"data values were changed to:",	
									paste(newList, collapse = ", ")), collapse = " "))
	}	
	if( inverse )
		inDat <- 1/inDat
        
	## Rank values according to ANOVA p value
	if( pCst || pCalc || pThresh < 1){
		metaDat$pVal <- findPval(sLab = sLab, inDat, type = 'p')
		
		if( pCst ){
			sidx <- order(metaDat$pVal)
			inDat <- inDat[sidx,]
			if(ncol( metaDat) > 1 )
				metaDat <- metaDat[sidx,]
			else
				metaDat <- data.frame( compound = metaDat[sidx,], stringsAsFactors = F)
		}		
		## Calculate pThreshold
		if( bonferroni ){
			pThresh <- .05 / nrow(inDat) 	
		}
	
		if( pThresh < 1 ){
			idxP <- which(metaDat$pVal < pThresh)
			inDat <- inDat[idxP,]
			if(ncol( metaDat) > 1 )
				metaDat <- metaDat[idxP,]
			else
				metaDat <- data.frame( compound = metaDat[idxP,], stringsAsFactors = F)	
		}
	}
	
	## Allow users to select sub groups of data
	if( usrSelect ){
		userChoice <- mySelect(names(inDat))
		if( length(userChoice) > 0 ){
			inDat <- inDat[,userChoice]
			sLab <- sLab[userChoice]	
		}
	}

	## Add metadata
	if( returnMeta ){
		inDat <- cbind( metaDat, inDat)
	}else{
		inDat <- data.frame( compound = row.names(inDat), inDat, stringsAsFactors = F)
		names(inDat) <- c('compound', sam)
	}
	invisible(inDat)	
}


## Default MET plotter
## This function takes a METen-format list and makes several types of plots
## inDat - METen format table
## rCst - Logical - TRUE clusters by row
## cCst - Logical - TRUE clusters by column
## pCst - Logical - TRUE ranks rows by P value 
## nCst - Logical - TRUE sorts samples by name
## wildMan - Logical - Allows users to select a smaller set of data to plot
## grid - Logical - TRUE plots the a sample grid over the heatmap
## pca - Logical - TRUE makes the PCA plots
## linePlot - Logical - TRUE makes a line plot
## stackPlot - Logical - TRUE makes a staked line plot
## heatMap - Logical - TRUE makes a classic heatmap
## barPlot - Logical - Plots bar figure, note this needs replicates see "by"
## zPlot - Logical - this plots data relative to a reference sample set
## scale - controls how data are normalized
##			must be one of the following: "log", "durbin", "row", "column", or "none"
##					'log' applies base 2 log transform
##          'durbin' applies durbin transform
##          'row' normalizes by row
##          'column' noralizes by column
## 					'zScore' normalizes relative to a control, sets equal variance
##          'none' does not apply any normalization
## c - integer - constant for durbin trasform (0 is ln)
## n - number of levels used in heatmap
## zlim - numeric vector of length 2 for the c(min, max) plotted value
##      if zlim is missing or NULL, then all data are plotted
## logScale - logical - TRUE uses log2 intervals and disable n
## byClass - Logical - TRUE agregates data according to values in the 'notes' column 
## excludeName - character, names of columns to be excluded from analysis
## na.rm - Logical - TRUE converts NA to 0, INF to positive max, 
##                   -INF to negative max
## maxDup - logical, finds the maximum record in cases of duplicated compounds
## avgRep - logical, returns the mean intensity for each group
## sdRep - logical, returns the standard deviation for each group
## 				- note, avgRep and sdRep cannont be both set to T
## inverse - logical, inverts the incoming data (1/data)
## pThresh - logical - excludes data above the pValue threshold
## bonferroni - logical - calculates the bonferroni-corrected pThresh
## returnPCA - logical - TRUE returnes PCA scores
## returnMeta - Logical, if T, will return all input meta data
##             note: this cannot be used with column/p clustering or thresholding
## rmZero - logical - removes rows with all zero intensities
## random - Logical, TRUE will generate random number distribution that match
##          the mean/sd for each compound
## zJitter - Logical, if true will offset the x scale for zPlots
## returns the dataframe used to plot the data

METPlot <- function(inDat, rCst = F, cCst = F, pCst = F, nCst = F, 
		usrSelect = F,grid = FALSE, sdRep = F, by = "\\.", excludeName = NULL, 
		scale = c("log", "durbin", "row", "column", "fold", "ratio", "zScore", "none"), 
		c = 1, n = 10, zlim, logScale = F, maxDup = F, avgRep = F, 
		byClass = F, na.rm  = T, sLab, dThresh,
		pca = F, linePlot = F, stackPlot = F, heatMap = F, barPlot = F,  zPlot = F, 
		vioPlot = F, ylim, xlim, 
		returnMeta = T, returnPCA = F, rmZero = F,zJitter = F, main = NULL, sub = NULL, xlab = NULL, ylab = NULL,  ...){
	
	## check format of peaklist
	if( !is.data.frame(inDat) || nrow(inDat) < 1 || ncol(inDat) < 2 
			|| !any(names(inDat) == 'compound'))
		stop("A dataframe with a 'compound' and at least one data column is required")
	
	## Set Scale argument
	scale <- if (missing(scale)) 
				"none"
			else match.arg(scale)
	
	## Make all compound names unique
	com <- inDat$compound
	if(any(is.na(com)) &&  any(names(inDat) == 'medMz') 
			&& any(names(inDat) == 'medRt')){
		com[which(is.na(com))] <- paste(inDat$medMz[which(is.na(com))], 
				inDat$medRt[which(is.na(com))], sep = "@")
	}
	
	for( i in unique(com) ){
		tCpd <- com[ com == i ]
		if( length(tCpd) == 1 )
			next
		com[com == i] <- paste(tCpd, 1:length(tCpd), sep = '_')
	}
	inDat$compound <- com
	
	
	## Remove duplicate metabolite records and find max for each
	if(maxDup)
		inDat <- findMaxRecord( inDat, excludeName = excludeName )
	
	## Combine any metadata fields to be excluded
	tDat <- metaSep(inDat, excludeName = excludeName)
	metaDat <- tDat[[1]]
	inDat <- tDat[[2]]
	
	## Order samples by name
	if (nCst){
		inDat <- inDat[,order(names(inDat))]
	}
	
	## hold sample names
	sam <- names(inDat) ## read tab header to get sample names
	rownames(inDat) <- metaDat$compound
	
	## Zero any data below the data threshold
	if( !missing(dThresh) )
		inDat[ inDat < dThresh ] <- 0
	
	## Remove all zero data
	if( any(!apply(inDat, 1, function(x) all(x == 0))) && rmZero ){
		inDat <- inDat[- which(apply(inDat, 1, function(x) all(x == 0))),]		
	}
	
	## Condenses metabolite signals according to compound classes
	if (byClass ){
		
		if( any(is.na(metaDat$note)) )
			stop("The note column must contain compound classes")
		
		cpdClass <- metaDat$note			
		out <- list()
		for(i in 1:ncol(inDat))	
			out[[i]] <- aggregate(inDat[,i], by = list(cpdClass), "sum")
		inDat <- data.frame(sapply(out, function(x) x[,2]))
		metaDat <- data.frame( compound = out[[1]][,1], stringsAsFactors = F)
		names(inDat) <- sam
		rownames(inDat) <- metaDat$compound
	}
	
	## Get sample labels
	if( (barPlot || scale == 'zScore' || pCst || stackPlot || vioPlot ||
				linePlot || zPlot) && missing (sLab))
		sLab <- sapply( strsplit(sam, split = by), function (x) x[1])		
	
	
	## Average sample replicates
	if( avgRep || sdRep ){
		
		tLab <- list(factor(
						sapply(strsplit(names(inDat), split = by), function (x) x[1])))
		
		if( avgRep )
			inDat <- aggregate(t(inDat), FUN = function(x) mean(x, na.rm = T), by =tLab)
		else
			inDat <- aggregate(t(inDat), FUN = function(x) sd(x, na.rm = T), by =tLab)
		
		
		uName <- levels(inDat[,1]) 
		inDat <- data.frame(t(inDat[,-1]))
		names(inDat) <- uName
		
		## Reorder colums to match input list
		sam <- unique(sapply(strsplit(sam, split = by), function (x) x[1]))
		inDat <- inDat[,c(which(names(inDat) == 'compound'), 
						match(sam, names(inDat)))]
	}
	
	
	## Scale data
	if( scale == 'log'){
		inDat <- log( inDat )	
	}else if( scale == 'durbin' ){
		inDat <- durbin(inDat, c = c )	
	}else if (scale == "row") {
		inDat <- sweep(inDat, 1L, rowMeans(inDat, na.rm = TRUE), check.margin = FALSE)
		sx <- apply(inDat, 1L, sd, na.rm = TRUE)
		inDat <- sweep(inDat, 1L, sx, "/", check.margin = FALSE)
	}else if (scale == "column") {
		inDat <- sweep(inDat, 2L, colMeans(inDat, na.rm = TRUE), check.margin = FALSE)
		sx <- apply(inDat, 2L, sd, na.rm = TRUE)
		inDat <- sweep(inDat, 2L, sx, "/", check.margin = FALSE)
	}else if (scale == "fold" || scale == "ratio" || scale == 'zScore'){
		userSel <- select.list( names( inDat), multiple = T, 
				title = "Select reference:" )
		if( !length(userSel) || !nzchar(userSel) )
			stop('No file selected')
		normCol <- match(userSel, names(inDat))
		inDat <- foldChange( inTable = inDat, normVec = normCol, scale = scale)	
	}
	
	## Remove NA and Inf data
	if( na.rm  ){
		printList <- newList <- NULL
		
		if( any(is.na(inDat)) ){
			inDat[is.na(inDat)] <- 0
			printList <- c(printList, "NA")
			newList <- c(newList, 0)
		}
		
		if( any(inDat == Inf)){
			inDat[ inDat == Inf ] <- range(unlist(inDat), finite = T)[2]	
			newList <- c(printList, "Inf")
			newList <- c(newList, max(inDat, na.rm = T))
		}
		
		if( any(inDat == -Inf)){
			inDat[ inDat == -Inf ] <- range(unlist(inDat), finite = T)[1]	
			printList <- c(printList, "-Inf")
			newList <- c(newList, min(inDat, na.rm = T))
		}
		if(!is.null(printList))
			print( paste(c(paste(printList, collapse = ", "), 
									"data values were changed to:",	
									paste(newList, collapse = ", ")), collapse = " "))
	}	        
		
	## Allow users to select sub groups of data
	if( usrSelect ){
		userChoice <- mySelect(names(inDat))
		if( length(userChoice) > 0 ){
			inDat <- inDat[,userChoice]
			sLab <- sLab[userChoice]	
		}
	}
	## Generate PCA plots
	if( pca ){
		outPCA <- pcaPlot( dat = inDat, samples = sam, compounds = metaDat$compound, 
				transform = FALSE)		
		title(main = main, sub = sub, xlab = xlab, ylab = ylab)
	}

	## Get user selections
	if( barPlot || linePlot || zPlot || stackPlot || vioPlot){
		if( nrow(metaDat) ==  1){
			usrRow <- 1
		}else{
			userSel <- select.list( metaDat$compound, 
					multiple = T, title = "Select compounds:" )
			if( !length(userSel) || !nzchar(userSel) ){
				print('No compounds selected')
				break
			}
			usrRow <- match(userSel, metaDat$compound )	
		} 	
	}
	
	## Generate bar plot	
	if( barPlot ){
		userSel <- select.list( unique(sLab), multiple = T )
		if( !length(userSel) || !nzchar(userSel) ){
			print('No files selected')
			break
		}		
		usrCol <- myMatch(userSel, sLab )
		plotBar(sLab = sLab[usrCol], y = inDat[usrRow,usrCol], 
				cpdNames = row.names(inDat)[usrRow], ylim = ylim, ...)
		invisible( inDat[usrRow,usrCol] ) 	
	}
	
	## Generate XY plot
	if(linePlot){
		plotXY(sLab = sLab, y = inDat[usrRow,], cpdNames = row.names(inDat)[usrRow],
				stack = F, xlim = xlim, ylim = ylim, ...)		
		title(main = main, sub = sub, xlab = xlab, ylab = ylab)	
	} 
	
	## Generate stacked XY plot
	if(stackPlot){
		plotXY(sLab = sLab, y = inDat[usrRow,], cpdNames = row.names(inDat)[usrRow],
				stack = T, xlim = xlim, ylim = ylim, ...)
		title(main = main, sub = sub, xlab = xlab, ylab = ylab)
		
	} 
	
	## Generate violin plot
	if( vioPlot ){
		
		plotXY(sLab = sLab, y = inDat[usrRow,], cpdNames = row.names(inDat)[usrRow],
				stack = F, vioPlot = T, xlim = xlim, ylim = ylim, ...)
		title(main = main, sub = sub, xlab = xlab, ylab = ylab)
		
	}
	
	## Generate Z plot
	if(zPlot){
		plotXY(sLab = sLab, y = inDat[usrRow,], cpdNames = row.names(inDat)[usrRow],
				zPlot = T, xlim = xlim, ylim = ylim, zJitter = zJitter, ...)
		if( scale == "zScore" )
			abline( h=c(-1.644852, 1.644852), lty = 2)
		title(main = main, sub = sub, xlab = xlab, ylab = ylab)
		
	} 

	## Generate heatmap
	if( heatMap ){
		lDat <- levelFun(inDat, n = n, zlim = zlim, logScale = logScale )
		out <- metImage( as.matrix(inDat[nrow(inDat):1,]), col = lDat$col, breaks = lDat$breaks, 
				rCst = rCst, cCst = cCst, grid = grid )	
		title(main = main, sub = sub, xlab = xlab, ylab = ylab)
		
		rowInd <- rev(length(out$rowInd)-out$rowInd +1)
		colInd <- rev(length(out$colInd)-out$colInd +1)
		inDat <- inDat[rowInd, colInd]

		if(ncol(metaDat) ==  1){
			tmp <- data.frame(metaDat[rowInd,], stringsAsFactors = F)
			names(tmp) <- names(metaDat)
			metaDat <- tmp
			
		}else
			metaDat <- metaDat[rowInd,]
		
		
	}
	## Return PCA scores
	if(returnPCA)
		return(outPCA)
	## Add metadata
	if( returnMeta ){
		inDat <- cbind( metaDat, inDat)
	}else{
		inDat <- data.frame( compound = row.names(inDat), inDat, stringsAsFactors = F)
		names(inDat) <- c('compound', sam)
	}
	invisible(inDat)	
}

## Function for clustering untargeted mz data
## mzTable - dataframe - METen format dataframe from untargeted csv METen output
## mzThresh - numeric - Mass tolerance in ppm
## rtThresh - numeric - Renention time tolerance in minutes
## pW - numeric vector, 0 to 1 - weight given to mass, 
##                      retention time, and correlation
## pThresh - numeric, 0 to 1 - overall probablility threshold
## adList - numeric vector - list of known adducts
## collapse - logical - T collapses groups to highest intensity signal
## clusterOnly - logical - T returns collapsed groups with multiple signals 
METCluster <- function( mzTable, mzThresh = 10,rtThresh = .5,pW = c(1,1,1), pThresh = 0.75,
		adList = c(0, 1.003355, 1.007825, 22.9897692809, 18.010565, 43.99038, 34.968853),
        collapse = FALSE, clusterOnly = FALSE){
	
	mzList <- mzTable$medMz
	rtList <- mzTable$medRt
	pW <- pW/sum(pW) 

	## Remove all standard METen meta data
	METNam <- c('label', 'metaGroupId', 'groupId', 'goodPeakCount', 'medMz', 'medRt', 
			'maxQuality', 'note', 'compound', 'compoundId', 'expectedRtDiff', 
			'ppmDiff', 'parent', 'adductName', 'isotopeLabel', 'formula', 'pVal')	
	mzData <- mzTable[,-na.omit(match(METNam, names(mzTable)))] ## remove the non numeric data
		
	## Find mz diff matrix 
	mzMat <- matrix(NA, ncol = length(mzList), nrow = length(mzList) )
	for( i in 1:length(mzList) ){
		mzMat[i,] <- abs(mzList[i] - mzList)
	}

	## Calculate mass prob matrix
	myMin <- function( x ){ return( min(x, na.rm = TRUE))}
	mzProb <- matrix(NA, ncol = ncol(mzMat), nrow = nrow(mzMat) )	
	for( i in adList){
		mzTmp <- abs(mzMat - i) 
		mzTmp <- (mzTmp / mzList) * 10^6
		
		## Keep the best mz match
		mzIdx <- which(mzTmp < mzThresh) 
		tmpDat <- cbind(mzTmp[mzIdx], mzProb[mzIdx])
		mzProb[mzIdx] <- apply( tmpDat, 1, myMin )	
	}
	
	mzProb <- (abs(mzProb - mzThresh)/mzThresh)
	mzProb[which(is.na(mzProb))] <- 0

	## Find rt diff matix
	rtMat <- matrix(NA, ncol = length(rtList), nrow = length(rtList) )
	for( i in 1:length(rtList) ){
		rtMat[i,] <- abs(rtList[i] - rtList)
	}
	
	## Create linear probability matix with a threshold
	rtMat[rtMat > rtThresh ] <- NA
	rtMat <- abs((rtMat - rtThresh)/rtThresh)
	rtMat[which(is.na(rtMat))] <- 0
	
	## Find covarance
	corMat <- abs(cor( t(mzData) ))
	
	
	## Make weighted prob matrix
	probMat <- (mzProb * pW[1]) + (rtMat * pW[2]) + (corMat * pW[3])
	
	## Find compound grouping
	group <- apply(probMat, 1, function(x) which(x > pThresh))
	
	oLen <- sapply(group, length)
	nLen <- 1
	while( any(oLen != nLen) ){
		oLen <- sapply(group, length)
		for( i in 1:length(group) ){
			#i == 2201
			group[[i]] <- sort(unique(unlist(
									group[sapply( group, 
													function(x) any(match(x, group[[i]]), na.rm = T) )])))
		}
		nLen <- sapply(group, length)
	}
	
	group <- unique( group )
	group <- unique(group[order(sapply(group, length), decreasing = T)])
	grp <- nGrp <-  rep( NA, length(unlist(group)))

	for( i in 1:length(group) ){
		grp[ group[[i]] ] <- i		
		nGrp[ group[[i]] ]  <- length(group[[i]])
	}

	## Reorder list
	mzTable <- cbind( grp, nGrp, gPar = NA, mzTable)[order(grp),]
	grp <- grp[order(grp)]
	
	## Find parent isotope for each group
	uGrp <- unique(grp)
	datMax <- as.vector(apply(mzTable[,grep('_', names(mzTable))], 1, sum))

	for( i in 1:length(uGrp) ){
		idx <- which( grp == uGrp[i] )
		if( any(mzTable[idx,]$nGrp == 1)) 
			next
		idx1 <- idx[ which.max(datMax[idx]) ]
		idx2 <- idx[ which.min(mzTable$medMz[idx]) ]
		
		if( idx1 == idx2 )
			mzTable[idx1,]$gPar <- 'Parent'
		else{
			mzTable[idx1,]$gPar <- 'MaxI'
		}
	}

	mzTable <- unique(mzTable)
	
	## Collapse groups to parent masses
	if( collapse ){
		sTab <- mzTable[ mzTable$nGrp > 1, ]
		sTab <- sTab[ !is.na(sTab$gPar), ]
		if( clusterOnly )
			mzTable <- sTab
		else
			mzTable <- rbind( sTab, mzTable[ -(which(mzTable$nGrp > 1)), ])		
		mzTable$nGrp <- 1			
	}
	
	return(mzTable)
}
