# functions to calculate recombination and coupling statistics 

library(rhdf5)
#library(scales)
library(ape)
library(spatstat)


# Barton 1983's coupling coefficient (theta = s/r)
# Kruuk 1999's summed coupling coefficient (phi = (L-1)s/r) 
# Le (effective nLoci)  (Le = s*/s)




calcRecomb <- function(nLoci, map, nChrom= 4,...){
	# calculate recombination rate from numberLoci, map, and number of chromosomes
	if(nLoci <= nChrom){
		recomb <- 0.5 
		return(recomb)}
	else {
		nNeighbors <- nLoci - 1
		nNeighborsSame <- nLoci - nChrom
        nNeighborsDiff <- nNeighbors - nNeighborsSame
        sameChromAvgDist <- map/nLoci
        diffChromDist <- 0.5
        recomb <- ((nNeighborsSame * sameChromAvgDist) + (nNeighborsDiff * diffChromDist)) / nNeighbors
		recomb[which(recomb < 0 )] <- NaN
		recomb[which(recomb > 0.5)] <- 0.5
		return(recomb)
	}
}

singleLocusEq <- function(s, m, singleS= T){
	# function to calculate the equilibrium frequencies for (vectors of) single loci  (with singleS = T for vector of s values, and F when using one single s value (e.g. sMax))
	peq <- (m * (0.5 + 0.75* s) - 0.125* s - (0.125* (sqrt(s^2 - (4* m* s^2) + (4* m^2) * ((2 + s)^2))))) / (-1* (0.25 + m)*s )
	clineWidth <- (peq - (1 - peq))
	if (singleS == T) {
		return(unlist(list(peq, abs(clineWidth))))
	}
	else {
		out <- list()
		out[[1]] <- peq
		out[[2]] <- abs(clineWidth)
		return(out)
	}
}

calcLe <- function(m, pBar, s) {
	# calculate effective number of loci (L_e), with pBar = afDiffs
	sStar <- (m * (pBar - 0.5)) / ((0.25 - 0.25* pBar)* pBar + m * (0.5 + pBar * (pBar - 1.5)))
	Le <- sStar / s
	outp <- list(sStar, Le, s, pBar)
	names(outp) <- c('sStar', 'Le', 's', 'pBar')
	return(outp)
}

# same as above, but outputs unlist(x)
calcLeOld <- function(m, pBar, s) {
	sStar <- (m * (pBar - 0.5)) / ((0.25 - 0.25* pBar)* pBar + m * (0.5 + pBar * (pBar - 1.5)))
	Le <- sStar / s
	return(unlist(list(sStar, Le, s)))
}


getFSTs <- function(fstSubGen, i) {
	# calculate average Fst values per generation and locType (mean neutral, mean selected, mean total)
	meanSel <- mean(fstSubGen[[i]]$Fst[which(fstSubGen[[i]]$locType == 1)])
	meanNeut <- mean(fstSubGen[[i]]$Fst[which(fstSubGen[[i]]$locType == 0)])
	total <- mean(fstSubGen[[i]]$Fst)
	return(list(meanNeut, meanSel, total))
}


calcMaxEffMig <- function(s, m, L) {
	# calculate maximum effective migration for a given s vector, mutation rate and number of (selected) loci in one generation
	# NOTE:  only for single population
	
	# calculate equilibrium frequencies with given s   (here: sBar or meanS) 
	peq <- singleLocusEq(s, m, singleS= F)[[1]]
	q <- 1 - peq

	# calc fitnesses and maxEffMig
	randomResFit <- (1 + (peq * s))^L
	randomImmFit <- (1 + (q * s))^L
	patchAvgFit <- (randomImmFit * m) + ((1 - m) * randomResFit)
	#
	maxEffMig <- (randomImmFit * m) / patchAvgFit
	return(list(peq, maxEffMig))
}


calcGWCtime <- function(em, maxEffMig, endAllopatry){
# calculate the time of genome-wide congealingi, i.e. (last time step where effective migration was > maxEffMig ) + 1
	nGen <- unique(em[,1]) 
	emVec <- em[,2]
	t <- length(nGen)
	if (emVec[t] >= 0.5 * maxEffMig[t]) {
		# the effective migration rate can still be considered close to random, thus gwcTime not reached
		gwcTime <- NA
		gwcReached <- 0
	}
	else {
		gwcIndex <- tail(which(emVec >= maxEffMig), n= 1) + 1 
		if (length(gwcIndex) < 1) {
			gwcTime <- endAllopatry
		}
		else {
			gwcTime <- nGen[gwcIndex]
		}
		gwcReached <- 1
	}
	out <- list(gwcTime, gwcReached)
	names(out) <- c('gwcTime', 'gwcReached')
	return(out)
}


geomean = function(x, na.rm=TRUE){
	# function to calculate the geometric mean of vector x
	exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}



calcMeanS <- function(sSel, gen, mutDist, sConst) { 
	# calculate mean s for each generation in a given run
	if (mutDist == 3) {      # constant s
		meanS <- rep(sConst, gen)
		return(meanS)
	}
	else {
		meanS <- lapply(sSel, geomean)
		return(unlist(meanS))
	}
}





ccStats.1 <- function(df, fst, afts, LDsel, LDneut, effMig, run, maf= 25e-4, nChrom= 4) {
	effMig <- effMig[,1:2]
	# function to calculate coupling/congealing stats for a given run 
	df <- as.data.frame(df)
	params <- df[which(df$run == run),]
	m = params$sd_move
	fstSpl <- split(fst, fst$nGen)
	aftsSpl <- split(afts, afts$nGen)
	gen <- length(aftsSpl)
	sConst <- params$deme0_constant_s
	mutDist <- params$mutation_distribution
	#
	nLociS <- list()
	recom <- list()
	phiObs <- list()
	phiTilde <- list()
	kruukPhi_sMax <- list()
	si <- list()
	afDiffS <- list()
	afDiffN <- list()
	sBar <- list()
	sMax <- list()
	pHat <- list()
	clineWsMax <- list()
	pBarAllS <- list()
	clineWallS <- list()
	sStLeS <- list()
	FSTs <- list()
	LDsell <- list()
	LDneutl <- list()
	meanS <- list()
	sSelList <- list()
	avgAFdiff <- list()
	#
	for (i in 1:gen){
		mafThresh <- which(fstSpl[[i]]$allele_frequencies >= maf)
		#fstCurGen <- fstSpl
		#afCurGen <- aftsSpl
		fstCurGen <- fstSpl[[i]][mafThresh,]
		afCurGen <- aftsSpl[[i]][mafThresh,]
		#nGen <- fstGen$nGen
		locType <- which(fstCurGen$locType == 1)
		p <- fstCurGen$allele_frequencies
		sSel <- fstCurGen$S_MAX0[locType]
		sSelList[[i]] <- sSel
		#map <- fstCurGen$MAP
		#mapS <- map[locType]
	#
		afDiffS[[i]] <- afCurGen$AFdiff[locType]
		afDiffN[[i]] <- afCurGen$AFdiff[-locType]
		avgAFdiff[[i]] <- mean(abs(afCurGen$AFdiff))

		FSTs[[i]] <- getFSTs(fstSpl, i)	
	#
		nLoci <- dim(fstCurGen)[1]     
		nLociS[[i]] <- length(sSel)
		#nLociSmaf <- length(sSel[mafThresh])
		recom[[i]] <- calcRecomb(nLoci, 0.02*params$total_map_length, nChrom)  
	#
		sMax[[i]] <- prod(1 + sSel) - 1   
		kruukPhi_sMax[[i]] <- (nLociS[[i]] - 1) * (sMax[[i]]/recom[[i]])     
	#	
		phiObs[[i]] <- (sum(sSel) / recom[[i]])
		#phiTilde[[i]] <- phiObs[[i]] / nLociS[[i]]
		sBar[[i]] <- (phiObs[[i]] / nLociS[[i]]) * recom[[i]]
	#
			#	
		# single locus exp. under full coupling (sMax)
		clineWidthSmax <- singleLocusEq(sMax[[i]], m, singleS= T)
		pHat[[i]] <- clineWidthSmax[1]
		clineWsMax[[i]] <- clineWidthSmax[2]
	#
		# equilibrium freq & cline width for all sSel
		clineWidthAllS <- singleLocusEq(sSel, m, singleS= F)
		pBarAllS[[i]] <- clineWidthAllS[1]
		clineWallS[[i]] <- clineWidthAllS[2]
		#pBarAvg[[i]] <- #  "avgAFdiff in the L counted loci" (= obs. freq. of avg allele in the favored deme, pBar 
		#pBarAvg <- mean(unlist(pBarAllS[[i]]))
		Le <- calcLeOld(m, avgAFdiff[[i]], sBar[[i]])
		sStLeS[[i]] <- as.data.frame(matrix(Le, ncol= 4, dimnames= list(NULL, c('sStar', 'Le', 's', 'pBar'))))
		
		#bartonTheta[[i]] <- sum(si[[i]])/rVec[[i]]

	   	# if all loci have constant s (sMax == s); expected value (~) 
		#PhiObs <- sum(fstSub$S_MAX1[Sloci]) / rV 
	}
	meanS <- calcMeanS(sSelList, mutDist, gen, sConst)

	LDsell <- LDsel[, c(1,4)]
	LDneutl <- LDneut[, c(1,4)]
	FSTout <- as.data.frame(matrix(unlist(FSTs),ncol=3, byrow=TRUE, dimnames= list(NULL, c('FSTneut', 'FSTsel', 'FSTtot'))))
	sStLeS <- do.call('rbind', sStLeS)
#
	# maxEffMig with sBar (phi/ (L * r))
	maxEffMigSbar <- calcMaxEffMig(unlist(sBar), m, unlist(lapply(fstSpl, length)))[[2]]
	gwcTimeSbar <- calcGWCtime(effMig, maxEffMigSbar, params$end_period_allopatry)
	# maxEffMig with meanS (geomean of sMax[i])
	maxEffMigMeanS <- calcMaxEffMig(meanS, m, unlist(lapply(fstSpl, length)))[[2]]
	gwcTimeMeanS <- calcGWCtime(effMig, maxEffMigMeanS, params$end_period_allopatry)
	
	#
	runOut <- list(afDiffS, afDiffN, avgAFdiff, FSTout, unlist(recom), unlist(sMax), unlist(kruukPhi_sMax), unlist(phiObs), meanS, unlist(sBar), unlist(pHat), unlist(clineWsMax), pBarAllS, clineWallS, sStLeS, LDneutl, LDsell, effMig, unlist(maxEffMigSbar), gwcTimeSbar, m, unlist(maxEffMigMeanS), gwcTimeMeanS)
	names(runOut) <- c('afDiffS', 'afDiffN', 'avgAFdiff', 'FSTs', 'recomb', 'sMax', 'kruukPhi_sMax', 'phiObs', 'meanS', 'sBar', 'pHat', 'clineWidthSmax', 'pBarAllS', 'clineWallS','sStarLeS', 'LDneut', 'LDsel', 'effMig', 'maxEffMigSbar', 'gwcTimeSbar', 'sd_move', 'maxEffMigMeanS', 'gwcTimeMeanS')
		return(runOut)
}


### 


readCCobj <- function(run, smParam, path, ...) {
	# function to read a single data set (as input for ccStats)
	path5 <- paste('/runs/', run, sep= '')
	fstTmp <- as.data.frame(h5read(paste(path, smParam, '/Fst_', smParam, 'T.h5', sep= ''), name= path5)[[1]])
	H5close()
	aftsTmp <- as.data.frame(h5read(paste(path, smParam, '/afts_', smParam, 'T.h5', sep= ''), name= path5)[[1]])
	colnames(fstTmp) <- c("nGen", "locusID", "Fst", "allele_frequencies", "S_MAX1", "S_MAX0", "chromosomeMembership", "MAP", "locType")
	colnames(aftsTmp) <- c("nGen", "locusID", "AFpatch0", "AFpatch1", "is_reversed_locus", "locType", "AF", "AFdiff")
	H5close()
	#
	LDselTmp <- h5read(paste(path, smParam, '/LDselAvg_', smParam, 'T.h5', sep= ''), name= path5)[[1]]
	H5close()
	LDneuTmp <- h5read(paste(path, smParam, '/LDneutAvg_', smParam, 'T.h5', sep= ''), name= path5)[[1]]
	H5close()
	effMig <- h5read(paste(path, smParam, '/effMig_', smParam, 'T.h5', sep= ''), name= path5)[[1]]
	colnames(effMig) <- c("nGen", "eme0", "eme1", "nVariableLoci", "nRes", "nImm", paste('V', seq(7, 26,1)))
	H5close()
	dXY <- h5read(paste(path, smParam, '/dXY_', smParam, 'T.h5', sep= ''), name= path5)[[1]]
	colnames(dXY) <- c('nGen', 'dXY', 'deme0', 'deme1')
	H5close()
	#
	out <- list(fstTmp, aftsTmp, LDselTmp, LDneuTmp, effMig, dXY)
	names(out) <- c('fst', 'afts', 'LDsel', 'LDneut', 'effMig', 'dXY')
	return(out)
}

readCCobjRude <- function(run, smParam, folder, path, ...) {
	# function to read a single data set (as input for ccStats) (on ruderalis.colorado.edu)
	path5 <- paste('/runs/', run, sep= '')
	fstTmp <- as.data.frame(h5read(paste(path, folder, '/Fst_', smParam, 'T.h5', sep= ''), name= path5)[[1]])
	H5close()
	aftsTmp <- as.data.frame(h5read(paste(path, folder, '/afts_', smParam, 'T.h5', sep= ''), name= path5)[[1]])
	colnames(fstTmp) <- c("nGen", "locusID", "Fst", "allele_frequencies", "S_MAX1", "S_MAX0", "chromosomeMembership", "MAP", "locType")
	colnames(aftsTmp) <- c("nGen", "locusID", "AFpatch0", "AFpatch1", "is_reversed_locus", "locType", "AF", "AFdiff")
	H5close()
	#
	LDselTmp <- h5read(paste(path, folder, '/LDselAvg_', smParam, 'T.h5', sep= ''), name= path5)[[1]]
	H5close()
	LDneuTmp <- h5read(paste(path, folder, '/LDneutAvg_', smParam, 'T.h5', sep= ''), name= path5)[[1]]
	H5close()
	effMig <- h5read(paste(path, folder, '/effMig_', smParam, 'T.h5', sep= ''), name= path5)[[1]]
	colnames(effMig) <- c("nGen", "eme0", "eme1", "nVariableLoci", "nRes", "nImm", paste('V', seq(7, 26,1)))
	H5close()
	dXY <- h5read(paste(path, folder, '/dXY_', smParam, 'T.h5', sep= ''), name= path5)[[1]]
	colnames(dXY) <- c('nGen', 'dXY', 'deme0', 'deme1')
	H5close()
	#
	out <- list(fstTmp, aftsTmp, LDselTmp, LDneuTmp, effMig, dXY)
	names(out) <- c('fst', 'afts', 'LDsel', 'LDneut', 'effMig', 'dXY')
	return(out)
}


ccStats.2 <- function(run, df, ccObj, maf= 25e-4, nChrom= 4) {    #fst, afts, LDsel, LDneut, effMig, run, maf= 25e-4, nChrom= 4) {
	# function to calculate coupling/congealing stats for a given run 
	fst <- ccObj$fst
	afts <- ccObj$afts
	LDsel <- ccObj$LDsel #emRed <- em[em$nGen %in% unique(fst$nGen),]
	#LDsel <- as.data.frame(ccObj$LDsel)[ccObj$LDsel[,1] %in% names(ccObj$phiObs),]
	LDneut <- ccObj$LDneut
	effMig <- ccObj$effMig[,1:2]
	params <- df[which(df$run == run),]
	m = params$sd_move
	fstSpl <- split(fst, fst$nGen)
	aftsSpl <- split(afts, afts$nGen)
	nLoci <- unlist(lapply(lapply(fstSpl, dim), '[', 1))
	gen <- length(aftsSpl)
	fstSplS <- lapply(fstSpl, subset, locType == 1)
	aftsSplS <- lapply(aftsSpl, subset, locType == 1)
	nLociS <- as.numeric(unlist(lapply(fstSplS, nrow)))
	p <- lapply(fstSplS, '[[', 4)
	sSel <- lapply(fstSplS, '[[', 6)
	# 
	sConst <- params$deme0_constant_s
	mutDist <- params$mutation_distribution
	#
	# calc Kruuk's phi and phiObs 
	PHIs <- calcPHIs(aftsSpl, fstSpl, maf= maf, mapL= params$total_map_length, nChrom= params$nchromosomes)     # NOTE: we use a MAF threshold here! 
	#
	nLoci <- as.data.frame(cbind(nLoci, nLociS, PHIs$nLoci, PHIs$nLociS))
	colnames(nLoci) <- c('nLoci', 'nLociS', 'nLocimaf', 'nLociSmaf')
	recomb <- PHIs$recomb
	#  afDiffs
	afDiff_s <- lapply(aftsSplS, '[[', 8)
	afDiff_n <- lapply(lapply(aftsSpl, subset, locType == 0), '[[', 8)
	avgAFdiff <- unlist(lapply(lapply(lapply(aftsSpl, '[[', 8), abs), mean))    # mean(abs(allAFdiffsPerGen))
	#
	#and using actual selection coefficients (sSel == S_MAX0)
	# calc single locus exp. under full coupling (sMax) 
	sMaxlist <- lapply(PHIs$sMax, singleLocusEq, m= m, singleS= T)
	pHatsMax <- unlist(lapply(sMaxlist, '[[', 1))
	cWsMax <- unlist(lapply(sMaxlist, '[[', 2))
	# calc single locus exp. using actual selection coefficients (sSel == S_MAX0)
	sBarlist <- lapply(lapply(PHIs$sBar, singleLocusEq, m= m, singleS= F), unlist)
	pHatsBar <- unlist(lapply(sBarlist, '[[', 1))
	cWsBar <- unlist(lapply(sBarlist, '[[', 2))
	
	# equilibrium freq & cline width for all sSel
	clineWidthAllS <- lapply(sSel, singleLocusEq, m= m, singleS= F)
	pBarAllS <- lapply(clineWidthAllS, '[[', 1)
	clineWallS <- lapply(clineWidthAllS, '[[', 2)

	# get Fst values (s, n & total) per generation
	fstSel <- unlist(lapply(lapply(lapply(fstSpl, function(x)x[x$locType == 1,]), '[[', 3), mean))
	fstNeut <- unlist(lapply(lapply(lapply(fstSpl, function(x)x[x$locType == 0,]), '[[', 3), mean))
	fstAll <- unlist(lapply(lapply(fstSpl, '[[', 3), mean))
	FSTs <- as.data.frame(cbind(fstSel, fstNeut, fstAll))
	names(FSTs) <- c('FSTsel', 'FSTneut', 'FSTtot')
	# get effective s and Le   
	sSLe <- calcLe(m, avgAFdiff, PHIs$sBar)
	sStarLeS <- as.data.frame(do.call(cbind, sSLe))
	names(sStarLeS) <- c('sStar', 'Le', 's', 'pBar')
	meanS <- calcMeanS(sSel, mutDist, gen, sConst)
	# 
	LDsell <- LDsel[, c(1,4)]
	LDneutl <- LDneut[, c(1,4)]
	#FSTout <- as.data.frame(matrix(unlist(FSTs),ncol=3, byrow=TRUE, dimnames= list(NULL, c('FSTneut', 'FSTsel', 'FSTtot'))))
	######sStLeS <- do.call('rbind', sStLeS)
	# maxEffMig with sBar (phi/ (L * r))
	maxEffMigSbar <- calcMaxEffMig(PHIs$sBar, m, unlist(lapply(fstSpl, length)))[[2]]
	gwcTimeSbar <- calcGWCtime(effMig, maxEffMigSbar, params$end_period_allopatry)
	# maxEffMig with meanS (geomean of sMax[i])
	maxEffMigMeanS <- calcMaxEffMig(meanS, m, unlist(lapply(fstSpl, length)))[[2]]
	gwcTimeMeanS <- calcGWCtime(effMig, maxEffMigMeanS, params$end_period_allopatry)
##### output
	out <- list(FSTs, LDsell, LDneutl, afDiff_s, afDiff_n, avgAFdiff, PHIs$sMax, PHIs$kphiSmax, pHatsMax, cWsMax, PHIs$phiObs, PHIs$sBar, pHatsBar, cWsBar, meanS, sStarLeS, m, effMig, unlist(maxEffMigSbar), gwcTimeSbar, unlist(maxEffMigMeanS), gwcTimeMeanS, clineWallS, pBarAllS, nLoci, maf, recomb)
	names(out) <- c('FSTs', 'LDsel', 'LDneut', 'afDiffS', 'afDiffN', 'avgAFdiffs', 'sMax', 'kphisMax', 'pHatsMax', 'cWsMax', 'phiObs', 'sBar', 'pHatsBar', 'cWsBar', 'meanS', 'sStarLeS', 'sd_move', 'effMig', 'maxEffMigSbar', 'gwcTimeSbar', 'maxEffMigMeanS', 'gwcTimeMeanS', 'cWallS', 'pBarAllS', 'nLoci', 'maf', 'recomb')
	return(out)
}


ccStats.2slim <- function(run, df, ccObj, maf= 25e-4, nChrom= 4) {    #fst, afts, LDsel, LDneut, effMig, run, maf= 25e-4, nChrom= 4) {
	# function to calculate coupling/congealing stats for a given run 
	fst <- ccObj$fst
	afts <- ccObj$afts
	LDsel <- ccObj$LDsel #emRed <- em[em$nGen %in% unique(fst$nGen),]
	#LDsel <- as.data.frame(ccObj$LDsel)[ccObj$LDsel[,1] %in% names(ccObj$phiObs),]
	LDneut <- ccObj$LDneut
	effMig <- ccObj$effMig[,1:2]
	params <- df[which(df$run == run),]
	m = params$sd_move
	fstSpl <- split(fst, fst$nGen)
	aftsSpl <- split(afts, afts$nGen)
	nLoci <- unlist(lapply(lapply(fstSpl, dim), '[', 1))
	gen <- length(aftsSpl)
	fstSplS <- lapply(fstSpl, subset, locType == 1)
	aftsSplS <- lapply(aftsSpl, subset, locType == 1)
	nLociS <- as.numeric(unlist(lapply(fstSplS, nrow)))
	p <- lapply(fstSplS, '[[', 4)
	sSel <- lapply(fstSplS, '[[', 6)
	# 
	sConst <- params$deme0_constant_s
	mutDist <- params$mutation_distribution
	#
	# calc Kruuk's phi and phiObs 
	PHIs <- calcPHIs(aftsSpl, fstSpl, maf= maf, mapL= params$total_map_length, nChrom= params$nchromosomes)     # NOTE: we use a MAF threshold here! 
	#
	nLoci <- as.data.frame(cbind(nLoci, nLociS, PHIs$nLoci, PHIs$nLociS))
	colnames(nLoci) <- c('nLoci', 'nLociS', 'nLocimaf', 'nLociSmaf')
	#  afDiffs
	afDiff_s <- lapply(aftsSplS, '[[', 8)
	afDiff_n <- lapply(lapply(aftsSpl, subset, locType == 0), '[[', 8)
	avgAFdiff <- unlist(lapply(lapply(lapply(aftsSpl, '[[', 8), abs), mean))    # mean(abs(allAFdiffsPerGen))
	
	#and using actual selection coefficients (sSel == S_MAX0)
	# calc single locus exp. under full coupling (sMax) 
	sMaxlist <- lapply(PHIs$sMax, singleLocusEq, m= m, singleS= T)
	pHatsMax <- unlist(lapply(sMaxlist, '[[', 1))
	cWsMax <- unlist(lapply(sMaxlist, '[[', 2))
	# calc single locus exp. using actual selection coefficients (sSel == S_MAX0)
	sBarlist <- lapply(lapply(PHIs$sBar, singleLocusEq, m= m, singleS= F), unlist)
	pHatsBar <- unlist(lapply(sBarlist, '[[', 1))
	cWsBar <- unlist(lapply(sBarlist, '[[', 2))
	
	# equilibrium freq & cline width for all sSel
	clineWidthAllS <- lapply(sSel, singleLocusEq, m= m, singleS= F)
	pBarAllS <- lapply(clineWidthAllS, '[[', 1)
	clineWallS <- lapply(clineWidthAllS, '[[', 2)


	#
##### output
	out <- list(afDiff_s, afDiff_n, avgAFdiff, PHIs$kphiSmax, pHatsMax, PHIs$phiObs, clineWallS, pBarAllS, clineWidthAllS)
				
				#m, effMig, unlist(maxEffMigSbar), gwcTimeSbar, unlist(maxEffMigMeanS), gwcTimeMeanS, clineWallS, pBarAllS, clineWidthAllS, nLoci, maf)
	#names(out) <- c('afDiffS', 'afDiffN', 'avgAFdiffs', 'sMax', 'kphisMax', 'pHatsMax', 'cWsMax', 'phiObs', 'sBar', 'pHatsBar', 'cWsBar', 'meanS', 'sStarLeS', 'sd_move', 'effMig', 'maxEffMigSbar', 'gwcTimeSbar', 'maxEffMigMeanS', 'gwcTimeMeanS', 'pBarAllS', 'cWallS'   		, 'nLoci', 'maf')
	names(out) <- c('afDiffS', 'afDiffN', 'avgAFdiffs', 'kphiSmax', 'pHatsMax', 'phiObs', 'pBarAllS', 'clineWallS')
	return(out)
}





plotStatic.2 <- function(ccObj, run, data, i= '', ...) {
	# function to create the per-run plots with LD, Fst, PHIs, Le and me
	close.screen(all.screens= T)
	par(oma=c(4.5,4.5,4.5,4.5), mar=c(4.0,4.0,4.0,4.0) + 0.2, mgp= c(3,1,0))
	#par(bg = "white") # erase.screen() will appear not to work if the background color is transparent 
	split.screen(c(2,4))  
	#mtext(run, outer = TRUE )
	mainList <- c("neutral sites", "selected sites")
	LDlist <- list(ccObj$LDneut, ccObj$LDsel)
	param <- data[which(data$run == run),]
	m <- param$sd_move
	s <- param$mean_s
	mutdist <- param$mutation_distribution  #[which(data$run == run)]
	tsFreq <- param$ts_sampling_frequency
	xLab <- paste('sampling times (tsFreq =', tsFreq, ')')
	#
	# LD
	screen(1)
	plot(1:length(ccObj$phiObs), LDlist[[1]][,2], xlim= c(0, length(ccObj$phiObs)+2), ylim= c(0, 1), main= 'LD', ylab= 'Avg LD', xlab= xLab, type= 'l', col= 'blue')
	points(1:length(ccObj$phiObs), LDlist[[2]][,2], col= 'red', type= 'l')
	legend('bottomright', legend= c('neutral', 'selected'), fill= c('blue', 'red'))
	#
	# Fst 
	screen(2)
	cols <- c('blue', 'red', 'black') #t(rainbow(3))
	plot(ccObj$FSTs$FSTneut, type= 'l', col= cols[1], ylim= c(0,1), xlim= c(0, length(ccObj$phiObs)+2), ylab= expression('avg F'[st]), main= expression('avg F'[st]), xlab= xLab)
	points(ccObj$FSTs$FSTsel, col= cols[2], type= 'l')
	points(ccObj$FSTs$FSTtot, col= cols[3], type= 'l', lty= 1)
	legend('bottomright', legend= c('neutral', 'selected', 'total'), fill= cols, cex= 0.8)
	#
	# PHIs per generation
	screen(3) #, new= F)
	plot(log10(ccObj$phiObs), ylab= expression(paste('log'[10]~ phi)), xlab= xLab, type= 'l', ylim= c(-2, max(log10(ccObj$kphisMax))+ 1), col= 'grey70')  # xlim=  c(0, length(ccObj$phiObs)+100)
	points(log10(ccObj$kphisMax), type= 'l', col= 'black')
	legend('topleft', legend= c(expression(paste(phi[Kruuk])), expression(paste(phi, ' ', bar(s)))), fill= c('black', 'grey70'), cex= 0.8)
	text(x= 0.8*length(ccObj$meanS), y= 0.6*max(log10(ccObj$kphisMax)), paste('s = ', s))
	text(x= 0.8*length(ccObj$meanS), y= 0.5*max(log10(ccObj$kphisMax)), paste('m = ', m))

	#
	# cline widths and PHIs
	screen(4)
	cWallS <- lapply(ccObj$cWallS, unlist)
	phiOncW <- mapply(rep, ccObj$phiObs, times= unlist(lapply(cWallS, length)))
	#yLim <- max(unlist(ccObj$clineWidthSmax))
	plot(log10(ccObj$kphisMax), ccObj$pHatsMax, type= 'l', ylim= c(-0.45,  1), xlim= c(-2, max(log10(ccObj$kphisMax))), xlab= expression(paste('log'[10], ' ', phi)), ylab= expression(paste('p'[i0]~'- p'[i1])), col= 'black')
	#abline(h= 0, lty= 3)
	#}
	points(log10(unlist(phiOncW)), unlist(cWallS), pch= '.', col= 'grey70') #type= 'l')
	points(log10(unlist(ccObj$phiObs)), unlist(lapply(lapply(ccObj$afDiffS, abs), mean)), pch= '.', cex= 1.4, col= 'red')
	points(log10(unlist(ccObj$phiObs)), unlist(lapply(lapply(ccObj$afDiffN, abs), mean)), pch= '.', cex= 1.4, col= 'blue')
	#points(log10(unlist(ccObj$phiObs)), unlist(ccObj$avgAFDiff), type= 'l', col= 'black')

	legend('bottomright', legend= c(expression(paste(phi[Kruuk], ' ~ ', 'peq '[sMax])), expression(paste(phi, ' ~ ', bar(p), ' all s')), expression(paste(phi, ' ~ avg p S')), expression(paste(phi, ' ~ avg p N'))), fill= c('black', 'grey70', 'red', 'blue'), cex= 0.8)
	#
	# plot Le and me
	screen(5)
	#plot(1, type="n", xlab= 'gen / 1e3', ylab="", ylim=c(-0.1, 0.1), xlim=c(0, (ccObj$sStarLeS)))
	emig <- ccObj$effMig[,2]
	emig[which(emig == 0)] <- 1e-10
	plot(log10(abs(ccObj$sStarLeS$Le)), type= 'l', xlab= xLab, ylab= expression(paste('log'[10], '.')) , ylim= c(-10, 7))#, ylim= c(min(log10(emig)), max(log10(abs(ccObj$sStarLeS$Le)))))      # ylab was expression(paste('log'[10]~'L'[e]))
	abline(v= which.min(log10(abs(ccObj$sStarLeS$Le))), lty= 2, col= 'black')
	points(log10(emig), col= 'blue', type= 'l')
	points(log10(ccObj$sBar), type= 'l', col= 'orange')
	points(log10(ccObj$sMax), type= 'l', col= 'red')
	abline(v= ccObj$gwcTimeMeanS$gwcTime/tsFreq, lty= 2, col= 'green')
	legend('bottomright', legend= c(expression('L'[e]), expression('m'[e]), 'gwcTime meanS', expression(bar(s)), 'sMax'), lty= c(1, 1, 2), col= c('black', 'blue', 'green', 'orange', 'red'), cex= 0.8)
	text(x= 0.6*length(ccObj$meanS), y= 0.6*max(log10(emig)), paste("gwc time = ", ccObj$gwcTimeMeanS$gwcTime))
	#
	# plot effMig  
	screen(6)
	plot(ccObj$maxEffMigSbar, col= 'orange', type= 'l', ylim= c(0, 0.25), ylab= 'migration rate', xlab= xLab)
	points(ccObj$effMig[,2], col= 'black', type= 'l')
	points(ccObj$maxEffMigMeanS, col= 'green', type= 'l')
	points(ccObj$recomb, type= 'l', col= 'cyan')

	abline(h= m, col= 'red', lty= 3)
	abline(h= s, col= 'blue', lty= 3)
	abline(v= ccObj$gwcTimeSbar$gwcTime/tsFreq, col= 'orange', lty= 2)
	abline(v= ccObj$gwcTimeMeanS$gwcTime/tsFreq, col= 'green', lty= 2)
	endAllo <- data$end_period_allopatry[which(data$run == run)]
	if (ccObj$gwcTimeMeanS$gwcTime == endAllo) {
		text(x= 0.3* length(ccObj$meanS), y= 1.1*max(ccObj$effMig[,2]), expression(paste('gwcTime == endPallo')))
	}
	text(x= 0.8* length(ccObj$meanS), y= 0.4*max(ccObj$effMig[,2]), paste('mutdist = ', mutdist)) 
	text(x= 0.3* length(ccObj$meanS), y= 1.3*max(ccObj$effMig[,2]), paste('endPallo = ', param$end_period_allopatry))
	cols <- c('orange', 'green', 'black', 'cyan', 'red', 'blue', 'orange', 'green')
	legend('topright', legend= c(expression(paste('maxEffMig ', bar(s))), 'maxEffMig meanS', expression(paste('m'[e0])), 'recomb.', 'm', 's', expression(paste('gwcTime ', bar(s))), 'gwcTime meanS'), col= cols, lty= c(1,1,1,3,3,2,2), pch= c(NA,NA,NA,NA,NA,NA,NA), pt.cex= 2, seg.len= 0.9, cex= 0.7)
	#
	screen(7)
	maflgnd <- paste('MAF (', ccObj$maf, ')', sep= '')
	plot(ccObj$nLoci$nLoci, type= 'l', ylab= 'numbers of loci', ylim= c(0, 1.1*max(ccObj$nLoci$nLoci)))
	points(ccObj$nLoci$nLociS, type= 'l', col= 'red')
	points(ccObj$nLoci$nLocimaf, type= 'l', col= 'blue')
	points(ccObj$nLoci$nLociSmaf, type= 'l', col= 'yellow')
	legend('topleft', legend= c('total', paste('total', maflgnd, sep= ' '), 'selected', paste('selected ',maflgnd, sep= ' ')), fill= c('black', 'blue', 'red', 'yellow'), cex= 0.8)

	#
	screen(8)
	plot(unlist(lapply(lapply(ccObj$afDiffS, abs), mean)), pch= '.', col= 'red', ylab= 'allele freq diffs', cex= 1.4)
	points(unlist(lapply(lapply(ccObj$afDiffN, abs), mean)), pch= '.', col= 'blue')
	legend('bottomright', legend= c('neutral', 'selected'), fill= c('blue', 'red'))
	#
	close.screen(all= T)
	title(paste(run, ' (',i, ')', sep= ''), outer= T)
}


wrapH5static <- function(data, setname, path= '/media/schimar/dapperdata/bu2s/h5/', maf= 25e-4, sleep= 0,...) {
	# function to read individual runs (from vector of runs), calculate CC and plotStatic
	#
	for (i in 1:dim(data)[1]){
		run <- data$run[i]
		path5 <- paste('/runs/', run, sep= '')
		#
		ccObjTmp <- readCCobj(run, setname, path)
		#ccTmp <- ccStats.2(data, ccObjTmp$fst, ccObjTmp$afts, ccObjTmp$LDsel, ccObjTmp$LDneut, ccObjTmp$effMig, run, maf= maf)
		ccTmp <- ccStats.2(run= run, df= df, ccObj= ccObjTmp, maf= maf)
		#
		#plotStatic.1(ccTmp, run, data)
		plotStatic.2(ccTmp, run, data, i= i)
		#
		Sys.sleep(sleep)
		H5close()
	}
}

wrapH5dynaGen <- function(data, setname, path= '/media/schimar/schimar2/bu2s/h5/', maf= 25e-4, sleep= 0, time= 1, static= F,...) {
	# function to read individual runs (from vector of runs), calculate CC and plotStatic
	#
	for (i in 1:dim(data)[1]){
		run <- data$run[i]
		path5 <- paste('/runs/', run, sep= '')
		#
		ccObjTmp <- readCCobj(run, setname, path)
		#ccTmp <- ccStats.2(data, ccObjTmp$fst, ccObjTmp$afts, ccObjTmp$LDsel, ccObjTmp$LDneut, ccObjTmp$effMig, run, maf= maf)
		ccTmp <- ccStats.2(run= run, df= df, ccObj= ccObjTmp, maf= maf)
		#
		fstSpl <- split(ccObjTmp$fst, ccObjTmp$fst$nGen)
		mIKgen <- calcMorIripK(fstSpl) 
		#plotStatic.1(ccTmp, run, data)
		plotDynaMorIgen(ccTmp, mIKgen[[1]], run, time= time, static= static, wait= 0.3)
		#
		Sys.sleep(sleep)
		H5close()
	}
}
# ccObj, morIdata, time= 1, static= F, wait= 0, ...)

wrapH5dynaBin <- function(data, setname, path= '/media/schimar/schimar2/bu2s/h5/', maf= 25e-4, sleep= 0, time= 1, static= F,...) {
	# function to read individual runs (from vector of runs), calculate CC and MorIbin, and then plot Fst vs MAP and MorIbin 
	#
	for (i in 1:dim(data)[1]){
		run <- data$run[i]
		path5 <- paste('/runs/', run, sep= '')
		#
		ccObjTmp <- readCCobj(run, setname, path)
		#ccTmp <- ccStats.2(data, ccObjTmp$fst, ccObjTmp$afts, ccObjTmp$LDsel, ccObjTmp$LDneut, ccObjTmp$effMig, run, maf= maf)
		ccTmp <- ccStats.2(run= run, df= df, ccObj= ccObjTmp, maf= maf)
		#
		fstSpl <- split(ccObjTmp$fst, ccObjTmp$fst$nGen)
		mIbin <- calcMorIbin(fstSpl) 
		plotFstMAPmorIbin(fstSpl, mIbin, time= 1, static= F, wait= 0.1, ...)
		
		
		#plotStatic.1(ccTmp, run, data)
		#plotDyna(ccTmp, mIKgen[[1]], run, time= time, static= static, wait= 0.3)
	
		#
		Sys.sleep(sleep)
		H5close()
	}
}

wrapH5FstMAPmIgen <- function(data, setname, path= '/media/schimar/schimar2/bu2s/h5/', maf= 25e-4, sleep= 0, time= 1, static= F,...) {
	# function to read individual runs (from vector of runs), calculate CC and MorIbin, and then plot Fst vs MAP and MorIbin 
	#
	for (i in 1:dim(data)[1]){
		run <- data$run[i]
		path5 <- paste('/runs/', run, sep= '')
		#
		ccObjTmp <- readCCobj(run, setname, path)
		#ccTmp <- ccStats.2(data, ccObjTmp$fst, ccObjTmp$afts, ccObjTmp$LDsel, ccObjTmp$LDneut, ccObjTmp$effMig, run, maf= maf)
		ccTmp <- ccStats.2(run= run, df= df, ccObj= ccObjTmp, maf= maf)
		#
		fstSpl <- split(ccObjTmp$fst, ccObjTmp$fst$nGen)
		mIKgen <- calcMorIripK(fstSpl) 
		plotFstMAPmorIgen(fstSpl, mIKgen, time= time, static= static, wait= 0.1, ...)
		
		
		#plotStatic.1(ccTmp, run, data)
		#plotDyna(ccTmp, mIKgen[[1]], run, time= time, static= static, wait= 0.3)
	
		#
		Sys.sleep(sleep)
		H5close()
	}
}



readIn <- function(run, data, setname, path= '/media/schimar/dapperdata/bu2s/h5/', ...) {
	# function to read single run and return ccStats output (so far, not really needed. refer to the two called functions instead(readCCobj & ccStats.2)
	#
	
	path5 <- paste('/runs/', run, sep= '')
	#
	ccObjTmp <- readCCobj(run, setname, path)
	ccTmp <- ccStats.1(data, ccObjTmp$fst, ccObjTmp$afts, ccObjTmp$LDsel, ccObjTmp$LDneut, ccObjTmp$effMig, run)
	#ccTmp <- ccStats.2(run, data, ccObjTmp)
	#
	return(ccTmp)
	H5close()
}

#r202971 <- readIn('Run202971', df, 'Sm')





### 
calcPHIs <- function(aftsSpl, fstSpl, maf= 25e-4, mapL= 100, nChrom= 4) {
	#afts <- subset(afts, fst$allele_frequencies >= maf)
	#fst <- subset(fst, fst$allele_frequencies >= maf)
	aftsSplmaf <- lapply(aftsSpl, function(x)x[x$AF >= 0.1,])
	fstSplmaf <- lapply(fstSpl, function(x)x[x$allele_frequencies >= 0.1,])
	nLoci <- unlist(lapply(lapply(fstSplmaf, dim), '[', 1))
	# 
	fstSplS <- lapply(fstSplmaf, subset, locType == 1)
	aftsSplS <- lapply(aftsSplmaf, subset, locType == 1)
	nLociS <- as.numeric(unlist(lapply(fstSplS, nrow)))
	p <- lapply(fstSplS, '[[', 4)
	sSel <- lapply(fstSplS, '[[', 6)
	# calc recombination 
	recomb <- unlist(lapply(nLociS, calcRecomb, map= 0.02*mapL, nChrom= nChrom))
	# calc sMax
	sMax <- unlist(lapply(sSel, function(x) prod(1 + x)-1))
	# calc Kruuk's phi 
	kphiSmax <- ((unlist(nLociS) - 1) * (sMax / recomb))
	#
	# phiObs & sBar
	phiObs <- unlist(lapply(sSel, sum)) / recomb
	sBar <- (phiObs / unlist(nLociS)) * recomb
	#
	out <- list(sMax, kphiSmax, phiObs, sBar, nLociS, nLoci, recomb)
	names(out) <- c('sMax', 'kphiSmax', 'phiObs', 'sBar', 'nLociS', 'nLoci', 'recomb')
	return(out)
}



xtractPhis <- function(data, setname, folder, path= '/media/schimar/FLAXMAN/h5/', maf= 25e-4, ...) {
	# function to read individual runs (from vector of runs), calculate CC and create new list (of length(data)) that contains phiObs and kphismax
	#
	runs <- list()
	#kphisMax <- list()
	for (i in 1:dim(data)[1]){
		run <- data$run[i]
		path5 <- paste('/runs/', run, sep= '')
		#
		ccObjTmp <- readCCobjRude(run, setname, folder, path)
		#ccTmp <- ccStats.2(data, ccObjTmp$fst, ccObjTmp$afts, ccObjTmp$LDsel, ccObjTmp$LDneut, ccObjTmp$effMig, run, maf= maf)
		ccTmp <- ccStats.2(run= run, df= df, ccObj= ccObjTmp, maf= maf)
		#ccTmp <- ccStats.2slim(run= run, df= df, ccObj= ccObjTmp, maf= maf)
		#
		avgAFdiffS <- unlist(lapply(lapply(ccTmp$afDiffS, abs), mean))
		avgAFdiffN <- unlist(lapply(lapply(ccTmp$afDiffN, abs), mean))
		cWallS <- lapply(ccTmp$cWallS, unlist)
		runs[[i]] <- list(ccTmp$phiObs, ccTmp$kphisMax, ccTmp$pHatsMax, avgAFdiffS, avgAFdiffN, cWallS, ccTmp$pBarAllS)   # 
		names(runs)[i] <- run
		names(runs[[i]]) <- c('phiObs', 'kphisMax', 'pHatsMax', 'afDiffS', 'afDiffN', 'cWallS', 'pBarAllS')
		#phiObs[[i]] <- ccTmp$phiObs
		#names(phiObs)[i] <- run
		#kphisMax[[i]] <- ccTmp$kphisMax
		#names(kphisMax)[i] <- run
	}
	#out <- list(phiObs, kphisMax)
	#names(out) <- c('phiObs', 'kphisMax')
	return(runs)
}		
# maybe write another function ('ccStats.2 abgespeckt') to calcPHIs and get afDiffs  (so it doesn't take as outlandischly long to get this...) 

phi2plot <- function(phis, ...) {
	# function to get the output for ggplot2 plots from the output of xtractPhis
	kps <- unlist(lapply(phis, '[[', 2))
	pHs <- unlist(lapply(phis, '[[', 3))
	
	phO <- unlist(lapply(phis, '[[', 1))
	pN <- unlist(lapply(phis, '[[', 5))
	pS <- unlist(lapply(phis, '[[', 4))
	
	cWs <- unlist(lapply(phis, '[[', 6))
	
	
	# get the 'phiOncW' for selected and neutral afDiffs
	phiOncW <- list()
	for (i in 1:length(phis)) {
		#k <- i + 800
		phiOncW[[i]] <- mapply(rep, phis[[1]]$phiObs, times= unlist(lapply(phis[[i]]$cWallS, length)))
	}
	# flatten 
	pcWallS <- unlist(lapply(phiOncW, unlist))

	### cbind phi and AFs/cline widths 
	# kruuk: 
	# 		(kps_pHs)
	kps_pHs <- as.data.frame(cbind(log10(kps), pHs, rep('C', length(kps))))
	names(kps_pHs) <- c('phi', 'p', 'group')
		
	# observed phis ~ avgAFdiff(S|N) 
	phiOAFd <- as.data.frame(cbind(phO, pN, pS))
	names(phiOAFd) <- c('phiObs', 'pN', 'pS')
	
	Sphi <- as.data.frame(cbind(as.numeric(log10(phiOAFd$phiObs)), as.numeric(phiOAFd$pS), rep('A', length(phiOAFd$pS))))	
	names(Sphi) <- c('phi', 'p', 'group')
	Nphi <- as.data.frame(cbind(as.numeric(log10(phiOAFd$phiObs)), as.numeric(phiOAFd$pN), rep('B', length(phiOAFd$pN))))	
	names(Nphi) <- c('phi', 'p', 'group')
	# obsPHI <- rbind(Sphi, Nphi)
	
	# expectations based on per-locus s:
	allS <- as.data.frame(cbind(log10(pcWallS), cWs, rep('D', length(cWs))))
	colnames(allS) <- c('phi', 'p', 'group')
	
	# all together (used for ggplot2)
	# allPhi <- as.data.frame(rbind(kps_pHs, Sphi, Nphi, allS))
	# colnames(allPhi) <- c('phi', 'p', 'group')
	allPhi <- list(kps_pHs, Sphi, Nphi, allS)
	names(allPhi) <- c('kphi', 'pS', 'pN', 'allS')
	return(allPhi)
}


xtractLe <- function(data, setname, folder, path= '/media/schimar/FLAXMAN/h5/', maf= 25e-4, ...) {
	# function to read individual runs (from vector of runs), calculate CC and create new list (of length(data)) that contains effMig, Le and gwcTime  
	#
	runs <- list()
	#kphisMax <- list()
	for (i in 1:dim(data)[1]){
		run <- data$run[i]
		path5 <- paste('/runs/', run, sep= '')
		#
		ccObjTmp <- readCCobjRude(run, setname, folder, path)
		#ccTmp <- ccStats.2(data, ccObjTmp$fst, ccObjTmp$afts, ccObjTmp$LDsel, ccObjTmp$LDneut, ccObjTmp$effMig, run, maf= maf)
		ccTmp <- ccStats.2(run= run, df= df, ccObj= ccObjTmp, maf= maf)
		runs[[i]] <- list(ccTmp$sStarLeS, ccTmp$effMig, ccTmp$maxEffMigMeanS, ccTmp$gwcTimeMeanS, data$ts_sampling_frequency[i])   # 
		names(runs)[i] <- run
		names(runs[[i]]) <- c('sStarLeS', 'effMig', 'maxEffMigMeanS', 'tsFreq')
		#phiObs[[i]] <- ccTmp$phiObs
		#names(phiObs)[i] <- run
		#kphisMax[[i]] <- ccTmp$kphisMax
		#names(kphisMax)[i] <- run
	}
	#out <- list(phiObs, kphisMax)
	#names(out) <- c('phiObs', 'kphisMax')
	return(runs)
}		


xtractLD <- function(data, setname, folder, path= '/media/schimar/FLAXMAN/h5/', maf= 25e-4, ...) {
	# function to read individual runs (from vector of runs), calculate CC and create new list (of length(data)) that contains effMig, Le and gwcTime  
	#
	runs <- list()
	#kphisMax <- list()
	for (i in 1:dim(data)[1]){
		run <- data$run[i]
		path5 <- paste('/runs/', run, sep= '')
		#
		ccObjTmp <- readCCobjRude(run, setname, folder, path)
		#ccTmp <- ccStats.2(data, ccObjTmp$fst, ccObjTmp$afts, ccObjTmp$LDsel, ccObjTmp$LDneut, ccObjTmp$effMig, run, maf= maf)
		#ccTmp <- ccStats.2(run= run, df= df, ccObj= ccObjTmp, maf= maf)
		#
		runs[[i]] <- list(ccObjTmp$LDneut, ccObjTmp$LDsel, ccObjTmp$dXY)   # 

		names(runs)[i] <- run
		names(runs[[i]]) <- c('LDneut', 'LDsel', 'dXY')
	}
	#out <- list(phiObs, kphisMax)
	#names(out) <- c('phiObs', 'kphisMax')
	return(runs)
}		


xtractIK <- function(data, setname, folder, path= '/home/schimar/FLAXMAN/h5/', maf= 25e-4, ...) {
	# function to read individual runs (from vector of runs), calculate CC and create new list (of length(data)) that contains effMig, Le and gwcTime  
	#
	runs <- list()
	#kphisMax <- list()
	for (i in 1:dim(data)[1]){
		run <- data$run[i]
		path5 <- paste('/runs/', run, sep= '')
		#
		ccObjTmp <- readCCobjRude(run, setname, folder, path)
		#ccTmp <- ccStats.2(data, ccObjTmp$fst, ccObjTmp$afts, ccObjTmp$LDsel, ccObjTmp$LDneut, ccObjTmp$effMig, run, maf= maf)
		#ccTmp <- ccStats.2(run= run, df= df, ccObj= ccObjTmp, maf= maf)
		#
		fstSpl <- split(ccObjTmp$fst, ccObjTmp$fst$nGen)
		#runs[[i]] <- list(ccObjTmp$LDneut, ccObjTmp$LDsel, ccObjTmp$dXY)   # 
		runs[[i]] <- calcMorIripK(fstSpl)

		names(runs)[i] <- run
		#names(runs[[i]]) <- c('LDneut', 'LDsel', 'dXY')
	}
	#out <- list(phiObs, kphisMax)
	#names(out) <- c('phiObs', 'kphisMax')
	return(runs)
}		


#############################################

		######## Moran's I ########

calcMorIripK <- function(fstspl, windowLen= 25, ...) {
# function to calculate Moran's I per generation and chromosome, with INPUT: fst data split by generation. 	
	mIgenAll <- list()
	KgenAll <- list()
	pointPatt <- list()
	for (i in 1:length(fstspl)) {
		cFst <- fstspl[[i]]
		cChromSpl <- split(cFst, cFst$chromosomeMembership)
		mIchrom <- as.data.frame(matrix(nrow= 4, ncol= 4, dimnames= list(seq(1, 4, 1), c('obs', 'exp', 'sd', 'pval'))))
		k4chrom <- list()
		ppchrom <- list()
		for (j in 1:length(cChromSpl)) {
			cChrom <- cChromSpl[[j]]
			ctFst <- cChrom$Fst
			ctMAP <- cChrom$MAP
			# Ripley's K-function
			pp <- ppp(ctMAP, rep(1, length(ctMAP)), window= owin(c(0, windowLen), c(0,1)))
			ppchrom[[j]] <- pp
			Kpp <- Kest(pp, correction= 'isotropic')
			k4chrom[[j]] <- Kpp
			# Moran's I
			ctDists.inv <- 1/as.matrix(dist(ctMAP))
			diag(ctDists.inv) <- 0
			ctDists.inv[is.infinite(ctDists.inv)] <- 0
			morI <- unlist(Moran.I(ctFst, ctDists.inv, scaled= T, alt= 'greater', na.rm= T))
			mIchrom[j,] <- morI
		}
		mIgenAll[[i]] <- mIchrom
		KgenAll[[i]] <- k4chrom
		pointPatt[[i]] <- ppchrom
	}
	return(list(mIgenAll, KgenAll, pointPatt))
}



calcMorIbin <- function(fstspl, nbin= 5, ...) {
	# function to calculate Moran's I per generation and chromosome and distance bin
	# INPUT: fst data split by generation and the number of distance bins (default 5). The functiuon outputs a list containing the Moran's I nested list, and the bin affiliation for each generation/chromosome/locus. 
	mIgen <- list()
	grouplist <- list()
	for (i in 1:length(fstspl)) {
		cFst <- fstspl[[i]]
		cChromSpl <- split(cFst, cFst$chromosomeMembership)
		mIchrom <- list()
		ctlist <- list()
		for (j in 1:length(cChromSpl)) {
			cChrom <- cChromSpl[[j]]
			dists <- dist(cChrom$MAP)
			hc <- hclust(dists, method= 'ward.D')
			ct <- cutree(hc, k= nbin)
			ctlist[[j]] <- ct
			mIbin <- as.data.frame(matrix(nrow= length(unique(ct)), ncol= 4, dimnames= list(seq(1, length(unique(ct)), 1), c('obs', 'exp', 'sd', 'pval'))))
			for (k in 1:length(unique(ct))) {
				ctFst <- cChrom$Fst[which(ct == k)]
				ctMAP <- cChrom$MAP[which(ct == k)]
				ctDists.inv <- 1/as.matrix(dist(ctMAP))
				diag(ctDists.inv) <- 0
				ctDists.inv[is.infinite(ctDists.inv)] <- 0
				morI <- Moran.I(ctFst, ctDists.inv, scaled= T, alt= 'greater', na.rm= T)
				#ctDists <- mat2listw(as.matrix(dist(ctMAP)))
				#morI <- moran.test(ctFst, ctDists, na.action= na.omit, zero.policy= T, rank= T, adjust.n= T)
				#morI <- moran(ctFst, ctDists, n= 1, S0= Szero(ctDists), NAOK= T)
				mIbin[k, ] <- unlist(morI)
	
			}
			
			mIchrom[[j]] <- mIbin
		}
		grouplist[[i]] <- ctlist
		mIgen[[i]] <- mIchrom
	}
	return(list(mIgen, grouplist))
}




########################
# plots for MorI

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
	if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
		stop("vectors must be same length")
	arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}



plotMorIbin <- function(morIgen, time= 1) {
	# function to plot Moran's I values per generation, chromosome and distance bin, with a list of length = nGen containing the Moran's I values (nested lists) and the starting time for the loop as input (nGen/tsFreq). 
	par(mfrow= c(1,4))
	
	for (i in time:length(morIgen)) {
		for (j in 1:length(morIgen[[i]])) {
			cMorI <- morIgen[[i]][[j]]
			barx <- barplot(cMorI$obs, ylim= c(-1, 1), pch= 20, main= paste("Moran's I, gen = ", i), names.arg= c(1,2,3,4,5), ylab= "Moran's I", xlab= 'distance bins')  
			sdI <- cMorI$sd
			sdI[which(sdI == Inf)] <- NA
			sdI[which(sdI == NaN)] <- NA
			error.bar(barx, cMorI$obs, sdI)
		}
		Sys.sleep(0.5)
	}
}



plotFstMAPmorIbin <- function(fstspl, morIbin, time= 1, static= F, wait= 0.1, ...) {
	# function to plot both the MAP vs Fst and Moran's I and Moran's I per generation, chromosome and distance bin
	# INPUT: fst data split by nGen, nested list of Moran's I values (per gen, chrom and bin) and the generation time to start the loop (with nGen/tsFreq).
	close.screen(all.screens= T)

	par(bg = "white") # erase.screen() will appear not to work if the background color is transparent 
	split.screen(c(2,4))  


	#par(mfrow= c(2,4))
	for(i in time:length(fstspl)) {
		chrom <- fstspl[[i]]$chromosomeMembership
		for (j in 1:length(unique(chrom))) {
			screen(j)
			cChrom <- fstspl[[i]][which(chrom == unique(chrom)[j]),]
			plot(cChrom$MAP, cChrom$Fst, col= as.factor(cChrom$locType), xlim= c(0, 25), ylim= c(0,1), main= paste('chrom ', j-1, ' gen = ', i), ylab= 'Fst', xlab= 'map', pch= 20)
			lines(cChrom$MAP, cChrom$Fst, col= 'grey70')
		
		#for (j in 1:length(morIbin[[i]])) {
			screen(j+4)
			cMorI <- morIbin[[1]][[i]][[j]]
			barx <- barplot(cMorI$obs, ylim= c(-1, 1), pch= 20, main= paste("Moran's I, gen = ", i), names.arg= c(1,2,3,4,5), ylab= "Moran's I", xlab= 'distance bins')      # col= grouplist[[i]][[j]]
			sdI <- cMorI$sd
			sdI[which(sdI == Inf)] <- NA
			sdI[which(sdI == NaN)] <- NA
			error.bar(barx, cMorI$obs, sdI)
		}
		if (static == T) {
				break
			}
		Sys.sleep(wait)

	}
}


plotFstMAPripKmorIbin <- function(fstspl, morIbin, mIKgen, time= 1, static= F, wait= 0.1, ...) {
	# function to plot both the MAP vs Fst, Ripley's K-function (with envelope)  and Moran's I per generation, chromosome and distance bin
	# INPUT: fst data split by nGen, nested list of Moran's I values (per gen, chrom and bin) and the generation time to start the loop (with nGen/tsFreq).

	close.screen(all.screens= T)
	par(bg = "white") # erase.screen() will appear not to work if the background color is transparent 
	split.screen(c(3,4))  

	for (i in time:length(fstspl)) {
		chrom <- fstspl[[i]]$chromosomeMembership
		cGmI <- mIKgen[[1]][[i]]
		cGpp <- mIKgen[[3]][[i]]
		cGmIbin <- mIbin[[1]][[i]]
		for (j in 1:length(cGpp)) {
			cCpp <- cGpp[[j]]
			cMorI <- cGmIbin[[j]]
			screen(j)
			cChrom <- fstspl[[i]][which(chrom == unique(chrom)[j]),]
			plot(cChrom$MAP, cChrom$Fst, col= as.factor(cChrom$locType), xlim= c(0, 25), ylim= c(0,1), main= paste('chrom ', j-1, ' gen = ', i), ylab= 'Fst', xlab= 'map', pch= 20)
			lines(cChrom$MAP, cChrom$Fst, col= 'grey70')
			screen(j+4)
			plot(envelope(cCpp, Kest, correction= 'isotropic'), main= paste('chrom', j-1))
			# 
			screen(j+8)
			barx <- barplot(cMorI$obs, ylim= c(-1, 1), pch= 20, main= paste("Moran's I, gen = ", i), names.arg= c(1,2,3,4,5), ylab= "Moran's I", xlab= 'distance bins')  
			sdI <- cMorI$sd
			sdI[which(sdI == Inf)] <- NA
			sdI[which(sdI == NaN)] <- NA
			error.bar(barx, cMorI$obs, sdI)
		}
		Sys.sleep(wait)
		if (static == T) {
				break
		}
	}
}


plotMorIripK <- function(mIKgenAll, time= 1) {
	# function to plot the single Moran's I values per generation and chromosome 
	# INPUT: output of calcMorI() and the generation time to start the loop (with nGen/tsFreq)
	close.screen(all.screens= T)

	par(bg = "white") # erase.screen() will appear not to work if the background color is transparent 
	split.screen(c(2,4))  


	for (i in time:length(mIgenAll)) {
		cMorIall <- mIKgenAll[[1]][[i]]
		ripKall <- mIKgenAll[[2]][[i]]
		barx <- barplot(cMorIall$obs, ylim= c(-1, 1), pch= 20, main= paste('Morans I; gen = ', i), names.arg= c(0,1,2,3), xlab= 'chromosome', ylab= "Moran's I")
		#lines(cMorIall$obs, col= 'grey70')
		error.bar(barx, cMorIall$obs, cMorIall$sd)
		Sys.sleep(0.6)
	}
}





# plotFstMorI


#############################################





	#lapply(fst, subset, allele_frequencies >= 25e-4)

#function(x) length(x[x<0])
#####################################################
# deprecated   (func to plot afDiffs per gen)


plotDynaMorIgen <- function(ccObj, morIdata, run, time= 1, static= F, wait= 0, ...) {
	# function to plot allele frequencies for selected and neutral sites per generation and Moran's I values (calcMorIripK[[1]]) per generation and chromosome, as well as LD and avg Fst values (for S & N sites)
	# INPUT: output of ccObj, calcMorI() and the generation time to start the loop (with nGen/tsFreq). static = T to stop after generation i (time). 

	close.screen(all.screens= T)
	#par(oma=c(4.5,4.5,4.5,4.5), mar=c(4.0,4.0,4.0,4.0) + 0.2, mgp= c(3,1,0))
	par(bg = "white") # erase.screen() will appear not to work if the background color is transparent 
	split.screen(c(2,3))  

	mainList <- c("neutral sites", "selected sites")
	#screen(6) #, new= F)
	#plot(unlist(ccObj$phiObs), xlim=  c(0, length(ccObj$phiObs)+2), ylim= c(0,10005), ylab= expression(paste(phi, ' obs.')), main= expression(paste(phi, ' obs.')), xlab= '', type= 'l')
	#
	#screen(4) 


	screen(2)
	cols <- c('blue', 'red', 'black') #t(rainbow(3))
	plot(ccObj$FSTs$FSTneut, type= 'l', col= cols[1], ylim= c(0,1), xlim= c(0, length(ccObj$phiObs)+2), ylab= expression('avg F'[st]), main= expression('avg F'[st]))
	points(ccObj$FSTs$FSTsel, col= cols[2], type= 'l')
	points(ccObj$FSTs$FSTtot, col= cols[3], type= 'l', lty= 1)
	legend('topleft', legend= c('neutral', 'selected', 'total'), fill= cols)
	#
	LDlist <- list(ccObj$LDneut, ccObj$LDsel)
	#for (k in 1:2) {
	#	screen(k+2) #, new= F)
	#	plot(1:length(ccObj$phiObs), LDlist[[k]][,2], xlim= c(0, length(ccObj$phiObs)+2), ylim= c(0,1), main= paste('LD', mainList[k]), ylab= 'Avg LD', xlab= '', type= 'l')
	#}
	screen(3)
	plot(1:length(ccObj$phiObs), LDlist[[1]][,2], xlim= c(0, length(ccObj$phiObs)+2), ylim= c(0,1), main= 'LD', ylab= 'Avg LD', xlab= '', type= 'l', col= 'blue')
	points(1:length(ccObj$phiObs), LDlist[[2]][,2], col= 'red', type= 'l')
	legend('topleft', legend= c('neutral', 'selected'), fill= c('blue', 'red'))
	
	screen(1)
	cWallS <- lapply(ccObj$cWallS, unlist)
	phiOncW <- mapply(rep, ccObj$phiObs, times= unlist(lapply(cWallS, length)))
	#yLim <- max(unlist(ccObj$clineWidthSmax))
	plot(log10(ccObj$kphisMax), ccObj$pHatsMax, type= 'l', ylim= c(-0.45,  1), xlim= c(-2, max(log10(ccObj$kphisMax))), xlab= expression(paste('log'[10], ' ', phi)), ylab= expression(paste('p'[i0]~'- p'[i1])), col= 'black')
	#abline(h= 0, lty= 3)
	#}
	points(log10(unlist(phiOncW)), unlist(cWallS), pch= '.', col= 'grey70') #type= 'l')
	points(log10(unlist(ccObj$phiObs)), unlist(lapply(lapply(ccObj$afDiffS, abs), mean)), pch= '.', cex= 1.4, col= 'red')
	points(log10(unlist(ccObj$phiObs)), unlist(lapply(lapply(ccObj$afDiffN, abs), mean)), pch= '.', cex= 1.4, col= 'blue')
	#points(log10(unlist(ccObj$phiObs)), unlist(ccObj$avgAFDiff), type= 'l', col= 'black')

	legend('bottomright', legend= c(expression(paste(phi[Kruuk], ' ~ ', 'peq '[sMax])), expression(paste(phi, ' ~ ', bar(p), ' all s')), expression(paste(phi, ' ~ avg p S')), expression(paste(phi, ' ~ avg p N'))), fill= c('black', 'grey70', 'red', 'blue'), cex= 0.8)

	for (i in time:length(ccObj$phiObs)) {									# 1174 
		afDifflist <- list(ccObj$afDiffN[[i]], ccObj$afDiffS[[i]])
		#LDlist <- list(ccObj$LDneut[1:i,], ccObj$LDsel[1:i,])
		for (j in 1:2){
			screen(j+3)
			hist(afDifflist[[j]], xlim= c(-1,1), ylim= c(0, 1000), main= paste('afDiffs', mainList[j]), xlab= 'allele freq diffs.', ylab= 'prop. of sites')
			if (j == 1) {
				text(0.5, 1000, paste('nGen=', i*1000))
			}
			screen(6)
			cMorIall <- morIdata[[i]]		# mIgen
			barx <- barplot(cMorIall$obs, ylim= c(-1, 1), pch= 20, main= paste('Morans I; gen = ', i), names.arg= c(0,1,2,3), xlab= 'chromosome', ylab= "Moran's I")
			error.bar(barx, cMorIall$obs, cMorIall$sd)

			#screen(j+2) #, new= F)
			#plot(1:i, LDlist[[j]][,2], xlim= c(0, length(ccObj$phiObs)+2), ylim= c(0,1), main= paste('LD', mainList[j]), ylab= 'Avg LD', xlab= '', type= 'l')
			
	
			#plot(LDlist[[j]]$V1, LDlist[[j]]$V4, xlim= c(0, 1174), y= c(0,1), type = 'l')
		}
				#screen(6) #, new= F)
		#plot(1:i, unlist(lapply(fstspl, mean))[1:i], xlim=  c(0, length(ccObj$phiObs)+2), ylim= c(0,1), ylab= expression('avg F'[st]), main= expression('avg F'[st]), xlab= '', type= 'l') 
		#close.screen(all= T)
		title(run, outer= T, line= -1)		# paste(run, ' (',i, ')', sep= '')
		if (static == T) {
				break
			}
		Sys.sleep(wait)
	}
}



plotFstMAPmorIgen <- function(fstspl, morIgen, time= 1, static= F, wait= 0.1, ...) {
	# function to plot both the MAP vs Fst and Moran's I and Moran's I per generation and chromosome 
	# INPUT: fst data split by nGen, nested list of Moran's I values (per gen, chrom and bin) and the generation time to start the loop (with nGen/tsFreq).
	close.screen(all.screens= T)

	par(bg = "white") # erase.screen() will appear not to work if the background color is transparent 
	split.screen(c(2,1))  # split into two rows, one column  (screens = 1,2)
	split.screen(c(1,4), screen= 1)   # split top row (screen 1) into 4 columns (screens = [3,4,5,6])
	#par(mfrow= c(2,4))
	for(i in time:length(fstspl)) {
		chrom <- fstspl[[i]]$chromosomeMembership
		for (j in 1:length(unique(chrom))) {
			screen(j+2)
			cChrom <- fstspl[[i]][which(chrom == unique(chrom)[j]),]
			plot(cChrom$MAP, cChrom$Fst, col= as.factor(cChrom$locType), xlim= c(0, 25), ylim= c(0,1), main= paste('chrom ', j-1, ' gen = ', i), ylab= 'Fst', xlab= 'map', pch= 20)
			lines(cChrom$MAP, cChrom$Fst, col= 'grey70')
		
		#for (j in 1:length(morIbin[[i]])) {
		}
		screen(2)
			cMorI <- morIgen[[1]][[i]]#[[j]]
			barx <- barplot(cMorI$obs, ylim= c(-1, 1), pch= 20, main= "Moran's I", names.arg= 0:3, ylab= "Moran's I", xlab= 'chrom')      # col= grouplist[[i]][[j]]
			sdI <- cMorI$sd
			sdI[which(sdI == Inf)] <- NA
			sdI[which(sdI == NaN)] <- NA
			error.bar(barx, cMorI$obs, sdI)

		if (static == T) {
				break
			}
		Sys.sleep(wait)

	}
}



plotStatic.1 <- function(ccObj, run, data, ...) {
	# function to create the per-run plots with LD, Fst, PHIs, Le and me
	close.screen(all.screens= T)
	par(oma=c(4.5,4.5,4.5,4.5), mar=c(4.0,4.0,4.0,4.0) + 0.2, mgp= c(3,1,0))
	#par(bg = "white") # erase.screen() will appear not to work if the background color is transparent 
	split.screen(c(3,2))  
	#mtext(run, outer = TRUE )
	mainList <- c("neutral sites", "selected sites")
	LDlist <- list(ccObj$LDneut, ccObj$LDsel)
	param <- data[which(data$run == run),]
	m <- param$sd_move
	s <- param$mean_s
	mutdist <- param$mutation_distribution  #[which(data$run == run)]
	tsFreq <- param$ts_sampling_frequency
	xLab <- paste('sampling times (tsFreq =', tsFreq, ')')
	#
	# LD
	screen(1)
	plot(1:length(ccObj$phiObs), LDlist[[1]][,2], xlim= c(0, length(ccObj$phiObs)+2), ylim= c(0, 1), main= 'LD', ylab= 'Avg LD', xlab= xLab, type= 'l', col= 'blue')
	points(1:length(ccObj$phiObs), LDlist[[2]][,2], col= 'red', type= 'l')
	legend('bottomright', legend= c('neutral', 'selected'), fill= c('blue', 'red'))
	#
	# Fst 
	screen(2)
	cols <- c('blue', 'red', 'black') #t(rainbow(3))
	plot(ccObj$FSTs$FSTneut, type= 'l', col= cols[1], ylim= c(0,1), xlim= c(0, length(ccObj$phiObs)+2), ylab= expression('avg F'[st]), main= expression('avg F'[st]), xlab= xLab)
	points(ccObj$FSTs$FSTsel, col= cols[2], type= 'l')
	points(ccObj$FSTs$FSTtot, col= cols[3], type= 'l', lty= 1)
	legend('bottomright', legend= c('neutral', 'selected', 'total'), fill= cols, cex= 0.8)
	#
	# PHIs per generation
	screen(3) #, new= F)
	plot(log10(ccObj$phiObs), ylab= expression(paste('log'[10]~ phi)), xlab= xLab, type= 'l', ylim= c(0, max(log10(ccObj$kruukPhi_sMax))+ 1), col= 'grey70')  # xlim=  c(0, length(ccObj$phiObs)+100)
	points(log10(ccObj$kruukPhi_sMax), type= 'l', col= 'black')
	legend('topleft', legend= c(expression(paste(phi[Kruuk])), expression(paste(phi, ' ', bar(s)))), fill= c('black', 'grey70'), cex= 0.8)
	text(x= 0.8*length(ccObj$meanS), y= 0.6*max(log10(ccObj$kruukPhi_sMax)), paste('s = ', s))
	text(x= 0.8*length(ccObj$meanS), y= 0.5*max(log10(ccObj$kruukPhi_sMax)), paste('m = ', m))

	#
	# cline widths and PHIs
	screen(4)
	cWallS <- lapply(ccObj$clineWallS, unlist)
	phiOncW <- mapply(rep, ccObj$phiObs, times= unlist(lapply(cWallS, length)))
	#yLim <- max(unlist(ccObj$clineWidthSmax))
	plot(log10(ccObj$kruukPhi_sMax), ccObj$pHat, type= 'l', ylim= c(-0.45,  1), xlim= c(0, max(log10(ccObj$kruukPhi_sMax))), xlab= expression(paste('log'[10], ' ', phi)), ylab= expression(paste('p'[i0]~'- p'[i1])), col= 'black')
	abline(h= 0, lty= 3)
	#}
	points(log10(unlist(phiOncW)), unlist(cWallS), pch= '.', col= 'grey70') #type= 'l')
	points(log10(unlist(ccObj$phiObs)), unlist(lapply(lapply(ccObj$afDiffS, abs), mean)), type= 'l', col= 'red')
	points(log10(unlist(ccObj$phiObs)), unlist(lapply(lapply(ccObj$afDiffN, abs), mean)), type= 'l', col= 'blue')
	#points(log10(unlist(ccObj$phiObs)), unlist(ccObj$avgAFDiff), type= 'l', col= 'black')

	legend('bottomright', legend= c(expression(paste(phi[Kruuk], ' ~ ', 'peq '[sMax])), expression(paste(phi, ' ~ pBar all s')), expression(paste(phi, ' ~ avg p S')), expression(paste(phi, ' ~ avg p N'))), fill= c('black', 'grey70', 'red', 'blue'), cex= 0.8)
	#
	# plot sStar, mean_s, and sBar
	screen(5)
	#plot(1, type="n", xlab= 'gen / 1e3', ylab="", ylim=c(-0.1, 0.1), xlim=c(0, (ccObj$sStarLeS)))
	emig <- ccObj$effMig[,2]
	emig[which(emig == 0)] <- 1e-10
	plot(log10(abs(ccObj$sStarLeS$Le)), type= 'l', xlab= xLab, ylab= expression(paste('log'[10]~'L'[e])) , ylim= c(-10, 3))#, ylim= c(min(log10(emig)), max(log10(abs(ccObj$sStarLeS$Le)))))
	points(log10(emig), col= 'blue', type= 'l')
	abline(v= ccObj$gwcTimeMeanS$gwcTime/tsFreq, lty= 2, col= 'green')
	legend('bottomright', legend= c(expression('L'[e]), expression('m'[e]), 'gwcTime meanS'), lty= c(1, 1, 2), col= c('black', 'blue', 'green'), cex= 0.8)
	text(x= 0.6*length(ccObj$meanS), y= 0.6*max(log10(emig)), paste("gwc time = ", ccObj$gwcTimeMeanS$gwcTime))
	#
	# plot effMig  
	screen(6)
	plot(ccObj$maxEffMigSbar, col= 'orange', type= 'l', ylim= c(0, 1.6*max(ccObj$effMig[,2])), ylab= 'migration rate', xlab= xLab)
	points(ccObj$effMig[,2], col= 'black', type= 'l')
	points(ccObj$maxEffMigMeanS, col= 'green', type= 'l')

	abline(h= m, col= 'red', lty= 3)
	abline(h= s, col= 'blue', lty= 3)
	abline(v= ccObj$gwcTimeSbar$gwcTime/tsFreq, col= 'orange', lty= 2)
	abline(v= ccObj$gwcTimeMeanS$gwcTime/tsFreq, col= 'green', lty= 2)
	endAllo <- data$end_period_allopatry[which(data$run == run)]
	if (ccObj$gwcTimeMeanS$gwcTime == endAllo) {
		text(x= 0.3* length(ccObj$meanS), y= 1.1*max(ccObj$effMig[,2]), expression(paste('gwcTime == endPallo')))
	}
	text(x= 0.8* length(ccObj$meanS), y= 0.4*max(ccObj$effMig[,2]), paste('mutdist = ', mutdist)) 
	text(x= 0.3* length(ccObj$meanS), y= 1.3*max(ccObj$effMig[,2]), paste('endPallo = ', param$end_period_allopatry))
	cols <- c('orange', 'green', 'black', 'red', 'blue', 'orange', 'green')
	legend('topright', legend= c(expression(paste('maxEffMig ', bar(s))), 'maxEffMig meanS', expression(paste('m'[e0])), 'm', 's', expression(paste('gwcTime ', bar(s))), 'gwcTime meanS'), col= cols, lty= c(1,1,1,3,3,2,2), pch= c(NA,NA,NA,NA,NA,NA,NA), pt.cex= 2, seg.len= 0.9, cex= 0.7)

	#
	close.screen(all= T)
	title(run, outer= T)
}

