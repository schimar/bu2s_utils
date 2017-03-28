# functions to calculate recombination and coupling statistics 

library(rhdf5)
library(scales)


# Barton 1983's coupling coefficient (theta = s/r)
# Kruuk 1999's summed coupling coefficient (phi = (L-1)s/r) 
# Le (effective nLoci)  (Le = s*/s)


# calc r

calcRecomb <- function(nLoci, map, nChrom= 4,...){
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

# function to calculate the expected cline width under full coupling (with sMax)
calcClineWidth <- function(s, m, sMax= T){
	pHat <- (m * (0.5 + 0.75* s) - 0.125* s - (0.125* sqrt(s^2 - 4* m* s^2 + (4* m^2) * (2 + s)^2))) / (-1* (0.25 + m)*s )
	clineWidth <- (pHat - (1 - pHat))
	if (sMax == T) {
		return(unlist(list(pHat, clineWidth)))
	}
	else {
		out <- list()
		out[[1]] <- pHat
		out[[2]] <- clineWidth
		return(out)
	}
}

# calculate effective number of loci (L_e), with pBar = afDiffs
calcLe <- function(m, pBar, s) {
	sStar <- (m * (pBar - 0.5)) / ((0.25 - 0.25* pBar)* pBar + m * (0.5 + pBar * (pBar - 1.5)))
	Le <- sStar / s
	return(unlist(list(sStar, Le)))
}


# calculate coupling/congealing stats (GWC time still pending)

coupCong <- function(df, fst, afts, run, maf= 25e-4, s= 1, nChrom= 4,...){
	#mafT <- which(afts$AF >= maf)moab
	params <- df[which(df$run == run),]
	m = params$sd_move
	fstSpl <- split(fst, fst$nGen)
	aftsSpl <- split(afts, afts$nGen)
	#
	#locType <- which(fstSpl[[i]]$locType == s)
	gen <- length(aftsSpl)
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
	#
	for (i in 1:gen){
		mafThresh <- which(fstSpl[[i]]$allele_frequencies >= maf)
		fstCurGen <- fstSpl[[i]][mafThresh,]
		afCurGen <- aftsSpl[[i]][mafThresh,]
		#nGen <- fstGen$nGen
		locType <- which(fstCurGen$locType == s)
		#locType <- fstSpl[[i]]$locType
		p <- fstCurGen$allele_frequencies
		sSel <- fstCurGen$S_MAX1[locType]
		#map <- fstCurGen$MAP
		#mapS <- map[locType]
	#
		#locID <- afCurGen$locusID
		#pPatches <- aftsSpl[[i]][,3:4]
		#aftslocType <- aftsSpl[[i]]$locType
		afDiffS[[i]] <- abs(afCurGen$AFdiff[locType])
		afDiffN[[i]] <- abs(afCurGen$AFdiff[-locType])
	#
		#nChrom <- df$nchromosomes[which(df$run == run)]   
		nLoci <- dim(fstCurGen)[1]     
		nLociS <- length(locType)
		recom[[i]] <- calcRecomb(nLoci, 0.02*params$total_map_length, nChrom)  
		sMax[[i]] <- prod(1 + sSel) - 1   
		kruukPhi_sMax[[i]] <- (nLociS - 1) * (sMax[[i]]/recom[[i]])     # better mean or max of rV ???
		phiObs[[i]] <- (sum(sSel) / recom[[i]])
		phiTilde[[i]] <- phiObs[[i]] / nLociS
		sBar[[i]] <- (phiObs[[i]] / nLociS) * recom[[i]]
		clineWidthSmax <- calcClineWidth(sMax[[i]], m, sMax= T)
		pHat[[i]] <- clineWidthSmax[1]
		clineWsMax[[i]] <- clineWidthSmax[2]
		#
		clineWidthAllS <- calcClineWidth(sSel, m, sMax= F)
		pBarAllS[[i]] <- clineWidthAllS[1]
		clineWallS[[i]] <- clineWidthAllS[2]
		
		# now calculate all pBar (per generation and locus), to get the full spread of allele frequencies (pBar) 
		# this will inform us about the expectations under no coupling 
		
		
		
		#sMax[[i]] <- prod(1 + aftsSpl[[i]]$AF) - 1 
		#bartonTheta[[i]] <- sum(si[[i]])/rVec[[i]]
	 

		#kruukPhi_sMax <- (nLociS - 1) * (sMax[[i]]/rV)   # if all loci have constant s (sMax == s); expected value (~) 

		#PhiObs <- sum(fstSub$S_MAX1[Sloci]) / rV 





				#kruukPhi[[i]] <- (nLociS - 1) * bartonTheta[[i]]
	}
	runOut <- list(afDiffS, afDiffN, recom, sMax, kruukPhi_sMax, phiObs, phiTilde, sBar, pHat, clineWsMax, pBarAllS, clineWallS)
	names(runOut) <- c('afDiffS', 'afDiffN', 'recomb', 'sMax', 'kruukPhi_sMax', 'phiObs', 'phiTil', 'sBar', 'pHat', 'clineWidthSmax', 'pBarAllS', 'clineWallS')
	return(runOut)
}



#####################################################
# plot function, that takes coupCong-list (cocol) as input

plotOcoco <- function(cocol, LD= LD, diff= dXY, ...) {
	par(mfrow= c(2,4))#, mar=c(1.5,2,1,1) + 0.1, oma= c(5,0,0,0), mgp= c(0,1,0))
	# plot afDiffs for S & N sites
	# LD with nearest selected site (for S & N sites)
	# log10 (N* m_e)  and time 
	# dXY by time (for all, selected & neutral sites)

	# coupling: 
	# | p_i1 - p_i0 |    by phi  (exp. no-coupling, exp. full coupling and observed) 
	
	
	plot(x)
}


#plotOstats <- function(LD 
