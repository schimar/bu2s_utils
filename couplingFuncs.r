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

# function to calculate the expected cline width under full coupling (with sMax == is.bool)
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
	return(unlist(list(sStar, Le, s)))
}


# calculate average Fst values per generation and locType (mean selected, mean neutral, mean total)
getFSTs <- function(fstSubGen, i) {
	meanSel <- mean(fstSubGen[[i]]$Fst[which(fstSubGen[[i]]$locType == 1)])
	meanNeut <- mean(fstSubGen[[i]]$Fst[which(fstSubGen[[i]]$locType == 0)])
	total <- mean(fstSubGen[[i]]$Fst)
	return(list(meanNeut, meanSel, total))
}



# calculate coupling/congealing stats (GWC time still pending)

coupCong <- function(df, fst, afts, LDsel, LDneut, run, maf= 25e-4, s= 1, nChrom= 4) {   # randProploci = T,...){
	params <- df[which(df$run == run),]
	m = params$sd_move
	fstSpl <- split(fst, fst$nGen)
	aftsSpl <- split(afts, afts$nGen)
	gen <- length(aftsSpl)
	#
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
	#
	for (i in 1:gen){
		mafThresh <- which(fstSpl[[i]]$allele_frequencies >= maf)
		fstCurGen <- fstSpl[[i]][mafThresh,]
		afCurGen <- aftsSpl[[i]][mafThresh,]
		#nGen <- fstGen$nGen
		locType <- which(fstCurGen$locType == s)
		p <- fstCurGen$allele_frequencies
		sSel <- fstCurGen$S_MAX1[locType]
		#map <- fstCurGen$MAP
		#mapS <- map[locType]
	#
		### if (randProploci == T) {
		#locID <- afCurGen$locusID
		#pPatches <- aftsSpl[[i]][,3:4]
		#aftslocType <- aftsSpl[[i]]$locType
		afDiffS[[i]] <- afCurGen$AFdiff[locType]
		afDiffN[[i]] <- afCurGen$AFdiff[-locType]
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

		# they're lumped together (sStar and Le)
		# we need sSel as well, I think
		clineWidthAllS <- calcClineWidth(sSel, m, sMax= F)
		pBarAllS[[i]] <- clineWidthAllS[1]
		clineWallS[[i]] <- clineWidthAllS[2]
		#pBarAvg[[i]] <- #  "avgAFdiff in the L counted loci" (= obs. freq. of avg allele in the favored deme, pBar 
			# pBar goes into calcLe, and m (from df) as well as s (sSel?)  
		Le <- calcLe(m, unlist(pBarAllS[[i]]), sSel)
		sStLeS[[i]] <- as.data.frame(matrix(Le, ncol= 3, dimnames= list(NULL, c('sStar', 'Le', 's'))))
		# they're both 
		FSTs[[i]] <- getFSTs(fstSpl, i)	
		
		#bartonTheta[[i]] <- sum(si[[i]])/rVec[[i]]

	   	# if all loci have constant s (sMax == s); expected value (~) 
		#PhiObs <- sum(fstSub$S_MAX1[Sloci]) / rV 
	
	}
	LDsell <- LDsel[, c(1,4)]
	LDneutl <- LDneut[, c(1,4)]
	FSTout <- as.data.frame(matrix(unlist(FSTs),ncol=3,byrow=TRUE, dimnames= list(NULL, c('FSTneut', 'FSTsel', 'FSTtot'))))
	runOut <- list(afDiffS, afDiffN, FSTout, recom, sMax, kruukPhi_sMax, phiObs, phiTilde, sBar, pHat, clineWsMax, pBarAllS, clineWallS, sStLeS, LDneutl, LDsell)
	names(runOut) <- c('afDiffS', 'afDiffN', 'FSTs', 'recomb', 'sMax', 'kruukPhi_sMax', 'phiObs', 'phiTil', 'sBar', 'pHat', 'clineWidthSmax', 'pBarAllS', 'clineWallS','sStarLeS', 'LDneut', 'LDsel')
	return(runOut)
}


### 

# function to read a single data set (as input for coupCong)

readCCobj <- function(run, smParam, path, ...) {
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
	out <- list(fstTmp, aftsTmp, LDselTmp, LDneuTmp)
	names(out) <- c('fst', 'afts', 'LDsel', 'LDneut')
	return(out)
}



#####################################################
# plot functions


plotGen <- function(ccObj, ...) {

	close.screen(all.screens= T)

	par(bg = "white") # erase.screen() will appear not to work if the background color is transparent 
	split.screen(c(3,2))  

	mainList <- c("neutral sites", "selected sites")
	screen(6) #, new= F)
	plot(unlist(ccObj$phiObs), xlim=  c(0, length(ccObj$phiObs)+2), ylim= c(0,10005), ylab= expression(paste(phi, ' obs.')), main= expression(paste(phi, ' obs.')), xlab= '', type= 'l')
	#
	screen(5)
	cols <- c('blue', 'red', 'black') #t(rainbow(3))
	plot(ccObj$FSTs$FSTneut, type= 'l', col= cols[1], ylim= c(0,1), xlim= c(0, length(ccObj$phiObs)+2), ylab= expression('avg F'[st]), main= expression('avg F'[st]))
	points(ccObj$FSTs$FSTsel, col= cols[2], type= 'l')
	points(ccObj$FSTs$FSTtot, col= cols[3], type= 'l', lty= 1)
	legend(200, 1, legend= c('neutral', 'selected', 'total'), fill= cols)
	#
	LDlist <- list(ccObj$LDneut, ccObj$LDsel)
	#for (k in 1:2) {
	#	screen(k+2) #, new= F)
	#	plot(1:length(ccObj$phiObs), LDlist[[k]][,2], xlim= c(0, length(ccObj$phiObs)+2), ylim= c(0,1), main= paste('LD', mainList[k]), ylab= 'Avg LD', xlab= '', type= 'l')
	#}
	screen(3)
	plot(1:length(ccObj$phiObs), LDlist[[1]][,2], xlim= c(0, length(ccObj$phiObs)+2), ylim= c(0,1), main= 'LD', ylab= 'Avg LD', xlab= '', type= 'l', col= 'blue')
	points(1:length(ccObj$phiObs), LDlist[[2]][,2], col= 'red', type= 'l')
	legend(200, 1, legend= c('neutral', 'selected'), fill= c('blue', 'red'))

	for (i in 1:length(ccObj$phiObs)) {									# 1174 
		afDifflist <- list(ccObj$afDiffN[[i]], ccObj$afDiffS[[i]])
		#LDlist <- list(ccObj$LDneut[1:i,], ccObj$LDsel[1:i,])
		for (j in 1:2){
			screen(j)
			hist(afDifflist[[j]], xlim= c(-1,1), ylim= c(0, 1000), main= paste('afDiffs', mainList[j]), xlab= 'allele freq diffs.', ylab= 'prop. of sites')
			if (j == 1) {
				text(0.5, 1000, paste('nGen=', i*1000))
			}
			#screen(j+2) #, new= F)
			#plot(1:i, LDlist[[j]][,2], xlim= c(0, length(ccObj$phiObs)+2), ylim= c(0,1), main= paste('LD', mainList[j]), ylab= 'Avg LD', xlab= '', type= 'l')
			
	
			#plot(LDlist[[j]]$V1, LDlist[[j]]$V4, xlim= c(0, 1174), y= c(0,1), type = 'l')
		}
				#screen(6) #, new= F)
		#plot(1:i, unlist(lapply(fstspl, mean))[1:i], xlim=  c(0, length(ccObj$phiObs)+2), ylim= c(0,1), ylab= expression('avg F'[st]), main= expression('avg F'[st]), xlab= '', type= 'l') 
	}
}


plotStatic <- function(ccObj, run,...) {

	close.screen(all.screens= T)
	par(oma=c(0,0,2,0))
	#par(bg = "white") # erase.screen() will appear not to work if the background color is transparent 
	split.screen(c(3,2))  
	#mtext(run, outer = TRUE )
	mainList <- c("neutral sites", "selected sites")
	LDlist <- list(ccObj$LDneut, ccObj$LDsel)
	#
	# LD
	screen(1)
	plot(1:length(ccObj$phiObs), LDlist[[1]][,2], xlim= c(0, length(ccObj$phiObs)+2), ylim= c(0, 1), main= 'LD', ylab= 'Avg LD', xlab= '', type= 'l', col= 'blue')
	points(1:length(ccObj$phiObs), LDlist[[2]][,2], col= 'red', type= 'l')
	legend('bottomright', legend= c('neutral', 'selected'), fill= c('blue', 'red'))
	#
	# Fst 
	screen(2)
	cols <- c('blue', 'red', 'black') #t(rainbow(3))
	plot(ccObj$FSTs$FSTneut, type= 'l', col= cols[1], ylim= c(0,1), xlim= c(0, length(ccObj$phiObs)+2), ylab= expression('avg F'[st]), main= expression('avg F'[st]))
	points(ccObj$FSTs$FSTsel, col= cols[2], type= 'l')
	points(ccObj$FSTs$FSTtot, col= cols[3], type= 'l', lty= 1)
	legend('bottomright', legend= c('neutral', 'selected', 'total'), fill= cols)
	#
	# PHIs per generation
	screen(3) #, new= F)
	plot(unlist(ccObj$phiObs), xlim=  c(0, length(ccObj$phiObs)+2), ylim= c(0, max(unlist(ccObj$phiObs))+ 200), ylab= expression(paste(phi, ' obs.')), main= expression(paste(phi, ' obs.')), xlab= '', type= 'l')
	points(unlist(ccObj$kruukPhi_sMax), type= 'l', col= 'grey70')
	legend('top', legend= c(expression(phi), expression(paste("Kruuk's ", phi))), fill= c('black', 'grey70'))
	#
	# cline widths and PHIs
	screen(4)
	cWallS <- lapply(ccObj$clineWallS, unlist)
	phiOncW <- mapply(rep, ccObj$phiObs, times= unlist(lapply(cWallS, length)))
	#yLim <- max(unlist(ccObj$clineWidthSmax))
	plot(log10(unlist(ccObj$kruukPhi_sMax)), unlist(ccObj$clineWidthSmax), type= 'l', ylim= c(-0.45,  1), xlab= expression(paste('log'[10], ' ', phi)), ylab= expression(paste('p'[i0]~'- p'[i1])), col= 'grey70')
	#}
	points(log10(unlist(phiOncW)), unlist(cWallS), pch= '.') #type= 'l')
	#
	# plot sStar, mean_s, and sBar
	screen(5)
	plot(1, type="n", xlab="", ylab="", xlim=c(0, length(ccObj$sStarLeS)+100), ylim=c(-0.05, 0.1))
	for (i in length(ccObj$sStarLeS)){
		points(ccObj$sStarLeS[[i]]$sStar, col= 'green', pch= 20, cex= 0.4) #, type= 'l')
		points(ccObj$sStarLeS[[i]]$s, col= 'grey50', pch= 20, cex= 0.4) # type= 'l')
		points(unlist(ccObj$sBar), col= 'orange', pch= 20, cex= 0.4)
	}
	legend('topleft', legend= c('mean s', expression(bar(s)), expression('s'^'*')), fill= c('grey50', 'orange', 'green'))

	#
	close.screen(all= T)
	title(run, outer= T)
	

}


# function to read individual runs (from vector of runs), calculate CC and plotStatic
wrapH5 <- function(path, data, setname, ...) {
	#pathTmp <- '/media/schimar/schimar2/bu2s/h5/'
	#set <- 'Sm'

	for (i in 1:dim(data)[1]){
		run <- data$run[i]
		path5 <- paste('/runs/', run, sep= '')
		#
		ccObjTmp <- readCCobj(run, setname, path)
		ccTmp <- coupCong(df, ccObjTmp$fst, ccObjTmp$afts, ccObjTmp$LDsel, ccObjTmp$LDneut, run)
		#
		plotStatic(ccTmp, run)
		#plotGen(ccTmp)
		#Sys.sleep(1.5)
		H5close()
	}
}


