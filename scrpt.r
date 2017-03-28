
##################################################################################################
##################################################################################################
library(rhdf5) 
library(scales)

####### 
#  coupling 

# Barton 1983's coupling coefficient (theta = s/r)
# Kruuk 1999's summed coupling coefficient (phi = (L-1)s/r) 
# Le (effective nLoci)  (Le = s*/s)

source("~/flaxmans/bu2s/bu2s_utils/couplingFuncs.r")

df <- read.table("~/flaxmans/bu2s/runs/paramsALL.txt", header= T, sep= '\t')
names(df) <- tolower(names(df))




# example run (215499) 

fst <- read.table("xRuns/Run215499/FSTtimeSeries.txt.bz2")
colnames(fst) <- c("nGen", "locusID", "Fst", "allele_frequencies", "S_MAX1", "S_MAX0", "chromosomeMembership", "MAP", "locType")


#plot(fst$nGen, fst$Fst, col= as.factor(fst$locType))


# afts

afts <- read.table("xRuns/Run215499/AlleleFreqTimeSeries.txt.bz2")
colnames(afts) <- c("nGen", "locusID", "AFpatch0", "AFpatch1", "is_reversed_locus", "locType", "AF", "AFdiff")


#plot(afts$nGen, afts$AFdiff, ylab= 'afDiff', xlab= 'nGen', pch= '.', cex= 2.5, col= alpha(afts$locType +1, 0.7))

r499 <- coupCong(df, fst, afts, 'Run215499')



###########################  hdf5  loop



# loop through the respective parameter file (df  ->  Sm, sm, sM)
paramSub <- Sm



for (i in 1:dim(paramSub)[1]){
	run <- paramSub$run[i]
	path <- paste('/runs/', run, sep= '')
	fstTmp <- h5read('/media/schimar/schimar2/bu2s/h5/Sm/Fst_SmT.h5', name= path)[[1]]
	aftsTmp <- h5read('/media/schimar/schimar2/bu2s/h5/Sm/Fst_SmT.h5', name= path)[[1]]
	#
	LDselTmp <- h5read('/media/schimar/schimar2/bu2s/h5/Sm/LDselAvg_SmT.h5', name= path)[[1]]
	LDneuTmp <- h5read('/media/schimar/schimar2/bu2s/h5/Sm/LDneutAvg_SmT.h5', name= path)[[1]]


	colnames(fstTmp) <- c("nGen", "locusID", "Fst", "allele_frequencies", "S_MAX1", "S_MAX0", "chromosomeMembership", "MAP", "locType")
	colnames(aftsTmp) <- c("nGen", "locusID", "AFpatch0", "AFpatch1", "is_reversed_locus", "locType", "AF", "AFdiff")
	# read LDneut & LDsel, dXY (or use Fst), effMig
	
	# calculate coupling stats
	# tmpCoup <- coupCong(df, run, fstTmp, aftsTmp, s= 1)
	
	# then plot stuff (func to be written)
}



	











###############  other files 

LDneut <- read.table('xRuns/Run215499/LDneutSitesAvg.txt.bz2')

LDsel <- read.table('xRuns/Run215499/LDselSitesAvg.txt.bz2')

par(mfrow = c(2,3))
plot(LDneut$V1, LDneut$V2)
plot(LDneut$V1, LDneut$V3)
plot(LDneut$V1, LDneut$V4)
#points(LDneut$V1, LDneut$V2, col= 'lightgreen')
#points(LDneut$V1, LDneut$V3, col= 'yellow')    # move along, nothing to see here...

#
plot(LDsel$V1, LDsel$V2)
plot(LDsel$V1, LDsel$V3)
plot(LDsel$V1, LDsel$V4)


# me
em <- read.table("xRuns/Run215499/EffectiveMigrationRates.txt.bz2")
em <- em[,1:18]
colnames(em)[1:6] <- c("nGen", "eme0", "eme1", "nVariableLoci", "nRes", "nImm")

plot(em$nGen, em$eme0)
points(em$nGen, em$eme1, col= 'orange')




# dXY

dXY <- read.table("xRuns/Run215499/DXYtimeSeries.csv", skip= 1, header= F, sep= ',')

plot(dXY$V1, dXY$V2)
points(dXY$V1, dXY$V3, col= 'lightgreen')
points(dXY$V1, dXY$V4, col= 'yellow')


# AFdiff
neutAF <- read.table("xRuns/Run215499/NeutralAlleleFrequencies.txt.bz2")

selAF <- read.table("xRuns/Run215499/SelectedAlleleFrequencies.txt.bz2")
plot(selAF$V1, selAF$V4)
points(neutAF$V1, neutAF$V4, col= 'red')

################## 


# param (e.g. Sm)    
plotOmat <- function(run, ...){
	par(mfrow= c(2,4))
	plot(LDneut$V1/1e5, LDneut$V4, ylab= 'avg LD neutral', xlab= 'nGen', pch= '.', cex= 2.5)
	plot(neutAF$V1/1e5, neutAF$V7, xlab= 'nGen', ylab= 'neutral - LD w/ nsn', pch= '.', cex= 2.5)
	plot(fst$V1/1e5, fst$V3, col= alpha(fst$V9+1, 0.7), xlab= 'Fst', ylab= 'nGen', pch= '.', cex= 2.5)
	plot(em$nGen/1e5, em$eme0, col= alpha('black', 0.4), xlab= 'eme', ylab= 'nGen', pch= '.', cex= 2.5)
	points(em$nGen/1e5, em$eme1, col= alpha('red', 0.7), pch= '.', cex= 2.5)
	
	
	# 2nd row
	plot(LDsel$V1/1e5, LDsel$V4, ylab= 'avg LD selected', xlab= 'nGen', pch= '.', cex= 2.5)
	plot(selAF$V1/1e5, selAF$V9, xlab= 'nGen', ylab= 'selected - LD w/ nsn', pch= '.', cex= 2.5)
	plot(dXY$V1/1e5, dXY$V2, xlab= 'dXY', ylab= 'nGen', pch= '.', cex= 2.5)
	points(dXY$V1/1e5, dXY$V3, col= 'lightgreen', pch= '.', cex= 2.5)
	points(dXY$V1/1e5, dXY$V4, col= 'yellow', pch= '.', cex= 2.5)
	plot(selAF$V1, selAF$V4, xlab= 'af diff', ylab= 'nGen', col= alpha('red', 0.7), pch= '.', cex= 2.5)
	points(neutAF$V1, neutAF$V4, col= alpha('black', 0.4), pch= '.', cex= 2.5)
	#
	Sys.sleep(2.5)
}




## Run215499
#plotOmat <- function(run, ...){
#	LDneut <- 
#	LDsel <- 
#	neutAF <- 
#	selAF <- 
#	fst <- 
#	em <- 
#	dXY <- dXYSm[[run]]$dXY
#	#
#	par(mfrow= c(2,4))#, mar=c(1.5,2,1,1) + 0.1, oma= c(5,0,0,0), mgp= c(0,1,0))
#	plot(LDneut$V1/1e5, LDneut$V4, ylab= 'avg LD neutral', xlab= 'nGen', pch= '.', cex= 2.5)
#	plot(neutAF$V1/1e5, neutAF$V7, xlab= 'nGen', ylab= 'neutral - LD w/ nsn', pch= '.', cex= 2.5)
#	plot(fst$V1/1e5, fst$V3, col= alpha(fst$V9+1, 0.7), ylab= 'Fst', xlab= 'nGen', pch= '.', cex= 2.5)
#	plot(em$nGen/1e5, em$eme0, col= alpha('black', 0.6), ylab= 'eme', xlab= 'nGen', pch= '.', cex= 2.5)
#	points(em$nGen/1e5, em$eme1, col= alpha('red', 0.7), pch= '.', cex= 2.5)
#	
#	
#	# 2nd row
#	plot(LDsel$V1/1e5, LDsel$V4, ylab= 'avg LD selected', xlab= 'nGen', pch= '.', cex= 2.5)
#	plot(selAF$V1/1e5, selAF$V9, xlab= 'nGen', ylab= 'selected - LD w/ nsn', pch= '.', cex= 2.5)
#	plot(dXY$V1/1e5, dXY$V2, ylab= 'dXY', xlab= 'nGen', pch= '.', cex= 2.5)
#	points(dXY$V1/1e5, dXY$V3, col= 'lightgreen', pch= '.', cex= 2.5)
#	points(dXY$V1/1e5, dXY$V4, col= 'yellow', pch= '.', cex= 2.5)
#	plot(selAF$V1/1e5, selAF$V4, ylab= 'af diff', xlab= 'nGen', col= alpha('red', 0.7), pch= '.', cex= 2.5)
#	points(neutAF$V1/1e5, neutAF$V4, col= alpha('black', 0.5), pch= '.', cex= 2.5)
#	
#	#
#	Sys.sleep(2.5)
#}

######################

# plots

plot(1, type="n", xlab="", ylab="", xlim=c(0, 20000), ylim=c(0, 1))
for (i in 1:1174){
	points(rep(r499$phiObs[[i]], length(r499$afDiffS[[i]])), r499$clineWallS[[i]][[1]])
}
#  expression(p[i0] - p[i1])

#plot(1, type="n", xlab="", ylab="", xlim=c(0, 20000), ylim=c(0, 1))
#for (i in 1:1174) {
# points(unlist(r499$phiObs), unlist(r499$pHat))


#for (i in 1:1174){
#	points(rep(r499$phiObs[[i]], length(r499$afDiffN[[i]])), r499$afDiffN[[i]], col= 'lightgreen')
#}

# 
#r499 <- coupCong(df, fst, afts, 'Run215499')

##
plot(unlist(r499$afDiffs), unlist(r499$clineWidthSmax), ylim= c(0,1))
points(unlist(r499$kruukPhi_sMax), unlist(r499$pHat), col= 'lightgreen')



# sBar vs sMax
plot(unlist(r499$sBar), ylim= c(0, 2.65e+6))
points(unlist(r499$sMax), col= 'lightgreen')

####     sBar vs mean_S  (= 0.005 here)   (especially early, we expect sBar > mean_S)
plot(unlist(r499$sBar), ylim= c(0, 0.03), ylab= expression(bar(s)), xlab= expression(paste("time (generations) x 10"^"3")))
abline(h= 0.005)


# phiObs vs Kruuk's phiSmax
plot(unlist(r499$phiObs), type= 'l', ylim= c(0, 2.6e+12))
points(unlist(r499$kruukPhi_sMax), type= 'l', col= 'cyan')

# a <- h5read('/media/schimar/schimar2/bu2s/h5/Sm/dXY_SmT.h5', name= '/runs/Run202971')
# a[[1]]


plot(unlist(r499$afDiffS[[1]]), unlist(r499$clineWallS[[1]][[1]]), ylim= c(-1,1))

	
	
#geomean = function(x, na.rm=TRUE){
#  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
#}

#	gmean <- geomean(afDiffs)
#		bartonTheta[[i]] <- gmean/recolist[[i]] 
#		kruukPhi[[i]] <- (nLoci - 1)* gmean / recolist[[i]]
#		#effnLoci <- getLe( # m, p, s)
#	}
#	#return(bartonTheta) #recolist)
#	return(kruukPhi)
##return(length(reco))
#}
#cc <- coupCong(fst, afts)


######################
# after chat w/ Sam on 16.03. 
# s[i]   (= prod(1+Si) of all selected loci (w/ MAF >= 25e-4))
#prod(1 + afts1[[1174]]$AF[afts1[[1174]]$AF >= 0.0025])


fst1 <- split(fst, fst$nGen)
afts1 <- split(afts, afts$nGen)

fstSub <- fst1[[1174]]
aftsSub <- afts1[[1174]]

maf <- which(aftsSub$AF >= 0.0025)
fstSub <- fstSub[maf,]
aftsSub <- aftsSub[maf,]

Sloci <- which(fstSub$locType == 1)
mapS <- fstSub$MAP[Sloci]
#lT <- which(fstSub$locType == 1)

rV <- calcRecomb(dim(fstSub)[1], 0.02*100, nChrom= 4)

sMax <- prod(1 + fstSub$S_MAX1[Sloci]) - 1 

nLociS <- length(Sloci)

kruukPhi_sMax <- (nLociS - 1) * sMax/rV   # if all loci have constant s (sMax == s); expected value (~) 

PhiObs <- sum(fstSub$S_MAX1[Sloci]) / rV 

#bartonT <- sum(siSub)/rV

L <- length(rV)
kruukP <- bartonT * (L - 1)




# s[i]
#si <- prod(1 + afts1$AF[maf])

# r
#fstSub <- fst1[maf,]
#rVec <- array()
#for (i in 1:dim(fstSub)[1]){
#	r <- calcRecomb(dim(fstSub)[1], fstSub$MAP[i], nChrom= 4)
#	rVec[i] <- r
#}

#prod(1 + afts1[[1174]]$AF[afts1[[1174]]$AF >= 0.0025])
# observed stats: (for selected and neutral) 

# phi (= (L-1)*(sum(s[i] at selected_loci)/r))

# theta (Barton83) (= sum(s[i] at selected_loci)/r)

# Le 
#effectiveS = FullSimplify[Solve[p1eq == p, s]] (* this is Barton 1983's "s*" *) 
#Plot[((s /. effectiveS) /. m -> 0.05), {p, 0.51, 0.8}]

#p1eq = FullSimplify[(p1 /. (sols[[3]]))]
# mathematica ReplaceAll func:
# {x, x^2, y, z} /. x -> 1
#Plot3D[p1eq, {m, 0.0001, 0.5}, {s, 0.0001, 0.5},  PlotRange -> {0.5, 1},  AxesLabel -> {{"m", "s", "p1"}, FontSize -> 14}]

# expected stats: 

# for (1) coupling and (2) no coupling 

##################
# function for effective nLoci

getLe <- function(m, p, s){
	effS <- (m*(p-0.5)/(0.25-0.25*p)*p + m*(0.5+p*(p-1.5)))
	le <- effS/s
	return(le)
}

getLe(0.01, 0.09, 0.001)




# nearest selected neighbors:
# LD
hist(selAF$V9)
hist(neutAF$V7)

plot(selAF$V1, selAF$V9, xlab= 'nGen', ylab= 'LD w/ nsn')
plot(neutAF$V1, neutAF$V7, xlab= 'nGen', ylab= 'LD w/ nsn')


#par(mfrow= c(2,3))
#	plot(LDneut$V1/1e5, LDneut$V4)
#	plot(fst$V1/1e5, fst$V3, col= as.factor(1-fst$V9))
#	plot(em$nGen/1e5, em$eme0, col= 'red')
#	points(em$nGen/1e5, em$eme1)#, col= '')
	
#	# 2nd row
#	plot(LDsel$V1/1e5, LDsel$V4)
#	plot(dXY$V1/1e5, dXY$V2)
#	points(dXY$V1/1e5, dXY$V3, col= 'lightgreen')
#	points(dXY$V1/1e5, dXY$V4, col= 'yellow')
#	plot(selAF$V1, selAF$V4)
#	points(neutAF$V1, neutAF$V4, col= 'red')


##################################################################################################
##################################################################################################

#df <- read.table("paramsNoOrg_v361.txt", header= T, sep= '\t')


#pH <- read.table("paramHeader.txt", sep= '\t')
#pH <- as.vector(pH[1,1:71])

#colnames(params) <- pH

#########


dfVar <- apply(df, 2, var)

# which(dfVar >0)
dfVar[dfVar != 0]

# subsets for h5 files:
x <- which(df$mean_s < df$sd_move)
# 8784 runs
#write.table(df$run[x], "sM.txt", sep= ' ', row.names= F, quote= F, col.names= F)

#
y <- which(df$mean_s > df$sd_move)
# 1800 runs
#write.table(df$run[y], "Sm.txt", sep= ' ', row.names= F, quote= F, col.names= F)

#
z <- which(df$mean_s == df$sd_move)
#write.table(df$run[z], "sm.txt", sep= ' ', row.names= F, quote= F, col.names= F)




################
# PCA on df
# dfPCA <- prcomp(df, center = T, scale.= F)

################

# subsets for different conditions (s & m) 

# weak s, strong m
which(df$mean_s == 0.005 & df$sd_move == 0.1)

# strong s, weak m
which(df$mean_s == 0.005 & df$sd_move == 0.1)

#  s, strong m
which(df$mean_s == 0.005 & df$sd_move == 0.1)






##################################################
##################################################

# HDF5 


library(rhdf5)

df <- read.table("~/flaxmans/bu2s/runs/paramsALL.txt", header= T, sep= '\t')
names(df) <- tolower(names(df))


########
# s > m 

# params for all 1600 runs
Sm <- df[which(df$sd_move < df$mean_s),]

#tmpFst <- h5ls('/media/schimar/schimar2/bu2s/h5/Sm/Fst_SmT.h5', recursive= T)
#fstSm <- h5read('/media/schimar/schimar2/bu2s/h5/Sm/Fst_SmT.h5', '/runs/')
dXYSm <- h5read('/media/schimar/schimar2/bu2s/h5/Sm/dXY_SmT.h5', name= '/runs/')

#tmpSelAF <- h5ls('SelAF_Sm.h5', recursive= T)

#selAF <- h5read('/media/schimar/dapperdata/bu2s/h5/Sm/SelAF_SmT.h5', '/runs/')



# Fst: 

# we need average Fst values (per nGen) for neutral and selected sites

# (bu2s_2h5fst.py)

#x <- h5read('fstAvgTst.h5', '/runs/')
x <- h5read('/media/schimar/schimar2/bu2s/h5/Sm/fstAvgPerLocType_Sm.h5', '/runs/')



### single run (215499)
x1 <- as.data.frame(x$Run215499$fstAvg,)
colnames(x1) <- c('locType', 'nGen', 'fst')
xx <- split(x1, x1$locType)

plot(xx[[1]]$nGen, xx[[1]]$fst)
points(xx[[2]]$nGen, xx[[2]]$fst, col= 'orange')

# more generic (if run names replaced with index)
plot(x$Run215500$fstAvg[,2], x$Run215500$fstAvg[,3], col= as.factor(x$Run215500$fstAvg[,1]))

#par(mfrow= c(1,2))
plot(1, type="n", xlab="", ylab="", xlim=c(0, 100000), ylim= c(0, 1))
for (i in 1:length(x)){
	points(x[[i]]$fstAvg[,2], x[[i]]$fstAvg[,3], col= as.factor(x[[i]]$fstAvg[,1]))#, add= T)
}


#######


#head(t(selAF[[1]]$SelAF))

plot(1, type="n", xlab="", ylab="", xlim=c(0, 100000), ylim=c(0, 1))
for( i in 1:length(dXYSm)){
	points(dXYSm[[i]]$dXY[,1], dXYSm[[i]]$dXY[,2], pch= '.', cex= 2.5) #, col= as.factor(dXYSm[i]$dXY))
}

order(names(dXYSm))

# subset based on paramsALL.txt 

#plus192 %in% dpl
#plus192[plus192 %in% dpl]


########################################   need to find a neat way to subset the runs with a list (subset) of df$runs (based on exact s&m values 

# figure out how the indices work in the hdf5 file!!! 
# then, you can read in separate chunks (or better: each run individually...) 
# h5f <- h5read('/media/schimar/schimar2/bu2s/h5/Sm/Fst_SmT.h5', name= '/runs/', index= list(1:2))

# 
# fstSm <- h5read('/media/schimar/schimar2/bu2s/h5/Sm/Fst_SmT.h5', name= '/runs/', index= list(1, NULL, NULL))




		
		
		
		
# write function !!


#lapply(list.df, function(x)x[x$B!=2,])


# this seems to work, but check for <NA> items in list (?!) 
dXYSm[Sm$run]


tmpFst <- h5read('/media/schimar/schimar2/bu2s/h5/Sm/fstAvgPerLocType_Sm.h5', '/runs/')
par(mfrow= c(2,4))
for (i in 1:dim(Sm)[1]){
	#tmpFst <- h5read('/media/schimar/schimar2/bu2s/h5/Sm/fstAvgPerLocType_Sm.h5')
	tmpFstA <- tmpFst[[i]]$fstAvg 
	plot(tmpFstA[,3], col= as.factor(tmpFstA[,1]))
}

# Fst for all Sm (1600) runs	
for (i in Sm$run){
	#print(dim(tmpFst[[i]]$fstAvg))
	plot(tmpFst[[i]]$fstAvg[,2], tmpFst[[i]]$fstAvg[,3], col= as.factor(tmpFst[[i]]$fstAvg[,1]))
}


#
for (i in Sm$run){
	#print(dim(tmpFst[[i]]$fstAvg))
	plot(tmpFst[[i]]$fstAvg[,2], tmpFst[[i]]$fstAvg[,3], col= as.factor(tmpFst[[i]]$fstAvg[,1]))
}





#############
# s < m  

# params for all 8784 runs
sM <- df[which(df$sd_move > df$mean_s),]













#############
# s == m  

# params for all 8784 runs
sm <- df[which(df$sd_move == df$mean_s),]



#h5run <- h5read("runs.hdf5", "/runs/")


# they all seem to be transposed, thus:
#for (i in 1:length(h5run)){
#	for (j in 1:length(h5run[[i]])){
#		h5run[[i]][[j]] <- as.data.frame(t(h5run[[i]][[j]]))
#	}
#	colnames(h5run[[i]]$Fst) <- c("totalGenerationsElapsed", "locusID", "Fst", "allele_frequencies", "S_MAX1", "S_MAX0", "chromosomeMembership", "MAP", "locType")
#	colnames(h5run[[i]]$Fix) <- c("totalGenerationsElapsed", "locusID", "MAP", "chromosomeMembership", "S_MAX1", "is_reversed_locus", "is_selected_locus")						  
#	colnames(h5run[[i]]$NeutAF) <- c("totalGenerationsElapsed", "locusID", "allele_frequencies", "AFdiff", "nsnDist", "nsnS", "nsnLD")
#	colnames(h5run[[i]]$SelAF) <- c("totalGenerationsElapsed", "locusID", "allele_frequencies", "AFdiff", "S_MAX1", "Fst", "nsnDist", "nsnS", "nsnLD")
#
#}


# no need to transpose anymore, data is read in transposed with h5py now
#for (i in 1:length(selAF)){
#	#for (j in 1:length(selAF[[i]]$selAF)){
#	selAF[[i]]$SelAF <- t(selAF[[i]]$SelAF)
#	}

##### split single df by locType   (NOTE:  check on how to subset with rhdf5 - potentially faster !!!) 
#r499split <- split(h5run$Run215499$Fst, h5run$Run215499$Fst$locType)
#r500split <- split(h5run$Run215500$Fst, h5run$Run215500$Fst$locType)


# 
#plot(r499split[[2]]$totalGenerationsElapsed, r499split[[2]]$Fst, col= 'grey20')
#points(r499split[[1]]$totalGenerationsElapsed, r499split[[1]]$Fst, col= 'grey70')



##########################################################
############   JAFSdata
#params <- read.table("../../fat_hd3tb/paramsNoOrg_v361.txt", header= T, sep= '\t')




# Run 215500 (last one on FAT_HD3TB)
# w/ numPerPatch2df.py written to r215500npp.txt

#files <- system("ls", intern= T)
## read JAFSdata
#jafsNames <- files[grep("JAFSdata", files)]
##jafs_nGen <- jafsNames[grep("[0-9]+", jafsNames]
#nGen <- as.numeric(substring(jafsNames, 9, nchar(as.character(jafsNames))-4))
#
## order/sort according to generation
#jafsNames <- jafsNames[order(nGen)]
#nGen <- sort(nGen)
#
##
#jafs215500 <- list()
#for (i in 1:length(jafsNames)){
#	curGen <- nGen[i]
#	jafs215500 <- read.csv(jafsNames[i], header= T)
#}



