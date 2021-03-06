
##################################################################################################
##################################################################################################
#library(rhdf5) 
#library(scales)

####### 
#  coupling 

# Barton 1983's coupling coefficient (theta = s/r)
# Kruuk 1999's summed coupling coefficient (phi = (L-1)s/r) 
# Le (effective nLoci)  (Le = s*/s)

#### pandapeter 
source("~/flaxmans/bu2s/bu2s_utils/couplingFuncs.r")

df <- read.table("~/flaxmans/bu2s/runs/paramsALL.txt", header= T, sep= '\t')
names(df) <- tolower(names(df))



# on ruderalis:
#source('~/schimar/bu2s/bu2s_utils/couplingFuncs.r')
#df  <- read.table("~/schimar/bu2s/runs/paramsALL.txt", header= T, sep= '\t')
#names(df) <- tolower(names(df))



#########################    simple SUBSETS 
# s > m  (1600 runs)
Sm <- df[which(df$sd_move < df$mean_s),]
Sm <- Sm[which(Sm$run != 'Run209578'),]      # something's up with dim of LDsel in this run... exclude!)

SmSub <- split(Sm, Sm$mutation_distribution)
#############
# s < m  (8784 runs)
sM <- df[which(df$sd_move > df$mean_s),]

sMsub <- split(sM, list(sM$sd_move, sM$mean_s), drop= T)
#############
# s == m  (2400 runs)

sm <- df[which(df$sd_move == df$mean_s),]

smSub <- split(sm, list(sm$sd_move, sm$mean_s), drop= T)

##
# sm$mutation_distribution[sm$run %in% smSub[[2]]$run]


############################################################################################################################
# example run (215500) 

fst <- read.table("xRuns/Run215500/FSTtimeSeries.txt.bz2")
colnames(fst) <- c("nGen", "locusID", "Fst", "allele_frequencies", "S_MAX1", "S_MAX0", "chromosomeMembership", "MAP", "locType")

# afts
afts <- read.table("xRuns/Run215500/AlleleFreqTimeSeries.txt.bz2")
colnames(afts) <- c("nGen", "locusID", "AFpatch0", "AFpatch1", "is_revesed_locus", "locType", "AF", "AFdiff")

# LDs
LDneut <- read.table('xRuns/Run215500/LDneutSitesAvg.txt.bz2')
LDsel <- read.table('xRuns/Run215500/LDselSitesAvg.txt.bz2')

# me
em <- read.table("xRuns/Run215500/EffectiveMigrationRates.txt.bz2")
em <- em[,1:18]
colnames(em)[1:6] <- c("nGen", "eme0", "eme1", "nVariableLoci", "nRes", "nImm")



#r499 <- ccStats(df, fst, afts, 'Run215500')
r499 <- ccStats.1(df, fst, afts, LDsel, LDneut, em, 'Run215500', maf= 25e-4)

plotStatic.1(r499, 'Run215500', df)


############################################################################################################################


#plot(afts$nGen, afts$AFdiff, ylab= 'afDiff', xlab= 'nGen', pch= '.', cex= 2.5, col= alpha(afts$locType +1, 0.7))

###########################################################################################################################
###########################  hdf5  loop

# test subset (from sM)
#rs <- paste('Run215', seq(499, 550, 1), sep= '')
#paramSub <- df[df$run %in% rs,]
paths <-  '/media/schimar/schimar2/bu2s/h5/'


paramSub <- sMsub[[2]][1:10,]
wrapH5static(paramSub, 'sM', path= paths)

wrapH5static(sM[1:100,], 'sM')


wrapH5static(paramSub, 'Sm')

wrapH5static(smSub[[1]], 'sm')

smSym <- sm[801:2400,][which(sm[801:2400,]$end_period_allopatry == -1),]
wrapH5static(smSym, 'sm', path= paths)



#Sm[which(Sm$run == 'Run203006'),which(dfVar != 0)]



# Sm 
# phis_Sm <- xtractPhis(Sm, 'Sm', maf= 0.025)
phis_Sm <- xtractPhis(Sm, setname= 'Sm', folder= 'Sm3', maf= 0.025)

# sm 
phis_sm <- xtractPhis(sm, 'sm', 'sm', maf= 0.025)		#path= '/media/schimar/FLAXMAN/h5/', 



# sM
#phis_sM <- xtractPhis(sM, setname= 'sM', folder= 'sM2', maf= 0.025)

# running on ruderalis as "phisM"
# resume with: screen -r phisM

phis_sM1 <- xtractPhis(sMsub[[1]], setname= 'sM', folder= 'sM2', maf= 0.025)
#### phis_sM2 <- xtractPhis(sMsub[[2]], setname= 'sM', folder= 'sM2', maf= 0.025)
phis_sM3 <- xtractPhis(sMsub[[3]], setname= 'sM', folder= 'sM2', maf= 0.025)
phis_sM4 <- xtractPhis(sMsub[[4]], setname= 'sM', folder= 'sM2', maf= 0.025) 		# path= '/media/schimar/schimar2/bu2s/h5/' 
phis_sM5 <- xtractPhis(sMsub[[5]], setname= 'sM', folder= 'sM2', maf= 0.025) 		# path= '/media/schimar/schimar2/bu2s/h5/' 
phis_sM6 <- xtractPhis(sMsub[[6]], setname= 'sM', folder= 'sM2', maf= 0.025) 		# path= '/media/schimar/schimar2/bu2s/h5/' 
phis_sM7 <- xtractPhis(sMsub[[7]], setname= 'sM', folder= 'sM2', maf= 0.025) 		# path= '/media/schimar/schimar2/bu2s/h5/' 
phis_sM8 <- xtractPhis(sMsub[[8]], setname= 'sM', folder= 'sM2', maf= 0.025) 		# path= '/media/schimar/schimar2/bu2s/h5/' 
phis_sM9 <- xtractPhis(sMsub[[9]], setname= 'sM', folder= 'sM2', maf= 0.025) 		# path= '/media/schimar/schimar2/bu2s/h5/' 



# xtractLe 

Le_Sm <- xtractLe(Sm, setname= 'Sm', folder= 'Sm', maf= 0.025, path= '/media/schimar/dapperdata/bu2s/h5/')
Le_sm <- xtractLe(sm, setname= 'sm', folder= 'sm', maf= 0.025)

Le_sM1 <- xtractLe(sMsub[[1]], setname= 'sM', folder= 'sM2', maf= 0.025)
Le_sM3 <- xtractLe(sMsub[[3]], setname= 'sM', folder= 'sM2', maf= 0.025)
Le_sM4 <- xtractLe(sMsub[[4]], setname= 'sM', folder= 'sM2', maf= 0.025)
Le_sM5 <- xtractLe(sMsub[[5]], setname= 'sM', folder= 'sM2', maf= 0.025)
Le_sM6 <- xtractLe(sMsub[[6]], setname= 'sM', folder= 'sM2', maf= 0.025)
Le_sM7 <- xtractLe(sMsub[[7]], setname= 'sM', folder= 'sM2', maf= 0.025)
Le_sM8 <- xtractLe(sMsub[[8]], setname= 'sM', folder= 'sM2', maf= 0.025)
Le_sM9 <- xtractLe(sMsub[[9]], setname= 'sM', folder= 'sM2', maf= 0.025)






# Sm
LeSm <- lapply(lapply(Le_Sm, '[[', 1), '[[', 2)
#sStSm <- lapply(lapply(Le_Sm, '[[', 1), '[[', 1)
minSm <- unlist(lapply(LeSm, min, na.rm= T))

# sm

Lesm <- lapply(lapply(Le_sm, '[[', 1), '[[', 2)

minsm <- unlist(lapply(Lesm, min, na.rm= T))


# sM

LesM1 <- lapply(lapply(Le_sM1, '[[', 1), '[[', 2)
minsM1 <- unlist(lapply(LesM1, min, na.rm= T))

LesM3 <- lapply(lapply(Le_sM3, '[[', 1), '[[', 2)
minsM3 <- unlist(lapply(LesM3, min, na.rm= T))

LesM4 <- lapply(lapply(Le_sM4, '[[', 1), '[[', 2)
minsM4 <- unlist(lapply(LesM4, min, na.rm= T))

LesM5 <- lapply(lapply(Le_sM5, '[[', 1), '[[', 2)
minsM5 <- unlist(lapply(LesM5, min, na.rm= T))

LesM6 <- lapply(lapply(Le_sM6, '[[', 1), '[[', 2)
minsM6 <- unlist(lapply(LesM6, min, na.rm= T))

LesM7 <- lapply(lapply(Le_sM7, '[[', 1), '[[', 2)
minsM7 <- unlist(lapply(LesM7, min, na.rm= T))

LesM8 <- lapply(lapply(Le_sM8, '[[', 1), '[[', 2)
minsM8 <- unlist(lapply(LesM8, min, na.rm= T))

LesM9 <- lapply(lapply(Le_sM9, '[[', 1), '[[', 2)
minsM9 <- unlist(lapply(LesM9, min, na.rm= T))

minsM <- c(minsM1, minsM3, minsM4, minsM5, minsM6, minsM7, minsM8, minsM9)

# 
par(mfrow= c(1,3))

plot(minSm,  type= 'l', main= 'Sm')
plot(minsm, type= 'l', main= 'sm')
plot(minsM, type= 'l', main= 'sM')

plot(minSm, type= 'l', main= 'Sm')
plot(minsm, type= 'l', main= 'sm')
plot(minsM[which(sM$mutation_distribution == 0)], type= 'l', main= 'sM')


LesmPallo <- Lesm[which(sm$end_period_allopatry != -1)]
LDsmPallo <- LDsm[which(sm$end_period_allopatry != -1)]


plot(LesmPallo[[751]], LDsmPallo[[751]]$LDneut[,4], xlim= c(-200, 10), ylim= c(0,1), type= 'n')
for (i in 751:length(LesmPallo)) {
	cLe <- LesmPallo[[i]]
	cLDn <- LDsmPallo[[i]]$LDneut[,4]
	cLDs <- LDsmPallo[[i]]$LDsel[,4]
	points(cLe, cLDs, type= 'l', col= 'red')
	points(cLe, cLDn, type= 'l', col= 'blue')
}



LeSmPallo <- LeSm[which(Sm$end_period_allopatry != -1)]
LDSmPallo <- LDSm[which(Sm$end_period_allopatry != -1)]


plot(LeSmPallo[[751]], LDSmPallo[[751]]$LDneut[,4], xlim= c(-200, 10), ylim= c(0,1), type= 'n')
for (i in 751:length(LeSmPallo)) {
	cLe <- LeSmPallo[[i]]
	cLDn <- LDSmPallo[[i]]$LDneut[,4]
	cLDs <- LDSmPallo[[i]]$LDsel[,4]
	points(cLe, cLDs, type= 'l', col= 'red')
	points(cLe, cLDn, type= 'l', col= 'blue')
}


LesMPallo <- LesM[which(sM$end_period_allopatry != -1)]
LDsMPallo <- LDsM[which(sM$end_period_allopatry != -1)]


plot(LesMPallo[[751]], LDsMPallo[[751]]$LDneut[,4], xlim= c(-200, 10), ylim= c(0,1), type= 'n')
for (i in 751:length(LesMPallo)) {
	cLe <- LesmPallo[[i]]
	cLDn <- LDsmPallo[[i]]$LDneut[,4]
	cLDs <- LDsmPallo[[i]]$LDsel[,4]
	points(cLe, cLDs, type= 'l', col= 'red')
	points(cLe, cLDn, type= 'l', col= 'blue')
}


#################################
# extract LD for each set of runs

LDSm <- xtractLD(Sm, setname= 'Sm', folder= 'Sm3')
LDsm <- xtractLD(sm, setname= 'sm', folder= 'sm')

LD_sM1 <- xtractLD(sMsub[[1]], setname= 'sM', folder= 'sM2')
LD_sM3 <- xtractLD(sMsub[[3]], setname= 'sM', folder= 'sM2')
LD_sM4 <- xtractLD(sMsub[[4]], setname= 'sM', folder= 'sM2')
LD_sM5 <- xtractLD(sMsub[[5]], setname= 'sM', folder= 'sM2')
LD_sM6 <- xtractLD(sMsub[[6]], setname= 'sM', folder= 'sM2')
LD_sM7 <- xtractLD(sMsub[[7]], setname= 'sM', folder= 'sM2')
LD_sM8 <- xtractLD(sMsub[[8]], setname= 'sM', folder= 'sM2')
LD_sM9 <- xtractLD(sMsub[[9]], setname= 'sM', folder= 'sM2')



plot(Lesm[[801]], LDsm[[801]]$LDsel[,4], type= 'n', xlim= c(-462, 10))
for (i in 801:1599) {
	points(Lesm[[i]], LDsm[[i]]$LDsel[,4], type= 'l')
}

# load the .RData for PHIs 

load('~/flaxmans/bu2s/runs/PHIs/sm/.RData')
load('~/flaxmans/bu2s/runs/PHIs/S_m/.RData')
load('~/flaxmans/bu2s/runs/PHIs/sM/.RData')

# on ruderalis
#load("~/schimar/bu2s/runs/Le/.RData")
#load("~/schimar/bu2s/runs/LD/.RData")
#load("~/schimar/bu2s/runs/PHIs/sm/.RData")
#load("~/schimar/bu2s/runs/PHIs/sM/.RData")
#load("~/schimar/bu2s/runs/PHIs/S_m/.RData")

#load("~/schimar/bu2s/runs/PHIs/sM/phis_sM1.RData")




# plot phis & afDiffS
# sm
plot(log10(phis_sm[[1]][[2]]), type= 'n', ylab= expression(phi)) #, ylim= c(-2, 1e+11)
for (i in 1:length(phis_sm[801:1600])) {
	phiO <- phis_sm[[i]][[1]]
	kphi <- phis_sm[[i]][[2]]
	afDiffS <- phis_sm[[i]][[3]]
	afDiffN <- phis_sm[[i]][[4]]
	points(log10(phis_sm[[i]][[2]]), col= 'black', type= 'l')   # kphisMax
	points(log10(phis_sm[[i]][[1]]), col= 'grey70', type= 'l')  # phiObs
}


plot(log10(phis_Sm[[1]]$kphisMax), type= 'n', ylab= expression(phi)) #, ylim= c(-2, 1e+11)
for (i in 1:length(phis_Sm)) {
	points(log10(phis_Sm[[i]]$kphisMax), col= 'black', type= 'l')
	points(log10(phis_Sm[[i]]$phiObs), col= 'grey70', type= 'l')
}


#
plot(log10(phis_sM[[1]]$kphisMax), type= 'n', ylab= expression(phi)) #, ylim= c(-2, 1e+11)
for (i in 1:length(phis_sM)) {
	points(log10(phis_sM[[i]]$kphisMax), col= 'black', type= 'l')
	points(log10(phis_sM[[i]]$phiObs), col= 'grey70', type= 'l')
}


plotPhis <- function(phis, ...) {
	plot(log10(phis[[1]]$kphisMax), phis[[1]]$pHatsMax, type= 'l', ylim= c(-0.45,  1), xlim= c(-2, max(log10(phis[[1]]$kphisMax))), xlab= expression(paste('log'[10], ' ', phi)), ylab= expression(paste('p'[i0]~'- p'[i1])), col= 'black')
	legend('bottomright', legend= c(expression(paste(phi[Kruuk], ' ~ ', 'peq '[sMax])), expression(paste(phi, ' ~ ', bar(p), ' all s')), expression(paste(phi, ' ~ avg p S')), expression(paste(phi, ' ~ avg p N'))), fill= c('black', 'grey70', 'red', 'blue'), cex= 0.8)
	#
	for (i in 1:length(phis)) {
		phiGen <- phis[[i]]
		cWallS <- phiGen$cWallS
		phiOncW <- mapply(rep, phiGen$phiObs, times= unlist(lapply(cWallS, length)))
		# 
		points(log10(phis[[1]]$kphisMax), phis[[1]]$pHatsMax, type= 'l')
		points(log10(unlist(phiOncW)), unlist(cWallS), pch= '.', col= 'grey70') #type= 'l')
		points(log10(unlist(phiGen$phiObs)), unlist(lapply(lapply(phiGen$afDiffS, abs), mean)), pch= '.', cex= 1.4, col= 'red')
		points(log10(unlist(phiGen$phiObs)), unlist(lapply(lapply(phiGen$afDiffN, abs), mean)), pch= '.', cex= 1.4, col= 'blue')
			}
}

# this is ok, but I should rather pull out the overall stats (mean & sd for phis across runs) 



dev.new()

par(mfrow= c(4,4))
plotPhis(phis_Sm[1:800])
plotPhis(phis_Sm[801:1600])
plotPhis(phis_sm[1:800])
plotPhis(phis_sm[801:1600])
plotPhis(phis_sm[1601:2400])

plotPhis(phis_sM1)

plotPhis(phis_sM3)
plotPhis(phis_sM4)
plotPhis(phis_sM5)
plotPhis(phis_sM6)
plotPhis(phis_sM7)
plotPhis(phis_sM8)
plotPhis(phis_sM9)




#for (i in length(phis_Sm)) {
	
############################ 
# plot sets of phis 

######################### manual run 
#load('sM_cWs_pcWallS.RData')


phis <- phis_sm[801:1600]


kps <- unlist(lapply(phis, '[[', 2))
pHs <- unlist(lapply(phis, '[[', 3))

phO <- unlist(lapply(phis, '[[', 1))
pN <- unlist(lapply(phis, '[[', 5))
pS <- unlist(lapply(phis, '[[', 4))

cWs <- unlist(lapply(phis, '[[', 6))
#

# get the 'phiOncW' for selected and neutral afDiffs
phiOncW <- list()
for (i in 1:length(phis)) {
	#k <- i + 800
	phiOncW[[i]] <- mapply(rep, phis[[i]]$phiObs, times= unlist(lapply(phis[[i]]$cWallS, length)))
}


# flatten 
pcWallS <- unlist(lapply(phiOncW, unlist))

gc()
rm(phis)
rm(phis_sm)

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
#allS <- as.data.frame(cbind(log10(pcWallS), cWs, rep('D', length(cWs))))
#colnames(allS) <- c('phi', 'p', 'group')

cWss <- cWs[c(TRUE, FALSE, FALSE)]
lpcWallSs <- log10(x$pcWallS[c(TRUE, FALSE, FALSE)])

##
matplot(pcWallS, cWs, type= 'p', pch= '.', cex= 1.0, col= 'grey80', ylim= c(0,1), xlim= c(-2, 4.5), xlab= expression(paste('log'[10], ' ', phi)), ylab= expression(paste('p'[i0]~'- p'[i1])), cex.lab= 1.4, cex.axis= 1.4)
matpoints(Sphi$phi, Sphi$p, pch= '.', col= 'firebrick1')
matpoints(kps_pHs$phi, kps_pHs$p, pch= '.', col= 'black', cex= 1.0)
matpoints(Nphi$phi, Nphi$p, pch= '.', col= 'deepskyblue1')
legend('topleft', legend= c(expression(paste(phi[Kruuk], ' ~ ', 'peq '[sMax])), expression(paste(phi, ' ~ ', bar(p), ' all s')), expression(paste(phi, ' ~ avg p S')), expression(paste(phi, ' ~ avg p N'))), fill= c('black', 'grey70', 'firebrick1', 'deepskyblue1'), cex= 0.8)

##################################################################

matplot(1:10, type= 'n', xlab= 'generations * 1e-3', ylab= expression('log'[10], ' (metric)'), cex.lab= 1.6, cex.axis= 1.4, xlim= c(0, 600), ylim= c(-0.6,0.6))

# create matrices with runs as rows and generations as columns
#kphis <- do.call(rbind, lapply(phis_Sm, '[[', 2))
#pHatsMax <- do.call(rbind, lapply(phis_Sm, '[[', 3))
#afDn <- do.call(rbind, lapply(phis_Sm, '[[', 5))
#afDs <- do.call(rbind, lapply(phis_Sm, '[[', 4))
#phiObs <- do.call(rbind, lapply(phis_Sm, '[[', 1))
#cWallS  <- do.call(rbind, lapply(phis_Sm, '[[', 6))
#
#


# all PHIs  (with calcPHIs.3) 

s500 <- ccStats.3('Run215500', df, cc500, maf= 0.025)


# phis per Gen

plot(log10(s500$kphisMax), type= 'l', col= 'green')
points(log10(s500$phiBarsMax), type= 'l', col= 'darkgreen')

points(log10(s500$phiBarMeanS), type= 'l', col= 'grey70')
points(log10(kphiMeanS), type= 'l')
legend('topleft', legend= c('kphi sMax', 'phiBar sMax', 'phiBar MeanS', 'kphi MeanS'), fill= c('green', 'darkgreen', 'grey70', 'black'), cex= 1.4)


# phis ~ AFs & peq

ccObj <- s500

plot(log10(ccObj$kphisMax), ccObj$pHatsMax, type= 'l', ylim= c(0,  1), xlim= c(-2, 4), xlab= expression(paste('log'[10], ' ', phi)), ylab= expression(paste('p'[i0]~'- p'[i1])), col= 'black')
	#abline(h= 0, lty= 3)
	#}
	points(log10(unlist(phiOncW)), unlist(cWallS), pch= '.', col= 'grey70') #type= 'l')
	#points(log10(unlist(ccObj$phiObs)), unlist(lapply(lapply(ccObj$afDiffS, abs), mean)), pch= '.', cex= 1.4, col= 'firebrick1')
	#points(log10(unlist(ccObj$phiObs)), unlist(lapply(lapply(ccObj$afDiffN, abs), mean)), pch= '.', cex= 1.4, col= 'deepskyblue1')
	points(log10(kphiMeanS), unlist(lapply(lapply(ccObj$afDiffS, abs), mean)), pch= '.', cex= 1.4, col= 'green')
	points(log10(kphiMeanS), unlist(lapply(lapply(ccObj$afDiffN, abs), mean)), pch= '.', cex= 1.4, col= 'purple')


# better write this shit to a file
write.table(allPhi, file= 'allPhi_Sm.txt', quote= F, col.names= T, row.names= F, sep= '\t')


# 
allPhi <- read.table('allPhi_Sm.txt', header= T, sep= '\t')
#####  finally some plotting 
library(ggplot2)

#
# call func 
#allPhi <- phi2plot(phis_Sm[801:900])
#
## 
#gg_allPhi <- ggplot(allPhi, aes(x= phi, y= p, colour= group)) + theme_bw() + coord_cartesian(xlim= c(-1.8, 4), ylim= c(0,1))
#
#gg_allPhi + stat_summary(aes(x= phi, y= p, colour= group), geom= 'ribbon', fun.ymin= 'min', fun.ymax= 'max', alpha= 0.7)


allPhiSm <- phi2plot(phis_Sm[801:1599])


allPhism <- phi2plot(phis_sm[801:1600])

allPhisM1 <- phi2plot(phis_sM1)


phis <- allPhisM1
#
kphi <- phis[which(phis$group == 'C'),]
pS <- phis[which(phis$group == 'A'),]
pN <- phis[which(phis$group == 'B'),]
allS <- phis[which(phis$group == 'D'),]




# R base plot 
#phis <- allPhism
#
##
matplot(allS$phi, allS$p, type= 'p', pch= '.', cex= 1.0, col= 'grey80', ylim= c(0,1), xlim= c(-2, 3.5), xlab= expression(paste('log'[10], ' ', phi)), ylab= expression(paste('p'[i0]~'- p'[i1])))
matpoints(pS$phi, pS$p, pch= '.', col= 'firebrick1')
matpoints(kphi$phi, kphi$p, pch= '.', col= 'black', cex= 1.0)
matpoints(pN$phi, pN$p, pch= '.', col= 'deepskyblue1')
legend('topleft', legend= c(expression(paste(phi[Kruuk], ' ~ ', 'peq '[sMax])), expression(paste(phi, ' ~ ', bar(p), ' all s')), expression(paste(phi, ' ~ avg p S')), expression(paste(phi, ' ~ avg p N'))), fill= c('black', 'grey70', 'firebrick1', 'deepskyblue1'), cex= 0.8)

#legend('topleft', legend= c(expression(paste(phi[Kruuk], ' ~ ', 'peq '[sMax])), expression(paste(phi, ' ~ ', bar(p), ' all s')), expression(paste(phi, ' ~ avg p S')), expression(paste(phi, ' ~ avg p N'))), fill= c('black', 'grey70', 'firebrick1', 'deepskyblue1'), cex= 0.8)

######
phis <- allPhism
#
kphi <- phis$kphi
pS <- phis$pS
pN <- phis$pN
allS <- phis$allS


#
plot(allS$phi, allS$p, type= 'p', pch= '.', cex= 1.0, col= 'grey80', ylim= c(0,1), xlim= c(-2, 3.5), xlab= expression(paste('log'[10], ' ', phi)), ylab= expression(paste('p'[i0]~'- p'[i1])))
points(pS$phi, pS$p, pch= '.', col= 'firebrick1')
points(kphi$phi, kphi$p, pch= '.', col= 'black')
points(pN$phi, pN$p, pch= '.', col= 'deepskyblue1')

legend('topleft', legend= c(expression(paste(phi[Kruuk], ' ~ ', 'peq '[sMax])), expression(paste(phi, ' ~ ', bar(p), ' all s')), expression(paste(phi, ' ~ avg p S')), expression(paste(phi, ' ~ avg p N'))), fill= c('black', 'grey70', 'firebrick1', 'deepskyblue1'), cex= 1.4)


#qplot(kphi$phi, kphi$p, geom='smooth', span =0.5)
############### test ggplot 

	#coord_cartesian(xlim= c(-2, 5), ylim= c(0,1)) +
	#geom_line(data= afDiffS)
	#theme_bw()

#ggkpH <- ggplot(kps_pHs, aes(x= lkps, y= pHs))  #+ coord_cartesian(xlim= c(-2, 5), ylim= c(0,1)) + theme_bw()

#ggkpH+ 
	#stat_summary(geom= 'ribbon', fun.ymin= 'min', fun.ymax= 'max', alpha= 0.7) + 
#		geom_point(aes(x= phO, y= pS, colour= 'red')) # phiObs ~ afdiffS 
	
ggAFDphiObs <- ggplot(phiOAFd, aes(x= log10(phiObs), y= pS))  + theme_bw() + coord_cartesian(xlim= c(-1.8, 4), ylim= c(0,1)) 	

ggAFDphiObs + stat_summary(geom= 'ribbon', aes(colour= 'red'), fun.ymin= 'min', fun.ymax= 'max', alpha= 0.7) +
	stat_summary(aes(x= log10(phO), y= pN, colour= 'blue'), geom= 'ribbon', fun.ymin= 'min', fun.ymax= 'max', alpha= 0.7) + 
	stat_summary(aes(x= log10(kps), y= pHs, colour= 'grey70'), geom= 'ribbon',  fun.ymin= 'min', fun.ymax= 'max', alpha= 0.7) 


	# geom_point(aes(x= log10(pcWallS), y= cWs, colour= 'grey20'))
	




#stat_summary(geom= 'ribbon', fun.ymin= 'min', fun.ymax= 'max', alpha= 0.5)


#ggplot(dat, aes(x=pos,y=value, colour=type)) +
#	stat_summary(geom="ribbon", fun.ymin="min", fun.ymax="max", aes(fill=type), alpha=0.3) +
#	theme_bw()

ggkpH +
	#geom_ribbon(aes(ymin= kps_pHs$pHs -0.02, ymax= kps_pHs$pHs +0.02), fill= 'grey70') +
	geom_line(aes(x= lkps, y= pHs))





	geom_smooth(method="gam", span=0.75, se=TRUE, alpha=0.3) +
  	theme_bw()




####
# Sm
#plot(log10(phis_Sm[[1]]$kphisMax), type= 'n', ylab= expression(phi)) #, ylim= c(-2, 1e+11)
plot(log10(cGenPHIs$kphisMax), cGenPHIs$pHatsMax, type= 'n') #, ylim= c(-0.45,  1), xlim= c(-2, max(log10(cGenPHIs$kphisMax))), xlab= expression(paste('log'[10], ' ', phi)), ylab= expression(paste('p'[i0]~'- p'[i1])), col= 'black')


for (i in 1:length(phis_Sm)) {
	cGenPHIs <- phis_Sm[[i]]
	cWallS <- lapply(cGenPHIs$clineWallS, unlist)
	phiOncW <- mapply(rep, cGenPHIs$phiObs, times= unlist(lapply(cWallS, length)))

	plot(log10(cGenPHIs$kphisMax), cGenPHIs$pHatsMax, type= 'l', ylim= c(-0.45,  1), xlim= c(-2, max(log10(cGenPHIs$kphisMax))), xlab= expression(paste('log'[10], ' ', phi)), ylab= expression(paste('p'[i0]~'- p'[i1])), col= 'black')
	#abline(h= 0, lty= 3)
	#}
	points(log10(unlist(phiOncW)), unlist(cWallS), pch= '.', col= 'grey70') #type= 'l')
	points(log10(unlist(cGenPHIs$phiObs)), unlist(lapply(lapply(cGenPHIs$afDiffS, abs), mean)), pch= '.', cex= 1.4, col= 'red')
	points(log10(unlist(cGenPHIs$phiObs)), unlist(lapply(lapply(cGenPHIs$afDiffN, abs), mean)), pch= '.', cex= 1.4, col= 'blue')
	}





###   plot phis 
plot(log10(phis_Sm[[801]]$phiObs), type= 'n')

for (i in 801:1599) {
	points(log10(phis_Sm[[i]]$phiObs), type= 'l', col= 'grey70')
	#points(log10(phis_Sm[[i]]$kphisMax), type= 'l', col= 'grey70')
}


##
plot(log10(phis_sm[[1]]$phiObs), type= 'n')

for (i in 1601:2400) {
	points(log10(phis_sm[[i]]$phiObs), type= 'l') #, col= 'grey70')
	#points(log10(phis_Sm[[i]]$kphisMax), type= 'l', col= 'grey70')
}


plot(log10(phis_sM5[[1]]$phiObs), type= 'n')

for (i in 1:800) {
	points(log10(phis_sM6[[i]]$phiObs), type= 'l', col= 'grey70')
	#points(log10(phis_Sm[[i]]$kphisMax), type= 'l', col= 'grey70')
}




#paste(pathTmp, set, "/afts_", set, 'T.h5', sep= '')

###########################################################################################################################
# ccStats2 

#c202971 <- readIn('Run202971', df, 'Sm', path= '/media/schimar/schimar2/bu2s/h5/')

# c500 <- readIn('Run215499', df, 'sM', path= '/media/schimar/schimar2/bu2s/h5/')

cc500 <- readCCobj('Run215500', 'sM', path= '/media/schimar/schimar2/bu2s/h5/')
s500 <- ccStats.2('Run215500', df, cc500, maf= 0.025)


# on ruderalis

cc500 <- readCCobjRude('Run215500', smParam= 'sM', folder= 'sM2', path= '/home/schimar/FLAXMAN/h5/')
s500 <- ccStats.2('Run215500', df, cc500, maf= 0.025)

#####  

# for poster: Run206171 (sm[801,])
cc171 <- readCCobj('Run206171', 'sm', path= '/media/schimar/dapperdata/bu2s/h5/') 	#paths)
s171 <- ccStats.2('Run206171', df, cc171, maf= 0.025)

cc210 <- readCCobj('Run210070', 'sm', path= '/media/schimar/dapperdata/bu2s/h5/')
s210 <- ccStats.2('Run210070', df, cc210, maf= 0.025)



# sMsub[[1]][1,]
cc571 <- readCCobj('Run212571', 'sM', path= '/media/schimar/dapperdata/bu2s/h5/')
s571 <- ccStats.2('Run212571', df, cc571, maf= 0.025)




#slim500 <- ccStats.2slim('Run215500', df, cc500, maf= 0.025)


# Run206976
#ccT <- readCCobj('Run206976', 'sm', path= '/media/schimar/schimar2/bu2s/h5/')
#rcT <- ccStats(df, ccT$fst, ccT$afts, ccT$LDsel, ccT$LDneut, ccT$effMig, 'Run206976', maf= 0.0)
#
## Run206976 (run with ts_sampling_freq = 194) 
##ccT578 <- readCCobj('Run209578', 'Sm', path= '/media/schimar/schimar2/bu2s/h5/')
##rcT578 <- ccStats.2('Run209578', df, ccT578, maf= 0.025)
##
##ccT577 <- readCCobj('Run209577', 'Sm', path= '/media/schimar/schimar2/bu2s/h5/')
##rcT577 <- ccStats.2('Run209577', df, ccT577, maf= 0.025)
##
##cc203124 <- readCCobj('Run203124', 'sM', path= '/media/schimar/schimar2/bu2s/h5/')
##rc203124 <- ccStats.2('Run203124', df, cc203124, maf= 0.025)
##
### in sm
##
cc207575 <- readCCobj('Run207575', 'sm', path= '/media/schimar/schimar2/bu2s/h5/')
##
rc207575 <- ccStats.2('Run207575', df, cc207575, maf= 0.025)


#lapply(list.df, function(x)x[x$B!=2,])

plotStatic.2(s500, 'Run215500', df)


fstSpl <- split(fst, fst$nGen)
aftsSpl <- split(afts, afts$nGen)

####
# plot MAP ~ Fst    (& dXY, LD, etc.)    NOTE:  not by chromosome !!
#####################################################

##### spatial autocorrelation 

###ct <- cut(frst$MAP, quantile(frst$MAP, probs = seq(0, 1, by = 1/5)), right= F, labels= 1:5, include.lowest= T)

##################################################### calc Moran's I for distance bins  (5)

# calc Moran's I per chromosome (no distance bins within) 
mIKgen <- calcMorIripK(fstSpl)

# now with distance bins (k = 5)
mIbin <- calcMorIbin(fstSpl)

#IKSm <- xtractIK(Sm, setname= 'Sm', folder= 'Sm3')
#IKsm <- xtractIK(sm, setname= 'sm', folder= 'sm')

### 
IK_sM1 <- xtractIK(sMsub[[1]][1:400,], setname= 'sM', folder= 'sM2', path= '/media/schimar/FLAXMAN/h5/')
IK_sM1.2 <- xtractIK(sMsub[[1]][1:400,], setname= 'sM', folder= 'sM2', path= '/media/schimar/FLAXMAN/h5/')



IK_sM3 <- xtractIK(sMsub[[3]], setname= 'sM', folder= 'sM2')
IK_sM4 <- xtractIK(sMsub[[4]], setname= 'sM', folder= 'sM2')
IK_sM5 <- xtractIK(sMsub[[5]], setname= 'sM', folder= 'sM2')

#gwc <- wrapGWCtime('/media/schimar/FLAXMAN/h5/', sm[801:1600,], setname= 'sm', folder= 'sm')


IK_sM6 <- xtractIK(sMsub[[6]], setname= 'sM', folder= 'sM2')
IK_sM7 <- xtractIK(sMsub[[7]], setname= 'sM', folder= 'sM2')
IK_sM8 <- xtractIK(sMsub[[8]], setname= 'sM', folder= 'sM2')
IK_sM9 <- xtractIK(sMsub[[9]], setname= 'sM', folder= 'sM2')


# get the absolute values and calc mean (for each chromosome) 
# then plot the 4 chromosomes' means across generations
#  do we see changes in values before/at/after GWC?) 

#### various wrapper functions  

wrapH5dynaGen(sm, 'sm', time= 80)
wrapH5dynaGen(sm[which(sm$end_period_allopatry != -1),], 'sm')

wrapH5dynaBin(sm[801:1600,], 'sm', time= 1)

wrapH5FstMAPmIgen(sm, 'sm', time= 1)


# plot Ripley's K    (isotropic correction, w/ 99 simulations of CSR (default "envelope()  ))
	
plotFstMAPripKmorIbin(fstSpl, mIbin, mIKgen, static= F, wait= 0.8, time= 600)

#plot(envelope(xpp, Kest))

# from ?Kest 

# 	   The estimates of K(r) are of the form

#           Kest(r) = (a/(n * (n-1))) * sum[i,j] I(d[i,j] <= r) e[i,j])   
#      
#      where a is the area of the window, n is the number of data points,
#      and the sum is taken over all ordered pairs of points i and j in
#      ‘X’.  Here d[i,j] is the distance between the two points, and
#      I(d[i,j] <= r) is the indicator that equals 1 if the distance is
#      less than or equal to r.  The term e[i,j] is the edge correction
#      weight (which depends on the choice of edge correction listed
#      above).

########  evol2017

plotFstMAPmorIgen(fstSpl, mIKgen, static= F, time= 600)


plotAFdiffs(r499$afDiffS, r499$afDiffN, wait= 0.1, time= 610)




######
blue <- '#00BFFF'
red <- '#FF3030'


plotAFdiffs <- function(xS, xN, time= 1, wait= 0.3){
	close.screen(all.screens= T)
	par(bg = "white")
	split.screen(c(2,1))
	for (i in time:1250){
		screen(1)
		hist(xS[[i]], xlim= c(-1,1), ylim= c(0, 1000), xlab= 'allele freq diffs.', ylab= '# of sites', main= i, col= red, cex.lab= 1.8, cex.axis= 1.8)
		screen(2)
		hist(xN[[i]], xlim= c(-1,1), ylim= c(0, 1000), xlab= 'allele freq diffs.', ylab= '# of sites', main= '', col= blue,cex.lab= 1.8, cex.axis= 1.8)
		Sys.sleep(wait)
	}
}


##############################################
# plot the Moran's I values 
plotFstMAPmorIbin(fstSpl, mIbin, static= F, time= 620)


plotFstMAPmorIgen(fstSpl, mIKgen, static= F, time= 1)
#plotMorIbin(mIgen, time= 630)

plotMorI(mIgen, time= 600)


plotDyna(s500, mIgen, static= F, time= 1)



##
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
stop("vectors must be same length")
arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

# barx <- barplot(y.means, names.arg=1:5,ylim=c(0,1.5), col="blue", axis.lty=1, xlab="Replicates", ylab="Value (arbitrary units)")
# error.bar(barx,y.means, 1.96*y.sd/10)
##




##############################################

# not done yet; how to arrange the plot ??
# NOTE: see plotDyna()
plotFstMAPmorI <- function(fstspl, morI, time= 1,...) {
	# function to plot both the MAP vs Fst and Moran's I and Moran's I per generation and chromosome
	# with INPUT: fst data split by nGen, the nested list of Moran's I values (per gen and chrom) and the generation time to start the loop (with nGen/tsFreq).
	
	par(mfrow= c(2,4))
	for(i in time:length(fstspl)) {
		chrom <- fstspl[[i]]$chromosomeMembership
		for (j in 1:length(unique(chrom))) {
			cChrom <- fstspl[[i]][which(chrom == unique(chrom)[j]),]
			plot(cChrom$MAP, cChrom$Fst, col= as.factor(cChrom$locType), xlim= c(0, 25), ylim= c(0,1), main= paste('chrom ', j-1, ' gen = ', i), ylab= 'Fst', xlab= 'map', pch= 20)
			lines(cChrom$MAP, cChrom$Fst, col= 'grey70')
		}
		cMorI <- morI
		barx <- barplot(cMorI$obs, ylim= c(-1, 1), pch= 20, main= paste('Morans I; gen = ', i), names.arg= c(0,1,2,3), xlab= 'chromosome', ylab= "Moran's I")
		#lines(cMorIall$obs, col= 'grey70')
		error.bar(barx, cMorIall$obs, cMorIall$sd)
	}
}


##########################

# get obs and sd from IKsm   
load("~/schimar/bu2s/runs/mIripK/Sm/.RData")


#Ism <- lapply(IKsm[801:1600], '[[', 1)

Ism <- lapply(IKSm[801:1599], '[[', 1)

Ism <- lapply(IK_sM1, '[[', 1)



Iobs <- list()
Isd <- list()
for (i in 1:length(Ism)) {
	Iobs[[i]] <- as.data.frame(do.call(rbind, lapply(Ism[[i]], '[[', 1)))
	Isd[[i]] <- as.data.frame(do.call(rbind, lapply(Ism[[i]], '[[', 3)))
}


## 
plot(Iobs[[1]][,1], type= 'l', ylim= c(-0.5, 0.5))
points(Iobs[[1]][,2], type= 'l', col= 'deepskyblue2')
points(Iobs[[1]][,3], type= 'l', col= 'firebrick3')
points(Iobs[[1]][,4], type= 'l', col= 'chartreuse1')

#############################################

chrom1sm <- lapply(Iobs, '[[', 1)

cols <- rainbow(500)
plot(chrom1sm[[1]], type= 'n', xlim= c(0, 675), ylim= c(-0.6, 0.7), xlab= 'generations / 1e3', ylab= "Moran's I", cex.axis= 1.6)
for (i in 1:length(chrom1sm)) {
	points(chrom1sm[[i]], type= 'l', col= cols[i])
}
abline(h= 0, lty= 2)


###
wrapGWCtime <- function(path, data, setname, folder, ...) {
	# function to read individual runs (from vector of runs), calculate CC and plotStatic
	gwcTimes <- list()
	for (i in 1:dim(data)[1]){
		run <- data$run[i]
		path5 <- paste('/runs/', run, sep= '')
		#
		ccObjTmp <- readCCobjRude(run, setname, folder, path)
		ccTmp <- ccStats.2(run, data, ccObjTmp, maf= 0.05)
		#
		gwcTimes[[i]] <- ccTmp$gwcTimeMeanS$gwcTime
		
	H5close()
	}
	return(gwcTimes)
}
gwc <- wrapGWCtime('/media/schimar/FLAXMAN/h5/', sm[801:1600,], setname= 'sm', folder= 'sm')
gwcSm <- wrapGWCtime('/media/schimar/FLAXMAN/h5/', Sm[801:1599,], setname= 'Sm', folder= 'Sm3')

gwcsM <- wrapGWCtime('/media/schimar/FLAXMAN/h5/', sMsub[[1]], setname= 'sM', folder= 'sM2')




#ccStats.2('Run215500', df, cc500, maf= 0.025)



#-----------------------------------------------------------------------------------------------------------------
#library(Hmisc) 

#arrows(x, avg-sdev, x, avg+sdev, length=0.05, angle=90, code=3)


#lapply(mIgen[[1]][[1]], '[[', 1)      # get the observed Moran's I value 
#
#lapply(mIgen[[1]][[1]], '[[', 3)	  # get the sd







dnearneigh(as.matrix(cbind(fst$MAP[1:115], rep(0, 115))), 0, 25)
# use knearneigh to get the d2 value??? 

# determine number of bins 
nclass.scott(fst$MAP[1:115])	   # 5
nclass.Sturges(fst$MAP[1:115])     # 8    (used in correlog in the ncf package)
nclass.FD(fst$MAP[1:115])   	   # 5

# create distance bins    (why not use kmeans clustering for that? this would give us the indeces for respective groups

x <- kmeans(fst$MAP[1:23], 5)      #, algo= 'Lloyd')

x$cluster 
# then loop through the group number and 

#dists <- mat2listw(as.matrix(dist(cChrom$MAP)))
#mTest <- moran.test(cChrom$Fst, dists, na.action= na.omit, zero.policy= T, rank= T, adjust.n= T)



# # keep separate screens open (with empty plot) 
# 
# for(i in 600:length(fstspl)) {
# 	close.screen(all.screens= T)
# 	#
# 	chrom <- fstspl[[i]]$chromosomeMembership
# 	#par(oma=c(4.5,4.5,4.5,4.5), mar=c(4.0,4.0,4.0,4.0) + 0.2, mgp= c(3,1,0))
# 	split.screen(c(1,4))
# 	screen(1)
# 	plot(0:10, seq(0,1, 0.1), xlim= c(0, 25), ylim= c(0,1), type= 'n', xlab= 'map location (cM)', ylab= expression('F'[ST])) 
# 	screen(2)
# 	plot(0:10, seq(0,1, 0.1), xlim= c(0, 25), ylim= c(0,1), type= 'n', xlab= 'map location (cM)', ylab= expression('F'[ST])) 
# 	screen(3)
# 	plot(0:10, seq(0,1, 0.1), xlim= c(0, 25), ylim= c(0,1), type= 'n', xlab= 'map location (cM)', ylab= expression('F'[ST])) 
# 	screen(4)
# 	plot(0:10, seq(0,1, 0.1), xlim= c(0, 25), ylim= c(0,1), type= 'n', xlab= 'map location (cM)', ylab= expression('F'[ST])) 
# 	#
# 	#screen(1)
# 	#points(fstspl[[i]]$MAP, fstspl[[i]]$Fst, col= as.factor(fstspl[[i]]$locType)
# 	nChrom <- unique(chrom)
# 	for (j in 1:length(nChrom)) {
# 		cChrom <- fstspl[[i]][which(chrom == unique(chrom)[j]),]
# 		screen(j)
# 		points(cChrom$MAP, cChrom$Fst, col= as.factor(cChrom$locType), pch= 20)
# 	}
# }




xLab <- ''


###
# plot numbers of loci (total, s, totalMAF, sMAF)
plot(s500$nLoci$nLoci, type= 'l')
points(s500$nLoci$nLociS, type= 'l', col= 'red')
points(s500$nLoci$nLocimaf, type= 'l', col= 'blue')
points(s500$nLoci$nLociSmaf, type= 'l', col= 'yellow')


### plot screen 3 with proper axes (log10-scale...)
plot(log10(s500$phiObs), log= 'y', ylab= expression(paste('log'[10]~ phi)), xlab= xLab, type= 'l', col= 'grey70', axes= F) 
points(log10(s500$kphisMax), type= 'l', col= 'black')
axis(1, at= c(0, 200, 400, 600, 800, 1000, 1200), cex.axis= 1, lwd= 1)

axis(2, at= seq(-2, 12, 2), labels= labels)   # not quite there yet!

#plot(d2, type ="b", log="y",axes=FALSE, ylim=c(1,10^7))
#axis(2, at=10^(0:6), labels=formatC(10^(0:6),format="f", digits=0),
#     cex.axis=0.8,las=2 )
#axis(1, at=1:25, cex.axis=.6)

# from Dan Standage's old blog  (not quite ...) 
ticks <- seq(0, 4, by=2)
labels <- sapply(ticks, function(i) as.expression(bquote(10^ .(i))))
axis(2, at=c(0, 1, 100), labels=labels)

#######################
# calculate GWC time 


mem <- calcMaxEffMig(r499$sBar, 0.05, unlist(lapply(fstspl, length)))[[2]]

gwc499 <- calcGWCtime(em, mem, -1)

# intersections (nGen) of em and afts&fst 

	# if not the same time points, then create emReduced (with only the times that are in fst&afts)

emRed <- em[em$nGen %in% unique(fst$nGen),]

# calcMeanS with geometric mean of s values (and then use sBar, to compare the two?)
	# maybe just calculate & output both? 

sBarfreqs <- singleLocusEq(unlist(r499$sBar), m, singleS= F)

# calcMaxEffMig (nLociS, sBarVector, m) 
	# calculate exp. single-locus equilibrium frequencies in favored deme (and in the non-favored deme)  (peqvec & qvec) (maybe single func)
	
	# calculate fitnesses (of random immigrant and random resident, respectively (maybe single func), patchAvg)
		# if equilibrium frequency == 1, then 1+s_i (nRes = fitness of random resident)  (since randomRes = (1 + (peq * s))^L )

		# randomImm = fitness of random immigrant (randomImm = (1 + (q * s))^L )

	# patchAvgFit <- (randomImm * m)+ ((1 - m) * randomRes)

	# maxEffMig(i) <- (randomImm * m) / patchAvgFit


# calcGWCtime  (calcGWCtime3 in makeDataSummary.m)
	# (emReduced, maxEffMig, df$END_PERIOD_ALLOPATRY)

	# if em[nGen] >= 0.5 * maxEffMig[nGen]      # still close to random
		# gwcTime <- NaN
		# gwc not reached
	# else 
		# find the last time where effMig was greater than the threshold for independent loci 
		# if gwcIndex < 1
			# gwcTime <- end_period_allopatry
		# else: 
			# gwcTime <- nGen[gwcIndex]

















#plot(1, type="n", xlab="", ylab="") #, xlim=c(0, 11000), ylim=c(-0.5, 0.5))


# plot sStar, mean_s, and sBar
plot(1, type="n", xlab="", ylab="", xlim=c(0, 1250), ylim=c(-0.2, 0.1))

points(r499$sStarLeS$sStar, col= 'green', pch= '.', cex= 1.5) #, type= 'l')
points(r499$sStarLeS$s, col= 'grey10', pch= '.', cex= 1.5) # type= 'l')
points(r499$sBar, col= 'orange', pch= '.', cex= 1.5)


legend('topleft', legend= c('mean s', expression(bar(s)), expression('s'^'*')), fill= c('black', 'orange', 'green'))

# set x&ylim values !
# summary(do.call("rbind", r499$sStarLeS))


# plot Le  & me    (on log10 scale)   (or keep as is in plotStatic())

plot(log10(abs(r499$sStarLeS$Le)), type= 'l', xlab= 'gen / 1e3', ylab= expression(paste('L'[e])), col= 'orange', ylim= c(-8, 4))
points(log10(r499$effMig[,2]), type= 'l', col= 'black')

# with log10(Le[zeros]) in effMig = 0:  
g <- log10(r499$sStarLeS$Le)
g[which(g == 'NaN')] <- 0

plot(g, type= 'l', xlab= 'gen / 1e3', ylab= expression(paste('L'[e])), col= 'orange', ylim= c(-8, 4))
points(log10(r499$effMig[,2]), type= 'l', col= 'black')

e  <- r499$effMig[,2]
e[which(e == 0)] <- 1e-5



points(log(unlist(r499$sBar)), type= 'l')

###############    (max) EffMig 
#
plot(r499$maxEffMigSbar, col= 'orange', type= 'l', ylim= c(0, (max(r499$effMig$eme0) + 0.8* max(r499$effMig$eme0))), ylab= 'migration rate', xlab= 'gen / 1e3' )
points(r499$effMig$eme0, col= 'black', type= 'l')
points(r499$maxEffMigMeanS, col= 'green', type= 'l')

abline(h= 0.05, col= 'grey70', lty= 3)
abline(v= r499$gwcTimeSbar$gwcTime/1e3, col= 'orange', lty= 2)
abline(v= r499$gwcTimeMeanS$gwcTime/1e3, col= 'green', lty= 2)

cols <- c('orange', 'green', 'black', 'grey70', 'orange', 'green')
legend('right', legend= c(expression(paste('maxEffMig ', bar(s))), 'maxEffMig meanS', expression(paste('m'[e0])), 'm', expression(paste('gwcTime ', bar(s))), 'gwcTime meanS'), col= cols, lty= c(1,1,1,3,2,2), pch= c(NA, NA, NA, NA, NA, NA), pt.cex= 2, seg.len= 0.9, cex= 1.1)


wrapEffMigPlot <- function(path, data, setname, ...) {
	# function to read individual runs (from vector of runs), calculate CC and plotStatic
	#
	for (i in 1:dim(data)[1]){
		run <- data$run[i]
		path5 <- paste('/runs/', run, sep= '')
		#
		ccObjTmp <- readCCobj(run, setname, path)
		ccTmp <- ccStats(data, ccObjTmp$fst, ccObjTmp$afts, ccObjTmp$LDsel, ccObjTmp$LDneut, ccObjTmp$effMig, run)
		#
		plot(ccTmp$maxEffMigSbar, col= 'orange', type= 'l', ylab= 'migration rate', xlab= 'gen / 1e3' )     #, ylim= c(0, (max(ccTmp$effMig[,2]) + 0.8* max(ccTmp$effMig[,2])))
		points(ccTmp$effMig[,2], col= 'black', type= 'l')
		points(ccTmp$maxEffMigMeanS, col= 'green', type= 'l')
		
		abline(h= 0.05, col= 'grey70', lty= 3)
		abline(v= ccTmp$gwcTimeSbar$gwcTime/1e3, col= 'orange', lty= 2)
		abline(v= ccTmp$gwcTimeMeanS$gwcTime/1e3, col= 'green', lty= 2)
		
		cols <- c('orange', 'green', 'black', 'grey70', 'orange', 'green')
		legend('right', legend= c(expression(paste('maxEffMig ', bar(s))), 'maxEffMig meanS', expression(paste('m'[e0])), 'm', expression(paste('gwcTime ', bar(s))), 'gwcTime meanS'), col= cols, lty= c(1,1,1,3,2,2), pch= c(NA, NA, NA, NA, NA, NA), pt.cex= 2, seg.len= 0.9, cex= 1.1)

	H5close()
	}
}
wrapEffMigPlot('/media/schimar/dapperdata/bu2s/h5/', sM[1:10,], 'sM')


# 
r202971 <- readIn('Run202971', df, 'Sm', path= '/media/schimar/schimar2/bu2s/h5/')

plotStatic(r202971, 'Run202971')







####################################################################################
######
# plot the PHIs ~ cline widths (full coupl, no coupling, and observed)


cWallS <- lapply(r499$clineWallS, unlist)
pBarallS <- lapply(r499$pBarAllS, unlist)
phiOncW <- mapply(rep, r499$phiObs, times= unlist(lapply(cWallS, length)))


####
# afDiffs (N & S) vs phiObs
phiOafDiffneut <- mapply(rep, r499$phiObs, times= unlist(lapply(r499$afDiffN, length)))
phiOafDiffsel <- mapply(rep, r499$phiObs, times= unlist(lapply(r499$afDiffS, length)))

plot(log10(unlist(r499$kruukPhi_sMax)), unlist(r499$pHat), type= 'l', ylim= c(0, 1))
points(log10(unlist(phiOncW)), unlist(pBarallS), pch= '.') #type= 'l')
#points(log10(unlist(phiOncW)), abs(unlist(cWallS)), pch= '.') #type= 'l')


# phiObs & afDiff (S & N)
points(log10(unlist(phiOafDiffneut)), abs(unlist(r499$afDiffN)), pch= '.', col= 'blue')
points(log10(unlist(phiOafDiffsel)), abs(unlist(r499$afDiffS)), pch= '.', col= 'red')

######################




#mydf[ sample( which(mydf$gender=='F'), round(0.2*length(which(mydf$gender=='F')))), ]


# initially set up plots

#dev.set(2)
#par(mfrow= c(1,1))
#
#dev.set(3)
#par(mfrow= c(1,1))

# 
for(i in 1:length(r499$clineWallS)){
	hist(r499$pBarAllS[[i]][[1]], xlim= c(-1,1))
}



# eventually m_e L_e and s* should go in here (or log10(N x m_e) ?)



close.screen(2, all.screens= T)


#dev.set(j+1)
		#hist(afDifflist[[j]], xlim= c(-1,1), ylim= c(0, 1400), main= mainList[j])
		text(0.5, 1300, paste('nGen=', i*1000))
		#hist(LDlist[[j]], main= paste('LD', mainList[j]))
		# plot LD 
		#dev.set(j)
		#hist(afDifflist[[j]], xlim= c(-1,1), ylim= c(0, 1400))


for (i in 1:1174){
	hist(r499$afDiffN[[i]], xlim= c(-1,1), ylim= c(0, 1400))
}


par(bg = "white") # erase.screen() will appear not to work
                  # if the background color is transparent 
                  # (as it is by default on most devices).
split.screen(c(2,2)) # split display into two screens

#split.screen(c(1,3), screen = 2) # now split the bottom half into 3
screen(1) # prepare screen 1 for output
plot(10:-10)
screen(4) # prepare screen 4 for output
plot(10:1)

cols <- c('blue', 'red', 'black') #t(rainbow(3))
plot(r499$FSTs$FSTneut, type= 'l', col= cols[1], ylim= c(0,1), xlim= c(0, length(r499$phiObs)+2))
points(r499$FSTs$FSTsel, col= cols[2], type= 'l')
points(r499$FSTs$FSTtot, col= cols[3], type= 'l', lty= 1)
legend(200, 1, legend= c('neutral', 'selected', 'total'), fill= cols)

plot(r499$phiObs, xlim=  c(0, length(r499$phiObs)+2), ylim= c(0,10005), ylab= expression(paste(phi, ' obs.')), main= expression(paste(phi, ' obs.')), xlab= '', type= 'l')



# multiple windows
#plot(1:1)
#dev.new()
#plot(2,2)
#dev.set(dev.prev()) # go back to first
#title(main="test dev 1")
#
#dev.set(dev.next()) # go to second
#title(main="test dev 2")


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
#r499 <- ccStats(df, fst, afts, 'Run215499')

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
#cc <- ccStats(fst, afts)


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



