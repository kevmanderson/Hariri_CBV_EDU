# Check for ICV outliers in GSP dataset
subs2check <- NULL

# base dir
hariri.dir <- '/Users/kevinanderson/PHD/PROJECTS/Hariri_ICV_EDU'

# Aggregated GSP Freesurfer 5.3 outputs
GSP.data   <- read.csv(paste0(hariri.dir, '/data/HARIRI_ICV_EDU_GSP_vals.txt'))

# GSP phenotype information with some FS 4.5 output
GSP.pheno <- read.csv(paste0(hariri.dir, '/data/GSP_MegaOverview_140124.csv'))
GSP.pheno <- GSP.pheno[GSP.pheno$Label %in% GSP.data$subj,]
GSP.pheno <- GSP.pheno[order(GSP.pheno$Label),]

# subset for white/non-hispanic latino subjects only 
wht.nonHL.idxs <- intersect(which(GSP.pheno$Race == 'WHT'), which(GSP.pheno$Ethn == 'NOT_HL'))
GSP.pheno.race <- GSP.pheno[wht.nonHL.idxs,] 


# Only examine subjects with passing T1/T2 images (determined by GSP QC procedures)
GSP.pheno.race.qc1 <- GSP.pheno.race[which(GSP.pheno.race$Anat_Trim == 'Keep'),]

# subset FS data 
GSP.data.race.qc1 <- GSP.data[which(GSP.data$subj %in% GSP.pheno.race.qc1$Label),]

# check sorting 
matches <- length(which(as.character(GSP.data.race.qc1$subj) == as.character(GSP.pheno.race.qc1$Label)))
if ( matches == length(GSP.data.race.qc1$subj) ){
  print('data good to go')
} else {
  print('data did not sort properly')
}

# quick rename
GSP.pheno <- GSP.pheno.race.qc1
GSP.data  <- GSP.data.race.qc1

# Identify subjects that have ICV values that are vary greatly from older FS 4.5 ICV values
icv.diffs <- GSP.pheno$ICV - GSP.data$IntraCranialVol
hist(icv.diffs)

# subjects with missing 4.5 ICV values
check.me <- as.character(GSP.data$subj[is.na(icv.diffs)])
check.out <- cbind(check.me, rep('noFS_4.5', length(check.me)))
subs2check <- check.out
#
mean.diff <- mean(icv.diffs, na.rm = TRUE)
sd.diff   <- sd(icv.diffs, na.rm = TRUE)
# subjects that are +/- 3 sd away from FS 4.5 values
top.outliers <- which(icv.diffs >= (mean.diff + (sd.diff*3)))
bot.outliers <- which(icv.diffs <= (mean.diff - (sd.diff*3)))
check.me <- as.character(GSP.data$subj[c(top.outliers, bot.outliers)])
check.out <- cbind(check.me, rep('FS_4.5_diff', length(check.me)))
subs2check <- rbind(subs2check, check.out)


# Control for covariates and then look for ICV outliers
GSP.pheno$ICV_5_3_HCP <- GSP.data$IntraCranialVol
residualized.ICV <- resid(lm(ICV_5_3_HCP ~ MF + Age_A + Scanner + Console + Coil, data=GSP.pheno))
# +/- two SDs for outliers
top.outliers.resid <- as.character(GSP.pheno$Label[which(residualized.ICV >= sd(residualized.ICV)*2)])
bot.outliers.resid <- as.character(GSP.pheno$Label[which(residualized.ICV <= -sd(residualized.ICV)*2)])

check.me <- c(top.outliers.resid, bot.outliers.resid)
check.out <- cbind(check.me, rep('ICV_outlier', length(check.me)))
subs2check <- rbind(subs2check, check.out)


# Control for covariates and then look for Ventricle outliers
ventricle.size <- GSP.data$BrainSegVol - GSP.data$BrainSegNotVent
GSP.pheno$Ventricle_5_3_HCP <- ventricle.size
residualized.Ventricle <- resid(lm(Ventricle_5_3_HCP ~ ICV_5_3_HCP + MF + Age_A + Scanner + Console + Coil, data=GSP.pheno))
# +/- two SDs for outliers
top.outliers.resid <- as.character(GSP.pheno$Label[which(residualized.Ventricle >= sd(residualized.Ventricle)*2)])
bot.outliers.resid <- as.character(GSP.pheno$Label[which(residualized.Ventricle <= -sd(residualized.Ventricle)*2)])

check.me <- c(top.outliers.resid, bot.outliers.resid)
check.out <- cbind(check.me, rep('Ventricles_outlier', length(check.me)))
subs2check <- rbind(subs2check, check.out)


# Control for covariates and then look for GM outliers
GSP.pheno$TotalGray_5_3_HCP <- GSP.data$TotalGray
residualized.TotalGray <- resid(lm(TotalGray_5_3_HCP ~ ICV_5_3_HCP + MF + Age_A + Scanner + Console + Coil, data=GSP.pheno))
# +/- two SDs for outliers
top.outliers.resid <- as.character(GSP.pheno$Label[which(residualized.TotalGray >= sd(residualized.TotalGray)*2)])
bot.outliers.resid <- as.character(GSP.pheno$Label[which(residualized.TotalGray <= -sd(residualized.TotalGray)*2)])

check.me <- c(top.outliers.resid, bot.outliers.resid)
check.out <- cbind(check.me, rep('GM_outlier', length(check.me)))
subs2check <- rbind(subs2check, check.out)




# Control for covariates and then look for WM outliers
GSP.pheno$TotalWM_5_3_HCP <- GSP.data$CorticalWhiteMatter
residualized.WM <- resid(lm(TotalWM_5_3_HCP ~ ICV_5_3_HCP + MF + Age_A + Scanner + Console + Coil, data=GSP.pheno))
# +/- two SDs for outliers
top.outliers.resid <- as.character(GSP.pheno$Label[which(residualized.WM >= sd(residualized.WM)*2)])
bot.outliers.resid <- as.character(GSP.pheno$Label[which(residualized.WM <= -sd(residualized.WM)*2)])

check.me <- c(top.outliers.resid, bot.outliers.resid)
check.out <- cbind(check.me, rep('WM_outlier', length(check.me)))
subs2check <- rbind(subs2check, check.out)


# write outlier subjects for manual QC checking
write.df <- as.data.frame(subs2check)
df$V3 <- 1
write.df <- spread(df, 'V2', 'V3')
write.df[is.na(write.df)] <- 0
write.table(x=write.df, file=paste0(hariri.dir, '/MANUAL_CHECK_outlier_subs.csv'), sep=',', quote=FALSE, row.names=FALSE, col.names=TRUE)

qc.info <- read.csv(file=paste0(hariri.dir, '/MANUAL_CHECK_outlier_subs_kma.csv'))
remove.subs <- qc.info$check.me[qc.info$PASS == FALSE]

GSP.pheno.qc <- GSP.pheno[which(!GSP.pheno$Label %in% remove.subs),]
GSP.data.qc <- GSP.data[which(!GSP.data$subj %in% remove.subs),]
write.table(x=GSP.data.qc$subj, file=paste0(hariri.dir, '/postQC_subs.csv'), sep=',', quote=FALSE, row.names=FALSE, col.names=FALSE)



##
## MANUAL QC happens here, subjects are marked as pass/fail
##

QC.visual.check <- read.csv(file=paste0(hariri.dir, '/outlier_subs_manual_check.csv'))
exclude.subs    <- as.character(QC.visual.check$Label[which(QC.visual.check$QC == 'FAIL')])

GSP.pass.QC <- GSP.data[!GSP.data$subj %in% exclude.subs,]
write.table(x=as.charactre(GSP.pass.QC$subj), file=paste0(hariri.dir, '/outlier_subs.csv'), quote=FALSE, row.names=FALSE, col.names=FALSE)
