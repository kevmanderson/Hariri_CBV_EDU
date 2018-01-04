
# Read data
# --------------
base.dir <- '/Users/kevinanderson/PHD/PROJECTS/FOR_HARIRI'

pheno.data.in <- read.csv(paste0(base.dir, '/GSP_MegaOverview_140124.csv')) # phenotype data with IQ + covariates
snp.scores    <- read.csv(paste0(base.dir, '/figures/EduYears_Main_SCORES_AT_ALL_THRESHOLDS.txt'), sep=' ') # Okbey PGS scores
brain.pheno   <- read.csv(paste0(base.dir, '/GSP_ICV_vals.txt'), sep=',') # FS derived values
eigen.vals    <- read.csv(paste0(base.dir, '/pca_n1735_imputed.eigenvec'), sep=' ', header=FALSE) # eigenvector PCs (10)
colnames(eigen.vals) <- c('Fam1', 'Fam2', paste('eigen',1:10, sep='_'))
# ---------------------
# ---------------------


# COMBINE SNP/EIGENVECTOR DATA
# ---------------------
if (length(which(eigen.vals$Fam1 == snp.scores$IID)) != length(snp.scores$IID)){
  print('check subject order before combining')
} else {
  print('all good')
}
# Combine PGS scores and eigenvectors
snp.scores <- cbind(snp.scores, eigen.vals)

# create residualized eigenvector values by regressing out PGS
for (i in 1:10){
  snp.scores[[paste0("eigen_resid_", i)]] <- scale(resid(lm(paste("eigen_", i," ~ pT_1.1",sep = ""), data=snp.scores)))
}
# ---------------------
# ---------------------



# COMBINE PHENOTYPE/BRAIN DATA
# ---------------------
pheno.data       <- pheno.data.in[which(pheno.data.in$Label %in% brain.pheno$subj),]
pheno.data.sort1 <- pheno.data[order(pheno.data$Label),] # order pheno data
brain.pheno.sort <- brain.pheno[order(brain.pheno$subj),] # order brain metric data
pheno.data.sort1 <- cbind(pheno.data.sort1, brain.pheno.sort)

# check that the row order of each dataset lines up before combining
if (length(as.character(pheno.data.sort1$Label) == brain.pheno.sort$subj) != length(brain.pheno.sort$subj)){
  print('check subject order before combining')
} else {
  print('all good')
}
# Combine phenotype/brain data
pheno.data.sort1 <- cbind(pheno.data.sort1, brain.pheno.sort)
# ---------------------
# ---------------------


# COMBINE ALL DATA
# ---------------------
pheno.data.sort       <- pheno.data.sort1[order(pheno.data.sort1$FamID),]
pheno.data.sort$FamID <- as.character(pheno.data.sort$FamID)
snp.scores.sort       <- snp.scores[order(snp.scores$IID),]
snp.scores.sort$IID   <- as.character(snp.scores.sort$IID)
if (length(which(pheno.data.sort$FamID == snp.scores.sort$IID)) != length(snp.scores.sort$IID)){
  print('check subject order before combining')
} else {
  print('all good')
}
pheno.data.sort <- cbind(pheno.data.sort, snp.scores.sort)
# get rid of 1 subject with 'Drop' GWAS indicator
pheno.data.sort <- pheno.data.sort[pheno.data.sort$GWAS_Trim == 'Pass',]
# ---------------------
# ---------------------


# RUN MODELS - ALL SUBJECTS
# ---------------------
# TBV predicting IQ
# Subjects with compliant Shipley IQ measures
pheno.data.ship <- pheno.data.sort[which(pheno.data.sort$Shipley_SM_Comp == 'Compliant'),]
summary(lm(scale(EstIQ_Ship) ~ scale(BrainSegNotVent) + MF + scale(Age_A) + Scanner + Console + Coil, data = pheno.data.ship))
# GSP-determined Anatomical cutoff (more stringent than TBV-only cutoff)
pheno.data.ship.gspAnat <- pheno.data.ship[which(pheno.data.ship$Anat_Trim == 'Keep'),]
summary(lm(scale(EstIQ_Ship) ~ scale(BrainSegNotVent) + MF + scale(Age_A) + Scanner + Console + Coil, data = pheno.data.ship.gspAnat))


# TBV PGS predicting IQ 
# (eigenvectors have already been scaled above)
summary(lm(scale(EstIQ_Ship) ~ scale(pT_1.1) + MF + scale(Age_A) + Scanner + Console + Coil + 
             eigen_resid_1 + eigen_resid_2 + eigen_resid_3 + eigen_resid_4 + eigen_resid_5 + 
             eigen_resid_6 + eigen_resid_7 + eigen_resid_8 + eigen_resid_9 + eigen_resid_10, data = pheno.data.ship))
# GSP-determined Anatomical cutoff (more stringent than TBV-only cutoff)
summary(lm(scale(EstIQ_Ship) ~ scale(pT_1.1) + MF + scale(Age_A) + Scanner + Console + Coil + 
             eigen_resid_1 + eigen_resid_2 + eigen_resid_3 + eigen_resid_4 + eigen_resid_5 + 
             eigen_resid_6 + eigen_resid_7 + eigen_resid_8 + eigen_resid_9 + eigen_resid_10, data = pheno.data.ship.gspAnat))


# TBV PGS/TBV mediation
summary(lm(scale(EstIQ_Ship) ~  scale(BrainSegNotVent) + scale(pT_1.1) + 
             MF + scale(Age_A) + Scanner + Console + Coil +
             eigen_resid_1 + eigen_resid_2 + eigen_resid_3 + eigen_resid_4 + eigen_resid_5 + 
             eigen_resid_6 + eigen_resid_7 + eigen_resid_8 + eigen_resid_9 + eigen_resid_10, data = pheno.data.ship))
# GSP-determined Anatomical cutoff (more stringent than TBV-only cutoff)
summary(lm(scale(EstIQ_Ship) ~  scale(BrainSegNotVent) + scale(pT_1.1) + 
             MF + scale(Age_A) + Scanner + Console + Coil +
             eigen_resid_1 + eigen_resid_2 + eigen_resid_3 + eigen_resid_4 + eigen_resid_5 + 
             eigen_resid_6 + eigen_resid_7 + eigen_resid_8 + eigen_resid_9 + eigen_resid_10, data = pheno.data.ship.gspAnat))


# PGS predicting TBV
summary(lm(scale(BrainSegNotVent) ~ scale(pT_1.1) + 
             MF + scale(Age_A) + Scanner + Console + Coil +
             eigen_resid_1 + eigen_resid_2 + eigen_resid_3 + eigen_resid_4 + eigen_resid_5 + 
             eigen_resid_6 + eigen_resid_7 + eigen_resid_8 + eigen_resid_9 + eigen_resid_10, data = pheno.data.ship))
# GSP-determined Anatomical cutoff (more stringent than TBV-only cutoff)
summary(lm(scale(BrainSegNotVent) ~ scale(pT_1.1) + 
             MF + scale(Age_A) + Scanner + Console + Coil +
             eigen_resid_1 + eigen_resid_2 + eigen_resid_3 + eigen_resid_4 + eigen_resid_5 + 
             eigen_resid_6 + eigen_resid_7 + eigen_resid_8 + eigen_resid_9 + eigen_resid_10, data = pheno.data.ship.gspAnat))


# PGS predicting TBV - with IQ covariate
summary(lm(scale(BrainSegNotVent) ~ scale(pT_1.1) + scale(EstIQ_Ship) +
             MF + scale(Age_A) + Scanner + Console + Coil +
             eigen_resid_1 + eigen_resid_2 + eigen_resid_3 + eigen_resid_4 + eigen_resid_5 + 
             eigen_resid_6 + eigen_resid_7 + eigen_resid_8 + eigen_resid_9 + eigen_resid_10, data = pheno.data.ship.gspAnat))
# ---------------------
# ---------------------



# RUN MODELS - 18-23yo SUBJECTS
# ---------------------
# Subset the data to examine 18-23 yos
young.data.ship         <- pheno.data.ship[intersect(which(pheno.data.ship$Age_A >= 18), which(pheno.data.ship$Age_A <= 23)),]
young.data.ship.gspAnat <- pheno.data.ship.gspAnat[intersect(which(pheno.data.ship.gspAnat$Age_A >= 18), which(pheno.data.ship.gspAnat$Age_A <= 23)),]


# TBV predicting IQ
summary(lm(scale(EstIQ_Ship) ~ scale(BrainSegNotVent) + MF + scale(Age_A) + Scanner + Console + Coil, data = young.data.ship))
# GSP-determined Anatomical cutoff (more stringent than TBV-only cutoff)
summary(lm(scale(EstIQ_Ship) ~ scale(BrainSegNotVent) + MF + scale(Age_A) + Scanner + Console + Coil, data = young.data.ship.gspAnat))


# PGS predicting IQ
summary(lm(scale(EstIQ_Ship) ~ scale(pT_1.1) + MF + scale(Age_A) + Scanner + Console + Coil + 
             eigen_resid_1 + eigen_resid_2 + eigen_resid_3 + eigen_resid_4 + eigen_resid_5 + 
             eigen_resid_6 + eigen_resid_7 + eigen_resid_8 + eigen_resid_9 + eigen_resid_10, data = young.data.ship))
# GSP-determined Anatomical cutoff (more stringent than TBV-only cutoff)
summary(lm(scale(EstIQ_Ship) ~ scale(pT_1.1) + MF + scale(Age_A) + Scanner + Console + Coil + 
             eigen_resid_1 + eigen_resid_2 + eigen_resid_3 + eigen_resid_4 + eigen_resid_5 + 
             eigen_resid_6 + eigen_resid_7 + eigen_resid_8 + eigen_resid_9 + eigen_resid_10, data = young.data.ship.gspAnat))


# PGS/TBV mediation predicting IQ
summary(lm(scale(EstIQ_Ship) ~ scale(pT_1.1) + scale(BrainSegNotVent) + 
             MF + scale(Age_A) + Scanner + Console + Coil +
             eigen_resid_1 + eigen_resid_2 + eigen_resid_3 + eigen_resid_4 + eigen_resid_5 + 
             eigen_resid_6 + eigen_resid_7 + eigen_resid_8 + eigen_resid_9 + eigen_resid_10, data = young.data.ship))
# GSP-determined Anatomical cutoff (more stringent than TBV-only cutoff)
summary(lm(scale(EstIQ_Ship) ~ scale(pT_1.1) + scale(BrainSegNotVent) + 
             MF + scale(Age_A) + Scanner + Console + Coil +
             eigen_resid_1 + eigen_resid_2 + eigen_resid_3 + eigen_resid_4 + eigen_resid_5 + 
             eigen_resid_6 + eigen_resid_7 + eigen_resid_8 + eigen_resid_9 + eigen_resid_10, data = young.data.ship.gspAnat))

# PGS predicting TBV
summary(lm(scale(BrainSegNotVent) ~ scale(pT_1.1) + 
             MF + scale(Age_A) + Scanner + Console + Coil +
             eigen_resid_1 + eigen_resid_2 + eigen_resid_3 + eigen_resid_4 + eigen_resid_5 + 
             eigen_resid_6 + eigen_resid_7 + eigen_resid_8 + eigen_resid_9 + eigen_resid_10, data = young.data.ship))
summary(lm(scale(BrainSegNotVent) ~ scale(pT_1.1) + 
             MF + scale(Age_A) + Scanner + Console + Coil +
             eigen_resid_1 + eigen_resid_2 + eigen_resid_3 + eigen_resid_4 + eigen_resid_5 + 
             eigen_resid_6 + eigen_resid_7 + eigen_resid_8 + eigen_resid_9 + eigen_resid_10, data = young.data.ship.gspAnat))

# PGS predicting TBV, with IQ covariate
summary(lm(scale(BrainSegNotVent) ~ scale(pT_1.1) + scale(EstIQ_Ship) +
             MF + scale(Age_A) + Scanner + Console + Coil +
             eigen_resid_1 + eigen_resid_2 + eigen_resid_3 + eigen_resid_4 + eigen_resid_5 + 
             eigen_resid_6 + eigen_resid_7 + eigen_resid_8 + eigen_resid_9 + eigen_resid_10, data = young.data))
# ---------------------
# ---------------------



# ---------------------
# BMI sanity check
# ---------------------
bmi.snp.scores    <- read.csv(paste0(base.dir, '/figures/GIANTBMI2015_SCORES_AT_ALL_THRESHOLDS.txt'), sep=' ') # Okbey PGS scores
colnames(bmi.snp.scores) <- paste0('bmi_', colnames(bmi.snp.scores))
bmi.snp.scores.sort <- bmi.snp.scores[order(bmi.snp.scores$bmi_IID),]
bmi.snp.scores.sort <- bmi.snp.scores.sort[bmi.snp.scores.sort$bmi_IID %in% pheno.data.sort$FamID,]
if (length(which(pheno.data.sort$FamID == bmi.snp.scores.sort$bmi_IID)) != length(bmi.snp.scores.sort$bmi_IID)){
  print('check subject order before combining')
} else {
  print('all good')
}
pheno.data.sort <- cbind(pheno.data.sort, bmi.snp.scores.sort)
delete.subs <- unique(c(which(is.na(pheno.data.sort$wt)), which(is.na(pheno.data.sort$ht)), 
                        which(pheno.data.sort$ht == 9999), which(pheno.data.sort$wt == 9999)))
pheno.data.bmi <- pheno.data.sort[-delete.subs,]

# convert lbs to kilograms
pheno.data.bmi$wt_kg <- pheno.data.bmi$wt * 0.453592
# convert inches to meters 
pheno.data.bmi$ht_m  <- pheno.data.bmi$ht * 0.0254
# calculate BMI
pheno.data.bmi$bmi <- pheno.data.bmi$wt_kg / (pheno.data.bmi$ht_m^2)

# create residualized eigenvector values by regressing out BMI-PGS
for (i in 1:10){
  pheno.data.bmi[[paste0("bmi_eigen_resid_", i)]] <- scale(resid(lm(paste("eigen_", i," ~ bmi_pT_1.1",sep = ""), data=pheno.data.bmi)))
}

# BMI-PGS predicting BMI
summary(lm(scale(bmi) ~ scale(bmi_pT_1.1) + 
             MF + scale(Age_A) + Scanner + Console + Coil +
             bmi_eigen_resid_1 + bmi_eigen_resid_2 + bmi_eigen_resid_3 + bmi_eigen_resid_4 + bmi_eigen_resid_5 + 
             bmi_eigen_resid_6 + bmi_eigen_resid_7 + bmi_eigen_resid_8 + bmi_eigen_resid_9 + bmi_eigen_resid_10, data = pheno.data.bmi))
# ---------------------
# ---------------------



# ---------------------
# Height sanity check
# ---------------------
hgt.snp.scores    <- read.csv(paste0(base.dir, '/figures/GIANTHEIGHT2014_SCORES_AT_ALL_THRESHOLDS.txt'), sep=' ') # Okbey PGS scores
colnames(hgt.snp.scores) <- paste0('hgt_', colnames(hgt.snp.scores))
hgt.snp.scores <- hgt.snp.scores[order(hgt.snp.scores$hgt_IID),]
hgt.snp.scores <- hgt.snp.scores[hgt.snp.scores$hgt_IID %in% pheno.data.sort$FamID,]
if (length(which(pheno.data.sort$FamID == hgt.snp.scores$hgt_IID)) != length(hgt.snp.scores$hgt_IID)){
  print('check subject order before combining')
} else {
  print('all good')
}
pheno.data.sort <- cbind(pheno.data.sort, hgt.snp.scores)
delete.subs <- unique(c(which(is.na(pheno.data.sort$ht)), which(pheno.data.sort$ht == 9999)))
pheno.data.hgt <- pheno.data.sort[-delete.subs,]


# create residualized eigenvector values by regressing out HGT-PGS
for (i in 1:10){
  pheno.data.hgt[[paste0("hgt_eigen_resid_", i)]] <- scale(resid(lm(paste("eigen_", i," ~ hgt_pT_1.1",sep = ""), data=pheno.data.hgt)))
}

# BMI-PGS predicting BMI
summary(lm(scale(ht) ~ scale(hgt_pT_1.1) + 
             MF + scale(Age_A) + Scanner + Console + Coil +
             hgt_eigen_resid_1 + hgt_eigen_resid_2 + hgt_eigen_resid_3 + hgt_eigen_resid_4 + hgt_eigen_resid_5 + 
             hgt_eigen_resid_6 + hgt_eigen_resid_7 + hgt_eigen_resid_8 + hgt_eigen_resid_9 + hgt_eigen_resid_10, data = pheno.data.hgt))
# ---------------------
# ---------------------







