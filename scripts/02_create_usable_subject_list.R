# R code

# extracted freesurfer values
usable.brains <- read.csv('/cluster/nexus/3/users/holmes/Anderson/PROJECTS/HARIRI_ICV_EDU/HARIRI_ICV_EDU_GSP_vals.txt', header=TRUE)

# GSP subject information
ref.table     <- read.csv('/cluster/nexus/3/users/holmes/Anderson/aholmes_2_3_2017_12_32_16.csv', header=TRUE)

# convert visit IDs to family IDs
fam.id.out    <- as.character(ref.table$FamID[which(ref.table$Label %in% usable.brains$subj)])
write.table(x=fam.id.out, '/cluster/nexus/3/users/holmes/Anderson/PROJECTS/HARIRI_ICV_EDU/GSP_usable_subs.txt', row.names=FALSE, col.names=FALSE, quote=FALSE)




# end
