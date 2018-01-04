library(data.table)
library(biomaRt)

grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")

attributes <- c("chromosome_name","start_position", "end_position","strand", "ensembl_gene_id", "hgnc_symbol")
results    <- getBM(attributes = attributes, filters='entrezgene', values='6750', mart=grch37)


gene.list.dir <- '/gpfs/milgram/project/holmes/GITHUB_REPOS/GENETICS_REPO/projects/SST_PVALB_gradients/gene_lists'
gene.list.dir <- '/Users/kevinanderson/PHD/PROJECTS/GENETICS_REPO/projects/SST_PVALB_gradients/gene_lists'
sst.list <- read.csv(paste0(gene.list.dir, '/sst_gtex_3cort_genes.txt'), header=F)


top.100 <- as.character(sst.list$V1[1:100])
results <- getBM(attributes=attributes, filters='ensembl_gene_id', values=top.100, mart=grch37)

gcta.dir <- '/gpfs/milgram/project/holmes/Open_Data/DATA_UKBIOBANK/REPOSITORY/GWAS/GCTA/OUT_imp_chr_v2_HaploRef_INFO_0.8_MAF_0.01_HWE_1e6_GENO_0.1_LD_r20.8_n8822/'
bim.file <- paste0(gcta.dir, 'ukb_imp_chrALL_v2_HaploRef_INFO_0.8_MAF_0.01_HWE_1e6_GENO_0.1_LDr2_0.8_n8822_pruned.bim')

bim.file <- paste0(gene.list.dir, '/ukb_imp_chrALL_v2_HaploRef_INFO_0.8_MAF_0.01_HWE_1e6_GENO_0.1_LDr2_0.8_n8822_pruned.bim')

bim.data <- fread(bim.file)

all.snp.dat <- NULL
for (row in 1:nrow(results) ){
  chr   <- results[row,]$chromosome_name
  start.pos <- results[row,]$start_position
  end.pos   <- results[row,]$end_position
  
  chr.match <- which(bim.data$V1 == chr)
  chr.data  <- bim.data[chr.match,]
  match.idxs <- intersect(which(chr.data$V4 >= start.pos-10000), which(chr.data$V4 <= end.pos+10000))
  snp.dat <- chr.data[match.idxs,]
  all.snp.dat <- rbind(all.snp.dat, snp.dat)
}

gtex.regions <- c('Amygdala',
                  'Hippocampus',
                  'Anterior_cingulate_cortex_BA24',
                  'Hypothalamus',
                  'Caudate_basal_ganglia',
                  'Nucleus_accumbens_basal_ganglia',
                  'Cerebellar_Hemisphere',	    
                  'Putamen_basal_ganglia',
                  'Cerebellum',
                  'Cortex',
                  'Substantia_nigra',
                  'Frontal_Cortex_BA9')

split_me <- function(x){ 
  tmp <- strsplit(as.character(x), '[.]')
  return(unlist(tmp)[1])
}
out         <- lapply(rownames(gtex.striatal.expr), split_me)
use_ensembl <- unlist(out)
# ------------------


gtex.dir <- '/gpfs/milgram/project/holmes/GITHUB_REPOS/GENETICS_REPO/data/GTEx/GTEx_Analysis_v7_eQTL'
all.gtex.dat <- NULL
for (reg in gtex.regions){
  print(reg)
  gtex.dat <- fread(paste0(gtex.dir, '/Brain_', reg, '.v7.signif_variant_gene_pairs.txt'))
  out <- lapply(gtex.dat$gene_id, split_me)
  ensembl.base <- unlist(out)
  gtex.dat$ensembl_base <- ensembl.base
  gtex.dat$region <- reg
  all.gtex.dat <- rbind(all.gtex.dat, gtex.dat)
}
paste0(gtex.dir, )

snp_mart = useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp", host="grch37.ensembl.org")
rs1333049 <- getBM(attributes=c('refsnp_id','refsnp_source','chr_name','chrom_start','chrom_end','minor_allele','minor_allele_freq','minor_allele_count','consequence_allele_string','ensembl_gene_stable_id','ensembl_transcript_stable_id'), filters = 'snp_filter', values ="rs1333049", mart = snp_mart)

snp.attributes <- c('refsnp_id','refsnp_source','chr_name','chrom_start','chrom_end','minor_allele','minor_allele_freq','minor_allele_count','consequence_allele_string','ensembl_gene_stable_id','ensembl_transcript_stable_id')

test <- getBM(attributes=snp.attributes, filters='chr_name', values = 22, mart = snp_mart)

rs.ref <- fread('/gpfs/milgram/project/holmes/GITHUB_REPOS/GENETICS_REPO/data/GTEx/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt')

