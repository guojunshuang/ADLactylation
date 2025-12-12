# exposure

library(data.table)
expo_rt1 <- read_exposure_data(
  filename = ".txt",
  sep = "\t",
  snp_col = "rsids",
  beta_col = "Beta",
  se_col = "SE",
  effect_allele_col = "effectAllele",
  other_allele_col = "otherAllele",
  pval_col = "Pval",
  eaf_col ="ImpMAF")
expo_rt <- expo_rt1[expo_rt1$pval.exposure < 5e-6,]
setwd("/home/data/t020511/")
get_plink_exe()
test <- ld_clump(
  
  dat = dplyr::tibble(rsid=expo_rt$SNP,
                      pval=expo_rt$pval.exposure,
                      id=expo_rt$exposure),
  clump_kb=10000, clump_r2=0.001,clump_p = 1,
  bfile = "./1kg/EUR", 
  plink_bin = get_plink_exe() 
)
library(dplyr)
exposure_local_ld <- expo_rt %>% 
  filter(SNP %in% test$rsid)
dim(exposure_local_ld)
write.table(exposure_local_ld,"_5e-6_clumped.txt",row.names = F,sep = "\t",quote = F)
# outcome
outcome_origin <- fread("finngen_R12_AD_AM_EXMORE.gz",
                        sep = "\t",header = T)
colnames(outcome_origin)[colnames(outcome_origin) == "rsids"] <- "SNP"
outcome_origin <- as.data.frame(outcome_origin)
common_snp <- merge(exposure_local_ld,outcome_origin,by="SNP")$SNP
length(common_snp)
outcome_data <- format_data(dat = outcome_origin,
                            type = "outcome", 
                            snps = common_snp, 
                            snp_col = "SNP",
                            beta_col = "beta",
                            se_col = "sebeta",
                            eaf_col = "af_alt",
                            effect_allele_col = "alt",
                            other_allele_col = "ref",
                            pval_col = "pval")
exposure_local_ld$id.exposure <- "4278_14_P4HB_Protein_disulfide_isomerase"
outcome_data$id.outcome <- "finngen_R12_AD_AM_EXMORE"
dat <- harmonise_data(exposure_local_ld, outcome_data)
# dat=dat[-c(5,8),]
res <- mr(dat)
res
mr_pleiotropy <- mr_pleiotropy_test(dat)
singlesnp_res<- mr_singlesnp(dat)
OR = generate_odds_ratios(mr_res <- mr(dat))
mr_heterogeneity <- mr_heterogeneity(dat) 
run_mr_presso <- run_mr_presso(dat,NbDistribution = 1000)

# outcome
vcf_outcome <- readVcf("~/MR/Lactate_CNS/Lactate_AD/ieu-b-2.vcf.gz")
outcome_origin <- gwasvcf_to_TwoSampleMR(vcf = vcf_outcome,type = "outcome")
common_snp <- merge(exposure_local_ld,outcome_origin,by="SNP")$SNP
length(common_snp)
outcome_data <- format_data(dat = outcome_origin,
                            type = "outcome", 
                            snps = common_snp, 
                            snp_col = "SNP",
                            beta_col = "beta.outcome",
                            se_col = "se.outcome",
                            eaf_col = "eaf.outcome",
                            effect_allele_col = "effect_allele.outcome",
                            other_allele_col = "other_allele.outcome",
                            pval_col = "pval.outcome",
                            samplesize_col = "samplesize.outcome",
                            ncase_col = "ncase.outcome",
                            ncontrol_col = "ncontrol.outcome",
                            id_col = "id.outcome")
exposure_local_ld$id.exposure <- "met-d-Lactate"
outcome_data$id.outcome <- "ieu-b-2"
dat <- harmonise_data(exposure_local_ld, outcome_data)
# dat=dat[-c(14),]
res <- mr(dat)
res
mr_pleiotropy <- mr_pleiotropy_test(dat)
singlesnp_res<- mr_singlesnp(dat)
OR = generate_odds_ratios(mr_res <- mr(dat))
mr_heterogeneity <- mr_heterogeneity(dat) 
run_mr_presso <- run_mr_presso(dat,NbDistribution = 1000)

# outcome
load("~/MR/Lactate_CNS/Lactate_AD/UK_AD_as_outcome.RData")
head(UKdata)
outcome_origin <- UKdata
common_snp <- merge(exposure_local_ld,outcome_origin,by="SNP")$SNP
length(common_snp)
outcome_data <- format_data(dat = outcome_origin,
                            type = "outcome", 
                            snps = common_snp, 
                            snp_col = "SNP",
                            beta_col = "beta.outcome",
                            se_col = "se.outcome",
                            eaf_col = "eaf.outcome",
                            effect_allele_col = "effect_allele.outcome",
                            other_allele_col = "other_allele.outcome",
                            pval_col = "pval.outcome",
                            samplesize_col = "samplesize.outcome",
                            ncase_col = "ncase.outcome",
                            ncontrol_col = "ncontrol.outcome",
                            id_col = "id.outcome")
exposure_local_ld$id.exposure <- "2524_56_HMGB1_HMG_1"
outcome_data$id.outcome <- "UKAD"
dat <- harmonise_data(exposure_local_ld, outcome_data)
res <- mr(dat)
mr_pleiotropy <- mr_pleiotropy_test(dat)
singlesnp_res<- mr_singlesnp(dat)
OR = generate_odds_ratios(mr_res <- mr(dat))
mr_heterogeneity <- mr_heterogeneity(dat) 
run_mr_presso <- run_mr_presso(dat,NbDistribution = 1000)
