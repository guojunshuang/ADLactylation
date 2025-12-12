rm(list=ls())
gc()
library(data.table)

input_folder <- "~/MR/Lactate_CNS/Lactylation_FTD"
output_folder <- "~/MR/Lactate_CNS/GSMR/Lactylation_FTD"

file_list <- list.files(path = input_folder, pattern = "harmonise\\.txt$", full.names = TRUE)

for (file_path in file_list) {
  
  harm_rt <- fread(file_path, sep = "\t", header = TRUE)

  harm_rt <- as.data.frame(harm_rt)
  
  harm_rt <- harm_rt[, !colnames(harm_rt) %in% "samplesize.outcome"]
  harm_rt <- as.data.frame(harm_rt)

  harm_rt$samplesize.exposure <- 35287
  
  gsmr_data1 <- data.frame(
    SNP = harm_rt$SNP,
    a1 = harm_rt$effect_allele.exposure,
    a2 = harm_rt$other_allele.exposure,
    a1_freq = harm_rt$eaf.exposure,
    bzx_n = harm_rt$samplesize.exposure,
    bzx = harm_rt$beta.exposure,            
    bzx_se = harm_rt$se.exposure,      
    bzx_pval = harm_rt$pval.exposure, 
    bzy = harm_rt$beta.outcome, 
    bzy_se = harm_rt$se.outcome,       
    bzy_pval = harm_rt$pval.outcome
  )

  gsmr_data <- gsmr_data1
  ldrho <- diag(nrow(gsmr_data))
  colnames(ldrho) <- rownames(ldrho) <- snp_coeff_id <- snpid <- as.character(gsmr_data$SNP)
  
  bzx <- gsmr_data$bzx       
  bzx_se <- gsmr_data$bzx_se      
  bzx_pval <- gsmr_data$bzx_pval
  bzy <- gsmr_data$bzy             
  bzy_se <- gsmr_data$bzy_se       
  bzy_pval <- gsmr_data$bzy_pval   
  n_ref <- gsmr_data$bzx_n
  
  thresh <- 1e-5 
  gwas_thresh <- thresh              
  multi_snps_heidi_thresh <- 0.01  
  nsnps_thresh <- 3               
  heidi_outlier_flag <- TRUE         
  ld_r2_thresh <- 0.05             
  ld_fdr_thresh <- 0.05            
  gsmr2_beta <- 1    
  
  gsmr_results <- gsmr(
    bzx, bzx_se, bzx_pval, bzy, bzy_se, bzy_pval, 
    ldrho, snp_coeff_id, n_ref, heidi_outlier_flag, 
    gwas_thresh, multi_snps_heidi_thresh, nsnps_thresh, 
    ld_r2_thresh, ld_fdr_thresh, gsmr2_beta
  )
  
  beta <- gsmr_results[["bxy"]]
  beta_lci <- gsmr_results[["bxy"]] - 1.96 * gsmr_results[["bxy_se"]]
  beta_uci <- gsmr_results[["bxy"]] + 1.96 * gsmr_results[["bxy_se"]]
  or <- exp(gsmr_results[["bxy"]])
  or_lci <- exp(gsmr_results[["bxy"]] - 1.96 * gsmr_results[["bxy_se"]])
  or_uci <- exp(gsmr_results[["bxy"]] + 1.96 * gsmr_results[["bxy_se"]])
  pvalue <- gsmr_results[["bxy_pval"]]
  
  result <- cbind(beta, beta_lci, beta_uci, or, or_lci, or_uci, pvalue)
  
  specific_marker <- sub("_harmonise.txt", "", basename(file_path))

  harm_rt_output_path <- file.path(output_folder, paste0(specific_marker, "_GSMR_harm_rt.txt"))
  result_output_path <- file.path(output_folder, paste0(specific_marker, "_GSMR.txt"))
  
  write.table(harm_rt, harm_rt_output_path, sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(result, result_output_path, sep = "\t", quote = FALSE, row.names = FALSE)
  
}
