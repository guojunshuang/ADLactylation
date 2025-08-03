rm(list=ls())
gc()
library(ggplot2)
library(foreach)
library(data.table)
library(gwasvcf)
library(gwasglue)
library(TwoSampleMR)
library(VariantAnnotation)
library(MRcML)

input_folder <- "~/MR/Lactate_CNS/Lactylation_AD"
output_folder <- "~/MR/Lactate_CNS/cMLMA/Lactylation_AD"

file_list <- list.files(path = input_folder, pattern = "harmonise\\.txt$", full.names = TRUE)


for (file_path in file_list) {
  
  tryCatch({

    harm_rt <- fread(file_path, sep = "\t", header = TRUE)
    
    harm_rt <- as.data.frame(harm_rt)
    
    harm_rt <- harm_rt[, !colnames(harm_rt) %in% "samplesize.outcome"]
    harm_rt <- as.data.frame(harm_rt)
    
    harm_rt$samplesize.exposure <- 35287
    
    MRcML_data1 <- data.frame(
      SNP = harm_rt$SNP,
      a1 = harm_rt$effect_allele.exposure,
      a2 = harm_rt$other_allele.exposure,
      a1_freq = harm_rt$eaf.exposure,
      bzx_n = harm_rt$samplesize.exposure,
      b_exp = harm_rt$beta.exposure,            
      se_exp = harm_rt$se.exposure,      
      bzx_pval = harm_rt$pval.exposure, 
      b_out = harm_rt$beta.outcome, 
      bzy_n = harm_rt$samplesize.outcome,
      se_out = harm_rt$se.outcome,       
      bzy_pval = harm_rt$pval.outcome
    )
    
    MRcML_data <- MRcML_data1
    
    cML_result <- mr_cML(
      MRcML_data$b_exp,
      MRcML_data$b_out,
      MRcML_data$se_exp,
      MRcML_data$se_out,
      n = harm_rt$samplesize.exposure,
      random_start = 100,
      random_seed = 1
    )
    
    specific_marker <- sub("_harmonise.txt", "", basename(file_path))
    
    result_output_path <- file.path(output_folder, paste0(specific_marker, "_cML.txt"))
    capture.output(cML_result, file = result_output_path)
    
  }, error = function(e) {

    message("Error processing file: ", file_path)
    message("Error details: ", e$message)
  })
}
