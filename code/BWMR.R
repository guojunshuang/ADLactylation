rm(list=ls())
gc()
library(data.table)
library(ggplot2)
setwd("~//BWMR")
source("BWMR_updated.R") 
input_folder <- "~/MR/Lactate_CNS/Lactylation_PD"
output_folder <- "~/MR/Lactate_CNS/BWMR/Lactylation_PD"

file_list <- list.files(input_folder, pattern = "harmonise.txt$")

for (file_name in file_list) {
  file_path <- file.path(input_folder, file_name)
  
  harm_rt <- fread(file_path, sep = "\t", header = TRUE)

  fit.BWMR <- BWMR(gammahat = harm_rt$beta.exposure,                 
                   Gammahat = harm_rt$beta.outcome,                 
                   sigmaX = harm_rt$se.exposure,                 
                   sigmaY = harm_rt$se.outcome) 

  lower95 <- fit.BWMR$beta - 1.96 * fit.BWMR$se_beta 
  upper95 <- fit.BWMR$beta + 1.96 * fit.BWMR$se_beta
  OR <- exp(fit.BWMR$beta)
  OR_lower95 <- exp(lower95)
  OR_upper95 <- exp(upper95)

  BWMR_results <- data.frame(
    beta = fit.BWMR$beta,
    se_beta = fit.BWMR$se_beta,
    P_value = fit.BWMR$P_value,
    lower95 = lower95,
    upper95 = upper95,
    OR = OR,
    OR_lower95 = OR_lower95,
    OR_upper95 = OR_upper95
  )

  specific_marker <- sub("_harmonise.txt", "", file_name)
  output_file <- file.path(output_folder, paste0("BWMR_", specific_marker, ".txt"))

  write.table(BWMR_results, output_file, row.names = FALSE, sep = "\t", quote = FALSE)
}
