## Combine individual .csv files into one dataset
sim_total=1000 #Number of simulations run
combined_csv_files=NULL
for(i in 1:sim_total)
{
  temp_data=read.csv(paste("results_boot_CI_HT_55", "_", i, ".csv", sep = ""))
  combined_csv_files=rbind(combined_csv_files,temp_data)
}


## Recreation of Coverage Table (code presented for reference simulation in Figure 7) ##

coverage_Wald <- sum(combined_csv_files$coverage_Wald * combined_csv_files$proceed_stg2_CC_2) / sum(combined_csv_files$proceed_stg2_CC_2)
coverage_percentile <- sum(combined_csv_files$coverage_percentile * combined_csv_files$proceed_stg2_CC_2) / sum(combined_csv_files$proceed_stg2_CC_2)
coverage_BC <- sum(combined_csv_files$coverage_BC * combined_csv_files$proceed_stg2_CC_2) / sum(combined_csv_files$proceed_stg2_CC_2)

## Recreation of Power Table (code presented for reference simulation in Figure 8) ##

# Conditional Power #

power_Wald <- sum(combined_csv_files$reject_Wald * combined_csv_files$proceed_stg2_CC_2) / sum(combined_csv_files$proceed_stg2_CC_2)
power_percentile <- sum(combined_csv_files$reject_percentile * combined_csv_files$proceed_stg2_CC_2) / sum(combined_csv_files$proceed_stg2_CC_2)
power_BC <- sum(combined_csv_files$reject_BC * combined_csv_files$proceed_stg2_CC_2) / sum(combined_csv_files$proceed_stg2_CC_2)

# Unconditional Power #

power_Wald <- sum(combined_csv_files$reject_Wald * combined_csv_files$proceed_stg2_CC_2) / nrow(combined_csv_files)
power_percentile <- sum(combined_csv_files$reject_percentile * combined_csv_files$proceed_stg2_CC_2) / nrow(combined_csv_files)
power_BC <- sum(combined_csv_files$reject_BC * combined_csv_files$proceed_stg2_CC_2) / nrow(combined_csv_files)