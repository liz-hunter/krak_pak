# Krak Pack Toolkit

### Toolkit built to parse kraken and bracken formatted reports 

### Dependencies:
- tidyverse
- taxonomizr

### Functions:

#### read_reports 
Reads in the reports, accounting for the specific naming scheme and whether there are replicates present
Usage: read_reports(directory, report_pattern)
Arguments: 
- directory 	file path for the directory of reports (required)
- report_pattern	regex formatted pattern for the report naming scheme (required)
- replicate 	boolean value, default false (optional)
- replicate_pattern	regex formatted pattern for the replicate naming scheme (req'd if replicates=TRUE)
Example: read_reports(directory, report_pattern = "_S[0-9]+_nt_report.txt", replicate = TRUE, replicate_pattern = "[a|b]")

#### overview_sum
Provide an overview of the top hits for each sample (can be used pre and post filtering)
Usage: overview_sum(reports, metric = c("readsI", "readsE", "abun", "minC", "minD") )
Arguments:
- grouping ("sample_ID" or "replicate")
- metric ("readsI", "readsE", "abun", "minC", "minD")
-

Example: overview_sum(reports, metric = "readsI", grouping="sample_ID", 1000, add_rank = FALSE)


#### filter_taxrank




#### filter_metrics


#### kraken_count 


#### kraken_visualiize
