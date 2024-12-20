library(tidyverse)
library(taxonomizr)

# Import kraken data from a directory and parse into a single dataframe, optionally accounting for replicates
# Requires directory path (string), report file pattern (regex), and optional replicate (boolean) + report pattern (regex)
read_reports <- function(directory, report_pattern, replicate = FALSE, replicate_pattern = NULL) {
  
  # Validate inputs
  if (missing(directory)) {
    stop("The 'directory' argument is required.")
  }
  if (missing(report_pattern)) {
    stop("The 'report_pattern' argument is required.")
  }
  if (replicate && is.null(replicate_pattern)) {
    stop("The 'replicate_pattern' argument is required when 'replicate = TRUE'.")
  }

  # Read in the reports from the directory
  myfiles <- list.files(directory, pattern = report_pattern, full.names = TRUE)
  
  # Process each file
  reports_list <- lapply(myfiles, function(file) {
    report <- read.table(file, header = FALSE, sep = "\t", quote = "", strip.white = TRUE)
    col_names <- c("kraken_abundance", "reads_inclusive", "reads_direct", 
                   "minimizer_count", "minimizers_distinct", "ID_level", "taxID", "ID")
    names(report) <- col_names
    
    # Add sample_ID column (based on file name)
    report$sample_ID <- gsub(report_pattern, "", basename(file)) # Use gsub pattern based on your needs
    
    return(report)
  })
  
  # Optionally clean up replicates using replicate pattern
  if (replicate) {
    reports_list <- lapply(reports_list, function(df) {
      df$replicate <- df$sample_ID
      df$replicate <- gsub(replicate_pattern, "", df$replicate) # Format the replicate names
      return(df)
    })
  }
  
  # Combine all reports into a single data frame
  reports <- bind_rows(reports_list)
  
  return(reports)
}

test <- read_reports(directory, report_pattern = "_S[0-9]+_nt_report.txt", replicate = TRUE, replicate_pattern = "[a|b]")

# Produce an overview of the highest abundance taxa for each sample_ID (or replicate)
# Can be used at any stage of filtering to get a sense of the data
# Requires reports dataframe and optional grouping variable, default is sample_ID
overview_sum <- function(reports, grouping = "sample_ID", metric = c("abun", "readsI", "readsE", "minC", "minD"), value, add_rank = TRUE) {
  
  # Validate the grouping variable
  valid_grouping <- c("sample_ID", "replicate")
  if (!grouping %in% valid_grouping) {
    stop("Invalid 'grouping' argument. Must be either: 'sample_ID' or 'replicate'.")
  }
  
  # Validate the metric
  metric <- match.arg(metric)
  
  # Define the corresponding column name based on the selected metric
  metric_column <- switch(metric,
                          abun = "kraken_abundance",
                          readsI = "reads_inclusive",
                          readsE = "reads_direct",
                          minC = "minimizer_count",
                          minD = "minimizers_distinct")
  
  # Check if the specified metric column exists in the reports data frame
  if (!metric_column %in% colnames(reports)) {
    stop(paste("The specified metric:", metric_column, "does not exist in the reports data frame."))
  }
  
  # Group by the selected grouping variable and retrieve the top taxa for the specified metric
  data_sum <- reports %>%
    group_by(across(all_of(grouping))) %>%
    slice_max(order_by = .data[[metric_column]], n = value, with_ties = FALSE) %>%  # Get top taxa by the specified metric
    arrange(desc(.data[[metric_column]]))  # Sort taxa by the specified metric
  
  # Add rank column if the add_rank parameter is TRUE
  if (add_rank) {
    data_sum <- data_sum %>% mutate(rank = row_number())  # Add rank column within each group
  }
  
  return(data_sum)
}

check <- overview_sum(test3, metric = "readsI", grouping="sample_ID", 1000, add_rank = FALSE)
#maybe default to report all above the specified threshold if value is not specified
#or maybe this is starting to overlap with the filtering_metric function

# Filter for desired taxonomic rank and optionally, a specific taxa 
# Add something to handle sub-ranks (S1, C22, etc. but caveat that it will cause over counting if multiple levels are pulled)
# REQUIRES A TAXONOMIZR SQL DB
filter_taxrank <- function(reports, rank = "genus", taxa = NULL, taxa_rank = NULL, sql_path = NULL) {
  
  # Validate inputs
  valid_ranks <- c("domain", "kingdom", "phylum", "class", "order", "family", "genus", "species")
  if (!rank %in% valid_ranks) {
    stop("Invalid 'rank' argument. Options are: 'domain', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'.")
  }
  if (!is.null(taxa) && is.null(taxa_rank)) {
    stop("The 'taxa_rank' argument is required when filtering by taxa.")
  }
  if (!is.null(taxa) && is.null(sql_path)) {
    stop("The 'sql_path' argument is required to retrieve taxonomy.")
  }
  
  # Map taxonomic rank to single-letter code (assumes 'ID_level' uses this convention)
  rank_sub <- toupper(substr(rank, 1, 1))
  
  # Filter the reports for the specified rank
  reports_filtered <- reports %>% filter(ID_level == rank_sub)
  
  # If no taxa filtering is requested, return the filtered reports
  if (is.null(taxa)) {
    return(reports_filtered)
  }
  
  # Retrieve unique taxIDs from the reports
  unique_taxIDs <- unique(reports_filtered$taxID)

  # Retrieve the taxonomic hierarchy for these unique taxIDs
  taxonomy <- as.data.frame(getTaxonomy(unique_taxIDs, sqlFile = sql_path, desiredTaxa = taxa_rank))
  taxonomy$taxID <- rownames(taxonomy) %>% as.integer()

  # Integrate taxIDs back into reports_filtered
  reports_filtered <- left_join(reports_filtered, taxonomy, by = "taxID")
  
  # Filter the reports for matching taxIDs & remove that column
  final_filtered_reports <- test3 %>% filter(.data[[taxa_rank]] %in% taxa)
  final_filtered_reports <- final_filtered_reports %>% select(-taxa_rank)
  
  return(final_filtered_reports)
}

test3 <- filter_taxrank(test, rank = "genus", taxa_rank = "kingdom", taxa = "Fungi", sql_path = "/Users/elizabeth.hunter/rworking/taxonomizr_20240126/accessionTaxa.sql")

#filter using a minimizer cutoff/ratio/min read count 
filter_metrics <- function(input_report, metric = c("readsI", "readsE", "minC", "minD", "abun"), value) {
  # Match the metric argument to ensure it's valid
  metric <- match.arg(metric)
  
  # Define the corresponding column names based on the metric selected
  metric_column <- switch(metric, 
                          readsI = "reads_inclusive", 
                          readsE = "reads_direct", 
                          minC = "minimizer_count", 
                          minD = "minimizers_distinct",
                          abun = "kraken_abundance")
  
  if (!metric_column %in% colnames(reports)) {
    stop(paste("The specified metric:", metric_column, "does not exist in the reports data frame."))
  }
  
  # Filter the input_report based on the specified metric and value
  filtered_report <- input_report %>%
    filter(.data[[metric_column]] > value)  # Adjust the condition as needed (>=, >, <=, etc.)
  
  return(filtered_report)
}

test4 <- filter_metrics(test3, metric = "abun", value = 0)

#create a count table for differential abundance analysis 
library(dplyr)
library(tidyr)

kraken_count <- function(kraken_data, metric = c("readsI", "readsE", "minC", "minD")) {
  
  # Validate the metric
  metric <- match.arg(metric)
  
  # Define the corresponding column name based on the selected metric
  metric_column <- switch(metric,
                          readsI = "reads_inclusive",
                          readsE = "reads_direct",
                          minC = "minimizer_count",
                          minD = "minimizers_distinct")
  
  # Check if the specified metric column exists in the kraken_data data frame
  if (!metric_column %in% colnames(kraken_data)) {
    stop(paste("The specified metric column:", metric_column, "does not exist in the kraken_data data frame."))
  }
  
  # Format the data: select relevant columns and pivot wider
  formatted_data <- kraken_data %>%
    select(sample_ID, taxID, all_of(metric_column)) %>%
    pivot_wider(names_from = taxID, values_from = all_of(metric_column)) %>%
    replace(is.na(.), NA) %>%
    column_to_rownames(var = "sample_ID")
  
  # Transpose the data frame
  formatted_data <- t(formatted_data)
  
  return(formatted_data)
}

formatted_kraken_data <- kraken_count(test4, metric = "readsI")


