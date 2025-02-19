# Krak Pak Toolkit

### Toolkit built to parse and filter kraken and bracken formatted reports 

### Dependencies:
- tidyverse
- taxonomizr

### Functions:

#### read_reports 
Reads in the reports, accounting for the specific naming scheme and whether there are replicates present

Usage: ```read_reports(directory, report_pattern, type = c("default", "bracken", "minimizers")```

Arguments: 
- directory: file path for the directory of reports (required)
- report_pattern: regex formatted pattern for the report naming scheme (required)
- type: type of report (required)
- replicate: boolean value, default false (optional)
- replicate_pattern: regex formatted pattern for the replicate naming scheme (req'd if replicates=TRUE)

Example: 
```read_reports(directory, report_pattern = "_S[0-9]+_nt_report.txt", replicate = TRUE, replicate_pattern = "[a|b]")```

#### overview_sum
Provide an overview of the top hits for each sample (can be used pre and post filtering)

Usage: ```overview_sum(reports, metric = c("readsI", "readsE", "abun", "minC", "minD"), 1000, add_rank=FALSE)```

Arguments:
- reports: the dataframe list of cleaned reports produced by read_reports 
- grouping: "sample_ID" or "replicate" (required)
- metric: metric to order by, such as reads inclusive, reads exclusive, minimizer count, or distinct minimizer count (required)
- value: the number of rows to report (required)
- add_rank: boolean value to add a column with the numerical rank (default=FALSE)

Example: 
```overview_sum(reports, metric = "readsI", grouping="sample_ID", 1000, add_rank = FALSE)```

#### filter_taxrank
Filters for a designated taxonomic group and optionally, a specific taxa (at any taxonomic level). Requires a pre-built taxonmizr database to filter for a specific taxa.

Usage: ```filter_taxrank(reports, rank = c("domain", "kingdom", "phylum", "class", "order", "family", "genus", "species"), taxa_rank = c("domain", "kingdom", "phylum", "class", "order", "family", "genus", "species"), taxa = "desired_taxa", sql_path = "accessionTaxa.sql")```

Arguments:
- reports: the dataframe list of cleaned reports produced by read_reports (required)
- rank: the taxonomic level to report back (required)
- taxa_rank: the taxonomic level to filter for (optional)
- taxa: the taxa to filter for (optional)
- sql_path: path to the taxonomizr sql database (optional)

Example: 
```filter_taxrank(reports, rank = "genus", taxa_rank = "kingdom", taxa = "Fungi", sql_path = "accessionTaxa.sql")```

#### filter_metrics
Filters the results to a user specified minimum based on the desired metric. 

Usage: ```filter_metrics(reports, metric = c("readsI", "readsE", "minC", "minD", "abun"), value)```

Arguments:
- reports: the dataframe list of cleaned reports produced by read_reports (required)
- metric: metric to order by, such as reads inclusive, reads exclusive, minimizer count, or distinct minimizer count (required)
- value: minimum integer value to filter by (required)

Example:
```filter_metrics(reports, metric = "abun", value = 10)```


#### kraken_count 
Create a count table for differential abundance analysis.

Usage: ```kraken_count(reports, metric = c("readsI", "readsE", "minC", "minD"))```

Arguments: 
- reports: the dataframe list of cleaned reports from read_reports, with desired filtering from filter_metrics or filter_taxrank (required)
- metric: metric to incorporate into the count table, such as reads inclusive, reads exclusive, minimizer count, or distinct minimizer count (required)

Example: 
```kraken_count(reports, metric = "readsI")```

#### compare_taxa_counts
Takes a list of named dataframes, and outputs a graph or a table with a summary of the number of taxa identified in each dataframe. This function expects that the same samples are present in each dataframe and will throw an error if this is not the case, but it will still run even if some samples are missing from some analyses. 

Usage: 
```compare_taxa_counts(df_list, output = c("table", "graph"))```

Arguments: 
- df_list: a list of named dataframes to summarise that include taxID and sample_ID columns (required)
- output: the output type desired, a table or a graph (required)

Example:
```
mydata <- list(
  "kraken" = kraken2,
  "bracken" = bracken2
)

mygraph <- compare_taxa_counts(mydata, output="graph")
```
