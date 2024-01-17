#!/bin/R
#The following is a script that is used after CIRIquant which calculates some circRNA stats for downstream use
#Version: 1.0
#Date: September 7th, 2023
#Publisher: Travis Haight
#Note: Used in command-line, see circ_stats_rstudio.R for a Rscript that can be used in Rstudio

#Check if required libraries are installed, and if not install and load
check_and_install_libraries <- function(library_list) {
  for (library_name in library_list) {
    if (!requireNamespace(library_name, quietly = TRUE)) {
      # If not installed, install the library
      install.packages(library_name, repos='https://cloud.r-project.org/')
    }
    # Load the library
    library(library_name, character.only = TRUE)
  }
}
libraries_to_load <- c("ggplot2", "gridExtra", "dplyr", "ggrepel", "tidyr") #list of required libraries
check_and_install_libraries(libraries_to_load)
#Read in command-line arguments
args <- commandArgs(trailingOnly = TRUE)
#Check that arguments were given
if (length(args) == 0) {
  cat("No arguments were provided.\n")
} else {
  compare_path <- args[1] 
  compare_table <- read.delim(file = compare_path, sep = "", header = TRUE) #First arg is the compare list
  project_dir <- args[2] #Path to project
  project <- args[3]
  compare_list <- compare_table$Subproject

  #Initialize the stats output
  headers <- c("ID", "number_downregulated", "number_upregulated", "number_common", "number_exonic_circRNA", "number_intronic_circRNA", "number_intergenic_circRNA", "number_antisense_circRNA")
  stats_output <- data.frame(matrix(ncol = length(headers), nrow = 0))
  #colnames(stats_output) <- headers
  cat("Calculating circRNA stats for the following project: ", project, "\n")
  for (compare in compare_list) {
    cat(compare, "\n")
    DE_path <- paste0(project_dir, "/", "DE_analysis/circular/", compare, "/", project, "_", compare, "_circRNA_DE.csv")
    bsj_path <- paste0(project_dir, "/", "DE_analysis/circular/", compare, "/", project, "_", compare, "_circRNA_info.csv")
    DE <- read.table(file =DE_path, sep = ',', header = TRUE)
    names(DE)[1] <- 'circ_id'
    #Remove Extra Columns
    DE_bsj <- read.table(file =bsj_path, sep = ',', header = TRUE)
    #Remove Extra Columns
    DE_bsj <- DE_bsj[,!grepl("gene_type",names(DE_bsj))]
    #Add Gene and type information
    DE <- full_join(DE, DE_bsj)
    #Add Gene and type information
    DE$diffexpressed <- "NO"
    #logFC > 0.6 and PValue < 0.05, set as "UP" 
    DE$diffexpressed[DE$logFC > 0.6 & DE$PValue <= 0.1] <- "UP"
    #logFC < -0.6 and PValue < 0.05, set as "DOWN"
    DE$diffexpressed[DE$logFC < -0.6 & DE$PValue <= 0.1] <- "DOWN"
    DE <- DE %>% separate(gene_id, c("GeneID","ExtraID"), sep = ',', extra = "merge", fill = "right")
    DE_exon <- DE
    DE_intron <- DE
    DE_antisense <- DE
    DE_intergenic <- DE
    DE_exon <- DE_exon[DE_exon$circ_type == 'exon',]
    DE_intron <- DE_intron[DE_intron$circ_type == 'intron',]
    DE_antisense <- DE_antisense[DE_antisense$circ_type == 'antisense',]
    DE_intergenic <- DE_intergenic[DE_intergenic$circ_type == 'intergenic',]
  
    #Determine the number of up and down regulated circRNA
    regulation_counts <- table(DE$diffexpressed)
    type_counts <- table(DE$circ_type)
  
    #Determine how many genes have transcripts that are both up and down regulated
    both <- DE %>%
      group_by(GeneID) %>%
        summarize(DOWN_count = sum(diffexpressed == "DOWN"),
              UP_count = sum(diffexpressed == "UP"))
    both_filter <- both %>%
      filter(DOWN_count > 0, UP_count > 0)
    both_count <- nrow(both_filter)
    down_count <- regulation_counts["DOWN"]
    up_count <- regulation_counts["UP"]
    exon_count <- type_counts["exon"]
    intron_count <- type_counts["intron"]
    intergenic_count <- type_counts["intergenic"]
    antisense_count <- type_counts["antisense"]
    new_row <- c(compare, down_count, up_count, both_count, exon_count, intron_count, intergenic_count, antisense_count)
    stats_output <- rbind(stats_output, new_row)
  }
  stats_file <- paste0(project_dir,"/", project, "_circRNA_stats.tsv")
  cat("Writing stats out to: ", stats_file, "\n")
  colnames(stats_output) <- headers
  write.table(stats_output, file = stats_file, sep = "\t", row.names = FALSE)
}




