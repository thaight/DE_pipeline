#!/bin/R
#The following is a script that is used after CIRIquant which annotates the output for graph making
#Version: 1.0
#Date: September 7th, 2023
#Publisher: Travis Haight
#Note: Used in command-line, see circRNA_DE_annotation_rstudio.R for a Rscript that can be used in Rstudio

#Check if required libraries are installed, and if not install and load
check_and_install_libraries <- function(library_list) {
  for (library_name in library_list) {
    if (!requireNamespace(library_name, quietly = TRUE)) {
      # If not installed, install the library
      install.packages(library_name, repos='https://cloud.r-project.org/', quiet = TRUE)
    }
    
    # Load the library
    library(library_name, character.only = TRUE, quietly = TRUE)
  }
}
libraries_to_load <- c("dplyr", "tidyr") #list of required libraries
check_and_install_libraries(libraries_to_load)
#Read in command-line arguments
args <- commandArgs(trailingOnly = TRUE)
#Check that arguments were given
if (length(args) == 0) {
  cat("No arguments were provided.\n")
} else {
  project_dir_DE <- args[1] #Path to project1
  project_dir_DE_tmp <- paste0(project_dir_DE, "/DE_analysis/circular")
  compare <- args[2] #Compare located within project1/DE_analysis/compare
  cat("Annotating circRNA expression on: ", compare, "\n", sep = "")
  input_DE <- paste0(project_dir_DE_tmp, "/", compare, "/", compare, "_circRNA_DE.csv")
  input_bsj <- paste0(project_dir_DE_tmp, "/", compare, "/", compare, "_circRNA_info.csv")
  DE_in <- read.table(file = input_DE, sep = ',', header = TRUE)
  bsj_in <- read.table(file= input_bsj, sep = ',', header = TRUE)
  names(DE_in)[1] <- 'circ_id'
  #Remove Extra Columns
  bsj_in <- bsj_in[,!grepl("gene_type",names(bsj_in))]
  #Add Gene and type information
  DE <- full_join(DE_in, bsj_in)
  # Sets DE as null to start
  DE$diffexpressed <- "NO"
  #logFC > 0.6 and PValue < 0.1, set as "UP" 
  DE$diffexpressed[DE$logFC > 0.6 & DE$PValue <= 0.1] <- "UP"
  #logFC < -0.6 and PValue < 0.1, set as "DOWN"
  DE$diffexpressed[DE$logFC < -0.6 & DE$PValue <= 0.1] <- "DOWN"
  #Set colours to match state of expression
  state <- c("blue", "red", "black")
  names(state) <- c("DOWN", "UP", "NO")
  #Split gene_id so that only 1 gene value is present for downstream labelling
  DE <- DE %>% separate(gene_id, c("GeneID","ExtraID"), sep = ',', extra = "merge", fill = "right")
  #Write out filtered results for any downstream use
  completed_name <- paste0(compare, "_analyzed.tsv")
  completed_path <- paste0(project_dir_DE_tmp, "/", compare, "/", completed_name)
  write.table(DE, file = completed_path, sep = "\t", quote = FALSE, row.names = FALSE)
}
