#!/bin/R
#The following is a script that is used to make chromosome density plots after running circRNA_DE_annotation
#Version: 1.0
#Date: September 7th, 2023
#Publisher: Travis Haight
#Note: Used in command-line, see chromosome_density_rscript.R for a Rscript that can be used from Rstudio

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
libraries_to_load <- c("ggplot2", "gridExtra", "dplyr", "ggrepel", "tidyr") #list of required libraries
check_and_install_libraries(libraries_to_load)

#Extract chromosome from circ_id, allows chromosome identification to be in separate column later
extract_chromosome <- function(circ_id) {
  chromosome <- sub("chr(\\w+):.*", "\\1", circ_id)
  if (chromosome %in% "X") {
    return("X")
  } else if (chromosome %in% "Y") {
    return("Y")
  } else if (chromosome %in% "M") {
    return("M")
  } else if (chromosome %in% "C") {
    return("C")
  } else {
    return(as.numeric(chromosome))
  }
}
#Assign integer values to not number chromosomes, this allows to sort X, Y, M later on for nicer looking figures
assign_chrom <- function(chromosome) {
  if (chromosome == "X") {
    return(max_chrom + 1)
  } else if (chromosome == "Y") {
    return(max_chrom + 2)
  } else if (chromosome == "M") {
    return(max_chrom + 3)
  } else if (chromosome == "C") {
    return(max_chrom + 3)
  } else {
    return(as.numeric(chromosome))
  }
}

#Read in command-line arguments
args <- commandArgs(trailingOnly = TRUE)
#Check that arguments were given
if (length(args) == 0) {
  cat("No arguments were provided.\n")
} else {
  project_dir_DE <- args[1] #Path to project1
  project_dir_DE_tmp <- paste0(project_dir_DE, "/DE_analysis/circular")
  compare <- args[2] #Compare located within project1/DE_analysis/compare
  logFC_max <- 0.6
  logFC_min <- logFC_max * -1
  cat("Creating chromosome density plots on: ", compare, "\n", sep = "")
  input_DE <- paste0(project_dir_DE_tmp, "/", compare, "/", compare, "_analyzed.tsv")
  DE_in <- read.table(file = input_DE, sep = '\t', header = TRUE)
  names(DE_in)[1] <- 'circ_id'
  #Extract chromosome
  DE_in$chrom <- sapply(DE_in$circ_id, extract_chromosome)
  #Filter out only those that are UP/DOWN regulated
  filter_DE_in <- DE_in %>% filter(diffexpressed %in% c("UP", "DOWN"))
  base_title <- paste0(project, "_", compare)
  #Total density plot
  max_chrom <- max(as.numeric(sub("chr(\\d+).*", "\\1", DE_in$chrom)[!DE_in$chrom %in% c("X", "Y", "M")]), na.rm = TRUE)
  #Assign X,Y,M chromosomes as integer values for temp sorting
  DE_in$chrom <- sapply(DE_in$chrom, assign_chrom)
  #Ensures thats the labels can grow with any chrom # organism
  x_axis_lab <- c(1:max_chrom, "X", "Y", "M", "C")
  #Create density plot
  total_title <- paste0(base_title, " total circRNA density")
  plot <- ggplot(DE_in, aes(x = factor(chrom, levels = c(paste0(1:max_chrom), (max_chrom + 1), (max_chrom + 2), (max_chrom + 3))), y = logFC, color = diffexpressed)) +
    geom_point() +
    scale_color_manual(values = c(UP = "red", DOWN = "blue")) +
    labs(x = "Chromosome", y = "logFC") +
    theme_minimal() +
    scale_x_discrete(labels = x_axis_lab) +
    ggtitle(total_title) + theme(plot.title = element_text(hjust = 0.5))
  #Filter density plot
  max_chrom <- max(as.numeric(sub("chr(\\d+).*", "\\1", filter_DE_in$chrom)[!filter_DE_in$chrom %in% c("X", "Y", "M")]), na.rm = TRUE)
  #Assign X,Y,M chromosomes as integer values for temp sorting
  filter_DE_in$chrom <- sapply(filter_DE_in$chrom, assign_chrom)
  #Ensures thats the labels can grow with any chrom # organism
  x_axis_lab <- c(1:max_chrom, "X", "Y", "M")
  #Create density plot
  filter_title <- paste0(base_title, " dsyregulated circRNA density")
  filter_plot <- ggplot(filter_DE_in, aes(x = factor(chrom, levels = c(paste0(1:max_chrom), (max_chrom + 1), (max_chrom + 2), (max_chrom + 3))), y = logFC, color = diffexpressed)) +
    geom_point() +
    scale_color_manual(values = c(UP = "red", DOWN = "blue")) +
    labs(x = "Chromosome", y = "logFC") +
    theme_minimal() +
    scale_x_discrete(labels = x_axis_lab) +
    ggtitle(filter_title) + theme(plot.title = element_text(hjust = 0.5)) +
    geom_rect(xmin = -Inf, xmax = Inf, ymin = logFC_min, ymax = logFC_max,
              fill = "grey", alpha = 0.2, color = "grey")
  #Create density plot Figures
  figures_path <- paste0(project_dir_DE, "/figures/")
  study_name <- paste0(compare, "_chromosome_density.pdf")
  filter_name <-  paste0(compare, "_filtered_chromosome_density.pdf")
  study_path <- paste0(figures_path, study_name)
  filter_path <- paste0(figures_path, filter_name)
  pdf(file = study_path, width = 12, height = 8)
  stress<- grid.arrange(plot)
  dev.off()
  pdf(file = filter_path, width = 12, height = 8)
  stress<- grid.arrange(filter_plot)
  dev.off()
}
