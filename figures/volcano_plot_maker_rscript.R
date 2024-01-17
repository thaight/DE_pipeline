#!/bin/R
#The following is a script that is used to make volcano plots after running circRNA_DE_annotation
#Version: 1.0
#Date: September 7th, 2023
#Publisher: Travis Haight
#Note: Used in command-line, see volcano_plot_maker_rstudio.R for a Rscript that can be used from Rstudio

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
#Read in command-line arguments
args <- commandArgs(trailingOnly = TRUE)
#Check that arguments were given
if (length(args) == 0) {
  cat("No arguments were provided.\n")
} else {
  project_dir_DE <- args[1] #Path to project1
  project_dir_DE_tmp <- paste0(project_dir_DE, "/DE_analysis/circular")
  compare <- args[2] #Compare located within project1/DE_analysis/compare
  cat("Creating volcano plots on the following project: ", compare, "\n", sep = "")
  input_DE <- paste0(project_dir_DE_tmp, "/", compare, "/", compare, "_analyzed.tsv")
  DE_in <- read.table(file = input_DE, sep = '\t', header = TRUE)
  names(DE_in)[1] <- 'circ_id'
  DE_exon <- DE_in
  DE_intron <- DE_in
  DE_antisense <- DE_in
  DE_intergenic <- DE_in
  DE_exon <- DE_exon[DE_exon$circ_type == 'exon',]
  DE_intron <- DE_intron[DE_intron$circ_type == 'intron',]
  DE_antisense <- DE_antisense[DE_antisense$circ_type == 'antisense',]
  DE_intergenic <- DE_intergenic[DE_intergenic$circ_type == 'intergenic',]
  state <- c("blue", "red", "black")
  names(state) <- c("DOWN", "UP", "NO")
  #Create volcano plot for the total study, for whatever reason cant group all the modifications under 1 ggplot call or else it results in blank image
  base_title <- paste0(project, "_", compare)
  DEP1 <- ggplot(data=DE_in, aes(x=logFC, y=-log10(PValue), col=diffexpressed)) + geom_point() + theme_minimal() + 
    geom_text(aes(label=ifelse(diffexpressed != "NO",as.character(GeneID),'')),hjust=0,vjust=0)
  DEP2 <- DEP1 + geom_vline(xintercept=c(-0.6, 0.6), col="red") + geom_hline(yintercept=-log10(0.1), col="red")
  DEP3 <- DEP2 + scale_colour_manual(values = state)
  base_title <- paste0(project, "_", compare)
  DEP3 <- DEP3 + ggtitle(paste0(base_title, " total circRNA expression")) + theme(plot.title = element_text(hjust = 0.5))
  #Create Volcano Plot for each circRNA type
  DEP1_exon <- ggplot(data=DE_exon, aes(x=logFC, y=-log10(PValue), col=diffexpressed)) + geom_point() + theme_minimal() + 
    geom_text(aes(label=ifelse(diffexpressed != "NO",as.character(GeneID),'')),hjust=0,vjust=0)
  DEP2_exon <- DEP1_exon + geom_vline(xintercept=c(-0.6, 0.6), col="red") + geom_hline(yintercept=-log10(0.1), col="red")
  DEP3_exon <- DEP2_exon + scale_colour_manual(values = state)
  DEP3_exon <- DEP3_exon + ggtitle("Exonic circRNA expression") + theme(plot.title = element_text(hjust = 0.5))
  DEP1_intron <- ggplot(data=DE_intron, aes(x=logFC, y=-log10(PValue), col=diffexpressed)) + geom_point() + theme_minimal() + 
    geom_text(aes(label=ifelse(diffexpressed != "NO",as.character(GeneID),'')),hjust=0,vjust=0)
  DEP2_intron <- DEP1_intron + geom_vline(xintercept=c(-0.6, 0.6), col="red") + geom_hline(yintercept=-log10(0.1), col="red")
  DEP3_intron <- DEP2_intron + scale_colour_manual(values = state)
  DEP3_intron <- DEP3_intron + ggtitle("Intronic circRNA expression") + theme(plot.title = element_text(hjust = 0.5))
  DEP1_antisense <- ggplot(data=DE_antisense, aes(x=logFC, y=-log10(PValue), col=diffexpressed)) + geom_point() + theme_minimal() + 
    geom_text(aes(label=ifelse(diffexpressed != "NO",as.character(GeneID),'')),hjust=0,vjust=0)
  DEP2_antisense <- DEP1_antisense + geom_vline(xintercept=c(-0.6, 0.6), col="red") + geom_hline(yintercept=-log10(0.1), col="red")
  DEP3_antisense <- DEP2_antisense + scale_colour_manual(values = state)
  DEP3_antisense <- DEP3_antisense + ggtitle("Antisense circRNA expression") + theme(plot.title = element_text(hjust = 0.5))
  DEP1_intergenic <- ggplot(data=DE_intergenic, aes(x=logFC, y=-log10(PValue), col=diffexpressed)) + geom_point() + theme_minimal() + 
    geom_text(aes(label=ifelse(diffexpressed != "NO",as.character(GeneID),'')),hjust=0,vjust=0)
  DEP2_intergenic <- DEP1_intergenic + geom_vline(xintercept=c(-0.6, 0.6), col="red") + geom_hline(yintercept=-log10(0.1), col="red")
  DEP3_intergenic <- DEP2_intergenic + scale_colour_manual(values = state)
  DEP3_intergenic <- DEP3_intergenic + ggtitle("Intergenic circRNA expression") + theme(plot.title = element_text(hjust = 0.5))
  #Create Volcano Plot Figures
  figures_path <- paste0(project_dir_DE, "/figures/")
  study_name <- paste0(compare, "_volcano.pdf")
  split_name <-  paste0(compare, "_splitbytype_volcano.pdf")
  study_path <- paste0(figures_path, study_name)
  split_path <- paste0(figures_path, split_name)
  pdf(file = study_path, width = 12, height = 8)
  stress<- grid.arrange(DEP3)
  dev.off()
  pdf(file = split_path, width = 12, height = 8)
  stress<- grid.arrange(arrangeGrob(DEP3_exon, DEP3_intron, ncol = 2), arrangeGrob(DEP3_intergenic, DEP3_antisense, ncol = 2), nrow = 2, top =  base_title)
  dev.off()
}
