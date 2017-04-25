# disable scientific notation
options(scipen=999)

# install required libraries
# uncomment if you don't have these installed
#install.packages("ape")
#install.packages("seqinr")
#install.packages("data.table")

# load needed libraries
library("ape")
library("seqinr")
library("data.table")

### SET WORKING DIRECTORY ###

setwd("/media/mlb/Seagate Backup Plus Drive/MLB_UCEs")

#trees_dir <- file.path("./uce-trees-test/")
#aln_dir <- file.path("./uce-upp-min-114-taxa-test/")
trees_dir <- file.path("./gene_trees/")
aln_dir <- file.path("./uce-upp-min-114-taxa/fasta/")

trees_files <- dir(path=trees_dir, pattern="RAxML_bipartitions.uce-*")
aln_files <- dir(path=aln_dir, pattern="*.fasta")

tree_regex <- "RAxML_bipartitions.(uce-[0-9]+)_alignment_masked.phylip-out"
aln_regex <- "(uce-[0-9]+)_alignment_masked.fasta"

### AVERAGE BOOTSTRAP SUPPORT ###

# This code should work with any Newick tree files
# that have some measure of support at nodes, 
# including PP from Bayesian analysis.

# define a function to calculate average support
Avg_support <- function(file) {
  
  # get UCE number from filename
  locus_no <- sub(tree_regex, "\\1", perl=TRUE, x=file)
  # read in tree
  tree <- read.tree(paste(trees_dir,file, sep=""))
  # store support values in a vector
  support <- c(as.numeric(tree$node.label))
  # calculate average support
  avg_supp <- mean(support, na.rm=T)
  return(c(locus_no,avg_supp))
  
}

# loop over all files
average_bootstrap <- lapply(trees_files, Avg_support)
average_bootstrap <- data.frame(matrix(unlist(average_bootstrap), nrow=(length(average_bootstrap)), byrow=T))
colnames(average_bootstrap) <- c("Locus", "Average_bootstrap")

write.csv(average_bootstrap, file="avg_boot_table.csv", quote=FALSE, row.names=FALSE)

### AVERAGE BRANCH LENGTHS ###

# This takes a Newick tree with branch lengths
# and returns the number of tips for each tree,
# and calculates average branch length.

Br_length.trees <- function(file) {
  
  # reads the phylogenetic tree
  tree <- read.tree(paste(trees_dir,file, sep=""))
  # gets number of tips
  no_tips <- length(tree$tip.label)
  # calculate avg branch length
  avg_br_length <- mean(tree$edge.length)
  # get UCE number from filename
  locus_no <- sub(tree_regex, "\\1", perl=TRUE, x=file)
  return(c(locus_no,avg_br_length))
  
}

# loop over all files
br_lengths <- lapply(trees_files, Br_length.trees)
br_lengths <- data.frame(matrix(unlist(br_lengths), nrow=(length(br_lengths)), byrow=T))
colnames(br_lengths) <- c("Locus", "Average_branch_length")

write.csv(br_lengths, file="avg_br_lengths_table.csv", quote=FALSE, row.names=FALSE)


### SATURATION ###

# The code below takes as input FASTA alignments and tree files with branch lengths
# and calculates regression slope and r-squared for each locus
# also saving regression plots for each locus 
# in a newly created 'Saturation_plots' folder

# define function to perform regression, 
# calculate R-squared, and save saturation plots
# needs fasta alignments and 
# corresponding raxml's tree files
Saturation <- function(seq, tree) {
  
  # read in alignment
  alignment <- read.alignment(file=seq, format="fasta")
  # read in tree
  tree <- read.tree(tree)
  # get locus number from filename
  # this has to be different regex from other functions used here
  # because file names here include directory
  locus_no <- sub(".+?(uce-[0-9]+)_alignment_masked.fasta", "\\1", perl=TRUE, x=seq)
  # matrix with pairwise identity
  mat <- dist.alignment(alignment, matrix="identity")
  # matrix with uncorrected p-distances
  p_mat <- mat*mat
  # mean p-distance
  avg_p_mat <- mean(p_mat, na.rm=TRUE)
  # make matrix of pairwise distances in branch lengths from the tree
  cophentr <- cophenetic(tree)  
  # store as matrix objects
  mat_mat <- as.matrix(mat)
  mat_p_mat <- as.matrix(p_mat)
  # order p-distance matrix by names
  mat_p_mat <- mat_p_mat[order(row.names(mat_p_mat)),order(row.names(mat_p_mat))]
  mat_co <- as.matrix(cophentr)
  # order pairwise distances matrix by names
  mat_co <- mat_co[order(row.names(mat_co)),order(row.names(mat_co))]
  # get lower triangulars of both matrices
  branch_dist <- mat_co[lower.tri(mat_co)]
  p_dist <- mat_p_mat[lower.tri(mat_p_mat)]
  # perform simple linear regression
  regress <- lm(p_dist ~ branch_dist)
  # get slope
  slope <- coef(regress)[2]
  # get r-squared
  Rsquared <- summary(regress)$r.squared
  
  # plot branch length pairwise distances on x
  # and uncorrected p-distances on y
  work_dir <- getwd()
  sat_dir <- file.path(work_dir, "saturation_plots/")
  
  # open png file
  png(file=paste(sat_dir, locus_no, "-saturation.png", sep=""), width=600, height=600)
  plot(branch_dist, p_dist)
  # add simple linear regression line
  abline(lm(p_dist ~ branch_dist), col="red")
  # give title as locus number and subtitle as tree length
  title(main=locus_no,sub=paste("Slope: ", round(slope, digits=3), " R-squared: ", round(Rsquared, digits=3), sep=""), cex.sub=1.25)
  # close png file
  dev.off()
  
  return(list(locus_no, avg_p_mat, slope, Rsquared))
  
}

# create a folder for saturation plots
dir.create("./saturation_plots")

# create a table with file names
files_table <- as.data.frame(cbind(paste(aln_dir, aln_files, sep=""), paste(trees_dir, trees_files, sep="")))

alns_table <- as.matrix(files_table$V1)
trees_table <- as.matrix(files_table$V2)

saturation_table <- t(mapply(Saturation, alns_table, trees_table))
saturation_table <- as.data.frame(saturation_table)
colnames(saturation_table) <- c("Locus","Avg_p_dist","Slope","R_squared")
row.names(saturation_table) <- NULL

write.csv(saturation_table, file="saturation_table.csv", quote=FALSE, row.names=FALSE)


### PLOTTING GENE TREES ### 


# Define a function that checks for monophyly;
# If the possible formicoid group is non-monophyletic
# unrooted tree will be plotted

checkMono <- function(tree,clade){
  taxa <- tree$tip.label
  for (taxon in clade) {
    if (!(taxon %in% taxa)) {
      clade <- clade[!clade==taxon]
    }
  }
  len <- length(clade)
  if (len <= 1) {
    # This is for Harpegnathos situation below
    # but should be turned to 'NA' if one is 
    # to check for a monophyly of a clade
    monophyly <- TRUE
  } else {
    monophyly <- is.monophyletic(tree,tips=clade)
  }
  return(monophyly)
}

# Then define a function to root with either Harpegnathos
# or some of the formicoids

rootNonDory <- function(tree) {
    formicoids <- c("Acromyrmex_echinatior_Genome","Atta_cephalotes_Genome","Camponotus_floridanus_Genome",
                                                "Cardiocondyla_obscurior_Genome","Linepithema_humile_Genome",
                                                "Pogonomyrmex_barbatus_Genome","Solenopsis_invicta_Genome",
                                                "Vollenhovia_emeryi_Genome")
    possible_outgroups <- list("Harpegnathos_saltator_Genome", formicoids)  
    for (outgroup in possible_outgroups) {
    taxa <- tree$tip.label
    # check for monophyly
    mono <- checkMono(tree, outgroup)
    if ((mono == TRUE) && (outgroup %in% taxa)) 
    {
      try(tree <- root(phy=tree,outgroup=outgroup,resolve.root=T))
      break
    }
  }
  return(tree)
}

# This will save png files of 600 x 600 pixel plots
# of all unrooted phylogenetic trees with tip labels
# and support values

Plot_trees <- function(file) {
  
  # reads the phylogenetic tree
  tree <- read.tree(paste(trees_dir,file, sep=""))
  rtree <- rootNonDory(tree)
  # extracts plot name (locus) from file name 
  plot_name <- sub(tree_regex, "\\1", perl=TRUE, x=file)
  # open png file
  png(file=paste(tree_plots_dir, plot_name, "-tree.png", sep=""), width=1000, height=1800)
  plot.phylo(rtree, show.node.label=T)
  # give title as locus number and subtitle as tree length
  title(main=plot_name)
  # close png file
  dev.off()
  
}

# create directory for tree plots
work_dir <- getwd()
dir.create("./tree_plots")
tree_plots_dir <- file.path(work_dir, "tree_plots/")

#loop over all files
lapply(trees_files, Plot_trees)

# putting together all the data
dfs <- list(average_bootstrap, br_lengths, saturation_table)

Multmerge <- function(dfs){
  datalist <- lapply(dfs, function(x){data.frame(x)})
  Reduce(function(x,y) {merge(x,y)}, datalist)
}

all_loci_stats <- Multmerge(edfs)
all_loci_stats <- data.frame(lapply(all_loci_stats, as.character), stringsAsFactors=FALSE)

write.csv(all_loci_stats, file="tree_stats_table.csv", quote=FALSE, row.names=FALSE)