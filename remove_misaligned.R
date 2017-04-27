# disable scientific notation
options(scipen=999)

# load needed libraries
library("ape")
library("seqinr")

# read in file names
fasta_dir <- file.path(".")
# always check if the pattern matched alignments properly
alignments <- dir(path=fasta_dir, pattern="*-gb")

## Get row means (taxon means)
## from a matrix of uncorrected p-distances
get_p_distance <- function(alignment) {
  # matrix with pairwise identity
  mat <- dist.alignment(alignment, matrix="identity")
  # matrix with uncorrected p-distances
  p_mat <- mat*mat
  # means for each row
  mns <- rowMeans(as.matrix(p_mat), na.rm=T)
  
  return(mns)
}

## Get names and p-distances of any taxa 
## that are 3 standard deviations from the mean 
get_taxa_above_3sd <- function(p_matrix_row_means) {
  # value 3 standard deviations above the mean
  p_above <- mean(p_matrix_row_means) + sd(p_matrix_row_means) * 3
  # names and p-distances for taxa above threshold
  above_3sd <- p_matrix_row_means[p_matrix_row_means > p_above]
  
  return(above_3sd)
}

## Remove taxa from alignment
## Given the alignment name
remove_outlier_taxa <- function(fasta) {
  # define taxa to retain
  taxa_retain <- c("Acromyrmex_echinatior_Genome",
	"Atta_cephalotes_Genome",
	"Camponotus_floridanus_Genome",
	"Cardiocondyla_obscurior_Genome",
	"Harpegnathos_saltator_Genome",
	"Linepithema_humile_Genome",
	"Pogonomyrmex_barbatus_Genome",
	"Solenopsis_invicta_Genome",
	"Vollenhovia_emeryi_Genome")
  # read in alignment
  alignment <- read.alignment(file=fasta, format="fasta")
  # get row means of p-distance matrix
  p_matrix_means <- get_p_distance(alignment)
  # get names of taxa 3 sd above mean
  taxa_remove <- names(get_taxa_above_3sd(p_matrix_means))
  # negation function
  '%ni%' <- Negate('%in%')
  # ignore taxa that are to be retained
  final_remove <- taxa_remove[taxa_remove %ni% taxa_retain]
  no_outliers_name <- paste("no_outliers", fasta, sep="-")
  
  if ( length(final_remove) > 0 ) {
    # collapse names into one string to print message
    t <- paste(final_remove, collapse=" ")
    print(paste(as.character(length(final_remove)), "taxa to be removed from alignment", fasta, ":", t, sep=" "))
    # construct amas command
    cmd <- paste("amas remove -i", fasta, "-f fasta -d aa -x", t, sep=" ")
    print(cmd)
    system(cmd)
  } else {
  	# if no outlier sequences found, rename processed file
  	# but print a message
  	print(paste("No outlier taxa found in file", fasta, sep=" "))
    cmd <- paste("cp ", fasta, " reduced_", fasta, "-out.fas", sep="")
	print(cmd) 
	system(cmd)   
  }
}


lapply(alignments, remove_outlier_taxa)