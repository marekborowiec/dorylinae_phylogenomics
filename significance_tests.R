# disable scientific notation
options(scipen=999)

library(ggplot2)
library(reshape)

# set working directory to where 'all-loci-min114-summary.txt' is
setwd("/media/mlb/Seagate Backup Plus Drive/MLB_UCEs/For_Dryad/Tabular_data/")

# get the data
all_data <- read.table("all_loci_summary.txt", header=T)

# high boot loci
high_boot <- read.csv("high_boot_loci.txt", header=F)

### SUBSETTING DATA ###

# loci that pass phylogeny-corrected test
Homogeneous <- all_data$Phylo_corrected_test >= 0.05
# check if this correctly identifies 379 
table(Homogeneous)["TRUE"]
# add a T/F homogeneous column to the data 
all_data <- cbind(all_data, Homogeneous)

# loci that were included in high bootstrap matrix
High_boot <- all_data$Locus %in% high_boot$V1
table(High_boot)["TRUE"]
all_data <- cbind(all_data, High_boot)

# subset data for slow-evolving loci
slow_evolving <- all_data$Rate_bin_no == 1
non_slow_evolving <- all_data$Rate_bin_no != 1
slow_subset <- all_data[slow_evolving, ]
non_slow_subset <-all_data[non_slow_evolving, ] 

high_boot <- all_data$High_boot == T
hboot_subset <- all_data[high_boot, ]

head(all_data)

### T-TESTS ###

# RCFV

slow_rcfv <- slow_subset$RCFV
non_slow_rcfv <- non_slow_subset$RCFV
hboot_rcfv <- hboot_subset$RCFV

# check if distribution of values is normal
plot(density(slow_rcfv))
plot(density(non_slow_rcfv))
plot(density(hboot_rcfv))

# check if variance ratios are between 0.5-2
slow_v <- var(slow_rcfv)
non_slow_v <- var(non_slow_rcfv)
hboot_v <- var(hboot_rcfv)

rat_slow_non <- slow_v/non_slow_v
rat_slow_hboot <- slow_v/hboot_v
rat_slow_non
rat_slow_hboot

t.test(non_slow_rcfv, slow_rcfv, alternative="greater", var.equal=T)
t.test(hboot_rcfv, slow_rcfv, alternative="greater", var.equal=T)

# Slope of regression

slow_sat <- slow_subset$Slope
non_slow_sat <- non_slow_subset$Slope
hboot_sat <- hboot_subset$Slope

# check if distribution of values is normal
plot(density(slow_sat))
plot(density(non_slow_sat))
plot(density(hboot_sat))

# check if variance ratios are between 0.5-2
slow_v_sat <- var(slow_sat)
non_slow_v_sat <- var(non_slow_sat)
hboot_v_sat <- var(hboot_sat)

rat_slow_non_sat <- slow_v_sat/non_slow_v_sat
rat_slow_non_sat

rat_slow_hboot_sat <- slow_v_sat/hboot_v_sat
rat_slow_non_sat

t.test(slow_sat, non_slow_sat, alternative="greater", var.equal=T)
t.test(slow_sat, hboot_sat, alternative="greater", var.equal=T)

### Plotting RCFV and saturation ###

# Get only data needed and transform by duplicating records for relevant loci
for_melting <- all_data[c("RCFV", "Slope", "Rate_bin_no", "Homogeneous", "High_boot")]
names(for_melting) <- c("RCFV", "Slope", "slow_evolving", "homogeneous", "high_bootstrap")
# Transform all rate bins except first into zeros
for_melting$slow_evolving[for_melting$slow_evolving != 1] <- 0
# Convert logical values into zeros and ones
for_melting$homogeneous <- as.numeric(as.logical(for_melting$homogeneous))
for_melting$high_bootstrap <- as.numeric(as.logical(for_melting$high_bootstrap))
# Reshape and multiply records
melted <- melt(for_melting, id.var=c("RCFV", "Slope"))
melted <- melted[melted$value != 0,]
names(melted) <- c("RCFV", "Regression_slope", "Loci", "value")

# open pdf device
cairo_pdf(filename="rcfv.pdf",    # create SVG file       
          width = 8,        
          height = 6,
          pointsize = 12)  

ggplot(melted, aes(x=Loci, y=RCFV, fill=Loci)) + 
  geom_boxplot() + guides(fill=F)

# close pdf device
dev.off()

# open pdf device
cairo_pdf(filename="saturation.pdf",    # create SVG file       
          width = 8,        
          height = 6,
          pointsize = 12)  

ggplot(melted, aes(x=Loci, y=Regression_slope, fill=Loci)) + 
  geom_boxplot() + guides(fill=F)

# close pdf device
dev.off()