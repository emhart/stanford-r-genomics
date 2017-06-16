# Simple script to parse a folder of E. coli VCF files
# and produce a simple plot of the proportion of variants per
# position along the chromosome
# Naupaka Zimmerman
# naupaka@gmail.com
# June 12, 2017

# Note: might need some fiddling to get vcfR to install
# because of dependencies of dependencies (should all be on CRAN though)

# See great tutorial and documentation for vcfR here:

# https://knausb.github.io/vcfR_documentation/index.html

# Knaus, Brian J., and Niklaus J. Grunwald. 2017. 
# VCFR: a package to manipulate and visualize variant 
# call format data in R. Molecular Ecology Resources 
# 17(1):44-53. http://dx.doi.org/10.1111/1755-0998.12549.

# install.packages("vcfR", dependencies = TRUE)
library("vcfR")
library("tidyr")
library("dplyr")

# Get list of vcf files in current working directory
# These are produced by running the script 
# run_variant_calling.sh 
# as described here:
# https://github.com/JasonJWilliamsNY/wrangling-genomics/blob/gh-pages/lessons/02-variant-calling-workflow.md

my_vcf_files <- list.files(path = "./",
                           pattern = "\\.vcf")

number_of_strains <- length(my_vcf_files)

# Initialize a 9 column empty matrix
all_vcf <- matrix(NA, 0, 9)

# Loop over files and append to bottom of matrix
# Add value for strain ID to first column each time
for (vcf_file in my_vcf_files) {
  
  # read vcf file using function from package "vcfR"
  my_vcf_in <- read.vcfR(vcf_file, 
                         verbose = FALSE)
  
  # Pull strain name out of first part fo filename
  strain_name <- strsplit(vcf_file, "_")[[1]][1]
  
  # Bind column of strain name to rest of data for that
  # strain, pulled out of @fix slot in vcf object
  all_this_vcf <- cbind(rep(strain_name, 
                            nrow(my_vcf_in@fix)), 
                        my_vcf_in@fix)
  
  # Add this newly labeled data to bottom of matrix
  all_vcf <- rbind(all_vcf, 
                   all_this_vcf)
}

# Give the strain name column a proper column name
colnames(all_vcf)[1] <- "STRAIN"

# Convert to df for use by dplyr
all_vcf <- as.data.frame(all_vcf)

## This data frame, "all_vcf", has a bunch of metadata in it
## and could be used for other R teaching stuff if you like
write.csv(all_vcf, "all_vcf.csv", row.names = FALSE)

# calculate proportion of strains that had 
# a variant in each position along the chromosome
# note: have to group POS (i.e. position) as a factor
# first, then mutate to numeric for later plotting
variants_each_position <- all_vcf %>%
  group_by(POS) %>%
  summarise(count = n()/number_of_strains) %>%
  mutate(POS = as.numeric(as.character(POS)))

# Plot simple figure showing basic results
plot(x = variants_each_position$POS, 
     y = variants_each_position$count,
     xlab = "Position on chromosome",
     ylab = "Porportion of strains with variation",
     main = "E. coli genome variation",
     ylim = c(0,1))
