## code to prepare `EYR` dataset goes here

# usethis::use_data(DATASET, overwrite = TRUE)

# Read raw data as genlight object
library(dartR)
dir <- "C:/Users/Diangie/My Drive/Projects/SwissArmyKnife/SCRIPTS/eyr/"
dart_file <- "Report_DYro21-6107_5_moreOrders_SNP_mapping_3_uniqIDs.csv"

EYR <- gl.read.dart(filename = paste(dir,
                                     dart_file, 
                                     sep = ""),
                    covfilename = paste(dir,
                                        "EYR2021_cov.NONAS.csv",
                                        sep = ""))

EYR  # 899 genotypes,  54,823 binary SNPs

# Remove individuals to not be used
EYR <- gl.drop.pop(EYR, pop.list = c("DoNotUse",
                                     "KeilorsTrack",
                                     "PileOfBricks"),
                   recalc = TRUE,
                   mono.rm = TRUE) 

EYR  # 782 genotypes,  53,324 binary SNPs

EYR@other$ind.metrics <- droplevels(EYR@other$ind.metrics)


# Change pops from ind.metrics
EYR@other$ind.metrics[EYR@other$ind.metrics$pop == "RailwayDam", "pop"] <- "Muckleford"
EYR@other$ind.metrics[EYR@other$ind.metrics$pop == "Sedgwick", "pop"] <- "Crusoe"

EYR@other$ind.metrics <- droplevels(EYR@other$ind.metrics)


# Change pops from @pops
EYR@pop[EYR@pop == "RailwayDam"] <- "Muckleford"
EYR@pop[EYR@pop == "Sedgwick"] <- "Crusoe"

EYR@pop <- droplevels(EYR@pop)


# Remove secondaries
# Keep only one SNP per contig to control for physical linkage. There are two 
# methods. Method "best" keeps SNP based on repeatability and avgPIC (average 
# Polymorphism Information Criteria, in that order). PIC works like this: for 
# example, a marker that reveals six alleles but one allele is found to be in 
# very high frequency, has a lower discriminatory capacity than another with six 
# alleles but with similar frequencies (Smith et al., 1997). Therefore, markers 
# with PIC values above 0.5 are more recommended to genetic studies while those 
# below 0.25 are not recommended (Lemos Serrote et al. 2020 Gene). 
# 
# This may introduce a bias in which SNPs with more equal allele freq (p = 0.5) 
# are more likely to be kept. So we use method "random".
EYR <- gl.filter.secondaries(EYR, method = "random")

EYR  # 782 genotypes,  35,663 binary SNPs

# Downsample loci to upload to dartR.sexlinked
EYR <- gl.keep.loc(EYR, 
                   loc.list = sample(EYR@loc.names, 
                                     size = 1000, 
                                     replace = FALSE))

EYR  # 782 genotypes,  1,000 binary SNPs

# Save as file in package
usethis::use_data(EYR)
