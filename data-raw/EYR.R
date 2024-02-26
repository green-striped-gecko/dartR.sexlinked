## code to prepare `EYR` dataset

# usethis::use_data(DATASET, overwrite = TRUE)

# Read raw data as genlight object
library(dartR)
dir <- "C:/Users/Diangie/My Drive/Projects/SwissArmyKnife/SCRIPTS/eyr/"
dart_file <- "Report_DYro21-6107_5_moreOrders_SNP_mapping_3_uniqIDs.csv"

data0 <- gl.read.dart(filename = paste(dir,
                                     dart_file, 
                                     sep = ""),
                    covfilename = paste(dir,
                                        "EYR2021_cov.NONAS.csv",
                                        sep = ""))

data0  # 899 genotypes,  54,823 binary SNPs

# Remove individuals to not be used
data <- gl.drop.pop(data0, pop.list = c("DoNotUse",
                                     "KeilorsTrack",
                                     "PileOfBricks"),
                   recalc = TRUE,
                   mono.rm = TRUE) 

data  # 782 genotypes,  53,324 binary SNPs

data@other$ind.metrics <- droplevels(data@other$ind.metrics)


# Change pops from ind.metrics
data@other$ind.metrics[data@other$ind.metrics$pop == "RailwayDam", "pop"] <- "Muckleford"
data@other$ind.metrics[data@other$ind.metrics$pop == "Sedgwick", "pop"] <- "Crusoe"

data@other$ind.metrics <- droplevels(data@other$ind.metrics)


# Change pops from @pops
data@pop[data@pop == "RailwayDam"] <- "Muckleford"
data@pop[data@pop == "Sedgwick"] <- "Crusoe"

data@pop <- droplevels(data@pop)


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
data <- gl.filter.secondaries(data, method = "random")

data  # 782 genotypes,  35,663 binary SNPs

# Identify autosomes in order to downsample them
data_sex <- filter.sex.linked(data, system = "zw")
# **FINISHED** Total of analyzed loci: 35663.
# Found 3722 sex-linked loci:
#   147 W-linked loci
# 2466 sex-biased loci
# 794 Z-linked loci
# 315 ZW gametologs.
# And 31941 autosomal loci.

aut <- sample(data_sex$autosomal@loc.names, 
              size = 850, 
              replace = FALSE)
W <- sample(data_sex$w.linked@loc.names, 
            size = 16, 
            replace = FALSE)
Z <- sample(data_sex$z.linked@loc.names, 
            size = 35, 
            replace = FALSE)
g <- sample(data_sex$gametolog@loc.names, 
            size = 20, 
            replace = FALSE)
bia <- sample(data_sex$sex.biased@loc.names, 
              size = 79, 
              replace = FALSE)


# Downsample loci to upload to dartR.sexlinked
EYR <- gl.keep.loc(data, 
                   loc.list = c(aut, 
                                Z,
                                bia,
                                W,
                                g))


EYR  # 782 genotypes,  1,000 binary SNPs,   size: 20.2 Mb


# Save as file in package
usethis::use_data(EYR, overwrite = TRUE)
