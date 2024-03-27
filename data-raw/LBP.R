## code to prepare `LBP` dataset goes here

# usethis::use_data(DATASET, overwrite = TRUE)

# Read raw data as genlight object
library(dartR)
dir <- "C:/Users/Diangie/My Drive/Projects/SwissArmyKnife/SCRIPTS/New/LBP/filtering_regime_using_filter.sex.linked/"

data0 <- gl.read.dart(filename = paste(dir,
                                     "LBP_DArT.csv", 
                                     sep = ""),
                    covfilename = paste(dir,
                                        "LBP_DArT_covariates.csv",
                                        sep = ""))

data0  # 400 genotypes,  10,149 binary SNPs

# Format sexes
data0$other$ind.metrics$sex <- gsub("Female", "F", data0$other$ind.metrics$sex)
data0$other$ind.metrics$sex <- gsub("Male", "M", data0$other$ind.metrics$sex)
data0$other$ind.metrics$sex <- gsub("Unknown", "", data0$other$ind.metrics$sex)

# Remove individuals for which silicoDArT sexes were wrong
remove <- c("Y153", "LM022", "LM024", "LM117", "LM028", "LM129", "LM065", 
            "LM111", "ARCHERON_2016_F", "YLB_O3E_GR70", "BBRN28_AD_F",
            "BBC4_AD_F", "LM020", "LM014", "LM023", "LM018", "LM116", "LM132", 
            "LM066", "LM121", "LM106", "YLB_M6B_BLX1_M", "YLB_O3E_BL21", 
            "YLB_L9C_GR27_M", "YLB_E3I_GRYRT_M", "DSN6A_JUV_M", "LM059","LM078")

data <- gl.keep.ind(data0, 
                   ind.list = data0@ind.names[!(data0@ind.names %in% remove)], 
                   mono.rm = TRUE)

data  # 376 genotypes,  9,508 binary SNPs

data@other$ind.metrics <- droplevels(data@other$ind.metrics)


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

data  # 376 genotypes,  8,436 binary SNPs

# Identify autosomes in order to downsample them
data_sex <- filter.sex.linked(data, system = "xy")
# **FINISHED** Total of analyzed loci: 8436.
# Found 73 sex-linked loci:
#   1 Y-linked loci
# 6 sex-biased loci
# 65 X-linked loci
# 1 XY gametologs.
# And 8363 autosomal loci.

aut <- sample(data_sex$autosomal@loc.names, 
              size = 927, 
              replace = FALSE)


# Downsample loci to upload to dartR.sexlinked
YTH <- gl.keep.loc(data, 
                   loc.list = c(aut, 
                                data_sex$x.linked@loc.names,    # all x-linked
                                data_sex$sex.biased@loc.names,  # all sex-biased
                                data_sex$y.linked@loc.names,    # all y-linked
                                data_sex$gametolog@loc.names))  # all gametologs

LBP  # 376 genotypes,  1,000 binary SNPs, 5.2 Mb

# Save as file in package
usethis::use_data(LBP, overwrite = TRUE)
