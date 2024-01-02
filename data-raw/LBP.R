## code to prepare `LBP` dataset goes here

# usethis::use_data(DATASET, overwrite = TRUE)

# Read raw data as genlight object
library(dartR)
dir <- "C:/Users/Diangie/My Drive/Projects/SwissArmyKnife/SCRIPTS/New/LBP/filtering_regime_using_filter.sex.linked/"

LBP <- gl.read.dart(filename = paste(dir,
                                     "LBP_DArT.csv", 
                                     sep = ""),
                    covfilename = paste(dir,
                                        "LBP_DArT_covariates.csv",
                                        sep = ""))

LBP  # 400 genotypes,  10,149 binary SNPs

# Format sexes
LBP$other$ind.metrics$sex <- gsub("Female", "F", LBP$other$ind.metrics$sex)
LBP$other$ind.metrics$sex <- gsub("Male", "M", LBP$other$ind.metrics$sex)
LBP$other$ind.metrics$sex <- gsub("Unknown", "", LBP$other$ind.metrics$sex)

# Remove individuals for which silicoDArT sexes were wrong
remove <- c("Y153", "LM022", "LM024", "LM117", "LM028", "LM129", "LM065", 
            "LM111", "ARCHERON_2016_F", "YLB_O3E_GR70", "BBRN28_AD_F",
            "BBC4_AD_F", "LM020", "LM014", "LM023", "LM018", "LM116", "LM132", 
            "LM066", "LM121", "LM106", "YLB_M6B_BLX1_M", "YLB_O3E_BL21", 
            "YLB_L9C_GR27_M", "YLB_E3I_GRYRT_M", "DSN6A_JUV_M", "LM059","LM078")

LBP <- gl.keep.ind(LBP, 
                   ind.list = LBP@ind.names[!(LBP@ind.names %in% remove)], 
                   mono.rm = TRUE)

LBP  # 376 genotypes,  9,508 binary SNPs

LBP@other$ind.metrics <- droplevels(LBP@other$ind.metrics)


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
LBP <- gl.filter.secondaries(LBP, method = "random")

LBP  # 376 genotypes,  8,436 binary SNPs

# Downsample loci to upload to dartR.sexlinked
LBP_test <- gl.keep.loc(LBP, 
                   loc.list = sample(LBP@loc.names, 
                                     size = 1000, 
                                     replace = FALSE))

LBP  # 782 genotypes,  1,000 binary SNPs

# Save as file in package
usethis::use_data(LBP)
