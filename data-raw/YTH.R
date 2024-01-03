## code to prepare `YTH` dataset goes here

# usethis::use_data(DATASET, overwrite = TRUE)

# Read raw data as genlight object
library(dartR)
dir <- "C:/Users/Diangie/My Drive/Projects/1.MateChoice_BreedingManagement/"
dart_file <- "Report-DLich19-4793/Report_DLich19-4793_2_moreOrders_SNP_2.csv"
cov_file <- "Analyses2020/1_2020.02.12_Sex_check_DArT_files/YTH_covariates_complete_dartR_repeated_geneticsexes.csv"

YTH <- gl.read.dart(filename = paste(dir,
                                     dart_file, 
                                     sep = ""),
                    covfilename = paste(dir,
                                        cov_file,
                                        sep = ""))

YTH  # 646 genotypes,  118,750 binary SNPs

# Remove repeated individuals
YTH <- gl.drop.ind(YTH, ind.list = c("043-04232.1",
                                   "043-00549.1",
                                   "043-00541.1",
                                   "043-04220.1",
                                   "043-04201.1"),
                  recalc = TRUE,
                  mono.rm = TRUE)

YTH  # 641 genotypes,  118,732 binary SNPs

YTH@other$ind.metrics <- droplevels(YTH@other$ind.metrics)


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
YTH <- gl.filter.secondaries(YTH, method = "random")

YTH  # 641 genotypes,  74,470 binary SNPs

# Downsample loci to upload to dartR.sexlinked
YTH_test <- gl.keep.loc(YTH, 
                   loc.list = sample(YTH@loc.names, 
                                     size = 1000, 
                                     replace = FALSE))

YTH  # 641 genotypes,  1,000 binary SNPs

# Save as file in package
usethis::use_data(YTH)
