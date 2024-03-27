## code to prepare `YTH` dataset goes here

# usethis::use_data(DATASET, overwrite = TRUE)

# Read raw data as genlight object
library(dartR)
dir <- "C:/Users/Diangie/My Drive/Projects/1.MateChoice_BreedingManagement/"
dart_file <- "Report-DLich19-4793/Report_DLich19-4793_2_moreOrders_SNP_2.csv"
cov_file <- "Analyses2020/1_2020.02.12_Sex_check_DArT_files/YTH_covariates_complete_dartR_repeated_geneticsexes.csv"

data0 <- gl.read.dart(filename = paste(dir,
                                     dart_file, 
                                     sep = ""),
                    covfilename = paste(dir,
                                        cov_file,
                                        sep = ""))

data0  # 646 genotypes,  118,750 binary SNPs

# Remove repeated individuals
data <- gl.drop.ind(data0, ind.list = c("043-04232.1",
                                   "043-00549.1",
                                   "043-00541.1",
                                   "043-04220.1",
                                   "043-04201.1",
                                   "A81754", 
                                   "043-04282",
                                   "042-07991", 
                                   "W41", 
                                   "W93",
                                   "A81196", 
                                   "043-04222", 
                                   "043-04285", 
                                   "043-12508", 
                                   "(042?) 59323"),
                  recalc = TRUE,
                  mono.rm = TRUE)

data  # 631 genotypes,  118,726 binary SNPs

# Remove individuals with uncertain population
data <- gl.drop.pop(data, pop.list = c("Gippslandicus/Melanops?",
                                       "Gippslandicus/Cassidix?",
                                       "Gippslandicus?",
                                       "Hybrid"),
                    recalc = TRUE,
                    mono.rm = TRUE)

data  # 623 genotypes,  118,583 SNPs

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

data  # 623 genotypes,  74,382 SNPs

# Identify autosomes in order to downsample them
data_sex <- filter.sex.linked(data, system = "zw")
# **FINISHED** Total of analyzed loci: 74382.
# Found 3204 sex-linked loci:
#   57 W-linked loci
# 2147 sex-biased loci
# 987 Z-linked loci
# 13 ZW gametologs.
# And 71178 autosomal loci.

aut <- sample(data_sex$autosomal@loc.names, 
              size = 500, 
              replace = FALSE)
Z <- sample(data_sex$z.linked@loc.names, 
              size = 200, 
              replace = FALSE)
bia <- sample(data_sex$sex.biased@loc.names, 
              size = 230, 
              replace = FALSE)


# Downsample loci to upload to dartR.sexlinked
YTH <- gl.keep.loc(data, 
                   loc.list = c(aut, 
                                Z,
                                bia,
                                data_sex$w.linked@loc.names,    # all w-linked
                                data_sex$gametolog@loc.names))  # all gametologs

YTH  # 623 genotypes,  1,000 binary SNPs, size 50 Mb


# Remove low quality individuals
gl.report.callrate(YTH, 
                   method = "ind")
YTH <- gl.filter.callrate(YTH, 
                          method = "ind", 
                          threshold = 0.6, 
                          mono.rm = TRUE)

YTH  # 609 genotypes,  994 SNPs , size: 49.9 Mb

# Save as file in package
usethis::use_data(YTH, overwrite = TRUE)









pca <- gl.pcoa(YTH)
gl.pcoa.plot(x = YTH, glPca = pca)
Ho <- gl.report.heterozygosity(YTH, method = "ind")
x <- merge(YTH@other$ind.metrics, Ho, 
           by.x = "id", by.y = "ind.name")
t.test(x[x$sex == "F", "Ho"], x[x$sex == "M", "Ho"])
# Welch Two Sample t-test
# 
# data:  x[x$sex == "F", "Ho"] and x[x$sex == "M", "Ho"]
# t = -52.817, df = 427.05, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.09005245 -0.08359053
# sample estimates:
#   mean of x  mean of y 
# 0.07292651 0.15974800 


YTH_fil <- filter.sex.linked(YTH, system = "zw")
Ho_fil <- gl.report.heterozygosity(YTH_fil$autosomal, method = "ind")
x_fil <- merge(YTH_fil$autosomal@other$ind.metrics, Ho_fil, 
           by.x = "id", by.y = "ind.name")
t.test(x_fil[x_fil$sex == "F", "Ho"], x_fil[x_fil$sex == "M", "Ho"])

x_fil1 <- x_fil[!x_fil$id %in% ex, ]

ex <- c("A81754", 
  "043-04282",
  "042-07991", 
  "W41", 
  "W93",
  "A81196", 
  "043-04222", 
  "043-04285", 
  "043-12508", 
  "(042?) 59323")
t.test(x_fil1[x_fil1$sex == "F", "Ho"], x_fil1[x_fil1$sex == "M", "Ho"])



