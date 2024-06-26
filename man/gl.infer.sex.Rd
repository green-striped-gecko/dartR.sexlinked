% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gl.infer.sex.r
\name{gl.infer.sex}
\alias{gl.infer.sex}
\title{Uses sex-linked loci to infer sex of individuals}
\usage{
gl.infer.sex(gl_sexlinked, system = NULL, seed = NULL)
}
\arguments{
\item{gl_sexlinked}{The output of function gl.keep.sexlinked (complete 
list with 5 elements). See explanation in "Details" section [required].}

\item{system}{String that declares the sex-determination system of the 
species: 'zw' or 'xy' [required].}

\item{seed}{User-defined integer for repeatability purposes. If not provided
by user, it is chosen randomly by default. See "Details" section.}
}
\value{
A dataframe.
}
\description{
This function uses the output of function gl.keep.sexlinked (list of 5 
elements) to infer the sex of all individuals in the dataset. 
It uses 3 types of sex-linked loci (W-/Y-linked, Z-/X-linked, and 
gametologs), assigns a preliminary genetic sex for each type of sex-linked
loci available, and outputs an agreed sex.

This function produces as output a dataframe with individuals in rows and 11 
columns.
}
\details{
Parameter \code{gl_sexlinked} must be the name of the output object (a 
list of 5 elements) produced by function \code{gl.keep.sexlinked}. Parameter 
\code{seed} must be an integer that will be used on the KMeans algorithm used
by the function. We highly recommend choosing the seed to guarantee 
repeatability.

Note that this function was created with the explicit intent that a human 
checks the evidence for the sex assignments that do NOT agree for all 
types of sex-linked loci (called "indefinite sex assignments" and denoted 
as "*M" or "*F" in the last column of dataframe output). This human can then 
use their criterion to validate these assignments.


\strong{ Function's output }

This function creates a dataframe with one row per individual and 11 columns: \itemize{
\item {id > Individuals' ID.}

\item {w.linked.sex or y.linked.sex > Sex inferred using w-linked or y-linked
loci.}
\item {#called > Number of W-linked or Y-linked loci for which the individual 
had a called genotype (cf. missing genotype).}
\item {#missing > Number of W-linked or Y-linked loci for which the 
individual had a missing genotype (cf. called genotype).} 

\item {z.linked.sex or x.linked.sex > Sex inferred using z-linked or x-linked
loci.}
\item {#Hom.z or #Hom.x > Number of z-linked or x-linked loci for which the 
individual is homozygous.}
\item {#Het.z or #Het.x > Number of z-linked or x-linked loci for which the 
individual is heterozygous.}

\item {gametolog.sex > Sex inferred using gametologs.} 
\item {#Hom.g > Number of gametologous loci for which the individual is 
homozygous.}
\item {#Het.g > Number of gametologous loci for which the individual is 
heterozygous.}

\item {agreed.sex > Agreed sex: 'F' or 'M' if all preliminary sex-assignments
match (i.e., definite sex assignment), and '*F' or '*M' if NOT all 
preliminary sex-assignments match (i.e., indefinite sex assignment).} 
}
}
\examples{
LBP_sexLinked <- gl.keep.sexlinked(x = LBP, system = "xy", plot.display = TRUE, ncores = 1)
inferred.sexes <- gl.infer.sex(gl_sexlinked = LBP_sexLinked, system = "xy", seed = 100)
inferred.sexes

}
\references{
\itemize{
\item Robledo‐Ruiz, D. A., Austin, L., Amos, J. N., Castrejón‐Figueroa, J.,
 Harley, D. K., Magrath, M. J., Sunnucks, P., & Pavlova, A. (2023). 
 Easy‐to‐use R functions to separate reduced‐representation genomic datasets
  into sex‐linked and autosomal loci, and conduct sex assignment. Molecular 
  Ecology Resources, 00, 1-21.
 }
}
\author{
Custodian: Diana Robledo-Ruiz  -- Post to
  \url{https://groups.google.com/d/forum/dartr}
}
